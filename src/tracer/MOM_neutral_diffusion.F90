!> A column-wise toolbox for implementing neutral (horizontal) diffusion
module MOM_neutral_diffusion

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_EOS, only : EOS_type, calculate_compress, calculate_density_derivs
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type

implicit none ; private

#include <MOM_memory.h>

public neutral_diffusion, neutral_diffusion_init, neutral_diffusion_end
public neutral_diffusion_calc_coeffs
public neutralDiffusionUnitTests

type, public :: neutral_diffusion_CS ; private
  integer :: nkp1   ! Number of interfaces for a column = nk + 1
  integer :: nkp1X2 ! Number of intersecting interfaces between columns = 2 * nkp1
end type neutral_diffusion_CS

! This include declares and sets the variable "version".
#include "version_variable.h"
character(len=40)  :: mod = "MOM_neutral_diffusion" ! This module's name.

contains

subroutine neutral_diffusion(nk, Dl, Dr, hl, hr, Tl, Tr, Sl, Sr, EOS, H_to_Pa)
  integer,             intent(in)    :: nk  !< Number of levels per column (assumed equal)
  real,                intent(in)    :: Dl  !< Bottom depth of left column (m)
  real,                intent(in)    :: Dr  !< Bottom depth of right column (m)
  real, dimension(nk), intent(in)    :: hl  !< Level thickness of left column (m or kg/m2)
  real, dimension(nk), intent(in)    :: hr  !< Level thickness of right column (m or kg/m2)
  real, dimension(nk), intent(in)    :: Tl  !< Potential temperature of left column (degC)
  real, dimension(nk), intent(in)    :: Tr  !< Potential temperature of right column (degC)
  real, dimension(nk), intent(in)    :: Sl  !< Salinity of left column (ppt)
  real, dimension(nk), intent(in)    :: Sr  !< Salinity of right column (ppt)
  type(EOS_type),      pointer       :: EOS !< Equation of state structure
  real,                intent(in)    :: H_to_Pa !< Converts h units to Pascals
  ! Local variables
  real, dimension(nk+1) :: Pil ! Interface pressures of left column (Pa)
  real, dimension(nk+1) :: Pir ! Interface pressures of right column (Pa)
  real, dimension(nk+1) :: Til !< Interface potential temperature of left column (degC)
  real, dimension(nk+1) :: Tir !< Interface potential temperature of right column (degC)
  real, dimension(nk+1) :: Sil !< Interface salinity of left column (ppt)
  real, dimension(nk+1) :: Sir !< Interface salinity of right column (ppt)
  real, dimension(nk+1) :: rhoil !< Interface densities of left column (kg/m3)
  real, dimension(nk+1) :: rhoir !< Interface densities of right column (kg/m3)
  real, dimension(nk+1) :: PirL !< Right column interface pressure projected on to left column (Pa)
  real, dimension(nk+1) :: PilR !< Right column interface pressure projected on to right column (Pa)
  real, dimension(2*nk+1) :: dPl !< Projected thickness of union layers in left column (Pa)
  real, dimension(2*nk+1) :: dPr !< Projected thickness of union layers in right column (Pa)
  integer, dimension(2*nk+1) :: Kl !< Left column indexes of projected union layers
  integer, dimension(2*nk+1) :: Kr !< Right column indexes of projected union layers
  real, dimension(nk+1) :: dRdTl !< Interface thermal expansion coefficient of left column (kg/m3/degC)
  real, dimension(nk+1) :: dRdSl !< Interface haline expandion coefficient of left column (kg/m3/ppt)
  real, dimension(nk+1) :: dRdTr !< Interface thermal expansion coefficient of right column (kg/m3/degC)
  real, dimension(nk+1) :: dRdSr !< Interface haline expandion coefficient of right column (kg/m3/ppt)
  integer :: k

  ! Find interface positions in meters for each column (ultimately needed for equation of state)
  Pil(1) = 0. ; Pir(1) = 0. ! Pressure loading from ice-shelves? -AJA
  do k = 1, nk
    Pil(k+1) = Pil(k) + H_to_Pa * hl(k)
    Pir(k+1) = Pir(k) + H_to_Pa * hr(k)
  enddo

  ! Create interface T and S for each column
  call interface_TS(nk, Tl, Sl, Til, Sil)
  call interface_TS(nk, Tr, Sr, Tir, Sir)

  ! For each interface on the left column, calculate the expansion coefficients
  ! that will be used to calculate potential density differences referenced to
  ! the left column interfaces.
  call calculate_density_derivs(Til, Sil, Pil, dRdTl, dRdSl, 1, nk+1, EOS)
  ! For each interface on the right column, calculate the expansion coefficients
  ! that will be used to calculate potential density differences referenced to
  ! the right column interfaces.
  call calculate_density_derivs(Tir, Sir, Pir, dRdTr, dRdSr, 1, nk+1, EOS)

  call find_neutral_surface_position(nk, Til, Sil, dRdTl, dRdSl, Tir, Sir, Pir, PilR)

end subroutine neutral_diffusion

!> Returns interface T and S for a column
subroutine interface_TS(nk, T, S, Ti, Si)
  integer,               intent(in)    :: nk !< Number of levels
  real, dimension(nk),   intent(in)    :: T  !< Layer potential temperature (degC)
  real, dimension(nk),   intent(in)    :: S  !< Layer salinity (ppt)
  real, dimension(nk+1), intent(inout) :: Ti !< Interface potential temperature (degC)
  real, dimension(nk+1), intent(inout) :: Si !< Interface salinity (ppt)
  ! Local variables
  integer :: k

  ! We use simple averaging for internal interfaces and piecewise-constant at the top and bottom.
  Ti(1) = T(1) ; Si(1) = S(1)
  do k = 2, nk
    Ti(k) = 0.5*( T(k-1) + T(k) )
    Si(k) = 0.5*( S(k-1) + S(k) )
  enddo
  Ti(nk+1) = T(nk) ; Si(nk+1) = S(nk)

end subroutine interface_TS

!> Returns position within right column of neutral surface corresponding to each left-column interface
subroutine find_neutral_surface_position(nk, Tl, Sl, dRdT, dRdS, Tr, Sr, Pr, Po)
  integer,               intent(in)    :: nk   !< Number of levels
  real, dimension(nk+1), intent(in)    :: Tl   !< Left-column interface potential temperature (degC)
  real, dimension(nk+1), intent(in)    :: Sl   !< Left-column interface salinity (ppt)
  real, dimension(nk+1), intent(in)    :: dRdT !< Left-column dRho/dT (kg/m3/degC)
  real, dimension(nk+1), intent(in)    :: dRdS !< Left-column dRho/dS (kg/m3/ppt)
  real, dimension(nk+1), intent(in)    :: Tr   !< Right-column interface potential temperature (degC)
  real, dimension(nk+1), intent(in)    :: Sr   !< Right-column interface salinity (ppt)
  real, dimension(nk+1), intent(in)    :: Pr   !< Right-column interface pressure (Pa)
  real, dimension(nk+1), intent(inout) :: Po   !< Position of neutral surface (Pa)
  ! Local variables
  integer kl, kr, krm1
  real dRkr, dRkrm1, Prm1, wghtU, wghtD, Pint

  ! Initialize variables for the search
  kr = 1
  Prm1 = Pr(1)

  ! For each left-column interface, search for the right layer whose interfaces (kr-1 and kr) bracket the left interface.
  do kl = 1, nk+1
    ! Search downward until the kr interface is denser than the kl interface
    krm1 = max(kr-1, 1)
    dRkrm1 = dRdT(kl) * ( Tr(krm1) - Tl(kl) ) + dRdS(kl) * ( Sr(krm1) - Sl(kl) ) !  Potential density difference, (kr-1) - kl
    dRkr = dRdT(kl) * ( Tr(kr) - Tl(kl) ) + dRdS(kl) * ( Sr(kr) - Sl(kl) ) !  Potential density difference, kr - kl
    kr_search: do while ( dRkr<=0. .and. (kr <= nk) )
      ! The potential density at kr is lighter than the kl interface so we increment kr
      Prm1 = Pr(kr)
      dRkrm1 = dRkr
      kr = kr + 1
      dRkr = dRdT(kl) * ( Tr(kr) - Tl(kl) ) + dRdS(kl) * ( Sr(kr) - Sl(kl) ) !  Potential density difference, kr - kl
    enddo kr_search
    ! At this point:
    !  1) If kr=1 we are at the surface and dRkr>0 but dRkrm1 could be either sign;
    !  2) In the interior, dRkrm1<=0. and dRkr>0;
    !  3) At the bottom,  dRkr<0.
    ! Note that dRkrm1 + dRkr >0 unless this next test fails.
    if ( dRkrm1 > dRkr ) stop 'find_neutral_surface_position: Houston, we have a problem! dRkrm1 >= dRkr'
    ! Linearly interpolate for the position between Pr(kr-1) and Pr(kr) where the density difference between right and
    ! left is zero; Pint = ( dRkr * Pr(kr) - dRkrm1 * Pr(kr-1) ) / ( dRdk - dRkrm1 ) .
    if ( dRkr - dRkrm1 > 0. ) then
      wghtU = -dRkrm1 / ( dRkr - dRkrm1 )
      wghtD = dRkr / ( dRkr - dRkrm1 )
      ! Ostensibly wghtU + wghtD = 1.
      ! The following is an accurate way to do Pint = wghtD * Pr(kr) + wghtU * Pr(kr-1) 
      if ( wghtU < 0.5 ) then
        Pint = Prm1 + max( wghtU, 0. ) * ( Pr(kr) - Prm1 )
      elseif ( wghtD < 0.5 ) then
        Pint = Pr(kr) + max( wghtD, 0. ) * ( Prm1 - Pr(kr) )
      else
        Pint = 0.5 * ( Prm1 + Pr(kr) )
      endif
    else ! dRkr - dRkrm1 = 0
      Pint = 0.5 * ( Prm1 + Pr(kr) )
    endif
    if ( Pint < Prm1 ) stop 'find_neutral_surface_position: Houston, we have a problem! Pint < Pr(kr-1)'
    if ( Pint > Pr(kr) ) stop 'find_neutral_surface_position: Houston, we have a problem! Pint > Pr(kr)'
    Po(kl) = Pint
  enddo

end subroutine find_neutral_surface_position

subroutine neutral_diffusion_init(Time, G, param_file, diag, CS)
  type(time_type), target,  intent(in)    :: Time
  type(ocean_grid_type),    intent(in)    :: G
  type(diag_ctrl), target,  intent(inout) :: diag
  type(param_file_type),    intent(in)    :: param_file
  type(neutral_diffusion_CS), pointer       :: CS
  character(len=256) :: mesg    ! Message for error messages.

  if (associated(CS)) then
    call MOM_error(WARNING, "neutral_diffusion_init called with associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
! call get_param(param_file, mod, "KHTR", CS%KhTr, &
!                "The background along-isopycnal tracer diffusivity.", &
!                units="m2 s-1", default=0.0)
end subroutine neutral_diffusion_init

!> Calculates remapping factors for u/v columns used to map adjoining columns to
!! shared coordinate space
subroutine neutral_diffusion_calc_coeffs(G, CS)
  type(ocean_grid_type),    intent(in)    :: G  !< Ocean grid structure
  type(neutral_diffusion_CS), pointer     :: CS !< Neutral slopes structure

end subroutine neutral_diffusion_calc_coeffs

subroutine neutral_diffusion_end(CS)
  type(neutral_diffusion_CS), pointer :: CS

  if (associated(CS)) deallocate(CS)

end subroutine neutral_diffusion_end

!> Returns true if unit tests of neutral_diffusion functions fail. Otherwise returns false.
logical function neutralDiffusionUnitTests()
  integer, parameter :: nk = 4
  real, dimension(nk+1) :: Til, Tir1, Tir2, Tir3, Tir4, Tir5, Tir6 ! Test interface temperatures
  real, dimension(nk+1) :: Sil, Sir1 ! Test interface salinities
  real, dimension(nk+1) :: Pil, Pir1 ! Test interface positions
  real, dimension(nk+1) :: dRdT, dRdS ! Test interface expansion coefficients
  real, dimension(nk+1) :: Po, Po1, Po2, Po3, Po4, Po5, Po6 ! Test positions
  data Pil / 0., 10., 20., 30., 40. /
  data Til / 10., 7.5, 5., 2.5, 0. /
  data Sil / 0., 0., 0., 0., 0. /
  data dRdT / -0.25, -0.25, -0.25, -0.25, -0.25 /
  data dRdS / 0.5, 0.5, 0.5, 0.5, 0.5 /
  ! Same grid, warmer values
  data Pir1 / 0., 10., 20., 30., 40. /
  data Tir1 / 12., 9.5, 7., 4.5, 2. /
  data Sir1 / 0., 0., 0., 0., 0. /
  data Sir1 / 0., 0., 0., 0., 0. /
  data Po1 / 8., 18., 28., 38., 40. / ! Correct answer
  ! Same grid, cooler values
  data Tir2 / 8., 5.5, 3., 0.5, -2. /
  data Po2 / 0., 2., 12., 22., 32. / ! Correct answer
  ! Same grid, all warmer values
  data Tir3 / 12., 12., 11., 11., 10.5 /
  data Po3 / 40., 40., 40., 40., 40. / ! Correct answer
  ! Same grid, all cooler values
  data Tir4 / -2., -2., -3., -3., -4. /
  data Po4 / 0., 0., 0., 0., 0. / ! Correct answer
  ! Same grid, spanning values
  data Tir5 / 15., 15., 15., -5., -5. /
  data Po5 / 22.5, 23.75, 25., 26.25, 27.5 / ! Correct answer
  ! Same grid, encompassed values
  data Tir6 / 7., 6.5, 6., 5.5, 5. /
  data Po6 / 0., 0., 40., 40., 40. / ! Correct answer
  integer :: k

  neutralDiffusionUnitTests = .false. ! Normally return false
  write(*,'(a)') '===== MOM_neutral_diffusion: neutralDiffusionUnitTests =================='

  call find_neutral_surface_position(nk, Til, Sil, dRdT, dRdS, Til, Sir1, Pir1, Po)
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. simpleTest(nk, Po, Pil, 'Project left->right 0, identical values')

  call find_neutral_surface_position(nk, Til, Sil, dRdT, dRdS, Tir1, Sir1, Pir1, Po)
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. simpleTest(nk, Po, Po1, 'Project left->right 1, slightly warmer values')

  call find_neutral_surface_position(nk, Til, Sil, dRdT, dRdS, Tir2, Sir1, Pir1, Po)
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. simpleTest(nk, Po, Po2, 'Project left->right 2, slightly cooler values')

  call find_neutral_surface_position(nk, Til, Sil, dRdT, dRdS, Tir3, Sir1, Pir1, Po)
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. simpleTest(nk, Po, Po3, 'Project left->right 3, all warmer values')

  call find_neutral_surface_position(nk, Til, Sil, dRdT, dRdS, Tir4, Sir1, Pir1, Po)
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. simpleTest(nk, Po, Po4, 'Project left->right 4, all cooler values')

  call find_neutral_surface_position(nk, Til, Sil, dRdT, dRdS, Tir5, Sir1, Pir1, Po)
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. simpleTest(nk, Po, Po5, 'Project left->right 5, spanning values')

  call find_neutral_surface_position(nk, Til, Sil, dRdT, dRdS, Tir6, Sir1, Pir1, Po)
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. simpleTest(nk, Po, Po6, 'Project left->right 6, encompassed values')

  write(*,'(a)') '=========================================================='

  contains

  !> Writes results to screen and return true is test fails
  logical function simpleTest(nk, Po, Ptrue, title)
    integer,          intent(in) :: nk !< Number of layers
    real,             intent(in) :: Po(nk+1) !< Calculated answer
    real,             intent(in) :: Ptrue(nk+1) !< True answer
    character(len=*), intent(in) :: title !< Title for messages
    ! Local variables
    integer :: k

    write(*,'(a)') title
    simpleTest = .false.
    do k = 1,nk+1
      if (Po(k) /= Ptrue(k)) simpleTest = .true.
      write(*,'(a,i2,2(x,a,f20.16),x,a,1pe22.15)') 'k=',k,'Po=',Po(k),'Ptrue=',Ptrue(k),'err=',Po(k)-Ptrue(k)
    enddo

  end function simpleTest

end function neutralDiffusionUnitTests

end module MOM_neutral_diffusion
