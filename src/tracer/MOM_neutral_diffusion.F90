!> A column-wise toolbox for implementing neutral (horizontal) diffusion
module MOM_neutral_diffusion

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_EOS, only : calculate_compress, EOS_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type

implicit none ; private

#include <MOM_memory.h>

public neutral_diffusion, neutral_diffusion_init, neutral_diffusion_end
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
  real, dimension(nk+1) :: Pil !< Interface pressures of left column (Pa)
  real, dimension(nk+1) :: Pir !< Interface pressures of right column (Pa)
  real, dimension(nk+1) :: Til !< Interface potential temperature of left column (degC)
  real, dimension(nk+1) :: Tir !< Interface potential temperature of right column (degC)
  real, dimension(nk+1) :: Sil !< Interface salinity of left column (ppt)
  real, dimension(nk+1) :: Sir !< Interface salinity of right column (ppt)
  real, dimension(nk+1) :: rhoil !< Interface densities of left column (kg/m3)
  real, dimension(nk+1) :: rhoir !< Interface densities of right column (kg/m3)
  real, dimension(nk+1) :: drdpil !< Interface compressibility of left column (kg/m3/pa)
  real, dimension(nk+1) :: drdpir !< Interface compressibility of right column (kg/m3/pa)
  real, dimension(nk+1) :: PirL !< Right column interface pressure projected on to left column (Pa)
  real, dimension(nk+1) :: PilR !< Right column interface pressure projected on to right column (Pa)
  real, dimension(2*nk+1) :: dPl !< Projected thickness of union layers in left column (Pa)
  real, dimension(2*nk+1) :: dPr !< Projected thickness of union layers in right column (Pa)
  integer, dimension(2*nk+1) :: Kl !< Left column indexes of projected union layers
  integer, dimension(2*nk+1) :: Kr !< Right column indexes of projected union layers
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
  
  ! Calculate interface densities using the interface pressure for each column (in-situ)
  ! and the compressibility d/dp rho.
  call calculate_compress(Til, Sil, Pil, rhoil, drdpil, 1, nk+1, EOS)
  call calculate_compress(Tir, Sir, Pir, rhoir, drdpir, 1, nk+1, EOS)

  ! For each interface of the right column, using the density of that interface referenced
  ! to that local pressure (isoneutral), find the position within the left column where the
  ! potential density referenced to the same pressure equals that of the right interface.
  call projected_positions(nk, rhoir, Pir, rhoil, Pil, drdpil, PirL)
  ! For each interface of the left column, using the density of that interface referenced
  ! to that local pressure (isoneutral), find the position within the right column where the
  ! potential density referenced to the same pressure equals that of the left interface.
  call projected_positions(nk, rhoil, Pil, rhoir, Pir, drdpir, PilR)

  ! Find corresponding indexes and thickness in the position space of the left column for
  ! each union layer
  call project_thicknesses(nk, Pil, PirL, dPl, Kl)
  ! Find corresponding indexes and thickness in the position space of the right column for
  ! each union layer
  call project_thicknesses(nk, Pir, PilR, dPr, Kr)

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

!> Returns position within left column of right column in-situ densities
subroutine projected_positions(nk, rhoir, Pir, rhoil, Pil, drdpil, Po)
  integer,               intent(in)    :: nk     !< Number of levels
  real, dimension(nk+1), intent(in)    :: rhoir  !< Interface in-situ density of right column (Pa)
  real, dimension(nk+1), intent(in)    :: Pir    !< Interface pressure of right column (Pa)
  real, dimension(nk+1), intent(in)    :: rhoil  !< Interface in-situ density of left column (Pa)
  real, dimension(nk+1), intent(in)    :: Pil    !< Interface pressure of left column (Pa)
  real, dimension(nk+1), intent(in)    :: drdpil !< Interface compressibility of left column (kg/m3/Pa)
  real, dimension(nk+1), intent(inout) :: Po     !< Projection of right interface positions on left column (Pa)
  ! Local variables
  real :: sigma_kl, sigma_kl_p1 ! Densities of left column referenced to a pressure from right column (kg/m3)
  real :: Po_last ! The last projected pressure from right side, used to avoid tangling. (Pa)
  real :: P_target, sigma_target, P, a, b
  integer :: kl, kr
  logical, parameter :: debug = .true.

  kl = 1
  Po_last = Pil(1)
  do kr = 1, nk+1
    ! We are looking for the position within the left column of a density referenced to a
    ! pressure from the right column
    P_target = Pir(kr)
    sigma_target = rhoir(kr)

    ! Search for a left column layer whose interface densities
    ! span rhoir(kr) when referenced to Pir(kr).
    search_for_kl: do while(.true.)
      !kl = min(nk+1, kl)

      ! Density of interface kl at reference pressure, pr(kr)
      ! The use of compressibility avoids the call
      !  call calculate_density( Til(kl), Sil(kl), Pir(kr), sigma_kl )
      sigma_kl = rhoil(kl) + drdpil(kl) * ( P_target - Pil(kl) )

      ! Density of interface kl+1 at reference pressure, pr(kr)
      ! The use of compressibility avoids the call
      !  call calculate_density( Til(kl), Sil(kl), Pir(kr), sigma_kl )
      sigma_kl_p1 = rhoil(kl+1) + drdpil(kl+1) * ( P_target - Pil(kl+1) )
                                               if (debug) write(*,'(2(x,a,i2),3(x,a,f10.3))') 'kr=',kr,'kl=',kl,'sigma_kl=',sigma_kl,'sigma_kl_p1=',sigma_kl_p1,'sigma_target=',sigma_target
                                               if (debug) write(*,'(12x,3(x,a,f10.3))') 'Pil(kl)=',Pil(kl),'Pil(kl+1)=',Pil(kl+1),'P_target=',P_target

      if (sigma_target <= sigma_kl) then
        ! The density we are searching for is lighter that the left column pair
        ! which should only happen at the surface when beginning the search
        P = Pil(kl) ! Would match the surface
                                               if (debug) write(*,'(4x,a)') 'lighter or match'
        exit search_for_kl
      else
        if (sigma_kl_p1 > sigma_kl) then
          ! The left column is stably stratified
          if (sigma_target > sigma_kl_p1) then
            ! The lower of the pair is lighter than the density we need so move kl down
            if (kl < nk) then
                                               if (debug) write(*,'(4x,a)') 'stable kl,kl+1 heavier than kl+1, increment'
              kl = kl + 1
              cycle search_for_kl
            else ! kl = nk+1
                                               if (debug) write(*,'(4x,a)') 'stable kl,kl+1 heavier than kl+1, at bottom'
              ! We have reached the bottom and so sigma_target is outside of the range of left densities
             !P = Pil(nk+1) ! Would match the bottom
              P = Pil(kl+1) ! Effectively creates a new layer below bottom
              exit search_for_kl
            endif
          elseif (sigma_target >= sigma_kl) then
                                               if (debug) write(*,'(4x,a)') 'stable kl,kl+1 heavier than kl, interpolate'
            ! The density we are searching for is bounded by sigma_kl < rho < sigma_kl_p1.
            ! Interpolate for the position of the right density between the left pair.
            ! a and b are weights such that a>=0, b>=0 and a + b = 1.
            a = ( sigma_target - sigma_kl ) / ( sigma_kl_p1 - sigma_kl )
            if (a<0.5) then ! To avoids overflows
              b = 1. - a
            else
              b = ( sigma_kl_p1 - sigma_target ) / ( sigma_kl_p1 - sigma_kl )
              a = 1. - b
            endif
            P = Pil(kl+1) * a + Pil(kl) * b
            exit search_for_kl
          else ! (sigma_target < sigma_kl)
            ! This should not happen?
            stop 'sigma_target < sigma_kl is impossible!'
          endif
        elseif (sigma_kl_p1 < sigma_kl) then
                                               if (debug) write(*,'(4x,a)') 'unstable kl,kl+1'
          ! The left column is unstably stratified
          P = Po_last ! Not sure what to do here???
          exit search_for_kl
        else ! (sigma_kl_p1 == sigma_kl)
          ! The left column is neutrally stratified
          if (sigma_target > sigma_kl_p1) then
                                               if (debug) write(*,'(4x,a)') 'neutral kl,kl+1, increment'
            ! The lower of the pair is lighter than the density we need so move kl down
            kl = kl + 1
            cycle search_for_kl
          elseif (sigma_target == sigma_kl) then
                                               if (debug) write(*,'(4x,a)') 'neutral kl,kl+1, at bottom'
            P = 0.5*( Pil(kl) + Pil(kl+1) ) ! No better choice but to split the layer?
            exit search_for_kl
          else ! (sigma_target < sigma_kl)
            ! This should not happen?
            stop 'sigma_target < sigma_kl and sigma_kl_p1 == sigma_kl is impossible!'
          endif
        endif
      endif ! if (sigma_target > sigma_kl)

    enddo search_for_kl
    P = max(Po_last, P) ! Avoid possible tangles
    Po(kr) = P
    Po_last = P
 
  enddo ! kr

end subroutine projected_positions

!> Returns the projected fractional thickness of union layers into the local column, along with indexes
subroutine project_thicknesses(nk, Pi, Pp, dP, Ks)
  integer,                    intent(in)    :: nk !< Number of levels
  real,    dimension(nk+1),   intent(in)    :: Pi !< Interface pressures of local column (Pa)
  real,    dimension(nk+1),   intent(in)    :: Pp !< Positions of another column interfaces projected onto this column (Pa)
  real,    dimension(2*nk+1), intent(inout) :: dP !< Projected thickness of union layer (Pa)
  integer, dimension(2*nk+1), intent(inout) :: Ks !< Local layer index in this column spanned union interfaces
  ! Local variables
  integer :: k ! Index of union of layers
  integer :: ki ! Index of candidate local interface
  integer :: kp ! Index of candidate projected position
  integer :: kip1 ! Index of next candidate local interface
  integer :: kpp1 ! Index of next candidate projected position
  real :: Ptop, Pbot ! Top and bottom positions of union layer
  real :: Ppkp, Piki, Ppkpp1, Pikip1
  logical :: move_forward_P, move_forward_L
  logical, parameter :: debug = .false.

  ki = 1 ; kp = 1
  Ptop = max(Pp(1), Pi(1)) ! Top of intersection of columns
  Pbot = min(Pp(nk+1), Pi(nk+1)) ! Bottom of intersection of columns
  do k=1, 2*nk+1 ! Loop over the union of layers
    ! A union layer can be bounded by 4 possible pairs of interfaces:
    !   P-P - Two projected interfaces that fall outside the range of the
    !         local column
    !   L-L - Two local interfaces that fall  outside the range of the
    !         projected column
    !   L-P - A local interface above and projected interface below
    !   P-L - A projected interface above and local interface below

    kpp1 = min(nk+1, kp+1) ! Next projected position index
    kip1 = min(nk+1, ki+1) ! Next local interface index
  ! dP(k) = max(0., Pbot - Ptop)   ! Non-overlapping layers have zero weight
    Ppkp = max(Ptop, Pp(kp))
    Piki = max(Ptop, Pi(ki))
    Ppkpp1 = min(Pbot, Pp(kpp1))
    Pikip1 = min(Pbot, Pi(kip1))
                                               if (debug) write(*,'(3(x,a,i2))') 'k=',k,'ki=',ki,'kp=',kp
                                               if (debug) write(*,'(4x,3(x,a,f10.3))') 'Pp(kp)=',Pp(kp),'Pi(ki)=',Pi(ki)
                                               if (debug) write(*,'(4x,3(x,a,f10.3))') 'Pp(kp+1)=',Pp(kpp1),'Pi(ki+1)=',Pi(kip1)
                                               if (debug) write(*,'(4x,3(x,a,f10.3))') 'Ppkp=',Ppkp,'Piki=',Piki
                                               if (debug) write(*,'(4x,3(x,a,f10.3))') 'Ppkp1=',Ppkpp1,'Piki1=',Pikip1

    move_forward_P = .false.
    move_forward_L = .false.
    if     (Pp(kpp1) < Pi(ki)) then ! case P-P
                                               if (debug) write(*,'(4x,2(x,a,i3))') 'P-P kp=',kp,'ki=',ki
      move_forward_P = .true.
        dP(k) = Ppkpp1 - Ppkp
        Ks(k) = ki-1
    elseif (Pp(kp) > Pi(kip1)) then ! case L-L
                                               if (debug) write(*,'(4x,2(x,a,i3))') 'L-L kp=',kp,'ki=',ki
      move_forward_L = .true.
        dP(k) = Pikip1 - Piki
        Ks(k) = ki
    elseif (Pp(kp) < Pi(ki)) then ! case P-L
                                               if (debug) write(*,'(4x,2(x,a,i3))') 'P-L kp=',kp,'ki=',ki
      move_forward_P = .true.
        dP(k) = Piki - Ppkp
        Ks(k) = ki-1
    elseif (Pp(kp) > Pi(ki)) then ! case L-P
                                               if (debug) write(*,'(4x,2(x,a,i3))') 'L-P kp=',kp,'ki=',ki
      move_forward_L = .true.
      dP(k) = Ppkp - Piki
      Ks(k) = ki
    elseif (Pp(kp) == Pi(ki)) then ! best choice depends on next interfaces
      if (Pp(kpp1) <= Pi(kip1)) then ! case L-P
                                               if (debug) write(*,'(4x,x,a,f10.3)') 'P==L L-P'
        move_forward_L = .true.
      else ! Pp(kpp1) >= Pi(kip1)  case P-L
                                               if (debug) write(*,'(4x,x,a,f10.3)') 'P==L P-L'
        move_forward_P = .true.
      endif
      dP(k) = 0.
    else
      stop 'I am thinking about how this might happen'
    endif
    dP(k) = max(0., dP(k))
    Ks(k) = max(1, Ks(k))
    Ks(k) = min(nk, Ks(k))
                                               if (debug) write(*,'(4x,3(x,a,f10.3))') 'dP(k)=',dP(k)

    if (move_forward_P) then
      if (kp < nk+1) then
        kp = kp + 1
      else ! No choice but to move ki
        ki = ki + 1
      endif
    elseif (move_forward_L) then
      if (ki < nk+1) then
        ki = ki + 1
      else ! No choice but to move kp
        kp = kp + 1
      endif
    else
      stop 'Neither L or P moved forward'
    endif

  enddo

end subroutine project_thicknesses

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

subroutine neutral_diffusion_end(CS)
  type(neutral_diffusion_CS), pointer :: CS

  if (associated(CS)) deallocate(CS)

end subroutine neutral_diffusion_end

!> Returns true if unit tests of neutral_diffusion functions fail. Otherwise returns false.
logical function neutralDiffusionUnitTests()
  integer, parameter :: nk = 4
  real, dimension(nk+1) :: zl, zr1, zr2, zr3 ! Test interface positions
  real, dimension(2*nk+1) :: dP, dP1, dP2 ! Test thicknesses
  integer, dimension(2*nk+1) :: Ko, Kl1, Kr1, Kl2, Kr2 ! Test indexes
  real, dimension(nk+1) :: rhol, rhor1, rhor3 ! Test interface densities
  real, dimension(nk+1) :: drdpl ! Test interface compressibility
  real, dimension(nk+1) :: Po ! Test positions
  data zl / 0., 10., 20., 30., 40. /
  data zr1 / -2., -1., 4., 42., 43. /
  data zr2 / 0., 0., 20., 20., 25. /
  data dP1 / 0., 0., 4., 6., 10., 10., 10., 0., 0. /
  data Kl1 / 1, 1, 1, 1, 2, 3, 4, 4, 4 /
  data Kr1 / 1, 2, 2, 3, 3, 3, 3, 3, 4 /
  data dP2 / 0., 0., 10., 10., 0., 0., 5., 0., 0. /
  data Kl2 / 1, 1, 1, 2, 3, 3, 3, 3, 4 /
  data Kr2 / 1, 1, 2, 2, 3, 3, 4, 4, 4 /
  data rhol / 1., 1.5, 2.5, 3.5, 4. /
  data drdpl / 0.1, 0.1, 0.1, 0.1, 0.1 /
  data rhor1 / 0.5, 1.25, 3.0, 3.5, 5.0 /
  data zr3 / 0., 0., 40., 40., 40. /
  data rhor3 / 1., 1.5, 2.5, 3.5, 4. /
  integer :: k

  neutralDiffusionUnitTests = .false. ! Normally return false
  write(*,*) '===== MOM_remapping: neutralDiffusionUnitTests =================='
  write(*,*) 'Project left->right 1'
  call project_thicknesses(nk, zl, zr1, dP, Ko)
  do k = 1,2*nk+1
    if (dP(k) /= dP1(k)) neutralDiffusionUnitTests = .true.
    if (Ko(k) /= Kl1(k)) neutralDiffusionUnitTests = .true.
    write(*,*) k, dP(k), Ko(k), .not.neutralDiffusionUnitTests
  enddo
  write(*,*) 'Project right->left 1'
  call project_thicknesses(nk, zr1, zl, dP, Ko)
  do k = 1,2*nk+1
    if (dP(k) /= dP1(k)) neutralDiffusionUnitTests = .true.
    if (Ko(k) /= Kr1(k)) neutralDiffusionUnitTests = .true.
    write(*,*) k, dP(k), Ko(k), .not.neutralDiffusionUnitTests
  enddo
  write(*,*) 'Project left->right 2'
  call project_thicknesses(nk, zl, zr2, dP, Ko)
  do k = 1,2*nk+1
    if (dP(k) /= dP2(k)) neutralDiffusionUnitTests = .true.
    if (Ko(k) /= Kl2(k)) neutralDiffusionUnitTests = .true.
    write(*,*) k, dP(k), Ko(k), .not.neutralDiffusionUnitTests
  enddo
  write(*,*) 'Project right->left 2'
  call project_thicknesses(nk, zr2, zl, dP, Ko)
  do k = 1,2*nk+1
    if (dP(k) /= dP2(k)) neutralDiffusionUnitTests = .true.
    if (Ko(k) /= Kr2(k)) neutralDiffusionUnitTests = .true.
    write(*,*) k, dP(k), Ko(k), .not.neutralDiffusionUnitTests
  enddo
  write(*,*) 'Find positions right->left, 0.0 compressibility'
  call projected_positions(nk, rhor1, zr1, rhol, zl, 0.*drdpl, Po)
  do k = 1,nk+1
    write(*,*) k, Po(k), .not.neutralDiffusionUnitTests
  enddo
  write(*,*) 'Find positions left->right, 0.0 compressibility 3'
  call projected_positions(nk, rhol, zl, rhor3, zr3, 0.*drdpl, Po)
  do k = 1,nk+1
    write(*,*) k, Po(k), .not.neutralDiffusionUnitTests
  enddo
  write(*,*) 'Find positions right->left, 0.0 compressibility 3'
  call projected_positions(nk, rhor3, zr3, rhol, zl, 0.*drdpl, Po)
  do k = 1,nk+1
    write(*,*) k, Po(k), .not.neutralDiffusionUnitTests
  enddo
  write(*,*) 'Find positions right->left, 0.25 compressibility'
  call projected_positions(nk, rhor1, zr1, rhol, zl, 0.25*drdpl, Po)
  do k = 1,nk+1
    write(*,*) k, Po(k), .not.neutralDiffusionUnitTests
  enddo
  write(*,*) 'Find positions left->right, 0.25 compressibility 3'
  call projected_positions(nk, rhol, zl, rhor3, zr3, 0.25*drdpl, Po)
  do k = 1,nk+1
    write(*,*) k, Po(k), .not.neutralDiffusionUnitTests
  enddo
  write(*,*) 'Find positions right->left, 0.25 compressibility 3'
  call projected_positions(nk, rhor3, zr3, rhol, zl, 0.25*drdpl, Po)
  do k = 1,nk+1
    write(*,*) k, Po(k), .not.neutralDiffusionUnitTests
  enddo

  write(*,*) '=========================================================='

end function neutralDiffusionUnitTests


end module MOM_neutral_diffusion
