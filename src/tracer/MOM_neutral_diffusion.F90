!> A column-wise toolbox for implementing neutral (horizontal) diffusion
module MOM_neutral_diffusion

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_EOS, only : EOS_type, calculate_compress, calculate_density_derivs
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_error_handler, only : MOM_get_verbosity
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_file_parser, only : openParameterBlock, closeParameterBlock
use MOM_grid, only : ocean_grid_type

implicit none ; private

#include <MOM_memory.h>

public neutral_diffusion, neutral_diffusion_init, neutral_diffusion_end
public neutral_diffusion_calc_coeffs
public neutralDiffusionUnitTests

type, public :: neutral_diffusion_CS ; private
  integer :: nkp1   ! Number of interfaces for a column = nk + 1
  integer :: nkp1X2 ! Number of intersecting interfaces between columns = 2 * nkp1

  real,    allocatable, dimension(:,:,:) :: uPoL ! Non-dimensional position with left layer uKoL-1, u-point
  real,    allocatable, dimension(:,:,:) :: uPoR ! Non-dimensional position with right layer uKoR-1, u-point
  integer, allocatable, dimension(:,:,:) :: uKoL ! Index of left interface corresponding to neutral surface, u-point
  integer, allocatable, dimension(:,:,:) :: uKoR ! Index of right interface corresponding to neutral surface, u-point
  real,    allocatable, dimension(:,:,:) :: uHeff ! Effective thickness at u-point (H units)
  real,    allocatable, dimension(:,:,:) :: vPoL ! Non-dimensional position with left layer uKoL-1, v-point
  real,    allocatable, dimension(:,:,:) :: vPoR ! Non-dimensional position with right layer uKoR-1, v-point
  integer, allocatable, dimension(:,:,:) :: vKoL ! Index of left interface corresponding to neutral surface, v-point
  integer, allocatable, dimension(:,:,:) :: vKoR ! Index of right interface corresponding to neutral surface, v-point
  real,    allocatable, dimension(:,:,:) :: vHeff ! Effective thickness at v-point (H units)

end type neutral_diffusion_CS

! This include declares and sets the variable "version".
#include "version_variable.h"
character(len=40)  :: mod = "MOM_neutral_diffusion" ! This module's name.

logical, parameter :: debug_this_module = .true.

contains

!> Read parameters and allocates control structure for neutral_diffusion module.
logical function neutral_diffusion_init(Time, G, param_file, diag, CS)
  type(time_type), target,    intent(in)    :: Time       !< Time structure
  type(ocean_grid_type),      intent(in)    :: G          !< Grid structure
  type(diag_ctrl), target,    intent(inout) :: diag       !< Diagnostics control structure
  type(param_file_type),      intent(in)    :: param_file !< Parameter file structure
  type(neutral_diffusion_CS), pointer       :: CS         !< Neutral diffusion control structure
  ! Local variables
  character(len=256) :: mesg    ! Message for error messages.

  neutral_diffusion_init = .false.
  if (associated(CS)) then
    call MOM_error(FATAL, "neutral_diffusion_init called with associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, &
       "This module implements neutral diffusion of tracers")
  call get_param(param_file, mod, "USE_NEUTRAL_DIFFUSION", neutral_diffusion_init, &
                 "If true, enables the neutral diffusion module.", &
                 default=.false.)
  call openParameterBlock(param_file,'NEUTRAL_DIFF')
! call get_param(param_file, mod, "KHTR", CS%KhTr, &
!                "The background along-isopycnal tracer diffusivity.", &
!                units="m2 s-1", default=0.0)
  call closeParameterBlock(param_file)

  if (.not.neutral_diffusion_init) then
    deallocate(CS)
    return
  endif
  
  ! U-points
  allocate(CS%uPoL(G%isd:G%ied,G%jsd:G%jed,2*G%ke+2)); CS%uPoL(G%isc-1:G%iec,G%jsc:G%jec,:) = 0.
  allocate(CS%uPoR(G%isd:G%ied,G%jsd:G%jed,2*G%ke+2)); CS%uPoR(G%isc-1:G%iec,G%jsc:G%jec,:) = 0.
  allocate(CS%uKoL(G%isd:G%ied,G%jsd:G%jed,2*G%ke+2)); CS%uKoL(G%isc-1:G%iec,G%jsc:G%jec,:) = 0
  allocate(CS%uKoR(G%isd:G%ied,G%jsd:G%jed,2*G%ke+2)); CS%uKoR(G%isc-1:G%iec,G%jsc:G%jec,:) = 0
  allocate(CS%uHeff(G%isd:G%ied,G%jsd:G%jed,2*G%ke+1)); CS%uHeff(G%isc-1:G%iec,G%jsc:G%jec,:) = 0
  ! V-points
  allocate(CS%vPoL(G%isd:G%ied,G%jsd:G%jed,2*G%ke+2)); CS%vPoL(G%isc:G%iec,G%jsc-1:G%jec,:) = 0.
  allocate(CS%vPoR(G%isd:G%ied,G%jsd:G%jed,2*G%ke+2)); CS%vPoR(G%isc:G%iec,G%jsc-1:G%jec,:) = 0.
  allocate(CS%vKoL(G%isd:G%ied,G%jsd:G%jed,2*G%ke+2)); CS%vKoL(G%isc:G%iec,G%jsc-1:G%jec,:) = 0
  allocate(CS%vKoR(G%isd:G%ied,G%jsd:G%jed,2*G%ke+2)); CS%vKoR(G%isc:G%iec,G%jsc-1:G%jec,:) = 0
  allocate(CS%vHeff(G%isd:G%ied,G%jsd:G%jed,2*G%ke+1)); CS%vHeff(G%isc:G%iec,G%jsc-1:G%jec,:) = 0

end function neutral_diffusion_init

!> Calculates remapping factors for u/v columns used to map adjoining columns to
!! a shared coordinate space.
subroutine neutral_diffusion_calc_coeffs(G, h, T, S, EOS, CS)
  type(ocean_grid_type),                 intent(in) :: G !< Ocean grid structure
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in) :: h !< Layer thickness (H units, m or Pa)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in) :: T !< Potential temperature (degC)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in) :: S !< Salinity (ppt)
  type(EOS_type),                        pointer    :: EOS !< Equation of state structure
  type(neutral_diffusion_CS),            pointer    :: CS !< Neutral diffusion constrol structure
  ! Local variables
  integer :: i, j, k
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: Tint ! Interface T (degC)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: Sint ! Interface S (ppt)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: Pint ! Interface pressure (Pa)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: dRdT ! Interface thermal expansion coefficient (kg/m3/degC)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: dRdS ! Interface haline expansion coefficient (kg/m3/ppt)

  do j = G%jsc-1, G%jec+1
    ! Interpolate state to interface
    do i = G%isc-1, G%iec+1
      call interface_TS(G%ke, T(i,j,:), S(i,j,:), Tint, Sint)
    enddo

    ! Caclulate interface properties
    Pint(:,j,1) = 0. ! Assume P=0 (Pa) at surface - needs correcting for atmospheric and ice loading - AJA
    do k = 1, G%ke+1
      call calculate_density_derivs(Tint(:,j,k), Sint(:,j,k), Pint(:,j,k), &
                                    dRdT(:,j,k), dRdS(:,j,k), G%isc-1, G%iec-G%isc+3, EOS)
      if (k<=G%ke) Pint(:,j,k+1) = Pint(:,j,k) + h(:,j,k) * G%H_to_Pa ! Pressure at next interface, k+1 (Pa)
    enddo
  enddo

  ! Neutral surface factors at U points
  do j = G%jsc, G%jec
    do I = G%isc-1, G%iec
      call find_neutral_surface_positions(G%ke, &
               Pint(i,j,k), Tint(i,j,k), Sint(i,j,k), dRdT(i,j,k), dRdS(i,j,k), &
               Pint(i+1,j,k), Tint(i+1,j,k), Sint(i+1,j,k), dRdT(i+1,j,k), dRdS(i+1,j,k), &
               CS%uPoL(I,j,k), CS%uPoR(I,j,k), CS%uKoL(I,j,k), CS%uKoR(I,j,k), CS%uhEff(I,j,k) )
    enddo
  enddo

  ! Neutral surface factors at V points
  do J = G%jsc-1, G%jec
    do i = G%isc, G%iec
      call find_neutral_surface_positions(G%ke, &
               Pint(i,j,k), Tint(i,j,k), Sint(i,j,k), dRdT(i,j,k), dRdS(i,j,k), &
               Pint(i,j+1,k), Tint(i,j+1,k), Sint(i,j+1,k), dRdT(i,j+1,k), dRdS(i,j+1,k), &
               CS%vPoL(i,J,k), CS%vPoR(i,J,k), CS%vKoL(i,J,k), CS%vKoR(i,J,k), CS%vhEff(i,J,k) )
    enddo
  enddo

end subroutine neutral_diffusion_calc_coeffs

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

  !call find_neutral_surface_position(nk, Til, Sil, dRdTl, dRdSl, Tir, Sir, Pir, PilR)

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
  ! NOTE: THIS IS A PLACEHOLDER FOR HIGHER ORDER INTERPOLATION T.B.I.
  Ti(1) = T(1) ; Si(1) = S(1)
  do k = 2, nk
    Ti(k) = 0.5*( T(k-1) + T(k) )
    Si(k) = 0.5*( S(k-1) + S(k) )
  enddo
  Ti(nk+1) = T(nk) ; Si(nk+1) = S(nk)

end subroutine interface_TS

!> Returns positions within left/right columns of combined interfaces
subroutine find_neutral_surface_positions(nk, Pl, Tl, Sl, dRdTl, dRdSl, Pr, Tr, Sr, dRdTr, dRdSr, PoL, PoR, KoL, KoR, hEff)
  integer,                    intent(in)    :: nk    !< Number of levels
  real, dimension(nk+1),      intent(in)    :: Pl    !< Left-column interface pressure (Pa)
  real, dimension(nk+1),      intent(in)    :: Tl    !< Left-column interface potential temperature (degC)
  real, dimension(nk+1),      intent(in)    :: Sl    !< Left-column interface salinity (ppt)
  real, dimension(nk+1),      intent(in)    :: dRdTl !< Left-column dRho/dT (kg/m3/degC)
  real, dimension(nk+1),      intent(in)    :: dRdSl !< Left-column dRho/dS (kg/m3/ppt)
  real, dimension(nk+1),      intent(in)    :: Pr    !< Right-column interface pressure (Pa)
  real, dimension(nk+1),      intent(in)    :: Tr    !< Right-column interface potential temperature (degC)
  real, dimension(nk+1),      intent(in)    :: Sr    !< Right-column interface salinity (ppt)
  real, dimension(nk+1),      intent(in)    :: dRdTr !< Left-column dRho/dT (kg/m3/degC)
  real, dimension(nk+1),      intent(in)    :: dRdSr !< Left-column dRho/dS (kg/m3/ppt)
  real, dimension(2*nk+2),    intent(inout) :: PoL   !< Position of neutral surface in left column (Pa)
  real, dimension(2*nk+2),    intent(inout) :: PoR   !< Position of neutral surface in right column (Pa)
  integer, dimension(2*nk+2), intent(inout) :: KoL   !< Index of first left interface below neutral surface.
  integer, dimension(2*nk+2), intent(inout) :: KoR   !< Index of first right interface below neutral surface.
  real, dimension(2*nk+1),    intent(inout) :: hEff  !< Effective thickness between two neutral surfaces (Pa).
  ! Local variables
  integer :: k_surface ! Index of neutral surface
  integer :: kl ! Index of left interface
  integer :: kr ! Index of right interface
  real :: dRdT, dRdS ! dRho/dT and dRho/dS for the neutral surface
  logical :: looking_left ! True if searching for the position of a right interface in the left column
  logical :: looking_right ! True if searching for the position of a left interface in the right column
  integer :: krm1, klm1
  real :: dRho, dRhoM1, hL, hR

  ! Initialize variables for the search
  kr = 1
  kl = 1

  ! Loop over each neutral surface, working from top to bottom
  neutral_surfaces: do k_surface = 1, 2*nk+2
    klm1 = max(kl-1, 1)
    krm1 = max(kr-1, 1)

    ! Potential density difference, kr - kl 
    dRho = 0.5 * ( ( dRdTr(kr) + dRdTl(kl) ) * ( Tr(kr) - Tl(kl) ) &
                 + ( dRdSr(kr) + dRdSl(kl) ) * ( Sr(kr) - Sl(kl) ) )
    ! Which column has the lighter surface for the current indexes, kr and kl
    if (k_surface<2*nk+2) then
      if (dRho < 0.) then
        looking_left = .true.
        looking_right = .false.
      elseif (dRho >= 0.) then
        looking_right = .true.
        looking_left = .false.
      endif
    else
      ! Note from AJA: This handles the last neutral surface that always sees the same values of kr and kl as
      ! the previous neutral surface. I do not like this kind of logic. Is there a better way to construct kr
      ! and kl from k_surface?
      looking_left = .not. looking_left
      looking_right = .not. looking_right
    endif
                                                     if (debug_this_module) write(0,*) 'k,kr,kl=',k_surface,kr,kl
 
    if (looking_left) then
      ! Interpolate for the neutral surface position within the left column, layer kl
                                                     if (debug_this_module) write(0,*) 'looking_left=',looking_left
      ! Potential density difference, rho(kl-1) - rho(kr) (should be negative)
      dRhoM1 = 0.5 * ( ( dRdTl(klm1) + dRdTr(kr) ) * ( Tl(klm1) - Tr(kr) ) &
                     + ( dRdSl(klm1) + dRdSr(kr) ) * ( Sl(klm1) - Sr(kr) ) )
      ! Potential density difference, rho(kl) - rho(kr) (will be positive)
      !dRho = 0.5 * ( ( dRdTl(kl) + dRdTr(kr) ) * ( Tl(kl) - Tr(kr) ) &
      !             + ( dRdSl(kl) + dRdSr(kr) ) * ( Sl(kl) - Sr(kr) ) )
      dRho = - dRho ! Re-use calculation above

      ! Because we a looking left, the right surface, kr, is lighter than kl and should be denser than kl-1
      ! unless we are still at the top of the left column (kl=1)
      if (dRhoM1 > 0.) then
       !if (kl>1) stop 'This should never happen: kl>1 and dRhoM1>=0.'
        PoL(k_surface) = 1.
      else
        ! Linearly interpolate for the position between Pl(kl-1) and Pl(kl) where the density difference
        ! between right and left is zero.
        PoL(k_surface) = interpolate_for_nondim_position( dRhoM1, Pl(klm1), dRho, Pl(kl) )
      endif
      PoR(k_surface) = 1.
      KoR(k_surface) = kr
      KoL(k_surface) = kl
                                                     if (debug_this_module) write(0,*) '  dRhoM1=',dRhoM1,' dRho=',dRho
                                                     if (debug_this_module) write(0,*) '  PoL(k)=',PoL(k_surface)
      if (kr <= nk) then
        kr = kr + 1
      else
        kl = min(kl + 1, nk+1)
      endif
                                                     if (debug_this_module) write(0,*) '  kr=',kr,' kl=',kl
    elseif (looking_right) then
      ! Interpolate for the neutral surface position within the right column, layer kr
                                                     if (debug_this_module) write(0,*) 'looking_right=',looking_right
      ! Potential density difference, rho(kr-1) - rho(kl) (should be negative)
      dRhoM1 = 0.5 * ( ( dRdTr(krm1) + dRdTl(kl) ) * ( Tr(krm1) - Tl(kl) ) &
                     + ( dRdSr(krm1) + dRdSl(kl) ) * ( Sr(krm1) - Sl(kl) ) )
      ! Potential density difference, rho(kr) - rho(kl) (will be positive)
      !dRho = 0.5 * ( ( dRdTr(kr) + dRdTl(kl) ) * ( Tr(kr) - Tl(kl) ) &
      !             + ( dRdSr(kr) + dRdSl(kl) ) * ( Sr(kr) - Sl(kl) ) )
      dRho = dRho ! Re-use calculation above

      ! Because we a looking right, the left surface, kl, is lighter than kr and should be denser than kr-1
      ! unless we are still at the top of the right column (kr=1)
      if (dRhoM1 > 0.) then
       !if (kr>1) stop 'This should never happen: kr>1 and dRhoM1>=0.'
        PoR(k_surface) = 1.
      else
        ! Linearly interpolate for the position between Pr(kr-1) and Pr(kr) where the density difference
        ! between right and left is zero.
        PoR(k_surface) = interpolate_for_nondim_position( dRhoM1, Pr(krm1), dRho, Pr(kr) )
      endif
      PoL(k_surface) = 1.
      KoL(k_surface) = kl
      KoR(k_surface) = kr
                                                     if (debug_this_module) write(0,*) '  dRhoM1=',dRhoM1,' dRho=',dRho
                                                     if (debug_this_module) write(0,*) '  PoR(k)=',PoR(k_surface)
      if (kl <= nk) then
        kl = kl + 1
      else
        kr = min(kr + 1, nk+1)
      endif
                                                     if (debug_this_module) write(0,*) '  kl=',kl,' kr=',kr
    else
      stop 'Else what?'
    endif

    ! Effective thickness
    ! NOTE: This would be better expressed in terms of the layers thicknesses rather
    ! than as differences of position - AJA
    if (k_surface>1) then
      hL = absolute_position(nk,Pl,KoL,PoL,k_surface) - absolute_position(nk,Pl,KoL,PoL,k_surface-1)
      hR = absolute_position(nk,Pr,KoR,PoR,k_surface) - absolute_position(nk,Pr,KoR,PoR,k_surface-1)
      if ( hL + hR > 0.) then
        hEff(k_surface-1) = 2. * hL * hR / ( hL + hR )
      else
        hEff(k_surface-1) = 0.
      endif
                                                     if (debug_this_module) write(0,*) '  hEff=',hEff(k_surface-1)
    endif

  enddo neutral_surfaces

  contains

  !> Converts non-dimensional positions within a layer to absolute positions (for debugging)
  real function absolute_position(n,Pint,Karr,NParr,k_surface)
    integer, intent(in) :: n            !< Number of levels
    real,    intent(in) :: Pint(n+1)    !< Position of interface (Pa)
    integer, intent(in) :: Karr(2*n+2)  !< Index of deeper 
    real,    intent(in) :: NParr(2*n+2) !< Non-dimensional position with layer Karr(:)-1
    ! Local variables
    integer :: k_surface, k, km1
  
    k = Karr(k_surface)
    km1 = max(1, k-1)
    absolute_position = Pint(km1) + NParr(k_surface) * ( Pint(k) - Pint(km1) )

  end function absolute_position

end subroutine find_neutral_surface_positions

!> Converts non-dimensional positions within a layer to absolute positions (for debugging)
function absolute_positions(n,Pint,Karr,NParr)
  integer, intent(in) :: n            !< Number of levels
  real,    intent(in) :: Pint(n+1)    !< Position of interface (Pa)
  integer, intent(in) :: Karr(2*n+2)  !< Index of deeper 
  real,    intent(in) :: NParr(2*n+2) !< Non-dimensional position with layer Karr(:)-1
  real,  dimension(2*n+2) :: absolute_positions ! Absolute positions (Pa)
  ! Local variables
  integer :: k_surface, k, km1

  do k_surface = 1, 2*n+2
    k = Karr(k_surface)
    km1 = max(1, k-1)
    absolute_positions(k_surface) = Pint(km1) + NParr(k_surface) * ( Pint(k) - Pint(km1) )
  enddo

end function absolute_positions

!> Returns the position between Pneg and Ppos where the interpolated density difference equals
!! zero: Pint = ( dRhoPos * Ppos - dRhoNeg * Pneg ) / ( dRhoPos - dRhoneg )
!! The result is always bounded to be between Pneg and Ppos.
real function interpolate_for_position(dRhoNeg, Pneg, dRhoPos, Ppos)
  real, intent(in) :: dRhoNeg !< Negative density difference
  real, intent(in) :: Pneg    !< Position of negative density difference
  real, intent(in) :: dRhoPos !< Positive density difference
  real, intent(in) :: Ppos    !< Position of positive density difference
  ! Local variables
  real :: wghtU, wghtD, Pint

  if (Ppos<Pneg) stop 'interpolate_for_position: Houston, we have a problem! Ppos<Pneg'
  if (dRhoNeg>dRhoPos) stop 'interpolate_for_position: Houston, we have a problem! dRhoNeg>dRhoPos'
  if (Ppos<=Pneg) then ! Handle vanished or inverted layers
    wghtU = 0.5
    wghtD = 0.5
    Pint = 0.5 * ( Pneg + Ppos )
  elseif ( dRhoPos - dRhoNeg > 0. ) then
    wghtU = -dRhoNeg / ( dRhoPos - dRhoNeg )
    wghtD = dRhoPos / ( dRhoPos - dRhoNeg )
    if ( wghtU < 0.5 ) then
      Pint = Pneg + max( wghtU, 0. ) * ( Ppos - Pneg )
    elseif ( wghtD < 0.5 ) then
      Pint = Ppos + max( wghtD, 0. ) * ( Pneg - Ppos )
    else
      Pint = 0.5 * ( Pneg + Ppos )
    endif
  elseif ( dRhoPos - dRhoNeg == 0) then
    if (dRhoNeg>0.) then
      wghtU = 0.
      wghtD = 1.
      Pint = Pneg
    elseif (dRhoNeg<0.) then
      wghtU = 1.
      wghtD = 0.
      Pint = Ppos
    else ! dRhoPos = dRhoNeg = 0
      wghtU = 0.5
      wghtD = 0.5
      Pint = 0.5 * ( Pneg + Ppos )
    endif
  else ! dRho - dRhoNeg < 0
    wghtU = 0.5
    wghtD = 0.5
    Pint = 0.5 * ( Pneg + Ppos )
  endif
  if ( Pint < Pneg ) stop 'interpolate_for_position: Houston, we have a problem! Pint < Pneg'
  if ( Pint > Ppos ) stop 'interpolate_for_position: Houston, we have a problem! Pint > Ppos'
  interpolate_for_position = Pint
end function interpolate_for_position

!> Returns the non-dimensnional position between Pneg and Ppos where the interpolated density difference equals zero.
!! The result is always bounded to be between 0 and 1.
real function interpolate_for_nondim_position(dRhoNeg, Pneg, dRhoPos, Ppos)
  real, intent(in) :: dRhoNeg !< Negative density difference
  real, intent(in) :: Pneg    !< Position of negative density difference
  real, intent(in) :: dRhoPos !< Positive density difference
  real, intent(in) :: Ppos    !< Position of positive density difference

  if (Ppos<Pneg) stop 'interpolate_for_position: Houston, we have a problem! Ppos<Pneg'
  if (dRhoNeg>dRhoPos) stop 'interpolate_for_position: Houston, we have a problem! dRhoNeg>dRhoPos'
  if (Ppos<=Pneg) then ! Handle vanished or inverted layers
    interpolate_for_nondim_position = 0.5
  elseif ( dRhoPos - dRhoNeg > 0. ) then
    interpolate_for_nondim_position = min( 1., max( 0., -dRhoNeg / ( dRhoPos - dRhoNeg ) ) )
  elseif ( dRhoPos - dRhoNeg == 0) then
    if (dRhoNeg>0.) then
      interpolate_for_nondim_position = 0.
    elseif (dRhoNeg<0.) then
      interpolate_for_nondim_position = 1.
    else ! dRhoPos = dRhoNeg = 0
      interpolate_for_nondim_position = 0.5
    endif
  else ! dRho - dRhoNeg < 0
    interpolate_for_nondim_position = 0.5
  endif
  if ( interpolate_for_nondim_position < 0. ) stop 'interpolate_for_position: Houston, we have a problem! Pint < Pneg'
  if ( interpolate_for_nondim_position > 1. ) stop 'interpolate_for_position: Houston, we have a problem! Pint > Ppos'
end function interpolate_for_nondim_position

!> Returns true if unit tests of neutral_diffusion functions fail. Otherwise returns false.
logical function neutralDiffusionUnitTests()
  integer, parameter :: nk = 4
  real, dimension(nk+1)   :: TiL, TiR1, TiR2, TiR4 ! Test interface temperatures
  real, dimension(nk+1)   :: SiL ! Test interface salinities
  real, dimension(nk+1)   :: PiL, PiR4 ! Test interface positions
  real, dimension(nk+1)   :: dRdT, dRdS ! Test interface expansion coefficients
  real, dimension(2*nk+2) :: PiLRo, PiRLo ! Test positions
  integer, dimension(2*nk+2) :: KoL, KoR ! Test indexes
  real, dimension(2*nk+1) :: hEff ! Test positions
  real, dimension(2*nk+2) :: pL0, pL1, pL2, pL3, pL4 ! Test positions
  real, dimension(2*nk+2) :: pR0, pR1, pR2, pR3, pR4 ! Test positions
  real, dimension(2*nk+1) :: hE0, hE1, hE2, hE3, hE4 ! Test positions
  integer, dimension(2*nk+2) :: kL0, kL1, kL2 ! Test indexes
  integer, dimension(2*nk+2) :: kR0, kR1, kR2 ! Test indexes
  ! Fixed left column values
  data PiL / 0., 10., 20., 30., 40. /
  data TiL / 10., 7.5, 5., 2.5, 0. /
  data SiL / 0., 0., 0., 0., 0. /
  data dRdT / -0.25, -0.25, -0.25, -0.25, -0.25 /
  data dRdS / 0.5, 0.5, 0.5, 0.5, 0.5 /
  ! Identical columns, answers
  data pL0 / 0., 0., 10., 10., 20., 20., 30., 30., 40., 40. /
  data pR0 / 0., 0., 10., 10., 20., 20., 30., 30., 40., 40. /
  data hE0 / 0., 10., 0., 10., 0., 10., 0., 10., 0. /
  data kL0 / 1, 1, 1, 2, 2, 3, 3, 4, 4, 5 /
  data kR0 / 1, 1, 2, 2, 3, 3, 4, 4, 5, 5 /
  ! Slightly warmer on right, answers
  data pL1 / 0., 0., 2., 10., 12., 20., 22., 30., 32., 40. /
  data pR1 / 0., 8., 10., 18., 20., 28., 30., 38., 40., 40. /
  data hE1 / 0., 2., 8., 2., 8., 2., 8., 2., 0. /
  data kL1 / 1, 1, 1, 2, 2, 3, 3, 4, 4, 5 /
  data kR1 / 1, 1, 2, 2, 3, 3, 4, 4, 5, 5 /
  ! Warmer on right, answers
  data pL2 / 0., 0., 0., 0., 0., 8., 20., 30., 40., 40. /
  data pR2 / 0., 10., 20., 30., 32., 40., 40., 40., 40., 40. /
  data hE2 / 0., 0., 0., 0., 8., 0., 0., 0., 0. /
  ! Strong stratification on right, input
  data TiR1 / 20., 10., 0., -7.5, -10. /
  ! Strong stratification on right, answers
  data pL3 / 0., 0., 0., 10., 20., 30., 40., 40., 40., 40. /
  data pR3 / 0., 10., 10., 12.5, 15., 17.5, 20., 20., 40., 40. /
  data hE3 / 0., 0., 4., 4., 4., 4., 0., 0., 0. /
  ! Weak stratification on right, input
  data TiR2 / 7.5, 6.875, 6.25, 5.625, 5.0 /
  ! Vanished layers on right, input
  data TiR4 / 10., 5., 5., 5., 0. /
  data PiR4 / 0., 20., 20., 20., 40. /
  ! Vanished layers on right, answers
  data pL4 / 0., 0., 10., 20., 20., 20., 20., 30., 40., 40. /
  data pR4 / 0., 0., 10., 20., 20., 20., 20., 30., 40., 40. /
  data hE4 / 0., 10., 10., 0., 0., 0., 10., 10., 0. /

  integer :: k
  logical :: verbosity

  verbosity = MOM_get_verbosity()

  neutralDiffusionUnitTests = .false. ! Normally return false
  write(*,'(a)') '===== MOM_neutral_diffusion: neutralDiffusionUnitTests =================='

  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifp(-1.0, 0.,  1.0, 1.0, 0.5, 'Check mid-point')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifp( 0.0, 0.,  1.0, 1.0, 0.0, 'Check bottom')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifp( 0.1, 0.,  1.1, 1.0, 0.0, 'Check below')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifp(-1.0, 0.,  0.0, 1.0, 1.0, 'Check top')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifp(-1.0, 0., -0.1, 1.0, 1.0, 'Check above')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifp(-1.0, 0.,  3.0, 1.0, 0.25, 'Check 1/4')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifp(-3.0, 0.,  1.0, 1.0, 0.75, 'Check 3/4')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifp( 1.0, 0.,  1.0, 1.0, 0.0, 'Check dRho=0 below')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifp(-1.0, 0., -1.0, 1.0, 1.0, 'Check dRho=0 above')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifp( 0.0, 0.,  0.0, 1.0, 0.5, 'Check dRho=0 mid')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifp(-2.0, .5,  5.0, 0.5, 0.5, 'Check dP=0')

  call find_neutral_surface_positions(nk, PiL, TiL, SiL, dRdt, dRdS, PiL, TiL, SiL, dRdT, dRdS, PiLRo, PiRLo, KoL, KoR, hEff)
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL, KoL, PiLRo), pL0, 'Identical columns, left positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL, KoR, PiRLo), pR0, 'Identical columns, right positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+1, hEff, hE0, 'Identical columns, thicknesses')
stop

  call find_neutral_surface_positions(nk, PiL, TiL, SiL, dRdt, dRdS, PiL+2., TiL, SiL, dRdT, dRdS, PiLRo, PiRLo, KoL, KoR, hEff)
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL, KoL, PiLRo), pL0, 'Same values raised on right, left positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL+2., KoR, PiRLo), pR0+2., 'Same values raised on right, right positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+1, hEff, hE0, 'Same values raised on right, thicknesses')

  call find_neutral_surface_positions(nk, PiL, TiL, SiL, dRdt, dRdS, PiL-2., TiL, SiL, dRdT, dRdS, PiLRo, PiRLo, KoL, KoR, hEff)
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL, KoL, PiLRo), pL0, 'Same values lowered on right, left positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL-2., KoR, PiRLo), pR0-2., 'Same values lowered on right, right positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+1, hEff, hE0, 'Same values lowered on right, thicknesses')

  call find_neutral_surface_positions(nk, PiL, TiL, SiL, dRdt, dRdS, PiL, TiL+2., SiL, dRdT, dRdS, PiLRo, PiRLo, KoL, KoR, hEff)
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL, KoL, PiLRo), pL1, 'Slightly warmer on right, left positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL, KoR, PiRLo), pR1, 'Slightly warmer on right, right positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+1, hEff, hE1, 'Slightly warmer on right, thicknesses')

  call find_neutral_surface_positions(nk, PiL, TiL, SiL, dRdt, dRdS, PiL, TiL-2., SiL, dRdT, dRdS, PiLRo, PiRLo, KoL, KoR, hEff)
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL, KoL, PiLRo), pR1, 'Slightly cooler on right, left positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL, KoR, PiRLo), pL1, 'Slightly cooler on right, right positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+1, hEff, hE1, 'Slightly cooler on right, thicknesses')

  call find_neutral_surface_positions(nk, PiL, TiL, SiL, dRdt, dRdS, PiL, TiL+8., SiL, dRdT, dRdS, PiLRo, PiRLo, KoL, KoR, hEff)
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL, KoL, PiLRo), pL2, 'Warmer on right, left positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL, KoR, PiRLo), pR2, 'Warmer on right, right positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+1, hEff, hE2, 'Warmer on right, thicknesses')

  call find_neutral_surface_positions(nk, PiL, TiL, SiL, dRdt, dRdS, PiL, TiR1, SiL, dRdT, dRdS, PiLRo, PiRLo, KoL, KoR, hEff)
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL, KoL, PiLRo), pL3, 'Strong stratification on right, left positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL, KoR, PiRLo), pR3, 'Strong stratification on right, right positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+1, hEff, hE3, 'Strong stratification on right, thicknesses')

  call find_neutral_surface_positions(nk, PiL, TiL, SiL, dRdt, dRdS, PiL, TiR2, SiL, dRdT, dRdS, PiLRo, PiRLo, KoL, KoR, hEff)
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL, KoL, PiLRo), pR3, 'Weak stratification on right, left positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL, KoR, PiRLo), pL3, 'Weak stratification on right, right positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+1, hEff, hE3, 'Weak stratification on right, thicknesses')

  call find_neutral_surface_positions(nk, PiL, TiL, SiL, dRdt, dRdS, PiR4, TiR4, SiL, dRdT, dRdS, PiLRo, PiRLo, KoL, KoR, hEff)
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL, KoL, PiLRo), pL4, 'Vanished layers on right, left positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiR4, KoR, PiRLo), pR4, 'Vanished layers on right, right positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+1, hEff, hE4, 'Vanished layers on right, thicknesses')

  call find_neutral_surface_positions(nk, PiR4, TiR4, SiL, dRdt, dRdS, PiL, TiL, SiL, dRdT, dRdS, PiLRo, PiRLo, KoL, KoR, hEff)
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiR4, KoL, PiLRo), pR4, 'Vanished layers on left, left positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL, KoR, PiRLo), pL4, 'Vanished layers on left, right positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+1, hEff, hE4, 'Vanished layers on left, thicknesses')
stop

  write(*,'(a)') '=========================================================='

  contains

  !> Returns true if test of interpolate_for_position() fails, and conditionally writes results to stream
  logical function test_ifp(rhoNeg, Pneg, rhoPos, Ppos, Ptrue, title)
    real,             intent(in) :: rhoNeg !< Lighter density (kg/m3)
    real,             intent(in) :: Pneg   !< Interface position of lighter density (pa)
    real,             intent(in) :: rhoPos !< Heavier density (kg/m3)
    real,             intent(in) :: Ppos   !< Interface position of heavier density (pa)
    real,             intent(in) :: Ptrue  !< True answer (Pa)
    character(len=*), intent(in) :: title !< Title for messages
    ! Local variables
    integer :: stdunit
    real :: Pret

    Pret = interpolate_for_position(rhoNeg, Pneg, rhoPos, Ppos)
    test_ifp = (Pret /= Ptrue) 

    if (test_ifp .or. verbosity>5) then
      stdunit = 6
      if (test_ifp.or.debug_this_module) stdunit = 0 ! In case of wrong results, write to error stream
      write(stdunit,'(a)') title
      if (test_ifp) then
        write(stdunit,'(4(x,a,f20.16),2(x,a,1pe22.15),x,a)') 'r1=',rhoNeg,'p1=',Pneg,'r2=',rhoPos,'p2=',Ppos,'pRet=',Pret,'pTrue=',Ptrue,'WRONG!'
      else
        write(stdunit,'(4(x,a,f20.16),2(x,a,1pe22.15))') 'r1=',rhoNeg,'p1=',Pneg,'r2=',rhoPos,'p2=',Ppos,'pRet=',Pret,'pTrue=',Ptrue
      endif
    endif

  end function test_ifp

  !> Returns true if comparison of Po and Ptrue fails, and conditionally writes results to stream
  logical function test_fnsp(nk, Po, Ptrue, title)
    integer,          intent(in) :: nk !< Number of layers
    real,             intent(in) :: Po(nk) !< Calculated answer
    real,             intent(in) :: Ptrue(nk) !< True answer
    character(len=*), intent(in) :: title !< Title for messages
    ! Local variables
    integer :: k, stdunit

    test_fnsp = .false.
    do k = 1,nk
      if (Po(k) /= Ptrue(k)) test_fnsp = .true.
    enddo

    if (test_fnsp .or. verbosity>5) then
      stdunit = 6
      if (test_fnsp.or.debug_this_module) stdunit = 0 ! In case of wrong results, write to error stream
      write(stdunit,'(a)') title
      do k = 1,nk
        if (Po(k) /= Ptrue(k)) then
          test_fnsp = .true.
          write(stdunit,'(a,i2,2(x,a,f20.16),x,a,1pe22.15,x,a)') 'k=',k,'Po=',Po(k),'Ptrue=',Ptrue(k),'err=',Po(k)-Ptrue(k),'WRONG!'
        else
          if (verbosity>5) &
            write(stdunit,'(a,i2,2(x,a,f20.16),x,a,1pe22.15)') 'k=',k,'Po=',Po(k),'Ptrue=',Ptrue(k),'err=',Po(k)-Ptrue(k)
        endif
      enddo
    endif

  end function test_fnsp

end function neutralDiffusionUnitTests

!> Deallocates neutral_diffusion control structure
subroutine neutral_diffusion_end(CS)
  type(neutral_diffusion_CS), pointer :: CS

  if (associated(CS)) deallocate(CS)

end subroutine neutral_diffusion_end

end module MOM_neutral_diffusion
