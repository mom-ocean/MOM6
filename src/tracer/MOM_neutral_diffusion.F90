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

logical, parameter :: debug_this_module = .false.

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
      call interface_scalar(G%ke, h(i,j,:), T(i,j,:), Tint)
      call interface_scalar(G%ke, h(i,j,:), S(i,j,:), Sint)
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
               Pint(i,j,:), Tint(i,j,:), Sint(i,j,:), dRdT(i,j,:), dRdS(i,j,:), &
               Pint(i+1,j,:), Tint(i+1,j,:), Sint(i+1,j,:), dRdT(i+1,j,:), dRdS(i+1,j,:), &
               CS%uPoL(I,j,:), CS%uPoR(I,j,:), CS%uKoL(I,j,:), CS%uKoR(I,j,:), CS%uhEff(I,j,:) )
    enddo
  enddo

  ! Neutral surface factors at V points
  do J = G%jsc-1, G%jec
    do i = G%isc, G%iec
      call find_neutral_surface_positions(G%ke, &
               Pint(i,j,:), Tint(i,j,:), Sint(i,j,:), dRdT(i,j,:), dRdS(i,j,:), &
               Pint(i,j+1,:), Tint(i,j+1,:), Sint(i,j+1,:), dRdT(i,j+1,:), dRdS(i,j+1,:), &
               CS%vPoL(i,J,:), CS%vPoR(i,J,:), CS%vKoL(i,J,:), CS%vKoR(i,J,:), CS%vhEff(i,J,:) )
    enddo
  enddo

end subroutine neutral_diffusion_calc_coeffs

subroutine neutral_diffusion()
  ! Local variables

end subroutine neutral_diffusion

!> Returns interface scalar, Si, for a column of layer values, S.
subroutine interface_scalar(nk, h, S, Si)
  integer,               intent(in)    :: nk !< Number of levels
  real, dimension(nk),   intent(in)    :: h  !< Layer thickness (H, m or Pa)
  real, dimension(nk),   intent(in)    :: S  !< Layer scalar (conc, e.g. ppt)
  real, dimension(nk+1), intent(inout) :: Si !< Interface scalar (conc, e.g. ppt)
  ! Local variables
  integer :: k
  real, dimension(nk) :: diff
  real :: Sb, Sa

  call PLM_diff(nk, h, S, 2, diff)
  Si(1) = S(1) - 0.5 * diff(1)
  do k = 2, nk
    ! Average of the two edge values (will be bounded and,
    ! when slopes are unlimited, notionally second-order accurate)
    Sa = S(k-1) + 0.5 * diff(k-1) ! Lower edge value of a PLM reconstruction for layer above
    Sb = S(k) - 0.5 * diff(k) ! Upper edge value of a PLM reconstruction for layer below
    Si(k) = 0.5 * ( Sa + Sb )
  enddo
  Si(nk+1) = S(nk) + 0.5 * diff(nk)

end subroutine interface_scalar

!> Returns PLM slopes for a column where the slopes are the difference in value across
!! each cell. Slopes in the first and last cell (top/bottom) are always zero.
subroutine PLM_diff(nk, h, S, c_method, diff)
  integer,             intent(in)    :: nk       !< Number of levels
  real, dimension(nk), intent(in)    :: h        !< Layer thickness (H, m or Pa)
  real, dimension(nk), intent(in)    :: S        !< Layer salinity (conc, e.g. ppt)
  integer,             intent(in)    :: c_method !< Method to use for the centered difference
  real, dimension(nk), intent(inout) :: diff     !< Scalar difference across layer (conc, e.g. ppt)
                                                 !! determined by the following values: 
                                                 !!   1. Second order finite difference (not recommended)
                                                 !!   2. Second order finite volume (used in original PPM)
                                                 !!   3. Finite-volume weighted least squares linear fit
                                                 !! \todo  The use of c_method to choose a scheme is inefficient
                                                 !! and should eventually be moved up the call tree.
  ! Local variables
  integer :: k
  real :: hkm1, hk, hkp1, Skm1, Sk, Skp1, diff_l, diff_r, diff_c

  do k = 2, nk-1
    hkm1 = h(k-1)
    hk = h(k)
    hkp1 = h(k+1)

    if ( ( hkp1 + hk ) * ( hkm1 + hk ) > 0.) then
      Skm1 = S(k-1)
      Sk = S(k)
      Skp1 = S(k+1)
      if (c_method==1) then
        ! Simple centered diff (from White)
        if ( hk + 0.5 * (hkm1 + hkp1) /= 0. ) then
          diff_c = ( Skp1 - Skm1 ) * ( hk / ( hk + 0.5 * (hkm1 + hkp1) ) )
        else
          diff_c = 0.
        endif
      elseif (c_method==2) then
        ! Second order accurate centered FV slope (from Colella and Woodward, JCP 1984)
        diff_c = fv_diff(hkm1, hk, hkp1, Skm1, Sk, Skp1)
      elseif (c_method==3) then
        ! Second order accurate finite-volume least squares slope
        diff_c = hk * fvlsq_slope(hkm1, hk, hkp1, Skm1, Sk, Skp1)
      endif
      ! Limit centered slope by twice the side differenced slopes
      diff_l = 2. * ( Sk - Skm1 )
      diff_r = 2. * ( Skp1 - Sk )
      if (diff_l * diff_r > 0.) then
        diff(k) = sign( min( abs(diff_l), abs(diff_c), abs(diff_r) ), diff_c )
      else
        diff(k) = 0. ! PCM for local extrema
      endif
    else
      diff(k) = 0. ! PCM next to vanished layers
    endif
  enddo
  ! PCM for top and bottom layer
 !diff(1) = 0.
 !diff(nk) = 0.
  ! Linear extrapolation for top and bottom interfaces
  diff(1) = ( S(2) - S(1) ) * 2. * ( h(1) / ( h(1) + h(2) ) )
  diff(nk) = S(nk) - S(nk-1) * 2. * ( h(nk) / ( h(nk-1) + h(nk) ) )

end subroutine PLM_diff

!> Returns the cell-centered second-order finite volume (unlimited PLM) slope
!! using three consecutive cell widths and average values. Slope is returned
!! as a difference across the central cell (i.e. units of scalar S).
real function fv_diff(hkm1, hk, hkp1, Skm1, Sk, Skp1)
  real, intent(in) :: hkm1 !< Left cell width
  real, intent(in) :: hk   !< Center cell width
  real, intent(in) :: hkp1 !< Right cell width
  real, intent(in) :: Skm1 !< Left cell average value
  real, intent(in) :: Sk   !< Center cell average value
  real, intent(in) :: Skp1 !< Right cell average value
  ! Local variables
  real :: h_sum, hp, hm

  h_sum = ( hkm1 + hkp1 ) + hk
  if (h_sum /= 0.) h_sum = 1./ h_sum
  hm =  hkm1 + hk
  if (hm /= 0.) hm = 1./ hm
  hp =  hkp1 + hk
  if (hp /= 0.) hp = 1./ hp
  fv_diff = ( hk * h_sum ) * &
            (   ( 2. * hkm1 + hk ) * hp * ( Skp1 - Sk ) &
              + ( 2. * hkp1 + hk ) * hm * ( Sk - Skm1 ) )
end function fv_diff

!> Returns the cell-centered second-order weigthed least squares slope
!! using three consecutive cell widths and average values. Slope is returned
!! as a gradient (i.e. units of scalar S over width units).
real function fvlsq_slope(hkm1, hk, hkp1, Skm1, Sk, Skp1)
  real, intent(in) :: hkm1 !< Left cell width
  real, intent(in) :: hk   !< Center cell width
  real, intent(in) :: hkp1 !< Right cell width
  real, intent(in) :: Skm1 !< Left cell average value
  real, intent(in) :: Sk   !< Center cell average value
  real, intent(in) :: Skp1 !< Right cell average value
  ! Local variables
  real :: xkm1, xkp1
  real :: h_sum, hx_sum, hxsq_sum, hxy_sum, hy_sum, det

  xkm1 = -0.5 * ( hk + hkm1 )
  xkp1 = 0.5 * ( hk + hkp1 )
  h_sum = ( hkm1 + hkp1 ) + hk
  hx_sum = hkm1*xkm1 + hkp1*xkp1
  hxsq_sum = hkm1*(xkm1**2) + hkp1*(xkp1**2)
  hxy_sum = hkm1*xkm1*Skm1 + hkp1*xkp1*Skp1
  hy_sum = ( hkm1*Skm1 + hkp1*Skp1 ) + hk*Sk
  det = h_sum * hxsq_sum - hx_sum**2
  if (det /= 0.) then
    !a = ( hxsq_sum * hy_sum - hx_sum*hxy_sum ) / det ! a would be mean of straight line fit
    fvlsq_slope = ( h_sum * hxy_sum - hx_sum*hy_sum ) / det ! Gradient of straight line fit
  else
    fvlsq_slope = 0. ! Adcroft's reciprocal rule
  endif
end function fvlsq_slope

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
  integer, dimension(2*nk+2), intent(inout) :: KoL   !< Index of first left interface below neutral surface
  integer, dimension(2*nk+2), intent(inout) :: KoR   !< Index of first right interface below neutral surface
  real, dimension(2*nk+1),    intent(inout) :: hEff  !< Effective thickness between two neutral surfaces (Pa)
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

!!! !> Returns the position between Pneg and Ppos where the interpolated density difference equals
!!! !! zero: Pint = ( dRhoPos * Ppos - dRhoNeg * Pneg ) / ( dRhoPos - dRhoneg )
!!! !! The result is always bounded to be between Pneg and Ppos.
!!! real function interpolate_for_position(dRhoNeg, Pneg, dRhoPos, Ppos)
!!!   real, intent(in) :: dRhoNeg !< Negative density difference
!!!   real, intent(in) :: Pneg    !< Position of negative density difference
!!!   real, intent(in) :: dRhoPos !< Positive density difference
!!!   real, intent(in) :: Ppos    !< Position of positive density difference
!!!   ! Local variables
!!!   real :: wghtU, wghtD, Pint
!!! 
!!!   if (Ppos<Pneg) stop 'interpolate_for_position: Houston, we have a problem! Ppos<Pneg'
!!!   if (dRhoNeg>dRhoPos) stop 'interpolate_for_position: Houston, we have a problem! dRhoNeg>dRhoPos'
!!!   if (Ppos<=Pneg) then ! Handle vanished or inverted layers
!!!     wghtU = 0.5
!!!     wghtD = 0.5
!!!     Pint = 0.5 * ( Pneg + Ppos )
!!!   elseif ( dRhoPos - dRhoNeg > 0. ) then
!!!     wghtU = -dRhoNeg / ( dRhoPos - dRhoNeg )
!!!     wghtD = dRhoPos / ( dRhoPos - dRhoNeg )
!!!     if ( wghtU < 0.5 ) then
!!!       Pint = Pneg + max( wghtU, 0. ) * ( Ppos - Pneg )
!!!     elseif ( wghtD < 0.5 ) then
!!!       Pint = Ppos + max( wghtD, 0. ) * ( Pneg - Ppos )
!!!     else
!!!       Pint = 0.5 * ( Pneg + Ppos )
!!!     endif
!!!   elseif ( dRhoPos - dRhoNeg == 0) then
!!!     if (dRhoNeg>0.) then
!!!       wghtU = 0.
!!!       wghtD = 1.
!!!       Pint = Pneg
!!!     elseif (dRhoNeg<0.) then
!!!       wghtU = 1.
!!!       wghtD = 0.
!!!       Pint = Ppos
!!!     else ! dRhoPos = dRhoNeg = 0
!!!       wghtU = 0.5
!!!       wghtD = 0.5
!!!       Pint = 0.5 * ( Pneg + Ppos )
!!!     endif
!!!   else ! dRho - dRhoNeg < 0
!!!     wghtU = 0.5
!!!     wghtD = 0.5
!!!     Pint = 0.5 * ( Pneg + Ppos )
!!!   endif
!!!   if ( Pint < Pneg ) stop 'interpolate_for_position: Houston, we have a problem! Pint < Pneg'
!!!   if ( Pint > Ppos ) stop 'interpolate_for_position: Houston, we have a problem! Pint > Ppos'
!!!   interpolate_for_position = Pint
!!! end function interpolate_for_position

!> Returns the non-dimensnional position between Pneg and Ppos where the interpolated density difference equals zero.
!! The result is always bounded to be between 0 and 1.
real function interpolate_for_nondim_position(dRhoNeg, Pneg, dRhoPos, Ppos)
  real, intent(in) :: dRhoNeg !< Negative density difference
  real, intent(in) :: Pneg    !< Position of negative density difference
  real, intent(in) :: dRhoPos !< Positive density difference
  real, intent(in) :: Ppos    !< Position of positive density difference

  if (Ppos<Pneg) stop 'interpolate_for_nondim_position: Houston, we have a problem! Ppos<Pneg'
  if (dRhoNeg>dRhoPos) stop 'interpolate_for_nondim_position: Houston, we have a problem! dRhoNeg>dRhoPos'
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
  if ( interpolate_for_nondim_position < 0. ) stop 'interpolate_for_nondim_position: Houston, we have a problem! Pint < Pneg'
  if ( interpolate_for_nondim_position > 1. ) stop 'interpolate_for_nondim_position: Houston, we have a problem! Pint > Ppos'
end function interpolate_for_nondim_position

!> Returns a single column of neutral diffusion fluxes of a tracer.
subroutine neutral_surface_flux(nk, Pl, Pr, Tl, Tr, PiL, PiR, KoL, KoR, hEff, Flx)
  integer,                    intent(in)    :: nk    !< Number of levels
  real, dimension(nk+1),      intent(in)    :: Pl    !< Left-column interface pressure (Pa)
  real, dimension(nk+1),      intent(in)    :: Pr    !< Right-column interface pressure (Pa)
  real, dimension(nk+1),      intent(in)    :: Tl    !< Left-column interface tracer (conc, e.g. degC)
  real, dimension(nk+1),      intent(in)    :: Tr    !< Right-column interface tracer (conc, e.g. degC)
  real, dimension(2*nk+2),    intent(in)    :: PiL   !< Position of neutral surface in left column (Pa)
  real, dimension(2*nk+2),    intent(in)    :: PiR   !< Position of neutral surface in right column (Pa)
  integer, dimension(2*nk+2), intent(in)    :: KoL   !< Index of first left interface below neutral surface
  integer, dimension(2*nk+2), intent(in)    :: KoR   !< Index of first right interface below neutral surface
  real, dimension(2*nk+1),    intent(in)    :: hEff  !< Effective thickness between two neutral surfaces (Pa)
  real, dimension(2*nk+1),    intent(inout) :: Flx   !< Flux of tracer between pairs of neutral layers (conc H)
  ! Local variables
  integer :: k_sublayer, kl, klm1, kr, krm1
  real :: T_right_top, T_right_bottom, T_left_top, T_left_bottom

  do k_sublayer = 1, 2*nk+1
 !  if (hEff(k_sublayer) == 0.) then
 !    Flx(k_sublayer) = 0.
 !  else

      kl = KoL(k_sublayer)
      klm1 = max(1, KoL(k_sublayer)-1)
      T_left_top = ( 1. - PiL(k_sublayer) ) * Tl(klm1) + PiL(k_sublayer) * Tl(kl)
!write(0,'(i3,2(i3,f8.2),2f8.2," left top")') k_sublayer,klm1,Tl(klm1),kl,Tl(kl),PiL(k_sublayer),T_left_top

      kl = KoL(k_sublayer+1)
      klm1 = max(1, KoL(k_sublayer+1)-1)
      T_left_bottom = ( 1. - PiL(k_sublayer+1) ) * Tl(klm1) + PiL(k_sublayer+1) * Tl(kl)
!write(0,'(i3,2(i3,f8.2),2f8.2," left bottom")') k_sublayer+1,klm1,Tl(klm1),kl,Tl(kl),PiL(k_sublayer+1),T_left_bottom

      kr = KoR(k_sublayer)
      krm1 = max(1, KoR(k_sublayer)-1)
      T_right_top = ( 1. - PiR(k_sublayer) ) * Tr(krm1) + PiR(k_sublayer) * Tr(kr)
!write(0,'(i3,2(i3,f8.2),2f8.2," right top")') k_sublayer,krm1,Tr(krm1),kr,Tr(kr),PiR(k_sublayer),T_right_top

      kr = KoR(k_sublayer+1)
      krm1 = max(1, KoR(k_sublayer+1)-1)
      T_right_bottom = ( 1. - PiR(k_sublayer+1) ) * Tr(krm1) + PiR(k_sublayer+1) * Tr(kr)
!write(0,'(i3,2(i3,f8.2),2f8.2," right bottom")') k_sublayer+1,krm1,Tr(krm1),kr,Tr(kr),PiR(k_sublayer+1),T_right_bottom

      Flx(k_sublayer) = 0.5 * ( ( T_right_top - T_left_top ) + ( T_right_bottom - T_left_bottom ) )
!write(0,'(i3,f8.3)') k_sublayer, Flx(k_sublayer) * hEff(k_sublayer)
 !  endif
  enddo

end subroutine neutral_surface_flux

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
  real, dimension(2*nk+1) :: Flx ! Test flux
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
  data kL0 / 1, 2, 2, 3, 3, 4, 4, 5, 5, 5 /
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

  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fv_diff(1.,1.,1., 0.,1.,2., 1., 'FV: Straight line on uniform grid')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fv_diff(1.,1.,0., 0.,4.,8., 7., 'FV: Vanished right cell')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fv_diff(0.,1.,1., 0.,4.,8., 7., 'FV: Vanished left cell')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fv_diff(1.,2.,4., 0.,3.,9., 4., 'FV: Stretched grid')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fv_diff(2.,0.,2., 0.,1.,2., 0., 'FV: Vanished middle cell')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fv_diff(0.,1.,0., 0.,1.,2., 2., 'FV: Vanished on both sides')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fv_diff(1.,0.,0., 0.,1.,2., 0., 'FV: Two vanished cell sides')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fv_diff(0.,0.,0., 0.,1.,2., 0., 'FV: All vanished cells')

  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fvlsq_slope(1.,1.,1., 0.,1.,2., 1., 'LSQ: Straight line on uniform grid')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fvlsq_slope(1.,1.,0., 0.,1.,2., 1., 'LSQ: Vanished right cell')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fvlsq_slope(0.,1.,1., 0.,1.,2., 1., 'LSQ: Vanished left cell')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fvlsq_slope(1.,2.,4., 0.,3.,9., 2., 'LSQ: Stretched grid')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fvlsq_slope(1.,0.,1., 0.,1.,2., 2., 'LSQ: Vanished middle cell')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fvlsq_slope(0.,1.,0., 0.,1.,2., 0., 'LSQ: Vanished on both sides')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fvlsq_slope(1.,0.,0., 0.,1.,2., 0., 'LSQ: Two vanished cell sides')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fvlsq_slope(0.,0.,0., 0.,1.,2., 0., 'LSQ: All vanished cells')

  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifndp(-1.0, 0.,  1.0, 1.0, 0.5, 'Check mid-point')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifndp( 0.0, 0.,  1.0, 1.0, 0.0, 'Check bottom')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifndp( 0.1, 0.,  1.1, 1.0, 0.0, 'Check below')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifndp(-1.0, 0.,  0.0, 1.0, 1.0, 'Check top')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifndp(-1.0, 0., -0.1, 1.0, 1.0, 'Check above')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifndp(-1.0, 0.,  3.0, 1.0, 0.25, 'Check 1/4')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifndp(-3.0, 0.,  1.0, 1.0, 0.75, 'Check 3/4')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifndp( 1.0, 0.,  1.0, 1.0, 0.0, 'Check dRho=0 below')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifndp(-1.0, 0., -1.0, 1.0, 1.0, 'Check dRho=0 above')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifndp( 0.0, 0.,  0.0, 1.0, 0.5, 'Check dRho=0 mid')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_ifndp(-2.0, .5,  5.0, 0.5, 0.5, 'Check dP=0')

  call find_neutral_surface_positions(nk, PiL, TiL, SiL, dRdt, dRdS, PiL, TiL, SiL, dRdT, dRdS, PiLRo, PiRLo, KoL, KoR, hEff)
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL, KoL, PiLRo), pL0, 'Identical columns, left positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+2, absolute_positions(nk, PiL, KoR, PiRLo), pR0, 'Identical columns, right positions')
  neutralDiffusionUnitTests = neutralDiffusionUnitTests .or. test_fnsp(2*nk+1, hEff, hE0, 'Identical columns, thicknesses')

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
  call neutral_surface_flux(nk, PiL, PiL, TiL, TiL+2., PiLRo, PiRLo, KoL, KoR, hEff, Flx)

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

  !> Returns true if a test of fv_diff() fails, and conditionally writes results to stream
  logical function test_fv_diff(hkm1, hk, hkp1, Skm1, Sk, Skp1, Ptrue, title)
    real,             intent(in) :: hkm1 !< Left cell width
    real,             intent(in) :: hk   !< Center cell width
    real,             intent(in) :: hkp1 !< Right cell width
    real,             intent(in) :: Skm1 !< Left cell average value
    real,             intent(in) :: Sk   !< Center cell average value
    real,             intent(in) :: Skp1 !< Right cell average value
    real,             intent(in) :: Ptrue  !< True answer (Pa)
    character(len=*), intent(in) :: title !< Title for messages
    ! Local variables
    integer :: stdunit
    real :: Pret

    Pret = fv_diff(hkm1, hk, hkp1, Skm1, Sk, Skp1)
    test_fv_diff = (Pret /= Ptrue) 

    if (test_fv_diff .or. verbosity>5) then
      stdunit = 6
      if (test_fv_diff.or.debug_this_module) stdunit = 0 ! In case of wrong results, write to error stream
      write(stdunit,'(a)') title
      if (test_fv_diff) then
        write(stdunit,'(2(x,a,f20.16),x,a)') 'pRet=',Pret,'pTrue=',Ptrue,'WRONG!'
      else
        write(stdunit,'(2(x,a,f20.16))') 'pRet=',Pret,'pTrue=',Ptrue
      endif
    endif

  end function test_fv_diff

  !> Returns true if a test of fvlsq_slope() fails, and conditionally writes results to stream
  logical function test_fvlsq_slope(hkm1, hk, hkp1, Skm1, Sk, Skp1, Ptrue, title)
    real,             intent(in) :: hkm1 !< Left cell width
    real,             intent(in) :: hk   !< Center cell width
    real,             intent(in) :: hkp1 !< Right cell width
    real,             intent(in) :: Skm1 !< Left cell average value
    real,             intent(in) :: Sk   !< Center cell average value
    real,             intent(in) :: Skp1 !< Right cell average value
    real,             intent(in) :: Ptrue  !< True answer (Pa)
    character(len=*), intent(in) :: title !< Title for messages
    ! Local variables
    integer :: stdunit
    real :: Pret

    Pret = fvlsq_slope(hkm1, hk, hkp1, Skm1, Sk, Skp1)
    test_fvlsq_slope = (Pret /= Ptrue) 

    if (test_fvlsq_slope .or. verbosity>5) then
      stdunit = 6
      if (test_fvlsq_slope.or.debug_this_module) stdunit = 0 ! In case of wrong results, write to error stream
      write(stdunit,'(a)') title
      if (test_fvlsq_slope) then
        write(stdunit,'(2(x,a,f20.16),x,a)') 'pRet=',Pret,'pTrue=',Ptrue,'WRONG!'
      else
        write(stdunit,'(2(x,a,f20.16))') 'pRet=',Pret,'pTrue=',Ptrue
      endif
    endif

  end function test_fvlsq_slope

  !> Returns true if a test of interpolate_for_nondim_position() fails, and conditionally writes results to stream
  logical function test_ifndp(rhoNeg, Pneg, rhoPos, Ppos, Ptrue, title)
    real,             intent(in) :: rhoNeg !< Lighter density (kg/m3)
    real,             intent(in) :: Pneg   !< Interface position of lighter density (pa)
    real,             intent(in) :: rhoPos !< Heavier density (kg/m3)
    real,             intent(in) :: Ppos   !< Interface position of heavier density (pa)
    real,             intent(in) :: Ptrue  !< True answer (Pa)
    character(len=*), intent(in) :: title !< Title for messages
    ! Local variables
    integer :: stdunit
    real :: Pret

    Pret = interpolate_for_nondim_position(rhoNeg, Pneg, rhoPos, Ppos)
    test_ifndp = (Pret /= Ptrue) 

    if (test_ifndp .or. verbosity>5) then
      stdunit = 6
      if (test_ifndp.or.debug_this_module) stdunit = 0 ! In case of wrong results, write to error stream
      write(stdunit,'(a)') title
      if (test_ifndp) then
        write(stdunit,'(4(x,a,f20.16),2(x,a,1pe22.15),x,a)') 'r1=',rhoNeg,'p1=',Pneg,'r2=',rhoPos,'p2=',Ppos,'pRet=',Pret,'pTrue=',Ptrue,'WRONG!'
      else
        write(stdunit,'(4(x,a,f20.16),2(x,a,1pe22.15))') 'r1=',rhoNeg,'p1=',Pneg,'r2=',rhoPos,'p2=',Ppos,'pRet=',Pret,'pTrue=',Ptrue
      endif
    endif

  end function test_ifndp

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
