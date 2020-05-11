!> Calculates and applies diffusive fluxes as a parameterization of lateral mixing (non-neutral) by
!! mesoscale eddies near the top and bottom (to be implemented) boundary layers of the ocean.

module MOM_lateral_boundary_diffusion

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,             only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,             only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_domains,               only : pass_var
use MOM_diag_mediator,         only : diag_ctrl, time_type
use MOM_diag_mediator,         only : post_data, register_diag_field
use MOM_error_handler,         only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_file_parser,           only : openParameterBlock, closeParameterBlock
use MOM_grid,                  only : ocean_grid_type
use MOM_remapping,             only : remapping_CS, initialize_remapping
use MOM_remapping,             only : extract_member_remapping_CS, build_reconstructions_1d
use MOM_remapping,             only : average_value_ppoly, remappingSchemesDoc, remappingDefaultScheme
use MOM_tracer_registry,       only : tracer_registry_type, tracer_type
use MOM_unit_scaling,          only : unit_scale_type
use MOM_verticalGrid,          only : verticalGrid_type
use MOM_CVMix_KPP,             only : KPP_get_BLD, KPP_CS
use MOM_energetic_PBL,         only : energetic_PBL_get_MLD, energetic_PBL_CS
use MOM_diabatic_driver,       only : diabatic_CS, extract_diabatic_member

use iso_fortran_env, only : stdout=>output_unit, stderr=>error_unit

implicit none ; private

public near_boundary_unit_tests, lateral_boundary_diffusion, lateral_boundary_diffusion_init
public boundary_k_range

! Private parameters to avoid doing string comparisons for bottom or top boundary layer
integer, public, parameter :: SURFACE = -1 !< Set a value that corresponds to the surface bopundary
integer, public, parameter :: BOTTOM  = 1  !< Set a value that corresponds to the bottom boundary
#include <MOM_memory.h>

!> Sets parameters for lateral boundary mixing module.
type, public :: lateral_boundary_diffusion_CS ; private
  integer :: method                                               !< Determine which of the three methods calculate
                                                                  !! and apply near boundary layer fluxes
                                                                  !! 1. Bulk-layer approach
                                                                  !! 2. Along layer
  integer :: deg                                                  !< Degree of polynomial reconstruction
  integer :: surface_boundary_scheme                              !< Which boundary layer scheme to use
                                                                  !! 1. ePBL; 2. KPP
  logical :: limiter                                              !< Controls wether a flux limiter is applied.
                                                                  !! Only valid when method = 1.
  type(remapping_CS)              :: remap_CS                     !< Control structure to hold remapping configuration
  type(KPP_CS),           pointer :: KPP_CSp => NULL()            !< KPP control structure needed to get BLD
  type(energetic_PBL_CS), pointer :: energetic_PBL_CSp => NULL()  !< ePBL control structure needed to get BLD
  type(diag_ctrl), pointer :: diag => NULL()                      !< A structure that is used to
                                                                  !! regulate the timing of diagnostic output.
end type lateral_boundary_diffusion_CS

! This include declares and sets the variable "version".
#include "version_variable.h"
character(len=40) :: mdl = "MOM_lateral_boundary_diffusion"       !< Name of this module

contains

!> Initialization routine that reads runtime parameters and sets up pointers to other control structures that might be
!! needed for lateral boundary diffusion.
logical function lateral_boundary_diffusion_init(Time, G, param_file, diag, diabatic_CSp, CS)
  type(time_type), target,          intent(in)    :: Time          !< Time structure
  type(ocean_grid_type),            intent(in)    :: G             !< Grid structure
  type(param_file_type),            intent(in)    :: param_file    !< Parameter file structure
  type(diag_ctrl), target,          intent(inout) :: diag          !< Diagnostics control structure
  type(diabatic_CS),                pointer       :: diabatic_CSp  !< KPP control structure needed to get BLD
  type(lateral_boundary_diffusion_CS), pointer    :: CS            !< Lateral boundary mixing control structure

  ! local variables
  character(len=80)  :: string  ! Temporary strings
  logical :: boundary_extrap

  if (ASSOCIATED(CS)) then
    call MOM_error(FATAL, "lateral_boundary_diffusion_init called with associated control structure.")
    return
  endif

  ! Log this module and master switch for turning it on/off
  call log_version(param_file, mdl, version, &
       "This module implements lateral diffusion of tracers near boundaries")
  call get_param(param_file, mdl, "USE_LATERAL_BOUNDARY_DIFFUSION", lateral_boundary_diffusion_init, &
                 "If true, enables the lateral boundary tracer's diffusion module.", &
                 default=.false.)

  if (.not. lateral_boundary_diffusion_init) then
    return
  endif

  allocate(CS)
  CS%diag => diag
  call extract_diabatic_member(diabatic_CSp, KPP_CSp=CS%KPP_CSp)
  call extract_diabatic_member(diabatic_CSp, energetic_PBL_CSp=CS%energetic_PBL_CSp)

  CS%surface_boundary_scheme = -1
  if ( .not. ASSOCIATED(CS%energetic_PBL_CSp) .and. .not. ASSOCIATED(CS%KPP_CSp) ) then
    call MOM_error(FATAL,"Lateral boundary diffusion is true, but no valid boundary layer scheme was found")
  endif

  ! Read all relevant parameters and write them to the model log.
  call get_param(param_file, mdl, "LATERAL_BOUNDARY_METHOD", CS%method, &
                 "Determine how to apply boundary lateral diffusion of tracers: \n"//&
                 "1. Bulk layer approach  \n"//&
                 "2. Along layer approach", default=1)
  if (CS%method == 1) then
    call get_param(param_file, mdl, "APPLY_LIMITER", CS%limiter, &
                   "If True, apply a flux limiter in the LBD. This is only available \n"//&
                   "when LATERAL_BOUNDARY_METHOD=1.", default=.false.)
  endif
  call get_param(param_file, mdl, "LBD_BOUNDARY_EXTRAP", boundary_extrap, &
                 "Use boundary extrapolation in LBD code", &
                 default=.false.)
  call get_param(param_file, mdl, "LBD_REMAPPING_SCHEME", string, &
                 "This sets the reconstruction scheme used "//&
                 "for vertical remapping for all variables. "//&
                 "It can be one of the following schemes: "//&
                 trim(remappingSchemesDoc), default=remappingDefaultScheme)
  call initialize_remapping( CS%remap_CS, string, boundary_extrapolation = boundary_extrap )
  call extract_member_remapping_CS(CS%remap_CS, degree=CS%deg)

end function lateral_boundary_diffusion_init

!> Driver routine for calculating lateral diffusive fluxes near the top and bottom boundaries.
!! Two different methods are available:
!! Method 1: lower order representation, calculate fluxes from bulk layer integrated quantities.
!! Method 2: more straight forward, diffusion is applied layer by layer using only information
!! from neighboring cells.
subroutine lateral_boundary_diffusion(G, GV, US, h, Coef_x, Coef_y, dt, Reg, CS)
  type(ocean_grid_type),                intent(inout) :: G   !< Grid type
  type(verticalGrid_type),              intent(in)    :: GV  !< ocean vertical grid structure
  type(unit_scale_type),                intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                                        intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G)),    intent(in)    :: Coef_x !< dt * Kh * dy / dx at u-points [L2 ~> m2]
  real, dimension(SZI_(G),SZJB_(G)),    intent(in)    :: Coef_y !< dt * Kh * dx / dy at v-points [L2 ~> m2]
  real,                                 intent(in)    :: dt     !< Tracer time step * I_numitts
                                                                !! (I_numitts in tracer_hordiff)
  type(tracer_registry_type),           pointer       :: Reg    !< Tracer registry
  type(lateral_boundary_diffusion_CS),  intent(in)    :: CS     !< Control structure for this module

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: hbl                           !< bnd. layer depth [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G),CS%deg+1) :: ppoly0_coefs !< Coefficients of polynomial
  real, dimension(SZI_(G),SZJ_(G),SZK_(G),2)        :: ppoly0_E     !< Edge values from reconstructions
  real, dimension(SZK_(G),CS%deg+1)                 :: ppoly_S      !< Slopes from reconstruction (placeholder)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: uFlx        !< Zonal flux of tracer [conc m^3]
  real, dimension(SZIB_(G),SZJ_(G))         :: uFLx_bulk   !< Total calculated bulk-layer u-flux for the tracer
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: vFlx        !< Meridional flux of tracer [conc m^3]
  real, dimension(SZI_(G),SZJB_(G))         :: vFlx_bulk   !< Total calculated bulk-layer v-flux for the tracer
  real, dimension(SZIB_(G),SZJ_(G))         :: uwork_2d    !< Layer summed u-flux transport
  real, dimension(SZI_(G),SZJB_(G))         :: vwork_2d    !< Layer summed v-flux transport
  real, dimension(SZI_(G),SZJ_(G),G%ke)     :: tendency    !< tendency array for diagn
  real, dimension(SZI_(G),SZJ_(G))          :: tendency_2d !< depth integrated content tendency for diagn
  type(tracer_type), pointer                :: Tracer => NULL() !< Pointer to the current tracer
  integer :: remap_method !< Reconstruction method
  integer :: i,j,k,m      !< indices to loop over
  real    :: Idt          !< inverse of the time step [s-1]

  Idt = 1./dt
  hbl(:,:) = 0.
  if (ASSOCIATED(CS%KPP_CSp)) call KPP_get_BLD(CS%KPP_CSp, hbl, G, US, m_to_BLD_units=GV%m_to_H)
  if (ASSOCIATED(CS%energetic_PBL_CSp)) &
    call energetic_PBL_get_MLD(CS%energetic_PBL_CSp, hbl, G, US, m_to_MLD_units=GV%m_to_H)

  call pass_var(hbl,G%Domain)
  do m = 1,Reg%ntr
    tracer => Reg%tr(m)

    ! for diagnostics
    if (tracer%id_lbdxy_conc > 0 .or. tracer%id_lbdxy_cont > 0 .or. tracer%id_lbdxy_cont_2d > 0) then
      tendency(:,:,:) = 0.0
    endif

    ! Interpolate state to interface
    do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
      call build_reconstructions_1d( CS%remap_CS, G%ke, h(i,j,:), tracer%t(i,j,:), ppoly0_coefs(i,j,:,:), &
                                     ppoly0_E(i,j,:,:), ppoly_S, remap_method, GV%H_subroundoff, GV%H_subroundoff)
    enddo ; enddo
    ! Diffusive fluxes in the i-direction
    uFlx(:,:,:) = 0.
    vFlx(:,:,:) = 0.
    uFlx_bulk(:,:) = 0.
    vFlx_bulk(:,:) = 0.

    ! Method #1
    if ( CS%method == 1 ) then
      do j=G%jsc,G%jec
        do i=G%isc-1,G%iec
          if (G%mask2dCu(I,j)>0.) then
            call fluxes_bulk_method(SURFACE, GV%ke, CS%deg, h(I,j,:), h(I+1,j,:), hbl(I,j), hbl(I+1,j), &
              G%areaT(I,j), G%areaT(I+1,j), tracer%t(I,j,:), tracer%t(I+1,j,:),                         &
              ppoly0_coefs(I,j,:,:), ppoly0_coefs(I+1,j,:,:), ppoly0_E(I,j,:,:),                        &
              ppoly0_E(I+1,j,:,:), remap_method, Coef_x(I,j), uFlx_bulk(I,j), uFlx(I,j,:), CS%limiter)
          endif
        enddo
      enddo
      do J=G%jsc-1,G%jec
        do i=G%isc,G%iec
          if (G%mask2dCv(i,J)>0.) then
            call fluxes_bulk_method(SURFACE, GV%ke, CS%deg, h(i,J,:), h(i,J+1,:), hbl(i,J), hbl(i,J+1), &
              G%areaT(i,J), G%areaT(i,J+1), tracer%t(i,J,:), tracer%t(i,J+1,:),                         &
              ppoly0_coefs(i,J,:,:), ppoly0_coefs(i,J+1,:,:), ppoly0_E(i,J,:,:),                        &
              ppoly0_E(i,J+1,:,:), remap_method, Coef_y(i,J), vFlx_bulk(i,J), vFlx(i,J,:), CS%limiter)
          endif
        enddo
      enddo
      ! Post tracer bulk diags
      if (tracer%id_lbd_bulk_dfx>0) call post_data(tracer%id_lbd_bulk_dfx, uFlx_bulk*Idt, CS%diag)
      if (tracer%id_lbd_bulk_dfy>0) call post_data(tracer%id_lbd_bulk_dfy, vFlx_bulk*Idt, CS%diag)

    ! Method #2
    elseif (CS%method == 2) then
      do j=G%jsc,G%jec
        do i=G%isc-1,G%iec
          if (G%mask2dCu(I,j)>0.) then
            call fluxes_layer_method(SURFACE, GV%ke, CS%deg, h(I,j,:), h(I+1,j,:), hbl(I,j), hbl(I+1,j), &
              G%areaT(I,j), G%areaT(I+1,j), tracer%t(I,j,:), tracer%t(I+1,j,:), ppoly0_coefs(I,j,:,:), &
              ppoly0_coefs(I+1,j,:,:), ppoly0_E(I,j,:,:), ppoly0_E(I+1,j,:,:), remap_method, Coef_x(I,j), uFlx(I,j,:))
          endif
        enddo
      enddo
      do J=G%jsc-1,G%jec
        do i=G%isc,G%iec
          if (G%mask2dCv(i,J)>0.) then
            call fluxes_layer_method(SURFACE, GV%ke, CS%deg, h(i,J,:), h(i,J+1,:), hbl(i,J), hbl(i,J+1), &
              G%areaT(i,J), G%areaT(i,J+1), tracer%t(i,J,:), tracer%t(i,J+1,:), ppoly0_coefs(i,J,:,:), &
              ppoly0_coefs(i,J+1,:,:), ppoly0_E(i,J,:,:), ppoly0_E(i,J+1,:,:), remap_method, Coef_y(i,J), vFlx(i,J,:))
          endif
        enddo
      enddo
    endif

    ! Update the tracer fluxes
    do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if (G%mask2dT(i,j)>0.) then
        tracer%t(i,j,k) = tracer%t(i,j,k) + (( (uFlx(I-1,j,k)-uFlx(I,j,k)) ) + ( (vFlx(i,J-1,k)-vFlx(i,J,k) ) ))* &
                          (G%IareaT(i,j)/( h(i,j,k) + GV%H_subroundoff))

        if (tracer%id_lbdxy_conc > 0  .or. tracer%id_lbdxy_cont > 0 .or. tracer%id_lbdxy_cont_2d > 0 ) then
          tendency(i,j,k) = (( (uFlx(I-1,j,k)-uFlx(I,j,k)) ) + ( (vFlx(i,J-1,k)-vFlx(i,J,k) ) )) * G%IareaT(i,j) * Idt
        endif

      endif
    enddo ; enddo ; enddo

    ! Post the tracer diagnostics
    if (tracer%id_lbd_dfx>0)      call post_data(tracer%id_lbd_dfx, uFlx*Idt, CS%diag)
    if (tracer%id_lbd_dfy>0)      call post_data(tracer%id_lbd_dfy, vFlx*Idt, CS%diag)
    if (tracer%id_lbd_dfx_2d>0) then
      uwork_2d(:,:) = 0.
      do k=1,GV%ke ; do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
        uwork_2d(I,j) = uwork_2d(I,j) + (uFlx(I,j,k) * Idt)
      enddo ; enddo ; enddo
      call post_data(tracer%id_lbd_dfx_2d, uwork_2d, CS%diag)
    endif

    if (tracer%id_lbd_dfy_2d>0) then
      vwork_2d(:,:) = 0.
      do k=1,GV%ke ; do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
        vwork_2d(i,J) = vwork_2d(i,J) + (vFlx(i,J,k) * Idt)
      enddo ; enddo ; enddo
      call post_data(tracer%id_lbd_dfy_2d, vwork_2d, CS%diag)
    endif

    ! post tendency of tracer content
    if (tracer%id_lbdxy_cont > 0) then
      call post_data(tracer%id_lbdxy_cont, tendency, CS%diag)
    endif

    ! post depth summed tendency for tracer content
    if (tracer%id_lbdxy_cont_2d > 0) then
      tendency_2d(:,:) = 0.
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        do k=1,GV%ke
          tendency_2d(i,j) = tendency_2d(i,j) + tendency(i,j,k)
        enddo
      enddo ; enddo
      call post_data(tracer%id_lbdxy_cont_2d, tendency_2d, CS%diag)
    endif

    ! post tendency of tracer concentration; this step must be
    ! done after posting tracer content tendency, since we alter
    ! the tendency array and its units.
    if (tracer%id_lbdxy_conc > 0) then
      do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
        tendency(i,j,k) =  tendency(i,j,k) / ( h(i,j,k) + GV%H_subroundoff )
      enddo ; enddo ; enddo
      call post_data(tracer%id_lbdxy_conc, tendency, CS%diag)
    endif

  enddo

end subroutine lateral_boundary_diffusion

!< Calculate bulk layer value of a scalar quantity as the thickness weighted average
real function bulk_average(boundary, nk, deg, h, hBLT, phi, ppoly0_E, ppoly0_coefs, method, k_top, zeta_top, k_bot, &
                           zeta_bot)
  integer             :: boundary          !< SURFACE or BOTTOM                                         [nondim]
  integer             :: nk                !< Number of layers                                          [nondim]
  integer             :: deg               !< Degree of polynomial                                      [nondim]
  real, dimension(nk) :: h                 !< Layer thicknesses                               [H ~> m or kg m-2]
  real                :: hBLT              !< Depth of the boundary layer                     [H ~> m or kg m-2]
  real, dimension(nk) :: phi               !< Scalar quantity
  real, dimension(nk,2)     :: ppoly0_E     !< Edge value of polynomial
  real, dimension(nk,deg+1) :: ppoly0_coefs !< Coefficients of polynomial
  integer                   :: method       !< Remapping scheme to use

  integer             :: k_top             !< Index of the first layer within the boundary
  real                :: zeta_top          !< Fraction of the layer encompassed by the bottom boundary layer
                                           !! (0 if none, 1. if all). For the surface, this is always 0. because
                                           !! integration starts at the surface                         [nondim]
  integer             :: k_bot             !< Index of the last layer within the boundary
  real                :: zeta_bot          !< Fraction of the layer encompassed by the surface boundary layer
                                           !! (0 if none, 1. if all). For the bottom boundary layer, this is always 1.
                                           !! because integration starts at the bottom                  [nondim]
  ! Local variables
  real    :: htot !< Running sum of the thicknesses (top to bottom) [H ~> m or kg m-2]
  integer :: k    !< k indice


  htot = 0.
  bulk_average = 0.
  if (hblt == 0.) return
  if (boundary == SURFACE) then
    htot = (h(k_bot) * zeta_bot)
    bulk_average = average_value_ppoly( nk, phi, ppoly0_E, ppoly0_coefs, method, k_bot, 0., zeta_bot) * htot
    do k = k_bot-1,1,-1
      bulk_average = bulk_average + phi(k)*h(k)
      htot = htot + h(k)
    enddo
  elseif (boundary == BOTTOM) then
    htot = (h(k_top) * zeta_top)
    ! (note 1-zeta_top because zeta_top is the fraction of the layer)
    bulk_average = average_value_ppoly( nk, phi, ppoly0_E, ppoly0_coefs, method, k_top, (1.-zeta_top), 1.) * htot
    do k = k_top+1,nk
      bulk_average = bulk_average + phi(k)*h(k)
      htot = htot + h(k)
    enddo
  else
    call MOM_error(FATAL, "bulk_average: a valid boundary type must be provided.")
  endif

  bulk_average = bulk_average / hBLT

end function bulk_average

!> Calculate the harmonic mean of two quantities
!! See \ref section_harmonic_mean.
real function harmonic_mean(h1,h2)
  real :: h1 !< Scalar quantity
  real :: h2 !< Scalar quantity
  if (h1 + h2 == 0.) then
    harmonic_mean = 0.
  else
    harmonic_mean = 2.*(h1*h2)/(h1+h2)
  endif
end function harmonic_mean

!> Find the k-index range corresponding to the layers that are within the boundary-layer region
subroutine boundary_k_range(boundary, nk, h, hbl, k_top, zeta_top, k_bot, zeta_bot)
  integer,             intent(in   ) :: boundary !< SURFACE or BOTTOM                       [nondim]
  integer,             intent(in   ) :: nk       !< Number of layers                        [nondim]
  real, dimension(nk), intent(in   ) :: h        !< Layer thicknesses of the column         [H ~> m or kg m-2]
  real,                intent(in   ) :: hbl      !< Thickness of the boundary layer         [H ~> m or kg m-2]
                                                 !! If surface, with respect to zbl_ref = 0.
                                                 !! If bottom, with respect to zbl_ref = SUM(h)
  integer,             intent(  out) :: k_top    !< Index of the first layer within the boundary
  real,                intent(  out) :: zeta_top !< Distance from the top of a layer to the intersection of the
                                                 !! top extent of the boundary layer (0 at top, 1 at bottom)  [nondim]
  integer,             intent(  out) :: k_bot    !< Index of the last layer within the boundary
  real,                intent(  out) :: zeta_bot !< Distance of the lower layer to the boundary layer depth
                                                 !! (0 at top, 1 at bottom)  [nondim]
  ! Local variables
  real :: htot ! Summed thickness [H ~> m or kg m-2]
  integer :: k
  ! Surface boundary layer
  if ( boundary == SURFACE ) then
    k_top = 1
    zeta_top = 0.
    htot = 0.
    k_bot = 1
    zeta_bot = 0.
    if (hbl == 0.) return
    if (hbl >= SUM(h(:))) then
      k_bot = nk
      zeta_bot = 1.
      return
    endif
    do k=1,nk
      htot = htot + h(k)
      if ( htot >= hbl) then
        k_bot = k
        zeta_bot = 1 - (htot - hbl)/h(k)
        return
      endif
    enddo
  ! Bottom boundary layer
  elseif ( boundary == BOTTOM ) then
    k_top = nk
    zeta_top = 1.
    k_bot = nk
    zeta_bot = 0.
    htot = 0.
    if (hbl == 0.) return
    if (hbl >= SUM(h(:))) then
      k_top = 1
      zeta_top = 1.
      return
    endif
    do k=nk,1,-1
      htot = htot + h(k)
      if (htot >= hbl) then
        k_top = k
        zeta_top = 1 - (htot - hbl)/h(k)
        return
      endif
    enddo
  else
    call MOM_error(FATAL,"Houston, we've had a problem in boundary_k_range")
  endif

end subroutine boundary_k_range


!> Calculate the lateral boundary diffusive fluxes using the layer by layer method.
!! See \ref section_method2
subroutine fluxes_layer_method(boundary, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, &
                              ppoly0_coefs_L, ppoly0_coefs_R, ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_layer)

  integer,                   intent(in   )       :: boundary !< Which boundary layer SURFACE or BOTTOM  [nondim]
  integer,                   intent(in   )       :: nk       !< Number of layers                        [nondim]
  integer,                   intent(in   )       :: deg      !< order of the polynomial reconstruction  [nondim]
  real, dimension(nk),       intent(in   )       :: h_L      !< Layer thickness (left)        [H ~> m or kg m-2]
  real, dimension(nk),       intent(in   )       :: h_R      !< Layer thickness (right)       [H ~> m or kg m-2]
  real,                      intent(in   )       :: hbl_L    !< Thickness of the boundary boundary
                                                                       !! layer (left)        [H ~> m or kg m-2]
  real,                      intent(in   )       :: hbl_R    !< Thickness of the boundary boundary
                                                             !! layer (right)                 [H ~> m or kg m-2]
  real,                      intent(in   )       :: area_L   !< Area of the horizontal grid (left)    [L2 ~> m2]
  real,                      intent(in   )       :: area_R   !< Area of the horizontal grid (right)   [L2 ~> m2]
  real, dimension(nk),       intent(in   )       :: phi_L    !< Tracer values (left)                    [conc]
  real, dimension(nk),       intent(in   )       :: phi_R    !< Tracer values (right)                   [conc]
  real, dimension(nk,deg+1), intent(in   )       :: ppoly0_coefs_L !< Tracer reconstruction (left)      [conc]
  real, dimension(nk,deg+1), intent(in   )       :: ppoly0_coefs_R !< Tracer reconstruction (right)     [conc]
  real, dimension(nk,2),     intent(in   )       :: ppoly0_E_L !< Polynomial edge values (left)         [ nondim ]
  real, dimension(nk,2),     intent(in   )       :: ppoly0_E_R !< Polynomial edge values (right)        [ nondim ]
  integer,                   intent(in   )       :: method   !< Method of polynomial integration        [ nondim ]
  real,                      intent(in   )       :: khtr_u   !< Horizontal diffusivities times delta t
                                                             !! at a velocity point [L2 ~> m2]
  real, dimension(nk),       intent(  out)       :: F_layer  !< Layerwise diffusive flux at U- or V-point
                                                             !! [H L2 conc ~> m3 conc]

  ! Local variables
  real, dimension(nk) :: h_means              !< Calculate the layer-wise harmonic means      [H ~> m or kg m-2]
  real                :: khtr_avg             !< Thickness-weighted diffusivity at the u-point        [m^2 s^-1]
                                              !! This is just to remind developers that khtr_avg should be
                                              !! computed once khtr is 3D.
  real                :: heff                 !< Harmonic mean of layer thicknesses           [H ~> m or kg m-2]
  real                :: inv_heff             !< Inverse of the harmonic mean of layer thicknesses
                                              !!  [H-1 ~> m-1 or m2 kg-1]
  real                :: phi_L_avg, phi_R_avg !< Bulk, thickness-weighted tracer averages (left and right column)
                                              !!                                                    [conc m^-3 ]
  real    :: htot                      !< Total column thickness [H ~> m or kg m-2]
  integer :: k, k_bot_min, k_top_max   !< k-indices, min and max for top and bottom, respectively
  integer :: k_top_L, k_bot_L          !< k-indices left
  integer :: k_top_R, k_bot_R          !< k-indices right
  real    :: zeta_top_L, zeta_top_R    !< distance from the top of a layer to the boundary
                                       !! layer depth                                     [nondim]
  real    :: zeta_bot_L, zeta_bot_R    !< distance from the bottom of a layer to the boundary
                                       !!layer depth                                      [nondim]
  real    :: h_work_L, h_work_R  !< dummy variables
  real    :: hbl_min             !< minimum BLD (left and right)                                   [m]

  F_layer(:) = 0.0
  if (hbl_L == 0. .or. hbl_R == 0.) then
    return
  endif

  ! Calculate vertical indices containing the boundary layer
  call boundary_k_range(boundary, nk, h_L, hbl_L, k_top_L, zeta_top_L, k_bot_L, zeta_bot_L)
  call boundary_k_range(boundary, nk, h_R, hbl_R, k_top_R, zeta_top_R, k_bot_R, zeta_bot_R)

  if (boundary == SURFACE) then
    k_bot_min = MIN(k_bot_L, k_bot_R)
    ! make sure left and right k indices span same range
    if (k_bot_min .ne. k_bot_L) then
      k_bot_L = k_bot_min
      zeta_bot_L = 1.0
    endif
    if (k_bot_min .ne. k_bot_R) then
      k_bot_R= k_bot_min
      zeta_bot_R = 1.0
    endif

    h_work_L = (h_L(k_bot_L) * zeta_bot_L)
    h_work_R = (h_R(k_bot_R) * zeta_bot_R)

    phi_L_avg = average_value_ppoly( nk, phi_L, ppoly0_E_L, ppoly0_coefs_L, method, k_bot_L, 0., zeta_bot_L)
    phi_R_avg = average_value_ppoly( nk, phi_R, ppoly0_E_R, ppoly0_coefs_R, method, k_bot_R, 0., zeta_bot_R)
    heff = harmonic_mean(h_work_L, h_work_R)
    ! tracer flux where the minimum BLD intersets layer
    ! GMM, khtr_avg should be computed once khtr is 3D
    F_layer(k_bot_min) = -(heff * khtr_u) * (phi_R_avg - phi_L_avg)

    do k = k_bot_min-1,1,-1
      heff = harmonic_mean(h_L(k), h_R(k))
      F_layer(k) = -(heff * khtr_u) * (phi_R(k) - phi_L(k))
    enddo
  endif

  if (boundary == BOTTOM) then
    k_top_max = MAX(k_top_L, k_top_R)
    ! make sure left and right k indices span same range
    if (k_top_max .ne. k_top_L) then
      k_top_L = k_top_max
      zeta_top_L = 1.0
    endif
    if (k_top_max .ne. k_top_R) then
      k_top_R= k_top_max
      zeta_top_R = 1.0
    endif

    h_work_L = (h_L(k_top_L) * zeta_top_L)
    h_work_R = (h_R(k_top_R) * zeta_top_R)

    phi_L_avg = average_value_ppoly( nk, phi_L, ppoly0_E_L, ppoly0_coefs_L, method, k_top_L, 1.0-zeta_top_L, 1.0)
    phi_R_avg = average_value_ppoly( nk, phi_R, ppoly0_E_R, ppoly0_coefs_R, method, k_top_R, 1.0-zeta_top_R, 1.0)
    heff = harmonic_mean(h_work_L, h_work_R)

    ! tracer flux where the minimum BLD intersets layer
    F_layer(k_top_max) = (-heff * khtr_u) * (phi_R_avg - phi_L_avg)

    do k = k_top_max+1,nk
      heff = harmonic_mean(h_L(k), h_R(k))
      F_layer(k) = -(heff * khtr_u) * (phi_R(k) - phi_L(k))
    enddo
  endif

end subroutine fluxes_layer_method

!> Apply the lateral boundary diffusive fluxes calculated from a 'bulk model'
!! See \ref section_method1
subroutine fluxes_bulk_method(boundary, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, ppoly0_coefs_L, &
                                   ppoly0_coefs_R, ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer, F_limit)

  integer,                   intent(in   )       :: boundary !< Which boundary layer SURFACE or BOTTOM  [nondim]
  integer,                   intent(in   )       :: nk       !< Number of layers                        [nondim]
  integer,                   intent(in   )       :: deg      !< order of the polynomial reconstruction  [nondim]
  real, dimension(nk),       intent(in   )       :: h_L      !< Layer thickness (left)          [H ~> m or kg m-2]
  real, dimension(nk),       intent(in   )       :: h_R      !< Layer thickness (right)         [H ~> m or kg m-2]
  real,                      intent(in   )       :: hbl_L    !< Thickness of the boundary boundary
                                                                       !! layer (left)          [H ~> m or kg m-2]
  real,                      intent(in   )       :: hbl_R    !< Thickness of the boundary boundary
                                                             !! layer (left)                    [H ~> m or kg m-2]
  real,                      intent(in   )       :: area_L   !< Area of the horizontal grid (left)      [L2 ~> m2]
  real,                      intent(in   )       :: area_R   !< Area of the horizontal grid (right)     [L2 ~> m2]
  real, dimension(nk),       intent(in   )       :: phi_L    !< Tracer values (left)                    [conc]
  real, dimension(nk),       intent(in   )       :: phi_R    !< Tracer values (right)                   [conc]
  real, dimension(nk,deg+1), intent(in   )       :: ppoly0_coefs_L !< Tracer reconstruction (left)      [conc]
  real, dimension(nk,deg+1), intent(in   )       :: ppoly0_coefs_R !< Tracer reconstruction (right)     [conc]
  real, dimension(nk,2),     intent(in   )       :: ppoly0_E_L !< Polynomial edge values (left)         [nondim]
  real, dimension(nk,2),     intent(in   )       :: ppoly0_E_R !< Polynomial edge values (right)        [nondim]
  integer,                   intent(in   )       :: method   !< Method of polynomial integration        [nondim]
  real,                      intent(in   )       :: khtr_u   !< Horizontal diffusivities times delta t
                                                             !! at a velocity point [L2 ~> m2]
  real,                      intent(  out)       :: F_bulk   !< The bulk mixed layer lateral flux
                                                             !! [H L2 conc ~> m3 conc]
  real, dimension(nk),       intent(  out)       :: F_layer  !< Layerwise diffusive flux at U- or V-point
                                                             !! [H L2 conc ~> m3 conc]
  logical, optional,         intent(in   )       :: F_limit  !< If True, apply a limiter

  ! Local variables
  real, dimension(nk) :: h_means              !< Calculate the layer-wise harmonic means           [H ~> m or kg m-2]
  real                :: khtr_avg             !< Thickness-weighted diffusivity at the u-point     [m^2 s^-1]
                                              !! This is just to remind developers that khtr_avg should be
                                              !! computed once khtr is 3D.
  real                :: heff                 !< Harmonic mean of layer thicknesses                [H ~> m or kg m-2]
  real                :: inv_heff             !< Inverse of the harmonic mean of layer thicknesses
                                              !! [H-1 ~> m-1 or m2 kg-1]
  real                :: phi_L_avg, phi_R_avg !< Bulk, thickness-weighted tracer averages (left and right column)
                                              !!                                                   [conc m^-3 ]
  real    :: htot                             !< Total column thickness [H ~> m or kg m-2]
  integer :: k, k_min, k_max                  !< k-indices, min and max for top and bottom, respectively
  integer :: k_top_L, k_bot_L                 !< k-indices left
  integer :: k_top_R, k_bot_R                 !< k-indices right
  real    :: zeta_top_L, zeta_top_R           !< distance from the top of a layer to the
                                              !! boundary layer   [nondim]
  real    :: zeta_bot_L, zeta_bot_R           !< distance from the bottom of a layer to the
                                              !! boundary layer   [nondim]
  real    :: h_work_L, h_work_R               !< dummy variables
  real    :: F_max                            !< The maximum amount of flux that can leave a
                                              !! cell  [m^3 conc]
  logical :: limiter                          !< True if flux limiter should be applied
  real    :: hfrac                            !< Layer fraction wrt sum of all layers [nondim]
  real    :: dphi                             !< tracer gradient                      [conc m^-3]

  if (hbl_L == 0. .or. hbl_R == 0.) then
    F_bulk = 0.
    F_layer(:) = 0.
    return
  endif

  limiter = .false.
  if (PRESENT(F_limit)) then
    limiter = F_limit
  endif

  ! Calculate vertical indices containing the boundary layer
  call boundary_k_range(boundary, nk, h_L, hbl_L, k_top_L, zeta_top_L, k_bot_L, zeta_bot_L)
  call boundary_k_range(boundary, nk, h_R, hbl_R, k_top_R, zeta_top_R, k_bot_R, zeta_bot_R)

  ! Calculate bulk averages of various quantities
  phi_L_avg  = bulk_average(boundary, nk, deg, h_L, hbl_L, phi_L, ppoly0_E_L, ppoly0_coefs_L, method, k_top_L, &
                            zeta_top_L, k_bot_L, zeta_bot_L)
  phi_R_avg  = bulk_average(boundary, nk, deg, h_R, hbl_R, phi_R, ppoly0_E_R, ppoly0_coefs_R, method, k_top_R, &
                            zeta_top_R, k_bot_R, zeta_bot_R)

  ! Calculate the 'bulk' diffusive flux from the bulk averaged quantities
  ! GMM, khtr_avg should be computed once khtr is 3D
  heff = harmonic_mean(hbl_L, hbl_R)
  F_bulk = -(khtr_u * heff) * (phi_R_avg - phi_L_avg)
  ! Calculate the layerwise sum of the vertical effective thickness. This is different than the heff calculated
  ! above, but is used as a way to decompose the fluxes onto the individual layers
  h_means(:) = 0.

  if (boundary == SURFACE) then
    k_min = MIN(k_bot_L, k_bot_R)

    ! left hand side
    if (k_bot_L == k_min) then
      h_work_L = h_L(k_min) * zeta_bot_L
    else
      h_work_L = h_L(k_min)
    endif

    ! right hand side
    if (k_bot_R == k_min) then
      h_work_R = h_R(k_min) * zeta_bot_R
    else
      h_work_R = h_R(k_min)
    endif

    h_means(k_min) = harmonic_mean(h_work_L,h_work_R)

    do k=1,k_min-1
      h_means(k) = harmonic_mean(h_L(k),h_R(k))
    enddo

  elseif (boundary == BOTTOM) then
    k_max = MAX(k_top_L, k_top_R)
    ! left hand side
    if (k_top_L == k_max) then
      h_work_L = h_L(k_max) * zeta_top_L
    else
      h_work_L = h_L(k_max)
    endif

    ! right hand side
    if (k_top_R == k_max) then
      h_work_R = h_R(k_max) * zeta_top_R
    else
      h_work_R = h_R(k_max)
    endif

    h_means(k_max) = harmonic_mean(h_work_L,h_work_R)

    do k=nk,k_max+1,-1
      h_means(k) = harmonic_mean(h_L(k),h_R(k))
    enddo
  endif

  if ( SUM(h_means) == 0. ) then
    return
 ! Decompose the bulk flux onto the individual layers
  else
    ! Initialize remaining thickness
    inv_heff = 1./SUM(h_means)
    do k=1,nk
      if (h_means(k) > 0.) then
        hfrac = h_means(k)*inv_heff
        F_layer(k) = F_bulk * hfrac

        if (limiter) then
          ! limit the flux to 0.2 of the tracer *gradient*
          ! Why 0.2?
          !  t=0         t=inf
          !   0           .2
          ! 0 1 0       .2.2.2
          !   0           .2
          !
          F_max = -0.2 * ((area_R*(phi_R(k)*h_R(k)))-(area_L*(phi_L(k)*h_R(k))))

          ! check if bulk flux (or F_layer) and F_max have same direction
          if ( SIGN(1.,F_bulk) == SIGN(1., F_max)) then
            ! Apply flux limiter calculated above
            if (F_max >= 0.) then
              F_layer(k) = MIN(F_layer(k),F_max)
            else
              F_layer(k) = MAX(F_layer(k),F_max)
            endif
          else
            ! do not apply a flux on this layer
            F_layer(k) = 0.
          endif
        else
          dphi = -(phi_R(k) - phi_L(k))
          if (.not. SIGN(1.,F_bulk) == SIGN(1., dphi)) then
            ! upgradient, do not apply a flux on this layer
            F_layer(k) = 0.
          endif
        endif ! limited
      endif
    enddo
  endif

end subroutine fluxes_bulk_method

!> Unit tests for near-boundary horizontal mixing
logical function near_boundary_unit_tests( verbose )
  logical,               intent(in) :: verbose !< If true, output additional information for debugging unit tests

  ! Local variables
  integer, parameter    :: nk = 2               ! Number of layers
  integer, parameter    :: deg = 1              ! Degree of reconstruction (linear here)
  integer, parameter    :: method = 1           ! Method used for integrating polynomials
  real, dimension(nk)   :: phi_L, phi_R         ! Tracer values (left and right column)             [ nondim m^-3 ]
  real, dimension(nk)   :: phi_L_avg, phi_R_avg ! Bulk, thickness-weighted tracer averages (left and right column)
  real, dimension(nk,deg+1) :: phi_pp_L, phi_pp_R   ! Coefficients for the linear pseudo-reconstructions
                                                !                                                   [ nondim m^-3 ]

  real, dimension(nk,2) :: ppoly0_E_L, ppoly0_E_R! Polynomial edge values (left and right)          [concentration]
  real, dimension(nk)   :: h_L, h_R             ! Layer thickness (left and right)                  [m]
  real                  :: khtr_u               ! Horizontal diffusivities at U-point               [m^2 s^-1]
  real                  :: hbl_L, hbl_R         ! Depth of the boundary layer (left and right)      [m]
  real                  :: F_bulk               ! Total diffusive flux across the U point           [nondim s^-1]
  real, dimension(nk)   :: F_layer              ! Diffusive flux within each layer at U-point       [nondim s^-1]
  real                  :: h_u, hblt_u          ! Thickness at the u-point                          [m]
  real                  :: khtr_avg             ! Thickness-weighted diffusivity at the u-point     [m^2 s^-1]
  real                  :: heff                 ! Harmonic mean of layer thicknesses                [m]
  real                  :: inv_heff             ! Inverse of the harmonic mean of layer thicknesses [m^[-1]
  character(len=120)    :: test_name            ! Title of the unit test
  integer               :: k_top                ! Index of cell containing top of boundary
  real                  :: zeta_top             ! Nondimension position
  integer               :: k_bot                ! Index of cell containing bottom of boundary
  real                  :: zeta_bot             ! Nondimension position
  real                  :: area_L,area_R        ! Area of grid cell [m^2]
  area_L = 1.; area_R = 1. ! Set to unity for all unit tests

  near_boundary_unit_tests = .false.

  ! Unit tests for boundary_k_range
  test_name = 'Surface boundary spans the entire top cell'
  h_L = (/5.,5./)
  call boundary_k_range(SURFACE, nk, h_L, 5., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 1, 1., test_name, verbose)

  test_name = 'Surface boundary spans the entire column'
  h_L = (/5.,5./)
  call boundary_k_range(SURFACE, nk, h_L, 10., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 2, 1., test_name, verbose)

  test_name = 'Bottom boundary spans the entire bottom cell'
  h_L = (/5.,5./)
  call boundary_k_range(BOTTOM, nk, h_L, 5., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 2, 1., 2, 0., test_name, verbose)

  test_name = 'Bottom boundary spans the entire column'
  h_L = (/5.,5./)
  call boundary_k_range(BOTTOM, nk, h_L, 10., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 1., 2, 0., test_name, verbose)

  test_name = 'Surface boundary intersects second layer'
  h_L = (/10.,10./)
  call boundary_k_range(SURFACE, nk, h_L, 17.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 2, 0.75, test_name, verbose)

  test_name = 'Surface boundary intersects first layer'
  h_L = (/10.,10./)
  call boundary_k_range(SURFACE, nk, h_L, 2.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 1, 0.25, test_name, verbose)

  test_name = 'Surface boundary is deeper than column thickness'
  h_L = (/10.,10./)
  call boundary_k_range(SURFACE, nk, h_L, 21.0, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 2, 1., test_name, verbose)

  test_name = 'Bottom boundary intersects first layer'
  h_L = (/10.,10./)
  call boundary_k_range(BOTTOM, nk, h_L, 17.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0.75, 2, 0., test_name, verbose)

  test_name = 'Bottom boundary intersects second layer'
  h_L = (/10.,10./)
  call boundary_k_range(BOTTOM, nk, h_L, 2.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 2, 0.25, 2, 0., test_name, verbose)

  ! All cases in this section have hbl which are equal to the column thicknesses
  test_name = 'Equal hbl and same layer thicknesses (gradient from right to left)'
  hbl_L = 10; hbl_R = 10
  h_L = (/5.,5./) ; h_R = (/5.,5./)
  phi_L = (/0.,0./) ; phi_R = (/1.,1./)
  phi_pp_L(1,1) = 0.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 1.; phi_pp_R(1,2) = 0.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 0.
  ppoly0_E_L(1,1) = 0.; ppoly0_E_L(1,2) = 0.
  ppoly0_E_L(2,1) = 0.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 1.; ppoly0_E_R(1,2) = 1.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 1.
  khtr_u = 1.
  ! Without limiter
  call fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, phi_pp_R, &
                                    ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, F_layer, (/-5.0,-5.0/) )

  ! same as above, but with limiter
  call fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, phi_pp_R, &
                                    ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer, .true.)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, F_layer, (/-1.0,-1.0/) )

  test_name = 'Equal hbl and same layer thicknesses (gradient from left to right)'
  hbl_L = 10.; hbl_R = 10.
  h_L = (/5.,5./) ; h_R = (/5.,5./)
  phi_L = (/1.,1./) ; phi_R = (/0.,0./)
  phi_pp_L(1,1) = 1.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 1.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 0.; phi_pp_R(1,2) = 0.
  phi_pp_R(2,1) = 0.; phi_pp_R(2,2) = 0.
  ppoly0_E_L(1,1) = 1.; ppoly0_E_L(1,2) = 1.
  ppoly0_E_L(2,1) = 1.; ppoly0_E_L(2,2) = 1.
  ppoly0_E_R(1,1) = 0.; ppoly0_E_R(1,2) = 0.
  ppoly0_E_R(2,1) = 0.; ppoly0_E_R(2,2) = 0.
  khtr_u = 1.
  call fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
                                    ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, F_layer, (/5.0,5.0/) )

  test_name = 'Equal hbl and same layer thicknesses (no gradient)'
  hbl_L = 10; hbl_R = 10
  h_L = (/5.,5./) ; h_R = (/5.,5./)
  phi_L = (/1.,1./) ; phi_R = (/1.,1./)
  phi_pp_L(1,1) = 1.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 1.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 1.; phi_pp_R(1,2) = 0.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 0.
  ppoly0_E_L(1,1) = 1.; ppoly0_E_L(1,2) = 1.
  ppoly0_E_L(2,1) = 1.; ppoly0_E_L(2,2) = 1.
  ppoly0_E_R(1,1) = 1.; ppoly0_E_R(1,2) = 1.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 1.
  khtr_u = 1.
  call fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
                                    ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, F_layer, (/0.0,0.0/) )

  test_name = 'Equal hbl and different layer thicknesses (gradient right to left)'
  hbl_L = 16.; hbl_R = 16.
  h_L = (/10.,6./) ; h_R = (/6.,10./)
  phi_L = (/0.,0./) ; phi_R = (/1.,1./)
  phi_pp_L(1,1) = 0.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 1.; phi_pp_R(1,2) = 0.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 0.
  ppoly0_E_L(1,1) = 0.; ppoly0_E_L(1,2) = 0.
  ppoly0_E_L(2,1) = 0.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 1.; ppoly0_E_R(1,2) = 1.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 1.
  khtr_u = 1.
  call fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
                                    ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, F_layer, (/-8.0,-8.0/) )

  test_name = 'Equal hbl and same layer thicknesses (diagonal tracer values)'
  hbl_L = 10.; hbl_R = 10.
  h_L = (/5.,5./) ; h_R = (/5.,5./)
  phi_L = (/1.,0./) ; phi_R = (/0.,1./)
  phi_pp_L(1,1) = 1.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 0.; phi_pp_R(1,2) = 0.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 0.
  ppoly0_E_L(1,1) = 1.; ppoly0_E_L(1,2) = 1.
  ppoly0_E_L(2,1) = 0.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 0.; ppoly0_E_R(1,2) = 0.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 1.
  khtr_u = 1.
  call fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
                                    ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, F_layer, (/0.0,0.0/) )

  test_name = 'Different hbl and different layer thicknesses (gradient from right to left)'
  hbl_L = 12; hbl_R = 20
  h_L = (/6.,6./) ; h_R = (/10.,10./)
  phi_L = (/0.,0./) ; phi_R = (/1.,1./)
  phi_pp_L(1,1) = 0.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 1.; phi_pp_R(1,2) = 0.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 0.
  ppoly0_E_L(1,1) = 0.; ppoly0_E_L(1,2) = 0.
  ppoly0_E_L(2,1) = 0.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 1.; ppoly0_E_R(1,2) = 1.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 1.
  khtr_u = 1.
  call fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
                                    ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, F_layer, (/-7.5,-7.5/) )

  ! Cases where hbl < column thickness (polynomial coefficients specified for pseudo-linear reconstruction)

  test_name = 'hbl < column thickness, hbl same, constant concentration each column'
  hbl_L = 2; hbl_R = 2
  h_L = (/1.,2./) ; h_R = (/1.,2./)
  phi_L = (/0.,0./) ; phi_R = (/1.,1./)
  phi_pp_L(1,1) = 0.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 1.; phi_pp_R(1,2) = 0.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 0.
  ppoly0_E_L(1,1) = 0.; ppoly0_E_L(1,2) = 0.
  ppoly0_E_L(2,1) = 0.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 1.; ppoly0_E_R(1,2) = 1.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 1.
  khtr_u = 1.
  call fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
                                    ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, F_layer, (/-1.,-1./) )

  test_name = 'hbl < column thickness, hbl same, linear profile right'
  hbl_L = 2; hbl_R = 2
  h_L = (/1.,2./) ; h_R = (/1.,2./)
  phi_L = (/0.,0./) ; phi_R = (/0.5,2./)
  phi_pp_L(1,1) = 0.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 0.; phi_pp_R(1,2) = 1.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 2.
  khtr_u = 1.
  ppoly0_E_L(1,1) = 0.; ppoly0_E_L(1,2) = 0.
  ppoly0_E_L(2,1) = 0.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 0.; ppoly0_E_R(1,2) = 1.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 3.
  call fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
                                    ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, F_layer, (/-1.,-1./) )

  test_name = 'hbl < column thickness, hbl same, linear profile right, khtr=2'
  hbl_L = 2; hbl_R = 2
  h_L = (/1.,2./) ; h_R = (/1.,2./)
  phi_L = (/0.,0./) ; phi_R = (/0.5,2./)
  phi_pp_L(1,1) = 0.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 0.; phi_pp_R(1,2) = 1.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 2.
  khtr_u = 2.
  ppoly0_E_L(1,1) = 0.; ppoly0_E_L(1,2) = 0.
  ppoly0_E_L(2,1) = 0.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 0.; ppoly0_E_R(1,2) = 1.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 3.
  call fluxes_layer_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, &
                                    phi_pp_R, ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_layer)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, F_layer, (/-1.,-3./) )

  ! unit tests for layer by layer method
  test_name = 'Different hbl and different column thicknesses (gradient from right to left)'
  hbl_L = 12; hbl_R = 20
  h_L = (/6.,6./) ; h_R = (/10.,10./)
  phi_L = (/0.,0./) ; phi_R = (/1.,1./)
  phi_pp_L(1,1) = 0.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 1.; phi_pp_R(1,2) = 0.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 0.
  ppoly0_E_L(1,1) = 0.; ppoly0_E_L(1,2) = 0.
  ppoly0_E_L(2,1) = 0.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 1.; ppoly0_E_R(1,2) = 1.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 1.
  khtr_u = 1.
  call fluxes_layer_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, &
                                    phi_pp_R, ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_layer)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, F_layer, (/-7.5,-7.5/) )

  test_name = 'Different hbl and different column thicknesses (linear profile right)'

  hbl_L = 15; hbl_R = 6
  h_L = (/10.,10./) ; h_R = (/12.,10./)
  phi_L = (/0.,0./) ; phi_R = (/1.,3./)
  phi_pp_L(1,1) = 0.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 0.; phi_pp_R(1,2) = 2.
  phi_pp_R(2,1) = 2.; phi_pp_R(2,2) = 2.
  ppoly0_E_L(1,1) = 0.; ppoly0_E_L(1,2) = 0.
  ppoly0_E_L(2,1) = 0.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 0.; ppoly0_E_R(1,2) = 2.
  ppoly0_E_R(2,1) = 2.; ppoly0_E_R(2,2) = 4.
  khtr_u = 1.
  call fluxes_layer_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, area_L, area_R, phi_L, phi_R, phi_pp_L, &
                                    phi_pp_R, ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_layer)
  near_boundary_unit_tests = near_boundary_unit_tests .or. &
                             test_layer_fluxes( verbose, nk, test_name, F_layer, (/-3.75,0.0/) )
end function near_boundary_unit_tests

!> Returns true if output of near-boundary unit tests does not match correct computed values
!! and conditionally writes results to stream
logical function test_layer_fluxes(verbose, nk, test_name, F_calc, F_ans)
  logical,                    intent(in) :: verbose   !< If true, write results to stdout
  character(len=80),          intent(in) :: test_name !< Brief description of the unit test
  integer,                    intent(in) :: nk        !< Number of layers
  real, dimension(nk),        intent(in) :: F_calc    !< Fluxes of the unitless tracer from the algorithm [s^-1]
  real, dimension(nk),        intent(in) :: F_ans     !< Fluxes of the unitless tracer calculated by hand [s^-1]
  ! Local variables
  integer :: k
  integer, parameter :: stdunit = stdout

  test_layer_fluxes = .false.
  do k=1,nk
    if ( F_calc(k) /= F_ans(k) ) then
      test_layer_fluxes = .true.
      write(stdunit,*) "MOM_lateral_boundary_diffusion, UNIT TEST FAILED: ", test_name
      write(stdunit,10) k, F_calc(k), F_ans(k)
      ! ### Once these unit tests are passing, and failures are caught properly,
      ! we will post failure notifications to both stdout and stderr.
      !write(stderr,*) "MOM_lateral_boundary_diffusion, UNIT TEST FAILED: ", test_name
      !write(stderr,10) k, F_calc(k), F_ans(k)
    elseif (verbose) then
      write(stdunit,10) k, F_calc(k), F_ans(k)
    endif
  enddo

10 format("Layer=",i3," F_calc=",f20.16," F_ans",f20.16)
end function test_layer_fluxes

!> Return true if output of unit tests for boundary_k_range does not match answers
logical function test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, k_top_ans, zeta_top_ans,&
                                       k_bot_ans, zeta_bot_ans, test_name, verbose)
  integer :: k_top               !< Index of cell containing top of boundary
  real    :: zeta_top            !< Nondimension position
  integer :: k_bot               !< Index of cell containing bottom of boundary
  real    :: zeta_bot            !< Nondimension position
  integer :: k_top_ans           !< Index of cell containing top of boundary
  real    :: zeta_top_ans        !< Nondimension position
  integer :: k_bot_ans           !< Index of cell containing bottom of boundary
  real    :: zeta_bot_ans        !< Nondimension position
  character(len=80) :: test_name !< Name of the unit test
  logical :: verbose             !< If true always print output

  integer, parameter :: stdunit = stdout

  test_boundary_k_range = k_top .ne. k_top_ans
  test_boundary_k_range = test_boundary_k_range .or. (zeta_top .ne. zeta_top_ans)
  test_boundary_k_range = test_boundary_k_range .or. (k_bot .ne. k_bot_ans)
  test_boundary_k_range = test_boundary_k_range .or. (zeta_bot .ne. zeta_bot_ans)

  if (test_boundary_k_range) write(stdunit,*) "UNIT TEST FAILED: ", test_name
  if (test_boundary_k_range .or. verbose) then
    write(stdunit,20) "k_top", k_top, "k_top_ans", k_top_ans
    write(stdunit,20) "k_bot", k_bot, "k_bot_ans", k_bot_ans
    write(stdunit,30) "zeta_top", zeta_top, "zeta_top_ans", zeta_top_ans
    write(stdunit,30) "zeta_bot", zeta_bot, "zeta_bot_ans", zeta_bot_ans
  endif

  20 format(A,"=",i3,X,A,"=",i3)
  30 format(A,"=",f20.16,X,A,"=",f20.16)


end function test_boundary_k_range

!> \namespace mom_lateral_boundary_diffusion
!!
!! \section section_LBD The Lateral Boundary Diffusion (LBD) framework
!!
!! The LBD framework accounts for the effects of diabatic mesoscale fluxes
!! within surface and bottom boundary layers. Unlike the equivalent adiabatic
!! fluxes, which is applied along neutral density surfaces, LBD is purely
!! horizontal.
!!
!! The bottom boundary layer fluxes remain to be implemented, although most
!! of the steps needed to do so have already been added and tested.
!!
!! Boundary lateral diffusion can be applied using one of the three methods:
!!
!! * [Method #1: Bulk layer](@ref section_method1) (default);
!! * [Method #2: Along layer](@ref section_method2);
!!
!! A brief summary of these methods is provided below.
!!
!! \subsection section_method1 Bulk layer approach (Method #1)
!!
!! Apply the lateral boundary diffusive fluxes calculated from a 'bulk model'.This
!! is a lower order representation (Kraus-Turner like approach) which assumes that
!! eddies are acting along well mixed layers (i.e., eddies do not know care about
!! vertical tracer gradients within the boundary layer).
!!
!! Step #1: compute vertical indices containing boundary layer (boundary_k_range).
!! For the TOP boundary layer, these are:
!!
!! k_top, k_bot, zeta_top, zeta_bot
!!
!! Step #2: compute bulk averages (thickness weighted) tracer averages (phi_L and phi_R),
!! then calculate the bulk diffusive flux (F_{bulk}):
!!
!! \f[ F_{bulk} = -KHTR \times h_{eff} \times (\phi_R - \phi_L),  \f]
!! where h_eff is the [harmonic mean](@ref section_harmonic_mean) of the boundary layer depth
!! in the left and right columns (\f[ HBL_L \f] and \f[ HBL_R \f], respectively).
!!
!! Step #3: decompose F_bulk onto individual layers:
!!
!! \f[ F_{layer}(k) = F_{bulk} \times h_{frac}(k) ,  \f]
!!
!! where h_{frac} is
!!
!! \f[ h_{frac}(k) = h_u(k) \times \frac{1}{\sum(h_u)}.  \f]
!!
!! h_u is the [harmonic mean](@ref section_harmonic_mean) of thicknesses at each layer.
!! Special care (layer reconstruction) must be taken at k_min = min(k_botL, k_bot_R).
!!
!! Step #4: limit the tracer flux so that 1) only down-gradient fluxes are applied,
!! and 2) the flux cannot be larger than F_max, which is defined using the tracer
!! gradient:
!!
!! \f[ F_{max} = -0.2 \times [(V_R(k) \times \phi_R(k)) - (V_L(k) \times \phi_L(k))],  \f]
!! where V is the cell volume. Why 0.2?
!!          t=0         t=inf
!!           0           .2
!!         0 1 0       .2.2.2
!!           0           .2
!!
!! \subsection section_method2 Along layer approach (Method #2)
!!
!! This is a more straight forward method where diffusion is applied layer by layer using
!! only information from neighboring cells.
!!
!! Step #1: compute vertical indices containing boundary layer (boundary_k_range).
!! For the TOP boundary layer, these are:
!!
!! k_top, k_bot, zeta_top, zeta_bot
!!
!! Step #2: calculate the diffusive flux at each layer:
!!
!! \f[ F_{k} = -KHTR \times h_{eff}(k) \times (\phi_R(k) - \phi_L(k)),  \f]
!! where h_eff is the [harmonic mean](@ref section_harmonic_mean) of the layer thickness
!! in the left and right columns. Special care (layer reconstruction) must be taken at
!! k_min = min(k_botL, k_bot_R). This method does not require a limiter since KHTR
!! is already limted based on a diffusive CFL condition prior to the call of this
!! module.
!!
!! \subsection section_harmonic_mean Harmonic Mean
!!
!! The harmonic mean (HM) betwen h1 and h2 is defined as:
!!
!! \f[ HM = \frac{2 \times h1 \times h2}{h1 + h2} \f]
!!
end module MOM_lateral_boundary_diffusion
