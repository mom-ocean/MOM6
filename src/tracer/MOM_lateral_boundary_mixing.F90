!> Calculate and apply diffusive fluxes as a parameterization of lateral mixing (non-neutral) by
!! mesoscale eddies near the top and bottom boundary layers of the ocean.
module MOM_lateral_boundary_mixing

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

implicit none ; private

public near_boundary_unit_tests, lateral_boundary_mixing, lateral_boundary_mixing_init

! Private parameters to avoid doing string comparisons for bottom or top boundary layer
integer, parameter :: SURFACE = -1 !< Set a value that corresponds to the surface bopundary
integer, parameter :: BOTTOM  = 1  !< Set a value that corresponds to the bottom boundary
#include <MOM_memory.h>

type, public :: lateral_boundary_mixing_CS ; private
  integer :: method                                               !< Determine which of the three methods calculate
                                                                  !! and apply near boundary layer fluxes
                                                                  !! 1. bulk-layer approach
                                                                  !! 2. Along layer
                                                                  !! 3. Decomposition onto pressure levels
  integer :: deg                                                  !< Degree of polynomial reconstruction
  integer :: surface_boundary_scheme                              !< Which boundary layer scheme to use
                                                                  !! 1. ePBL; 2. KPP
  type(remapping_CS)              :: remap_CS                     !< Control structure to hold remapping configuration
  type(KPP_CS),           pointer :: KPP_CSp => NULL()            !< KPP control structure needed to get BLD
  type(energetic_PBL_CS), pointer :: energetic_PBL_CSp => NULL()  !< ePBL control structure needed to get MLD
  type(diag_ctrl), pointer :: diag => NULL()                      !< A structure that is used to
                                                                  !! regulate the timing of diagnostic output.
end type lateral_boundary_mixing_CS

! This include declares and sets the variable "version".
#include "version_variable.h"
character(len=40) :: mdl = "MOM_lateral_boundary_mixing"

contains

!> Initialization routine that reads runtime parameters and sets up pointers to other control structures that might be
!! needed for lateral boundary mixing
logical function lateral_boundary_mixing_init(Time, G, param_file, diag, diabatic_CSp, CS)
  type(time_type), target,          intent(in)    :: Time       !< Time structure
  type(ocean_grid_type),            intent(in)    :: G          !< Grid structure
  type(param_file_type),            intent(in)    :: param_file !< Parameter file structure
  type(diag_ctrl), target,          intent(inout) :: diag       !< Diagnostics control structure
  type(diabatic_CS),                pointer       :: diabatic_CSp !< KPP control structure needed to get BLD
  type(lateral_boundary_mixing_CS), pointer       :: CS         !< Lateral boundary mixing control structure

  character(len=80)  :: string  ! Temporary strings
  logical :: boundary_extrap

  if (ASSOCIATED(CS)) then
    call MOM_error(FATAL, "lateral_boundary_mixing_init called with associated control structure.")
    return
  endif

  ! Log this module and master switch for turning it on/off
  call log_version(param_file, mdl, version, &
       "This module implements lateral boundary  mixing of tracers")
  call get_param(param_file, mdl, "USE_LATERAL_BOUNDARY_MIXING", lateral_boundary_mixing_init, &
                 "If true, enables the lateral boundary mixing module.", &
                 default=.false.)

  if (.not. lateral_boundary_mixing_init) then
    return
  endif

  allocate(CS)
  CS%diag => diag
  call extract_diabatic_member(diabatic_CSp, KPP_CSp=CS%KPP_CSp)
  call extract_diabatic_member(diabatic_CSp, energetic_PBL_CSp=CS%energetic_PBL_CSp)

  CS%surface_boundary_scheme = -1
  if ( .not. ASSOCIATED(CS%energetic_PBL_CSp) .and. .not. ASSOCIATED(CS%KPP_CSp) ) then
    call MOM_error(FATAL,"Lateral boundary mixing is true, but no valid boundary layer scheme was found")
  endif

  ! Read all relevant parameters and write them to the model log.
  call get_param(param_file, mdl, "LATERAL_BOUNDARY_METHOD", CS%method, &
                 "Determine how to apply near-boundary lateral mixing of tracers"//&
                 "1. Bulk layer approach"//&
                 "2. Along layer approach"//&
                 "3. Decomposition on to pressure levels", default=1)
  call get_param(param_file, mdl, "LBM_BOUNDARY_EXTRAP", boundary_extrap, &
                 "Use boundary extrapolation in LBM code", &
                 default=.false.)
  call get_param(param_file, mdl, "LBM_REMAPPING_SCHEME", string, &
                 "This sets the reconstruction scheme used "//&
                 "for vertical remapping for all variables. "//&
                 "It can be one of the following schemes: "//&
                 trim(remappingSchemesDoc), default=remappingDefaultScheme)
  call initialize_remapping( CS%remap_CS, string, boundary_extrapolation = boundary_extrap )
  call extract_member_remapping_CS(CS%remap_CS, degree=CS%deg)

end function lateral_boundary_mixing_init

!> Driver routine for calculating lateral diffusive fluxes near the top and bottom boundaries. Two different methods
!! Method 1: Calculate fluxes from bulk layer integrated quantities
subroutine lateral_boundary_mixing(G, GV, US, h, Coef_x, Coef_y, dt, Reg, CS)
  type(ocean_grid_type),                intent(inout) :: G       !< Grid type
  type(verticalGrid_type),              intent(in)    :: GV      !< ocean vertical grid structure
  type(unit_scale_type),            intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                                        intent(in)    :: h       !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G)),    intent(in)    :: Coef_x !< dt * Kh * dy / dx at u-points [m2]
  real, dimension(SZI_(G),SZJB_(G)),    intent(in)    :: Coef_y !< dt * Kh * dx / dy at v-points [m2]
  real,                                 intent(in)    :: dt     !< Tracer time step * I_numitts
                                                                !! (I_numitts in tracer_hordiff)
  type(tracer_registry_type),           pointer       :: Reg    !< Tracer registry
  type(lateral_boundary_mixing_CS),     intent(in)    :: CS      !< Control structure for this module
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: hbl   !< bnd. layer depth [m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G),CS%deg+1) :: ppoly0_coefs !< Coefficients of polynomial
  real, dimension(SZI_(G),SZJ_(G),SZK_(G),2)        :: ppoly0_E     !< Edge values from reconstructions
  real, dimension(SZK_(G),CS%deg+1)                    :: ppoly_S      !< Slopes from reconstruction (placeholder)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: uFlx        ! Zonal flux of tracer [H conc ~> m conc or conc kg m-2]
  real, dimension(SZI_(G),SZJ_(G))          :: uFLx_bulk   ! Total calculated bulk-layer u-flux for the tracer
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: vFlx        ! Meridional flux of tracer
  real, dimension(SZI_(G),SZJB_(G))         :: vFlx_bulk   ! Total calculated bulk-layer v-flux for the tracer
  type(tracer_type), pointer                :: Tracer => NULL() ! Pointer to the current tracer
  integer :: remap_method !< Reconstruction method
  integer :: i,j,k,m

  hbl(:,:) = 0.
  if (ASSOCIATED(CS%KPP_CSp)) call KPP_get_BLD(CS%KPP_CSp, hbl, G)
  if (ASSOCIATED(CS%energetic_PBL_CSp)) call energetic_PBL_get_MLD(CS%energetic_PBL_CSp, hbl, G, US)

  call pass_var(hbl,G%Domain)

  do m = 1,Reg%ntr
    tracer => Reg%tr(m)
    do j = G%jsc-1, G%jec+1
      ! Interpolate state to interface
      do i = G%isc-1, G%iec+1
          call build_reconstructions_1d( CS%remap_CS, G%ke, h(i,j,:), tracer%t(i,j,:), ppoly0_coefs(i,j,:,:), &
                                         ppoly0_E(i,j,:,:), ppoly_S, remap_method, GV%H_subroundoff, GV%H_subroundoff)
      enddo
    enddo
    ! Diffusive fluxes in the i-direction
    uFlx(:,:,:) = 0.
    vFlx(:,:,:) = 0.
    uFlx_bulk(:,:) = 0.
    vFlx_bulk(:,:) = 0.
    if ( CS%method == 1 ) then
      do j=G%jsc,G%jec
        do i=G%isc-1,G%iec
          if (G%mask2dCu(I,j)>0.) then
            call layer_fluxes_bulk_method(SURFACE, GV%ke, CS%deg, h(i,j,:), h(i+1,j,:), hbl(i,j), hbl(i+1,j), &
              tracer%t(i,j,:), tracer%t(i+1,j,:), ppoly0_coefs(i,j,:,:), ppoly0_coefs(i+1,j,:,:), ppoly0_E(i,j,:,:), &
              ppoly0_E(i+1,j,:,:), remap_method, Coef_x(I,j), uFlx_bulk(I,j), uFlx(I,j,:))
          endif
        enddo
      enddo
      do J=G%jsc-1,G%jec
        do i=G%isc,G%iec
          if (G%mask2dCv(i,J)>0.) then
            call layer_fluxes_bulk_method(SURFACE, GV%ke, CS%deg, h(i,J,:), h(i,J+1,:), hbl(i,J), hbl(i,J+1), &
              tracer%t(i,J,:), tracer%t(i,J+1,:), ppoly0_coefs(i,J,:,:), ppoly0_coefs(i,J+1,:,:), ppoly0_E(i,J,:,:), &
              ppoly0_E(i,J+1,:,:), remap_method, Coef_y(i,J), vFlx_bulk(i,J), vFlx(i,J,:))
          endif
        enddo
      enddo
    endif

    ! Update the tracer fluxes
    do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if (G%mask2dT(i,j)>0.) then
        tracer%t(i,j,k) = tracer%t(i,j,k) + (( (uFlx(I-1,j,k)-uFlx(I,j,k)) ) + ( (vFlx(i,J-1,k)-vFlx(i,J,k) ) ))*(G%IareaT(i,j)/( h(i,j,k) + GV%H_subroundoff))
      endif
    enddo ; enddo ; enddo

    ! Post the tracer diagnostics
    if (tracer%id_lbm_bulk_dfx>0) call post_data(tracer%id_lbm_bulk_dfx, uFlx_bulk, CS%diag)
    if (tracer%id_lbm_bulk_dfy>0) call post_data(tracer%id_lbm_bulk_dfy, vFlx_bulk, CS%diag)
    if (tracer%id_lbm_dfx>0)      call post_data(tracer%id_lbm_dfx, uFlx, CS%diag)
    if (tracer%id_lbm_dfy>0)      call post_data(tracer%id_lbm_dfy, vFlx, CS%diag)

  enddo

end subroutine lateral_boundary_mixing

!< Calculate bulk layer value of a scalar quantity as the thickness weighted average
real function bulk_average(boundary, nk, deg, h, hBLT, phi, ppoly0_E, ppoly0_coefs, method, k_top, zeta_top, k_bot, zeta_bot)
  integer             :: boundary          !< SURFACE or BOTTOM                                         [nondim]
  integer             :: nk                !< Number of layers                                          [nondim]
  integer             :: deg               !< Degree of polynomial                                      [nondim]
  real, dimension(nk) :: h                 !< Layer thicknesses                                         [m]
  real                :: hBLT              !< Depth of the mixing layer                                 [m]
  real, dimension(nk) :: phi               !< Scalar quantity
  real, dimension(nk,2)    :: ppoly0_E(:,:)     !< Edge value of polynomial
  real, dimension(nk,deg+1) :: ppoly0_coefs(:,:) !< Coefficients of polynomial
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
  real    :: htot ! Running sum of the thicknesses (top to bottom)
  integer :: k


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
  real, dimension(nk), intent(in   ) :: h        !< Layer thicknesses of the coluymn        [m]
  real,                intent(in   ) :: hbl      !< Thickness of the boundary layer         [m]
                                                 !! If surface, with respect to zbl_ref = 0.
                                                 !! If bottom, with respect to zbl_ref = SUM(h)
  integer,             intent(  out) :: k_top    !< Index of the first layer within the boundary
  real,                intent(  out) :: zeta_top !< Distance from the top of a layer to the intersection of the
                                                 !! top extent of the boundary layer (0 at top, 1 at bottom)  [nondim]
  integer,             intent(  out) :: k_bot    !< Index of the last layer within the boundary
  real,                intent(  out) :: zeta_bot !< Distance of the lower layer to the boundary layer depth
                                                 !! (0 at top, 1 at bottom)  [nondim]
  ! Local variables
  real :: htot
  integer :: k
  ! Surface boundary layer
  if ( boundary == SURFACE ) then
    k_top = 1
    zeta_top = 0.
    htot = 0.
    k_bot = 1
    zeta_bot = 0.
    if (hbl == 0.) return
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
    zeta_bot = 1.
    htot = 0.
    if (hbl == 0.) return
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

!> Calculate the near-boundary diffusive fluxes calculated from a 'bulk model'
subroutine layer_fluxes_bulk_method(boundary, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R, ppoly0_coefs_L, &
                                   ppoly0_coefs_R, ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  integer,                   intent(in   )       :: boundary !< Which boundary layer SURFACE or BOTTOM  [nondim]
  integer,                   intent(in   )       :: nk       !< Number of layers                        [nondim]
  integer,                   intent(in   )       :: deg      !< order of the polynomial reconstruction  [nondim]
  real, dimension(nk),       intent(in   )       :: h_L      !< Layer thickness (left)                  [m]
  real, dimension(nk),       intent(in   )       :: h_R      !< Layer thickness (right)                 [m]
  real,                      intent(in   )       :: hbl_L    !< Thickness of the boundary boundary
                                                                       !! layer (left)                  [m]
  real,                      intent(in   )       :: hbl_R    !< Thickness of the boundary boundary
                                                             !! layer (left)                            [m]
  real, dimension(nk),       intent(in   )       :: phi_L    !< Tracer values (left)                    [ nondim m^-3 ]
  real, dimension(nk),       intent(in   )       :: phi_R    !< Tracer values (right)                   [ nondim m^-3 ]
  real, dimension(nk,deg+1), intent(in   )       :: ppoly0_coefs_L !< Tracer reconstruction (left)      [ nondim m^-3 ]
  real, dimension(nk,deg+1), intent(in   )       :: ppoly0_coefs_R !< Tracer reconstruction (right)     [ nondim m^-3 ]
  real, dimension(nk,2),     intent(in   )       :: ppoly0_E_L !< Polynomial edge values (left)         [ nondim ]
  real, dimension(nk,2),     intent(in   )       :: ppoly0_E_R !< Polynomial edge values (right)        [ nondim ]
  integer,                   intent(in   )       :: method   !< Method of polynomial integration        [ nondim ]
  real,                      intent(in   )       :: khtr_u   !< Horizontal diffusivities at U-point     [m^2 s^-1]
  real,                      intent(  out)       :: F_bulk   !< The bulk mixed layer lateral flux       [trunit s^-1]
  real, dimension(nk),       intent(  out)       :: F_layer  !< Layerwise diffusive flux at U-point     [trunit s^-1]
  ! Local variables
  real, dimension(nk) :: h_means              ! Calculate the layer-wise harmonic means           [m]
  real, dimension(nk) :: h_u                  ! Thickness at the u-point                          [m]
  real                :: hbl_u                ! Boundary layer Thickness at the u-point           [m]
  real                :: khtr_avg             ! Thickness-weighted diffusivity at the u-point     [m^2 s^-1]
  real                :: heff                 ! Harmonic mean of layer thicknesses                [m]
  real                :: inv_heff             ! Inverse of the harmonic mean of layer thicknesses [m^[-1]
  real                :: phi_L_avg, phi_R_avg ! Bulk, thickness-weighted tracer averages (left and right column)
                                              !                                                   [trunit m^-3 ]
  real    :: htot                 ! Total column thickness [m]
  integer :: k, k_min, k_max
  integer :: k_top_L, k_bot_L, k_top_u
  integer :: k_top_R, k_bot_R, k_bot_u
  real    :: zeta_top_L, zeta_top_R, zeta_top_u
  real    :: zeta_bot_L, zeta_bot_R, zeta_bot_u
  real    :: h_work_L, h_work_R  ! dummy variables
  if (hbl_L == 0. .or. hbl_R == 0.) then
    F_bulk = 0.
    F_layer(:) = 0.
    return
  endif
  ! Calculate vertical indices containing the boundary layer
  call boundary_k_range(boundary, nk, h_L, hbl_L, k_top_L, zeta_top_L, k_bot_L, zeta_bot_L)
  call boundary_k_range(boundary, nk, h_R, hbl_R, k_top_R, zeta_top_R, k_bot_R, zeta_bot_R)
  ! Calculate bulk averages of various quantities
  phi_L_avg  = bulk_average(boundary, nk, deg, h_L, hbl_L, phi_L, ppoly0_E_L, ppoly0_coefs_L, method, k_top_L, zeta_top_L,&
                            k_bot_L, zeta_bot_L)
  phi_R_avg  = bulk_average(boundary, nk, deg, h_R, hbl_R, phi_R, ppoly0_E_R, ppoly0_coefs_R, method, k_top_R, zeta_top_R,&
                            k_bot_R, zeta_bot_R)
  do k=1,nk
    h_u(k) = 0.5 * (h_L(k) + h_R(k))
  enddo
  hbl_u = 0.5*(hbl_L + hbl_R)
  call boundary_k_range(boundary, nk, h_u, hbl_u, k_top_u, zeta_top_u, k_bot_u, zeta_bot_u)

  ! Calculate the 'bulk' diffusive flux from the bulk averaged quantities
  heff = harmonic_mean(hbl_L, hbl_R)
  F_bulk = -(khtr_u * heff) * (phi_R_avg - phi_L_avg)
  if (F_bulk .ne. F_bulk) print *, khtr_avg, heff, phi_R_avg, phi_L_avg, hbl_L, hbl_R
  ! Calculate the layerwise sum of the vertical effective thickness. This is different than the heff calculated
  ! above, but is used as a way to decompose decompose the fluxes onto the individual layers
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
  endif

  if (boundary == BOTTOM) then
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
  else
    inv_heff = 1./SUM(h_means)
    ! Decompose the bulk flux onto the individual layers
    do k=1,nk
      if ( SIGN(1.,F_bulk) == SIGN(1., -(phi_R(k)-phi_L(k))) ) then
        F_layer(k) = F_bulk * (h_means(k)*inv_heff)
      else
        F_layer(k) = 0.
      endif
    enddo
  endif

end subroutine layer_fluxes_bulk_method

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
  real                  :: hbl_L, hbl_R       ! Depth of the boundary layer (left and right)      [m]
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
  near_boundary_unit_tests = .false.

  ! Unit tests for boundary_k_range
  test_name = 'Surface boundary spans the entire top cell'
  h_L = (/5.,5./)
  call boundary_k_range(SURFACE, nk, h_L, 5., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 1, 1., test_name, verbose)

  test_name = 'Surface boundary spans the entire column'
  h_L = (/5.,5./)
  call boundary_k_range(SURFACE, nk, h_L, 10., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 2, 1., test_name, verbose)

  test_name = 'Bottom boundary spans the entire bottom cell'
  h_L = (/5.,5./)
  call boundary_k_range(BOTTOM, nk, h_L, 5., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 2, 0., 2, 1., test_name, verbose)

  test_name = 'Bottom boundary spans the entire column'
  h_L = (/5.,5./)
  call boundary_k_range(BOTTOM, nk, h_L, 10., k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 2, 1., test_name, verbose)

  test_name = 'Surface boundary intersects second layer'
  h_L = (/10.,10./)
  call boundary_k_range(SURFACE, nk, h_L, 17.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 2, 0.75, test_name, verbose)

  test_name = 'Surface boundary intersects first layer'
  h_L = (/10.,10./)
  call boundary_k_range(SURFACE, nk, h_L, 2.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0., 1, 0.25, test_name, verbose)

  test_name = 'Bottom boundary intersects first layer'
  h_L = (/10.,10./)
  call boundary_k_range(BOTTOM, nk, h_L, 17.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 1, 0.75, 2, 1., test_name, verbose)

  test_name = 'Bottom boundary intersects second layer'
  h_L = (/10.,10./)
  call boundary_k_range(BOTTOM, nk, h_L, 2.5, k_top, zeta_top, k_bot, zeta_bot)
  near_boundary_unit_tests = test_boundary_k_range(k_top, zeta_top, k_bot, zeta_bot, 2, 0.25, 2, 1., test_name, verbose)

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
  call layer_fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
                                    ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  near_boundary_unit_tests = test_layer_fluxes( verbose, nk, test_name, F_layer, (/-5.0,-5.0/) )

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
  call layer_fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
                                    ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  near_boundary_unit_tests = test_layer_fluxes( verbose, nk, test_name, F_layer, (/5.0,5.0/) )

  test_name = 'Equal hbl and same layer thicknesses (no gradient)'
  hbl_L = 10; hbl_R = 10
  h_L = (/5.,5./) ; h_R = (/5.,5./)
  phi_L = (/1.,1./) ; phi_R = (/1.,1./)
  phi_pp_L(1,1) = 1.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 1.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 1.; phi_pp_R(1,2) = 0.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 0.
  ppoly0_E_L(1,1) = 1.; ppoly0_E_L(1,2) = 0.
  ppoly0_E_L(2,1) = 1.; ppoly0_E_L(2,2) = 0.
  ppoly0_E_R(1,1) = 1.; ppoly0_E_R(1,2) = 1.
  ppoly0_E_R(2,1) = 1.; ppoly0_E_R(2,2) = 1.
  khtr_u = 1.
  call layer_fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
                                    ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  near_boundary_unit_tests = test_layer_fluxes( verbose, nk, test_name, F_layer, (/0.0,0.0/) )

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
  call layer_fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
                                    ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  near_boundary_unit_tests = test_layer_fluxes( verbose, nk, test_name, F_layer, (/-8.0,-8.0/) )

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
  call layer_fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
                                    ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  near_boundary_unit_tests = test_layer_fluxes( verbose, nk, test_name, F_layer, (/0.0,0.0/) )

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
  call layer_fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
                                    ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  near_boundary_unit_tests = test_layer_fluxes( verbose, nk, test_name, F_layer, (/-7.5,-7.5/) )

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
  call layer_fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
                                    ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  near_boundary_unit_tests = test_layer_fluxes( verbose, nk, test_name, F_layer, (/-7.5,-7.5/) )

  ! Cases where hbl < column thickness (polynomial coefficients specified for pseudo-linear reconstruction)

  test_name = 'hbl < column thickness, hbl same, constant concentration each column'
  hbl_L = 2; hbl_R = 2
  h_L = (/1.,2./) ; h_R = (/1.,2./)
  phi_L = (/0.,0./) ; phi_R = (/1.,1./)
  phi_pp_L(1,1) = 0.; phi_pp_L(1,2) = 0.
  phi_pp_L(2,1) = 0.; phi_pp_L(2,2) = 0.
  phi_pp_R(1,1) = 1.; phi_pp_R(1,2) = 0.
  phi_pp_R(2,1) = 1.; phi_pp_R(2,2) = 0.
  khtr_u = 1.
  call layer_fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
                                    ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  near_boundary_unit_tests = test_layer_fluxes( verbose, nk, test_name, F_layer, (/-1.,-1./) )

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
  call layer_fluxes_bulk_method(SURFACE, nk, deg, h_L, h_R, hbl_L, hbl_R, phi_L, phi_R, phi_pp_L, phi_pp_R,&
                                    ppoly0_E_L, ppoly0_E_R, method, khtr_u, F_bulk, F_layer)
  near_boundary_unit_tests = test_layer_fluxes( verbose, nk, test_name, F_layer, (/-1.,-1./) )
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
  integer, parameter :: stdunit = 6

  test_layer_fluxes = .false.
  do k=1,nk
    if ( F_calc(k) /= F_ans(k) ) then
      test_layer_fluxes = .true.
      write(stdunit,*) "UNIT TEST FAILED: ", test_name
      write(stdunit,10) k, F_calc(k), F_ans(k)
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

  integer, parameter :: stdunit = 6

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
end module MOM_lateral_boundary_mixing
