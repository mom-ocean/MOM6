!> Calculates horizontal viscosity and viscous stresses
module MOM_hor_visc

! This file is part of MOM6. See LICENSE.md for the license.
use MOM_checksums,             only : hchksum, Bchksum, uvchksum
use MOM_coms,                  only : min_across_PEs
use MOM_diag_mediator,         only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator,         only : post_product_u, post_product_sum_u
use MOM_diag_mediator,         only : post_product_v, post_product_sum_v
use MOM_diag_mediator,         only : diag_ctrl, time_type
use MOM_domains,               only : pass_var, CORNER, pass_vector, AGRID, BGRID_NE
use MOM_domains,               only : To_All, Scalar_Pair
use MOM_error_handler,         only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_grid,                  only : ocean_grid_type
use MOM_interface_heights,     only : thickness_to_dz
use MOM_lateral_mixing_coeffs, only : VarMix_CS, calc_QG_slopes, calc_QG_Leith_viscosity
use MOM_barotropic,            only : barotropic_CS, barotropic_get_tav
use MOM_thickness_diffuse,     only : thickness_diffuse_CS, thickness_diffuse_get_KH
use MOM_io,                    only : MOM_read_data, slasher
use MOM_MEKE_types,            only : MEKE_type
use MOM_open_boundary,         only : ocean_OBC_type, OBC_DIRECTION_E, OBC_DIRECTION_W
use MOM_open_boundary,         only : OBC_DIRECTION_N, OBC_DIRECTION_S, OBC_NONE
use MOM_unit_scaling,          only : unit_scale_type
use MOM_verticalGrid,          only : verticalGrid_type
use MOM_variables,             only : accel_diag_ptrs, thermo_var_ptrs
use MOM_Zanna_Bolton,          only : ZB2020_lateral_stress, ZB2020_init, ZB2020_end
use MOM_Zanna_Bolton,          only : ZB2020_CS, ZB2020_copy_gradient_and_thickness

implicit none ; private

#include <MOM_memory.h>

public horizontal_viscosity, hor_visc_init, hor_visc_end

!> Control structure for horizontal viscosity
type, public :: hor_visc_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  logical :: Laplacian       !< Use a Laplacian horizontal viscosity if true.
  logical :: biharmonic      !< Use a biharmonic horizontal viscosity if true.
  logical :: debug           !< If true, write verbose checksums for debugging purposes.
  logical :: no_slip         !< If true, no slip boundary conditions are used.
                             !! Otherwise free slip boundary conditions are assumed.
                             !! The implementation of the free slip boundary
                             !! conditions on a C-grid is much cleaner than the
                             !! no slip boundary conditions. The use of free slip
                             !! b.c.s is strongly encouraged. The no slip b.c.s
                             !! are not implemented with the biharmonic viscosity.
  logical :: bound_Kh        !< If true, the Laplacian coefficient is locally
                             !! limited to guarantee stability.
  logical :: better_bound_Kh !< If true, use a more careful bounding of the
                             !! Laplacian viscosity to guarantee stability.
  logical :: bound_Ah        !< If true, the biharmonic coefficient is locally
                             !! limited to guarantee stability.
  logical :: better_bound_Ah !< If true, use a more careful bounding of the
                             !! biharmonic viscosity to guarantee stability.
  real    :: Re_Ah           !! If nonzero, the biharmonic coefficient is scaled
                             !< so that the biharmonic Reynolds number is equal to this [nondim].
  real    :: bound_coef      !< The nondimensional coefficient of the ratio of
                             !! the viscosity bounds to the theoretical maximum
                             !! for stability without considering other terms [nondim].
                             !! The default is 0.8.
  logical :: backscatter_underbound !< If true, the bounds on the biharmonic viscosity are allowed
                             !! to increase where the Laplacian viscosity is negative (due to
                             !! backscatter parameterizations) beyond the largest timestep-dependent
                             !! stable values of biharmonic viscosity when no Laplacian viscosity is
                             !! applied.  The default is true for historical reasons, but this option
                             !! probably should not be used as it can lead to numerical instabilities.
  logical :: Smagorinsky_Kh  !< If true, use Smagorinsky nonlinear eddy
                             !! viscosity. KH is the background value.
  logical :: Smagorinsky_Ah  !< If true, use a biharmonic form of Smagorinsky
                             !! nonlinear eddy viscosity. AH is the background.
  logical :: Leith_Kh        !< If true, use 2D Leith nonlinear eddy
                             !! viscosity. KH is the background value.
  logical :: Modified_Leith  !< If true, use extra component of Leith viscosity
                             !! to damp divergent flow. To use, still set Leith_Kh=.TRUE.
  logical :: use_beta_in_Leith !< If true, includes the beta term in the Leith viscosity
  logical :: Leith_Ah        !< If true, use a biharmonic form of 2D Leith
                             !! nonlinear eddy viscosity. AH is the background.
  logical :: use_Leithy      !< If true, use a biharmonic form of 2D Leith
                             !! nonlinear eddy viscosity with harmonic backscatter.
                             !! Ah is the background. Leithy = Leith+E
  real    :: c_K             !< Fraction of energy dissipated by the biharmonic term
                             !! that gets backscattered in the Leith+E scheme. [nondim]
  logical :: smooth_Ah       !< If true (default), then Ah and m_leithy are smoothed.
                             !! This smoothing requires a lot of blocking communication.
  logical :: use_QG_Leith_visc    !< If true, use QG Leith nonlinear eddy viscosity.
                             !! KH is the background value.
  logical :: bound_Coriolis  !< If true & SMAGORINSKY_AH is used, the biharmonic
                             !! viscosity is modified to include a term that
                             !! scales quadratically with the velocity shears.
  logical :: use_Kh_bg_2d    !< Read 2d background viscosity from a file.
  logical :: Kh_bg_2d_bug    !< If true, retain an answer-changing horizontal indexing bug
                             !! in setting the corner-point viscosities when USE_KH_BG_2D=True.
  real    :: Kh_bg_min       !< The minimum value allowed for Laplacian horizontal
                             !! viscosity [L2 T-1 ~> m2 s-1]. The default is 0.0.
  logical :: FrictWork_bug    !< If true, retain an answer-changing bug in calculating FrictWork,
                             !! which cancels the h in thickness flux and the h at velocity point.
  logical :: use_land_mask   !< Use the land mask for the computation of thicknesses
                             !! at velocity locations. This eliminates the dependence on
                             !! arbitrary values over land or outside of the domain.
                             !! Default is False to maintain answers with legacy experiments
                             !! but should be changed to True for new experiments.
  logical :: anisotropic     !< If true, allow anisotropic component to the viscosity.
  logical :: add_LES_viscosity!< If true, adds the viscosity from Smagorinsky and Leith to
                             !! the background viscosity instead of taking the maximum.
  real    :: Kh_aniso        !< The anisotropic viscosity [L2 T-1 ~> m2 s-1].
  logical :: dynamic_aniso   !< If true, the anisotropic viscosity is recomputed as a function
                             !! of state. This is set depending on ANISOTROPIC_MODE.
  logical :: res_scale_MEKE  !< If true, the viscosity contribution from MEKE is scaled by
                             !! the resolution function.
  logical :: use_GME         !< If true, use GME backscatter scheme.
  integer :: answer_date     !< The vintage of the order of arithmetic and expressions in the
                             !! horizontal viscosity calculations.  Values below 20190101 recover
                             !! the answers from the end of 2018, while higher values use updated
                             !! and more robust forms of the same expressions.
  real    :: GME_h0          !< The strength of GME tapers quadratically to zero when the bathymetric
                             !! total water column thickness is less than GME_H0 [H ~> m or kg m-2]
  real    :: GME_efficiency  !< The nondimensional prefactor multiplying the GME coefficient [nondim]
  real    :: GME_limiter     !< The absolute maximum value the GME coefficient is allowed to take [L2 T-1 ~> m2 s-1].
  real    :: min_grid_Kh     !< Minimum horizontal Laplacian viscosity used to
                             !! limit the grid Reynolds number [L2 T-1 ~> m2 s-1]
  real    :: min_grid_Ah     !< Minimun horizontal biharmonic viscosity used to
                             !! limit grid Reynolds number [L4 T-1 ~> m4 s-1]
  logical :: use_cont_thick  !< If true, thickness at velocity points adopts h[uv] in BT_cont from continuity solver.
  type(ZB2020_CS) :: ZB2020  !< Zanna-Bolton 2020 control structure.
  logical :: use_ZB2020      !< If true, use Zanna-Bolton 2020 parameterization.

  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: Kh_bg_xx
                      !< The background Laplacian viscosity at h points [L2 T-1 ~> m2 s-1].
                      !! The actual viscosity may be the larger of this
                      !! viscosity and the Smagorinsky and Leith viscosities.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: Kh_bg_2d
                      !< The background Laplacian viscosity at h points [L2 T-1 ~> m2 s-1].
                      !! The actual viscosity may be the larger of this
                      !! viscosity and the Smagorinsky and Leith viscosities.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: Ah_bg_xx
                      !< The background biharmonic viscosity at h points [L4 T-1 ~> m4 s-1].
                      !! The actual viscosity may be the larger of this
                      !! viscosity and the Smagorinsky and Leith viscosities.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: reduction_xx
                      !< The amount by which stresses through h points are reduced
                      !! due to partial barriers [nondim].
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    Kh_Max_xx,      & !< The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1].
    Ah_Max_xx,      & !< The maximum permitted biharmonic viscosity [L4 T-1 ~> m4 s-1].
    n1n2_h,         & !< Factor n1*n2 in the anisotropic direction tensor at h-points [nondim]
    n1n1_m_n2n2_h,  & !< Factor n1**2-n2**2 in the anisotropic direction tensor at h-points [nondim]
    grid_sp_h2,     & !< Harmonic mean of the squares of the grid [L2 ~> m2]
    grid_sp_h3        !< Harmonic mean of the squares of the grid^(3/2) [L3 ~> m3]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: Kh_bg_xy
                      !< The background Laplacian viscosity at q points [L2 T-1 ~> m2 s-1].
                      !! The actual viscosity may be the larger of this
                      !! viscosity and the Smagorinsky and Leith viscosities.
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: Ah_bg_xy
                      !< The background biharmonic viscosity at q points [L4 T-1 ~> m4 s-1].
                      !! The actual viscosity may be the larger of this
                      !! viscosity and the Smagorinsky and Leith viscosities.
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: reduction_xy
                      !< The amount by which stresses through q points are reduced
                      !! due to partial barriers [nondim].
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    Kh_Max_xy,      & !< The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1].
    Ah_Max_xy,      & !< The maximum permitted biharmonic viscosity [L4 T-1 ~> m4 s-1].
    n1n2_q,         & !< Factor n1*n2 in the anisotropic direction tensor at q-points [nondim]
    n1n1_m_n2n2_q     !< Factor n1**2-n2**2 in the anisotropic direction tensor at q-points [nondim]

  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    dx2h,           & !< Pre-calculated dx^2 at h points [L2 ~> m2]
    dy2h,           & !< Pre-calculated dy^2 at h points [L2 ~> m2]
    dx_dyT,         & !< Pre-calculated dx/dy at h points [nondim]
    dy_dxT,         & !< Pre-calculated dy/dx at h points [nondim]
    m_const_leithy, & !< Pre-calculated .5*sqrt(c_K)*max{dx,dy} [L ~> m]
    m_leithy_max      !< Pre-calculated 4./max(dx,dy)^2 at h points [L-2 ~> m-2]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    dx2q,    & !< Pre-calculated dx^2 at q points [L2 ~> m2]
    dy2q,    & !< Pre-calculated dy^2 at q points [L2 ~> m2]
    dx_dyBu, & !< Pre-calculated dx/dy at q points [nondim]
    dy_dxBu    !< Pre-calculated dy/dx at q points [nondim]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: &
    Idx2dyCu, & !< 1/(dx^2 dy) at u points [L-3 ~> m-3]
    Idxdy2u     !< 1/(dx dy^2) at u points [L-3 ~> m-3]
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: &
    Idx2dyCv, & !< 1/(dx^2 dy) at v points [L-3 ~> m-3]
    Idxdy2v     !< 1/(dx dy^2) at v points [L-3 ~> m-3]

  ! The following variables are precalculated time-invariant combinations of
  ! parameters and metric terms.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    Laplac2_const_xx, & !< Laplacian  metric-dependent constants [L2 ~> m2]
    Biharm6_const_xx, & !< Biharmonic metric-dependent constants [L6 ~> m6]
    Laplac3_const_xx, & !< Laplacian  metric-dependent constants [L3 ~> m3]
    Biharm_const_xx,  & !< Biharmonic metric-dependent constants [L4 ~> m4]
    Biharm_const2_xx, & !< Biharmonic metric-dependent constants [T L4 ~> s m4]
    Re_Ah_const_xx      !< Biharmonic metric-dependent constants [L3 ~> m3]

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    Laplac2_const_xy, & !< Laplacian  metric-dependent constants [L2 ~> m2]
    Biharm6_const_xy, & !< Biharmonic metric-dependent constants [L6 ~> m6]
    Laplac3_const_xy, & !< Laplacian  metric-dependent constants [L3 ~> m3]
    Biharm_const_xy,  & !< Biharmonic metric-dependent constants [L4 ~> m4]
    Biharm_const2_xy, & !< Biharmonic metric-dependent constants [T L4 ~> s m4]
    Re_Ah_const_xy      !< Biharmonic metric-dependent constants [L3 ~> m3]

  type(diag_ctrl), pointer :: diag => NULL() !< structure to regulate diagnostics

  ! real, allocatable :: hf_diffu(:,:,:)  ! Zonal horizontal viscous acceleleration times
  !                                       ! fractional thickness [L T-2 ~> m s-2].
  ! real, allocatable :: hf_diffv(:,:,:)  ! Meridional horizontal viscous acceleleration times
  !                                       ! fractional thickness [L T-2 ~> m s-2].
  ! 3D diagnostics hf_diffu(diffv) are commented because there is no clarity on proper remapping grid option.
  ! The code is retained for debugging purposes in the future.

  integer :: num_smooth_gme !< number of smoothing passes for the GME fluxes.
  !>@{
  !! Diagnostic id
  integer :: id_grid_Re_Ah = -1, id_grid_Re_Kh   = -1
  integer :: id_diffu     = -1, id_diffv         = -1
  ! integer :: id_hf_diffu  = -1, id_hf_diffv      = -1
  integer :: id_h_diffu  = -1, id_h_diffv      = -1
  integer :: id_hf_diffu_2d = -1, id_hf_diffv_2d = -1
  integer :: id_intz_diffu_2d = -1, id_intz_diffv_2d = -1
  integer :: id_diffu_visc_rem = -1, id_diffv_visc_rem = -1
  integer :: id_Ah_h      = -1, id_Ah_q          = -1
  integer :: id_Kh_h      = -1, id_Kh_q          = -1
  integer :: id_GME_coeff_h = -1, id_GME_coeff_q = -1
  integer :: id_dudx_bt = -1, id_dvdy_bt = -1
  integer :: id_dudy_bt = -1, id_dvdx_bt = -1
  integer :: id_vort_xy_q = -1, id_div_xx_h      = -1
  integer :: id_sh_xy_q = -1,    id_sh_xx_h      = -1
  integer :: id_FrictWork = -1, id_FrictWorkIntz = -1
  integer :: id_FrictWork_GME = -1
  integer :: id_normstress = -1, id_shearstress = -1
  !>@}

end type hor_visc_CS

contains

!> Calculates the acceleration due to the horizontal viscosity.
!!
!! A combination of biharmonic and Laplacian forms can be used. The coefficient
!! may either be a constant or a shear-dependent form. The biharmonic is
!! determined by twice taking the divergence of an appropriately defined stress
!! tensor. The Laplacian is determined by doing so once.
!!
!! To work, the following fields must be set outside of the usual
!! is:ie range before this subroutine is called:
!!   u(is-2:ie+2,js-2:je+2)
!!   v(is-2:ie+2,js-2:je+2)
!!   h(is-1:ie+1,js-1:je+1) or up to h(is-2:ie+2,js-2:je+2) with some Leith options.
subroutine horizontal_viscosity(u, v, h, uh, vh, diffu, diffv, MEKE, VarMix, G, GV, US, &
                                CS, tv, dt, OBC, BT, TD, ADp, hu_cont, hv_cont)
  type(ocean_grid_type),         intent(in)  :: G      !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)  :: GV     !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in)  :: u      !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in)  :: v      !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                 intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in)  :: uh      !< The zonal volume transport [H L2 T-1 ~> m3 s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in)  :: vh      !< The meridional volume transport [H L2 T-1 ~> m3 s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(out) :: diffu  !< Zonal acceleration due to convergence of
                                                       !! along-coordinate stress tensor [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(out) :: diffv  !< Meridional acceleration due to convergence
                                                       !! of along-coordinate stress tensor [L T-2 ~> m s-2].
  type(MEKE_type),               intent(inout) :: MEKE !< MEKE fields
                                                       !! related to Mesoscale Eddy Kinetic Energy.
  type(VarMix_CS),               intent(inout) :: VarMix !< Variable mixing control structure
  type(unit_scale_type),         intent(in)    :: US   !< A dimensional unit scaling type
  type(hor_visc_CS),             intent(inout) :: CS   !< Horizontal viscosity control structure
  type(thermo_var_ptrs),         intent(in)    :: tv   !< A structure pointing to various
                                                       !! thermodynamic variables
  real,                          intent(in)    :: dt   !< Time increment [T ~> s]
  type(ocean_OBC_type), optional, pointer      :: OBC  !< Pointer to an open boundary condition type
  type(barotropic_CS), optional, intent(in)    :: BT   !< Barotropic control structure
  type(thickness_diffuse_CS), optional, intent(in) :: TD !< Thickness diffusion control structure
  type(accel_diag_ptrs), optional, intent(in)  :: ADp  !< Acceleration diagnostics
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                        optional, intent(in) :: hu_cont !< Layer thickness at u-points [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                        optional, intent(in) :: hv_cont !< Layer thickness at v-points [H ~> m or kg m-2].

  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    Del2u, &      ! The u-component of the Laplacian of velocity [L-1 T-1 ~> m-1 s-1]
    h_u, &        ! Thickness interpolated to u points [H ~> m or kg m-2].
    vort_xy_dy, & ! y-derivative of vertical vorticity (d/dy(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
    vort_xy_dy_smooth, & ! y-derivative of smoothed vertical vorticity [L-1 T-1 ~> m-1 s-1]
    div_xx_dx, &  ! x-derivative of horizontal divergence (d/dx(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
    ubtav         ! zonal barotropic velocity averaged over a baroclinic time-step [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G)) :: &
    Del2v, &      ! The v-component of the Laplacian of velocity [L-1 T-1 ~> m-1 s-1]
    h_v, &        ! Thickness interpolated to v points [H ~> m or kg m-2].
    vort_xy_dx, & ! x-derivative of vertical vorticity (d/dx(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
    vort_xy_dx_smooth, & ! x-derivative of smoothed vertical vorticity [L-1 T-1 ~> m-1 s-1]
    div_xx_dy, &  ! y-derivative of horizontal divergence (d/dy(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
    vbtav         ! meridional barotropic velocity averaged over a baroclinic time-step [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G)) :: &
    dudx_bt, dvdy_bt, & ! components in the barotropic horizontal tension [T-1 ~> s-1]
    div_xx, &     ! Estimate of horizontal divergence at h-points [T-1 ~> s-1]
    sh_xx, &      ! horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
    sh_xx_smooth, & ! horizontal tension from smoothed velocity including metric terms [T-1 ~> s-1]
    sh_xx_bt, &   ! barotropic horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
    str_xx,&      ! str_xx is the diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but
                  ! at some points in the code it is not yet layer integrated, so is in [L2 T-2 ~> m2 s-2].
    str_xx_GME,&  ! smoothed diagonal term in the stress tensor from GME [L2 T-2 ~> m2 s-2]
    bhstr_xx, &   ! A copy of str_xx that only contains the biharmonic contribution [H L2 T-2 ~> m3 s-2 or kg s-2]
    FrictWorkIntz, & ! depth integrated energy dissipated by lateral friction [R L2 T-3 ~> W m-2]
    grad_vort_mag_h, & ! Magnitude of vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
    grad_vort_mag_h_2d, & ! Magnitude of 2d vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
    grad_div_mag_h, &     ! Magnitude of divergence gradient at h-points [L-1 T-1 ~> m-1 s-1]
    dudx, dvdy, &    ! components in the horizontal tension [T-1 ~> s-1]
    dudx_smooth, dvdy_smooth, & ! components in the horizontal tension from smoothed velocity [T-1 ~> s-1]
    GME_effic_h, &  ! The filtered efficiency of the GME terms at h points [nondim]
    m_leithy, &   ! Kh=m_leithy*Ah in Leith+E parameterization [L-2 ~> m-2]
    Ah_sq, &      ! The square of the biharmonic viscosity [L8 T-2 ~> m8 s-2]
    htot          ! The total thickness of all layers [H ~> m or kg m-2]
  real :: Del2vort_h ! Laplacian of vorticity at h-points [L-2 T-1 ~> m-2 s-1]
  real :: grad_vel_mag_bt_h ! Magnitude of the barotropic velocity gradient tensor squared at h-points [T-2 ~> s-2]
  real :: boundary_mask_h ! A mask that zeroes out cells with at least one land edge [nondim]

  real, dimension(SZIB_(G),SZJB_(G)) :: &
    dvdx, dudy, & ! components in the shearing strain [T-1 ~> s-1]
    dvdx_smooth, dudy_smooth, & ! components in the shearing strain from smoothed velocity [T-1 ~> s-1]
    dDel2vdx, dDel2udy, & ! Components in the biharmonic equivalent of the shearing strain [L-2 T-1 ~> m-2 s-1]
    dvdx_bt, dudy_bt,   & ! components in the barotropic shearing strain [T-1 ~> s-1]
    sh_xy,  &     ! horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
    sh_xy_smooth,  & ! horizontal shearing strain from smoothed velocity including metric terms [T-1 ~> s-1]
    sh_xy_bt, &   ! barotropic horizontal shearing strain (du/dy + dv/dx) inc. metric terms [T-1 ~> s-1]
    str_xy, &     ! str_xy is the cross term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but
                  ! at some points in the code it is not yet layer integrated, so is in [L2 T-2 ~> m2 s-2].
    str_xy_GME, & ! smoothed cross term in the stress tensor from GME [L2 T-2 ~> m2 s-2]
    bhstr_xy, &   ! A copy of str_xy that only contains the biharmonic contribution [H L2 T-2 ~> m3 s-2 or kg s-2]
    vort_xy, &    ! Vertical vorticity (dv/dx - du/dy) including metric terms [T-1 ~> s-1]
    vort_xy_smooth, & ! Vertical vorticity including metric terms, smoothed [T-1 ~> s-1]
    grad_vort_mag_q, & ! Magnitude of vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
    grad_vort_mag_q_2d, & ! Magnitude of 2d vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
    Del2vort_q, & ! Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
    grad_div_mag_q, &  ! Magnitude of divergence gradient at q-points [L-1 T-1 ~> m-1 s-1]
    hq, &          ! harmonic mean of the harmonic means of the u- & v point thicknesses [H ~> m or kg m-2]
                   ! This form guarantees that hq/hu < 4.
    GME_effic_q    ! The filtered efficiency of the GME terms at q points [nondim]
  real :: grad_vel_mag_bt_q ! Magnitude of the barotropic velocity gradient tensor squared at q-points [T-2 ~> s-2]
  real :: boundary_mask_q ! A mask that zeroes out cells with at least one land edge [nondim]

  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)) :: &
    Ah_q, &      ! biharmonic viscosity at corner points [L4 T-1 ~> m4 s-1]
    Kh_q, &      ! Laplacian viscosity at corner points [L2 T-1 ~> m2 s-1]
    vort_xy_q, & ! vertical vorticity at corner points [T-1 ~> s-1]
    sh_xy_q,   & ! horizontal shearing strain at corner points [T-1 ~> s-1]
    GME_coeff_q, &  !< GME coeff. at q-points [L2 T-1 ~> m2 s-1]
    ShSt         ! A diagnostic array of shear stress [T-1 ~> s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1) :: &
    KH_u_GME, &  !< Isopycnal height diffusivities in u-columns [L2 T-1 ~> m2 s-1]
    slope_x      !< Isopycnal slope in i-direction [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1) :: &
    KH_v_GME, &  !< Isopycnal height diffusivities in v-columns [L2 T-1 ~> m2 s-1]
    slope_y      !< Isopycnal slope in j-direction [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    Ah_h, &          ! biharmonic viscosity at thickness points [L4 T-1 ~> m4 s-1]
    Kh_h, &          ! Laplacian viscosity at thickness points [L2 T-1 ~> m2 s-1]
    dz, &            ! Height change across layers [Z ~> m]
    FrictWork, &     ! work done by MKE dissipation mechanisms [R L2 T-3 ~> W m-2]
    FrictWork_GME, & ! work done by GME [R L2 T-3 ~> W m-2]
    div_xx_h,      & ! horizontal divergence [T-1 ~> s-1]
    sh_xx_h,       & ! horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
    NoSt             ! A diagnostic array of normal stress [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    grid_Re_Kh, &    ! Grid Reynolds number for Laplacian horizontal viscosity at h points [nondim]
    grid_Re_Ah, &    ! Grid Reynolds number for Biharmonic horizontal viscosity at h points [nondim]
    GME_coeff_h      ! GME coefficient at h-points [L2 T-1 ~> m2 s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: &
    u_smooth         ! Zonal velocity, smoothed with a spatial low-pass filter [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: &
    v_smooth         ! Meridional velocity, smoothed with a spatial low-pass filter [L T-1 ~> m s-1]
  real :: AhSm       ! Smagorinsky biharmonic viscosity [L4 T-1 ~> m4 s-1]
  real :: AhLth      ! 2D Leith biharmonic viscosity [L4 T-1 ~> m4 s-1]
  real :: AhLthy     ! 2D Leith+E biharmonic viscosity [L4 T-1 ~> m4 s-1]
  real :: Shear_mag_bc  ! Shear_mag value in backscatter [T-1 ~> s-1]
  real :: sh_xx_sq   ! Square of tension (sh_xx) [T-2 ~> s-2]
  real :: sh_xy_sq   ! Square of shearing strain (sh_xy) [T-2 ~> s-2]
  real :: h2uq, h2vq ! temporary variables [H2 ~> m2 or kg2 m-4].
  real :: hu, hv     ! Thicknesses interpolated by arithmetic means to corner
                     ! points; these are first interpolated to u or v velocity
                     ! points where masks are applied [H ~> m or kg m-2].
  real :: h_arith_q  ! The arithmetic mean total thickness at q points [H ~> m or kg m-2]
  real :: I_GME_h0   ! The inverse of GME tapering scale [H-1 ~> m-1 or m2 kg-1]
  real :: h_neglect  ! thickness so small it can be lost in roundoff and so neglected [H ~> m or kg m-2]
  real :: h_neglect3 ! h_neglect^3 [H3 ~> m3 or kg3 m-6]
  real :: h_min      ! Minimum h at the 4 neighboring velocity points [H ~> m]
  real :: Kh_max_here ! The local maximum Laplacian viscosity for stability [L2 T-1 ~> m2 s-1]
  real :: RoScl     ! The scaling function for MEKE source term [nondim]
  real :: FatH      ! abs(f) at h-point for MEKE source term [T-1 ~> s-1]
  real :: local_strain ! Local variable for interpolating computed strain rates [T-1 ~> s-1].
  real :: meke_res_fn ! A copy of the resolution scaling factor if being applied to MEKE [nondim]. Otherwise = 1.
  real :: GME_coeff ! The GME (negative) viscosity coefficient [L2 T-1 ~> m2 s-1]
  real :: DY_dxBu   ! Ratio of meridional over zonal grid spacing at vertices [nondim]
  real :: DX_dyBu   ! Ratio of zonal over meridional grid spacing at vertices [nondim]
  real :: Sh_F_pow  ! The ratio of shear over the absolute value of f raised to some power and rescaled [nondim]
  real :: backscat_subround ! The ratio of f over Shear_mag that is so small that the backscatter
                    ! calculation gives the same value as if f were 0 [nondim].
  real :: KE        ! Local kinetic energy [L2 T-2 ~> m2 s-2]
  real :: d_del2u   ! dy-weighted Laplacian(u) diff in x [L-2 T-1 ~> m-2 s-1]
  real :: d_del2v   ! dx-weighted Laplacian(v) diff in y [L-2 T-1 ~> m-2 s-1]
  real :: d_str     ! Stress tensor update [L2 T-2 ~> m2 s-2]
  real :: grad_vort ! Vorticity gradient magnitude [L-1 T-1 ~> m-1 s-1]
  real :: grad_vort_qg ! QG-based vorticity gradient magnitude [L-1 T-1 ~> m-1 s-1]
  real :: grid_Kh   ! Laplacian viscosity bound by grid [L2 T-1 ~> m2 s-1]
  real :: grid_Ah   ! Biharmonic viscosity bound by grid [L4 T-1 ~> m4 s-1]

  logical :: rescale_Kh, legacy_bound
  logical :: find_FrictWork
  logical :: apply_OBC = .false.
  logical :: use_MEKE_Ku
  logical :: use_MEKE_Au
  logical :: use_cont_huv
  integer :: is_vort, ie_vort, js_vort, je_vort  ! Loop ranges for vorticity terms
  integer :: is_Kh, ie_Kh, js_Kh, je_Kh  ! Loop ranges for thickness point viscosities
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k, n
  real :: inv_PI3, inv_PI2, inv_PI6 ! Powers of the inverse of pi [nondim]

  ! Fields evaluated on active layers, used for constructing 3D stress fields
  ! NOTE: The position of these declarations can impact performance, due to the
  !   very large number of stack arrays in this function.  Move with caution!
  ! NOTE: Several of these are declared with the memory extent of q-points, but the
  !   same arrays are also used at h-points to reduce the memory footprint of this
  !   module, so they should never be used in halo point or checksum calls.
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    Ah, &           ! biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
    Kh, &           ! Laplacian  viscosity (h or q) [L2 T-1 ~> m2 s-1]
    Shear_mag, &    ! magnitude of the shear (h or q) [T-1 ~> s-1]
    vert_vort_mag, &  ! magnitude of the vertical vorticity gradient (h or q) [L-1 T-1 ~> m-1 s-1]
    vert_vort_mag_smooth, &  ! magnitude of gradient of smoothed vertical vorticity (h or q) [L-1 T-1 ~> m-1 s-1]
    hrat_min, &     ! h_min divided by the thickness at the stress point (h or q) [nondim]
    visc_bound_rem  ! fraction of overall viscous bounds that remain to be applied (h or q) [nondim]

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  h_neglect  = GV%H_subroundoff
  !h_neglect3 = h_neglect**3
  h_neglect3 = h_neglect*h_neglect*h_neglect
  inv_PI3 = 1.0/((4.0*atan(1.0))**3)
  inv_PI2 = 1.0/((4.0*atan(1.0))**2)
  inv_PI6 = inv_PI3 * inv_PI3

  m_leithy(:,:) = 0.0 ! Initialize

  if (present(OBC)) then ; if (associated(OBC)) then ; if (OBC%OBC_pe) then
    apply_OBC = OBC%Flather_u_BCs_exist_globally .or. OBC%Flather_v_BCs_exist_globally
    apply_OBC = .true.
  endif ; endif ; endif

  if (.not.CS%initialized) call MOM_error(FATAL, &
         "MOM_hor_visc: Module must be initialized before it is used.")

  if (.not.(CS%Laplacian .or. CS%biharmonic)) return

  find_FrictWork = (CS%id_FrictWork > 0)
  if (CS%id_FrictWorkIntz > 0) find_FrictWork = .true.

  if (allocated(MEKE%mom_src)) find_FrictWork = .true.
  backscat_subround = 0.0
  if (find_FrictWork .and. allocated(MEKE%mom_src) .and. (MEKE%backscatter_Ro_c > 0.0) .and. &
      (MEKE%backscatter_Ro_Pow /= 0.0)) &
    backscat_subround = (1.0e-16/MEKE%backscatter_Ro_c)**(1.0/MEKE%backscatter_Ro_Pow)

  ! Toggle whether to use a Laplacian viscosity derived from MEKE
  use_MEKE_Ku = allocated(MEKE%Ku)
  use_MEKE_Au = allocated(MEKE%Au)

  use_cont_huv = CS%use_cont_thick .and. present(hu_cont) .and. present(hv_cont)

  rescale_Kh = .false.
  if (VarMix%use_variable_mixing) then
    rescale_Kh = VarMix%Resoln_scaled_Kh
    if ((rescale_Kh .or. CS%res_scale_MEKE) &
        .and. (.not. allocated(VarMix%Res_fn_h) .or. .not. allocated(VarMix%Res_fn_q))) &
      call MOM_error(FATAL, "MOM_hor_visc: VarMix%Res_fn_h and VarMix%Res_fn_q "//&
                     "both need to be associated with Resoln_scaled_Kh or RES_SCALE_MEKE_VISC.")
  elseif (CS%res_scale_MEKE) then
    call MOM_error(FATAL, "MOM_hor_visc: VarMix needs to be associated if "//&
                          "RES_SCALE_MEKE_VISC is True.")
  endif

  ! Set the halo sizes used for the thickness-point viscosities.
  if (CS%use_Leithy) then
    js_Kh = js-1 ; je_Kh = je+1 ; is_Kh = is-1 ; ie_Kh = ie+1
  else
    js_Kh = Jsq ; je_Kh = je+1 ; is_Kh = Isq ; ie_Kh = ie+1
  endif

  ! Set the halo sizes used for the vorticity calculations.
  if ((CS%Leith_Kh) .or. (CS%Leith_Ah) .or. (CS%use_Leithy)) then
    js_vort = js_Kh-2 ; je_vort = Jeq+2 ; is_vort = is_Kh-2 ; ie_vort = Ieq+2
    if ((G%isc-G%isd < 3) .or. (G%isc-G%isd < 3)) call MOM_error(FATAL, &
          "The minimum halo size is 3 when a Leith viscosity is being used.")
  else
    js_vort = js-2 ; je_vort = Jeq+1 ; is_vort = is-2 ; ie_vort = Ieq+1
  endif

  legacy_bound = (CS%Smagorinsky_Kh .or. CS%Leith_Kh) .and. &
                 (CS%bound_Kh .and. .not.CS%better_bound_Kh)

  if (CS%use_GME) then

    ! Initialize diagnostic arrays with zeros
    GME_coeff_h(:,:,:) = 0.0
    GME_coeff_q(:,:,:) = 0.0
    str_xx_GME(:,:) = 0.0
    str_xy_GME(:,:) = 0.0

    ! Get barotropic velocities and their gradients
    call barotropic_get_tav(BT, ubtav, vbtav, G, US)

    call pass_vector(ubtav, vbtav, G%Domain)
    call pass_var(h, G%domain, halo=2)

    ! Calculate the barotropic horizontal tension
    do j=js-2,je+2 ; do i=is-2,ie+2
      dudx_bt(i,j) = CS%DY_dxT(i,j)*(G%IdyCu(I,j) * ubtav(I,j) - &
                                     G%IdyCu(I-1,j) * ubtav(I-1,j))
      dvdy_bt(i,j) = CS%DX_dyT(i,j)*(G%IdxCv(i,J) * vbtav(i,J) - &
                                     G%IdxCv(i,J-1) * vbtav(i,J-1))
    enddo ; enddo
    do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
      sh_xx_bt(i,j) = dudx_bt(i,j) - dvdy_bt(i,j)
    enddo ; enddo

    ! Components for the barotropic shearing strain
    do J=Jsq-2,Jeq+2 ; do I=Isq-2,Ieq+2
      dvdx_bt(I,J) = CS%DY_dxBu(I,J)*(vbtav(i+1,J)*G%IdyCv(i+1,J) &
                                    - vbtav(i,J)*G%IdyCv(i,J))
      dudy_bt(I,J) = CS%DX_dyBu(I,J)*(ubtav(I,j+1)*G%IdxCu(I,j+1) &
                                    - ubtav(I,j)*G%IdxCu(I,j))
    enddo ; enddo

    if (CS%no_slip) then
      do J=js-2,je+1 ; do I=is-2,ie+1
        sh_xy_bt(I,J) = (2.0-G%mask2dBu(I,J)) * ( dvdx_bt(I,J) + dudy_bt(I,J) )
      enddo ; enddo
    else
      do J=js-2,je+1 ; do I=is-2,ie+1
        sh_xy_bt(I,J) = G%mask2dBu(I,J) * ( dvdx_bt(I,J) + dudy_bt(I,J) )
      enddo ; enddo
    endif

    do j=js-2,je+2 ; do i=is-2,ie+2
      htot(i,j) = 0.0
    enddo ; enddo
    do k=1,nz ; do j=js-2,je+2 ; do i=is-2,ie+2
      htot(i,j) = htot(i,j) + h(i,j,k)
    enddo ; enddo ; enddo

    I_GME_h0 = 1.0 / CS%GME_h0
    do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
      boundary_mask_h = (G%mask2dCu(I,j) * G%mask2dCu(I-1,j)) * (G%mask2dCv(i,J) * G%mask2dCv(i,J-1))
      grad_vel_mag_bt_h = G%mask2dT(I,J) * boundary_mask_h * (dudx_bt(i,j)**2 + dvdy_bt(i,j)**2 + &
            (0.25*((dvdx_bt(I,J)+dvdx_bt(I-1,J-1)) + (dvdx_bt(I,J-1)+dvdx_bt(I-1,J))))**2 + &
            (0.25*((dudy_bt(I,J)+dudy_bt(I-1,J-1)) + (dudy_bt(I,J-1)+dudy_bt(I-1,J))))**2)
      ! Probably the following test could be simplified to
      ! if (boundary_mask_h * G%mask2dT(I,J) > 0.0) then
      if (grad_vel_mag_bt_h > 0.0) then
        GME_effic_h(i,j) = CS%GME_efficiency * G%mask2dT(I,J) * (MIN(htot(i,j) * I_GME_h0, 1.0)**2)
      else
        GME_effic_h(i,j) = 0.0
      endif
    enddo ; enddo

    do J=js-2,je+1 ; do I=is-2,ie+1
      boundary_mask_q = (G%mask2dCv(i,J) * G%mask2dCv(i+1,J)) * (G%mask2dCu(I,j) * G%mask2dCu(I,j+1))
      grad_vel_mag_bt_q = G%mask2dBu(I,J) * boundary_mask_q * (dvdx_bt(I,J)**2 + dudy_bt(I,J)**2 + &
            (0.25*((dudx_bt(i,j)+dudx_bt(i+1,j+1)) + (dudx_bt(i,j+1)+dudx_bt(i+1,j))))**2 + &
            (0.25*((dvdy_bt(i,j)+dvdy_bt(i+1,j+1)) + (dvdy_bt(i,j+1)+dvdy_bt(i+1,j))))**2)
      ! Probably the following test could be simplified to
      ! if (boundary_mask_q * G%mask2dBu(I,J) > 0.0) then
      if (grad_vel_mag_bt_q > 0.0) then
        h_arith_q = 0.25 * ((htot(i,j) + htot(i+1,j+1)) + (htot(i+1,j) + htot(i,j+1)))
        GME_effic_q(I,J) = CS%GME_efficiency * G%mask2dBu(I,J) * (MIN(h_arith_q * I_GME_h0, 1.0)**2)
      else
        GME_effic_q(I,J) = 0.0
      endif
    enddo ; enddo

    call thickness_diffuse_get_KH(TD, KH_u_GME, KH_v_GME, G, GV)

    call pass_vector(KH_u_GME, KH_v_GME, G%domain, To_All+Scalar_Pair)

    if (CS%debug) &
      call uvchksum("GME KH[u,v]_GME", KH_u_GME, KH_v_GME, G%HI, haloshift=2, unscale=US%L_to_m**2*US%s_to_T)

  endif ! use_GME

  if (CS%use_Leithy) then
    ! Smooth the velocity. Right now it happens twice. In the future
    ! one might make the number of smoothing cycles a user-specified parameter
    do k=1,nz
      ! One call applies the filter twice
      u_smooth(:,:,k) = u(:,:,k)
      v_smooth(:,:,k) = v(:,:,k)
      call smooth_x9_uv(G, u_smooth(:,:,k), v_smooth(:,:,k), zero_land=.false.)
    enddo
    call pass_vector(u_smooth, v_smooth, G%Domain)
  endif

  if (CS%use_QG_Leith_visc .and. ((CS%Leith_Kh) .or. (CS%Leith_Ah))) then
    call thickness_to_dz(h, tv, dz, G, GV, US, halo_size=2)
    ! Calculate isopycnal slopes that will be used for some forms of viscosity.
    call calc_QG_slopes(h, tv, dt, G, GV, US, slope_x, slope_y, VarMix, OBC)
    ! If the following halo update is added, the calculations in calc_QG_slopes could work on just
    ! the computational domains, and some halo updates outside of this routine could be smaller.
    ! call pass_vector(slope_x, slope_y, G%Domain, halo=2)
  endif

  !$OMP parallel do default(none) &
  !$OMP shared( &
  !$OMP   CS, G, GV, US, OBC, VarMix, MEKE, u, v, h, uh, vh, &
  !$OMP   is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, &
  !$OMP   is_vort, ie_vort, js_vort, je_vort, &
  !$OMP   is_Kh, ie_Kh, js_Kh, je_Kh, &
  !$OMP   apply_OBC, rescale_Kh, legacy_bound, find_FrictWork, &
  !$OMP   use_MEKE_Ku, use_MEKE_Au, u_smooth, v_smooth, use_cont_huv, slope_x, slope_y, dz, &
  !$OMP   backscat_subround, GME_effic_h, GME_effic_q, &
  !$OMP   h_neglect, h_neglect3, inv_PI3, inv_PI6, &
  !$OMP   diffu, diffv, Kh_h, Kh_q, Ah_h, Ah_q, FrictWork, FrictWork_GME, &
  !$OMP   div_xx_h, sh_xx_h, vort_xy_q, sh_xy_q, GME_coeff_h, GME_coeff_q, &
  !$OMP   KH_u_GME, KH_v_GME, grid_Re_Kh, grid_Re_Ah, NoSt, ShSt, hu_cont, hv_cont &
  !$OMP ) &
  !$OMP private( &
  !$OMP   i, j, k, n, &
  !$OMP   dudx, dudy, dvdx, dvdy, sh_xx, sh_xy, h_u, h_v, &
  !$OMP   Del2u, Del2v, DY_dxBu, DX_dyBu, sh_xx_bt, sh_xy_bt, &
  !$OMP   str_xx, str_xy, bhstr_xx, bhstr_xy, str_xx_GME, str_xy_GME, &
  !$OMP   vort_xy, vort_xy_dx, vort_xy_dy, div_xx, div_xx_dx, div_xx_dy, &
  !$OMP   grad_div_mag_h, grad_div_mag_q, grad_vort_mag_h, grad_vort_mag_q, &
  !$OMP   grad_vort, grad_vort_qg, grad_vort_mag_h_2d, grad_vort_mag_q_2d, &
  !$OMP   sh_xx_sq, sh_xy_sq, meke_res_fn, Shear_mag, Shear_mag_bc, vert_vort_mag, &
  !$OMP   h_min, hrat_min, visc_bound_rem, Kh_max_here, &
  !$OMP   grid_Ah, grid_Kh, d_Del2u, d_Del2v, d_str, &
  !$OMP   Kh, Ah, AhSm, AhLth, local_strain, Sh_F_pow, &
  !$OMP   dDel2vdx, dDel2udy, Del2vort_q, Del2vort_h, KE, &
  !$OMP   h2uq, h2vq, hu, hv, hq, FatH, RoScl, GME_coeff, &
  !$OMP   dudx_smooth, dudy_smooth, dvdx_smooth, dvdy_smooth, &
  !$OMP   vort_xy_smooth, vort_xy_dx_smooth, vort_xy_dy_smooth, &
  !$OMP   sh_xx_smooth, sh_xy_smooth, &
  !$OMP   vert_vort_mag_smooth, m_leithy, Ah_sq, AhLthy &
  !$OMP )
  do k=1,nz

    ! The following are the forms of the horizontal tension and horizontal
    ! shearing strain advocated by Smagorinsky (1993) and discussed in
    ! Griffies and Hallberg (2000).

    ! NOTE: There is a ~1% speedup when the tension and shearing loops below
    !   are fused (presumably due to shared access of Id[xy]C[uv]).  However,
    !   this breaks the center/vertex index case convention, and also evaluates
    !   the dudx and dvdy terms beyond their valid bounds.
    ! TODO: Explore methods for retaining both the syntax and speedup.

    ! Calculate horizontal tension
    do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
      dudx(i,j) = CS%DY_dxT(i,j)*(G%IdyCu(I,j) * u(I,j,k) - &
                                  G%IdyCu(I-1,j) * u(I-1,j,k))
      dvdy(i,j) = CS%DX_dyT(i,j)*(G%IdxCv(i,J) * v(i,J,k) - &
                                  G%IdxCv(i,J-1) * v(i,J-1,k))
      sh_xx(i,j) = dudx(i,j) - dvdy(i,j)
    enddo ; enddo

    ! Components for the shearing strain
    do J=js_vort,je_vort ; do I=is_vort,ie_vort
      dvdx(I,J) = CS%DY_dxBu(I,J)*(v(i+1,J,k)*G%IdyCv(i+1,J) - v(i,J,k)*G%IdyCv(i,J))
      dudy(I,J) = CS%DX_dyBu(I,J)*(u(I,j+1,k)*G%IdxCu(I,j+1) - u(I,j,k)*G%IdxCu(I,j))
    enddo ; enddo

    if (CS%use_Leithy) then
      ! Calculate horizontal tension from smoothed velocity
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        dudx_smooth(i,j) = CS%DY_dxT(i,j)*(G%IdyCu(I,j) * u_smooth(I,j,k) - &
                                           G%IdyCu(I-1,j) * u_smooth(I-1,j,k))
        dvdy_smooth(i,j) = CS%DX_dyT(i,j)*(G%IdxCv(i,J) * v_smooth(i,J,k) - &
                                           G%IdxCv(i,J-1) * v_smooth(i,J-1,k))
        sh_xx_smooth(i,j) = dudx_smooth(i,j) - dvdy_smooth(i,j)
      enddo ; enddo

      ! Components for the shearing strain from smoothed velocity
      do J=js_Kh-1,je_Kh ; do I=is_Kh-1,ie_Kh
        dvdx_smooth(I,J) = CS%DY_dxBu(I,J) * &
                         (v_smooth(i+1,J,k)*G%IdyCv(i+1,J) - v_smooth(i,J,k)*G%IdyCv(i,J))
        dudy_smooth(I,J) = CS%DX_dyBu(I,J) * &
                         (u_smooth(I,j+1,k)*G%IdxCu(I,j+1) - u_smooth(I,j,k)*G%IdxCu(I,j))
      enddo ; enddo
    endif ! use Leith+E

    if (CS%id_normstress > 0) then
      do j=js,je ; do i=is,ie
        NoSt(i,j,k) = sh_xx(i,j)
      enddo ; enddo
    endif

    ! Interpolate the thicknesses to velocity points.
    ! The extra wide halos are to accommodate the cross-corner-point projections
    ! in OBCs, which are not ordinarily be necessary, and might not be necessary
    ! even with OBCs if the accelerations are zeroed at OBC points, in which
    ! case the j-loop for h_u could collapse to j=js=1,je+1. -RWH
    if (CS%use_land_mask) then
      do j=js-2,je+2 ; do I=is-2,Ieq+1
        h_u(I,j) = 0.5 * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i+1,j)*h(i+1,j,k))
      enddo ; enddo
      do J=js-2,Jeq+1 ; do i=is-2,ie+2
        h_v(i,J) = 0.5 * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i,j+1)*h(i,j+1,k))
      enddo ; enddo
    else
      do j=js-2,je+2 ; do I=is-2,Ieq+1
        h_u(I,j) = 0.5 * (h(i,j,k) + h(i+1,j,k))
      enddo ; enddo
      do J=js-2,Jeq+1 ; do i=is-2,ie+2
        h_v(i,J) = 0.5 * (h(i,j,k) + h(i,j+1,k))
      enddo ; enddo
    endif

    ! The following should obviously be combined with the previous block if adopted.
    if (use_cont_huv) then
      do j=js-2,je+2 ; do I=Isq-1,Ieq+1
        h_u(I,j) = hu_cont(I,j,k)
      enddo ; enddo
      do J=Jsq-1,Jeq+1 ; do i=is-2,ie+2
        h_v(i,J) = hv_cont(i,J,k)
      enddo ; enddo
    endif

    ! Adjust contributions to shearing strain and interpolated values of
    ! thicknesses on open boundaries.
    if (apply_OBC) then ; do n=1,OBC%number_of_segments
      J = OBC%segment(n)%HI%JsdB ; I = OBC%segment(n)%HI%IsdB
      if (OBC%zero_strain .or. OBC%freeslip_strain .or. OBC%computed_strain) then
        if (OBC%segment(n)%is_N_or_S .and. (J >= Js_vort) .and. (J <= Je_vort)) then
          do I = max(OBC%segment(n)%HI%IsdB,Is_vort), min(OBC%segment(n)%HI%IedB,Ie_vort)
            if (OBC%zero_strain) then
              dvdx(I,J) = 0. ; dudy(I,J) = 0.
            elseif (OBC%freeslip_strain) then
              dudy(I,J) = 0.
            elseif (OBC%computed_strain) then
              if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
                dudy(I,J) = 2.0*CS%DX_dyBu(I,J)* &
                            (OBC%segment(n)%tangential_vel(I,J,k) - u(I,j,k))*G%IdxCu(I,j)
              else
                dudy(I,J) = 2.0*CS%DX_dyBu(I,J)* &
                            (u(I,j+1,k) - OBC%segment(n)%tangential_vel(I,J,k))*G%IdxCu(I,j+1)
              endif
            elseif (OBC%specified_strain) then
              if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
                dudy(I,J) = CS%DX_dyBu(I,J)*OBC%segment(n)%tangential_grad(I,J,k)*G%IdxCu(I,j)*G%dxBu(I,J)
              else
                dudy(I,J) = CS%DX_dyBu(I,J)*OBC%segment(n)%tangential_grad(I,J,k)*G%IdxCu(I,j+1)*G%dxBu(I,J)
              endif
            endif
            if (CS%use_Leithy) then
              dvdx_smooth(I,J) = dvdx(I,J)
              dudy_smooth(I,J) = dudy(I,J)
            endif
          enddo
        elseif (OBC%segment(n)%is_E_or_W .and. (I >= is_vort) .and. (I <= ie_vort)) then
          do J = max(OBC%segment(n)%HI%JsdB,js_vort), min(OBC%segment(n)%HI%JedB,je_vort)
            if (OBC%zero_strain) then
              dvdx(I,J) = 0. ; dudy(I,J) = 0.
            elseif (OBC%freeslip_strain) then
              dvdx(I,J) = 0.
            elseif (OBC%computed_strain) then
              if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
                dvdx(I,J) = 2.0*CS%DY_dxBu(I,J)* &
                            (OBC%segment(n)%tangential_vel(I,J,k) - v(i,J,k))*G%IdyCv(i,J)
              else
                dvdx(I,J) = 2.0*CS%DY_dxBu(I,J)* &
                            (v(i+1,J,k) - OBC%segment(n)%tangential_vel(I,J,k))*G%IdyCv(i+1,J)
              endif
            elseif (OBC%specified_strain) then
              if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
                dvdx(I,J) = CS%DY_dxBu(I,J)*OBC%segment(n)%tangential_grad(I,J,k)*G%IdyCv(i,J)*G%dxBu(I,J)
              else
                dvdx(I,J) = CS%DY_dxBu(I,J)*OBC%segment(n)%tangential_grad(I,J,k)*G%IdyCv(i+1,J)*G%dxBu(I,J)
              endif
            endif
            if (CS%use_Leithy) then
              dvdx_smooth(I,J) = dvdx(I,J)
              dudy_smooth(I,J) = dudy(I,J)
            endif
          enddo
        endif
      endif

      if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
        ! There are extra wide halos here to accommodate the cross-corner-point
        ! OBC projections, but they might not be necessary if the accelerations
        ! are always zeroed out at OBC points, in which case the i-loop below
        ! becomes do i=is-1,ie+1. -RWH
        if ((J >= js-2) .and. (J <= Jeq+1)) then
          do i = max(is-2,OBC%segment(n)%HI%isd), min(ie+2,OBC%segment(n)%HI%ied)
            h_v(i,J) = h(i,j,k)
          enddo
        endif
      elseif (OBC%segment(n)%direction == OBC_DIRECTION_S) then
        if ((J >= js-2) .and. (J <= Jeq+1)) then
          do i = max(is-2,OBC%segment(n)%HI%isd), min(ie+2,OBC%segment(n)%HI%ied)
            h_v(i,J) = h(i,j+1,k)
          enddo
        endif
      elseif (OBC%segment(n)%direction == OBC_DIRECTION_E) then
        if ((I >= is-2) .and. (I <= Ieq+1)) then
          do j = max(js-2,OBC%segment(n)%HI%jsd), min(je+2,OBC%segment(n)%HI%jed)
            h_u(I,j) = h(i,j,k)
          enddo
        endif
      elseif (OBC%segment(n)%direction == OBC_DIRECTION_W) then
        if ((I >= is-2) .and. (I <= Ieq+1)) then
          do j = max(js-2,OBC%segment(n)%HI%jsd), min(je+2,OBC%segment(n)%HI%jed)
            h_u(I,j) = h(i+1,j,k)
          enddo
        endif
      endif
    enddo ; endif
    ! Now project thicknesses across corner points on OBCs.
    if (apply_OBC) then ; do n=1,OBC%number_of_segments
      J = OBC%segment(n)%HI%JsdB ; I = OBC%segment(n)%HI%IsdB
      if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
        if ((J >= js-2) .and. (J <= je)) then
          do I = max(is-2,OBC%segment(n)%HI%IsdB), min(Ieq+1,OBC%segment(n)%HI%IedB)
            h_u(I,j+1) = h_u(I,j)
          enddo
        endif
      elseif (OBC%segment(n)%direction == OBC_DIRECTION_S) then
        if ((J >= js-1) .and. (J <= je+1)) then
          do I = max(is-2,OBC%segment(n)%HI%isd), min(Ieq+1,OBC%segment(n)%HI%ied)
            h_u(I,j) = h_u(I,j+1)
          enddo
        endif
      elseif (OBC%segment(n)%direction == OBC_DIRECTION_E) then
        if ((I >= is-2) .and. (I <= ie)) then
          do J = max(js-2,OBC%segment(n)%HI%jsd), min(Jeq+1,OBC%segment(n)%HI%jed)
            h_v(i+1,J) = h_v(i,J)
          enddo
        endif
      elseif (OBC%segment(n)%direction == OBC_DIRECTION_W) then
        if ((I >= is-1) .and. (I <= ie+1)) then
          do J = max(js-2,OBC%segment(n)%HI%jsd), min(Jeq+1,OBC%segment(n)%HI%jed)
            h_v(i,J) = h_v(i+1,J)
          enddo
        endif
      endif
    enddo ; endif

    ! Shearing strain (including no-slip boundary conditions at the 2-D land-sea mask).
    ! dudy and dvdx include modifications at OBCs from above.
    if (CS%no_slip) then
      do J=js-2,Jeq+1 ; do I=is-2,Ieq+1
        sh_xy(I,J) = (2.0-G%mask2dBu(I,J)) * ( dvdx(I,J) + dudy(I,J) )
        if (CS%id_shearstress > 0) ShSt(I,J,k) = sh_xy(I,J)
      enddo ; enddo
    else
      do J=js-2,Jeq+1 ; do I=is-2,Ieq+1
        sh_xy(I,J) = G%mask2dBu(I,J) * ( dvdx(I,J) + dudy(I,J) )
        if (CS%id_shearstress > 0) ShSt(I,J,k) = sh_xy(I,J)
      enddo ; enddo
    endif

    if (CS%use_Leithy) then
      ! Shearing strain (including no-slip boundary conditions at the 2-D land-sea mask).
      ! dudy_smooth and dvdx_smooth do not (yet) include modifications at OBCs from above.
      if (CS%no_slip) then
        do J=js-1,Jeq ; do I=is-1,Ieq
          sh_xy_smooth(I,J) = (2.0-G%mask2dBu(I,J)) * ( dvdx_smooth(I,J) + dudy_smooth(I,J) )
        enddo ; enddo
      else
        do J=js-1,Jeq ; do I=is-1,Ieq
          sh_xy_smooth(I,J) = G%mask2dBu(I,J) * ( dvdx_smooth(I,J) + dudy_smooth(I,J) )
        enddo ; enddo
      endif
    endif ! use Leith+E

    !  Evaluate Del2u = x.Div(Grad u) and Del2v = y.Div( Grad u)
    if (CS%biharmonic) then
      do j=js-1,Jeq+1 ; do I=Isq-1,Ieq+1
        Del2u(I,j) = CS%Idxdy2u(I,j)*(CS%dy2h(i+1,j)*sh_xx(i+1,j) - CS%dy2h(i,j)*sh_xx(i,j)) + &
                     CS%Idx2dyCu(I,j)*(CS%dx2q(I,J)*sh_xy(I,J) - CS%dx2q(I,J-1)*sh_xy(I,J-1))
      enddo ; enddo
      do J=Jsq-1,Jeq+1 ; do i=is-1,Ieq+1
        Del2v(i,J) = CS%Idxdy2v(i,J)*(CS%dy2q(I,J)*sh_xy(I,J) - CS%dy2q(I-1,J)*sh_xy(I-1,J)) - &
                     CS%Idx2dyCv(i,J)*(CS%dx2h(i,j+1)*sh_xx(i,j+1) - CS%dx2h(i,j)*sh_xx(i,j))
      enddo ; enddo
      if (apply_OBC) then ; if (OBC%zero_biharmonic) then
        do n=1,OBC%number_of_segments
          I = OBC%segment(n)%HI%IsdB ; J = OBC%segment(n)%HI%JsdB
          if (OBC%segment(n)%is_N_or_S .and. (J >= Jsq-1) .and. (J <= Jeq+1)) then
            do I=OBC%segment(n)%HI%isd,OBC%segment(n)%HI%ied
              Del2v(i,J) = 0.
            enddo
          elseif (OBC%segment(n)%is_E_or_W .and. (I >= Isq-1) .and. (I <= Ieq+1)) then
            do j=OBC%segment(n)%HI%jsd,OBC%segment(n)%HI%jed
              Del2u(I,j) = 0.
            enddo
          endif
        enddo
      endif ; endif
    endif

    ! Vorticity
    if ((CS%Leith_Kh) .or. (CS%Leith_Ah) .or. (CS%use_Leithy) .or. (CS%id_vort_xy_q>0)) then
      if (CS%no_slip) then
        do J=js_vort,je_vort ; do I=is_vort,ie_vort
          vort_xy(I,J) = (2.0-G%mask2dBu(I,J)) * ( dvdx(I,J) - dudy(I,J) )
        enddo ; enddo
      else
        do J=js_vort,je_vort ; do I=is_vort,ie_vort
          vort_xy(I,J) = G%mask2dBu(I,J) * ( dvdx(I,J) - dudy(I,J) )
        enddo ; enddo
      endif
    endif

    if (CS%use_Leithy) then
      if (CS%no_slip) then
        do J=js_Kh-1,je_Kh ; do I=is_Kh-1,ie_Kh
          vort_xy_smooth(I,J) = (2.0-G%mask2dBu(I,J)) * ( dvdx_smooth(I,J) - dudy_smooth(I,J) )
        enddo ; enddo
      else
        do J=js_Kh-1,je_Kh ; do I=is_Kh-1,ie_Kh
          vort_xy_smooth(I,J) = G%mask2dBu(I,J) * ( dvdx_smooth(I,J) - dudy_smooth(I,J) )
        enddo ; enddo
      endif
    endif


    if ((CS%Leith_Kh) .or. (CS%Leith_Ah) .or. (CS%use_Leithy)) then

      ! Vorticity gradient
      do J=js-2,je_Kh ; do i=is_Kh-1,ie_Kh+1
        DY_dxBu = G%dyBu(I,J) * G%IdxBu(I,J)
        vort_xy_dx(i,J) = DY_dxBu * (vort_xy(I,J) * G%IdyCu(I,j) - vort_xy(I-1,J) * G%IdyCu(I-1,j))
      enddo ; enddo

      do j=js_Kh-1,je_Kh+1 ; do I=is-2,ie_Kh
        DX_dyBu = G%dxBu(I,J) * G%IdyBu(I,J)
        vort_xy_dy(I,j) = DX_dyBu * (vort_xy(I,J) * G%IdxCv(i,J) - vort_xy(I,J-1) * G%IdxCv(i,J-1))
      enddo ; enddo

      if (CS%use_Leithy) then
        ! Gradient of smoothed vorticity
        do J=js_Kh-1,je_Kh ; do i=is_Kh,ie_Kh
          DY_dxBu = G%dyBu(I,J) * G%IdxBu(I,J)
          vort_xy_dx_smooth(i,J) = DY_dxBu * &
                      (vort_xy_smooth(I,J) * G%IdyCu(I,j) - vort_xy_smooth(I-1,J) * G%IdyCu(I-1,j))
        enddo ; enddo

        do j=js_Kh,je_Kh ; do I=is_Kh-1,ie_Kh
          DX_dyBu = G%dxBu(I,J) * G%IdyBu(I,J)
          vort_xy_dy_smooth(I,j) = DX_dyBu * &
                      (vort_xy_smooth(I,J) * G%IdxCv(i,J) - vort_xy_smooth(I,J-1) * G%IdxCv(i,J-1))
        enddo ; enddo
      endif ! If Leithy

      ! Laplacian of vorticity
      ! if (CS%Leith_Ah .or. CS%use_Leithy) then
      do J=js_Kh-1,je_Kh ; do I=is_Kh-1,ie_Kh
        DY_dxBu = G%dyBu(I,J) * G%IdxBu(I,J)
        DX_dyBu = G%dxBu(I,J) * G%IdyBu(I,J)

        Del2vort_q(I,J) = DY_dxBu * (vort_xy_dx(i+1,J) * G%IdyCv(i+1,J) - vort_xy_dx(i,J) * G%IdyCv(i,J)) + &
                          DX_dyBu * (vort_xy_dy(I,j+1) * G%IdyCu(I,j+1) - vort_xy_dy(I,j) * G%IdyCu(I,j))
      enddo ; enddo
      ! endif

      if (CS%modified_Leith) then

        ! Divergence
        do j=js_Kh-1,je_Kh+1 ; do i=is_Kh-1,ie_Kh+1
          div_xx(i,j) = dudx(i,j) + dvdy(i,j)
        enddo ; enddo

        ! Divergence gradient
        do j=js-1,je+1 ; do I=is_Kh-1,ie_Kh
          div_xx_dx(I,j) = G%IdxCu(I,j)*(div_xx(i+1,j) - div_xx(i,j))
        enddo ; enddo
        do J=js_Kh-1,je_Kh ; do i=is-1,ie+1
          div_xx_dy(i,J) = G%IdyCv(i,J)*(div_xx(i,j+1) - div_xx(i,j))
        enddo ; enddo

        ! Magnitude of divergence gradient
        do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
          grad_div_mag_h(i,j) = sqrt((0.5*(div_xx_dx(I,j) + div_xx_dx(I-1,j)))**2 + &
                                     (0.5*(div_xx_dy(i,J) + div_xx_dy(i,J-1)))**2)
        enddo ; enddo
        do J=js-1,Jeq ; do I=is-1,Ieq
          grad_div_mag_q(I,J) = sqrt((0.5*(div_xx_dx(I,j) + div_xx_dx(I,j+1)))**2 + &
                                     (0.5*(div_xx_dy(i,J) + div_xx_dy(i+1,J)))**2)
        enddo ; enddo

      else

        do j=js-1,je+1 ; do I=is_Kh-1,ie_Kh
          div_xx_dx(I,j) = 0.0
        enddo ; enddo
        do J=js_Kh-1,je_Kh ; do i=is-1,ie+1
          div_xx_dy(i,J) = 0.0
        enddo ; enddo
        do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
          grad_div_mag_h(i,j) = 0.0
        enddo ; enddo
        do J=js-1,Jeq ; do I=is-1,Ieq
          grad_div_mag_q(I,J) = 0.0
        enddo ; enddo

      endif ! CS%modified_Leith

      ! Add in beta for the Leith viscosity
      if (CS%use_beta_in_Leith) then
        do J=js-2,Jeq+1 ; do i=is-1,ie+1
          vort_xy_dx(i,J) = vort_xy_dx(i,J) + 0.5 * ( G%dF_dx(i,j) + G%dF_dx(i,j+1))
        enddo ; enddo
        do j=js-1,je+1 ; do I=is-2,Ieq+1
          vort_xy_dy(I,j) = vort_xy_dy(I,j) + 0.5 * ( G%dF_dy(i,j) + G%dF_dy(i+1,j))
        enddo ; enddo
      endif ! CS%use_beta_in_Leith

      if (CS%use_QG_Leith_visc) then

        do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
          grad_vort_mag_h_2d(i,j) = SQRT((0.5*(vort_xy_dx(i,J) + vort_xy_dx(i,J-1)))**2 + &
                                         (0.5*(vort_xy_dy(I,j) + vort_xy_dy(I-1,j)))**2 )
        enddo ; enddo
        do J=js-1,Jeq ; do I=is-1,Ieq
          grad_vort_mag_q_2d(I,J) = SQRT((0.5*(vort_xy_dx(i,J) + vort_xy_dx(i+1,J)))**2 + &
                                         (0.5*(vort_xy_dy(I,j) + vort_xy_dy(I,j+1)))**2 )
        enddo ; enddo

        ! This accumulates terms, some of which are in VarMix.
        call calc_QG_Leith_viscosity(VarMix, G, GV, US, h, dz, k, div_xx_dx, div_xx_dy, &
                                     slope_x, slope_y, vort_xy_dx, vort_xy_dy)

      endif

      do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
        grad_vort_mag_h(i,j) = SQRT((0.5*(vort_xy_dx(i,J) + vort_xy_dx(i,J-1)))**2 + &
                                    (0.5*(vort_xy_dy(I,j) + vort_xy_dy(I-1,j)))**2 )
      enddo ; enddo
      do J=js-1,Jeq ; do I=is-1,Ieq
        grad_vort_mag_q(I,J) = SQRT((0.5*(vort_xy_dx(i,J) + vort_xy_dx(i+1,J)))**2 + &
                                    (0.5*(vort_xy_dy(I,j) + vort_xy_dy(I,j+1)))**2 )
      enddo ; enddo

      if (CS%use_Leithy) then
        do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
          vert_vort_mag_smooth(i,j) = SQRT((0.5*(vort_xy_dx_smooth(i,J) + &
                                                 vort_xy_dx_smooth(i,J-1)))**2 + &
                                           (0.5*(vort_xy_dy_smooth(I,j) + &
                                                 vort_xy_dy_smooth(I-1,j)))**2 )
        enddo ; enddo
      endif ! Leithy

    endif ! CS%Leith_Kh

    if ((CS%Smagorinsky_Kh) .or. (CS%Smagorinsky_Ah)) then
      do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
        sh_xx_sq = sh_xx(i,j)**2
        sh_xy_sq = 0.25 * ( (sh_xy(I-1,J-1)**2 + sh_xy(I,J)**2) &
                          + (sh_xy(I-1,J)**2 + sh_xy(I,J-1)**2) )
        Shear_mag(i,j) = sqrt(sh_xx_sq + sh_xy_sq)
      enddo ; enddo
    endif

    if (CS%better_bound_Ah .or. CS%better_bound_Kh) then
      do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
        h_min = min(h_u(I,j), h_u(I-1,j), h_v(i,J), h_v(i,J-1))
        hrat_min(i,j) = min(1.0, h_min / (h(i,j,k) + h_neglect))
      enddo ; enddo
    endif

    if (CS%Laplacian) then
      ! Determine the Laplacian viscosity at h points, using the
      ! largest value from several parameterizations. Also get
      ! the Laplacian component of str_xx.

      if ((CS%Leith_Kh) .or. (CS%Leith_Ah) .or. (CS%use_Leithy)) then
        if (CS%use_QG_Leith_visc) then
          do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
            grad_vort = grad_vort_mag_h(i,j) + grad_div_mag_h(i,j)
            grad_vort_qg = 3. * grad_vort_mag_h_2d(i,j)
            vert_vort_mag(i,j) = min(grad_vort, grad_vort_qg)
          enddo ; enddo
        else
          do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
            vert_vort_mag(i,j) = grad_vort_mag_h(i,j) + grad_div_mag_h(i,j)
          enddo ; enddo
        endif
      endif

      ! Static (pre-computed) background viscosity
      do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
        Kh(i,j) = CS%Kh_bg_xx(i,j)
      enddo ; enddo

      ! NOTE: The following do-block can be decomposed and vectorized after the
      !   stack size has been reduced.
      do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
        if (CS%add_LES_viscosity) then
          if (CS%Smagorinsky_Kh) &
            Kh(i,j) = Kh(i,j) + CS%Laplac2_const_xx(i,j) * Shear_mag(i,j)
          if (CS%Leith_Kh) &
            Kh(i,j) = Kh(i,j) + CS%Laplac3_const_xx(i,j) * vert_vort_mag(i,j) * inv_PI3
        else
          if (CS%Smagorinsky_Kh) &
            Kh(i,j) = max(Kh(i,j), CS%Laplac2_const_xx(i,j) * Shear_mag(i,j))
          if (CS%Leith_Kh) &
            Kh(i,j) = max(Kh(i,j), CS%Laplac3_const_xx(i,j) * vert_vort_mag(i,j) * inv_PI3)
        endif
      enddo ; enddo

      ! All viscosity contributions above are subject to resolution scaling

      if (rescale_Kh) then
        do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
          Kh(i,j) = VarMix%Res_fn_h(i,j) * Kh(i,j)
        enddo ; enddo
      endif

      if (legacy_bound) then
        ! Older method of bounding for stability
        do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
          Kh(i,j) = min(Kh(i,j), CS%Kh_Max_xx(i,j))
        enddo ; enddo
      endif

      ! Place a floor on the viscosity, if desired.
      do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
        Kh(i,j) = max(Kh(i,j), CS%Kh_bg_min)
      enddo ; enddo

      if (use_MEKE_Ku) then
        ! *Add* the MEKE contribution (which might be negative)
        if (CS%res_scale_MEKE) then
          do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
            Kh(i,j) = Kh(i,j) + MEKE%Ku(i,j) * VarMix%Res_fn_h(i,j)
          enddo ; enddo
        else
          do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
            Kh(i,j) = Kh(i,j) + MEKE%Ku(i,j)
          enddo ; enddo
        endif
      endif

      if (CS%anisotropic) then
        do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
          ! *Add* the tension component of anisotropic viscosity
          Kh(i,j) = Kh(i,j) + CS%Kh_aniso * (1. - CS%n1n2_h(i,j)**2)
        enddo ; enddo
      endif

      ! Newer method of bounding for stability
      if ((CS%better_bound_Kh) .and. (CS%better_bound_Ah)) then
        do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
          visc_bound_rem(i,j) = 1.0
          Kh_max_here = hrat_min(i,j) * CS%Kh_Max_xx(i,j)
          if (Kh(i,j) >= Kh_max_here) then
            visc_bound_rem(i,j) = 0.0
            Kh(i,j) = Kh_max_here
          elseif ((Kh(i,j) > 0.0) .or. (CS%backscatter_underbound .and. (Kh_max_here > 0.0))) then
            visc_bound_rem(i,j) = 1.0 - Kh(i,j) / Kh_max_here
          endif
        enddo ; enddo
      elseif (CS%better_bound_Kh) then
        do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
          Kh(i,j) = min(Kh(i,j), hrat_min(i,j) * CS%Kh_Max_xx(i,j))
        enddo ; enddo
      endif

      ! In Leith+E parameterization Kh is computed after Ah in the biharmonic loop.
      ! The harmonic component of str_xx is added in the biharmonic loop.
      if (CS%use_Leithy) then
        do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
          Kh(i,j) = 0.
        enddo ; enddo
      endif

      if (CS%id_Kh_h>0 .or. CS%debug) then
        do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
          Kh_h(i,j,k) = Kh(i,j)
        enddo ; enddo
      endif

      if (CS%id_grid_Re_Kh>0) then
        do j=js,je ; do i=is,ie
          KE = 0.125*((u(I,j,k)+u(I-1,j,k))**2 + (v(i,J,k)+v(i,J-1,k))**2)
          grid_Kh = max(Kh(i,j), CS%min_grid_Kh)
          grid_Re_Kh(i,j,k) = (sqrt(KE) * sqrt(CS%grid_sp_h2(i,j))) / grid_Kh
        enddo ; enddo
      endif

      if (CS%id_div_xx_h>0) then
        do j=js,je ; do i=is,ie
          div_xx_h(i,j,k) = dudx(i,j) + dvdy(i,j)
        enddo ; enddo
      endif

      if (CS%id_sh_xx_h>0) then
        do j=js,je ; do i=is,ie
          sh_xx_h(i,j,k) = sh_xx(i,j)
        enddo ; enddo
      endif

      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        str_xx(i,j) = -Kh(i,j) * sh_xx(i,j)
      enddo ; enddo
    else
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        str_xx(i,j) = 0.0
      enddo ; enddo
    endif ! Get Kh at h points and get Laplacian component of str_xx

    if (CS%anisotropic) then
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        ! Shearing-strain averaged to h-points
        local_strain = 0.25 * ( (sh_xy(I,J) + sh_xy(I-1,J-1)) + (sh_xy(I-1,J) + sh_xy(I,J-1)) )
        ! *Add* the shear-strain contribution to the xx-component of stress
        str_xx(i,j) = str_xx(i,j) - CS%Kh_aniso * CS%n1n2_h(i,j) * CS%n1n1_m_n2n2_h(i,j) * local_strain
      enddo ; enddo
    endif

    if (CS%biharmonic) then
      ! Determine the biharmonic viscosity at h points, using the
      ! largest value from several parameterizations. Also get the
      ! biharmonic component of str_xx.
      do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
        Ah(i,j) = CS%Ah_bg_xx(i,j)
      enddo ; enddo

      if ((CS%Smagorinsky_Ah) .or. (CS%Leith_Ah) .or. (CS%use_Leithy)) then
        if (CS%Smagorinsky_Ah) then
          if (CS%bound_Coriolis) then
           do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
              AhSm = Shear_mag(i,j) * (CS%Biharm_const_xx(i,j) &
                  + CS%Biharm_const2_xx(i,j) * Shear_mag(i,j) &
              )
              Ah(i,j) = max(Ah(i,j), AhSm)
            enddo ; enddo
          else
            do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
              AhSm = CS%Biharm_const_xx(i,j) * Shear_mag(i,j)
              Ah(i,j) = max(Ah(i,j), AhSm)
            enddo ; enddo
          endif
        endif

        if (CS%Leith_Ah) then
          do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
            Del2vort_h = 0.25 * ((Del2vort_q(I,J) + Del2vort_q(I-1,J-1)) + &
                                 (Del2vort_q(I-1,J) + Del2vort_q(I,J-1)))
            AhLth = CS%Biharm6_const_xx(i,j) * abs(Del2vort_h) * inv_PI6
            Ah(i,j) = max(Ah(i,j), AhLth)
          enddo ; enddo
        endif

        if (CS%use_Leithy) then
          ! Get m_leithy
          if (CS%smooth_Ah) m_leithy(:,:) = 0.0 ! This is here to initialize domain edge halo values.
          do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
            Del2vort_h = 0.25 * ((Del2vort_q(I,J) + Del2vort_q(I-1,J-1)) + &
                                 (Del2vort_q(I-1,J) + Del2vort_q(I,J-1)))
            AhLth  = CS%Biharm6_const_xx(i,j) * inv_PI6 * abs(Del2vort_h)
            if (AhLth <= CS%Ah_bg_xx(i,j)) then
              m_leithy(i,j) = 0.0
            else
              if ((CS%m_const_leithy(i,j)*vert_vort_mag(i,j)) < abs(vort_xy_smooth(i,j))) then
                m_leithy(i,j) = CS%c_K * (vert_vort_mag(i,j) / vort_xy_smooth(i,j))**2
              else
                m_leithy(i,j) = CS%m_leithy_max(i,j)
              endif
            endif
          enddo ; enddo

          if (CS%smooth_Ah) then
            ! Smooth m_leithy.  A single call smoothes twice.
            call pass_var(m_leithy, G%Domain, halo=2)
            call smooth_x9_h(G, m_leithy, zero_land=.true.)
            call pass_var(m_leithy, G%Domain)
          endif
          ! Get Ah
          do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
            Del2vort_h = 0.25 * ((Del2vort_q(I,J) + Del2vort_q(I-1,J-1)) + &
                                 (Del2vort_q(I-1,J) + Del2vort_q(I,J-1)))
            AhLthy = CS%Biharm6_const_xx(i,j) * inv_PI6 * &
                    sqrt(max(0.,Del2vort_h**2 - m_leithy(i,j)*vert_vort_mag_smooth(i,j)**2))
            Ah(i,j) = max(CS%Ah_bg_xx(i,j), AhLthy)
          enddo ; enddo
          if (CS%smooth_Ah) then
            ! Smooth Ah before applying upper bound.  Square Ah, then smooth, then take its square root.
            Ah_sq(:,:) = 0.0 ! This is here to initialize domain edge halo values.
            do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
              Ah_sq(i,j) = Ah(i,j)**2
            enddo ; enddo
            call pass_var(Ah_sq, G%Domain, halo=2)
            ! A single call smoothes twice.
            call smooth_x9_h(G, Ah_sq, zero_land=.false.)
            call pass_var(Ah_sq, G%Domain)
            do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
              Ah_h(i,j,k) = max(CS%Ah_bg_xx(i,j), sqrt(max(0., Ah_sq(i,j))))
              Ah(i,j)     = Ah_h(i,j,k)
            enddo ; enddo
          else
            do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
              Ah_h(i,j,k) = Ah(i,j)
            enddo ; enddo
          endif
        endif

        if (CS%bound_Ah .and. .not. CS%better_bound_Ah) then
          do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
            Ah(i,j) = min(Ah(i,j), CS%Ah_Max_xx(i,j))
          enddo ; enddo
        endif
      endif ! Smagorinsky_Ah or Leith_Ah or Leith+E

      if (use_MEKE_Au) then
        ! *Add* the MEKE contribution
        do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
          Ah(i,j) = Ah(i,j) + MEKE%Au(i,j)
        enddo ; enddo
      endif

      if (CS%Re_Ah > 0.0) then
        do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
          KE = 0.125*((u(I,j,k)+u(I-1,j,k))**2 + (v(i,J,k)+v(i,J-1,k))**2)
          Ah(i,j) = sqrt(KE) * CS%Re_Ah_const_xx(i,j)
        enddo ; enddo
      endif

      if (CS%better_bound_Ah) then
        if (CS%better_bound_Kh) then
          do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
            Ah(i,j) = min(Ah(i,j), visc_bound_rem(i,j) * hrat_min(i,j) * CS%Ah_Max_xx(i,j))
          enddo ; enddo
        else
          do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
            Ah(i,j) = min(Ah(i,j), hrat_min(i,j) * CS%Ah_Max_xx(i,j))
          enddo ; enddo
        endif
      endif

      if ((CS%id_Ah_h>0) .or. CS%debug .or. CS%use_Leithy) then
        do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
          Ah_h(i,j,k) = Ah(i,j)
        enddo ; enddo
      endif

      if (CS%use_Leithy) then
        ! Compute Leith+E Kh after bounds have been applied to Ah
        ! and after it has been smoothed. Kh = -m_leithy * Ah
        do j=js_Kh,je_Kh ; do i=is_Kh,ie_Kh
          Kh(i,j) = -m_leithy(i,j) * Ah(i,j)
          Kh_h(i,j,k) = Kh(i,j)
        enddo ; enddo
      endif

      if (CS%id_grid_Re_Ah>0) then
        do j=js,je ; do i=is,ie
          KE = 0.125 * ((u(I,j,k) + u(I-1,j,k))**2 + (v(i,J,k) + v(i,J-1,k))**2)
          grid_Ah = max(Ah(i,j), CS%min_grid_Ah)
          grid_Re_Ah(i,j,k) = (sqrt(KE) * CS%grid_sp_h3(i,j)) / grid_Ah
        enddo ; enddo
      endif

      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        d_del2u = G%IdyCu(I,j) * Del2u(I,j) - G%IdyCu(I-1,j) * Del2u(I-1,j)
        d_del2v = G%IdxCv(i,J) * Del2v(i,J) - G%IdxCv(i,J-1) * Del2v(i,J-1)
        d_str = Ah(i,j) * (CS%DY_dxT(i,j) * d_del2u - CS%DX_dyT(i,j) * d_del2v)

        str_xx(i,j) = str_xx(i,j) + d_str

        if (CS%use_Leithy) str_xx(i,j) = str_xx(i,j) - Kh(i,j) * sh_xx_smooth(i,j)

        ! Keep a copy of the biharmonic contribution for backscatter parameterization
        bhstr_xx(i,j) = d_str * (h(i,j,k) * CS%reduction_xx(i,j))
      enddo ; enddo
    endif ! Get biharmonic coefficient at h points and biharmonic part of str_xx

    if (CS%biharmonic) then
      ! Gradient of Laplacian, for use in bi-harmonic term
      do J=js-1,Jeq ; do I=is-1,Ieq
        dDel2vdx(I,J) = CS%DY_dxBu(I,J)*(Del2v(i+1,J)*G%IdyCv(i+1,J) - Del2v(i,J)*G%IdyCv(i,J))
        dDel2udy(I,J) = CS%DX_dyBu(I,J)*(Del2u(I,j+1)*G%IdxCu(I,j+1) - Del2u(I,j)*G%IdxCu(I,j))
      enddo ; enddo
      ! Adjust contributions to shearing strain on open boundaries.
      if (apply_OBC) then ; if (OBC%zero_strain .or. OBC%freeslip_strain) then
        do n=1,OBC%number_of_segments
          J = OBC%segment(n)%HI%JsdB ; I = OBC%segment(n)%HI%IsdB
          if (OBC%segment(n)%is_N_or_S .and. (J >= js-1) .and. (J <= Jeq)) then
            do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
              if (OBC%zero_strain) then
                dDel2vdx(I,J) = 0. ; dDel2udy(I,J) = 0.
              elseif (OBC%freeslip_strain) then
                dDel2udy(I,J) = 0.
              endif
            enddo
          elseif (OBC%segment(n)%is_E_or_W .and. (I >= is-1) .and. (I <= Ieq)) then
            do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
              if (OBC%zero_strain) then
                dDel2vdx(I,J) = 0. ; dDel2udy(I,J) = 0.
              elseif (OBC%freeslip_strain) then
                dDel2vdx(I,J) = 0.
              endif
            enddo
          endif
        enddo
      endif ; endif
    endif

    meke_res_fn = 1.

    if ((CS%Smagorinsky_Kh) .or. (CS%Smagorinsky_Ah)) then
      do J=js-1,Jeq ; do I=is-1,Ieq
        sh_xy_sq = sh_xy(I,J)**2
        sh_xx_sq = 0.25 * ( (sh_xx(i,j)**2 + sh_xx(i+1,j+1)**2) &
                          + (sh_xx(i,j+1)**2 + sh_xx(i+1,j)**2) )
        Shear_mag(I,J) = sqrt(sh_xy_sq + sh_xx_sq)
      enddo ; enddo
    endif

    do J=js-1,Jeq ; do I=is-1,Ieq
      h2uq = 4.0 * (h_u(I,j) * h_u(I,j+1))
      h2vq = 4.0 * (h_v(i,J) * h_v(i+1,J))
      hq(I,J) = (2.0 * (h2uq * h2vq)) &
          / (h_neglect3 + (h2uq + h2vq) * ((h_u(I,j) + h_u(I,j+1)) + (h_v(i,J) + h_v(i+1,J))))
    enddo ; enddo

    if (CS%better_bound_Ah .or. CS%better_bound_Kh) then
      do J=js-1,Jeq ; do I=is-1,Ieq
        h_min = min(h_u(I,j), h_u(I,j+1), h_v(i,J), h_v(i+1,J))
        hrat_min(I,J) = min(1.0, h_min / (hq(I,J) + h_neglect))
      enddo ; enddo

    endif

    if (CS%no_slip) then
      do J=js-1,Jeq ; do I=is-1,Ieq
        if (CS%no_slip .and. (G%mask2dBu(I,J) < 0.5)) then
          if ((G%mask2dCu(I,j) + G%mask2dCu(I,j+1)) + &
              (G%mask2dCv(i,J) + G%mask2dCv(i+1,J)) > 0.0) then
            ! This is a coastal vorticity point, so modify hq and hrat_min.

            hu = G%mask2dCu(I,j) * h_u(I,j) + G%mask2dCu(I,j+1) * h_u(I,j+1)
            hv = G%mask2dCv(i,J) * h_v(i,J) + G%mask2dCv(i+1,J) * h_v(i+1,J)
            if ((G%mask2dCu(I,j) + G%mask2dCu(I,j+1)) * &
                (G%mask2dCv(i,J) + G%mask2dCv(i+1,J)) == 0.0) then
              ! Only one of hu and hv is nonzero, so just add them.
              hq(I,J) = hu + hv
              hrat_min(I,J) = 1.0
            else
              ! Both hu and hv are nonzero, so take the harmonic mean.
              hq(I,J) = 2.0 * (hu * hv) / ((hu + hv) + h_neglect)
              hrat_min(I,J) = min(1.0, min(hu, hv) / (hq(I,J) + h_neglect) )
            endif
          endif
        endif
      enddo ; enddo
    endif

    ! Pass the velocity gradients and thickness to ZB2020
    if (CS%use_ZB2020) then
      call ZB2020_copy_gradient_and_thickness( &
           sh_xx, sh_xy, vort_xy,              &
           hq,                                 &
           G, GV, CS%ZB2020, k)
    endif

    if (CS%Laplacian) then
      ! Determine the Laplacian viscosity at q points, using the
      ! largest value from several parameterizations. Also get the
      ! Laplacian component of str_xy.

      if ((CS%Leith_Kh) .or. (CS%Leith_Ah)) then
        if (CS%use_QG_Leith_visc) then
          do J=js-1,Jeq ; do I=is-1,Ieq
            grad_vort = grad_vort_mag_q(I,J) + grad_div_mag_q(I,J)
            grad_vort_qg = 3. * grad_vort_mag_q_2d(I,J)
            vert_vort_mag(I,J) = min(grad_vort, grad_vort_qg)
          enddo ; enddo
        else
          do J=js-1,Jeq ; do I=is-1,Ieq
            vert_vort_mag(I,J) = grad_vort_mag_q(I,J) + grad_div_mag_q(I,J)
          enddo ; enddo
        endif
      endif

      ! Static (pre-computed) background viscosity
      do J=js-1,Jeq ; do I=is-1,Ieq
        Kh(I,J) = CS%Kh_bg_xy(I,J)
      enddo ; enddo

      if (CS%Smagorinsky_Kh) then
        if (CS%add_LES_viscosity) then
          do J=js-1,Jeq ; do I=is-1,Ieq
            Kh(I,J) = Kh(I,J) + CS%Laplac2_const_xy(I,J) * Shear_mag(I,J)
          enddo ; enddo
        else
          do J=js-1,Jeq ; do I=is-1,Ieq
            Kh(I,J) = max(Kh(I,J), CS%Laplac2_const_xy(I,J) * Shear_mag(I,J) )
          enddo ; enddo
        endif
      endif

      if (CS%Leith_Kh) then
        if (CS%add_LES_viscosity) then
          do J=js-1,Jeq ; do I=is-1,Ieq
            Kh(I,J) = Kh(I,J) + CS%Laplac3_const_xy(I,J) * vert_vort_mag(I,J) * inv_PI3 ! Is this right? -AJA
          enddo ; enddo
        else
          do J=js-1,Jeq ; do I=is-1,Ieq
            Kh(I,J) = max(Kh(I,J), CS%Laplac3_const_xy(I,J) * vert_vort_mag(I,J) * inv_PI3)
          enddo ; enddo
        endif
      endif

      ! All viscosity contributions above are subject to resolution scaling

      ! NOTE: The following do-block can be decomposed and vectorized after the
      !   stack size has been reduced.
      do J=js-1,Jeq ; do I=is-1,Ieq
        if (rescale_Kh) &
          Kh(I,J) = VarMix%Res_fn_q(I,J) * Kh(I,J)

        if (CS%res_scale_MEKE) &
          meke_res_fn = VarMix%Res_fn_q(I,J)

        ! Older method of bounding for stability
        if (legacy_bound) &
          Kh(I,J) = min(Kh(I,J), CS%Kh_Max_xy(I,J))

        Kh(I,J) = max(Kh(I,J), CS%Kh_bg_min) ! Place a floor on the viscosity, if desired.

        if (use_MEKE_Ku) then
          ! *Add* the MEKE contribution (might be negative)
          Kh(I,J) = Kh(I,J) + 0.25*( (MEKE%Ku(i,j) + MEKE%Ku(i+1,j+1)) + &
                           (MEKE%Ku(i+1,j) + MEKE%Ku(i,j+1)) ) * meke_res_fn
        endif

        if (CS%anisotropic) &
          ! *Add* the shear component of anisotropic viscosity
          Kh(I,J) = Kh(I,J) + CS%Kh_aniso * CS%n1n2_q(I,J)**2

        ! Newer method of bounding for stability
        if ((CS%better_bound_Kh) .and. (CS%better_bound_Ah)) then
          visc_bound_rem(I,J) = 1.0
          Kh_max_here = hrat_min(I,J) * CS%Kh_Max_xy(I,J)
          if (Kh(I,J) >= Kh_max_here) then
            visc_bound_rem(I,J) = 0.0
            Kh(I,J) = Kh_max_here
          elseif ((Kh(I,J) > 0.0) .or. (CS%backscatter_underbound .and. (Kh_max_here > 0.0))) then
            visc_bound_rem(I,J) = 1.0 - Kh(I,J) / Kh_max_here
          endif
        elseif (CS%better_bound_Kh) then
          Kh(I,J) = min(Kh(I,J), hrat_min(I,J) * CS%Kh_Max_xy(I,J))
        endif

        ! Leith+E doesn't recompute Kh at q points, it just interpolates it from h to q points
        if (CS%use_Leithy) then
          Kh(I,J) = 0.25 * ((Kh_h(i,j,k) + Kh_h(i+1,j+1,k)) + (Kh_h(i,j+1,k) + Kh_h(i+1,j,k)))
        end if

        if (CS%id_Kh_q>0 .or. CS%debug) &
          Kh_q(I,J,k) = Kh(I,J)

        if (CS%id_vort_xy_q>0) &
          vort_xy_q(I,J,k) = vort_xy(I,J)

        if (CS%id_sh_xy_q>0) &
          sh_xy_q(I,J,k) = sh_xy(I,J)
      enddo ; enddo

      if ( .not. CS%use_Leithy) then
        do J=js-1,Jeq ; do I=is-1,Ieq
          str_xy(I,J) = -Kh(I,J) * sh_xy(I,J)
        enddo ; enddo
      else
        do J=js-1,Jeq ; do I=is-1,Ieq
          str_xy(I,J) = -Kh(I,J) * sh_xy_smooth(I,J)
        enddo ; enddo
      endif
    else
      do J=js-1,Jeq ; do I=is-1,Ieq
        str_xy(I,J) = 0.
      enddo ; enddo
    endif ! get harmonic coefficient Kh at q points and harmonic part of str_xy

    if (CS%anisotropic) then
      do J=js-1,Jeq ; do I=is-1,Ieq
        ! Horizontal-tension averaged to q-points
        local_strain = 0.25 * ( (sh_xx(i,j) + sh_xx(i+1,j+1)) + (sh_xx(i+1,j) + sh_xx(i,j+1)) )
        ! *Add* the tension contribution to the xy-component of stress
        str_xy(I,J) = str_xy(I,J) - CS%Kh_aniso * CS%n1n2_q(I,J) * CS%n1n1_m_n2n2_q(I,J) * local_strain
      enddo ; enddo
    endif

    if (CS%biharmonic) then
      ! Determine the biharmonic viscosity at q points, using the
      ! largest value from several parameterizations. Also get the
      ! biharmonic component of str_xy.
      do J=js-1,Jeq ; do I=is-1,Ieq
        Ah(I,J) = CS%Ah_bg_xy(I,J)
      enddo ; enddo

      if (CS%Smagorinsky_Ah .or. CS%Leith_Ah) then
        if (CS%Smagorinsky_Ah) then
          if (CS%bound_Coriolis) then
            do J=js-1,Jeq ; do I=is-1,Ieq
              AhSm = Shear_mag(I,J) * (CS%Biharm_const_xy(I,J) &
                  + CS%Biharm_const2_xy(I,J) * Shear_mag(I,J) &
              )
              Ah(I,J) = max(Ah(I,J), AhSm)
            enddo ; enddo
          else
            do J=js-1,Jeq ; do I=is-1,Ieq
              AhSm = CS%Biharm_const_xy(I,J) * Shear_mag(I,J)
              Ah(I,J) = max(Ah(I,J), AhSm)
            enddo ; enddo
          endif
        endif

        if (CS%Leith_Ah) then
          do J=js-1,Jeq ; do I=is-1,Ieq
            AhLth = CS%Biharm6_const_xy(I,J) * abs(Del2vort_q(I,J)) * inv_PI6
            Ah(I,J) = max(Ah(I,J), AhLth)
          enddo ; enddo
        endif

        if (CS%bound_Ah .and. .not.CS%better_bound_Ah) then
          do J=js-1,Jeq ; do I=is-1,Ieq
            Ah(I,J) = min(Ah(I,J), CS%Ah_Max_xy(I,J))
          enddo ; enddo
        endif
      endif ! Smagorinsky_Ah or Leith_Ah

      if (use_MEKE_Au) then
        ! *Add* the MEKE contribution
        do J=js-1,Jeq ; do I=is-1,Ieq
          Ah(I,J) = Ah(I,J) + 0.25 * ( &
              (MEKE%Au(i,j) + MEKE%Au(i+1,j+1)) + (MEKE%Au(i+1,j) + MEKE%Au(i,j+1)) &
          )
        enddo ; enddo
      endif

      if (CS%Re_Ah > 0.0) then
        do J=js-1,Jeq ; do I=is-1,Ieq
          KE = 0.125 * ((u(I,j,k) + u(I,j+1,k))**2 + (v(i,J,k) + v(i+1,J,k))**2)
          Ah(I,J) = sqrt(KE) * CS%Re_Ah_const_xy(I,J)
        enddo ; enddo
      endif

      if (CS%better_bound_Ah) then
        if (CS%better_bound_Kh) then
          do J=js-1,Jeq ; do I=is-1,Ieq
            Ah(I,J) = min(Ah(I,J), visc_bound_rem(I,J) * hrat_min(I,J) * CS%Ah_Max_xy(I,J))
          enddo ; enddo
        else
          do J=js-1,Jeq ; do I=is-1,Ieq
            Ah(I,J) = min(Ah(I,J), hrat_min(I,J) * CS%Ah_Max_xy(I,J))
          enddo ; enddo
        endif
      endif

      ! Leith+E doesn't recompute Ah at q points, it just interpolates it from h to q points
      if (CS%use_Leithy) then
        do J=js-1,Jeq ; do I=is-1,Ieq
          Ah(I,J) = 0.25 * ((Ah_h(i,j,k) + Ah_h(i+1,j+1,k)) + (Ah_h(i,j+1,k) + Ah_h(i+1,j,k)))
        enddo ; enddo
      end if

      if (CS%id_Ah_q>0 .or. CS%debug) then
        do J=js-1,Jeq ; do I=is-1,Ieq
          Ah_q(I,J,k) = Ah(I,J)
        enddo ; enddo
      endif

      ! Again, need to initialize str_xy as if its biharmonic
      do J=js-1,Jeq ; do I=is-1,Ieq
        d_str = Ah(I,J) * (dDel2vdx(I,J) + dDel2udy(I,J))

        str_xy(I,J) = str_xy(I,J) + d_str

        ! Keep a copy of the biharmonic contribution for backscatter parameterization
        bhstr_xy(I,J) = d_str * (hq(I,J) * G%mask2dBu(I,J) * CS%reduction_xy(I,J))
      enddo ; enddo
    endif ! Get Ah at q points and biharmonic part of str_xy

    if (CS%use_GME) then
      ! The wider halo here is to permit one pass of smoothing without a halo update.
      do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
        GME_coeff = GME_effic_h(i,j) * 0.25 * &
            ((KH_u_GME(I,j,k)+KH_u_GME(I-1,j,k)) + (KH_v_GME(i,J,k)+KH_v_GME(i,J-1,k)))
        GME_coeff = MIN(GME_coeff, CS%GME_limiter)

        if ((CS%id_GME_coeff_h>0) .or. find_FrictWork) GME_coeff_h(i,j,k) = GME_coeff
        str_xx_GME(i,j) = GME_coeff * sh_xx_bt(i,j)
      enddo ; enddo

      ! The wider halo here is to permit one pass of smoothing without a halo update.
      do J=js-2,je+1 ; do I=is-2,ie+1
        GME_coeff = GME_effic_q(I,J) * 0.25 * &
            ((KH_u_GME(I,j,k)+KH_u_GME(I,j+1,k)) + (KH_v_GME(i,J,k)+KH_v_GME(i+1,J,k)))
        GME_coeff = MIN(GME_coeff, CS%GME_limiter)

        if (CS%id_GME_coeff_q>0) GME_coeff_q(I,J,k) = GME_coeff
        str_xy_GME(I,J) = GME_coeff * sh_xy_bt(I,J)
      enddo ; enddo

      ! Applying GME diagonal term.  This is linear and the arguments can be rescaled.
      call smooth_GME(CS, G, GME_flux_h=str_xx_GME)
      call smooth_GME(CS, G, GME_flux_q=str_xy_GME)

      ! This changes the units of str_xx from [L2 T-2 ~> m2 s-2] to [H L2 T-2 ~> m3 s-2 or kg s-2].
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        str_xx(i,j) = (str_xx(i,j) + str_xx_GME(i,j)) * (h(i,j,k) * CS%reduction_xx(i,j))
      enddo ; enddo

      ! This adds in GME and changes the units of str_xx from [L2 T-2 ~> m2 s-2] to [H L2 T-2 ~> m3 s-2 or kg s-2].
      if (CS%no_slip) then
        do J=js-1,Jeq ; do I=is-1,Ieq
          str_xy(I,J) = (str_xy(I,J) + str_xy_GME(I,J)) * (hq(I,J) * CS%reduction_xy(I,J))
        enddo ; enddo
      else
        do J=js-1,Jeq ; do I=is-1,Ieq
          str_xy(I,J) = (str_xy(I,J) + str_xy_GME(I,J)) * (hq(I,J) * G%mask2dBu(I,J) * CS%reduction_xy(I,J))
        enddo ; enddo
      endif

    else ! .not. use_GME
      ! This changes the units of str_xx from [L2 T-2 ~> m2 s-2] to [H L2 T-2 ~> m3 s-2 or kg s-2].
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        str_xx(i,j) = str_xx(i,j) * (h(i,j,k) * CS%reduction_xx(i,j))
      enddo ; enddo

      ! This changes the units of str_xy from [L2 T-2 ~> m2 s-2] to [H L2 T-2 ~> m3 s-2 or kg s-2].
      if (CS%no_slip) then
        do J=js-1,Jeq ; do I=is-1,Ieq
          str_xy(I,J) = str_xy(I,J) * (hq(I,J) * CS%reduction_xy(I,J))
        enddo ; enddo
      else
        do J=js-1,Jeq ; do I=is-1,Ieq
          str_xy(I,J) = str_xy(I,J) * (hq(I,J) * G%mask2dBu(I,J) * CS%reduction_xy(I,J))
        enddo ; enddo
      endif
    endif ! use_GME

    ! Evaluate 1/h x.Div(h Grad u) or the biharmonic equivalent.
    do j=js,je ; do I=Isq,Ieq
      diffu(I,j,k) = ((G%IdyCu(I,j)*(CS%dy2h(i,j)*str_xx(i,j) - CS%dy2h(i+1,j)*str_xx(i+1,j)) + &
                       G%IdxCu(I,j)*(CS%dx2q(I,J-1)*str_xy(I,J-1) - CS%dx2q(I,J)*str_xy(I,J))) * &
                     G%IareaCu(I,j)) / (h_u(I,j) + h_neglect)
    enddo ; enddo

    if (apply_OBC) then
      ! This is not the right boundary condition. If all the masking of tendencies are done
      ! correctly later then eliminating this block should not change answers.
      do n=1,OBC%number_of_segments
        if (OBC%segment(n)%is_E_or_W) then
          I = OBC%segment(n)%HI%IsdB
          do j=OBC%segment(n)%HI%jsd,OBC%segment(n)%HI%jed
            diffu(I,j,k) = 0.
          enddo
        endif
      enddo
    endif

    ! Evaluate 1/h y.Div(h Grad u) or the biharmonic equivalent.
    do J=Jsq,Jeq ; do i=is,ie
      diffv(i,J,k) = ((G%IdyCv(i,J)*(CS%dy2q(I-1,J)*str_xy(I-1,J) - CS%dy2q(I,J)*str_xy(I,J)) - &
                       G%IdxCv(i,J)*(CS%dx2h(i,j)*str_xx(i,j) - CS%dx2h(i,j+1)*str_xx(i,j+1))) * &
                     G%IareaCv(i,J)) / (h_v(i,J) + h_neglect)
    enddo ; enddo

    if (apply_OBC) then
      ! This is not the right boundary condition. If all the masking of tendencies are done
      ! correctly later then eliminating this block should not change answers.
      do n=1,OBC%number_of_segments
        if (OBC%segment(n)%is_N_or_S) then
          J = OBC%segment(n)%HI%JsdB
          do i=OBC%segment(n)%HI%isd,OBC%segment(n)%HI%ied
            diffv(i,J,k) = 0.
          enddo
        endif
      enddo
    endif

    if (find_FrictWork) then
      if (CS%FrictWork_bug) then ; do j=js,je ; do i=is,ie
      ! Diagnose   str_xx*d_x u - str_yy*d_y v + str_xy*(d_y u + d_x v)
      ! This is the old formulation that includes energy diffusion
        FrictWork(i,j,k) = GV%H_to_RZ * ( &
              (str_xx(i,j) * (u(I,j,k)-u(I-1,j,k))*G%IdxT(i,j)    &
             - str_xx(i,j) * (v(i,J,k)-v(i,J-1,k))*G%IdyT(i,j))   &
          + 0.25*((str_xy(I,J) *                                  &
                   ((u(I,j+1,k)-u(I,j,k))*G%IdyBu(I,J)            &
                  + (v(i+1,J,k)-v(i,J,k))*G%IdxBu(I,J))           &
                 + str_xy(I-1,J-1) *                              &
                   ((u(I-1,j,k)-u(I-1,j-1,k))*G%IdyBu(I-1,J-1)    &
                  + (v(i,J-1,k)-v(i-1,J-1,k))*G%IdxBu(I-1,J-1)) ) &
                + (str_xy(I-1,J) *                                &
                   ((u(I-1,j+1,k)-u(I-1,j,k))*G%IdyBu(I-1,J)      &
                  + (v(i,J,k)-v(i-1,J,k))*G%IdxBu(I-1,J))         &
                 + str_xy(I,J-1) *                                &
                   ((u(I,j,k)-u(I,j-1,k))*G%IdyBu(I,J-1)          &
                  + (v(i+1,J-1,k)-v(i,J-1,k))*G%IdxBu(I,J-1)) ) ) )
          enddo ; enddo
      else ; do j=js,je ; do i=is,ie
        FrictWork(i,j,k) = GV%H_to_RZ * G%IareaT(i,j) * ( &
            ((str_xx(i,j)*CS%dy2h(i,j) * ( &
                  (uh(I,j,k)*G%dxCu(I,j)*G%IdyCu(I,j)*G%IareaCu(I,j)/(h_u(I,j)+h_neglect)) &
                - (uh(I-1,j,k)*G%dxCu(I-1,j)*G%IdyCu(I-1,j)*G%IareaCu(I-1,j)/(h_u(I-1,j)+h_neglect)) ) ) &
           - (str_xx(i,j)*CS%dx2h(i,j) * ( &
                  (vh(i,J,k)*G%dyCv(i,J)*G%IdxCv(i,J)*G%IareaCv(i,J)/(h_v(i,J)+h_neglect)) &
                - (vh(i,J-1,k)*G%dyCv(i,J-1)*G%IdxCv(i,J-1)*G%IareaCv(i,J-1)/(h_v(i,J-1)+h_neglect)) ) )) &
       + (0.25*(((str_xy(I,J)*(                                     &
                     (CS%dx2q(I,J)*((uh(I,j+1,k)*G%IareaCu(I,j+1)/(h_u(I,j+1)+h_neglect)) &
                                  - (uh(I,j,k)*G%IareaCu(I,j)/(h_u(I,j)+h_neglect))))            &
                   + (CS%dy2q(I,J)*((vh(i+1,J,k)*G%IareaCv(i+1,J)/(h_v(i+1,J)+h_neglect)) &
                                  - (vh(i,J,k)*G%IareaCv(i,J)/(h_v(i,J)+h_neglect)))) ))          &
                +(str_xy(I-1,J-1)*(                                 &
                     (CS%dx2q(I-1,J-1)*((uh(I-1,j,k)*G%IareaCu(I-1,j)/(h_u(I-1,j)+h_neglect)) &
                                      - (uh(I-1,j-1,k)*G%IareaCu(I-1,j-1)/(h_u(I-1,j-1)+h_neglect))))    &
                   + (CS%dy2q(I-1,J-1)*((vh(i,J-1,k)*G%IareaCv(i,J-1)/(h_v(i,J-1)+h_neglect)) &
                                      - (vh(i-1,J-1,k)*G%IareaCv(i-1,J-1)/(h_v(i-1,J-1)+h_neglect)))) )) ) &
               +((str_xy(I-1,J)*(                                   &
                     (CS%dx2q(I-1,J)*((uh(I-1,j+1,k)*G%IareaCu(I-1,j+1)/(h_u(I-1,j+1)+h_neglect)) &
                                    - (uh(I-1,j,k)*G%IareaCu(I-1,j)/(h_u(I-1,j)+h_neglect))))      &
                   + (CS%dy2q(I-1,J)*((vh(i,J,k)*G%IareaCv(i,J)/(h_v(i,J)+h_neglect)) &
                                    - (vh(i-1,J,k)*G%IareaCv(i-1,J)/(h_v(i-1,J)+h_neglect)))) ))        &
                +(str_xy(I,J-1)*(                                   &
                     (CS%dx2q(I,J-1)*((uh(I,j,k)*G%IareaCu(I,j)/(h_u(I,j)+h_neglect)) &
                                    - (uh(I,j-1,k)*G%IareaCu(I,j-1)/(h_u(I,j-1)+h_neglect))))          &
                   + (CS%dy2q(I,J-1)*((vh(i+1,J-1,k)*G%IareaCv(i+1,J-1)/(h_v(i+1,J-1)+h_neglect)) &
                                    - (vh(i,J-1,k)*G%IareaCv(i,J-1)/(h_v(i,J-1)+h_neglect)))) )) ) )) )
      enddo ; enddo ; endif
    endif


    if (CS%use_GME) then
      if (CS%FrictWork_bug) then ; do j=js,je ; do i=is,ie
      ! Diagnose   str_xx_GME*d_x u - str_yy_GME*d_y v + str_xy_GME*(d_y u + d_x v)
      ! This is the old formulation that includes energy diffusion
        FrictWork_GME(i,j,k) = GV%H_to_RZ * ( &
              (str_xx_GME(i,j)*(u(I,j,k)-u(I-1,j,k))*G%IdxT(i,j)     &
             - str_xx_GME(i,j)*(v(i,J,k)-v(i,J-1,k))*G%IdyT(i,j))    &
            + 0.25*((str_xy_GME(I,J) *                               &
                     ((u(I,j+1,k)-u(I,j,k))*G%IdyBu(I,J)             &
                    + (v(i+1,J,k)-v(i,J,k))*G%IdxBu(I,J))            &
                   + str_xy_GME(I-1,J-1) *                           &
                     ((u(I-1,j,k)-u(I-1,j-1,k))*G%IdyBu(I-1,J-1)     &
                    + (v(i,J-1,k)-v(i-1,J-1,k))*G%IdxBu(I-1,J-1)) )  &
                  + (str_xy_GME(I-1,J) *                             &
                     ((u(I-1,j+1,k)-u(I-1,j,k))*G%IdyBu(I-1,J)       &
                    + (v(i,J,k)-v(i-1,J,k))*G%IdxBu(I-1,J))          &
                   + str_xy_GME(I,J-1) *                             &
                     ((u(I,j,k)-u(I,j-1,k))*G%IdyBu(I,J-1)           &
                    + (v(i+1,J-1,k)-v(i,J-1,k))*G%IdxBu(I,J-1)) ) ) )
        enddo ; enddo
      else ; do j=js,je ; do i=is,ie
        FrictWork_GME(i,j,k) = GV%H_to_RZ * G%IareaT(i,j) * ( &
            ((str_xx_GME(i,j)*CS%dy2h(i,j) * ( &
                  (uh(I,j,k)*G%dxCu(I,j)*G%IdyCu(I,j)*G%IareaCu(I,j)/(h_u(I,j)+h_neglect)) &
                - (uh(I-1,j,k)*G%dxCu(I-1,j)*G%IdyCu(I-1,j)*G%IareaCu(I-1,j)/(h_u(I-1,j)+h_neglect)) ) ) &
           - (str_xx_GME(i,j)*CS%dx2h(i,j) * ( &
                  (vh(i,J,k)*G%dyCv(i,J)*G%IdxCv(i,J)*G%IareaCv(i,J)/(h_v(i,J)+h_neglect)) &
                - (vh(i,J-1,k)*G%dyCv(i,J-1)*G%IdxCv(i,J-1)*G%IareaCv(i,J-1)/(h_v(i,J-1)+h_neglect)) ) )) &
       + (0.25*(((str_xy_GME(I,J)*(                                     &
                     (CS%dx2q(I,J)*((uh(I,j+1,k)*G%IareaCu(I,j+1)/(h_u(I,j+1)+h_neglect)) &
                                  - (uh(I,j,k)*G%IareaCu(I,j)/(h_u(I,j)+h_neglect))))            &
                   + (CS%dy2q(I,J)*((vh(i+1,J,k)*G%IareaCv(i+1,J)/(h_v(i+1,J)+h_neglect)) &
                                  - (vh(i,J,k)*G%IareaCv(i,J)/(h_v(i,J)+h_neglect)))) ))          &
                +(str_xy_GME(I-1,J-1)*(                                 &
                     (CS%dx2q(I-1,J-1)*((uh(I-1,j,k)*G%IareaCu(I-1,j)/(h_u(I-1,j)+h_neglect)) &
                                      - (uh(I-1,j-1,k)*G%IareaCu(I-1,j-1)/(h_u(I-1,j-1)+h_neglect))))    &
                   + (CS%dy2q(I-1,J-1)*((vh(i,J-1,k)*G%IareaCv(i,J-1)/(h_v(i,J-1)+h_neglect)) &
                                      - (vh(i-1,J-1,k)*G%IareaCv(i-1,J-1)/(h_v(i-1,J-1)+h_neglect)))) )) ) &
               +((str_xy_GME(I-1,J)*(                                   &
                     (CS%dx2q(I-1,J)*((uh(I-1,j+1,k)*G%IareaCu(I-1,j+1)/(h_u(I-1,j+1)+h_neglect)) &
                                    - (uh(I-1,j,k)*G%IareaCu(I-1,j)/(h_u(I-1,j)+h_neglect))))      &
                   + (CS%dy2q(I-1,J)*((vh(i,J,k)*G%IareaCv(i,J)/(h_v(i,J)+h_neglect)) &
                                    - (vh(i-1,J,k)*G%IareaCv(i-1,J)/(h_v(i-1,J)+h_neglect)))) ))        &
                +(str_xy_GME(I,J-1)*(                                   &
                     (CS%dx2q(I,J-1)*((uh(I,j,k)*G%IareaCu(I,j)/(h_u(I,j)+h_neglect)) &
                                    - (uh(I,j-1,k)*G%IareaCu(I,j-1)/(h_u(I,j-1)+h_neglect))))          &
                   + (CS%dy2q(I,J-1)*((vh(i+1,J-1,k)*G%IareaCv(i+1,J-1)/(h_v(i+1,J-1)+h_neglect)) &
                                    - (vh(i,J-1,k)*G%IareaCv(i,J-1)/(h_v(i,J-1)+h_neglect)))) )) ) )) )

      enddo ; enddo ; endif
    endif

    ! Make a similar calculation as for FrictWork above but accumulating into
    ! the vertically integrated MEKE source term, and adjusting for any
    ! energy loss seen as a reduction in the (biharmonic) frictional source term.
    if (find_FrictWork .and. allocated(MEKE%mom_src)) then
      if (k==1) then
        do j=js,je ; do i=is,ie
          MEKE%mom_src(i,j) = 0.
        enddo ; enddo
        if (allocated(MEKE%GME_snk)) then
          do j=js,je ; do i=is,ie
            MEKE%GME_snk(i,j) = 0.
          enddo ; enddo
        endif
      endif
      if (MEKE%backscatter_Ro_c /= 0.) then
        do j=js,je ; do i=is,ie
          FatH = 0.25*( (abs(G%CoriolisBu(I-1,J-1)) + abs(G%CoriolisBu(I,J))) + &
                        (abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J-1))) )
          Shear_mag_bc = sqrt(sh_xx(i,j) * sh_xx(i,j) + &
            0.25*((sh_xy(I-1,J-1)*sh_xy(I-1,J-1) + sh_xy(I,J)*sh_xy(I,J)) + &
                  (sh_xy(I-1,J)*sh_xy(I-1,J) + sh_xy(I,J-1)*sh_xy(I,J-1))))
          if (CS%answer_date > 20190101) then
            FatH = (US%s_to_T*FatH)**MEKE%backscatter_Ro_pow ! f^n
            ! Note the hard-coded dimensional constant in the following line that can not
            ! be rescaled for dimensional consistency.
            Shear_mag_bc = (((US%s_to_T * Shear_mag_bc)**MEKE%backscatter_Ro_pow) + 1.e-30) &
                        * MEKE%backscatter_Ro_c ! c * D^n
            ! The Rossby number function is g(Ro) = 1/(1+c.Ro^n)
            ! RoScl = 1 - g(Ro)
            RoScl = Shear_mag_bc / (FatH + Shear_mag_bc) ! = 1 - f^n/(f^n+c*D^n)
          else
            if (FatH <= backscat_subround*Shear_mag_bc) then
              RoScl = 1.0
            else
              Sh_F_pow = MEKE%backscatter_Ro_c * (Shear_mag_bc / FatH)**MEKE%backscatter_Ro_pow
              RoScl = Sh_F_pow / (1.0 + Sh_F_pow) ! = 1 - f^n/(f^n+c*D^n)
            endif
          endif

          MEKE%mom_src(i,j) = MEKE%mom_src(i,j) + GV%H_to_RZ * ( &
                ((str_xx(i,j)-RoScl*bhstr_xx(i,j))*(u(I,j,k)-u(I-1,j,k))*G%IdxT(i,j)  &
                -(str_xx(i,j)-RoScl*bhstr_xx(i,j))*(v(i,J,k)-v(i,J-1,k))*G%IdyT(i,j)) &
              + 0.25*(((str_xy(I,J)-RoScl*bhstr_xy(I,J)) *                            &
                       ((u(I,j+1,k)-u(I,j,k))*G%IdyBu(I,J)                            &
                      + (v(i+1,J,k)-v(i,J,k))*G%IdxBu(I,J) )                          &
                     + (str_xy(I-1,J-1)-RoScl*bhstr_xy(I-1,J-1)) *                    &
                       ((u(I-1,j,k)-u(I-1,j-1,k))*G%IdyBu(I-1,J-1)                    &
                      + (v(i,J-1,k)-v(i-1,J-1,k))*G%IdxBu(I-1,J-1)) )                 &
                    + ((str_xy(I-1,J)-RoScl*bhstr_xy(I-1,J)) *                        &
                       ((u(I-1,j+1,k)-u(I-1,j,k))*G%IdyBu(I-1,J)                      &
                      + (v(i,J,k)-v(i-1,J,k))*G%IdxBu(I-1,J))                         &
                     + (str_xy(I,J-1)-RoScl*bhstr_xy(I,J-1)) *                        &
                       ((u(I,j,k)-u(I,j-1,k))*G%IdyBu(I,J-1)                          &
                      + (v(i+1,J-1,k)-v(i,J-1,k))*G%IdxBu(I,J-1)) ) ) )
        enddo ; enddo
      endif ! MEKE%backscatter_Ro_c

      do j=js,je ; do i=is,ie
        MEKE%mom_src(i,j) = MEKE%mom_src(i,j) + FrictWork(i,j,k)
      enddo ; enddo

      if (CS%use_GME .and. allocated(MEKE%GME_snk)) then
        do j=js,je ; do i=is,ie
          MEKE%GME_snk(i,j) = MEKE%GME_snk(i,j) + FrictWork_GME(i,j,k)
        enddo ; enddo
      endif

    endif ! find_FrictWork and associated(mom_src)

  enddo ! end of k loop

  ! Offer fields for diagnostic averaging.
  if (CS%id_normstress > 0) call post_data(CS%id_normstress, NoSt, CS%diag)
  if (CS%id_shearstress > 0) call post_data(CS%id_shearstress, ShSt, CS%diag)
  if (CS%id_diffu>0)     call post_data(CS%id_diffu, diffu, CS%diag)
  if (CS%id_diffv>0)     call post_data(CS%id_diffv, diffv, CS%diag)
  if (CS%id_FrictWork>0) call post_data(CS%id_FrictWork, FrictWork, CS%diag)
  if (CS%id_Ah_h>0)      call post_data(CS%id_Ah_h, Ah_h, CS%diag)
  if (CS%id_grid_Re_Ah>0) call post_data(CS%id_grid_Re_Ah, grid_Re_Ah, CS%diag)
  if (CS%id_div_xx_h>0)  call post_data(CS%id_div_xx_h, div_xx_h, CS%diag)
  if (CS%id_vort_xy_q>0) call post_data(CS%id_vort_xy_q, vort_xy_q, CS%diag)
  if (CS%id_sh_xx_h>0)   call post_data(CS%id_sh_xx_h, sh_xx_h, CS%diag)
  if (CS%id_sh_xy_q>0)   call post_data(CS%id_sh_xy_q, sh_xy_q, CS%diag)
  if (CS%id_Ah_q>0)      call post_data(CS%id_Ah_q, Ah_q, CS%diag)
  if (CS%id_Kh_h>0)      call post_data(CS%id_Kh_h, Kh_h, CS%diag)
  if (CS%id_grid_Re_Kh>0) call post_data(CS%id_grid_Re_Kh, grid_Re_Kh, CS%diag)
  if (CS%id_Kh_q>0)      call post_data(CS%id_Kh_q, Kh_q, CS%diag)
  if (CS%use_GME) then  ! post barotropic tension and strain
    if (CS%id_GME_coeff_h > 0) call post_data(CS%id_GME_coeff_h, GME_coeff_h, CS%diag)
    if (CS%id_GME_coeff_q > 0) call post_data(CS%id_GME_coeff_q, GME_coeff_q, CS%diag)
    if (CS%id_FrictWork_GME>0) call post_data(CS%id_FrictWork_GME, FrictWork_GME, CS%diag)
    if (CS%id_dudx_bt > 0) call post_data(CS%id_dudx_bt, dudx_bt, CS%diag)
    if (CS%id_dvdy_bt > 0) call post_data(CS%id_dvdy_bt, dvdy_bt, CS%diag)
    if (CS%id_dudy_bt > 0) call post_data(CS%id_dudy_bt, dudy_bt, CS%diag)
    if (CS%id_dvdx_bt > 0) call post_data(CS%id_dvdx_bt, dvdx_bt, CS%diag)
  endif

  if (CS%debug) then
    if (CS%Laplacian) then
      call hchksum(Kh_h, "Kh_h", G%HI, haloshift=0, unscale=US%L_to_m**2*US%s_to_T)
      call Bchksum(Kh_q, "Kh_q", G%HI, haloshift=0, unscale=US%L_to_m**2*US%s_to_T)
    endif
    if (CS%biharmonic) call hchksum(Ah_h, "Ah_h", G%HI, haloshift=0, unscale=US%L_to_m**4*US%s_to_T)
    if (CS%biharmonic) call Bchksum(Ah_q, "Ah_q", G%HI, haloshift=0, unscale=US%L_to_m**4*US%s_to_T)
  endif

  if (CS%id_FrictWorkIntz > 0) then
    do j=js,je
      do i=is,ie ; FrictWorkIntz(i,j) = FrictWork(i,j,1) ; enddo
      do k=2,nz ; do i=is,ie
        FrictWorkIntz(i,j) = FrictWorkIntz(i,j) + FrictWork(i,j,k)
      enddo ; enddo
    enddo
    call post_data(CS%id_FrictWorkIntz, FrictWorkIntz, CS%diag)
  endif

  if (present(ADp)) then
    ! Diagnostics of the fractional thicknesses times momentum budget terms
    ! 3D diagnostics of hf_diffu(diffv) are commented because there is no clarity on proper remapping grid option.
    ! The code is retained for debugging purposes in the future.
    !if (CS%id_hf_diffu > 0) call post_product_u(CS%id_hf_diffu, diffu, ADp%diag_hfrac_u, G, nz, CS%diag)
    !if (CS%id_hf_diffv > 0) call post_product_v(CS%id_hf_diffv, diffv, ADp%diag_hfrac_v, G, nz, CS%diag)

    ! Diagnostics for thickness-weighted vertically averaged momentum budget terms
    if (CS%id_hf_diffu_2d > 0) call post_product_sum_u(CS%id_hf_diffu_2d, diffu, ADp%diag_hfrac_u, G, nz, CS%diag)
    if (CS%id_hf_diffv_2d > 0) call post_product_sum_v(CS%id_hf_diffv_2d, diffv, ADp%diag_hfrac_v, G, nz, CS%diag)

    ! Diagnostics for the vertical sum of layer thickness x momentum budget terms
    if (CS%id_intz_diffu_2d > 0) call post_product_sum_u(CS%id_intz_diffu_2d, diffu, ADp%diag_hu, G, nz, CS%diag)
    if (CS%id_intz_diffv_2d > 0) call post_product_sum_v(CS%id_intz_diffv_2d, diffv, ADp%diag_hv, G, nz, CS%diag)

    ! Diagnostics for thickness x momentum budget terms
    if (CS%id_h_diffu > 0) call post_product_u(CS%id_h_diffu, diffu, ADp%diag_hu, G, nz, CS%diag)
    if (CS%id_h_diffv > 0) call post_product_v(CS%id_h_diffv, diffv, ADp%diag_hv, G, nz, CS%diag)

    ! Diagnostics for momentum budget terms multiplied by visc_rem_[uv],
    if (CS%id_diffu_visc_rem > 0) call post_product_u(CS%id_diffu_visc_rem, diffu, ADp%visc_rem_u, G, nz, CS%diag)
    if (CS%id_diffv_visc_rem > 0) call post_product_v(CS%id_diffv_visc_rem, diffv, ADp%visc_rem_v, G, nz, CS%diag)
  endif

  if (CS%use_ZB2020) then
    call ZB2020_lateral_stress(u, v, h, diffu, diffv, G, GV, CS%ZB2020, &
                               CS%dx2h, CS%dy2h, CS%dx2q, CS%dy2q)
  endif

end subroutine horizontal_viscosity

!> Allocates space for and calculates static variables used by horizontal_viscosity().
!! hor_visc_init calculates and stores the values of a number of metric functions that
!! are used in horizontal_viscosity().
subroutine hor_visc_init(Time, G, GV, US, param_file, diag, CS, ADp)
  type(time_type),         intent(in)    :: Time !< Current model time.
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< Structure to regulate diagnostic output.
  type(hor_visc_CS),       intent(inout) :: CS   !< Horizontal viscosity control structure
  type(accel_diag_ptrs), intent(in), optional :: ADp !< Acceleration diagnostics

  ! u0v is the Laplacian sensitivities to the v velocities at u points, with u0u, v0u, and v0v defined analogously.
  real, dimension(SZIB_(G),SZJ_(G)) :: u0u, u0v ! Laplacian sensitivities at u points [L-2 ~> m-2]
  real, dimension(SZI_(G),SZJB_(G)) :: v0u, v0v ! Laplacian sensitivities at v points [L-2 ~> m-2]
  real :: grid_sp_h2       ! Harmonic mean of the squares of the grid [L2 ~> m2]
  real :: grid_sp_h3       ! Harmonic mean of the squares of the grid^(3/2) [L3 ~> m3]
  real :: grid_sp_q2       ! spacings at h and q points [L2 ~> m2]
  real :: grid_sp_q3       ! spacings at h and q points^(3/2) [L3 ~> m3]
  real :: min_grid_sp_h2   ! Minimum value of grid_sp_h2 [L2 ~> m2]
  real :: min_grid_sp_h4   ! Minimum value of grid_sp_h2**2 [L4 ~> m4]
  real :: Kh_Limit         ! A coefficient [T-1 ~> s-1] used, along with the
                           ! grid spacing, to limit Laplacian viscosity.
  real :: fmax             ! maximum absolute value of f at the four
                           ! vorticity points around a thickness point [T-1 ~> s-1]
  real :: BoundCorConst    ! A constant used when using viscosity to bound the Coriolis accelerations
                           ! [T2 L-2 ~> s2 m-2]
  real :: Ah_Limit         ! coefficient [T-1 ~> s-1] used, along with the
                           ! grid spacing, to limit biharmonic viscosity
  real :: Kh               ! Lapacian horizontal viscosity [L2 T-1 ~> m2 s-1]
  real :: Ah               ! biharmonic horizontal viscosity [L4 T-1 ~> m4 s-1]
  real :: Kh_vel_scale     ! this speed [L T-1 ~> m s-1] times grid spacing gives Laplacian viscosity
  real :: Ah_vel_scale     ! this speed [L T-1 ~> m s-1] times grid spacing cubed gives biharmonic viscosity
  real :: Ah_time_scale    ! damping time-scale for biharmonic visc [T ~> s]
  real :: Smag_Lap_const   ! nondimensional Laplacian Smagorinsky constant [nondim]
  real :: Smag_bi_const    ! nondimensional biharmonic Smagorinsky constant [nondim]
  real :: Leith_Lap_const  ! nondimensional Laplacian Leith constant [nondim]
  real :: Leith_bi_const   ! nondimensional biharmonic Leith constant [nondim]
  real :: dt               ! The dynamics time step [T ~> s]
  real :: Idt              ! The inverse of dt [T-1 ~> s-1]
  real :: denom            ! work variable; the denominator of a fraction [L-2 ~> m-2] or [L-4 ~> m-4]
  real :: maxvel           ! largest permitted velocity components [L T-1 ~> m s-1]
  real :: bound_Cor_vel    ! grid-scale velocity variations at which value
                           ! the quadratically varying biharmonic viscosity
                           ! balances Coriolis acceleration [L T-1 ~> m s-1]
  real :: Kh_sin_lat       ! Amplitude of latitudinally dependent viscosity [L2 T-1 ~> m2 s-1]
  real :: Kh_pwr_of_sine   ! Power used to raise sin(lat) when using Kh_sin_lat [nondim]
  logical :: bound_Cor_def ! parameter setting of BOUND_CORIOLIS
  logical :: split         ! If true, use the split time stepping scheme.
                           ! If false and USE_GME = True, issue a FATAL error.
  logical :: use_MEKE      ! If true, the MEKE parameterization is in use.
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags
  character(len=200) :: inputdir, filename ! Input file names and paths
  character(len=80) ::  Kh_var ! Input variable names
  real    :: deg2rad       ! Converts degrees to radians [radians degree-1]
  real    :: slat_fn       ! sin(lat)**Kh_pwr_of_sine [nondim]
  real    :: aniso_grid_dir(2) ! Vector (n1,n2) for anisotropic direction [nondim]
  integer :: aniso_mode    ! Selects the mode for setting the anisotropic direction
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  integer :: i, j
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_hor_visc"  ! module name
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec ; nz = GV%ke
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ! init control structure
  call ZB2020_init(Time, G, GV, US, param_file, diag, CS%ZB2020, CS%use_ZB2020)

  CS%initialized = .true.

  CS%diag => diag
  ! Read parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  ! All parameters are read in all cases to enable parameter spelling checks.
  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
  call get_param(param_file, mdl, "HOR_VISC_ANSWER_DATE", CS%answer_date, &
                 "The vintage of the order of arithmetic and expressions in the horizontal "//&
                 "viscosity calculations.  Values below 20190101 recover the answers from the "//&
                 "end of 2018, while higher values use updated and more robust forms of the "//&
                 "same expressions.", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
  if (.not.GV%Boussinesq) CS%answer_date = max(CS%answer_date, 20230701)

  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false.)
  call get_param(param_file, mdl, "USE_CONT_THICKNESS", CS%use_cont_thick, &
                 "If true, use thickness at velocity points from continuity solver. This option "//&
                 "currently only works with split mode.", default=.false.)
  call get_param(param_file, mdl, "LAPLACIAN", CS%Laplacian, &
                 "If true, use a Laplacian horizontal viscosity.", &
                 default=.false.)

  call get_param(param_file, mdl, "KH", Kh,                      &
                 "The background Laplacian horizontal viscosity.", &
                 units="m2 s-1", default=0.0, scale=US%m_to_L**2*US%T_to_s, &
                 do_not_log=.not.CS%Laplacian)
  call get_param(param_file, mdl, "KH_BG_MIN", CS%Kh_bg_min, &
                 "The minimum value allowed for Laplacian horizontal viscosity, KH.", &
                 units="m2 s-1", default=0.0, scale=US%m_to_L**2*US%T_to_s, &
                 do_not_log=.not.CS%Laplacian)
  call get_param(param_file, mdl, "KH_VEL_SCALE", Kh_vel_scale, &
                 "The velocity scale which is multiplied by the grid "//&
                 "spacing to calculate the Laplacian viscosity. "//&
                 "The final viscosity is the largest of this scaled "//&
                 "viscosity, the Smagorinsky and Leith viscosities, and KH.", &
                 units="m s-1", default=0.0, scale=US%m_s_to_L_T, &
                 do_not_log=.not.CS%Laplacian)
  call get_param(param_file, mdl, "KH_SIN_LAT", Kh_sin_lat, &
                 "The amplitude of a latitudinally-dependent background "//&
                 "viscosity of the form KH_SIN_LAT*(SIN(LAT)**KH_PWR_OF_SINE).", &
                 units="m2 s-1", default=0.0, scale=US%m_to_L**2*US%T_to_s, &
                 do_not_log=.not.CS%Laplacian)
  call get_param(param_file, mdl, "KH_PWR_OF_SINE", Kh_pwr_of_sine, &
                 "The power used to raise SIN(LAT) when using a latitudinally "//&
                 "dependent background viscosity.", &
                 units="nondim", default=4.0, &
                 do_not_log=.not.(CS%Laplacian .and. (Kh_sin_lat>0.)) )
  call get_param(param_file, mdl, "SMAGORINSKY_KH", CS%Smagorinsky_Kh, &
                 "If true, use a Smagorinsky nonlinear eddy viscosity.", &
                 default=.false., do_not_log=.not.CS%Laplacian)
  if (.not.CS%Laplacian) CS%Smagorinsky_Kh = .false.
  call get_param(param_file, mdl, "SMAG_LAP_CONST", Smag_Lap_const, &
                 "The nondimensional Laplacian Smagorinsky constant, "//&
                 "often 0.15.", units="nondim", default=0.0, &
                 fail_if_missing=CS%Smagorinsky_Kh, do_not_log=.not.CS%Smagorinsky_Kh)
  call get_param(param_file, mdl, "LEITH_KH", CS%Leith_Kh, &
                 "If true, use a Leith nonlinear eddy viscosity.", &
                 default=.false., do_not_log=.not.CS%Laplacian)
  if (.not.CS%Laplacian) CS%Leith_Kh = .false.
  call get_param(param_file, mdl, "LEITH_LAP_CONST", Leith_Lap_const, &
                 "The nondimensional Laplacian Leith constant, "//&
                 "often set to 1.0", units="nondim", default=0.0, &
                  fail_if_missing=CS%Leith_Kh, do_not_log=.not.CS%Leith_Kh)
  call get_param(param_file, mdl, "USE_MEKE", use_MEKE, &
                 default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "RES_SCALE_MEKE_VISC", CS%res_scale_MEKE, &
                 "If true, the viscosity contribution from MEKE is scaled by "//&
                 "the resolution function.", default=.false., &
                 do_not_log=.not.(CS%Laplacian.and.use_MEKE))
  if (.not.(CS%Laplacian.and.use_MEKE)) CS%res_scale_MEKE = .false.

  call get_param(param_file, mdl, "BOUND_KH", CS%bound_Kh, &
                 "If true, the Laplacian coefficient is locally limited "//&
                 "to be stable.", default=.true., do_not_log=.not.CS%Laplacian)
  call get_param(param_file, mdl, "BETTER_BOUND_KH", CS%better_bound_Kh, &
                 "If true, the Laplacian coefficient is locally limited "//&
                 "to be stable with a better bounding than just BOUND_KH.", &
                 default=CS%bound_Kh, do_not_log=.not.CS%Laplacian)
  if (.not.CS%Laplacian) CS%bound_Kh = .false.
  if (.not.CS%Laplacian) CS%better_bound_Kh = .false.
  call get_param(param_file, mdl, "ANISOTROPIC_VISCOSITY", CS%anisotropic, &
                 "If true, allow anistropic viscosity in the Laplacian "//&
                 "horizontal viscosity.", default=.false., &
                 do_not_log=.not.CS%Laplacian)
  if (.not.CS%Laplacian) CS%anisotropic = .false. ! This replicates the prior code, but is it intended?
  call get_param(param_file, mdl, "ADD_LES_VISCOSITY", CS%add_LES_viscosity, &
                 "If true, adds the viscosity from Smagorinsky and Leith to the "//&
                 "background viscosity instead of taking the maximum.", default=.false., &
                 do_not_log=.not.CS%Laplacian)

  call get_param(param_file, mdl, "KH_ANISO", CS%Kh_aniso, &
                 "The background Laplacian anisotropic horizontal viscosity.", &
                 units="m2 s-1", default=0.0, scale=US%m_to_L**2*US%T_to_s, &
                 do_not_log=.not.CS%anisotropic)
  call get_param(param_file, mdl, "ANISOTROPIC_MODE", aniso_mode, &
                 "Selects the mode for setting the direction of anisotropy.\n"//&
                 "\t 0 - Points along the grid i-direction.\n"//&
                 "\t 1 - Points towards East.\n"//&
                 "\t 2 - Points along the flow direction, U/|U|.", &
                 default=0, do_not_log=.not.CS%anisotropic)
  if (aniso_mode == 0) then
    call get_param(param_file, mdl, "ANISO_GRID_DIR", aniso_grid_dir, &
                 "The vector pointing in the direction of anisotropy for horizontal viscosity. "//&
                 "n1,n2 are the i,j components relative to the grid.", &
                 units="nondim", fail_if_missing=CS%anisotropic, do_not_log=.not.CS%anisotropic)
  elseif (aniso_mode == 1) then
    call get_param(param_file, mdl, "ANISO_GRID_DIR", aniso_grid_dir, &
                 "The vector pointing in the direction of anisotropy for horizontal viscosity. "//&
                 "n1,n2 are the i,j components relative to the spherical coordinates.", &
                 units="nondim", fail_if_missing=CS%anisotropic, do_not_log=.not.CS%anisotropic)
  else
    call get_param(param_file, mdl, "ANISO_GRID_DIR", aniso_grid_dir, &
                 "The vector pointing in the direction of anisotropy for horizontal viscosity.", &
                 units="nondim", fail_if_missing=.false., do_not_log=.true.)
  endif

  call get_param(param_file, mdl, "BIHARMONIC", CS%biharmonic, &
                 "If true, use a biharmonic horizontal viscosity. "//&
                 "BIHARMONIC may be used with LAPLACIAN.", &
                 default=.true.)
  call get_param(param_file, mdl, "AH", Ah, &
                 "The background biharmonic horizontal viscosity.", &
                 units="m4 s-1", default=0.0, scale=US%m_to_L**4*US%T_to_s, &
                 do_not_log=.not.CS%biharmonic)
  call get_param(param_file, mdl, "AH_VEL_SCALE", Ah_vel_scale, &
                 "The velocity scale which is multiplied by the cube of "//&
                 "the grid spacing to calculate the biharmonic viscosity. "//&
                 "The final viscosity is the largest of this scaled "//&
                 "viscosity, the Smagorinsky and Leith viscosities, and AH.", &
                 units="m s-1", default=0.0, scale=US%m_s_to_L_T, do_not_log=.not.CS%biharmonic)
  call get_param(param_file, mdl, "AH_TIME_SCALE", Ah_time_scale, &
                 "A time scale whose inverse is multiplied by the fourth "//&
                 "power of the grid spacing to calculate biharmonic viscosity. "//&
                 "The final viscosity is the largest of all viscosity "//&
                 "formulations in use. 0.0 means that it's not used.", &
                 units="s", default=0.0, scale=US%s_to_T, do_not_log=.not.CS%biharmonic)
  call get_param(param_file, mdl, "SMAGORINSKY_AH", CS%Smagorinsky_Ah, &
                 "If true, use a biharmonic Smagorinsky nonlinear eddy "//&
                 "viscosity.", default=.false., do_not_log=.not.CS%biharmonic)
  if (.not.CS%biharmonic) CS%Smagorinsky_Ah = .false.
  call get_param(param_file, mdl, "LEITH_AH", CS%Leith_Ah, &
                 "If true, use a biharmonic Leith nonlinear eddy "//&
                 "viscosity.", default=.false., do_not_log=.not.CS%biharmonic)
  if (.not.CS%biharmonic) CS%Leith_Ah = .false.
  call get_param(param_file, mdl, "USE_LEITHY", CS%use_Leithy, &
                 "If true, use a biharmonic Leith nonlinear eddy "//&
                 "viscosity together with a harmonic backscatter.", &
                 default=.false.)
  call get_param(param_file, mdl, "BOUND_AH", CS%bound_Ah, &
                 "If true, the biharmonic coefficient is locally limited "//&
                 "to be stable.", default=.true., do_not_log=.not.CS%biharmonic)
  call get_param(param_file, mdl, "BETTER_BOUND_AH", CS%better_bound_Ah, &
                 "If true, the biharmonic coefficient is locally limited "//&
                 "to be stable with a better bounding than just BOUND_AH.", &
                 default=CS%bound_Ah, do_not_log=.not.CS%biharmonic)
  if (.not.CS%biharmonic) CS%bound_Ah = .false.
  if (.not.CS%biharmonic) CS%better_bound_Ah = .false.
  call get_param(param_file, mdl, "RE_AH", CS%Re_Ah, &
                 "If nonzero, the biharmonic coefficient is scaled "//&
                 "so that the biharmonic Reynolds number is equal to this.", &
                 units="nondim", default=0.0, do_not_log=.not.CS%biharmonic)

  call get_param(param_file, mdl, "BACKSCATTER_UNDERBOUND", CS%backscatter_underbound, &
                 "If true, the bounds on the biharmonic viscosity are allowed to "//&
                 "increase where the Laplacian viscosity is negative (due to backscatter "//&
                 "parameterizations) beyond the largest timestep-dependent stable values of "//&
                 "biharmonic viscosity when no Laplacian viscosity is applied.  The default "//&
                 "is true for historical reasons, but this option probably should not be used "//&
                 "because it can contribute to numerical instabilities.", &
                 default=.true., do_not_log=.not.((CS%better_bound_Kh).and.(CS%better_bound_Ah)))
                 !### The default for BACKSCATTER_UNDERBOUND should be false.

  call get_param(param_file, mdl, "SMAG_BI_CONST",Smag_bi_const, &
                 "The nondimensional biharmonic Smagorinsky constant, "//&
                 "typically 0.015 - 0.06.", units="nondim", default=0.0, &
                 fail_if_missing=CS%Smagorinsky_Ah, do_not_log=.not.CS%Smagorinsky_Ah)

  call get_param(param_file, mdl, "USE_BETA_IN_LEITH", CS%use_beta_in_Leith, &
                 "If true, include the beta term in the Leith nonlinear eddy viscosity.", &
                 default=CS%Leith_Kh, do_not_log=.not.(CS%Leith_Kh .or. CS%Leith_Ah) )
  call get_param(param_file, mdl, "MODIFIED_LEITH", CS%modified_Leith, &
                 "If true, add a term to Leith viscosity which is "//&
                 "proportional to the gradient of divergence.", &
                 default=.false., do_not_log=.not.(CS%Leith_Kh .or. CS%Leith_Ah) )
  call get_param(param_file, mdl, "USE_QG_LEITH_VISC", CS%use_QG_Leith_visc, &
                 "If true, use QG Leith nonlinear eddy viscosity.", &
                 default=.false., do_not_log=.not.(CS%Leith_Kh .or. CS%Leith_Ah) )
!  if (CS%use_QG_Leith_visc) then
!    call MOM_error(FATAL, "USE_QG_LEITH_VISC=True activates code that is a work-in-progress and "//&
!          "should not be used until a number of bugs are fixed.  Specifically it does not "//&
!          "reproduce across PE count or layout, and may use arrays that have not been properly "//&
!          "set or allocated.  See github.com/mom-ocean/MOM6/issues/1590 for a discussion.")
!  endif
  if (CS%use_QG_Leith_visc .and. .not. (CS%Leith_Kh .or. CS%Leith_Ah) ) then
    call MOM_error(FATAL, "MOM_hor_visc.F90, hor_visc_init:"//&
                 "LEITH_KH or LEITH_AH must be True when USE_QG_LEITH_VISC=True.")
  endif

  call get_param(param_file, mdl, "BOUND_CORIOLIS", bound_Cor_def, default=.false.)
  call get_param(param_file, mdl, "BOUND_CORIOLIS_BIHARM", CS%bound_Coriolis, &
                 "If true use a viscosity that increases with the square "//&
                 "of the velocity shears, so that the resulting viscous "//&
                 "drag is of comparable magnitude to the Coriolis terms "//&
                 "when the velocity differences between adjacent grid "//&
                 "points is 0.5*BOUND_CORIOLIS_VEL.  The default is the "//&
                 "value of BOUND_CORIOLIS (or false).", default=bound_Cor_def, &
                 do_not_log=.not.CS%Smagorinsky_Ah)
  if (.not.CS%Smagorinsky_Ah) CS%bound_Coriolis = .false.
  call get_param(param_file, mdl, "MAXVEL", maxvel, &
                 units="m s-1", default=3.0e8, scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "BOUND_CORIOLIS_VEL", bound_Cor_vel, &
                 "The velocity scale at which BOUND_CORIOLIS_BIHARM causes "//&
                 "the biharmonic drag to have comparable magnitude to the "//&
                 "Coriolis acceleration.  The default is set by MAXVEL.", &
                 units="m s-1", default=maxvel*US%L_T_to_m_s, scale=US%m_s_to_L_T, &
                 do_not_log=.not.(CS%Smagorinsky_Ah .and. CS%bound_Coriolis))
  call get_param(param_file, mdl, "LEITH_BI_CONST", Leith_bi_const, &
                 "The nondimensional biharmonic Leith constant, "//&
                 "typical values are thus far undetermined.", units="nondim", default=0.0, &
                 fail_if_missing=(CS%Leith_Ah .or. CS%use_Leithy), &
                 do_not_log=.not.(CS%Leith_Ah .or. CS%use_Leithy))
  call get_param(param_file, mdl, "USE_LAND_MASK_FOR_HVISC", CS%use_land_mask, &
                 "If true, use the land mask for the computation of thicknesses "//&
                 "at velocity locations. This eliminates the dependence on arbitrary "//&
                 "values over land or outside of the domain.", default=.true.)
  call get_param(param_file, mdl, "HORVISC_BOUND_COEF", CS%bound_coef, &
                 "The nondimensional coefficient of the ratio of the "//&
                 "viscosity bounds to the theoretical maximum for "//&
                 "stability without considering other terms.", units="nondim", &
                 default=0.8, do_not_log=.not.(CS%better_bound_Ah .or. CS%better_bound_Kh))
  call get_param(param_file, mdl, "NOSLIP", CS%no_slip, &
                 "If true, no slip boundary conditions are used; otherwise "//&
                 "free slip boundary conditions are assumed. The "//&
                 "implementation of the free slip BCs on a C-grid is much "//&
                 "cleaner than the no slip BCs. The use of free slip BCs "//&
                 "is strongly encouraged, and no slip BCs are not used with "//&
                 "the biharmonic viscosity.", default=.false.)
  call get_param(param_file, mdl, "USE_KH_BG_2D", CS%use_Kh_bg_2d, &
                 "If true, read a file containing 2-d background harmonic "//&
                 "viscosities. The final viscosity is the maximum of the other "//&
                 "terms and this background value.", default=.false., do_not_log=.not.CS%Laplacian)
  if (.not.CS%Laplacian) CS%use_Kh_bg_2d = .false.
  call get_param(param_file, mdl, "KH_BG_2D_BUG", CS%Kh_bg_2d_bug, &
                 "If true, retain an answer-changing horizontal indexing bug in setting "//&
                 "the corner-point viscosities when USE_KH_BG_2D=True.  This is "//&
                 "not recommended.", default=.false., do_not_log=.not.CS%use_Kh_bg_2d)
  call get_param(param_file, mdl, "FRICTWORK_BUG", CS%FrictWork_bug, &
                 "If true, retain an answer-changing bug in calculating "//&
                 "the FrictWork, which cancels the h in thickness flux and the h at velocity point. This is"//&
                 "not recommended.", default=.true.)

  call get_param(param_file, mdl, "USE_GME", CS%use_GME, &
                 "If true, use the GM+E backscatter scheme in association \n"//&
                 "with the Gent and McWilliams parameterization.", default=.false.)
  call get_param(param_file, mdl, "SPLIT", split, &
                 "Use the split time stepping if true.", default=.true., do_not_log=.true.)
  if (CS%use_Leithy) then
    if (.not.(CS%biharmonic .and. CS%Laplacian)) then
                   call MOM_error(FATAL, "MOM_hor_visc.F90, hor_visc_init: "//&
                   "LAPLACIAN and BIHARMONIC must both be True when USE_LEITHY=True.")
    endif
  endif
  call get_param(param_file, mdl, "LEITHY_CK", CS%c_K, &
                 "Fraction of biharmonic dissipation that gets backscattered, "//&
                 "in Leith+E.", units="nondim", default=1.0, do_not_log=.not.CS%use_Leithy)
  call get_param(param_file, mdl, "SMOOTH_AH", CS%smooth_Ah, &
                 "If true, Ah and m_leithy are smoothed within Leith+E.  This requires "//&
                 "lots of blocking communications, which can be expensive", &
                 default=.true., do_not_log=.not.CS%use_Leithy)

  if (CS%use_GME .and. .not.split) call MOM_error(FATAL,"ERROR: Currently, USE_GME = True "// &
                                           "cannot be used with SPLIT=False.")

  if (CS%use_GME) then
    call get_param(param_file, mdl, "GME_NUM_SMOOTHINGS", CS%num_smooth_gme, &
                   "Number of smoothing passes for the GME fluxes.", &
                   default=1)
    call get_param(param_file, mdl, "GME_H0", CS%GME_h0, &
                   "The strength of GME tapers quadratically to zero when the bathymetric "//&
                   "depth is shallower than GME_H0.", &
                   units="m", scale=GV%m_to_H, default=1000.0)
    call get_param(param_file, mdl, "GME_EFFICIENCY", CS%GME_efficiency, &
                   "The nondimensional prefactor multiplying the GME coefficient.", &
                   units="nondim", default=1.0)
    call get_param(param_file, mdl, "GME_LIMITER", CS%GME_limiter, &
                   "The absolute maximum value the GME coefficient is allowed to take.", &
                   units="m2 s-1", scale=US%m_to_L**2*US%T_to_s, default=1.0e7)
  endif

  if (CS%Laplacian .or. CS%biharmonic) then
    call get_param(param_file, mdl, "DT", dt, &
                 "The (baroclinic) dynamics time step.", units="s", scale=US%s_to_T, &
                 fail_if_missing=.true.)
    Idt = 1.0 / dt
  endif
  if (CS%no_slip .and. CS%biharmonic) &
    call MOM_error(FATAL,"ERROR: NOSLIP and BIHARMONIC cannot be defined "// &
                         "at the same time in MOM.")
  if (.not.(CS%Laplacian .or. CS%biharmonic)) then
    ! Only issue inviscid warning if not in single column mode (usually 2x2 domain)
    if ( max(G%domain%niglobal, G%domain%njglobal)>2 ) call MOM_error(WARNING, &
      "hor_visc_init:  It is usually a very bad idea not to use either "//&
      "LAPLACIAN or BIHARMONIC viscosity.")
    return ! We are not using either Laplacian or Bi-harmonic lateral viscosity
  endif
  deg2rad = atan(1.0) / 45.
  ALLOC_(CS%dx2h(isd:ied,jsd:jed))        ; CS%dx2h(:,:)    = 0.0
  ALLOC_(CS%dy2h(isd:ied,jsd:jed))        ; CS%dy2h(:,:)    = 0.0
  ALLOC_(CS%dx2q(IsdB:IedB,JsdB:JedB))    ; CS%dx2q(:,:)    = 0.0
  ALLOC_(CS%dy2q(IsdB:IedB,JsdB:JedB))    ; CS%dy2q(:,:)    = 0.0
  ALLOC_(CS%dx_dyT(isd:ied,jsd:jed))      ; CS%dx_dyT(:,:)  = 0.0
  ALLOC_(CS%dy_dxT(isd:ied,jsd:jed))      ; CS%dy_dxT(:,:)  = 0.0
  ALLOC_(CS%dx_dyBu(IsdB:IedB,JsdB:JedB)) ; CS%dx_dyBu(:,:) = 0.0
  ALLOC_(CS%dy_dxBu(IsdB:IedB,JsdB:JedB)) ; CS%dy_dxBu(:,:) = 0.0
  if (CS%Laplacian) then
    ALLOC_(CS%grid_sp_h2(isd:ied,jsd:jed))   ; CS%grid_sp_h2(:,:) = 0.0
    ALLOC_(CS%Kh_bg_xx(isd:ied,jsd:jed))     ; CS%Kh_bg_xx(:,:) = 0.0
    ALLOC_(CS%Kh_bg_xy(IsdB:IedB,JsdB:JedB)) ; CS%Kh_bg_xy(:,:) = 0.0
    if (CS%bound_Kh .or. CS%better_bound_Kh) then
      ALLOC_(CS%Kh_Max_xx(Isd:Ied,Jsd:Jed)) ; CS%Kh_Max_xx(:,:) = 0.0
      ALLOC_(CS%Kh_Max_xy(IsdB:IedB,JsdB:JedB)) ; CS%Kh_Max_xy(:,:) = 0.0
    endif
    if (CS%Smagorinsky_Kh) then
      ALLOC_(CS%Laplac2_const_xx(isd:ied,jsd:jed))     ; CS%Laplac2_const_xx(:,:) = 0.0
      ALLOC_(CS%Laplac2_const_xy(IsdB:IedB,JsdB:JedB)) ; CS%Laplac2_const_xy(:,:) = 0.0
    endif
    if (CS%Leith_Kh) then
      ALLOC_(CS%Laplac3_const_xx(isd:ied,jsd:jed)) ; CS%Laplac3_const_xx(:,:) = 0.0
      ALLOC_(CS%Laplac3_const_xy(IsdB:IedB,JsdB:JedB)) ; CS%Laplac3_const_xy(:,:) = 0.0
    endif
  endif
  ALLOC_(CS%reduction_xx(isd:ied,jsd:jed))     ; CS%reduction_xx(:,:) = 0.0
  ALLOC_(CS%reduction_xy(IsdB:IedB,JsdB:JedB)) ; CS%reduction_xy(:,:) = 0.0

  CS%dynamic_aniso = .false.
  if (CS%anisotropic) then
    ALLOC_(CS%n1n2_h(isd:ied,jsd:jed)) ; CS%n1n2_h(:,:) = 0.0
    ALLOC_(CS%n1n1_m_n2n2_h(isd:ied,jsd:jed)) ; CS%n1n1_m_n2n2_h(:,:) = 0.0
    ALLOC_(CS%n1n2_q(IsdB:IedB,JsdB:JedB)) ; CS%n1n2_q(:,:) = 0.0
    ALLOC_(CS%n1n1_m_n2n2_q(IsdB:IedB,JsdB:JedB)) ; CS%n1n1_m_n2n2_q(:,:) = 0.0
    select case (aniso_mode)
      case (0)
        call align_aniso_tensor_to_grid(CS, aniso_grid_dir(1), aniso_grid_dir(2))
      case (1)
      ! call align_aniso_tensor_to_grid(CS, aniso_grid_dir(1), aniso_grid_dir(2))
      case (2)
        CS%dynamic_aniso = .true.
      case default
        call MOM_error(FATAL, "MOM_hor_visc: "//&
             "Runtime parameter ANISOTROPIC_MODE is out of range.")
    end select
  endif

  call get_param(param_file, mdl, "KH_BG_2D_FILENAME", filename, &
                 'The filename containing a 2d map of "Kh".', &
                 default='KH_background_2d.nc', do_not_log=.not.CS%use_Kh_bg_2d)
  call get_param(param_file, mdl, "KH_BG_2D_VARNAME", Kh_var, &
                 'The name in the input file of the horizontal viscosity variable.', &
                 default='Kh', do_not_log=.not.CS%use_Kh_bg_2d)

  if (CS%use_Kh_bg_2d) then
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    ALLOC_(CS%Kh_bg_2d(isd:ied,jsd:jed))     ; CS%Kh_bg_2d(:,:) = 0.0
    call MOM_read_data(trim(inputdir)//trim(filename), Kh_var, CS%Kh_bg_2d, &
                       G%domain, timelevel=1, scale=US%m_to_L**2*US%T_to_s)
    call pass_var(CS%Kh_bg_2d, G%domain)
  endif
  if (CS%biharmonic) then
    ALLOC_(CS%Idx2dyCu(IsdB:IedB,jsd:jed)) ; CS%Idx2dyCu(:,:) = 0.0
    ALLOC_(CS%Idx2dyCv(isd:ied,JsdB:JedB)) ; CS%Idx2dyCv(:,:) = 0.0
    ALLOC_(CS%Idxdy2u(IsdB:IedB,jsd:jed))  ; CS%Idxdy2u(:,:)  = 0.0
    ALLOC_(CS%Idxdy2v(isd:ied,JsdB:JedB))  ; CS%Idxdy2v(:,:)  = 0.0
    ALLOC_(CS%Ah_bg_xx(isd:ied,jsd:jed))     ; CS%Ah_bg_xx(:,:) = 0.0
    ALLOC_(CS%Ah_bg_xy(IsdB:IedB,JsdB:JedB)) ; CS%Ah_bg_xy(:,:) = 0.0
    ALLOC_(CS%grid_sp_h3(isd:ied,jsd:jed))   ; CS%grid_sp_h3(:,:) = 0.0
    if (CS%bound_Ah .or. CS%better_bound_Ah) then
      ALLOC_(CS%Ah_Max_xx(isd:ied,jsd:jed))     ; CS%Ah_Max_xx(:,:) = 0.0
      ALLOC_(CS%Ah_Max_xy(IsdB:IedB,JsdB:JedB)) ; CS%Ah_Max_xy(:,:) = 0.0
    endif
    if (CS%Smagorinsky_Ah) then
      ALLOC_(CS%Biharm_const_xx(isd:ied,jsd:jed))     ; CS%Biharm_const_xx(:,:) = 0.0
      ALLOC_(CS%Biharm_const_xy(IsdB:IedB,JsdB:JedB)) ; CS%Biharm_const_xy(:,:) = 0.0
      if (CS%bound_Coriolis) then
        ALLOC_(CS%Biharm_const2_xx(isd:ied,jsd:jed))     ; CS%Biharm_const2_xx(:,:) = 0.0
        ALLOC_(CS%Biharm_const2_xy(IsdB:IedB,JsdB:JedB)) ; CS%Biharm_const2_xy(:,:) = 0.0
      endif
    endif
    if ((CS%Leith_Ah) .or. (CS%use_Leithy)) then
      ALLOC_(CS%biharm6_const_xx(isd:ied,jsd:jed)) ; CS%biharm6_const_xx(:,:) = 0.0
      ALLOC_(CS%biharm6_const_xy(IsdB:IedB,JsdB:JedB)) ; CS%biharm6_const_xy(:,:) = 0.0
    endif
    if (CS%use_Leithy) then
      ALLOC_(CS%m_const_leithy(isd:ied,jsd:jed)) ; CS%m_const_leithy(:,:) = 0.0
      ALLOC_(CS%m_leithy_max(isd:ied,jsd:jed)) ; CS%m_leithy_max(:,:) = 0.0
    endif
    if (CS%Re_Ah > 0.0) then
      ALLOC_(CS%Re_Ah_const_xx(isd:ied,jsd:jed)); CS%Re_Ah_const_xx(:,:) = 0.0
      ALLOC_(CS%Re_Ah_const_xy(IsdB:IedB,JsdB:JedB)); CS%Re_Ah_const_xy(:,:) = 0.0
    endif
  endif
  do J=js-2,Jeq+1 ; do I=is-2,Ieq+1
    CS%dx2q(I,J) = G%dxBu(I,J)*G%dxBu(I,J) ; CS%dy2q(I,J) = G%dyBu(I,J)*G%dyBu(I,J)
    CS%DX_dyBu(I,J) = G%dxBu(I,J)*G%IdyBu(I,J) ; CS%DY_dxBu(I,J) = G%dyBu(I,J)*G%IdxBu(I,J)
  enddo ; enddo
  do j=js-2,Jeq+2 ; do i=is-2,Ieq+2
    CS%dx2h(i,j) = G%dxT(i,j)*G%dxT(i,j) ; CS%dy2h(i,j) = G%dyT(i,j)*G%dyT(i,j)
    CS%DX_dyT(i,j) = G%dxT(i,j)*G%IdyT(i,j) ; CS%DY_dxT(i,j) = G%dyT(i,j)*G%IdxT(i,j)
  enddo ; enddo
  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    CS%reduction_xx(i,j) = 1.0
    if ((G%dy_Cu(I,j) > 0.0) .and. (G%dy_Cu(I,j) < G%dyCu(I,j)) .and. &
        (G%dy_Cu(I,j) < G%dyCu(I,j) * CS%reduction_xx(i,j))) &
      CS%reduction_xx(i,j) = G%dy_Cu(I,j) / (G%dyCu(I,j))
    if ((G%dy_Cu(I-1,j) > 0.0) .and. (G%dy_Cu(I-1,j) < G%dyCu(I-1,j)) .and. &
        (G%dy_Cu(I-1,j) < G%dyCu(I-1,j) * CS%reduction_xx(i,j))) &
      CS%reduction_xx(i,j) = G%dy_Cu(I-1,j) / (G%dyCu(I-1,j))
    if ((G%dx_Cv(i,J) > 0.0) .and. (G%dx_Cv(i,J) < G%dxCv(i,J)) .and. &
        (G%dx_Cv(i,J) < G%dxCv(i,J) * CS%reduction_xx(i,j))) &
      CS%reduction_xx(i,j) = G%dx_Cv(i,J) / (G%dxCv(i,J))
    if ((G%dx_Cv(i,J-1) > 0.0) .and. (G%dx_Cv(i,J-1) < G%dxCv(i,J-1)) .and. &
        (G%dx_Cv(i,J-1) < G%dxCv(i,J-1) * CS%reduction_xx(i,j))) &
      CS%reduction_xx(i,j) = G%dx_Cv(i,J-1) / (G%dxCv(i,J-1))
  enddo ; enddo
  do J=js-1,Jeq ; do I=is-1,Ieq
    CS%reduction_xy(I,J) = 1.0
    if ((G%dy_Cu(I,j) > 0.0) .and. (G%dy_Cu(I,j) < G%dyCu(I,j)) .and. &
        (G%dy_Cu(I,j) < G%dyCu(I,j) * CS%reduction_xy(I,J))) &
      CS%reduction_xy(I,J) = G%dy_Cu(I,j) / (G%dyCu(I,j))
    if ((G%dy_Cu(I,j+1) > 0.0) .and. (G%dy_Cu(I,j+1) < G%dyCu(I,j+1)) .and. &
        (G%dy_Cu(I,j+1) < G%dyCu(I,j+1) * CS%reduction_xy(I,J))) &
      CS%reduction_xy(I,J) = G%dy_Cu(I,j+1) / (G%dyCu(I,j+1))
    if ((G%dx_Cv(i,J) > 0.0) .and. (G%dx_Cv(i,J) < G%dxCv(i,J)) .and. &
        (G%dx_Cv(i,J) < G%dxCv(i,J) * CS%reduction_xy(I,J))) &
      CS%reduction_xy(I,J) = G%dx_Cv(i,J) / (G%dxCv(i,J))
    if ((G%dx_Cv(i+1,J) > 0.0) .and. (G%dx_Cv(i+1,J) < G%dxCv(i+1,J)) .and. &
        (G%dx_Cv(i+1,J) < G%dxCv(i+1,J) * CS%reduction_xy(I,J))) &
      CS%reduction_xy(I,J) = G%dx_Cv(i+1,J) / (G%dxCv(i+1,J))
  enddo ; enddo
  if (CS%Laplacian) then
   ! The 0.3 below was 0.4 in MOM1.10.  The change in hq requires
   ! this to be less than 1/3, rather than 1/2 as before.
    if (CS%bound_Kh .or. CS%bound_Ah) Kh_Limit = 0.3 / (dt*4.0)
    ! Calculate and store the background viscosity at h-points

    min_grid_sp_h2 = huge(1.)
    do j=js-1,Jeq+1 ; do i=is-1,Ieq+1
      ! Static factors in the Smagorinsky and Leith schemes
      grid_sp_h2 = (2.0*CS%dx2h(i,j)*CS%dy2h(i,j)) / (CS%dx2h(i,j) + CS%dy2h(i,j))
      CS%grid_sp_h2(i,j) = grid_sp_h2
      grid_sp_h3 = grid_sp_h2*sqrt(grid_sp_h2)
      if (CS%Smagorinsky_Kh) CS%Laplac2_const_xx(i,j) = Smag_Lap_const * grid_sp_h2
      if (CS%Leith_Kh)       CS%Laplac3_const_xx(i,j) = Leith_Lap_const * grid_sp_h3
      ! Maximum of constant background and MICOM viscosity
      CS%Kh_bg_xx(i,j) = MAX(Kh, Kh_vel_scale * sqrt(grid_sp_h2))
      ! Use the larger of the above and values read from a file
      if (CS%use_Kh_bg_2d) CS%Kh_bg_xx(i,j) = MAX(CS%Kh_bg_2d(i,j), CS%Kh_bg_xx(i,j))
      ! Use the larger of the above and a function of sin(latitude)
      if (Kh_sin_lat>0.) then
        slat_fn = abs( sin( deg2rad * G%geoLatT(i,j) ) ) ** Kh_pwr_of_sine
        CS%Kh_bg_xx(i,j) = MAX(Kh_sin_lat * slat_fn, CS%Kh_bg_xx(i,j))
      endif
      if (CS%bound_Kh .and. .not.CS%better_bound_Kh) then
        ! Limit the background viscosity to be numerically stable
        CS%Kh_Max_xx(i,j) = Kh_Limit * grid_sp_h2
        CS%Kh_bg_xx(i,j) = MIN(CS%Kh_bg_xx(i,j), CS%Kh_Max_xx(i,j))
      endif
      min_grid_sp_h2 = min(grid_sp_h2, min_grid_sp_h2)
    enddo ; enddo
    call min_across_PEs(min_grid_sp_h2)

    ! Calculate and store the background viscosity at q-points
    do J=js-1,Jeq ; do I=is-1,Ieq
      ! Static factors in the Smagorinsky and Leith schemes
      grid_sp_q2 = (2.0*CS%dx2q(I,J)*CS%dy2q(I,J)) / (CS%dx2q(I,J) + CS%dy2q(I,J))
      grid_sp_q3 = grid_sp_q2*sqrt(grid_sp_q2)
      if (CS%Smagorinsky_Kh) CS%Laplac2_const_xy(I,J) = Smag_Lap_const * grid_sp_q2
      if (CS%Leith_Kh)       CS%Laplac3_const_xy(I,J) = Leith_Lap_const * grid_sp_q3
      ! Maximum of constant background and MICOM viscosity
      CS%Kh_bg_xy(I,J) = MAX(Kh, Kh_vel_scale * sqrt(grid_sp_q2))
      ! Use the larger of the above and values read from a file
      if (CS%use_Kh_bg_2d) then
        if (CS%Kh_bg_2d_bug) then
          ! This option is unambiguously wrong but is needed to recover old answers
          CS%Kh_bg_xy(I,J) = MAX(CS%Kh_bg_2d(i,j), CS%Kh_bg_xy(I,J))
        else
          CS%Kh_bg_xy(I,J) = MAX(CS%Kh_bg_xy(I,J), &
              0.25*((CS%Kh_bg_2d(i,j) + CS%Kh_bg_2d(i+1,j+1)) + &
                    (CS%Kh_bg_2d(i+1,j) + CS%Kh_bg_2d(i,j+1))) )
        endif
      endif

      ! Use the larger of the above and a function of sin(latitude)
      if (Kh_sin_lat>0.) then
        slat_fn = abs( sin( deg2rad * G%geoLatBu(I,J) ) ) ** Kh_pwr_of_sine
        CS%Kh_bg_xy(I,J) = MAX(Kh_sin_lat * slat_fn, CS%Kh_bg_xy(I,J))
      endif
      if (CS%bound_Kh .and. .not.CS%better_bound_Kh) then
        ! Limit the background viscosity to be numerically stable
        CS%Kh_Max_xy(I,J) = Kh_Limit * grid_sp_q2
        CS%Kh_bg_xy(I,J) = MIN(CS%Kh_bg_xy(I,J), CS%Kh_Max_xy(I,J))
      endif
    enddo ; enddo
  endif
  if (CS%biharmonic) then
    do j=js-1,Jeq+1 ; do I=is-2,Ieq+1
      CS%Idx2dyCu(I,j) = (G%IdxCu(I,j)*G%IdxCu(I,j)) * G%IdyCu(I,j)
      CS%Idxdy2u(I,j) = G%IdxCu(I,j) * (G%IdyCu(I,j)*G%IdyCu(I,j))
    enddo ; enddo
    do J=js-2,Jeq+1 ; do i=is-1,Ieq+1
      CS%Idx2dyCv(i,J) = (G%IdxCv(i,J)*G%IdxCv(i,J)) * G%IdyCv(i,J)
      CS%Idxdy2v(i,J) = G%IdxCv(i,J) * (G%IdyCv(i,J)*G%IdyCv(i,J))
    enddo ; enddo
    CS%Ah_bg_xy(:,:) = 0.0
    ! The 0.3 below was 0.4 in HIM 1.10.  The change in hq requires
    ! this to be less than 1/3, rather than 1/2 as before.
    if (CS%better_bound_Ah .or. CS%bound_Ah) Ah_Limit = 0.3 / (dt*64.0)
    if (CS%Smagorinsky_Ah .and. CS%bound_Coriolis) &
      BoundCorConst = 1.0 / (5.0*(bound_Cor_vel*bound_Cor_vel))

    min_grid_sp_h4 = huge(1.)
    do j=js-1,Jeq+1 ; do i=is-1,Ieq+1
      grid_sp_h2 = (2.0*CS%dx2h(i,j)*CS%dy2h(i,j)) / (CS%dx2h(i,j)+CS%dy2h(i,j))
      grid_sp_h3 = grid_sp_h2*sqrt(grid_sp_h2)
      CS%grid_sp_h3(i,j) = grid_sp_h3
      if (CS%Smagorinsky_Ah) then
        CS%Biharm_const_xx(i,j) = Smag_bi_const * (grid_sp_h2 * grid_sp_h2)
        if (CS%bound_Coriolis) then
          fmax = MAX(abs(G%CoriolisBu(I-1,J-1)), abs(G%CoriolisBu(I,J-1)), &
                     abs(G%CoriolisBu(I-1,J)),   abs(G%CoriolisBu(I,J)))
          CS%Biharm_const2_xx(i,j) = (grid_sp_h2 * grid_sp_h2 * grid_sp_h2) * &
                                     (fmax * BoundCorConst)
        endif
      endif
      if (CS%Leith_Ah) then
         CS%biharm6_const_xx(i,j) = Leith_bi_const * (grid_sp_h3 * grid_sp_h3)
      endif
      if (CS%use_Leithy) then
         CS%biharm6_const_xx(i,j) = Leith_bi_const * max(G%dxT(i,j),G%dyT(i,j))**6
         CS%m_const_leithy(i,j) = 0.5 * sqrt(CS%c_K) * max(G%dxT(i,j),G%dyT(i,j))
         CS%m_leithy_max(i,j) = 4. / max(G%dxT(i,j),G%dyT(i,j))**2
      endif
      CS%Ah_bg_xx(i,j) = MAX(Ah, Ah_vel_scale * grid_sp_h2 * sqrt(grid_sp_h2))
      if (CS%Re_Ah > 0.0) CS%Re_Ah_const_xx(i,j) = grid_sp_h3 / CS%Re_Ah
      if (Ah_time_scale > 0.) CS%Ah_bg_xx(i,j) = &
            MAX(CS%Ah_bg_xx(i,j), (grid_sp_h2 * grid_sp_h2) / Ah_time_scale)
      if (CS%bound_Ah .and. .not.CS%better_bound_Ah) then
        CS%Ah_Max_xx(i,j) = Ah_Limit * (grid_sp_h2 * grid_sp_h2)
        CS%Ah_bg_xx(i,j) = MIN(CS%Ah_bg_xx(i,j), CS%Ah_Max_xx(i,j))
      endif
      min_grid_sp_h4 = min(grid_sp_h2**2, min_grid_sp_h4)
    enddo ; enddo
    call min_across_PEs(min_grid_sp_h4)

    do J=js-1,Jeq ; do I=is-1,Ieq
      grid_sp_q2 = (2.0*CS%dx2q(I,J)*CS%dy2q(I,J)) / (CS%dx2q(I,J)+CS%dy2q(I,J))
      grid_sp_q3 = grid_sp_q2*sqrt(grid_sp_q2)
      if (CS%Smagorinsky_Ah) then
        CS%Biharm_const_xy(I,J) = Smag_bi_const * (grid_sp_q2 * grid_sp_q2)
        if (CS%bound_Coriolis) then
          CS%Biharm_const2_xy(I,J) = (grid_sp_q2 * grid_sp_q2 * grid_sp_q2) * &
                                     (abs(G%CoriolisBu(I,J)) * BoundCorConst)
        endif
      endif
      if ((CS%Leith_Ah) .or. (CS%use_Leithy))then
         CS%biharm6_const_xy(I,J) = Leith_bi_const * (grid_sp_q3 * grid_sp_q3)
      endif
      CS%Ah_bg_xy(I,J) = MAX(Ah, Ah_vel_scale * grid_sp_q2 * sqrt(grid_sp_q2))
      if (CS%Re_Ah > 0.0) CS%Re_Ah_const_xy(i,j) = grid_sp_q3 / CS%Re_Ah
      if (Ah_time_scale > 0.) CS%Ah_bg_xy(i,j) = &
           MAX(CS%Ah_bg_xy(i,j), (grid_sp_q2 * grid_sp_q2) / Ah_time_scale)
      if (CS%bound_Ah .and. .not.CS%better_bound_Ah) then
        CS%Ah_Max_xy(I,J) = Ah_Limit * (grid_sp_q2 * grid_sp_q2)
        CS%Ah_bg_xy(I,J) = MIN(CS%Ah_bg_xy(I,J), CS%Ah_Max_xy(I,J))
      endif
    enddo ; enddo
  endif
  ! The Laplacian bounds should avoid overshoots when CS%bound_coef < 1.
  if (CS%Laplacian .and. CS%better_bound_Kh) then
    do j=js-1,Jeq+1 ; do i=is-1,Ieq+1
      denom = max( &
         (CS%dy2h(i,j) * CS%DY_dxT(i,j) * (G%IdyCu(I,j) + G%IdyCu(I-1,j)) * &
          max(G%IdyCu(I,j)*G%IareaCu(I,j), G%IdyCu(I-1,j)*G%IareaCu(I-1,j)) ), &
         (CS%dx2h(i,j) * CS%DX_dyT(i,j) * (G%IdxCv(i,J) + G%IdxCv(i,J-1)) * &
          max(G%IdxCv(i,J)*G%IareaCv(i,J), G%IdxCv(i,J-1)*G%IareaCv(i,J-1)) ) )
      CS%Kh_Max_xx(i,j) = 0.0
      if (denom > 0.0) &
        CS%Kh_Max_xx(i,j) = CS%bound_coef * 0.25 * Idt / denom
    enddo ; enddo
    do J=js-1,Jeq ; do I=is-1,Ieq
      denom = max( &
         (CS%dx2q(I,J) * CS%DX_dyBu(I,J) * (G%IdxCu(I,j+1) + G%IdxCu(I,j)) * &
          max(G%IdxCu(I,j)*G%IareaCu(I,j), G%IdxCu(I,j+1)*G%IareaCu(I,j+1)) ), &
         (CS%dy2q(I,J) * CS%DY_dxBu(I,J) * (G%IdyCv(i+1,J) + G%IdyCv(i,J)) * &
          max(G%IdyCv(i,J)*G%IareaCv(i,J), G%IdyCv(i+1,J)*G%IareaCv(i+1,J)) ) )
      CS%Kh_Max_xy(I,J) = 0.0
      if (denom > 0.0) &
        CS%Kh_Max_xy(I,J) = CS%bound_coef * 0.25 * Idt / denom
    enddo ; enddo
    if (CS%debug) then
      call hchksum(CS%Kh_Max_xx, "Kh_Max_xx", G%HI, haloshift=0, unscale=US%L_to_m**2*US%s_to_T)
      call Bchksum(CS%Kh_Max_xy, "Kh_Max_xy", G%HI, haloshift=0, unscale=US%L_to_m**2*US%s_to_T)
    endif
  endif
  ! The biharmonic bounds should avoid overshoots when CS%bound_coef < 0.5, but
  ! empirically work for CS%bound_coef <~ 1.0
  if (CS%biharmonic .and. CS%better_bound_Ah) then
    do j=js-1,Jeq+1 ; do I=is-2,Ieq+1
      u0u(I,j) = (CS%Idxdy2u(I,j)*(CS%dy2h(i+1,j)*CS%DY_dxT(i+1,j)*(G%IdyCu(I+1,j) + G%IdyCu(I,j))   + &
                                   CS%dy2h(i,j) * CS%DY_dxT(i,j) * (G%IdyCu(I,j) + G%IdyCu(I-1,j)) ) + &
                 CS%Idx2dyCu(I,j)*(CS%dx2q(I,J) * CS%DX_dyBu(I,J) * (G%IdxCu(I,j+1) + G%IdxCu(I,j)) + &
                                   CS%dx2q(I,J-1)*CS%DX_dyBu(I,J-1)*(G%IdxCu(I,j) + G%IdxCu(I,j-1)) ) )
      u0v(I,j) = (CS%Idxdy2u(I,j)*(CS%dy2h(i+1,j)*CS%DX_dyT(i+1,j)*(G%IdxCv(i+1,J) + G%IdxCv(i+1,J-1)) + &
                                   CS%dy2h(i,j) * CS%DX_dyT(i,j) * (G%IdxCv(i,J) + G%IdxCv(i,J-1)) )   + &
                 CS%Idx2dyCu(I,j)*(CS%dx2q(I,J) * CS%DY_dxBu(I,J) * (G%IdyCv(i+1,J) + G%IdyCv(i,J))   + &
                                   CS%dx2q(I,J-1)*CS%DY_dxBu(I,J-1)*(G%IdyCv(i+1,J-1) + G%IdyCv(i,J-1)) ) )
    enddo ; enddo
    do J=js-2,Jeq+1 ; do i=is-1,Ieq+1
      v0u(i,J) = (CS%Idxdy2v(i,J)*(CS%dy2q(I,J) * CS%DX_dyBu(I,J) * (G%IdxCu(I,j+1) + G%IdxCu(I,j))       + &
                                   CS%dy2q(I-1,J)*CS%DX_dyBu(I-1,J)*(G%IdxCu(I-1,j+1) + G%IdxCu(I-1,j)) ) + &
                 CS%Idx2dyCv(i,J)*(CS%dx2h(i,j+1)*CS%DY_dxT(i,j+1)*(G%IdyCu(I,j+1) + G%IdyCu(I-1,j+1))   + &
                                   CS%dx2h(i,j) * CS%DY_dxT(i,j) * (G%IdyCu(I,j) + G%IdyCu(I-1,j)) ) )
      v0v(i,J) = (CS%Idxdy2v(i,J)*(CS%dy2q(I,J) * CS%DY_dxBu(I,J) * (G%IdyCv(i+1,J) + G%IdyCv(i,J))   + &
                                   CS%dy2q(I-1,J)*CS%DY_dxBu(I-1,J)*(G%IdyCv(i,J) + G%IdyCv(i-1,J)) ) + &
                 CS%Idx2dyCv(i,J)*(CS%dx2h(i,j+1)*CS%DX_dyT(i,j+1)*(G%IdxCv(i,J+1) + G%IdxCv(i,J))   + &
                                   CS%dx2h(i,j) * CS%DX_dyT(i,j) * (G%IdxCv(i,J) + G%IdxCv(i,J-1)) ) )
    enddo ; enddo
    do j=js-1,Jeq+1 ; do i=is-1,Ieq+1
      denom = max( &
         (CS%dy2h(i,j) * &
          (CS%DY_dxT(i,j)*(G%IdyCu(I,j)*u0u(I,j) + G%IdyCu(I-1,j)*u0u(I-1,j))  + &
           CS%DX_dyT(i,j)*(G%IdxCv(i,J)*v0u(i,J) + G%IdxCv(i,J-1)*v0u(i,J-1))) * &
          max(G%IdyCu(I,j)*G%IareaCu(I,j), G%IdyCu(I-1,j)*G%IareaCu(I-1,j)) ),   &
         (CS%dx2h(i,j) * &
          (CS%DY_dxT(i,j)*(G%IdyCu(I,j)*u0v(I,j) + G%IdyCu(I-1,j)*u0v(I-1,j))  + &
           CS%DX_dyT(i,j)*(G%IdxCv(i,J)*v0v(i,J) + G%IdxCv(i,J-1)*v0v(i,J-1))) * &
          max(G%IdxCv(i,J)*G%IareaCv(i,J), G%IdxCv(i,J-1)*G%IareaCv(i,J-1)) ) )
      CS%Ah_Max_xx(I,J) = 0.0
      if (denom > 0.0) &
        CS%Ah_Max_xx(I,J) = CS%bound_coef * 0.5 * Idt / denom
    enddo ; enddo
    do J=js-1,Jeq ; do I=is-1,Ieq
      denom = max( &
         (CS%dx2q(I,J) * &
          (CS%DX_dyBu(I,J)*(u0u(I,j+1)*G%IdxCu(I,j+1) + u0u(I,j)*G%IdxCu(I,j))  + &
           CS%DY_dxBu(I,J)*(v0u(i+1,J)*G%IdyCv(i+1,J) + v0u(i,J)*G%IdyCv(i,J))) * &
          max(G%IdxCu(I,j)*G%IareaCu(I,j), G%IdxCu(I,j+1)*G%IareaCu(I,j+1)) ),    &
         (CS%dy2q(I,J) * &
          (CS%DX_dyBu(I,J)*(u0v(I,j+1)*G%IdxCu(I,j+1) + u0v(I,j)*G%IdxCu(I,j))  + &
           CS%DY_dxBu(I,J)*(v0v(i+1,J)*G%IdyCv(i+1,J) + v0v(i,J)*G%IdyCv(i,J))) * &
          max(G%IdyCv(i,J)*G%IareaCv(i,J), G%IdyCv(i+1,J)*G%IareaCv(i+1,J)) ) )
      CS%Ah_Max_xy(I,J) = 0.0
      if (denom > 0.0) &
        CS%Ah_Max_xy(I,J) = CS%bound_coef * 0.5 * Idt / denom
    enddo ; enddo
    if (CS%debug) then
      call hchksum(CS%Ah_Max_xx, "Ah_Max_xx", G%HI, haloshift=0, unscale=US%L_to_m**4*US%s_to_T)
      call Bchksum(CS%Ah_Max_xy, "Ah_Max_xy", G%HI, haloshift=0, unscale=US%L_to_m**4*US%s_to_T)
    endif
  endif
  ! Register fields for output from this module.
  CS%id_normstress = register_diag_field('ocean_model', 'NoSt', diag%axesTL, Time, &
      'Normal Stress', 's-1', conversion=US%s_to_T)
  CS%id_shearstress = register_diag_field('ocean_model', 'ShSt', diag%axesBL, Time, &
      'Shear Stress', 's-1', conversion=US%s_to_T)
  CS%id_diffu = register_diag_field('ocean_model', 'diffu', diag%axesCuL, Time, &
      'Zonal Acceleration from Horizontal Viscosity', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_diffv = register_diag_field('ocean_model', 'diffv', diag%axesCvL, Time, &
      'Meridional Acceleration from Horizontal Viscosity', 'm s-2', conversion=US%L_T2_to_m_s2)

  !CS%id_hf_diffu = register_diag_field('ocean_model', 'hf_diffu', diag%axesCuL, Time, &
  !    'Fractional Thickness-weighted Zonal Acceleration from Horizontal Viscosity', &
  !    'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if ((CS%id_hf_diffu > 0) .and. (present(ADp))) then
  !  call safe_alloc_alloc(CS%hf_diffu,G%IsdB,G%IedB,G%jsd,G%jed,GV%ke)
  !  call safe_alloc_ptr(ADp%diag_hfrac_u,G%IsdB,G%IedB,G%jsd,G%jed,GV%ke)
  !endif

  !CS%id_hf_diffv = register_diag_field('ocean_model', 'hf_diffv', diag%axesCvL, Time, &
  !    'Fractional Thickness-weighted Meridional Acceleration from Horizontal Viscosity', &
  !    'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if ((CS%id_hf_diffv > 0) .and. (present(ADp))) then
  !  call safe_alloc_alloc(CS%hf_diffv,G%isd,G%ied,G%JsdB,G%JedB,GV%ke)
  !  call safe_alloc_ptr(ADp%diag_hfrac_v,G%isd,G%ied,G%JsdB,G%JedB,GV%ke)
  !endif

  CS%id_hf_diffu_2d = register_diag_field('ocean_model', 'hf_diffu_2d', diag%axesCu1, Time, &
      'Depth-sum Fractional Thickness-weighted Zonal Acceleration from Horizontal Viscosity', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if ((CS%id_hf_diffu_2d > 0) .and. (present(ADp))) then
    call safe_alloc_ptr(ADp%diag_hfrac_u,G%IsdB,G%IedB,G%jsd,G%jed,GV%ke)
  endif

  CS%id_hf_diffv_2d = register_diag_field('ocean_model', 'hf_diffv_2d', diag%axesCv1, Time, &
      'Depth-sum Fractional Thickness-weighted Meridional Acceleration from Horizontal Viscosity', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if ((CS%id_hf_diffv_2d > 0) .and. (present(ADp))) then
    call safe_alloc_ptr(ADp%diag_hfrac_v,G%isd,G%ied,G%JsdB,G%JedB,GV%ke)
  endif

  CS%id_h_diffu = register_diag_field('ocean_model', 'h_diffu', diag%axesCuL, Time, &
      'Thickness Multiplied Zonal Acceleration from Horizontal Viscosity', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if ((CS%id_h_diffu > 0) .and. (present(ADp))) then
    call safe_alloc_ptr(ADp%diag_hu,G%IsdB,G%IedB,G%jsd,G%jed,GV%ke)
  endif

  CS%id_h_diffv = register_diag_field('ocean_model', 'h_diffv', diag%axesCvL, Time, &
      'Thickness Multiplied Meridional Acceleration from Horizontal Viscosity', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if ((CS%id_h_diffv > 0) .and. (present(ADp))) then
    call safe_alloc_ptr(ADp%diag_hv,G%isd,G%ied,G%JsdB,G%JedB,GV%ke)
  endif

  CS%id_intz_diffu_2d = register_diag_field('ocean_model', 'intz_diffu_2d', diag%axesCu1, Time, &
      'Depth-integral of Zonal Acceleration from Horizontal Viscosity', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if ((CS%id_intz_diffu_2d > 0) .and. (present(ADp))) then
    call safe_alloc_ptr(ADp%diag_hu,G%IsdB,G%IedB,G%jsd,G%jed,GV%ke)
  endif

  CS%id_intz_diffv_2d = register_diag_field('ocean_model', 'intz_diffv_2d', diag%axesCv1, Time, &
      'Depth-integral of Meridional Acceleration from Horizontal Viscosity', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if ((CS%id_intz_diffv_2d > 0) .and. (present(ADp))) then
    call safe_alloc_ptr(ADp%diag_hv,G%isd,G%ied,G%JsdB,G%JedB,GV%ke)
  endif

  CS%id_diffu_visc_rem = register_diag_field('ocean_model', 'diffu_visc_rem', diag%axesCuL, Time, &
      'Zonal Acceleration from Horizontal Viscosity multiplied by viscous remnant', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if ((CS%id_diffu_visc_rem > 0) .and. (present(ADp))) then
    call safe_alloc_ptr(ADp%visc_rem_u,G%IsdB,G%IedB,G%jsd,G%jed,GV%ke)
  endif

  CS%id_diffv_visc_rem = register_diag_field('ocean_model', 'diffv_visc_rem', diag%axesCvL, Time, &
      'Meridional Acceleration from Horizontal Viscosity multiplied by viscous remnant', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if ((CS%id_diffv_visc_rem > 0) .and. (present(ADp))) then
    call safe_alloc_ptr(ADp%visc_rem_v,G%isd,G%ied,G%JsdB,G%JedB,GV%ke)
  endif

  if (CS%biharmonic) then
    CS%id_Ah_h = register_diag_field('ocean_model', 'Ahh', diag%axesTL, Time,    &
        'Biharmonic Horizontal Viscosity at h Points', 'm4 s-1', conversion=US%L_to_m**4*US%s_to_T, &
        cmor_field_name='difmxybo',                                             &
        cmor_long_name='Ocean lateral biharmonic viscosity',                     &
        cmor_standard_name='ocean_momentum_xy_biharmonic_diffusivity')
    CS%id_Ah_q = register_diag_field('ocean_model', 'Ahq', diag%axesBL, Time, &
        'Biharmonic Horizontal Viscosity at q Points', 'm4 s-1', conversion=US%L_to_m**4*US%s_to_T)
    CS%id_grid_Re_Ah = register_diag_field('ocean_model', 'grid_Re_Ah', diag%axesTL, Time, &
        'Grid Reynolds number for the Biharmonic horizontal viscosity at h points', 'nondim')

    if (CS%id_grid_Re_Ah > 0) &
      ! Compute the smallest biharmonic viscosity capable of modifying the
      ! velocity at floating point precision.
      CS%min_grid_Ah = spacing(1.) * min_grid_sp_h4 * Idt
  endif
  if (CS%Laplacian) then
    CS%id_Kh_h = register_diag_field('ocean_model', 'Khh', diag%axesTL, Time,   &
        'Laplacian Horizontal Viscosity at h Points', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T, &
        cmor_field_name='difmxylo',                                             &
        cmor_long_name='Ocean lateral Laplacian viscosity',                     &
        cmor_standard_name='ocean_momentum_xy_laplacian_diffusivity')
    CS%id_Kh_q = register_diag_field('ocean_model', 'Khq', diag%axesBL, Time, &
        'Laplacian Horizontal Viscosity at q Points', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
    CS%id_grid_Re_Kh = register_diag_field('ocean_model', 'grid_Re_Kh', diag%axesTL, Time, &
        'Grid Reynolds number for the Laplacian horizontal viscosity at h points', 'nondim')
    CS%id_vort_xy_q = register_diag_field('ocean_model', 'vort_xy_q', diag%axesBL, Time, &
      'Vertical vorticity at q Points', 's-1', conversion=US%s_to_T)
    CS%id_div_xx_h = register_diag_field('ocean_model', 'div_xx_h', diag%axesTL, Time, &
      'Horizontal divergence at h Points', 's-1', conversion=US%s_to_T)
    CS%id_sh_xy_q = register_diag_field('ocean_model', 'sh_xy_q', diag%axesBL, Time, &
      'Shearing strain at q Points', 's-1', conversion=US%s_to_T)
    CS%id_sh_xx_h = register_diag_field('ocean_model', 'sh_xx_h', diag%axesTL, Time, &
      'Horizontal tension at h Points', 's-1', conversion=US%s_to_T)

    if (CS%id_grid_Re_Kh > 0) &
      ! Compute a smallest Laplacian viscosity capable of modifying the
      ! velocity at floating point precision.
      CS%min_grid_Kh = spacing(1.) * min_grid_sp_h2 * Idt
  endif
  if (CS%use_GME) then
    CS%id_dudx_bt = register_diag_field('ocean_model', 'dudx_bt', diag%axesT1, Time, &
        'Zonal component of the barotropic shearing strain at h points', 's-1', &
        conversion=US%s_to_T)
    CS%id_dudy_bt = register_diag_field('ocean_model', 'dudy_bt', diag%axesB1, Time, &
        'Zonal component of the barotropic shearing strain at q points', 's-1', &
        conversion=US%s_to_T)
    CS%id_dvdy_bt = register_diag_field('ocean_model', 'dvdy_bt', diag%axesT1, Time, &
        'Meridional component of the barotropic shearing strain at h points', 's-1', &
        conversion=US%s_to_T)
    CS%id_dvdx_bt = register_diag_field('ocean_model', 'dvdx_bt', diag%axesB1, Time, &
        'Meridional component of the barotropic shearing strain at q points', 's-1', &
        conversion=US%s_to_T)
    CS%id_GME_coeff_h = register_diag_field('ocean_model', 'GME_coeff_h', diag%axesTL, Time, &
        'GME coefficient at h Points', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
    CS%id_GME_coeff_q = register_diag_field('ocean_model', 'GME_coeff_q', diag%axesBL, Time, &
        'GME coefficient at q Points', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
    CS%id_FrictWork_GME = register_diag_field('ocean_model','FrictWork_GME',diag%axesTL,Time,&
        'Integral work done by lateral friction terms in GME (excluding diffusion of energy)', &
        'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  endif
  CS%id_FrictWork = register_diag_field('ocean_model','FrictWork',diag%axesTL,Time,&
      'Integral work done by lateral friction terms. If GME is turned on, this '//&
      'includes the GME contribution.', &
      'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  CS%id_FrictWorkIntz = register_diag_field('ocean_model','FrictWorkIntz',diag%axesT1,Time,      &
      'Depth integrated work done by lateral friction', &
      'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2, &
      cmor_field_name='dispkexyfo',                                                              &
      cmor_long_name='Depth integrated ocean kinetic energy dissipation due to lateral friction',&
      cmor_standard_name='ocean_kinetic_energy_dissipation_per_unit_area_due_to_xy_friction')

end subroutine hor_visc_init

!> Calculates factors in the anisotropic orientation tensor to be align with the grid.
!! With n1=1 and n2=0, this recovers the approach of Large et al, 2001.
subroutine align_aniso_tensor_to_grid(CS, n1, n2)
  type(hor_visc_CS), intent(inout) :: CS !< Control structure for horizontal viscosity
  real,              intent(in) :: n1 !< i-component of direction vector [nondim]
  real,              intent(in) :: n2 !< j-component of direction vector [nondim]
  ! Local variables
  real :: recip_n2_norm ! The inverse of the squared magnitude of n1 and n2 [nondim]
  ! For normalizing n=(n1,n2) in case arguments are not a unit vector
  recip_n2_norm = n1**2 + n2**2
  if (recip_n2_norm > 0.) recip_n2_norm = 1. / recip_n2_norm
  CS%n1n2_h(:,:) = 2. * ( n1 * n2 ) * recip_n2_norm
  CS%n1n2_q(:,:) = 2. * ( n1 * n2 ) * recip_n2_norm
  CS%n1n1_m_n2n2_h(:,:) = ( n1 * n1 - n2 * n2 ) * recip_n2_norm
  CS%n1n1_m_n2n2_q(:,:) = ( n1 * n1 - n2 * n2 ) * recip_n2_norm
end subroutine align_aniso_tensor_to_grid

!> Apply a 1-1-4-1-1 Laplacian filter one time on GME diffusive flux to reduce any
!! horizontal two-grid-point noise
subroutine smooth_GME(CS, G, GME_flux_h, GME_flux_q)
  type(hor_visc_CS),                            intent(in)    :: CS        !< Control structure
  type(ocean_grid_type),                        intent(in)    :: G         !< Ocean grid
  real, dimension(SZI_(G),SZJ_(G)),   optional, intent(inout) :: GME_flux_h!< GME diffusive flux
                                                              !! at h points [L2 T-2 ~> m2 s-2]
  real, dimension(SZIB_(G),SZJB_(G)), optional, intent(inout) :: GME_flux_q!< GME diffusive flux
                                                              !! at q points [L2 T-2 ~> m2 s-2]
  ! local variables
  real, dimension(SZI_(G),SZJ_(G)) :: GME_flux_h_original ! The previous value of GME_flux_h [L2 T-2 ~> m2 s-2]
  real, dimension(SZIB_(G),SZJB_(G)) :: GME_flux_q_original ! The previous value of GME_flux_q [L2 T-2 ~> m2 s-2]
  real :: wc, ww, we, wn, ws ! averaging weights for smoothing [nondim]
  integer :: i, j, s, halosz
  integer :: xh, xq  ! The number of valid extra halo points for h and q points.
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  xh = 0 ; xq = 0

  do s=1,CS%num_smooth_gme
    if (present(GME_flux_h)) then
      if (xh < 0) then
        ! Update halos if needed, but avoid doing so more often than is needed.
        halosz = min(G%isc-G%isd, G%jsc-G%jsd, 2+CS%num_smooth_gme-s)
        call pass_var(GME_flux_h, G%Domain, halo=halosz)
        xh = halosz - 2
      endif
      GME_flux_h_original(:,:) = GME_flux_h(:,:)
      ! apply smoothing on GME
      do j=Jsq-xh,Jeq+1+xh ; do i=Isq-xh,Ieq+1+xh
        ! skip land points
        if (G%mask2dT(i,j)==0.) cycle
        ! compute weights
        ww = 0.125 * G%mask2dT(i-1,j)
        we = 0.125 * G%mask2dT(i+1,j)
        ws = 0.125 * G%mask2dT(i,j-1)
        wn = 0.125 * G%mask2dT(i,j+1)
        wc = 1.0 - ((ww+we)+(wn+ws))
        GME_flux_h(i,j) =  wc * GME_flux_h_original(i,j)   &
                         + ((ww * GME_flux_h_original(i-1,j) + we * GME_flux_h_original(i+1,j)) &
                          + (ws * GME_flux_h_original(i,j-1) + wn * GME_flux_h_original(i,j+1)))
      enddo ; enddo
      xh = xh - 1
    endif
    if (present(GME_flux_q)) then
      if (xq < 0) then
        ! Update halos if needed, but avoid doing so more often than is needed.
        halosz = min(G%isc-G%isd, G%jsc-G%jsd, 2+CS%num_smooth_gme-s)
        call pass_var(GME_flux_q, G%Domain, position=CORNER, complete=.true., halo=halosz)
        xq = halosz - 2
      endif
      GME_flux_q_original(:,:) = GME_flux_q(:,:)
      ! apply smoothing on GME
      do J=js-1-xq,je+xq ; do I=is-1-xq,ie+xq
        ! skip land points
        if (G%mask2dBu(I,J)==0.) cycle
        ! compute weights
        ww = 0.125 * G%mask2dBu(I-1,J)
        we = 0.125 * G%mask2dBu(I+1,J)
        ws = 0.125 * G%mask2dBu(I,J-1)
        wn = 0.125 * G%mask2dBu(I,J+1)
        wc = 1.0 - ((ww+we)+(wn+ws))
        GME_flux_q(I,J) =  wc * GME_flux_q_original(I,J)   &
                         + ((ww * GME_flux_q_original(I-1,J) + we * GME_flux_q_original(I+1,J)) &
                          + (ws * GME_flux_q_original(I,J-1) + wn * GME_flux_q_original(I,J+1)))
      enddo ; enddo
      xq = xq - 1
    endif
  enddo ! s-loop
end subroutine smooth_GME

!> Apply a 9-point smoothing filter twice to a field staggered at a thickness point to reduce
!! horizontal two-grid-point noise.
!! Note that this subroutine does not conserve mass, so don't use it in situations where you
!! need conservation.  Also note that it assumes that the input field has valid values in the
!! first two halo points upon entry.
subroutine smooth_x9_h(G, field_h, zero_land)
  type(ocean_grid_type),            intent(in)    :: G         !< Ocean grid
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: field_h   !< h-point field to be smoothed [arbitrary]
  logical,                optional, intent(in)    :: zero_land !< If present and false, return the average
                                                               !! of the surrounding ocean points when
                                                               !! smoothing, otherwise use a value of 0 for
                                                               !! land points and include them in the averages.

  ! Local variables
  real :: fh_prev(SZI_(G),SZJ_(G))  ! The value of the h-point field at the previous iteration [arbitrary]
  real :: Iwts             ! The inverse of the sum of the weights [nondim]
  logical :: zero_land_val ! The value of the zero_land optional argument or .true. if it is absent.
  integer :: i, j, s, is, ie, js, je

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec

  zero_land_val = .true. ; if (present(zero_land)) zero_land_val = zero_land

  do s=1,0,-1
    fh_prev(:,:) = field_h(:,:)
    ! apply smoothing on field_h using rotationally symmetric expressions.
    do j=js-s,je+s ; do i=is-s,ie+s ; if (G%mask2dT(i,j) > 0.0) then
      Iwts = 0.0625
      if (.not. zero_land_val) &
        Iwts = 1.0 / ( (4.0*G%mask2dT(i,j) + &
                        ( 2.0*((G%mask2dT(i-1,j) + G%mask2dT(i+1,j)) + &
                               (G%mask2dT(i,j-1) + G%mask2dT(i,j+1))) + &
                         ((G%mask2dT(i-1,j-1) + G%mask2dT(i+1,j+1)) + &
                          (G%mask2dT(i-1,j+1) + G%mask2dT(i+1,j-1))) ) ) + 1.0e-16 )
      field_h(i,j) = Iwts * ( 4.0*G%mask2dT(i,j) * fh_prev(i,j) &
                            + (2.0*((G%mask2dT(i-1,j) * fh_prev(i-1,j) + G%mask2dT(i+1,j) * fh_prev(i+1,j)) + &
                                    (G%mask2dT(i,j-1) * fh_prev(i,j-1) + G%mask2dT(i,j+1) * fh_prev(i,j+1))) &
                              + ((G%mask2dT(i-1,j-1) * fh_prev(i-1,j-1) + G%mask2dT(i+1,j+1) * fh_prev(i+1,j+1)) + &
                                 (G%mask2dT(i-1,j+1) * fh_prev(i-1,j+1) + G%mask2dT(i+1,j-1) * fh_prev(i-1,j-1))) ))
    endif ; enddo ; enddo
  enddo

end subroutine smooth_x9_h

!> Apply a 9-point smoothing filter twice to a pair of velocity components to reduce
!! horizontal two-grid-point noise.
!! Note that this subroutine does not conserve angular momentum, so don't use it
!! in situations where you need conservation.  Also note that it assumes that the
!! input fields have valid values in the first two halo points upon entry.
subroutine smooth_x9_uv(G, field_u, field_v, zero_land)
  type(ocean_grid_type),             intent(in)    :: G         !< Ocean grid
  real, dimension(SZIB_(G),SZJ_(G)), intent(inout) :: field_u   !< u-point field to be smoothed [arbitrary]
  real, dimension(SZI_(G),SZJB_(G)), intent(inout) :: field_v   !< v-point field to be smoothed [arbitrary]
  logical,                 optional, intent(in)    :: zero_land !< If present and false, return the average
                                                                !! of the surrounding ocean points when
                                                                !! smoothing, otherwise use a value of 0 for
                                                                !! land points and include them in the averages.

  ! Local variables.
  real :: fu_prev(SZIB_(G),SZJ_(G))  ! The value of the u-point field at the previous iteration [arbitrary]
  real :: fv_prev(SZI_(G),SZJB_(G))  ! The value of the v-point field at the previous iteration [arbitrary]
  real :: Iwts             ! The inverse of the sum of the weights [nondim]
  logical :: zero_land_val ! The value of the zero_land optional argument or .true. if it is absent.
  integer :: i, j, s, is, ie, js, je, Isq, Ieq, Jsq, Jeq

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  zero_land_val = .true. ; if (present(zero_land)) zero_land_val = zero_land

  do s=1,0,-1
    fu_prev(:,:) = field_u(:,:)
    ! apply smoothing on field_u using the original non-rotationally symmetric expressions.
    do j=js-s,je+s ; do I=Isq-s,Ieq+s ; if (G%mask2dCu(I,j) > 0.0) then
      Iwts = 0.0625
      if (.not. zero_land_val) &
        Iwts = 1.0 / ( (4.0*G%mask2dCu(I,j) + &
                        ( 2.0*((G%mask2dCu(I-1,j) + G%mask2dCu(I+1,j)) + &
                               (G%mask2dCu(I,j-1) + G%mask2dCu(I,j+1))) + &
                         ((G%mask2dCu(I-1,j-1) + G%mask2dCu(I+1,j+1)) + &
                          (G%mask2dCu(I-1,j+1) + G%mask2dCu(I+1,j-1))) ) ) + 1.0e-16 )
      field_u(I,j) = Iwts * ( 4.0*G%mask2dCu(I,j) * fu_prev(I,j) &
                            + (2.0*((G%mask2dCu(I-1,j) * fu_prev(I-1,j) + G%mask2dCu(I+1,j) * fu_prev(I+1,j)) + &
                                    (G%mask2dCu(I,j-1) * fu_prev(I,j-1) + G%mask2dCu(I,j+1) * fu_prev(I,j+1))) &
                              + ((G%mask2dCu(I-1,j-1) * fu_prev(I-1,j-1) + G%mask2dCu(I+1,j+1) * fu_prev(I+1,j+1)) + &
                                 (G%mask2dCu(I-1,j+1) * fu_prev(I-1,j+1) + G%mask2dCu(I+1,j-1) * fu_prev(I-1,j-1))) ))
    endif ; enddo ; enddo

    fv_prev(:,:) = field_v(:,:)
    ! apply smoothing on field_v using the original non-rotationally symmetric expressions.
    do J=Jsq-s,Jeq+s ; do i=is-s,ie+s ; if (G%mask2dCv(i,J) > 0.0) then
      Iwts = 0.0625
      if (.not. zero_land_val) &
        Iwts = 1.0 / ( (4.0*G%mask2dCv(i,J) + &
                        ( 2.0*((G%mask2dCv(i-1,J) + G%mask2dCv(i+1,J)) + &
                               (G%mask2dCv(i,J-1) + G%mask2dCv(i,J+1))) + &
                         ((G%mask2dCv(i-1,J-1) + G%mask2dCv(i+1,J+1)) + &
                          (G%mask2dCv(i-1,J+1) + G%mask2dCv(i+1,J-1))) ) ) + 1.0e-16 )
      field_v(i,J) = Iwts * ( 4.0*G%mask2dCv(i,J) * fv_prev(i,J) &
                            + (2.0*((G%mask2dCv(i-1,J) * fv_prev(i-1,J) + G%mask2dCv(i+1,J) * fv_prev(i+1,J)) + &
                                    (G%mask2dCv(i,J-1) * fv_prev(i,J-1) + G%mask2dCv(i,J+1) * fv_prev(i,J+1))) &
                              + ((G%mask2dCv(i-1,J-1) * fv_prev(i-1,J-1) + G%mask2dCv(i+1,J+1) * fv_prev(i+1,J+1)) + &
                                 (G%mask2dCv(i-1,J+1) * fv_prev(i-1,J+1) + G%mask2dCv(i+1,J-1) * fv_prev(i-1,J-1))) ))
    endif ; enddo ; enddo
  enddo

end subroutine smooth_x9_uv

!> Deallocates any variables allocated in hor_visc_init.
subroutine hor_visc_end(CS)
  type(hor_visc_CS), intent(inout) :: CS !< Horizontal viscosity control structure
  if (CS%Laplacian .or. CS%biharmonic) then
    DEALLOC_(CS%dx2h) ; DEALLOC_(CS%dx2q) ; DEALLOC_(CS%dy2h) ; DEALLOC_(CS%dy2q)
    DEALLOC_(CS%dx_dyT) ; DEALLOC_(CS%dy_dxT) ; DEALLOC_(CS%dx_dyBu) ; DEALLOC_(CS%dy_dxBu)
    DEALLOC_(CS%reduction_xx) ; DEALLOC_(CS%reduction_xy)
  endif
  if (CS%Laplacian) then
    DEALLOC_(CS%Kh_bg_xx) ; DEALLOC_(CS%Kh_bg_xy)
    DEALLOC_(CS%grid_sp_h2)
    if (CS%bound_Kh) then
      DEALLOC_(CS%Kh_Max_xx) ; DEALLOC_(CS%Kh_Max_xy)
    endif
    if (CS%Smagorinsky_Kh) then
      DEALLOC_(CS%Laplac2_const_xx) ; DEALLOC_(CS%Laplac2_const_xy)
    endif
    if (CS%Leith_Kh) then
      DEALLOC_(CS%Laplac3_const_xx) ; DEALLOC_(CS%Laplac3_const_xy)
    endif
  endif
  if (CS%biharmonic) then
    DEALLOC_(CS%grid_sp_h3)
    DEALLOC_(CS%Idx2dyCu) ; DEALLOC_(CS%Idx2dyCv)
    DEALLOC_(CS%Idxdy2u) ; DEALLOC_(CS%Idxdy2v)
    DEALLOC_(CS%Ah_bg_xx) ; DEALLOC_(CS%Ah_bg_xy)
    if (CS%bound_Ah) then
      DEALLOC_(CS%Ah_Max_xx) ; DEALLOC_(CS%Ah_Max_xy)
    endif
    if (CS%Smagorinsky_Ah) then
      DEALLOC_(CS%Biharm_const_xx) ; DEALLOC_(CS%Biharm_const_xy)
    endif
    if ((CS%Leith_Ah) .or. (CS%use_Leithy)) then
      DEALLOC_(CS%Biharm6_const_xx) ; DEALLOC_(CS%Biharm6_const_xy)
    endif
    if (CS%use_Leithy) then
      DEALLOC_(CS%m_const_leithy)
      DEALLOC_(CS%m_leithy_max)
    endif
    if (CS%Re_Ah > 0.0) then
      DEALLOC_(CS%Re_Ah_const_xx) ; DEALLOC_(CS%Re_Ah_const_xy)
    endif
  endif
  if (CS%anisotropic) then
    DEALLOC_(CS%n1n2_h)
    DEALLOC_(CS%n1n2_q)
    DEALLOC_(CS%n1n1_m_n2n2_h)
    DEALLOC_(CS%n1n1_m_n2n2_q)
  endif

  if (CS%use_ZB2020) then
    call ZB2020_end(CS%ZB2020)
  endif

end subroutine hor_visc_end
!> \namespace mom_hor_visc
!!
!! This module contains the subroutine horizontal_viscosity() that calculates the
!! effects of horizontal viscosity, including parameterizations of the value of
!! the viscosity itself. horizontal_viscosity() calculates the acceleration due to
!! some combination of a biharmonic viscosity and a Laplacian viscosity. Either or
!! both may use a coefficient that depends on the shear and strain of the flow.
!! All metric terms are retained. The Laplacian is calculated as the divergence of
!! a stress tensor, using the form suggested by Smagorinsky (1993). The biharmonic
!! is calculated by twice applying the divergence of the stress tensor that is
!! used to calculate the Laplacian, but without the dependence on thickness in the
!! first pass. This form permits a variable viscosity, and indicates no
!! acceleration for either resting fluid or solid body rotation.
!!
!! The form of the viscous accelerations is discussed extensively in Griffies and
!! Hallberg (2000), and the implementation here follows that discussion closely.
!! We use the notation of Smith and McWilliams (2003) with the exception that the
!! isotropic viscosity is \f$\kappa_h\f$.
!!
!! \section section_horizontal_viscosity Horizontal viscosity in MOM
!!
!! In general, the horizontal stress tensor can be written as
!! \f[
!! {\bf \sigma} =
!! \begin{pmatrix}
!! \frac{1}{2} \left( \sigma_D + \sigma_T \right) & \frac{1}{2} \sigma_S \\\\
!! \frac{1}{2} \sigma_S & \frac{1}{2} \left( \sigma_D - \sigma_T \right)
!! \end{pmatrix}
!! \f]
!! where \f$\sigma_D\f$, \f$\sigma_T\f$ and \f$\sigma_S\f$ are stresses associated with
!! invariant factors in the strain-rate tensor. For a Newtonian fluid, the stress
!! tensor is usually linearly related to the strain-rate tensor. The horizontal
!! strain-rate tensor is
!! \f[
!! \dot{\bf e} =
!! \begin{pmatrix}
!! \frac{1}{2} \left( \dot{e}_D + \dot{e}_T \right) & \frac{1}{2} \dot{e}_S \\\\
!! \frac{1}{2} \dot{e}_S & \frac{1}{2} \left( \dot{e}_D - \dot{e}_T \right)
!! \end{pmatrix}
!! \f]
!! where \f$\dot{e}_D = \partial_x u + \partial_y v\f$ is the horizontal divergence,
!! \f$\dot{e}_T = \partial_x u - \partial_y v\f$ is the horizontal tension, and
!! \f$\dot{e}_S = \partial_y u + \partial_x v\f$ is the horizontal shear strain.
!!
!! The trace of the stress tensor, \f$tr(\bf \sigma) = \sigma_D\f$, is usually
!! absorbed into the pressure and only the deviatoric stress tensor considered.
!! From here on, we drop \f$\sigma_D\f$. The trace of the strain tensor, \f$tr(\bf e) =
!! \dot{e}_D\f$ is non-zero for horizontally divergent flow but only enters the
!! stress tensor through \f$\sigma_D\f$ and so we will drop \f$\sigma_D\f$ from
!! calculations of the strain tensor in the code. Therefore the horizontal stress
!! tensor can be considered to be
!! \f[
!! {\bf \sigma} =
!! \begin{pmatrix}
!! \frac{1}{2} \sigma_T & \frac{1}{2} \sigma_S \\\\
!! \frac{1}{2} \sigma_S & - \frac{1}{2} \sigma_T
!! \end{pmatrix}
!! .\f]
!!
!! The stresses above are linearly related to the strain through a viscosity
!! coefficient, \f$\kappa_h\f$:
!! \f{eqnarray*}{
!! \sigma_T & = & 2 \kappa_h \dot{e}_T \\\\
!! \sigma_S & = & 2 \kappa_h \dot{e}_S
!! .
!! \f}
!!
!! The viscosity \f$\kappa_h\f$ may either be a constant or variable. For example,
!! \f$\kappa_h\f$ may vary with the shear, as proposed by Smagorinsky (1993).
!!
!! The accelerations resulting form the divergence of the stress tensor are
!! \f{eqnarray*}{
!! \hat{\bf x} \cdot \left( \nabla \cdot {\bf \sigma} \right)
!! & = &
!! \partial_x \left( \frac{1}{2} \sigma_T \right)
!! + \partial_y \left( \frac{1}{2} \sigma_S \right)
!! \\\\
!! & = &
!! \partial_x \left( \kappa_h \dot{e}_T \right)
!! + \partial_y \left( \kappa_h \dot{e}_S \right)
!! \\\\
!! \hat{\bf y} \cdot \left( \nabla \cdot {\bf \sigma} \right)
!! & = &
!! \partial_x \left( \frac{1}{2} \sigma_S \right)
!! + \partial_y \left( - \frac{1}{2} \sigma_T \right)
!! \\\\
!! & = &
!! \partial_x \left( \kappa_h \dot{e}_S \right)
!! + \partial_y \left( - \kappa_h \dot{e}_T \right)
!! .
!! \f}
!!
!! The form of the Laplacian viscosity in general coordinates is:
!! \f{eqnarray*}{
!! \hat{\bf x} \cdot \left( \nabla \cdot \sigma \right)
!! & = &
!! \frac{1}{h} \left[ \partial_x \left( \kappa_h h \dot{e}_T \right)
!! + \partial_y \left( \kappa_h h \dot{e}_S \right) \right]
!! \\\\
!! \hat{\bf y} \cdot \left( \nabla \cdot \sigma \right)
!! & = &
!! \frac{1}{h} \left[ \partial_x \left( \kappa_h h \dot{e}_S \right)
!! - \partial_y \left( \kappa_h h \dot{e}_T \right) \right]
!! .
!! \f}
!!
!! \subsection section_laplacian_viscosity_coefficient Laplacian viscosity coefficient
!!
!! The horizontal viscosity coefficient, \f$\kappa_h\f$, can have multiple components.
!! The isotropic components are:
!!   - A uniform background component, \f$\kappa_{bg}\f$.
!!   - A constant but spatially variable 2D map, \f$\kappa_{2d}(x,y)\f$.
!!   - A ''MICOM'' viscosity, \f$U_\nu \Delta(x,y)\f$, which uses a constant
!! velocity scale, \f$U_\nu\f$ and a measure of the grid-spacing \f$\Delta(x,y)^2 =
!! \frac{2 \Delta x^2 \Delta y^2}{\Delta x^2 + \Delta y^2}\f$.
!!   - A function of
!! latitude, \f$\kappa_{\phi}(x,y) = \kappa_{\pi/2} |\sin(\phi)|^n\f$.
!!   - A dynamic Smagorinsky viscosity, \f$\kappa_{Sm}(x,y,t) = C_{Sm} \Delta^2 \sqrt{\dot{e}_T^2 + \dot{e}_S^2}\f$.
!!   - A dynamic Leith viscosity, \f$\kappa_{Lth}(x,y,t) =
!!                                    C_{Lth} \Delta^3 \sqrt{|\nabla \zeta|^2 + |\nabla \dot{e}_D|^2}\f$.
!!
!! A maximum stable viscosity, \f$\kappa_{max}(x,y)\f$ is calculated based on the
!! grid-spacing and time-step and used to clip calculated viscosities.
!!
!! The static components of \f$\kappa_h\f$ are first combined as follows:
!! \f[
!! \kappa_{static} = \min \left[ \max\left(
!! \kappa_{bg},
!! U_\nu \Delta(x,y),
!! \kappa_{2d}(x,y),
!! \kappa_\phi(x,y)
!! \right)
!! , \kappa_{max}(x,y) \right]
!! \f]
!! and stored in the module control structure as variables <code>Kh_bg_xx</code> and
!! <code>Kh_bg_xy</code> for the tension (h-points) and shear (q-points) components
!! respectively.
!!
!! The full viscosity includes the dynamic components as follows:
!! \f[
!! \kappa_h(x,y,t) = r(\Delta,L_d)
!! \max \left( \kappa_{static}, \kappa_{Sm}, \kappa_{Lth} \right)
!! \f]
!! where \f$r(\Delta,L_d)\f$ is a resolution function.
!!
!! The dynamic Smagorinsky and Leith viscosity schemes are exclusive with each
!! other.
!!
!! \subsection section_viscous_boundary_conditions Viscous boundary conditions
!!
!! Free slip boundary conditions have been coded, although no slip boundary
!! conditions can be used with the Laplacian viscosity based on the 2D land-sea
!! mask. For a western boundary, for example, the boundary conditions with the
!! biharmonic operator would be written as:
!! \f[
!!   \partial_x v = 0 ; \partial_x^3 v = 0 ; u = 0 ; \partial_x^2 u = 0 ,
!! \f]
!! while for a Laplacian operator, they are simply
!! \f[
!!   \partial_x v = 0 ; u = 0 .
!! \f]
!! These boundary conditions are largely dictated by the use of an Arakawa
!! C-grid and by the varying layer thickness.
!!
!! \subsection section_anisotropic_viscosity Anisotropic viscosity
!!
!! Large et al., 2001, proposed enhancing viscosity in a particular direction and the
!! approach was generalized in Smith and McWilliams, 2003. We use the second form of their
!! two coefficient anisotropic viscosity (section 4.3). We also replace their
!! \f$A^\prime\f$ and $D$ such that \f$2A^\prime = 2 \kappa_h + D\f$ and
!! \f$\kappa_a = D\f$ so that \f$\kappa_h\f$ can be considered the isotropic
!! viscosity and \f$\kappa_a=D\f$ can be consider the anisotropic viscosity. The
!! direction of anisotropy is defined by a unit vector \f$\hat{\bf
!! n}=(n_1,n_2)\f$.
!!
!! The contributions to the stress tensor are
!! \f[
!! \begin{pmatrix}
!! \sigma_T \\\\ \sigma_S
!! \end{pmatrix}
!! =
!! \left[
!! \begin{pmatrix}
!! 2 \kappa_h + \kappa_a & 0 \\\\
!! 0 & 2 \kappa_h
!! \end{pmatrix}
!! + 2 \kappa_a n_1 n_2
!! \begin{pmatrix}
!! - 2 n_1 n_2 & n_1^2 - n_2^2 \\\\
!! n_1^2 - n_2^2 & 2 n_1 n_2
!! \end{pmatrix}
!! \right]
!! \begin{pmatrix}
!! \dot{e}_T \\\\ \dot{e}_S
!! \end{pmatrix}
!! \f]
!! Dissipation of kinetic energy requires \f$\kappa_h \geq 0\f$ and \f$2 \kappa_h + \kappa_a \geq 0\f$.
!! Note that when anisotropy is aligned with the x-direction, \f$n_1 = \pm 1\f$, then
!! \f$n_2 = 0\f$ and the cross terms vanish. The accelerations in this aligned limit
!! with constant coefficients become
!! \f{eqnarray*}{
!! \hat{\bf x} \cdot \nabla \cdot {\bf \sigma}
!! & = &
!! \partial_x \left( \left( \kappa_h + \frac{1}{2} \kappa_a \right) \dot{e}_T \right)
!! + \partial_y \left( \kappa_h \dot{e}_S \right)
!! \\\\
!! & = &
!! \left( \kappa_h + \kappa_a \right) \partial_{xx} u
!! + \kappa_h \partial_{yy} u
!! - \frac{1}{2} \kappa_a \partial_x \left( \partial_x u + \partial_y v \right)
!! \\\\
!! \hat{\bf y} \cdot \nabla \cdot {\bf \sigma}
!! & = &
!! \partial_x \left( \kappa_h \dot{e}_S \right)
!! - \partial_y \left( \left( \kappa_h + \frac{1}{2} \kappa_a \right) \dot{e}_T \right)
!! \\\\
!! & = &
!! \kappa_h \partial_{xx} v
!! + \left( \kappa_h + \kappa_a \right) \partial_{yy} v
!! - \frac{1}{2} \kappa_a \partial_y \left( \partial_x u + \partial_y v \right)
!! \f}
!! which has contributions akin to a negative divergence damping (a divergence
!! enhancement?) but which is weaker than the enhanced tension terms by half.
!!
!! \subsection section_viscous_discretization Discretization
!!
!! The horizontal tension, \f$\dot{e}_T\f$, is stored in variable <code>sh_xx</code> and
!! discretized as
!! \f[
!! \dot{e}_T
!! = \frac{\Delta y}{\Delta x} \delta_i \left( \frac{1}{\Delta y} u \right)
!! - \frac{\Delta x}{\Delta y} \delta_j \left( \frac{1}{\Delta x} v \right)
!! .
!! \f]
!! The horizontal divergent strain, \f$\dot{e}_D\f$, is stored in variable
!! <code>div_xx</code> and discretized as
!! \f[
!! \dot{e}_D
!! = \frac{1}{h A} \left( \delta_i \left( \overline{h}^i \Delta y \, u \right)
!! + \delta_j \left( \overline{h}^j\Delta x \, v \right) \right)
!! .
!! \f]
!! Note that for expediency this is the exact discretization used in the
!! continuity equation.
!!
!! The horizontal shear strain, \f$\dot{e}_S\f$, is stored in variable <code>sh_xy</code>
!! and discretized as
!! \f[
!! \dot{e}_S = v_x + u_y
!! \f]
!! where
!! \f{align*}{
!! v_x &= \frac{\Delta y}{\Delta x} \delta_i \left( \frac{1}{\Delta y} v \right) \\\\
!! u_y &= \frac{\Delta x}{\Delta y} \delta_j \left( \frac{1}{\Delta x} u \right)
!! \f}
!! which are calculated separately so that no-slip or free-slip boundary
!! conditions can be applied to \f$v_x\f$ and \f$u_y\f$ where appropriate.
!!
!! The tendency for the x-component of the divergence of stress is stored in
!! variable <code>diffu</code> and discretized as
!! \f[
!! \hat{\bf x} \cdot \left( \nabla \cdot {\bf \sigma} \right) =
!! \frac{1}{A \overline{h}^i} \left(
!! \frac{1}{\Delta y} \delta_i \left( h \Delta y^2 \kappa_h \dot{e}_T \right) +
!! \frac{1}{\Delta x} \delta_j \left( \tilde{h}^{ij} \Delta x^2 \kappa_h \dot{e}_S \right)
!! \right)
!! .
!! \f]
!!
!! The tendency for the y-component of the divergence of stress is stored in
!! variable <code>diffv</code> and discretized as
!! \f[
!! \hat{\bf y} \cdot \left( \nabla \cdot {\bf \sigma} \right) =
!! \frac{1}{A \overline{h}^j} \left(
!! \frac{1}{\Delta y} \delta_i \left( \tilde{h}^{ij} \Delta y^2 A_M \dot{e}_S \right)
!! - \frac{1}{\Delta x} \delta_j \left( h \Delta x^2 A_M \dot{e}_T \right)
!! \right)
!! .
!! \f]
!!
!! \subsection section_viscous_refs References
!!
!! Griffies, S.M., and Hallberg, R.W., 2000: Biharmonic friction with a
!! Smagorinsky-like viscosity for use in large-scale eddy-permitting ocean models.
!! Monthly Weather Review, 128(8), 2935-2946.
!! https://doi.org/10.1175/1520-0493(2000)128%3C2935:BFWASL%3E2.0.CO;2
!!
!! Large, W.G., Danabasoglu, G., McWilliams, J.C., Gent, P.R. and Bryan, F.O.,
!! 2001: Equatorial circulation of a global ocean climate model with
!! anisotropic horizontal viscosity.
!! Journal of Physical Oceanography, 31(2), pp.518-536.
!! https://doi.org/10.1175/1520-0485(2001)031%3C0518:ECOAGO%3E2.0.CO;2
!!
!! Smagorinsky, J., 1993: Some historical remarks on the use of nonlinear
!! viscosities. Large eddy simulation of complex engineering and geophysical
!! flows, 1, 69-106.
!!
!! Smith, R.D., and McWilliams, J.C., 2003: Anisotropic horizontal viscosity for
!! ocean models. Ocean Modelling, 5(2), 129-156.
!! https://doi.org/10.1016/S1463-5003(02)00016-1
end module MOM_hor_visc
