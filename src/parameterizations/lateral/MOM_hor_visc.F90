!> Calculates horizontal viscosity and viscous stresses
module MOM_hor_visc

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_checksums,             only : hchksum, Bchksum
use MOM_coms,                  only : min_across_PEs
use MOM_diag_mediator,         only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator,         only : diag_ctrl, time_type
use MOM_domains,               only : pass_var, CORNER, pass_vector, AGRID, BGRID_NE
use MOM_error_handler,         only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_grid,                  only : ocean_grid_type
use MOM_lateral_mixing_coeffs, only : VarMix_CS, calc_QG_Leith_viscosity
use MOM_barotropic,            only : barotropic_CS, barotropic_get_tav
use MOM_thickness_diffuse,     only : thickness_diffuse_CS, thickness_diffuse_get_KH
use MOM_io,                    only : MOM_read_data, slasher
use MOM_MEKE_types,            only : MEKE_type
use MOM_open_boundary,         only : ocean_OBC_type, OBC_DIRECTION_E, OBC_DIRECTION_W
use MOM_open_boundary,         only : OBC_DIRECTION_N, OBC_DIRECTION_S, OBC_NONE
use MOM_unit_scaling,          only : unit_scale_type
use MOM_verticalGrid,          only : verticalGrid_type
use MOM_variables,             only : accel_diag_ptrs

implicit none ; private

#include <MOM_memory.h>

public horizontal_viscosity, hor_visc_init, hor_visc_end

!> Control structure for horizontal viscosity
type, public :: hor_visc_CS ; private
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
                             !< so that the biharmonic Reynolds number is equal to this.
  real    :: bound_coef      !< The nondimensional coefficient of the ratio of
                             !! the viscosity bounds to the theoretical maximum
                             !! for stability without considering other terms [nondim].
                             !! The default is 0.8.
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
  logical :: use_QG_Leith_visc    !< If true, use QG Leith nonlinear eddy viscosity.
                             !! KH is the background value.
  logical :: bound_Coriolis  !< If true & SMAGORINSKY_AH is used, the biharmonic
                             !! viscosity is modified to include a term that
                             !! scales quadratically with the velocity shears.
  logical :: use_Kh_bg_2d    !< Read 2d background viscosity from a file.
  real    :: Kh_bg_min       !< The minimum value allowed for Laplacian horizontal
                             !! viscosity [L2 T-1 ~> m2 s-1]. The default is 0.0.
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
  logical :: answers_2018    !< If true, use the order of arithmetic and expressions that recover the
                             !! answers from the end of 2018.  Otherwise, use updated and more robust
                             !! forms of the same expressions.
  real    :: GME_h0          !< The strength of GME tapers quadratically to zero when the bathymetric
                             !! depth is shallower than GME_H0 [Z ~> m]
  real    :: GME_efficiency  !< The nondimensional prefactor multiplying the GME coefficient [nondim]
  real    :: GME_limiter     !< The absolute maximum value the GME coefficient is allowed to take [L2 T-1 ~> m2 s-1].
  real    :: min_grid_Kh     !< Minimum horizontal Laplacian viscosity used to
                             !! limit the grid Reynolds number [L2 T-1 ~> m2 s-1]
  real    :: min_grid_Ah     !< Minimun horizontal biharmonic viscosity used to
                             !! limit grid Reynolds number [L4 T-1 ~> m4 s-1]

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
    n1n2_h,         & !< Factor n1*n2 in the anisotropic direction tensor at h-points
    n1n1_m_n2n2_h,  & !< Factor n1**2-n2**2 in the anisotropic direction tensor at h-points
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
    n1n2_q,         & !< Factor n1*n2 in the anisotropic direction tensor at q-points
    n1n1_m_n2n2_q     !< Factor n1**2-n2**2 in the anisotropic direction tensor at q-points

  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    dx2h,   & !< Pre-calculated dx^2 at h points [L2 ~> m2]
    dy2h,   & !< Pre-calculated dy^2 at h points [L2 ~> m2]
    dx_dyT, & !< Pre-calculated dx/dy at h points [nondim]
    dy_dxT    !< Pre-calculated dy/dx at h points [nondim]
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

  ! real, pointer :: hf_diffu(:,:,:)   => NULL() ! Zonal hor. visc. accel. x fract. thickness [L T-2 ~> m s-2].
  ! real, pointer :: hf_diffv(:,:,:)   => NULL() ! Merdional hor. visc. accel. x fract. thickness [L T-2 ~> m s-2].
  ! 3D diagnostics hf_diffu(diffv) are commented because there is no clarity on proper remapping grid option.
  ! The code is retained for degugging purposes in the future.

  !>@{
  !! Diagnostic id
  integer :: id_grid_Re_Ah = -1, id_grid_Re_Kh   = -1
  integer :: id_diffu     = -1, id_diffv         = -1
  ! integer :: id_hf_diffu  = -1, id_hf_diffv      = -1
  integer :: id_hf_diffu_2d = -1, id_hf_diffv_2d = -1
  integer :: id_Ah_h      = -1, id_Ah_q          = -1
  integer :: id_Kh_h      = -1, id_Kh_q          = -1
  integer :: id_GME_coeff_h = -1, id_GME_coeff_q = -1
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
!!   u[is-2:ie+2,js-2:je+2]
!!   v[is-2:ie+2,js-2:je+2]
!!   h[is-1:ie+1,js-1:je+1]
subroutine horizontal_viscosity(u, v, h, diffu, diffv, MEKE, VarMix, G, GV, US, &
                                CS, OBC, BT, TD, ADp)
  type(ocean_grid_type),         intent(in)  :: G      !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)  :: GV     !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                                 intent(in)  :: u      !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                                 intent(in)  :: v      !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  &
                                 intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                                 intent(out) :: diffu  !< Zonal acceleration due to convergence of
                                                       !! along-coordinate stress tensor [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                                 intent(out) :: diffv  !< Meridional acceleration due to convergence
                                                       !! of along-coordinate stress tensor [L T-2 ~> m s-2].
  type(MEKE_type),               pointer     :: MEKE   !< Pointer to a structure containing fields
                                                       !! related to Mesoscale Eddy Kinetic Energy.
  type(VarMix_CS),               pointer     :: VarMix !< Pointer to a structure with fields that
                                                       !! specify the spatially variable viscosities
  type(unit_scale_type),         intent(in)  :: US     !< A dimensional unit scaling type
  type(hor_visc_CS),             pointer     :: CS     !< Control structure returned by a previous
                                                       !! call to hor_visc_init.
  type(ocean_OBC_type), optional, pointer    :: OBC    !< Pointer to an open boundary condition type
  type(barotropic_CS),  optional, pointer    :: BT     !< Pointer to a structure containing
                                                       !! barotropic velocities.
  type(thickness_diffuse_CS), optional, pointer :: TD  !< Pointer to a structure containing
                                                       !! thickness diffusivities.
  type(accel_diag_ptrs), optional, pointer :: ADp      !< Acceleration diagnostic pointers

  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    Del2u, &      ! The u-compontent of the Laplacian of velocity [L-1 T-1 ~> m-1 s-1]
    h_u, &        ! Thickness interpolated to u points [H ~> m or kg m-2].
    vort_xy_dy, & ! y-derivative of vertical vorticity (d/dy(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
    div_xx_dx, &  ! x-derivative of horizontal divergence (d/dx(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
    ubtav         ! zonal barotropic vel. ave. over baroclinic time-step [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G)) :: &
    Del2v, &      ! The v-compontent of the Laplacian of velocity [L-1 T-1 ~> m-1 s-1]
    h_v, &        ! Thickness interpolated to v points [H ~> m or kg m-2].
    vort_xy_dx, & ! x-derivative of vertical vorticity (d/dx(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
    div_xx_dy, &  ! y-derivative of horizontal divergence (d/dy(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
    vbtav         ! meridional barotropic vel. ave. over baroclinic time-step [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G)) :: &
    dudx_bt, dvdy_bt, & ! components in the barotropic horizontal tension [T-1 ~> s-1]
    div_xx, &     ! Estimate of horizontal divergence at h-points [T-1 ~> s-1]
    sh_xx, &      ! horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
    sh_xx_bt, &   ! barotropic horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
    str_xx,&      ! str_xx is the diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2]
    str_xx_GME,&  ! smoothed diagonal term in the stress tensor from GME [H L2 T-2 ~> m3 s-2 or kg s-2]
    bhstr_xx, &   ! A copy of str_xx that only contains the biharmonic contribution [H L2 T-2 ~> m3 s-2 or kg s-2]
    FrictWorkIntz, & ! depth integrated energy dissipated by lateral friction [R L2 T-3 ~> W m-2]
    ! Leith_Kh_h, & ! Leith Laplacian viscosity at h-points [L2 T-1 ~> m2 s-1]
    ! Leith_Ah_h, & ! Leith bi-harmonic viscosity at h-points [L4 T-1 ~> m4 s-1]
    grad_vort_mag_h, & ! Magnitude of vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
    grad_vort_mag_h_2d, & ! Magnitude of 2d vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
    Del2vort_h, & ! Laplacian of vorticity at h-points [L-2 T-1 ~> m-2 s-1]
    grad_div_mag_h, &     ! Magnitude of divergence gradient at h-points [L-1 T-1 ~> m-1 s-1]
    dudx, dvdy, &    ! components in the horizontal tension [T-1 ~> s-1]
    grad_vel_mag_h, & ! Magnitude of the velocity gradient tensor squared at h-points [T-2 ~> s-2]
    grad_vel_mag_bt_h, & ! Magnitude of the barotropic velocity gradient tensor squared at h-points [T-2 ~> s-2]
    grad_d2vel_mag_h, & ! Magnitude of the Laplacian of the velocity vector, squared [L-2 T-2 ~> m-2 s-2]
    boundary_mask_h ! A mask that zeroes out cells with at least one land edge [nondim]

  real, allocatable, dimension(:,:) :: hf_diffu_2d ! Depth sum of hf_diffu [L T-2 ~> m s-2]
  real, allocatable, dimension(:,:) :: hf_diffv_2d ! Depth sum of hf_diffv [L T-2 ~> m s-2]

  real, dimension(SZIB_(G),SZJB_(G)) :: &
    dvdx, dudy, & ! components in the shearing strain [T-1 ~> s-1]
    dDel2vdx, dDel2udy, & ! Components in the biharmonic equivalent of the shearing strain [L-2 T-1 ~> m-2 s-1]
    dvdx_bt, dudy_bt,   & ! components in the barotropic shearing strain [T-1 ~> s-1]
    sh_xy,  &     ! horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
    sh_xy_bt, &   ! barotropic horizontal shearing strain (du/dy + dv/dx) inc. metric terms [T-1 ~> s-1]
    str_xy, &     ! str_xy is the cross term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2]
    str_xy_GME, & ! smoothed cross term in the stress tensor from GME [H L2 T-2 ~> m3 s-2 or kg s-2]
    bhstr_xy, &   ! A copy of str_xy that only contains the biharmonic contribution [H L2 T-2 ~> m3 s-2 or kg s-2]
    vort_xy, & ! Vertical vorticity (dv/dx - du/dy) including metric terms [T-1 ~> s-1]
    Leith_Kh_q, & ! Leith Laplacian viscosity at q-points [L2 T-1 ~> m2 s-1]
    Leith_Ah_q, & ! Leith bi-harmonic viscosity at q-points [L4 T-1 ~> m4 s-1]
    grad_vort_mag_q, & ! Magnitude of vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
    grad_vort_mag_q_2d, & ! Magnitude of 2d vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
    Del2vort_q, & ! Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
    grad_div_mag_q, &  ! Magnitude of divergence gradient at q-points [L-1 T-1 ~> m-1 s-1]
    grad_vel_mag_q, &  ! Magnitude of the velocity gradient tensor squared at q-points [T-2 ~> s-2]
    hq, &  ! harmonic mean of the harmonic means of the u- & v point thicknesses [H ~> m or kg m-2]
           ! This form guarantees that hq/hu < 4.
    grad_vel_mag_bt_q, &  ! Magnitude of the barotropic velocity gradient tensor squared at q-points [T-2 ~> s-2]
    boundary_mask_q ! A mask that zeroes out cells with at least one land edge [nondim]

  real, dimension(SZIB_(G),SZJB_(G),SZK_(G)) :: &
    Ah_q, &      ! biharmonic viscosity at corner points [L4 T-1 ~> m4 s-1]
    Kh_q, &      ! Laplacian viscosity at corner points [L2 T-1 ~> m2 s-1]
    vort_xy_q, & ! vertical vorticity at corner points [T-1 ~> s-1]
    sh_xy_q,   & ! horizontal shearing strain at corner points [T-1 ~> s-1]
    GME_coeff_q, &  !< GME coeff. at q-points [L2 T-1 ~> m2 s-1]
    max_diss_rate_q, & ! maximum possible energy dissipated by lateral friction [L2 T-3 ~> m2 s-3]
    ShSt         ! A diagnostic array of shear stress [T-1 ~> s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)+1) :: &
    KH_u_GME  !< interface height diffusivities in u-columns [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)+1) :: &
    KH_v_GME  !< interface height diffusivities in v-columns [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    Ah_h, &          ! biharmonic viscosity at thickness points [L4 T-1 ~> m4 s-1]
    Kh_h, &          ! Laplacian viscosity at thickness points [L2 T-1 ~> m2 s-1]
    max_diss_rate_h, & ! maximum possible energy dissipated by lateral friction [L2 T-3 ~> m2 s-3]
    FrictWork, &     ! work done by MKE dissipation mechanisms [R L2 T-3 ~> W m-2]
    FrictWork_GME, & ! work done by GME [R L2 T-3 ~> W m-2]
    div_xx_h,      & ! horizontal divergence [T-1 ~> s-1]
    sh_xx_h,       & ! horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
    NoSt             ! A diagnostic array of normal stress [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    grid_Re_Kh, &    !< Grid Reynolds number for Laplacian horizontal viscosity at h points [nondim]
    grid_Re_Ah, &    !< Grid Reynolds number for Biharmonic horizontal viscosity at h points [nondim]
    GME_coeff_h      !< GME coeff. at h-points [L2 T-1 ~> m2 s-1]
  real :: Ah         ! biharmonic viscosity [L4 T-1 ~> m4 s-1]
  real :: Kh         ! Laplacian  viscosity [L2 T-1 ~> m2 s-1]
  real :: AhSm       ! Smagorinsky biharmonic viscosity [L4 T-1 ~> m4 s-1]
  real :: AhLth      ! 2D Leith biharmonic viscosity [L4 T-1 ~> m4 s-1]
  real :: mod_Leith  ! nondimensional coefficient for divergence part of modified Leith
                     ! viscosity. Here set equal to nondimensional Laplacian Leith constant.
                     ! This is set equal to zero if modified Leith is not used.
  real :: Shear_mag  ! magnitude of the shear [T-1 ~> s-1]
  real :: vert_vort_mag ! magnitude of the vertical vorticity gradient [L-1 T-1 ~> m-1 s-1]
  real :: h2uq, h2vq ! temporary variables [H2 ~> m2 or kg2 m-4].
  real :: hu, hv     ! Thicknesses interpolated by arithmetic means to corner
                     ! points; these are first interpolated to u or v velocity
                     ! points where masks are applied [H ~> m or kg m-2].
  real :: h_neglect  ! thickness so small it can be lost in roundoff and so neglected [H ~> m or kg m-2]
  real :: h_neglect3 ! h_neglect^3 [H3 ~> m3 or kg3 m-6]
  real :: hrat_min   ! minimum thicknesses at the 4 neighboring
                     ! velocity points divided by the thickness at the stress
                     ! point (h or q point) [nondim]
  real :: visc_bound_rem ! fraction of overall viscous bounds that
                         ! remain to be applied [nondim]
  real :: Kh_scale  ! A factor between 0 and 1 by which the horizontal
                    ! Laplacian viscosity is rescaled [nondim]
  real :: RoScl     ! The scaling function for MEKE source term [nondim]
  real :: FatH      ! abs(f) at h-point for MEKE source term [T-1 ~> s-1]
  real :: local_strain ! Local variable for interpolating computed strain rates [T-1 ~> s-1].
  real :: meke_res_fn ! A copy of the resolution scaling factor if being applied to MEKE. Otherwise =1.
  real :: GME_coeff ! The GME (negative) viscosity coefficient [L2 T-1 ~> m2 s-1]
  real :: GME_coeff_limiter ! Maximum permitted value of the GME coefficient [L2 T-1 ~> m2 s-1]
  real :: FWfrac    ! Fraction of maximum theoretical energy transfer to use when scaling GME coefficient [nondim]
  real :: DY_dxBu   ! Ratio of meridional over zonal grid spacing at vertices [nondim]
  real :: DX_dyBu   ! Ratio of zonal over meridiononal grid spacing at vertices [nondim]
  real :: DY_dxCv   ! Ratio of meridional over zonal grid spacing at faces [nondim]
  real :: DX_dyCu   ! Ratio of zonal over meridional grid spacing at faces [nondim]
  real :: Sh_F_pow  ! The ratio of shear over the absolute value of f raised to some power and rescaled [nondim]
  real :: backscat_subround ! The ratio of f over Shear_mag that is so small that the backscatter
                    ! calculation gives the same value as if f were 0 [nondim].
  real :: H0_GME    ! Depth used to scale down GME coefficient in shallow areas [Z ~> m]
  real :: KE        ! Local kinetic energy [L2 T-2 ~> m2 s-2]

  logical :: rescale_Kh, legacy_bound
  logical :: find_FrictWork
  logical :: apply_OBC = .false.
  logical :: use_MEKE_Ku
  logical :: use_MEKE_Au
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k, n
  real :: inv_PI3, inv_PI2, inv_PI6
  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  h_neglect  = GV%H_subroundoff
  h_neglect3 = h_neglect**3
  inv_PI3 = 1.0/((4.0*atan(1.0))**3)
  inv_PI2 = 1.0/((4.0*atan(1.0))**2)
  inv_PI6 = inv_PI3 * inv_PI3

  Ah_h(:,:,:) = 0.0
  Kh_h(:,:,:) = 0.0

  if (present(OBC)) then ; if (associated(OBC)) then ; if (OBC%OBC_pe) then
    apply_OBC = OBC%Flather_u_BCs_exist_globally .or. OBC%Flather_v_BCs_exist_globally
    apply_OBC = .true.
  endif ; endif ; endif

  if (.not.associated(CS)) call MOM_error(FATAL, &
         "MOM_hor_visc: Module must be initialized before it is used.")
  if (.not.(CS%Laplacian .or. CS%biharmonic)) return

  find_FrictWork = (CS%id_FrictWork > 0)
  if (CS%id_FrictWorkIntz > 0) find_FrictWork = .true.
  if (associated(MEKE)) then
    if (associated(MEKE%mom_src)) find_FrictWork = .true.
    backscat_subround = 0.0
    if (find_FrictWork .and. associated(MEKE%mom_src) .and. (MEKE%backscatter_Ro_c > 0.0) .and. &
        (MEKE%backscatter_Ro_Pow /= 0.0)) &
      backscat_subround = (1.0e-16/MEKE%backscatter_Ro_c)**(1.0/MEKE%backscatter_Ro_Pow)
  endif

  rescale_Kh = .false.
  if (associated(VarMix)) then
    rescale_Kh = VarMix%Resoln_scaled_Kh
    if (rescale_Kh .and. &
    (.not.associated(VarMix%Res_fn_h) .or. .not.associated(VarMix%Res_fn_q))) &
      call MOM_error(FATAL, "MOM_hor_visc: VarMix%Res_fn_h and " //&
        "VarMix%Res_fn_q both need to be associated with Resoln_scaled_Kh.")
  endif
  legacy_bound = (CS%Smagorinsky_Kh .or. CS%Leith_Kh) .and. &
                 (CS%bound_Kh .and. .not.CS%better_bound_Kh)

  ! Toggle whether to use a Laplacian viscosity derived from MEKE
  use_MEKE_Ku = associated(MEKE%Ku)
  use_MEKE_Au = associated(MEKE%Au)

  if (CS%use_GME) then
    do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
      boundary_mask_h(i,j) = (G%mask2dCu(I,j) * G%mask2dCv(i,J) * G%mask2dCu(I-1,j) * G%mask2dCv(i,J-1))
    enddo ; enddo

    do J=js-2,Jeq+1 ; do I=is-2,Ieq+1
      boundary_mask_q(I,J) = (G%mask2dCv(i,J) * G%mask2dCv(i+1,J) * G%mask2dCu(I,j) * G%mask2dCu(I,j-1))
    enddo; enddo

    ! initialize diag. array with zeros
    GME_coeff_h(:,:,:) = 0.0
    GME_coeff_q(:,:,:) = 0.0
    str_xx_GME(:,:) = 0.0
    str_xy_GME(:,:) = 0.0

    ! Get barotropic velocities and their gradients
    call barotropic_get_tav(BT, ubtav, vbtav, G, US)
    call pass_vector(ubtav, vbtav, G%Domain)

    do j=js-1,je+2 ; do i=is-1,ie+2
      dudx_bt(i,j) = CS%DY_dxT(i,j)*(G%IdyCu(I,j) * ubtav(I,j) - &
                                     G%IdyCu(I-1,j) * ubtav(I-1,j))
      dvdy_bt(i,j) = CS%DX_dyT(i,j)*(G%IdxCv(i,J) * vbtav(i,J) - &
                                     G%IdxCv(i,J-1) * vbtav(i,J-1))
    enddo; enddo

    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      sh_xx_bt(i,j) = dudx_bt(i,j) - dvdy_bt(i,j)
    enddo ; enddo

    ! Components for the barotropic shearing strain
    do J=js-2,Jeq+1 ; do I=is-2,Ieq+1
      dvdx_bt(I,J) = CS%DY_dxBu(I,J)*(vbtav(i+1,J)*G%IdyCv(i+1,J) &
                                    - vbtav(i,J)*G%IdyCv(i,J))
      dudy_bt(I,J) = CS%DX_dyBu(I,J)*(ubtav(I,j+1)*G%IdxCu(I,j+1) &
                                    - ubtav(I,j)*G%IdxCu(I,j))
    enddo ; enddo

    call pass_vector(dudx_bt, dvdy_bt, G%Domain, stagger=BGRID_NE)
    call pass_vector(dvdx_bt, dudy_bt, G%Domain, stagger=AGRID)

    if (CS%no_slip) then
      do J=js-1,Jeq ; do I=is-1,Ieq
        sh_xy_bt(I,J) = (2.0-G%mask2dBu(I,J)) * ( dvdx_bt(I,J) + dudy_bt(I,J) )
      enddo ; enddo
    else
      do J=js-1,Jeq ; do I=is-1,Ieq
        sh_xy_bt(I,J) = G%mask2dBu(I,J) * ( dvdx_bt(I,J) + dudy_bt(I,J) )
      enddo ; enddo
    endif

    do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
      grad_vel_mag_bt_h(i,j) = boundary_mask_h(i,j) * (dudx_bt(i,j)**2 + dvdy_bt(i,j)**2 + &
            (0.25*((dvdx_bt(I,J)+dvdx_bt(I-1,J-1))+(dvdx_bt(I,J-1)+dvdx_bt(I-1,J))))**2 + &
            (0.25*((dudy_bt(I,J)+dudy_bt(I-1,J-1))+(dudy_bt(I,J-1)+dudy_bt(I-1,J))))**2)
    enddo ; enddo

    do J=js-2,Jeq+1 ; do I=is-2,Ieq+1
      grad_vel_mag_bt_q(I,J) = boundary_mask_q(I,J) * (dvdx_bt(i,j)**2 + dudy_bt(i,j)**2 + &
            (0.25*((dudx_bt(i,j)+dudx_bt(i+1,j+1))+(dudx_bt(i,j+1)+dudx_bt(i+1,j))))**2 + &
            (0.25*((dvdy_bt(i,j)+dvdy_bt(i+1,j+1))+(dvdy_bt(i,j+1)+dvdy_bt(i+1,j))))**2)
    enddo ; enddo

  endif ! use_GME

  !$OMP parallel do default(none) &
  !$OMP shared( &
  !$OMP   CS, G, GV, US, OBC, VarMix, MEKE, u, v, h, &
  !$OMP   is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, &
  !$OMP   apply_OBC, rescale_Kh, legacy_bound, find_FrictWork, &
  !$OMP   use_MEKE_Ku, use_MEKE_Au, boundary_mask_h, boundary_mask_q, &
  !$OMP   backscat_subround, GME_coeff_limiter, &
  !$OMP   h_neglect, h_neglect3, FWfrac, inv_PI3, inv_PI6, H0_GME, &
  !$OMP   diffu, diffv, max_diss_rate_h, max_diss_rate_q, &
  !$OMP   Kh_h, Kh_q, Ah_h, Ah_q, FrictWork, FrictWork_GME, &
  !$OMP   div_xx_h, sh_xx_h, vort_xy_q, sh_xy_q, GME_coeff_h, GME_coeff_q, &
  !$OMP   TD, KH_u_GME, KH_v_GME, grid_Re_Kh, grid_Re_Ah, NoSt, ShSt &
  !$OMP ) &
  !$OMP private( &
  !$OMP   i, j, k, n, &
  !$OMP   dudx, dudy, dvdx, dvdy, sh_xx, sh_xy, h_u, h_v, &
  !$OMP   Del2u, Del2v, DY_dxBu, DX_dyBu, sh_xx_bt, sh_xy_bt, &
  !$OMP   str_xx, str_xy, bhstr_xx, bhstr_xy, str_xx_GME, str_xy_GME, &
  !$OMP   vort_xy, vort_xy_dx, vort_xy_dy, div_xx, div_xx_dx, div_xx_dy, &
  !$OMP   grad_div_mag_h, grad_div_mag_q, grad_vort_mag_h, grad_vort_mag_q, &
  !$OMP   grad_vort_mag_h_2d, grad_vort_mag_q_2d, &
  !$OMP   grad_vel_mag_h, grad_vel_mag_q, &
  !$OMP   grad_vel_mag_bt_h, grad_vel_mag_bt_q, grad_d2vel_mag_h, &
  !$OMP   meke_res_fn, Shear_mag, vert_vort_mag, hrat_min, visc_bound_rem, &
  !$OMP   Kh, Ah, AhSm, AhLth, local_strain, Sh_F_pow, &
  !$OMP   dDel2vdx, dDel2udy, DY_dxCv, DX_dyCu, Del2vort_q, Del2vort_h, KE, &
  !$OMP   h2uq, h2vq, hu, hv, hq, FatH, RoScl, GME_coeff &
  !$OMP )
  do k=1,nz

    ! The following are the forms of the horizontal tension and horizontal
    ! shearing strain advocated by Smagorinsky (1993) and discussed in
    ! Griffies and Hallberg (2000).

    ! Calculate horizontal tension
    do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
      dudx(i,j) = CS%DY_dxT(i,j)*(G%IdyCu(I,j) * u(I,j,k) - &
                                  G%IdyCu(I-1,j) * u(I-1,j,k))
      dvdy(i,j) = CS%DX_dyT(i,j)*(G%IdxCv(i,J) * v(i,J,k) - &
                                  G%IdxCv(i,J-1) * v(i,J-1,k))
      sh_xx(i,j) = dudx(i,j) - dvdy(i,j)
      if (CS%id_normstress > 0) NoSt(i,j,k) = sh_xx(i,j)
    enddo ; enddo

    ! Components for the shearing strain
    do J=Jsq-2,Jeq+2 ; do I=Isq-2,Ieq+2
      dvdx(I,J) = CS%DY_dxBu(I,J)*(v(i+1,J,k)*G%IdyCv(i+1,J) - v(i,J,k)*G%IdyCv(i,J))
      dudy(I,J) = CS%DX_dyBu(I,J)*(u(I,j+1,k)*G%IdxCu(I,j+1) - u(I,j,k)*G%IdxCu(I,j))
    enddo ; enddo

    ! Interpolate the thicknesses to velocity points.
    ! The extra wide halos are to accommodate the cross-corner-point projections
    ! in OBCs, which are not ordinarily be necessary, and might not be necessary
    ! even with OBCs if the accelerations are zeroed at OBC points, in which
    ! case the j-loop for h_u could collapse to j=js=1,je+1. -RWH
    if (CS%use_land_mask) then
      do j=js-2,je+2 ; do I=Isq-1,Ieq+1
        h_u(I,j) = 0.5 * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i+1,j)*h(i+1,j,k))
      enddo ; enddo
      do J=Jsq-1,Jeq+1 ; do i=is-2,ie+2
        h_v(i,J) = 0.5 * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i,j+1)*h(i,j+1,k))
      enddo ; enddo
    else
      do j=js-2,je+2 ; do I=Isq-1,Ieq+1
        h_u(I,j) = 0.5 * (h(i,j,k) + h(i+1,j,k))
      enddo ; enddo
      do J=Jsq-1,Jeq+1 ; do i=is-2,ie+2
        h_v(i,J) = 0.5 * (h(i,j,k) + h(i,j+1,k))
      enddo ; enddo
    endif

    ! Adjust contributions to shearing strain and interpolated values of
    ! thicknesses on open boundaries.
    if (apply_OBC) then ; do n=1,OBC%number_of_segments
      J = OBC%segment(n)%HI%JsdB ; I = OBC%segment(n)%HI%IsdB
      if (OBC%zero_strain .or. OBC%freeslip_strain .or. OBC%computed_strain) then
        if (OBC%segment(n)%is_N_or_S .and. (J >= js-2) .and. (J <= Jeq+1)) then
          do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
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
          enddo
        elseif (OBC%segment(n)%is_E_or_W .and. (I >= is-2) .and. (I <= Ieq+1)) then
          do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
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
          enddo
        endif
      endif


      if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
        ! There are extra wide halos here to accommodate the cross-corner-point
        ! OBC projections, but they might not be necessary if the accelerations
        ! are always zeroed out at OBC points, in which case the i-loop below
        ! becomes do i=is-1,ie+1. -RWH
        if ((J >= Jsq-1) .and. (J <= Jeq+1)) then
          do i = max(is-2,OBC%segment(n)%HI%isd), min(ie+2,OBC%segment(n)%HI%ied)
            h_v(i,J) = h(i,j,k)
          enddo
        endif
      elseif (OBC%segment(n)%direction == OBC_DIRECTION_S) then
        if ((J >= Jsq-1) .and. (J <= Jeq+1)) then
          do i = max(is-2,OBC%segment(n)%HI%isd), min(ie+2,OBC%segment(n)%HI%ied)
            h_v(i,J) = h(i,j+1,k)
          enddo
        endif
      elseif (OBC%segment(n)%direction == OBC_DIRECTION_E) then
        if ((I >= Isq-1) .and. (I <= Ieq+1)) then
          do j = max(js-2,OBC%segment(n)%HI%jsd), min(je+2,OBC%segment(n)%HI%jed)
            h_u(I,j) = h(i,j,k)
          enddo
        endif
      elseif (OBC%segment(n)%direction == OBC_DIRECTION_W) then
        if ((I >= Isq-1) .and. (I <= Ieq+1)) then
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
          do I = max(Isq-1,OBC%segment(n)%HI%IsdB), min(Ieq+1,OBC%segment(n)%HI%IedB)
            h_u(I,j+1) = h_u(I,j)
          enddo
        endif
      elseif (OBC%segment(n)%direction == OBC_DIRECTION_S) then
        if ((J >= js-1) .and. (J <= je+1)) then
          do I = max(Isq-1,OBC%segment(n)%HI%isd), min(Ieq+1,OBC%segment(n)%HI%ied)
            h_u(I,j) = h_u(I,j+1)
          enddo
        endif
      elseif (OBC%segment(n)%direction == OBC_DIRECTION_E) then
        if ((I >= is-2) .and. (I <= ie)) then
          do J = max(Jsq-1,OBC%segment(n)%HI%jsd), min(Jeq+1,OBC%segment(n)%HI%jed)
            h_v(i+1,J) = h_v(i,J)
          enddo
        endif
      elseif (OBC%segment(n)%direction == OBC_DIRECTION_W) then
        if ((I >= is-1) .and. (I <= ie+1)) then
          do J = max(Jsq-1,OBC%segment(n)%HI%jsd), min(Jeq+1,OBC%segment(n)%HI%jed)
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
      if (apply_OBC) then; if (OBC%zero_biharmonic) then
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
      endif; endif
    endif

    ! Vorticity
    if (CS%no_slip) then
      do J=Jsq-2,Jeq+2 ; do I=Isq-2,Ieq+2
        vort_xy(I,J) = (2.0-G%mask2dBu(I,J)) * ( dvdx(I,J) - dudy(I,J) )
      enddo ; enddo
    else
      do J=Jsq-2,Jeq+2 ; do I=Isq-2,Ieq+2
        vort_xy(I,J) = G%mask2dBu(I,J) * ( dvdx(I,J) - dudy(I,J) )
      enddo ; enddo
    endif

    ! Divergence
    do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
      div_xx(i,j) = dudx(i,j) + dvdy(i,j)
    enddo ; enddo

    if ((CS%Leith_Kh) .or. (CS%Leith_Ah)) then

      ! Vorticity gradient
      do J=Jsq-1,Jeq+1 ; do i=Isq-1,Ieq+2
        DY_dxBu = G%dyBu(I,J) * G%IdxBu(I,J)
        vort_xy_dx(i,J) = DY_dxBu * (vort_xy(I,J) * G%IdyCu(I,j) - vort_xy(I-1,J) * G%IdyCu(I-1,j))
      enddo ; enddo

      do j=Jsq-1,Jeq+2 ; do I=Isq-1,Ieq+1
        DX_dyBu = G%dxBu(I,J) * G%IdyBu(I,J)
        vort_xy_dy(I,j) = DX_dyBu * (vort_xy(I,J) * G%IdxCv(i,J) - vort_xy(I,J-1) * G%IdxCv(i,J-1))
      enddo ; enddo

      ! Laplacian of vorticity
      do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
        DY_dxBu = G%dyBu(I,J) * G%IdxBu(I,J)
        DX_dyBu = G%dxBu(I,J) * G%IdyBu(I,J)

        Del2vort_q(I,J) = DY_dxBu * (vort_xy_dx(i+1,J) * G%IdyCv(i+1,J) - vort_xy_dx(i,J) * G%IdyCv(i,J)) + &
                        DX_dyBu * (vort_xy_dy(I,j+1) * G%IdyCu(I,j+1) - vort_xy_dy(I,j) * G%IdyCu(I,j))
      enddo ; enddo
      do J=Jsq,Jeq+1 ; do I=Isq,Ieq+1
        Del2vort_h(i,j) = 0.25*(Del2vort_q(I,J) + Del2vort_q(I-1,J) + Del2vort_q(I,J-1) + Del2vort_q(I-1,J-1))
      enddo ; enddo

      if (CS%modified_Leith) then

        ! Divergence gradient
        do j=Jsq-1,Jeq+2 ; do I=Isq-1,Ieq+1
          div_xx_dx(I,j) = G%IdxCu(I,j)*(div_xx(i+1,j) - div_xx(i,j))
        enddo ; enddo
        do J=Jsq-1,Jeq+1 ; do i=Isq-1,Ieq+2
          div_xx_dy(i,J) = G%IdyCv(i,J)*(div_xx(i,j+1) - div_xx(i,j))
        enddo ; enddo

        ! Magnitude of divergence gradient
        do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
          grad_div_mag_h(i,j) =sqrt((0.5*(div_xx_dx(I,j) + div_xx_dx(I-1,j)))**2 + &
          (0.5 * (div_xx_dy(i,J) + div_xx_dy(i,J-1)))**2)
        enddo ; enddo
        !do J=js-1,Jeq ; do I=is-1,Ieq
        do j=Jsq-1,Jeq+1 ; do i=Isq-1,Ieq+1
          grad_div_mag_q(I,J) =sqrt((0.5*(div_xx_dx(I,j) + div_xx_dx(I,j+1)))**2 + &
          (0.5 * (div_xx_dy(i,J) + div_xx_dy(i+1,J)))**2)
        enddo ; enddo

      else

        do j=Jsq-1,Jeq+2 ; do I=is-2,Ieq+1
          div_xx_dx(I,j) = 0.0
        enddo ; enddo
        do J=Jsq-1,Jeq+1 ; do i=Isq-1,Ieq+2
          div_xx_dy(i,J) = 0.0
        enddo ; enddo
        do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
          grad_div_mag_h(i,j) = 0.0
        enddo ; enddo
        do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
          grad_div_mag_q(I,J) = 0.0
        enddo ; enddo

      endif ! CS%modified_Leith

      ! Add in beta for the Leith viscosity
      if (CS%use_beta_in_Leith) then
        do J=js-2,Jeq+1 ; do i=is-1,Ieq+1
          vort_xy_dx(i,J) = vort_xy_dx(i,J) + 0.5 * ( G%dF_dx(i,j) + G%dF_dx(i,j+1))
        enddo ; enddo
        do j=js-1,Jeq+1 ; do I=is-2,Ieq+1
          vort_xy_dy(I,j) = vort_xy_dy(I,j) + 0.5 * ( G%dF_dy(i,j) + G%dF_dy(i+1,j))
        enddo ; enddo
      endif ! CS%use_beta_in_Leith

      if (CS%use_QG_Leith_visc) then

        do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
          grad_vort_mag_h_2d(i,j) = SQRT((0.5*(vort_xy_dx(i,J) + vort_xy_dx(i,J-1)))**2 + &
                                         (0.5*(vort_xy_dy(I,j) + vort_xy_dy(I-1,j)))**2 )
        enddo ; enddo
        do J=js-1,Jeq ; do I=is-1,Ieq
          grad_vort_mag_q_2d(I,J) = SQRT((0.5*(vort_xy_dx(i,J) + vort_xy_dx(i+1,J)))**2 + &
                                         (0.5*(vort_xy_dy(I,j) + vort_xy_dy(I,j+1)))**2 )
        enddo ; enddo

        ! This accumulates terms, some of which are in VarMix, so rescaling can not be done here.
        call calc_QG_Leith_viscosity(VarMix, G, GV, US, h, k, div_xx_dx, div_xx_dy, &
                                     vort_xy_dx, vort_xy_dy)

      endif

      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        grad_vort_mag_h(i,j) = SQRT((0.5*(vort_xy_dx(i,J) + vort_xy_dx(i,J-1)))**2 + &
                                    (0.5*(vort_xy_dy(I,j) + vort_xy_dy(I-1,j)))**2 )
      enddo ; enddo
      do J=js-1,Jeq ; do I=is-1,Ieq
        grad_vort_mag_q(I,J) = SQRT((0.5*(vort_xy_dx(i,J) + vort_xy_dx(i+1,J)))**2 + &
                                    (0.5*(vort_xy_dy(I,j) + vort_xy_dy(I,j+1)))**2 )
      enddo ; enddo

    endif ! CS%Leith_Kh

    meke_res_fn = 1.

    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      if ((CS%Smagorinsky_Kh) .or. (CS%Smagorinsky_Ah)) then
        Shear_mag = sqrt(sh_xx(i,j)*sh_xx(i,j) + &
          0.25*((sh_xy(I-1,J-1)*sh_xy(I-1,J-1) + sh_xy(I,J)*sh_xy(I,J)) + &
                (sh_xy(I-1,J)*sh_xy(I-1,J) + sh_xy(I,J-1)*sh_xy(I,J-1))))
      endif
      if ((CS%Leith_Kh) .or. (CS%Leith_Ah)) then
        if (CS%use_QG_Leith_visc) then
          vert_vort_mag = MIN(grad_vort_mag_h(i,j) + grad_div_mag_h(i,j),3.*grad_vort_mag_h_2d(i,j))
        else
          vert_vort_mag = (grad_vort_mag_h(i,j) + grad_div_mag_h(i,j))
        endif
      endif
      if (CS%better_bound_Ah .or. CS%better_bound_Kh) then
        hrat_min = min(1.0, min(h_u(I,j), h_u(I-1,j), h_v(i,J), h_v(i,J-1)) / &
                            (h(i,j,k) + h_neglect) )
        visc_bound_rem = 1.0
      endif

      if (CS%Laplacian) then
        ! Determine the Laplacian viscosity at h points, using the
        ! largest value from several parameterizations.
        Kh = CS%Kh_bg_xx(i,j) ! Static (pre-computed) background viscosity
        if (CS%add_LES_viscosity) then
          if (CS%Smagorinsky_Kh) Kh = Kh + CS%Laplac2_const_xx(i,j) * Shear_mag
          if (CS%Leith_Kh) Kh = Kh + CS%Laplac3_const_xx(i,j) * vert_vort_mag*inv_PI3
        else
          if (CS%Smagorinsky_Kh) Kh = max( Kh, CS%Laplac2_const_xx(i,j) * Shear_mag )
          if (CS%Leith_Kh) Kh = max( Kh, CS%Laplac3_const_xx(i,j) * vert_vort_mag*inv_PI3)
        endif
        ! All viscosity contributions above are subject to resolution scaling
        if (rescale_Kh) Kh = VarMix%Res_fn_h(i,j) * Kh
        if (CS%res_scale_MEKE) meke_res_fn = VarMix%Res_fn_h(i,j)
        ! Older method of bounding for stability
        if (legacy_bound) Kh = min(Kh, CS%Kh_Max_xx(i,j))
        Kh = max( Kh, CS%Kh_bg_min ) ! Place a floor on the viscosity, if desired.
        if (use_MEKE_Ku) &
          Kh = Kh + MEKE%Ku(i,j) * meke_res_fn ! *Add* the MEKE contribution (might be negative)
        if (CS%anisotropic) Kh = Kh + CS%Kh_aniso * ( 1. - CS%n1n2_h(i,j)**2 ) ! *Add* the tension component
                                                                               ! of anisotropic viscosity

        ! Newer method of bounding for stability
        if (CS%better_bound_Kh) then
          if (Kh >= hrat_min*CS%Kh_Max_xx(i,j)) then
            visc_bound_rem = 0.0
            Kh = hrat_min*CS%Kh_Max_xx(i,j)
          else
            visc_bound_rem = 1.0 - Kh / (hrat_min*CS%Kh_Max_xx(i,j))
          endif
        endif

        if ((CS%id_Kh_h>0) .or. find_FrictWork .or. CS%debug) Kh_h(i,j,k) = Kh

        if (CS%id_grid_Re_Kh>0) then
          KE = 0.125*((u(I,j,k)+u(I-1,j,k))**2 + (v(i,J,k)+v(i,J-1,k))**2)
          grid_Re_Kh(i,j,k) = (sqrt(KE) * sqrt(CS%grid_sp_h2(i,j))) &
              / max(Kh, CS%min_grid_Kh)
        endif

        if (CS%id_div_xx_h>0) div_xx_h(i,j,k) = div_xx(i,j)
        if (CS%id_sh_xx_h>0) sh_xx_h(i,j,k) = sh_xx(i,j)

        str_xx(i,j) = -Kh * sh_xx(i,j)
      else   ! not Laplacian
        str_xx(i,j) = 0.0
      endif ! Laplacian

      if (CS%anisotropic) then
        ! Shearing-strain averaged to h-points
        local_strain = 0.25 * ( (sh_xy(I,J) + sh_xy(I-1,J-1)) + (sh_xy(I-1,J) + sh_xy(I,J-1)) )
        ! *Add* the shear-strain contribution to the xx-component of stress
        str_xx(i,j) = str_xx(i,j) - CS%Kh_aniso * CS%n1n2_h(i,j) * CS%n1n1_m_n2n2_h(i,j) * local_strain
      endif

      if (CS%biharmonic) then
        ! Determine the biharmonic viscosity at h points, using the
        ! largest value from several parameterizations.
        AhSm = 0.0; AhLth = 0.0
        if ((CS%Smagorinsky_Ah) .or. (CS%Leith_Ah)) then
          if (CS%Smagorinsky_Ah) then
            if (CS%bound_Coriolis) then
              AhSm = Shear_mag * (CS%Biharm_const_xx(i,j) + &
                                  CS%Biharm_const2_xx(i,j)*Shear_mag)
            else
              AhSm = CS%Biharm_const_xx(i,j) * Shear_mag
            endif
          endif
          if (CS%Leith_Ah) AhLth = CS%Biharm6_const_xx(i,j) * abs(Del2vort_h(i,j)) * inv_PI6
          Ah = MAX(MAX(CS%Ah_bg_xx(i,j), AhSm), AhLth)
          if (CS%bound_Ah .and. .not.CS%better_bound_Ah) &
            Ah = MIN(Ah, CS%Ah_Max_xx(i,j))
        else
          Ah = CS%Ah_bg_xx(i,j)
        endif ! Smagorinsky_Ah or Leith_Ah

        if (use_MEKE_Au) Ah = Ah + MEKE%Au(i,j) ! *Add* the MEKE contribution

        if (CS%Re_Ah > 0.0) then
          KE = 0.125*((u(I,j,k)+u(I-1,j,k))**2 + (v(i,J,k)+v(i,J-1,k))**2)
          Ah = sqrt(KE) * CS%Re_Ah_const_xx(i,j)
        endif

        if (CS%better_bound_Ah) then
          Ah = MIN(Ah, visc_bound_rem*hrat_min*CS%Ah_Max_xx(i,j))
        endif

        if ((CS%id_Ah_h>0) .or. find_FrictWork .or. CS%debug) Ah_h(i,j,k) = Ah

        if (CS%id_grid_Re_Ah>0) then
          KE = 0.125*((u(I,j,k)+u(I-1,j,k))**2 + (v(i,J,k)+v(i,J-1,k))**2)
          grid_Re_Ah(i,j,k) = (sqrt(KE) * CS%grid_sp_h3(i,j)) &
              / max(Ah, CS%min_grid_Ah)
        endif

        str_xx(i,j) = str_xx(i,j) + Ah * &
          (CS%DY_dxT(i,j) * (G%IdyCu(I,j)*Del2u(I,j) - G%IdyCu(I-1,j)*Del2u(I-1,j)) - &
           CS%DX_dyT(i,j) * (G%IdxCv(i,J)*Del2v(i,J) - G%IdxCv(i,J-1)*Del2v(i,J-1)))

        ! Keep a copy of the biharmonic contribution for backscatter parameterization
        bhstr_xx(i,j) = Ah * &
          (CS%DY_dxT(i,j) * (G%IdyCu(I,j)*Del2u(I,j) - G%IdyCu(I-1,j)*Del2u(I-1,j)) - &
           CS%DX_dyT(i,j) * (G%IdxCv(i,J)*Del2v(i,J) - G%IdxCv(i,J-1)*Del2v(i,J-1)))
        bhstr_xx(i,j) = bhstr_xx(i,j) * (h(i,j,k) * CS%reduction_xx(i,j))

      endif  ! biharmonic

    enddo ; enddo

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

    do J=js-1,Jeq ; do I=is-1,Ieq
      if ((CS%Smagorinsky_Kh) .or. (CS%Smagorinsky_Ah)) then
        Shear_mag = sqrt(sh_xy(I,J)*sh_xy(I,J) + &
            0.25*((sh_xx(i,j)*sh_xx(i,j) + sh_xx(i+1,j+1)*sh_xx(i+1,j+1)) + &
                  (sh_xx(i,j+1)*sh_xx(i,j+1) + sh_xx(i+1,j)*sh_xx(i+1,j))))
      endif
      if ((CS%Leith_Kh) .or. (CS%Leith_Ah)) then
        if (CS%use_QG_Leith_visc) then
          vert_vort_mag = MIN(grad_vort_mag_q(I,J) + grad_div_mag_q(I,J), 3.*grad_vort_mag_q_2d(I,J))
        else
          vert_vort_mag = (grad_vort_mag_q(I,J) + grad_div_mag_q(I,J))
        endif
      endif
      h2uq = 4.0 * h_u(I,j) * h_u(I,j+1)
      h2vq = 4.0 * h_v(i,J) * h_v(i+1,J)
      hq(I,J) = 2.0 * h2uq * h2vq / (h_neglect3 + (h2uq + h2vq) * &
              ((h_u(I,j) + h_u(I,j+1)) + (h_v(i,J) + h_v(i+1,J))))

      if (CS%better_bound_Ah .or. CS%better_bound_Kh) then
        hrat_min = min(1.0, min(h_u(I,j), h_u(I,j+1), h_v(i,J), h_v(i+1,J)) / &
                            (hq(I,J) + h_neglect) )
        visc_bound_rem = 1.0
      endif

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
            hrat_min = 1.0
          else
            ! Both hu and hv are nonzero, so take the harmonic mean.
            hq(I,J) = 2.0 * (hu * hv) / ((hu + hv) + h_neglect)
            hrat_min = min(1.0, min(hu, hv) / (hq(I,J) + h_neglect) )
          endif
        endif
      endif

      if (CS%Laplacian) then
        ! Determine the Laplacian viscosity at q points, using the
        ! largest value from several parameterizations.
        Kh = CS%Kh_bg_xy(i,j) ! Static (pre-computed) background viscosity
        if (CS%add_LES_viscosity) then
          if (CS%Smagorinsky_Kh) Kh = Kh + CS%Laplac2_const_xx(i,j) * Shear_mag
          if (CS%Leith_Kh) Kh = Kh + CS%Laplac3_const_xx(i,j) * vert_vort_mag*inv_PI3
        else
          if (CS%Smagorinsky_Kh) Kh = max( Kh, CS%Laplac2_const_xy(I,J) * Shear_mag )
          if (CS%Leith_Kh) Kh = max( Kh, CS%Laplac3_const_xy(I,J) * vert_vort_mag*inv_PI3)
        endif
        ! All viscosity contributions above are subject to resolution scaling
        if (rescale_Kh) Kh = VarMix%Res_fn_q(i,j) * Kh
        if (CS%res_scale_MEKE) meke_res_fn = VarMix%Res_fn_q(i,j)
        ! Older method of bounding for stability
        if (legacy_bound) Kh = min(Kh, CS%Kh_Max_xy(i,j))
        Kh = max( Kh, CS%Kh_bg_min ) ! Place a floor on the viscosity, if desired.
        if (use_MEKE_Ku) then ! *Add* the MEKE contribution (might be negative)
          Kh = Kh + 0.25*( (MEKE%Ku(i,j) + MEKE%Ku(i+1,j+1)) + &
                           (MEKE%Ku(i+1,j) + MEKE%Ku(i,j+1)) ) * meke_res_fn
        endif
        ! Older method of bounding for stability
        if (CS%anisotropic) Kh = Kh + CS%Kh_aniso * CS%n1n2_q(I,J)**2 ! *Add* the shear component
                                                                      ! of anisotropic viscosity

        ! Newer method of bounding for stability
        if (CS%better_bound_Kh) then
          if (Kh >= hrat_min*CS%Kh_Max_xy(I,J)) then
            visc_bound_rem = 0.0
            Kh = hrat_min*CS%Kh_Max_xy(I,J)
          elseif (CS%Kh_Max_xy(I,J)>0.) then
            visc_bound_rem = 1.0 - Kh / (hrat_min*CS%Kh_Max_xy(I,J))
          endif
        endif

        if (CS%id_Kh_q>0 .or. CS%debug) Kh_q(I,J,k) = Kh
        if (CS%id_vort_xy_q>0) vort_xy_q(I,J,k) = vort_xy(I,J)
        if (CS%id_sh_xy_q>0) sh_xy_q(I,J,k) = sh_xy(I,J)

        str_xy(I,J) = -Kh * sh_xy(I,J)
      else   ! not Laplacian
        str_xy(I,J) = 0.0
      endif ! Laplacian

      if (CS%anisotropic) then
        ! Horizontal-tension averaged to q-points
        local_strain = 0.25 * ( (sh_xx(i,j) + sh_xx(i+1,j+1)) + (sh_xx(i+1,j) + sh_xx(i,j+1)) )
        ! *Add* the tension contribution to the xy-component of stress
        str_xy(I,J) = str_xy(I,J) - CS%Kh_aniso * CS%n1n2_q(i,j) * CS%n1n1_m_n2n2_q(i,j) * local_strain
      endif

      if (CS%biharmonic) then
      ! Determine the biharmonic viscosity at q points, using the
      ! largest value from several parameterizations.
        AhSm = 0.0 ; AhLth = 0.0
        if (CS%Smagorinsky_Ah .or. CS%Leith_Ah) then
          if (CS%Smagorinsky_Ah) then
            if (CS%bound_Coriolis) then
              AhSm = Shear_mag * (CS%Biharm_const_xy(I,J) + &
                                  CS%Biharm_const2_xy(I,J)*Shear_mag)
            else
              AhSm = CS%Biharm_const_xy(I,J) * Shear_mag
            endif
          endif
          if (CS%Leith_Ah) AhLth = CS%Biharm6_const_xy(I,J) * abs(Del2vort_q(I,J)) * inv_PI6
          Ah = MAX(MAX(CS%Ah_bg_xy(I,J), AhSm), AhLth)
          if (CS%bound_Ah .and. .not.CS%better_bound_Ah) &
            Ah = MIN(Ah, CS%Ah_Max_xy(I,J))
        else
          Ah = CS%Ah_bg_xy(I,J)
        endif ! Smagorinsky_Ah or Leith_Ah

        if (use_MEKE_Au) then ! *Add* the MEKE contribution
          Ah = Ah + 0.25*( (MEKE%Au(i,j) + MEKE%Au(i+1,j+1)) +  &
                           (MEKE%Au(i+1,j) + MEKE%Au(i,j+1)) )
        endif

        if (CS%Re_Ah > 0.0) then
          KE = 0.125*((u(I,j,k)+u(I,j+1,k))**2 + (v(i,J,k)+v(i+1,J,k))**2)
          Ah = sqrt(KE) * CS%Re_Ah_const_xy(i,j)
        endif

        if (CS%better_bound_Ah) then
          Ah = MIN(Ah, visc_bound_rem*hrat_min*CS%Ah_Max_xy(I,J))
        endif

        if (CS%id_Ah_q>0 .or. CS%debug) Ah_q(I,J,k) = Ah

        str_xy(I,J) = str_xy(I,J) + Ah * ( dDel2vdx(I,J) + dDel2udy(I,J) )

        ! Keep a copy of the biharmonic contribution for backscatter parameterization
        bhstr_xy(I,J) = Ah * ( dDel2vdx(I,J) + dDel2udy(I,J) ) * &
                        (hq(I,J) * G%mask2dBu(I,J) * CS%reduction_xy(I,J))

      endif  ! biharmonic

    enddo ; enddo

    if (CS%use_GME) then
      call thickness_diffuse_get_KH(TD, KH_u_GME, KH_v_GME, G)
      call pass_vector(KH_u_GME, KH_v_GME, G%Domain)

      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        if (grad_vel_mag_bt_h(i,j)>0) then
          GME_coeff = CS%GME_efficiency * (MIN(G%bathyT(i,j)/CS%GME_h0,1.0)**2) * &
            (0.25*(KH_u_GME(I,j,k)+KH_u_GME(I-1,j,k)+KH_v_GME(i,J,k)+KH_v_GME(i,J-1,k)))
        else
          GME_coeff = 0.0
        endif

        ! apply mask
        GME_coeff = GME_coeff * boundary_mask_h(i,j)
        GME_coeff = MIN(GME_coeff, CS%GME_limiter)

        if ((CS%id_GME_coeff_h>0) .or. find_FrictWork) GME_coeff_h(i,j,k) = GME_coeff
        str_xx_GME(i,j) = GME_coeff * sh_xx_bt(i,j)

      enddo ; enddo

      do J=js-1,Jeq ; do I=is-1,Ieq
        if (grad_vel_mag_bt_q(i,j)>0) then
          GME_coeff = CS%GME_efficiency * (MIN(G%bathyT(i,j)/CS%GME_h0,1.0)**2) * &
              (0.25*(KH_u_GME(I,j,k)+KH_u_GME(I,j+1,k)+KH_v_GME(i,J,k)+KH_v_GME(i+1,J,k)))
        else
          GME_coeff = 0.0
        endif

        ! apply mask
        GME_coeff = GME_coeff * boundary_mask_q(I,J)
        GME_coeff = MIN(GME_coeff, CS%GME_limiter)

        if (CS%id_GME_coeff_q>0) GME_coeff_q(I,J,k) = GME_coeff
        str_xy_GME(I,J) = GME_coeff * sh_xy_bt(I,J)

      enddo ; enddo

      ! Applying GME diagonal term.  This is linear and the arguments can be rescaled.
      call smooth_GME(CS,G,GME_flux_h=str_xx_GME)
      call smooth_GME(CS,G,GME_flux_q=str_xy_GME)

      do J=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        str_xx(i,j) = (str_xx(i,j) + str_xx_GME(i,j)) * (h(i,j,k) * CS%reduction_xx(i,j))
      enddo ; enddo

      do J=js-1,Jeq ; do I=is-1,Ieq
        ! GME is applied below
        if (CS%no_slip) then
          str_xy(I,J) = (str_xy(I,J) + str_xy_GME(I,J)) * (hq(I,J) * CS%reduction_xy(I,J))
        else
          str_xy(I,J) = (str_xy(I,J) + str_xy_GME(I,J)) * (hq(I,J) * G%mask2dBu(I,J) * CS%reduction_xy(I,J))
        endif
      enddo ; enddo

      if (associated(MEKE%GME_snk)) then
        do j=js,je ; do i=is,ie
          FrictWork_GME(i,j,k) = GME_coeff_h(i,j,k) * h(i,j,k) * GV%H_to_kg_m2 * grad_vel_mag_bt_h(i,j)
        enddo ; enddo
      endif

    else ! use_GME
      do J=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        str_xx(i,j) = str_xx(i,j) * (h(i,j,k) * CS%reduction_xx(i,j))
      enddo ; enddo

      do J=js-1,Jeq ; do I=is-1,Ieq
        if (CS%no_slip) then
          str_xy(I,J) = str_xy(I,J) * (hq(I,J) * CS%reduction_xy(I,J))
        else
          str_xy(I,J) = str_xy(I,J) * (hq(I,J) * G%mask2dBu(I,J) * CS%reduction_xy(I,J))
        endif
      enddo ; enddo

    endif ! use_GME

    ! Evaluate 1/h x.Div(h Grad u) or the biharmonic equivalent.
    do j=js,je ; do I=Isq,Ieq
      diffu(I,j,k) = ((G%IdyCu(I,j)*(CS%dy2h(i,j) *str_xx(i,j) - &
                                     CS%dy2h(i+1,j)*str_xx(i+1,j)) + &
                       G%IdxCu(I,j)*(CS%dx2q(I,J-1)*str_xy(I,J-1) - &
                                     CS%dx2q(I,J) *str_xy(I,J))) * &
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
      diffv(i,J,k) = ((G%IdyCv(i,J)*(CS%dy2q(I-1,J)*str_xy(I-1,J) - &
                                    CS%dy2q(I,J) *str_xy(I,J)) - &
                       G%IdxCv(i,J)*(CS%dx2h(i,j) *str_xx(i,j) - &
                                    CS%dx2h(i,j+1)*str_xx(i,j+1))) * &
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

    if (find_FrictWork) then ; do j=js,je ; do i=is,ie
      ! Diagnose   str_xx*d_x u - str_yy*d_y v + str_xy*(d_y u + d_x v)
      ! This is the old formulation that includes energy diffusion
      FrictWork(i,j,k) = GV%H_to_RZ * ( &
              (str_xx(i,j)*(u(I,j,k)-u(I-1,j,k))*G%IdxT(i,j)     &
              -str_xx(i,j)*(v(i,J,k)-v(i,J-1,k))*G%IdyT(i,j))    &
       +0.25*((str_xy(I,J)*(                                     &
                   (u(I,j+1,k)-u(I,j,k))*G%IdyBu(I,J)            &
                  +(v(i+1,J,k)-v(i,J,k))*G%IdxBu(I,J) )          &
              +str_xy(I-1,J-1)*(                                 &
                   (u(I-1,j,k)-u(I-1,j-1,k))*G%IdyBu(I-1,J-1)    &
                  +(v(i,J-1,k)-v(i-1,J-1,k))*G%IdxBu(I-1,J-1) )) &
             +(str_xy(I-1,J)*(                                   &
                   (u(I-1,j+1,k)-u(I-1,j,k))*G%IdyBu(I-1,J)      &
                  +(v(i,J,k)-v(i-1,J,k))*G%IdxBu(I-1,J) )        &
              +str_xy(I,J-1)*(                                   &
                   (u(I,j,k)-u(I,j-1,k))*G%IdyBu(I,J-1)          &
                  +(v(i+1,J-1,k)-v(i,J-1,k))*G%IdxBu(I,J-1) )) ) )
    enddo ; enddo ; endif

    ! Make a similar calculation as for FrictWork above but accumulating into
    ! the vertically integrated MEKE source term, and adjusting for any
    ! energy loss seen as a reduction in the (biharmonic) frictional source term.
    if (find_FrictWork .and. associated(MEKE)) then ; if (associated(MEKE%mom_src)) then
      if (k==1) then
        do j=js,je ; do i=is,ie
          MEKE%mom_src(i,j) = 0.
        enddo ; enddo
        if (associated(MEKE%GME_snk)) then
          do j=js,je ; do i=is,ie
            MEKE%GME_snk(i,j) = 0.
          enddo ; enddo
        endif
      endif
      if (MEKE%backscatter_Ro_c /= 0.) then
        do j=js,je ; do i=is,ie
          FatH = 0.25*( (abs(G%CoriolisBu(I-1,J-1)) + abs(G%CoriolisBu(I,J))) + &
                        (abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J-1))) )
          Shear_mag = sqrt(sh_xx(i,j)*sh_xx(i,j) + &
            0.25*((sh_xy(I-1,J-1)*sh_xy(I-1,J-1) + sh_xy(I,J)*sh_xy(I,J)) + &
                  (sh_xy(I-1,J)*sh_xy(I-1,J) + sh_xy(I,J-1)*sh_xy(I,J-1))))
          if (CS%answers_2018) then
            FatH = (US%s_to_T*FatH)**MEKE%backscatter_Ro_pow ! f^n
            ! Note the hard-coded dimensional constant in the following line that can not
            ! be rescaled for dimensional consistency.
            Shear_mag = ( ( (US%s_to_T*Shear_mag)**MEKE%backscatter_Ro_pow ) + 1.e-30 ) &
                        * MEKE%backscatter_Ro_c ! c * D^n
            ! The Rossby number function is g(Ro) = 1/(1+c.Ro^n)
            ! RoScl = 1 - g(Ro)
            RoScl = Shear_mag / ( FatH + Shear_mag ) ! = 1 - f^n/(f^n+c*D^n)
          else
            if (FatH <= backscat_subround*Shear_mag) then
              RoScl = 1.0
            else
              Sh_F_pow = MEKE%backscatter_Ro_c * (Shear_mag / FatH)**MEKE%backscatter_Ro_pow
              RoScl = Sh_F_pow / (1.0 + Sh_F_pow) ! = 1 - f^n/(f^n+c*D^n)
            endif
          endif


          MEKE%mom_src(i,j) = MEKE%mom_src(i,j) + GV%H_to_RZ * ( &
                ((str_xx(i,j)-RoScl*bhstr_xx(i,j))*(u(I,j,k)-u(I-1,j,k))*G%IdxT(i,j)  &
                -(str_xx(i,j)-RoScl*bhstr_xx(i,j))*(v(i,J,k)-v(i,J-1,k))*G%IdyT(i,j)) &
         +0.25*(((str_xy(I,J)-RoScl*bhstr_xy(I,J))*(                                  &
                     (u(I,j+1,k)-u(I,j,k))*G%IdyBu(I,J)                               &
                    +(v(i+1,J,k)-v(i,J,k))*G%IdxBu(I,J) )                             &
                +(str_xy(I-1,J-1)-RoScl*bhstr_xy(I-1,J-1))*(                          &
                     (u(I-1,j,k)-u(I-1,j-1,k))*G%IdyBu(I-1,J-1)                       &
                    +(v(i,J-1,k)-v(i-1,J-1,k))*G%IdxBu(I-1,J-1) ))                    &
                +((str_xy(I-1,J)-RoScl*bhstr_xy(I-1,J))*(                             &
                     (u(I-1,j+1,k)-u(I-1,j,k))*G%IdyBu(I-1,J)                         &
                    +(v(i,J,k)-v(i-1,J,k))*G%IdxBu(I-1,J) )                           &
                +(str_xy(I,J-1)-RoScl*bhstr_xy(I,J-1))*(                              &
                     (u(I,j,k)-u(I,j-1,k))*G%IdyBu(I,J-1)                             &
                    +(v(i+1,J-1,k)-v(i,J-1,k))*G%IdxBu(I,J-1) )) ) )
        enddo ; enddo
      endif ! MEKE%backscatter_Ro_c

      do j=js,je ; do i=is,ie
        MEKE%mom_src(i,j) = MEKE%mom_src(i,j) + FrictWork(i,j,k)
      enddo ; enddo

      if (CS%use_GME .and. associated(MEKE)) then ; if (associated(MEKE%GME_snk)) then
        do j=js,je ; do i=is,ie
          MEKE%GME_snk(i,j) = MEKE%GME_snk(i,j) + FrictWork_GME(i,j,k)
        enddo ; enddo
      endif ; endif

    endif ; endif ! find_FrictWork and associated(mom_src)

  enddo ! end of k loop

  ! Offer fields for diagnostic averaging.
  if (CS%id_normstress > 0) call post_data(CS%id_normstress, NoSt, CS%diag)
  if (CS%id_shearstress > 0) call post_data(CS%id_shearstress, ShSt, CS%diag)
  if (CS%id_diffu>0)     call post_data(CS%id_diffu, diffu, CS%diag)
  if (CS%id_diffv>0)     call post_data(CS%id_diffv, diffv, CS%diag)
  if (CS%id_FrictWork>0) call post_data(CS%id_FrictWork, FrictWork, CS%diag)
  if (CS%id_FrictWork_GME>0) call post_data(CS%id_FrictWork_GME, FrictWork_GME, CS%diag)
  if (CS%id_Ah_h>0)      call post_data(CS%id_Ah_h, Ah_h, CS%diag)
  if (CS%id_grid_Re_Ah>0) call post_data(CS%id_grid_Re_Ah, grid_Re_Ah, CS%diag)
  if (CS%id_div_xx_h>0)  call post_data(CS%id_div_xx_h, div_xx_h, CS%diag)
  if (CS%id_vort_xy_q>0) call post_data(CS%id_vort_xy_q, vort_xy_q, CS%diag)
  if (CS%id_sh_xx_h>0)  call post_data(CS%id_sh_xx_h, sh_xx_h, CS%diag)
  if (CS%id_sh_xy_q>0) call post_data(CS%id_sh_xy_q, sh_xy_q, CS%diag)
  if (CS%id_Ah_q>0)      call post_data(CS%id_Ah_q, Ah_q, CS%diag)
  if (CS%id_Kh_h>0)      call post_data(CS%id_Kh_h, Kh_h, CS%diag)
  if (CS%id_grid_Re_Kh>0) call post_data(CS%id_grid_Re_Kh, grid_Re_Kh, CS%diag)
  if (CS%id_Kh_q>0)      call post_data(CS%id_Kh_q, Kh_q, CS%diag)
  if (CS%id_GME_coeff_h > 0)  call post_data(CS%id_GME_coeff_h, GME_coeff_h, CS%diag)
  if (CS%id_GME_coeff_q > 0)  call post_data(CS%id_GME_coeff_q, GME_coeff_q, CS%diag)

  if (CS%debug) then
    if (CS%Laplacian) then
      call hchksum(Kh_h, "Kh_h", G%HI, haloshift=0, scale=US%L_to_m**2*US%s_to_T)
      call Bchksum(Kh_q, "Kh_q", G%HI, haloshift=0, scale=US%L_to_m**2*US%s_to_T)
    endif
    if (CS%biharmonic) call hchksum(Ah_h, "Ah_h", G%HI, haloshift=0, scale=US%L_to_m**4*US%s_to_T)
    if (CS%biharmonic) call Bchksum(Ah_q, "Ah_q", G%HI, haloshift=0, scale=US%L_to_m**4*US%s_to_T)
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

  ! Diagnostics for terms multiplied by fractional thicknesses

  ! 3D diagnostics hf_diffu(diffv) are commented because there is no clarity on proper remapping grid option.
  ! The code is retained for degugging purposes in the future.
  !if (present(ADp) .and. (CS%id_hf_diffu > 0)) then
  !  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
  !    CS%hf_diffu(I,j,k) = diffu(I,j,k) * ADp%diag_hfrac_u(I,j,k)
  !  enddo ; enddo ; enddo
  !  call post_data(CS%id_hf_diffu, CS%hf_diffu, CS%diag)
  !endif
  !if (present(ADp) .and. (CS%id_hf_diffv > 0)) then
  !  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
  !    CS%hf_diffv(i,J,k) = diffv(i,J,k) * ADp%diag_hfrac_v(i,J,k)
  !  enddo ; enddo ; enddo
  !  call post_data(CS%id_hf_diffv, CS%hf_diffv, CS%diag)
  !endif
  if (present(ADp) .and. (CS%id_hf_diffu_2d > 0)) then
    allocate(hf_diffu_2d(G%IsdB:G%IedB,G%jsd:G%jed))
    hf_diffu_2d(:,:) = 0.0
    do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      hf_diffu_2d(I,j) = hf_diffu_2d(I,j) + diffu(I,j,k) * ADp%diag_hfrac_u(I,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_hf_diffu_2d, hf_diffu_2d, CS%diag)
    deallocate(hf_diffu_2d)
  endif
  if (present(ADp) .and. (CS%id_hf_diffv_2d > 0)) then
    allocate(hf_diffv_2d(G%isd:G%ied,G%JsdB:G%JedB))
    hf_diffv_2d(:,:) = 0.0
    do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      hf_diffv_2d(i,J) = hf_diffv_2d(i,J) + diffv(i,J,k) * ADp%diag_hfrac_v(i,J,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_hf_diffv_2d, hf_diffv_2d, CS%diag)
    deallocate(hf_diffv_2d)
  endif

end subroutine horizontal_viscosity

!> Allocates space for and calculates static variables used by horizontal_viscosity().
!! hor_visc_init calculates and stores the values of a number of metric functions that
!! are used in horizontal_viscosity().
subroutine hor_visc_init(Time, G, US, param_file, diag, CS, MEKE, ADp)
  type(time_type),         intent(in)    :: Time !< Current model time.
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< Structure to regulate diagnostic output.
  type(hor_visc_CS), pointer             :: CS   !< Pointer to the control structure for this module
  type(MEKE_type), pointer               :: MEKE !< MEKE data
  type(accel_diag_ptrs), optional, pointer :: ADp !< Acceleration diagnostic pointers
  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G)) :: u0u, u0v
  real, dimension(SZI_(G),SZJB_(G)) :: v0u, v0v
                ! u0v is the Laplacian sensitivities to the v velocities
                ! at u points [L-2 ~> m-2], with u0u, v0u, and v0v defined similarly.
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
  real :: Kh_vel_scale     ! this speed [L T-1 ~> m s-1] times grid spacing gives Lap visc
  real :: Ah_vel_scale     ! this speed [L T-1 ~> m s-1] times grid spacing cubed gives bih visc
  real :: Ah_time_scale    ! damping time-scale for biharmonic visc [T ~> s]
  real :: Smag_Lap_const   ! nondimensional Laplacian Smagorinsky constant
  real :: Smag_bi_const    ! nondimensional biharmonic Smagorinsky constant
  real :: Leith_Lap_const  ! nondimensional Laplacian Leith constant
  real :: Leith_bi_const   ! nondimensional biharmonic Leith constant
  real :: dt               ! The dynamics time step [T ~> s]
  real :: Idt              ! The inverse of dt [T-1 ~> s-1]
  real :: denom            ! work variable; the denominator of a fraction
  real :: maxvel           ! largest permitted velocity components [m s-1]
  real :: bound_Cor_vel    ! grid-scale velocity variations at which value
                           ! the quadratically varying biharmonic viscosity
                           ! balances Coriolis acceleration [L T-1 ~> m s-1]
  real :: Kh_sin_lat       ! Amplitude of latitudinally dependent viscosity [L2 T-1 ~> m2 s-1]
  real :: Kh_pwr_of_sine   ! Power used to raise sin(lat) when using Kh_sin_lat
  logical :: bound_Cor_def ! parameter setting of BOUND_CORIOLIS
  logical :: get_all       ! If true, read and log all parameters, regardless of
                           ! whether they are used, to enable spell-checking of
                           ! valid parameters.
  logical :: split         ! If true, use the split time stepping scheme.
                           ! If false and USE_GME = True, issue a FATAL error.
  logical :: default_2018_answers
  character(len=64) :: inputdir, filename
  real    :: deg2rad       ! Converts degrees to radians
  real    :: slat_fn       ! sin(lat)**Kh_pwr_of_sine
  real    :: aniso_grid_dir(2) ! Vector (n1,n2) for anisotropic direction
  integer :: aniso_mode    ! Selects the mode for setting the anisotropic direction
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  integer :: i, j
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_hor_visc"  ! module name
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec ; nz = G%ke
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  if (associated(CS)) then
    call MOM_error(WARNING, "hor_visc_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)
  CS%diag => diag
  ! Read parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  !   It is not clear whether these initialization lines are needed for the
  ! cases where the corresponding parameters are not read.
  CS%bound_Kh = .false. ; CS%better_bound_Kh = .false. ; CS%Smagorinsky_Kh = .false. ; CS%Leith_Kh = .false.
  CS%bound_Ah = .false. ; CS%better_bound_Ah = .false. ; CS%Smagorinsky_Ah = .false. ; CS%Leith_Ah = .false.
  CS%use_QG_Leith_visc = .false.
  CS%bound_Coriolis = .false.
  CS%Modified_Leith = .false.
  CS%anisotropic = .false.
  CS%dynamic_aniso = .false.
  Kh = 0.0 ; Ah = 0.0
  !   If GET_ALL_PARAMS is true, all parameters are read in all cases to enable
  ! parameter spelling checks.
  call get_param(param_file, mdl, "GET_ALL_PARAMS", get_all, default=.false.)
  call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false.)
  call get_param(param_file, mdl, "HOR_VISC_2018_ANSWERS", CS%answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)
  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false.)
  call get_param(param_file, mdl, "LAPLACIAN", CS%Laplacian, &
                 "If true, use a Laplacian horizontal viscosity.", &
                 default=.false.)
  if (CS%Laplacian .or. get_all) then
    call get_param(param_file, mdl, "KH", Kh,                      &
                 "The background Laplacian horizontal viscosity.", &
                 units = "m2 s-1", default=0.0, scale=US%m_to_L**2*US%T_to_s)
    call get_param(param_file, mdl, "KH_BG_MIN", CS%Kh_bg_min, &
                 "The minimum value allowed for Laplacian horizontal viscosity, KH.", &
                 units = "m2 s-1",  default=0.0, scale=US%m_to_L**2*US%T_to_s)
    call get_param(param_file, mdl, "KH_VEL_SCALE", Kh_vel_scale, &
                 "The velocity scale which is multiplied by the grid "//&
                 "spacing to calculate the Laplacian viscosity. "//&
                 "The final viscosity is the largest of this scaled "//&
                 "viscosity, the Smagorinsky and Leith viscosities, and KH.", &
                 units="m s-1", default=0.0, scale=US%m_s_to_L_T)
    call get_param(param_file, mdl, "KH_SIN_LAT", Kh_sin_lat, &
                 "The amplitude of a latitudinally-dependent background "//&
                 "viscosity of the form KH_SIN_LAT*(SIN(LAT)**KH_PWR_OF_SINE).", &
                 units = "m2 s-1",  default=0.0, scale=US%m_to_L**2*US%T_to_s)
    if (Kh_sin_lat>0. .or. get_all) &
      call get_param(param_file, mdl, "KH_PWR_OF_SINE", Kh_pwr_of_sine, &
                 "The power used to raise SIN(LAT) when using a latitudinally "//&
                 "dependent background viscosity.", &
                 units = "nondim",  default=4.0)
    call get_param(param_file, mdl, "SMAGORINSKY_KH", CS%Smagorinsky_Kh, &
                 "If true, use a Smagorinsky nonlinear eddy viscosity.", &
                 default=.false.)
    if (CS%Smagorinsky_Kh .or. get_all) &
      call get_param(param_file, mdl, "SMAG_LAP_CONST", Smag_Lap_const, &
                 "The nondimensional Laplacian Smagorinsky constant, "//&
                 "often 0.15.", units="nondim", default=0.0, &
                  fail_if_missing = CS%Smagorinsky_Kh)
    call get_param(param_file, mdl, "LEITH_KH", CS%Leith_Kh, &
                 "If true, use a Leith nonlinear eddy viscosity.", &
                 default=.false.)
    call get_param(param_file, mdl, "MODIFIED_LEITH", CS%Modified_Leith, &
                 "If true, add a term to Leith viscosity which is "//&
                 "proportional to the gradient of divergence.", &
                 default=.false.)
    call get_param(param_file, mdl, "RES_SCALE_MEKE_VISC", CS%res_scale_MEKE, &
                 "If true, the viscosity contribution from MEKE is scaled by "//&
                 "the resolution function.", default=.false.)
    if (CS%Leith_Kh .or. get_all) then
      call get_param(param_file, mdl, "LEITH_LAP_CONST", Leith_Lap_const, &
                 "The nondimensional Laplacian Leith constant, "//&
                 "often set to 1.0", units="nondim", default=0.0, &
                  fail_if_missing = CS%Leith_Kh)
      call get_param(param_file, mdl, "USE_QG_LEITH_VISC", CS%use_QG_Leith_visc, &
                 "If true, use QG Leith nonlinear eddy viscosity.", &
                 default=.false.)
      if (CS%use_QG_Leith_visc .and. .not. CS%Leith_Kh) call MOM_error(FATAL, &
                 "MOM_lateral_mixing_coeffs.F90, VarMix_init:"//&
                 "LEITH_KH must be True when USE_QG_LEITH_VISC=True.")
    endif
    if (CS%Leith_Kh .or. CS%Leith_Ah .or. get_all) then
      call get_param(param_file, mdl, "USE_BETA_IN_LEITH", CS%use_beta_in_Leith, &
                 "If true, include the beta term in the Leith nonlinear eddy viscosity.", &
                 default=CS%Leith_Kh)
      call get_param(param_file, mdl, "MODIFIED_LEITH", CS%modified_Leith, &
                 "If true, add a term to Leith viscosity which is \n"//&
                 "proportional to the gradient of divergence.", &
                 default=.false.)
    endif
    call get_param(param_file, mdl, "BOUND_KH", CS%bound_Kh, &
                 "If true, the Laplacian coefficient is locally limited "//&
                 "to be stable.", default=.true.)
    call get_param(param_file, mdl, "BETTER_BOUND_KH", CS%better_bound_Kh, &
                 "If true, the Laplacian coefficient is locally limited "//&
                 "to be stable with a better bounding than just BOUND_KH.", &
                 default=CS%bound_Kh)
    call get_param(param_file, mdl, "ANISOTROPIC_VISCOSITY", CS%anisotropic, &
                 "If true, allow anistropic viscosity in the Laplacian "//&
                 "horizontal viscosity.", default=.false.)
    call get_param(param_file, mdl, "ADD_LES_VISCOSITY", CS%add_LES_viscosity, &
                 "If true, adds the viscosity from Smagorinsky and Leith to the "//&
                 "background viscosity instead of taking the maximum.", default=.false.)
  endif
  if (CS%anisotropic .or. get_all) then
    call get_param(param_file, mdl, "KH_ANISO", CS%Kh_aniso, &
                 "The background Laplacian anisotropic horizontal viscosity.", &
                 units = "m2 s-1", default=0.0, scale=US%m_to_L**2*US%T_to_s)
    call get_param(param_file, mdl, "ANISOTROPIC_MODE", aniso_mode, &
                 "Selects the mode for setting the direction of anistropy.\n"//&
                 "\t 0 - Points along the grid i-direction.\n"//&
                 "\t 1 - Points towards East.\n"//&
                 "\t 2 - Points along the flow direction, U/|U|.", &
                 default=0)
    select case (aniso_mode)
      case (0)
        call get_param(param_file, mdl, "ANISO_GRID_DIR", aniso_grid_dir, &
                 "The vector pointing in the direction of anistropy for "//&
                 "horizont viscosity. n1,n2 are the i,j components relative "//&
                 "to the grid.", units = "nondim", fail_if_missing=.true.)
      case (1)
        call get_param(param_file, mdl, "ANISO_GRID_DIR", aniso_grid_dir, &
                 "The vector pointing in the direction of anistropy for "//&
                 "horizont viscosity. n1,n2 are the i,j components relative "//&
                 "to the spherical coordinates.", units = "nondim", fail_if_missing=.true.)
    end select
  endif
  call get_param(param_file, mdl, "BIHARMONIC", CS%biharmonic, &
                 "If true, use a biharmonic horizontal viscosity. "//&
                 "BIHARMONIC may be used with LAPLACIAN.", &
                 default=.true.)
  if (CS%biharmonic .or. get_all) then
    call get_param(param_file, mdl, "AH", Ah, &
                 "The background biharmonic horizontal viscosity.", &
                 units = "m4 s-1", default=0.0, scale=US%m_to_L**4*US%T_to_s)
    call get_param(param_file, mdl, "AH_VEL_SCALE", Ah_vel_scale, &
                 "The velocity scale which is multiplied by the cube of "//&
                 "the grid spacing to calculate the biharmonic viscosity. "//&
                 "The final viscosity is the largest of this scaled "//&
                 "viscosity, the Smagorinsky and Leith viscosities, and AH.", &
                 units="m s-1", default=0.0, scale=US%m_s_to_L_T)
    call get_param(param_file, mdl, "AH_TIME_SCALE", Ah_time_scale, &
                 "A time scale whose inverse is multiplied by the fourth "//&
                 "power of the grid spacing to calculate biharmonic viscosity. "//&
                 "The final viscosity is the largest of all viscosity "//&
                 "formulations in use. 0.0 means that it's not used.", &
                 units="s", default=0.0, scale=US%s_to_T)
    call get_param(param_file, mdl, "SMAGORINSKY_AH", CS%Smagorinsky_Ah, &
                 "If true, use a biharmonic Smagorinsky nonlinear eddy "//&
                 "viscosity.", default=.false.)
    call get_param(param_file, mdl, "LEITH_AH", CS%Leith_Ah, &
                 "If true, use a biharmonic Leith nonlinear eddy "//&
                 "viscosity.", default=.false.)
    call get_param(param_file, mdl, "BOUND_AH", CS%bound_Ah, &
                 "If true, the biharmonic coefficient is locally limited "//&
                 "to be stable.", default=.true.)
    call get_param(param_file, mdl, "BETTER_BOUND_AH", CS%better_bound_Ah, &
                 "If true, the biharmonic coefficient is locally limited "//&
                 "to be stable with a better bounding than just BOUND_AH.", &
                 default=CS%bound_Ah)
    call get_param(param_file, mdl, "RE_AH", CS%Re_Ah, &
                 "If nonzero, the biharmonic coefficient is scaled "//&
                 "so that the biharmonic Reynolds number is equal to this.", &
                 units="nondim", default=0.0)

    if (CS%Smagorinsky_Ah .or. get_all) then
      call get_param(param_file, mdl, "SMAG_BI_CONST",Smag_bi_const, &
                 "The nondimensional biharmonic Smagorinsky constant, "//&
                 "typically 0.015 - 0.06.", units="nondim", default=0.0, &
                 fail_if_missing = CS%Smagorinsky_Ah)
      call get_param(param_file, mdl, "BOUND_CORIOLIS", bound_Cor_def, default=.false.)
      call get_param(param_file, mdl, "BOUND_CORIOLIS_BIHARM", CS%bound_Coriolis, &
                 "If true use a viscosity that increases with the square "//&
                 "of the velocity shears, so that the resulting viscous "//&
                 "drag is of comparable magnitude to the Coriolis terms "//&
                 "when the velocity differences between adjacent grid "//&
                 "points is 0.5*BOUND_CORIOLIS_VEL.  The default is the "//&
                 "value of BOUND_CORIOLIS (or false).", default=bound_Cor_def)
      if (CS%bound_Coriolis .or. get_all) then
        call get_param(param_file, mdl, "MAXVEL", maxvel, default=3.0e8)
        bound_Cor_vel = maxvel
        call get_param(param_file, mdl, "BOUND_CORIOLIS_VEL", bound_Cor_vel, &
                 "The velocity scale at which BOUND_CORIOLIS_BIHARM causes "//&
                 "the biharmonic drag to have comparable magnitude to the "//&
                 "Coriolis acceleration.  The default is set by MAXVEL.", &
                 units="m s-1", default=maxvel, scale=US%m_s_to_L_T)
      endif
    endif
    if (CS%Leith_Ah .or. get_all) &
      call get_param(param_file, mdl, "LEITH_BI_CONST", Leith_bi_const, &
                 "The nondimensional biharmonic Leith constant, "//&
                 "typical values are thus far undetermined.", units="nondim", default=0.0, &
                 fail_if_missing = CS%Leith_Ah)
  endif
  call get_param(param_file, mdl, "USE_LAND_MASK_FOR_HVISC", CS%use_land_mask, &
                 "If true, use Use the land mask for the computation of thicknesses "//&
                 "at velocity locations. This eliminates the dependence on arbitrary "//&
                 "values over land or outside of the domain.", default=.true.)
  if (CS%better_bound_Ah .or. CS%better_bound_Kh .or. get_all) &
    call get_param(param_file, mdl, "HORVISC_BOUND_COEF", CS%bound_coef, &
                 "The nondimensional coefficient of the ratio of the "//&
                 "viscosity bounds to the theoretical maximum for "//&
                 "stability without considering other terms.", units="nondim", &
                 default=0.8)
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
                 "terms and this background value.", default=.false.)

  call get_param(param_file, mdl, "USE_GME", CS%use_GME, &
                 "If true, use the GM+E backscatter scheme in association \n"//&
                 "with the Gent and McWilliams parameterization.", default=.false.)
  if (CS%use_GME) then
    call get_param(param_file, mdl, "SPLIT", split, &
                 "Use the split time stepping if true.", default=.true., &
                  do_not_log=.true.)
    if (.not. split) call MOM_error(FATAL,"ERROR: Currently, USE_GME = True "// &
                                           "cannot be used with SPLIT=False.")
    call get_param(param_file, mdl, "GME_H0", CS%GME_h0, &
                 "The strength of GME tapers quadratically to zero when the bathymetric "//&
                 "depth is shallower than GME_H0.", units="m", scale=US%m_to_Z, &
                 default=1000.0)
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
  if (CS%use_Kh_bg_2d) then
    ALLOC_(CS%Kh_bg_2d(isd:ied,jsd:jed))     ; CS%Kh_bg_2d(:,:) = 0.0
    call get_param(param_file, mdl, "KH_BG_2D_FILENAME", filename, &
                 'The filename containing a 2d map of "Kh".', &
                 default='KH_background_2d.nc')
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    call MOM_read_data(trim(inputdir)//trim(filename), 'Kh', CS%Kh_bg_2d, &
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
    if (CS%Leith_Ah) then
        ALLOC_(CS%biharm6_const_xx(isd:ied,jsd:jed)) ; CS%biharm6_const_xx(:,:) = 0.0
        ALLOC_(CS%biharm6_const_xy(IsdB:IedB,JsdB:JedB)) ; CS%biharm6_const_xy(:,:) = 0.0
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
  do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
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
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
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
        CS%Kh_bg_xy(I,J) = MAX(CS%Kh_bg_xy(I,J), &
            0.25*((CS%Kh_bg_2d(i,j) + CS%Kh_bg_2d(i+1,j+1)) + &
                  (CS%Kh_bg_2d(i+1,j) + CS%Kh_bg_2d(i,j+1))) )
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
    do j=js-1,Jeq+1 ; do I=Isq-1,Ieq+1
      CS%Idx2dyCu(I,j) = (G%IdxCu(I,j)*G%IdxCu(I,j)) * G%IdyCu(I,j)
      CS%Idxdy2u(I,j) = G%IdxCu(I,j) * (G%IdyCu(I,j)*G%IdyCu(I,j))
    enddo ; enddo
    do J=Jsq-1,Jeq+1 ; do i=is-1,Ieq+1
      CS%Idx2dyCv(i,J) = (G%IdxCv(i,J)*G%IdxCv(i,J)) * G%IdyCv(i,J)
      CS%Idxdy2v(i,J) = G%IdxCv(i,J) * (G%IdyCv(i,J)*G%IdyCv(i,J))
    enddo ; enddo
    CS%Ah_bg_xy(:,:) = 0.0
   ! The 0.3 below was 0.4 in MOM1.10.  The change in hq requires
   ! this to be less than 1/3, rather than 1/2 as before.
    if (CS%better_bound_Ah .or. CS%bound_Ah) Ah_Limit = 0.3 / (dt*64.0)
    if (CS%Smagorinsky_Ah .and. CS%bound_Coriolis) &
      BoundCorConst = 1.0 / (5.0*(bound_Cor_vel*bound_Cor_vel))

    min_grid_sp_h4 = huge(1.)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
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
      if (CS%Leith_Ah) then
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
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
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
      call hchksum(CS%Kh_Max_xx, "Kh_Max_xx", G%HI, haloshift=0, scale=US%L_to_m**2*US%s_to_T)
      call Bchksum(CS%Kh_Max_xy, "Kh_Max_xy", G%HI, haloshift=0, scale=US%L_to_m**2*US%s_to_T)
    endif
  endif
  ! The biharmonic bounds should avoid overshoots when CS%bound_coef < 0.5, but
  ! empirically work for CS%bound_coef <~ 1.0
  if (CS%biharmonic .and. CS%better_bound_Ah) then
    do j=js-1,Jeq+1 ; do I=Isq-1,Ieq+1
      u0u(I,j) = (CS%Idxdy2u(I,j)*(CS%dy2h(i+1,j)*CS%DY_dxT(i+1,j)*(G%IdyCu(I+1,j) + G%IdyCu(I,j))   + &
                                   CS%dy2h(i,j) * CS%DY_dxT(i,j) * (G%IdyCu(I,j) + G%IdyCu(I-1,j)) ) + &
                 CS%Idx2dyCu(I,j)*(CS%dx2q(I,J) * CS%DX_dyBu(I,J) * (G%IdxCu(I,j+1) + G%IdxCu(I,j)) + &
                                   CS%dx2q(I,J-1)*CS%DX_dyBu(I,J-1)*(G%IdxCu(I,j) + G%IdxCu(I,j-1)) ) )
      u0v(I,j) = (CS%Idxdy2u(I,j)*(CS%dy2h(i+1,j)*CS%DX_dyT(i+1,j)*(G%IdxCv(i+1,J) + G%IdxCv(i+1,J-1)) + &
                                   CS%dy2h(i,j) * CS%DX_dyT(i,j) * (G%IdxCv(i,J) + G%IdxCv(i,J-1)) )   + &
                 CS%Idx2dyCu(I,j)*(CS%dx2q(I,J) * CS%DY_dxBu(I,J) * (G%IdyCv(i+1,J) + G%IdyCv(i,J))   + &
                                   CS%dx2q(I,J-1)*CS%DY_dxBu(I,J-1)*(G%IdyCv(i+1,J-1) + G%IdyCv(i,J-1)) ) )
    enddo ; enddo
    do J=Jsq-1,Jeq+1 ; do i=is-1,Ieq+1
      v0u(i,J) = (CS%Idxdy2v(i,J)*(CS%dy2q(I,J) * CS%DX_dyBu(I,J) * (G%IdxCu(I,j+1) + G%IdxCu(I,j))       + &
                                   CS%dy2q(I-1,J)*CS%DX_dyBu(I-1,J)*(G%IdxCu(I-1,j+1) + G%IdxCu(I-1,j)) ) + &
                 CS%Idx2dyCv(i,J)*(CS%dx2h(i,j+1)*CS%DY_dxT(i,j+1)*(G%IdyCu(I,j+1) + G%IdyCu(I-1,j+1))   + &
                                   CS%dx2h(i,j) * CS%DY_dxT(i,j) * (G%IdyCu(I,j) + G%IdyCu(I-1,j)) ) )
      v0v(i,J) = (CS%Idxdy2v(i,J)*(CS%dy2q(I,J) * CS%DY_dxBu(I,J) * (G%IdyCv(i+1,J) + G%IdyCv(i,J))   + &
                                   CS%dy2q(I-1,J)*CS%DY_dxBu(I-1,J)*(G%IdyCv(i,J) + G%IdyCv(i-1,J)) ) + &
                 CS%Idx2dyCv(i,J)*(CS%dx2h(i,j+1)*CS%DX_dyT(i,j+1)*(G%IdxCv(i,J+1) + G%IdxCv(i,J))   + &
                                   CS%dx2h(i,j) * CS%DX_dyT(i,j) * (G%IdxCv(i,J) + G%IdxCv(i,J-1)) ) )
    enddo ; enddo
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
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
      call hchksum(CS%Ah_Max_xx, "Ah_Max_xx", G%HI, haloshift=0, scale=US%L_to_m**4*US%s_to_T)
      call Bchksum(CS%Ah_Max_xy, "Ah_Max_xy", G%HI, haloshift=0, scale=US%L_to_m**4*US%s_to_T)
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
  !    'Fractional Thickness-weighted Zonal Acceleration from Horizontal Viscosity', 'm s-2', &
  !    v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if ((CS%id_hf_diffu > 0) .and. (present(ADp))) then
  !  call safe_alloc_ptr(CS%hf_diffu,G%IsdB,G%IedB,G%jsd,G%jed,G%ke)
  !  call safe_alloc_ptr(ADp%diag_hfrac_u,G%IsdB,G%IedB,G%jsd,G%jed,G%ke)
  !endif

  !CS%id_hf_diffv = register_diag_field('ocean_model', 'hf_diffv', diag%axesCvL, Time, &
  !    'Fractional Thickness-weighted Meridional Acceleration from Horizontal Viscosity', 'm s-2', &
  !    v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if ((CS%id_hf_diffv > 0) .and. (present(ADp))) then
  !  call safe_alloc_ptr(CS%hf_diffv,G%isd,G%ied,G%JsdB,G%JedB,G%ke)
  !  call safe_alloc_ptr(ADp%diag_hfrac_v,G%isd,G%ied,G%JsdB,G%JedB,G%ke)
  !endif

  CS%id_hf_diffu_2d = register_diag_field('ocean_model', 'hf_diffu_2d', diag%axesCu1, Time, &
      'Depth-sum Fractional Thickness-weighted Zonal Acceleration from Horizontal Viscosity', 'm s-2', &
      conversion=US%L_T2_to_m_s2)
  if ((CS%id_hf_diffu_2d > 0) .and. (present(ADp))) then
    call safe_alloc_ptr(ADp%diag_hfrac_u,G%IsdB,G%IedB,G%jsd,G%jed,G%ke)
  endif

  CS%id_hf_diffv_2d = register_diag_field('ocean_model', 'hf_diffv_2d', diag%axesCv1, Time, &
      'Depth-sum Fractional Thickness-weighted Meridional Acceleration from Horizontal Viscosity', 'm s-2', &
      conversion=US%L_T2_to_m_s2)
  if ((CS%id_hf_diffv_2d > 0) .and. (present(ADp))) then
    call safe_alloc_ptr(ADp%diag_hfrac_v,G%isd,G%ied,G%JsdB,G%JedB,G%ke)
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
      CS%id_GME_coeff_h = register_diag_field('ocean_model', 'GME_coeff_h', diag%axesTL, Time, &
        'GME coefficient at h Points', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
      CS%id_GME_coeff_q = register_diag_field('ocean_model', 'GME_coeff_q', diag%axesBL, Time, &
        'GME coefficient at q Points', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
      CS%id_FrictWork_GME = register_diag_field('ocean_model','FrictWork_GME',diag%axesTL,Time,&
      'Integral work done by lateral friction terms in GME (excluding diffusion of energy)', &
      'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  endif
  CS%id_FrictWork = register_diag_field('ocean_model','FrictWork',diag%axesTL,Time,&
      'Integral work done by lateral friction terms', &
      'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  CS%id_FrictWorkIntz = register_diag_field('ocean_model','FrictWorkIntz',diag%axesT1,Time,      &
      'Depth integrated work done by lateral friction', &
      'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2, &
      cmor_field_name='dispkexyfo',                                                              &
      cmor_long_name='Depth integrated ocean kinetic energy dissipation due to lateral friction',&
      cmor_standard_name='ocean_kinetic_energy_dissipation_per_unit_area_due_to_xy_friction')
  if (CS%Laplacian .or. get_all) then
  endif
end subroutine hor_visc_init
!> Calculates factors in the anisotropic orientation tensor to be align with the grid.
!! With n1=1 and n2=0, this recovers the approach of Large et al, 2001.
subroutine align_aniso_tensor_to_grid(CS, n1, n2)
  type(hor_visc_CS), pointer :: CS !< Control structure for horizontal viscosity
  real,              intent(in) :: n1 !< i-component of direction vector [nondim]
  real,              intent(in) :: n2 !< j-component of direction vector [nondim]
  ! Local variables
  real :: recip_n2_norm
  ! For normalizing n=(n1,n2) in case arguments are not a unit vector
  recip_n2_norm = n1**2 + n2**2
  if (recip_n2_norm > 0.) recip_n2_norm = 1./recip_n2_norm
  CS%n1n2_h(:,:) = 2. * ( n1 * n2 ) * recip_n2_norm
  CS%n1n2_q(:,:) = 2. * ( n1 * n2 ) * recip_n2_norm
  CS%n1n1_m_n2n2_h(:,:) = ( n1 * n1 - n2 * n2 ) * recip_n2_norm
  CS%n1n1_m_n2n2_q(:,:) = ( n1 * n1 - n2 * n2 ) * recip_n2_norm
end subroutine align_aniso_tensor_to_grid
!> Apply a 1-1-4-1-1 Laplacian filter one time on GME diffusive flux to reduce any
!! horizontal two-grid-point noise
subroutine smooth_GME(CS,G,GME_flux_h,GME_flux_q)
  ! Arguments
  type(hor_visc_CS),                            pointer       :: CS        !< Control structure
  type(ocean_grid_type),                        intent(in)    :: G         !< Ocean grid
  real, dimension(SZI_(G),SZJ_(G)),   optional, intent(inout) :: GME_flux_h!< GME diffusive flux
                                                              !! at h points
  real, dimension(SZIB_(G),SZJB_(G)), optional, intent(inout) :: GME_flux_q!< GME diffusive flux
                                                              !! at q points
  ! local variables
  real, dimension(SZI_(G),SZJ_(G)) :: GME_flux_h_original
  real, dimension(SZIB_(G),SZJB_(G)) :: GME_flux_q_original
  real :: wc, ww, we, wn, ws ! averaging weights for smoothing
  integer :: i, j, k, s
  do s=1,1
    ! Update halos
    if (present(GME_flux_h)) then
      call pass_var(GME_flux_h, G%Domain)
      GME_flux_h_original = GME_flux_h
      ! apply smoothing on GME
      do j = G%jsc, G%jec
        do i = G%isc, G%iec
          ! skip land points
          if (G%mask2dT(i,j)==0.) cycle
          ! compute weights
          ww = 0.125 * G%mask2dT(i-1,j)
          we = 0.125 * G%mask2dT(i+1,j)
          ws = 0.125 * G%mask2dT(i,j-1)
          wn = 0.125 * G%mask2dT(i,j+1)
          wc = 1.0 - (ww+we+wn+ws)
          GME_flux_h(i,j) =  wc * GME_flux_h_original(i,j)   &
                           + ww * GME_flux_h_original(i-1,j) &
                           + we * GME_flux_h_original(i+1,j) &
                           + ws * GME_flux_h_original(i,j-1) &
                           + wn * GME_flux_h_original(i,j+1)
      enddo; enddo
    endif
    ! Update halos
    if (present(GME_flux_q)) then
      call pass_var(GME_flux_q, G%Domain, position=CORNER, complete=.true.)
      GME_flux_q_original = GME_flux_q
      ! apply smoothing on GME
      do J = G%JscB, G%JecB
        do I = G%IscB, G%IecB
          ! skip land points
          if (G%mask2dBu(I,J)==0.) cycle
          ! compute weights
          ww = 0.125 * G%mask2dBu(I-1,J)
          we = 0.125 * G%mask2dBu(I+1,J)
          ws = 0.125 * G%mask2dBu(I,J-1)
          wn = 0.125 * G%mask2dBu(I,J+1)
          wc = 1.0 - (ww+we+wn+ws)
          GME_flux_q(I,J) =  wc * GME_flux_q_original(I,J)   &
                           + ww * GME_flux_q_original(I-1,J) &
                           + we * GME_flux_q_original(I+1,J) &
                           + ws * GME_flux_q_original(I,J-1) &
                           + wn * GME_flux_q_original(I,J+1)
      enddo; enddo
    endif
  enddo ! s-loop
end subroutine smooth_GME
!> Deallocates any variables allocated in hor_visc_init.
subroutine hor_visc_end(CS)
  type(hor_visc_CS), pointer :: CS !< The control structure returned by a
                                   !! previous call to hor_visc_init.
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
      DEALLOC_(CS%Biharm6_const_xx) ; DEALLOC_(CS%Biharm6_const_xy)
    endif
    if (CS%Leith_Ah) then
      DEALLOC_(CS%Biharm_const_xx) ; DEALLOC_(CS%Biharm_const_xy)
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
  deallocate(CS)
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
!! + \partial_y \left( \frac{1}{2} \sigma_T \right)
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
