!> Baropotric solver
module MOM_barotropic

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging, only : hchksum, uvchksum
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, query_averaging_enabled, register_diag_field
use MOM_diag_mediator, only : safe_alloc_ptr, diag_ctrl, enable_averaging
use MOM_domains, only : min_across_PEs, clone_MOM_domain, pass_vector
use MOM_domains, only : To_All, Scalar_Pair, AGRID, CORNER, MOM_domain_type
use MOM_domains, only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains, only : start_group_pass, complete_group_pass, pass_var
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : mech_forcing
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_io, only : vardesc, var_desc, MOM_read_data, slasher
use MOM_open_boundary, only : ocean_OBC_type, OBC_SIMPLE, OBC_NONE, open_boundary_query
use MOM_open_boundary, only : OBC_DIRECTION_E, OBC_DIRECTION_W
use MOM_open_boundary, only : OBC_DIRECTION_N, OBC_DIRECTION_S, OBC_segment_type
use MOM_restart, only : register_restart_field, register_restart_pair
use MOM_restart, only : query_initialized, MOM_restart_CS
use MOM_tidal_forcing, only : tidal_forcing_sensitivity, tidal_forcing_CS
use MOM_time_manager, only : time_type, real_to_time, operator(+), operator(-)
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : BT_cont_type, alloc_bt_cont_type
use MOM_verticalGrid, only : verticalGrid_type
use MOM_variables, only : accel_diag_ptrs

implicit none ; private

#include <MOM_memory.h>
#ifdef STATIC_MEMORY_
#  ifndef BTHALO_
#    define BTHALO_ 0
#  endif
#  define WHALOI_ MAX(BTHALO_-NIHALO_,0)
#  define WHALOJ_ MAX(BTHALO_-NJHALO_,0)
#  define NIMEMW_   1-WHALOI_:NIMEM_+WHALOI_
#  define NJMEMW_   1-WHALOJ_:NJMEM_+WHALOJ_
#  define NIMEMBW_  -WHALOI_:NIMEM_+WHALOI_
#  define NJMEMBW_  -WHALOJ_:NJMEM_+WHALOJ_
#  define SZIW_(G)  NIMEMW_
#  define SZJW_(G)  NJMEMW_
#  define SZIBW_(G) NIMEMBW_
#  define SZJBW_(G) NJMEMBW_
#else
#  define NIMEMW_   :
#  define NJMEMW_   :
#  define NIMEMBW_  :
#  define NJMEMBW_  :
#  define SZIW_(G)  G%isdw:G%iedw
#  define SZJW_(G)  G%jsdw:G%jedw
#  define SZIBW_(G) G%isdw-1:G%iedw
#  define SZJBW_(G) G%jsdw-1:G%jedw
#endif

public btcalc, bt_mass_source, btstep, barotropic_init, barotropic_end
public register_barotropic_restarts, set_dtbt, barotropic_get_tav

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The barotropic stepping open boundary condition type
type, private :: BT_OBC_type
  real, dimension(:,:), pointer :: Cg_u => NULL()  !< The external wave speed at u-points [L T-1 ~> m s-1].
  real, dimension(:,:), pointer :: Cg_v => NULL()  !< The external wave speed at u-points [L T-1 ~> m s-1].
  real, dimension(:,:), pointer :: H_u => NULL()   !< The total thickness at the u-points [H ~> m or kg m-2].
  real, dimension(:,:), pointer :: H_v => NULL()   !< The total thickness at the v-points [H ~> m or kg m-2].
  real, dimension(:,:), pointer :: uhbt => NULL()  !< The zonal barotropic thickness fluxes specified
                                     !! for open boundary conditions (if any) [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(:,:), pointer :: vhbt => NULL()  !< The meridional barotropic thickness fluxes specified
                                     !! for open boundary conditions (if any) [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(:,:), pointer :: ubt_outer => NULL() !< The zonal velocities just outside the domain,
                                     !! as set by the open boundary conditions [L T-1 ~> m s-1].
  real, dimension(:,:), pointer :: vbt_outer => NULL() !< The meridional velocities just outside the domain,
                                     !! as set by the open boundary conditions [L T-1 ~> m s-1].
  real, dimension(:,:), pointer :: eta_outer_u => NULL() !< The surface height outside of the domain
                                     !! at a u-point with an open boundary condition [H ~> m or kg m-2].
  real, dimension(:,:), pointer :: eta_outer_v => NULL() !< The surface height outside of the domain
                                     !! at a v-point with an open boundary condition [H ~> m or kg m-2].
  logical :: apply_u_OBCs !< True if this PE has an open boundary at a u-point.
  logical :: apply_v_OBCs !< True if this PE has an open boundary at a v-point.
  !>@{ Index ranges for the open boundary conditions
  integer :: is_u_obc, ie_u_obc, js_u_obc, je_u_obc
  integer :: is_v_obc, ie_v_obc, js_v_obc, je_v_obc
  !>@}
  logical :: is_alloced = .false. !< True if BT_OBC is in use and has been allocated

  type(group_pass_type) :: pass_uv   !< Structure for group halo pass
  type(group_pass_type) :: pass_uhvh !< Structure for group halo pass
  type(group_pass_type) :: pass_h    !< Structure for group halo pass
  type(group_pass_type) :: pass_cg   !< Structure for group halo pass
  type(group_pass_type) :: pass_eta_outer  !< Structure for group halo pass
end type BT_OBC_type

!> The barotropic stepping control stucture
type, public :: barotropic_CS ; private
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: frhatu
          !< The fraction of the total column thickness interpolated to u grid points in each layer [nondim].
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: frhatv
          !< The fraction of the total column thickness interpolated to v grid points in each layer [nondim].
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: IDatu
          !< Inverse of the basin depth at u grid points [Z-1 ~> m-1].
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: lin_drag_u
          !< A spatially varying linear drag coefficient acting on the zonal barotropic flow
          !! [H T-1 ~> m s-1 or kg m-2 s-1].
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: ubt_IC
          !< The barotropic solvers estimate of the zonal velocity that will be the initial
          !! condition for the next call to btstep [L T-1 ~> m s-1].
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: ubtav
          !< The barotropic zonal velocity averaged over the baroclinic time step [L T-1 ~> m s-1].
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: IDatv
          !< Inverse of the basin depth at v grid points [Z-1 ~> m-1].
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: lin_drag_v
          !< A spatially varying linear drag coefficient acting on the zonal barotropic flow
          !! [H T-1 ~> m s-1 or kg m-2 s-1].
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: vbt_IC
          !< The barotropic solvers estimate of the zonal velocity that will be the initial
          !! condition for the next call to btstep [L T-1 ~> m s-1].
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: vbtav
          !< The barotropic meridional velocity averaged over the  baroclinic time step [L T-1 ~> m s-1].
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: eta_cor
          !< The difference between the free surface height from the barotropic calculation and the sum
          !! of the layer thicknesses. This difference is imposed as a forcing term in the barotropic
          !! calculation over a baroclinic timestep [H ~> m or kg m-2].
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: eta_cor_bound
          !< A limit on the rate at which eta_cor can be applied while avoiding instability
          !! [H T-1 ~> m s-1 or kg m-2 s-1]. This is only used if CS%bound_BT_corr is true.
  real ALLOCABLE_, dimension(NIMEMW_,NJMEMW_) :: &
    ua_polarity, &  !< Test vector components for checking grid polarity.
    va_polarity, &  !< Test vector components for checking grid polarity.
    bathyT          !< A copy of bathyT (ocean bottom depth) with wide halos [Z ~> m]
  real ALLOCABLE_, dimension(NIMEMW_,NJMEMW_) :: IareaT
                    !<   This is a copy of G%IareaT with wide halos, but will
                    !! still utilize the macro IareaT when referenced, [L-2 ~> m-2].
  real ALLOCABLE_, dimension(NIMEMBW_,NJMEMW_) :: &
    D_u_Cor, &      !<   A simply averaged depth at u points [Z ~> m].
    dy_Cu, &        !<   A copy of G%dy_Cu with wide halos [L ~> m].
    IdxCu           !<   A copy of G%IdxCu with wide halos [L-1 ~> m-1].
  real ALLOCABLE_, dimension(NIMEMW_,NJMEMBW_) :: &
    D_v_Cor, &      !<   A simply averaged depth at v points [Z ~> m].
    dx_Cv, &        !<   A copy of G%dx_Cv with wide halos [L ~> m].
    IdyCv           !<   A copy of G%IdyCv with wide halos [L-1 ~> m-1].
  real ALLOCABLE_, dimension(NIMEMBW_,NJMEMBW_) :: &
    q_D             !< f / D at PV points [Z-1 T-1 ~> m-1 s-1].

  real, dimension(:,:,:), pointer :: frhatu1 => NULL() !< Predictor step values of frhatu stored for diagnostics.
  real, dimension(:,:,:), pointer :: frhatv1 => NULL() !< Predictor step values of frhatv stored for diagnostics.

  type(BT_OBC_type) :: BT_OBC !< A structure with all of this modules fields
                              !! for applying open boundary conditions.

  real    :: dtbt            !< The barotropic time step [T ~> s].
  real    :: dtbt_fraction   !<   The fraction of the maximum time-step that
                             !! should used.  The default is 0.98.
  real    :: dtbt_max        !<   The maximum stable barotropic time step [T ~> s].
  real    :: dt_bt_filter    !<   The time-scale over which the barotropic mode solutions are
                             !! filtered [T ~> s] if positive, or as a fraction of DT if
                             !! negative [nondim].  This can never be taken to be longer than 2*dt.
                             !! Set this to 0 to apply no filtering.
  integer :: nstep_last = 0  !< The number of barotropic timesteps per baroclinic
                             !! time step the last time btstep was called.
  real    :: bebt            !< A nondimensional number, from 0 to 1, that
                             !! determines the gravity wave time stepping scheme.
                             !! 0.0 gives a forward-backward scheme, while 1.0
                             !! give backward Euler. In practice, bebt should be
                             !! of order 0.2 or greater.
  logical :: split           !< If true, use the split time stepping scheme.
  logical :: bound_BT_corr   !< If true, the magnitude of the fake mass source
                             !! in the barotropic equation that drives the two
                             !! estimates of the free surface height toward each
                             !! other is bounded to avoid driving corrective
                             !! velocities that exceed MAXCFL_BT_CONT.
  logical :: gradual_BT_ICs  !< If true, adjust the initial conditions for the
                             !! barotropic solver to the values from the layered
                             !! solution over a whole timestep instead of
                             !! instantly.  This is a decent approximation to the
                             !! inclusion of sum(u dh_dt) while also correcting
                             !! for truncation errors.
  logical :: Sadourny        !< If true, the Coriolis terms are discretized
                             !! with Sadourny's energy conserving scheme,
                             !! otherwise the Arakawa & Hsu scheme is used.  If
                             !! the deformation radius is not resolved Sadourny's
                             !! scheme should probably be used.
  logical :: integral_bt_cont !< If true, use the time-integrated velocity over the barotropic steps
                             !! to determine the integrated transports used to update the continuity
                             !! equation.  Otherwise the transports are the sum of the transports
                             !! based on ]a series of instantaneous velocities and the BT_CONT_TYPE
                             !! for transports.  This is only valid if a BT_CONT_TYPE is used.
  logical :: Nonlinear_continuity !< If true, the barotropic continuity equation
                             !! uses the full ocean thickness for transport.
  integer :: Nonlin_cont_update_period !< The number of barotropic time steps
                             !! between updates to the face area, or 0 only to
                             !! update at the start of a call to btstep.  The
                             !! default is 1.
  logical :: BT_project_velocity !< If true, step the barotropic velocity first
                             !! and project out the velocity tendency by 1+BEBT
                             !! when calculating the transport.  The default
                             !! (false) is to use a predictor continuity step to
                             !! find the pressure field, and then do a corrector
                             !! continuity step using a weighted average of the
                             !! old and new velocities, with weights of (1-BEBT) and BEBT.
  logical :: nonlin_stress   !< If true, use the full depth of the ocean at the start of the
                             !! barotropic step when calculating the surface stress contribution to
                             !! the barotropic acclerations.  Otherwise use the depth based on bathyT.
  real    :: BT_Coriolis_scale !< A factor by which the barotropic Coriolis acceleration anomaly
                             !! terms are scaled.
  logical :: answers_2018    !< If true, use expressions for the barotropic solver that recover
                             !! the answers from the end of 2018.  Otherwise, use more efficient
                             !! or general expressions.

  logical :: dynamic_psurf   !< If true, add a dynamic pressure due to a viscous
                             !! ice shelf, for instance.
  real    :: Dmin_dyn_psurf  !< The minimum depth to use in limiting the size
                             !! of the dynamic surface pressure for stability [Z ~> m].
  real    :: ice_strength_length  !< The length scale at which the damping rate
                             !! due to the ice strength should be the same as if
                             !! a Laplacian were applied [L ~> m].
  real    :: const_dyn_psurf !< The constant that scales the dynamic surface
                             !! pressure [nondim].  Stable values are < ~1.0.
                             !! The default is 0.9.
  logical :: tides           !< If true, apply tidal momentum forcing.
  real    :: G_extra         !< A nondimensional factor by which gtot is enhanced.
  integer :: hvel_scheme     !< An integer indicating how the thicknesses at
                             !! velocity points are calculated. Valid values are
                             !! given by the parameters defined below:
                             !!   HARMONIC, ARITHMETIC, HYBRID, and FROM_BT_CONT
  logical :: strong_drag     !< If true, use a stronger estimate of the retarding
                             !! effects of strong bottom drag.
  logical :: linear_wave_drag  !< If true, apply a linear drag to the barotropic
                             !! velocities, using rates set by lin_drag_u & _v
                             !! divided by the depth of the ocean.
  logical :: linearized_BT_PV  !< If true, the PV and interface thicknesses used
                             !! in the barotropic Coriolis calculation is time
                             !! invariant and linearized.
  logical :: use_wide_halos  !< If true, use wide halos and march in during the
                             !! barotropic time stepping for efficiency.
  logical :: clip_velocity   !< If true, limit any velocity components that are
                             !! are large enough for a CFL number to exceed
                             !! CFL_trunc.  This should only be used as a
                             !! desperate debugging measure.
  logical :: debug           !< If true, write verbose checksums for debugging purposes.
  logical :: debug_bt        !< If true, write verbose checksums for debugging purposes.
  real    :: vel_underflow   !< Velocity components smaller than vel_underflow
                             !! are set to 0 [L T-1 ~> m s-1].
  real    :: maxvel          !< Velocity components greater than maxvel are
                             !! truncated to maxvel [L T-1 ~> m s-1].
  real    :: CFL_trunc       !< If clip_velocity is true, velocity components will
                             !! be truncated when they are large enough that the
                             !! corresponding CFL number exceeds this value, nondim.
  real    :: maxCFL_BT_cont  !< The maximum permitted CFL number associated with the
                             !! barotropic accelerations from the summed velocities
                             !! times the time-derivatives of thicknesses.  The
                             !! default is 0.1, and there will probably be real
                             !! problems if this were set close to 1.
  logical :: BT_cont_bounds  !< If true, use the BT_cont_type variables to set limits
                             !! on the magnitude of the corrective mass fluxes.
  logical :: visc_rem_u_uh0  !< If true, use the viscous remnants when estimating
                             !! the barotropic velocities that were used to
                             !! calculate uh0 and vh0.  False is probably the
                             !! better choice.
  logical :: adjust_BT_cont  !< If true, adjust the curve fit to the BT_cont type
                             !! that is used by the barotropic solver to match the
                             !! transport about which the flow is being linearized.
  logical :: use_old_coriolis_bracket_bug !< If True, use an order of operations
                             !! that is not bitwise rotationally symmetric in the
                             !! meridional Coriolis term of the barotropic solver.
  type(time_type), pointer :: Time  => NULL() !< A pointer to the ocean models clock.
  type(diag_ctrl), pointer :: diag => NULL()  !< A structure that is used to regulate
                             !! the timing of diagnostic output.
  type(MOM_domain_type), pointer :: BT_Domain => NULL()  !< Barotropic MOM domain
  type(hor_index_type), pointer :: debug_BT_HI => NULL() !< debugging copy of horizontal index_type
  type(tidal_forcing_CS), pointer :: tides_CSp => NULL() !< Control structure for tides
  logical :: module_is_initialized = .false.  !< If true, module has been initialized

  integer :: isdw !< The lower i-memory limit for the wide halo arrays.
  integer :: iedw !< The upper i-memory limit for the wide halo arrays.
  integer :: jsdw !< The lower j-memory limit for the wide halo arrays.
  integer :: jedw !< The upper j-memory limit for the wide halo arrays.

  type(group_pass_type) :: pass_q_DCor !< Handle for a group halo pass
  type(group_pass_type) :: pass_gtot !< Handle for a group halo pass
  type(group_pass_type) :: pass_tmp_uv !< Handle for a group halo pass
  type(group_pass_type) :: pass_eta_bt_rem !< Handle for a group halo pass
  type(group_pass_type) :: pass_force_hbt0_Cor_ref !< Handle for a group halo pass
  type(group_pass_type) :: pass_Dat_uv !< Handle for a group halo pass
  type(group_pass_type) :: pass_eta_ubt !< Handle for a group halo pass
  type(group_pass_type) :: pass_etaav !< Handle for a group halo pass
  type(group_pass_type) :: pass_ubt_Cor !< Handle for a group halo pass
  type(group_pass_type) :: pass_ubta_uhbta !< Handle for a group halo pass
  type(group_pass_type) :: pass_e_anom !< Handle for a group halo pass

  !>@{ Diagnostic IDs
  integer :: id_PFu_bt = -1, id_PFv_bt = -1, id_Coru_bt = -1, id_Corv_bt = -1
  integer :: id_ubtforce = -1, id_vbtforce = -1, id_uaccel = -1, id_vaccel = -1
  integer :: id_visc_rem_u = -1, id_visc_rem_v = -1, id_eta_cor = -1
  integer :: id_ubt = -1, id_vbt = -1, id_eta_bt = -1, id_ubtav = -1, id_vbtav = -1
  integer :: id_ubt_st = -1, id_vbt_st = -1, id_eta_st = -1
  integer :: id_ubtdt = -1, id_vbtdt = -1
  integer :: id_ubt_hifreq = -1, id_vbt_hifreq = -1, id_eta_hifreq = -1
  integer :: id_uhbt_hifreq = -1, id_vhbt_hifreq = -1, id_eta_pred_hifreq = -1
  integer :: id_gtotn = -1, id_gtots = -1, id_gtote = -1, id_gtotw = -1
  integer :: id_uhbt = -1, id_frhatu = -1, id_vhbt = -1, id_frhatv = -1
  integer :: id_frhatu1 = -1, id_frhatv1 = -1

  integer :: id_BTC_FA_u_EE = -1, id_BTC_FA_u_E0 = -1, id_BTC_FA_u_W0 = -1, id_BTC_FA_u_WW = -1
  integer :: id_BTC_ubt_EE = -1, id_BTC_ubt_WW = -1
  integer :: id_BTC_FA_v_NN = -1, id_BTC_FA_v_N0 = -1, id_BTC_FA_v_S0 = -1, id_BTC_FA_v_SS = -1
  integer :: id_BTC_vbt_NN = -1, id_BTC_vbt_SS = -1
  integer :: id_BTC_FA_u_rat0 = -1, id_BTC_FA_v_rat0 = -1, id_BTC_FA_h_rat0 = -1
  integer :: id_uhbt0 = -1, id_vhbt0 = -1
  !>@}

end type barotropic_CS

!> A desciption of the functional dependence of transport at a u-point
type, private :: local_BT_cont_u_type
  real :: FA_u_EE !< The effective open face area for zonal barotropic transport
                  !! drawing from locations far to the east [H L ~> m2 or kg m-1].
  real :: FA_u_E0 !< The effective open face area for zonal barotropic transport
                  !! drawing from nearby to the east [H L ~> m2 or kg m-1].
  real :: FA_u_W0 !< The effective open face area for zonal barotropic transport
                  !! drawing from nearby to the west [H L ~> m2 or kg m-1].
  real :: FA_u_WW !< The effective open face area for zonal barotropic transport
                  !! drawing from locations far to the west [H L ~> m2 or kg m-1].
  real :: uBT_WW  !< uBT_WW is the barotropic velocity [L T-1 ~> m s-1], or with INTEGRAL_BT_CONTINUITY
                  !! the time-integrated barotropic velocity [L ~> m], beyond which the marginal
                  !! open face area is FA_u_WW.  uBT_WW must be non-negative.
  real :: uBT_EE  !< uBT_EE is a barotropic velocity [L T-1 ~> m s-1], or with INTEGRAL_BT_CONTINUITY
                  !! the time-integrated barotropic velocity [L ~> m], beyond which the marginal
                  !! open face area is FA_u_EE. uBT_EE must be non-positive.
  real :: uh_crvW !< The curvature of face area with velocity for flow from the west [H T2 L-1 ~> s2 or kg s2 m-3]
                  !! or [H L-1 ~> 1 or kg m-3] with INTEGRAL_BT_CONTINUITY.
  real :: uh_crvE !< The curvature of face area with velocity for flow from the east [H T2 L-1 ~> s2 or kg s2 m-3]
                  !! or [H L-1 ~> 1 or kg m-3] with INTEGRAL_BT_CONTINUITY.
  real :: uh_WW   !< The zonal transport when ubt=ubt_WW [H L2 T-1 ~> m3 s-1 or kg s-1], or the equivalent
                  !! time-integrated transport with INTEGRAL_BT_CONTINUITY [H L2 ~> m3 or kg].
  real :: uh_EE   !< The zonal transport when ubt=ubt_EE [H L2 T-1 ~> m3 s-1 or kg s-1], or the equivalent
                  !! time-integrated transport with INTEGRAL_BT_CONTINUITY [H L2 ~> m3 or kg].
end type local_BT_cont_u_type

!> A desciption of the functional dependence of transport at a v-point
type, private :: local_BT_cont_v_type
  real :: FA_v_NN !< The effective open face area for meridional barotropic transport
                  !! drawing from locations far to the north [H L ~> m2 or kg m-1].
  real :: FA_v_N0 !< The effective open face area for meridional barotropic transport
                  !! drawing from nearby to the north [H L ~> m2 or kg m-1].
  real :: FA_v_S0 !< The effective open face area for meridional barotropic transport
                  !! drawing from nearby to the south [H L ~> m2 or kg m-1].
  real :: FA_v_SS !< The effective open face area for meridional barotropic transport
                  !! drawing from locations far to the south [H L ~> m2 or kg m-1].
  real :: vBT_SS  !< vBT_SS is the barotropic velocity [L T-1 ~> m s-1], or with INTEGRAL_BT_CONTINUITY
                  !! the time-integrated barotropic velocity [L ~> m], beyond which the marginal
                  !! open face area is FA_v_SS. vBT_SS must be non-negative.
  real :: vBT_NN  !< vBT_NN is the barotropic velocity [L T-1 ~> m s-1], or with INTEGRAL_BT_CONTINUITY
                  !! the time-integrated barotropic velocity [L ~> m], beyond which the marginal
                  !! open face area is FA_v_NN.  vBT_NN must be non-positive.
  real :: vh_crvS !< The curvature of face area with velocity for flow from the south [H T2 L-1 ~> s2 or kg s2 m-3]
                  !! or [H L-1 ~> 1 or kg m-3] with INTEGRAL_BT_CONTINUITY.
  real :: vh_crvN !< The curvature of face area with velocity for flow from the north [H T2 L-1 ~> s2 or kg s2 m-3]
                  !! or [H L-1 ~> 1 or kg m-3] with INTEGRAL_BT_CONTINUITY.
  real :: vh_SS   !< The meridional transport when vbt=vbt_SS [H L2 T-1 ~> m3 s-1 or kg s-1], or the equivalent
                  !! time-integrated transport with INTEGRAL_BT_CONTINUITY [H L2 ~> m3 or kg].
  real :: vh_NN   !< The meridional transport when vbt=vbt_NN [H L2 T-1 ~> m3 s-1 or kg s-1], or the equivalent
                  !! time-integrated transport with INTEGRAL_BT_CONTINUITY [H L2 ~> m3 or kg].
end type local_BT_cont_v_type

!> A container for passing around active tracer point memory limits
type, private :: memory_size_type
  !>@{ Currently active memory limits
  integer :: isdw, iedw, jsdw, jedw ! The memory limits of the wide halo arrays.
  !>@}
end type memory_size_type

!>@{ CPU time clock IDs
integer :: id_clock_sync=-1, id_clock_calc=-1
integer :: id_clock_calc_pre=-1, id_clock_calc_post=-1
integer :: id_clock_pass_step=-1, id_clock_pass_pre=-1, id_clock_pass_post=-1
!>@}

!>@{ Enumeration values for various schemes
integer, parameter :: HARMONIC        = 1
integer, parameter :: ARITHMETIC      = 2
integer, parameter :: HYBRID          = 3
integer, parameter :: FROM_BT_CONT    = 4
integer, parameter :: HYBRID_BT_CONT  = 5
character*(20), parameter :: HYBRID_STRING = "HYBRID"
character*(20), parameter :: HARMONIC_STRING = "HARMONIC"
character*(20), parameter :: ARITHMETIC_STRING = "ARITHMETIC"
character*(20), parameter :: BT_CONT_STRING = "FROM_BT_CONT"
!>@}

contains

!> This subroutine time steps the barotropic equations explicitly.
!! For gravity waves, anything between a forwards-backwards scheme
!! and a simulated backwards Euler scheme is used, with bebt between
!! 0.0 and 1.0 determining the scheme.  In practice, bebt must be of
!! order 0.2 or greater.  A forwards-backwards treatment of the
!! Coriolis terms is always used.
subroutine btstep(U_in, V_in, eta_in, dt, bc_accel_u, bc_accel_v, forces, pbce, &
                  eta_PF_in, U_Cor, V_Cor, accel_layer_u, accel_layer_v, &
                  eta_out, uhbtav, vhbtav, G, GV, US, CS, &
                  visc_rem_u, visc_rem_v, etaav, ADp, OBC, BT_cont, eta_PF_start, &
                  taux_bot, tauy_bot, uh0, vh0, u_uh0, v_vh0)
  type(ocean_grid_type),                   intent(inout) :: G       !< The ocean's grid structure.
  type(verticalGrid_type),                   intent(in)  :: GV      !< The ocean's vertical grid structure.
  type(unit_scale_type),                     intent(in)  :: US      !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)  :: U_in    !< The initial (3-D) zonal
                                                                    !! velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)  :: V_in    !< The initial (3-D) meridional
                                                                    !! velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G)),          intent(in)  :: eta_in  !< The initial barotropic free surface height
                                                         !! anomaly or column mass anomaly [H ~> m or kg m-2].
  real,                                      intent(in)  :: dt      !< The time increment to integrate over [T ~> s].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)  :: bc_accel_u !< The zonal baroclinic accelerations,
                                                                       !! [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)  :: bc_accel_v !< The meridional baroclinic accelerations,
                                                                       !! [L T-2 ~> m s-2].
  type(mech_forcing),                        intent(in)  :: forces     !< A structure with the driving mechanical forces
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: pbce       !< The baroclinic pressure anomaly in each layer
                                                         !! due to free surface height anomalies
                                                         !! [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
  real, dimension(SZI_(G),SZJ_(G)),          intent(in)  :: eta_PF_in  !< The 2-D eta field (either SSH anomaly or
                                                         !! column mass anomaly) that was used to calculate the input
                                                         !! pressure gradient accelerations (or its final value if
                                                         !! eta_PF_start is provided [H ~> m or kg m-2].
                                                         !! Note: eta_in, pbce, and eta_PF_in must have up-to-date
                                                         !! values in the first point of their halos.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)  :: U_Cor      !< The (3-D) zonal velocities used to
                                                         !! calculate the Coriolis terms in bc_accel_u [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)  :: V_Cor      !< The (3-D) meridional velocities used to
                                                         !! calculate the Coriolis terms in bc_accel_u [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out) :: accel_layer_u !< The zonal acceleration of each layer due
                                                         !! to the barotropic calculation [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: accel_layer_v !< The meridional acceleration of each layer
                                                         !! due to the barotropic calculation [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJ_(G)),          intent(out) :: eta_out       !< The final barotropic free surface
                                                         !! height anomaly or column mass anomaly [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G)),         intent(out) :: uhbtav        !< the barotropic zonal volume or mass
                                                         !! fluxes averaged through the barotropic steps
                                                         !! [H L2 T-1 ~> m3 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)),         intent(out) :: vhbtav        !< the barotropic meridional volume or mass
                                                         !! fluxes averaged through the barotropic steps
                                                         !! [H L2 T-1 ~> m3 or kg s-1].
  type(barotropic_CS),                       pointer     :: CS            !< The control structure returned by a
                                                         !! previous call to barotropic_init.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)  :: visc_rem_u    !< Both the fraction of the momentum
                                                         !! originally in a layer that remains after a time-step of
                                                         !! viscosity, and the fraction of a time-step's worth of a
                                                         !! barotropic acceleration that a layer experiences after
                                                         !! viscosity is applied, in the zonal direction. Nondimensional
                                                         !! between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)  :: visc_rem_v    !< Ditto for meridional direction [nondim].
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(out) :: etaav        !< The free surface height or column mass
                                                         !! averaged over the barotropic integration [H ~> m or kg m-2].
  type(accel_diag_ptrs),               optional, pointer :: ADp          !< Acceleration diagnostic pointers
  type(ocean_OBC_type),                optional, pointer :: OBC          !< The open boundary condition structure.
  type(BT_cont_type),                  optional, pointer :: BT_cont      !< A structure with elements that describe
                                                         !! the effective open face areas as a function of barotropic
                                                         !! flow.
  real, dimension(:,:),                optional, pointer :: eta_PF_start !< The eta field consistent with the pressure
                                                         !! gradient at the start of the barotropic stepping
                                                         !! [H ~> m or kg m-2].
  real, dimension(:,:),                optional, pointer :: taux_bot     !< The zonal bottom frictional stress from
                                                         !! ocean to the seafloor [R L Z T-2 ~> Pa].
  real, dimension(:,:),                optional, pointer :: tauy_bot     !< The meridional bottom frictional stress
                                                         !! from ocean to the seafloor [R L Z T-2 ~> Pa].
  real, dimension(:,:,:),              optional, pointer :: uh0     !< The zonal layer transports at reference
                                                                    !! velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(:,:,:),              optional, pointer :: u_uh0   !< The velocities used to calculate
                                                                    !! uh0 [L T-1 ~> m s-1]
  real, dimension(:,:,:),              optional, pointer :: vh0     !< The zonal layer transports at reference
                                                                    !! velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(:,:,:),              optional, pointer :: v_vh0   !< The velocities used to calculate
                                                                    !! vh0 [L T-1 ~> m s-1]

  ! Local variables
  real :: ubt_Cor(SZIB_(G),SZJ_(G)) ! The barotropic velocities that had been
  real :: vbt_Cor(SZI_(G),SZJB_(G)) ! used to calculate the input Coriolis
                                    ! terms [L T-1 ~> m s-1].
  real :: wt_u(SZIB_(G),SZJ_(G),SZK_(GV)) ! wt_u and wt_v are the
  real :: wt_v(SZI_(G),SZJB_(G),SZK_(GV)) ! normalized weights to
                ! be used in calculating barotropic velocities, possibly with
                ! sums less than one due to viscous losses.  Nondimensional.
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    av_rem_u, &   ! The weighted average of visc_rem_u, nondimensional.
    tmp_u, &      ! A temporary array at u points.
    ubt_st, &     ! The zonal barotropic velocity at the start of timestep [L T-1 ~> m s-1].
    ubt_dt        ! The zonal barotropic velocity tendency [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G)) :: &
    av_rem_v, &   ! The weighted average of visc_rem_v, nondimensional.
    tmp_v, &      ! A temporary array at v points.
    vbt_st, &     ! The meridional barotropic velocity at the start of timestep [L T-1 ~> m s-1].
    vbt_dt        ! The meridional barotropic velocity tendency [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    tmp_h, &      ! A temporary array at h points.
    e_anom        ! The anomaly in the sea surface height or column mass
                  ! averaged between the beginning and end of the time step,
                  ! relative to eta_PF, with SAL effects included [H ~> m or kg m-2].

  ! These are always allocated with symmetric memory and wide halos.
  real :: q(SZIBW_(CS),SZJBW_(CS))  ! A pseudo potential vorticity [T-1 Z-1 ~> s-1 m-1]
                  ! or [T-1 H-1 ~> s-1 m-1 or m2 s-1 kg-1]
  real, dimension(SZIBW_(CS),SZJW_(CS)) :: &
    ubt, &        ! The zonal barotropic velocity [L T-1 ~> m s-1].
    bt_rem_u, &   ! The fraction of the barotropic zonal velocity that remains
                  ! after a time step, the remainder being lost to bottom drag.
                  ! bt_rem_u is a nondimensional number between 0 and 1.
    BT_force_u, & ! The vertical average of all of the u-accelerations that are
                  ! not explicitly included in the barotropic equation [L T-2 ~> m s-2].
    u_accel_bt, & ! The difference between the zonal acceleration from the
                  ! barotropic calculation and BT_force_u [L T-2 ~> m s-2].
    uhbt, &       ! The zonal barotropic thickness fluxes [H L2 T-1 ~> m3 s-1 or kg s-1].
    uhbt0, &      ! The difference between the sum of the layer zonal thickness
                  ! fluxes and the barotropic thickness flux using the same
                  ! velocity [H L2 T-1 ~> m3 s-1 or kg s-1].
    ubt_old, &    ! The starting value of ubt in a barotropic step [L T-1 ~> m s-1].
    ubt_first, &  ! The starting value of ubt in a series of barotropic steps [L T-1 ~> m s-1].
    ubt_sum, &    ! The sum of ubt over the time steps [L T-1 ~> m s-1].
    ubt_int, &    ! The running time integral of ubt over the time steps [L ~> m].
    uhbt_sum, &   ! The sum of uhbt over the time steps [H L2 T-1 ~> m3 s-1 or kg s-1].
    uhbt_int, &   ! The running time integral of uhbt over the time steps [H L2  ~> m3].
    ubt_wtd, &    ! A weighted sum used to find the filtered final ubt [L T-1 ~> m s-1].
    ubt_trans, &  ! The latest value of ubt used for a transport [L T-1 ~> m s-1].
    azon, bzon, & ! _zon & _mer are the values of the Coriolis force which
    czon, dzon, & ! are applied to the neighboring values of vbtav & ubtav,
    amer, bmer, & ! respectively to get the barotropic inertial rotation
    cmer, dmer, & ! [T-1 ~> s-1].
    Cor_u, &      ! The zonal Coriolis acceleration [L T-2 ~> m s-2].
    Cor_ref_u, &  ! The zonal barotropic Coriolis acceleration due
                  ! to the reference velocities [L T-2 ~> m s-2].
    PFu, &        ! The zonal pressure force acceleration [L T-2 ~> m s-2].
    Rayleigh_u, & ! A Rayleigh drag timescale operating at u-points [T-1 ~> s-1].
    PFu_bt_sum, & ! The summed zonal barotropic pressure gradient force [L T-2 ~> m s-2].
    Coru_bt_sum, & ! The summed zonal barotropic Coriolis acceleration [L T-2 ~> m s-2].
    DCor_u, &     ! An averaged depth or total thickness at u points [Z ~> m] or [H ~> m or kg m-2].
    Datu          ! Basin depth at u-velocity grid points times the y-grid
                  ! spacing [H L ~> m2 or kg m-1].
  real, dimension(SZIW_(CS),SZJBW_(CS)) :: &
    vbt, &        ! The meridional barotropic velocity [L T-1 ~> m s-1].
    bt_rem_v, &   ! The fraction of the barotropic meridional velocity that
                  ! remains after a time step, the rest being lost to bottom
                  ! drag.  bt_rem_v is a nondimensional number between 0 and 1.
    BT_force_v, & ! The vertical average of all of the v-accelerations that are
                  ! not explicitly included in the barotropic equation [L T-2 ~> m s-2].
    v_accel_bt, & ! The difference between the meridional acceleration from the
                  ! barotropic calculation and BT_force_v [L T-2 ~> m s-2].
    vhbt, &       ! The meridional barotropic thickness fluxes [H L2 T-1 ~> m3 s-1 or kg s-1].
    vhbt0, &      ! The difference between the sum of the layer meridional
                  ! thickness fluxes and the barotropic thickness flux using
                  ! the same velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
    vbt_old, &    ! The starting value of vbt in a barotropic step [L T-1 ~> m s-1].
    vbt_first, &  ! The starting value of ubt in a series of barotropic steps [L T-1 ~> m s-1].
    vbt_sum, &    ! The sum of vbt over the time steps [L T-1 ~> m s-1].
    vbt_int, &    ! The running time integral of vbt over the time steps [L ~> m].
    vhbt_sum, &   ! The sum of vhbt over the time steps [H L2 T-1 ~> m3 s-1 or kg s-1].
    vhbt_int, &   ! The running time integral of vhbt over the time steps [H L2  ~> m3].
    vbt_wtd, &    ! A weighted sum used to find the filtered final vbt [L T-1 ~> m s-1].
    vbt_trans, &  ! The latest value of vbt used for a transport [L T-1 ~> m s-1].
    Cor_v, &      ! The meridional Coriolis acceleration [L T-2 ~> m s-2].
    Cor_ref_v, &  ! The meridional barotropic Coriolis acceleration due
                  ! to the reference velocities [L T-2 ~> m s-2].
    PFv, &        ! The meridional pressure force acceleration [L T-2 ~> m s-2].
    Rayleigh_v, & ! A Rayleigh drag timescale operating at v-points [T-1 ~> s-1].
    PFv_bt_sum, & ! The summed meridional barotropic pressure gradient force,
                  ! [L T-2 ~> m s-2].
    Corv_bt_sum, & ! The summed meridional barotropic Coriolis acceleration,
                  ! [L T-2 ~> m s-2].
    DCor_v, &     ! An averaged depth or total thickness at v points [Z ~> m] or [H ~> m or kg m-2].
    Datv          ! Basin depth at v-velocity grid points times the x-grid
                  ! spacing [H L ~> m2 or kg m-1].
  real, target, dimension(SZIW_(CS),SZJW_(CS)) :: &
    eta, &        ! The barotropic free surface height anomaly or column mass
                  ! anomaly [H ~> m or kg m-2]
    eta_pred      ! A predictor value of eta [H ~> m or kg m-2] like eta.
  real, dimension(:,:), pointer :: &
    eta_PF_BT     ! A pointer to the eta array (either eta or eta_pred) that
                  ! determines the barotropic pressure force [H ~> m or kg m-2]
  real, dimension(SZIW_(CS),SZJW_(CS)) :: &
    eta_sum, &    ! eta summed across the timesteps [H ~> m or kg m-2].
    eta_wtd, &    ! A weighted estimate used to calculate eta_out [H ~> m or kg m-2].
    eta_IC, &     ! A local copy of the initial 2-D eta field (eta_in) [H ~> m or kg m-2]
    eta_PF, &     ! A local copy of the 2-D eta field (either SSH anomaly or
                  ! column mass anomaly) that was used to calculate the input
                  ! pressure gradient accelerations [H ~> m or kg m-2].
    eta_PF_1, &   ! The initial value of eta_PF, when interp_eta_PF is
                  ! true [H ~> m or kg m-2].
    d_eta_PF, &   ! The change in eta_PF over the barotropic time stepping when
                  ! interp_eta_PF is true [H ~> m or kg m-2].
    gtot_E, &     ! gtot_X is the effective total reduced gravity used to relate
    gtot_W, &     ! free surface height deviations to pressure forces (including
    gtot_N, &     ! GFS and baroclinic  contributions) in the barotropic momentum
    gtot_S, &     ! equations half a grid-point in the X-direction (X is N, S, E, or W)
                  ! from the thickness point [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
                  ! (See Hallberg, J Comp Phys 1997 for a discussion.)
    eta_src, &    ! The source of eta per barotropic timestep [H ~> m or kg m-2].
    dyn_coef_eta, & ! The coefficient relating the changes in eta to the
                  ! dynamic surface pressure under rigid ice
                  ! [L2 T-2 H-1 ~> m s-2 or m4 s-2 kg-1].
    p_surf_dyn    ! A dynamic surface pressure under rigid ice [L2 T-2 ~> m2 s-2].
  type(local_BT_cont_u_type), dimension(SZIBW_(CS),SZJW_(CS)) :: &
    BTCL_u        ! A repackaged version of the u-point information in BT_cont.
  type(local_BT_cont_v_type), dimension(SZIW_(CS),SZJBW_(CS)) :: &
    BTCL_v        ! A repackaged version of the v-point information in BT_cont.
  ! End of wide-sized variables.

  real, dimension(SZIBW_(CS),SZJW_(CS)) :: &
    ubt_prev, ubt_sum_prev, ubt_wtd_prev, & ! Previous velocities stored for OBCs [L T-1 ~> m s-1]
    uhbt_prev, uhbt_sum_prev, & ! Previous transports stored for OBCs [L2 H T-1 ~> m3 s-1]
    ubt_int_prev, & ! Previous value of time-integrated velocity stored for OBCs [L ~> m]
    uhbt_int_prev   ! Previous value of time-integrated transport stored for OBCs [L2 H ~> m3]
  real, dimension(SZIW_(CS),SZJBW_(CS)) :: &
    vbt_prev, vbt_sum_prev, vbt_wtd_prev, & ! Previous velocities stored for OBCs [L T-1 ~> m s-1]
    vhbt_prev, vhbt_sum_prev, & ! Previous transports stored for OBCs [L2 H T-1 ~> m3 s-1]
    vbt_int_prev, & ! Previous value of time-integrated velocity stored for OBCs [L ~> m]
    vhbt_int_prev   ! Previous value of time-integrated transport stored for OBCs [L2 H ~> m3]
  real :: mass_to_Z   ! The depth unit conversion divided by the mean density (Rho0) [Z m-1 R-1 ~> m3 kg-1].
  real :: mass_accel_to_Z ! The inverse of the mean density (Rho0) [R-1 ~> m3 kg-1].
  real :: visc_rem    ! A work variable that may equal visc_rem_[uv].  Nondim.
  real :: vel_prev    ! The previous velocity [L T-1 ~> m s-1].
  real :: dtbt        ! The barotropic time step [T ~> s].
  real :: bebt        ! A copy of CS%bebt [nondim].
  real :: be_proj     ! The fractional amount by which velocities are projected
                      ! when project_velocity is true. For now be_proj is set
                      ! to equal bebt, as they have similar roles and meanings.
  real :: Idt         ! The inverse of dt [T-1 ~> s-1].
  real :: det_de      ! The partial derivative due to self-attraction and loading
                      ! of the reference geopotential with the sea surface height.
                      ! This is typically ~0.09 or less.
  real :: dgeo_de     ! The constant of proportionality between geopotential and
                      ! sea surface height.  It is a nondimensional number of
                      ! order 1.  For stability, this may be made larger
                      ! than the physical problem would suggest.
  real :: Instep      ! The inverse of the number of barotropic time steps to take.
  real :: wt_end      ! The weighting of the final value of eta_PF [nondim]
  integer :: nstep    ! The number of barotropic time steps to take.
  type(time_type) :: &
    time_bt_start, &  ! The starting time of the barotropic steps.
    time_step_end, &  ! The end time of a barotropic step.
    time_end_in       ! The end time for diagnostics when this routine started.
  real :: time_int_in ! The diagnostics' time interval when this routine started.
  real :: Htot_avg    ! The average total thickness of the tracer columns adjacent to a
                      ! velocity point [H ~> m or kg m-2]
  logical :: do_hifreq_output  ! If true, output occurs every barotropic step.
  logical :: use_BT_cont, do_ave, find_etaav, find_PF, find_Cor
  logical :: integral_BT_cont ! If true, update the barotropic continuity equation directly
                      ! from the initial condition using the time-integrated barotropic velocity.
  logical :: ice_is_rigid, nonblock_setup, interp_eta_PF
  logical :: project_velocity, add_uh0

  real :: dyn_coef_max ! The maximum stable value of dyn_coef_eta
                      ! [L2 T-2 H-1 ~> m s-2 or m4 s-2 kg-1].
  real :: ice_strength = 0.0  ! The effective strength of the ice [L2 Z-1 T-2 ~> m s-2].
  real :: Idt_max2    ! The squared inverse of the local maximum stable
                      ! barotropic time step [T-2 ~> s-2].
  real :: H_min_dyn   ! The minimum depth to use in limiting the size of the
                      ! dynamic surface pressure for stability [H ~> m or kg m-2].
  real :: H_eff_dx2   ! The effective total thickness divided by the grid spacing
                      ! squared [H L-2 ~> m-1 or kg m-4].
  real :: u_max_cor, v_max_cor ! The maximum corrective velocities [L T-1 ~> m s-1].
  real :: uint_cor, vint_cor ! The maximum time-integrated corrective velocities [L ~> m].
  real :: Htot        ! The total thickness [H ~> m or kg m-2].
  real :: eta_cor_max ! The maximum fluid that can be added as a correction to eta [H ~> m or kg m-2].
  real :: accel_underflow ! An acceleration that is so small it should be zeroed out [L T-2 ~> m s-2].
  real :: h_neglect            ! A thickness that is so small it is usually lost
                               ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: Idtbt       ! The inverse of the barotropic time step [T-1 ~> s-1]

  real, allocatable, dimension(:) :: wt_vel, wt_eta, wt_accel, wt_trans, wt_accel2
  real :: sum_wt_vel, sum_wt_eta, sum_wt_accel, sum_wt_trans
  real :: I_sum_wt_vel, I_sum_wt_eta, I_sum_wt_accel, I_sum_wt_trans
  real :: dt_filt     ! The half-width of the barotropic filter [T ~> s].
  real :: trans_wt1, trans_wt2 ! The weights used to compute ubt_trans and vbt_trans
  integer :: nfilter

  logical :: apply_OBCs, apply_OBC_flather, apply_OBC_open
  type(memory_size_type) :: MS
  character(len=200) :: mesg
  integer :: isv, iev, jsv, jev ! The valid array size at the end of a step.
  integer :: stencil  ! The stencil size of the algorithm, often 1 or 2.
  integer :: isvf, ievf, jsvf, jevf, num_cycles
  integer :: i, j, k, n
  integer :: is, ie, js, je, nz, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  integer :: ioff, joff
  integer :: l_seg

  if (.not.associated(CS)) call MOM_error(FATAL, &
      "btstep: Module MOM_barotropic must be initialized before it is used.")
  if (.not.CS%split) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  MS%isdw = CS%isdw ; MS%iedw = CS%iedw ; MS%jsdw = CS%jsdw ; MS%jedw = CS%jedw
  h_neglect = GV%H_subroundoff

  Idt = 1.0 / dt
  accel_underflow = CS%vel_underflow * Idt

  use_BT_cont = .false.
  if (present(BT_cont)) use_BT_cont = (associated(BT_cont))
  integral_BT_cont = use_BT_cont .and. CS%integral_BT_cont

  interp_eta_PF = .false.
  if (present(eta_PF_start)) interp_eta_PF = (associated(eta_PF_start))

  project_velocity = CS%BT_project_velocity

  ! Figure out the fullest arrays that could be updated.
  stencil = 1
  if ((.not.use_BT_cont) .and. CS%Nonlinear_continuity .and. &
      (CS%Nonlin_cont_update_period > 0)) stencil = 2

  do_ave = query_averaging_enabled(CS%diag)
  find_etaav = present(etaav)
  find_PF = (do_ave .and. ((CS%id_PFu_bt > 0) .or. (CS%id_PFv_bt > 0)))
  find_Cor = (do_ave .and. ((CS%id_Coru_bt > 0) .or. (CS%id_Corv_bt > 0)))

  add_uh0 = .false.
  if (present(uh0)) add_uh0 = associated(uh0)
  if (add_uh0 .and. .not.(present(vh0) .and. present(u_uh0) .and. &
                          present(v_vh0))) call MOM_error(FATAL, &
      "btstep: vh0, u_uh0, and v_vh0 must be present if uh0 is used.")
  if (add_uh0 .and. .not.(associated(vh0) .and. associated(u_uh0) .and. &
                          associated(v_vh0))) call MOM_error(FATAL, &
      "btstep: vh0, u_uh0, and v_vh0 must be associated if uh0 is used.")

  ! This can be changed to try to optimize the performance.
  nonblock_setup = G%nonblocking_updates

  if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)

  apply_OBCs = .false. ; CS%BT_OBC%apply_u_OBCs = .false. ; CS%BT_OBC%apply_v_OBCs = .false.
  apply_OBC_open = .false.
  apply_OBC_flather = .false.
  if (present(OBC)) then ; if (associated(OBC)) then
    CS%BT_OBC%apply_u_OBCs = OBC%open_u_BCs_exist_globally .or. OBC%specified_u_BCs_exist_globally
    CS%BT_OBC%apply_v_OBCs = OBC%open_v_BCs_exist_globally .or. OBC%specified_v_BCs_exist_globally
    apply_OBC_flather = open_boundary_query(OBC, apply_Flather_OBC=.true.)
    apply_OBC_open = open_boundary_query(OBC, apply_open_OBC=.true.)
    apply_OBCs = open_boundary_query(OBC, apply_specified_OBC=.true.) .or. &
           apply_OBC_flather .or. apply_OBC_open

    if (apply_OBC_flather .and. .not.GV%Boussinesq) call MOM_error(FATAL, &
      "btstep: Flather open boundary conditions have not yet been "// &
      "implemented for a non-Boussinesq model.")
  endif ; endif

  num_cycles = 1
  if (CS%use_wide_halos) &
    num_cycles = min((is-CS%isdw) / stencil, (js-CS%jsdw) / stencil)
  isvf = is - (num_cycles-1)*stencil ; ievf = ie + (num_cycles-1)*stencil
  jsvf = js - (num_cycles-1)*stencil ; jevf = je + (num_cycles-1)*stencil

  nstep = CEILING(dt/CS%dtbt - 0.0001)
  if (is_root_PE() .and. (nstep /= CS%nstep_last)) then
    write(mesg,'("btstep is using a dynamic barotropic timestep of ", ES12.6, &
               & " seconds, max ", ES12.6, ".")') (US%T_to_s*dt/nstep), US%T_to_s*CS%dtbt_max
    call MOM_mesg(mesg, 3)
  endif
  CS%nstep_last = nstep

  ! Set the actual barotropic time step.
  Instep = 1.0 / real(nstep)
  dtbt = dt * Instep
  Idtbt = 1.0 / dtbt
  bebt = CS%bebt
  be_proj = CS%bebt
  mass_accel_to_Z = 1.0 / GV%Rho0
  mass_to_Z = US%m_to_Z / GV%Rho0

  !--- setup the weight when computing vbt_trans and ubt_trans
  if (project_velocity) then
    trans_wt1 = (1.0 + be_proj); trans_wt2 = -be_proj
  else
    trans_wt1 = bebt ;           trans_wt2 = (1.0-bebt)
  endif

  do_hifreq_output = .false.
  if ((CS%id_ubt_hifreq > 0) .or. (CS%id_vbt_hifreq > 0) .or. &
      (CS%id_eta_hifreq > 0) .or. (CS%id_eta_pred_hifreq > 0) .or. &
      (CS%id_uhbt_hifreq > 0) .or. (CS%id_vhbt_hifreq > 0)) then
    do_hifreq_output = query_averaging_enabled(CS%diag, time_int_in, time_end_in)
    if (do_hifreq_output) &
      time_bt_start = time_end_in - real_to_time(US%T_to_s*dt)
  endif

!--- begin setup for group halo update
  if (id_clock_pass_pre > 0) call cpu_clock_begin(id_clock_pass_pre)
  if (.not. CS%linearized_BT_PV) then
    call create_group_pass(CS%pass_q_DCor, q, CS%BT_Domain, To_All, position=CORNER)
    call create_group_pass(CS%pass_q_DCor, DCor_u, DCor_v, CS%BT_Domain, &
         To_All+Scalar_Pair)
  endif
  if ((Isq > is-1) .or. (Jsq > js-1)) &
    call create_group_pass(CS%pass_tmp_uv, tmp_u, tmp_v, G%Domain)
  call create_group_pass(CS%pass_gtot, gtot_E, gtot_N, CS%BT_Domain, &
       To_All+Scalar_Pair, AGRID)
  call create_group_pass(CS%pass_gtot, gtot_W, gtot_S, CS%BT_Domain, &
       To_All+Scalar_Pair, AGRID)

  if (CS%dynamic_psurf) &
    call create_group_pass(CS%pass_eta_bt_rem, dyn_coef_eta, CS%BT_Domain)
  if (interp_eta_PF) then
    call create_group_pass(CS%pass_eta_bt_rem, eta_PF_1, CS%BT_Domain)
    call create_group_pass(CS%pass_eta_bt_rem, d_eta_PF, CS%BT_Domain)
  else
    call create_group_pass(CS%pass_eta_bt_rem, eta_PF, CS%BT_Domain)
  endif
  if (integral_BT_cont) &
    call create_group_pass(CS%pass_eta_bt_rem, eta_IC, CS%BT_Domain)
  call create_group_pass(CS%pass_eta_bt_rem, eta_src, CS%BT_Domain)
  ! The following halo updates are not needed without wide halos.  RWH
  ! We do need them after all.
! if (ievf > ie) then
    call create_group_pass(CS%pass_eta_bt_rem, bt_rem_u, bt_rem_v, &
                      CS%BT_Domain, To_All+Scalar_Pair)
    if (CS%linear_wave_drag) &
      call create_group_pass(CS%pass_eta_bt_rem, Rayleigh_u, Rayleigh_v, &
                      CS%BT_Domain, To_All+Scalar_Pair)
! endif
  ! The following halo update is not needed without wide halos.  RWH
  if (((G%isd > CS%isdw) .or. (G%jsd > CS%jsdw)) .or. (Isq <= is-1) .or. (Jsq <= js-1)) &
    call create_group_pass(CS%pass_force_hbt0_Cor_ref, BT_force_u, BT_force_v, CS%BT_Domain)
  if (add_uh0) call create_group_pass(CS%pass_force_hbt0_Cor_ref, uhbt0, vhbt0, CS%BT_Domain)
  call create_group_pass(CS%pass_force_hbt0_Cor_ref, Cor_ref_u, Cor_ref_v, CS%BT_Domain)
  if (.not. use_BT_cont) then
    call create_group_pass(CS%pass_Dat_uv, Datu, Datv, CS%BT_Domain, To_All+Scalar_Pair)
  endif
  call create_group_pass(CS%pass_eta_ubt, eta, CS%BT_Domain)
  call create_group_pass(CS%pass_eta_ubt, ubt, vbt, CS%BT_Domain)
  if (integral_BT_cont) then
    call create_group_pass(CS%pass_eta_ubt, ubt_int, vbt_int, CS%BT_Domain)
    ! This is only needed with integral_BT_cont, OBCs and multiple barotropic steps between halo updates.
    if (apply_OBC_open) &
      call create_group_pass(CS%pass_eta_ubt, uhbt_int, vhbt_int, CS%BT_Domain)
  endif

  call create_group_pass(CS%pass_ubt_Cor, ubt_Cor, vbt_Cor, G%Domain)
  ! These passes occur at the end of the routine, as data is being readied to
  ! share with the main part of the MOM6 code.
  if (find_etaav) then
    call create_group_pass(CS%pass_etaav, etaav, G%Domain)
  endif
  call create_group_pass(CS%pass_e_anom, e_anom, G%Domain)
  call create_group_pass(CS%pass_ubta_uhbta, CS%ubtav, CS%vbtav, G%Domain)
  call create_group_pass(CS%pass_ubta_uhbta, uhbtav, vhbtav, G%Domain)

  if (id_clock_pass_pre > 0) call cpu_clock_end(id_clock_pass_pre)
!--- end setup for group halo update

!   Calculate the constant coefficients for the Coriolis force terms in the
! barotropic momentum equations.  This has to be done quite early to start
! the halo update that needs to be completed before the next calculations.
  if (CS%linearized_BT_PV) then
    !$OMP parallel do default(shared)
    do J=jsvf-2,jevf+1 ; do I=isvf-2,ievf+1
      q(I,J) = CS%q_D(I,j)
    enddo ; enddo
    !$OMP parallel do default(shared)
    do j=jsvf-1,jevf+1 ; do I=isvf-2,ievf+1
      DCor_u(I,j) = CS%D_u_Cor(I,j)
    enddo ; enddo
    !$OMP parallel do default(shared)
    do J=jsvf-2,jevf+1 ; do i=isvf-1,ievf+1
      DCor_v(i,J) = CS%D_v_Cor(i,J)
    enddo ; enddo
  else
    q(:,:) = 0.0 ; DCor_u(:,:) = 0.0 ; DCor_v(:,:) = 0.0
    if (GV%Boussinesq) then
      !$OMP parallel do default(shared)
      do j=js,je ; do I=is-1,ie
        DCor_u(I,j) = 0.5 * (max(GV%Z_to_H*G%bathyT(i+1,j) + eta_in(i+1,j), 0.0) + &
                             max(GV%Z_to_H*G%bathyT(i,j) + eta_in(i,j), 0.0) )
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=js-1,je ; do i=is,ie
        DCor_v(i,J) = 0.5 * (max(GV%Z_to_H*G%bathyT(i,j+1) + eta_in(i+1,j), 0.0) + &
                             max(GV%Z_to_H*G%bathyT(i,j) + eta_in(i,j), 0.0) )
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=js-1,je ; do I=is-1,ie
        q(I,J) = 0.25 * (CS%BT_Coriolis_scale * G%CoriolisBu(I,J)) * &
             ((G%areaT(i,j) + G%areaT(i+1,j+1)) + (G%areaT(i+1,j) + G%areaT(i,j+1))) / &
             (max((G%areaT(i,j) * max(GV%Z_to_H*G%bathyT(i,j) + eta_in(i,j), 0.0) + &
               G%areaT(i+1,j+1) * max(GV%Z_to_H*G%bathyT(i+1,j+1) + eta_in(i+1,j+1), 0.0)) + &
              (G%areaT(i+1,j) * max(GV%Z_to_H*G%bathyT(i+1,j) + eta_in(i+1,j), 0.0) + &
               G%areaT(i,j+1) * max(GV%Z_to_H*G%bathyT(i,j+1) + eta_in(i,j+1), 0.0)), h_neglect) )
      enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=js,je ; do I=is-1,ie
        DCor_u(I,j) = 0.5 * (eta_in(i+1,j) + eta_in(i,j))
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=js-1,je ; do i=is,ie
        DCor_v(i,J) = 0.5 * (eta_in(i,j+1) + eta_in(i,j))
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=js-1,je ; do I=is-1,ie
        q(I,J) = 0.25 * (CS%BT_Coriolis_scale * G%CoriolisBu(I,J)) * &
             ((G%areaT(i,j) + G%areaT(i+1,j+1)) + (G%areaT(i+1,j) + G%areaT(i,j+1))) / &
             (max((G%areaT(i,j) * eta_in(i,j) + G%areaT(i+1,j+1) * eta_in(i+1,j+1)) + &
                  (G%areaT(i+1,j) * eta_in(i+1,j) + G%areaT(i,j+1) * eta_in(i,j+1)), h_neglect) )
      enddo ; enddo
    endif

    ! With very wide halos, q and D need to be calculated on the available data
    ! domain and then updated onto the full computational domain.
    ! These calculations can be done almost immediately, but the halo updates
    ! must be done before the [abcd]mer and [abcd]zon are calculated.
    if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
    if (nonblock_setup) then
      call start_group_pass(CS%pass_q_DCor, CS%BT_Domain, clock=id_clock_pass_pre)
    else
      call do_group_pass(CS%pass_q_DCor, CS%BT_Domain, clock=id_clock_pass_pre)
    endif
    if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)
  endif

  ! Zero out various wide-halo arrays.
  !$OMP parallel do default(shared)
  do j=CS%jsdw,CS%jedw ; do i=CS%isdw,CS%iedw
    gtot_E(i,j) = 0.0 ; gtot_W(i,j) = 0.0
    gtot_N(i,j) = 0.0 ; gtot_S(i,j) = 0.0
    eta(i,j) = 0.0
    eta_PF(i,j) = 0.0
    if (interp_eta_PF) then
      eta_PF_1(i,j) = 0.0 ; d_eta_PF(i,j) = 0.0
    endif
    if (integral_BT_cont) then
      eta_IC(i,j) = 0.0
    endif
    p_surf_dyn(i,j) = 0.0
    if (CS%dynamic_psurf) dyn_coef_eta(i,j) = 0.0
  enddo ; enddo
  !   The halo regions of various arrays need to be initialized to
  ! non-NaNs in case the neighboring domains are not part of the ocean.
  ! Otherwise a halo update later on fills in the correct values.
  !$OMP parallel do default(shared)
  do j=CS%jsdw,CS%jedw ; do I=CS%isdw-1,CS%iedw
    Cor_ref_u(I,j) = 0.0 ; BT_force_u(I,j) = 0.0 ; ubt(I,j) = 0.0
    Datu(I,j) = 0.0 ; bt_rem_u(I,j) = 0.0 ; uhbt0(I,j) = 0.0
  enddo ; enddo
  !$OMP parallel do default(shared)
  do J=CS%jsdw-1,CS%jedw ; do i=CS%isdw,CS%iedw
    Cor_ref_v(i,J) = 0.0 ; BT_force_v(i,J) = 0.0 ; vbt(i,J) = 0.0
    Datv(i,J) = 0.0 ; bt_rem_v(i,J) = 0.0 ; vhbt0(i,J) = 0.0
  enddo ; enddo

  if (CS%linear_wave_drag) then
    !$OMP parallel do default(shared)
    do j=CS%jsdw,CS%jedw ; do I=CS%isdw-1,CS%iedw
      Rayleigh_u(I,j) = 0.0
    enddo ; enddo
    !$OMP parallel do default(shared)
    do J=CS%jsdw-1,CS%jedw ; do i=CS%isdw,CS%iedw
      Rayleigh_v(i,J) = 0.0
    enddo ; enddo
  endif

  ! Copy input arrays into their wide-halo counterparts.
  if (interp_eta_PF) then
    !$OMP parallel do default(shared)
    do j=G%jsd,G%jed ; do i=G%isd,G%ied ! Was "do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1" but doing so breaks OBC. Not sure why?
      eta(i,j) = eta_in(i,j)
      eta_PF_1(i,j) = eta_PF_start(i,j)
      d_eta_PF(i,j) = eta_PF_in(i,j) - eta_PF_start(i,j)
    enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do j=G%jsd,G%jed ; do i=G%isd,G%ied !: Was "do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1" but doing so breaks OBC. Not sure why?
      eta(i,j) = eta_in(i,j)
      eta_PF(i,j) = eta_PF_in(i,j)
    enddo ; enddo
  endif
  if (integral_BT_cont) then
    !$OMP parallel do default(shared)
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      eta_IC(i,j) = eta_in(i,j)
    enddo ; enddo
  endif

  !$OMP parallel do default(shared) private(visc_rem)
  do k=1,nz ; do j=js,je ; do I=is-1,ie
    ! rem needs greater than visc_rem_u and 1-Instep/visc_rem_u.
    ! The 0.5 below is just for safety.
    if (visc_rem_u(I,j,k) <= 0.0) then ; visc_rem = 0.0
    elseif (visc_rem_u(I,j,k) >= 1.0) then ; visc_rem = 1.0
    elseif (visc_rem_u(I,j,k)**2 > visc_rem_u(I,j,k) - 0.5*Instep) then
      visc_rem = visc_rem_u(I,j,k)
    else ; visc_rem = 1.0 - 0.5*Instep/visc_rem_u(I,j,k) ; endif
    wt_u(I,j,k) = CS%frhatu(I,j,k) * visc_rem
  enddo ; enddo ; enddo
  !$OMP parallel do default(shared) private(visc_rem)
  do k=1,nz ; do J=js-1,je ; do i=is,ie
    ! rem needs greater than visc_rem_v and 1-Instep/visc_rem_v.
    if (visc_rem_v(i,J,k) <= 0.0) then ; visc_rem = 0.0
    elseif (visc_rem_v(i,J,k) >= 1.0) then ; visc_rem = 1.0
    elseif (visc_rem_v(i,J,k)**2 > visc_rem_v(i,J,k) - 0.5*Instep) then
      visc_rem = visc_rem_v(i,J,k)
    else ; visc_rem = 1.0 - 0.5*Instep/visc_rem_v(i,J,k) ; endif
    wt_v(i,J,k) = CS%frhatv(i,J,k) * visc_rem
  enddo ; enddo ; enddo

  !   Use u_Cor and v_Cor as the reference values for the Coriolis terms,
  ! including the viscous remnant.
  !$OMP parallel do default(shared)
  do j=js-1,je+1 ; do I=is-1,ie ; ubt_Cor(I,j) = 0.0 ; enddo ; enddo
  !$OMP parallel do default(shared)
  do J=js-1,je ; do i=is-1,ie+1 ; vbt_Cor(i,J) = 0.0 ; enddo ; enddo
  !$OMP parallel do default(shared)
  do j=js,je ; do k=1,nz ; do I=is-1,ie
    ubt_Cor(I,j) = ubt_Cor(I,j) + wt_u(I,j,k) * U_Cor(I,j,k)
  enddo ; enddo ; enddo
  !$OMP parallel do default(shared)
  do J=js-1,je ; do k=1,nz ; do i=is,ie
    vbt_Cor(i,J) = vbt_Cor(i,J) + wt_v(i,J,k) * V_Cor(i,J,k)
  enddo ; enddo ; enddo

  ! The gtot arrays are the effective layer-weighted reduced gravities for
  ! accelerations across the various faces, with names for the relative
  ! locations of the faces to the pressure point.  They will have their halos
  ! updated later on.
  !$OMP parallel do default(shared)
  do j=js,je
    do k=1,nz ; do I=is-1,ie
      gtot_E(i,j)   = gtot_E(i,j)   + pbce(i,j,k)   * wt_u(I,j,k)
      gtot_W(i+1,j) = gtot_W(i+1,j) + pbce(i+1,j,k) * wt_u(I,j,k)
    enddo ; enddo
  enddo
  !$OMP parallel do default(shared)
  do J=js-1,je
    do k=1,nz ; do i=is,ie
      gtot_N(i,j)   = gtot_N(i,j)   + pbce(i,j,k)   * wt_v(i,J,k)
      gtot_S(i,j+1) = gtot_S(i,j+1) + pbce(i,j+1,k) * wt_v(i,J,k)
    enddo ; enddo
  enddo

  if (CS%tides) then
    call tidal_forcing_sensitivity(G, CS%tides_CSp, det_de)
    dgeo_de = 1.0 + det_de + CS%G_extra
  else
    dgeo_de = 1.0 + CS%G_extra
  endif

  if (nonblock_setup .and. .not.CS%linearized_BT_PV) then
    if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
    call complete_group_pass(CS%pass_q_DCor, CS%BT_Domain, clock=id_clock_pass_pre)
    if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)
  endif

  ! Calculate the open areas at the velocity points.
  ! The halo updates are needed before Datu is first used, either in set_up_BT_OBC or ubt_Cor.
  if (integral_BT_cont) then
    call set_local_BT_cont_types(BT_cont, BTCL_u, BTCL_v, G, US, MS, CS%BT_Domain, 1+ievf-ie, dt_baroclinic=dt)
  elseif (use_BT_cont) then
    call set_local_BT_cont_types(BT_cont, BTCL_u, BTCL_v, G, US, MS, CS%BT_Domain, 1+ievf-ie)
  else
    if (CS%Nonlinear_continuity) then
      call find_face_areas(Datu, Datv, G, GV, US, CS, MS, eta, 1)
    else
      call find_face_areas(Datu, Datv, G, GV, US, CS, MS, halo=1)
    endif
  endif

  ! Set up fields related to the open boundary conditions.
  if (apply_OBCs) then
    call set_up_BT_OBC(OBC, eta, CS%BT_OBC, CS%BT_Domain, G, GV, US, MS, ievf-ie, use_BT_cont, &
                       integral_BT_cont, dt, Datu, Datv, BTCL_u, BTCL_v)
  endif

  ! Determine the difference between the sum of the layer fluxes and the
  ! barotropic fluxes found from the same input velocities.
  if (add_uh0) then
    !$OMP parallel do default(shared)
    do j=js,je ; do I=is-1,ie ; uhbt(I,j) = 0.0 ; ubt(I,j) = 0.0 ; enddo ; enddo
    !$OMP parallel do default(shared)
    do J=js-1,je ; do i=is,ie ; vhbt(i,J) = 0.0 ; vbt(i,J) = 0.0 ; enddo ; enddo
    if (CS%visc_rem_u_uh0) then
      !$OMP parallel do default(shared)
      do j=js,je ; do k=1,nz ; do I=is-1,ie
        uhbt(I,j) = uhbt(I,j) + uh0(I,j,k)
        ubt(I,j) = ubt(I,j) + wt_u(I,j,k) * u_uh0(I,j,k)
      enddo ; enddo ; enddo
      !$OMP parallel do default(shared)
      do J=js-1,je ; do k=1,nz ; do i=is,ie
        vhbt(i,J) = vhbt(i,J) + vh0(i,J,k)
        vbt(i,J) = vbt(i,J) + wt_v(i,J,k) * v_vh0(i,J,k)
      enddo ; enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=js,je ; do k=1,nz ; do I=is-1,ie
        uhbt(I,j) = uhbt(I,j) + uh0(I,j,k)
        ubt(I,j) = ubt(I,j) + CS%frhatu(I,j,k) * u_uh0(I,j,k)
      enddo ; enddo ; enddo
      !$OMP parallel do default(shared)
      do J=js-1,je ; do k=1,nz ; do i=is,ie
        vhbt(i,J) = vhbt(i,J) + vh0(i,J,k)
        vbt(i,J) = vbt(i,J) + CS%frhatv(i,J,k) * v_vh0(i,J,k)
      enddo ; enddo ; enddo
    endif
    if ((use_BT_cont .or. integral_BT_cont) .and. CS%adjust_BT_cont) then
      ! Use the additional input transports to broaden the fits
      ! over which the bt_cont_type applies.

      ! Fill in the halo data for ubt, vbt, uhbt, and vhbt.
      if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
      if (id_clock_pass_pre > 0) call cpu_clock_begin(id_clock_pass_pre)
      call pass_vector(ubt, vbt, CS%BT_Domain, complete=.false., halo=1+ievf-ie)
      call pass_vector(uhbt, vhbt, CS%BT_Domain, complete=.true., halo=1+ievf-ie)
      if (id_clock_pass_pre > 0) call cpu_clock_end(id_clock_pass_pre)
      if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)

      if (integral_BT_cont) then
        call adjust_local_BT_cont_types(ubt, uhbt, vbt, vhbt, BTCL_u, BTCL_v, &
                                        G, US, MS, halo=1+ievf-ie, dt_baroclinic=dt)
      else
        call adjust_local_BT_cont_types(ubt, uhbt, vbt, vhbt, BTCL_u, BTCL_v, &
                                        G, US, MS, halo=1+ievf-ie)
      endif
    endif
    if (integral_BT_cont) then
      !$OMP parallel do default(shared)
      do j=js,je ; do I=is-1,ie
        uhbt0(I,j) = uhbt(I,j) - find_uhbt(dt*ubt(I,j), BTCL_u(I,j)) * Idt
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=js-1,je ; do i=is,ie
        vhbt0(i,J) = vhbt(i,J) - find_vhbt(dt*vbt(i,J), BTCL_v(i,J)) * Idt
      enddo ; enddo
    elseif (use_BT_cont) then
      !$OMP parallel do default(shared)
      do j=js,je ; do I=is-1,ie
        uhbt0(I,j) = uhbt(I,j) - find_uhbt(ubt(I,j), BTCL_u(I,j))
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=js-1,je ; do i=is,ie
        vhbt0(i,J) = vhbt(i,J) - find_vhbt(vbt(i,J), BTCL_v(i,J))
      enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=js,je ; do I=is-1,ie
        uhbt0(I,j) = uhbt(I,j) - Datu(I,j)*ubt(I,j)
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=js-1,je ; do i=is,ie
        vhbt0(i,J) = vhbt(i,J) - Datv(i,J)*vbt(i,J)
      enddo ; enddo
    endif
    if (CS%BT_OBC%apply_u_OBCs) then  ! zero out pressure force across boundary
      !$OMP parallel do default(shared)
      do j=js,je ; do I=is-1,ie ; if (OBC%segnum_u(I,j) /= OBC_NONE) then
        uhbt0(I,j) = 0.0
      endif ; enddo ; enddo
    endif
    if (CS%BT_OBC%apply_v_OBCs) then  ! zero out PF across boundary
      !$OMP parallel do default(shared)
      do J=js-1,je ; do i=is,ie ; if (OBC%segnum_v(i,J) /= OBC_NONE) then
        vhbt0(i,J) = 0.0
      endif ; enddo ; enddo
    endif
  endif

! Calculate the initial barotropic velocities from the layer's velocities.
  if (integral_BT_cont) then
    !$OMP parallel do default(shared)
    do j=jsvf-1,jevf+1 ; do I=isvf-2,ievf+1
      ubt(I,j) = 0.0 ; uhbt(I,j) = 0.0 ; u_accel_bt(I,j) = 0.0
      ubt_int(I,j) = 0.0 ; uhbt_int(I,j) = 0.0
    enddo ; enddo
    !$OMP parallel do default(shared)
    do J=jsvf-2,jevf+1 ; do i=isvf-1,ievf+1
      vbt(i,J) = 0.0 ; vhbt(i,J) = 0.0 ; v_accel_bt(i,J) = 0.0
      vbt_int(i,J) = 0.0 ; vhbt_int(i,J) = 0.0
    enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do j=jsvf-1,jevf+1 ; do I=isvf-2,ievf+1
      ubt(I,j) = 0.0 ; uhbt(I,j) = 0.0 ; u_accel_bt(I,j) = 0.0
    enddo ; enddo
    !$OMP parallel do default(shared)
    do J=jsvf-2,jevf+1 ; do i=isvf-1,ievf+1
      vbt(i,J) = 0.0 ; vhbt(i,J) = 0.0 ; v_accel_bt(i,J) = 0.0
    enddo ; enddo
  endif
  !$OMP parallel do default(shared)
  do j=js,je ; do k=1,nz ; do I=is-1,ie
    ubt(I,j) = ubt(I,j) + wt_u(I,j,k) * U_in(I,j,k)
  enddo ; enddo ; enddo
  !$OMP parallel do default(shared)
  do J=js-1,je ; do k=1,nz ; do i=is,ie
    vbt(i,J) = vbt(i,J) + wt_v(i,J,k) * V_in(i,J,k)
  enddo ; enddo ;  enddo
  !$OMP parallel do default(shared)
  do j=js,je ; do I=is-1,ie
    if (abs(ubt(I,j)) < CS%vel_underflow) ubt(I,j) = 0.0
  enddo ; enddo
  !$OMP parallel do default(shared)
  do J=js-1,je ; do i=is,ie
    if (abs(vbt(i,J)) < CS%vel_underflow) vbt(i,J) = 0.0
  enddo ; enddo

  if (apply_OBCs) then
    ubt_first(:,:) = ubt(:,:) ; vbt_first(:,:) = vbt(:,:)
  endif

!   Here the vertical average accelerations due to the Coriolis, advective,
! pressure gradient and horizontal viscous terms in the layer momentum
! equations are calculated.  These will be used to determine the difference
! between the accelerations due to the average of the layer equations and the
! barotropic calculation.

  !$OMP parallel do default(shared)
  do j=js,je ; do I=is-1,ie ; if (G%mask2dCu(I,j) > 0.0) then
    if (CS%nonlin_stress) then
      if (GV%Boussinesq) then
        Htot_avg = 0.5*(max(CS%bathyT(i,j)*GV%Z_to_H + eta(i,j), 0.0) + &
                        max(CS%bathyT(i+1,j)*GV%Z_to_H + eta(i+1,j), 0.0))
      else
        Htot_avg = 0.5*(eta(i,j) + eta(i+1,j))
      endif
      if (Htot_avg*CS%dy_Cu(I,j) <= 0.0) then
        CS%IDatu(I,j) = 0.0
      elseif (integral_BT_cont) then
        CS%IDatu(I,j) = CS%dy_Cu(I,j) / (max(find_duhbt_du(ubt(I,j)*dt, BTCL_u(I,j)), &
                                             CS%dy_Cu(I,j)*Htot_avg) )
      elseif (use_BT_cont) then ! Reconsider the max and whether there should be some scaling.
        CS%IDatu(I,j) = CS%dy_Cu(I,j) / (max(find_duhbt_du(ubt(I,j), BTCL_u(I,j)), &
                                             CS%dy_Cu(I,j)*Htot_avg) )
      else
        CS%IDatu(I,j) = 1.0 / Htot_avg
      endif
    endif

    BT_force_u(I,j) = forces%taux(I,j) * mass_accel_to_Z * CS%IDatu(I,j)*visc_rem_u(I,j,1)
  else
    BT_force_u(I,j) = 0.0
  endif ; enddo ; enddo
  !$OMP parallel do default(shared)
  do J=js-1,je ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.0) then
    if (CS%nonlin_stress) then
      if (GV%Boussinesq) then
        Htot_avg = 0.5*(max(CS%bathyT(i,j)*GV%Z_to_H + eta(i,j), 0.0) + &
                        max(CS%bathyT(i,j+1)*GV%Z_to_H + eta(i,j+1), 0.0))
      else
        Htot_avg = 0.5*(eta(i,j) + eta(i,j+1))
      endif
      if (Htot_avg*CS%dx_Cv(i,J) <= 0.0) then
        CS%IDatv(i,J) = 0.0
      elseif (integral_BT_cont) then
        CS%IDatv(i,J) = CS%dx_Cv(i,J) / (max(find_dvhbt_dv(vbt(i,J)*dt, BTCL_v(i,J)), &
                                             CS%dx_Cv(i,J)*Htot_avg) )
      elseif (use_BT_cont) then ! Reconsider the max and whether there should be some scaling.
        CS%IDatv(i,J) = CS%dx_Cv(i,J) / (max(find_dvhbt_dv(vbt(i,J), BTCL_v(i,J)), &
                                             CS%dx_Cv(i,J)*Htot_avg) )
      else
        CS%IDatv(i,J) = 1.0 / Htot_avg
      endif
    endif

    BT_force_v(i,J) = forces%tauy(i,J) * mass_accel_to_Z * CS%IDatv(i,J)*visc_rem_v(i,J,1)
  else
    BT_force_v(i,J) = 0.0
  endif ; enddo ; enddo
  if (present(taux_bot) .and. present(tauy_bot)) then
    if (associated(taux_bot) .and. associated(tauy_bot)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do I=is-1,ie ; if (G%mask2dCu(I,j) > 0.0) then
        BT_force_u(I,j) = BT_force_u(I,j) - taux_bot(I,j) * mass_to_Z  * CS%IDatu(I,j)
      endif ; enddo ; enddo
      !$OMP parallel do default(shared)
      do J=js-1,je ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.0) then
        BT_force_v(i,J) = BT_force_v(i,J) - tauy_bot(i,J) * mass_to_Z  * CS%IDatv(i,J)
      endif ; enddo ; enddo
    endif
  endif

  ! bc_accel_u & bc_accel_v are only available on the potentially
  ! non-symmetric computational domain.
  !$OMP parallel do default(shared)
  do j=js,je ; do k=1,nz ; do I=Isq,Ieq
    BT_force_u(I,j) = BT_force_u(I,j) + wt_u(I,j,k) * bc_accel_u(I,j,k)
  enddo ; enddo ; enddo
  !$OMP parallel do default(shared)
  do J=Jsq,Jeq ; do k=1,nz ; do i=is,ie
    BT_force_v(i,J) = BT_force_v(i,J) + wt_v(i,J,k) * bc_accel_v(i,J,k)
  enddo ; enddo ; enddo

  if (CS%gradual_BT_ICs) then
    !$OMP parallel do default(shared)
    do j=js,je ; do I=is-1,ie
      BT_force_u(I,j) = BT_force_u(I,j) + (ubt(I,j) - CS%ubt_IC(I,j)) * Idt
      ubt(I,j) = CS%ubt_IC(I,j)
      if (abs(ubt(I,j)) < CS%vel_underflow) ubt(I,j) = 0.0
    enddo ; enddo
    !$OMP parallel do default(shared)
    do J=js-1,je ; do i=is,ie
      BT_force_v(i,J) = BT_force_v(i,J) + (vbt(i,J) - CS%vbt_IC(i,J)) * Idt
      vbt(i,J) = CS%vbt_IC(i,J)
      if (abs(vbt(i,J)) < CS%vel_underflow) vbt(i,J) = 0.0
    enddo ; enddo
  endif

  if ((Isq > is-1) .or. (Jsq > js-1)) then
    ! Non-symmetric memory is being used, so the edge values need to be
    ! filled in with a halo update of a non-symmetric array.
    if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
    if (id_clock_pass_pre > 0) call cpu_clock_begin(id_clock_pass_pre)
    tmp_u(:,:) = 0.0 ; tmp_v(:,:) = 0.0
    do j=js,je ; do I=Isq,Ieq ; tmp_u(I,j) = BT_force_u(I,j) ; enddo ; enddo
    do J=Jsq,Jeq ; do i=is,ie ; tmp_v(i,J) = BT_force_v(i,J) ; enddo ; enddo
    if (nonblock_setup) then
      call start_group_pass(CS%pass_tmp_uv, G%Domain)
    else
      call do_group_pass(CS%pass_tmp_uv, G%Domain)
      do j=jsd,jed ; do I=IsdB,IedB ; BT_force_u(I,j) = tmp_u(I,j) ; enddo ; enddo
      do J=JsdB,JedB ; do i=isd,ied ; BT_force_v(i,J) = tmp_v(i,J) ; enddo ; enddo
    endif
    if (id_clock_pass_pre > 0) call cpu_clock_end(id_clock_pass_pre)
    if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)
  endif

  if (nonblock_setup) then
    if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
    if (id_clock_pass_pre > 0) call cpu_clock_begin(id_clock_pass_pre)
    call start_group_pass(CS%pass_gtot, CS%BT_Domain)
    call start_group_pass(CS%pass_ubt_Cor, G%Domain)
    if (id_clock_pass_pre > 0) call cpu_clock_end(id_clock_pass_pre)
    if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)
  endif

  ! Determine the weighted Coriolis parameters for the neighboring velocities.
  !$OMP parallel do default(shared)
  do J=jsvf-1,jevf ; do i=isvf-1,ievf+1
    if (CS%Sadourny) then
      amer(I-1,j) = DCor_u(I-1,j) * q(I-1,J)
      bmer(I,j) = DCor_u(I,j) * q(I,J)
      cmer(I,j+1) = DCor_u(I,j+1) * q(I,J)
      dmer(I-1,j+1) = DCor_u(I-1,j+1) * q(I-1,J)
    else
      amer(I-1,j) = DCor_u(I-1,j) * &
                    ((q(I,J) + q(I-1,J-1)) + q(I-1,J)) / 3.0
      bmer(I,j) = DCor_u(I,j) * &
                  (q(I,J) + (q(I-1,J) + q(I,J-1))) / 3.0
      cmer(I,j+1) = DCor_u(I,j+1) * &
                    (q(I,J) + (q(I-1,J) + q(I,J+1))) / 3.0
      dmer(I-1,j+1) = DCor_u(I-1,j+1) * &
                      ((q(I,J) + q(I-1,J+1)) + q(I-1,J)) / 3.0
    endif
  enddo ; enddo

  !$OMP parallel do default(shared)
  do j=jsvf-1,jevf+1 ; do I=isvf-1,ievf
    if (CS%Sadourny) then
      azon(I,j) = DCor_v(i+1,J) * q(I,J)
      bzon(I,j) = DCor_v(i,J) * q(I,J)
      czon(I,j) = DCor_v(i,J-1) * q(I,J-1)
      dzon(I,j) = DCor_v(i+1,J-1) * q(I,J-1)
    else
      azon(I,j) = DCor_v(i+1,J) * &
                  (q(I,J) + (q(I+1,J) + q(I,J-1))) / 3.0
      bzon(I,j) = DCor_v(i,J) * &
                  (q(I,J) + (q(I-1,J) + q(I,J-1))) / 3.0
      czon(I,j) = DCor_v(i,J-1) * &
                  ((q(I,J) + q(I-1,J-1)) + q(I,J-1)) / 3.0
      dzon(I,j) = DCor_v(i+1,J-1) * &
                  ((q(I,J) + q(I+1,J-1)) + q(I,J-1)) / 3.0
    endif
  enddo ; enddo

! Complete the previously initiated message passing.
  if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
  if (id_clock_pass_pre > 0) call cpu_clock_begin(id_clock_pass_pre)
  if (nonblock_setup) then
    if ((Isq > is-1) .or. (Jsq > js-1)) then
      call complete_group_pass(CS%pass_tmp_uv, G%Domain)
      do j=jsd,jed ; do I=IsdB,IedB ; BT_force_u(I,j) = tmp_u(I,j) ; enddo ; enddo
      do J=JsdB,JedB ; do i=isd,ied ; BT_force_v(i,J) = tmp_v(i,J) ; enddo ; enddo
    endif
    call complete_group_pass(CS%pass_gtot, CS%BT_Domain)
    call complete_group_pass(CS%pass_ubt_Cor, G%Domain)
  else
    call do_group_pass(CS%pass_gtot, CS%BT_Domain)
    call do_group_pass(CS%pass_ubt_Cor, G%Domain)
  endif
  ! The various elements of gtot are positive definite but directional, so use
  ! the polarity arrays to sort out when the directions have shifted.
  do j=jsvf-1,jevf+1 ; do i=isvf-1,ievf+1
    if (CS%ua_polarity(i,j) < 0.0) call swap(gtot_E(i,j), gtot_W(i,j))
    if (CS%va_polarity(i,j) < 0.0) call swap(gtot_N(i,j), gtot_S(i,j))
  enddo ; enddo

  !$OMP parallel do default(shared)
  do j=js,je ; do I=is-1,ie
    Cor_ref_u(I,j) =  &
        ((azon(I,j) * vbt_Cor(i+1,j) + czon(I,j) * vbt_Cor(i  ,j-1)) + &
         (bzon(I,j) * vbt_Cor(i  ,j) + dzon(I,j) * vbt_Cor(i+1,j-1)))
  enddo ; enddo
  !$OMP parallel do default(shared)
  do J=js-1,je ; do i=is,ie
    Cor_ref_v(i,J) = -1.0 * &
        ((amer(I-1,j) * ubt_Cor(I-1,j) + cmer(I  ,j+1) * ubt_Cor(I  ,j+1)) + &
         (bmer(I  ,j) * ubt_Cor(I  ,j) + dmer(I-1,j+1) * ubt_Cor(I-1,j+1)))
  enddo ; enddo

  ! Now start new halo updates.
  if (nonblock_setup) then
    if (.not.use_BT_cont) &
      call start_group_pass(CS%pass_Dat_uv, CS%BT_Domain)

    ! The following halo update is not needed without wide halos.  RWH
    call start_group_pass(CS%pass_force_hbt0_Cor_ref, CS%BT_Domain)
  endif
  if (id_clock_pass_pre > 0) call cpu_clock_end(id_clock_pass_pre)
  if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)
  !$OMP parallel default(shared) private(u_max_cor,uint_cor,v_max_cor,vint_cor,eta_cor_max,Htot)
  !$OMP do
  do j=js-1,je+1 ; do I=is-1,ie ; av_rem_u(I,j) = 0.0 ; enddo ; enddo
  !$OMP do
  do J=js-1,je ; do i=is-1,ie+1 ; av_rem_v(i,J) = 0.0 ; enddo ; enddo
  !$OMP do
  do j=js,je ; do k=1,nz ; do I=is-1,ie
    av_rem_u(I,j) = av_rem_u(I,j) + CS%frhatu(I,j,k) * visc_rem_u(I,j,k)
  enddo ; enddo ; enddo
  !$OMP do
  do J=js-1,je ; do k=1,nz ; do i=is,ie
    av_rem_v(i,J) = av_rem_v(i,J) + CS%frhatv(i,J,k) * visc_rem_v(i,J,k)
  enddo ; enddo ; enddo
  if (CS%strong_drag) then
    !$OMP do
    do j=js,je ; do I=is-1,ie
      bt_rem_u(I,j) = G%mask2dCu(I,j) * &
         ((nstep * av_rem_u(I,j)) / (1.0 + (nstep-1)*av_rem_u(I,j)))
    enddo ; enddo
    !$OMP do
    do J=js-1,je ; do i=is,ie
      bt_rem_v(i,J) = G%mask2dCv(i,J) * &
         ((nstep * av_rem_v(i,J)) / (1.0 + (nstep-1)*av_rem_v(i,J)))
    enddo ; enddo
  else
    !$OMP do
    do j=js,je ; do I=is-1,ie
      bt_rem_u(I,j) = 0.0
      if (G%mask2dCu(I,j) * av_rem_u(I,j) > 0.0) &
        bt_rem_u(I,j) = G%mask2dCu(I,j) * (av_rem_u(I,j)**Instep)
    enddo ; enddo
    !$OMP do
    do J=js-1,je ; do i=is,ie
      bt_rem_v(i,J) = 0.0
      if (G%mask2dCv(i,J) * av_rem_v(i,J) > 0.0) &
        bt_rem_v(i,J) = G%mask2dCv(i,J) * (av_rem_v(i,J)**Instep)
    enddo ; enddo
  endif
  if (CS%linear_wave_drag) then
    !$OMP do
    do j=js,je ; do I=is-1,ie ; if (CS%lin_drag_u(I,j) > 0.0) then
      Htot = 0.5 * (eta(i,j) + eta(i+1,j))
      if (GV%Boussinesq) &
        Htot = Htot + 0.5*GV%Z_to_H * (CS%bathyT(i,j) + CS%bathyT(i+1,j))
      bt_rem_u(I,j) = bt_rem_u(I,j) * (Htot / (Htot + CS%lin_drag_u(I,j) * dtbt))

      Rayleigh_u(I,j) = CS%lin_drag_u(I,j) / Htot
    endif ; enddo ; enddo
    !$OMP do
    do J=js-1,je ; do i=is,ie ; if (CS%lin_drag_v(i,J) > 0.0) then
      Htot = 0.5 * (eta(i,j) + eta(i,j+1))
      if (GV%Boussinesq) &
        Htot = Htot + 0.5*GV%Z_to_H * (CS%bathyT(i,j) + CS%bathyT(i,j+1))
      bt_rem_v(i,J) = bt_rem_v(i,J) * (Htot / (Htot + CS%lin_drag_v(i,J) * dtbt))

      Rayleigh_v(i,J) = CS%lin_drag_v(i,J) / Htot
    endif ; enddo ; enddo
  endif

  ! Zero out the arrays for various time-averaged quantities.
  if (find_etaav) then
    !$OMP do
    do j=jsvf-1,jevf+1 ; do i=isvf-1,ievf+1
      eta_sum(i,j) = 0.0 ; eta_wtd(i,j) = 0.0
    enddo ; enddo
  else
    !$OMP do
    do j=jsvf-1,jevf+1 ; do i=isvf-1,ievf+1
      eta_wtd(i,j) = 0.0
    enddo ; enddo
  endif
  !$OMP do
  do j=jsvf-1,jevf+1 ; do I=isvf-1,ievf
    ubt_sum(I,j) = 0.0 ; uhbt_sum(I,j) = 0.0
    PFu_bt_sum(I,j) = 0.0 ; Coru_bt_sum(I,j) = 0.0
    ubt_wtd(I,j) = 0.0 ; ubt_trans(I,j) = 0.0
  enddo ; enddo
  !$OMP do
  do J=jsvf-1,jevf ; do i=isvf-1,ievf+1
    vbt_sum(i,J) = 0.0 ; vhbt_sum(i,J) = 0.0
    PFv_bt_sum(i,J) = 0.0 ; Corv_bt_sum(i,J) = 0.0
    vbt_wtd(i,J) = 0.0 ; vbt_trans(i,J) = 0.0
  enddo ; enddo

  ! Set the mass source, after first initializing the halos to 0.
  !$OMP do
  do j=jsvf-1,jevf+1 ; do i=isvf-1,ievf+1 ; eta_src(i,j) = 0.0 ; enddo ; enddo
  if (CS%bound_BT_corr) then ; if ((use_BT_Cont.or.integral_BT_cont) .and. CS%BT_cont_bounds) then
    do j=js,je ; do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
      if (CS%eta_cor(i,j) > 0.0) then
        !   Limit the source (outward) correction to be a fraction the mass that
        ! can be transported out of the cell by velocities with a CFL number of CFL_cor.
        if (integral_BT_cont) then
          uint_cor = G%dxT(i,j) * CS%maxCFL_BT_cont
          vint_cor = G%dyT(i,j) * CS%maxCFL_BT_cont
          eta_cor_max = (CS%IareaT(i,j) * &
                   (((find_uhbt(uint_cor, BTCL_u(I,j)) + dt*uhbt0(I,j)) - &
                     (find_uhbt(-uint_cor, BTCL_u(I-1,j)) + dt*uhbt0(I-1,j))) + &
                    ((find_vhbt(vint_cor, BTCL_v(i,J)) + dt*vhbt0(i,J)) - &
                     (find_vhbt(-vint_cor, BTCL_v(i,J-1)) + dt*vhbt0(i,J-1))) ))
        else ! (use_BT_Cont) then
          u_max_cor = G%dxT(i,j) * (CS%maxCFL_BT_cont*Idt)
          v_max_cor = G%dyT(i,j) * (CS%maxCFL_BT_cont*Idt)
          eta_cor_max = dt * (CS%IareaT(i,j) * &
                   (((find_uhbt(u_max_cor, BTCL_u(I,j)) + uhbt0(I,j)) - &
                     (find_uhbt(-u_max_cor, BTCL_u(I-1,j)) + uhbt0(I-1,j))) + &
                    ((find_vhbt(v_max_cor, BTCL_v(i,J)) + vhbt0(i,J)) - &
                     (find_vhbt(-v_max_cor, BTCL_v(i,J-1)) + vhbt0(i,J-1))) ))
        endif
        CS%eta_cor(i,j) = min(CS%eta_cor(i,j), max(0.0, eta_cor_max))
      else
        ! Limit the sink (inward) correction to the amount of mass that is already inside the cell.
        Htot = eta(i,j)
        if (GV%Boussinesq) Htot = CS%bathyT(i,j)*GV%Z_to_H + eta(i,j)

        CS%eta_cor(i,j) = max(CS%eta_cor(i,j), -max(0.0,Htot))
      endif
    endif ; enddo ; enddo
  else ; do j=js,je ; do i=is,ie
    if (abs(CS%eta_cor(i,j)) > dt*CS%eta_cor_bound(i,j)) &
      CS%eta_cor(i,j) = sign(dt*CS%eta_cor_bound(i,j), CS%eta_cor(i,j))
  enddo ; enddo ; endif ; endif
  !$OMP do
  do j=js,je ; do i=is,ie
    eta_src(i,j) = G%mask2dT(i,j) * (Instep * CS%eta_cor(i,j))
  enddo ; enddo
!$OMP end parallel

  if (CS%dynamic_psurf) then
    ice_is_rigid = (associated(forces%rigidity_ice_u) .and. &
                    associated(forces%rigidity_ice_v))
    H_min_dyn = GV%Z_to_H * CS%Dmin_dyn_psurf
    if (ice_is_rigid .and. use_BT_cont) &
      call BT_cont_to_face_areas(BT_cont, Datu, Datv, G, US, MS, 0, .true.)
    if (ice_is_rigid) then
      !$OMP parallel do default(shared) private(Idt_max2,H_eff_dx2,dyn_coef_max,ice_strength)
      do j=js,je ; do i=is,ie
      ! First determine the maximum stable value for dyn_coef_eta.

      !   This estimate of the maximum stable time step is pretty accurate for
      ! gravity waves, but it is a conservative estimate since it ignores the
      ! stabilizing effect of the bottom drag.
      Idt_max2 = 0.5 * (dgeo_de * (1.0 + 2.0*bebt)) * (G%IareaT(i,j) * &
            ((gtot_E(i,j) * (Datu(I,j)*G%IdxCu(I,j)) + &
              gtot_W(i,j) * (Datu(I-1,j)*G%IdxCu(I-1,j))) + &
             (gtot_N(i,j) * (Datv(i,J)*G%IdyCv(i,J)) + &
              gtot_S(i,j) * (Datv(i,J-1)*G%IdyCv(i,J-1)))) + &
            ((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
             (G%CoriolisBu(I-1,J)**2 + G%CoriolisBu(I,J-1)**2)) * CS%BT_Coriolis_scale**2 )
      H_eff_dx2 = max(H_min_dyn * ((G%IdxT(i,j))**2 + (G%IdyT(i,j))**2), &
                      G%IareaT(i,j) * &
                        ((Datu(I,j)*G%IdxCu(I,j) + Datu(I-1,j)*G%IdxCu(I-1,j)) + &
                         (Datv(i,J)*G%IdyCv(i,J) + Datv(i,J-1)*G%IdyCv(i,J-1)) ) )
      dyn_coef_max = CS%const_dyn_psurf * max(0.0, 1.0 - dtbt**2 * Idt_max2) / &
                     (dtbt**2 * H_eff_dx2)

      ! ice_strength has units of [L2 Z-1 T-2 ~> m s-2]. rigidity_ice_[uv] has units of [L4 Z-1 T-1 ~> m3 s-1].
      ice_strength = ((forces%rigidity_ice_u(I,j) + forces%rigidity_ice_u(I-1,j)) + &
                      (forces%rigidity_ice_v(i,J) + forces%rigidity_ice_v(i,J-1))) / &
                      (CS%ice_strength_length**2 * dtbt)

      ! Units of dyn_coef: [L2 T-2 H-1 ~> m s-2 or m4 s-2 kg-1]
      dyn_coef_eta(i,j) = min(dyn_coef_max, ice_strength * GV%H_to_Z)
    enddo ; enddo ; endif
  endif

  if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
  if (id_clock_pass_pre > 0) call cpu_clock_begin(id_clock_pass_pre)
  if (nonblock_setup) then
    call start_group_pass(CS%pass_eta_bt_rem, CS%BT_Domain)
    ! The following halo update is not needed without wide halos.  RWH
  else
    call do_group_pass(CS%pass_eta_bt_rem, CS%BT_Domain)
    if (.not.use_BT_cont) &
      call do_group_pass(CS%pass_Dat_uv, CS%BT_Domain)
    call do_group_pass(CS%pass_force_hbt0_Cor_ref, CS%BT_Domain)
  endif
  if (id_clock_pass_pre > 0) call cpu_clock_end(id_clock_pass_pre)
  if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)

  ! Complete all of the outstanding halo updates.
  if (nonblock_setup) then
    if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
    if (id_clock_pass_pre > 0) call cpu_clock_begin(id_clock_pass_pre)

    if (.not.use_BT_cont) call complete_group_pass(CS%pass_Dat_uv, CS%BT_Domain)
    call complete_group_pass(CS%pass_force_hbt0_Cor_ref, CS%BT_Domain)
    call complete_group_pass(CS%pass_eta_bt_rem, CS%BT_Domain)

    if (id_clock_pass_pre > 0) call cpu_clock_end(id_clock_pass_pre)
    if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)
  endif

  if (CS%debug) then
    call uvchksum("BT [uv]hbt", uhbt, vhbt, CS%debug_BT_HI, haloshift=0, &
                  scale=US%s_to_T*US%L_to_m**2*GV%H_to_m)
    call uvchksum("BT Initial [uv]bt", ubt, vbt, CS%debug_BT_HI, haloshift=0, scale=US%L_T_to_m_s)
    call hchksum(eta, "BT Initial eta", CS%debug_BT_HI, haloshift=0, scale=GV%H_to_m)
    call uvchksum("BT BT_force_[uv]", BT_force_u, BT_force_v, &
                  CS%debug_BT_HI, haloshift=0, scale=US%L_T2_to_m_s2)
    if (interp_eta_PF) then
      call hchksum(eta_PF_1, "BT eta_PF_1",CS%debug_BT_HI,haloshift=0, scale=GV%H_to_m)
      call hchksum(d_eta_PF, "BT d_eta_PF",CS%debug_BT_HI,haloshift=0, scale=GV%H_to_m)
    else
      call hchksum(eta_PF, "BT eta_PF",CS%debug_BT_HI,haloshift=0, scale=GV%H_to_m)
      call hchksum(eta_PF_in, "BT eta_PF_in",G%HI,haloshift=0, scale=GV%H_to_m)
    endif
    call uvchksum("BT Cor_ref_[uv]", Cor_ref_u, Cor_ref_v, CS%debug_BT_HI, haloshift=0, scale=US%L_T2_to_m_s2)
    call uvchksum("BT [uv]hbt0", uhbt0, vhbt0, CS%debug_BT_HI, haloshift=0, &
                  scale=US%L_to_m**2*US%s_to_T*GV%H_to_m)
    if (.not. use_BT_cont) then
      call uvchksum("BT Dat[uv]", Datu, Datv, CS%debug_BT_HI, haloshift=1, scale=US%L_to_m*GV%H_to_m)
    endif
    call uvchksum("BT wt_[uv]", wt_u, wt_v, G%HI, haloshift=0, &
                  symmetric=.true., omit_corners=.true., scalar_pair=.true.)
    call uvchksum("BT frhat[uv]", CS%frhatu, CS%frhatv, G%HI, haloshift=0, &
                  symmetric=.true., omit_corners=.true., scalar_pair=.true.)
    call uvchksum("BT bc_accel_[uv]", bc_accel_u, bc_accel_v, G%HI, haloshift=0, scale=US%L_T2_to_m_s2)
    call uvchksum("BT IDat[uv]", CS%IDatu, CS%IDatv, G%HI, haloshift=0, &
                  scale=US%m_to_Z, scalar_pair=.true.)
    call uvchksum("BT visc_rem_[uv]", visc_rem_u, visc_rem_v, G%HI, &
                  haloshift=1, scalar_pair=.true.)
  endif

  if (CS%id_ubtdt > 0) then
    do j=js-1,je+1 ; do I=is-1,ie
      ubt_st(I,j) = ubt(I,j)
    enddo ; enddo
  endif
  if (CS%id_vbtdt > 0) then
    do J=js-1,je ; do i=is-1,ie+1
      vbt_st(i,J) = vbt(i,J)
    enddo ; enddo
  endif

  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_eta_st > 0) call post_data(CS%id_eta_st, eta(isd:ied,jsd:jed), CS%diag)
    if (CS%id_ubt_st > 0) call post_data(CS%id_ubt_st, ubt(IsdB:IedB,jsd:jed), CS%diag)
    if (CS%id_vbt_st > 0) call post_data(CS%id_vbt_st, vbt(isd:ied,JsdB:JedB), CS%diag)
  endif

  if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
  if (id_clock_calc > 0) call cpu_clock_begin(id_clock_calc)

  if (project_velocity) then ; eta_PF_BT => eta ; else ; eta_PF_BT => eta_pred ; endif

  if (CS%dt_bt_filter >= 0.0) then
    dt_filt = 0.5 * max(0.0, min(CS%dt_bt_filter, 2.0*dt))
  else
    dt_filt = 0.5 * max(0.0, dt * min(-CS%dt_bt_filter, 2.0))
  endif
  nfilter = ceiling(dt_filt / dtbt)

  if (nstep+nfilter==0 ) call MOM_error(FATAL, &
      "btstep: number of barotropic step (nstep+nfilter) is 0")

  ! Set up the normalized weights for the filtered velocity.
  sum_wt_vel = 0.0 ; sum_wt_eta = 0.0 ; sum_wt_accel = 0.0 ; sum_wt_trans = 0.0
  allocate(wt_vel(nstep+nfilter)) ; allocate(wt_eta(nstep+nfilter))
  allocate(wt_trans(nstep+nfilter+1)) ; allocate(wt_accel(nstep+nfilter+1))
  allocate(wt_accel2(nstep+nfilter+1))
  do n=1,nstep+nfilter
    ! Modify this to use a different filter...

    ! This is a filter that ramps down linearly over a time dt_filt.
    if ( (n==nstep) .or. (dt_filt - abs(n-nstep)*dtbt >= 0.0)) then
      wt_vel(n) = 1.0  ; wt_eta(n) = 1.0
    elseif (dtbt + dt_filt - abs(n-nstep)*dtbt > 0.0) then
      wt_vel(n) = 1.0 + (dt_filt / dtbt) - abs(n-nstep) ; wt_eta(n) = wt_vel(n)
    else
      wt_vel(n) = 0.0  ; wt_eta(n) = 0.0
    endif
    ! This is a simple stepfunction filter.
    ! if (n < nstep-nfilter) then ; wt_vel(n) = 0.0 ; else ; wt_vel(n) = 1.0 ; endif
    ! wt_eta(n) = wt_vel(n)

    ! The rest should not be changed.
    sum_wt_vel = sum_wt_vel + wt_vel(n) ; sum_wt_eta = sum_wt_eta + wt_eta(n)
  enddo
  wt_trans(nstep+nfilter+1) = 0.0 ; wt_accel(nstep+nfilter+1) = 0.0
  do n=nstep+nfilter,1,-1
    wt_trans(n) = wt_trans(n+1) + wt_eta(n)
    wt_accel(n) = wt_accel(n+1) + wt_vel(n)
    sum_wt_accel = sum_wt_accel + wt_accel(n) ; sum_wt_trans = sum_wt_trans + wt_trans(n)
  enddo
  ! Normalize the weights.
  I_sum_wt_vel = 1.0 / sum_wt_vel ; I_sum_wt_accel = 1.0 / sum_wt_accel
  I_sum_wt_eta = 1.0 / sum_wt_eta ; I_sum_wt_trans = 1.0 / sum_wt_trans
  do n=1,nstep+nfilter
    wt_vel(n) = wt_vel(n) * I_sum_wt_vel
    if (CS%answers_2018) then
      wt_accel2(n) = wt_accel(n)
     ! wt_trans(n) = wt_trans(n) * I_sum_wt_trans
    else
      wt_accel2(n) = wt_accel(n) * I_sum_wt_accel
      wt_trans(n) = wt_trans(n) * I_sum_wt_trans
    endif
    wt_accel(n) = wt_accel(n) * I_sum_wt_accel
    wt_eta(n) = wt_eta(n) * I_sum_wt_eta
  enddo

  sum_wt_vel = 0.0 ; sum_wt_eta = 0.0 ; sum_wt_accel = 0.0 ; sum_wt_trans = 0.0

  ! The following loop contains all of the time steps.
  isv=is ; iev=ie ; jsv=js ; jev=je
  do n=1,nstep+nfilter

    sum_wt_vel = sum_wt_vel + wt_vel(n)
    sum_wt_eta = sum_wt_eta + wt_eta(n)
    sum_wt_accel = sum_wt_accel + wt_accel2(n)
    sum_wt_trans = sum_wt_trans + wt_trans(n)

    if (CS%clip_velocity) then
      do j=jsv,jev ; do I=isv-1,iev
        if ((ubt(I,j) * (dt * G%dy_Cu(I,j))) * G%IareaT(i+1,j) < -CS%CFL_trunc) then
          ! Add some error reporting later.
          ubt(I,j) = (-0.95*CS%CFL_trunc) * (G%areaT(i+1,j) / (dt * G%dy_Cu(I,j)))
        elseif ((ubt(I,j) * (dt * G%dy_Cu(I,j))) * G%IareaT(i,j) > CS%CFL_trunc) then
          ! Add some error reporting later.
          ubt(I,j) = (0.95*CS%CFL_trunc) * (G%areaT(i,j) / (dt * G%dy_Cu(I,j)))
        endif
      enddo ; enddo
      do J=jsv-1,jev ; do i=isv,iev
        if ((vbt(i,J) * (dt * G%dx_Cv(i,J))) * G%IareaT(i,j+1) < -CS%CFL_trunc) then
          ! Add some error reporting later.
          vbt(i,J) = (-0.9*CS%CFL_trunc) * (G%areaT(i,j+1) / (dt * G%dx_Cv(i,J)))
        elseif ((vbt(i,J) * (dt * G%dx_Cv(i,J))) * G%IareaT(i,j) > CS%CFL_trunc) then
          ! Add some error reporting later.
          vbt(i,J) = (0.9*CS%CFL_trunc) * (G%areaT(i,j) / (dt * G%dx_Cv(i,J)))
        endif
      enddo ; enddo
    endif

    if ((iev - stencil < ie) .or. (jev - stencil < je)) then
      if (id_clock_calc > 0) call cpu_clock_end(id_clock_calc)
      call do_group_pass(CS%pass_eta_ubt, CS%BT_Domain, clock=id_clock_pass_step)
      isv = isvf ; iev = ievf ; jsv = jsvf ; jev = jevf
      if (id_clock_calc > 0) call cpu_clock_begin(id_clock_calc)
    else
      isv = isv+stencil ; iev = iev-stencil
      jsv = jsv+stencil ; jev = jev-stencil
    endif

    if ((.not.use_BT_cont) .and. CS%Nonlinear_continuity .and. &
        (CS%Nonlin_cont_update_period > 0)) then
      if ((n>1) .and. (mod(n-1,CS%Nonlin_cont_update_period) == 0)) &
        call find_face_areas(Datu, Datv, G, GV, US, CS, MS, eta, 1+iev-ie)
    endif

    if (integral_BT_cont) then
      !$OMP parallel do default(shared)
      do j=jsv-1,jev+1 ; do I=isv-2,iev+1
        ubt_int_prev(I,j) = ubt_int(I,j) ; uhbt_int_prev(I,j) = uhbt_int(I,j)
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=jsv-2,jev+1 ; do i=isv-1,iev+1
        vbt_int_prev(i,J) = vbt_int(i,J) ; vhbt_int_prev(i,J) = vhbt_int(i,J)
      enddo ; enddo
    endif

    !$OMP parallel default(shared) private(vel_prev, ioff, joff)
    if (CS%dynamic_psurf .or. .not.project_velocity) then
      if (integral_BT_cont) then
        !$OMP do
        do j=jsv-1,jev+1 ; do I=isv-2,iev+1
          uhbt_int(I,j) = find_uhbt(ubt_int(I,j) + dtbt*ubt(I,j), BTCL_u(I,j)) + n*dtbt*uhbt0(I,j)
        enddo ; enddo
        !$OMP end do nowait
        !$OMP do
        do J=jsv-2,jev+1 ; do i=isv-1,iev+1
          vhbt_int(i,J) = find_vhbt(vbt_int(i,J) + dtbt*vbt(i,J), BTCL_v(i,J)) + n*dtbt*vhbt0(i,J)
        enddo ; enddo
        !$OMP do
        do j=jsv-1,jev+1 ; do i=isv-1,iev+1
          eta_pred(i,j) = (eta_IC(i,j) + n*eta_src(i,j)) + CS%IareaT(i,j) * &
                     ((uhbt_int(I-1,j) - uhbt_int(I,j)) + (vhbt_int(i,J-1) - vhbt_int(i,J)))
        enddo ; enddo
      elseif (use_BT_cont) then
        !$OMP do
        do j=jsv-1,jev+1 ; do I=isv-2,iev+1
          uhbt(I,j) = find_uhbt(ubt(I,j), BTCL_u(I,j)) + uhbt0(I,j)
        enddo ; enddo
        !$OMP do
        do J=jsv-2,jev+1 ; do i=isv-1,iev+1
          vhbt(i,J) = find_vhbt(vbt(i,J), BTCL_v(i,J)) + vhbt0(i,J)
        enddo ; enddo
        !$OMP do
        do j=jsv-1,jev+1 ; do i=isv-1,iev+1
          eta_pred(i,j) = (eta(i,j) + eta_src(i,j)) + (dtbt * CS%IareaT(i,j)) * &
                     ((uhbt(I-1,j) - uhbt(I,j)) + (vhbt(i,J-1) - vhbt(i,J)))
        enddo ; enddo
      else
        !$OMP do
        do j=jsv-1,jev+1 ; do i=isv-1,iev+1
          eta_pred(i,j) = (eta(i,j) + eta_src(i,j)) + (dtbt * CS%IareaT(i,j)) * &
              (((Datu(I-1,j)*ubt(I-1,j) + uhbt0(I-1,j)) - &
                (Datu(I,j)*ubt(I,j) + uhbt0(I,j))) + &
               ((Datv(i,J-1)*vbt(i,J-1) + vhbt0(i,J-1)) - &
                (Datv(i,J)*vbt(i,J) + vhbt0(i,J))))
        enddo ; enddo
      endif

      if (CS%dynamic_psurf) then
        !$OMP do
        do j=jsv-1,jev+1 ; do i=isv-1,iev+1
          p_surf_dyn(i,j) = dyn_coef_eta(i,j) * (eta_pred(i,j) - eta(i,j))
        enddo ; enddo
      endif
    endif

    ! Recall that just outside the do n loop, there is code like...
    !  eta_PF_BT => eta_pred ; if (project_velocity) eta_PF_BT => eta

    if (find_etaav) then
      !$OMP do
      do j=js,je ; do i=is,ie
        eta_sum(i,j) = eta_sum(i,j) + wt_accel2(n) * eta_PF_BT(i,j)
      enddo ; enddo
      !$OMP end do nowait
    endif

    if (interp_eta_PF) then
      wt_end = n*Instep  ! This could be (n-0.5)*Instep.
      !$OMP do
      do j=jsv-1,jev+1 ; do i=isv-1,iev+1
        eta_PF(i,j) = eta_PF_1(i,j) + wt_end*d_eta_PF(i,j)
      enddo ; enddo
    endif

    if (apply_OBC_flather .or. apply_OBC_open) then
      !$OMP do
      do j=jsv,jev ; do I=isv-2,iev+1
        ubt_old(I,j) = ubt(I,j)
      enddo ; enddo
      !$OMP do
      do J=jsv-2,jev+1 ; do i=isv,iev
        vbt_old(i,J) = vbt(i,J)
      enddo ; enddo
    endif

    if (apply_OBCs) then
      if (MOD(n+G%first_direction,2)==1) then
        ioff = 1; joff = 0
      else
        ioff = 0; joff = 1
      endif

      if (CS%BT_OBC%apply_u_OBCs) then  ! save the old value of ubt and uhbt
        !$OMP do
        do j=jsv-joff,jev+joff ; do I=isv-1,iev
          ubt_prev(I,j) = ubt(I,j) ; uhbt_prev(I,j) = uhbt(I,j)
          ubt_sum_prev(I,j) = ubt_sum(I,j) ; uhbt_sum_prev(I,j) = uhbt_sum(I,j) ; ubt_wtd_prev(I,j) = ubt_wtd(I,j)
        enddo ; enddo
      endif

      if (CS%BT_OBC%apply_v_OBCs) then  ! save the old value of vbt and vhbt
        !$OMP do
        do J=jsv-1,jev ; do i=isv-ioff,iev+ioff
          vbt_prev(i,J) = vbt(i,J) ; vhbt_prev(i,J) = vhbt(i,J)
          vbt_sum_prev(i,J) = vbt_sum(i,J) ; vhbt_sum_prev(i,J) = vhbt_sum(i,J) ; vbt_wtd_prev(i,J) = vbt_wtd(i,J)
        enddo ; enddo
      endif
    endif

    if (MOD(n+G%first_direction,2)==1) then
      ! On odd-steps, update v first.
      !$OMP do schedule(static)
      do J=jsv-1,jev ; do i=isv-1,iev+1
        Cor_v(i,J) = -1.0*((amer(I-1,j) * ubt(I-1,j) + cmer(I,j+1) * ubt(I,j+1)) + &
               (bmer(I,j) * ubt(I,j) + dmer(I-1,j+1) * ubt(I-1,j+1))) - Cor_ref_v(i,J)
        PFv(i,J) = ((eta_PF_BT(i,j)-eta_PF(i,j))*gtot_N(i,j) - &
                     (eta_PF_BT(i,j+1)-eta_PF(i,j+1))*gtot_S(i,j+1)) * &
                   dgeo_de * CS%IdyCv(i,J)
      enddo ; enddo
      !$OMP end do nowait
      if (CS%dynamic_psurf) then
        !$OMP do schedule(static)
        do J=jsv-1,jev ; do i=isv-1,iev+1
          PFv(i,J) = PFv(i,J) + (p_surf_dyn(i,j) - p_surf_dyn(i,j+1)) * CS%IdyCv(i,J)
        enddo ; enddo
        !$OMP end do nowait
      endif

      if (CS%BT_OBC%apply_v_OBCs) then  ! zero out PF across boundary
        !$OMP do schedule(static)
        do J=jsv-1,jev ; do i=isv-1,iev+1 ; if (OBC%segnum_v(i,J) /= OBC_NONE) then
          PFv(i,J) = 0.0
        endif ; enddo ; enddo
        !$OMP end do nowait
      endif

      !$OMP do schedule(static)
      do J=jsv-1,jev ; do i=isv-1,iev+1
        vel_prev = vbt(i,J)
        vbt(i,J) = bt_rem_v(i,J) * (vbt(i,J) + &
             dtbt * ((BT_force_v(i,J) + Cor_v(i,J)) + PFv(i,J)))
        if (abs(vbt(i,J)) < CS%vel_underflow) vbt(i,J) = 0.0
        vbt_trans(i,J) = trans_wt1*vbt(i,J) + trans_wt2*vel_prev
      enddo ; enddo

      if (CS%linear_wave_drag) then
        !$OMP do schedule(static)
        do J=jsv-1,jev ; do i=isv-1,iev+1
          v_accel_bt(i,J) = v_accel_bt(i,J) + wt_accel(n) * &
              ((Cor_v(i,J) + PFv(i,J)) - vbt(i,J)*Rayleigh_v(i,J))
        enddo ; enddo
      else
        !$OMP do schedule(static)
        do J=jsv-1,jev ; do i=isv-1,iev+1
          v_accel_bt(i,J) = v_accel_bt(i,J) + wt_accel(n) * (Cor_v(i,J) + PFv(i,J))
        enddo ; enddo
      endif

      if (integral_BT_cont) then
        !$OMP do schedule(static)
        do J=jsv-1,jev ; do i=isv-1,iev+1
          vbt_int(i,J) = vbt_int(i,J) + dtbt * vbt_trans(i,J)
          vhbt_int(i,J) = find_vhbt(vbt_int(i,J), BTCL_v(i,J)) + n*dtbt*vhbt0(i,J)
          ! Estimate the mass flux within a single timestep to take the filtered average.
          vhbt(i,J) = (vhbt_int(i,J) - vhbt_int_prev(i,J)) * Idtbt
        enddo ; enddo
      elseif (use_BT_cont) then
        !$OMP do schedule(static)
        do J=jsv-1,jev ; do i=isv-1,iev+1
          vhbt(i,J) = find_vhbt(vbt_trans(i,J), BTCL_v(i,J)) + vhbt0(i,J)
        enddo ; enddo
        !$OMP end do nowait
      else
        !$OMP do schedule(static)
        do J=jsv-1,jev ; do i=isv-1,iev+1
          vhbt(i,J) = Datv(i,J)*vbt_trans(i,J) + vhbt0(i,J)
        enddo ; enddo
        !$OMP end do nowait
      endif
      if (CS%BT_OBC%apply_v_OBCs) then  ! copy back the value for v-points on the boundary.
        !$OMP do schedule(static)
        do J=jsv-1,jev ; do i=isv-1,iev+1 ; if (OBC%segnum_v(i,J) /= OBC_NONE) then
          vbt(i,J) = vbt_prev(i,J) ; vhbt(i,J) = vhbt_prev(i,J)
        endif ; enddo ; enddo
      endif
      ! Now update the zonal velocity.
      !$OMP do schedule(static)
      do j=jsv,jev ; do I=isv-1,iev
        Cor_u(I,j) = ((azon(I,j) * vbt(i+1,J) + czon(I,j) * vbt(i,J-1)) + &
                      (bzon(I,j) * vbt(i,J) + dzon(I,j) * vbt(i+1,J-1))) - &
                     Cor_ref_u(I,j)
        PFu(I,j) = ((eta_PF_BT(i,j)-eta_PF(i,j))*gtot_E(i,j) - &
                     (eta_PF_BT(i+1,j)-eta_PF(i+1,j))*gtot_W(i+1,j)) * &
                    dgeo_de * CS%IdxCu(I,j)
      enddo ; enddo
      !$OMP end do nowait

      if (CS%dynamic_psurf) then
        !$OMP do schedule(static)
        do j=jsv,jev ; do I=isv-1,iev
          PFu(I,j) = PFu(I,j) + (p_surf_dyn(i,j) - p_surf_dyn(i+1,j)) * CS%IdxCu(I,j)
        enddo ; enddo
        !$OMP end do nowait
      endif

      if (CS%BT_OBC%apply_u_OBCs) then  ! zero out pressure force across boundary
        !$OMP do schedule(static)
        do j=jsv,jev ; do I=isv-1,iev ; if (OBC%segnum_u(I,j) /= OBC_NONE) then
          PFu(I,j) = 0.0
        endif ; enddo ; enddo
        !$OMP end do nowait
      endif

      !$OMP do schedule(static)
      do j=jsv,jev ; do I=isv-1,iev
        vel_prev = ubt(I,j)
        ubt(I,j) = bt_rem_u(I,j) * (ubt(I,j) + &
             dtbt * ((BT_force_u(I,j) + Cor_u(I,j)) + PFu(I,j)))
        if (abs(ubt(I,j)) < CS%vel_underflow) ubt(I,j) = 0.0
        ubt_trans(I,j) = trans_wt1*ubt(I,j) + trans_wt2*vel_prev
      enddo ; enddo
      !$OMP end do nowait

      if (CS%linear_wave_drag) then
        !$OMP do schedule(static)
        do j=jsv,jev ; do I=isv-1,iev
          u_accel_bt(I,j) = u_accel_bt(I,j) + wt_accel(n) * &
             ((Cor_u(I,j) + PFu(I,j)) - ubt(I,j)*Rayleigh_u(I,j))
        enddo ; enddo
        !$OMP end do nowait
      else
        !$OMP do schedule(static)
        do j=jsv,jev ; do I=isv-1,iev
          u_accel_bt(I,j) = u_accel_bt(I,j) + wt_accel(n) * (Cor_u(I,j) + PFu(I,j))
        enddo ; enddo
        !$OMP end do nowait
      endif

      if (integral_BT_cont) then
        !$OMP do schedule(static)
        do j=jsv,jev ; do I=isv-1,iev
          ubt_int(I,j) = ubt_int(I,j) + dtbt * ubt_trans(I,j)
          uhbt_int(I,j) = find_uhbt(ubt_int(I,j), BTCL_u(I,j)) + n*dtbt*uhbt0(I,j)
          ! Estimate the mass flux within a single timestep to take the filtered average.
          uhbt(I,j) = (uhbt_int(I,j) - uhbt_int_prev(I,j)) * Idtbt
        enddo ; enddo
      elseif (use_BT_cont) then
        !$OMP do schedule(static)
        do j=jsv,jev ; do I=isv-1,iev
          uhbt(I,j) = find_uhbt(ubt_trans(I,j), BTCL_u(I,j)) + uhbt0(I,j)
        enddo ; enddo
      else
        !$OMP do schedule(static)
        do j=jsv,jev ; do I=isv-1,iev
          uhbt(I,j) = Datu(I,j)*ubt_trans(I,j) + uhbt0(I,j)
        enddo ; enddo
      endif
      if (CS%BT_OBC%apply_u_OBCs) then  ! copy back the value for u-points on the boundary.
        !$OMP do schedule(static)
        do j=jsv,jev ; do I=isv-1,iev ; if (OBC%segnum_u(I,j) /= OBC_NONE) then
          ubt(I,j) = ubt_prev(I,j) ; uhbt(I,j) = uhbt_prev(I,j)
        endif ; enddo ; enddo
      endif
    else
      ! On even steps, update u first.
      !$OMP do schedule(static)
      do j=jsv-1,jev+1 ; do I=isv-1,iev
        Cor_u(I,j) = ((azon(I,j) * vbt(i+1,J) + czon(I,j) * vbt(i,J-1)) + &
                      (bzon(I,j) * vbt(i,J) +  dzon(I,j) * vbt(i+1,J-1))) - &
                     Cor_ref_u(I,j)
        PFu(I,j) = ((eta_PF_BT(i,j)-eta_PF(i,j))*gtot_E(i,j) - &
                     (eta_PF_BT(i+1,j)-eta_PF(i+1,j))*gtot_W(i+1,j)) * &
                     dgeo_de * CS%IdxCu(I,j)
      enddo ; enddo
      !$OMP end do nowait

      if (CS%dynamic_psurf) then
        !$OMP do schedule(static)
        do j=jsv-1,jev+1 ; do I=isv-1,iev
          PFu(I,j) = PFu(I,j) + (p_surf_dyn(i,j) - p_surf_dyn(i+1,j)) * CS%IdxCu(I,j)
        enddo ; enddo
        !$OMP end do nowait
      endif

      if (CS%BT_OBC%apply_u_OBCs) then  ! zero out pressure force across boundary
        !$OMP do schedule(static)
        do j=jsv,jev ; do I=isv-1,iev ; if (OBC%segnum_u(I,j) /= OBC_NONE) then
          PFu(I,j) = 0.0
        endif ; enddo ; enddo
      endif

      !$OMP do schedule(static)
      do j=jsv-1,jev+1 ; do I=isv-1,iev
        vel_prev = ubt(I,j)
        ubt(I,j) = bt_rem_u(I,j) * (ubt(I,j) + &
             dtbt * ((BT_force_u(I,j) + Cor_u(I,j)) + PFu(I,j)))
        if (abs(ubt(I,j)) < CS%vel_underflow) ubt(I,j) = 0.0
        ubt_trans(I,j) = trans_wt1*ubt(I,j) + trans_wt2*vel_prev
      enddo ; enddo

      if (CS%linear_wave_drag) then
        !$OMP do schedule(static)
        do j=jsv-1,jev+1 ; do I=isv-1,iev
          u_accel_bt(I,j) = u_accel_bt(I,j) + wt_accel(n) * &
              ((Cor_u(I,j) + PFu(I,j)) - ubt(I,j)*Rayleigh_u(I,j))
        enddo ; enddo
      else
        !$OMP do schedule(static)
        do j=jsv-1,jev+1 ; do I=isv-1,iev
          u_accel_bt(I,j) = u_accel_bt(I,j) + wt_accel(n) * (Cor_u(I,j) + PFu(I,j))
        enddo ; enddo
      endif

      if (integral_BT_cont) then
        !$OMP do schedule(static)
        do j=jsv-1,jev+1 ; do I=isv-1,iev
          ubt_int(I,j) = ubt_int(I,j) + dtbt * ubt_trans(I,j)
          uhbt_int(I,j) = find_uhbt(ubt_int(I,j), BTCL_u(I,j)) + n*dtbt*uhbt0(I,j)
          ! Estimate the mass flux within a single timestep to take the filtered average.
          uhbt(I,j) = (uhbt_int(I,j) - uhbt_int_prev(I,j)) * Idtbt
        enddo ; enddo
      elseif (use_BT_cont) then
        !$OMP do schedule(static)
        do j=jsv-1,jev+1 ; do I=isv-1,iev
          uhbt(I,j) = find_uhbt(ubt_trans(I,j), BTCL_u(I,j)) + uhbt0(I,j)
        enddo ; enddo
        !$OMP end do nowait
      else
        !$OMP do schedule(static)
        do j=jsv-1,jev+1 ; do I=isv-1,iev
          uhbt(I,j) = Datu(I,j)*ubt_trans(I,j) + uhbt0(I,j)
        enddo ; enddo
        !$OMP end do nowait
      endif
      if (CS%BT_OBC%apply_u_OBCs) then  ! copy back the value for u-points on the boundary.
        !$OMP do schedule(static)
        do j=jsv-1,jev+1 ; do I=isv-1,iev ; if (OBC%segnum_u(I,j) /= OBC_NONE) then
          ubt(I,j) = ubt_prev(I,j) ; uhbt(I,j) = uhbt_prev(I,j)
        endif ; enddo ; enddo
      endif

      ! Now update the meridional velocity.
      if (CS%use_old_coriolis_bracket_bug) then
        !$OMP do schedule(static)
        do J=jsv-1,jev ; do i=isv,iev
          Cor_v(i,J) = -1.0*((amer(I-1,j) * ubt(I-1,j) + bmer(I,j) * ubt(I,j)) + &
                  (cmer(I,j+1) * ubt(I,j+1) + dmer(I-1,j+1) * ubt(I-1,j+1))) - Cor_ref_v(i,J)
          PFv(i,J) = ((eta_PF_BT(i,j)-eta_PF(i,j))*gtot_N(i,j) - &
                       (eta_PF_BT(i,j+1)-eta_PF(i,j+1))*gtot_S(i,j+1)) * &
                      dgeo_de * CS%IdyCv(i,J)
        enddo ; enddo
        !$OMP end do nowait
      else
        !$OMP do schedule(static)
        do J=jsv-1,jev ; do i=isv,iev
          Cor_v(i,J) = -1.0*((amer(I-1,j) * ubt(I-1,j) + cmer(I,j+1) * ubt(I,j+1)) + &
                  (bmer(I,j) * ubt(I,j) + dmer(I-1,j+1) * ubt(I-1,j+1))) - Cor_ref_v(i,J)
          PFv(i,J) = ((eta_PF_BT(i,j)-eta_PF(i,j))*gtot_N(i,j) - &
                       (eta_PF_BT(i,j+1)-eta_PF(i,j+1))*gtot_S(i,j+1)) * &
                      dgeo_de * CS%IdyCv(i,J)
        enddo ; enddo
        !$OMP end do nowait
      endif

      if (CS%dynamic_psurf) then
        !$OMP do schedule(static)
        do J=jsv-1,jev ; do i=isv,iev
          PFv(i,J) = PFv(i,J) + (p_surf_dyn(i,j) - p_surf_dyn(i,j+1)) * CS%IdyCv(i,J)
        enddo ; enddo
        !$OMP end do nowait
      endif

      if (CS%BT_OBC%apply_v_OBCs) then  ! zero out PF across boundary
        !$OMP do schedule(static)
        do J=jsv-1,jev ; do i=isv-1,iev+1 ; if (OBC%segnum_v(i,J) /= OBC_NONE) then
          PFv(i,J) = 0.0
        endif ; enddo ; enddo
      endif

      !$OMP do schedule(static)
      do J=jsv-1,jev ; do i=isv,iev
        vel_prev = vbt(i,J)
        vbt(i,J) = bt_rem_v(i,J) * (vbt(i,J) + &
             dtbt * ((BT_force_v(i,J) + Cor_v(i,J)) + PFv(i,J)))
        if (abs(vbt(i,J)) < CS%vel_underflow) vbt(i,J) = 0.0
        vbt_trans(i,J) = trans_wt1*vbt(i,J) + trans_wt2*vel_prev
      enddo ; enddo
      !$OMP end do nowait

      if (CS%linear_wave_drag) then
        !$OMP do schedule(static)
        do J=jsv-1,jev ; do i=isv,iev
          v_accel_bt(i,J) = v_accel_bt(i,J) + wt_accel(n) * &
             ((Cor_v(i,J) + PFv(i,J)) - vbt(i,J)*Rayleigh_v(i,J))
        enddo ; enddo
        !$OMP end do nowait
      else
        !$OMP do schedule(static)
        do J=jsv-1,jev ; do i=isv,iev
          v_accel_bt(i,J) = v_accel_bt(i,J) + wt_accel(n) * (Cor_v(i,J) + PFv(i,J))
        enddo ; enddo
        !$OMP end do nowait
      endif

      if (integral_BT_cont) then
        !$OMP do schedule(static)
        do J=jsv-1,jev ; do i=isv,iev
          vbt_int(i,J) = vbt_int(i,J) + dtbt * vbt_trans(i,J)
          vhbt_int(i,J) = find_vhbt(vbt_int(i,J), BTCL_v(i,J)) + n*dtbt*vhbt0(i,J)
          ! Estimate the mass flux within a single timestep to take the filtered average.
          vhbt(i,J) = (vhbt_int(i,J) - vhbt_int_prev(i,J)) * Idtbt
        enddo ; enddo
      elseif (use_BT_cont) then
        !$OMP do schedule(static)
        do J=jsv-1,jev ; do i=isv,iev
          vhbt(i,J) = find_vhbt(vbt_trans(i,J), BTCL_v(i,J)) + vhbt0(i,J)
        enddo ; enddo
      else
        !$OMP do schedule(static)
        do J=jsv-1,jev ; do i=isv,iev
          vhbt(i,J) = Datv(i,J)*vbt_trans(i,J) + vhbt0(i,J)
        enddo ; enddo
      endif
      if (CS%BT_OBC%apply_v_OBCs) then  ! copy back the value for v-points on the boundary.
        !$OMP do schedule(static)
        do J=jsv-1,jev ; do i=isv,iev ; if (OBC%segnum_v(i,J) /= OBC_NONE) then
          vbt(i,J) = vbt_prev(i,J); vhbt(i,J) = vhbt_prev(i,J)
        endif ; enddo ; enddo
      endif
    endif

    ! This might need to be moved outside of the OMP do loop directives.
    if (CS%debug_bt) then
      write(mesg,'("BT vel update ",I4)') n
      call uvchksum(trim(mesg)//" PF[uv]", PFu, PFv, CS%debug_BT_HI, haloshift=iev-ie, &
                    scale=US%L_T_to_m_s*US%s_to_T)
      call uvchksum(trim(mesg)//" Cor_[uv]", Cor_u, Cor_v, CS%debug_BT_HI, haloshift=iev-ie, &
                    scale=US%L_T_to_m_s*US%s_to_T)
      call uvchksum(trim(mesg)//" BT_force_[uv]", BT_force_u, BT_force_v, CS%debug_BT_HI, haloshift=iev-ie, &
                    scale=US%L_T_to_m_s*US%s_to_T)
      call uvchksum(trim(mesg)//" BT_rem_[uv]", BT_rem_u, BT_rem_v, CS%debug_BT_HI, &
                    haloshift=iev-ie, scalar_pair=.true.)
      call uvchksum(trim(mesg)//" [uv]bt", ubt, vbt, CS%debug_BT_HI, haloshift=iev-ie, &
                    scale=US%L_T_to_m_s)
      call uvchksum(trim(mesg)//" [uv]bt_trans", ubt_trans, vbt_trans, CS%debug_BT_HI, haloshift=iev-ie, &
                    scale=US%L_T_to_m_s)
      call uvchksum(trim(mesg)//" [uv]hbt", uhbt, vhbt, CS%debug_BT_HI, haloshift=iev-ie, &
                    scale=US%s_to_T*US%L_to_m**2*GV%H_to_m)
      if (integral_BT_cont) &
        call uvchksum(trim(mesg)//" [uv]hbt_int", uhbt_int, vhbt_int, CS%debug_BT_HI, haloshift=iev-ie, &
                      scale=US%L_to_m**2*GV%H_to_m)
    endif

    if (find_PF) then
      !$OMP do
      do j=js,je ; do I=is-1,ie
        PFu_bt_sum(I,j)  = PFu_bt_sum(I,j) + wt_accel2(n) * PFu(I,j)
      enddo ; enddo
      !$OMP end do nowait
      !$OMP do
      do J=js-1,je ; do i=is,ie
        PFv_bt_sum(i,J)  = PFv_bt_sum(i,J) + wt_accel2(n) * PFv(i,J)
      enddo ; enddo
      !$OMP end do nowait
    endif
    if (find_Cor) then
      !$OMP do
      do j=js,je ; do I=is-1,ie
        Coru_bt_sum(I,j) = Coru_bt_sum(I,j) + wt_accel2(n) * Cor_u(I,j)
      enddo ; enddo
      !$OMP end do nowait
      !$OMP do
      do J=js-1,je ; do i=is,ie
        Corv_bt_sum(i,J) = Corv_bt_sum(i,J) + wt_accel2(n) * Cor_v(i,J)
      enddo ; enddo
      !$OMP end do nowait
    endif

    !$OMP do
    do j=js,je ; do I=is-1,ie
      ubt_sum(I,j) = ubt_sum(I,j) + wt_trans(n) * ubt_trans(I,j)
      uhbt_sum(I,j) = uhbt_sum(I,j) + wt_trans(n) * uhbt(I,j)
      ubt_wtd(I,j) = ubt_wtd(I,j) + wt_vel(n) * ubt(I,j)
    enddo ; enddo
    !$OMP end do nowait
    !$OMP do
    do J=js-1,je ; do i=is,ie
      vbt_sum(i,J) = vbt_sum(i,J) + wt_trans(n) * vbt_trans(i,J)
      vhbt_sum(i,J) = vhbt_sum(i,J) + wt_trans(n) * vhbt(i,J)
      vbt_wtd(i,J) = vbt_wtd(i,J) + wt_vel(n) * vbt(i,J)
    enddo ; enddo
    !$OMP end do nowait

    if (apply_OBCs) then

      !$OMP single
      call apply_velocity_OBCs(OBC, ubt, vbt, uhbt, vhbt, &
             ubt_trans, vbt_trans, eta, ubt_old, vbt_old, CS%BT_OBC, &
             G, MS, US, iev-ie, dtbt, bebt, use_BT_cont, integral_BT_cont, &
             n*dtbt, Datu, Datv, BTCL_u, BTCL_v, uhbt0, vhbt0, &
             ubt_int_prev, vbt_int_prev, uhbt_int_prev, vhbt_int_prev)
      !$OMP end single

      if (CS%BT_OBC%apply_u_OBCs) then
        !$OMP do
        do j=js,je ; do I=is-1,ie
          if (OBC%segnum_u(I,j) /= OBC_NONE) then
            ! Update the summed and integrated quantities from the saved previous values.
            ubt_sum(I,j) = ubt_sum_prev(I,j) + wt_trans(n) * ubt_trans(I,j)
            uhbt_sum(I,j) = uhbt_sum_prev(I,j) + wt_trans(n) * uhbt(I,j)
            ubt_wtd(I,j) = ubt_wtd_prev(I,j) + wt_vel(n) * ubt(I,j)
            if (integral_BT_cont) then
              uhbt_int(I,j) = uhbt_int_prev(I,j) + dtbt * uhbt(I,j)
              ubt_int(I,j) = ubt_int_prev(I,j) + dtbt * ubt_trans(I,j)
            endif
          endif
        enddo ; enddo
      endif
      if (CS%BT_OBC%apply_v_OBCs) then
        !$OMP do
        do J=js-1,je ; do i=is,ie
          if (OBC%segnum_v(i,J) /= OBC_NONE) then
            ! Update the summed and integrated quantities from the saved previous values.
            vbt_sum(i,J) = vbt_sum_prev(i,J) + wt_trans(n) * vbt_trans(i,J)
            vhbt_sum(i,J) = vhbt_sum_prev(i,J) + wt_trans(n) * vhbt(i,J)
            vbt_wtd(i,J) = vbt_wtd_prev(i,J) + wt_vel(n) * vbt(i,J)
            if (integral_BT_cont) then
              vbt_int(i,J) = vbt_int_prev(i,J) + dtbt * vbt_trans(i,J)
              vhbt_int(i,J) = vhbt_int_prev(i,J) + dtbt * vhbt(i,J)
            endif
          endif
        enddo ; enddo
      endif
    endif

    if (CS%debug_bt) then
      call uvchksum("BT [uv]hbt just after OBC", uhbt, vhbt, CS%debug_BT_HI, haloshift=iev-ie, &
                    scale=US%s_to_T*US%L_to_m**2*GV%H_to_m)
      if (integral_BT_cont) &
        call uvchksum("BT [uv]hbt_int just after OBC", uhbt_int, vhbt_int, CS%debug_BT_HI, &
                      haloshift=iev-ie, scale=US%L_to_m**2*GV%H_to_m)
    endif

    if (integral_BT_cont) then
      !$OMP do
      do j=jsv,jev ; do i=isv,iev
        eta(i,j) = (eta_IC(i,j) + n*eta_src(i,j)) + CS%IareaT(i,j) * &
                   ((uhbt_int(I-1,j) - uhbt_int(I,j)) + (vhbt_int(i,J-1) - vhbt_int(i,J)))
        eta_wtd(i,j) = eta_wtd(i,j) + eta(i,j) * wt_eta(n)
      enddo ; enddo
    else
      !$OMP do
      do j=jsv,jev ; do i=isv,iev
        eta(i,j) = (eta(i,j) + eta_src(i,j)) + (dtbt * CS%IareaT(i,j)) * &
                   ((uhbt(I-1,j) - uhbt(I,j)) + (vhbt(i,J-1) - vhbt(i,J)))
        eta_wtd(i,j) = eta_wtd(i,j) + eta(i,j) * wt_eta(n)
      enddo ; enddo
    endif
    !$OMP end parallel

    if (do_hifreq_output) then
      time_step_end = time_bt_start + real_to_time(n*US%T_to_s*dtbt)
      call enable_averaging(US%T_to_s*dtbt, time_step_end, CS%diag)
      if (CS%id_ubt_hifreq > 0) call post_data(CS%id_ubt_hifreq, ubt(IsdB:IedB,jsd:jed), CS%diag)
      if (CS%id_vbt_hifreq > 0) call post_data(CS%id_vbt_hifreq, vbt(isd:ied,JsdB:JedB), CS%diag)
      if (CS%id_eta_hifreq > 0) call post_data(CS%id_eta_hifreq, eta(isd:ied,jsd:jed), CS%diag)
      if (CS%id_uhbt_hifreq > 0) call post_data(CS%id_uhbt_hifreq, uhbt(IsdB:IedB,jsd:jed), CS%diag)
      if (CS%id_vhbt_hifreq > 0) call post_data(CS%id_vhbt_hifreq, vhbt(isd:ied,JsdB:JedB), CS%diag)
      if (CS%id_eta_pred_hifreq > 0) call post_data(CS%id_eta_pred_hifreq, eta_PF_BT(isd:ied,jsd:jed), CS%diag)
    endif

    if (CS%debug_bt) then
      write(mesg,'("BT step ",I4)') n
      call uvchksum(trim(mesg)//" [uv]bt", ubt, vbt, CS%debug_BT_HI, haloshift=iev-ie, &
                    scale=US%L_T_to_m_s)
      call hchksum(eta, trim(mesg)//" eta", CS%debug_BT_HI, haloshift=iev-ie, scale=GV%H_to_m)
    endif

    if (GV%Boussinesq) then
      do j=js,je ; do i=is,ie
        if (eta(i,j) < -GV%Z_to_H*G%bathyT(i,j)) &
          call MOM_error(WARNING, "btstep: eta has dropped below bathyT.")
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        if (eta(i,j) < 0.0) &
          call MOM_error(WARNING, "btstep: negative eta in a non-Boussinesq barotropic solver.")
      enddo ; enddo
    endif

  enddo ! end of do n=1,ntimestep
  if (id_clock_calc > 0) call cpu_clock_end(id_clock_calc)
  if (id_clock_calc_post > 0) call cpu_clock_begin(id_clock_calc_post)

  ! Reset the time information in the diag type.
  if (do_hifreq_output) call enable_averaging(time_int_in, time_end_in, CS%diag)

  if (CS%answers_2018) then
    I_sum_wt_vel = 1.0 / sum_wt_vel ; I_sum_wt_eta = 1.0 / sum_wt_eta
    I_sum_wt_accel = 1.0 / sum_wt_accel ; I_sum_wt_trans = 1.0 / sum_wt_trans
  else
    I_sum_wt_vel = 1.0 ; I_sum_wt_eta = 1.0 ; I_sum_wt_accel = 1.0 ; I_sum_wt_trans = 1.0
  endif

  if (find_etaav) then ; do j=js,je ; do i=is,ie
    etaav(i,j) = eta_sum(i,j) * I_sum_wt_accel
  enddo ; enddo ; endif
  do j=js-1,je+1 ; do i=is-1,ie+1 ; e_anom(i,j) = 0.0 ; enddo ; enddo
  if (interp_eta_PF) then
    do j=js,je ; do i=is,ie
      e_anom(i,j) = dgeo_de * (0.5 * (eta(i,j) + eta_in(i,j)) - &
                               (eta_PF_1(i,j) + 0.5*d_eta_PF(i,j)))
    enddo ; enddo
  else
    do j=js,je ; do i=is,ie
      e_anom(i,j) = dgeo_de * (0.5 * (eta(i,j) + eta_in(i,j)) - eta_PF(i,j))
    enddo ; enddo
  endif
  if (apply_OBCs) then
    !!! Not safe for wide halos...
    if (CS%BT_OBC%apply_u_OBCs) then  ! copy back the value for u-points on the boundary.
      !GOMP parallel do default(shared)
      do j=js,je ; do I=is-1,ie
        l_seg = OBC%segnum_u(I,j)
        if (l_seg == OBC_NONE) cycle

        if (OBC%segment(l_seg)%direction == OBC_DIRECTION_E) then
          e_anom(i+1,j) = e_anom(i,j)
        elseif (OBC%segment(l_seg)%direction == OBC_DIRECTION_W) then
          e_anom(i,j) = e_anom(i+1,j)
        endif
      enddo ; enddo
    endif

    if (CS%BT_OBC%apply_v_OBCs) then  ! copy back the value for v-points on the boundary.
      !GOMP parallel do default(shared)
      do J=js-1,je ; do I=is,ie
        l_seg = OBC%segnum_v(i,J)
        if (l_seg == OBC_NONE) cycle

        if (OBC%segment(l_seg)%direction == OBC_DIRECTION_N) then
          e_anom(i,j+1) = e_anom(i,j)
        elseif (OBC%segment(l_seg)%direction == OBC_DIRECTION_S) then
          e_anom(i,j) = e_anom(i,j+1)
        endif
      enddo ; enddo
    endif
  endif

  ! It is possible that eta_out and eta_in are the same.
  do j=js,je ; do i=is,ie
    eta_out(i,j) = eta_wtd(i,j) * I_sum_wt_eta
  enddo ; enddo

  if (id_clock_calc_post > 0) call cpu_clock_end(id_clock_calc_post)
  if (id_clock_pass_post > 0) call cpu_clock_begin(id_clock_pass_post)
  if (G%nonblocking_updates) then
    call start_group_pass(CS%pass_e_anom, G%Domain)
  else
    if (find_etaav) call do_group_pass(CS%pass_etaav, G%Domain)
    call do_group_pass(CS%pass_e_anom, G%Domain)
  endif
  if (id_clock_pass_post > 0) call cpu_clock_end(id_clock_pass_post)
  if (id_clock_calc_post > 0) call cpu_clock_begin(id_clock_calc_post)

  if (CS%answers_2018) then
    do j=js,je ; do I=is-1,ie
      CS%ubtav(I,j) = ubt_sum(I,j) * I_sum_wt_trans
      uhbtav(I,j) = uhbt_sum(I,j) * I_sum_wt_trans
      ubt_wtd(I,j) = ubt_wtd(I,j) * I_sum_wt_vel
    enddo ; enddo

    do J=js-1,je ; do i=is,ie
      CS%vbtav(i,J) = vbt_sum(i,J) * I_sum_wt_trans
      vhbtav(i,J) = vhbt_sum(i,J) * I_sum_wt_trans
      vbt_wtd(i,J) = vbt_wtd(i,J) * I_sum_wt_vel
    enddo ; enddo
  else
    do j=js,je ; do I=is-1,ie
      CS%ubtav(I,j) = ubt_sum(I,j)
      uhbtav(I,j) = uhbt_sum(I,j)
    enddo ; enddo

    do J=js-1,je ; do i=is,ie
      CS%vbtav(i,J) = vbt_sum(i,J)
      vhbtav(i,J) = vhbt_sum(i,J)
    enddo ; enddo
  endif


  if (id_clock_calc_post > 0) call cpu_clock_end(id_clock_calc_post)
  if (id_clock_pass_post > 0) call cpu_clock_begin(id_clock_pass_post)
  if (G%nonblocking_updates) then
    call complete_group_pass(CS%pass_e_anom, G%Domain)
    if (find_etaav) call start_group_pass(CS%pass_etaav, G%Domain)
    call start_group_pass(CS%pass_ubta_uhbta, G%DoMain)
  else
    call do_group_pass(CS%pass_ubta_uhbta, G%Domain)
  endif
  if (id_clock_pass_post > 0) call cpu_clock_end(id_clock_pass_post)
  if (id_clock_calc_post > 0) call cpu_clock_begin(id_clock_calc_post)

  ! Now calculate each layer's accelerations.
  !$OMP parallel do default(shared)
  do k=1,nz
    do j=js,je ; do I=is-1,ie
      accel_layer_u(I,j,k) = (u_accel_bt(I,j) - &
           ((pbce(i+1,j,k) - gtot_W(i+1,j)) * e_anom(i+1,j) - &
            (pbce(i,j,k) - gtot_E(i,j)) * e_anom(i,j)) * CS%IdxCu(I,j) )
      if (abs(accel_layer_u(I,j,k)) < accel_underflow) accel_layer_u(I,j,k) = 0.0
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      accel_layer_v(i,J,k) = (v_accel_bt(i,J) - &
           ((pbce(i,j+1,k) - gtot_S(i,j+1)) * e_anom(i,j+1) - &
            (pbce(i,j,k) - gtot_N(i,j)) * e_anom(i,j)) * CS%IdyCv(i,J) )
      if (abs(accel_layer_v(i,J,k)) < accel_underflow) accel_layer_v(i,J,k) = 0.0
    enddo ; enddo
  enddo

  if (apply_OBCs) then
    ! Correct the accelerations at OBC velocity points, but only in the
    ! symmetric-memory computational domain, not in the wide halo regions.
    if (CS%BT_OBC%apply_u_OBCs) then ; do j=js,je ; do I=is-1,ie
      if (OBC%segnum_u(I,j) /= OBC_NONE) then
        u_accel_bt(I,j) = (ubt_wtd(I,j) - ubt_first(I,j)) / dt
        do k=1,nz ; accel_layer_u(I,j,k) = u_accel_bt(I,j) ; enddo
      endif
    enddo ; enddo ; endif
    if (CS%BT_OBC%apply_v_OBCs) then ; do J=js-1,je ; do i=is,ie
      if (OBC%segnum_v(i,J) /= OBC_NONE) then
        v_accel_bt(i,J) = (vbt_wtd(i,J) - vbt_first(i,J)) / dt
        do k=1,nz ; accel_layer_v(i,J,k) = v_accel_bt(i,J) ; enddo
      endif
    enddo ; enddo ; endif
  endif

  if (id_clock_calc_post > 0) call cpu_clock_end(id_clock_calc_post)

  ! Calculate diagnostic quantities.
  if (query_averaging_enabled(CS%diag)) then

    if (CS%gradual_BT_ICs) then
      do j=js,je ; do I=is-1,ie ; CS%ubt_IC(I,j) = ubt_wtd(I,j) ; enddo ; enddo
      do J=js-1,je ; do i=is,ie ; CS%vbt_IC(i,J) = vbt_wtd(i,J) ; enddo ; enddo
    endif

!  Offer various barotropic terms for averaging.
    if (CS%id_PFu_bt > 0) then
      do j=js,je ; do I=is-1,ie
        PFu_bt_sum(I,j) = PFu_bt_sum(I,j) * I_sum_wt_accel
      enddo ; enddo
      call post_data(CS%id_PFu_bt, PFu_bt_sum(IsdB:IedB,jsd:jed), CS%diag)
    endif
    if (CS%id_PFv_bt > 0) then
      do J=js-1,je ; do i=is,ie
        PFv_bt_sum(i,J) = PFv_bt_sum(i,J) * I_sum_wt_accel
      enddo ; enddo
      call post_data(CS%id_PFv_bt, PFv_bt_sum(isd:ied,JsdB:JedB), CS%diag)
    endif
    if (CS%id_Coru_bt > 0) then
      do j=js,je ; do I=is-1,ie
        Coru_bt_sum(I,j) = Coru_bt_sum(I,j) * I_sum_wt_accel
      enddo ; enddo
      call post_data(CS%id_Coru_bt, Coru_bt_sum(IsdB:IedB,jsd:jed), CS%diag)
    endif
    if (CS%id_Corv_bt > 0) then
      do J=js-1,je ; do i=is,ie
        Corv_bt_sum(i,J) = Corv_bt_sum(i,J) * I_sum_wt_accel
      enddo ; enddo
      call post_data(CS%id_Corv_bt, Corv_bt_sum(isd:ied,JsdB:JedB), CS%diag)
    endif
    if (CS%id_ubtdt > 0) then
      do j=js,je ; do I=is-1,ie
        ubt_dt(I,j) = (ubt_wtd(I,j) - ubt_st(I,j))*Idt
      enddo ; enddo
      call post_data(CS%id_ubtdt, ubt_dt(IsdB:IedB,jsd:jed), CS%diag)
    endif
    if (CS%id_vbtdt > 0) then
      do J=js-1,je ; do i=is,ie
        vbt_dt(i,J) = (vbt_wtd(i,J) - vbt_st(i,J))*Idt
      enddo ; enddo
      call post_data(CS%id_vbtdt, vbt_dt(isd:ied,JsdB:JedB), CS%diag)
    endif

    if (CS%id_ubtforce > 0) call post_data(CS%id_ubtforce, BT_force_u(IsdB:IedB,jsd:jed), CS%diag)
    if (CS%id_vbtforce > 0) call post_data(CS%id_vbtforce, BT_force_v(isd:ied,JsdB:JedB), CS%diag)
    if (CS%id_uaccel > 0) call post_data(CS%id_uaccel, u_accel_bt(IsdB:IedB,jsd:jed), CS%diag)
    if (CS%id_vaccel > 0) call post_data(CS%id_vaccel, v_accel_bt(isd:ied,JsdB:JedB), CS%diag)

    if (CS%id_eta_cor > 0) call post_data(CS%id_eta_cor, CS%eta_cor, CS%diag)
    if (CS%id_eta_bt > 0) call post_data(CS%id_eta_bt, eta_out, CS%diag)
    if (CS%id_gtotn > 0) call post_data(CS%id_gtotn, gtot_N(isd:ied,jsd:jed), CS%diag)
    if (CS%id_gtots > 0) call post_data(CS%id_gtots, gtot_S(isd:ied,jsd:jed), CS%diag)
    if (CS%id_gtote > 0) call post_data(CS%id_gtote, gtot_E(isd:ied,jsd:jed), CS%diag)
    if (CS%id_gtotw > 0) call post_data(CS%id_gtotw, gtot_W(isd:ied,jsd:jed), CS%diag)
    if (CS%id_ubt > 0) call post_data(CS%id_ubt, ubt_wtd(IsdB:IedB,jsd:jed), CS%diag)
    if (CS%id_vbt > 0) call post_data(CS%id_vbt, vbt_wtd(isd:ied,JsdB:JedB), CS%diag)
    if (CS%id_ubtav > 0) call post_data(CS%id_ubtav, CS%ubtav, CS%diag)
    if (CS%id_vbtav > 0) call post_data(CS%id_vbtav, CS%vbtav, CS%diag)
    if (CS%id_visc_rem_u > 0) call post_data(CS%id_visc_rem_u, visc_rem_u, CS%diag)
    if (CS%id_visc_rem_v > 0) call post_data(CS%id_visc_rem_v, visc_rem_v, CS%diag)

    if (CS%id_frhatu > 0) call post_data(CS%id_frhatu, CS%frhatu, CS%diag)
    if (CS%id_uhbt > 0) call post_data(CS%id_uhbt, uhbtav, CS%diag)
    if (CS%id_frhatv > 0) call post_data(CS%id_frhatv, CS%frhatv, CS%diag)
    if (CS%id_vhbt > 0) call post_data(CS%id_vhbt, vhbtav, CS%diag)
    if (CS%id_uhbt0 > 0) call post_data(CS%id_uhbt0, uhbt0(IsdB:IedB,jsd:jed), CS%diag)
    if (CS%id_vhbt0 > 0) call post_data(CS%id_vhbt0, vhbt0(isd:ied,JsdB:JedB), CS%diag)

    if (CS%id_frhatu1 > 0) call post_data(CS%id_frhatu1, CS%frhatu1, CS%diag)
    if (CS%id_frhatv1 > 0) call post_data(CS%id_frhatv1, CS%frhatv1, CS%diag)

    if (use_BT_cont) then
      if (CS%id_BTC_FA_u_EE > 0) call post_data(CS%id_BTC_FA_u_EE, BT_cont%FA_u_EE, CS%diag)
      if (CS%id_BTC_FA_u_E0 > 0) call post_data(CS%id_BTC_FA_u_E0, BT_cont%FA_u_E0, CS%diag)
      if (CS%id_BTC_FA_u_W0 > 0) call post_data(CS%id_BTC_FA_u_W0, BT_cont%FA_u_W0, CS%diag)
      if (CS%id_BTC_FA_u_WW > 0) call post_data(CS%id_BTC_FA_u_WW, BT_cont%FA_u_WW, CS%diag)
      if (CS%id_BTC_uBT_EE > 0) call post_data(CS%id_BTC_uBT_EE, BT_cont%uBT_EE, CS%diag)
      if (CS%id_BTC_uBT_WW > 0) call post_data(CS%id_BTC_uBT_WW, BT_cont%uBT_WW, CS%diag)
      if (CS%id_BTC_FA_u_rat0 > 0) then
        tmp_u(:,:) = 0.0
        do j=js,je ; do I=is-1,ie
          if ((G%mask2dCu(I,j) > 0.0) .and. (BT_cont%FA_u_W0(I,j) > 0.0)) then
            tmp_u(I,j) = (BT_cont%FA_u_E0(I,j)/ BT_cont%FA_u_W0(I,j))
          else
            tmp_u(I,j) = 1.0
          endif
        enddo ; enddo
        call post_data(CS%id_BTC_FA_u_rat0, tmp_u, CS%diag)
      endif
      if (CS%id_BTC_FA_v_NN > 0) call post_data(CS%id_BTC_FA_v_NN, BT_cont%FA_v_NN, CS%diag)
      if (CS%id_BTC_FA_v_N0 > 0) call post_data(CS%id_BTC_FA_v_N0, BT_cont%FA_v_N0, CS%diag)
      if (CS%id_BTC_FA_v_S0 > 0) call post_data(CS%id_BTC_FA_v_S0, BT_cont%FA_v_S0, CS%diag)
      if (CS%id_BTC_FA_v_SS > 0) call post_data(CS%id_BTC_FA_v_SS, BT_cont%FA_v_SS, CS%diag)
      if (CS%id_BTC_vBT_NN > 0) call post_data(CS%id_BTC_vBT_NN, BT_cont%vBT_NN, CS%diag)
      if (CS%id_BTC_vBT_SS > 0) call post_data(CS%id_BTC_vBT_SS, BT_cont%vBT_SS, CS%diag)
      if (CS%id_BTC_FA_v_rat0 > 0) then
        tmp_v(:,:) = 0.0
        do J=js-1,je ; do i=is,ie
          if ((G%mask2dCv(i,J) > 0.0) .and. (BT_cont%FA_v_S0(i,J) > 0.0)) then
            tmp_v(i,J) = (BT_cont%FA_v_N0(i,J)/ BT_cont%FA_v_S0(i,J))
          else
            tmp_v(i,J) = 1.0
          endif
        enddo ; enddo
        call post_data(CS%id_BTC_FA_v_rat0, tmp_v, CS%diag)
      endif
      if (CS%id_BTC_FA_h_rat0 > 0) then
        tmp_h(:,:) = 0.0
        do j=js,je ; do i=is,ie
          tmp_h(i,j) = 1.0
          if ((G%mask2dCu(I,j) > 0.0) .and. (BT_cont%FA_u_W0(I,j) > 0.0) .and. (BT_cont%FA_u_E0(I,j) > 0.0)) then
            if (BT_cont%FA_u_W0(I,j) > BT_cont%FA_u_E0(I,j)) then
              tmp_h(i,j) = max(tmp_h(i,j), (BT_cont%FA_u_W0(I,j)/ BT_cont%FA_u_E0(I,j)))
            else
              tmp_h(i,j) = max(tmp_h(i,j), (BT_cont%FA_u_E0(I,j)/ BT_cont%FA_u_W0(I,j)))
            endif
          endif
          if ((G%mask2dCu(I-1,j) > 0.0) .and. (BT_cont%FA_u_W0(I-1,j) > 0.0) .and. (BT_cont%FA_u_E0(I-1,j) > 0.0)) then
            if (BT_cont%FA_u_W0(I-1,j) > BT_cont%FA_u_E0(I-1,j)) then
              tmp_h(i,j) = max(tmp_h(i,j), (BT_cont%FA_u_W0(I-1,j)/ BT_cont%FA_u_E0(I-1,j)))
            else
              tmp_h(i,j) = max(tmp_h(i,j), (BT_cont%FA_u_E0(I-1,j)/ BT_cont%FA_u_W0(I-1,j)))
            endif
          endif
          if ((G%mask2dCv(i,J) > 0.0) .and. (BT_cont%FA_v_S0(i,J) > 0.0) .and. (BT_cont%FA_v_N0(i,J) > 0.0)) then
            if (BT_cont%FA_v_S0(i,J) > BT_cont%FA_v_N0(i,J)) then
              tmp_h(i,j) = max(tmp_h(i,j), (BT_cont%FA_v_S0(i,J)/ BT_cont%FA_v_N0(i,J)))
            else
              tmp_h(i,j) = max(tmp_h(i,j), (BT_cont%FA_v_N0(i,J)/ BT_cont%FA_v_S0(i,J)))
            endif
          endif
          if ((G%mask2dCv(i,J-1) > 0.0) .and. (BT_cont%FA_v_S0(i,J-1) > 0.0) .and. (BT_cont%FA_v_N0(i,J-1) > 0.0)) then
            if (BT_cont%FA_v_S0(i,J-1) > BT_cont%FA_v_N0(i,J-1)) then
              tmp_h(i,j) = max(tmp_h(i,j), (BT_cont%FA_v_S0(i,J-1)/ BT_cont%FA_v_N0(i,J-1)))
            else
              tmp_h(i,j) = max(tmp_h(i,j), (BT_cont%FA_v_N0(i,J-1)/ BT_cont%FA_v_S0(i,J-1)))
            endif
          endif
        enddo ; enddo
        call post_data(CS%id_BTC_FA_h_rat0, tmp_h, CS%diag)
      endif
    endif
  else
    if (CS%id_frhatu1 > 0) CS%frhatu1(:,:,:) = CS%frhatu(:,:,:)
    if (CS%id_frhatv1 > 0) CS%frhatv1(:,:,:) = CS%frhatv(:,:,:)
  endif

  if ((present(ADp)) .and. (associated(ADp%diag_hfrac_u))) then
    do k=1,nz ; do j=js,je ; do I=is-1,ie
      ADp%diag_hfrac_u(I,j,k) = CS%frhatu(I,j,k)
    enddo ; enddo ; enddo
  endif
  if ((present(ADp)) .and. (associated(ADp%diag_hfrac_v))) then
    do k=1,nz ; do J=js-1,je ; do i=is,ie
      ADp%diag_hfrac_v(i,J,k) = CS%frhatv(i,J,k)
    enddo ; enddo ; enddo
  endif

  if ((present(ADp)) .and. (present(BT_cont)) .and. (associated(ADp%diag_hu))) then
    do k=1,nz ; do j=js,je ; do I=is-1,ie
      ADp%diag_hu(I,j,k) = BT_cont%h_u(I,j,k)
    enddo ; enddo ; enddo
  endif
  if ((present(ADp)) .and. (present(BT_cont)) .and. (associated(ADp%diag_hv))) then
    do k=1,nz ; do J=js-1,je ; do i=is,ie
      ADp%diag_hv(i,J,k) = BT_cont%h_v(i,J,k)
    enddo ; enddo ; enddo
  endif

  if (G%nonblocking_updates) then
    if (find_etaav) call complete_group_pass(CS%pass_etaav, G%Domain)
    call complete_group_pass(CS%pass_ubta_uhbta, G%Domain)
  endif

end subroutine btstep

!> This subroutine automatically determines an optimal value for dtbt based
!! on some state of the ocean.
subroutine set_dtbt(G, GV, US, CS, eta, pbce, BT_cont, gtot_est, SSH_add)
  type(ocean_grid_type),        intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type),      intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),        intent(in)    :: US   !< A dimensional unit scaling type
  type(barotropic_CS),          pointer       :: CS   !< Barotropic control structure.
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: eta  !< The barotropic free surface
                                                      !! height anomaly or column mass anomaly [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), optional, intent(in) :: pbce  !< The baroclinic pressure
                                                      !! anomaly in each layer due to free surface
                                                      !! height anomalies [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
  type(BT_cont_type), optional, pointer       :: BT_cont  !< A structure with elements that describe
                                                      !! the effective open face areas as a
                                                      !! function of barotropic flow.
  real,               optional, intent(in)    :: gtot_est !< An estimate of the total gravitational
                                                      !! acceleration [L2 Z-1 T-2 ~> m s-2].
  real,               optional, intent(in)    :: SSH_add  !< An additional contribution to SSH to
                                                      !! provide a margin of error when
                                                      !! calculating the external wave speed [Z ~> m].

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    gtot_E, &     ! gtot_X is the effective total reduced gravity used to relate
    gtot_W, &     ! free surface height deviations to pressure forces (including
    gtot_N, &     ! GFS and baroclinic  contributions) in the barotropic momentum
    gtot_S        ! equations half a grid-point in the X-direction (X is N, S, E, or W)
                  ! from the thickness point [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
                  ! (See Hallberg, J Comp Phys 1997 for a discussion.)
  real, dimension(SZIBS_(G),SZJ_(G)) :: &
    Datu          ! Basin depth at u-velocity grid points times the y-grid
                  ! spacing [H L ~> m2 or kg m-1].
  real, dimension(SZI_(G),SZJBS_(G)) :: &
    Datv          ! Basin depth at v-velocity grid points times the x-grid
                  ! spacing [H L ~> m2 or kg m-1].
  real :: det_de  ! The partial derivative due to self-attraction and loading
                  ! of the reference geopotential with the sea surface height [nondim].
                  ! This is typically ~0.09 or less.
  real :: dgeo_de ! The constant of proportionality between geopotential and
                  ! sea surface height [nondim].  It is a nondimensional number of
                  ! order 1.  For stability, this may be made larger
                  ! than physical problem would suggest.
  real :: add_SSH ! An additional contribution to SSH to provide a margin of error
                  ! when calculating the external wave speed [Z ~> m].
  real :: min_max_dt2 ! The square of the minimum value of the largest stable barotropic
                      ! timesteps [T2 ~> s2]
  real :: dtbt_max    ! The maximum barotropic timestep [T ~> s]
  real :: Idt_max2    ! The squared inverse of the local maximum stable
                      ! barotropic time step [T-2 ~> s-2].
  logical :: use_BT_cont
  type(memory_size_type) :: MS

  character(len=200) :: mesg
  integer :: i, j, k, is, ie, js, je, nz

  if (.not.associated(CS)) call MOM_error(FATAL, &
      "set_dtbt: Module MOM_barotropic must be initialized before it is used.")
  if (.not.CS%split) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  MS%isdw = G%isd ; MS%iedw = G%ied ; MS%jsdw = G%jsd ; MS%jedw = G%jed

  if (.not.(present(pbce) .or. present(gtot_est))) call MOM_error(FATAL, &
      "set_dtbt: Either pbce or gtot_est must be present.")

  add_SSH = 0.0 ; if (present(SSH_add)) add_SSH = SSH_add

  use_BT_cont = .false.
  if (present(BT_cont)) use_BT_cont = (associated(BT_cont))

  if (use_BT_cont) then
    call BT_cont_to_face_areas(BT_cont, Datu, Datv, G, US, MS, 0, .true.)
  elseif (CS%Nonlinear_continuity .and. present(eta)) then
    call find_face_areas(Datu, Datv, G, GV, US, CS, MS, eta=eta, halo=0)
  else
    call find_face_areas(Datu, Datv, G, GV, US, CS, MS, halo=0, add_max=add_SSH)
  endif

  det_de = 0.0
  if (CS%tides) call tidal_forcing_sensitivity(G, CS%tides_CSp, det_de)
  dgeo_de = 1.0 + max(0.0, det_de + CS%G_extra)
  if (present(pbce)) then
    do j=js,je ; do i=is,ie
      gtot_E(i,j) = 0.0 ; gtot_W(i,j) = 0.0
      gtot_N(i,j) = 0.0 ; gtot_S(i,j) = 0.0
    enddo ; enddo
    do k=1,nz ; do j=js,je ; do i=is,ie
      gtot_E(i,j) = gtot_E(i,j) + pbce(i,j,k) * CS%frhatu(I,j,k)
      gtot_W(i,j) = gtot_W(i,j) + pbce(i,j,k) * CS%frhatu(I-1,j,k)
      gtot_N(i,j) = gtot_N(i,j) + pbce(i,j,k) * CS%frhatv(i,J,k)
      gtot_S(i,j) = gtot_S(i,j) + pbce(i,j,k) * CS%frhatv(i,J-1,k)
    enddo ; enddo ; enddo
  else
    do j=js,je ; do i=is,ie
      gtot_E(i,j) = gtot_est * GV%H_to_Z ; gtot_W(i,j) = gtot_est * GV%H_to_Z
      gtot_N(i,j) = gtot_est * GV%H_to_Z ; gtot_S(i,j) = gtot_est * GV%H_to_Z
    enddo ; enddo
  endif

  min_max_dt2 = 1.0e38*US%s_to_T**2  ! A huge value for the permissible timestep squared.
  do j=js,je ; do i=is,ie
    !   This is pretty accurate for gravity waves, but it is a conservative
    ! estimate since it ignores the stabilizing effect of the bottom drag.
    Idt_max2 = 0.5 * (1.0 + 2.0*CS%bebt) * (G%IareaT(i,j) * &
      ((gtot_E(i,j)*Datu(I,j)*G%IdxCu(I,j) + gtot_W(i,j)*Datu(I-1,j)*G%IdxCu(I-1,j)) + &
       (gtot_N(i,j)*Datv(i,J)*G%IdyCv(i,J) + gtot_S(i,j)*Datv(i,J-1)*G%IdyCv(i,J-1))) + &
      ((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
       (G%CoriolisBu(I-1,J)**2 + G%CoriolisBu(I,J-1)**2)) * CS%BT_Coriolis_scale**2 )
    if (Idt_max2 * min_max_dt2 > 1.0) min_max_dt2 = 1.0 / Idt_max2
  enddo ; enddo
  dtbt_max = sqrt(min_max_dt2 / dgeo_de)
  if (id_clock_sync > 0) call cpu_clock_begin(id_clock_sync)
  call min_across_PEs(dtbt_max)
  if (id_clock_sync > 0) call cpu_clock_end(id_clock_sync)

  CS%dtbt = CS%dtbt_fraction * dtbt_max
  CS%dtbt_max = dtbt_max
end subroutine set_dtbt

!> The following 4 subroutines apply the open boundary conditions.
!! This subroutine applies the open boundary conditions on barotropic
!! velocities and mass transports, as developed by Mehmet Ilicak.
subroutine apply_velocity_OBCs(OBC, ubt, vbt, uhbt, vhbt, ubt_trans, vbt_trans, eta, &
                               ubt_old, vbt_old, BT_OBC, G, MS, US, halo, dtbt, bebt, &
                               use_BT_cont, integral_BT_cont, dt_elapsed, Datu, Datv, &
                               BTCL_u, BTCL_v, uhbt0, vhbt0, ubt_int, vbt_int, uhbt_int, vhbt_int)
  type(ocean_OBC_type),                  pointer       :: OBC     !< An associated pointer to an OBC type.
  type(ocean_grid_type),                 intent(inout) :: G       !< The ocean's grid structure.
  type(memory_size_type),                intent(in)    :: MS      !< A type that describes the memory sizes of
                                                                  !! the argument arrays.
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(inout) :: ubt     !< the zonal barotropic velocity [L T-1 ~> m s-1].
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(inout) :: uhbt    !< the zonal barotropic transport
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(inout) :: ubt_trans !< The zonal barotropic velocity used in
                                                                  !! transport [L T-1 ~> m s-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(inout) :: vbt     !< The meridional barotropic velocity
                                                                  !! [L T-1 ~> m s-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(inout) :: vhbt    !< the meridional barotropic transport
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(inout) :: vbt_trans !< the meridional BT velocity used in
                                                                  !! transports [L T-1 ~> m s-1].
  real, dimension(SZIW_(MS),SZJW_(MS)),  intent(in)    :: eta     !< The barotropic free surface height anomaly or
                                                                  !! column mass anomaly [H ~> m or kg m-2].
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(in)    :: ubt_old !< The starting value of ubt in a barotropic
                                                                  !! step [L T-1 ~> m s-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(in)    :: vbt_old !< The starting value of vbt in a barotropic
                                                                  !! step [L T-1 ~> m s-1].
  type(BT_OBC_type),                     intent(in)    :: BT_OBC  !< A structure with the private barotropic arrays
                                                                  !! related to the open boundary conditions,
                                                                  !! set by set_up_BT_OBC.
  type(unit_scale_type),                 intent(in)    :: US      !< A dimensional unit scaling type
  integer,                               intent(in)    :: halo    !< The extra halo size to use here.
  real,                                  intent(in)    :: dtbt    !< The time step [T ~> s].
  real,                                  intent(in)    :: bebt    !< The fractional weighting of the future velocity
                                                                  !! in determining the transport.
  logical,                               intent(in)    :: use_BT_cont !< If true, use the BT_cont_types to calculate
                                                                  !! transports.
  logical,                               intent(in)    :: integral_BT_cont !< If true, update the barotropic continuity
                                                                  !! equation directly from the initial condition
                                                                  !! using the time-integrated barotropic velocity.
  real,                                  intent(in)    :: dt_elapsed !< The amount of time in the barotropic stepping
                                                                  !! that will have elapsed [T ~> s].
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(in)    :: Datu    !< A fixed estimate of the face areas at u points
                                                                  !! [H L ~> m2 or kg m-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(in)    :: Datv    !< A fixed estimate of the face areas at v points
                                                                  !! [H L ~> m2 or kg m-1].
  type(local_BT_cont_u_type), dimension(SZIBW_(MS),SZJW_(MS)), intent(in) :: BTCL_u !< Structure of information used
                                                                  !! for a dynamic estimate of the face areas at
                                                                  !! u-points.
  type(local_BT_cont_v_type), dimension(SZIW_(MS),SZJBW_(MS)), intent(in) :: BTCL_v !< Structure of information used
                                                                  !! for a dynamic estimate of the face areas at
                                                                  !! v-points.
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(in)    :: uhbt0   !< A correction to the zonal transport so that
                                                                  !! the barotropic functions agree with the sum
                                                                  !! of the layer transports
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(in)    :: vhbt0   !< A correction to the meridional transport so that
                                                                  !! the barotropic functions agree with the sum
                                                                  !! of the layer transports
                                                                  !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(in) :: ubt_int    !< The time-integrated zonal barotropic
                                                                  !! velocity before this update [L T-1 ~> m s-1].
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(in) :: uhbt_int   !< The time-integrated zonal barotropic
                                                                  !! transport [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(in) :: vbt_int    !< The time-integrated meridional barotropic
                                                                  !! velocity before this update [L T-1 ~> m s-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(in) :: vhbt_int   !< The time-integrated meridional barotropic
                                                                  !! transport [H L2 T-1 ~> m3 s-1 or kg s-1].

  ! Local variables
  real :: vel_prev    ! The previous velocity [L T-1 ~> m s-1].
  real :: vel_trans   ! The combination of the previous and current velocity
                      ! that does the mass transport [L T-1 ~> m s-1].
  real :: H_u         ! The total thickness at the u-point [H ~> m or kg m-2].
  real :: H_v         ! The total thickness at the v-point [H ~> m or kg m-2].
  real :: cfl         ! The CFL number at the point in question [nondim]
  real :: u_inlet     ! The zonal inflow velocity [L T-1 ~> m s-1]
  real :: v_inlet     ! The meridional inflow velocity [L T-1 ~> m s-1]
  real :: uhbt_int_new ! The updated time-integrated zonal transport [H L2 ~> m3]
  real :: vhbt_int_new ! The updated time-integrated meridional transport [H L2 ~> m3]
  real :: h_in        ! The inflow thickess [H ~> m or kg m-2].
  real :: cff, Cx, Cy, tau
  real :: dhdt, dhdx, dhdy
  real :: Idtbt       ! The inverse of the barotropic time step [T-1 ~> s-1]
  integer :: i, j, is, ie, js, je
  real, dimension(SZIB_(G),SZJB_(G)) :: grad
  real, parameter :: eps = 1.0e-20
  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo

  if (.not.(BT_OBC%apply_u_OBCs .or. BT_OBC%apply_v_OBCs)) return

  Idtbt = 1.0 / dtbt

  if (BT_OBC%apply_u_OBCs) then
    do j=js,je ; do I=is-1,ie ; if (OBC%segnum_u(I,j) /= OBC_NONE) then
      if (OBC%segment(OBC%segnum_u(I,j))%specified) then
        uhbt(I,j) = BT_OBC%uhbt(I,j)
        ubt(I,j) = BT_OBC%ubt_outer(I,j)
        vel_trans = ubt(I,j)
      elseif (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_E) then
        if (OBC%segment(OBC%segnum_u(I,j))%Flather) then
          cfl = dtbt * BT_OBC%Cg_u(I,j) * G%IdxCu(I,j) ! CFL
          u_inlet = cfl*ubt_old(I-1,j) + (1.0-cfl)*ubt_old(I,j)  ! Valid for cfl<1
          h_in = eta(i,j) + (0.5-cfl)*(eta(i,j)-eta(i-1,j))      ! internal
          H_u = BT_OBC%H_u(I,j)
          vel_prev = ubt(I,j)
          ubt(I,j) = 0.5*((u_inlet + BT_OBC%ubt_outer(I,j)) + &
              (BT_OBC%Cg_u(I,j)/H_u) * (h_in-BT_OBC%eta_outer_u(I,j)))
          vel_trans = (1.0-bebt)*vel_prev + bebt*ubt(I,j)
        elseif (OBC%segment(OBC%segnum_u(I,j))%gradient) then
          ubt(I,j) = ubt(I-1,j)
          vel_trans = ubt(I,j)
        endif
      elseif (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_W) then
        if (OBC%segment(OBC%segnum_u(I,j))%Flather) then
          cfl = dtbt * BT_OBC%Cg_u(I,j) * G%IdxCu(I,j) ! CFL
          u_inlet = cfl*ubt_old(I+1,j) + (1.0-cfl)*ubt_old(I,j)  ! Valid for cfl<1
          h_in = eta(i+1,j) + (0.5-cfl)*(eta(i+1,j)-eta(i+2,j))  ! external

          H_u = BT_OBC%H_u(I,j)
          vel_prev = ubt(I,j)
          ubt(I,j) = 0.5*((u_inlet + BT_OBC%ubt_outer(I,j)) + &
              (BT_OBC%Cg_u(I,j)/H_u) * (BT_OBC%eta_outer_u(I,j)-h_in))

          vel_trans = (1.0-bebt)*vel_prev + bebt*ubt(I,j)
        elseif (OBC%segment(OBC%segnum_u(I,j))%gradient) then
          ubt(I,j) = ubt(I+1,j)
          vel_trans = ubt(I,j)
        endif
      endif

      if (.not. OBC%segment(OBC%segnum_u(I,j))%specified) then
        if (integral_BT_cont) then
          uhbt_int_new = find_uhbt(ubt_int(I,j) + dtbt*vel_trans, BTCL_u(I,j)) + &
                         dt_elapsed*uhbt0(I,j)
          uhbt(I,j) = (uhbt_int_new - uhbt_int(I,j)) * Idtbt
        elseif (use_BT_cont) then
          uhbt(I,j) = find_uhbt(vel_trans, BTCL_u(I,j)) + uhbt0(I,j)
        else
          uhbt(I,j) = Datu(I,j)*vel_trans + uhbt0(I,j)
        endif
      endif

      ubt_trans(I,j) = vel_trans
    endif ; enddo ; enddo
  endif

  if (BT_OBC%apply_v_OBCs) then
    do J=js-1,je ; do i=is,ie ; if (OBC%segnum_v(i,J) /= OBC_NONE) then
      if (OBC%segment(OBC%segnum_v(i,J))%specified) then
        vhbt(i,J) = BT_OBC%vhbt(i,J)
        vbt(i,J) = BT_OBC%vbt_outer(i,J)
        vel_trans = vbt(i,J)
      elseif (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_N) then
        if (OBC%segment(OBC%segnum_v(i,J))%Flather) then
          cfl = dtbt * BT_OBC%Cg_v(i,J) * G%IdyCv(i,J) ! CFL
          v_inlet = cfl*vbt_old(i,J-1) + (1.0-cfl)*vbt_old(i,J)  ! Valid for cfl<1
          h_in = eta(i,j) + (0.5-cfl)*(eta(i,j)-eta(i,j-1))      ! internal

          H_v = BT_OBC%H_v(i,J)
          vel_prev = vbt(i,J)
          vbt(i,J) = 0.5*((v_inlet + BT_OBC%vbt_outer(i,J)) + &
              (BT_OBC%Cg_v(i,J)/H_v) * (h_in-BT_OBC%eta_outer_v(i,J)))

          vel_trans = (1.0-bebt)*vel_prev + bebt*vbt(i,J)
        elseif (OBC%segment(OBC%segnum_v(i,J))%gradient) then
          vbt(i,J) = vbt(i,J-1)
          vel_trans = vbt(i,J)
        endif
      elseif (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_S) then
        if (OBC%segment(OBC%segnum_v(i,J))%Flather) then
          cfl = dtbt * BT_OBC%Cg_v(i,J) * G%IdyCv(i,J) ! CFL
          v_inlet = cfl*vbt_old(i,J+1) + (1.0-cfl)*vbt_old(i,J)  ! Valid for cfl <1
          h_in = eta(i,j+1) + (0.5-cfl)*(eta(i,j+1)-eta(i,j+2))  ! internal

          H_v = BT_OBC%H_v(i,J)
          vel_prev = vbt(i,J)
          vbt(i,J) = 0.5*((v_inlet + BT_OBC%vbt_outer(i,J)) + &
              (BT_OBC%Cg_v(i,J)/H_v) * (BT_OBC%eta_outer_v(i,J)-h_in))

          vel_trans = (1.0-bebt)*vel_prev + bebt*vbt(i,J)
        elseif (OBC%segment(OBC%segnum_v(i,J))%gradient) then
          vbt(i,J) = vbt(i,J+1)
          vel_trans = vbt(i,J)
        endif
      endif

      if (.not. OBC%segment(OBC%segnum_v(i,J))%specified) then
        if (integral_BT_cont) then
          vhbt_int_new = find_vhbt(vbt_int(i,J) + dtbt*vel_trans, BTCL_v(i,J)) + &
                         dt_elapsed*vhbt0(i,J)
          vhbt(i,J) = (vhbt_int_new - vhbt_int(i,J)) * Idtbt
        elseif (use_BT_cont) then
          vhbt(i,J) = find_vhbt(vel_trans, BTCL_v(i,J)) + vhbt0(i,J)
        else
          vhbt(i,J) = vel_trans*Datv(i,J) + vhbt0(i,J)
        endif
      endif

      vbt_trans(i,J) = vel_trans
    endif ; enddo ; enddo
  endif

end subroutine apply_velocity_OBCs

!> This subroutine sets up the private structure used to apply the open
!! boundary conditions, as developed by Mehmet Ilicak.
subroutine set_up_BT_OBC(OBC, eta, BT_OBC, BT_Domain, G, GV, US, MS, halo, use_BT_cont, &
                         integral_BT_cont, dt_baroclinic, Datu, Datv, BTCL_u, BTCL_v)
  type(ocean_OBC_type),                  pointer       :: OBC    !< An associated pointer to an OBC type.
  type(memory_size_type),                intent(in)    :: MS     !< A type that describes the memory sizes of the
                                                                 !! argument arrays.
  real, dimension(SZIW_(MS),SZJW_(MS)),  intent(in)    :: eta    !< The barotropic free surface height anomaly or
                                                                 !! column mass anomaly [H ~> m or kg m-2].
  type(BT_OBC_type),                     intent(inout) :: BT_OBC !< A structure with the private barotropic arrays
                                                                 !! related to the open boundary conditions,
                                                                 !! set by set_up_BT_OBC.
  type(MOM_domain_type),                 intent(inout) :: BT_Domain !< MOM_domain_type associated with wide arrays
  type(ocean_grid_type),                 intent(inout) :: G      !< The ocean's grid structure.
  type(verticalGrid_type),               intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),                 intent(in)    :: US     !< A dimensional unit scaling type
  integer,                               intent(in)    :: halo   !< The extra halo size to use here.
  logical,                               intent(in)    :: use_BT_cont !< If true, use the BT_cont_types to calculate
                                                                 !! transports.
  logical,                               intent(in)    :: integral_BT_cont !< If true, update the barotropic continuity
                                                                 !! equation directly from the initial condition
                                                                 !! using the time-integrated barotropic velocity.
  real,                                  intent(in)    :: dt_baroclinic !< The baroclinic timestep for this cycle of
                                                                 !! updates to the barotropic solver [T ~> s]
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(in)    :: Datu   !< A fixed estimate of the face areas at u points
                                                                 !! [H L ~> m2 or kg m-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(in)    :: Datv   !< A fixed estimate of the face areas at v points
                                                                 !! [H L ~> m2 or kg m-1].
  type(local_BT_cont_u_type), dimension(SZIBW_(MS),SZJW_(MS)), intent(in) :: BTCL_u !< Structure of information used
                                                                 !! for a dynamic estimate of the face areas at
                                                                 !! u-points.
  type(local_BT_cont_v_type), dimension(SZIW_(MS),SZJBW_(MS)), intent(in) :: BTCL_v !< Structure of information used
                                                                 !! for a dynamic estimate of the face areas at
                                                                 !! v-points.

  ! Local variables
  real :: I_dt      ! The inverse of the time interval of this call [T-1 ~> s-1].
  integer :: i, j, k, is, ie, js, je, n, nz, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  integer :: isdw, iedw, jsdw, jedw
  logical :: OBC_used
  type(OBC_segment_type), pointer  :: segment !< Open boundary segment

  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  isdw = MS%isdw ; iedw = MS%iedw ; jsdw = MS%jsdw ; jedw = MS%jedw

  I_dt = 1.0 / dt_baroclinic

  if ((isdw < isd) .or. (jsdw < jsd)) then
    call MOM_error(FATAL, "set_up_BT_OBC: Open boundary conditions are not "//&
                           "yet fully implemented with wide barotropic halos.")
  endif

  if (.not. BT_OBC%is_alloced) then
    allocate(BT_OBC%Cg_u(isdw-1:iedw,jsdw:jedw))        ; BT_OBC%Cg_u(:,:) = 0.0
    allocate(BT_OBC%H_u(isdw-1:iedw,jsdw:jedw))         ; BT_OBC%H_u(:,:) = 0.0
    allocate(BT_OBC%uhbt(isdw-1:iedw,jsdw:jedw))        ; BT_OBC%uhbt(:,:) = 0.0
    allocate(BT_OBC%ubt_outer(isdw-1:iedw,jsdw:jedw))   ; BT_OBC%ubt_outer(:,:) = 0.0
    allocate(BT_OBC%eta_outer_u(isdw-1:iedw,jsdw:jedw)) ; BT_OBC%eta_outer_u(:,:) = 0.0

    allocate(BT_OBC%Cg_v(isdw:iedw,jsdw-1:jedw))        ; BT_OBC%Cg_v(:,:) = 0.0
    allocate(BT_OBC%H_v(isdw:iedw,jsdw-1:jedw))         ; BT_OBC%H_v(:,:) = 0.0
    allocate(BT_OBC%vhbt(isdw:iedw,jsdw-1:jedw))        ; BT_OBC%vhbt(:,:) = 0.0
    allocate(BT_OBC%vbt_outer(isdw:iedw,jsdw-1:jedw))   ; BT_OBC%vbt_outer(:,:) = 0.0
    allocate(BT_OBC%eta_outer_v(isdw:iedw,jsdw-1:jedw)) ; BT_OBC%eta_outer_v(:,:)=0.0
    BT_OBC%is_alloced = .true.
    call create_group_pass(BT_OBC%pass_uv, BT_OBC%ubt_outer, BT_OBC%vbt_outer, BT_Domain)
    call create_group_pass(BT_OBC%pass_uhvh, BT_OBC%uhbt, BT_OBC%vhbt, BT_Domain)
    call create_group_pass(BT_OBC%pass_eta_outer, BT_OBC%eta_outer_u, BT_OBC%eta_outer_v, BT_Domain,To_All+Scalar_Pair)
    call create_group_pass(BT_OBC%pass_h, BT_OBC%H_u, BT_OBC%H_v, BT_Domain,To_All+Scalar_Pair)
    call create_group_pass(BT_OBC%pass_cg, BT_OBC%Cg_u, BT_OBC%Cg_v, BT_Domain,To_All+Scalar_Pair)
  endif

  if (BT_OBC%apply_u_OBCs) then
    if (OBC%specified_u_BCs_exist_globally) then
      do n = 1, OBC%number_of_segments
        segment => OBC%segment(n)
        if (segment%is_E_or_W .and. segment%specified) then
          do j=segment%HI%jsd,segment%HI%jed ; do I=segment%HI%IsdB,segment%HI%IedB
            BT_OBC%uhbt(I,j) = 0.
          enddo ; enddo
          do k=1,nz ; do j=segment%HI%jsd,segment%HI%jed ; do I=segment%HI%IsdB,segment%HI%IedB
            BT_OBC%uhbt(I,j) = BT_OBC%uhbt(I,j) + segment%normal_trans(I,j,k)
          enddo ; enddo ; enddo
        endif
      enddo
    endif
    do j=js,je ; do I=is-1,ie ; if (OBC%segnum_u(I,j) /= OBC_NONE) then
      ! Can this go in segment loop above? Is loop above wrong for wide halos??
      if (OBC%segment(OBC%segnum_u(I,j))%specified) then
        if (integral_BT_cont) then
          BT_OBC%ubt_outer(I,j) = uhbt_to_ubt(BT_OBC%uhbt(I,j)*dt_baroclinic, BTCL_u(I,j)) * I_dt
        elseif (use_BT_cont) then
          BT_OBC%ubt_outer(I,j) = uhbt_to_ubt(BT_OBC%uhbt(I,j), BTCL_u(I,j))
        else
          if (Datu(I,j) > 0.0) BT_OBC%ubt_outer(I,j) = BT_OBC%uhbt(I,j) / Datu(I,j)
        endif
      else  ! This is assuming Flather as only other option
        if (GV%Boussinesq) then
          if (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_E) then
            BT_OBC%H_u(I,j) = G%bathyT(i,j)*GV%Z_to_H + eta(i,j)
          elseif (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_W) then
            BT_OBC%H_u(I,j) = G%bathyT(i+1,j)*GV%Z_to_H + eta(i+1,j)
          endif
        else
          if (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_E) then
            BT_OBC%H_u(I,j) = eta(i,j)
          elseif (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_W) then
            BT_OBC%H_u(I,j) = eta(i+1,j)
          endif
        endif
        BT_OBC%Cg_u(I,j) = SQRT(GV%g_prime(1) * GV%H_to_Z*BT_OBC%H_u(i,j))
      endif
    endif ; enddo ; enddo
    if (OBC%Flather_u_BCs_exist_globally) then
      do n = 1, OBC%number_of_segments
        segment => OBC%segment(n)
        if (segment%is_E_or_W .and. segment%Flather) then
          do j=segment%HI%jsd,segment%HI%jed ; do I=segment%HI%IsdB,segment%HI%IedB
            BT_OBC%ubt_outer(I,j) = segment%normal_vel_bt(I,j)
            BT_OBC%eta_outer_u(I,j) = segment%eta(I,j)
          enddo ; enddo
        endif
      enddo
    endif
  endif

  if (BT_OBC%apply_v_OBCs) then
    if (OBC%specified_v_BCs_exist_globally) then
      do n = 1, OBC%number_of_segments
        segment => OBC%segment(n)
        if (segment%is_N_or_S .and. segment%specified) then
          do J=segment%HI%JsdB,segment%HI%JedB ; do i=segment%HI%isd,segment%HI%ied
            BT_OBC%vhbt(i,J) = 0.
          enddo ; enddo
          do k=1,nz ; do J=segment%HI%JsdB,segment%HI%JedB ; do i=segment%HI%isd,segment%HI%ied
            BT_OBC%vhbt(i,J) = BT_OBC%vhbt(i,J) + segment%normal_trans(i,J,k)
          enddo ; enddo ; enddo
        endif
      enddo
    endif
    do J=js-1,je ; do i=is,ie ; if (OBC%segnum_v(i,J) /= OBC_NONE) then
      ! Can this go in segment loop above? Is loop above wrong for wide halos??
      if (OBC%segment(OBC%segnum_v(i,J))%specified) then
        if (integral_BT_cont) then
          BT_OBC%vbt_outer(i,J) = vhbt_to_vbt(BT_OBC%vhbt(i,J)*dt_baroclinic, BTCL_v(i,J)) * I_dt
        elseif (use_BT_cont) then
          BT_OBC%vbt_outer(i,J) = vhbt_to_vbt(BT_OBC%vhbt(i,J), BTCL_v(i,J))
        else
          if (Datv(i,J) > 0.0) BT_OBC%vbt_outer(i,J) = BT_OBC%vhbt(i,J) / Datv(i,J)
        endif
      else  ! This is assuming Flather as only other option
        if (GV%Boussinesq) then
          if (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_N) then
            BT_OBC%H_v(i,J) = G%bathyT(i,j)*GV%Z_to_H + eta(i,j)
          elseif (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_S) then
            BT_OBC%H_v(i,J) = G%bathyT(i,j+1)*GV%Z_to_H + eta(i,j+1)
          endif
        else
          if (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_N) then
            BT_OBC%H_v(i,J) = eta(i,j)
          elseif (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_S) then
            BT_OBC%H_v(i,J) = eta(i,j+1)
          endif
        endif
        BT_OBC%Cg_v(i,J) = SQRT(GV%g_prime(1) * GV%H_to_Z*BT_OBC%H_v(i,J))
      endif
    endif ; enddo ; enddo
    if (OBC%Flather_v_BCs_exist_globally) then
      do n = 1, OBC%number_of_segments
        segment => OBC%segment(n)
        if (segment%is_N_or_S .and. segment%Flather) then
          do J=segment%HI%JsdB,segment%HI%JedB ; do i=segment%HI%isd,segment%HI%ied
            BT_OBC%vbt_outer(i,J) = segment%normal_vel_bt(i,J)
            BT_OBC%eta_outer_v(i,J) = segment%eta(i,J)
          enddo ; enddo
        endif
      enddo
    endif
  endif

  call do_group_pass(BT_OBC%pass_uv, BT_Domain)
  call do_group_pass(BT_OBC%pass_uhvh, BT_Domain)
  call do_group_pass(BT_OBC%pass_eta_outer, BT_Domain)
  call do_group_pass(BT_OBC%pass_h, BT_Domain)
  call do_group_pass(BT_OBC%pass_cg, BT_Domain)

end subroutine set_up_BT_OBC

!> Clean up the BT_OBC memory.
subroutine destroy_BT_OBC(BT_OBC)
  type(BT_OBC_type), intent(inout) :: BT_OBC !< A structure with the private barotropic arrays
                                             !! related to the open boundary conditions,
                                             !! set by set_up_BT_OBC.

  if (BT_OBC%is_alloced) then
    deallocate(BT_OBC%Cg_u)
    deallocate(BT_OBC%H_u)
    deallocate(BT_OBC%uhbt)
    deallocate(BT_OBC%ubt_outer)
    deallocate(BT_OBC%eta_outer_u)

    deallocate(BT_OBC%Cg_v)
    deallocate(BT_OBC%H_v)
    deallocate(BT_OBC%vhbt)
    deallocate(BT_OBC%vbt_outer)
    deallocate(BT_OBC%eta_outer_v)
    BT_OBC%is_alloced = .false.
  endif
end subroutine destroy_BT_OBC

!> btcalc calculates the barotropic velocities from the full velocity and
!! thickness fields, determines the fraction of the total water column in each
!! layer at velocity points, and determines a corrective fictitious mass source
!! that will drive the barotropic estimate of the free surface height toward the
!! baroclinic estimate.
subroutine btcalc(h, G, GV, CS, h_u, h_v, may_use_default, OBC)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2].
  type(barotropic_CS),     pointer       :: CS   !< The control structure returned by a previous
                                                 !! call to barotropic_init.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: h_u  !< The specified thicknesses at u-points [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 optional, intent(in)    :: h_v  !< The specified thicknesses at v-points [H ~> m or kg m-2].
  logical,       optional, intent(in)    :: may_use_default !< An optional logical argument
                                                 !! to indicate that the default velocity point
                                                 !! thicknesses may be used for this particular
                                                 !! calculation, even though the setting of
                                                 !! CS%hvel_scheme would usually require that h_u
                                                 !! and h_v be passed in.
  type(ocean_OBC_type), optional, pointer :: OBC !< Open boundary control structure.

  ! Local variables
  real :: hatutot(SZIB_(G))    ! The sum of the layer thicknesses interpolated to u points [H ~> m or kg m-2].
  real :: hatvtot(SZI_(G))     ! The sum of the layer thicknesses interpolated to v points [H ~> m or kg m-2].
  real :: Ihatutot(SZIB_(G))   ! Ihatutot is the inverse of hatutot [H-1 ~> m-1 or m2 kg-1].
  real :: Ihatvtot(SZI_(G))    ! Ihatvtot is the inverse of hatvtot [H-1 ~> m-1 or m2 kg-1].
  real :: h_arith              ! The arithmetic mean thickness [H ~> m or kg m-2].
  real :: h_harm               ! The harmonic mean thicknesses [H ~> m or kg m-2].
  real :: h_neglect            ! A thickness that is so small it is usually lost
                               ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: wt_arith             ! The nondimensional weight for the arithmetic mean thickness.
                               ! The harmonic mean uses a weight of (1 - wt_arith).
  real :: Rh                   ! A ratio of summed thicknesses, nondim.
  real :: e_u(SZIB_(G),SZK_(GV)+1) !   The interface heights at u-velocity and
  real :: e_v(SZI_(G),SZK_(GV)+1)  ! v-velocity points [H ~> m or kg m-2].
  real :: D_shallow_u(SZI_(G)) ! The shallower of the adjacent depths [H ~> m or kg m-2].
  real :: D_shallow_v(SZIB_(G))! The shallower of the adjacent depths [H ~> m or kg m-2].
  real :: htot                 ! The sum of the layer thicknesses [H ~> m or kg m-2].
  real :: Ihtot                ! The inverse of htot [H-1 ~> m-1 or m2 kg-1].

  logical :: use_default, test_dflt, apply_OBCs
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, i, j, k
  integer :: iss, ies, n

!    This section interpolates thicknesses onto u & v grid points with the
! second order accurate estimate h = 2*(h+ * h-)/(h+ + h-).
  if (.not.associated(CS)) call MOM_error(FATAL, &
      "btcalc: Module MOM_barotropic must be initialized before it is used.")
  if (.not.CS%split) return

  use_default = .false.
  test_dflt = .false. ; if (present(may_use_default)) test_dflt = may_use_default

  if (test_dflt) then
    if (.not.((present(h_u) .and. present(h_v)) .or. &
              (CS%hvel_scheme == HARMONIC) .or. (CS%hvel_scheme == HYBRID) .or.&
              (CS%hvel_scheme == ARITHMETIC))) use_default = .true.
  else
    if (.not.((present(h_u) .and. present(h_v)) .or. &
              (CS%hvel_scheme == HARMONIC) .or. (CS%hvel_scheme == HYBRID) .or.&
              (CS%hvel_scheme == ARITHMETIC))) call MOM_error(FATAL, &
        "btcalc: Inconsistent settings of optional arguments and hvel_scheme.")
  endif

  apply_OBCs = .false.
  if (present(OBC)) then ; if (associated(OBC)) then ; if (OBC%OBC_pe) then
    ! Some open boundary condition points might be in this processor's symmetric
    ! computational domain.
    apply_OBCs = (OBC%number_of_segments > 0)
  endif ; endif ; endif

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  h_neglect = GV%H_subroundoff

  !   This estimates the fractional thickness of each layer at the velocity
  ! points, using a harmonic mean estimate.
!$OMP parallel do default(none) shared(is,ie,js,je,nz,h_u,CS,h_neglect,h,use_default,G,GV) &
!$OMP                          private(hatutot,Ihatutot,e_u,D_shallow_u,h_arith,h_harm,wt_arith)

  do j=js,je
    if (present(h_u)) then
      do I=is-1,ie ; hatutot(I) = h_u(I,j,1) ; enddo
      do k=2,nz ; do I=is-1,ie
        hatutot(I) = hatutot(I) + h_u(I,j,k)
      enddo ; enddo
      do I=is-1,ie ; Ihatutot(I) = G%mask2dCu(I,j) / (hatutot(I) + h_neglect) ; enddo
      do k=1,nz ; do I=is-1,ie
        CS%frhatu(I,j,k) = h_u(I,j,k) * Ihatutot(I)
      enddo ; enddo
    else
      if (CS%hvel_scheme == ARITHMETIC) then
        do I=is-1,ie
          CS%frhatu(I,j,1) = 0.5 * (h(i+1,j,1) + h(i,j,1))
          hatutot(I) = CS%frhatu(I,j,1)
        enddo
        do k=2,nz ; do I=is-1,ie
          CS%frhatu(I,j,k) = 0.5 * (h(i+1,j,k) + h(i,j,k))
          hatutot(I) = hatutot(I) + CS%frhatu(I,j,k)
        enddo ; enddo
      elseif (CS%hvel_scheme == HYBRID .or. use_default) then
        do I=is-1,ie
          e_u(I,nz+1) = -0.5 * GV%Z_to_H * (G%bathyT(i+1,j) + G%bathyT(i,j))
          D_shallow_u(I) = -GV%Z_to_H * min(G%bathyT(i+1,j), G%bathyT(i,j))
          hatutot(I) = 0.0
        enddo
        do k=nz,1,-1 ; do I=is-1,ie
          e_u(I,K) = e_u(I,K+1) + 0.5 * (h(i+1,j,k) + h(i,j,k))
          h_arith = 0.5 * (h(i+1,j,k) + h(i,j,k))
          if (e_u(I,K+1) >= D_shallow_u(I)) then
            CS%frhatu(I,j,k) = h_arith
          else
            h_harm = (h(i+1,j,k) * h(i,j,k)) / (h_arith + h_neglect)
            if (e_u(I,K) <= D_shallow_u(I)) then
              CS%frhatu(I,j,k) = h_harm
            else
              wt_arith = (e_u(I,K) - D_shallow_u(I)) / (h_arith + h_neglect)
              CS%frhatu(I,j,k) = wt_arith*h_arith + (1.0-wt_arith)*h_harm
            endif
          endif
          hatutot(I) = hatutot(I) + CS%frhatu(I,j,k)
        enddo ; enddo
      elseif (CS%hvel_scheme == HARMONIC) then
        do I=is-1,ie
          CS%frhatu(I,j,1) = 2.0*(h(i+1,j,1) * h(i,j,1)) / &
                             ((h(i+1,j,1) + h(i,j,1)) + h_neglect)
          hatutot(I) = CS%frhatu(I,j,1)
        enddo
        do k=2,nz ; do I=is-1,ie
          CS%frhatu(I,j,k) = 2.0*(h(i+1,j,k) * h(i,j,k)) / &
                             ((h(i+1,j,k) + h(i,j,k)) + h_neglect)
          hatutot(I) = hatutot(I) + CS%frhatu(I,j,k)
        enddo ; enddo
      endif
      do I=is-1,ie ; Ihatutot(I) = G%mask2dCu(I,j) / (hatutot(I) + h_neglect) ; enddo
      do k=1,nz ; do I=is-1,ie
        CS%frhatu(I,j,k) = CS%frhatu(I,j,k) * Ihatutot(I)
      enddo ; enddo
    endif
  enddo

!$OMP parallel do default(none) shared(is,ie,js,je,nz,CS,G,GV,h_v,h_neglect,h,use_default) &
!$OMP                          private(hatvtot,Ihatvtot,e_v,D_shallow_v,h_arith,h_harm,wt_arith)
  do J=js-1,je
    if (present(h_v)) then
      do i=is,ie ; hatvtot(i) = h_v(i,J,1) ; enddo
      do k=2,nz ; do i=is,ie
        hatvtot(i) = hatvtot(i) + h_v(i,J,k)
      enddo ; enddo
      do i=is,ie ; Ihatvtot(i) = G%mask2dCv(i,J) / (hatvtot(i) + h_neglect) ; enddo
      do k=1,nz ; do i=is,ie
        CS%frhatv(i,J,k) = h_v(i,J,k) * Ihatvtot(i)
      enddo ; enddo
    else
      if (CS%hvel_scheme == ARITHMETIC) then
        do i=is,ie
          CS%frhatv(i,J,1) = 0.5 * (h(i,j+1,1) + h(i,j,1))
          hatvtot(i) = CS%frhatv(i,J,1)
        enddo
        do k=2,nz ; do i=is,ie
          CS%frhatv(i,J,k) = 0.5 * (h(i,j+1,k) + h(i,j,k))
          hatvtot(i) = hatvtot(i) + CS%frhatv(i,J,k)
        enddo ; enddo
      elseif (CS%hvel_scheme == HYBRID .or. use_default) then
        do i=is,ie
          e_v(i,nz+1) = -0.5 * GV%Z_to_H * (G%bathyT(i,j+1) + G%bathyT(i,j))
          D_shallow_v(I) = -GV%Z_to_H * min(G%bathyT(i,j+1), G%bathyT(i,j))
          hatvtot(I) = 0.0
        enddo
        do k=nz,1,-1 ; do i=is,ie
          e_v(i,K) = e_v(i,K+1) + 0.5 * (h(i,j+1,k) + h(i,j,k))
          h_arith = 0.5 * (h(i,j+1,k) + h(i,j,k))
          if (e_v(i,K+1) >= D_shallow_v(i)) then
            CS%frhatv(i,J,k) = h_arith
          else
            h_harm = (h(i,j+1,k) * h(i,j,k)) / (h_arith + h_neglect)
            if (e_v(i,K) <= D_shallow_v(i)) then
              CS%frhatv(i,J,k) = h_harm
            else
              wt_arith = (e_v(i,K) - D_shallow_v(i)) / (h_arith + h_neglect)
              CS%frhatv(i,J,k) = wt_arith*h_arith + (1.0-wt_arith)*h_harm
            endif
          endif
          hatvtot(i) = hatvtot(i) + CS%frhatv(i,J,k)
        enddo ; enddo
      elseif (CS%hvel_scheme == HARMONIC) then
        do i=is,ie
          CS%frhatv(i,J,1) = 2.0*(h(i,j+1,1) * h(i,j,1)) / &
                             ((h(i,j+1,1) + h(i,j,1)) + h_neglect)
          hatvtot(i) = CS%frhatv(i,J,1)
        enddo
        do k=2,nz ; do i=is,ie
          CS%frhatv(i,J,k) = 2.0*(h(i,j+1,k) * h(i,j,k)) / &
                             ((h(i,j+1,k) + h(i,j,k)) + h_neglect)
          hatvtot(i) = hatvtot(i) + CS%frhatv(i,J,k)
        enddo ; enddo
      endif
      do i=is,ie ; Ihatvtot(i) = G%mask2dCv(i,J) / (hatvtot(i) + h_neglect) ; enddo
      do k=1,nz ; do i=is,ie
        CS%frhatv(i,J,k) = CS%frhatv(i,J,k) * Ihatvtot(i)
      enddo ; enddo
    endif
  enddo

  if (apply_OBCs) then ; do n=1,OBC%number_of_segments ! Test for segment type?
    if (.not. OBC%segment(n)%on_pe) cycle
    if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
      J = OBC%segment(n)%HI%JsdB
      if ((J >= js-1) .and. (J <= je)) then
        iss = max(is,OBC%segment(n)%HI%isd) ; ies = min(ie,OBC%segment(n)%HI%ied)
        do i=iss,ies ; hatvtot(i) = h(i,j,1) ; enddo
        do k=2,nz ; do i=iss,ies
          hatvtot(i) = hatvtot(i) + h(i,j,k)
        enddo ; enddo
        do i=iss,ies
          Ihatvtot(i) = G%mask2dCv(i,J) / (hatvtot(i) + h_neglect)
        enddo
        do k=1,nz ; do i=iss,ies
          CS%frhatv(i,J,k) = h(i,j,k) * Ihatvtot(i)
        enddo ; enddo
      endif
    elseif (OBC%segment(n)%direction == OBC_DIRECTION_S) then
      J = OBC%segment(n)%HI%JsdB
      if ((J >= js-1) .and. (J <= je)) then
        iss = max(is,OBC%segment(n)%HI%isd) ; ies = min(ie,OBC%segment(n)%HI%ied)
        do i=iss,ies ; hatvtot(i) = h(i,j+1,1) ; enddo
        do k=2,nz ; do i=iss,ies
          hatvtot(i) = hatvtot(i) + h(i,j+1,k)
        enddo ; enddo
        do i=iss,ies
          Ihatvtot(i) = G%mask2dCv(i,J) / (hatvtot(i) + h_neglect)
        enddo
        do k=1,nz ; do i=iss,ies
          CS%frhatv(i,J,k) = h(i,j+1,k) * Ihatvtot(i)
        enddo ; enddo
      endif
    elseif (OBC%segment(n)%direction == OBC_DIRECTION_E) then
      I = OBC%segment(n)%HI%IsdB
      if ((I >= is-1) .and. (I <= ie)) then
        do j = max(js,OBC%segment(n)%HI%jsd), min(je,OBC%segment(n)%HI%jed)
          htot = h(i,j,1)
          do k=2,nz ; htot = htot + h(i,j,k) ; enddo
          Ihtot = G%mask2dCu(I,j) / (htot + h_neglect)
          do k=1,nz ; CS%frhatu(I,j,k) = h(i,j,k) * Ihtot ; enddo
        enddo
      endif
    elseif (OBC%segment(n)%direction == OBC_DIRECTION_W) then
      I = OBC%segment(n)%HI%IsdB
      if ((I >= is-1) .and. (I <= ie)) then
        do j = max(js,OBC%segment(n)%HI%jsd), min(je,OBC%segment(n)%HI%jed)
          htot = h(i+1,j,1)
          do k=2,nz ; htot = htot + h(i+1,j,k) ; enddo
          Ihtot = G%mask2dCu(I,j) / (htot + h_neglect)
          do k=1,nz ; CS%frhatu(I,j,k) = h(i+1,j,k) * Ihtot ; enddo
        enddo
      endif
    else
      call MOM_error(fatal, "btcalc encountered and OBC segment of indeterminate direction.")
    endif
  enddo ; endif

  if (CS%debug) then
    call uvchksum("btcalc frhat[uv]", CS%frhatu, CS%frhatv, G%HI, &
                  haloshift=0, symmetric=.true., omit_corners=.true., &
                  scalar_pair=.true.)
    if (present(h_u) .and. present(h_v)) &
      call uvchksum("btcalc h_[uv]", h_u, h_v, G%HI, haloshift=0, &
                    symmetric=.true., omit_corners=.true., scale=GV%H_to_m, &
                    scalar_pair=.true.)
    call hchksum(h, "btcalc h",G%HI, haloshift=1, scale=GV%H_to_m)
  endif

end subroutine btcalc

!> The function find_uhbt determines the zonal transport for a given velocity, or with
!! INTEGRAL_BT_CONT=True it determines the time-integrated zonal transport for a given
!! time-integrated velocity.
function find_uhbt(u, BTC) result(uhbt)
  real, intent(in) :: u    !< The local zonal velocity [L T-1 ~> m s-1] or time integrated velocity [L ~> m]
  type(local_BT_cont_u_type), intent(in) :: BTC !< A structure containing various fields that
                           !! allow the barotropic transports to be calculated consistently
                           !! with the layers' continuity equations.  The dimensions of some
                           !! of the elements in this type vary depending on INTEGRAL_BT_CONT.

  real :: uhbt !< The zonal barotropic transport [L2 H T-1 ~> m3 s-1] or time integrated transport [L2 H ~> m3]

  if (u == 0.0) then
    uhbt = 0.0
  elseif (u < BTC%uBT_EE) then
    uhbt = (u - BTC%uBT_EE) * BTC%FA_u_EE + BTC%uh_EE
  elseif (u < 0.0) then
    uhbt = u * (BTC%FA_u_E0 + BTC%uh_crvE * u**2)
  elseif (u <= BTC%uBT_WW) then
    uhbt = u * (BTC%FA_u_W0 + BTC%uh_crvW * u**2)
  else ! (u > BTC%uBT_WW)
    uhbt = (u - BTC%uBT_WW) * BTC%FA_u_WW + BTC%uh_WW
  endif

end function find_uhbt

!> The function find_duhbt_du determines the marginal zonal face area for a given velocity, or
!! with INTEGRAL_BT_CONT=True for a given time-integrated velocity.
function find_duhbt_du(u, BTC) result(duhbt_du)
  real, intent(in) :: u    !< The local zonal velocity [L T-1 ~> m s-1] or time integrated velocity [L ~> m]
  type(local_BT_cont_u_type), intent(in) :: BTC !< A structure containing various fields that
                           !! allow the barotropic transports to be calculated consistently
                           !! with the layers' continuity equations.  The dimensions of some
                           !! of the elements in this type vary depending on INTEGRAL_BT_CONT.
  real :: duhbt_du !< The zonal barotropic face area [L H ~> m2]

  if (u == 0.0) then
    duhbt_du = 0.5*(BTC%FA_u_E0 + BTC%FA_u_W0)  ! Note the potential discontinuity here.
  elseif (u < BTC%uBT_EE) then
    duhbt_du = BTC%FA_u_EE
  elseif (u < 0.0) then
    duhbt_du = (BTC%FA_u_E0 + 3.0*BTC%uh_crvE * u**2)
  elseif (u <= BTC%uBT_WW) then
    duhbt_du = (BTC%FA_u_W0 + 3.0*BTC%uh_crvW * u**2)
  else ! (u > BTC%uBT_WW)
    duhbt_du = BTC%FA_u_WW
  endif

end function find_duhbt_du

!> This function inverts the transport function to determine the barotopic
!! velocity that is consistent with a given transport, or if INTEGRAL_BT_CONT=True
!! this finds the time-integrated velocity that is consistent with a time-integrated transport.
function uhbt_to_ubt(uhbt, BTC, guess) result(ubt)
  real, intent(in) :: uhbt                      !< The barotropic zonal transport that should be inverted for,
                                                !! [H L2 T-1 ~> m3 s-1 or kg s-1] or the time-integrated
                                                !! transport [H L2 ~> m3 or kg].
  type(local_BT_cont_u_type), intent(in) :: BTC !< A structure containing various fields that allow the
                                                !! barotropic transports to be calculated consistently with the
                                                !! layers' continuity equations.  The dimensions of some
                                                !! of the elements in this type vary depending on INTEGRAL_BT_CONT.
  real, optional, intent(in) :: guess           !< A guess at what ubt will be [L T-1 ~> m s-1] or [L ~> m].
                                                !! The result is not allowed to be dramatically larger than guess.
  real :: ubt                                   !< The result - The velocity that gives uhbt transport [L T-1 ~> m s-1]
                                                !! or the time-integrated velocity [L ~> m].

  ! Local variables
  real :: ubt_min, ubt_max       ! Bounding values of vbt [L T-1 ~> m s-1] or [L ~> m]
  real :: uhbt_err               ! The transport error [H L2 T-1 ~> m3 s-1 or kg s-1] or [H L2 ~> m3 or kg].
  real :: derr_du                ! The change in transport error with vbt, i.e. the face area [H L ~> m2 or kg m-1].
  real :: uherr_min, uherr_max   ! The bounding values of the transport error [H L2 T-1 ~> m3 s-1 or kg s-1]
                                 ! or [H L2 ~> m3 or kg].
  real, parameter :: tol = 1.0e-10 ! A fractional match tolerance [nondim]
  real :: dvel  ! Temporary variable used in the limiting the velocity [L T-1 ~> m s-1] or [L ~> m].
  real :: vsr   ! Temporary variable used in the limiting the velocity [nondim].
  real, parameter :: vs1 = 1.25  ! Nondimensional parameters used in limiting
  real, parameter :: vs2 = 2.0   ! the velocity, starting at vs1, with the
                                 ! maximum increase of vs2, both nondim.
  integer :: itt, max_itt = 20

  ! Find the value of ubt that gives uhbt.
  if (uhbt == 0.0) then
    ubt = 0.0
  elseif (uhbt < BTC%uh_EE) then
    ubt = BTC%uBT_EE + (uhbt - BTC%uh_EE) / BTC%FA_u_EE
  elseif (uhbt < 0.0) then
    ! Iterate to convergence with Newton's method (when bounded) and the
    ! false position method otherwise.  ubt will be negative.
    ubt_min = BTC%uBT_EE ; uherr_min = BTC%uh_EE - uhbt
    ubt_max = 0.0 ; uherr_max = -uhbt
    ! Use a false-position method first guess.
    ubt = BTC%uBT_EE * (uhbt / BTC%uh_EE)
    do itt = 1, max_itt
      uhbt_err = ubt * (BTC%FA_u_E0 + BTC%uh_crvE * ubt**2) - uhbt

      if (abs(uhbt_err) < tol*abs(uhbt)) exit
      if (uhbt_err > 0.0) then ; ubt_max = ubt ; uherr_max = uhbt_err ; endif
      if (uhbt_err < 0.0) then ; ubt_min = ubt ; uherr_min = uhbt_err ; endif

      derr_du = BTC%FA_u_E0 + 3.0 * BTC%uh_crvE * ubt**2
      if ((uhbt_err >= derr_du*(ubt - ubt_min)) .or. &
          (-uhbt_err >= derr_du*(ubt_max - ubt)) .or. (derr_du <= 0.0)) then
        ! Use a false-position method guess.
        ubt = ubt_max + (ubt_min-ubt_max) * (uherr_max / (uherr_max-uherr_min))
      else ! Use Newton's method.
        ubt = ubt - uhbt_err / derr_du
        if (abs(uhbt_err) < (0.01*tol)*abs(ubt_min*derr_du)) exit
      endif
    enddo
  elseif (uhbt <= BTC%uh_WW) then
    ! Iterate to convergence with Newton's method.  ubt will be positive.
    ubt_min = 0.0 ; uherr_min = -uhbt
    ubt_max = BTC%uBT_WW ; uherr_max = BTC%uh_WW - uhbt
    ! Use a false-position method first guess.
    ubt = BTC%uBT_WW * (uhbt / BTC%uh_WW)
    do itt = 1, max_itt
      uhbt_err = ubt * (BTC%FA_u_W0 + BTC%uh_crvW * ubt**2) - uhbt

      if (abs(uhbt_err) < tol*abs(uhbt)) exit
      if (uhbt_err > 0.0) then ; ubt_max = ubt ; uherr_max = uhbt_err ; endif
      if (uhbt_err < 0.0) then ; ubt_min = ubt ; uherr_min = uhbt_err ; endif

      derr_du = BTC%FA_u_W0 + 3.0 * BTC%uh_crvW * ubt**2
      if ((uhbt_err >= derr_du*(ubt - ubt_min)) .or. &
          (-uhbt_err >= derr_du*(ubt_max - ubt)) .or. (derr_du <= 0.0)) then
        ! Use a false-position method guess.
        ubt = ubt_min + (ubt_max-ubt_min) * (-uherr_min / (uherr_max-uherr_min))
      else ! Use Newton's method.
        ubt = ubt - uhbt_err / derr_du
        if (abs(uhbt_err) < (0.01*tol)*(ubt_max*derr_du)) exit
      endif
    enddo
  else ! (uhbt > BTC%uh_WW)
    ubt = BTC%uBT_WW + (uhbt - BTC%uh_WW) / BTC%FA_u_WW
  endif

  if (present(guess)) then
    dvel = abs(ubt) - vs1*abs(guess)
    if (dvel > 0.0) then ! Limit the velocity
      if (dvel < 40.0 * (abs(guess)*(vs2-vs1)) ) then
        vsr = vs2 - (vs2-vs1) * exp(-dvel / (abs(guess)*(vs2-vs1)))
      else  ! The exp is less than 4e-18 anyway in this case, so neglect it.
        vsr = vs2
      endif
      ubt = SIGN(vsr * guess, ubt)
    endif
  endif

end function uhbt_to_ubt

!> The function find_vhbt determines the meridional transport for a given velocity, or with
!! INTEGRAL_BT_CONT=True it determines the time-integrated meridional transport for a given
!! time-integrated velocity.
function find_vhbt(v, BTC) result(vhbt)
  real, intent(in) :: v    !< The local meridional velocity [L T-1 ~> m s-1] or time integrated velocity [L ~> m]
  type(local_BT_cont_v_type), intent(in) :: BTC !< A structure containing various fields that
                           !! allow the barotropic transports to be calculated consistently
                           !! with the layers' continuity equations.  The dimensions of some
                           !! of the elements in this type vary depending on INTEGRAL_BT_CONT.
  real :: vhbt !< The meridional barotropic transport [L2 H T-1 ~> m3 s-1] or time integrated transport [L2 H ~> m3]

  if (v == 0.0) then
    vhbt = 0.0
  elseif (v < BTC%vBT_NN) then
    vhbt = (v - BTC%vBT_NN) * BTC%FA_v_NN + BTC%vh_NN
  elseif (v < 0.0) then
    vhbt = v * (BTC%FA_v_N0 + BTC%vh_crvN * v**2)
  elseif (v <= BTC%vBT_SS) then
    vhbt = v * (BTC%FA_v_S0 + BTC%vh_crvS * v**2)
  else ! (v > BTC%vBT_SS)
    vhbt = (v - BTC%vBT_SS) * BTC%FA_v_SS + BTC%vh_SS
  endif

end function find_vhbt

!> The function find_dvhbt_dv determines the marginal meridional face area for a given velocity, or
!! with INTEGRAL_BT_CONT=True for a given time-integrated velocity.
function find_dvhbt_dv(v, BTC) result(dvhbt_dv)
  real, intent(in) :: v    !< The local meridional velocity [L T-1 ~> m s-1] or time integrated velocity [L ~> m]
  type(local_BT_cont_v_type), intent(in) :: BTC !< A structure containing various fields that
                           !! allow the barotropic transports to be calculated consistently
                           !! with the layers' continuity equations.  The dimensions of some
                           !! of the elements in this type vary depending on INTEGRAL_BT_CONT.
  real :: dvhbt_dv !< The meridional barotropic face area [L H ~> m2]

  if (v == 0.0) then
    dvhbt_dv = 0.5*(BTC%FA_v_N0 + BTC%FA_v_S0)  ! Note the potential discontinuity here.
  elseif (v < BTC%vBT_NN) then
    dvhbt_dv = BTC%FA_v_NN
  elseif (v < 0.0) then
    dvhbt_dv = BTC%FA_v_N0 + 3.0*BTC%vh_crvN * v**2
  elseif (v <= BTC%vBT_SS) then
    dvhbt_dv = BTC%FA_v_S0 + 3.0*BTC%vh_crvS * v**2
  else ! (v > BTC%vBT_SS)
    dvhbt_dv = BTC%FA_v_SS
  endif

end function find_dvhbt_dv

!> This function inverts the transport function to determine the barotopic
!! velocity that is consistent with a given transport, or if INTEGRAL_BT_CONT=True
!! this finds the time-integrated velocity that is consistent with a time-integrated transport.
function vhbt_to_vbt(vhbt, BTC, guess) result(vbt)
  real, intent(in) :: vhbt                      !< The barotropic meridional transport that should be
                                                !! inverted for [H L2 T-1 ~> m3 s-1 or kg s-1] or the
                                                !! time-integrated transport [H L2 ~> m3 or kg].
  type(local_BT_cont_v_type), intent(in) :: BTC !< A structure containing various fields that allow the
                                                !! barotropic transports to be calculated consistently
                                                !! with the layers' continuity equations.  The dimensions of some
                                                !! of the elements in this type vary depending on INTEGRAL_BT_CONT.
  real, optional, intent(in) :: guess           !< A guess at what vbt will be [L T-1 ~> m s-1] or [L ~> m].
                                                !! The result is not allowed to be dramatically larger than guess.
  real :: vbt                                   !< The result - The velocity that gives vhbt transport [L T-1 ~> m s-1]
                                                !! or the time-integrated velocity [L ~> m].

  ! Local variables
  real :: vbt_min, vbt_max       ! Bounding values of vbt [L T-1 ~> m s-1] or [L ~> m]
  real :: vhbt_err               ! The transport error [H L2 T-1 ~> m3 s-1 or kg s-1] or [H L2 ~> m3 or kg].
  real :: derr_dv                ! The change in transport error with vbt, i.e. the face area [H L ~> m2 or kg m-1].
  real :: vherr_min, vherr_max   ! The bounding values of the transport error [H L2 T-1 ~> m3 s-1 or kg s-1]
                                 ! or [H L2 ~> m3 or kg].
  real, parameter :: tol = 1.0e-10 ! A fractional match tolerance [nondim]
  real :: dvel  ! Temporary variable used in the limiting the velocity [L T-1 ~> m s-1] or [L ~> m].
  real :: vsr   ! Temporary variable used in the limiting the velocity [nondim].
  real, parameter :: vs1 = 1.25  ! Nondimensional parameters used in limiting
  real, parameter :: vs2 = 2.0   ! the velocity, starting at vs1, with the
                                 ! maximum increase of vs2, both nondim.
  integer :: itt, max_itt = 20

  ! Find the value of vbt that gives vhbt.
  if (vhbt == 0.0) then
    vbt = 0.0
  elseif (vhbt < BTC%vh_NN) then
    vbt = BTC%vBT_NN + (vhbt - BTC%vh_NN) / BTC%FA_v_NN
  elseif (vhbt < 0.0) then
    ! Iterate to convergence with Newton's method (when bounded) and the
    ! false position method otherwise.  vbt will be negative.
    vbt_min = BTC%vBT_NN ; vherr_min = BTC%vh_NN - vhbt
    vbt_max = 0.0 ; vherr_max = -vhbt
    ! Use a false-position method first guess.
    vbt = BTC%vBT_NN * (vhbt / BTC%vh_NN)
    do itt = 1, max_itt
      vhbt_err = vbt * (BTC%FA_v_N0 + BTC%vh_crvN * vbt**2) - vhbt

      if (abs(vhbt_err) < tol*abs(vhbt)) exit
      if (vhbt_err > 0.0) then ; vbt_max = vbt ; vherr_max = vhbt_err ; endif
      if (vhbt_err < 0.0) then ; vbt_min = vbt ; vherr_min = vhbt_err ; endif

      derr_dv = BTC%FA_v_N0 + 3.0 * BTC%vh_crvN * vbt**2
      if ((vhbt_err >= derr_dv*(vbt - vbt_min)) .or. &
          (-vhbt_err >= derr_dv*(vbt_max - vbt)) .or. (derr_dv <= 0.0)) then
        ! Use a false-position method guess.
        vbt = vbt_max + (vbt_min-vbt_max) * (vherr_max / (vherr_max-vherr_min))
      else ! Use Newton's method.
        vbt = vbt - vhbt_err / derr_dv
        if (abs(vhbt_err) < (0.01*tol)*abs(derr_dv*vbt_min)) exit
      endif
    enddo
  elseif (vhbt <= BTC%vh_SS) then
    ! Iterate to convergence with Newton's method.  vbt will be positive.
    vbt_min = 0.0 ; vherr_min = -vhbt
    vbt_max = BTC%vBT_SS ; vherr_max = BTC%vh_SS - vhbt
    ! Use a false-position method first guess.
    vbt = BTC%vBT_SS * (vhbt / BTC%vh_SS)
    do itt = 1, max_itt
      vhbt_err = vbt * (BTC%FA_v_S0 + BTC%vh_crvS * vbt**2) - vhbt

      if (abs(vhbt_err) < tol*abs(vhbt)) exit
      if (vhbt_err > 0.0) then ; vbt_max = vbt ; vherr_max = vhbt_err ; endif
      if (vhbt_err < 0.0) then ; vbt_min = vbt ; vherr_min = vhbt_err ; endif

      derr_dv = BTC%FA_v_S0 + 3.0 * BTC%vh_crvS * vbt**2
      if ((vhbt_err >= derr_dv*(vbt - vbt_min)) .or. &
          (-vhbt_err >= derr_dv*(vbt_max - vbt)) .or. (derr_dv <= 0.0)) then
        ! Use a false-position method guess.
        vbt = vbt_min + (vbt_max-vbt_min) * (-vherr_min / (vherr_max-vherr_min))
      else ! Use Newton's method.
        vbt = vbt - vhbt_err / derr_dv
        if (abs(vhbt_err) < (0.01*tol)*(vbt_max*derr_dv)) exit
      endif
    enddo
  else ! (vhbt > BTC%vh_SS)
    vbt = BTC%vBT_SS + (vhbt - BTC%vh_SS) / BTC%FA_v_SS
  endif

  if (present(guess)) then
    dvel = abs(vbt) - vs1*abs(guess)
    if (dvel > 0.0) then ! Limit the velocity
      if (dvel < 40.0 * (abs(guess)*(vs2-vs1)) ) then
        vsr = vs2 - (vs2-vs1) * exp(-dvel / (abs(guess)*(vs2-vs1)))
      else  ! The exp is less than 4e-18 anyway in this case, so neglect it.
        vsr = vs2
      endif
      vbt = SIGN(guess * vsr, vbt)
    endif
  endif

end function vhbt_to_vbt

!> This subroutine sets up reordered versions of the BT_cont type in the
!! local_BT_cont types, which have wide halos properly filled in.
subroutine set_local_BT_cont_types(BT_cont, BTCL_u, BTCL_v, G, US, MS, BT_Domain, halo, dt_baroclinic)
  type(BT_cont_type),                                    intent(inout) :: BT_cont    !< The BT_cont_type input to the
                                                                                     !! barotropic solver.
  type(memory_size_type),                                intent(in)    :: MS         !< A type that describes the
                                                                                     !! memory sizes of the argument
                                                                                     !! arrays.
  type(local_BT_cont_u_type), dimension(SZIBW_(MS),SZJW_(MS)), intent(out) :: BTCL_u !< A structure with the u
                                                                                     !! information from BT_cont.
  type(local_BT_cont_v_type), dimension(SZIW_(MS),SZJBW_(MS)), intent(out) :: BTCL_v !< A structure with the v
                                                                                     !! information from BT_cont.
  type(ocean_grid_type),                                 intent(in)    :: G          !< The ocean's grid structure.
  type(unit_scale_type),                                 intent(in)    :: US         !< A dimensional unit scaling type
  type(MOM_domain_type),                                 intent(inout) :: BT_Domain  !< The domain to use for updating
                                                                                     !! the halos of wide arrays.
  integer,                                     optional, intent(in)    :: halo       !< The extra halo size to use here.
  real,                                        optional, intent(in)    :: dt_baroclinic !< The baroclinic time step
                                                                                     !! [T ~> s], which is provided if
                                                                                     !! INTEGRAL_BT_CONTINUITY is true.

  ! Local variables
  real, dimension(SZIBW_(MS),SZJW_(MS)) :: &
    u_polarity, &      ! An array used to test for halo update polarity [nondim]
    uBT_EE, uBT_WW, &  ! Zonal velocities at which the form of the fit changes [L T-1 ~> m s-1]
    FA_u_EE, FA_u_E0, FA_u_W0, FA_u_WW ! Zonal face areas [H L ~> m2 or kg m-1]
  real, dimension(SZIW_(MS),SZJBW_(MS)) :: &
    v_polarity, &      ! An array used to test for halo update polarity [nondim]
    vBT_NN, vBT_SS, &  ! Meridional velocities at which the form of the fit changes [L T-1 ~> m s-1]
    FA_v_NN, FA_v_N0, FA_v_S0, FA_v_SS ! Meridional face areas [H L ~> m2 or kg m-1]
  real :: dt ! The baroclinic timestep [T ~> s] or 1.0 [nondim]
  real, parameter :: C1_3 = 1.0/3.0
  integer :: i, j, is, ie, js, je, hs

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  hs = 1 ; if (present(halo)) hs = max(halo,0)
  dt = 1.0 ; if (present(dt_baroclinic)) dt = dt_baroclinic

  ! Copy the BT_cont arrays into symmetric, potentially wide haloed arrays.
!$OMP parallel default(none) shared(is,ie,js,je,hs,u_polarity,uBT_EE,uBT_WW,FA_u_EE, &
!$OMP                               FA_u_E0,FA_u_W0,FA_u_WW,v_polarity,vBT_NN,vBT_SS,&
!$OMP                               FA_v_NN,FA_v_N0,FA_v_S0,FA_v_SS,BT_cont )
!$OMP do
  do j=js-hs,je+hs ; do i=is-hs-1,ie+hs
    u_polarity(i,j) = 1.0
    uBT_EE(i,j) = 0.0 ; uBT_WW(i,j) = 0.0
    FA_u_EE(i,j) = 0.0 ; FA_u_E0(i,j) = 0.0 ; FA_u_W0(i,j) = 0.0 ; FA_u_WW(i,j) = 0.0
  enddo ; enddo
!$OMP do
  do j=js-hs-1,je+hs ; do i=is-hs,ie+hs
    v_polarity(i,j) = 1.0
    vBT_NN(i,j) = 0.0 ; vBT_SS(i,j) = 0.0
    FA_v_NN(i,j) = 0.0 ; FA_v_N0(i,j) = 0.0 ; FA_v_S0(i,j) = 0.0 ; FA_v_SS(i,j) = 0.0
  enddo ; enddo
!$OMP do
  do j=js,je ; do I=is-1,ie
    uBT_EE(I,j) = BT_cont%uBT_EE(I,j) ; uBT_WW(I,j) = BT_cont%uBT_WW(I,j)
    FA_u_EE(I,j) = BT_cont%FA_u_EE(I,j) ; FA_u_E0(I,j) = BT_cont%FA_u_E0(I,j)
    FA_u_W0(I,j) = BT_cont%FA_u_W0(I,j) ; FA_u_WW(I,j) = BT_cont%FA_u_WW(I,j)
  enddo ; enddo
!$OMP do
  do J=js-1,je ; do i=is,ie
    vBT_NN(i,J) = BT_cont%vBT_NN(i,J) ; vBT_SS(i,J) = BT_cont%vBT_SS(i,J)
    FA_v_NN(i,J) = BT_cont%FA_v_NN(i,J) ; FA_v_N0(i,J) = BT_cont%FA_v_N0(i,J)
    FA_v_S0(i,J) = BT_cont%FA_v_S0(i,J) ; FA_v_SS(i,J) = BT_cont%FA_v_SS(i,J)
  enddo ; enddo
!$OMP end parallel

  if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
  if (id_clock_pass_pre > 0) call cpu_clock_begin(id_clock_pass_pre)
!--- begin setup for group halo update
  call create_group_pass(BT_cont%pass_polarity_BT, u_polarity, v_polarity, BT_Domain)
  call create_group_pass(BT_cont%pass_polarity_BT, uBT_EE, vBT_NN, BT_Domain)
  call create_group_pass(BT_cont%pass_polarity_BT, uBT_WW, vBT_SS, BT_Domain)

  call create_group_pass(BT_cont%pass_FA_uv, FA_u_EE, FA_v_NN, BT_Domain, To_All+Scalar_Pair)
  call create_group_pass(BT_cont%pass_FA_uv, FA_u_E0, FA_v_N0, BT_Domain, To_All+Scalar_Pair)
  call create_group_pass(BT_cont%pass_FA_uv, FA_u_W0, FA_v_S0, BT_Domain, To_All+Scalar_Pair)
  call create_group_pass(BT_cont%pass_FA_uv, FA_u_WW, FA_v_SS, BT_Domain, To_All+Scalar_Pair)
!--- end setup for group halo update
  ! Do halo updates on BT_cont.
  call do_group_pass(BT_cont%pass_polarity_BT, BT_Domain)
  call do_group_pass(BT_cont%pass_FA_uv, BT_Domain)
  if (id_clock_pass_pre > 0) call cpu_clock_end(id_clock_pass_pre)
  if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)

  !$OMP parallel default(shared)
  !$OMP do
  do j=js-hs,je+hs ; do I=is-hs-1,ie+hs
    BTCL_u(I,j)%FA_u_EE = FA_u_EE(I,j) ; BTCL_u(I,j)%FA_u_E0 = FA_u_E0(I,j)
    BTCL_u(I,j)%FA_u_W0 = FA_u_W0(I,j) ; BTCL_u(I,j)%FA_u_WW = FA_u_WW(I,j)
    BTCL_u(I,j)%uBT_EE = dt*uBT_EE(I,j)   ; BTCL_u(I,j)%uBT_WW = dt*uBT_WW(I,j)
    ! Check for reversed polarity in the tripolar halo regions.
    if (u_polarity(I,j) < 0.0) then
      call swap(BTCL_u(I,j)%FA_u_EE, BTCL_u(I,j)%FA_u_WW)
      call swap(BTCL_u(I,j)%FA_u_E0, BTCL_u(I,j)%FA_u_W0)
      call swap(BTCL_u(I,j)%uBT_EE,  BTCL_u(I,j)%uBT_WW)
    endif

    BTCL_u(I,j)%uh_EE = BTCL_u(I,j)%uBT_EE * &
        (C1_3 * (2.0*BTCL_u(I,j)%FA_u_E0 + BTCL_u(I,j)%FA_u_EE))
    BTCL_u(I,j)%uh_WW = BTCL_u(I,j)%uBT_WW * &
        (C1_3 * (2.0*BTCL_u(I,j)%FA_u_W0 + BTCL_u(I,j)%FA_u_WW))

    BTCL_u(I,j)%uh_crvE = 0.0 ; BTCL_u(I,j)%uh_crvW = 0.0
    if (abs(BTCL_u(I,j)%uBT_WW) > 0.0) BTCL_u(I,j)%uh_crvW = &
      (C1_3 * (BTCL_u(I,j)%FA_u_WW - BTCL_u(I,j)%FA_u_W0)) / BTCL_u(I,j)%uBT_WW**2
    if (abs(BTCL_u(I,j)%uBT_EE) > 0.0) BTCL_u(I,j)%uh_crvE = &
      (C1_3 * (BTCL_u(I,j)%FA_u_EE - BTCL_u(I,j)%FA_u_E0)) / BTCL_u(I,j)%uBT_EE**2
  enddo ; enddo
  !$OMP do
  do J=js-hs-1,je+hs ; do i=is-hs,ie+hs
    BTCL_v(i,J)%FA_v_NN = FA_v_NN(i,J) ; BTCL_v(i,J)%FA_v_N0 = FA_v_N0(i,J)
    BTCL_v(i,J)%FA_v_S0 = FA_v_S0(i,J) ; BTCL_v(i,J)%FA_v_SS = FA_v_SS(i,J)
    BTCL_v(i,J)%vBT_NN = dt*vBT_NN(i,J)   ; BTCL_v(i,J)%vBT_SS = dt*vBT_SS(i,J)
    ! Check for reversed polarity in the tripolar halo regions.
    if (v_polarity(i,J) < 0.0) then
      call swap(BTCL_v(i,J)%FA_v_NN, BTCL_v(i,J)%FA_v_SS)
      call swap(BTCL_v(i,J)%FA_v_N0, BTCL_v(i,J)%FA_v_S0)
      call swap(BTCL_v(i,J)%vBT_NN,  BTCL_v(i,J)%vBT_SS)
    endif

    BTCL_v(i,J)%vh_NN = BTCL_v(i,J)%vBT_NN * &
        (C1_3 * (2.0*BTCL_v(i,J)%FA_v_N0 + BTCL_v(i,J)%FA_v_NN))
    BTCL_v(i,J)%vh_SS = BTCL_v(i,J)%vBT_SS * &
        (C1_3 * (2.0*BTCL_v(i,J)%FA_v_S0 + BTCL_v(i,J)%FA_v_SS))

    BTCL_v(i,J)%vh_crvN = 0.0 ; BTCL_v(i,J)%vh_crvS = 0.0
    if (abs(BTCL_v(i,J)%vBT_SS) > 0.0) BTCL_v(i,J)%vh_crvS = &
      (C1_3 * (BTCL_v(i,J)%FA_v_SS - BTCL_v(i,J)%FA_v_S0)) / BTCL_v(i,J)%vBT_SS**2
    if (abs(BTCL_v(i,J)%vBT_NN) > 0.0) BTCL_v(i,J)%vh_crvN = &
      (C1_3 * (BTCL_v(i,J)%FA_v_NN - BTCL_v(i,J)%FA_v_N0)) / BTCL_v(i,J)%vBT_NN**2
  enddo ; enddo
  !$OMP end parallel
end subroutine set_local_BT_cont_types


!> Adjust_local_BT_cont_types expands the range of velocities with a cubic curve
!! translating velocities into transports to match the inital values of velocities and
!! summed transports when the velocities are larger than the first guesses of the cubic
!! transition velocities used to set up the local_BT_cont types.
subroutine adjust_local_BT_cont_types(ubt, uhbt, vbt, vhbt, BTCL_u, BTCL_v, &
                                      G, US, MS, halo, dt_baroclinic)
  type(memory_size_type), intent(in)  :: MS   !< A type that describes the memory sizes of the argument arrays.
  real, dimension(SZIBW_(MS),SZJW_(MS)), &
                          intent(in)  :: ubt  !< The linearization zonal barotropic velocity [L T-1 ~> m s-1].
  real, dimension(SZIBW_(MS),SZJW_(MS)), &
                          intent(in)  :: uhbt !< The linearization zonal barotropic transport
                                              !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), &
                          intent(in)  :: vbt  !< The linearization meridional barotropic velocity [L T-1 ~> m s-1].
  real, dimension(SZIW_(MS),SZJBW_(MS)), &
                          intent(in)  :: vhbt !< The linearization meridional barotropic transport
                                              !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(local_BT_cont_u_type), dimension(SZIBW_(MS),SZJW_(MS)), &
                          intent(out) :: BTCL_u !< A structure with the u information from BT_cont.
  type(local_BT_cont_v_type), dimension(SZIW_(MS),SZJBW_(MS)), &
                          intent(out) :: BTCL_v !< A structure with the v information from BT_cont.
  type(ocean_grid_type),  intent(in)  :: G    !< The ocean's grid structure.
  type(unit_scale_type),  intent(in)  :: US   !< A dimensional unit scaling type
  integer,      optional, intent(in)  :: halo !< The extra halo size to use here.
  real,         optional, intent(in)  :: dt_baroclinic !< The baroclinic time step [T ~> s], which is
                                                       !! provided if INTEGRAL_BT_CONTINUITY is true.

  ! Local variables
  real, dimension(SZIBW_(MS),SZJW_(MS)) :: &
    u_polarity, uBT_EE, uBT_WW, FA_u_EE, FA_u_E0, FA_u_W0, FA_u_WW
  real, dimension(SZIW_(MS),SZJBW_(MS)) :: &
    v_polarity, vBT_NN, vBT_SS, FA_v_NN, FA_v_N0, FA_v_S0, FA_v_SS
  real :: dt ! The baroclinic timestep [T ~> s] or 1.0 [nondim]
  real, parameter :: C1_3 = 1.0/3.0
  integer :: i, j, is, ie, js, je, hs

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  hs = 1 ; if (present(halo)) hs = max(halo,0)
  dt = 1.0 ; if (present(dt_baroclinic)) dt = dt_baroclinic

  !$OMP parallel do default(shared)
  do j=js-hs,je+hs ; do I=is-hs-1,ie+hs
    if ((dt*ubt(I,j) > BTCL_u(I,j)%uBT_WW) .and. (dt*uhbt(I,j) > BTCL_u(I,j)%uh_WW)) then
      ! Expand the cubic fit to use this new point.  ubt is negative.
      BTCL_u(I,j)%ubt_WW = dt * ubt(I,j)
      if (3.0*uhbt(I,j) < 2.0*ubt(I,j) * BTCL_u(I,j)%FA_u_W0) then
        ! No further bounding is needed.
        BTCL_u(I,j)%uh_crvW = (uhbt(I,j) - ubt(I,j) * BTCL_u(I,j)%FA_u_W0) / (dt**2 * ubt(I,j)**3)
      else ! This should not happen often!
        BTCL_u(I,j)%FA_u_W0 = 1.5*uhbt(I,j) / ubt(I,j)
        BTCL_u(I,j)%uh_crvW = -0.5*uhbt(I,j) / (dt**2 * ubt(I,j)**3)
      endif
      BTCL_u(I,j)%uh_WW = dt * uhbt(I,j)
      ! I don't know whether this is helpful.
!     BTCL_u(I,j)%FA_u_WW = min(BTCL_u(I,j)%FA_u_WW, uhbt(I,j) / ubt(I,j))
    elseif ((dt*ubt(I,j) < BTCL_u(I,j)%uBT_EE) .and. (dt*uhbt(I,j) < BTCL_u(I,j)%uh_EE)) then
      ! Expand the cubic fit to use this new point.  ubt is negative.
      BTCL_u(I,j)%ubt_EE = dt * ubt(I,j)
      if (3.0*uhbt(I,j) < 2.0*ubt(I,j) * BTCL_u(I,j)%FA_u_E0) then
        ! No further bounding is needed.
        BTCL_u(I,j)%uh_crvE = (uhbt(I,j) - ubt(I,j) * BTCL_u(I,j)%FA_u_E0) / (dt**2 * ubt(I,j)**3)
      else ! This should not happen often!
        BTCL_u(I,j)%FA_u_E0 = 1.5*uhbt(I,j) / ubt(I,j)
        BTCL_u(I,j)%uh_crvE = -0.5*uhbt(I,j) / (dt**2 * ubt(I,j)**3)
      endif
      BTCL_u(I,j)%uh_EE = dt * uhbt(I,j)
      ! I don't know whether this is helpful.
!     BTCL_u(I,j)%FA_u_EE = min(BTCL_u(I,j)%FA_u_EE, uhbt(I,j) / ubt(I,j))
    endif
  enddo ; enddo
  !$OMP parallel do default(shared)
  do J=js-hs-1,je+hs ; do i=is-hs,ie+hs
    if ((dt*vbt(i,J) > BTCL_v(i,J)%vBT_SS) .and. (dt*vhbt(i,J) > BTCL_v(i,J)%vh_SS)) then
      ! Expand the cubic fit to use this new point.  vbt is negative.
      BTCL_v(i,J)%vbt_SS = dt * vbt(i,J)
      if (3.0*vhbt(i,J) < 2.0*vbt(i,J) * BTCL_v(i,J)%FA_v_S0) then
        ! No further bounding is needed.
        BTCL_v(i,J)%vh_crvS = (vhbt(i,J) - vbt(i,J) * BTCL_v(i,J)%FA_v_S0) /  (dt**2 * vbt(i,J)**3)
      else ! This should not happen often!
        BTCL_v(i,J)%FA_v_S0 = 1.5*vhbt(i,J) / (vbt(i,J))
        BTCL_v(i,J)%vh_crvS = -0.5*vhbt(i,J) /  (dt**2 * vbt(i,J)**3)
      endif
      BTCL_v(i,J)%vh_SS = dt * vhbt(i,J)
      ! I don't know whether this is helpful.
!     BTCL_v(i,J)%FA_v_SS = min(BTCL_v(i,J)%FA_v_SS, vhbt(i,J) / vbt(i,J))
    elseif ((dt*vbt(i,J) < BTCL_v(i,J)%vBT_NN) .and. (dt*vhbt(i,J) < BTCL_v(i,J)%vh_NN)) then
      ! Expand the cubic fit to use this new point.  vbt is negative.
      BTCL_v(i,J)%vbt_NN = dt * vbt(i,J)
      if (3.0*vhbt(i,J) < 2.0*vbt(i,J) * BTCL_v(i,J)%FA_v_N0) then
        ! No further bounding is needed.
        BTCL_v(i,J)%vh_crvN = (vhbt(i,J) - vbt(i,J) * BTCL_v(i,J)%FA_v_N0) /  (dt**2 * vbt(i,J)**3)
      else ! This should not happen often!
        BTCL_v(i,J)%FA_v_N0 = 1.5*vhbt(i,J) / (vbt(i,J))
        BTCL_v(i,J)%vh_crvN = -0.5*vhbt(i,J) /  (dt**2 * vbt(i,J)**3)
      endif
      BTCL_v(i,J)%vh_NN = dt * vhbt(i,J)
      ! I don't know whether this is helpful.
!     BTCL_v(i,J)%FA_v_NN = min(BTCL_v(i,J)%FA_v_NN, vhbt(i,J) / vbt(i,J))
    endif
  enddo ; enddo

end subroutine adjust_local_BT_cont_types

!> This subroutine uses the BTCL types to find typical or maximum face
!! areas, which can then be used for finding wave speeds, etc.
subroutine BT_cont_to_face_areas(BT_cont, Datu, Datv, G, US, MS, halo, maximize)
  type(BT_cont_type),     intent(inout) :: BT_cont    !< The BT_cont_type input to the
                                                      !! barotropic solver.
  type(memory_size_type), intent(in)    :: MS         !< A type that describes the memory
                                                      !! sizes of the argument arrays.
  real, dimension(MS%isdw-1:MS%iedw,MS%jsdw:MS%jedw), &
                          intent(out)   :: Datu       !< The effective zonal face area [H L ~> m2 or kg m-1].
  real, dimension(MS%isdw:MS%iedw,MS%jsdw-1:MS%jedw), &
                          intent(out)   :: Datv       !< The effective meridional face area [H L ~> m2 or kg m-1].
  type(ocean_grid_type),  intent(in)    :: G          !< The ocean's grid structure.
  type(unit_scale_type),  intent(in)    :: US   !< A dimensional unit scaling type
  integer,      optional, intent(in)    :: halo       !< The extra halo size to use here.
  logical,      optional, intent(in)    :: maximize   !< If present and true, find the
                                                      !! maximum face area for any velocity.

  ! Local variables
  logical :: find_max
  integer :: i, j, is, ie, js, je, hs
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  hs = 1 ; if (present(halo)) hs = max(halo,0)
  find_max = .false. ; if (present(maximize)) find_max = maximize

  if (find_max) then
    do j=js-hs,je+hs ; do I=is-1-hs,ie+hs
      Datu(I,j) = max(BT_cont%FA_u_EE(I,j), BT_cont%FA_u_E0(I,j), &
                      BT_cont%FA_u_W0(I,j), BT_cont%FA_u_WW(I,j))
    enddo ; enddo
    do J=js-1-hs,je+hs ; do i=is-hs,ie+hs
      Datv(i,J) = max(BT_cont%FA_v_NN(i,J), BT_cont%FA_v_N0(i,J), &
                      BT_cont%FA_v_S0(i,J), BT_cont%FA_v_SS(i,J))
    enddo ; enddo
  else
    do j=js-hs,je+hs ; do I=is-1-hs,ie+hs
      Datu(I,j) = 0.5 * (BT_cont%FA_u_E0(I,j) + BT_cont%FA_u_W0(I,j))
    enddo ; enddo
    do J=js-1-hs,je+hs ; do i=is-hs,ie+hs
      Datv(i,J) = 0.5 * (BT_cont%FA_v_N0(i,J) + BT_cont%FA_v_S0(i,J))
    enddo ; enddo
  endif

end subroutine BT_cont_to_face_areas

!> Swap the values of two real variables
subroutine swap(a,b)
  real, intent(inout) :: a !< The first variable to be swapped.
  real, intent(inout) :: b !< The second variable to be swapped.
  real :: tmp
  tmp = a ; a = b ; b = tmp
end subroutine swap

!> This subroutine determines the open face areas of cells for calculating
!! the barotropic transport.
subroutine find_face_areas(Datu, Datv, G, GV, US, CS, MS, eta, halo, add_max)
  type(memory_size_type),  intent(in) :: MS    !< A type that describes the memory sizes of the argument arrays.
  real, dimension(MS%isdw-1:MS%iedw,MS%jsdw:MS%jedw), &
                           intent(out) :: Datu !< The open zonal face area [H L ~> m2 or kg m-1].
  real, dimension(MS%isdw:MS%iedw,MS%jsdw-1:MS%jedw), &
                           intent(out) :: Datv !< The open meridional face area [H L ~> m2 or kg m-1].
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  type(barotropic_CS),     pointer     :: CS   !< The control structure returned by a previous
                                               !! call to barotropic_init.
  real, dimension(MS%isdw:MS%iedw,MS%jsdw:MS%jedw), &
                 optional, intent(in)  :: eta  !< The barotropic free surface height anomaly
                                               !! or column mass anomaly [H ~> m or kg m-2].
  integer,       optional, intent(in)  :: halo !< The halo size to use, default = 1.
  real,          optional, intent(in)  :: add_max !< A value to add to the maximum depth (used
                                               !! to overestimate the external wave speed) [Z ~> m].

  ! Local variables
  real :: H1, H2      ! Temporary total thicknesses [H ~> m or kg m-2].
  integer :: i, j, is, ie, js, je, hs
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  hs = 1 ; if (present(halo)) hs = max(halo,0)

!$OMP parallel default(none) shared(is,ie,js,je,hs,eta,GV,CS,Datu,Datv,add_max) &
!$OMP                       private(H1,H2)
  if (present(eta)) then
    ! The use of harmonic mean thicknesses ensure positive definiteness.
    if (GV%Boussinesq) then
!$OMP do
      do j=js-hs,je+hs ; do I=is-1-hs,ie+hs
        H1 = CS%bathyT(i,j)*GV%Z_to_H + eta(i,j) ; H2 = CS%bathyT(i+1,j)*GV%Z_to_H + eta(i+1,j)
        Datu(I,j) = 0.0 ; if ((H1 > 0.0) .and. (H2 > 0.0)) &
        Datu(I,j) = CS%dy_Cu(I,j) * (2.0 * H1 * H2) / (H1 + H2)
!       Datu(I,j) = CS%dy_Cu(I,j) * 0.5 * (H1 + H2)
      enddo ; enddo
!$OMP do
      do J=js-1-hs,je+hs ; do i=is-hs,ie+hs
        H1 = CS%bathyT(i,j)*GV%Z_to_H + eta(i,j) ; H2 = CS%bathyT(i,j+1)*GV%Z_to_H + eta(i,j+1)
        Datv(i,J) = 0.0 ; if ((H1 > 0.0) .and. (H2 > 0.0)) &
        Datv(i,J) = CS%dx_Cv(i,J) * (2.0 * H1 * H2) / (H1 + H2)
!       Datv(i,J) = CS%dy_v(i,J) * 0.5 * (H1 + H2)
      enddo ; enddo
    else
!$OMP do
      do j=js-hs,je+hs ; do I=is-1-hs,ie+hs
        Datu(I,j) = 0.0 ; if ((eta(i,j) > 0.0) .and. (eta(i+1,j) > 0.0)) &
        Datu(I,j) = CS%dy_Cu(I,j) * (2.0 * eta(i,j) * eta(i+1,j)) / &
                                  (eta(i,j) + eta(i+1,j))
        ! Datu(I,j) = CS%dy_Cu(I,j) * 0.5 * (eta(i,j) + eta(i+1,j))
      enddo ; enddo
!$OMP do
      do J=js-1-hs,je+hs ; do i=is-hs,ie+hs
        Datv(i,J) = 0.0 ; if ((eta(i,j) > 0.0) .and. (eta(i,j+1) > 0.0)) &
        Datv(i,J) = CS%dx_Cv(i,J) * (2.0 * eta(i,j) * eta(i,j+1)) / &
                                  (eta(i,j) + eta(i,j+1))
        ! Datv(i,J) = CS%dy_v(i,J) * 0.5 * (eta(i,j) + eta(i,j+1))
      enddo ; enddo
    endif
  elseif (present(add_max)) then
!$OMP do
    do j=js-hs,je+hs ; do I=is-1-hs,ie+hs
      Datu(I,j) = CS%dy_Cu(I,j) * GV%Z_to_H * &
                 (max(CS%bathyT(i+1,j), CS%bathyT(i,j)) + add_max)
    enddo ; enddo
!$OMP do
    do J=js-1-hs,je+hs ; do i=is-hs,ie+hs
      Datv(i,J) = CS%dx_Cv(i,J) * GV%Z_to_H * &
                 (max(CS%bathyT(i,j+1), CS%bathyT(i,j)) + add_max)
    enddo ; enddo
  else
!$OMP do
    do j=js-hs,je+hs ; do I=is-1-hs,ie+hs
      Datu(I, j) = 0.0
      !Would be "if (G%mask2dCu(I,j)>0.) &" is G was valid on BT domain
      if (CS%bathyT(i+1,j)+CS%bathyT(i,j)>0.) &
        Datu(I,j) = 2.0*CS%dy_Cu(I,j) * GV%Z_to_H * &
                  (CS%bathyT(i+1,j) * CS%bathyT(i,j)) / &
                  (CS%bathyT(i+1,j) + CS%bathyT(i,j))
    enddo ; enddo
!$OMP do
    do J=js-1-hs,je+hs ; do i=is-hs,ie+hs
      Datv(i, J) = 0.0
      !Would be "if (G%mask2dCv(i,J)>0.) &" is G was valid on BT domain
      if (CS%bathyT(i,j+1)+CS%bathyT(i,j)>0.) &
        Datv(i,J) = 2.0*CS%dx_Cv(i,J) * GV%Z_to_H * &
                  (CS%bathyT(i,j+1) * CS%bathyT(i,j)) / &
                  (CS%bathyT(i,j+1) + CS%bathyT(i,j))
    enddo ; enddo
  endif
!$OMP end parallel

end subroutine find_face_areas

!> bt_mass_source determines the appropriately limited mass source for
!! the barotropic solver, along with a corrective fictitious mass source that
!! will drive the barotropic estimate of the free surface height toward the
!! baroclinic estimate.
subroutine bt_mass_source(h, eta, set_cor, G, GV, CS)
  type(ocean_grid_type),              intent(in) :: G        !< The ocean's grid structure.
  type(verticalGrid_type),            intent(in) :: GV       !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h  !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),   intent(in) :: eta      !< The free surface height that is to be
                                                             !! corrected [H ~> m or kg m-2].
  logical,                            intent(in) :: set_cor  !< A flag to indicate whether to set the corrective
                                                             !! fluxes (and update the slowly varying part of eta_cor)
                                                             !! (.true.) or whether to incrementally update the
                                                             !! corrective fluxes.
  type(barotropic_CS),                pointer    :: CS       !< The control structure returned by a previous call
                                                             !! to barotropic_init.

  ! Local variables
  real :: h_tot(SZI_(G))      ! The sum of the layer thicknesses [H ~> m or kg m-2].
  real :: eta_h(SZI_(G))      ! The free surface height determined from
                              ! the sum of the layer thicknesses [H ~> m or kg m-2].
  real :: d_eta               ! The difference between estimates of the total
                              ! thicknesses [H ~> m or kg m-2].
  integer :: is, ie, js, je, nz, i, j, k

  if (.not.associated(CS)) call MOM_error(FATAL, "bt_mass_source: "// &
        "Module MOM_barotropic must be initialized before it is used.")
  if (.not.CS%split) return

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  !$OMP parallel do default(shared) private(eta_h,h_tot,d_eta)
  do j=js,je
    do i=is,ie ; h_tot(i) = h(i,j,1) ; enddo
    if (GV%Boussinesq) then
      do i=is,ie ; eta_h(i) = h(i,j,1) - G%bathyT(i,j)*GV%Z_to_H ; enddo
    else
      do i=is,ie ; eta_h(i) = h(i,j,1) ; enddo
    endif
    do k=2,nz ; do i=is,ie
      eta_h(i) = eta_h(i) + h(i,j,k)
      h_tot(i) = h_tot(i) + h(i,j,k)
    enddo ; enddo

    if (set_cor) then
      do i=is,ie
        d_eta = eta_h(i) - eta(i,j)
        CS%eta_cor(i,j) = d_eta
      enddo
    else
      do i=is,ie
        d_eta = eta_h(i) - eta(i,j)
        CS%eta_cor(i,j) = CS%eta_cor(i,j) + d_eta
      enddo
    endif
  enddo

end subroutine bt_mass_source

!> barotropic_init initializes a number of time-invariant fields used in the
!! barotropic calculation and initializes any barotropic fields that have not
!! already been initialized.
subroutine barotropic_init(u, v, h, eta, Time, G, GV, US, param_file, diag, CS, &
                           restart_CS, calc_dtbt, BT_cont, tides_CSp)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)    :: eta  !< Free surface height or column mass anomaly
                                                 !! [Z ~> m] or [H ~> kg m-2].
  type(time_type), target, intent(in)    :: Time !< The current model time.
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic
                                                 !! output.
  type(barotropic_CS),     pointer       :: CS   !< A pointer to the control structure for this module
                                                 !! that is set in register_barotropic_restarts.
  type(MOM_restart_CS),    pointer       :: restart_CS !< A pointer to the restart control structure.
  logical,                 intent(out)   :: calc_dtbt  !< If true, the barotropic time step must
                                                 !! be recalculated before stepping.
  type(BT_cont_type), optional, &
                           pointer       :: BT_cont    !< A structure with elements that describe the
                                                 !! effective open face areas as a function of
                                                 !! barotropic flow.
  type(tidal_forcing_CS), optional, &
                           pointer       :: tides_CSp  !< A pointer to the control structure of the
                                                 !! tide module.

! This include declares and sets the variable "version".
#include "version_variable.h"
  ! Local variables
  character(len=40)  :: mdl = "MOM_barotropic"  ! This module's name.
  real :: Datu(SZIBS_(G),SZJ_(G))   ! Zonal open face area [H L ~> m2 or kg m-1].
  real :: Datv(SZI_(G),SZJBS_(G))   ! Meridional open face area [H L ~> m2 or kg m-1].
  real :: gtot_estimate ! Summed GV%g_prime [L2 Z-1 T-2 ~> m s-2], to give an upper-bound estimate for pbce.
  real :: SSH_extra     ! An estimate of how much higher SSH might get, for use
                        ! in calculating the safe external wave speed [Z ~> m].
  real :: dtbt_input    ! The input value of DTBT, [nondim] if negative or [s] if positive.
  real :: dtbt_tmp      ! A temporary copy of CS%dtbt read from a restart file [T ~> s]
  real :: wave_drag_scale ! A scaling factor for the barotropic linear wave drag
                          ! piston velocities.
  character(len=200) :: inputdir       ! The directory in which to find input files.
  character(len=200) :: wave_drag_file ! The file from which to read the wave
                                       ! drag piston velocity.
  character(len=80)  :: wave_drag_var  ! The wave drag piston velocity variable
                                       ! name in wave_drag_file.
  real :: vel_rescale ! A rescaling factor for horizontal velocity from the representation in
                      ! a restart file to the internal representation in this run.
  real :: uH_rescale  ! A rescaling factor for thickness transports from the representation in
                      ! a restart file to the internal representation in this run.
  real :: mean_SL     ! The mean sea level that is used along with the bathymetry to estimate the
                      ! geometry when LINEARIZED_BT_CORIOLIS is true or BT_NONLIN_STRESS is false [Z ~> m].
  real, allocatable, dimension(:,:) :: lin_drag_h
  type(memory_size_type) :: MS
  type(group_pass_type) :: pass_static_data, pass_q_D_Cor
  type(group_pass_type) :: pass_bt_hbt_btav, pass_a_polarity
  logical :: default_2018_answers ! The default setting for the various 2018_ANSWERS flags.
  logical :: apply_bt_drag, use_BT_cont_type
  character(len=48) :: thickness_units, flux_units
  character*(40) :: hvel_str
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  integer :: isdw, iedw, jsdw, jedw
  integer :: i, j, k
  integer :: wd_halos(2), bt_halo_sz
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  MS%isdw = G%isd ; MS%iedw = G%ied ; MS%jsdw = G%jsd ; MS%jedw = G%jed

  if (CS%module_is_initialized) then
    call MOM_error(WARNING, "barotropic_init called with a control structure "// &
                            "that has already been initialized.")
    return
  endif
  CS%module_is_initialized = .true.

  CS%diag => diag ; CS%Time => Time
  if (present(tides_CSp)) then
    if (associated(tides_CSp)) CS%tides_CSp => tides_CSp
  endif

  ! Read all relevant parameters and write them to the model log.
  call get_param(param_file, mdl, "SPLIT", CS%split, default=.true., do_not_log=.true.)
  call log_version(param_file, mdl, version, "", log_to_all=.true., layout=CS%split, &
                   debugging=CS%split, all_default=.not.CS%split)
  call get_param(param_file, mdl, "SPLIT", CS%split, &
                 "Use the split time stepping if true.", default=.true.)
  if (.not.CS%split) return

  call get_param(param_file, mdl, "USE_BT_CONT_TYPE", use_BT_cont_type, &
                 "If true, use a structure with elements that describe "//&
                 "effective face areas from the summed continuity solver "//&
                 "as a function the barotropic flow in coupling between "//&
                 "the barotropic and baroclinic flow.  This is only used "//&
                 "if SPLIT is true.", default=.true.)
  call get_param(param_file, mdl, "INTEGRAL_BT_CONTINUITY", CS%integral_bt_cont, &
                 "If true, use the time-integrated velocity over the barotropic steps "//&
                 "to determine the integrated transports used to update the continuity "//&
                 "equation.  Otherwise the transports are the sum of the transports based on "//&
                 "a series of instantaneous velocities and the BT_CONT_TYPE for transports.  "//&
                 "This is only valid if USE_BT_CONT_TYPE = True.", &
                 default=.false., do_not_log=.not.use_BT_cont_type)
  call get_param(param_file, mdl, "BOUND_BT_CORRECTION", CS%bound_BT_corr, &
                 "If true, the corrective pseudo mass-fluxes into the "//&
                 "barotropic solver are limited to values that require "//&
                 "less than maxCFL_BT_cont to be accommodated.",default=.false.)
  call get_param(param_file, mdl, "BT_CONT_CORR_BOUNDS", CS%BT_cont_bounds, &
                 "If true, and BOUND_BT_CORRECTION is true, use the "//&
                 "BT_cont_type variables to set limits determined by "//&
                 "MAXCFL_BT_CONT on the CFL number of the velocities "//&
                 "that are likely to be driven by the corrective mass fluxes.", &
                 default=.true., do_not_log=.not.CS%bound_BT_corr)
  call get_param(param_file, mdl, "ADJUST_BT_CONT", CS%adjust_BT_cont, &
                 "If true, adjust the curve fit to the BT_cont type "//&
                 "that is used by the barotropic solver to match the "//&
                 "transport about which the flow is being linearized.", &
                 default=.false., do_not_log=.not.use_BT_cont_type)
  call get_param(param_file, mdl, "GRADUAL_BT_ICS", CS%gradual_BT_ICs, &
                 "If true, adjust the initial conditions for the "//&
                 "barotropic solver to the values from the layered "//&
                 "solution over a whole timestep instead of instantly. "//&
                 "This is a decent approximation to the inclusion of "//&
                 "sum(u dh_dt) while also correcting for truncation errors.", &
                 default=.false.)
  call get_param(param_file, mdl, "BT_USE_VISC_REM_U_UH0", CS%visc_rem_u_uh0, &
                 "If true, use the viscous remnants when estimating the "//&
                 "barotropic velocities that were used to calculate uh0 "//&
                 "and vh0.  False is probably the better choice.", default=.false.)
  call get_param(param_file, mdl, "BT_USE_WIDE_HALOS", CS%use_wide_halos, &
                 "If true, use wide halos and march in during the "//&
                 "barotropic time stepping for efficiency.", default=.true., &
                 layoutParam=.true.)
  call get_param(param_file, mdl, "BTHALO", bt_halo_sz, &
                 "The minimum halo size for the barotropic solver.", default=0, &
                 layoutParam=.true.)
#ifdef STATIC_MEMORY_
  if ((bt_halo_sz > 0) .and. (bt_halo_sz /= BTHALO_)) call MOM_error(FATAL, &
      "barotropic_init: Run-time values of BTHALO must agree with the "//&
      "macro BTHALO_ with STATIC_MEMORY_.")
  wd_halos(1) = WHALOI_+NIHALO_ ; wd_halos(2) = WHALOJ_+NJHALO_
#else
  wd_halos(1) = bt_halo_sz; wd_halos(2) =  bt_halo_sz
#endif
  call log_param(param_file, mdl, "!BT x-halo", wd_halos(1), &
                 "The barotropic x-halo size that is actually used.", &
                 layoutParam=.true.)
  call log_param(param_file, mdl, "!BT y-halo", wd_halos(2), &
                 "The barotropic y-halo size that is actually used.", &
                 layoutParam=.true.)

  call get_param(param_file, mdl, "NONLINEAR_BT_CONTINUITY", CS%Nonlinear_continuity, &
                 "If true, use nonlinear transports in the barotropic "//&
                 "continuity equation.  This does not apply if "//&
                 "USE_BT_CONT_TYPE is true.", default=.false., do_not_log=use_BT_cont_type)
  call get_param(param_file, mdl, "NONLIN_BT_CONT_UPDATE_PERIOD", CS%Nonlin_cont_update_period, &
                 "If NONLINEAR_BT_CONTINUITY is true, this is the number "//&
                 "of barotropic time steps between updates to the face "//&
                 "areas, or 0 to update only before the barotropic stepping.", &
                 units="nondim", default=1, do_not_log=.not.CS%Nonlinear_continuity)

  call get_param(param_file, mdl, "BT_PROJECT_VELOCITY", CS%BT_project_velocity,&
                 "If true, step the barotropic velocity first and project "//&
                 "out the velocity tendency by 1+BEBT when calculating the "//&
                 "transport.  The default (false) is to use a predictor "//&
                 "continuity step to find the pressure field, and then "//&
                 "to do a corrector continuity step using a weighted "//&
                 "average of the old and new velocities, with weights "//&
                 "of (1-BEBT) and BEBT.", default=.false.)
  call get_param(param_file, mdl, "BT_NONLIN_STRESS", CS%nonlin_stress, &
                 "If true, use the full depth of the ocean at the start of the barotropic "//&
                 "step when calculating the surface stress contribution to the barotropic "//&
                 "acclerations.  Otherwise use the depth based on bathyT.", default=.false.)

  call get_param(param_file, mdl, "DYNAMIC_SURFACE_PRESSURE", CS%dynamic_psurf, &
                 "If true, add a dynamic pressure due to a viscous ice "//&
                 "shelf, for instance.", default=.false.)
  call get_param(param_file, mdl, "ICE_LENGTH_DYN_PSURF", CS%ice_strength_length, &
                 "The length scale at which the Rayleigh damping rate due "//&
                 "to the ice strength should be the same as if a Laplacian "//&
                 "were applied, if DYNAMIC_SURFACE_PRESSURE is true.", &
                 units="m", default=1.0e4, scale=US%m_to_L, do_not_log=.not.CS%dynamic_psurf)
  call get_param(param_file, mdl, "DEPTH_MIN_DYN_PSURF", CS%Dmin_dyn_psurf, &
                 "The minimum depth to use in limiting the size of the "//&
                 "dynamic surface pressure for stability, if "//&
                 "DYNAMIC_SURFACE_PRESSURE is true..", &
                 units="m", default=1.0e-6, scale=US%m_to_Z, do_not_log=.not.CS%dynamic_psurf)
  call get_param(param_file, mdl, "CONST_DYN_PSURF", CS%const_dyn_psurf, &
                 "The constant that scales the dynamic surface pressure, "//&
                 "if DYNAMIC_SURFACE_PRESSURE is true.  Stable values "//&
                 "are < ~1.0.", units="nondim", default=0.9, do_not_log=.not.CS%dynamic_psurf)

  call get_param(param_file, mdl, "BT_CORIOLIS_SCALE", CS%BT_Coriolis_scale, &
                 "A factor by which the barotropic Coriolis anomaly terms are scaled.", &
                 units="nondim", default=1.0)
  call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false.)
  call get_param(param_file, mdl, "BAROTROPIC_2018_ANSWERS", CS%answers_2018, &
                 "If true, use expressions for the barotropic solver that recover the answers "//&
                 "from the end of 2018.  Otherwise, use more efficient or general expressions.", &
                 default=default_2018_answers)

  call get_param(param_file, mdl, "TIDES", CS%tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)
  call get_param(param_file, mdl, "SADOURNY", CS%Sadourny, &
                 "If true, the Coriolis terms are discretized with the "//&
                 "Sadourny (1975) energy conserving scheme, otherwise "//&
                 "the Arakawa & Hsu scheme is used.  If the internal "//&
                 "deformation radius is not resolved, the Sadourny scheme "//&
                 "should probably be used.", default=.true.)

  call get_param(param_file, mdl, "BT_THICK_SCHEME", hvel_str, &
                 "A string describing the scheme that is used to set the "//&
                 "open face areas used for barotropic transport and the "//&
                 "relative weights of the accelerations. Valid values are:\n"//&
                 "\t ARITHMETIC - arithmetic mean layer thicknesses \n"//&
                 "\t HARMONIC - harmonic mean layer thicknesses \n"//&
                 "\t HYBRID (the default) - use arithmetic means for \n"//&
                 "\t    layers above the shallowest bottom, the harmonic \n"//&
                 "\t    mean for layers below, and a weighted average for \n"//&
                 "\t    layers that straddle that depth \n"//&
                 "\t FROM_BT_CONT - use the average thicknesses kept \n"//&
                 "\t    in the h_u and h_v fields of the BT_cont_type", &
                 default=BT_CONT_STRING)
  select case (hvel_str)
    case (HYBRID_STRING) ; CS%hvel_scheme = HYBRID
    case (HARMONIC_STRING) ; CS%hvel_scheme = HARMONIC
    case (ARITHMETIC_STRING) ; CS%hvel_scheme = ARITHMETIC
    case (BT_CONT_STRING) ; CS%hvel_scheme = FROM_BT_CONT
    case default
      call MOM_mesg('barotropic_init: BT_THICK_SCHEME ="'//trim(hvel_str)//'"', 0)
      call MOM_error(FATAL, "barotropic_init: Unrecognized setting "// &
            "#define BT_THICK_SCHEME "//trim(hvel_str)//" found in input file.")
  end select
  if ((CS%hvel_scheme == FROM_BT_CONT) .and. .not.use_BT_cont_type) &
    call MOM_error(FATAL, "barotropic_init: BT_THICK_SCHEME FROM_BT_CONT "//&
                           "can only be used if USE_BT_CONT_TYPE is defined.")

  call get_param(param_file, mdl, "BT_STRONG_DRAG", CS%strong_drag, &
                 "If true, use a stronger estimate of the retarding "//&
                 "effects of strong bottom drag, by making it implicit "//&
                 "with the barotropic time-step instead of implicit with "//&
                 "the baroclinic time-step and dividing by the number of "//&
                 "barotropic steps.", default=.false.)
  call get_param(param_file, mdl, "BT_LINEAR_WAVE_DRAG", CS%linear_wave_drag, &
                 "If true, apply a linear drag to the barotropic velocities, "//&
                 "using rates set by lin_drag_u & _v divided by the depth of "//&
                 "the ocean.  This was introduced to facilitate tide modeling.", &
                 default=.false.)
  call get_param(param_file, mdl, "BT_WAVE_DRAG_FILE", wave_drag_file, &
                 "The name of the file with the barotropic linear wave drag "//&
                 "piston velocities.", default="", do_not_log=.not.CS%linear_wave_drag)
  call get_param(param_file, mdl, "BT_WAVE_DRAG_VAR", wave_drag_var, &
                 "The name of the variable in BT_WAVE_DRAG_FILE with the "//&
                 "barotropic linear wave drag piston velocities at h points.", &
                 default="rH", do_not_log=.not.CS%linear_wave_drag)
  call get_param(param_file, mdl, "BT_WAVE_DRAG_SCALE", wave_drag_scale, &
                 "A scaling factor for the barotropic linear wave drag "//&
                 "piston velocities.", default=1.0, units="nondim", &
                 do_not_log=.not.CS%linear_wave_drag)

  call get_param(param_file, mdl, "CLIP_BT_VELOCITY", CS%clip_velocity, &
                 "If true, limit any velocity components that exceed "//&
                 "CFL_TRUNCATE.  This should only be used as a desperate "//&
                 "debugging measure.", default=.false.)
  call get_param(param_file, mdl, "CFL_TRUNCATE", CS%CFL_trunc, &
                 "The value of the CFL number that will cause velocity "//&
                 "components to be truncated; instability can occur past 0.5.", &
                 units="nondim", default=0.5, do_not_log=.not.CS%clip_velocity)
  call get_param(param_file, mdl, "MAXVEL", CS%maxvel, &
                 "The maximum velocity allowed before the velocity "//&
                 "components are truncated.", units="m s-1", default=3.0e8, scale=US%m_s_to_L_T, &
                 do_not_log=.not.CS%clip_velocity)
  call get_param(param_file, mdl, "MAXCFL_BT_CONT", CS%maxCFL_BT_cont, &
                 "The maximum permitted CFL number associated with the "//&
                 "barotropic accelerations from the summed velocities "//&
                 "times the time-derivatives of thicknesses.", units="nondim", &
                 default=0.25)
  call get_param(param_file, mdl, "VEL_UNDERFLOW", CS%vel_underflow, &
                 "A negligibly small velocity magnitude below which velocity "//&
                 "components are set to 0.  A reasonable value might be "//&
                 "1e-30 m/s, which is less than an Angstrom divided by "//&
                 "the age of the universe.", units="m s-1", default=0.0, scale=US%m_s_to_L_T)

  call get_param(param_file, mdl, "DT_BT_FILTER", CS%dt_bt_filter, &
                 "A time-scale over which the barotropic mode solutions "//&
                 "are filtered, in seconds if positive, or as a fraction "//&
                 "of DT if negative. When used this can never be taken to "//&
                 "be longer than 2*dt.  Set this to 0 to apply no filtering.", &
                 units="sec or nondim", default=-0.25)
  if (CS%dt_bt_filter > 0.0) CS%dt_bt_filter = US%s_to_T*CS%dt_bt_filter
  call get_param(param_file, mdl, "G_BT_EXTRA", CS%G_extra, &
                 "A nondimensional factor by which gtot is enhanced.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "SSH_EXTRA", SSH_extra, &
                 "An estimate of how much higher SSH might get, for use "//&
                 "in calculating the safe external wave speed. The "//&
                 "default is the minimum of 10 m or 5% of MAXIMUM_DEPTH.", &
                 units="m", default=min(10.0,0.05*G%max_depth*US%Z_to_m), scale=US%m_to_Z)

  call get_param(param_file, mdl, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_BT", CS%debug_bt, &
                 "If true, write out verbose debugging data within the "//&
                 "barotropic time-stepping loop. The data volume can be "//&
                 "quite large if this is true.", default=CS%debug, &
                 debuggingParam=.true.)

  call get_param(param_file, mdl, "LINEARIZED_BT_CORIOLIS", CS%linearized_BT_PV, &
                 "If true use the bottom depth instead of the total water column thickness "//&
                 "in the barotropic Coriolis term calculations.", default=.true.)
  call get_param(param_file, mdl, "BEBT", CS%bebt, &
                 "BEBT determines whether the barotropic time stepping "//&
                 "uses the forward-backward time-stepping scheme or a "//&
                 "backward Euler scheme. BEBT is valid in the range from "//&
                 "0 (for a forward-backward treatment of nonrotating "//&
                 "gravity waves) to 1 (for a backward Euler treatment). "//&
                 "In practice, BEBT must be greater than about 0.05.", &
                 units="nondim", default=0.1)
  call get_param(param_file, mdl, "DTBT", dtbt_input, &
                 "The barotropic time step, in s. DTBT is only used with "//&
                 "the split explicit time stepping. To set the time step "//&
                 "automatically based the maximum stable value use 0, or "//&
                 "a negative value gives the fraction of the stable value. "//&
                 "Setting DTBT to 0 is the same as setting it to -0.98. "//&
                 "The value of DTBT that will actually be used is an "//&
                 "integer fraction of DT, rounding down.", units="s or nondim",&
                 default = -0.98)
  call get_param(param_file, mdl, "BT_USE_OLD_CORIOLIS_BRACKET_BUG", &
                 CS%use_old_coriolis_bracket_bug , &
                 "If True, use an order of operations that is not bitwise "//&
                 "rotationally symmetric in the meridional Coriolis term of "//&
                 "the barotropic solver.", default=.false.)

  ! Initialize a version of the MOM domain that is specific to the barotropic solver.
  call clone_MOM_domain(G%Domain, CS%BT_Domain, min_halo=wd_halos, symmetric=.true.)
#ifdef STATIC_MEMORY_
  if (wd_halos(1) /= WHALOI_+NIHALO_) call MOM_error(FATAL, "barotropic_init: "//&
          "Barotropic x-halo sizes are incorrectly resized with STATIC_MEMORY_.")
  if (wd_halos(2) /= WHALOJ_+NJHALO_) call MOM_error(FATAL, "barotropic_init: "//&
          "Barotropic y-halo sizes are incorrectly resized with STATIC_MEMORY_.")
#else
  if (bt_halo_sz > 0) then
    if (wd_halos(1) > bt_halo_sz) &
      call MOM_mesg("barotropic_init: barotropic x-halo size increased.", 3)
    if (wd_halos(2) > bt_halo_sz) &
      call MOM_mesg("barotropic_init: barotropic y-halo size increased.", 3)
  endif
#endif

  CS%isdw = G%isc-wd_halos(1) ; CS%iedw = G%iec+wd_halos(1)
  CS%jsdw = G%jsc-wd_halos(2) ; CS%jedw = G%jec+wd_halos(2)
  isdw = CS%isdw ; iedw = CS%iedw ; jsdw = CS%jsdw ; jedw = CS%jedw

  ALLOC_(CS%frhatu(IsdB:IedB,jsd:jed,nz)) ; ALLOC_(CS%frhatv(isd:ied,JsdB:JedB,nz))
  ALLOC_(CS%eta_cor(isd:ied,jsd:jed))
  if (CS%bound_BT_corr) then
    ALLOC_(CS%eta_cor_bound(isd:ied,jsd:jed)) ; CS%eta_cor_bound(:,:) = 0.0
  endif
  ALLOC_(CS%IDatu(IsdB:IedB,jsd:jed)) ; ALLOC_(CS%IDatv(isd:ied,JsdB:JedB))

  ALLOC_(CS%ua_polarity(isdw:iedw,jsdw:jedw))
  ALLOC_(CS%va_polarity(isdw:iedw,jsdw:jedw))

  CS%frhatu(:,:,:) = 0.0 ; CS%frhatv(:,:,:) = 0.0
  CS%eta_cor(:,:) = 0.0
  CS%IDatu(:,:) = 0.0 ; CS%IDatv(:,:) = 0.0

  CS%ua_polarity(:,:) = 1.0 ; CS%va_polarity(:,:) = 1.0
  call create_group_pass(pass_a_polarity, CS%ua_polarity, CS%va_polarity, CS%BT_domain, To_All, AGRID)
  call do_group_pass(pass_a_polarity, CS%BT_domain)

  if (use_BT_cont_type) &
    call alloc_BT_cont_type(BT_cont, G, GV, (CS%hvel_scheme == FROM_BT_CONT))

  if (CS%debug) then ! Make a local copy of loop ranges for chksum calls
    allocate(CS%debug_BT_HI)
    CS%debug_BT_HI%isc=G%isc
    CS%debug_BT_HI%iec=G%iec
    CS%debug_BT_HI%jsc=G%jsc
    CS%debug_BT_HI%jec=G%jec
    CS%debug_BT_HI%IscB=G%isc-1
    CS%debug_BT_HI%IecB=G%iec
    CS%debug_BT_HI%JscB=G%jsc-1
    CS%debug_BT_HI%JecB=G%jec
    CS%debug_BT_HI%isd=CS%isdw
    CS%debug_BT_HI%ied=CS%iedw
    CS%debug_BT_HI%jsd=CS%jsdw
    CS%debug_BT_HI%jed=CS%jedw
    CS%debug_BT_HI%IsdB=CS%isdw-1
    CS%debug_BT_HI%IedB=CS%iedw
    CS%debug_BT_HI%JsdB=CS%jsdw-1
    CS%debug_BT_HI%JedB=CS%jedw
    CS%debug_BT_HI%turns = G%HI%turns
  endif

  ! IareaT, IdxCu, and IdyCv need to be allocated with wide halos.
  ALLOC_(CS%IareaT(CS%isdw:CS%iedw,CS%jsdw:CS%jedw)) ; CS%IareaT(:,:) = 0.0
  ALLOC_(CS%bathyT(CS%isdw:CS%iedw,CS%jsdw:CS%jedw)) ; CS%bathyT(:,:) = GV%Angstrom_m !### Change to 0.0?
  ALLOC_(CS%IdxCu(CS%isdw-1:CS%iedw,CS%jsdw:CS%jedw)) ; CS%IdxCu(:,:) = 0.0
  ALLOC_(CS%IdyCv(CS%isdw:CS%iedw,CS%jsdw-1:CS%jedw)) ; CS%IdyCv(:,:) = 0.0
  ALLOC_(CS%dy_Cu(CS%isdw-1:CS%iedw,CS%jsdw:CS%jedw)) ; CS%dy_Cu(:,:) = 0.0
  ALLOC_(CS%dx_Cv(CS%isdw:CS%iedw,CS%jsdw-1:CS%jedw)) ; CS%dx_Cv(:,:) = 0.0
  do j=G%jsd,G%jed ; do i=G%isd,G%ied
    CS%IareaT(i,j) = G%IareaT(i,j)
    CS%bathyT(i,j) = G%bathyT(i,j)
  enddo ; enddo

  ! Note: G%IdxCu & G%IdyCv may be valid for a smaller extent than CS%IdxCu & CS%IdyCv, even without
  !   wide halos.
  do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB
    CS%IdxCu(I,j) = G%IdxCu(I,j) ; CS%dy_Cu(I,j) = G%dy_Cu(I,j)
  enddo ; enddo
  do J=G%JsdB,G%JedB ; do i=G%isd,G%ied
    CS%IdyCv(i,J) = G%IdyCv(i,J) ; CS%dx_Cv(i,J) = G%dx_Cv(i,J)
  enddo ; enddo
  call create_group_pass(pass_static_data, CS%IareaT, CS%BT_domain, To_All)
  call create_group_pass(pass_static_data, CS%bathyT, CS%BT_domain, To_All)
  call create_group_pass(pass_static_data, CS%IdxCu, CS%IdyCv, CS%BT_domain, To_All+Scalar_Pair)
  call create_group_pass(pass_static_data, CS%dy_Cu, CS%dx_Cv, CS%BT_domain, To_All+Scalar_Pair)
  call do_group_pass(pass_static_data, CS%BT_domain)

  if (CS%linearized_BT_PV) then
    ALLOC_(CS%q_D(CS%isdw-1:CS%iedw,CS%jsdw-1:CS%jedw))
    ALLOC_(CS%D_u_Cor(CS%isdw-1:CS%iedw,CS%jsdw:CS%jedw))
    ALLOC_(CS%D_v_Cor(CS%isdw:CS%iedw,CS%jsdw-1:CS%jedw))
    CS%q_D(:,:) = 0.0 ; CS%D_u_Cor(:,:) = 0.0 ; CS%D_v_Cor(:,:) = 0.0

    Mean_SL = G%Z_ref
    do j=js,je ; do I=is-1,ie
      CS%D_u_Cor(I,j) = 0.5 * (max(Mean_SL+G%bathyT(i+1,j),0.0) + max(Mean_SL+G%bathyT(i,j),0.0))
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      CS%D_v_Cor(i,J) = 0.5 * (max(Mean_SL+G%bathyT(i,j+1),0.0) + max(Mean_SL+G%bathyT(i,j),0.0))
    enddo ; enddo
    do J=js-1,je ; do I=is-1,ie
      if (G%mask2dT(i,j)+G%mask2dT(i,j+1)+G%mask2dT(i+1,j)+G%mask2dT(i+1,j+1)>0.) then
        CS%q_D(I,J) = 0.25 * (CS%BT_Coriolis_scale * G%CoriolisBu(I,J)) * &
           ((G%areaT(i,j) + G%areaT(i+1,j+1)) + (G%areaT(i+1,j) + G%areaT(i,j+1))) / &
           (max(((G%areaT(i,j) * max(Mean_SL+G%bathyT(i,j),0.0) + &
                  G%areaT(i+1,j+1) * max(Mean_SL+G%bathyT(i+1,j+1),0.0)) + &
                 (G%areaT(i+1,j) * max(Mean_SL+G%bathyT(i+1,j),0.0) + &
                  G%areaT(i,j+1) * max(Mean_SL+G%bathyT(i,j+1),0.0))), GV%H_to_Z*GV%H_subroundoff) )
      else ! All four h points are masked out so q_D(I,J) will is meaningless
        CS%q_D(I,J) = 0.
      endif
    enddo ; enddo
    ! With very wide halos, q and D need to be calculated on the available data
    ! domain and then updated onto the full computational domain.
    call create_group_pass(pass_q_D_Cor, CS%q_D, CS%BT_Domain, To_All, position=CORNER)
    call create_group_pass(pass_q_D_Cor, CS%D_u_Cor, CS%D_v_Cor, CS%BT_Domain, &
                           To_All+Scalar_Pair)
    call do_group_pass(pass_q_D_Cor, CS%BT_Domain)
  endif

  if (CS%linear_wave_drag) then
    ALLOC_(CS%lin_drag_u(IsdB:IedB,jsd:jed)) ; CS%lin_drag_u(:,:) = 0.0
    ALLOC_(CS%lin_drag_v(isd:ied,JsdB:JedB)) ; CS%lin_drag_v(:,:) = 0.0

    if (len_trim(wave_drag_file) > 0) then
      inputdir = "." ;  call get_param(param_file, mdl, "INPUTDIR", inputdir)
      wave_drag_file = trim(slasher(inputdir))//trim(wave_drag_file)
      call log_param(param_file, mdl, "INPUTDIR/BT_WAVE_DRAG_FILE", wave_drag_file)

      allocate(lin_drag_h(isd:ied,jsd:jed)) ; lin_drag_h(:,:) = 0.0

      call MOM_read_data(wave_drag_file, wave_drag_var, lin_drag_h, G%Domain, scale=US%m_to_Z*US%T_to_s)
      call pass_var(lin_drag_h, G%Domain)
      do j=js,je ; do I=is-1,ie
        CS%lin_drag_u(I,j) = (GV%Z_to_H * wave_drag_scale) * &
           0.5 * (lin_drag_h(i,j) + lin_drag_h(i+1,j))
      enddo ; enddo
      do J=js-1,je ; do i=is,ie
        CS%lin_drag_v(i,J) = (GV%Z_to_H * wave_drag_scale) * &
            0.5 * (lin_drag_h(i,j) + lin_drag_h(i,j+1))
      enddo ; enddo
      deallocate(lin_drag_h)
    endif
  endif

  CS%dtbt_fraction = 0.98 ; if (dtbt_input < 0.0) CS%dtbt_fraction = -dtbt_input

  dtbt_tmp = -1.0
  if (query_initialized(CS%dtbt, "DTBT", restart_CS)) then
    dtbt_tmp = CS%dtbt
    if ((US%s_to_T_restart /= 0.0) .and. (US%s_to_T_restart /= US%s_to_T)) &
      dtbt_tmp = (US%s_to_T / US%s_to_T_restart) * CS%dtbt
  endif

  ! Estimate the maximum stable barotropic time step.
  gtot_estimate = 0.0
  do k=1,GV%ke ; gtot_estimate = gtot_estimate + GV%g_prime(K) ; enddo
  call set_dtbt(G, GV, US, CS, gtot_est=gtot_estimate, SSH_add=SSH_extra)

  if (dtbt_input > 0.0) then
    CS%dtbt = US%s_to_T * dtbt_input
  elseif (dtbt_tmp > 0.0) then
    CS%dtbt = dtbt_tmp
  endif
  if ((dtbt_tmp > 0.0) .and. (dtbt_input > 0.0)) calc_dtbt = .false.

  call log_param(param_file, mdl, "DTBT as used", CS%dtbt*US%T_to_s)
  call log_param(param_file, mdl, "estimated maximum DTBT", CS%dtbt_max*US%T_to_s)

  ! ubtav and vbtav, and perhaps ubt_IC and vbt_IC, are allocated and
  ! initialized in register_barotropic_restarts.

  if (GV%Boussinesq) then
    thickness_units = "m" ; flux_units = "m3 s-1"
  else
    thickness_units = "kg m-2" ; flux_units = "kg s-1"
  endif

  CS%id_PFu_bt = register_diag_field('ocean_model', 'PFuBT', diag%axesCu1, Time, &
      'Zonal Anomalous Barotropic Pressure Force Force Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_PFv_bt = register_diag_field('ocean_model', 'PFvBT', diag%axesCv1, Time, &
      'Meridional Anomalous Barotropic Pressure Force Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_Coru_bt = register_diag_field('ocean_model', 'CoruBT', diag%axesCu1, Time, &
      'Zonal Barotropic Coriolis Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_Corv_bt = register_diag_field('ocean_model', 'CorvBT', diag%axesCv1, Time, &
      'Meridional Barotropic Coriolis Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_uaccel = register_diag_field('ocean_model', 'u_accel_bt', diag%axesCu1, Time, &
      'Barotropic zonal acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_vaccel = register_diag_field('ocean_model', 'v_accel_bt', diag%axesCv1, Time, &
      'Barotropic meridional acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_ubtforce = register_diag_field('ocean_model', 'ubtforce', diag%axesCu1, Time, &
      'Barotropic zonal acceleration from baroclinic terms', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_vbtforce = register_diag_field('ocean_model', 'vbtforce', diag%axesCv1, Time, &
      'Barotropic meridional acceleration from baroclinic terms', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_ubtdt = register_diag_field('ocean_model', 'ubt_dt', diag%axesCu1, Time, &
      'Barotropic zonal acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_vbtdt = register_diag_field('ocean_model', 'vbt_dt', diag%axesCv1, Time, &
      'Barotropic meridional acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_eta_bt = register_diag_field('ocean_model', 'eta_bt', diag%axesT1, Time, &
      'Barotropic end SSH', thickness_units, conversion=GV%H_to_m)
  CS%id_ubt = register_diag_field('ocean_model', 'ubt', diag%axesCu1, Time, &
      'Barotropic end zonal velocity', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vbt = register_diag_field('ocean_model', 'vbt', diag%axesCv1, Time, &
      'Barotropic end meridional velocity', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_eta_st = register_diag_field('ocean_model', 'eta_st', diag%axesT1, Time, &
      'Barotropic start SSH', thickness_units, conversion=GV%H_to_m)
  CS%id_ubt_st = register_diag_field('ocean_model', 'ubt_st', diag%axesCu1, Time, &
      'Barotropic start zonal velocity', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vbt_st = register_diag_field('ocean_model', 'vbt_st', diag%axesCv1, Time, &
      'Barotropic start meridional velocity', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_ubtav = register_diag_field('ocean_model', 'ubtav', diag%axesCu1, Time, &
      'Barotropic time-average zonal velocity', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vbtav = register_diag_field('ocean_model', 'vbtav', diag%axesCv1, Time, &
      'Barotropic time-average meridional velocity', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_eta_cor = register_diag_field('ocean_model', 'eta_cor', diag%axesT1, Time, &
      'Corrective mass flux', 'm s-1', conversion=GV%H_to_m)
  CS%id_visc_rem_u = register_diag_field('ocean_model', 'visc_rem_u', diag%axesCuL, Time, &
      'Viscous remnant at u', 'nondim')
  CS%id_visc_rem_v = register_diag_field('ocean_model', 'visc_rem_v', diag%axesCvL, Time, &
      'Viscous remnant at v', 'nondim')
  CS%id_gtotn = register_diag_field('ocean_model', 'gtot_n', diag%axesT1, Time, &
      'gtot to North', 'm s-2', conversion=GV%m_to_H*(US%L_T_to_m_s**2))
  CS%id_gtots = register_diag_field('ocean_model', 'gtot_s', diag%axesT1, Time, &
      'gtot to South', 'm s-2', conversion=GV%m_to_H*(US%L_T_to_m_s**2))
  CS%id_gtote = register_diag_field('ocean_model', 'gtot_e', diag%axesT1, Time, &
      'gtot to East', 'm s-2', conversion=GV%m_to_H*(US%L_T_to_m_s**2))
  CS%id_gtotw = register_diag_field('ocean_model', 'gtot_w', diag%axesT1, Time, &
      'gtot to West', 'm s-2', conversion=GV%m_to_H*(US%L_T_to_m_s**2))
  CS%id_eta_hifreq = register_diag_field('ocean_model', 'eta_hifreq', diag%axesT1, Time, &
      'High Frequency Barotropic SSH', thickness_units, conversion=GV%H_to_m)
  CS%id_ubt_hifreq = register_diag_field('ocean_model', 'ubt_hifreq', diag%axesCu1, Time, &
      'High Frequency Barotropic zonal velocity', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vbt_hifreq = register_diag_field('ocean_model', 'vbt_hifreq', diag%axesCv1, Time, &
      'High Frequency Barotropic meridional velocity', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_eta_pred_hifreq = register_diag_field('ocean_model', 'eta_pred_hifreq', diag%axesT1, Time, &
      'High Frequency Predictor Barotropic SSH', thickness_units, &
      conversion=GV%H_to_m)
  CS%id_uhbt_hifreq = register_diag_field('ocean_model', 'uhbt_hifreq', diag%axesCu1, Time, &
      'High Frequency Barotropic zonal transport', 'm3 s-1', &
      conversion=GV%H_to_m*US%L_to_m*US%L_T_to_m_s)
  CS%id_vhbt_hifreq = register_diag_field('ocean_model', 'vhbt_hifreq', diag%axesCv1, Time, &
      'High Frequency Barotropic meridional transport', 'm3 s-1', &
      conversion=GV%H_to_m*US%L_to_m*US%L_T_to_m_s)
  CS%id_frhatu = register_diag_field('ocean_model', 'frhatu', diag%axesCuL, Time, &
      'Fractional thickness of layers in u-columns', 'nondim')
  CS%id_frhatv = register_diag_field('ocean_model', 'frhatv', diag%axesCvL, Time, &
      'Fractional thickness of layers in v-columns', 'nondim')
  CS%id_frhatu1 = register_diag_field('ocean_model', 'frhatu1', diag%axesCuL, Time, &
      'Predictor Fractional thickness of layers in u-columns', 'nondim')
  CS%id_frhatv1 = register_diag_field('ocean_model', 'frhatv1', diag%axesCvL, Time, &
      'Predictor Fractional thickness of layers in v-columns', 'nondim')
  CS%id_uhbt = register_diag_field('ocean_model', 'uhbt', diag%axesCu1, Time, &
      'Barotropic zonal transport averaged over a baroclinic step', 'm3 s-1', &
      conversion=GV%H_to_m*US%L_to_m*US%L_T_to_m_s)
  CS%id_vhbt = register_diag_field('ocean_model', 'vhbt', diag%axesCv1, Time, &
      'Barotropic meridional transport averaged over a baroclinic step', 'm3 s-1', &
      conversion=GV%H_to_m*US%L_to_m*US%L_T_to_m_s)

  if (use_BT_cont_type) then
    CS%id_BTC_FA_u_EE = register_diag_field('ocean_model', 'BTC_FA_u_EE', diag%axesCu1, Time, &
        'BTCont type far east face area', 'm2', conversion=US%L_to_m*GV%H_to_m)
    CS%id_BTC_FA_u_E0 = register_diag_field('ocean_model', 'BTC_FA_u_E0', diag%axesCu1, Time, &
        'BTCont type near east face area', 'm2', conversion=US%L_to_m*GV%H_to_m)
    CS%id_BTC_FA_u_WW = register_diag_field('ocean_model', 'BTC_FA_u_WW', diag%axesCu1, Time, &
        'BTCont type far west face area', 'm2', conversion=US%L_to_m*GV%H_to_m)
    CS%id_BTC_FA_u_W0 = register_diag_field('ocean_model', 'BTC_FA_u_W0', diag%axesCu1, Time, &
        'BTCont type near west face area', 'm2', conversion=US%L_to_m*GV%H_to_m)
    CS%id_BTC_ubt_EE = register_diag_field('ocean_model', 'BTC_ubt_EE', diag%axesCu1, Time, &
        'BTCont type far east velocity', 'm s-1', conversion=US%L_T_to_m_s)
    CS%id_BTC_ubt_WW = register_diag_field('ocean_model', 'BTC_ubt_WW', diag%axesCu1, Time, &
        'BTCont type far west velocity', 'm s-1', conversion=US%L_T_to_m_s)
    ! This is a specialized diagnostic that is not being made widely available (yet).
    ! CS%id_BTC_FA_u_rat0 = register_diag_field('ocean_model', 'BTC_FA_u_rat0', diag%axesCu1, Time, &
    !     'BTCont type ratio of near east and west face areas', 'nondim')
    CS%id_BTC_FA_v_NN = register_diag_field('ocean_model', 'BTC_FA_v_NN', diag%axesCv1, Time, &
        'BTCont type far north face area', 'm2', conversion=US%L_to_m*GV%H_to_m)
    CS%id_BTC_FA_v_N0 = register_diag_field('ocean_model', 'BTC_FA_v_N0', diag%axesCv1, Time, &
        'BTCont type near north face area', 'm2', conversion=US%L_to_m*GV%H_to_m)
    CS%id_BTC_FA_v_SS = register_diag_field('ocean_model', 'BTC_FA_v_SS', diag%axesCv1, Time, &
        'BTCont type far south face area', 'm2', conversion=US%L_to_m*GV%H_to_m)
    CS%id_BTC_FA_v_S0 = register_diag_field('ocean_model', 'BTC_FA_v_S0', diag%axesCv1, Time, &
        'BTCont type near south face area', 'm2', conversion=US%L_to_m*GV%H_to_m)
    CS%id_BTC_vbt_NN = register_diag_field('ocean_model', 'BTC_vbt_NN', diag%axesCv1, Time, &
        'BTCont type far north velocity', 'm s-1', conversion=US%L_T_to_m_s)
    CS%id_BTC_vbt_SS = register_diag_field('ocean_model', 'BTC_vbt_SS', diag%axesCv1, Time, &
        'BTCont type far south velocity', 'm s-1', conversion=US%L_T_to_m_s)
    ! This is a specialized diagnostic that is not being made widely available (yet).
    ! CS%id_BTC_FA_v_rat0 = register_diag_field('ocean_model', 'BTC_FA_v_rat0', diag%axesCv1, Time, &
    !     'BTCont type ratio of near north and south face areas', 'nondim')
    ! CS%id_BTC_FA_h_rat0 = register_diag_field('ocean_model', 'BTC_FA_h_rat0', diag%axesT1, Time, &
    !     'BTCont type maximum ratios of near face areas around cells', 'nondim')
  endif
  CS%id_uhbt0 = register_diag_field('ocean_model', 'uhbt0', diag%axesCu1, Time, &
      'Barotropic zonal transport difference', 'm3 s-1', conversion=GV%H_to_m*US%L_to_m**2*US%s_to_T)
  CS%id_vhbt0 = register_diag_field('ocean_model', 'vhbt0', diag%axesCv1, Time, &
      'Barotropic meridional transport difference', 'm3 s-1', conversion=GV%H_to_m*US%L_to_m**2*US%s_to_T)

  if (CS%id_frhatu1 > 0) call safe_alloc_ptr(CS%frhatu1, IsdB,IedB,jsd,jed,nz)
  if (CS%id_frhatv1 > 0) call safe_alloc_ptr(CS%frhatv1, isd,ied,JsdB,JedB,nz)

  if (.NOT.query_initialized(CS%ubtav,"ubtav",restart_CS) .or. &
      .NOT.query_initialized(CS%vbtav,"vbtav",restart_CS)) then
    call btcalc(h, G, GV, CS, may_use_default=.true.)
    CS%ubtav(:,:) = 0.0 ; CS%vbtav(:,:) = 0.0
    do k=1,nz ; do j=js,je ; do I=is-1,ie
      CS%ubtav(I,j) = CS%ubtav(I,j) + CS%frhatu(I,j,k) * u(I,j,k)
    enddo ; enddo ; enddo
    do k=1,nz ; do J=js-1,je ; do i=is,ie
      CS%vbtav(i,J) = CS%vbtav(i,J) + CS%frhatv(i,J,k) * v(i,J,k)
    enddo ; enddo ; enddo
  elseif ((US%s_to_T_restart*US%m_to_L_restart /= 0.0) .and. &
          (US%m_to_L*US%s_to_T_restart) /= (US%m_to_L_restart*US%s_to_T)) then
    vel_rescale = (US%m_to_L*US%s_to_T_restart) / (US%m_to_L_restart*US%s_to_T)
    do j=js,je ; do I=is-1,ie ; CS%ubtav(I,j) = vel_rescale * CS%ubtav(I,j) ; enddo ; enddo
    do J=js-1,je ; do i=is,ie ; CS%vbtav(i,J) = vel_rescale * CS%vbtav(i,J) ; enddo ; enddo
  endif

  if (CS%gradual_BT_ICs) then
    if (.NOT.query_initialized(CS%ubt_IC,"ubt_IC",restart_CS) .or. &
        .NOT.query_initialized(CS%vbt_IC,"vbt_IC",restart_CS)) then
      do j=js,je ; do I=is-1,ie ; CS%ubt_IC(I,j) = CS%ubtav(I,j) ; enddo ; enddo
      do J=js-1,je ; do i=is,ie ; CS%vbt_IC(i,J) = CS%vbtav(i,J) ; enddo ; enddo
    elseif ((US%s_to_T_restart*US%m_to_L_restart /= 0.0) .and. &
            (US%m_to_L*US%s_to_T_restart) /= (US%m_to_L_restart*US%s_to_T)) then
      vel_rescale = (US%m_to_L*US%s_to_T_restart) / (US%m_to_L_restart*US%s_to_T)
      do j=js,je ; do I=is-1,ie ; CS%ubt_IC(I,j) = vel_rescale * CS%ubt_IC(I,j) ; enddo ; enddo
      do J=js-1,je ; do i=is,ie ; CS%vbt_IC(i,J) = vel_rescale * CS%vbt_IC(i,J) ; enddo ; enddo
    endif
  endif
!   Calculate other constants which are used for btstep.

  if (.not.CS%nonlin_stress) then
    Mean_SL = G%Z_ref
    do j=js,je ; do I=is-1,ie
      if (G%mask2dCu(I,j)>0.) then
        CS%IDatu(I,j) = G%mask2dCu(I,j) * 2.0 / ((G%bathyT(i+1,j) + G%bathyT(i,j)) + 2.0*Mean_SL)
      else ! Both neighboring H points are masked out so IDatu(I,j) is meaningless
        CS%IDatu(I,j) = 0.
      endif
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      if (G%mask2dCv(i,J)>0.) then
        CS%IDatv(i,J) = G%mask2dCv(i,J) * 2.0 / ((G%bathyT(i,j+1) + G%bathyT(i,j)) + 2.0*Mean_SL)
      else ! Both neighboring H points are masked out so IDatv(i,J) is meaningless
        CS%IDatv(i,J) = 0.
      endif
    enddo ; enddo
  endif

  call find_face_areas(Datu, Datv, G, GV, US, CS, MS, halo=1)
  if ((CS%bound_BT_corr) .and. .not.(use_BT_Cont_type .and. CS%BT_cont_bounds)) then
    ! This is not used in most test cases.  Were it ever to become more widely used, consider
    ! replacing maxvel with min(G%dxT(i,j),G%dyT(i,j)) * (CS%maxCFL_BT_cont*Idt) .
    do j=js,je ; do i=is,ie
      CS%eta_cor_bound(i,j) = G%IareaT(i,j) * 0.1 * CS%maxvel * &
         ((Datu(I-1,j) + Datu(I,j)) + (Datv(i,J) + Datv(i,J-1)))
    enddo ; enddo
  endif

  if (CS%gradual_BT_ICs) &
    call create_group_pass(pass_bt_hbt_btav, CS%ubt_IC, CS%vbt_IC, G%Domain)
  call create_group_pass(pass_bt_hbt_btav, CS%ubtav, CS%vbtav, G%Domain)
  call do_group_pass(pass_bt_hbt_btav, G%Domain)

!  id_clock_pass = cpu_clock_id('(Ocean BT halo updates)', grain=CLOCK_ROUTINE)
  id_clock_calc_pre  = cpu_clock_id('(Ocean BT pre-calcs only)', grain=CLOCK_ROUTINE)
  id_clock_pass_pre = cpu_clock_id('(Ocean BT pre-step halo updates)', grain=CLOCK_ROUTINE)
  id_clock_calc = cpu_clock_id('(Ocean BT stepping calcs only)', grain=CLOCK_ROUTINE)
  id_clock_pass_step = cpu_clock_id('(Ocean BT stepping halo updates)', grain=CLOCK_ROUTINE)
  id_clock_calc_post = cpu_clock_id('(Ocean BT post-calcs only)', grain=CLOCK_ROUTINE)
  id_clock_pass_post = cpu_clock_id('(Ocean BT post-step halo updates)', grain=CLOCK_ROUTINE)
  if (dtbt_input <= 0.0) &
    id_clock_sync = cpu_clock_id('(Ocean BT global synch)', grain=CLOCK_ROUTINE)

end subroutine barotropic_init

!> Copies ubtav and vbtav from private type into arrays
subroutine barotropic_get_tav(CS, ubtav, vbtav, G, US)
  type(barotropic_CS),               pointer       :: CS    !< Control structure for this module
  type(ocean_grid_type),             intent(in)    :: G     !< Grid structure
  real, dimension(SZIB_(G),SZJ_(G)), intent(inout) :: ubtav !< Zonal barotropic velocity averaged
                                                            !! over a baroclinic timestep [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G)), intent(inout) :: vbtav !< Meridional barotropic velocity averaged
                                                            !! over a baroclinic timestep [L T-1 ~> m s-1]
  type(unit_scale_type),             intent(in)    :: US    !< A dimensional unit scaling type
  ! Local variables
  integer :: i,j

  do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
    ubtav(I,j) = CS%ubtav(I,j)
  enddo ; enddo

  do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
    vbtav(i,J) = CS%vbtav(i,J)
  enddo ; enddo

end subroutine barotropic_get_tav


!> Clean up the barotropic control structure.
subroutine barotropic_end(CS)
  type(barotropic_CS), pointer :: CS  !< Control structure to clear out.
  DEALLOC_(CS%frhatu)   ; DEALLOC_(CS%frhatv)
  DEALLOC_(CS%IDatu)    ; DEALLOC_(CS%IDatv)
  DEALLOC_(CS%ubtav)    ; DEALLOC_(CS%vbtav)
  DEALLOC_(CS%eta_cor)
  DEALLOC_(CS%ua_polarity) ; DEALLOC_(CS%va_polarity)
  if (CS%bound_BT_corr) then
    DEALLOC_(CS%eta_cor_bound)
  endif

  call destroy_BT_OBC(CS%BT_OBC)

  deallocate(CS)
end subroutine barotropic_end

!> This subroutine is used to register any fields from MOM_barotropic.F90
!! that should be written to or read from the restart file.
subroutine register_barotropic_restarts(HI, GV, param_file, CS, restart_CS)
  type(hor_index_type),    intent(in) :: HI         !< A horizontal index type structure.
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters.
  type(barotropic_CS),     pointer    :: CS         !< A pointer that is set to point to the control
                                                    !! structure for this module.
  type(verticalGrid_type), intent(in) :: GV         !< The ocean's vertical grid structure.
  type(MOM_restart_CS),    pointer    :: restart_CS !< A pointer to the restart control structure.

  ! Local variables
  type(vardesc) :: vd(3)
  character(len=40)  :: mdl = "MOM_barotropic"  ! This module's name.
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed
  IsdB = HI%IsdB ; IedB = HI%IedB ; JsdB = HI%JsdB ; JedB = HI%JedB

  if (associated(CS)) then
    call MOM_error(WARNING, "register_barotropic_restarts called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  call get_param(param_file, mdl, "GRADUAL_BT_ICS", CS%gradual_BT_ICs, &
                 "If true, adjust the initial conditions for the "//&
                 "barotropic solver to the values from the layered "//&
                 "solution over a whole timestep instead of instantly. "//&
                 "This is a decent approximation to the inclusion of "//&
                 "sum(u dh_dt) while also correcting for truncation errors.", &
                 default=.false., do_not_log=.true.)

  ALLOC_(CS%ubtav(IsdB:IedB,jsd:jed))      ; CS%ubtav(:,:) = 0.0
  ALLOC_(CS%vbtav(isd:ied,JsdB:JedB))      ; CS%vbtav(:,:) = 0.0
  if (CS%gradual_BT_ICs) then
    ALLOC_(CS%ubt_IC(IsdB:IedB,jsd:jed))     ; CS%ubt_IC(:,:) = 0.0
    ALLOC_(CS%vbt_IC(isd:ied,JsdB:JedB))     ; CS%vbt_IC(:,:) = 0.0
  endif

  vd(2) = var_desc("ubtav","m s-1","Time mean barotropic zonal velocity", &
                hor_grid='u', z_grid='1')
  vd(3) = var_desc("vbtav","m s-1","Time mean barotropic meridional velocity",&
                hor_grid='v', z_grid='1')
  call register_restart_pair(CS%ubtav, CS%vbtav, vd(2), vd(3), .false., restart_CS)

  if (CS%gradual_BT_ICs) then
    vd(2) = var_desc("ubt_IC", "m s-1", &
                longname="Next initial condition for the barotropic zonal velocity", &
                hor_grid='u', z_grid='1')
    vd(3) = var_desc("vbt_IC", "m s-1", &
                longname="Next initial condition for the barotropic meridional velocity",&
                hor_grid='v', z_grid='1')
    call register_restart_pair(CS%ubt_IC, CS%vbt_IC, vd(2), vd(3), .false., restart_CS)
  endif


  call register_restart_field(CS%dtbt, "DTBT", .false., restart_CS, &
                              longname="Barotropic timestep", units="seconds")

end subroutine register_barotropic_restarts

!> \namespace mom_barotropic
!!
!!  By Robert Hallberg, April 1994 - January 2007
!!
!!    This program contains the subroutines that time steps the
!!  linearized barotropic equations.  btstep is used to actually
!!  time step the barotropic equations, and contains most of the
!!  substance of this module.
!!
!!    btstep uses a forwards-backwards based scheme to time step
!!  the barotropic equations, returning the layers' accelerations due
!!  to the barotropic changes in the ocean state, the final free
!!  surface height (or column mass), and the volume (or mass) fluxes
!!  summed through the layers and averaged over the baroclinic time
!!  step.  As input, btstep takes the initial 3-D velocities, the
!!  inital free surface height, the 3-D accelerations of the layers,
!!  and the external forcing.  Everything in btstep is cast in terms
!!  of anomalies, so if everything is in balance, there is explicitly
!!  no acceleration due to btstep.
!!
!!    The spatial discretization of the continuity equation is second
!!  order accurate.  A flux conservative form is used to guarantee
!!  global conservation of volume.  The spatial discretization of the
!!  momentum equation is second order accurate.  The Coriolis force
!!  is written in a form which does not contribute to the energy
!!  tendency and which conserves linearized potential vorticity, f/D.
!!  These terms are exactly removed from the baroclinic momentum
!!  equations, so the linearization of vorticity advection will not
!!  degrade the overall solution.
!!
!!    btcalc calculates the fractional thickness of each layer at the
!!  velocity points, for later use in calculating the barotropic
!!  velocities and the averaged accelerations.  Harmonic mean
!!  thicknesses (i.e. 2*h_L*h_R/(h_L + h_R)) are used to avoid overly
!!  strong weighting of overly thin layers.  This may later be relaxed
!!  to use thicknesses determined from the continuity equations.
!!
!!    bt_mass_source determines the real mass sources for the
!!  barotropic solver, along with the corrective pseudo-fluxes that
!!  keep the barotropic and baroclinic estimates of the free surface
!!  height close to each other.  Given the layer thicknesses and the
!!  free surface height that correspond to each other, it calculates
!!  a corrective mass source that is added to the barotropic continuity*
!!  equation, and optionally adjusts a slowly varying correction rate.
!!  Newer algorithmic changes have deemphasized the need for this, but
!!  it is still here to add net water sources to the barotropic solver.*
!!
!!    barotropic_init allocates and initializes any barotropic arrays
!!  that have not been read from a restart file, reads parameters from
!!  the inputfile, and sets up diagnostic fields.
!!
!!    barotropic_end deallocates anything allocated in barotropic_init
!!  or register_barotropic_restarts.
!!
!!    register_barotropic_restarts is used to indicate any fields that
!!  are private to the barotropic solver that need to be included in
!!  the restart files, and to ensure that they are read.

end module MOM_barotropic
