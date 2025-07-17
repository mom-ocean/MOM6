!> Time step the adiabatic dynamic core of MOM using RK2 method.
module MOM_dynamics_split_RK2

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_variables,    only : vertvisc_type, thermo_var_ptrs, porous_barrier_type
use MOM_variables,    only : BT_cont_type, alloc_bt_cont_type, dealloc_bt_cont_type
use MOM_variables,    only : accel_diag_ptrs, ocean_internal_state, cont_diag_ptrs
use MOM_forcing_type, only : mech_forcing

use MOM_checksum_packages, only : MOM_thermo_chksum, MOM_state_chksum, MOM_accel_chksum
use MOM_cpu_clock,         only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,         only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT
use MOM_cpu_clock,         only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diabatic_driver,   only : diabatic_CS, extract_diabatic_member
use MOM_diag_mediator,     only : diag_mediator_init, enable_averages
use MOM_diag_mediator,     only : disable_averaging, post_data, safe_alloc_ptr
use MOM_diag_mediator,     only : post_product_u, post_product_sum_u
use MOM_diag_mediator,     only : post_product_v, post_product_sum_v
use MOM_diag_mediator,     only : register_diag_field, register_static_field
use MOM_diag_mediator,     only : set_diag_mediator_grid, diag_ctrl, diag_update_remap_grids
use MOM_domains,           only : To_South, To_West, To_All, CGRID_NE, SCALAR_PAIR
use MOM_domains,           only : To_North, To_East, Omit_Corners
use MOM_domains,           only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,           only : start_group_pass, complete_group_pass, pass_var, pass_vector
use MOM_debugging,         only : hchksum, uvchksum, query_debugging_checks
use MOM_error_handler,     only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_error_handler,     only : MOM_set_verbosity, callTree_showQuery
use MOM_error_handler,     only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,       only : get_param, log_version, param_file_type
use MOM_get_input,         only : directories
use MOM_io,                only : vardesc, var_desc, EAST_FACE, NORTH_FACE
use MOM_restart,           only : register_restart_field, register_restart_pair
use MOM_restart,           only : query_initialized, set_initialized, save_restart
use MOM_restart,           only : only_read_from_restarts
use MOM_restart,           only : restart_init, is_new_run, MOM_restart_CS
use MOM_time_manager,      only : time_type, time_type_to_real, operator(+)
use MOM_time_manager,      only : operator(-), operator(>), operator(*), operator(/)

use MOM_ALE,                   only : ALE_CS, ALE_remap_velocities
use MOM_barotropic,            only : barotropic_init, btstep, btcalc, bt_mass_source
use MOM_barotropic,            only : register_barotropic_restarts, set_dtbt, barotropic_CS
use MOM_barotropic,            only : barotropic_end
use MOM_boundary_update,       only : update_OBC_data, update_OBC_CS
use MOM_continuity,            only : continuity, continuity_CS
use MOM_continuity,            only : continuity_init, continuity_stencil
use MOM_CoriolisAdv,           only : CorAdCalc, CoriolisAdv_CS
use MOM_CoriolisAdv,           only : CoriolisAdv_init, CoriolisAdv_end
use MOM_CVMix_KPP,             only : KPP_get_BLD, KPP_CS
use MOM_debugging,             only : check_redundant
use MOM_energetic_PBL,         only : energetic_PBL_get_MLD, energetic_PBL_CS
use MOM_grid,                  only : ocean_grid_type
use MOM_harmonic_analysis,     only : harmonic_analysis_CS
use MOM_hor_index,             only : hor_index_type
use MOM_hor_visc,              only : horizontal_viscosity, hor_visc_CS, hor_visc_vel_stencil
use MOM_hor_visc,              only : hor_visc_init, hor_visc_end
use MOM_interface_heights,     only : thickness_to_dz, find_col_avg_SpV
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use MOM_MEKE_types,            only : MEKE_type
use MOM_open_boundary,         only : ocean_OBC_type, radiation_open_bdry_conds
use MOM_open_boundary,         only : open_boundary_zero_normal_flow, open_boundary_query
use MOM_open_boundary,         only : open_boundary_test_extern_h, update_OBC_ramp
use MOM_PressureForce,         only : PressureForce, PressureForce_CS
use MOM_PressureForce,         only : PressureForce_init
use MOM_set_visc,              only : set_viscous_ML, set_visc_CS
use MOM_stochastics,           only : stochastic_CS
use MOM_thickness_diffuse,     only : thickness_diffuse_CS
use MOM_self_attr_load,        only : SAL_CS
use MOM_self_attr_load,        only : SAL_init, SAL_end
use MOM_tidal_forcing,         only : tidal_forcing_CS
use MOM_tidal_forcing,         only : tidal_forcing_init, tidal_forcing_end
use MOM_unit_scaling,          only : unit_scale_type
use MOM_vert_friction,         only : vertvisc, vertvisc_coef, vertvisc_remnant
use MOM_vert_friction,         only : vertvisc_init, vertvisc_end, vertvisc_CS
use MOM_vert_friction,         only : updateCFLtruncationValue, vertFPmix
use MOM_verticalGrid,          only : verticalGrid_type, get_thickness_units
use MOM_verticalGrid,          only : get_flux_units, get_tr_flux_units
use MOM_wave_interface,        only : wave_parameters_CS, Stokes_PGF

implicit none ; private

#include <MOM_memory.h>

!> MOM_dynamics_split_RK2 module control structure
type, public :: MOM_dyn_split_RK2_CS ; private
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
    CAu, &    !< CAu = f*v - u.grad(u) [L T-2 ~> m s-2]
    CAu_pred, & !< The predictor step value of CAu = f*v - u.grad(u) [L T-2 ~> m s-2]
    PFu, &    !< PFu = -dM/dx [L T-2 ~> m s-2]
    PFu_Stokes, & !< PFu_Stokes = -d/dx int_r (u_L*duS/dr) [L T-2 ~> m s-2]
    diffu     !< Zonal acceleration due to convergence of the along-isopycnal stress tensor [L T-2 ~> m s-2]

  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
    CAv, &    !< CAv = -f*u - u.grad(v) [L T-2 ~> m s-2]
    CAv_pred, & !< The predictor step value of CAv = -f*u - u.grad(v) [L T-2 ~> m s-2]
    PFv, &    !< PFv = -dM/dy [L T-2 ~> m s-2]
    PFv_Stokes, & !< PFv_Stokes = -d/dy int_r (v_L*dvS/dr) [L T-2 ~> m s-2]
    diffv     !< Meridional acceleration due to convergence of the along-isopycnal stress tensor [L T-2 ~> m s-2]

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: visc_rem_u
              !< Both the fraction of the zonal momentum originally in a
              !! layer that remains after a time-step of viscosity, and the
              !! fraction of a time-step worth of a barotropic acceleration
              !! that a layer experiences after viscosity is applied [nondim].
              !! Nondimensional between 0 (at the bottom) and 1 (far above).
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: u_accel_bt
              !< The zonal layer accelerations due to the difference between
              !! the barotropic accelerations and the baroclinic accelerations
              !! that were fed into the barotopic calculation [L T-2 ~> m s-2]
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: visc_rem_v
              !< Both the fraction of the meridional momentum originally in
              !! a layer that remains after a time-step of viscosity, and the
              !! fraction of a time-step worth of a barotropic acceleration
              !! that a layer experiences after viscosity is applied [nondim].
              !! Nondimensional between 0 (at the bottom) and 1 (far above).
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: v_accel_bt
              !< The meridional layer accelerations due to the difference between
              !! the barotropic accelerations and the baroclinic accelerations
              !! that were fed into the barotopic calculation [L T-2 ~> m s-2]

  ! The following variables are only used with the split time stepping scheme.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_)             :: eta    !< Instantaneous free surface height (in Boussinesq
                                                                  !! mode) or column mass anomaly (in non-Boussinesq
                                                                  !! mode) [H ~> m or kg m-2]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: u_av   !< layer x-velocity with vertical mean replaced by
                                                                  !! time-mean barotropic velocity over a baroclinic
                                                                  !! timestep [L T-1 ~> m s-1]
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: v_av   !< layer y-velocity with vertical mean replaced by
                                                                  !! time-mean barotropic velocity over a baroclinic
                                                                  !! timestep [L T-1 ~> m s-1]
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_)      :: h_av   !< arithmetic mean of two successive layer
                                                                  !! thicknesses [H ~> m or kg m-2]
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_)             :: eta_PF !< instantaneous SSH used in calculating PFu and
                                                                  !! PFv [H ~> m or kg m-2]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_)        :: uhbt   !< average x-volume or mass flux determined by the
                                                                  !! barotropic solver [H L2 T-1 ~> m3 s-1 or kg s-1].
                                                                  !! uhbt is roughly equal to the vertical sum of uh.
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_)        :: vhbt   !< average y-volume or mass flux determined by the
                                                                  !! barotropic solver [H L2 T-1 ~> m3 s-1 or kg s-1].
                                                                  !! vhbt is roughly equal to vertical sum of vh.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_)      :: pbce   !< pbce times eta gives the baroclinic pressure
                                                                  !! anomaly in each layer due to free surface height
                                                                  !! anomalies [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2].
  type(KPP_CS),           pointer :: KPP_CSp => NULL() !< KPP control structure needed to ge
  type(energetic_PBL_CS), pointer :: energetic_PBL_CSp => NULL()  !< ePBL control structure

  real, pointer, dimension(:,:) :: taux_bot => NULL()  !< frictional x-bottom stress from the ocean
                                                       !! to the seafloor [R L Z T-2 ~> Pa]
  real, pointer, dimension(:,:) :: tauy_bot => NULL()  !< frictional y-bottom stress from the ocean
                                                       !! to the seafloor [R L Z T-2 ~> Pa]
  type(BT_cont_type), pointer   :: BT_cont  => NULL()  !<  A structure with elements that describe the
                                                       !! effective summed open face areas as a function
                                                       !! of barotropic flow.

  ! This is to allow the previous, velocity-based coupling with between the
  ! baroclinic and barotropic modes.
  logical :: BT_use_layer_fluxes  !< If true, use the summed layered fluxes plus
                                  !! an adjustment due to a changed barotropic
                                  !! velocity in the barotropic continuity equation.
  logical :: split_bottom_stress  !< If true, provide the bottom stress
                                  !! calculated by the vertical viscosity to the
                                  !! barotropic solver.
  logical :: dtbt_use_bt_cont     !< If true, use BT_cont to calculate DTBT.
  logical :: store_CAu            !< If true, store the Coriolis and advective accelerations at the
                                  !! end of the timestep for use in the next predictor step.
  logical :: CAu_pred_stored      !< If true, the Coriolis and advective accelerations at the
                                  !! end of the timestep have been stored for use in the next
                                  !! predictor step.  This is used to accomodate various generations
                                  !! of restart files.
  logical :: calculate_SAL        !< If true, calculate self-attraction and loading.
  logical :: use_tides            !< If true, tidal forcing is enabled.
  logical :: remap_aux            !< If true, apply ALE remapping to all of the auxiliary 3-D
                                  !! variables that are needed to reproduce across restarts,
                                  !! similarly to what is done with the primary state variables.

  real    :: be      !< A nondimensional number from 0.5 to 1 that controls
                     !! the backward weighting of the time stepping scheme [nondim]
  real    :: begw    !< A nondimensional number from 0 to 1 that controls
                     !! the extent to which the treatment of gravity waves
                     !! is forward-backward (0) or simulated backward
                     !! Euler (1) [nondim].  0 is often used.
  real    :: Cemp_NL   !< Empirical coefficient of non-local momentum mixing [nondim]
  logical :: debug     !< If true, write verbose checksums for debugging purposes.
  logical :: debug_OBC !< If true, do additional calls resetting values to help debug the correctness
                       !! of the open boundary condition code.
  logical :: fpmix     !< If true, add non-local momentum flux increments and diffuse down the Eulerian gradient.
  logical :: module_is_initialized = .false. !< Record whether this module has been initialized.
  logical :: visc_rem_dt_bug = .true. !< If true, recover a bug that uses dt_pred rather than dt for vertvisc_rem
                                      !! at the end of predictor.

  !>@{ Diagnostic IDs
  integer :: id_uh     = -1, id_vh     = -1
  integer :: id_umo    = -1, id_vmo    = -1
  integer :: id_umo_2d = -1, id_vmo_2d = -1
  integer :: id_PFu    = -1, id_PFv    = -1
  integer :: id_CAu    = -1, id_CAv    = -1
  integer :: id_ueffA = -1, id_veffA = -1
  ! integer :: id_hf_PFu    = -1, id_hf_PFv    = -1
  integer :: id_h_PFu    = -1, id_h_PFv    = -1
  integer :: id_hf_PFu_2d = -1, id_hf_PFv_2d = -1
  integer :: id_intz_PFu_2d = -1, id_intz_PFv_2d = -1
  integer :: id_PFu_visc_rem = -1, id_PFv_visc_rem = -1
  ! integer :: id_hf_CAu    = -1, id_hf_CAv    = -1
  integer :: id_h_CAu    = -1, id_h_CAv    = -1
  integer :: id_hf_CAu_2d = -1, id_hf_CAv_2d = -1
  integer :: id_intz_CAu_2d = -1, id_intz_CAv_2d = -1
  integer :: id_CAu_visc_rem = -1, id_CAv_visc_rem = -1
  integer :: id_deta_dt = -1

  ! Split scheme only.
  integer :: id_uav        = -1, id_vav        = -1
  integer :: id_u_BT_accel = -1, id_v_BT_accel = -1
  ! integer :: id_hf_u_BT_accel    = -1, id_hf_v_BT_accel    = -1
  integer :: id_h_u_BT_accel    = -1, id_h_v_BT_accel    = -1
  integer :: id_hf_u_BT_accel_2d = -1, id_hf_v_BT_accel_2d = -1
  integer :: id_intz_u_BT_accel_2d = -1, id_intz_v_BT_accel_2d = -1
  integer :: id_u_BT_accel_visc_rem    = -1, id_v_BT_accel_visc_rem    = -1
  !>@}

  type(diag_ctrl), pointer       :: diag => NULL() !< A structure that is used to regulate the
                                         !! timing of diagnostic output.
  type(accel_diag_ptrs), pointer :: ADp => NULL()  !< A structure pointing to the various
                                         !! accelerations in the momentum equations,
                                         !! which can later be used to calculate
                                         !! derived diagnostics like energy budgets.
  type(accel_diag_ptrs), pointer :: AD_pred => NULL() !< A structure pointing to the various
                                         !! predictor step accelerations in the momentum equations,
                                         !! which can be used to debug truncations.
  type(cont_diag_ptrs), pointer  :: CDp => NULL()  !< A structure with pointers to various
                                         !! terms in the continuity equations,
                                         !! which can later be used to calculate
                                         !! derived diagnostics like energy budgets.

  ! The remainder of the structure points to child subroutines' control structures.
  !> A pointer to the horizontal viscosity control structure
  type(hor_visc_CS) :: hor_visc
  !> A pointer to the continuity control structure
  type(continuity_CS) :: continuity_CSp
  !> The CoriolisAdv control structure
  type(CoriolisAdv_CS) :: CoriolisAdv
  !> A pointer to the PressureForce control structure
  type(PressureForce_CS) :: PressureForce_CSp
  !> A pointer to a structure containing interface height diffusivities
  type(vertvisc_CS),      pointer :: vertvisc_CSp      => NULL()
  !> A pointer to the set_visc control structure
  type(set_visc_CS),      pointer :: set_visc_CSp      => NULL()
  !> A pointer to the barotropic stepping control structure
  type(barotropic_CS) :: barotropic_CSp
  !> A pointer to the SAL control structure
  type(SAL_CS) :: SAL_CSp
  !> A pointer to the tidal forcing control structure
  type(tidal_forcing_CS) :: tides_CSp
  !> A pointer to the harmonic analysis control structure
  type(harmonic_analysis_CS) :: HA_CSp
  !> A pointer to the ALE control structure.
  type(ALE_CS), pointer :: ALE_CSp => NULL()

  type(ocean_OBC_type),   pointer :: OBC => NULL() !< A pointer to an open boundary
     !! condition type that specifies whether, where, and  what open boundary
     !! conditions are used.  If no open BCs are used, this pointer stays
     !! nullified.  Flather OBCs use open boundary_CS as well.
  !> A pointer to the update_OBC control structure
  type(update_OBC_CS),    pointer :: update_OBC_CSp => NULL()

  type(group_pass_type) :: pass_eta  !< Structure for group halo pass
  type(group_pass_type) :: pass_visc_rem  !< Structure for group halo pass
  type(group_pass_type) :: pass_uvp  !< Structure for group halo pass
  type(group_pass_type) :: pass_hp_uv  !< Structure for group halo pass
  type(group_pass_type) :: pass_uv  !< Structure for group halo pass
  type(group_pass_type) :: pass_h  !< Structure for group halo pass
  type(group_pass_type) :: pass_av_uvh  !< Structure for group halo pass

end type MOM_dyn_split_RK2_CS


public step_MOM_dyn_split_RK2
public register_restarts_dyn_split_RK2
public initialize_dyn_split_RK2
public remap_dyn_split_RK2_aux_vars
public init_dyn_split_RK2_diabatic
public end_dyn_split_RK2

!>@{ CPU time clock IDs
integer :: id_clock_Cor, id_clock_pres, id_clock_vertvisc
integer :: id_clock_horvisc, id_clock_mom_update
integer :: id_clock_continuity, id_clock_thick_diff
integer :: id_clock_btstep, id_clock_btcalc, id_clock_btforce
integer :: id_clock_pass, id_clock_pass_init
!>@}

contains

!> RK2 splitting for time stepping MOM adiabatic dynamics
subroutine step_MOM_dyn_split_RK2(u_inst, v_inst, h, tv, visc, Time_local, dt, forces, &
                                  p_surf_begin, p_surf_end, uh, vh, uhtr, vhtr, eta_av, G, GV, US, CS, &
                                  calc_dtbt, VarMix, MEKE, thickness_diffuse_CSp, pbv, STOCH, Waves)
  type(ocean_grid_type),             intent(inout) :: G            !< Ocean grid structure
  type(verticalGrid_type),           intent(in)    :: GV           !< Ocean vertical grid structure
  type(unit_scale_type),             intent(in)    :: US           !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                             target, intent(inout) :: u_inst       !< Instantaneous zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                             target, intent(inout) :: v_inst       !< Instantaneous meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                     intent(inout) :: h            !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),             intent(in)    :: tv           !< Thermodynamic type
  type(vertvisc_type),               intent(inout) :: visc         !< Vertical visc, bottom drag, and related
  type(time_type),                   intent(in)    :: Time_local   !< Model time at end of time step
  real,                              intent(in)    :: dt           !< Baroclinic dynamics time step [T ~> s]
  type(mech_forcing),                intent(in)    :: forces       !< A structure with the driving mechanical forces
  real, dimension(:,:),              pointer       :: p_surf_begin !< Surface pressure at the start of this dynamic
                                                                   !! time step [R L2 T-2 ~> Pa]
  real, dimension(:,:),              pointer       :: p_surf_end   !< Surface pressure at the end of this dynamic
                                                                   !! time step [R L2 T-2 ~> Pa]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                             target, intent(inout) :: uh           !< Zonal volume or mass transport
                                                                   !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                             target, intent(inout) :: vh           !< Meridional volume or mass transport
                                                                   !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                     intent(inout) :: uhtr         !< Accumulated zonal volume or mass transport
                                                                   !! since last tracer advection [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                     intent(inout) :: vhtr         !< Accumulated meridional volume or mass transport
                                                                   !! since last tracer advection [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G)),  intent(out)   :: eta_av       !< Free surface height or column mass
                                                                   !! averaged over time step [H ~> m or kg m-2]
  type(MOM_dyn_split_RK2_CS),        pointer       :: CS           !< Module control structure
  logical,                           intent(in)    :: calc_dtbt    !< If true, recalculate the barotropic time step
  type(VarMix_CS),                   intent(inout) :: VarMix       !< Variable mixing control structure
  type(MEKE_type),                   intent(inout) :: MEKE         !< MEKE fields
  type(thickness_diffuse_CS),        intent(inout) :: thickness_diffuse_CSp !< Pointer to a structure containing
                                                                   !! interface height diffusivities
  type(porous_barrier_type),         intent(in)    :: pbv          !< porous barrier fractional cell metrics
  type(stochastic_CS),     optional, intent(inout) :: STOCH        !< Stochastic control structure
  type(wave_parameters_CS), optional, pointer      :: Waves        !< A pointer to a structure containing
                                                                   !! fields related to the surface wave conditions

  ! local variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: up  ! Predicted zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: vp  ! Predicted meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))  :: hp  ! Predicted thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))  :: dz  ! Distance between the interfaces around a layer [Z ~> m]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: ueffA   ! Effective Area of U-Faces [H L ~> m2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: veffA   ! Effective Area of V-Faces [H L ~> m2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: u_bc_accel ! The summed zonal baroclinic accelerations
                                                        ! of each layer calculated by the non-barotropic
                                                        ! part of the model [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: v_bc_accel ! The summed meridional baroclinic accelerations
                                                        ! of each layer calculated by the non-barotropic
                                                        ! part of the model [L T-2 ~> m s-2]

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), target :: uh_in ! The zonal mass transports that would be
                                ! obtained using the initial velocities [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), target :: vh_in ! The meridional mass transports that would be
                                ! obtained using the initial velocities [H L2 T-1 ~> m3 s-1 or kg s-1]

  real, dimension(SZI_(G),SZJ_(G)) :: eta_pred ! The predictor value of the free surface height
                                               ! or column mass [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)) :: SpV_avg  ! The column averaged specific volume [R-1 ~> m3 kg-1]
  real, dimension(SZI_(G),SZJ_(G)) :: deta_dt  ! A diagnostic of the time derivative of the free surface
                                               ! height or column mass [H T-1 ~> m s-1 or kg m-2 s-1]

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: u_old_rad_OBC ! The starting zonal velocities, which are
                                ! saved for use in the radiation open boundary condition code [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: v_old_rad_OBC ! The starting meridional velocities, which are
                                ! saved for use in the radiation open boundary condition code [L T-1 ~> m s-1]

  ! GMM, TODO: make these allocatable?
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: uold ! u-velocity before vert_visc is applied, for fpmix
                                                     !                                      [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: vold ! v-velocity before vert_visc is applied, for fpmix
                                                     !                                      [L T-1 ~> m s-1]
  real :: pres_to_eta ! A factor that converts pressures to the units of eta
                      ! [H T2 R-1 L-2 ~> m Pa-1 or kg m-2 Pa-1]
  real, pointer, dimension(:,:) :: &
    p_surf => NULL(), &         ! A pointer to the surface pressure [R L2 T-2 ~> Pa]
    eta_PF_start => NULL(), &   ! The value of eta that corresponds to the starting pressure
                                ! for the barotropic solver [H ~> m or kg m-2]
    taux_bot => NULL(), &       ! A pointer to the zonal bottom stress in some cases [R L Z T-2 ~> Pa]
    tauy_bot => NULL(), &       ! A pointer to the meridional bottom stress in some cases [R L Z T-2 ~> Pa]
    ! This pointer is just used as shorthand for CS%eta.
    eta => NULL()               ! A pointer to the instantaneous free surface height (in Boussinesq
                                ! mode) or column mass anomaly (in non-Boussinesq mode) [H ~> m or kg m-2]

  real, pointer, dimension(:,:,:) :: &
    ! These pointers are used to alter which fields are passed to btstep with various options:
    u_ptr => NULL(), &   ! A pointer to a zonal velocity [L T-1 ~> m s-1]
    v_ptr => NULL(), &   ! A pointer to a meridional velocity [L T-1 ~> m s-1]
    uh_ptr => NULL(), &  ! A pointer to a zonal volume or mass transport [H L2 T-1 ~> m3 s-1 or kg s-1]
    vh_ptr => NULL(), &  ! A pointer to a meridional volume or mass transport [H L2 T-1 ~> m3 s-1 or kg s-1]
    ! These pointers are just used as shorthand for CS%u_av, CS%v_av, and CS%h_av.
    u_av, & ! The zonal velocity time-averaged over a time step [L T-1 ~> m s-1].
    v_av, & ! The meridional velocity time-averaged over a time step [L T-1 ~> m s-1].
    h_av    ! The layer thickness time-averaged over a time step [H ~> m or kg m-2].

  real, dimension(SZI_(G),SZJ_(G)) :: hbl       ! Boundary layer depth from Cvmix [H ~> m or kg m-2]
  real :: dt_pred   ! The time step for the predictor part of the baroclinic time stepping [T ~> s].
  real :: Idt_bc    ! Inverse of the baroclinic timestep [T-1 ~> s-1]
  logical :: dyn_p_surf
  logical :: debug_redundant  ! If true, check redundant values on PE boundaries when debugging
  logical :: BT_cont_BT_thick ! If true, use the BT_cont_type to estimate the
                              ! relative weightings of the layers in calculating
                              ! the barotropic accelerations.
  logical :: Use_Stokes_PGF ! If true, add Stokes PGF to hydrostatic PGF
  !---For group halo pass
  logical :: showCallTree, sym
  logical :: lFPpost        ! Used to only post diagnostics in vertFPmix when fpmix=true and
                            ! in the  corrector step (not the predict)
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: cont_stencil, obc_stencil, vel_stencil

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  u_av => CS%u_av ; v_av => CS%v_av ; h_av => CS%h_av ; eta => CS%eta

  Idt_bc = 1.0 / dt

  sym = G%Domain%symmetric  ! switch to include symmetric domain in checksums

  showCallTree = callTree_showQuery()
  if (showCallTree) call callTree_enter("step_MOM_dyn_split_RK2(), MOM_dynamics_split_RK2.F90")

  !$OMP parallel do default(shared)
  do k=1,nz
    do j=G%jsd,G%jed   ; do i=G%isdB,G%iedB ;  up(i,j,k) = 0.0 ; enddo ; enddo
    do j=G%jsdB,G%jedB ; do i=G%isd,G%ied   ;  vp(i,j,k) = 0.0 ; enddo ; enddo
    do j=G%jsd,G%jed   ; do i=G%isd,G%ied   ;  hp(i,j,k) = h(i,j,k) ; enddo ; enddo
  enddo

  ! Update CFL truncation value as function of time
  call updateCFLtruncationValue(Time_local, CS%vertvisc_CSp, US)

  if (CS%debug) then
    call query_debugging_checks(do_redundant=debug_redundant)
    call MOM_state_chksum("Start predictor ", u_inst, v_inst, h, uh, vh, G, GV, US, symmetric=sym)
    if (debug_redundant) then
      call check_redundant("Start predictor u ", u_inst, v_inst, G, unscale=US%L_T_to_m_s)
      call check_redundant("Start predictor uh ", uh, vh, G, unscale=GV%H_to_MKS*US%L_to_m**2*US%s_to_T)
    endif
  endif

  dyn_p_surf = associated(p_surf_begin) .and. associated(p_surf_end)
  if (dyn_p_surf) then
    p_surf => p_surf_end
    call safe_alloc_ptr(eta_PF_start,G%isd,G%ied,G%jsd,G%jed)
    eta_PF_start(:,:) = 0.0
  else
    p_surf => forces%p_surf
  endif

  if (associated(CS%OBC)) then
    if (CS%debug_OBC) call open_boundary_test_extern_h(G, GV, CS%OBC, h)

    ! Update OBC ramp value as function of time
    call update_OBC_ramp(Time_local, CS%OBC, US)

    do k=1,nz ; do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB
      u_old_rad_OBC(I,j,k) = u_av(I,j,k)
    enddo ; enddo ; enddo
    do k=1,nz ; do J=G%JsdB,G%JedB ; do i=G%isd,G%ied
      v_old_rad_OBC(i,J,k) = v_av(i,J,k)
    enddo ; enddo ; enddo
  endif

  BT_cont_BT_thick = .false.
  if (associated(CS%BT_cont)) BT_cont_BT_thick = &
    (allocated(CS%BT_cont%h_u) .and. allocated(CS%BT_cont%h_v))

  if (CS%split_bottom_stress) then
    taux_bot => CS%taux_bot ; tauy_bot => CS%tauy_bot
  endif

  !--- begin set up for group halo pass

  cont_stencil = continuity_stencil(CS%continuity_CSp)
  obc_stencil = 2
  if (associated(CS%OBC)) then
    if (CS%OBC%oblique_BCs_exist_globally) obc_stencil = 3
  endif
  vel_stencil = max(2, obc_stencil, hor_visc_vel_stencil(CS%hor_visc))
  call cpu_clock_begin(id_clock_pass)
  call create_group_pass(CS%pass_eta, eta, G%Domain, halo=1)
  call create_group_pass(CS%pass_visc_rem, CS%visc_rem_u, CS%visc_rem_v, G%Domain, &
                         To_All+SCALAR_PAIR, CGRID_NE, halo=max(1,cont_stencil))
  call create_group_pass(CS%pass_uvp, up, vp, G%Domain, halo=max(1,cont_stencil))
  call create_group_pass(CS%pass_hp_uv, hp, G%Domain, halo=2)
  call create_group_pass(CS%pass_hp_uv, u_av, v_av, G%Domain, halo=vel_stencil)
  call create_group_pass(CS%pass_hp_uv, uh(:,:,:), vh(:,:,:), G%Domain, halo=vel_stencil)

  call create_group_pass(CS%pass_uv, u_inst, v_inst, G%Domain, halo=max(2,cont_stencil))
  call create_group_pass(CS%pass_h, h, G%Domain, halo=max(2,cont_stencil))
  call create_group_pass(CS%pass_av_uvh, u_av, v_av, G%Domain, halo=vel_stencil)
  call create_group_pass(CS%pass_av_uvh, uh(:,:,:), vh(:,:,:), G%Domain, halo=vel_stencil)
  call cpu_clock_end(id_clock_pass)
  !--- end set up for group halo pass

! PFu = d/dx M(h,T,S)
! pbce = dM/deta
  if (CS%begw == 0.0) call enable_averages(dt, Time_local, CS%diag)
  call cpu_clock_begin(id_clock_pres)
  call PressureForce(h, tv, CS%PFu, CS%PFv, G, GV, US, CS%PressureForce_CSp, &
                     CS%ALE_CSp, CS%ADp, p_surf, CS%pbce, CS%eta_PF)
  if (dyn_p_surf) then
    pres_to_eta = 1.0 / (GV%g_Earth * GV%H_to_RZ)
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta_PF_start(i,j) = CS%eta_PF(i,j) - pres_to_eta * (p_surf_begin(i,j) - p_surf_end(i,j))
    enddo ; enddo
  endif
  ! Stokes shear force contribution to pressure gradient
  Use_Stokes_PGF = present(Waves)
  if (Use_Stokes_PGF) then
    Use_Stokes_PGF = associated(Waves)
    if (Use_Stokes_PGF) Use_Stokes_PGF = Waves%Stokes_PGF
    if (Use_Stokes_PGF) then
      call thickness_to_dz(h, tv, dz, G, GV, US, halo_size=1)
      call Stokes_PGF(G, GV, US, dz, u_inst, v_inst, CS%PFu_Stokes, CS%PFv_Stokes, Waves)

      ! We are adding Stokes_PGF to hydrostatic PGF here.  The diag PFu/PFv
      ! will therefore report the sum total PGF and we avoid other
      ! modifications in the code.  The PFu_Stokes is output within the waves routines.
      if (.not.Waves%Passive_Stokes_PGF) then
        do k=1,nz
          do j=js,je ; do I=Isq,Ieq
            CS%PFu(I,j,k) = CS%PFu(I,j,k) + CS%PFu_Stokes(I,j,k)
          enddo ; enddo
        enddo
        do k=1,nz
          do J=Jsq,Jeq ; do i=is,ie
             CS%PFv(i,J,k) = CS%PFv(i,J,k) + CS%PFv_Stokes(i,J,k)
          enddo ; enddo
        enddo
      endif
    endif
  endif
  call cpu_clock_end(id_clock_pres)
  call disable_averaging(CS%diag)
  if (showCallTree) call callTree_wayPoint("done with PressureForce (step_MOM_dyn_split_RK2)")

  if (associated(CS%OBC)) then ; if (CS%OBC%update_OBC) then
    call update_OBC_data(CS%OBC, G, GV, US, tv, h, CS%update_OBC_CSp, Time_local)
  endif ; endif
  if (associated(CS%OBC) .and. CS%debug_OBC) &
    call open_boundary_zero_normal_flow(CS%OBC, G, GV, CS%PFu, CS%PFv)

  if (G%nonblocking_updates) &
    call start_group_pass(CS%pass_eta, G%Domain, clock=id_clock_pass)

! CAu = -(f+zeta_av)/h_av vh + d/dx KE_av
  if (.not.CS%CAu_pred_stored) then
    ! Calculate a predictor-step estimate of the Coriolis and momentum advection terms,
    ! if it was not already stored from the end of the previous time step.
    call cpu_clock_begin(id_clock_Cor)
    call CorAdCalc(u_av, v_av, h_av, uh, vh, CS%CAu_pred, CS%CAv_pred, CS%OBC, CS%AD_pred, &
                   G, GV, US, CS%CoriolisAdv, pbv, Waves=Waves)
    call cpu_clock_end(id_clock_Cor)
    if (showCallTree) call callTree_wayPoint("done with CorAdCalc (step_MOM_dyn_split_RK2)")
  endif

! u_bc_accel = CAu + PFu + diffu(u[n-1])
  call cpu_clock_begin(id_clock_btforce)
  !$OMP parallel do default(shared)
  do k=1,nz
    do j=js,je ; do I=Isq,Ieq
      u_bc_accel(I,j,k) = (CS%CAu_pred(I,j,k) + CS%PFu(I,j,k)) + CS%diffu(I,j,k)
    enddo ; enddo
    do J=Jsq,Jeq ; do i=is,ie
      v_bc_accel(i,J,k) = (CS%CAv_pred(i,J,k) + CS%PFv(i,J,k)) + CS%diffv(i,J,k)
    enddo ; enddo
  enddo
  if (associated(CS%OBC)) then
    call open_boundary_zero_normal_flow(CS%OBC, G, GV, u_bc_accel, v_bc_accel)
  endif
  call cpu_clock_end(id_clock_btforce)

  if (CS%debug) then
    call MOM_accel_chksum("pre-btstep accel", CS%CAu_pred, CS%CAv_pred, CS%PFu, CS%PFv, &
                          CS%diffu, CS%diffv, G, GV, US, CS%pbce, u_bc_accel, v_bc_accel, &
                          symmetric=sym)
    if (debug_redundant) then
      call check_redundant("pre-btstep CS%CA ", CS%CAu_pred, CS%CAv_pred, G, unscale=US%L_T2_to_m_s2)
      call check_redundant("pre-btstep CS%PF ", CS%PFu, CS%PFv, G, unscale=US%L_T2_to_m_s2)
      call check_redundant("pre-btstep CS%diff ", CS%diffu, CS%diffv, G, unscale=US%L_T2_to_m_s2)
      call check_redundant("pre-btstep u_bc_accel ", u_bc_accel, v_bc_accel, G, unscale=US%L_T2_to_m_s2)
    endif
  endif

  call cpu_clock_begin(id_clock_vertvisc)
  !$OMP parallel do default(shared)
  do k=1,nz
    do j=js,je ; do I=Isq,Ieq
      up(I,j,k) = G%mask2dCu(I,j) * (u_inst(I,j,k) + dt * u_bc_accel(I,j,k))
    enddo ; enddo
    do J=Jsq,Jeq ; do i=is,ie
      vp(i,J,k) = G%mask2dCv(i,J) * (v_inst(i,J,k) + dt * v_bc_accel(i,J,k))
    enddo ; enddo
  enddo

  call enable_averages(dt, Time_local, CS%diag)
  call set_viscous_ML(u_inst, v_inst, h, tv, forces, visc, dt, G, GV, US, CS%set_visc_CSp)
  call disable_averaging(CS%diag)

  if (CS%debug) then
    call uvchksum("before vertvisc: up", up, vp, G%HI, haloshift=0, symmetric=sym, unscale=US%L_T_to_m_s)
  endif
  call thickness_to_dz(h, tv, dz, G, GV, US, halo_size=1)
  call vertvisc_coef(up, vp, h, dz, forces, visc, tv, dt, G, GV, US, CS%vertvisc_CSp, CS%OBC, VarMix)
  call vertvisc_remnant(visc, CS%visc_rem_u, CS%visc_rem_v, dt, G, GV, US, CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)
  if (showCallTree) call callTree_wayPoint("done with vertvisc_coef (step_MOM_dyn_split_RK2)")


  call cpu_clock_begin(id_clock_pass)
  if (G%nonblocking_updates) then
    call complete_group_pass(CS%pass_eta, G%Domain)
    call start_group_pass(CS%pass_visc_rem, G%Domain)
  else
    call do_group_pass(CS%pass_eta, G%Domain)
    call do_group_pass(CS%pass_visc_rem, G%Domain)
  endif
  call cpu_clock_end(id_clock_pass)

  call cpu_clock_begin(id_clock_btcalc)
  ! Calculate the relative layer weights for determining barotropic quantities.
  if (.not.BT_cont_BT_thick) &
    call btcalc(h, G, GV, CS%barotropic_CSp, OBC=CS%OBC)
  call bt_mass_source(h, eta, .true., G, GV, CS%barotropic_CSp)

  SpV_avg(:,:) = 0.0
  if ((.not.GV%Boussinesq) .and. associated(CS%OBC)) then
    ! Determine the column average specific volume if it is needed due to the
    ! use of Flather open boundary conditions in non-Boussinesq mode.
    if (open_boundary_query(CS%OBC, apply_Flather_OBC=.true.)) &
      call find_col_avg_SpV(h, SpV_avg, tv, G, GV, US)
  endif
  call cpu_clock_end(id_clock_btcalc)

  if (G%nonblocking_updates) &
    call complete_group_pass(CS%pass_visc_rem, G%Domain, clock=id_clock_pass)

! u_accel_bt = layer accelerations due to barotropic solver
  if (associated(CS%BT_cont) .or. CS%BT_use_layer_fluxes) then
    call cpu_clock_begin(id_clock_continuity)
    call continuity(u_inst, v_inst, h, hp, uh_in, vh_in, dt, G, GV, US, CS%continuity_CSp, CS%OBC, pbv, &
                    visc_rem_u=CS%visc_rem_u, visc_rem_v=CS%visc_rem_v, BT_cont=CS%BT_cont)
    call cpu_clock_end(id_clock_continuity)
    if (BT_cont_BT_thick) then
      call btcalc(h, G, GV, CS%barotropic_CSp, CS%BT_cont%h_u, CS%BT_cont%h_v, &
                  OBC=CS%OBC)
    endif
    if (showCallTree) call callTree_wayPoint("done with continuity[BT_cont] (step_MOM_dyn_split_RK2)")
  endif

  if (CS%BT_use_layer_fluxes) then
    uh_ptr => uh_in ; vh_ptr => vh_in; u_ptr => u_inst ; v_ptr => v_inst
  endif

  call cpu_clock_begin(id_clock_btstep)
  if (calc_dtbt) then
    if (CS%dtbt_use_bt_cont .and. associated(CS%BT_cont)) then
      call set_dtbt(G, GV, US, CS%barotropic_CSp, CS%pbce, BT_cont=CS%BT_cont)
    else
      ! In the following call, eta is only used when NONLINEAR_BT_CONTINUITY is True. Otherwise, dtbt is effectively
      ! calculated with eta=0. Note that NONLINEAR_BT_CONTINUITY is False if BT_CONT is used, which is the default.
      call set_dtbt(G, GV, US, CS%barotropic_CSp, CS%pbce, eta=eta)
    endif
  endif
  if (showCallTree) call callTree_enter("btstep(), MOM_barotropic.F90")
  ! This is the predictor step call to btstep.
  ! The CS%ADp argument here stores the weights for certain integrated diagnostics.
  call btstep(u_inst, v_inst, eta, dt, u_bc_accel, v_bc_accel, forces, CS%pbce, CS%eta_PF, u_av, v_av, &
              CS%u_accel_bt, CS%v_accel_bt, eta_pred, CS%uhbt, CS%vhbt, G, GV, US, &
              CS%barotropic_CSp, CS%visc_rem_u, CS%visc_rem_v, SpV_avg, CS%ADp, CS%OBC, CS%BT_cont, &
              eta_PF_start, taux_bot, tauy_bot, uh_ptr, vh_ptr, u_ptr, v_ptr)
  if (showCallTree) call callTree_leave("btstep()")
  call cpu_clock_end(id_clock_btstep)

! up = u + dt_pred*( u_bc_accel + u_accel_bt )
  dt_pred = dt * CS%be
  call cpu_clock_begin(id_clock_mom_update)

  !$OMP parallel do default(shared)
  do k=1,nz
    do J=Jsq,Jeq ; do i=is,ie
      vp(i,J,k) = G%mask2dCv(i,J) * (v_inst(i,J,k) + dt_pred * &
                      (v_bc_accel(i,J,k) + CS%v_accel_bt(i,J,k)))
    enddo ; enddo
    do j=js,je ; do I=Isq,Ieq
      up(I,j,k) = G%mask2dCu(I,j) * (u_inst(I,j,k) + dt_pred * &
                      (u_bc_accel(I,j,k) + CS%u_accel_bt(I,j,k)))
    enddo ; enddo
  enddo
  call cpu_clock_end(id_clock_mom_update)

  if (CS%debug) then
    call uvchksum("Predictor 1 [uv]", up, vp, G%HI, haloshift=0, symmetric=sym, unscale=US%L_T_to_m_s)
    call hchksum(h, "Predictor 1 h", G%HI, haloshift=1, unscale=GV%H_to_MKS)
    call uvchksum("Predictor 1 [uv]h", uh, vh, G%HI,haloshift=2, &
                  symmetric=sym, unscale=GV%H_to_MKS*US%L_to_m**2*US%s_to_T)
!   call MOM_state_chksum("Predictor 1", up, vp, h, uh, vh, G, GV, US, haloshift=1)
    call MOM_accel_chksum("Predictor accel", CS%CAu_pred, CS%CAv_pred, CS%PFu, CS%PFv, &
             CS%diffu, CS%diffv, G, GV, US, CS%pbce, CS%u_accel_bt, CS%v_accel_bt, symmetric=sym)
    call MOM_state_chksum("Predictor 1 init", u_inst, v_inst, h, uh, vh, G, GV, US, haloshift=1, &
                          symmetric=sym)
    if (debug_redundant) then
      call check_redundant("Predictor 1 up", up, vp, G, unscale=US%L_T_to_m_s)
      call check_redundant("Predictor 1 uh", uh, vh, G, unscale=GV%H_to_MKS*US%L_to_m**2*US%s_to_T)
    endif
  endif

! up <- up + dt_pred d/dz visc d/dz up
! u_av  <- u_av  + dt_pred d/dz visc d/dz u_av
  call cpu_clock_begin(id_clock_vertvisc)
  if (CS%debug) then
    call uvchksum("0 before vertvisc: [uv]p", up, vp, G%HI,haloshift=0, symmetric=sym, unscale=US%L_T_to_m_s)
  endif

  if (CS%fpmix) then
    uold(:,:,:) = 0.0
    vold(:,:,:) = 0.0
    do k = 1, nz
      do j = js , je
        do I = Isq, Ieq
          uold(I,j,k)   = up(I,j,k)
        enddo
      enddo
      do J = Jsq, Jeq
        do i = is, ie
          vold(i,J,k)   = vp(i,J,k)
        enddo
      enddo
    enddo
  endif

  call thickness_to_dz(h, tv, dz, G, GV, US, halo_size=1)
  call vertvisc_coef(up, vp, h, dz, forces, visc, tv, dt_pred, G, GV, US, CS%vertvisc_CSp, &
                     CS%OBC, VarMix)

  if (CS%fpmix) then
    hbl(:,:) = 0.0
    if (ASSOCIATED(CS%KPP_CSp)) call KPP_get_BLD(CS%KPP_CSp, hbl, G, US, m_to_BLD_units=GV%m_to_H)
    if (ASSOCIATED(CS%energetic_PBL_CSp)) &
      call energetic_PBL_get_MLD(CS%energetic_PBL_CSp, hbl, G, US, m_to_MLD_units=GV%m_to_H)

    ! lFPpost must be false in the predictor step to avoid averaging into the diagnostics
    lFPpost = .false.
    call vertFPmix(up, vp, uold, vold, hbl, h, forces, dt_pred, lFPpost, CS%Cemp_NL,  &
                   G, GV, US, CS%vertvisc_CSp, CS%OBC, waves=waves)
    call vertvisc(up, vp, h, forces, visc, dt_pred, CS%OBC, CS%AD_pred, CS%CDp, G, &
                  GV, US, CS%vertvisc_CSp, CS%taux_bot, CS%tauy_bot, fpmix=CS%fpmix, waves=waves)
  else
    call vertvisc(up, vp, h, forces, visc, dt_pred, CS%OBC, CS%AD_pred, CS%CDp, G, &
                  GV, US, CS%vertvisc_CSp, CS%taux_bot, CS%tauy_bot, waves=waves)
  endif

  if (showCallTree) call callTree_wayPoint("done with vertvisc (step_MOM_dyn_split_RK2)")
  if (G%nonblocking_updates) then
    call cpu_clock_end(id_clock_vertvisc)
    call start_group_pass(CS%pass_uvp, G%Domain, clock=id_clock_pass)
    call cpu_clock_begin(id_clock_vertvisc)
  endif
  if (CS%visc_rem_dt_bug) then
    call vertvisc_remnant(visc, CS%visc_rem_u, CS%visc_rem_v, dt_pred, G, GV, US, CS%vertvisc_CSp)
  else
    call vertvisc_remnant(visc, CS%visc_rem_u, CS%visc_rem_v, dt, G, GV, US, CS%vertvisc_CSp)
  endif
  call cpu_clock_end(id_clock_vertvisc)

  call do_group_pass(CS%pass_visc_rem, G%Domain, clock=id_clock_pass)
  if (G%nonblocking_updates) then
    call complete_group_pass(CS%pass_uvp, G%Domain, clock=id_clock_pass)
  else
    call do_group_pass(CS%pass_uvp, G%Domain, clock=id_clock_pass)
  endif

  ! uh = u_av * h
  ! hp = h + dt * div . uh
  call cpu_clock_begin(id_clock_continuity)
  call continuity(up, vp, h, hp, uh, vh, dt, G, GV, US, CS%continuity_CSp, CS%OBC, pbv, &
                  uhbt=CS%uhbt, vhbt=CS%vhbt, visc_rem_u=CS%visc_rem_u, visc_rem_v=CS%visc_rem_v, &
                  u_cor=u_av, v_cor=v_av, BT_cont=CS%BT_cont)
  call cpu_clock_end(id_clock_continuity)
  if (showCallTree) call callTree_wayPoint("done with continuity (step_MOM_dyn_split_RK2)")

  call do_group_pass(CS%pass_hp_uv, G%Domain, clock=id_clock_pass)

  if (associated(CS%OBC)) then

    if (CS%debug) &
      call uvchksum("Pre OBC avg [uv]", u_av, v_av, G%HI, haloshift=1, symmetric=sym, unscale=US%L_T_to_m_s)

    call radiation_open_bdry_conds(CS%OBC, u_av, u_old_rad_OBC, v_av, v_old_rad_OBC, G, GV, US, dt_pred)

    if (CS%debug) &
      call uvchksum("Post OBC avg [uv]", u_av, v_av, G%HI, haloshift=1, symmetric=sym, unscale=US%L_T_to_m_s)

    ! These should be done with a pass that excludes uh & vh.
!   call do_group_pass(CS%pass_hp_uv, G%Domain, clock=id_clock_pass)
  endif

  if (G%nonblocking_updates) then
    call start_group_pass(CS%pass_av_uvh, G%Domain, clock=id_clock_pass)
  endif

  ! h_av = (h + hp)/2
  !$OMP parallel do default(shared)
  do k=1,nz ; do j=js-2,je+2 ; do i=is-2,ie+2
    h_av(i,j,k) = 0.5*(h(i,j,k) + hp(i,j,k))
  enddo ; enddo ; enddo

  ! The correction phase of the time step starts here.
  call enable_averages(dt, Time_local, CS%diag)

  ! Calculate a revised estimate of the free-surface height correction to be
  ! used in the next call to btstep.  This call is at this point so that
  ! hp can be changed if CS%begw /= 0.
  ! eta_cor = ...                 (hidden inside CS%barotropic_CSp)
  call cpu_clock_begin(id_clock_btcalc)
  call bt_mass_source(hp, eta_pred, .false., G, GV, CS%barotropic_CSp)
  call cpu_clock_end(id_clock_btcalc)

  if (CS%begw /= 0.0) then
    ! hp <- (1-begw)*h_in + begw*hp
    ! Back up hp to the value it would have had after a time-step of
    ! begw*dt.  hp is not used again until recalculated by continuity.
    !$OMP parallel do default(shared)
    do k=1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
      hp(i,j,k) = (1.0-CS%begw)*h(i,j,k) + CS%begw*hp(i,j,k)
    enddo ; enddo ; enddo

    ! PFu = d/dx M(hp,T,S)
    ! pbce = dM/deta
    call cpu_clock_begin(id_clock_pres)
    call PressureForce(hp, tv, CS%PFu, CS%PFv, G, GV, US, CS%PressureForce_CSp, &
                       CS%ALE_CSp, CS%ADp, p_surf, CS%pbce, CS%eta_PF)
    ! Stokes shear force contribution to pressure gradient
    Use_Stokes_PGF = present(Waves)
    if (Use_Stokes_PGF) then
      Use_Stokes_PGF = associated(Waves)
      if (Use_Stokes_PGF) Use_Stokes_PGF = Waves%Stokes_PGF
      if (Use_Stokes_PGF) then
        call thickness_to_dz(h, tv, dz, G, GV, US, halo_size=1)
        call Stokes_PGF(G, GV, US, dz, u_inst, v_inst, CS%PFu_Stokes, CS%PFv_Stokes, Waves)
        if (.not.Waves%Passive_Stokes_PGF) then
          do k=1,nz
            do j=js,je ; do I=Isq,Ieq
              CS%PFu(I,j,k) = CS%PFu(I,j,k) + CS%PFu_Stokes(I,j,k)
            enddo ; enddo
          enddo
          do k=1,nz
            do J=Jsq,Jeq ; do i=is,ie
              CS%PFv(i,J,k) = CS%PFv(i,J,k) + CS%PFv_Stokes(i,J,k)
            enddo ; enddo
          enddo
        endif
      endif
    endif
    call cpu_clock_end(id_clock_pres)
    if (showCallTree) call callTree_wayPoint("done with PressureForce[hp=(1-b).h+b.h] (step_MOM_dyn_split_RK2)")
  endif

  if (G%nonblocking_updates) &
    call complete_group_pass(CS%pass_av_uvh, G%Domain, clock=id_clock_pass)

  if (BT_cont_BT_thick) then
    call btcalc(h, G, GV, CS%barotropic_CSp, CS%BT_cont%h_u, CS%BT_cont%h_v, &
                OBC=CS%OBC)
    if (showCallTree) call callTree_wayPoint("done with btcalc[BT_cont_BT_thick] (step_MOM_dyn_split_RK2)")
  endif

  if (CS%debug) then
    call MOM_state_chksum("Predictor ", up, vp, hp, uh, vh, G, GV, US, symmetric=sym)
    call uvchksum("Predictor avg [uv]", u_av, v_av, G%HI, haloshift=1, symmetric=sym, unscale=US%L_T_to_m_s)
    call hchksum(h_av, "Predictor avg h", G%HI, haloshift=2, unscale=GV%H_to_MKS)
  ! call MOM_state_chksum("Predictor avg ", u_av, v_av, h_av, uh, vh, G, GV, US)
    if (debug_redundant) then
      call check_redundant("Predictor up ", up, vp, G, unscale=US%L_T_to_m_s)
      call check_redundant("Predictor uh ", uh, vh, G, unscale=GV%H_to_MKS*US%L_to_m**2*US%s_to_T)
    endif
  endif

! diffu = horizontal viscosity terms (u_av)
  call cpu_clock_begin(id_clock_horvisc)
  call horizontal_viscosity(u_av, v_av, h_av, uh, vh, CS%diffu, CS%diffv, &
                            MEKE, Varmix, G, GV, US, CS%hor_visc, tv, dt, &
                            OBC=CS%OBC, BT=CS%barotropic_CSp, TD=thickness_diffuse_CSp, &
                            ADp=CS%ADp, hu_cont=CS%BT_cont%h_u, hv_cont=CS%BT_cont%h_v, STOCH=STOCH)
  call cpu_clock_end(id_clock_horvisc)
  if (showCallTree) call callTree_wayPoint("done with horizontal_viscosity (step_MOM_dyn_split_RK2)")

! CAu = -(f+zeta_av)/h_av vh + d/dx KE_av
  call cpu_clock_begin(id_clock_Cor)
  call CorAdCalc(u_av, v_av, h_av, uh, vh, CS%CAu, CS%CAv, CS%OBC, CS%ADp, &
                 G, GV, US, CS%CoriolisAdv, pbv, Waves=Waves)
  call cpu_clock_end(id_clock_Cor)
  if (showCallTree) call callTree_wayPoint("done with CorAdCalc (step_MOM_dyn_split_RK2)")

! Calculate the momentum forcing terms for the barotropic equations.

! u_bc_accel = CAu + PFu + diffu(u[n-1])
  call cpu_clock_begin(id_clock_btforce)
  !$OMP parallel do default(shared)
  do k=1,nz
    do j=js,je ; do I=Isq,Ieq
      u_bc_accel(I,j,k) = (CS%Cau(I,j,k) + CS%PFu(I,j,k)) + CS%diffu(I,j,k)
    enddo ; enddo
    do J=Jsq,Jeq ; do i=is,ie
      v_bc_accel(i,J,k) = (CS%Cav(i,J,k) + CS%PFv(i,J,k)) + CS%diffv(i,J,k)
    enddo ; enddo
  enddo
  if (associated(CS%OBC)) then
    call open_boundary_zero_normal_flow(CS%OBC, G, GV, u_bc_accel, v_bc_accel)
  endif
  call cpu_clock_end(id_clock_btforce)

  if (CS%debug) then
    call MOM_accel_chksum("corr pre-btstep accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv, &
                          CS%diffu, CS%diffv, G, GV, US, CS%pbce, u_bc_accel, v_bc_accel, &
                          symmetric=sym)
    if (debug_redundant) then
      call check_redundant("corr pre-btstep CS%CA ", CS%CAu, CS%CAv, G, unscale=US%L_T2_to_m_s2)
      call check_redundant("corr pre-btstep CS%PF ", CS%PFu, CS%PFv, G, unscale=US%L_T2_to_m_s2)
      call check_redundant("corr pre-btstep CS%diff ", CS%diffu, CS%diffv, G, unscale=US%L_T2_to_m_s2)
      call check_redundant("corr pre-btstep u_bc_accel ", u_bc_accel, v_bc_accel, G, unscale=US%L_T2_to_m_s2)
    endif
  endif

  ! u_accel_bt = layer accelerations due to barotropic solver
  ! pbce = dM/deta
  call cpu_clock_begin(id_clock_btstep)
  if (CS%BT_use_layer_fluxes) then
    uh_ptr => uh ; vh_ptr => vh ; u_ptr => u_av ; v_ptr => v_av
  endif

  if (showCallTree) call callTree_enter("btstep(), MOM_barotropic.F90")
  ! This is the corrector step call to btstep.
  call btstep(u_inst, v_inst, eta, dt, u_bc_accel, v_bc_accel, forces, CS%pbce, CS%eta_PF, u_av, v_av, &
              CS%u_accel_bt, CS%v_accel_bt, eta_pred, CS%uhbt, CS%vhbt, G, GV, US, &
              CS%barotropic_CSp, CS%visc_rem_u, CS%visc_rem_v, SpV_avg, CS%ADp, CS%OBC, CS%BT_cont, &
              eta_PF_start, taux_bot, tauy_bot, uh_ptr, vh_ptr, u_ptr, v_ptr, etaav=eta_av)
  if (CS%id_deta_dt>0) then
    do j=js,je ; do i=is,ie ; deta_dt(i,j) = (eta_pred(i,j) - eta(i,j))*Idt_bc ; enddo ; enddo
  endif
  do j=js,je ; do i=is,ie ; eta(i,j) = eta_pred(i,j) ; enddo ; enddo

  call cpu_clock_end(id_clock_btstep)
  if (showCallTree) call callTree_leave("btstep()")

  if (CS%debug .and. debug_redundant) then
    call check_redundant("u_accel_bt ", CS%u_accel_bt, CS%v_accel_bt, G, unscale=US%L_T2_to_m_s2)
  endif

  ! u = u + dt*( u_bc_accel + u_accel_bt )
  call cpu_clock_begin(id_clock_mom_update)
  !$OMP parallel do default(shared)
  do k=1,nz
    do j=js,je ; do I=Isq,Ieq
      u_inst(I,j,k) = G%mask2dCu(I,j) * (u_inst(I,j,k) + dt * &
                      (u_bc_accel(I,j,k) + CS%u_accel_bt(I,j,k)))
    enddo ; enddo
    do J=Jsq,Jeq ; do i=is,ie
      v_inst(i,J,k) = G%mask2dCv(i,J) * (v_inst(i,J,k) + dt * &
                      (v_bc_accel(i,J,k) + CS%v_accel_bt(i,J,k)))
    enddo ; enddo
  enddo
  call cpu_clock_end(id_clock_mom_update)

  if (CS%debug) then
    call uvchksum("Corrector 1 [uv]", u_inst, v_inst, G%HI, haloshift=0, symmetric=sym, unscale=US%L_T_to_m_s)
    call hchksum(h, "Corrector 1 h", G%HI, haloshift=1, unscale=GV%H_to_MKS)
    call uvchksum("Corrector 1 [uv]h", uh, vh, G%HI, haloshift=2, &
                  symmetric=sym, unscale=GV%H_to_MKS*US%L_to_m**2*US%s_to_T)
  ! call MOM_state_chksum("Corrector 1", u_inst, v_inst, h, uh, vh, G, GV, US, haloshift=1)
    call MOM_accel_chksum("Corrector accel", CS%CAu, CS%CAv, CS%PFu, CS%PFv, &
                          CS%diffu, CS%diffv, G, GV, US, CS%pbce, CS%u_accel_bt, CS%v_accel_bt, &
                          symmetric=sym)
  endif

  ! u <- u + dt d/dz visc d/dz u
  ! u_av <- u_av + dt d/dz visc d/dz u_av
  call cpu_clock_begin(id_clock_vertvisc)

  if (CS%fpmix) then
    uold(:,:,:) = 0.0
    vold(:,:,:) = 0.0
    do k = 1, nz
      do j = js , je
        do I = Isq, Ieq
          uold(I,j,k)   = u_inst(I,j,k)
        enddo
      enddo
      do J = Jsq, Jeq
        do i = is, ie
          vold(i,J,k)   = v_inst(i,J,k)
        enddo
      enddo
    enddo
  endif

  call thickness_to_dz(h, tv, dz, G, GV, US, halo_size=1)
  call vertvisc_coef(u_inst, v_inst, h, dz, forces, visc, tv, dt, G, GV, US, CS%vertvisc_CSp, CS%OBC, VarMix)

  if (CS%fpmix) then
    lFPpost = .true.
    call vertFPmix(u_inst, v_inst, uold, vold, hbl, h, forces, dt, lFPpost, CS%Cemp_NL, &
                   G, GV, US, CS%vertvisc_CSp, CS%OBC, Waves=Waves)
    call vertvisc(u_inst, v_inst, h, forces, visc, dt, CS%OBC, CS%ADp, CS%CDp, G, GV, US, &
         CS%vertvisc_CSp, CS%taux_bot, CS%tauy_bot, fpmix=CS%fpmix, waves=waves)

  else
    call vertvisc(u_inst, v_inst, h, forces, visc, dt, CS%OBC, CS%ADp, CS%CDp, G, GV, US, &
                  CS%vertvisc_CSp, CS%taux_bot, CS%tauy_bot, waves=waves)
  endif

  if (G%nonblocking_updates) then
    call cpu_clock_end(id_clock_vertvisc)
    call start_group_pass(CS%pass_uv, G%Domain, clock=id_clock_pass)
    call cpu_clock_begin(id_clock_vertvisc)
  endif
  call vertvisc_remnant(visc, CS%visc_rem_u, CS%visc_rem_v, dt, G, GV, US, CS%vertvisc_CSp)
  call cpu_clock_end(id_clock_vertvisc)
  if (showCallTree) call callTree_wayPoint("done with vertvisc (step_MOM_dyn_split_RK2)")

! Later, h_av = (h_in + h_out)/2, but for now use h_av to store h_in.
  !$OMP parallel do default(shared)
  do k=1,nz ; do j=js-2,je+2 ; do i=is-2,ie+2
    h_av(i,j,k) = h(i,j,k)
  enddo ; enddo ; enddo

  call do_group_pass(CS%pass_visc_rem, G%Domain, clock=id_clock_pass)
  if (G%nonblocking_updates) then
    call complete_group_pass(CS%pass_uv, G%Domain, clock=id_clock_pass)
  else
    call do_group_pass(CS%pass_uv, G%Domain, clock=id_clock_pass)
  endif

  ! uh = u_av * h
  ! h  = h + dt * div . uh
  ! u_av and v_av adjusted so their mass transports match uhbt and vhbt.
  call cpu_clock_begin(id_clock_continuity)
  call continuity(u_inst, v_inst, h, h, uh, vh, dt, G, GV, US, CS%continuity_CSp, CS%OBC, pbv, &
                  uhbt=CS%uhbt, vhbt=CS%vhbt, visc_rem_u=CS%visc_rem_u, visc_rem_v=CS%visc_rem_v, &
                  u_cor=u_av, v_cor=v_av)
  call cpu_clock_end(id_clock_continuity)
  call do_group_pass(CS%pass_h, G%Domain, clock=id_clock_pass)
  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  call diag_update_remap_grids(CS%diag)
  if (showCallTree) call callTree_wayPoint("done with continuity (step_MOM_dyn_split_RK2)")

  if (G%nonblocking_updates) then
    call start_group_pass(CS%pass_av_uvh, G%Domain, clock=id_clock_pass)
  else
    call do_group_pass(CS%pass_av_uvh, G%domain, clock=id_clock_pass)
  endif

  if (associated(CS%OBC)) then
    !### I suspect that there is a bug here when u_inst is compared with a previous value of u_av
    ! to estimate the dominant outward group velocity, but a fix is not available yet.
    call radiation_open_bdry_conds(CS%OBC, u_inst, u_old_rad_OBC, v_inst, v_old_rad_OBC, G, GV, US, dt)
  endif

! h_av = (h_in + h_out)/2 . Going in to this line, h_av = h_in.
  !$OMP parallel do default(shared)
  do k=1,nz ; do j=js-2,je+2 ; do i=is-2,ie+2
    h_av(i,j,k) = 0.5*(h_av(i,j,k) + h(i,j,k))
  enddo ; enddo ; enddo

  if (G%nonblocking_updates) &
    call complete_group_pass(CS%pass_av_uvh, G%Domain, clock=id_clock_pass)

  !$OMP parallel do default(shared)
  do k=1,nz
    do j=js-2,je+2 ; do I=Isq-2,Ieq+2
      uhtr(I,j,k) = uhtr(I,j,k) + uh(I,j,k)*dt
    enddo ; enddo
    do J=Jsq-2,Jeq+2 ; do i=is-2,ie+2
      vhtr(i,J,k) = vhtr(i,J,k) + vh(i,J,k)*dt
    enddo ; enddo
  enddo

  if (CS%store_CAu) then
    ! Calculate a predictor-step estimate of the Coriolis and momentum advection terms
    ! for use in the next time step, possibly after it has been vertically remapped.
    call cpu_clock_begin(id_clock_Cor)
    call disable_averaging(CS%diag)  ! These calculations should not be used for diagnostics.
    ! CAu = -(f+zeta_av)/h_av vh + d/dx KE_av
    call CorAdCalc(u_av, v_av, h_av, uh, vh, CS%CAu_pred, CS%CAv_pred, CS%OBC, CS%AD_pred, &
                   G, GV, US, CS%CoriolisAdv, pbv, Waves=Waves)
    CS%CAu_pred_stored = .true.
    call enable_averages(dt, Time_local, CS%diag) ! Reenable the averaging
    call cpu_clock_end(id_clock_Cor)
    if (showCallTree) call callTree_wayPoint("done with CorAdCalc (step_MOM_dyn_split_RK2)")
  else
    CS%CAu_pred_stored = .false.
  endif

  ! The time-averaged free surface height has already been set by the last call to btstep.

  ! Deallocate this memory to avoid a memory leak. ### We should revisit how this array is declared. -RWH
  if (dyn_p_surf .and. associated(eta_PF_start)) deallocate(eta_PF_start)

  !  Here various terms used in to update the momentum equations are
  !  offered for time averaging.
  if (CS%id_PFu > 0) call post_data(CS%id_PFu, CS%PFu, CS%diag)
  if (CS%id_PFv > 0) call post_data(CS%id_PFv, CS%PFv, CS%diag)
  if (CS%id_CAu > 0) call post_data(CS%id_CAu, CS%CAu, CS%diag)
  if (CS%id_CAv > 0) call post_data(CS%id_CAv, CS%CAv, CS%diag)

  ! Here the thickness fluxes are offered for time averaging.
  if (CS%id_uh         > 0) call post_data(CS%id_uh,  uh,                   CS%diag)
  if (CS%id_vh         > 0) call post_data(CS%id_vh,  vh,                   CS%diag)
  if (CS%id_uav        > 0) call post_data(CS%id_uav, u_av,                 CS%diag)
  if (CS%id_vav        > 0) call post_data(CS%id_vav, v_av,                 CS%diag)
  if (CS%id_u_BT_accel > 0) call post_data(CS%id_u_BT_accel, CS%u_accel_bt, CS%diag)
  if (CS%id_v_BT_accel > 0) call post_data(CS%id_v_BT_accel, CS%v_accel_bt, CS%diag)

  ! Calculate effective areas and post data
  if (CS%id_ueffA > 0) then
    ueffA(:,:,:) = 0
    do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      if (abs(up(I,j,k)) > 0.) ueffA(I,j,k) = uh(I,j,k) / up(I,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_ueffA, ueffA, CS%diag)
  endif

  if (CS%id_veffA > 0) then
    veffA(:,:,:) = 0
    do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      if (abs(vp(i,J,k)) > 0.) veffA(i,J,k) = vh(i,J,k) / vp(i,J,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_veffA, veffA, CS%diag)
  endif

  ! Diagnostics of the fractional thicknesses times momentum budget terms
  ! 3D diagnostics hf_PFu etc. are commented because there is no clarity on proper remapping grid option.
  ! The code is retained for debugging purposes in the future.
  !if (CS%id_hf_PFu > 0) call post_product_u(CS%id_hf_PFu, CS%PFu, CS%ADp%diag_hfrac_u, G, nz, CS%diag)
  !if (CS%id_hf_PFv > 0) call post_product_v(CS%id_hf_PFv, CS%PFv, CS%ADp%diag_hfrac_v, G, nz, CS%diag)
  !if (CS%id_hf_CAu > 0) call post_product_u(CS%id_hf_CAu, CS%CAu, CS%ADp%diag_hfrac_u, G, nz, CS%diag)
  !if (CS%id_hf_CAv > 0) call post_product_v(CS%id_hf_CAv, CS%CAv, CS%ADp%diag_hfrac_v, G, nz, CS%diag)
  !if (CS%id_hf_u_BT_accel > 0) &
  !  call post_product_u(CS%id_hf_u_BT_accel, CS%u_accel_bt, CS%ADp%diag_hfrac_u, G, nz, CS%diag)
  !if (CS%id_hf_v_BT_accel > 0) &
  !  call post_product_v(CS%id_hf_v_BT_accel, CS%v_accel_bt, CS%ADp%diag_hfrac_v, G, nz, CS%diag)

  ! Diagnostics for the vertical sum of layer thickness x prssure force accelerations
  if (CS%id_intz_PFu_2d > 0) call post_product_sum_u(CS%id_intz_PFu_2d, CS%PFu, CS%ADp%diag_hu, G, nz, CS%diag)
  if (CS%id_intz_PFv_2d > 0) call post_product_sum_v(CS%id_intz_PFv_2d, CS%PFv, CS%ADp%diag_hv, G, nz, CS%diag)

  ! Diagnostics for thickness-weighted vertically averaged prssure force accelerations
  if (CS%id_hf_PFu_2d > 0) call post_product_sum_u(CS%id_hf_PFu_2d, CS%PFu, CS%ADp%diag_hfrac_u, G, nz, CS%diag)
  if (CS%id_hf_PFv_2d > 0) call post_product_sum_v(CS%id_hf_PFv_2d, CS%PFv, CS%ADp%diag_hfrac_v, G, nz, CS%diag)

  ! Diagnostics for thickness x prssure force accelerations
  if (CS%id_h_PFu > 0) call post_product_u(CS%id_h_PFu, CS%PFu, CS%ADp%diag_hu, G, nz, CS%diag)
  if (CS%id_h_PFv > 0) call post_product_v(CS%id_h_PFv, CS%PFv, CS%ADp%diag_hv, G, nz, CS%diag)

  ! Diagnostics of Coriolis acceleratations
  if (CS%id_intz_CAu_2d > 0) call post_product_sum_u(CS%id_intz_CAu_2d, CS%CAu, CS%ADp%diag_hu, G, nz, CS%diag)
  if (CS%id_intz_CAv_2d > 0) call post_product_sum_v(CS%id_intz_CAv_2d, CS%CAv, CS%ADp%diag_hv, G, nz, CS%diag)
  if (CS%id_hf_CAu_2d > 0) call post_product_sum_u(CS%id_hf_CAu_2d, CS%CAu, CS%ADp%diag_hfrac_u, G, nz, CS%diag)
  if (CS%id_hf_CAv_2d > 0) call post_product_sum_v(CS%id_hf_CAv_2d, CS%CAv, CS%ADp%diag_hfrac_v, G, nz, CS%diag)
  if (CS%id_h_CAu > 0) call post_product_u(CS%id_h_CAu, CS%CAu, CS%ADp%diag_hu, G, nz, CS%diag)
  if (CS%id_h_CAv > 0) call post_product_v(CS%id_h_CAv, CS%CAv, CS%ADp%diag_hv, G, nz, CS%diag)

  ! Diagnostics of barotropic solver acceleratations
  if (CS%id_intz_u_BT_accel_2d > 0) &
    call post_product_sum_u(CS%id_intz_u_BT_accel_2d, CS%u_accel_bt, CS%ADp%diag_hu, G, nz, CS%diag)
  if (CS%id_intz_v_BT_accel_2d > 0) &
    call post_product_sum_v(CS%id_intz_v_BT_accel_2d, CS%v_accel_bt, CS%ADp%diag_hv, G, nz, CS%diag)
  if (CS%id_hf_u_BT_accel_2d > 0) &
    call post_product_sum_u(CS%id_hf_u_BT_accel_2d, CS%u_accel_bt, CS%ADp%diag_hfrac_u, G, nz, CS%diag)
  if (CS%id_hf_v_BT_accel_2d > 0) &
    call post_product_sum_v(CS%id_hf_v_BT_accel_2d, CS%v_accel_bt, CS%ADp%diag_hfrac_v, G, nz, CS%diag)
  if (CS%id_h_u_BT_accel > 0) &
    call post_product_u(CS%id_h_u_BT_accel, CS%u_accel_bt, CS%ADp%diag_hu, G, nz, CS%diag)
  if (CS%id_h_v_BT_accel > 0) &
    call post_product_v(CS%id_h_v_BT_accel, CS%v_accel_bt, CS%ADp%diag_hv, G, nz, CS%diag)

  ! Diagnostics for momentum budget terms multiplied by visc_rem_[uv],
  if (CS%id_PFu_visc_rem > 0) call post_product_u(CS%id_PFu_visc_rem, CS%PFu, CS%ADp%visc_rem_u, G, nz, CS%diag)
  if (CS%id_PFv_visc_rem > 0) call post_product_v(CS%id_PFv_visc_rem, CS%PFv, CS%ADp%visc_rem_v, G, nz, CS%diag)
  if (CS%id_CAu_visc_rem > 0) call post_product_u(CS%id_CAu_visc_rem, CS%CAu, CS%ADp%visc_rem_u, G, nz, CS%diag)
  if (CS%id_CAv_visc_rem > 0) call post_product_v(CS%id_CAv_visc_rem, CS%CAv, CS%ADp%visc_rem_v, G, nz, CS%diag)
  if (CS%id_u_BT_accel_visc_rem > 0) &
    call post_product_u(CS%id_u_BT_accel_visc_rem, CS%u_accel_bt, CS%ADp%visc_rem_u, G, nz, CS%diag)
  if (CS%id_v_BT_accel_visc_rem > 0) &
    call post_product_v(CS%id_v_BT_accel_visc_rem, CS%v_accel_bt, CS%ADp%visc_rem_v, G, nz, CS%diag)

  ! Diagnostics related to changes in eta
  if (CS%id_deta_dt > 0) call post_data(CS%id_deta_dt, deta_dt, CS%diag)

  if (CS%debug) then
    call MOM_state_chksum("Corrector ", u_inst, v_inst, h, uh, vh, G, GV, US, symmetric=sym)
    call uvchksum("Corrector avg [uv]", u_av, v_av, G%HI, haloshift=1, symmetric=sym, unscale=US%L_T_to_m_s)
    call hchksum(h_av, "Corrector avg h", G%HI, haloshift=1, unscale=GV%H_to_MKS)
 !  call MOM_state_chksum("Corrector avg ", u_av, v_av, h_av, uh, vh, G, GV, US)
  endif

  if (showCallTree) call callTree_leave("step_MOM_dyn_split_RK2()")

end subroutine step_MOM_dyn_split_RK2

!> This subroutine sets up any auxiliary restart variables that are specific
!! to the split-explicit time stepping scheme.  All variables registered here should
!! have the ability to be recreated if they are not present in a restart file.
subroutine register_restarts_dyn_split_RK2(HI, GV, US, param_file, CS, restart_CS, uh, vh)
  type(hor_index_type),          intent(in)    :: HI         !< Horizontal index structure
  type(verticalGrid_type),       intent(in)    :: GV         !< ocean vertical grid structure
  type(unit_scale_type),         intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),         intent(in)    :: param_file !< parameter file
  type(MOM_dyn_split_RK2_CS),    pointer       :: CS         !< module control structure
  type(MOM_restart_CS),          intent(inout) :: restart_CS !< MOM restart control structure
  real, dimension(SZIB_(HI),SZJ_(HI),SZK_(GV)), &
                         target, intent(inout) :: uh !< zonal volume or mass transport [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(HI),SZJB_(HI),SZK_(GV)), &
                         target, intent(inout) :: vh !< merid volume or mass transport [H L2 T-1 ~> m3 s-1 or kg s-1]

  character(len=40) :: mdl = "MOM_dynamics_split_RK2" ! This module's name.
  type(vardesc)     :: vd(2)
  character(len=48) :: thickness_units, flux_units

  integer :: isd, ied, jsd, jed, nz, IsdB, IedB, JsdB, JedB

  isd  = HI%isd  ; ied  = HI%ied  ; jsd  = HI%jsd  ; jed  = HI%jed ; nz = GV%ke
  IsdB = HI%IsdB ; IedB = HI%IedB ; JsdB = HI%JsdB ; JedB = HI%JedB

  ! This is where a control structure specific to this module would be allocated.
  if (associated(CS)) then
    call MOM_error(WARNING, "register_restarts_dyn_split_RK2 called with an associated "// &
                             "control structure.")
    return
  endif
  allocate(CS)

  ALLOC_(CS%diffu(IsdB:IedB,jsd:jed,nz)) ; CS%diffu(:,:,:) = 0.0
  ALLOC_(CS%diffv(isd:ied,JsdB:JedB,nz)) ; CS%diffv(:,:,:) = 0.0
  ALLOC_(CS%CAu(IsdB:IedB,jsd:jed,nz))   ; CS%CAu(:,:,:)   = 0.0
  ALLOC_(CS%CAv(isd:ied,JsdB:JedB,nz))   ; CS%CAv(:,:,:)   = 0.0
  ALLOC_(CS%CAu_pred(IsdB:IedB,jsd:jed,nz)) ; CS%CAu_pred(:,:,:)   = 0.0
  ALLOC_(CS%CAv_pred(isd:ied,JsdB:JedB,nz)) ; CS%CAv_pred(:,:,:)   = 0.0
  ALLOC_(CS%PFu(IsdB:IedB,jsd:jed,nz))   ; CS%PFu(:,:,:)   = 0.0
  ALLOC_(CS%PFv(isd:ied,JsdB:JedB,nz))   ; CS%PFv(:,:,:)   = 0.0

  ALLOC_(CS%eta(isd:ied,jsd:jed))       ; CS%eta(:,:)    = 0.0
  ALLOC_(CS%u_av(IsdB:IedB,jsd:jed,nz)) ; CS%u_av(:,:,:) = 0.0
  ALLOC_(CS%v_av(isd:ied,JsdB:JedB,nz)) ; CS%v_av(:,:,:) = 0.0
  ALLOC_(CS%h_av(isd:ied,jsd:jed,nz))   ; CS%h_av(:,:,:) = GV%Angstrom_H

  thickness_units = get_thickness_units(GV)
  flux_units = get_flux_units(GV)

  call get_param(param_file, mdl, "STORE_CORIOLIS_ACCEL", CS%store_CAu, &
                 "If true, calculate the Coriolis accelerations at the end of each "//&
                 "timestep for use in the predictor step of the next split RK2 timestep.", &
                 default=.true., do_not_log=.true.)

  if (GV%Boussinesq) then
    call register_restart_field(CS%eta, "sfc", .false., restart_CS, &
             longname="Free surface Height", units=thickness_units, conversion=GV%H_to_mks)
  else
    call register_restart_field(CS%eta, "p_bot", .false., restart_CS, &
             longname="Bottom Pressure", units=thickness_units, conversion=GV%H_to_mks)
  endif

  ! These are needed, either to calculate CAu and CAv or to calculate the velocity anomalies in
  ! the barotropic solver's Coriolis terms.
  vd(1) = var_desc("u2", "m s-1", "Auxiliary Zonal velocity", 'u', 'L')
  vd(2) = var_desc("v2", "m s-1", "Auxiliary Meridional velocity", 'v', 'L')
  call register_restart_pair(CS%u_av, CS%v_av, vd(1), vd(2), .false., restart_CS, &
                             conversion=US%L_T_to_m_s)

  if (CS%store_CAu) then
    vd(1) = var_desc("CAu", "m s-2", "Zonal Coriolis and advactive acceleration", 'u', 'L')
    vd(2) = var_desc("CAv", "m s-2", "Meridional Coriolis and advactive  acceleration", 'v', 'L')
    call register_restart_pair(CS%CAu_pred, CS%CAv_pred, vd(1), vd(2), .false., restart_CS, &
                               conversion=US%L_T2_to_m_s2)
  else
    call register_restart_field(CS%h_av, "h2", .false., restart_CS, &
           longname="Auxiliary Layer Thickness", units=thickness_units, conversion=GV%H_to_mks)

    vd(1) = var_desc("uh", flux_units, "Zonal thickness flux", 'u', 'L')
    vd(2) = var_desc("vh", flux_units, "Meridional thickness flux", 'v', 'L')
    call register_restart_pair(uh, vh, vd(1), vd(2), .false., restart_CS, &
                               conversion=GV%H_to_MKS*US%L_to_m**2*US%s_to_T)
  endif

  vd(1) = var_desc("diffu", "m s-2", "Zonal horizontal viscous acceleration", 'u', 'L')
  vd(2) = var_desc("diffv", "m s-2", "Meridional horizontal viscous acceleration", 'v', 'L')
  call register_restart_pair(CS%diffu, CS%diffv, vd(1), vd(2), .false., restart_CS, &
                             conversion=US%L_T2_to_m_s2)

  call register_barotropic_restarts(HI, GV, US, param_file, CS%barotropic_CSp, restart_CS)

end subroutine register_restarts_dyn_split_RK2

!> This subroutine does remapping for the auxiliary restart variables that are used
!! with the split RK2 time stepping scheme.
subroutine remap_dyn_split_RK2_aux_vars(G, GV, CS, h_old_u, h_old_v, h_new_u, h_new_v, ALE_CSp)
  type(ocean_grid_type),            intent(inout) :: G        !< ocean grid structure
  type(verticalGrid_type),          intent(in)    :: GV       !< ocean vertical grid structure
  type(MOM_dyn_split_RK2_CS),       pointer       :: CS       !< module control structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: h_old_u  !< Source grid thickness at zonal
                                                              !! velocity points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                    intent(in)    :: h_old_v  !< Source grid thickness at meridional
                                                              !! velocity points [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: h_new_u  !< Destination grid thickness at zonal
                                                              !! velocity points [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                    intent(in)    :: h_new_v  !< Destination grid thickness at meridional
                                                              !! velocity points [H ~> m or kg m-2]
  type(ALE_CS),                     pointer       :: ALE_CSp  !< ALE control structure to use when remapping

  if (.not.CS%remap_aux) return

  if (CS%store_CAu) then
    call ALE_remap_velocities(ALE_CSp, G, GV, h_old_u, h_old_v, h_new_u, h_new_v, CS%u_av, CS%v_av)
    call pass_vector(CS%u_av, CS%v_av, G%Domain, complete=.false.)
    call ALE_remap_velocities(ALE_CSp, G, GV, h_old_u, h_old_v, h_new_u, h_new_v, CS%CAu_pred, CS%CAv_pred)
    call pass_vector(CS%CAu_pred, CS%CAv_pred, G%Domain, complete=.true.)
  endif

  call ALE_remap_velocities(ALE_CSp, G, GV, h_old_u, h_old_v, h_new_u, h_new_v, CS%diffu, CS%diffv)

end subroutine remap_dyn_split_RK2_aux_vars

!> Initializes aspects of the dyn_split_RK2 that depend on diabatic processes.
!! Needed when BLDs are used in the dynamics.
subroutine init_dyn_split_RK2_diabatic(diabatic_CSp, CS)
  type(diabatic_CS),                intent(in) :: diabatic_CSp !< diabatic structure
  type(MOM_dyn_split_RK2_CS),       pointer    :: CS           !< module control structure

  call extract_diabatic_member(diabatic_CSp, KPP_CSp=CS%KPP_CSp)
  call extract_diabatic_member(diabatic_CSp, energetic_PBL_CSp=CS%energetic_PBL_CSp)

end subroutine init_dyn_split_RK2_diabatic

!> This subroutine initializes all of the variables that are used by this
!! dynamic core, including diagnostics and the cpu clocks.
subroutine initialize_dyn_split_RK2(u, v, h, tv, uh, vh, eta, Time, G, GV, US, param_file, &
                      diag, CS, HA_CSp, restart_CS, dt, Accel_diag, Cont_diag, MIS, &
                      VarMix, MEKE, thickness_diffuse_CSp,                  &
                      OBC, update_OBC_CSp, ALE_CSp, set_visc, &
                      visc, dirs, ntrunc, pbv, calc_dtbt, cont_stencil)
  type(ocean_grid_type),            intent(inout) :: G          !< ocean grid structure
  type(verticalGrid_type),          intent(in)    :: GV         !< ocean vertical grid structure
  type(unit_scale_type),            intent(in)    :: US         !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                    intent(inout) :: u          !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                    intent(inout) :: v          !< merid velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                                    intent(inout) :: h          !< layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),            intent(in)    :: tv         !< Thermodynamic type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            target, intent(inout) :: uh    !< zonal volume/mass transport [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                            target, intent(inout) :: vh    !< merid volume/mass transport [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: eta        !< free surface height or column mass [H ~> m or kg m-2]
  type(time_type),          target, intent(in)    :: Time       !< current model time
  type(param_file_type),            intent(in)    :: param_file !< parameter file for parsing
  type(diag_ctrl),          target, intent(inout) :: diag       !< to control diagnostics
  type(MOM_dyn_split_RK2_CS),       pointer       :: CS         !< module control structure
  type(harmonic_analysis_CS),       pointer       :: HA_CSp     !< A pointer to the control structure of the
                                                                !! harmonic analysis module
  type(MOM_restart_CS),             intent(inout) :: restart_CS !< MOM restart control structure
  real,                             intent(in)    :: dt         !< time step [T ~> s]
  type(accel_diag_ptrs),    target, intent(inout) :: Accel_diag !< points to momentum equation terms for
                                                                !! budget analysis
  type(cont_diag_ptrs),     target, intent(inout) :: Cont_diag  !< points to terms in continuity equation
  type(ocean_internal_state),       intent(inout) :: MIS        !< "MOM6 internal state" used to pass
                                                                !! diagnostic pointers
  type(VarMix_CS),                  intent(inout) :: VarMix     !< points to spatially variable viscosities
  type(MEKE_type),                  intent(inout) :: MEKE       !< MEKE fields
  type(thickness_diffuse_CS),       intent(inout) :: thickness_diffuse_CSp !< Pointer to the control structure
                                                                !! used for the isopycnal height diffusive transport.
  type(ocean_OBC_type),             pointer       :: OBC        !< points to OBC related fields
  type(update_OBC_CS),              pointer       :: update_OBC_CSp !< points to OBC update related fields
  type(ALE_CS),                     pointer       :: ALE_CSp    !< points to ALE control structure
  type(set_visc_CS),        target, intent(in)    :: set_visc   !< set_visc control structure
  type(vertvisc_type),              intent(inout) :: visc       !< vertical viscosities, bottom drag, and related
  type(directories),                intent(in)    :: dirs       !< contains directory paths
  integer, target,                  intent(inout) :: ntrunc     !< A target for the variable that records
                                                                !! the number of times the velocity is
                                                                !! truncated (this should be 0).
  logical,                          intent(out)   :: calc_dtbt  !< If true, recalculate the barotropic time step
  type(porous_barrier_type),        intent(in)    :: pbv        !< porous barrier fractional cell metrics
  integer,                          intent(out)   :: cont_stencil !< The stencil for thickness
                                                                !! from the continuity solver.

  ! local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_tmp ! A temporary copy of the layer thicknesses [H ~> m or kg m-2]
  character(len=40) :: mdl = "MOM_dynamics_split_RK2" ! This module's name.
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=48) :: thickness_units, flux_units, eta_rest_name
  type(group_pass_type) :: pass_av_h_uvh
  logical :: debug_truncations
  logical :: read_uv, read_h2
  logical :: enable_bugs  ! If true, the defaults for recently added bug-fix flags are set to
                          ! recreate the bugs, or if false bugs are only used if actively selected.
  logical :: visc_rem_bug ! Stores the value of runtime paramter VISC_REM_BUG.

  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec ; nz = GV%ke
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(CS)) call MOM_error(FATAL, &
      "initialize_dyn_split_RK2 called with an unassociated control structure.")
  if (CS%module_is_initialized) then
    call MOM_error(WARNING, "initialize_dyn_split_RK2 called with a control "// &
                            "structure that has already been initialized.")
    return
  endif
  CS%module_is_initialized = .true.

  CS%diag => diag

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "TIDES", CS%use_tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)
  call get_param(param_file, mdl, "CALCULATE_SAL", CS%calculate_SAL, &
                 "If true, calculate self-attraction and loading.", default=CS%use_tides)
  call get_param(param_file, mdl, "BE", CS%be, &
                 "If SPLIT is true, BE determines the relative weighting "//&
                 "of a  2nd-order Runga-Kutta baroclinic time stepping "//&
                 "scheme (0.5) and a backward Euler scheme (1) that is "//&
                 "used for the Coriolis and inertial terms.  BE may be "//&
                 "from 0.5 to 1, but instability may occur near 0.5. "//&
                 "BE is also applicable if SPLIT is false and USE_RK2 "//&
                 "is true.", units="nondim", default=0.6)
  call get_param(param_file, mdl, "BEGW", CS%begw, &
                 "If SPLIT is true, BEGW is a number from 0 to 1 that "//&
                 "controls the extent to which the treatment of gravity "//&
                 "waves is forward-backward (0) or simulated backward "//&
                 "Euler (1).  0 is almost always used. "//&
                 "If SPLIT is false and USE_RK2 is true, BEGW can be "//&
                 "between 0 and 0.5 to damp gravity waves.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "SET_DTBT_USE_BT_CONT", CS%dtbt_use_bt_cont, &
                 "If true, use BT_CONT to calculate DTBT if possible.", default=.false.)
  call get_param(param_file, mdl, "SPLIT_BOTTOM_STRESS", CS%split_bottom_stress, &
                 "If true, provide the bottom stress calculated by the "//&
                 "vertical viscosity to the barotropic solver.", default=.false.)
  call get_param(param_file, mdl, "BT_USE_LAYER_FLUXES", CS%BT_use_layer_fluxes, &
                 "If true, use the summed layered fluxes plus an "//&
                 "adjustment due to the change in the barotropic velocity "//&
                 "in the barotropic continuity equation.", default=.true.)
  call get_param(param_file, mdl, "STORE_CORIOLIS_ACCEL", CS%store_CAu, &
                 "If true, calculate the Coriolis accelerations at the end of each "//&
                 "timestep for use in the predictor step of the next split RK2 timestep.", &
                 default=.true.)
  call get_param(param_file, mdl, "FPMIX", CS%fpmix, &
                 "If true, add non-local momentum flux increments and diffuse down the Eulerian gradient.", &
                 default=.false.)
  if (CS%fpmix) then
    call get_param(param_file, "MOM", "CEMP_NL", CS%Cemp_NL, &
                 "Empirical coefficient of non-local momentum mixing.", &
                 units="nondim", default=3.6)
  endif
  call get_param(param_file, mdl, "REMAP_AUXILIARY_VARS", CS%remap_aux, &
                 "If true, apply ALE remapping to all of the auxiliary 3-dimensional "//&
                 "variables that are needed to reproduce across restarts, similarly to "//&
                 "what is already being done with the primary state variables.  "//&
                 "The default should be changed to true.", default=.false., do_not_log=.true.)
  if (CS%remap_aux .and. .not.CS%store_CAu) call MOM_error(FATAL, &
      "REMAP_AUXILIARY_VARS requires that STORE_CORIOLIS_ACCEL = True.")
  call get_param(param_file, mdl, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "OBC_DEBUGGING_TESTS", CS%debug_OBC, &
                 "If true, do additional calls resetting certain values to help verify the "//&
                 "correctness of the open boundary condition code.", &
                 default=.false., old_name="DEBUG_OBC", debuggingParam=.true., do_not_log=.true.)
  call get_param(param_file, mdl, "DEBUG_TRUNCATIONS", debug_truncations, &
                 default=.false.)
  call get_param(param_file, mdl, "ENABLE_BUGS_BY_DEFAULT", enable_bugs, &
                 default=.true., do_not_log=.true.)  ! This is logged from MOM.F90.
  call get_param(param_file, mdl, "VISC_REM_BUG", visc_rem_bug, &
                 "If true, visc_rem_[uv] in split mode is incorrectly calculated or accounted "//&
                 "for in two places. This parameter controls the defaults of two individual "//&
                 "flags, VISC_REM_TIMESTEP_BUG in MOM_dynamics_split_RK2(b) and "//&
                 "VISC_REM_BT_WEIGHT_BUG in MOM_barotropic.", default=.false.)
  call get_param(param_file, mdl, "VISC_REM_TIMESTEP_BUG", CS%visc_rem_dt_bug, &
                 "If true, recover a bug that uses dt_pred rather than dt in "//&
                 "vertvisc_remnant() at the end of predictor stage for the following "//&
                 "continuity() and btstep() calls in the corrector step. Default of this flag "//&
                 "is set by VISC_REM_BUG", default=visc_rem_bug)

  allocate(CS%taux_bot(IsdB:IedB,jsd:jed), source=0.0)
  allocate(CS%tauy_bot(isd:ied,JsdB:JedB), source=0.0)

  ALLOC_(CS%uhbt(IsdB:IedB,jsd:jed))          ; CS%uhbt(:,:)         = 0.0
  ALLOC_(CS%vhbt(isd:ied,JsdB:JedB))          ; CS%vhbt(:,:)         = 0.0
  ALLOC_(CS%visc_rem_u(IsdB:IedB,jsd:jed,nz)) ; CS%visc_rem_u(:,:,:) = 0.0
  ALLOC_(CS%visc_rem_v(isd:ied,JsdB:JedB,nz)) ; CS%visc_rem_v(:,:,:) = 0.0
  ALLOC_(CS%eta_PF(isd:ied,jsd:jed))          ; CS%eta_PF(:,:)       = 0.0
  ALLOC_(CS%pbce(isd:ied,jsd:jed,nz))         ; CS%pbce(:,:,:)       = 0.0

  ALLOC_(CS%u_accel_bt(IsdB:IedB,jsd:jed,nz)) ; CS%u_accel_bt(:,:,:) = 0.0
  ALLOC_(CS%v_accel_bt(isd:ied,JsdB:JedB,nz)) ; CS%v_accel_bt(:,:,:) = 0.0
  ALLOC_(CS%PFu_Stokes(IsdB:IedB,jsd:jed,nz)) ; CS%PFu_Stokes(:,:,:) = 0.0
  ALLOC_(CS%PFv_Stokes(isd:ied,JsdB:JedB,nz)) ; CS%PFv_Stokes(:,:,:) = 0.0

  MIS%diffu      => CS%diffu
  MIS%diffv      => CS%diffv
  MIS%PFu        => CS%PFu
  MIS%PFv        => CS%PFv
  MIS%CAu        => CS%CAu
  MIS%CAv        => CS%CAv
  MIS%pbce       => CS%pbce
  MIS%u_accel_bt => CS%u_accel_bt
  MIS%v_accel_bt => CS%v_accel_bt
  MIS%u_av       => CS%u_av
  MIS%v_av       => CS%v_av

  CS%ADp           => Accel_diag
  CS%CDp           => Cont_diag
  Accel_diag%diffu => CS%diffu
  Accel_diag%diffv => CS%diffv
  Accel_diag%PFu   => CS%PFu
  Accel_diag%PFv   => CS%PFv
  Accel_diag%CAu   => CS%CAu
  Accel_diag%CAv   => CS%CAv
  Accel_diag%u_accel_bt => CS%u_accel_bt
  Accel_diag%v_accel_bt => CS%v_accel_bt

  allocate(CS%AD_pred)
  CS%AD_pred%diffu => CS%diffu
  CS%AD_pred%diffv => CS%diffv
  CS%AD_pred%PFu   => CS%PFu
  CS%AD_pred%PFv   => CS%PFv
  CS%AD_pred%CAu   => CS%CAu_pred
  CS%AD_pred%CAv   => CS%CAv_pred
  CS%AD_pred%u_accel_bt => CS%u_accel_bt
  CS%AD_pred%v_accel_bt => CS%v_accel_bt

!  Accel_diag%pbce => CS%pbce
!  Accel_diag%u_accel_bt => CS%u_accel_bt ; Accel_diag%v_accel_bt => CS%v_accel_bt
!  Accel_diag%u_av => CS%u_av ; Accel_diag%v_av => CS%v_av

  id_clock_pass_init = cpu_clock_id('(Ocean init message passing)', grain=CLOCK_ROUTINE)

  call continuity_init(Time, G, GV, US, param_file, diag, CS%continuity_CSp)
  cont_stencil = continuity_stencil(CS%continuity_CSp)
  call CoriolisAdv_init(Time, G, GV, US, param_file, diag, CS%ADp, CS%CoriolisAdv)
  if (CS%calculate_SAL) call SAL_init(h, tv, G, GV, US, param_file, CS%SAL_CSp, restart_CS)
  if (CS%use_tides) then
    call tidal_forcing_init(Time, G, US, param_file, CS%tides_CSp, CS%HA_CSp)
    HA_CSp => CS%HA_CSp
  else
    HA_CSp => NULL()
  endif
  call PressureForce_init(Time, G, GV, US, param_file, diag, CS%PressureForce_CSp, CS%ADp, &
                          CS%SAL_CSp, CS%tides_CSp)
  call hor_visc_init(Time, G, GV, US, param_file, diag, CS%hor_visc, ADp=CS%ADp)
  call vertvisc_init(MIS, Time, G, GV, US, param_file, diag, CS%ADp, dirs, &
                     ntrunc, CS%vertvisc_CSp, CS%fpmix)
  CS%set_visc_CSp => set_visc
  call updateCFLtruncationValue(Time, CS%vertvisc_CSp, US, activate=is_new_run(restart_CS) )

  if (associated(ALE_CSp)) CS%ALE_CSp => ALE_CSp
  if (associated(OBC)) then
    CS%OBC => OBC
    if (OBC%ramp) call update_OBC_ramp(Time, CS%OBC, US, activate=is_new_run(restart_CS) )
  endif
  if (associated(update_OBC_CSp)) CS%update_OBC_CSp => update_OBC_CSp

  eta_rest_name = "sfc" ; if (.not.GV%Boussinesq) eta_rest_name = "p_bot"
  if (.not. query_initialized(CS%eta, trim(eta_rest_name), restart_CS)) then
    ! Estimate eta based on the layer thicknesses - h.  With the Boussinesq
    ! approximation, eta is the free surface height anomaly, while without it
    ! eta is the mass of ocean per unit area.  eta always has the same
    ! dimensions as h, either m or kg m-3.
    !   CS%eta(:,:) = 0.0 already from initialization.
    if (GV%Boussinesq) then
      do j=js,je ; do i=is,ie ; CS%eta(i,j) = -GV%Z_to_H * G%bathyT(i,j) ; enddo ; enddo
    endif
    do k=1,nz ; do j=js,je ; do i=is,ie
      CS%eta(i,j) = CS%eta(i,j) + h(i,j,k)
    enddo ; enddo ; enddo
    call set_initialized(CS%eta, trim(eta_rest_name), restart_CS)
  endif
  ! Copy eta into an output array.
  do j=js,je ; do i=is,ie ; eta(i,j) = CS%eta(i,j) ; enddo ; enddo

  call barotropic_init(u, v, h, Time, G, GV, US, param_file, diag, &
                       CS%barotropic_CSp, restart_CS, calc_dtbt, CS%BT_cont, &
                       CS%OBC, CS%SAL_CSp, HA_CSp)

  if (.not. query_initialized(CS%diffu, "diffu", restart_CS) .or. &
      .not. query_initialized(CS%diffv, "diffv", restart_CS)) then
    call horizontal_viscosity(u, v, h, uh, vh, CS%diffu, CS%diffv, MEKE, VarMix, G, GV, US, CS%hor_visc, &
                              tv, dt, OBC=CS%OBC, BT=CS%barotropic_CSp, TD=thickness_diffuse_CSp, &
                              hu_cont=CS%BT_cont%h_u, hv_cont=CS%BT_cont%h_v)
    call set_initialized(CS%diffu, "diffu", restart_CS)
    call set_initialized(CS%diffv, "diffv", restart_CS)
  endif

  if (.not. query_initialized(CS%u_av, "u2", restart_CS) .or. &
      .not. query_initialized(CS%v_av, "v2", restart_CS)) then
    do k=1,nz ; do j=jsd,jed ; do I=IsdB,IedB ; CS%u_av(I,j,k) = u(I,j,k) ; enddo ; enddo ; enddo
    do k=1,nz ; do J=JsdB,JedB ; do i=isd,ied ; CS%v_av(i,J,k) = v(i,J,k) ; enddo ; enddo ; enddo
    call set_initialized(CS%u_av, "u2", restart_CS)
    call set_initialized(CS%v_av, "v2", restart_CS)
  endif

  if (CS%store_CAu) then
    if (query_initialized(CS%CAu_pred, "CAu", restart_CS) .and. &
        query_initialized(CS%CAv_pred, "CAv", restart_CS)) then
      CS%CAu_pred_stored = .true.
    else
      call only_read_from_restarts(uh, vh, 'uh', 'vh', G, restart_CS, stagger=CGRID_NE, &
                   filename=dirs%input_filename, directory=dirs%restart_input_dir, &
                   success=read_uv, scale=US%m_to_L**2*US%T_to_s/GV%H_to_mks)
      call only_read_from_restarts('h2', CS%h_av, G, restart_CS, &
                   filename=dirs%input_filename, directory=dirs%restart_input_dir, &
                   success=read_h2, scale=1.0/GV%H_to_mks)
      if (read_uv .and. read_h2) then
        call pass_var(CS%h_av, G%Domain, clock=id_clock_pass_init)
      else
        do k=1,nz ; do j=jsd,jed ; do i=isd,ied ; h_tmp(i,j,k) = h(i,j,k) ; enddo ; enddo ; enddo
        call continuity(CS%u_av, CS%v_av, h, h_tmp, uh, vh, dt, G, GV, US, CS%continuity_CSp, CS%OBC, pbv)
        call pass_var(h_tmp, G%Domain, clock=id_clock_pass_init)
        do k=1,nz ; do j=jsd,jed ; do i=isd,ied
          CS%h_av(i,j,k) = 0.5*(h(i,j,k) + h_tmp(i,j,k))
        enddo ; enddo ; enddo
      endif
      call pass_vector(CS%u_av, CS%v_av, G%Domain, halo=2, clock=id_clock_pass_init, complete=.false.)
      call pass_vector(uh, vh, G%Domain, halo=2, clock=id_clock_pass_init, complete=.true.)
      call CorAdCalc(CS%u_av, CS%v_av, CS%h_av, uh, vh, CS%CAu_pred, CS%CAv_pred, CS%OBC, CS%ADp, &
                     G, GV, US, CS%CoriolisAdv, pbv) !, Waves=Waves)
      CS%CAu_pred_stored = .true.
    endif
  else
    CS%CAu_pred_stored = .false.
    ! This call is just here to initialize uh and vh.
    if (.not. query_initialized(uh, "uh", restart_CS) .or. &
        .not. query_initialized(vh, "vh", restart_CS)) then
      do k=1,nz ; do j=jsd,jed ; do i=isd,ied ; h_tmp(i,j,k) = h(i,j,k) ; enddo ; enddo ; enddo
      call continuity(u, v, h, h_tmp, uh, vh, dt, G, GV, US, CS%continuity_CSp, CS%OBC, pbv)
      call pass_var(h_tmp, G%Domain, clock=id_clock_pass_init)
      do k=1,nz ; do j=jsd,jed ; do i=isd,ied
        CS%h_av(i,j,k) = 0.5*(h(i,j,k) + h_tmp(i,j,k))
      enddo ; enddo ; enddo
      call set_initialized(uh, "uh", restart_CS)
      call set_initialized(vh, "vh", restart_CS)
      call set_initialized(CS%h_av, "h2", restart_CS)
      ! Try reading the CAu and CAv fields from the restart file, in case this restart file is
      ! using a newer format.
      call only_read_from_restarts(CS%CAu_pred, CS%CAv_pred, "CAu", "CAv", G, restart_CS, &
                   stagger=CGRID_NE, filename=dirs%input_filename, directory=dirs%restart_input_dir, &
                   success=read_uv, scale=US%m_s_to_L_T*US%T_to_s)
      CS%CAu_pred_stored = read_uv
    else
      if (.not. query_initialized(CS%h_av, "h2", restart_CS)) then
        CS%h_av(:,:,:) = h(:,:,:)
        call set_initialized(CS%h_av, "h2", restart_CS)
      endif
    endif
  endif
  call cpu_clock_begin(id_clock_pass_init)
  call create_group_pass(pass_av_h_uvh, CS%u_av, CS%v_av, G%Domain, halo=2)
  if (CS%CAu_pred_stored) then
    call create_group_pass(pass_av_h_uvh, CS%CAu_pred, CS%CAv_pred, G%Domain, halo=2)
  else
    call create_group_pass(pass_av_h_uvh, CS%h_av, G%Domain, halo=2)
    call create_group_pass(pass_av_h_uvh, uh, vh, G%Domain, halo=2)
  endif
  call do_group_pass(pass_av_h_uvh, G%Domain)
  call cpu_clock_end(id_clock_pass_init)

  flux_units = get_flux_units(GV)
  thickness_units = get_thickness_units(GV)
  CS%id_uh = register_diag_field('ocean_model', 'uh', diag%axesCuL, Time, &
      'Zonal Thickness Flux', flux_units, conversion=GV%H_to_MKS*US%L_to_m**2*US%s_to_T, &
      y_cell_method='sum', v_extensive=.true.)
  CS%id_vh = register_diag_field('ocean_model', 'vh', diag%axesCvL, Time, &
      'Meridional Thickness Flux', flux_units, conversion=GV%H_to_MKS*US%L_to_m**2*US%s_to_T, &
      x_cell_method='sum', v_extensive=.true.)

  CS%id_CAu = register_diag_field('ocean_model', 'CAu', diag%axesCuL, Time, &
      'Zonal Coriolis and Advective Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_CAv = register_diag_field('ocean_model', 'CAv', diag%axesCvL, Time, &
      'Meridional Coriolis and Advective Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_PFu = register_diag_field('ocean_model', 'PFu', diag%axesCuL, Time, &
      'Zonal Pressure Force Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_PFv = register_diag_field('ocean_model', 'PFv', diag%axesCvL, Time, &
      'Meridional Pressure Force Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_ueffA = register_diag_field('ocean_model', 'ueffA', diag%axesCuL, Time, &
       'Effective U-Face Area', 'm^2', conversion=GV%H_to_m*US%L_to_m, &
       y_cell_method='sum', v_extensive=.true.)
  CS%id_veffA = register_diag_field('ocean_model', 'veffA', diag%axesCvL, Time, &
       'Effective V-Face Area', 'm^2', conversion=GV%H_to_m*US%L_to_m, &
       x_cell_method='sum', v_extensive=.true.)
  if (GV%Boussinesq) then
    CS%id_deta_dt = register_diag_field('ocean_model', 'deta_dt', diag%axesT1, Time, &
      'Barotropic SSH tendency due to dynamics', trim(thickness_units)//' s-1', conversion=GV%H_to_MKS*US%s_to_T)
  else
    CS%id_deta_dt = register_diag_field('ocean_model', 'deta_dt', diag%axesT1, Time, &
      'Barotropic column-mass tendency due to dynamics', trim(thickness_units)//' s-1', &
      conversion=GV%H_to_mks*US%s_to_T)
  endif

  !CS%id_hf_PFu = register_diag_field('ocean_model', 'hf_PFu', diag%axesCuL, Time, &
  !    'Fractional Thickness-weighted Zonal Pressure Force Acceleration', &
  !    'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if (CS%id_hf_PFu > 0) call safe_alloc_ptr(CS%ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)

  !CS%id_hf_PFv = register_diag_field('ocean_model', 'hf_PFv', diag%axesCvL, Time, &
  !    'Fractional Thickness-weighted Meridional Pressure Force Acceleration', &
  !    'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if (CS%id_hf_PFv > 0) call safe_alloc_ptr(CS%ADp%diag_hfrac_v,isd,ied,JsdB,JedB,nz)

  !CS%id_hf_CAu = register_diag_field('ocean_model', 'hf_CAu', diag%axesCuL, Time, &
  !    'Fractional Thickness-weighted Zonal Coriolis and Advective Acceleration', &
  !    'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if (CS%id_hf_CAu > 0) call safe_alloc_ptr(CS%ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)

  !CS%id_hf_CAv = register_diag_field('ocean_model', 'hf_CAv', diag%axesCvL, Time, &
  !    'Fractional Thickness-weighted Meridional Coriolis and Advective Acceleration', &
  !    'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if (CS%id_hf_CAv > 0) call safe_alloc_ptr(CS%ADp%diag_hfrac_v,isd,ied,JsdB,JedB,nz)

  CS%id_hf_PFu_2d = register_diag_field('ocean_model', 'hf_PFu_2d', diag%axesCu1, Time, &
      'Depth-sum Fractional Thickness-weighted Zonal Pressure Force Acceleration', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_hf_PFu_2d > 0) call safe_alloc_ptr(CS%ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)

  CS%id_hf_PFv_2d = register_diag_field('ocean_model', 'hf_PFv_2d', diag%axesCv1, Time, &
      'Depth-sum Fractional Thickness-weighted Meridional Pressure Force Acceleration', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_hf_PFv_2d > 0) call safe_alloc_ptr(CS%ADp%diag_hfrac_v,isd,ied,JsdB,JedB,nz)

  CS%id_h_PFu = register_diag_field('ocean_model', 'h_PFu', diag%axesCuL, Time, &
      'Thickness Multiplied Zonal Pressure Force Acceleration', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_PFu > 0) call safe_alloc_ptr(CS%ADp%diag_hu,IsdB,IedB,jsd,jed,nz)

  CS%id_h_PFv = register_diag_field('ocean_model', 'h_PFv', diag%axesCvL, Time, &
      'Thickness Multiplied Meridional Pressure Force Acceleration', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_PFv > 0) call safe_alloc_ptr(CS%ADp%diag_hv,isd,ied,JsdB,JedB,nz)

  CS%id_intz_PFu_2d = register_diag_field('ocean_model', 'intz_PFu_2d', diag%axesCu1, Time, &
      'Depth-integral of Zonal Pressure Force Acceleration', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_intz_PFu_2d > 0) call safe_alloc_ptr(CS%ADp%diag_hu,IsdB,IedB,jsd,jed,nz)

  CS%id_intz_PFv_2d = register_diag_field('ocean_model', 'intz_PFv_2d', diag%axesCv1, Time, &
      'Depth-integral of Meridional Pressure Force Acceleration', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_intz_PFv_2d > 0) call safe_alloc_ptr(CS%ADp%diag_hv,isd,ied,JsdB,JedB,nz)

  CS%id_hf_CAu_2d = register_diag_field('ocean_model', 'hf_CAu_2d', diag%axesCu1, Time, &
      'Depth-sum Fractional Thickness-weighted Zonal Coriolis and Advective Acceleration', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_hf_CAu_2d > 0) call safe_alloc_ptr(CS%ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)

  CS%id_hf_CAv_2d = register_diag_field('ocean_model', 'hf_CAv_2d', diag%axesCv1, Time, &
      'Depth-sum Fractional Thickness-weighted Meridional Coriolis and Advective Acceleration', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_hf_CAv_2d > 0) call safe_alloc_ptr(CS%ADp%diag_hfrac_v,isd,ied,JsdB,JedB,nz)

  CS%id_h_CAu = register_diag_field('ocean_model', 'h_CAu', diag%axesCuL, Time, &
      'Thickness Multiplied Zonal Coriolis and Advective Acceleration', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_CAu > 0) call safe_alloc_ptr(CS%ADp%diag_hu,IsdB,IedB,jsd,jed,nz)

  CS%id_h_CAv = register_diag_field('ocean_model', 'h_CAv', diag%axesCvL, Time, &
      'Thickness Multiplied Meridional Coriolis and Advective Acceleration', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_CAv > 0) call safe_alloc_ptr(CS%ADp%diag_hv,isd,ied,JsdB,JedB,nz)

  CS%id_intz_CAu_2d = register_diag_field('ocean_model', 'intz_CAu_2d', diag%axesCu1, Time, &
      'Depth-integral of Zonal Coriolis and Advective Acceleration', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_intz_CAu_2d > 0) call safe_alloc_ptr(CS%ADp%diag_hu,IsdB,IedB,jsd,jed,nz)

  CS%id_intz_CAv_2d = register_diag_field('ocean_model', 'intz_CAv_2d', diag%axesCv1, Time, &
      'Depth-integral of Meridional Coriolis and Advective Acceleration', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_intz_CAv_2d > 0) call safe_alloc_ptr(CS%ADp%diag_hv,isd,ied,JsdB,JedB,nz)

  CS%id_uav = register_diag_field('ocean_model', 'uav', diag%axesCuL, Time, &
      'Barotropic-step Averaged Zonal Velocity', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vav = register_diag_field('ocean_model', 'vav', diag%axesCvL, Time, &
      'Barotropic-step Averaged Meridional Velocity', 'm s-1', conversion=US%L_T_to_m_s)

  CS%id_u_BT_accel = register_diag_field('ocean_model', 'u_BT_accel', diag%axesCuL, Time, &
    'Barotropic Anomaly Zonal Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_v_BT_accel = register_diag_field('ocean_model', 'v_BT_accel', diag%axesCvL, Time, &
    'Barotropic Anomaly Meridional Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)

  !CS%id_hf_u_BT_accel = register_diag_field('ocean_model', 'hf_u_BT_accel', diag%axesCuL, Time, &
  !    'Fractional Thickness-weighted Barotropic Anomaly Zonal Acceleration', &
  !    'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if (CS%id_hf_u_BT_accel > 0) call safe_alloc_ptr(CS%ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)

  !CS%id_hf_v_BT_accel = register_diag_field('ocean_model', 'hf_v_BT_accel', diag%axesCvL, Time, &
  !    'Fractional Thickness-weighted Barotropic Anomaly Meridional Acceleration', &
  !    'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if (CS%id_hf_v_BT_accel > 0) call safe_alloc_ptr(CS%ADp%diag_hfrac_v,isd,ied,JsdB,JedB,nz)

  CS%id_hf_u_BT_accel_2d = register_diag_field('ocean_model', 'hf_u_BT_accel_2d', diag%axesCu1, Time, &
      'Depth-sum Fractional Thickness-weighted Barotropic Anomaly Zonal Acceleration', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_hf_u_BT_accel_2d > 0) call safe_alloc_ptr(CS%ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)

  CS%id_hf_v_BT_accel_2d = register_diag_field('ocean_model', 'hf_v_BT_accel_2d', diag%axesCv1, Time, &
      'Depth-sum Fractional Thickness-weighted Barotropic Anomaly Meridional Acceleration', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_hf_v_BT_accel_2d > 0) call safe_alloc_ptr(CS%ADp%diag_hfrac_v,isd,ied,JsdB,JedB,nz)

  CS%id_h_u_BT_accel = register_diag_field('ocean_model', 'h_u_BT_accel', diag%axesCuL, Time, &
      'Thickness Multiplied Barotropic Anomaly Zonal Acceleration', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_u_BT_accel > 0) call safe_alloc_ptr(CS%ADp%diag_hu,IsdB,IedB,jsd,jed,nz)

  CS%id_h_v_BT_accel = register_diag_field('ocean_model', 'h_v_BT_accel', diag%axesCvL, Time, &
      'Thickness Multiplied Barotropic Anomaly Meridional Acceleration', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_v_BT_accel > 0) call safe_alloc_ptr(CS%ADp%diag_hv,isd,ied,JsdB,JedB,nz)

  CS%id_intz_u_BT_accel_2d = register_diag_field('ocean_model', 'intz_u_BT_accel_2d', diag%axesCu1, Time, &
      'Depth-integral of Barotropic Anomaly Zonal Acceleration', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_intz_u_BT_accel_2d > 0) call safe_alloc_ptr(CS%ADp%diag_hu,IsdB,IedB,jsd,jed,nz)

  CS%id_intz_v_BT_accel_2d = register_diag_field('ocean_model', 'intz_v_BT_accel_2d', diag%axesCv1, Time, &
      'Depth-integral of Barotropic Anomaly Meridional Acceleration', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_intz_v_BT_accel_2d > 0) call safe_alloc_ptr(CS%ADp%diag_hv,isd,ied,JsdB,JedB,nz)

  CS%id_PFu_visc_rem = register_diag_field('ocean_model', 'PFu_visc_rem', diag%axesCuL, Time, &
      'Zonal Pressure Force Acceleration multiplied by the viscous remnant', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_PFu_visc_rem > 0) call safe_alloc_ptr(CS%ADp%visc_rem_u,IsdB,IedB,jsd,jed,nz)
  CS%id_PFv_visc_rem = register_diag_field('ocean_model', 'PFv_visc_rem', diag%axesCvL, Time, &
      'Meridional Pressure Force Acceleration multiplied by the viscous remnant', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_PFv_visc_rem > 0) call safe_alloc_ptr(CS%ADp%visc_rem_v,isd,ied,JsdB,JedB,nz)

  CS%id_CAu_visc_rem = register_diag_field('ocean_model', 'CAu_visc_rem', diag%axesCuL, Time, &
      'Zonal Coriolis and Advective Acceleration multiplied by the viscous remnant', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_CAu_visc_rem > 0) call safe_alloc_ptr(CS%ADp%visc_rem_u,IsdB,IedB,jsd,jed,nz)
  CS%id_CAv_visc_rem = register_diag_field('ocean_model', 'CAv_visc_rem', diag%axesCvL, Time, &
      'Meridional Coriolis and Advective Acceleration multiplied by the viscous remnant', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_CAv_visc_rem > 0) call safe_alloc_ptr(CS%ADp%visc_rem_v,isd,ied,JsdB,JedB,nz)

  CS%id_u_BT_accel_visc_rem = register_diag_field('ocean_model', 'u_BT_accel_visc_rem', diag%axesCuL, Time, &
      'Barotropic Anomaly Zonal Acceleration multiplied by the viscous remnant', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_u_BT_accel_visc_rem > 0) call safe_alloc_ptr(CS%ADp%visc_rem_u,IsdB,IedB,jsd,jed,nz)
  CS%id_v_BT_accel_visc_rem = register_diag_field('ocean_model', 'v_BT_accel_visc_rem', diag%axesCvL, Time, &
      'Barotropic Anomaly Meridional Acceleration multiplied by the viscous remnant', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_v_BT_accel_visc_rem > 0) call safe_alloc_ptr(CS%ADp%visc_rem_v,isd,ied,JsdB,JedB,nz)

  id_clock_Cor        = cpu_clock_id('(Ocean Coriolis & mom advection)', grain=CLOCK_MODULE)
  id_clock_continuity = cpu_clock_id('(Ocean continuity equation)',      grain=CLOCK_MODULE)
  id_clock_pres       = cpu_clock_id('(Ocean pressure force)',           grain=CLOCK_MODULE)
  id_clock_vertvisc   = cpu_clock_id('(Ocean vertical viscosity)',       grain=CLOCK_MODULE)
  id_clock_horvisc    = cpu_clock_id('(Ocean horizontal viscosity)',     grain=CLOCK_MODULE)
  id_clock_mom_update = cpu_clock_id('(Ocean momentum increments)',      grain=CLOCK_MODULE)
  id_clock_pass       = cpu_clock_id('(Ocean message passing)',          grain=CLOCK_MODULE)
  id_clock_btcalc     = cpu_clock_id('(Ocean barotropic mode calc)',     grain=CLOCK_MODULE)
  id_clock_btstep     = cpu_clock_id('(Ocean barotropic mode stepping)', grain=CLOCK_MODULE)
  id_clock_btforce    = cpu_clock_id('(Ocean barotropic forcing calc)',  grain=CLOCK_MODULE)

end subroutine initialize_dyn_split_RK2


!> Close the dyn_split_RK2 module
subroutine end_dyn_split_RK2(CS)
  type(MOM_dyn_split_RK2_CS), pointer :: CS  !< module control structure

  call barotropic_end(CS%barotropic_CSp)

  call vertvisc_end(CS%vertvisc_CSp)
  deallocate(CS%vertvisc_CSp)

  call hor_visc_end(CS%hor_visc)
  if (CS%calculate_SAL) call SAL_end(CS%SAL_CSp)
  if (CS%use_tides) call tidal_forcing_end(CS%tides_CSp)
  call CoriolisAdv_end(CS%CoriolisAdv)

  DEALLOC_(CS%diffu) ; DEALLOC_(CS%diffv)
  DEALLOC_(CS%CAu)   ; DEALLOC_(CS%CAv)
  DEALLOC_(CS%CAu_pred) ; DEALLOC_(CS%CAv_pred)
  DEALLOC_(CS%PFu)   ; DEALLOC_(CS%PFv)

  if (associated(CS%taux_bot)) deallocate(CS%taux_bot)
  if (associated(CS%tauy_bot)) deallocate(CS%tauy_bot)
  DEALLOC_(CS%uhbt) ; DEALLOC_(CS%vhbt)
  DEALLOC_(CS%u_accel_bt) ; DEALLOC_(CS%v_accel_bt)
  DEALLOC_(CS%visc_rem_u) ; DEALLOC_(CS%visc_rem_v)

  DEALLOC_(CS%eta) ; DEALLOC_(CS%eta_PF) ; DEALLOC_(CS%pbce)
  DEALLOC_(CS%h_av) ; DEALLOC_(CS%u_av) ; DEALLOC_(CS%v_av)

  call dealloc_BT_cont_type(CS%BT_cont)
  deallocate(CS%AD_pred)

  deallocate(CS)
end subroutine end_dyn_split_RK2


!> \namespace mom_dynamics_split_rk2
!!
!!  This file time steps the adiabatic dynamic core by splitting
!!  between baroclinic and barotropic modes. It uses a pseudo-second order
!!  Runge-Kutta time stepping scheme for the baroclinic momentum
!!  equation and a forward-backward coupling between the baroclinic
!!  momentum and continuity equations.  This split time-stepping
!!  scheme is described in detail in Hallberg (JCP, 1997). Additional
!!  issues related to exact tracer conservation and how to
!!  ensure consistency between the barotropic and layered estimates
!!  of the free surface height are described in Hallberg and
!!  Adcroft (Ocean Modelling, 2009).  This was the time stepping code
!!  that is used for most GOLD applications, including GFDL's ESM2G
!!  Earth system model, and all of the examples provided with the
!!  MOM code (although several of these solutions are routinely
!!  verified by comparison with the slower unsplit schemes).
!!
!!  The subroutine step_MOM_dyn_split_RK2 actually does the time
!!  stepping, while register_restarts_dyn_split_RK2 sets the fields
!!  that are found in a full restart file with this scheme, and
!!  initialize_dyn_split_RK2 initializes the cpu clocks that are
!!  used in this module.  For largely historical reasons, this module
!!  does not have its own control structure, but shares the same
!!  control structure with MOM.F90 and the other MOM_dynamics_...
!!  modules.

end module MOM_dynamics_split_RK2
