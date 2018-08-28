!> This is the main routine for MOM
module MOM

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_variables, only : vertvisc_type
use MOM_open_boundary, only : ocean_OBC_type

! A Structure with pointers to forcing fields to drive MOM;
! all fluxes are positive downward.
use MOM_forcing_type, only : forcing, mech_forcing

use MOM_variables, only: accel_diag_ptrs, cont_diag_ptrs, ocean_internal_state

! A structure containing pointers to various fields
! to describe surface state, and will be returned
! to the calling program.
use MOM_variables, only : surface

! A structure containing pointers to an assortment of
! thermodynamic fields, including potential/Conservative
! temperature, salinity and mixed layer density.
use MOM_variables, only: thermo_var_ptrs

! Infrastructure modules
use MOM_debugging,            only : MOM_debugging_init, hchksum, uvchksum
use MOM_checksum_packages,    only : MOM_thermo_chksum, MOM_state_chksum, MOM_accel_chksum
use MOM_cpu_clock,            only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,            only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT
use MOM_cpu_clock,            only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_coms,                 only : reproducing_sum
use MOM_coord_initialization, only : MOM_initialize_coord
use MOM_diag_mediator,        only : diag_mediator_init, enable_averaging
use MOM_diag_mediator,        only : diag_mediator_infrastructure_init
use MOM_diag_mediator,        only : diag_register_area_ids
use MOM_diag_mediator,        only : diag_associate_volume_cell_measure
use MOM_diag_mediator,        only : diag_set_state_ptrs, diag_update_remap_grids
use MOM_diag_mediator,        only : disable_averaging, post_data, safe_alloc_ptr
use MOM_diag_mediator,        only : register_diag_field, register_static_field
use MOM_diag_mediator,        only : register_scalar_field, get_diag_time_end
use MOM_diag_mediator,        only : set_axes_info, diag_ctrl, diag_masks_set
use MOM_diag_mediator,        only : set_masks_for_axes
use MOM_domains,              only : MOM_domains_init, clone_MOM_domain
use MOM_domains,              only : sum_across_PEs, pass_var, pass_vector
use MOM_domains,              only : To_North, To_East, To_South, To_West
use MOM_domains,              only : To_All, Omit_corners, CGRID_NE, SCALAR_PAIR
use MOM_domains,              only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,              only : start_group_pass, complete_group_pass, Omit_Corners
use MOM_error_handler,        only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler,        only : MOM_set_verbosity, callTree_showQuery
use MOM_error_handler,        only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,          only : read_param, get_param, log_version, param_file_type
use MOM_fixed_initialization, only : MOM_initialize_fixed
use MOM_forcing_type,         only : MOM_forcing_chksum
use MOM_get_input,            only : Get_MOM_Input, directories
use MOM_io,                   only : MOM_io_init, vardesc, var_desc
use MOM_io,                   only : slasher, file_exists, read_data
use MOM_obsolete_params,      only : find_obsolete_params
use MOM_restart,              only : register_restart_field, query_initialized, save_restart
use MOM_restart,              only : restart_init, is_new_run, MOM_restart_CS
use MOM_spatial_means,        only : global_area_mean, global_area_integral
use MOM_state_initialization, only : MOM_initialize_state
use MOM_time_manager,         only : time_type, set_time, time_type_to_real, operator(+)
use MOM_time_manager,         only : operator(-), operator(>), operator(*), operator(/)
use MOM_time_manager,         only : increment_date
use MOM_unit_tests,           only : unit_tests
use coupler_types_mod,        only : coupler_type_send_data, coupler_1d_bc_type, coupler_type_spawn

! MOM core modules
use MOM_ALE,                   only : ALE_init, ALE_end, ALE_main, ALE_CS, adjustGridForIntegrity
use MOM_ALE,                   only : ALE_getCoordinate, ALE_getCoordinateUnits, ALE_writeCoordinateFile
use MOM_ALE,                   only : ALE_updateVerticalGridType, ALE_remap_init_conds, ALE_register_diags
use MOM_boundary_update,       only : call_OBC_register, OBC_register_end, update_OBC_CS
use MOM_continuity,            only : continuity, continuity_init, continuity_CS, continuity_stencil
use MOM_CoriolisAdv,           only : CorAdCalc, CoriolisAdv_init, CoriolisAdv_CS
use MOM_diabatic_driver,       only : diabatic, diabatic_driver_init, diabatic_CS
use MOM_diabatic_driver,       only : adiabatic, adiabatic_driver_init, diabatic_driver_end
use MOM_diagnostics,           only : calculate_diagnostic_fields, MOM_diagnostics_init
use MOM_diagnostics,           only : diagnostics_CS
use MOM_diag_to_Z,             only : calculate_Z_diag_fields, calculate_Z_transport
use MOM_diag_to_Z,             only : MOM_diag_to_Z_init, register_Z_tracer, diag_to_Z_CS
use MOM_diag_to_Z,             only : MOM_diag_to_Z_end
use MOM_dynamics_unsplit,      only : step_MOM_dyn_unsplit, register_restarts_dyn_unsplit
use MOM_dynamics_unsplit,      only : initialize_dyn_unsplit, end_dyn_unsplit
use MOM_dynamics_unsplit,      only : MOM_dyn_unsplit_CS
use MOM_dynamics_split_RK2,    only : step_MOM_dyn_split_RK2, register_restarts_dyn_split_RK2
use MOM_dynamics_split_RK2,    only : initialize_dyn_split_RK2, end_dyn_split_RK2
use MOM_dynamics_split_RK2,    only : MOM_dyn_split_RK2_CS
use MOM_dynamics_unsplit_RK2,  only : step_MOM_dyn_unsplit_RK2, register_restarts_dyn_unsplit_RK2
use MOM_dynamics_unsplit_RK2,  only : initialize_dyn_unsplit_RK2, end_dyn_unsplit_RK2
use MOM_dynamics_unsplit_RK2,  only : MOM_dyn_unsplit_RK2_CS
use MOM_dyn_horgrid,           only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_EOS,                   only : EOS_init
use MOM_EOS,                   only : gsw_sp_from_sr, gsw_pt_from_ct
use MOM_EOS,                   only : calculate_density
use MOM_debugging,             only : check_redundant
use MOM_grid,                  only : ocean_grid_type, set_first_direction
use MOM_grid,                  only : MOM_grid_init, MOM_grid_end
use MOM_hor_index,             only : hor_index_type, hor_index_init
use MOM_hor_visc,              only : horizontal_viscosity, hor_visc_init
use MOM_interface_heights,     only : find_eta
use MOM_lateral_mixing_coeffs, only : calc_slope_functions, VarMix_init
use MOM_lateral_mixing_coeffs, only : calc_resoln_function, VarMix_CS
use MOM_MEKE,                  only : MEKE_init, MEKE_alloc_register_restart, step_forward_MEKE, MEKE_CS
use MOM_MEKE_types,            only : MEKE_type
use MOM_mixed_layer_restrat,   only : mixedlayer_restrat, mixedlayer_restrat_init, mixedlayer_restrat_CS
use MOM_mixed_layer_restrat,   only : mixedlayer_restrat_register_restarts
use MOM_neutral_diffusion,     only : neutral_diffusion_CS, neutral_diffusion_diag_init
use MOM_obsolete_diagnostics,  only : register_obsolete_diagnostics
use MOM_open_boundary,         only : OBC_registry_type, register_temp_salt_segments
use MOM_PressureForce,         only : PressureForce, PressureForce_init, PressureForce_CS
use MOM_set_visc,              only : set_viscous_BBL, set_viscous_ML, set_visc_init
use MOM_set_visc,              only : set_visc_register_restarts, set_visc_CS
use MOM_sponge,                only : init_sponge_diags, sponge_CS
use MOM_ALE_sponge,            only : init_ALE_sponge_diags, ALE_sponge_CS
use MOM_thickness_diffuse,     only : thickness_diffuse, thickness_diffuse_init, thickness_diffuse_CS
use MOM_tidal_forcing,         only : tidal_forcing_init, tidal_forcing_CS
use MOM_tracer_advect,         only : advect_tracer, tracer_advect_init
use MOM_tracer_advect,         only : tracer_advect_end, tracer_advect_CS
use MOM_tracer_hor_diff,       only : tracer_hordiff, tracer_hor_diff_init
use MOM_tracer_hor_diff,       only : tracer_hor_diff_end, tracer_hor_diff_CS
use MOM_tracer_registry,       only : register_tracer, tracer_registry_init
use MOM_tracer_registry,       only : add_tracer_diagnostics, tracer_registry_type
use MOM_tracer_registry,       only : lock_tracer_registry, tracer_registry_end
use MOM_tracer_flow_control,   only : call_tracer_register, tracer_flow_control_CS
use MOM_tracer_flow_control,   only : tracer_flow_control_init, call_tracer_surface_state
use MOM_tracer_flow_control,   only : tracer_flow_control_end
use MOM_transcribe_grid,       only : copy_dyngrid_to_MOM_grid, copy_MOM_grid_to_dyngrid
use MOM_vert_friction,         only : vertvisc, vertvisc_remnant
use MOM_vert_friction,         only : vertvisc_limit_vel, vertvisc_init
use MOM_verticalGrid,          only : verticalGrid_type, verticalGridInit, verticalGridEnd
use MOM_verticalGrid,          only : get_thickness_units, get_flux_units, get_tr_flux_units
! Offline modules
use MOM_offline_main,          only : offline_transport_CS, offline_transport_init, update_offline_fields
use MOM_offline_main,          only : insert_offline_main, extract_offline_main, post_offline_convergence_diags
use MOM_offline_main,          only : register_diags_offline_transport, offline_advection_ale
use MOM_offline_main,          only : offline_redistribute_residual, offline_diabatic_ale
use MOM_offline_main,          only : offline_fw_fluxes_into_ocean, offline_fw_fluxes_out_ocean
use MOM_offline_main,          only : offline_advection_layer, offline_transport_end
use MOM_ALE,                   only : ale_offline_tracer_final, ALE_main_offline


implicit none ; private

#include <MOM_memory.h>

!> Control structure for this module
type, public :: MOM_control_struct
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: &
    h, &      !< layer thickness (m or kg/m2 (H))
    T, &      !< potential temperature (degrees C)
    S         !< salinity (ppt)
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
    u,  &     !< zonal velocity component (m/s)
    uh, &     !< uh = u * h * dy at u grid points (m3/s or kg/s)
    uhtr      !< accumulated zonal thickness fluxes to advect tracers (m3 or kg)
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
    v,  &     !< meridional velocity (m/s)
    vh, &     !< vh = v * h * dx at v grid points (m3/s or kg/s)
    vhtr      !< accumulated meridional thickness fluxes to advect tracers (m3 or kg)
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    ave_ssh   !< time-averaged (ave over baroclinic time steps) sea surface height (meter)
  real, pointer, dimension(:,:)   :: Hml => NULL() !< active mixed layer depth, in m
  real, pointer, dimension(:,:,:) :: &
    u_prev => NULL(), &  !< previous value of u stored for diagnostics
    v_prev => NULL()     !< previous value of v stored for diagnostics

  type(ocean_grid_type) :: G       !< structure containing metrics and grid info
  type(verticalGrid_type), pointer :: GV => NULL() !< structure containing vertical grid info
  type(thermo_var_ptrs) :: tv      !< structure containing pointers to available
                                   !! thermodynamic fields
  type(diag_ctrl)       :: diag    !< structure to regulate diagnostic output timing
  type(vertvisc_type)   :: visc    !< structure containing vertical viscosities,
                                   !! bottom drag viscosities, and related fields
  type(MEKE_type), pointer :: MEKE => NULL()  !<  structure containing fields
                                   !! related to the Mesoscale Eddy Kinetic Energy
  type(accel_diag_ptrs) :: ADp     !< structure containing pointers to accelerations,
                                   !! for derived diagnostics (e.g., energy budgets)
  type(cont_diag_ptrs)  :: CDp     !< structure containing pointers continuity equation
                                   !! terms, for derived diagnostics (e.g., energy budgets)

  logical :: split                   !< If true, use the split time stepping scheme.
  logical :: use_RK2                 !< If true, use RK2 instead of RK3 in unsplit mode
                                     !! (i.e., no split between barotropic and baroclinic).
  logical :: adiabatic               !< If true, then no diapycnal mass fluxes, with no calls
                                     !! to routines to calculate or apply diapycnal fluxes.
  logical :: use_temperature         !< If true, temp and saln used as state variables.
  logical :: calc_rho_for_sea_lev    !< If true, calculate rho to convert pressure to sea level
  logical :: use_frazil              !< If true, liquid seawater freezes if temp below freezing,
                                     !! with accumulated heat deficit returned to surface ocean.
  logical :: bound_salinity          !< If true, salt is added to keep salinity above
                                     !! a minimum value, and the deficit is reported.
  logical :: bulkmixedlayer          !< If true, a refined bulk mixed layer scheme is used
                                     !! with nkml sublayers and nkbl buffer layer.
  logical :: diabatic_first          !< If true, apply diabatic and thermodynamic
                                     !! processes before time stepping the dynamics.
  logical :: use_conT_absS           !< If true, , the prognostics T&S are the conservative temperature
                                     !! and absolute salinity. Care should be taken to convert them
                                     !! to potential temperature and practical salinity before
                                     !! exchanging them with the coupler and/or reporting T&S diagnostics.
  logical :: thickness_diffuse       !< If true, diffuse interface height w/ a diffusivity KHTH.
  logical :: thickness_diffuse_first !< If true, diffuse thickness before dynamics.
  logical :: mixedlayer_restrat      !< If true, use submesoscale mixed layer restratifying scheme.
  logical :: useMEKE                 !< If true, call the MEKE parameterization.
  logical :: debug                   !< If true, write verbose checksums for debugging purposes.
  logical :: debug_truncations       !< If true, turn on diagnostics useful for debugging truncations.
  logical :: use_ALE_algorithm       !< If true, use the ALE algorithm rather than layered
                                     !! isopycnal/stacked shallow water mode. This logical is
                                     !! set by calling the function useRegridding() from the
                                     !! MOM_regridding module.
  logical :: do_dynamics             !< If false, does not call step_MOM_dyn_*. This is an
                                     !! undocumented run-time flag that is fragile.
  logical :: offline_tracer_mode = .false.
                                     !< If true, step_offline() is called instead of step_MOM().
                                     !! This is intended for running MOM6 in offline tracer mode
  logical :: advect_TS               !< If false, then no horizontal advection of temperature
                                     !! and salnity is performed
  real    :: dt                      !< (baroclinic) dynamics time step (seconds)
  real    :: dt_therm                !< thermodynamics time step (seconds)
  logical :: thermo_spans_coupling   !< If true, thermodynamic and tracer time
                                     !! steps can span multiple coupled time steps.
  real    :: t_dyn_rel_adv           !< The time of the dynamics relative to tracer
                                     !! advection and lateral mixing (in seconds), or
                                     !! equivalently the elapsed time since advectively
                                     !! updating the tracers.  t_dyn_rel_adv is invariably
                                     !! positive and may span multiple coupling timesteps.
  real    :: t_dyn_rel_thermo        !< The time of the dynamics relative to diabatic
                                     !! processes and remapping (in seconds).  t_dyn_rel_thermo
                                     !! can be negative or positive depending on whether
                                     !! the diabatic processes are applied before or after
                                     !! the dynamics and may span multiple coupling timesteps.
  real    :: t_dyn_rel_diag          !< The time of the diagnostics relative to diabatic
                                     !! processes and remapping (in seconds).  t_dyn_rel_diag
                                     !! is always positive, since the diagnostics must lag.
  type(time_type) :: Z_diag_interval !< amount of time between calculating Z-space diagnostics
  type(time_type) :: Z_diag_time     !< next time to compute Z-space diagnostics
  type(time_type), pointer :: Time   !< pointer to ocean clock
  real :: rel_time = 0.0             !< relative time (sec) since start of current execution
  real :: dtbt_reset_period          !< The time interval in seconds between dynamic
                                     !! recalculation of the barotropic time step.  If
                                     !! this is negative, it is never calculated, and
                                     !! if it is 0, it is calculated every step.

  logical :: interp_p_surf           !< If true, linearly interpolate surface pressure
                                     !! over the coupling time step, using specified value
                                     !! at the end of the coupling step. False by default.
  logical :: p_surf_prev_set         !< If true, p_surf_prev has been properly set from
                                     !! a previous time-step or the ocean restart file.
                                     !! This is only valid when interp_p_surf is true.

  real :: Hmix                       !< Diagnostic mixed layer thickness (meter) when
                                     !! bulk mixed layer is not used.
  real :: Hmix_UV                    !< Depth scale over which to average surface flow to
                                     !! feedback to the coupler/driver (m).
                                     !! bulk mixed layer is not used.
  real :: missing=-1.0e34            !< missing data value for masked fields

  ! Flags needed to reach between start and finish phases of initialization
  logical :: write_IC                !< If true, then the initial conditions will be written to file
  character(len=120) :: IC_file      !< A file into which the initial conditions are
                                     !! written in a new run if SAVE_INITIAL_CONDS is true.

  integer :: ntrunc                  !< number u,v truncations since last call to write_energy
  logical :: check_bad_surface_vals  !< If true, scan surface state for ridiculous values.
  real    :: bad_val_ssh_max         !< Maximum SSH before triggering bad value message
  real    :: bad_val_sst_max         !< Maximum SST before triggering bad value message
  real    :: bad_val_sst_min         !< Minimum SST before triggering bad value message
  real    :: bad_val_sss_max         !< Maximum SSS before triggering bad value message
  real    :: bad_val_column_thickness!< Minimum column thickness before triggering bad value message

  real, pointer, dimension(:,:) :: &
    p_surf_prev  => NULL(), & !< surface pressure (Pa) at end  previous call to step_MOM
    p_surf_begin => NULL(), & !< surface pressure (Pa) at start of step_MOM_dyn_...
    p_surf_end   => NULL()    !< surface pressure (Pa) at end   of step_MOM_dyn_...

  type(vardesc) :: &
    vd_T, &   !< vardesc array describing potential temperature
    vd_S      !< vardesc array describing salinity

  real, pointer, dimension(:,:,:) :: &  !< diagnostic arrays of advective/diffusive tracer fluxes
    T_adx => NULL(), T_ady => NULL(), T_diffx => NULL(), T_diffy => NULL(), &
    S_adx => NULL(), S_ady => NULL(), S_diffx => NULL(), S_diffy => NULL()

  real, pointer, dimension(:,:) :: & !< diagnostic arrays of vertically integrated advective/diffusive fluxes
    T_adx_2d => NULL(), T_ady_2d => NULL(), T_diffx_2d => NULL(), T_diffy_2d => NULL(), &
    S_adx_2d => NULL(), S_ady_2d => NULL(), S_diffx_2d => NULL(), S_diffy_2d => NULL()

  real, pointer, dimension(:,:,:) :: &  !< diagnostic arrays for advection tendencies and total tendencies
    T_advection_xy => NULL(), S_advection_xy => NULL(),  &
    T_prev         => NULL(), S_prev         => NULL(),  &
    Th_prev        => NULL(), Sh_prev        => NULL()

  real, pointer, dimension(:,:,:) :: & !< diagnostic arrays for variance decay through ALE
    T_squared => NULL(), S_squared => NULL()

  logical :: tendency_diagnostics = .false.

  ! diagnostic ids

  ! 3-d state fields
  integer :: id_u  = -1
  integer :: id_v  = -1
  integer :: id_h  = -1
  integer :: id_T  = -1
  integer :: id_S  = -1
  integer :: id_Tcon  = -1
  integer :: id_Sabs  = -1

  ! 2-d surface and bottom fields
  integer :: id_zos      = -1
  integer :: id_zossq    = -1
  integer :: id_volo     = -1
  integer :: id_ssh      = -1
  integer :: id_ssh_ga   = -1
  integer :: id_sst      = -1
  integer :: id_sst_sq   = -1
  integer :: id_sss      = -1
  integer :: id_sss_sq   = -1
  integer :: id_ssu      = -1
  integer :: id_ssv      = -1
  integer :: id_speed    = -1
  integer :: id_ssh_inst = -1
  integer :: id_tob      = -1
  integer :: id_sob      = -1
  integer :: id_sstcon   = -1
  integer :: id_sssabs   = -1

  ! heat and salt flux fields
  integer :: id_fraz         = -1
  integer :: id_salt_deficit = -1
  integer :: id_Heat_PmE     = -1
  integer :: id_intern_heat  = -1

  ! transport of temperature and salinity
  integer :: id_Tadx      = -1
  integer :: id_Tady      = -1
  integer :: id_Tdiffx    = -1
  integer :: id_Tdiffy    = -1
  integer :: id_Sadx      = -1
  integer :: id_Sady      = -1
  integer :: id_Sdiffx    = -1
  integer :: id_Sdiffy    = -1
  integer :: id_Tadx_2d   = -1
  integer :: id_Tady_2d   = -1
  integer :: id_Tdiffx_2d = -1
  integer :: id_Tdiffy_2d = -1
  integer :: id_Sadx_2d   = -1
  integer :: id_Sady_2d   = -1
  integer :: id_Sdiffx_2d = -1
  integer :: id_Sdiffy_2d = -1

  ! tendencies for temp/heat and saln/salt
  integer :: id_T_advection_xy    = -1
  integer :: id_T_advection_xy_2d = -1
  integer :: id_T_tendency        = -1
  integer :: id_Th_tendency       = -1
  integer :: id_Th_tendency_2d    = -1
  integer :: id_S_advection_xy    = -1
  integer :: id_S_advection_xy_2d = -1
  integer :: id_S_tendency        = -1
  integer :: id_Sh_tendency       = -1
  integer :: id_Sh_tendency_2d    = -1

  ! variance decay for temp and heat
  integer :: id_T_vardec = -1
  integer :: id_S_vardec = -1

  ! diagnostic for fields prior to applying diapycnal physics
  integer :: id_u_predia = -1
  integer :: id_v_predia = -1
  integer :: id_h_predia = -1
  integer :: id_T_predia = -1
  integer :: id_S_predia = -1
  integer :: id_e_predia = -1
  integer :: id_u_preale = -1
  integer :: id_v_preale = -1
  integer :: id_h_preale = -1
  integer :: id_T_preale = -1
  integer :: id_S_preale = -1
  integer :: id_e_preale = -1

  ! Diagnostics for tracer horizontal transport
  integer :: id_uhtr = -1, id_umo = -1, id_umo_2d = 1
  integer :: id_vhtr = -1, id_vmo = -1, id_vmo_2d = 1

  ! The remainder provides pointers to child module control structures.
  type(MOM_dyn_unsplit_CS),      pointer :: dyn_unsplit_CSp      => NULL()
  type(MOM_dyn_unsplit_RK2_CS),  pointer :: dyn_unsplit_RK2_CSp  => NULL()
  type(MOM_dyn_split_RK2_CS),    pointer :: dyn_split_RK2_CSp    => NULL()

  type(set_visc_CS),             pointer :: set_visc_CSp           => NULL()
  type(diabatic_CS),             pointer :: diabatic_CSp           => NULL()
  type(thickness_diffuse_CS),    pointer :: thickness_diffuse_CSp  => NULL()
  type(mixedlayer_restrat_CS),   pointer :: mixedlayer_restrat_CSp => NULL()
  type(MEKE_CS),                 pointer :: MEKE_CSp               => NULL()
  type(VarMix_CS),               pointer :: VarMix                 => NULL()
  type(tracer_registry_type),    pointer :: tracer_Reg             => NULL()
  type(tracer_advect_CS),        pointer :: tracer_adv_CSp         => NULL()
  type(tracer_hor_diff_CS),      pointer :: tracer_diff_CSp        => NULL()
  type(neutral_diffusion_CS),    pointer :: neutral_diffusion_CSp  => NULL()
  type(tracer_flow_control_CS),  pointer :: tracer_flow_CSp        => NULL()
  type(diagnostics_CS),          pointer :: diagnostics_CSp        => NULL()
  type(diag_to_Z_CS),            pointer :: diag_to_Z_CSp          => NULL()
  type(MOM_restart_CS),          pointer :: restart_CSp            => NULL()
  type(update_OBC_CS),           pointer :: update_OBC_CSp         => NULL()
  type(ocean_OBC_type),          pointer :: OBC                    => NULL()
  type(sponge_CS),               pointer :: sponge_CSp             => NULL()
  type(ALE_sponge_CS),           pointer :: ALE_sponge_CSp         => NULL()
  type(ALE_CS),                  pointer :: ALE_CSp                => NULL()
  type(offline_transport_CS),    pointer :: offline_CSp            => NULL()

  ! These are used for group halo updates.
  type(group_pass_type) :: pass_tau_ustar_psurf
  type(group_pass_type) :: pass_ray
  type(group_pass_type) :: pass_bbl_thick_kv_bbl
  type(group_pass_type) :: pass_T_S_h
  type(group_pass_type) :: pass_T_S
  type(group_pass_type) :: pass_kv_turb
  type(group_pass_type) :: pass_uv_T_S_h
  type(group_pass_type) :: pass_ssh

end type MOM_control_struct

public initialize_MOM
public finish_MOM_initialization
public step_MOM
public step_offline
public MOM_end
public allocate_surface_state
public calculate_surface_state

integer :: id_clock_ocean
integer :: id_clock_dynamics
integer :: id_clock_thermo
integer :: id_clock_tracer
integer :: id_clock_diabatic
integer :: id_clock_continuity  ! also in dynamics s/r
integer :: id_clock_thick_diff
integer :: id_clock_BBL_visc
integer :: id_clock_ml_restrat
integer :: id_clock_diagnostics
integer :: id_clock_Z_diag
integer :: id_clock_init
integer :: id_clock_MOM_init
integer :: id_clock_pass       ! also in dynamics d/r
integer :: id_clock_pass_init  ! also in dynamics d/r
integer :: id_clock_ALE
integer :: id_clock_other
integer :: id_clock_offline_tracer

contains


!> This subroutine orchestrates the time stepping of MOM.  The adiabatic
!! dynamics are stepped by calls to one of the step_MOM_dyn_...routines.
!! The action of lateral processes on tracers occur in calls to
!! advect_tracer and tracer_hordiff.  Vertical mixing and possibly remapping
!! occur inside of diabatic.
subroutine step_MOM(forces, fluxes, sfc_state, Time_start, time_interval, CS)
  type(mech_forcing), intent(in)     :: forces        !< A structure with the driving mechanical forces
  type(forcing),    intent(inout)    :: fluxes        !< pointers to forcing fields
  type(surface),    intent(inout)    :: sfc_state     !< surface ocean state
  type(time_type),  intent(in)       :: Time_start    !< starting time of a segment, as a time type
  real,             intent(in)       :: time_interval !< time interval covered by this run segment, in s.
  type(MOM_control_struct), pointer  :: CS            !< control structure from initialize_MOM

  ! local
  type(ocean_grid_type), pointer :: G ! pointer to a structure containing
                                      ! metrics and related information
  type(verticalGrid_type),  pointer :: GV => NULL()
  integer, save :: nt_debug = 1 ! running number of iterations, for debugging only.
  integer       :: ntstep ! time steps between tracer updates or diabatic forcing
  integer       :: n_max  ! number of steps to take in this call

  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, n
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  real :: dt              ! baroclinic time step (sec)
  real :: dtth            ! time step for thickness diffusion (sec)
  real :: dtdia           ! time step for diabatic processes (sec)
  real :: dt_therm        ! a limited and quantized version of CS%dt_therm (sec)
  real :: dtbt_reset_time ! value of CS%rel_time when DTBT was last calculated (sec)

  real :: mass_src_time   ! The amount of time for the surface mass source from
                          ! precipitation-evaporation, rivers, etc., that should
                          ! be applied to the start of the barotropic solver to
                          ! avoid generating tsunamis, in s.  This is negative
                          ! if the precipation has already been applied to the
                          ! layers, and positive if it will be applied later.

  real :: wt_end, wt_beg
  real :: bbl_time_int    ! The amount of time over which the calculated BBL
                          ! properties will apply, for use in diagnostics.

  logical :: calc_dtbt                 ! Indicates whether the dynamically adjusted
                                       ! barotropic time step needs to be updated.
  logical :: do_advection              ! If true, it is time to advect tracers.
  logical :: do_calc_bbl               ! If true, calculate the boundary layer properties.
  logical :: thermo_does_span_coupling ! If true, thermodynamic forcing spans
                                       ! multiple dynamic timesteps.
  real, dimension(SZI_(CS%G),SZJ_(CS%G)) :: &
    eta_av, &   ! average sea surface height or column mass over a timestep (meter or kg/m2)
    ssh         ! sea surface height based on eta_av (meter or kg/m2)

  real, pointer, dimension(:,:,:) :: &
    u, & ! u : zonal velocity component (m/s)
    v, & ! v : meridional velocity component (m/s)
    h    ! h : layer thickness (meter (Bouss) or kg/m2 (non-Bouss))

  ! Store the layer thicknesses before any changes by the dynamics. This is necessary for remapped
  ! mass transport diagnostics
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)) :: h_pre_dyn
  real :: tot_wt_ssh, Itot_wt_ssh

  type(time_type) :: Time_local
  logical :: showCallTree
  ! These are used for group halo passes.
  logical :: do_pass_kv_turb, do_pass_Ray, do_pass_kv_bbl_thick

  G => CS%G ; GV => CS%GV
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec ; nz = G%ke
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  u => CS%u ; v => CS%v ; h => CS%h

  call cpu_clock_begin(id_clock_ocean)
  call cpu_clock_begin(id_clock_other)

  if (CS%debug) then
    call MOM_state_chksum("Beginning of step_MOM ", u, v, h, CS%uh, CS%vh, G, GV)
    call hchksum(CS%h,"CS%h beginning of step_MOM",G%HI, scale=GV%H_to_m)
  endif

  showCallTree = callTree_showQuery()
  if (showCallTree) call callTree_enter("step_MOM(), MOM.F90")

  ! First determine the time step that is consistent with this call.
  ! It is anticipated that the time step will almost always coincide
  ! with dt. In addition, ntstep is determined, subject to the constraint
  ! that ntstep cannot exceed n_max.
  if (time_interval <= CS%dt) then
    n_max = 1
  else
    n_max = ceiling(time_interval/CS%dt - 0.001)
  endif

  dt    = time_interval / real(n_max)
  dtdia = 0.0
  thermo_does_span_coupling = (CS%thermo_spans_coupling .and. &
                              (CS%dt_therm > 1.5*time_interval))
  if (thermo_does_span_coupling) then
    ! Set dt_therm to be an integer multiple of the coupling time step.
    dt_therm = time_interval * floor(CS%dt_therm / time_interval + 0.001)
    ntstep = floor(dt_therm/dt + 0.001)
  else
    ntstep = MAX(1,MIN(n_max,floor(CS%dt_therm/dt + 0.001)))
    dt_therm = dt*ntstep
  endif

  if (.not.ASSOCIATED(forces%p_surf)) CS%interp_p_surf = .false.

  !---------- Begin setup for group halo pass

  call cpu_clock_begin(id_clock_pass)
  call create_group_pass(CS%pass_tau_ustar_psurf, forces%taux, forces%tauy, G%Domain)
  if (ASSOCIATED(forces%ustar)) &
    call create_group_pass(CS%pass_tau_ustar_psurf, forces%ustar, G%Domain)
  if (ASSOCIATED(forces%p_surf)) &
    call create_group_pass(CS%pass_tau_ustar_psurf, forces%p_surf, G%Domain)

  do_pass_Ray = .FALSE. ; do_pass_kv_bbl_thick = .FALSE.
  if (.not.G%Domain%symmetric) then
    if (associated(CS%visc%Ray_u) .and. associated(CS%visc%Ray_v)) then
      call create_group_pass(CS%pass_ray, CS%visc%Ray_u, CS%visc%Ray_v, G%Domain, &
                             To_North+To_East+SCALAR_PAIR+Omit_corners, CGRID_NE, halo=1)
      do_pass_Ray = .TRUE.
    endif
    if (associated(CS%visc%bbl_thick_u) .and. associated(CS%visc%bbl_thick_v)) then
      call create_group_pass(CS%pass_bbl_thick_kv_bbl, CS%visc%bbl_thick_u, &
                             CS%visc%bbl_thick_v, G%Domain, &
                             To_North+To_East+SCALAR_PAIR+Omit_corners, CGRID_NE, halo=1)
      do_pass_kv_bbl_thick = .TRUE.
    endif
    if (associated(CS%visc%kv_bbl_u) .and. associated(CS%visc%kv_bbl_v)) then
      call create_group_pass(CS%pass_bbl_thick_kv_bbl, CS%visc%kv_bbl_u, &
                             CS%visc%kv_bbl_v, G%Domain, &
                             To_North+To_East+SCALAR_PAIR+Omit_corners, CGRID_NE, halo=1)
      do_pass_kv_bbl_thick = .TRUE.
    endif
  endif
  do_pass_kv_turb = associated(CS%visc%Kv_turb)
  if (associated(CS%visc%Kv_turb)) &
    call create_group_pass(CS%pass_kv_turb, CS%visc%Kv_turb, G%Domain, To_All+Omit_Corners, halo=1)

  if (.not.CS%adiabatic .AND. CS%use_ALE_algorithm ) then
    if (CS%use_temperature) then
      call create_group_pass(CS%pass_T_S_h, CS%tv%T, G%Domain, To_All+Omit_Corners, halo=1)
      call create_group_pass(CS%pass_T_S_h, CS%tv%S, G%Domain, To_All+Omit_Corners, halo=1)
    endif
    call create_group_pass(CS%pass_T_S_h, h, G%Domain, To_All+Omit_Corners, halo=1)
  endif

  if ((CS%adiabatic .OR. CS%diabatic_first) .AND. CS%use_temperature) then
    call create_group_pass(CS%pass_T_S, CS%tv%T, G%Domain, To_All+Omit_Corners, halo=1)
    call create_group_pass(CS%pass_T_S, CS%tv%S, G%Domain, To_All+Omit_Corners, halo=1)
  endif

  !---------- End setup for group halo pass

  if (G%nonblocking_updates) then
    call start_group_pass(CS%pass_tau_ustar_psurf, G%Domain)
  else
    call do_group_pass(CS%pass_tau_ustar_psurf, G%Domain)
  endif
  call cpu_clock_end(id_clock_pass)

  if (ASSOCIATED(CS%tv%frazil))        CS%tv%frazil(:,:)        = 0.0
  if (ASSOCIATED(CS%tv%salt_deficit))  CS%tv%salt_deficit(:,:)  = 0.0
  if (ASSOCIATED(CS%tv%TempxPmE))      CS%tv%TempxPmE(:,:)      = 0.0
  if (ASSOCIATED(CS%tv%internal_heat)) CS%tv%internal_heat(:,:) = 0.0

  CS%rel_time = 0.0

  tot_wt_ssh = 0.0
  do j=js,je ; do i=is,ie ; CS%ave_ssh(i,j) = 0.0 ; ssh(i,j) = CS%missing; enddo ; enddo

  if (associated(CS%VarMix)) then
    call enable_averaging(time_interval, Time_start+set_time(int(time_interval)), &
                          CS%diag)
    call calc_resoln_function(h, CS%tv, G, GV, CS%VarMix)
    call disable_averaging(CS%diag)
  endif

  if (G%nonblocking_updates) &
    call complete_group_pass(CS%pass_tau_ustar_psurf, G%Domain, clock=id_clock_pass)

  if (CS%interp_p_surf) then
    if (.not.ASSOCIATED(CS%p_surf_end))   allocate(CS%p_surf_end(isd:ied,jsd:jed))
    if (.not.ASSOCIATED(CS%p_surf_begin)) allocate(CS%p_surf_begin(isd:ied,jsd:jed))
    if (.not.CS%p_surf_prev_set) then
      do j=jsd,jed ; do i=isd,ied
        CS%p_surf_prev(i,j) = forces%p_surf(i,j)
      enddo ; enddo
      CS%p_surf_prev_set = .true.
    endif
  else
    CS%p_surf_end  => forces%p_surf
  endif

  if (CS%debug) then
    call MOM_state_chksum("Before steps ", u, v, h, CS%uh, CS%vh, G, GV)
    call MOM_forcing_chksum("Before steps", fluxes, G, haloshift=0)
    call check_redundant("Before steps ", u, v, G)
    call check_redundant("Before steps ", forces%taux, forces%tauy, G)
  endif
  call cpu_clock_end(id_clock_other)

  do n=1,n_max

    nt_debug = nt_debug + 1

    ! Set the universally visible time to the middle of the time step
    CS%Time = Time_start + set_time(int(floor(CS%rel_time+0.5*dt+0.5)))
    CS%rel_time = CS%rel_time + dt

    ! Set the local time to the end of the time step.
    Time_local = Time_start + set_time(int(floor(CS%rel_time+0.5)))
    if (showCallTree) call callTree_enter("DT cycles (step_MOM) n=",n)

    !===========================================================================
    ! This is the first place where the diabatic processes and remapping could occur.
    if (CS%diabatic_first .and. (CS%t_dyn_rel_adv==0.0)) then ! do thermodynamics.

      if (thermo_does_span_coupling) then
        dtdia = dt_therm
        if ((fluxes%dt_buoy_accum > 0.0) .and. (dtdia > time_interval) .and. &
            (abs(fluxes%dt_buoy_accum - dtdia) > 1e-6*dtdia)) then
          call MOM_error(FATAL, "step_MOM: Mismatch between long thermodynamic "//&
            "timestep and time over which buoyancy fluxes have been accumulated.")
        endif
        call MOM_error(FATAL, "MOM is not yet set up to have restarts that work "//&
          "with THERMO_SPANS_COUPLING and DIABATIC_FIRST.")
      else
        dtdia = dt*min(ntstep,n_max-(n-1))
      endif

      ! The end-time of the diagnostic interval needs to be set ahead if there
      ! are multiple dynamic time steps worth of thermodynamics applied here.
      call enable_averaging(dtdia, Time_local + &
                                   set_time(int(floor(dtdia-dt+0.5))), CS%diag)

      if (CS%debug) then
        call uvchksum("Pre set_viscous_BBL [uv]", u, v, G%HI, haloshift=1)
        call hchksum(h,"Pre set_viscous_BBL h", G%HI, haloshift=1, scale=GV%H_to_m)
        if (associated(CS%tv%T)) call hchksum(CS%tv%T, "Pre set_viscous_BBL T", G%HI, haloshift=1)
        if (associated(CS%tv%S)) call hchksum(CS%tv%S, "Pre set_viscous_BBL S", G%HI, haloshift=1)
      endif

      !   Calculate the BBL properties and store them inside visc (u,h).
      ! This is here so that CS%visc is updated before diabatic() when
      ! DIABATIC_FIRST=True. Otherwise diabatic() is called after the dynamics
      ! and set_viscous_BBL is called as a part of the dynamic stepping.
      call cpu_clock_begin(id_clock_BBL_visc)
      call set_viscous_BBL(u, v, h, CS%tv, CS%visc, G, GV, CS%set_visc_CSp)
      call cpu_clock_end(id_clock_BBL_visc)

      if (do_pass_Ray) call do_group_pass(CS%pass_ray, G%Domain, clock=id_clock_pass)
      if (do_pass_kv_bbl_thick) &
        call do_group_pass(CS%pass_bbl_thick_kv_bbl, G%Domain, clock=id_clock_pass)
      if (showCallTree) call callTree_wayPoint("done with set_viscous_BBL (diabatic_first)")

      call cpu_clock_begin(id_clock_thermo)

      ! Apply diabatic forcing, do mixing, and regrid.
      call step_MOM_thermo(CS, G, GV, u, v, h, CS%tv, fluxes, dtdia)
      do_pass_kv_turb = associated(CS%visc%Kv_turb)

      ! The diabatic processes are now ahead of the dynamics by dtdia.
      CS%t_dyn_rel_thermo = -dtdia
      if (showCallTree) call callTree_waypoint("finished diabatic_first (step_MOM)")

      call disable_averaging(CS%diag)
      call cpu_clock_end(id_clock_thermo)

    endif ! end of block "(CS%diabatic_first .and. (CS%t_dyn_rel_adv==0.0))"

    !===========================================================================
    ! This is the start of the dynamics stepping part of the algorithm.

    call cpu_clock_begin(id_clock_dynamics)
    call disable_averaging(CS%diag)

    if ((CS%t_dyn_rel_adv == 0.0) .and. CS%thickness_diffuse .and. CS%thickness_diffuse_first) then
      if (thermo_does_span_coupling) then
        dtth = dt_therm
      else
        dtth = dt*min(ntstep,n_max-n+1)
      endif

      call enable_averaging(dtth,Time_local+set_time(int(floor(dtth-dt+0.5))), CS%diag)
      call cpu_clock_begin(id_clock_thick_diff)
      if (associated(CS%VarMix)) &
        call calc_slope_functions(h, CS%tv, dt, G, GV, CS%VarMix)
      call thickness_diffuse(h, CS%uhtr, CS%vhtr, CS%tv, dtth, G, GV, &
                             CS%MEKE, CS%VarMix, CS%CDp, CS%thickness_diffuse_CSp)
      call cpu_clock_end(id_clock_thick_diff)
      call pass_var(h, G%Domain, clock=id_clock_pass) !###, halo=max(2,cont_stensil))
      call disable_averaging(CS%diag)
      if (showCallTree) call callTree_waypoint("finished thickness_diffuse_first (step_MOM)")

      ! Whenever thickness changes let the diag manager know, target grids
      ! for vertical remapping may need to be regenerated.
      call diag_update_remap_grids(CS%diag)
    endif

    ! The bottom boundary layer properties are out-of-date and need to be
    ! recalculated.  This always occurs at the start of a coupling time
    ! step because the externally prescribed stresses may have changed.
    do_calc_bbl = ((CS%t_dyn_rel_adv == 0.0) .or. (n==1))
    if (do_calc_bbl) then
      ! Calculate the BBL properties and store them inside visc (u,h).
      call cpu_clock_begin(id_clock_BBL_visc)
      bbl_time_int = max(dt, min(dt_therm - CS%t_dyn_rel_adv, dt*(1+n_max-n)) )
      call enable_averaging(bbl_time_int, &
                Time_local+set_time(int(bbl_time_int-dt+0.5)), CS%diag)
      call set_viscous_BBL(u, v, h, CS%tv, CS%visc, G, GV, CS%set_visc_CSp)
      call disable_averaging(CS%diag)
      call cpu_clock_end(id_clock_BBL_visc)
      if (showCallTree) call callTree_wayPoint("done with set_viscous_BBL (step_MOM)")
    endif

    if (do_pass_kv_turb) &
      call do_group_pass(CS%pass_kv_turb, G%Domain, clock=id_clock_pass)
    do_pass_kv_turb = .false.

    if (do_calc_bbl) then
      if (G%nonblocking_updates) then
        if (do_pass_Ray) &
          call start_group_pass(CS%pass_Ray, G%Domain, clock=id_clock_pass)
        if (do_pass_kv_bbl_thick) &
          call start_group_pass(CS%pass_bbl_thick_kv_bbl, G%Domain, clock=id_clock_pass)
        ! do_calc_bbl will be set to .false. when the message passing is complete.
      else
        if (do_pass_Ray) &
          call do_group_pass(CS%pass_Ray, G%Domain, clock=id_clock_pass)
        if (do_pass_kv_bbl_thick) &
          call do_group_pass(CS%pass_bbl_thick_kv_bbl, G%Domain, clock=id_clock_pass)
      endif
    endif

    if (CS%interp_p_surf) then
      wt_end = real(n) / real(n_max)
      wt_beg = real(n-1) / real(n_max)
      do j=jsd,jed ; do i=isd,ied
        CS%p_surf_end(i,j) = wt_end * forces%p_surf(i,j) + &
                        (1.0-wt_end) * CS%p_surf_prev(i,j)
        CS%p_surf_begin(i,j) = wt_beg * forces%p_surf(i,j) + &
                        (1.0-wt_beg) * CS%p_surf_prev(i,j)
      enddo ; enddo
    endif

    ! The original velocities might be stored for debugging.
    if (associated(CS%u_prev) .and. associated(CS%v_prev)) then
      do k=1,nz ; do j=jsd,jed ; do I=IsdB,IedB
        CS%u_prev(I,j,k) = u(I,j,k)
      enddo ; enddo ; enddo
      do k=1,nz ; do J=JsdB,JedB ; do i=isd,ied
        CS%v_prev(I,j,k) = u(I,j,k)
      enddo ; enddo ; enddo
    endif

    ! Store pre-dynamics layer thicknesses so that mass fluxes are remapped correctly
    do k=1,nz ; do j=jsd,jed ; do i=isd,ied
      h_pre_dyn(i,j,k) = h(i,j,k)
    enddo ; enddo ; enddo

    if (G%nonblocking_updates) then ; if (do_calc_bbl) then
      if (do_pass_Ray) &
        call complete_group_pass(CS%pass_Ray, G%Domain, clock=id_clock_pass)
      if (do_pass_kv_bbl_thick) &
        call complete_group_pass(CS%pass_bbl_thick_kv_bbl, G%Domain, clock=id_clock_pass)
    endif ; endif

    if (CS%do_dynamics .and. CS%split) then !--------------------------- start SPLIT
      ! This section uses a split time stepping scheme for the dynamic equations,
      ! basically the stacked shallow water equations with viscosity.

      calc_dtbt = .false.
      if ((CS%dtbt_reset_period >= 0.0) .and. &
          ((n==1) .or. (CS%dtbt_reset_period == 0.0) .or. &
           (CS%rel_time >= dtbt_reset_time + 0.999*CS%dtbt_reset_period))) then
        calc_dtbt = .true.
        dtbt_reset_time = CS%rel_time
      endif

      mass_src_time = CS%t_dyn_rel_thermo
      call step_MOM_dyn_split_RK2(u, v, h, CS%tv, CS%visc, &
                  Time_local, dt, forces, CS%p_surf_begin, CS%p_surf_end, &
                  mass_src_time, dt_therm, CS%uh, CS%vh, CS%uhtr, CS%vhtr, &
                  eta_av, G, GV, CS%dyn_split_RK2_CSp, calc_dtbt, CS%VarMix, CS%MEKE)
      if (showCallTree) call callTree_waypoint("finished step_MOM_dyn_split (step_MOM)")

    elseif (CS%do_dynamics) then ! --------------------------------------------------- not SPLIT
      !   This section uses an unsplit stepping scheme for the dynamic
      ! equations; basically the stacked shallow water equations with viscosity.
      ! Because the time step is limited by CFL restrictions on the external
      ! gravity waves, the unsplit is usually much less efficient that the split
      ! approaches. But because of its simplicity, the unsplit method is very
      ! useful for debugging purposes.

      if (CS%use_RK2) then
        call step_MOM_dyn_unsplit_RK2(u, v, h, CS%tv, CS%visc, Time_local, dt, forces, &
                 CS%p_surf_begin, CS%p_surf_end, CS%uh, CS%vh, CS%uhtr, CS%vhtr, &
                 eta_av, G, GV, CS%dyn_unsplit_RK2_CSp, CS%VarMix, CS%MEKE)
      else
        call step_MOM_dyn_unsplit(u, v, h, CS%tv, CS%visc, Time_local, dt, forces, &
                 CS%p_surf_begin, CS%p_surf_end, CS%uh, CS%vh, CS%uhtr, CS%vhtr, &
                 eta_av, G, GV, CS%dyn_unsplit_CSp, CS%VarMix, CS%MEKE)
      endif
      if (showCallTree) call callTree_waypoint("finished step_MOM_dyn_unsplit (step_MOM)")

    endif ! -------------------------------------------------- end SPLIT


    if (CS%thickness_diffuse .and. .not.CS%thickness_diffuse_first) then
      call cpu_clock_begin(id_clock_thick_diff)

      if (CS%debug) call hchksum(h,"Pre-thickness_diffuse h", G%HI, haloshift=0, scale=GV%H_to_m)

      if (associated(CS%VarMix)) &
        call calc_slope_functions(h, CS%tv, dt, G, GV, CS%VarMix)
      call thickness_diffuse(h, CS%uhtr, CS%vhtr, CS%tv, dt, G, GV, &
                             CS%MEKE, CS%VarMix, CS%CDp, CS%thickness_diffuse_CSp)

      if (CS%debug) call hchksum(h,"Post-thickness_diffuse h", G%HI, haloshift=1, scale=GV%H_to_m)
      call cpu_clock_end(id_clock_thick_diff)
      call pass_var(h, G%Domain, clock=id_clock_pass) !###, halo=max(2,cont_stensil))
      if (showCallTree) call callTree_waypoint("finished thickness_diffuse (step_MOM)")
    endif

    ! apply the submesoscale mixed layer restratification parameterization
    if (CS%mixedlayer_restrat) then
      if (CS%debug) then
        call hchksum(h,"Pre-mixedlayer_restrat h", G%HI, haloshift=1, scale=GV%H_to_m)
        call uvchksum("Pre-mixedlayer_restrat uhtr", &
                      CS%uhtr, CS%vhtr, G%HI, haloshift=0)
      endif
      call cpu_clock_begin(id_clock_ml_restrat)
      call mixedlayer_restrat(h, CS%uhtr, CS%vhtr, CS%tv, forces, dt, CS%visc%MLD, &
                              CS%VarMix, G, GV, CS%mixedlayer_restrat_CSp)
      call cpu_clock_end(id_clock_ml_restrat)
      call pass_var(h, G%Domain, clock=id_clock_pass) !###, halo=max(2,cont_stensil))
      if (CS%debug) then
        call hchksum(h,"Post-mixedlayer_restrat h", G%HI, haloshift=1, scale=GV%H_to_m)
        call uvchksum("Post-mixedlayer_restrat [uv]htr", &
                      CS%uhtr, CS%vhtr, G%HI, haloshift=0)
      endif
    endif

    ! Whenever thickness changes let the diag manager know, target grids
    ! for vertical remapping may need to be regenerated.
    call diag_update_remap_grids(CS%diag)

    if (CS%useMEKE) call step_forward_MEKE(CS%MEKE, h, CS%VarMix%SN_u, CS%VarMix%SN_v, &
                                           CS%visc, dt, G, GV, CS%MEKE_CSp, CS%uhtr, CS%vhtr)
    call disable_averaging(CS%diag)

    ! Advance the dynamics time by dt.
    CS%t_dyn_rel_adv = CS%t_dyn_rel_adv + dt
    CS%t_dyn_rel_thermo = CS%t_dyn_rel_thermo + dt
    CS%t_dyn_rel_diag = CS%t_dyn_rel_diag + dt

    call cpu_clock_end(id_clock_dynamics)

    !===========================================================================
    ! This is the start of the tracer advection part of the algorithm.

    if (thermo_does_span_coupling) then
      do_advection = (CS%t_dyn_rel_adv + 0.5*dt > dt_therm)
    else
      do_advection = ((MOD(n,ntstep) == 0) .or. (n==n_max))
    endif

    if (do_advection) then ! Do advective transport and lateral tracer mixing.

      if (CS%debug) then
        call cpu_clock_begin(id_clock_other)
        call uvchksum("Pre-advection [uv]", u, v, G%HI, haloshift=2)
        call hchksum(h,"Pre-advection h", G%HI, haloshift=1, scale=GV%H_to_m)
        call uvchksum("Pre-advection uhtr", CS%uhtr, CS%vhtr, G%HI, &
                      haloshift=0, scale=GV%H_to_m)
      ! call MOM_state_chksum("Pre-advection ", u, v, &
      !                       h, CS%uhtr, CS%vhtr, G, GV, haloshift=1)
          if (associated(CS%tv%T)) call hchksum(CS%tv%T, "Pre-advection T", G%HI, haloshift=1)
          if (associated(CS%tv%S)) call hchksum(CS%tv%S, "Pre-advection S", G%HI, haloshift=1)
          if (associated(CS%tv%frazil)) call hchksum(CS%tv%frazil, &
                         "Pre-advection frazil", G%HI, haloshift=0)
          if (associated(CS%tv%salt_deficit)) call hchksum(CS%tv%salt_deficit, &
                         "Pre-advection salt deficit", G%HI, haloshift=0)
      ! call MOM_thermo_chksum("Pre-advection ", CS%tv, G)
        call check_redundant("Pre-advection ", u, v, G)
        call cpu_clock_end(id_clock_other)
      endif

      call cpu_clock_begin(id_clock_thermo) ; call cpu_clock_begin(id_clock_tracer)
      call enable_averaging(CS%t_dyn_rel_adv, Time_local, CS%diag)

      call advect_tracer(h, CS%uhtr, CS%vhtr, CS%OBC, CS%t_dyn_rel_adv, G, GV, &
                         CS%tracer_adv_CSp, CS%tracer_Reg)
      call tracer_hordiff(h, CS%t_dyn_rel_adv, CS%MEKE, CS%VarMix, G, GV, &
                          CS%tracer_diff_CSp, CS%tracer_Reg, CS%tv)
      if (showCallTree) call callTree_waypoint("finished tracer advection/diffusion (step_MOM)")
      call cpu_clock_end(id_clock_tracer) ; call cpu_clock_end(id_clock_thermo)

      call cpu_clock_begin(id_clock_other) ; call cpu_clock_begin(id_clock_diagnostics)
      call post_transport_diagnostics(G, GV, CS, CS%diag, CS%t_dyn_rel_adv, h, h_pre_dyn)
      ! Rebuild the remap grids now that we've posted the fields which rely on thicknesses
      ! from before the dynamics calls
      call diag_update_remap_grids(CS%diag)

      call disable_averaging(CS%diag)
      call cpu_clock_end(id_clock_diagnostics) ; call cpu_clock_end(id_clock_other)

      ! Reset the accumulated transports to 0 and record that the dynamics
      ! and advective times now agree.
      call cpu_clock_begin(id_clock_thermo) ; call cpu_clock_begin(id_clock_tracer)
      CS%uhtr(:,:,:) = 0.0
      CS%vhtr(:,:,:) = 0.0
      CS%t_dyn_rel_adv = 0.0
      call cpu_clock_end(id_clock_tracer) ; call cpu_clock_end(id_clock_thermo)

      if (CS%diabatic_first .and. CS%use_temperature) then
        ! Temperature and salinity need halo updates because they will be used
        ! in the dynamics before they are changed again.
        call do_group_pass(CS%pass_T_S, G%Domain, clock=id_clock_pass)
      endif

    endif

    !===========================================================================
    ! This is the second place where the diabatic processes and remapping could occur.
    if (CS%t_dyn_rel_adv == 0.0) then
      call cpu_clock_begin(id_clock_thermo)

      if (.not.CS%diabatic_first) then
        dtdia = CS%t_dyn_rel_thermo
        if (thermo_does_span_coupling .and. (abs(dt_therm - dtdia) > 1e-6*dt_therm)) then
          call MOM_error(FATAL, "step_MOM: Mismatch between dt_therm and dtdia "//&
                         "before call to diabatic.")
        endif

        call enable_averaging(CS%t_dyn_rel_thermo, Time_local, CS%diag)

        ! Apply diabatic forcing, do mixing, and regrid.
        call step_MOM_thermo(CS, G, GV, u, v, h, CS%tv, fluxes, dtdia)
        do_pass_kv_turb = associated(CS%visc%Kv_turb)

        call disable_averaging(CS%diag)

      endif

      if (CS%diabatic_first .and. abs(CS%t_dyn_rel_thermo) > 1e-6*dt) call MOM_error(FATAL, &
              "step_MOM: Mismatch between the dynamics and diabatic times "//&
              "with DIABATIC_FIRST.")
      ! Record that the dynamics and diabatic processes are synchronized.
      CS%t_dyn_rel_thermo = 0.0
      call cpu_clock_end(id_clock_thermo)
    endif

    call cpu_clock_begin(id_clock_dynamics)

    ! Determining the time-average sea surface height is part of the algorithm.
    ! This may be eta_av if Boussinesq, or need to be diagnosed if not.
    tot_wt_ssh = tot_wt_ssh + dt
    call find_eta(h, CS%tv, GV%g_Earth, G, GV, ssh, eta_av)
    do j=js,je ; do i=is,ie
      CS%ave_ssh(i,j) = CS%ave_ssh(i,j) + dt*ssh(i,j)
    enddo ; enddo
    call cpu_clock_end(id_clock_dynamics)

    !===========================================================================
    ! Calculate diagnostics at the end of the time step.
    call cpu_clock_begin(id_clock_other) ; call cpu_clock_begin(id_clock_diagnostics)

    call enable_averaging(dt, Time_local, CS%diag)
    ! These diagnostics are available every time step.
    if (CS%id_u > 0) call post_data(CS%id_u, u, CS%diag)
    if (CS%id_v > 0) call post_data(CS%id_v, v, CS%diag)
    if (CS%id_h > 0) call post_data(CS%id_h, h, CS%diag)
    if (CS%id_ssh_inst > 0) call post_data(CS%id_ssh_inst, ssh, CS%diag)
    call disable_averaging(CS%diag)

    if (CS%t_dyn_rel_adv == 0.0) then
      ! Diagnostics that require the complete state to be up-to-date can be calculated.

      call enable_averaging(CS%t_dyn_rel_diag, Time_local, CS%diag)
      call calculate_diagnostic_fields(u, v, h, CS%uh, CS%vh, CS%tv, CS%ADp, &
                          CS%CDp, fluxes, CS%t_dyn_rel_diag, G, GV, CS%diagnostics_CSp)
      call post_TS_diagnostics(CS, G, GV, CS%tv, CS%diag, CS%t_dyn_rel_diag)
      if (showCallTree) call callTree_waypoint("finished calculate_diagnostic_fields (step_MOM)")
      call disable_averaging(CS%diag)
      CS%t_dyn_rel_diag = 0.0

      call cpu_clock_begin(id_clock_Z_diag)
      if (Time_local + set_time(int(0.5*dt_therm)) > CS%Z_diag_time) then
        call enable_averaging(real(time_type_to_real(CS%Z_diag_interval)), &
                              CS%Z_diag_time, CS%diag)
        call calculate_Z_diag_fields(u, v, h, ssh, fluxes%frac_shelf_h, &
                                     G, GV, CS%diag_to_Z_CSp)
        CS%Z_diag_time = CS%Z_diag_time + CS%Z_diag_interval
        call disable_averaging(CS%diag)
        if (showCallTree) call callTree_waypoint("finished calculate_Z_diag_fields (step_MOM)")
      endif
      call cpu_clock_end(id_clock_Z_diag)
    endif
    call cpu_clock_end(id_clock_diagnostics) ; call cpu_clock_end(id_clock_other)

    if (showCallTree) call callTree_leave("DT cycles (step_MOM)")

  enddo ! complete the n loop

  call cpu_clock_begin(id_clock_other)

  Itot_wt_ssh = 1.0/tot_wt_ssh
  do j=js,je ; do i=is,ie
    CS%ave_ssh(i,j) = CS%ave_ssh(i,j)*Itot_wt_ssh
    ssh(i,j) = CS%ave_ssh(i,j)
  enddo ; enddo
  call adjust_ssh_for_p_atm(CS, G, GV, CS%ave_ssh, forces%p_surf_SSH)

  if (CS%interp_p_surf) then ; do j=jsd,jed ; do i=isd,ied
    CS%p_surf_prev(i,j) = forces%p_surf(i,j)
  enddo ; enddo ; endif

  if (showCallTree) call callTree_waypoint("calling calculate_surface_state (step_MOM)")
  call calculate_surface_state(sfc_state, u, v, h, CS%ave_ssh, G, GV, CS)

  ! Do diagnostics that only occur at the end of a complete forcing step.
  call cpu_clock_begin(id_clock_diagnostics)
  call enable_averaging(dt*n_max, Time_local, CS%diag)
  call post_integrated_diagnostics(CS, G, GV, CS%diag, dt*n_max, CS%tv, ssh, fluxes)
  call post_surface_diagnostics(CS, G, CS%diag, sfc_state)
  call disable_averaging(CS%diag)
  call cpu_clock_end(id_clock_diagnostics)

  call cpu_clock_end(id_clock_other)

  if (showCallTree) call callTree_leave("step_MOM()")
  call cpu_clock_end(id_clock_ocean)

end subroutine step_MOM

!> MOM_step_thermo orchestrates the thermodynamic time stepping and vertical
!! remapping, via calls to diabatic (or adiabatic) and ALE_main.
subroutine step_MOM_thermo(CS, G, GV, u, v, h, tv, fluxes, dtdia)
  type(MOM_control_struct), intent(inout) :: CS     !< control structure
  type(ocean_grid_type),    intent(inout) :: G      !< ocean grid structure
  type(verticalGrid_type),  intent(inout) :: GV     !< ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                            intent(inout) :: u      !< zonal velocity (m/s)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                            intent(inout) :: v      !< meridional velocity (m/s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  &
                            intent(inout) :: h      !< layer thickness (m or kg/m2)
  type(thermo_var_ptrs),    intent(inout) :: tv     !< A structure pointing to various thermodynamic variables
  type(forcing),            intent(inout) :: fluxes !< pointers to forcing fields
  real,                     intent(in)    :: dtdia  !< The time interval over which to advance, in s

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: eta_predia, eta_preale
  integer :: i, j, k, is, ie, js, je, nz! , Isq, Ieq, Jsq, Jeq, n
  logical :: use_ice_shelf ! Needed for selecting the right ALE interface.
  logical :: showCallTree

  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec ; nz = G%ke
  showCallTree = callTree_showQuery()
  if (showCallTree) call callTree_enter("step_MOM_thermo(), MOM.F90")

  use_ice_shelf = .false.
  if (associated(fluxes%frac_shelf_h)) use_ice_shelf = .true.

  if (.not.CS%adiabatic) then

    if (CS%debug) then
      call uvchksum("Pre-diabatic [uv]", u, v, G%HI, haloshift=2)
      call hchksum(h,"Pre-diabatic h", G%HI, haloshift=1, scale=GV%H_to_m)
      call uvchksum("Pre-diabatic [uv]h", CS%uhtr, CS%vhtr, G%HI, &
                    haloshift=0, scale=GV%H_to_m)
    ! call MOM_state_chksum("Pre-diabatic ",u, v, h, CS%uhtr, CS%vhtr, G, GV)
      call MOM_thermo_chksum("Pre-diabatic ", CS%tv, G,haloshift=0)
      call check_redundant("Pre-diabatic ", u, v, G)
      call MOM_forcing_chksum("Pre-diabatic", fluxes, G, haloshift=0)
    endif

    if (CS%id_u_predia > 0) call post_data(CS%id_u_predia, u, CS%diag)
    if (CS%id_v_predia > 0) call post_data(CS%id_v_predia, v, CS%diag)
    if (CS%id_h_predia > 0) call post_data(CS%id_h_predia, h, CS%diag)
    if (CS%id_T_predia > 0) call post_data(CS%id_T_predia, CS%tv%T, CS%diag)
    if (CS%id_S_predia > 0) call post_data(CS%id_S_predia, CS%tv%S, CS%diag)
    if (CS%id_e_predia > 0) then
      call find_eta(h, CS%tv, GV%g_Earth, G, GV, eta_predia)
      call post_data(CS%id_e_predia, eta_predia, CS%diag)
    endif

    call cpu_clock_begin(id_clock_diabatic)
    call diabatic(u, v, h, tv, CS%Hml, fluxes, CS%visc, CS%ADp, CS%CDp, &
                  dtdia, G, GV, CS%diabatic_CSp)
    fluxes%fluxes_used = .true.
    call cpu_clock_end(id_clock_diabatic)

    if (CS%id_u_preale > 0) call post_data(CS%id_u_preale, u,    CS%diag)
    if (CS%id_v_preale > 0) call post_data(CS%id_v_preale, v,    CS%diag)
    if (CS%id_h_preale > 0) call post_data(CS%id_h_preale, h,    CS%diag)
    if (CS%id_T_preale > 0) call post_data(CS%id_T_preale, tv%T, CS%diag)
    if (CS%id_S_preale > 0) call post_data(CS%id_S_preale, tv%S, CS%diag)
    if (CS%id_e_preale > 0) then
      call find_eta(h, tv, GV%g_Earth, G, GV, eta_preale)
      call post_data(CS%id_e_preale, eta_preale, CS%diag)
    endif

    if (showCallTree) call callTree_waypoint("finished diabatic (step_MOM_thermo)")

    ! Regridding/remapping is done here, at end of thermodynamics time step
    ! (that may comprise several dynamical time steps)
    ! The routine 'ALE_main' can be found in 'MOM_ALE.F90'.
    if ( CS%use_ALE_algorithm ) then
!         call pass_vector(u, v, G%Domain)
      call do_group_pass(CS%pass_T_S_h, G%Domain)

      ! update squared quantities
      if (associated(CS%S_squared)) then ; do k=1,nz ; do j=js,je ; do i=is,ie
        CS%S_squared(i,j,k) = tv%S(i,j,k)**2
      enddo ; enddo ; enddo ; endif
      if (associated(CS%T_squared)) then ; do k=1,nz ; do j=js,je ; do i=is,ie
        CS%T_squared(i,j,k) = tv%T(i,j,k)**2
      enddo ; enddo ; enddo ; endif

      if (CS%debug) then
        call MOM_state_chksum("Pre-ALE ", u, v, h, CS%uh, CS%vh, G, GV)
        call hchksum(tv%T,"Pre-ALE T", G%HI, haloshift=1)
        call hchksum(tv%S,"Pre-ALE S", G%HI, haloshift=1)
        call check_redundant("Pre-ALE ", u, v, G)
      endif
      call cpu_clock_begin(id_clock_ALE)
      if (use_ice_shelf) then
        call ALE_main(G, GV, h, u, v, tv, CS%tracer_Reg, CS%ALE_CSp, dtdia, &
                      fluxes%frac_shelf_h)
      else
        call ALE_main(G, GV, h, u, v, tv, CS%tracer_Reg, CS%ALE_CSp, dtdia)
      endif

      if (showCallTree) call callTree_waypoint("finished ALE_main (step_MOM_thermo)")
      call cpu_clock_end(id_clock_ALE)
    endif   ! endif for the block "if ( CS%use_ALE_algorithm )"

    call do_group_pass(CS%pass_uv_T_S_h, G%Domain, clock=id_clock_pass)

    if (CS%debug .and. CS%use_ALE_algorithm) then
      call MOM_state_chksum("Post-ALE ", u, v, h, CS%uh, CS%vh, G, GV)
      call hchksum(tv%T, "Post-ALE T", G%HI, haloshift=1)
      call hchksum(tv%S, "Post-ALE S", G%HI, haloshift=1)
      call check_redundant("Post-ALE ", u, v, G)
    endif

    ! Whenever thickness changes let the diag manager know, target grids
    ! for vertical remapping may need to be regenerated. This needs to
    ! happen after the H update and before the next post_data.
    call diag_update_remap_grids(CS%diag)

    call post_diags_TS_vardec(G, CS, dtdia)

    if (CS%debug) then
      call uvchksum("Post-diabatic u", u, v, G%HI, haloshift=2)
      call hchksum(h, "Post-diabatic h", G%HI, haloshift=1, scale=GV%H_to_m)
      call uvchksum("Post-diabatic [uv]h", CS%uhtr, CS%vhtr, G%HI, &
                    haloshift=0, scale=GV%H_to_m)
    ! call MOM_state_chksum("Post-diabatic ", u, v, &
    !                       h, CS%uhtr, CS%vhtr, G, GV, haloshift=1)
      if (associated(tv%T)) call hchksum(tv%T, "Post-diabatic T", G%HI, haloshift=1)
      if (associated(tv%S)) call hchksum(tv%S, "Post-diabatic S", G%HI, haloshift=1)
      if (associated(tv%frazil)) call hchksum(tv%frazil, &
                               "Post-diabatic frazil", G%HI, haloshift=0)
      if (associated(tv%salt_deficit)) call hchksum(tv%salt_deficit, &
                               "Post-diabatic salt deficit", G%HI, haloshift=0)
    ! call MOM_thermo_chksum("Post-diabatic ", tv, G)
      call check_redundant("Post-diabatic ", u, v, G)
    endif

  else   ! complement of "if (.not.CS%adiabatic)"

    call cpu_clock_begin(id_clock_diabatic)
    call adiabatic(h, tv, fluxes, dtdia, G, GV, CS%diabatic_CSp)
    fluxes%fluxes_used = .true.
    call cpu_clock_end(id_clock_diabatic)

    if (CS%use_temperature) then
      call do_group_pass(CS%pass_T_S, G%Domain, clock=id_clock_pass)
      if (CS%debug) then
        if (associated(tv%T)) call hchksum(tv%T, "Post-diabatic T", G%HI, haloshift=1)
        if (associated(tv%S)) call hchksum(tv%S, "Post-diabatic S", G%HI, haloshift=1)
      endif
    endif

  endif   ! endif for the block "if (.not.CS%adiabatic)"

  if (showCallTree) call callTree_leave("step_MOM_thermo(), MOM.F90")

end subroutine step_MOM_thermo


!> step_offline is the main driver for running tracers offline in MOM6. This has been primarily
!! developed with ALE configurations in mind. Some work has been done in isopycnal configuration, but
!! the work is very preliminary. Some more detail about this capability along with some of the subroutines
!! called here can be found in tracers/MOM_offline_control.F90
subroutine step_offline(forces, fluxes, sfc_state, Time_start, time_interval, CS)
  type(mech_forcing), intent(in)     :: forces        !< A structure with the driving mechanical forces
  type(forcing),    intent(inout)    :: fluxes        !< pointers to forcing fields
  type(surface),    intent(inout)    :: sfc_state     !< surface ocean state
  type(time_type),  intent(in)       :: Time_start    !< starting time of a segment, as a time type
  real,             intent(in)       :: time_interval !< time interval
  type(MOM_control_struct), pointer  :: CS            !< control structure from initialize_MOM

  ! Local pointers
  type(ocean_grid_type),      pointer :: G  => NULL() ! Pointer to a structure containing
                                                      ! metrics and related information
  type(verticalGrid_type),    pointer :: GV => NULL() ! Pointer to structure containing information
                                                      ! about the vertical grid

  logical :: first_iter    !< True if this is the first time step_offline has been called in a given interval
  logical :: last_iter     !< True if this is the last time step_tracer is to be called in an offline interval
  logical :: do_vertical   !< If enough time has elapsed, do the diabatic tracer sources/sinks
  logical :: adv_converged !< True if all the horizontal fluxes have been used

  integer :: dt_offline, dt_offline_vertical
  logical :: skip_diffusion
  integer :: id_eta_diff_end

  integer, pointer :: accumulated_time
  integer :: i,j,k
  integer :: is, ie, js, je, isd, ied, jsd, jed

  ! 3D pointers
  real, dimension(:,:,:), pointer   :: &
    uhtr, vhtr, &
    eatr, ebtr, &
    h_end

  ! 2D Array for diagnostics
  real, dimension(SZI_(CS%G),SZJ_(CS%G)) :: eta_pre, eta_end
  type(time_type) :: Time_end    ! End time of a segment, as a time type

  ! Grid-related pointer assignments
  G => CS%G
  GV => CS%GV

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec
  isd = G%isd  ; ied = G%ied  ; jsd = G%jsd  ; jed = G%jed

  call cpu_clock_begin(id_clock_offline_tracer)
  call extract_offline_main(CS%offline_CSp, uhtr, vhtr, eatr, ebtr, h_end, accumulated_time, &
                            dt_offline, dt_offline_vertical, skip_diffusion)
  Time_end = increment_date(Time_start, seconds=floor(time_interval+0.001))
  call enable_averaging(time_interval, Time_end, CS%diag)

  ! Check to see if this is the first iteration of the offline interval
  if(accumulated_time==0) then
    first_iter = .true.
  else ! This is probably unnecessary but is used to guard against unwanted behavior
    first_iter = .false.
  endif

  ! Check to see if vertical tracer functions should  be done
  if ( mod(accumulated_time, dt_offline_vertical) == 0 ) then
    do_vertical = .true.
  else
    do_vertical = .false.
  endif

  ! Increment the amount of time elapsed since last read and check if it's time to roll around
  accumulated_time = mod(accumulated_time + int(time_interval), dt_offline)
  if(accumulated_time==0) then
    last_iter = .true.
  else
    last_iter = .false.
  endif

  if(CS%use_ALE_algorithm) then
    ! If this is the first iteration in the offline timestep, then we need to read in fields and
    ! perform the main advection.
    if (first_iter) then
      if(is_root_pe()) print *, "Reading in new offline fields"
      ! Read in new transport and other fields
      ! call update_transport_from_files(G, GV, CS%offline_CSp, h_end, eatr, ebtr, uhtr, vhtr, &
      !     CS%tv%T, CS%tv%S, fluxes, CS%use_ALE_algorithm)
      ! call update_transport_from_arrays(CS%offline_CSp)
      call update_offline_fields(CS%offline_CSp, CS%h, fluxes, CS%use_ALE_algorithm)

      ! Apply any fluxes into the ocean
      call offline_fw_fluxes_into_ocean(G, GV, CS%offline_CSp, fluxes, CS%h)

      if (.not.CS%diabatic_first) then
        call offline_advection_ale(fluxes, Time_start, time_interval, CS%offline_CSp, id_clock_ALE, &
            CS%h, uhtr, vhtr, converged=adv_converged)

        ! Redistribute any remaining transport
        call offline_redistribute_residual(CS%offline_CSp, CS%h, uhtr, vhtr, adv_converged)

        ! Perform offline diffusion if requested
        if (.not. skip_diffusion) then
          if (associated(CS%VarMix)) then
            call pass_var(CS%h,G%Domain)
            call calc_resoln_function(CS%h, CS%tv, G, GV, CS%VarMix)
            call calc_slope_functions(CS%h, CS%tv, REAL(dt_offline), G, GV, CS%VarMix)
          endif
          call tracer_hordiff(CS%h, REAL(dt_offline), CS%MEKE, CS%VarMix, G, GV, &
              CS%tracer_diff_CSp, CS%tracer_Reg, CS%tv)
        endif
      endif
    endif
    ! The functions related to column physics of tracers is performed separately in ALE mode
    if (do_vertical) then
      call offline_diabatic_ale(fluxes, Time_start, Time_end, CS%offline_CSp, CS%h, eatr, ebtr)
    endif

    ! Last thing that needs to be done is the final ALE remapping
    if(last_iter) then
      if (CS%diabatic_first) then
        call offline_advection_ale(fluxes, Time_start, time_interval, CS%offline_CSp, id_clock_ALE, &
            CS%h, uhtr, vhtr, converged=adv_converged)

        ! Redistribute any remaining transport and perform the remaining advection
        call offline_redistribute_residual(CS%offline_CSp, CS%h, uhtr, vhtr, adv_converged)
                ! Perform offline diffusion if requested
        if (.not. skip_diffusion) then
          if (associated(CS%VarMix)) then
            call pass_var(CS%h,G%Domain)
            call calc_resoln_function(CS%h, CS%tv, G, GV, CS%VarMix)
            call calc_slope_functions(CS%h, CS%tv, REAL(dt_offline), G, GV, CS%VarMix)
          endif
          call tracer_hordiff(CS%h, REAL(dt_offline), CS%MEKE, CS%VarMix, G, GV, &
              CS%tracer_diff_CSp, CS%tracer_Reg, CS%tv)
        endif
      endif

      if(is_root_pe()) print *, "Last iteration of offline interval"

      ! Apply freshwater fluxes out of the ocean
      call offline_fw_fluxes_out_ocean(G, GV, CS%offline_CSp, fluxes, CS%h)
      ! These diagnostic can be used to identify which grid points did not converge within
      ! the specified number of advection sub iterations
      call post_offline_convergence_diags(CS%offline_CSp, CS%h, h_end, uhtr, vhtr)

      ! Call ALE one last time to make sure that tracers are remapped onto the layer thicknesses
      ! stored from the forward run
      call cpu_clock_begin(id_clock_ALE)
      call ALE_offline_tracer_final( G, GV, CS%h, CS%tv, h_end, CS%tracer_Reg, CS%ALE_CSp)
      call cpu_clock_end(id_clock_ALE)
      call pass_var(CS%h, G%Domain)
    endif
  else ! NON-ALE MODE...NOT WELL TESTED
    call MOM_error(WARNING, &
        "Offline tracer mode in non-ALE configuration has not been thoroughly tested")
    ! Note that for the layer mode case, the calls to tracer sources and sinks is embedded in
    ! main_offline_advection_layer. Warning: this may not be appropriate for tracers that
    ! exchange with the atmosphere
    if(time_interval .NE. dt_offline) then
      call MOM_error(FATAL, &
          "For offline tracer mode in a non-ALE configuration, dt_offline must equal time_interval")
    endif
    call update_offline_fields(CS%offline_CSp, CS%h, fluxes, CS%use_ALE_algorithm)
    call offline_advection_layer(fluxes, Time_start, time_interval, CS%offline_CSp, &
        CS%h, eatr, ebtr, uhtr, vhtr)
    ! Perform offline diffusion if requested
    if (.not. skip_diffusion) then
      call tracer_hordiff(h_end, REAL(dt_offline), CS%MEKE, CS%VarMix, G, GV, &
        CS%tracer_diff_CSp, CS%tracer_Reg, CS%tv)
    endif

    CS%h = h_end

    call pass_var(CS%tv%T, G%Domain)
    call pass_var(CS%tv%S, G%Domain)
    call pass_var(CS%h, G%Domain)

  endif

  call adjust_ssh_for_p_atm(CS, G, GV, CS%ave_ssh, forces%p_surf_SSH)
  call calculate_surface_state(sfc_state, CS%u, CS%v, CS%h, CS%ave_ssh, G, GV, CS)

  call disable_averaging(CS%diag)
  call pass_var(CS%tv%T,G%Domain)
  call pass_var(CS%tv%S,G%Domain)
  call pass_var(CS%h,G%Domain)

  fluxes%fluxes_used = .true.

  call cpu_clock_end(id_clock_offline_tracer)

end subroutine step_offline

!> This subroutine initializes MOM.
subroutine initialize_MOM(Time, param_file, dirs, CS, Time_in, offline_tracer_mode)
  type(time_type), target,   intent(inout) :: Time        !< model time, set in this routine
  type(param_file_type),     intent(out)   :: param_file  !< structure indicating paramater file to parse
  type(directories),         intent(out)   :: dirs        !< structure with directory paths
  type(MOM_control_struct),  pointer       :: CS          !< pointer set in this routine to MOM control structure
  type(time_type), optional, intent(in)    :: Time_in     !< time passed to MOM_initialize_state when
                                                          !! model is not being started from a restart file
  logical,         optional, intent(out)   :: offline_tracer_mode !< True if tracers are being run offline

  ! local
  type(ocean_grid_type),  pointer :: G => NULL() ! A pointer to a structure with metrics and related
  type(hor_index_type)            :: HI  !  A hor_index_type for array extents
  type(verticalGrid_type), pointer :: GV => NULL()
  type(dyn_horgrid_type), pointer :: dG => NULL()
  type(diag_ctrl),        pointer :: diag

  character(len=4), parameter :: vers_num = 'v2.0'

! This include declares and sets the variable "version".
#include "version_variable.h"

  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB
  real    :: dtbt
  real    :: Z_diag_int  ! minimum interval between calc depth-space diagnostics (sec)

  real, allocatable, dimension(:,:,:) :: e   ! interface heights (meter)
  real, allocatable, dimension(:,:)   :: eta ! free surface height (m) or bottom press (Pa)
  real, allocatable, dimension(:,:)   :: area_shelf_h ! area occupied by ice shelf
  real, dimension(:,:), allocatable, target  :: frac_shelf_h ! fraction of total area occupied by ice shelf
  real, dimension(:,:), pointer :: shelf_area
  type(MOM_restart_CS),  pointer      :: restart_CSp_tmp => NULL()
  type(group_pass_type) :: tmp_pass_uv_T_S_h

  real    :: default_val       ! default value for a parameter
  logical :: write_geom_files  ! If true, write out the grid geometry files.
  logical :: new_sim
  logical :: use_geothermal    ! If true, apply geothermal heating.
  logical :: use_EOS           ! If true, density calculated from T & S using an equation of state.
  logical :: symmetric         ! If true, use symmetric memory allocation.
  logical :: save_IC           ! If true, save the initial conditions.
  logical :: do_unit_tests     ! If true, call unit tests.
  logical :: test_grid_copy = .false.
  logical :: use_ice_shelf     ! Needed for ALE
  logical :: global_indexing   ! If true use global horizontal index values instead
                               ! of having the data domain on each processor start at 1.
  logical :: bathy_at_vel      ! If true, also define bathymetric fields at the
                               ! the velocity points.
  integer :: first_direction   ! An integer that indicates which direction is to be
                               ! updated first in directionally split parts of the
                               ! calculation.  This can be altered during the course
                               ! of the run via calls to set_first_direction.
  integer :: nkml, nkbl, verbosity, write_geom
  integer :: dynamics_stencil  ! The computational stencil for the calculations
                               ! in the dynamic core.

  type(time_type)                 :: Start_time
  type(ocean_internal_state)      :: MOM_internal_state
  character(len=200) :: area_varname, ice_shelf_file, inputdir, filename

  if (associated(CS)) then
    call MOM_error(WARNING, "initialize_MOM called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  if (test_grid_copy) then ; allocate(G)
  else ; G => CS%G ; endif

  CS%Time => Time

  id_clock_init = cpu_clock_id('Ocean Initialization', grain=CLOCK_SUBCOMPONENT)
  call cpu_clock_begin(id_clock_init)

  Start_time = Time ; if (present(Time_in)) Start_time = Time_in

  call Get_MOM_Input(param_file, dirs)

  verbosity = 2 ; call read_param(param_file, "VERBOSITY", verbosity)
  call MOM_set_verbosity(verbosity)
  call callTree_enter("initialize_MOM(), MOM.F90")

  call find_obsolete_params(param_file)

  ! Read relevant parameters and write them to the model log.
  call log_version(param_file, "MOM", version, "")
  call get_param(param_file, "MOM", "VERBOSITY", verbosity,  &
                 "Integer controlling level of messaging\n" // &
                 "\t0 = Only FATAL messages\n" // &
                 "\t2 = Only FATAL, WARNING, NOTE [default]\n" // &
                 "\t9 = All)", default=2)
  call get_param(param_file, "MOM", "DO_UNIT_TESTS", do_unit_tests, &
                 "If True, exercises unit tests at model start up.", &
                 default=.false.)
  if (do_unit_tests) then
    call unit_tests(verbosity)
  endif

  call get_param(param_file, "MOM", "SPLIT", CS%split, &
                 "Use the split time stepping if true.", default=.true.)
  if (CS%split) then
    CS%use_RK2 = .false.
  else
    call get_param(param_file, "MOM", "USE_RK2", CS%use_RK2, &
                 "If true, use RK2 instead of RK3 in the unsplit time stepping.", &
                 default=.false.)
  endif

  call get_param(param_file, "MOM", "CALC_RHO_FOR_SEA_LEVEL", CS%calc_rho_for_sea_lev, &
                 "If true, the in-situ density is used to calculate the\n"//&
                 "effective sea level that is returned to the coupler. If false,\n"//&
                 "the Boussinesq parameter RHO_0 is used.", default=.false.)
  call get_param(param_file, "MOM", "ENABLE_THERMODYNAMICS", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state \n"//&
                 "variables.", default=.true.)
  call get_param(param_file, "MOM", "USE_EOS", use_EOS, &
                 "If true,  density is calculated from temperature and \n"//&
                 "salinity with an equation of state.  If USE_EOS is \n"//&
                 "true, ENABLE_THERMODYNAMICS must be true as well.", &
                 default=CS%use_temperature)
  call get_param(param_file, "MOM", "DIABATIC_FIRST", CS%diabatic_first, &
                 "If true, apply diabatic and thermodynamic processes, \n"//&
                 "including buoyancy forcing and mass gain or loss, \n"//&
                 "before stepping the dynamics forward.", default=.false.)
  call get_param(param_file, "MOM", "USE_CONTEMP_ABSSAL", CS%use_conT_absS, &
                 "If true, , the prognostics T&S are the conservative temperature \n"//&
                 "and absolute salinity. Care should be taken to convert them \n"//&
                 "to potential temperature and practical salinity before  \n"//&
                 "exchanging them with the coupler and/or reporting T&S diagnostics. \n"&
                 , default=.false.)
  call get_param(param_file, "MOM", "ADIABATIC", CS%adiabatic, &
                 "There are no diapycnal mass fluxes if ADIABATIC is \n"//&
                 "true. This assumes that KD = KDML = 0.0 and that \n"//&
                 "there is no buoyancy forcing, but makes the model \n"//&
                 "faster by eliminating subroutine calls.", default=.false.)
  call get_param(param_file, "MOM", "DO_DYNAMICS", CS%do_dynamics, &
                 "If False, skips the dynamics calls that update u & v, as well as\n"//&
                 "the gravity wave adjustment to h. This is a fragile feature and\n"//&
                 "thus undocumented.", default=.true., do_not_log=.true. )
  call get_param(param_file, "MOM", "ADVECT_TS", CS%advect_TS , &
                 "If True, advect temperature and salinity horizontally\n"//&
                 "If False, T/S are registered for advection.\n"//&
                 "This is intended only to be used in offline tracer mode.", &
                 "and is by default false in that case", &
                 do_not_log = .true., default=.true. )
  if (present(offline_tracer_mode)) then ! Only read this parameter in solo mode
    call get_param(param_file, "MOM", "OFFLINE_TRACER_MODE", CS%offline_tracer_mode, &
                 "If true, barotropic and baroclinic dynamics, thermodynamics\n"//&
                 "are all bypassed with all the fields necessary to integrate\n"//&
                 "the tracer advection and diffusion equation are read in from\n"//&
                 "files stored from a previous integration of the prognostic model.\n"//&
                 "NOTE: This option only used in the ocean_solo_driver.", default=.false.)
    if(CS%offline_tracer_mode) then
      call get_param(param_file, "MOM", "ADVECT_TS", CS%advect_TS , &
                   "If True, advect temperature and salinity horizontally\n"//&
                   "If False, T/S are registered for advection.\n"//&
                   "This is intended only to be used in offline tracer mode."//&
                   "and is by default false in that case", &
                   default=.false. )
    endif
  endif
  call get_param(param_file, "MOM", "USE_REGRIDDING", CS%use_ALE_algorithm , &
                 "If True, use the ALE algorithm (regridding/remapping).\n"//&
                 "If False, use the layered isopycnal algorithm.", default=.false. )
  call get_param(param_file, "MOM", "BULKMIXEDLAYER", CS%bulkmixedlayer, &
                 "If true, use a Kraus-Turner-like bulk mixed layer \n"//&
                 "with transitional buffer layers.  Layers 1 through  \n"//&
                 "NKML+NKBL have variable densities. There must be at \n"//&
                 "least NKML+NKBL+1 layers if BULKMIXEDLAYER is true. \n"//&
                 "BULKMIXEDLAYER can not be used with USE_REGRIDDING. \n"//&
                 "The default is influenced by ENABLE_THERMODYNAMICS.", &
                 default=CS%use_temperature .and. .not.CS%use_ALE_algorithm)
  call get_param(param_file, "MOM", "THICKNESSDIFFUSE", CS%thickness_diffuse, &
                 "If true, interface heights are diffused with a \n"//&
                 "coefficient of KHTH.", default=.false.)
  call get_param(param_file, "MOM",  "THICKNESSDIFFUSE_FIRST", &
                                      CS%thickness_diffuse_first, &
                 "If true, do thickness diffusion before dynamics.\n"//&
                 "This is only used if THICKNESSDIFFUSE is true.", &
                 default=.false.)
  if (.not.CS%thickness_diffuse) CS%thickness_diffuse_first = .false.
  call get_param(param_file, "MOM", "BATHYMETRY_AT_VEL", bathy_at_vel, &
                 "If true, there are separate values for the basin depths \n"//&
                 "at velocity points.  Otherwise the effects of topography \n"//&
                 "are entirely determined from thickness points.", &
                 default=.false.)

  call get_param(param_file, "MOM", "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", default=.false.)
  call get_param(param_file, "MOM", "DEBUG_TRUNCATIONS", CS%debug_truncations, &
                 "If true, calculate all diagnostics that are useful for \n"//&
                 "debugging truncations.", default=.false.)

  call get_param(param_file, "MOM", "DT", CS%dt, &
                 "The (baroclinic) dynamics time step.  The time-step that \n"//&
                 "is actually used will be an integer fraction of the \n"//&
                 "forcing time-step (DT_FORCING in ocean-only mode or the \n"//&
                 "coupling timestep in coupled mode.)", units="s", &
                 fail_if_missing=.true.)
  call get_param(param_file, "MOM", "DT_THERM", CS%dt_therm, &
                 "The thermodynamic and tracer advection time step. \n"//&
                 "Ideally DT_THERM should be an integer multiple of DT \n"//&
                 "and less than the forcing or coupling time-step, unless \n"//&
                 "THERMO_SPANS_COUPLING is true, in which case DT_THERM \n"//&
                 "can be an integer multiple of the coupling timestep.  By \n"//&
                 "default DT_THERM is set to DT.", units="s", default=CS%dt)
  call get_param(param_file, "MOM", "THERMO_SPANS_COUPLING", CS%thermo_spans_coupling, &
                 "If true, the MOM will take thermodynamic and tracer \n"//&
                 "timesteps that can be longer than the coupling timestep. \n"//&
                 "The actual thermodynamic timestep that is used in this \n"//&
                 "case is the largest integer multiple of the coupling \n"//&
                 "timestep that is less than or equal to DT_THERM.", default=.false.)

  if (.not.CS%bulkmixedlayer) then
    call get_param(param_file, "MOM", "HMIX_SFC_PROP", CS%Hmix, &
                 "If BULKMIXEDLAYER is false, HMIX_SFC_PROP is the depth \n"//&
                 "over which to average to find surface properties like \n"//&
                 "SST and SSS or density (but not surface velocities).", &
                 units="m", default=1.0)
    call get_param(param_file, "MOM", "HMIX_UV_SFC_PROP", CS%Hmix_UV, &
                 "If BULKMIXEDLAYER is false, HMIX_UV_SFC_PROP is the depth\n"//&
                 "over which to average to find surface flow properties,\n"//&
                 "SSU, SSV. A non-positive value indicates no averaging.", &
                 units="m", default=0.)
  endif
  call get_param(param_file, "MOM", "MIN_Z_DIAG_INTERVAL", Z_diag_int, &
                 "The minimum amount of time in seconds between \n"//&
                 "calculations of depth-space diagnostics. Making this \n"//&
                 "larger than DT_THERM reduces the  performance penalty \n"//&
                 "of regridding to depth online.", units="s", default=0.0)
  call get_param(param_file, "MOM", "INTERPOLATE_P_SURF", CS%interp_p_surf, &
                 "If true, linearly interpolate the surface pressure \n"//&
                 "over the coupling time step, using the specified value \n"//&
                 "at the end of the step.", default=.false.)

  if (CS%split) then
    call get_param(param_file, "MOM", "DTBT", dtbt, default=-0.98)
    default_val = CS%dt_therm ; if (dtbt > 0.0) default_val = -1.0
    CS%dtbt_reset_period = -1.0
    call get_param(param_file, "MOM", "DTBT_RESET_PERIOD", CS%dtbt_reset_period, &
                 "The period between recalculations of DTBT (if DTBT <= 0). \n"//&
                 "If DTBT_RESET_PERIOD is negative, DTBT is set based \n"//&
                 "only on information available at initialization.  If \n"//&
                 "dynamic, DTBT will be set at least every forcing time \n"//&
                 "step, and if 0, every dynamics time step.  The default is \n"//&
                 "set by DT_THERM.  This is only used if SPLIT is true.", &
                 units="s", default=default_val, do_not_read=(dtbt > 0.0))
  endif

  ! This is here in case these values are used inappropriately.
  CS%use_frazil = .false. ; CS%bound_salinity = .false. ; CS%tv%P_Ref = 2.0e7
  if (CS%use_temperature) then
    call get_param(param_file, "MOM", "FRAZIL", CS%use_frazil, &
                 "If true, water freezes if it gets too cold, and the \n"//&
                 "the accumulated heat deficit is returned in the \n"//&
                 "surface state.  FRAZIL is only used if \n"//&
                 "ENABLE_THERMODYNAMICS is true.", default=.false.)
    call get_param(param_file, "MOM", "DO_GEOTHERMAL", use_geothermal, &
                 "If true, apply geothermal heating.", default=.false.)
    call get_param(param_file, "MOM", "BOUND_SALINITY", CS%bound_salinity, &
                 "If true, limit salinity to being positive. (The sea-ice \n"//&
                 "model may ask for more salt than is available and \n"//&
                 "drive the salinity negative otherwise.)", default=.false.)
    call get_param(param_file, "MOM", "C_P", CS%tv%C_p, &
                 "The heat capacity of sea water, approximated as a \n"//&
                 "constant. This is only used if ENABLE_THERMODYNAMICS is \n"//&
                 "true. The default value is from the TEOS-10 definition \n"//&
                 "of conservative temperature.", units="J kg-1 K-1", &
                 default=3991.86795711963)
  endif
  if (use_EOS) call get_param(param_file, "MOM", "P_REF", CS%tv%P_Ref, &
                 "The pressure that is used for calculating the coordinate \n"//&
                 "density.  (1 Pa = 1e4 dbar, so 2e7 is commonly used.) \n"//&
                 "This is only used if USE_EOS and ENABLE_THERMODYNAMICS \n"//&
                 "are true.", units="Pa", default=2.0e7)

  if (CS%bulkmixedlayer) then
    call get_param(param_file, "MOM", "NKML", nkml, &
                 "The number of sublayers within the mixed layer if \n"//&
                 "BULKMIXEDLAYER is true.", units="nondim", default=2)
    call get_param(param_file, "MOM", "NKBL", nkbl, &
                 "The number of layers that are used as variable density \n"//&
                 "buffer layers if BULKMIXEDLAYER is true.", units="nondim", &
                 default=2)
  endif

  call get_param(param_file, "MOM", "GLOBAL_INDEXING", global_indexing, &
                 "If true, use a global lateral indexing convention, so \n"//&
                 "that corresponding points on different processors have \n"//&
                 "the same index. This does not work with static memory.", &
                 default=.false., layoutParam=.true.)
#ifdef STATIC_MEMORY_
  if (global_indexing) call MOM_error(FATAL, "initialize_MOM: "//&
       "GLOBAL_INDEXING can not be true with STATIC_MEMORY.")
#endif
  call get_param(param_file, "MOM", "FIRST_DIRECTION", first_direction, &
                 "An integer that indicates which direction goes first \n"//&
                 "in parts of the code that use directionally split \n"//&
                 "updates, with even numbers (or 0) used for x- first \n"//&
                 "and odd numbers used for y-first.", default=0)

  call get_param(param_file, "MOM", "CHECK_BAD_SURFACE_VALS", &
                                     CS%check_bad_surface_vals, &
                 "If true, check the surface state for ridiculous values.", &
                 default=.false.)
  if (CS%check_bad_surface_vals) then
    call get_param(param_file, "MOM", "BAD_VAL_SSH_MAX", CS%bad_val_ssh_max, &
                 "The value of SSH above which a bad value message is \n"//&
                 "triggered, if CHECK_BAD_SURFACE_VALS is true.", units="m", &
                 default=20.0)
    call get_param(param_file, "MOM", "BAD_VAL_SSS_MAX", CS%bad_val_sss_max, &
                 "The value of SSS above which a bad value message is \n"//&
                 "triggered, if CHECK_BAD_SURFACE_VALS is true.", units="PPT", &
                 default=45.0)
    call get_param(param_file, "MOM", "BAD_VAL_SST_MAX", CS%bad_val_sst_max, &
                 "The value of SST above which a bad value message is \n"//&
                 "triggered, if CHECK_BAD_SURFACE_VALS is true.", &
                 units="deg C", default=45.0)
    call get_param(param_file, "MOM", "BAD_VAL_SST_MIN", CS%bad_val_sst_min, &
                 "The value of SST below which a bad value message is \n"//&
                 "triggered, if CHECK_BAD_SURFACE_VALS is true.", &
                 units="deg C", default=-2.1)
    call get_param(param_file, "MOM", "BAD_VAL_COLUMN_THICKNESS", CS%bad_val_column_thickness, &
         "The value of column thickness below which a bad value message is \n"//&
         "triggered, if CHECK_BAD_SURFACE_VALS is true.", units="m", &
                          default=0.0)
  endif

  call get_param(param_file, "MOM", "SAVE_INITIAL_CONDS", save_IC, &
                 "If true, write the initial conditions to a file given \n"//&
                 "by IC_OUTPUT_FILE.", default=.false.)
  call get_param(param_file, "MOM", "IC_OUTPUT_FILE", CS%IC_file, &
                 "The file into which to write the initial conditions.", &
                 default="MOM_IC")
  call get_param(param_file, "MOM", "WRITE_GEOM", write_geom, &
                 "If =0, never write the geometry and vertical grid files.\n"//&
                 "If =1, write the geometry and vertical grid files only for\n"//&
                 "a new simulation. If =2, always write the geometry and\n"//&
                 "vertical grid files. Other values are invalid.", default=1)
  if (write_geom<0 .or. write_geom>2) call MOM_error(FATAL,"MOM: "//&
         "WRITE_GEOM must be equal to 0, 1 or 2.")
  write_geom_files = ((write_geom==2) .or. ((write_geom==1) .and. &
     ((dirs%input_filename(1:1)=='n') .and. (LEN_TRIM(dirs%input_filename)==1))))
! If the restart file type had been initialized, this could become:
!  write_geom_files = ((write_geom==2) .or. &
!                      ((write_geom==1) .and. is_new_run(CS%restart_CSp)))

  ! Check for inconsistent parameter settings.
  if (CS%use_ALE_algorithm .and. CS%bulkmixedlayer) call MOM_error(FATAL, &
    "MOM: BULKMIXEDLAYER can not currently be used with the ALE algorithm.")
  if (CS%use_ALE_algorithm .and. .not.CS%use_temperature) call MOM_error(FATAL, &
     "MOM: At this time, USE_EOS should be True when using the ALE algorithm.")
  if (CS%adiabatic .and. CS%use_temperature) call MOM_error(WARNING, &
    "MOM: ADIABATIC and ENABLE_THERMODYNAMICS both defined is usually unwise.")
  if (use_EOS .and. .not.CS%use_temperature) call MOM_error(FATAL, &
    "MOM: ENABLE_THERMODYNAMICS must be defined to use USE_EOS.")
  if (CS%adiabatic .and. CS%bulkmixedlayer) call MOM_error(FATAL, &
    "MOM: ADIABATIC and BULKMIXEDLAYER can not both be defined.")
  if (CS%bulkmixedlayer .and. .not.use_EOS) call MOM_error(FATAL, &
      "initialize_MOM: A bulk mixed layer can only be used with T & S as "//&
      "state variables. Add USE_EOS = True to MOM_input.")

  call get_param(param_file, 'MOM', "ICE_SHELF", use_ice_shelf, default=.false., do_not_log=.true.)
  if (use_ice_shelf) then
     inputdir = "." ;  call get_param(param_file, 'MOM', "INPUTDIR", inputdir)
     inputdir = slasher(inputdir)
     call get_param(param_file, 'MOM', "ICE_THICKNESS_FILE", ice_shelf_file, &
                    "The file from which the ice bathymetry and area are read.", &
                    fail_if_missing=.true.)
     call get_param(param_file, 'MOM', "ICE_AREA_VARNAME", area_varname, &
                    "The name of the area variable in ICE_THICKNESS_FILE.", &
                    fail_if_missing=.true.)
  endif


  call callTree_waypoint("MOM parameters read (initialize_MOM)")

  ! Set up the model domain and grids.
#ifdef SYMMETRIC_MEMORY_
  symmetric = .true.
#else
  symmetric = .false.
#endif
#ifdef STATIC_MEMORY_
  call MOM_domains_init(G%domain, param_file, symmetric=symmetric, &
            static_memory=.true., NIHALO=NIHALO_, NJHALO=NJHALO_, &
            NIGLOBAL=NIGLOBAL_, NJGLOBAL=NJGLOBAL_, NIPROC=NIPROC_, &
            NJPROC=NJPROC_)
#else
  call MOM_domains_init(G%domain, param_file, symmetric=symmetric)
#endif
  call callTree_waypoint("domains initialized (initialize_MOM)")

  call MOM_debugging_init(param_file)
  call diag_mediator_infrastructure_init()
  call MOM_io_init(param_file)

  call hor_index_init(G%Domain, HI, param_file, &
                      local_indexing=.not.global_indexing)

  call create_dyn_horgrid(dG, HI, bathymetry_at_vel=bathy_at_vel)
  call clone_MOM_domain(G%Domain, dG%Domain)

  call verticalGridInit( param_file, CS%GV )
  GV => CS%GV
!  dG%g_Earth = GV%g_Earth

  ! Allocate the auxiliary non-symmetric domain for debugging or I/O purposes.
  if (CS%debug .or. dG%symmetric) &
    call clone_MOM_domain(dG%Domain, dG%Domain_aux, symmetric=.false.)

  call callTree_waypoint("grids initialized (initialize_MOM)")


  call MOM_timing_init(CS)

  ! Allocate initialize time-invariant MOM variables.
  call MOM_initialize_fixed(dG, CS%OBC, param_file, write_geom_files, dirs%output_directory)
  call callTree_waypoint("returned from MOM_initialize_fixed() (initialize_MOM)")
  if (associated(CS%OBC)) call call_OBC_register(param_file, CS%update_OBC_CSp, CS%OBC)

  call tracer_registry_init(param_file, CS%tracer_Reg)

  ! Allocate and initialize space for the primary time-varying MOM variables.
  is   = dG%isc   ; ie   = dG%iec  ; js   = dG%jsc  ; je   = dG%jec ; nz = GV%ke
  isd  = dG%isd   ; ied  = dG%ied  ; jsd  = dG%jsd  ; jed  = dG%jed
  IsdB = dG%IsdB  ; IedB = dG%IedB ; JsdB = dG%JsdB ; JedB = dG%JedB
  ALLOC_(CS%u(IsdB:IedB,jsd:jed,nz))   ; CS%u(:,:,:) = 0.0
  ALLOC_(CS%v(isd:ied,JsdB:JedB,nz))   ; CS%v(:,:,:) = 0.0
  ALLOC_(CS%h(isd:ied,jsd:jed,nz))     ; CS%h(:,:,:) = GV%Angstrom
  ALLOC_(CS%uh(IsdB:IedB,jsd:jed,nz))  ; CS%uh(:,:,:) = 0.0
  ALLOC_(CS%vh(isd:ied,JsdB:JedB,nz))  ; CS%vh(:,:,:) = 0.0
  if (CS%use_temperature) then
    ALLOC_(CS%T(isd:ied,jsd:jed,nz))   ; CS%T(:,:,:) = 0.0
    ALLOC_(CS%S(isd:ied,jsd:jed,nz))   ; CS%S(:,:,:) = 0.0
    CS%tv%T => CS%T ; CS%tv%S => CS%S
    CS%vd_T = var_desc(name="T",units="degC",longname="Potential Temperature", &
                       cmor_field_name="thetao",                               &
                       conversion=CS%tv%C_p)
    CS%vd_S = var_desc(name="S",units="psu",longname="Salinity",&
                       cmor_field_name="so",                    &
                       conversion=0.001)
    if(CS%advect_TS) then
      call register_tracer(CS%tv%T, CS%vd_T, param_file, dG%HI, GV, CS%tracer_Reg, CS%vd_T)
      call register_tracer(CS%tv%S, CS%vd_S, param_file, dG%HI, GV, CS%tracer_Reg, CS%vd_S)
    endif
    if (associated(CS%OBC)) &
      call register_temp_salt_segments(GV, CS%OBC, CS%tv, CS%vd_T, CS%vd_S, param_file)
  endif
  if (CS%use_frazil) then
    allocate(CS%tv%frazil(isd:ied,jsd:jed)) ; CS%tv%frazil(:,:) = 0.0
  endif
  if (CS%bound_salinity) then
    allocate(CS%tv%salt_deficit(isd:ied,jsd:jed)) ; CS%tv%salt_deficit(:,:)=0.0
  endif

  if (CS%bulkmixedlayer .or. CS%use_temperature) then
    allocate(CS%Hml(isd:ied,jsd:jed)) ; CS%Hml(:,:) = 0.0
  endif

  if (CS%bulkmixedlayer) then
    GV%nkml = nkml ; GV%nk_rho_varies = nkml + nkbl
  else
    GV%nkml = 0 ; GV%nk_rho_varies = 0
  endif
  if (CS%use_ALE_algorithm) then
    call get_param(param_file, "MOM", "NK_RHO_VARIES", GV%nk_rho_varies, default=0) ! Will default to nz later... -AJA
  endif

  ALLOC_(CS%uhtr(IsdB:IedB,jsd:jed,nz)) ; CS%uhtr(:,:,:) = 0.0
  ALLOC_(CS%vhtr(isd:ied,JsdB:JedB,nz)) ; CS%vhtr(:,:,:) = 0.0
  CS%t_dyn_rel_adv = 0.0 ; CS%t_dyn_rel_thermo = 0.0 ; CS%t_dyn_rel_diag = 0.0

  if (CS%debug_truncations) then
    allocate(CS%u_prev(IsdB:IedB,jsd:jed,nz)) ; CS%u_prev(:,:,:) = 0.0
    allocate(CS%v_prev(isd:ied,JsdB:JedB,nz)) ; CS%v_prev(:,:,:) = 0.0
  endif

  MOM_internal_state%u => CS%u ; MOM_internal_state%v => CS%v
  MOM_internal_state%h => CS%h
  MOM_internal_state%uh => CS%uh ; MOM_internal_state%vh => CS%vh
  if (CS%use_temperature) then
    MOM_internal_state%T => CS%T ; MOM_internal_state%S => CS%S
  endif

  CS%CDp%uh => CS%uh ; CS%CDp%vh => CS%vh

  if (CS%interp_p_surf) then
    allocate(CS%p_surf_prev(isd:ied,jsd:jed)) ; CS%p_surf_prev(:,:) = 0.0
  endif

  ALLOC_(CS%ave_ssh(isd:ied,jsd:jed)) ; CS%ave_ssh(:,:) = 0.0

  ! Use the Wright equation of state by default, unless otherwise specified
  ! Note: this line and the following block ought to be in a separate
  ! initialization routine for tv.
  if (use_EOS) call EOS_init(param_file, CS%tv%eqn_of_state)
  if (CS%use_temperature) then
    allocate(CS%tv%TempxPmE(isd:ied,jsd:jed))
    CS%tv%TempxPmE(:,:) = 0.0
    if (use_geothermal) then
      allocate(CS%tv%internal_heat(isd:ied,jsd:jed))
      CS%tv%internal_heat(:,:) = 0.0
    endif
  endif
  call callTree_waypoint("state variables allocated (initialize_MOM)")

  ! Set the fields that are needed for bitwise identical restarting
  ! the time stepping scheme.
  call restart_init(param_file, CS%restart_CSp)
  call set_restart_fields(GV, param_file, CS)
  if (CS%split) then
    call register_restarts_dyn_split_RK2(dG%HI, GV, param_file, &
             CS%dyn_split_RK2_CSp, CS%restart_CSp, CS%uh, CS%vh)
  elseif (CS%use_RK2) then
    call register_restarts_dyn_unsplit_RK2(dG%HI, GV, param_file, &
           CS%dyn_unsplit_RK2_CSp, CS%restart_CSp)
  else
    call register_restarts_dyn_unsplit(dG%HI, GV, param_file, &
           CS%dyn_unsplit_CSp, CS%restart_CSp)
  endif

  ! This subroutine calls user-specified tracer registration routines.
  ! Additional calls can be added to MOM_tracer_flow_control.F90.
  call call_tracer_register(dG%HI, GV, param_file, CS%tracer_flow_CSp, &
                            CS%tracer_Reg, CS%restart_CSp)

  call MEKE_alloc_register_restart(dG%HI, param_file, CS%MEKE, CS%restart_CSp)
  call set_visc_register_restarts(dG%HI, GV, param_file, CS%visc, CS%restart_CSp)
  call mixedlayer_restrat_register_restarts(dG%HI, param_file, CS%mixedlayer_restrat_CSp, CS%restart_CSp)

  call callTree_waypoint("restart registration complete (initialize_MOM)")

  ! Initialize dynamically evolving fields, perhaps from restart files.
  call cpu_clock_begin(id_clock_MOM_init)
  call MOM_initialize_coord(GV, param_file, write_geom_files, &
                            dirs%output_directory, CS%tv, dG%max_depth)
  call callTree_waypoint("returned from MOM_initialize_coord() (initialize_MOM)")

  if (CS%use_ALE_algorithm) then
    call ALE_init(param_file, GV, dG%max_depth, CS%ALE_CSp)
    call callTree_waypoint("returned from ALE_init() (initialize_MOM)")
  endif

  !   Shift from using the temporary dynamic grid type to using the final
  ! (potentially static) ocean-specific grid type.
  !   The next line would be needed if G%Domain had not already been init'd above:
  !     call clone_MOM_domain(dG%Domain, G%Domain)
  call MOM_grid_init(G, param_file, HI, bathymetry_at_vel=bathy_at_vel)
  call copy_dyngrid_to_MOM_grid(dG, G)
  call destroy_dyn_horgrid(dG)

  ! Set a few remaining fields that are specific to the ocean grid type.
  call set_first_direction(G, first_direction)
  ! Allocate the auxiliary non-symmetric domain for debugging or I/O purposes.
  if (CS%debug .or. G%symmetric) &
    call clone_MOM_domain(G%Domain, G%Domain_aux, symmetric=.false.)
  ! Copy common variables from the vertical grid to the horizontal grid.
  ! Consider removing this later?
  G%ke = GV%ke ; G%g_Earth = GV%g_Earth

  call MOM_initialize_state(CS%u, CS%v, CS%h, CS%tv, Time, G, GV, param_file, &
                            dirs, CS%restart_CSp, CS%ALE_CSp, CS%tracer_Reg, &
                            CS%sponge_CSp, CS%ALE_sponge_CSp, CS%OBC, Time_in)
  call cpu_clock_end(id_clock_MOM_init)
  call callTree_waypoint("returned from MOM_initialize_state() (initialize_MOM)")

  ! From this point, there may be pointers being set, so the final grid type
  ! that will persist throughout the run has to be used.

  if (test_grid_copy) then
    !  Copy the data from the temporary grid to the dyn_hor_grid to CS%G.
    call create_dyn_horgrid(dG, G%HI)
    call clone_MOM_domain(G%Domain, dG%Domain)

    call clone_MOM_domain(G%Domain, CS%G%Domain)
    call MOM_grid_init(CS%G, param_file)

    call copy_MOM_grid_to_dyngrid(G, dg)
    call copy_dyngrid_to_MOM_grid(dg, CS%G)

    call destroy_dyn_horgrid(dG)
    call MOM_grid_end(G) ; deallocate(G)

    G => CS%G
    if (CS%debug .or. CS%G%symmetric) &
      call clone_MOM_domain(CS%G%Domain, CS%G%Domain_aux, symmetric=.false.)
    G%ke = GV%ke ; G%g_Earth = GV%g_Earth
  endif


  ! At this point, all user-modified initialization code has been called.  The
  ! remainder of this subroutine is controlled by the parameters that have
  ! have already been set.

  if (ALE_remap_init_conds(CS%ALE_CSp) .and. .not. query_initialized(CS%h,"h",CS%restart_CSp)) then
    ! This block is controlled by the ALE parameter REMAP_AFTER_INITIALIZATION.
    ! \todo This block exists for legacy reasons and we should phase it out of
    ! all examples. !###
    if (CS%debug) then
      call uvchksum("Pre ALE adjust init cond [uv]", &
                    CS%u, CS%v, G%HI, haloshift=1)
      call hchksum(CS%h,"Pre ALE adjust init cond h", G%HI, haloshift=1, scale=GV%H_to_m)
    endif
    call callTree_waypoint("Calling adjustGridForIntegrity() to remap initial conditions (initialize_MOM)")
    call adjustGridForIntegrity(CS%ALE_CSp, G, GV, CS%h )
    call callTree_waypoint("Calling ALE_main() to remap initial conditions (initialize_MOM)")
    if (use_ice_shelf) then
      filename = trim(inputdir)//trim(ice_shelf_file)
      if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
        "MOM: Unable to open "//trim(filename))

      allocate(area_shelf_h(isd:ied,jsd:jed))
      allocate(frac_shelf_h(isd:ied,jsd:jed))
      call read_data(filename,trim(area_varname),area_shelf_h,domain=G%Domain%mpp_domain)
      ! initialize frac_shelf_h with zeros (open water everywhere)
      frac_shelf_h(:,:) = 0.0
      ! compute fractional ice shelf coverage of h
      do j=jsd,jed ; do i=isd,ied
        if (G%areaT(i,j) > 0.0) &
          frac_shelf_h(i,j) = area_shelf_h(i,j) / G%areaT(i,j)
      enddo ; enddo
      ! pass to the pointer
      shelf_area => frac_shelf_h
      call ALE_main(G, GV, CS%h, CS%u, CS%v, CS%tv, CS%tracer_Reg, CS%ALE_CSp, &
                    frac_shelf_h = shelf_area)
    else
      call ALE_main( G, GV, CS%h, CS%u, CS%v, CS%tv, CS%tracer_Reg, CS%ALE_CSp )
    endif

    call cpu_clock_begin(id_clock_pass_init)
    call create_group_pass(tmp_pass_uv_T_S_h, CS%u, CS%v, G%Domain)
    if (CS%use_temperature) then
      call create_group_pass(tmp_pass_uv_T_S_h, CS%tv%T, G%Domain, halo=1)
      call create_group_pass(tmp_pass_uv_T_S_h, CS%tv%S, G%Domain, halo=1)
    endif
    call create_group_pass(tmp_pass_uv_T_S_h, CS%h, G%Domain, halo=1)
    call do_group_pass(tmp_pass_uv_T_S_h, G%Domain)
    call cpu_clock_end(id_clock_pass_init)

    if (CS%debug) then
      call uvchksum("Post ALE adjust init cond [uv]", CS%u, CS%v, G%HI, haloshift=1)
      call hchksum(CS%h, "Post ALE adjust init cond h", G%HI, haloshift=1, scale=GV%H_to_m)
    endif
  endif
  if ( CS%use_ALE_algorithm ) call ALE_updateVerticalGridType( CS%ALE_CSp, GV )

  diag => CS%diag
  ! Initialize the diag mediator.
  call diag_mediator_init(G, GV%ke, param_file, diag, doc_file_dir=dirs%output_directory)

  ! Initialize the diagnostics masks for native arrays.
  ! This step has to be done after call to MOM_initialize_state
  ! and before MOM_diagnostics_init
  call diag_masks_set(G, GV%ke, diag)

  ! Set up pointers within diag mediator control structure,
  ! this needs to occur _after_ CS%h etc. have been allocated.
  call diag_set_state_ptrs(CS%h, CS%T, CS%S, CS%tv%eqn_of_state, diag)

  ! This call sets up the diagnostic axes. These are needed,
  ! e.g. to generate the target grids below.
  call set_axes_info(G, GV, param_file, diag)

  ! Whenever thickness/T/S changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  ! FIXME: are h, T, S updated at the same time? Review these for T, S updates.
  call diag_update_remap_grids(diag)

  ! Calculate masks for diagnostics arrays in non-native coordinates
  ! This step has to be done after set_axes_info() because the axes needed
  ! to be configured, and after diag_update_remap_grids() because the grids
  ! must be defined.
  call set_masks_for_axes(G, diag)

  ! Diagnose static fields AND associate areas/volumes with axes
  call write_static_fields(G, CS%diag)
  call callTree_waypoint("static fields written (initialize_MOM)")

  ! Register the volume cell measure (must be one of first diagnostics)
  call register_cell_measure(G, CS%diag, Time)

  call cpu_clock_begin(id_clock_MOM_init)
  if (CS%use_ALE_algorithm) then
    call ALE_writeCoordinateFile( CS%ALE_CSp, GV, dirs%output_directory )
  endif
  call cpu_clock_end(id_clock_MOM_init)
  call callTree_waypoint("ALE initialized (initialize_MOM)")

  CS%useMEKE = MEKE_init(Time, G, param_file, diag, CS%MEKE_CSp, CS%MEKE, CS%restart_CSp)

  call VarMix_init(Time, G, param_file, diag, CS%VarMix)
  call set_visc_init(Time, G, GV, param_file, diag, CS%visc, CS%set_visc_CSp,CS%OBC)
  if (CS%split) then
    allocate(eta(SZI_(G),SZJ_(G))) ; eta(:,:) = 0.0
    call initialize_dyn_split_RK2(CS%u, CS%v, CS%h, CS%uh, CS%vh, eta, Time, &
              G, GV, param_file, diag, CS%dyn_split_RK2_CSp, CS%restart_CSp, &
              CS%dt, CS%ADp, CS%CDp, MOM_internal_state, CS%VarMix, CS%MEKE, &
              CS%OBC, CS%update_OBC_CSp, CS%ALE_CSp, CS%set_visc_CSp,        &
              CS%visc, dirs, CS%ntrunc)
  elseif (CS%use_RK2) then
    call initialize_dyn_unsplit_RK2(CS%u, CS%v, CS%h, Time, G, GV,         &
            param_file, diag, CS%dyn_unsplit_RK2_CSp, CS%restart_CSp,      &
            CS%ADp, CS%CDp, MOM_internal_state, CS%OBC, CS%update_OBC_CSp, &
            CS%ALE_CSp, CS%set_visc_CSp, CS%visc, dirs, CS%ntrunc)
  else
    call initialize_dyn_unsplit(CS%u, CS%v, CS%h, Time, G, GV,             &
            param_file, diag, CS%dyn_unsplit_CSp, CS%restart_CSp,          &
            CS%ADp, CS%CDp, MOM_internal_state, CS%OBC, CS%update_OBC_CSp, &
            CS%ALE_CSp, CS%set_visc_CSp, CS%visc, dirs, CS%ntrunc)
  endif
  call callTree_waypoint("dynamics initialized (initialize_MOM)")

  call thickness_diffuse_init(Time, G, GV, param_file, diag, CS%CDp, CS%thickness_diffuse_CSp)
  CS%mixedlayer_restrat = mixedlayer_restrat_init(Time, G, GV, param_file, diag, &
                                                  CS%mixedlayer_restrat_CSp)
  if (CS%mixedlayer_restrat) then
    if (.not.(CS%bulkmixedlayer .or. CS%use_ALE_algorithm)) &
      call MOM_error(FATAL, "MOM: MIXEDLAYER_RESTRAT true requires a boundary layer scheme.")
    ! When DIABATIC_FIRST=False and using CS%visc%ML in mixedlayer_restrat we need to update after a restart
    if (.not. CS%diabatic_first .and. associated(CS%visc%MLD)) &
      call pass_var(CS%visc%MLD, G%domain, halo=1)
  endif

  call MOM_diagnostics_init(MOM_internal_state, CS%ADp, CS%CDp, Time, G, GV, &
                            param_file, diag, CS%diagnostics_CSp)

  CS%Z_diag_interval = set_time(int((CS%dt_therm) * &
       max(1,floor(0.01 + Z_diag_int/(CS%dt_therm)))))
  call MOM_diag_to_Z_init(Time, G, GV, param_file, diag, CS%diag_to_Z_CSp)
  CS%Z_diag_time = Start_time + CS%Z_diag_interval * (1 + &
    ((Time + set_time(int(CS%dt_therm))) - Start_time) / CS%Z_diag_interval)

  if (associated(CS%sponge_CSp)) &
    call init_sponge_diags(Time, G, diag, CS%sponge_CSp)

  if (associated(CS%ALE_sponge_CSp)) &
    call init_ALE_sponge_diags(Time, G, diag, CS%ALE_sponge_CSp)

  if (CS%adiabatic) then
    call adiabatic_driver_init(Time, G, param_file, diag, CS%diabatic_CSp, &
                               CS%tracer_flow_CSp, CS%diag_to_Z_CSp)
  else
    call diabatic_driver_init(Time, G, GV, param_file, CS%use_ALE_algorithm, diag,     &
                              CS%ADp, CS%CDp, CS%diabatic_CSp, CS%tracer_flow_CSp, &
                              CS%sponge_CSp, CS%ALE_sponge_CSp, CS%diag_to_Z_CSp)
  endif

  call tracer_advect_init(Time, G, param_file, diag, CS%tracer_adv_CSp)
  call tracer_hor_diff_init(Time, G, param_file, diag, CS%tracer_diff_CSp, CS%neutral_diffusion_CSp)

  if (CS%use_ALE_algorithm) &
    call register_diags_TS_vardec(Time, G%HI, GV, param_file, CS)

  call lock_tracer_registry(CS%tracer_Reg)
  call callTree_waypoint("tracer registry now locked (initialize_MOM)")

  ! now register some diagnostics since tracer registry is locked
  call register_diags(Time, G, GV, CS, CS%ADp, CS%tv%C_p)
  call register_diags_TS_tendency(Time, G, CS)
  if (CS%use_ALE_algorithm) then
    call ALE_register_diags(Time, G, diag, CS%tv%C_p, CS%tracer_Reg, CS%ALE_CSp)
  endif



  ! If need a diagnostic field, then would have been allocated in register_diags.
  if (CS%use_temperature) then
    if(CS%advect_TS) then
      call add_tracer_diagnostics("T", CS%tracer_Reg, CS%T_adx, CS%T_ady, &
                        CS%T_diffx, CS%T_diffy, CS%T_adx_2d, CS%T_ady_2d, &
                        CS%T_diffx_2d, CS%T_diffy_2d, CS%T_advection_xy)
      call add_tracer_diagnostics("S", CS%tracer_Reg, CS%S_adx, CS%S_ady, &
                        CS%S_diffx, CS%S_diffy, CS%S_adx_2d, CS%S_ady_2d, &
                        CS%S_diffx_2d, CS%S_diffy_2d, CS%S_advection_xy)
    endif
    call register_Z_tracer(CS%tv%T, "temp", "Potential Temperature", "degC", Time,   &
                      G, CS%diag_to_Z_CSp, cmor_field_name="thetao",                 &
                      cmor_standard_name="sea_water_potential_temperature",          &
                      cmor_long_name ="Sea Water Potential Temperature")
    call register_Z_tracer(CS%tv%S, "salt", "Salinity", "psu", Time,               &
                      G, CS%diag_to_Z_CSp, cmor_field_name="so",                   &
                      cmor_standard_name="sea_water_salinity",                     &
                      cmor_long_name ="Sea Water Salinity")
  endif

  ! This subroutine initializes any tracer packages.
  new_sim = is_new_run(CS%restart_CSp)
  call tracer_flow_control_init(.not.new_sim, Time, G, GV, CS%h, param_file, &
             CS%diag, CS%OBC, CS%tracer_flow_CSp, CS%sponge_CSp, &
             CS%ALE_sponge_CSp, CS%diag_to_Z_CSp, CS%tv)


  ! If running in offline tracer mode, initialize the necessary control structure and
  ! parameters
  if(present(offline_tracer_mode)) offline_tracer_mode=CS%offline_tracer_mode

  if(CS%offline_tracer_mode) then
    ! Setup some initial parameterizations and also assign some of the subtypes
    call offline_transport_init(param_file, CS%offline_CSp, CS%diabatic_CSp, G, GV)
    call insert_offline_main( CS=CS%offline_CSp, ALE_CSp=CS%ALE_CSp, diabatic_CSp=CS%diabatic_CSp, &
                              diag=CS%diag, OBC=CS%OBC, tracer_adv_CSp=CS%tracer_adv_CSp,              &
                              tracer_flow_CSp=CS%tracer_flow_CSp, tracer_Reg=CS%tracer_Reg,            &
                              tv=CS%tv, x_before_y = (MOD(first_direction,2)==0), debug=CS%debug )
    call register_diags_offline_transport(Time, CS%diag, CS%offline_CSp)
  endif

  !--- set up group pass for u,v,T,S and h. pass_uv_T_S_h also is used in step_MOM
  call cpu_clock_begin(id_clock_pass_init)
  dynamics_stencil = min(3, G%Domain%nihalo, G%Domain%njhalo)
  call create_group_pass(CS%pass_uv_T_S_h, CS%u, CS%v, G%Domain, halo=dynamics_stencil)
  if (CS%use_temperature) then
    call create_group_pass(CS%pass_uv_T_S_h, CS%tv%T, G%Domain, halo=dynamics_stencil)
    call create_group_pass(CS%pass_uv_T_S_h, CS%tv%S, G%Domain, halo=dynamics_stencil)
  endif
  call create_group_pass(CS%pass_uv_T_S_h, CS%h, G%Domain, halo=dynamics_stencil)

  call do_group_pass(CS%pass_uv_T_S_h, G%Domain)
  call cpu_clock_end(id_clock_pass_init)

  call register_obsolete_diagnostics(param_file, CS%diag)
  call neutral_diffusion_diag_init(Time, G, diag, CS%tv%C_p, CS%tracer_Reg, CS%neutral_diffusion_CSp)

  if (CS%use_frazil) then
    if (.not.query_initialized(CS%tv%frazil,"frazil",CS%restart_CSp)) &
      CS%tv%frazil(:,:) = 0.0
  endif

  if (CS%interp_p_surf) then
    CS%p_surf_prev_set = &
      query_initialized(CS%p_surf_prev,"p_surf_prev",CS%restart_CSp)

    if (CS%p_surf_prev_set) call pass_var(CS%p_surf_prev, G%domain)
  endif

  if (.not.query_initialized(CS%ave_ssh,"ave_ssh",CS%restart_CSp)) then
    if (CS%split) then
      call find_eta(CS%h, CS%tv, GV%g_Earth, G, GV, CS%ave_ssh, eta)
    else
      call find_eta(CS%h, CS%tv, GV%g_Earth, G, GV, CS%ave_ssh)
    endif
  endif
  if (CS%split) deallocate(eta)

  ! Flag whether to save initial conditions in finish_MOM_initialization() or not.
  CS%write_IC = save_IC .and. &
                .not.((dirs%input_filename(1:1) == 'r') .and. &
                      (LEN_TRIM(dirs%input_filename) == 1))

  call callTree_leave("initialize_MOM()")
  call cpu_clock_end(id_clock_init)

end subroutine initialize_MOM

!> This subroutine finishes initializing MOM and writes out the initial conditions.
subroutine finish_MOM_initialization(Time, dirs, CS, fluxes)
  type(time_type),           intent(in)    :: Time        !< model time, used in this routine
  type(directories),         intent(in)    :: dirs        !< structure with directory paths
  type(MOM_control_struct),  pointer       :: CS          !< pointer set in this routine to MOM control structure
  type(forcing),             intent(inout) :: fluxes      !< pointers to forcing fields
  ! Local variables
  type(ocean_grid_type), pointer :: G => NULL()
  type(verticalGrid_type), pointer :: GV => NULL()
  type(MOM_restart_CS), pointer :: restart_CSp_tmp => NULL()
  real, allocatable :: z_interface(:,:,:) ! Interface heights (meter)
  real, allocatable :: eta(:,:) ! Interface heights (meter)
  type(vardesc) :: vd

  call cpu_clock_begin(id_clock_init)
  call callTree_enter("finish_MOM_initialization()")

  ! Pointers for convenience
  G => CS%G ; GV => CS%GV

  ! Write initial conditions
  if (CS%write_IC) then
    allocate(restart_CSp_tmp)
    restart_CSp_tmp = CS%restart_CSp
    allocate(z_interface(SZI_(G),SZJ_(G),SZK_(G)+1))
    call find_eta(CS%h, CS%tv, GV%g_Earth, G, GV, z_interface)
    vd = var_desc("eta","meter","Interface heights",z_grid='i')
    call register_restart_field(z_interface, vd, .true., restart_CSp_tmp)

    call save_restart(dirs%output_directory, Time, G, &
                      restart_CSp_tmp, filename=CS%IC_file, GV=GV)
    deallocate(z_interface)
    deallocate(restart_CSp_tmp)
  endif

  call callTree_leave("finish_MOM_initialization()")
  call cpu_clock_end(id_clock_init)

end subroutine finish_MOM_initialization

!> Register the diagnostics
subroutine register_diags(Time, G, GV, CS, ADp, C_p)
  type(time_type),           intent(in)    :: Time  !< current model time
  type(ocean_grid_type),     intent(inout) :: G     !< ocean grid structu
  type(verticalGrid_type),   intent(inout) :: GV    !< ocean vertical grid structure
  type(MOM_control_struct),  pointer       :: CS    !< control structure set up by initialize_MOM
  type(accel_diag_ptrs),     intent(inout) :: ADp   !< structure pointing to accelerations in momentum equation
  real,                      intent(in)    :: C_p   !< Heat capacity used in conversion to watts

  real :: conv2watt
  character(len=48) :: thickness_units, flux_units, S_flux_units
  type(diag_ctrl), pointer :: diag
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed ; nz = G%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  diag => CS%diag

  thickness_units = get_thickness_units(GV)
  flux_units      = get_flux_units(GV)
  S_flux_units    = get_tr_flux_units(GV, "psu")
  conv2watt       = GV%H_to_kg_m2 * C_p

  !Initialize the diagnostics mask arrays.
  !This has to be done after MOM_initialize_state call.
  !call diag_masks_set(G, CS%missing)

  CS%id_u = register_diag_field('ocean_model', 'u', diag%axesCuL, Time,              &
      'Zonal velocity', 'm s-1', cmor_field_name='uo', &
      cmor_standard_name='sea_water_x_velocity', cmor_long_name='Sea Water X Velocity')
  CS%id_v = register_diag_field('ocean_model', 'v', diag%axesCvL, Time,                  &
      'Meridional velocity', 'm s-1', cmor_field_name='vo', &
      cmor_standard_name='sea_water_y_velocity', cmor_long_name='Sea Water Y Velocity')
  CS%id_h = register_diag_field('ocean_model', 'h', diag%axesTL, Time, &
      'Layer Thickness', thickness_units, v_extensive=.true.)

  CS%id_volo = register_scalar_field('ocean_model', 'volo', Time, diag,&
      long_name='Total volume of liquid ocean', units='m3',            &
      standard_name='sea_water_volume')
  CS%id_zos = register_diag_field('ocean_model', 'zos', diag%axesT1, Time,&
      standard_name = 'sea_surface_height_above_geoid',                   &
      long_name= 'Sea surface height above geoid', units='m', missing_value=CS%missing)
  CS%id_zossq = register_diag_field('ocean_model', 'zossq', diag%axesT1, Time,&
      standard_name='square_of_sea_surface_height_above_geoid',             &
      long_name='Square of sea surface height above geoid', units='m2', missing_value=CS%missing)
  CS%id_ssh = register_diag_field('ocean_model', 'SSH', diag%axesT1, Time, &
      'Sea Surface Height', 'm', CS%missing)
  CS%id_ssh_ga = register_scalar_field('ocean_model', 'ssh_ga', Time, diag,&
      long_name='Area averaged sea surface height', units='m',            &
      standard_name='area_averaged_sea_surface_height')
  CS%id_ssh_inst = register_diag_field('ocean_model', 'SSH_inst', diag%axesT1, Time, &
      'Instantaneous Sea Surface Height', 'm', CS%missing)
  CS%id_ssu = register_diag_field('ocean_model', 'SSU', diag%axesCu1, Time, &
      'Sea Surface Zonal Velocity', 'm s-1', CS%missing)
  CS%id_ssv = register_diag_field('ocean_model', 'SSV', diag%axesCv1, Time, &
      'Sea Surface Meridional Velocity', 'm s-1', CS%missing)
  CS%id_speed = register_diag_field('ocean_model', 'speed', diag%axesT1, Time, &
      'Sea Surface Speed', 'm s-1', CS%missing)

  if (CS%use_temperature) then
    CS%id_T = register_diag_field('ocean_model', 'temp', diag%axesTL, Time, &
        'Potential Temperature', 'degC',                                    &
         cmor_field_name="thetao",                                          &
         cmor_standard_name="sea_water_potential_temperature",              &
         cmor_long_name ="Sea Water Potential Temperature")
    CS%id_S = register_diag_field('ocean_model', 'salt', diag%axesTL, Time, &
        long_name='Salinity', units='psu', cmor_field_name='so',            &
        cmor_long_name='Sea Water Salinity',                                &
        cmor_standard_name='sea_water_salinity')
    CS%id_tob = register_diag_field('ocean_model','tob', diag%axesT1, Time,          &
        long_name='Sea Water Potential Temperature at Sea Floor',                    &
        standard_name='sea_water_potential_temperature_at_sea_floor', units='degC')
    CS%id_sob = register_diag_field('ocean_model','sob',diag%axesT1, Time,           &
        long_name='Sea Water Salinity at Sea Floor',                                 &
        standard_name='sea_water_salinity_at_sea_floor', units='psu')
    CS%id_sst = register_diag_field('ocean_model', 'SST', diag%axesT1, Time,     &
        'Sea Surface Temperature', 'degC', CS%missing, cmor_field_name='tos', &
        cmor_long_name='Sea Surface Temperature',                                &
        cmor_standard_name='sea_surface_temperature')
    CS%id_sst_sq = register_diag_field('ocean_model', 'SST_sq', diag%axesT1, Time, &
        'Sea Surface Temperature Squared', 'degC2', CS%missing, cmor_field_name='tossq', &
        cmor_long_name='Square of Sea Surface Temperature ',                      &
        cmor_standard_name='square_of_sea_surface_temperature')
    CS%id_sss = register_diag_field('ocean_model', 'SSS', diag%axesT1, Time, &
        'Sea Surface Salinity', 'psu', CS%missing, cmor_field_name='sos', &
        cmor_long_name='Sea Surface Salinity',                            &
        cmor_standard_name='sea_surface_salinity')
    CS%id_sss_sq = register_diag_field('ocean_model', 'SSS_sq', diag%axesT1, Time, &
        'Sea Surface Salinity Squared', 'psu', CS%missing, cmor_field_name='sossq', &
        cmor_long_name='Square of Sea Surface Salinity ',                     &
        cmor_standard_name='square_of_sea_surface_salinity')
    if (CS%use_conT_absS) then
      CS%id_Tcon = register_diag_field('ocean_model', 'contemp', diag%axesTL, Time, &
          'Conservative Temperature', 'Celsius')
      CS%id_Sabs = register_diag_field('ocean_model', 'abssalt', diag%axesTL, Time, &
          long_name='Absolute Salinity', units='g kg-1')
      CS%id_sstcon = register_diag_field('ocean_model', 'conSST', diag%axesT1, Time,     &
          'Sea Surface Conservative Temperature', 'Celsius', CS%missing)
      CS%id_sssabs = register_diag_field('ocean_model', 'absSSS', diag%axesT1, Time,     &
          'Sea Surface Absolute Salinity', 'g kg-1', CS%missing)
    endif
  endif

  if (CS%use_temperature .and. CS%use_frazil) then
    CS%id_fraz = register_diag_field('ocean_model', 'frazil', diag%axesT1, Time,                         &
          'Heat from frazil formation', 'W m-2', cmor_field_name='hfsifrazil',                    &
          cmor_standard_name='heat_flux_into_sea_water_due_to_frazil_ice_formation', &
          cmor_long_name='Heat Flux into Sea Water due to Frazil Ice Formation')
  endif

  CS%id_salt_deficit = register_diag_field('ocean_model', 'salt_deficit', diag%axesT1, Time, &
         'Salt sink in ocean due to ice flux', 'psu m-2 s-1')
  CS%id_Heat_PmE = register_diag_field('ocean_model', 'Heat_PmE', diag%axesT1, Time, &
         'Heat flux into ocean from mass flux into ocean', 'W m-2')
  CS%id_intern_heat = register_diag_field('ocean_model', 'internal_heat', diag%axesT1, Time,&
         'Heat flux into ocean from geothermal or other internal sources', 'W m-2')


  ! lateral heat advective and diffusive fluxes
  CS%id_Tadx = register_diag_field('ocean_model', 'T_adx', diag%axesCuL, Time,          &
      'Advective (by residual mean) Zonal Flux of Potential Temperature', 'W m-2', &
      v_extensive = .true., conversion = conv2watt)
  CS%id_Tady = register_diag_field('ocean_model', 'T_ady', diag%axesCvL, Time,               &
      'Advective (by residual mean) Meridional Flux of Potential Temperature', 'W m-2', &
      v_extensive = .true., conversion = conv2watt)
  CS%id_Tdiffx = register_diag_field('ocean_model', 'T_diffx', diag%axesCuL, Time,        &
      'Diffusive Zonal Flux of Potential Temperature', 'W m-2',                      &
      v_extensive = .true., conversion = conv2watt)
  CS%id_Tdiffy = register_diag_field('ocean_model', 'T_diffy', diag%axesCvL, Time, &
      'Diffusive Meridional Flux of Potential Temperature', 'W m-2',          &
      v_extensive = .true., conversion = conv2watt)
  if (CS%id_Tadx   > 0) call safe_alloc_ptr(CS%T_adx,IsdB,IedB,jsd,jed,nz)
  if (CS%id_Tady   > 0) call safe_alloc_ptr(CS%T_ady,isd,ied,JsdB,JedB,nz)
  if (CS%id_Tdiffx > 0) call safe_alloc_ptr(CS%T_diffx,IsdB,IedB,jsd,jed,nz)
  if (CS%id_Tdiffy > 0) call safe_alloc_ptr(CS%T_diffy,isd,ied,JsdB,JedB,nz)


  ! lateral salt advective and diffusive fluxes
  CS%id_Sadx = register_diag_field('ocean_model', 'S_adx', diag%axesCuL, Time, &
      'Advective (by residual mean) Zonal Flux of Salinity', S_flux_units, v_extensive = .true.)
  CS%id_Sady = register_diag_field('ocean_model', 'S_ady', diag%axesCvL, Time, &
      'Advective (by residual mean) Meridional Flux of Salinity', S_flux_units, v_extensive = .true.)
  CS%id_Sdiffx = register_diag_field('ocean_model', 'S_diffx', diag%axesCuL, Time, &
      'Diffusive Zonal Flux of Salinity', S_flux_units, v_extensive = .true.)
  CS%id_Sdiffy = register_diag_field('ocean_model', 'S_diffy', diag%axesCvL, Time, &
      'Diffusive Meridional Flux of Salinity', S_flux_units, v_extensive = .true.)
  if (CS%id_Sadx   > 0) call safe_alloc_ptr(CS%S_adx,IsdB,IedB,jsd,jed,nz)
  if (CS%id_Sady   > 0) call safe_alloc_ptr(CS%S_ady,isd,ied,JsdB,JedB,nz)
  if (CS%id_Sdiffx > 0) call safe_alloc_ptr(CS%S_diffx,IsdB,IedB,jsd,jed,nz)
  if (CS%id_Sdiffy > 0) call safe_alloc_ptr(CS%S_diffy,isd,ied,JsdB,JedB,nz)


  ! vertically integrated lateral heat advective and diffusive fluxes
  CS%id_Tadx_2d = register_diag_field('ocean_model', 'T_adx_2d', diag%axesCu1, Time, &
      'Vertically Integrated Advective Zonal Flux of Potential Temperature', 'W m-2', conversion = conv2watt)
  CS%id_Tady_2d = register_diag_field('ocean_model', 'T_ady_2d', diag%axesCv1, Time, &
      'Vertically Integrated Advective Meridional Flux of Potential Temperature', 'W m-2', conversion = conv2watt)
  CS%id_Tdiffx_2d = register_diag_field('ocean_model', 'T_diffx_2d', diag%axesCu1, Time, &
      'Vertically Integrated Diffusive Zonal Flux of Potential Temperature', 'W m-2', conversion = conv2watt)
  CS%id_Tdiffy_2d = register_diag_field('ocean_model', 'T_diffy_2d', diag%axesCv1, Time, &
      'Vertically Integrated Diffusive Meridional Flux of Potential Temperature', 'W m-2', conversion = conv2watt)
  if (CS%id_Tadx_2d   > 0) call safe_alloc_ptr(CS%T_adx_2d,IsdB,IedB,jsd,jed)
  if (CS%id_Tady_2d   > 0) call safe_alloc_ptr(CS%T_ady_2d,isd,ied,JsdB,JedB)
  if (CS%id_Tdiffx_2d > 0) call safe_alloc_ptr(CS%T_diffx_2d,IsdB,IedB,jsd,jed)
  if (CS%id_Tdiffy_2d > 0) call safe_alloc_ptr(CS%T_diffy_2d,isd,ied,JsdB,JedB)

  ! vertically integrated lateral salt advective and diffusive fluxes
  CS%id_Sadx_2d = register_diag_field('ocean_model', 'S_adx_2d', diag%axesCu1, Time, &
      'Vertically Integrated Advective Zonal Flux of Salinity', S_flux_units)
  CS%id_Sady_2d = register_diag_field('ocean_model', 'S_ady_2d', diag%axesCv1, Time, &
      'Vertically Integrated Advective Meridional Flux of Salinity', S_flux_units)
  CS%id_Sdiffx_2d = register_diag_field('ocean_model', 'S_diffx_2d', diag%axesCu1, Time, &
      'Vertically Integrated Diffusive Zonal Flux of Salinity', S_flux_units)
  CS%id_Sdiffy_2d = register_diag_field('ocean_model', 'S_diffy_2d', diag%axesCv1, Time, &
      'Vertically Integrated Diffusive Meridional Flux of Salinity', S_flux_units)
  if (CS%id_Sadx_2d   > 0) call safe_alloc_ptr(CS%S_adx_2d,IsdB,IedB,jsd,jed)
  if (CS%id_Sady_2d   > 0) call safe_alloc_ptr(CS%S_ady_2d,isd,ied,JsdB,JedB)
  if (CS%id_Sdiffx_2d > 0) call safe_alloc_ptr(CS%S_diffx_2d,IsdB,IedB,jsd,jed)
  if (CS%id_Sdiffy_2d > 0) call safe_alloc_ptr(CS%S_diffy_2d,isd,ied,JsdB,JedB)


  if (CS%debug_truncations) then
    call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
    if (.not.CS%adiabatic) then
      call safe_alloc_ptr(ADp%du_dt_dia,IsdB,IedB,jsd,jed,nz)
      call safe_alloc_ptr(ADp%dv_dt_dia,isd,ied,JsdB,JedB,nz)
    endif
  endif

  ! diagnostics for values prior to diabatic and prior to ALE
  CS%id_u_predia = register_diag_field('ocean_model', 'u_predia', diag%axesCuL, Time, &
      'Zonal velocity before diabatic forcing', 'm s-1')
  CS%id_v_predia = register_diag_field('ocean_model', 'v_predia', diag%axesCvL, Time, &
      'Meridional velocity before diabatic forcing', 'm s-1')
  CS%id_h_predia = register_diag_field('ocean_model', 'h_predia', diag%axesTL, Time, &
      'Layer Thickness before diabatic forcing', thickness_units, v_extensive=.true.)
  CS%id_e_predia = register_diag_field('ocean_model', 'e_predia', diag%axesTi, Time, &
      'Interface Heights before diabatic forcing', 'm')
  if (.not. CS%adiabatic) then
    CS%id_u_preale = register_diag_field('ocean_model', 'u_preale', diag%axesCuL, Time, &
        'Zonal velocity before remapping', 'm s-1')
    CS%id_v_preale = register_diag_field('ocean_model', 'v_preale', diag%axesCvL, Time, &
        'Meridional velocity before remapping', 'm s-1')
    CS%id_h_preale = register_diag_field('ocean_model', 'h_preale', diag%axesTL, Time, &
        'Layer Thickness before remapping', thickness_units, v_extensive=.true.)
    CS%id_T_preale = register_diag_field('ocean_model', 'T_preale', diag%axesTL, Time, &
        'Temperature before remapping', 'degC')
    CS%id_S_preale = register_diag_field('ocean_model', 'S_preale', diag%axesTL, Time, &
        'Salinity before remapping', 'psu')
    CS%id_e_preale = register_diag_field('ocean_model', 'e_preale', diag%axesTi, Time, &
        'Interface Heights before remapping', 'm')
  endif

  if (CS%use_temperature) then
    CS%id_T_predia = register_diag_field('ocean_model', 'temp_predia', diag%axesTL, Time, &
        'Potential Temperature', 'degC')
    CS%id_S_predia = register_diag_field('ocean_model', 'salt_predia', diag%axesTL, Time, &
        'Salinity', 'psu')
  endif

  ! Diagnostics related to tracer transport
  CS%id_uhtr = register_diag_field('ocean_model', 'uhtr', diag%axesCuL, Time, &
      'Accumulated zonal thickness fluxes to advect tracers', 'kg', &
      y_cell_method='sum', v_extensive=.true.)
  CS%id_vhtr = register_diag_field('ocean_model', 'vhtr', diag%axesCvL, Time, &
      'Accumulated meridional thickness fluxes to advect tracers', 'kg', &
      x_cell_method='sum', v_extensive=.true.)
  CS%id_umo = register_diag_field('ocean_model', 'umo', &
      diag%axesCuL, Time, 'Ocean Mass X Transport', 'kg s-1', &
      standard_name='ocean_mass_x_transport', y_cell_method='sum', v_extensive=.true.)
  CS%id_vmo = register_diag_field('ocean_model', 'vmo', &
      diag%axesCvL, Time, 'Ocean Mass Y Transport', 'kg s-1', &
      standard_name='ocean_mass_y_transport', x_cell_method='sum', v_extensive=.true.)
  CS%id_umo_2d = register_diag_field('ocean_model', 'umo_2d', &
      diag%axesCu1, Time, 'Ocean Mass X Transport Vertical Sum', 'kg s-1', &
      standard_name='ocean_mass_x_transport_vertical_sum', y_cell_method='sum')
  CS%id_vmo_2d = register_diag_field('ocean_model', 'vmo_2d', &
      diag%axesCv1, Time, 'Ocean Mass Y Transport Vertical Sum', 'kg s-1', &
      standard_name='ocean_mass_y_transport_vertical_sum', x_cell_method='sum')

end subroutine register_diags


!> Initialize diagnostics for temp/heat and saln/salt tendencies.
subroutine register_diags_TS_tendency(Time, G, CS)
  type(time_type),           intent(in)    :: Time  !< current model time
  type(ocean_grid_type),     intent(inout) :: G     !< ocean grid structure
  type(MOM_control_struct),  pointer       :: CS    !< control structure set up by initialize_MOM

  type(diag_ctrl), pointer :: diag
  integer :: i, j, k
  integer :: isd, ied, jsd, jed, nz
  integer :: is, ie, js, je

  diag => CS%diag
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed ; nz = G%ke
  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec


  ! heat tendencies from lateral advection
  CS%id_T_advection_xy = register_diag_field('ocean_model', 'T_advection_xy', diag%axesTL, Time, &
      'Horizontal convergence of residual mean heat advective fluxes', 'W m-2',v_extensive=.true.)
  CS%id_T_advection_xy_2d = register_diag_field('ocean_model', 'T_advection_xy_2d', diag%axesT1, Time,&
      'Vertical sum of horizontal convergence of residual mean heat advective fluxes', 'W m-2')
  if (CS%id_T_advection_xy > 0 .or. CS%id_T_advection_xy_2d > 0) then
    call safe_alloc_ptr(CS%T_advection_xy,isd,ied,jsd,jed,nz)
    CS%tendency_diagnostics = .true.
  endif

  ! net temperature and heat tendencies
  CS%id_T_tendency = register_diag_field('ocean_model', 'T_tendency', diag%axesTL, Time, &
      'Net time tendency for temperature', 'degC s-1')
  CS%id_Th_tendency = register_diag_field('ocean_model', 'Th_tendency', diag%axesTL, Time,        &
      'Net time tendency for heat', 'W m-2',                                                      &
      cmor_field_name="opottemptend",                                                             &
      cmor_standard_name="tendency_of_sea_water_potential_temperature_expressed_as_heat_content", &
      cmor_long_name ="Tendency of Sea Water Potential Temperature Expressed as Heat Content",    &
      v_extensive=.true.)
  CS%id_Th_tendency_2d = register_diag_field('ocean_model', 'Th_tendency_2d', diag%axesT1, Time,              &
      'Vertical sum of net time tendency for heat', 'W m-2',                                                  &
      cmor_field_name="opottemptend_2d",                                                                      &
      cmor_standard_name="tendency_of_sea_water_potential_temperature_expressed_as_heat_content_vertical_sum",&
      cmor_long_name ="Tendency of Sea Water Potential Temperature Expressed as Heat Content Vertical Sum")
  if (CS%id_T_tendency > 0) then
    CS%tendency_diagnostics = .true.
    call safe_alloc_ptr(CS%T_prev,isd,ied,jsd,jed,nz)
    do k=1,nz ; do j=js,je ; do i=is,ie
      CS%T_prev(i,j,k) = CS%tv%T(i,j,k)
    enddo ; enddo ; enddo
  endif
  if (CS%id_Th_tendency > 0 .or. CS%id_Th_tendency_2d > 0) then
    CS%tendency_diagnostics = .true.
    call safe_alloc_ptr(CS%Th_prev,isd,ied,jsd,jed,nz)
    do k=1,nz ; do j=js,je ; do i=is,ie
      CS%Th_prev(i,j,k) = CS%tv%T(i,j,k) * CS%h(i,j,k)
    enddo ; enddo ; enddo
  endif


  ! salt tendencies from lateral advection
  CS%id_S_advection_xy = register_diag_field('ocean_model', 'S_advection_xy', diag%axesTL, Time, &
      'Horizontal convergence of residual mean salt advective fluxes', 'kg m-2 s-1', v_extensive=.true.)
  CS%id_S_advection_xy_2d = register_diag_field('ocean_model', 'S_advection_xy_2d', diag%axesT1, Time,&
      'Vertical sum of horizontal convergence of residual mean salt advective fluxes', 'kg m-2 s-1')
  if (CS%id_S_advection_xy > 0 .or. CS%id_S_advection_xy_2d > 0) then
    call safe_alloc_ptr(CS%S_advection_xy,isd,ied,jsd,jed,nz)
    CS%tendency_diagnostics = .true.
  endif

  ! net salinity and salt tendencies
  CS%id_S_tendency = register_diag_field('ocean_model', 'S_tendency', diag%axesTL, Time, &
      'Net time tendency for salinity', 'psu s-1')
  CS%id_Sh_tendency = register_diag_field('ocean_model', 'Sh_tendency', diag%axesTL, Time,&
      'Net time tendency for salt', 'kg m-2 s-1',                                         &
      cmor_field_name="osalttend",                                                        &
      cmor_standard_name="tendency_of_sea_water_salinity_expressed_as_salt_content",      &
      cmor_long_name ="Tendency of Sea Water Salinity Expressed as Salt Content",         &
      v_extensive=.true.)
  CS%id_Sh_tendency_2d = register_diag_field('ocean_model', 'Sh_tendency_2d', diag%axesT1, Time, &
      'Vertical sum of net time tendency for salt', 'kg m-2 s-1',                                &
      cmor_field_name="osalttend_2d",                                                            &
      cmor_standard_name="tendency_of_sea_water_salinity_expressed_as_salt_content_vertical_sum",&
      cmor_long_name ="Tendency of Sea Water Salinity Expressed as Salt Content Vertical Sum")
  if (CS%id_S_tendency > 0) then
    CS%tendency_diagnostics = .true.
    call safe_alloc_ptr(CS%S_prev,isd,ied,jsd,jed,nz)
    do k=1,nz ; do j=js,je ; do i=is,ie
      CS%S_prev(i,j,k) = CS%tv%S(i,j,k)
    enddo ; enddo ; enddo
  endif
  if (CS%id_Sh_tendency > 0 .or. CS%id_Sh_tendency_2d > 0) then
    CS%tendency_diagnostics = .true.
    call safe_alloc_ptr(CS%Sh_prev,isd,ied,jsd,jed,nz)
    do k=1,nz ; do j=js,je ; do i=is,ie
      CS%Sh_prev(i,j,k) = CS%tv%S(i,j,k) * CS%h(i,j,k)
    enddo ; enddo ; enddo
  endif

end subroutine register_diags_TS_tendency


!> Initialize diagnostics for the variance decay of temp/salt
!! across regridding/remapping
subroutine register_diags_TS_vardec(Time, HI, GV, param_file, CS)
  type(time_type),         intent(in) :: Time     !< current model time
  type(hor_index_type),    intent(in) :: HI       !< horizontal index type
  type(verticalGrid_type), intent(in) :: GV       !< ocean vertical grid structure
  type(param_file_type),   intent(in) :: param_file !< parameter file
  type(MOM_control_struct), pointer :: CS   !< control structure for MOM

  integer :: isd, ied, jsd, jed, nz
  type(vardesc) :: vd_tmp
  type(diag_ctrl), pointer :: diag

  diag => CS%diag
  isd  = HI%isd  ; ied  = HI%ied  ; jsd  = HI%jsd  ; jed  = HI%jed ; nz = GV%ke

  ! variancy decay through ALE operation
  CS%id_T_vardec = register_diag_field('ocean_model', 'T_vardec', diag%axesTL, Time, &
      'ALE variance decay for temperature', 'degC2 s-1')
  if (CS%id_T_vardec > 0) then
    call safe_alloc_ptr(CS%T_squared,isd,ied,jsd,jed,nz)
    CS%T_squared(:,:,:) = 0.

    vd_tmp = var_desc(name="T2", units="degC2", longname="Squared Potential Temperature")
    call register_tracer(CS%T_squared, vd_tmp, param_file, HI, GV, CS%tracer_reg)
  endif

  CS%id_S_vardec = register_diag_field('ocean_model', 'S_vardec', diag%axesTL, Time, &
      'ALE variance decay for salinity', 'psu2 s-1')
  if (CS%id_S_vardec > 0) then
    call safe_alloc_ptr(CS%S_squared,isd,ied,jsd,jed,nz)
    CS%S_squared(:,:,:) = 0.

    vd_tmp = var_desc(name="S2", units="psu2", longname="Squared Salinity")
    call register_tracer(CS%S_squared, vd_tmp, param_file, HI, GV, CS%tracer_reg)
  endif

end subroutine register_diags_TS_vardec

!> This subroutine sets up clock IDs for timing various subroutines.
subroutine MOM_timing_init(CS)
  type(MOM_control_struct), intent(in) :: CS  !< control structure set up by initialize_MOM.

 id_clock_ocean    = cpu_clock_id('Ocean', grain=CLOCK_COMPONENT)
 id_clock_dynamics = cpu_clock_id('Ocean dynamics', grain=CLOCK_SUBCOMPONENT)
 id_clock_thermo   = cpu_clock_id('Ocean thermodynamics and tracers', grain=CLOCK_SUBCOMPONENT)
 id_clock_other    = cpu_clock_id('Ocean Other', grain=CLOCK_SUBCOMPONENT)
 id_clock_tracer   = cpu_clock_id('(Ocean tracer advection)', grain=CLOCK_MODULE_DRIVER)
 if (.not.CS%adiabatic) &
   id_clock_diabatic = cpu_clock_id('(Ocean diabatic driver)', grain=CLOCK_MODULE_DRIVER)

 id_clock_continuity = cpu_clock_id('(Ocean continuity equation *)', grain=CLOCK_MODULE)
 id_clock_BBL_visc = cpu_clock_id('(Ocean set BBL viscosity)', grain=CLOCK_MODULE)
 id_clock_pass       = cpu_clock_id('(Ocean message passing *)', grain=CLOCK_MODULE)
 id_clock_MOM_init   = cpu_clock_id('(Ocean MOM_initialize_state)', grain=CLOCK_MODULE)
 id_clock_pass_init  = cpu_clock_id('(Ocean init message passing *)', grain=CLOCK_ROUTINE)
 if (CS%thickness_diffuse) &
   id_clock_thick_diff = cpu_clock_id('(Ocean thickness diffusion *)', grain=CLOCK_MODULE)
!if (CS%mixedlayer_restrat) &
   id_clock_ml_restrat = cpu_clock_id('(Ocean mixed layer restrat)', grain=CLOCK_MODULE)
 id_clock_diagnostics  = cpu_clock_id('(Ocean collective diagnostics)', grain=CLOCK_MODULE)
 id_clock_Z_diag       = cpu_clock_id('(Ocean Z-space diagnostics)', grain=CLOCK_MODULE)
 id_clock_ALE          = cpu_clock_id('(Ocean ALE)', grain=CLOCK_MODULE)
 if(CS%offline_tracer_mode) then
  id_clock_offline_tracer = cpu_clock_id('Ocean offline tracers', grain=CLOCK_SUBCOMPONENT)
 endif

end subroutine MOM_timing_init

!> This routine posts diagnostics of the transports, including the subgridscale
!! contributions.
subroutine post_transport_diagnostics(G, GV, CS, diag, dt_trans, h, h_pre_dyn)
  type(ocean_grid_type),    intent(inout) :: G   !< ocean grid structure
  type(verticalGrid_type),  intent(in)    :: GV  !< ocean vertical grid structure
  type(MOM_control_struct), intent(in)    :: CS  !< control structure
  type(diag_ctrl),          intent(inout) :: diag !< regulates diagnostic output
  real                    , intent(in)    :: dt_trans !< total time step associated with the transports, in s.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                            intent(in)    :: h   !< The updated layer thicknesses, in H
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                            intent(in)    :: h_pre_dyn !< The thickness before the transports, in H.

  real, dimension(SZIB_(G), SZJ_(G)) :: umo2d ! Diagnostics of integrated mass transport, in kg s-1
  real, dimension(SZI_(G), SZJB_(G)) :: vmo2d ! Diagnostics of integrated mass transport, in kg s-1
  real, dimension(SZIB_(G), SZJ_(G), SZK_(G)) :: umo ! Diagnostics of layer mass transport, in kg s-1
  real, dimension(SZI_(G), SZJB_(G), SZK_(G)) :: vmo ! Diagnostics of layer mass transport, in kg s-1
  real :: H_to_kg_m2_dt   ! A conversion factor from accumulated transports to fluxes, in kg m-2 H-1 s-1.
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call cpu_clock_begin(id_clock_Z_diag)
  call calculate_Z_transport(CS%uhtr, CS%vhtr, h, dt_trans, G, GV, &
                             CS%diag_to_Z_CSp)
  call cpu_clock_end(id_clock_Z_diag)

  ! Post mass transports, including SGS
  ! Build the remap grids using the layer thicknesses from before the dynamics
  call diag_update_remap_grids(diag, alt_h = h_pre_dyn)

  H_to_kg_m2_dt = GV%H_to_kg_m2 / dt_trans
  if (CS%id_umo_2d > 0) then
    umo2d(:,:) = 0.0
    do k=1,nz ; do j=js,je ; do I=is-1,ie
      umo2d(I,j) = umo2d(I,j) + CS%uhtr(I,j,k) * H_to_kg_m2_dt
    enddo ; enddo ; enddo
    call post_data(CS%id_umo_2d, umo2d, diag)
  endif
  if (CS%id_umo > 0) then
    ! Convert to kg/s. Modifying the array for diagnostics is allowed here since it is set to zero immediately below
    do k=1,nz ; do j=js,je ; do I=is-1,ie
      umo(I,j,k) =  CS%uhtr(I,j,k) * H_to_kg_m2_dt
    enddo ; enddo ; enddo
    call post_data(CS%id_umo, umo, diag, alt_h = h_pre_dyn)
  endif
  if (CS%id_vmo_2d > 0) then
    vmo2d(:,:) = 0.0
    do k=1,nz ; do J=js-1,je ; do i=is,ie
      vmo2d(i,J) = vmo2d(i,J) + CS%vhtr(i,J,k) * H_to_kg_m2_dt
    enddo ; enddo ; enddo
    call post_data(CS%id_vmo_2d, vmo2d, diag)
  endif
  if (CS%id_vmo > 0) then
    ! Convert to kg/s. Modifying the array for diagnostics is allowed here since it is set to zero immediately below
    do k=1,nz ; do J=js-1,je ; do i=is,ie
      vmo(i,J,k) = CS%vhtr(i,J,k) * H_to_kg_m2_dt
    enddo ; enddo ; enddo
    call post_data(CS%id_vmo, vmo, diag, alt_h = h_pre_dyn)
  endif

  if (CS%id_uhtr > 0) call post_data(CS%id_uhtr, CS%uhtr, diag, alt_h = h_pre_dyn)
  if (CS%id_vhtr > 0) call post_data(CS%id_vhtr, CS%vhtr, diag, alt_h = h_pre_dyn)

end subroutine post_transport_diagnostics

!> Post diagnostics of temperatures and salinities, their fluxes, and tendencies.
subroutine post_TS_diagnostics(CS, G, GV, tv, diag, dt)
  type(MOM_control_struct), intent(inout) :: CS  !< control structure
  type(ocean_grid_type),    intent(in)    :: G   !< ocean grid structure
  type(verticalGrid_type),  intent(in)    :: GV  !< ocean vertical grid structure
  type(thermo_var_ptrs),    intent(in)    :: tv  !< A structure pointing to various thermodynamic variables
  type(diag_ctrl),          intent(in)    :: diag !< regulates diagnostic output
  real,                     intent(in)    :: dt  !< total time step for T,S update

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: potTemp, pracSal !TEOS10 Diagnostics
  real    :: work3d(SZI_(G),SZJ_(G),SZK_(G))
  real    :: work2d(SZI_(G),SZJ_(G))
  real    :: Idt, ppt2mks
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.NOT. CS%use_conT_absS) then
    !Internal T&S variables are assumed to be potential&practical
    if (CS%id_T > 0) call post_data(CS%id_T, tv%T, diag)
    if (CS%id_S > 0) call post_data(CS%id_S, tv%S, diag)

    if (CS%id_tob > 0) call post_data(CS%id_tob, tv%T(:,:,G%ke), diag, mask=G%mask2dT)
    if (CS%id_sob > 0) call post_data(CS%id_sob, tv%S(:,:,G%ke), diag, mask=G%mask2dT)
  else
    !Internal T&S variables are assumed to be conservative&absolute
    if (CS%id_Tcon > 0) call post_data(CS%id_Tcon, tv%T, diag)
    if (CS%id_Sabs > 0) call post_data(CS%id_Sabs, tv%S, diag)
    !Using TEOS-10 function calls convert T&S diagnostics
    !from conservative temp to potential temp and
    !from absolute salinity to practical salinity
    do k=1,nz ; do j=js,je ; do i=is,ie
      pracSal(i,j,k) = gsw_sp_from_sr(tv%S(i,j,k))
      potTemp(i,j,k) = gsw_pt_from_ct(tv%S(i,j,k),tv%T(i,j,k))
    enddo; enddo ; enddo
    if (CS%id_T > 0) call post_data(CS%id_T, potTemp, diag)
    if (CS%id_S > 0) call post_data(CS%id_S, pracSal, diag)
    if (CS%id_tob > 0) call post_data(CS%id_tob, potTemp(:,:,G%ke), diag, mask=G%mask2dT)
    if (CS%id_sob > 0) call post_data(CS%id_sob, pracSal(:,:,G%ke), diag, mask=G%mask2dT)
  endif

  if (CS%id_Tadx   > 0) call post_data(CS%id_Tadx,   CS%T_adx,   diag)
  if (CS%id_Tady   > 0) call post_data(CS%id_Tady,   CS%T_ady,   diag)
  if (CS%id_Tdiffx > 0) call post_data(CS%id_Tdiffx, CS%T_diffx, diag)
  if (CS%id_Tdiffy > 0) call post_data(CS%id_Tdiffy, CS%T_diffy, diag)

  if (CS%id_Sadx   > 0) call post_data(CS%id_Sadx,   CS%S_adx,   diag)
  if (CS%id_Sady   > 0) call post_data(CS%id_Sady,   CS%S_ady,   diag)
  if (CS%id_Sdiffx > 0) call post_data(CS%id_Sdiffx, CS%S_diffx, diag)
  if (CS%id_Sdiffy > 0) call post_data(CS%id_Sdiffy, CS%S_diffy, diag)

  if (CS%id_Tadx_2d   > 0) call post_data(CS%id_Tadx_2d,   CS%T_adx_2d,   diag)
  if (CS%id_Tady_2d   > 0) call post_data(CS%id_Tady_2d,   CS%T_ady_2d,   diag)
  if (CS%id_Tdiffx_2d > 0) call post_data(CS%id_Tdiffx_2d, CS%T_diffx_2d, diag)
  if (CS%id_Tdiffy_2d > 0) call post_data(CS%id_Tdiffy_2d, CS%T_diffy_2d, diag)

  if (CS%id_Sadx_2d   > 0) call post_data(CS%id_Sadx_2d,   CS%S_adx_2d,   diag)
  if (CS%id_Sady_2d   > 0) call post_data(CS%id_Sady_2d,   CS%S_ady_2d,   diag)
  if (CS%id_Sdiffx_2d > 0) call post_data(CS%id_Sdiffx_2d, CS%S_diffx_2d, diag)
  if (CS%id_Sdiffy_2d > 0) call post_data(CS%id_Sdiffy_2d, CS%S_diffy_2d, diag)

  if(.not. CS%tendency_diagnostics) return

  Idt = 0.; if (dt/=0.) Idt = 1.0 / dt ! The "if" is in case the diagnostic is called for a zero length interval
  ppt2mks       = 0.001
  work3d(:,:,:) = 0.0
  work2d(:,:)   = 0.0

  ! Diagnose tendency of heat from convergence of lateral advective,
  ! fluxes, where advective transport arises from residual mean velocity.
  if (CS%id_T_advection_xy > 0 .or. CS%id_T_advection_xy_2d > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work3d(i,j,k) = CS%T_advection_xy(i,j,k) * GV%H_to_kg_m2 * tv%C_p
    enddo ; enddo ; enddo
    if (CS%id_T_advection_xy    > 0) call post_data(CS%id_T_advection_xy, work3d, diag)
    if (CS%id_T_advection_xy_2d > 0) then
      do j=js,je ; do i=is,ie
        work2d(i,j) = 0.0
        do k=1,nz
          work2d(i,j) = work2d(i,j) + work3d(i,j,k)
        enddo
      enddo ; enddo
      call post_data(CS%id_T_advection_xy_2d, work2d, diag)
    endif
  endif

  ! Diagnose tendency of salt from convergence of lateral advective
  ! fluxes, where advective transport arises from residual mean velocity.
  if (CS%id_S_advection_xy > 0 .or. CS%id_S_advection_xy_2d > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work3d(i,j,k) = CS%S_advection_xy(i,j,k) * GV%H_to_kg_m2 * ppt2mks
    enddo ; enddo ; enddo
    if (CS%id_S_advection_xy    > 0) call post_data(CS%id_S_advection_xy, work3d, diag)
    if (CS%id_S_advection_xy_2d > 0) then
      do j=js,je ; do i=is,ie
        work2d(i,j) = 0.0
        do k=1,nz
          work2d(i,j) = work2d(i,j) + work3d(i,j,k)
        enddo
      enddo ; enddo
      call post_data(CS%id_S_advection_xy_2d, work2d, diag)
    endif
  endif

  ! diagnose net tendency for temperature over a time step and update T_prev
  if (CS%id_T_tendency > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work3d(i,j,k)    = (tv%T(i,j,k) - CS%T_prev(i,j,k))*Idt
      CS%T_prev(i,j,k) =  tv%T(i,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_T_tendency, work3d, diag)
  endif

  ! diagnose net tendency for salinity over a time step and update S_prev
  if (CS%id_S_tendency > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work3d(i,j,k)    = (tv%S(i,j,k) - CS%S_prev(i,j,k))*Idt
      CS%S_prev(i,j,k) =  tv%S(i,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_S_tendency, work3d, diag)
  endif

  ! diagnose net tendency for heat content of a grid cell over a time step and update Th_prev
  if (CS%id_Th_tendency > 0 .or. CS%id_Th_tendency_2d > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work3d(i,j,k)     = (tv%T(i,j,k)*CS%h(i,j,k) - CS%Th_prev(i,j,k)) * Idt * GV%H_to_kg_m2 * tv%C_p
      CS%Th_prev(i,j,k) =  tv%T(i,j,k)*CS%h(i,j,k)
    enddo ; enddo ; enddo
    if (CS%id_Th_tendency    > 0) call post_data(CS%id_Th_tendency, work3d, diag)
    if (CS%id_Th_tendency_2d > 0) then
      do j=js,je ; do i=is,ie
        work2d(i,j) = 0.0
        do k=1,nz
          work2d(i,j) = work2d(i,j) + work3d(i,j,k)
        enddo
      enddo ; enddo
      call post_data(CS%id_Th_tendency_2d, work2d, diag)
    endif
  endif

  ! diagnose net tendency for salt content of a grid cell over a time step and update Sh_prev
  if (CS%id_Sh_tendency > 0 .or. CS%id_Sh_tendency_2d > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work3d(i,j,k)     = (tv%S(i,j,k)*CS%h(i,j,k) - CS%Sh_prev(i,j,k)) * Idt * GV%H_to_kg_m2 * ppt2mks
      CS%Sh_prev(i,j,k) =  tv%S(i,j,k)*CS%h(i,j,k)
    enddo ; enddo ; enddo
    if (CS%id_Sh_tendency    > 0) call post_data(CS%id_Sh_tendency, work3d, diag)
    if (CS%id_Sh_tendency_2d > 0) then
      do j=js,je ; do i=is,ie
        work2d(i,j) = 0.0
        do k=1,nz
          work2d(i,j) = work2d(i,j) + work3d(i,j,k)
        enddo
      enddo ; enddo
      call post_data(CS%id_Sh_tendency_2d, work2d, diag)
    endif
  endif

end subroutine post_TS_diagnostics

!> Calculate and post variance decay diagnostics for temp/salt
subroutine post_diags_TS_vardec(G, CS, dt)
  type(ocean_grid_type),    intent(in) :: G    !< ocean grid structure
  type(MOM_control_struct), intent(in) :: CS   !< control structure
  real, intent(in) :: dt                       !< total time step

  real :: work(SZI_(G),SZJ_(G),SZK_(G))
  real :: Idt
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  Idt = 0.; if (dt/=0.) Idt = 1.0 / dt ! The "if" is in case the diagnostic is called for a zero length interval

  if (CS%id_T_vardec > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work(i,j,k) = (CS%T_squared(i,j,k) - CS%tv%T(i,j,k)**2) * Idt
    enddo ; enddo ; enddo
    call post_data(CS%id_T_vardec, work, CS%diag)
  endif

  if (CS%id_S_vardec > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work(i,j,k) = (CS%S_squared(i,j,k) - CS%tv%S(i,j,k)**2) * Idt
    enddo ; enddo ; enddo
    call post_data(CS%id_S_vardec, work, CS%diag)
  endif
end subroutine post_diags_TS_vardec

!> This routine posts diagnostics of various integrated quantities.
subroutine post_integrated_diagnostics(CS, G, GV, diag, dt_int, tv, ssh, fluxes)
  type(MOM_control_struct), intent(in) :: CS  !< control structure
  type(ocean_grid_type),    intent(in) :: G   !< ocean grid structure
  type(verticalGrid_type),  intent(in) :: GV  !< ocean vertical grid structure
  type(diag_ctrl),          intent(in) :: diag  !< regulates diagnostic output
  real,                     intent(in) :: dt_int  !< total time step associated with these diagnostics, in s.
  type(thermo_var_ptrs),    intent(in) :: tv  !< A structure pointing to various thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G)), &
                            intent(in) :: ssh !< Time mean surface height without
                                              !! corrections for ice displacement(m)
  type(forcing),            intent(in) :: fluxes  !< pointers to forcing fields

  real, allocatable, dimension(:,:) :: &
    tmp,              & ! temporary 2d field
    zos,              & ! dynamic sea lev (zero area mean) from inverse-barometer adjusted ssh (meter)
    zossq,            & ! square of zos (m^2)
    sfc_speed,        & ! sea surface speed at h-points (m/s)
    frazil_ave,       & ! average frazil heat flux required to keep temp above freezing (W/m2)
    salt_deficit_ave, & ! average salt flux required to keep salinity above 0.01ppt (gSalt m-2 s-1)
    Heat_PmE_ave,     & ! average effective heat flux into the ocean due to
                        ! the exchange of water with other components, times the
                        ! heat capacity of water, in W m-2.
    intern_heat_ave     ! avg heat flux into ocean from geothermal or
                        ! other internal heat sources (W/m2)
  real :: I_time_int    ! The inverse of the time interval in s-1.
  real :: zos_area_mean, volo, ssh_ga
  integer :: i, j, k, is, ie, js, je, nz! , Isq, Ieq, Jsq, Jeq

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  ! area mean SSH
  if (CS%id_ssh_ga > 0) then
    ssh_ga = global_area_mean(ssh, G)
    call post_data(CS%id_ssh_ga, ssh_ga, diag)
  endif

  I_time_int = 1.0 / dt_int
  if (CS%id_ssh > 0) &
    call post_data(CS%id_ssh, ssh, diag, mask=G%mask2dT)

  ! post the dynamic sea level, zos, and zossq.
  ! zos is ave_ssh with sea ice inverse barometer removed,
  ! and with zero global area mean.
  if(CS%id_zos > 0 .or. CS%id_zossq > 0) then
     allocate(zos(G%isd:G%ied,G%jsd:G%jed))
     zos(:,:) = 0.0
     do j=js,je ; do i=is,ie
       zos(i,j) = ssh(i,j)
     enddo ; enddo
     if (ASSOCIATED(fluxes%p_surf)) then
       do j=js,je ; do i=is,ie
         zos(i,j) = zos(i,j) + G%mask2dT(i,j)*fluxes%p_surf(i,j) / &
                              (GV%Rho0 * GV%g_Earth)
       enddo ; enddo
     endif
     zos_area_mean = global_area_mean(zos, G)
     do j=js,je ; do i=is,ie
       zos(i,j) = zos(i,j) - G%mask2dT(i,j)*zos_area_mean
     enddo ; enddo
     if(CS%id_zos > 0) then
       call post_data(CS%id_zos, zos, diag, mask=G%mask2dT)
     endif
     if(CS%id_zossq > 0) then
       allocate(zossq(G%isd:G%ied,G%jsd:G%jed))
       zossq(:,:) = 0.0
       do j=js,je ; do i=is,ie
         zossq(i,j) = zos(i,j)*zos(i,j)
       enddo ; enddo
       call post_data(CS%id_zossq, zossq, diag, mask=G%mask2dT)
       deallocate(zossq)
     endif
     deallocate(zos)
  endif

  ! post total volume of the liquid ocean
  if(CS%id_volo > 0) then
    allocate(tmp(G%isd:G%ied,G%jsd:G%jed))
    do j=js,je ; do i=is,ie
      tmp(i,j) = G%mask2dT(i,j)*(ssh(i,j) + G%bathyT(i,j))
    enddo ; enddo
    volo = global_area_integral(tmp, G)
    call post_data(CS%id_volo, volo, diag)
    deallocate(tmp)
  endif

  ! post frazil
  if (ASSOCIATED(tv%frazil) .and. (CS%id_fraz > 0)) then
    allocate(frazil_ave(G%isd:G%ied,G%jsd:G%jed))
    do j=js,je ; do i=is,ie
      frazil_ave(i,j) = tv%frazil(i,j) * I_time_int
    enddo ; enddo
    call post_data(CS%id_fraz, frazil_ave, diag, mask=G%mask2dT)
    deallocate(frazil_ave)
  endif

  ! post the salt deficit
  if (ASSOCIATED(tv%salt_deficit) .and. (CS%id_salt_deficit > 0)) then
    allocate(salt_deficit_ave(G%isd:G%ied,G%jsd:G%jed))
    do j=js,je ; do i=is,ie
      salt_deficit_ave(i,j) = tv%salt_deficit(i,j) * I_time_int
    enddo ; enddo
    call post_data(CS%id_salt_deficit, salt_deficit_ave, diag, mask=G%mask2dT)
    deallocate(salt_deficit_ave)
  endif

  ! post temperature of P-E+R
  if (ASSOCIATED(tv%TempxPmE) .and. (CS%id_Heat_PmE > 0)) then
    allocate(Heat_PmE_ave(G%isd:G%ied,G%jsd:G%jed))
    do j=js,je ; do i=is,ie
      Heat_PmE_ave(i,j) = tv%TempxPmE(i,j) * (tv%C_p * I_time_int)
    enddo ; enddo
    call post_data(CS%id_Heat_PmE, Heat_PmE_ave, diag, mask=G%mask2dT)
    deallocate(Heat_PmE_ave)
  endif

  ! post geothermal heating or internal heat source/sinks
  if (ASSOCIATED(tv%internal_heat) .and. (CS%id_intern_heat > 0)) then
    allocate(intern_heat_ave(G%isd:G%ied,G%jsd:G%jed))
    do j=js,je ; do i=is,ie
      intern_heat_ave(i,j) = tv%internal_heat(i,j) * (tv%C_p * I_time_int)
    enddo ; enddo
    call post_data(CS%id_intern_heat, intern_heat_ave, diag, mask=G%mask2dT)
    deallocate(intern_heat_ave)
  endif

end subroutine post_integrated_diagnostics

!> This routine posts diagnostics of various ocean surface quantities.
subroutine post_surface_diagnostics(CS, G, diag, sfc_state)
  type(MOM_control_struct), intent(in)    :: CS  !< control structure
  type(ocean_grid_type),    intent(in)    :: G   !< ocean grid structure
  type(diag_ctrl),          intent(in)    :: diag  !< regulates diagnostic output
  type(surface),            intent(in)    :: sfc_state !< ocean surface state

  real, dimension(SZI_(G),SZJ_(G)) :: &
    potTemp, &  ! TEOS10 potential temperature (deg C)
    pracSal, &  ! TEOS10 practical salinity
    SST_sq, &   ! Surface temperature squared, in degC^2
    SSS_sq, &   ! Surface salinity squared, in salnity units^2
    sfc_speed   ! sea surface speed at h-points (m/s)

  integer :: i, j, is, ie, js, je

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (.NOT.CS%use_conT_absS) then
    !Internal T&S variables are assumed to be potential&practical
    if (CS%id_sst > 0) call post_data(CS%id_sst, sfc_state%SST, diag, mask=G%mask2dT)
    if (CS%id_sss > 0) call post_data(CS%id_sss, sfc_state%SSS, diag, mask=G%mask2dT)
  else
    !Internal T&S variables are assumed to be conservative&absolute
    if (CS%id_sstcon > 0) call post_data(CS%id_sstcon, sfc_state%SST, diag, mask=G%mask2dT)
    if (CS%id_sssabs > 0) call post_data(CS%id_sssabs, sfc_state%SSS, diag, mask=G%mask2dT)
    !Using TEOS-10 function calls convert T&S diagnostics
    !from conservative temp to potential temp and
    !from absolute salinity to practical salinity
    do j=js,je ; do i=is,ie
      pracSal(i,j) = gsw_sp_from_sr(sfc_state%SSS(i,j))
      potTemp(i,j) = gsw_pt_from_ct(sfc_state%SSS(i,j),sfc_state%SST(i,j))
    enddo ; enddo
    if (CS%id_sst > 0) call post_data(CS%id_sst, potTemp, diag, mask=G%mask2dT)
    if (CS%id_sss > 0) call post_data(CS%id_sss, pracSal, diag, mask=G%mask2dT)
  endif

  if (CS%id_sst_sq > 0) then
    do j=js,je ; do i=is,ie
      SST_sq(i,j) = sfc_state%SST(i,j)*sfc_state%SST(i,j)
    enddo ; enddo
    call post_data(CS%id_sst_sq, SST_sq, diag, mask=G%mask2dT)
  endif
  if (CS%id_sss_sq > 0) then
    do j=js,je ; do i=is,ie
      SSS_sq(i,j) = sfc_state%SSS(i,j)*sfc_state%SSS(i,j)
    enddo ; enddo
    call post_data(CS%id_sss_sq, SSS_sq, diag, mask=G%mask2dT)
  endif

  if (CS%id_ssu > 0) &
    call post_data(CS%id_ssu, sfc_state%u, diag, mask=G%mask2dCu)
  if (CS%id_ssv > 0) &
    call post_data(CS%id_ssv, sfc_state%v, diag, mask=G%mask2dCv)

  if (CS%id_speed > 0) then
    do j=js,je ; do i=is,ie
      sfc_speed(i,j) = sqrt(0.5*(sfc_state%u(I-1,j)**2 + sfc_state%u(I,j)**2) + &
                            0.5*(sfc_state%v(i,J-1)**2 + sfc_state%v(i,J)**2))
    enddo ; enddo
    call post_data(CS%id_speed, sfc_speed, diag, mask=G%mask2dT)
  endif

  call coupler_type_send_data(sfc_state%tr_fields, get_diag_time_end(diag))

end subroutine post_surface_diagnostics

!> Sets a handle inside diagnostics mediator to associate 3d cell measures
subroutine register_cell_measure(G, diag, Time)
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean grid structure
  type(diag_ctrl), target, intent(inout) :: diag !< Regulates diagnostic output
  type(time_type),         intent(in)    :: Time !< Model time
  ! Local variables
  integer :: id
  id = register_diag_field('ocean_model', 'volcello', diag%axesTL, &
                           Time, 'Ocean grid-cell volume', 'm3', &
                           standard_name='ocean_volume', v_extensive=.true., &
                           x_cell_method='sum', y_cell_method='sum')
  call diag_associate_volume_cell_measure(diag, id)

end subroutine register_cell_measure

!> Offers the static fields in the ocean grid type
!! for output via the diag_manager.
subroutine write_static_fields(G, diag)
  type(ocean_grid_type),   intent(in)    :: G      !< ocean grid structure
  type(diag_ctrl), target, intent(inout) :: diag   !< regulates diagnostic output
  ! Local variables
  real    :: tmp_h(SZI_(G),SZJ_(G))
  integer :: id, i, j

  id = register_static_field('ocean_model', 'geolat', diag%axesT1, &
        'Latitude of tracer (T) points', 'degrees_north')
  if (id > 0) call post_data(id, G%geoLatT, diag, .true.)

  id = register_static_field('ocean_model', 'geolon', diag%axesT1, &
        'Longitude of tracer (T) points', 'degrees_east')
  if (id > 0) call post_data(id, G%geoLonT, diag, .true.)

  id = register_static_field('ocean_model', 'geolat_c', diag%axesB1, &
        'Latitude of corner (Bu) points', 'degrees_north', interp_method='none')
  if (id > 0) call post_data(id, G%geoLatBu, diag, .true.)

  id = register_static_field('ocean_model', 'geolon_c', diag%axesB1, &
        'Longitude of corner (Bu) points', 'degrees_east', interp_method='none')
  if (id > 0) call post_data(id, G%geoLonBu, diag, .true.)

  id = register_static_field('ocean_model', 'geolat_v', diag%axesCv1, &
        'Latitude of meridional velocity (Cv) points', 'degrees_north', interp_method='none')
  if (id > 0) call post_data(id, G%geoLatCv, diag, .true.)

  id = register_static_field('ocean_model', 'geolon_v', diag%axesCv1, &
        'Longitude of meridional velocity (Cv) points', 'degrees_east', interp_method='none')
  if (id > 0) call post_data(id, G%geoLonCv, diag, .true.)

  id = register_static_field('ocean_model', 'geolat_u', diag%axesCu1, &
        'Latitude of zonal velocity (Cu) points', 'degrees_north', interp_method='none')
  if (id > 0) call post_data(id, G%geoLatCu, diag, .true.)

  id = register_static_field('ocean_model', 'geolon_u', diag%axesCu1, &
        'Longitude of zonal velocity (Cu) points', 'degrees_east', interp_method='none')
  if (id > 0) call post_data(id, G%geoLonCu, diag, .true.)

  id = register_static_field('ocean_model', 'area_t', diag%axesT1,   &
        'Surface area of tracer (T) cells', 'm2',                    &
        cmor_field_name='areacello', cmor_standard_name='cell_area', &
        cmor_long_name='Ocean Grid-Cell Area',      &
        x_cell_method='sum', y_cell_method='sum', area_cell_method='sum')
  if (id > 0) then
    call post_data(id, G%areaT, diag, .true.)
    call diag_register_area_ids(diag, id_area_t=id)
  endif

  id = register_static_field('ocean_model', 'area_u', diag%axesCu1,     &
        'Surface area of x-direction flow (U) cells', 'm2',             &
        cmor_field_name='areacello_cu', cmor_standard_name='cell_area', &
        cmor_long_name='Ocean Grid-Cell Area',         &
        x_cell_method='sum', y_cell_method='sum', area_cell_method='sum')
  if (id > 0) then
    call post_data(id, G%areaCu, diag, .true.)
  endif

  id = register_static_field('ocean_model', 'area_v', diag%axesCv1,     &
        'Surface area of y-direction flow (V) cells', 'm2',             &
        cmor_field_name='areacello_cv', cmor_standard_name='cell_area', &
        cmor_long_name='Ocean Grid-Cell Area',         &
        x_cell_method='sum', y_cell_method='sum', area_cell_method='sum')
  if (id > 0) then
    call post_data(id, G%areaCv, diag, .true.)
  endif

  id = register_static_field('ocean_model', 'area_q', diag%axesB1,      &
        'Surface area of B-grid flow (Q) cells', 'm2',                  &
        cmor_field_name='areacello_bu', cmor_standard_name='cell_area', &
        cmor_long_name='Ocean Grid-Cell Area',         &
        x_cell_method='sum', y_cell_method='sum', area_cell_method='sum')
  if (id > 0) then
    call post_data(id, G%areaBu, diag, .true.)
  endif

  id = register_static_field('ocean_model', 'depth_ocean', diag%axesT1,  &
        'Depth of the ocean at tracer points', 'm',                      &
        standard_name='sea_floor_depth_below_geoid',                     &
        cmor_field_name='deptho', cmor_long_name='Sea Floor Depth',      &
        cmor_standard_name='sea_floor_depth_below_geoid',&
        area=diag%axesT1%id_area, &
        x_cell_method='mean', y_cell_method='mean', area_cell_method='mean')
  if (id > 0) call post_data(id, G%bathyT, diag, .true., mask=G%mask2dT)

  id = register_static_field('ocean_model', 'wet', diag%axesT1, &
        '0 if land, 1 if ocean at tracer points', 'none', area=diag%axesT1%id_area)
  if (id > 0) call post_data(id, G%mask2dT, diag, .true.)

  id = register_static_field('ocean_model', 'wet_c', diag%axesB1, &
        '0 if land, 1 if ocean at corner (Bu) points', 'none', interp_method='none')
  if (id > 0) call post_data(id, G%mask2dBu, diag, .true.)

  id = register_static_field('ocean_model', 'wet_u', diag%axesCu1, &
        '0 if land, 1 if ocean at zonal velocity (Cu) points', 'none', interp_method='none')
  if (id > 0) call post_data(id, G%mask2dCu, diag, .true.)

  id = register_static_field('ocean_model', 'wet_v', diag%axesCv1, &
        '0 if land, 1 if ocean at meridional velocity (Cv) points', 'none', interp_method='none')
  if (id > 0) call post_data(id, G%mask2dCv, diag, .true.)

  id = register_static_field('ocean_model', 'Coriolis', diag%axesB1, &
        'Coriolis parameter at corner (Bu) points', 's-1', interp_method='none')
  if (id > 0) call post_data(id, G%CoriolisBu, diag, .true.)

  id = register_static_field('ocean_model', 'dxt', diag%axesT1, &
        'Delta(x) at thickness/tracer points (meter)', 'm', interp_method='none')
  if (id > 0) call post_data(id, G%dxt, diag, .true.)

  id = register_static_field('ocean_model', 'dyt', diag%axesT1, &
        'Delta(y) at thickness/tracer points (meter)', 'm', interp_method='none')
  if (id > 0) call post_data(id, G%dyt, diag, .true.)

  id = register_static_field('ocean_model', 'dxCu', diag%axesCu1, &
        'Delta(x) at u points (meter)', 'm', interp_method='none')
  if (id > 0) call post_data(id, G%dxCu, diag, .true.)

  id = register_static_field('ocean_model', 'dyCu', diag%axesCu1, &
        'Delta(y) at u points (meter)', 'm', interp_method='none')
  if (id > 0) call post_data(id, G%dyCu, diag, .true.)

  id = register_static_field('ocean_model', 'dxCv', diag%axesCv1, &
        'Delta(x) at v points (meter)', 'm', interp_method='none')
  if (id > 0) call post_data(id, G%dxCv, diag, .true.)

  id = register_static_field('ocean_model', 'dyCv', diag%axesCv1, &
        'Delta(y) at v points (meter)', 'm', interp_method='none')
  if (id > 0) call post_data(id, G%dyCv, diag, .true.)

  ! This static diagnostic is from CF 1.8, and is the fraction of a cell
  ! covered by ocean, given as a percentage (poorly named).
  id = register_static_field('ocean_model', 'area_t_percent', diag%axesT1, &
        'Percentage of cell area covered by ocean', '%', &
        cmor_field_name='sftof', cmor_standard_name='SeaAreaFraction', &
        cmor_long_name='Sea Area Fraction', &
        x_cell_method='mean', y_cell_method='mean', area_cell_method='mean')
  if (id > 0) then
    tmp_h(:,:) = 0.
    tmp_h(G%isc:G%iec,G%jsc:G%jec) = 100. * G%mask2dT(G%isc:G%iec,G%jsc:G%jec)
    call post_data(id, tmp_h, diag, .true.)
  endif

end subroutine write_static_fields


!> Set the fields that are needed for bitwise identical restarting
!! the time stepping scheme.  In addition to those specified here
!! directly, there may be fields related to the forcing or to the
!! barotropic solver that are needed; these are specified in sub-
!! routines that are called from this one.
!!
!! This routine should be altered if there are any changes to the
!! time stepping scheme.  The CHECK_RESTART facility may be used to
!! confirm that all needed restart fields have been included.
subroutine set_restart_fields(GV, param_file, CS)
  type(verticalGrid_type),  intent(inout) :: GV         !< ocean vertical grid structure
  type(param_file_type),    intent(in) :: param_file    !< opened file for parsing to get parameters
  type(MOM_control_struct), intent(in) :: CS            !< control structure set up by inialize_MOM
  ! Local variables
  logical :: use_ice_shelf ! Needed to determine whether to add CS%Hml to restarts
  type(vardesc) :: vd
  character(len=48) :: thickness_units, flux_units

  call get_param(param_file, '', "ICE_SHELF", use_ice_shelf, default=.false., do_not_log=.true.)

  thickness_units = get_thickness_units(GV)
  flux_units = get_flux_units(GV)

  if (CS%use_temperature) then
    vd = var_desc("Temp","degC","Potential Temperature")
    call register_restart_field(CS%tv%T, vd, .true., CS%restart_CSp)

    vd = var_desc("Salt","PPT","Salinity")
    call register_restart_field(CS%tv%S, vd, .true., CS%restart_CSp)
  endif

  vd = var_desc("h",thickness_units,"Layer Thickness")
  call register_restart_field(CS%h, vd, .true., CS%restart_CSp)

  vd = var_desc("u","m s-1","Zonal velocity",'u','L')
  call register_restart_field(CS%u, vd, .true., CS%restart_CSp)

  vd = var_desc("v","m s-1","Meridional velocity",'v','L')
  call register_restart_field(CS%v, vd, .true., CS%restart_CSp)

  if (CS%use_frazil) then
    vd = var_desc("frazil","J m-2","Frazil heat flux into ocean",'h','1')
    call register_restart_field(CS%tv%frazil, vd, .false., CS%restart_CSp)
  endif

  if (CS%interp_p_surf) then
    vd = var_desc("p_surf_prev","Pa","Previous ocean surface pressure",'h','1')
    call register_restart_field(CS%p_surf_prev, vd, .false., CS%restart_CSp)
  endif

  vd = var_desc("ave_ssh","meter","Time average sea surface height",'h','1')
  call register_restart_field(CS%ave_ssh, vd, .false., CS%restart_CSp)

  ! hML is needed when using the ice shelf module
  if (use_ice_shelf .and. associated(CS%Hml)) then
     vd = var_desc("hML","meter","Mixed layer thickness",'h','1')
     call register_restart_field(CS%Hml, vd, .false., CS%restart_CSp)
  endif

end subroutine set_restart_fields

!> This subroutine applies a correction to the sea surface height to compensate
!! for the atmospheric pressure (the inverse barometer).
subroutine adjust_ssh_for_p_atm(CS, G, GV, ssh, p_atm)
  type(MOM_control_struct),          intent(in)    :: CS     !< control structure
  type(ocean_grid_type),             intent(in)    :: G      !< ocean grid structure
  type(verticalGrid_type),           intent(in)    :: GV     !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: ssh    !< time mean surface height (m)
  real, dimension(:,:),    optional, pointer       :: p_atm  !< atmospheric pressure (Pascal)

  real :: Rho_conv    ! The density used to convert surface pressure to
                      ! a corrected effective SSH, in kg m-3.
  real :: IgR0        ! The SSH conversion factor from Pa to m.
  integer :: i, j, is, ie, js, je

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec
  if (ASSOCIATED(p_atm)) then
    ! Correct the output sea surface height for the contribution from the
    ! atmospheric pressure
    do j=js,je ; do i=is,ie
      if ((ASSOCIATED(CS%tv%eqn_of_state)) .and. (CS%calc_rho_for_sea_lev)) then
        call calculate_density(CS%tv%T(i,j,1), CS%tv%S(i,j,1), p_atm(i,j)/2.0, &
                               Rho_conv, CS%tv%eqn_of_state)
      else
        Rho_conv=GV%Rho0
      endif
      IgR0 = 1.0 / (Rho_conv * GV%g_Earth)
      ssh(i,j) = ssh(i,j) + p_atm(i,j) * IgR0
    enddo ; enddo
  endif

end subroutine adjust_ssh_for_p_atm

!> This subroutine allocates the fields for the surface (return) properties of
!! the ocean model.  Unused fields are unallocated.
subroutine allocate_surface_state(sfc_state, G, use_temperature, do_integrals, &
                                  gas_fields_ocn)
  type(ocean_grid_type), intent(in)    :: G                !< ocean grid structure
  type(surface),         intent(inout) :: sfc_state        !< ocean surface state type to be allocated.
  logical,     optional, intent(in)    :: use_temperature  !< If true, allocate the space for thermodynamic variables.
  logical,     optional, intent(in)    :: do_integrals     !< If true, allocate the space for vertically integrated fields.
  type(coupler_1d_bc_type), &
               optional, intent(in)    :: gas_fields_ocn   !< If present, this type describes the ocean
                                              !! ocean and surface-ice fields that will participate
                                              !! in the calculation of additional gas or other
                                              !! tracer fluxes, and can be used to spawn related
                                              !! internal variables in the ice model.

  logical :: use_temp, alloc_integ
  integer :: is, ie, js, je, isd, ied, jsd, jed
  integer :: isdB, iedB, jsdB, jedB

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isdB = G%isdB ; iedB = G%iedB; jsdB = G%jsdB ; jedB = G%jedB

  use_temp = .true. ; if (present(use_temperature)) use_temp = use_temperature
  alloc_integ = .true. ; if (present(do_integrals)) alloc_integ = do_integrals

  if (sfc_state%arrays_allocated) return

  if (use_temp) then
    allocate(sfc_state%SST(isd:ied,jsd:jed)) ; sfc_state%SST(:,:) = 0.0
    allocate(sfc_state%SSS(isd:ied,jsd:jed)) ; sfc_state%SSS(:,:) = 0.0
  else
    allocate(sfc_state%sfc_density(isd:ied,jsd:jed)) ; sfc_state%sfc_density(:,:) = 0.0
  endif
  allocate(sfc_state%sea_lev(isd:ied,jsd:jed)) ; sfc_state%sea_lev(:,:) = 0.0
  allocate(sfc_state%Hml(isd:ied,jsd:jed)) ; sfc_state%Hml(:,:) = 0.0
  allocate(sfc_state%u(IsdB:IedB,jsd:jed)) ; sfc_state%u(:,:) = 0.0
  allocate(sfc_state%v(isd:ied,JsdB:JedB)) ; sfc_state%v(:,:) = 0.0

  if (alloc_integ) then
    ! Allocate structures for the vertically integrated ocean_mass, ocean_heat,
    ! and ocean_salt.
    allocate(sfc_state%ocean_mass(isd:ied,jsd:jed)) ; sfc_state%ocean_mass(:,:) = 0.0
    if (use_temp) then
      allocate(sfc_state%ocean_heat(isd:ied,jsd:jed)) ; sfc_state%ocean_heat(:,:) = 0.0
      allocate(sfc_state%ocean_salt(isd:ied,jsd:jed)) ; sfc_state%ocean_salt(:,:) = 0.0
    endif
    allocate(sfc_state%salt_deficit(isd:ied,jsd:jed)) ; sfc_state%salt_deficit(:,:) = 0.0
  endif

  if (present(gas_fields_ocn)) &
    call coupler_type_spawn(gas_fields_ocn, sfc_state%tr_fields, &
                            (/is,is,ie,ie/), (/js,js,je,je/), as_needed=.true.)

  sfc_state%arrays_allocated = .true.

end subroutine allocate_surface_state

!> This subroutine sets the surface (return) properties of the ocean
!! model by setting the appropriate fields in state.  Unused fields
!! are set to NULL or are unallocated.
subroutine calculate_surface_state(sfc_state, u, v, h, ssh, G, GV, CS)
  type(ocean_grid_type),                     intent(inout) :: G      !< ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV     !< ocean vertical grid structure
  type(surface),                             intent(inout) :: sfc_state !< ocean surface state
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u      !< zonal velocity (m/s)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v      !< meridional velocity (m/s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h      !< layer thickness (m or kg/m2)
  real, dimension(SZI_(G),SZJ_(G)),          intent(in)    :: ssh    !< time mean surface height (m)
  type(MOM_control_struct),                  intent(inout) :: CS     !< control structure

  ! local
  real :: depth(SZI_(G))              ! distance from the surface (meter)
  real :: depth_ml                    ! depth over which to average to
                                      ! determine mixed layer properties (meter)
  real :: dh                          ! thickness of a layer within mixed layer (meter)
  real :: mass                        ! mass per unit area of a layer (kg/m2)

  real :: hu, hv
  integer :: i, j, k, is, ie, js, je, nz, numberOfErrors
  integer :: isd, ied, jsd, jed
  integer :: iscB, iecB, jscB, jecB, isdB, iedB, jsdB, jedB
  logical :: localError
  character(240) :: msg

  call callTree_enter("calculate_surface_state(), MOM.F90")
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  iscB = G%iscB ; iecB = G%iecB; jscB = G%jscB ; jecB = G%jecB
  isdB = G%isdB ; iedB = G%iedB; jsdB = G%jsdB ; jedB = G%jedB

  if (.not.sfc_state%arrays_allocated) then
    !  Consider using a run-time flag to determine whether to do the vertical
    ! integrals, since the 3-d sums are not negligible in cost.
    call allocate_surface_state(sfc_state, G, CS%use_temperature, do_integrals=.true.)
  endif
  sfc_state%frazil => CS%tv%frazil
  sfc_state%TempxPmE => CS%tv%TempxPmE
  sfc_state%internal_heat => CS%tv%internal_heat
  if (associated(CS%visc%taux_shelf)) sfc_state%taux_shelf => CS%visc%taux_shelf
  if (associated(CS%visc%tauy_shelf)) sfc_state%tauy_shelf => CS%visc%tauy_shelf

  do j=js,je ; do i=is,ie
    sfc_state%sea_lev(i,j) = ssh(i,j)
  enddo ; enddo

  if (CS%bulkmixedlayer) then
    if (CS%use_temperature) then ; do j=js,je ; do i=is,ie
      sfc_state%SST(i,j) = CS%tv%T(i,j,1)
      sfc_state%SSS(i,j) = CS%tv%S(i,j,1)
    enddo ; enddo ; endif
    do j=js,je ; do I=IscB,IecB
      sfc_state%u(I,j) = u(I,j,1)
    enddo ; enddo
    do J=JscB,JecB ; do i=is,ie
      sfc_state%v(i,J) = v(i,J,1)
    enddo ; enddo

    if (associated(CS%Hml)) then ; do j=js,je ; do i=is,ie
      sfc_state%Hml(i,j) = CS%Hml(i,j)
    enddo ; enddo ; endif
  else

    depth_ml = CS%Hmix
  !   Determine the mean tracer properties of the uppermost depth_ml fluid.
    !$OMP parallel do default(shared) private(depth,dh)
    do j=js,je
      do i=is,ie
        depth(i) = 0.0
        if (CS%use_temperature) then
          sfc_state%SST(i,j) = 0.0 ; sfc_state%SSS(i,j) = 0.0
        else
          sfc_state%sfc_density(i,j) = 0.0
        endif
      enddo

      do k=1,nz ; do i=is,ie
        if (depth(i) + h(i,j,k)*GV%H_to_m < depth_ml) then
          dh = h(i,j,k)*GV%H_to_m
        elseif (depth(i) < depth_ml) then
          dh = depth_ml - depth(i)
        else
          dh = 0.0
        endif
        if (CS%use_temperature) then
          sfc_state%SST(i,j) = sfc_state%SST(i,j) + dh * CS%tv%T(i,j,k)
          sfc_state%SSS(i,j) = sfc_state%SSS(i,j) + dh * CS%tv%S(i,j,k)
        else
          sfc_state%sfc_density(i,j) = sfc_state%sfc_density(i,j) + dh * GV%Rlay(k)
        endif
        depth(i) = depth(i) + dh
      enddo ; enddo
  ! Calculate the average properties of the mixed layer depth.
      do i=is,ie
        if (depth(i) < GV%H_subroundoff*GV%H_to_m) &
            depth(i) = GV%H_subroundoff*GV%H_to_m
        if (CS%use_temperature) then
          sfc_state%SST(i,j) = sfc_state%SST(i,j) / depth(i)
          sfc_state%SSS(i,j) = sfc_state%SSS(i,j) / depth(i)
        else
          sfc_state%sfc_density(i,j) = sfc_state%sfc_density(i,j) / depth(i)
        endif
        sfc_state%Hml(i,j) = depth(i)
      enddo
    enddo ! end of j loop

!   Determine the mean velocities in the uppermost depth_ml fluid.
    if (CS%Hmix_UV>0.) then
      depth_ml = CS%Hmix_UV
      !$OMP parallel do default(shared) private(depth,dh,hv)
      do J=jscB,jecB
        do i=is,ie
          depth(i) = 0.0
          sfc_state%v(i,J) = 0.0
        enddo
        do k=1,nz ; do i=is,ie
          hv = 0.5 * (h(i,j,k) + h(i,j+1,k)) * GV%H_to_m
          if (depth(i) + hv < depth_ml) then
            dh = hv
          elseif (depth(i) < depth_ml) then
            dh = depth_ml - depth(i)
          else
            dh = 0.0
          endif
          sfc_state%v(i,J) = sfc_state%v(i,J) + dh * v(i,J,k)
          depth(i) = depth(i) + dh
        enddo ; enddo
        ! Calculate the average properties of the mixed layer depth.
        do i=is,ie
          if (depth(i) < GV%H_subroundoff*GV%H_to_m) &
              depth(i) = GV%H_subroundoff*GV%H_to_m
          sfc_state%v(i,J) = sfc_state%v(i,J) / depth(i)
        enddo
      enddo ! end of j loop

      !$OMP parallel do default(shared) private(depth,dh,hu)
      do j=js,je
        do I=iscB,iecB
          depth(I) = 0.0
          sfc_state%u(I,j) = 0.0
        enddo
        do k=1,nz ; do I=iscB,iecB
          hu = 0.5 * (h(i,j,k) + h(i+1,j,k)) * GV%H_to_m
          if (depth(i) + hu < depth_ml) then
            dh = hu
          elseif (depth(I) < depth_ml) then
            dh = depth_ml - depth(I)
          else
            dh = 0.0
          endif
          sfc_state%u(I,j) = sfc_state%u(I,j) + dh * u(I,j,k)
          depth(I) = depth(I) + dh
        enddo ; enddo
        ! Calculate the average properties of the mixed layer depth.
        do I=iscB,iecB
          if (depth(I) < GV%H_subroundoff*GV%H_to_m) &
              depth(I) = GV%H_subroundoff*GV%H_to_m
          sfc_state%u(I,j) = sfc_state%u(I,j) / depth(I)
        enddo
      enddo ! end of j loop
    else ! Hmix_UV<=0.
      do j=js,je ; do I=IscB,IecB
        sfc_state%u(I,j) = u(I,j,1)
      enddo ; enddo
      do J=JscB,JecB ; do i=is,ie
        sfc_state%v(i,J) = v(i,J,1)
      enddo ; enddo
    endif
  endif                                              ! end BULKMIXEDLAYER

  if (allocated(sfc_state%salt_deficit) .and. associated(CS%tv%salt_deficit)) then
    !$OMP parallel do default(shared)
    do j=js,je ; do i=is,ie
      ! Convert from gSalt to kgSalt
      sfc_state%salt_deficit(i,j) = 1000.0 * CS%tv%salt_deficit(i,j)
    enddo ; enddo
  endif

  if (allocated(sfc_state%ocean_mass) .and. allocated(sfc_state%ocean_heat) .and. &
      allocated(sfc_state%ocean_salt)) then
    !$OMP parallel do default(shared)
    do j=js,je ; do i=is,ie
      sfc_state%ocean_mass(i,j) = 0.0
      sfc_state%ocean_heat(i,j) = 0.0 ; sfc_state%ocean_salt(i,j) = 0.0
    enddo ; enddo
    !$OMP parallel do default(shared) private(mass)
    do j=js,je ; do k=1,nz; do i=is,ie
      mass = GV%H_to_kg_m2*h(i,j,k)
      sfc_state%ocean_mass(i,j) = sfc_state%ocean_mass(i,j) + mass
      sfc_state%ocean_heat(i,j) = sfc_state%ocean_heat(i,j) + mass*CS%tv%T(i,j,k)
      sfc_state%ocean_salt(i,j) = sfc_state%ocean_salt(i,j) + &
                              mass * (1.0e-3*CS%tv%S(i,j,k))
    enddo ; enddo ; enddo
  else
    if (allocated(sfc_state%ocean_mass)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie ; sfc_state%ocean_mass(i,j) = 0.0 ; enddo ; enddo
      !$OMP parallel do default(shared)
      do j=js,je ; do k=1,nz ; do i=is,ie
        sfc_state%ocean_mass(i,j) = sfc_state%ocean_mass(i,j) + GV%H_to_kg_m2*h(i,j,k)
      enddo ; enddo ; enddo
    endif
    if (allocated(sfc_state%ocean_heat)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie ; sfc_state%ocean_heat(i,j) = 0.0 ; enddo ; enddo
      !$OMP parallel do default(shared) private(mass)
      do j=js,je ; do k=1,nz ; do i=is,ie
        mass = GV%H_to_kg_m2*h(i,j,k)
        sfc_state%ocean_heat(i,j) = sfc_state%ocean_heat(i,j) + mass*CS%tv%T(i,j,k)
      enddo ; enddo ; enddo
    endif
    if (allocated(sfc_state%ocean_salt)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie ; sfc_state%ocean_salt(i,j) = 0.0 ; enddo ; enddo
      !$OMP parallel do default(shared) private(mass)
      do j=js,je ; do k=1,nz ; do i=is,ie
        mass = GV%H_to_kg_m2*h(i,j,k)
        sfc_state%ocean_salt(i,j) = sfc_state%ocean_salt(i,j) + &
                                mass * (1.0e-3*CS%tv%S(i,j,k))
      enddo ; enddo ; enddo
    endif
  endif

  if (associated(CS%tracer_flow_CSp)) then
    call call_tracer_surface_state(sfc_state, h, G, CS%tracer_flow_CSp)
  endif

  if (CS%check_bad_surface_vals) then
    numberOfErrors=0 ! count number of errors
    do j=js,je; do i=is,ie
      if (G%mask2dT(i,j)>0.) then
        localError = sfc_state%sea_lev(i,j)<=-G%bathyT(i,j)       &
                .or. sfc_state%sea_lev(i,j)>= CS%bad_val_ssh_max  &
                .or. sfc_state%sea_lev(i,j)<=-CS%bad_val_ssh_max  &
                .or. sfc_state%sea_lev(i,j)+G%bathyT(i,j) < CS%bad_val_column_thickness
        if (CS%use_temperature) localError = localError &
                .or. sfc_state%SSS(i,j)<0.                        &
                .or. sfc_state%SSS(i,j)>=CS%bad_val_sss_max       &
                .or. sfc_state%SST(i,j)< CS%bad_val_sst_min       &
                .or. sfc_state%SST(i,j)>=CS%bad_val_sst_max
        if (localError) then
          numberOfErrors=numberOfErrors+1
          if (numberOfErrors<9) then ! Only report details for the first few errors
            if (CS%use_temperature) then
              write(msg(1:240),'(2(a,i4,x),2(a,f8.3,x),8(a,es11.4,x))') &
                'Extreme surface sfc_state detected: i=',i,'j=',j, &
                'x=',G%geoLonT(i,j),'y=',G%geoLatT(i,j), &
                'D=',G%bathyT(i,j),                      &
                'SSH=',sfc_state%sea_lev(i,j),           &
                'SST=',sfc_state%SST(i,j),               &
                'SSS=',sfc_state%SSS(i,j),               &
                'U-=',sfc_state%u(I-1,j),                &
                'U+=',sfc_state%u(I,j),                  &
                'V-=',sfc_state%v(i,J-1),                &
                'V+=',sfc_state%v(i,J)
            else
              write(msg(1:240),'(2(a,i4,x),2(a,f8.3,x),6(a,es11.4))') &
                'Extreme surface sfc_state detected: i=',i,'j=',j, &
                'x=',G%geoLonT(i,j),'y=',G%geoLatT(i,j), &
                'D=',G%bathyT(i,j),                      &
                'SSH=',sfc_state%sea_lev(i,j),           &
                'U-=',sfc_state%u(I-1,j),                &
                'U+=',sfc_state%u(I,j),                  &
                'V-=',sfc_state%v(i,J-1),                &
                'V+=',sfc_state%v(i,J)
            endif
            call MOM_error(WARNING, trim(msg), all_print=.true.)
          elseif (numberOfErrors==9) then ! Indicate once that there are more errors
            call MOM_error(WARNING, 'There were more unreported extreme events!', all_print=.true.)
          endif ! numberOfErrors
        endif ! localError
      endif ! mask2dT
    enddo; enddo
    call sum_across_PEs(numberOfErrors)
    if (numberOfErrors>0) then
      write(msg(1:240),'(3(a,i9,x))') 'There were a total of ',numberOfErrors, &
          'locations detected with extreme surface values!'
      call MOM_error(FATAL, trim(msg))
    endif
  endif

  call callTree_leave("calculate_surface_sfc_state()")
end subroutine calculate_surface_state


!> End of model
subroutine MOM_end(CS)
  type(MOM_control_struct), pointer :: CS   !< MOM control structure

  if (CS%use_ALE_algorithm) then
    call ALE_end(CS%ALE_CSp)
  endif

  DEALLOC_(CS%u) ; DEALLOC_(CS%v) ; DEALLOC_(CS%h)
  DEALLOC_(CS%uh) ; DEALLOC_(CS%vh)

  if (CS%use_temperature) then
    DEALLOC_(CS%T) ; CS%tv%T => NULL() ; DEALLOC_(CS%S) ; CS%tv%S => NULL()
  endif
  if (associated(CS%tv%frazil)) deallocate(CS%tv%frazil)
  if (associated(CS%tv%salt_deficit)) deallocate(CS%tv%salt_deficit)
  if (associated(CS%Hml)) deallocate(CS%Hml)

  call tracer_advect_end(CS%tracer_adv_CSp)
  call tracer_hor_diff_end(CS%tracer_diff_CSp)
  call tracer_registry_end(CS%tracer_Reg)
  call tracer_flow_control_end(CS%tracer_flow_CSp)

  if(CS%offline_tracer_mode) then
    call offline_transport_end(CS%offline_CSp)
  endif

  DEALLOC_(CS%uhtr) ; DEALLOC_(CS%vhtr)
  if (CS%split) then
    call end_dyn_split_RK2(CS%dyn_split_RK2_CSp)
  elseif (CS%use_RK2) then
    call end_dyn_unsplit_RK2(CS%dyn_unsplit_RK2_CSp)
  else
    call end_dyn_unsplit(CS%dyn_unsplit_CSp)
  endif
  DEALLOC_(CS%ave_ssh)
  if (associated(CS%update_OBC_CSp)) call OBC_register_end(CS%update_OBC_CSp)

  call verticalGridEnd(CS%GV)
  call MOM_grid_end(CS%G)

  deallocate(CS)

end subroutine MOM_end

!> \namespace mom
!!
!! Modular Ocean Model (MOM) Version 6.0 (MOM6)
!!
!! \authors Alistair Adcroft, Robert Hallberg, and Stephen Griffies
!!
!!  Additional contributions from:
!!    * Whit Anderson
!!    * Brian Arbic
!!    * Will Cooke
!!    * Anand Gnanadesikan
!!    * Matthew Harrison
!!    * Mehmet Ilicak
!!    * Laura Jackson
!!    * Jasmine John
!!    * John Krasting
!!    * Zhi Liang
!!    * Bonnie Samuels
!!    * Harper Simmons
!!    * Laurent White
!!    * Niki Zadeh
!!
!!  MOM ice-shelf code was developed by
!!  * Daniel Goldberg
!!  * Robert Hallberg
!!  * Chris Little
!!  * Olga Sergienko
!!
!!  \section section_overview Overview of MOM
!!
!!  This program (MOM) simulates the ocean by numerically solving
!!  the hydrostatic primitive equations in generalized Lagrangian
!!  vertical coordinates, typically tracking stretched pressure (p*)
!!  surfaces or following isopycnals in the ocean's interior, and
!!  general orthogonal horizontal coordinates. Unlike earlier versions
!!  of MOM, in MOM6 these equations are horizontally discretized on an
!!  Arakawa C-grid.  (It remains to be seen whether a B-grid dynamic
!!  core will be revived in MOM6 at a later date; for now applications
!!  requiring a B-grid discretization should use MOM5.1.)  MOM6 offers
!!  a range of options for the physical parameterizations, from those
!!  most appropriate to highly idealized models for geophysical fluid
!!  dynamics studies to a rich suite of processes appropriate for
!!  realistic ocean simulations.  The thermodynamic options typically
!!  use conservative temperature and preformed salinity as conservative
!!  state variables and a full nonlinear equation of state, but there
!!  are also idealized adiabatic configurations of the model that use
!!  fixed density layers.  Version 6.0 of MOM continues in the long
!!  tradition of a commitment to climate-quality ocean simulations
!!  embodied in previous versions of MOM, even as it draws extensively
!!  on the lessons learned in the development of the Generalized Ocean
!!  Layered Dynamics (GOLD) ocean model, which was also primarily
!!  developed at NOAA/GFDL.  MOM has also benefited tremendously from
!!  the FMS infrastructure, which it utilizes and shares with other
!!  component models developed at NOAA/GFDL.
!!
!!    When run is isopycnal-coordinate mode, the uppermost few layers
!!  are often used to describe a bulk mixed layer, including the
!!  effects of penetrating shortwave radiation.  Either a split-
!!  explicit time stepping scheme or a non-split scheme may be used
!!  for the dynamics, while the time stepping may be split (and use
!!  different numbers of steps to cover the same interval) for the
!!  forcing, the thermodynamics, and for the dynamics.  Most of the
!!  numerics are second order accurate in space.  MOM can run with an
!!  absurdly thin minimum layer thickness. A variety of non-isopycnal
!!  vertical coordinate options are under development, but all exploit
!!  the advantages of a Lagrangian vertical coordinate, as discussed
!!  in detail by Adcroft and Hallberg (Ocean Modelling, 2006).
!!
!!    Details of the numerics and physical parameterizations are
!!  provided in the appropriate source files.  All of the available
!!  options are selected at run-time by parsing the input files,
!!  usually MOM_input and MOM_override, and the options choices are
!!  then documented for each run in MOM_param_docs.
!!
!!    MOM6 integrates the equations forward in time in three distinct
!!  phases.  In one phase, the dynamic equations for the velocities
!!  and layer thicknesses are advanced, capturing the propagation of
!!  external and internal inertia-gravity waves, Rossby waves, and
!!  other strictly adiabatic processes, including lateral stresses,
!!  vertical viscosity and momentum forcing, and interface height
!!  diffusion (commonly called Gent-McWilliams diffusion in depth-
!!  coordinate models).  In the second phase, all tracers are advected
!!  and diffused along the layers.  The third phase applies diabatic
!!  processes, vertical mixing of water properties, and perhaps
!!  vertical remapping to cause the layers to track the desired
!!  vertical coordinate.
!!
!!    The present file (MOM.F90) orchestrates the main time stepping
!!  loops. One time integration option for the dynamics uses a split
!!  explicit time stepping scheme to rapidly step the barotropic
!!  pressure and velocity fields. The barotropic velocities are
!!  averaged over the baroclinic time step before they are used to
!!  advect thickness and determine the baroclinic accelerations.  As
!!  described in Hallberg and Adcroft (2009), a barotropic correction
!!  is applied to the time-mean layer velocities to ensure that the
!!  sum of the layer transports agrees with the time-mean barotropic
!!  transport, thereby ensuring that the estimates of the free surface
!!  from the sum of the layer thicknesses agrees with the final free
!!  surface height as calculated by the barotropic solver.  The
!!  barotropic and baroclinic velocities are kept consistent by
!!  recalculating the barotropic velocities from the baroclinic
!!  transports each time step. This scheme is described in Hallberg,
!!  1997, J. Comp. Phys. 135, 54-65 and in Hallberg and Adcroft, 2009,
!!  Ocean Modelling, 29, 15-26.
!!
!!    The other time integration options use non-split time stepping
!!  schemes based on the 3-step third order Runge-Kutta scheme
!!  described in Matsuno, 1966, J. Met. Soc. Japan, 44, 85-88, or on
!!  a two-step quasi-2nd order Runge-Kutta scheme.  These are much
!!  slower than the split time-stepping scheme, but they are useful
!!  for providing a more robust solution for debugging cases where the
!!  more complicated split time-stepping scheme may be giving suspect
!!  solutions.
!!
!!    There are a range of closure options available.  Horizontal
!!  velocities are subject to a combination of horizontal biharmonic
!!  and Laplacian friction (based on a stress tensor formalism) and a
!!  vertical Fickian viscosity (perhaps using the kinematic viscosity
!!  of water).  The horizontal viscosities may be constant, spatially
!!  varying or may be dynamically calculated using Smagorinsky's
!!  approach.  A diapycnal diffusion of density and thermodynamic
!!  quantities is also allowed, but not required, as is horizontal
!!  diffusion of interface heights (akin to the Gent-McWilliams
!!  closure of geopotential coordinate models).  The diapycnal mixing
!!  may use a fixed diffusivity or it may use the shear Richardson
!!  number dependent closure, like that described in Jackson et al.
!!  (JPO, 2008).  When there is diapycnal diffusion, it applies to
!!  momentum as well. As this is in addition to the vertical viscosity,
!!  the vertical Prandtl always exceeds 1.  A refined bulk-mixed layer
!!  is often used to describe the planetary boundary layer in realistic
!!  ocean simulations.
!!
!!    MOM has a number of noteworthy debugging capabilities.
!!  Excessively large velocities are truncated and MOM will stop
!!  itself after a number of such instances to keep the model from
!!  crashing altogether.  This is useful in diagnosing failures,
!!  or (by accepting some truncations) it may be useful for getting
!!  the model past the adjustment from an ill-balanced initial
!!  condition.  In addition, all of the accelerations in the columns
!!  with excessively large velocities may be directed to a text file.
!!  Parallelization errors may be diagnosed using the DEBUG option,
!!  which causes extensive checksums to be written out along with
!!  comments indicating where in the algorithm the sums originate and
!!  what variable is being summed.  The point where these checksums
!!  differ between runs is usually a good indication of where in the
!!  code the problem lies.  All of the test cases provided with MOM
!!  are routinely tested to ensure that they give bitwise identical
!!  results regardless of the domain decomposition, or whether they
!!  use static or dynamic memory allocation.
!!
!!  \section section_structure Structure of MOM
!!
!!  About 115 other files of source code and 4 header files comprise
!!  the MOM code, although there are several hundred more files that
!!  make up the FMS infrastructure upon which MOM is built.  Each of
!!  the MOM files contains comments documenting what it does, and
!!  most of the file names are fairly self-evident. In addition, all
!!  subroutines and data types are referenced via a module use, only
!!  statement, and the module names are consistent with the file names,
!!  so it is not too hard to find the source file for a subroutine.
!!
!!    The typical MOM directory tree is as follows:
!!
!! \verbatim
!!        ../MOM
!!        |-- config_src
!!        |   |-- coupled_driver
!!        |   |-- dynamic
!!        |   `-- solo_driver
!!        |-- examples
!!        |   |-- CM2G
!!        |   |-- ...
!!        |   `-- torus_advection_test
!!        `-- src
!!            |-- core
!!            |-- diagnostics
!!            |-- equation_of_state
!!            |-- framework
!!            |-- ice_shelf
!!            |-- initialization
!!            |-- parameterizations
!!            |   |-- lateral
!!            |   `-- vertical
!!            |-- tracer
!!            `-- user
!! \endverbatim
!!
!!  Rather than describing each file here, each directory contents
!!  will be described to give a broad overview of the MOM code
!!  structure.
!!
!!    The directories under config_src contain files that are used for
!!  configuring the code, for instance for coupled or ocean-only runs.
!!  Only one or two of these directories are used in compiling any,
!!  particular run.
!!
!!  * config_src/coupled_driver:
!!    The files here are used to couple MOM as a component in a larger
!!    run driven by the FMS coupler.  This includes code that converts
!!    various forcing fields into the code structures and flux and unit
!!    conventions used by MOM, and converts the MOM surface fields
!!    back to the forms used by other FMS components.
!!
!!  * config_src/dynamic:
!!    The only file here is the version of MOM_memory.h that is used
!!    for dynamic memory configurations of MOM.
!!
!!  * config_src/solo_driver:
!!    The files here are include the _main driver that is used when
!!    MOM is configured as an ocean-only model, as well as the files
!!    that specify the surface forcing in this configuration.
!!
!!    The directories under examples provide a large number of working
!!  configurations of MOM, along with reference solutions for several
!!  different compilers on GFDL's latest large computer.  The versions
!!  of MOM_memory.h in these directories need not be used if dynamic
!!  memory allocation is desired, and the answers should be unchanged.
!!
!!    The directories under src contain most of the MOM files.  These
!!  files are used in every configuration using MOM.
!!
!!  * src/core:
!!    The files here constitute the MOM dynamic core.  This directory
!!    also includes files with the types that describe the model's
!!    lateral grid and have defined types that are shared across
!!    various MOM modules to allow for more succinct and flexible
!!    subroutine argument lists.
!!
!!  * src/diagnostics:
!!    The files here calculate various diagnostics that are anciliary
!!    to the model itself.  While most of these diagnostics do not
!!    directly affect the model's solution, there are some, like the
!!    calculation of the deformation radius, that are used in some
!!    of the process parameterizations.
!!
!!  * src/equation_of_state:
!!    These files describe the physical properties of sea-water,
!!    including both the equation of state and when it freezes.
!!
!!  * src/framework:
!!    These files provide infrastructure utilities for MOM.  Many are
!!    simply wrappers for capabilities provided by FMS, although others
!!    provide capabilities (like the file_parser) that are unique to
!!    MOM. When MOM is adapted to use a modeling infrastructure
!!    distinct from FMS, most of the required changes are in this
!!    directory.
!!
!!  * src/initialization:
!!    These are the files that are used to initialize the MOM grid
!!    or provide the initial physical state for MOM.  These files are
!!    not intended to be modified, but provide a means for calling
!!    user-specific initialization code like the examples in src/user.
!!
!!  * src/parameterizations/lateral:
!!    These files implement a number of quasi-lateral (along-layer)
!!    process parameterizations, including lateral viscosities,
!!    parameterizations of eddy effects, and the calculation of tidal
!!    forcing.
!!
!!  * src/parameterizations/vertical:
!!    These files implement a number of vertical mixing or diabatic
!!    processes, including the effects of vertical viscosity and
!!    code to parameterize the planetary boundary layer.  There is a
!!    separate driver that orchestrates this portion of the algorithm,
!!    and there is a diversity of parameterizations to be found here.
!!
!!  * src/tracer:
!!    These files handle the lateral transport and diffusion of
!!    tracers, or are the code to implement various passive tracer
!!    packages.  Additional tracer packages are readily accomodated.
!!
!!  * src/user:
!!    These are either stub routines that a user could use to change
!!    the model's initial conditions or forcing, or are examples that
!!    implement specific test cases.  These files can easily  be hand
!!    edited to create new analytically specified configurations.
!!
!!
!!  Most simulations can be set up by modifying only the files
!!  MOM_input, and possibly one or two of the files in src/user.
!!  In addition, the diag_table (MOM_diag_table) will commonly be
!!  modified to tailor the output to the needs of the question at
!!  hand.  The FMS utility mkmf works with a file called path_names
!!  to build an appropriate makefile, and path_names should be edited
!!  to reflect the actual location of the desired source code.
!!
!!
!!  There are 3 publicly visible subroutines in this file (MOM.F90).
!!  * step_MOM steps MOM over a specified interval of time.
!!  * MOM_initialize calls initialize and does other initialization
!!    that does not warrant user modification.
!!  * calculate_surface_state determines the surface (bulk mixed layer
!!    if traditional isoycnal vertical coordinate) properties of the
!!    current model state and packages pointers to these fields into an
!!    exported structure.
!!
!!    The remaining subroutines in this file (src/core/MOM.F90) are:
!!  * find_total_transport determines the barotropic mass transport.
!!  * register_diags registers many diagnostic fields for the dynamic
!!    solver, or of the main model variables.
!!  * MOM_timing_init initializes various CPU time clocks.
!!  * write_static_fields writes out various time-invariant fields.
!!  * set_restart_fields is used to specify those fields that are
!!    written to and read from the restart file.
!!
!!  \section section_heat_budget Diagnosing MOM heat budget
!!
!!  Here are some example heat budgets for the ALE version of MOM6.
!!
!!  \subsection subsection_2d_heat_budget Depth integrated heat budget
!!
!!  Depth integrated heat budget diagnostic for MOM.
!!
!! * OPOTTEMPTEND_2d = T_ADVECTION_XY_2d + OPOTTEMPPMDIFF_2d + HFDS + HFGEOU
!!
!! * T_ADVECTION_XY_2d = horizontal advection
!! * OPOTTEMPPMDIFF_2d = neutral diffusion
!! * HFDS              = net surface boundary heat flux
!! * HFGEOU            = geothermal heat flux
!!
!! * HFDS = net surface boundary heat flux entering the ocean
!!        = rsntds + rlntds + hfls + hfss + heat_pme + hfsifrazil
!!
!! * More heat flux cross-checks
!!   * hfds     = net_heat_coupler + hfsifrazil + heat_pme
!!   * heat_pme = heat_content_surfwater
!!              = heat_content_massin + heat_content_massout
!!              = heat_content_fprec + heat_content_cond + heat_content_vprec
!!               + hfrunoffds + hfevapds + hfrainds
!!
!!  \subsection subsection_3d_heat_budget Depth integrated heat budget
!!
!!  Here is an example 3d heat budget diagnostic for MOM.
!!
!! * OPOTTEMPTEND = T_ADVECTION_XY + TH_TENDENCY_VERT_REMAP + OPOTTEMPDIFF + OPOTTEMPPMDIFF
!!                + BOUNDARY_FORCING_HEAT_TENDENCY + FRAZIL_HEAT_TENDENCY
!!
!! * OPOTTEMPTEND                   = net tendency of heat as diagnosed in MOM.F90
!! * T_ADVECTION_XY                 = heating of a cell from lateral advection
!! * TH_TENDENCY_VERT_REMAP         = heating of a cell from vertical remapping
!! * OPOTTEMPDIFF                   = heating of a cell from diabatic diffusion
!! * OPOTTEMPPMDIFF                 = heating of a cell from neutral diffusion
!! * BOUNDARY_FORCING_HEAT_TENDENCY = heating of cell from boundary fluxes
!! * FRAZIL_HEAT_TENDENCY           = heating of cell from frazil
!!
!! * TH_TENDENCY_VERT_REMAP has zero vertical sum, as it redistributes heat in vertical.
!!
!! * OPOTTEMPDIFF has zero vertical sum, as it redistributes heat in the vertical.
!!
!! * BOUNDARY_FORCING_HEAT_TENDENCY generally has 3d structure, with k > 1 contributions from
!!   penetrative shortwave, and from other fluxes for the case when layers are tiny, in which
!!   case MOM6 partitions tendencies into k > 1 layers.
!!
!! * FRAZIL_HEAT_TENDENCY generally has 3d structure, since MOM6 frazil calculation checks the
!!   full ocean column.
!!
!! * FRAZIL_HEAT_TENDENCY[k=\@sum] = HFSIFRAZIL = column integrated frazil heating.
!!
!! * HFDS = FRAZIL_HEAT_TENDENCY[k=\@sum] + BOUNDARY_FORCING_HEAT_TENDENCY[k=\@sum]
!!
!!  Here is an example 2d heat budget (depth summed) diagnostic for MOM.
!!
!! * OPOTTEMPTEND_2d = T_ADVECTION_XY_2d + OPOTTEMPPMDIFF_2d + HFDS
!!
!!
!!  Here is an example 3d salt budget diagnostic for MOM.
!!
!! * OSALTTEND = S_ADVECTION_XY + SH_TENDENCY_VERT_REMAP + OSALTDIFF + OSALTPMDIFF
!!                + BOUNDARY_FORCING_SALT_TENDENCY
!!
!! * OSALTTEND                      = net tendency of salt as diagnosed in MOM.F90
!! * S_ADVECTION_XY                 = salt convergence to cell from lateral advection
!! * SH_TENDENCY_VERT_REMAP         = salt convergence to cell from vertical remapping
!! * OSALTDIFF                      = salt convergence to cell from diabatic diffusion
!! * OSALTPMDIFF                    = salt convergence to cell from neutral diffusion
!! * BOUNDARY_FORCING_SALT_TENDENCY = salt convergence to cell from boundary fluxes
!!
!! * SH_TENDENCY_VERT_REMAP has zero vertical sum, as it redistributes salt in vertical.
!!
!! * OSALTDIFF has zero vertical sum, as it redistributes salt in the vertical.
!!
!! * BOUNDARY_FORCING_SALT_TENDENCY generally has 3d structure, with k > 1 contributions from
!!   the case when layers are tiny, in which case MOM6 partitions tendencies into k > 1 layers.
!!
!! * SFDSI = BOUNDARY_FORCING_SALT_TENDENCY[k=\@sum]
!!
!!  Here is an example 2d salt budget (depth summed) diagnostic for MOM.
!!
!! * OSALTTEND_2d = S_ADVECTION_XY_2d + OSALTPMDIFF_2d + SFDSI (+ SALT_FLUX_RESTORE)
!!
!!
!!
end module MOM
