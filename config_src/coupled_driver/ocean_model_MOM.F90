!> Top-level module for the MOM6 ocean model in coupled mode.
module ocean_model_mod

! This file is part of MOM6. See LICENSE.md for the license.

! This is the top level module for the MOM6 ocean model.  It contains routines
! for initialization, termination and update of ocean model state.  This
! particular version wraps all of the calls for MOM6 in the calls that had
! been used for MOM4.
!
! This code is a stop-gap wrapper of the MOM6 code to enable it to be called
! in the same way as MOM4.

use MOM, only : initialize_MOM, step_MOM, MOM_control_struct, MOM_end
use MOM, only : extract_surface_state, allocate_surface_state, finish_MOM_initialization
use MOM, only : get_MOM_state_elements, MOM_state_is_synchronized
use MOM, only : get_ocean_stocks, step_offline
use MOM_constants, only : CELSIUS_KELVIN_OFFSET, hlf
use MOM_diag_mediator, only : diag_ctrl, enable_averaging, disable_averaging
use MOM_diag_mediator, only : diag_mediator_close_registration, diag_mediator_end
use MOM_domains, only : pass_var, pass_vector, AGRID, BGRID_NE, CGRID_NE
use MOM_domains, only : TO_ALL, Omit_Corners
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave
use MOM_file_parser, only : get_param, log_version, close_param_file, param_file_type
use MOM_forcing_type, only : forcing, mech_forcing, allocate_forcing_type
use MOM_forcing_type, only : fluxes_accumulate, get_net_mass_forcing
use MOM_forcing_type, only : copy_back_forcing_fields
use MOM_forcing_type, only : forcing_diagnostics, mech_forcing_diags
use MOM_get_input, only : Get_MOM_Input, directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : close_file, file_exists, read_data, write_version_number
use MOM_marine_ice, only : iceberg_forces, iceberg_fluxes, marine_ice_init, marine_ice_CS
use MOM_restart, only : MOM_restart_CS, save_restart
use MOM_string_functions, only : uppercase
use MOM_surface_forcing_gfdl, only : surface_forcing_init, convert_IOB_to_fluxes
use MOM_surface_forcing_gfdl, only : convert_IOB_to_forces, ice_ocn_bnd_type_chksum
use MOM_surface_forcing_gfdl, only : ice_ocean_boundary_type, surface_forcing_CS
use MOM_surface_forcing_gfdl, only : forcing_save_restart
use MOM_time_manager, only : time_type, operator(>), operator(+), operator(-)
use MOM_time_manager, only : operator(*), operator(/), operator(/=)
use MOM_time_manager, only : operator(<=), operator(>=), operator(<)
use MOM_time_manager, only : real_to_time, time_type_to_real
use MOM_tracer_flow_control, only : call_tracer_register, tracer_flow_control_init
use MOM_tracer_flow_control, only : call_tracer_flux_init
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface
use MOM_verticalGrid, only : verticalGrid_type
use MOM_ice_shelf, only : initialize_ice_shelf, shelf_calc_flux, ice_shelf_CS
use MOM_ice_shelf, only : add_shelf_forces, ice_shelf_end, ice_shelf_save_restart
use coupler_types_mod, only : coupler_1d_bc_type, coupler_2d_bc_type
use coupler_types_mod, only : coupler_type_spawn, coupler_type_write_chksums
use coupler_types_mod, only : coupler_type_initialized, coupler_type_copy_data
use coupler_types_mod, only : coupler_type_set_diags, coupler_type_send_data
use mpp_domains_mod, only : domain2d, mpp_get_layout, mpp_get_global_domain
use mpp_domains_mod, only : mpp_define_domains, mpp_get_compute_domain, mpp_get_data_domain
use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux
use fms_mod, only : stdout
use mpp_mod, only : mpp_chksum
use MOM_EOS, only : gsw_sp_from_sr, gsw_pt_from_ct
use MOM_wave_interface, only: wave_parameters_CS, MOM_wave_interface_init
use MOM_wave_interface, only: MOM_wave_interface_init_lite, Update_Surface_Waves

#include <MOM_memory.h>

#ifdef _USE_GENERIC_TRACER
use MOM_generic_tracer, only : MOM_generic_tracer_fluxes_accumulate
#endif

implicit none ; private

public ocean_model_init, ocean_model_end, update_ocean_model
public ocean_model_save_restart, Ocean_stock_pe
public ice_ocean_boundary_type
public ocean_model_init_sfc, ocean_model_flux_init
public ocean_model_restart
public ice_ocn_bnd_type_chksum
public ocean_public_type_chksum
public ocean_model_data_get
public get_ocean_grid
public ocean_model_get_UV_surf

!> This interface extracts a named scalar field or array from the ocean surface or public type
interface ocean_model_data_get
  module procedure ocean_model_data1D_get
  module procedure ocean_model_data2D_get
end interface


!> This type is used for communication with other components via the FMS coupler.
!! The element names and types can be changed only with great deliberation, hence
!! the persistnce of things like the cutsy element name "avg_kount".
type, public ::  ocean_public_type
  type(domain2d) :: Domain    !< The domain for the surface fields.
  logical :: is_ocean_pe      !< .true. on processors that run the ocean model.
  character(len=32) :: instance_name = '' !< A name that can be used to identify
                                 !! this instance of an ocean model, for example
                                 !! in ensembles when writing messages.
  integer, pointer, dimension(:) :: pelist => NULL()   !< The list of ocean PEs.
  logical, pointer, dimension(:,:) :: maskmap =>NULL() !< A pointer to an array
                    !! indicating which logical processors are actually used for
                    !! the ocean code. The other logical processors would be all
                    !! land points and are not assigned to actual processors.
                    !! This need not be assigned if all logical processors are used.

  integer :: stagger = -999   !< The staggering relative to the tracer points
                    !! points of the two velocity components. Valid entries
                    !! include AGRID, BGRID_NE, CGRID_NE, BGRID_SW, and CGRID_SW,
                    !! corresponding to the community-standard Arakawa notation.
                    !! (These are named integers taken from mpp_parameter_mod.)
                    !! Following MOM5, stagger is BGRID_NE by default when the
                    !! ocean is initialized, but here it is set to -999 so that
                    !! a global max across ocean and non-ocean processors can be
                    !! used to determine its value.
  real, pointer, dimension(:,:)  :: &
    t_surf => NULL(), & !< SST on t-cell (degrees Kelvin)
    s_surf => NULL(), & !< SSS on t-cell (psu)
    u_surf => NULL(), & !< i-velocity at the locations indicated by stagger [m s-1].
    v_surf => NULL(), & !< j-velocity at the locations indicated by stagger [m s-1].
    sea_lev => NULL(), & !< Sea level in m after correction for surface pressure,
                        !! i.e. dzt(1) + eta_t + patm/rho0/grav [m]
    frazil =>NULL(), &  !< Accumulated heating [J m-2] from frazil
                        !! formation in the ocean.
    area => NULL()      !< cell area of the ocean surface [m2].
  type(coupler_2d_bc_type) :: fields    !< A structure that may contain named
                                        !! arrays of tracer-related surface fields.
  integer                  :: avg_kount !< A count of contributions to running
                                        !! sums, used externally by the FMS coupler
                                        !! for accumulating averages of this type.
  integer, dimension(2)    :: axes = 0  !< Axis numbers that are available
                                        !! for I/O using this surface data.
end type ocean_public_type


!> The ocean_state_type contains all information about the state of the ocean,
!! with a format that is private so it can be readily changed without disrupting
!! other coupled components.
type, public :: ocean_state_type ; private
  ! This type is private, and can therefore vary between different ocean models.
  logical :: is_ocean_PE = .false.  !< True if this is an ocean PE.
  type(time_type) :: Time     !< The ocean model's time and master clock.
  type(time_type) :: Time_dyn !< The ocean model's time for the dynamics.  Time and Time_dyn
                              !! should be the same after a full time step.
  integer :: Restart_control  !< An integer that is bit-tested to determine whether
                              !! incremental restart files are saved and whether they
                              !! have a time stamped name.  +1 (bit 0) for generic
                              !! files and +2 (bit 1) for time-stamped files.  A
                              !! restart file is saved at the end of a run segment
                              !! unless Restart_control is negative.

  integer :: nstep = 0        !< The number of calls to update_ocean that update the dynamics.
  integer :: nstep_thermo = 0 !< The number of calls to update_ocean that update the thermodynamics.
  logical :: use_ice_shelf    !< If true, the ice shelf model is enabled.
  logical :: use_waves        !< If true use wave coupling.

  logical :: icebergs_alter_ocean !< If true, the icebergs can change ocean the
                              !! ocean dynamics and forcing fluxes.
  real :: press_to_z          !< A conversion factor between pressure and ocean
                              !! depth in m, usually 1/(rho_0*g) [m Pa-1].
  real :: C_p                 !< The heat capacity of seawater [J degC-1 kg-1].
  logical :: offline_tracer_mode = .false. !< If false, use the model in prognostic mode
                              !! with the barotropic and baroclinic dynamics, thermodynamics,
                              !! etc. stepped forward integrated in time.
                              !! If true, all of the above are bypassed with all
                              !! fields necessary to integrate only the tracer advection
                              !! and diffusion equation read in from files stored from
                              !! a previous integration of the prognostic model.

  logical :: single_step_call !< If true, advance the state of MOM with a single
                              !! step including both dynamics and thermodynamics.
                              !! If false, the two phases are advanced with
                              !! separate calls. The default is true.
  ! The following 3 variables are only used here if single_step_call is false.
  real    :: dt               !< (baroclinic) dynamics time step [s]
  real    :: dt_therm         !< thermodynamics time step [s]
  logical :: thermo_spans_coupling !< If true, thermodynamic and tracer time
                              !! steps can span multiple coupled time steps.
  logical :: diabatic_first   !< If true, apply diabatic and thermodynamic
                              !! processes before time stepping the dynamics.

  type(directories) :: dirs   !< A structure containing several relevant directory paths.
  type(mech_forcing) :: forces !< A structure with the driving mechanical surface forces
  type(forcing)   :: fluxes   !< A structure containing pointers to
                              !! the thermodynamic ocean forcing fields.
  type(forcing)   :: flux_tmp !< A secondary structure containing pointers to the
                              !! ocean forcing fields for when multiple coupled
                              !! timesteps are taken per thermodynamic step.
  type(surface)   :: sfc_state !< A structure containing pointers to
                              !! the ocean surface state fields.
  type(ocean_grid_type), pointer :: &
    grid => NULL()            !< A pointer to a grid structure containing metrics
                              !! and related information.
  type(verticalGrid_type), pointer :: &
    GV => NULL()              !< A pointer to a structure containing information
                              !! about the vertical grid.
  type(unit_scale_type), pointer :: &
    US => NULL()              !< A pointer to a structure containing dimensional
                              !! unit scaling factors.
  type(MOM_control_struct), pointer :: &
    MOM_CSp => NULL()         !< A pointer to the MOM control structure
  type(ice_shelf_CS), pointer :: &
    Ice_shelf_CSp => NULL()   !< A pointer to the control structure for the
                              !! ice shelf model that couples with MOM6.  This
                              !! is null if there is no ice shelf.
  type(marine_ice_CS), pointer :: &
    marine_ice_CSp => NULL()  !< A pointer to the control structure for the
                              !! marine ice effects module.
  type(wave_parameters_cs), pointer :: &
    Waves !< A structure containing pointers to the surface wave fields
  type(surface_forcing_CS), pointer :: &
    forcing_CSp => NULL()     !< A pointer to the MOM forcing control structure
  type(MOM_restart_CS), pointer :: &
    restart_CSp => NULL()     !< A pointer set to the restart control structure
                              !! that will be used for MOM restart files.
  type(diag_ctrl), pointer :: &
    diag => NULL()            !< A pointer to the diagnostic regulatory structure
end type ocean_state_type

contains

!> ocean_model_init initializes the ocean model, including registering fields
!! for restarts and reading restart files if appropriate.
!!
!!   This subroutine initializes both the ocean state and the ocean surface type.
!! Because of the way that indicies and domains are handled, Ocean_sfc must have
!! been used in a previous call to initialize_ocean_type.
subroutine ocean_model_init(Ocean_sfc, OS, Time_init, Time_in, wind_stagger, gas_fields_ocn)
  type(ocean_public_type), target, &
                       intent(inout) :: Ocean_sfc !< A structure containing various publicly
                                !! visible ocean surface properties after initialization,
                                !! the data in this type is intent out.
  type(ocean_state_type), pointer    :: OS        !< A structure whose internal
                                !! contents are private to ocean_model_mod that may be used to
                                !! contain all information about the ocean's interior state.
  type(time_type),     intent(in)    :: Time_init !< The start time for the coupled model's calendar
  type(time_type),     intent(in)    :: Time_in   !< The time at which to initialize the ocean model.
  integer, optional,   intent(in)    :: wind_stagger !< If present, the staggering of the winds that are
                                                     !! being provided in calls to update_ocean_model
  type(coupler_1d_bc_type), &
             optional, intent(in)    :: gas_fields_ocn !< If present, this type describes the
                                              !! ocean and surface-ice fields that will participate
                                              !! in the calculation of additional gas or other
                                              !! tracer fluxes, and can be used to spawn related
                                              !! internal variables in the ice model.
  ! Local variables
  real :: Rho0        ! The Boussinesq ocean density [kg m-3].
  real :: G_Earth     ! The gravitational acceleration [m s-2].
  real :: HFrz        !< If HFrz > 0 (m), melt potential will be computed.
                      !! The actual depth over which melt potential is computed will
                      !! min(HFrz, OBLD), where OBLD is the boundary layer depth.
                      !! If HFrz <= 0 (default), melt potential will not be computed.
  logical :: use_melt_pot!< If true, allocate melt_potential array

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "ocean_model_init"  ! This module's name.
  character(len=48)  :: stagger ! A string indicating the staggering locations for the
                                ! surface velocities returned to the coupler.
  type(param_file_type) :: param_file !< A structure to parse for run-time parameters
  logical :: use_temperature ! If true, temperature and salinity are state variables.

  call callTree_enter("ocean_model_init(), ocean_model_MOM.F90")
  if (associated(OS)) then
    call MOM_error(WARNING, "ocean_model_init called with an associated "// &
                   "ocean_state_type structure. Model is already initialized.")
    return
  endif
  allocate(OS)

  OS%is_ocean_pe = Ocean_sfc%is_ocean_pe
  if (.not.OS%is_ocean_pe) return

  OS%Time = Time_in ; OS%Time_dyn = Time_in
  call initialize_MOM(OS%Time, Time_init, param_file, OS%dirs, OS%MOM_CSp, &
                      OS%restart_CSp, Time_in, offline_tracer_mode=OS%offline_tracer_mode, &
                      diag_ptr=OS%diag, count_calls=.true.)
  call get_MOM_state_elements(OS%MOM_CSp, G=OS%grid, GV=OS%GV, US=OS%US, C_p=OS%C_p, &
                              C_p_scaled=OS%fluxes%C_p, use_temp=use_temperature)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "SINGLE_STEPPING_CALL", OS%single_step_call, &
                 "If true, advance the state of MOM with a single step "//&
                 "including both dynamics and thermodynamics.  If false, "//&
                 "the two phases are advanced with separate calls.", default=.true.)
  call get_param(param_file, mdl, "DT", OS%dt, &
                 "The (baroclinic) dynamics time step.  The time-step that "//&
                 "is actually used will be an integer fraction of the "//&
                 "forcing time-step.", units="s", fail_if_missing=.true.)
  call get_param(param_file, mdl, "DT_THERM", OS%dt_therm, &
                 "The thermodynamic and tracer advection time step. "//&
                 "Ideally DT_THERM should be an integer multiple of DT "//&
                 "and less than the forcing or coupling time-step, unless "//&
                 "THERMO_SPANS_COUPLING is true, in which case DT_THERM "//&
                 "can be an integer multiple of the coupling timestep.  By "//&
                 "default DT_THERM is set to DT.", units="s", default=OS%dt)
  call get_param(param_file, "MOM", "THERMO_SPANS_COUPLING", OS%thermo_spans_coupling, &
                 "If true, the MOM will take thermodynamic and tracer "//&
                 "timesteps that can be longer than the coupling timestep. "//&
                 "The actual thermodynamic timestep that is used in this "//&
                 "case is the largest integer multiple of the coupling "//&
                 "timestep that is less than or equal to DT_THERM.", default=.false.)
  call get_param(param_file, mdl, "DIABATIC_FIRST", OS%diabatic_first, &
                 "If true, apply diabatic and thermodynamic processes, "//&
                 "including buoyancy forcing and mass gain or loss, "//&
                 "before stepping the dynamics forward.", default=.false.)

  call get_param(param_file, mdl, "RESTART_CONTROL", OS%Restart_control, &
                 "An integer whose bits encode which restart files are "//&
                 "written. Add 2 (bit 1) for a time-stamped file, and odd "//&
                 "(bit 0) for a non-time-stamped file.  A restart file "//&
                 "will be saved at the end of the run segment for any "//&
                 "non-negative value.", default=1)
  call get_param(param_file, mdl, "OCEAN_SURFACE_STAGGER", stagger, &
                 "A case-insensitive character string to indicate the "//&
                 "staggering of the surface velocity field that is "//&
                 "returned to the coupler.  Valid values include "//&
                 "'A', 'B', or 'C'.", default="C")
  if (uppercase(stagger(1:1)) == 'A') then ; Ocean_sfc%stagger = AGRID
  elseif (uppercase(stagger(1:1)) == 'B') then ; Ocean_sfc%stagger = BGRID_NE
  elseif (uppercase(stagger(1:1)) == 'C') then ; Ocean_sfc%stagger = CGRID_NE
  else ; call MOM_error(FATAL,"ocean_model_init: OCEAN_SURFACE_STAGGER = "// &
                        trim(stagger)//" is invalid.") ; endif

  call get_param(param_file, mdl, "RHO_0", Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINSEQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0)
  call get_param(param_file, mdl, "G_EARTH", G_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80)

  call get_param(param_file, mdl, "ICE_SHELF",  OS%use_ice_shelf, &
                 "If true, enables the ice shelf model.", default=.false.)

  call get_param(param_file, mdl, "ICEBERGS_APPLY_RIGID_BOUNDARY",  OS%icebergs_alter_ocean, &
                 "If true, allows icebergs to change boundary condition felt by ocean", default=.false.)

  OS%press_to_z = 1.0/(Rho0*G_Earth)

  !   Consider using a run-time flag to determine whether to do the diagnostic
  ! vertical integrals, since the related 3-d sums are not negligible in cost.
  call get_param(param_file, mdl, "HFREEZE", HFrz, &
                 "If HFREEZE > 0, melt potential will be computed. The actual depth "//&
                 "over which melt potential is computed will be min(HFREEZE, OBLD), "//&
                 "where OBLD is the boundary layer depth. If HFREEZE <= 0 (default), "//&
                 "melt potential will not be computed.", units="m", default=-1.0, do_not_log=.true.)

  if (HFrz .gt. 0.0) then
    use_melt_pot=.true.
  else
    use_melt_pot=.false.
  endif

  call allocate_surface_state(OS%sfc_state, OS%grid, use_temperature, do_integrals=.true., &
                              gas_fields_ocn=gas_fields_ocn, use_meltpot=use_melt_pot)

  if (present(wind_stagger)) then
    call surface_forcing_init(Time_in, OS%grid, OS%US, param_file, OS%diag, &
                              OS%forcing_CSp, wind_stagger)
  else
    call surface_forcing_init(Time_in, OS%grid, OS%US, param_file, OS%diag, &
                              OS%forcing_CSp)
  endif

  if (OS%use_ice_shelf)  then
    call initialize_ice_shelf(param_file, OS%grid, OS%Time, OS%ice_shelf_CSp, &
                              OS%diag, OS%forces, OS%fluxes)
  endif
  if (OS%icebergs_alter_ocean)  then
    call marine_ice_init(OS%Time, OS%grid, param_file, OS%diag, OS%marine_ice_CSp)
    if (.not. OS%use_ice_shelf) &
      call allocate_forcing_type(OS%grid, OS%fluxes, shelf=.true.)
  endif

  call get_param(param_file, mdl, "USE_WAVES", OS%Use_Waves, &
       "If true, enables surface wave modules.", default=.false.)
  if (OS%use_waves) then
    call MOM_wave_interface_init(OS%Time, OS%grid, OS%GV, OS%US, param_file, OS%Waves, OS%diag)
  else
    call MOM_wave_interface_init_lite(param_file)
  endif

  if (associated(OS%grid%Domain%maskmap)) then
    call initialize_ocean_public_type(OS%grid%Domain%mpp_domain, Ocean_sfc, &
                                      OS%diag, maskmap=OS%grid%Domain%maskmap, &
                                      gas_fields_ocn=gas_fields_ocn)
  else
    call initialize_ocean_public_type(OS%grid%Domain%mpp_domain, Ocean_sfc, &
                                      OS%diag, gas_fields_ocn=gas_fields_ocn)
  endif

  ! This call can only occur here if the coupler_bc_type variables have been
  ! initialized already using the information from gas_fields_ocn.
  if (present(gas_fields_ocn)) then
    call coupler_type_set_diags(Ocean_sfc%fields, "ocean_sfc", &
                                Ocean_sfc%axes(1:2), Time_in)

    call extract_surface_state(OS%MOM_CSp, OS%sfc_state)

    call convert_state_to_ocean_type(OS%sfc_state, Ocean_sfc, OS%grid, OS%US)

  endif

  call close_param_file(param_file)
  call diag_mediator_close_registration(OS%diag)

  call MOM_mesg('==== Completed MOM6 Coupled Initialization ====', 2)

  call callTree_leave("ocean_model_init(")
end subroutine ocean_model_init

!> update_ocean_model uses the forcing in Ice_ocean_boundary to advance the
!! ocean model's state from the input value of Ocean_state (which must be for
!! time time_start_update) for a time interval of Ocean_coupling_time_step,
!! returning the publicly visible ocean surface properties in Ocean_sfc and
!! storing the new ocean properties in Ocean_state.
subroutine update_ocean_model(Ice_ocean_boundary, OS, Ocean_sfc, time_start_update, &
                              Ocean_coupling_time_step, update_dyn, update_thermo, &
                              Ocn_fluxes_used, start_cycle, end_cycle, cycle_length)
  type(ice_ocean_boundary_type), &
                     intent(in)    :: Ice_ocean_boundary !< A structure containing the various
                                              !! forcing fields coming from the ice and atmosphere.
  type(ocean_state_type), &
                     pointer       :: OS      !< A pointer to a private structure containing the
                                              !! internal ocean state.
  type(ocean_public_type), &
                     intent(inout) :: Ocean_sfc !< A structure containing all the publicly visible
                                              !! ocean surface fields after a coupling time step.
                                              !! The data in this type is intent out.
  type(time_type),   intent(in)    :: time_start_update  !< The time at the beginning of the update step.
  type(time_type),   intent(in)    :: Ocean_coupling_time_step !< The amount of time over which to
                                              !! advance the ocean.
  logical, optional, intent(in)    :: update_dyn !< If present and false, do not do updates
                                              !! due to the ocean dynamics.
  logical, optional, intent(in)    :: update_thermo !< If present and false, do not do updates
                                              !! due to the ocean thermodynamics or remapping.
  logical, optional, intent(in)    :: Ocn_fluxes_used !< If present, this indicates whether the
                                              !! cumulative thermodynamic fluxes from the ocean,
                                              !! like frazil, have been used and should be reset.
  logical, optional, intent(in)    :: start_cycle !< This indicates whether this call is to be
                                              !! treated as the first call to step_MOM in a
                                              !! time-stepping cycle; missing is like true.
  logical, optional, intent(in)    :: end_cycle   !< This indicates whether this call is to be
                                              !! treated as the last call to step_MOM in a
                                              !! time-stepping cycle; missing is like true.
  real,    optional, intent(in)    :: cycle_length !< The duration of a coupled time stepping cycle [s].

  ! Local variables
  type(time_type) :: Time_seg_start ! Stores the dynamic or thermodynamic ocean model time at the
                            ! start of this call to allow step_MOM to temporarily change the time
                            ! as seen by internal modules.
  type(time_type) :: Time_thermo_start ! Stores the ocean model thermodynamics time at the start of
                            ! this call to allow step_MOM to temporarily change the time as seen by
                            ! internal modules.
  type(time_type) :: Time1  ! The value of the ocean model's time at the start of a call to step_MOM.
  integer :: index_bnds(4)  ! The computational domain index bounds in the ice-ocean boundary type.
  real :: weight            ! Flux accumulation weight of the current fluxes.
  real :: dt_coupling       ! The coupling time step [s].
  real :: dt_therm          ! A limited and quantized version of OS%dt_therm [s].
  real :: dt_dyn            ! The dynamics time step [s].
  real :: dtdia             ! The diabatic time step [s].
  real :: t_elapsed_seg     ! The elapsed time in this update segment [s].
  integer :: n              ! The internal iteration counter.
  integer :: nts            ! The number of baroclinic dynamics time steps in a thermodynamic step.
  integer :: n_max          ! The number of calls to step_MOM dynamics in this call to update_ocean_model.
  integer :: n_last_thermo  ! The iteration number the last time thermodynamics were updated.
  logical :: thermo_does_span_coupling ! If true, thermodynamic forcing spans multiple dynamic timesteps.
  logical :: do_dyn         ! If true, step the ocean dynamics and transport.
  logical :: do_thermo      ! If true, step the ocean thermodynamics.
  logical :: step_thermo    ! If true, take a thermodynamic step.
  integer :: is, ie, js, je

  call callTree_enter("update_ocean_model(), ocean_model_MOM.F90")
  dt_coupling = time_type_to_real(Ocean_coupling_time_step)

  if (.not.associated(OS)) then
    call MOM_error(FATAL, "update_ocean_model called with an unassociated "// &
                    "ocean_state_type structure. ocean_model_init must be "//  &
                    "called first to allocate this structure.")
    return
  endif

  do_dyn = .true. ; if (present(update_dyn)) do_dyn = update_dyn
  do_thermo = .true. ; if (present(update_thermo)) do_thermo = update_thermo

  if (do_thermo .and. (time_start_update /= OS%Time)) &
    call MOM_error(WARNING, "update_ocean_model: internal clock does not "//&
                            "agree with time_start_update argument.")
  if (do_dyn .and. (time_start_update /= OS%Time_dyn)) &
    call MOM_error(WARNING, "update_ocean_model: internal dynamics clock does not "//&
                            "agree with time_start_update argument.")

  if (.not.(do_dyn .or. do_thermo)) call MOM_error(FATAL, &
      "update_ocean_model called without updating either dynamics or thermodynamics.")
  if (do_dyn .and. do_thermo .and. (OS%Time /= OS%Time_dyn)) call MOM_error(FATAL, &
      "update_ocean_model called to update both dynamics and thermodynamics with inconsistent clocks.")

  ! This is benign but not necessary if ocean_model_init_sfc was called or if
  ! OS%sfc_state%tr_fields was spawned in ocean_model_init.  Consider removing it.
  is = OS%grid%isc ; ie = OS%grid%iec ; js = OS%grid%jsc ; je = OS%grid%jec
  call coupler_type_spawn(Ocean_sfc%fields, OS%sfc_state%tr_fields, &
                          (/is,is,ie,ie/), (/js,js,je,je/), as_needed=.true.)

  ! Translate Ice_ocean_boundary into fluxes and forces.
  call mpp_get_compute_domain(Ocean_sfc%Domain, index_bnds(1), index_bnds(2), &
                              index_bnds(3), index_bnds(4))

  if (do_dyn) then
    call convert_IOB_to_forces(Ice_ocean_boundary, OS%forces, index_bnds, OS%Time_dyn, OS%grid, OS%US, &
                               OS%forcing_CSp, dt_forcing=dt_coupling, reset_avg=OS%fluxes%fluxes_used)
    if (OS%use_ice_shelf) &
      call add_shelf_forces(OS%grid, OS%US, OS%Ice_shelf_CSp, OS%forces)
    if (OS%icebergs_alter_ocean) &
      call iceberg_forces(OS%grid, OS%forces, OS%use_ice_shelf, &
                          OS%sfc_state, dt_coupling, OS%marine_ice_CSp)
  endif

  if (do_thermo) then
    if (OS%fluxes%fluxes_used) then
      call convert_IOB_to_fluxes(Ice_ocean_boundary, OS%fluxes, index_bnds, OS%Time, dt_coupling, &
                                 OS%grid, OS%US, OS%forcing_CSp, OS%sfc_state)

      ! Add ice shelf fluxes
      if (OS%use_ice_shelf) &
        call shelf_calc_flux(OS%sfc_state, OS%fluxes, OS%Time, dt_coupling, OS%Ice_shelf_CSp)
      if (OS%icebergs_alter_ocean) &
        call iceberg_fluxes(OS%grid, OS%US, OS%fluxes, OS%use_ice_shelf, &
                            OS%sfc_state, dt_coupling, OS%marine_ice_CSp)

#ifdef _USE_GENERIC_TRACER
      call enable_averaging(dt_coupling, OS%Time + Ocean_coupling_time_step, OS%diag) !Is this needed?
      call MOM_generic_tracer_fluxes_accumulate(OS%fluxes, 1.0) ! Here weight=1, so just store the current fluxes
      call disable_averaging(OS%diag)
#endif
    else
      ! The previous fluxes have not been used yet, so translate the input fluxes
      ! into a temporary type and then accumulate them in about 20 lines.
      OS%flux_tmp%C_p = OS%fluxes%C_p
      call convert_IOB_to_fluxes(Ice_ocean_boundary, OS%flux_tmp, index_bnds, OS%Time, dt_coupling, &
                                 OS%grid, OS%US, OS%forcing_CSp, OS%sfc_state)

      if (OS%use_ice_shelf) &
        call shelf_calc_flux(OS%sfc_state, OS%flux_tmp, OS%Time, dt_coupling, OS%Ice_shelf_CSp)
      if (OS%icebergs_alter_ocean) &
        call iceberg_fluxes(OS%grid, OS%US, OS%flux_tmp, OS%use_ice_shelf, &
                            OS%sfc_state, dt_coupling, OS%marine_ice_CSp)

      call fluxes_accumulate(OS%flux_tmp, OS%fluxes, OS%grid, weight)
#ifdef _USE_GENERIC_TRACER
       ! Incorporate the current tracer fluxes into the running averages
      call MOM_generic_tracer_fluxes_accumulate(OS%flux_tmp, weight)
#endif
    endif
  endif

  ! The net mass forcing is not currently used in the MOM6 dynamics solvers, so this is may be unnecessary.
  if (do_dyn .and. associated(OS%forces%net_mass_src) .and. .not.OS%forces%net_mass_src_set) &
    call get_net_mass_forcing(OS%fluxes, OS%grid, OS%US, OS%forces%net_mass_src)

  if (OS%use_waves .and. do_thermo) then
    ! For now, the waves are only updated on the thermodynamics steps, because that is where
    ! the wave intensities are actually used to drive mixing.  At some point, the wave updates
    ! might also need to become a part of the ocean dynamics, according to B. Reichl.
    call Update_Surface_Waves(OS%grid, OS%GV, OS%US, OS%time, ocean_coupling_time_step, OS%waves)
  endif

  if ((OS%nstep==0) .and. (OS%nstep_thermo==0)) then ! This is the first call to update_ocean_model.
    call finish_MOM_initialization(OS%Time, OS%dirs, OS%MOM_CSp, OS%restart_CSp)
  endif

  Time_thermo_start = OS%Time
  Time_seg_start = OS%Time ; if (do_dyn) Time_seg_start = OS%Time_dyn
  Time1 = Time_seg_start

  if (OS%offline_tracer_mode .and. do_thermo) then
    call step_offline(OS%forces, OS%fluxes, OS%sfc_state, Time1, dt_coupling, OS%MOM_CSp)
  elseif ((.not.do_thermo) .or. (.not.do_dyn)) then
    ! The call sequence is being orchestrated from outside of update_ocean_model.
    call step_MOM(OS%forces, OS%fluxes, OS%sfc_state, Time1, dt_coupling, OS%MOM_CSp, &
                  Waves=OS%Waves, do_dynamics=do_dyn, do_thermodynamics=do_thermo, &
                  start_cycle=start_cycle, end_cycle=end_cycle, cycle_length=cycle_length, &
                  reset_therm=Ocn_fluxes_used)
  elseif (OS%single_step_call) then
    call step_MOM(OS%forces, OS%fluxes, OS%sfc_state, Time1, dt_coupling, OS%MOM_CSp, Waves=OS%Waves)
  else  ! Step both the dynamics and thermodynamics with separate calls.
    n_max = 1 ; if (dt_coupling > OS%dt) n_max = ceiling(dt_coupling/OS%dt - 0.001)
    dt_dyn = dt_coupling / real(n_max)
    thermo_does_span_coupling = (OS%thermo_spans_coupling .and. (OS%dt_therm > 1.5*dt_coupling))

    if (thermo_does_span_coupling) then
      dt_therm = dt_coupling * floor(OS%dt_therm / dt_coupling + 0.001)
      nts = floor(dt_therm/dt_dyn + 0.001)
    else
      nts = MAX(1,MIN(n_max,floor(OS%dt_therm/dt_dyn + 0.001)))
      n_last_thermo = 0
    endif

    Time1 = Time_seg_start ; t_elapsed_seg = 0.0
    do n=1,n_max
      if (OS%diabatic_first) then
        if (thermo_does_span_coupling) call MOM_error(FATAL, &
            "MOM is not yet set up to have restarts that work with "//&
            "THERMO_SPANS_COUPLING and DIABATIC_FIRST.")
        if (modulo(n-1,nts)==0) then
          dtdia = dt_dyn*min(nts,n_max-(n-1))
          call step_MOM(OS%forces, OS%fluxes, OS%sfc_state, Time1, dtdia, OS%MOM_CSp, &
                        Waves=OS%Waves, do_dynamics=.false., do_thermodynamics=.true., &
                        start_cycle=(n==1), end_cycle=.false., cycle_length=dt_coupling)
        endif

        call step_MOM(OS%forces, OS%fluxes, OS%sfc_state, Time1, dt_dyn, OS%MOM_CSp, &
                      Waves=OS%Waves, do_dynamics=.true., do_thermodynamics=.false., &
                      start_cycle=.false., end_cycle=(n==n_max), cycle_length=dt_coupling)
      else
        call step_MOM(OS%forces, OS%fluxes, OS%sfc_state, Time1, dt_dyn, OS%MOM_CSp, &
                      Waves=OS%Waves, do_dynamics=.true., do_thermodynamics=.false., &
                      start_cycle=(n==1), end_cycle=.false., cycle_length=dt_coupling)

        step_thermo = .false.
        if (thermo_does_span_coupling) then
          dtdia = dt_therm
          step_thermo = MOM_state_is_synchronized(OS%MOM_CSp, adv_dyn=.true.)
        elseif ((modulo(n,nts)==0) .or. (n==n_max)) then
          dtdia = dt_dyn*(n - n_last_thermo)
          n_last_thermo = n
          step_thermo = .true.
        endif

        if (step_thermo) then
          ! Back up Time1 to the start of the thermodynamic segment.
          Time1 = Time1 - real_to_time(dtdia - dt_dyn)
          call step_MOM(OS%forces, OS%fluxes, OS%sfc_state, Time1, dtdia, OS%MOM_CSp, &
                        Waves=OS%Waves, do_dynamics=.false., do_thermodynamics=.true., &
                        start_cycle=.false., end_cycle=(n==n_max), cycle_length=dt_coupling)
        endif
      endif

      t_elapsed_seg = t_elapsed_seg + dt_dyn
      Time1 = Time_seg_start + real_to_time(t_elapsed_seg)
    enddo
  endif

  if (do_dyn) OS%Time_dyn = Time_seg_start + Ocean_coupling_time_step
  if (do_dyn) OS%nstep = OS%nstep + 1
  OS%Time = Time_thermo_start  ! Reset the clock to compensate for shared pointers.
  if (do_thermo) OS%Time = OS%Time + Ocean_coupling_time_step
  if (do_thermo) OS%nstep_thermo = OS%nstep_thermo + 1

  if (do_dyn) then
    call mech_forcing_diags(OS%forces, dt_coupling, OS%grid, OS%Time_dyn, OS%diag, OS%forcing_CSp%handles)
  endif

  if (OS%fluxes%fluxes_used .and. do_thermo) then
    call forcing_diagnostics(OS%fluxes, OS%sfc_state, OS%grid, OS%US, OS%Time, OS%diag, OS%forcing_CSp%handles)
  endif

! Translate state into Ocean.
!  call convert_state_to_ocean_type(OS%sfc_state, Ocean_sfc, OS%grid, OS%US, &
!                                   Ice_ocean_boundary%p, OS%press_to_z)
  call convert_state_to_ocean_type(OS%sfc_state, Ocean_sfc, OS%grid, OS%US)
  Time1 = OS%Time ; if (do_dyn) Time1 = OS%Time_dyn
  call coupler_type_send_data(Ocean_sfc%fields, Time1)

  call callTree_leave("update_ocean_model()")
end subroutine update_ocean_model

!> This subroutine writes out the ocean model restart file.
subroutine ocean_model_restart(OS, timestamp)
  type(ocean_state_type),     pointer    :: OS !< A pointer to the structure containing the
                                               !! internal ocean state being saved to a restart file
  character(len=*), optional, intent(in) :: timestamp !< An optional timestamp string that should be
                                               !! prepended to the file name. (Currently this is unused.)

  if (.not.MOM_state_is_synchronized(OS%MOM_CSp)) &
      call MOM_error(WARNING, "End of MOM_main reached with inconsistent "//&
         "dynamics and advective times.  Additional restart fields "//&
         "that have not been coded yet would be required for reproducibility.")
  if (.not.OS%fluxes%fluxes_used) call MOM_error(FATAL, "ocean_model_restart "//&
      "was called with unused buoyancy fluxes.  For conservation, the ocean "//&
      "restart files can only be created after the buoyancy forcing is applied.")

  if (BTEST(OS%Restart_control,1)) then
    call save_restart(OS%dirs%restart_output_dir, OS%Time, OS%grid, &
                      OS%restart_CSp, .true., GV=OS%GV)
    call forcing_save_restart(OS%forcing_CSp, OS%grid, OS%Time, &
                              OS%dirs%restart_output_dir, .true.)
    if (OS%use_ice_shelf) then
      call ice_shelf_save_restart(OS%Ice_shelf_CSp, OS%Time, OS%dirs%restart_output_dir, .true.)
    endif
  endif
  if (BTEST(OS%Restart_control,0)) then
    call save_restart(OS%dirs%restart_output_dir, OS%Time, OS%grid, &
                      OS%restart_CSp, GV=OS%GV)
    call forcing_save_restart(OS%forcing_CSp, OS%grid, OS%Time, &
                              OS%dirs%restart_output_dir)
    if (OS%use_ice_shelf) then
      call ice_shelf_save_restart(OS%Ice_shelf_CSp, OS%Time, OS%dirs%restart_output_dir)
    endif
  endif

end subroutine ocean_model_restart
! </SUBROUTINE> NAME="ocean_model_restart"

!> ocean_model_end terminates the model run, saving the ocean state in a restart
!! and deallocating any data associated with the ocean.
subroutine ocean_model_end(Ocean_sfc, Ocean_state, Time)
  type(ocean_public_type), intent(inout) :: Ocean_sfc   !< An ocean_public_type structure that is
                                                        !! to be deallocated upon termination.
  type(ocean_state_type),  pointer       :: Ocean_state !< A pointer to the structure containing
                                                        !! the internal ocean state to be deallocated
                                                        !! upon termination.
  type(time_type),         intent(in)    :: Time        !< The model time, used for writing restarts.

  call ocean_model_save_restart(Ocean_state, Time)
  call diag_mediator_end(Time, Ocean_state%diag)
  call MOM_end(Ocean_state%MOM_CSp)
  if (Ocean_state%use_ice_shelf) call ice_shelf_end(Ocean_state%Ice_shelf_CSp)
end subroutine ocean_model_end

!> ocean_model_save_restart causes restart files associated with the ocean to be
!! written out.
subroutine ocean_model_save_restart(OS, Time, directory, filename_suffix)
  type(ocean_state_type),     pointer    :: OS  !< A pointer to the structure containing the
                                                !! internal ocean state (in).
  type(time_type),            intent(in) :: Time !< The model time at this call, needed for mpp_write calls.
  character(len=*), optional, intent(in) :: directory  !<  An optional directory into which to
                                                !! write these restart files.
  character(len=*), optional, intent(in) :: filename_suffix !< An optional suffix (e.g., a time-stamp)
                                                !! to append to the restart file names.
! Note: This is a new routine - it will need to exist for the new incremental
!   checkpointing.  It will also be called by ocean_model_end, giving the same
!   restart behavior as now in FMS.
  character(len=200) :: restart_dir

  if (.not.MOM_state_is_synchronized(OS%MOM_CSp)) &
    call MOM_error(WARNING, "ocean_model_save_restart called with inconsistent "//&
         "dynamics and advective times.  Additional restart fields "//&
         "that have not been coded yet would be required for reproducibility.")
  if (.not.OS%fluxes%fluxes_used) call MOM_error(FATAL, "ocean_model_save_restart "//&
       "was called with unused buoyancy fluxes.  For conservation, the ocean "//&
       "restart files can only be created after the buoyancy forcing is applied.")

  if (present(directory)) then ; restart_dir = directory
  else ; restart_dir = OS%dirs%restart_output_dir ; endif

  call save_restart(restart_dir, Time, OS%grid, OS%restart_CSp, GV=OS%GV)

  call forcing_save_restart(OS%forcing_CSp, OS%grid, Time, restart_dir)

  if (OS%use_ice_shelf) then
    call ice_shelf_save_restart(OS%Ice_shelf_CSp, OS%Time, OS%dirs%restart_output_dir)
  endif

end subroutine ocean_model_save_restart

!> Initialize the public ocean type
subroutine initialize_ocean_public_type(input_domain, Ocean_sfc, diag, maskmap, &
                                        gas_fields_ocn)
  type(domain2D),          intent(in)    :: input_domain !< The ocean model domain description
  type(ocean_public_type), intent(inout) :: Ocean_sfc !< A structure containing various publicly
                                !! visible ocean surface properties after initialization, whose
                                !! elements are allocated here.
  type(diag_ctrl),         intent(in)    :: diag  !< A structure that regulates diagnsotic output
  logical, dimension(:,:), &
                 optional, intent(in)    :: maskmap !< A mask indicating which virtual processors
                                              !! are actually in use.  If missing, all are used.
  type(coupler_1d_bc_type), &
                 optional, intent(in)    :: gas_fields_ocn !< If present, this type describes the
                                              !! ocean and surface-ice fields that will participate
                                              !! in the calculation of additional gas or other
                                              !! tracer fluxes.

  integer :: xsz, ysz, layout(2)
  ! ice-ocean-boundary fields are always allocated using absolute indicies
  ! and have no halos.
  integer :: isc, iec, jsc, jec

  call mpp_get_layout(input_domain,layout)
  call mpp_get_global_domain(input_domain, xsize=xsz, ysize=ysz)
  if (PRESENT(maskmap)) then
     call mpp_define_domains((/1,xsz,1,ysz/),layout,Ocean_sfc%Domain, maskmap=maskmap)
  else
     call mpp_define_domains((/1,xsz,1,ysz/),layout,Ocean_sfc%Domain)
  endif
  call mpp_get_compute_domain(Ocean_sfc%Domain, isc, iec, jsc, jec)

  allocate ( Ocean_sfc%t_surf (isc:iec,jsc:jec), &
             Ocean_sfc%s_surf (isc:iec,jsc:jec), &
             Ocean_sfc%u_surf (isc:iec,jsc:jec), &
             Ocean_sfc%v_surf (isc:iec,jsc:jec), &
             Ocean_sfc%sea_lev(isc:iec,jsc:jec), &
             Ocean_sfc%area   (isc:iec,jsc:jec), &
             Ocean_sfc%frazil (isc:iec,jsc:jec))

  Ocean_sfc%t_surf(:,:)  = 0.0  ! time averaged sst (Kelvin) passed to atmosphere/ice model
  Ocean_sfc%s_surf(:,:)  = 0.0  ! time averaged sss (psu) passed to atmosphere/ice models
  Ocean_sfc%u_surf(:,:)  = 0.0  ! time averaged u-current (m/sec) passed to atmosphere/ice models
  Ocean_sfc%v_surf(:,:)  = 0.0  ! time averaged v-current (m/sec)  passed to atmosphere/ice models
  Ocean_sfc%sea_lev(:,:) = 0.0  ! time averaged thickness of top model grid cell (m) plus patm/rho0/grav
  Ocean_sfc%frazil(:,:)  = 0.0  ! time accumulated frazil (J/m^2) passed to ice model
  Ocean_sfc%area(:,:)    = 0.0
  Ocean_sfc%axes    = diag%axesT1%handles !diag axes to be used by coupler tracer flux diagnostics

  if (present(gas_fields_ocn)) then
    call coupler_type_spawn(gas_fields_ocn, Ocean_sfc%fields, (/isc,isc,iec,iec/), &
                              (/jsc,jsc,jec,jec/), suffix = '_ocn', as_needed=.true.)
  endif

end subroutine initialize_ocean_public_type

!> This subroutine translates the coupler's ocean_data_type into MOM's
!! surface state variable.  This may eventually be folded into the MOM
!! code that calculates the surface state in the first place.
!! Note the offset in the arrays because the ocean_data_type has no
!! halo points in its arrays and always uses absolute indicies.
subroutine convert_state_to_ocean_type(sfc_state, Ocean_sfc, G, US, patm, press_to_z)
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                               !! describe the surface state of the ocean.
  type(ocean_public_type), &
                 target, intent(inout) :: Ocean_sfc !< A structure containing various publicly
                                               !! visible ocean surface fields, whose elements
                                               !! have their data set here.
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  real,        optional, intent(in)    :: patm(:,:)  !< The pressure at the ocean surface [Pa].
  real,        optional, intent(in)    :: press_to_z !< A conversion factor between pressure and
                                               !! ocean depth in m, usually 1/(rho_0*g) [m Pa-1].
  ! Local variables
  real :: IgR0
  character(len=48)  :: val_str
  integer :: isc_bnd, iec_bnd, jsc_bnd, jec_bnd
  integer :: i, j, i0, j0, is, ie, js, je

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  call pass_vector(sfc_state%u, sfc_state%v, G%Domain)

  call mpp_get_compute_domain(Ocean_sfc%Domain, isc_bnd, iec_bnd, &
                              jsc_bnd, jec_bnd)
  if (present(patm)) then
    ! Check that the inidicies in patm are (isc_bnd:iec_bnd,jsc_bnd:jec_bnd).
    if (.not.present(press_to_z)) call MOM_error(FATAL, &
        'convert_state_to_ocean_type: press_to_z must be present if patm is.')
  endif

  i0 = is - isc_bnd ; j0 = js - jsc_bnd
  if (sfc_state%T_is_conT) then
    ! Convert the surface T from conservative T to potential T.
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%t_surf(i,j) = gsw_pt_from_ct(sfc_state%SSS(i+i0,j+j0), &
                               sfc_state%SST(i+i0,j+j0)) + CELSIUS_KELVIN_OFFSET
    enddo ; enddo
  else
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%t_surf(i,j) = sfc_state%SST(i+i0,j+j0) + CELSIUS_KELVIN_OFFSET
    enddo ; enddo
  endif
  if (sfc_state%S_is_absS) then
    ! Convert the surface S from absolute salinity to practical salinity.
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%s_surf(i,j) = gsw_sp_from_sr(sfc_state%SSS(i+i0,j+j0))
    enddo ; enddo
  else
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%s_surf(i,j) = sfc_state%SSS(i+i0,j+j0)
    enddo ; enddo
  endif

  if (present(patm)) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%sea_lev(i,j) = US%Z_to_m * sfc_state%sea_lev(i+i0,j+j0) + patm(i,j) * press_to_z
      Ocean_sfc%area(i,j) = US%L_to_m**2 * G%areaT(i+i0,j+j0)
    enddo ; enddo
  else
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%sea_lev(i,j) = US%Z_to_m * sfc_state%sea_lev(i+i0,j+j0)
      Ocean_sfc%area(i,j) = US%L_to_m**2 * G%areaT(i+i0,j+j0)
    enddo ; enddo
  endif

  if (allocated(sfc_state%frazil)) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%frazil(i,j) = US%Q_to_J_kg*US%RZ_to_kg_m2 * sfc_state%frazil(i+i0,j+j0)
    enddo ; enddo
  endif

  if (Ocean_sfc%stagger == AGRID) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%u_surf(i,j) = G%mask2dT(i+i0,j+j0) * US%L_T_to_m_s * &
                0.5*(sfc_state%u(I+i0,j+j0)+sfc_state%u(I-1+i0,j+j0))
      Ocean_sfc%v_surf(i,j) = G%mask2dT(i+i0,j+j0) * US%L_T_to_m_s * &
                0.5*(sfc_state%v(i+i0,J+j0)+sfc_state%v(i+i0,J-1+j0))
    enddo ; enddo
  elseif (Ocean_sfc%stagger == BGRID_NE) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%u_surf(i,j) = G%mask2dBu(I+i0,J+j0) * US%L_T_to_m_s * &
                0.5*(sfc_state%u(I+i0,j+j0)+sfc_state%u(I+i0,j+j0+1))
      Ocean_sfc%v_surf(i,j) = G%mask2dBu(I+i0,J+j0) * US%L_T_to_m_s * &
                0.5*(sfc_state%v(i+i0,J+j0)+sfc_state%v(i+i0+1,J+j0))
    enddo ; enddo
  elseif (Ocean_sfc%stagger == CGRID_NE) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%u_surf(i,j) = G%mask2dCu(I+i0,j+j0)*US%L_T_to_m_s * sfc_state%u(I+i0,j+j0)
      Ocean_sfc%v_surf(i,j) = G%mask2dCv(i+i0,J+j0)*US%L_T_to_m_s * sfc_state%v(i+i0,J+j0)
    enddo ; enddo
  else
    write(val_str, '(I8)') Ocean_sfc%stagger
    call MOM_error(FATAL, "convert_state_to_ocean_type: "//&
      "Ocean_sfc%stagger has the unrecognized value of "//trim(val_str))
  endif

  if (coupler_type_initialized(sfc_state%tr_fields)) then
    if (.not.coupler_type_initialized(Ocean_sfc%fields)) then
      call MOM_error(FATAL, "convert_state_to_ocean_type: "//&
               "Ocean_sfc%fields has not been initialized.")
    endif
    call coupler_type_copy_data(sfc_state%tr_fields, Ocean_sfc%fields)
  endif

end subroutine convert_state_to_ocean_type

!>   This subroutine extracts the surface properties from the ocean's internal
!! state and stores them in the ocean type returned to the calling ice model.
!! It has to be separate from the ocean_initialization call because the coupler
!! module allocates the space for some of these variables.
subroutine ocean_model_init_sfc(OS, Ocean_sfc)
  type(ocean_state_type),  pointer       :: OS  !< The structure with the complete ocean state
  type(ocean_public_type), intent(inout) :: Ocean_sfc !< A structure containing various publicly
                                !! visible ocean surface properties after initialization, whose
                                !! elements have their data set here.
  integer :: is, ie, js, je

  is = OS%grid%isc ; ie = OS%grid%iec ; js = OS%grid%jsc ; je = OS%grid%jec
  call coupler_type_spawn(Ocean_sfc%fields, OS%sfc_state%tr_fields, &
                          (/is,is,ie,ie/), (/js,js,je,je/), as_needed=.true.)

  call extract_surface_state(OS%MOM_CSp, OS%sfc_state)

  call convert_state_to_ocean_type(OS%sfc_state, Ocean_sfc, OS%grid, OS%US)

end subroutine ocean_model_init_sfc

!> ocean_model_flux_init is used to initialize properties of the air-sea fluxes
!! as determined by various run-time parameters.  It can be called from
!! non-ocean PEs, or PEs that have not yet been initialzed, and it can safely
!! be called multiple times.
subroutine ocean_model_flux_init(OS, verbosity)
  type(ocean_state_type), optional, pointer :: OS  !< An optional pointer to the ocean state,
                                             !! used to figure out if this is an ocean PE that
                                             !! has already been initialized.
  integer, optional, intent(in) :: verbosity !< A 0-9 integer indicating a level of verbosity.

  logical :: OS_is_set
  integer :: verbose

  OS_is_set = .false. ; if (present(OS)) OS_is_set = associated(OS)

  ! Use this to control the verbosity of output; consider rethinking this logic later.
  verbose = 5 ; if (OS_is_set) verbose = 3
  if (present(verbosity)) verbose = verbosity

  call call_tracer_flux_init(verbosity=verbose)

end subroutine ocean_model_flux_init

!> Ocean_stock_pe - returns the integrated stocks of heat, water, etc. for conservation checks.
!!   Because of the way FMS is coded, only the root PE has the integrated amount,
!!   while all other PEs get 0.
subroutine Ocean_stock_pe(OS, index, value, time_index)
  use stock_constants_mod, only : ISTOCK_WATER, ISTOCK_HEAT,ISTOCK_SALT
  type(ocean_state_type), pointer     :: OS         !< A structure containing the internal ocean state.
                                                    !! The data in OS is intent in.
  integer,                intent(in)  :: index      !< The stock index for the quantity of interest.
  real,                   intent(out) :: value      !< Sum returned for the conservation quantity of interest.
  integer,      optional, intent(in)  :: time_index !< An unused optional argument, present only for
                                                    !! interfacial compatibility with other models.
! Arguments: OS - A structure containing the internal ocean state.
!  (in)      index - Index of conservation quantity of interest.
!  (in)      value -  Sum returned for the conservation quantity of interest.
!  (in,opt)  time_index - Index for time level to use if this is necessary.

  real :: salt

  value = 0.0
  if (.not.associated(OS)) return
  if (.not.OS%is_ocean_pe) return

  select case (index)
    case (ISTOCK_WATER)  ! Return the mass of fresh water in the ocean in kg.
      if (OS%GV%Boussinesq) then
        call get_ocean_stocks(OS%MOM_CSp, mass=value, on_PE_only=.true.)
      else  ! In non-Boussinesq mode, the mass of salt needs to be subtracted.
        call get_ocean_stocks(OS%MOM_CSp, mass=value, salt=salt, on_PE_only=.true.)
        value = value - salt
      endif
    case (ISTOCK_HEAT)  ! Return the heat content of the ocean in J.
      call get_ocean_stocks(OS%MOM_CSp, heat=value, on_PE_only=.true.)
    case (ISTOCK_SALT)  ! Return the mass of the salt in the ocean in kg.
       call get_ocean_stocks(OS%MOM_CSp, salt=value, on_PE_only=.true.)
    case default ; value = 0.0
  end select
  ! If the FMS coupler is changed so that Ocean_stock_PE is only called on
  ! ocean PEs, uncomment the following and eliminate the on_PE_only flags above.
  !  if (.not.is_root_pe()) value = 0.0

end subroutine Ocean_stock_pe

!> This subroutine extracts a named 2-D field from the ocean surface or public type
subroutine ocean_model_data2D_get(OS, Ocean, name, array2D, isc, jsc)
  use MOM_constants, only : CELSIUS_KELVIN_OFFSET
  type(ocean_state_type),     pointer    :: OS    !< A pointer to the structure containing the
                                                  !! internal ocean state (intent in).
  type(ocean_public_type),    intent(in) :: Ocean !< A structure containing various publicly
                                                  !! visible ocean surface fields.
  character(len=*)          , intent(in) :: name  !< The name of the field to extract
  real, dimension(isc:,jsc:), intent(out):: array2D !< The values of the named field, it must
                                                  !! cover only the computational domain
  integer                   , intent(in) :: isc   !< The starting i-index of array2D
  integer                   , intent(in) :: jsc   !< The starting j-index of array2D

  integer :: g_isc, g_iec, g_jsc, g_jec,g_isd, g_ied, g_jsd, g_jed, i, j

  if (.not.associated(OS)) return
  if (.not.OS%is_ocean_pe) return

! The problem is %areaT is on MOM domain but Ice_Ocean_Boundary%... is on mpp domain.
! We want to return the MOM data on the mpp (compute) domain
! Get MOM domain extents
  call mpp_get_compute_domain(OS%grid%Domain%mpp_domain, g_isc, g_iec, g_jsc, g_jec)
  call mpp_get_data_domain   (OS%grid%Domain%mpp_domain, g_isd, g_ied, g_jsd, g_jed)

  g_isc = g_isc-g_isd+1 ; g_iec = g_iec-g_isd+1 ; g_jsc = g_jsc-g_jsd+1 ; g_jec = g_jec-g_jsd+1


  select case(name)
  case('area')
     array2D(isc:,jsc:) = OS%US%L_to_m**2*OS%grid%areaT(g_isc:g_iec,g_jsc:g_jec)
  case('mask')
     array2D(isc:,jsc:) = OS%grid%mask2dT(g_isc:g_iec,g_jsc:g_jec)
!OR same result
!     do j=g_jsc,g_jec ; do i=g_isc,g_iec
!        array2D(isc+i-g_isc,jsc+j-g_jsc) = OS%grid%mask2dT(i,j)
!     enddo ; enddo
  case('t_surf')
     array2D(isc:,jsc:) = Ocean%t_surf(isc:,jsc:)-CELSIUS_KELVIN_OFFSET
  case('t_pme')
     array2D(isc:,jsc:) = Ocean%t_surf(isc:,jsc:)-CELSIUS_KELVIN_OFFSET
  case('t_runoff')
     array2D(isc:,jsc:) = Ocean%t_surf(isc:,jsc:)-CELSIUS_KELVIN_OFFSET
  case('t_calving')
     array2D(isc:,jsc:) = Ocean%t_surf(isc:,jsc:)-CELSIUS_KELVIN_OFFSET
  case('btfHeat')
     array2D(isc:,jsc:) = 0
  case('cos_rot')
     array2D(isc:,jsc:) = OS%grid%cos_rot(g_isc:g_iec,g_jsc:g_jec) ! =1
  case('sin_rot')
     array2D(isc:,jsc:) = OS%grid%sin_rot(g_isc:g_iec,g_jsc:g_jec) ! =0
  case('s_surf')
     array2D(isc:,jsc:) = Ocean%s_surf(isc:,jsc:)
  case('sea_lev')
     array2D(isc:,jsc:) = Ocean%sea_lev(isc:,jsc:)
  case('frazil')
     array2D(isc:,jsc:) = Ocean%frazil(isc:,jsc:)
  case default
     call MOM_error(FATAL,'get_ocean_grid_data2D: unknown argument name='//name)
  end select

end subroutine ocean_model_data2D_get

!> This subroutine extracts a named scalar field from the ocean surface or public type
subroutine ocean_model_data1D_get(OS, Ocean, name, value)
  type(ocean_state_type),     pointer    :: OS    !< A pointer to the structure containing the
                                                  !! internal ocean state (intent in).
  type(ocean_public_type),    intent(in) :: Ocean !< A structure containing various publicly
                                                  !! visible ocean surface fields.
  character(len=*)          , intent(in) :: name  !< The name of the field to extract
  real                      , intent(out):: value !< The value of the named field

  if (.not.associated(OS)) return
  if (.not.OS%is_ocean_pe) return

  select case(name)
  case('c_p')
    value = OS%C_p
  case default
    call MOM_error(FATAL,'get_ocean_grid_data1D: unknown argument name='//name)
  end select

end subroutine ocean_model_data1D_get

!> Write out FMS-format checsums on fields from the ocean surface state
subroutine ocean_public_type_chksum(id, timestep, ocn)

  character(len=*),        intent(in) :: id  !< An identifying string for this call
  integer,                 intent(in) :: timestep !< The number of elapsed timesteps
  type(ocean_public_type), intent(in) :: ocn !< A structure containing various publicly
                                             !! visible ocean surface fields.
  integer :: n, m, outunit

  outunit = stdout()

  write(outunit,*) "BEGIN CHECKSUM(ocean_type):: ", id, timestep
  write(outunit,100) 'ocean%t_surf   ',mpp_chksum(ocn%t_surf )
  write(outunit,100) 'ocean%s_surf   ',mpp_chksum(ocn%s_surf )
  write(outunit,100) 'ocean%u_surf   ',mpp_chksum(ocn%u_surf )
  write(outunit,100) 'ocean%v_surf   ',mpp_chksum(ocn%v_surf )
  write(outunit,100) 'ocean%sea_lev  ',mpp_chksum(ocn%sea_lev)
  write(outunit,100) 'ocean%frazil   ',mpp_chksum(ocn%frazil )

  call coupler_type_write_chksums(ocn%fields, outunit, 'ocean%')
100 FORMAT("   CHECKSUM::",A20," = ",Z20)

end subroutine ocean_public_type_chksum

!> This subroutine gives a handle to the grid from ocean state
subroutine get_ocean_grid(OS, Gridp)
  ! Obtain the ocean grid.
  type(ocean_state_type) :: OS              !< A structure containing the
                                            !! internal ocean state
  type(ocean_grid_type) , pointer :: Gridp  !< The ocean's grid structure

  Gridp => OS%grid
  return
end subroutine get_ocean_grid

!> This subroutine extracts a named (u- or v-) 2-D surface current from ocean internal state
subroutine ocean_model_get_UV_surf(OS, Ocean, name, array2D, isc, jsc)

  type(ocean_state_type),     pointer    :: OS    !< A pointer to the structure containing the
                                                  !! internal ocean state (intent in).
  type(ocean_public_type),    intent(in) :: Ocean !< A structure containing various publicly
                                                  !! visible ocean surface fields.
  character(len=*)          , intent(in) :: name  !< The name of the current (ua or va) to extract
  real, dimension(isc:,jsc:), intent(out):: array2D !< The values of the named field, it must
                                                  !! cover only the computational domain
  integer                   , intent(in) :: isc   !< The starting i-index of array2D
  integer                   , intent(in) :: jsc   !< The starting j-index of array2D

  type(ocean_grid_type) , pointer :: G            !< The ocean's grid structure
  type(surface),          pointer :: sfc_state    !< A structure containing fields that
                                                  !! describe the surface state of the ocean.

  integer :: isc_bnd, iec_bnd, jsc_bnd, jec_bnd
  integer :: i, j, i0, j0
  integer :: is, ie, js, je

  if (.not.associated(OS)) return
  if (.not.OS%is_ocean_pe) return

  G => OS%grid
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  call mpp_get_compute_domain(Ocean%Domain, isc_bnd, iec_bnd, &
                              jsc_bnd, jec_bnd)

  i0 = is - isc_bnd ; j0 = js - jsc_bnd

  sfc_state => OS%sfc_state

  select case(name)
  case('ua')
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      array2D(i,j) = G%mask2dT(i+i0,j+j0) * &
                0.5*(sfc_state%u(I+i0,j+j0)+sfc_state%u(I-1+i0,j+j0))
    enddo ; enddo
  case('va')
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      array2D(i,j) = G%mask2dT(i+i0,j+j0) * &
                0.5*(sfc_state%v(i+i0,J+j0)+sfc_state%v(i+i0,J-1+j0))
    enddo ; enddo
  case('ub')
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      array2D(i,j) = G%mask2dBu(I+i0,J+j0) * &
                0.5*(sfc_state%u(I+i0,j+j0)+sfc_state%u(I+i0,j+j0+1))
    enddo ; enddo
  case('vb')
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      array2D(i,j) = G%mask2dBu(I+i0,J+j0) * &
                0.5*(sfc_state%v(i+i0,J+j0)+sfc_state%v(i+i0+1,J+j0))
    enddo ; enddo
  case default
     call MOM_error(FATAL,'ocean_model_get_UV_surf: unknown argument name='//name)
  end select

end subroutine ocean_model_get_UV_surf

end module ocean_model_mod
