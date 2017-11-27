module ocean_model_mod

! This file is part of MOM6. See LICENSE.md for the license.

!-----------------------------------------------------------------------
!
! This is the top level module for the MOM6 ocean model.  It contains routines
! for initialization, termination and update of ocean model state.  This
! particular version wraps all of the calls for MOM6 in the calls that had
! been used for MOM4.
!
! <CONTACT EMAIL="Robert.Hallberg@noaa.gov"> Robert Hallberg
! </CONTACT>
!
!<OVERVIEW>
! This code is a stop-gap wrapper of the MOM6 code to enable it to be called
! in the same way as MOM4.
!</OVERVIEW>

use MOM, only : initialize_MOM, step_MOM, MOM_control_struct, MOM_end
use MOM, only : calculate_surface_state, allocate_surface_state, finish_MOM_initialization
use MOM, only : step_offline
use MOM_constants, only : CELSIUS_KELVIN_OFFSET, hlf
use MOM_diag_mediator, only : diag_ctrl, enable_averaging, disable_averaging
use MOM_diag_mediator, only : diag_mediator_close_registration, diag_mediator_end
use MOM_domains, only : pass_vector, AGRID, BGRID_NE, CGRID_NE
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave
use MOM_file_parser, only : get_param, log_version, close_param_file, param_file_type
use MOM_forcing_type, only : forcing, mech_forcing
use MOM_forcing_type, only : forcing_accumulate, copy_common_forcing_fields
use MOM_forcing_type, only : copy_back_forcing_fields, set_net_mass_forcing
use MOM_forcing_type, only : set_derived_forcing_fields
use MOM_forcing_type, only : forcing_diagnostics, mech_forcing_diags
use MOM_get_input, only : Get_MOM_Input, directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : close_file, file_exists, read_data, write_version_number
use MOM_restart, only : save_restart
use MOM_sum_output, only : write_energy, accumulate_net_input
use MOM_sum_output, only : MOM_sum_output_init, sum_output_CS
use MOM_string_functions, only : uppercase
use MOM_surface_forcing, only : surface_forcing_init, convert_IOB_to_fluxes
use MOM_surface_forcing, only : ice_ocn_bnd_type_chksum
use MOM_surface_forcing, only : ice_ocean_boundary_type, surface_forcing_CS
use MOM_surface_forcing, only : forcing_save_restart
use MOM_time_manager, only : time_type, get_time, set_time, operator(>)
use MOM_time_manager, only : operator(+), operator(-), operator(*), operator(/)
use MOM_time_manager, only : operator(/=), operator(<=), operator(>=)
use MOM_time_manager, only : operator(<), real_to_time_type, time_type_to_real
use MOM_tracer_flow_control, only : call_tracer_register, tracer_flow_control_init
use MOM_tracer_flow_control, only : call_tracer_flux_init
use MOM_variables, only : surface
use MOM_verticalGrid, only : verticalGrid_type
use MOM_ice_shelf, only : initialize_ice_shelf, shelf_calc_flux, ice_shelf_CS
use MOM_ice_shelf, only : ice_shelf_end, ice_shelf_save_restart
use coupler_types_mod, only : coupler_1d_bc_type, coupler_2d_bc_type
use coupler_types_mod, only : coupler_type_spawn, coupler_type_write_chksums
use coupler_types_mod, only : coupler_type_initialized, coupler_type_copy_data
use coupler_types_mod, only : coupler_type_set_diags, coupler_type_send_data
use mpp_domains_mod, only : domain2d, mpp_get_layout, mpp_get_global_domain
use mpp_domains_mod, only : mpp_define_domains, mpp_get_compute_domain, mpp_get_data_domain
use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux
use MOM_forcing_type, only : allocate_forcing_type
use fms_mod, only : stdout
use mpp_mod, only : mpp_chksum
use MOM_domains, only : pass_var, pass_vector, TO_ALL, CGRID_NE, BGRID_NE
use MOM_EOS, only : gsw_sp_from_sr, gsw_pt_from_ct

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
    u_surf => NULL(), & !< i-velocity at the locations indicated by stagger, m/s.
    v_surf => NULL(), & !< j-velocity at the locations indicated by stagger, m/s.
    sea_lev => NULL(), & !< Sea level in m after correction for surface pressure,
                        !! i.e. dzt(1) + eta_t + patm/rho0/grav (m)
    frazil =>NULL(), &  !< Accumulated heating (in Joules/m^2) from frazil
                        !! formation in the ocean.
    area => NULL()      !< cell area of the ocean surface, in m2.
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
  type(time_type) :: Time    ! The ocean model's time and master clock.
  integer :: Restart_control ! An integer that is bit-tested to determine whether
                             ! incremental restart files are saved and whether they
                             ! have a time stamped name.  +1 (bit 0) for generic
                             ! files and +2 (bit 1) for time-stamped files.  A
                             ! restart file is saved at the end of a run segment
                             ! unless Restart_control is negative.

  type(time_type) :: energysavedays            ! The interval between writing the energies
                                               ! and other integral quantities of the run.
  type(time_type) :: energysavedays_geometric  ! The starting interval for computing a geometric
                                               ! progression of time deltas between calls to
                                               ! write_energy. This interval will increase by a factor of 2.
                                               ! after each call to write_energy.
  logical         :: energysave_geometric      ! Logical to control whether calls to write_energy should
                                               ! follow a geometric progression
  type(time_type) :: write_energy_time         ! The next time to write to the energy file.
  type(time_type) :: geometric_end_time        ! Time at which to stop the geometric progression
                                               ! of calls to write_energy and revert to the standard
                                               ! energysavedays interval

  integer :: nstep = 0        ! The number of calls to update_ocean.
  logical :: use_ice_shelf    ! If true, the ice shelf model is enabled.

  ! Many of the following variables do not appear to belong here. -RWH
  logical :: icebergs_apply_rigid_boundary  ! If true, the icebergs can change ocean bd condition.
  real :: kv_iceberg          ! The viscosity of the icebergs in m2/s (for ice rigidity)
  real :: berg_area_threshold ! Fraction of grid cell which iceberg must occupy
                              !so that fluxes below are set to zero. (0.5 is a
                              !good value to use. Not applied for negative values.
  real :: latent_heat_fusion  ! Latent heat of fusion
  real :: density_iceberg     ! A typical density of icebergs in kg/m3 (for ice rigidity)

  type(ice_shelf_CS), pointer :: Ice_shelf_CSp => NULL()
  logical :: restore_salinity ! If true, the coupled MOM driver adds a term to
                              ! restore salinity to a specified value.
  logical :: restore_temp     ! If true, the coupled MOM driver adds a term to
                              ! restore sst to a specified value.
  real :: press_to_z          ! A conversion factor between pressure and ocean
                              ! depth in m, usually 1/(rho_0*g), in m Pa-1.
  real :: C_p                 ! The heat capacity of seawater, in J K-1 kg-1.

  type(directories) :: dirs   ! A structure containing several relevant directory paths.
  type(mech_forcing) :: forces !< A structure with the driving mechanical surface forces
  type(forcing)   :: fluxes   ! A structure containing pointers to
                              ! the thermodynamic ocean forcing fields.
  type(forcing)   :: flux_tmp ! A secondary structure containing pointers to the
                              ! ocean forcing fields for when multiple coupled
                              ! timesteps are taken per thermodynamic step.
  type(surface)   :: sfc_state ! A structure containing pointers to
                              ! the ocean surface state fields.
  type(ocean_grid_type), pointer :: grid => NULL() ! A pointer to a grid structure
                              ! containing metrics and related information.
  type(verticalGrid_type), pointer :: GV => NULL() ! A pointer to a vertical grid
                              ! structure containing metrics and related information.
  type(MOM_control_struct), pointer :: MOM_CSp => NULL()
  type(surface_forcing_CS), pointer :: forcing_CSp => NULL()
  type(sum_output_CS),      pointer :: sum_output_CSp => NULL()
end type ocean_state_type

contains

!=======================================================================
! <SUBROUTINE NAME="ocean_model_init">
!
! <DESCRIPTION>
! Initialize the ocean model.
! </DESCRIPTION>

!> ocean_model_init initializes the ocean model, including registering fields
!! for restarts and reading restart files if appropriate.
subroutine ocean_model_init(Ocean_sfc, OS, Time_init, Time_in, gas_fields_ocn)
  type(ocean_public_type), target, &
                       intent(inout) :: Ocean_sfc !< A structure containing various
                                !! publicly visible ocean surface properties after initialization,
                                !! the data in this type is intent(out).
  type(ocean_state_type), pointer    :: OS        !< A structure whose internal
                                !! contents are private to ocean_model_mod that may be used to
                                !! contain all information about the ocean's interior state.
  type(time_type),     intent(in)    :: Time_init !< The start time for the coupled model's calendar
  type(time_type),     intent(in)    :: Time_in   !< The time at which to initialize the ocean model.
  type(coupler_1d_bc_type), &
             optional, intent(in)    :: gas_fields_ocn !< If present, this type describes the
                                              !! ocean and surface-ice fields that will participate
                                              !! in the calculation of additional gas or other
                                              !! tracer fluxes, and can be used to spawn related
                                              !! internal variables in the ice model.

!   This subroutine initializes both the ocean state and the ocean surface type.
! Because of the way that indicies and domains are handled, Ocean_sfc must have
! been used in a previous call to initialize_ocean_type.

! Arguments: Ocean_sfc - A structure containing various publicly visible ocean
!                    surface properties after initialization, this is intent(out).
!  (out,private) OS - A structure whose internal contents are private
!                    to ocean_model_mod that may be used to contain all
!                    information about the ocean's interior state.
!  (in)      Time_init - The start time for the coupled model's calendar.
!  (in)      Time_in - The time at which to initialize the ocean model.
  real :: Time_unit   ! The time unit in seconds for ENERGYSAVEDAYS.
  real :: Rho0        ! The Boussinesq ocean density, in kg m-3.
  real :: G_Earth     ! The gravitational acceleration in m s-2.
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "ocean_model_init"  ! This module's name.
  character(len=48)  :: stagger
  integer :: secs, days
  type(param_file_type) :: param_file !< A structure to parse for run-time parameters
  logical :: offline_tracer_mode
  type(time_type) :: dt_geometric, dt_savedays, dt_from_base

  call callTree_enter("ocean_model_init(), ocean_model_MOM.F90")
  if (associated(OS)) then
    call MOM_error(WARNING, "ocean_model_init called with an associated "// &
                    "ocean_state_type structure. Model is already initialized.")
    return
  endif
  allocate(OS)

  OS%is_ocean_pe = Ocean_sfc%is_ocean_pe
  if (.not.OS%is_ocean_pe) return

  OS%Time = Time_in
  call initialize_MOM(OS%Time, param_file, OS%dirs, OS%MOM_CSp, Time_in, &
      offline_tracer_mode=offline_tracer_mode)
  OS%grid => OS%MOM_CSp%G ; OS%GV => OS%MOM_CSp%GV
  OS%C_p = OS%MOM_CSp%tv%C_p
  OS%fluxes%C_p = OS%MOM_CSp%tv%C_p

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "RESTART_CONTROL", OS%Restart_control, &
                 "An integer whose bits encode which restart files are \n"//&
                 "written. Add 2 (bit 1) for a time-stamped file, and odd \n"//&
                 "(bit 0) for a non-time-stamped file.  A restart file \n"//&
                 "will be saved at the end of the run segment for any \n"//&
                 "non-negative value.", default=1)

  call get_param(param_file, mdl, "TIMEUNIT", Time_unit, &
                 "The time unit for ENERGYSAVEDAYS.", &
                 units="s", default=86400.0)
  call get_param(param_file, mdl, "ENERGYSAVEDAYS",OS%energysavedays, &
                 "The interval in units of TIMEUNIT between saves of the \n"//&
                 "energies of the run and other globally summed diagnostics.",&
                 default=set_time(0,days=1), timeunit=Time_unit)
  call get_param(param_file, mdl, "ENERGYSAVEDAYS_GEOMETRIC",OS%energysavedays_geometric, &
                 "The starting interval in units of TIMEUNIT for the first call \n"//&
                 "to save the energies of the run and other globally summed diagnostics. \n"//&
                 "The interval increases by a factor of 2. after each call to write_energy.",&
                 default=set_time(seconds=0), timeunit=Time_unit)

  if ((time_type_to_real(OS%energysavedays_geometric) > 0.) .and. &
     (OS%energysavedays_geometric < OS%energysavedays)) then
         OS%energysave_geometric = .true.
  else
         OS%energysave_geometric = .false.
  endif

  call get_param(param_file, mdl, "OCEAN_SURFACE_STAGGER", stagger, &
                 "A case-insensitive character string to indicate the \n"//&
                 "staggering of the surface velocity field that is \n"//&
                 "returned to the coupler.  Valid values include \n"//&
                 "'A', 'B', or 'C'.", default="C")
  if (uppercase(stagger(1:1)) == 'A') then ; Ocean_sfc%stagger = AGRID
  elseif (uppercase(stagger(1:1)) == 'B') then ; Ocean_sfc%stagger = BGRID_NE
  elseif (uppercase(stagger(1:1)) == 'C') then ; Ocean_sfc%stagger = CGRID_NE
  else ; call MOM_error(FATAL,"ocean_model_init: OCEAN_SURFACE_STAGGER = "// &
                        trim(stagger)//" is invalid.") ; endif

  call get_param(param_file, mdl, "RESTORE_SALINITY",OS%restore_salinity, &
                 "If true, the coupled driver will add a globally-balanced \n"//&
                 "fresh-water flux that drives sea-surface salinity \n"//&
                 "toward specified values.", default=.false.)
  call get_param(param_file, mdl, "RESTORE_TEMPERATURE",OS%restore_temp, &
                 "If true, the coupled driver will add a  \n"//&
                 "heat flux that drives sea-surface temperauture \n"//&
                 "toward specified values.", default=.false.)
  call get_param(param_file, mdl, "RHO_0", Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to \n"//&
                 "calculate accelerations and the mass for conservation \n"//&
                 "properties, or with BOUSSINSEQ false to convert some \n"//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0)
  call get_param(param_file, mdl, "G_EARTH", G_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80)

  call get_param(param_file, mdl, "ICE_SHELF",  OS%use_ice_shelf, &
                 "If true, enables the ice shelf model.", default=.false.)

  call get_param(param_file, mdl, "ICEBERGS_APPLY_RIGID_BOUNDARY",  OS%icebergs_apply_rigid_boundary, &
                 "If true, allows icebergs to change boundary condition felt by ocean", default=.false.)

  if (OS%icebergs_apply_rigid_boundary) then
    call get_param(param_file, mdl, "KV_ICEBERG",  OS%kv_iceberg, &
                 "The viscosity of the icebergs",  units="m2 s-1",default=1.0e10)
    call get_param(param_file, mdl, "DENSITY_ICEBERGS",  OS%density_iceberg, &
                 "A typical density of icebergs.", units="kg m-3", default=917.0)
    call get_param(param_file, mdl, "LATENT_HEAT_FUSION", OS%latent_heat_fusion, &
                 "The latent heat of fusion.", units="J/kg", default=hlf)
    call get_param(param_file, mdl, "BERG_AREA_THRESHOLD", OS%berg_area_threshold, &
                 "Fraction of grid cell which iceberg must occupy, so that fluxes \n"//&
                 "below berg are set to zero. Not applied for negative \n"//&
                 " values.", units="non-dim", default=-1.0)
  endif

  OS%press_to_z = 1.0/(Rho0*G_Earth)

  !   Consider using a run-time flag to determine whether to do the diagnostic
  ! vertical integrals, since the related 3-d sums are not negligible in cost.
  call allocate_surface_state(OS%sfc_state, OS%grid, OS%MOM_CSp%use_temperature, &
                              do_integrals=.true., gas_fields_ocn=gas_fields_ocn)

  call surface_forcing_init(Time_in, OS%grid, param_file, OS%MOM_CSp%diag, &
                            OS%forcing_CSp, OS%restore_salinity, OS%restore_temp)

  if (OS%use_ice_shelf)  then
    call initialize_ice_shelf(param_file, OS%grid, OS%Time, OS%ice_shelf_CSp, &
                              OS%MOM_CSp%diag, OS%forces, OS%fluxes)
  endif
  if (OS%icebergs_apply_rigid_boundary)  then
    !call allocate_forcing_type(OS%grid, OS%fluxes, iceberg=.true.)
    !This assumes that the iceshelf and ocean are on the same grid. I hope this is true
    if (.not. OS%use_ice_shelf) &
      call allocate_forcing_type(OS%grid, OS%fluxes, shelf=.true.)
  endif

  call MOM_sum_output_init(OS%grid, param_file, OS%dirs%output_directory, &
                            OS%MOM_CSp%ntrunc, Time_init, OS%sum_output_CSp)

  ! This call has been moved into the first call to update_ocean_model.
  !  call write_energy(OS%MOM_CSp%u, OS%MOM_CSp%v, OS%MOM_CSp%h, OS%MOM_CSp%tv, &
  !             OS%Time, 0, OS%grid, OS%GV, OS%sum_output_CSp, OS%MOM_CSp%tracer_flow_CSp)

  ! write_energy_time is the next integral multiple of energysavedays.
  if (OS%energysave_geometric) then
    if (OS%energysavedays_geometric < OS%energysavedays) then
      OS%write_energy_time = OS%Time + OS%energysavedays_geometric
      OS%geometric_end_time = Time_init + OS%energysavedays * &
       (1 + (OS%Time - Time_init) / OS%energysavedays)
    else
      OS%write_energy_time = Time_init + OS%energysavedays * &
        (1 + (OS%Time - Time_init) / OS%energysavedays)
    endif
  else
    OS%write_energy_time = Time_init + OS%energysavedays * &
      (1 + (OS%Time - Time_init) / OS%energysavedays)
  endif

  if (ASSOCIATED(OS%grid%Domain%maskmap)) then
    call initialize_ocean_public_type(OS%grid%Domain%mpp_domain, Ocean_sfc, &
                                      OS%MOM_CSp%diag, maskmap=OS%grid%Domain%maskmap, &
                                      gas_fields_ocn=gas_fields_ocn)
  else
    call initialize_ocean_public_type(OS%grid%Domain%mpp_domain, Ocean_sfc, &
                                      OS%MOM_CSp%diag, gas_fields_ocn=gas_fields_ocn)
  endif

  ! This call can only occur here if the coupler_bc_type variables have been
  ! initialized already using the information from gas_fields_ocn.
  if (present(gas_fields_ocn)) then
    call coupler_type_set_diags(Ocean_sfc%fields, "ocean_sfc", &
                                Ocean_sfc%axes(1:2), Time_in)

    call calculate_surface_state(OS%sfc_state, OS%MOM_CSp%u, &
             OS%MOM_CSp%v, OS%MOM_CSp%h, OS%MOM_CSp%ave_ssh,&
             OS%grid, OS%GV, OS%MOM_CSp)

    call convert_state_to_ocean_type(OS%sfc_state, Ocean_sfc, OS%grid, &
                                     OS%MOM_CSp%use_conT_absS)

  endif

  call close_param_file(param_file)
  call diag_mediator_close_registration(OS%MOM_CSp%diag)

  if (is_root_pe()) &
    write(*,'(/12x,a/)') '======== COMPLETED MOM INITIALIZATION ========'

  call callTree_leave("ocean_model_init(")
end subroutine ocean_model_init
! </SUBROUTINE> NAME="ocean_model_init"


!=======================================================================
! <SUBROUTINE NAME="update_ocean_model">
!
! <DESCRIPTION>
! Update in time the ocean model fields.  This code wraps the call to step_MOM
! with MOM4's call.
! </DESCRIPTION>
!

!> update_ocean_model uses the forcing in Ice_ocean_boundary to advance the
!! ocean model's state from the input value of Ocean_state (which must be for
!! time time_start_update) for a time interval of Ocean_coupling_time_step,
!! returning the publicly visible ocean surface properties in Ocean_sfc and
!! storing the new ocean properties in Ocean_state.
subroutine update_ocean_model(Ice_ocean_boundary, OS, Ocean_sfc, &
                              time_start_update, Ocean_coupling_time_step)
  type(ice_ocean_boundary_type), &
                           intent(in)    :: Ice_ocean_boundary !< A structure containing the
                                                    !! various forcing fields coming from the ice.
  type(ocean_state_type),  pointer       :: OS      !< A pointer to a private structure containing
                                                    !! the internal ocean state.
  type(ocean_public_type), intent(inout) :: Ocean_sfc !< A structure containing all the
                                                    !! publicly visible ocean surface fields after
                                                    !! a coupling time step.  The data in this type is
                                                    !! intent out.
  type(time_type),         intent(in)    :: time_start_update  !< The time at the beginning of the update step.
  type(time_type),         intent(in)    :: Ocean_coupling_time_step !< The amount of time over
                                                    !! which to advance the ocean.
!   This subroutine uses the forcing in Ice_ocean_boundary to advance the
! ocean model's state from the input value of Ocean_state (which must be for
! time time_start_update) for a time interval of Ocean_coupling_time_step,
! returning the publicly visible ocean surface properties in Ocean_sfc and
! storing the new ocean properties in Ocean_state.

! Arguments: Ice_ocean_boundary - A structure containing the various forcing
!                                 fields coming from the ice. It is intent in.
!  (inout)   Ocean_state - A structure containing the internal ocean state.
!  (out)     Ocean_sfc - A structure containing all the publicly visible ocean
!                        surface fields after a coupling time step.
!  (in)      time_start_update - The time at the beginning of the update step.
!  (in)      Ocean_coupling_time_step - The amount of time over which to advance
!                                       the ocean.

! Note: although several types are declared intent(inout), this is to allow for
!   the possibility of halo updates and to keep previously allocated memory.
!   In practice, Ice_ocean_boundary is intent in, Ocean_state is private to
!   this module and intent inout, and Ocean_sfc is intent out.
  type(time_type) :: Master_time ! This allows step_MOM to temporarily change
                                 ! the time that is seen by internal modules.
  type(time_type) :: Time1       ! The value of the ocean model's time at the
                                 ! start of a call to step_MOM.
  integer :: index_bnds(4)       ! The computational domain index bounds in the
                                 ! ice-ocean boundary type.
  real :: weight            ! Flux accumulation weight
  real :: time_step         ! The time step of a call to step_MOM in seconds.
  integer :: secs, days
  integer :: is, ie, js, je
  type(time_type) :: write_energy_time_geometric

  call callTree_enter("update_ocean_model(), ocean_model_MOM.F90")
  call get_time(Ocean_coupling_time_step, secs, days)
  time_step = 86400.0*real(days) + real(secs)

  if (time_start_update /= OS%Time) then
    call MOM_error(WARNING, "update_ocean_model: internal clock does not "//&
                            "agree with time_start_update argument.")
  endif
  if (.not.associated(OS)) then
    call MOM_error(FATAL, "update_ocean_model called with an unassociated "// &
                    "ocean_state_type structure. ocean_model_init must be "//  &
                    "called first to allocate this structure.")
    return
  endif

  ! This is benign but not necessary if ocean_model_init_sfc was called or if
  ! OS%sfc_state%tr_fields was spawned in ocean_model_init.  Consider removing it.
  is = OS%grid%isc ; ie = OS%grid%iec ; js = OS%grid%jsc ; je = OS%grid%jec
  call coupler_type_spawn(Ocean_sfc%fields, OS%sfc_state%tr_fields, &
                          (/is,is,ie,ie/), (/js,js,je,je/), as_needed=.true.)

  ! Translate Ice_ocean_boundary into fluxes.
  call mpp_get_compute_domain(Ocean_sfc%Domain, index_bnds(1), index_bnds(2), &
                              index_bnds(3), index_bnds(4))

  weight = 1.0

  if (OS%fluxes%fluxes_used) then
    call enable_averaging(time_step, OS%Time + Ocean_coupling_time_step, OS%MOM_CSp%diag) ! Needed to allow diagnostics in convert_IOB
    call convert_IOB_to_fluxes(Ice_ocean_boundary, OS%forces, OS%fluxes, index_bnds, OS%Time, &
                               OS%grid, OS%forcing_CSp, OS%sfc_state, OS%restore_salinity,OS%restore_temp)

    ! Add ice shelf fluxes
    if (OS%use_ice_shelf) then
      call shelf_calc_flux(OS%sfc_state, OS%forces, OS%fluxes, OS%Time, time_step, OS%Ice_shelf_CSp)
    endif
    if (OS%icebergs_apply_rigid_boundary)  then
      !This assumes that the iceshelf and ocean are on the same grid. I hope this is true
      call add_berg_flux_to_shelf(OS%grid, OS%forces, OS%fluxes, OS%use_ice_shelf, &
                    OS%density_iceberg,OS%kv_iceberg, OS%latent_heat_fusion, OS%sfc_state, &
                    time_step, OS%berg_area_threshold)
    endif

    ! Fields that exist in both the forcing and mech_forcing types must be copied.
    call copy_common_forcing_fields(OS%forces, OS%fluxes, OS%grid)

#ifdef _USE_GENERIC_TRACER
    call MOM_generic_tracer_fluxes_accumulate(OS%fluxes, weight) !here weight=1, just saving the current fluxes
#endif
    ! Indicate that there are new unused fluxes.
    OS%fluxes%fluxes_used = .false.
    OS%fluxes%dt_buoy_accum = time_step
  else
    OS%flux_tmp%C_p = OS%fluxes%C_p
    call convert_IOB_to_fluxes(Ice_ocean_boundary, OS%forces, OS%flux_tmp, index_bnds, OS%Time, &
                               OS%grid, OS%forcing_CSp, OS%sfc_state, OS%restore_salinity,OS%restore_temp)
    if (OS%use_ice_shelf) then
      call shelf_calc_flux(OS%sfc_state, OS%forces, OS%flux_tmp, OS%Time, time_step, OS%Ice_shelf_CSp)
    endif
    if (OS%icebergs_apply_rigid_boundary)  then
     !This assumes that the iceshelf and ocean are on the same grid. I hope this is true
     call add_berg_flux_to_shelf(OS%grid, OS%forces, OS%flux_tmp, OS%use_ice_shelf, OS%density_iceberg, &
            OS%kv_iceberg, OS%latent_heat_fusion, OS%sfc_state, time_step, OS%berg_area_threshold)
    endif

    call forcing_accumulate(OS%flux_tmp, OS%forces, OS%fluxes, time_step, OS%grid, weight)
    ! Some of the fields that exist in both the forcing and mech_forcing types
    ! are time-averages must be copied back to the forces type.
    call copy_back_forcing_fields(OS%fluxes, OS%forces, OS%grid)

#ifdef _USE_GENERIC_TRACER
    call MOM_generic_tracer_fluxes_accumulate(OS%flux_tmp, weight) !weight of the current flux in the running average
#endif
  endif
  call set_derived_forcing_fields(OS%forces, OS%fluxes, OS%grid, OS%GV%Rho0)
  call set_net_mass_forcing(OS%fluxes, OS%forces, OS%grid)

  if (OS%nstep==0) then
    call finish_MOM_initialization(OS%Time, OS%dirs, OS%MOM_CSp, OS%fluxes)

    call write_energy(OS%MOM_CSp%u, OS%MOM_CSp%v, OS%MOM_CSp%h, OS%MOM_CSp%tv, &
                      OS%Time, 0, OS%grid, OS%GV, OS%sum_output_CSp, &
                      OS%MOM_CSp%tracer_flow_CSp)
  endif

  call disable_averaging(OS%MOM_CSp%diag)
  Master_time = OS%Time ; Time1 = OS%Time

  if(OS%MOM_Csp%offline_tracer_mode) then
    call step_offline(OS%forces, OS%fluxes, OS%sfc_state, Time1, time_step, OS%MOM_CSp)
  else
    call step_MOM(OS%forces, OS%fluxes, OS%sfc_state, Time1, time_step, OS%MOM_CSp)
  endif

  OS%Time = Master_time + Ocean_coupling_time_step
  OS%nstep = OS%nstep + 1

  call enable_averaging(time_step, OS%Time, OS%MOM_CSp%diag)
  call mech_forcing_diags(OS%forces, OS%fluxes, time_step, OS%grid, &
                          OS%MOM_CSp%diag, OS%forcing_CSp%handles)
  call disable_averaging(OS%MOM_CSp%diag)

  if (OS%fluxes%fluxes_used) then
    call enable_averaging(OS%fluxes%dt_buoy_accum, OS%Time, OS%MOM_CSp%diag)
    call forcing_diagnostics(OS%fluxes, OS%sfc_state, OS%fluxes%dt_buoy_accum, &
                             OS%grid, OS%MOM_CSp%diag, OS%forcing_CSp%handles)
    call accumulate_net_input(OS%fluxes, OS%sfc_state, OS%fluxes%dt_buoy_accum, &
                              OS%grid, OS%sum_output_CSp)
    call disable_averaging(OS%MOM_CSp%diag)
  endif

!  See if it is time to write out the energy.

  if (OS%energysave_geometric) then
    if ((OS%Time + ((Ocean_coupling_time_step)/2) > OS%geometric_end_time) .and. &
        (OS%MOM_CSp%t_dyn_rel_adv==0.0)) then
        call write_energy(OS%MOM_CSp%u, OS%MOM_CSp%v, OS%MOM_CSp%h, OS%MOM_CSp%tv, &
                        OS%Time, OS%nstep, OS%grid, OS%GV, OS%sum_output_CSp, &
                        OS%MOM_CSp%tracer_flow_CSp)
        OS%write_energy_time = OS%geometric_end_time + OS%energysavedays
        OS%energysave_geometric = .false.  ! stop geometric progression
    endif
  endif

  if ((OS%Time + ((Ocean_coupling_time_step)/2) > OS%write_energy_time) .and. &
      (OS%MOM_CSp%t_dyn_rel_adv==0.0)) then
    call write_energy(OS%MOM_CSp%u, OS%MOM_CSp%v, OS%MOM_CSp%h, OS%MOM_CSp%tv, &
                      OS%Time, OS%nstep, OS%grid, OS%GV, OS%sum_output_CSp, &
                      OS%MOM_CSp%tracer_flow_CSp)
    if (OS%energysave_geometric) then
        OS%energysavedays_geometric = OS%energysavedays_geometric*2
        OS%write_energy_time = OS%write_energy_time + OS%energysavedays_geometric
    else
      OS%write_energy_time = OS%write_energy_time + OS%energysavedays
    endif
  endif



! Translate state into Ocean.
!  call convert_state_to_ocean_type(OS%sfc_state, Ocean_sfc, OS%grid, &
!                                   Ice_ocean_boundary%p, OS%press_to_z)
  call convert_state_to_ocean_type(OS%sfc_state, Ocean_sfc, OS%grid, &
                                   OS%MOM_CSp%use_conT_absS)
  call coupler_type_send_data(Ocean_sfc%fields, OS%Time)

  call callTree_leave("update_ocean_model()")
end subroutine update_ocean_model
! </SUBROUTINE> NAME="update_ocean_model"

!=======================================================================
! <SUBROUTINE NAME="ocean_model_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments:
!   timestamp (optional, intent(in)) : A character string that represents the model time,
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix.
! </DESCRIPTION>
!

subroutine add_berg_flux_to_shelf(G, forces, fluxes, use_ice_shelf, density_ice, kv_ice, &
                                  latent_heat_fusion, sfc_state, time_step, berg_area_threshold)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
  type(mech_forcing),    intent(inout) :: forces  !< A structure with the driving mechanical forces
  type(forcing),         intent(inout) :: fluxes
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  logical,               intent(in)    :: use_ice_shelf  !< If true, this configuration uses ice shelves.
  real, intent(in) :: kv_ice       ! The viscosity of ice, in m2 s-1.
  real, intent(in) :: density_ice  ! A typical density of ice, in kg m-3.
  real, intent(in) :: latent_heat_fusion   ! The latent heat of fusion, in J kg-1.
  real, intent(in) :: time_step   ! The latent heat of fusion, in J kg-1.
  real, intent(in) :: berg_area_threshold  !Area threshold for zero'ing fluxes bellow iceberg
! Arguments:
!  (in)      fluxes - A structure of surface fluxes that may be used.
!  (in)      G - The ocean's grid structure.
  real :: fraz          ! refreezing rate in kg m-2 s-1
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed
  !This routine adds iceberg data to the ice shelf data (if ice shelf is used)
  !which can then be used to change the top of ocean boundary condition used in
  !the ocean model. This routine is taken from the add_shelf_flux subroutine
  !within the ice shelf model.

  if (.not. (((associated(fluxes%frac_shelf_h) .and. associated(forces%frac_shelf_u)) &
    .and.(associated(forces%frac_shelf_v) .and. associated(fluxes%ustar_shelf)))&
    .and.(associated(forces%rigidity_ice_u) .and. associated(forces%rigidity_ice_v)))) return

  if (.not.(associated(fluxes%area_berg) .and. associated(fluxes%ustar_berg) .and. &
            associated(fluxes%mass_berg) ) ) return

  if (.not. use_ice_shelf) then
    fluxes%frac_shelf_h(:,:) = 0.
    forces%frac_shelf_u(:,:) = 0.
    forces%frac_shelf_v(:,:) = 0.
    fluxes%ustar_shelf(:,:) = 0.
    forces%rigidity_ice_u(:,:) = 0.
    forces%rigidity_ice_v(:,:) = 0.
  endif

  do j=jsd,jed ; do i=isd,ied
    if (G%areaT(i,j) > 0.0) &
      fluxes%frac_shelf_h(i,j) = fluxes%frac_shelf_h(i,j) +  fluxes%area_berg(i,j)
      fluxes%ustar_shelf(i,j)  = fluxes%ustar_shelf(i,j)  +  fluxes%ustar_berg(i,j)
  enddo ; enddo
  do j=jsd,jed ; do I=isd,ied-1
    forces%frac_shelf_u(I,j) = 0.0
    if ((G%areaT(i,j) + G%areaT(i+1,j) > 0.0)) & ! .and. (G%dxdy_u(I,j) > 0.0)) &
      forces%frac_shelf_u(I,j) = forces%frac_shelf_u(I,j) + &
          (((fluxes%area_berg(i,j)*G%areaT(i,j)) + (fluxes%area_berg(i+1,j)*G%areaT(i+1,j))) / &
           (G%areaT(i,j) + G%areaT(i+1,j)) )
    !### Either the min here or the max below must be wrong, but is either right? -RWH
    forces%rigidity_ice_u(I,j) = forces%rigidity_ice_u(I,j) +((kv_ice / density_ice) * &
                                  min(fluxes%mass_berg(i,j), fluxes%mass_berg(i+1,j)))
  enddo ; enddo
  do J=jsd,jed-1 ; do i=isd,ied
    forces%frac_shelf_v(i,J) = 0.0
    if ((G%areaT(i,j) + G%areaT(i,j+1) > 0.0)) & ! .and. (G%dxdy_v(i,J) > 0.0)) &
      forces%frac_shelf_v(i,J) = forces%frac_shelf_v(i,J) + &
          (((fluxes%area_berg(i,j)*G%areaT(i,j)) + (fluxes%area_berg(i,j+1)*G%areaT(i,j+1))) / &
           (G%areaT(i,j) + G%areaT(i,j+1)) )
    !### Either the max here or the min above must be wrong, but is either right? -RWH
    forces%rigidity_ice_v(i,J) = forces%rigidity_ice_v(i,J) +((kv_ice / density_ice) * &
                                  max(fluxes%mass_berg(i,j), fluxes%mass_berg(i,j+1)))
  enddo ; enddo
  call pass_vector(forces%frac_shelf_u, forces%frac_shelf_v, G%domain, TO_ALL, CGRID_NE)

  !Zero'ing out other fluxes under the tabular icebergs
  if (berg_area_threshold >= 0.) then
    do j=jsd,jed ; do i=isd,ied
      if (fluxes%frac_shelf_h(i,j) > berg_area_threshold) then  !Only applying for ice shelf covering most of cell

        if (associated(fluxes%sw)) fluxes%sw(i,j) = 0.0
        if (associated(fluxes%lw)) fluxes%lw(i,j) = 0.0
        if (associated(fluxes%latent)) fluxes%latent(i,j) = 0.0
        if (associated(fluxes%evap)) fluxes%evap(i,j) = 0.0

        ! Add frazil formation diagnosed by the ocean model (J m-2) in the
        ! form of surface layer evaporation (kg m-2 s-1). Update lprec in the
        ! control structure for diagnostic purposes.

        if (associated(sfc_state%frazil)) then
          fraz = sfc_state%frazil(i,j) / time_step / latent_heat_fusion
          if (associated(fluxes%evap)) fluxes%evap(i,j) = fluxes%evap(i,j) - fraz
          !CS%lprec(i,j)=CS%lprec(i,j) - fraz
          sfc_state%frazil(i,j) = 0.0
        endif

        !Alon: Should these be set to zero too?
        if (associated(fluxes%sens)) fluxes%sens(i,j) = 0.0
        if (associated(fluxes%salt_flux)) fluxes%salt_flux(i,j) = 0.0
        if (associated(fluxes%lprec)) fluxes%lprec(i,j) = 0.0
      endif
    enddo ; enddo
  endif

end subroutine add_berg_flux_to_shelf

subroutine ocean_model_restart(OS, timestamp)
  type(ocean_state_type),        pointer :: OS
  character(len=*), intent(in), optional :: timestamp

  if (OS%MOM_CSp%t_dyn_rel_adv > 0.0) call MOM_error(WARNING, "End of MOM_main reached "//&
       "with inconsistent dynamics and advective times.  Additional restart fields "//&
       "that have not been coded yet would be required for reproducibility.")
  if (.not.OS%fluxes%fluxes_used) call MOM_error(FATAL, "ocean_model_restart "//&
      "was called with unused buoyancy fluxes.  For conservation, the ocean "//&
      "restart files can only be created after the buoyancy forcing is applied.")

  if (BTEST(OS%Restart_control,1)) then
    call save_restart(OS%dirs%restart_output_dir, OS%Time, OS%grid, &
                      OS%MOM_CSp%restart_CSp, .true., GV=OS%GV)
    call forcing_save_restart(OS%forcing_CSp, OS%grid, OS%Time, &
                              OS%dirs%restart_output_dir, .true.)
    if (OS%use_ice_shelf) then
      call ice_shelf_save_restart(OS%Ice_shelf_CSp, OS%Time, OS%dirs%restart_output_dir, .true.)
    endif
  endif
  if (BTEST(OS%Restart_control,0)) then
    call save_restart(OS%dirs%restart_output_dir, OS%Time, OS%grid, &
                      OS%MOM_CSp%restart_CSp, GV=OS%GV)
    call forcing_save_restart(OS%forcing_CSp, OS%grid, OS%Time, &
                              OS%dirs%restart_output_dir)
    if (OS%use_ice_shelf) then
      call ice_shelf_save_restart(OS%Ice_shelf_CSp, OS%Time, OS%dirs%restart_output_dir)
    endif
  endif

end subroutine ocean_model_restart
! </SUBROUTINE> NAME="ocean_model_restart"

!=======================================================================
! <SUBROUTINE NAME="ocean_model_end">
!
! <DESCRIPTION>
! Close down the ocean model
! </DESCRIPTION>

!> ocean_model_end terminates the model run, saving the ocean state in a restart
!! and deallocating any data associated with the ocean.
subroutine ocean_model_end(Ocean_sfc, Ocean_state, Time)
  type(ocean_public_type), intent(inout) :: Ocean_sfc   !< An ocean_public_type structure that is
                                                        !! to be deallocated upon termination.
  type(ocean_state_type),  pointer       :: Ocean_state !< A pointer to the structure containing
                                                        !! the internal ocean state to be deallocated
                                                        !! upon termination.
  type(time_type),         intent(in)    :: Time        !< The model time, used for writing restarts.

!   This subroutine terminates the model run, saving the ocean state in a
! restart file and deallocating any data associated with the ocean.

! Arguments: Ocean_sfc - An ocean_public_type structure that is to be
!                        deallocated upon termination.
!  (inout)   Ocean_state - A pointer to the structure containing the internal
!                          ocean state to be deallocated upon termination.
!  (in)      Time - The model time, used for writing restarts.

  call ocean_model_save_restart(Ocean_state, Time)
  call diag_mediator_end(Time, Ocean_state%MOM_CSp%diag)
  call MOM_end(Ocean_state%MOM_CSp)
  if (Ocean_state%use_ice_shelf) call ice_shelf_end(Ocean_state%Ice_shelf_CSp)
end subroutine ocean_model_end
! </SUBROUTINE> NAME="ocean_model_end"


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
! Arguments: Ocean_state - A structure containing the internal ocean state (in).
!  (in)      Time - The model time at this call.  This is needed for mpp_write calls.
!  (in, opt) directory - An optional directory into which to write these restart files.
!  (in, opt) filename_suffix - An optional suffix (e.g., a time-stamp) to append
!                              to the restart file names.

! Note: This is a new routine - it will need to exist for the new incremental
!   checkpointing.  It will also be called by ocean_model_end, giving the same
!   restart behavior as now in FMS.
  character(len=200) :: restart_dir

  if (OS%MOM_CSp%t_dyn_rel_adv > 0.0) call MOM_error(WARNING, "End of MOM_main reached "//&
       "with inconsistent dynamics and advective times.  Additional restart fields "//&
       "that have not been coded yet would be required for reproducibility.")
  if (.not.OS%fluxes%fluxes_used) call MOM_error(FATAL, "ocean_model_save_restart "//&
       "was called with unused buoyancy fluxes.  For conservation, the ocean "//&
       "restart files can only be created after the buoyancy forcing is applied.")

  if (present(directory)) then ; restart_dir = directory
  else ; restart_dir = OS%dirs%restart_output_dir ; endif

  call save_restart(restart_dir, Time, OS%grid, OS%MOM_CSp%restart_CSp, GV=OS%GV)

  call forcing_save_restart(OS%forcing_CSp, OS%grid, Time, restart_dir)

  if (OS%use_ice_shelf) then
    call ice_shelf_save_restart(OS%Ice_shelf_CSp, OS%Time, OS%dirs%restart_output_dir)
  endif

end subroutine ocean_model_save_restart

!=======================================================================

subroutine initialize_ocean_public_type(input_domain, Ocean_sfc, diag, maskmap, &
                                        gas_fields_ocn)
  type(domain2D), intent(in)             :: input_domain
  type(ocean_public_type), intent(inout) :: Ocean_sfc
  type(diag_ctrl), intent(in)            :: diag
  logical, intent(in), optional          :: maskmap(:,:)
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
  if(PRESENT(maskmap)) then
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

  Ocean_sfc%t_surf  = 0.0  ! time averaged sst (Kelvin) passed to atmosphere/ice model
  Ocean_sfc%s_surf  = 0.0  ! time averaged sss (psu) passed to atmosphere/ice models
  Ocean_sfc%u_surf  = 0.0  ! time averaged u-current (m/sec) passed to atmosphere/ice models
  Ocean_sfc%v_surf  = 0.0  ! time averaged v-current (m/sec)  passed to atmosphere/ice models
  Ocean_sfc%sea_lev = 0.0  ! time averaged thickness of top model grid cell (m) plus patm/rho0/grav
  Ocean_sfc%frazil  = 0.0  ! time accumulated frazil (J/m^2) passed to ice model
  Ocean_sfc%area    = 0.0
  Ocean_sfc%axes    = diag%axesT1%handles !diag axes to be used by coupler tracer flux diagnostics

  if (present(gas_fields_ocn)) then
    call coupler_type_spawn(gas_fields_ocn, Ocean_sfc%fields, (/isc,isc,iec,iec/), &
                              (/jsc,jsc,jec,jec/), suffix = '_ocn', as_needed=.true.)
  endif

end subroutine initialize_ocean_public_type

subroutine convert_state_to_ocean_type(sfc_state, Ocean_sfc, G, use_conT_absS, &
                                       patm, press_to_z)
  type(surface),           intent(inout) :: sfc_state !< A structure containing fields that
                                                      !! describe the surface state of the ocean.
  type(ocean_public_type), &
                   target, intent(inout) :: Ocean_sfc
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure
  logical,                 intent(in)    :: use_conT_absS
  real,          optional, intent(in)    :: patm(:,:)
  real,          optional, intent(in)    :: press_to_z
! This subroutine translates the coupler's ocean_data_type into MOM's
! surface state variable.  This may eventually be folded into the MOM
! code that calculates the surface state in the first place.
! Note the offset in the arrays because the ocean_data_type has no
! halo points in its arrays and always uses absolute indicies.
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
  if (use_conT_absS) then
    !If directed convert the surface T&S from conservative T to potential T and
    !from absolute (reference) salinity to practical salinity
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%s_surf(i,j) = gsw_sp_from_sr(sfc_state%SSS(i+i0,j+j0))
      Ocean_sfc%t_surf(i,j) = gsw_pt_from_ct(sfc_state%SSS(i+i0,j+j0), &
                               sfc_state%SST(i+i0,j+j0)) + CELSIUS_KELVIN_OFFSET
    enddo ; enddo
  else
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%t_surf(i,j) = sfc_state%SST(i+i0,j+j0) + CELSIUS_KELVIN_OFFSET
      Ocean_sfc%s_surf(i,j) = sfc_state%SSS(i+i0,j+j0)
    enddo ; enddo
  endif

  if (present(patm)) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%sea_lev(i,j) = sfc_state%sea_lev(i+i0,j+j0) + patm(i,j) * press_to_z
      Ocean_sfc%area(i,j) = G%areaT(i+i0,j+j0)
    enddo ; enddo
  else
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%sea_lev(i,j) = sfc_state%sea_lev(i+i0,j+j0)
      Ocean_sfc%area(i,j) = G%areaT(i+i0,j+j0)
    enddo ; enddo
  endif

  if (associated(sfc_state%frazil)) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%frazil(i,j) = sfc_state%frazil(i+i0,j+j0)
    enddo ; enddo
  endif

  if (Ocean_sfc%stagger == AGRID) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%u_surf(i,j) = G%mask2dT(i+i0,j+j0) * &
                0.5*(sfc_state%u(I+i0,j+j0)+sfc_state%u(I-1+i0,j+j0))
      Ocean_sfc%v_surf(i,j) = G%mask2dT(i+i0,j+j0) * &
                0.5*(sfc_state%v(i+i0,J+j0)+sfc_state%v(i+i0,J-1+j0))
    enddo ; enddo
  elseif (Ocean_sfc%stagger == BGRID_NE) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%u_surf(i,j) = G%mask2dBu(I+i0,J+j0) * &
                0.5*(sfc_state%u(I+i0,j+j0)+sfc_state%u(I+i0,j+j0+1))
      Ocean_sfc%v_surf(i,j) = G%mask2dBu(I+i0,J+j0) * &
                0.5*(sfc_state%v(i+i0,J+j0)+sfc_state%v(i+i0+1,J+j0))
    enddo ; enddo
  elseif (Ocean_sfc%stagger == CGRID_NE) then
    do j=jsc_bnd,jec_bnd ; do i=isc_bnd,iec_bnd
      Ocean_sfc%u_surf(i,j) = G%mask2dCu(I+i0,j+j0)*sfc_state%u(I+i0,j+j0)
      Ocean_sfc%v_surf(i,j) = G%mask2dCv(i+i0,J+j0)*sfc_state%v(i+i0,J+j0)
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


!=======================================================================
! <SUBROUTINE NAME="ocean_model_init_sfc">
!
! <DESCRIPTION>
!   This subroutine extracts the surface properties from the ocean's internal
! state and stores them in the ocean type returned to the calling ice model.
! It has to be separate from the ocean_initialization call because the coupler
! module allocates the space for some of these variables.
! </DESCRIPTION>

subroutine ocean_model_init_sfc(OS, Ocean_sfc)
  type(ocean_state_type),  pointer       :: OS
  type(ocean_public_type), intent(inout) :: Ocean_sfc

  integer :: is, ie, js, je

  is = OS%grid%isc ; ie = OS%grid%iec ; js = OS%grid%jsc ; je = OS%grid%jec
  call coupler_type_spawn(Ocean_sfc%fields, OS%sfc_state%tr_fields, &
                          (/is,is,ie,ie/), (/js,js,je,je/), as_needed=.true.)

  call calculate_surface_state(OS%sfc_state, OS%MOM_CSp%u, &
           OS%MOM_CSp%v, OS%MOM_CSp%h, OS%MOM_CSp%ave_ssh,&
           OS%grid, OS%GV, OS%MOM_CSp)

  call convert_state_to_ocean_type(OS%sfc_state, Ocean_sfc, OS%grid, &
                                   OS%MOM_CSp%use_conT_absS)

end subroutine ocean_model_init_sfc
! </SUBROUTINE NAME="ocean_model_init_sfc">

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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Ocean_stock_pe - returns stocks of heat, water, etc. for conservation checks.!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Ocean_stock_pe - returns the integrated stocks of heat, water, etc. for conservation checks.
subroutine Ocean_stock_pe(OS, index, value, time_index)
  use stock_constants_mod, only : ISTOCK_WATER, ISTOCK_HEAT,ISTOCK_SALT
  type(ocean_state_type), pointer     :: OS         !< A structure containing the internal ocean state.
                                                    !! The data in OS is intent(in).
  integer,                intent(in)  :: index      !< The stock index for the quantity of interest.
  real,                   intent(out) :: value      !< Sum returned for the conservation quantity of interest.
  integer,      optional, intent(in)  :: time_index !< An unused optional argument, present only for
                                                    !! interfacial compatibility with other models.
! Arguments: OS - A structure containing the internal ocean state.
!  (in)      index - Index of conservation quantity of interest.
!  (in)      value -  Sum returned for the conservation quantity of interest.
!  (in,opt)  time_index - Index for time level to use if this is necessary.

  real :: to_heat, to_mass, to_salt, PSU_to_kg ! Conversion constants.
  integer :: i, j, k, is, ie, js, je, nz, ind

  value = 0.0
  if (.not.associated(OS)) return
  if (.not.OS%is_ocean_pe) return

  is = OS%grid%isc ; ie = OS%grid%iec
  js = OS%grid%jsc ; je = OS%grid%jec ; nz = OS%grid%ke

  select case (index)
    case (ISTOCK_WATER)
      ! Return the mass of fresh water in the ocean on this PE in kg.
      to_mass = OS%GV%H_to_kg_m2
      if (OS%GV%Boussinesq) then
        do k=1,nz ; do j=js,je ; do i=is,ie ; if (OS%grid%mask2dT(i,j) > 0.5) then
          value = value + to_mass*(OS%MOM_CSp%h(i,j,k) * OS%grid%areaT(i,j))
        endif ; enddo ; enddo ; enddo
      else
        ! In non-Boussinesq mode, the mass of salt needs to be subtracted.
        PSU_to_kg = 1.0e-3
        do k=1,nz ; do j=js,je ; do i=is,ie ; if (OS%grid%mask2dT(i,j) > 0.5) then
          value = value + to_mass * ((1.0 - PSU_to_kg*OS%MOM_CSp%tv%S(i,j,k))*&
                                  (OS%MOM_CSp%h(i,j,k) * OS%grid%areaT(i,j)))
        endif ; enddo ; enddo ; enddo
      endif
    case (ISTOCK_HEAT)
      ! Return the heat content of the ocean on this PE in J.
      to_heat = OS%GV%H_to_kg_m2 * OS%C_p
      do k=1,nz ; do j=js,je ; do i=is,ie ; if (OS%grid%mask2dT(i,j) > 0.5) then
        value = value + (to_heat * OS%MOM_CSp%tv%T(i,j,k)) * &
                        (OS%MOM_CSp%h(i,j,k)*OS%grid%areaT(i,j))
      endif ; enddo ; enddo ; enddo
    case (ISTOCK_SALT)
      ! Return the mass of the salt in the ocean on this PE in kg.
      ! The 1000 converts salinity in PSU to salt in kg kg-1.
      to_salt = OS%GV%H_to_kg_m2 / 1000.0
      do k=1,nz ; do j=js,je ; do i=is,ie ; if (OS%grid%mask2dT(i,j) > 0.5) then
        value = value + (to_salt * OS%MOM_CSp%tv%S(i,j,k)) * &
                        (OS%MOM_CSp%h(i,j,k)*OS%grid%areaT(i,j))
      endif ; enddo ; enddo ; enddo
    case default ; value = 0.0
  end select

end subroutine Ocean_stock_pe

subroutine ocean_model_data2D_get(OS,Ocean, name, array2D,isc,jsc)
  use MOM_constants, only : CELSIUS_KELVIN_OFFSET
  type(ocean_state_type),     pointer    :: OS
  type(ocean_public_type),    intent(in) :: Ocean
  character(len=*)          , intent(in) :: name
  real, dimension(isc:,jsc:), intent(out):: array2D
  integer                   , intent(in) :: isc,jsc

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
     array2D(isc:,jsc:) = OS%grid%areaT(g_isc:g_iec,g_jsc:g_jec)
  case('mask')
     array2D(isc:,jsc:) = OS%grid%mask2dT(g_isc:g_iec,g_jsc:g_jec)
!OR same result
!     do j=g_jsc,g_jec; do i=g_isc,g_iec
!        array2D(isc+i-g_isc,jsc+j-g_jsc) = OS%grid%mask2dT(i,j)
!     enddo; enddo
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
  case default
     call MOM_error(FATAL,'get_ocean_grid_data2D: unknown argument name='//name)
  end select


end subroutine ocean_model_data2D_get

subroutine ocean_model_data1D_get(OS,Ocean, name, value)
  type(ocean_state_type),     pointer    :: OS
  type(ocean_public_type),    intent(in) :: Ocean
  character(len=*)          , intent(in) :: name
  real                      , intent(out):: value

  if (.not.associated(OS)) return
  if (.not.OS%is_ocean_pe) return

  select case(name)
  case('c_p')
     value = OS%C_p
  case default
     call MOM_error(FATAL,'get_ocean_grid_data1D: unknown argument name='//name)
  end select


end subroutine ocean_model_data1D_get

subroutine ocean_public_type_chksum(id, timestep, ocn)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ocean_public_type), intent(in) :: ocn
    integer ::   n,m, outunit

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

end module ocean_model_mod
