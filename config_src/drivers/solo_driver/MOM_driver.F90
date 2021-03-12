program MOM_main

! This file is part of MOM6. See LICENSE.md for the license.

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*                  The Modular Ocean Model, version 6                 *
!*                               MOM6                                  *
!*                                                                     *
!*  By Alistair Adcroft, Stephen Griffies and Robert Hallberg          *
!*                                                                     *
!*    This file is the ocean-only driver for Version 6 of the Modular  *
!*  Ocean Model (MOM).  A separate ocean interface for use with        *
!*  coupled models is provided in ocean_model_MOM.F90.   These two     *
!*  drivers are kept in separate directories for convenience of code   *
!*  selection during compiling.  This file orchestrates the calls to   *
!*  the MOM initialization routines, to the subroutine that steps      *
!*  the model, and coordinates the output and saving restarts.  A      *
!*  description of all of the files that constitute MOM is found in    *
!*  the comments at the beginning of MOM.F90.  The arguments of each   *
!*  subroutine are described where the subroutine is defined.          *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

  use MOM_cpu_clock,       only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
  use MOM_cpu_clock,       only : CLOCK_COMPONENT
  use MOM_diag_mediator,   only : enable_averaging, disable_averaging, diag_mediator_end
  use MOM_diag_mediator,   only : diag_ctrl, diag_mediator_close_registration
  use MOM,                 only : initialize_MOM, step_MOM, MOM_control_struct, MOM_end
  use MOM,                 only : extract_surface_state, finish_MOM_initialization
  use MOM,                 only : get_MOM_state_elements, MOM_state_is_synchronized
  use MOM,                 only : step_offline
  use MOM_coms,            only : Set_PElist
  use MOM_domains,         only : MOM_infra_init, MOM_infra_end, set_MOM_thread_affinity
  use MOM_ensemble_manager, only : ensemble_manager_init, get_ensemble_size
  use MOM_ensemble_manager, only : ensemble_pelist_setup
  use MOM_error_handler,   only : MOM_error, MOM_mesg, WARNING, FATAL, is_root_pe
  use MOM_error_handler,   only : callTree_enter, callTree_leave, callTree_waypoint
  use MOM_file_parser,     only : read_param, get_param, log_param, log_version, param_file_type
  use MOM_file_parser,     only : close_param_file
  use MOM_forcing_type,    only : forcing, mech_forcing, forcing_diagnostics
  use MOM_forcing_type,    only : mech_forcing_diags, MOM_forcing_chksum, MOM_mech_forcing_chksum
  use MOM_get_input,       only : get_MOM_input, directories
  use MOM_grid,            only : ocean_grid_type
  use MOM_ice_shelf,       only : initialize_ice_shelf, ice_shelf_end, ice_shelf_CS
  use MOM_ice_shelf,       only : shelf_calc_flux, add_shelf_forces, ice_shelf_save_restart
  use MOM_ice_shelf,       only : initialize_ice_shelf_fluxes, initialize_ice_shelf_forces
  use MOM_interpolate,     only : time_interp_external_init
  use MOM_io,              only : file_exists, open_ASCII_file, close_file
  use MOM_io,              only : check_nml_error, io_infra_init, io_infra_end
  use MOM_io,              only : APPEND_FILE, READONLY_FILE
  use MOM_restart,         only : MOM_restart_CS, save_restart
  use MOM_string_functions,only : uppercase
  use MOM_surface_forcing, only : set_forcing, forcing_save_restart
  use MOM_surface_forcing, only : surface_forcing_init, surface_forcing_CS
  use MOM_time_manager,    only : time_type, set_date, get_date, real_to_time, time_type_to_real
  use MOM_time_manager,    only : operator(+), operator(-), operator(*), operator(/)
  use MOM_time_manager,    only : operator(>), operator(<), operator(>=)
  use MOM_time_manager,    only : increment_date, set_calendar_type, month_name
  use MOM_time_manager,    only : JULIAN, GREGORIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
  use MOM_tracer_flow_control, only : tracer_flow_control_CS
  use MOM_unit_scaling,    only : unit_scale_type
  use MOM_variables,       only : surface
  use MOM_verticalGrid,    only : verticalGrid_type
  use MOM_wave_interface,  only : wave_parameters_CS, MOM_wave_interface_init
  use MOM_wave_interface,  only : MOM_wave_interface_init_lite, Update_Surface_Waves
  use MOM_write_cputime,   only : write_cputime, MOM_write_cputime_init
  use MOM_write_cputime,   only : write_cputime_start_clock, write_cputime_CS

  implicit none

#include <MOM_memory.h>

  ! A structure with the driving mechanical surface forces
  type(mech_forcing) :: forces
  ! A structure containing pointers to the thermodynamic forcing fields
  ! at the ocean surface.
  type(forcing) :: fluxes
  ! A structure containing pointers to the ocean surface state fields.
  type(surface) :: sfc_state

  ! A pointer to a structure containing metrics and related information.
  type(ocean_grid_type), pointer :: grid => NULL()
  type(verticalGrid_type), pointer :: GV => NULL()
  ! A pointer to a structure containing dimensional unit scaling factors.
  type(unit_scale_type), pointer :: US => NULL()

  ! If .true., use the ice shelf model for part of the domain.
  logical :: use_ice_shelf = .false.

  ! If .true., use surface wave coupling
  logical :: use_waves = .false.

  ! This is .true. if incremental restart files may be saved.
  logical :: permit_incr_restart = .true.

  integer :: ns

  ! nmax is the number of iterations after which to stop so that the
  ! simulation does not exceed its CPU time limit.  nmax is determined by
  ! evaluating the CPU time used between successive calls to write_cputime.
  ! Initially it is set to be very large.
  integer :: nmax=2000000000

  ! A structure containing several relevant directory paths.
  type(directories) :: dirs

  ! A suite of time types for use by MOM
  type(time_type), target :: Time       ! A copy of the ocean model's time.
                                        ! Other modules can set pointers to this and
                                        ! change it to manage diagnostics.
  type(time_type) :: Master_Time        ! The ocean model's master clock. No other
                                        ! modules are ever given access to this.
  type(time_type) :: Time1              ! The value of the ocean model's time at the
                                        ! start of a call to step_MOM.
  type(time_type) :: Start_time         ! The start time of the simulation.
  type(time_type) :: segment_start_time ! The start time of this run segment.
  type(time_type) :: Time_end           ! End time for the segment or experiment.
  type(time_type) :: restart_time       ! The next time to write restart files.
  type(time_type) :: Time_step_ocean    ! A time_type version of dt_forcing.

  real    :: elapsed_time = 0.0   ! Elapsed time in this run  [s].
  logical :: elapsed_time_master  ! If true, elapsed time is used to set the
                                  ! model's master clock (Time).  This is needed
                                  ! if Time_step_ocean is not an exact
                                  ! representation of dt_forcing.
  real :: dt_forcing              ! The coupling time step [s].
  real :: dt                      ! The nominal baroclinic dynamics time step [s].
  real :: dt_off                  ! Offline time step [s].
  integer :: ntstep               ! The number of baroclinic dynamics time steps
                                  ! within dt_forcing.
  real :: dt_therm                ! The thermodynamic timestep [s]
  real :: dt_dyn                  ! The actual dynamic timestep used [s].  The value of dt_dyn is
                                  ! chosen so that dt_forcing is an integer multiple of dt_dyn.
  real :: dtdia                   ! The diabatic timestep [s]
  real :: t_elapsed_seg           ! The elapsed time in this run segment [s]
  integer :: n, n_max, nts, n_last_thermo
  logical :: diabatic_first, single_step_call
  type(time_type) :: Time2, time_chg

  integer :: Restart_control    ! An integer that is bit-tested to determine whether
                                ! incremental restart files are saved and whether they
                                ! have a time stamped name.  +1 (bit 0) for generic
                                ! files and +2 (bit 1) for time-stamped files.  A
                                ! restart file is saved at the end of a run segment
                                ! unless Restart_control is negative.

  real            :: Time_unit       ! The time unit for the following input fields [s].
  type(time_type) :: restint         ! The time between saves of the restart file.
  type(time_type) :: daymax          ! The final day of the simulation.

  integer :: CPU_steps          ! The number of steps between writing CPU time.
  integer :: date_init(6)=0                ! The start date of the whole simulation.
  integer :: date(6)=-1                    ! Possibly the start date of this run segment.
  integer :: years=0, months=0, days=0     ! These may determine the segment run
  integer :: hours=0, minutes=0, seconds=0 ! length, if read from a namelist.
  integer :: yr, mon, day, hr, mins, sec   ! Temp variables for writing the date.
  type(param_file_type) :: param_file      ! The structure indicating the file(s)
                                           ! containing all run-time parameters.
  character(len=9)  :: month
  character(len=16) :: calendar = 'julian'
  integer :: calendar_type=-1

  integer :: unit, io_status, ierr
  integer :: ensemble_size, nPEs_per, ensemble_info(6)

  integer, dimension(0) :: atm_PElist, land_PElist, ice_PElist
  integer, dimension(:), allocatable :: ocean_PElist
  logical :: unit_in_use
  integer :: initClock, mainClock, termClock

  logical :: debug               ! If true, write verbose checksums for debugging purposes.
  logical :: offline_tracer_mode ! If false, use the model in prognostic mode where
                                 ! the barotropic and baroclinic dynamics, thermodynamics,
                                 ! etc. are stepped forward integrated in time.
                                 ! If true, then all of the above are bypassed with all
                                 ! fields necessary to integrate only the tracer advection
                                 ! and diffusion equation are read in from files stored from
                                 ! a previous integration of the prognostic model

  type(MOM_control_struct),  pointer :: MOM_CSp => NULL()
  !> A pointer to the tracer flow control structure.
  type(tracer_flow_control_CS), pointer :: &
    tracer_flow_CSp => NULL()  !< A pointer to the tracer flow control structure
  type(surface_forcing_CS),  pointer :: surface_forcing_CSp => NULL()
  type(write_cputime_CS),    pointer :: write_CPU_CSp => NULL()
  type(ice_shelf_CS),        pointer :: ice_shelf_CSp => NULL()
  type(wave_parameters_cs),  pointer :: waves_CSp => NULL()
  type(MOM_restart_CS),      pointer :: &
    restart_CSp => NULL()     !< A pointer to the restart control structure
                              !! that will be used for MOM restart files.
  type(diag_ctrl),           pointer :: &
       diag => NULL()         !< A pointer to the diagnostic regulatory structure
  !-----------------------------------------------------------------------

  character(len=4), parameter :: vers_num = 'v2.0'
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod_name = "MOM_main (MOM_driver)" ! This module's name.

  integer :: ocean_nthreads = 1
  logical :: use_hyper_thread = .false.
  integer :: omp_get_num_threads,omp_get_thread_num
  namelist /ocean_solo_nml/ date_init, calendar, months, days, hours, minutes, seconds,&
                            ocean_nthreads, use_hyper_thread

  !=====================================================================

  call write_cputime_start_clock(write_CPU_CSp)

  call MOM_infra_init() ; call io_infra_init()

  !allocate(forces,fluxes,sfc_state)

  ! Initialize the ensemble manager.  If there are no settings for ensemble_size
  ! in input.nml(ensemble.nml), these should not do anything.  In coupled
  ! configurations, this all occurs in the external driver.
  call ensemble_manager_init() ; ensemble_info(:) =  get_ensemble_size()
  ensemble_size=ensemble_info(1) ; nPEs_per=ensemble_info(2)
  if (ensemble_size > 1) then ! There are multiple ensemble members.
    allocate(ocean_pelist(nPEs_per))
    call ensemble_pelist_setup(.true., 0, nPEs_per, 0, 0, atm_pelist, ocean_pelist, &
                               land_pelist, ice_pelist)
    call Set_PElist(ocean_pelist)
    deallocate(ocean_pelist)
  endif

  ! These clocks are on the global pelist.
  initClock = cpu_clock_id( 'Initialization' )
  mainClock = cpu_clock_id( 'Main loop' )
  termClock = cpu_clock_id( 'Termination' )
  call cpu_clock_begin(initClock)

  call MOM_mesg('======== Model being driven by MOM_driver ========', 2)
  call callTree_waypoint("Program MOM_main, MOM_driver.F90")

  if (file_exists('input.nml')) then
    ! Provide for namelist specification of the run length and calendar data.
    call open_ASCII_file(unit, 'input.nml', action=READONLY_FILE)
    read(unit, ocean_solo_nml, iostat=io_status)
    call close_file(unit)
    ierr = check_nml_error(io_status,'ocean_solo_nml')
    if (years+months+days+hours+minutes+seconds > 0) then
      if (is_root_pe()) write(*,ocean_solo_nml)
    endif
  endif

  ! This call sets the number and affinity of threads with openMP.
  !$  call set_MOM_thread_affinity(ocean_nthreads, use_hyper_thread)

  ! This call is required to initiate dirs%restart_input_dir for ocean_solo.res
  ! The contents of dirs will be reread in initialize_MOM.
  call get_MOM_input(dirs=dirs)

  ! Read ocean_solo restart, which can override settings from the namelist.
  if (file_exists(trim(dirs%restart_input_dir)//'ocean_solo.res')) then
    call open_ASCII_file(unit, trim(dirs%restart_input_dir)//'ocean_solo.res', action=READONLY_FILE)
    read(unit,*) calendar_type
    read(unit,*) date_init
    read(unit,*) date
    call close_file(unit)
  else
    calendar = uppercase(calendar)
    if (calendar(1:6) == 'JULIAN') then ;        calendar_type = JULIAN
    elseif (calendar(1:9) == 'GREGORIAN') then ; calendar_type = GREGORIAN
    elseif (calendar(1:6) == 'NOLEAP') then ;    calendar_type = NOLEAP
    elseif (calendar(1:10)=='THIRTY_DAY') then ; calendar_type = THIRTY_DAY_MONTHS
    elseif (calendar(1:11)=='NO_CALENDAR') then; calendar_type = NO_CALENDAR
    elseif (calendar(1:1) /= ' ') then
      call MOM_error(FATAL,'MOM_driver: Invalid namelist value '//trim(calendar)//' for calendar')
    else
      call MOM_error(FATAL,'MOM_driver: No namelist value for calendar')
    endif
  endif
  call set_calendar_type(calendar_type)


  if (sum(date_init) > 0) then
    Start_time = set_date(date_init(1), date_init(2), date_init(3), &
                          date_init(4), date_init(5), date_init(6))
  else
    Start_time = real_to_time(0.0)
  endif

  call time_interp_external_init()

  if (sum(date) >= 0) then
    ! In this case, the segment starts at a time fixed by ocean_solo.res
    segment_start_time = set_date(date(1), date(2), date(3), date(4), date(5), date(6))
    Time = segment_start_time
  else
    ! In this case, the segment starts at a time read from the MOM restart file
    ! or left as Start_time by MOM_initialize.
    Time = Start_time
  endif

  ! Call initialize MOM with an optional Ice Shelf CS which, if present triggers
  ! initialization of ice shelf parameters and arrays.
  if (sum(date) >= 0) then
    call initialize_MOM(Time, Start_time, param_file, dirs, MOM_CSp, restart_CSp, &
                        segment_start_time, offline_tracer_mode=offline_tracer_mode, &
                        diag_ptr=diag, tracer_flow_CSp=tracer_flow_CSp, ice_shelf_CSp=ice_shelf_CSp)
  else
    call initialize_MOM(Time, Start_time, param_file, dirs, MOM_CSp, restart_CSp, &
                        offline_tracer_mode=offline_tracer_mode, diag_ptr=diag, &
                        tracer_flow_CSp=tracer_flow_CSp, ice_shelf_CSp=ice_shelf_CSp)
  endif

  call get_MOM_state_elements(MOM_CSp, G=grid, GV=GV, US=US, C_p_scaled=fluxes%C_p)
  Master_Time = Time
  use_ice_shelf = associated(ice_shelf_CSp)

  if (use_ice_shelf) then
    ! These arrays are not initialized in most solo cases, but are needed
    ! when using an ice shelf
    call initialize_ice_shelf_fluxes(ice_shelf_CSp, grid, US, fluxes)
    call initialize_ice_shelf_forces(ice_shelf_CSp, grid, US, forces)
  endif


  call callTree_waypoint("done initialize_MOM")

  call extract_surface_state(MOM_CSp, sfc_state)

  call surface_forcing_init(Time, grid, US, param_file, diag, &
                            surface_forcing_CSp, tracer_flow_CSp)
  call callTree_waypoint("done surface_forcing_init")


  call get_param(param_file,mod_name, "USE_WAVES", Use_Waves, &
       "If true, enables surface wave modules.",default=.false.)
  if (use_waves) then
    call MOM_wave_interface_init(Time, grid, GV, US, param_file, Waves_CSp, diag)
  else
    call MOM_wave_interface_init_lite(param_file)
  endif

  segment_start_time = Time
  elapsed_time = 0.0

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod_name, version, "")
  call get_param(param_file, mod_name, "DT", dt, fail_if_missing=.true.)
  call get_param(param_file, mod_name, "DT_FORCING", dt_forcing, &
                 "The time step for changing forcing, coupling with other "//&
                 "components, or potentially writing certain diagnostics. "//&
                 "The default value is given by DT.", units="s", default=dt)
  if (offline_tracer_mode) then
    call get_param(param_file, mod_name, "DT_OFFLINE", dt_forcing, &
                   "Time step for the offline time step")
    dt = dt_forcing
  endif
  ntstep = MAX(1,ceiling(dt_forcing/dt - 0.001))

  Time_step_ocean = real_to_time(dt_forcing)
  elapsed_time_master = (abs(dt_forcing - time_type_to_real(Time_step_ocean)) > 1.0e-12*dt_forcing)
  if (elapsed_time_master) &
    call MOM_mesg("Using real elapsed time for the master clock.", 2)

  ! Determine the segment end time, either from the namelist file or parsed input file.
  call get_param(param_file, mod_name, "TIMEUNIT", Time_unit, &
                 "The time unit for DAYMAX, ENERGYSAVEDAYS, and RESTINT.", &
                 units="s", default=86400.0)
  if (years+months+days+hours+minutes+seconds > 0) then
    Time_end = increment_date(Time, years, months, days, hours, minutes, seconds)
    call MOM_mesg('Segment run length determined from ocean_solo_nml.', 2)
    call get_param(param_file, mod_name, "DAYMAX", daymax, timeunit=Time_unit, &
                   default=Time_end, do_not_log=.true.)
    call log_param(param_file, mod_name, "DAYMAX", daymax, &
                 "The final time of the whole simulation, in units of "//&
                 "TIMEUNIT seconds.  This also sets the potential end "//&
                 "time of the present run segment if the end time is "//&
                 "not set via ocean_solo_nml in input.nml.", &
                 timeunit=Time_unit)
  else
    call get_param(param_file, mod_name, "DAYMAX", daymax, &
                 "The final time of the whole simulation, in units of "//&
                 "TIMEUNIT seconds.  This also sets the potential end "//&
                 "time of the present run segment if the end time is "//&
                 "not set via ocean_solo_nml in input.nml.", &
                 timeunit=Time_unit, fail_if_missing=.true.)
    Time_end = daymax
  endif

  call get_param(param_file, mod_name, "SINGLE_STEPPING_CALL", single_step_call, &
                 "If true, advance the state of MOM with a single step "//&
                 "including both dynamics and thermodynamics.  If false "//&
                 "the two phases are advanced with separate calls.", default=.true.)
  call get_param(param_file, mod_name, "DT_THERM", dt_therm, &
                 "The thermodynamic and tracer advection time step. "//&
                 "Ideally DT_THERM should be an integer multiple of DT "//&
                 "and less than the forcing or coupling time-step, unless "//&
                 "THERMO_SPANS_COUPLING is true, in which case DT_THERM "//&
                 "can be an integer multiple of the coupling timestep.  By "//&
                 "default DT_THERM is set to DT.", units="s", default=dt)
  call get_param(param_file, mod_name, "DIABATIC_FIRST", diabatic_first, &
                 "If true, apply diabatic and thermodynamic processes, "//&
                 "including buoyancy forcing and mass gain or loss, "//&
                 "before stepping the dynamics forward.", default=.false.)


  if (Time >= Time_end) call MOM_error(FATAL, &
    "MOM_driver: The run has been started at or after the end time of the run.")

  call get_param(param_file, mod_name, "RESTART_CONTROL", Restart_control, &
                 "An integer whose bits encode which restart files are "//&
                 "written. Add 2 (bit 1) for a time-stamped file, and odd "//&
                 "(bit 0) for a non-time-stamped file. A non-time-stamped "//&
                 "restart file is saved at the end of the run segment "//&
                 "for any non-negative value.", default=1)
  call get_param(param_file, mod_name, "RESTINT", restint, &
                 "The interval between saves of the restart file in units "//&
                 "of TIMEUNIT.  Use 0 (the default) to not save "//&
                 "incremental restart files at all.", default=real_to_time(0.0), &
                 timeunit=Time_unit)
  call get_param(param_file, mod_name, "WRITE_CPU_STEPS", cpu_steps, &
                 "The number of coupled timesteps between writing the cpu "//&
                 "time. If this is not positive, do not check cpu time, and "//&
                 "the segment run-length can not be set via an elapsed CPU time.", &
                 default=1000)
  call get_param(param_file, "MOM", "DEBUG", debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)

  call log_param(param_file, mod_name, "ELAPSED TIME AS MASTER", elapsed_time_master)

  if (cpu_steps > 0) &
    call MOM_write_cputime_init(param_file, dirs%output_directory, Start_time, &
                                write_CPU_CSp)

  ! Close the param_file.  No further parsing of input is possible after this.
  call close_param_file(param_file)
  call diag_mediator_close_registration(diag)

  ! Write out a time stamp file.
  if (is_root_pe() .and. (calendar_type /= NO_CALENDAR)) then
    call open_ASCII_file(unit, 'time_stamp.out', action=APPEND_FILE)
    call get_date(Time, date(1), date(2), date(3), date(4), date(5), date(6))
    month = month_name(date(2))
    write(unit,'(6i4,2x,a3)') date, month(1:3)
    call get_date(Time_end, date(1), date(2), date(3), date(4), date(5), date(6))
    month = month_name(date(2))
    write(unit,'(6i4,2x,a3)') date, month(1:3)
    call close_file(unit)
  endif

  if (cpu_steps > 0) call write_cputime(Time, 0, write_CPU_CSp)

  if (((.not.BTEST(Restart_control,1)) .and. (.not.BTEST(Restart_control,0))) &
      .or. (Restart_control < 0)) permit_incr_restart = .false.

  if (restint > real_to_time(0.0)) then
    ! restart_time is the next integral multiple of restint.
    restart_time = Start_time + restint * &
        (1 + ((Time + Time_step_ocean) - Start_time) / restint)
  else
    ! Set the time so late that there is no intermediate restart.
    restart_time = Time_end + Time_step_ocean
    permit_incr_restart = .false.
  endif

  call cpu_clock_end(initClock) !end initialization

  call cpu_clock_begin(mainClock) !begin main loop

  ns = 1
  do while ((ns < nmax) .and. (Time < Time_end))
    call callTree_enter("Main loop, MOM_driver.F90",ns)

    ! Set the forcing for the next steps.
    if (.not. offline_tracer_mode) then
        call set_forcing(sfc_state, forces, fluxes, Time, Time_step_ocean, grid, US, &
                     surface_forcing_CSp)
    endif
    if (debug) then
      call MOM_mech_forcing_chksum("After set forcing", forces, grid, US, haloshift=0)
      call MOM_forcing_chksum("After set forcing", fluxes, grid, US, haloshift=0)
    endif

    if (use_ice_shelf) then
      call shelf_calc_flux(sfc_state, fluxes, Time, dt_forcing, ice_shelf_CSp)
      call add_shelf_forces(grid, US, Ice_shelf_CSp, forces, external_call=.true.)
    endif
    fluxes%fluxes_used = .false.
    fluxes%dt_buoy_accum = US%s_to_T*dt_forcing

    if (use_waves) then
      call Update_Surface_Waves(grid, GV, US, time, time_step_ocean, waves_csp)
    endif

    if (ns==1) then
      call finish_MOM_initialization(Time, dirs, MOM_CSp, restart_CSp)
    endif

    ! This call steps the model over a time dt_forcing.
    Time1 = Master_Time ; Time = Master_Time
    if (offline_tracer_mode) then
      call step_offline(forces, fluxes, sfc_state, Time1, dt_forcing, MOM_CSp)
    elseif (single_step_call) then
      call step_MOM(forces, fluxes, sfc_state, Time1, dt_forcing, MOM_CSp, Waves=Waves_CSP)
    else
      n_max = 1 ; if (dt_forcing > dt) n_max = ceiling(dt_forcing/dt - 0.001)
      dt_dyn = dt_forcing / real(n_max)

      nts = MAX(1,MIN(n_max,floor(dt_therm/dt_dyn + 0.001)))
      n_last_thermo = 0

      Time2 = Time1 ; t_elapsed_seg = 0.0
      do n=1,n_max
        if (diabatic_first) then
          if (modulo(n-1,nts)==0) then
            dtdia = dt_dyn*min(ntstep,n_max-(n-1))
            call step_MOM(forces, fluxes, sfc_state, Time2, dtdia, MOM_CSp, &
                          do_dynamics=.false., do_thermodynamics=.true., &
                          start_cycle=(n==1), end_cycle=.false., cycle_length=dt_forcing)
          endif

          call step_MOM(forces, fluxes, sfc_state, Time2, dt_dyn, MOM_CSp, &
                        do_dynamics=.true., do_thermodynamics=.false., &
                        start_cycle=.false., end_cycle=(n==n_max), cycle_length=dt_forcing)
        else
          call step_MOM(forces, fluxes, sfc_state, Time2, dt_dyn, MOM_CSp, &
                        do_dynamics=.true., do_thermodynamics=.false., &
                        start_cycle=(n==1), end_cycle=.false., cycle_length=dt_forcing)

          if ((modulo(n,nts)==0) .or. (n==n_max)) then
            dtdia = dt_dyn*(n - n_last_thermo)
            ! Back up Time2 to the start of the thermodynamic segment.
            if (n > n_last_thermo+1) &
              Time2 = Time2 - real_to_time(dtdia - dt_dyn)
            call step_MOM(forces, fluxes, sfc_state, Time2, dtdia, MOM_CSp, &
                          do_dynamics=.false., do_thermodynamics=.true., &
                          start_cycle=.false., end_cycle=(n==n_max), cycle_length=dt_forcing)
            n_last_thermo = n
          endif
        endif

        t_elapsed_seg = t_elapsed_seg + dt_dyn
        Time2 = Time1 + real_to_time(t_elapsed_seg)
      enddo
    endif

!   Time = Time + Time_step_ocean
!   This is here to enable fractional-second time steps.
    elapsed_time = elapsed_time + dt_forcing
    if (elapsed_time > 2e9) then
      ! This is here to ensure that the conversion from a real to an integer can be accurately
      ! represented in long runs (longer than ~63 years). It will also ensure that elapsed time
      ! does not lose resolution of order the timetype's resolution, provided that the timestep and
      ! tick are larger than 10-5 seconds.  If a clock with a finer resolution is used, a smaller
      ! value would be required.
      time_chg = real_to_time(elapsed_time)
      segment_start_time = segment_start_time + time_chg
      elapsed_time = elapsed_time - time_type_to_real(time_chg)
    endif
    if (elapsed_time_master) then
      Master_Time = segment_start_time + real_to_time(elapsed_time)
    else
      Master_Time = Master_Time + Time_step_ocean
    endif
    Time = Master_Time

    if (cpu_steps > 0) then ; if (MOD(ns, cpu_steps) == 0) then
      call write_cputime(Time, ns+ntstep-1, write_CPU_CSp, nmax)
    endif ; endif

    call mech_forcing_diags(forces, dt_forcing, grid, Time, diag, surface_forcing_CSp%handles)

    if (.not. offline_tracer_mode) then
      if (fluxes%fluxes_used) then
        call forcing_diagnostics(fluxes, sfc_state, grid, US, Time, &
                                 diag, surface_forcing_CSp%handles)
      else
        call MOM_error(FATAL, "The solo MOM_driver is not yet set up to handle "//&
               "thermodynamic time steps that are longer than the coupling timestep.")
      endif
    endif

!  See if it is time to write out a restart file - timestamped or not.
    if ((permit_incr_restart) .and. (fluxes%fluxes_used) .and. &
        (Time + (Time_step_ocean/2) > restart_time)) then
      if (BTEST(Restart_control,1)) then
        call save_restart(dirs%restart_output_dir, Time, grid, &
                          restart_CSp, .true., GV=GV)
        call forcing_save_restart(surface_forcing_CSp, grid, Time, &
                            dirs%restart_output_dir, .true.)
        if (use_ice_shelf) call ice_shelf_save_restart(ice_shelf_CSp, Time, &
                                    dirs%restart_output_dir, .true.)
      endif
      if (BTEST(Restart_control,0)) then
        call save_restart(dirs%restart_output_dir, Time, grid, &
                          restart_CSp, GV=GV)
        call forcing_save_restart(surface_forcing_CSp, grid, Time, &
                            dirs%restart_output_dir)
        if (use_ice_shelf) call ice_shelf_save_restart(ice_shelf_CSp, Time, &
                                    dirs%restart_output_dir)
      endif
      restart_time = restart_time + restint
    endif

    ns = ns + ntstep
    call callTree_leave("Main loop")
  enddo

  call cpu_clock_end(mainClock)
  call cpu_clock_begin(termClock)
  if (Restart_control>=0) then
    if (.not.MOM_state_is_synchronized(MOM_CSp)) &
      call MOM_error(WARNING, "End of MOM_main reached with inconsistent "//&
         "dynamics and advective times.  Additional restart fields "//&
         "that have not been coded yet would be required for reproducibility.")
    if (.not.fluxes%fluxes_used .and. .not.offline_tracer_mode) call MOM_error(FATAL, &
         "End of MOM_main reached with unused buoyancy fluxes. "//&
         "For conservation, the ocean restart files can only be "//&
         "created after the buoyancy forcing is applied.")

    call save_restart(dirs%restart_output_dir, Time, grid, restart_CSp, GV=GV)
    if (use_ice_shelf) call ice_shelf_save_restart(ice_shelf_CSp, Time, &
                                dirs%restart_output_dir)
    ! Write ocean solo restart file.
    if (is_root_pe()) then
      call open_ASCII_file(unit, trim(dirs%restart_output_dir)//'ocean_solo.res')
      write(unit, '(i6,8x,a)') calendar_type, &
            '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'

      call get_date(Start_time, yr, mon, day, hr, mins, sec)
      write(unit, '(6i6,8x,a)') yr, mon, day, hr, mins, sec, &
            'Model start time:   year, month, day, hour, minute, second'
      call get_date(Time, yr, mon, day, hr, mins, sec)
      write(unit, '(6i6,8x,a)') yr, mon, day, hr, mins, sec, &
            'Current model time: year, month, day, hour, minute, second'
      call close_file(unit)
    endif
  endif

  if (is_root_pe()) then
    do unit=10,1967
      INQUIRE(unit,OPENED=unit_in_use)
      if (.not.unit_in_use) exit
    enddo
    open(unit,FILE="exitcode",FORM="FORMATTED",STATUS="REPLACE",action="WRITE")
    if (Time < daymax) then
      write(unit,*) 9
    else
      write(unit,*) 0
    endif
    close(unit)
  endif

  call callTree_waypoint("End MOM_main")
  if (use_ice_shelf) call ice_shelf_end(ice_shelf_CSp)
  call diag_mediator_end(Time, diag, end_diag_manager=.true.)
  if (cpu_steps > 0) call write_cputime(Time, ns-1, write_CPU_CSp, call_end=.true.)
  call cpu_clock_end(termClock)

  call io_infra_end ; call MOM_infra_end

  call MOM_end(MOM_CSp)

end program MOM_main
