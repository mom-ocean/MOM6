program Shelf_main

! This file is part of MOM6. See LICENSE.md for the license.

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Daniel Goldberg, Olga Sergienko, and Robert Hallberg            *
!*                                                                     *
!*    This file is the driver for the stand-alone ice-sheet model that *
!*  is under development at GFDL.  When used in a mode that is coupled *
!*  with an ocean model or a full coupled model, a different driver    *
!*  will be used. This file orchestrates the calls to the appropriate  *
!*  initialization routines, to the subroutine that steps the model,   *
!*  and coordinates the saving and reading of restarts.                *
!*  A description of all of the files that constitute this ice shelf   *
!*  component is found in the comments at the beginning of             *
!*  MOM_ice_shelf.F90.  The arguments of each subroutine are           *
!*  described where the subroutine is defined.                         *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h  *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

  use MOM_cpu_clock,       only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
  use MOM_cpu_clock,       only : CLOCK_COMPONENT
  use MOM_debugging,       only : MOM_debugging_init
  use MOM_diag_mediator,   only : diag_mediator_init, diag_mediator_infrastructure_init
  use MOM_diag_mediator,   only : enable_averaging, disable_averaging, diag_mediator_end
  use MOM_diag_mediator,   only : diag_ctrl, diag_mediator_close_registration
  use MOM_domains,         only : MOM_infra_init, MOM_infra_end
  use MOM_domains,         only : MOM_domains_init, clone_MOM_domain, pass_var
  use MOM_dyn_horgrid,     only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
  use MOM_error_handler,   only : MOM_error, MOM_mesg, WARNING, FATAL, is_root_pe
  use MOM_error_handler,   only : callTree_enter, callTree_leave, callTree_waypoint
  use MOM_file_parser,     only : read_param, get_param, log_param, log_version, param_file_type
  use MOM_file_parser,     only : close_param_file
  use MOM_fixed_initialization, only : MOM_initialize_fixed
  use MOM_get_input,       only : Get_MOM_Input, directories
  use MOM_grid,            only : ocean_grid_type, MOM_grid_init, MOM_grid_end
  use MOM_hor_index,       only : hor_index_type, hor_index_init
  use MOM_io,              only : MOM_io_init, file_exists, open_ASCII_file, close_file
  use MOM_io,              only : check_nml_error, io_infra_init, io_infra_end
  use MOM_io,              only : APPEND_FILE, READONLY_FILE, SINGLE_FILE
  use MOM_open_boundary,   only : ocean_OBC_type
  use MOM_restart,         only : save_restart
  use MOM_string_functions,only : uppercase
  use MOM_time_manager,    only : time_type, set_date, get_date
  use MOM_time_manager,    only : real_to_time, time_type_to_real
  use MOM_time_manager,    only : operator(+), operator(-), operator(*), operator(/)
  use MOM_time_manager,    only : operator(>), operator(<), operator(>=)
  use MOM_time_manager,    only : increment_date, set_calendar_type, month_name
  use MOM_time_manager,    only : JULIAN, GREGORIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
  use MOM_transcribe_grid, only : copy_dyngrid_to_MOM_grid, copy_MOM_grid_to_dyngrid
  use MOM_unit_scaling,    only : unit_scale_type, unit_scaling_init
  use MOM_verticalGrid,    only : verticalGrid_type, verticalGridInit, verticalGridEnd
  use MOM_write_cputime,   only : write_cputime, MOM_write_cputime_init
  use MOM_write_cputime,   only : write_cputime_start_clock, write_cputime_CS

  use MOM_ice_shelf, only : initialize_ice_shelf, ice_shelf_end, ice_shelf_CS
  use MOM_ice_shelf, only : ice_shelf_save_restart, solo_step_ice_shelf

  implicit none

#include <MOM_memory.h>

  logical :: use_ice_shelf = .false. ! If .true., use the ice shelf model for
                                  ! part of the domain.

  ! This is .true. if incremental restart files may be saved.
  logical :: permit_incr_restart = .true.

  integer :: ns     ! Running number of external timesteps.
  integer :: ns_ice ! Running number of internal timesteps in solo_step_ice_shelf.

  ! nmax is the number of iterations after which to stop so that the simulation does not exceed its
  ! CPU time limit.  nmax is determined by evaluating the CPU time used between successive calls to
  ! write_cputime.  Initially it is set to be very large.
  integer :: nmax=2000000000

  ! A structure containing several relevant directory paths.
  type(directories) :: dirs

  ! A suite of time types for use by the solo ice model.
  type(time_type), target :: Time       ! A copy of the model's time.
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
  type(time_type) :: Time_step_shelf    ! A time_type version of time_step.
  type(time_type) :: time_chg           ! An amount of time to adjust the segment_start_time
                                        ! and elapsed time to avoid roundoff problems.

  real    :: elapsed_time = 0.0   ! Elapsed time in this run  [s].

  logical :: elapsed_time_master  ! If true, elapsed time is used to set the
                                  ! model's master clock (Time).  This is needed
                                  ! if Time_step_shelf is not an exact
                                  ! representation of time_step.
  real :: time_step               ! The time step [s]

  ! A pointer to a structure containing metrics and related information.
  type(ocean_grid_type), pointer :: ocn_grid

  type(dyn_horgrid_type), pointer :: dG => NULL()   ! A dynamic version of the horizontal grid
  type(hor_index_type),   pointer :: HI => NULL()   ! A hor_index_type for array extents
  type(verticalGrid_type), pointer :: GV => NULL()  ! Pointer to the ocean vertical grid structure

  !> Pointer to the MOM open boundary condition type
  type(ocean_OBC_type),          pointer :: OBC => NULL()

  ! A pointer to a structure containing dimensional unit scaling factors.
  type(unit_scale_type), pointer :: US

  type(diag_ctrl), pointer :: &
    diag => NULL()              ! A pointer to the diagnostic regulatory structure

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
  logical :: symmetric

  logical :: unit_in_use
  integer :: initClock, mainClock, termClock

  type(write_cputime_CS),    pointer :: write_CPU_CSp => NULL()
  type(ice_shelf_CS),        pointer :: ice_shelf_CSp => NULL()
  !-----------------------------------------------------------------------

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mod_name = "SHELF_main (ice_shelf_driver)" ! This module's name.

  namelist /ice_solo_nml/ date_init, calendar, months, days, hours, minutes, seconds

  !=====================================================================

  call write_cputime_start_clock(write_CPU_CSp)

  call MOM_infra_init() ; call io_infra_init()

  ! These clocks are on the global pelist.
  initClock = cpu_clock_id( 'Initialization' )
  mainClock = cpu_clock_id( 'Main loop' )
  termClock = cpu_clock_id( 'Termination' )
  call cpu_clock_begin(initClock)

  call MOM_mesg('======== Model being driven by ice_shelf_driver ========', 2)
  call callTree_waypoint("Program Shelf_main, ice_shelf_driver.F90")

  if (file_exists('input.nml')) then
    ! Provide for namelist specification of the run length and calendar data.
    call open_ASCII_file(unit, 'input.nml', action=READONLY_FILE)
    read(unit, ice_solo_nml, iostat=io_status)
    call close_file(unit)
    ierr = check_nml_error(io_status,'ice_solo_nml')
    if (years+months+days+hours+minutes+seconds > 0) then
      if (is_root_pe()) write(*,ice_solo_nml)
    endif
  endif

  ! Read ocean_solo restart, which can override settings from the namelist.
  if (file_exists(trim(dirs%restart_input_dir)//'ice_solo.res')) then
    call open_ASCII_file(unit, trim(dirs%restart_input_dir)//'ice_solo.res', action=READONLY_FILE)
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
      call MOM_error(FATAL,'Shelf_driver: Invalid namelist value '//trim(calendar)//' for calendar')
    else
      call MOM_error(FATAL,'Shelf_driver: No namelist value for calendar')
    endif
  endif
  call set_calendar_type(calendar_type)


  if (sum(date_init) > 0) then
    Start_time = set_date(date_init(1),date_init(2), date_init(3), &
         date_init(4),date_init(5),date_init(6))
  else
    Start_time = real_to_time(0.0)
  endif

  call Get_MOM_Input(param_file, dirs)
  ! Determining the internal unit scaling factors for this run.
  call unit_scaling_init(param_file, US)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod_name, version, "")

  call get_param(param_file, mod_name, "ICE_SHELF", use_ice_shelf, &
                 "If true, call the code to apply an ice shelf model over "//&
                 "some of the domain.", default=.false.)

  if (.not.use_ice_shelf) call MOM_error(FATAL, "Shelf_driver: Run stops unless ICE_SHELF is true.")

  call get_param(param_file, mod_name, "ICE_VELOCITY_TIMESTEP", time_step, &
                 "The time step for changing forcing, coupling with other "//&
                 "components, or potentially writing certain diagnostics.", &
                 units="s", fail_if_missing=.true.)

  if (sum(date) >= 0) then
    ! In this case, the segment starts at a time fixed by ocean_solo.res
    segment_start_time = set_date(date(1),date(2),date(3),date(4),date(5),date(6))
    Time = segment_start_time
  else
    ! In this case, the segment starts at Start_time.
    Time = Start_time
  endif

  ! This is the start of the code that is the counterpart of MOM_initialization.
  call callTree_waypoint("Start of ice shelf initialization.")

  call MOM_debugging_init(param_file)
  call diag_mediator_infrastructure_init()
  call MOM_io_init(param_file)

  ! Set up the ocean model domain and grid; the ice model grid is set in initialize_ice_shelf,
  ! but the grids have strong commonalities in this configuration, and the ocean grid is required
  ! to set up the diag mediator control structure.
  call MOM_domains_init(ocn_grid%domain, param_file)
  call hor_index_init(ocn_grid%Domain, HI, param_file)
  call create_dyn_horgrid(dG, HI)
  call clone_MOM_domain(ocn_grid%Domain, dG%Domain)

  ! Initialize the ocean grid and topography.
  call MOM_initialize_fixed(dG, US, OBC, param_file, .true., dirs%output_directory)
  call MOM_grid_init(ocn_grid, param_file, US, HI)
  call copy_dyngrid_to_MOM_grid(dG, ocn_grid, US)
  call destroy_dyn_horgrid(dG)

  ! Initialize the diag mediator.  The ocean's vertical grid is not really used here, but at
  ! present the interface to diag_mediator_init assumes the presence of ocean-specific information.
  call verticalGridInit(param_file, GV, US)
  call diag_mediator_init(ocn_grid, GV, US, GV%ke, param_file, diag, doc_file_dir=dirs%output_directory)

  call callTree_waypoint("returned from diag_mediator_init()")

  call initialize_ice_shelf(param_file, ocn_grid, Time, ice_shelf_CSp, diag)

  ! This is the end of the code that is the counterpart of MOM_initialization.
  call callTree_waypoint("End of ice shelf initialization.")

  Master_Time = Time
!   grid => ice_shelf_CSp%grid

  segment_start_time = Time
  elapsed_time = 0.0

  Time_step_shelf = real_to_time(time_step)
  elapsed_time_master = (abs(time_step - time_type_to_real(Time_step_shelf)) > 1.0e-12*time_step)
  if (elapsed_time_master) &
    call MOM_mesg("Using real elapsed time for the master clock.", 2)

  ! Determine the segment end time, either from the namelist file or parsed input file.
  call get_param(param_file, mod_name, "TIMEUNIT", Time_unit, &
                 "The time unit for DAYMAX and RESTINT.", &
                 units="s", default=86400.0)
  if (years+months+days+hours+minutes+seconds > 0) then
    Time_end = increment_date(Time, years, months, days, hours, minutes, seconds)
    call MOM_mesg('Segment run length determined from ice_solo_nml.', 2)
    call get_param(param_file, mod_name, "DAYMAX", daymax, &
                 "The final time of the whole simulation, in units of "//&
                 "TIMEUNIT seconds.  This also sets the potential end "//&
                 "time of the present run segment if the end time is "//&
                 "not set (as it was here) via ocean_solo_nml in input.nml.", &
                 timeunit=Time_unit, default=Time_end)
  else
    call get_param(param_file, mod_name, "DAYMAX", daymax, &
                 "The final time of the whole simulation, in units of "//&
                 "TIMEUNIT seconds.  This also sets the potential end "//&
                 "time of the present run segment if the end time is "//&
                 "not set via ocean_solo_nml in input.nml.", &
                 timeunit=Time_unit, fail_if_missing=.true.)
    Time_end = daymax
  endif

  if (Time >= Time_end) call MOM_error(FATAL, &
    "Shelf_driver: The run has been started at or after the end time of the run.")

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
        (1 + ((Time + Time_step_shelf) - Start_time) / restint)
  else
    ! Set the time so late that there is no intermediate restart.
    restart_time = Time_end + Time_step_shelf
    permit_incr_restart = .false.
  endif

  call cpu_clock_end(initClock) !end initialization

  call cpu_clock_begin(mainClock) !begin main loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MAIN LOOP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ns = 1 ; ns_ice = 1
  do while ((ns < nmax) .and. (Time < Time_end))
    call callTree_enter("Main loop, Shelf_driver.F90", ns)

    ! This call steps the model over a time time_step.
    Time1 = Master_Time ; Time = Master_Time
    call solo_step_ice_shelf(ice_shelf_CSp, Time_step_shelf, ns_ice, Time)

!   Time = Time + Time_step_shelf
!   This is here to enable fractional-second time steps.
    elapsed_time = elapsed_time + time_step
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
      Master_Time = Master_Time + Time_step_shelf
    endif
    Time = Master_Time

    if (cpu_steps > 0) then ; if (MOD(ns, cpu_steps) == 0) then
      call write_cputime(Time, ns, write_CPU_CSp, nmax)
    endif ; endif

!  See if it is time to write out a restart file - timestamped or not.
    if ((permit_incr_restart) .and. (Time + (Time_step_shelf/2) > restart_time)) then
      if (BTEST(Restart_control,1)) then
        call ice_shelf_save_restart(ice_shelf_CSp, Time, dirs%restart_output_dir, .true.)
      endif
      if (BTEST(Restart_control,0)) then
        call ice_shelf_save_restart(ice_shelf_CSp, Time, dirs%restart_output_dir)
      endif
      restart_time = restart_time + restint
    endif

    ns = ns + 1
    call callTree_leave("Main loop")
  enddo

  call cpu_clock_end(mainClock)
  call cpu_clock_begin(termClock)
  if (Restart_control>=0) then
    call ice_shelf_save_restart(ice_shelf_CSp, Time, &
                                dirs%restart_output_dir)

    ! Write ice shelf solo restart file.
    if (is_root_pe())then
      call open_ASCII_file(unit, trim(dirs%restart_output_dir)//'shelf.res')
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

  call callTree_waypoint("End Shelf_main")
  call diag_mediator_end(Time, diag, end_diag_manager=.true.)
  if (cpu_steps > 0) call write_cputime(Time, ns-1, write_CPU_CSp, call_end=.true.)
  call cpu_clock_end(termClock)

  call io_infra_end ; call MOM_infra_end

  call ice_shelf_end(ice_shelf_CSp)

end program Shelf_main
