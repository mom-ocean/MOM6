!> A module to monitor the overall CPU time used by MOM6 and project when to stop the model
module MOM_write_cputime

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms, only : sum_across_PEs, num_pes
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, is_root_pe
use MOM_io, only : open_file, close_file, APPEND_FILE, ASCII_FILE, WRITEONLY_FILE
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_time_manager, only : time_type, get_time, operator(>)

implicit none ; private

public write_cputime, MOM_write_cputime_init, MOM_write_cputime_end, write_cputime_start_clock

!-----------------------------------------------------------------------

integer :: CLOCKS_PER_SEC = 1000 !< The number of clock cycles per second, used by the system clock
integer :: MAX_TICKS      = 1000 !< The number of ticks per second, used by the system clock

!> A control structure that regulates the writing of CPU time
type, public :: write_cputime_CS ; private
  real :: maxcpu                !<   The maximum amount of cpu time per processor
                                !! for which MOM should run before saving a restart
                                !! file and quiting with a return value that
                                !! indicates that further execution is required to
                                !! complete the simulation, in wall-clock seconds.
  type(time_type) :: Start_time !< The start time of the simulation.
                                !! Start_time is set in MOM_initialization.F90
  real :: startup_cputime       !< The CPU time used in the startup phase of the model.
  real :: prev_cputime = 0.0    !< The last measured CPU time.
  real :: dn_dcpu_min = -1.0    !< The minimum derivative of timestep with CPU time.
  real :: cputime2 = 0.0        !< The accumulated cpu time.
  integer :: previous_calls = 0 !< The number of times write_CPUtime has been called.
  integer :: prev_n = 0         !< The value of n from the last call.
  integer :: fileCPU_ascii= -1  !< The unit number of the CPU time file.
  character(len=200) :: CPUfile !< The name of the CPU time file.
end type write_cputime_CS

contains

!> Evaluate the CPU time returned by SYSTEM_CLOCK at the start of a run
subroutine write_cputime_start_clock(CS)
  type(write_cputime_CS), pointer :: CS !< The control structure set up by a previous
                                        !! call to MOM_write_cputime_init.
  integer :: new_cputime   ! The CPU time returned by SYSTEM_CLOCK
  if (.not.associated(CS)) allocate(CS)

  call SYSTEM_CLOCK(new_cputime, CLOCKS_PER_SEC, MAX_TICKS)
  CS%prev_cputime = new_cputime
end subroutine write_cputime_start_clock

!> Initialize the MOM_write_cputime module.
subroutine MOM_write_cputime_init(param_file, directory, Input_start_time, CS)
  type(param_file_type),  intent(in) :: param_file !< A structure to parse for run-time parameters
  character(len=*),       intent(in) :: directory  !< The directory where the CPU time file goes.
  type(time_type),        intent(in) :: Input_start_time !< The start model time of the simulation.
  type(write_cputime_CS), pointer    :: CS         !< A pointer that may be set to point to the
                                                   !! control structure for this module.

  ! Local variables
  integer :: new_cputime   ! The CPU time returned by SYSTEM_CLOCK
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = 'MOM_write_cputime'  ! This module's name.
  logical :: all_default   ! If true, all parameters are using their default values.

  if (.not.associated(CS)) then
    allocate(CS)
    call SYSTEM_CLOCK(new_cputime, CLOCKS_PER_SEC, MAX_TICKS)
    CS%prev_cputime = new_cputime
  endif

  ! Read all relevant parameters and write them to the model log.

  ! Determine whether all paramters are set to their default values.
  call get_param(param_file, mdl, "MAXCPU", CS%maxcpu, default=-1.0, do_not_log=.true.)
  call get_param(param_file, mdl, "CPU_TIME_FILE", CS%CPUfile, default="CPU_stats", do_not_log=.true.)
  all_default = (CS%maxcpu == -1.0) .and. (trim(CS%CPUfile) == trim("CPU_stats"))

  call log_version(param_file, mdl, version, "", all_default=all_default)
  call get_param(param_file, mdl, "MAXCPU", CS%maxcpu, &
                 "The maximum amount of cpu time per processor for which "//&
                 "MOM should run before saving a restart file and "//&
                 "quitting with a return value that indicates that a "//&
                 "further run is required to complete the simulation. "//&
                 "If automatic restarts are not desired, use a negative "//&
                 "value for MAXCPU.  MAXCPU has units of wall-clock "//&
                 "seconds, so the actual CPU time used is larger by a "//&
                 "factor of the number of processors used.", &
                 units="wall-clock seconds", default=-1.0)
  call get_param(param_file, mdl, "CPU_TIME_FILE", CS%CPUfile, &
                 "The file into which CPU time is written.",default="CPU_stats")
  CS%CPUfile = trim(directory)//trim(CS%CPUfile)
  call log_param(param_file, mdl, "directory/CPU_TIME_FILE", CS%CPUfile)
#ifdef STATSLABEL
  CS%CPUfile = trim(CS%CPUfile)//"."//trim(adjustl(STATSLABEL))
#endif

  CS%Start_time = Input_start_time

end subroutine MOM_write_cputime_init

!> Close the MOM_write_cputime module.
subroutine MOM_write_cputime_end(CS)
  type(write_cputime_CS), pointer    :: CS    !< The control structure set up by a previous
                                              !! call to MOM_write_cputime_init.

  if (.not.associated(CS)) return

  ! Flush and close the output files.
  if (is_root_pe() .and. CS%fileCPU_ascii > 0) then
    flush(CS%fileCPU_ascii)
    call close_file(CS%fileCPU_ascii)
  endif

  deallocate(CS)

end subroutine MOM_write_cputime_end

!> This subroutine assesses how much CPU time the model has taken and determines how long the model
!! should be run before it saves a restart file and stops itself.  Optionally this may also be used
!! to trigger this module's end routine.
subroutine write_cputime(day, n, CS, nmax, call_end)
  type(time_type),        intent(inout) :: day  !< The current model time.
  integer,                intent(in)    :: n    !< The time step number of the current execution.
  type(write_cputime_CS), pointer       :: CS   !< The control structure set up by a previous
                                                !! call to MOM_write_cputime_init.
  integer,      optional, intent(inout) :: nmax !< The number of iterations after which to stop so
                                                !! that the simulation will not run out of CPU time.
  logical,      optional, intent(in)    :: call_end !< If true, also call MOM_write_cputime_end.

  ! Local variables
  real    :: d_cputime     ! The change in CPU time since the last call
                           ! this subroutine.
  integer :: new_cputime   ! The CPU time returned by SYSTEM_CLOCK
  real    :: reday         ! A real version of day.
  character(len=256) :: mesg  ! The text of an error message
  integer :: start_of_day, num_days

  if (.not.associated(CS)) call MOM_error(FATAL, &
         "write_energy: Module must be initialized before it is used.")

  call SYSTEM_CLOCK(new_cputime, CLOCKS_PER_SEC, MAX_TICKS)
!   The following lines extract useful information even if the clock has rolled
! over, assuming a 32-bit SYSTEM_CLOCK.  With more bits, rollover is essentially
! impossible. Negative fluctuations of less than 10 seconds are not interpreted
! as the clock rolling over.  This should be unnecessary but is sometimes needed
! on the GFDL SGI/O3k.
  if (new_cputime < CS%prev_cputime-(10.0*CLOCKS_PER_SEC)) then
    d_cputime = new_cputime - CS%prev_cputime + MAX_TICKS
  else
    d_cputime = new_cputime - CS%prev_cputime
  endif

  call sum_across_PEs(d_cputime)
  if (CS%previous_calls == 0) CS%startup_cputime = d_cputime

  CS%cputime2 = CS%cputime2 + d_cputime

  if ((CS%previous_calls >= 1) .and. (CS%maxcpu > 0.0)) then
    ! Determine the slowest rate at which time steps are executed.
    if ((n > CS%prev_n) .and. (d_cputime > 0.0) .and. &
        ((CS%dn_dcpu_min*d_cputime < (n - CS%prev_n)) .or. &
         (CS%dn_dcpu_min < 0.0))) &
      CS%dn_dcpu_min = (n - CS%prev_n) / d_cputime
    if (present(nmax) .and. (CS%dn_dcpu_min >= 0.0)) then
      ! Have the model stop itself after 95% of the CPU time has been used.
      nmax = n + INT( CS%dn_dcpu_min * &
          (0.95*CS%maxcpu * REAL(num_pes())*CLOCKS_PER_SEC - &
           (CS%startup_cputime + CS%cputime2)) )
!     write(mesg,*) "Resetting nmax to ",nmax," at day",reday
!     call MOM_mesg(mesg)
    endif
  endif
  CS%prev_cputime = new_cputime ; CS%prev_n = n

  call get_time(day, start_of_day, num_days)
  reday = REAL(num_days)+ (REAL(start_of_day)/86400.0)

  !  Reopen or create a text output file.
  if ((CS%previous_calls == 0) .and. (is_root_pe())) then
    if (day > CS%Start_time) then
      call open_file(CS%fileCPU_ascii, trim(CS%CPUfile), &
                     action=APPEND_FILE, form=ASCII_FILE, nohdrs=.true.)
    else
      call open_file(CS%fileCPU_ascii, trim(CS%CPUfile), &
                     action=WRITEONLY_FILE, form=ASCII_FILE, nohdrs=.true.)
    endif
  endif

  if (is_root_pe()) then
    if (CS%previous_calls == 0) then
      write(CS%fileCPU_ascii, &
        '("Startup CPU time: ", F12.3, " sec summed across", I5, " PEs.")') &
                            (CS%startup_cputime / CLOCKS_PER_SEC), num_pes()
      write(CS%fileCPU_ascii,*)"        Day, Step number,     CPU time, CPU time change"
    endif
    write(CS%fileCPU_ascii,'(F12.3,", "I11,", ", F12.3,", ", F12.3)') &
           reday, n, (CS%cputime2 / real(CLOCKS_PER_SEC)), &
           d_cputime / real(CLOCKS_PER_SEC)

    flush(CS%fileCPU_ascii)
  endif
  CS%previous_calls = CS%previous_calls + 1

  if (present(call_end)) then
    if (call_end) call MOM_write_cputime_end(CS)
  endif

end subroutine write_cputime

!> \namespace mom_write_cputime
!!
!!  By Robert Hallberg, May 2006.
!!
!!    This file contains the subroutine (write_cputime) that writes
!!  the summed CPU time across all processors to an output file. In
!!  addition, write_cputime estimates how many more time steps can be
!!  taken before 95% of the available CPU time is used, so that the
!!  model can be checkpointed at that time.

end module MOM_write_cputime
