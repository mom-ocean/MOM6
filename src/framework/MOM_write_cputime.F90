module MOM_write_cputime
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, May 2006.                                      *
!*                                                                     *
!*    This file contains the subroutine (write_cputime) that writes    *
!*  the summed CPU time across all processors to an output file. In    *
!*  addition, write_cputime estimates how many more time steps can be  *
!*  taken before 95% of the available CPU time is used, so that the    *
!*  model can be checkpointed at that time.                            *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_coms, only : sum_across_PEs, pe_here, num_pes
use MOM_error_handler, only : MOM_error, FATAL, is_root_pe
use MOM_io, only : open_file, APPEND_FILE, ASCII_FILE, WRITEONLY_FILE
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_time_manager, only : time_type, get_time, operator(>)

implicit none ; private

public write_cputime, MOM_write_cputime_init, write_cputime_start_clock

!-----------------------------------------------------------------------

integer :: CLOCKS_PER_SEC = 1000
integer :: MAX_TICKS      = 1000

type, public :: write_cputime_CS ; private
  real :: maxcpu                !   The maximum amount of cpu time per processor
                                ! for which MOM should run before saving a restart
                                ! file and quiting with a return value that
                                ! indicates that further execution is required to
                                ! complete the simulation, in wall-clock seconds.
  type(time_type) :: Start_time ! The start time of the simulation.
                                ! Start_time is set in MOM_initialization.F90
  real :: startup_cputime       ! The CPU time used in the startup phase of the model.
  real :: prev_cputime = 0.0    ! The last measured CPU time.
  real :: dn_dcpu_min = -1.0    ! The minimum derivative of timestep with CPU time.
  real :: cputime2 = 0.0        ! The accumulated cpu time.
  integer :: previous_calls = 0 ! The number of times write_CPUtime has been called.
  integer :: prev_n = 0         ! The value of n from the last call.
  integer :: fileCPU_ascii      ! The unit number of the CPU time file.
  character(len=200) :: CPUfile ! The name of the CPU time file.
end type write_cputime_CS

contains

subroutine write_cputime_start_clock(CS)
  type(write_cputime_CS), pointer :: CS
! Argument:  CS - A pointer that is set to point to the control structure
!                 for this module
  integer :: new_cputime   ! The CPU time returned by SYSTEM_CLOCK
  if (.not.associated(CS)) allocate(CS)

  call SYSTEM_CLOCK(new_cputime, CLOCKS_PER_SEC, MAX_TICKS)
  CS%prev_cputime = new_cputime
end subroutine write_cputime_start_clock

subroutine MOM_write_cputime_init(param_file, directory, Input_start_time, CS)
  type(param_file_type),  intent(in) :: param_file
  character(len=*),       intent(in) :: directory
  type(time_type),        intent(in) :: Input_start_time
  type(write_cputime_CS), pointer    :: CS
! Arguments: param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      directory - The directory where the energy file goes.
!  (in)      Input_start_time - The start time of the simulation.
!  (in/out)  CS - A pointer that may be set to point to the control structure
!                 for this module.
  integer :: new_cputime   ! The CPU time returned by SYSTEM_CLOCK
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = 'MOM_write_cputime'  ! This module's name.

  if (.not.associated(CS)) then
    allocate(CS)
    call SYSTEM_CLOCK(new_cputime, CLOCKS_PER_SEC, MAX_TICKS)
    CS%prev_cputime = new_cputime
  endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "MAXCPU", CS%maxcpu, &
                 "The maximum amount of cpu time per processor for which \n"//&
                 "MOM should run before saving a restart file and \n"//&
                 "quitting with a return value that indicates that a \n"//&
                 "further run is required to complete the simulation. \n"//&
                 "If automatic restarts are not desired, use a negative \n"//&
                 "value for MAXCPU.  MAXCPU has units of wall-clock \n"//&
                 "seconds, so the actual CPU time used is larger by a \n"//&
                 "factor of the number of processors used.", &
                 units="wall-clock seconds", default=-1.0)
  call get_param(param_file, mod, "CPU_TIME_FILE", CS%CPUfile, &
                 "The file into which CPU time is written.",default="CPU_stats")
  CS%CPUfile = trim(directory)//trim(CS%CPUfile)
  call log_param(param_file, mod, "directory/CPU_TIME_FILE", CS%CPUfile)
#ifdef STATSLABEL
  CS%CPUfile = trim(CS%CPUfile)//"."//trim(adjustl(STATSLABEL))
#endif

  CS%Start_time = Input_start_time

end subroutine MOM_write_cputime_init

subroutine write_cputime(day, n, nmax, CS)
  type(time_type),                     intent(inout) :: day
  integer,                             intent(in)    :: n
  integer,                             intent(inout) :: nmax
  type(write_cputime_CS),              pointer       :: CS
!  This subroutine assesses how much CPU time the model has
! taken and determines how long the model should be run before it
! saves a restart file and stops itself.

! Arguments: day - The current model time.
!  (in)      n - The time step number of the current execution.
!  (out)     nmax - The number of iterations after which to stop so
!                   that the simulation will not run out of CPU time.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 MOM_write_cputime_init.
  real    :: d_cputime     ! The change in CPU time since the last call
                           ! this subroutine.
  integer :: new_cputime   ! The CPU time returned by SYSTEM_CLOCK
  real    :: reday         ! A real version of day.
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
    if (CS%dn_dcpu_min >= 0.0) then
      ! Have the model stop itself after 95% of the CPU time has been used.
      nmax = n + INT( CS%dn_dcpu_min * &
          (0.95*CS%maxcpu * REAL(num_pes())*CLOCKS_PER_SEC - &
           (CS%startup_cputime + CS%cputime2)) )
!     if (is_root_pe() ) then
!       write(*,*) "Resetting nmax to ",nmax," at day",reday
!     endif
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
  endif
  CS%previous_calls = CS%previous_calls + 1

end subroutine write_cputime

end module MOM_write_cputime
