!> Wraps the MPP cpu clock functions
!!
!! The functions and constants should be accessed via mom_cpu_clock
module MOM_cpu_clock_infra

! This file is part of MOM6. See LICENSE.md for the license.

! These interfaces and constants from MPP/FMS will not be directly exposed outside of this module
use fms_mod, only : clock_flag_default
use mpp_mod, only : mpp_clock_begin
use mpp_mod, only : mpp_clock_end, mpp_clock_id
use mpp_mod, only : MPP_CLOCK_COMPONENT => CLOCK_COMPONENT
use mpp_mod, only : MPP_CLOCK_SUBCOMPONENT => CLOCK_SUBCOMPONENT
use mpp_mod, only : MPP_CLOCK_MODULE_DRIVER => CLOCK_MODULE_DRIVER
use mpp_mod, only : MPP_CLOCK_MODULE => CLOCK_MODULE
use mpp_mod, only : MPP_CLOCK_ROUTINE => CLOCK_ROUTINE
use mpp_mod, only : MPP_CLOCK_LOOP => CLOCK_LOOP
use mpp_mod, only : MPP_CLOCK_INFRA => CLOCK_INFRA

implicit none ; private

! Public entities
public :: cpu_clock_id, cpu_clock_begin, cpu_clock_end
public :: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE_DRIVER, CLOCK_MODULE
public :: CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA

!> A granularity value to passed to cpu_clock_id() to indicate the clock is for a
!! component, e.g. the entire MOM6 model
integer, parameter :: CLOCK_COMPONENT = MPP_CLOCK_COMPONENT

!> A granularity value to passed to cpu_clock_id() to indicate the clock is for a
!! sub-component, e.g. dynamics or thermodynamics
integer, parameter :: CLOCK_SUBCOMPONENT = MPP_CLOCK_SUBCOMPONENT

!> A granularity value to passed to cpu_clock_id() to indicate the clock is for a
!! module driver, e.g. a routine that calls multiple other routines
integer, parameter :: CLOCK_MODULE_DRIVER = MPP_CLOCK_MODULE_DRIVER

!> A granularity value to passed to cpu_clock_id() to indicate the clock is for a
!! module, e.g. the main entry routine for a module
integer, parameter :: CLOCK_MODULE = MPP_CLOCK_MODULE

!> A granularity value to passed to cpu_clock_id() to indicate the clock is for a
!! subroutine or function
integer, parameter :: CLOCK_ROUTINE = MPP_CLOCK_ROUTINE

!> A granularity value to passed to cpu_clock_id() to indicate the clock is for a
!! section with in a routine, e.g. around a loop
integer, parameter :: CLOCK_LOOP = MPP_CLOCK_LOOP

!> A granularity value to passed to cpu_clock_id() to indicate the clock is for an
!! infrastructure operation, e.g. a halo update
integer, parameter :: CLOCK_INFRA = MPP_CLOCK_INFRA

contains

!> Turns on clock with handle "id"
subroutine cpu_clock_begin(id)
  integer, intent(in) :: id !< Handle for clock

  call mpp_clock_begin(id)

end subroutine cpu_clock_begin

!> Turns off clock with handle "id"
subroutine cpu_clock_end(id)
  integer, intent(in) :: id !< Handle for clock

  call mpp_clock_end(id)

end subroutine cpu_clock_end

!> Returns the integer handle for a named CPU clock.
integer function cpu_clock_id( name, synchro_flag, grain )
  character(len=*),  intent(in) :: name  !< The unique name of the CPU clock
  integer, optional, intent(in) :: synchro_flag !< An integer flag that controls whether the PEs
                                       !! are synchronized before the cpu clocks start counting.
                                       !! Synchronization occurs before the start of a clock if this
                                       !! is odd, while additional (expensive) statistics can set
                                       !! for other values. If absent, the default is taken from the
                                       !! settings for FMS.
  integer, optional, intent(in) :: grain !< The timing granularity for this clock, usually set to
                                       !! the values of CLOCK_COMPONENT, CLOCK_ROUTINE, CLOCK_LOOP, etc.

  if (present(synchro_flag)) then
    cpu_clock_id = mpp_clock_id(name, flags=synchro_flag, grain=grain)
  else
    cpu_clock_id = mpp_clock_id(name, flags=clock_flag_default, grain=grain)
  endif

end function cpu_clock_id

end module MOM_cpu_clock_infra
