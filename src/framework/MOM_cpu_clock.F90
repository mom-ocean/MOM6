!> Wraps the MPP cpu clock functions
module MOM_cpu_clock

! This file is part of MOM6. See LICENSE.md for the license.

use fms_mod, only : clock_flag_default
use mpp_mod, only : cpu_clock_begin => mpp_clock_begin
use mpp_mod, only : cpu_clock_end => mpp_clock_end, mpp_clock_id
use mpp_mod, only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE_DRIVER
use mpp_mod, only : CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA
use mpp_mod, only : CLOCK_SYNC => MPP_CLOCK_SYNC

implicit none ; private

public :: cpu_clock_id, cpu_clock_begin, cpu_clock_end
public :: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE_DRIVER, CLOCK_MODULE
public :: CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA, CLOCK_SYNC

contains

!> cpu_clock_id returns the integer handle for a named CPU clock.
function cpu_clock_id( name, synchro_flag, grain )
  character(len=*),  intent(in) :: name  !< The unique name of the CPU clock
  integer, intent(in), optional :: synchro_flag !< An integer flag that controls whether the PEs
                                       !! are synchronized before the cpu clocks start counting.
                                       !! Synchronization occurs before the start of a clock if this
                                       !! is odd, while additional (expensive) statistics can set
                                       !! for other values. If absent, the default is taken from the
                                       !! settings for FMS.
  integer, intent(in), optional :: grain !< The timing granularity for this clock, usually set to
                                       !! the values of CLOCK_COMPONENT, CLOCK_ROUTINE, CLOCK_LOOP, etc.
  integer                       :: cpu_clock_id !< The integer CPU clock handle.

  if (present(synchro_flag)) then
    cpu_clock_id = mpp_clock_id(name, flags=synchro_flag, grain=grain)
  else
    cpu_clock_id = mpp_clock_id(name, flags=clock_flag_default, grain=grain)
  endif

end function cpu_clock_id

end module MOM_cpu_clock
