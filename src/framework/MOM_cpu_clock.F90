!> Wraps the MPP cpu clock functions
module MOM_cpu_clock

! This file is part of MOM6. See LICENSE.md for the license.

use mpp_mod, only : cpu_clock_begin => mpp_clock_begin
use mpp_mod, only : cpu_clock_end => mpp_clock_end, cpu_clock_id => mpp_clock_id
use mpp_mod, only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE_DRIVER
use mpp_mod, only : CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA
use mpp_mod, only : CLOCK_SYNC => MPP_CLOCK_SYNC

implicit none ; private

public :: cpu_clock_id, cpu_clock_begin, cpu_clock_end
public :: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE_DRIVER, CLOCK_MODULE
public :: CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA, CLOCK_SYNC

end module MOM_cpu_clock
