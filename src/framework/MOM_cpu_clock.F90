!> Provides cpu clock functions
module MOM_cpu_clock

! This file is part of MOM6. See LICENSE.md for the license.

! These interfaces and constants from MPP/FMS will not be directly exposed outside of this module
use MOM_cpu_clock_infra, only : cpu_clock_begin
use MOM_cpu_clock_infra, only : cpu_clock_end
use MOM_cpu_clock_infra, only : cpu_clock_id
use MOM_cpu_clock_infra, only : CLOCK_COMPONENT
use MOM_cpu_clock_infra, only : CLOCK_SUBCOMPONENT
use MOM_cpu_clock_infra, only : CLOCK_MODULE_DRIVER
use MOM_cpu_clock_infra, only : CLOCK_MODULE
use MOM_cpu_clock_infra, only : CLOCK_ROUTINE
use MOM_cpu_clock_infra, only : CLOCK_LOOP
use MOM_cpu_clock_infra, only : CLOCK_INFRA

implicit none ; private

!> Public functions:
!> mom_cpu_clock_infra::cpu_clock_id, mom_cpu_clock_infra::cpu_clock_begin, mom_cpu_clock_infra::cpu_clock_end
public :: cpu_clock_id, cpu_clock_begin, cpu_clock_end

!> Public constants:
!> mom_cpu_clock_infra::clock_component, mom_cpu_clock_infra::clock_subcomponent
!> mom_cpu_clock_infra::clock_module_driver, mom_cpu_clock_infra::clock_module_driver
public :: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE_DRIVER, CLOCK_MODULE
!> mom_cpu_clock_infra::clock_routine, mom_cpu_clock_infra::clock_loop
!> mom_cpu_clock_infra::clock_infra
public :: CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA

end module MOM_cpu_clock

!> \namespace mom_cpu_clock
!!
!! APIs are defined and implemented in mom_cpu_clock_infra.
