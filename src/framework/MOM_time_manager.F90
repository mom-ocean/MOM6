!> Wraps the FMS time manager functions
module MOM_time_manager

! This file is part of MOM6. See LICENSE.md for the license.

use time_manager_mod, only : time_type, get_time, set_time
use time_manager_mod, only : time_type_to_real, real_to_time_type
use time_manager_mod, only : operator(+), operator(-), operator(*), operator(/)
use time_manager_mod, only : operator(>), operator(<), operator(>=), operator(<=)
use time_manager_mod, only : operator(==), operator(/=), operator(//)
use time_manager_mod, only : set_ticks_per_second , get_ticks_per_second
use time_manager_mod, only : get_date, set_date, increment_date
use time_manager_mod, only : days_in_month, month_name
use time_manager_mod, only : set_calendar_type, get_calendar_type
use time_manager_mod, only : JULIAN, NOLEAP, THIRTY_DAY_MONTHS, GREGORIAN
use time_manager_mod, only : NO_CALENDAR
use time_interp_external_mod, only : init_external_field, time_interp_external, time_interp_external_init
use time_interp_external_mod, only : get_external_field_size
use time_interp_external_mod, only : get_external_field_axes, get_external_field_missing

implicit none ; private

public :: time_type, get_time, set_time, time_type_to_real, real_to_time_type
public :: set_ticks_per_second , get_ticks_per_second
public :: operator(+), operator(-), operator(*), operator(/)
public :: operator(>), operator(<), operator(>=), operator(<=)
public :: operator(==), operator(/=), operator(//)
public :: get_date, set_date, increment_date, month_name, days_in_month
public :: JULIAN, NOLEAP, THIRTY_DAY_MONTHS, GREGORIAN, NO_CALENDAR
public :: set_calendar_type, get_calendar_type
public :: init_external_field
public :: time_interp_external
public :: time_interp_external_init
public :: get_external_field_size
public :: get_external_field_axes
public :: get_external_field_missing

end module MOM_time_manager
