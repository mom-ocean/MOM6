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

implicit none ; private

public :: time_type, get_time, set_time
public :: time_type_to_real, real_to_time_type, real_to_time
public :: set_ticks_per_second, get_ticks_per_second
public :: operator(+), operator(-), operator(*), operator(/)
public :: operator(>), operator(<), operator(>=), operator(<=)
public :: operator(==), operator(/=), operator(//)
public :: get_date, set_date, increment_date, month_name, days_in_month
public :: JULIAN, NOLEAP, THIRTY_DAY_MONTHS, GREGORIAN, NO_CALENDAR
public :: set_calendar_type, get_calendar_type

contains

!> Returns a time_type version of a real time in seconds, using an alternate implementation to the
!! FMS function real_to_time_type that is accurate over a larger range of input values.  With 32 bit
!! signed integers, this version should work over the entire valid range (2^31 days or ~5.8835
!! million years) of time_types, whereas the standard version in the FMS time_manager stops working
!! for conversions of times greater than 2^31 seconds, or ~68.1 years.
type(time_type) function real_to_time(x, err_msg)
!  type(time_type)  :: real_to_time !< The output time as a time_type
  real,                       intent(in)  :: x       !< The input time in real seconds.
  character(len=*), optional, intent(out) :: err_msg !< An optional returned error message.

  ! Local variables
  integer          :: seconds, days, ticks
  real             :: real_subsecond_remainder

  days = floor(x/86400.)
  seconds = floor(x - 86400.*days)
  real_subsecond_remainder = x - (days*86400. + seconds)
  ticks = nint(real_subsecond_remainder * get_ticks_per_second())

  real_to_time = set_time(seconds=seconds, days=days, ticks=ticks, err_msg=err_msg)
end function real_to_time

end module MOM_time_manager
