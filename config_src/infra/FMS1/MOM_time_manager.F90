! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Wraps the FMS time manager functions
module MOM_time_manager

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

! FMS re-exports
public :: time_type, get_time, set_time
public :: time_type_to_real, real_to_time_type
public :: set_ticks_per_second, get_ticks_per_second
public :: operator(+), operator(-), operator(*), operator(/)
public :: operator(>), operator(<), operator(>=), operator(<=)
public :: operator(==), operator(/=), operator(//)
public :: get_date, set_date, increment_date, month_name, days_in_month
public :: JULIAN, NOLEAP, THIRTY_DAY_MONTHS, GREGORIAN, NO_CALENDAR
public :: set_calendar_type, get_calendar_type
! Module functions
public :: real_to_time, time_minus_signed, time_to_real

contains

!> Returns a time_type version of a real time in seconds, using an alternate implementation to the
!! FMS function real_to_time_type that is accurate over a larger range of input values.  With 32 bit
!! signed integers, this version should work over the entire valid range (2^31 days or ~5.8835
!! million years) of time_types, whereas the standard version in the FMS time_manager stops working
!! for conversions of times greater than 2^31 seconds, or ~68.1 years.
type(time_type) function real_to_time(time_in, err_msg, unscale)
!  type(time_type)  :: real_to_time !< The output time as a time_type
  real,                       intent(in)  :: time_in !< The input time in [s] or [T ~> s]
  character(len=*), optional, intent(out) :: err_msg !< An optional returned error message.
  real,             optional, intent(in)  :: unscale !< A scaling factor that the input time is
                                                     !! multiplied by, often in [s T-1 ~> nondim]

  ! Local variables
  real             :: x  ! The time in real seconds [s]
  real             :: real_subsecond_remainder ! The fractional seconds from time_in [s]
  integer          :: seconds, days, ticks

  x = time_in ; if (present(unscale)) x = unscale*time_in

  days = floor(x/86400.)
  seconds = floor(x - 86400.*days)
  real_subsecond_remainder = x - (days*86400. + seconds)
  ticks = nint(real_subsecond_remainder * get_ticks_per_second())

  real_to_time = set_time(seconds=seconds, days=days, ticks=ticks, err_msg=err_msg)

end function real_to_time

!> Returns the real number of seconds encoded in the time type [s] or the rescaled
!! amount of time in other units, often [T ~> s]
real function time_to_real(time, scale)
  type(time_type), intent(in) :: time   !< The time to be converted.
  real, optional,  intent(in) :: scale  !< A scaling factor that returned the time is
                                        !! multiplied by, often in [T s-1 ~> nondim]

  time_to_real = time_type_to_real(time)
  if (present(scale)) time_to_real = scale * time_to_real

end function time_to_real

!> Returns a real number representing time_a - time_b in [s] or [T ~> s] if scale is present.
!! The FMS - operator for time types returns a new time type representing
!! a difference that is always >= 0.
!! In contrast, this function returns a negative real number if time_b > time_a,
!! and a positive real otherwise, as would be expected for subtraction.
real function time_minus_signed(time_a, time_b, scale)
  type(time_type), intent(in) :: time_a, time_b !< Two times for calculating time_a - time_b
  real, optional,  intent(in) :: scale          !< A scaling factor that returned the time is
                                                !! multiplied by, often in [T s-1 ~> nondim]

  ! Local variables
  real :: abs_diff ! The absolute value of the difference in times [s] or [T ~> s]

  ! Do FMS time subtraction, which will always be >= 0,
  ! and convert to a real number.
  abs_diff = time_to_real(time_a - time_b, scale)

  ! Add the sign back by comparing time_a and time_b
  time_minus_signed = merge(abs_diff, -abs_diff, time_a >= time_b)

end function time_minus_signed

end module MOM_time_manager
