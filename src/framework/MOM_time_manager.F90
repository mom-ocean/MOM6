module MOM_time_manager
!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*    This file is a part of MOM.  See MOM.F90 for licensing.          *
!*                                                                     *
!*  By R. Hallberg November 2005                                       *
!*                                                                     *
!*    This module wraps the FMS time manager module.                   *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use time_manager_mod, only : time_type, get_time, set_time
use time_manager_mod, only : time_type_to_real, real_to_time_type
use time_manager_mod, only : operator(+), operator(-), operator(*), operator(/)
use time_manager_mod, only : operator(>), operator(<), operator(>=), operator(<=)
use time_manager_mod, only : operator(==), operator(/=), operator(//)
use time_manager_mod, only : set_ticks_per_second , get_ticks_per_second
use time_manager_mod, only : get_date, set_date, increment_date
use time_manager_mod, only : days_in_month, month_name
use time_manager_mod, only : set_calendar_type, get_calendar_type
use time_manager_mod, only : JULIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR

implicit none ; private

public :: time_type, get_time, set_time, time_type_to_real, real_to_time_type
public :: set_ticks_per_second , get_ticks_per_second
public :: operator(+), operator(-), operator(*), operator(/)
public :: operator(>), operator(<), operator(>=), operator(<=)
public :: operator(==), operator(/=), operator(//)
public :: get_date, set_date, increment_date, month_name, days_in_month
public :: JULIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
public :: set_calendar_type, get_calendar_type

end module MOM_time_manager
