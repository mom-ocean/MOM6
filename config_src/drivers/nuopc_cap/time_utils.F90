!> Set of time utilities for converting between FMS and ESMF time type.
module time_utils_mod

! FMS
use fms_mod,            only: uppercase
use mpp_mod,            only: mpp_error, FATAL
use time_manager_mod,   only: time_type, set_time, set_date, get_date
use time_manager_mod,   only: GREGORIAN, JULIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
use time_manager_mod,   only: fms_get_calendar_type => get_calendar_type
! ESMF
use ESMF,               only: ESMF_CALKIND_FLAG, ESMF_CALKIND_GREGORIAN
use ESMF,               only: ESMF_CALKIND_JULIAN, ESMF_CALKIND_NOLEAP
use ESMF,               only: ESMF_CALKIND_360DAY, ESMF_CALKIND_NOCALENDAR
use ESMF,               only: ESMF_Time, ESMF_TimeGet, ESMF_LogFoundError
use ESMF,               only: ESMF_LOGERR_PASSTHRU,ESMF_TimeInterval
use ESMF,               only: ESMF_TimeIntervalGet, ESMF_TimeSet, ESMF_SUCCESS
use MOM_cap_methods,    only: ChkErr

implicit none; private

!> Converts calendar from FMS to ESMF format
interface fms2esmf_cal
  module procedure fms2esmf_cal_c
  module procedure fms2esmf_cal_i
end interface fms2esmf_cal

!> Converts time from FMS to ESMF format
interface esmf2fms_time
  module procedure esmf2fms_time_t
  module procedure esmf2fms_timestep
end interface esmf2fms_time

public fms2esmf_cal
public esmf2fms_time
public fms2esmf_time
public string_to_date

character(len=*),parameter :: u_FILE_u = &
     __FILE__

contains

!> Sets fms2esmf_cal_c to the corresponding ESMF calendar type
function fms2esmf_cal_c(calendar)
  type(ESMF_CALKIND_FLAG)            :: fms2esmf_cal_c !< ESMF calendar type
  character(len=*), intent(in)       :: calendar       !< Type of calendar

  select case( uppercase(trim(calendar)) )
  case( 'GREGORIAN' )
    fms2esmf_cal_c = ESMF_CALKIND_GREGORIAN
  case( 'JULIAN' )
    fms2esmf_cal_c = ESMF_CALKIND_JULIAN
  case( 'NOLEAP' )
    fms2esmf_cal_c = ESMF_CALKIND_NOLEAP
  case( 'THIRTY_DAY' )
    fms2esmf_cal_c = ESMF_CALKIND_360DAY
  case( 'NO_CALENDAR' )
    fms2esmf_cal_c = ESMF_CALKIND_NOCALENDAR
  case default
    call mpp_error(FATAL, &
    'ocean_solo: ocean_solo_nml entry calendar must be one of GREGORIAN|JULIAN|NOLEAP|THIRTY_DAY|NO_CALENDAR.' )
  end select
end function fms2esmf_cal_c

!> Sets fms2esmf_cal_i to the corresponding ESMF calendar type
function fms2esmf_cal_i(calendar)
  type(ESMF_CALKIND_FLAG)            :: fms2esmf_cal_i !< ESMF calendar structure
  integer, intent(in)                :: calendar       !< Type of calendar

  select case(calendar)
    case(THIRTY_DAY_MONTHS)
      fms2esmf_cal_i = ESMF_CALKIND_360DAY
    case(GREGORIAN)
      fms2esmf_cal_i = ESMF_CALKIND_GREGORIAN
    case(JULIAN)
      fms2esmf_cal_i = ESMF_CALKIND_JULIAN
    case(NOLEAP)
      fms2esmf_cal_i = ESMF_CALKIND_NOLEAP
    case(NO_CALENDAR)
      fms2esmf_cal_i = ESMF_CALKIND_NOCALENDAR
  end select
end function fms2esmf_cal_i

!> Converts date from ESMF format to FMS format.
function esmf2fms_time_t(time)
  type(Time_type)                    :: esmf2fms_time_t !< FMS time structure
  type(ESMF_Time), intent(in)        :: time            !< ESMF time structure

  ! Local Variables
  integer                            :: yy, mm, dd, h, m, s
  type(ESMF_CALKIND_FLAG)            :: calkind

  integer                            :: rc

  call ESMF_TimeGet(time, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, &
      calkindflag=calkind, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  esmf2fms_time_t = set_date(yy, mm, dd, h, m, s)

end function esmf2fms_time_t

!> Converts time-interval from ESMF format to FMS format.
function esmf2fms_timestep(timestep)
  type(Time_type)                    :: esmf2fms_timestep !< FMS time structure
  type(ESMF_TimeInterval), intent(in):: timestep          !< time-interval following
                                                          !! ESMF format [s]
  ! Local Variables
  integer                            :: s
  type(ESMF_CALKIND_FLAG)            :: calkind

  integer                            :: rc

  call ESMF_TimeIntervalGet(timestep, s=s, calkindflag=calkind, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  esmf2fms_timestep = set_time(s, 0)

end function esmf2fms_timestep

!> Converts date from FMS format to ESMF format.
function fms2esmf_time(time, calkind)
  type(ESMF_Time)                    :: fms2esmf_time !< ESMF time structure
  type(time_type), intent(in)        :: time          !< FMS time structure
  type(ESMF_CALKIND_FLAG), intent(in), optional :: calkind !< ESMF calendar structure

  ! Local Variables
  integer                            :: yy, mm, d, h, m, s
  type(ESMF_CALKIND_FLAG)            :: l_calkind

  integer                            :: rc

  if(present(calkind)) then
    l_calkind = calkind
  else
    l_calkind = fms2esmf_cal(fms_get_calendar_type())
  endif

  call get_date(time, yy, mm, d, h, m, s)

  call ESMF_TimeSet(fms2esmf_time, yy=yy, mm=mm, d=d, h=h, m=m, s=s, &
      calkindflag=l_calkind, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

end function fms2esmf_time

!> Converts a string (I4.4,I2.2,I2.2,".",I2.2,I2.2,I2.2) that represents
!! yr, mon, day, hr, min, sec to a FMS data format.
function string_to_date(string, rc)
  character(len=15), intent(in)           :: string        !< String representing a date
  integer, intent(out), optional          :: rc            !< ESMF error handler
  type(time_type)                         :: string_to_date!< FMS time structure

  ! Local variables
  integer                                 :: yr,mon,day,hr,min,sec

  if(present(rc)) rc = ESMF_SUCCESS

  read(string, '(I4.4,I2.2,I2.2,".",I2.2,I2.2,I2.2)') yr, mon, day, hr, min, sec
  string_to_date = set_date(yr, mon, day, hr, min, sec)

end function string_to_date

end module time_utils_mod
