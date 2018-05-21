module time_utils_mod

  use fms_mod,                  only: uppercase
  use mpp_mod,                  only: mpp_error, FATAL
  use time_manager_mod,         only: time_type, set_time, set_date, get_date
  use time_manager_mod,         only: GREGORIAN, JULIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
  use time_manager_mod,         only: fms_get_calendar_type => get_calendar_type
  use ESMF
  
  implicit none
  private

  !-------------------- interface blocks ---------------------
  interface fms2esmf_cal
    module procedure fms2esmf_cal_c
    module procedure fms2esmf_cal_i
  end interface fms2esmf_cal
  interface esmf2fms_time
    module procedure esmf2fms_time_t
    module procedure esmf2fms_timestep
  end interface esmf2fms_time

  public fms2esmf_cal
  public esmf2fms_time
  public fms2esmf_time
  public string_to_date

  contains

  !-------------------- module code ---------------------

  function fms2esmf_cal_c(calendar)
!    ! Return Value:
    type(ESMF_CALKIND_FLAG)            :: fms2esmf_cal_c
!    ! Arguments:
    character(len=*), intent(in)       :: calendar

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

  function fms2esmf_cal_i(calendar)
!    ! Return Value:
    type(ESMF_CALKIND_FLAG)            :: fms2esmf_cal_i
!    ! Arguments:
    integer, intent(in)                :: calendar

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

  function esmf2fms_time_t(time)
    ! Return Value
    type(Time_type)                    :: esmf2fms_time_t
    ! Input Arguments
    type(ESMF_Time), intent(in)        :: time
    ! Local Variables
    integer                            :: yy, mm, dd, h, m, s
    type(ESMF_CALKIND_FLAG)            :: calkind

    integer                            :: rc

    call ESMF_TimeGet(time, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, &
        calkindflag=calkind, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    esmf2fms_time_t = Set_date(yy, mm, dd, h, m, s)

  end function esmf2fms_time_t

  function esmf2fms_timestep(timestep)
    ! Return Value
    type(Time_type)                    :: esmf2fms_timestep
    ! Input Arguments
    type(ESMF_TimeInterval), intent(in):: timestep
    ! Local Variables
    integer                            :: s
    type(ESMF_CALKIND_FLAG)            :: calkind

    integer                            :: rc

    call ESMF_TimeIntervalGet(timestep, s=s, calkindflag=calkind, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    esmf2fms_timestep = set_time(s, 0)

  end function esmf2fms_timestep

  function fms2esmf_time(time, calkind)
    ! Return Value
    type(ESMF_Time)                    :: fms2esmf_time
    ! Input Arguments
    type(Time_type), intent(in)        :: time
    type(ESMF_CALKIND_FLAG), intent(in), optional :: calkind
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
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end function fms2esmf_time

  function string_to_date(string, rc)
    character(len=15), intent(in)           :: string
    integer, intent(out), optional          :: rc
    type(time_type)                         :: string_to_date

    integer                                 :: yr,mon,day,hr,min,sec
    
    if(present(rc)) rc = ESMF_SUCCESS

    read(string, '(I4.4,I2.2,I2.2,".",I2.2,I2.2,I2.2)') yr, mon, day, hr, min, sec
    string_to_date = set_date(yr, mon, day, hr, min, sec)

  end function string_to_date

end module time_utils_mod
