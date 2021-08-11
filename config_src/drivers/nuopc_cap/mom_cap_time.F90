!> This was originally share code in CIME, but required CIME as a
!! dependency to build the MOM cap.  The options here for setting
!! a restart alarm are useful for all caps, so a second step is to
!! determine if/how these could be offered more generally in a
!! shared library.  For now we really want the MOM cap to only
!! depend on MOM and ESMF/NUOPC.
module MOM_cap_time

! !USES:
use ESMF                  , only : ESMF_Time, ESMF_Clock, ESMF_Calendar, ESMF_Alarm
use ESMF                  , only : ESMF_TimeGet, ESMF_TimeSet
use ESMF                  , only : ESMF_TimeInterval, ESMF_TimeIntervalSet
use ESMF                  , only : ESMF_ClockGet, ESMF_AlarmCreate
use ESMF                  , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO
use ESMF                  , only : ESMF_LogSetError, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
use ESMF                  , only : ESMF_RC_ARG_BAD
use ESMF                  , only : operator(<), operator(/=), operator(+), operator(-), operator(*) , operator(>=)
use ESMF                  , only : operator(<=), operator(>), operator(==)
use MOM_cap_methods       , only : ChkErr

implicit none; private

public  :: AlarmInit  ! initialize an alarm

private :: TimeInit
private :: date2ymd

! Clock and alarm options
character(len=*), private, parameter :: &
     optNONE           = "none"      , &
     optNever          = "never"     , &
     optNSteps         = "nsteps"    , &
     optNStep          = "nstep"     , &
     optNSeconds       = "nseconds"  , &
     optNSecond        = "nsecond"   , &
     optNMinutes       = "nminutes"  , &
     optNMinute        = "nminute"   , &
     optNHours         = "nhours"    , &
     optNHour          = "nhour"     , &
     optNDays          = "ndays"     , &
     optNDay           = "nday"      , &
     optNMonths        = "nmonths"   , &
     optNMonth         = "nmonth"    , &
     optNYears         = "nyears"    , &
     optNYear          = "nyear"     , &
     optMonthly        = "monthly"   , &
     optYearly         = "yearly"    , &
     optDate           = "date"      , &
     optIfdays0        = "ifdays0"   , &
     optGLCCouplingPeriod = "glc_coupling_period"

! Module data
integer, parameter          :: SecPerDay = 86400 ! Seconds per day
character(len=*), parameter :: u_FILE_u = &
     __FILE__

contains

!> Setup an alarm in a clock. The ringtime sent to AlarmCreate
!! MUST be the next alarm time.  If you send an arbitrary but
!! proper ringtime from the past and the ring interval, the alarm
!! will always go off on the next clock advance and this will cause
!! serious problems. Even if it makes sense to initialize an alarm
!! with some reference time and the alarm interval, that reference
!! time has to be advance forward to be >= the current time.
!! In the logic below  we set an appropriate "NextAlarm" and then
!! we make sure to advance it properly based on the ring interval.
subroutine AlarmInit( clock, alarm, option, &
     opt_n, opt_ymd, opt_tod, RefTime, alarmname, rc)
  type(ESMF_Clock)            , intent(inout) :: clock     !< ESMF clock
  type(ESMF_Alarm)            , intent(inout) :: alarm     !< ESMF alarm
  character(len=*)            , intent(in)    :: option    !< alarm option
  integer          , optional , intent(in)    :: opt_n     !< alarm freq
  integer          , optional , intent(in)    :: opt_ymd   !< alarm ymd
  integer          , optional , intent(in)    :: opt_tod   !< alarm tod (sec)
  type(ESMF_Time)  , optional , intent(in)    :: RefTime   !< ref time
  character(len=*) , optional , intent(in)    :: alarmname !< alarm name
  integer                     , intent(inout) :: rc        !< Return code

  ! local variables
  type(ESMF_Calendar)     :: cal              ! calendar
  integer                 :: lymd             ! local ymd
  integer                 :: ltod             ! local tod
  integer                 :: cyy,cmm,cdd,csec ! time info
  integer                 :: nyy,nmm,ndd,nsec ! time info
  character(len=64)       :: lalarmname       ! local alarm name
  logical                 :: update_nextalarm ! update next alarm
  type(ESMF_Time)         :: CurrTime         ! Current Time
  type(ESMF_Time)         :: NextAlarm        ! Next restart alarm time
  type(ESMF_TimeInterval) :: AlarmInterval    ! Alarm interval
  character(len=*), parameter :: subname = '(AlarmInit): '
  !-------------------------------------------------------------------------------

  rc = ESMF_SUCCESS

  lalarmname = 'alarm_unknown'
  if (present(alarmname)) lalarmname = trim(alarmname)
  ltod = 0
  if (present(opt_tod)) ltod = opt_tod
  lymd = -1
  if (present(opt_ymd)) lymd = opt_ymd

  ! verify parameters
  if (trim(option) == optNSteps    .or. trim(option) == optNStep   .or. &
      trim(option) == optNSeconds  .or. trim(option) == optNSecond .or. &
      trim(option) == optNMinutes  .or. trim(option) == optNMinute .or. &
      trim(option) == optNHours    .or. trim(option) == optNHour   .or. &
      trim(option) == optNDays     .or. trim(option) == optNDay    .or. &
      trim(option) == optNMonths   .or. trim(option) == optNMonth  .or. &
      trim(option) == optNYears    .or. trim(option) == optNYear   .or. &
      trim(option) == optIfdays0) then
     if (.not. present(opt_n)) then
        call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
             msg=subname//trim(option)//' requires opt_n', &
             line=__LINE__, &
             file=__FILE__, rcToReturn=rc)
        return
     endif
     if (opt_n <= 0) then
        call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
             msg=subname//trim(option)//' invalid opt_n', &
             line=__LINE__, &
             file=__FILE__, rcToReturn=rc)
        return
     endif
  endif

  call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_TimeGet(CurrTime, yy=cyy, mm=cmm, dd=cdd, s=csec, rc=rc )
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_TimeGet(CurrTime, yy=nyy, mm=nmm, dd=ndd, s=nsec, rc=rc )
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! initial guess of next alarm, this will be updated below
  if (present(RefTime)) then
     NextAlarm = RefTime
  else
     NextAlarm = CurrTime
  endif

  ! Determine calendar
  call ESMF_ClockGet(clock, calendar=cal, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! Determine inputs for call to create alarm
  selectcase (trim(option))

  case (optNONE, optNever)
     call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=cal, rc=rc )
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     update_nextalarm  = .false.

  case (optDate)
     if (.not. present(opt_ymd)) then
        call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
             msg=subname//trim(option)//' requires opt_ymd', &
             line=__LINE__, &
             file=__FILE__, rcToReturn=rc)
        return
     endif
     if (lymd < 0 .or. ltod < 0) then
        call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
             msg=subname//trim(option)//'opt_ymd, opt_tod invalid', &
             line=__LINE__, &
             file=__FILE__, rcToReturn=rc)
        return
     endif
     call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     call TimeInit(NextAlarm, lymd, cal, tod=ltod, desc="optDate", rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     update_nextalarm  = .false.

  case (optIfdays0)
     if (.not. present(opt_ymd)) then
        call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
             msg=subname//trim(option)//' requires opt_ymd', &
             line=__LINE__, &
             file=__FILE__, rcToReturn=rc)
        return
     endif
     call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=opt_n, s=0, calendar=cal, rc=rc )
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     update_nextalarm  = .true.

  case (optNSteps, optNStep)
     call ESMF_ClockGet(clock, TimeStep=AlarmInterval, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     AlarmInterval = AlarmInterval * opt_n
     update_nextalarm  = .true.

  case (optNSeconds, optNSecond)
     call ESMF_TimeIntervalSet(AlarmInterval, s=1, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     AlarmInterval = AlarmInterval * opt_n
     update_nextalarm  = .true.

  case (optNMinutes, optNMinute)
     call ESMF_TimeIntervalSet(AlarmInterval, s=60, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     AlarmInterval = AlarmInterval * opt_n
     update_nextalarm  = .true.

  case (optNHours, optNHour)
     call ESMF_TimeIntervalSet(AlarmInterval, s=3600, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     AlarmInterval = AlarmInterval * opt_n
     update_nextalarm  = .true.

  case (optNDays, optNDay)
     call ESMF_TimeIntervalSet(AlarmInterval, d=1, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     AlarmInterval = AlarmInterval * opt_n
     update_nextalarm  = .true.

  case (optNMonths, optNMonth)
     call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     AlarmInterval = AlarmInterval * opt_n
     update_nextalarm  = .true.

  case (optMonthly)
     call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=1, s=0, calendar=cal, rc=rc )
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     update_nextalarm  = .true.

  case (optNYears, optNYear)
     call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     AlarmInterval = AlarmInterval * opt_n
     update_nextalarm  = .true.

  case (optYearly)
     call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     call ESMF_TimeSet( NextAlarm, yy=cyy, mm=1, dd=1, s=0, calendar=cal, rc=rc )
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     update_nextalarm  = .true.

  case default
     call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
          msg=subname//' unknown option: '//trim(option), &
          line=__LINE__, &
          file=__FILE__, rcToReturn=rc)
     return

  end select

  ! --------------------------------------------------------------------------------
  ! --- AlarmInterval and NextAlarm should be set ---
  ! --------------------------------------------------------------------------------

  ! --- advance Next Alarm so it won't ring on first timestep for
  ! --- most options above. go back one alarminterval just to be careful

  if (update_nextalarm) then
     NextAlarm = NextAlarm - AlarmInterval
     do while (NextAlarm <= CurrTime)
        NextAlarm = NextAlarm + AlarmInterval
     enddo
  endif

  alarm = ESMF_AlarmCreate( name=lalarmname, clock=clock, ringTime=NextAlarm, ringInterval=AlarmInterval, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

end subroutine AlarmInit

!> Creates the ESMF_Time object corresponding to the given input time,
!! given in YMD (Year Month Day) and TOD (Time-of-day) format. Sets
!! the time by an integer as YYYYMMDD and integer seconds in the day.
subroutine TimeInit( Time, ymd, cal, tod, desc, logunit, rc)
  type(ESMF_Time)     , intent(inout)         :: Time !< ESMF time
  integer             , intent(in)            :: ymd  !< year, month, day YYYYMMDD
  type(ESMF_Calendar) , intent(in)            :: cal  !< ESMF calendar
  integer             , intent(in),  optional :: tod  !< time of day in [sec]
  character(len=*)    , intent(in),  optional :: desc !< description of time to set
  integer             , intent(in),  optional :: logunit!< Unit for stdout output
  integer             , intent(out), optional :: rc   !< Return code

  ! local varaibles
  integer                     :: yr, mon, day ! Year, month, day as integers
  integer                     :: ltod         ! local tod
  character(len=256)          :: ldesc        ! local desc
  character(len=*), parameter :: subname = '(TimeInit) '
  !-------------------------------------------------------------------------------

  ltod = 0
  if (present(tod)) ltod = tod
  ldesc = ''
  if (present(desc)) ldesc = desc

  if ( (ymd < 0) .or. (ltod < 0) .or. (ltod > SecPerDay) )then
     if (present(logunit)) then
        write(logunit,*) subname//': ERROR yymmdd is a negative number or '// &
             'time-of-day out of bounds', ymd, ltod
     endif
     call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
          msg=subname//' yymmdd is negative or time-of-day out of bounds ', &
          line=__LINE__, &
          file=__FILE__, rcToReturn=rc)
     return
  endif

  call date2ymd (ymd,yr,mon,day)

  call ESMF_TimeSet( Time, yy=yr, mm=mon, dd=day, s=ltod, calendar=cal, rc=rc )
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

end subroutine TimeInit

!> Converts a coded-date (yyyymmdd) into calendar year,month,day.
subroutine date2ymd (date, year, month, day)
  integer, intent(in)  :: date             !< coded-date (yyyymmdd)
  integer, intent(out) :: year,month,day   !< calendar year,month,day

  ! local variables
  integer :: tdate   ! temporary date
  character(*),parameter :: subName = "(date2ymd)"
  !-------------------------------------------------------------------------------

  tdate = abs(date)
  year = int(tdate/10000)
  if (date < 0) then
     year = -year
  endif
  month = int( mod(tdate,10000)/  100)
  day = mod(tdate,  100)

end subroutine date2ymd

end module MOM_cap_time
