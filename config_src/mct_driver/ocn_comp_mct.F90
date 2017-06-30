module ocn_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: ocn_comp_mct
! !INTERFACE:

! !DESCRIPTION:
!  This is the main driver for the MOM6 in CIME
!
! !REVISION HISTORY:
!
! !USES:
   use esmf
   use seq_cdata_mod
   use mct_mod

  ! From MOM6
  use ocean_model_mod, only: ocean_state_type, ocean_public_type
  use ocean_model_mod, only: ocean_model_init
  use MOM_time_manager, only: time_type, set_date, set_calendar_type, THIRTY_DAY_MONTHS
  use MOM_domains, only: MOM_infra_init
!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: ocn_init_mct
  public :: ocn_run_mct
  public :: ocn_final_mct
  private                              ! By default make data private

!
! ! PUBLIC DATA:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
! !PRIVATE MODULE FUNCTIONS:

!
! !PRIVATE MODULE VARIABLES
 type(ocean_state_type), pointer  :: ocn_state => NULL()   ! Private state of ocean
 type(ocean_public_type), pointer :: ocn_surface => NULL() ! Public surface state of ocean

!=======================================================================

contains

!***********************************************************************
!BOP
!
! !IROUTINE: ocn_init_mct
!
! !INTERFACE:
  subroutine ocn_init_mct( EClock, cdata_o, x2o_o, o2x_o, NLFilename )
!
! !DESCRIPTION:
! Initialize POP 
!
! !INPUT/OUTPUT PARAMETERS:

  type(ESMF_Clock)            , intent(inout) :: EClock
  type(seq_cdata)             , intent(inout) :: cdata_o
  type(mct_aVect)             , intent(inout) :: x2o_o, o2x_o
  character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!EOP
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
  type(time_type)     :: time_init ! Start time of coupled model's calendar
  type(time_type)     :: time_in   ! Start time for ocean model at initialization
  type(ESMF_time)     :: current_time
  integer             :: year, month, day, hour, minute, seconds, rc
  character(len=128)  :: errMsg
  integer             :: mpicom

  mpicom = cdata_o%mpicom

  call MOM_infra_init(mpicom)

  call ESMF_ClockGet(EClock, currTime=current_time, rc=rc) 

  call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)

  call set_calendar_type(THIRTY_DAY_MONTHS)

  time_init = set_date(year, month, day, hour, minute, seconds, err_msg=errMsg)

  allocate(ocn_surface)

  !call ocean_model_init(ocn_surface, ocn_state, time_init, time_in)

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_init_mct

!***********************************************************************
!BOP
!
! !IROUTINE: ocn_run_mct
!
! !INTERFACE:
  subroutine ocn_run_mct( EClock, cdata_o, x2o_o, o2x_o)
!
! !DESCRIPTION:
! Run POP for a coupling interval
!
! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o

!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!EOP
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_run_mct

!***********************************************************************
!BOP
!
! !IROUTINE: ocn_final_mct
!
! !INTERFACE:
  subroutine ocn_final_mct( EClock, cdata_o, x2o_o, o2x_o)
!
! !DESCRIPTION:
! Finalize POP
!
! !USES:
! !ARGUMENTS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o
!
! !REVISION HISTORY:
! Author: Fei Liu
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

  end subroutine ocn_final_mct


end module ocn_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
