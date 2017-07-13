module ocn_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: ocn_comp_mct
! !INTERFACE:

! !DESCRIPTION:
!  This is the main driver for MOM6 in CIME
!
! !REVISION HISTORY:
!
! !USES:
  use esmf
  use seq_cdata_mod
  use mct_mod
  use seq_flds_mod,       only: seq_flds_x2o_fields,            &
                                seq_flds_o2x_fields
  use seq_infodata_mod,   only: seq_infodata_type,              &
                                seq_infodata_GetData,           &
                                seq_infodata_start_type_start,  &
                                seq_infodata_start_type_cont,   &
                                seq_infodata_start_type_brnch
  use seq_comm_mct,       only: seq_comm_name, seq_comm_inst, seq_comm_suffix
  use perf_mod,           only: t_startf, t_stopf
   
 
  ! From MOM6
  use ocean_model_mod,    only: ocean_state_type, ocean_public_type
  use ocean_model_mod,    only: ocean_model_init
  use MOM_time_manager,   only: time_type, set_date, set_calendar_type, NOLEAP
  use MOM_domains,        only: MOM_infra_init, num_pes, root_pe, pe_here
  use coupler_indices,    only: coupler_indices_init
  use ocn_import_export,  only: SBUFF_SUM

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
  private :: ocn_SetGSMap_mct
  private :: ocn_domain_mct

! !PRIVATE MODULE VARIABLES
  type(ocean_state_type), pointer  :: ocn_state => NULL()   ! Private state of ocean
  type(ocean_public_type), pointer :: ocn_surface => NULL() ! Public surface state of ocean

  type(seq_infodata_type), pointer :: &
      infodata

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
  character(len=384)  :: runid
  character(len=384)  :: runtype
  character(len=32)   :: starttype          ! infodata start type
  integer             :: mpicom_ocn
  integer             :: npes, pe0
  integer             :: i
  integer             :: lsize, nsend, nrecv

  ! mct variables (these are local for now)
  integer                   :: MOM_MCT_ID
  type(mct_gsMap), pointer  :: MOM_MCT_gsMap    ! 2d, points to cdata
  type(mct_gGrid), pointer  :: MOM_MCT_dom      ! 2d, points to cdata
  type(mct_gsMap)           :: MOM_MCT_gsMap3d  ! for 3d streams, local
  type(mct_gGrid)           :: MOM_MCT_dom3d    ! for 3d streams, local

  ! instance control vars (these are local for now)
  integer(kind=4)     :: inst_index
  character(len=16)   :: inst_name
  character(len=16)   :: inst_suffix

  !!!DANGER!!!: change the following vars with the corresponding MOM6 vars
  integer :: km=62 ! number of vertical levels 

!-----------------------------------------------------------------------

  ! set (actually, get from mct) the cdata pointers:
  call seq_cdata_setptrs(cdata_o, id=MOM_MCT_ID, mpicom=mpicom_ocn, infodata=infodata)

  !---------------------------------------------------------------------
  ! Initialize the model run 
  !---------------------------------------------------------------------

  call coupler_indices_init()

  call seq_infodata_GetData( infodata, case_name=runid )

  call seq_infodata_GetData( infodata, start_type=starttype)

  if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
     runtype = "initial"
  else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
     runtype = "continue"
  else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
     runtype = "branch"
  else
     write(*,*) 'ocn_comp_mct ERROR: unknown starttype'
     call exit(0)
  end if

  ! instance control

  inst_name   = seq_comm_name(MOM_MCT_ID)
  inst_index  = seq_comm_inst(MOM_MCT_ID)
  inst_suffix = seq_comm_suffix(MOM_MCT_ID)

  !---------------------------------------------------------------------
  ! Initialize MOM6 
  !---------------------------------------------------------------------

  call t_startf('MOM_init')

  call MOM_infra_init(mpicom_ocn)

  call ESMF_ClockGet(EClock, currTime=current_time, rc=rc) 
  call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
  call set_calendar_type(NOLEAP)  !TODO: confirm this

  time_init = set_date(year, month, day, hour, minute, seconds, err_msg=errMsg)
  time_in = set_date(year, month, day, hour, minute, seconds, err_msg=errMsg)

  npes = num_pes()
  pe0 = root_pe()

  allocate(ocn_surface)
  ocn_surface%is_ocean_PE = .true.
  allocate(ocn_surface%pelist(npes))
  ocn_surface%pelist(:) = (/(i,i=pe0,pe0+npes)/)

  ! initialize the MOM6 model
  call ocean_model_init(ocn_surface, ocn_state, time_init, time_in)

  call t_stopf('MOM_init')

  !---------------------------------------------------------------------
  ! Initialize MCT attribute vectors and indices 
  !---------------------------------------------------------------------

  call t_startf('MOM_mct_init')

  ! Set mct global seg maps:

  call ocn_SetGSMap_mct(mpicom_ocn, MOM_MCT_ID, MOM_MCT_GSMap, MOM_MCT_GSMap3d) 
  lsize = mct_gsMap_lsize(MOM_MCT_gsmap, mpicom_ocn) 

  ! Initialize mct ocn domain (needs ocn initialization info)

  call ocn_domain_mct(lsize, MOM_MCT_gsmap, MOM_MCT_dom)
  call ocn_domain_mct(lsize*km, MOM_MCT_gsmap3d, MOM_MCT_dom3d)

  ! Inialize mct attribute vectors
  
  ! Initialize the mct attribute vector x2o_o, given Attribute list and length:
  call mct_aVect_init(x2o_o, rList=seq_flds_x2o_fields, lsize=lsize)
  ! set the mct attribute vector x2o_o to zero:
  call mct_aVect_zero(x2o_o)
    
  ! Initialize the mct attribute vector o2x_o, given Attribute list and length:
  call mct_aVect_init(o2x_o, rList=seq_flds_o2x_fields, lsize=lsize)
  ! set the mct attribute vector o2x_o to zero:
  call mct_aVect_zero(o2x_o)

  nsend = mct_avect_nRattr(o2x_o)
  nrecv = mct_avect_nRattr(x2o_o)    
  !allocate (SBUFF_SUM(nx_block,ny_block,max_blocks_clinic,nsend))















  call t_stopf('MOM_mct_init')



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


!***********************************************************************
!BOP
!IROUTINE: ocn_SetGSMap_mct
! !INTERFACE:

  subroutine ocn_SetGSMap_mct(mpicom_ocn, MOM_MCT_ID, gsMap_ocn, gsMap3d_ocn)

! !DESCRIPTION:
!  This routine mct global seg maps for the MOM decomposition
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

    implicit none
    integer        , intent(in)    :: mpicom_ocn
    integer        , intent(in)    :: MOM_MCT_ID
    type(mct_gsMap), intent(inout) :: gsMap_ocn
    type(mct_gsMap), intent(inout) :: gsMap3d_ocn

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


 
!-----------------------------------------------------------------------
!EOC
  
  end subroutine ocn_SetGSMap_mct


!***********************************************************************
!BOP
! !IROUTINE: ocn_domain_mct
! !INTERFACE:

 subroutine ocn_domain_mct( lsize, gsMap_ocn, dom_ocn)

! !DESCRIPTION:
!  This routine mct global seg maps for the pop decomposition
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT/OUTPUT PARAMETERS:

    implicit none
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(in)    :: gsMap_ocn
    type(mct_ggrid), intent(inout) :: dom_ocn

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_domain_mct


end module ocn_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
