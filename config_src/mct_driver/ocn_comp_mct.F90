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
  use ESMF,                only: ESMF_clock, ESMF_time, ESMF_timeInterval
  use ESMF,                only: ESMF_ClockGet, ESMF_TimeGet, ESMF_TimeIntervalGet
  use seq_cdata_mod,       only: seq_cdata
  use seq_cdata_mod,       only: seq_cdata_setptrs
  use mct_mod,             only: mct_gsMap, mct_gsmap_init, mct_gsMap_lsize, mct_gsmap_orderedpoints
  use mct_mod,             only: mct_aVect, mct_aVect_init, mct_aVect_zero, mct_aVect_nRattr
  use mct_mod,             only: mct_gGrid, mct_gGrid_init, mct_gGrid_importRAttr, mct_gGrid_importIAttr
  use seq_flds_mod,        only: seq_flds_x2o_fields,            &
                                seq_flds_o2x_fields,            &
                                SEQ_FLDS_DOM_COORD,             &
                                SEQ_FLDS_DOM_other
  use seq_infodata_mod,    only: seq_infodata_type,              &
                                seq_infodata_GetData,           &
                                seq_infodata_start_type_start,  &
                                seq_infodata_start_type_cont,   &
                                seq_infodata_start_type_brnch,  &
                                seq_infodata_PutData
  use seq_comm_mct,        only: seq_comm_name, seq_comm_inst, seq_comm_suffix
  use seq_timemgr_mod,     only: seq_timemgr_EClockGetData, seq_timemgr_RestartAlarmIsOn
  use perf_mod,            only: t_startf, t_stopf
  use shr_kind_mod,        only: SHR_KIND_R8


  ! From MOM6
  use ocean_model_mod,    only: ocean_state_type, ocean_public_type, ocean_model_init_sfc
  use ocean_model_mod,    only: ocean_model_init, get_state_pointers
  use ocean_model_mod,    only: ice_ocean_boundary_type, update_ocean_model
  use MOM_domains,        only: MOM_infra_init, num_pes, root_pe, pe_here
  use MOM_grid,           only: ocean_grid_type, get_global_grid_size
  use MOM_variables,      only: surface
  use MOM_error_handler,  only: MOM_error, FATAL, is_root_pe
  use MOM_time_manager,   only: time_type, set_date, set_calendar_type, NOLEAP
  use coupler_indices,    only: coupler_indices_init, cpl_indices
  use coupler_indices,    only: ocn_export

!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: ocn_init_mct
  public :: ocn_run_mct
  public :: ocn_final_mct
  private                              ! By default make data private
  logical, parameter :: debug=.true.

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
  type MCT_MOM_Data
    type(ocean_state_type), pointer  :: ocn_state => NULL()   !< Private state of ocean
    type(ocean_public_type), pointer :: ocn_public => NULL()  !< Public state of ocean
    type(ocean_grid_type), pointer   :: grid => NULL()        !< A pointer to a grid structure
    type(surface), pointer           :: ocn_surface => NULL() !< A pointer to the ocean surface state
    type(ice_ocean_boundary_type)    :: ice_ocean_boundary    !< A pointer to the ice ocean boundary type
    type(seq_infodata_type), pointer :: infodata

    type(cpl_indices), public :: ind !< Variable IDs

  end type
  type(MCT_MOM_Data) :: glb

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

  type(ESMF_Clock),             intent(inout) :: EClock  !< Time and time step ? \todo Why must this be intent(inout)?
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
  type(ESMF_timeInterval) :: time_interval
  integer             :: year, month, day, hour, minute, seconds, seconds_n, seconds_d, rc
  character(len=384)  :: runid
  character(len=384)  :: runtype
  character(len=32)   :: starttype          ! infodata start type
  integer             :: mpicom_ocn
  integer             :: npes, pe0
  integer             :: i, errorCode
  integer             :: lsize, nsend, nrecv
  logical             :: ldiag_cpl = .false.
  integer             :: ni, nj
  integer             :: isc, iec, jsc, jec     !< Indices for the start and end of the domain
                                                !! in the x and y dir., respectively.
  ! mct variables (these are local for now)
  integer                   :: MOM_MCT_ID
  type(mct_gsMap), pointer  :: MOM_MCT_gsMap => NULL() ! 2d, points to cdata
  type(mct_gGrid), pointer  :: MOM_MCT_dom => NULL()   ! 2d, points to cdata
  type(mct_gsMap)           :: MOM_MCT_gsMap3d  ! for 3d streams, local
  type(mct_gGrid)           :: MOM_MCT_dom3d    ! for 3d streams, local

  ! time management
  integer                   :: ocn_cpl_dt
  real (kind=8)             :: mom_cpl_dt
  real (kind=8), parameter  ::        &
      seconds_in_minute =    60.0d0, &
      seconds_in_hour   =  3600.0d0, &
      seconds_in_day    = 86400.0d0, &
      minutes_in_hour   =    60.0d0


  ! instance control vars (these are local for now)
  integer(kind=4)     :: inst_index
  character(len=16)   :: inst_name
  character(len=16)   :: inst_suffix

  !!!DANGER!!!: change the following vars with the corresponding MOM6 vars
  integer :: km=1 ! number of vertical levels
  integer :: nx_block=0, ny_block=0 ! size of block domain in x,y dir including ghost cells
  integer :: max_blocks_clinic=0 !max number of blocks per processor in each distribution
  integer :: ncouple_per_day = 48
  logical :: lsend_precip_fact ! if T,send precip_fact to cpl for use in fw balance
                               ! (partially-coupled option)
  character(len=128) :: err_msg

!-----------------------------------------------------------------------

  ! set (actually, get from mct) the cdata pointers:
  call seq_cdata_setptrs(cdata_o, id=MOM_MCT_ID, mpicom=mpicom_ocn, &
                         gsMap=MOM_MCT_gsMap, dom=MOM_MCT_dom, infodata=glb%infodata)

  !---------------------------------------------------------------------
  ! Initialize the model run
  !---------------------------------------------------------------------

  call coupler_indices_init(glb%ind)

  call seq_infodata_GetData( glb%infodata, case_name=runid )

  call seq_infodata_GetData( glb%infodata, start_type=starttype)

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


  time_init = set_date(year, month, day, hour, minute, seconds, err_msg=err_msg)
  time_in = set_date(year, month, day, hour, minute, seconds, err_msg=err_msg)

  ! Debugging clocks
  if (debug .and. is_root_pe()) then
    write(6,*) 'ocn_init_mct, current time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, StartTime=current_time, rc=rc)
    call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(6,*) 'ocn_init_mct, start time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, StopTime=current_time, rc=rc)
    call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(6,*) 'ocn_init_mct, stop time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, PrevTime=current_time, rc=rc)
    call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(6,*) 'ocn_init_mct, previous time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, TimeStep=time_interval, rc=rc)
    call ESMF_TimeIntervalGet(time_interval, yy=year, mm=month, d=day, s=seconds, sn=seconds_n, sd=seconds_d, rc=rc)
    write(6,*) 'ocn_init_mct, time step: y,m,d-',year,month,day,'s,sn,sd=',seconds,seconds_n,seconds_d
  endif

  npes = num_pes()
  pe0 = root_pe()

  allocate(glb%ocn_public)
  glb%ocn_public%is_ocean_PE = .true.
  allocate(glb%ocn_public%pelist(npes))
  glb%ocn_public%pelist(:) = (/(i,i=pe0,pe0+npes)/)
  ! \todo Set other bits of glb$ocn_public

  ! Initialize the MOM6 model
  call ocean_model_init(glb%ocn_public, glb%ocn_state, time_init, time_in)

  ! Initialize ocn_state%state out of sight
  call ocean_model_init_sfc(glb%ocn_state, glb%ocn_public)

  ! store pointers to components inside MOM
  call get_state_pointers(glb%ocn_state, grid=glb%grid)

  call t_stopf('MOM_init')


  !---------------------------------------------------------------------
  ! Initialize MCT attribute vectors and indices
  !---------------------------------------------------------------------

  call t_startf('MOM_mct_init')

  if (debug .and. root_pe().eq.pe_here()) print *, "calling ocn_SetGSMap_mct"

  ! Set mct global seg maps:

  call ocn_SetGSMap_mct(mpicom_ocn, MOM_MCT_ID, MOM_MCT_GSMap, MOM_MCT_GSMap3d)
  lsize = mct_gsMap_lsize(MOM_MCT_gsmap, mpicom_ocn)

  ! Initialize mct ocn domain (needs ocn initialization info)

  if (debug .and. root_pe().eq.pe_here()) print *, "calling ocn_domain_mct"
  call ocn_domain_mct(lsize, MOM_MCT_gsmap, MOM_MCT_dom)
  call ocn_domain_mct(lsize*km, MOM_MCT_gsmap3d, MOM_MCT_dom3d) !TODO: this is not used

  ! Inialize mct attribute vectors

  if (debug .and. root_pe().eq.pe_here()) print *, "calling mct_avect_init a"

  ! Initialize the mct attribute vector x2o_o, given Attribute list and length:
  call mct_aVect_init(x2o_o, rList=seq_flds_x2o_fields, lsize=lsize)
  ! set the mct attribute vector x2o_o to zero:
  call mct_aVect_zero(x2o_o)

  if (debug .and. root_pe().eq.pe_here()) print *, "calling mct_avect_init b"

  ! Initialize the mct attribute vector o2x_o, given Attribute list and length:
  call mct_aVect_init(o2x_o, rList=seq_flds_o2x_fields, lsize=lsize)
  ! set the mct attribute vector o2x_o to zero:
  call mct_aVect_zero(o2x_o)

  ! allocate send buffer
  nsend = mct_avect_nRattr(o2x_o)
  nrecv = mct_avect_nRattr(x2o_o)

  ! initialize necessary coupling info

  if (debug .and. root_pe().eq.pe_here()) print *, "calling seq_timemgr_eclockgetdata"

  call seq_timemgr_EClockGetData(EClock, dtime=ocn_cpl_dt)

  ! \todo Need interface to get dt from MOM6
  mom_cpl_dt = seconds_in_day / ncouple_per_day
  if (mom_cpl_dt /= ocn_cpl_dt) then
     write(*,*) 'ERROR pop_cpl_dt and ocn_cpl_dt must be identical'
     call exit(0)
  end if

  ! send initial state to driver

  !TODO:
  ! if ( lsend_precip_fact )  then
  !    call seq_infodata_PutData( infodata, precip_fact=precip_fact)
  ! end if


  if (debug .and. root_pe().eq.pe_here()) print *, "calling ocn_export"
  call ocn_export(glb%ind, glb%ocn_public, glb%grid, o2x_o%rattr)

  call t_stopf('MOM_mct_init')

  if (debug .and. root_pe().eq.pe_here()) print *, "calling get_state_pointers"

  ! Size of global domain
  call get_global_grid_size(glb%grid, ni, nj)

  ! allocate ice_ocean_boundary
  isc = glb%grid%isc; iec = glb%grid%iec;  
  jsc = glb%grid%jsc; jec = glb%grid%jec;  
  allocate(glb%ice_ocean_boundary%u_flux(isc:iec,jsc:jec));          glb%ice_ocean_boundary%u_flux(:,:) = 0.0
  allocate(glb%ice_ocean_boundary%v_flux(isc:iec,jsc:jec));          glb%ice_ocean_boundary%v_flux(:,:) = 0.0
  allocate(glb%ice_ocean_boundary%t_flux(isc:iec,jsc:jec));          glb%ice_ocean_boundary%t_flux(:,:) = 0.0
  allocate(glb%ice_ocean_boundary%q_flux(isc:iec,jsc:jec));          glb%ice_ocean_boundary%q_flux(:,:) = 0.0
  allocate(glb%ice_ocean_boundary%salt_flux(isc:iec,jsc:jec));       glb%ice_ocean_boundary%salt_flux(:,:) = 0.0
  allocate(glb%ice_ocean_boundary%lw_flux(isc:iec,jsc:jec));         glb%ice_ocean_boundary%lw_flux(:,:) = 0.0
  allocate(glb%ice_ocean_boundary%sw_flux_vis_dir(isc:iec,jsc:jec)); glb%ice_ocean_boundary%sw_flux_vis_dir(:,:) = 0.0
  allocate(glb%ice_ocean_boundary%sw_flux_vis_dif(isc:iec,jsc:jec)); glb%ice_ocean_boundary%sw_flux_vis_dif(:,:) = 0.0
  allocate(glb%ice_ocean_boundary%sw_flux_nir_dir(isc:iec,jsc:jec)); glb%ice_ocean_boundary%sw_flux_nir_dir(:,:) = 0.0
  allocate(glb%ice_ocean_boundary%sw_flux_nir_dif(isc:iec,jsc:jec)); glb%ice_ocean_boundary%sw_flux_nir_dif(:,:) = 0.0
  allocate(glb%ice_ocean_boundary%lprec(isc:iec,jsc:jec));           glb%ice_ocean_boundary%lprec(:,:) = 0.0
  allocate(glb%ice_ocean_boundary%fprec(isc:iec,jsc:jec));           glb%ice_ocean_boundary%fprec(:,:) = 0.0
  allocate(glb%ice_ocean_boundary%runoff(isc:iec,jsc:jec));          glb%ice_ocean_boundary%runoff(:,:) = 0.0
  allocate(glb%ice_ocean_boundary%calving(isc:iec,jsc:jec));         glb%ice_ocean_boundary%calving(:,:) = 0.0
  allocate(glb%ice_ocean_boundary%runoff_hflx(isc:iec,jsc:jec));     glb%ice_ocean_boundary%runoff_hflx(:,:) = 0.0
  allocate(glb%ice_ocean_boundary%calving_hflx(isc:iec,jsc:jec));    glb%ice_ocean_boundary%calving_hflx(:,:) = 0.0
  allocate(glb%ice_ocean_boundary%p(isc:iec,jsc:jec));               glb%ice_ocean_boundary%p(:,:) = 0.0
  allocate(glb%ice_ocean_boundary%mi(isc:iec,jsc:jec));              glb%ice_ocean_boundary%mi(:,:) = 0.0
 
 
  if (debug .and. root_pe().eq.pe_here()) print *, "calling seq_infodata_putdata"

   call seq_infodata_PutData( glb%infodata, &
        ocn_nx = ni , ocn_ny = nj)
   call seq_infodata_PutData( glb%infodata, &
        ocn_prognostic=.true., ocnrof_prognostic=.true.)


  if (debug .and. root_pe().eq.pe_here()) print *, "leaving ocean_init_mct"

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_init_mct

  !> Step forward ocean model for coupling interval
  subroutine ocn_run_mct( EClock, cdata_o, x2o_o, o2x_o)
  type(ESMF_Clock), intent(inout) :: EClock  !< Time and time step ? \todo Why must this be intent(inout)?
  type(seq_cdata),  intent(inout) :: cdata_o
  type(mct_aVect),  intent(inout) :: x2o_o
  type(mct_aVect),  intent(inout) :: o2x_o
  ! Local variables
  type(ESMF_time) :: current_time
  type(ESMF_timeInterval) :: time_interval
  integer :: year, month, day, hour, minute, seconds, seconds_n, seconds_d, rc
  logical :: write_restart_at_eod
  type(time_type) :: time_start ! Start of coupled time interval to pass to MOM6
  type(time_type) :: coupling_timestep ! Coupled time interval to pass to MOM6
  character(len=128) :: err_msg

  ! Translate the current time (start of coupling interval)
  call ESMF_ClockGet(EClock, currTime=current_time, rc=rc)
  call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
  time_start = set_date(year, month, day, hour, minute, seconds, err_msg=err_msg)

  ! Debugging clocks
  if (debug .and. is_root_pe()) then
    write(6,*) 'ocn_run_mct, current time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, StartTime=current_time, rc=rc)
    call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(6,*) 'ocn_run_mct, start time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, StopTime=current_time, rc=rc)
    call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(6,*) 'ocn_run_mct, stop time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, PrevTime=current_time, rc=rc)
    call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(6,*) 'ocn_run_mct, previous time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
  endif

  ! Translate the coupling time interval
  call ESMF_ClockGet(EClock, TimeStep=time_interval, rc=rc)
  call ESMF_TimeIntervalGet(time_interval, yy=year, mm=month, d=day, s=seconds, sn=seconds_n, sd=seconds_d, rc=rc)
  time_start = set_date(year, month, day, 0, 0, seconds, err_msg=err_msg)
  if (debug .and. is_root_pe()) then
    write(6,*) 'ocn_run_mct, time step: y,m,d-',year,month,day,'s,sn,sd=',seconds,seconds_n,seconds_d
  endif

  ! set (actually, get from mct) the cdata pointers:
  ! \todo this was done in _init_, is it needed again. Does this infodata need to be in glb%?
  call seq_cdata_setptrs(cdata_o, infodata=glb%infodata)

  ! Check alarms for flag to write restart at end of day
  write_restart_at_eod = seq_timemgr_RestartAlarmIsOn(EClock)
  ! \todo Let MOM6 know to write restart...
  if (debug .and. is_root_pe()) write(6,*) 'ocn_run_mct, write_restart_at_eod=', write_restart_at_eod

  ! fill ice ocean boundary
  call fill_ice_ocean_bnd(glb%ice_ocean_boundary, glb%grid, x2o_o%rattr, glb%ind)
  if (debug .and. is_root_pe()) write(6,*) 'fill_ice_ocean_bnd'

!  call update_ocean_model(glb%ice_ocean_boundary, glb%ocn_state, glb%ocn_public, &
!                          time_start, coupling_timestep)

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


!> This routine mct global seg maps for the MOM decomposition
!!
!! \todo Find out if we should only provide indirect indexing for ocean points and not land.
subroutine ocn_SetGSMap_mct(mpicom_ocn, MOM_MCT_ID, gsMap_ocn, gsMap3d_ocn)
  integer,         intent(in)    :: mpicom_ocn  !< MPI communicator
  integer,         intent(in)    :: MOM_MCT_ID  !< MCT component ID
  type(mct_gsMap), intent(inout) :: gsMap_ocn   !< MCT global segment map for 2d data
  type(mct_gsMap), intent(inout) :: gsMap3d_ocn !< MCT global segment map for 3d data
  ! Local variables
  integer :: lsize ! Local size of indirect indexing array
  integer :: i, j, k ! Local indices
  integer :: ni, nj ! Declared sizes of h-point arrays
  integer :: ig, jg ! Global indices
  type(ocean_grid_type), pointer :: grid => NULL() ! A pointer to a grid structure
  integer, allocatable :: gindex(:) ! Indirect indices

  grid => glb%grid ! for convenience
  if (.not. associated(grid)) call MOM_error(FATAL, 'ocn_comp_mct.F90, ocn_SetGSMap_mct():' // &
      'grid returned from get_state_pointers() was not associated!')

  ! Size of computational domain
  lsize = ( grid%iec - grid%isc + 1 ) * ( grid%jec - grid%jsc + 1 )

  ! Size of global domain
  call get_global_grid_size(grid, ni, nj)

  ! Create indirect indices for the computational domain
  allocate( gindex( lsize ) )

  ! Set indirect indices in gindex
  k = 0
  do j = grid%jsc, grid%jec
    jg = j + grid%jdg_offset ! TODO: check this calculation
    do i = grid%isc, grid%iec
      ig = i + grid%idg_offset ! TODO: check this calculation
      k = k + 1 ! Increment position within gindex
      gindex(k) = ni * ( jg - 1 ) + ig
    enddo
  enddo

  ! Tell MCT how to indirectly index into the 2d buffer
  call mct_gsMap_init( gsMap_ocn, gindex, mpicom_ocn, MOM_MCT_ID, lsize, ni * nj)

  deallocate( gindex )

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

! Local Variables
  integer, parameter              :: SHR_REAL_R8 = selected_real_kind(12)
  integer, pointer                :: idata(:)
  integer                         :: i,j,k
  real(kind=SHR_REAL_R8), pointer :: data(:)
  real(kind=SHR_REAL_R8)          :: m2_to_rad2
  type(ocean_grid_type), pointer :: grid => NULL() ! A pointer to a grid structure

  grid => glb%grid ! for convenience

  ! set coords to lat and lon, and areas to rad^2
  call mct_gGrid_init(GGrid=dom_ocn, CoordChars=trim(seq_flds_dom_coord), &
                      OtherChars=trim(seq_flds_dom_other), lsize=lsize )

  call mct_avect_zero(dom_ocn%data)
  allocate(data(lsize))

  ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
  k = pe_here()
  call mct_gsMap_orderedPoints(gsMap_ocn, k, idata)
  call mct_gGrid_importIAttr(dom_ocn,'GlobGridNum',idata,lsize)

  !initialization
  data(:) = -9999.0
  call mct_gGrid_importRAttr(dom_ocn,"lat"  ,data,lsize)
  call mct_gGrid_importRAttr(dom_ocn,"lon"  ,data,lsize)
  call mct_gGrid_importRAttr(dom_ocn,"area" ,data,lsize)
  call mct_gGrid_importRAttr(dom_ocn,"aream",data,lsize)
  data(:) = 0.0
  call mct_gGrid_importRAttr(dom_ocn,"mask",data,lsize)
  call mct_gGrid_importRAttr(dom_ocn,"frac",data,lsize)

  k = 0
  do j = grid%jsc, grid%jec
    do i = grid%isc, grid%iec
      k = k + 1 ! Increment position within gindex
      data(k) = grid%geoLonT(i,j)
    enddo
  enddo
  call mct_gGrid_importRattr(dom_ocn,"lon",data,lsize)

  k = 0
  do j = grid%jsc, grid%jec
    do i = grid%isc, grid%iec
      k = k + 1 ! Increment position within gindex
      data(k) = grid%geoLatT(i,j)
    enddo
  enddo
  call mct_gGrid_importRattr(dom_ocn,"lat",data,lsize)

  k = 0
  m2_to_rad2 = 1./grid%Rad_Earth**2
  do j = grid%jsc, grid%jec
    do i = grid%isc, grid%iec
      k = k + 1 ! Increment position within gindex
      data(k) = grid%AreaT(i,j) * m2_to_rad2
    enddo
  enddo
  call mct_gGrid_importRattr(dom_ocn,"area",data,lsize)

  k = 0
  do j = grid%jsc, grid%jec
    do i = grid%isc, grid%iec
      k = k + 1 ! Increment position within gindex
      data(k) = grid%mask2dT(i,j)
    enddo
  enddo
  call mct_gGrid_importRattr(dom_ocn,"mask",data,lsize)
  call mct_gGrid_importRattr(dom_ocn,"frac",data,lsize)

  deallocate(data)
  deallocate(idata)

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
