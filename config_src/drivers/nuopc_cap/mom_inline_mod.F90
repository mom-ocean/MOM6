!> This module contains a set of subroutines that enables inline CDEPS capability

module mom_inline_mod

use NUOPC            , only: NUOPC_CompAttributeGet
use ESMF             , only: ESMF_GridComp, ESMF_Mesh
use ESMF             , only: ESMF_Clock, ESMF_Time, ESMF_TimeGet, ESMF_ClockGet
use ESMF             , only: ESMF_KIND_R8, ESMF_SUCCESS, ESMF_LogFoundError
use ESMF             , only: ESMF_LOGERR_PASSTHRU, ESMF_LOGMSG_INFO, ESMF_LOGWRITE
use ESMF             , only: ESMF_END_ABORT, ESMF_Finalize, ESMF_MAXSTR
use dshr_mod         , only: dshr_pio_init
use dshr_strdata_mod , only: shr_strdata_type, shr_strdata_print
use dshr_strdata_mod , only: shr_strdata_init_from_inline
use dshr_strdata_mod , only: shr_strdata_advance
use dshr_methods_mod , only: dshr_fldbun_getfldptr, dshr_fldbun_Field_diagnose
use dshr_stream_mod  , only: shr_stream_init_from_esmfconfig
use MOM_cap_methods  , only: ChkErr

implicit none
private

public mom_inline_init
public mom_inline_run

type(shr_strdata_type), allocatable :: sdat(:)

integer                    :: logunit      ! the logunit on the root task
character(len=ESMF_MAXSTR) :: stream_name  ! generic identifier

character(len=*), parameter :: u_FILE_u =  __FILE__
contains

!===============================================================================
subroutine mom_inline_init(gcomp, model_clock, model_mesh, mytask, rc)
  type(ESMF_GridComp)    , intent(in)  :: gcomp        !< ESMF_GridComp object
  type(ESMF_Clock)       , intent(in)  :: model_clock  !< ESMF_Clock object
  type(ESMF_Mesh)        , intent(in)  :: model_mesh   !< ESMF mesh
  integer                , intent(in)  :: mytask       !< the current task
  integer                , intent(out) :: rc           !< Return code

  ! local variables
  logical :: isPresent, isSet
  integer :: ns, l
  integer :: nstreams, nvars
  type(shr_strdata_type) :: sdatconfig !< stream data from config (xml or esmf), one or more streams

  character(len=ESMF_MAXSTR) :: value, streamfilename
  character(len=ESMF_MAXSTR), allocatable :: filelist(:)
  character(len=ESMF_MAXSTR), allocatable :: filevars(:,:)

  character(len=*), parameter  :: subname='(mom_inline_init)'
  !----------------------------------------------------------------------

  rc = ESMF_SUCCESS

  call NUOPC_CompAttributeGet(gcomp, name="streamfilename", value=value, isPresent=isPresent, isSet=isSet, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (isPresent .and. isSet) then
    streamfilename = value
  else
    call ESMF_LogWrite(trim(subname)//': streamfilename must be provided', ESMF_LOGMSG_INFO)
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
    return
  endif

#ifndef CESMCOUPLED
  if (mytask == 0) then
    open (newunit=logunit, file='log.mom6.cdeps')
  else
    logunit = 6
  endif

  ! CMEPS Init PIO
  call dshr_pio_init(gcomp, sdatconfig, logunit, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) return

  ! read the available stream definitions, each data stream is one or more data_files
  ! which have the same spatial and temporal coordinates
  call shr_stream_init_from_esmfconfig(trim(streamfilename), sdatconfig%stream, logunit, &
       sdatconfig%pio_subsystem, sdatconfig%io_type, sdatconfig%io_format, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) return
#else
  !TODO: CESM logunit, configuration via xml etc
  !call shr_stream_init_from_xml(trim(streamfilename) ....
#endif

  nstreams = size(sdatconfig%stream)
  ! allocate stream data type
  if (.not. allocated(sdat)) allocate(sdat(nstreams))

  ! set the model clock and mesh
  sdat(:)%model_clock = model_clock
  sdat(:)%model_mesh = model_mesh

  ! loop over streams and initialize
  do ns = 1, nstreams
    sdat(ns)%pio_subsystem => sdatconfig%pio_subsystem
    sdat(ns)%io_type = sdatconfig%io_type
    sdat(ns)%io_format = sdatconfig%io_format

    allocate(filelist(sdatconfig%stream(ns)%nfiles))
    allocate(filevars(sdatconfig%stream(ns)%nvars,2))

    ! fill stream info
    do l = 1, sdatconfig%stream(ns)%nfiles
      filelist(l) = trim(sdatconfig%stream(ns)%file(l)%name)
    enddo
    do l = 1, sdatconfig%stream(ns)%nvars
      filevars(l,1) = trim(sdatconfig%stream(ns)%varlist(l)%nameinfile)
      filevars(l,2) = trim(sdatconfig%stream(ns)%varlist(l)%nameinmodel)
    enddo

    write(stream_name,fmt='(a,i2.2)') 'stream_', ns
    call shr_strdata_init_from_inline(sdat(ns),                         &
         my_task             = mytask,                                  &
         logunit             = logunit,                                 &
         compname            = 'OCN',                                   &
         model_clock         = sdat(ns)%model_clock,                    &
         model_mesh          = sdat(ns)%model_mesh,                     &
         stream_name         = trim(stream_name),                       &
         stream_meshfile     = trim(sdatconfig%stream(ns)%meshfile),    &
         stream_filenames    = filelist,                                &
         stream_yearFirst    = sdatconfig%stream(ns)%yearFirst,         &
         stream_yearLast     = sdatconfig%stream(ns)%yearLast,          &
         stream_yearAlign    = sdatconfig%stream(ns)%yearAlign,         &
         stream_fldlistFile  = filevars(:,1),                           &
         stream_fldListModel = filevars(:,2),                           &
         stream_lev_dimname  = trim(sdatconfig%stream(ns)%lev_dimname), &
         stream_mapalgo      = trim(sdatconfig%stream(ns)%mapalgo),     &
         stream_offset       = sdatconfig%stream(ns)%offset,            &
         stream_taxmode      = trim(sdatconfig%stream(ns)%taxmode),     &
         stream_dtlimit      = sdatconfig%stream(ns)%dtlimit,           &
         stream_tintalgo     = trim(sdatconfig%stream(ns)%tInterpAlgo), &
         stream_src_mask     = sdatconfig%stream(ns)%src_mask_val,      &
         stream_dst_mask     = sdatconfig%stream(ns)%dst_mask_val,      &
         rc                  = rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    deallocate(filelist)
    deallocate(filevars)
  enddo

end subroutine mom_inline_init
!===============================================================================
subroutine mom_inline_run(clock, ocean_public, ocean_grid, ice_ocean_boundary, dbug, rc)
  use MOM_ocean_model_nuopc,     only: ocean_public_type
  use MOM_surface_forcing_nuopc, only: ice_ocean_boundary_type
  use MOM_grid,                  only: ocean_grid_type
  use mpp_domains_mod,           only: mpp_get_compute_domain

  type(ESMF_Clock) ,              intent(in)    :: clock              !< ESMF_Clock object
  type(ocean_public_type)       , intent(in)    :: ocean_public       !< Ocean surface state
  type(ocean_grid_type)         , intent(in)    :: ocean_grid         !< Ocean model grid
  type(ice_ocean_boundary_type) , intent(inout) :: ice_ocean_boundary !< Ocean boundary forcing
  integer ,                       intent(in)    :: dbug               !< Integer debug flag
  integer ,                       intent(out)   :: rc                 !< Return code

  ! local variables
  type(ESMF_Time)             :: date
  integer                     :: nstreams, nflds
  integer                     :: ns,nf,n,i,j
  integer                     :: isc, iec, jsc, jec
  integer                     :: year    ! year (0, ...) for nstep+1
  integer                     :: mon     ! month (1, ..., 12) for nstep+1
  integer                     :: day     ! day of month (1, ..., 31) for nstep+1
  integer                     :: sec     ! seconds into current date for nstep+1
  integer                     :: mcdate  ! Current model date (yyyymmdd)
  character(len=ESMF_MAXSTR)  :: fldname
  real(ESMF_KIND_R8), pointer :: dataPtr1d(:)
  character(len=*), parameter :: subname='(mom_inline_run)'
  !-----------------------------------------------------------------------

  rc = ESMF_SUCCESS

  ! The following are global indices without halos
  call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)

  ! Current model date
  call ESMF_ClockGet( clock, currTime=date, rc=rc )
  if (chkerr(rc,__LINE__,u_FILE_u)) return
  call ESMF_TimeGet(date, yy=year, mm=mon, dd=day, s=sec, rc=rc)
  if (chkerr(rc,__LINE__,u_FILE_u)) return
  mcdate = year*10000 + mon*100 + day

  nstreams = size(sdat)
  ! Advance the streams
  do ns = 1,nstreams
    write(stream_name,fmt='(a,i2.2)') 'stream_', ns
    call shr_strdata_advance(sdat(ns), ymd=mcdate, tod=sec, logunit=logunit, istr=trim(stream_name),rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    nflds = size(sdat(ns)%pstrm(1)%fldlist_model)
    do nf = 1,nflds
      fldname = trim(sdat(ns)%pstrm(1)%fldlist_model(nf))

      if (fldname == 'lrunoff') then
        ! Get pointer for stream data that is time and spatially interpolated to model time and grid
        call dshr_fldbun_getFldPtr(sdat(ns)%pstrm(1)%fldbun_model, trim(fldname), dataPtr1d, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return

        n = 0
        do j = jsc,jec
          do i = isc,iec
            n = n + 1
            ice_ocean_boundary%lrunoff(i,j)  = dataPtr1d(n)
          enddo
        enddo
      endif

      if (dbug > 1) then
        call dshr_fldbun_Field_diagnose(sdat(ns)%pstrm(1)%fldbun_model, trim(fldname), 'inline_run ', rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return
      endif
    enddo !nf
  enddo !ns

end subroutine mom_inline_run
end module mom_inline_mod
