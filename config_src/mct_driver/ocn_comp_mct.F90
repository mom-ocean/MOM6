!> This is the main driver for MOM6 in CIME
module ocn_comp_mct

! This file is part of MOM6. See LICENSE.md for the license.

! mct modules
use ESMF,                only: ESMF_clock, ESMF_time, ESMF_timeInterval
use ESMF,                only: ESMF_ClockGet, ESMF_TimeGet, ESMF_TimeIntervalGet
use seq_cdata_mod,       only: seq_cdata, seq_cdata_setptrs
use seq_flds_mod,        only: seq_flds_x2o_fields, seq_flds_o2x_fields
use mct_mod,             only: mct_gsMap, mct_gsmap_init, mct_gsMap_lsize, &
                               mct_gsmap_orderedpoints
use mct_mod,             only: mct_aVect, mct_aVect_init, mct_aVect_zero, &
                               mct_aVect_nRattr
use mct_mod,             only: mct_gGrid, mct_gGrid_init, mct_gGrid_importRAttr, &
                               mct_gGrid_importIAttr
use seq_infodata_mod,    only: seq_infodata_type, seq_infodata_GetData, &
                               seq_infodata_start_type_start, seq_infodata_start_type_cont,   &
                               seq_infodata_start_type_brnch, seq_infodata_PutData
use seq_comm_mct,        only: seq_comm_name, seq_comm_inst, seq_comm_suffix
use seq_timemgr_mod,     only: seq_timemgr_EClockGetData, seq_timemgr_RestartAlarmIsOn
use perf_mod,            only: t_startf, t_stopf
use shr_file_mod,        only: shr_file_getUnit, shr_file_freeUnit, shr_file_setIO, &
                               shr_file_getLogUnit, shr_file_getLogLevel, &
                               shr_file_setLogUnit, shr_file_setLogLevel

! MOM6 modules
use MOM,                  only: extract_surface_state
use MOM_variables,        only: surface
use MOM_domains,          only: MOM_infra_init
use MOM_restart,          only: save_restart
use MOM_ice_shelf,        only: ice_shelf_save_restart
use MOM_domains,          only: num_pes, root_pe, pe_here
use MOM_grid,             only: ocean_grid_type, get_global_grid_size
use MOM_error_handler,    only: MOM_error, FATAL, is_root_pe, WARNING
use MOM_time_manager,     only: time_type, set_date, set_time, set_calendar_type, NOLEAP
use MOM_time_manager,     only: operator(+), operator(-), operator(*), operator(/)
use MOM_time_manager,     only: operator(==), operator(/=), operator(>), get_time
use MOM_file_parser,      only: get_param, log_version, param_file_type, close_param_file
use MOM_get_input,        only: Get_MOM_Input, directories
use MOM_EOS,              only: gsw_sp_from_sr, gsw_pt_from_ct
use MOM_constants,        only: CELSIUS_KELVIN_OFFSET
use MOM_domains,          only: AGRID, BGRID_NE, CGRID_NE, pass_vector
use mpp_domains_mod,      only: mpp_get_compute_domain

! Previously inlined - now in separate modules
use MOM_ocean_model_mct,     only: ocean_public_type, ocean_state_type
use MOM_ocean_model_mct,     only: ocean_model_init , update_ocean_model, ocean_model_end
use MOM_ocean_model_mct,     only: convert_state_to_ocean_type
use MOM_surface_forcing_mct, only: surface_forcing_CS, forcing_save_restart, ice_ocean_boundary_type
use ocn_cap_methods,      only: ocn_import, ocn_export

! FMS modules
use time_interp_external_mod, only : time_interp_external

! MCT indices structure and import and export routines that access mom data
use ocn_cpl_indices,   only : cpl_indices_type, cpl_indices_init

! GFDL coupler modules
use coupler_types_mod,   only : coupler_type_spawn
use coupler_types_mod,   only : coupler_type_initialized, coupler_type_copy_data

! By default make data private
implicit none; private

#include <MOM_memory.h>

! Public member functions
public :: ocn_init_mct
public :: ocn_run_mct
public :: ocn_final_mct

! Private member functions
private :: ocn_SetGSMap_mct
private :: ocn_domain_mct
private :: get_runtype
private :: ocean_model_init_sfc

! Flag for debugging
logical, parameter :: debug=.true.

!> Control structure for this module
type MCT_MOM_Data
  type(ocean_state_type),  pointer :: ocn_state => NULL()  !< The private state of ocean
  type(ocean_public_type), pointer :: ocn_public => NULL() !< The public state of ocean
  type(ocean_grid_type),   pointer :: grid => NULL()       !< The grid structure
  type(seq_infodata_type), pointer :: infodata             !< The input info type
  type(cpl_indices_type)           :: ind                  !< Variable IDs
  logical                          :: sw_decomp            !< Controls whether shortwave is decomposed into 4 components
  real                             :: c1, c2, c3, c4       !< Coeffs. used in the shortwave decomposition  i/o
  integer                          :: stdout               !< standard output unit. (by default, points to ocn.log.* )
  character(len=384)               :: pointer_filename     !< Name of the ascii file that contains the path
                                                           !! and filename of the latest restart file.
end type MCT_MOM_Data

type(MCT_MOM_Data)            :: glb !< global structure
type(ice_ocean_boundary_type) :: ice_ocean_boundary

!=======================================================================
contains
!=======================================================================

!> This subroutine initializes MOM6.
subroutine ocn_init_mct( EClock, cdata_o, x2o_o, o2x_o, NLFilename )
  type(ESMF_Clock),             intent(inout) :: EClock      !< Time and time step ? \todo Why must this
                                                             !! be intent(inout)?
  type(seq_cdata)             , intent(inout) :: cdata_o     !< Input parameters
  type(mct_aVect)             , intent(inout) :: x2o_o       !< Fluxes from coupler to ocean, computed by ocean
  type(mct_aVect)             , intent(inout) :: o2x_o       !< Fluxes from ocean to coupler, computed by ocean
  character(len=*), optional  , intent(in)    :: NLFilename  !< Namelist filename

  !  local variable
  type(time_type)         :: time0                !< Start time of coupled model's calendar.
  type(time_type)         :: time_start           !< The time at which to initialize the ocean model
  type(ESMF_time)         :: time_var             !< ESMF_time variable to query time
  type(ESMF_time)         :: time_in_ESMF         !< Initial time for ocean
  type(ESMF_timeInterval) :: ocn_cpl_interval     !< Ocean coupling interval
  integer                 :: ncouple_per_day
  integer                 :: year, month, day, hour, minute, seconds, seconds_n, seconds_d, rc
  character(len=240)      :: runid                !< Run ID
  character(len=32)       :: runtype              !< Run type
  character(len=240)      :: restartfile          !< Path/Name of restart file
  integer                 :: nu                   !< i/o unit to read pointer file
  character(len=240)      :: restart_pointer_file !< File name for restart pointer file
  character(len=240)      :: restartpath          !< Path of the restart file
  integer                 :: mpicom_ocn           !< MPI ocn communicator
  integer                 :: npes, pe0            !< # of processors and current processor
  integer                 :: i, errorCode
  integer                 :: lsize, nsend, nrecv
  logical                 :: ldiag_cpl = .false.
  integer                 :: isc, iec, jsc, jec, ni, nj !< Indices for the start and end of the domain
                                                        !! in the x and y dir., respectively.
  ! runtime params
  type(param_file_type) :: param_file           !< A structure to parse for run-time parameters
  type(directories)     :: dirs_tmp             !< A structure containing several relevant directory paths
  character(len=40)     :: mdl = "ocn_comp_mct" !< This module's name.

  ! mct variables (these are local for now)
  integer                   :: MOM_MCT_ID
  type(mct_gsMap), pointer  :: MOM_MCT_gsMap => NULL() !< 2d, points to cdata
  type(mct_gGrid), pointer  :: MOM_MCT_dom => NULL()   !< 2d, points to cdata
  type(mct_gsMap)           :: MOM_MCT_gsMap3d         !< for 3d streams, local
  type(mct_gGrid)           :: MOM_MCT_dom3d           !< for 3d streams, local

  ! time management
  integer                   :: ocn_cpl_dt   !< one ocn coupling interval in seconds. (to be received from cesm)
  real (kind=8)             :: mom_cpl_dt   !< one ocn coupling interval in seconds. (internal)
  real (kind=8), parameter  ::        &
      seconds_in_minute =    60.0d0, &
      seconds_in_hour   =  3600.0d0, &
      seconds_in_day    = 86400.0d0, &
      minutes_in_hour   =    60.0d0

  character(len=99) :: ocn_modelio_name !< ocn model input namelist filename
  integer           :: shrlogunit       !< original log file unit
  integer           :: shrloglev        !< original log level

  integer(kind=4)     :: inst_index     !< instance control vars (these are local for now)
  character(len=16)   :: inst_name
  character(len=16)   :: inst_suffix

  ! TODO: Change the following vars with the corresponding MOM6 vars
  integer :: km=1                   !< Number of vertical levels
  !logical :: lsend_precip_fact      !< If T,send precip_fact to cpl for use in fw balance
                                    !! (partially-coupled option)
  character(len=128) :: err_msg     !< Error message

  ! set the cdata pointers:
  call seq_cdata_setptrs(cdata_o, id=MOM_MCT_ID, mpicom=mpicom_ocn, &
                         gsMap=MOM_MCT_gsMap, dom=MOM_MCT_dom, infodata=glb%infodata)

  ! Determine attribute vector indices
  call cpl_indices_init(glb%ind)

  call seq_infodata_GetData( glb%infodata, case_name=runid )

  ! instance control
  inst_name   = seq_comm_name(MOM_MCT_ID)
  inst_index  = seq_comm_inst(MOM_MCT_ID)
  inst_suffix = seq_comm_suffix(MOM_MCT_ID)

  call t_startf('MOM_init')

  ! Initialize MOM6 comm
  call MOM_infra_init(mpicom_ocn)

  ! initialize ocn log file
  if (is_root_pe()) then

    ! get original log file properties
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)

    glb%stdout = shr_file_getUnit() ! get an unused unit number

    ! open the ocn_modelio.nml file and then open a log file associated with stdout
    ocn_modelio_name = 'ocn_modelio.nml' // trim(inst_suffix)
    call shr_file_setIO(ocn_modelio_name,glb%stdout)

    !  set the shr log io unit number
    call shr_file_setLogUnit(glb%stdout)
  end if

  call set_calendar_type(NOLEAP)  !TODO: confirm this

  ! Get start time
  call ESMF_ClockGet(EClock, StartTime=time_var, rc=rc)
  call ESMF_TimeGet(time_var, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
  time0 = set_date(year, month, day, hour, minute, seconds, err_msg=err_msg)

  ! Get current time
  call ESMF_ClockGet(EClock, currTime=time_var, rc=rc)
  call ESMF_TimeGet(time_var, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
  time_start = set_date(year, month, day, hour, minute, seconds, err_msg=err_msg)

  ! Debugging clocks
  if (debug .and. is_root_pe()) then
    write(glb%stdout,*) 'ocn_init_mct, current time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds

    call ESMF_ClockGet(EClock, StartTime=time_var, rc=rc)
    call ESMF_TimeGet(time_var, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(glb%stdout,*) 'ocn_init_mct, start time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds

    call ESMF_ClockGet(EClock, StopTime=time_var, rc=rc)
    call ESMF_TimeGet(time_var, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(glb%stdout,*) 'ocn_init_mct, stop time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds

    call ESMF_ClockGet(EClock, PrevTime=time_var, rc=rc)
    call ESMF_TimeGet(time_var, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(glb%stdout,*) 'ocn_init_mct, previous time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds

    call ESMF_ClockGet(EClock, TimeStep=ocn_cpl_interval, rc=rc)
    call ESMF_TimeIntervalGet(ocn_cpl_interval, yy=year, mm=month, d=day, s=seconds, sn=seconds_n, sd=seconds_d, rc=rc)
    write(glb%stdout,*) 'ocn_init_mct, time step: y,m,d-',year,month,day,'s,sn,sd=',seconds,seconds_n,seconds_d
  endif

  npes = num_pes()
  pe0 = root_pe()

  allocate(glb%ocn_public)
  glb%ocn_public%is_ocean_PE = .true.

  allocate(glb%ocn_public%pelist(npes))
  glb%ocn_public%pelist(:) = (/(i,i=pe0,pe0+npes)/)
  ! \todo Set other bits of glb$ocn_public

  ! This include declares and sets the variable "version".
  ! read useful runtime params
  call get_MOM_Input(param_file, dirs_tmp, check_params=.false.)
  !call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "POINTER_FILENAME", glb%pointer_filename, &
                 "Name of the ascii file that contains the path and filename of" // &
                 " the latest restart file.", default='rpointer.ocn')

  call get_param(param_file, mdl, "SW_DECOMP", glb%sw_decomp, &
                 "If True, read coeffs c1, c2, c3 and c4 and decompose" // &
                 "the net shortwave radiation (SW) into four components:\n" // &
                 "visible, direct shortwave  = c1 * SW \n" // &
                 "visible, diffuse shortwave = c2 * SW \n" // &
                 "near-IR, direct shortwave  = c3 * SW \n" // &
                 "near-IR, diffuse shortwave = c4 * SW", default=.true.)

  if (glb%sw_decomp) then
    call get_param(param_file, mdl, "SW_c1", glb%c1, &
                  "Coeff. used to convert net shortwave rad. into "//&
                  "visible, direct shortwave.", units="nondim", default=0.285)

    call get_param(param_file, mdl, "SW_c2", glb%c2, &
                  "Coeff. used to convert net shortwave rad. into "//&
                  "visible, diffuse shortwave.", units="nondim", default=0.285)

    call get_param(param_file, mdl, "SW_c3", glb%c3, &
                  "Coeff. used to convert net shortwave rad. into "//&
                  "near-IR, direct shortwave.", units="nondim", default=0.215)

    call get_param(param_file, mdl, "SW_c4", glb%c4, &
                  "Coeff. used to convert net shortwave rad. into "//&
                  "near-IR, diffuse shortwave.", units="nondim", default=0.215)
  else
    glb%c1 = 0.0; glb%c2 = 0.0; glb%c3 = 0.0; glb%c4 = 0.0
  endif

  ! Close param file before it gets opened by ocean_model_init again.
  call close_param_file(param_file)

  ! Initialize the MOM6 model
  runtype = get_runtype()
  if (runtype == "initial") then
    ! startup (new run) - 'n' is needed below since we don't specify input_filename in input.nml
    call ocean_model_init(glb%ocn_public, glb%ocn_state, time0, time_start, input_restart_file = 'n')
  else  ! hybrid or branch or continuos runs
    ! get output path root
    call seq_infodata_GetData( glb%infodata, outPathRoot=restartpath )
    ! read name of restart file in the pointer file
    nu = shr_file_getUnit()
    restart_pointer_file = trim(glb%pointer_filename)
    if (is_root_pe()) write(glb%stdout,*) 'Reading ocn pointer file: ',restart_pointer_file
    open(nu, file=restart_pointer_file, form='formatted', status='unknown')
    read(nu,'(a)') restartfile
    close(nu)
    !restartfile = trim(restartpath) // trim(restartfile)
    if (is_root_pe()) then
      write(glb%stdout,*) 'Reading restart file: ',trim(restartfile)
    end if
    call shr_file_freeUnit(nu)
    call ocean_model_init(glb%ocn_public, glb%ocn_state, time0, time_start, input_restart_file=trim(restartfile))
  endif
  if (is_root_pe()) then
    write(glb%stdout,'(/12x,a/)') '======== COMPLETED MOM INITIALIZATION ========'
  end if

  ! Initialize ocn_state%sfc_state out of sight
  call ocean_model_init_sfc(glb%ocn_state, glb%ocn_public)

  ! Store pointers to components inside MOM
  glb%grid => glb%ocn_state%grid

  ! Allocate IOB data type (needs to be called after glb%grid is set)
  !write(6,*)'DEBUG: isc,iec,jsc,jec= ',glb%grid%isc, glb%grid%iec, glb%grid%jsc, glb%grid%jec
  call IOB_allocate(ice_ocean_boundary, glb%grid%isc, glb%grid%iec, glb%grid%jsc, glb%grid%jec)

  call t_stopf('MOM_init')

  ! Initialize MCT attribute vectors and indices
  call t_startf('MOM_mct_init')

  if (debug .and. root_pe().eq.pe_here()) print *, "calling ocn_SetGSMap_mct"

  ! Set mct global seg maps:

  call ocn_SetGSMap_mct(mpicom_ocn, MOM_MCT_ID, MOM_MCT_GSMap, MOM_MCT_GSMap3d)
  lsize = mct_gsMap_lsize(MOM_MCT_gsmap, mpicom_ocn)

  ! Initialize mct ocn domain (needs ocn initialization info)

  if (debug .and. root_pe().eq.pe_here()) print *, "calling ocn_domain_mct"
  call ocn_domain_mct(lsize, MOM_MCT_gsmap, MOM_MCT_dom)
  !call ocn_domain_mct(lsize*km, MOM_MCT_gsmap3d, MOM_MCT_dom3d) !TODO: this is not used

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
  ncouple_per_day = seconds_in_day / ocn_cpl_dt
  mom_cpl_dt = seconds_in_day / ncouple_per_day
  if (mom_cpl_dt /= ocn_cpl_dt) then
    write(glb%stdout,*) 'ERROR mom_cpl_dt and ocn_cpl_dt must be identical'
    call exit(0)
  end if

  ! send initial state to driver

  !TODO:
  ! if ( lsend_precip_fact )  then
  !    call seq_infodata_PutData( infodata, precip_fact=precip_fact)
  ! end if

  if (debug .and. root_pe().eq.pe_here()) print *, "calling ocn_export"
  call ocn_export(glb%ind, glb%ocn_public, glb%grid, o2x_o%rattr, mom_cpl_dt, ncouple_per_day)

  call t_stopf('MOM_mct_init')

  ! Size of global domain
  call get_global_grid_size(glb%grid, ni, nj)

  if (debug .and. root_pe().eq.pe_here()) print *, "calling seq_infodata_putdata"

  call seq_infodata_PutData( glb%infodata, &
       ocn_nx = ni , ocn_ny = nj)
  call seq_infodata_PutData( glb%infodata, &
       ocn_prognostic=.true., ocnrof_prognostic=.true.)

  if (debug .and. root_pe().eq.pe_here()) print *, "leaving ocean_init_mct"

  ! Reset shr logging to original values
  if (is_root_pe()) then
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
  end if

end subroutine ocn_init_mct

!=======================================================================

!> Step forward ocean model for coupling interval
subroutine ocn_run_mct( EClock, cdata_o, x2o_o, o2x_o)
  type(ESMF_Clock), intent(inout) :: EClock  !< Time and time step ? \todo Why must this be intent(inout)?
  type(seq_cdata),  intent(inout) :: cdata_o !< Input parameters
  type(mct_aVect),  intent(inout) :: x2o_o   !< Fluxes from coupler to ocean, computed by ocean
  type(mct_aVect),  intent(inout) :: o2x_o   !< Fluxes from ocean to coupler, computed by ocean
  ! Local variables
  type(ESMF_time) :: time_var                 !< ESMF_time variable to query time
  type(ESMF_timeInterval) :: ocn_cpl_interval !< The length of one ocean coupling interval
  integer :: year, month, day, hour, minute, seconds, seconds_n, seconds_d, rc
  logical :: write_restart_at_eod      !< Controls if restart files must be written
  logical :: debug=.false.
  type(time_type) :: time_start        !< Start of coupled time interval to pass to MOM6
  type(time_type) :: coupling_timestep !< Coupled time interval to pass to MOM6
  character(len=128) :: err_msg        !< Error message
  character(len=32)  :: timestamp      !< Name of intermediate restart file
  character(len=384) :: restartname    !< The restart file name (no dir)
  character(len=384) :: restart_pointer_file !< File name for restart pointer file
  character(len=384) :: runid                !< Run ID
  character(len=32)  :: runtype              !< Run type
  integer            :: nu                   !< i/o unit to write pointer file
  integer            :: shrlogunit ! original log file unit
  integer            :: shrloglev  ! original log level
  logical, save      :: firstCall = .true.
  real (kind=8), parameter  ::  seconds_in_day = 86400.0 !< number of seconds in one day
  integer                   :: ocn_cpl_dt   !< one ocn coupling interval in seconds. (to be received from cesm)
  real (kind=8)             :: mom_cpl_dt   !< one ocn coupling interval in seconds. (internal)
  integer                   :: ncouple_per_day !< number of ocean coupled call in one day (non-dim)

  ! reset shr logging to ocn log file:
  if (is_root_pe()) then
    call shr_file_getLogUnit(shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit(glb%stdout)
  endif

  ! Query the beginning time of the current coupling interval
  call ESMF_ClockGet(EClock, PrevTime=time_var, rc=rc)
  call ESMF_TimeGet(time_var, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
  time_start = set_date(year, month, day, hour, minute, seconds, err_msg=err_msg)

  ! Query the coupling interval duration
  call ESMF_ClockGet(EClock, TimeStep=ocn_cpl_interval, rc=rc)
  call ESMF_TimeIntervalGet(ocn_cpl_interval, yy=year, mm=month, d=day, s=seconds, sn=seconds_n, sd=seconds_d, rc=rc)
  coupling_timestep = set_time(seconds, days=day, err_msg=err_msg)

  call seq_timemgr_EClockGetData(EClock, dtime=ocn_cpl_dt)
  ncouple_per_day = seconds_in_day / ocn_cpl_dt
  mom_cpl_dt = seconds_in_day / ncouple_per_day

  ! The following if-block is to correct monthly mean outputs:
  ! With this change, MOM6 starts at the same date as the other components, and runs for the same
  ! duration as other components, unlike POP, which would have one missing interval due to ocean
  ! lag. MOM6 accounts for this lag by doubling the duration of the first coupling interval.
  if (firstCall) then

    runtype = get_runtype()
    if (runtype /= "continue" .and. runtype /= "branch") then

      if (debug .and. is_root_pe()) then
        write(glb%stdout,*) 'doubling first interval duration!'
      endif

      ! shift back the start time by one coupling interval (to align the start time with other components)
      time_start = time_start-coupling_timestep
      ! double the first coupling interval (to account for the missing coupling interval to due to lag)
      coupling_timestep = coupling_timestep*2
    end if

    firstCall = .false.
  end if

  ! Debugging clocks
  if (debug .and. is_root_pe()) then
    call ESMF_ClockGet(EClock, CurrTime=time_var, rc=rc)
    call ESMF_TimeGet(time_var, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(glb%stdout,*) 'ocn_run_mct, current time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, StartTime=time_var, rc=rc)
    call ESMF_TimeGet(time_var, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(glb%stdout,*) 'ocn_run_mct, start time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, StopTime=time_var, rc=rc)
    call ESMF_TimeGet(time_var, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(glb%stdout,*) 'ocn_run_mct, stop time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, PrevTime=time_var, rc=rc)
    call ESMF_TimeGet(time_var, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    write(glb%stdout,*) 'ocn_run_mct, previous time: y,m,d-',year,month,day,'h,m,s=',hour,minute,seconds
    call ESMF_ClockGet(EClock, TimeStep=ocn_cpl_interval, rc=rc)
    call ESMF_TimeIntervalGet(ocn_cpl_interval, yy=year, mm=month, d=day, s=seconds, sn=seconds_n, sd=seconds_d, rc=rc)
    write(glb%stdout,*) 'ocn_init_mct, time step: y,m,d-',year,month,day,'s,sn,sd=',seconds,seconds_n,seconds_d
  endif

  ! set the cdata pointers:
  ! \todo this was done in _init_, is it needed again. Does this infodata need to be in glb%?
  ! GMM, check if  this is needed!
  call seq_cdata_setptrs(cdata_o, infodata=glb%infodata)

  ! Translate import fields to ice_ocean_boundary
  !TODO: make this an input variable
  !glb%sw_decomp = .false.
  !END TODO:
  if (glb%sw_decomp) then
    call ocn_import(x2o_o%rattr, glb%ind,  glb%grid, Ice_ocean_boundary, glb%ocn_public, glb%stdout, Eclock, &
          c1=glb%c1, c2=glb%c2, c3=glb%c3, c4=glb%c4)
  else
    call ocn_import(x2o_o%rattr, glb%ind,  glb%grid, Ice_ocean_boundary, glb%ocn_public, glb%stdout, Eclock )
  end if

  ! Update internal ocean
  call update_ocean_model(ice_ocean_boundary, glb%ocn_state, glb%ocn_public, time_start, coupling_timestep)

  ! Return export state to driver
  call ocn_export(glb%ind, glb%ocn_public, glb%grid, o2x_o%rattr, mom_cpl_dt, ncouple_per_day)

  !--- write out intermediate restart file when needed.
  ! Check alarms for flag to write restart at end of day
  write_restart_at_eod = seq_timemgr_RestartAlarmIsOn(EClock)
  if (debug .and. is_root_pe()) write(glb%stdout,*) 'ocn_run_mct, write_restart_at_eod=', write_restart_at_eod

  if (write_restart_at_eod) then
    ! case name
    call seq_infodata_GetData( glb%infodata, case_name=runid )
    ! add time stamp to the restart filename
    call ESMF_ClockGet(EClock, CurrTime=time_var, rc=rc)
    call ESMF_TimeGet(time_var, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    seconds = seconds + hour*3600 + minute*60
    write(restartname,'(A,".mom6.r.",I4.4,"-",I2.2,"-",I2.2,"-",I5.5)') trim(runid), year, month, day, seconds

    call save_restart(glb%ocn_state%dirs%restart_output_dir, glb%ocn_state%Time, glb%grid, &
                      glb%ocn_state%restart_CSp, .false., filename=restartname, GV=glb%ocn_state%GV)

    ! write name of restart file in the rpointer file
    nu = shr_file_getUnit()
    if (is_root_pe()) then
      restart_pointer_file = trim(glb%pointer_filename)
      open(nu, file=restart_pointer_file, form='formatted', status='unknown')
      write(nu,'(a)') trim(restartname) //'.nc'
      close(nu)
      write(glb%stdout,*) 'ocn restart pointer file written: ',trim(restartname)
    endif
    call shr_file_freeUnit(nu)

    ! Is this needed?
    call forcing_save_restart(glb%ocn_state%forcing_CSp, glb%grid, glb%ocn_state%Time, &
                              glb%ocn_state%dirs%restart_output_dir, .true.)

    ! Once we start using the ice shelf module, the following will be needed
    if (glb%ocn_state%use_ice_shelf) then
      call ice_shelf_save_restart(glb%ocn_state%Ice_shelf_CSp, glb%ocn_state%Time, &
                                  glb%ocn_state%dirs%restart_output_dir, .true.)
    endif

  endif

  ! reset shr logging to original values
  if (is_root_pe()) then
    call shr_file_setLogUnit(shrlogunit)
    call shr_file_setLogLevel(shrloglev)
  endif

end subroutine ocn_run_mct

!=======================================================================

!> Finalizes MOM6
!!
!! \todo This needs to be done here.
subroutine ocn_final_mct( EClock, cdata_o, x2o_o, o2x_o)
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o  !< Fluxes from coupler to ocean, computed by ocean
    type(mct_aVect)             , intent(inout) :: o2x_o  !< Fluxes from ocean to coupler, computed by ocean

    call ocean_model_end(glb%ocn_public, glb%ocn_state, glb%ocn_state%Time)

end subroutine ocn_final_mct

!=======================================================================

!> Sets mct global segment maps for the MOM decomposition.
!!
!! \todo Find out if we should only provide indirect indexing for ocean points and not land.
subroutine ocn_SetGSMap_mct(mpicom_ocn, MOM_MCT_ID, gsMap_ocn, gsMap3d_ocn)
  integer,         intent(in)    :: mpicom_ocn  !< MPI communicator
  integer,         intent(in)    :: MOM_MCT_ID  !< MCT component ID
  type(mct_gsMap), intent(inout) :: gsMap_ocn   !< MCT global segment map for 2d data
  type(mct_gsMap), intent(inout) :: gsMap3d_ocn !< MCT global segment map for 3d data

  ! Local variables
  integer                        :: lsize          !< Local size of indirect indexing array
  integer                        :: i, j, k        !< Local indices
  integer                        :: ni, nj         !< Declared sizes of h-point arrays
  integer                        :: ig, jg         !< Global indices
  type(ocean_grid_type), pointer :: grid => NULL() !< A pointer to a grid structure
  integer, allocatable           :: gindex(:)      !< Indirect indices

  grid => glb%grid ! for convenience
  if (.not. associated(grid)) call MOM_error(FATAL, 'ocn_comp_mct.F90, ocn_SetGSMap_mct():' // &
       'grid is not associated!')

  ! Size of computational domain
  lsize = ( grid%iec - grid%isc + 1 ) * ( grid%jec - grid%jsc + 1 )

  ! Size of global domain
  call get_global_grid_size(grid, ni, nj)

  ! Create indirect indices for the computational domain
  allocate(gindex(lsize))

  ! Set indirect indices in gindex
  k = 0
  do j = grid%jsc, grid%jec
    jg = j + grid%jdg_offset ! TODO: check this calculation
    do i = grid%isc, grid%iec
      ig = i + grid%idg_offset ! TODO: check this calculation
      k = k + 1 ! Increment position within gindex
      gindex(k) = ni * (jg - 1) + ig
    enddo
  enddo

  ! Tell MCT how to indirectly index into the 2d buffer
  call mct_gsMap_init(gsMap_ocn, gindex, mpicom_ocn, MOM_MCT_ID, lsize, ni * nj)

  deallocate(gindex)

end subroutine ocn_SetGSMap_mct

!=======================================================================

!> Sets MCT global segment maps for the MOM6 decomposition
subroutine ocn_domain_mct( lsize, gsMap_ocn, dom_ocn)
  integer        , intent(in)    :: lsize      !< Size of attr. vector
  type(mct_gsMap), intent(in)    :: gsMap_ocn  !< MCT global segment map for 2d data
  type(mct_ggrid), intent(inout) :: dom_ocn    !< WHAT IS THIS?

  ! Local Variables
  integer, parameter              :: SHR_REAL_R8 = selected_real_kind(12)
  integer, pointer                :: idata(:)
  integer                         :: i,j,k
  real(kind=SHR_REAL_R8), pointer :: data(:)
  real(kind=SHR_REAL_R8)          :: L2_to_rad2
  type(ocean_grid_type), pointer :: grid => NULL() ! A pointer to a grid structure

  grid => glb%grid ! for convenience

  ! set coords to lat and lon, and areas to rad^2
  call mct_gGrid_init(GGrid=dom_ocn, CoordChars='lat:lon:hgt', OtherChars='area:aream:mask:frac', lsize=lsize )

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
  L2_to_rad2 = grid%US%L_to_m**2 / grid%Rad_Earth**2
  do j = grid%jsc, grid%jec
    do i = grid%isc, grid%iec
      k = k + 1 ! Increment position within gindex
      data(k) = grid%AreaT(i,j) * L2_to_rad2
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

end subroutine ocn_domain_mct

!=======================================================================

!> Returns the CESM run type
character(32) function get_runtype()
  character(len=32)   :: starttype         !< infodata start type

  call seq_infodata_GetData( glb%infodata, start_type=starttype)

    if (   trim(starttype) == trim(seq_infodata_start_type_start)) then
     get_runtype = "initial"
  else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
     get_runtype = "continue"
  else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
     get_runtype = "branch"
  else
     write(glb%stdout,*) 'ocn_comp_mct ERROR: unknown starttype'
     call exit(0)
  end if
  return

end function

!=======================================================================

!> It has to be separate from the ocean_initialization call because the coupler
!! module allocates the space for some of these variables.
subroutine ocean_model_init_sfc(OS, Ocean_sfc)
  type(ocean_state_type),  pointer       :: OS
  type(ocean_public_type), intent(inout) :: Ocean_sfc

  integer :: is, ie, js, je

  is = OS%grid%isc ; ie = OS%grid%iec ; js = OS%grid%jsc ; je = OS%grid%jec
  call coupler_type_spawn(Ocean_sfc%fields, OS%sfc_state%tr_fields, &
                          (/is,is,ie,ie/), (/js,js,je,je/), as_needed=.true.)

  call extract_surface_state(OS%MOM_CSp, OS%sfc_state)

  call convert_state_to_ocean_type(OS%sfc_state, Ocean_sfc, OS%grid, OS%US)

end subroutine ocean_model_init_sfc

!=======================================================================

!> \namespace ocn_comp_mct
!!
!! \section section_ocn_import Fluxes imported from the coupler (MCT) to MOM6
!! The following summarizes the mismatches between MCT and MOM6 in terms
!! of ice ocean fluxes.
!!
!! Redundancies:
!! x2o_Faxa_prec = x2o_Faxa_rain + x2o_Faxa_snow
!!
!! Variables whose units and sign  **could not** be verified so far:
!! x2o_Foxx_rofl
!! x2o_Foxx_rof
!!
!! Variables in MOM6 fluxes that are **NOT** filled by the coupler:
!! ustar_berg, frictional velocity beneath icebergs [m s-1]
!! area_berg, area covered by icebergs(m2/m2)
!! mass_berg, mass of icebergs(kg/m2)
!! runoff_hflx, heat content of liquid runoff (W/m2)
!! calving_hflx, heat content of frozen runoff (W/m2)
!! mi, mass of ice (kg/m2)
!!
!! Variables in the coupler that are **NOT** used in MOM6 (i.e., no corresponding field in fluxes):
!! x2o_Si_ifrac, fractional ice wrt ocean
!! x2o_So_duu10n, 10m wind speed squared (m^2/s^2)
!! x2o_Sa_co2prog, bottom atm level prognostic CO2
!! x2o_Sa_co2diag, bottom atm level diagnostic CO2
!!
!! \TODO Langmuir related fields:
!! surface Stokes drift, x-comp. (x2o_Sw_ustokes)
!! surface Stokes drift, y-comp. (x2o_Sw_vstokes)
!! wave model langmuir multiplier (x2o_Sw_lamult)
!!
!! \TODO Biogeochemistry:
!! x2o_Fioi_bcpho, Black Carbon hydrophobic release from sea ice component
!! x2o_Fioi_bcphi, Black Carbon hydrophilic release from sea ice component
!! x2o_Fioi_flxdst, Dust release from sea ice component
!! x2o_Faxa_bcphidry, Black Carbon hydrophilic dry deposition
!! x2o_Faxa_bcphodry, Black Carbon hydrophobic dry deposition
!! x2o_Faxa_bcphiwet, Black Carbon hydrophobic wet deposition
!! x2o_Faxa_ocphidry, Organic Carbon hydrophilic dry deposition
!! x2o_Faxa_ocphodry, Organic Carbon hydrophobic dry deposition
!! x2o_Faxa_ocphiwet, Organic Carbon hydrophilic dry deposition
!! x2o_Faxa_dstwet, Sizes 1 to 4 dust - wet deposition
!! x2o_Faxa_dstdry, Sizes 1 to 4 dust - dry deposition
!!
!! \section section_ocn_export Fluxes exported from MOM6 to the coupler (MCT)
!!
!! Variables that are currently being exported:
!!
!! Surface temperature (Kelvin)
!! Surface salinity (psu)
!! Surface eastward velocity [m s-1]
!! Surface northward velocity [m s-1]
!! Zonal slope in the sea surface height
!! Meridional slope in the sea surface height
!!
!! \TODO Variables that **are not** currently being exported:
!!
!! Boundary layer depth
!! CO2
!! DMS

!> Allocates ice-ocean boundary type containers and sets to 0.
subroutine IOB_allocate(IOB, isc, iec, jsc, jec)
  type(ice_ocean_boundary_type), intent(inout)    :: IOB    !< An ice-ocean boundary type with fluxes to drive
  integer, intent(in) :: isc, iec, jsc, jec                 !< The ocean's local grid size

  allocate ( IOB% rofl_flux (isc:iec,jsc:jec),       &
             IOB% rofi_flux (isc:iec,jsc:jec),       &
             IOB% u_flux (isc:iec,jsc:jec),          &
             IOB% v_flux (isc:iec,jsc:jec),          &
             IOB% t_flux (isc:iec,jsc:jec),          &
             IOB% seaice_melt_heat (isc:iec,jsc:jec),&
             IOB% seaice_melt (isc:iec,jsc:jec),     &
             IOB% q_flux (isc:iec,jsc:jec),          &
             IOB% salt_flux (isc:iec,jsc:jec),       &
             IOB% lw_flux (isc:iec,jsc:jec),         &
             IOB% sw_flux_vis_dir (isc:iec,jsc:jec), &
             IOB% sw_flux_vis_dif (isc:iec,jsc:jec), &
             IOB% sw_flux_nir_dir (isc:iec,jsc:jec), &
             IOB% sw_flux_nir_dif (isc:iec,jsc:jec), &
             IOB% lprec (isc:iec,jsc:jec),           &
             IOB% fprec (isc:iec,jsc:jec),           &
             IOB% ustar_berg (isc:iec,jsc:jec),      &
             IOB% area_berg (isc:iec,jsc:jec),       &
             IOB% mass_berg (isc:iec,jsc:jec),       &
             IOB% calving (isc:iec,jsc:jec),         &
             IOB% runoff_hflx (isc:iec,jsc:jec),     &
             IOB% calving_hflx (isc:iec,jsc:jec),    &
             IOB% mi (isc:iec,jsc:jec),              &
             IOB% p (isc:iec,jsc:jec))

  IOB%rofl_flux        = 0.0
  IOB%rofi_flux        = 0.0
  IOB%u_flux           = 0.0
  IOB%v_flux           = 0.0
  IOB%t_flux           = 0.0
  IOB%seaice_melt_heat = 0.0
  IOB%seaice_melt      = 0.0
  IOB%q_flux           = 0.0
  IOB%salt_flux        = 0.0
  IOB%lw_flux          = 0.0
  IOB%sw_flux_vis_dir  = 0.0
  IOB%sw_flux_vis_dif  = 0.0
  IOB%sw_flux_nir_dir  = 0.0
  IOB%sw_flux_nir_dif  = 0.0
  IOB%lprec            = 0.0
  IOB%fprec            = 0.0
  IOB%ustar_berg       = 0.0
  IOB%area_berg        = 0.0
  IOB%mass_berg        = 0.0
  IOB%calving          = 0.0
  IOB%runoff_hflx      = 0.0
  IOB%calving_hflx     = 0.0
  IOB%mi               = 0.0
  IOB%p                = 0.0

end subroutine IOB_allocate

end module ocn_comp_mct
