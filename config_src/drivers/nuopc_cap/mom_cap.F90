!> This module contains a set of subroutines that are required by NUOPC.

module MOM_cap_mod

use constants_mod,            only: constants_init
use diag_manager_mod,         only: diag_manager_init, diag_manager_end
use field_manager_mod,        only: field_manager_init, field_manager_end
use fms_mod,                  only: fms_init, fms_end, open_namelist_file, check_nml_error
use fms_mod,                  only: close_file, file_exist, uppercase
use fms_io_mod,               only: fms_io_exit
use mpp_domains_mod,          only: domain2d, mpp_get_compute_domain, mpp_get_compute_domains
use mpp_domains_mod,          only: mpp_get_ntile_count, mpp_get_pelist, mpp_get_global_domain
use mpp_domains_mod,          only: mpp_get_domain_npes
use mpp_io_mod,               only: mpp_open, MPP_RDONLY, MPP_ASCII, MPP_OVERWR, MPP_APPEND, mpp_close, MPP_SINGLE
use mpp_mod,                  only: stdlog, stdout, mpp_root_pe, mpp_clock_id
use mpp_mod,                  only: mpp_clock_begin, mpp_clock_end, MPP_CLOCK_SYNC
use mpp_mod,                  only: MPP_CLOCK_DETAILED, CLOCK_COMPONENT, MAXPES
use time_manager_mod,         only: set_calendar_type, time_type, increment_date
use time_manager_mod,         only: set_time, set_date, get_time, get_date, month_name
use time_manager_mod,         only: GREGORIAN, JULIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
use time_manager_mod,         only: operator( <= ), operator( < ), operator( >= )
use time_manager_mod,         only: operator( + ),  operator( - ), operator( / )
use time_manager_mod,         only: operator( * ), operator( /= ), operator( > )
use time_manager_mod,         only: date_to_string
use time_manager_mod,         only: fms_get_calendar_type => get_calendar_type
use MOM_domains,              only: MOM_infra_init, num_pes, root_pe, pe_here
use MOM_file_parser,          only: get_param, log_version, param_file_type, close_param_file
use MOM_get_input,            only: get_MOM_input, directories
use MOM_domains,              only: pass_var
use MOM_error_handler,        only: MOM_error, FATAL, is_root_pe
use MOM_ocean_model_nuopc,    only: ice_ocean_boundary_type
use MOM_grid,                 only: ocean_grid_type, get_global_grid_size
use MOM_ocean_model_nuopc,    only: ocean_model_restart, ocean_public_type, ocean_state_type
use MOM_ocean_model_nuopc,    only: ocean_model_init_sfc
use MOM_ocean_model_nuopc,    only: ocean_model_init, update_ocean_model, ocean_model_end
use MOM_ocean_model_nuopc,    only: get_ocean_grid, get_eps_omesh
use MOM_cap_time,             only: AlarmInit
use MOM_cap_methods,          only: mom_import, mom_export, mom_set_geomtype, state_diagnose
use MOM_cap_methods,          only: ChkErr
#ifdef CESMCOUPLED
use shr_file_mod,             only: shr_file_setLogUnit, shr_file_getLogUnit
#endif
use time_utils_mod,           only: esmf2fms_time

use, intrinsic :: iso_fortran_env, only: output_unit

use ESMF,  only: ESMF_ClockAdvance, ESMF_ClockGet, ESMF_ClockPrint
use ESMF,  only: ESMF_ClockGetAlarm, ESMF_ClockGetNextTime, ESMF_ClockAdvance
use ESMF,  only: ESMF_ClockSet, ESMF_Clock, ESMF_GeomType_Flag, ESMF_LOGMSG_INFO
use ESMF,  only: ESMF_Grid, ESMF_GridCreate, ESMF_GridAddCoord
use ESMF,  only: ESMF_GridGetCoord, ESMF_GridAddItem, ESMF_GridGetItem
use ESMF,  only: ESMF_GridComp, ESMF_GridCompSetEntryPoint, ESMF_GridCompGet
use ESMF,  only: ESMF_LogFoundError, ESMF_LogWrite, ESMF_LogSetError
use ESMF,  only: ESMF_LOGERR_PASSTHRU, ESMF_KIND_R8, ESMF_RC_VAL_WRONG
use ESMF,  only: ESMF_GEOMTYPE_MESH, ESMF_GEOMTYPE_GRID, ESMF_SUCCESS
use ESMF,  only: ESMF_METHOD_INITIALIZE, ESMF_MethodRemove, ESMF_State
use ESMF,  only: ESMF_LOGMSG_INFO, ESMF_RC_ARG_BAD, ESMF_VM, ESMF_Time
use ESMF,  only: ESMF_TimeInterval, ESMF_MAXSTR, ESMF_VMGetCurrent
use ESMF,  only: ESMF_VMGet, ESMF_TimeGet, ESMF_TimeIntervalGet, ESMF_MeshGet
use ESMF,  only: ESMF_MethodExecute, ESMF_Mesh, ESMF_DeLayout, ESMF_Distgrid
use ESMF,  only: ESMF_DistGridConnection, ESMF_StateItem_Flag, ESMF_KIND_I4
use ESMF,  only: ESMF_KIND_I8, ESMF_FAILURE, ESMF_DistGridCreate, ESMF_MeshCreate
use ESMF,  only: ESMF_FILEFORMAT_ESMFMESH, ESMF_DELayoutCreate, ESMF_DistGridConnectionSet
use ESMF,  only: ESMF_DistGridGet, ESMF_STAGGERLOC_CORNER, ESMF_GRIDITEM_MASK
use ESMF,  only: ESMF_TYPEKIND_I4, ESMF_TYPEKIND_R8, ESMF_STAGGERLOC_CENTER
use ESMF,  only: ESMF_GRIDITEM_AREA, ESMF_Field, ESMF_ALARM, ESMF_VMLogMemInfo
use ESMF,  only: ESMF_AlarmIsRinging, ESMF_AlarmRingerOff, ESMF_StateRemove
use ESMF,  only: ESMF_FieldCreate, ESMF_LOGMSG_ERROR, ESMF_LOGMSG_WARNING
use ESMF,  only: ESMF_COORDSYS_SPH_DEG, ESMF_GridCreate, ESMF_INDEX_DELOCAL
use ESMF,  only: ESMF_MESHLOC_ELEMENT, ESMF_RC_VAL_OUTOFRANGE, ESMF_StateGet
use ESMF,  only: ESMF_TimePrint, ESMF_AlarmSet, ESMF_FieldGet, ESMF_Array
use ESMF,  only: ESMF_ArrayCreate
use ESMF,  only: ESMF_RC_FILE_OPEN, ESMF_RC_FILE_READ, ESMF_RC_FILE_WRITE
use ESMF,  only: ESMF_VMBroadcast
use ESMF,  only: ESMF_AlarmCreate, ESMF_ClockGetAlarmList, ESMF_AlarmList_Flag
use ESMF,  only: ESMF_AlarmGet, ESMF_AlarmIsCreated, ESMF_ALARMLIST_ALL, ESMF_AlarmIsEnabled
use ESMF,  only: ESMF_STATEITEM_NOTFOUND, ESMF_FieldWrite
use ESMF,  only: ESMF_END_ABORT, ESMF_Finalize
use ESMF,  only: operator(==), operator(/=), operator(+), operator(-)

! TODO ESMF_GridCompGetInternalState does not have an explicit Fortran interface.
!! Model does not compile with "use ESMF,  only: ESMF_GridCompGetInternalState"
!! Is this okay?

use NUOPC,       only: NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
use NUOPC,       only: NUOPC_CompFilterPhaseMap, NUOPC_CompAttributeGet, NUOPC_CompAttributeAdd
use NUOPC,       only: NUOPC_Advertise, NUOPC_SetAttribute, NUOPC_IsUpdated, NUOPC_Write
use NUOPC,       only: NUOPC_IsConnected, NUOPC_Realize, NUOPC_CompAttributeSet
use NUOPC_Model, only: NUOPC_ModelGet
use NUOPC_Model, only: model_routine_SS           => SetServices
use NUOPC_Model, only: model_label_Advance        => label_Advance
use NUOPC_Model, only: model_label_DataInitialize => label_DataInitialize
use NUOPC_Model, only: model_label_SetRunClock    => label_SetRunClock
use NUOPC_Model, only: model_label_Finalize       => label_Finalize
use NUOPC_Model, only: SetVM

implicit none; private

public SetServices
public SetVM

!> Internal state type with pointers to three types defined by MOM.
type ocean_internalstate_type
  type(ocean_public_type),       pointer :: ocean_public_type_ptr
  type(ocean_state_type),        pointer :: ocean_state_type_ptr
  type(ice_ocean_boundary_type), pointer :: ice_ocean_boundary_type_ptr
end type

!>  Wrapper-derived type required to associate an internal state instance
!! with the ESMF/NUOPC component
type ocean_internalstate_wrapper
  type(ocean_internalstate_type), pointer :: ptr
end type

!> Contains field information
type fld_list_type
  character(len=64) :: stdname
  character(len=64) :: shortname
  character(len=64) :: transferOffer
end type fld_list_type

integer,parameter    :: fldsMax = 100
integer              :: fldsToOcn_num = 0
type (fld_list_type) :: fldsToOcn(fldsMax)
integer              :: fldsFrOcn_num = 0
type (fld_list_type) :: fldsFrOcn(fldsMax)

integer              :: dbug = 0
integer              :: import_slice = 1
integer              :: export_slice = 1
character(len=256)   :: tmpstr
logical              :: write_diagnostics = .false.
logical              :: overwrite_timeslice = .false.
character(len=32)    :: runtype  !< run type
integer              :: logunit  !< stdout logging unit number
logical              :: profile_memory = .true.
logical              :: grid_attach_area = .false.
logical              :: use_coldstart = .true.
logical              :: use_mommesh = .true.
character(len=128)   :: scalar_field_name = ''
integer              :: scalar_field_count = 0
integer              :: scalar_field_idx_grid_nx = 0
integer              :: scalar_field_idx_grid_ny = 0
character(len=*),parameter :: u_FILE_u = &
     __FILE__

#ifdef CESMCOUPLED
logical :: cesm_coupled = .true.
type(ESMF_GeomType_Flag) :: geomtype = ESMF_GEOMTYPE_MESH
#else
logical :: cesm_coupled = .false.
type(ESMF_GeomType_Flag) :: geomtype
#endif
character(len=8) :: restart_mode = 'alarms'

contains

!> NUOPC SetService method is the only public entry point.
!! SetServices registers all of the user-provided subroutines
!! in the module with the NUOPC layer.
!!
!! @param gcomp an ESMF_GridComp object
!! @param rc return code
subroutine SetServices(gcomp, rc)

  type(ESMF_GridComp)  :: gcomp !< an ESMF_GridComp object
  integer, intent(out) :: rc    !< return code

  ! local variables
  character(len=*),parameter  :: subname='(MOM_cap:SetServices)'

  rc = ESMF_SUCCESS

  ! the NUOPC model component will register the generic methods
  call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! switching to IPD versions
  call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
    userRoutine=InitializeP0, phase=0, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! set entry point for methods that require specific implementation
  call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
    phaseLabelList=(/"IPDv03p1"/), userRoutine=InitializeAdvertise, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
    phaseLabelList=(/"IPDv03p3"/), userRoutine=InitializeRealize, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !------------------
  ! attach specializing method(s)
  !------------------

  call NUOPC_CompSpecialize(gcomp, specLabel=model_label_DataInitialize, &
    specRoutine=DataInitialize, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
    specRoutine=ModelAdvance, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
       specRoutine=ModelSetRunClock, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
    specRoutine=ocean_model_finalize, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

end subroutine SetServices

!> First initialize subroutine called by NUOPC.  The purpose
!! is to set which version of the Initialize Phase Definition (IPD)
!! to use.
!!
!! For this MOM cap, we are using IPDv01.
!!
!! @param gcomp an ESMF_GridComp object
!! @param importState an ESMF_State object for import fields
!! @param exportState an ESMF_State object for export fields
!! @param clock an ESMF_Clock object
!! @param rc return code
subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
  type(ESMF_GridComp)   :: gcomp                    !< ESMF_GridComp object
  type(ESMF_State)      :: importState, exportState !< ESMF_State object for
                                                    !! import/export fields
  type(ESMF_Clock)      :: clock                    !< ESMF_Clock object
  integer, intent(out)  :: rc                       !< return code

  ! local variables
  logical                     :: isPresent, isSet
  integer                     :: iostat
  character(len=64)           :: value, logmsg
  character(len=*),parameter  :: subname='(MOM_cap:InitializeP0)'

  rc = ESMF_SUCCESS

  ! Switch to IPDv03 by filtering all other phaseMap entries
  call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
       acceptStringList=(/"IPDv03p"/), rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  write_diagnostics = .false.
  call NUOPC_CompAttributeGet(gcomp, name="DumpFields", value=value, &
       isPresent=isPresent, isSet=isSet, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (isPresent .and. isSet) write_diagnostics=(trim(value)=="true")

  write(logmsg,*) write_diagnostics
  call ESMF_LogWrite('MOM_cap:DumpFields = '//trim(logmsg), ESMF_LOGMSG_INFO)

  overwrite_timeslice = .false.
  call NUOPC_CompAttributeGet(gcomp, name="OverwriteSlice", value=value, &
       isPresent=isPresent, isSet=isSet, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (isPresent .and. isSet) overwrite_timeslice=(trim(value)=="true")
  write(logmsg,*) overwrite_timeslice
  call ESMF_LogWrite('MOM_cap:OverwriteSlice = '//trim(logmsg), ESMF_LOGMSG_INFO)

  profile_memory = .false.
  call NUOPC_CompAttributeGet(gcomp, name="ProfileMemory", value=value, &
       isPresent=isPresent, isSet=isSet, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (isPresent .and. isSet) profile_memory=(trim(value)=="true")
  write(logmsg,*) profile_memory
  call ESMF_LogWrite('MOM_cap:ProfileMemory = '//trim(logmsg), ESMF_LOGMSG_INFO)

  grid_attach_area = .false.
  call NUOPC_CompAttributeGet(gcomp, name="GridAttachArea", value=value, &
       isPresent=isPresent, isSet=isSet, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (isPresent .and. isSet) grid_attach_area=(trim(value)=="true")
  write(logmsg,*) grid_attach_area
  call ESMF_LogWrite('MOM_cap:GridAttachArea = '//trim(logmsg), ESMF_LOGMSG_INFO)

  call NUOPC_CompAttributeGet(gcomp, name='dbug_flag', value=value, isPresent=isPresent, isSet=isSet, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (isPresent .and. isSet) then
   read(value,*) dbug
  end if
  write(logmsg,'(i6)') dbug
  call ESMF_LogWrite('MOM_cap:dbug = '//trim(logmsg), ESMF_LOGMSG_INFO)

  scalar_field_name = ""
  call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=value, &
       isPresent=isPresent, isSet=isSet, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (isPresent .and. isSet) then
     scalar_field_name = trim(value)
     call ESMF_LogWrite('MOM_cap:ScalarFieldName = '//trim(scalar_field_name), ESMF_LOGMSG_INFO)
  endif

  scalar_field_count = 0
  call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=value, &
       isPresent=isPresent, isSet=isSet, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (isPresent .and. isSet) then
     read(value, *, iostat=iostat) scalar_field_count
     if (iostat /= 0) then
       call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg=subname//": ScalarFieldCount not an integer: "//trim(value), &
            line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
     endif
     write(logmsg,*) scalar_field_count
     call ESMF_LogWrite('MOM_cap:ScalarFieldCount = '//trim(logmsg), ESMF_LOGMSG_INFO)
  endif

  scalar_field_idx_grid_nx = 0
  call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=value, &
       isPresent=isPresent, isSet=isSet, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (isPresent .and. isSet) then
     read(value, *, iostat=iostat) scalar_field_idx_grid_nx
     if (iostat /= 0) then
        call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
             msg=subname//": ScalarFieldIdxGridNX not an integer: "//trim(value), &
             line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
     endif
     write(logmsg,*) scalar_field_idx_grid_nx
     call ESMF_LogWrite('MOM_cap:ScalarFieldIdxGridNX = '//trim(logmsg), ESMF_LOGMSG_INFO)
  endif

  scalar_field_idx_grid_ny = 0
  call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=value, &
       isPresent=isPresent, isSet=isSet, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (isPresent .and. isSet) then
     read(value, *, iostat=iostat) scalar_field_idx_grid_ny
     if (iostat /= 0) then
        call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
             msg=subname//": ScalarFieldIdxGridNY not an integer: "//trim(value), &
             line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
     endif
     write(logmsg,*) scalar_field_idx_grid_ny
     call ESMF_LogWrite('MOM_cap:ScalarFieldIdxGridNY = '//trim(logmsg), ESMF_LOGMSG_INFO)
  endif

  use_coldstart = .true.
  call NUOPC_CompAttributeGet(gcomp, name="use_coldstart", value=value, &
       isPresent=isPresent, isSet=isSet, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (isPresent .and. isSet) use_coldstart=(trim(value)=="true")
  write(logmsg,*) use_coldstart
  call ESMF_LogWrite('MOM_cap:use_coldstart = '//trim(logmsg), ESMF_LOGMSG_INFO)

  use_mommesh = .true.
  call NUOPC_CompAttributeGet(gcomp, name="use_mommesh", value=value, &
       isPresent=isPresent, isSet=isSet, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (isPresent .and. isSet) use_mommesh=(trim(value)=="true")
  write(logmsg,*) use_mommesh
  call ESMF_LogWrite('MOM_cap:use_mommesh = '//trim(logmsg), ESMF_LOGMSG_INFO)

  if(use_mommesh)then
    geomtype = ESMF_GEOMTYPE_MESH
    call NUOPC_CompAttributeGet(gcomp, name='mesh_ocn', isPresent=isPresent, isSet=isSet, rc=rc)
      if (.not. isPresent .and. .not. isSet) then
        call ESMF_LogWrite('geomtype set to mesh but mesh_ocn is not specified', ESMF_LOGMSG_INFO)
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
     endif
  else
    geomtype = ESMF_GEOMTYPE_GRID
  endif

end subroutine

!> Called by NUOPC to advertise import and export fields.  "Advertise"
!! simply means that the standard names of all import and export
!! fields are supplied.  The NUOPC layer uses these to match fields
!! between components in the coupled system.
!!
!! @param gcomp an ESMF_GridComp object
!! @param importState an ESMF_State object for import fields
!! @param exportState an ESMF_State object for export fields
!! @param clock an ESMF_Clock object
!! @param rc return code
subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
  type(ESMF_GridComp)            :: gcomp                    !< ESMF_GridComp object
  type(ESMF_State)               :: importState, exportState !< ESMF_State object for
                                                             !! import/export fields
  type(ESMF_Clock)               :: clock                    !< ESMF_Clock object
  integer, intent(out)           :: rc                       !< return code

  ! local variables
  type(ESMF_VM)                          :: vm
  type(ESMF_Time)                        :: MyTime
  type(ESMF_TimeInterval)                :: TINT
  type (ocean_public_type),      pointer :: ocean_public => NULL()
  type (ocean_state_type),       pointer :: ocean_state => NULL()
  type(ice_ocean_boundary_type), pointer :: Ice_ocean_boundary => NULL()
  type(ocean_internalstate_wrapper)      :: ocean_internalstate
  type(ocean_grid_type),         pointer :: ocean_grid => NULL()
  type(directories)                      :: dirs
  type(time_type)                        :: Run_len      !< length of experiment
  type(time_type)                        :: time0        !< Start time of coupled model's calendar.
  type(time_type)                        :: time_start   !< The time at which to initialize the ocean model
  type(time_type)                        :: Time_restart
  type(time_type)                        :: DT
  integer                                :: DT_OCEAN
  integer                                :: isc,iec,jsc,jec
  integer                                :: year=0, month=0, day=0, hour=0, minute=0, second=0
  integer                                :: mpi_comm_mom
  integer                                :: i,n
  character(len=256)                     :: stdname, shortname
  character(len=32)                      :: starttype            ! model start type
  character(len=512)                     :: diro
  character(len=512)                     :: logfile
  character(ESMF_MAXSTR)                 :: cvalue
  logical                                :: isPresent, isPresentDiro, isPresentLogfile, isSet
  logical                                :: existflag
  integer                                :: userRc
  integer                                :: localPet
  integer                                :: iostat
  integer                                :: readunit
  character(len=512)                     :: restartfile          ! Path/Name of restart file
  character(len=2048)                    :: restartfiles         ! Path/Name of restart files
                                                                 ! (same as restartfile if single restart file)
  character(len=*), parameter            :: subname='(MOM_cap:InitializeAdvertise)'
  character(len=32)                      :: calendar
!--------------------------------

  rc = ESMF_SUCCESS

  call ESMF_LogWrite(subname//' enter', ESMF_LOGMSG_INFO)

  allocate(Ice_ocean_boundary)
  !allocate(ocean_state) ! ocean_model_init allocate this pointer
  allocate(ocean_public)
  allocate(ocean_internalstate%ptr)
  ocean_internalstate%ptr%ice_ocean_boundary_type_ptr => Ice_ocean_boundary
  ocean_internalstate%ptr%ocean_public_type_ptr       => ocean_public
  ocean_internalstate%ptr%ocean_state_type_ptr        => ocean_state

  call ESMF_VMGetCurrent(vm, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_VMGet(VM, mpiCommunicator=mpi_comm_mom, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_ClockGet(CLOCK, currTIME=MyTime, TimeStep=TINT,  RC=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_TimeGet (MyTime, YY=YEAR, MM=MONTH, DD=DAY, H=HOUR, M=MINUTE, S=SECOND, RC=rc )
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  CALL ESMF_TimeIntervalGet(TINT, S=DT_OCEAN, RC=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !TODO: next two lines not present in NCAR
  call fms_init(mpi_comm_mom)
  call constants_init
  call field_manager_init

  ! determine the calendar
  if (cesm_coupled) then
     call NUOPC_CompAttributeGet(gcomp, name="calendar", value=cvalue, &
          isPresent=isPresent, isSet=isSet, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     if (isPresent .and. isSet) then
        read(cvalue,*) calendar
        select case (trim(calendar))
           case ("NO_LEAP")
              call set_calendar_type (NOLEAP)
           case ("GREGORIAN")
              call set_calendar_type (GREGORIAN)
           case default
              call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
                 msg=subname//": Calendar not supported in MOM6: "//trim(calendar), &
                 line=__LINE__, file=__FILE__, rcToReturn=rc)
           end select
     else
        call set_calendar_type (NOLEAP)
     endif

  else
     call set_calendar_type (JULIAN)
  endif

  call diag_manager_init

  ! this ocean connector will be driven at set interval
  DT = set_time (DT_OCEAN, 0)
  ! get current time
  time_start = set_date (YEAR,MONTH,DAY,HOUR,MINUTE,SECOND)

  if (is_root_pe()) then
    write(logunit,*) subname//'current time: y,m,d-',year,month,day,'h,m,s=',hour,minute,second
  endif

  ! get start/reference time
  call ESMF_ClockGet(CLOCK, refTime=MyTime, RC=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_TimeGet (MyTime, YY=YEAR, MM=MONTH, DD=DAY, H=HOUR, M=MINUTE, S=SECOND, RC=rc )
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  time0 = set_date (YEAR,MONTH,DAY,HOUR,MINUTE,SECOND)

  if (is_root_pe()) then
    write(logunit,*) subname//'start time: y,m,d-',year,month,day,'h,m,s=',hour,minute,second
  endif

  ! rsd need to figure out how to get this without share code
  !call shr_nuopc_get_component_instance(gcomp, inst_suffix, inst_index)
  !inst_name = "OCN"//trim(inst_suffix)

  ! reset shr logging to my log file
  if (is_root_pe()) then
     call NUOPC_CompAttributeGet(gcomp, name="diro", &
          isPresent=isPresentDiro, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     call NUOPC_CompAttributeGet(gcomp, name="logfile", &
          isPresent=isPresentLogfile, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     if (isPresentDiro .and. isPresentLogfile) then
          call NUOPC_CompAttributeGet(gcomp, name="diro", value=diro, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call NUOPC_CompAttributeGet(gcomp, name="logfile", value=logfile, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          open(newunit=logunit,file=trim(diro)//"/"//trim(logfile))
     else
        logunit = output_unit
     endif
  else
     logunit = output_unit
  endif

  starttype = ""
  call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, &
       isPresent=isPresent, isSet=isSet, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  if (isPresent .and. isSet) then
     read(cvalue,*) starttype
  else
     call ESMF_LogWrite('MOM_cap:start_type unset - using input.nml for restart option', &
          ESMF_LOGMSG_INFO)
  endif

  runtype = ""
  if (trim(starttype) == trim('startup')) then
     runtype = "initial"
  else if (trim(starttype) == trim('continue') ) then
     runtype = "continue"
  else if (trim(starttype) == trim('branch')) then
     runtype = "continue"
  else if (len_trim(starttype) > 0) then
     call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
          msg=subname//": unknown starttype - "//trim(starttype), &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
     return
  endif

  if (len_trim(runtype) > 0) then
     call ESMF_LogWrite('MOM_cap:startup = '//trim(runtype), ESMF_LOGMSG_INFO)
  endif

  restartfile = ""; restartfiles = ""
  if (runtype == "initial") then
    if (cesm_coupled) then
      restartfiles = "n"
    else
      call get_MOM_input(dirs=dirs)
      restartfiles = dirs%input_filename(1:1)
    endif
    call ESMF_LogWrite('MOM_cap:restartfile = '//trim(restartfiles), ESMF_LOGMSG_INFO)

  else if (runtype == "continue") then ! hybrid or branch or continuos runs

     if (cesm_coupled) then
        call ESMF_LogWrite('MOM_cap: restart requested, using rpointer.ocn', ESMF_LOGMSG_WARNING)
        call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        call ESMF_VMGet(vm, localPet=localPet, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        if (localPet == 0) then
           ! this hard coded for rpointer.ocn right now
            open(newunit=readunit, file='rpointer.ocn', form='formatted', status='old', iostat=iostat)
            if (iostat /= 0) then
                call ESMF_LogSetError(ESMF_RC_FILE_OPEN, msg=subname//' ERROR opening rpointer.ocn', &
                     line=__LINE__, file=u_FILE_u, rcToReturn=rc)
                return
            endif
            do
              read(readunit,'(a)', iostat=iostat) restartfile
              if (iostat /= 0) then
                if (len(trim(restartfiles))>1 .and. iostat<0) then
                  exit ! done reading restart files list.
                else
                   call ESMF_LogSetError(ESMF_RC_FILE_READ, msg=subname//' ERROR reading rpointer.ocn', &
                     line=__LINE__, file=u_FILE_u, rcToReturn=rc)
                   return
                endif
              endif
              ! check if the length of restartfiles variable is sufficient:
              if (len(restartfiles)-len(trim(restartfiles)) < len(trim(restartfile))) then
                call MOM_error(FATAL, "Restart file name(s) too long.")
              endif
              restartfiles = trim(restartfiles) // " " // trim(restartfile)
            enddo
            close(readunit)
         endif
         ! broadcast attribute set on master task to all tasks
         call ESMF_VMBroadcast(vm, restartfiles, count=len(restartfiles), rootPet=0, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      else
         call ESMF_LogWrite('MOM_cap: restart requested, use input.nml', ESMF_LOGMSG_WARNING)
       endif

  endif

  ocean_public%is_ocean_pe = .true.
  call ocean_model_init(ocean_public, ocean_state, time0, time_start, input_restart_file=trim(restartfiles))

  call ocean_model_init_sfc(ocean_state, ocean_public)

  call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)

  allocate ( Ice_ocean_boundary% u_flux (isc:iec,jsc:jec),          &
             Ice_ocean_boundary% v_flux (isc:iec,jsc:jec),          &
             Ice_ocean_boundary% t_flux (isc:iec,jsc:jec),          &
             Ice_ocean_boundary% q_flux (isc:iec,jsc:jec),          &
             Ice_ocean_boundary% salt_flux (isc:iec,jsc:jec),       &
             Ice_ocean_boundary% lw_flux (isc:iec,jsc:jec),         &
             Ice_ocean_boundary% sw_flux_vis_dir (isc:iec,jsc:jec), &
             Ice_ocean_boundary% sw_flux_vis_dif (isc:iec,jsc:jec), &
             Ice_ocean_boundary% sw_flux_nir_dir (isc:iec,jsc:jec), &
             Ice_ocean_boundary% sw_flux_nir_dif (isc:iec,jsc:jec), &
             Ice_ocean_boundary% lprec (isc:iec,jsc:jec),           &
             Ice_ocean_boundary% fprec (isc:iec,jsc:jec),           &
             Ice_ocean_boundary% seaice_melt_heat (isc:iec,jsc:jec),&
             Ice_ocean_boundary% seaice_melt (isc:iec,jsc:jec),     &
             Ice_ocean_boundary% mi (isc:iec,jsc:jec),              &
             Ice_ocean_boundary% p (isc:iec,jsc:jec),               &
             Ice_ocean_boundary% lrunoff_hflx (isc:iec,jsc:jec),    &
             Ice_ocean_boundary% frunoff_hflx (isc:iec,jsc:jec),    &
             Ice_ocean_boundary% lrunoff (isc:iec,jsc:jec),       &
             Ice_ocean_boundary% frunoff (isc:iec,jsc:jec))

  Ice_ocean_boundary%u_flux          = 0.0
  Ice_ocean_boundary%v_flux          = 0.0
  Ice_ocean_boundary%t_flux          = 0.0
  Ice_ocean_boundary%q_flux          = 0.0
  Ice_ocean_boundary%salt_flux       = 0.0
  Ice_ocean_boundary%lw_flux         = 0.0
  Ice_ocean_boundary%sw_flux_vis_dir = 0.0
  Ice_ocean_boundary%sw_flux_vis_dif = 0.0
  Ice_ocean_boundary%sw_flux_nir_dir = 0.0
  Ice_ocean_boundary%sw_flux_nir_dif = 0.0
  Ice_ocean_boundary%lprec           = 0.0
  Ice_ocean_boundary%fprec           = 0.0
  Ice_ocean_boundary%seaice_melt     = 0.0
  Ice_ocean_boundary%seaice_melt_heat= 0.0
  Ice_ocean_boundary%mi              = 0.0
  Ice_ocean_boundary%p               = 0.0
  Ice_ocean_boundary%lrunoff_hflx    = 0.0
  Ice_ocean_boundary%frunoff_hflx    = 0.0
  Ice_ocean_boundary%lrunoff         = 0.0
  Ice_ocean_boundary%frunoff         = 0.0

  if (ocean_state%use_waves) then
    Ice_ocean_boundary%num_stk_bands=ocean_state%Waves%NumBands
    allocate ( Ice_ocean_boundary% ustk0 (isc:iec,jsc:jec),         &
               Ice_ocean_boundary% vstk0 (isc:iec,jsc:jec),         &
               Ice_ocean_boundary% ustkb (isc:iec,jsc:jec,Ice_ocean_boundary%num_stk_bands), &
               Ice_ocean_boundary% vstkb (isc:iec,jsc:jec,Ice_ocean_boundary%num_stk_bands), &
               Ice_ocean_boundary%stk_wavenumbers (Ice_ocean_boundary%num_stk_bands))
    Ice_ocean_boundary%ustk0           = 0.0
    Ice_ocean_boundary%vstk0           = 0.0
    Ice_ocean_boundary%stk_wavenumbers = ocean_state%Waves%WaveNum_Cen
    Ice_ocean_boundary%ustkb           = 0.0
    Ice_ocean_boundary%vstkb           = 0.0
  endif

  ocean_internalstate%ptr%ocean_state_type_ptr => ocean_state
  call ESMF_GridCompSetInternalState(gcomp, ocean_internalstate, rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  if (len_trim(scalar_field_name) > 0) then
   call fld_list_add(fldsToOcn_num, fldsToOcn, trim(scalar_field_name), "will_provide")
   call fld_list_add(fldsFrOcn_num, fldsFrOcn, trim(scalar_field_name), "will_provide")
  end if

  if (cesm_coupled) then
    !call fld_list_add(fldsToOcn_num, fldsToOcn, "Sw_lamult"                 , "will provide")
    !call fld_list_add(fldsToOcn_num, fldsToOcn, "Sw_ustokes"                , "will provide")
    !call fld_list_add(fldsToOcn_num, fldsToOcn, "Sw_vstokes"                , "will provide")
    !call fld_list_add(fldsToOcn_num, fldsToOcn, "Sw_hstokes"                , "will provide")
    !call fld_list_add(fldsToOcn_num, fldsToOcn, "Fioi_melth"                , "will provide")
    !call fld_list_add(fldsToOcn_num, fldsToOcn, "Fioi_meltw"                , "will provide")
    !call fld_list_add(fldsFrOcn_num, fldsFrOcn, "So_fswpen"                 , "will provide")
  else
    !call fld_list_add(fldsToOcn_num, fldsToOcn, "mass_of_overlying_sea_ice" , "will provide")
    !call fld_list_add(fldsFrOcn_num, fldsFrOcn, "sea_lev"                   , "will provide")
  endif

  !--------- import fields -------------
  call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_salt_rate"             , "will provide") ! from ice
  call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_zonal_moment_flx"      , "will provide")
  call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_merid_moment_flx"      , "will provide")
  call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_sensi_heat_flx"        , "will provide")
  call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_evap_rate"             , "will provide")
  call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_lw_flx"            , "will provide")
  call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_sw_vis_dir_flx"    , "will provide")
  call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_sw_vis_dif_flx"    , "will provide")
  call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_sw_ir_dir_flx"     , "will provide")
  call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_sw_ir_dif_flx"     , "will provide")
  call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_prec_rate"             , "will provide")
  call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_fprec_rate"            , "will provide")
  call fld_list_add(fldsToOcn_num, fldsToOcn, "inst_pres_height_surface"   , "will provide")
  call fld_list_add(fldsToOcn_num, fldsToOcn, "Foxx_rofl"                  , "will provide") !-> liquid runoff
  call fld_list_add(fldsToOcn_num, fldsToOcn, "Foxx_rofi"                  , "will provide") !-> ice runoff
  call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_fresh_water_to_ocean_rate", "will provide")
  call fld_list_add(fldsToOcn_num, fldsToOcn, "net_heat_flx_to_ocn"        , "will provide")
  !These are not currently used and changing requires a nuopc dictionary change
  !call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_runoff_heat_flx"        , "will provide")
  !call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_calving_heat_flx"       , "will provide")
  if (ocean_state%use_waves) then
    if (Ice_ocean_boundary%num_stk_bands > 3) then
      call MOM_error(FATAL, "Number of Stokes Bands > 3, NUOPC cap not set up for this")
    endif
    call fld_list_add(fldsToOcn_num, fldsToOcn, "eastward_partitioned_stokes_drift_1" , "will provide")
    call fld_list_add(fldsToOcn_num, fldsToOcn, "northward_partitioned_stokes_drift_1", "will provide")
    call fld_list_add(fldsToOcn_num, fldsToOcn, "eastward_partitioned_stokes_drift_2" , "will provide")
    call fld_list_add(fldsToOcn_num, fldsToOcn, "northward_partitioned_stokes_drift_2", "will provide")
    call fld_list_add(fldsToOcn_num, fldsToOcn, "eastward_partitioned_stokes_drift_3" , "will provide")
    call fld_list_add(fldsToOcn_num, fldsToOcn, "northward_partitioned_stokes_drift_3", "will provide")
  endif

  !--------- export fields -------------
  call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocean_mask"                 , "will provide")
  call fld_list_add(fldsFrOcn_num, fldsFrOcn, "sea_surface_temperature"    , "will provide")
  call fld_list_add(fldsFrOcn_num, fldsFrOcn, "s_surf"                     , "will provide")
  call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocn_current_zonal"          , "will provide")
  call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocn_current_merid"          , "will provide")
  call fld_list_add(fldsFrOcn_num, fldsFrOcn, "sea_surface_slope_zonal"    , "will provide")
  call fld_list_add(fldsFrOcn_num, fldsFrOcn, "sea_surface_slope_merid"    , "will provide")
  call fld_list_add(fldsFrOcn_num, fldsFrOcn, "freezing_melting_potential" , "will provide")
  call fld_list_add(fldsFrOcn_num, fldsFrOcn, "So_bldepth"                 , "will provide")

  do n = 1,fldsToOcn_num
    call NUOPC_Advertise(importState, standardName=fldsToOcn(n)%stdname, name=fldsToOcn(n)%shortname, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  enddo

  do n = 1,fldsFrOcn_num
    call NUOPC_Advertise(exportState, standardName=fldsFrOcn(n)%stdname, name=fldsFrOcn(n)%shortname, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  enddo

end subroutine InitializeAdvertise

!> Called by NUOPC to realize import and export fields.  "Realizing" a field
!! means that its grid has been defined and an ESMF_Field object has been
!! created and put into the import or export State.
!!
!! @param gcomp an ESMF_GridComp object
!! @param importState an ESMF_State object for import fields
!! @param exportState an ESMF_State object for export fields
!! @param clock an ESMF_Clock object
!! @param rc return code
subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
  type(ESMF_GridComp)  :: gcomp                    !< ESMF_GridComp object
  type(ESMF_State)     :: importState, exportState !< ESMF_State object for
                                                   !! import/export fields
  type(ESMF_Clock)     :: clock                    !< ESMF_Clock object
  integer, intent(out) :: rc                       !< return code

  ! Local Variables
  type(ESMF_VM)                              :: vm
  type(ESMF_Grid)                            :: gridIn, gridOut
  type(ESMF_Mesh)                            :: Emesh, EmeshTemp
  type(ESMF_DeLayout)                        :: delayout
  type(ESMF_Distgrid)                        :: Distgrid
  type(ESMF_DistGridConnection), allocatable :: connectionList(:)
  type(ESMF_StateItem_Flag)                  :: itemFlag
  type (ocean_public_type),      pointer     :: ocean_public   => NULL()
  type (ocean_state_type),       pointer     :: ocean_state => NULL()
  type(ice_ocean_boundary_type), pointer     :: Ice_ocean_boundary => NULL()
  type(ocean_grid_type)        , pointer     :: ocean_grid
  type(ocean_internalstate_wrapper)          :: ocean_internalstate
  integer                                    :: npet, ntiles
  integer                                    :: nxg, nyg, cnt
  integer                                    :: isc,iec,jsc,jec
  integer, allocatable                       :: xb(:),xe(:),yb(:),ye(:),pe(:)
  integer, allocatable                       :: deBlockList(:,:,:)
  integer, allocatable                       :: petMap(:)
  integer, allocatable                       :: deLabelList(:)
  integer, allocatable                       :: indexList(:)
  integer                                    :: ioff, joff
  integer                                    :: i, j, n, i1, j1, n1, jlast
  integer                                    :: lbnd1,ubnd1,lbnd2,ubnd2
  integer                                    :: lbnd3,ubnd3,lbnd4,ubnd4
  integer                                    :: nblocks_tot
  logical                                    :: found
  integer(ESMF_KIND_I4), pointer             :: dataPtr_mask(:,:)
  real(ESMF_KIND_R8), pointer                :: dataPtr_area(:,:)
  real(ESMF_KIND_R8), pointer                :: dataPtr_xcen(:,:)
  real(ESMF_KIND_R8), pointer                :: dataPtr_ycen(:,:)
  real(ESMF_KIND_R8), pointer                :: dataPtr_xcor(:,:)
  real(ESMF_KIND_R8), pointer                :: dataPtr_ycor(:,:)
  integer                                    :: mpicom
  integer                                    :: localPet
  integer                                    :: lsize
  integer                                    :: ig,jg, ni,nj,k
  integer, allocatable                       :: gindex(:) ! global index space
  character(len=128)                         :: fldname
  character(len=256)                         :: cvalue
  character(len=256)                         :: frmt    ! format specifier for several error msgs
  character(len=512)                         :: err_msg ! error messages
  character(len=*), parameter                :: subname='(MOM_cap:InitializeRealize)'
  integer                         :: spatialDim
  integer                         :: numOwnedElements
  type(ESMF_Array)                :: elemMaskArray
  real(ESMF_KIND_R8)    , pointer :: ownedElemCoords(:)
  real(ESMF_KIND_R8)    , pointer :: lat(:), latMesh(:)
  real(ESMF_KIND_R8)    , pointer :: lon(:), lonMesh(:)
  integer(ESMF_KIND_I4) , pointer :: mask(:), maskMesh(:)
  real(ESMF_KIND_R8)              :: diff_lon, diff_lat
  real                            :: eps_omesh
  !--------------------------------

  rc = ESMF_SUCCESS

  call shr_file_setLogUnit (logunit)

  !----------------------------------------------------------------------------
  ! Get pointers to ocean internal state
  !----------------------------------------------------------------------------

  call ESMF_GridCompGetInternalState(gcomp, ocean_internalstate, rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  Ice_ocean_boundary => ocean_internalstate%ptr%ice_ocean_boundary_type_ptr
  ocean_public       => ocean_internalstate%ptr%ocean_public_type_ptr
  ocean_state        => ocean_internalstate%ptr%ocean_state_type_ptr

  !----------------------------------------------------------------------------
  ! Get mpi information
  !----------------------------------------------------------------------------

  call ESMF_VMGetCurrent(vm, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_VMGet(vm, petCount=npet, mpiCommunicator=mpicom, localPet=localPet, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !---------------------------------
  ! global mom grid size
  !---------------------------------

  call mpp_get_global_domain(ocean_public%domain, xsize=nxg, ysize=nyg)
  write(tmpstr,'(a,2i6)') subname//' nxg,nyg = ',nxg,nyg
  call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

  !---------------------------------
  ! number of tiles per PET, assumed to be 1, and number of pes (tiles) total
  !---------------------------------

  ntiles=mpp_get_ntile_count(ocean_public%domain) ! this is tiles on this pe
  if (ntiles /= 1) then
    rc = ESMF_FAILURE
    call ESMF_LogWrite(subname//' ntiles must be 1', ESMF_LOGMSG_ERROR)
  endif
  ntiles=mpp_get_domain_npes(ocean_public%domain)
  write(tmpstr,'(a,1i6)') subname//' ntiles = ',ntiles
  call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

  !---------------------------------
  ! get start and end indices of each tile and their PET
  !---------------------------------

  allocate(xb(ntiles),xe(ntiles),yb(ntiles),ye(ntiles),pe(ntiles))
  call mpp_get_compute_domains(ocean_public%domain, xbegin=xb, xend=xe, ybegin=yb, yend=ye)
  call mpp_get_pelist(ocean_public%domain, pe)
  if (dbug > 1) then
     do n = 1,ntiles
        write(tmpstr,'(a,6i6)') subname//' tiles ',n,pe(n),xb(n),xe(n),yb(n),ye(n)
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
     enddo
  endif

  !---------------------------------
  ! Create either a grid or a mesh
  !---------------------------------

   !Get the ocean grid and sizes of global and computational domains
   call get_ocean_grid(ocean_state, ocean_grid)

  if (geomtype == ESMF_GEOMTYPE_MESH) then

     !---------------------------------
     ! Create a MOM6 mesh
     !---------------------------------

     call get_global_grid_size(ocean_grid, ni, nj)
     lsize = ( ocean_grid%iec - ocean_grid%isc + 1 ) * ( ocean_grid%jec - ocean_grid%jsc + 1 )

     ! Create the global index space for the computational domain
     allocate(gindex(lsize))
     k = 0
     do j = ocean_grid%jsc, ocean_grid%jec
        jg = j + ocean_grid%jdg_offset
        do i = ocean_grid%isc, ocean_grid%iec
           ig = i + ocean_grid%idg_offset
           k = k + 1 ! Increment position within gindex
           gindex(k) = ni * (jg - 1) + ig
        enddo
     enddo

     DistGrid = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     ! read in the mesh
     call NUOPC_CompAttributeGet(gcomp, name='mesh_ocn', value=cvalue, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     EMeshTemp = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     if (localPet == 0) then
        write(logunit,*)'mesh file for mom6 domain is ',trim(cvalue)
     endif

     ! recreate the mesh using the above distGrid
     EMesh = ESMF_MeshCreate(EMeshTemp, elementDistgrid=Distgrid, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     ! Check for consistency of lat, lon and mask between mesh and mom6 grid
     call ESMF_MeshGet(Emesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     allocate(ownedElemCoords(spatialDim*numOwnedElements))
     allocate(lonMesh(numOwnedElements), lon(numOwnedElements))
     allocate(latMesh(numOwnedElements), lat(numOwnedElements))
     allocate(maskMesh(numOwnedElements), mask(numOwnedElements))

     call ESMF_MeshGet(Emesh, ownedElemCoords=ownedElemCoords, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     do n = 1,numOwnedElements
        lonMesh(n) = ownedElemCoords(2*n-1)
        latMesh(n) = ownedElemCoords(2*n)
     end do

     elemMaskArray = ESMF_ArrayCreate(Distgrid, maskMesh, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     call ESMF_MeshGet(Emesh, elemMaskArray=elemMaskArray, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)
     n = 0
     do j = jsc, jec
       jg = j + ocean_grid%jsc - jsc
       do i = isc, iec
         ig = i + ocean_grid%isc - isc
         n = n+1
         mask(n) = ocean_grid%mask2dT(ig,jg)
         lon(n)  = ocean_grid%geolonT(ig,jg)
         lat(n)  = ocean_grid%geolatT(ig,jg)
       end do
     end do

     eps_omesh = get_eps_omesh(ocean_state)
     do n = 1,numOwnedElements
       diff_lon = abs(mod(lonMesh(n) - lon(n),360.0))
       if (diff_lon > eps_omesh) then
         frmt = "('ERROR: Difference between ESMF Mesh and MOM6 domain coords is "//&
                "greater than parameter EPS_OMESH. n, lonMesh(n), lon(n), diff_lon, "//&
                "EPS_OMESH= ',i8,2(f21.13,3x),2(d21.5))"
         write(err_msg, frmt)n,lonMesh(n),lon(n), diff_lon, eps_omesh
         call MOM_error(FATAL, err_msg)
       end if
       diff_lat = abs(latMesh(n) - lat(n))
       if (diff_lat > eps_omesh) then
         frmt = "('ERROR: Difference between ESMF Mesh and MOM6 domain coords is"//&
                "greater than parameter EPS_OMESH. n, latMesh(n), lat(n), diff_lat, "//&
                "EPS_OMESH= ',i8,2(f21.13,3x),2(d21.5))"
         write(err_msg, frmt)n,latMesh(n),lat(n), diff_lat, eps_omesh
         call MOM_error(FATAL, err_msg)
        end if
        if (abs(maskMesh(n) - mask(n)) > 0) then
          frmt = "('ERROR: ESMF mesh and MOM6 domain masks are inconsistent! - "//&
                 "MOM n, maskMesh(n), mask(n) = ',3(i8,2x))"
          write(err_msg, frmt)n,maskMesh(n),mask(n)
          call MOM_error(FATAL, err_msg)
        end if
     end do

     deallocate(ownedElemCoords)
     deallocate(lonMesh , lon )
     deallocate(latMesh , lat )
     deallocate(maskMesh, mask)
     ! realize the import and export fields using the mesh
     call MOM_RealizeFields(importState, fldsToOcn_num, fldsToOcn, "Ocn import", mesh=Emesh, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     call MOM_RealizeFields(exportState, fldsFrOcn_num, fldsFrOcn, "Ocn export", mesh=Emesh, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

  else if (geomtype == ESMF_GEOMTYPE_GRID) then

     !---------------------------------
     ! create a MOM6 grid
     !---------------------------------

     ! generate delayout and dist_grid

     allocate(deBlockList(2,2,ntiles))
     allocate(petMap(ntiles))
     allocate(deLabelList(ntiles))

     do n = 1, ntiles
       deLabelList(n) = n
       deBlockList(1,1,n) = xb(n)
       deBlockList(1,2,n) = xe(n)
       deBlockList(2,1,n) = yb(n)
       deBlockList(2,2,n) = ye(n)
       petMap(n) = pe(n)
       ! write(tmpstr,'(a,3i8)') subname//' iglo = ',n,deBlockList(1,1,n),deBlockList(1,2,n)
       ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
       ! write(tmpstr,'(a,3i8)') subname//' jglo = ',n,deBlockList(2,1,n),deBlockList(2,2,n)
       ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
       ! write(tmpstr,'(a,2i8)') subname//' pe  = ',n,petMap(n)
       ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
       !--- assume a tile with starting index of 1 has an equivalent wraparound tile on the other side
     enddo

     delayout = ESMF_DELayoutCreate(petMap, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     ! rsd this assumes tripole grid, but sometimes in CESM a bipole
     ! grid is used -- need to introduce conditional logic here

     allocate(connectionList(2))

     ! bipolar boundary condition at top row: nyg
     call ESMF_DistGridConnectionSet(connectionList(1), tileIndexA=1, &
          tileIndexB=1, positionVector=(/nxg+1, 2*nyg+1/), &
          orientationVector=(/-1, -2/), rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     ! periodic boundary condition along first dimension
     call ESMF_DistGridConnectionSet(connectionList(2), tileIndexA=1, &
          tileIndexB=1, positionVector=(/nxg, 0/), rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     distgrid = ESMF_DistGridCreate(minIndex=(/1,1/), maxIndex=(/nxg,nyg/), &
          !        indexflag = ESMF_INDEX_DELOCAL, &
          deBlockList=deBlockList, &
          !        deLabelList=deLabelList, &
          delayout=delayout, &
          connectionList=connectionList, &
          rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     deallocate(xb,xe,yb,ye,pe)
     deallocate(connectionList)
     deallocate(deLabelList)
     deallocate(deBlockList)
     deallocate(petMap)

     call ESMF_DistGridGet(distgrid=distgrid, localDE=0, elementCount=cnt, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     allocate(indexList(cnt))
     write(tmpstr,'(a,i8)') subname//' distgrid cnt= ',cnt
     call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

     call ESMF_DistGridGet(distgrid=distgrid, localDE=0, seqIndexList=indexList, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     write(tmpstr,'(a,4i8)') subname//' distgrid list= ',&
          indexList(1),indexList(cnt),minval(indexList), maxval(indexList)
     call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

     deallocate(IndexList)

     ! create grid

     gridIn = ESMF_GridCreate(distgrid=distgrid, &
          gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/0,1/), &
          coordSys = ESMF_COORDSYS_SPH_DEG, &
          rc = rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     call ESMF_GridAddCoord(gridIn, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     call ESMF_GridAddCoord(gridIn, staggerLoc=ESMF_STAGGERLOC_CORNER, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     call ESMF_GridAddItem(gridIn, itemFlag=ESMF_GRIDITEM_MASK, itemTypeKind=ESMF_TYPEKIND_I4, &
          staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     ! Attach area to the Grid optionally. By default the cell areas are computed.
     if(grid_attach_area) then
        call ESMF_GridAddItem(gridIn, itemFlag=ESMF_GRIDITEM_AREA, itemTypeKind=ESMF_TYPEKIND_R8, &
             staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
     endif

     call ESMF_GridGetCoord(gridIn, coordDim=1, &
          staggerloc=ESMF_STAGGERLOC_CENTER, &
          farrayPtr=dataPtr_xcen, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     call ESMF_GridGetCoord(gridIn, coordDim=2, &
          staggerloc=ESMF_STAGGERLOC_CENTER, &
          farrayPtr=dataPtr_ycen, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     call ESMF_GridGetCoord(gridIn, coordDim=1, &
          staggerloc=ESMF_STAGGERLOC_CORNER, &
          farrayPtr=dataPtr_xcor, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     call ESMF_GridGetCoord(gridIn, coordDim=2, &
          staggerloc=ESMF_STAGGERLOC_CORNER, &
          farrayPtr=dataPtr_ycor, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     call ESMF_GridGetItem(gridIn, itemflag=ESMF_GRIDITEM_MASK, &
          staggerloc=ESMF_STAGGERLOC_CENTER, &
          farrayPtr=dataPtr_mask, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     if(grid_attach_area) then
       call ESMF_GridGetItem(gridIn, itemflag=ESMF_GRIDITEM_AREA, &
             staggerloc=ESMF_STAGGERLOC_CENTER, &
             farrayPtr=dataPtr_area, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
     endif

     ! load up area, mask, center and corner values
     ! area, mask, and centers should be same size in mom and esmf grid
     ! corner points may not be, need to offset corner points by 1 in i and j
     ! retrieve these values directly from ocean_grid, which contains halo
     ! values for j=0 and wrap-around in i. on tripole seam, decomposition
     ! domains are 1 larger in j; to load corner values need to loop one extra row

     call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)

     lbnd1 = lbound(dataPtr_mask,1)
     ubnd1 = ubound(dataPtr_mask,1)
     lbnd2 = lbound(dataPtr_mask,2)
     ubnd2 = ubound(dataPtr_mask,2)

     lbnd3 = lbound(dataPtr_xcor,1)
     ubnd3 = ubound(dataPtr_xcor,1)
     lbnd4 = lbound(dataPtr_xcor,2)
     ubnd4 = ubound(dataPtr_xcor,2)

     write(tmpstr,*) subname//' iscjsc = ',isc,iec,jsc,jec
     call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

     write(tmpstr,*) subname//' lbub12 = ',lbnd1,ubnd1,lbnd2,ubnd2
     call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

     write(tmpstr,*) subname//' lbub34 = ',lbnd3,ubnd3,lbnd4,ubnd4
     call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

     if (iec-isc /= ubnd1-lbnd1 .or. jec-jsc /= ubnd2-lbnd2) then
        call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
             msg=SUBNAME//": fld and grid do not have the same size.", &
             line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
     endif

     do j = jsc, jec
       j1 = j + lbnd2 - jsc
       jg = j + ocean_grid%jsc - jsc
       do i = isc, iec
         i1 = i + lbnd1 - isc
         ig = i + ocean_grid%isc - isc
         dataPtr_mask(i1,j1)  = ocean_grid%mask2dT(ig,jg)
         dataPtr_xcen(i1,j1)  = ocean_grid%geolonT(ig,jg)
         dataPtr_ycen(i1,j1)  = ocean_grid%geolatT(ig,jg)
         if(grid_attach_area) then
           dataPtr_area(i1,j1) = ocean_grid%US%L_to_m**2 * ocean_grid%areaT(ig,jg)
         endif
       enddo
     enddo

     jlast = jec
     if(jec == nyg)jlast = jec+1

     do j = jsc, jlast
       j1 = j + lbnd4 - jsc
       jg = j + ocean_grid%jsc - jsc - 1
       do i = isc, iec
         i1 = i + lbnd3 - isc
         ig = i + ocean_grid%isc - isc - 1
         dataPtr_xcor(i1,j1)  = ocean_grid%geolonBu(ig,jg)
         dataPtr_ycor(i1,j1)  = ocean_grid%geolatBu(ig,jg)
       enddo
     enddo

     write(tmpstr,*) subname//' mask = ',minval(dataPtr_mask),maxval(dataPtr_mask)
     call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

     if(grid_attach_area) then
        write(tmpstr,*) subname//' area = ',minval(dataPtr_area),maxval(dataPtr_area)
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
     endif

     write(tmpstr,*) subname//' xcen = ',minval(dataPtr_xcen),maxval(dataPtr_xcen)
     call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

     write(tmpstr,*) subname//' ycen = ',minval(dataPtr_ycen),maxval(dataPtr_ycen)
     call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

     write(tmpstr,*) subname//' xcor = ',minval(dataPtr_xcor),maxval(dataPtr_xcor)
     call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

     write(tmpstr,*) subname//' ycor = ',minval(dataPtr_ycor),maxval(dataPtr_ycor)
     call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

     gridOut = gridIn ! for now out same as in

     call MOM_RealizeFields(importState, fldsToOcn_num, fldsToOcn, "Ocn import", grid=gridIn, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     call MOM_RealizeFields(exportState, fldsFrOcn_num, fldsFrOcn, "Ocn export", grid=gridOut, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

  endif

  !---------------------------------
  ! set scalar data in export state
  !---------------------------------

  if (len_trim(scalar_field_name) > 0) then
     call State_SetScalar(real(nxg,ESMF_KIND_R8),scalar_field_idx_grid_nx, exportState, localPet, &
         scalar_field_name, scalar_field_count, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     call State_SetScalar(real(nyg,ESMF_KIND_R8),scalar_field_idx_grid_ny, exportState, localPet, &
          scalar_field_name, scalar_field_count, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
  endif

  !---------------------------------
  ! Set module variable geomtype in MOM_cap_methods
  !---------------------------------
  call mom_set_geomtype(geomtype)

  !---------------------------------
  ! write out diagnostics
  !---------------------------------

  !call NUOPC_Write(exportState, fileNamePrefix='post_realize_field_ocn_export_', &
  !     timeslice=1, relaxedFlag=.true., rc=rc)
  !if (ChkErr(rc,__LINE__,u_FILE_u)) return

end subroutine InitializeRealize

!> TODO
!!
!! @param gcomp an ESMF_GridComp object
!! @param rc return code
subroutine DataInitialize(gcomp, rc)
  type(ESMF_GridComp)  :: gcomp !< ESMF_GridComp object
  integer, intent(out) :: rc    !< return code

  ! local variables
  type(ESMF_Clock)                       :: clock
  type(ESMF_State)                       :: importState, exportState
  type(ESMF_Time)                        :: currTime
  type(ESMF_TimeInterval)                :: timeStep
  type(ESMF_StateItem_Flag)              :: itemType
  type (ocean_public_type),      pointer :: ocean_public       => NULL()
  type (ocean_state_type),       pointer :: ocean_state        => NULL()
  type(ice_ocean_boundary_type), pointer :: Ice_ocean_boundary => NULL()
  type(ocean_internalstate_wrapper)      :: ocean_internalstate
  type(ocean_grid_type), pointer         :: ocean_grid
  character(240)                         :: msgString
  character(240)                         :: fldname
  character(240)                         :: timestr
  integer                                :: fieldCount, n
  type(ESMF_Field)                       :: field
  character(len=64),allocatable          :: fieldNameList(:)
  character(len=*),parameter  :: subname='(MOM_cap:DataInitialize)'
  !--------------------------------

  ! query the Component for its clock, importState and exportState
  call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, exportState=exportState, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  call ESMF_TimeGet(currTime,          timestring=timestr, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_GridCompGetInternalState(gcomp, ocean_internalstate, rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  Ice_ocean_boundary => ocean_internalstate%ptr%ice_ocean_boundary_type_ptr
  ocean_public       => ocean_internalstate%ptr%ocean_public_type_ptr
  ocean_state        => ocean_internalstate%ptr%ocean_state_type_ptr
  call get_ocean_grid(ocean_state, ocean_grid)

  call mom_export(ocean_public, ocean_grid, ocean_state, exportState, clock, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  allocate(fieldNameList(fieldCount))
  call ESMF_StateGet(exportState, itemNameList=fieldNameList, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  do n=1, fieldCount
    call ESMF_StateGet(exportState, itemName=fieldNameList(n), field=field, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_SetAttribute(field, name="Updated", value="true", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  enddo
  deallocate(fieldNameList)

  ! check whether all Fields in the exportState are "Updated"
  if (NUOPC_IsUpdated(exportState)) then
    call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="true", rc=rc)
    call ESMF_LogWrite("MOM6 - Initialize-Data-Dependency SATISFIED!!!", ESMF_LOGMSG_INFO)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  endif

  if(write_diagnostics) then
     do n = 1,fldsFrOcn_num
      fldname = fldsFrOcn(n)%shortname
      call ESMF_StateGet(exportState, itemName=trim(fldname), itemType=itemType, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      if (itemType /= ESMF_STATEITEM_NOTFOUND) then
        call ESMF_StateGet(exportState, itemName=trim(fldname), field=field, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        call ESMF_FieldWrite(field, fileName='field_init_ocn_export_'//trim(timestr)//'.nc', &
          timeslice=1, overwrite=overwrite_timeslice, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
     enddo
  endif

end subroutine DataInitialize

!> Called by NUOPC to advance the model a single timestep.
!!
!! @param gcomp an ESMF_GridComp object
!! @param rc return code
subroutine ModelAdvance(gcomp, rc)
  type(ESMF_GridComp)                    :: gcomp !< ESMF_GridComp object
  integer, intent(out)                   :: rc    !< return code

  ! local variables
  integer                                :: userRc
  logical                                :: existflag, isPresent, isSet
  logical                                :: do_advance = .true.
  type(ESMF_Clock)                       :: clock!< ESMF Clock class definition
  type(ESMF_Alarm)                       :: restart_alarm, stop_alarm
  type(ESMF_State)                       :: importState, exportState
  type(ESMF_Time)                        :: currTime
  type(ESMF_TimeInterval)                :: timeStep
  type(ESMF_Time)                        :: startTime
  type(ESMF_TimeInterval)                :: time_elapsed
  integer(ESMF_KIND_I8)                  :: n_interval, time_elapsed_sec
  type(ESMF_Field)                       :: lfield
  type(ESMF_StateItem_Flag)              :: itemType
  character(len=64)                      :: timestamp
  type (ocean_public_type),      pointer :: ocean_public       => NULL()
  type (ocean_state_type),       pointer :: ocean_state        => NULL()
  type(ice_ocean_boundary_type), pointer :: Ice_ocean_boundary => NULL()
  type(ocean_internalstate_wrapper)      :: ocean_internalstate
  type(ocean_grid_type)        , pointer :: ocean_grid
  type(time_type)                        :: Time
  type(time_type)                        :: Time_step_coupled
  type(time_type)                        :: Time_restart_current
  integer                                :: dth, dtm, dts
  integer                                :: nc
  type(ESMF_Time)                        :: MyTime
  integer                                :: seconds, day, year, month, hour, minute
  character(ESMF_MAXSTR)                 :: restartname, cvalue
  character(240)                         :: msgString
  character(ESMF_MAXSTR)                 :: casename
  integer                                :: iostat
  integer                                :: writeunit
  integer                                :: localPet
  type(ESMF_VM)                          :: vm
  integer                                :: n, i
  character(240)                         :: import_timestr, export_timestr
  character(len=128)                     :: fldname
  character(len=*),parameter             :: subname='(MOM_cap:ModelAdvance)'
  character(len=8)                       :: suffix
  integer                                :: num_rest_files

  rc = ESMF_SUCCESS
  if(profile_memory) call ESMF_VMLogMemInfo("Entering MOM Model_ADVANCE: ")

  call shr_file_setLogUnit (logunit)

  ! query the Component for its clock, importState and exportState
  call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
    exportState=exportState, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep

  call ESMF_ClockPrint(clock, options="currTime", &
    preString="------>Advancing OCN from: ", unit=msgString, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  call ESMF_LogWrite(subname//trim(msgString), ESMF_LOGMSG_INFO)

  call ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, &
    timeStep=timeStep, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_TimePrint(currTime + timeStep, &
    preString="--------------------------------> to: ", unit=msgString, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

  call ESMF_TimeGet(currTime,          timestring=import_timestr, rc=rc)
  call ESMF_TimeGet(currTime+timestep, timestring=export_timestr, rc=rc)

  Time_step_coupled = esmf2fms_time(timeStep)
  Time = esmf2fms_time(currTime)

  !---------------
  ! Apply ocean lag for startup runs:
  !---------------

  if (cesm_coupled .or. (.not.use_coldstart)) then
    if (trim(runtype) == "initial") then

      ! Do not call MOM6 timestepping routine if the first cpl tstep of a startup run
      if (currTime == startTime) then
        call ESMF_LogWrite("MOM6 - Skipping the first coupling timestep", ESMF_LOGMSG_INFO)
        do_advance = .false.
      else
        do_advance = .true.
      endif

      if (do_advance) then
         ! If the second cpl tstep of a startup run, step back a cpl tstep and advance for two cpl tsteps
         if (currTime == startTime + timeStep) then
            call ESMF_LogWrite("MOM6 - Stepping back one coupling timestep", ESMF_LOGMSG_INFO)
            Time = esmf2fms_time(currTime-timeStep) ! i.e., startTime

            call ESMF_LogWrite("MOM6 - doubling the coupling timestep", ESMF_LOGMSG_INFO)
            Time_step_coupled = 2 * esmf2fms_time(timeStep)
         endif
      end if

    endif
  endif

  if (do_advance) then

     call ESMF_GridCompGetInternalState(gcomp, ocean_internalstate, rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     Ice_ocean_boundary => ocean_internalstate%ptr%ice_ocean_boundary_type_ptr
     ocean_public       => ocean_internalstate%ptr%ocean_public_type_ptr
     ocean_state        => ocean_internalstate%ptr%ocean_state_type_ptr

     !---------------
     ! Write diagnostics for import
     !---------------

     if (write_diagnostics) then
      do n = 1,fldsToOcn_num
       fldname = fldsToOcn(n)%shortname
       call ESMF_StateGet(importState, itemName=trim(fldname), itemType=itemType, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (itemType /= ESMF_STATEITEM_NOTFOUND) then
         call ESMF_StateGet(importState, itemName=trim(fldname), field=lfield, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return

         call ESMF_FieldWrite(lfield, fileName='field_ocn_import_'//trim(import_timestr)//'.nc', &
           timeslice=1, overwrite=overwrite_timeslice, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
      enddo
     endif

     if (dbug > 0) then
       call state_diagnose(importState,subname//':IS ',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
     end if

     !---------------
     ! Get ocean grid
     !---------------

     call get_ocean_grid(ocean_state, ocean_grid)

     !---------------
     ! Import data
     !---------------

     call mom_import(ocean_public, ocean_grid, importState, ice_ocean_boundary, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     !---------------
     ! Update MOM6
     !---------------

     if(profile_memory) call ESMF_VMLogMemInfo("Entering MOM update_ocean_model: ")
     call update_ocean_model(Ice_ocean_boundary, ocean_state, ocean_public, Time, Time_step_coupled)
     if(profile_memory) call ESMF_VMLogMemInfo("Leaving MOM update_ocean_model: ")

     !---------------
     ! Export Data
     !---------------

     call mom_export(ocean_public, ocean_grid, ocean_state, exportState, clock, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     if (dbug > 0) then
       call state_diagnose(exportState,subname//':ES ',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
     end if
  endif

  !---------------
  ! Get the stop alarm
  !---------------

   call ESMF_ClockGetAlarm(clock, alarmname='stop_alarm', alarm=stop_alarm, rc=rc)
   if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !---------------
  ! If restart alarm exists and is ringing - write restart file
  !---------------

  if (restart_mode == 'alarms') then
     call ESMF_ClockGetAlarm(clock, alarmname='restart_alarm', alarm=restart_alarm, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     if (ESMF_AlarmIsRinging(restart_alarm, rc=rc)) then
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     ! turn off the alarm
     call ESMF_AlarmRingerOff(restart_alarm, rc=rc )
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ! determine restart filename
     call ESMF_ClockGetNextTime(clock, MyTime, rc=rc)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return
     call ESMF_TimeGet (MyTime, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc )
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     if (cesm_coupled) then
        call NUOPC_CompAttributeGet(gcomp, name='case_name', value=casename, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        call ESMF_VMGet(vm, localPet=localPet, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        write(restartname,'(A,".mom6.r.",I4.4,"-",I2.2,"-",I2.2,"-",I5.5)') &
             trim(casename), year, month, day, seconds
        call ESMF_LogWrite("MOM_cap: Writing restart :  "//trim(restartname), ESMF_LOGMSG_INFO)
        ! write restart file(s)
        call ocean_model_restart(ocean_state, restartname=restartname, num_rest_files=num_rest_files)
        if (localPet == 0) then
           ! Write name of restart file in the rpointer file - this is currently hard-coded for the ocean
           open(newunit=writeunit, file='rpointer.ocn', form='formatted', status='unknown', iostat=iostat)
           if (iostat /= 0) then
              call ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
                   msg=subname//' ERROR opening rpointer.ocn', line=__LINE__, file=u_FILE_u, rcToReturn=rc)
              return
           endif
           write(writeunit,'(a)') trim(restartname)//'.nc'

           if (num_rest_files > 1) then
              ! append i.th restart file name to rpointer
              do i=1, num_rest_files-1
                if (i < 10) then
                  write(suffix,'("_",I1)') i
                else
                  write(suffix,'("_",I2)') i
                endif
                write(writeunit,'(a)') trim(restartname) // trim(suffix) // '.nc'
              enddo
           endif
           close(writeunit)
        endif
     else  ! not cesm_coupled
        ! write the final restart without a timestamp
        if (ESMF_AlarmIsRinging(stop_alarm, rc=rc)) then
           write(restartname,'(A)')"MOM.res"
        else
           write(restartname,'(A,I4.4,"-",I2.2,"-",I2.2,"-",I2.2,"-",I2.2,"-",I2.2)') &
                "MOM.res.", year, month, day, hour, minute, seconds
        endif
        call ESMF_LogWrite("MOM_cap: Writing restart :  "//trim(restartname), ESMF_LOGMSG_INFO)

        ! write restart file(s)
        call ocean_model_restart(ocean_state, restartname=restartname)
     endif

     if (is_root_pe()) then
       write(logunit,*) subname//' writing restart file ',trim(restartname)
     endif
   endif
  end if ! restart_mode

  !---------------
  ! Write diagnostics
  !---------------

  if (write_diagnostics) then
     do n = 1,fldsFrOcn_num
      fldname = fldsFrOcn(n)%shortname
      call ESMF_StateGet(exportState, itemName=trim(fldname), itemType=itemType, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      if (itemType /= ESMF_STATEITEM_NOTFOUND) then
        call ESMF_StateGet(exportState, itemName=trim(fldname), field=lfield, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        call ESMF_FieldWrite(lfield, fileName='field_ocn_export_'//trim(export_timestr)//'.nc', &
          timeslice=1, overwrite=overwrite_timeslice, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
     enddo
  endif

  if(profile_memory) call ESMF_VMLogMemInfo("Leaving MOM Model_ADVANCE: ")

end subroutine ModelAdvance


subroutine ModelSetRunClock(gcomp, rc)
  type(ESMF_GridComp)  :: gcomp
  integer, intent(out) :: rc

  ! local variables
  type(ESMF_Clock)         :: mclock, dclock
  type(ESMF_Time)          :: mcurrtime, dcurrtime
  type(ESMF_Time)          :: mstoptime, dstoptime
  type(ESMF_TimeInterval)  :: mtimestep, dtimestep
  character(len=128)       :: mtimestring, dtimestring
  character(len=256)       :: cvalue
  character(len=256)       :: restart_option ! Restart option units
  integer                  :: restart_n      ! Number until restart interval
  integer                  :: restart_ymd    ! Restart date (YYYYMMDD)
  type(ESMF_Alarm)         :: restart_alarm
  type(ESMF_Alarm)         :: stop_alarm
  logical                  :: isPresent, isSet
  logical                  :: first_time = .true.
  character(len=*),parameter :: subname='MOM_cap:(ModelSetRunClock) '
  character(len=256)       :: timestr
  !--------------------------------

  rc = ESMF_SUCCESS

  ! query the Component for its clock, importState and exportState
  call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, &
                     stopTime=dstoptime, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !--------------------------------
  ! check that the current time in the model and driver are the same
  !--------------------------------

  if (mcurrtime /= dcurrtime) then
    call ESMF_TimeGet(dcurrtime, timeString=dtimestring, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(mcurrtime, timeString=mtimestring, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
         msg=subname//": ERROR in time consistency: "//trim(dtimestring)//" != "//trim(mtimestring),  &
         line=__LINE__, file=__FILE__, rcToReturn=rc)
    return
  endif

  !--------------------------------
  ! force model clock currtime and timestep to match driver and set stoptime
  !--------------------------------

  mstoptime = mcurrtime + dtimestep

  call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  if (first_time) then
     !--------------------------------
     ! set restart alarm
     !--------------------------------

     ! defaults
     restart_n = 0
     restart_ymd = 0

     if (cesm_coupled) then

        call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        ! If restart_option is set then must also have set either restart_n or restart_ymd
        call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, &
                isPresent=isPresent, isSet=isSet, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        if (isPresent .and. isSet) then
           read(cvalue,*) restart_n
        endif
        call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, &
             isPresent=isPresent, isSet=isSet, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        if (isPresent .and. isSet) then
           read(cvalue,*) restart_ymd
        endif
        if (restart_n == 0 .and. restart_ymd == 0) then
           call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
                msg=subname//": ERROR both restart_n and restart_ymd are zero for restart_option set ",  &
                line=__LINE__, file=__FILE__, rcToReturn=rc)
           return
        endif
        call ESMF_LogWrite(subname//" Set restart option = "//restart_option, ESMF_LOGMSG_INFO)

     else
        call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, &
             isPresent=isPresent, isSet=isSet, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        ! If restart_n is set and non-zero, then restart_option must be available from config
        if (isPresent .and. isSet) then
          call ESMF_LogWrite(subname//" Restart_n = "//trim(cvalue), ESMF_LOGMSG_INFO)
          read(cvalue,*) restart_n
          if(restart_n /= 0)then
            call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=cvalue, &
                 isPresent=isPresent, isSet=isSet, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            if (isPresent .and. isSet) then
              read(cvalue,*) restart_option
              call ESMF_LogWrite(subname//" Restart_option = "//restart_option, &
                   ESMF_LOGMSG_INFO)
            else
              call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
                   msg=subname//": ERROR both restart_n and restart_option must be set ",  &
                   line=__LINE__, file=__FILE__, rcToReturn=rc)
              return
            endif
            ! not used in nems
            call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, &
                 isPresent=isPresent, isSet=isSet, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            if (isPresent .and. isSet) then
               read(cvalue,*) restart_ymd
               call ESMF_LogWrite(subname//" Restart_ymd = "//trim(cvalue), ESMF_LOGMSG_INFO)
            endif
          else
            ! restart_n is zero, restarts will be written at finalize only (no alarm control)
            restart_mode = 'no_alarms'
            call ESMF_LogWrite(subname//" Restarts will be written at finalize only", ESMF_LOGMSG_INFO)
          endif
        endif
     endif

     if (restart_mode == 'alarms') then
        call AlarmInit(mclock, &
             alarm   = restart_alarm,         &
             option  = trim(restart_option),  &
             opt_n   = restart_n,             &
             opt_ymd = restart_ymd,           &
             RefTime = mcurrTime,             &
             alarmname = 'restart_alarm', rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        call ESMF_LogWrite(subname//" Restart alarm is Created and Set", ESMF_LOGMSG_INFO)
     end if

     ! create a 1-shot alarm at the driver stop time
     stop_alarm = ESMF_AlarmCreate(mclock, ringtime=dstopTime, name = "stop_alarm", rc=rc)
     call ESMF_LogWrite(subname//" Create Stop alarm", ESMF_LOGMSG_INFO)
     if (ChkErr(rc,__LINE__,u_FILE_u)) return

     call ESMF_TimeGet(dstoptime, timestring=timestr, rc=rc)
     call ESMF_LogWrite("Stop Alarm will ring at : "//trim(timestr), ESMF_LOGMSG_INFO)

     first_time = .false.

  endif

  !--------------------------------
  ! Advance model clock to trigger alarms then reset model clock back to currtime
  !--------------------------------

  call ESMF_ClockAdvance(mclock,rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

end subroutine ModelSetRunClock

!===============================================================================

!> Called by NUOPC at the end of the run to clean up.
!!
!! @param gcomp an ESMF_GridComp object
!! @param rc return code
subroutine ocean_model_finalize(gcomp, rc)

  type(ESMF_GridComp)  :: gcomp !< ESMF_GridComp object
  integer, intent(out) :: rc    !< return code

  ! local variables
  type (ocean_public_type),      pointer :: ocean_public
  type (ocean_state_type),       pointer :: ocean_state
  type(ocean_internalstate_wrapper)      :: ocean_internalstate
  type(TIME_TYPE)                        :: Time
  type(ESMF_Clock)                       :: clock
  type(ESMF_Time)                        :: currTime
  type(ESMF_Alarm), allocatable          :: alarmList(:)
  integer                                :: alarmCount
  character(len=64)                      :: timestamp
  logical                                :: write_restart
  character(len=*),parameter  :: subname='(MOM_cap:ocean_model_finalize)'

  write(*,*) 'MOM: --- finalize called ---'
  rc = ESMF_SUCCESS

  call ESMF_GridCompGetInternalState(gcomp, ocean_internalstate, rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  ocean_public => ocean_internalstate%ptr%ocean_public_type_ptr
  ocean_state  => ocean_internalstate%ptr%ocean_state_type_ptr

  call NUOPC_ModelGet(gcomp, modelClock=clock, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return
  Time = esmf2fms_time(currTime)

  ! Do not write a restart unless mode is no_alarms
  if (restart_mode == 'no_alarms') then
     write_restart = .true.
  else
     write_restart = .false.
  end if
  if (write_restart)call ESMF_LogWrite("No Restart Alarm, writing restart at Finalize ", &
                         ESMF_LOGMSG_INFO)

  call ocean_model_end(ocean_public, ocean_State, Time, write_restart=write_restart)
  call field_manager_end()

  call fms_io_exit()
  call fms_end()

  write(*,*) 'MOM: --- completed ---'

end subroutine ocean_model_finalize


!> Set scalar data from state for a particula name
subroutine State_SetScalar(value, scalar_id, State, mytask, scalar_name, scalar_count,  rc)
  real(ESMF_KIND_R8),intent(in)     :: value
  integer,           intent(in)     :: scalar_id
  type(ESMF_State),  intent(inout)  :: State
  integer,           intent(in)     :: mytask
  character(len=*),  intent(in)     :: scalar_name
  integer,           intent(in)     :: scalar_count
  integer,           intent(inout)  :: rc           !< return code

  ! local variables
  type(ESMF_Field)                :: field
  real(ESMF_KIND_R8), pointer     :: farrayptr(:,:)
  character(len=*), parameter     :: subname='(MOM_cap:State_SetScalar)'
  !--------------------------------------------------------

  rc = ESMF_SUCCESS

  call ESMF_StateGet(State, itemName=trim(scalar_name), field=field, rc=rc)
  if (ChkErr(rc,__LINE__,u_FILE_u)) return

  if (mytask == 0) then
    call ESMF_FieldGet(field, farrayPtr=farrayptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (scalar_id < 0 .or. scalar_id > scalar_count) then
       call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg=subname//": ERROR in scalar_id", line=__LINE__, file=__FILE__, rcToReturn=rc)
       return
    endif

    farrayptr(scalar_id,1) = value
  endif

end subroutine State_SetScalar

!> Realize the import and export fields using either a grid or a mesh.
subroutine MOM_RealizeFields(state, nfields, field_defs, tag, grid, mesh, rc)
  type(ESMF_State)    , intent(inout)        :: state !< ESMF_State object for
                                                      !! import/export fields.
  integer             , intent(in)           :: nfields !< Number of fields.
  type(fld_list_type) , intent(inout)        :: field_defs(:) !< Structure with field's
                                                              !! information.
  character(len=*)    , intent(in)           :: tag !< Import or export.
  type(ESMF_Grid)     , intent(in), optional :: grid!< ESMF grid.
  type(ESMF_Mesh)     , intent(in), optional :: mesh!< ESMF mesh.
  integer             , intent(inout)        :: rc  !< Return code.

  ! local variables
  integer                     :: i
  type(ESMF_Field)            :: field
  real(ESMF_KIND_R8), pointer :: fldptr1d(:)   ! for mesh
  real(ESMF_KIND_R8), pointer :: fldptr2d(:,:) ! for grid
  character(len=*),parameter  :: subname='(MOM_cap:MOM_RealizeFields)'
  !--------------------------------------------------------

  rc = ESMF_SUCCESS

  do i = 1, nfields

    if (NUOPC_IsConnected(state, fieldName=field_defs(i)%shortname)) then

      if (field_defs(i)%shortname == scalar_field_name) then

        call ESMF_LogWrite(subname // tag // " Field "// trim(field_defs(i)%stdname) // " is connected on root pe.", &
          ESMF_LOGMSG_INFO)

        call SetScalarField(field, rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

      else

        call ESMF_LogWrite(subname // tag // " Field "// trim(field_defs(i)%stdname) // " is connected.", &
          ESMF_LOGMSG_INFO)

        if (present(grid)) then

           field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, indexflag=ESMF_INDEX_DELOCAL, &
                name=field_defs(i)%shortname, rc=rc)
           if (ChkErr(rc,__LINE__,u_FILE_u)) return

           ! initialize fldptr to zero
           call ESMF_FieldGet(field, farrayPtr=fldptr2d, rc=rc)
           if (ChkErr(rc,__LINE__,u_FILE_u)) return
           fldptr2d(:,:) = 0.0

        else if (present(mesh)) then

           field = ESMF_FieldCreate(mesh=mesh, typekind=ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
                name=field_defs(i)%shortname, rc=rc)
           if (ChkErr(rc,__LINE__,u_FILE_u)) return

           ! initialize fldptr to zero
           call ESMF_FieldGet(field, farrayPtr=fldptr1d, rc=rc)
           if (ChkErr(rc,__LINE__,u_FILE_u)) return
           fldptr1d(:) = 0.0

        endif

      endif

      ! Realize connected field
      call NUOPC_Realize(state, field=field, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

    else ! field is not connected

      call ESMF_LogWrite(subname // tag // " Field "// trim(field_defs(i)%stdname) // " is not connected.", &
        ESMF_LOGMSG_INFO)

      ! remove a not connected Field from State
      call ESMF_StateRemove(state, (/field_defs(i)%shortname/), rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

    endif

  enddo

contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  subroutine SetScalarField(field, rc)

    ! create a field with scalar data on the root pe
    type(ESMF_Field), intent(inout)  :: field
    integer,          intent(inout)  :: rc

    ! local variables
    type(ESMF_Distgrid) :: distgrid
    type(ESMF_Grid)     :: grid
    character(len=*), parameter :: subname='(MOM_cap:SetScalarField)'

    rc = ESMF_SUCCESS

    ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
    distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    grid = ESMF_GridCreate(distgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! num of scalar values
    field = ESMF_FieldCreate(name=trim(scalar_field_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
         ungriddedLBound=(/1/), ungriddedUBound=(/scalar_field_count/), gridToFieldMap=(/2/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine SetScalarField

end subroutine MOM_RealizeFields

!===============================================================================

!> Set up list of field information
subroutine fld_list_add(num, fldlist, stdname, transferOffer, shortname)
  integer,                    intent(inout) :: num
  type(fld_list_type),        intent(inout) :: fldlist(:)
  character(len=*),           intent(in)    :: stdname
  character(len=*),           intent(in)    :: transferOffer
  character(len=*), optional, intent(in)    :: shortname

  ! local variables
  integer :: rc
  character(len=*), parameter :: subname='(MOM_cap:fld_list_add)'

  ! fill in the new entry
  num = num + 1
  if (num > fldsMax) then
     call ESMF_LogSetError(ESMF_RC_VAL_OUTOFRANGE, &
          msg=trim(subname)//": ERROR number of field exceeded fldsMax: "//trim(stdname), &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
     return
  endif

  fldlist(num)%stdname        = trim(stdname)
  if (present(shortname)) then
     fldlist(num)%shortname   = trim(shortname)
  else
     fldlist(num)%shortname   = trim(stdname)
  endif
  fldlist(num)%transferOffer  = trim(transferOffer)

end subroutine fld_list_add


#ifndef CESMCOUPLED
subroutine shr_file_setLogUnit(nunit)
  integer, intent(in) :: nunit
  ! do nothing for this stub - its just here to replace
  ! having cppdefs in the main program
end subroutine shr_file_setLogUnit

subroutine shr_file_getLogUnit(nunit)
  integer, intent(in) :: nunit
  ! do nothing for this stub - its just here to replace
  ! having cppdefs in the main program
end subroutine shr_file_getLogUnit
#endif

!>
!! @page nuopc_cap NUOPC Cap
!! @author Fei Liu (fei.liu@gmail.com)
!! @date 5/10/13 Original documentation
!! @author Rocky Dunlap (rocky.dunlap@noaa.gov)
!! @date 1/12/17 Moved to doxygen
!! @date 2/28/19 Rewrote for unified cap
!! @tableofcontents
!!
!! @section Overview Overview
!!
!! **This MOM cap has been tested with MOM6.**
!!
!! This document describes the MOM NUOPC "cap", which is a light weight software layer that is
!! required when the [MOM ocean model](https://github.com/NOAA-GFDL/MOM6/tree/dev/master)
!! is used in [National Unified Operation Prediction Capability]
!! (http://www.earthsystemcog.org/projects/nuopc) (NUOPC) coupled systems. Also see the
!! [MOM wiki](https://github.com/NOAA-GFDL/MOM6-Examples/wiki) for more documentation.
!!
!! NUOPC is a software layer built on top of the [Earth System Modeling
!! Framework] (https://www.earthsystemcog.org/projects/esmf) (ESMF).
!! ESMF is a high-performance modeling framework that provides
!! data structures, interfaces, and operations suited for building coupled models
!! from a set of components. NUOPC refines the capabilities of ESMF by providing
!! a more precise definition of what it means for a model to be a component and
!! how components should interact and share data in a coupled system. The NUOPC
!! Layer software is designed to work with typical high-performance models in the
!! Earth sciences domain, most of which are written in Fortran and are based on a
!! distributed memory model of parallelism (MPI).
!!
!! A NUOPC "cap" is a Fortran module that serves as the interface to a model
!! when it's used in a NUOPC-based coupled system.
!! The term "cap" is used because it is a light weight software layer that sits on top
!! of model code, making calls into it and exposing model data structures in a
!! standard way.
!!
!! The MOM cap package includes the cap code itself (MOM_cap.F90, MOM_cap_methods.F90
!! and MOM_cap_time.F90), a set of time utilities (time_utils.F90) for converting between ESMF and FMS
!! time type and two modules MOM_ocean_model_nuopc.F90 and MOM_surface_forcing_nuopc.F90. MOM_surface_forcing_nuopc.F90
!! converts the input ESMF data (import data) to a MOM-specific data type (surface_forcing_CS).
!! MOM_ocean_model_nuopc.F90 contains routines for initialization, update and finalization of the ocean model state.
!!
!! @subsection CapSubroutines Cap Subroutines
!!
!! The MOM cap modules contains a set of subroutines that are required
!! by NUOPC.  These subroutines are called by the NUOPC infrastructure according
!! to a predefined calling sequence.  Some subroutines are called during
!! initialization of the coupled system, some during the run of the coupled
!! system, and some during finalization of the coupled system.
!!
!! The initialization sequence is the most complex and is governed by the NUOPC technical rules.
!! Details about the initialization sequence can be found in the [NUOPC Reference Manual]
!! (http://www.earthsystemmodeling.org/esmf_releases/last_built/NUOPC_refdoc/).
!! The cap requires beta snapshot ESMF v8.0.0bs16 or later.
!!
!! The following table summarizes the NUOPC-required subroutines that appear in the
!! MOM cap.  The "Phase" column says whether the subroutine is called during the
!! initialization, run, or finalize part of the coupled system run.
!!
!!<table>
!!<tr><th> Phase  <th> MOM Cap Subroutine  <th> Description
!!<tr>
!!  <td> Init
!!  <td> [InitializeP0] (@ref MOM_cap_mod::initializep0)
!!  <td> Sets the Initialize Phase Definition (IPD) version to use
!!<tr>
!!  <td> Init
!!  <td> [InitializeAdvertise] (@ref MOM_cap_mod::initializeadvertise)
!!  <td> Advertises standard names of import and export fields
!!<tr>
!!  <td> Init
!!  <td> [InitializeRealize] (@ref MOM_cap_mod::initializerealize)
!!  <td> Creates an ESMF_Grid or ESMF_Mesh as well as ESMF_Fields for import and export fields
!!<tr>
!!  <td> Run
!!  <td> [ModelAdvance] (@ref MOM_cap_mod::modeladvance)
!!  <td> Advances the model by a timestep
!!<tr>
!!  <td> Final
!!  <td> [Finalize] (@ref MOM_cap_mod::ocean_model_finalize)
!!  <td> Cleans up
!!</table>
!!
!!
!! @section UnderlyingModelInterfaces Underlying Model Interfaces
!!
!!
!! @subsection DomainCreation Domain Creation
!!
!! The cap can accomodate a MOM tripolar grid which is represented either as a 2D `ESMF_Grid` or
!! as a 1D `ESMF_Mesh`. Other MOM grids (e.g. a bipolar grid) can be represented as a 1d `ESMF_Mesh` only.
!! Coupling fields are placed on either the `ESMF_Grid` or `ESMF_Mesh`.
!! Note that for either the `ESMF_Grid` or `ESMF_Mesh` representation, the fields are translated into
!! a 2D MOM specific surface boundary type and the distinction between the two is no longer there.
!! Calls related to creating the grid are located in the [InitializeRealize]
!! (@ref MOM_cap_mod::initializerealize) subroutine, which is called by the NUOPC infrastructure
!! during the intialization sequence.
!!
!! The cap determines parameters for setting up the grid by calling subroutines in the
!! `mpp_domains_mod` module. The global domain size is determined by calling `mpp_get_global_domain()`.
!! A check is in place to ensure that there is only a single tile in the domain (the
!! cap is currently limited to one tile; multi-tile mosaics are not supported).  The
!! decomposition across processors is determined via calls to `mpp_get_compute_domains()`
!! (to retrieve decomposition block indices) and `mpp_get_pelist()` (to determine how
!! blocks are assigned to processors).
!!
!! The `ESMF_Grid` is created in several steps:
!!  - an `ESMF_DELayout` is created based on the pelist from MOM
!!  - an `ESMF_DistGrid` is created over the global index space. Connections are set
!!    up so that the index space is periodic in the first dimension and has a
!!    fold at the top for the bipole. The decompostion blocks are also passed in
!!    along with the `ESMF_DELayout` mentioned above.
!!  - an `ESMF_Grid` is then created by passing in the above `ESMF_DistGrid`.
!!  - masks, areas, center (tlat, tlon), and corner (ulat, ulon) coordinates are then added to the `ESMF_Grid`
!!    by retrieving those fields from the MOM datatype `ocean_grid` elements.
!!
!! The `ESMF_Mesh` is also created in several steps:
!!   - the target mesh is generated offline.
!!   - a temporary mesh is created from an input file specified by the config variable `mesh_ocn`.
!!     the mesh has a distribution that is automatically generated by ESMF when reading in the mesh
!!   - an `ESMF_DistGrid` is created from the global index space for the computational domain.
!!   - the final `ESMF_Mesh` is then created by distributing the temporary mesh using the created `ESMF_DistGrid`.
!!
!!
!! @subsection Initialization Initialization
!!
!! During the [InitializeAdvertise] (@ref MOM_cap_mod::initializeadvertise) phase, calls are
!! made to MOM's native initialization subroutines, including `fms_init()`, `constants_init()`,
!! `field_manager_init()`, `diag_manager_init()`, and `set_calendar_type()`.  The MPI communicator
!! is pulled in through the ESMF VM object for the MOM component. The dt and start time are set
!! from parameters from the incoming ESMF clock with calls to `set_time()` and `set_date().`
!!
!!
!! @subsection Run Run
!!
!! The [ModelAdvance] (@ref MOM_cap_mod::modeladvance) subroutine is called by the NUOPC
!! infrastructure when it's time for MOM to advance in time. During this subroutine, there is a
!! call into the MOM update routine:
!!
!!      call update_ocean_model(Ice_ocean_boundary, Ocean_state, Ocean_public, Time, Time_step_coupled)
!!
!! Priori to the call to `update_ocean_model()`, the cap performs these steps
!! - the `Time` and `Time_step_coupled` parameters, based on FMS types, are derived from the incoming ESMF clock
!! - diagnostics are optionally written to files `field_ocn_import_*`, one for each import field
!! - mom_import is called and translates to the ESMF input data to a MOM specific data type
!!    - momentum flux vectors are rotated to internal grid
!!
!! After the call to `update_ocean_model()`, the cap performs these steps:
!! - mom_export is called
!!   - the `ocean_mask` export is set to match that of the internal MOM mask
!!   - the `freezing_melting_potential` export is converted from J m-2 to W m-2 by dividing by the coupling interval
!!   - vector rotations are applied to the `ocean_current_zonal` and `ocean_current_merid` exports, back to lat-lon grid
!!   - diagnostics are optionally written to files `field_ocn_export_*`, one for each export field
!!   - optionally, a call is made to `ocean_model_restart()` at the interval `restart_interval`
!!
!! @subsubsection VectorRotations Vector Rotations
!!
!! Vector rotations are applied to incoming momentum fluxes (from regular lat-lon to tripolar grid) and
!! outgoing ocean currents (from tripolar to regular lat-lon). The rotation angles are provided
!! from the native MOM grid by a call to `get_ocean_grid(Ocean_grid)`.
!! The cosine and sine of the rotation angle are:
!!
!!     ocean_grid%cos_rot(i,j)
!!     ocean_grid%sin_rot(i,j)
!!
!! The rotation of momentum flux from regular lat-lon to tripolar is:
!! \f[
!! \begin{bmatrix}
!! \tau_x' \\
!! \tau_y'
!! \end{bmatrix} =
!! \begin{bmatrix}
!!  cos \theta   & sin \theta \\
!!  -sin \theta  & cos \theta
!! \end{bmatrix} *
!! \begin{bmatrix}
!! \tau_x \\
!! \tau_y
!! \end{bmatrix}
!! \f]
!!
!! The rotation of ocean current from tripolar to regular lat-lon is:
!! \f[
!! \begin{bmatrix}
!! u' \\
!! v'
!! \end{bmatrix} =
!! \begin{bmatrix}
!!  cos \theta   & -sin \theta \\
!!  sin \theta  & cos \theta
!! \end{bmatrix} *
!! \begin{bmatrix}
!! u \\
!! v
!! \end{bmatrix}
!! \f]
!! @subsection Finalization Finalization
!!
!! NUOPC infrastructure calls [ocean_model_finalize] (@ref MOM_cap_mod::ocean_model_finalize)
!! at the end of the run. This subroutine is a hook to call into MOM's native shutdown
!! procedures:
!!
!!     call ocean_model_end (ocean_public, ocean_State, Time)
!!     call diag_manager_end(Time )
!!     call field_manager_end
!!     call fms_io_exit
!!     call fms_end
!!
!! @section ModelFields Model Fields
!!
!! The following tables list the import and export fields currently set up in the MOM cap.
!!
!! @subsection ImportFields Import Fields
!!
!! <table>
!! <tr>
!!     <th>Standard Name</td>
!!     <th>Units</td>
!!     <th>Model Variable</td>
!!     <th>Description</td>
!!     <th>Notes</td>
!! <tr>
!!     <td>inst_pres_height_surface</td>
!!     <td>Pa</td>
!!     <td>p</td>
!!     <td>pressure of overlying sea ice and atmosphere</td>
!!     <td></td>
!! </tr>
!! <tr>
!!     <td>mass_of_overlying_sea_ice</td>
!!     <td>kg</td>
!!     <td>mi</td>
!!     <td>mass of overlying sea ice</td>
!!     <td></td>
!! </tr>
!! <tr>
!!     <td>seaice_melt_heat</td>
!!     <td>W m-2</td>
!!     <td>seaice_melt_heat</td>
!!     <td>sea ice and snow melt heat flux</td>
!!     <td></td>
!! </tr>
!! <tr>
!!     <td>seaice_melt</td>
!!     <td>kg m-2 s-1</td>
!!     <td>seaice_melt</td>
!!     <td>water flux due to sea ice and snow melting</td>
!!     <td></td>
!! </tr>
!! <tr>
!!     <td>mean_calving_heat_flx</td>
!!     <td>W m-2</td>
!!     <td>calving_hflx</td>
!!     <td>heat flux, relative to 0C, of frozen land water into ocean</td>
!!     <td></td>
!! </tr>
!! <tr>
!!     <td>mean_calving_rate</td>
!!     <td>kg m-2 s-1</td>
!!     <td>calving</td>
!!     <td>mass flux of frozen runoff</td>
!!     <td></td>
!! </tr>
!! <tr>
!!     <td>mean_evap_rate</td>
!!     <td>kg m-2 s-1</td>
!!     <td>q_flux</td>
!!     <td>specific humidity flux</td>
!!     <td></td>
!! </tr>
!! <tr>
!!     <td>mean_fprec_rate</td>
!!     <td>kg m-2 s-1</td>
!!     <td>fprec</td>
!!     <td>mass flux of frozen precip</td>
!!     <td></td>
!! </tr>
!! <tr>
!!     <td>mean_merid_moment_flx</td>
!!     <td>Pa</td>
!!     <td>v_flux</td>
!!     <td>j-directed wind stress into ocean</td>
!!     <td>[vector rotation] (@ref VectorRotations) applied - lat-lon to tripolar</td>
!! </tr>
!! <tr>
!!     <td>mean_net_lw_flx</td>
!!     <td>W m-2</td>
!!     <td>lw_flux</td>
!!     <td>long wave radiation</td>
!!     <td></td>
!! </tr>
!! <tr>
!!     <td>mean_net_sw_ir_dif_flx</td>
!!     <td>W m-2</td>
!!     <td>sw_flux_nir_dif</td>
!!     <td>diffuse near IR shortwave radiation</td>
!!     <td></td>
!! </tr>
!! <tr>
!!     <td>mean_net_sw_ir_dir_flx</td>
!!     <td>W m-2</td>
!!     <td>sw_flux_nir_dir</td>
!!     <td>direct near IR shortwave radiation</td>
!!     <td></td>
!! </tr>
!! <tr>
!!     <td>mean_net_sw_vis_dif_flx</td>
!!     <td>W m-2</td>
!!     <td>sw_flux_vis_dif</td>
!!     <td>diffuse visible shortware radiation</td>
!!     <td></td>
!! </tr>
!! <tr>
!!     <td>mean_net_sw_vis_dir_flx</td>
!!     <td>W m-2</td>
!!     <td>sw_flux_vis_dir</td>
!!     <td>direct visible shortware radiation</td>
!!     <td></td>
!! </tr>
!! <tr>
!!     <td>mean_prec_rate</td>
!!     <td>kg m-2 s-1</td>
!!     <td>lprec</td>
!!     <td>mass flux of liquid precip</td>
!!     <td></td>
!! </tr>
!! <tr>
!!     <td>mean_runoff_heat_flx</td>
!!     <td>W m-2</td>
!!     <td>runoff_hflx</td>
!!     <td>heat flux, relative to 0C, of liquid land water into ocean</td>
!!     <td></td>
!! </tr>
!! <tr>
!!     <td>mean_runoff_rate</td>
!!     <td>kg m-2 s-1</td>
!!     <td>runoff</td>
!!     <td>mass flux of liquid runoff</td>
!!     <td></td>
!! </tr>
!! <tr>
!!     <td>mean_salt_rate</td>
!!     <td>kg m-2 s-1</td>
!!     <td>salt_flux</td>
!!     <td>salt flux</td>
!!     <td></td>
!! </tr>
!! <tr>
!!     <td>mean_sensi_heat_flx</td>
!!     <td>W m-2</td>
!!     <td>t_flux</td>
!!     <td>sensible heat flux into ocean</td>
!!     <td></td>
!! </tr>
!! <tr>
!!     <td>mean_zonal_moment_flx</td>
!!     <td>Pa</td>
!!     <td>u_flux</td>
!!     <td>i-directed wind stress into ocean</td>
!!     <td>[vector rotation] (@ref VectorRotations) applied - lat-lon to tripolar</td>
!! </tr>
!! </table>
!!
!! @subsection ExportField Export Fields
!!
!! Export fields are populated from the `ocean_public` parameter (type `ocean_public_type`)
!! after the call to `update_ocean_model()`.
!!
!! <table>
!!   <tr>
!!       <th>Standard Name</th>
!!       <th>Units</th>
!!       <th>Model Variable</th>
!!       <th>Description</th>
!!       <th>Notes</th>
!!   </tr>
!!   <tr>
!!       <td>freezing_melting_potential</td>
!!       <td>W m-2</td>
!!       <td>combination of frazil and melt_potential</td>
!!       <td>cap converts model units (J m-2) to (W m-2) for export</td>
!!       <td></td>
!!   </tr>
!!   <tr>
!!       <td>ocean_mask</td>
!!       <td></td>
!!       <td></td>
!!       <td>ocean mask</td>
!!       <td></td>
!!   </tr>
!!   <tr>
!!       <td>ocn_current_merid</td>
!!       <td>m s-1</td>
!!       <td>v_surf</td>
!!       <td>j-directed surface velocity on u-cell</td>
!!       <td>[vector rotation] (@ref VectorRotations) applied - tripolar to lat-lon</td>
!!   </tr>
!!   <tr>
!!       <td>ocn_current_zonal</td>
!!       <td>m s-1</td>
!!       <td>u_surf</td>
!!       <td>i-directed surface velocity on u-cell</td>
!!       <td>[vector rotation] (@ref VectorRotations) applied - tripolar to lat-lon</td>
!!   </tr>
!!   <tr>
!!       <td>s_surf</td>
!!       <td>psu</td>
!!       <td>s_surf</td>
!!       <td>sea surface salinity on t-cell</td>
!!       <td></td>
!!   </tr>
!!   <tr>
!!       <td>sea_surface_temperature</td>
!!       <td>K</td>
!!       <td>t_surf</td>
!!       <td>sea surface temperature on t-cell</td>
!!       <td></td>
!!   </tr>
!!   <tr>
!!       <td>sea_surface_slope_zonal</td>
!!       <td>unitless</td>
!!       <td>created from ssh</td>
!!       <td>sea surface zonal slope</td>
!!       <td></td>
!!   </tr>
!!   <tr>
!!       <td>sea_surface_slope_merid</td>
!!       <td>unitless</td>
!!       <td>created from ssh</td>
!!       <td>sea surface meridional slope</td>
!!       <td></td>
!!   </tr>
!!   <tr>
!!       <td>so_bldepth</td>
!!       <td>m</td>
!!       <td>obld</td>
!!       <td>ocean surface boundary layer depth</td>
!!       <td></td>
!!   </tr>
!! </table>
!!
!! @subsection MemoryManagement Memory Management
!!
!! The MOM cap has an internal state type with pointers to three
!! types defined by MOM. There is also a small wrapper derived type
!! required to associate an internal state instance
!! with the ESMF/NUOPC component:
!!
!!     type ocean_internalstate_type
!!        type(ocean_public_type),       pointer :: ocean_public_type_ptr
!!        type(ocean_state_type),        pointer :: ocean_state_type_ptr
!!        type(ice_ocean_boundary_type), pointer :: ice_ocean_boundary_type_ptr
!!     end type
!!
!!     type ocean_internalstate_wrapper
!!        type(ocean_internalstate_type), pointer :: ptr
!!     end type
!!
!! The member of type `ocean_public_type` stores ocean surface fields used during the coupling.
!! The member of type `ocean_state_type` is required by the ocean driver,
!! although its internals are private (not to be used by the coupling directly).
!! This type is passed to the ocean init and update routines
!! so that it can maintain state there if desired.
!! The member of type `ice_ocean_boundary_type` is populated by this cap
!! with incoming coupling fields from other components. These three derived types are allocated during the
!! [InitializeAdvertise] (@ref MOM_cap_mod::initializeadvertise) phase.  Also during that
!! phase, the `ice_ocean_boundary` type members are all allocated using bounds retrieved
!! from `mpp_get_compute_domain()`.
!!
!! During the [InitializeRealize] (@ref MOM_cap_mod::initializerealize) phase,
!! `ESMF_Field`s are created for each of the coupling fields in the `ice_ocean_boundary`
!! and `ocean_public_type` members of the internal state. These fields directly reference into the members of
!! the `ice_ocean_boundary` and `ocean_public_type` so that memory-to-memory copies are not required to move
!! data from the cap's import and export states to the memory areas used internally
!! by MOM.
!!
!! @subsection IO I/O
!!
!! The cap can optionally output coupling fields for diagnostic purposes if the ESMF attribute
!! "DumpFields" has been set to "true". In this case the cap will write out NetCDF files
!! with names "field_ocn_import_<fieldname>.nc" and "field_ocn_export_<fieldname>.nc".
!! Additionally, calls will be made to the cap subroutine [dumpMomInternal]
!! (@ref MOM_cap_mod::dumpmominternal) to write out model internal fields to files
!! named "field_ocn_internal_<fieldname>.nc".  In all cases these NetCDF files will
!! contain a time series of field data.
!!
!! @section RuntimeConfiguration Runtime Configuration
!!
!! At runtime, the MOM cap can be configured with several options provided
!! as ESMF attributes.  Attributes can be set in the cap by the NUOPC Driver
!! above this cap, or in some systems ESMF attributes are set by
!! reading in from a configuration file.  The available attributes are:
!!
!! * `DumpFields` - when set to "true", write out diagnostic NetCDF files for import/export/internal fields
!! * `ProfileMemory` - when set to "true", write out memory usage information to the ESMF log files; this
!!    information is written when entering and leaving the [ModelAdvance]
!!    (@ref MOM_cap_mod::modeladvance) subroutine and before and after the call to
!!   `update_ocean_model()`.
!! * `restart_interval` - integer number of seconds indicating the interval at
!!    which to call `ocean_model_restart()`; no restarts written if set to 0

end module MOM_cap_mod
