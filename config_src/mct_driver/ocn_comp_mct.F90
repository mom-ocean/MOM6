!> This is the main driver for MOM6 in CIME
module ocn_comp_mct

! This file is part of MOM6. See LICENSE.md for the license.

! mct modules
use ESMF,                only: ESMF_clock, ESMF_time, ESMF_timeInterval, ESMF_TimeInc
use ESMF,                only: ESMF_ClockGet, ESMF_TimeGet, ESMF_TimeIntervalGet
use seq_cdata_mod,       only: seq_cdata, seq_cdata_setptrs
use seq_flds_mod,        only: ice_ncat, seq_flds_i2o_per_cat
use mct_mod,             only: mct_gsMap, mct_gsmap_init, mct_gsMap_lsize, &
                               mct_gsmap_orderedpoints
use mct_mod,             only: mct_aVect, mct_aVect_init, mct_aVect_zero, &
                               mct_aVect_nRattr
use mct_mod,             only: mct_gGrid, mct_gGrid_init, mct_gGrid_importRAttr, &
                               mct_gGrid_importIAttr
use mct_mod,             only : mct_avect_indexra, mct_aVect_clean
use seq_flds_mod,        only: seq_flds_x2o_fields, seq_flds_o2x_fields, seq_flds_dom_coord, &
                               seq_flds_dom_other
use seq_infodata_mod,    only: seq_infodata_type, seq_infodata_GetData, &
                               seq_infodata_start_type_start, seq_infodata_start_type_cont,   &
                               seq_infodata_start_type_brnch, seq_infodata_PutData
use seq_comm_mct,        only: seq_comm_name, seq_comm_inst, seq_comm_suffix
use seq_timemgr_mod,     only: seq_timemgr_EClockGetData, seq_timemgr_RestartAlarmIsOn
use perf_mod,            only: t_startf, t_stopf
use shr_kind_mod,        only: shr_kind_r8

! MOM6 modules
use ocean_model_mod,    only: ocean_state_type, ocean_public_type, ocean_model_init_sfc
use ocean_model_mod,    only: ocean_model_init, get_state_pointers, ocean_model_restart
use ocean_model_mod,    only: ice_ocean_boundary_type, update_ocean_model
use MOM_domains,        only: MOM_infra_init, num_pes, root_pe, pe_here
use MOM_domains,        only : pass_var, AGRID
use MOM_grid,           only: ocean_grid_type, get_global_grid_size
use MOM_variables,      only: surface
use MOM_error_handler,  only: MOM_error, FATAL, is_root_pe
use MOM_time_manager,   only: time_type, set_date, set_time, set_calendar_type, NOLEAP
use MOM_file_parser,    only: get_param, log_version, param_file_type
use MOM_get_input,      only: Get_MOM_Input, directories

! By default make data private
implicit none; private
  ! Public member functions
  public :: ocn_init_mct
  public :: ocn_run_mct
  public :: ocn_final_mct
  ! Flag for debugging
  logical, parameter :: debug=.true.

  !> Structure with MCT attribute vectors and indices
  type cpl_indices

    ! ocean to coupler
    integer :: o2x_So_t         !< Surface potential temperature (deg C)
    integer :: o2x_So_u         !< Surface zonal velocity (m/s)
    integer :: o2x_So_v         !< Surface meridional velocity (m/s)
    integer :: o2x_So_s         !< Surface salinity (PSU)
    integer :: o2x_So_dhdx      !< Zonal slope in the sea surface height
    integer :: o2x_So_dhdy      !< Meridional lope in the sea surface height
    integer :: o2x_So_bldepth   !< Boundary layer depth (m)
    integer :: o2x_Fioo_q       !< Heat flux?
    integer :: o2x_Faoo_fco2_ocn!< CO2 flux
    integer :: o2x_Faoo_fdms_ocn!< DMS flux

    ! coupler to ocean
    integer :: x2o_Si_ifrac        !< Fractional ice wrt ocean
    integer :: x2o_So_duu10n       !< 10m wind speed squared (m^2/s^2)
    integer :: x2o_Sa_pslv         !< Sea-level pressure (Pa)
    integer :: x2o_Sa_co2prog      !< Bottom atm level prognostic CO2
    integer :: x2o_Sa_co2diag      !< Bottom atm level diagnostic CO2
    integer :: x2o_Sw_lamult       !< Wave model langmuir multiplier
    integer :: x2o_Sw_ustokes      !< Surface Stokes drift, x-component
    integer :: x2o_Sw_vstokes      !< Surface Stokes drift, y-component
    integer :: x2o_Foxx_taux       !< Zonal wind stress (W/m2)
    integer :: x2o_Foxx_tauy       !< Meridonal wind stress (W/m2)
    integer :: x2o_Foxx_swnet      !< Net short-wave heat flux (W/m2)
    integer :: x2o_Foxx_sen        !< Sensible heat flux (W/m2)
    integer :: x2o_Foxx_lat        !< Latent heat flux  (W/m2)
    integer :: x2o_Foxx_lwup       !< Longwave radiation, up (W/m2)
    integer :: x2o_Faxa_lwdn       !< Longwave radiation, down (W/m2)
    integer :: x2o_Fioi_melth      !< Heat flux from snow & ice melt (W/m2)
    integer :: x2o_Fioi_meltw      !< Snow melt flux (kg/m2/s)
    integer :: x2o_Fioi_bcpho      !< Black Carbon hydrophobic release
                                   !! from sea ice component
    integer :: x2o_Fioi_bcphi      !< Black Carbon hydrophilic release from
                                   !! sea ice component
    integer :: x2o_Fioi_flxdst     !< Dust release from sea ice component
    integer :: x2o_Fioi_salt       !< Salt flux    (kg(salt)/m2/s)
    integer :: x2o_Foxx_evap       !< Evaporation flux  (kg/m2/s)
    integer :: x2o_Faxa_prec       !< Total precipitation flux (kg/m2/s)
    integer :: x2o_Faxa_snow       !< Water flux due to snow (kg/m2/s)
    integer :: x2o_Faxa_rain       !< Water flux due to rain (kg/m2/s)
    integer :: x2o_Faxa_bcphidry   !< Black   Carbon hydrophilic dry deposition
    integer :: x2o_Faxa_bcphodry   !< Black   Carbon hydrophobic dry deposition
    integer :: x2o_Faxa_bcphiwet   !< Black   Carbon hydrophilic wet deposition
    integer :: x2o_Faxa_ocphidry   !< Organic Carbon hydrophilic dry deposition
    integer :: x2o_Faxa_ocphodry   !< Organic Carbon hydrophobic dry deposition
    integer :: x2o_Faxa_ocphiwet   !< Organic Carbon hydrophilic dry deposition
    integer :: x2o_Faxa_dstwet1    !< Size 1 dust -- wet deposition
    integer :: x2o_Faxa_dstwet2    !< Size 2 dust -- wet deposition
    integer :: x2o_Faxa_dstwet3    !< Size 3 dust -- wet deposition
    integer :: x2o_Faxa_dstwet4    !< Size 4 dust -- wet deposition
    integer :: x2o_Faxa_dstdry1    !< Size 1 dust -- dry deposition
    integer :: x2o_Faxa_dstdry2    !< Size 2 dust -- dry deposition
    integer :: x2o_Faxa_dstdry3    !< Size 3 dust -- dry deposition
    integer :: x2o_Faxa_dstdry4    !< Size 4 dust -- dry deposition
    integer :: x2o_Foxx_rofl       !< River runoff flux (kg/m2/s)
    integer :: x2o_Foxx_rofi       !< Ice runoff flux (kg/m2/s)

    ! optional per thickness category fields

    integer, dimension(:), allocatable :: x2o_frac_col !< Fraction of ocean cell,
                                                       !! per column
    integer, dimension(:), allocatable :: x2o_fracr_col!< Fraction of ocean cell used
                                                       !! in radiation computations,
                                                       !! per column
    integer, dimension(:), allocatable :: x2o_qsw_fracr_col !< qsw * fracr, per column

  end type cpl_indices

  !> Control structure for this module
  type MCT_MOM_Data

    type(ocean_state_type), pointer  :: ocn_state => NULL()   !< The private state of ocean
    type(ocean_public_type), pointer :: ocn_public => NULL()  !< The public state of ocean
    type(ocean_grid_type), pointer   :: grid => NULL()        !< The grid structure
    type(surface), pointer           :: ocn_surface => NULL() !< The ocean surface state
    type(ice_ocean_boundary_type)    :: ice_ocean_boundary    !< The ice ocean boundary type
    type(seq_infodata_type), pointer :: infodata              !< The input info type
    type(cpl_indices), public        :: ind                   !< Variable IDs
    ! runtime params
    logical :: sw_decomp   !< Controls whether shortwave is decomposed into four components
    real :: c1, c2, c3, c4 !< Coeffs. used in the shortwave decomposition

  end type MCT_MOM_Data

  type(MCT_MOM_Data) :: glb                                   !< global structure

contains

!> Initializes MOM6
subroutine ocn_init_mct( EClock, cdata_o, x2o_o, o2x_o, NLFilename )
  type(ESMF_Clock),             intent(inout) :: EClock      !< Time and time step ? \todo Why must this
                                                             !! be intent(inout)?
  type(seq_cdata)             , intent(inout) :: cdata_o     !< Input parameters
  type(mct_aVect)             , intent(inout) :: x2o_o       !< Fluxes from coupler to ocean, computed by ocean
  type(mct_aVect)             , intent(inout) :: o2x_o       !< Fluxes from ocean to coupler, computed by ocean
  character(len=*), optional  , intent(in)    :: NLFilename  !< Namelist filename

  !  local variables
  type(time_type)     :: time_init         !< Start time of coupled model's calendar
  type(time_type)     :: time_in           !< Time at the beginning of the first ocn coupling interval
  type(ESMF_time)     :: current_time      !< Current time
  type(ESMF_time)     :: time_in_ESMF      !< Time after first ocean coupling interval (of type ESMF_time)
  type(ESMF_timeInterval) :: ocn_cpl_interval !< Ocean coupling interval
  integer             :: year, month, day, hour, minute, seconds, seconds_n, seconds_d, rc
  character(len=384)  :: runid             !< Run ID
  character(len=384)  :: runtype           !< Run type
  character(len=32)   :: starttype         !< infodata start type
  integer             :: mpicom_ocn        !< MPI ocn communicator
  integer             :: npes, pe0         !< # of processors and current processor
  integer             :: i, errorCode
  integer             :: lsize, nsend, nrecv
  logical             :: ldiag_cpl = .false.
  integer             :: isc, iec, jsc, jec, ni, nj     !< Indices for the start and end of the domain
                                                        !! in the x and y dir., respectively.
  ! runi-time params
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


  ! instance control vars (these are local for now)
  integer(kind=4)     :: inst_index
  character(len=16)   :: inst_name
  character(len=16)   :: inst_suffix

  !!!DANGER!!!: change the following vars with the corresponding MOM6 vars
  integer :: km=1                   !< Number of vertical levels
  integer :: nx_block=0, ny_block=0 !< Size of block domain in x,y dir including ghost cells
  integer :: max_blocks_clinic=0    !< Max. number of blocks per processor in each distribution
  integer :: ncouple_per_day
  logical :: lsend_precip_fact      !< If T,send precip_fact to cpl for use in fw balance
                                    !! (partially-coupled option)
  character(len=128) :: err_msg     !< Error message

  ! set (actually, get from mct) the cdata pointers:
  call seq_cdata_setptrs(cdata_o, id=MOM_MCT_ID, mpicom=mpicom_ocn, &
                         gsMap=MOM_MCT_gsMap, dom=MOM_MCT_dom, infodata=glb%infodata)

  ! Determine attribute vector indices
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

  call t_startf('MOM_init')

  ! Initialize MOM6 comm
  call MOM_infra_init(mpicom_ocn)

  call set_calendar_type(NOLEAP)  !TODO: confirm this

  ! Get the ESMF clock instance (assigned by CESM for MOM6)
  call ESMF_ClockGet(EClock, currTime=current_time, rc=rc)

  ! Get the initial CESM time
  call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
  time_init = set_date(year, month, day, hour, minute, seconds, err_msg=err_msg)

  ! Compute time_in: time at the beginning of the first ocn coupling interval
  !   (In CESM, ocean component is lagged by one ocean coupling interval, so:
  !   time_in = time_init + ocn_cpl_interval )
  call ESMF_ClockGet(EClock, TimeStep=ocn_cpl_interval, rc=rc)
  time_in_ESMF = ESMF_TimeInc(current_time, ocn_cpl_interval)
  call ESMF_TimeGet(time_in_ESMF, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
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
    call ESMF_ClockGet(EClock, TimeStep=ocn_cpl_interval, rc=rc)
    call ESMF_TimeIntervalGet(ocn_cpl_interval, yy=year, mm=month, d=day, s=seconds, sn=seconds_n, sd=seconds_d, rc=rc)
    write(6,*) 'ocn_init_mct, time step: y,m,d-',year,month,day,'s,sn,sd=',seconds,seconds_n,seconds_d
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
  call get_param(param_file, mdl, "SW_DECOMP", glb%sw_decomp, &
                 "If True, read coeffs c1, c2, c3 and c4 and decompose" // &
                 "the net shortwave radiation (SW) into four components:\n" // &
                 "visible, direct shortwave  = c1 * SW \n" // &
                 "visible, diffuse shortwave = c2 * SW \n" // &
                 "near-IR, direct shortwave  = c3 * SW \n" // &
                 "near-IR, diffuse shortwave = c4 * SW", default=.true.)
  if (glb%sw_decomp) then
    call get_param(param_file, mdl, "SW_c1", glb%c1, &
                  "Coeff. used to convert net shortwave rad. into \n"//&
                  "visible, direct shortwave.", units="nondim", default=0.285)
    call get_param(param_file, mdl, "SW_c2", glb%c2, &
                  "Coeff. used to convert net shortwave rad. into \n"//&
                  "visible, diffuse shortwave.", units="nondim", default=0.285)
    call get_param(param_file, mdl, "SW_c3", glb%c3, &
                  "Coeff. used to convert net shortwave rad. into \n"//&
                  "near-IR, direct shortwave.", units="nondim", default=0.215)
    call get_param(param_file, mdl, "SW_c4", glb%c4, &
                  "Coeff. used to convert net shortwave rad. into \n"//&
                  "near-IR, diffuse shortwave.", units="nondim", default=0.215)
  else
    glb%c1 = 0.0; glb%c2 = 0.0; glb%c3 = 0.0; glb%c4 = 0.0
  endif

  ! Initialize the MOM6 model
  call ocean_model_init(glb%ocn_public, glb%ocn_state, time_init, time_in)

  ! Initialize ocn_state%state out of sight
  call ocean_model_init_sfc(glb%ocn_state, glb%ocn_public)

  ! store pointers to components inside MOM
  call get_state_pointers(glb%ocn_state, grid=glb%grid)

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
  ncouple_per_day = seconds_in_day / ocn_cpl_dt
  mom_cpl_dt = seconds_in_day / ncouple_per_day
  if (mom_cpl_dt /= ocn_cpl_dt) then
     write(*,*) 'ERROR mom_cpl_dt and ocn_cpl_dt must be identical'
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

  ! Allocate ice_ocean_boundary using global indexing without halos
  isc = glb%grid%isc + glb%grid%idg_offset
  iec = glb%grid%iec + glb%grid%idg_offset
  jsc = glb%grid%jsc + glb%grid%jdg_offset
  jec = glb%grid%jec + glb%grid%jdg_offset
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

end subroutine ocn_init_mct

!> Determines attribute vector indices
subroutine coupler_indices_init(ind)

  type(cpl_indices), intent(inout) :: ind !< Structure with coupler indices
                                          !! and vectors

  ! Local Variables
  type(mct_aVect) :: o2x      !< Array with ocean to coupler data
  type(mct_aVect) :: x2o      !< Array with coupler to ocean data

  integer          :: ncat    !< Thickness category index
  character(len=2) :: cncat   !< Character version of ncat
  integer          :: ncol    !< Column index
  integer          :: mcog_ncols !< Number of ice thickness categories?
  integer          :: lmcog_flds_sent !< Used to convert per thickness
                                      !! category fields?

  ! create temporary attribute vectors
  call mct_aVect_init(x2o, rList=seq_flds_x2o_fields, lsize=1)
  call mct_aVect_init(o2x, rList=seq_flds_o2x_fields, lsize=1)

  ! ocean to coupler
  ind%o2x_So_t          = mct_avect_indexra(o2x,'So_t')
  ind%o2x_So_u          = mct_avect_indexra(o2x,'So_u')
  ind%o2x_So_v          = mct_avect_indexra(o2x,'So_v')
  ind%o2x_So_s          = mct_avect_indexra(o2x,'So_s')
  ind%o2x_So_dhdx       = mct_avect_indexra(o2x,'So_dhdx')
  ind%o2x_So_dhdy       = mct_avect_indexra(o2x,'So_dhdy')
  ! QL, 150526, to wav, boundary layer depth
  ind%o2x_So_bldepth    = mct_avect_indexra(o2x,'So_bldepth')
  ind%o2x_Fioo_q        = mct_avect_indexra(o2x,'Fioo_q')
  ind%o2x_Faoo_fco2_ocn = mct_avect_indexra(o2x,'Faoo_fco2_ocn',perrWith='quiet')
  ind%o2x_Faoo_fdms_ocn = mct_avect_indexra(o2x,'Faoo_fdms_ocn',perrWith='quiet')

  ! coupler to ocean
  ind%x2o_Si_ifrac      = mct_avect_indexra(x2o,'Si_ifrac')
  ind%x2o_Sa_pslv       = mct_avect_indexra(x2o,'Sa_pslv')
  ind%x2o_So_duu10n     = mct_avect_indexra(x2o,'So_duu10n')
  ! QL, 150526, from wav
  ind%x2o_Sw_lamult     = mct_avect_indexra(x2o,'Sw_lamult')
  ind%x2o_Sw_ustokes    = mct_avect_indexra(x2o,'Sw_ustokes')
  ind%x2o_Sw_vstokes    = mct_avect_indexra(x2o,'Sw_vstokes')
  ind%x2o_Foxx_tauy     = mct_avect_indexra(x2o,'Foxx_tauy')
  ind%x2o_Foxx_taux     = mct_avect_indexra(x2o,'Foxx_taux')
  ind%x2o_Foxx_swnet    = mct_avect_indexra(x2o,'Foxx_swnet')
  ind%x2o_Foxx_lat      = mct_avect_indexra(x2o,'Foxx_lat')
  ind%x2o_Foxx_sen      = mct_avect_indexra(x2o,'Foxx_sen')
  ind%x2o_Foxx_lwup     = mct_avect_indexra(x2o,'Foxx_lwup')
  ind%x2o_Faxa_lwdn     = mct_avect_indexra(x2o,'Faxa_lwdn')
  ind%x2o_Fioi_melth    = mct_avect_indexra(x2o,'Fioi_melth')
  ind%x2o_Fioi_meltw    = mct_avect_indexra(x2o,'Fioi_meltw')
  ind%x2o_Fioi_salt     = mct_avect_indexra(x2o,'Fioi_salt')
  ind%x2o_Fioi_bcpho    = mct_avect_indexra(x2o,'Fioi_bcpho')
  ind%x2o_Fioi_bcphi    = mct_avect_indexra(x2o,'Fioi_bcphi')
  ind%x2o_Fioi_flxdst   = mct_avect_indexra(x2o,'Fioi_flxdst')
  ind%x2o_Faxa_prec     = mct_avect_indexra(x2o,'Faxa_prec')
  ind%x2o_Faxa_snow     = mct_avect_indexra(x2o,'Faxa_snow')
  ind%x2o_Faxa_rain     = mct_avect_indexra(x2o,'Faxa_rain')
  ind%x2o_Foxx_evap     = mct_avect_indexra(x2o,'Foxx_evap')
  ind%x2o_Foxx_rofl     = mct_avect_indexra(x2o,'Foxx_rofl')
  ind%x2o_Foxx_rofi     = mct_avect_indexra(x2o,'Foxx_rofi')
  ind%x2o_Faxa_bcphidry = mct_avect_indexra(x2o,'Faxa_bcphidry')
  ind%x2o_Faxa_bcphodry = mct_avect_indexra(x2o,'Faxa_bcphodry')
  ind%x2o_Faxa_bcphiwet = mct_avect_indexra(x2o,'Faxa_bcphiwet')
  ind%x2o_Faxa_ocphidry = mct_avect_indexra(x2o,'Faxa_ocphidry')
  ind%x2o_Faxa_ocphodry = mct_avect_indexra(x2o,'Faxa_ocphodry')
  ind%x2o_Faxa_ocphiwet = mct_avect_indexra(x2o,'Faxa_ocphiwet')
  ind%x2o_Faxa_dstdry1  = mct_avect_indexra(x2o,'Faxa_dstdry1')
  ind%x2o_Faxa_dstdry2  = mct_avect_indexra(x2o,'Faxa_dstdry2')
  ind%x2o_Faxa_dstdry3  = mct_avect_indexra(x2o,'Faxa_dstdry3')
  ind%x2o_Faxa_dstdry4  = mct_avect_indexra(x2o,'Faxa_dstdry4')
  ind%x2o_Faxa_dstwet1  = mct_avect_indexra(x2o,'Faxa_dstwet1')
  ind%x2o_Faxa_dstwet2  = mct_avect_indexra(x2o,'Faxa_dstwet2')
  ind%x2o_Faxa_dstwet3  = mct_avect_indexra(x2o,'Faxa_dstwet3')
  ind%x2o_Faxa_dstwet4  = mct_avect_indexra(x2o,'Faxa_dstwet4')
  ind%x2o_Sa_co2prog    = mct_avect_indexra(x2o,'Sa_co2prog',perrWith='quiet')
  ind%x2o_Sa_co2diag    = mct_avect_indexra(x2o,'Sa_co2diag',perrWith='quiet')
  ! optional per thickness category fields

  ! convert cpl indices to mcog column indices
  ! this implementation only handles columns due to ice thickness categories
  lmcog_flds_sent = seq_flds_i2o_per_cat

  if (seq_flds_i2o_per_cat) then
    mcog_ncols = ice_ncat+1
    allocate(ind%x2o_frac_col(mcog_ncols))
    allocate(ind%x2o_fracr_col(mcog_ncols))
    allocate(ind%x2o_qsw_fracr_col(mcog_ncols))
    ncol = 1
    ind%x2o_frac_col(ncol)        = mct_avect_indexra(x2o,'Sf_afrac')
    ind%x2o_fracr_col(ncol)       = mct_avect_indexra(x2o,'Sf_afracr')
    ind%x2o_qsw_fracr_col(ncol)   = mct_avect_indexra(x2o,'Foxx_swnet_afracr')

    do ncat = 1, ice_ncat
      write(cncat,'(i2.2)') ncat
      ncol = ncat+1
      ind%x2o_frac_col(ncol)      = mct_avect_indexra(x2o,'Si_ifrac_'//cncat)
      ind%x2o_fracr_col(ncol)     = ind%x2o_frac_col(ncol)
      ind%x2o_qsw_fracr_col(ncol) = mct_avect_indexra(x2o,'PFioi_swpen_ifrac_'//cncat)
    enddo
  else
    mcog_ncols = 1
  endif

  call mct_aVect_clean(x2o)
  call mct_aVect_clean(o2x)

end subroutine coupler_indices_init

!> Maps outgoing ocean data to MCT buffer
subroutine ocn_export(ind, ocn_public, grid, o2x)
  type(cpl_indices),       intent(inout) :: ind        !< Structure with coupler
                                                       !! indices and vectors
  type(ocean_public_type), intent(in)    :: ocn_public !< Ocean surface state
  type(ocean_grid_type),   intent(in)    :: grid       !< Ocean model grid
  real(kind=8),            intent(inout) :: o2x(:,:)   !< MCT outgoing bugger
  ! Local variables
  real, dimension(grid%isd:grid%ied,grid%jsd:grid%jed) :: ssh !< Local copy of sea_lev with updated halo
  integer :: i, j, n, ig, jg  !< Grid indices
  real :: slp_L, slp_R, slp_C, slope, u_min, u_max

  ! Copy from ocn_public to o2x. ocn_public uses global indexing with no halos.
  ! The mask comes from "grid" that uses the usual MOM domain that has halos
  ! and does not use global indexing.
  n = 0
  do j=grid%jsc, grid%jec
    jg = j + grid%jdg_offset
    do i=grid%isc,grid%iec
      n = n+1
      ig = i + grid%idg_offset
      ! surface temperature in Kelvin
      o2x(ind%o2x_So_t, n) = ocn_public%t_surf(ig,jg)  * grid%mask2dT(i,j)
      o2x(ind%o2x_So_s, n) = ocn_public%s_surf(ig,jg) * grid%mask2dT(i,j)
      o2x(ind%o2x_So_u, n) = ocn_public%u_surf(ig,jg) * grid%mask2dT(i,j)
      o2x(ind%o2x_So_v, n) = ocn_public%v_surf(ig,jg) * grid%mask2dT(i,j)
      ! Make a copy of ssh in order to do a halo update. We use the usual MOM domain
      ! in order to update halos. i.e. does not use global indexing.
      ssh(i,j) = ocn_public%sea_lev(ig,jg)
    end do
  end do

  ! Update halo of ssh so we can calculate gradients
  call pass_var(ssh, grid%domain)


  ! d/dx ssh
  n = 0
  do j=grid%jsc, grid%jec ; do i=grid%isc,grid%iec
    n = n+1
    ! This is a simple second-order difference
    ! o2x(ind%o2x_So_dhdx, n) = 0.5 * (ssh(i+1,j) + ssh(i-1,j)) * grid%IdxT(i,j) * grid%mask2dT(i,j)
    ! This is a PLM slope which might be less prone to the A-grid null mode
    slp_L = (ssh(i,j) - ssh(i-1,j)) * grid%mask2dCu(I-1,j)
    !if (grid%mask2dCu(I-1,j)==0.) slp_L = 0.
    slp_R = (ssh(i+1,j) - ssh(i,j)) * grid%mask2dCu(I,j)
    !if (grid%mask2dCu(I,j)==0.) slp_R = 0.
    slp_C = 0.5 * (slp_L + slp_R)
    if ( (slp_L * slp_R) > 0.0 ) then
      ! This limits the slope so that the edge values are bounded by the
      ! two cell averages spanning the edge.
      u_min = min( ssh(i-1,j), ssh(i,j), ssh(i+1,j) )
      u_max = max( ssh(i-1,j), ssh(i,j), ssh(i+1,j) )
      slope = sign( min( abs(slp_C), 2.*min( ssh(i,j) - u_min, u_max - ssh(i,j) ) ), slp_C )
    else
      ! Extrema in the mean values require a PCM reconstruction avoid generating
      ! larger extreme values.
      slope = 0.0
    end if
    o2x(ind%o2x_So_dhdx, n) = slope * grid%IdxT(i,j) * grid%mask2dT(i,j)
  end do; end do

  ! d/dy ssh
  do j=grid%jsc, grid%jec ; do i=grid%isc,grid%iec
    ! This is a simple second-order difference
    ! o2x(ind%o2x_So_dhdy, n) = 0.5 * (ssh(i,j+1) + ssh(i,j-1)) * grid%IdyT(i,j) * grid%mask2dT(i,j)
    ! This is a PLM slope which might be less prone to the A-grid null mode
    slp_L = ssh(i,j) - ssh(i,j-1)
    slp_R = ssh(i,j+1) - ssh(i,j)
    slp_L=0.; slp_R=0.
    slp_C = 0.5 * (slp_L + slp_R)
    if ((slp_L * slp_R) > 0.0) then
      ! This limits the slope so that the edge values are bounded by the
      ! two cell averages spanning the edge.
      u_min = min( ssh(i,j-1), ssh(i,j), ssh(i,j+1) )
      u_max = max( ssh(i,j-1), ssh(i,j), ssh(i,j+1) )
      slope = sign( min( abs(slp_C), 2.*min( ssh(i,j) - u_min, u_max - ssh(i,j) ) ), slp_C )
    else
      ! Extrema in the mean values require a PCM reconstruction avoid generating
      ! larger extreme values.
      slope = 0.0
    end if
    o2x(ind%o2x_So_dhdy, n) = slope * grid%IdyT(i,j) * grid%mask2dT(i,j)
  end do; end do

end subroutine ocn_export

!> Fills the ice ocean boundary type
subroutine fill_ice_ocean_bnd(ice_ocean_boundary, grid, x2o_o, ind, &
                             sw_decomp, c1, c2, c3, c4)

  type(ice_ocean_boundary_type), intent(inout) :: ice_ocean_boundary !< A type for the ice ocean boundary
  type(ocean_grid_type),  intent(in)    :: grid      !<
  !type(mct_aVect),       intent(in)    :: x2o_o
  real(kind=8),           intent(in)    :: x2o_o(:,:)
  type(cpl_indices),      intent(inout) :: ind
  logical,                intent(in)    :: sw_decomp !< controls if shortwave is
                                                     !!decomposed into four components
  real(kind=8),           intent(in), optional :: c1, c2, c3, c4 !< Coeffs. used in the shortwave decomposition

  ! local variables
  integer               :: i, j, k, ig, jg !< grid indices

  ! /TODO The following comments summarizes the mismatches between MCT and MOM6 in terms
  ! of ice ocean fluxes.

  ! Redundancies:
  ! x2o_Foxx_lat, x2o_Faxa_prec, x2o_Faxa_evap - **latent from MOM6 and coupler are different**
  ! x2o_Faxa_prec = x2o_Faxa_rain + x2o_Faxa_snow

  ! Variables that **could not** be verified so far:
  ! x2o_Foxx_rofl
  ! x2o_Foxx_rof

  ! Variables in IOB that are **NOT** filled by the coupler:
  ! ustar_berg, frictional velocity beneath icebergs (m/s)
  ! area_berg, area covered by icebergs(m2/m2)
  ! mass_berg, mass of icebergs(kg/m2)
  ! runoff_hflx, heat content of liquid runoff (W/m2)
  ! calving_hflx, heat content of frozen runoff (W/m2)
  ! mi, mass of ice (kg/m2)

  ! Variables in the coupler that are **NOT** used in MOM6 (i.e., no corresponding field in IOB):
  ! x2o_Fioi_melth, heat flux from snow & ice melt (W/m2)
  ! x2o_Fioi_meltw, snow melt flux (kg/m2/s)
  ! x2o_Si_ifrac, fractional ice wrt ocean
  ! x2o_So_duu10n, 10m wind speed squared (m^2/s^2)
  ! x2o_Sa_co2prog, bottom atm level prognostic CO2
  ! x2o_Sa_co2diag, bottom atm level diagnostic CO2
  ! x2o_Foxx_lat, latent heat flux  (W/m2)

  ! Langmuir related fields:
  ! surface Stokes drift, x-comp. (x2o_Sw_ustokes)
  ! surface Stokes drift, y-comp. (x2o_Sw_vstokes)
  ! wave model langmuir multiplier (x2o_Sw_lamult)

  ! Biogeochemistry:
  ! x2o_Fioi_bcpho, Black Carbon hydrophobic release from sea ice component
  ! x2o_Fioi_bcphi, Black Carbon hydrophilic release from sea ice component
  ! x2o_Fioi_flxdst, Dust release from sea ice component
  ! x2o_Faxa_bcphidry, Black Carbon hydrophilic dry deposition
  ! x2o_Faxa_bcphodry, Black Carbon hydrophobic dry deposition
  ! x2o_Faxa_bcphiwet, Black Carbon hydrophobic wet deposition
  ! x2o_Faxa_ocphidry, Organic Carbon hydrophilic dry deposition
  ! x2o_Faxa_ocphodry, Organic Carbon hydrophobic dry deposition
  ! x2o_Faxa_ocphiwet, Organic Carbon hydrophilic dry deposition
  ! x2o_Faxa_dstwet, Sizes 1 to 4 dust - wet deposition
  ! x2o_Faxa_dstdry, Sizes 1 to 4 dust - dry depositi

  k = 0
  do j = grid%jsc, grid%jec
    jg = j + grid%jdg_offset
    do i = grid%isc, grid%iec
      k = k + 1 ! Increment position within gindex
      ig = i + grid%idg_offset
      ! sea-level pressure (Pa)
      ice_ocean_boundary%p(ig,jg) = x2o_o(ind%x2o_Sa_pslv,k)
      ! zonal wind stress (taux)
      ice_ocean_boundary%u_flux(ig,jg) = x2o_o(ind%x2o_Foxx_taux,k)
      ! meridional wind stress (tauy)
      ice_ocean_boundary%v_flux(ig,jg) = x2o_o(ind%x2o_Foxx_tauy,k)
      ! sensible heat flux (W/m2)
      ice_ocean_boundary%t_flux(ig,jg) = -x2o_o(ind%x2o_Foxx_sen,k)
      ! salt flux
      ice_ocean_boundary%salt_flux(ig,jg) = -x2o_o(ind%x2o_Fioi_salt,k)
      ! heat content from frozen runoff
      ice_ocean_boundary%calving_hflx(ig,jg) = 0.0
      ! snow melt flux
      !ice_ocean_boundary%fprec(ig,jg) = x2o_o(ind%x2o_Fioi_meltw,k)
      ! river runoff flux
      ice_ocean_boundary%runoff(ig,jg) = x2o_o(ind%x2o_Foxx_rofl,k)
      ! ice runoff flux
      ice_ocean_boundary%calving(ig,jg) = x2o_o(ind%x2o_Foxx_rofi,k)
      ! liquid precipitation (rain)
      ice_ocean_boundary%lprec(ig,jg) = x2o_o(ind%x2o_Faxa_rain,k)
      ! frozen precipitation (snow)
      ice_ocean_boundary%fprec(ig,jg) = x2o_o(ind%x2o_Faxa_snow,k)
      ! evaporation flux, MOM6 calls q_flux specific humidity (kg/m2/s)
      ice_ocean_boundary%q_flux(ig,jg) = -x2o_o(ind%x2o_Foxx_evap,k)
      ! longwave radiation, sum up and down (W/m2)
      ice_ocean_boundary%lw_flux(ig,jg) = x2o_o(ind%x2o_Faxa_lwdn,k) + x2o_o(ind%x2o_Foxx_lwup,k)
      if (sw_decomp) then
        ! Use runtime coefficients to decompose net short-wave heat flux into 4 components
        ! 1) visible, direct shortwave (W/m2)
        ice_ocean_boundary%sw_flux_vis_dir(ig,jg) = x2o_o(ind%x2o_Foxx_swnet,k)*c1
        ! 2) visible, diffuse shortwave (W/m2)
        ice_ocean_boundary%sw_flux_vis_dif(ig,jg) = x2o_o(ind%x2o_Foxx_swnet,k)*c2
        ! 3) near-IR, direct shortwave (W/m2)
        ice_ocean_boundary%sw_flux_nir_dir(ig,jg) = x2o_o(ind%x2o_Foxx_swnet,k)*c3
        ! 4) near-IR, diffuse shortwave (W/m2)
        ice_ocean_boundary%sw_flux_nir_dif(ig,jg) = x2o_o(ind%x2o_Foxx_swnet,k)*c4
      else
        call MOM_error(FATAL,"fill_ice_ocean_bnd: this option has not been implemented yet."// &
                       "Shortwave must be decomposed using coeffs. c1, c2, c3, c4.");
      endif
    enddo
  enddo

  ice_ocean_boundary%wind_stagger = AGRID

end subroutine fill_ice_ocean_bnd

!> Step forward ocean model for coupling interval
subroutine ocn_run_mct( EClock, cdata_o, x2o_o, o2x_o)
  type(ESMF_Clock), intent(inout) :: EClock  !< Time and time step ? \todo Why must this be intent(inout)?
  type(seq_cdata),  intent(inout) :: cdata_o !< Input parameters
  type(mct_aVect),  intent(inout) :: x2o_o   !< Fluxes from coupler to ocean, computed by ocean
  type(mct_aVect),  intent(inout) :: o2x_o   !< Fluxes from ocean to coupler, computed by ocean
  ! Local variables
  type(ESMF_time) :: current_time             !< Time to be reached at the end of ocean cpl interval
  type(ESMF_time) :: time_start_ESMF          !< Time at the start of the coupling interval
  type(ESMF_timeInterval) :: ocn_cpl_interval !< The length of one ocean coupling interval
  integer :: year, month, day, hour, minute, seconds, seconds_n, seconds_d, rc
  logical :: write_restart_at_eod
  type(time_type) :: time_start        !< Start of coupled time interval to pass to MOM6
  type(time_type) :: coupling_timestep !< Coupled time interval to pass to MOM6
  character(len=128) :: err_msg
  character(len=32)  :: timestamp      !< Name of intermediate restart file

  ! Compute the time at the start of this coupling interval
  call ESMF_ClockGet(EClock, PrevTime=time_start_ESMF, rc=rc)
  call ESMF_TimeGet(time_start_ESMF, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
  time_start = set_date(year, month, day, hour, minute, seconds, err_msg=err_msg)

  ! Debugging clocks
  if (debug .and. is_root_pe()) then
    call ESMF_ClockGet(EClock, CurrTime=current_time, rc=rc)
    call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
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
    call ESMF_ClockGet(EClock, TimeStep=ocn_cpl_interval, rc=rc)
    call ESMF_TimeIntervalGet(ocn_cpl_interval, yy=year, mm=month, d=day, s=seconds, sn=seconds_n, sd=seconds_d, rc=rc)
    write(6,*) 'ocn_init_mct, time step: y,m,d-',year,month,day,'s,sn,sd=',seconds,seconds_n,seconds_d
  endif

  ! Translate the coupling time interval
  call ESMF_ClockGet(EClock, TimeStep=ocn_cpl_interval, rc=rc)
  call ESMF_TimeIntervalGet(ocn_cpl_interval, yy=year, mm=month, d=day, s=seconds, sn=seconds_n, sd=seconds_d, rc=rc)
  coupling_timestep = set_time(seconds, days=day, err_msg=err_msg)

  ! set (actually, get from mct) the cdata pointers:
  ! \todo this was done in _init_, is it needed again. Does this infodata need to be in glb%?
  call seq_cdata_setptrs(cdata_o, infodata=glb%infodata)

  ! fill ice ocean boundary
  call fill_ice_ocean_bnd(glb%ice_ocean_boundary, glb%grid, x2o_o%rattr, glb%ind, glb%sw_decomp, &
                          glb%c1, glb%c2, glb%c3, glb%c4)
  if (debug .and. is_root_pe()) write(6,*) 'fill_ice_ocean_bnd'

  call update_ocean_model(glb%ice_ocean_boundary, glb%ocn_state, glb%ocn_public, &
                          time_start, coupling_timestep)

  !--- write out intermediate restart file when needed.
  ! Check alarms for flag to write restart at end of day
  write_restart_at_eod = seq_timemgr_RestartAlarmIsOn(EClock)
  if (debug .and. is_root_pe()) write(6,*) 'ocn_run_mct, write_restart_at_eod=', write_restart_at_eod

  if (write_restart_at_eod) then
    !timestamp = date_to_string(EClock)
    ! \todo add time stamp to ocean_model_restart
    !call ocean_model_restart(glb%ocn_state, timestamp)
    call ocean_model_restart(glb%ocn_state)
  endif

end subroutine ocn_run_mct

!> Finalizes MOM6
!!
!! \todo This needs to be done here.
subroutine ocn_final_mct( EClock, cdata_o, x2o_o, o2x_o)
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o  !< Fluxes from coupler to ocean, computed by ocean
    type(mct_aVect)             , intent(inout) :: o2x_o  !< Fluxes from ocean to coupler, computed by ocean

end subroutine ocn_final_mct


!> Sets mct global segment maps for the MOM decomposition.
!!
!! \todo Find out if we should only provide indirect indexing for ocean points and not land.
subroutine ocn_SetGSMap_mct(mpicom_ocn, MOM_MCT_ID, gsMap_ocn, gsMap3d_ocn)
  integer,         intent(in)    :: mpicom_ocn  !< MPI communicator
  integer,         intent(in)    :: MOM_MCT_ID  !< MCT component ID
  type(mct_gsMap), intent(inout) :: gsMap_ocn   !< MCT global segment map for 2d data
  type(mct_gsMap), intent(inout) :: gsMap3d_ocn !< MCT global segment map for 3d data
  ! Local variables
  integer :: lsize   !< Local size of indirect indexing array
  integer :: i, j, k !< Local indices
  integer :: ni, nj  !< Declared sizes of h-point arrays
  integer :: ig, jg  !< Global indices
  type(ocean_grid_type), pointer :: grid => NULL() !< A pointer to a grid structure
  integer, allocatable           :: gindex(:)      !< Indirect indices

  grid => glb%grid ! for convenience
  if (.not. associated(grid)) call MOM_error(FATAL, 'ocn_comp_mct.F90, ocn_SetGSMap_mct():' // &
      'grid returned from get_state_pointers() was not associated!')

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

end subroutine ocn_domain_mct

end module ocn_comp_mct
