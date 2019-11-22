!> Interfaces for MOM6 ensembles and data assimilation.
module MOM_oda_driver_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the top-level module for MOM6 ocean data assimilation.
! It can be used to gather an ensemble of ocean states
! before calling ensemble filter routines which calculate
! increments based on cross-ensemble co-variance. It can also
! be used to compare gridded model state variables to in-situ
! observations without applying DA incrementa.
!
! init_oda:  Initialize the ODA module
! set_analysis_time : update time for performing next analysis
! set_prior: Store prior model state
! oda: call to filter
! get_posterior : returns posterior increments (or full state) for the current ensemble member
!
! Authors: Matthew.Harrison@noaa.gov
!          Feiyu.Liu@noaa.gov and
!          Tony.Rosati@noaa.gov
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This file is part of MOM6. see LICENSE.md for the license.
use fms_mod, only : open_namelist_file, close_file, check_nml_error
use fms_mod, only : error_mesg, FATAL
use mpp_mod, only : stdout, stdlog, mpp_error, npes=>mpp_npes,pe=>mpp_pe
use mpp_mod, only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : set_current_pelist => mpp_set_current_pelist
use mpp_mod, only : set_root_pe => mpp_set_root_pe
use mpp_mod, only : mpp_sync_self, mpp_sum, get_pelist=>mpp_get_current_pelist, mpp_root_pe
use mpp_mod, only : set_stack_size=>mpp_set_stack_size, broadcast=>mpp_broadcast
use mpp_io_mod, only : io_set_stack_size=>mpp_io_set_stack_size
use mpp_io_mod, only : MPP_SINGLE,MPP_MULTI
use mpp_domains_mod, only : domain2d, mpp_global_field, mpp_update_domains
use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain
use mpp_domains_mod, only : mpp_redistribute, mpp_broadcast_domain
use mpp_domains_mod, only : set_domains_stack_size=>mpp_domains_set_stack_size
use MOM_diag_mediator, only : register_diag_field, post_data, diag_update_remap_grids
use MOM_diag_mediator, only : diag_ctrl, enable_averaging, disable_averaging
use MOM_time_manager, only : init_external_field, get_external_field_size, time_interp_external_init
use MOM_time_manager, only : time_type, decrement_time, increment_time
use MOM_time_manager, only : get_date, operator(>=),operator(/=),operator(==),operator(<)
use MOM_horizontal_regridding, only : horiz_interp_and_extrap_tracer, myStats
use ensemble_manager_mod, only : get_ensemble_id, get_ensemble_size
use ensemble_manager_mod, only : get_ensemble_pelist, get_ensemble_filter_pelist
use MOM_constants, only : seconds_per_hour
! ODA Modules
use ocean_da_types_mod, only : grid_type, ocean_profile_type, ocean_control_struct
use ocean_da_core_mod, only : ocean_da_core_init, get_profiles, kd_root
use ocean_da_types_mod, only : TEMP_ID, SALT_ID
use ocean_da_types_mod, only : ODA_PFL, ODA_XBT, ODA_MRB, ODA_OISST
#ifdef ENABLE_FILTER
use eakf_oda_mod, only : ensemble_filter
#endif
use write_ocean_obs_mod, only : open_profile_file
use write_ocean_obs_mod, only : write_profile,close_profile_file
#ifdef ENABLE_FILTER
use kdtree, only : kd_root !# JEDI
#endif
! MOM Modules
use MOM_io, only : slasher, MOM_read_data
use MOM_diag_mediator, only : diag_ctrl, set_axes_info
use MOM_error_handler, only : FATAL, WARNING, MOM_error, MOM_mesg, is_root_pe
use MOM_get_input, only : get_MOM_input, directories
use MOM_grid, only : ocean_grid_type, MOM_grid_init
use MOM_grid_initialize, only : set_grid_metrics
use MOM_hor_index, only : hor_index_type, hor_index_init
use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_transcribe_grid, only : copy_dyngrid_to_MOM_grid, copy_MOM_grid_to_dyngrid
use MOM_fixed_initialization, only : MOM_initialize_fixed, MOM_initialize_topography
use MOM_coord_initialization, only : MOM_initialize_coord
use MOM_file_parser, only : read_param, get_param, param_file_type
use MOM_string_functions, only : lowercase
use MOM_ALE, only : ALE_CS, ALE_initThicknessToCoord, ALE_init, ALE_updateVerticalGridType
use MOM_domains, only : MOM_domains_init, MOM_domain_type, clone_MOM_domain
use MOM_remapping, only : remapping_CS, initialize_remapping, remapping_core_h
use MOM_regridding, only : regridding_CS, initialize_regridding
use MOM_regridding, only : regridding_main, set_regrid_params
use MOM_unit_scaling, only : unit_scale_type, unit_scaling_init
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type, verticalGridInit

implicit none ; private

public :: init_oda, oda_end, set_prior_tracer, get_posterior_tracer
public :: set_analysis_time, oda, save_obs_diff, apply_oda_tracer_increments

#include <MOM_memory.h>

!> Control structure that contains tracer ids for temperature and salinity
type :: INC_CS
   integer :: fldno = 0 !< integer field identification
   integer :: T_id !< integer id for temperature
   integer :: S_id !< integer id for salinity
end type INC_CS

!> Control structure that contains a transpose of the ocean state across ensemble members.
type, public :: ODA_CS ; private
  type(ocean_control_struct), pointer :: Ocean_prior=> NULL()     !< ensemble ocean prior states in DA space
  type(ocean_control_struct), pointer :: Ocean_posterior=> NULL() !< ensemble ocean posterior states
                                                                  !! or increments to prior in DA space
  type(ocean_control_struct), pointer :: Ocean_increment=> NULL() !< ensemble increments in DA space
  integer :: nk                                                   !< number of vertical layers used for DA
  type(ocean_grid_type), pointer :: Grid => NULL()                !< MOM6 grid type and decomposition for the DA
  type(ocean_grid_type), pointer :: G => NULL()                   !< MOM6 grid type and decomposition for the model
  type(ptr_mpp_domain), pointer, dimension(:) :: domains => NULL()!< Pointer to mpp_domain objects
                                                                  !! for ensemble members
  type(verticalGrid_type), pointer :: GV => NULL() !< vertical grid for DA
  type(unit_scale_type), pointer :: &
    US => NULL()    !< structure containing various unit conversion factors for DA
  type(domain2d), pointer :: mpp_domain => NULL()   !< Pointer to a mpp domain object for DA
  type(grid_type), pointer :: oda_grid              !< local tracer grid
  real, pointer, dimension(:,:,:) :: h => NULL()    !<layer thicknesses [H ~> m or kg m-2] for DA
  type(thermo_var_ptrs), pointer :: tv => NULL()    !< pointer to thermodynamic variables
  type(thermo_var_ptrs), pointer :: tv_bc => NULL() !< pointer for thermodynamic variable bias adjustment
  integer :: ni          !< global i-direction grid size
  integer :: nj          !< global j-direction grid size
  logical :: reentrant_x !< grid is reentrant in the x direction
  logical :: reentrant_y !< grid is reentrant in the y direction
  logical :: tripolar_N !< grid is folded at its north edge
  logical :: symmetric !< Values at C-grid locations are symmetric
  integer :: assim_method !< Method: NO_FILTER or EAKF
  integer :: ensemble_size !< Size of the ensemble
  integer :: ensemble_id = 0 !< id of the current ensemble member
  integer, pointer, dimension(:,:) :: ensemble_pelist !< PE list for ensemble members
  integer, pointer, dimension(:) :: filter_pelist !< PE list for ensemble members
  integer :: assim_frequency !< analysis interval in hours
  ! Profiles local to the analysis domain
  type(ocean_profile_type), pointer :: Profiles => NULL() !< pointer to linked list of all available profiles
  type(ocean_profile_type), pointer :: CProfiles => NULL()!< pointer to linked list of current profiles
  type(kd_root), pointer :: kdroot => NULL() !< A structure for storing nearest neighbors
  type(ALE_CS), pointer :: ALE_CS=>NULL() !< ALE control structure for DA
  logical :: use_ALE_algorithm !< true is using ALE remapping
  type(regridding_CS) :: regridCS !< ALE control structure for regridding
  type(remapping_CS) :: remapCS !< ALE control structure for remapping
  type(time_type) :: Time !< Current Analysis time
  type(diag_ctrl), pointer :: diag_cs => NULL() !<Diagnostics control structure
  logical :: do_bias_correction !< setting this to true will enable bias adjustment
  real :: correction_multiplier !< non-dimensional factor for bias adjustment [nodim]
  logical :: write_obs !< setting this to true enables writing of profile misfits to file
  integer :: id_inc_t, id_inc_s !< integer ids for temperature and salinity increments
  type(INC_CS) :: INC_CS !< control structure containing field indices
end type ODA_CS

!> A structure with a pointer to a domain2d, to allow for the creation of arrays of pointers.
type :: ptr_mpp_domain
  type(domain2d), pointer :: mpp_domain => NULL() !< pointer to an mpp domain2d
end type ptr_mpp_domain


  integer, parameter :: NO_FILTER = 0 !< No filter adjustments
  integer, parameter :: EAKF = 1 !< Use interfaces to GFDL ensemble adjustment filter
  integer :: id_clock_oda_init !< timer id for oda_init
  integer :: id_clock_oda_prior !< timer id for oda_prior
  integer :: id_clock_oda_filter !< timer id for oda_filter
  integer :: id_clock_oda_posterior !< timer id for oda_get_posterior
  integer :: id_clock_bias_correction !< timer id for bias correction
  integer :: id_clock_apply_increments !< timer id for oda_apply increments
  integer :: temp_fid, salt_fid !< profile file handles for temperature and salinity

contains
!! initialize analysis grid and ODA-related variables
!! information for all ensemble members
  subroutine init_oda(Time, G, GV, US, diag_CS, CS)

    type(time_type), intent(in) :: Time !< The current model time.
    type(ocean_grid_type), pointer :: G !< domain and grid information for ocean model
    type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
    type(unit_scale_type), pointer :: US !< unit scaling type
    type(diag_ctrl), target, intent(inout) :: diag_CS !< diagnostic control structure
    type(ODA_CS), pointer :: CS                       !< ocean DA control structure

 ! Local variables
    type(thermo_var_ptrs) :: tv_dummy
    type(dyn_horgrid_type), pointer :: dG=> NULL()
    type(hor_index_type), pointer :: HI=> NULL()
    type(directories) :: dirs
    type(grid_type), pointer :: T_grid => NULL() !< global tracer grid
    type(param_file_type) :: PF
    integer :: n, m, k, i, j, nk
    integer :: is,ie,js,je,isd,ied,jsd,jed
    integer :: stdout_unit, ioun, ierr, io_status
    character(len=32) :: assim_method
    integer :: npes_pm, ens_info(6), ni, nj
    character(len=128) :: mesg
    character(len=32) :: fldnam
    character(len=30) :: coord_mode
    character(len=200) :: inputdir, basin_file
    logical :: reentrant_x, reentrant_y, tripolar_N, symmetric
    character(len=80) :: remap_scheme
    character(len=80) :: inc_file
    integer, dimension(4) :: fld_sz
    real :: missing_value

    !---- namelist with default values
    logical :: write_obs = .false.
    character(len=80) :: obs_file
    logical :: do_bias_correction = .false.
    character(len=80) :: bias_correction_file
    real :: correction_multiplier = 1.0
    namelist /oda_init_nml/ write_obs, obs_file, do_bias_correction, bias_correction_file, correction_multiplier

    if (associated(CS)) call mpp_error(FATAL,'Calling oda_init with associated control structure')
    allocate(CS)

    ioun = open_namelist_file()
    read(UNIT=ioun, NML=oda_init_nml, IOSTAT=io_status)
    ierr = check_nml_error(io_status,'oda_init_nml')
    call close_file(ioun)
    CS%do_bias_correction = do_bias_correction
    CS%correction_multiplier = correction_multiplier
    CS%write_obs = write_obs

    !! Use ens0 parameters, which is set up solely for the analysis grid
    call get_MOM_input(PF,dirs,ensemble_num=0)
    call get_param(PF, "MOM", "ASSIM_METHOD", assim_method,  &
            "String which determines the data assimilation method" // &
            "Valid methods are: \'EAKF\' and \'NO_ASSIM\'", default='NO_ASSIM')
    call get_param(PF, "MOM", "ASSIM_FREQUENCY", CS%assim_frequency,  &
            "data assimilation frequency in hours")
    call get_param(PF, "MOM", "USE_REGRIDDING", CS%use_ALE_algorithm , &
            "If True, use the ALE algorithm (regridding/remapping).\n"//&
            "If False, use the layered isopycnal algorithm.", default=.false. )
    call get_param(PF, "MOM", "REENTRANT_X", CS%reentrant_x, &
            "If true, the domain is zonally reentrant.", default=.true.)
    call get_param(PF, "MOM", "REENTRANT_Y", CS%reentrant_y, &
            "If true, the domain is meridionally reentrant.", default=.false.)
    call get_param(PF,"MOM", "TRIPOLAR_N", CS%tripolar_N, &
            "Use tripolar connectivity at the northern edge of the \n"//&
            "domain.  With TRIPOLAR_N, NIGLOBAL must be even.", default=.false.)
    call get_param(PF,"MOM", "NIGLOBAL", CS%ni, &
            "The total number of thickness grid points in the \n"//&
            "x-direction in the physical domain.")
    call get_param(PF,"MOM", "NJGLOBAL", CS%nj, &
            "The total number of thickness grid points in the \n"//&
            "y-direction in the physical domain.")
    call get_param(PF, 'MOM', "INPUTDIR", inputdir)
    call get_param(PF, "MOM", "REMAPPING_SCHEME", remap_scheme, default="PPM_H4")
    inputdir = slasher(inputdir)

    select case(lowercase(trim(assim_method)))
    case('eakf')
        CS%assim_method = EAKF
    case('none')
        CS%assim_method = NO_FILTER
    case default
        call mpp_error(FATAL,'Invalid assimilation method provided')
    end select

    ens_info = get_ensemble_size()
    CS%ensemble_size = ens_info(1)
    npes_pm=ens_info(3)
    CS%ensemble_id = get_ensemble_id()
    allocate(CS%ensemble_pelist(CS%ensemble_size,npes_pm))
    allocate(CS%filter_pelist(CS%ensemble_size*npes_pm))
    call get_ensemble_pelist(CS%ensemble_pelist,'ocean')
    call get_ensemble_filter_pelist(CS%filter_pelist,'ocean')
    id_clock_apply_increments = mpp_clock_id('(ODA applying increments)')
    id_clock_bias_correction = mpp_clock_id('(Bias correction through increments)')
    !! Switch to global ocean pelist
    call set_current_pelist(CS%filter_pelist)
    if(is_root_pe()) print *, 'Initialize ODA'

    id_clock_oda_init = mpp_clock_id('(ODA initialization)')
    id_clock_oda_prior = mpp_clock_id('(ODA setting prior)')
    id_clock_oda_filter = mpp_clock_id('(ODA filter computation)')
    id_clock_oda_posterior = mpp_clock_id('(ODA getting posterior)')
    call mpp_clock_begin(id_clock_oda_init)

    !! set up and broadcast ensemble domains to enable redistribution later
    allocate(CS%domains(CS%ensemble_size))
    CS%domains(CS%ensemble_id)%mpp_domain => G%Domain%mpp_domain
    do n=1,CS%ensemble_size
      if(.not. associated(CS%domains(n)%mpp_domain)) allocate(CS%domains(n)%mpp_domain)
      call mpp_broadcast_domain(CS%domains(n)%mpp_domain)
    enddo

    CS%G => G
    allocate(CS%Grid)
    !! params NIHALO_ODA, NJHALO_ODA set the DA halo size
    call MOM_domains_init(CS%Grid%Domain,PF,param_suffix='_ODA')
    allocate(HI)
    call hor_index_init(CS%Grid%Domain, HI, PF)
    call verticalGridInit( PF, CS%GV, US )
    allocate(dG)
    call create_dyn_horgrid(dG,HI)
    call clone_MOM_domain(CS%Grid%Domain, dG%Domain,symmetric=.false.)
    call set_grid_metrics(dG,PF)
    call MOM_initialize_topography(dg%bathyT,dG%max_depth,dG,PF,US)
    call MOM_initialize_coord(CS%GV, US, PF, .false., &
         dirs%output_directory, tv_dummy, dG%max_depth)
    call ALE_init(PF, CS%GV, US, dG%max_depth, CS%ALE_CS)
    call MOM_grid_init(CS%Grid, PF)
    call ALE_updateVerticalGridType(CS%ALE_CS,CS%GV)
    call copy_dyngrid_to_MOM_grid(dG, CS%Grid, US)
    CS%mpp_domain => CS%Grid%Domain%mpp_domain
    CS%Grid%ke = CS%GV%ke
    CS%nk = CS%GV%ke

    ! initialize storage for prior and posterior
    allocate(CS%Ocean_prior)
    call init_ocean_ensemble(CS%Ocean_prior,CS%Grid,CS%GV,CS%ensemble_size)
    allocate(CS%Ocean_posterior)
    call init_ocean_ensemble(CS%Ocean_posterior,CS%Grid,CS%GV,CS%ensemble_size)
    allocate(CS%Ocean_increment)
    call init_ocean_ensemble(CS%Ocean_increment,CS%Grid,CS%GV,CS%ensemble_size)

    call get_param(PF, 'oda_driver', "REGRIDDING_COORDINATE_MODE", coord_mode, &
         "Coordinate mode for vertical regridding.", &
         default="ZSTAR", fail_if_missing=.false.)
    call initialize_regridding(CS%regridCS, CS%GV, US, dG%max_depth,PF,'oda_driver',coord_mode,'','')
    call initialize_remapping(CS%remapCS,remap_scheme)
    call set_regrid_params(CS%regridCS, min_thickness=0.)

    ! get domain indices from model domain decomposition
    isd = G%isd; ied = G%ied; jsd = G%jsd; jed = G%jed
    if(.not. associated(CS%h)) then
        allocate(CS%h(isd:ied,jsd:jed,CS%GV%ke)); CS%h(:,:,:) = CS%GV%Angstrom_m
        ! assign thicknesses
        call ALE_initThicknessToCoord(CS%ALE_CS,G,CS%GV,CS%h)
    endif
    allocate(CS%tv)     ! storage for increment
    ! increments are stored in z* and model domain decomposition
    allocate(CS%tv%T(isd:ied,jsd:jed,CS%GV%ke)); CS%tv%T(:,:,:)=0.0
    allocate(CS%tv%S(isd:ied,jsd:jed,CS%GV%ke)); CS%tv%S(:,:,:)=0.0

    ! get domain indices from analysis domain decomposition
    isd = CS%Grid%isd; ied = CS%Grid%ied; jsd = CS%Grid%jsd; jed = CS%Grid%jed
    allocate(CS%oda_grid) ! local grid information needed from analysis
    CS%oda_grid%x => CS%Grid%geolonT
    CS%oda_grid%y => CS%Grid%geolatT

    ! get basin flag from file
    call get_param(PF, 'oda_driver', "BASIN_FILE", basin_file, &
            "A file in which to find the basin masks, in variable 'basin'.", &
            default="basin.nc")
    basin_file = trim(inputdir) // trim(basin_file)
    allocate(CS%oda_grid%basin_mask(isd:ied,jsd:jed))
    CS%oda_grid%basin_mask(:,:) = 0.0
    call MOM_read_data(basin_file,'basin',CS%oda_grid%basin_mask,CS%Grid%domain, timelevel=1)

    ! set up diag variables for analysis increments
    CS%diag_cs => diag_CS
    CS%id_inc_t=register_diag_field('ocean_model','temp_increment',diag_CS%axesTL,&
            Time,'ocean potential temperature increments','degC')
    CS%id_inc_s=register_diag_field('ocean_model','salt_increment',diag_CS%axesTL,&
            Time,'ocean salinity increments','psu')

    !!  get global grid information from ocean model needed for ODA initialization
    call set_up_global_tgrid(T_grid, CS, G)

    call ocean_da_core_init(CS%mpp_domain, T_grid, CS%Profiles, Time)

    !! Set the initial assimilation time
    CS%Time=Time
    !CS%Time=increment_time(Time,CS%assim_frequency*seconds_per_hour)

    deallocate(T_grid)
    !deallocate(h)
    call mpp_clock_end(id_clock_oda_init)
    !! switch back to ensemble member pelist
    call set_current_pelist(CS%ensemble_pelist(CS%ensemble_id,:))

    if (CS%do_bias_correction) then
      call time_interp_external_init()
      inc_file = trim(inputdir) // trim(bias_correction_file)
      CS%INC_CS%T_id = init_external_field(inc_file, "temp_increment", &
              correct_leap_year_inconsistency=.true.,verbose=.true.,domain=G%Domain%mpp_domain)
      CS%INC_CS%S_id = init_external_field(inc_file, "salt_increment", &
              correct_leap_year_inconsistency=.true.,verbose=.true.,domain=G%Domain%mpp_domain)
      fld_sz = get_external_field_size(CS%INC_CS%T_id)
      CS%INC_CS%fldno = 2
      if (CS%nk .ne. fld_sz(3)) call mpp_error(FATAL,'Increment levels /= ODA levels')

      isd = G%isd; ied = G%ied; jsd = G%jsd; jed = G%jed
      allocate(CS%tv_bc)     ! storage for increment
      allocate(CS%tv_bc%T(isd:ied,jsd:jed,CS%GV%ke)); CS%tv_bc%T(:,:,:)=0.0
      allocate(CS%tv_bc%S(isd:ied,jsd:jed,CS%GV%ke)); CS%tv_bc%S(:,:,:)=0.0
    endif

    if (CS%write_obs) then
       temp_fid = open_profile_file("temp_"//trim(obs_file))
       salt_fid = open_profile_file("salt_"//trim(obs_file))
    end if

end subroutine init_oda

!> Copy ensemble member tracers to ensemble vector.
subroutine set_prior_tracer(Time, G, GV, h, tv, CS)
    type(time_type), intent(in)    :: Time                       !< The current model time
    type(ocean_grid_type), pointer :: G                          !< domain and grid information for ocean model
    type(verticalGrid_type),               intent(in)    :: GV   !< The ocean's vertical grid structure
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h    !< Layer thicknesses, in H (usually m or kg m-2)
    type(thermo_var_ptrs),                 intent(in) :: tv      !< A structure pointing to various thermodynamic variables
    type(ODA_CS), pointer :: CS                                  !< ocean DA control structure

    integer :: i, j, m, n, ss
    integer :: isc, iec, jsc, jec
    logical :: used
    real, dimension(SZI_(G),SZJ_(G),SZK_(CS%Grid)) :: T
    real, dimension(SZI_(G),SZJ_(G),SZK_(CS%Grid)) :: S

    !! return if not time for analysis
    if (Time < CS%Time .or. .not. associated(CS) .or. CS%assim_method .eq. NO_FILTER) return

    if (.not. ASSOCIATED(CS%Grid)) call MOM_ERROR(FATAL,'ODA_CS ensemble horizontal grid not associated')
    if (.not. ASSOCIATED(CS%GV)) call MOM_ERROR(FATAL,'ODA_CS ensemble vertical grid not associated')

    !! switch to global pelist
    call set_current_pelist(CS%filter_pelist)
    call mpp_clock_begin(id_clock_oda_prior)
    if(is_root_pe()) print *, 'Setting prior'

    T = 0.0; S = 0.0
    isc=G%isc; iec=G%iec; jsc=G%jsc; jec=G%jec
    do j=jsc,jec; do i=isc,iec
      call remapping_core_h(CS%remapCS, GV%ke, h(i,j,:), tv%T(i,j,:), &
           CS%nk, CS%h(i,j,:), T(i,j,:))
      call remapping_core_h(CS%remapCS, GV%ke, h(i,j,:), tv%S(i,j,:), &
           CS%nk, CS%h(i,j,:), S(i,j,:))
    enddo; enddo

    do m=1,CS%ensemble_size
      call mpp_redistribute(CS%domains(m)%mpp_domain, T,&
           CS%mpp_domain, CS%Ocean_prior%T(:,:,:,m), complete=.true.)
      call mpp_redistribute(CS%domains(m)%mpp_domain, S,&
           CS%mpp_domain, CS%Ocean_prior%S(:,:,:,m), complete=.true.)
    enddo

    do m=1,CS%ensemble_size
      call mpp_update_domains(CS%Ocean_prior%T(:,:,:,m), CS%mpp_domain)
      call mpp_update_domains(CS%Ocean_prior%S(:,:,:,m), CS%mpp_domain)
    enddo

    call mpp_clock_end(id_clock_oda_prior)
    !! switch back to ensemble member pelist
    call set_current_pelist(CS%ensemble_pelist(CS%ensemble_id,:))

    return

end subroutine set_prior_tracer

!> Returns posterior adjustments or full state
!!Note that only those PEs associated with an ensemble member receive data
subroutine get_posterior_tracer(Time, CS, h, tv, increment)
    type(time_type), intent(in) :: Time               !< the current model time
    type(ODA_CS), pointer :: CS                       !< ocean DA control structure
    real, dimension(:,:,:), pointer, optional :: h    !< Layer thicknesses, in H (usually m or kg m-2)
    type(thermo_var_ptrs), pointer, optional :: tv    !< A structure pointing to various thermodynamic variables
    logical, optional, intent(in) :: increment        !< logical flag to determine if tracer increment or full value are returned

    integer :: i, j, m
    logical :: used, get_inc

    !! return if not analysis time (retain pointers for h and tv)
    if (Time < CS%Time .or. CS%assim_method .eq. NO_FILTER) return

    !! switch to global pelist
    call set_current_pelist(CS%filter_pelist)

    call mpp_clock_begin(id_clock_oda_posterior)
    if(is_root_pe()) print *, 'Getting posterior'
    if( present(h) .and. .not. associated(h) ) h => CS%h ! Get analysis thickness

    !! Calculate and redistribute increments to CS%tv right after assimilation
    !! Retain CS%tv to calculate increments for IAU updates CS%tv_inc otherwise
    get_inc = .true.
    if(present(increment)) get_inc = increment

    if(get_inc) then
      CS%Ocean_increment%T = CS%Ocean_posterior%T - CS%Ocean_prior%T
      CS%Ocean_increment%S = CS%Ocean_posterior%S - CS%Ocean_prior%S
    endif

    do m=1,CS%ensemble_size
      if(get_inc) then
        call mpp_redistribute(CS%mpp_domain, CS%Ocean_increment%T(:,:,:,m), &
                CS%domains(m)%mpp_domain, CS%tv%T, complete=.true.)
        call mpp_redistribute(CS%mpp_domain, CS%Ocean_increment%S(:,:,:,m), &
                CS%domains(m)%mpp_domain, CS%tv%S, complete=.true.)
      else
        call mpp_redistribute(CS%mpp_domain, CS%Ocean_posterior%T(:,:,:,m), &
                CS%domains(m)%mpp_domain, CS%tv%T, complete=.true.)
        call mpp_redistribute(CS%mpp_domain, CS%Ocean_posterior%S(:,:,:,m), &
                CS%domains(m)%mpp_domain, CS%tv%S, complete=.true.)
        if(present(tv)) tv => CS%tv
      endif
    end do

    call mpp_clock_end(id_clock_oda_posterior)
    !! switch back to ensemble member pelist
    call set_current_pelist(CS%ensemble_pelist(CS%ensemble_id,:))

    call mpp_update_domains(CS%tv%T, CS%domains(CS%ensemble_id)%mpp_domain)
    call mpp_update_domains(CS%tv%S, CS%domains(CS%ensemble_id)%mpp_domain)

    CS%tv%T = CS%tv%T / (CS%assim_frequency * seconds_per_hour)
    CS%tv%S = CS%tv%S / (CS%assim_frequency * seconds_per_hour)

end subroutine get_posterior_tracer

  subroutine get_bias_correction_tracer(Time, CS)
    type(time_type), intent(in) :: Time !< the current model time
    type(ODA_CS), pointer :: CS !< ocean DA control structure

    integer :: i,j,k
    real, allocatable, dimension(:,:,:) :: T_bias, S_bias
    real, allocatable, dimension(:,:,:) :: mask_z
    real, allocatable, dimension(:), target :: z_in, z_edges_in
    real :: missing_value
    integer,dimension(3) :: fld_sz

    if(is_root_pe()) print *, 'Getting bias correction'

    call mpp_clock_begin(id_clock_bias_correction)
    call horiz_interp_and_extrap_tracer(CS%INC_CS%T_id,Time,1.0,CS%G,T_bias,&
            mask_z,z_in,z_edges_in,missing_value,.true.,.false.,.false.,.true.)
    call horiz_interp_and_extrap_tracer(CS%INC_CS%S_id,Time,1.0,CS%G,S_bias,&
            mask_z,z_in,z_edges_in,missing_value,.true.,.false.,.false.,.true.)

    fld_sz=shape(T_bias)
    do i=1,fld_sz(1)
       do j=1,fld_sz(2)
          do k=1,fld_sz(3)
             if (T_bias(i,j,k) .gt. 1.0E-3) T_bias(i,j,k) = 0.0
             if (S_bias(i,j,k) .gt. 1.0E-3) S_bias(i,j,k) = 0.0
          enddo
       enddo
    enddo

    CS%tv_bc%T = T_bias * CS%correction_multiplier
    CS%tv_bc%S = S_bias * CS%correction_multiplier

    call mpp_update_domains(CS%tv_bc%T, CS%domains(CS%ensemble_id)%mpp_domain)
    call mpp_update_domains(CS%tv_bc%S, CS%domains(CS%ensemble_id)%mpp_domain)

    call mpp_clock_end(id_clock_bias_correction)

  end subroutine get_bias_correction_tracer

!> Gather observations and sall ODA routines
subroutine oda(Time, CS)
    type(time_type), intent(in) :: Time   !< time type
    type(ODA_CS), pointer :: CS           !< ocean DA control structure

    integer :: i, j
    integer :: m
    integer :: yr, mon, day, hr, min, sec

    if ( Time >= CS%Time .and. associated(CS) ) then

      if ( CS%assim_method > NO_FILTER) then

        !! switch to global pelist
        call set_current_pelist(CS%filter_pelist)
        call mpp_clock_begin(id_clock_oda_filter)

        !! get profiles for current assimilation step
        call get_profiles(Time, CS%Profiles, CS%CProfiles)
#ifdef ENABLE_FILTER
        call ensemble_filter(CS%Ocean_prior, CS%Ocean_posterior, CS%CProfiles, CS%kdroot, CS%mpp_domain, CS%oda_grid)
#endif
        call mpp_clock_end(id_clock_oda_filter)

        if (CS%write_obs) call save_obs_diff(CS%CProfiles)
        !! switch back to ensemble member pelist
        call set_current_pelist(CS%ensemble_pelist(CS%ensemble_id,:))

        call get_posterior_tracer(Time, CS, increment=.true.)

      endif

      if (CS%do_bias_correction) call get_bias_correction_tracer(Time, CS)

    end if

    return
end subroutine oda

!> Finalize DA module
subroutine oda_end(CS)
  type(ODA_CS), intent(inout) :: CS !< the ocean DA control structure
  call close_profile_file(temp_fid)
  call close_profile_file(salt_fid)
end subroutine oda_end

!> Initialize DA module
subroutine init_ocean_ensemble(CS,Grid,GV,ens_size)
  type(ocean_control_struct), pointer :: CS !< Pointer to ODA control structure
  type(ocean_grid_type), pointer :: Grid !< Pointer to ocean analysis grid
  type(verticalGrid_type), pointer :: GV !< Pointer to DA vertical grid
  integer, intent(in) :: ens_size !< ensemble size

  integer :: n,is,ie,js,je,nk

  nk=GV%ke
  is=Grid%isd;ie=Grid%ied
  js=Grid%jsd;je=Grid%jed
  CS%ensemble_size=ens_size
  allocate(CS%T(is:ie,js:je,nk,ens_size));CS%T(:,:,:,:)=0.0
  allocate(CS%S(is:ie,js:je,nk,ens_size));CS%S(:,:,:,:)=0.0

  return
end subroutine init_ocean_ensemble

!> Set the next analysis time
subroutine set_analysis_time(Time,CS)
  type(time_type), intent(in) :: Time !< the current model time
  type(ODA_CS), pointer, intent(inout) :: CS !< the DA control structure

  character(len=160) :: mesg  ! The text of an error message
  integer :: yr, mon, day, hr, min, sec

  if (Time >= CS%Time) then
    CS%Time=increment_time(CS%Time,CS%assim_frequency*int(seconds_per_hour))

    call get_date(Time, yr, mon, day, hr, min, sec)
    write(mesg,*) 'Model Time: ', yr, mon, day, hr, min, sec
    call MOM_mesg("set_analysis_time: "//trim(mesg))
    call get_date(CS%time, yr, mon, day, hr, min, sec)
    write(mesg,*) 'Assimilation Time: ', yr, mon, day, hr, min, sec
    call MOM_mesg("set_analysis_time: "//trim(mesg))
  endif
  if (CS%Time < Time) then
    call MOM_error(FATAL, " set_analysis_time: " // &
         "assimilation interval appears to be shorter than " // &
         "the model timestep")
  endif
  return

end subroutine set_analysis_time

!> Write observation differences to a file
subroutine save_obs_diff(profiles)
  type(ocean_profile_type), pointer :: profiles      !< pointer to profile type
  type(ocean_profile_type), pointer :: Prof=>NULL()
  Prof=>profiles

  do while (associated(Prof))
    if(Prof%compute .and. Prof%inst_type .eq. ODA_PFL) then
       if(Prof%variable .eq. TEMP_ID) then
          call write_profile(temp_fid,Prof)
       elseif(Prof%variable .eq. SALT_ID) then
          call write_profile(salt_fid,Prof)
       endif
    endif
    Prof=>Prof%cnext
  enddo

  return
end subroutine save_obs_diff


!> Apply increments to tracers
subroutine apply_oda_tracer_increments(dt,Time_end,G,tv,h,CS)
  real,                     intent(in)    :: dt !< The tracer timestep [s]
  type(time_type), intent(in) :: Time_end       !< Time at the end of the interval
  type(ocean_grid_type),    intent(in)    :: G  !< ocean grid structure
  type(thermo_var_ptrs),    intent(inout) :: tv !< A structure pointing to various thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  &
                            intent(in)    :: h  !< layer thickness [H ~> m or kg m-2]
  type(ODA_CS), pointer,     intent(inout) :: CS !< the data assimilation structure

    !! local variables
    integer :: yr, mon, day, hr, min, sec
    integer :: i, j
    integer :: isc, iec, jsc, jec, k
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: T_inc
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: S_inc
    real, dimension(SZI_(G),SZJ_(G),SZK_(CS%Grid)) :: T
    real, dimension(SZI_(G),SZJ_(G),SZK_(CS%Grid)) :: S
    real :: missing_value

    if (.not. associated(CS)) return
    if (CS%assim_method .eq. NO_FILTER .and. (.not. CS%do_bias_correction)) return

    call mpp_clock_begin(id_clock_apply_increments)

    T_inc = 0.0; S_inc = 0.0; T = 0.0; S = 0.0
    if (CS%assim_method > 0 ) then
      T = T + CS%tv%T
      S = S + CS%tv%S
    endif
    if (CS%do_bias_correction ) then
      T = T + CS%tv_bc%T
      S = S + CS%tv_bc%S
    endif

    isc=G%isc; iec=G%iec; jsc=G%jsc; jec=G%jec
    do j=jsc,jec; do i=isc,iec
      call remapping_core_h(CS%remapCS, CS%nk, CS%h(i,j,:), T(i,j,:), &
              G%ke, h(i,j,:), T_inc(i,j,:))
      call remapping_core_h(CS%remapCS, CS%nk, CS%h(i,j,:), S(i,j,:), &
              G%ke, h(i,j,:), S_inc(i,j,:))
    enddo; enddo

    !missing_value = get_external_field_missing(CS%INC_CS%T_id)
    !do k = 1,G%ke
       !call myStats(T_inc(:,:,k),missing_value,1,G%ied,1,G%jed,k,'Applied T increments')
    !enddo
    !do k = 1,G%ke
       !call myStats(S_inc(:,:,k),missing_value,1,G%ied,1,G%jed,k,'Applied S increments')
    !enddo

    call mpp_update_domains(T_inc, G%Domain%mpp_domain)
    call mpp_update_domains(S_inc, G%Domain%mpp_domain)

    tv%T(isc:iec,jsc:jec,:)=tv%T(isc:iec,jsc:jec,:)+T_inc(isc:iec,jsc:jec,:)*dt
    tv%S(isc:iec,jsc:jec,:)=tv%S(isc:iec,jsc:jec,:)+S_inc(isc:iec,jsc:jec,:)*dt

    call mpp_update_domains(tv%T, G%Domain%mpp_domain)
    call mpp_update_domains(tv%S, G%Domain%mpp_domain)

    call enable_averaging(dt, Time_end, CS%diag_cs)
    if (CS%id_inc_t > 0) call post_data(CS%id_inc_t, T_inc, CS%diag_cs)
    if (CS%id_inc_s > 0) call post_data(CS%id_inc_s, S_inc, CS%diag_cs)
    call disable_averaging(CS%diag_cs)

    call diag_update_remap_grids(CS%diag_cs)
    call mpp_clock_end(id_clock_apply_increments)

    return

end subroutine apply_oda_tracer_increments

  subroutine set_up_global_tgrid(T_grid, CS, G)
    type(grid_type), pointer :: T_grid      !< global tracer grid
    type(ODA_CS), pointer, intent(in) :: CS !< ocean DA control structure
    type(ocean_grid_type), pointer :: G     !< domain and grid information for ocean model

    ! local variables
    real, dimension(:,:), allocatable :: global2D, global2D_old
    integer :: i, j, k

    if(.not. associated(T_grid)) allocate(T_grid)
    T_grid%ni = CS%ni
    T_grid%nj = CS%nj
    T_grid%nk = CS%nk

    allocate(T_grid%x(CS%ni,CS%nj))
    call mpp_global_field(CS%mpp_domain, CS%Grid%geolonT, T_grid%x)
    allocate(T_grid%y(CS%ni,CS%nj))
    call mpp_global_field(CS%mpp_domain, CS%Grid%geolatT, T_grid%y)

    allocate(T_grid%basin_mask(CS%ni,CS%nj))
    call mpp_global_field(CS%mpp_domain, CS%oda_grid%basin_mask, T_grid%basin_mask)

    allocate(T_grid%bathyT(CS%ni,CS%nj))
    call mpp_global_field(CS%domains(CS%ensemble_id)%mpp_domain, G%bathyT, T_grid%bathyT)

    allocate(T_grid%mask(CS%ni,CS%nj,CS%nk));   T_grid%mask(:,:,:) = 0.0
    allocate(T_grid%z(CS%ni,CS%nj,CS%nk));      T_grid%z(:,:,:) = 0.0
    allocate(global2D(CS%ni,CS%nj))
    allocate(global2D_old(CS%ni,CS%nj))

    do k = 1, CS%nk
      call mpp_global_field(CS%domains(CS%ensemble_id)%mpp_domain, CS%h(:,:,k), global2D)
      do i=1, CS%ni; do j=1, CS%nj
        if ( global2D(i,j) > 1 ) then
          T_grid%mask(i,j,k) = 1.0
        end if
      end do; end do
      if (k .eq. 1) then
        T_grid%z(:,:,k) = global2D/2
      else
        T_grid%z(:,:,k) = T_grid%z(:,:,k-1) + (global2D + global2D_old)/2
      end if
      global2D_old = global2D
    end do

    deallocate(global2D)
    deallocate(global2D_old)
  end subroutine set_up_global_tgrid

!> \namespace MOM_oda_driver_mod
!!
!! \section section_ODA The Ocean data assimilation (DA) and Ensemble Framework
!!
!! The DA framework implements ensemble capability in MOM6.   Currently, this framework
!! is enabled using the cpp directive ENSEMBLE_OCEAN.  The ensembles need to be generated
!! at the level of the calling routine for oda_init or above. The ensemble instances may
!! exist on overlapping or non-overlapping processors. The ensemble information is accessed
!! via the FMS ensemble manager. An independent PE layout is used to gather (prior) ensemble
!! member information where this information is stored in the ODA control structure.  This
!! module was developed in collaboration with Feiyu Lu and Tony Rosati in the GFDL prediction
!! group for use in their coupled ensemble framework. These interfaces should be suitable for
!! interfacing MOM6 to other data assimilation packages as well.


end module MOM_oda_driver_mod
