!> Interfaces for MOM6 ensembles and data assimilation.
module MOM_oda_driver_mod

! This file is part of MOM6. see LICENSE.md for the license.

! MOM infrastructure
use MOM_coms, only : PE_here, num_PEs
use MOM_coms, only : set_PElist, set_rootPE, Get_PElist, broadcast
use MOM_domains, only : domain2d, global_field, get_domain_extent
use MOM_domains, only : pass_var, redistribute_array, broadcast_domain
use MOM_diag_mediator, only : register_diag_field, diag_axis_init, post_data
use MOM_diag_mediator, only : enable_averaging, disable_averaging
use MOM_diag_mediator, only : diag_update_remap_grids
use MOM_ensemble_manager, only : get_ensemble_id, get_ensemble_size
use MOM_ensemble_manager, only : get_ensemble_pelist, get_ensemble_filter_pelist
use MOM_error_handler, only : stdout, stdlog, MOM_error
use MOM_io, only : SINGLE_FILE
use MOM_interp_infra, only : init_extern_field, get_external_field_info
use MOM_interp_infra, only : time_interp_extern
use MOM_interpolate, only : external_field
use MOM_remapping,    only : remappingSchemesDoc
use MOM_time_manager, only : time_type, real_to_time, get_date
use MOM_time_manager, only : operator(+), operator(>=), operator(/=)
use MOM_time_manager, only : operator(==), operator(<)
use MOM_cpu_clock, only : cpu_clock_begin, cpu_clock_end, cpu_clock_id
use MOM_horizontal_regridding, only : horiz_interp_and_extrap_tracer
! ODA Modules
use ocean_da_types_mod, only : grid_type, ocean_profile_type, ocean_control_struct
use ocean_da_core_mod, only : ocean_da_core_init, get_profiles
!This preprocessing directive enables the SPEAR online ensemble data assimilation
!configuration. Existing community based APIs for data assimilation are currently
!called offline for forecast applications using information read from a MOM6 state file.
!The SPEAR configuration (https://doi.org/10.1029/2020MS002149) calculated increments
!efficiently online. A community-based set of APIs should be implemented in place
!of the CPP directive when this is available.
#ifdef ENABLE_ECDA
use eakf_oda_mod, only : ensemble_filter
#endif
use kdtree, only : kd_root !# A kd-tree object using JEDI APIs
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
use MOM_fixed_initialization, only : MOM_initialize_topography
use MOM_coord_initialization, only : MOM_initialize_coord
use MOM_file_parser, only : read_param, get_param, param_file_type
use MOM_string_functions, only : lowercase
use MOM_ALE, only : ALE_CS, ALE_initThicknessToCoord, ALE_init, ALE_updateVerticalGridType
use MOM_domains, only : MOM_domains_init, MOM_domain_type, clone_MOM_domain
use MOM_remapping, only : remapping_CS, initialize_remapping, remapping_core_h
use MOM_regridding, only : regridding_CS, initialize_regridding
use MOM_regridding, only : regridding_main, set_regrid_params, set_h_neglect
use MOM_unit_scaling, only : unit_scale_type, unit_scaling_init
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type, verticalGridInit

implicit none ; private

public :: init_oda, oda_end, set_prior_tracer, get_posterior_tracer
public :: set_analysis_time, oda, apply_oda_tracer_increments

!>@{ CPU time clock ID
integer :: id_clock_oda_init
integer :: id_clock_bias_adjustment
integer :: id_clock_apply_increments
!>@}

#include <MOM_memory.h>

!> A structure with a pointer to a domain2d, to allow for the creation of arrays of pointers.
type :: ptr_mpp_domain
  type(domain2d), pointer :: mpp_domain => NULL() !< pointer to a domain2d
end type ptr_mpp_domain

!> A structure containing integer handles for bias adjustment of tracers
type :: INC_CS
  integer :: fldno = 0 !< The number of tracers
  type(external_field) :: T  !< The handle for the temperature file
  type(external_field) :: S  !< The handle for the salinity file
end type INC_CS

!> Control structure that contains a transpose of the ocean state across ensemble members.
type, public :: ODA_CS ; private
  type(ocean_control_struct), pointer :: Ocean_prior=> NULL() !< ensemble ocean prior states in DA space
  type(ocean_control_struct), pointer :: Ocean_posterior=> NULL() !< ensemble ocean posterior states
                                                                  !! or increments to prior in DA space
  type(ocean_control_struct), pointer :: Ocean_increment=> NULL() !< A separate structure for
                                                                  !! increment diagnostics
  integer :: nk !< number of vertical layers used for DA
  type(ocean_grid_type), pointer :: Grid => NULL() !< MOM6 grid type and decomposition for the DA
  type(ocean_grid_type), pointer :: G => NULL() !< MOM6 grid type and decomposition for the model
  type(MOM_domain_type), pointer, dimension(:) :: domains => NULL() !< Pointer to mpp_domain objects
                                                                       !! for ensemble members
  type(verticalGrid_type), pointer :: GV => NULL() !< vertical grid for DA
  type(unit_scale_type), pointer :: &
    US => NULL()    !< structure containing various unit conversion factors for DA

  type(domain2d), pointer :: mpp_domain => NULL() !< Pointer to a mpp domain object for DA
  type(grid_type), pointer :: oda_grid !< local tracer grid
  real, pointer, dimension(:,:,:) :: h => NULL() !<layer thicknesses [H ~> m or kg m-2] for DA
  real, pointer, dimension(:,:,:) :: T_tend => NULL() !<layer temperature tendency from DA [C T-1 ~> degC s-1]
  real, pointer, dimension(:,:,:) :: S_tend => NULL() !<layer salinity tendency from DA [S T-1 ~> ppt s-1]
  real, pointer, dimension(:,:,:) :: T_bc_tend => NULL() !< The layer temperature tendency due
                                                         !! to bias adjustment [C T-1 ~> degC s-1]
  real, pointer, dimension(:,:,:) :: S_bc_tend => NULL() !< The layer salinity tendency due
                                                         !! to bias adjustment [S T-1 ~> ppt s-1]
  integer :: ni          !< global i-direction grid size
  integer :: nj          !< global j-direction grid size
  logical :: reentrant_x !< grid is reentrant in the x direction
  logical :: reentrant_y !< grid is reentrant in the y direction
  logical :: tripolar_N !< grid is folded at its north edge
  logical :: symmetric !< Values at C-grid locations are symmetric
  logical :: use_basin_mask !< If true, use a basin file to delineate weakly coupled ocean basins
  logical :: do_bias_adjustment !< If true, use spatio-temporally varying climatological tendency
                                !! adjustment for Temperature and Salinity
  real :: bias_adjustment_multiplier !< A scaling for the bias adjustment [nondim]
  integer :: assim_method !< Method: NO_ASSIM,EAKF_ASSIM or OI_ASSIM
  integer :: ensemble_size !< Size of the ensemble
  integer :: ensemble_id = 0 !< id of the current ensemble member
  integer, pointer, dimension(:,:) :: ensemble_pelist !< PE list for ensemble members
  integer, pointer, dimension(:) :: filter_pelist !< PE list for ensemble members
  real :: assim_interval !< analysis interval [ T ~> s]
  ! Profiles local to the analysis domain
  type(ocean_profile_type), pointer :: Profiles => NULL() !< pointer to linked list of all available profiles
  type(ocean_profile_type), pointer :: CProfiles => NULL()!< pointer to linked list of current profiles
  type(kd_root), pointer :: kdroot => NULL() !< A structure for storing nearest neighbors
  type(ALE_CS), pointer :: ALE_CS=>NULL() !< ALE control structure for DA
  logical :: use_ALE_algorithm !< true is using ALE remapping
  type(regridding_CS) :: regridCS !< ALE control structure for regridding
  type(remapping_CS) :: remapCS !< ALE control structure for remapping
  type(time_type) :: Time !< Current Analysis time
  type(diag_ctrl), pointer :: diag_cs=> NULL() !<Pointer to diagnostics control structure
  type(INC_CS) :: INC_CS !< A Structure containing integer file handles for bias adjustment
  integer :: id_inc_t !< A diagnostic handle for the temperature climatological adjustment
  integer :: id_inc_s !< A diagnostic handle for the salinity climatological adjustment
  integer :: answer_date    !< The vintage of the order of arithmetic and expressions in the
                            !! remapping invoked by the ODA driver.  Values below 20190101 recover
                            !! the answers from the end of 2018, while higher values use updated
                            !! and more robust forms of the same expressions.
  logical :: reproduce_2018_nmme !< true if reproducing older NMME answers.
end type ODA_CS


!>@{  DA parameters
integer, parameter :: NO_ASSIM = 0, OI_ASSIM=1, EAKF_ASSIM=2
!>@}
character(len=40)  :: mdl = "MOM_oda_driver" !< This module's name.

contains

!> initialize First_guess (prior) and Analysis grid
!! information for all ensemble members
subroutine init_oda(Time, G, GV, US, diag_CS, CS)

  type(time_type), intent(in) :: Time !< The current model time.
  type(ocean_grid_type), pointer :: G !< domain and grid information for ocean model
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(diag_ctrl), target, intent(inout) :: diag_CS !< A pointer to a diagnostic control structure
  type(ODA_CS), pointer, intent(inout) :: CS  !< The DA control structure

! Local variables
  type(thermo_var_ptrs) :: tv_dummy
  type(dyn_horgrid_type), pointer :: dG=> NULL()
  type(hor_index_type), pointer :: HI=> NULL()
  type(directories) :: dirs

  type(grid_type), pointer :: T_grid !< global tracer grid
  type(param_file_type) :: PF
  integer :: n
  integer :: isd, ied, jsd, jed
  integer :: is_oda, ie_oda, js_oda, je_oda
  integer :: isd_oda, ied_oda, jsd_oda, jed_oda
  integer, dimension(4) :: fld_sz
  character(len=32) :: assim_method
  integer :: npes_pm, ens_info(6)
  character(len=30) :: coord_mode
  character(len=200) :: inputdir, basin_file
  character(len=80) :: basin_var
  character(len=80) :: remap_scheme
  character(len=80) :: bias_correction_file, inc_file
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  logical :: om4_remap_via_sub_cells ! If true, use the OM4 remapping algorithm
  real :: h_neglect, h_neglect_edge                 ! small thicknesses [H ~> m or kg m-2]

  if (associated(CS)) call MOM_error(FATAL, 'Calling oda_init with associated control structure')
  allocate(CS)

  id_clock_oda_init=cpu_clock_id('(ODA initialization)')
  call cpu_clock_begin(id_clock_oda_init)

! Use ens1 parameters , this could be changed at a later time
! if it were desirable to have alternate parameters, e.g. for the grid
! for the analysis
  call get_MOM_input(PF,dirs,ensemble_num=0)
  call unit_scaling_init(PF, CS%US)

  call get_param(PF, mdl, "ASSIM_METHOD", assim_method,  &
       "String which determines the data assimilation method "//&
       "Valid methods are: \'EAKF\',\'OI\', and \'NO_ASSIM\'", default='NO_ASSIM')
  call get_param(PF, mdl, "ASSIM_INTERVAL", CS%assim_interval,  &
       "data assimilation update interval in hours",default=-1.0,units="hours",scale=3600.*US%s_to_T)
  if (CS%assim_interval < 0.) then
     call get_param(PF, mdl, "ASSIM_FREQUENCY", CS%assim_interval,  &
          "data assimilation update  in hours. This parameter name will \n"//&
          "be deprecated in the future. ASSIM_INTERVAL should be used instead.",default=-1.0, &
          units="hours",scale=3600.*US%s_to_T)
  endif

  call get_param(PF, mdl, "USE_REGRIDDING", CS%use_ALE_algorithm , &
                "If True, use the ALE algorithm (regridding/remapping).\n"//&
                "If False, use the layered isopycnal algorithm.", default=.false. )
  call get_param(PF, mdl, "REENTRANT_X", CS%reentrant_x, &
       "If true, the domain is zonally reentrant.", default=.true.)
  call get_param(PF, mdl, "REENTRANT_Y", CS%reentrant_y, &
       "If true, the domain is meridionally reentrant.", &
       default=.false.)
  call get_param(PF, mdl, "TRIPOLAR_N", CS%tripolar_N, &
       "Use tripolar connectivity at the northern edge of the "//&
       "domain.  With TRIPOLAR_N, NIGLOBAL must be even.", &
       default=.false.)
  call get_param(PF, mdl, "APPLY_TRACER_TENDENCY_ADJUSTMENT", CS%do_bias_adjustment, &
       "If true, add a spatio-temporally varying climatological adjustment "//&
       "to temperature and salinity.", &
       default=.false.)
  if (CS%do_bias_adjustment) then
    call get_param(PF, mdl, "TRACER_ADJUSTMENT_FACTOR", CS%bias_adjustment_multiplier, &
       "A multiplicative scaling factor for the climatological tracer tendency adjustment ", &
       units="nondim", default=1.0)
  endif
  call get_param(PF, mdl, "USE_BASIN_MASK", CS%use_basin_mask, &
       "If true, add a basin mask to delineate weakly connected "//&
       "ocean basins for the purpose of data assimilation.", &
       default=.false.)

  call get_param(PF, mdl, "NIGLOBAL", CS%ni, &
       "The total number of thickness grid points in the "//&
       "x-direction in the physical domain.")
  call get_param(PF, mdl, "NJGLOBAL", CS%nj, &
       "The total number of thickness grid points in the "//&
       "y-direction in the physical domain.")
  call get_param(PF, mdl, "INPUTDIR", inputdir)
  call get_param(PF, mdl, "ODA_REMAPPING_SCHEME", remap_scheme, &
                 "This sets the reconstruction scheme used "//&
                 "for vertical remapping for all ODA variables. "//&
                 "It can be one of the following schemes: "//&
                 trim(remappingSchemesDoc), default="PPM_H4")
  call get_param(PF, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
  call get_param(PF, mdl, "ODA_ANSWER_DATE", CS%answer_date, &
               "The vintage of the order of arithmetic and expressions used by the ODA driver "//&
               "Values below 20190101 recover the answers from the end of 2018, while higher "//&
               "values use updated and more robust forms of the same expressions.", &
               default=default_answer_date, do_not_log=.not.GV%Boussinesq)
  if (.not.GV%Boussinesq) CS%answer_date = max(CS%answer_date, 20230701)

  call get_param(PF, mdl, "REPRODUCE_2018_NMME_ANSWERS", CS%reproduce_2018_nmme, &
               "Logical flag needed to reproduce older NMME forecast answers."//&
               "True gives old answers, the default of false gives different answers.", &
               default=.false.)

  inputdir = slasher(inputdir)

  select case(lowercase(trim(assim_method)))
    case('eakf')
      CS%assim_method = EAKF_ASSIM
    case('oi')
      CS%assim_method = OI_ASSIM
    case('no_assim')
      CS%assim_method = NO_ASSIM
    case default
      call MOM_error(FATAL, "Invalid assimilation method provided")
  end select

  ens_info = get_ensemble_size()
  CS%ensemble_size = ens_info(1)
  npes_pm=ens_info(3)
  CS%ensemble_id = get_ensemble_id()
  !! Switch to global pelist
  allocate(CS%ensemble_pelist(CS%ensemble_size,npes_pm))
  allocate(CS%filter_pelist(CS%ensemble_size*npes_pm))
  call get_ensemble_pelist(CS%ensemble_pelist, 'ocean')
  call get_ensemble_filter_pelist(CS%filter_pelist, 'ocean')

  call set_PElist(CS%filter_pelist)

  allocate(CS%domains(CS%ensemble_size))
  CS%domains(CS%ensemble_id)%mpp_domain => G%Domain%mpp_domain ! this should go away
  do n=1,CS%ensemble_size
    if (.not. associated(CS%domains(n)%mpp_domain)) allocate(CS%domains(n)%mpp_domain)
    call set_rootPE(CS%ensemble_pelist(n,1)) ! this line is not in Feiyu's version (needed?)
    call broadcast_domain(CS%domains(n)%mpp_domain)
  enddo
  call set_rootPE(CS%filter_pelist(1)) ! this line is not in Feiyu's version (needed?)
  CS%G => G
  allocate(CS%Grid)
  ! params NIHALO_ODA, NJHALO_ODA set the DA halo size
  call MOM_domains_init(CS%Grid%Domain, PF, param_suffix='_ODA', US=CS%US)
  allocate(HI)
  call hor_index_init(CS%Grid%Domain, HI, PF)
  call verticalGridInit( PF, CS%GV, CS%US )
  allocate(dG)
  call create_dyn_horgrid(dG, HI)
  call clone_MOM_domain(CS%Grid%Domain, dG%Domain,symmetric=.false.)
  call set_grid_metrics(dG, PF, CS%US)
  call MOM_initialize_topography(dG%bathyT, dG%max_depth, dG, PF, CS%US)
  call MOM_initialize_coord(CS%GV, CS%US, PF, tv_dummy, dG%max_depth)
  call ALE_init(PF, CS%GV, CS%US, dG%max_depth, CS%ALE_CS)
  call MOM_grid_init(CS%Grid, PF, global_indexing=.false.)
  call ALE_updateVerticalGridType(CS%ALE_CS, CS%GV)
  call copy_dyngrid_to_MOM_grid(dG, CS%Grid, CS%US)
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
  call get_param(PF, mdl, "REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                 do_not_log=.true., default=.true.)
  call get_param(PF, mdl, "ODA_REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
       "If true, use the OM4 remapping-via-subcells algorithm for ODA. "//&
       "See REMAPPING_USE_OM4_SUBCELLS for more details. "//&
       "We recommend setting this option to false.", default=om4_remap_via_sub_cells)
  call initialize_regridding(CS%regridCS, CS%GV, CS%US, dG%max_depth,PF,'oda_driver',coord_mode,'','')

  h_neglect = set_h_neglect(GV, CS%answer_date, h_neglect_edge)
  call initialize_remapping(CS%remapCS, remap_scheme, om4_remap_via_sub_cells=om4_remap_via_sub_cells, &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge, answer_date=CS%answer_date)
  call set_regrid_params(CS%regridCS, min_thickness=0.)
  isd = G%isd; ied = G%ied; jsd = G%jsd; jed = G%jed

  ! breaking with the MOM6 convention and using global indices
  !call get_domain_extent(G%Domain,is,ie,js,je,isd,ied,jsd,jed,&
  !                       isg,ieg,jsg,jeg,idg_offset,jdg_offset,symmetric)
  !isd=isd+idg_offset; ied=ied+idg_offset ! using global indexing within the DA module
  !jsd=jsd+jdg_offset; jed=jed+jdg_offset ! TODO:  switch to local indexing? (mjh)

  if (.not. associated(CS%h)) then
    allocate(CS%h(isd:ied,jsd:jed,CS%GV%ke), source=CS%GV%Angstrom_H)
    ! assign thicknesses
    call ALE_initThicknessToCoord(CS%ALE_CS, G, CS%GV, CS%h)
  endif

  allocate(CS%T_tend(isd:ied,jsd:jed,CS%GV%ke), source=0.0)
  allocate(CS%S_tend(isd:ied,jsd:jed,CS%GV%ke), source=0.0)
!  call set_axes_info(CS%Grid, CS%GV, CS%US, PF, CS%diag_cs, set_vertical=.true.) ! missing in Feiyu's fork
  allocate(CS%oda_grid)
  CS%oda_grid%x => CS%Grid%geolonT
  CS%oda_grid%y => CS%Grid%geolatT


  if (CS%use_basin_mask) then
    call get_param(PF, 'oda_driver', "BASIN_FILE", basin_file, &
          "A file in which to find the basin masks.", default="basin.nc")
    basin_file = trim(inputdir) // trim(basin_file)
    call get_param(PF, 'oda_driver', "BASIN_VAR", basin_var, &
          "The basin mask variable in BASIN_FILE.", default="basin")
    ! Need different data domain indices for the ODA ensemble basin mask.
    call get_domain_extent(CS%Grid%Domain, is_oda, ie_oda, js_oda, je_oda, isd_oda, ied_oda, jsd_oda, jed_oda)
    allocate(CS%oda_grid%basin_mask(isd_oda:ied_oda,jsd_oda:jed_oda), source=0.0)
    call MOM_read_data(basin_file, basin_var, CS%oda_grid%basin_mask, CS%Grid%domain, timelevel=1)
  endif

  ! set up diag variables for analysis increments
  CS%diag_CS => diag_CS
  CS%id_inc_t = register_diag_field('ocean_model', 'temp_increment', diag_CS%axesTL, &
       Time, 'ocean potential temperature increments', 'degC', conversion=US%C_to_degC)
  CS%id_inc_s = register_diag_field('ocean_model', 'salt_increment', diag_CS%axesTL, &
       Time, 'ocean salinity increments', 'psu', conversion=US%S_to_ppt)

  !!  get global grid information from ocean model needed for ODA initialization
  T_grid => NULL()
  call set_up_global_tgrid(T_grid, CS, G)

  call ocean_da_core_init(CS%mpp_domain, T_grid, CS%Profiles, Time)
  deallocate(T_grid)
  CS%Time = Time
  !! switch back to ensemble member pelist
  call set_PElist(CS%ensemble_pelist(CS%ensemble_id,:))

  if (CS%do_bias_adjustment) then
    call get_param(PF, mdl, "TEMP_SALT_ADJUSTMENT_FILE", bias_correction_file,  &
                "The name of the file containing temperature and salinity "//&
                "tendency adjustments", default='temp_salt_adjustment.nc')

    inc_file = trim(inputdir) // trim(bias_correction_file)
    CS%INC_CS%T = init_extern_field(inc_file, "temp_increment", &
          correct_leap_year_inconsistency=.true.,verbose=.true.,domain=G%Domain%mpp_domain)
    CS%INC_CS%S = init_extern_field(inc_file, "salt_increment", &
          correct_leap_year_inconsistency=.true.,verbose=.true.,domain=G%Domain%mpp_domain)
    call get_external_field_info(CS%INC_CS%T, size=fld_sz)
    CS%INC_CS%fldno = 2
    if (CS%nk /= fld_sz(3)) call MOM_error(FATAL,'Increment levels /= ODA levels')

    allocate(CS%T_bc_tend(G%isd:G%ied,G%jsd:G%jed,CS%GV%ke), source=0.0)
    allocate(CS%S_bc_tend(G%isd:G%ied,G%jsd:G%jed,CS%GV%ke), source=0.0)
  endif

  call cpu_clock_end(id_clock_oda_init)

!  if (CS%write_obs) then
!    temp_fid = open_profile_file("temp_"//trim(obs_file))
!    salt_fid = open_profile_file("salt_"//trim(obs_file))
!  end if

end subroutine init_oda

!> Copy ensemble member tracers to ensemble vector.
subroutine set_prior_tracer(Time, G, GV, h, tv, CS)
  type(time_type), intent(in)    :: Time !< The current model time
  type(ocean_grid_type), pointer :: G !< domain and grid information for ocean model
  type(verticalGrid_type),               intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                 intent(in) :: tv   !< A structure pointing to various thermodynamic variables

  type(ODA_CS), pointer :: CS !< ocean DA control structure
  real, dimension(SZI_(G),SZJ_(G),CS%nk) :: T  ! Temperature on the analysis grid [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),CS%nk) :: S  ! Salinity on the analysis grid [S ~> ppt]
  integer :: i, j, m
  integer :: isc, iec, jsc, jec

  ! return if not time for analysis
  if (Time < CS%Time) return

  if (.not. associated(CS%Grid)) call MOM_ERROR(FATAL,'ODA_CS ensemble horizontal grid not associated')
  if (.not. associated(CS%GV)) call MOM_ERROR(FATAL,'ODA_CS ensemble vertical grid not associated')

  !! switch to global pelist
  call set_PElist(CS%filter_pelist)
  !call MOM_mesg('Setting prior')

  ! computational domain for the analysis grid
  isc=CS%Grid%isc;iec=CS%Grid%iec;jsc=CS%Grid%jsc;jec=CS%Grid%jec
  ! array extents for the ensemble member
  !call get_domain_extent(CS%domains(CS%ensemble_id),is,ie,js,je,isd,ied,jsd,jed,&
  !     isg,ieg,jsg,jeg,idg_offset,jdg_offset,symmetric)
  ! remap temperature and salinity from the ensemble member to the analysis grid
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    call remapping_core_h(CS%remapCS, GV%ke, h(i,j,:), tv%T(i,j,:), &
                          CS%nk, CS%h(i,j,:), T(i,j,:))
    call remapping_core_h(CS%remapCS, GV%ke, h(i,j,:), tv%S(i,j,:), &
                          CS%nk, CS%h(i,j,:), S(i,j,:))
  enddo ; enddo
  ! cast ensemble members to the analysis domain
  do m=1,CS%ensemble_size
    call redistribute_array(CS%domains(m)%mpp_domain, T,&
         CS%mpp_domain, CS%Ocean_prior%T(:,:,:,m), complete=.true.)
    call redistribute_array(CS%domains(m)%mpp_domain, S,&
         CS%mpp_domain, CS%Ocean_prior%S(:,:,:,m), complete=.true.)
  enddo

  do m=1,CS%ensemble_size
    call pass_var(CS%Ocean_prior%T(:,:,:,m),CS%Grid%domain)
    call pass_var(CS%Ocean_prior%S(:,:,:,m),CS%Grid%domain)
  enddo

  !! switch back to ensemble member pelist
  call set_PElist(CS%ensemble_pelist(CS%ensemble_id,:))

  return

end subroutine set_prior_tracer

!> Returns posterior adjustments or full state
!!Note that only those PEs associated with an ensemble member receive data
subroutine get_posterior_tracer(Time, CS, increment)
  type(time_type), intent(in) :: Time !< the current model time
  type(ODA_CS), pointer :: CS !< ocean DA control structure
  logical, optional, intent(in) :: increment !< True if returning increment only

  type(ocean_control_struct), pointer :: Ocean_increment=>NULL()
  integer :: m
  logical :: get_inc


  ! return if not analysis time (retain pointers for h and tv)
  if (Time < CS%Time .or. CS%assim_method == NO_ASSIM) return


  !! switch to global pelist
  call set_PElist(CS%filter_pelist)
  call MOM_mesg('Getting posterior')

  !! Calculate and redistribute increments to CS%tv right after assimilation
  !! Retain CS%tv to calculate increments for IAU updates CS%tv_inc otherwise
  get_inc = .true.
  if (present(increment)) get_inc = increment

  if (get_inc) then
    CS%Ocean_increment%T = CS%Ocean_posterior%T - CS%Ocean_prior%T
    CS%Ocean_increment%S = CS%Ocean_posterior%S - CS%Ocean_prior%S
  endif
  ! It may be necessary to check whether the increment and ocean state have the
  ! same dimensionally rescaled units.
  do m=1,CS%ensemble_size
    if (get_inc) then
      call redistribute_array(CS%mpp_domain, CS%Ocean_increment%T(:,:,:,m),&
           CS%domains(m)%mpp_domain, CS%T_tend, complete=.true.)
      call redistribute_array(CS%mpp_domain, CS%Ocean_increment%S(:,:,:,m),&
           CS%domains(m)%mpp_domain, CS%S_tend, complete=.true.)
    else
      call redistribute_array(CS%mpp_domain, CS%Ocean_posterior%T(:,:,:,m),&
           CS%domains(m)%mpp_domain, CS%T_tend, complete=.true.)
      call redistribute_array(CS%mpp_domain, CS%Ocean_posterior%S(:,:,:,m),&
           CS%domains(m)%mpp_domain, CS%S_tend, complete=.true.)
    endif
  enddo


  !! switch back to ensemble member pelist
  call set_PElist(CS%ensemble_pelist(CS%ensemble_id,:))

  call pass_var(CS%T_tend,CS%domains(CS%ensemble_id))
  call pass_var(CS%S_tend,CS%domains(CS%ensemble_id))

  !convert to a tendency (degC or PSU per second)
  CS%T_tend = CS%T_tend / (CS%assim_interval)
  CS%S_tend = CS%S_tend / (CS%assim_interval)


end subroutine get_posterior_tracer

!> Gather observations and call ODA routines
subroutine oda(Time, CS)
  type(time_type), intent(in) :: Time !< the current model time
  type(oda_CS), pointer :: CS !< A pointer the ocean DA control structure

  if ( Time >= CS%Time ) then

    !! switch to global pelist
    call set_PElist(CS%filter_pelist)
    call get_profiles(Time, CS%Profiles, CS%CProfiles)
#ifdef ENABLE_ECDA
    call ensemble_filter(CS%Ocean_prior, CS%Ocean_posterior, CS%CProfiles, CS%kdroot, CS%mpp_domain, CS%oda_grid)
#endif
    !! switch back to ensemble member pelist
    call set_PElist(CS%ensemble_pelist(CS%ensemble_id,:))
    call get_posterior_tracer(Time, CS, increment=.true.)
    if (CS%do_bias_adjustment) call get_bias_correction_tracer(Time, CS%US, CS)

  endif

  return
end subroutine oda

subroutine get_bias_correction_tracer(Time, US, CS)
  type(time_type), intent(in) :: Time !< the current model time
  type(unit_scale_type), intent(in) :: US !< A dimensional unit scaling type
  type(ODA_CS), pointer :: CS !< ocean DA control structure

  ! Local variables
  real, allocatable, dimension(:,:,:) :: T_bias ! Estimated temperature tendency bias [C T-1 ~> degC s-1]
  real, allocatable, dimension(:,:,:) :: S_bias ! Estimated salinity tendency bias [S T-1 ~> ppt s-1]
  real, allocatable, dimension(:,:,:) :: valid_flag ! Valid value flag on the horizontal model grid
                                                    ! and input-file vertical levels [nondim]
  real, allocatable, dimension(:), target :: z_in       ! Cell center depths for input data [Z ~> m]
  real, allocatable, dimension(:), target :: z_edges_in ! Cell edge depths for input data [Z ~> m]
  real :: missing_value ! A value indicating that there is no valid input data at this point [CU ~> conc]
  integer, dimension(3) :: fld_sz
  integer :: i,j,k


  call cpu_clock_begin(id_clock_bias_adjustment)
  call horiz_interp_and_extrap_tracer(CS%INC_CS%T, Time, CS%G, T_bias, &
            valid_flag, z_in, z_edges_in, missing_value, scale=US%degC_to_C*US%s_to_T, spongeOngrid=.true., &
            answer_date=CS%answer_date)
  call horiz_interp_and_extrap_tracer(CS%INC_CS%S, Time, CS%G, S_bias, &
            valid_flag, z_in, z_edges_in, missing_value, scale=US%ppt_to_S*US%s_to_T, spongeOngrid=.true., &
            answer_date=CS%answer_date)

  ! This should be replaced to use mask_z instead of the following lines
  ! which are intended to zero land values using an arbitrary limit.
  fld_sz=shape(T_bias)
  if (CS%reproduce_2018_nmme) then
    do i=1,fld_sz(1)
      do j=1,fld_sz(2)
        do k=1,fld_sz(3)
          ! The following two lines are needed for backward compatibility for NMME answers (2018 vintage)
          ! These were implemented to catch missing values, so large values are excluded.
          if (T_bias(i,j,k) > 1.0E-3*US%degC_to_C) T_bias(i,j,k) = 0.0
          if (S_bias(i,j,k) > 1.0E-3*US%ppt_to_S) S_bias(i,j,k) = 0.0
        enddo
      enddo
    enddo
  else
    do i=1,fld_sz(1)
      do j=1,fld_sz(2)
        do k=1,fld_sz(3)
          if (valid_flag(i,j,k)==0.) then
            T_bias(i,j,k)=0.0
            S_bias(i,j,k)=0.0
          endif
        enddo
      enddo
    enddo
  endif

  CS%T_bc_tend = T_bias * CS%bias_adjustment_multiplier
  CS%S_bc_tend = S_bias * CS%bias_adjustment_multiplier

  call pass_var(CS%T_bc_tend, CS%domains(CS%ensemble_id))
  call pass_var(CS%S_bc_tend, CS%domains(CS%ensemble_id))

  call cpu_clock_end(id_clock_bias_adjustment)

end subroutine get_bias_correction_tracer

!> Finalize DA module
subroutine oda_end(CS)
  type(ODA_CS), intent(inout) :: CS !< the ocean DA control structure

end subroutine oda_end

!> Initialize DA module
subroutine init_ocean_ensemble(CS,Grid,GV,ens_size)
  type(ocean_control_struct), pointer :: CS !< Pointer to ODA control structure
  type(ocean_grid_type), pointer :: Grid !< Pointer to ocean analysis grid
  type(verticalGrid_type), pointer :: GV !< Pointer to DA vertical grid
  integer, intent(in) :: ens_size !< ensemble size

  integer :: is, ie, js, je, nk

  nk=GV%ke
  is=Grid%isd;ie=Grid%ied
  js=Grid%jsd;je=Grid%jed
  CS%ensemble_size=ens_size
  allocate(CS%T(is:ie,js:je,nk,ens_size))
  allocate(CS%S(is:ie,js:je,nk,ens_size))
  allocate(CS%SSH(is:ie,js:je,ens_size))
!  allocate(CS%id_t(ens_size), source=-1)
!  allocate(CS%id_s(ens_size), source=-1)
!  allocate(CS%U(is:ie,js:je,nk,ens_size))
!  allocate(CS%V(is:ie,js:je,nk,ens_size))
!  allocate(CS%id_u(ens_size), source=-1)
!  allocate(CS%id_v(ens_size), source=-1)
!  allocate(CS%id_ssh(ens_size), source=-1)

  return
end subroutine init_ocean_ensemble

!> Set the next analysis time
subroutine set_analysis_time(Time,CS)
  type(time_type), intent(in) :: Time !< the current model time
  type(ODA_CS), pointer, intent(inout) :: CS !< the DA control structure

  character(len=160) :: mesg  ! The text of an error message
  integer :: yr, mon, day, hr, min, sec

  if (Time >= CS%Time) then
    ! increment the analysis time to the next step
    CS%Time = CS%Time + real_to_time(CS%US%T_to_s*(CS%assim_interval))

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


!> Apply increments to tracers
subroutine apply_oda_tracer_increments(dt, Time_end, G, GV, tv, h, CS)
  real,                     intent(in)    :: dt !< The tracer timestep [T ~> s]
  type(time_type), intent(in)             :: Time_end !< Time at the end of the interval
  type(ocean_grid_type),    intent(in)    :: G  !< ocean grid structure
  type(verticalGrid_type),  intent(in)    :: GV !< The ocean's vertical grid structure
  type(thermo_var_ptrs),    intent(inout) :: tv !< A structure pointing to various thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h  !< layer thickness [H ~> m or kg m-2]
  type(ODA_CS), pointer                   :: CS !< the data assimilation structure

  !! local variables
  integer :: i, j
  integer :: isc, iec, jsc, jec
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: T_tend_inc !< an adjustment to the temperature
                                                    !! tendency [C T-1 ~> degC s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: S_tend_inc !< an adjustment to the salinity
                                                    !! tendency [S T-1 ~> ppt s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(CS%Grid)) :: T_tend !< The temperature tendency adjustment from
                                                           !! DA [C T-1 ~> degC s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(CS%Grid)) :: S_tend !< The salinity tendency adjustment from DA
                                                          !! [S T-1 ~> ppt s-1]

  if (.not. associated(CS)) return
  if (CS%assim_method == NO_ASSIM .and. (.not. CS%do_bias_adjustment)) return

  call cpu_clock_begin(id_clock_apply_increments)

  T_tend_inc(:,:,:) = 0.0; S_tend_inc(:,:,:) = 0.0; T_tend(:,:,:) = 0.0; S_tend(:,:,:) = 0.0
  if (CS%assim_method > 0 ) then
    T_tend = T_tend + CS%T_tend
    S_tend = S_tend + CS%S_tend
  endif
  if (CS%do_bias_adjustment ) then
    T_tend = T_tend + CS%T_bc_tend
    S_tend = S_tend + CS%S_bc_tend
  endif

  isc=G%isc; iec=G%iec; jsc=G%jsc; jec=G%jec
  do j=jsc,jec; do i=isc,iec
    call remapping_core_h(CS%remapCS, CS%nk, CS%h(i,j,:), T_tend(i,j,:), &
                          G%ke, h(i,j,:), T_tend_inc(i,j,:))
    call remapping_core_h(CS%remapCS, CS%nk, CS%h(i,j,:), S_tend(i,j,:), &
                          G%ke, h(i,j,:), S_tend_inc(i,j,:))
  enddo; enddo


  call pass_var(T_tend_inc, G%Domain)
  call pass_var(S_tend_inc, G%Domain)

  tv%T(isc:iec,jsc:jec,:) = tv%T(isc:iec,jsc:jec,:) + T_tend_inc(isc:iec,jsc:jec,:)*dt
  tv%S(isc:iec,jsc:jec,:) = tv%S(isc:iec,jsc:jec,:) + S_tend_inc(isc:iec,jsc:jec,:)*dt

  call pass_var(tv%T, G%Domain)
  call pass_var(tv%S, G%Domain)

  call enable_averaging(dt, Time_end, CS%diag_CS)
  if (CS%id_inc_t > 0) call post_data(CS%id_inc_t, T_tend_inc, CS%diag_CS)
  if (CS%id_inc_s > 0) call post_data(CS%id_inc_s, S_tend_inc, CS%diag_CS)
  call disable_averaging(CS%diag_CS)

  call diag_update_remap_grids(CS%diag_CS)
  call cpu_clock_end(id_clock_apply_increments)


end subroutine apply_oda_tracer_increments

!> Set up the grid of thicknesses at tracer points throughout the global domain
  subroutine set_up_global_tgrid(T_grid, CS, G)
    type(grid_type), pointer :: T_grid !< global tracer grid
    type(ODA_CS), pointer, intent(in) :: CS !< A pointer to DA control structure.
    type(ocean_grid_type), pointer :: G !< domain and grid information for ocean model

    ! local variables
    real, dimension(:,:), allocatable :: &
      global2D, &  ! A layer thickness in the entire global domain [H ~> m or kg m-2]
      global2D_old ! The thickness of the layer above the one in global2D in the entire
                   ! global domain [H ~> m or kg m-2]
    integer :: i, j, k

    !    get global grid information from ocean_model
    T_grid=>NULL()
    !if (associated(T_grid)) call MOM_error(FATAL,'MOM_oda_driver:set_up_global_tgrid called with associated T_grid')

    allocate(T_grid)
    T_grid%ni = CS%ni
    T_grid%nj = CS%nj
    T_grid%nk = CS%nk
    allocate(T_grid%x(CS%ni,CS%nj))
    allocate(T_grid%y(CS%ni,CS%nj))
    allocate(T_grid%bathyT(CS%ni,CS%nj))
    call global_field(CS%mpp_domain, CS%Grid%geolonT, T_grid%x)
    call global_field(CS%mpp_domain, CS%Grid%geolatT, T_grid%y)
    call global_field(CS%domains(CS%ensemble_id)%mpp_domain, G%bathyT, T_grid%bathyT)
    if (CS%use_basin_mask) then
      allocate(T_grid%basin_mask(CS%ni,CS%nj))
      call global_field(CS%mpp_domain, CS%oda_grid%basin_mask, T_grid%basin_mask)
    endif
    allocate(T_grid%mask(CS%ni,CS%nj,CS%nk), source=0.0)
    allocate(T_grid%z(CS%ni,CS%nj,CS%nk), source=0.0)
    allocate(global2D(CS%ni,CS%nj))
    allocate(global2D_old(CS%ni,CS%nj))

    do k = 1, CS%nk
      call global_field(G%Domain%mpp_domain, CS%h(:,:,k), global2D)
      do i=1,CS%ni ; do j=1,CS%nj
        ! ###Does the next line need to be revised?  Perhaps it should be
        ! if ( global2D(i,j) > 1.0*GV%H_to_m ) then
        if ( global2D(i,j) > 1 ) then
           T_grid%mask(i,j,k) = 1.0
        endif
      enddo; enddo
      if (k == 1) then
         T_grid%z(:,:,k) = global2D/2
      else
         T_grid%z(:,:,k) = T_grid%z(:,:,k-1) + (global2D + global2D_old)/2
      endif
      global2D_old = global2D
    enddo

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
