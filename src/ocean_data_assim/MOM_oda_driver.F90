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

  use fms_mod, only : open_namelist_file, close_file, check_nml_error
  use fms_mod, only : error_mesg, FATAL
  use mpp_mod, only : stdout, stdlog, mpp_error, npes=>mpp_npes,pe=>mpp_pe
  use mpp_mod, only : set_current_pelist => mpp_set_current_pelist
  use mpp_mod, only : set_root_pe => mpp_set_root_pe
  use mpp_mod, only : mpp_sync_self, mpp_sum, get_pelist=>mpp_get_current_pelist, mpp_root_pe
  use mpp_mod, only : set_stack_size=>mpp_set_stack_size, broadcast=>mpp_broadcast
  use mpp_mod, only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use mpp_io_mod, only : io_set_stack_size=>mpp_io_set_stack_size
  use mpp_io_mod, only : MPP_SINGLE,MPP_MULTI
  use mpp_domains_mod, only : domain2d, mpp_global_field
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain
  use mpp_domains_mod, only : mpp_redistribute, mpp_broadcast_domain
  use mpp_domains_mod, only : set_domains_stack_size=>mpp_domains_set_stack_size
  use ensemble_manager_mod, only : get_ensemble_id, get_ensemble_size
  use ensemble_manager_mod, only : get_ensemble_pelist, get_ensemble_filter_pelist
  use time_manager_mod, only : time_type, decrement_time, increment_time
  use time_manager_mod, only : get_date, get_time, operator(>=),operator(/=),operator(==),operator(<)
  use constants_mod, only : radius, epsln
  ! ODA Modules
  use oda_types_mod, only : grid_type, ocean_profile_type, ocean_control_struct
  use oda_core_mod, only : oda_core_init, get_profiles
#ifdef ENABLE_ECDA
  use eakf_oda_mod, only : ensemble_filter
#endif
  use write_ocean_data_mod, only : open_profile_file
  use write_ocean_data_mod, only : write_profile,close_profile_file
  use kdtree, only : kd_root !# JEDI
  ! MOM Modules
  use MOM_io, only : slasher, MOM_read_data
  use MOM_error_handler, only : FATAL, WARNING, MOM_error, is_root_pe
  use MOM_get_input, only : get_MOM_input, directories
  use MOM_variables, only : thermo_var_ptrs
  use MOM_grid, only : ocean_grid_type, MOM_grid_init
  use MOM_grid_initialize, only : set_grid_metrics
  use MOM_hor_index, only : hor_index_type, hor_index_init
  use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
  use MOM_transcribe_grid, only : copy_dyngrid_to_MOM_grid, copy_MOM_grid_to_dyngrid
  use MOM_fixed_initialization, only : MOM_initialize_fixed, MOM_initialize_topography
  use MOM_coord_initialization, only : MOM_initialize_coord
  use MOM_verticalGrid, only : verticalGrid_type, verticalGridInit
  use MOM_file_parser, only : read_param, get_param, param_file_type
  use MOM_string_functions, only : lowercase
  use MOM_ALE, only : ALE_CS, ALE_initThicknessToCoord, ALE_init, ALE_updateVerticalGridType
  use MOM_domains, only : MOM_domains_init, MOM_domain_type, clone_MOM_domain
  use MOM_remapping, only : remapping_CS, initialize_remapping, remapping_core_h
  use MOM_regridding, only : regridding_CS, initialize_regridding
  use MOM_regridding, only : regridding_main, set_regrid_params

  implicit none
  private

  public :: init_oda, oda_end, set_prior_tracer, get_posterior_tracer
  public :: set_analysis_time, oda, save_obs_diff, apply_oda_tracer_increments

#include <MOM_memory.h>

  type, public :: ODA_CS; private
     type(ocean_control_struct), pointer :: Ocean_prior=> NULL() !< ensemble ocean prior states in DA space
     type(ocean_control_struct), pointer :: Ocean_posterior=> NULL() !< ensemble ocean posterior states
                                                                     !! or increments to prior in DA space
     integer :: nk !< number of vertical layers used for DA
     type(ocean_grid_type), pointer :: Grid => NULL() !< MOM6 grid type and decomposition for the DA
     type(pointer_mpp_domain), pointer, dimension(:) :: domains => NULL() !< Pointer to mpp_domain objects
                                                                          !! for ensemble members
     type(verticalGrid_type), pointer :: GV => NULL() !< vertical grid for DA
     type(domain2d), pointer :: mpp_domain => NULL() !< Pointer to a mpp domain object for DA
     type(grid_type), pointer :: oda_grid !< local tracer grid
     real, pointer, dimension(:,:,:) :: h => NULL() !<layer thicknesses (m or kg/m2) for DA
     type(thermo_var_ptrs), pointer :: tv => NULL() !< pointer to thermodynamic variables
     integer :: ni, nj !< global grid size
     logical :: reentrant_x !< grid is reentrant in the x direction
     logical :: reentrant_y !< grid is reentrant in the y direction
     logical :: tripolar_N !< grid is folded at its north edge
     logical :: symmetric !< Values at C-grid locations are symmetric
     integer :: assim_method !< Method: NO_ASSIM,EAKF_ASSIM or OI_ASSIM
     integer :: ensemble_size !< Size of the ensemble
     integer :: ensemble_id = 0 !< id of the current ensemble member
     integer, pointer, dimension(:,:) :: ensemble_pelist !< PE list for ensemble members
     integer, pointer, dimension(:) :: filter_pelist !< PE list for ensemble members
     integer :: assim_frequency !< analysis interval in hours
     ! Profiles local to the analysis domain
     type(ocean_profile_type), pointer :: Profiles => NULL() !< pointer to linked list of all available profiles
     type(ocean_profile_type), pointer :: CProfiles => NULL()!< pointer to linked list of current profiles
     type(kd_root), pointer :: kdroot
     type(ALE_CS), pointer :: ALE_CS=>NULL() !< ALE control structure for DA
     logical :: use_ALE_algorithm !< true is using ALE remapping
     type(regridding_CS) :: regridCS !< ALE control structure for regridding
     type(remapping_CS) :: remapCS !< ALE control structure for remapping
     type(time_type) :: Time !< Current Analysis time
  end type ODA_CS

  type :: pointer_mpp_domain
     type(domain2d), pointer :: mpp_domain => NULL()
  end type pointer_mpp_domain

  real, dimension(:,:,:), allocatable :: T, S
  type(ocean_control_struct), pointer :: Ocean_increment=>NULL()

  integer, parameter :: NO_ASSIM = 0, OI_ASSIM=1, EAKF_ASSIM=2

contains

!V initialize First_guess (prior) and Analysis grid
!! information for all ensemble members
!!
  subroutine init_oda(Time, G, GV, CS)

    type(time_type), intent(in) :: Time !< The current model time.
    type(ocean_grid_type), pointer :: G !< domain and grid information for ocean model
    type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
    type(ODA_CS), pointer, intent(inout) :: CS

 ! Local variables
    type(thermo_var_ptrs) :: tv_dummy
    type(dyn_horgrid_type), pointer :: dG=> NULL()
    type(hor_index_type), pointer :: HI=> NULL()
    type(directories) :: dirs

    type(grid_type), pointer :: T_grid !< global tracer grid
    real, dimension(:,:), allocatable :: global2D, global2D_old
    real, dimension(:), allocatable :: lon1D, lat1D, glon1D, glat1D
    type(param_file_type) :: PF
    integer :: n, m, k, i, j, nk, id_oda_init
    integer :: is,ie,js,je,isd,ied,jsd,jed
    integer :: stdout_unit
    character(len=32) :: assim_method
    integer :: npes_pm, ens_info(6), ni, nj
    character(len=128) :: mesg
    character(len=32) :: fldnam
    character(len=30) :: coord_mode
    character(len=200) :: inputdir, basin_file
    logical :: reentrant_x, reentrant_y, tripolar_N, symmetric

    if (associated(CS)) call mpp_error(FATAL,'Calling oda_init with associated control structure')
    allocate(CS)
! Use ens1 parameters , this could be changed at a later time
! if it were desirable to have alternate parameters, e.g. for the grid
! for the analysis
    call get_MOM_input(PF,dirs,ensemble_num=0)
    call get_param(PF, "MOM", "ASSIM_METHOD", assim_method,  &
         "String which determines the data assimilation method" // &
         "Valid methods are: \'EAKF\',\'OI\', and \'NO_ASSIM\'", default='NO_ASSIM')
    call get_param(PF, "MOM", "ASSIM_FREQUENCY", CS%assim_frequency,  &
         "data assimilation frequency in hours")
    call get_param(PF, "MOM", "USE_REGRIDDING", CS%use_ALE_algorithm , &
                  "If True, use the ALE algorithm (regridding/remapping).\n"//&
                  "If False, use the layered isopycnal algorithm.", default=.false. )
    call get_param(PF, "MOM", "REENTRANT_X", CS%reentrant_x, &
         "If true, the domain is zonally reentrant.", default=.true.)
    call get_param(PF, "MOM", "REENTRANT_Y", CS%reentrant_y, &
         "If true, the domain is meridionally reentrant.", &
         default=.false.)
    call get_param(PF,"MOM", "TRIPOLAR_N", CS%tripolar_N, &
         "Use tripolar connectivity at the northern edge of the \n"//&
         "domain.  With TRIPOLAR_N, NIGLOBAL must be even.", &
         default=.false.)
    call get_param(PF,"MOM", "NIGLOBAL", CS%ni, &
         "The total number of thickness grid points in the \n"//&
         "x-direction in the physical domain.")
    call get_param(PF,"MOM", "NJGLOBAL", CS%nj, &
         "The total number of thickness grid points in the \n"//&
         "y-direction in the physical domain.")
    call get_param(PF, 'MOM', "INPUTDIR", inputdir)
    inputdir = slasher(inputdir)

    select case(lowercase(trim(assim_method)))
    case('eakf')
        CS%assim_method = EAKF_ASSIM
    case('oi')
       CS%assim_method = OI_ASSIM
    case('no_assim')
        CS%assim_method = NO_ASSIM
    case default
        call mpp_error(FATAL,'Invalid assimilation method provided')
    end select

    ens_info = get_ensemble_size()
    CS%ensemble_size = ens_info(1)
    npes_pm=ens_info(3)
    CS%ensemble_id = get_ensemble_id()
    !! Switch to global pelist
    allocate(CS%ensemble_pelist(CS%ensemble_size,npes_pm))
    allocate(CS%filter_pelist(CS%ensemble_size*npes_pm))
    call get_ensemble_pelist(CS%ensemble_pelist,'ocean')
    call get_ensemble_filter_pelist(CS%filter_pelist,'ocean')

    call set_current_pelist(CS%filter_pelist)
    if(is_root_pe()) print *, 'Initialize ODA'

    id_oda_init = mpp_clock_id('(ODA initialization computation)')
    call mpp_clock_begin(id_oda_init)

    allocate(CS%domains(CS%ensemble_size))
    CS%domains(CS%ensemble_id)%mpp_domain => G%Domain%mpp_domain
    do n=1,CS%ensemble_size
      if(.not. associated(CS%domains(n)%mpp_domain)) allocate(CS%domains(n)%mpp_domain)
      call mpp_broadcast_domain(CS%domains(n)%mpp_domain)
    enddo

    allocate(CS%Grid)
    ! params NIHALO_ODA, NJHALO_ODA set the DA halo size
    call MOM_domains_init(CS%Grid%Domain,PF,param_suffix='_ODA')
    allocate(HI)
    call hor_index_init(CS%Grid%Domain, HI, PF)
    call verticalGridInit( PF, CS%GV )
    allocate(dG)
    call create_dyn_horgrid(dG,HI)
    call clone_MOM_domain(CS%Grid%Domain, dG%Domain,symmetric=.false.)
    call set_grid_metrics(dG,PF)
    call MOM_initialize_topography(dg%bathyT,dG%max_depth,dG,PF)
    call MOM_initialize_coord(CS%GV, PF, .false., &
         dirs%output_directory, tv_dummy, dG%max_depth)
    call ALE_init(PF, CS%GV, dG%max_depth, CS%ALE_CS)
    call MOM_grid_init(CS%Grid, PF)
    call ALE_updateVerticalGridType(CS%ALE_CS,CS%GV)
    call copy_dyngrid_to_MOM_grid(dG, CS%Grid)
    CS%mpp_domain => CS%Grid%Domain%mpp_domain
    CS%Grid%ke = CS%GV%ke
    CS%nk = CS%GV%ke
    ! initialize storage for prior and posterior
    allocate(CS%Ocean_prior)
    call init_ocean_ensemble(CS%Ocean_prior,CS%Grid,CS%GV,CS%ensemble_size)
    allocate(CS%Ocean_posterior)
    call init_ocean_ensemble(CS%Ocean_posterior,CS%Grid,CS%GV,CS%ensemble_size)
    allocate(CS%tv)

    call get_param(PF, 'oda_driver', "REGRIDDING_COORDINATE_MODE", coord_mode, &
         "Coordinate mode for vertical regridding.", &
         default="ZSTAR", fail_if_missing=.false.)
    call initialize_regridding(CS%regridCS, CS%GV, dG%max_depth,PF,'oda_driver',coord_mode,'','')
    call initialize_remapping(CS%remapCS,'PCM')
    call set_regrid_params(CS%regridCS, min_thickness=0.)
    isd = G%isd; ied = G%ied; jsd = G%jsd; jed = G%jed
    if(.not. associated(CS%h)) then
        allocate(CS%h(isd:ied,jsd:jed,CS%GV%ke)); CS%h(:,:,:) = CS%GV%Angstrom
! assign thicknesses
        call ALE_initThicknessToCoord(CS%ALE_CS,G,CS%GV,CS%h)
    endif
    allocate(CS%tv%T(isd:ied,jsd:jed,CS%GV%ke)); CS%tv%T(:,:,:)=0.0
    allocate(CS%tv%S(isd:ied,jsd:jed,CS%GV%ke)); CS%tv%S(:,:,:)=0.0
    allocate(T(isd:ied,jsd:jed,CS%nk)); T = 0.0
    allocate(S(isd:ied,jsd:jed,CS%nk)); S = 0.0

    isd = CS%Grid%isd; ied = CS%Grid%ied; jsd = CS%Grid%jsd; jed = CS%Grid%jed
    allocate(CS%oda_grid)
    CS%oda_grid%x => CS%Grid%geolonT
    CS%oda_grid%y => CS%Grid%geolatT

    call get_param(PF, 'oda_driver', "BASIN_FILE", basin_file, &
            "A file in which to find the basin masks, in variable 'basin'.", &
            default="basin.nc")
    basin_file = trim(inputdir) // trim(basin_file)
    allocate(CS%oda_grid%basin_mask(isd:ied,jsd:jed))
    CS%oda_grid%basin_mask(:,:) = 0.0
    call MOM_read_data(basin_file,'basin',CS%oda_grid%basin_mask,CS%Grid%domain, timelevel=1)

!    get global grid information from ocean_model
    allocate(T_grid)
    allocate(T_grid%x(CS%ni,CS%nj))
    allocate(T_grid%y(CS%ni,CS%nj))
    allocate(T_grid%basin_mask(CS%ni,CS%nj))
    call mpp_global_field(CS%mpp_domain, CS%Grid%geolonT, T_grid%x)
    call mpp_global_field(CS%mpp_domain, CS%Grid%geolatT, T_grid%y)
    call mpp_global_field(CS%mpp_domain, CS%oda_grid%basin_mask, T_grid%basin_mask)
    T_grid%ni = CS%ni
    T_grid%nj = CS%nj
    T_grid%nk = CS%nk
    allocate(T_grid%mask(CS%ni,CS%nj,CS%nk))
    allocate(T_grid%z(CS%ni,CS%nj,CS%nk))
    allocate(global2D(CS%ni,CS%nj))
    allocate(global2D_old(CS%ni,CS%nj))
    T_grid%mask(:,:,:) = 0.0
    T_grid%z(:,:,:) = 0.0
 
    do k = 1, CS%nk
      call mpp_global_field(G%Domain%mpp_domain, CS%h(:,:,k), global2D)
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
 
    call oda_core_init(CS%mpp_domain, T_grid, CS%Profiles, Time)
 
    !CS%Time=Time
    CS%Time=increment_time(Time,CS%assim_frequency*3600)
    !! switch back to ensemble member pelist
    call mpp_clock_end(id_oda_init)
    call set_current_pelist(CS%ensemble_pelist(CS%ensemble_id,:))
  end subroutine init_oda

  subroutine set_prior_tracer(Time, G, GV, h, tv, CS)
    type(time_type), intent(in)    :: Time !< The current model time
    type(ocean_grid_type), pointer :: G !< domain and grid information for ocean model
    type(verticalGrid_type),               intent(in)    :: GV   !< The ocean's vertical grid structure
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h    !< Layer thicknesses, in H (usually m or kg m-2)
    type(thermo_var_ptrs),                 intent(in) :: tv   !< A structure pointing to various thermodynamic variables
    type(ODA_CS), pointer :: CS !< ocean DA control structure

    integer :: i, j, m, n, ss
    integer :: isc, iec, jsc, jec
    integer :: id, id_oda_prior
    logical :: used

    ! return if not time for analysis
    if (Time < CS%Time) return

    if (.not. ASSOCIATED(CS%Grid)) call MOM_ERROR(FATAL,'ODA_CS ensemble horizontal grid not associated')
    if (.not. ASSOCIATED(CS%GV)) call MOM_ERROR(FATAL,'ODA_CS ensemble vertical grid not associated')

    call set_current_pelist(CS%filter_pelist)
    id_oda_prior = mpp_clock_id('(ODA setting prior)')
    call mpp_clock_begin(id_oda_prior)

    !! switch to global pelist
    if(is_root_pe()) print *, 'Setting prior'

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

    !! switch back to ensemble member pelist
    call mpp_clock_end(id_oda_prior)
    call set_current_pelist(CS%ensemble_pelist(CS%ensemble_id,:))

    return

  end subroutine set_prior_tracer

  !> Returns posterior adjustments or full state
  !!Note that only those PEs associated with an ensemble member receive data
  subroutine get_posterior_tracer(Time, CS, G, GV, h, tv, increment)
    type(time_type), intent(in) :: Time !< the current model time
    type(ODA_CS), pointer :: CS !< ocean DA control structure
    type(ocean_grid_type), pointer :: G !< domain and grid information for ocean model
    type(verticalGrid_type),               intent(in)    :: GV   !< The ocean's vertical grid structure
    real, dimension(:,:,:), pointer :: h    !< Layer thicknesses, in H (usually m or kg m-2)
    type(thermo_var_ptrs), pointer :: tv   !< A structure pointing to various thermodynamic variables

    logical, optional, intent(in) :: increment

    integer :: i, j, m, id_oda_posterior
    logical :: used, get_inc

    ! return if not analysis time (retain pointers for h and tv)
    if (Time < CS%Time) return

    !! switch to global pelist
    call set_current_pelist(CS%filter_pelist)

    id_oda_posterior = mpp_clock_id('(ODA getting posterior)')
    call mpp_clock_begin(id_oda_posterior)
    if(is_root_pe()) print *, 'Getting posterior'

    get_inc = .true.
    if(present(increment)) get_inc = increment

    if(get_inc) then
      if( .not. associated(Ocean_increment) ) then
        allocate(Ocean_increment)
        call init_ocean_ensemble(Ocean_increment,CS%Grid,CS%GV,CS%ensemble_size)
      endif
      Ocean_increment%T = CS%Ocean_posterior%T - CS%Ocean_prior%T
      Ocean_increment%S = CS%Ocean_posterior%S - CS%Ocean_prior%S
    endif
    do m=1,CS%ensemble_size
      if(get_inc) then
        call mpp_redistribute(CS%mpp_domain, Ocean_increment%T(:,:,:,m), &
                CS%domains(m)%mpp_domain, CS%tv%T, complete=.true.)
        call mpp_redistribute(CS%mpp_domain, Ocean_increment%S(:,:,:,m), &
                CS%domains(m)%mpp_domain, CS%tv%S, complete=.true.)
      else
        call mpp_redistribute(CS%mpp_domain, CS%Ocean_posterior%T(:,:,:,m), &
                CS%domains(m)%mpp_domain, CS%tv%T, complete=.true.)
        call mpp_redistribute(CS%mpp_domain, CS%Ocean_posterior%S(:,:,:,m), &
                CS%domains(m)%mpp_domain, CS%tv%S, complete=.true.)
      endif

    end do

    tv => CS%tv
    h => CS%h
    !! switch back to ensemble member pelist
    call mpp_clock_end(id_oda_posterior)
    call set_current_pelist(CS%ensemble_pelist(CS%ensemble_id,:))

   end subroutine get_posterior_tracer

  subroutine oda(Time, CS)
    type(time_type), intent(in) :: Time
    type(oda_CS), intent(inout) :: CS

    integer :: i, j
    integer :: m
    integer :: yr, mon, day, hr, min, sec

    if ( Time >= CS%Time ) then

      !! switch to global pelist
      call set_current_pelist(CS%filter_pelist)

      call get_profiles(Time, CS%Profiles, CS%CProfiles)
#ifdef ENABLE_ECDA
      call ensemble_filter(CS%Ocean_prior, CS%Ocean_posterior, CS%CProfiles, CS%kdroot, CS%mpp_domain, CS%oda_grid)
#endif

      !! switch back to ensemble member pelist
      call set_current_pelist(CS%ensemble_pelist(CS%ensemble_id,:))

    end if

    return
  end subroutine oda

  subroutine oda_end(CS)
    type(ODA_CS), intent(inout) :: CS

  end subroutine oda_end

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
    allocate(CS%T(is:ie,js:je,nk,ens_size))
    allocate(CS%S(is:ie,js:je,nk,ens_size))
    !allocate(CS%SSH(is:ie,js:je,ens_size))
    !allocate(CS%id_t(ens_size));CS%id_t(:)=-1
    !allocate(CS%id_s(ens_size));CS%id_s(:)=-1
    !allocate(CS%U(is:ie,js:je,nk,ens_size))
    !allocate(CS%V(is:ie,js:je,nk,ens_size))
    !allocate(CS%id_u(ens_size));CS%id_u(:)=-1
    !allocate(CS%id_v(ens_size));CS%id_v(:)=-1
    !allocate(CS%id_ssh(ens_size));CS%id_ssh(:)=-1

    return
  end subroutine init_ocean_ensemble

  subroutine set_analysis_time(Time,CS)
    type(time_type), intent(in) :: Time
    type(ODA_CS), pointer, intent(inout) :: CS

    integer :: yr, mon, day, hr, min, sec

    if (Time >= CS%Time) then
      CS%Time=increment_time(CS%Time,CS%assim_frequency*3600)

      call get_date(Time, yr, mon, day, hr, min, sec)
      if(pe() .eq. mpp_root_pe()) print *, 'Model Time: ', yr, mon, day, hr, min, sec
      call get_date(CS%time, yr, mon, day, hr, min, sec)
      if(pe() .eq. mpp_root_pe()) print *, 'Assimilation Time: ', yr, mon, day, hr, min, sec
    endif
    if (CS%Time < Time) then
        call MOM_error(FATAL, " set_analysis_time: " // &
             "assimilation interval appears to be shorter than " // &
             "the model timestep")
    endif
    
    return

  end subroutine set_analysis_time

  subroutine save_obs_diff(filename,CS)
    character(len=*), intent(in) :: filename
    type(ODA_CS), pointer, intent(in) :: CS

    integer :: fid ! profile file handle
    type(ocean_profile_type), pointer :: Prof=>NULL()

    fid = open_profile_file(trim(filename), nvar=2, thread=MPP_SINGLE, fset=MPP_SINGLE)
    Prof=>CS%CProfiles

    !! switch to global pelist
    !call set_current_pelist(CS%filter_pelist)

    do while (associated(Prof))
      call write_profile(fid,Prof)
      Prof=>Prof%cnext
    enddo
    call close_profile_file(fid)

    !! switch back to ensemble member pelist
    !call set_current_pelist(CS%ensemble_pelist(CS%ensemble_id,:))

    return
  end subroutine save_obs_diff

  subroutine apply_oda_tracer_increments(dt,G,tv,h,CS)
    real, intent(in) :: dt ! the tracer timestep (seconds)
    type(ocean_grid_type),    intent(in)    :: G      !< ocean grid structure
    type(thermo_var_ptrs),    intent(inout) :: tv     !< A structure pointing to various thermodynamic variables
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  &
                              intent(in)    :: h      !< layer thickness (m or kg/m2)
    type(ODA_CS),              intent(inout) :: CS     !< the data assimilation structure

  end subroutine apply_oda_tracer_increments
end module MOM_oda_driver_mod
