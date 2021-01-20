!> Interfaces for MOM6 ensembles and data assimilation.
module MOM_oda_driver_mod

! This file is part of MOM6. see LICENSE.md for the license.

! MOM infrastructure
use MOM_coms, only : PE_here, num_PEs
use MOM_coms, only : set_PElist, set_rootPE, Get_PElist, broadcast
use MOM_domains, only : domain2d, global_field, get_domain_extent
use MOM_domains, only : pass_var, redistribute_array, broadcast_domain
use MOM_diag_mediator, only : register_diag_field, diag_axis_init, post_data
use MOM_ensemble_manager, only : get_ensemble_id, get_ensemble_size
use MOM_ensemble_manager, only : get_ensemble_pelist, get_ensemble_filter_pelist
use MOM_error_handler, only : stdout, stdlog, MOM_error
use MOM_io, only : SINGLE_FILE
use MOM_time_manager, only : time_type, real_to_time, get_date
use MOM_time_manager, only : operator(+), operator(>=), operator(/=)
use MOM_time_manager, only : operator(==), operator(<)

! ODA Modules
use ocean_da_types_mod, only : grid_type, ocean_profile_type, ocean_control_struct
use ocean_da_core_mod, only : ocean_da_core_init, get_profiles
#ifdef ENABLE_ECDA
use eakf_oda_mod, only : ensemble_filter
#endif
use write_ocean_obs_mod, only : open_profile_file
use write_ocean_obs_mod, only : write_profile,close_profile_file
use kdtree, only : kd_root !# JEDI
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

!> A structure with a pointer to a domain2d, to allow for the creation of arrays of pointers.
type :: ptr_mpp_domain
  type(domain2d), pointer :: mpp_domain => NULL() !< pointer to a domain2d
end type ptr_mpp_domain

!> Control structure that contains a transpose of the ocean state across ensemble members.
type, public :: ODA_CS ; private
  type(ocean_control_struct), pointer :: Ocean_prior=> NULL() !< ensemble ocean prior states in DA space
  type(ocean_control_struct), pointer :: Ocean_posterior=> NULL() !< ensemble ocean posterior states
                                                                  !! or increments to prior in DA space
  integer :: nk !< number of vertical layers used for DA
  type(ocean_grid_type), pointer :: Grid => NULL() !< MOM6 grid type and decomposition for the DA
  type(MOM_domain_type), pointer, dimension(:) :: domains => NULL() !< Pointer to mpp_domain objects
                                                                       !! for ensemble members
  type(verticalGrid_type), pointer :: GV => NULL() !< vertical grid for DA
  type(unit_scale_type), pointer :: &
    US => NULL()    !< structure containing various unit conversion factors for DA

  type(domain2d), pointer :: mpp_domain => NULL() !< Pointer to a mpp domain object for DA
  type(grid_type), pointer :: oda_grid !< local tracer grid
  real, pointer, dimension(:,:,:) :: h => NULL() !<layer thicknesses [H ~> m or kg m-2] for DA
  type(thermo_var_ptrs), pointer :: tv => NULL() !< pointer to thermodynamic variables
  integer :: ni          !< global i-direction grid size
  integer :: nj          !< global j-direction grid size
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
  type(kd_root), pointer :: kdroot => NULL() !< A structure for storing nearest neighbors
  type(ALE_CS), pointer :: ALE_CS=>NULL() !< ALE control structure for DA
  logical :: use_ALE_algorithm !< true is using ALE remapping
  type(regridding_CS) :: regridCS !< ALE control structure for regridding
  type(remapping_CS) :: remapCS !< ALE control structure for remapping
  type(time_type) :: Time !< Current Analysis time
  type(diag_ctrl) :: diag_cs !<Diagnostics control structure
end type ODA_CS


!>@{  DA parameters
integer, parameter :: NO_ASSIM = 0, OI_ASSIM=1, EAKF_ASSIM=2
!>@}

contains

!> initialize First_guess (prior) and Analysis grid
!! information for all ensemble members
subroutine init_oda(Time, G, GV, CS)

  type(time_type), intent(in) :: Time !< The current model time.
  type(ocean_grid_type), pointer :: G !< domain and grid information for ocean model
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  type(ODA_CS), pointer, intent(inout) :: CS  !< The DA control structure

! Local variables
  type(thermo_var_ptrs) :: tv_dummy
  type(dyn_horgrid_type), pointer :: dG=> NULL()
  type(hor_index_type), pointer :: HI=> NULL()
  type(directories) :: dirs

  type(grid_type), pointer :: T_grid !< global tracer grid
  real, dimension(:,:), allocatable :: global2D, global2D_old
  real, dimension(:), allocatable :: lon1D, lat1D, glon1D, glat1D
  type(param_file_type) :: PF
  integer :: n, m, k, i, j, nk
  integer :: is,ie,js,je,isd,ied,jsd,jed
  integer :: isg,ieg,jsg,jeg
  integer :: idg_offset, jdg_offset
  integer :: stdout_unit
  character(len=32) :: assim_method
  integer :: npes_pm, ens_info(6), ni, nj
  character(len=128) :: mesg
  character(len=32) :: fldnam
  character(len=30) :: coord_mode
  character(len=200) :: inputdir, basin_file
  logical :: reentrant_x, reentrant_y, tripolar_N, symmetric

  if (associated(CS)) call MOM_error(FATAL, 'Calling oda_init with associated control structure')
  allocate(CS)
! Use ens1 parameters , this could be changed at a later time
! if it were desirable to have alternate parameters, e.g. for the grid
! for the analysis
  call get_MOM_input(PF,dirs,ensemble_num=0)

  call unit_scaling_init(PF, CS%US)

  call get_param(PF, "MOM", "ASSIM_METHOD", assim_method,  &
       "String which determines the data assimilation method "//&
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
       "Use tripolar connectivity at the northern edge of the "//&
       "domain.  With TRIPOLAR_N, NIGLOBAL must be even.", &
       default=.false.)
  call get_param(PF,"MOM", "NIGLOBAL", CS%ni, &
       "The total number of thickness grid points in the "//&
       "x-direction in the physical domain.")
  call get_param(PF,"MOM", "NJGLOBAL", CS%nj, &
       "The total number of thickness grid points in the "//&
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
  CS%domains(CS%ensemble_id)%mpp_domain => G%Domain%mpp_domain
  do n=1,CS%ensemble_size
    if (.not. associated(CS%domains(n)%mpp_domain)) allocate(CS%domains(n)%mpp_domain)
    call set_rootPE(CS%ensemble_pelist(n,1))
    call broadcast_domain(CS%domains(n)%mpp_domain)
  enddo
  call set_rootPE(CS%filter_pelist(1))
  allocate(CS%Grid)
  ! params NIHALO_ODA, NJHALO_ODA set the DA halo size
  call MOM_domains_init(CS%Grid%Domain,PF,param_suffix='_ODA')
  allocate(HI)
  call hor_index_init(CS%Grid%Domain, HI, PF, &
                      local_indexing=.false.)  ! Use global indexing for DA
  call verticalGridInit( PF, CS%GV, CS%US )
  allocate(dG)
  call create_dyn_horgrid(dG, HI)
  call clone_MOM_domain(CS%Grid%Domain, dG%Domain,symmetric=.false.)
  call set_grid_metrics(dG,PF)
  call MOM_initialize_topography(dg%bathyT,dG%max_depth,dG,PF)
  call MOM_initialize_coord(CS%GV, CS%US, PF, .false., &
           dirs%output_directory, tv_dummy, dG%max_depth)
  call ALE_init(PF, CS%GV, CS%US, dG%max_depth, CS%ALE_CS)
  call MOM_grid_init(CS%Grid, PF, global_indexing=.true.)
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
  allocate(CS%tv)

  call get_param(PF, 'oda_driver', "REGRIDDING_COORDINATE_MODE", coord_mode, &
       "Coordinate mode for vertical regridding.", &
       default="ZSTAR", fail_if_missing=.false.)
  call initialize_regridding(CS%regridCS, CS%GV, CS%US, dG%max_depth,PF,'oda_driver',coord_mode,'','')
  call initialize_remapping(CS%remapCS,'PLM')
  call set_regrid_params(CS%regridCS, min_thickness=0.)
  ! breaking with the MOM6 convention and using global indices
  call get_domain_extent(G%Domain,is,ie,js,je,isd,ied,jsd,jed,&
                         isg,ieg,jsg,jeg,idg_offset,jdg_offset,symmetric)
  isd=isd+idg_offset; ied=ied+idg_offset
  jsd=jsd+jdg_offset; jed=jed+jdg_offset
  !call mpp_get_data_domain(G%Domain%mpp_domain,isd,ied,jsd,jed)
  if (.not. associated(CS%h)) then
    allocate(CS%h(isd:ied,jsd:jed,CS%GV%ke)); CS%h(:,:,:)=0.0
    ! assign thicknesses
    call ALE_initThicknessToCoord(CS%ALE_CS,G,CS%GV,CS%h)
  endif
  allocate(CS%tv%T(isd:ied,jsd:jed,CS%GV%ke)); CS%tv%T(:,:,:)=0.0
  allocate(CS%tv%S(isd:ied,jsd:jed,CS%GV%ke)); CS%tv%S(:,:,:)=0.0
  call set_axes_info(CS%Grid, CS%GV, CS%US, PF, CS%diag_cs, set_vertical=.true.)
  ! get domain extents for the analysis grid and use global indexing
  !call get_domain_extent(CS%Grid%Domain,is,ie,js,je,isd,ied,jsd,jed,&
  !                       isg,ieg,jsg,jeg,idg_offset,jdg_offset,symmetric)
  !isd=isd+idg_offset; ied=ied+idg_offset
  !jsd=jsd+jdg_offset; jed=jed+jdg_offset
  !call mpp_get_data_domain(CS%mpp_domain,isd,ied,jsd,jed)
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
  call global_field(CS%mpp_domain, CS%Grid%geolonT, T_grid%x)
  call global_field(CS%mpp_domain, CS%Grid%geolatT, T_grid%y)
  call global_field(CS%mpp_domain, CS%oda_grid%basin_mask, T_grid%basin_mask)
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
    call global_field(G%Domain%mpp_domain, CS%h(:,:,k), global2D)
    do i=1,CS%ni ; do j=1,CS%nj
      if ( global2D(i,j) > 1 ) then
        T_grid%mask(i,j,k) = 1.0
      endif
    enddo ; enddo
    if (k == 1) then
      T_grid%z(:,:,k) = global2D/2
    else
      T_grid%z(:,:,k) = T_grid%z(:,:,k-1) + (global2D + global2D_old)/2
    endif
    global2D_old = global2D
  enddo

  call ocean_da_core_init(CS%mpp_domain, T_grid, CS%Profiles, Time)

  CS%Time=Time
  !! switch back to ensemble member pelist
  call set_PElist(CS%ensemble_pelist(CS%ensemble_id,:))
end subroutine init_oda

!> Copy ensemble member tracers to ensemble vector.
subroutine set_prior_tracer(Time, G, GV, h, tv, CS)
  type(time_type), intent(in)    :: Time !< The current model time
  type(ocean_grid_type), pointer :: G !< domain and grid information for ocean model
  type(verticalGrid_type),               intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                 intent(in) :: tv   !< A structure pointing to various thermodynamic variables

  type(ODA_CS), pointer :: CS !< ocean DA control structure
  real, dimension(SZI_(G),SZJ_(G),CS%nk) :: T, S
  type(ocean_grid_type), pointer :: Grid=>NULL()
  integer :: i,j, m, n, ss
  integer :: is, ie, js, je
  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: isg, ieg, jsg, jeg, idg_offset, jdg_offset
  integer :: id
  logical :: used, symmetric

  ! return if not time for analysis
  if (Time < CS%Time) return

  if (.not. associated(CS%Grid)) call MOM_ERROR(FATAL,'ODA_CS ensemble horizontal grid not associated')
  if (.not. associated(CS%GV)) call MOM_ERROR(FATAL,'ODA_CS ensemble vertical grid not associated')

  !! switch to global pelist
  call set_PElist(CS%filter_pelist)
  call MOM_mesg('Setting prior')

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
subroutine get_posterior_tracer(Time, CS, h, tv, increment)
  type(time_type), intent(in) :: Time !< the current model time
  type(ODA_CS), pointer :: CS !< ocean DA control structure
  real, dimension(:,:,:), pointer :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs), pointer :: tv   !< A structure pointing to various thermodynamic variables
  logical, optional, intent(in) :: increment !< True if returning increment only

  type(ocean_control_struct), pointer :: Ocean_increment=>NULL()
  integer :: i, j, m
  logical :: used, get_inc

  ! return if not analysis time (retain pointers for h and tv)
  if (Time < CS%Time) return


  !! switch to global pelist
  call set_PElist(CS%filter_pelist)
  call MOM_mesg('Getting posterior')

  get_inc = .true.
  if (present(increment)) get_inc = increment

  if (get_inc) then
    allocate(Ocean_increment)
    call init_ocean_ensemble(Ocean_increment,CS%Grid,CS%GV,CS%ensemble_size)
    Ocean_increment%T = CS%Ocean_posterior%T - CS%Ocean_prior%T
    Ocean_increment%S = CS%Ocean_posterior%S - CS%Ocean_prior%S
  endif
  do m=1,CS%ensemble_size
    if (get_inc) then
      call redistribute_array(CS%mpp_domain, Ocean_increment%T(:,:,:,m),&
           CS%domains(m)%mpp_domain, CS%tv%T, complete=.true.)
      call redistribute_array(CS%mpp_domain, Ocean_increment%S(:,:,:,m),&
           CS%domains(m)%mpp_domain, CS%tv%S, complete=.true.)
    else
      call redistribute_array(CS%mpp_domain, CS%Ocean_posterior%T(:,:,:,m),&
           CS%domains(m)%mpp_domain, CS%tv%T, complete=.true.)
      call redistribute_array(CS%mpp_domain, CS%Ocean_posterior%S(:,:,:,m),&
           CS%domains(m)%mpp_domain, CS%tv%S, complete=.true.)
    endif
  enddo

  tv => CS%tv
  h => CS%h
  !! switch back to ensemble member pelist
  call set_PElist(CS%ensemble_pelist(CS%ensemble_id,:))

end subroutine get_posterior_tracer

!> Gather observations and call ODA routines
subroutine oda(Time, CS)
  type(time_type), intent(in) :: Time !< the current model time
  type(oda_CS), intent(inout) :: CS !< the ocean DA control structure

  integer :: i, j
  integer :: m
  integer :: yr, mon, day, hr, min, sec

  if ( Time >= CS%Time ) then

    !! switch to global pelist
    call set_PElist(CS%filter_pelist)

    call get_profiles(Time, CS%Profiles, CS%CProfiles)
#ifdef ENABLE_ECDA
    call ensemble_filter(CS%Ocean_prior, CS%Ocean_posterior, CS%CProfiles, CS%kdroot, CS%mpp_domain, CS%oda_grid)
#endif

    !! switch back to ensemble member pelist
    call set_PElist(CS%ensemble_pelist(CS%ensemble_id,:))

  endif

  return
end subroutine oda

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

  integer :: n,is,ie,js,je,nk

  nk=GV%ke
  is=Grid%isd;ie=Grid%ied
  js=Grid%jsd;je=Grid%jed
  CS%ensemble_size=ens_size
  allocate(CS%T(is:ie,js:je,nk,ens_size))
  allocate(CS%S(is:ie,js:je,nk,ens_size))
  allocate(CS%SSH(is:ie,js:je,ens_size))
!  allocate(CS%id_t(ens_size));CS%id_t(:)=-1
!  allocate(CS%id_s(ens_size));CS%id_s(:)=-1
!  allocate(CS%U(is:ie,js:je,nk,ens_size))
!  allocate(CS%V(is:ie,js:je,nk,ens_size))
!  allocate(CS%id_u(ens_size));CS%id_u(:)=-1
!  allocate(CS%id_v(ens_size));CS%id_v(:)=-1
!  allocate(CS%id_ssh(ens_size));CS%id_ssh(:)=-1

  return
end subroutine init_ocean_ensemble

!> Set the next analysis time
subroutine set_analysis_time(Time,CS)
  type(time_type), intent(in) :: Time !< the current model time
  type(ODA_CS), pointer, intent(inout) :: CS !< the DA control structure

  character(len=160) :: mesg  ! The text of an error message
  integer :: yr, mon, day, hr, min, sec

  if (Time >= CS%Time) then
    ! increment the analysis time to the next step converting to seconds
    CS%Time = CS%Time + real_to_time(CS%US%T_to_s*(CS%assim_frequency*3600.))

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
subroutine save_obs_diff(filename,CS)
  character(len=*), intent(in) :: filename !< name of output file
  type(ODA_CS), pointer, intent(in) :: CS  !< pointer to DA control structure

  integer :: fid ! profile file handle
  type(ocean_profile_type), pointer :: Prof=>NULL()

  fid = open_profile_file(trim(filename), nvar=2, thread=SINGLE_FILE, fset=SINGLE_FILE)
  Prof=>CS%CProfiles

  !! switch to global pelist
  !call set_PElist(CS%filter_pelist)

  do while (associated(Prof))
    call write_profile(fid,Prof)
    Prof=>Prof%cnext
  enddo
  call close_profile_file(fid)

  !! switch back to ensemble member pelist
  !call set_PElist(CS%ensemble_pelist(CS%ensemble_id,:))

  return
end subroutine save_obs_diff


!> Apply increments to tracers
subroutine apply_oda_tracer_increments(dt, G, GV, tv, h, CS)
  real,                     intent(in)    :: dt !< The tracer timestep [s]
  type(ocean_grid_type),    intent(in)    :: G  !< ocean grid structure
  type(verticalGrid_type),  intent(in)    :: GV !< The ocean's vertical grid structure
  type(thermo_var_ptrs),    intent(inout) :: tv !< A structure pointing to various thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h  !< layer thickness [H ~> m or kg m-2]
  type(ODA_CS),             intent(inout) :: CS !< the data assimilation structure

end subroutine apply_oda_tracer_increments

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
