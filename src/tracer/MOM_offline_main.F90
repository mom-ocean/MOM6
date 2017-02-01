module MOM_offline_main

use mpp_domains_mod,          only : CENTER, CORNER, NORTH, EAST

use MOM_ALE,                  only : ALE_CS, ALE_main_offline, ALE_offline_tracer_final
use MOM_debugging,            only : hchksum, uchksum, vchksum
use MOM_cpu_clock,            only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_diabatic_aux,         only : diabatic_aux_CS
use MOM_diabatic_driver,      only : diabatic_CS
use MOM_diabatic_aux,         only : tridiagTS
use MOM_diag_mediator,        only : diag_ctrl, post_data, register_diag_field
use MOM_domains,              only : sum_across_PEs, pass_var, pass_vector
use MOM_error_handler,        only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler,        only : callTree_enter, callTree_leave
use MOM_file_parser,          only : read_param, get_param, log_version, param_file_type
use MOM_forcing_type,         only : forcing
use MOM_grid,                 only : ocean_grid_type
use MOM_io,                   only : read_data
use MOM_offline_aux,          only : update_h_horizontal_flux, update_h_vertical_flux, limit_mass_flux_3d
use MOM_offline_aux,          only : distribute_residual_uh_barotropic, distribute_residual_vh_barotropic
use MOM_offline_aux,          only : distribute_residual_uh_upwards, distribute_residual_vh_upwards
use MOM_offline_aux,          only : offline_add_diurnal_sw
use MOM_opacity,              only : set_opacity
use MOM_open_boundary,        only : ocean_OBC_type
use MOM_time_manager,         only : time_type
use MOM_tracer_advect,        only : tracer_advect_CS, advect_tracer
use MOM_tracer_diabatic,      only : applyTracerBoundaryFluxesInOut
use MOM_tracer_flow_control,  only : tracer_flow_control_CS, call_tracer_column_fns
use MOM_tracer_registry,      only : tracer_registry_type
use MOM_variables,            only : thermo_var_ptrs
use MOM_verticalGrid,         only : verticalGrid_type

implicit none ; private

#include "MOM_memory.h"
#include "version_variable.h"

type, public :: offline_transport_CS

  !> Pointers to relevant fields from the main MOM control structure
  type(ALE_CS),                  pointer :: ALE_CSp                => NULL()
  type(diabatic_CS),             pointer :: diabatic_CSp           => NULL()
  type(diag_ctrl),               pointer :: diag                   => NULL()
  type(ocean_OBC_type),          pointer :: OBC                    => NULL()
  type(tracer_advect_CS),        pointer :: tracer_adv_CSp         => NULL()
  type(tracer_flow_control_CS),  pointer :: tracer_flow_CSp        => NULL()
  type(tracer_registry_type),    pointer :: tracer_Reg             => NULL()
  type(thermo_var_ptrs),         pointer :: tv                     => NULL()
  type(ocean_grid_type),         pointer :: G                      => NULL()
  type(verticalGrid_type),       pointer :: GV                     => NULL()

  !> Variables related to reading in fields from online run
  integer :: start_index  ! Timelevel to start
  integer :: iter_no      ! Timelevel to start
  integer :: numtime      ! How many timelevels in the input fields
  integer :: accumulated_time ! Length of time accumulatedi n the current offline interval
  integer :: &            ! Index of each of the variables to be read in
    ridx_sum = -1, &      ! Separate indices for each variabile if they are
    ridx_snap = -1        ! setoff from each other in time
  character(len=200) :: offlinedir  ! Directory where offline fields are stored
  character(len=200) :: & !         ! Names of input files
    snap_file,  &
    sum_file,   &
    mean_file
  character(len=20)  :: redistribute_method
  logical :: fields_are_offset ! True if the time-averaged fields and snapshot fields are
                               ! offset by one time level
  logical :: x_before_y        ! Which horizontal direction is advected first
  logical :: print_adv_offline ! Prints out some updates each advection sub interation
  logical :: skip_diffusion    ! Skips horizontal diffusion of tracers
  logical :: read_sw           ! Read in averaged values for shortwave radiation
  logical :: diurnal_sw        ! Adds a synthetic diurnal cycle on shortwave radiation
  logical :: debug
  !> Variables controlling some of the numerical considerations of offline transport
  integer           ::  num_off_iter
  real              ::  dt_offline ! Timestep used for offline tracers
  real              ::  dt_offline_vertical ! Timestep used for calls to tracer vertical physics
  real              ::  max_off_cfl=0.5 ! Hardcoded for now, only used in non-ALE mode
  real              ::  evap_CFL_limit, minimum_forcing_depth

  !> Diagnostic manager IDs for some fields that may be of interest when doing offline transport
  integer :: &
    id_uhr = -1, &
    id_vhr = -1, &
    id_ear = -1, &
    id_ebr = -1, &
    id_hr = -1,  &
    id_uhr_redist = -1, &
    id_vhr_redist = -1, &
    id_h_redist = -1, &
    id_eta_diff = -1

  !> Variables that may need to be stored between calls to step_MOM
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_)       :: uhtr
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_)       :: vhtr

  ! Fields at T-point
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: &
      eatr,     &  ! Amount of fluid entrained from the layer above within
                   ! one time step  (m for Bouss, kg/m^2 for non-Bouss)
      ebtr         ! Amount of fluid entrained from the layer below within
                   ! one time step  (m for Bouss, kg/m^2 for non-Bouss)
  ! Work arrays for temperature and salinity
  ! Arrays for temperature and salinity
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: &
      temp_mean, salt_mean, &
      h_end
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
      netMassIn, netMassOut


end type offline_transport_CS

public offline_advection_ale
public offline_redistribute_residual
public offline_diabatic_ale
public offline_advection_layer
public transport_by_files
public offline_transport_init
public register_diags_offline_transport

contains

subroutine offline_advection_ale(fluxes, Time_start, time_interval, CS, id_clock_ALE,&
      h_pre, uhtr, vhtr, converged)
  type(forcing),    intent(inout)    :: fluxes        !< pointers to forcing fields
  type(time_type),  intent(in)       :: Time_start    !< starting time of a segment, as a time type
  real,             intent(in)       :: time_interval !< time interval
  type(offline_transport_CS), pointer  :: CS            !< control structure from initialize_MOM
  integer :: id_clock_ALE !< ID for the ALE clock
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)),  intent(inout)  :: h_pre !< layer thicknesses before advection
  real, dimension(SZIB_(CS%G),SZJ_(CS%G),SZK_(CS%G)), intent(inout)  :: uhtr  !< Zonal mass transport
  real, dimension(SZI_(CS%G),SZJB_(CS%G),SZK_(CS%G)), intent(inout)  :: vhtr  !< Meridional mass transport
  logical,                                            intent(  out)  :: converged

  ! Local pointers
  type(ocean_grid_type),      pointer :: G  => NULL() ! Pointer to a structure containing
                                                      ! metrics and related information
  type(verticalGrid_type),    pointer :: GV => NULL() ! Pointer to structure containing information
                                                      ! about the vertical grid
  ! Work arrays for mass transports
  real, dimension(SZIB_(CS%G),SZJ_(CS%G),SZK_(CS%G))   :: uhtr_sub
  ! Meridional mass transports
  real, dimension(SZI_(CS%G),SZJB_(CS%G),SZK_(CS%G))   :: vhtr_sub

  real :: sum_abs_fluxes, sum_u, sum_v  ! Used to keep track of how close to convergence we are

  ! Variables used to keep track of layer thicknesses at various points in the code
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)) :: &
      h_new, &
      h_vol
  ! Fields for eta_diff diagnostic
  real, dimension(SZI_(CS%G),SZJ_(CS%G))         :: eta_pre, eta_end
  integer                                        :: niter, iter
  real                                           :: Inum_iter
  integer :: i, j, k, m, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: isv, iev, jsv, jev ! The valid range of the indices.
  integer :: IsdB, IedB, JsdB, JedB
  logical :: z_first, x_before_y
  real :: evap_CFL_limit, minimum_forcing_depth, dt_iter, dt_offline

  ! Grid-related pointer assignments
  G => CS%G
  GV => CS%GV

  x_before_y = CS%x_before_y

  ! Initialize some shorthand variables from other structures
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  dt_offline = CS%dt_offline
  evap_CFL_limit = CS%evap_CFL_limit
  minimum_forcing_depth = CS%minimum_forcing_depth

  niter = CS%num_off_iter
  Inum_iter = 1./real(niter)
  dt_iter = dt_offline*Inum_iter

  ! Initialize working arrays
  h_new(:,:,:) = 0.0
  h_vol(:,:,:) = 0.0
  uhtr_sub(:,:,:) = 0.0
  vhtr_sub(:,:,:) = 0.0

  ! Initialize logicals
  converged = .false.

  ! Tracers are transported using the stored mass fluxes. Where possible, operators are Strang-split around
  ! the call to
  ! 1)  Using the layer thicknesses and tracer concentrations from the previous timestep,
  !     half of the accumulated vertical mixing (eatr and ebtr) is applied in the call to tracer_column_fns.
  !     For tracers whose source/sink terms need dt, this value is set to 1/2 dt_offline
  ! 2)  Half of the accumulated surface freshwater fluxes are applied
  !! START ITERATION
  ! 3)  Accumulated mass fluxes are used to do horizontal transport. The number of iterations used in
  !     advect_tracer is limited to 2 (e.g x->y->x->y). The remaining mass fluxes are stored for later use
  !     and resulting layer thicknesses fed into the next step
  ! 4)  Tracers and the h-grid are regridded and remapped in a call to ALE. This allows for layers which might
  !     'vanish' because of horizontal mass transport to be 'reinflated'
  ! 5)  Check that transport is done if the remaining mass fluxes equals 0 or if the max number of iterations
  !     has been reached
  !! END ITERATION
  ! 6)  Repeat steps 1 and 2
  ! 7)  Force a remapping to the stored layer thicknesses that correspond to the snapshot of the online model
  ! 8)  Reset T/S and h to their stored snapshotted values to prevent model drift

  ! Copy over the horizontal mass fluxes from the total mass fluxes
  do k=1,nz ; do j=jsd,jed ; do i=isdB,iedB
    uhtr_sub(i,j,k) = uhtr(i,j,k)
  enddo ; enddo ; enddo
  do k=1,nz ; do j=jsdB,jedB ; do i=isd,ied
    vhtr_sub(i,j,k) = vhtr(i,j,k)
  enddo ; enddo ; enddo
  call pass_vector(uhtr_sub,vhtr_sub,G%Domain)

  if(CS%debug) then
    call hchksum(h_pre,"h_pre before transport",G%HI)
    call uchksum(uhtr_sub,"uhtr_sub before transport",G%HI)
    call vchksum(vhtr_sub,"vhtr_sub before transport",G%HI)
  endif

  ! This loop does essentially a flux-limited, nonlinear advection scheme until all mass fluxes
  ! are used. ALE is done after the horizontal advection.
  do iter=1,CS%num_off_iter

    do k=1,nz ; do j=jsd,jed ; do i=isd,ied
      h_vol(i,j,k) = h_pre(i,j,k)*G%areaT(i,j)
    enddo ; enddo ; enddo

    call advect_tracer(h_pre, uhtr_sub, vhtr_sub, CS%OBC, CS%dt_offline, G, GV, &
        CS%tracer_adv_CSp, CS%tracer_Reg, h_vol, max_iter_in=2, &
        uhr_out=uhtr, vhr_out=vhtr, h_out=h_new, x_first_in=x_before_y)

    ! Switch the direction every iteration
    x_before_y = .not. x_before_y

    ! Update the new layer thicknesses after one round of advection has happened
    do k=1,nz ; do j=jsd,jed ; do i=isd,ied
      h_pre(i,j,k) = h_new(i,j,k)/G%areaT(i,j)
    enddo ; enddo ; enddo

    if(CS%debug) then
      call hchksum(h_pre,"h_pre before ALE",G%HI)
    endif

    ! Do ALE remapping/regridding to allow for more advection to occur in the next iteration
    call cpu_clock_begin(id_clock_ALE)
    call ALE_main_offline(G, GV, h_pre, CS%tv, &
        CS%tracer_Reg, CS%ALE_CSp, CS%dt_offline)
    call cpu_clock_end(id_clock_ALE)

    call pass_var(h_pre, G%Domain)
    if(CS%debug) then
      call hchksum(h_pre,"h_pre after ALE",G%HI)
    endif

    ! Copy over remaining mass transports
    do k=1,nz ; do j=jsd,jed ; do i=isdB,iedB
      uhtr_sub(i,j,k) = uhtr(i,j,k)
    enddo ; enddo ; enddo
    do k=1,nz ; do j=jsdB,jedB ; do i=isd,ied
      vhtr_sub(i,j,k) = vhtr(i,j,k)
    enddo ; enddo ; enddo
    call pass_vector(uhtr_sub,vhtr_sub,G%Domain)

    sum_u = 0.0
    do k=1,nz; do j=js,je ; do i=is-1,ie
      sum_u = sum_u + abs(uhtr_sub(i,j,k))
    enddo; enddo; enddo
    sum_v = 0.0
    do k=1,nz; do j=js-1,je; do i=is,ie
      sum_v = sum_v + abs(vhtr_sub(i,j,k))
    enddo; enddo ; enddo

    call sum_across_PEs(sum_u)
    call sum_across_PEs(sum_v)
    if(CS%print_adv_offline .and. is_root_pe()) &
        print *, "Remaining transport: u", sum_u, "v", sum_v

    if(sum_u+sum_v==0.0) then
      if(is_root_pe()) print *, "Converged after iteration", iter
      converged = .true.
      exit
    else
      converged=.false.
    endif
  enddo

  ! Make sure that uhtr and vhtr halos are updated
  call pass_vector(uhtr,vhtr,G%Domain)

  if(CS%debug) then
    call hchksum(h_pre,"h after offline_advection_ale",G%HI)
    call uchksum(uhtr,"uhtr after offline_advection_ale",G%HI)
    call vchksum(vhtr,"vhtr after offline_advection_ale",G%HI)
  endif

end subroutine offline_advection_ale

subroutine offline_redistribute_residual(CS, h_pre, h_end, uhtr, vhtr, converged)
  type(offline_transport_CS), pointer  :: CS            !< control structure from initialize_MOM

  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)),  intent(inout)  :: h_pre !< layer thicknesses before advection
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)),  intent(in)     :: h_end !< target layer thicknesses
  real, dimension(SZIB_(CS%G),SZJ_(CS%G),SZK_(CS%G)), intent(inout)  :: uhtr  !< Zonal mass transport
  real, dimension(SZI_(CS%G),SZJB_(CS%G),SZK_(CS%G)), intent(inout)  :: vhtr  !< Meridional mass transport
  logical,                                            intent(in)     :: converged

  type(ocean_grid_type),      pointer :: G  => NULL() ! Pointer to a structure containing
                                                      ! metrics and related information
  type(verticalGrid_type),    pointer :: GV => NULL() ! Pointer to structure containing information
                                                      ! about the vertical grid
  logical :: x_before_y
  ! Variables used to keep track of layer thicknesses at various points in the code
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)) :: &
      h_new, &
      h_vol

  ! Used to calculate the eta_diff diagnostic
  real, dimension(SZI_(CS%G),SZJ_(CS%G)) :: &
      eta_pre, &
      eta_end
  real, dimension(SZIB_(CS%G),SZJ_(CS%G),SZK_(CS%G)) :: uhr  !< Zonal mass transport
  real, dimension(SZI_(CS%G),SZJB_(CS%G),SZK_(CS%G)) :: vhr  !< Meridional mass transport

  integer :: i, j, k, m, is, ie, js, je, isd, ied, jsd, jed, nz

  ! Assign grid pointers
  G  => CS%G
  GV => CS%GV

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  x_before_y = CS%x_before_y

  uhr(:,:,:) = uhtr
  vhr(:,:,:) = vhtr

  do k=1,nz ; do j=jsd,jed ; do i=isd,ied
    h_vol(i,j,k) = h_pre(i,j,k)*G%areaT(i,j)
  enddo ; enddo ; enddo
  call pass_var(h_vol,G%Domain)
  call pass_vector(uhtr,vhtr,G%Domain)

  if (CS%debug) then
    call hchksum(h_pre,"h_pre before upwards redistribute",G%HI, haloshift = 1)
    call hchksum(h_vol,"h_vol before upwards redistribute",G%HI)
    call uchksum(uhtr,"uhtr before upwards redistribute",G%HI)
    call vchksum(vhtr,"vhtr before upwards redistribute",G%HI)
  endif

  ! These are used to find out how much will be redistributed in this routine
  if (CS%id_h_redist>0) call post_data(CS%id_h_redist, h_pre, CS%diag)
  if (CS%id_uhr_redist>0) call post_data(CS%id_uhr_redist, uhtr, CS%diag)
  if (CS%id_vhr_redist>0) call post_data(CS%id_vhr_redist, vhtr, CS%diag)

  ! First try to distribute the residual upwards and advect
  if(x_before_y) then
    call distribute_residual_uh_upwards(G, GV, h_pre, uhtr)
    call pass_var(h_pre,G%Domain)
    call distribute_residual_vh_upwards(G, GV, h_pre, vhtr)
    call pass_var(h_pre,G%Domain)
  else
    call distribute_residual_vh_upwards(G, GV, h_pre, vhtr)
    call pass_var(h_pre,G%Domain)
    call distribute_residual_uh_upwards(G, GV, h_pre, uhtr)
    call pass_var(h_pre,G%Domain)
  endif
  call pass_vector(uhtr,vhtr,G%Domain)

  if (CS%debug) then
    call hchksum(h_vol,"h_vol after upwards redistribute",G%HI,haloshift = 1)
    call uchksum(uhtr,"uh after upwards redistribute",G%HI)
    call vchksum(vhtr,"vh after upwards redistribute",G%HI)
  endif
  call advect_tracer(h_pre, uhtr, vhtr, CS%OBC, CS%dt_offline, G, GV, &
        CS%tracer_adv_CSp, CS%tracer_Reg, h_vol, max_iter_in=2, &
        h_out=h_new, uhr_out=uhr, vhr_out=vhr, x_first_in=x_before_y)

  do k=1,nz ; do j=jsd,jed ; do i=isd,ied
    h_vol(i,j,k) = h_new(i,j,k)
    h_new(i,j,k) = h_vol(i,j,k)/G%areaT(i,j)
  enddo ; enddo ; enddo

  call pass_var(h_vol,G%Domain)
  call pass_var(h_new,G%Domain)
  call pass_vector(uhr,vhr,G%Domain)

  if (CS%debug) then
    call hchksum(h_vol,"h_vol before barotropic redistribute",G%HI)
    call hchksum(h_new,"h_new before barotropic redistribute",G%HI)
    call uchksum(uhr,"uhr before barotropic redistribute",G%HI)
    call vchksum(vhr,"vhr before barotropic redistribute",G%HI)
  endif

  ! Then check if there's any transport left and if so, distribute it equally
  ! throughout the rest of the water column
  if(x_before_y) then
    call distribute_residual_uh_barotropic(G, GV, h_new, uhr)
    call pass_var(h_new,G%Domain)
    call distribute_residual_vh_barotropic(G, GV, h_new, vhr)
    call pass_var(h_new,G%Domain)
  else
    call distribute_residual_vh_barotropic(G, GV, h_new, vhr)
    call pass_var(h_new,G%Domain)
    call distribute_residual_uh_barotropic(G, GV, h_new, uhr)
    call pass_var(h_new,G%Domain)
  endif
  call pass_vector(uhr,vhr,G%Domain)
  if (CS%debug) then
    call hchksum(h_vol,"h_vol after barotropic redistribute",G%HI)
    call uchksum(uhr,"uhr after barotropic redistribute",G%HI)
    call vchksum(vhr,"vhr after barotropic redistribute",G%HI)
  endif

  call advect_tracer(h_new, uhr, vhr, CS%OBC, CS%dt_offline, G, GV, &
      CS%tracer_adv_CSp, CS%tracer_Reg, h_vol, max_iter_in=2, &
      h_out=h_pre, uhr_out=uhtr, vhr_out=vhtr, x_first_in=x_before_y)

  if (CS%debug) then
    call hchksum(h_pre,"h_pre after advection barotropic redistribute",G%HI)
    call uchksum(uhtr,"uhtr after advection barotropic redistribute",G%HI)
    call vchksum(vhtr,"vhtr after advection barotropic redistribute",G%HI)

    do k=1,nz ; do j=js,je ; do i=is-1,ie
      if( abs(uhtr(i,j,k))>0.0 ) print *, "Remaining uhtr i, j, k: ", uhr(i,j,k), i, j, k
    enddo ; enddo ; enddo

    do k=1,nz ; do j=js-1,je ; do i=is,ie
      if( abs(vhtr(i,j,k))>0.0 ) print *, "Remaining vhtr i, j, k: ", vhr(i,j,k), i, j, k
    enddo ; enddo ; enddo

  endif

  do k=1,nz ; do j=jsd,jed ; do i=isd,ied
    h_pre(i,j,k) = h_pre(i,j,k)/G%areaT(i,j)
  enddo ; enddo ; enddo
  call pass_var(h_pre,G%Domain)

  ! This diagnostic can be used to identify which grid points did not converge within
  ! the specified number of advection sub iterations
  if(CS%id_eta_diff>0) then
    eta_pre(:,:) = 0.0
    eta_end(:,:) = 0.0
    do k=1,nz ; do j=jsd,jed ; do i=isd,ied
      if(h_pre(i,j,k)>GV%Angstrom) eta_pre(i,j) = eta_pre(i,j)+h_pre(i,j,k)
      if(h_end(i,j,k)>GV%Angstrom) eta_end(i,j) = eta_end(i,j)+h_end(i,j,k)
    enddo ; enddo; enddo

    call post_data(CS%id_eta_diff,eta_pre-eta_end,CS%diag)

  endif

  if(CS%id_uhr>0) call post_data(CS%id_uhr,uhtr,CS%diag)
  if(CS%id_vhr>0) call post_data(CS%id_vhr,vhtr,CS%diag)

end subroutine offline_redistribute_residual

subroutine offline_diabatic_ale(fluxes, Time_start, Time_end, dt, CS, h_pre, eatr, ebtr)

  type(forcing),    intent(inout)      :: fluxes        !< pointers to forcing fields
  type(time_type),  intent(in)         :: Time_start    !< starting time of a segment, as a time type
  type(time_type),  intent(in)         :: Time_end      !< time interval
  real,             intent(in)         :: dt            !< Time step to be used for column functions
  type(offline_transport_CS), pointer  :: CS            !< control structure from initialize_MOM
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)), intent(inout) :: h_pre
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)), intent(in)    :: eatr !< Entrainment from layer above
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)), intent(in)    :: ebtr !< Entrainment from layer below
  real, dimension(SZI_(CS%G),SZJ_(CS%G))    :: sw, sw_vis, sw_nir !< Save old value of shortwave radiation

  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)) :: zero_3dh
  zero_3dh(:,:,:) = 0.0

  if(CS%debug) then
    call hchksum(h_pre,"h_pre before offline_diabatic_ale",CS%G%HI)
    call hchksum(eatr,"eatr before offline_diabatic_ale",CS%G%HI)
    call hchksum(ebtr,"ebtr before offline_diabatic_ale",CS%G%HI)
  endif

  if(CS%diurnal_SW .and. CS%read_sw) then
    sw(:,:) = fluxes%sw
    sw_vis(:,:) = fluxes%sw_vis_dir
    sw_nir(:,:) = fluxes%sw_nir_dir
    call offline_add_diurnal_SW(fluxes, CS%G, Time_start, Time_end)
  endif

  if (associated(CS%diabatic_CSp%optics)) &
    call set_opacity(CS%diabatic_CSp%optics, fluxes, CS%G, CS%GV, CS%diabatic_CSp%opacity_CSp)

  ! Note that first two arguments are identical, because in ALE mode, there is no change in layer thickness
  ! because of eatr and ebtr

  call call_tracer_column_fns(h_pre, h_pre, eatr, ebtr, &
      fluxes, dt, CS%G, CS%GV, CS%tv, &
      CS%diabatic_CSp%optics, CS%tracer_flow_CSp, CS%debug, &
      evap_CFL_limit=CS%evap_CFL_limit, &
      minimum_forcing_depth=CS%minimum_forcing_depth)
  ! This next line is called to calculate new layer thicknesses based on the freshwater fluxes
  call applyTracerBoundaryFluxesInOut(CS%G, CS%GV, zero_3dh, dt, fluxes, h_pre, &
      CS%evap_CFL_limit, CS%minimum_forcing_depth)

  ! The mixing of T/S is probably not appropriate because time averages are being used
  ! call triDiagTS(CS%G, CS%GV, CS%G%isc, CS%G%iec, CS%G%jsc, CS%G%jec, h_pre, eatr, ebtr, CS%tv%T, CS%tv%S)

  if(CS%diurnal_SW .and. CS%read_sw) then
    fluxes%sw(:,:) = sw
    fluxes%sw_vis_dir(:,:) = sw_vis
    fluxes%sw_nir_dir(:,:) = sw_nir
  endif

  if(CS%debug) then
    call hchksum(h_pre,"h_pre after offline_diabatic_ale",CS%G%HI)
    call hchksum(eatr,"eatr after offline_diabatic_ale",CS%G%HI)
    call hchksum(ebtr,"ebtr after offline_diabatic_ale",CS%G%HI)
  endif

end subroutine offline_diabatic_ale

subroutine offline_advection_layer(fluxes, Time_start, time_interval, CS, h_pre, eatr, ebtr, uhtr, vhtr)
  type(forcing),    intent(inout)    :: fluxes        !< pointers to forcing fields
  type(time_type),  intent(in)       :: Time_start    !< starting time of a segment, as a time type
  real,             intent(in)       :: time_interval !< time interval
  type(offline_transport_CS), pointer  :: CS            !< control structure from initialize_MOM
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)),  intent(inout) :: h_pre !< layer thicknesses before advection
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)),  intent(inout) :: eatr !< Entrainment from layer above
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)),  intent(inout) :: ebtr !< Entrainment from layer below
  real, dimension(SZIB_(CS%G),SZJ_(CS%G),SZK_(CS%G)), intent(inout) :: uhtr  !< Zonal mass transport
  real, dimension(SZI_(CS%G),SZJB_(CS%G),SZK_(CS%G)), intent(inout) :: vhtr  !< Meridional mass transport
  ! Local pointers
  type(ocean_grid_type),      pointer :: G  => NULL() ! Pointer to a structure containing
                                                      ! metrics and related information
  type(verticalGrid_type),    pointer :: GV => NULL() ! Pointer to structure containing information
                                                      ! about the vertical grid
  ! Zonal mass transports
  real, dimension(SZIB_(CS%G),SZJ_(CS%G),SZK_(CS%G))   :: uhtr_sub
  ! Meridional mass transports
  real, dimension(SZI_(CS%G),SZJB_(CS%G),SZK_(CS%G))   :: vhtr_sub

  real :: sum_abs_fluxes, sum_u, sum_v  ! Used to keep track of how close to convergence we are
  real :: dt_offline ! Shorthand variables from offline CS

  ! Local variables
  ! Vertical diffusion related variables
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)) :: &
      eatr_sub, &
      ebtr_sub
  ! Variables used to keep track of layer thicknesses at various points in the code
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)) :: &
      h_new, &
      h_vol
  ! Work arrays for temperature and salinity
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)) :: &
      temp_old, salt_old, &
      temp_mean, salt_mean, &
      zero_3dh     !
  integer                                        :: niter, iter
  real                                           :: Inum_iter, dt_iter
  logical                                        :: converged
  integer :: i, j, k, m, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: isv, iev, jsv, jev ! The valid range of the indices.
  integer :: IsdB, IedB, JsdB, JedB
  logical :: z_first, x_before_y

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  do iter=1,CS%num_off_iter

    do k = 1, nz ; do j=js-1,je+1 ; do i=is-1,ie+1
      eatr_sub(i,j,k) = eatr(i,j,k)
      ebtr_sub(i,j,k) = ebtr(i,j,k)
    enddo; enddo ; enddo

    do k = 1, nz ; do j=js-1,je+1 ; do i=is-2,ie+1
      uhtr_sub(I,j,k) = uhtr(I,j,k)
    enddo; enddo ; enddo

    do k = 1, nz ; do j=js-2,je+1 ; do i=is-1,ie+1
      vhtr_sub(i,J,k) = vhtr(i,J,k)
    enddo; enddo ; enddo


    ! Calculate 3d mass transports to be used in this iteration
    call limit_mass_flux_3d(G, GV, uhtr_sub, vhtr_sub, eatr_sub, ebtr_sub, h_pre, &
        CS%max_off_cfl)

    if (z_first) then
      ! First do vertical advection
      call update_h_vertical_flux(G, GV, eatr_sub, ebtr_sub, h_pre, h_new)
      call call_tracer_column_fns(h_pre, h_new, eatr_sub, ebtr_sub, &
          fluxes, dt_iter, G, GV, CS%tv, CS%diabatic_CSp%optics, CS%tracer_flow_CSp, CS%debug)
      ! We are now done with the vertical mass transports, so now h_new is h_sub
      do k = 1, nz ; do j=js-1,je+1 ; do i=is-1,ie+1
        h_pre(i,j,k) = h_new(i,j,k)
      enddo ; enddo ; enddo
      call pass_var(h_pre,G%Domain)

      ! Second zonal and meridional advection
      call update_h_horizontal_flux(G, GV, uhtr_sub, vhtr_sub, h_pre, h_new)
      do k = 1, nz ; do i = is-1, ie+1 ; do j=js-1, je+1
        h_vol(i,j,k) = h_pre(i,j,k)*G%areaT(i,j)
      enddo; enddo; enddo
      call advect_tracer(h_pre, uhtr_sub, vhtr_sub, CS%OBC, dt_iter, G, GV, &
          CS%tracer_adv_CSp, CS%tracer_Reg, h_vol, max_iter_in=30, x_first_in=x_before_y)

      ! Done with horizontal so now h_pre should be h_new
      do k = 1, nz ; do i=is-1,ie+1 ; do j=js-1,je+1
          h_pre(i,j,k) = h_new(i,j,k)
      enddo ; enddo ; enddo

    endif

    if (.not. z_first) then

      ! First zonal and meridional advection
      call update_h_horizontal_flux(G, GV, uhtr_sub, vhtr_sub, h_pre, h_new)
      do k = 1, nz ; do i = is-1, ie+1 ; do j=js-1, je+1
        h_vol(i,j,k) = h_pre(i,j,k)*G%areaT(i,j)
      enddo; enddo; enddo
      call advect_tracer(h_pre, uhtr_sub, vhtr_sub, CS%OBC, dt_iter, G, GV, &
          CS%tracer_adv_CSp, CS%tracer_Reg, h_vol, max_iter_in=30, x_first_in=x_before_y)

      ! Done with horizontal so now h_pre should be h_new
      do k = 1, nz ; do i=is-1,ie+1 ; do j=js-1,je+1
          h_pre(i,j,k) = h_new(i,j,k)
      enddo ; enddo ; enddo

      ! Second vertical advection
      call update_h_vertical_flux(G, GV, eatr_sub, ebtr_sub, h_pre, h_new)
      call call_tracer_column_fns(h_pre, h_new, eatr_sub, ebtr_sub, &
          fluxes, dt_iter, G, GV, CS%tv, CS%diabatic_CSp%optics, CS%tracer_flow_CSp, CS%debug)
      ! We are now done with the vertical mass transports, so now h_new is h_sub
      do k = 1, nz ; do i=is-1,ie+1 ; do j=js-1,je+1
        h_pre(i,j,k) = h_new(i,j,k)
      enddo ; enddo ; enddo

    endif

    ! Update remaining transports
    do k = 1, nz ; do j=js-1,je+1 ; do i=is-1,ie+1
      eatr(i,j,k) = eatr(i,j,k) - eatr_sub(i,j,k)
      ebtr(i,j,k) = ebtr(i,j,k) - ebtr_sub(i,j,k)
    enddo; enddo ; enddo

    do k = 1, nz ; do j=js-1,je+1 ; do i=is-2,ie+1
      uhtr(I,j,k) = uhtr(I,j,k) - uhtr_sub(I,j,k)
    enddo; enddo ; enddo

    do k = 1, nz ; do j=js-2,je+1 ; do i=is-1,ie+1
      vhtr(i,J,k) = vhtr(i,J,k) - vhtr_sub(i,J,k)
    enddo; enddo ; enddo

    call pass_var(eatr,G%Domain)
    call pass_var(ebtr,G%Domain)
    call pass_var(h_pre,G%Domain)
    call pass_vector(uhtr,vhtr,G%Domain)
  !
    ! Calculate how close we are to converging by summing the remaining fluxes at each point
    sum_abs_fluxes = 0.0
    sum_u = 0.0
    sum_v = 0.0
    do k=1,nz; do j=js,je; do i=is,ie
      sum_u = sum_u + abs(uhtr(I-1,j,k))+abs(uhtr(I,j,k))
      sum_v = sum_v + abs(vhtr(i,J-1,k))+abs(vhtr(I,J,k))
      sum_abs_fluxes = sum_abs_fluxes + abs(eatr(i,j,k)) + abs(ebtr(i,j,k)) + abs(uhtr(I-1,j,k)) + &
          abs(uhtr(I,j,k)) + abs(vhtr(i,J-1,k)) + abs(vhtr(i,J,k))
    enddo; enddo; enddo
    call sum_across_PEs(sum_abs_fluxes)

    print *, "Remaining u-flux, v-flux:", sum_u, sum_v
    if (sum_abs_fluxes==0) then
      print *, 'Converged after iteration', iter
      exit
    endif

    ! Switch order of Strang split every iteration
    z_first = .not. z_first
    x_before_y = .not. x_before_y
  end do

end subroutine offline_advection_layer

!> Controls the reading in 3d mass fluxes, diffusive fluxes, and other fields stored
!! in a previous integration of the online model
subroutine transport_by_files(G, GV, CS, h_end, eatr, ebtr, uhtr, vhtr, &
    temp_mean, salt_mean, fluxes, do_ale_in)

  type(ocean_grid_type),                     intent(inout)    :: G
  type(verticalGrid_type),                   intent(inout)    :: GV
  type(offline_transport_CS),                intent(inout)    :: CS
  logical, optional                                           :: do_ale_in

  !! Mandatory variables
  ! Fields at U-points
  !  3D
  ! Zonal mass transports
  real, dimension(SZIB_(CS%G),SZJ_(CS%G),SZK_(CS%G))   :: uhtr
  ! Meridional mass transports
  real, dimension(SZI_(CS%G),SZJB_(CS%G),SZK_(CS%G))   :: vhtr

  ! Vertical diffusion related variables
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)) :: &
    h_end, &
    eatr, ebtr, &
    temp_mean, salt_mean
  type(forcing)                                               :: fluxes
  logical                                                     :: do_ale
  integer :: i, j, k, is, ie, js, je, nz
  real                                                        :: Initer_vert
  do_ale = .false.;
  if (present(do_ale_in) ) do_ale = do_ale_in

  is   = G%isc   ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec ; nz = GV%ke

  call callTree_enter("transport_by_files, MOM_offline_control.F90")

  uhtr(:,:,:) = 0.0
  vhtr(:,:,:) = 0.0
  eatr(:,:,:) = 0.0
  ebtr(:,:,:) = 0.0
  h_end(:,:,:) = 0.0
  temp_mean(:,:,:) = 0.0
  salt_mean(:,:,:) = 0.0

  !! Time-summed fields
  ! U-grid
  call read_data(CS%sum_file, 'uhtr_sum',     uhtr,domain=G%Domain%mpp_domain, &
    timelevel=CS%ridx_sum,position=EAST)
  ! V-grid
  call read_data(CS%sum_file, 'vhtr_sum',     vhtr, domain=G%Domain%mpp_domain, &
    timelevel=CS%ridx_sum,position=NORTH)
  ! T-grid
  call read_data(CS%sum_file, 'ea_sum',   eatr, domain=G%Domain%mpp_domain, &
    timelevel=CS%ridx_sum,position=CENTER)
  call read_data(CS%sum_file, 'eb_sum',   ebtr, domain=G%Domain%mpp_domain, &
    timelevel=CS%ridx_sum,position=CENTER)

  !! Time-averaged fields
  call read_data(CS%mean_file, 'temp',   temp_mean, domain=G%Domain%mpp_domain, &
    timelevel=CS%ridx_sum,position=CENTER)
  call read_data(CS%mean_file, 'salt',   salt_mean, domain=G%Domain%mpp_domain, &
    timelevel=CS%ridx_sum,position=CENTER)

  !! Read snapshot fields (end of time interval timestamp)
  call read_data(CS%snap_file, 'h_end', h_end, domain=G%Domain%mpp_domain, &
    timelevel=CS%ridx_snap,position=CENTER)

  ! This block makes sure that the fluxes control structure, which may not be used in the solo_driver,
  ! contains netMassIn and netMassOut which is necessary for the applyTracerBoundaryFluxesInOut routine
  if (do_ale) then
    if (.not. ASSOCIATED(fluxes%netMassOut)) then
      ALLOC_(fluxes%netMassOut(G%isd:G%ied,G%jsd:G%jed))
    endif
    if (.not. ASSOCIATED(fluxes%netMassIn)) then
      ALLOC_(fluxes%netMassIn(G%isd:G%ied,G%jsd:G%jed))
    endif

    CS%netMassOut(:,:) = 0.0
    CS%netMassIn(:,:) = 0.0
    call read_data(CS%sum_file,'massout_flux_sum',CS%netMassOut, domain=G%Domain%mpp_domain, &
        timelevel=CS%ridx_sum)
    call read_data(CS%sum_file,'massin_flux_sum', CS%netMassIn,  domain=G%Domain%mpp_domain, &
        timelevel=CS%ridx_sum)

    do j=js,je ; do i=is,ie
      if(G%mask2dT(i,j)<1.0) then
        CS%netMassOut(i,j) = 0.0
        CS%netMassIn(i,j) = 0.0
      endif
    enddo ; enddo

  endif

  if(CS%read_sw) then

    ! Shortwave radiation is only needed for offline mode with biogeochemistry.
    ! Need to double check, but set_opacity seems to only need the sum of the diffuse and
    ! direct fluxes in the visible and near-infrared bands. For convenience, we store the
    ! sum of the direct and diffuse fluxes in the 'dir' field and set the 'dif' fields to zero
    call read_data(CS%mean_file,'sw_vis',fluxes%sw_vis_dir, domain=G%Domain%mpp_domain, &
        timelevel=CS%ridx_sum)
    call read_data(CS%mean_file,'sw_nir',fluxes%sw_nir_dir, domain=G%Domain%mpp_domain, &
        timelevel=CS%ridx_sum)
    fluxes%sw_vis_dir(:,:) = fluxes%sw_vis_dir(:,:)*0.5
    fluxes%sw_vis_dif(:,:) = fluxes%sw_vis_dir
    fluxes%sw_nir_dir(:,:) = fluxes%sw_nir_dir(:,:)*0.5
    fluxes%sw_nir_dif(:,:) = fluxes%sw_nir_dir
    fluxes%sw = fluxes%sw_vis_dir + fluxes%sw_vis_dif + fluxes%sw_nir_dir + fluxes%sw_nir_dif
    do j=js,je ; do i=is,ie
      if(G%mask2dT(i,j)<1.0) then
        fluxes%sw(i,j) = 0.0
        fluxes%sw_vis_dir(i,j) = 0.0
        fluxes%sw_nir_dir(i,j) = 0.0
        fluxes%sw_vis_dif(i,j) = 0.0
        fluxes%sw_nir_dif(i,j) = 0.0
      endif
    enddo ; enddo
    call pass_var(fluxes%sw,G%Domain)
    call pass_var(fluxes%sw_vis_dir,G%Domain)
    call pass_var(fluxes%sw_vis_dif,G%Domain)
    call pass_var(fluxes%sw_nir_dir,G%Domain)
    call pass_var(fluxes%sw_nir_dif,G%Domain)
  endif

  ! Apply masks at T, U, and V points
  do k=1,nz ; do j=js,je ; do i=is,ie
    if(G%mask2dT(i,j)<1.0) then
      h_end(i,j,k) = GV%Angstrom
      eatr(i,j,k) = 0.0
      ebtr(i,j,k) = 0.0
      temp_mean(i,j,k) = 0.0
      salt_mean(i,j,k) = 0.0
    endif

  enddo; enddo ; enddo

  do k=1,nz ; do J=js-1,je ; do i=is,ie
    if(G%mask2dCv(i,J)<1.0) then
      vhtr(i,J,k) = 0.0
    endif
  enddo; enddo ; enddo

  do k=1,nz ; do j=js,je ; do I=is-1,ie
    if(G%mask2dCu(I,j)<1.0) then
      uhtr(I,j,k) = 0.0
    endif
  enddo; enddo ; enddo


  !! Make sure all halos have been updated
  ! Vector fields
  call pass_vector(uhtr, vhtr, G%Domain)

  ! Scalar fields
  call pass_var(h_end, G%Domain)
  call pass_var(eatr, G%Domain)
  call pass_var(ebtr, G%Domain)
  call pass_var(temp_mean, G%Domain)
  call pass_var(salt_mean, G%Domain)

  if (do_ale) then
    call pass_var(CS%netMassOut,G%Domain)
    call pass_var(CS%netMassIn,G%Domain)
  endif

  ! Update the read indices
  CS%ridx_snap = next_modulo_time(CS%ridx_snap,CS%numtime)
  CS%ridx_sum = next_modulo_time(CS%ridx_sum,CS%numtime)

  if(CS%debug) then
    call hchksum(h_end,"h_end after transport_by_file",G%HI)
    call hchksum(eatr,"eatr after transport_by_file",G%HI)
    call hchksum(ebtr,"ebtr after transport_by_file",G%HI)
    call uchksum(uhtr,"uhtr after transport_by_file",G%HI)
    call vchksum(vhtr,"vhtr after transport_by_file",G%HI)
  endif

  call callTree_leave("transport_by_file")

end subroutine transport_by_files

!> Calculates the next timelevel to read from the input fields. This allows the 'looping'
!! of the fields
function next_modulo_time(inidx, numtime)
  ! Returns the next time interval to be read
  integer                 :: numtime              ! Number of time levels in input fields
  integer                 :: inidx                ! The current time index

  integer                 :: read_index           ! The index in the input files that corresponds
                                                  ! to the current timestep

  integer                 :: next_modulo_time

  read_index = mod(inidx+1,numtime)
  if (read_index < 0)  read_index = inidx-read_index
  if (read_index == 0) read_index = numtime

  next_modulo_time = read_index

end function next_modulo_time

!> Initialize additional diagnostics required for offline tracer transport
subroutine register_diags_offline_transport(Time, diag, CS)

  type(offline_transport_CS), pointer :: CS         !< control structure for MOM
  type(time_type), intent(in) :: Time               !< current model time
  type(diag_ctrl)             :: diag


  ! U-cell fields
  CS%id_uhr = register_diag_field('ocean_model', 'uhr', diag%axesCuL, Time, &
    'Zonal thickness fluxes remaining at end of timestep', 'kg')
  CS%id_uhr_redist = register_diag_field('ocean_model', 'uhr_redist', diag%axesCuL, Time, &
    'Zonal thickness fluxes to be redistributed vertically', 'kg')

  ! V-cell fields
  CS%id_vhr = register_diag_field('ocean_model', 'vhr', diag%axesCvL, Time, &
    'Meridional thickness fluxes remaining at end of timestep', 'kg')
  CS%id_vhr_redist = register_diag_field('ocean_model', 'vhr_redist', diag%axesCvL, Time, &
    'Meridional thickness to be redistributed vertically', 'kg')

  ! T-cell fields
  CS%id_hr  = register_diag_field('ocean_model', 'hdiff', diag%axesTL, Time, &
    'Difference between the stored and calculated layer thickness', 'm')
  CS%id_ear  = register_diag_field('ocean_model', 'ear', diag%axesTL, Time, &
    'Remaining thickness entrained from above', 'm')
  CS%id_ebr  = register_diag_field('ocean_model', 'ebr', diag%axesTL, Time, &
    'Remaining thickness entrained from below', 'm')
  CS%id_eta_diff = register_diag_field('ocean_model','eta_diff', diag%axesT1, Time, &
    'Difference in total water column height from online and offline','m')
  CS%id_h_redist = register_diag_field('ocean_model','h_redist', diag%axesTL, Time, &
    'Layer thicknesses before redistribution of mass fluxes','m')

end subroutine register_diags_offline_transport

!> Initializes the control structure for offline transport and reads in some of the
! run time parameters from MOM_input
subroutine offline_transport_init(param_file, CS, diabatic_aux_CSp, G, GV)

  type(param_file_type),               intent(in)     :: param_file
  type(offline_transport_CS), pointer, intent(inout)  :: CS
  type(diabatic_aux_CS),      pointer, intent(in)     :: diabatic_aux_CSp
  type(ocean_grid_type),      pointer, intent(in)     :: G
  type(verticalGrid_type),    pointer, intent(in)     :: GV

  character(len=40)                               :: mod = "offline_transport"

  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB

  is   = G%isc   ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec ; nz = GV%ke
  isd  = G%isd   ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB  ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call callTree_enter("offline_transport_init, MOM_offline_control.F90")

  if (associated(CS)) then
    call MOM_error(WARNING, "offline_transport_init called with an associated "// &
      "control structure.")
    return
  endif
  allocate(CS)
  call log_version(param_file,mod,version, &
    "This module allows for tracers to be run offline")

  ! Parse MOM_input for offline control
  call get_param(param_file, mod, "OFFLINEDIR", CS%offlinedir, &
    "Input directory where the offline fields can be found", default=" ")
  call get_param(param_file, mod, "OFF_SUM_FILE", CS%sum_file, &
    "Filename where the accumulated fields can be found", default = " ")
  call get_param(param_file, mod, "OFF_SNAP_FILE", CS%snap_file, &
    "Filename where snapshot fields can be found",default=" ")
  call get_param(param_file, mod, "OFF_MEAN_FILE", CS%mean_file, &
    "Filename where averaged fields can be found",default=" ")
  call get_param(param_file, mod, "START_INDEX", CS%start_index, &
    "Which time index to start from", default=1)
  call get_param(param_file, mod, "NUMTIME", CS%numtime, &
    "Number of timelevels in offline input files", default=0)
  call get_param(param_file, mod, "FIELDS_ARE_OFFSET", CS%fields_are_offset, &
    "True if the time-averaged fields and snapshot fields\n"//&
    "are offset by one time level", default=.false.)
  call get_param(param_file, mod, "REDISTRIBUTE_METHOD", CS%redistribute_method, &
    "Redistributes any remaining horizontal fluxes throughout\n"//&
    "the rest of water column. Options are 'barotropic' which\n"//&
    "evenly distributes flux throughout the entire water column,\n"//&
    "'upwards' which adds the maximum of the remaining flux in\n"//&
    "each layer above, and 'none' which does no redistribution", &
    default='barotropic')
  call get_param(param_file, mod, "NUM_OFF_ITER", CS%num_off_iter, &
    "Number of iterations to subdivide the offline tracer advection and diffusion" )
  call get_param(param_file, mod, "DT_OFFLINE", CS%dt_offline, &
    "Length of the offline timestep")
  call get_param(param_file, mod, "DT_OFFLINE_VERTICAL", CS%dt_offline_vertical, &
    "Length of the offline timestep for tracer column sources/sinks")
  call get_param(param_file, mod, "PRINT_ADV_OFFLINE", CS%print_adv_offline, &
    "Print diagnostic output every advection subiteration",default=.false.)
  call get_param(param_file, mod, "SKIP_DIFFUSION_OFFLINE", CS%skip_diffusion, &
    "Do not do horizontal diffusion",default=.false.)
  call get_param(param_file, mod, "READ_SW", CS%read_sw, &
    "Read in shortwave radiation field instead of using values from the coupler"//&
    "when in offline tracer mode",default=.false.)
  call get_param(param_file, mod, "OFFLINE_ADD_DIURNAL_SW", CS%diurnal_sw, &
    "Adds a synthetic diurnal cycle in the same way that the ice model would have"//&
    "when time-averaged fields of shortwave radiation are read in", default=.true.)

  ! Concatenate offline directory and file names
  CS%snap_file = trim(CS%offlinedir)//trim(CS%snap_file)
  CS%mean_file = trim(CS%offlinedir)//trim(CS%mean_file)
  CS%sum_file = trim(CS%offlinedir)//trim(CS%sum_file)

  ! Set the accumulated time to zero
  CS%accumulated_time = 0
  ! Set the starting read index for time-averaged and snapshotted fields
  CS%ridx_sum = CS%start_index
  if(CS%fields_are_offset) CS%ridx_snap = next_modulo_time(CS%start_index,CS%numtime)
  if(.not. CS%fields_are_offset) CS%ridx_snap = CS%start_index

  ! Copy over parameters from other control structures
  if(associated(diabatic_aux_CSp)) then
    CS%evap_CFL_limit = diabatic_aux_CSp%evap_CFL_limit
    CS%minimum_forcing_depth = diabatic_aux_CSp%minimum_forcing_depth
  endif

  ! Grid pointer assignments
  CS%G  => G
  CS%GV => GV

  ! Allocate arrays
  ALLOC_(CS%uhtr(IsdB:IedB,jsd:jed,nz))   ; CS%uhtr(:,:,:) = 0.0
  ALLOC_(CS%vhtr(isd:ied,JsdB:JedB,nz))   ; CS%vhtr(:,:,:) = 0.0
  ALLOC_(CS%eatr(isd:ied,jsd:jed,nz))          ; CS%eatr(:,:,:) = 0.0
  ALLOC_(CS%ebtr(isd:ied,jsd:jed,nz))          ; CS%ebtr(:,:,:) = 0.0
  ALLOC_(CS%temp_mean(isd:ied,jsd:jed,nz))     ; CS%temp_mean(:,:,:) = 0.0
  ALLOC_(CS%salt_mean(isd:ied,jsd:jed,nz))     ; CS%salt_mean(:,:,:) = 0.0
  ALLOC_(CS%h_end(isd:ied,jsd:jed,nz))         ; CS%h_end(:,:,:) = 0.0
  ALLOC_(CS%netMassOut(G%isd:G%ied,G%jsd:G%jed)) ; CS%netMassOut(:,:) = 0.0
  ALLOC_(CS%netMassIn(G%isd:G%ied,G%jsd:G%jed))  ; CS%netMassIn(:,:) = 0.0

  call callTree_leave("offline_transport_init")


end subroutine offline_transport_init

!> \namespace mom_offline_main
!! \section offline_overview Offline Tracer Transport in MOM6
!!  'Offline tracer modeling' uses physical fields (e.g. mass transports and layer thicknesses) saved
!!  from a previous integration of the physical model to transport passive tracers. These fields are
!!  accumulated or averaged over a period of time (in this test case, 1 day) and used to integrate
!!  portions of the MOM6 code base that handle the 3d advection and diffusion of passive tracers.
!!
!!  The distribution of tracers in the ocean modeled offline should not be expected to match an online
!!  simulation. Accumulating transports over more than one online model timestep implicitly assumes
!!  homogeneity over that time period and essentially aliases over processes that occur with higher
!!  frequency. For example, consider the case of a surface boundary layer with a strong diurnal cycle.
!!  An offline simulation with a 1 day timestep, captures the net transport into or out of that layer,
!!  but not the exact cycling. This effective aliasing may also complicate online model configurations
!!  which strongly-eddying regions. In this case, the offline model timestep must be limited to some
!!  fraction of the eddy correlation timescale. Lastly, the nonlinear advection scheme which applies
!!  limited mass-transports over a sequence of iterations means that tracers are not transported along
!!  exactly the same path as they are in the online model.
!!
!!  This capability has currently targeted the Baltic_ALE_z test case, though some work has also been
!!  done with the OM4 1/2 degree configuration. Work is ongoing to develop recommendations and best
!!  practices for investigators seeking to use MOM6 for offline tracer modeling.
!!
!!  \section offline_technical Implementation of offline routine in MOM6
!!
!!  The subroutine step_tracers that coordinates this can be found in MOM.F90 and is only called
!!  using the solo ocean driver. This is to avoid issues with coupling to other climate components
!!  that may be relying on fluxes from the ocean to be coupled more often than the offline time step.
!!  Other routines related to offline tracer modeling can be found in tracers/MOM_offline_control.F90
!!
!!  As can also be seen in the comments for the step_tracers subroutine, an offline time step
!!  comprises the following steps:
!!        -#  Using the layer thicknesses and tracer concentrations from the previous timestep,
!!            half of the accumulated vertical mixing (eatr and ebtr) is applied in the call to
!!            tracer_column_fns.
!!            For tracers whose source/sink terms need dt, this value is set to 1/2 dt_offline
!!        -#  Half of the accumulated surface freshwater fluxes are applied
!!        START ITERATION
!!        -#  Accumulated mass fluxes are used to do horizontal transport. The number of iterations
!!            used in advect_tracer is limited to 2 (e.g x->y->x->y). The remaining mass fluxes are
!!            stored for later use and resulting layer thicknesses fed into the next step
!!        -#  Tracers and the h-grid are regridded and remapped in a call to ALE. This allows for
!!            layers which might 'vanish' because of horizontal mass transport to be 'reinflated'
!!            and essentially allows for the vertical transport of tracers
!!        -#  Check that transport is done if the remaining mass fluxes equals 0 or if the max
!!            number of iterations has been reached
!!        END ITERATION
!!        -#  Repeat steps 1 and 2
!!        -#  Redistribute any residual mass fluxes that remain after the advection iterations
!!            in a barotropic manner, progressively upward through the water column.
!!        -#  Force a remapping to the stored layer thicknesses that correspond to the snapshot of
!!            the online model at the end of an accumulation interval
!!        -#  Reset T/S and h to their stored snapshotted values to prevent model drift
!!
!!  \section  offline_evaluation Evaluating the utility of an offline tracer model
!!  How well an offline tracer model can be used as an alternative to integrating tracers online
!!  with the prognostic model must be evaluated for each application. This efficacy may be related
!!  to the native coordinate of the online model, to the length of the offline timestep, and to the
!!  behavior of the tracer itself.
!!
!!  A framework for formally regression testing the offline capability still needs to be developed.
!!  However, as a simple way of testing whether the offline model is nominally behaving as expected,
!!  the total inventory of the advection test tracers (tr1, tr2, etc.) should be conserved between
!!  time steps except for the last 4 decimal places. As a general guideline, an offline timestep of
!!  5 days or less.
!!
!!  \section offline_parameters Runtime parameters for offline tracers
!!    - OFFLINEDIR:    Input directory where the offline fields can be found
!!    - OFF_SUM_FILE:  Filename where the accumulated fields can be found (e.g. horizontal mass transports)
!!    - OFF_SNAP_FILE: Filename where snapshot fields can be found (e.g. end of timestep layer thickness)
!!    - START_INDEX:   Which timelevel of the input files to read first
!!    - NUMTIME:       How many timelevels to read before 'looping' back to 1
!!    - FIELDS_ARE_OFFSET: True if the time-averaged fields and snapshot fields are offset by one
!!                        time level, probably not needed
!!    -NUM_OFF_ITER:  Maximum number of iterations to do for the nonlinear advection scheme
!!    -REDISTRIBUTE_METHOD: Redistributes any remaining horizontal fluxes throughout the rest of water column.
!!                          Options are 'barotropic' which "evenly distributes flux throughout the entire water
!!                          column,'upwards' which adds the maximum of the remaining flux in each layer above,
!!                          and 'none' which does no redistribution"

end module MOM_offline_main
