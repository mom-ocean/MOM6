module MOM_offline_main
  
implicit none

type, public :: offline_transport_CS

  !> Pointers to relevant fields from the main MOM control structure
  type(ALE_CS),                  pointer :: ALE_CSp                => NULL()
  type(diabatic_CS),             pointer :: diabatic_CSp           => NULL()
  type(diag_ctrl),               pointer :: diag                   => NULL()
  type(ocean_OBC_type),          pointer :: OBC                    => NULL()
  type(tracer_advect_CS),        pointer :: tracer_adv_CSp         => NULL()
  type(tracer_flow_control_CS),  pointer :: tracer_flow_CSp        => NULL()
  type(tracer_registry_type),    pointer :: tracer_Reg             => NULL()
  type(thermo_var_ptrs)          pointer :: tv                     => NULL()
  
  !> Variables related to reading in fields from online run
  integer :: start_index  ! Timelevel to start
  integer :: numtime      ! How many timelevels in the input fields
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
  logical :: print_adv_offline ! Prints out some updates each advection sub interation
  logical :: skip_diffusion    ! Skips horizontal diffusion of tracers
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

end type offline_transport_CS

contains

subroutine offline_advection_ale(fluxes, state, Time_start, time_interval, CS, h_pre, h_end, uhtr, vhtr)
  type(forcing),    intent(inout)    :: fluxes        !< pointers to forcing fields
  type(surface),    intent(inout)    :: state         !< surface ocean state
  type(time_type),  intent(in)       :: Time_start    !< starting time of a segment, as a time type
  real,             intent(in)       :: time_interval !< time interval
  type(MOM_control_struct), pointer  :: CS            !< control structure from initialize_MOM
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)),  intent(in)     :: h_pre !< layer thicknesses before advection
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)),  intent(in)     :: h_end!< target layer thicknesses
  real, dimension(SZIB_(CS%G),SZJ_(CS%G),SZK_(CS%G)), intent(inout)  :: uhtr  !< Zonal mass transport
  real, dimension(SZI_(CS%G),SZJB_(CS%G),SZK_(CS%G)), intent(inout)  :: vhtr  !< Meridional mass transport
  
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
  real :: dt_offline

  ! Variables used to keep track of layer thicknesses at various points in the code    
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)) :: &    
      h_new, &        
      h_end, &
      h_vol, &
      h_pre
  ! Fields for eta_diff diagnostic
  real, dimension(SZI_(CS%G),SZJ_(CS%G))         :: eta_pre, eta_end    
  integer                                        :: niter, iter
  real                                           :: Inum_iter, dt_iter
  logical                                        :: converged
  integer :: i, j, k, m, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: isv, iev, jsv, jev ! The valid range of the indices.
  integer :: IsdB, IedB, JsdB, JedB
  logical :: z_first, x_before_y
  integer :: niter_vert

  ! Fail out if offline_tracer_mode is not true
  if (.not.CS%offline_tracer_mode) call MOM_error(FATAL,"OFFLINE_TRACER_MODE=False when calling step_tracers")

  ! Grid-related pointer assignments
  G => CS%G
  GV => CS%GV

  ! Initialize some shorthand variables from other structures
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  dt_offline = CS%offline_CSp%dt_offline
  evap_CFL_limit = CS%offline_CSp%evap_CFL_limit
  minimum_forcing_depth = CS%offline_CSp%minimum_forcing_depth
  niter_vert = CEILING(dt_offline/CS%offline_CSp%dt_offline_vertical)

  niter = CS%offline_CSp%num_off_iter
  Inum_iter = 1./real(niter)
  dt_iter = dt_offline*Inum_iter

  ! Initialize working arrays
  h_pre(:,:,:) = GV%Angstrom
  h_new(:,:,:) = GV%Angstrom
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

  ! Convert flux rates into explicit mass/height of freshwater flux. Also note, that
  ! fluxes are halved because diabatic processes are split before and after advection

  ! Copy over the horizontal mass fluxes from the total mass fluxes
  do k=1,nz ; do j=jsd,jed ; do i=isdB,iedB
    uhtr_sub(i,j,k) = uhtr(i,j,k)
  enddo ; enddo ; enddo
  do k=1,nz ; do j=jsdB,jedB ; do i=isd,ied
    vhtr_sub(i,j,k) = vhtr(i,j,k)
  enddo ; enddo ; enddo
  call pass_vector(uhtr_sub,vhtr_sub,G%Domain)
  
  if(CS%debug) then
    call uchksum(uhtr_sub,"uhtr_sub before transport",G%HI)
    call vchksum(vhtr_sub,"vhtr_sub before transport",G%HI)
    call hchksum(h_pre,"h_pre before transport",G%HI)
  endif

  ! This loop does essentially a flux-limited, nonlinear advection scheme until all mass fluxes
  ! are used. ALE is done after the horizontal advection.
  do iter=1,CS%offline_CSp%num_off_iter

    do k=1,nz ; do j=jsd,jed ; do i=isd,ied
      h_vol(i,j,k) = h_pre(i,j,k)*G%areaT(i,j)
    enddo ; enddo ; enddo

    call advect_tracer(h_pre, uhtr_sub, vhtr_sub, CS%OBC, dt_iter, G, GV, &
        CS%tracer_adv_CSp, CS%tracer_Reg, h_vol, max_iter_in=1, &
        uhr_out=uhtr, vhr_out=vhtr, h_out=h_new, x_first_in=x_before_y)
        
    ! Switch the direction every iteration
    x_before_y = .not. x_before_y

    ! Update the new layer thicknesses after one round of advection has happened
    do k=1,nz ; do j=jsd,jed ; do i=isd,ied
      h_pre(i,j,k) = h_new(i,j,k)/G%areaT(i,j)
    enddo ; enddo ; enddo
    call pass_var(h_pre, G%Domain)

    if(CS%debug) then
      call hchksum(h_pre,"h_pre before ALE",G%HI)
    endif
    
    ! Do ALE remapping/regridding to allow for more advection to occur in the next iteration
    call cpu_clock_begin(id_clock_ALE)
    call ALE_main_offline(G, GV, h_pre, CS%tv, &
        CS%tracer_Reg, CS%ALE_CSp, CS%offline_CSp%dt_offline)
    call cpu_clock_end(id_clock_ALE)

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
    if(CS%offline_CSp%print_adv_offline .and. is_root_pe()) &
        print *, "Remaining transport: u", sum_u, "v", sum_v

    if(sum_u+sum_v==0.0) then
      if(is_root_pe()) print *, "Converged after iteration", iter
      converged = .true.
      exit
    else
      converged=.false.
    endif
  enddo               

  ! This diagnostic can be used to identify which grid points did not converge within
  ! the specified number of advection sub iterations
  if(CS%offline_CSp%id_eta_diff>0) then
    eta_pre(:,:) = 0.0
    eta_end(:,:) = 0.0
    do k=1,nz ; do j=jsd,jed ; do i=isd,ied
      if(h_pre(i,j,k)>GV%Angstrom) eta_pre(i,j) = eta_pre(i,j)+h_pre(i,j,k)
      if(h_end(i,j,k)>GV%Angstrom) eta_end(i,j) = eta_end(i,j)+h_end(i,j,k)
    enddo ; enddo; enddo

    call post_data(CS%offline_CSp%id_eta_diff,eta_pre-eta_end,CS%diag)

  endif

  if (.not. converged) then

    do k=1,nz ; do j=jsd,jed ; do i=isd,ied
      h_vol(i,j,k) = h_pre(i,j,k)*G%areaT(i,j)
    enddo ; enddo ; enddo
    call pass_var(h_vol,G%Domain)

    if (CS%debug) then
      call hchksum(h_pre,"h_pre after before redistribute",G%HI)
      call uchksum(uhtr_sub,"uhtr_sub before redistribute",G%HI)
      call vchksum(vhtr_sub,"vhtr_sub before redistribute",G%HI)
    endif

    ! These are used to find out how much will be redistributed in this routine
    if (CS%offline_CSp%id_h_redist>0) call post_data(CS%offline_CSp%id_h_redist, h_pre, CS%diag)
    if (CS%offline_CSp%id_uhr_redist>0) call post_data(CS%offline_CSp%id_uhr_redist, uhtr, CS%diag)
    if (CS%offline_CSp%id_vhr_redist>0) call post_data(CS%offline_CSp%id_vhr_redist, vhtr, CS%diag)

    ! Determine how to actually redistribute the remaining fluxes
    select case (CS%offline_CSp%redistribute_method)
      case ('barotropic')
        if (x_before_y) then
          call distribute_residual_uh_barotropic(G, GV, h_pre, uhtr_sub)
          call distribute_residual_vh_barotropic(G, GV, h_pre, vhtr_sub)
        else
          call distribute_residual_vh_barotropic(G, GV, h_pre, vhtr_sub)
          call distribute_residual_uh_barotropic(G, GV, h_pre, uhtr_sub)
        endif 
        call advect_tracer(h_pre, uhtr_sub, vhtr_sub, CS%OBC, dt_iter, G, GV, &
            CS%tracer_adv_CSp, CS%tracer_Reg, h_vol, max_iter_in=1, &
            uhr_out=uhtr, vhr_out=vhtr, h_out=h_new, x_first_in=x_before_y)

      case ('upwards')
        if (x_before_y) then
          call distribute_residual_uh_upwards(G, GV, h_pre, uhtr_sub)
          call distribute_residual_vh_upwards(G, GV, h_pre, vhtr_sub)
        else
          call distribute_residual_vh_upwards(G, GV, h_pre, vhtr_sub)
          call distribute_residual_uh_upwards(G, GV, h_pre, uhtr_sub)
        endif
        call advect_tracer(h_pre, uhtr_sub, vhtr_sub, CS%OBC, dt_iter, G, GV, &
              CS%tracer_adv_CSp, CS%tracer_Reg, h_vol, max_iter_in=1, &
              uhr_out=uhtr, vhr_out=vhtr, h_out=h_new, x_first_in=x_before_y)
      case ('none')
        call MOM_error(WARNING,"Offline advection did not converge")

      case default
        call MOM_error(FATAL,"Unrecognized REDISTRIBUTE_METHOD")
    end select

    if (CS%debug) then
      call hchksum(h_pre,"h_pre after after redistribute",G%HI)
      call uchksum(uhtr_sub,"uhtr_sub after redistribute",G%HI)
      call vchksum(vhtr_sub,"vhtr_sub after redistribute",G%HI)
    endif

    do k=1,nz ; do j=jsd,jed ; do i=isd,ied
      h_pre(i,j,k) = h_new(i,j,k)/G%areaT(i,j)
    enddo ; enddo ; enddo
    call pass_var(h_pre,G%Domain)

  endif

  ! Call ALE one last time to make sure that tracers are remapped onto the layer thicknesses
  ! stored from the forward run
  call cpu_clock_begin(id_clock_ALE)
  call ALE_offline_tracer_final( G, GV, h_pre, h_end, CS%tracer_Reg, CS%ALE_CSp)
  call cpu_clock_end(id_clock_ALE)        
    
end subroutine offline_advection_ale

subroutine offline_diabatic_ale(fluxes, state, Time_start, time_interval, CS, h_pre, eatr, ebtr)

  type(forcing),    intent(inout)    :: fluxes        !< pointers to forcing fields
  type(surface),    intent(inout)    :: state         !< surface ocean state
  type(time_type),  intent(in)       :: Time_start    !< starting time of a segment, as a time type
  real,             intent(in)       :: time_interval !< time interval
  type(MOM_control_struct), pointer  :: CS            !< control structure from initialize_MOM
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)), intent(inout) :: h_pre
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)), intent(inout) :: eatr !< Entrainment from layer above
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%G)), intent(inout) :: ebtr !< Entrainment from layer below
  
  if (associated(CS%diabatic_CSp%optics)) &
    call set_opacity(CS%diabatic_CSp%optics, fluxes, G, GV, CS%diabatic_CSp%opacity_CSp)

  call offline_add_diurnal_SW(fluxes, G, Time_start, Time_end)
  ! Note that first two arguments are identical, because in ALE mode, there is no change in layer thickness
  ! because of eatr and ebtr
  call call_tracer_column_fns(h_pre, h_pre, eatr, ebtr, &
      fluxes, CS%offline_CSp%dt_offline_vertical, G, GV, CS%tv, &
      CS%diabatic_CSp%optics, CS%tracer_flow_CSp, CS%debug, &
      evap_CFL_limit=evap_CFL_limit, &
      minimum_forcing_depth=minimum_forcing_depth)
  call applyTracerBoundaryFluxesInOut(G, GV, zero_3dh, CS%offline_CSp%dt_offline_vertical, fluxes, h_pre, &
      evap_CFL_limit, minimum_forcing_depth)

end subroutine offline_diabatic_ale

subroutine offline_advection_layer(fluxes, state, Time_start, time_interval, CS, h_pre, eatr, ebtr, uhtr, vhtr)
  type(forcing),    intent(inout)    :: fluxes        !< pointers to forcing fields
  type(surface),    intent(inout)    :: state         !< surface ocean state
  type(time_type),  intent(in)       :: Time_start    !< starting time of a segment, as a time type
  real,             intent(in)       :: time_interval !< time interval
  type(MOM_control_struct), pointer  :: CS            !< control structure from initialize_MOM
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
  ! Zonal diffusive transport
  real, dimension(SZIB_(CS%G),SZJ_(CS%G))              :: khdt_x
  ! Meridional mass transports
  real, dimension(SZI_(CS%G),SZJB_(CS%G),SZK_(CS%G))   :: vhtr_sub
  ! Meridional diffusive transports
  real, dimension(SZI_(CS%G),SZJB_(CS%G))              :: khdt_y

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
      h_vol, &
      h_pre
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
  integer :: niter_vert
  
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  
  do iter=1,CS%offline_CSp%num_off_iter

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
        CS%offline_CSp%max_off_cfl)

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

end module MOM_offline_main
