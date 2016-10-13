!> Contains routines related to offline transport of tracers
module MOM_offline_transport
! This file is part of MOM6. See LICENSE.md for the license.

  use data_override_mod,    only : data_override_init, data_override
  use MOM_time_manager,     only : time_type
  use MOM_domains,          only : pass_var, pass_vector, To_All
  use MOM_error_handler,    only : callTree_enter, callTree_leave, MOM_error, FATAL, WARNING, is_root_pe
  use MOM_grid,             only : ocean_grid_type
  use MOM_verticalGrid,     only : verticalGrid_type
  use MOM_io,               only : read_data
  use MOM_file_parser,      only : get_param, log_version, param_file_type
  use MOM_diag_mediator,    only : diag_ctrl, register_diag_field
  use mpp_domains_mod,      only : CENTER, CORNER, NORTH, EAST
  use MOM_variables,        only : vertvisc_type
  use MOM_forcing_type,     only : forcing
  use MOM_shortwave_abs,    only : optics_type
  use MOM_diag_mediator,    only : post_data
  use MOM_forcing_type,     only : forcing
  use MOM_diabatic_driver,  only : diabatic_CS

  implicit none

#include <MOM_memory.h>

  type, public :: offline_transport_CS

    !> Variables related to reading in fields from online run
    integer :: start_index  ! Timelevel to start
    integer :: numtime      ! How many timelevels in the input fields
    integer :: &            ! Index of each of the variables to be read in
      ridx_sum = -1, &      ! Separate indices for each variabile if they are
      ridx_snap = -1        ! setoff from each other in time
    character(len=200) :: offlinedir  ! Directory where offline fields are stored
    character(len=200) :: & !         ! Names of input files
      snap_file,  &
      sum_file
    logical :: fields_are_offset ! True if the time-averaged fields and snapshot fields are
                                 ! offset by one time level
    
    !> Variables controlling some of the numerical considerations of offline transport
    integer           ::  num_off_iter
    real              ::  dt_offline ! Timestep used for offline tracers
    real              ::  max_off_cfl=0.5 ! Hardcoded for now, only used in non-ALE mode
    real              ::  evap_CFL_limit, minimum_forcing_depth

    !> Diagnostic manager IDs for use in the online model of additional fields necessary
    !> for offline tracer modeling
    integer :: &
      id_uhtr_preadv = -1, &
      id_vhtr_preadv = -1, &
    !> Diagnostic manager IDs for some fields that may be of interest when doing offline transport  
      id_uhr = -1, &
      id_vhr = -1, &
      id_ear = -1, &
      id_ebr = -1, &
      id_hr = -1

  end type offline_transport_CS

#include "MOM_memory.h"
#include "version_variable.h"
  public offline_transport_init
  public transport_by_files
  public register_diags_offline_transport
  public update_h_horizontal_flux
  public update_h_vertical_flux
  public limit_mass_flux_3d

contains
  !> Controls the reading in 3d mass fluxes, diffusive fluxes, and other fields stored
  !! in a previous integration of the online model
  subroutine transport_by_files(G, GV, CS, h_end, eatr, ebtr, uhtr, vhtr, khdt_x, khdt_y, &
      temp, salt, fluxes, do_ale_in)

    type(ocean_grid_type),                     intent(inout)    :: G
    type(verticalGrid_type),                   intent(inout)    :: GV
    type(offline_transport_CS),                intent(inout)    :: CS
    logical, optional                                           :: do_ale_in

    !! Mandatory variables
    ! Fields at U-points
    !  3D
    real, dimension(SZIB_(G),SZJ_(G),SZK_(G))                   :: uhtr
    !  2D
    real, dimension(SZIB_(G),SZJ_(G))                           :: khdt_x
    ! Fields at V-points
    !  3D
    real, dimension(SZI_(G),SZJB_(G),SZK_(G))                   :: vhtr
    !  2D
    real, dimension(SZI_(G),SZJB_(G))                           :: khdt_y
    ! Fields at T-point
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
      h_end, &
      eatr, ebtr, &
      temp, salt
    type(forcing)                                               :: fluxes
    logical                                                     :: do_ale
    integer :: i, j, k, is, ie, js, je, nz

    do_ale = .false.;
    if (present(do_ale_in) ) do_ale = do_ale_in

    is   = G%isc   ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec ; nz = GV%ke


    call callTree_enter("transport_by_files, MOM_offline_control.F90")


    !! Time-summed fields
    ! U-grid
    call read_data(CS%sum_file, 'uhtr_sum',     uhtr,domain=G%Domain%mpp_domain, &
      timelevel=CS%ridx_sum,position=EAST)
    call read_data(CS%sum_file, 'khdt_x_sum', khdt_x,domain=G%Domain%mpp_domain, &
      timelevel=CS%ridx_sum,position=EAST)
    ! V-grid
    call read_data(CS%sum_file, 'vhtr_sum',     vhtr, domain=G%Domain%mpp_domain, &
      timelevel=CS%ridx_sum,position=NORTH)
    call read_data(CS%sum_file, 'khdt_y_sum', khdt_y, domain=G%Domain%mpp_domain, &
      timelevel=CS%ridx_sum,position=NORTH)
    ! T-grid
    call read_data(CS%sum_file, 'ea_sum',   eatr, domain=G%Domain%mpp_domain, &
      timelevel=CS%ridx_sum,position=CENTER)
    call read_data(CS%sum_file, 'eb_sum',   ebtr, domain=G%Domain%mpp_domain, &
      timelevel=CS%ridx_sum,position=CENTER)

    !! Time-averaged fields
    call read_data(CS%snap_file, 'temp',   temp, domain=G%Domain%mpp_domain, &
      timelevel=CS%ridx_sum,position=CENTER)
    call read_data(CS%snap_file, 'salt',   salt, domain=G%Domain%mpp_domain, &
      timelevel=CS%ridx_sum,position=CENTER)

    !! Read snapshot fields (end of time interval timestamp)
    call read_data(CS%snap_file, 'h_end', h_end, domain=G%Domain%mpp_domain, &
      timelevel=CS%ridx_snap,position=CENTER)

    ! Apply masks at T, U, and V points
    do k=1,nz ; do j=js,je ; do i=is,ie
      if(G%mask2dT(i,j)<1.0) then
        h_end(i,j,k) = GV%Angstrom
        eatr(i,j,k) = 0.0
        ebtr(i,j,k) = 0.0
      endif
    enddo; enddo ; enddo

    do k=1,nz ; do j=js-1,je ; do i=is,ie
      if(G%mask2dCv(i,j)<1.0) then
        khdt_y(i,j) = 0.0
        vhtr(i,j,k) = 0.0
      endif
    enddo; enddo ; enddo

    do k=1,nz ; do j=js,je ; do i=is-1,ie
      if(G%mask2dCu(i,j)<1.0) then
        khdt_x(i,j) = 0.0
        uhtr(i,j,k) = 0.0
      endif
    enddo; enddo ; enddo

    ! This block makes sure that the fluxes control structure, which may not be used in the solo_driver,
    ! contains netMassIn and netMassOut which is necessary for the applyTracerBoundaryFluxesInOut routine
    if (do_ale) then
      if (.not. ASSOCIATED(fluxes%netMassOut)) then
        ALLOCATE(fluxes%netMassOut(G%isd:G%ied,G%jsd:G%jed))
        fluxes%netMassOut(:,:) = 0.0
      endif
      if (.not. ASSOCIATED(fluxes%netMassIn)) then
        ALLOCATE(fluxes%netMassIn(G%isd:G%ied,G%jsd:G%jed))
        fluxes%netMassIn(:,:) = 0.0
      endif
      
      call read_data(CS%sum_file,'massout_flux_sum',fluxes%netMassOut, domain=G%Domain%mpp_domain, &
          timelevel=CS%ridx_snap,position=center)
      call read_data(CS%sum_file,'massin_flux_sum', fluxes%netMassIn,  domain=G%Domain%mpp_domain, &
          timelevel=CS%ridx_snap,position=center)
    endif

    !! Make sure all halos have been updated
    ! Vector fields
    call pass_vector(uhtr, vhtr, G%Domain)
    call pass_vector(khdt_x, khdt_y, G%Domain)

    ! Scalar fields
    call pass_var(h_end, G%Domain)
    call pass_var(eatr, G%Domain)
    call pass_var(ebtr, G%Domain)
    call pass_var(temp, G%Domain)
    call pass_var(salt, G%Domain)

    if (do_ale) then
      call pass_var(fluxes%netMassOut,G%Domain)
      call pass_var(fluxes%netMassIn,G%Domain)
    endif

    ! Update the read indices
    CS%ridx_snap = next_modulo_time(CS%ridx_snap,CS%numtime)
    CS%ridx_sum = next_modulo_time(CS%ridx_sum,CS%numtime)

    call callTree_leave("transport_by_file")

  end subroutine transport_by_files

  !> Initialize additional diagnostics required for offline tracer transport
  subroutine register_diags_offline_transport(Time, diag, CS)

    type(offline_transport_CS), pointer :: CS         !< control structure for MOM
    type(time_type), intent(in) :: Time             !< current model time
    type(diag_ctrl)             :: diag


    ! U-cell fields
    CS%id_uhr = register_diag_field('ocean_model', 'uhr', diag%axesCuL, Time, &
      'Zonal thickness fluxes remaining at end of timestep', 'kg')

    ! V-cell fields
    CS%id_vhr = register_diag_field('ocean_model', 'vhr', diag%axesCvL, Time, &
      'Meridional thickness fluxes remaining at end of timestep', 'kg')

    ! T-cell fields
    CS%id_hr  = register_diag_field('ocean_model', 'hdiff', diag%axesTL, Time, &
      'Difference between the stored and calculated layer thickness', 'm')
    CS%id_ear  = register_diag_field('ocean_model', 'ear', diag%axesTL, Time, &
      'Remaining thickness entrained from above', 'm')
    CS%id_ebr  = register_diag_field('ocean_model', 'ebr', diag%axesTL, Time, &
      'Remaining thickness entrained from below', 'm')

  end subroutine register_diags_offline_transport

  ! Initializes the control structure for offline transport and reads in some of the
  ! run time parameters from MOM_input
  subroutine offline_transport_init(param_file, CS, diabatic_CSp, G, GV)

    type(param_file_type),               intent(in)     :: param_file
    type(offline_transport_CS), pointer, intent(inout)  :: CS
    type(diabatic_CS),          pointer, intent(in)     :: diabatic_CSp
    type(ocean_grid_type),               intent(in)     :: G
    type(verticalGrid_type),             intent(in)     :: GV

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
    call get_param(param_file, mod, "START_INDEX", CS%start_index, &
      "Which time index to start from", default=1)
    call get_param(param_file, mod, "NUMTIME", CS%numtime, &
      "Number of timelevels in offline input files", default=0)
    call get_param(param_file, mod, "FIELDS_ARE_OFFSET", CS%fields_are_offset, &
      "True if the time-averaged fields and snapshot fields are offset by one time level", &
      default=.false.)
    call get_param(param_file, mod, "NUM_OFF_ITER", CS%num_off_iter, &
      "Number of iterations to subdivide the offline tracer advection and diffusion" )
    call get_param(param_file, mod, "DT_OFFLINE", CS%dt_offline, &
      "Length of the offline timestep")

    ! Concatenate offline directory and file names
    CS%snap_file = trim(CS%offlinedir)//trim(CS%snap_file)
    CS%sum_file = trim(CS%offlinedir)//trim(CS%sum_file)

    ! Set the starting read index for time-averaged and snapshotted fields
    CS%ridx_sum = CS%start_index
    if(CS%fields_are_offset) CS%ridx_snap = next_modulo_time(CS%start_index,CS%numtime)
    if(.not. CS%fields_are_offset) CS%ridx_snap = CS%start_index
    
    ! Copy over parameters from other control structures
    CS%evap_CFL_limit = diabatic_CSp%diabatic_aux_CSp%evap_CFL_limit
    CS%minimum_forcing_depth = diabatic_CSp%diabatic_aux_CSp%minimum_forcing_depth

    call callTree_leave("offline_transport_init")

  end subroutine offline_transport_init
  
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

  !> This updates thickness based on the convergence of horizontal mass fluxes
  !! NOTE: Only used in non-ALE mode
  subroutine update_h_horizontal_flux(G, GV, uhtr, vhtr, h_pre, h_new)
    type(ocean_grid_type),    pointer                           :: G
    type(verticalGrid_type),  pointer                           :: GV
    real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)       :: uhtr
    real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)       :: vhtr
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(in)       :: h_pre
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)    :: h_new

    ! Local variables
    integer :: i, j, k, m, is, ie, js, je, nz
    ! Set index-related variables for fields on T-grid
    is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke

    do k = 1, nz
      do i=is-1,ie+1 ; do j=js-1,je+1

        h_new(i,j,k) = max(0.0, G%areaT(i,j)*h_pre(i,j,k) + &
          ((uhtr(I-1,j,k) - uhtr(I,j,k)) + (vhtr(i,J-1,k) - vhtr(i,J,k))))

        ! In the case that the layer is now dramatically thinner than it was previously,
        ! add a bit of mass to avoid truncation errors.  This will lead to
        ! non-conservation of tracers
        h_new(i,j,k) = h_new(i,j,k) + &
          max(GV%Angstrom, 1.0e-13*h_new(i,j,k) - G%areaT(i,j)*h_pre(i,j,k))

        ! Convert back to thickness
        h_new(i,j,k) = h_new(i,j,k)/G%areaT(i,j)

      enddo ; enddo
    enddo
  end subroutine update_h_horizontal_flux

  !> Updates layer thicknesses due to vertical mass transports
  !! NOTE: Only used in non-ALE configuration
  subroutine update_h_vertical_flux(G, GV, ea, eb, h_pre, h_new)
    type(ocean_grid_type),    pointer                           :: G
    type(verticalGrid_type),  pointer                           :: GV
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(in)       :: ea
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(in)       :: eb
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(in)       :: h_pre
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)    :: h_new

    ! Local variables
    integer :: i, j, k, m, is, ie, js, je, nz
    ! Set index-related variables for fields on T-grid
    is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke

    ! Update h_new with convergence of vertical mass transports
    do j=js-1,je+1
      do i=is-1,ie+1

        ! Top layer
        h_new(i,j,1) = max(0.0, h_pre(i,j,1) + (eb(i,j,1) - ea(i,j,2) + ea(i,j,1) ))
        h_new(i,j,1) = h_new(i,j,1) + &
            max(0.0, 1.0e-13*h_new(i,j,1) - h_pre(i,j,1))

        ! Bottom layer
!        h_new(i,j,nz) = h_pre(i,j,nz) + (ea(i,j,nz) - eb(i,j,nz-1)+eb(i,j,nz))
        h_new(i,j,nz) = max(0.0, h_pre(i,j,nz) + (ea(i,j,nz) - eb(i,j,nz-1)+eb(i,j,nz)))
        h_new(i,j,nz) = h_new(i,j,nz) + &
            max(0.0, 1.0e-13*h_new(i,j,nz) - h_pre(i,j,nz))

      enddo

      ! Interior layers
      do k=2,nz-1 ; do i=is-1,ie+1

        h_new(i,j,k) = max(0.0, h_pre(i,j,k) + ((ea(i,j,k) - eb(i,j,k-1)) + &
            (eb(i,j,k) - ea(i,j,k+1))))
        h_new(i,j,k) = h_new(i,j,k) + &
          max(0.0, 1.0e-13*h_new(i,j,k) - h_pre(i,j,k))

      enddo ; enddo

    enddo

  end subroutine update_h_vertical_flux

  !> This routine limits the mass fluxes so that the a layer cannot be completely depleted.
  !! NOTE: Only used in non-ALE mode
  subroutine limit_mass_flux_3d(G, GV, uh, vh, ea, eb, h_pre, max_off_cfl)
    type(ocean_grid_type),    pointer                           :: G
    type(verticalGrid_type),  pointer                           :: GV
    real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout)    :: uh
    real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout)    :: vh
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)    :: ea
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)    :: eb
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(in)       :: h_pre
    real,                                      intent(in)       :: max_off_cfl

    ! Local variables
    integer :: i, j, k, m, is, ie, js, je, nz
    real, dimension(SZI_(G),SZJ_(G),SZK_(G))                    :: top_flux, bottom_flux
    real                                                        :: pos_flux, hvol, h_neglect, scale_factor


    ! In this subroutine, fluxes out of the box are scaled away if they deplete
    ! the layer, note that we define the positive direction as flux out of the box.
    ! Hence, uh(I-1) is multipled by negative one, but uh(I) is not

    ! Set index-related variables for fields on T-grid
    is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke

    ! Calculate top and bottom fluxes from ea and eb. Note the explicit negative signs
    ! to enforce the positive out convention
    k = 1
    do j=js-1,je+1 ; do i=is-1,ie+1
      top_flux(i,j,k) = -ea(i,j,k)
      bottom_flux(i,j,k) = -(eb(i,j,k)-ea(i,j,k+1))
    enddo ; enddo

    do k=2, nz-1 ; do j=js-1,je+1 ; do i=is-1,ie+1
      top_flux(i,j,k) = -(ea(i,j,k)-eb(i,j,k-1))
      bottom_flux(i,j,k) = -(eb(i,j,k)-ea(i,j,k+1))
    enddo ; enddo ; enddo

    k=nz
    do j=js-1,je+1 ; do i=is-1,ie+1
      top_flux(i,j,k) = -(ea(i,j,k)-eb(i,j,k-1))
      bottom_flux(i,j,k) = -eb(i,j,k)
    enddo ; enddo


    ! Calculate sum of positive fluxes (negatives applied to enforce convention)
    ! in a given cell and scale it back if it would deplete a layer
    do k = 1, nz ; do j=js-1,je+1 ; do i=is-1,ie+1

      hvol = h_pre(i,j,k)*G%areaT(i,j)
      pos_flux  = max(0.0,-uh(I-1,j,k)) + max(0.0, -vh(i,J-1,k)) + &
        max(0.0, uh(I,j,k)) + max(0.0, vh(i,J,k)) + &
        max(0.0, top_flux(i,j,k)*G%areaT(i,j)) + max(0.0, bottom_flux(i,j,k)*G%areaT(i,j))

      if (pos_flux>hvol .and. pos_flux>0.0) then
        scale_factor = ( hvol )/pos_flux*max_off_cfl
      else ! Don't scale
        scale_factor = 1.0
      endif

      ! Scale horizontal fluxes
      if (-uh(I-1,j,k)>0) uh(I-1,j,k) = uh(I-1,j,k)*scale_factor
      if (uh(I,j,k)>0)    uh(I,j,k)   = uh(I,j,k)*scale_factor
      if (-vh(i,J-1,k)>0) vh(i,J-1,k) = vh(i,J-1,k)*scale_factor
      if (vh(i,J,k)>0)    vh(i,J,k)   = vh(i,J,k)*scale_factor

      if (k>1 .and. k<nz) then
      ! Scale interior layers
        if(top_flux(i,j,k)>0.0) then
          ea(i,j,k) = ea(i,j,k)*scale_factor
          eb(i,j,k-1) = eb(i,j,k-1)*scale_factor
        endif
        if(bottom_flux(i,j,k)>0.0) then
          eb(i,j,k) = eb(i,j,k)*scale_factor
          ea(i,j,k+1) = ea(i,j,k+1)*scale_factor
        endif
      ! Scale top layer
      elseif (k==1) then
        if(top_flux(i,j,k)>0.0)    ea(i,j,k) = ea(i,j,k)*scale_factor
        if(bottom_flux(i,j,k)>0.0) then
          eb(i,j,k)   = eb(i,j,k)*scale_factor
          ea(i,j,k+1) = ea(i,j,k+1)*scale_factor
        endif
      ! Scale bottom layer
      elseif (k==nz) then
        if(top_flux(i,j,k)>0.0) then
          ea(i,j,k)   = ea(i,j,k)*scale_factor
          eb(i,j,k-1) = eb(i,j,k-1)*scale_factor
        endif
        if (bottom_flux(i,j,k)>0.0) eb(i,j,k)=eb(i,j,k)*scale_factor
      endif
    enddo ; enddo ; enddo

  end subroutine limit_mass_flux_3d
  
  !> In the case where offline advection has failed to converge. Redistribute the flux
  !! into remainder of the water column
  subroutine redistribute_residual(G, GV, h, uhtr, vhtr)
    type(ocean_grid_type),    pointer                           :: G
    type(verticalGrid_type),  pointer                           :: GV
    real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout)    :: uh
    real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout)    :: vh
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)    :: ea
    
  end subroutine distribute_residual_upwards
  
!> \namespace mom_offline_transport
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
!!        -#  Force a remapping to the stored layer thicknesses that correspond to the snapshot of
!!            the online model at the end of an accumulation interval
!!        -3  Reset T/S and h to their stored snapshotted values to prevent model drift
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
  
end module MOM_offline_transport

