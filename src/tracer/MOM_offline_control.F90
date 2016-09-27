!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the impliec warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Andrew Shao 2016                                                *
!*                                                                     *
!*  The subroutines here allow MOM6 to be run in a so-called 'offline' *
!*  mode ostensibly for the purpose of modeling tracers. Instead of    *
!*  calculating u, v, and h prognostically, these fields are read in   *
!*  at regular intervals which have been saved from a previous         *
!*  integration of MOM6.                                               *
!*                                                                     *
!*  Users are warned that the usual diagnostics (i.e. conservation of  *
!*  mass) cannot be expected to be replicated to the same accuracy as  *
!*  the online model because some information is loss due to the       *
!*  averaging and snapshotting involved with saving offline files. The *
!*  responsibility lies on the user that the loss of accuracy is       *
!*  acceptable for their application.
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h  *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

module MOM_offline_transport


  use data_override_mod,  only : data_override_init, data_override
  use MOM_time_manager,   only : time_type
  use MOM_domains,        only : pass_var, pass_vector, To_All
  use MOM_error_handler,  only : callTree_enter, callTree_leave, MOM_error, FATAL, WARNING, is_root_pe
  use MOM_grid,           only : ocean_grid_type
  use MOM_verticalGrid,   only : verticalGrid_type
  use MOM_io,             only : read_data
  use MOM_file_parser,    only : get_param, log_version, param_file_type
  use MOM_diag_mediator,  only : diag_ctrl, register_diag_field
  use mpp_domains_mod,    only : CENTER, CORNER, NORTH, EAST
  use MOM_variables,      only : vertvisc_type
  use MOM_forcing_type,   only : forcing
  use MOM_shortwave_abs,  only : optics_type
  use MOM_diag_mediator,  only : post_data

  implicit none

#include <MOM_memory.h>

  type, public :: offline_transport_CS

    integer :: start_index  ! Timelevel to start
    integer :: numtime     ! How many timelevels in the input fields

    integer :: &            ! Index of each of the variables to be read in
      ridx_mean = -1, &     ! Separate indices for each variabile if they are
      ridx_snap = -1      ! setoff from each other in time


    character(len=200) :: offlinedir  ! Directory where offline fields are stored
    character(len=200) :: & ! Names
      mean_file,  &
      snap_file,  &
      sum_file,   &
      preale_file

    logical :: fields_are_offset ! True if the time-averaged fields and snapshot fields are
                                 ! offset by one time level

    real :: max_off_cfl
    ! These fields for preale are allocatable because they are not necessary for all runs
    real, allocatable, dimension(NIMEM_,NJMEM_,NKMEM_)      :: &
      T_preale, &
      S_preale, &
      h_preale
    real, allocatable, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
      u_preale
    real, allocatable, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
      v_preale

    real              ::  dt_offline ! Timestep used for offline tracers

    integer           ::  num_off_iter

    integer :: &
      id_uhtr_preadv = -1, &
      id_vhtr_preadv = -1, &
      id_temp_preadv = -1, &
      id_salt_preadv = -1, &
      id_uhr = -1, &
      id_vhr = -1, &
      id_ear = -1, &
      id_ebr = -1, &
      id_hr = -1

  end type offline_transport_CS

#include "MOM_memory.h"
#include "version_variable.h"
  public offline_transport_init
  public post_advection_fields
  public transport_by_files
  public register_diags_offline_transport
  public update_h_horizontal_flux
  public update_h_vertical_flux
  public limit_mass_flux_3d

contains

    ! Called right before tracer_advect call in MOM.F90 to ensure that all terms
    ! in the tracer advection routine are the same online and offline
  subroutine post_advection_fields( G, CS, diag, h_adv, uhtr, vhtr, temp, salt )

    type(ocean_grid_type),                     intent(in)       :: G
    type(offline_transport_CS),                intent(in)       :: CS
    type(diag_ctrl),                           intent(inout)    :: diag
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)       :: h_adv, temp, salt
    real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)       :: uhtr   !< accumulated volume/mass flux through zonal face (m3 or kg)
    real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)       :: vhtr   !< accumulated volume/mass flux through merid face (m3 or kg)

    real, dimension(SZI_(G),SZJ_(G),SZK_(G))                    :: write_all_3dt
    real, dimension(SZIB_(G),SZJ_(G),SZK_(G))                   :: write_all_3du
    real, dimension(SZI_(G),SZJB_(G),SZK_(G))                   :: write_all_3dv


    write_all_3dt(:,:,:) = 1.
    write_all_3du(:,:,:) = 1.
    write_all_3dv(:,:,:) = 1.


    if (CS%id_uhtr_preadv>0)   call post_data(CS%id_uhtr_preadv,  uhtr, diag )
    if (CS%id_vhtr_preadv>0)   call post_data(CS%id_vhtr_preadv,  vhtr, diag )
    if (CS%id_temp_preadv>0)   call post_data(CS%id_temp_preadv,  temp, diag )
    if (CS%id_salt_preadv>0)   call post_data(CS%id_salt_preadv,  salt, diag )

!    if (CS%id_uhtr_preadv>0)   call post_data(CS%id_uhtr_preadv,  uhtr, diag, mask = write_all_3du )
!    if (CS%id_vhtr_preadv>0)   call post_data(CS%id_vhtr_preadv,  vhtr, diag, mask = write_all_3dv )
!    if (CS%id_temp_preadv>0)   call post_data(CS%id_temp_preadv,  temp, diag, mask = write_all_3dt )
!    if (CS%id_salt_preadv>0)   call post_data(CS%id_salt_preadv,  salt, diag, mask = write_all_3dt )

  end subroutine post_advection_fields

  subroutine transport_by_files(G, GV, CS, h_end, eatr, ebtr, uhtr, vhtr, khdt_x, khdt_y, &
    temp, salt, do_ale_in)
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
    logical                                                     :: do_ale
    integer :: i, j, k, is, ie, js, je, nz

    do_ale = .false.;
    if (present(do_ale_in) ) do_ale = do_ale_in

    is   = G%isc   ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec ; nz = GV%ke


    call callTree_enter("transport_by_files, MOM_offline_control.F90")


    !! Time-summed fields
    ! U-grid
    call read_data(CS%sum_file, 'uhtr_preadv_sum',     uhtr,domain=G%Domain%mpp_domain, &
      timelevel=CS%ridx_mean,position=EAST)
    call read_data(CS%sum_file, 'khdt_x_sum', khdt_x,domain=G%Domain%mpp_domain, &
      timelevel=CS%ridx_mean,position=EAST)
    ! V-grid
    call read_data(CS%sum_file, 'vhtr_preadv_sum',     vhtr, domain=G%Domain%mpp_domain, &
      timelevel=CS%ridx_mean,position=NORTH)
    call read_data(CS%sum_file, 'khdt_y_sum', khdt_y, domain=G%Domain%mpp_domain, &
      timelevel=CS%ridx_mean,position=NORTH)
    ! T-grid
    call read_data(CS%sum_file, 'ea_sum',   eatr, domain=G%Domain%mpp_domain, &
      timelevel=CS%ridx_mean,position=CENTER)
    call read_data(CS%sum_file, 'eb_sum',   ebtr, domain=G%Domain%mpp_domain, &
      timelevel=CS%ridx_mean,position=CENTER)

    !! Time-averaged fields
    call read_data(CS%snap_file, 'temp_preadv',   temp, domain=G%Domain%mpp_domain, &
      timelevel=CS%ridx_mean,position=CENTER)
    call read_data(CS%snap_file, 'salt_preadv',   salt, domain=G%Domain%mpp_domain, &
      timelevel=CS%ridx_mean,position=CENTER)

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


    if (do_ale) then
      CS%h_preale = GV%Angstrom
      CS%T_preale = 0.0
      CS%S_preale = 0.0
      CS%u_preale = 0.0
      CS%v_preale = 0.0
      call read_data(CS%preale_file, 'h_preale',   CS%h_preale, domain=G%Domain%mpp_domain, &
        timelevel=CS%ridx_snap,position=CENTER)
      call read_data(CS%preale_file, 'T_preale',   CS%T_preale, domain=G%Domain%mpp_domain, &
        timelevel=CS%ridx_mean,position=CENTER)
      call read_data(CS%preale_file, 'S_preale',   CS%S_preale, domain=G%Domain%mpp_domain, &
        timelevel=CS%ridx_mean,position=CENTER)
      call read_data(CS%preale_file, 'u_preale',   CS%u_preale, domain=G%Domain%mpp_domain, &
        timelevel=CS%ridx_mean,position=EAST)
      call read_data(CS%preale_file, 'v_preale',   CS%v_preale, domain=G%Domain%mpp_domain, &
        timelevel=CS%ridx_mean,position=NORTH)

!      do k=1,nz ; do j=js-1,je ; do i=is-1,ie
!        if (G%mask2dCu(i,j)<1.0) then
!          CS%u_preale(I,j,k) = 0.0
!        endif
!        if (G%mask2dCv(i,j)<1.0) then
!          CS%v_preale(I,j,k) = 0.0
!        endif
!        if (G%mask2dT(i,j)<1.0) then
!          CS%h_preale(i,j,k) = GV%Angstrom
!        endif
!      enddo; enddo; enddo

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

      call pass_vector(CS%u_preale,CS%v_preale,G%Domain)
      call pass_var(CS%h_preale, G%Domain)
      call pass_var(CS%T_preale, G%Domain)
      call pass_var(CS%S_preale, G%Domain)


    endif

    ! Update the read indices
    CS%ridx_snap = next_modulo_time(CS%ridx_snap,CS%numtime)
    CS%ridx_mean = next_modulo_time(CS%ridx_mean,CS%numtime)

    call callTree_leave("transport_by_file")

  end subroutine transport_by_files

  !> Initialize additional diagnostics required for offline tracer transport
  subroutine register_diags_offline_transport(Time, diag, CS)

    type(offline_transport_CS), pointer :: CS         !< control structure for MOM
    type(time_type), intent(in) :: Time             !< current model time
    type(diag_ctrl)             :: diag


    ! U-cell fields
    CS%id_uhtr_preadv = register_diag_field('ocean_model', 'uhtr_preadv', diag%axesCuL, Time, &
      'Accumulated zonal thickness fluxes to advect tracers', 'kg')
    CS%id_uhr = register_diag_field('ocean_model', 'uhr', diag%axesCuL, Time, &
      'Zonal thickness fluxes remaining at end of timestep', 'kg')

    ! V-cell fields
    CS%id_vhtr_preadv = register_diag_field('ocean_model', 'vhtr_preadv', diag%axesCvL, Time, &
      'Accumulated meridional thickness fluxes to advect tracers', 'kg')
    CS%id_vhr = register_diag_field('ocean_model', 'vhr', diag%axesCvL, Time, &
      'Meridional thickness fluxes remaining at end of timestep', 'kg')

    ! T-cell fields
    CS%id_temp_preadv  = register_diag_field('ocean_model', 'temp_preadv', diag%axesTL, Time, &
      'Temperature prior to advection', 'C')
    CS%id_salt_preadv  = register_diag_field('ocean_model', 'salt_preadv', diag%axesTL, Time, &
      'Salinity prior to advection', 'S')
    CS%id_hr  = register_diag_field('ocean_model', 'hdiff', diag%axesTL, Time, &
      'Difference between the stored and calculated layer thickness', 'm')
    CS%id_ear  = register_diag_field('ocean_model', 'ear', diag%axesTL, Time, &
      'Remaining thickness entrained from above', 'm')
    CS%id_ebr  = register_diag_field('ocean_model', 'ebr', diag%axesTL, Time, &
      'Remaining thickness entrained from below', 'm')

  end subroutine register_diags_offline_transport

  subroutine offline_transport_init(param_file, CS, do_ale, G, GV)

    type(param_file_type)              , intent(in) :: param_file
    type(offline_transport_CS), pointer, intent(inout) :: CS
    logical                            , intent(in) :: do_ale
    type(ocean_grid_type),               intent(in) :: G
    type(verticalGrid_type),             intent(in) :: GV

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
    call get_param(param_file, mod, "OFF_MEAN_FILE", CS%mean_file, &
      "Filename where time-averaged fields are fund can be found", default=" ")
    call get_param(param_file, mod, "OFF_SUM_FILE", CS%sum_file, &
      "Filename where the accumulated fields can be found", default = " ")
    call get_param(param_file, mod, "OFF_SNAP_FILE", CS%snap_file, &
      "Filename where snapshot fields can be found",default=" ")
    call get_param(param_file, mod, "OFF_PREALE_FILE", CS%preale_file, &
      "Filename where the preale T, S, u, v, and h fields are found",default=" ")
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
    call get_param(param_file, "MOM_mixed_layer", "MAX_OFF_CFL", CS%max_off_cfl, &
      "Maximum CFL when advection is done offline. This should be less than 1 \n", &
      units="nondim", default=0.9)

    ! Concatenate offline directory and file names
    CS%mean_file = trim(CS%offlinedir)//trim(CS%mean_file)
    CS%snap_file = trim(CS%offlinedir)//trim(CS%snap_file)
    CS%sum_file = trim(CS%offlinedir)//trim(CS%sum_file)
    CS%preale_file = trim(CS%offlinedir)//trim(CS%preale_file)

    ! Set the starting read index for time-averaged and snapshotted fields
    CS%ridx_mean = CS%start_index
    if(CS%fields_are_offset) CS%ridx_snap = next_modulo_time(CS%start_index,CS%numtime)
    if(.not. CS%fields_are_offset) CS%ridx_snap = CS%start_index

    if (do_ale) then
      ALLOC_(CS%u_preale(IsdB:IedB,jsd:jed,nz))   ; CS%u_preale(:,:,:) = 0.0
      ALLOC_(CS%v_preale(isd:ied,JsdB:JedB,nz))   ; CS%v_preale(:,:,:) = 0.0
      ALLOC_(CS%h_preale(isd:ied,jsd:jed,nz))     ; CS%h_preale(:,:,:) = GV%Angstrom
      ALLOC_(CS%T_preale(isd:ied,jsd:jed,nz))     ; CS%T_preale(:,:,:) = 0.0
      ALLOC_(CS%S_preale(isd:ied,jsd:jed,nz))     ; CS%S_preale(:,:,:) = 0.0
    endif

    call callTree_leave("offline_transport_init")

  end subroutine offline_transport_init

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

  subroutine limit_mass_flux_ordered_3d(G, GV, uh, vh, ea, eb, h_in, h_end, max_off_cfl, z_first, x_before_y)
    type(ocean_grid_type),    pointer                           :: G
    type(verticalGrid_type),  pointer                           :: GV
    real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout)    :: uh
    real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout)    :: vh
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)    :: ea
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)    :: eb
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(in)       :: h_in
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)       :: h_end
    real,                                      intent(in)       :: max_off_cfl
    logical,                                   intent(in)       :: z_first, x_before_y

    ! Local variables
    integer :: i, j, k, m, is, ie, js, je, nz
    real, dimension(SZI_(G),SZJ_(G),SZK_(G))                    :: h_budget ! Tracks how much thickness
                                                                            ! remains for other fluxes
    integer                                                     :: flux_order = -1
    ! In this subroutine, fluxes out of the box are scaled down if they deplete
    ! the layer. Here the positive direction is defined as flux out of the box as opposed to the
    ! typical, strictly upwind convention. Hence, uh(I-1) is multipled by negative one,
    ! but uh(I) is not. This routine differs from limit_mass_flux_3d because in this case,
    ! the ordering of direction matters. While this is more aggressive than the other routine which,
    ! Scales fluxes if they would deplete the layer (independent of any convergence within an
    ! iteration), this routine should still maintain a CFL less than 1
    ! Because horizontal transport must always be together (i.e. cannot do x->z->y),
    ! four cases are considered)
    !   1: z -> x -> y
    !   2: z -> y -> x
    !   3: x -> y -> z
    !   4: y -> x -> z

    ! Set index-related variables for fields on T-grid
    is  = G%isd ; ie  = G%ied ; js  = G%jsd ; je  = G%jed ; nz = GV%ke
    ! Copy layer thicknesses into a working array for this subroutine
    do k=1,nz ; do j=js,je ; do i=is,ie
      h_budget(i,j,k) = h_in(i,j,k)*G%areaT(i,j)
    enddo ; enddo ; enddo

    ! Set the flux order (corresponding to one of the four cases described previously)
    if (z_first .and. x_before_y)                 flux_order = 1
    if (z_first .and. (.not. x_before_y))         flux_order = 2
    if ((.not. z_first) .and. x_before_y)         flux_order = 3
    if ((.not. z_first) .and. (.not. x_before_y)) flux_order = 4

    select case (flux_order)
      case (1) ! z -> x -> y
        ! Check first to see if either the top or bottom flux would deplete the layer
        !call flux_limiter_vertical(G, GV, h_budget, ea, eb, max_off_cfl)
        call flux_limiter_u(G, GV, h_budget, uh, max_off_cfl)
        call flux_limiter_v(G, GV, h_budget, vh, max_off_cfl)
      case (2) ! z -> y -> x
        !call flux_limiter_vertical(G, GV, h_budget, ea, eb, max_off_cfl)
        call flux_limiter_v(G, GV, h_budget, vh, max_off_cfl)
        call flux_limiter_u(G, GV, h_budget, uh, max_off_cfl)
      case (3) ! x -> y -> z
        call flux_limiter_u(G, GV, h_budget, uh, max_off_cfl)
        call flux_limiter_v(G, GV, h_budget, vh, max_off_cfl)
!        call flux_limiter_vertical(G, GV, h_budget, ea, eb, max_off_cfl)
      case (4) ! y -> x -> z
        call flux_limiter_v(G, GV, h_budget, vh, max_off_cfl)
        call flux_limiter_u(G, GV, h_budget, uh, max_off_cfl)
!        call flux_limiter_vertical(G, GV, h_budget, ea, eb, max_off_cfl)
      case default
        call MOM_error(FATAL, "Invalid choice of flux_order")
    end select
    
    do k=1,nz ; do j=js,je ; do i=is,ie
      h_end(i,j,k) = h_budget(i,j,k)/G%areaT(i,j)
      if (h_end(i,j,k)<0.0 ) then
        print *, "i,j,k,h,",i,j,k,h_end(i,j,k)
      endif
    enddo ; enddo ; enddo

  end subroutine limit_mass_flux_ordered_3d

  subroutine flux_limiter_vertical(G, GV, h, ea, eb, max_off_cfl)
    type(ocean_grid_type),    pointer                           :: G
    type(verticalGrid_type),  pointer                           :: GV
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)    :: ea
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)    :: eb
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)       :: h
    real                                                        :: max_off_cfl
    ! Limits how much the a layer can be depleted in the vertical direction
    real, dimension(SZI_(G),SZK_(G))                            :: ea2d, eb2d
    real, dimension(SZI_(G),SZK_(G))                            :: h2d, scale
    real, dimension(SZI_(G),SZK_(G))                            :: top_flux, bottom_flux
    real                                                        :: total_out_flux, h_budget
    integer :: i, j, k, m, is, ie, js, je, nz

    ! Set index-related variables for fields on T-grid
    is  = G%isd ; ie  = G%ied ; js  = G%jsd ; je  = G%jed ; nz = GV%ke

    do j=js,je
      do k=1,nz ; do i=is,ie
        ea2d(i,k) = ea(i,j,k)
        eb2d(i,k) = eb(i,j,k)
        h2d(i,k) = h(i,j,k)
        scale(i,k) = 1.0
      enddo ; enddo;

      k=1 ! Top layer
      do i=is,ie
        top_flux(i,k) = -ea2d(i,k)
        bottom_flux(i,k) = -(eb2d(i,k)-ea2d(i,k+1))
      enddo
      ! Interior layers
      do k=2, nz-1 ; do i=is,ie
        top_flux(i,k) = -(ea2d(i,k)-eb2d(i,k-1))
        bottom_flux(i,k) = -(eb2d(i,k)-ea2d(i,k+1))
      enddo ; enddo
      k=nz ! Bottom layer
      do i=is,ie
        top_flux(i,k) = -(ea2d(i,k)-eb2d(i,k-1))
        bottom_flux(i,k) = -eb2d(i,k)
      enddo

      do k=1,nz ; do i=is,ie
        h_budget = h2d(i,k)*max_off_cfl ! How much the layer can be depleted in any given step
                                        ! based on the specified max CFL
        total_out_flux = (max(0.0,top_flux(i,k)) + max(0.0, bottom_flux(i,k)))*G%areaT(i,j)
        if (total_out_flux>h_budget) scale(i,k) = h_budget/total_out_flux
        if (scale(j,k)>1.0) call MOM_error(FATAL, "scale(j,k) is larger than 1")
      enddo ; enddo

      k=1
      do i=is,ie
        if(top_flux(i,k)>0.0) then
          ea2d(i,k) = ea2d(i,k)*scale(i,k)
        endif
        if(bottom_flux(i,k)>0.0) then
          ea2d(i,k+1) = ea2d(i,k+1)*scale(i,k)
          eb2d(i,k)  = eb2d(i,k)*scale(i,k)
        endif
      enddo
      ! Interior layers
      do k=2, nz-1 ; do i=is,ie
        if(top_flux(i,k)>0.0) then
          ea2d(i,k)   = ea2d(i,k)*scale(i,k)
          eb2d(i,k-1) = eb2d(i,k-1)*scale(i,k)
        endif
        if(bottom_flux(i,k)>0.0) then
          ea2d(i,k+1) = ea2d(i,k+1)*scale(i,k)
          eb2d(i,k)   = eb2d(i,k)*scale(i,k)
        endif
      enddo; enddo;
      k=nz
      do i=is,ie
        if(top_flux(i,k)>0.0) then
          ea2d(i,k)   = ea2d(i,k)*scale(i,k)
          eb2d(i,k-1) = eb2d(i,k-1)*scale(i,k)
        endif
        if(bottom_flux(i,k)>0.0) then
          eb2d(i,k)   = eb2d(i,k)*scale(i,k)
        endif
      enddo

      ! Update h with new scaled fluxes
      k=1 ! Top layer
      do i=is,ie
        top_flux(i,k) = -ea2d(i,k)
        bottom_flux(i,k) = -(eb2d(i,k)-ea2d(i,k+1))
      enddo
      ! Interior layers
      do k=2, nz-1 ; do i=is,ie
        top_flux(i,k) = -(ea2d(i,k)-eb2d(i,k-1))
        bottom_flux(i,k) = -(eb2d(i,k)-ea2d(i,k+1))
      enddo ; enddo
      k=nz ! Bottom layer
      do i=is,ie
        top_flux(i,k) = -(ea2d(i,k)-eb2d(i,k-1))
        bottom_flux(i,k) = -eb2d(i,k)
      enddo

      do k=1,nz ; do i=is,ie
        h(i,j,k)  = h2d(i,k) - (top_flux(i,k)+bottom_flux(i,k))*G%areaT(i,j)
        ea(i,j,k) = ea2d(i,k)
        eb(i,j,k) = eb2d(i,k)
      enddo; enddo

    enddo
  end subroutine flux_limiter_vertical

  subroutine flux_limiter_u(G, GV, h, uh, max_off_cfl)
    type(ocean_grid_type),    pointer                           :: G
    type(verticalGrid_type),  pointer                           :: GV
    real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout)    :: uh
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)    :: h
    real                                                        :: max_off_cfl
    ! Limits how much the a layer can be depleted in the vertical direction
    real, dimension(SZIB_(G),SZK_(G))                           :: uh2d
    real                                                        :: hup, hlos, min_h, h_remain
    integer :: i, j, k, m, is, ie, js, je, nz

    min_h= 0.1*GV%Angstrom
    ! Set index-related variables for fields on T-grid
    is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke

    do j=js,je
      do k=1,nz ; do i=is-2,ie
        uh2d(I,k) = uh(I,j,k)        
      enddo ; enddo;
      
      do k=1,nz ; do i=is-1,ie
        if(uh2d(I,k)<0.0) then
          hup = h(i+1,j,k) - min_h*G%areaT(i+1,j)
          hlos = max(0.0, uh2d(I+1,k))
          if (( ((hup-hlos)+uh2d(I,k))<0.0) .and. &
              ((0.5*hup + uh2d(I,k))<0.0)) then
                uh2d(I,k) = min(-0.5*hup,-hup+hlos,0.0)
          endif
        else
          hup = h(i,j,k) - G%areaT(i,j)*min_h
          hlos = max(0.0,-uh2d(I-1,k))
          if ((((hup-hlos)-uh2d(I,k))<0.0) .and. &
              ((0.5*hup-uh2d(I,k))<0.0))  then
                uh2d(I,k) = max(0.5*hup,hup-hlos,0.0)
          endif
        endif  
      enddo ; enddo

      do k=1,nz
        do i=is-1,ie 
          uh(I,j,k) = uh2d(I,k)
        enddo
        do i=is,ie 
          h(i,j,k) = h(i,j,k) - (uh2d(I,k) - uh2d(I-1,k))
        enddo
      enddo
      
    enddo  

  end subroutine flux_limiter_u

  subroutine flux_limiter_v(G, GV, h, vh, max_off_cfl)
    type(ocean_grid_type),    pointer                           :: G
    type(verticalGrid_type),  pointer                           :: GV
    real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout)    :: vh
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)    :: h
    real,                                      intent(in)       :: max_off_cfl
    ! Limits how much the a layer can be depleted in the vertical direction
    real, dimension(SZJB_(G),SZK_(G))                           :: vh2d
    real, dimension(SZJ_(G),SZK_(G))                            :: h2d
    real                                                        :: hup, hlos, min_h, h_remain
    integer :: i, j, k, m, is, ie, js, je, nz

    ! Set index-related variables for fields on T-grid
    is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke
    min_h= 0.1*GV%Angstrom
    do i=is,ie
      do k=1,nz ; do j=js-2,je
        vh2d(J,k) = vh(i,J,k)        
      enddo ; enddo;
      do k=1,nz ; do j=js-1,je
        if(vh2d(J,k)<0.0) then
          hup = h(i,j+1,k)-G%areaT(i,j+1)*min_h
          hlos = max(0.0,vh2d(J+1,k))
          if ((((hup-hlos)+vh2d(J,k))<0.0) .and. &
              ((0.5*hup+vh2d(J,k))<0.0)) then
              vh2d(J,k) = min(-0.5*hup,-hup+hlos,0.0)
          endif
        else
          hup = h(i,j,k) -G%areaT(i,j)*min_h
          hlos = max(0.0,-vh2d(J-1,k))
          if ((((hup-hlos)-vh2d(J,k))<0.0) .and. &
              ((0.5*hup - vh2d(J,k))<0.0)) then
              vh2d(J,k) = max(0.5*hup,hup-hlos,0.0)
          endif
        endif
      enddo ; enddo
      
      do k=1,nz
        do j=js-1,je 
          vh(i,J,k) = vh2d(J,k)
        enddo
        do j=js,je 
          h(i,j,k) = h(i,j,k) - (vh2d(J,k) - vh2d(J-1,k))
        enddo
      enddo
    enddo
  end subroutine flux_limiter_v


end module MOM_offline_transport
