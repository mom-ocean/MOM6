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

    use data_override_mod, only : data_override_init, data_override
    use MOM_time_manager,  only : time_type
    use MOM_domains,       only : pass_var, pass_vector, To_All
    use MOM_error_handler, only : callTree_enter, callTree_leave, MOM_error, WARNING, is_root_pe
    use MOM_grid,          only : ocean_grid_type
    use MOM_io,            only : read_data
    use MOM_file_parser,   only : get_param, log_version, param_file_type
    use MOM_diag_mediator, only : diag_ctrl, register_diag_field
    use mpp_domains_mod,      only : CENTER, CORNER, NORTH, EAST
    implicit none


    type, public :: offline_transport_CS

        integer :: total_counter ! How many total timesteps have been taken since
                                 ! the start of the run

        integer :: start_index  ! Timelevel to start
        integer :: numtime     ! How many timelevels in the input fields

        integer :: &            ! Index of each of the variables to be read in
            ridx_mean = -1, &     ! Separate indices for each variabile if they are
            ridx_snap = -1      ! setoff from each other in time


        character(len=200) :: offlinedir  ! Directory where offline fields are stored
        character(len=200) :: & ! Names
            transport_file,   &
            h_file

        logical :: fields_are_offset ! True if the time-averaged fields and snapshot fields are
                                     ! offset by one time level

        integer :: &
            id_h = -1, &
            id_u = -1, id_v = -1, &
            id_uh = -1, id_vh = -1, &
            id_uhtr = -1, id_vhtr = -1, &
            id_eta = -1

    end type offline_transport_CS

#include "MOM_memory.h"

    public next_modulo_time
    public offline_transport_init
    public transport_by_files
    public transport_by_data_override

contains

    function next_modulo_time(inidx, total_counter)
        ! Returns the next time interval to be read
        integer                 :: total_counter        ! How many times advect_tracer has been called
        integer                 :: inidx                ! Number of time levels in the input files

        integer                 :: read_index           ! The index in the input files that corresponds
                                                        ! to the current timestep

        integer                 :: next_modulo_time

        read_index = mod(inidx+1,total_counter)
        if (read_index < 0)  read_index = inidx-read_index
        if (read_index == 0) read_index = 1

        next_modulo_time = read_index

    end function next_modulo_time

    !> Initialize additional diagnostics required for offline tracer transport
    subroutine register_diags_offline_transport(Time, diag, CS)

        type(offline_transport_CS), pointer :: CS         !< control structure for MOM
        type(time_type), intent(in) :: Time             !< current model time
        type(diag_ctrl)             :: diag

        CS%id_uhtr = register_diag_field('ocean_model', 'uhtr_off', diag%axesCuL, Time, &
            'Accumulated zonal thickness fluxes to advect tracers', 'kg')
        CS%id_vhtr = register_diag_field('ocean_model', 'vhtr_off', diag%axesCvL, Time, &
            'Accumulated meridional thickness fluxes to advect tracers', 'kg')
        CS%id_uh = register_diag_field('ocean_model', 'uh_off', diag%axesCuL, Time, &
            'Accumulated meridional thickness fluxes to advect tracers', 'kg')
        CS%id_vh = register_diag_field('ocean_model', 'vh_off', diag%axesCvL, Time, &
            'Accumulated meridional thickness fluxes to advect tracers', 'kg')
        CS%id_u = register_diag_field('ocean_model', 'u_off', diag%axesCuL, Time, &
            'Accumulated meridional thickness fluxes to advect tracers', 'kg')
        CS%id_v = register_diag_field('ocean_model', 'v_off', diag%axesCvL, Time, &
            'Accumulated meridional thickness fluxes to advect tracers', 'kg')
        CS%id_h = register_diag_field('ocean_model', 'h_off', diag%axesTL, Time, &
            'Accumulated meridional thickness fluxes to advect tracers', 'kg')
        CS%id_eta = register_diag_field('ocean_model', 'eta_av', diag%axesT1, Time, &
            'Accumulated meridional thickness fluxes to advect tracers', 'kg')
        CS%id_eta = register_diag_field('ocean_model', 'p_surf_begin', diag%axesT1, Time, &
            'Accumulated meridional thickness fluxes to advect tracers', 'kg')

    end subroutine register_diags_offline_transport

    subroutine offline_transport_init(param_file, CS)

        type(param_file_type)               :: param_file
        type(offline_transport_CS), pointer :: CS

        character(len=40)                   :: mod = "offline_transport"


        call callTree_enter("offline_transport_init, MOM_offline_control.F90")

        if (associated(CS)) then
            call MOM_error(WARNING, "offline_transport_init called with an associated "// &
                "control structure.")
            return
        endif
        allocate(CS)

        CS%total_counter = 0;
        call get_param(param_file, mod, "OFFLINEDIR", CS%offlinedir, &
            "Input directory where the offline fields can be found", fail_if_missing=.true.)
        call get_param(param_file, mod, "TRANSPORT_FILE", CS%transport_file, &
            "Filename where uhtr, vhtr, u, v fields can be found", default="offline_transport.nc")
        call get_param(param_file, mod, "H_FILE", CS%h_file, &
            "Filename where the h field can be found", default="offline_h.nc")
        call get_param(param_file, mod, "START_INDEX", CS%start_index, &
            "Which time index to start from", fail_if_missing=.true.)
        call get_param(param_file, mod, "NUMTIME", CS%numtime, &
            "Number of timelevels in offline input files", fail_if_missing=.true.)
        call get_param(param_file, mod, "FIELDS_ARE_OFFSET", CS%fields_are_offset, &
            "True if the time-averaged fields and snapshot fields are offset by one time level", &
            default=.false.)

        CS%transport_file = trim(CS%offlinedir)//trim(CS%transport_file)
        CS%h_file = trim(CS%offlinedir)//trim(CS%h_file)


        ! Set the starting read index for time-averaged and snapshotted fields
        CS%ridx_mean = CS%start_index
        if(CS%fields_are_offset) CS%ridx_snap = next_modulo_time(CS%start_index,CS%numtime)
        if(.not. CS%fields_are_offset) CS%ridx_snap = CS%start_index

        call callTree_leave("offline_transport_init")

    end subroutine offline_transport_init

    subroutine transport_by_files(G, CS, angstrom, u, v, uh, vh, uhtr, vhtr, h , eta_av, missing)
        type(ocean_grid_type)                    , intent(inout)       :: G
        type(offline_transport_CS)               , intent(inout)       :: CS
        real                                                        :: angstrom
        real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout)    :: u
        real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout)    :: v
        real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout)    :: uh
        real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout)    :: vh
        real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout)    :: uhtr  !< accumulated volume/mass flux through zonal face (m3 or kg)
        real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout)    :: vhtr  !< accumulated volume/mass flux through merid face (m3 or kg)
        real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)    :: h     !< layer thickness after advection (m or kg m-2)
        real, dimension(SZI_(G),SZJ_(G)) , intent(inout)            :: eta_av
        real                                                        :: missing

        integer :: i, j, k

        call callTree_enter("transport_by_files, MOM_offline_control.F90")

        if ( is_root_pe() ) print *, "Read index: ", CS%ridx_mean

        ! Read time-averaged fields (middle of time interval timestamp)
        call read_data(CS%transport_file, 'u', u(:,:,:),domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_mean,position=EAST)
        call read_data(CS%transport_file, 'v', v(:,:,:),domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_mean,position=NORTH)
!
        call read_data(CS%transport_file, 'uh', uh(:,:,:),domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_mean,position=EAST)
        call read_data(CS%transport_file, 'vh', vh(:,:,:),domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_mean,position=NORTH)


        call read_data(CS%transport_file, 'uhtr', uhtr(:,:,:),domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_mean,position=EAST)
        call read_data(CS%transport_file, 'vhtr', vhtr(:,:,:),domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_mean,position=NORTH)
!
!        call read_data(CS%transport_file, 'eta_av', eta_av(:,:), domain=G%Domain%mpp_domain, &
!            timelevel=CS%ridx_mean,position=CENTER)
!
!        ! Read snapshot fields (end of time interval timestamp)
        call read_data(CS%transport_file, 'h', h(:,:,:),domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_mean,position=CENTER)

        ! Apply masks sinice read_data doesn't account for missing values?!
        do k = 1,G%ke
!
!            ! Fields on T-cell
            do j=G%jsd, G%jed ; do i=G%isd, G%ied
                if( h(i,j,k)<0.0 )  h(i,j,k) = angstrom
!                if(eta_av(i,j).EQ.missing) eta_av = 0
            enddo ; enddo
!
!
            ! Fields on U-Grid
            do j=G%jsd, G%jed ; do i=G%isdb, G%iedb
                if(uhtr(i,j,k) .EQ. missing)  uhtr(i,j,k) = 0
                if(uh(i,j,k)   .EQ. missing)  uh(i,j,k) = 0
                if(u(i,j,k)    .EQ. missing)  u(i,j,k) = 0
            enddo ; enddo

            ! Fields on V-Grid
            do j=G%jsdb, G%jedb ; do i=G%isd, G%ied
                if(vhtr(i,j,k) .EQ. missing)  vhtr(i,j,k) = 0
                if(vh(i,j,k)   .EQ. missing)  vh(i,j,k) = 0
                if(v(i,j,k)    .EQ. missing)  v(i,j,k) = 0
            enddo ; enddo

!
!
        enddo
!
!        ! Make sure all halos have been updated
        call pass_vector(uhtr, vhtr, G%Domain)
        call pass_vector(uh, vh, G%Domain)
        call pass_vector(u, v, G%Domain)
        call pass_var(h,G%Domain)
        call pass_var(eta_av, G%Domain)

        ! Update the read indices
        CS%ridx_snap = next_modulo_time(CS%ridx_snap,CS%numtime)
        CS%ridx_mean = next_modulo_time(CS%ridx_mean,CS%numtime)



        call callTree_leave("transport_by_file")

    end subroutine transport_by_files

    subroutine transport_by_data_override(G, day, u, v, uhtr, vhtr, h)
        type(time_type)                          , intent(in)       :: day   !< Current model time
        type(ocean_grid_type)                    , intent(inout)       :: G
        real, dimension(SZI_(G),SZJ_(G),SZK_(G)) , intent(inout)    :: h     !< layer thickness after advection (m or kg m-2)
        real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout)    :: u
        real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout)    :: v
        real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout)    :: uhtr  !< accumulated volume/mass flux through zonal face (m3 or kg)
        real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout)    :: vhtr  !< accumulated volume/mass flux through merid face (m3 or kg)

        ! This subroutine sets the surface wind stresses

        ! Arguments:
        !            state  = structure describing ocean surface state
        !  (out)     fluxes = structure with pointers to forcing fields; unused have NULL ptrs
        !  (in)      day    = time of the fluxes
        !  (in)      G      = ocean grid structure
        !  (in)      CS     = pointer to control struct returned by previous surface_forcing_init call

        integer :: i, j, is_in, ie_in, js_in, je_in

        call callTree_enter("ocean_transport_by_data_override, MOM_offline_control.F90")

        is_in = G%isc - G%isd + 1
        ie_in = G%iec - G%isd + 1
        js_in = G%jsc - G%jsd + 1
        je_in = G%jec - G%jsd + 1

        call data_override('OCN', 'uhtr', uhtr, day, is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)
        call data_override('OCN', 'vhtr', vhtr, day, is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)
        call data_override('OCN', 'u', u, day, is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)
        call data_override('OCN', 'v', v, day, is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)
        call data_override('OCN', 'h', h, day, is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in)

        !        call data_override('OCN', 'uhtr', uhtr, day)
        !        call data_override('OCN', 'vhtr', vhtr, day)
        !        call data_override('OCN', 'u', u, day)
        !        call data_override('OCN', 'v', v, day)
        !        call data_override('OCN', 'h', h, day)

        call pass_vector(uhtr, vhtr, G%Domain)
        call pass_vector(u, v, G%Domain)
        call pass_var(h, G%Domain)
        call callTree_leave("transport_by_data_override")

    end subroutine transport_by_data_override

end module MOM_offline_transport
