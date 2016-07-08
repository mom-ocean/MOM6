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
    use MOM_error_handler,  only : callTree_enter, callTree_leave, MOM_error, WARNING, is_root_pe
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
            transport_file,   &
            h_file, &
            ts_file, &
            preale_file

        logical :: fields_are_offset ! True if the time-averaged fields and snapshot fields are
                                     ! offset by one time level

        ! These fields for preale are allocatable because they are not necessary for all runs
        real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_)      :: &
            T_preale, &
            S_preale, &
            h_preale
        real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
            u_preale
        real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
            v_preale

        integer :: &
            id_h_new = -1, &
            id_h_old = -1, &
            id_h_adv = -1, &
            id_uhtr = -1, &
            id_vhtr = -1, &
            id_eatr = -1, &
            id_ebtr = -1, &
            id_temp = -1, &
            id_salt = -1

    end type offline_transport_CS

#include "MOM_memory.h"
    public offline_transport_init
    public post_diabatic_fields
    public post_advection_fields

contains

    ! Called from call_tracer_column_fns to make sure that all the terms in the
    ! diabatic driver routine are the same online and offline
    subroutine post_diabatic_fields( G, CS, diag, dt, h_old, h_new, eatr, ebtr )

        type(ocean_grid_type),                      intent(in)    :: G
        type(offline_transport_CS),                 intent(in)       :: CS
        type(diag_ctrl),        target,             intent(inout)    :: diag
        real,                                       intent(in)       :: dt
        real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)       :: h_old, h_new, &
            eatr, ebtr

        real                                                        :: Idt
        real, dimension(SZI_(G),SZJ_(G),SZK_(G))                    :: write_all_3dt
        write_all_3dt = 1.

        if (CS%id_h_old>0)  call post_data(CS%id_h_old, h_old,  diag, mask = write_all_3dt )
        if (CS%id_h_new>0)  call post_data(CS%id_h_new, h_new,  diag, mask = write_all_3dt )
        if (CS%id_eatr>0)   call post_data(CS%id_eatr,  eatr,   diag, mask = write_all_3dt)
        if (CS%id_ebtr>0)   call post_data(CS%id_ebtr,  ebtr,   diag, mask = write_all_3dt)

    end subroutine post_diabatic_fields

    ! Called right before tracer_advect call in MOM.F90 to ensure that all terms
    ! in the tracer advection routine are the same online and offline
    subroutine post_advection_fields( G, CS, diag, dt, h_adv, uhtr, vhtr, temp, salt )

        type(ocean_grid_type),                     intent(in)       :: G
        type(offline_transport_CS),                intent(in)       :: CS
        type(diag_ctrl),        target,            intent(inout)    :: diag
        real,                                      intent(in)       :: dt
        real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)       :: h_adv, temp, salt
        real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)       :: uhtr   !< accumulated volume/mass flux through zonal face (m3 or kg)
        real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)       :: vhtr   !< accumulated volume/mass flux through merid face (m3 or kg)

        real, dimension(SZI_(G),SZJ_(G),SZK_(G))                    :: write_all_3dt
        real, dimension(SZIB_(G),SZJ_(G),SZK_(G))                   :: write_all_3du
        real, dimension(SZI_(G),SZJB_(G),SZK_(G))                   :: write_all_3dv


        write_all_3dt = 1.
        write_all_3du = 1.
        write_all_3dv = 1.

        if (CS%id_h_adv>0)  call post_data(CS%id_h_adv, h_adv,  diag, mask = write_all_3dt )
        if (CS%id_uhtr>0)   call post_data(CS%id_uhtr,  uhtr,   diag, mask = write_all_3du )
        if (CS%id_vhtr>0)   call post_data(CS%id_vhtr,  vhtr,   diag, mask = write_all_3dv )
        if (CS%id_temp>0)   call post_data(CS%id_temp,  temp,   diag, mask = write_all_3dt )
        if (CS%id_salt>0)   call post_data(CS%id_salt,  salt,   diag, mask = write_all_3dt )

    end subroutine post_advection_fields

    subroutine transport_by_files(G, CS, h_old, h_new, h_adv, h_end, eatr, ebtr, uhtr, vhtr, khdt_x, khdt_y, &
        temp, salt, dt, fluxes, optics, do_ale_in)
        type(ocean_grid_type),                     intent(inout)    :: G
        type(offline_transport_CS),                intent(inout)    :: CS
        type(forcing),                             intent(inout)    :: fluxes
        type(optics_type),                         intent(inout)    :: optics
        real,                                      intent(in)       :: dt
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
            h_new, h_old, h_adv, h_end, &
            eatr, ebtr, &
            temp, salt
        logical                                                     :: do_ale

        do_ale = .false.;
        if (present(do_ale_in) ) do_ale = do_ale_in


        call callTree_enter("transport_by_files, MOM_offline_control.F90")


        !! Time-averaged fields
        ! U-cell fields
        call read_data(CS%transport_file, 'uhtr_sum',     uhtr,domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_mean,position=EAST)
        call read_data(CS%transport_file, 'khdt_x_sum', khdt_x,domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_mean,position=EAST)

        ! V-cell fields
        call read_data(CS%transport_file, 'vhtr_sum',     vhtr, domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_mean,position=NORTH)
        call read_data(CS%transport_file, 'khdt_y_sum', khdt_y, domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_mean,position=NORTH)

        ! T-cell fields
        call read_data(CS%transport_file, 'eatr_sum',   eatr, domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_mean,position=CENTER)
        call read_data(CS%transport_file, 'ebtr_sum',   ebtr, domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_mean,position=CENTER)
        call read_data(CS%ts_file, 'temp',   temp, domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_mean,position=CENTER)
        call read_data(CS%ts_file, 'salt',   salt, domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_mean,position=CENTER)

        !! Read snapshot fields (end of time interval timestamp)
        ! T-cell fields
        call read_data(CS%h_file, 'h_new', h_new, domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_snap,position=CENTER)
        call read_data(CS%h_file, 'h_old', h_old, domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_snap,position=CENTER)
        call read_data(CS%h_file, 'h_adv', h_adv, domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_snap,position=CENTER)
        call read_data(CS%h_file, 'h_end', h_end, domain=G%Domain%mpp_domain, &
            timelevel=CS%ridx_snap,position=CENTER)

        if (do_ale) then
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

        endif



        ! Convert all transport from time-averages to total amounts
!        uhtr = uhtr * dt
!        vhtr = vhtr * dt
!        eatr = eatr * dt
!        ebtr = ebtr * dt
!        khdt_x = khdt_x * dt
!        khdt_y = khdt_y * dt

        !! Make sure all halos have been updated
        ! Vector fields
        call pass_vector(uhtr, vhtr, G%Domain)
        call pass_vector(khdt_x, khdt_y, G%Domain)

        ! Scalar fields
        call pass_var(h_adv, G%Domain)
        call pass_var(h_old, G%Domain)
        call pass_var(h_new, G%Domain)
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
        CS%id_uhtr = register_diag_field('ocean_model', 'uh_off', diag%axesCuL, Time, &
            'Accumulated zonal thickness fluxes to advect tracers', 'kg')

        ! V-cell fields
        CS%id_vhtr = register_diag_field('ocean_model', 'vh_off', diag%axesCvL, Time, &
            'Accumulated meridional thickness fluxes to advect tracers', 'kg')

        ! T-cell fields
        CS%id_h_adv = register_diag_field('ocean_model', 'h_adv',    diag%axesTL, Time, &
            'Layer Thickness prior to advection', 'm')
        CS%id_h_old = register_diag_field('ocean_model', 'h_old',    diag%axesTL, Time, &
            'Layer Thickness before diabatic', 'm')
        CS%id_h_new = register_diag_field('ocean_model', 'h_new',    diag%axesTL, Time, &
            'Layer Thickness after diabatic', 'm')
        CS%id_eatr  = register_diag_field('ocean_model', 'eatr_off', diag%axesTL, Time, &
            'Entrainment from layer above', 'kg')
        CS%id_ebtr  = register_diag_field('ocean_model', 'ebtr_off', diag%axesTL, Time, &
            'Entrainment from layer below', 'kg')
        CS%id_temp  = register_diag_field('ocean_model', 'temp_off', diag%axesTL, Time, &
            'Temperature prior to advection', 'C')
        CS%id_salt  = register_diag_field('ocean_model', 'salt_off', diag%axesTL, Time, &
            'Salinity prior to advection', 'S')



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

        ! Parse MOM_input for offline control
        call get_param(param_file, mod, "OFFLINEDIR", CS%offlinedir, &
            "Input directory where the offline fields can be found", fail_if_missing=.true.)
        call get_param(param_file, mod, "TRANSPORT_FILE", CS%transport_file, &
            "Filename where uhtr, vhtr, ea, eb fields can be found")
        call get_param(param_file, mod, "H_FILE", CS%h_file, &
            "Filename where the h fields can be found")
        call get_param(param_file, mod, "TS_FILE", CS%ts_file, &
            "Filename where the temperature and salinity fields can be found")
        call get_param(param_file, mod, "PREALE_FILE", CS%preale_file, &
            "Filename where the preale T, S, u, v, and h fields are found")
        call get_param(param_file, mod, "START_INDEX", CS%start_index, &
            "Which time index to start from", fail_if_missing=.true.)
        call get_param(param_file, mod, "NUMTIME", CS%numtime, &
            "Number of timelevels in offline input files", fail_if_missing=.true.)
        call get_param(param_file, mod, "FIELDS_ARE_OFFSET", CS%fields_are_offset, &
            "True if the time-averaged fields and snapshot fields are offset by one time level", &
            default=.false.)

        ! Concatenate offline directory and file names
        CS%transport_file = trim(CS%offlinedir)//trim(CS%transport_file)
        CS%h_file = trim(CS%offlinedir)//trim(CS%h_file)
        CS%ts_file = trim(CS%offlinedir)//trim(CS%ts_file)
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


end module MOM_offline_transport
