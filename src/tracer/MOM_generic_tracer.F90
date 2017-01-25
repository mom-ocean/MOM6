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
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************
!----------------------------------------------------------------
! <CONTACT EMAIL="Niki.Zadeh@noaa.gov"> Niki Zadeh
! </CONTACT>
!
! <REVIEWER EMAIL="William.Cooke@noaa.gov"> William Cooke
! </REVIEWER>
!
! <OVERVIEW>
!  This module drives the generic version of tracers TOPAZ and CFC
! </OVERVIEW>
!----------------------------------------------------------------

#include <MOM_memory.h>

module MOM_generic_tracer

#ifdef _USE_GENERIC_TRACER
#include <fms_platform.h>

  use mpp_mod,        only: stdout, mpp_error, FATAL,WARNING,NOTE
  use field_manager_mod, only: fm_get_index,fm_string_len

  use generic_tracer, only: generic_tracer_register, generic_tracer_get_diag_list
  use generic_tracer, only: generic_tracer_init, generic_tracer_source, generic_tracer_register_diag
  use generic_tracer, only: generic_tracer_coupler_get, generic_tracer_coupler_set
  use generic_tracer, only: generic_tracer_end, generic_tracer_get_list, do_generic_tracer
  use generic_tracer, only: generic_tracer_update_from_bottom,generic_tracer_vertdiff_G
  use generic_tracer, only: generic_tracer_coupler_accumulate

  use g_tracer_utils,   only: g_tracer_get_name,g_tracer_set_values,g_tracer_set_common,g_tracer_get_common
  use g_tracer_utils,   only: g_tracer_get_next,g_tracer_type,g_tracer_is_prog,g_tracer_flux_init
  use g_tracer_utils,   only: g_tracer_send_diag,g_tracer_get_values
  use g_tracer_utils,   only: g_tracer_get_pointer,g_tracer_get_alias,g_diag_type,g_tracer_set_csdiag

  use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
  use MOM_diag_mediator, only : diag_ctrl, get_diag_time_end
  use MOM_diag_to_Z, only : register_Z_tracer, diag_to_Z_CS
  use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE, is_root_pe
  use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
  use MOM_forcing_type, only : forcing, optics_type
  use MOM_grid, only : ocean_grid_type
  use MOM_hor_index, only : hor_index_type
  use MOM_io, only : file_exists, read_data, slasher, vardesc, var_desc
  use MOM_restart, only : register_restart_field, query_initialized, MOM_restart_CS
  use MOM_spatial_means, only : global_area_mean
  use MOM_sponge, only : set_up_sponge_field, sponge_CS
  use MOM_ALE_sponge, only : set_up_ALE_sponge_field, ALE_sponge_CS
  use MOM_time_manager, only : time_type, get_time, set_time
  use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
  use MOM_tracer_registry, only : register_tracer, tracer_registry_type
  use MOM_tracer_registry, only : add_tracer_diagnostics, add_tracer_OBC_values
  use MOM_tracer_Z_init, only : tracer_Z_init
  use MOM_tracer_initialization_from_Z, only : MOM_initialize_tracer_from_Z
  use MOM_variables, only : surface, thermo_var_ptrs
  use MOM_open_boundary, only : ocean_OBC_type
  use MOM_verticalGrid, only : verticalGrid_type


  implicit none ; private
  logical :: g_registered = .false.

  public register_MOM_generic_tracer, initialize_MOM_generic_tracer
  public MOM_generic_tracer_column_physics, MOM_generic_tracer_surface_state
  public end_MOM_generic_tracer, MOM_generic_tracer_get
  public MOM_generic_tracer_stock
  public MOM_generic_flux_init
  public MOM_generic_tracer_min_max
  public MOM_generic_tracer_fluxes_accumulate

  type, public :: MOM_generic_tracer_CS ; private
     character(len = 200) :: IC_file ! The file in which the generic tracer initial values can
                       ! be found, or an empty string for internal initialization.
     logical :: Z_IC_file ! If true, the generic_tracer IC_file is in Z-space.  The default is false.
     real :: tracer_IC_val = 0.0    ! The initial value assigned to tracers.
     real :: tracer_land_val = -1.0 ! The values of tracers used where  land is masked out.
     logical :: tracers_may_reinit  ! If true, tracers may go through the
                              ! initialization code if they are not found in the
                              ! restart files.

     type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
     type(MOM_restart_CS), pointer :: restart_CSp => NULL()

     !   The following pointer will be directed to the first element of the
     ! linked list of generic tracers.
     type(g_tracer_type), pointer :: g_tracer_list => NULL()
     !   The following pointer will be directed to the first element of the
     ! linked list of generic diagnostics fields that must be Z registered by MOM.
     type(g_diag_type), pointer :: g_diag_list => NULL()

  end type MOM_generic_tracer_CS

! This include declares and sets the variable "version".
#include "version_variable.h"

contains

  ! <SUBROUTINE NAME="register_MOM_generic_tracer">
  !  <OVERVIEW>
  !   Initialize phase I: Add the generic tracers
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine:
  !     Initializes the generic tracer packages and adds their tracers to the list
  !     Adds the tracers in the list of generic tracers to the set of MOM tracers (i.e., MOM-register them)
  !     Register these tracers for restart
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call register_MOM_generic_tracer(G, param_file, CS, diag, tr_Reg, restart_CS)
  !  </TEMPLATE>
  ! </SUBROUTINE>

  function register_MOM_generic_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
    type(hor_index_type),       intent(in)   :: HI
    type(verticalGrid_type),    intent(in) :: GV
    type(param_file_type), intent(in)   :: param_file
    type(MOM_generic_tracer_CS),   pointer      :: CS
    type(tracer_registry_type), pointer     :: tr_Reg
    type(MOM_restart_CS),   pointer     :: restart_CS
    ! This subroutine is used to register tracer fields and subroutines
    ! to be used with MOM.
    ! Arguments: G - The ocean's grid structure.
    !  (in)      param_file - A structure indicating the open file to parse for
    !                         model parameter values.
    !  (in/out)  CS - A pointer that is set to point to the control structure
    !                 for this module
    !  (in)      diag - A structure that is used to regulate diagnostic output.
    !  (in/out)  tr_Reg - A pointer that is set to point to the control structure
    !                  for the tracer advection and diffusion module.
    !  (in)      restart_CS - A pointer to the restart control structure.
    logical :: register_MOM_generic_tracer


    character(len=fm_string_len), parameter :: sub_name = 'register_MOM_generic_tracer'
    character(len=200) :: inputdir ! The directory where NetCDF input files are.
    ! These can be overridden later in via the field manager?

    integer :: ntau, k,i,j,axes(3)
    type(g_tracer_type), pointer      :: g_tracer,g_tracer_next
    character(len=fm_string_len)      :: g_tracer_name,longname,units
    real, dimension(:,:,:,:), pointer   :: tr_field
    real, dimension(:,:,:), pointer     :: tr_ptr
    real, dimension(HI%isd:HI%ied, HI%jsd:HI%jed,GV%ke)         :: grid_tmask
    integer, dimension(HI%isd:HI%ied, HI%jsd:HI%jed)           :: grid_kmt
    type(vardesc) :: vdesc

    register_MOM_generic_tracer = .false.
    if (associated(CS)) then
       call mpp_error(WARNING, "register_MOM_generic_tracer called with an "// &
            "associated control structure.")
       return
    endif
    allocate(CS)


    !Register all the generic tracers used and create the list of them.
    !This can be called by ALL PE's. No array fields allocated.
    if (.not. g_registered) then
       call generic_tracer_register
       g_registered = .true.
    endif


  ! Read all relevant parameters and write them to the model log.
    call log_version(param_file, sub_name, version, "")
    call get_param(param_file, sub_name, "GENERIC_TRACER_IC_FILE", CS%IC_file, &
                 "The file in which the generic trcer initial values can \n"//&
                 "be found, or an empty string for internal initialization.", &
                 default=" ")
    if ((len_trim(CS%IC_file) > 0) .and. (scan(CS%IC_file,'/') == 0)) then
      ! Add the directory if CS%IC_file is not already a complete path.
      call get_param(param_file, sub_name, "INPUTDIR", inputdir, default=".")
      CS%IC_file = trim(slasher(inputdir))//trim(CS%IC_file)
      call log_param(param_file, sub_name, "INPUTDIR/GENERIC_TRACER_IC_FILE", CS%IC_file)
    endif
    call get_param(param_file, sub_name, "GENERIC_TRACER_IC_FILE_IS_Z", CS%Z_IC_file, &
                 "If true, GENERIC_TRACER_IC_FILE is in depth space, not \n"//&
                 "layer space.",default=.false.)
    call get_param(param_file, sub_name, "TRACERS_MAY_REINIT", CS%tracers_may_reinit, &
                 "If true, tracers may go through the initialization code \n"//&
                 "if they are not found in the restart files.  Otherwise \n"//&
                 "it is a fatal error if tracers are not found in the \n"//&
                 "restart files of a restarted run.", default=.false.)

    CS%restart_CSp => restart_CS


    ntau=1 ! MOM needs the fields at only one time step


    !   At this point G%mask2dT and CS%diag%axesTL are not allocated.
    ! postpone diag_registeration to initialize_MOM_generic_tracer

    !Fields cannot be diag registered as they are allocated and have to registered later.
    grid_tmask(:,:,:) = 0.0
    grid_kmt(:,:) = 0.0
    axes(:) = -1

    !
    ! Initialize all generic tracers
    !
    call generic_tracer_init(HI%isc,HI%iec,HI%jsc,HI%jec,HI%isd,HI%ied,HI%jsd,HI%jed,&
         GV%ke,ntau,axes,grid_tmask,grid_kmt,set_time(0,0))


    !
    ! MOM-register the generic tracers
    !

    !Get the tracer list
    call generic_tracer_get_list(CS%g_tracer_list)
    if(.NOT. associated(CS%g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")
    ! For each tracer name get its T_prog index and get its fields

    g_tracer=>CS%g_tracer_list
    do
       call g_tracer_get_alias(g_tracer,g_tracer_name)

       call g_tracer_get_pointer(g_tracer,g_tracer_name,'field',tr_field)
       call g_tracer_get_values(g_tracer,g_tracer_name,'longname', longname)
       call g_tracer_get_values(g_tracer,g_tracer_name,'units',units )

       !nnz: Hard coded stuff. Need get/set routines
       vdesc = var_desc(g_tracer_name, units, longname, &
                        caller="MOM_generic_tracer")
       !!nnz: MOM field is 3D. Does this affect performance? Need it be override field?
       tr_ptr => tr_field(:,:,:,1)
       ! Register tracer for restart file.
       ! mandatory field in restart file is set to .false.
       ! 2008/12/08 jgj: change default to true, so all fields must be present in restart.
       ! 2010/02/04 jgj: if tracers_may_reinit is true, tracers may go through
       ! initialization code if not found in restart
       call register_restart_field(tr_ptr, vdesc, .not.CS%tracers_may_reinit, restart_CS)

       ! Register prognastic tracer for horizontal advection & diffusion. Note
       ! that because the generic tracer code uses only a temporary copy of
       ! the vardesc type, a pointer to this type can not be set as a target
       ! for register_tracer to use.
       if (g_tracer_is_prog(g_tracer)) &
         call register_tracer(tr_ptr, vdesc, param_file, HI, GV, tr_Reg)

       !traverse the linked list till hit NULL
       call g_tracer_get_next(g_tracer, g_tracer_next)
       if(.NOT. associated(g_tracer_next)) exit
       g_tracer=>g_tracer_next

    enddo

    register_MOM_generic_tracer = .true.
  end function register_MOM_generic_tracer

  ! <SUBROUTINE NAME="initialize_MOM_generic_tracer">
  !  <OVERVIEW>
  !   Initialize phase II:  Initialize required variables for generic tracers
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   There are some steps of initialization that cannot be done in register_MOM_generic_tracer
  !   This is the place and time to do them:
  !       Set the grid mask and initial time for all generic tracers.
  !       Diag_register them.
  !       Z_diag_register them.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call initialize_MOM_generic_tracer(restart, day, G, h, OBC, CS, sponge_CSp,ALE_sponge_CSp, diag_to_Z_CSp)
 ! </SUBROUTINE>
  subroutine initialize_MOM_generic_tracer(restart, day, G, GV, h, param_file, diag, OBC, CS, &
                                          sponge_CSp, ALE_sponge_CSp,diag_to_Z_CSp)
    logical,                               intent(in) :: restart
    type(time_type), target,               intent(in) :: day
    type(ocean_grid_type),                 intent(inout) :: G
    type(verticalGrid_type),               intent(in) :: GV
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h
    type(param_file_type),                 intent(in) :: param_file
    type(diag_ctrl),               target, intent(in) :: diag
    type(ocean_OBC_type),                  pointer    :: OBC
    type(MOM_generic_tracer_CS),           pointer    :: CS
    type(sponge_CS),                       pointer    :: sponge_CSp
    type(ALE_sponge_CS),                   pointer    :: ALE_sponge_CSp
    type(diag_to_Z_CS),                    pointer    :: diag_to_Z_CSp
    !   This subroutine initializes the NTR tracer fields in tr(:,:,:,:)
    ! and it sets up the tracer output.

    ! Arguments: restart - .true. if the fields have already been read from
    !                     a restart file.
    !  (in)      day - Time of the start of the run.
    !  (in)      G - The ocean's grid structure.
    !  (in)      GV - The ocean's vertical grid structure.
    !  (in)      h - Layer thickness, in m or kg m-2.
    !  (in)      OBC - This open boundary condition type specifies whether, where,
    !                  and what open boundary conditions are used.
    !  (in/out)  CS - The control structure returned by a previous call to
    !                 register_MOM_generic_tracer.
    !  (in/out)  sponge_CSp - A pointer to the control structure for the sponges, if
    !                         they are in use.  Otherwise this may be unassociated.
    !  (in/out)  ALE_sponge_CSp - A pointer to the control structure for the ALE
    !                 sponges, if they are in use.  Otherwise this may be unassociated
    !  (in/out)  diag_to_Z_Csp - A pointer to the control structure for diagnostics
    !                            in depth space.

    character(len=fm_string_len), parameter :: sub_name = 'initialize_MOM_generic_tracer'
    logical :: OK
    integer :: i, j, k, isc, iec, jsc, jec, nk
    type(g_tracer_type), pointer    :: g_tracer,g_tracer_next
    type(g_diag_type)  , pointer    :: g_diag,g_diag_next
    character(len=fm_string_len)      :: g_tracer_name, longname, units
    real, dimension(:,:,:,:), pointer   :: tr_field
    real, dimension(:,:,:), pointer     :: tr_ptr
    real,    dimension(G%isd:G%ied, G%jsd:G%jed,1:G%ke) :: grid_tmask
    integer, dimension(G%isd:G%ied, G%jsd:G%jed)        :: grid_kmt

    !! 2010/02/04  Add code to re-initialize Generic Tracers if needed during a model simulation
    !! By default, restart cpio should not contain a Generic Tracer IC file and step below will be skipped.
    !! Ideally, the generic tracer IC file should have the tracers on Z levels.

    isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; nk = G%ke

    CS%diag=>diag
    !Get the tracer list
    if(.NOT. associated(CS%g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")
    !For each tracer name get its  fields
    g_tracer=>CS%g_tracer_list

    do
      if(INDEX(CS%IC_file, '_NULL_') .ne. 0) then
         call MOM_error(WARNING,"The name of the IC_file "//trim(CS%IC_file)//&
                              " indicates no MOM initialization was asked for the generic tracers."//&
                              "Bypassing the MOM initialization of ALL generic tracers!")
         exit
      endif
      call g_tracer_get_alias(g_tracer,g_tracer_name)
      call g_tracer_get_pointer(g_tracer,g_tracer_name,'field',tr_field)
      tr_ptr => tr_field(:,:,:,1)

      if (.not.restart .or. (CS%tracers_may_reinit .and. &
          .not.query_initialized(tr_ptr, g_tracer_name, CS%restart_CSp))) then

       if(g_tracer%requires_src_info ) then
         call MOM_error(NOTE,"initialize_MOM_generic_tracer: "//&
                             "initializing generic tracer "//trim(g_tracer_name)//&
                             " using MOM_initialize_tracer_from_Z ")

         call MOM_initialize_tracer_from_Z(h, tr_ptr, G, GV, param_file,                       &
                                src_file = g_tracer%src_file,                              &
                                src_var_nam = g_tracer%src_var_name,                       &
                                src_var_unit_conversion = g_tracer%src_var_unit_conversion,&
                                src_var_record = g_tracer%src_var_record,                  &
                                src_var_gridspec = g_tracer%src_var_gridspec               )

         !Check/apply the bounds for each g_tracer
         do k=1,nk ; do j=jsc,jec ; do i=isc,iec
            if(tr_ptr(i,j,k) .ne. CS%tracer_land_val) then
              if(tr_ptr(i,j,k) .lt. g_tracer%src_var_valid_min) tr_ptr(i,j,k) = g_tracer%src_var_valid_min
              !Jasmin does not want to apply the maximum for now
              !if(tr_ptr(i,j,k) .gt. g_tracer%src_var_valid_max) tr_ptr(i,j,k) = g_tracer%src_var_valid_max
            endif
         enddo; enddo ; enddo

         !jgj: Reset CASED to 0 below K=1
         if(trim(g_tracer_name) .eq. 'cased') then
            do k=2,nk ; do j=jsc,jec ; do i=isc,iec
               if(tr_ptr(i,j,k) .ne. CS%tracer_land_val) then
                 tr_ptr(i,j,k) = 0.0
               endif
            enddo; enddo ; enddo
         endif

       else !Do it old way if the tracer is not registered to start from a specific source file.
            !This path should be deprecated if all generic tracers are required to start from specified sources.
        if (len_trim(CS%IC_file) > 0) then
        !  Read the tracer concentrations from a netcdf file.
          if (.not.file_exists(CS%IC_file)) call MOM_error(FATAL, &
                  "initialize_MOM_Generic_tracer: Unable to open "//CS%IC_file)
          if (CS%Z_IC_file) then
            OK = tracer_Z_init(tr_ptr, h, CS%IC_file, g_tracer_name, G)
            if (.not.OK) then
              OK = tracer_Z_init(tr_ptr, h, CS%IC_file, trim(g_tracer_name), G)
              if (.not.OK) call MOM_error(FATAL,"initialize_MOM_Generic_tracer: "//&
                      "Unable to read "//trim(g_tracer_name)//" from "//&
                      trim(CS%IC_file)//".")
            endif
            call MOM_error(NOTE,"initialize_MOM_generic_tracer: "//&
                            "initialized generic tracer "//trim(g_tracer_name)//&
                            " using Generic Tracer File on Z: "//CS%IC_file)
          else
            ! native grid
            call MOM_error(NOTE,"initialize_MOM_generic_tracer: "//&
                  "Using Generic Tracer IC file on native grid "//trim(CS%IC_file)//&
                  " for tracer "//trim(g_tracer_name))
            call read_data(CS%IC_file, trim(g_tracer_name), tr_ptr, domain=G%Domain%mpp_domain)
          endif
        else
          call MOM_error(FATAL,"initialize_MOM_generic_tracer: "//&
                  "check Generic Tracer IC filename "//trim(CS%IC_file)//".")
        endif

       endif
      endif

      !traverse the linked list till hit NULL
      call g_tracer_get_next(g_tracer, g_tracer_next)
      if(.NOT. associated(g_tracer_next)) exit
      g_tracer=>g_tracer_next
    enddo
    !! end section to re-initialize generic tracers


    !Now we can reset the grid mask, axes and time to their true values
    !Note that grid_tmask must be set correctly on the data domain boundary
    !so that coast mask can be deduced from it.
    grid_tmask(:,:,:) = 0.0
    grid_kmt(:,:) = 0
    do j = G%jsd, G%jed ; do i = G%isd, G%ied
       if (G%mask2dT(i,j) .gt. 0) then
          grid_tmask(i,j,:) = 1.0
          grid_kmt(i,j) = G%ke ! Tell the code that a layer thicker than 1m is the bottom layer.
       endif
    enddo ; enddo
    call g_tracer_set_common(G%isc,G%iec,G%jsc,G%jec,G%isd,G%ied,G%jsd,G%jed,&
                             GV%ke,1,CS%diag%axesTL%handles,grid_tmask,grid_kmt,day)

    ! Register generic tracer modules diagnostics

#ifdef _USE_MOM6_DIAG
    call g_tracer_set_csdiag(CS%diag)
#endif
    call generic_tracer_register_diag()
#ifdef _USE_MOM6_DIAG
    call g_tracer_set_csdiag(CS%diag)
#endif


    ! Register Z diagnostic output.
    !Get the tracer list
    if(.NOT. associated(CS%g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")
    !For each tracer name get its  fields
    g_tracer=>CS%g_tracer_list
    do
       call g_tracer_get_alias(g_tracer,g_tracer_name)

       call g_tracer_get_pointer(g_tracer,g_tracer_name,'field',tr_field)
       tr_ptr => tr_field(:,:,:,1)
       call g_tracer_get_values(g_tracer,g_tracer_name,'longname', longname)
       call g_tracer_get_values(g_tracer,g_tracer_name,'units',units )

       call register_Z_tracer(tr_ptr, trim(g_tracer_name),longname , units, &
            day, G, diag_to_Z_CSp)

       !traverse the linked list till hit NULL
       call g_tracer_get_next(g_tracer, g_tracer_next)
       if(.NOT. associated(g_tracer_next)) exit
       g_tracer=>g_tracer_next

    enddo

    !For each special diagnostics name get its  fields
    !Get the diag list
    call generic_tracer_get_diag_list(CS%g_diag_list)
    if(associated(CS%g_diag_list)) then
       g_diag=>CS%g_diag_list
       do
          if(g_diag%Z_diag .ne. 0) &
               call register_Z_tracer(g_diag%field_ptr, trim(g_diag%name),g_diag%longname , g_diag%units, &
               day, G, diag_to_Z_CSp)

          !traverse the linked list till hit NULL
          g_diag=>g_diag%next
          if(.NOT. associated(g_diag)) exit

       enddo
    endif

  end subroutine initialize_MOM_generic_tracer

  ! <SUBROUTINE NAME="MOM_generic_tracer_column_physics">
  !  <OVERVIEW>
  !   Column physics for generic tracers.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine does:
  !       Get the coupler values for generic tracers that exchange with atmosphere
  !       Update generic tracer concentration fields from sources and sinks.
  !       Vertically diffuse generic tracer concentration fields.
  !       Update generic tracers from bottom and their bottom reservoir.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call MOM_generic_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, CS, tv, optics)
  !  </TEMPLATE>
  ! </SUBROUTINE>

  subroutine MOM_generic_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, CS, tv, optics, &
        evap_CFL_limit, minimum_forcing_depth)
    type(ocean_grid_type),                 intent(in) :: G
    type(verticalGrid_type),               intent(in) :: GV
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h_old, h_new, ea, eb
    type(forcing),                         intent(in) :: fluxes
    real,                                  intent(in) :: dt
    type(MOM_generic_tracer_CS),           pointer    :: CS
    type(thermo_var_ptrs),                 intent(in) :: tv
    type(optics_type),                     intent(in) :: optics
    real,                        optional,intent(in)  :: evap_CFL_limit
    real,                        optional,intent(in)  :: minimum_forcing_depth
    !   This subroutine applies diapycnal diffusion and any other column
    ! tracer physics or chemistry to the tracers from this file.
    ! CFCs are relatively simple, as they are passive tracers. with only a surface
    ! flux as a source.

    ! Arguments: h_old -  Layer thickness before entrainment, in m or kg m-2.
    !  (in)      h_new -  Layer thickness after entrainment, in m or kg m-2.
    !  (in)      ea - an array to which the amount of fluid entrained
    !                 from the layer above during this call will be
    !                 added, in m or kg m-2.
    !  (in)      eb - an array to which the amount of fluid entrained
    !                 from the layer below during this call will be
    !                 added, in m or kg m-2.
    !  (in)      fluxes - A structure containing pointers to any possible
    !                     forcing fields.  Unused fields have NULL ptrs.
    !  (in)      dt - The amount of time covered by this call, in s.
    !  (in)      G - The ocean's grid structure.
    !  (in)      GV - The ocean's vertical grid structure.
    !  (in)      CS - The control structure returned by a previous call to
    !                 register_MOM_generic_tracer.
    !  (in)      evap_CFL_limit - Limits how much water can be fluxed out of the top layer
    !                             Stored previously in diabatic CS.
    !  (in)      minimum_forcing_depth - The smallest depth over which fluxes can be applied
    !                             Stored previously in diabatic CS.
    !
    ! The arguments to this subroutine are redundant in that
    !     h_new[k] = h_old[k] + ea[k] - eb[k-1] + eb[k] - ea[k+1]
    character(len=fm_string_len), parameter :: sub_name = 'MOM_generic_tracer_column_physics'

    type(g_tracer_type), pointer  :: g_tracer, g_tracer_next
    character(len=fm_string_len)  :: g_tracer_name
    real, dimension(:,:), pointer :: stf_array,trunoff_array,runoff_tracer_flux_array

    real :: surface_field(SZI_(G),SZJ_(G))
    real :: sosga

    real, dimension(G%isd:G%ied,G%jsd:G%jed,G%ke) :: rho_dzt, dzt
    real, dimension(G%isd:G%ied,G%jsd:G%jed)      :: hblt_depth
    real, dimension(SZI_(G),SZJ_(G),SZK_(G))      :: h_work
    integer :: i, j, k, isc, iec, jsc, jec, nk

    isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; nk = G%ke

    !Get the tracer list
    if(.NOT. associated(CS%g_tracer_list)) call mpp_error(FATAL,&
         trim(sub_name)//": No tracer in the list.")

#ifdef _USE_MOM6_DIAG
    call g_tracer_set_csdiag(CS%diag)
#endif

    !
    !Extract the tracer surface fields from coupler and update tracer fields from sources
    !
    !call generic_tracer_coupler_get(fluxes%tr_fluxes)
    !Niki: This is moved out to ocean_model_MOM.F90 because if dt_therm>dt_cpld we need to average
    !      the fluxes without coming into this subroutine.
    !      MOM5 has to modified to conform.

    !
    !Add contribution of river to surface flux
    !
    g_tracer=>CS%g_tracer_list
    do
       if(_ALLOCATED(g_tracer%trunoff)) then
          call g_tracer_get_alias(g_tracer,g_tracer_name)
          call g_tracer_get_pointer(g_tracer,g_tracer_name,'stf',   stf_array)
          call g_tracer_get_pointer(g_tracer,g_tracer_name,'trunoff',trunoff_array)
          call g_tracer_get_pointer(g_tracer,g_tracer_name,'runoff_tracer_flux',runoff_tracer_flux_array)
          !nnz: Why is fluxes%river = 0?
          runoff_tracer_flux_array = trunoff_array * fluxes%lrunoff
          stf_array = stf_array + runoff_tracer_flux_array
       endif

       !traverse the linked list till hit NULL
       call g_tracer_get_next(g_tracer, g_tracer_next)
       if(.NOT. associated(g_tracer_next)) exit
       g_tracer=>g_tracer_next

    enddo

    !
    !Prepare input arrays for source update
    !

    rho_dzt(:,:,:) = GV%H_to_kg_m2 * GV%Angstrom
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
      rho_dzt(i,j,k) = GV%H_to_kg_m2 * h_old(i,j,k)
    enddo; enddo ; enddo !}

    dzt(:,:,:) = 1.0
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
      dzt(i,j,k) = GV%H_to_m * h_old(i,j,k)
    enddo; enddo ; enddo !}


    ! Boussinesq model
    hblt_depth(:,:) = GV%H_to_m * GV%Angstrom
    do j=jsc,jec ; do i=isc,iec ;
      hblt_depth(i,j) = GV%H_to_m * h_old(i,j,1)
    enddo; enddo
    do k=2,GV%nkml ; do j=jsc,jec ; do i=isc,iec
      hblt_depth(i,j) = hblt_depth(i,j) + GV%H_to_m * h_old(i,j,k)
    enddo; enddo ; enddo


    do j=jsc,jec ; do i=isc,iec
       surface_field(i,j) = tv%S(i,j,1)
    enddo ; enddo
    sosga = global_area_mean(surface_field, G)

    !
    !Calculate tendencies (i.e., field changes at dt) from the sources / sinks
    !

    call generic_tracer_source(tv%T,tv%S,rho_dzt,dzt,hblt_depth,G%isd,G%jsd,1,dt,&
         G%areaT,get_diag_time_end(CS%diag),&
         optics%nbands, optics%max_wavelength_band, optics%sw_pen_band, optics%opacity_band, sosga=sosga)

    ! This uses applyTracerBoundaryFluxesInOut to handle the change in tracer due to freshwater fluxes
    ! usually in ALE mode
    if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
      g_tracer=>CS%g_tracer_list
      do
        if (g_tracer_is_prog(g_tracer)) then
          do k=1,nk ;do j=jsc,jec ; do i=isc,iec
            h_work(i,j,k) = h_old(i,j,k)
          enddo ; enddo ; enddo;
          call applyTracerBoundaryFluxesInOut(G, GV, g_tracer%field(:,:,:,1), dt, fluxes, h_work, &
              evap_CFL_limit, minimum_forcing_depth)
        endif

         !traverse the linked list till hit NULL
         call g_tracer_get_next(g_tracer, g_tracer_next)
        if(.NOT. associated(g_tracer_next)) exit
        g_tracer=>g_tracer_next
      enddo
    endif

    !
    !Update Tr(n)%field from explicit vertical diffusion
    !
    ! Use a tridiagonal solver to determine the concentrations after the
    ! surface source is applied and diapycnal advection and diffusion occurs.
    if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
      call generic_tracer_vertdiff_G(h_work, ea, eb, dt, GV%kg_m2_to_H, GV%m_to_H, 1) !Last arg is tau which is always 1 for MOM
    else
      call generic_tracer_vertdiff_G(h_old, ea, eb, dt, GV%kg_m2_to_H, GV%m_to_H, 1) !Last arg is tau which is always 1 for MOM
    endif

    ! Update bottom fields after vertical processes

    call generic_tracer_update_from_bottom(dt, 1, get_diag_time_end(CS%diag)) !Second arg is tau which is always 1 for MOM

    !Output diagnostics via diag_manager for all generic tracers and their fluxes
    call g_tracer_send_diag(CS%g_tracer_list, get_diag_time_end(CS%diag), tau=1)
#ifdef _USE_MOM6_DIAG
    call g_tracer_set_csdiag(CS%diag)
#endif


  end subroutine MOM_generic_tracer_column_physics

  ! <SUBROUTINE NAME="MOM_generic_tracer_stock">
  !  <OVERVIEW>
  !   Calculate the mass-weighted integral of tracer concentrations.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !     This subroutine calculates mass-weighted integral on the PE either
  !   of all available tracer concentrations, or of a tracer that is
  !   being requested specifically, returning the number of stocks it has
  !   calculated.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   ns = MOM_generic_tracer_stock(h, stocks, G, CS, names, units, stock_index)
  !  </TEMPLATE>
  ! </SUBROUTINE>

  function MOM_generic_tracer_stock(h, stocks, G, GV, CS, names, units, stock_index)
    type(ocean_grid_type),              intent(in)    :: G
    type(verticalGrid_type),            intent(in)    :: GV
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h
    real, dimension(:),                 intent(out)   :: stocks
    type(MOM_generic_tracer_CS),        pointer       :: CS
    character(len=*), dimension(:),     intent(out)   :: names
    character(len=*), dimension(:),     intent(out)   :: units
    integer, optional,                  intent(in)    :: stock_index
    integer                                           :: MOM_generic_tracer_stock
  ! This function calculates the mass-weighted integral of all tracer stocks,
  ! returning the number of stocks it has calculated.  If the stock_index
  ! is present, only the stock corresponding to that coded index is returned.

  ! Arguments: h - Layer thickness, in m or kg m-2.
  !  (out)     stocks - the mass-weighted integrated amount of each tracer,
  !                     in kg times concentration units.
  !  (in)      G - The ocean's grid structure.
  !  (in)      GV - The ocean's vertical grid structure.
  !  (in)      CS - The control structure returned by a previous call to
  !                 register_MOM_generic_tracer.
  !  (out)     names - the names of the stocks calculated.
  !  (out)     units - the units of the stocks calculated.
  !  (in,opt)  stock_index - the coded index of a specific stock being sought.
  ! Return value: the number of stocks calculated here.
    type(g_tracer_type), pointer  :: g_tracer, g_tracer_next
    real, dimension(:,:,:,:), pointer   :: tr_field
    real, dimension(:,:,:), pointer     :: tr_ptr
    character(len=fm_string_len), parameter :: sub_name = 'MOM_generic_tracer_stock'

    integer :: i, j, k, is, ie, js, je, nz, m
    is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

    MOM_generic_tracer_stock = 0
    if (.not.associated(CS)) return

    if (present(stock_index)) then ; if (stock_index > 0) then
      ! Check whether this stock is available from this routine.

      ! No stocks from this routine are being checked yet.  Return 0.
      return
    endif ; endif

    if(.NOT. associated(CS%g_tracer_list)) return ! No stocks.

    m=1 ; g_tracer=>CS%g_tracer_list
    do
      call g_tracer_get_alias(g_tracer,names(m))
      call g_tracer_get_values(g_tracer,names(m),'units',units(m))
      units(m) = trim(units(m))//" kg"
      call g_tracer_get_pointer(g_tracer,names(m),'field',tr_field)

      stocks(m) = 0.0
      tr_ptr => tr_field(:,:,:,1)
      do k=1,nz ; do j=js,je ; do i=is,ie
        stocks(m) = stocks(m) + tr_ptr(i,j,k) * &
                               (G%mask2dT(i,j) * G%areaT(i,j) * h(i,j,k))
      enddo ; enddo ; enddo
      stocks(m) = GV%H_to_kg_m2 * stocks(m)

      !traverse the linked list till hit NULL
      call g_tracer_get_next(g_tracer, g_tracer_next)
      if(.NOT. associated(g_tracer_next)) exit
      g_tracer=>g_tracer_next
      m = m+1
    enddo

    MOM_generic_tracer_stock = m

  end function MOM_generic_tracer_stock

  ! <SUBROUTINE NAME="MOM_generic_tracer_min_max">
  !  <OVERVIEW>
  !   Find the min and max of tracer concentrations.
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !     This subroutine find the global min and max of either
  !   of all available tracer concentrations, or of a tracer that is
  !   being requested specifically, returning the number of tracers it has gone through.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   ns = MOM_generic_tracer_min_max(do_minmax, gmin, gmax, igmin, jgmin, kgmin, igmax, jgmax, kgmax , G, CS, names, units)
  !  </TEMPLATE>
  ! </SUBROUTINE>

  function MOM_generic_tracer_min_max(ind_start, got_minmax, gmin, gmax, xgmin, ygmin, zgmin, xgmax, ygmax, zgmax , G, CS, names, units)
    use mpp_utilities_mod, only: mpp_array_global_min_max
    integer,                            intent(in)    :: ind_start
    logical, dimension(:),              intent(out)   :: got_minmax
    real, dimension(:),                 intent(out)   :: gmin,gmax
    real, dimension(:),                 intent(out)   :: xgmin, ygmin, zgmin, xgmax, ygmax, zgmax
    type(ocean_grid_type),              intent(in)    :: G
    type(MOM_generic_tracer_CS),       pointer       :: CS
    character(len=*), dimension(:),     intent(out)   :: names
    character(len=*), dimension(:),     intent(out)   :: units
    integer                                           :: MOM_generic_tracer_min_max
  ! This function calculates the mass-weighted integral of all tracer stocks,
  ! returning the number of stocks it has calculated.  If the stock_index
  ! is present, only the stock corresponding to that coded index is returned.

  ! Arguments: h - Layer thickness, in m or kg m-2.
  !  (out)     gmin , gmax - the global minimum and maximum of each tracer,
  !                     in kg times concentration units.
  !  (in)      G - The ocean's grid structure.
  !  (in)      CS - The control structure returned by a previous call to
  !                 register_MOM_generic_tracer.
  !  (out)     names - the names of the stocks calculated.
  !  (out)     units - the units of the stocks calculated.
  !  (in,opt)  trace_index - the coded index of a specific tracer being sought.
  ! Return value: the number of tracers done here.

    type(g_tracer_type), pointer  :: g_tracer, g_tracer_next
    real, dimension(:,:,:,:), pointer   :: tr_field
    real, dimension(:,:,:), pointer     :: tr_ptr
    character(len=fm_string_len), parameter :: sub_name = 'MOM_generic_tracer_min_max'

    real, dimension(:,:,:),pointer :: grid_tmask
    integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau

    integer :: i, j, k, is, ie, js, je, nz, m
    real, allocatable, dimension(:) :: geo_z

    is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

    MOM_generic_tracer_min_max = 0
    if (.not.associated(CS)) return

    if(.NOT. associated(CS%g_tracer_list)) return ! No stocks.


    call g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,grid_tmask=grid_tmask)

    !  Because the use of a simple z-coordinate can not be assumed, simply
    ! use the layer index as the vertical label.
    allocate(geo_z(nk))
    do k=1,nk ; geo_z(k) = real(k) ; enddo


    m=ind_start ; g_tracer=>CS%g_tracer_list
    do
      call g_tracer_get_alias(g_tracer,names(m))
      call g_tracer_get_values(g_tracer,names(m),'units',units(m))
      units(m) = trim(units(m))//" kg"
      call g_tracer_get_pointer(g_tracer,names(m),'field',tr_field)

      gmin(m) = -1.0
      gmax(m) = -1.0

      tr_ptr => tr_field(:,:,:,1)


      call mpp_array_global_min_max(tr_ptr, grid_tmask,isd,jsd,isc,iec,jsc,jec,nk , gmin(m), gmax(m), &
                                    G%geoLonT,G%geoLatT,geo_z,xgmin(m), ygmin(m), zgmin(m), xgmax(m), ygmax(m), zgmax(m))

      got_minmax(m) = .true.


      !traverse the linked list till hit NULL
      call g_tracer_get_next(g_tracer, g_tracer_next)
      if(.NOT. associated(g_tracer_next)) exit
      g_tracer=>g_tracer_next
      m = m+1
    enddo

    MOM_generic_tracer_min_max = m

  end function MOM_generic_tracer_min_max


  ! <SUBROUTINE NAME="MOM_generic_tracer_surface_state">
  !  <OVERVIEW>
  !   Calculate the surface state and set coupler values
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   This subroutine calculates the surface state and set coupler values for
  !   those generic tracers that havd flux exchange with atmosphere.
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call MOM_generic_tracer_surface_state(state, h, G, CS)
  !  </TEMPLATE>
  ! </SUBROUTINE>

  subroutine MOM_generic_tracer_surface_state(state, h, G, CS)
    type(ocean_grid_type),                 intent(in) :: G
    type(surface),                         intent(inout) :: state
    real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h
    type(MOM_generic_tracer_CS),           pointer    :: CS
    !   This subroutine sets up the fields that the coupler needs to calculate the
    ! CFC fluxes between the ocean and atmosphere.
    ! Arguments: state - A structure containing fields that describe the
    !                    surface state of the ocean.
    !  (in)      h - Layer thickness, in m.
    !  (in)      G - The ocean's grid structure.
    !  (in)      CS - The control structure returned by a previous call to
    !                 register_MOM_generic_tracer.

    real :: sosga

    character(len=fm_string_len), parameter :: sub_name = 'MOM_generic_tracer_surface_state'
    real, dimension(G%isd:G%ied,G%jsd:G%jed,1:G%ke,1) :: rho0
    type(g_tracer_type), pointer :: g_tracer

    !Set coupler values
    !nnz: fake rho0
    rho0=1.0

    sosga = global_area_mean(state%SSS, G)

    call generic_tracer_coupler_set(state%tr_fields,&
         ST=state%SST,&
         SS=state%SSS,&
         rho=rho0,& !nnz: required for MOM5 and previous versions.
         ilb=G%isd, jlb=G%jsd,&
         tau=1,sosga=sosga,model_time=get_diag_time_end(CS%diag))

    !Output diagnostics via diag_manager for all tracers in this module
!    if(.NOT. associated(CS%g_tracer_list)) call mpp_error(FATAL, trim(sub_name)//&
!         "No tracer in the list.")
!    call g_tracer_send_diag(CS%g_tracer_list, get_diag_time_end(CS%diag), tau=1)
    !Niki: The problem with calling diagnostic outputs here is that this subroutine is called every dt_cpld
    !      hence if dt_therm > dt_cpld we get output (and contribution to the mean) at times that tracers
    !      had not been updated.
    !      Moving this to the end of column physics subrotuine fixes this issue.

  end subroutine MOM_generic_tracer_surface_state

!ALL PE subroutine on Ocean!  Due to otpm design the fluxes should be initialized like this on ALL PE's!
  subroutine MOM_generic_flux_init

    integer :: ind
    character(len=fm_string_len)   :: g_tracer_name,longname, package,units,old_package,file_in,file_out
    real :: const_init_value
    character(len=fm_string_len), parameter :: sub_name = 'MOM_generic_flux_init'
    type(g_tracer_type), pointer :: g_tracer_list,g_tracer,g_tracer_next

    if (.not. g_registered) then
       call generic_tracer_register
       g_registered = .true.
    endif

    call generic_tracer_get_list(g_tracer_list)
    if(.NOT. associated(g_tracer_list)) then
       call mpp_error(WARNING, trim(sub_name)// ": No generic tracer in the list.")
       return
    endif

    g_tracer=>g_tracer_list
    do

       call g_tracer_flux_init(g_tracer)

       !traverse the linked list till hit NULL
       call g_tracer_get_next(g_tracer, g_tracer_next)
       if(.NOT. associated(g_tracer_next)) exit
       g_tracer=>g_tracer_next

    enddo

  end subroutine MOM_generic_flux_init

  subroutine MOM_generic_tracer_fluxes_accumulate(flux_tmp, weight)
    type(forcing),         intent(in)    :: flux_tmp
    real,                  intent(in)    :: weight

   call generic_tracer_coupler_accumulate(flux_tmp%tr_fluxes, weight)

  end subroutine MOM_generic_tracer_fluxes_accumulate

  subroutine MOM_generic_tracer_get(name,member,array, CS)
    character(len=*),         intent(in)  :: name
    character(len=*),         intent(in)  :: member
    real, dimension(:,:,:),   intent(out) :: array
    type(MOM_generic_tracer_CS), pointer :: CS
    !  (in)      CS - The control structure returned by a previous call to
    !                 register_MOM_generic_tracer.

    real, dimension(:,:,:),   pointer :: array_ptr
    character(len=fm_string_len), parameter :: sub_name = 'MOM_generic_tracer_get'

    call g_tracer_get_pointer(CS%g_tracer_list,name,member,array_ptr)
    array(:,:,:) = array_ptr(:,:,:)

  end subroutine MOM_generic_tracer_get

  ! <SUBROUTINE NAME="end_MOM_generic_tracer">
  !  <OVERVIEW>
  !   Ends the generic tracer module
  !  </OVERVIEW>
  !  <DESCRIPTION>
  !   Call the end for generic tracer module and deallocate all temp arrays
  !  </DESCRIPTION>
  !  <TEMPLATE>
  !   call end_MOM_generic_tracer(CS)
  !  </TEMPLATE>
  ! </SUBROUTINE>

  subroutine end_MOM_generic_tracer(CS)
    type(MOM_generic_tracer_CS), pointer :: CS
    !   This subroutine deallocates the memory owned by this module.
    ! Argument: CS - The control structure returned by a previous call to
    !                register_MOM_generic_tracer.
    integer :: m

    call generic_tracer_end

    if (associated(CS)) then
       deallocate(CS)
    endif
  end subroutine end_MOM_generic_tracer

#endif /* _USE_GENERIC_TRACER */
end module MOM_generic_tracer
