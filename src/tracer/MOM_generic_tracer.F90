!> Drives the generic version of tracers TOPAZ and CFC and other GFDL BGC components
module MOM_generic_tracer

! This file is part of MOM6. See LICENSE.md for the license.

#include <MOM_memory.h>

! The following macro is usually defined in <fms_platform.h> but since MOM6 should not directly
! include files from FMS we replicate the macro lines here:
#ifdef NO_F2000
#define _ALLOCATED associated
#else
#define _ALLOCATED allocated
#endif

  ! ### These imports should not reach into FMS directly ###
  use field_manager_mod, only: fm_string_len

  use generic_tracer, only: generic_tracer_register, generic_tracer_get_diag_list
  use generic_tracer, only: generic_tracer_init, generic_tracer_source, generic_tracer_register_diag
  use generic_tracer, only: generic_tracer_coupler_get, generic_tracer_coupler_set
  use generic_tracer, only: generic_tracer_end, generic_tracer_get_list, do_generic_tracer
  use generic_tracer, only: generic_tracer_update_from_bottom,generic_tracer_vertdiff_G
  use generic_tracer, only: generic_tracer_coupler_accumulate

  use g_tracer_utils,   only: g_tracer_get_name,g_tracer_set_values,g_tracer_set_common,g_tracer_get_common
  use g_tracer_utils,   only: g_tracer_get_next,g_tracer_type,g_tracer_is_prog,g_tracer_flux_init
  use g_tracer_utils,   only: g_tracer_send_diag,g_tracer_get_values
  use g_tracer_utils,   only: g_tracer_get_pointer,g_tracer_get_alias,g_tracer_set_csdiag

  use MOM_ALE_sponge, only : set_up_ALE_sponge_field, ALE_sponge_CS
  use MOM_coms, only : max_across_PEs, min_across_PEs, PE_here
  use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
  use MOM_diag_mediator, only : diag_ctrl, get_diag_time_end
  use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE, is_root_pe
  use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
  use MOM_forcing_type, only : forcing, optics_type
  use MOM_grid, only : ocean_grid_type
  use MOM_hor_index, only : hor_index_type
  use MOM_io, only : file_exists, MOM_read_data, slasher
  use MOM_open_boundary, only : ocean_OBC_type
  use MOM_restart, only : register_restart_field, query_initialized, MOM_restart_CS
  use MOM_spatial_means, only : global_area_mean
  use MOM_sponge, only : set_up_sponge_field, sponge_CS
  use MOM_time_manager, only : time_type, set_time
  use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
  use MOM_tracer_registry, only : register_tracer, tracer_registry_type
  use MOM_tracer_Z_init, only : tracer_Z_init
  use MOM_tracer_initialization_from_Z, only : MOM_initialize_tracer_from_Z
  use MOM_unit_scaling, only : unit_scale_type
  use MOM_variables, only : surface, thermo_var_ptrs
  use MOM_verticalGrid, only : verticalGrid_type


  implicit none ; private

  !> An state hidden in module data that is very much not allowed in MOM6
  ! ### This needs to be fixed
  logical :: g_registered = .false.

  public register_MOM_generic_tracer, initialize_MOM_generic_tracer
  public MOM_generic_tracer_column_physics, MOM_generic_tracer_surface_state
  public end_MOM_generic_tracer, MOM_generic_tracer_get
  public MOM_generic_tracer_stock
  public MOM_generic_flux_init
  public MOM_generic_tracer_min_max
  public MOM_generic_tracer_fluxes_accumulate

  !> Control structure for generic tracers
  type, public :: MOM_generic_tracer_CS ; private
    character(len = 200) :: IC_file !< The file in which the generic tracer initial values can
                                    !! be found, or an empty string for internal initialization.
    logical :: Z_IC_file !< If true, the generic_tracer IC_file is in Z-space.  The default is false.
    real :: tracer_IC_val = 0.0    !< The initial value assigned to tracers.
    real :: tracer_land_val = -1.0 !< The values of tracers used where  land is masked out.
    logical :: tracers_may_reinit  !< If true, tracers may go through the
                                   !! initialization code if they are not found in the restart files.

    type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                               !! regulate the timing of diagnostic output.
    type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< Restart control structure

    !> Pointer to the first element of the linked list of generic tracers.
    type(g_tracer_type), pointer :: g_tracer_list => NULL()

    integer :: H_to_m !< Auxiliary to access GV%H_to_m in routines that do not have access to GV

  end type MOM_generic_tracer_CS

! This include declares and sets the variable "version".
#include "version_variable.h"

contains

  !> Initializes the generic tracer packages and adds their tracers to the list
  !! Adds the tracers in the list of generic tracers to the set of MOM tracers (i.e., MOM-register them)
  !! Register these tracers for restart
  function register_MOM_generic_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
    type(hor_index_type),       intent(in)   :: HI         !< Horizontal index ranges
    type(verticalGrid_type),    intent(in)   :: GV         !< The ocean's vertical grid structure
    type(param_file_type),      intent(in)   :: param_file !< A structure to parse for run-time parameters
    type(MOM_generic_tracer_CS), pointer     :: CS         !< Pointer to the control structure for this module
    type(tracer_registry_type), pointer      :: tr_Reg     !< Pointer to the control structure for the tracer
                                                           !! advection and diffusion module.
    type(MOM_restart_CS),       pointer      :: restart_CS !< Pointer to the restart control structure.

! Local variables
    logical :: register_MOM_generic_tracer

    character(len=128), parameter :: sub_name = 'register_MOM_generic_tracer'
    character(len=200) :: inputdir ! The directory where NetCDF input files are.
    ! These can be overridden later in via the field manager?

    integer :: ntau, k,i,j,axes(3)
    type(g_tracer_type), pointer      :: g_tracer,g_tracer_next
    character(len=fm_string_len)      :: g_tracer_name,longname,units
    real, dimension(:,:,:,:), pointer   :: tr_field
    real, dimension(:,:,:), pointer     :: tr_ptr
    real, dimension(HI%isd:HI%ied, HI%jsd:HI%jed,GV%ke)         :: grid_tmask
    integer, dimension(HI%isd:HI%ied, HI%jsd:HI%jed)           :: grid_kmt

    register_MOM_generic_tracer = .false.
    if (associated(CS)) then
      call MOM_error(WARNING, "register_MOM_generic_tracer called with an "// &
            "associated control structure.")
      return
    endif
    allocate(CS)


    !Register all the generic tracers used and create the list of them.
    !This can be called by ALL PE's. No array fields allocated.
    if (.not. g_registered) then
       call generic_tracer_register()
       g_registered = .true.
    endif


  ! Read all relevant parameters and write them to the model log.
    call log_version(param_file, sub_name, version, "")
    call get_param(param_file, sub_name, "GENERIC_TRACER_IC_FILE", CS%IC_file, &
                 "The file in which the generic trcer initial values can "//&
                 "be found, or an empty string for internal initialization.", &
                 default=" ")
    if ((len_trim(CS%IC_file) > 0) .and. (scan(CS%IC_file,'/') == 0)) then
      ! Add the directory if CS%IC_file is not already a complete path.
      call get_param(param_file, sub_name, "INPUTDIR", inputdir, default=".")
      CS%IC_file = trim(slasher(inputdir))//trim(CS%IC_file)
      call log_param(param_file, sub_name, "INPUTDIR/GENERIC_TRACER_IC_FILE", CS%IC_file)
    endif
    call get_param(param_file, sub_name, "GENERIC_TRACER_IC_FILE_IS_Z", CS%Z_IC_file, &
                 "If true, GENERIC_TRACER_IC_FILE is in depth space, not "//&
                 "layer space.",default=.false.)
    call get_param(param_file, sub_name, "TRACERS_MAY_REINIT", CS%tracers_may_reinit, &
                 "If true, tracers may go through the initialization code "//&
                 "if they are not found in the restart files.  Otherwise "//&
                 "it is a fatal error if tracers are not found in the "//&
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
    if (.NOT. associated(CS%g_tracer_list)) call MOM_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")
    ! For each tracer name get its T_prog index and get its fields

    g_tracer=>CS%g_tracer_list
    do
       call g_tracer_get_alias(g_tracer,g_tracer_name)

       call g_tracer_get_pointer(g_tracer,g_tracer_name,'field',tr_field)
       call g_tracer_get_values(g_tracer,g_tracer_name,'longname', longname)
       call g_tracer_get_values(g_tracer,g_tracer_name,'units',units )

       !!nnz: MOM field is 3D. Does this affect performance? Need it be override field?
       tr_ptr => tr_field(:,:,:,1)
       ! Register prognastic tracer for horizontal advection, diffusion, and restarts.
       if (g_tracer_is_prog(g_tracer)) then
         call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, &
                              name=g_tracer_name, longname=longname, units=units, &
                              registry_diags=.false., &   !### CHANGE TO TRUE?
                              restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit)
       else
         call register_restart_field(tr_ptr, g_tracer_name, .not.CS%tracers_may_reinit, &
                                     restart_CS, longname=longname, units=units)
       endif

       !traverse the linked list till hit NULL
       call g_tracer_get_next(g_tracer, g_tracer_next)
       if (.NOT. associated(g_tracer_next)) exit
       g_tracer=>g_tracer_next

    enddo

    register_MOM_generic_tracer = .true.
  end function register_MOM_generic_tracer

  !>  Initialize phase II:  Initialize required variables for generic tracers
  !!  There are some steps of initialization that cannot be done in register_MOM_generic_tracer
  !!  This is the place and time to do them:
  !!      Set the grid mask and initial time for all generic tracers.
  !!      Diag_register them.
  !!      Z_diag_register them.
  !!
  !!   This subroutine initializes the NTR tracer fields in tr(:,:,:,:)
  !! and it sets up the tracer output.
  subroutine initialize_MOM_generic_tracer(restart, day, G, GV, US, h, param_file, diag, OBC, CS, &
                                          sponge_CSp, ALE_sponge_CSp)
    logical,                               intent(in) :: restart !< .true. if the fields have already been
                                                                 !! read from a restart file.
    type(time_type), target,               intent(in) :: day     !< Time of the start of the run.
    type(ocean_grid_type),                 intent(inout) :: G    !< The ocean's grid structure
    type(verticalGrid_type),               intent(in)    :: GV   !< The ocean's vertical grid structure
    type(unit_scale_type),                 intent(in)    :: US   !< A dimensional unit scaling type
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
    type(param_file_type),                 intent(in) :: param_file !< A structure to parse for run-time parameters
    type(diag_ctrl),               target, intent(in) :: diag    !< Regulates diagnostic output.
    type(ocean_OBC_type),                  pointer    :: OBC     !< This open boundary condition type specifies whether,
                                                                 !! where, and what open boundary conditions are used.
    type(MOM_generic_tracer_CS),           pointer    :: CS      !< Pointer to the control structure for this module.
    type(sponge_CS),                       pointer    :: sponge_CSp !< Pointer to the control structure for the sponges.
    type(ALE_sponge_CS),                   pointer    :: ALE_sponge_CSp !< Pointer  to the control structure for the
                                                                 !! ALE sponges.

    character(len=128), parameter :: sub_name = 'initialize_MOM_generic_tracer'
    logical :: OK
    integer :: i, j, k, isc, iec, jsc, jec, nk
    type(g_tracer_type), pointer    :: g_tracer,g_tracer_next
    character(len=fm_string_len)      :: g_tracer_name
    real, dimension(:,:,:,:), pointer   :: tr_field
    real, dimension(:,:,:), pointer     :: tr_ptr
    real,    dimension(G%isd:G%ied, G%jsd:G%jed, 1:GV%ke) :: grid_tmask
    integer, dimension(G%isd:G%ied, G%jsd:G%jed)          :: grid_kmt

    !! 2010/02/04  Add code to re-initialize Generic Tracers if needed during a model simulation
    !! By default, restart cpio should not contain a Generic Tracer IC file and step below will be skipped.
    !! Ideally, the generic tracer IC file should have the tracers on Z levels.

    isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; nk = GV%ke

    CS%diag=>diag
    !Get the tracer list
    if (.NOT. associated(CS%g_tracer_list)) call MOM_error(FATAL, trim(sub_name)//&
         ": No tracer in the list.")
    !For each tracer name get its  fields
    g_tracer=>CS%g_tracer_list

    do
      if (INDEX(CS%IC_file, '_NULL_') /= 0) then
        call MOM_error(WARNING, "The name of the IC_file "//trim(CS%IC_file)//&
                              " indicates no MOM initialization was asked for the generic tracers."//&
                              "Bypassing the MOM initialization of ALL generic tracers!")
        exit
      endif
      call g_tracer_get_alias(g_tracer,g_tracer_name)
      call g_tracer_get_pointer(g_tracer,g_tracer_name,'field',tr_field)
      tr_ptr => tr_field(:,:,:,1)

      if (.not.restart .or. (CS%tracers_may_reinit .and. &
          .not.query_initialized(tr_ptr, g_tracer_name, CS%restart_CSp))) then

       if (g_tracer%requires_src_info ) then
         call MOM_error(NOTE,"initialize_MOM_generic_tracer: "//&
                             "initializing generic tracer "//trim(g_tracer_name)//&
                             " using MOM_initialize_tracer_from_Z ")

         call MOM_initialize_tracer_from_Z(h, tr_ptr, G, GV, US, param_file,               &
                                src_file = g_tracer%src_file,                              &
                                src_var_nam = g_tracer%src_var_name,                       &
                                src_var_unit_conversion = g_tracer%src_var_unit_conversion,&
                                src_var_record = g_tracer%src_var_record,                  &
                                src_var_gridspec = g_tracer%src_var_gridspec               )

         !Check/apply the bounds for each g_tracer
         do k=1,nk ; do j=jsc,jec ; do i=isc,iec
           if (tr_ptr(i,j,k) /= CS%tracer_land_val) then
             if (tr_ptr(i,j,k) < g_tracer%src_var_valid_min) tr_ptr(i,j,k) = g_tracer%src_var_valid_min
             !Jasmin does not want to apply the maximum for now
             !if (tr_ptr(i,j,k) > g_tracer%src_var_valid_max) tr_ptr(i,j,k) = g_tracer%src_var_valid_max
           endif
         enddo ; enddo ; enddo

         !jgj: Reset CASED to 0 below K=1
         if ( (trim(g_tracer_name) == 'cased') .or. (trim(g_tracer_name) == 'ca13csed') ) then
           do k=2,nk ; do j=jsc,jec ; do i=isc,iec
             if (tr_ptr(i,j,k) /= CS%tracer_land_val) then
               tr_ptr(i,j,k) = 0.0
             endif
           enddo ; enddo ; enddo
         endif
       elseif(.not. g_tracer%requires_restart) then
         !Do nothing for this tracer, it is initialized by the tracer package
          call MOM_error(NOTE,"initialize_MOM_generic_tracer: "//&
                            "skip initialization of generic tracer "//trim(g_tracer_name))
       else !Do it old way if the tracer is not registered to start from a specific source file.
            !This path should be deprecated if all generic tracers are required to start from specified sources.
        if (len_trim(CS%IC_file) > 0) then
        !  Read the tracer concentrations from a netcdf file.
          if (.not.file_exists(CS%IC_file)) call MOM_error(FATAL, &
                  "initialize_MOM_Generic_tracer: Unable to open "//CS%IC_file)
          if (CS%Z_IC_file) then
            OK = tracer_Z_init(tr_ptr, h, CS%IC_file, g_tracer_name, G, GV, US)
            if (.not.OK) then
              OK = tracer_Z_init(tr_ptr, h, CS%IC_file, trim(g_tracer_name), G, GV, US)
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
            call MOM_read_data(CS%IC_file, trim(g_tracer_name), tr_ptr, G%Domain)
          endif
        else
          call MOM_error(FATAL,"initialize_MOM_generic_tracer: "//&
                  "check Generic Tracer IC filename "//trim(CS%IC_file)//&
                  " for tracer "//trim(g_tracer_name))
        endif

       endif
      endif

      !traverse the linked list till hit NULL
      call g_tracer_get_next(g_tracer, g_tracer_next)
      if (.NOT. associated(g_tracer_next)) exit
      g_tracer=>g_tracer_next
    enddo
    !! end section to re-initialize generic tracers


    !Now we can reset the grid mask, axes and time to their true values
    !Note that grid_tmask must be set correctly on the data domain boundary
    !so that coast mask can be deduced from it.
    grid_tmask(:,:,:) = 0.0
    grid_kmt(:,:) = 0
    do j = G%jsd, G%jed ; do i = G%isd, G%ied
      if (G%mask2dT(i,j) > 0) then
        grid_tmask(i,j,:) = 1.0
        grid_kmt(i,j) = GV%ke ! Tell the code that a layer thicker than 1m is the bottom layer.
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

    CS%H_to_m = GV%H_to_m

  end subroutine initialize_MOM_generic_tracer

  !>  Column physics for generic tracers.
  !!      Get the coupler values for generic tracers that exchange with atmosphere
  !!      Update generic tracer concentration fields from sources and sinks.
  !!      Vertically diffuse generic tracer concentration fields.
  !!      Update generic tracers from bottom and their bottom reservoir.
  !!
  !!   This subroutine applies diapycnal diffusion and any other column
  !! tracer physics or chemistry to the tracers from this file.
  !! CFCs are relatively simple, as they are passive tracers. with only a surface
  !! flux as a source.
  subroutine MOM_generic_tracer_column_physics(h_old, h_new, ea, eb, fluxes, Hml, dt, G, GV, CS, tv, optics, &
        evap_CFL_limit, minimum_forcing_depth)
    type(ocean_grid_type),   intent(in) :: G     !< The ocean's grid structure
    type(verticalGrid_type), intent(in) :: GV    !< The ocean's vertical grid structure
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(in) :: ea    !< The amount of fluid entrained from the layer
                                                 !! above during this call [H ~> m or kg m-2].
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(in) :: eb    !< The amount of fluid entrained from the layer
                                                 !! below during this call [H ~> m or kg m-2].
    type(forcing),           intent(in) :: fluxes !< A structure containing pointers to thermodynamic
                                                 !! and tracer forcing fields.
    real, dimension(SZI_(G),SZJ_(G)), intent(in) :: Hml  !< Mixed layer depth [Z ~> m]
    real,                    intent(in) :: dt    !< The amount of time covered by this call [s]
    type(MOM_generic_tracer_CS), pointer :: CS   !< Pointer to the control structure for this module.
    type(thermo_var_ptrs),   intent(in) :: tv    !< A structure pointing to various thermodynamic variables
    type(optics_type),       intent(in) :: optics !< The structure containing optical properties.
    real,          optional, intent(in) :: evap_CFL_limit !< Limits how much water can be fluxed out of
                                                 !! the top layer Stored previously in diabatic CS.
    real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which fluxes
                                                 !!  can be applied [H ~> m or kg m-2]
                                                 !   Stored previously in diabatic CS.
    ! The arguments to this subroutine are redundant in that
    !     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)

    ! Local variables
    character(len=128), parameter :: sub_name = 'MOM_generic_tracer_column_physics'

    type(g_tracer_type), pointer  :: g_tracer, g_tracer_next
    character(len=fm_string_len)  :: g_tracer_name
    real, dimension(:,:), pointer :: stf_array,trunoff_array,runoff_tracer_flux_array

    real :: surface_field(SZI_(G),SZJ_(G))
    real :: dz_ml(SZI_(G),SZJ_(G))  ! The mixed layer depth in the MKS units used for generic tracers [m]
    real :: sosga

    real, dimension(G%isd:G%ied,G%jsd:G%jed,GV%ke) :: rho_dzt, dzt
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV))      :: h_work
    integer :: i, j, k, isc, iec, jsc, jec, nk

    isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; nk = GV%ke

    !Get the tracer list
    if (.NOT. associated(CS%g_tracer_list)) call MOM_error(FATAL,&
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
       if (_ALLOCATED(g_tracer%trunoff)) then
          call g_tracer_get_alias(g_tracer,g_tracer_name)
          call g_tracer_get_pointer(g_tracer,g_tracer_name,'stf',   stf_array)
          call g_tracer_get_pointer(g_tracer,g_tracer_name,'trunoff',trunoff_array)
          call g_tracer_get_pointer(g_tracer,g_tracer_name,'runoff_tracer_flux',runoff_tracer_flux_array)
          !nnz: Why is fluxes%river = 0?
          runoff_tracer_flux_array(:,:) = trunoff_array(:,:) * &
                   G%US%RZ_T_to_kg_m2s*fluxes%lrunoff(:,:)
          stf_array = stf_array + runoff_tracer_flux_array
       endif

       !traverse the linked list till hit NULL
       call g_tracer_get_next(g_tracer, g_tracer_next)
       if (.NOT. associated(g_tracer_next)) exit
       g_tracer=>g_tracer_next

    enddo

    !
    !Prepare input arrays for source update
    !

    rho_dzt(:,:,:) = GV%H_to_kg_m2 * GV%Angstrom_H
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
      rho_dzt(i,j,k) = GV%H_to_kg_m2 * h_old(i,j,k)
    enddo ; enddo ; enddo !}

    dzt(:,:,:) = 1.0
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec  !{
      dzt(i,j,k) = GV%H_to_m * h_old(i,j,k)
    enddo ; enddo ; enddo !}
    dz_ml(:,:) = 0.0
    do j=jsc,jec ; do i=isc,iec
      surface_field(i,j) = tv%S(i,j,1)
      dz_ml(i,j) = G%US%Z_to_m * Hml(i,j)
    enddo ; enddo
    sosga = global_area_mean(surface_field, G)

    !
    !Calculate tendencies (i.e., field changes at dt) from the sources / sinks
    !
    if ((G%US%L_to_m == 1.0) .and. (G%US%RZ_to_kg_m2 == 1.0) .and. (G%US%s_to_T == 1.0)) then
      ! Avoid unnecessary copies when no unit conversion is needed.
      call generic_tracer_source(tv%T, tv%S, rho_dzt, dzt, dz_ml, G%isd, G%jsd, 1, dt, &
               G%areaT, get_diag_time_end(CS%diag), &
               optics%nbands, optics%max_wavelength_band, optics%sw_pen_band, optics%opacity_band, &
               internal_heat=tv%internal_heat, frunoff=fluxes%frunoff, sosga=sosga)
    else
      call generic_tracer_source(tv%T, tv%S, rho_dzt, dzt, dz_ml, G%isd, G%jsd, 1, dt, &
               G%US%L_to_m**2*G%areaT(:,:), get_diag_time_end(CS%diag), &
               optics%nbands, optics%max_wavelength_band, optics%sw_pen_band, optics%opacity_band, &
               internal_heat=G%US%RZ_to_kg_m2*tv%internal_heat(:,:), &
               frunoff=G%US%RZ_T_to_kg_m2s*fluxes%frunoff(:,:), sosga=sosga)
    endif

    ! This uses applyTracerBoundaryFluxesInOut to handle the change in tracer due to freshwater fluxes
    ! usually in ALE mode
    if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
      g_tracer=>CS%g_tracer_list
      do
        if (g_tracer_is_prog(g_tracer)) then
          do k=1,nk ;do j=jsc,jec ; do i=isc,iec
            h_work(i,j,k) = h_old(i,j,k)
          enddo ; enddo ; enddo
          call applyTracerBoundaryFluxesInOut(G, GV, g_tracer%field(:,:,:,1), G%US%s_to_T*dt, &
                            fluxes, h_work, evap_CFL_limit, minimum_forcing_depth)
        endif

         !traverse the linked list till hit NULL
         call g_tracer_get_next(g_tracer, g_tracer_next)
        if (.NOT. associated(g_tracer_next)) exit
        g_tracer=>g_tracer_next
      enddo
    endif

    !
    !Update Tr(n)%field from explicit vertical diffusion
    !
    ! Use a tridiagonal solver to determine the concentrations after the
    ! surface source is applied and diapycnal advection and diffusion occurs.
    if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
      ! Last arg is tau which is always 1 for MOM6
      call generic_tracer_vertdiff_G(h_work, ea, eb, dt, GV%kg_m2_to_H, GV%m_to_H, 1)
    else
      ! Last arg is tau which is always 1 for MOM6
      call generic_tracer_vertdiff_G(h_old, ea, eb, dt, GV%kg_m2_to_H, GV%m_to_H, 1)
    endif

    ! Update bottom fields after vertical processes

    ! Second arg is tau which is always 1 for MOM6
    call generic_tracer_update_from_bottom(dt, 1, get_diag_time_end(CS%diag))

    !Output diagnostics via diag_manager for all generic tracers and their fluxes
    call g_tracer_send_diag(CS%g_tracer_list, get_diag_time_end(CS%diag), tau=1)
#ifdef _USE_MOM6_DIAG
    call g_tracer_set_csdiag(CS%diag)
#endif

  end subroutine MOM_generic_tracer_column_physics

  !> This subroutine calculates mass-weighted integral on the PE either
  !! of all available tracer concentrations, or of a tracer that is
  !! being requested specifically, returning the number of stocks it has
  !! calculated. If the stock_index is present, only the stock corresponding
  !! to that coded index is returned.
  function MOM_generic_tracer_stock(h, stocks, G, GV, CS, names, units, stock_index)
    type(ocean_grid_type),              intent(in)    :: G    !< The ocean's grid structure
    type(verticalGrid_type),            intent(in)    :: GV   !< The ocean's vertical grid structure
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< Layer thicknesses [H ~> m or kg m-2]
    real, dimension(:),                 intent(out)   :: stocks !< The mass-weighted integrated amount of each
                                                                !! tracer, in kg times concentration units [kg conc].
    type(MOM_generic_tracer_CS),        pointer       :: CS     !< Pointer to the control structure for this module.
    character(len=*), dimension(:),     intent(out)   :: names  !< The names of the stocks calculated.
    character(len=*), dimension(:),     intent(out)   :: units  !< The units of the stocks calculated.
    integer, optional,                  intent(in)    :: stock_index !< The coded index of a specific stock
                                                                     !! being sought.
    integer                                           :: MOM_generic_tracer_stock !< Return value, the
                                                                     !! number of stocks calculated here.

    ! Local variables
    real :: stock_scale ! The dimensional scaling factor to convert stocks to kg [kg H-1 L-2 ~> kg m-3 or nondim]
    type(g_tracer_type), pointer  :: g_tracer, g_tracer_next
    real, dimension(:,:,:,:), pointer   :: tr_field
    real, dimension(:,:,:), pointer     :: tr_ptr
    character(len=128), parameter :: sub_name = 'MOM_generic_tracer_stock'

    integer :: i, j, k, is, ie, js, je, nz, m
    is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

    MOM_generic_tracer_stock = 0
    if (.not.associated(CS)) return

    if (present(stock_index)) then ; if (stock_index > 0) then
      ! Check whether this stock is available from this routine.

      ! No stocks from this routine are being checked yet.  Return 0.
      return
    endif ; endif

    if (.NOT. associated(CS%g_tracer_list)) return ! No stocks.

    stock_scale = G%US%L_to_m**2 * GV%H_to_kg_m2
    m=1 ; g_tracer=>CS%g_tracer_list
    do
      call g_tracer_get_alias(g_tracer,names(m))
      call g_tracer_get_values(g_tracer,names(m),'units',units(m))
      units(m) = trim(units(m))//" kg"
      call g_tracer_get_pointer(g_tracer,names(m),'field',tr_field)

      stocks(m) = 0.0
      tr_ptr => tr_field(:,:,:,1)
      do k=1,nz ; do j=js,je ; do i=is,ie
        stocks(m) = stocks(m) + tr_ptr(i,j,k) * (G%mask2dT(i,j) * G%areaT(i,j) * h(i,j,k))
      enddo ; enddo ; enddo
      stocks(m) = stock_scale * stocks(m)

      !traverse the linked list till hit NULL
      call g_tracer_get_next(g_tracer, g_tracer_next)
      if (.NOT. associated(g_tracer_next)) exit
      g_tracer=>g_tracer_next
      m = m+1
    enddo

    MOM_generic_tracer_stock = m

  end function MOM_generic_tracer_stock

  !> This subroutine find the global min and max of either of all
  !! available tracer concentrations, or of a tracer that is being
  !! requested specifically, returning the number of tracers it has gone through.
  function MOM_generic_tracer_min_max(ind_start, got_minmax, gmin, gmax, xgmin, ygmin, zgmin, &
                                      xgmax, ygmax, zgmax , G, CS, names, units)
    integer,                        intent(in)    :: ind_start !< The index of the tracer to start with
    logical, dimension(:),          intent(out)   :: got_minmax !< Indicates whether the global min and
                                                            !! max are found for each tracer
    real, dimension(:),             intent(out)   :: gmin   !< Global minimum of each tracer, in kg
                                                            !! times concentration units.
    real, dimension(:),             intent(out)   :: gmax   !< Global maximum of each tracer, in kg
                                                            !! times concentration units.
    real, dimension(:),             intent(out)   :: xgmin  !< The x-position of the global minimum
    real, dimension(:),             intent(out)   :: ygmin  !< The y-position of the global minimum
    real, dimension(:),             intent(out)   :: zgmin  !< The z-position of the global minimum
    real, dimension(:),             intent(out)   :: xgmax  !< The x-position of the global maximum
    real, dimension(:),             intent(out)   :: ygmax  !< The y-position of the global maximum
    real, dimension(:),             intent(out)   :: zgmax  !< The z-position of the global maximum
    type(ocean_grid_type),          intent(in)    :: G      !< The ocean's grid structure
    type(MOM_generic_tracer_CS),    pointer       :: CS     !< Pointer to the control structure for this module.
    character(len=*), dimension(:), intent(out)   :: names  !< The names of the stocks calculated.
    character(len=*), dimension(:), intent(out)   :: units  !< The units of the stocks calculated.
    integer                                       :: MOM_generic_tracer_min_max !< Return value, the
                                                            !! number of tracers done here.

! Local variables
    type(g_tracer_type), pointer  :: g_tracer, g_tracer_next
    real, dimension(:,:,:,:), pointer   :: tr_field
    real, dimension(:,:,:), pointer     :: tr_ptr
    character(len=128), parameter :: sub_name = 'MOM_generic_tracer_min_max'

    real, dimension(:,:,:),pointer :: grid_tmask
    integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau

    integer :: i, j, k, is, ie, js, je, m
    real, allocatable, dimension(:) :: geo_z

    is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

    MOM_generic_tracer_min_max = 0
    if (.not.associated(CS)) return

    if (.NOT. associated(CS%g_tracer_list)) return ! No stocks.


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

      call array_global_min_max(tr_ptr, grid_tmask, isd, jsd, isc, iec, jsc, jec, nk, gmin(m), gmax(m), &
                                    G%geoLonT, G%geoLatT, geo_z, xgmin(m), ygmin(m), zgmin(m), &
                                    xgmax(m), ygmax(m), zgmax(m))

      got_minmax(m) = .true.

      !traverse the linked list till hit NULL
      call g_tracer_get_next(g_tracer, g_tracer_next)
      if (.NOT. associated(g_tracer_next)) exit
      g_tracer=>g_tracer_next
      m = m+1
    enddo

    MOM_generic_tracer_min_max = m

  end function MOM_generic_tracer_min_max

  !> Find the global maximum and minimum of a tracer array and return the locations of the extrema.
  subroutine array_global_min_max(tr_array, tmask, isd, jsd, isc, iec, jsc, jec, nk, g_min, g_max, &
                                  geo_x, geo_y, geo_z, xgmin, ygmin, zgmin, xgmax, ygmax, zgmax)
    integer,                      intent(in)  :: isd   !< The starting data domain i-index
    integer,                      intent(in)  :: jsd   !< The starting data domain j-index
    real, dimension(isd:,jsd:,:), intent(in)  :: tr_array !< The tracer array to search for extrema
    real, dimension(isd:,jsd:,:), intent(in)  :: tmask !< A mask that is 0 for points to exclude
    integer,                      intent(in)  :: isc   !< The starting compute domain i-index
    integer,                      intent(in)  :: iec   !< The ending compute domain i-index
    integer,                      intent(in)  :: jsc   !< The starting compute domain j-index
    integer,                      intent(in)  :: jec   !< The ending compute domain j-index
    integer,                      intent(in)  :: nk    !< The number of vertical levels
    real,                         intent(out) :: g_min !< The global minimum of tr_array
    real,                         intent(out) :: g_max !< The global maximum of tr_array
    real, dimension(isd:,jsd:),   intent(in)  :: geo_x !< The geographic x-positions of points
    real, dimension(isd:,jsd:),   intent(in)  :: geo_y !< The geographic y-positions of points
    real, dimension(:),           intent(in)  :: geo_z !< The vertical pseudo-positions of points
    real,                         intent(out) :: xgmin !< The x-position of the global minimum
    real,                         intent(out) :: ygmin !< The y-position of the global minimum
    real,                         intent(out) :: zgmin !< The z-position of the global minimum
    real,                         intent(out) :: xgmax !< The x-position of the global maximum
    real,                         intent(out) :: ygmax !< The y-position of the global maximum
    real,                         intent(out) :: zgmax !< The z-position of the global maximum

    ! This subroutine is an exact transcription (bugs and all) of mpp_array_global_min_max()
    ! from the version in FMS/mpp/mpp_utilities.F90, but with some whitespace changes to match
    ! MOM6 code styles and to use infrastructure routines via the MOM6 framework code, and with
    ! added comments to document its arguments.i

    !### The obvious problems with this routine as currently written include:
    !  1. It does not return exactly the maximum and minimum values.
    !  2. The reported maximum and minimum are dependent on PE count and layout.
    !  3. For all-zero arrays, the reported maxima scale with the PE_count
    !  4. For arrays with a large enough offset or scaling, so that the magnitude of values exceed
    !     1e10, the values it returns are simply wrong.
    !  5. The results do not scale appropriately if the argument is rescaled.
    !  6. The extrema and locations are not rotationally invariant.
    !  7. It is inefficient because it uses 8 blocking global reduction calls when it could use just 2 or 3.

    ! Local variables
    real    :: tmax, tmin   ! Maximum and minimum tracer values, in the same units as tr_array
    real    :: tmax0, tmin0 ! First-guest values of tmax and tmin.
    integer :: itmax, jtmax, ktmax, itmin, jtmin, ktmin
    integer :: igmax, jgmax, kgmax, igmin, jgmin, kgmin
    real    :: fudge ! A factor that is close to 1 that is used to find the location of the extrema.

     ! arrays to enable vectorization
    integer :: iminarr(3), imaxarr(3)

    !### These dimensional constant values mean that the results can not be guaranteed to be rescalable.
    g_min = -88888888888.0 ; g_max = -999999999.0
    tmax = -1.e10 ; tmin = 1.e10
    itmax = 0 ; jtmax = 0 ; ktmax = 0
    itmin = 0 ; jtmin = 0 ; ktmin = 0

    if (ANY(tmask(isc:iec,jsc:jec,:) > 0.)) then
      ! Vectorized using maxloc() and minloc() intrinsic functions by Russell.Fiedler@csiro.au.
      iminarr = minloc(tr_array(isc:iec,jsc:jec,:), (tmask(isc:iec,jsc:jec,:) > 0.))
      imaxarr = maxloc(tr_array(isc:iec,jsc:jec,:), (tmask(isc:iec,jsc:jec,:) > 0.))
      itmin = iminarr(1)+isc-1
      jtmin = iminarr(2)+jsc-1
      ktmin = iminarr(3)
      itmax = imaxarr(1)+isc-1
      jtmax = imaxarr(2)+jsc-1
      ktmax = imaxarr(3)
      tmin = tr_array(itmin,jtmin,ktmin)
      tmax = tr_array(itmax,jtmax,ktmax)
    end if

    ! use "fudge" to distinguish processors when tracer extreme is independent of processor
    !### This fudge factor is not independent of PE layout, and while it mostly works for finding
    !    a positive maximum or a negative minimum, it could miss the true extrema in the opposite
    !    cases, for which the fudge factor should be slightly reduced.  The fudge factor should
    !    be based on global index-space conventions, which are decomposition invariant, and
    !    not the PE-number!
    fudge = 1.0 + 1.e-12*real(PE_here() )
    tmax = tmax*fudge
    tmin = tmin*fudge
    if (tmax == 0.0) then
      tmax = tmax + 1.e-12*real(PE_here() )
    endif
    if (tmin == 0.0) then
      tmin = tmin + 1.e-12*real(PE_here() )
    endif

    tmax0 = tmax ; tmin0 = tmin

    call max_across_PEs(tmax)
    call min_across_PEs(tmin)

    g_max = tmax
    g_min = tmin

    ! Now find the location of the global extrema.
    !
    ! Note that the fudge factor above guarantees that the location of max (min) is uinque,
    ! since tmax0 (tmin0) has slightly different values on each processor.
    ! Otherwise, the function tr_array(i,j,k) could be equal to global max (min) at more
    ! than one point in space and this would be a much more difficult problem to solve.
    !
    !-999 on all current PE's
    xgmax = -999. ; ygmax = -999. ; zgmax = -999.
    xgmin = -999. ; ygmin = -999. ; zgmin = -999.

    if (tmax0 == tmax) then !This happens ONLY on ONE processor because of fudge factor above.
      xgmax = geo_x(itmax,jtmax)
      ygmax = geo_y(itmax,jtmax)
      zgmax = geo_z(ktmax)
    endif

    !### These three calls and the three calls that follow in about 10 lines should be combined
    !    into a single call for efficiency.
    call max_across_PEs(xgmax)
    call max_across_PEs(ygmax)
    call max_across_PEs(zgmax)

    if (tmin0 == tmin) then !This happens ONLY on ONE processor because of fudge factor above.
      xgmin = geo_x(itmin,jtmin)
      ygmin = geo_y(itmin,jtmin)
      zgmin = geo_z(ktmin)
    endif

    call max_across_PEs(xgmin)
    call max_across_PEs(ygmin)
    call max_across_PEs(zgmin)

  end subroutine array_global_min_max

  !> This subroutine calculates the surface state and sets coupler values for
  !! those generic tracers that have flux exchange with atmosphere.
  !!
  !! This subroutine sets up the fields that the coupler needs to calculate the
  !! CFC fluxes between the ocean and atmosphere.
  subroutine MOM_generic_tracer_surface_state(sfc_state, h, G, GV, CS)
    type(ocean_grid_type),                 intent(in)    :: G    !< The ocean's grid structure
    type(verticalGrid_type),               intent(in)    :: GV   !< The ocean's vertical grid structure
    type(surface),                         intent(inout) :: sfc_state !< A structure containing fields that
                                                                 !! describe the surface state of the ocean.
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
    type(MOM_generic_tracer_CS),           pointer       :: CS   !< Pointer to the control structure for this module.

! Local variables
    real :: sosga

    character(len=128), parameter :: sub_name = 'MOM_generic_tracer_surface_state'
    real, dimension(G%isd:G%ied,G%jsd:G%jed,1:GV%ke,1) :: rho0
    real, dimension(G%isd:G%ied,G%jsd:G%jed,1:GV%ke) ::  dzt
    type(g_tracer_type), pointer :: g_tracer

    !Set coupler values
    !nnz: fake rho0
    rho0=1.0

    dzt(:,:,:) = CS%H_to_m * h(:,:,:)

    sosga = global_area_mean(sfc_state%SSS, G)

    call generic_tracer_coupler_set(sfc_state%tr_fields,&
         ST=sfc_state%SST,&
         SS=sfc_state%SSS,&
         rho=rho0,& !nnz: required for MOM5 and previous versions.
         ilb=G%isd, jlb=G%jsd,&
         dzt=dzt,& !This is needed for the Mocsy method of carbonate system vars
         tau=1,sosga=sosga,model_time=get_diag_time_end(CS%diag))

    !Output diagnostics via diag_manager for all tracers in this module
!    if (.NOT. associated(CS%g_tracer_list)) call MOM_error(FATAL, trim(sub_name)//&
!         "No tracer in the list.")
!    call g_tracer_send_diag(CS%g_tracer_list, get_diag_time_end(CS%diag), tau=1)
    !Niki: The problem with calling diagnostic outputs here is that this subroutine is called every dt_cpld
    !      hence if dt_therm > dt_cpld we get output (and contribution to the mean) at times that tracers
    !      had not been updated.
    !      Moving this to the end of column physics subrotuine fixes this issue.

  end subroutine MOM_generic_tracer_surface_state

!ALL PE subroutine on Ocean!  Due to otpm design the fluxes should be initialized like this on ALL PE's!
  subroutine MOM_generic_flux_init(verbosity)
    integer, optional, intent(in) :: verbosity  !< A 0-9 integer indicating a level of verbosity.

    integer :: ind
    character(len=fm_string_len)   :: g_tracer_name,longname, package,units,old_package,file_in,file_out
    real :: const_init_value
    character(len=128), parameter :: sub_name = 'MOM_generic_flux_init'
    type(g_tracer_type), pointer :: g_tracer_list,g_tracer,g_tracer_next

    if (.not. g_registered) then
      call generic_tracer_register()
      g_registered = .true.
    endif

    call generic_tracer_get_list(g_tracer_list)
    if (.NOT. associated(g_tracer_list)) then
      call MOM_error(WARNING, trim(sub_name)// ": No generic tracer in the list.")
      return
    endif

    g_tracer=>g_tracer_list
    do

      call g_tracer_flux_init(g_tracer) !, verbosity=verbosity) !### Add this after ocean shared is updated.

      ! traverse the linked list till hit NULL
      call g_tracer_get_next(g_tracer, g_tracer_next)
      if (.NOT. associated(g_tracer_next)) exit
      g_tracer=>g_tracer_next

    enddo

  end subroutine MOM_generic_flux_init

  subroutine MOM_generic_tracer_fluxes_accumulate(flux_tmp, weight)
    type(forcing), intent(in)    :: flux_tmp  !< A structure containing pointers to
                                              !! thermodynamic and tracer forcing fields.
    real,          intent(in)    :: weight    !< A weight for accumulating this flux

    call generic_tracer_coupler_accumulate(flux_tmp%tr_fluxes, weight)

  end subroutine MOM_generic_tracer_fluxes_accumulate

  !> Copy the requested tracer into an array.
  subroutine MOM_generic_tracer_get(name,member,array, CS)
    character(len=*),         intent(in)  :: name   !< Name of requested tracer.
    character(len=*),         intent(in)  :: member !< The tracer element to return.
    real, dimension(:,:,:),   intent(out) :: array  !< Array filled by this routine.
    type(MOM_generic_tracer_CS), pointer :: CS   !< Pointer to the control structure for this module.

    real, dimension(:,:,:),   pointer :: array_ptr
    character(len=128), parameter :: sub_name = 'MOM_generic_tracer_get'

    call g_tracer_get_pointer(CS%g_tracer_list,name,member,array_ptr)
    array(:,:,:) = array_ptr(:,:,:)

  end subroutine MOM_generic_tracer_get

  !> This subroutine deallocates the memory owned by this module.
  subroutine end_MOM_generic_tracer(CS)
    type(MOM_generic_tracer_CS), pointer :: CS   !< Pointer to the control structure for this module.

    call generic_tracer_end()

    if (associated(CS)) then
      deallocate(CS)
    endif
  end subroutine end_MOM_generic_tracer

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

end module MOM_generic_tracer
