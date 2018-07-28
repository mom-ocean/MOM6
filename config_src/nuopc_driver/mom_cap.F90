!>
!! @mainpage MOM NUOPC Cap
!! @author Fei Liu (fei.liu@gmail.com)
!! @date 5/10/13 Original documentation
!! @author Rocky Dunlap (rocky.dunlap@noaa.gov)
!! @date 1/12/17 Moved to doxygen
!!
!! @tableofcontents
!!
!! @section Overview Overview
!!
!! **This MOM cap has been tested with MOM6.**
!!
!! This document describes the MOM "cap", which is a small software layer that is
!! required when the [MOM ocean model] (http://mom-ocean.org/web)
!! is used in [National Unified Operation Prediction Capability]
!! (http://www.earthsystemcog.org/projects/nuopc) (NUOPC) coupled systems.
!! The NUOPC Layer is a software layer built on top of the [Earth System Modeling
!! Framework] (https://www.earthsystemcog.org/projects/esmf) (ESMF).
!! ESMF is a high-performance modeling framework that provides
!! data structures, interfaces, and operations suited for building coupled models
!! from a set of components. NUOPC refines the capabilities of ESMF by providing
!! a more precise definition of what it means for a model to be a component and
!! how components should interact and share data in a coupled system. The NUOPC
!! Layer software is designed to work with typical high-performance models in the
!! Earth sciences domain, most of which are written in Fortran and are based on a
!! distributed memory model of parallelism (MPI).
!! A NUOPC "cap" is a Fortran module that serves as the interface to a model
!! when it's used in a NUOPC-based coupled system.
!! The term "cap" is used because it is a small software layer that sits on top
!! of model code, making calls into it and exposing model data structures in a
!! standard way. For more information about creating NUOPC caps in general, please
!! see the [Building a NUOPC Model]
!! (http://www.earthsystemmodeling.org/esmf_releases/non_public/ESMF_7_0_0/NUOPC_howtodoc/)
!! how-to document.
!!
!! The MOM cap package includes the cap itself (mom_cap.F90, a Fortran module), a
!! set of time utilities (time_utils.F90) for converting between ESMF and FMS
!! time types, and two makefiles. Also included are self-describing dependency
!! makefile fragments (mom.mk and mom.mk.template), although these can be generated
!! by the makefiles for specific installations of the MOM cap.
!!
!! @subsection CapSubroutines Cap Subroutines
!!
!! The MOM cap Fortran module contains a set of subroutines that are required
!! by NUOPC.  These subroutines are called by the NUOPC infrastructure according
!! to a predefined calling sequence.  Some subroutines are called during
!! initialization of the coupled system, some during the run of the coupled
!! system, and some during finalization of the coupled system.  The initialization
!! sequence is the most complex and is governed by the NUOPC technical rules.
!! Details about the initialization sequence can be found in the [NUOPC Reference Manual]
!! (http://www.earthsystemmodeling.org/esmf_releases/non_public/ESMF_7_0_0/NUOPC_refdoc/node3.html
!!  #SECTION00034000000000000000).
!!
!! A particularly important part of the NUOPC intialization sequence is to establish
!! field connections between models.  Simply put, a field connection is established
!! when a field output by one model can be consumed by another.  As an example, the
!! MOM model is able to accept a precipitation rate when coupled to an atmosphere
!! model.  In this case a field connection will be established between the precipitation
!! rate exported from the atmosphere and the precipitation rate imported into the
!! MOM model.  Because models may uses different variable names for physical
!! quantities, NUOPC relies on a set of standard names and a built-in, extensible
!! standard name dictionary to match fields between models.  More information about
!! the use of standard names can be found in the [NUOPC Reference Manual]
!! (http://www.earthsystemmodeling.org/esmf_releases/non_public/ESMF_7_0_0/NUOPC_refdoc/node3.html
!!  #SECTION00032000000000000000).
!!
!! Two key initialization phases that appear in every NUOPC cap, including this MOM
!! cap are the field "advertise" and field "realize" phases.  *Advertise* is a special
!! NUOPC term that refers to a model participating in a coupled system
!! providing a list of standard names of required import fields and available export
!! fields.  In other words, each model will advertise to the other models which physical fields
!! it needs and which fields it can provide when coupled. NUOPC compares all of the advertised
!! standard names and creates a set of unidirectional links, each from one export field
!! in a model to one import field in another model.  When these connections have been established,
!! all models in the coupled system need to provide a description of their geographic
!! grid (e.g., lat-lon, tri-polar, cubed sphere, etc.) and allocate their connected
!! fields on that grid.  In NUOPC terms, this is refered to as *realizing* a set of
!! fields.  NUOPC relies on ESMF data types for this, such as the [ESMF_Grid]
!! (http://www.earthsystemmodeling.org/esmf_releases/public/last/ESMF_refdoc/node5.html#SECTION05080000000000000000)
!! type, which describes logically rectangular grids and the [ESMF_Field]
!! (http://www.earthsystemmodeling.org/esmf_releases/public/last/ESMF_refdoc/node5.html#SECTION05030000000000000000)
!! type, which wraps a models data arrays and provides basic metadata. Because ESMF supports
!! interpolation between different grids (sometimes called "regridding" or "grid remapping"),
!! it is not necessary that models share a grid.  As you will see below
!! the *advertise* and *realize* phases each have a subroutine in the HYCOM cap.
!!
!! The following table summarizes the NUOPC-required subroutines that appear in the
!! MOM cap.  The "Phase" column says whether the subroutine is called during the
!! initialization, run, or finalize part of the coupled system run.
!!
!! Phase    | MOM Cap Subroutine                                                 |  Description
!! ---------|--------------------------------------------------------------------|--------------------------------------
!! Init     | [InitializeP0] (@ref mom_cap_mod::initializep0)                    | Sets the Initialize Phase Definition
!!                                                                               |  (IPD) version to use
!! Init     | [InitializeAdvertise] (@ref mom_cap_mod::initializeadvertise)      | Advertises standard names of import
!!                                                                               |  and export fields
!! Init     | [InitializeRealize] (@ref mom_cap_mod::initializerealize)          | Creates an ESMF_Grid for the MOM grid
!!                                                                               |  as well as ESMF_Fields for import
!!                                                                               |  and export fields
!! Run      | [ModelAdvance] (@ref mom_cap_mod::modeladvance)                    | Advances the model by a timestep
!! Final    | [Finalize] (@ref mom_cap_mod::ocean_model_finalize)                | Cleans up
!!
!! @section UnderlyingModelInterfaces Underlying Model Interfaces
!!
!!
!! @subsection DomainCreation Domain Creation
!!
!! The MOM tripolar grid is represented as a 2D `ESMF_Grid` and coupling fields are placed
!! on this grid. Calls related to creating the grid are located in the [InitializeRealize]
!! (@ref mom_cap_mod::initializerealize) subroutine, which is called by the NUOPC infrastructure
!! during the intialization sequence.
!!
!! The cap determines parameters for setting up the grid by calling subroutines in the
!! `mpp_domains_mod` module. The global domain size is determined by calling `mpp_get_global_domain()`.
!! A check is in place to ensure that there is only a single tile in the domain (the
!! cap is currently limited to one tile; multi-tile mosaics are not supported).  The
!! decomposition across processors is determined via calls to `mpp_get_compute_domains()`
!! (to retrieve decomposition block indices) and `mpp_get_pelist()` (to determine how
!! blocks are assigned to processors).
!!
!! The grid is created in several steps:
!!  - an `ESMF_DELayout` is created based on the pelist from MOM
!!  - an `ESMF_DistGrid` is created over the global index space. Connections are set
!!    up so that the index space is periodic in the first dimension and has a
!!    fold at the top for the bipole. The decompostion blocks are also passed in
!!    along with the `ESMF_DELayout` mentioned above.
!!  - an `ESMF_Grid` is then created by passing in the above `ESMF_DistGrid`.
!!
!! Masks, areas, center (tlat, tlon), and corner (ulat, ulon) coordinates are then added to the `ESMF_Grid`
!! by retrieving those fields from MOM with calls to `ocean_model_data_get()`.
!!
!! @subsection Initialization Initialization
!!
!! During the [InitializeAdvertise] (@ref mom_cap_mod::initializeadvertise) phase, calls are
!! made to MOM's native initialization subroutines, including `fms_init()`, `constants_init()`,
!! `field_manager_init()`, `diag_manager_init()`, and `set_calendar_type()`.  The MPI communicator
!! is pulled in through the ESMF VM object for the MOM component. The dt and start time are set
!! from parameters from the incoming ESMF clock with calls to `set_time()` and `set_date().`
!!
!!
!! @subsection Run Run
!!
!! The [ModelAdvance] (@ref mom_cap_mod::modeladvance) subroutine is called by the NUOPC
!! infrastructure when it's time for MOM to advance in time. During this subroutine, there is a
!! call into the MOM update routine:
!!
!!      call update_ocean_model(Ice_ocean_boundary, Ocean_state, Ocean_public, Time, Time_step_coupled)
!!
!! Prior to this call, the cap performs a few steps:
!! - the `Time` and `Time_step_coupled` parameters, based on FMS types, are derived from the incoming ESMF clock
!! - there are calls to two stubs: `ice_ocn_bnd_from_data()` and `external_coupler_sbc_before()` - these are currently
!!   inactive, but may be modified to read in import data from file or from an external coupler
!! - diagnostics are optionally written to files `field_ocn_import_*`, one for each import field
!! - import fields are prepared:
!!    - the sign is reversed on `mean_evap_rate` and `mean_sensi_heat_flux`
!!    - momentum flux vectors are rotated to internal grid
!! - optionally, a call is made to `ocean_model_restart()` at the interval `restart_interval`
!!
!! After the call to `update_ocean_model()`, the cap performs these steps:
!! - the `ocean_mask` export is set to match that of the internal MOM mask
!! - the `freezing_melting_potential` export is converted from J m-2 to W m-2 by dividing by the coupling interval
!! - vector rotations are applied to the `ocean_current_zonal` and `ocean_current_merid` exports, back to lat-lon grid
!! - diagnostics are optionally written to files `field_ocn_export_*`, one for each export field
!! - a call is made to `external_coupler_sbc_after()` to update exports from an external coupler (currently an inactive
!!    stub)
!!
!! @subsubsection VectorRotations Vector Rotations
!!
!! Vector rotations are applied to incoming momentum fluxes (from regular lat-lon to tripolar grid) and
!! outgoing ocean currents (from tripolar to regular lat-lon). The rotation angles are provided
!! from the native MOM grid by a call to `get_ocean_grid(Ocean_grid)`.
!! The cosine and sine of the rotation angle are:
!!
!!     ocean_grid%cos_rot(i,j)
!!     ocean_grid%sin_rot(i,j)
!!
!! The rotation of momentum flux from regular lat-lon to tripolar is:
!! \f[
!! \begin{bmatrix}
!! \tau_x' \\
!! \tau_y'
!! \end{bmatrix} =
!! \begin{bmatrix}
!!  cos \theta   & sin \theta \\
!!  -sin \theta  & cos \theta
!! \end{bmatrix} *
!! \begin{bmatrix}
!! \tau_x \\
!! \tau_y
!! \end{bmatrix}
!! \f]
!!
!! The rotation of ocean current from tripolar to regular lat-lon is:
!! \f[
!! \begin{bmatrix}
!! u' \\
!! v'
!! \end{bmatrix} =
!! \begin{bmatrix}
!!  cos \theta   & -sin \theta \\
!!  sin \theta  & cos \theta
!! \end{bmatrix} *
!! \begin{bmatrix}
!! u \\
!! v
!! \end{bmatrix}
!! \f]
!! @subsection Finalization Finalization
!!
!! NUOPC infrastructure calls [ocean_model_finalize] (@ref mom_cap_mod::ocean_model_finalize)
!! at the end of the run. This subroutine is a hook to call into MOM's native shutdown
!! procedures:
!!
!!     call ocean_model_end (ocean_public, ocean_State, Time)
!!     call diag_manager_end(Time )
!!     call field_manager_end
!!     call fms_io_exit
!!     call fms_end
!!
!! @section ModelFields Model Fields
!!
!! The following tables list the import and export fields currently set up in the MOM cap.
!!
!! @subsection ImportFields Import Fields
!!
!! Standard Name             | Units      | Model Variable  | Description                      | Notes
!! --------------------------|------------|-----------------|---------------------------------------|-------------------
!! inst_pres_height_surface  | Pa         | p               | pressure of overlying sea ice and atmosphere
!! mass_of_overlying_sea_ice | kg         | mi              | mass of overlying sea ice        | |
!! mean_calving_heat_flx     | W m-2      | calving_hflx    | heat flux, relative to 0C, of frozen land water into ocean
!! mean_calving_rate         | kg m-2 s-1 | calving         | mass flux of frozen runoff       | |
!! mean_evap_rate            | kg m-2 s-1 | q_flux          | specific humidity flux           | sign reversed (- evap)
!! mean_fprec_rate           | kg m-2 s-1 | fprec           | mass flux of frozen precip       | |
!! mean_merid_moment_flx     | Pa         | v_flux          | j-directed wind stress into ocean
!!                                              | [vector rotation] (@ref VectorRotations) applied - lat-lon to tripolar
!! mean_net_lw_flx           | W m-2      | lw_flux         | long wave radiation              | |
!! mean_net_sw_ir_dif_flx    | W m-2      | sw_flux_nir_dif | diffuse near IR shortwave radiation| |
!! mean_net_sw_ir_dir_flx    | W m-2      | sw_flux_nir_dir | direct near IR shortwave radiation| |
!! mean_net_sw_vis_dif_flx   | W m-2      | sw_flux_vis_dif | diffuse visible shortware radiation| |
!! mean_net_sw_vis_dir_flx   | W m-2      | sw_flux_vis_dir | direct visible shortware radiation| |
!! mean_prec_rate            | kg m-2 s-1 | lprec           | mass flux of liquid precip       | |
!! mean_runoff_heat_flx      | W m-2      | runoff_hflx     | heat flux, relative to 0C, of liquid land water into ocean
!! mean_runoff_rate          | kg m-2 s-1 | runoff          | mass flux of liquid runoff       | |
!! mean_salt_rate            | kg m-2 s-1 | salt_flux       | salt flux                        | |
!! mean_sensi_heat_flx       | W m-2      | t_flux          | sensible heat flux into ocean    | sign reversed (- sensi)
!! mean_zonal_moment_flx     | Pa         | u_flux          | j-directed wind stress into ocean
!!                                              | [vector rotation] (@ref VectorRotations) applied - lat-lon to tripolar
!!
!!
!! @subsection ExportField Export Fields
!!
!! Export fields are populated from the `ocean_public` parameter (type `ocean_public_type`)
!! after the call to `update_ocean_model()`.
!!
!! Standard Name              | Units | Model Variable | Description                               | Notes
!! ---------------------------|-------|----------------|-------------------------------------------|--------------------
!! freezing_melting_potential | W m-2 | frazil         | accumulated heating from frazil formation
!!                                              | cap converts model units (J m-2) to (W m-2) for export
!! ocean_mask                 |       |                | ocean mask                                | |
!! ocn_current_merid          | m s-1 | v_surf         | j-directed surface velocity on u-cell
!!                                              | [vector rotation] (@ref VectorRotations) applied - tripolar to lat-lon
!! ocn_current_zonal          | m s-1 | u_surf         | i-directed surface velocity on u-cell
!!                                              | [vector rotation] (@ref VectorRotations) applied - tripolar to lat-lon
!! s_surf                     | psu   | s_surf         | sea surface salinity on t-cell            | |
!! sea_lev                    | m     | sea_lev        | sea level
!!                                              | model computation is eta_t + patm/(rho0*grav) - eta_geoid - eta_tide
!! sea_surface_temperature    | K     | t_surf         | sea surface temperature on t-cell         | |
!!
!! @subsection MemoryManagement Memory Management
!!
!! The MOM cap has an internal state type with pointers to three
!! types defined by MOM. There is also a small wrapper derived type
!! required to associate an internal state instance
!! with the ESMF/NUOPC component:
!!
!!     type ocean_internalstate_type
!!        type(ocean_public_type),       pointer :: ocean_public_type_ptr
!!        type(ocean_state_type),        pointer :: ocean_state_type_ptr
!!        type(ice_ocean_boundary_type), pointer :: ice_ocean_boundary_type_ptr
!!     end type
!!
!!     type ocean_internalstate_wrapper
!!        type(ocean_internalstate_type), pointer :: ptr
!!     end type
!!
!! The member of type `ocean_public_type` stores ocean surface fields used during the coupling.
!! The member of type `ocean_state_type` is required by the ocean driver,
!! although its internals are private (not to be used by the coupling directly).
!! This type is passed to the ocean init and update routines
!! so that it can maintain state there if desired.
!! The member of type `ice_ocean_boundary_type` is populated by this cap
!! with incoming coupling fields from other components. These three derived types are allocated during the
!! [InitializeAdvertise] (@ref mom_cap_mod::initializeadvertise) phase.  Also during that
!! phase, the `ice_ocean_boundary` type members are all allocated using bounds retrieved
!! from `mpp_get_compute_domain()`.
!!
!! During the [InitializeRealize] (@ref mom_cap_mod::initializerealize) phase,
!! `ESMF_Field`s are created for each of the coupling fields in the `ice_ocean_boundary`
!! and `ocean_public_type` members of the internal state. These fields directly reference into the members of
!! the `ice_ocean_boundary` and `ocean_public_type` so that memory-to-memory copies are not required to move
!! data from the cap's import and export states to the memory areas used internally
!! by MOM.
!!
!! @subsection IO I/O
!!
!! The cap can optionally output coupling fields for diagnostic purposes if the ESMF attribute
!! "DumpFields" has been set to "true". In this case the cap will write out NetCDF files
!! with names "field_ocn_import_<fieldname>.nc" and "field_ocn_export_<fieldname>.nc".
!! Additionally, calls will be made to the cap subroutine [dumpMomInternal]
!! (@ref mom_cap_mod::dumpmominternal) to write out model internal fields to files
!! named "field_ocn_internal_<fieldname>.nc".  In all cases these NetCDF files will
!! contain a time series of field data.
!!
!! @section BuildingAndInstalling Building and Installing
!!
!! There are two makefiles included with the MOM cap, makefile and makefile.nuopc.
!! The makefile.nuopc file is intended to be used within another build system, such
!! as the NEMSAppBuilder.  The regular makefile can be used generally for building
!! and installing the cap.  Two variables must be customized at the top:
!! - `INSTALLDIR` - where to copy the cap library and dependent libraries
!! - `NEMSMOMDIR` - location of the MOM library and FMS library
!!
!! To install run:
!!    $ make install
!!
!! A makefile fragment, mom.mk, will also be copied into the directory. The fragment
!! defines several variables that can be used by another build system to include the
!! MOM cap and its dependencies.
!!
!! @subsection Dependencies Dependencies
!!
!! The MOM cap is dependent on the MOM library itself (lib_ocean.a) and the FMS
!! library (lib_FMS.a).
!!
!! @section RuntimeConfiguration Runtime Configuration
!!
!! At runtime, the MOM cap can be configured with several options provided
!! as ESMF attributes.  Attributes can be set in the cap by the NUOPC Driver
!! above this cap, or in some systems (e.g., NEMS) attributes are set by
!! reading in from a configuration file.  The available attributes are:
!!
!! * `DumpFields` - when set to "true", write out diagnostic NetCDF files for import/export/internal fields
!! * `ProfileMemory` - when set to "true", write out memory usage information to the ESMF log files; this
!!   information is written when entering and leaving the [ModelAdvance]
!!   (@ref mom_cap_mod::modeladvance) subroutine and before and after the call to
!!   `update_ocean_model()`.
!! * `restart_interval` - integer number of seconds indicating the interval at
!!   which to call `ocean_model_restart()`; no restarts written if set to 0
!! * `GridAttachArea` - when set to "true", this option indicates that MOM grid attaches cell area
!!   using internal values computed in MOM. The default value is "false", grid cell area will
!!   be computed in ESMF.
!!
!!
!! @section Repository
!! The MOM NUOPC cap is maintained in a GitHub repository:
!! https://github.com/feiliuesmf/nems_mom_cap
!!
!! @section References
!!
!! - [MOM Home Page] (http://mom-ocean.org/web)
!!
!!
module mom_cap_mod
  use constants_mod,            only: constants_init
  use data_override_mod,        only: data_override_init, data_override
  use diag_manager_mod,         only: diag_manager_init, diag_manager_end
  use field_manager_mod,        only: field_manager_init, field_manager_end
  use fms_mod,                  only: fms_init, fms_end, open_namelist_file, check_nml_error
  use fms_mod,                  only: close_file, file_exist, uppercase
  use fms_io_mod,               only: fms_io_exit
  use mpp_domains_mod,          only: domain2d, mpp_get_compute_domain, mpp_get_compute_domains
  use mpp_domains_mod,          only: mpp_get_ntile_count, mpp_get_pelist, mpp_get_global_domain
  use mpp_domains_mod,          only: mpp_get_domain_npes, mpp_global_field
  use mpp_io_mod,               only: mpp_open, MPP_RDONLY, MPP_ASCII, MPP_OVERWR, MPP_APPEND, mpp_close, MPP_SINGLE
  use mpp_mod,                  only: input_nml_file, mpp_error, FATAL, NOTE, mpp_pe, mpp_npes, mpp_set_current_pelist
  use mpp_mod,                  only: stdlog, stdout, mpp_root_pe, mpp_clock_id
  use mpp_mod,                  only: mpp_clock_begin, mpp_clock_end, MPP_CLOCK_SYNC
  use mpp_mod,                  only: MPP_CLOCK_DETAILED, CLOCK_COMPONENT, MAXPES
  use time_interp_external_mod, only: time_interp_external_init
  use time_manager_mod,         only: set_calendar_type, time_type, increment_date
  use time_manager_mod,         only: set_time, set_date, get_time, get_date, month_name
  use time_manager_mod,         only: GREGORIAN, JULIAN, NOLEAP, THIRTY_DAY_MONTHS, NO_CALENDAR
  use time_manager_mod,         only: operator( <= ), operator( < ), operator( >= )
  use time_manager_mod,         only: operator( + ),  operator( - ), operator( / )
  use time_manager_mod,         only: operator( * ), operator( /= ), operator( > )
  use time_manager_mod,         only: date_to_string
  use time_manager_mod,         only: fms_get_calendar_type => get_calendar_type
  use MOM_domains,              only: MOM_infra_init, num_pes, root_pe, pe_here
  use MOM_surface_forcing,      only: IOB_allocate
  use MOM_file_parser,          only: get_param, log_version, param_file_type, close_param_file
  use MOM_get_input,            only: Get_MOM_Input, directories
  use MOM_domains,              only: pass_var
  use MOM_error_handler,        only: is_root_pe
  use MOM_ocean_model,          only: ice_ocean_boundary_type
  use MOM_grid,                 only: ocean_grid_type
  use MOM_ocean_model,          only: ocean_model_restart, ocean_public_type, ocean_state_type
  use MOM_ocean_model,          only: ocean_model_data_get, ocean_model_init_sfc
  use MOM_ocean_model,          only: ocean_model_init, update_ocean_model, ocean_model_end, get_ocean_grid
#ifdef CESMCOUPLED
  use mom_cap_methods,          only: mom_import, mom_export
  use esmFlds,                  only: flds_scalar_name, flds_scalar_num
  use esmFlds,                  only: flds_scalar_index_nx, flds_scalar_index_ny
  use shr_file_mod,             only: shr_file_getUnit, shr_file_freeUnit
  use shr_file_mod,             only: shr_file_getLogUnit, shr_file_getLogLevel
  use shr_file_mod,             only: shr_file_setLogUnit, shr_file_setLogLevel
  use shr_nuopc_time_mod,       only: shr_nuopc_time_alarmInit
#endif

  use ESMF                      ! TODO: only: ...
  use NUOPC                     ! TODO: only: ...
  use NUOPC_Model, &            ! TODO: only: ...
    model_routine_SS           => SetServices, &
    model_label_DataInitialize => label_DataInitialize, &
    model_label_Advance        => label_Advance, &
#ifdef CESMCOUPLED
    model_label_SetRunClock    => label_SetRunClock, &
#endif
    model_label_Finalize       => label_Finalize

  use time_utils_mod,           only: esmf2fms_time

  implicit none
  private

  public SetServices

  type ocean_internalstate_type
    type(ocean_public_type),       pointer :: ocean_public_type_ptr
    type(ocean_state_type),        pointer :: ocean_state_type_ptr
    type(ice_ocean_boundary_type), pointer :: ice_ocean_boundary_type_ptr
  end type

  type ocean_internalstate_wrapper
    type(ocean_internalstate_type), pointer :: ptr
  end type

  type fld_list_type
    character(len=64) :: stdname
    character(len=64) :: shortname
    character(len=64) :: transferOffer
    logical           :: assoc    ! is the farrayPtr associated with internal data
    real(ESMF_KIND_R8), dimension(:,:), pointer :: farrayPtr
  end type fld_list_type

  integer,parameter    :: fldsMax = 100
  integer              :: fldsToOcn_num = 0
  type (fld_list_type) :: fldsToOcn(fldsMax)
  integer              :: fldsFrOcn_num = 0
  type (fld_list_type) :: fldsFrOcn(fldsMax)

  integer              :: debug = 0
  integer              :: import_slice = 1
  integer              :: export_slice = 1
  character(len=256)   :: tmpstr
  integer              :: dbrc
  type(ESMF_Grid)      :: mom_grid_i

#ifdef CESMCOUPLED
  logical                 :: write_diagnostics = .false.
  integer                 :: logunit  ! stdout logging unit number
  character(len=32)       :: runtype  ! run type
#else
  logical                 :: write_diagnostics = .true.
#endif
  logical                 :: profile_memory = .true.
  logical                 :: grid_attach_area = .false.
  integer(ESMF_KIND_I8)   :: restart_interval
  logical                 :: sw_decomp
  character(len=*),parameter :: u_file_u = &
       __FILE__

contains

  !===============================================================================
  !> NUOPC SetService method is the only public entry point.
  !! SetServices registers all of the user-provided subroutines
  !! in the module with the NUOPC layer.
  !!
  !! @param gcomp an ESMF_GridComp object
  !! @param rc return code
  subroutine SetServices(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname='(mom_cap:SetServices)'

    rc = ESMF_SUCCESS

    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !------------------
    ! attach specializing method(s)
    !------------------

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_DataInitialize, &
      specRoutine=DataInitialize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

#ifdef CESMCOUPLED
    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
      specRoutine=ocean_model_finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine SetServices

  !===============================================================================

  !> First initialize subroutine called by NUOPC.  The purpose
  !! is to set which version of the Initialize Phase Definition (IPD)
  !! to use.
  !!
  !! For this MOM cap, we are using IPDv01.
  !!
  !! @param gcomp an ESMF_GridComp object
  !! @param importState an ESMF_State object for import fields
  !! @param exportState an ESMF_State object for export fields
  !! @param clock an ESMF_Clock object
  !! @param rc return code
  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc

    character(len=10)                         :: value
    character(len=*),parameter  :: subname='(mom_cap:InitializeP0)'

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_AttributeGet(gcomp, name="DumpFields", value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    !    write_diagnostics=(trim(value)=="true")
    call ESMF_LogWrite('MOM_CAP:DumpFields = '//trim(value), ESMF_LOGMSG_INFO, rc=dbrc)

    call ESMF_AttributeGet(gcomp, name="ProfileMemory", value=value, defaultValue="true", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    profile_memory=(trim(value)/="false")
    call ESMF_LogWrite('MOM_CAP:ProfileMemory = '//trim(value), ESMF_LOGMSG_INFO, rc=dbrc)

#ifndef CESMCOUPLED
    ! Retrieve restart_interval in (seconds)
    ! A restart_interval value of 0 means no restart will be written.
    call ESMF_AttributeGet(gcomp, name="restart_interval", value=value, defaultValue="0", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    restart_interval = ESMF_UtilString2Int(value, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if(restart_interval < 0) then
      call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
        msg="MOM_CAP: OCN attribute: restart_interval cannot be negative.", &
        line=__LINE__, &
        file=__FILE__, rcToReturn=rc)
      return
    endif
    call ESMF_LogWrite('MOM_CAP:restart_interval = '//trim(value), ESMF_LOGMSG_INFO, rc=dbrc)
#endif

    call ESMF_AttributeGet(gcomp, name="GridAttachArea", value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    grid_attach_area=(trim(value)=="true")
    call ESMF_LogWrite('MOM_CAP:GridAttachArea = '//trim(value), ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine

  !===============================================================================

  !> Called by NUOPC to advertise import and export fields.  "Advertise"
  !! simply means that the standard names of all import and export
  !! fields are supplied.  The NUOPC layer uses these to match fields
  !! between components in the coupled system.
  !!
  !! @param gcomp an ESMF_GridComp object
  !! @param importState an ESMF_State object for import fields
  !! @param exportState an ESMF_State object for export fields
  !! @param clock an ESMF_Clock object
  !! @param rc return code
  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)                    :: gcomp
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Clock)                       :: clock
    integer, intent(out)                   :: rc

    type(ESMF_VM)                          :: vm
    type(ESMF_Time)                        :: MyTime
    type(ESMF_TimeInterval)                :: TINT
    type (ocean_public_type),      pointer :: ocean_public => NULL()
    type (ocean_state_type),       pointer :: ocean_state => NULL()
    type(ice_ocean_boundary_type), pointer :: Ice_ocean_boundary => NULL()
    type(ocean_internalstate_wrapper)      :: ocean_internalstate
    type(ocean_grid_type),         pointer :: ocean_grid => NULL()
    type(time_type)                        :: Run_len      ! length of experiment
    type(time_type)                        :: Time
    type(time_type)                        :: Time_restart
    type(time_type)                        :: DT
    integer                                :: DT_OCEAN
    integer                                :: isc,iec,jsc,jec
    integer                                :: dt_cpld  = 86400
    integer                                :: year=0, month=0, day=0, hour=0, minute=0, second=0
    integer                                :: mpi_comm_mom
    integer                                :: i,n
    character(80)                          :: stdname, shortname
#ifdef CESMCOUPLED
    integer                                :: nflds
    character(len=32)                      :: starttype            ! model start type
    character(len=512)                     :: diro
    character(len=512)                     :: logfile
    character(len=64)                      :: cvalue
    integer                                :: shrlogunit           ! original log unit
    integer                                :: shrloglev            ! original log level
    integer                                :: inst_index           ! number of current instance (ie. 1)
    character(len=16)                      :: inst_name            ! fullname of current instance (ie. "lnd_0001")
    character(len=16)                      :: inst_suffix = ""     ! char string associated with instance
                                                                   ! (ie. "_0001" or "")
    logical                                :: isPresent
    character(len=384)                     :: restart_pointer_file ! File name for restart pointer file
    character(len=384)                     :: restartfile          ! Path/Name of restart file
    character(len=384)                     :: restartname          ! The restart file name (no dir)
    integer                                :: nu                   ! i/o unit to read pointer file
#endif
    character(len=*),parameter             :: subname='(mom_cap:InitializeAdvertise)'
    !--------------------------------

    rc = ESMF_SUCCESS

    allocate(Ice_ocean_boundary)
    !allocate(ocean_state) ! ocean_model_init allocate this pointer
    allocate(ocean_public)
    allocate(ocean_internalstate%ptr)
    ocean_internalstate%ptr%ice_ocean_boundary_type_ptr => Ice_ocean_boundary
    ocean_internalstate%ptr%ocean_public_type_ptr       => ocean_public
    ocean_internalstate%ptr%ocean_state_type_ptr        => ocean_state

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(VM, mpiCommunicator=mpi_comm_mom, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(CLOCK, currTIME=MyTime, TimeStep=TINT,  RC=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimeGet (MyTime, YY=YEAR, MM=MONTH, DD=DAY, H=HOUR, M=MINUTE, S=SECOND, RC=rc )
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    CALL ESMF_TimeIntervalGet(TINT, S=DT_OCEAN, RC=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call fms_init(mpi_comm_mom)
    call constants_init
    call field_manager_init
    call set_calendar_type (JULIAN)
    call diag_manager_init

    ! this ocean connector will be driven at set interval
    dt_cpld = DT_OCEAN
    DT = set_time (DT_OCEAN, 0)
    Time = set_date (YEAR,MONTH,DAY,HOUR,MINUTE,SECOND)

#ifdef CESMCOUPLED

    ! determine instance information
    call NUOPC_CompAttributeGet(gcomp, name="inst_name", value=inst_name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_CompAttributeGet(gcomp, name="inst_index", value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    read(cvalue,*) inst_index

    call ESMF_AttributeGet(gcomp, name="inst_suffix", isPresent=isPresent, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (isPresent) then
      call NUOPC_CompAttributeGet(gcomp, name="inst_suffix", value=inst_suffix, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    else
       inst_suffix = ''
    end if

    ! reset shr logging to my log file
    if (is_root_pe()) then
      call NUOPC_CompAttributeGet(gcomp, name="diro", value=diro, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call NUOPC_CompAttributeGet(gcomp, name="logfile", value=logfile, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      logunit = shr_file_getUnit()
      open(logunit,file=trim(diro)//"/"//trim(logfile))
    else
      logunit = 6
    endif

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogLevel(max(shrloglev,1))
    call shr_file_setLogUnit (logunit)

    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    read(cvalue,*) starttype

    if (trim(starttype) == trim('startup')) then
      runtype = "initial"
    else if (trim(starttype) == trim('continue') ) then
      runtype = "continue"
    else if (trim(starttype) == trim('branch')) then
      runtype = "continue"
    else
      call ESMF_LogWrite(subname//' ERROR: unknown starttype '//trim(starttype), ESMF_LOGMSG_ERROR, rc=dbrc)
      rc = ESMF_FAILURE
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end if

    if (runtype == "initial") then

      ! startup (new run) - 'n' is needed below since we don't specify input_filename in input.nml
      ocean_public%is_ocean_pe = .true.
      call ocean_model_init(ocean_public, ocean_state, Time, Time, input_restart_file = 'n')


    else  ! hybrid or branch or continuos runs

      ! read name of restart file in the pointer file
      nu = shr_file_getUnit()
      restart_pointer_file = 'rpointer.ocn'
      if (is_root_pe()) then
        write(logunit,*) 'Reading ocn pointer file: ',restart_pointer_file
      end if
      open(nu, file=restart_pointer_file, form='formatted', status='unknown')
      read(nu,'(a)') restartfile
      close(nu)

      ! initialize from restart file
      if (is_root_pe()) then
        write(logunit,*) 'Reading restart file: ',trim(restartfile)
      end if
      call shr_file_freeUnit(nu)

      ocean_public%is_ocean_pe = .true.
      call ocean_model_init(ocean_public, ocean_state, Time, Time, input_restart_file=trim(restartfile))

    end if
    call ocean_model_init_sfc(ocean_state, ocean_public)

#else

    ocean_public%is_ocean_pe = .true.
    call ocean_model_init(ocean_public, ocean_state, Time, Time)

#endif

    !tcx tcraig This results in errors in CESM with help from Alper
    ! FATAL error "MPP_OPEN: error in OPEN for data_table"
    ! The subroutine data_override_init shouldn't be called because ALLOW_FLUX_ADJUSTMENTS is set to FALSE
    !tcx    call data_override_init(ocean_domain_in = ocean_public%domain)

    call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)
    call IOB_allocate(ice_ocean_boundary, isc, iec, jsc, jec)

    call external_coupler_sbc_init(ocean_public%domain, dt_cpld, Run_len)

    ocean_internalstate%ptr%ocean_state_type_ptr => ocean_state
    call ESMF_GridCompSetInternalState(gcomp, ocean_internalstate, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

#ifdef CESMCOUPLED

    !--------- import fields -------------
    call fld_list_add(fldsToOcn_num, fldsToOcn, trim(flds_scalar_name), "will_provide") ! not in EMC
    call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_rain"     , "will provide") ! -> mean_prec_rat
    call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_snow"     , "will provide") ! -> mean_fprec_rate
    call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_lwdn"     , "will provide")
    call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_swndr"    , "will provide") ! -> mean_net_sw_ir_dif_flx
    call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_swvdr"    , "will provide") ! -> mean_net_sw_vis_dir_flx
    call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_swndf"    , "will provide") ! -> mean_net_sw_ir_dir_flx 
    call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_swvdf"    , "will provide") ! -> mean_net_sw_vis_dif_flx
    call fld_list_add(fldsToOcn_num, fldsToOcn, "Foxx_taux"     , "will provide") ! -> mean_zonal_moment_flx
    call fld_list_add(fldsToOcn_num, fldsToOcn, "Foxx_tauy"     , "will provide") ! -> mean_merid_moment_flx
    call fld_list_add(fldsToOcn_num, fldsToOcn, "Foxx_sen"      , "will provide") ! -> mean_sensi_heat_flx
    call fld_list_add(fldsToOcn_num, fldsToOcn, "Foxx_lat"      , "will provide")
    call fld_list_add(fldsToOcn_num, fldsToOcn, "Foxx_lwup"     , "will provide")
    call fld_list_add(fldsToOcn_num, fldsToOcn, "Foxx_evap"     , "will provide") ! -> mean_evap_rate
    call fld_list_add(fldsToOcn_num, fldsToOcn, "Fioi_salt"     , "will provide") ! -> mean_salt_rate
    call fld_list_add(fldsToOcn_num, fldsToOcn, "Foxx_rofl"     , "will provide")
    call fld_list_add(fldsToOcn_num, fldsToOcn, "Foxx_rofi"     , "will provide")
    call fld_list_add(fldsToOcn_num, fldsToOcn, "Sa_pslv"       , "will provide") ! -> inst_pres_height_surface
    
    ! EMC fields not used
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_runoff_rate"           , "will provide") ! for CESM rofl + rofi
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_lw_flx"            , "will provide") ! for CESM lwup + lwdn
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_calving_rate"          , "will provide") ! not in CESM
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_runoff_heat_flx"       , "will provide") ! not in CESM
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_calving_heat_flx"      , "will provide") ! not in CESM
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "mass_of_overlying_sea_ice"  , "will provide") ! not in CESM

    ! CESM currently not used
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Sw_lamult"     , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Sw_ustokes"    , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Sw_vstokes"    , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Sw_hstokes"    , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Si_ifrac"      , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Fioi_melth"    , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Fioi_meltw"    , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_prec"     , "will provide") 
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_bcphidry" , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_bcphodry" , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_bcphiwet" , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_ocphidry" , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_ocphodry" , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_ocphiwet" , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_dstwet1"  , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_dstwet2"  , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_dstwet3"  , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_dstwet4"  , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_dstdry1"  , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_dstdry2"  , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_dstdry3"  , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Faxa_dstdry4"  , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Fioi_bcphi"    , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Fioi_bcpho"    , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "Fioi_flxdst"   , "will provide")
    ! call fld_list_add(fldsToOcn_num, fldsToOcn, "So_duu10n"     , "will provide")

    ! Optional CESM fields currently not used
    ! call NUOPC_CompAttributeGet(gcomp, name='flds_co2a', value=cvalue, rc=rc)
    ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ! read(cvalue,*) flds_co2a
    ! call ESMF_LogWrite('flds_co2a = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)
    ! call NUOPC_CompAttributeGet(gcomp, name='flds_co2c', value=cvalue, rc=rc)
    ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ! read(cvalue,*) flds_co2c
    ! call ESMF_LogWrite('flds_co2c = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)
    ! if (flds_co2a .or. flds_co2c) then
    !    call fld_list_add(fldsToOcn_num, fldsToOcn, "Sa_co2prog"    , "will provide")
    !    call fld_list_add(fldsToOcn_num, fldsToOcn, "Sa_co2diag"    , "will provide")
    ! end if
    ! call NUOPC_CompAttributeGet(gcomp, name='ice_ncat', value=cvalue, rc=rc)
    ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ! read(cvalue,*) ice_ncat
    ! call ESMF_LogWrite('ice_ncat = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)
    ! call NUOPC_CompAttributeGet(gcomp, name='flds_i2o_per_cat', value=cvalue, rc=rc)
    ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ! read(cvalue,*) flds_i2o_per_cat
    ! call ESMF_LogWrite('flds_i2o_per_cat = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)
    ! if (flds_i2o_per_cat) then
    !    do num = 1, ice_ncat
    !       name = 'Si_ifrac_' // cnum
    !       call fld_list_add(fldsToOcn_num, fldsToOcn, trim(name), "will provide")
    !       name = 'PFioi_swpen_ifrac_' // cnum
    !       call fld_list_add(fldsToOcn_num, fldsToOcn, trim(name), "will provide")
    !    end do
    !    call fld_list_add(fldsToOcn_num, fldsToOcn, "Sf_afrac"      , "will provide")
    !    call fld_list_add(fldsToOcn_num, fldsToOcn, "Sf_afracr"     , "will provide")
    !    call fld_list_add(fldsToOcn_num, fldsToOcn, "Foxx_swnet_afracr", "will provide")
    ! end if
    ! do n = 1,shr_string_listGetNum(ndep_fields)
    !    call shr_string_listGetName(ndep_fields, n, name)
    !    call fld_list_add(fldsToOcn_num, fldsToOcn, trim(name), "will provide")
    ! end do

    !--------- export fields -------------
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, trim(flds_scalar_name), "will_provide") ! not in EMC
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "So_omask"      , "will provide") ! -> ocean_mask
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "So_t"          , "will provide") ! -> sea_surface_temperature
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "So_s"          , "will provide") ! -> s_surf
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "So_u"          , "will provide") ! -> ocn_current_zonal
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "So_v"          , "will provide") ! -> ocn_current_merid
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "So_dhdx"       , "will provide") ! not in EMC
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "So_dhdy"       , "will provide") ! not in EMC
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "So_bldepth"    , "will provide") ! not in EMC 
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "Fioo_q"        , "will provide") ! not in EMC

    ! EMC fields not used
    ! call fld_list_add(fldsFrOcn_num, fldsFrOcn, "sea_lev"                    , "will provide") ! not in CESM
    ! call fld_list_add(fldsFrOcn_num, fldsFrOcn, "freezing_melting_potential" , "will provide") ! not in CESM

    ! Optional CESM fields currently not used
    ! call fld_list_add(fldsFrOcn_num, fldsFrOcn, "So_fswpen"     , "will provide") ! not in EMC
    ! if (flds_co2c) then
    !    call fld_list_add(fldsToOcn_num, fldsFrOcn, "Faoo_fco2_ocn" , "will provide") 
    ! end if


#else

    ! This sets pointers of the fldsToOcn to the iceocean_boundary_type
    ! Don't point directly into mom data YET (last field is optional in interface)
    ! instead, create space for the field when it's "realized".

    !--------- import fields -------------
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_zonal_moment_flx"      , "will provide",&
                      data=Ice_ocean_boundary%u_flux)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_merid_moment_flx"      , "will provide",&
                      data=Ice_ocean_boundary%v_flux)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_sensi_heat_flx"        , "will provide",&
                      data=Ice_ocean_boundary%t_flux)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_evap_rate"             , "will provide",&
                      data=Ice_ocean_boundary%q_flux)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_salt_rate"             , "will provide",&
                      data=Ice_ocean_boundary%salt_flux)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_lw_flx"            , "will provide",&
                      data=Ice_ocean_boundary%lw_flux  )
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_sw_vis_dir_flx"    , "will provide",&
                      data=Ice_ocean_boundary%sw_flux_vis_dir)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_sw_vis_dif_flx"    , "will provide",&
                      data=Ice_ocean_boundary%sw_flux_vis_dif)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_sw_ir_dir_flx"     , "will provide",&
                      data=Ice_ocean_boundary%sw_flux_nir_dir)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_sw_ir_dif_flx"     , "will provide",&
                      data=Ice_ocean_boundary%sw_flux_nir_dif)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_prec_rate"             , "will provide",&
                      data=Ice_ocean_boundary%lprec  )
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_fprec_rate"            , "will provide",&
                      data=Ice_ocean_boundary%fprec  )
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_runoff_rate"           , "will provide",&
                      data=Ice_ocean_boundary%runoff )
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_calving_rate"          , "will provide",&
                      data=Ice_ocean_boundary%calving)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_runoff_heat_flx"       , "will provide",&
                      data=Ice_ocean_boundary%runoff_hflx )
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_calving_heat_flx"      , "will provide",&
                      data=Ice_ocean_boundary%calving_hflx)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "inst_pres_height_surface"   , "will provide",&
                      data=Ice_ocean_boundary%p )
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mass_of_overlying_sea_ice"  , "will provide",&
                      data=Ice_ocean_boundary%mi)

    !--------- export fields -------------
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocean_mask"                 , "will provide")
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "sea_surface_temperature"    , "will provide",&
                      data=ocean_public%t_surf)
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "s_surf"                     , "will provide",&
                      data=ocean_public%s_surf )
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocn_current_zonal"          , "will provide",&
                      data=ocean_public%u_surf )
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocn_current_merid"          , "will provide",&
                      data=ocean_public%v_surf )
   !call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocn_current_idir"           , "will provide")
   !call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocn_current_jdir"           , "will provide")
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "sea_lev"                    , "will provide",&
                      data=ocean_public%sea_lev)
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "freezing_melting_potential" , "will provide",&
                      data=ocean_public%frazil)

#endif

    do n = 1,fldsToOcn_num
      call NUOPC_Advertise(importState, standardName=fldsToOcn(n)%stdname, name=fldsToOcn(n)%shortname, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    enddo

    do n = 1,fldsFrOcn_num
      call NUOPC_Advertise(exportState, standardName=fldsFrOcn(n)%stdname, name=fldsFrOcn(n)%shortname, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    enddo

    ! When running mom6 solo, the rotation angles are not computed internally
    ! in MOM6. We need to calculate cos and sin of rotational angle for MOM6;
    ! the values are stored in ocean_internalstate%ptr%ocean_grid_ptr%cos_rot and sin_rot
    ! The rotation angles are retrieved during run time to rotate incoming
    ! and outgoing vectors
    ! call calculate_rot_angle(ocean_state, ocean_public)
    ! tcraig, this is handled fine internally and if not, then later call this
    ! call initialize_grid_rotation_angle(ocean_grid, PF)

    write(*,*) '----- MOM initialization phase Advertise completed'

  end subroutine InitializeAdvertise

  !===============================================================================

  !> Called by NUOPC to realize import and export fields.  "Realizing" a field
  !! means that its grid has been defined and an ESMF_Field object has been
  !! created and put into the import or export State.
  !!
  !! @param gcomp an ESMF_GridComp object
  !! @param importState an ESMF_State object for import fields
  !! @param exportState an ESMF_State object for export fields
  !! @param clock an ESMF_Clock object
  !! @param rc return code
  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local Variables
    type(ESMF_VM)                              :: vm
    type(ESMF_Grid)                            :: gridIn
    type(ESMF_Grid)                            :: gridOut
    type(ESMF_DeLayout)                        :: delayout
    type(ESMF_Distgrid)                        :: Distgrid
    type(ESMF_DistGridConnection), allocatable :: connectionList(:)
    type (ocean_public_type),      pointer     :: ocean_public   => NULL()
    type (ocean_state_type),       pointer     :: ocean_state => NULL()
    type(ice_ocean_boundary_type), pointer     :: Ice_ocean_boundary => NULL()
    type(ocean_internalstate_wrapper)          :: ocean_internalstate
    integer                                    :: npet, ntiles
    integer                                    :: nxg, nyg, cnt
    integer                                    :: isc,iec,jsc,jec
    integer, allocatable                       :: xb(:),xe(:),yb(:),ye(:),pe(:)
    integer, allocatable                       :: deBlockList(:,:,:)
    integer, allocatable                       :: petMap(:)
    integer, allocatable                       :: deLabelList(:)
    integer, allocatable                       :: indexList(:)
    integer                                    :: ioff, joff
    integer                                    :: i, j, n, i1, j1, n1, icount
    integer                                    :: lbnd1,ubnd1,lbnd2,ubnd2
    integer                                    :: lbnd3,ubnd3,lbnd4,ubnd4
    integer                                    :: nblocks_tot
    logical                                    :: found
    real(ESMF_KIND_R8), allocatable            :: ofld(:,:), gfld(:,:)
    real(ESMF_KIND_R8), pointer                :: t_surf(:,:)
    integer(ESMF_KIND_I4), pointer             :: dataPtr_mask(:,:)
    real(ESMF_KIND_R8), pointer                :: dataPtr_area(:,:)
    real(ESMF_KIND_R8), pointer                :: dataPtr_xcen(:,:)
    real(ESMF_KIND_R8), pointer                :: dataPtr_ycen(:,:)
    real(ESMF_KIND_R8), pointer                :: dataPtr_xcor(:,:)
    real(ESMF_KIND_R8), pointer                :: dataPtr_ycor(:,:)
    type(ESMF_Field)                           :: field_t_surf
    integer                                    :: mpicom
    integer                                    :: localPet
#ifdef CESMCOUPLED
    integer                                    :: shrlogunit       ! original log unit
    integer                                    :: shrloglev        ! original log level
    integer                                    :: inst_index       ! number of current instance (ie. 1)
    character(len=16)                          :: inst_name        ! fullname of current instance (ie. "lnd_0001")
    character(len=16)                          :: inst_suffix = "" ! char string associated with instance
                                                                   ! (ie. "_0001" or "")
    character(len=64)                          :: cvalue
    logical                                    :: isPresent
#endif
    character(len=*),parameter  :: subname='(mom_cap:InitializeRealize)'
    !--------------------------------

    rc = ESMF_SUCCESS

    !----------------------------------------------------------------------------
    ! Get pointers to ocean internal state
    !----------------------------------------------------------------------------

    call ESMF_GridCompGetInternalState(gcomp, ocean_internalstate, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    Ice_ocean_boundary => ocean_internalstate%ptr%ice_ocean_boundary_type_ptr
    ocean_public       => ocean_internalstate%ptr%ocean_public_type_ptr
    ocean_state        => ocean_internalstate%ptr%ocean_state_type_ptr

    !----------------------------------------------------------------------------
    ! Get mpi information
    !----------------------------------------------------------------------------

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm, petCount=npet, mpiCommunicator=mpicom, localPet=localPet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !---------------------------------
    ! global mom grid size
    !---------------------------------

    call mpp_get_global_domain(ocean_public%domain, xsize=nxg, ysize=nyg)
    write(tmpstr,'(a,2i6)') subname//' nxg,nyg = ',nxg,nyg
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    !---------------------------------
    ! number of tiles per PET, assumed to be 1, and number of pes (tiles) total
    !---------------------------------

    ntiles=mpp_get_ntile_count(ocean_public%domain) ! this is tiles on this pe
    if (ntiles /= 1) then
      rc = ESMF_FAILURE
      call ESMF_LogWrite(subname//' ntiles must be 1', ESMF_LOGMSG_ERROR, rc=dbrc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    ntiles=mpp_get_domain_npes(ocean_public%domain)
    write(tmpstr,'(a,1i6)') subname//' ntiles = ',ntiles
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    !---------------------------------
    ! get start and end indices of each tile and their PET
    !---------------------------------

    allocate(xb(ntiles),xe(ntiles),yb(ntiles),ye(ntiles),pe(ntiles))
    call mpp_get_compute_domains(ocean_public%domain, xbegin=xb, xend=xe, ybegin=yb, yend=ye)
    call mpp_get_pelist(ocean_public%domain, pe)
    if (debug > 0) then
      do n = 1,ntiles
        write(tmpstr,'(a,6i6)') subname//' tiles ',n,pe(n),xb(n),xe(n),yb(n),ye(n)
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
      enddo
    end if

    !---------------------------------
    ! create delayout and distgrid
    !---------------------------------

    allocate(deBlockList(2,2,ntiles))
    allocate(petMap(ntiles))
    allocate(deLabelList(ntiles))

    do n = 1, ntiles
      deLabelList(n) = n
      deBlockList(1,1,n) = xb(n)
      deBlockList(1,2,n) = xe(n)
      deBlockList(2,1,n) = yb(n)
      deBlockList(2,2,n) = ye(n)
      petMap(n) = pe(n)
      ! write(tmpstr,'(a,3i8)') subname//' iglo = ',n,deBlockList(1,1,n),deBlockList(1,2,n)
      ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
      ! write(tmpstr,'(a,3i8)') subname//' jglo = ',n,deBlockList(2,1,n),deBlockList(2,2,n)
      ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
      ! write(tmpstr,'(a,2i8)') subname//' pe  = ',n,petMap(n)
      ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
      !--- assume a tile with starting index of 1 has an equivalent wraparound tile on the other side
    enddo

    delayout = ESMF_DELayoutCreate(petMap, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    allocate(connectionList(2))

    ! bipolar boundary condition at top row: nyg
    call ESMF_DistGridConnectionSet(connectionList(1), tileIndexA=1, &
      tileIndexB=1, positionVector=(/nxg+1, 2*nyg+1/), &
      orientationVector=(/-1, -2/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! periodic boundary condition along first dimension
    call ESMF_DistGridConnectionSet(connectionList(2), tileIndexA=1, &
      tileIndexB=1, positionVector=(/nxg, 0/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    distgrid = ESMF_DistGridCreate(minIndex=(/1,1/), maxIndex=(/nxg,nyg/), &
!        indexflag = ESMF_INDEX_DELOCAL, &
        deBlockList=deBlockList, &
!        deLabelList=deLabelList, &
        delayout=delayout, &
        connectionList=connectionList, &
        rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    deallocate(xb,xe,yb,ye,pe)
    deallocate(connectionList)
    deallocate(deLabelList)
    deallocate(deBlockList)
    deallocate(petMap)

    call ESMF_DistGridGet(distgrid=distgrid, localDE=0, elementCount=cnt, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    allocate(indexList(cnt))
    write(tmpstr,'(a,i8)') subname//' distgrid cnt= ',cnt
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call ESMF_DistGridGet(distgrid=distgrid, localDE=0, seqIndexList=indexList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write(tmpstr,'(a,4i8)') subname//' distgrid list= ',&
         indexList(1),indexList(cnt),minval(indexList), maxval(indexList)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    deallocate(IndexList)

    !---------------------------------
    ! create grid
    !---------------------------------

    gridIn = ESMF_GridCreate(distgrid=distgrid, &
      gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/0,1/), &
      coordSys = ESMF_COORDSYS_SPH_DEG, &
      rc = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    mom_grid_i = gridIn

    call ESMF_GridAddCoord(gridIn, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridAddCoord(gridIn, staggerLoc=ESMF_STAGGERLOC_CORNER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridAddItem(gridIn, itemFlag=ESMF_GRIDITEM_MASK, itemTypeKind=ESMF_TYPEKIND_I4, &
       staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Attach area to the Grid optionally. By default the cell areas are computed.
    if(grid_attach_area) then
      call ESMF_GridAddItem(gridIn, itemFlag=ESMF_GRIDITEM_AREA, itemTypeKind=ESMF_TYPEKIND_R8, &
        staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif

    call ESMF_GridGetCoord(gridIn, coordDim=1, &
        staggerloc=ESMF_STAGGERLOC_CENTER, &
        farrayPtr=dataPtr_xcen, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetCoord(gridIn, coordDim=2, &
        staggerloc=ESMF_STAGGERLOC_CENTER, &
        farrayPtr=dataPtr_ycen, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridGetCoord(gridIn, coordDim=1, &
        staggerloc=ESMF_STAGGERLOC_CORNER, &
        farrayPtr=dataPtr_xcor, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetCoord(gridIn, coordDim=2, &
        staggerloc=ESMF_STAGGERLOC_CORNER, &
        farrayPtr=dataPtr_ycor, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridGetItem(gridIn, itemflag=ESMF_GRIDITEM_MASK, &
        staggerloc=ESMF_STAGGERLOC_CENTER, &
        farrayPtr=dataPtr_mask, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if(grid_attach_area) then
      call ESMF_GridGetItem(gridIn, itemflag=ESMF_GRIDITEM_AREA, &
          staggerloc=ESMF_STAGGERLOC_CENTER, &
          farrayPtr=dataPtr_area, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif

    !---------------------------------
    ! load up area, mask, center and corner values
    ! area, mask, and centers should be same size in mom and esmf grid
    ! corner points may not be, need to offset corner points by 1 in i and j
    !   for esmf and also need to "make up" j=1 values.  use wraparound in i
    !---------------------------------

    call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)

    lbnd1 = lbound(dataPtr_mask,1)
    ubnd1 = ubound(dataPtr_mask,1)
    lbnd2 = lbound(dataPtr_mask,2)
    ubnd2 = ubound(dataPtr_mask,2)

    lbnd3 = lbound(dataPtr_xcor,1)
    ubnd3 = ubound(dataPtr_xcor,1)
    lbnd4 = lbound(dataPtr_xcor,2)
    ubnd4 = ubound(dataPtr_xcor,2)

    write(tmpstr,*) subname//' iscjsc = ',isc,iec,jsc,jec
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    write(tmpstr,*) subname//' lbub12 = ',lbnd1,ubnd1,lbnd2,ubnd2
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    write(tmpstr,*) subname//' lbub34 = ',lbnd3,ubnd3,lbnd4,ubnd4
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    if (iec-isc /= ubnd1-lbnd1 .or. jec-jsc /= ubnd2-lbnd2) then
      rc=ESMF_FAILURE
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg=SUBNAME//": fld and grid do not have the same size.", &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      return  ! bail out
    endif

    allocate(ofld(isc:iec,jsc:jec))
    allocate(gfld(nxg,nyg))

    call ocean_model_data_get(ocean_state, ocean_public, 'mask', ofld, isc, jsc)
    write(tmpstr,*) subname//' ofld mask = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(ocean_public%domain, ofld, gfld)
    write(tmpstr,*) subname//' gfld mask = ',minval(gfld),maxval(gfld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    do j = lbnd2, ubnd2
    do i = lbnd1, ubnd1
      j1 = j - lbnd2 + jsc
      i1 = i - lbnd1 + isc
      dataPtr_mask(i,j) = nint(ofld(i1,j1))
    enddo
    enddo

    if(grid_attach_area) then
      call ocean_model_data_get(ocean_state, ocean_public, 'area', ofld, isc, jsc)
      write(tmpstr,*) subname//' ofld area = ',minval(ofld),maxval(ofld)
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
      call mpp_global_field(ocean_public%domain, ofld, gfld)
      write(tmpstr,*) subname//' gfld area = ',minval(gfld),maxval(gfld)
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
      do j = lbnd2, ubnd2
      do i = lbnd1, ubnd1
        j1 = j - lbnd2 + jsc
        i1 = i - lbnd1 + isc
        dataPtr_area(i,j) = ofld(i1,j1)
      enddo
      enddo
    endif

    call ocean_model_data_get(ocean_state, ocean_public, 'tlon', ofld, isc, jsc)
    write(tmpstr,*) subname//' ofld xt = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(ocean_public%domain, ofld, gfld)
    write(tmpstr,*) subname//' gfld xt = ',minval(gfld),maxval(gfld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    do j = lbnd2, ubnd2
    do i = lbnd1, ubnd1
      j1 = j - lbnd2 + jsc
      i1 = i - lbnd1 + isc
      dataPtr_xcen(i,j) = ofld(i1,j1)
      dataPtr_xcen(i,j) = mod(dataPtr_xcen(i,j)+720.0_ESMF_KIND_R8,360.0_ESMF_KIND_R8)
    enddo
    enddo

    call ocean_model_data_get(ocean_state, ocean_public, 'tlat', ofld, isc, jsc)
    write(tmpstr,*) subname//' ofld yt = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(ocean_public%domain, ofld, gfld)
    write(tmpstr,*) subname//' gfld yt = ',minval(gfld),maxval(gfld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    do j = lbnd2, ubnd2
    do i = lbnd1, ubnd1
      j1 = j - lbnd2 + jsc
      i1 = i - lbnd1 + isc
      dataPtr_ycen(i,j) = ofld(i1,j1)
    enddo
    enddo

    call ocean_model_data_get(ocean_state, ocean_public, 'geoLonBu', ofld, isc, jsc)
    write(tmpstr,*) subname//' ofld xu = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(ocean_public%domain, ofld, gfld)
    write(tmpstr,*) subname//' gfld xu = ',minval(gfld),maxval(gfld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    do j = lbnd4, ubnd4
    do i = lbnd3, ubnd3
      j1 = j - lbnd4 + jsc - 1
      i1 = mod(i - lbnd3 + isc - 2 + nxg, nxg) + 1
      if (j1 == 0) then
        dataPtr_xcor(i,j) = 2*gfld(i1,1) - gfld(i1,2)
!        if (dataPtr_xcor(i,j)-dataPtr_xcen(i,j) > 180.) dataPtr_xcor(i,j) = dataPtr_xcor(i,j) - 360.
!        if (dataPtr_xcor(i,j)-dataPtr_xcen(i,j) < 180.) dataPtr_xcor(i,j) = dataPtr_xcor(i,j) + 360.
      elseif (j1 >= 1 .and. j1 <= nyg) then
        dataPtr_xcor(i,j) = gfld(i1,j1)
      else
        rc=ESMF_FAILURE
        call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
          msg=SUBNAME//": error in xu j1.", &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return  ! bail out
      endif
      dataPtr_xcor(i,j) = mod(dataPtr_xcor(i,j)+720.0_ESMF_KIND_R8,360.0_ESMF_KIND_R8)
      ! write(tmpstr,*) subname//' ijfld xu = ',i,i1,j,j1,dataPtr_xcor(i,j)
      ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    enddo
    enddo

    ! MOM6 runs on C-Grid.
    call ocean_model_data_get(ocean_state, ocean_public, 'geoLatBu', ofld, isc, jsc)
    write(tmpstr,*) subname//' ofld yu = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(ocean_public%domain, ofld, gfld)
    write(tmpstr,*) subname//' gfld yu = ',minval(gfld),maxval(gfld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    do j = lbnd4, ubnd4
    do i = lbnd3, ubnd3
      j1 = j - lbnd4 + jsc - 1
      i1 = mod(i - lbnd3 + isc - 2 + nxg, nxg) + 1
      if (j1 == 0) then
        dataPtr_ycor(i,j) = 2*gfld(i1,1) - gfld(i1,2)
      elseif (j1 >= 1 .and. j1 <= nyg) then
        dataPtr_ycor(i,j) = gfld(i1,j1)
      else
        rc=ESMF_FAILURE
        call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
          msg=SUBNAME//": error in yu j1.", &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return  ! bail out
      endif
      ! write(tmpstr,*) subname//' ijfld yu = ',i,i1,j,j1,dataPtr_ycor(i,j)
      ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    enddo
    enddo

    write(tmpstr,*) subname//' mask = ',minval(dataPtr_mask),maxval(dataPtr_mask)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    if(grid_attach_area) then
      write(tmpstr,*) subname//' area = ',minval(dataPtr_area),maxval(dataPtr_area)
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    write(tmpstr,*) subname//' xcen = ',minval(dataPtr_xcen),maxval(dataPtr_xcen)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    write(tmpstr,*) subname//' ycen = ',minval(dataPtr_ycen),maxval(dataPtr_ycen)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    write(tmpstr,*) subname//' xcor = ',minval(dataPtr_xcor),maxval(dataPtr_xcor)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    write(tmpstr,*) subname//' ycor = ',minval(dataPtr_ycor),maxval(dataPtr_ycor)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    deallocate(gfld)

    gridOut = gridIn ! for now out same as in

    !---------------------------------
    ! realize fields on grid
    !---------------------------------


    call MOM_RealizeFields(importState, gridIn , fldsToOcn_num, fldsToOcn, "Ocn import", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call MOM_RealizeFields(exportState, gridOut, fldsFrOcn_num, fldsFrOcn, "Ocn export", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

#ifdef CESMCOUPLED
    call State_SetScalar(dble(nxg),flds_scalar_index_nx, exportState, localPet, &
        flds_scalar_name, flds_scalar_num, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    call State_SetScalar(dble(nyg),flds_scalar_index_ny, exportState, localPet, &
        flds_scalar_name, flds_scalar_num, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
#endif

    call ESMF_StateGet(exportState, itemSearch="sea_surface_temperature", itemCount=icount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Do sst initialization if it's part of export state
    if(icount /= 0) then

      call ESMF_StateGet(exportState, itemName='sea_surface_temperature', field=field_t_surf, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_FieldGet(field_t_surf, localDe=0, farrayPtr=t_surf, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call ocean_model_data_get(ocean_state, ocean_public, 'mask', ofld, isc, jsc)

      lbnd1 = lbound(t_surf,1)
      ubnd1 = ubound(t_surf,1)
      lbnd2 = lbound(t_surf,2)
      ubnd2 = ubound(t_surf,2)

      do j = lbnd2, ubnd2
      do i = lbnd1, ubnd1
        j1 = j - lbnd2 + jsc
        i1 = i - lbnd1 + isc
        if (ofld(i1,j1) == 0.) t_surf(i,j) = 0.0
      enddo
      enddo

      deallocate(ofld)
    endif

! tcraig, turn this off for now, have issues with overwriting failures
!    call NUOPC_Write(exportState, fileNamePrefix='init_field_ocn_export_', &
!      timeslice=1, relaxedFlag=.true., rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    write(*,*) '----- MOM initialization phase Realize completed'

  end subroutine InitializeRealize

  !===============================================================================

  subroutine DataInitialize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)                       :: clock
    type(ESMF_State)                       :: importState, exportState
    type (ocean_public_type),      pointer :: ocean_public       => NULL()
    type (ocean_state_type),       pointer :: ocean_state        => NULL()
    type(ice_ocean_boundary_type), pointer :: Ice_ocean_boundary => NULL()
    type(ocean_internalstate_wrapper)      :: ocean_internalstate
    type(ocean_grid_type), pointer         :: ocean_grid
    character(240)                         :: msgString
    integer                                :: fieldCount, n
    type(ESMF_Field)                       :: field
    character(len=64),allocatable          :: fieldNameList(:)
    character(len=*),parameter  :: subname='(mom_cap:DataInitialize)'

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridCompGetInternalState(gcomp, ocean_internalstate, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    Ice_ocean_boundary => ocean_internalstate%ptr%ice_ocean_boundary_type_ptr
    ocean_public       => ocean_internalstate%ptr%ocean_public_type_ptr
    ocean_state        => ocean_internalstate%ptr%ocean_state_type_ptr
    call get_ocean_grid(ocean_state, ocean_grid)

    call mom_export(ocean_public, ocean_grid, exportState, logunit, clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    allocate(fieldNameList(fieldCount))
    call ESMF_StateGet(exportState, itemNameList=fieldNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    do n=1, fieldCount
      call ESMF_StateGet(exportState, itemName=fieldNameList(n), field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call NUOPC_SetAttribute(field, name="Updated", value="true", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end do
    deallocate(fieldNameList)

    ! check whether all Fields in the exportState are "Updated"
    if (NUOPC_IsUpdated(exportState)) then
      call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="true", rc=rc)

      call ESMF_LogWrite("MOM6 - Initialize-Data-Dependency SATISFIED!!!", ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end if

    if(write_diagnostics) then
      call NUOPC_Write(exportState, fileNamePrefix='field_init_ocn_export_', &
        timeslice=import_slice, relaxedFlag=.true., rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif

  end subroutine DataInitialize

  !===============================================================================

  !> Called by NUOPC to advance the model a single timestep.
  !!
  !! @param gcomp an ESMF_GridComp object
  !! @param rc return code
  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)                    :: gcomp
    integer, intent(out)                   :: rc

    ! local variables
    type(ESMF_Clock)                       :: clock
    type(ESMF_Alarm)                       :: alarm
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Time)                        :: currTime
    type(ESMF_TimeInterval)                :: timeStep
    type(ESMF_Time)                        :: startTime
    type(ESMF_TimeInterval)                :: time_elapsed
    integer(ESMF_KIND_I8)                  :: n_interval, time_elapsed_sec
    character(len=64)                      :: timestamp
    type (ocean_public_type),      pointer :: ocean_public       => NULL()
    type (ocean_state_type),       pointer :: ocean_state        => NULL()
    type(ice_ocean_boundary_type), pointer :: Ice_ocean_boundary => NULL()
    type(ocean_internalstate_wrapper)      :: ocean_internalstate
    type(time_type)                        :: Time
    type(time_type)                        :: Time_step_coupled
    type(time_type)                        :: Time_restart_current
    integer                                :: dth, dtm, dts, dt_cpld  = 86400
    integer                                :: isc,iec,jsc,jec,lbnd1,ubnd1,lbnd2,ubnd2
    integer                                :: i,j,i1,j1
    integer                                :: nc
#ifdef CESMCOUPLED
    type(ESMF_Time)                        :: MyTime
    integer                                :: seconds, day, year, month, hour, minute
    integer                                :: shrlogunit       ! original log unit
    integer                                :: shrloglev        ! original log level
    integer                                :: logunit          ! i/o unit for stdout
    integer                                :: nu               ! i/o unit to write pointer file
    character(ESMF_MAXSTR)                 :: cvalue
    character(ESMF_MAXSTR)                 :: runid            ! Run ID
    character(len=384)                     :: restartname      ! restart file name (no dir)
    character(len=384)                     :: restart_pointer_file ! file name for restart pointer file
#else
    real(ESMF_KIND_R8), allocatable        :: ofld(:,:), ocz(:,:), ocm(:,:)
    real(ESMF_KIND_R8), allocatable        :: mmmf(:,:), mzmf(:,:)
    real(ESMF_KIND_R8), pointer            :: dataPtr_mask(:,:)
    real(ESMF_KIND_R8), pointer            :: dataPtr_mmmf(:,:)
    real(ESMF_KIND_R8), pointer            :: dataPtr_mzmf(:,:)
    real(ESMF_KIND_R8), pointer            :: dataPtr_ocz(:,:)
    real(ESMF_KIND_R8), pointer            :: dataPtr_ocm(:,:)
    real(ESMF_KIND_R8), pointer            :: dataPtr_frazil(:,:)
    real(ESMF_KIND_R8), pointer            :: dataPtr_evap(:,:)
    real(ESMF_KIND_R8), pointer            :: dataPtr_sensi(:,:)
#endif
    type(ocean_grid_type), pointer         :: ocean_grid
    character(240)                         :: msgString
    character(len=*),parameter             :: subname='(mom_cap:ModelAdvance)'
    !--------------------------------

    rc = ESMF_SUCCESS
    if(profile_memory) call ESMF_VMLogMemInfo("Entering MOM Model_ADVANCE: ")

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridCompGetInternalState(gcomp, ocean_internalstate, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    Ice_ocean_boundary => ocean_internalstate%ptr%ice_ocean_boundary_type_ptr
    ocean_public       => ocean_internalstate%ptr%ocean_public_type_ptr
    ocean_state        => ocean_internalstate%ptr%ocean_state_type_ptr

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep

    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing OCN from: ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(subname//trim(msgString), ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, &
      timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimePrint(currTime + timeStep, &
      preString="--------------------------------> to: ", &
      unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimeIntervalGet(timeStep, h=dth, m=dtm, s=dts, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    Time = esmf2fms_time(currTime)
    Time_step_coupled = esmf2fms_time(timeStep)
    dt_cpld = dth*3600+dtm*60+dts

    call ice_ocn_bnd_from_data(Ice_ocean_boundary, Time, Time_step_coupled)

    call external_coupler_sbc_before(Ice_ocean_boundary, ocean_public, nc, dt_cpld )

    if(write_diagnostics) then
      call NUOPC_Write(importState, fileNamePrefix='field_ocn_import_', &
        timeslice=import_slice, relaxedFlag=.true., rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      import_slice = import_slice + 1
    endif

    ! rotate the lat/lon wind vector (CW) onto local tripolar coordinate system

    call mpp_get_compute_domain(ocean_public%domain, isc, iec, jsc, jec)

    call get_ocean_grid(ocean_state, ocean_grid)

#ifdef CESMCOUPLED
    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogLevel(max(shrloglev,1))
    call shr_file_setLogUnit (logunit)

    call mom_import(ocean_public, ocean_grid, importState, ice_ocean_boundary, logunit, runtype, clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#else
    call State_getFldPtr(exportState,'ocean_mask',dataPtr_mask,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    lbnd1 = lbound(dataPtr_mask,1)
    ubnd1 = ubound(dataPtr_mask,1)
    lbnd2 = lbound(dataPtr_mask,2)
    ubnd2 = ubound(dataPtr_mask,2)

    call State_getFldPtr(importState,'mean_zonal_moment_flx',dataPtr_mzmf,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,'mean_merid_moment_flx',dataPtr_mmmf,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,'mean_evap_rate',dataPtr_evap,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(importState,'mean_sensi_heat_flx',dataPtr_sensi,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    dataPtr_evap = - dataPtr_evap
    dataPtr_sensi = - dataPtr_sensi

    allocate(mzmf(lbnd1:ubnd1,lbnd2:ubnd2))
    allocate(mmmf(lbnd1:ubnd1,lbnd2:ubnd2))
    do j  = lbnd2, ubnd2
      do i = lbnd1, ubnd1
        j1 = j - lbnd2 + jsc  ! work around local vs global indexing
        i1 = i - lbnd1 + isc
        mzmf(i,j) = ocean_grid%cos_rot(i1,j1)*dataPtr_mzmf(i,j) &
                  + ocean_grid%sin_rot(i1,j1)*dataPtr_mmmf(i,j)
        mmmf(i,j) = ocean_grid%cos_rot(i1,j1)*dataPtr_mmmf(i,j) &
                  - ocean_grid%sin_rot(i1,j1)*dataPtr_mzmf(i,j)
      enddo
    enddo
    dataPtr_mzmf = mzmf
    dataPtr_mmmf = mmmf
    deallocate(mzmf, mmmf)

    !Optionally write restart files when currTime-startTime is integer multiples of restart_interval
    if (restart_interval > 0 ) then
      time_elapsed = currTime - startTime
      call ESMF_TimeIntervalGet(time_elapsed, s_i8=time_elapsed_sec, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      n_interval = time_elapsed_sec / restart_interval
      if ((n_interval .gt. 0) .and. (n_interval*restart_interval == time_elapsed_sec)) then
          time_restart_current = esmf2fms_time(currTime)
          timestamp = date_to_string(time_restart_current)
          call ESMF_LogWrite("MOM: Writing restart at "//trim(timestamp), ESMF_LOGMSG_INFO, rc=dbrc)
          write(*,*) 'calling ocean_model_restart'
          call ocean_model_restart(ocean_state, timestamp)
      endif
    endif
#endif

    ! Update MOM6

    if(profile_memory) call ESMF_VMLogMemInfo("Entering MOM update_ocean_model: ")
    call update_ocean_model(Ice_ocean_boundary, ocean_state, ocean_public, Time, Time_step_coupled)
    if(profile_memory) call ESMF_VMLogMemInfo("Leaving MOM update_ocean_model: ")

#ifdef CESMCOUPLED

    call mom_export(ocean_public, ocean_grid, exportState, logunit, clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! If restart alarm is ringing - write restart file
    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call ESMF_AlarmRingerOff( alarm, rc=rc )
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! determine restart filename
      ! Need to use next time step since clock is not advanced until the end of the time interval
      call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      read(cvalue,*) runid

      call ESMF_ClockGetNextTime(clock, MyTime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call ESMF_TimeGet (MyTime, yy=year, mm=month, dd=day, s=seconds, rc=rc )
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      write(restartname,'(A,".mom6.r.",I4.4,"-",I2.2,"-",I2.2,"-",I5.5)') trim(runid), year, month, day, seconds

      ! write name of restart file in the rpointer file
      nu = shr_file_getUnit()
      if (is_root_pe()) then
        restart_pointer_file = 'rpointer.ocn'
        open(nu, file=restart_pointer_file, form='formatted', status='unknown')
        write(nu,'(a)') trim(restartname) //'.nc'
        close(nu)
        write(logunit,*) 'ocn restart pointer file written: ',trim(restartname)
      endif
      call shr_file_freeUnit(nu)

      ! write restart file(s)
      call ocean_model_restart(ocean_state, restartname=restartname)

      if (is_root_pe()) then
        write(logunit,*) subname//' writing restart file ',trim(restartname)
      end if
    endif

    ! reset shr logging to my original values
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#else

    allocate(ofld(isc:iec,jsc:jec))

    call ocean_model_data_get(ocean_state, ocean_public, 'mask', ofld, isc, jsc)
    do j = lbnd2, ubnd2
    do i = lbnd1, ubnd1
      j1 = j - lbnd2 + jsc
      i1 = i - lbnd1 + isc
      dataPtr_mask(i,j) = nint(ofld(i1,j1))
    enddo
    enddo
    deallocate(ofld)

    ! Now rotate ocn current from tripolar grid back to lat/lon grid (CCW)
    allocate(ocz(lbnd1:ubnd1,lbnd2:ubnd2))
    allocate(ocm(lbnd1:ubnd1,lbnd2:ubnd2))

    call State_getFldPtr(exportState,'ocn_current_zonal',dataPtr_ocz,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(exportState,'ocn_current_merid',dataPtr_ocm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(exportState,'freezing_melting_potential',dataPtr_frazil,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    dataPtr_frazil = dataPtr_frazil/dt_cpld !convert from J/m^2 to W/m^2 for CICE coupling

    ocz = dataPtr_ocz
    ocm = dataPtr_ocm
    do j  = lbnd2, ubnd2
      do i = lbnd1, ubnd1
        j1 = j - lbnd2 + jsc  ! work around local vs global indexing
        i1 = i - lbnd1 + isc
        dataPtr_ocz(i,j) = ocean_grid%cos_rot(i1,j1)*ocz(i,j) &
                         - ocean_grid%sin_rot(i1,j1)*ocm(i,j)
        dataPtr_ocm(i,j) = ocean_grid%cos_rot(i1,j1)*ocm(i,j) &
                         + ocean_grid%sin_rot(i1,j1)*ocz(i,j)
      enddo
    enddo
    deallocate(ocz, ocm)

    !call ESMF_LogWrite("B4 writing diags", dataPtr_model_data_get(ocean_state, ocean_public, 'mask', ofld, isc, jsc))
    do j = lbnd2, ubnd2
    do i = lbnd1, ubnd1
      j1 = j - lbnd2 + jsc
      i1 = i - lbnd1 + isc
      dataPtr_mask(i,j) = nint(ofld(i1,j1))
    enddo
    enddo
    deallocate(ofld)

    ! Now rotate ocn current from tripolar grid back to lat/lon grid (CCW)
    allocate(ocz(lbnd1:ubnd1,lbnd2:ubnd2))
    allocate(ocm(lbnd1:ubnd1,lbnd2:ubnd2))

    call State_getFldPtr(exportState,'ocn_current_zonal',dataPtr_ocz,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(exportState,'ocn_current_merid',dataPtr_ocm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call State_getFldPtr(exportState,'freezing_melting_potential',dataPtr_frazil,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    dataPtr_frazil = dataPtr_frazil/dt_cpld !convert from J/m^2 to W/m^2 for CICE coupling

    ocz = dataPtr_ocz
    ocm = dataPtr_ocm
    do j  = lbnd2, ubnd2
      do i = lbnd1, ubnd1
        j1 = j - lbnd2 + jsc  ! work around local vs global indexing
        i1 = i - lbnd1 + isc
        dataPtr_ocz(i,j) = ocean_grid%cos_rot(i1,j1)*ocz(i,j) &
                         - ocean_grid%sin_rot(i1,j1)*ocm(i,j)
        dataPtr_ocm(i,j) = ocean_grid%cos_rot(i1,j1)*ocm(i,j) &
                         + ocean_grid%sin_rot(i1,j1)*ocz(i,j)
      enddo
    enddo
    deallocate(ocz, ocm)

#endif

    if (write_diagnostics) then
      call NUOPC_Write(exportState, fileNamePrefix='field_ocn_export_', &
        timeslice=export_slice, relaxedFlag=.true., rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      export_slice = export_slice + 1
    endif

    !call ESMF_LogWrite("Before calling sbc forcing", ESMF_LOGMSG_INFO, rc=rc)
    call external_coupler_sbc_after(Ice_ocean_boundary, ocean_public, nc, dt_cpld )
    !call ESMF_LogWrite("Before dumpMomInternal", ESMF_LOGMSG_INFO, rc=rc)
    !write(*,*) 'MOM: --- run phase called ---'

    if(profile_memory) call ESMF_VMLogMemInfo("Leaving MOM Model_ADVANCE: ")

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=128)       :: mtimestring, dtimestring
    character(len=256)       :: cvalue
    character(len=256)       :: restart_option ! Restart option units
    integer                  :: restart_n      ! Number until restart interval
    integer                  :: restart_ymd    ! Restart date (YYYYMMDD)
    type(ESMF_ALARM)         :: restart_alarm
    logical                  :: first_time = .true.
    character(len=*),parameter :: subname='mom_cap:(ModelSetRunClock) '
    !--------------------------------

    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !--------------------------------
    ! check that the current time in the model and driver are the same
    !--------------------------------

    if (mcurrtime /= dcurrtime) then
      call ESMF_TimeGet(dcurrtime, timeString=dtimestring, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call ESMF_TimeGet(mcurrtime, timeString=mtimestring, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call ESMF_LogWrite(subname//" ERROR in time consistency; "//trim(dtimestring)//" ne "//trim(mtimestring),  &
          ESMF_LOGMSG_ERROR, rc=dbrc)
      rc = ESMF_FAILURE
      return
    endif

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (first_time) then
       !--------------------------------                                                                                 
       ! set restart alarm
       !--------------------------------                                                                                 
       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       call NUOPC_CompAttributeGet(gcomp,  name="restart_n", value=cvalue, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
       read(cvalue,*) restart_ymd

       call shr_nuopc_time_alarmInit(mclock, &
            alarm   = restart_alarm,         &
            option  = trim(restart_option),  &
            opt_n   = restart_n,             &
            opt_ymd = restart_ymd,           &
            RefTime = mcurrTime,             &
            alarmname = 'alarm_restart', rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
       first_time = .false.
    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine ModelSetRunClock

  !===============================================================================

  !> Called by NUOPC at the end of the run to clean up.
  !!
  !! @param gcomp an ESMF_GridComp object
  !! @param rc return code
  subroutine ocean_model_finalize(gcomp, rc)

    ! input arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type (ocean_public_type),      pointer :: ocean_public
    type (ocean_state_type),       pointer :: ocean_state
    type(ocean_internalstate_wrapper)      :: ocean_internalstate
    type(TIME_TYPE)                        :: Time
    type(ESMF_Clock)                       :: clock
    type(ESMF_Time)                        :: currTime
    character(len=64)                      :: timestamp
    character(len=*),parameter  :: subname='(mom_cap:ocean_model_finalize)'

    write(*,*) 'MOM: --- finalize called ---'
    rc = ESMF_SUCCESS

    call ESMF_GridCompGetInternalState(gcomp, ocean_internalstate, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ocean_public => ocean_internalstate%ptr%ocean_public_type_ptr
    ocean_state  => ocean_internalstate%ptr%ocean_state_type_ptr

    call NUOPC_ModelGet(gcomp, modelClock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    Time = esmf2fms_time(currTime)

#ifdef CESMCOUPLED
    call ocean_model_end (ocean_public, ocean_State, Time, write_restart=.false.)
#else
    call ocean_model_end (ocean_public, ocean_State, Time, write_restart=.true.)
#endif
    call field_manager_end

    call fms_io_exit
    call fms_end

    write(*,*) 'MOM: --- completed ---'

  end subroutine ocean_model_finalize

  !====================================================================

  ! get forcing data from data_overide
  subroutine ice_ocn_bnd_from_data(x, Time, Time_step_coupled)

    type (ice_ocean_boundary_type) :: x
    type(Time_type), intent(in)    :: Time, Time_step_coupled

    type(Time_type)                :: Time_next
    character(len=*),parameter  :: subname='(mom_cap:ice_ocn_bnd_from_data)'

    Time_next = Time + Time_step_coupled

    !call data_override('OCN', 't_flux',          x%t_flux         , Time_next)
    !call data_override('OCN', 'u_flux',          x%u_flux         , Time_next)
    !call data_override('OCN', 'v_flux',          x%v_flux         , Time_next)
    !call data_override('OCN', 'q_flux',          x%q_flux         , Time_next)
    !call data_override('OCN', 'salt_flux',       x%salt_flux      , Time_next)
    !call data_override('OCN', 'lw_flux',         x%lw_flux        , Time_next)
    !call data_override('OCN', 'sw_flux_vis_dir', x%sw_flux_vis_dir, Time_next)
    !call data_override('OCN', 'sw_flux_vis_dif', x%sw_flux_vis_dif, Time_next)
    !call data_override('OCN', 'sw_flux_nir_dir', x%sw_flux_nir_dir, Time_next)
    !call data_override('OCN', 'sw_flux_nir_dif', x%sw_flux_nir_dif, Time_next)
    !call data_override('OCN', 'lprec',           x%lprec          , Time_next)
    !call data_override('OCN', 'fprec',           x%fprec          , Time_next)
    !call data_override('OCN', 'runoff',          x%runoff         , Time_next)
    !call data_override('OCN', 'calving',         x%calving        , Time_next)
    !call data_override('OCN', 'p',               x%p              , Time_next)

  end subroutine ice_ocn_bnd_from_data

!-----------------------------------------------------------------------------------------
!
! Subroutines  for enabling coupling to external programs through a third party coupler
! such as OASIS/PRISM.
! If no external coupler then these will mostly be dummy routines.
! These routines can also serve as spots to call other user defined routines
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------

! Dummy subroutines.

  subroutine external_coupler_mpi_init(mom_local_communicator, external_initialization)
  implicit none
  integer, intent(out) :: mom_local_communicator
  logical, intent(out) :: external_initialization
  external_initialization = .false.
  mom_local_communicator = -100         ! Is there mpp_undefined parameter corresponding to MPI_UNDEFINED?
                                        ! probably wouldn't need logical flag.
  return
  end subroutine external_coupler_mpi_init

!-----------------------------------------------------------------------------------------
  subroutine external_coupler_sbc_init(Dom, dt_cpld, Run_len)
  implicit none
  type(domain2d) :: Dom
  integer :: dt_cpld
  type(time_type) :: Run_len
  return
  end  subroutine external_coupler_sbc_init

  subroutine external_coupler_sbc_before(Ice_ocean_boundary, ocean_public, nsteps, dt_cpld )
  implicit none
  type (ice_ocean_boundary_type), intent(INOUT) :: Ice_ocean_boundary
  type (ocean_public_type) , intent(INOUT)        :: ocean_public
  integer , intent(IN)                       :: nsteps, dt_cpld
  return
  end subroutine external_coupler_sbc_before


  subroutine external_coupler_sbc_after(Ice_ocean_boundary, ocean_public, nsteps, dt_cpld )
  type (ice_ocean_boundary_type) :: Ice_ocean_boundary
  type (ocean_public_type)         :: ocean_public
  integer                        :: nsteps, dt_cpld
  return
  end subroutine external_coupler_sbc_after

  subroutine external_coupler_restart( dt_cpld, num_cpld_calls )
  implicit none
  integer, intent(in)               :: dt_cpld, num_cpld_calls
  return
  end subroutine external_coupler_restart

  subroutine external_coupler_exit
  return
  end subroutine external_coupler_exit

!-----------------------------------------------------------------------------------------
  subroutine external_coupler_mpi_exit(mom_local_communicator, external_initialization)
  implicit none
  integer, intent(in) :: mom_local_communicator
  logical, intent(in) :: external_initialization
  return
  end subroutine external_coupler_mpi_exit
!-----------------------------------------------------------------------------------------
    subroutine writeSliceFields(state, filename_prefix, slice, rc)
      type(ESMF_State)                :: state
      character(len=*)                :: filename_prefix
      integer                         :: slice
      integer, intent(out), optional  :: rc

      integer                         :: n, nfields
      type(ESMF_Field)                :: field
      type(ESMF_StateItem_Flag)       :: itemType
      character(len=40)               :: fileName
      character(len=64),allocatable   :: fieldNameList(:)
      character(len=*),parameter :: subname='(mom_cap:writeSliceFields)'

      if (present(rc)) rc = ESMF_SUCCESS

      if (ESMF_IO_PIO_PRESENT .and. &
        (ESMF_IO_NETCDF_PRESENT .or. ESMF_IO_PNETCDF_PRESENT)) then

        call ESMF_StateGet(state, itemCount=nfields, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        allocate(fieldNameList(nfields))
        call ESMF_StateGet(state, itemNameList=fieldNameList, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        do n=1, size(fieldNameList)
          call ESMF_StateGet(state, itemName=fieldNameList(n), &
            itemType=itemType, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          if (itemType /= ESMF_STATEITEM_NOTFOUND) then
            ! field is available in the state
            call ESMF_StateGet(state, itemName=fieldNameList(n), field=field, &
              rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              return  ! bail out
            ! -> output to file
            write (fileName,"(A)") &
              filename_prefix//trim(fieldNameList(n))//".nc"
            call ESMF_FieldWrite(field, fileName=trim(fileName), &
              timeslice=slice, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__)) &
              call ESMF_Finalize(endflag=ESMF_END_ABORT)
          endif
        enddo

        deallocate(fieldNameList)

      endif

    end subroutine writeSliceFields

  !-----------------------------------------------------------------------------

  subroutine State_GetFldPtr(ST, fldname, fldptr, rc)
    type(ESMF_State)   , intent(in)            :: ST
    character(len=*)   , intent(in)            :: fldname
    real(ESMF_KIND_R8) , pointer, intent(in)   :: fldptr(:,:)
    integer            , intent(out), optional :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer          :: lrc
    character(len=*),parameter :: subname='(mom_cap:State_GetFldPtr)'

    call ESMF_StateGet(ST, itemName=trim(fldname), field=lfield, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (present(rc)) rc = lrc

  end subroutine State_GetFldPtr

  !-----------------------------------------------------------------------------

  subroutine State_SetScalar(value, scalar_id, State, mytask, scalar_name, scalar_num,  rc)
    ! ----------------------------------------------
    ! Set scalar data from State for a particular name
    ! ----------------------------------------------
    real(ESMF_KIND_R8),intent(in)     :: value
    integer,           intent(in)     :: scalar_id
    type(ESMF_State),  intent(inout)  :: State
    integer,           intent(in)     :: mytask
    character(len=*),  intent(in)     :: scalar_name   
    integer,           intent(in)     :: scalar_num 
    integer,           intent(inout)  :: rc

    ! local variables
    integer                         :: ierr, len
    type(ESMF_Field)                :: field
    real(ESMF_KIND_R8), pointer     :: farrayptr(:,:)
    character(len=*), parameter     :: subname='(mom_cap:State_SetScalar)'

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(scalar_name), field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    if (mytask == 0) then
      call ESMF_FieldGet(field, farrayPtr=farrayptr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      if (scalar_id < 0 .or. scalar_id > scalar_num) then
         call ESMF_LogWrite(trim(subname)//": ERROR in scalar_id", &
              ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u, rc=dbrc)
        rc = ESMF_FAILURE
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
      endif

      farrayptr(1,scalar_id) = value
    endif

  end subroutine State_SetScalar

  !-----------------------------------------------------------------------------

  subroutine MOM_RealizeFields(state, grid, nfields, field_defs, tag, rc)

    type(ESMF_State)    , intent(inout) :: state
    type(ESMF_Grid)     , intent(in)    :: grid
    integer             , intent(in)    :: nfields
    type(fld_list_type) , intent(inout) :: field_defs(:)
    character(len=*)    , intent(in)    :: tag
    integer             , intent(inout) :: rc

    integer          :: i
    type(ESMF_Field) :: field
    integer          :: npet, nx, ny, pet
    integer          :: elb(2), eub(2), clb(2), cub(2), tlb(2), tub(2)
    type(ESMF_VM)    :: vm
    character(len=*),parameter  :: subname='(mom_cap:MOM_RealizeFields)'

    rc = ESMF_SUCCESS

    do i = 1, nfields

      if (NUOPC_IsConnected(state, fieldName=field_defs(i)%shortname)) then

        if (field_defs(i)%shortname == flds_scalar_name) then
          call ESMF_LogWrite(subname // tag // " Field "// trim(field_defs(i)%stdname) // " is connected on root pe.", &
            ESMF_LOGMSG_INFO, &
            line=__LINE__, &
            file=__FILE__, &
            rc=dbrc)
          call SetScalarField(field, rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        elseif (field_defs(i)%assoc) then
          call ESMF_LogWrite(subname // tag // " Field "// trim(field_defs(i)%stdname)&
            // " is connected and associated.", &
            ESMF_LOGMSG_INFO, &
            line=__LINE__, &
            file=__FILE__, &
            rc=dbrc)
          write(tmpstr,'(a,4i12)') subname//trim(tag)//' Field '//trim(field_defs(i)%shortname)//':', &
            lbound(field_defs(i)%farrayPtr,1), ubound(field_defs(i)%farrayPtr,1), &
            lbound(field_defs(i)%farrayPtr,2), ubound(field_defs(i)%farrayPtr,2)
          call ESMF_LogWrite(tmpstr, ESMF_LOGMSG_INFO, rc=dbrc)
          field = ESMF_FieldCreate(grid=grid, &
            farray=field_defs(i)%farrayPtr, indexflag=ESMF_INDEX_DELOCAL, &
           !farray=field_defs(i)%farrayPtr, indexflag=ESMF_INDEX_GLOBAL, &
            name=field_defs(i)%shortname, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        else
          call ESMF_LogWrite(subname // tag // " Field "// trim(field_defs(i)%stdname) // " is connected.", &
            ESMF_LOGMSG_INFO, &
            line=__LINE__, &
            file=__FILE__, &
            rc=dbrc)
          field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, indexflag=ESMF_INDEX_DELOCAL, &
            name=field_defs(i)%shortname, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        endif

        call NUOPC_Realize(state, field=field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        !        call ESMF_FieldPrint(field=field, rc=rc)
        !        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        !          line=__LINE__, &
        !          file=__FILE__)) &
        !          return  ! bail out
      else
        call ESMF_LogWrite(subname // tag // " Field "// trim(field_defs(i)%stdname) // " is not connected.", &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=dbrc)
        ! TODO: Initialize the value in the pointer to 0 after proper restart is setup
        ! if(associated(field_defs(i)%farrayPtr) ) field_defs(i)%farrayPtr = 0.0
        ! remove a not connected Field from State
        call ESMF_StateRemove(state, (/field_defs(i)%shortname/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

    enddo

  end subroutine MOM_RealizeFields

  !-----------------------------------------------------------------------------

  subroutine SetScalarField(field, rc)
    ! ----------------------------------------------
    ! create a field with scalar data on the root pe
    ! ----------------------------------------------
    type(ESMF_Field), intent(inout)  :: field
    integer,          intent(inout)  :: rc

    ! local variables
    type(ESMF_Distgrid) :: distgrid
    type(ESMF_Grid)     :: grid
    character(len=*), parameter :: subname='(mom_cap:SetScalarField)'

    rc = ESMF_SUCCESS

    ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
    distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    grid = ESMF_GridCreate(distgrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, &
      typekind=ESMF_TYPEKIND_R8, &
      ungriddedLBound=(/1/), &
      ungriddedUBound=(/flds_scalar_num/), & ! num of scalar values
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine SetScalarField

  !-----------------------------------------------------------------------------

  subroutine fld_list_add(num, fldlist, stdname, transferOffer, data, shortname)
    ! ----------------------------------------------
    ! Set up a list of field information
    ! ----------------------------------------------
    integer,             intent(inout)  :: num
    type(fld_list_type), intent(inout)  :: fldlist(:)
    character(len=*),    intent(in)     :: stdname
    character(len=*),    intent(in)     :: transferOffer
    real(ESMF_KIND_R8), dimension(:,:), optional, target :: data
    character(len=*),    intent(in),optional :: shortname

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(mom_cap:fld_list_add)'

    ! fill in the new entry

    num = num + 1
    if (num > fldsMax) then
      call ESMF_LogWrite(trim(subname)//": ERROR num gt fldsMax "//trim(stdname), &
        ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=dbrc)
      return
    endif

    fldlist(num)%stdname        = trim(stdname)
    if (present(shortname)) then
       fldlist(num)%shortname   = trim(shortname)
    else
       fldlist(num)%shortname   = trim(stdname)
    endif
    fldlist(num)%transferOffer  = trim(transferOffer)
    if (present(data)) then
      fldlist(num)%assoc        = .true.
      ! The following sets up the data pointer that will be used in the realize call
      fldlist(num)%farrayPtr    => data
    else
      fldlist(num)%assoc        = .false.
    endif

  end subroutine fld_list_add

  !-----------------------------------------------------------------------------

#if (1 == 0)
  subroutine calculate_rot_angle(OS, OSFC)
    type(ocean_state_type), intent(in)    :: OS
    type(ocean_public_type), intent(in)   :: OSFC

    integer                               :: i,j,ishift,jshift,ilb,iub,jlb,jub
    real                                  :: angle, lon_scale
    type(ocean_grid_type), pointer        :: grid

    call get_ocean_grid(OS, grid)

    !print *, 'lbound: ', lbound(grid%geoLatT), lbound(grid%geoLonT), lbound(grid%sin_rot)
    !print *, 'ubound: ', ubound(grid%geoLatT), ubound(grid%geoLonT), ubound(grid%sin_rot)

    !print *, minval(grid%geoLatT), maxval(grid%geoLatT)
    !print *, minval(grid%geoLonT), maxval(grid%geoLonT)
    !print *, grid%isc, grid%jsc, grid%iec, grid%jec

    !
    ! The bounds isc:iec goes from 5-104, isc-ishift:iec-ishift goes from 1:100
    !
    call mpp_get_compute_domain(OSFC%Domain, ilb, iub, jlb, jub)
    ishift = ilb-grid%isc
    jshift = jlb-grid%jsc
    !print *, 'ilb, iub, jlb, jub', ilb, iub, jlb, jub, ishift, jshift
    !print *, 'sizes', iub-ilb, jub-jlb, grid%iec-grid%isc, grid%jec-grid%jsc
!    allocate(grid%sin_rot(ilb:iub, jlb:jub))
!    allocate(grid%cos_rot(ilb:iub, jlb:jub))

    ! loop 5-104
    do j=grid%jsc,grid%jec ; do i=grid%isc,grid%iec
      lon_scale    = cos((grid%geoLatBu(I-1,J-1) + grid%geoLatBu(I,J-1  ) + &
                          grid%geoLatBu(I-1,J) + grid%geoLatBu(I,J)) * atan(1.0)/180)
      angle        = atan2((grid%geoLonBu(I-1,J) + grid%geoLonBu(I,J) - &
                            grid%geoLonBu(I-1,J-1) - grid%geoLonBu(I,J-1))*lon_scale, &
                            grid%geoLatBu(I-1,J) + grid%geoLatBu(I,J) - &
                            grid%geoLatBu(I-1,J-1) - grid%geoLatBu(I,J-1) )
      grid%sin_rot(i+ishift,j+jshift) = sin(angle) ! angle is the clockwise angle from lat/lon to ocean
      grid%cos_rot(i+ishift,j+jshift) = cos(angle) ! grid (e.g. angle of ocean "north" from true north)
    enddo ; enddo
    !print *, minval(grid%sin_rot), maxval(grid%sin_rot)
    !print *, minval(grid%cos_rot), maxval(grid%cos_rot)

  end subroutine
#endif

end module mom_cap_mod
