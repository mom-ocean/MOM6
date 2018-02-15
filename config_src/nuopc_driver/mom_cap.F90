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
!! **This MOM cap has been tested with MOM5 and MOM6.**
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
!! (http://www.earthsystemmodeling.org/esmf_releases/non_public/ESMF_7_0_0/NUOPC_refdoc/node3.html#SECTION00034000000000000000).
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
!! (http://www.earthsystemmodeling.org/esmf_releases/non_public/ESMF_7_0_0/NUOPC_refdoc/node3.html#SECTION00032000000000000000).
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
!! ---------|--------------------------------------------------------------------|-------------------------------------------------------------
!! Init     | [InitializeP0] (@ref mom_cap_mod::initializep0)                    | Sets the Initialize Phase Definition (IPD) version to use
!! Init     | [InitializeAdvertise] (@ref mom_cap_mod::initializeadvertise)      | Advertises standard names of import and export fields
!! Init     | [InitializeRealize] (@ref mom_cap_mod::initializerealize)          | Creates an ESMF_Grid for the MOM grid as well as ESMF_Fields for import and export fields
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
!!      call update_ocean_model(Ice_ocean_boundary, Ocean_state, Ocean_sfc, Time, Time_step_coupled)
!!
!! Prior to this call, the cap performs a few steps:
!! - the `Time` and `Time_step_coupled` parameters, which are based on FMS types, are derived from the incoming ESMF clock
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
!! - a call is made to `external_coupler_sbc_after()` to update exports from an external coupler (currently an inactive stub)
!! - calls are made to `dumpMomInternal()` to write files `field_ocn_internal_*` for all internal fields (both import and export)
!!
!! @subsubsection VectorRotations Vector Rotations
!!
!! Vector rotations are applied to incoming momentum fluxes (from regular lat-lon to tripolar grid) and
!! outgoing ocean currents (from tripolar to regular lat-lon). The rotation angles are provided
!! from the native MOM grid by a call to `get_ocean_grid(Ocean_grid)`.
!! The cosine and sine of the rotation angle are:
!!
!!     Ocean_grid%cos_rot(i,j)
!!     Ocean_grid%sin_rot(i,j)
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
!!     call ocean_model_end (Ocean_sfc, Ocean_State, Time)
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
!! Standard Name                     | Units      | Model Variable  | Description                                   | Notes
!! ----------------------------------|------------|-----------------|-----------------------------------------------|--------------------------------------
!! inst_pres_height_surface          | Pa         | p               | pressure of overlying sea ice and atmosphere  | |
!! mass_of_overlying_sea_ice         | kg         | mi              | mass of overlying sea ice                     | |
!! mean_calving_heat_flx             | W m-2      | calving_hflx    | heat flux, relative to 0C, of frozen land water into ocean | |
!! mean_calving_rate                 | kg m-2 s-1 | calving         | mass flux of frozen runoff                    | |
!! mean_evap_rate                    | kg m-2 s-1 | q_flux          | specific humidity flux                        | sign reversed (- evap)
!! mean_fprec_rate                   | kg m-2 s-1 | fprec           | mass flux of frozen precip                    | |
!! mean_merid_moment_flx             | Pa         | v_flux          | j-directed wind stress into ocean             | [vector rotation] (@ref VectorRotations) applied - lat-lon to tripolar
!! mean_net_lw_flx                   | W m-2      | lw_flux         | long wave radiation                           | |
!! mean_net_sw_ir_dif_flx            | W m-2      | sw_flux_nir_dif | diffuse near IR shortwave radiation           | |
!! mean_net_sw_ir_dir_flx            | W m-2      | sw_flux_nir_dir | direct near IR shortwave radiation            | |
!! mean_net_sw_vis_dif_flx           | W m-2      | sw_flux_vis_dif | diffuse visible shortware radiation           | |
!! mean_net_sw_vis_dir_flx           | W m-2      | sw_flux_vis_dir | direct visible shortware radiation            | |
!! mean_prec_rate                    | kg m-2 s-1 | lprec           | mass flux of liquid precip                    | |
!! mean_runoff_heat_flx              | W m-2      | runoff_hflx     | heat flux, relative to 0C, of liquid land water into ocean | |
!! mean_runoff_rate                  | kg m-2 s-1 | runoff          | mass flux of liquid runoff                    | |
!! mean_salt_rate                    | kg m-2 s-1 | salt_flux       | salt flux                                     | |
!! mean_sensi_heat_flx               | W m-2      | t_flux          | sensible heat flux into ocean                 | sign reversed (- sensi)
!! mean_zonal_moment_flx             | Pa         | u_flux          | j-directed wind stress into ocean             | [vector rotation] (@ref VectorRotations) applied - lat-lon to tripolar
!!
!!
!! @subsection ExportField Export Fields
!!
!! Export fields are populated from the `ocean_sfc` parameter (type `ocean_public_type`)
!! after the call to `update_ocean_model()`.
!!
!! Standard Name                     | Units      | Model Variable  | Description                               | Notes
!! ----------------------------------|------------|-----------------|-------------------------------------------|---------------------------------------------------------------------
!! freezing_melting_potential        | W m-2      | frazil          | accumulated heating from frazil formation | cap converts model units (J m-2) to (W m-2) for export
!! ocean_mask                        |            |                 | ocean mask                                | |
!! ocn_current_merid                 | m s-1      | v_surf          | j-directed surface velocity on u-cell     | [vector rotation] (@ref VectorRotations) applied - tripolar to lat-lon
!! ocn_current_zonal                 | m s-1      | u_surf          | i-directed surface velocity on u-cell     | [vector rotation] (@ref VectorRotations) applied - tripolar to lat-lon
!! s_surf                            | psu        | s_surf          | sea surface salinity on t-cell            | |
!! sea_lev                           | m          | sea_lev         | sea level                                 | model computation is eta_t + patm/(rho0*grav) - eta_geoid - eta_tide
!! sea_surface_temperature           | K          | t_surf          | sea surface temperature on t-cell         | |
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
!! * `OceanSolo` - when set to "true", this option indicates that MOM is being run
!!   uncoupled; in this case the vector rotations and other data manipulations
!!   on import fields are skipped
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

  use ocean_model_mod,          only: ocean_model_restart, ocean_public_type, ocean_state_type
  use ocean_model_mod,          only: ocean_model_data_get
  use ocean_model_mod,          only: ocean_model_init , update_ocean_model, ocean_model_end, get_ocean_grid
#ifdef MOM6_CAP
  use ocean_model_mod,          only: ice_ocean_boundary_type
  use MOM_grid,                 only: ocean_grid_type
#else
  use ocean_types_mod,          only: ice_ocean_boundary_type, ocean_grid_type
#endif

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS      => SetServices, &
    model_label_Advance   => label_Advance, &
    model_label_Finalize  => label_Finalize

  use time_utils_mod

  implicit none
  private
  public SetServices

  type ocean_internalstate_type
    type(ocean_public_type),       pointer :: ocean_public_type_ptr
    type(ocean_state_type),        pointer :: ocean_state_type_ptr
    type(ice_ocean_boundary_type), pointer :: ice_ocean_boundary_type_ptr
    type(ocean_grid_type),         pointer :: ocean_grid_ptr
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

  integer,parameter :: fldsMax = 100
  integer :: fldsToOcn_num = 0
  type (fld_list_type) :: fldsToOcn(fldsMax)
  integer :: fldsFrOcn_num = 0
  type (fld_list_type) :: fldsFrOcn(fldsMax)

  integer   :: import_slice = 1
  integer   :: export_slice = 1
  character(len=256) :: tmpstr
  integer   :: dbrc

  type(ESMF_Grid), save   :: mom_grid_i
  logical                 :: write_diagnostics = .true.
  logical                 :: profile_memory = .true.
  logical                 :: ocean_solo = .true.
  logical                 :: grid_attach_area = .false.
  integer(ESMF_KIND_I8)   :: restart_interval

  contains
  !-----------------------------------------------------------------------
  !------------------- Solo Ocean code starts here -----------------------
  !-----------------------------------------------------------------------

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
    
    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
      specRoutine=ocean_model_finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine SetServices

  !-----------------------------------------------------------------------------

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

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_AttributeGet(gcomp, name="DumpFields", value=value, defaultValue="true", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write_diagnostics=(trim(value)=="true")
    call ESMF_LogWrite('MOM_CAP:DumpFields = '//trim(value), ESMF_LOGMSG_INFO, rc=dbrc)  

    call ESMF_AttributeGet(gcomp, name="ProfileMemory", value=value, defaultValue="true", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    profile_memory=(trim(value)/="false")
    call ESMF_LogWrite('MOM_CAP:ProfileMemory = '//trim(value), ESMF_LOGMSG_INFO, rc=dbrc)  

    call ESMF_AttributeGet(gcomp, name="OceanSolo", value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ocean_solo=(trim(value)=="true")
    call ESMF_LogWrite('MOM_CAP:OceanSolo = '//trim(value), ESMF_LOGMSG_INFO, rc=dbrc)  

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

    call ESMF_AttributeGet(gcomp, name="GridAttachArea", value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    grid_attach_area=(trim(value)=="true")
    call ESMF_LogWrite('MOM_CAP:GridAttachArea = '//trim(value), ESMF_LOGMSG_INFO, rc=dbrc)  
    
  end subroutine
  
  !-----------------------------------------------------------------------------

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
    
    type (ocean_public_type),      pointer :: Ocean_sfc   => NULL()
    type (ocean_state_type),       pointer :: Ocean_state => NULL()
    type(ice_ocean_boundary_type), pointer :: Ice_ocean_boundary => NULL()
    type(ocean_internalstate_wrapper)      :: ocean_internalstate

    type(time_type)                        :: Run_len      ! length of experiment 
    type(time_type)                        :: Time        
    type(time_type)                        :: Time_restart
    type(time_type)                        :: DT
    integer                                :: DT_OCEAN
    integer                                :: isc,iec,jsc,jec
    integer                                :: dt_cpld  = 86400
    integer                                :: year=0, month=0, day=0, hour=0, minute=0, second=0
    integer                                :: mpi_comm_mom

    type(ESMF_Grid)                        :: gridIn
    type(ESMF_Grid)                        :: gridOut

    integer                                :: npet, npet_x, npet_y
    character(len=*),parameter  :: subname='(mom_cap:InitializeAdvertise)'

    rc = ESMF_SUCCESS

    allocate(Ice_ocean_boundary)
    !allocate(Ocean_state) ! ocean_model_init allocate this pointer
    allocate(Ocean_sfc)
    allocate(ocean_internalstate%ptr)
    ocean_internalstate%ptr%ice_ocean_boundary_type_ptr => Ice_ocean_boundary
    ocean_internalstate%ptr%ocean_public_type_ptr       => Ocean_sfc
    ocean_internalstate%ptr%ocean_state_type_ptr        => Ocean_state

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

    call ESMF_TimeGet (MyTime,                    &
                       YY=YEAR, MM=MONTH, DD=DAY, &
                       H=HOUR,    M =MINUTE,    S =SECOND,  &
                                        RC=rc )
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
    call set_calendar_type (JULIAN                )
    call diag_manager_init
    ! this ocean connector will be driven at set interval
    dt_cpld = DT_OCEAN
    DT = set_time (DT_OCEAN, 0)         
    Time = set_date (YEAR,MONTH,DAY,HOUR,MINUTE,SECOND)

    Ocean_sfc%is_ocean_pe = .true.
    call ocean_model_init(Ocean_sfc, Ocean_state, Time, Time)
    call data_override_init(Ocean_domain_in = Ocean_sfc%domain)
    call mpp_get_compute_domain(Ocean_sfc%domain, isc, iec, jsc, jec)

    allocate ( Ice_ocean_boundary% u_flux (isc:iec,jsc:jec),          &
               Ice_ocean_boundary% v_flux (isc:iec,jsc:jec),          &
               Ice_ocean_boundary% t_flux (isc:iec,jsc:jec),          &
               Ice_ocean_boundary% q_flux (isc:iec,jsc:jec),          &
               Ice_ocean_boundary% salt_flux (isc:iec,jsc:jec),       &
               Ice_ocean_boundary% lw_flux (isc:iec,jsc:jec),         &
               Ice_ocean_boundary% sw_flux_vis_dir (isc:iec,jsc:jec), &
               Ice_ocean_boundary% sw_flux_vis_dif (isc:iec,jsc:jec), &
               Ice_ocean_boundary% sw_flux_nir_dir (isc:iec,jsc:jec), &
               Ice_ocean_boundary% sw_flux_nir_dif (isc:iec,jsc:jec), &
               Ice_ocean_boundary% lprec (isc:iec,jsc:jec),           &
               Ice_ocean_boundary% fprec (isc:iec,jsc:jec),           &
               Ice_ocean_boundary% runoff (isc:iec,jsc:jec),          &
               Ice_ocean_boundary% calving (isc:iec,jsc:jec),         &
               Ice_ocean_boundary% runoff_hflx (isc:iec,jsc:jec),     &
               Ice_ocean_boundary% calving_hflx (isc:iec,jsc:jec),    &
               Ice_ocean_boundary% mi (isc:iec,jsc:jec),              &
               Ice_ocean_boundary% p (isc:iec,jsc:jec))

    Ice_ocean_boundary%u_flux          = 0.0
    Ice_ocean_boundary%v_flux          = 0.0
    Ice_ocean_boundary%t_flux          = 0.0
    Ice_ocean_boundary%q_flux          = 0.0
    Ice_ocean_boundary%salt_flux       = 0.0
    Ice_ocean_boundary%lw_flux         = 0.0
    Ice_ocean_boundary%sw_flux_vis_dir = 0.0
    Ice_ocean_boundary%sw_flux_vis_dif = 0.0
    Ice_ocean_boundary%sw_flux_nir_dir = 0.0
    Ice_ocean_boundary%sw_flux_nir_dif = 0.0
    Ice_ocean_boundary%lprec           = 0.0
    Ice_ocean_boundary%fprec           = 0.0
    Ice_ocean_boundary%runoff          = 0.0
    Ice_ocean_boundary%calving         = 0.0
    Ice_ocean_boundary%runoff_hflx     = 0.0
    Ice_ocean_boundary%calving_hflx    = 0.0
    Ice_ocean_boundary%mi              = 0.0
    Ice_ocean_boundary%p               = 0.0

    call external_coupler_sbc_init(Ocean_sfc%domain, dt_cpld, Run_len)

    ocean_internalstate%ptr%ocean_state_type_ptr => Ocean_state
    call ESMF_GridCompSetInternalState(gcomp, ocean_internalstate, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call MOM_FieldsSetup(ice_ocean_boundary, ocean_sfc)

    call MOM_AdvertiseFields(importState, fldsToOcn_num, fldsToOcn, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call MOM_AdvertiseFields(exportState, fldsFrOcn_num, fldsFrOcn, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

#ifdef MOM6_CAP
    ! When running mom6 solo, the rotation angles are not computed internally
    ! in MOM6. We need to 
    ! calculate cos and sin of rotational angle for MOM6; the values
    ! are stored in ocean_internalstate%ptr%ocean_grid_ptr%cos_rot and sin_rot
    ! The rotation angles are retrieved during run time to rotate incoming
    ! and outgoing vectors
    !
    call calculate_rot_angle(Ocean_state, ocean_sfc, &
      ocean_internalstate%ptr%ocean_grid_ptr)
#endif

    write(*,*) '----- MOM initialization phase Advertise completed'

  end subroutine InitializeAdvertise
  
  !-----------------------------------------------------------------------------

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
    type(ESMF_VM)                          :: vm
    type(ESMF_Grid)                        :: gridIn
    type(ESMF_Grid)                        :: gridOut
    type(ESMF_DeLayout)                    :: delayout
    type(ESMF_Distgrid)                    :: Distgrid
    type(ESMF_DistGridConnection), allocatable :: connectionList(:)
    type (ocean_public_type),      pointer :: Ocean_sfc   => NULL()
    type (ocean_state_type),       pointer :: Ocean_state => NULL()
    type(ice_ocean_boundary_type), pointer :: Ice_ocean_boundary => NULL()
    type(ocean_internalstate_wrapper)      :: ocean_internalstate
    integer                                :: npet, ntiles
    integer                                :: nxg, nyg, cnt
    integer                                :: isc,iec,jsc,jec
    integer, allocatable                   :: xb(:),xe(:),yb(:),ye(:),pe(:)
    integer, allocatable                   :: deBlockList(:,:,:), &
                                              petMap(:),deLabelList(:), &
                                              indexList(:)
    integer                                :: ioff, joff
    integer                                :: i, j, n, i1, j1, n1, icount
    integer                                :: lbnd1,ubnd1,lbnd2,ubnd2
    integer                                :: lbnd3,ubnd3,lbnd4,ubnd4
    integer                                :: nblocks_tot
    logical                                :: found
    real(ESMF_KIND_R8), allocatable        :: ofld(:,:), gfld(:,:)
    real(ESMF_KIND_R8), pointer            :: t_surf(:,:)
    integer(ESMF_KIND_I4), pointer         :: dataPtr_mask(:,:)
    real(ESMF_KIND_R8), pointer            :: dataPtr_area(:,:)
    real(ESMF_KIND_R8), pointer            :: dataPtr_xcen(:,:)
    real(ESMF_KIND_R8), pointer            :: dataPtr_ycen(:,:)
    real(ESMF_KIND_R8), pointer            :: dataPtr_xcor(:,:)
    real(ESMF_KIND_R8), pointer            :: dataPtr_ycor(:,:)
    type(ESMF_Field)                       :: field_t_surf
    character(len=*),parameter  :: subname='(mom_cap:InitializeRealize)'
    
    rc = ESMF_SUCCESS

    call ESMF_GridCompGetInternalState(gcomp, ocean_internalstate, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    Ice_ocean_boundary => ocean_internalstate%ptr%ice_ocean_boundary_type_ptr
    Ocean_sfc          => ocean_internalstate%ptr%ocean_public_type_ptr
    Ocean_state        => ocean_internalstate%ptr%ocean_state_type_ptr

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm, petCount=npet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !---------------------------------
    ! global mom grid size
    !---------------------------------

    call mpp_get_global_domain(Ocean_sfc%domain, xsize=nxg, ysize=nyg)
    write(tmpstr,'(a,2i6)') subname//' nxg,nyg = ',nxg,nyg
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)  

    !---------------------------------
    ! number of tiles per PET, assumed to be 1, and number of pes (tiles) total
    !---------------------------------

    ntiles=mpp_get_ntile_count(Ocean_sfc%domain) ! this is tiles on this pe
    if (ntiles /= 1) then
      rc = ESMF_FAILURE
      call ESMF_LogWrite(subname//' ntiles must be 1', ESMF_LOGMSG_ERROR, rc=dbrc)  
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    ntiles=mpp_get_domain_npes(Ocean_sfc%domain)
    write(tmpstr,'(a,1i6)') subname//' ntiles = ',ntiles
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)  

    !---------------------------------
    ! get start and end indices of each tile and their PET
    !---------------------------------

    allocate(xb(ntiles),xe(ntiles),yb(ntiles),ye(ntiles),pe(ntiles))
    call mpp_get_compute_domains(Ocean_sfc%domain, xbegin=xb, xend=xe, ybegin=yb, yend=ye)
    call mpp_get_pelist(Ocean_sfc%domain, pe)
    do n = 1,ntiles
      write(tmpstr,'(a,6i6)') subname//' tiles ',n,pe(n),xb(n),xe(n),yb(n),ye(n)
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)  
    enddo

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

    call mpp_get_compute_domain(Ocean_sfc%domain, isc, iec, jsc, jec)

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

    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'mask', ofld, isc, jsc)
    write(tmpstr,*) subname//' ofld mask = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(Ocean_sfc%domain, ofld, gfld)
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
      call ocean_model_data_get(Ocean_state, Ocean_sfc, 'area', ofld, isc, jsc)
      write(tmpstr,*) subname//' ofld area = ',minval(ofld),maxval(ofld)
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
      call mpp_global_field(Ocean_sfc%domain, ofld, gfld)
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

    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'tlon', ofld, isc, jsc)
    write(tmpstr,*) subname//' ofld xt = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(Ocean_sfc%domain, ofld, gfld)
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

    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'tlat', ofld, isc, jsc)
    write(tmpstr,*) subname//' ofld yt = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(Ocean_sfc%domain, ofld, gfld)
    write(tmpstr,*) subname//' gfld yt = ',minval(gfld),maxval(gfld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    do j = lbnd2, ubnd2
    do i = lbnd1, ubnd1
       j1 = j - lbnd2 + jsc
       i1 = i - lbnd1 + isc
       dataPtr_ycen(i,j) = ofld(i1,j1)
    enddo
    enddo

#ifdef MOM5_CAP
    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'ulon', ofld, isc, jsc)
#endif

#ifdef MOM6_CAP
    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'geoLonBu', ofld, isc, jsc)
#endif
    write(tmpstr,*) subname//' ofld xu = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(Ocean_sfc%domain, ofld, gfld)
    write(tmpstr,*) subname//' gfld xu = ',minval(gfld),maxval(gfld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    do j = lbnd4, ubnd4
    do i = lbnd3, ubnd3
       j1 = j - lbnd4 + jsc - 1
       i1 = mod(i - lbnd3 + isc - 2 + nxg, nxg) + 1
       if (j1 == 0) then
          dataPtr_xcor(i,j) = 2*gfld(i1,1) - gfld(i1,2)
!          if (dataPtr_xcor(i,j)-dataPtr_xcen(i,j) > 180.) dataPtr_xcor(i,j) = dataPtr_xcor(i,j) - 360.
!          if (dataPtr_xcor(i,j)-dataPtr_xcen(i,j) < 180.) dataPtr_xcor(i,j) = dataPtr_xcor(i,j) + 360.
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

! The corner latitude values are treated differently because MOM5 runs on B-Grid while
! MOM6 runs on C-Grid.
#ifdef MOM5_CAP
    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'ulat', ofld, isc, jsc)
#endif

#ifdef MOM6_CAP
    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'geoLatBu', ofld, isc, jsc)
#endif

    write(tmpstr,*) subname//' ofld yu = ',minval(ofld),maxval(ofld)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call mpp_global_field(Ocean_sfc%domain, ofld, gfld)
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

      call ocean_model_data_get(Ocean_state, Ocean_sfc, 'mask', ofld, isc, jsc)

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

    call NUOPC_Write(exportState, fileNamePrefix='init_field_ocn_export_', &
      timeslice=1, relaxedFlag=.true., rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(*,*) '----- MOM initialization phase Realize completed'

  end subroutine InitializeRealize
  
  !> Called by NUOPC to advance the model a single timestep.
  !!  
  !! @param gcomp an ESMF_GridComp object
  !! @param rc return code
  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)                    :: gcomp
    integer, intent(out)                   :: rc
    
    ! local variables
    type(ESMF_Clock)                       :: clock
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Time)                        :: currTime
    type(ESMF_TimeInterval)                :: timeStep
    type(ESMF_Time)                        :: startTime
    type(ESMF_TimeInterval)                :: time_elapsed
    integer(ESMF_KIND_I8)                  :: n_interval, time_elapsed_sec
    character(len=64)                      :: timestamp

    type (ocean_public_type),      pointer :: Ocean_sfc          => NULL()
    type (ocean_state_type),       pointer :: Ocean_state        => NULL()
    type(ice_ocean_boundary_type), pointer :: Ice_ocean_boundary => NULL()
    type(ocean_internalstate_wrapper)      :: ocean_internalstate

    ! define some time types 
    type(time_type)                        :: Time        
    type(time_type)                        :: Time_step_coupled
    type(time_type)                        :: Time_restart_current

    integer :: dth, dtm, dts, dt_cpld  = 86400
    integer :: isc,iec,jsc,jec,lbnd1,ubnd1,lbnd2,ubnd2
    integer :: i,j,i1,j1
    real(ESMF_KIND_R8), allocatable        :: ofld(:,:), ocz(:,:), ocm(:,:)
    real(ESMF_KIND_R8), allocatable        :: mmmf(:,:), mzmf(:,:)
    integer :: nc
    real(ESMF_KIND_R8), pointer :: dataPtr_mask(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_mmmf(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_mzmf(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ocz(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ocm(:,:) 
    real(ESMF_KIND_R8), pointer :: dataPtr_frazil(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_evap(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_sensi(:,:)
    type(ocean_grid_type), pointer :: Ocean_grid
    character(240)              :: msgString
    character(len=*),parameter  :: subname='(mom_cap:ModelAdvance)'

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
    Ocean_sfc          => ocean_internalstate%ptr%ocean_public_type_ptr
    Ocean_state        => ocean_internalstate%ptr%ocean_state_type_ptr

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing OCN from: ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
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
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
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

    call external_coupler_sbc_before(Ice_ocean_boundary, Ocean_sfc, nc, dt_cpld )

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

    call mpp_get_compute_domain(Ocean_sfc%domain, isc, iec, jsc, jec)

   if(.not. ocean_solo) then
    call State_getFldPtr(exportState,'ocean_mask',dataPtr_mask,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    lbnd1 = lbound(dataPtr_mask,1)
    ubnd1 = ubound(dataPtr_mask,1)
    lbnd2 = lbound(dataPtr_mask,2)
    ubnd2 = ubound(dataPtr_mask,2)

#ifdef MOM5_CAP
    call get_ocean_grid(Ocean_grid)
#endif
#ifdef MOM6_CAP
    Ocean_grid => ocean_internalstate%ptr%ocean_grid_ptr
#endif

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
        mzmf(i,j) = Ocean_grid%cos_rot(i1,j1)*dataPtr_mzmf(i,j) &
                  + Ocean_grid%sin_rot(i1,j1)*dataPtr_mmmf(i,j)
        mmmf(i,j) = Ocean_grid%cos_rot(i1,j1)*dataPtr_mmmf(i,j) &
                  - Ocean_grid%sin_rot(i1,j1)*dataPtr_mzmf(i,j)
      enddo
    enddo
    dataPtr_mzmf = mzmf
    dataPtr_mmmf = mmmf
    deallocate(mzmf, mmmf)
   endif  ! not ocean_solo

    !Optionally write restart files when currTime-startTime is integer multiples of restart_interval
    if(restart_interval > 0 ) then
      time_elapsed = currTime - startTime
      call ESMF_TimeIntervalGet(time_elapsed, s_i8=time_elapsed_sec, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      n_interval = time_elapsed_sec / restart_interval
      if((n_interval .gt. 0) .and. (n_interval*restart_interval == time_elapsed_sec)) then
          time_restart_current = esmf2fms_time(currTime)
          timestamp = date_to_string(time_restart_current)
          call ESMF_LogWrite("MOM: Writing restart at "//trim(timestamp), ESMF_LOGMSG_INFO, rc=dbrc)
          write(*,*) 'calling ocean_model_restart'
          call ocean_model_restart(Ocean_state, timestamp)
      endif
    endif

    if(profile_memory) call ESMF_VMLogMemInfo("Entering MOM update_ocean_model: ")
    call update_ocean_model(Ice_ocean_boundary, Ocean_state, Ocean_sfc, Time, Time_step_coupled)
    if(profile_memory) call ESMF_VMLogMemInfo("Leaving MOM update_ocean_model: ")

   if(.not. ocean_solo) then
    allocate(ofld(isc:iec,jsc:jec))

    call ocean_model_data_get(Ocean_state, Ocean_sfc, 'mask', ofld, isc, jsc)
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
        dataPtr_ocz(i,j) = Ocean_grid%cos_rot(i1,j1)*ocz(i,j) &
                         - Ocean_grid%sin_rot(i1,j1)*ocm(i,j)
        dataPtr_ocm(i,j) = Ocean_grid%cos_rot(i1,j1)*ocm(i,j) &
                         + Ocean_grid%sin_rot(i1,j1)*ocz(i,j)
      enddo
    enddo
    deallocate(ocz, ocm)
   endif  ! not ocean_solo

    call ESMF_LogWrite("Before writing diagnostics", ESMF_LOGMSG_INFO, rc=rc)
    if(write_diagnostics) then
      call NUOPC_Write(exportState, fileNamePrefix='field_ocn_export_', &
        timeslice=export_slice, relaxedFlag=.true., rc=rc) 
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      export_slice = export_slice + 1
    endif

    call ESMF_LogWrite("Before calling sbc forcing", ESMF_LOGMSG_INFO, rc=rc)
    call external_coupler_sbc_after(Ice_ocean_boundary, Ocean_sfc, nc, dt_cpld )

    call ESMF_LogWrite("Before dumpMomInternal", ESMF_LOGMSG_INFO, rc=rc)
    !write(*,*) 'MOM: --- run phase called ---'
    call dumpMomInternal(mom_grid_i, import_slice, "mean_zonal_moment_flx", "will provide", Ice_ocean_boundary%u_flux)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_merid_moment_flx", "will provide", Ice_ocean_boundary%v_flux)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_sensi_heat_flx"  , "will provide", Ice_ocean_boundary%t_flux)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_evap_rate"       , "will provide", Ice_ocean_boundary%q_flux)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_salt_rate"       , "will provide", Ice_ocean_boundary%salt_flux)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_net_lw_flx"      , "will provide", Ice_ocean_boundary%lw_flux  )
    call dumpMomInternal(mom_grid_i, import_slice, "mean_net_sw_vis_dir_flx", "will provide", Ice_ocean_boundary%sw_flux_vis_dir)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_net_sw_vis_dif_flx", "will provide", Ice_ocean_boundary%sw_flux_vis_dif)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_net_sw_ir_dir_flx" , "will provide", Ice_ocean_boundary%sw_flux_nir_dir)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_net_sw_ir_dif_flx" , "will provide", Ice_ocean_boundary%sw_flux_nir_dif)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_prec_rate"       , "will provide", Ice_ocean_boundary%lprec  )
    call dumpMomInternal(mom_grid_i, import_slice, "mean_fprec_rate"      , "will provide", Ice_ocean_boundary%fprec  )
    call dumpMomInternal(mom_grid_i, import_slice, "mean_runoff_rate"     , "will provide", Ice_ocean_boundary%runoff )
    call dumpMomInternal(mom_grid_i, import_slice, "mean_calving_rate"    , "will provide", Ice_ocean_boundary%calving)
    call dumpMomInternal(mom_grid_i, import_slice, "mean_runoff_heat_flx" , "will provide", Ice_ocean_boundary%runoff_hflx )
    call dumpMomInternal(mom_grid_i, import_slice, "mean_calving_heat_flx", "will provide", Ice_ocean_boundary%calving_hflx)
    call dumpMomInternal(mom_grid_i, import_slice, "inst_pres_height_surface" , "will provide", Ice_ocean_boundary%p )
    call dumpMomInternal(mom_grid_i, import_slice, "mass_of_overlying_sea_ice", "will provide", Ice_ocean_boundary%mi)

!--------- export fields -------------

    call dumpMomInternal(mom_grid_i, export_slice, "ocean_mask", "will provide", dataPtr_mask)
    call dumpMomInternal(mom_grid_i, export_slice, "sea_surface_temperature", "will provide", Ocean_sfc%t_surf)
    call dumpMomInternal(mom_grid_i, export_slice, "s_surf"    , "will provide", Ocean_sfc%s_surf )
    call dumpMomInternal(mom_grid_i, export_slice, "ocn_current_zonal", "will provide", Ocean_sfc%u_surf )
    call dumpMomInternal(mom_grid_i, export_slice, "ocn_current_merid", "will provide", Ocean_sfc%v_surf )
    call dumpMomInternal(mom_grid_i, export_slice, "sea_lev"   , "will provide", Ocean_sfc%sea_lev)

    if(profile_memory) call ESMF_VMLogMemInfo("Leaving MOM Model_ADVANCE: ")
  end subroutine ModelAdvance

  !> Called by NUOPC at the end of the run to clean up.
  !!
  !! @param gcomp an ESMF_GridComp object
  !! @param rc return code
  subroutine ocean_model_finalize(gcomp, rc)

    ! input arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type (ocean_public_type),      pointer :: Ocean_sfc          
    type (ocean_state_type),       pointer :: Ocean_state
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

    Ocean_sfc          => ocean_internalstate%ptr%ocean_public_type_ptr
    Ocean_state        => ocean_internalstate%ptr%ocean_state_type_ptr

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

    call ocean_model_end (Ocean_sfc, Ocean_State, Time)
    call diag_manager_end(Time )
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

  subroutine external_coupler_sbc_before(Ice_ocean_boundary, Ocean_sfc, nsteps, dt_cpld )
  implicit none
  type (ice_ocean_boundary_type), intent(INOUT) :: Ice_ocean_boundary
  type (ocean_public_type) , intent(INOUT)        :: Ocean_sfc
  integer , intent(IN)                       :: nsteps, dt_cpld
  return
  end subroutine external_coupler_sbc_before


  subroutine external_coupler_sbc_after(Ice_ocean_boundary, Ocean_sfc, nsteps, dt_cpld )
  type (ice_ocean_boundary_type) :: Ice_ocean_boundary
  type (ocean_public_type)         :: Ocean_sfc
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
    type(ESMF_State), intent(in) :: ST
    character(len=*), intent(in) :: fldname
    real(ESMF_KIND_R8), pointer, intent(in) :: fldptr(:,:)
    integer, intent(out), optional :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
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
  subroutine MOM_AdvertiseFields(state, nfields, field_defs, rc)

    type(ESMF_State), intent(inout)             :: state
    integer,intent(in)                          :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    integer, intent(inout)                      :: rc

    integer                                     :: i
    character(len=*),parameter  :: subname='(mom_cap:MOM_AdvertiseFields)'

    rc = ESMF_SUCCESS

    do i = 1, nfields

      call NUOPC_Advertise(state, &
        standardName=field_defs(i)%stdname, &
        name=field_defs(i)%shortname, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    enddo

  end subroutine MOM_AdvertiseFields

  !-----------------------------------------------------------------------------

  subroutine MOM_RealizeFields(state, grid, nfields, field_defs, tag, rc)

    type(ESMF_State), intent(inout)             :: state
    type(ESMF_Grid), intent(in)                 :: grid
    integer, intent(in)                         :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    character(len=*), intent(in)                :: tag
    integer, intent(inout)                      :: rc

    integer                                     :: i
    type(ESMF_Field)                            :: field
    integer                                     :: npet, nx, ny, pet, elb(2), eub(2), clb(2), cub(2), tlb(2), tub(2)
    type(ESMF_VM)                               :: vm
    character(len=*),parameter  :: subname='(mom_cap:MOM_RealizeFields)'
 
    rc = ESMF_SUCCESS

    do i = 1, nfields

      if (field_defs(i)%assoc) then
        write(tmpstr, *) subname, tag, ' Field ', field_defs(i)%shortname, ':', &
          lbound(field_defs(i)%farrayPtr,1), ubound(field_defs(i)%farrayPtr,1), &
          lbound(field_defs(i)%farrayPtr,2), ubound(field_defs(i)%farrayPtr,2)
        call ESMF_LogWrite(tmpstr, ESMF_LOGMSG_INFO, rc=dbrc)
        field = ESMF_FieldCreate(grid=grid, &
          farray=field_defs(i)%farrayPtr, indexflag=ESMF_INDEX_DELOCAL, &
!          farray=field_defs(i)%farrayPtr, indexflag=ESMF_INDEX_GLOBAL, &
          name=field_defs(i)%shortname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      else
        field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, indexflag=ESMF_INDEX_DELOCAL, &
          name=field_defs(i)%shortname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (NUOPC_IsConnected(state, fieldName=field_defs(i)%shortname)) then
        call NUOPC_Realize(state, field=field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is connected.", &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=dbrc)
!        call ESMF_FieldPrint(field=field, rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out
      else
        call ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is not connected.", &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=dbrc)
        ! TODO: Initialize the value in the pointer to 0 after proper restart is setup
        !if(associated(field_defs(i)%farrayPtr) ) field_defs(i)%farrayPtr = 0.0
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

  subroutine MOM_FieldsSetup(ice_ocean_boundary,ocean_sfc)
    type(ice_ocean_boundary_type), intent(in)   :: Ice_ocean_boundary
    type(ocean_public_type), intent(in)         :: Ocean_sfc
    character(len=*),parameter  :: subname='(mom_cap:MOM_FieldsSetup)'

  !!! fld_list_add(num, fldlist, stdname, transferOffer, data(optional), shortname(optional))

!--------- import fields -------------

! tcraig, don't point directly into mom data YET (last field is optional in interface)
! instead, create space for the field when it's "realized".
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_zonal_moment_flx", "will provide", data=Ice_ocean_boundary%u_flux)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_merid_moment_flx", "will provide", data=Ice_ocean_boundary%v_flux)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_sensi_heat_flx"  , "will provide", data=Ice_ocean_boundary%t_flux)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_evap_rate"       , "will provide", data=Ice_ocean_boundary%q_flux)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_salt_rate"       , "will provide", data=Ice_ocean_boundary%salt_flux)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_lw_flx"      , "will provide", data=Ice_ocean_boundary%lw_flux  )
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_sw_vis_dir_flx", "will provide", data=Ice_ocean_boundary%sw_flux_vis_dir)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_sw_vis_dif_flx", "will provide", data=Ice_ocean_boundary%sw_flux_vis_dif)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_sw_ir_dir_flx" , "will provide", data=Ice_ocean_boundary%sw_flux_nir_dir)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_net_sw_ir_dif_flx" , "will provide", data=Ice_ocean_boundary%sw_flux_nir_dif)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_prec_rate"       , "will provide", data=Ice_ocean_boundary%lprec  )
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_fprec_rate"      , "will provide", data=Ice_ocean_boundary%fprec  )
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_runoff_rate"     , "will provide", data=Ice_ocean_boundary%runoff )
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_calving_rate"    , "will provide", data=Ice_ocean_boundary%calving)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_runoff_heat_flx" , "will provide", data=Ice_ocean_boundary%runoff_hflx )
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mean_calving_heat_flx", "will provide", data=Ice_ocean_boundary%calving_hflx)
    call fld_list_add(fldsToOcn_num, fldsToOcn, "inst_pres_height_surface" , "will provide", data=Ice_ocean_boundary%p )
    call fld_list_add(fldsToOcn_num, fldsToOcn, "mass_of_overlying_sea_ice", "will provide", data=Ice_ocean_boundary%mi)

!--------- export fields -------------

    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocean_mask", "will provide")
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "sea_surface_temperature", "will provide", data=Ocean_sfc%t_surf)
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "s_surf"    , "will provide", data=Ocean_sfc%s_surf )
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocn_current_zonal", "will provide", data=Ocean_sfc%u_surf )
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocn_current_merid", "will provide", data=Ocean_sfc%v_surf )
!    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocn_current_idir", "will provide")
!    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "ocn_current_jdir", "will provide")
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "sea_lev"   , "will provide", data=Ocean_sfc%sea_lev)
    call fld_list_add(fldsFrOcn_num, fldsFrOcn, "freezing_melting_potential"   , "will provide", data=Ocean_sfc%frazil)

  end subroutine MOM_FieldsSetup

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
      fldlist(num)%farrayPtr    => data
    else
      fldlist(num)%assoc        = .false.
    endif

  end subroutine fld_list_add

  subroutine dumpMomInternal(grid, slice, stdname, nop, farray)

    type(ESMF_Grid)          :: grid
    integer, intent(in)      :: slice
    character(len=*)         :: stdname
    character(len=*)         :: nop
    real(ESMF_KIND_R8), dimension(:,:), target :: farray

    type(ESMF_Field)         :: field
    real(ESMF_KIND_R8), dimension(:,:), pointer  :: f2d
    integer                  :: rc

#ifdef MOM6_CAP
    return
#endif

    if(.not. write_diagnostics) return ! nop in production mode
    if(ocean_solo) return ! do not dump internal fields in ocean solo mode

    field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, &
      indexflag=ESMF_INDEX_DELOCAL, &
      name=stdname, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_FieldGet(field, farrayPtr=f2d, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    f2d(:,:) = farray(:,:)

    call ESMF_FieldWrite(field, fileName='field_ocn_internal_'//trim(stdname)//'.nc', &
      timeslice=slice, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_FieldDestroy(field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine

#ifdef MOM6_CAP
  subroutine calculate_rot_angle(OS, OSFC, OG)
    type(ocean_state_type), intent(in)    :: OS
    type(ocean_public_type), intent(in)   :: OSFC
    type(ocean_grid_type), pointer        :: OG

    integer                               :: i,j,ishift,jshift,ilb,iub,jlb,jub
    real                                  :: angle, lon_scale
    type(ocean_grid_type), pointer        :: G

    call get_ocean_grid(OS, G)

    !print *, 'lbound: ', lbound(G%geoLatT), lbound(G%geoLonT), lbound(G%sin_rot)
    !print *, 'ubound: ', ubound(G%geoLatT), ubound(G%geoLonT), ubound(G%sin_rot)

    !print *, minval(G%geoLatT), maxval(G%geoLatT)
    !print *, minval(G%geoLonT), maxval(G%geoLonT)
    !print *, G%isc, G%jsc, G%iec, G%jec

    !
    ! The bounds isc:iec goes from 5-104, isc-ishift:iec-ishift goes from 1:100
    !
    call mpp_get_compute_domain(OSFC%Domain, ilb, iub, jlb, jub)
    ishift = ilb-G%isc
    jshift = jlb-G%jsc
    !print *, 'ilb, iub, jlb, jub', ilb, iub, jlb, jub, ishift, jshift
    !print *, 'sizes', iub-ilb, jub-jlb, G%iec-G%isc, G%jec-G%jsc
    allocate(OG)
    allocate(OG%sin_rot(ilb:iub, jlb:jub))
    allocate(OG%cos_rot(ilb:iub, jlb:jub))

    ! loop 5-104
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      lon_scale    = cos((G%geoLatBu(I-1,J-1) + G%geoLatBu(I,J-1  ) + &
                          G%geoLatBu(I-1,J) + G%geoLatBu(I,J)) * atan(1.0)/180)
      angle        = atan2((G%geoLonBu(I-1,J) + G%geoLonBu(I,J) - &
                            G%geoLonBu(I-1,J-1) - G%geoLonBu(I,J-1))*lon_scale, &
                            G%geoLatBu(I-1,J) + G%geoLatBu(I,J) - &
                            G%geoLatBu(I-1,J-1) - G%geoLatBu(I,J-1) )
      OG%sin_rot(i+ishift,j+jshift) = sin(angle) ! angle is the clockwise angle from lat/lon to ocean
      OG%cos_rot(i+ishift,j+jshift) = cos(angle) ! grid (e.g. angle of ocean "north" from true north)
    enddo ; enddo
    !print *, minval(OG%sin_rot), maxval(OG%sin_rot)
    !print *, minval(OG%cos_rot), maxval(OG%cos_rot)

  end subroutine
#endif


end module mom_cap_mod
