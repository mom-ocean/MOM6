!> Controls where open boundary conditions are applied
module MOM_open_boundary

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms,                 only : sum_across_PEs
use MOM_cpu_clock,            only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator,        only : diag_ctrl, time_type
use MOM_domains,              only : pass_var, pass_vector
use MOM_domains,              only : To_All, SCALAR_PAIR, CGRID_NE
use MOM_error_handler,        only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,          only : get_param, log_version, param_file_type, log_param
use MOM_grid,                 only : ocean_grid_type, hor_index_type
use MOM_dyn_horgrid,          only : dyn_horgrid_type
use MOM_io,                   only : EAST_FACE, NORTH_FACE
use MOM_io,                   only : slasher, read_data, field_size, SINGLE_FILE
use MOM_io,                   only : vardesc, query_vardesc, var_desc
use MOM_restart,              only : register_restart_field, query_initialized, MOM_restart_CS
use MOM_obsolete_params,      only : obsolete_logical, obsolete_int, obsolete_real, obsolete_char
use MOM_string_functions,     only : extract_word, remove_spaces
use MOM_tracer_registry,      only : tracer_type, tracer_registry_type, tracer_name_lookup
use time_interp_external_mod, only : init_external_field, time_interp_external
use time_interp_external_mod, only : time_interp_external_init
use MOM_remapping,            only : remappingSchemesDoc, remappingDefaultScheme, remapping_CS
use MOM_remapping,            only : initialize_remapping, remapping_core_h, end_remapping
use MOM_regridding,           only : regridding_CS
use MOM_unit_scaling,         only : unit_scale_type
use MOM_variables,            only : thermo_var_ptrs
use MOM_verticalGrid,         only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public open_boundary_apply_normal_flow
public open_boundary_config
public open_boundary_init
public open_boundary_query
public open_boundary_end
public open_boundary_impose_normal_slope
public open_boundary_impose_land_mask
public radiation_open_bdry_conds
public set_tracer_data
public update_OBC_segment_data
public open_boundary_test_extern_uv
public open_boundary_test_extern_h
public open_boundary_zero_normal_flow
public register_OBC, OBC_registry_init
public register_file_OBC, file_OBC_end
public segment_tracer_registry_init
public segment_tracer_registry_end
public register_segment_tracer
public register_temp_salt_segments
public fill_temp_salt_segments
public open_boundary_register_restarts
public update_segment_tracer_reservoirs

integer, parameter, public :: OBC_NONE = 0      !< Indicates the use of no open boundary
integer, parameter, public :: OBC_SIMPLE = 1    !< Indicates the use of a simple inflow open boundary
integer, parameter, public :: OBC_WALL = 2      !< Indicates the use of a closed wall
integer, parameter, public :: OBC_FLATHER =  3  !< Indicates the use of a Flather open boundary
integer, parameter, public :: OBC_RADIATION = 4 !< Indicates the use of a radiation open boundary
integer, parameter, public :: OBC_DIRECTION_N = 100 !< Indicates the boundary is an effective northern boundary
integer, parameter, public :: OBC_DIRECTION_S = 200 !< Indicates the boundary is an effective southern boundary
integer, parameter, public :: OBC_DIRECTION_E = 300 !< Indicates the boundary is an effective eastern boundary
integer, parameter, public :: OBC_DIRECTION_W = 400 !< Indicates the boundary is an effective western boundary
integer, parameter         :: MAX_OBC_FIELDS = 100  !< Maximum number of data fields needed for OBC segments

!> Open boundary segment data from files (mostly).
type, public :: OBC_segment_data_type
  integer :: fid                                !< handle from FMS associated with segment data on disk
  integer :: fid_dz                             !< handle from FMS associated with segment thicknesses on disk
  character(len=8)                :: name       !< a name identifier for the segment data
  real, pointer, dimension(:,:,:) :: buffer_src=>NULL() !< buffer for segment data located at cell faces
                                                !! and on the original vertical grid
  integer                         :: nk_src     !< Number of vertical levels in the source data
  real, dimension(:,:,:), pointer :: dz_src=>NULL()     !< vertical grid cell spacing of the incoming segment data [m]
  real, dimension(:,:,:), pointer :: buffer_dst=>NULL() !< buffer src data remapped to the target vertical grid
  real, dimension(:,:), pointer   :: bt_vel=>NULL()     !< barotropic velocity [L T-1 ~> m s-1]
  real                            :: value              !< constant value if fid is equal to -1
end type OBC_segment_data_type

!> Tracer on OBC segment data structure, for putting into a segment tracer registry.
type, public :: OBC_segment_tracer_type
  real, dimension(:,:,:), pointer :: t          => NULL()  !< tracer concentration array
  real                            :: OBC_inflow_conc = 0.0 !< tracer concentration for generic inflows
  character(len=32)               :: name                  !< tracer name used for error messages
  type(tracer_type), pointer      :: Tr         => NULL()  !< metadata describing the tracer
  real, dimension(:,:,:), pointer :: tres       => NULL()  !< tracer reservoir array
  logical                         :: is_initialized        !< reservoir values have been set when True
end type OBC_segment_tracer_type

!> Registry type for tracers on segments
type, public :: segment_tracer_registry_type
  integer                       :: ntseg = 0         !< number of registered tracer segments
  type(OBC_segment_tracer_type) :: Tr(MAX_FIELDS_)   !< array of registered tracers
  logical                       :: locked = .false.  !< New tracers may be registered if locked=.false.
                                                     !! When locked=.true.,no more tracers can be registered.
                                                     !! Not sure who should lock it or when...
end type segment_tracer_registry_type

!> Open boundary segment data structure.
type, public :: OBC_segment_type
  logical :: Flather        !< If true, applies Flather + Chapman radiation of barotropic gravity waves.
  logical :: radiation      !< If true, 1D Orlanksi radiation boundary conditions are applied.
                            !! If False, a gradient condition is applied.
  logical :: radiation_tan  !< If true, 1D Orlanksi radiation boundary conditions are applied to
                            !! tangential flows.
  logical :: radiation_grad !< If true, 1D Orlanksi radiation boundary conditions are applied to
                            !! dudv and dvdx.
  logical :: oblique        !< Oblique waves supported at radiation boundary.
  logical :: oblique_tan    !< If true, 2D radiation boundary conditions are applied to
                            !! tangential flows.
  logical :: oblique_grad   !< If true, 2D radiation boundary conditions are applied to
                            !! dudv and dvdx.
  logical :: nudged         !< Optional supplement to radiation boundary.
  logical :: nudged_tan     !< Optional supplement to nudge tangential velocity.
  logical :: nudged_grad    !< Optional supplement to nudge normal gradient of tangential velocity.
  logical :: specified      !< Boundary normal velocity fixed to external value.
  logical :: specified_tan  !< Boundary tangential velocity fixed to external value.
  logical :: specified_grad !< Boundary gradient of tangential velocity fixed to external value.
  logical :: open           !< Boundary is open for continuity solver.
  logical :: gradient       !< Zero gradient at boundary.
  logical :: values_needed  !< Whether or not any external OBC fields are needed.
  logical :: u_values_needed!< Whether or not external u OBC fields are needed.
  logical :: v_values_needed!< Whether or not external v OBC fields are needed.
  logical :: t_values_needed!< Whether or not external T OBC fields are needed.
  logical :: s_values_needed!< Whether or not external S OBC fields are needed.
  logical :: z_values_needed!< Whether or not external zeta OBC fields are needed.
  logical :: g_values_needed!< Whether or not external gradient OBC fields are needed.
  integer :: direction      !< Boundary faces one of the four directions.
  logical :: is_N_or_S      !< True if the OB is facing North or South and exists on this PE.
  logical :: is_E_or_W      !< True if the OB is facing East or West and exists on this PE.
  logical :: is_E_or_W_2    !< True if the OB is facing East or West anywhere.
  type(OBC_segment_data_type), pointer, dimension(:) :: field=>NULL()   !<  OBC data
  integer :: num_fields     !< number of OBC data fields (e.g. u_normal,u_parallel and eta for Flather)
  character(len=32), pointer, dimension(:) :: field_names=>NULL() !< field names for this segment
  integer :: Is_obc         !< i-indices of boundary segment.
  integer :: Ie_obc         !< i-indices of boundary segment.
  integer :: Js_obc         !< j-indices of boundary segment.
  integer :: Je_obc         !< j-indices of boundary segment.
  real :: Velocity_nudging_timescale_in  !< Nudging timescale on inflow [T ~> s].
  real :: Velocity_nudging_timescale_out !< Nudging timescale on outflow [T ~> s].
  logical :: on_pe          !< true if segment is located in the computational domain
  logical :: temp_segment_data_exists !< true if temperature data arrays are present
  logical :: salt_segment_data_exists !< true if salinity data arrays are present
  real, pointer, dimension(:,:)   :: Cg=>NULL()     !< The external gravity wave speed [L T-1 ~> m s-1]
                                                    !! at OBC-points.
  real, pointer, dimension(:,:)   :: Htot=>NULL()   !< The total column thickness [H ~> m or kg m-2] at OBC-points.
  real, pointer, dimension(:,:,:) :: h=>NULL()      !< The cell thickness [H ~> m or kg m-2] at OBC-points.
  real, pointer, dimension(:,:,:) :: normal_vel=>NULL()     !< The layer velocity normal to the OB
                                                            !! segment [L T-1 ~> m s-1].
  real, pointer, dimension(:,:,:) :: tangential_vel=>NULL() !< The layer velocity tangential to the
                                                            !! OB segment [L T-1 ~> m s-1].
  real, pointer, dimension(:,:,:) :: tangential_grad=>NULL() !< The gradient of the velocity tangential
                                                            !! to the OB segment [T-1 ~> s-1].
  real, pointer, dimension(:,:,:) :: normal_trans=>NULL()   !< The layer transport normal to the OB
                                                            !! segment [H L2 T-1 ~> m3 s-1].
  real, pointer, dimension(:,:)   :: normal_vel_bt=>NULL()  !< The barotropic velocity normal to
                                                            !! the OB segment [L T-1 ~> m s-1].
  real, pointer, dimension(:,:)   :: eta=>NULL()            !< The sea-surface elevation along the segment [m].
  real, pointer, dimension(:,:,:) :: grad_normal=>NULL()    !< The gradient of the normal flow along the
                                                            !! segment times the grid spacing [L T-1 ~> m s-1]
  real, pointer, dimension(:,:,:) :: grad_tan=>NULL()       !< The gradient of the tangential flow along the
                                                            !! segment times the grid spacing [L T-1 ~> m s-1]
  real, pointer, dimension(:,:,:) :: grad_gradient=>NULL()  !< The gradient of the gradient of tangential flow along
                                                            !! the segment times the grid spacing [T-1 ~> s-1]
  real, pointer, dimension(:,:,:) :: rx_norm_rad=>NULL()    !< The previous normal phase speed use for EW radiation
                                                            !! OBC, in grid points per timestep [nondim]
  real, pointer, dimension(:,:,:) :: ry_norm_rad=>NULL()    !< The previous normal phase speed use for NS radiation
                                                            !! OBC, in grid points per timestep [nondim]
  real, pointer, dimension(:,:,:) :: rx_norm_obl=>NULL()    !< The previous normal radiation coefficient for EW
                                                            !! oblique OBCs [L2 T-2 ~> m2 s-2]
  real, pointer, dimension(:,:,:) :: ry_norm_obl=>NULL()    !< The previous normal radiation coefficient for NS
                                                            !! oblique OBCs [L2 T-2 ~> m2 s-2]
  real, pointer, dimension(:,:,:) :: cff_normal=>NULL()     !< The denominator for oblique radiation
                                                            !! for normal velocity [L2 T-2 ~> m2 s-2]
  real, pointer, dimension(:,:,:) :: nudged_normal_vel=>NULL() !< The layer velocity normal to the OB segment
                                                            !! that values should be nudged towards [L T-1 ~> m s-1].
  real, pointer, dimension(:,:,:) :: nudged_tangential_vel=>NULL() !< The layer velocity tangential to the OB segment
                                                            !! that values should be nudged towards [L T-1 ~> m s-1].
  real, pointer, dimension(:,:,:) :: nudged_tangential_grad=>NULL() !< The layer dvdx or dudy towards which nudging
                                                            !! can occur [T-1 ~> s-1].
  type(segment_tracer_registry_type), pointer  :: tr_Reg=> NULL()!< A pointer to the tracer registry for the segment.
  type(hor_index_type) :: HI !< Horizontal index ranges
  real :: Tr_InvLscale_out                                  !< An effective inverse length scale for restoring
                                                            !! the tracer concentration in a ficticious
                                                            !! reservior towards interior values when flow
                                                            !! is exiting the domain [L-1 ~> m-1]
  real :: Tr_InvLscale_in                                   !< An effective inverse length scale for restoring
                                                            !! the tracer concentration towards an externally
                                                            !! imposed value when flow is entering [L-1 ~> m-1]
end type OBC_segment_type

!> Open-boundary data
type, public :: ocean_OBC_type
  integer :: number_of_segments = 0                   !< The number of open-boundary segments.
  integer :: ke = 0                                   !< The number of model layers
  logical :: open_u_BCs_exist_globally = .false.      !< True if any zonal velocity points
                                                      !! in the global domain use open BCs.
  logical :: open_v_BCs_exist_globally = .false.      !< True if any meridional velocity points
                                                      !! in the global domain use open BCs.
  logical :: Flather_u_BCs_exist_globally = .false.   !< True if any zonal velocity points
                                                      !! in the global domain use Flather BCs.
  logical :: Flather_v_BCs_exist_globally = .false.   !< True if any meridional velocity points
                                                      !! in the global domain use Flather BCs.
  logical :: oblique_BCs_exist_globally = .false.     !< True if any velocity points
                                                      !! in the global domain use oblique BCs.
  logical :: nudged_u_BCs_exist_globally = .false.    !< True if any velocity points in the
                                                      !! global domain use nudged BCs.
  logical :: nudged_v_BCs_exist_globally = .false.    !< True if any velocity points in the
                                                      !! global domain use nudged BCs.
  logical :: specified_u_BCs_exist_globally = .false. !< True if any zonal velocity points
                                                      !! in the global domain use specified BCs.
  logical :: specified_v_BCs_exist_globally = .false. !< True if any meridional velocity points
                                                      !! in the global domain use specified BCs.
  logical :: radiation_BCs_exist_globally = .false.   !< True if radiations BCs are in use anywhere.
  logical :: user_BCs_set_globally = .false.          !< True if any OBC_USER_CONFIG is set
                                                      !! for input from user directory.
  logical :: update_OBC = .false.                     !< Is OBC data time-dependent
  logical :: needs_IO_for_data = .false.              !< Is any i/o needed for OBCs
  logical :: zero_vorticity = .false.                 !< If True, sets relative vorticity to zero on open boundaries.
  logical :: freeslip_vorticity = .false.             !< If True, sets normal gradient of tangential velocity to zero
                                                      !! in the relative vorticity on open boundaries.
  logical :: computed_vorticity = .false.             !< If True, uses external data for tangential velocity
                                                      !! in the relative vorticity on open boundaries.
  logical :: specified_vorticity = .false.            !< If True, uses external data for tangential velocity
                                                      !! gradients in the relative vorticity on open boundaries.
  logical :: zero_strain = .false.                    !< If True, sets strain to zero on open boundaries.
  logical :: freeslip_strain = .false.                !< If True, sets normal gradient of tangential velocity to zero
                                                      !! in the strain on open boundaries.
  logical :: computed_strain = .false.                !< If True, uses external data for tangential velocity to compute
                                                      !! normal gradient in the strain on open boundaries.
  logical :: specified_strain = .false.               !< If True, uses external data for tangential velocity gradients
                                                      !! to compute strain on open boundaries.
  logical :: zero_biharmonic = .false.                !< If True, zeros the Laplacian of flow on open boundaries for
                                                      !! use in the biharmonic viscosity term.
  logical :: brushcutter_mode = .false.               !< If True, read data on supergrid.
  real :: g_Earth                                     !< The gravitational acceleration [m s-2].
  logical, pointer, dimension(:) :: &
                   tracer_x_reservoirs_used => NULL() !< Dimensioned by the number of tracers, set globally,
                                                      !! true for those with x reservoirs (needed for restarts).
  logical, pointer, dimension(:) :: &
                   tracer_y_reservoirs_used => NULL() !< Dimensioned by the number of tracers, set globally,
                                                      !! true for those with y reservoirs (needed for restarts).
  integer                       :: ntr = 0            !< number of tracers
  ! Properties of the segments used.
  type(OBC_segment_type), pointer, dimension(:) :: &
    segment => NULL()   !< List of segment objects.
  ! Which segment object describes the current point.
  integer, pointer, dimension(:,:) :: &
    segnum_u => NULL(), &   !< Segment number of u-points.
    segnum_v => NULL()      !< Segment number of v-points.

  ! The following parameters are used in the baroclinic radiation code:
  real :: gamma_uv !< The relative weighting for the baroclinic radiation
                   !! velocities (or speed of characteristics) at the
                   !! new time level (1) or the running mean (0) for velocities.
                   !! Valid values range from 0 to 1, with a default of 0.3.
  real :: rx_max   !< The maximum magnitude of the baroclinic radiation
                   !! velocity (or speed of characteristics) [m s-1].  The
                   !! default value is 10 m s-1.
  logical :: OBC_pe !< Is there an open boundary on this tile?
  type(remapping_CS),      pointer :: remap_CS   !< ALE remapping control structure for segments only
  type(OBC_registry_type), pointer :: OBC_Reg => NULL()  !< Registry type for boundaries
  real, pointer, dimension(:,:,:) :: &
    rx_normal => NULL(), & !< Array storage for normal phase speed for EW radiation OBCs in units of
                           !! grid points per timestep [nondim]
    ry_normal => NULL(), & !< Array storage for normal phase speed for NS radiation OBCs in units of
                           !! grid points per timestep [nondim]
    rx_oblique => NULL(), & !< Array storage for oblique boundary condition restarts [L2 T-2 ~> m2 s-2]
    ry_oblique => NULL(), & !< Array storage for oblique boundary condition restarts [L2 T-2 ~> m2 s-2]
    cff_normal => NULL()   !< Array storage for oblique boundary condition restarts [L2 T-2 ~> m2 s-2]
  real, pointer, dimension(:,:,:,:) :: &
    tres_x => NULL(), & !< Array storage of tracer reservoirs for restarts [conc L ~> conc m]
    tres_y => NULL()    !< Array storage of tracer reservoirs for restarts [conc L ~> conc m]
  real :: silly_h  !< A silly value of thickness outside of the domain that can be used to test
                   !! the independence of the OBCs to this external data [H ~> m or kg m-2].
  real :: silly_u  !< A silly value of velocity outside of the domain that can be used to test
                   !! the independence of the OBCs to this external data [L T-1 ~> m s-1].
end type ocean_OBC_type

!> Control structure for open boundaries that read from files.
!! Probably lots to update here.
type, public :: file_OBC_CS ; private
  real :: tide_flow = 3.0e6         !< Placeholder for now...
end type file_OBC_CS

!> Type to carry something (what??) for the OBC registry.
type, public :: OBC_struct_type
  character(len=32)               :: name             !< OBC name used for error messages
end type OBC_struct_type

!> Type to carry basic OBC information needed for updating values.
type, public :: OBC_registry_type
  integer               :: nobc = 0          !< number of registered open boundary types.
  type(OBC_struct_type) :: OB(MAX_FIELDS_)   !< array of registered boundary types.
  logical               :: locked = .false.  !< New OBC types may be registered if locked=.false.
                                             !! When locked=.true.,no more boundaries can be registered.
end type OBC_registry_type

integer :: id_clock_pass !< A CPU time clock

character(len=40)  :: mdl = "MOM_open_boundary" !< This module's name.
! This include declares and sets the variable "version".
#include "version_variable.h"

contains

!> Enables OBC module and reads configuration parameters
!> This routine is called from MOM_initialize_fixed which
!> occurs before the initialization of the vertical coordinate
!> and ALE_init.  Therefore segment data are not fully initialized
!> here. The remainder of the segment data are initialized in a
!> later call to update_open_boundary_data

subroutine open_boundary_config(G, US, param_file, OBC)
  type(dyn_horgrid_type),  intent(inout) :: G   !< Ocean grid structure
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file handle
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundary control structure
  ! Local variables
  integer :: l ! For looping over segments
  logical :: debug_OBC, debug, mask_outside, reentrant_x, reentrant_y
  character(len=15) :: segment_param_str ! The run-time parameter name for each segment
  character(len=100) :: segment_str      ! The contents (rhs) for parameter "segment_param_str"
  character(len=200) :: config1          ! String for OBC_USER_CONFIG
  real               :: Lscale_in, Lscale_out ! parameters controlling tracer values at the boundaries [L ~> m]
  allocate(OBC)

  call log_version(param_file, mdl, version, &
                 "Controls where open boundaries are located, what kind of boundary condition "//&
                 "to impose, and what data to apply, if any.")
  call get_param(param_file, mdl, "OBC_NUMBER_OF_SEGMENTS", OBC%number_of_segments, &
                 "The number of open boundary segments.", &
                 default=0)
  call get_param(param_file, mdl, "G_EARTH", OBC%g_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80)
  call get_param(param_file, mdl, "OBC_USER_CONFIG", config1, &
                 "A string that sets how the open boundary conditions are "//&
                 " configured: \n", default="none", do_not_log=.true.)
  call get_param(param_file, mdl, "NK", OBC%ke, &
                 "The number of model layers", default=0, do_not_log=.true.)

  if (config1 /= "none" .and. config1 /= "dyed_obcs") OBC%user_BCs_set_globally = .true.

  if (OBC%number_of_segments > 0) then
    call get_param(param_file, mdl, "OBC_ZERO_VORTICITY", OBC%zero_vorticity, &
         "If true, sets relative vorticity to zero on open boundaries.", &
         default=.false.)
    call get_param(param_file, mdl, "OBC_FREESLIP_VORTICITY", OBC%freeslip_vorticity, &
         "If true, sets the normal gradient of tangential velocity to "//&
         "zero in the relative vorticity on open boundaries. This cannot "//&
         "be true if another OBC_XXX_VORTICITY option is True.", default=.true.)
    call get_param(param_file, mdl, "OBC_COMPUTED_VORTICITY", OBC%computed_vorticity, &
         "If true, uses the external values of tangential velocity "//&
         "in the relative vorticity on open boundaries. This cannot "//&
         "be true if another OBC_XXX_VORTICITY option is True.", default=.false.)
    call get_param(param_file, mdl, "OBC_SPECIFIED_VORTICITY", OBC%specified_vorticity, &
         "If true, uses the external values of tangential velocity "//&
         "in the relative vorticity on open boundaries. This cannot "//&
         "be true if another OBC_XXX_VORTICITY option is True.", default=.false.)
    if ((OBC%zero_vorticity .and. OBC%freeslip_vorticity) .or.  &
        (OBC%zero_vorticity .and. OBC%computed_vorticity) .or.  &
        (OBC%zero_vorticity .and. OBC%specified_vorticity) .or.  &
        (OBC%freeslip_vorticity .and. OBC%computed_vorticity) .or.  &
        (OBC%freeslip_vorticity .and. OBC%specified_vorticity) .or.  &
        (OBC%computed_vorticity .and. OBC%specified_vorticity))  &
         call MOM_error(FATAL, "MOM_open_boundary.F90, open_boundary_config:\n"//&
         "Only one of OBC_ZERO_VORTICITY, OBC_FREESLIP_VORTICITY, OBC_COMPUTED_VORTICITY\n"//&
         "and OBC_IMPORTED_VORTICITY can be True at once.")
    call get_param(param_file, mdl, "OBC_ZERO_STRAIN", OBC%zero_strain, &
         "If true, sets the strain used in the stress tensor to zero on open boundaries.", &
         default=.false.)
    call get_param(param_file, mdl, "OBC_FREESLIP_STRAIN", OBC%freeslip_strain, &
         "If true, sets the normal gradient of tangential velocity to "//&
         "zero in the strain use in the stress tensor on open boundaries. This cannot "//&
         "be true if another OBC_XXX_STRAIN option is True.", default=.true.)
    call get_param(param_file, mdl, "OBC_COMPUTED_STRAIN", OBC%computed_strain, &
         "If true, sets the normal gradient of tangential velocity to "//&
         "zero in the strain use in the stress tensor on open boundaries. This cannot "//&
         "be true if another OBC_XXX_STRAIN option is True.", default=.false.)
    call get_param(param_file, mdl, "OBC_SPECIFIED_STRAIN", OBC%specified_strain, &
         "If true, sets the normal gradient of tangential velocity to "//&
         "zero in the strain use in the stress tensor on open boundaries. This cannot "//&
         "be true if another OBC_XXX_STRAIN option is True.", default=.false.)
    if ((OBC%zero_strain .and. OBC%freeslip_strain) .or.  &
        (OBC%zero_strain .and. OBC%computed_strain) .or.  &
        (OBC%zero_strain .and. OBC%specified_strain) .or.  &
        (OBC%freeslip_strain .and. OBC%computed_strain) .or.  &
        (OBC%freeslip_strain .and. OBC%specified_strain) .or.  &
        (OBC%computed_strain .and. OBC%specified_strain))  &
         call MOM_error(FATAL, "MOM_open_boundary.F90, open_boundary_config: \n"//&
         "Only one of OBC_ZERO_STRAIN, OBC_FREESLIP_STRAIN, OBC_COMPUTED_STRAIN \n"//&
         "and OBC_IMPORTED_STRAIN can be True at once.")
    call get_param(param_file, mdl, "OBC_ZERO_BIHARMONIC", OBC%zero_biharmonic, &
         "If true, zeros the Laplacian of flow on open boundaries in the biharmonic "//&
         "viscosity term.", default=.false.)
    call get_param(param_file, mdl, "MASK_OUTSIDE_OBCS", mask_outside, &
         "If true, set the areas outside open boundaries to be land.", &
         default=.false.)
    call get_param(param_file, mdl, "DEBUG", debug, default=.false.)
    call get_param(param_file, mdl, "DEBUG_OBC", debug_OBC, default=.false.)
    if (debug_OBC .or. debug) &
      call log_param(param_file, mdl, "DEBUG_OBC", debug_OBC, &
                 "If true, do additional calls to help debug the performance "//&
                 "of the open boundary condition code.", default=.false., &
                 debuggingParam=.true.)

    call get_param(param_file, mdl, "OBC_SILLY_THICK", OBC%silly_h, &
                 "A silly value of thicknesses used outside of open boundary "//&
                 "conditions for debugging.", units="m", default=0.0, scale=US%m_to_Z, &
                 do_not_log=.not.debug_OBC, debuggingParam=.true.)
    call get_param(param_file, mdl, "OBC_SILLY_VEL", OBC%silly_u, &
                 "A silly value of velocities used outside of open boundary "//&
                 "conditions for debugging.", units="m/s", default=0.0, scale=US%m_s_to_L_T, &
                 do_not_log=.not.debug_OBC, debuggingParam=.true.)
    reentrant_x = .false.
    call get_param(param_file, mdl, "REENTRANT_X", reentrant_x, default=.true.)
    reentrant_y = .false.
    call get_param(param_file, mdl, "REENTRANT_Y", reentrant_y, default=.false.)

    ! Allocate everything
    ! Note the 0-segment is needed when %segnum_u/v(:,:) = 0
    allocate(OBC%segment(0:OBC%number_of_segments))
    do l=0,OBC%number_of_segments
      OBC%segment(l)%Flather = .false.
      OBC%segment(l)%radiation = .false.
      OBC%segment(l)%radiation_tan = .false.
      OBC%segment(l)%radiation_grad = .false.
      OBC%segment(l)%oblique = .false.
      OBC%segment(l)%oblique_tan = .false.
      OBC%segment(l)%oblique_grad = .false.
      OBC%segment(l)%nudged = .false.
      OBC%segment(l)%nudged_tan = .false.
      OBC%segment(l)%nudged_grad = .false.
      OBC%segment(l)%specified = .false.
      OBC%segment(l)%specified_tan = .false.
      OBC%segment(l)%specified_grad = .false.
      OBC%segment(l)%open = .false.
      OBC%segment(l)%gradient = .false.
      OBC%segment(l)%values_needed = .false.
      OBC%segment(l)%u_values_needed = .false.
      OBC%segment(l)%v_values_needed = .false.
      OBC%segment(l)%t_values_needed = .false.
      OBC%segment(l)%s_values_needed = .false.
      OBC%segment(l)%z_values_needed = .false.
      OBC%segment(l)%g_values_needed = .false.
      OBC%segment(l)%direction = OBC_NONE
      OBC%segment(l)%is_N_or_S = .false.
      OBC%segment(l)%is_E_or_W = .false.
      OBC%segment(l)%is_E_or_W_2 = .false.
      OBC%segment(l)%Velocity_nudging_timescale_in = 0.0
      OBC%segment(l)%Velocity_nudging_timescale_out = 0.0
      OBC%segment(l)%num_fields = 0
    enddo
    allocate(OBC%segnum_u(G%IsdB:G%IedB,G%jsd:G%jed)) ; OBC%segnum_u(:,:) = OBC_NONE
    allocate(OBC%segnum_v(G%isd:G%ied,G%JsdB:G%JedB)) ; OBC%segnum_v(:,:) = OBC_NONE

    do l = 1, OBC%number_of_segments
      write(segment_param_str(1:15),"('OBC_SEGMENT_',i3.3)") l
      call get_param(param_file, mdl, segment_param_str, segment_str, &
           "Documentation needs to be dynamic?????", &
           fail_if_missing=.true.)
      segment_str = remove_spaces(segment_str)
      if (segment_str(1:2) == 'I=') then
        call setup_u_point_obc(OBC, G, US, segment_str, l, param_file, reentrant_y)
      elseif (segment_str(1:2) == 'J=') then
        call setup_v_point_obc(OBC, G, US, segment_str, l, param_file, reentrant_x)
      else
        call MOM_error(FATAL, "MOM_open_boundary.F90, open_boundary_config: "//&
             "Unable to interpret "//segment_param_str//" = "//trim(segment_str))
      endif
    enddo

    !    if (open_boundary_query(OBC, needs_ext_seg_data=.true.)) &
    call initialize_segment_data(G, OBC, param_file)

    if (open_boundary_query(OBC, apply_open_OBC=.true.)) then
      call get_param(param_file, mdl, "OBC_RADIATION_MAX", OBC%rx_max, &
                   "The maximum magnitude of the baroclinic radiation "//&
                   "velocity (or speed of characteristics).  This is only "//&
                   "used if one of the open boundary segments is using Orlanski.", &
                   units="m s-1", default=10.0)
      call get_param(param_file, mdl, "OBC_RAD_VEL_WT", OBC%gamma_uv, &
                   "The relative weighting for the baroclinic radiation "//&
                   "velocities (or speed of characteristics) at the new "//&
                   "time level (1) or the running mean (0) for velocities. "//&
                   "Valid values range from 0 to 1. This is only used if "//&
                   "one of the open boundary segments is using Orlanski.", &
                   units="nondim", default=0.3)
    endif

    Lscale_in = 0.
    Lscale_out = 0.
    if (open_boundary_query(OBC, apply_open_OBC=.true.)) then
      call get_param(param_file, mdl, "OBC_TRACER_RESERVOIR_LENGTH_SCALE_OUT ", Lscale_out, &
                 "An effective length scale for restoring the tracer concentration "//&
                 "at the boundaries to externally imposed values when the flow "//&
                 "is exiting the domain.", units="m", default=0.0, scale=US%m_to_L)

      call get_param(param_file, mdl, "OBC_TRACER_RESERVOIR_LENGTH_SCALE_IN ", Lscale_in, &
                 "An effective length scale for restoring the tracer concentration "//&
                 "at the boundaries to values from the interior when the flow "//&
                 "is entering the domain.", units="m", default=0.0, scale=US%m_to_L)
    endif

    if (mask_outside) call mask_outside_OBCs(G, US, param_file, OBC)

    ! All tracers are using the same restoring length scale for now, but we may want to make this
    ! tracer-specific in the future for example, in cases where certain tracers are poorly constrained
    ! by data while others are well constrained - MJH.
    do l = 1, OBC%number_of_segments
      OBC%segment(l)%Tr_InvLscale_in = 0.0
      if (Lscale_in>0.) OBC%segment(l)%Tr_InvLscale_in =  1.0/Lscale_in
      OBC%segment(l)%Tr_InvLscale_out = 0.0
      if (Lscale_out>0.) OBC%segment(l)%Tr_InvLscale_out =  1.0/Lscale_out
    enddo

  endif ! OBC%number_of_segments > 0

    ! Safety check
  if ((OBC%open_u_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally) .and. &
       .not.G%symmetric ) call MOM_error(FATAL, &
       "MOM_open_boundary, open_boundary_config: "//&
       "Symmetric memory must be used when using Flather OBCs.")

  if (.not.(OBC%specified_u_BCs_exist_globally .or. OBC%specified_v_BCs_exist_globally .or. &
              OBC%open_u_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally)) then
    ! No open boundaries have been requested
    call open_boundary_dealloc(OBC)
  else
    ! Need this for ocean_only mode boundary interpolation.
    call time_interp_external_init()
  endif

end subroutine open_boundary_config

!> Allocate space for reading OBC data from files. It sets up the required vertical
!! remapping. In the process, it does funky stuff with the MPI processes.
subroutine initialize_segment_data(G, OBC, PF)
  use mpp_mod, only : mpp_pe, mpp_set_current_pelist, mpp_get_current_pelist,mpp_npes

  type(dyn_horgrid_type), intent(in)    :: G   !< Ocean grid structure
  type(ocean_OBC_type),   intent(inout) :: OBC !< Open boundary control structure
  type(param_file_type),  intent(in)    :: PF  !< Parameter file handle

  integer :: n,m,num_fields
  character(len=256) :: segstr, filename
  character(len=20)  :: segnam, suffix
  character(len=32)  :: varnam, fieldname
  real               :: value
  character(len=32), dimension(MAX_OBC_FIELDS) :: fields  ! segment field names
  character(len=128) :: inputdir
  type(OBC_segment_type), pointer :: segment => NULL() ! pointer to segment type list
  character(len=32)  :: remappingScheme
  character(len=256) :: mesg    ! Message for error messages.
  logical :: check_reconstruction, check_remapping, force_bounds_in_subcell
  integer, dimension(4) :: siz,siz2
  integer :: is, ie, js, je
  integer :: isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB
  integer, dimension(:), allocatable :: saved_pelist
  integer :: current_pe
  integer, dimension(1) :: single_pelist
  !will be able to dynamically switch between sub-sampling refined grid data or model grid

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  ! There is a problem with the order of the OBC initialization
  ! with respect to ALE_init. Currently handling this by copying the
  ! param file so that I can use it later in step_MOM in order to finish
  ! initializing segments on the first step.

  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)

  call get_param(PF, mdl, "REMAPPING_SCHEME", remappingScheme, &
          "This sets the reconstruction scheme used "//&
          "for vertical remapping for all variables. "//&
          "It can be one of the following schemes: \n"//&
          trim(remappingSchemesDoc), default=remappingDefaultScheme,do_not_log=.true.)
  call get_param(PF, mdl, "FATAL_CHECK_RECONSTRUCTIONS", check_reconstruction, &
          "If true, cell-by-cell reconstructions are checked for "//&
          "consistency and if non-monotonicity or an inconsistency is "//&
          "detected then a FATAL error is issued.", default=.false.,do_not_log=.true.)
  call get_param(PF, mdl, "FATAL_CHECK_REMAPPING", check_remapping, &
          "If true, the results of remapping are checked for "//&
          "conservation and new extrema and if an inconsistency is "//&
          "detected then a FATAL error is issued.", default=.false.,do_not_log=.true.)
  call get_param(PF, mdl, "REMAP_BOUND_INTERMEDIATE_VALUES", force_bounds_in_subcell, &
          "If true, the values on the intermediate grid used for remapping "//&
          "are forced to be bounded, which might not be the case due to "//&
          "round off.", default=.false.,do_not_log=.true.)
    call get_param(PF, mdl, "BRUSHCUTTER_MODE", OBC%brushcutter_mode, &
         "If true, read external OBC data on the supergrid.", &
         default=.false.)

  allocate(OBC%remap_CS)
  call initialize_remapping(OBC%remap_CS, remappingScheme, boundary_extrapolation = .false., &
       check_reconstruction=check_reconstruction, &
       check_remapping=check_remapping, force_bounds_in_subcell=force_bounds_in_subcell)

  if (OBC%user_BCs_set_globally) return

  ! Try this here just for the documentation. It is repeated below.
  do n=1, OBC%number_of_segments
    segment => OBC%segment(n)
    write(segnam,"('OBC_SEGMENT_',i3.3,'_DATA')") n
    call get_param(PF, mdl, segnam, segstr, 'OBC segment docs')
  enddo

  !< temporarily disable communication in order to read segment data independently

  allocate(saved_pelist(0:mpp_npes()-1))
  call mpp_get_current_pelist(saved_pelist)
  current_pe = mpp_pe()
  single_pelist(1) = current_pe
  call mpp_set_current_pelist(single_pelist)

  do n=1, OBC%number_of_segments
    segment => OBC%segment(n)
    if (.not. segment%values_needed) cycle

    write(segnam,"('OBC_SEGMENT_',i3.3,'_DATA')") n
    write(suffix,"('_segment_',i3.3)") n
    ! needs documentation !!  Yet, unsafe for now, causes grief for
    ! MOM_parameter_docs in circle_obcs on two processes.
!   call get_param(PF, mdl, segnam, segstr, 'xyz')
    ! Clear out any old values
    segstr = ''
    call get_param(PF, mdl, segnam, segstr)
    if (segstr == '') then
      write(mesg,'("No OBC_SEGMENT_XXX_DATA string for OBC segment ",I3)') n
      call MOM_error(FATAL, mesg)
    endif

    call parse_segment_data_str(trim(segstr), fields=fields, num_fields=num_fields)
    if (num_fields == 0) then
      call MOM_mesg('initialize_segment_data: num_fields = 0')
      cycle ! cycle to next segment
    endif

    allocate(segment%field(num_fields))
    segment%num_fields = num_fields

    segment%temp_segment_data_exists=.false.
    segment%salt_segment_data_exists=.false.
!!
! CODE HERE FOR OTHER OPTIONS (CLAMPED, NUDGED,..)
!!

    isd = segment%HI%isd ; ied = segment%HI%ied
    jsd = segment%HI%jsd ; jed = segment%HI%jed
    IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
    JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB

    do m=1,num_fields
      call parse_segment_data_str(trim(segstr), var=trim(fields(m)), value=value, filenam=filename, fieldnam=fieldname)
      if (trim(filename) /= 'none') then
        OBC%update_OBC = .true. ! Data is assumed to be time-dependent if we are reading from file
        OBC%needs_IO_for_data = .true. ! At least one segment is using I/O for OBC data
!       segment%values_needed = .true. ! Indicates that i/o will be needed for this segment
        segment%field(m)%name = trim(fields(m))
        if (segment%field(m)%name == 'TEMP') then
           segment%temp_segment_data_exists=.true.
           segment%t_values_needed = .false.
        endif
        if (segment%field(m)%name == 'SALT') then
           segment%salt_segment_data_exists=.true.
           segment%s_values_needed = .false.
        endif
        filename = trim(inputdir)//trim(filename)
        fieldname = trim(fieldname)//trim(suffix)
        call field_size(filename,fieldname,siz,no_domain=.true.)
!       if (siz(4) == 1) segment%values_needed = .false.
        if (segment%on_pe) then
          if (OBC%brushcutter_mode .and. (modulo(siz(1),2) == 0 .or. modulo(siz(2),2) == 0)) then
            call MOM_error(FATAL,'segment data are not on the supergrid')
          endif
          siz2(1)=1

          if (siz(1)>1) then
            if (OBC%brushcutter_mode) then
              siz2(1)=(siz(1)-1)/2
            else
              siz2(1)=siz(1)
            endif
          endif
          siz2(2)=1
          if (siz(2)>1) then
            if (OBC%brushcutter_mode) then
              siz2(2)=(siz(2)-1)/2
            else
              siz2(2)=siz(2)
            endif
          endif
          siz2(3)=siz(3)

          if (segment%is_E_or_W) then
            if (segment%field(m)%name == 'V') then
              allocate(segment%field(m)%buffer_src(IsdB:IedB,JsdB:JedB,siz2(3)))
              segment%v_values_needed = .false.
            else if (segment%field(m)%name == 'DVDX') then
              allocate(segment%field(m)%buffer_src(IsdB:IedB,JsdB:JedB,siz2(3)))
              segment%g_values_needed = .false.
            else
              allocate(segment%field(m)%buffer_src(IsdB:IedB,jsd:jed,siz2(3)))
              if (segment%field(m)%name == 'U') then
                segment%u_values_needed = .false.
              else if (segment%field(m)%name == 'SSH') then
                segment%z_values_needed = .false.
              else if (segment%field(m)%name == 'TEMP') then
                segment%t_values_needed = .false.
              else if (segment%field(m)%name == 'SALT') then
                segment%s_values_needed = .false.
              endif
            endif
          else
            if (segment%field(m)%name == 'U') then
              allocate(segment%field(m)%buffer_src(IsdB:IedB,JsdB:JedB,siz2(3)))
              segment%u_values_needed = .false.
            else if (segment%field(m)%name == 'DUDY') then
              allocate(segment%field(m)%buffer_src(IsdB:IedB,JsdB:JedB,siz2(3)))
              segment%g_values_needed = .false.
            else
              allocate(segment%field(m)%buffer_src(isd:ied,JsdB:JedB,siz2(3)))
              if (segment%field(m)%name == 'V') then
                segment%v_values_needed = .false.
              else if (segment%field(m)%name == 'SSH') then
                segment%z_values_needed = .false.
              else if (segment%field(m)%name == 'TEMP') then
                segment%t_values_needed = .false.
              else if (segment%field(m)%name == 'SALT') then
                segment%s_values_needed = .false.
              endif
            endif
          endif
          segment%field(m)%buffer_src(:,:,:)=0.0
          segment%field(m)%fid = init_external_field(trim(filename),&
               trim(fieldname),ignore_axis_atts=.true.,threading=SINGLE_FILE)
          if (siz(3) > 1) then
            fieldname = 'dz_'//trim(fieldname)
            call field_size(filename,fieldname,siz,no_domain=.true.)
            if (segment%is_E_or_W) then
              if (segment%field(m)%name == 'V' .or. segment%field(m)%name == 'DVDX') then
                allocate(segment%field(m)%dz_src(IsdB:IedB,JsdB:JedB,siz(3)))
              else
                allocate(segment%field(m)%dz_src(IsdB:IedB,jsd:jed,siz(3)))
              endif
            else
              if (segment%field(m)%name == 'U' .or. segment%field(m)%name == 'DUDY') then
                allocate(segment%field(m)%dz_src(IsdB:IedB,JsdB:JedB,siz(3)))
              else
                allocate(segment%field(m)%dz_src(isd:ied,JsdB:JedB,siz(3)))
              endif
            endif
            segment%field(m)%dz_src(:,:,:)=0.0
            segment%field(m)%nk_src=siz(3)
            segment%field(m)%fid_dz = init_external_field(trim(filename),trim(fieldname),&
                       ignore_axis_atts=.true.,threading=SINGLE_FILE)
          else
            segment%field(m)%nk_src=1
          endif
        endif
      else
        segment%field(m)%fid = -1
        segment%field(m)%value = value
        segment%field(m)%name = trim(fields(m))
        if (segment%field(m)%name == 'U') then
          segment%u_values_needed = .false.
        elseif (segment%field(m)%name == 'V') then
          segment%v_values_needed = .false.
        elseif (segment%field(m)%name == 'SSH') then
          segment%z_values_needed = .false.
        elseif (segment%field(m)%name == 'TEMP') then
          segment%t_values_needed = .false.
        elseif (segment%field(m)%name == 'SALT') then
          segment%s_values_needed = .false.
        elseif (segment%field(m)%name == 'DVDX' .or. segment%field(m)%name == 'DUDY') then
          segment%g_values_needed = .false.
        endif
      endif
    enddo
    if (segment%u_values_needed .or. segment%v_values_needed .or. &
        segment%t_values_needed .or. segment%s_values_needed .or. &
        segment%z_values_needed .or. segment%g_values_needed) then
      write(mesg,'("Values needed for OBC segment ",I3)') n
      call MOM_error(FATAL, mesg)
    endif
  enddo

  call mpp_set_current_pelist(saved_pelist)

end subroutine initialize_segment_data

!> Define indices for segment and store in hor_index_type
!> using global segment bounds corresponding to q-points
subroutine setup_segment_indices(G, seg, Is_obc, Ie_obc, Js_obc, Je_obc)
  type(dyn_horgrid_type), intent(in) :: G !< grid type
  type(OBC_segment_type), intent(inout) :: seg  !< Open boundary segment
  integer, intent(in) :: Is_obc !< Q-point global i-index of start of segment
  integer, intent(in) :: Ie_obc !< Q-point global i-index of end of segment
  integer, intent(in) :: Js_obc !< Q-point global j-index of start of segment
  integer, intent(in) :: Je_obc !< Q-point global j-index of end of segment
  ! Local variables
  integer :: Isg,Ieg,Jsg,Jeg

  ! Isg, Ieg will be I*_obc in global space
  if (Ie_obc<Is_obc) then
    Isg=Ie_obc;Ieg=Is_obc
  else
    Isg=Is_obc;Ieg=Ie_obc
  endif
  if (Je_obc<Js_obc) then
    Jsg=Je_obc;Jeg=Js_obc
  else
    Jsg=Js_obc;Jeg=Je_obc
  endif

  ! Global space I*_obc but sorted
  seg%HI%IsgB = Isg ; seg%HI%IegB = Ieg
  seg%HI%isg = Isg+1 ; seg%HI%ieg = Ieg
  seg%HI%JsgB = Jsg ; seg%HI%JegB = Jeg
  seg%HI%jsg = Jsg+1 ; seg%HI%Jeg = Jeg

  ! Move into local index space
  Isg = Isg - G%idg_offset
  Jsg = Jsg - G%jdg_offset
  Ieg = Ieg - G%idg_offset
  Jeg = Jeg - G%jdg_offset

  ! This is the i-extent of the segment on this PE.
  ! The values are nonsense if the segment is not on this PE.
  seg%HI%IsdB = min( max(Isg, G%HI%IsdB), G%HI%IedB)
  seg%HI%IedB = min( max(Ieg, G%HI%IsdB), G%HI%IedB)
  seg%HI%isd = min( max(Isg+1, G%HI%isd), G%HI%ied)
  seg%HI%ied = min( max(Ieg, G%HI%isd), G%HI%ied)
  seg%HI%IscB = min( max(Isg, G%HI%IscB), G%HI%IecB)
  seg%HI%IecB = min( max(Ieg, G%HI%IscB), G%HI%IecB)
  seg%HI%isc = min( max(Isg+1, G%HI%isc), G%HI%iec)
  seg%HI%iec = min( max(Ieg, G%HI%isc), G%HI%iec)

  ! This is the j-extent of the segment on this PE.
  ! The values are nonsense if the segment is not on this PE.
  seg%HI%JsdB = min( max(Jsg, G%HI%JsdB), G%HI%JedB)
  seg%HI%JedB = min( max(Jeg, G%HI%JsdB), G%HI%JedB)
  seg%HI%jsd = min( max(Jsg+1, G%HI%jsd), G%HI%jed)
  seg%HI%jed = min( max(Jeg, G%HI%jsd), G%HI%jed)
  seg%HI%JscB = min( max(Jsg, G%HI%JscB), G%HI%JecB)
  seg%HI%JecB = min( max(Jeg, G%HI%JscB), G%HI%JecB)
  seg%HI%jsc = min( max(Jsg+1, G%HI%jsc), G%HI%jec)
  seg%HI%jec = min( max(Jeg, G%HI%jsc), G%HI%jec)

end subroutine setup_segment_indices

!> Parse an OBC_SEGMENT_%%% string starting with "I=" and configure placement and type of OBC accordingly
subroutine setup_u_point_obc(OBC, G, US, segment_str, l_seg, PF, reentrant_y)
  type(ocean_OBC_type),    pointer    :: OBC !< Open boundary control structure
  type(dyn_horgrid_type),  intent(in) :: G   !< Ocean grid structure
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  character(len=*),        intent(in) :: segment_str !< A string in form of "I=%,J=%:%,string"
  integer,                 intent(in) :: l_seg !< which segment is this?
  type(param_file_type), intent(in)   :: PF  !< Parameter file handle
  logical, intent(in)                 :: reentrant_y !< is the domain reentrant in y?
  ! Local variables
  integer :: I_obc, Js_obc, Je_obc ! Position of segment in global index space
  integer :: j, a_loop
  character(len=32) :: action_str(8)
  character(len=128) :: segment_param_str
  real, allocatable, dimension(:)  :: tnudge
  ! This returns the global indices for the segment
  call parse_segment_str(G%ieg, G%jeg, segment_str, I_obc, Js_obc, Je_obc, action_str, reentrant_y)

  call setup_segment_indices(G, OBC%segment(l_seg),I_obc,I_obc,Js_obc,Je_obc)

  I_obc = I_obc - G%idg_offset ! Convert to local tile indices on this tile
  Js_obc = Js_obc - G%jdg_offset ! Convert to local tile indices on this tile
  Je_obc = Je_obc - G%jdg_offset ! Convert to local tile indices on this tile

  if (Je_obc>Js_obc) then
    OBC%segment(l_seg)%direction = OBC_DIRECTION_E
  elseif (Je_obc<Js_obc) then
    OBC%segment(l_seg)%direction = OBC_DIRECTION_W
    j=js_obc;js_obc=je_obc;je_obc=j
  endif

  OBC%segment(l_seg)%on_pe = .false.

  do a_loop = 1,8 ! up to 8 options available
    if (len_trim(action_str(a_loop)) == 0) then
      cycle
    elseif (trim(action_str(a_loop)) == 'FLATHER') then
      OBC%segment(l_seg)%Flather = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%Flather_u_BCs_exist_globally = .true.
      OBC%open_u_BCs_exist_globally = .true.
      OBC%segment%z_values_needed = .true.
      OBC%segment%u_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'ORLANSKI') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%open_u_BCs_exist_globally = .true.
      OBC%radiation_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'ORLANSKI_TAN') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%segment(l_seg)%radiation_tan = .true.
      OBC%radiation_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'ORLANSKI_GRAD') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%segment(l_seg)%radiation_grad = .true.
    elseif (trim(action_str(a_loop)) == 'OBLIQUE') then
      OBC%segment(l_seg)%oblique = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%oblique_BCs_exist_globally = .true.
      OBC%open_u_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'OBLIQUE_TAN') then
      OBC%segment(l_seg)%oblique = .true.
      OBC%segment(l_seg)%oblique_tan = .true.
      OBC%oblique_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'OBLIQUE_GRAD') then
      OBC%segment(l_seg)%oblique = .true.
      OBC%segment(l_seg)%oblique_grad = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED') then
      OBC%segment(l_seg)%nudged = .true.
      OBC%nudged_u_BCs_exist_globally = .true.
      OBC%segment%u_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED_TAN') then
      OBC%segment(l_seg)%nudged_tan = .true.
      OBC%nudged_u_BCs_exist_globally = .true.
      OBC%segment%v_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED_GRAD') then
      OBC%segment(l_seg)%nudged_grad = .true.
      OBC%segment%g_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'GRADIENT') then
      OBC%segment(l_seg)%gradient = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%open_u_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE') then
      OBC%segment(l_seg)%specified = .true.
      OBC%specified_u_BCs_exist_globally = .true. ! This avoids deallocation
      OBC%segment%u_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE_TAN') then
      OBC%segment(l_seg)%specified_tan = .true.
      OBC%segment%v_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE_GRAD') then
      OBC%segment(l_seg)%specified_grad = .true.
      OBC%segment%g_values_needed = .true.
    else
      call MOM_error(FATAL, "MOM_open_boundary.F90, setup_u_point_obc: "//&
                     "String '"//trim(action_str(a_loop))//"' not understood.")
    endif
    if (OBC%segment(l_seg)%nudged .or. OBC%segment(l_seg)%nudged_tan) then
      write(segment_param_str(1:43),"('OBC_SEGMENT_',i3.3,'_VELOCITY_NUDGING_TIMESCALES')") l_seg
      allocate(tnudge(2))
      call get_param(PF, mdl, segment_param_str(1:43), tnudge, &
                     "Timescales in days for nudging along a segment, "//&
                     "for inflow, then outflow. Setting both to zero should "//&
                     "behave like SIMPLE obcs for the baroclinic velocities.", &
                     fail_if_missing=.true., default=0., units="days", scale=86400.0*US%s_to_T)
      OBC%segment(l_seg)%Velocity_nudging_timescale_in = tnudge(1)
      OBC%segment(l_seg)%Velocity_nudging_timescale_out = tnudge(2)
      deallocate(tnudge)
    endif

  enddo ! a_loop

  OBC%segment(l_seg)%is_E_or_W_2 = .true.

  if (I_obc<=G%HI%IsdB+1 .or. I_obc>=G%HI%IedB-1) return ! Boundary is not on tile
  if (Je_obc<=G%HI%JsdB .or. Js_obc>=G%HI%JedB) return ! Segment is not on tile

  OBC%segment(l_seg)%on_pe = .true.
  OBC%segment(l_seg)%is_E_or_W = .true.

  do j=G%HI%jsd, G%HI%jed
    if (j>Js_obc .and. j<=Je_obc) then
      OBC%segnum_u(I_obc,j) = l_seg
    endif
  enddo
  OBC%segment(l_seg)%Is_obc = I_obc
  OBC%segment(l_seg)%Ie_obc = I_obc
  OBC%segment(l_seg)%Js_obc = Js_obc
  OBC%segment(l_seg)%Je_obc = Je_obc
  call allocate_OBC_segment_data(OBC, OBC%segment(l_seg))

  if (OBC%segment(l_seg)%oblique .and.  OBC%segment(l_seg)%radiation) &
         call MOM_error(FATAL, "MOM_open_boundary.F90, setup_u_point_obc: \n"//&
         "Orlanski and Oblique OBC options cannot be used together on one segment.")

  if (OBC%segment(l_seg)%u_values_needed .or. OBC%segment(l_seg)%v_values_needed .or. &
      OBC%segment(l_seg)%t_values_needed .or. OBC%segment(l_seg)%s_values_needed .or. &
      OBC%segment(l_seg)%z_values_needed .or. OBC%segment(l_seg)%g_values_needed) &
    OBC%segment(l_seg)%values_needed = .true.
end subroutine setup_u_point_obc

!> Parse an OBC_SEGMENT_%%% string starting with "J=" and configure placement and type of OBC accordingly
subroutine setup_v_point_obc(OBC, G, US, segment_str, l_seg, PF, reentrant_x)
  type(ocean_OBC_type),    pointer    :: OBC !< Open boundary control structure
  type(dyn_horgrid_type),  intent(in) :: G   !< Ocean grid structure
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  character(len=*),        intent(in) :: segment_str !< A string in form of "J=%,I=%:%,string"
  integer,                 intent(in) :: l_seg !< which segment is this?
  type(param_file_type),   intent(in) :: PF  !< Parameter file handle
  logical, intent(in)                 :: reentrant_x !< is the domain reentrant in x?
  ! Local variables
  integer :: J_obc, Is_obc, Ie_obc ! Position of segment in global index space
  integer :: i, a_loop
  character(len=32) :: action_str(8)
  character(len=128) :: segment_param_str
  real, allocatable, dimension(:)  :: tnudge

  ! This returns the global indices for the segment
  call parse_segment_str(G%ieg, G%jeg, segment_str, J_obc, Is_obc, Ie_obc, action_str, reentrant_x)

  call setup_segment_indices(G, OBC%segment(l_seg),Is_obc,Ie_obc,J_obc,J_obc)

  J_obc = J_obc - G%jdg_offset ! Convert to local tile indices on this tile
  Is_obc = Is_obc - G%idg_offset ! Convert to local tile indices on this tile
  Ie_obc = Ie_obc - G%idg_offset ! Convert to local tile indices on this tile

  if (Ie_obc>Is_obc) then
     OBC%segment(l_seg)%direction = OBC_DIRECTION_S
  elseif (Ie_obc<Is_obc) then
     OBC%segment(l_seg)%direction = OBC_DIRECTION_N
     i=Is_obc;Is_obc=Ie_obc;Ie_obc=i
  endif

  OBC%segment(l_seg)%on_pe = .false.

  do a_loop = 1,8
    if (len_trim(action_str(a_loop)) == 0) then
      cycle
    elseif (trim(action_str(a_loop)) == 'FLATHER') then
      OBC%segment(l_seg)%Flather = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%Flather_v_BCs_exist_globally = .true.
      OBC%open_v_BCs_exist_globally = .true.
      OBC%segment%z_values_needed = .true.
      OBC%segment%v_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'ORLANSKI') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%open_v_BCs_exist_globally = .true.
      OBC%radiation_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'ORLANSKI_TAN') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%segment(l_seg)%radiation_tan = .true.
      OBC%radiation_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'ORLANSKI_GRAD') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%segment(l_seg)%radiation_grad = .true.
    elseif (trim(action_str(a_loop)) == 'OBLIQUE') then
      OBC%segment(l_seg)%oblique = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%oblique_BCs_exist_globally = .true.
      OBC%open_v_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'OBLIQUE_TAN') then
      OBC%segment(l_seg)%oblique = .true.
      OBC%segment(l_seg)%oblique_tan = .true.
      OBC%oblique_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'OBLIQUE_GRAD') then
      OBC%segment(l_seg)%oblique = .true.
      OBC%segment(l_seg)%oblique_grad = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED') then
      OBC%segment(l_seg)%nudged = .true.
      OBC%nudged_v_BCs_exist_globally = .true.
      OBC%segment%v_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED_TAN') then
      OBC%segment(l_seg)%nudged_tan = .true.
      OBC%nudged_v_BCs_exist_globally = .true.
      OBC%segment%u_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED_GRAD') then
      OBC%segment(l_seg)%nudged_grad = .true.
      OBC%segment%g_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'GRADIENT') then
      OBC%segment(l_seg)%gradient = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%open_v_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE') then
      OBC%segment(l_seg)%specified = .true.
      OBC%specified_v_BCs_exist_globally = .true. ! This avoids deallocation
      OBC%segment%v_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE_TAN') then
      OBC%segment(l_seg)%specified_tan = .true.
      OBC%segment%u_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE_GRAD') then
      OBC%segment(l_seg)%specified_grad = .true.
      OBC%segment%g_values_needed = .true.
    else
      call MOM_error(FATAL, "MOM_open_boundary.F90, setup_v_point_obc: "//&
                     "String '"//trim(action_str(a_loop))//"' not understood.")
    endif
    if (OBC%segment(l_seg)%nudged .or. OBC%segment(l_seg)%nudged_tan) then
      write(segment_param_str(1:43),"('OBC_SEGMENT_',i3.3,'_VELOCITY_NUDGING_TIMESCALES')") l_seg
      allocate(tnudge(2))
      call get_param(PF, mdl, segment_param_str(1:43), tnudge, &
                     "Timescales in days for nudging along a segment, "//&
                     "for inflow, then outflow. Setting both to zero should "//&
                     "behave like SIMPLE obcs for the baroclinic velocities.", &
                     fail_if_missing=.true., default=0., units="days", scale=86400.0*US%s_to_T)
      OBC%segment(l_seg)%Velocity_nudging_timescale_in = tnudge(1)
      OBC%segment(l_seg)%Velocity_nudging_timescale_out = tnudge(2)
      deallocate(tnudge)
    endif

  enddo ! a_loop

  if (J_obc<=G%HI%JsdB+1 .or. J_obc>=G%HI%JedB-1) return ! Boundary is not on tile
  if (Ie_obc<=G%HI%IsdB .or. Is_obc>=G%HI%IedB) return ! Segment is not on tile

  OBC%segment(l_seg)%on_pe = .true.
  OBC%segment(l_seg)%is_N_or_S = .true.

  do i=G%HI%isd, G%HI%ied
    if (i>Is_obc .and. i<=Ie_obc) then
      OBC%segnum_v(i,J_obc) = l_seg
    endif
  enddo
  OBC%segment(l_seg)%Is_obc = Is_obc
  OBC%segment(l_seg)%Ie_obc = Ie_obc
  OBC%segment(l_seg)%Js_obc = J_obc
  OBC%segment(l_seg)%Je_obc = J_obc
  call allocate_OBC_segment_data(OBC, OBC%segment(l_seg))

  if (OBC%segment(l_seg)%oblique .and.  OBC%segment(l_seg)%radiation) &
         call MOM_error(FATAL, "MOM_open_boundary.F90, setup_v_point_obc: \n"//&
         "Orlanski and Oblique OBC options cannot be used together on one segment.")

  if (OBC%segment(l_seg)%u_values_needed .or. OBC%segment(l_seg)%v_values_needed .or. &
      OBC%segment(l_seg)%t_values_needed .or. OBC%segment(l_seg)%s_values_needed .or. &
      OBC%segment(l_seg)%z_values_needed .or. OBC%segment(l_seg)%g_values_needed) &
    OBC%segment(l_seg)%values_needed = .true.
end subroutine setup_v_point_obc

!> Parse an OBC_SEGMENT_%%% string
subroutine parse_segment_str(ni_global, nj_global, segment_str, l, m, n, action_str, reentrant)
  integer,          intent(in)  :: ni_global !< Number of h-points in zonal direction
  integer,          intent(in)  :: nj_global !< Number of h-points in meridional direction
  character(len=*), intent(in)  :: segment_str !< A string in form of "I=l,J=m:n,string" or "J=l,I=m,n,string"
  integer,          intent(out) :: l !< The value of I=l, if segment_str begins with I=l, or the value of J=l
  integer,          intent(out) :: m !< The value of J=m, if segment_str begins with I=, or the value of I=m
  integer,          intent(out) :: n !< The value of J=n, if segment_str begins with I=, or the value of I=n
  character(len=*), intent(out) :: action_str(:) !< The "string" part of segment_str
  logical,          intent(in)  :: reentrant !< is domain reentrant in relevant direction?
  ! Local variables
  character(len=24) :: word1, word2, m_word, n_word !< Words delineated by commas in a string in form of
                                                    !! "I=%,J=%:%,string"
  integer :: l_max !< Either ni_global or nj_global, depending on whether segment_str begins with "I=" or "J="
  integer :: mn_max !< Either nj_global or ni_global, depending on whether segment_str begins with "I=" or "J="
  integer :: j
  integer, parameter :: halo = 10

  ! Process first word which will started with either 'I=' or 'J='
  word1 = extract_word(segment_str,',',1)
  word2 = extract_word(segment_str,',',2)
  if (word1(1:2)=='I=') then
    l_max = ni_global
    mn_max = nj_global
    if (.not. (word2(1:2)=='J=')) call MOM_error(FATAL, "MOM_open_boundary.F90, parse_segment_str: "//&
                     "Second word of string '"//trim(segment_str)//"' must start with 'J='.")
  elseif (word1(1:2)=='J=') then ! Note that the file_parser uniformaly expands "=" to " = "
    l_max = nj_global
    mn_max = ni_global
    if (.not. (word2(1:2)=='I=')) call MOM_error(FATAL, "MOM_open_boundary.F90, parse_segment_str: "//&
                     "Second word of string '"//trim(segment_str)//"' must start with 'I='.")
  else
    call MOM_error(FATAL, "MOM_open_boundary.F90, parse_segment_str"//&
                   "String '"//segment_str//"' must start with 'I=' or 'J='.")
  endif

  ! Read l
  l = interpret_int_expr( word1(3:24), l_max )
  if (l<0 .or. l>l_max) then
    call MOM_error(FATAL, "MOM_open_boundary.F90, parse_segment_str: "//&
                   "First value from string '"//trim(segment_str)//"' is outside of the physical domain.")
  endif

  ! Read m
  m_word = extract_word(word2(3:24),':',1)
  m = interpret_int_expr( m_word, mn_max )
  if (reentrant) then
    if (m<-halo .or. m>mn_max+halo) then
      call MOM_error(FATAL, "MOM_open_boundary.F90, parse_segment_str: "//&
                     "Beginning of range in string '"//trim(segment_str)//"' is outside of the physical domain.")
    endif
  else
    if (m<-1 .or. m>mn_max+1) then
      call MOM_error(FATAL, "MOM_open_boundary.F90, parse_segment_str: "//&
                     "Beginning of range in string '"//trim(segment_str)//"' is outside of the physical domain.")
    endif
  endif

  ! Read n
  n_word = extract_word(word2(3:24),':',2)
  n = interpret_int_expr( n_word, mn_max )
  if (reentrant) then
    if (n<-halo .or. n>mn_max+halo) then
      call MOM_error(FATAL, "MOM_open_boundary.F90, parse_segment_str: "//&
                     "End of range in string '"//trim(segment_str)//"' is outside of the physical domain.")
    endif
  else
    if (n<-1 .or. n>mn_max+1) then
      call MOM_error(FATAL, "MOM_open_boundary.F90, parse_segment_str: "//&
                     "End of range in string '"//trim(segment_str)//"' is outside of the physical domain.")
    endif
  endif

  if (abs(n-m)==0) then
    call MOM_error(FATAL, "MOM_open_boundary.F90, parse_segment_str: "//&
                   "Range in string '"//trim(segment_str)//"' must span one cell.")
  endif

  ! Type of open boundary condition
  do j = 1, size(action_str)
    action_str(j) = extract_word(segment_str,',',2+j)
  enddo

  contains

  ! Returns integer value interpreted from string in form of %I, N or N+-%I
  integer function interpret_int_expr(string, imax)
    character(len=*), intent(in) :: string !< Integer in form or %I, N or N-%I
    integer,          intent(in) :: imax !< Value to replace 'N' with
    ! Local variables
    integer slen

    slen = len_trim(string)
    if (slen==0) call MOM_error(FATAL, "MOM_open_boundary.F90, parse_segment_str"//&
                                "Parsed string was empty!")
    if (len_trim(string)==1 .and. string(1:1)=='N') then
      interpret_int_expr = imax
    elseif (string(1:1)=='N') then
      if (string(2:2)=='+') then
        read(string(3:slen),*,err=911) interpret_int_expr
        interpret_int_expr = imax + interpret_int_expr
      elseif (string(2:2)=='-') then
        read(string(3:slen),*,err=911) interpret_int_expr
        interpret_int_expr = imax - interpret_int_expr
      endif
    else
      read(string(1:slen),*,err=911) interpret_int_expr
    endif
    return
    911 call MOM_error(FATAL, "MOM_open_boundary.F90, parse_segment_str"//&
                       "Problem reading value from string '"//trim(string)//"'.")
  end function interpret_int_expr
end subroutine parse_segment_str

!> Parse an OBC_SEGMENT_%%%_DATA string
 subroutine parse_segment_data_str(segment_str, var, value, filenam, fieldnam, fields, num_fields, debug )
   character(len=*),           intent(in)   :: segment_str !< A string in form of
                                                          !! "VAR1=file:foo1.nc(varnam1),VAR2=file:foo2.nc(varnam2),..."
   character(len=*), optional, intent(in)   :: var        !< The name of the variable for which parameters are needed
   character(len=*), optional, intent(out)  :: filenam    !< The name of the input file if using "file" method
   character(len=*), optional, intent(out)  :: fieldnam   !< The name of the variable in the input file if using
                                                          !! "file" method
   real,             optional, intent(out)  :: value      !< A constant value if using the "value" method
   character(len=*), dimension(MAX_OBC_FIELDS), &
                     optional, intent(out)  :: fields     !< List of fieldnames for each segment
   integer, optional, intent(out)           :: num_fields !< The number of fields in the segment data
   logical, optional, intent(in)            :: debug      !< If present and true, write verbose debugging messages
   ! Local variables
   character(len=128) :: word1, word2, word3, method
   integer :: lword, nfields, n, m
   logical :: continue,dbg
   character(len=32), dimension(MAX_OBC_FIELDS) :: flds

   nfields=0
   continue=.true.
   dbg=.false.
   if (PRESENT(debug)) dbg=debug

   do while (continue)
      word1 = extract_word(segment_str,',',nfields+1)
      if (trim(word1) == '') exit
      nfields=nfields+1
      word2 = extract_word(word1,'=',1)
      flds(nfields) = trim(word2)
   enddo

   if (PRESENT(fields)) then
     do n=1,nfields
       fields(n) = flds(n)
     enddo
   endif

   if (PRESENT(num_fields)) then
      num_fields=nfields
      return
   endif

   m=0
   if (PRESENT(var)) then
     do n=1,nfields
       if (trim(var)==trim(flds(n))) then
          m=n
          exit
       endif
     enddo
     if (m==0) then
        call abort()
     endif

    ! Process first word which will start with the fieldname
     word3 = extract_word(segment_str,',',m)
     word1 = extract_word(word3,':',1)
!     if (trim(word1) == '') exit
     word2 = extract_word(word1,'=',1)
     if (trim(word2) == trim(var)) then
        method=trim(extract_word(word1,'=',2))
        lword=len_trim(method)
        if (method(lword-3:lword) == 'file') then
           ! raise an error id filename/fieldname not in argument list
           word1 = extract_word(word3,':',2)
           filenam = extract_word(word1,'(',1)
           fieldnam = extract_word(word1,'(',2)
           lword=len_trim(fieldnam)
           fieldnam = fieldnam(1:lword-1)  ! remove trailing parenth
           value=-999.
        elseif (method(lword-4:lword) == 'value') then
           filenam = 'none'
           fieldnam = 'none'
           word1 = extract_word(word3,':',2)
           lword=len_trim(word1)
           read(word1(1:lword),*,end=986,err=987) value
        endif
      endif
    endif

   return
 986 call MOM_error(FATAL,'End of record while parsing segment data specification! '//trim(segment_str))
 987 call MOM_error(FATAL,'Error while parsing segment data specification! '//trim(segment_str))

 end subroutine parse_segment_data_str


!> Parse all the OBC_SEGMENT_%%%_DATA strings again
!! to see which need tracer reservoirs (all pes need to know).
 subroutine parse_for_tracer_reservoirs(OBC, PF, use_temperature)
  type(ocean_OBC_type),   intent(inout) :: OBC !< Open boundary control structure
  type(param_file_type),  intent(in)    :: PF  !< Parameter file handle
  logical,                intent(in) :: use_temperature !< If true, T and S are used

  ! Local variables
  integer :: n,m,num_fields
  character(len=256) :: segstr, filename
  character(len=20)  :: segnam, suffix
  character(len=32)  :: varnam, fieldname
  real               :: value
  character(len=32), dimension(MAX_OBC_FIELDS) :: fields  ! segment field names
  type(OBC_segment_type), pointer :: segment => NULL() ! pointer to segment type list
  character(len=256) :: mesg    ! Message for error messages.

  do n=1, OBC%number_of_segments
    segment => OBC%segment(n)
    write(segnam,"('OBC_SEGMENT_',i3.3,'_DATA')") n
    write(suffix,"('_segment_',i3.3)") n
    ! Clear out any old values
    segstr = ''
    call get_param(PF, mdl, segnam, segstr)
    if (segstr == '') cycle

    call parse_segment_data_str(trim(segstr), fields=fields, num_fields=num_fields)
    if (num_fields == 0) cycle

    ! At this point, just search for TEMP and SALT as tracers 1 and 2.
    do m=1,num_fields
      call parse_segment_data_str(trim(segstr), var=trim(fields(m)), value=value, filenam=filename, fieldnam=fieldname)
      if (trim(filename) /= 'none') then
        if (fields(m) == 'TEMP') then
          if (segment%is_E_or_W_2) then
            OBC%tracer_x_reservoirs_used(1) = .true.
          else
            OBC%tracer_y_reservoirs_used(1) = .true.
          endif
        endif
        if (fields(m) == 'SALT') then
          if (segment%is_E_or_W_2) then
            OBC%tracer_x_reservoirs_used(2) = .true.
          else
            OBC%tracer_y_reservoirs_used(2) = .true.
          endif
        endif
      endif
    enddo
    ! Alternately, set first two to true if use_temperature is true
    if (use_temperature) then
      if (segment%is_E_or_W_2) then
        OBC%tracer_x_reservoirs_used(1) = .true.
        OBC%tracer_x_reservoirs_used(2) = .true.
      else
        OBC%tracer_y_reservoirs_used(1) = .true.
        OBC%tracer_y_reservoirs_used(2) = .true.
      endif
    endif
  enddo

  return

end subroutine parse_for_tracer_reservoirs

!> Parse an OBC_SEGMENT_%%%_PARAMS string
 subroutine parse_segment_param_real(segment_str, var, param_value, debug )
   character(len=*),  intent(in)  :: segment_str !< A string in form of
                                                 !! "VAR1=file:foo1.nc(varnam1),VAR2=file:foo2.nc(varnam2),..."
   character(len=*),  intent(in)  :: var         !< The name of the variable for which parameters are needed
   real,              intent(out) :: param_value !< The value of the parameter
   logical, optional, intent(in)  :: debug       !< If present and true, write verbose debugging messages
   ! Local variables
   character(len=128) :: word1, word2, word3, method
   integer :: lword, nfields, n, m
   logical :: continue,dbg
   character(len=32), dimension(MAX_OBC_FIELDS) :: flds

   nfields=0
   continue=.true.
   dbg=.false.
   if (PRESENT(debug)) dbg=debug

   do while (continue)
      word1 = extract_word(segment_str,',',nfields+1)
      if (trim(word1) == '') exit
      nfields=nfields+1
      word2 = extract_word(word1,'=',1)
      flds(nfields) = trim(word2)
   enddo

   ! if (PRESENT(fields)) then
   !   do n=1,nfields
   !     fields(n) = flds(n)
   !   enddo
   ! endif

   ! if (PRESENT(num_fields)) then
   !    num_fields=nfields
   !    return
   ! endif

   m=0
!   if (PRESENT(var)) then
     do n=1,nfields
       if (trim(var)==trim(flds(n))) then
          m=n
          exit
       endif
     enddo
     if (m==0) then
        call abort()
     endif

    ! Process first word which will start with the fieldname
     word3 = extract_word(segment_str,',',m)
!     word1 = extract_word(word3,':',1)
!     if (trim(word1) == '') exit
     word2 = extract_word(word1,'=',1)
     if (trim(word2) == trim(var)) then
        method=trim(extract_word(word1,'=',2))
        lword=len_trim(method)
        read(method(1:lword),*,err=987) param_value
        ! if (method(lword-3:lword) == 'file') then
        !    ! raise an error id filename/fieldname not in argument list
        !    word1 = extract_word(word3,':',2)
        !    filenam = extract_word(word1,'(',1)
        !    fieldnam = extract_word(word1,'(',2)
        !    lword=len_trim(fieldnam)
        !    fieldnam = fieldnam(1:lword-1)  ! remove trailing parenth
        !    value=-999.
        ! elseif (method(lword-4:lword) == 'value') then
        !    filenam = 'none'
        !    fieldnam = 'none'
        !    word1 = extract_word(word3,':',2)
        !    lword=len_trim(word1)
        !    read(word1(1:lword),*,end=986,err=987) value
        ! endif
      endif
!    endif

   return
 986 call MOM_error(FATAL,'End of record while parsing segment data specification! '//trim(segment_str))
 987 call MOM_error(FATAL,'Error while parsing segment parameter specification! '//trim(segment_str))

 end subroutine parse_segment_param_real

!> Initialize open boundary control structure and do any necessary rescaling of OBC
!! fields that have been read from a restart file.
subroutine open_boundary_init(G, GV, US, param_file, OBC, restart_CSp)
  type(ocean_grid_type),   intent(in) :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV  !< Container for vertical grid information
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  type(param_file_type),   intent(in) :: param_file !< Parameter file handle
  type(ocean_OBC_type),    pointer    :: OBC !< Open boundary control structure
  type(MOM_restart_CS),    pointer    :: restart_CSp !< Restart structure, data intent(inout)

  ! Local variables
  real :: vel2_rescale ! A rescaling factor for squared velocities from the representation in
                       ! a restart file to the internal representation in this run.
  integer :: i, j, k, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) return

  id_clock_pass = cpu_clock_id('(Ocean OBC halo updates)', grain=CLOCK_ROUTINE)

  ! The rx_normal and ry_normal arrays used with radiation OBCs are currently in units of grid
  ! points per timestep, but if this were to be corrected to [L T-1 ~> m s-1] or [T-1 ~> s-1] to
  ! permit timesteps to change between calls to the OBC code, the following would be needed:
!  if ( OBC%radiation_BCs_exist_globally .and. (US%s_to_T_restart * US%m_to_L_restart /= 0.0) .and. &
!       ((US%m_to_L * US%s_to_T_restart) /= (US%m_to_L_restart * US%s_to_T)) ) then
!    vel_rescale = (US%m_to_L * US%s_to_T_restart) /  (US%m_to_L_restart * US%s_to_T)
!    if (query_initialized(OBC%rx_normal, "rx_normal", restart_CSp)) then
!      do k=1,nz ; do j=jsd,jed ; do I=IsdB,IedB
!        OBC%rx_normal(I,j,k) = vel_rescale * OBC%rx_normal(I,j,k)
!      enddo ; enddo ; enddo
!    endif
!    if (query_initialized(OBC%ry_normal, "ry_normal", restart_CSp)) then
!      do k=1,nz ; do J=JsdB,JedB ; do i=isd,ied
!        OBC%ry_normal(i,J,k) = vel_rescale * OBC%ry_normal(i,J,k)
!      enddo ; enddo ; enddo
!    endif
!  endif

  ! The oblique boundary condition terms have units of [L2 T-2 ~> m2 s-2] and may need to be rescaled.
  if ( OBC%oblique_BCs_exist_globally .and. (US%s_to_T_restart * US%m_to_L_restart /= 0.0) .and. &
       ((US%m_to_L * US%s_to_T_restart) /= (US%m_to_L_restart * US%s_to_T)) ) then
    vel2_rescale = (US%m_to_L * US%s_to_T_restart)**2 /  (US%m_to_L_restart * US%s_to_T)**2
    if (query_initialized(OBC%rx_oblique, "rx_oblique", restart_CSp)) then
      do k=1,nz ; do j=jsd,jed ; do I=IsdB,IedB
        OBC%rx_oblique(I,j,k) = vel2_rescale * OBC%rx_oblique(I,j,k)
      enddo ; enddo ; enddo
    endif
    if (query_initialized(OBC%ry_oblique, "ry_oblique", restart_CSp)) then
      do k=1,nz ; do J=JsdB,JedB ; do i=isd,ied
        OBC%ry_oblique(i,J,k) = vel2_rescale * OBC%ry_oblique(i,J,k)
      enddo ; enddo ; enddo
    endif
    if (query_initialized(OBC%cff_normal, "cff_normal", restart_CSp)) then
      do k=1,nz ; do J=JsdB,JedB ; do I=IsdB,IedB
        OBC%cff_normal(I,J,k) = vel2_rescale * OBC%cff_normal(I,J,k)
      enddo ; enddo ; enddo
    endif
   endif

end subroutine open_boundary_init

logical function open_boundary_query(OBC, apply_open_OBC, apply_specified_OBC, apply_Flather_OBC, &
                                     apply_nudged_OBC, needs_ext_seg_data)
  type(ocean_OBC_type), pointer    :: OBC !< Open boundary control structure
  logical, optional,    intent(in) :: apply_open_OBC      !< Returns True if open_*_BCs_exist_globally is true
  logical, optional,    intent(in) :: apply_specified_OBC !< Returns True if specified_*_BCs_exist_globally is true
  logical, optional,    intent(in) :: apply_Flather_OBC   !< Returns True if Flather_*_BCs_exist_globally is true
  logical, optional,    intent(in) :: apply_nudged_OBC    !< Returns True if nudged_*_BCs_exist_globally is true
  logical, optional,    intent(in) :: needs_ext_seg_data  !< Returns True if external segment data needed
  open_boundary_query = .false.
  if (.not. associated(OBC)) return
  if (present(apply_open_OBC)) open_boundary_query = OBC%open_u_BCs_exist_globally .or. &
                                                     OBC%open_v_BCs_exist_globally
  if (present(apply_specified_OBC)) open_boundary_query = OBC%specified_u_BCs_exist_globally .or. &
                                                          OBC%specified_v_BCs_exist_globally
  if (present(apply_Flather_OBC)) open_boundary_query = OBC%Flather_u_BCs_exist_globally .or. &
                                                        OBC%Flather_v_BCs_exist_globally
  if (present(apply_nudged_OBC)) open_boundary_query = OBC%nudged_u_BCs_exist_globally .or. &
                                                       OBC%nudged_v_BCs_exist_globally
  if (present(needs_ext_seg_data)) open_boundary_query = OBC%needs_IO_for_data

end function open_boundary_query

!> Deallocate open boundary data
subroutine open_boundary_dealloc(OBC)
  type(ocean_OBC_type), pointer :: OBC !< Open boundary control structure
  type(OBC_segment_type), pointer :: segment => NULL()
  integer :: n

  if (.not. associated(OBC)) return

  do n=1, OBC%number_of_segments
    segment => OBC%segment(n)
    call deallocate_OBC_segment_data(OBC, segment)
  enddo
  if (associated(OBC%segment)) deallocate(OBC%segment)
  if (associated(OBC%segnum_u)) deallocate(OBC%segnum_u)
  if (associated(OBC%segnum_v)) deallocate(OBC%segnum_v)
  if (associated(OBC%rx_normal)) deallocate(OBC%rx_normal)
  if (associated(OBC%ry_normal)) deallocate(OBC%ry_normal)
  if (associated(OBC%rx_oblique)) deallocate(OBC%rx_oblique)
  if (associated(OBC%ry_oblique)) deallocate(OBC%ry_oblique)
  if (associated(OBC%cff_normal)) deallocate(OBC%cff_normal)
  if (associated(OBC%tres_x)) deallocate(OBC%tres_x)
  if (associated(OBC%tres_y)) deallocate(OBC%tres_y)
  deallocate(OBC)
end subroutine open_boundary_dealloc

!> Close open boundary data
subroutine open_boundary_end(OBC)
  type(ocean_OBC_type), pointer :: OBC !< Open boundary control structure
  call open_boundary_dealloc(OBC)
end subroutine open_boundary_end

!> Sets the slope of bathymetry normal to an open bounndary to zero.
subroutine open_boundary_impose_normal_slope(OBC, G, depth)
  type(ocean_OBC_type),             pointer       :: OBC !< Open boundary control structure
  type(dyn_horgrid_type),           intent(in)    :: G !< Ocean grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: depth !< Bathymetry at h-points
  ! Local variables
  integer :: i, j, n
  type(OBC_segment_type), pointer :: segment => NULL()

  if (.not.associated(OBC)) return

  if (.not.(OBC%specified_u_BCs_exist_globally .or. OBC%specified_v_BCs_exist_globally .or. &
              OBC%open_u_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally)) &
    return

  do n=1,OBC%number_of_segments
    segment=>OBC%segment(n)
    if (.not. segment%on_pe) cycle
    if (segment%direction == OBC_DIRECTION_E) then
      I=segment%HI%IsdB
      do j=segment%HI%jsd,segment%HI%jed
        depth(i+1,j) = depth(i,j)
      enddo
    elseif (segment%direction == OBC_DIRECTION_W) then
      I=segment%HI%IsdB
      do j=segment%HI%jsd,segment%HI%jed
        depth(i,j) = depth(i+1,j)
      enddo
    elseif (segment%direction == OBC_DIRECTION_N) then
      J=segment%HI%JsdB
      do i=segment%HI%isd,segment%HI%ied
        depth(i,j+1) = depth(i,j)
      enddo
    elseif (segment%direction == OBC_DIRECTION_S) then
      J=segment%HI%JsdB
      do i=segment%HI%isd,segment%HI%ied
        depth(i,j) = depth(i,j+1)
      enddo
    endif
  enddo

end subroutine open_boundary_impose_normal_slope

!> Reconcile masks and open boundaries, deallocate OBC on PEs where it is not needed.
!! Also adjust u- and v-point cell area on specified open boundaries and mask all
!! points outside open boundaries.
subroutine open_boundary_impose_land_mask(OBC, G, areaCu, areaCv, US)
  type(ocean_OBC_type),              pointer       :: OBC !< Open boundary control structure
  type(dyn_horgrid_type),            intent(inout) :: G   !< Ocean grid structure
  type(unit_scale_type),             intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G)), intent(inout) :: areaCu !< Area of a u-cell [L2 ~> m2]
  real, dimension(SZI_(G),SZJB_(G)), intent(inout) :: areaCv !< Area of a u-cell [L2 ~> m2]
  ! Local variables
  integer :: i, j, n
  type(OBC_segment_type), pointer :: segment => NULL()
  logical :: any_U, any_V

  if (.not.associated(OBC)) return

  do n=1,OBC%number_of_segments
    segment=>OBC%segment(n)
    if (.not. segment%on_pe) cycle
    if (segment%is_E_or_W) then
      ! Sweep along u-segments and delete the OBC for blocked points.
      ! Also, mask all points outside.
      I=segment%HI%IsdB
      do j=segment%HI%jsd,segment%HI%jed
        if (G%mask2dCu(I,j) == 0) OBC%segnum_u(I,j) = OBC_NONE
        if (segment%direction == OBC_DIRECTION_W) then
          G%mask2dT(i,j) = 0
        else
          G%mask2dT(i+1,j) = 0
        endif
      enddo
      do J=segment%HI%JsdB+1,segment%HI%JedB-1
        if (segment%direction == OBC_DIRECTION_W) then
          G%mask2dCv(i,J) = 0
        else
          G%mask2dCv(i+1,J) = 0
        endif
      enddo
    else
      ! Sweep along v-segments and delete the OBC for blocked points.
      J=segment%HI%JsdB
      do i=segment%HI%isd,segment%HI%ied
        if (G%mask2dCv(i,J) == 0) OBC%segnum_v(i,J) = OBC_NONE
        if (segment%direction == OBC_DIRECTION_S) then
          G%mask2dT(i,j) = 0
        else
          G%mask2dT(i,j+1) = 0
        endif
      enddo
      do I=segment%HI%IsdB+1,segment%HI%IedB-1
        if (segment%direction == OBC_DIRECTION_S) then
          G%mask2dCu(I,j) = 0
        else
          G%mask2dCu(I,j+1) = 0
        endif
      enddo
    endif
  enddo

  do n=1,OBC%number_of_segments
    segment=>OBC%segment(n)
    if (.not. segment%on_pe .or. .not. segment%specified) cycle
    if (segment%is_E_or_W) then
      ! Sweep along u-segments and for %specified BC points reset the u-point area which was masked out
      I=segment%HI%IsdB
      do j=segment%HI%jsd,segment%HI%jed
        if (segment%direction == OBC_DIRECTION_E) then
          areaCu(I,j) = G%areaT(i,j)   ! Both of these are in [L2]
        else   ! West
          areaCu(I,j) = G%areaT(i+1,j) ! Both of these are in [L2]
        endif
      enddo
    else
      ! Sweep along v-segments and for %specified BC points reset the v-point area which was masked out
      J=segment%HI%JsdB
      do i=segment%HI%isd,segment%HI%ied
        if (segment%direction == OBC_DIRECTION_S) then
          areaCv(i,J) = G%areaT(i,j+1) ! Both of these are in [L2]
        else      ! North
          areaCu(i,J) = G%areaT(i,j)   ! Both of these are in [L2]
        endif
      enddo
    endif
  enddo

  ! G%mask2du will be open wherever bathymetry allows it.
  ! Bathymetry outside of the open boundary was adjusted to match
  ! the bathymetry inside so these points will be open unless the
  ! bathymetry inside the boundary was too shallow and flagged as land.
  any_U = .false.
  any_V = .false.
  do n=1,OBC%number_of_segments
    segment=>OBC%segment(n)
    if (.not. segment%on_pe) cycle
    if (segment%is_E_or_W) then
      I=segment%HI%IsdB
      do j=segment%HI%jsd,segment%HI%jed
        if (OBC%segnum_u(I,j) /= OBC_NONE) any_U = .true.
      enddo
    else
      J=segment%HI%JsdB
      do i=segment%HI%isd,segment%HI%ied
        if (OBC%segnum_v(i,J) /= OBC_NONE) any_V = .true.
      enddo
    endif
  enddo

  OBC%OBC_pe = .true.
  if (.not.(any_U .or. any_V)) OBC%OBC_pe = .false.

end subroutine open_boundary_impose_land_mask

!> Make sure the OBC tracer reservoirs are initialized.
subroutine setup_OBC_tracer_reservoirs(G, OBC)
  type(ocean_grid_type),      intent(inout) :: G          !< Ocean grid structure
  type(ocean_OBC_type),       pointer       :: OBC !< Open boundary control structure
  ! Local variables
  type(OBC_segment_type), pointer :: segment => NULL()
  integer :: i, j, k, m, n

  do n=1,OBC%number_of_segments
    segment=>OBC%segment(n)
    if (associated(segment%tr_Reg)) then
      if (segment%is_E_or_W) then
        I = segment%HI%IsdB
        do m=1,OBC%ntr
          if (associated(segment%tr_Reg%Tr(m)%tres)) then
            do k=1,G%ke
              do j=segment%HI%jsd,segment%HI%jed
                OBC%tres_x(I,j,k,m) = segment%tr_Reg%Tr(m)%t(i,j,k)
              enddo
            enddo
          endif
        enddo
      else
        J = segment%HI%JsdB
        do m=1,OBC%ntr
          if (associated(segment%tr_Reg%Tr(m)%tres)) then
            do k=1,G%ke
              do i=segment%HI%isd,segment%HI%ied
                OBC%tres_y(i,J,k,m) = segment%tr_Reg%Tr(m)%t(i,J,k)
              enddo
            enddo
          endif
        enddo
      endif
    endif
  enddo

end subroutine setup_OBC_tracer_reservoirs

!> Apply radiation conditions to 3D  u,v at open boundaries
subroutine radiation_open_bdry_conds(OBC, u_new, u_old, v_new, v_old, G, US, dt)
  type(ocean_grid_type),                     intent(inout) :: G !< Ocean grid structure
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary control structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: u_new !< On exit, new u values on open boundaries
                                                                    !! On entry, the old time-level v but including
                                                                    !! barotropic accelerations [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u_old !< Original unadjusted u [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: v_new !< On exit, new v values on open boundaries.
                                                                    !! On entry, the old time-level v but including
                                                                    !! barotropic accelerations [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v_old !< Original unadjusted v  [L T-1 ~> m s-1]
  type(unit_scale_type),                     intent(in)    :: US    !< A dimensional unit scaling type
  real,                                      intent(in)    :: dt    !< Appropriate timestep [T ~> s]
  ! Local variables
  real :: dhdt, dhdx, dhdy  ! One-point differences in time or space [L T-1 ~> m s-1]
  real :: gamma_u, gamma_2  ! Fractional weightings of new values [nondim]
  real :: tau            ! A local nudging timescale [T ~> s]
  real :: rx_max, ry_max ! coefficients for radiation [nondim] or [L2 T-2 ~> m2 s-2]
  real :: rx_new, rx_avg ! coefficients for radiation [nondim] or [L2 T-2 ~> m2 s-2]
  real :: ry_new, ry_avg ! coefficients for radiation [nondim] or [L2 T-2 ~> m2 s-2]
  real :: cff_new, cff_avg ! denominator in oblique [L2 T-2 ~> m2 s-2]
  real, allocatable, dimension(:,:,:) :: &
    rx_tang_rad, & ! The phase speed at u-points for tangential oblique OBCs
                   ! in units of grid points per timestep [nondim]
    ry_tang_rad, & ! The phase speed at v-points for tangential oblique OBCs
                   ! in units of grid points per timestep [nondim]
    rx_tang_obl, & ! The x-coefficient for tangential oblique OBCs [L2 T-2 ~> m2 s-2]
    ry_tang_obl, & ! The y-coefficient for tangential oblique OBCs [L2 T-2 ~> m2 s-2]
    cff_tangential ! The denominator for tangential oblique OBCs [L2 T-2 ~> m2 s-2]
  real :: eps      ! A small velocity squared [L2 T-2 ~> m2 s-2]
  type(OBC_segment_type), pointer :: segment => NULL()
  integer :: i, j, k, is, ie, js, je, m, nz, n
  integer :: is_obc, ie_obc, js_obc, je_obc

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not.associated(OBC)) return

  if (.not.(OBC%open_u_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally)) &
    return

  eps = 1.0e-20*US%m_s_to_L_T**2

  !! Copy previously calculated phase velocity from global arrays into segments
  !! This is terribly inefficient and temporary solution for continuity across restarts
  !! and needs to be revisited in the future.
  if (OBC%gamma_uv < 1.0) then
    do n=1,OBC%number_of_segments
      segment=>OBC%segment(n)
      if (.not. segment%on_pe) cycle
      if (segment%is_E_or_W .and. segment%radiation) then
        do k=1,G%ke
          I=segment%HI%IsdB
          do j=segment%HI%jsd,segment%HI%jed
            segment%rx_norm_rad(I,j,k) = OBC%rx_normal(I,j,k)
          enddo
        enddo
      elseif (segment%is_N_or_S .and. segment%radiation) then
        do k=1,G%ke
          J=segment%HI%JsdB
          do i=segment%HI%isd,segment%HI%ied
            segment%ry_norm_rad(i,J,k) = OBC%ry_normal(i,J,k)
          enddo
        enddo
      endif
      if (segment%is_E_or_W .and. segment%oblique) then
        do k=1,G%ke
          I=segment%HI%IsdB
          do j=segment%HI%jsd,segment%HI%jed
            segment%rx_norm_obl(I,j,k) = OBC%rx_oblique(I,j,k)
            segment%ry_norm_obl(I,j,k) = OBC%ry_oblique(I,j,k)
            segment%cff_normal(I,j,k) = OBC%cff_normal(I,j,k)
          enddo
        enddo
      elseif (segment%is_N_or_S .and. segment%oblique) then
        do k=1,G%ke
          J=segment%HI%JsdB
          do i=segment%HI%isd,segment%HI%ied
            segment%rx_norm_obl(i,J,k) = OBC%rx_oblique(i,J,k)
            segment%ry_norm_obl(i,J,k) = OBC%ry_oblique(i,J,k)
            segment%cff_normal(i,J,k) = OBC%cff_normal(i,J,k)
          enddo
        enddo
      endif
    enddo
  endif

  ! Now tracers (if any)
  do n=1,OBC%number_of_segments
    segment=>OBC%segment(n)
    if (associated(segment%tr_Reg)) then
      if (segment%is_E_or_W) then
        I = segment%HI%IsdB
        do m=1,OBC%ntr
          if (associated(segment%tr_Reg%Tr(m)%tres)) then
            do k=1,G%ke
              do j=segment%HI%jsd,segment%HI%jed
                segment%tr_Reg%Tr(m)%tres(I,j,k) = OBC%tres_x(I,j,k,m)
              enddo
            enddo
          endif
        enddo
      else
        J = segment%HI%JsdB
        do m=1,OBC%ntr
          if (associated(segment%tr_Reg%Tr(m)%tres)) then
            do k=1,G%ke
              do i=segment%HI%isd,segment%HI%ied
                segment%tr_Reg%Tr(m)%tres(i,J,k) = OBC%tres_y(i,J,k,m)
              enddo
            enddo
          endif
        enddo
      endif
    endif
  enddo

  gamma_u = OBC%gamma_uv
  rx_max = OBC%rx_max ; ry_max = OBC%rx_max
  do n=1,OBC%number_of_segments
     segment=>OBC%segment(n)
     if (.not. segment%on_pe) cycle
     if (segment%oblique) call gradient_at_q_points(G, segment, u_new(:,:,:), v_new(:,:,:))
     if (segment%direction == OBC_DIRECTION_E) then
       I=segment%HI%IsdB
       if (I<G%HI%IscB) cycle
       do k=1,nz ;  do j=segment%HI%jsd,segment%HI%jed
         if (segment%radiation) then
           dhdt = (u_old(I-1,j,k) - u_new(I-1,j,k)) !old-new
           dhdx = (u_new(I-1,j,k) - u_new(I-2,j,k)) !in new time backward sasha for I-1
           rx_new = 0.0
           if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max) ! outward phase speed
           if (gamma_u < 1.0) then
             rx_avg = (1.0-gamma_u)*segment%rx_norm_rad(I,j,k) + gamma_u*rx_new
           else
             rx_avg = rx_new
           endif
           segment%rx_norm_rad(I,j,k) = rx_avg
           ! The new boundary value is interpolated between future interior
           ! value, u_new(I-1) and past boundary value but with barotropic
           ! accelerations, u_new(I).
           segment%normal_vel(I,j,k) = (u_new(I,j,k) + rx_avg*u_new(I-1,j,k)) / (1.0+rx_avg)
           ! Copy restart fields into 3-d arrays. This is an inefficient and temporary issues
           ! implemented as a work-around to limitations in restart capability
           if (gamma_u < 1.0) then
             OBC%rx_normal(I,j,k) = segment%rx_norm_rad(I,j,k)
           endif
         elseif (segment%oblique) then
           dhdt = (u_old(I-1,j,k) - u_new(I-1,j,k)) !old-new
           dhdx = (u_new(I-1,j,k) - u_new(I-2,j,k)) !in new time backward sasha for I-1
           if (dhdt*(segment%grad_normal(J,1,k) + segment%grad_normal(J-1,1,k)) > 0.0) then
             dhdy = segment%grad_normal(J-1,1,k)
           elseif (dhdt*(segment%grad_normal(J,1,k) + segment%grad_normal(J-1,1,k)) == 0.0) then
             dhdy = 0.0
           else
             dhdy = segment%grad_normal(J,1,k)
           endif
           if (dhdt*dhdx < 0.0) dhdt = 0.0
           rx_new = dhdt*dhdx
           cff_new = max(dhdx*dhdx + dhdy*dhdy, eps)
           ry_new = min(cff_new,max(dhdt*dhdy,-cff_new))
           if (gamma_u < 1.0) then
             rx_avg = (1.0-gamma_u)*segment%rx_norm_obl(I,j,k) + gamma_u*rx_new
             ry_avg = (1.0-gamma_u)*segment%ry_norm_obl(i,J,k) + gamma_u*ry_new
             cff_avg = (1.0-gamma_u)*segment%cff_normal(i,J,k) + gamma_u*cff_new
           else
             rx_avg = rx_new
             ry_avg = ry_new
             cff_avg = cff_new
           endif
           segment%rx_norm_obl(I,j,k) = rx_avg
           segment%ry_norm_obl(i,J,k) = ry_avg
           segment%cff_normal(i,J,k) = cff_avg
           segment%normal_vel(I,j,k) = ((cff_avg*u_new(I,j,k) + rx_avg*u_new(I-1,j,k)) - &
                              (max(ry_avg,0.0)*segment%grad_normal(J-1,2,k) + &
                               min(ry_avg,0.0)*segment%grad_normal(J,2,k))) / &
                            (cff_avg + rx_avg)
           if (gamma_u < 1.0) then
             ! Copy restart fields into 3-d arrays. This is an inefficient and temporary
             ! implementation as a work-around to limitations in restart capability
             OBC%rx_oblique(I,j,k) = segment%rx_norm_obl(I,j,k)
             OBC%ry_oblique(i,J,k) = segment%ry_norm_obl(i,J,k)
             OBC%cff_normal(I,j,k) = segment%cff_normal(I,j,k)
           endif
         elseif (segment%gradient) then
           segment%normal_vel(I,j,k) = u_new(I-1,j,k)
         endif
         if ((segment%radiation .or. segment%oblique) .and. segment%nudged) then
           ! dhdt gets set to 0 on inflow in oblique case
           if (dhdt*dhdx <= 0.0) then
             tau = segment%Velocity_nudging_timescale_in
           else
             tau = segment%Velocity_nudging_timescale_out
           endif
           gamma_2 = dt / (tau + dt)
           segment%normal_vel(I,j,k) = (1.0 - gamma_2) * segment%normal_vel(I,j,k) + &
                                 gamma_2 * segment%nudged_normal_vel(I,j,k)
         endif
       enddo ; enddo
       if (segment%radiation_tan .or. segment%radiation_grad) then
         I=segment%HI%IsdB
         allocate(rx_tang_rad(segment%HI%IsdB:segment%HI%IedB,segment%HI%JsdB:segment%HI%JedB,nz))
         do k=1,nz
           if (gamma_u < 1.0) then
             rx_tang_rad(I,segment%HI%JsdB,k) = segment%rx_norm_rad(I,segment%HI%jsd,k)
             rx_tang_rad(I,segment%HI%JedB,k) = segment%rx_norm_rad(I,segment%HI%jed,k)
             do J=segment%HI%JsdB+1,segment%HI%JedB-1
               rx_tang_rad(I,J,k) = 0.5*(segment%rx_norm_rad(I,j,k) + segment%rx_norm_rad(I,j+1,k))
             enddo
           else
             do J=segment%HI%JsdB,segment%HI%JedB
               dhdt = v_old(i,J,k)-v_new(i,J,k)   !old-new
               dhdx = v_new(i,J,k)-v_new(i-1,J,k) !in new time backward sasha for I-1
               rx_tang_rad(I,J,k) = 0.0
               if (dhdt*dhdx > 0.0) rx_tang_rad(I,J,k) = min( (dhdt/dhdx), rx_max) ! outward phase speed
             enddo
           endif
         enddo
         if (segment%radiation_tan) then
           do k=1,nz ;  do J=segment%HI%JsdB,segment%HI%JedB
             rx_avg = rx_tang_rad(I,J,k)
             segment%tangential_vel(I,J,k) = (v_new(i,J,k) + rx_avg*v_new(i-1,J,k)) / (1.0+rx_avg)
           enddo ; enddo
         endif
         if (segment%nudged_tan) then
           do k=1,nz ; do J=segment%HI%JsdB,segment%HI%JedB
             ! dhdt gets set to 0 on inflow in oblique case
             if (rx_tang_rad(I,J,k) <= 0.0) then
               tau = segment%Velocity_nudging_timescale_in
             else
               tau = segment%Velocity_nudging_timescale_out
             endif
             gamma_2 = dt / (tau + dt)
             segment%tangential_vel(I,J,k) = (1.0 - gamma_2) * segment%tangential_vel(I,J,k) + &
                                 gamma_2 * segment%nudged_tangential_vel(I,J,k)
           enddo ; enddo
         endif
         if (segment%radiation_grad) then
           Js_obc = max(segment%HI%JsdB,G%jsd+1)
           Je_obc = min(segment%HI%JedB,G%jed-1)
           do k=1,nz ; do J=Js_obc,Je_obc
             rx_avg = rx_tang_rad(I,J,k)
!            if (G%mask2dCu(I-1,j) > 0.0 .and. G%mask2dCu(I-1,j+1) > 0.0) then
!              rx_avg = 0.5*(u_new(I-1,j,k) + u_new(I-1,j+1,k)) * dt * G%IdxBu(I-1,J)
!            elseif (G%mask2dCu(I-1,j) > 0.0) then
!              rx_avg = u_new(I-1,j,k) * dt * G%IdxBu(I-1,J)
!            elseif (G%mask2dCu(I-1,j+1) > 0.0) then
!              rx_avg = u_new(I-1,j+1,k) * dt * G%IdxBu(I-1,J)
!            else
!              rx_avg = 0.0
!            endif
             segment%tangential_grad(I,J,k) = ((v_new(i,J,k) - v_new(i-1,J,k))*G%IdxBu(I-1,J) + &
                               rx_avg*(v_new(i-1,J,k) - v_new(i-2,J,k))*G%IdxBu(I-2,J)) / (1.0+rx_avg)
           enddo ; enddo
         endif
         if (segment%nudged_grad) then
           do k=1,nz ; do J=segment%HI%JsdB,segment%HI%JedB
             ! dhdt gets set to 0 on inflow in oblique case
             if (rx_tang_rad(I,J,k) <= 0.0) then
               tau = segment%Velocity_nudging_timescale_in
             else
               tau = segment%Velocity_nudging_timescale_out
             endif
             gamma_2 = dt / (tau + dt)
             segment%tangential_grad(I,J,k) = (1.0 - gamma_2) * segment%tangential_grad(I,J,k) + &
                                 gamma_2 * segment%nudged_tangential_grad(I,J,k)
           enddo ; enddo
         endif
         deallocate(rx_tang_rad)
       endif
       if (segment%oblique_tan .or. segment%oblique_grad) then
         I=segment%HI%IsdB
         allocate(rx_tang_obl(segment%HI%IsdB:segment%HI%IedB,segment%HI%JsdB:segment%HI%JedB,nz))
         allocate(ry_tang_obl(segment%HI%IsdB:segment%HI%IedB,segment%HI%JsdB:segment%HI%JedB,nz))
         allocate(cff_tangential(segment%HI%IsdB:segment%HI%IedB,segment%HI%JsdB:segment%HI%JedB,nz))
         do k=1,nz
           if (gamma_u < 1.0) then
             rx_tang_obl(I,segment%HI%JsdB,k) = segment%rx_norm_obl(I,segment%HI%jsd,k)
             rx_tang_obl(I,segment%HI%JedB,k) = segment%rx_norm_obl(I,segment%HI%jed,k)
             ry_tang_obl(I,segment%HI%JsdB,k) = segment%ry_norm_obl(I,segment%HI%jsd,k)
             ry_tang_obl(I,segment%HI%JedB,k) = segment%ry_norm_obl(I,segment%HI%jed,k)
             cff_tangential(I,segment%HI%JsdB,k) = segment%cff_normal(I,segment%HI%jsd,k)
             cff_tangential(I,segment%HI%JedB,k) = segment%cff_normal(I,segment%HI%jed,k)
             do J=segment%HI%JsdB+1,segment%HI%JedB-1
               rx_tang_obl(I,J,k) = 0.5*(segment%rx_norm_obl(I,j,k) + segment%rx_norm_obl(I,j+1,k))
               ry_tang_obl(I,J,k) = 0.5*(segment%ry_norm_obl(I,j,k) + segment%ry_norm_obl(I,j+1,k))
               cff_tangential(I,J,k) = 0.5*(segment%cff_normal(I,j,k) + segment%cff_normal(I,j+1,k))
             enddo
           else
             do J=segment%HI%JsdB,segment%HI%JedB
               dhdt = v_old(i,J,k)-v_new(i,J,k)   !old-new
               dhdx = v_new(i,J,k)-v_new(i-1,J,k) !in new time backward sasha for I-1
               if (dhdt*(segment%grad_tan(j,1,k) + segment%grad_tan(j+1,1,k)) > 0.0) then
                 dhdy = segment%grad_tan(j,1,k)
               elseif (dhdt*(segment%grad_tan(j,1,k) + segment%grad_tan(j+1,1,k)) == 0.0) then
                 dhdy = 0.0
               else
                 dhdy = segment%grad_tan(j+1,1,k)
               endif
               if (dhdt*dhdx < 0.0) dhdt = 0.0
               rx_new = dhdt*dhdx
               cff_new = max(dhdx*dhdx + dhdy*dhdy, eps)
               ry_new = min(cff_new,max(dhdt*dhdy,-cff_new))
               rx_tang_obl(I,j,k) = rx_new
               ry_tang_obl(i,J,k) = ry_new
               cff_tangential(i,J,k) = cff_new
             enddo
           endif
         enddo
         if (segment%oblique_tan) then
           do k=1,nz ;  do J=segment%HI%JsdB,segment%HI%JedB
             rx_avg = rx_tang_obl(I,J,k)
             ry_avg = ry_tang_obl(I,J,k)
             cff_avg = cff_tangential(I,J,k)
             segment%tangential_vel(I,J,k) = ((cff_avg*v_new(i,J,k) + rx_avg*v_new(i-1,J,k)) - &
                                              (max(ry_avg,0.0)*segment%grad_tan(j,2,k) + &
                                               min(ry_avg,0.0)*segment%grad_tan(j+1,2,k))) / &
                                             (cff_avg + rx_avg)
           enddo ; enddo
         endif
         if (segment%nudged_tan) then
           do k=1,nz ; do J=segment%HI%JsdB,segment%HI%JedB
             ! dhdt gets set to 0 on inflow in oblique case
             if (rx_tang_obl(I,J,k) <= 0.0) then
               tau = segment%Velocity_nudging_timescale_in
             else
               tau = segment%Velocity_nudging_timescale_out
             endif
             gamma_2 = dt / (tau + dt)
             segment%tangential_vel(I,J,k) = (1.0 - gamma_2) * segment%tangential_vel(I,J,k) + &
                                 gamma_2 * segment%nudged_tangential_vel(I,J,k)
           enddo ; enddo
         endif
         if (segment%oblique_grad) then
           Js_obc = max(segment%HI%JsdB,G%jsd+1)
           Je_obc = min(segment%HI%JedB,G%jed-1)
           do k=1,nz ;  do J=segment%HI%JsdB+1,segment%HI%JedB-1
             rx_avg = rx_tang_obl(I,J,k)
             ry_avg = ry_tang_obl(I,J,k)
             cff_avg = cff_tangential(I,J,k)
             segment%tangential_grad(I,J,k) =  &
                 ((cff_avg*(v_new(i,J,k)  - v_new(i-1,J,k))*G%IdxBu(I-1,J) + &
                   rx_avg*(v_new(i-1,J,k) - v_new(i-2,J,k))*G%IdxBu(I-2,J)) - &
                  (max(ry_avg,0.0)*segment%grad_gradient(J,2,k) + &
                   min(ry_avg,0.0)*segment%grad_gradient(J+1,2,k)) ) / &
                 (cff_avg + rx_avg)
           enddo ; enddo
         endif
         if (segment%nudged_grad) then
           do k=1,nz ; do J=segment%HI%JsdB,segment%HI%JedB
             ! dhdt gets set to 0 on inflow in oblique case
             if (rx_tang_obl(I,J,k) <= 0.0) then
               tau = segment%Velocity_nudging_timescale_in
             else
               tau = segment%Velocity_nudging_timescale_out
             endif
             gamma_2 = dt / (tau + dt)
             segment%tangential_grad(I,J,k) = (1.0 - gamma_2) * segment%tangential_grad(I,J,k) + &
                                 gamma_2 * segment%nudged_tangential_grad(I,J,k)
           enddo ; enddo
         endif
         deallocate(rx_tang_obl)
         deallocate(ry_tang_obl)
         deallocate(cff_tangential)
       endif
     endif

     if (segment%direction == OBC_DIRECTION_W) then
       I=segment%HI%IsdB
       if (I>G%HI%IecB) cycle
       do k=1,nz ; do j=segment%HI%jsd,segment%HI%jed
         if (segment%radiation) then
           dhdt = (u_old(I+1,j,k) - u_new(I+1,j,k)) !old-new
           dhdx = (u_new(I+1,j,k) - u_new(I+2,j,k)) !in new time forward sasha for I+1
           rx_new = 0.0
           if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
           if (gamma_u < 1.0) then
             rx_avg = (1.0-gamma_u)*segment%rx_norm_rad(I,j,k) + gamma_u*rx_new
           else
             rx_avg = rx_new
           endif
           segment%rx_norm_rad(I,j,k) = rx_avg
           ! The new boundary value is interpolated between future interior
           ! value, u_new(I+1) and past boundary value but with barotropic
           ! accelerations, u_new(I).
           segment%normal_vel(I,j,k) = (u_new(I,j,k) + rx_avg*u_new(I+1,j,k)) / (1.0+rx_avg)
           if (gamma_u < 1.0) then
             ! Copy restart fields into 3-d arrays. This is an inefficient and temporary issues
             ! implemented as a work-around to limitations in restart capability
             OBC%rx_normal(I,j,k) = segment%rx_norm_rad(I,j,k)
           endif
         elseif (segment%oblique) then
           dhdt = (u_old(I+1,j,k) - u_new(I+1,j,k)) !old-new
           dhdx = (u_new(I+1,j,k) - u_new(I+2,j,k)) !in new time forward sasha for I+1
           if (dhdt*(segment%grad_normal(J,1,k) + segment%grad_normal(J-1,1,k)) > 0.0) then
             dhdy = segment%grad_normal(J-1,1,k)
           elseif (dhdt*(segment%grad_normal(J,1,k) + segment%grad_normal(J-1,1,k)) == 0.0) then
             dhdy = 0.0
           else
             dhdy = segment%grad_normal(J,1,k)
           endif
           if (dhdt*dhdx < 0.0) dhdt = 0.0

           rx_new = dhdt*dhdx
           cff_new = max(dhdx*dhdx + dhdy*dhdy, eps)
           ry_new = min(cff_new,max(dhdt*dhdy,-cff_new))
           if (gamma_u < 1.0) then
             rx_avg = (1.0-gamma_u)*segment%rx_norm_obl(I,j,k) + gamma_u*rx_new
             ry_avg = (1.0-gamma_u)*segment%ry_norm_obl(i,J,k) + gamma_u*ry_new
             cff_avg = (1.0-gamma_u)*segment%cff_normal(I,j,k) + gamma_u*cff_new
           else
             rx_avg = rx_new
             ry_avg = ry_new
             cff_avg = cff_new
           endif
           segment%rx_norm_obl(I,j,k) = rx_avg
           segment%ry_norm_obl(i,J,k) = ry_avg
           segment%cff_normal(i,J,k) = cff_avg
           segment%normal_vel(I,j,k) = ((cff_avg*u_new(I,j,k) + rx_avg*u_new(I+1,j,k)) - &
                                        (max(ry_avg,0.0)*segment%grad_normal(J-1,2,k) + &
                                         min(ry_avg,0.0)*segment%grad_normal(J,2,k))) / &
                                       (cff_avg + rx_avg)
           if (gamma_u < 1.0) then
             ! Copy restart fields into 3-d arrays. This is an inefficient and temporary issues
             ! implemented as a work-around to limitations in restart capability
             OBC%rx_oblique(I,j,k) = segment%rx_norm_obl(I,j,k)
             OBC%ry_oblique(i,J,k) = segment%ry_norm_obl(i,J,k)
             OBC%cff_normal(I,j,k) = segment%cff_normal(I,j,k)
           endif
         elseif (segment%gradient) then
           segment%normal_vel(I,j,k) = u_new(I+1,j,k)
         endif
         if ((segment%radiation .or. segment%oblique) .and. segment%nudged) then
           ! dhdt gets set to 0. on inflow in oblique case
           if (dhdt*dhdx <= 0.0) then
             tau = segment%Velocity_nudging_timescale_in
           else
             tau = segment%Velocity_nudging_timescale_out
           endif
           gamma_2 = dt / (tau + dt)
           segment%normal_vel(I,j,k) = (1.0 - gamma_2) * segment%normal_vel(I,j,k) + &
                                 gamma_2 * segment%nudged_normal_vel(I,j,k)
         endif
       enddo ; enddo
       if (segment%radiation_tan .or. segment%radiation_grad) then
         I=segment%HI%IsdB
         allocate(rx_tang_rad(segment%HI%IsdB:segment%HI%IedB,segment%HI%JsdB:segment%HI%JedB,nz))
         do k=1,nz
           if (gamma_u < 1.0) then
             rx_tang_rad(I,segment%HI%JsdB,k) = segment%rx_norm_rad(I,segment%HI%jsd,k)
             rx_tang_rad(I,segment%HI%JedB,k) = segment%rx_norm_rad(I,segment%HI%jed,k)
             do J=segment%HI%JsdB+1,segment%HI%JedB-1
               rx_tang_rad(I,J,k) = 0.5*(segment%rx_norm_rad(I,j,k) + segment%rx_norm_rad(I,j+1,k))
             enddo
           else
             do J=segment%HI%JsdB,segment%HI%JedB
               dhdt = v_old(i+1,J,k)-v_new(i+1,J,k)   !old-new
               dhdx = v_new(i+1,J,k)-v_new(i+2,J,k) !in new time backward sasha for I-1
               rx_tang_rad(I,J,k) = 0.0
               if (dhdt*dhdx > 0.0) rx_tang_rad(I,J,k) = min( (dhdt/dhdx), rx_max) ! outward phase speed
             enddo
           endif
         enddo
         if (segment%radiation_tan) then
           do k=1,nz ;  do J=segment%HI%JsdB,segment%HI%JedB
             rx_avg = rx_tang_rad(I,J,k)
             segment%tangential_vel(I,J,k) = (v_new(i+1,J,k) + rx_avg*v_new(i+2,J,k)) / (1.0+rx_avg)
           enddo ; enddo
         endif
         if (segment%nudged_tan) then
           do k=1,nz ; do J=segment%HI%JsdB,segment%HI%JedB
             ! dhdt gets set to 0 on inflow in oblique case
             if (rx_tang_rad(I,J,k) <= 0.0) then
               tau = segment%Velocity_nudging_timescale_in
             else
               tau = segment%Velocity_nudging_timescale_out
             endif
             gamma_2 = dt / (tau + dt)
             segment%tangential_vel(I,J,k) = (1.0 - gamma_2) * segment%tangential_vel(I,J,k) + &
                                 gamma_2 * segment%nudged_tangential_vel(I,J,k)
           enddo ; enddo
         endif
         if (segment%radiation_grad) then
           Js_obc = max(segment%HI%JsdB,G%jsd+1)
           Je_obc = min(segment%HI%JedB,G%jed-1)
           do k=1,nz ;  do J=Js_obc,Je_obc
             rx_avg = rx_tang_rad(I,J,k)
!            if (G%mask2dCu(I+1,j) > 0.0 .and. G%mask2dCu(I+1,j+1) > 0.0) then
!              rx_avg = 0.5*(u_new(I+1,j,k) + u_new(I+1,j+1,k)) * dt * G%IdxBu(I+1,J)
!            elseif (G%mask2dCu(I+1,j) > 0.0) then
!              rx_avg = u_new(I+1,j,k) * dt * G%IdxBu(I+1,J)
!            elseif (G%mask2dCu(I+1,j+1) > 0.0) then
!              rx_avg = u_new(I+1,j+1,k) * dt * G%IdxBu(I+1,J)
!            else
!              rx_avg = 0.0
!            endif
             segment%tangential_grad(I,J,k) = ((v_new(i+2,J,k) - v_new(i+1,J,k))*G%IdxBu(I+1,J) + &
                               rx_avg*(v_new(i+3,J,k) - v_new(i+2,J,k))*G%IdxBu(I+2,J)) / (1.0+rx_avg)
           enddo ; enddo
         endif
         if (segment%nudged_grad) then
           do k=1,nz ; do J=segment%HI%JsdB,segment%HI%JedB
             ! dhdt gets set to 0 on inflow in oblique case
             if (rx_tang_rad(I,J,k) <= 0.0) then
               tau = segment%Velocity_nudging_timescale_in
             else
               tau = segment%Velocity_nudging_timescale_out
             endif
             gamma_2 = dt / (tau + dt)
             segment%tangential_grad(I,J,k) = (1.0 - gamma_2) * segment%tangential_grad(I,J,k) + &
                                 gamma_2 * segment%nudged_tangential_grad(I,J,k)
           enddo ; enddo
         endif
         deallocate(rx_tang_rad)
       endif
       if (segment%oblique_tan .or. segment%oblique_grad) then
         I=segment%HI%IsdB
         allocate(rx_tang_obl(segment%HI%IsdB:segment%HI%IedB,segment%HI%JsdB:segment%HI%JedB,nz))
         allocate(ry_tang_obl(segment%HI%IsdB:segment%HI%IedB,segment%HI%JsdB:segment%HI%JedB,nz))
         allocate(cff_tangential(segment%HI%IsdB:segment%HI%IedB,segment%HI%JsdB:segment%HI%JedB,nz))
         do k=1,nz
           if (gamma_u < 1.0) then
             rx_tang_obl(I,segment%HI%JsdB,k) = segment%rx_norm_obl(I,segment%HI%jsd,k)
             rx_tang_obl(I,segment%HI%JedB,k) = segment%rx_norm_obl(I,segment%HI%jed,k)
             ry_tang_obl(I,segment%HI%JsdB,k) = segment%ry_norm_obl(I,segment%HI%jsd,k)
             ry_tang_obl(I,segment%HI%JedB,k) = segment%ry_norm_obl(I,segment%HI%jed,k)
             cff_tangential(I,segment%HI%JsdB,k) = segment%cff_normal(I,segment%HI%jsd,k)
             cff_tangential(I,segment%HI%JedB,k) = segment%cff_normal(I,segment%HI%jed,k)
             do J=segment%HI%JsdB+1,segment%HI%JedB-1
               rx_tang_obl(I,J,k) = 0.5*(segment%rx_norm_obl(I,j,k) + segment%rx_norm_obl(I,j+1,k))
               ry_tang_obl(I,J,k) = 0.5*(segment%ry_norm_obl(I,j,k) + segment%ry_norm_obl(I,j+1,k))
               cff_tangential(I,J,k) = 0.5*(segment%cff_normal(I,j,k) + segment%cff_normal(I,j+1,k))
             enddo
           else
             do J=segment%HI%JsdB,segment%HI%JedB
               dhdt = v_old(i+1,J,k)-v_new(i+1,J,k)   !old-new
               dhdx = v_new(i+1,J,k)-v_new(i+2,J,k) !in new time backward sasha for I-1
               if (dhdt*(segment%grad_tan(j,1,k) + segment%grad_tan(j+1,1,k)) > 0.0) then
                 dhdy = segment%grad_tan(j,1,k)
               elseif (dhdt*(segment%grad_tan(j,1,k) + segment%grad_tan(j+1,1,k)) == 0.0) then
                 dhdy = 0.0
               else
                 dhdy = segment%grad_tan(j+1,1,k)
               endif
               if (dhdt*dhdx < 0.0) dhdt = 0.0
               rx_new = dhdt*dhdx
               cff_new = max(dhdx*dhdx + dhdy*dhdy, eps)
               ry_new = min(cff_new,max(dhdt*dhdy,-cff_new))
               rx_tang_obl(I,j,k) = rx_new
               ry_tang_obl(i,J,k) = ry_new
               cff_tangential(i,J,k) = cff_new
             enddo
           endif
         enddo
         if (segment%oblique_tan) then
           do k=1,nz ;  do J=segment%HI%JsdB,segment%HI%JedB
             rx_avg = rx_tang_obl(I,J,k)
             ry_avg = ry_tang_obl(I,J,k)
             cff_avg = cff_tangential(I,J,k)
             segment%tangential_vel(I,J,k) = ((cff_avg*v_new(i+1,J,k) + rx_avg*v_new(i+2,J,k)) - &
                                              (max(ry_avg,0.0)*segment%grad_tan(j,2,k) + &
                                               min(ry_avg,0.0)*segment%grad_tan(j+1,2,k))) / &
                                             (cff_avg + rx_avg)
           enddo ; enddo
         endif
         if (segment%nudged_tan) then
           do k=1,nz ; do J=segment%HI%JsdB,segment%HI%JedB
             ! dhdt gets set to 0 on inflow in oblique case
             if (rx_tang_obl(I,J,k) <= 0.0) then
               tau = segment%Velocity_nudging_timescale_in
             else
               tau = segment%Velocity_nudging_timescale_out
             endif
             gamma_2 = dt / (tau + dt)
             segment%tangential_vel(I,J,k) = (1.0 - gamma_2) * segment%tangential_vel(I,J,k) + &
                                 gamma_2 * segment%nudged_tangential_vel(I,J,k)
           enddo ; enddo
         endif
         if (segment%oblique_grad) then
           Js_obc = max(segment%HI%JsdB,G%jsd+1)
           Je_obc = min(segment%HI%JedB,G%jed-1)
           do k=1,nz ;  do J=segment%HI%JsdB+1,segment%HI%JedB-1
             rx_avg = rx_tang_obl(I,J,k)
             ry_avg = ry_tang_obl(I,J,k)
             cff_avg = cff_tangential(I,J,k)
             segment%tangential_grad(I,J,k) = &
                 ((cff_avg*(v_new(i+2,J,k) - v_new(i+1,J,k))*G%IdxBu(I+1,J) + &
                    rx_avg*(v_new(i+3,J,k) - v_new(i+2,J,k))*G%IdxBu(I+2,J)) - &
                  (max(ry_avg,0.0)*segment%grad_gradient(J,2,k) + &
                   min(ry_avg,0.0)*segment%grad_gradient(J+1,2,k))) / &
                 (cff_avg + rx_avg)
           enddo ; enddo
         endif
         if (segment%nudged_grad) then
           do k=1,nz ; do J=segment%HI%JsdB,segment%HI%JedB
             ! dhdt gets set to 0 on inflow in oblique case
             if (rx_tang_obl(I,J,k) <= 0.0) then
               tau = segment%Velocity_nudging_timescale_in
             else
               tau = segment%Velocity_nudging_timescale_out
             endif
             gamma_2 = dt / (tau + dt)
             segment%tangential_grad(I,J,k) = (1.0 - gamma_2) * segment%tangential_grad(I,J,k) + &
                                 gamma_2 * segment%nudged_tangential_grad(I,J,k)
           enddo ; enddo
         endif
         deallocate(rx_tang_obl)
         deallocate(ry_tang_obl)
         deallocate(cff_tangential)
       endif
     endif

     if (segment%direction == OBC_DIRECTION_N) then
       J=segment%HI%JsdB
       if (J<G%HI%JscB) cycle
       do k=1,nz ;  do i=segment%HI%isd,segment%HI%ied
         if (segment%radiation) then
           dhdt = (v_old(i,J-1,k) - v_new(i,J-1,k)) !old-new
           dhdy = (v_new(i,J-1,k) - v_new(i,J-2,k)) !in new time backward sasha for J-1
           ry_new = 0.0
           if (dhdt*dhdy > 0.0) ry_new = min( (dhdt/dhdy), ry_max)
           if (gamma_u < 1.0) then
             ry_avg = (1.0-gamma_u)*segment%ry_norm_rad(I,j,k) + gamma_u*ry_new
           else
             ry_avg = ry_new
           endif
           segment%ry_norm_rad(i,J,k) = ry_avg
           ! The new boundary value is interpolated between future interior
           ! value, v_new(J-1) and past boundary value but with barotropic
           ! accelerations, v_new(J).
           segment%normal_vel(i,J,k) = (v_new(i,J,k) + ry_avg*v_new(i,J-1,k)) / (1.0+ry_avg)
           if (gamma_u < 1.0) then
             ! Copy restart fields into 3-d arrays. This is an inefficient and temporary issues
             ! implemented as a work-around to limitations in restart capability
             OBC%ry_normal(i,J,k) = segment%ry_norm_rad(i,J,k)
           endif
         elseif (segment%oblique) then
           dhdt = (v_old(i,J-1,k) - v_new(i,J-1,k)) !old-new
           dhdy = (v_new(i,J-1,k) - v_new(i,J-2,k)) !in new time backward sasha for J-1
           if (dhdt*(segment%grad_normal(I,1,k) + segment%grad_normal(I-1,1,k)) > 0.0) then
             dhdx = segment%grad_normal(I-1,1,k)
           elseif (dhdt*(segment%grad_normal(I,1,k) + segment%grad_normal(I-1,1,k)) == 0.0) then
             dhdx = 0.0
           else
             dhdx = segment%grad_normal(I,1,k)
           endif
           if (dhdt*dhdy < 0.0) dhdt = 0.0
           ry_new = dhdt*dhdy
           cff_new = max(dhdx*dhdx + dhdy*dhdy, eps)
           rx_new = min(cff_new,max(dhdt*dhdx,-cff_new))
           if (gamma_u < 1.0) then
             rx_avg = (1.0-gamma_u)*segment%rx_norm_obl(I,j,k) + gamma_u*rx_new
             ry_avg = (1.0-gamma_u)*segment%ry_norm_obl(i,J,k) + gamma_u*ry_new
             cff_avg = (1.0-gamma_u)*segment%cff_normal(i,J,k) + gamma_u*cff_new
           else
             rx_avg = rx_new
             ry_avg = ry_new
             cff_avg = cff_new
           endif
           segment%rx_norm_obl(I,j,k) = rx_avg
           segment%ry_norm_obl(i,J,k) = ry_avg
           segment%cff_normal(i,J,k) = cff_avg
           segment%normal_vel(i,J,k) = ((cff_avg*v_new(i,J,k) + ry_avg*v_new(i,J-1,k)) - &
                                        (max(rx_avg,0.0)*segment%grad_normal(I-1,2,k) +&
                                         min(rx_avg,0.0)*segment%grad_normal(I,2,k))) / &
                                       (cff_avg + ry_avg)
           if (gamma_u < 1.0) then
             ! Copy restart fields into 3-d arrays. This is an inefficient and temporary issues
             ! implemented as a work-around to limitations in restart capability
             OBC%rx_oblique(I,j,k) = segment%rx_norm_obl(I,j,k)
             OBC%ry_oblique(i,J,k) = segment%ry_norm_obl(i,J,k)
             OBC%cff_normal(i,J,k) = segment%cff_normal(i,J,k)
           endif
         elseif (segment%gradient) then
           segment%normal_vel(i,J,k) = v_new(i,J-1,k)
         endif
         if ((segment%radiation .or. segment%oblique) .and. segment%nudged) then
           ! dhdt gets set to 0 on inflow in oblique case
           if (dhdt*dhdy <= 0.0) then
             tau = segment%Velocity_nudging_timescale_in
           else
             tau = segment%Velocity_nudging_timescale_out
           endif
           gamma_2 = dt / (tau + dt)
           segment%normal_vel(i,J,k) = (1.0 - gamma_2) * segment%normal_vel(i,J,k) + &
                                 gamma_2 * segment%nudged_normal_vel(i,J,k)
         endif
       enddo ; enddo
       if (segment%radiation_tan .or. segment%radiation_grad) then
         J=segment%HI%JsdB
         allocate(ry_tang_rad(segment%HI%IsdB:segment%HI%IedB,segment%HI%JsdB:segment%HI%JedB,nz))
         do k=1,nz
           if (gamma_u < 1.0) then
             ry_tang_rad(segment%HI%IsdB,J,k) = segment%ry_norm_rad(segment%HI%isd,J,k)
             ry_tang_rad(segment%HI%IedB,J,k) = segment%ry_norm_rad(segment%HI%ied,J,k)
             do I=segment%HI%IsdB+1,segment%HI%IedB-1
               ry_tang_rad(I,J,k) = 0.5*(segment%ry_norm_rad(i,J,k) + segment%ry_norm_rad(i+1,J,k))
             enddo
           else
             do I=segment%HI%IsdB,segment%HI%IedB
               dhdt = u_old(I,j-1,k)-u_new(I,j-1,k)   !old-new
               dhdy = u_new(I,j-1,k)-u_new(I,j-2,k) !in new time backward sasha for I-1
               ry_tang_rad(I,J,k) = 0.0
               if (dhdt*dhdy > 0.0) ry_tang_rad(I,J,k) = min( (dhdt/dhdy), rx_max) ! outward phase speed
             enddo
           endif
         enddo
         if (segment%radiation_tan) then
           do k=1,nz ;  do I=segment%HI%IsdB,segment%HI%IedB
             ry_avg = ry_tang_rad(I,J,k)
             segment%tangential_vel(I,J,k) = (u_new(I,j,k) + ry_avg*u_new(I,j-1,k)) / (1.0+ry_avg)
           enddo ; enddo
         endif
         if (segment%nudged_tan) then
           do k=1,nz ; do I=segment%HI%IsdB,segment%HI%IedB
             ! dhdt gets set to 0 on inflow in oblique case
             if (ry_tang_rad(I,J,k) <= 0.0) then
               tau = segment%Velocity_nudging_timescale_in
             else
               tau = segment%Velocity_nudging_timescale_out
             endif
             gamma_2 = dt / (tau + dt)
             segment%tangential_vel(I,J,k) = (1.0 - gamma_2) * segment%tangential_vel(I,J,k) + &
                                 gamma_2 * segment%nudged_tangential_vel(I,J,k)
           enddo ; enddo
         endif
         if (segment%radiation_grad) then
           Is_obc = max(segment%HI%IsdB,G%isd+1)
           Ie_obc = min(segment%HI%IedB,G%ied-1)
           do k=1,nz ;  do I=Is_obc,Ie_obc
             ry_avg = ry_tang_rad(I,J,k)
!            if (G%mask2dCv(i,J-1) > 0.0 .and. G%mask2dCv(i+1,J-1) > 0.0) then
!              ry_avg = 0.5*(v_new(i,J-1,k) + v_new(i+1,J-1,k) * dt * G%IdyBu(I,J-1))
!            elseif (G%mask2dCv(i,J-1) > 0.0) then
!              ry_avg = v_new(i,J-1,k) * dt *G%IdyBu(I,J-1)
!            elseif (G%mask2dCv(i+1,J-1) > 0.0) then
!              ry_avg = v_new(i+1,J-1,k) * dt *G%IdyBu(I,J-1)
!            else
!              ry_avg = 0.0
!            endif
             segment%tangential_grad(I,J,k) = ((u_new(I,j,k) - u_new(I,j-1,k))*G%IdyBu(I,J-1) + &
                               ry_avg*(u_new(I,j-1,k) - u_new(I,j-2,k))*G%IdyBu(I,J-2)) / (1.0+ry_avg)
           enddo ; enddo
         endif
         if (segment%nudged_grad) then
           do k=1,nz ; do I=segment%HI%IsdB,segment%HI%IedB
             ! dhdt gets set to 0 on inflow in oblique case
             if (ry_tang_rad(I,J,k) <= 0.0) then
               tau = segment%Velocity_nudging_timescale_in
             else
               tau = segment%Velocity_nudging_timescale_out
             endif
             gamma_2 = dt / (tau + dt)
             segment%tangential_grad(I,J,k) = (1.0 - gamma_2) * segment%tangential_grad(I,J,k) + &
                                 gamma_2 * segment%nudged_tangential_grad(I,J,k)
           enddo ; enddo
         endif
         deallocate(ry_tang_rad)
       endif
       if (segment%oblique_tan .or. segment%oblique_grad) then
         J=segment%HI%JsdB
         allocate(rx_tang_obl(segment%HI%IsdB:segment%HI%IedB,segment%HI%JsdB:segment%HI%JedB,nz))
         allocate(ry_tang_obl(segment%HI%IsdB:segment%HI%IedB,segment%HI%JsdB:segment%HI%JedB,nz))
         allocate(cff_tangential(segment%HI%IsdB:segment%HI%IedB,segment%HI%JsdB:segment%HI%JedB,nz))
         do k=1,nz
           if (gamma_u < 1.0) then
             rx_tang_obl(segment%HI%IsdB,J,k) = segment%rx_norm_obl(segment%HI%isd,J,k)
             rx_tang_obl(segment%HI%IedB,J,k) = segment%rx_norm_obl(segment%HI%ied,J,k)
             ry_tang_obl(segment%HI%IsdB,J,k) = segment%ry_norm_obl(segment%HI%isd,J,k)
             ry_tang_obl(segment%HI%IedB,J,k) = segment%ry_norm_obl(segment%HI%ied,J,k)
             cff_tangential(segment%HI%IsdB,J,k) = segment%cff_normal(segment%HI%isd,J,k)
             cff_tangential(segment%HI%IedB,J,k) = segment%cff_normal(segment%HI%ied,J,k)
             do I=segment%HI%IsdB+1,segment%HI%IedB-1
               rx_tang_obl(I,J,k) = 0.5*(segment%rx_norm_obl(i,J,k) + segment%rx_norm_obl(i+1,J,k))
               ry_tang_obl(I,J,k) = 0.5*(segment%ry_norm_obl(i,J,k) + segment%ry_norm_obl(i+1,J,k))
               cff_tangential(I,J,k) = 0.5*(segment%cff_normal(i,J,k) + segment%cff_normal(i+1,J,k))
             enddo
           else
             do I=segment%HI%IsdB,segment%HI%IedB
               dhdt = u_old(I,j,k)-u_new(I,j,k)   !old-new
               dhdy = u_new(I,j,k)-u_new(I,j-1,k) !in new time backward sasha for I-1
               if (dhdt*(segment%grad_tan(i,1,k) + segment%grad_tan(i+1,1,k)) > 0.0) then
                 dhdx = segment%grad_tan(i,1,k)
               elseif (dhdt*(segment%grad_tan(i,1,k) + segment%grad_tan(i+1,1,k)) == 0.0) then
                 dhdx = 0.0
               else
                 dhdx = segment%grad_tan(i+1,1,k)
               endif
               if (dhdt*dhdy < 0.0) dhdt = 0.0
               ry_new = dhdt*dhdy
               cff_new = max(dhdx*dhdx + dhdy*dhdy, eps)
               rx_new = min(cff_new,max(dhdt*dhdx,-cff_new))
               rx_tang_obl(I,j,k) = rx_new
               ry_tang_obl(i,J,k) = ry_new
               cff_tangential(i,J,k) = cff_new
             enddo
           endif
         enddo
         if (segment%oblique_tan) then
           do k=1,nz ;  do I=segment%HI%IsdB,segment%HI%IedB
             rx_avg = rx_tang_obl(I,J,k)
             ry_avg = ry_tang_obl(I,J,k)
             cff_avg = cff_tangential(I,J,k)
             segment%tangential_vel(I,J,k) = ((cff_avg*u_new(I,j,k) + ry_avg*u_new(I,j-1,k)) - &
                                              (max(rx_avg,0.0)*segment%grad_tan(i,2,k) + &
                                               min(rx_avg,0.0)*segment%grad_tan(i+1,2,k))) / &
                                             (cff_avg + ry_avg)
           enddo ; enddo
         endif
         if (segment%nudged_tan) then
           do k=1,nz ; do I=segment%HI%IsdB,segment%HI%IedB
             ! dhdt gets set to 0 on inflow in oblique case
             if (ry_tang_obl(I,J,k) <= 0.0) then
               tau = segment%Velocity_nudging_timescale_in
             else
               tau = segment%Velocity_nudging_timescale_out
             endif
             gamma_2 = dt / (tau + dt)
             segment%tangential_vel(I,J,k) = (1.0 - gamma_2) * segment%tangential_vel(I,J,k) + &
                                 gamma_2 * segment%nudged_tangential_vel(I,J,k)
           enddo ; enddo
         endif
         if (segment%oblique_grad) then
           Is_obc = max(segment%HI%IsdB,G%isd+1)
           Ie_obc = min(segment%HI%IedB,G%ied-1)
           do k=1,nz ;  do I=segment%HI%IsdB+1,segment%HI%IedB-1
             rx_avg = rx_tang_obl(I,J,k)
             ry_avg = ry_tang_obl(I,J,k)
             cff_avg = cff_tangential(I,J,k)
             segment%tangential_grad(I,J,k) =  &
                 ((cff_avg*(u_new(I,j,k)   - u_new(I,j-1,k))*G%IdyBu(I,J-1) + &
                    ry_avg*(u_new(I,j-1,k) - u_new(I,j-2,k))*G%IdyBu(I,J-2)) - &
                                  (max(rx_avg,0.0)*segment%grad_gradient(I,2,k) + &
                                   min(rx_avg,0.0)*segment%grad_gradient(I+1,2,k))) / &
                                 (cff_avg + ry_avg)
           enddo ; enddo
         endif
         if (segment%nudged_grad) then
           do k=1,nz ; do I=segment%HI%IsdB,segment%HI%IedB
             ! dhdt gets set to 0 on inflow in oblique case
             if (ry_tang_obl(I,J,k) <= 0.0) then
               tau = segment%Velocity_nudging_timescale_in
             else
               tau = segment%Velocity_nudging_timescale_out
             endif
             gamma_2 = dt / (tau + dt)
             segment%tangential_grad(I,J,k) = (1.0 - gamma_2) * segment%tangential_grad(I,J,k) + &
                                 gamma_2 * segment%nudged_tangential_grad(I,J,k)
           enddo ; enddo
         endif
         deallocate(rx_tang_obl)
         deallocate(ry_tang_obl)
         deallocate(cff_tangential)
       endif
     endif

     if (segment%direction == OBC_DIRECTION_S) then
       J=segment%HI%JsdB
       if (J>G%HI%JecB) cycle
       do k=1,nz ;  do i=segment%HI%isd,segment%HI%ied
         if (segment%radiation) then
           dhdt = (v_old(i,J+1,k) - v_new(i,J+1,k)) !old-new
           dhdy = (v_new(i,J+1,k) - v_new(i,J+2,k)) !in new time backward sasha for J-1
           ry_new = 0.0
           if (dhdt*dhdy > 0.0) ry_new = min( (dhdt/dhdy), ry_max)
           if (gamma_u < 1.0) then
             ry_avg = (1.0-gamma_u)*segment%ry_norm_rad(I,j,k) + gamma_u*ry_new
           else
             ry_avg = ry_new
           endif
           segment%ry_norm_rad(i,J,k) = ry_avg
           ! The new boundary value is interpolated between future interior
           ! value, v_new(J+1) and past boundary value but with barotropic
           ! accelerations, v_new(J).
           segment%normal_vel(i,J,k) = (v_new(i,J,k) + ry_avg*v_new(i,J+1,k)) / (1.0+ry_avg)
           if (gamma_u < 1.0) then
             ! Copy restart fields into 3-d arrays. This is an inefficient and temporary issues
             ! implemented as a work-around to limitations in restart capability
             OBC%ry_normal(i,J,k) = segment%ry_norm_rad(i,J,k)
           endif
         elseif (segment%oblique) then
           dhdt = (v_old(i,J+1,k) - v_new(i,J+1,k)) !old-new
           dhdy = (v_new(i,J+1,k) - v_new(i,J+2,k)) !in new time backward sasha for J-1
           if (dhdt*(segment%grad_normal(I,1,k) + segment%grad_normal(I-1,1,k)) > 0.0) then
             dhdx = segment%grad_normal(I-1,1,k)
           elseif (dhdt*(segment%grad_normal(I,1,k) + segment%grad_normal(I-1,1,k)) == 0.0) then
             dhdx = 0.0
           else
             dhdx = segment%grad_normal(I,1,k)
           endif
           if (dhdt*dhdy < 0.0) dhdt = 0.0

           ry_new = dhdt*dhdy
           cff_new = max(dhdx*dhdx + dhdy*dhdy, eps)
           rx_new = min(cff_new,max(dhdt*dhdx,-cff_new))
           if (gamma_u < 1.0) then
             rx_avg = (1.0-gamma_u)*segment%rx_norm_obl(I,j,k) + gamma_u*rx_new
             ry_avg = (1.0-gamma_u)*segment%ry_norm_obl(i,J,k) + gamma_u*ry_new
             cff_avg = (1.0-gamma_u)*segment%cff_normal(i,J,k) + gamma_u*cff_new
           else
             rx_avg = rx_new
             ry_avg = ry_new
             cff_avg = cff_new
           endif
           segment%rx_norm_obl(I,j,k) = rx_avg
           segment%ry_norm_obl(i,J,k) = ry_avg
           segment%cff_normal(i,J,k) = cff_avg
           segment%normal_vel(i,J,k) = ((cff_avg*v_new(i,J,k) + ry_avg*v_new(i,J+1,k)) - &
                                        (max(rx_avg,0.0)*segment%grad_normal(I-1,2,k) + &
                                         min(rx_avg,0.0)*segment%grad_normal(I,2,k))) / &
                                       (cff_avg + ry_avg)
           if (gamma_u < 1.0) then
             ! Copy restart fields into 3-d arrays. This is an inefficient and temporary issues
             ! implemented as a work-around to limitations in restart capability
             OBC%rx_oblique(I,j,k) = segment%rx_norm_obl(I,j,k)
             OBC%ry_oblique(i,J,k) = segment%ry_norm_obl(i,J,k)
             OBC%cff_normal(i,J,k) = segment%cff_normal(i,J,k)
           endif
         elseif (segment%gradient) then
           segment%normal_vel(i,J,k) = v_new(i,J+1,k)
         endif
         if ((segment%radiation .or. segment%oblique) .and. segment%nudged) then
           ! dhdt gets set to 0 on inflow in oblique case
           if (dhdt*dhdy <= 0.0) then
             tau = segment%Velocity_nudging_timescale_in
           else
             tau = segment%Velocity_nudging_timescale_out
           endif
           gamma_2 = dt / (tau + dt)
           segment%normal_vel(i,J,k) = (1.0 - gamma_2) * segment%normal_vel(i,J,k) + &
                                 gamma_2 * segment%nudged_normal_vel(i,J,k)
         endif
       enddo ; enddo
       if (segment%radiation_tan .or. segment%radiation_grad) then
         J=segment%HI%JsdB
         allocate(ry_tang_rad(segment%HI%IsdB:segment%HI%IedB,segment%HI%JsdB:segment%HI%JedB,nz))
         do k=1,nz
           if (gamma_u < 1.0) then
             ry_tang_rad(segment%HI%IsdB,J,k) = segment%ry_norm_rad(segment%HI%isd,J,k)
             ry_tang_rad(segment%HI%IedB,J,k) = segment%ry_norm_rad(segment%HI%ied,J,k)
             do I=segment%HI%IsdB+1,segment%HI%IedB-1
               ry_tang_rad(I,J,k) = 0.5*(segment%ry_norm_rad(i,J,k) + segment%ry_norm_rad(i+1,J,k))
             enddo
           else
             do I=segment%HI%IsdB,segment%HI%IedB
               dhdt = u_old(I,j+1,k)-u_new(I,j+1,k)   !old-new
               dhdy = u_new(I,j+1,k)-u_new(I,j+2,k) !in new time backward sasha for I-1
               ry_tang_rad(I,J,k) = 0.0
               if (dhdt*dhdy > 0.0) ry_tang_rad(I,J,k) = min( (dhdt/dhdy), rx_max) ! outward phase speed
             enddo
           endif
         enddo
         if (segment%radiation_tan) then
           do k=1,nz ;  do I=segment%HI%IsdB,segment%HI%IedB
             ry_avg = ry_tang_rad(I,J,k)
             segment%tangential_vel(I,J,k) = (u_new(I,j+1,k) + ry_avg*u_new(I,j+2,k)) / (1.0+ry_avg)
           enddo ; enddo
         endif
         if (segment%nudged_tan) then
           do k=1,nz ; do I=segment%HI%IsdB,segment%HI%IedB
             ! dhdt gets set to 0 on inflow in oblique case
             if (ry_tang_rad(I,J,k) <= 0.0) then
               tau = segment%Velocity_nudging_timescale_in
             else
               tau = segment%Velocity_nudging_timescale_out
             endif
             gamma_2 = dt / (tau + dt)
             segment%tangential_vel(I,J,k) = (1.0 - gamma_2) * segment%tangential_vel(I,J,k) + &
                                 gamma_2 * segment%nudged_tangential_vel(I,J,k)
           enddo ; enddo
         endif
         if (segment%radiation_grad) then
           Is_obc = max(segment%HI%IsdB,G%isd+1)
           Ie_obc = min(segment%HI%IedB,G%ied-1)
           do k=1,nz ;  do I=Is_obc,Ie_obc
             ry_avg = ry_tang_rad(I,J,k)
!            if (G%mask2dCv(i,J+1) > 0.0 .and. G%mask2dCv(i+1,J+1) > 0.0) then
!              ry_avg = 0.5*(v_new(i,J+1,k) + v_new(i+1,J+1,k)) * dt * G%IdyBu(I,J+1)
!            elseif (G%mask2dCv(i,J+1) > 0.0) then
!              ry_avg = v_new(i,J+1,k) * dt * G%IdyBu(I,J+1)
!            elseif (G%mask2dCv(i+1,J+1) > 0.0) then
!              ry_avg = v_new(i+1,J+1,k) * dt * G%IdyBu(I,J+1)
!            else
!              ry_avg = 0.0
!            endif
             segment%tangential_grad(I,J,k) = ((u_new(I,j+2,k) - u_new(I,j+1,k))*G%IdyBu(I,J+1) + &
                               ry_avg*(u_new(I,j+3,k) - u_new(I,j+2,k))*G%IdyBu(I,J+2)) / (1.0+ry_avg)
           enddo ; enddo
         endif
         if (segment%nudged_grad) then
           do k=1,nz ; do I=segment%HI%IsdB,segment%HI%IedB
             ! dhdt gets set to 0 on inflow in oblique case
             if (ry_tang_rad(I,J,k) <= 0.0) then
               tau = segment%Velocity_nudging_timescale_in
             else
               tau = segment%Velocity_nudging_timescale_out
             endif
             gamma_2 = dt / (tau + dt)
             segment%tangential_grad(I,J,k) = (1.0 - gamma_2) * segment%tangential_grad(I,J,k) + &
                                 gamma_2 * segment%nudged_tangential_grad(I,J,k)
           enddo ; enddo
         endif
         deallocate(ry_tang_rad)
       endif
       if (segment%oblique_tan .or. segment%oblique_grad) then
         J=segment%HI%JsdB
         allocate(rx_tang_obl(segment%HI%IsdB:segment%HI%IedB,segment%HI%JsdB:segment%HI%JedB,nz))
         allocate(ry_tang_obl(segment%HI%IsdB:segment%HI%IedB,segment%HI%JsdB:segment%HI%JedB,nz))
         allocate(cff_tangential(segment%HI%IsdB:segment%HI%IedB,segment%HI%JsdB:segment%HI%JedB,nz))
         do k=1,nz
           if (gamma_u < 1.0) then
             rx_tang_obl(segment%HI%IsdB,J,k) = segment%rx_norm_obl(segment%HI%isd,J,k)
             rx_tang_obl(segment%HI%IedB,J,k) = segment%rx_norm_obl(segment%HI%ied,J,k)
             ry_tang_obl(segment%HI%IsdB,J,k) = segment%ry_norm_obl(segment%HI%isd,J,k)
             ry_tang_obl(segment%HI%IedB,J,k) = segment%ry_norm_obl(segment%HI%ied,J,k)
             cff_tangential(segment%HI%IsdB,J,k) = segment%cff_normal(segment%HI%isd,J,k)
             cff_tangential(segment%HI%IedB,J,k) = segment%cff_normal(segment%HI%ied,J,k)
             do I=segment%HI%IsdB+1,segment%HI%IedB-1
               rx_tang_obl(I,J,k) = 0.5*(segment%rx_norm_obl(i,J,k) + segment%rx_norm_obl(i+1,J,k))
               ry_tang_obl(I,J,k) = 0.5*(segment%ry_norm_obl(i,J,k) + segment%ry_norm_obl(i+1,J,k))
               cff_tangential(I,J,k) = 0.5*(segment%cff_normal(i,J,k) + segment%cff_normal(i+1,J,k))
             enddo
           else
             do I=segment%HI%IsdB,segment%HI%IedB
               dhdt = u_old(I,j+1,k)-u_new(I,j+1,k)   !old-new
               dhdy = u_new(I,j+1,k)-u_new(I,j+2,k) !in new time backward sasha for I-1
               if (dhdt*(segment%grad_tan(i,1,k) + segment%grad_tan(i+1,1,k)) > 0.0) then
                 dhdx = segment%grad_tan(i,1,k)
               elseif (dhdt*(segment%grad_tan(i,1,k) + segment%grad_tan(i+1,1,k)) == 0.0) then
                 dhdx = 0.0
               else
                 dhdx = segment%grad_tan(i+1,1,k)
               endif
               if (dhdt*dhdy < 0.0) dhdt = 0.0
               ry_new = dhdt*dhdy
               cff_new = max(dhdx*dhdx + dhdy*dhdy, eps)
               rx_new = min(cff_new,max(dhdt*dhdx,-cff_new))
               rx_tang_obl(I,j,k) = rx_new
               ry_tang_obl(i,J,k) = ry_new
               cff_tangential(i,J,k) = cff_new
             enddo
           endif
         enddo
         if (segment%oblique_tan) then
           do k=1,nz ;  do I=segment%HI%IsdB,segment%HI%IedB
             rx_avg = rx_tang_obl(I,J,k)
             ry_avg = ry_tang_obl(I,J,k)
             cff_avg = cff_tangential(I,J,k)
             segment%tangential_vel(I,J,k) = ((cff_avg*u_new(I,j+1,k) + ry_avg*u_new(I,j+2,k)) - &
                                              (max(rx_avg,0.0)*segment%grad_tan(i,2,k) + &
                                               min(rx_avg,0.0)*segment%grad_tan(i+1,2,k)) ) / &
                                             (cff_avg + ry_avg)
           enddo ; enddo
         endif
         if (segment%nudged_tan) then
           do k=1,nz ; do I=segment%HI%IsdB,segment%HI%IedB
             ! dhdt gets set to 0 on inflow in oblique case
             if (ry_tang_obl(I,J,k) <= 0.0) then
               tau = segment%Velocity_nudging_timescale_in
             else
               tau = segment%Velocity_nudging_timescale_out
             endif
             gamma_2 = dt / (tau + dt)
             segment%tangential_vel(I,J,k) = (1.0 - gamma_2) * segment%tangential_vel(I,J,k) + &
                                             gamma_2 * segment%nudged_tangential_vel(I,J,k)
           enddo ; enddo
         endif
         if (segment%oblique_grad) then
           Is_obc = max(segment%HI%IsdB,G%isd+1)
           Ie_obc = min(segment%HI%IedB,G%ied-1)
           do k=1,nz ;  do I=segment%HI%IsdB+1,segment%HI%IedB-1
             rx_avg = rx_tang_obl(I,J,k)
             ry_avg = ry_tang_obl(I,J,k)
             cff_avg = cff_tangential(I,J,k)
             segment%tangential_grad(I,J,k) = &
                 ((cff_avg*(u_new(I,j+2,k) - u_new(I,j+1,k))*G%IdyBu(I,J+1) + &
                    ry_avg*(u_new(I,j+3,k) - u_new(I,j+2,k))*G%IdyBu(I,J+2)) - &
                                  (max(rx_avg,0.0)*segment%grad_gradient(i,2,k) + &
                                   min(rx_avg,0.0)*segment%grad_gradient(i+1,2,k))) / &
                                 (cff_avg + ry_avg)
           enddo ; enddo
         endif
         if (segment%nudged_grad) then
           do k=1,nz ; do J=segment%HI%JsdB,segment%HI%JedB
             ! dhdt gets set to 0 on inflow in oblique case
             if (ry_tang_obl(I,J,k) <= 0.0) then
               tau = segment%Velocity_nudging_timescale_in
             else
               tau = segment%Velocity_nudging_timescale_out
             endif
             gamma_2 = dt / (tau + dt)
             segment%tangential_grad(I,J,k) = (1.0 - gamma_2) * segment%tangential_grad(I,J,k) + &
                                 gamma_2 * segment%nudged_tangential_grad(I,J,k)
           enddo ; enddo
         endif
         deallocate(rx_tang_obl)
         deallocate(ry_tang_obl)
         deallocate(cff_tangential)
       endif
     endif
  enddo

  ! Actually update u_new, v_new
  call open_boundary_apply_normal_flow(OBC, G, u_new, v_new)

  call pass_vector(u_new, v_new, G%Domain, clock=id_clock_pass)

end subroutine radiation_open_bdry_conds

!> Applies OBC values stored in segments to 3d u,v fields
subroutine open_boundary_apply_normal_flow(OBC, G, u, v)
  ! Arguments
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary control structure
  type(ocean_grid_type),                     intent(inout) :: G   !< Ocean grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: u   !< u field to update on open
                                                                  !! boundaries [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: v   !< v field to update on open
                                                                  !! boundaries [L T-1 ~> m s-1]
  ! Local variables
  integer :: i, j, k, n
  type(OBC_segment_type), pointer :: segment => NULL()

  if (.not.associated(OBC)) return ! Bail out if OBC is not available

  do n=1,OBC%number_of_segments
    segment => OBC%segment(n)
    if (.not. segment%on_pe) then
      cycle
    elseif (segment%radiation .or. segment%oblique .or. segment%gradient) then
      if (segment%is_E_or_W) then
        I=segment%HI%IsdB
        do k=1,G%ke ;  do j=segment%HI%jsd,segment%HI%jed
          u(I,j,k) = segment%normal_vel(I,j,k)
        enddo ; enddo
      elseif (segment%is_N_or_S) then
        J=segment%HI%JsdB
        do k=1,G%ke ;  do i=segment%HI%isd,segment%HI%ied
          v(i,J,k) = segment%normal_vel(i,J,k)
        enddo ; enddo
      endif
    endif
  enddo

end subroutine open_boundary_apply_normal_flow

!> Applies zero values to 3d u,v fields on OBC segments
subroutine open_boundary_zero_normal_flow(OBC, G, u, v)
  ! Arguments
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary control structure
  type(ocean_grid_type),                     intent(inout) :: G   !< Ocean grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: u   !< u field to update on open boundaries
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: v   !< v field to update on open boundaries
  ! Local variables
  integer :: i, j, k, n
  type(OBC_segment_type), pointer :: segment => NULL()

  if (.not.associated(OBC)) return ! Bail out if OBC is not available

  do n=1,OBC%number_of_segments
    segment => OBC%segment(n)
    if (.not. segment%on_pe) then
      cycle
    elseif (segment%is_E_or_W) then
      I=segment%HI%IsdB
      do k=1,G%ke ;  do j=segment%HI%jsd,segment%HI%jed
        u(I,j,k) = 0.
      enddo ; enddo
    elseif (segment%is_N_or_S) then
      J=segment%HI%JsdB
      do k=1,G%ke ;  do i=segment%HI%isd,segment%HI%ied
        v(i,J,k) = 0.
      enddo ; enddo
    endif
  enddo

end subroutine open_boundary_zero_normal_flow

!> Calculate the tangential gradient of the normal flow at the boundary q-points.
subroutine gradient_at_q_points(G, segment, uvel, vvel)
  type(ocean_grid_type), intent(in) :: G !< Ocean grid structure
  type(OBC_segment_type), pointer :: segment !< OBC segment structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: uvel !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: vvel !< meridional velocity [L T-1 ~> m s-1]
  integer :: i,j,k

  if (.not. segment%on_pe) return

  if (segment%is_E_or_W) then
    if (segment%direction == OBC_DIRECTION_E) then
      I=segment%HI%isdB
      do k=1,G%ke
        do J=max(segment%HI%JsdB, G%HI%JsdB+1),min(segment%HI%JedB, G%HI%JedB-1)
          segment%grad_normal(J,1,k) = (uvel(I-1,j+1,k)-uvel(I-1,j,k)) * G%mask2dBu(I-1,J)
          segment%grad_normal(J,2,k) = (uvel(I,j+1,k)-uvel(I,j,k)) * G%mask2dBu(I,J)
        enddo
      enddo
      if (segment%oblique_tan) then
        do k=1,G%ke
          do J=max(segment%HI%jsd-1, G%HI%jsd),min(segment%HI%jed+1, G%HI%jed)
            segment%grad_tan(j,1,k) = (vvel(i-1,J,k)-vvel(i-1,J-1,k)) * G%mask2dT(i-1,j)
            segment%grad_tan(j,2,k) = (vvel(i,J,k)-vvel(i,J-1,k)) * G%mask2dT(i,j)
          enddo
        enddo
      endif
      if (segment%oblique_grad) then
        do k=1,G%ke
          do J=max(segment%HI%jsd, G%HI%jsd+1),min(segment%HI%jed, G%HI%jed-1)
            segment%grad_gradient(j,1,k) = (((vvel(i-1,J,k) - vvel(i-2,J,k))*G%IdxBu(I-2,J)) - &
                 (vvel(i-1,J-1,k) - vvel(i-2,J-1,k))*G%IdxBu(I-2,J-1)) * G%mask2dCu(I-2,j)
            segment%grad_gradient(j,2,k) = (((vvel(i,J,k) - vvel(i-1,J,k))*G%IdxBu(I-1,J)) - &
                 (vvel(i,J-1,k) - vvel(i-1,J-1,k))*G%IdxBu(I-1,J-1)) * G%mask2dCu(I-1,j)
          enddo
        enddo
      endif
    else ! western segment
      I=segment%HI%isdB
      do k=1,G%ke
        do J=max(segment%HI%JsdB, G%HI%JsdB+1),min(segment%HI%JedB, G%HI%JedB-1)
          segment%grad_normal(J,1,k) = (uvel(I+1,j+1,k)-uvel(I+1,j,k)) * G%mask2dBu(I+1,J)
          segment%grad_normal(J,2,k) = (uvel(I,j+1,k)-uvel(I,j,k)) * G%mask2dBu(I,J)
        enddo
      enddo
      if (segment%oblique_tan) then
        do k=1,G%ke
          do J=max(segment%HI%jsd-1, G%HI%jsd),min(segment%HI%jed+1, G%HI%jed)
            segment%grad_tan(j,1,k) = (vvel(i+2,J,k)-vvel(i+2,J-1,k)) * G%mask2dT(i+2,j)
            segment%grad_tan(j,2,k) = (vvel(i+1,J,k)-vvel(i+1,J-1,k)) * G%mask2dT(i+1,j)
          enddo
        enddo
      endif
      if (segment%oblique_grad) then
        do k=1,G%ke
          do J=max(segment%HI%jsd, G%HI%jsd+1),min(segment%HI%jed, G%HI%jed-1)
            segment%grad_gradient(j,1,k) = (((vvel(i+3,J,k) - vvel(i+2,J,k))*G%IdxBu(I+2,J)) - &
                 (vvel(i+3,J-1,k) - vvel(i+2,J-1,k))*G%IdxBu(I+2,J-1)) * G%mask2dCu(I+2,j)
            segment%grad_gradient(j,2,k) = (((vvel(i+2,J,k) - vvel(i+1,J,k))*G%IdxBu(I+1,J)) - &
                 (vvel(i+2,J-1,k) - vvel(i+1,J-1,k))*G%IdxBu(I+1,J-1)) * G%mask2dCu(I+1,j)
          enddo
        enddo
      endif
    endif
  elseif (segment%is_N_or_S) then
    if (segment%direction == OBC_DIRECTION_N) then
      J=segment%HI%jsdB
      do k=1,G%ke
        do I=max(segment%HI%IsdB, G%HI%IsdB+1),min(segment%HI%IedB, G%HI%IedB-1)
          segment%grad_normal(I,1,k) = (vvel(i+1,J-1,k)-vvel(i,J-1,k)) * G%mask2dBu(I,J-1)
          segment%grad_normal(I,2,k) = (vvel(i+1,J,k)-vvel(i,J,k)) * G%mask2dBu(I,J)
        enddo
      enddo
      if (segment%oblique_tan) then
        do k=1,G%ke
          do I=max(segment%HI%isd-1, G%HI%isd),min(segment%HI%ied+1, G%HI%ied)
            segment%grad_tan(i,1,k) = (uvel(I,j-1,k)-uvel(I-1,j-1,k)) * G%mask2dT(i,j-1)
            segment%grad_tan(i,2,k) = (uvel(I,j,k)-uvel(I-1,j,k)) * G%mask2dT(i,j)
          enddo
        enddo
      endif
      if (segment%oblique_grad) then
        do k=1,G%ke
          do I=max(segment%HI%isd, G%HI%isd+1),min(segment%HI%ied, G%HI%ied-1)
            segment%grad_gradient(i,1,k) = (((uvel(I,j-1,k) - uvel(I,j-2,k))*G%IdyBu(I,J-2)) - &
                 (uvel(I-1,j-1,k) - uvel(I-1,j-2,k))*G%IdyBu(I-1,J-2)) * G%mask2dCv(i,J-2)
            segment%grad_gradient(i,2,k) = (((uvel(I,j,k) - uvel(I,j-1,k))*G%IdyBu(I,J-1)) - &
                 (uvel(I-1,j,k) - uvel(I-1,j-1,k))*G%IdyBu(I-1,J-1)) * G%mask2dCv(i,J-1)
          enddo
        enddo
      endif
    else ! south segment
      J=segment%HI%jsdB
      do k=1,G%ke
        do I=max(segment%HI%IsdB, G%HI%IsdB+1),min(segment%HI%IedB, G%HI%IedB-1)
          segment%grad_normal(I,1,k) = (vvel(i+1,J+1,k)-vvel(i,J+1,k)) * G%mask2dBu(I,J+1)
          segment%grad_normal(I,2,k) = (vvel(i+1,J,k)-vvel(i,J,k)) * G%mask2dBu(I,J)
        enddo
      enddo
      if (segment%oblique_tan) then
        do k=1,G%ke
          do I=max(segment%HI%isd-1, G%HI%isd),min(segment%HI%ied+1, G%HI%ied)
            segment%grad_tan(i,1,k) = (uvel(I,j+2,k)-uvel(I-1,j+2,k)) * G%mask2dT(i,j+2)
            segment%grad_tan(i,2,k) = (uvel(I,j+1,k)-uvel(I-1,j+1,k)) * G%mask2dT(i,j+1)
          enddo
        enddo
      endif
      if (segment%oblique_grad) then
        do k=1,G%ke
          do I=max(segment%HI%isd, G%HI%isd+1),min(segment%HI%ied, G%HI%ied-1)
            segment%grad_gradient(i,1,k) = (((uvel(I,j+3,k) - uvel(I,j+2,k))*G%IdyBu(I,J+2)) - &
                 (uvel(I-1,j+3,k) - uvel(I-1,j+2,k))*G%IdyBu(I-1,J+2)) * G%mask2dCv(i,J+2)
            segment%grad_gradient(i,2,k) = (((uvel(I,j+2,k) - uvel(I,j+1,k))*G%IdyBu(I,J+1)) - &
                 (uvel(I-1,j+2,k) - uvel(I-1,j+1,k))*G%IdyBu(I-1,J+1)) * G%mask2dCv(i,J+1)
          enddo
        enddo
      endif
    endif
  endif

end subroutine gradient_at_q_points


!> Sets the initial values of the tracer open boundary conditions.
!! Redoing this elsewhere.
subroutine set_tracer_data(OBC, tv, h, G, PF, tracer_Reg)
  type(ocean_grid_type),                     intent(inout) :: G !< Ocean grid structure
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary structure
  type(thermo_var_ptrs),                     intent(inout) :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(inout) :: h !< Thickness
  type(param_file_type),                     intent(in)    :: PF !< Parameter file handle
  type(tracer_registry_type),                pointer       :: tracer_Reg !< Tracer registry
  ! Local variables
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, nz, n
  integer :: isd_off, jsd_off
  integer :: IsdB, IedB, JsdB, JedB
  type(OBC_segment_type), pointer :: segment  => NULL() ! pointer to segment type list
  character(len=40)  :: mdl = "set_tracer_data" ! This subroutine's name.
  character(len=200) :: filename, OBC_file, inputdir ! Strings for file/path

  real :: temp_u(G%domain%niglobal+1,G%domain%njglobal)
  real :: temp_v(G%domain%niglobal,G%domain%njglobal+1)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ! For now, there are no radiation conditions applied to the thicknesses, since
  ! the thicknesses might not be physically motivated.  Instead, sponges should be
  ! used to enforce the near-boundary layer structure.

  if (associated(tv%T)) then

    call pass_var(tv%T, G%Domain)
    call pass_var(tv%S, G%Domain)

    do n=1,OBC%number_of_segments
      segment => OBC%segment(n)
      if (.not. segment%on_pe) cycle

      if (segment%direction == OBC_DIRECTION_E) then
        I=segment%HI%IsdB
        do k=1,G%ke ;  do j=segment%HI%jsd,segment%HI%jed
          tv%T(i+1,j,k) = tv%T(i,j,k) ; tv%S(i+1,j,k) = tv%S(i,j,k)
        enddo ; enddo
      elseif (segment%direction == OBC_DIRECTION_W) then
        I=segment%HI%IsdB
        do k=1,G%ke ;  do j=segment%HI%jsd,segment%HI%jed
          tv%T(i,j,k) = tv%T(i+1,j,k) ; tv%S(i,j,k) = tv%S(i+1,j,k)
        enddo ; enddo
      elseif (segment%direction == OBC_DIRECTION_N) then
        J=segment%HI%JsdB
        do k=1,G%ke ;  do i=segment%HI%isd,segment%HI%ied
          tv%T(i,j+1,k) = tv%T(i,j,k) ; tv%S(i,j+1,k) = tv%S(i,j,k)
        enddo ; enddo
      elseif (segment%direction == OBC_DIRECTION_S) then
        J=segment%HI%JsdB
        do k=1,G%ke ;  do i=segment%HI%isd,segment%HI%ied
          tv%T(i,j,k) = tv%T(i,j+1,k) ; tv%S(i,j,k) = tv%S(i,j+1,k)
        enddo ; enddo
      endif
    enddo
  endif

end subroutine set_tracer_data

!> Needs documentation
function lookup_seg_field(OBC_seg,field)
  type(OBC_segment_type), pointer :: OBC_seg !< OBC segment
  character(len=32), intent(in) :: field !< The field name
  integer :: lookup_seg_field
  ! Local variables
  integer :: n

  lookup_seg_field=-1
  do n=1,OBC_seg%num_fields
   if (trim(field) == OBC_seg%field(n)%name) then
     lookup_seg_field=n
     return
   endif
  enddo

end function lookup_seg_field


!> Allocate segment data fields
subroutine allocate_OBC_segment_data(OBC, segment)
  type(ocean_OBC_type),   pointer       :: OBC     !< Open boundary structure
  type(OBC_segment_type), intent(inout) :: segment !< Open boundary segment
  ! Local variables
  integer :: isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB
  integer :: IscB, IecB, JscB, JecB
  character(len=40)  :: mdl = "allocate_OBC_segment_data" ! This subroutine's name.

  isd = segment%HI%isd ; ied = segment%HI%ied
  jsd = segment%HI%jsd ; jed = segment%HI%jed
  IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
  JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB
  IscB = segment%HI%IscB ; IecB = segment%HI%IecB
  JscB = segment%HI%JscB ; JecB = segment%HI%JecB

  if (.not. segment%on_pe) return

  if (segment%is_E_or_W) then
    ! If these are just Flather, change update_OBC_segment_data accordingly
    allocate(segment%Cg(IsdB:IedB,jsd:jed));                  segment%Cg(:,:)=0.
    allocate(segment%Htot(IsdB:IedB,jsd:jed));                segment%Htot(:,:)=0.0
    allocate(segment%h(IsdB:IedB,jsd:jed,OBC%ke));            segment%h(:,:,:)=0.0
    allocate(segment%eta(IsdB:IedB,jsd:jed));                 segment%eta(:,:)=0.0
    if (segment%radiation) then
      allocate(segment%rx_norm_rad(IsdB:IedB,jsd:jed,OBC%ke));  segment%rx_norm_rad(:,:,:)=0.0
    endif
    allocate(segment%normal_vel(IsdB:IedB,jsd:jed,OBC%ke));   segment%normal_vel(:,:,:)=0.0
    allocate(segment%normal_vel_bt(IsdB:IedB,jsd:jed));       segment%normal_vel_bt(:,:)=0.0
    allocate(segment%normal_trans(IsdB:IedB,jsd:jed,OBC%ke)); segment%normal_trans(:,:,:)=0.0
    if (segment%nudged) then
      allocate(segment%nudged_normal_vel(IsdB:IedB,jsd:jed,OBC%ke)); segment%nudged_normal_vel(:,:,:)=0.0
    endif
    if (segment%radiation_tan .or. segment%nudged_tan .or. segment%specified_tan .or. &
        segment%oblique_tan .or. OBC%computed_vorticity .or. OBC%computed_strain) then
      allocate(segment%tangential_vel(IsdB:IedB,JsdB:JedB,OBC%ke)); segment%tangential_vel(:,:,:)=0.0
    endif
    if (segment%nudged_tan) then
      allocate(segment%nudged_tangential_vel(IsdB:IedB,JsdB:JedB,OBC%ke)); segment%nudged_tangential_vel(:,:,:)=0.0
    endif
    if (segment%nudged_grad) then
      allocate(segment%nudged_tangential_grad(IsdB:IedB,JsdB:JedB,OBC%ke)); segment%nudged_tangential_grad(:,:,:)=0.0
    endif
    if (OBC%specified_vorticity .or. OBC%specified_strain .or. segment%radiation_grad .or. &
              segment%oblique_grad .or. segment%specified_grad) then
      allocate(segment%tangential_grad(IsdB:IedB,JsdB:JedB,OBC%ke)); segment%tangential_grad(:,:,:)=0.0
    endif
    if (segment%oblique) then
      allocate(segment%grad_normal(JsdB:JedB,2,OBC%ke));      segment%grad_normal(:,:,:) = 0.0
      allocate(segment%rx_norm_obl(IsdB:IedB,jsd:jed,OBC%ke));  segment%rx_norm_obl(:,:,:)=0.0
      allocate(segment%ry_norm_obl(IsdB:IedB,jsd:jed,OBC%ke));  segment%ry_norm_obl(:,:,:)=0.0
      allocate(segment%cff_normal(IsdB:IedB,jsd:jed,OBC%ke)); segment%cff_normal(:,:,:)=0.0
    endif
    if (segment%oblique_tan) then
      allocate(segment%grad_tan(jsd-1:jed+1,2,OBC%ke));           segment%grad_tan(:,:,:) = 0.0
    endif
    if (segment%oblique_grad) then
      allocate(segment%grad_gradient(jsd:jed,2,OBC%ke));      segment%grad_gradient(:,:,:) = 0.0
    endif
  endif

  if (segment%is_N_or_S) then
    ! If these are just Flather, change update_OBC_segment_data accordingly
    allocate(segment%Cg(isd:ied,JsdB:JedB));                  segment%Cg(:,:)=0.
    allocate(segment%Htot(isd:ied,JsdB:JedB));                segment%Htot(:,:)=0.0
    allocate(segment%h(isd:ied,JsdB:JedB,OBC%ke));            segment%h(:,:,:)=0.0
    allocate(segment%eta(isd:ied,JsdB:JedB));                 segment%eta(:,:)=0.0
    if (segment%radiation) then
      allocate(segment%ry_norm_rad(isd:ied,JsdB:JedB,OBC%ke));  segment%ry_norm_rad(:,:,:)=0.0
    endif
    allocate(segment%normal_vel(isd:ied,JsdB:JedB,OBC%ke));   segment%normal_vel(:,:,:)=0.0
    allocate(segment%normal_vel_bt(isd:ied,JsdB:JedB));       segment%normal_vel_bt(:,:)=0.0
    allocate(segment%normal_trans(isd:ied,JsdB:JedB,OBC%ke)); segment%normal_trans(:,:,:)=0.0
    if (segment%nudged) then
      allocate(segment%nudged_normal_vel(isd:ied,JsdB:JedB,OBC%ke)); segment%nudged_normal_vel(:,:,:)=0.0
    endif
    if (segment%radiation_tan .or. segment%nudged_tan .or. segment%specified_tan .or. &
        segment%oblique_tan .or. OBC%computed_vorticity .or. OBC%computed_strain) then
      allocate(segment%tangential_vel(IsdB:IedB,JsdB:JedB,OBC%ke)); segment%tangential_vel(:,:,:)=0.0
    endif
    if (segment%nudged_tan) then
      allocate(segment%nudged_tangential_vel(IsdB:IedB,JsdB:JedB,OBC%ke)); segment%nudged_tangential_vel(:,:,:)=0.0
    endif
    if (segment%nudged_grad) then
      allocate(segment%nudged_tangential_grad(IsdB:IedB,JsdB:JedB,OBC%ke)); segment%nudged_tangential_grad(:,:,:)=0.0
    endif
    if (OBC%specified_vorticity .or. OBC%specified_strain .or. segment%radiation_grad .or. &
              segment%oblique_grad .or. segment%specified_grad) then
      allocate(segment%tangential_grad(IsdB:IedB,JsdB:JedB,OBC%ke)); segment%tangential_grad(:,:,:)=0.0
    endif
    if (segment%oblique) then
      allocate(segment%grad_normal(IsdB:IedB,2,OBC%ke));      segment%grad_normal(:,:,:) = 0.0
      allocate(segment%rx_norm_obl(isd:ied,JsdB:JedB,OBC%ke));  segment%rx_norm_obl(:,:,:)=0.0
      allocate(segment%ry_norm_obl(isd:ied,JsdB:JedB,OBC%ke));  segment%ry_norm_obl(:,:,:)=0.0
      allocate(segment%cff_normal(isd:ied,JsdB:JedB,OBC%ke)); segment%cff_normal(:,:,:)=0.0
    endif
    if (segment%oblique_tan) then
      allocate(segment%grad_tan(isd-1:ied+1,2,OBC%ke));           segment%grad_tan(:,:,:) = 0.0
    endif
    if (segment%oblique_grad) then
      allocate(segment%grad_gradient(isd:ied,2,OBC%ke));      segment%grad_gradient(:,:,:) = 0.0
    endif
  endif

end subroutine allocate_OBC_segment_data

!> Deallocate segment data fields
subroutine deallocate_OBC_segment_data(OBC, segment)
  type(ocean_OBC_type),   pointer       :: OBC     !< Open boundary structure
  type(OBC_segment_type), intent(inout) :: segment !< Open boundary segment
  ! Local variables
  character(len=40)  :: mdl = "deallocate_OBC_segment_data" ! This subroutine's name.

  if (.not. segment%on_pe) return

  if (associated (segment%Cg)) deallocate(segment%Cg)
  if (associated (segment%Htot)) deallocate(segment%Htot)
  if (associated (segment%h)) deallocate(segment%h)
  if (associated (segment%eta)) deallocate(segment%eta)
  if (associated (segment%rx_norm_rad)) deallocate(segment%rx_norm_rad)
  if (associated (segment%ry_norm_rad)) deallocate(segment%ry_norm_rad)
  if (associated (segment%rx_norm_obl)) deallocate(segment%rx_norm_obl)
  if (associated (segment%ry_norm_obl)) deallocate(segment%ry_norm_obl)
  if (associated (segment%cff_normal)) deallocate(segment%cff_normal)
  if (associated (segment%grad_normal)) deallocate(segment%grad_normal)
  if (associated (segment%grad_tan)) deallocate(segment%grad_tan)
  if (associated (segment%grad_gradient)) deallocate(segment%grad_gradient)
  if (associated (segment%normal_vel)) deallocate(segment%normal_vel)
  if (associated (segment%normal_vel_bt)) deallocate(segment%normal_vel_bt)
  if (associated (segment%normal_trans)) deallocate(segment%normal_trans)
  if (associated (segment%nudged_normal_vel)) deallocate(segment%nudged_normal_vel)
  if (associated (segment%tangential_vel)) deallocate(segment%tangential_vel)
  if (associated (segment%nudged_tangential_vel)) deallocate(segment%nudged_tangential_vel)
  if (associated (segment%nudged_tangential_grad)) deallocate(segment%nudged_tangential_grad)
  if (associated (segment%tangential_grad)) deallocate(segment%tangential_grad)
  if (associated (segment%tr_Reg)) call segment_tracer_registry_end(segment%tr_Reg)


end subroutine deallocate_OBC_segment_data

!> Set tangential velocities outside of open boundaries to silly values
!! (used for checking the interior state is independent of values outside
!! of the domain).
subroutine open_boundary_test_extern_uv(G, OBC, u, v)
  type(ocean_grid_type),                     intent(in)    :: G !< Ocean grid structure
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary structure
  real, dimension(SZIB_(G),SZJ_(G), SZK_(G)),intent(inout) :: u !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G), SZK_(G)),intent(inout) :: v !< Meridional velocity [L T-1 ~> m s-1]
  ! Local variables
  integer :: i, j, k, n

  if (.not. associated(OBC)) return

  do n = 1, OBC%number_of_segments
    do k = 1, G%ke
      if (OBC%segment(n)%is_N_or_S) then
        J = OBC%segment(n)%HI%JsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
          do I = OBC%segment(n)%HI%IsdB, OBC%segment(n)%HI%IedB
            u(I,j+1,k) = OBC%silly_u
          enddo
        else
          do I = OBC%segment(n)%HI%IsdB, OBC%segment(n)%HI%IedB
            u(I,j,k) = OBC%silly_u
          enddo
        endif
      elseif (OBC%segment(n)%is_E_or_W) then
        I = OBC%segment(n)%HI%IsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
          do J = OBC%segment(n)%HI%JsdB, OBC%segment(n)%HI%JedB
            v(i+1,J,k) = OBC%silly_u
          enddo
        else
          do J = OBC%segment(n)%HI%JsdB, OBC%segment(n)%HI%JedB
            v(i,J,k) = OBC%silly_u
          enddo
        endif
      endif
    enddo
  enddo

end subroutine open_boundary_test_extern_uv

!> Set thicknesses outside of open boundaries to silly values
!! (used for checking the interior state is independent of values outside
!! of the domain).
subroutine open_boundary_test_extern_h(G, GV, OBC, h)
  type(ocean_grid_type),                     intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV  !<  Ocean vertical grid structure
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)),intent(inout) :: h   !< Layer thickness [H ~> m or kg m-2]
  ! Local variables
  real :: silly_h  ! A silly thickness for testing [H ~> m or kg m-2]
  integer :: i, j, k, n

  if (.not. associated(OBC)) return

  silly_h = GV%Z_to_H*OBC%silly_h

  do n = 1, OBC%number_of_segments
    do k = 1, GV%ke
      if (OBC%segment(n)%is_N_or_S) then
        J = OBC%segment(n)%HI%JsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
          do i = OBC%segment(n)%HI%isd, OBC%segment(n)%HI%ied
            h(i,j+1,k) = silly_h
          enddo
        else
          do i = OBC%segment(n)%HI%isd, OBC%segment(n)%HI%ied
            h(i,j,k) = silly_h
          enddo
        endif
      elseif (OBC%segment(n)%is_E_or_W) then
        I = OBC%segment(n)%HI%IsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
          do j = OBC%segment(n)%HI%jsd, OBC%segment(n)%HI%jed
            h(i+1,j,k) = silly_h
          enddo
        else
          do j = OBC%segment(n)%HI%jsd, OBC%segment(n)%HI%jed
            h(i,j,k) = silly_h
          enddo
        endif
      endif
    enddo
  enddo

end subroutine open_boundary_test_extern_h

!> Update the OBC values on the segments.
subroutine update_OBC_segment_data(G, GV, US, OBC, tv, h, Time)
  type(ocean_grid_type),                     intent(in)    :: G    !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV   !<  Ocean vertical grid structure
  type(unit_scale_type),                     intent(in)    :: US   !< A dimensional unit scaling type
  type(ocean_OBC_type),                      pointer       :: OBC  !< Open boundary structure
  type(thermo_var_ptrs),                     intent(in)    :: tv   !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(inout) :: h    !< Thickness [m]
  type(time_type),                           intent(in)    :: Time !< Model time
  ! Local variables
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB, n, m, nz
  character(len=40)  :: mdl = "set_OBC_segment_data" ! This subroutine's name.
  character(len=200) :: filename, OBC_file, inputdir ! Strings for file/path
  type(OBC_segment_type), pointer :: segment => NULL()
  integer, dimension(4) :: siz,siz2
  real :: sumh ! column sum of thicknesses [m]
  integer :: ni_seg, nj_seg  ! number of src gridpoints along the segments
  integer :: i2, j2          ! indices for referencing local domain array
  integer :: is_obc, ie_obc, js_obc, je_obc  ! segment indices within local domain
  integer :: ishift, jshift  ! offsets for staggered locations
  real, dimension(:,:), pointer :: seg_vel => NULL()  ! pointer to segment velocity array
  real, dimension(:,:), pointer :: seg_trans => NULL()  ! pointer to segment transport array
  real, dimension(:,:,:), allocatable :: tmp_buffer
  real, dimension(:), allocatable :: h_stack
  integer :: is_obc2, js_obc2
  real :: net_H_src, net_H_int, scl_fac
  real, pointer, dimension(:,:)   :: normal_trans_bt=>NULL() ! barotropic transport

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  nz=G%ke

  if (.not. associated(OBC)) return

  do n = 1, OBC%number_of_segments
    segment => OBC%segment(n)

    if (.not. segment%on_pe) cycle ! continue to next segment if not in computational domain

    ni_seg = segment%ie_obc-segment%is_obc+1
    nj_seg = segment%je_obc-segment%js_obc+1
    is_obc = max(segment%is_obc,isd-1)
    ie_obc = min(segment%ie_obc,ied)
    js_obc = max(segment%js_obc,jsd-1)
    je_obc = min(segment%je_obc,jed)

! Calculate auxiliary fields at staggered locations.
! Segment indices are on q points:
!
!       |-----------|------------|-----------|-----------|  J_obc
!     Is_obc                                          Ie_obc
!
! i2 has to start at Is_obc+1 and end at Ie_obc.
! j2 is J_obc and jshift has to be +1 at both the north and south.

     ! calculate auxiliary fields at staggered locations
    ishift=0;jshift=0
    if (segment%is_E_or_W) then
      allocate(normal_trans_bt(segment%HI%IsdB:segment%HI%IedB,segment%HI%jsd:segment%HI%jed))
      normal_trans_bt(:,:)=0.0
      if (segment%direction == OBC_DIRECTION_W) ishift=1
      I=segment%HI%IsdB
      do j=segment%HI%jsd,segment%HI%jed
        segment%Cg(I,j) = sqrt(GV%g_prime(1)*G%bathyT(i+ishift,j))
        segment%Htot(I,j)=0.0
        do k=1,G%ke
          segment%h(I,j,k) = h(i+ishift,j,k)
          segment%Htot(I,j)=segment%Htot(I,j)+segment%h(I,j,k)
        enddo
      enddo
    else! (segment%direction == OBC_DIRECTION_N .or. segment%direction == OBC_DIRECTION_S)
      allocate(normal_trans_bt(segment%HI%isd:segment%HI%ied,segment%HI%JsdB:segment%HI%JedB))
      normal_trans_bt(:,:)=0.0
      if (segment%direction == OBC_DIRECTION_S) jshift=1
      J=segment%HI%JsdB
      do i=segment%HI%isd,segment%HI%ied
        segment%Cg(i,J) = sqrt(GV%g_prime(1)*G%bathyT(i,j+jshift))
        segment%Htot(i,J)=0.0
        do k=1,G%ke
          segment%h(i,J,k) = h(i,j+jshift,k)
          segment%Htot(i,J)=segment%Htot(i,J)+segment%h(i,J,k)
        enddo
      enddo
    endif

    allocate(h_stack(G%ke))
    h_stack(:) = 0.0
    do m = 1,segment%num_fields
      if (segment%field(m)%fid > 0) then
        siz(1)=size(segment%field(m)%buffer_src,1)
        siz(2)=size(segment%field(m)%buffer_src,2)
        siz(3)=size(segment%field(m)%buffer_src,3)
        if (.not.associated(segment%field(m)%buffer_dst)) then
          if (siz(3) /= segment%field(m)%nk_src) call MOM_error(FATAL,'nk_src inconsistency')
          if (segment%field(m)%nk_src > 1) then
            if (segment%is_E_or_W) then
              if (segment%field(m)%name == 'V' .or. segment%field(m)%name == 'DVDX') then
                allocate(segment%field(m)%buffer_dst(is_obc:ie_obc,js_obc:je_obc,G%ke))
              else
                allocate(segment%field(m)%buffer_dst(is_obc:ie_obc,js_obc+1:je_obc,G%ke))
              endif
              if (segment%field(m)%name == 'U') then
                allocate(segment%field(m)%bt_vel(is_obc:ie_obc,js_obc+1:je_obc))
                segment%field(m)%bt_vel(:,:)=0.0
              endif
            else
              if (segment%field(m)%name == 'U' .or. segment%field(m)%name == 'DUDY') then
                allocate(segment%field(m)%buffer_dst(is_obc:ie_obc,js_obc:je_obc,G%ke))
              else
                allocate(segment%field(m)%buffer_dst(is_obc+1:ie_obc,js_obc:je_obc,G%ke))
              endif
              if (segment%field(m)%name == 'V') then
                allocate(segment%field(m)%bt_vel(is_obc+1:ie_obc,js_obc:je_obc))
                segment%field(m)%bt_vel(:,:)=0.0
              endif
            endif
          else
            if (segment%is_E_or_W) then
              if (segment%field(m)%name == 'V' .or. segment%field(m)%name == 'DVDX') then
                allocate(segment%field(m)%buffer_dst(is_obc:ie_obc,js_obc:je_obc,1))
              else
                allocate(segment%field(m)%buffer_dst(is_obc:ie_obc,js_obc+1:je_obc,1))
              endif
              if (segment%field(m)%name == 'U') then
                allocate(segment%field(m)%bt_vel(is_obc:ie_obc,js_obc+1:je_obc))
                segment%field(m)%bt_vel(:,:)=0.0
              endif
            else
              if (segment%field(m)%name == 'U' .or. segment%field(m)%name == 'DUDY') then
                allocate(segment%field(m)%buffer_dst(is_obc:ie_obc,js_obc:je_obc,1))
              else
                allocate(segment%field(m)%buffer_dst(is_obc+1:ie_obc,js_obc:je_obc,1))
              endif
              if (segment%field(m)%name == 'V') then
                allocate(segment%field(m)%bt_vel(is_obc+1:ie_obc,js_obc:je_obc))
                segment%field(m)%bt_vel(:,:)=0.0
              endif
            endif
          endif
          segment%field(m)%buffer_dst(:,:,:)=0.0
        endif
        ! read source data interpolated to the current model time
        if (siz(1)==1) then
          if (OBC%brushcutter_mode) then
            allocate(tmp_buffer(1,nj_seg*2-1,segment%field(m)%nk_src))  ! segment data is currrently on supergrid
          else
            allocate(tmp_buffer(1,nj_seg,segment%field(m)%nk_src))  ! segment data is currrently on supergrid
          endif
        else
          if (OBC%brushcutter_mode) then
            allocate(tmp_buffer(ni_seg*2-1,1,segment%field(m)%nk_src))  ! segment data is currrently on supergrid
          else
            allocate(tmp_buffer(ni_seg,1,segment%field(m)%nk_src))  ! segment data is currrently on supergrid
          endif
        endif

        call time_interp_external(segment%field(m)%fid,Time, tmp_buffer)
        if (OBC%brushcutter_mode) then
          if (segment%is_E_or_W) then
            if (segment%field(m)%name == 'V' .or. segment%field(m)%name == 'DVDX') then
              segment%field(m)%buffer_src(is_obc,:,:) = &
                  tmp_buffer(1,2*(js_obc+G%jdg_offset)+1:2*(je_obc+G%jdg_offset)+1:2,:)
            else
              segment%field(m)%buffer_src(is_obc,:,:) = &
                  tmp_buffer(1,2*(js_obc+G%jdg_offset)+1:2*(je_obc+G%jdg_offset):2,:)
            endif
          else
            if (segment%field(m)%name == 'U' .or. segment%field(m)%name == 'DUDY') then
              segment%field(m)%buffer_src(:,js_obc,:) = &
                  tmp_buffer(2*(is_obc+G%idg_offset)+1:2*(ie_obc+G%idg_offset)+1:2,1,:)
            else
              segment%field(m)%buffer_src(:,js_obc,:) = &
                  tmp_buffer(2*(is_obc+G%idg_offset)+1:2*(ie_obc+G%idg_offset):2,1,:)
            endif
          endif
        else
          if (segment%is_E_or_W) then
            if (segment%field(m)%name == 'V' .or. segment%field(m)%name == 'DVDX') then
              segment%field(m)%buffer_src(is_obc,:,:)=tmp_buffer(1,js_obc+G%jdg_offset+1:je_obc+G%jdg_offset+1,:)
            else
              segment%field(m)%buffer_src(is_obc,:,:)=tmp_buffer(1,js_obc+G%jdg_offset+1:je_obc+G%jdg_offset,:)
            endif
          else
            if (segment%field(m)%name == 'U' .or. segment%field(m)%name == 'DUDY') then
              segment%field(m)%buffer_src(:,js_obc,:)=tmp_buffer(is_obc+G%idg_offset+1:ie_obc+G%idg_offset+1,1,:)
            else
              segment%field(m)%buffer_src(:,js_obc,:)=tmp_buffer(is_obc+G%idg_offset+1:ie_obc+G%idg_offset,1,:)
            endif
          endif
        endif
        if (segment%field(m)%nk_src > 1) then
          call time_interp_external(segment%field(m)%fid_dz,Time, tmp_buffer)
          if (OBC%brushcutter_mode) then
            if (segment%is_E_or_W) then
              if (segment%field(m)%name == 'V' .or. segment%field(m)%name == 'DVDX') then
                segment%field(m)%dz_src(is_obc,:,:) = &
                    tmp_buffer(1,2*(js_obc+G%jdg_offset)+1:2*(je_obc+G%jdg_offset)+1:2,:)
              else
                segment%field(m)%dz_src(is_obc,:,:) = &
                    tmp_buffer(1,2*(js_obc+G%jdg_offset)+1:2*(je_obc+G%jdg_offset):2,:)
              endif
            else
              if (segment%field(m)%name == 'U' .or. segment%field(m)%name == 'DUDY') then
                segment%field(m)%dz_src(:,js_obc,:) = &
                    tmp_buffer(2*(is_obc+G%idg_offset)+1:2*(ie_obc+G%idg_offset)+1:2,1,:)
              else
                segment%field(m)%dz_src(:,js_obc,:) = &
                    tmp_buffer(2*(is_obc+G%idg_offset)+1:2*(ie_obc+G%idg_offset):2,1,:)
              endif
            endif
          else
            if (segment%is_E_or_W) then
              if (segment%field(m)%name == 'V' .or. segment%field(m)%name == 'DVDX') then
                segment%field(m)%dz_src(is_obc,:,:)=tmp_buffer(1,js_obc+G%jdg_offset+1:je_obc+G%jdg_offset+1,:)
              else
                segment%field(m)%dz_src(is_obc,:,:)=tmp_buffer(1,js_obc+G%jdg_offset+1:je_obc+G%jdg_offset,:)
              endif
            else
              if (segment%field(m)%name == 'U' .or. segment%field(m)%name == 'DUDY') then
                segment%field(m)%dz_src(:,js_obc,:)=tmp_buffer(is_obc+G%idg_offset+1:ie_obc+G%idg_offset+1,1,:)
              else
                segment%field(m)%dz_src(:,js_obc,:)=tmp_buffer(is_obc+G%idg_offset+1:ie_obc+G%idg_offset,1,:)
              endif
            endif
          endif

          call adjustSegmentEtaToFitBathymetry(G,GV,US,segment,m)

          if (segment%is_E_or_W) then
            ishift=1
            if (segment%direction == OBC_DIRECTION_E) ishift=0
            I=is_obc
            if (segment%field(m)%name == 'V' .or. segment%field(m)%name == 'DVDX') then
              ! Do q points for the whole segment
              do J=max(js_obc,jsd),min(je_obc,jed-1)
                ! Using the h remapping approach
                ! Pretty sure we need to check for source/target grid consistency here
                segment%field(m)%buffer_dst(I,J,:)=0.0  ! initialize remap destination buffer
                if (G%mask2dCu(I,j)>0. .and. G%mask2dCu(I,j+1)>0.) then
                  h_stack(:) = 0.5*(h(i+ishift,j,:) + h(i+ishift,j+1,:))
                  call remapping_core_h(OBC%remap_CS, &
                       segment%field(m)%nk_src,segment%field(m)%dz_src(I,J,:), &
                       segment%field(m)%buffer_src(I,J,:), &
                       G%ke, h_stack, segment%field(m)%buffer_dst(I,J,:))
                elseif (G%mask2dCu(I,j)>0.) then
                  h_stack(:) = h(i+ishift,j,:)
                  call remapping_core_h(OBC%remap_CS, &
                       segment%field(m)%nk_src,segment%field(m)%dz_src(I,J,:), &
                       segment%field(m)%buffer_src(I,J,:), &
                       G%ke, h_stack, segment%field(m)%buffer_dst(I,J,:))
                elseif (G%mask2dCu(I,j+1)>0.) then
                  h_stack(:) = h(i+ishift,j+1,:)
                  call remapping_core_h(OBC%remap_CS, &
                       segment%field(m)%nk_src,segment%field(m)%dz_src(I,j,:), &
                       segment%field(m)%buffer_src(I,J,:), &
                       G%ke, h_stack, segment%field(m)%buffer_dst(I,J,:))
                endif
              enddo
            else
              do j=js_obc+1,je_obc
                ! Using the h remapping approach
                ! Pretty sure we need to check for source/target grid consistency here
                segment%field(m)%buffer_dst(I,j,:)=0.0  ! initialize remap destination buffer
                if (G%mask2dCu(I,j)>0.) then
                  net_H_src = sum( segment%field(m)%dz_src(I,j,:) )
                  net_H_int = sum( h(i+ishift,j,:) )
                  scl_fac = net_H_int / net_H_src
                  call remapping_core_h(OBC%remap_CS, &
                       segment%field(m)%nk_src, scl_fac*segment%field(m)%dz_src(I,j,:), &
                       segment%field(m)%buffer_src(I,j,:), &
                       G%ke, h(i+ishift,j,:), segment%field(m)%buffer_dst(I,j,:))
                endif
              enddo
            endif
          else
            jshift=1
            if (segment%direction == OBC_DIRECTION_N) jshift=0
            J=js_obc
            if (segment%field(m)%name == 'U' .or. segment%field(m)%name == 'DUDY') then
              ! Do q points for the whole segment
              do I=max(is_obc,isd),min(ie_obc,ied-1)
                segment%field(m)%buffer_dst(I,J,:)=0.0  ! initialize remap destination buffer
                if (G%mask2dCv(i,J)>0. .and. G%mask2dCv(i+1,J)>0.) then
              ! Using the h remapping approach
              ! Pretty sure we need to check for source/target grid consistency here
                  h_stack(:) = 0.5*(h(i,j+jshift,:) + h(i+1,j+jshift,:))
                  call remapping_core_h(OBC%remap_CS, &
                       segment%field(m)%nk_src,segment%field(m)%dz_src(I,J,:), &
                       segment%field(m)%buffer_src(I,J,:), &
                       G%ke, h_stack, segment%field(m)%buffer_dst(I,J,:))
                elseif (G%mask2dCv(i,J)>0.) then
                  h_stack(:) = h(i,j+jshift,:)
                  call remapping_core_h(OBC%remap_CS, &
                       segment%field(m)%nk_src,segment%field(m)%dz_src(I,J,:), &
                       segment%field(m)%buffer_src(I,J,:), &
                       G%ke, h_stack, segment%field(m)%buffer_dst(I,J,:))
                elseif (G%mask2dCv(i+1,J)>0.) then
                  h_stack(:) = h(i+1,j+jshift,:)
                  call remapping_core_h(OBC%remap_CS, &
                       segment%field(m)%nk_src,segment%field(m)%dz_src(I,J,:), &
                       segment%field(m)%buffer_src(I,J,:), &
                       G%ke, h_stack, segment%field(m)%buffer_dst(I,J,:))
                endif
              enddo
            else
              do i=is_obc+1,ie_obc
              ! Using the h remapping approach
              ! Pretty sure we need to check for source/target grid consistency here
                segment%field(m)%buffer_dst(i,J,:)=0.0  ! initialize remap destination buffer
                if (G%mask2dCv(i,J)>0.) then
                  net_H_src = sum( segment%field(m)%dz_src(i,J,:) )
                  net_H_int = sum( h(i,j+jshift,:) )
                  scl_fac = net_H_int / net_H_src
                  call remapping_core_h(OBC%remap_CS, &
                       segment%field(m)%nk_src,segment%field(m)%dz_src(i,J,:), &
                       segment%field(m)%buffer_src(i,J,:), &
                       G%ke, h(i,j+jshift,:), segment%field(m)%buffer_dst(i,J,:))
                endif
              enddo
            endif
          endif
        else  ! 2d data
          segment%field(m)%buffer_dst(:,:,1) = segment%field(m)%buffer_src(:,:,1)  ! initialize remap destination buffer
        endif
        deallocate(tmp_buffer)
      else ! fid <= 0 (Uniform value)
        if (.not. associated(segment%field(m)%buffer_dst)) then
          if (segment%is_E_or_W) then
            if (segment%field(m)%name == 'V') then
              allocate(segment%field(m)%buffer_dst(is_obc:ie_obc,js_obc:je_obc,G%ke))
              allocate(segment%field(m)%bt_vel(is_obc:ie_obc,js_obc:je_obc))
            elseif (segment%field(m)%name == 'U') then
              allocate(segment%field(m)%buffer_dst(is_obc:ie_obc,js_obc+1:je_obc,G%ke))
              allocate(segment%field(m)%bt_vel(is_obc:ie_obc,js_obc+1:je_obc))
            elseif (segment%field(m)%name == 'DVDX') then
              allocate(segment%field(m)%buffer_dst(is_obc:ie_obc,js_obc:je_obc,G%ke))
            elseif (segment%field(m)%name == 'SSH') then
              allocate(segment%field(m)%buffer_dst(is_obc:ie_obc,js_obc:je_obc,1))
            else
              allocate(segment%field(m)%buffer_dst(is_obc:ie_obc,js_obc+1:je_obc,G%ke))
            endif
          else
            if (segment%field(m)%name == 'U') then
              allocate(segment%field(m)%buffer_dst(is_obc:ie_obc,js_obc:je_obc,G%ke))
              allocate(segment%field(m)%bt_vel(is_obc:ie_obc,js_obc:je_obc))
            elseif (segment%field(m)%name == 'V') then
              allocate(segment%field(m)%buffer_dst(is_obc+1:ie_obc,js_obc:je_obc,G%ke))
              allocate(segment%field(m)%bt_vel(is_obc+1:ie_obc,js_obc:je_obc))
            elseif (segment%field(m)%name == 'DUDY') then
              allocate(segment%field(m)%buffer_dst(is_obc:ie_obc,js_obc:je_obc,G%ke))
            elseif (segment%field(m)%name == 'SSH') then
              allocate(segment%field(m)%buffer_dst(is_obc:ie_obc,js_obc:je_obc,1))
            else
              allocate(segment%field(m)%buffer_dst(is_obc+1:ie_obc,js_obc:je_obc,G%ke))
            endif
          endif
          segment%field(m)%buffer_dst(:,:,:) = segment%field(m)%value
          if (trim(segment%field(m)%name) == 'U' .or. trim(segment%field(m)%name) == 'V') then
            segment%field(m)%bt_vel(:,:) = segment%field(m)%value
          endif
        endif
      endif

      if (segment%field(m)%fid>0) then
        ! calculate external BT velocity and transport if needed
        if (trim(segment%field(m)%name) == 'U' .or. trim(segment%field(m)%name) == 'V') then
          if (trim(segment%field(m)%name) == 'U' .and. segment%is_E_or_W) then
            I=is_obc
            do j=js_obc+1,je_obc
              normal_trans_bt(I,j) = 0.0
              do k=1,G%ke
                segment%normal_vel(I,j,k) = US%m_s_to_L_T*segment%field(m)%buffer_dst(I,j,k)
                segment%normal_trans(I,j,k) = US%m_s_to_L_T*segment%field(m)%buffer_dst(I,j,k)*segment%h(I,j,k) * &
                          G%dyCu(I,j)
                normal_trans_bt(I,j) = normal_trans_bt(I,j) + segment%normal_trans(I,j,k)
              enddo
              segment%normal_vel_bt(I,j) = normal_trans_bt(I,j) / (max(segment%Htot(I,j),1.e-12) * G%dyCu(I,j))
              if (associated(segment%nudged_normal_vel)) segment%nudged_normal_vel(I,j,:) = segment%normal_vel(I,j,:)
            enddo
          elseif (trim(segment%field(m)%name) == 'V' .and. segment%is_N_or_S) then
            J=js_obc
            do i=is_obc+1,ie_obc
              normal_trans_bt(i,J) = 0.0
              do k=1,G%ke
                segment%normal_vel(i,J,k) = US%m_s_to_L_T*segment%field(m)%buffer_dst(i,J,k)
                segment%normal_trans(i,J,k) = US%m_s_to_L_T*segment%field(m)%buffer_dst(i,J,k)*segment%h(i,J,k) * &
                          G%dxCv(i,J)
                normal_trans_bt(i,J) = normal_trans_bt(i,J) + segment%normal_trans(i,J,k)
              enddo
              segment%normal_vel_bt(i,J) = normal_trans_bt(i,J) / (max(segment%Htot(i,J),1.e-12) * G%dxCv(i,J))
              if (associated(segment%nudged_normal_vel)) segment%nudged_normal_vel(i,J,:) = segment%normal_vel(i,J,:)
            enddo
          elseif (trim(segment%field(m)%name) == 'V' .and. segment%is_E_or_W .and. &
                  associated(segment%tangential_vel)) then
            I=is_obc
            do J=js_obc,je_obc
              do k=1,G%ke
                segment%tangential_vel(I,J,k) = US%m_s_to_L_T*segment%field(m)%buffer_dst(I,J,k)
              enddo
              if (associated(segment%nudged_tangential_vel)) &
                segment%nudged_tangential_vel(I,J,:) = segment%tangential_vel(I,J,:)
            enddo
          elseif (trim(segment%field(m)%name) == 'U' .and. segment%is_N_or_S .and. &
                  associated(segment%tangential_vel)) then
            J=js_obc
            do I=is_obc,ie_obc
              do k=1,G%ke
                segment%tangential_vel(I,J,k) = US%m_s_to_L_T*segment%field(m)%buffer_dst(I,J,k)
              enddo
              if (associated(segment%nudged_tangential_vel)) &
                segment%nudged_tangential_vel(I,J,:) = segment%tangential_vel(I,J,:)
            enddo
          endif
        elseif (trim(segment%field(m)%name) == 'DVDX' .and. segment%is_E_or_W .and. &
                associated(segment%tangential_grad)) then
          I=is_obc
          do J=js_obc,je_obc
            do k=1,G%ke
              segment%tangential_grad(I,J,k) = US%T_to_s*segment%field(m)%buffer_dst(I,J,k)
              if (associated(segment%nudged_tangential_grad)) &
                segment%nudged_tangential_grad(I,J,:) = segment%tangential_grad(I,J,:)
            enddo
          enddo
        elseif (trim(segment%field(m)%name) == 'DUDY' .and. segment%is_N_or_S .and. &
                associated(segment%tangential_grad)) then
          J=js_obc
          do I=is_obc,ie_obc
            do k=1,G%ke
              segment%tangential_grad(I,J,k) = US%T_to_s*segment%field(m)%buffer_dst(I,J,k)
              if (associated(segment%nudged_tangential_grad)) &
                segment%nudged_tangential_grad(I,J,:) = segment%tangential_grad(I,J,:)
            enddo
          enddo
        endif
      endif

      ! from this point on, data are entirely on segments - will
      ! write all segment loops as 2d loops.
      if (segment%is_E_or_W) then
        js_obc2 = js_obc+1
        is_obc2 = is_obc
      else
        js_obc2 = js_obc
        is_obc2 = is_obc+1
      endif
      if (segment%is_N_or_S) then
        is_obc2 = is_obc+1
        js_obc2 = js_obc
      else
        is_obc2 = is_obc
        js_obc2 = js_obc+1
      endif

      if (trim(segment%field(m)%name) == 'SSH') then
        do j=js_obc2,je_obc
          do i=is_obc2,ie_obc
            segment%eta(i,j) = segment%field(m)%buffer_dst(i,j,1)
          enddo
        enddo
      endif

      if (trim(segment%field(m)%name) == 'TEMP') then
        if (associated(segment%field(m)%buffer_dst)) then
          do k=1,nz; do j=js_obc2, je_obc; do i=is_obc2,ie_obc
            segment%tr_Reg%Tr(1)%t(i,j,k) = segment%field(m)%buffer_dst(i,j,k)
          enddo ; enddo ; enddo
          if (.not. segment%tr_Reg%Tr(1)%is_initialized) then
            ! if the tracer reservoir has not yet been initialized, then set to external value.
            do k=1,nz; do j=js_obc2, je_obc; do i=is_obc2,ie_obc
              segment%tr_Reg%Tr(1)%tres(i,j,k) = segment%tr_Reg%Tr(1)%t(i,j,k)
            enddo ; enddo ; enddo
            segment%tr_Reg%Tr(1)%is_initialized=.true.
          endif
        else
          segment%tr_Reg%Tr(1)%OBC_inflow_conc = segment%field(m)%value
        endif
      elseif (trim(segment%field(m)%name) == 'SALT') then
        if (associated(segment%field(m)%buffer_dst)) then
          do k=1,nz; do j=js_obc2, je_obc; do i=is_obc2,ie_obc
            segment%tr_Reg%Tr(2)%t(i,j,k) = segment%field(m)%buffer_dst(i,j,k)
          enddo ; enddo ; enddo
          if (.not. segment%tr_Reg%Tr(2)%is_initialized) then
            !if the tracer reservoir has not yet been initialized, then set to external value.
            do k=1,nz; do j=js_obc2, je_obc; do i=is_obc2,ie_obc
              segment%tr_Reg%Tr(2)%tres(i,j,k) = segment%tr_Reg%Tr(2)%t(i,j,k)
            enddo ; enddo ; enddo
            segment%tr_Reg%Tr(2)%is_initialized=.true.
          endif
        else
          segment%tr_Reg%Tr(2)%OBC_inflow_conc = segment%field(m)%value
        endif
      endif

    enddo ! end field loop
    deallocate(h_stack)
    deallocate(normal_trans_bt)

  enddo ! end segment loop

end subroutine update_OBC_segment_data

!> register open boundary objects for boundary updates.
subroutine register_OBC(name, param_file, Reg)
  character(len=32),     intent(in)  :: name        !< OBC name used for error messages
  type(param_file_type), intent(in)  :: param_file  !< file to parse for  model parameter values
  type(OBC_registry_type), pointer   :: Reg         !< pointer to the tracer registry
  integer :: nobc
  character(len=256) :: mesg    ! Message for error messages.

  if (.not. associated(Reg)) call OBC_registry_init(param_file, Reg)

  if (Reg%nobc>=MAX_FIELDS_) then
    write(mesg,'("Increase MAX_FIELDS_ in MOM_memory.h to at least ",I3," to allow for &
        &all the open boundaries being registered via register_OBC.")') Reg%nobc+1
    call MOM_error(FATAL,"MOM register_tracer: "//mesg)
  endif
  Reg%nobc = Reg%nobc + 1
  nobc     = Reg%nobc

  Reg%OB(nobc)%name = name

  if (Reg%locked) call MOM_error(FATAL, &
      "MOM register_tracer was called for variable "//trim(Reg%OB(nobc)%name)//&
      " with a locked tracer registry.")

end subroutine register_OBC

!> This routine include declares and sets the variable "version".
subroutine OBC_registry_init(param_file, Reg)
  type(param_file_type),   intent(in) :: param_file !< open file to parse for model parameters
  type(OBC_registry_type), pointer    :: Reg        !< pointer to OBC registry

  integer, save :: init_calls = 0

#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_open_boundary" ! This module's name.
  character(len=256) :: mesg    ! Message for error messages.

  if (.not.associated(Reg)) then ; allocate(Reg)
  else ; return ; endif

  ! Read all relevant parameters and write them to the model log.
! call log_version(param_file, mdl,s version, "")

  init_calls = init_calls + 1
  if (init_calls > 1) then
    write(mesg,'("OBC_registry_init called ",I3, &
      &" times with different registry pointers.")') init_calls
    if (is_root_pe()) call MOM_error(WARNING,"MOM_open_boundary"//mesg)
  endif

end subroutine OBC_registry_init

!> Add file to OBC registry.
function register_file_OBC(param_file, CS, OBC_Reg)
  type(param_file_type),    intent(in) :: param_file !< parameter file.
  type(file_OBC_CS),        pointer    :: CS         !< file control structure.
  type(OBC_registry_type),  pointer    :: OBC_Reg    !< OBC registry.
  logical                              :: register_file_OBC
  character(len=32)  :: casename = "OBC file"        !< This case's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "register_file_OBC called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  ! Register the file for boundary updates.
  call register_OBC(casename, param_file, OBC_Reg)
  register_file_OBC = .true.

end function register_file_OBC

!> Clean up the file OBC from registry.
subroutine file_OBC_end(CS)
  type(file_OBC_CS), pointer    :: CS   !< OBC file control structure.

  if (associated(CS)) then
    deallocate(CS)
  endif
end subroutine file_OBC_end

!> Initialize the segment tracer registry.
subroutine segment_tracer_registry_init(param_file, segment)
  type(param_file_type),      intent(in)      :: param_file !< open file to parse for model parameters
  type(OBC_segment_type), intent(inout)       :: segment    !<  the segment

  integer, save :: init_calls = 0

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "segment_tracer_registry_init" ! This routine's name.
  character(len=256) :: mesg    ! Message for error messages.

  if (.not.associated(segment%tr_Reg)) then
    allocate(segment%tr_Reg)
  else
    return
  endif

  init_calls = init_calls + 1

  ! Read all relevant parameters and write them to the model log.
  if (init_calls == 1) call log_version(param_file, mdl, version, "")

! Need to call once per segment with tracers...
! if (init_calls > 1) then
!   write(mesg,'("segment_tracer_registry_init called ",I3, &
!     &" times with different registry pointers.")') init_calls
!   if (is_root_pe()) call MOM_error(WARNING,"MOM_tracer"//mesg)
! endif

end subroutine segment_tracer_registry_init

subroutine register_segment_tracer(tr_ptr, param_file, GV, segment, &
                                   OBC_scalar, OBC_array)
  type(verticalGrid_type), intent(in)   :: GV         !< ocean vertical grid structure
  type(tracer_type), target             :: tr_ptr     !< A target that can be used to set a pointer to the
                                                      !! stored value of tr. This target must be
                                                      !! an enduring part of the control structure,
                                                      !! because the tracer registry will use this memory,
                                                      !! but it also means that any updates to this
                                                      !! structure in the calling module will be
                                                      !! available subsequently to the tracer registry.
  type(param_file_type), intent(in)     :: param_file !< file to parse for  model parameter values
  type(OBC_segment_type), intent(inout) :: segment    !< current segment data structure
  real, optional, intent(in)            :: OBC_scalar !< If present, use scalar value for segment tracer
                                                      !! inflow concentration.
  logical, optional, intent(in)         :: OBC_array  !< If true, use array values for segment tracer
                                                      !! inflow concentration.


! Local variables
  integer :: ntseg
  integer :: isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB
  character(len=256) :: mesg    ! Message for error messages.

  call segment_tracer_registry_init(param_file, segment)

  if (segment%tr_Reg%ntseg>=MAX_FIELDS_) then
    write(mesg,'("Increase MAX_FIELDS_ in MOM_memory.h to at least ",I3," to allow for &
        &all the tracers being registered via register_tracer.")') segment%tr_Reg%ntseg+1
    call MOM_error(FATAL,"MOM register_tracer: "//mesg)
  endif
  segment%tr_Reg%ntseg = segment%tr_Reg%ntseg + 1
  ntseg     = segment%tr_Reg%ntseg

  isd = segment%HI%isd ; ied = segment%HI%ied
  jsd = segment%HI%jsd ; jed = segment%HI%jed
  IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
  JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB

  segment%tr_Reg%Tr(ntseg)%Tr => tr_ptr
  segment%tr_Reg%Tr(ntseg)%name = tr_ptr%name

  if (segment%tr_Reg%locked) call MOM_error(FATAL, &
      "MOM register_tracer was called for variable "//trim(segment%tr_Reg%Tr(ntseg)%name)//&
      " with a locked tracer registry.")

  if (present(OBC_scalar)) segment%tr_Reg%Tr(ntseg)%OBC_inflow_conc = OBC_scalar ! initialize tracer value later
  if (present(OBC_array)) then
    if (segment%is_E_or_W) then
      allocate(segment%tr_Reg%Tr(ntseg)%t(IsdB:IedB,jsd:jed,1:GV%ke));segment%tr_Reg%Tr(ntseg)%t(:,:,:)=0.0
      allocate(segment%tr_Reg%Tr(ntseg)%tres(IsdB:IedB,jsd:jed,1:GV%ke));segment%tr_Reg%Tr(ntseg)%tres(:,:,:)=0.0
      segment%tr_Reg%Tr(ntseg)%is_initialized=.false.
    elseif (segment%is_N_or_S) then
      allocate(segment%tr_Reg%Tr(ntseg)%t(isd:ied,JsdB:JedB,1:GV%ke));segment%tr_Reg%Tr(ntseg)%t(:,:,:)=0.0
      allocate(segment%tr_Reg%Tr(ntseg)%tres(isd:ied,JsdB:JedB,1:GV%ke));segment%tr_Reg%Tr(ntseg)%tres(:,:,:)=0.0
      segment%tr_Reg%Tr(ntseg)%is_initialized=.false.
    endif
  endif

end subroutine register_segment_tracer

!> Clean up the segment tracer registry.
subroutine segment_tracer_registry_end(Reg)
  type(segment_tracer_registry_type), pointer :: Reg        !< pointer to tracer registry

! Local variables
  integer n

  if (associated(Reg)) then
    do n = 1, Reg%ntseg
      if (associated(Reg%Tr(n)%t)) deallocate(Reg%Tr(n)%t)
    enddo
    deallocate(Reg)
  endif
end subroutine segment_tracer_registry_end

subroutine register_temp_salt_segments(GV, OBC, tr_Reg, param_file)
  type(verticalGrid_type),    intent(in)    :: GV         !< ocean vertical grid structure
  type(ocean_OBC_type),       pointer       :: OBC        !< Open boundary structure
  type(tracer_registry_type), pointer       :: tr_Reg     !< Tracer registry
  type(param_file_type),      intent(in)    :: param_file !< file to parse for  model parameter values

! Local variables
  integer :: isd, ied, IsdB, IedB, jsd, jed, JsdB, JedB, nz, nf
  integer :: i, j, k, n
  character(len=32)  :: name
  type(OBC_segment_type), pointer :: segment => NULL() ! pointer to segment type list
  type(tracer_type), pointer      :: tr_ptr => NULL()

  if (.not. associated(OBC)) return

  do n=1, OBC%number_of_segments
    segment=>OBC%segment(n)
    if (.not. segment%on_pe) cycle

    if (associated(segment%tr_Reg)) &
         call MOM_error(FATAL,"register_temp_salt_segments: tracer array was previously allocated")

    name = 'temp'
    call tracer_name_lookup(tr_Reg, tr_ptr, name)
    call register_segment_tracer(tr_ptr, param_file, GV, segment, &
                                 OBC_array=segment%temp_segment_data_exists)
    name = 'salt'
    call tracer_name_lookup(tr_Reg, tr_ptr, name)
    call register_segment_tracer(tr_ptr, param_file, GV, segment, &
                                 OBC_array=segment%salt_segment_data_exists)
  enddo

end subroutine register_temp_salt_segments

subroutine fill_temp_salt_segments(G, OBC, tv)
  type(ocean_grid_type),      intent(inout) :: G          !< Ocean grid structure
  type(ocean_OBC_type),       pointer       :: OBC        !< Open boundary structure
  type(thermo_var_ptrs),      intent(inout) :: tv         !< Thermodynamics structure

! Local variables
  integer :: isd, ied, IsdB, IedB, jsd, jed, JsdB, JedB, n, nz
  integer :: i, j, k
  type(OBC_segment_type), pointer :: segment => NULL() ! pointer to segment type list

  if (.not. associated(OBC)) return
  if (.not. associated(tv%T) .and. associated(tv%S)) return
  ! Both temperature and salinity fields

  call pass_var(tv%T, G%Domain)
  call pass_var(tv%S, G%Domain)

  nz = G%ke

  do n=1, OBC%number_of_segments
    segment => OBC%segment(n)
    if (.not. segment%on_pe) cycle

    isd = segment%HI%isd ; ied = segment%HI%ied
    jsd = segment%HI%jsd ; jed = segment%HI%jed
    IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
    JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB

    ! Fill with T and S values
    if (segment%is_E_or_W) then
      I=segment%HI%IsdB
      do k=1,nz ; do j=segment%HI%jsd,segment%HI%jed
        if (segment%direction == OBC_DIRECTION_W) then
          segment%tr_Reg%Tr(1)%t(I,j,k) = tv%T(i+1,j,k)
          segment%tr_Reg%Tr(2)%t(I,j,k) = tv%S(i+1,j,k)
        else
          segment%tr_Reg%Tr(1)%t(I,j,k) = tv%T(i,j,k)
          segment%tr_Reg%Tr(2)%t(I,j,k) = tv%S(i,j,k)
        endif
      enddo ; enddo
    else
      J=segment%HI%JsdB
      do k=1,nz ; do i=segment%HI%isd,segment%HI%ied
        if (segment%direction == OBC_DIRECTION_S) then
          segment%tr_Reg%Tr(1)%t(i,J,k) = tv%T(i,j+1,k)
          segment%tr_Reg%Tr(2)%t(i,J,k) = tv%S(i,j+1,k)
        else
          segment%tr_Reg%Tr(1)%t(i,J,k) = tv%T(i,j,k)
          segment%tr_Reg%Tr(2)%t(i,J,k) = tv%S(i,j,k)
        endif
      enddo ; enddo
    endif
    segment%tr_Reg%Tr(1)%tres(:,:,:) = segment%tr_Reg%Tr(1)%t(:,:,:)
    segment%tr_Reg%Tr(2)%tres(:,:,:) = segment%tr_Reg%Tr(2)%t(:,:,:)
  enddo
  call setup_OBC_tracer_reservoirs(G, OBC)
end subroutine fill_temp_salt_segments

!> Find the region outside of all open boundary segments and
!! make sure it is set to land mask. Gonna need to know global land
!! mask as well to get it right...
subroutine mask_outside_OBCs(G, US, param_file, OBC)
  type(dyn_horgrid_type),       intent(inout) :: G          !< Ocean grid structure
  type(param_file_type),        intent(in)    :: param_file !< Parameter file handle
  type(ocean_OBC_type),         pointer       :: OBC        !< Open boundary structure
  type(unit_scale_type),        intent(in)    :: US         !< A dimensional unit scaling type

  ! Local variables
  integer :: isd, ied, IsdB, IedB, jsd, jed, JsdB, JedB, n
  integer :: i, j
  logical :: fatal_error = .False.
  real    :: min_depth
  integer, parameter :: cin = 3, cout = 4, cland = -1, cedge = -2
  character(len=256) :: mesg    ! Message for error messages.
  type(OBC_segment_type), pointer :: segment => NULL() ! pointer to segment type list
  real, allocatable, dimension(:,:) :: color, color2  ! For sorting inside from outside,
                                                      ! two different ways

  if (.not. associated(OBC)) return

  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 units="m", default=0.0, scale=US%m_to_Z, do_not_log=.true.)

  allocate(color(G%isd:G%ied, G%jsd:G%jed)) ; color = 0
  allocate(color2(G%isd:G%ied, G%jsd:G%jed)) ; color2 = 0


  ! Paint a frame around the outside.
  do j=G%jsd,G%jed
    color(G%isd,j) = cedge
    color(G%ied,j) = cedge
    color2(G%isd,j) = cedge
    color2(G%ied,j) = cedge
  enddo
  do i=G%isd,G%ied
    color(i,G%jsd) = cedge
    color(i,G%jed) = cedge
    color2(i,G%jsd) = cedge
    color2(i,G%jed) = cedge
  enddo

  ! Set color to cland in the land. Note that this is before the land
  ! mask has been initialized, set mask values based on depth.
  do j=G%jsd,G%jed
    do i=G%isd,G%ied
      if (G%bathyT(i,j) <= min_depth) color(i,j) = cland
      if (G%bathyT(i,j) <= min_depth) color2(i,j) = cland
    enddo
  enddo

  do j=G%jsd,G%jed ; do i=G%IsdB+1,G%IedB-1
    if (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_W) then
      if (color(i,j) == 0.0) color(i,j) = cout
      if (color(i+1,j) == 0.0) color(i+1,j) = cin
    elseif (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_E) then
      if (color(i,j) == 0.0) color(i,j) = cin
      if (color(i+1,j) == 0.0) color(i+1,j) = cout
    endif
  enddo ; enddo
  do J=G%JsdB+1,G%JedB-1 ; do i=G%isd,G%ied
    if (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_S) then
      if (color(i,j) == 0.0) color(i,j) = cout
      if (color(i,j+1) == 0.0) color(i,j+1) = cin
    elseif (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_N) then
      if (color(i,j) == 0.0) color(i,j) = cin
      if (color(i,j+1) == 0.0) color(i,j+1) = cout
    endif
  enddo ; enddo

  do J=G%JsdB+1,G%JedB-1 ; do i=G%isd,G%ied
    if (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_S) then
      if (color2(i,j) == 0.0) color2(i,j) = cout
      if (color2(i,j+1) == 0.0) color2(i,j+1) = cin
    elseif (OBC%segment(OBC%segnum_v(i,J))%direction == OBC_DIRECTION_N) then
      if (color2(i,j) == 0.0) color2(i,j) = cin
      if (color2(i,j+1) == 0.0) color2(i,j+1) = cout
    endif
  enddo ; enddo
  do j=G%jsd,G%jed ; do i=G%IsdB+1,G%IedB-1
    if (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_W) then
      if (color2(i,j) == 0.0) color2(i,j) = cout
      if (color2(i+1,j) == 0.0) color2(i+1,j) = cin
    elseif (OBC%segment(OBC%segnum_u(I,j))%direction == OBC_DIRECTION_E) then
      if (color2(i,j) == 0.0) color2(i,j) = cin
      if (color2(i+1,j) == 0.0) color2(i+1,j) = cout
    endif
  enddo ; enddo

  ! Do the flood fill until there are no more uncolored cells.
  call flood_fill(G, color, cin, cout, cland)
  call flood_fill2(G, color2, cin, cout, cland)

  ! Use the color to set outside to min_depth on this process.
  do j=G%jsd,G%jed ; do i=G%isd,G%ied
    if (color(i,j) /= color2(i,j)) then
      fatal_error = .True.
      write(mesg,'("MOM_open_boundary: problem with OBC segments specification at ",I5,",",I5," during\n", &
          "the masking of the outside grid points.")') i, j
      call MOM_error(WARNING,"MOM register_tracer: "//mesg, all_print=.true.)
    endif
    if (color(i,j) == cout) G%bathyT(i,j) = min_depth
  enddo ; enddo
  if (fatal_error) call MOM_error(FATAL, &
      "MOM_open_boundary: inconsistent OBC segments.")

  deallocate(color)
  deallocate(color2)
end subroutine mask_outside_OBCs

!> flood the cin, cout values
subroutine flood_fill(G, color, cin, cout, cland)
  type(dyn_horgrid_type), intent(inout) :: G      !< Ocean grid structure
  real, dimension(:,:),   intent(inout) :: color  !< For sorting inside from outside
  integer, intent(in) :: cin    !< color for inside the domain
  integer, intent(in) :: cout   !< color for outside the domain
  integer, intent(in) :: cland  !< color for inside the land mask

! Local variables
  integer :: i, j, ncount

  ncount = 1
  do while (ncount > 0)
    ncount = 0
    do j=G%jsd+1,G%jed-1
      do i=G%isd+1,G%ied-1
        if (color(i,j) == 0.0 .and. color(i-1,j) > 0.0) then
          color(i,j) = color(i-1,j)
          ncount = ncount + 1
        endif
        if (color(i,j) == 0.0 .and. color(i+1,j) > 0.0) then
          color(i,j) = color(i+1,j)
          ncount = ncount + 1
        endif
        if (color(i,j) == 0.0 .and. color(i,j-1) > 0.0) then
          color(i,j) = color(i,j-1)
          ncount = ncount + 1
        endif
        if (color(i,j) == 0.0 .and. color(i,j+1) > 0.0) then
          color(i,j) = color(i,j+1)
          ncount = ncount + 1
        endif
      enddo
    enddo
    do j=G%jed-1,G%jsd+1,-1
      do i=G%ied-1,G%isd+1,-1
        if (color(i,j) == 0.0 .and. color(i-1,j) > 0.0) then
          color(i,j) = color(i-1,j)
          ncount = ncount + 1
        endif
        if (color(i,j) == 0.0 .and. color(i+1,j) > 0.0) then
          color(i,j) = color(i+1,j)
          ncount = ncount + 1
        endif
        if (color(i,j) == 0.0 .and. color(i,j-1) > 0.0) then
          color(i,j) = color(i,j-1)
          ncount = ncount + 1
        endif
        if (color(i,j) == 0.0 .and. color(i,j+1) > 0.0) then
          color(i,j) = color(i,j+1)
          ncount = ncount + 1
        endif
      enddo
    enddo
    call pass_var(color, G%Domain)
    call sum_across_PEs(ncount)
  enddo

end subroutine flood_fill

!> flood the cin, cout values
subroutine flood_fill2(G, color, cin, cout, cland)
  type(dyn_horgrid_type), intent(inout) :: G       !< Ocean grid structure
  real, dimension(:,:),   intent(inout) :: color   !< For sorting inside from outside
  integer, intent(in) :: cin    !< color for inside the domain
  integer, intent(in) :: cout   !< color for outside the domain
  integer, intent(in) :: cland  !< color for inside the land mask

! Local variables
  integer :: i, j, ncount

  ncount = 1
  do while (ncount > 0)
    ncount = 0
    do i=G%isd+1,G%ied-1
      do j=G%jsd+1,G%jed-1
        if (color(i,j) == 0.0 .and. color(i-1,j) > 0.0) then
          color(i,j) = color(i-1,j)
          ncount = ncount + 1
        endif
        if (color(i,j) == 0.0 .and. color(i+1,j) > 0.0) then
          color(i,j) = color(i+1,j)
          ncount = ncount + 1
        endif
        if (color(i,j) == 0.0 .and. color(i,j-1) > 0.0) then
          color(i,j) = color(i,j-1)
          ncount = ncount + 1
        endif
        if (color(i,j) == 0.0 .and. color(i,j+1) > 0.0) then
          color(i,j) = color(i,j+1)
          ncount = ncount + 1
        endif
      enddo
    enddo
    do i=G%ied-1,G%isd+1,-1
      do j=G%jed-1,G%jsd+1,-1
        if (color(i,j) == 0.0 .and. color(i-1,j) > 0.0) then
          color(i,j) = color(i-1,j)
          ncount = ncount + 1
        endif
        if (color(i,j) == 0.0 .and. color(i+1,j) > 0.0) then
          color(i,j) = color(i+1,j)
          ncount = ncount + 1
        endif
        if (color(i,j) == 0.0 .and. color(i,j-1) > 0.0) then
          color(i,j) = color(i,j-1)
          ncount = ncount + 1
        endif
        if (color(i,j) == 0.0 .and. color(i,j+1) > 0.0) then
          color(i,j) = color(i,j+1)
          ncount = ncount + 1
        endif
      enddo
    enddo
    call pass_var(color, G%Domain)
    call sum_across_PEs(ncount)
  enddo

end subroutine flood_fill2

!> Register OBC segment data for restarts
subroutine open_boundary_register_restarts(HI, GV, OBC, Reg, param_file, restart_CSp, &
                                           use_temperature)
  type(hor_index_type),    intent(in) :: HI !< Horizontal indices
  type(verticalGrid_type), pointer    :: GV !< Container for vertical grid information
  type(ocean_OBC_type),    pointer    :: OBC !< OBC data structure, data intent(inout)
  type(tracer_registry_type), pointer :: Reg !< pointer to tracer registry
  type(param_file_type),   intent(in) :: param_file !< Parameter file handle
  type(MOM_restart_CS),    pointer    :: restart_CSp !< Restart structure, data intent(inout)
  logical,                 intent(in) :: use_temperature !< If true, T and S are used
  ! Local variables
  type(vardesc) :: vd
  integer       :: m, n
  character(len=100) :: mesg
  type(OBC_segment_type), pointer :: segment=>NULL()

  if (.not. associated(OBC)) &
       call MOM_error(FATAL, "open_boundary_register_restarts: Called with "//&
                      "uninitialized OBC control structure")

  if (associated(OBC%rx_normal) .or. associated(OBC%ry_normal) .or. &
      associated(OBC%rx_oblique) .or. associated(OBC%ry_oblique) .or. associated(OBC%cff_normal)) &
       call MOM_error(FATAL, "open_boundary_register_restarts: Restart "//&
                      "arrays were previously allocated")

  if (associated(OBC%tres_x) .or. associated(OBC%tres_y)) &
       call MOM_error(FATAL, "open_boundary_register_restarts: Restart "//&
                      "arrays were previously allocated")

  ! *** This is a temporary work around for restarts with OBC segments.
  ! This implementation uses 3D arrays solely for restarts. We need
  ! to be able to add 2D ( x,z or y,z ) data to restarts to avoid using
  ! so much memory and disk space. ***
  if (OBC%radiation_BCs_exist_globally) then
    allocate(OBC%rx_normal(HI%isdB:HI%iedB,HI%jsd:HI%jed,GV%ke))
    OBC%rx_normal(:,:,:) = 0.0
    vd = var_desc("rx_normal", "m s-1", "Normal Phase Speed for EW radiation OBCs", 'u', 'L')
    call register_restart_field(OBC%rx_normal, vd, .false., restart_CSp)
    allocate(OBC%ry_normal(HI%isd:HI%ied,HI%jsdB:HI%jedB,GV%ke))
    OBC%ry_normal(:,:,:) = 0.0
    vd = var_desc("ry_normal", "m s-1", "Normal Phase Speed for NS radiation OBCs", 'v', 'L')
    call register_restart_field(OBC%ry_normal, vd, .false., restart_CSp)
  endif
  if (OBC%oblique_BCs_exist_globally) then
    allocate(OBC%rx_oblique(HI%isdB:HI%iedB,HI%jsd:HI%jed,GV%ke))
    OBC%rx_oblique(:,:,:) = 0.0
    vd = var_desc("rx_oblique", "m2 s-2", "Radiation Speed Squared for EW oblique OBCs", 'u', 'L')
    call register_restart_field(OBC%rx_oblique, vd, .false., restart_CSp)
    allocate(OBC%ry_oblique(HI%isd:HI%ied,HI%jsdB:HI%jedB,GV%ke))
    OBC%ry_oblique(:,:,:) = 0.0
    vd = var_desc("ry_oblique", "m2 s-2", "Radiation Speed Squared for NS oblique OBCs", 'v', 'L')
    call register_restart_field(OBC%ry_oblique, vd, .false., restart_CSp)
    allocate(OBC%cff_normal(HI%IsdB:HI%IedB,HI%jsdB:HI%jedB,GV%ke))
    OBC%cff_normal(:,:,:) = 0.0
    vd = var_desc("cff_normal", "m2 s-2", "denominator for oblique OBCs", 'q', 'L')
    call register_restart_field(OBC%cff_normal, vd, .false., restart_CSp)
  endif

  if (Reg%ntr == 0) return
  if (.not. associated(OBC%tracer_x_reservoirs_used)) then
    OBC%ntr = Reg%ntr
    allocate(OBC%tracer_x_reservoirs_used(Reg%ntr))
    allocate(OBC%tracer_y_reservoirs_used(Reg%ntr))
    OBC%tracer_x_reservoirs_used(:) = .false.
    OBC%tracer_y_reservoirs_used(:) = .false.
    call parse_for_tracer_reservoirs(OBC, param_file, use_temperature)
  else
    ! This would be coming from user code such as DOME.
    if (OBC%ntr /= Reg%ntr) then
!        call MOM_error(FATAL, "open_boundary_regiser_restarts: Inconsistent value for ntr")
      write(mesg,'("Inconsistent values for ntr ", I8," and ",I8,".")') OBC%ntr, Reg%ntr
      call MOM_error(WARNING, 'open_boundary_register_restarts: '//mesg)
    endif
  endif

  ! Still painfully inefficient, now in four dimensions.
  if (any(OBC%tracer_x_reservoirs_used)) then
    allocate(OBC%tres_x(HI%isdB:HI%iedB,HI%jsd:HI%jed,GV%ke,OBC%ntr))
    OBC%tres_x(:,:,:,:) = 0.0
    do m=1,OBC%ntr
      if (OBC%tracer_x_reservoirs_used(m)) then
        write(mesg,'("tres_x_",I3.3)') m
        vd = var_desc(mesg,"Conc", "Tracer concentration for EW OBCs",'u','L')
        call register_restart_field(OBC%tres_x(:,:,:,m), vd, .false., restart_CSp)
      endif
    enddo
  endif
  if (any(OBC%tracer_y_reservoirs_used)) then
    allocate(OBC%tres_y(HI%isd:HI%ied,HI%jsdB:HI%jedB,GV%ke,OBC%ntr))
    OBC%tres_y(:,:,:,:) = 0.0
    do m=1,OBC%ntr
      if (OBC%tracer_y_reservoirs_used(m)) then
        write(mesg,'("tres_y_",I3.3)') m
        vd = var_desc(mesg,"Conc", "Tracer concentration for NS OBCs",'v','L')
        call register_restart_field(OBC%tres_y(:,:,:,m), vd, .false., restart_CSp)
      endif
    enddo
  endif

end subroutine open_boundary_register_restarts

!> Update the OBC tracer reservoirs after the tracers have been updated.
subroutine update_segment_tracer_reservoirs(G, GV, uhr, vhr, h, OBC, dt, Reg)
  type(ocean_grid_type),                     intent(in) :: G   !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in) :: GV  !<  Ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in) :: uhr !< accumulated volume/mass flux through
                                                               !! the zonal face [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in) :: vhr !< accumulated volume/mass flux through
                                                               !! the meridional face [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in) :: h   !< layer thickness after advection
                                                               !! [H ~> m or kg m-2]
  type(ocean_OBC_type),                      pointer    :: OBC !< Open boundary structure
  real,                                      intent(in) :: dt  !< time increment [T ~> s]
  type(tracer_registry_type),                pointer    :: Reg !< pointer to tracer registry
  ! Local variables
  type(OBC_segment_type), pointer :: segment=>NULL()
  real :: u_L_in, u_L_out ! The zonal distance moved in or out of a cell [L ~> m]
  real :: v_L_in, v_L_out ! The meridional distance moved in or out of a cell [L ~> m]
  real :: fac1            ! The denominator of the expression for tracer updates [nondim]
  integer :: i, j, k, m, n, ntr, nz
  integer :: ishift, idir, jshift, jdir

  nz = GV%ke
  ntr = Reg%ntr
  if (associated(OBC)) then ; if (OBC%OBC_pe) then ; do n=1,OBC%number_of_segments
    segment=>OBC%segment(n)
    if (.not. associated(segment%tr_Reg)) cycle
    if (segment%is_E_or_W) then
      do j=segment%HI%jsd,segment%HI%jed
        I = segment%HI%IsdB
        ! ishift+I corresponds to the nearest interior tracer cell index
        ! idir switches the sign of the flow so that positive is into the reservoir
        if (segment%direction == OBC_DIRECTION_W) then
          ishift = 1 ; idir = -1
        else
          ishift = 0 ; idir = 1
        endif
        ! Update the reservoir tracer concentration implicitly using a Backward-Euler timestep
        do m=1,ntr ; if (associated(segment%tr_Reg%Tr(m)%tres)) then ; do k=1,nz
          u_L_out = max(0.0, (idir*uhr(I,j,k))*segment%Tr_InvLscale_out / (h(i+ishift,j,k)*G%dyCu(I,j)))
          u_L_in  = min(0.0, (idir*uhr(I,j,k))*segment%Tr_InvLscale_in  / (h(i+ishift,j,k)*G%dyCu(I,j)))
          fac1 = 1.0 + (u_L_out-u_L_in)
          segment%tr_Reg%Tr(m)%tres(I,j,k) = (1.0/fac1)*(segment%tr_Reg%Tr(m)%tres(I,j,k) + &
                            (u_L_out*Reg%Tr(m)%t(I+ishift,j,k) - &
                             u_L_in*segment%tr_Reg%Tr(m)%t(I,j,k)))
          if (associated(OBC%tres_x)) OBC%tres_x(I,j,k,m) = segment%tr_Reg%Tr(m)%tres(I,j,k)
        enddo ; endif ; enddo
      enddo
    else
      do i=segment%HI%isd,segment%HI%ied
        J = segment%HI%JsdB
        ! jshift+J corresponds to the nearest interior tracer cell index
        ! jdir switches the sign of the flow so that positive is into the reservoir
        if (segment%direction == OBC_DIRECTION_S) then
          jshift = 1 ; jdir = -1
        else
          jshift = 0 ; jdir = 1
        endif
        ! Update the reservoir tracer concentration implicitly using a Backward-Euler timestep
        do m=1,ntr ; if (associated(segment%tr_Reg%Tr(m)%tres)) then ; do k=1,nz
          v_L_out = max(0.0, (jdir*vhr(i,J,k))*segment%Tr_InvLscale_out / (h(i,j+jshift,k)*G%dxCv(i,J)))
          v_L_in  = min(0.0, (jdir*vhr(i,J,k))*segment%Tr_InvLscale_in  / (h(i,j+jshift,k)*G%dxCv(i,J)))
          fac1 = 1.0 + (v_L_out-v_L_in)
          segment%tr_Reg%Tr(m)%tres(i,J,k) = (1.0/fac1)*(segment%tr_Reg%Tr(m)%tres(i,J,k) + &
                            (v_L_out*Reg%Tr(m)%t(i,J+jshift,k) - &
                             v_L_in*segment%tr_Reg%Tr(m)%t(i,J,k)))
          if (associated(OBC%tres_y)) OBC%tres_y(i,J,k,m) = segment%tr_Reg%Tr(m)%tres(i,J,k)
        enddo ; endif ; enddo
      enddo
    endif
  enddo ; endif ; endif

end subroutine update_segment_tracer_reservoirs

!> Adjust interface heights to fit the bathymetry and diagnose layer thickness.
!!
!! If the bottom most interface is below the topography then the bottom-most
!! layers are contracted to GV%Angstrom_m.
!! If the bottom most interface is above the topography then the entire column
!! is dilated (expanded) to fill the void.
!!   @remark{There is a (hard-wired) "tolerance" parameter such that the
!! criteria for adjustment must equal or exceed 10cm.}
subroutine adjustSegmentEtaToFitBathymetry(G, GV, US, segment,fld)
  type(ocean_grid_type),                      intent(in)    :: G   !< The ocean's grid structure
  type(verticalGrid_type),                    intent(in)    :: GV  !< The ocean's vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US  !< A dimensional unit scaling type
  type(OBC_segment_type),                     intent(inout) :: segment !< pointer to segment type
  integer,                                    intent(in)    :: fld  !< field index to adjust thickness
  ! Local variables
  integer :: i, j, k, is, ie, js, je, nz, contractions, dilations
  integer :: n
  real, allocatable, dimension(:,:,:) :: eta ! Segment source data interface heights, [Z -> m]
  real :: hTolerance = 0.1 !<  Tolerance to exceed adjustment criteria [Z ~> m]
  real :: hTmp, eTmp, dilate
  character(len=100) :: mesg

  hTolerance = 0.1*US%m_to_Z

  nz = size(segment%field(fld)%dz_src,3)

  if (segment%is_E_or_W) then
     ! segment thicknesses are defined at cell face centers.
     is = segment%HI%isdB ; ie = segment%HI%iedB
     js = segment%HI%jsd ; je = segment%HI%jed
  else
     is = segment%HI%isd ; ie = segment%HI%ied
     js = segment%HI%jsdB ; je = segment%HI%jedB
  endif
  allocate(eta(is:ie,js:je,nz+1))
  contractions=0; dilations=0
  do j=js,je ; do i=is,ie
    eta(i,j,1)=0.0   ! segment data are assumed to be located on a static grid
    ! For remapping calls, the entire column will be dilated
    ! by a factor equal to the ratio of the sum of the geopotential referenced
    ! source data thicknesses, and the current model thicknesses. This could be
    ! an issue to be addressed, for instance if we are placing open boundaries
    ! under ice shelf cavities.
    do k=2,nz+1
      eta(i,j,k)=eta(i,j,k-1)-segment%field(fld)%dz_src(i,j,k-1)
    enddo
    ! The normal slope at the boundary is zero by a
    ! previous call to open_boundary_impose_normal_slope
    do k=nz+1,1,-1
      if (-eta(i,j,k) > segment%Htot(i,j) + hTolerance) then
         eta(i,j,k) = -segment%Htot(i,j)
         contractions = contractions + 1
      endif
    enddo

    do k=1,nz
      ! Collapse layers to thinnest possible if the thickness less than
      ! the thinnest possible (or negative).
      if (eta(i,j,K) < (eta(i,j,K+1) + GV%Angstrom_Z)) then
         eta(i,j,K) = eta(i,j,K+1) + GV%Angstrom_Z
         segment%field(fld)%dz_src(i,j,k) = GV%Angstrom_Z
      else
         segment%field(fld)%dz_src(i,j,k) = (eta(i,j,K) - eta(i,j,K+1))
      endif
    enddo

    !   The whole column is dilated to accommodate deeper topography than
    ! the bathymetry would indicate.
    if (-eta(i,j,nz+1) < segment%Htot(i,j) - hTolerance) then
       dilations = dilations + 1
       ! expand bottom-most cell only
       eta(i,j,nz+1) = -segment%Htot(i,j)
       segment%field(fld)%dz_src(i,j,nz)= eta(i,j,nz)-eta(i,j,nz+1)
       ! if (eta(i,j,1) <= eta(i,j,nz+1)) then
       !   do k=1,nz ; segment%field(fld)%dz_src(i,j,k) = (eta(i,j,1) + G%bathyT(i,j)) / real(nz) ; enddo
       ! else
       !   dilate = (eta(i,j,1) + G%bathyT(i,j)) / (eta(i,j,1) - eta(i,j,nz+1))
       !   do k=1,nz ; segment%field(fld)%dz_src(i,j,k) = segment%field(fld)%dz_src(i,j,k) * dilate ; enddo
       ! endif
       !do k=nz,2,-1 ; eta(i,j,K) = eta(i,j,K+1) + segment%field(fld)%dz_src(i,j,k) ; enddo
    endif
    ! Now convert thicknesses to units of H.
    do k=1,nz
      segment%field(fld)%dz_src(i,j,k) = segment%field(fld)%dz_src(i,j,k)*GV%Z_to_H
    enddo
  enddo; enddo

  ! can not do communication call here since only PEs on the current segment are here

  ! call sum_across_PEs(contractions)
  ! if ((contractions > 0) .and. (is_root_pe())) then
  !    write(mesg,'("Thickness OBCs were contracted ",'// &
  !         '"to fit topography in ",I8," places.")') contractions
  !    call MOM_error(WARNING, 'adjustEtaToFitBathymetry: '//mesg)
  ! endif
  ! call sum_across_PEs(dilations)
  ! if ((dilations > 0) .and. (is_root_pe())) then
  !    write(mesg,'("Thickness OBCs were dilated ",'// &
  !         '"to fit topography in ",I8," places.")') dilations
  !    call MOM_error(WARNING, 'adjustEtaToFitBathymetry: '//mesg)
  ! endif
  deallocate(eta)



end subroutine adjustSegmentEtaToFitBathymetry

!> \namespace mom_open_boundary
!! This module implements some aspects of internal open boundary
!! conditions in MOM.
!!
!! A small fragment of the grid is shown below:
!!
!!    j+1  x ^ x ^ x   At x:  q, CoriolisBu
!!    j+1  > o > o >   At ^:  v, tauy
!!    j    x ^ x ^ x   At >:  u, taux
!!    j    > o > o >   At o:  h, bathyT, buoy, tr, T, S, Rml, ustar
!!    j-1  x ^ x ^ x
!!        i-1  i  i+1  At x & ^:
!!           i  i+1    At > & o:
!!
!! The boundaries always run through q grid points (x).

end module MOM_open_boundary
