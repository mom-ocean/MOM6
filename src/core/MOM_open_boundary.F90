!> Controls where open boundary conditions are applied
module MOM_open_boundary

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_array_transform,      only : rotate_array, rotate_array_pair
use MOM_coms,                 only : sum_across_PEs, Set_PElist, Get_PElist, PE_here, num_PEs
use MOM_cpu_clock,            only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_debugging,            only : hchksum, uvchksum, chksum
use MOM_diag_mediator,        only : diag_ctrl, time_type
use MOM_domains,              only : pass_var, pass_vector
use MOM_domains,              only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,              only : To_All, EAST_FACE, NORTH_FACE, SCALAR_PAIR, CGRID_NE, CORNER
use MOM_dyn_horgrid,          only : dyn_horgrid_type
use MOM_error_handler,        only : MOM_mesg, MOM_error, FATAL, WARNING, NOTE, is_root_pe
use MOM_file_parser,          only : get_param, log_version, param_file_type, log_param
use MOM_grid,                 only : ocean_grid_type, hor_index_type
use MOM_interface_heights,    only : thickness_to_dz
use MOM_interpolate,          only : init_external_field, time_interp_external, time_interp_external_init
use MOM_interpolate,          only : external_field
use MOM_io,                   only : slasher, field_size, file_exists, stderr, SINGLE_FILE
use MOM_io,                   only : vardesc, query_vardesc, var_desc
use MOM_obsolete_params,      only : obsolete_logical, obsolete_int, obsolete_real, obsolete_char
use MOM_regridding,           only : regridding_CS
use MOM_remapping,            only : remappingSchemesDoc, remappingDefaultScheme, remapping_CS
use MOM_remapping,            only : initialize_remapping, remapping_core_h, end_remapping
use MOM_restart,              only : register_restart_field, register_restart_pair
use MOM_restart,              only : query_initialized, set_initialized, MOM_restart_CS
use MOM_string_functions,     only : extract_word, remove_spaces, uppercase, lowercase
use MOM_tidal_forcing,        only : astro_longitudes, astro_longitudes_init, eq_phase, nodal_fu, tidal_frequency
use MOM_time_manager,         only : set_date, time_type, time_type_to_real, operator(-)
use MOM_tracer_registry,      only : tracer_type, tracer_registry_type, tracer_name_lookup
use MOM_unit_scaling,         only : unit_scale_type
use MOM_variables,            only : thermo_var_ptrs
use MOM_verticalGrid,         only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public open_boundary_apply_normal_flow
public open_boundary_config
public open_boundary_setup_vert
public open_boundary_halo_update
public open_boundary_query
public open_boundary_end
public open_boundary_impose_normal_slope
public open_boundary_impose_land_mask
public radiation_open_bdry_conds
public update_OBC_segment_data
public open_boundary_test_extern_uv
public open_boundary_test_extern_h
public open_boundary_zero_normal_flow
public parse_segment_str
public parse_segment_manifest_str
public parse_segment_data_str
public register_OBC, OBC_registry_init
public register_file_OBC, file_OBC_end
public segment_tracer_registry_init
public segment_tracer_registry_end
public register_segment_tracer
public register_temp_salt_segments
public register_obgc_segments
public fill_temp_salt_segments
public fill_obgc_segments
public set_obgc_segments_props
public setup_OBC_tracer_reservoirs
public open_boundary_register_restarts
public update_segment_tracer_reservoirs
public set_initialized_OBC_tracer_reservoirs
public update_OBC_ramp
public remap_OBC_fields
public rotate_OBC_config
public rotate_OBC_segment_fields, rotate_OBC_segment_tracer_registry
public rotate_OBC_segment_direction
public write_OBC_info, chksum_OBC_segments
public initialize_segment_data
public flood_fill
public flood_fill2

integer, parameter, public :: OBC_NONE = 0      !< Indicates the use of no open boundary
integer, parameter, public :: OBC_DIRECTION_N = 100 !< Indicates the boundary is an effective northern boundary
integer, parameter, public :: OBC_DIRECTION_S = 200 !< Indicates the boundary is an effective southern boundary
integer, parameter, public :: OBC_DIRECTION_E = 300 !< Indicates the boundary is an effective eastern boundary
integer, parameter, public :: OBC_DIRECTION_W = 400 !< Indicates the boundary is an effective western boundary
integer, parameter         :: MAX_OBC_FIELDS = 100  !< Maximum number of data fields needed for OBC segments

!> Open boundary segment data from files (mostly).
type, public :: OBC_segment_data_type
  type(external_field) :: handle            !< handle from FMS associated with segment data on disk
  type(external_field) :: dz_handle         !< handle from FMS associated with segment thicknesses on disk
  logical           :: use_IO = .false.     !< True if segment data is based on file input
  character(len=32) :: name                 !< A name identifier for the segment data.  When there is grid
                                            !! rotation, this is the name on the rotated internal grid.
  character(len=8)  :: genre                !< an identifier for the segment data
  logical           :: on_face              !< If true, this field is discretized on the OBC segment
                                            !! (velocity-point) faces, or if false it as the vorticiy points
  real              :: scale                !< A scaling factor for converting input data to
                                            !! the internal units of this field.  For salinity this would
                                            !! be in units of [S ppt-1 ~> 1]
  real, allocatable :: buffer_src(:,:,:)    !< buffer for segment data located at cell faces and on
                                            !! the original vertical grid in the internally scaled
                                            !! units for the field in question, such as [L T-1 ~> m s-1]
                                            !! for a velocity or [S ~> ppt] for salinity.
  integer           :: nk_src               !< Number of vertical levels in the source data
  real, allocatable :: dz_src(:,:,:)        !< vertical grid cell spacing of the incoming segment
                                            !! data in [Z ~> m].
  real, allocatable :: buffer_dst(:,:,:)    !< buffer src data remapped to the target vertical grid
                                            !! in the internally scaled units for the field in
                                            !! question, such as [L T-1 ~> m s-1] for a velocity or
                                            !! [S ~> ppt] for salinity.
  real              :: value                !< A constant value for the inflow concentration if not read
                                            !! from file, in the internal units of a field, such as [S ~> ppt]
                                            !! for salinity.
  real              :: resrv_lfac_in = 1.   !< The reservoir inverse length scale factor for the inward
                                            !! direction per field [nondim].  The general 1/Lscale_in is
                                            !! multiplied by this factor for a specific tracer.
  real              :: resrv_lfac_out= 1.   !< The reservoir inverse length scale factor for the outward
                                            !! direction per field [nondim].  The general 1/Lscale_out is
                                            !! multiplied by this factor for a specific tracer.
end type OBC_segment_data_type

!> Tracer on OBC segment data structure, for putting into a segment tracer registry.
type, public :: OBC_segment_tracer_type
  real, allocatable          :: t(:,:,:)              !< tracer concentration array in rescaled units,
                                                      !! like [S ~> ppt] for salinity.
  real                       :: OBC_inflow_conc = 0.0 !< tracer concentration for generic inflows in rescaled units,
                                                      !! like [S ~> ppt] for salinity.
  character(len=32)          :: name                  !< tracer name used for error messages
  type(tracer_type), pointer :: Tr => NULL()          !< metadata describing the tracer
  real, allocatable          :: tres(:,:,:)           !< tracer reservoir array in rescaled units,
                                                      !! like [S ~> ppt] for salinity.
  real                       :: scale                 !< A scaling factor for converting the units of input
                                                      !! data, like [S ppt-1 ~> 1] for salinity.
  logical                    :: is_initialized        !< reservoir values have been set when True
  integer                    :: ntr_index = -1        !< index of segment tracer in the global tracer registry
  integer                    :: fd_index = -1         !< index of segment tracer in the input fields
end type OBC_segment_tracer_type

!> Registry type for tracers on segments
type, public :: segment_tracer_registry_type
  integer                       :: ntseg = 0         !< number of registered tracer segments
  type(OBC_segment_tracer_type) :: Tr(MAX_FIELDS_)   !< array of registered tracers
  logical                       :: locked = .false.  !< New tracers may be registered if locked=.false.
                                                     !! When locked=.true.,no more tracers can be registered.
                                                     !! Not sure who should lock it or when...
end type segment_tracer_registry_type

!> Open boundary segment data structure.  Unless otherwise noted, 2-d and 3-d arrays are discretized
!! at the same position as normal velocity points in the middle of the OBC segments.
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
  logical :: open           !< Boundary is open for continuity solver, and there are no other
                            !! parameterized mass fluxes at the open boundary.
  logical :: gradient       !< Zero gradient at boundary.
  logical :: values_needed  !< Whether or not any external OBC fields are needed.
  logical :: u_values_needed      !< Whether or not external u OBC fields are needed.
  logical :: uamp_values_needed   !< Whether or not external u amplitude OBC fields are needed.
  logical :: uphase_values_needed !< Whether or not external u phase OBC fields are needed.
  logical :: v_values_needed      !< Whether or not external v OBC fields are needed.
  logical :: vamp_values_needed   !< Whether or not external v amplitude OBC fields are needed.
  logical :: vphase_values_needed !< Whether or not external v phase OBC fields are needed.
  logical :: t_values_needed!< Whether or not external T OBC fields are needed.
  logical :: s_values_needed!< Whether or not external S OBC fields are needed.
  logical :: z_values_needed!< Whether or not external zeta OBC fields are needed.
  logical :: zamp_values_needed   !< Whether or not external zeta amplitude OBC fields are needed.
  logical :: zphase_values_needed !< Whether or not external zeta phase OBC fields are needed.
  logical :: g_values_needed!< Whether or not external gradient OBC fields are needed.
  integer :: direction      !< Boundary faces one of the four directions.
  logical :: is_N_or_S      !< True if the OB is facing North or South and exists on this PE.
  logical :: is_E_or_W      !< True if the OB is facing East or West and exists on this PE.
  logical :: is_E_or_W_2    !< True if the OB is facing East or West anywhere.
  type(OBC_segment_data_type), pointer :: field(:) => NULL()  !< OBC data
  integer :: num_fields     !< number of OBC data fields (e.g. u_normal,u_parallel and eta for Flather)
  integer :: Is_obc         !< Starting local i-index of boundary segment, this may be outside of the local PE.
  integer :: Ie_obc         !< Ending local i-index of boundary segment, this may be outside of the local PE.
  integer :: Js_obc         !< Starting local j-index of boundary segment, this may be outside of the local PE.
  integer :: Je_obc         !< Ending local j-index of boundary segment, this may be outside of the local PE.
  integer :: uamp_index     !< Save where uamp is in segment%field.
  integer :: uphase_index   !< Save where uphase is in segment%field.
  integer :: vamp_index     !< Save where vamp is in segment%field.
  integer :: vphase_index   !< Save where vphase is in segment%field.
  integer :: zamp_index     !< Save where zamp is in segment%field.
  integer :: zphase_index   !< Save where zphase is in segment%field.
  real :: Velocity_nudging_timescale_in  !< Nudging timescale on inflow [T ~> s].
  real :: Velocity_nudging_timescale_out !< Nudging timescale on outflow [T ~> s].
  logical :: on_pe          !< true if any portion of the segment is located in this PE's data domain
  logical :: temp_segment_data_exists !< true if temperature data arrays are present
  logical :: salt_segment_data_exists !< true if salinity data arrays are present
  real, allocatable :: Cg(:,:)  !< The external gravity wave speed [L T-1 ~> m s-1]
                                !! at OBC-points.
  real, allocatable :: Htot(:,:)  !< The total column thickness [H ~> m or kg m-2] at OBC-points.
  real, allocatable :: dZtot(:,:) !< The total column vertical extent [Z ~> m] at OBC segment faces.
  real, allocatable :: h(:,:,:)   !< The cell thickness [H ~> m or kg m-2] at OBC segment faces
  real, allocatable :: normal_vel(:,:,:)      !< The layer velocity normal to the OB
                                              !! segment [L T-1 ~> m s-1].
  real, allocatable :: tangential_vel(:,:,:)  !< The layer velocity tangential to the OB segment
                                              !! [L T-1 ~> m s-1], discretized at the corner points.
  real, allocatable :: tangential_grad(:,:,:) !< The gradient of the velocity tangential to the OB
                                              !! segment [T-1 ~> s-1], discretized at the corner points.
  real, allocatable :: normal_trans(:,:,:)    !< The layer transport normal to the OB
                                              !! segment [H L2 T-1 ~> m3 s-1].
  real, allocatable :: normal_vel_bt(:,:)     !< The barotropic velocity normal to
                                              !! the OB segment [L T-1 ~> m s-1].
  real, allocatable :: SSH(:,:)               !< The sea-surface elevation along the
                                              !! segment [Z ~> m].
  real, allocatable :: grad_normal(:,:,:)     !< The gradient of the normal flow along the
                                              !! segment times the grid spacing [L T-1 ~> m s-1],
                                              !! with the first index being the corner-point index
                                              !! along the segment, and the second index being 1 (for
                                              !! values one point into the domain) or 2 (for values
                                              !! along the OBC itself)
  real, allocatable :: grad_tan(:,:,:)        !< The gradient of the tangential flow along the
                                              !! segment times the grid spacing [L T-1 ~> m s-1], with the
                                              !! first index being the velocity/tracer point index along the
                                              !! segment, and the second being 1 for the value 1.5 points
                                              !! inside the domain and 2 for the value half a point
                                              !! inside the domain.
  real, allocatable :: grad_gradient(:,:,:)   !< The gradient normal to the segment of the gradient
                                              !! tangetial to the segment of tangential flow along the segment
                                              !! times the grid spacing [T-1 ~> s-1], with the first
                                              !! index being the velocity/tracer point index along the segment,
                                              !! and the second being 1 for the value 2 points into the domain
                                              !! and 2 for the value 1 point into the domain.
  real, allocatable :: rx_norm_rad(:,:,:)     !< The previous normal phase speed use for EW radiation
                                              !! OBC, in grid points per timestep [nondim]
  real, allocatable :: ry_norm_rad(:,:,:)     !< The previous normal phase speed use for NS radiation
                                              !! OBC, in grid points per timestep [nondim]
  real, allocatable :: rx_norm_obl(:,:,:)     !< The previous x-direction normalized radiation coefficient
                                              !! for either EW or NS oblique OBCs [L2 T-2 ~> m2 s-2]
  real, allocatable :: ry_norm_obl(:,:,:)     !< The previous y-direction normalized radiation coefficient
                                              !! for either EW or NS oblique OBCs [L2 T-2 ~> m2 s-2]
  real, allocatable :: cff_normal(:,:,:)      !< The denominator for oblique radiation of the normal
                                              !! velocity [L2 T-2 ~> m2 s-2]
  real, allocatable :: nudged_normal_vel(:,:,:) !< The layer velocity normal to the OB segment
                                              !! that values should be nudged towards [L T-1 ~> m s-1].
  real, allocatable :: nudged_tangential_vel(:,:,:) !< The layer velocity tangential to the OB segment
                                              !! that values should be nudged towards [L T-1 ~> m s-1],
                                              !! discretized at the corner (PV) points.
  real, allocatable :: nudged_tangential_grad(:,:,:)  !< The layer dvdx or dudy towards which nudging
                                              !! can occur [T-1 ~> s-1].
  type(segment_tracer_registry_type), pointer  :: tr_Reg=> NULL()!< A pointer to the tracer registry for the segment.
  type(hor_index_type) :: HI !< Horizontal index ranges
  real :: Tr_InvLscale_out                                  !< An effective inverse length scale for restoring
                                                            !! the tracer concentration in a fictitious
                                                            !! reservoir towards interior values when flow
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
  logical :: update_OBC_seg_data = .false.            !< Is it the time for OBC segment data update for fields that
                                                      !! require less frequent update
  logical :: needs_IO_for_data = .false.              !< Is any i/o needed for OBCs on the current PE
  logical :: any_needs_IO_for_data = .false.          !< Is any i/o needed for OBCs globally
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
  logical, allocatable :: tracer_x_reservoirs_used(:) !< Dimensioned by the number of tracers, set globally,
                                                      !! true for those with x reservoirs (needed for restarts).
  logical, allocatable :: tracer_y_reservoirs_used(:) !< Dimensioned by the number of tracers, set globally,
                                                      !! true for those with y reservoirs (needed for restarts).
  integer                       :: ntr = 0            !< number of tracers
  integer :: n_tide_constituents = 0                  !< Number of tidal constituents to add to the boundary.
  logical :: add_tide_constituents = .false.          !< If true, add tidal constituents to the boundary elevation
                                                      !! and velocity. Will be set to true if n_tide_constituents > 0.
  character(len=2), allocatable, dimension(:) :: tide_names  !< Names of tidal constituents to add to the boundary data.
  real, allocatable, dimension(:) :: tide_frequencies !< Angular frequencies of chosen tidal
                                                      !! constituents [rad T-1 ~> rad s-1].
  real, allocatable, dimension(:) :: tide_eq_phases   !< Equilibrium phases of chosen tidal constituents [rad].
  real, allocatable, dimension(:) :: tide_fn          !< Amplitude modulation of boundary tides by nodal cycle [nondim].
  real, allocatable, dimension(:) :: tide_un          !< Phase modulation of boundary tides by nodal cycle [rad].
  logical :: add_eq_phase = .false.                   !< If true, add the equilibrium phase argument
                                                      !! to the specified boundary tidal phase.
  logical :: add_nodal_terms = .false.                !< If true, insert terms for the 18.6 year modulation when
                                                      !! calculating tidal boundary conditions.
  type(time_type) :: time_ref                         !< Reference date (t = 0) for tidal forcing.
  type(astro_longitudes) :: tidal_longitudes          !< Lunar and solar longitudes used to calculate tidal forcing.
  ! Properties of the segments used.
  type(OBC_segment_type), allocatable :: segment(:)   !< List of segment objects.
  ! Which segment object describes the current point.
  integer, allocatable :: segnum_u(:,:) !< The absolute value gives the segment number of any OBCs at u-points,
                                        !! while the sign indicates whether they are Eastern (> 0) or Western (< 0)
                                        !! OBCs, with 0 for velocities that are not on an OBC.
  integer, allocatable :: segnum_v(:,:) !< The absolute value gives the segment number of any OBCs at v-points,
                                        !! while the sign indicates whether they are Northern (> 0) or Southern (< 0)
                                        !! OBCs, with 0 for velocities that are not on an OBC.
  ! Keep the OBC segment properties for external BGC tracers
  type(external_tracers_segments_props), pointer :: obgc_segments_props => NULL() !< obgc segment properties
  integer :: num_obgc_tracers = 0       !< The total number of obgc tracers

  ! The following parameters are used in the baroclinic radiation code:
  real :: gamma_uv !< The relative weighting for the baroclinic radiation
                   !! velocities (or speed of characteristics) at the
                   !! new time level (1) or the running mean (0) for velocities [nondim].
                   !! Valid values range from 0 to 1, with a default of 0.3.
  real :: rx_max   !< The maximum magnitude of the baroclinic radiation velocity (or speed of
                   !! characteristics) in units of grid points per timestep [nondim].
  logical :: OBC_pe !< Is there an open boundary on this tile?
  logical :: u_OBCs_on_PE   !< True if there are any u-point OBCs on this PE, including in its halos.
  logical :: v_OBCs_on_PE   !< True if there are any v-point OBCs on this PE, including in its halos.
  logical :: v_N_OBCs_on_PE !< True if there are any northern v-point OBCs on this PE, including in its halos.
  logical :: v_S_OBCs_on_PE !< True if there are any southern v-point OBCs on this PE, including in its halos.
  logical :: u_E_OBCs_on_PE !< True if there are any eastern u-point OBCs on this PE, including in its halos.
  logical :: u_W_OBCs_on_PE !< True if there are any western u-point OBCs on this PE, including in its halos.
  !>@{ Index ranges on the local PE for the open boundary conditions in various directions
  integer :: Is_u_W_obc, Ie_u_W_obc, js_u_W_obc, je_u_W_obc
  integer :: Is_u_E_obc, Ie_u_E_obc, js_u_E_obc, je_u_E_obc
  integer :: is_v_S_obc, ie_v_S_obc, Js_v_S_obc, Je_v_S_obc
  integer :: is_v_N_obc, ie_v_N_obc, Js_v_N_obc, Je_v_N_obc
  !>@}
  type(remapping_CS), pointer :: remap_z_CS => NULL() !< ALE remapping control structure for
                                                      !! z-space data on segments
  type(remapping_CS), pointer :: remap_h_CS => NULL() !< ALE remapping control structure for
                                                      !! thickness-based fields on segments
  type(OBC_registry_type), pointer :: OBC_Reg => NULL()  !< Registry type for boundaries
  real, allocatable :: rx_normal(:,:,:)     !< Array storage for normal phase speed for EW radiation OBCs
                                            !! in units of grid points per timestep [nondim]
  real, allocatable :: ry_normal(:,:,:)     !< Array storage for normal phase speed for NS radiation OBCs
                                            !! in units of grid points per timestep [nondim]
  real, allocatable :: rx_oblique_u(:,:,:)  !< X-direction oblique boundary condition radiation speeds
                                            !! squared at u points for restarts [L2 T-2 ~> m2 s-2]
  real, allocatable :: ry_oblique_u(:,:,:)  !< Y-direction oblique boundary condition radiation speeds
                                            !! squared at u points for restarts [L2 T-2 ~> m2 s-2]
  real, allocatable :: rx_oblique_v(:,:,:)  !< X-direction oblique boundary condition radiation speeds
                                            !! squared at v points for restarts [L2 T-2 ~> m2 s-2]
  real, allocatable :: ry_oblique_v(:,:,:)  !< Y-direction oblique boundary condition radiation speeds
                                            !! squared at v points for restarts [L2 T-2 ~> m2 s-2]
  real, allocatable :: cff_normal_u(:,:,:)  !< Denominator for normalizing EW oblique boundary condition
                                            !! radiation rates at u points for restarts [L2 T-2 ~> m2 s-2]
  real, allocatable :: cff_normal_v(:,:,:)  !< Denominator for normalizing NS oblique boundary condition
                                            !! radiation rates at v points for restarts [L2 T-2 ~> m2 s-2]
  real, allocatable :: tres_x(:,:,:,:)      !< Array storage of tracer reservoirs for restarts,
                                            !! in unscaled units [conc]
  real, allocatable :: tres_y(:,:,:,:)      !< Array storage of tracer reservoirs for restarts,
                                            !! in unscaled units [conc]
  logical :: debug                         !< If true, write verbose checksums for debugging purposes.
  integer :: nk_OBC_debug = 0              !< The number of layers of OBC segment data to write out
                                           !! in full when DEBUG_OBCS is true.
  real :: silly_h  !< A silly value of thickness outside of the domain that can be used to test
                   !! the independence of the OBCs to this external data [Z ~> m].
  real :: silly_u  !< A silly value of velocity outside of the domain that can be used to test
                   !! the independence of the OBCs to this external data [L T-1 ~> m s-1].
  logical :: ramp = .false.                 !< If True, ramp from zero to the external values for SSH.
  logical :: ramping_is_activated = .false. !< True if the ramping has been initialized
  real :: ramp_timescale                    !< If ramp is True, use this timescale for ramping [T ~> s].
  real :: trunc_ramp_time                   !< If ramp is True, time after which ramp is done [T ~> s].
  real :: ramp_value                        !< If ramp is True, where we are on the ramp from
                                            !! zero to one [nondim].
  type(time_type) :: ramp_start_time        !< Time when model was started.
  integer :: remap_answer_date  !< The vintage of the order of arithmetic and expressions to use
                                !! for remapping.  Values below 20190101 recover the remapping
                                !! answers from 2018, while higher values use more robust
                                !! forms of the same remapping expressions.
  logical :: check_reconstruction !< Flag for remapping to run checks on reconstruction
  logical :: check_remapping      !< Flag for remapping to run internal checks
  logical :: force_bounds_in_subcell !< Flag for remapping to hide overshoot using bounds
  logical :: om4_remap_via_sub_cells !< If true, use the OM4 remapping algorithm
  character(40) :: remappingScheme !< String selecting the vertical remapping scheme
  type(group_pass_type) :: pass_oblique  !< Structure for group halo pass
  logical :: exterior_OBC_bug   !< If true, use incorrect form of tracers exterior to OBCs.
  logical :: hor_index_bug      !< If true, recover set of a horizontal indexing bugs in the OBC code.
  logical :: reservoir_init_bug !< If true, set the OBC tracer reservoirs at the startup of a new
                                !! run from the interior tracer concentrations regardless of
                                !! properties that may be explicitly specified for the reservoir
                                !! concentrations.
end type ocean_OBC_type

!> Control structure for open boundaries that read from files.
!! Probably lots to update here.
type, public :: file_OBC_CS ; private
  logical :: OBC_file_used = .false.     !< Placeholder for now to avoid an empty type.
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

!> Type to carry OBC information needed for setting segments for OBGC tracers
type, private :: external_tracers_segments_props
   type(external_tracers_segments_props), pointer :: next => NULL() !< pointer to the next node
   character(len=128) :: tracer_name      !< tracer name
   character(len=128) :: tracer_src_file  !< tracer source file for BC
   character(len=128) :: tracer_src_field !< name of the field in source file to extract BC
   real               :: lfac_in  !< multiplicative factor for inbound  tracer reservoir length scale [nondim]
   real               :: lfac_out !< multiplicative factor for outbound tracer reservoir length scale [nondim]
end type external_tracers_segments_props
integer :: id_clock_pass !< A CPU time clock

character(len=40)  :: mdl = "MOM_open_boundary" !< This module's name.

contains

!> Enables OBC module and reads configuration parameters
!! This routine is called from MOM_initialize_fixed which
!! occurs before the initialization of the vertical coordinate
!! and ALE_init.  Therefore segment data are not fully initialized
!! here. The remainder of the segment data are initialized in a
!! later call to update_open_boundary_data
subroutine open_boundary_config(G, US, param_file, OBC)
  type(dyn_horgrid_type),  intent(inout) :: G   !< Ocean grid structure
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file handle
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundary control structure

  ! Local variables
  integer :: l ! For looping over segments
  logical :: debug, mask_outside, reentrant_x, reentrant_y
  character(len=15) :: segment_param_str ! The run-time parameter name for each segment
  character(len=1024) :: segment_str      ! The contents (rhs) for parameter "segment_param_str"
  character(len=200) :: config1          ! String for OBC_USER_CONFIG
  real               :: Lscale_in, Lscale_out ! parameters controlling tracer values at the boundaries [L ~> m]
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  logical :: check_remapping, force_bounds_in_subcell
  logical :: enable_bugs     ! If true, the defaults for recently added bug-fix flags are set to
                             ! recreate the bugs, or if false bugs are only used if actively selected.
  logical :: debugging_tests ! If true, do additional calls resetting values to help debug the performance
                             ! of the open boundary condition code.
  logical :: om4_remap_via_sub_cells ! If true, use the OM4 remapping algorithm
  ! This include declares and sets the variable "version".
# include "version_variable.h"

  allocate(OBC)

  call get_param(param_file, mdl, "OBC_NUMBER_OF_SEGMENTS", OBC%number_of_segments, &
                 default=0, do_not_log=.true.)
  call log_version(param_file, mdl, version, &
                 "Controls where open boundaries are located, what kind of boundary condition "//&
                 "to impose, and what data to apply, if any.", &
                 all_default=(OBC%number_of_segments<=0))
  call get_param(param_file, mdl, "OBC_NUMBER_OF_SEGMENTS", OBC%number_of_segments, &
                 "The number of open boundary segments.", &
                 default=0)
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
    call get_param(param_file, mdl, "RAMP_OBCS", OBC%ramp, &
         "If true, ramps from zero to the external values over time, with"//&
         "a ramping timescale given by RAMP_TIMESCALE. Ramping SSH only so far", &
         default=.false.)
    call get_param(param_file, mdl, "OBC_RAMP_TIMESCALE", OBC%ramp_timescale, &
         "If RAMP_OBCS is true, this sets the ramping timescale.", &
         units="days", default=1.0, scale=86400.0*US%s_to_T)
    call get_param(param_file, mdl, "OBC_TIDE_N_CONSTITUENTS", OBC%n_tide_constituents, &
         "Number of tidal constituents being added to the open boundary.", &
         default=0)

    if (OBC%n_tide_constituents > 0) then
      OBC%add_tide_constituents = .true.
    else
      OBC%add_tide_constituents = .false.
    endif

    call get_param(param_file, mdl, "DEBUG", debug, default=.false.)
    call get_param(param_file, mdl, "DEBUG_OBCS", OBC%debug, &
                 "If true, do additional calls to help debug the performance "//&
                 "of the open boundary condition code.", &
                 default=.false., debuggingParam=.true.)
    if (OBC%debug .and. (num_PEs() > 1)) &
      call MOM_error(FATAL, "DEBUG_OBCS = True is currently only supported for single PE runs.")
    call get_param(param_file, mdl, "OBC_DEBUGGING_TESTS", debugging_tests, &
                 "If true, do additional calls resetting certain values to help verify the correctness "//&
                 "of the open boundary condition code.", &
                 default=.false., old_name="DEBUG_OBC", debuggingParam=.true.)
    call get_param(param_file, mdl, "NK_OBC_DEBUG", OBC%nk_OBC_debug, &
                 "The number of layers of OBC segment data to write out in full "//&
                 "when DEBUG_OBCS is true.", &
                 default=0, debuggingParam=.true., do_not_log=.not.OBC%debug)

    call get_param(param_file, mdl, "OBC_SILLY_THICK", OBC%silly_h, &
                 "A silly value of thicknesses used outside of open boundary "//&
                 "conditions for debugging.", units="m", default=0.0, scale=US%m_to_Z, &
                 do_not_log=.not.debugging_tests, debuggingParam=.true.)
    call get_param(param_file, mdl, "OBC_SILLY_VEL", OBC%silly_u, &
                 "A silly value of velocities used outside of open boundary "//&
                 "conditions for debugging.", units="m/s", default=0.0, scale=US%m_s_to_L_T, &
                 do_not_log=.not.debugging_tests, debuggingParam=.true.)
    call get_param(param_file, mdl, "ENABLE_BUGS_BY_DEFAULT", enable_bugs, &
                 default=.true., do_not_log=.true.)  ! This is logged from MOM.F90.
    call get_param(param_file, mdl, "EXTERIOR_OBC_BUG", OBC%exterior_OBC_bug, &
                 "If true, recover a bug in barotropic solver and other routines when "//&
                 "boundary contitions interior to the domain are used.", &
                 default=enable_bugs)
    call get_param(param_file, mdl, "OBC_HOR_INDEXING_BUG", OBC%hor_index_bug, &
                 "If true, recover set of a horizontal indexing bugs in the OBC code.", &
                 default=enable_bugs)
    call get_param(param_file, mdl, "OBC_RESERVOIR_INIT_BUG", OBC%reservoir_init_bug, &
                 "If true, set the OBC tracer reservoirs at the startup of a new run from the "//&
                 "interior tracer concentrations regardless of properties that may be explicitly "//&
                 "specified for the reservoir concentrations.", default=enable_bugs, do_not_log=.true.)
    reentrant_x = .false.
    call get_param(param_file, mdl, "REENTRANT_X", reentrant_x, default=.true.)
    reentrant_y = .false.
    call get_param(param_file, mdl, "REENTRANT_Y", reentrant_y, default=.false.)

    ! Allocate everything
    allocate(OBC%segment(1:OBC%number_of_segments))
    do l=1,OBC%number_of_segments
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
      OBC%segment(l)%uamp_values_needed = OBC%add_tide_constituents
      OBC%segment(l)%uphase_values_needed = OBC%add_tide_constituents
      OBC%segment(l)%v_values_needed = .false.
      OBC%segment(l)%vamp_values_needed = OBC%add_tide_constituents
      OBC%segment(l)%vphase_values_needed = OBC%add_tide_constituents
      OBC%segment(l)%t_values_needed = .false.
      OBC%segment(l)%s_values_needed = .false.
      OBC%segment(l)%z_values_needed = .false.
      OBC%segment(l)%zamp_values_needed = OBC%add_tide_constituents
      OBC%segment(l)%zphase_values_needed = OBC%add_tide_constituents
      OBC%segment(l)%g_values_needed = .false.
      OBC%segment(l)%direction = OBC_NONE
      OBC%segment(l)%is_N_or_S = .false.
      OBC%segment(l)%is_E_or_W = .false.
      OBC%segment(l)%is_E_or_W_2 = .false.
      OBC%segment(l)%Velocity_nudging_timescale_in = 0.0
      OBC%segment(l)%Velocity_nudging_timescale_out = 0.0
      OBC%segment(l)%num_fields = 0
    enddo
    allocate(OBC%segnum_u(G%IsdB:G%IedB,G%jsd:G%jed), source=0)
    allocate(OBC%segnum_v(G%isd:G%ied,G%JsdB:G%JedB), source=0)
    OBC%u_OBCs_on_PE = .false.
    OBC%v_OBCs_on_PE = .false.

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
    ! Set arrays indicating the segment number and segment direction, and also store the
    ! range of indices within which various orientations of OBCs can be found on this PE.
    call set_segnum_signs(OBC, G)

    ! Moved this earlier because time_interp_external_init needs to be called
    ! before anything that uses time_interp_external (such as initialize_segment_data)
    if (OBC%specified_u_BCs_exist_globally .or. OBC%specified_v_BCs_exist_globally .or. &
      OBC%open_u_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally) then
      ! Need this for ocean_only mode boundary interpolation.
      call time_interp_external_init()
    endif
    !    if (open_boundary_query(OBC, needs_ext_seg_data=.true.)) &
 !   call initialize_segment_data(G, OBC, param_file)

    if (open_boundary_query(OBC, apply_open_OBC=.true.)) then
      call get_param(param_file, mdl, "OBC_RADIATION_MAX", OBC%rx_max, &
                   "The maximum magnitude of the baroclinic radiation velocity (or speed of "//&
                   "characteristics), in gridpoints per timestep.  This is only "//&
                   "used if one of the open boundary segments is using Orlanski.", &
                   units="nondim", default=1.0)
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

    call get_param(param_file, mdl, "REMAPPING_SCHEME", OBC%remappingScheme, &
          default=remappingDefaultScheme, do_not_log=.true.)
    call get_param(param_file, mdl, "OBC_REMAPPING_SCHEME", OBC%remappingScheme, &
          "This sets the reconstruction scheme used "//&
          "for OBC vertical remapping for all variables. "//&
          "It can be one of the following schemes: \n"//&
          trim(remappingSchemesDoc), default=OBC%remappingScheme)
    call get_param(param_file, mdl, "FATAL_CHECK_RECONSTRUCTIONS", OBC%check_reconstruction, &
          "If true, cell-by-cell reconstructions are checked for "//&
          "consistency and if non-monotonicity or an inconsistency is "//&
          "detected then a FATAL error is issued.", default=.false.,do_not_log=.true.)
    call get_param(param_file, mdl, "FATAL_CHECK_REMAPPING", OBC%check_remapping, &
          "If true, the results of remapping are checked for "//&
          "conservation and new extrema and if an inconsistency is "//&
          "detected then a FATAL error is issued.", default=.false.,do_not_log=.true.)
    call get_param(param_file, mdl, "BRUSHCUTTER_MODE", OBC%brushcutter_mode, &
         "If true, read external OBC data on the supergrid.", &
         default=.false.)
    call get_param(param_file, mdl, "REMAP_BOUND_INTERMEDIATE_VALUES", OBC%force_bounds_in_subcell, &
          "If true, the values on the intermediate grid used for remapping "//&
          "are forced to be bounded, which might not be the case due to "//&
          "round off.", default=.false.,do_not_log=.true.)
    call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
    call get_param(param_file, mdl, "REMAPPING_ANSWER_DATE", OBC%remap_answer_date, &
                 "The vintage of the expressions and order of arithmetic to use for remapping.  "//&
                 "Values below 20190101 result in the use of older, less accurate expressions "//&
                 "that were in use at the end of 2018.  Higher values result in the use of more "//&
                 "robust and accurate forms of mathematically equivalent expressions.", &
                 default=default_answer_date)
    call get_param(param_file, mdl, "REMAPPING_USE_OM4_SUBCELLS", OBC%om4_remap_via_sub_cells, &
                   do_not_log=.true., default=.true.)

    call get_param(param_file, mdl, "OBC_REMAPPING_USE_OM4_SUBCELLS", OBC%om4_remap_via_sub_cells, &
                 "If true, use the OM4 remapping-via-subcells algorithm for neutral diffusion. "//&
                 "See REMAPPING_USE_OM4_SUBCELLS for more details. "//&
                 "We recommend setting this option to false.", default=OBC%om4_remap_via_sub_cells)

  endif ! OBC%number_of_segments > 0

    ! Safety check
  if ((OBC%open_u_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally) .and. &
       .not.G%symmetric ) call MOM_error(FATAL, &
       "MOM_open_boundary, open_boundary_config: "//&
       "Symmetric memory must be used when using Flather OBCs.")
  ! Need to do this last, because it depends on time_interp_external_init having already been called
  if (OBC%add_tide_constituents) then
    call initialize_obc_tides(OBC, US, param_file)
    ! Tide update is done within update_OBC_segment_data, so this should be true if tides are included.
    OBC%update_OBC = .true.
  endif

  if (.not.(OBC%specified_u_BCs_exist_globally .or. OBC%specified_v_BCs_exist_globally .or. &
              OBC%open_u_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally)) then
    ! No open boundaries have been requested
    call open_boundary_dealloc(OBC)
  endif

end subroutine open_boundary_config

!> Setup vertical remapping for open boundaries
subroutine open_boundary_setup_vert(GV, US, OBC)
  type(verticalGrid_type), intent(in)    :: GV  !< Container for vertical grid information
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundary control structure

  ! Local variables
  real :: dz_neglect, dz_neglect_edge ! Small thicknesses in vertical height units [Z ~> m]

  if (associated(OBC)) then
    if (OBC%number_of_segments > 0) then
      if (GV%Boussinesq .and. (OBC%remap_answer_date < 20190101)) then
        dz_neglect = US%m_to_Z * 1.0e-30 ; dz_neglect_edge = US%m_to_Z * 1.0e-10
      elseif (GV%semi_Boussinesq .and. (OBC%remap_answer_date < 20190101)) then
        dz_neglect = GV%kg_m2_to_H*GV%H_to_Z * 1.0e-30 ; dz_neglect_edge = GV%kg_m2_to_H*GV%H_to_Z * 1.0e-10
      else
        dz_neglect = GV%dZ_subroundoff ; dz_neglect_edge = GV%dZ_subroundoff
      endif
      allocate(OBC%remap_z_CS)
      call initialize_remapping(OBC%remap_z_CS, OBC%remappingScheme, boundary_extrapolation=.false., &
                 check_reconstruction=OBC%check_reconstruction, check_remapping=OBC%check_remapping, &
                 om4_remap_via_sub_cells=OBC%om4_remap_via_sub_cells, &
                 force_bounds_in_subcell=OBC%force_bounds_in_subcell, answer_date=OBC%remap_answer_date, &
                 h_neglect=dz_neglect, h_neglect_edge=dz_neglect_edge)
      allocate(OBC%remap_h_CS)
      call initialize_remapping(OBC%remap_h_CS, OBC%remappingScheme, boundary_extrapolation=.false., &
                 check_reconstruction=OBC%check_reconstruction, check_remapping=OBC%check_remapping, &
                 om4_remap_via_sub_cells=OBC%om4_remap_via_sub_cells, &
                 force_bounds_in_subcell=OBC%force_bounds_in_subcell, answer_date=OBC%remap_answer_date, &
                 h_neglect=GV%H_subroundoff, h_neglect_edge=GV%H_subroundoff)
    endif
  endif

end subroutine open_boundary_setup_vert

!> Allocate space for reading OBC data from files. It sets up the required vertical
!! remapping. In the process, it does funky stuff with the MPI processes.
subroutine initialize_segment_data(G, GV, US, OBC, PF)
  type(ocean_grid_type),        intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type),      intent(in)    :: GV  !< Container for vertical grid information
  type(unit_scale_type),        intent(in)    :: US  !< A dimensional unit scaling type
  type(ocean_OBC_type), target, intent(inout) :: OBC !< Open boundary control structure
  type(param_file_type),        intent(in)    :: PF  !< Parameter file handle

  integer :: n, m, num_fields, mm
  character(len=1024) :: segstr
  character(len=256) :: filename
  character(len=20)  :: segnam, suffix
  character(len=32)  :: fieldname
  real               :: value  ! A value that is parsed from the segment data string [various units]
  character(len=32), dimension(MAX_OBC_FIELDS) :: fields  ! segment field names
  character(len=128) :: inputdir
  type(OBC_segment_type), pointer :: segment => NULL() ! pointer to segment type list
  character(len=256) :: mesg    ! Message for error messages.
  integer, dimension(4) :: siz
  integer :: is, ie, js, je
  integer :: isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB
  integer, dimension(:), allocatable :: saved_pelist
  integer :: current_pe
  integer, dimension(1) :: single_pelist
  type(external_tracers_segments_props), pointer :: obgc_segments_props_list =>NULL()
  !will be able to dynamically switch between sub-sampling refined grid data or model grid
  integer :: IO_needs(2) ! Sums to determine global OBC data use and update patterns.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  ! There is a problem with the order of the OBC initialization
  ! with respect to ALE_init. Currently handling this by copying the
  ! param file so that I can use it later in step_MOM in order to finish
  ! initializing segments on the first step.

  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)

  if (OBC%user_BCs_set_globally) return

  ! Try this here just for the documentation. It is repeated below.
  do n=1, OBC%number_of_segments
    write(segnam,"('OBC_SEGMENT_',i3.3,'_DATA')") n
    call get_param(PF, mdl, segnam, segstr, 'OBC segment docs')
  enddo

  !< temporarily disable communication in order to read segment data independently

  allocate(saved_pelist(0:num_PEs()-1))
  call Get_PElist(saved_pelist)
  current_pe = PE_here()
  single_pelist(1) = current_pe
  call Set_PElist(single_pelist)

  do n=1, OBC%number_of_segments
    segment => OBC%segment(n)
    ! segment%values_needed is only true if this segment is on the local PE and some values need to be read.
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

    call parse_segment_manifest_str(trim(segstr), num_fields, fields)
    if (num_fields == 0) then
      call MOM_mesg('initialize_segment_data: num_fields = 0')
      cycle ! cycle to next segment
    endif

    !There are OBC%num_obgc_tracers obgc tracers are there that are not listed in param file
    segment%num_fields = num_fields + OBC%num_obgc_tracers
    allocate(segment%field(segment%num_fields))

    segment%temp_segment_data_exists = .false.
    segment%salt_segment_data_exists = .false.
!!
! CODE HERE FOR OTHER OPTIONS (CLAMPED, NUDGED,..)
!!

    isd = segment%HI%isd ; ied = segment%HI%ied
    jsd = segment%HI%jsd ; jed = segment%HI%jed
    IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
    JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB

    obgc_segments_props_list => OBC%obgc_segments_props !pointer to the head node

    do m=1,segment%num_fields
      if (m <= num_fields) then
        ! These are tracers with segments specified in MOM6 style override files
        call parse_segment_data_str(trim(segstr), m, trim(fields(m)), value, filename, fieldname)
        segment%field(m)%genre = ''
      else
        ! These are obgc tracers with segments specified by external modules.
        ! Set a flag so that these can be distinguished from native tracers as they may need
        ! extra steps for preparation and handling.
        segment%field(m)%genre = 'obgc'
        ! Query the obgc segment properties by traversing the linkedlist
        call get_obgc_segments_props(obgc_segments_props_list, fields(m), filename, fieldname, &
                                     segment%field(m)%resrv_lfac_in, segment%field(m)%resrv_lfac_out)
        ! Make sure the obgc tracer is not specified in the MOM6 param file too.
        do mm=1,num_fields
          if (trim(fields(m)) == trim(fields(mm))) then
            if (is_root_pe()) &
              call MOM_error(FATAL,"MOM_open_boundary:initialize_segment_data(): obgc tracer " //trim(fields(m))// &
                               " appears in OBC_SEGMENT_XXX_DATA string in MOM6 param file. This is not supported!")
          endif
        enddo
      endif

      segment%field(m)%name = trim(fields(m))
      ! The scale factor for tracers may also be set in register_segment_tracer, and a constant input
      ! value is rescaled there.
      segment%field(m)%scale = scale_factor_from_name(fields(m), GV, US, segment%tr_Reg)
      segment%field(m)%on_face = field_is_on_face(fields(m), segment%is_E_or_W)

      if (trim(filename) /= 'none') then
        OBC%update_OBC = .true. ! Data is assumed to be time-dependent if we are reading from file
        OBC%needs_IO_for_data = .true. ! At least one segment is using I/O for OBC data
!       segment%values_needed = .true. ! Indicates that i/o will be needed for this segment
        segment%field(m)%use_IO = .true.

        filename = trim(inputdir)//trim(filename)
        fieldname = trim(fieldname)//trim(suffix)
        call field_size(filename, fieldname, siz, no_domain=.true.)
!       if (siz(4) == 1) segment%values_needed = .false.

        if (.not.file_exists(filename)) &
          call MOM_error(FATAL," Unable to open OBC file " // trim(filename))

        if (OBC%brushcutter_mode .and. (modulo(siz(1),2) == 0 .or. modulo(siz(2),2) == 0)) then
          write(mesg,'("Brushcutter mode sizes ", I6, I6)') siz(1), siz(2)
          call MOM_error(WARNING, mesg // " " // trim(filename) // " " // trim(fieldname))
          call MOM_error(FATAL,'segment data are not on the supergrid')
        endif

        if (.not.segment%field(m)%on_face) then
          allocate(segment%field(m)%buffer_src(IsdB:IedB,JsdB:JedB,siz(3)), source=0.0)
        elseif (segment%is_E_or_W) then
          allocate(segment%field(m)%buffer_src(IsdB:IedB,jsd:jed,siz(3)), source=0.0)
        else
          allocate(segment%field(m)%buffer_src(isd:ied,JsdB:JedB,siz(3)), source=0.0)
        endif

        segment%field(m)%handle = init_external_field(trim(filename), trim(fieldname), &
                  ignore_axis_atts=.true., threading=SINGLE_FILE)
        if (siz(3) > 1) then
          if ((index(segment%field(m)%name, 'phase') > 0) .or. (index(segment%field(m)%name, 'amp') > 0)) then
            ! siz(3) is constituent for tidal variables
            call field_size(filename, 'constituent', siz, no_domain=.true.)
            ! expect third dimension to be number of constituents in MOM_input
            if (siz(3) /= OBC%n_tide_constituents .and. OBC%add_tide_constituents) then
              call MOM_error(FATAL, 'Number of constituents in input data is not '//&
                  'the same as the number specified')
            endif
          else
            ! siz(3) is depth for everything else
            fieldname = 'dz_'//trim(fieldname)
            call field_size(filename, fieldname, siz, no_domain=.true.)

            if (.not.segment%field(m)%on_face) then
              allocate(segment%field(m)%dz_src(IsdB:IedB,JsdB:JedB,siz(3)), source=0.0)
            elseif (segment%is_E_or_W) then
              allocate(segment%field(m)%dz_src(IsdB:IedB,jsd:jed,siz(3)), source=0.0)
            else
              allocate(segment%field(m)%dz_src(isd:ied,JsdB:JedB,siz(3)), source=0.0)
            endif
            segment%field(m)%dz_handle = init_external_field(trim(filename), trim(fieldname), &
                      ignore_axis_atts=.true., threading=SINGLE_FILE)
          endif
          segment%field(m)%nk_src = siz(3)
        else
          segment%field(m)%nk_src = 1
        endif

        if (segment%field(m)%name == 'TEMP') segment%temp_segment_data_exists = .true.
        if (segment%field(m)%name == 'SALT') segment%salt_segment_data_exists = .true.

      else  ! This data is not being read from a file.
        segment%field(m)%value = segment%field(m)%scale * value
        segment%field(m)%use_IO = .false.

        ! Check if this is a tidal field. If so, the number
        ! of expected constituents must be 1.
        if ((index(segment%field(m)%name, 'phase') > 0) .or. (index(segment%field(m)%name, 'amp') > 0)) then
          if (OBC%n_tide_constituents > 1 .and. OBC%add_tide_constituents) then
            call MOM_error(FATAL, 'Only one constituent is supported when specifying '//&
                'tidal boundary conditions by value rather than file.')
          endif
        endif
      endif

      ! Check on which values this field is providing.
      if (segment%field(m)%name == 'TEMP') segment%t_values_needed = .false.
      if (segment%field(m)%name == 'SALT') segment%s_values_needed = .false.
      if (segment%field(m)%name == 'U') segment%u_values_needed = .false.
      if (segment%field(m)%name == 'V') segment%v_values_needed = .false.
      if (segment%field(m)%name == 'SSH') segment%z_values_needed = .false.
      if ((segment%is_N_or_S .and. segment%field(m)%name == 'DUDY') .or. &
          (segment%is_E_or_W .and. segment%field(m)%name == 'DVDX')) segment%g_values_needed = .false.
      if (segment%field(m)%name == 'Uamp') segment%uamp_values_needed = .false.
      if (segment%field(m)%name == 'Uphase') segment%uphase_values_needed = .false.
      if (segment%field(m)%name == 'Vamp') segment%vamp_values_needed = .false.
      if (segment%field(m)%name == 'Vphase') segment%vphase_values_needed = .false.
      if (segment%field(m)%name == 'SSHamp') segment%zamp_values_needed = .false.
      if (segment%field(m)%name == 'SSHphase') segment%zphase_values_needed = .false.

      ! Store the field number for later retrievals.
      if (segment%field(m)%name == 'Uamp') segment%uamp_index = m
      if (segment%field(m)%name == 'Uphase') segment%uphase_index = m
      if (segment%field(m)%name == 'Vamp') segment%vamp_index = m
      if (segment%field(m)%name == 'Vphase') segment%vphase_index = m
      if (segment%field(m)%name == 'SSHamp') segment%zamp_index = m
      if (segment%field(m)%name == 'SSHphase') segment%zphase_index = m

    enddo

    ! Check for any values that have not been provided.
    if (segment%u_values_needed .or. segment%uamp_values_needed .or. segment%uphase_values_needed .or. &
        segment%v_values_needed .or. segment%vamp_values_needed .or. segment%vphase_values_needed .or. &
        segment%t_values_needed .or. segment%s_values_needed .or. segment%g_values_needed .or. &
        segment%z_values_needed .or. segment%zamp_values_needed .or. segment%zphase_values_needed ) then
      write(mesg,'("Values needed for OBC segment ",I3)') n
      call MOM_error(FATAL, mesg)
    endif
  enddo

  call Set_PElist(saved_pelist)

  ! Determine global IO data requirement patterns.
  IO_needs(1) = 0 ; if (OBC%needs_IO_for_data) IO_needs(1) = 1
  IO_needs(2) = 0 ; if (OBC%update_OBC) IO_needs(2) = 1
  call sum_across_PES(IO_needs, 2)
  OBC%any_needs_IO_for_data = (IO_needs(1) > 0)
  OBC%update_OBC = (IO_needs(2) > 0)

end subroutine initialize_segment_data

!> Determine whether a particular field is descretized at the normal-velocity faces of an open
!! boundary condition segment.
logical function field_is_on_face(name, is_E_or_W)
  character(len=*), intent(in) :: name       !< The OBC segment data name to interpret
  logical,          intent(in) :: is_E_or_W  !< This is true for an eastern or western open boundary condition

  field_is_on_face = .true.
  if (is_E_or_W) then
    if ((name == 'V') .or. (name == 'Vamp') .or. (name == 'Vphase') .or. (name == 'DVDX')) &
      field_is_on_face = .false.
  else
    if ((name == 'U') .or. (name == 'Uamp') .or. (name == 'Uphase') .or. (name == 'DUDY')) &
      field_is_on_face = .false.
  endif
end function field_is_on_face

!> Determine based on its name whether a particular field a barotropic tidal field, for which the
!! third dimension is the tidal constituent rather than a vertical axis
logical function field_is_tidal(name)
  character(len=*), intent(in) :: name       !< The OBC segment data name to interpret

  field_is_tidal = ((index(name, 'phase') > 0) .or. (index(name, 'amp') > 0))
end function field_is_tidal

!> This subroutine sets the sign of the OBC%segnum_u and OBC%segnum_v arrays to indicate the
!! direction of the faces - positive for logically eastern or northern OBCs and neagative
!! for logically western or southern OBCs, or zero on non-OBC points.  Also store information
!! about which orientations of OBCs ar on this PE and the range of indices within which the
!! various orientations of OBCs can be found on this PE.
subroutine set_segnum_signs(OBC, G)
  type(ocean_OBC_type),   intent(inout) :: OBC !< Open boundary control structure, perhaps on a rotated grid.
  type(dyn_horgrid_type), intent(in)    :: G   !< Ocean grid structure used by OBC

  integer :: i, j

  OBC%u_OBCs_on_PE = .false. ; OBC%v_OBCs_on_PE = .false.
  do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB
    OBC%segnum_u(I,j) = abs(OBC%segnum_u(I,j))
    if (abs(OBC%segnum_u(I,j)) > 0) then
      OBC%u_OBCs_on_PE = .true.
      if (OBC%segment(abs(OBC%segnum_u(I,j)))%direction == OBC_DIRECTION_W) &
        OBC%segnum_u(I,j) = -abs(OBC%segnum_u(I,j))
    endif
  enddo ; enddo
  do J=G%JsdB,G%JedB ; do i=G%isd,G%ied
    OBC%segnum_v(i,J) = abs(OBC%segnum_v(i,J))
    if (abs(OBC%segnum_v(i,J)) > 0) then
      OBC%v_OBCs_on_PE = .true.
      if (OBC%segment(abs(OBC%segnum_v(i,J)))%direction == OBC_DIRECTION_S) &
        OBC%segnum_v(i,J) = -abs(OBC%segnum_v(i,J))
    endif
  enddo ; enddo

  ! Determine the maximum and minimum index range for various directions of OBC points on this PE
  ! by first setting these one point outside of the wrong side of the domain.
  OBC%Is_u_W_obc = G%IedB + 1 ; OBC%Ie_u_W_obc = G%IsdB - 1
  OBC%js_u_W_obc = G%jed + 1 ; OBC%je_u_W_obc = G%jsd - 1
  OBC%Is_u_E_obc = G%IedB + 1 ; OBC%Ie_u_E_obc = G%IsdB - 1
  OBC%js_u_E_obc = G%jed + 1 ; OBC%je_u_E_obc = G%jsd - 1
  OBC%is_v_S_obc = G%ied + 1 ; OBC%ie_v_S_obc = G%isd - 1
  OBC%Js_v_S_obc = G%JedB + 1 ; OBC%Je_v_S_obc = G%JsdB - 1
  OBC%is_v_N_obc = G%ied + 1 ; OBC%ie_v_N_obc = G%isd - 1
  OBC%Js_v_N_obc = G%JedB + 1 ; OBC%Je_v_N_obc = G%JsdB - 1
  OBC%v_N_OBCs_on_PE = .false. ; OBC%v_S_OBCs_on_PE = .false.
  OBC%u_E_OBCs_on_PE = .false. ; OBC%u_W_OBCs_on_PE = .false.
  ! Note that the loop ranges are reduced because outward facing OBCs can not be applied at edge points.
  do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB-1
    if (OBC%segnum_u(I,j) < 0) then ! This point has OBC_DIRECTION_W.
      OBC%Is_u_W_obc = min(I, OBC%Is_u_W_obc) ; OBC%Ie_u_W_obc = max(I, OBC%Ie_u_W_obc)
      OBC%js_u_W_obc = min(j, OBC%js_u_W_obc) ; OBC%je_u_W_obc = max(j, OBC%je_u_W_obc)
      OBC%u_W_OBCs_on_PE = .true.
    endif
  enddo ; enddo
  do j=G%jsd,G%jed ; do I=G%IsdB+1,G%IedB
    if (OBC%segnum_u(I,j) > 0) then ! This point has OBC_DIRECTION_E.
      OBC%Is_u_E_obc = min(I, OBC%Is_u_E_obc) ; OBC%Ie_u_E_obc = max(I, OBC%Ie_u_E_obc)
      OBC%js_u_E_obc = min(j, OBC%js_u_E_obc) ; OBC%je_u_E_obc = max(j, OBC%je_u_E_obc)
      OBC%u_E_OBCs_on_PE = .true.
    endif
  enddo ; enddo
  do J=G%JsdB,G%JedB-1 ; do i=G%isd,G%ied
    if (OBC%segnum_v(i,J) < 0)  then ! This point has OBC_DIRECTION_S.
      OBC%is_v_S_obc = min(i, OBC%is_v_S_obc) ; OBC%ie_v_S_obc = max(i, OBC%ie_v_S_obc)
      OBC%Js_v_S_obc = min(J, OBC%Js_v_S_obc) ; OBC%Je_v_S_obc = max(J, OBC%Je_v_S_obc)
      OBC%v_S_OBCs_on_PE = .true.
    endif
  enddo ; enddo
  do J=G%JsdB+1,G%JedB ; do i=G%isd,G%ied
    if (OBC%segnum_v(i,J) > 0) then ! This point has OBC_DIRECTION_N.
      OBC%is_v_N_obc = min(i, OBC%is_v_N_obc) ; OBC%ie_v_N_obc = max(i, OBC%ie_v_N_obc)
      OBC%Js_v_N_obc = min(J, OBC%Js_v_N_obc) ; OBC%Je_v_N_obc = max(J, OBC%Je_v_N_obc)
      OBC%v_N_OBCs_on_PE = .true.
    endif
  enddo ; enddo

end subroutine set_segnum_signs

!> Return an appropriate dimensional scaling factor for input data based on an OBC segment data
!! name [various ~> 1], or 1 for tracers or other fields that do not match one of the specified names.
!! Note that calls to register_segment_tracer can come before or after calls to scale_factor_from_name.

real function scale_factor_from_name(name, GV, US, Tr_Reg)
  character(len=*),        intent(in) :: name  !< The OBC segment data name to interpret
  type(verticalGrid_type), intent(in) :: GV  !< Container for vertical grid information
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  type(segment_tracer_registry_type), pointer :: Tr_Reg  !< pointer to tracer registry for this segment

  integer :: m

  select case (trim(name))
    case ('U') ; scale_factor_from_name = US%m_s_to_L_T
    case ('V') ; scale_factor_from_name = US%m_s_to_L_T
    case ('Uamp') ; scale_factor_from_name = US%m_s_to_L_T
    case ('Vamp') ; scale_factor_from_name = US%m_s_to_L_T
    case ('DVDX') ; scale_factor_from_name = US%T_to_s
    case ('DUDY') ; scale_factor_from_name = US%T_to_s
    case ('SSH') ; scale_factor_from_name = US%m_to_Z
    case ('SSHamp') ; scale_factor_from_name = US%m_to_Z
    case default ; scale_factor_from_name = 1.0
  end select

  if (associated(Tr_Reg) .and. (scale_factor_from_name == 1.0)) then
    ! Check for name matches with previously registered tracers.
    do m=1,Tr_Reg%ntseg
      if (uppercase(name) == uppercase(Tr_Reg%Tr(m)%name)) then
        scale_factor_from_name = Tr_Reg%Tr(m)%scale
        exit
      endif
    enddo
  endif

end function scale_factor_from_name

!> Initize parameters and fields related to the specification of tides at open boundaries.
subroutine initialize_obc_tides(OBC, US, param_file)
  type(ocean_OBC_type), intent(inout) :: OBC  !< Open boundary control structure
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type), intent(in) :: param_file !< Parameter file handle
  integer, dimension(3) :: tide_ref_date      !< Reference date (t = 0) for tidal forcing (year, month, day).
  integer, dimension(3) :: nodal_ref_date     !< Date to calculate nodal modulation for (year, month, day).
  character(len=50) :: tide_constituent_str   !< List of tidal constituents to include on boundary.
  type(astro_longitudes) :: nodal_longitudes  !< Solar and lunar longitudes for tidal forcing
  type(time_type) :: nodal_time               !< Model time to calculate nodal modulation for.
  integer :: c                                !< Index to tidal constituent.
  logical :: tides                            !< True if astronomical tides are also used.

  call get_param(param_file, mdl, "OBC_TIDE_CONSTITUENTS", tide_constituent_str, &
      "Names of tidal constituents being added to the open boundaries.", &
      fail_if_missing=.true.)

  call get_param(param_file, mdl, "TIDES", tides, &
      "If true, apply tidal momentum forcing.", default=.false., do_not_log=.true.)

  call get_param(param_file, mdl, "TIDE_USE_EQ_PHASE", OBC%add_eq_phase, &
      "If true, add the equilibrium phase argument to the specified tidal phases.", &
      old_name="OBC_TIDE_ADD_EQ_PHASE", default=.false., do_not_log=tides)

  call get_param(param_file, mdl, "TIDE_ADD_NODAL", OBC%add_nodal_terms, &
      "If true, include 18.6 year nodal modulation in the boundary tidal forcing.", &
      old_name="OBC_TIDE_ADD_NODAL", default=.false., do_not_log=tides)

  call get_param(param_file, mdl, "TIDE_REF_DATE", tide_ref_date, &
      "Reference date to use for tidal calculations and equilibrium phase.", &
      old_name="OBC_TIDE_REF_DATE", defaults=(/0, 0, 0/), do_not_log=tides)

  call get_param(param_file, mdl, "TIDE_NODAL_REF_DATE", nodal_ref_date, &
      "Fixed reference date to use for nodal modulation of boundary tides.", &
      old_name="OBC_TIDE_NODAL_REF_DATE", defaults=(/0, 0, 0/), do_not_log=tides)

  if (.not. OBC%add_eq_phase) then
    ! If equilibrium phase argument is not added, the input phases
    ! should already be relative to the reference time.
    call MOM_mesg('OBC tidal phases will *not* be corrected with equilibrium arguments.')
  endif

  allocate(OBC%tide_names(OBC%n_tide_constituents))
  read(tide_constituent_str, *) OBC%tide_names

  ! Set reference time (t = 0) for boundary tidal forcing.
  OBC%time_ref = set_date(tide_ref_date(1), tide_ref_date(2), tide_ref_date(3), 0, 0, 0)

  ! Find relevant lunar and solar longitudes at the reference time
  if (OBC%add_eq_phase) call astro_longitudes_init(OBC%time_ref, OBC%tidal_longitudes)

  ! If the nodal correction is based on a different time, initialize that.
  ! Otherwise, it can use N from the time reference.
  if (OBC%add_nodal_terms) then
    if (sum(nodal_ref_date) /= 0) then
      ! A reference date was provided for the nodal correction
      nodal_time = set_date(nodal_ref_date(1), nodal_ref_date(2), nodal_ref_date(3), 0, 0, 0)
      call astro_longitudes_init(nodal_time, nodal_longitudes)
    elseif (OBC%add_eq_phase) then
      ! Astronomical longitudes were already calculated for use in equilibrium phases,
      ! so use nodal longitude from that.
      nodal_longitudes = OBC%tidal_longitudes
    else
      ! Tidal reference time is a required parameter, so calculate the longitudes from that.
      call astro_longitudes_init(OBC%time_ref, nodal_longitudes)
    endif
  endif

  allocate(OBC%tide_frequencies(OBC%n_tide_constituents))
  allocate(OBC%tide_eq_phases(OBC%n_tide_constituents))
  allocate(OBC%tide_fn(OBC%n_tide_constituents))
  allocate(OBC%tide_un(OBC%n_tide_constituents))

  do c=1,OBC%n_tide_constituents
    ! If tidal frequency is overridden by setting TIDE_*_FREQ, use that, otherwise use the
    ! default realistic frequency for this constituent.
    call get_param(param_file, mdl, "TIDE_"//trim(OBC%tide_names(c))//"_FREQ", OBC%tide_frequencies(c), &
        "Frequency of the "//trim(OBC%tide_names(c))//" tidal constituent. "//&
        "This is only used if TIDES and TIDE_"//trim(OBC%tide_names(c))// &
        " are true, or if OBC_TIDE_N_CONSTITUENTS > 0 and "//trim(OBC%tide_names(c))//&
        " is in OBC_TIDE_CONSTITUENTS.", &
        units="rad s-1", default=tidal_frequency(trim(OBC%tide_names(c))), scale=US%T_to_s)

    ! Find equilibrium phase if needed
    if (OBC%add_eq_phase) then
      OBC%tide_eq_phases(c) = eq_phase(trim(OBC%tide_names(c)), OBC%tidal_longitudes)
    else
      OBC%tide_eq_phases(c) = 0.0
    endif

    ! Find nodal corrections if needed
    if (OBC%add_nodal_terms) then
      call nodal_fu(trim(OBC%tide_names(c)), nodal_longitudes%N, OBC%tide_fn(c), OBC%tide_un(c))
    else
      OBC%tide_fn(c) = 1.0
      OBC%tide_un(c) = 0.0
    endif
  enddo
end subroutine initialize_obc_tides

!> Define indices for segment and store in hor_index_type
!! using global segment bounds corresponding to q-points
subroutine setup_segment_indices(G, seg, Is_obc, Ie_obc, Js_obc, Je_obc)
  type(dyn_horgrid_type), intent(in) :: G !< grid type
  type(OBC_segment_type), intent(inout) :: seg  !< Open boundary segment
  integer, intent(in) :: Is_obc !< Q-point global i-index of start of segment
  integer, intent(in) :: Ie_obc !< Q-point global i-index of end of segment
  integer, intent(in) :: Js_obc !< Q-point global j-index of start of segment
  integer, intent(in) :: Je_obc !< Q-point global j-index of end of segment
  ! Local variables
  integer :: IsgB, IegB, JsgB, JegB  ! Global corner point indices at the ends of the OBC segments
  integer :: isg, ieg, jsg, jeg

  ! Isg, Ieg will be I*_obc in global space
  if (Ie_obc < Is_obc) then
    IsgB = Ie_obc
    IegB = Is_obc
  else
    IsgB = Is_obc
    IegB = Ie_obc
  endif

  if (Je_obc < Js_obc) then
    JsgB = Je_obc
    JegB = Js_obc
  else
    JsgB = Js_obc
    JegB = Je_obc
  endif

  ! NOTE: h-points are defined along the interior of the segment q-points.
  !   For a given segment and its start and end index pairs, [IJ][se]gB, the
  !   h-cell corresponding to this pair are shown in the figure below.
  !
  ! x-x----------------x-x
  ! | |        N       | |
  ! x-x   W         E  x-x
  !   |        S         |
  ! x-x----------------x-x
  ! | |                | |
  ! x-x                x-x
  !
  ! For segment points on the west and south, h-point indices are incremented
  ! in order to move to the interior cell.

  if (Is_obc > Ie_obc) then
    ! Northern boundary
    isg = IsgB + 1
    jsg = JsgB
    ieg = IegB
    jeg = JegB
  endif

  if (Is_obc < Ie_obc) then
    ! Southern boundary
    isg = IsgB + 1
    jsg = JsgB + 1
    ieg = IegB
    jeg = JegB + 1
  endif

  if (Js_obc < Je_obc) then
    ! Eastern boundary
    isg = IsgB
    jsg = JsgB + 1
    ieg = IegB
    jeg = JegB
  endif

  if (Js_obc > Je_obc) then
    ! Western boundary
    isg = IsgB + 1
    jsg = JsgB + 1
    ieg = IegB + 1
    jeg = JegB
  endif

  ! Global space I*_obc but sorted
  seg%HI%IsgB = IsgB
  seg%HI%JegB = JegB
  seg%HI%IegB = IegB
  seg%HI%JsgB = JsgB

  seg%HI%isg = isg
  seg%HI%jsg = jsg
  seg%HI%ieg = ieg
  seg%HI%jeg = jeg

  ! Move into local index space
  IsgB = IsgB - G%idg_offset
  JsgB = JsgB - G%jdg_offset
  IegB = IegB - G%idg_offset
  JegB = JegB - G%jdg_offset

  isg = isg - G%idg_offset
  jsg = jsg - G%jdg_offset
  ieg = ieg - G%idg_offset
  jeg = jeg - G%jdg_offset

  ! This is the i-extent of the segment on this PE.
  ! The values are nonsense if the segment is not on this PE.
  seg%HI%IsdB = min(max(IsgB, G%HI%IsdB), G%HI%IedB)
  seg%HI%IedB = min(max(IegB, G%HI%IsdB), G%HI%IedB)
  seg%HI%isd = min(max(isg, G%HI%isd), G%HI%ied)
  seg%HI%ied = min(max(ieg, G%HI%isd), G%HI%ied)
  seg%HI%IscB = min(max(IsgB, G%HI%IscB), G%HI%IecB)
  seg%HI%IecB = min(max(IegB, G%HI%IscB), G%HI%IecB)
  seg%HI%isc = min(max(isg, G%HI%isc), G%HI%iec)
  seg%HI%iec = min(max(ieg, G%HI%isc), G%HI%iec)

  ! This is the j-extent of the segment on this PE.
  ! The values are nonsense if the segment is not on this PE.
  seg%HI%JsdB = min(max(JsgB, G%HI%JsdB), G%HI%JedB)
  seg%HI%JedB = min(max(JegB, G%HI%JsdB), G%HI%JedB)
  seg%HI%jsd = min(max(jsg, G%HI%jsd), G%HI%jed)
  seg%HI%jed = min(max(jeg, G%HI%jsd), G%HI%jed)
  seg%HI%JscB = min(max(JsgB, G%HI%JscB), G%HI%JecB)
  seg%HI%JecB = min(max(JegB, G%HI%JscB), G%HI%JecB)
  seg%HI%jsc = min(max(jsg, G%HI%jsc), G%HI%jec)
  seg%HI%jec = min(max(jeg, G%HI%jsc), G%HI%jec)

end subroutine setup_segment_indices

!> Parse an OBC_SEGMENT_%%% string starting with "I=" and configure placement and type of OBC accordingly
subroutine setup_u_point_obc(OBC, G, US, segment_str, l_seg, PF, reentrant_y)
  type(ocean_OBC_type),    intent(inout) :: OBC !< Open boundary control structure
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
  real, allocatable, dimension(:)  :: tnudge ! Nudging timescales [T ~> s]
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
    j = js_obc ; js_obc = je_obc ; je_obc = j
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
      OBC%segment(l_seg)%z_values_needed = .true.
      OBC%segment(l_seg)%u_values_needed = .true.
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
      OBC%segment(l_seg)%u_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED_TAN') then
      OBC%segment(l_seg)%nudged_tan = .true.
      OBC%nudged_u_BCs_exist_globally = .true.
      OBC%segment(l_seg)%v_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED_GRAD') then
      OBC%segment(l_seg)%nudged_grad = .true.
      OBC%segment(l_seg)%g_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'GRADIENT') then
      OBC%segment(l_seg)%gradient = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%open_u_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE') then
      OBC%segment(l_seg)%specified = .true.
      OBC%specified_u_BCs_exist_globally = .true. ! This avoids deallocation
      OBC%segment(l_seg)%u_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE_TAN') then
      OBC%segment(l_seg)%specified_tan = .true.
      OBC%segment(l_seg)%v_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE_GRAD') then
      OBC%segment(l_seg)%specified_grad = .true.
      OBC%segment(l_seg)%g_values_needed = .true.
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
                     fail_if_missing=.true., units="days", scale=86400.0*US%s_to_T)
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
      if (OBC%segment(l_seg)%direction == OBC_DIRECTION_W) OBC%segnum_u(I_obc,j) = -l_seg
      OBC%u_OBCs_on_PE = .true.
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
  type(ocean_OBC_type),    intent(inout) :: OBC !< Open boundary control structure
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
  real, allocatable, dimension(:)  :: tnudge ! Nudging timescales [T ~> s]

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
    i = Is_obc ; Is_obc = Ie_obc ; Ie_obc = i
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
      OBC%segment(l_seg)%z_values_needed = .true.
      OBC%segment(l_seg)%v_values_needed = .true.
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
      OBC%segment(l_seg)%v_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED_TAN') then
      OBC%segment(l_seg)%nudged_tan = .true.
      OBC%nudged_v_BCs_exist_globally = .true.
      OBC%segment(l_seg)%u_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED_GRAD') then
      OBC%segment(l_seg)%nudged_grad = .true.
      OBC%segment(l_seg)%g_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'GRADIENT') then
      OBC%segment(l_seg)%gradient = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%open_v_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE') then
      OBC%segment(l_seg)%specified = .true.
      OBC%specified_v_BCs_exist_globally = .true. ! This avoids deallocation
      OBC%segment(l_seg)%v_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE_TAN') then
      OBC%segment(l_seg)%specified_tan = .true.
      OBC%segment(l_seg)%u_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE_GRAD') then
      OBC%segment(l_seg)%specified_grad = .true.
      OBC%segment(l_seg)%g_values_needed = .true.
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
                     fail_if_missing=.true., units="days", scale=86400.0*US%s_to_T)
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
      if (OBC%segment(l_seg)%direction == OBC_DIRECTION_S) OBC%segnum_v(i,J_obc) = -l_seg
      OBC%v_OBCs_on_PE = .true.
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
  elseif (word1(1:2)=='J=') then ! Note that the file_parser uniformly expands "=" to " = "
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


!> Parse an OBC_SEGMENT_%%%_DATA string and determine its fields
subroutine parse_segment_manifest_str(segment_str, num_fields, fields)
  character(len=*), intent(in) :: segment_str   !< A string in form of
                                        !< "VAR1=file:foo1.nc(varnam1),VAR2=file:foo2.nc(varnam2),..."
  integer, intent(out) :: num_fields    !< The number of fields in the segment data
  character(len=*), dimension(MAX_OBC_FIELDS), intent(out) :: fields
                                        !< List of fieldnames for each segment

  ! Local variables
  character(len=128) :: word1, word2

  num_fields = 0
  do
    word1 = extract_word(segment_str, ',', num_fields+1)
    if (trim(word1) == '') exit
    num_fields = num_fields + 1
    word2 = extract_word(word1, '=', 1)
    fields(num_fields) = trim(word2)
  enddo
end subroutine parse_segment_manifest_str


!> Parse an OBC_SEGMENT_%%%_DATA string
subroutine parse_segment_data_str(segment_str, idx, var, value, filename, fieldname)
  character(len=*), intent(in) :: segment_str   !< A string in form of
      !! "VAR1=file:foo1.nc(varnam1),VAR2=file:foo2.nc(varnam2),..."
  integer, intent(in) :: idx                    !< Index of segment_str record
  character(len=*), intent(in) :: var           !< The name of the variable for which parameters are needed
  character(len=*), intent(out) :: filename     !< The name of the input file if using "file" method
  character(len=*), intent(out) :: fieldname    !< The name of the variable in the input file if using
                                                !! "file" method
  real, optional, intent(out)  :: value         !< A constant value if using the "value" method in various
                                                !! units but without the internal rescaling [various units]

  ! Local variables
  character(len=128) :: word1, word2, word3, method
  integer :: lword

  ! Process first word which will start with the fieldname
  word3 = extract_word(segment_str, ',', idx)
  word1 = extract_word(word3, ':', 1)
  !if (trim(word1) == '') exit
  word2 = extract_word(word1, '=', 1)
  if (trim(word2) == trim(var)) then
    method = trim(extract_word(word1, '=', 2))
    lword = len_trim(method)
    if (method(lword-3:lword) == 'file') then
      ! raise an error id filename/fieldname not in argument list
      word1 = extract_word(word3, ':', 2)
      filename = extract_word(word1, '(', 1)
      fieldname = extract_word(word1, '(', 2)
      lword = len_trim(fieldname)
      fieldname = fieldname(1:lword-1)  ! remove trailing parenth
      value = -999.
    elseif (method(lword-4:lword) == 'value') then
      filename = 'none'
      fieldname = 'none'
      word1 = extract_word(word3, ':', 2)
      lword = len_trim(word1)
      read(word1(1:lword), *, end=986, err=987) value
    endif
  endif

  return
986 call MOM_error(FATAL,'End of record while parsing segment data specification! '//trim(segment_str))
987 call MOM_error(FATAL,'Error while parsing segment data specification! '//trim(segment_str))
end subroutine parse_segment_data_str


!> Parse all the OBC_SEGMENT_%%%_DATA strings again
!! to see which need tracer reservoirs (all pes need to know).
subroutine parse_for_tracer_reservoirs(OBC, PF, use_temperature)
  type(ocean_OBC_type), target, intent(inout) :: OBC !< Open boundary control structure
  type(param_file_type),  intent(in) :: PF  !< Parameter file handle
  logical,                intent(in) :: use_temperature !< If true, T and S are used

  ! Local variables
  integer :: n,m,num_fields,na
  character(len=1024) :: segstr
  character(len=256) :: filename
  character(len=20)  :: segnam, suffix
  character(len=32)  :: fieldname
  real               :: value  ! A value that is parsed from the segment data string [various units]
  character(len=32), dimension(MAX_OBC_FIELDS) :: fields  ! segment field names
  type(OBC_segment_type), pointer :: segment => NULL() ! pointer to segment type list

  do n=1, OBC%number_of_segments
    segment => OBC%segment(n)
    write(segnam,"('OBC_SEGMENT_',i3.3,'_DATA')") n
    write(suffix,"('_segment_',i3.3)") n
    ! Clear out any old values
    segstr = ''
    call get_param(PF, mdl, segnam, segstr)
    if (segstr == '') cycle

    call parse_segment_manifest_str(trim(segstr), num_fields, fields)
    if (num_fields == 0) cycle

    ! At this point, just search for TEMP and SALT as tracers 1 and 2.
    do m=1,num_fields
      call parse_segment_data_str(trim(segstr), m, trim(fields(m)), value, filename, fieldname)
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
    !Add reservoirs for external/obgc tracers
    !There is a diconnect in the above logic between tracer index and reservoir index.
    !It arbitarily assigns reservoir indexes 1&2 to tracers T&S,
    !So we need to start from reservoir index for non-native tracers from 3, hence na=2 below.
    !num_fields is the number of vars in segstr (6 of them now,   U,V,SSH,TEMP,SALT,dye)
    !but OBC%tracer_x_reservoirs_used is allocated to size Reg%ntr, which is the total number of tracers
    na=2 !number of native MOM6 tracers (T&S) with reservoirs
    do m=1,OBC%num_obgc_tracers
       !This logic assumes all external tarcers need a reservoir
       !The segments for tracers are not initialized yet (that happens later in initialize_segment_data())
       !so we cannot query to determine if this tracer needs a reservoir.
      if (segment%is_E_or_W_2) then
        OBC%tracer_x_reservoirs_used(m+na) = .true.
      else
        OBC%tracer_y_reservoirs_used(m+na) = .true.
      endif
    enddo
  enddo

  return

end subroutine parse_for_tracer_reservoirs

!> Do any necessary halo updates on OBC-related fields.
subroutine open_boundary_halo_update(G, OBC)
  type(ocean_grid_type),   intent(in) :: G   !< Ocean grid structure
  type(ocean_OBC_type),    pointer    :: OBC !< Open boundary control structure

  ! Local variables
  integer :: m

  if (.not.associated(OBC)) return

  id_clock_pass = cpu_clock_id('(Ocean OBC halo updates)', grain=CLOCK_ROUTINE)
  if (OBC%radiation_BCs_exist_globally) call pass_vector(OBC%rx_normal, OBC%ry_normal, G%Domain, &
                     To_All+Scalar_Pair)
  if (OBC%oblique_BCs_exist_globally) then
!   call pass_vector(OBC%rx_oblique_u, OBC%ry_oblique_v, G%Domain, To_All+Scalar_Pair)
!   call pass_vector(OBC%ry_oblique_u, OBC%rx_oblique_v, G%Domain, To_All+Scalar_Pair)
!   call pass_vector(OBC%cff_normal_u, OBC%cff_normal_v, G%Domain, To_All+Scalar_Pair)
    call create_group_pass(OBC%pass_oblique, OBC%rx_oblique_u, OBC%ry_oblique_v, G%Domain, To_All+Scalar_Pair)
    call create_group_pass(OBC%pass_oblique, OBC%ry_oblique_u, OBC%rx_oblique_v, G%Domain, To_All+Scalar_Pair)
    call create_group_pass(OBC%pass_oblique, OBC%cff_normal_u, OBC%cff_normal_v, G%Domain, To_All+Scalar_Pair)
    call do_group_pass(OBC%pass_oblique, G%Domain)
  endif
  if (allocated(OBC%tres_x) .and. allocated(OBC%tres_y)) then
    do m=1,OBC%ntr
      call pass_vector(OBC%tres_x(:,:,:,m), OBC%tres_y(:,:,:,m), G%Domain, To_All+Scalar_Pair)
    enddo
  elseif (allocated(OBC%tres_x)) then
    do m=1,OBC%ntr
      call pass_var(OBC%tres_x(:,:,:,m), G%Domain, position=EAST_FACE)
    enddo
  elseif (allocated(OBC%tres_y)) then
    do m=1,OBC%ntr
      call pass_var(OBC%tres_y(:,:,:,m), G%Domain, position=NORTH_FACE)
    enddo
  endif

end subroutine open_boundary_halo_update

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
  if (present(needs_ext_seg_data)) open_boundary_query = OBC%any_needs_IO_for_data

end function open_boundary_query

!> Deallocate open boundary data
subroutine open_boundary_dealloc(OBC)
  type(ocean_OBC_type), pointer :: OBC !< Open boundary control structure
  type(OBC_segment_type), pointer :: segment => NULL()
  integer :: n

  if (.not. associated(OBC)) return

  do n=1, OBC%number_of_segments
    segment => OBC%segment(n)
    call deallocate_OBC_segment_data(segment)
  enddo
  if (allocated(OBC%segment)) deallocate(OBC%segment)
  if (allocated(OBC%segnum_u)) deallocate(OBC%segnum_u)
  if (allocated(OBC%segnum_v)) deallocate(OBC%segnum_v)
  if (allocated(OBC%rx_normal)) deallocate(OBC%rx_normal)
  if (allocated(OBC%ry_normal)) deallocate(OBC%ry_normal)
  if (allocated(OBC%rx_oblique_u)) deallocate(OBC%rx_oblique_u)
  if (allocated(OBC%ry_oblique_u)) deallocate(OBC%ry_oblique_u)
  if (allocated(OBC%rx_oblique_v)) deallocate(OBC%rx_oblique_v)
  if (allocated(OBC%ry_oblique_v)) deallocate(OBC%ry_oblique_v)
  if (allocated(OBC%cff_normal_u)) deallocate(OBC%cff_normal_u)
  if (allocated(OBC%cff_normal_v)) deallocate(OBC%cff_normal_v)
  if (allocated(OBC%tres_x)) deallocate(OBC%tres_x)
  if (allocated(OBC%tres_y)) deallocate(OBC%tres_y)
  if (associated(OBC%remap_z_CS)) deallocate(OBC%remap_z_CS)
  if (associated(OBC%remap_h_CS)) deallocate(OBC%remap_h_CS)
  deallocate(OBC)
end subroutine open_boundary_dealloc

!> Close open boundary data
subroutine open_boundary_end(OBC)
  type(ocean_OBC_type), pointer :: OBC !< Open boundary control structure
  call open_boundary_dealloc(OBC)
end subroutine open_boundary_end

!> Sets the slope of bathymetry normal to an open boundary to zero.
subroutine open_boundary_impose_normal_slope(OBC, G, depth)
  type(ocean_OBC_type),             pointer       :: OBC   !< Open boundary control structure
  type(dyn_horgrid_type),           intent(in)    :: G     !< Ocean grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: depth !< Bathymetry at h-points, in [Z ~> m] or other units
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
    segment => OBC%segment(n)
    if (.not. segment%on_pe) cycle
    if (segment%is_E_or_W) then
      ! Sweep along u-segments and delete the OBC for blocked points.
      ! Also, mask all points outside.
      I=segment%HI%IsdB
      do j=segment%HI%jsd,segment%HI%jed
        if (G%mask2dCu(I,j) == 0) OBC%segnum_u(I,j) = 0
        if (segment%direction == OBC_DIRECTION_W) then
          G%mask2dT(i,j) = 0.0
        else
          G%mask2dT(i+1,j) = 0.0
        endif
      enddo
      do J=segment%HI%JsdB+1,segment%HI%JedB-1
        if (segment%direction == OBC_DIRECTION_W) then
          G%mask2dCv(i,J) = 0 ; G%OBCmaskCv(i,J) = 0.0
        else
          G%mask2dCv(i+1,J) = 0.0 ; G%OBCmaskCv(i+1,J) = 0.0
        endif
      enddo
    else
      ! Sweep along v-segments and delete the OBC for blocked points.
      J=segment%HI%JsdB
      do i=segment%HI%isd,segment%HI%ied
        if (G%mask2dCv(i,J) == 0) OBC%segnum_v(i,J) = 0
        if (segment%direction == OBC_DIRECTION_S) then
          G%mask2dT(i,j) = 0.0
        else
          G%mask2dT(i,j+1) = 0.0
        endif
      enddo
      do I=segment%HI%IsdB+1,segment%HI%IedB-1
        if (segment%direction == OBC_DIRECTION_S) then
          G%mask2dCu(I,j) = 0.0 ; G%OBCmaskCu(I,j) = 0.0
        else
          G%mask2dCu(I,j+1) = 0.0 ; G%OBCmaskCu(I,j+1) = 0.0
        endif
      enddo
    endif
  enddo

  do n=1,OBC%number_of_segments
    segment=>OBC%segment(n)
    if (.not. (segment%on_pe .and. segment%open)) cycle
    ! Set the OBCmask values to help eliminate certain terms at u- or v- OBC points.
    if (segment%is_E_or_W) then
      I=segment%HI%IsdB
      do j=segment%HI%jsd,segment%HI%jed
        G%OBCmaskCu(I,j) = 0.0
      enddo
    else
      J=segment%HI%JsdB
      do i=segment%HI%isd,segment%HI%ied
        G%OBCmaskCv(i,J) = 0.0
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
          areaCu(I,j) = G%areaT(i,j)   ! Both of these are in [L2 ~> m2]
        else   ! West
          areaCu(I,j) = G%areaT(i+1,j) ! Both of these are in [L2 ~> m2]
        endif
      enddo
    else
      ! Sweep along v-segments and for %specified BC points reset the v-point area which was masked out
      J=segment%HI%JsdB
      do i=segment%HI%isd,segment%HI%ied
        if (segment%direction == OBC_DIRECTION_S) then
          areaCv(i,J) = G%areaT(i,j+1) ! Both of these are in [L2 ~> m2]
        else      ! North
          areaCv(i,J) = G%areaT(i,j)   ! Both of these are in [L2 ~> m2]
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
        if (OBC%segnum_u(I,j) /= 0) any_U = .true.
      enddo
    else
      J=segment%HI%JsdB
      do i=segment%HI%isd,segment%HI%ied
        if (OBC%segnum_v(i,J) /= 0) any_V = .true.
      enddo
    endif
  enddo

  OBC%u_OBCs_on_PE = any_U
  OBC%v_OBCs_on_PE = any_V
  OBC%OBC_pe = (any_U .or. any_V)

end subroutine open_boundary_impose_land_mask

!> Initialize the tracer reservoirs values, perhaps only if they have not been set via a restart file.
subroutine setup_OBC_tracer_reservoirs(G, GV, OBC, restart_CS)
  type(ocean_grid_type),          intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type),        intent(in)    :: GV  !< The ocean's vertical grid structure
  type(ocean_OBC_type), target,   intent(inout) :: OBC !< Open boundary control structure
  type(MOM_restart_CS), optional, intent(in)    :: restart_CS !< MOM restart control structure

  ! Local variables
  type(OBC_segment_type), pointer :: segment => NULL()
  real :: I_scale         ! The inverse of the scaling factor for the tracers.
                          ! For salinity the units would be [ppt S-1 ~> 1]
  logical :: set_tres_x, set_tres_y
  character(len=12) :: x_var_name, y_var_name
  integer :: i, j, k, m, n

  do m=1,OBC%ntr

    set_tres_x = allocated(OBC%tres_x) .and. OBC%tracer_x_reservoirs_used(m)
    set_tres_y = allocated(OBC%tres_y) .and. OBC%tracer_y_reservoirs_used(m)

    if (present(restart_CS)) then
      ! Set the names of the reservoirs for this tracer in the restart file, and inquire whether
      ! they have been initialized
      if (modulo(G%HI%turns, 2) == 0) then
        write(x_var_name,'("tres_x_",I3.3)') m
        write(y_var_name,'("tres_y_",I3.3)') m
      else
        write(x_var_name,'("tres_y_",I3.3)') m
        write(y_var_name,'("tres_x_",I3.3)') m
      endif
      if (set_tres_x) set_tres_x = .not.query_initialized(OBC%tres_x, x_var_name, restart_CS)
      if (set_tres_y) set_tres_y = .not.query_initialized(OBC%tres_y, y_var_name, restart_CS)
    endif

    do n=1,OBC%number_of_segments
      segment => OBC%segment(n)
      if (associated(segment%tr_Reg)) then ; if (allocated(segment%tr_Reg%Tr(m)%tres)) then
        I_scale = 1.0 ; if (segment%tr_Reg%Tr(m)%scale /= 0.0) I_scale = 1.0 / segment%tr_Reg%Tr(m)%scale

        if (segment%is_E_or_W .and. set_tres_x) then
          I = segment%HI%IsdB
          if (segment%tr_Reg%Tr(m)%is_initialized) then
            do k=1,GV%ke ; do j=segment%HI%jsd,segment%HI%jed
              OBC%tres_x(I,j,k,m) = I_scale * segment%tr_Reg%Tr(m)%tres(i,j,k)
            enddo ; enddo
          else
            do k=1,GV%ke ; do j=segment%HI%jsd,segment%HI%jed
              OBC%tres_x(I,j,k,m) = I_scale * segment%tr_Reg%Tr(m)%t(i,j,k)
            enddo ; enddo
          endif
        elseif (segment%is_N_or_S .and. set_tres_y) then
          J = segment%HI%JsdB
          if (segment%tr_Reg%Tr(m)%is_initialized) then
            do k=1,GV%ke ; do i=segment%HI%isd,segment%HI%ied
              OBC%tres_y(i,J,k,m) = I_scale * segment%tr_Reg%Tr(m)%tres(i,J,k)
            enddo ; enddo
          else
            do k=1,GV%ke ; do i=segment%HI%isd,segment%HI%ied
              OBC%tres_y(i,J,k,m) = I_scale * segment%tr_Reg%Tr(m)%t(i,J,k)
            enddo ; enddo
          endif
        endif
      endif ; endif
    enddo
  enddo

end subroutine setup_OBC_tracer_reservoirs

!> Record that the tracer reservoirs have been initialized so that their values are not reset later.
subroutine set_initialized_OBC_tracer_reservoirs(G, OBC, restart_CS)
  type(ocean_grid_type),          intent(in)    :: G   !< Ocean grid structure
  type(ocean_OBC_type),           intent(in)    :: OBC !< Open boundary control structure
  type(MOM_restart_CS),           intent(inout) :: restart_CS !< MOM restart control structure
  character(len=12) :: x_var_name, y_var_name
  integer :: i, j, k, m, n

  do m=1,OBC%ntr
    ! Set the names of the reservoirs for this tracer in the restart file
    if (modulo(G%HI%turns, 2) == 0) then
      write(x_var_name,'("tres_x_",I3.3)') m
      write(y_var_name,'("tres_y_",I3.3)') m
    else
      write(x_var_name,'("tres_y_",I3.3)') m
      write(y_var_name,'("tres_x_",I3.3)') m
    endif

    if (OBC%tracer_x_reservoirs_used(m)) call set_initialized(OBC%tres_x, x_var_name, restart_CS)
    if (OBC%tracer_y_reservoirs_used(m)) call set_initialized(OBC%tres_y, y_var_name, restart_CS)
  enddo

end subroutine set_initialized_OBC_tracer_reservoirs


!> Apply radiation conditions to 3D u,v at open boundaries
subroutine radiation_open_bdry_conds(OBC, u_new, u_old, v_new, v_old, G, GV, US, dt)
  type(ocean_grid_type),                      intent(inout) :: G     !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV    !< The ocean's vertical grid structure
  type(ocean_OBC_type),                       pointer       :: OBC   !< Open boundary control structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u_new !< On exit, new u values on open boundaries
                                                                     !! On entry, the old time-level v but including
                                                                     !! barotropic accelerations [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: u_old !< Original unadjusted u [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v_new !< On exit, new v values on open boundaries.
                                                                     !! On entry, the old time-level v but including
                                                                     !! barotropic accelerations [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: v_old !< Original unadjusted v  [L T-1 ~> m s-1]
  type(unit_scale_type),                      intent(in)    :: US    !< A dimensional unit scaling type
  real,                                       intent(in)    :: dt    !< Appropriate timestep [T ~> s]
  ! Local variables
  real :: dhdt, dhdx, dhdy  ! One-point differences in time or space [L T-1 ~> m s-1]
  real :: gamma_u, gamma_2  ! Fractional weightings of new values [nondim]
  real :: tau            ! A local nudging timescale [T ~> s]
  real :: rx_max, ry_max ! coefficients for radiation [nondim]
  real :: rx_new, rx_avg ! coefficients for radiation [nondim] or [L2 T-2 ~> m2 s-2]
  real :: ry_new, ry_avg ! coefficients for radiation [nondim] or [L2 T-2 ~> m2 s-2]
  real :: cff_new, cff_avg ! denominator in oblique [L2 T-2 ~> m2 s-2]
  real, allocatable, dimension(:,:,:) :: &
    rx_tang_rad, & ! The phase speed at u-points for tangential oblique OBCs
                   ! in units of grid points per timestep [nondim],
                   ! discretized at the corner (PV) points.
    ry_tang_rad, & ! The phase speed at v-points for tangential oblique OBCs
                   ! in units of grid points per timestep [nondim],
                   ! discretized at the corner (PV) points.
    rx_tang_obl, & ! The x-coefficient for tangential oblique OBCs [L2 T-2 ~> m2 s-2],
                   ! discretized at the corner (PV) points.
    ry_tang_obl, & ! The y-coefficient for tangential oblique OBCs [L2 T-2 ~> m2 s-2],
                   ! discretized at the corner (PV) points.
    cff_tangential ! The denominator for tangential oblique OBCs [L2 T-2 ~> m2 s-2],
                   ! discretized at the corner (PV) points.
  real :: eps      ! A small velocity squared [L2 T-2 ~> m2 s-2]
  type(OBC_segment_type), pointer :: segment => NULL()
  integer :: i, j, k, is, ie, js, je, m, nz, n
  integer :: is_obc, ie_obc, js_obc, je_obc
  logical :: sym
  character(len=3) :: var_num

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(OBC)) return

  if (.not.(OBC%open_u_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally)) &
    return

  if (OBC%debug) call chksum_OBC_segments(OBC, G, GV, US, OBC%nk_OBC_debug)


  eps = 1.0e-20*US%m_s_to_L_T**2

  !! Copy previously calculated phase velocity from global arrays into segments
  !! This is terribly inefficient and temporary solution for continuity across restarts
  !! and needs to be revisited in the future.
  if (OBC%gamma_uv < 1.0) then
    do n=1,OBC%number_of_segments
      segment=>OBC%segment(n)
      if (.not. segment%on_pe) cycle
      if (segment%is_E_or_W .and. segment%radiation) then
        do k=1,GV%ke
          I=segment%HI%IsdB
          do j=segment%HI%jsd,segment%HI%jed
            segment%rx_norm_rad(I,j,k) = OBC%rx_normal(I,j,k)
          enddo
        enddo
      elseif (segment%is_N_or_S .and. segment%radiation) then
        do k=1,GV%ke
          J=segment%HI%JsdB
          do i=segment%HI%isd,segment%HI%ied
            segment%ry_norm_rad(i,J,k) = OBC%ry_normal(i,J,k)
          enddo
        enddo
      endif
      if (segment%is_E_or_W .and. segment%oblique) then
        do k=1,GV%ke
          I=segment%HI%IsdB
          do j=segment%HI%jsd,segment%HI%jed
            segment%rx_norm_obl(I,j,k) = OBC%rx_oblique_u(I,j,k)
            segment%ry_norm_obl(I,j,k) = OBC%ry_oblique_u(I,j,k)
            segment%cff_normal(I,j,k) = OBC%cff_normal_u(I,j,k)
          enddo
        enddo
      elseif (segment%is_N_or_S .and. segment%oblique) then
        do k=1,GV%ke
          J=segment%HI%JsdB
          do i=segment%HI%isd,segment%HI%ied
            segment%rx_norm_obl(i,J,k) = OBC%rx_oblique_v(i,J,k)
            segment%ry_norm_obl(i,J,k) = OBC%ry_oblique_v(i,J,k)
            segment%cff_normal(i,J,k) = OBC%cff_normal_v(i,J,k)
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
          if (allocated(segment%tr_Reg%Tr(m)%tres)) then
            do k=1,GV%ke
              do j=segment%HI%jsd,segment%HI%jed
                segment%tr_Reg%Tr(m)%tres(I,j,k) = segment%tr_Reg%Tr(m)%scale * OBC%tres_x(I,j,k,m)
              enddo
            enddo
          endif
        enddo
      else
        J = segment%HI%JsdB
        do m=1,OBC%ntr
          if (allocated(segment%tr_Reg%Tr(m)%tres)) then
            do k=1,GV%ke
              do i=segment%HI%isd,segment%HI%ied
                segment%tr_Reg%Tr(m)%tres(i,J,k) = segment%tr_Reg%Tr(m)%scale * OBC%tres_y(i,J,k,m)
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
    if (segment%oblique) call gradient_at_q_points(G, GV, segment, u_new(:,:,:), v_new(:,:,:))
    if (segment%direction == OBC_DIRECTION_E) then
      I=segment%HI%IsdB
      if (I<G%HI%IscB) cycle
      do k=1,nz ;  do j=segment%HI%jsd,segment%HI%jed
        if (segment%radiation) then
          dhdt = (u_old(I-1,j,k) - u_new(I-1,j,k)) !old-new
          dhdx = (u_new(I-1,j,k) - u_new(I-2,j,k)) !in new time backward sashay for I-1
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
          dhdx = (u_new(I-1,j,k) - u_new(I-2,j,k)) !in new time backward sashay for I-1
          if (dhdt*(segment%grad_normal(J,1,k) + segment%grad_normal(J-1,1,k)) > 0.0) then
            dhdy = segment%grad_normal(J-1,1,k)
          elseif (dhdt*(segment%grad_normal(J,1,k) + segment%grad_normal(J-1,1,k)) == 0.0) then
            dhdy = 0.0
          else
            dhdy = segment%grad_normal(J,1,k)
          endif
          if (dhdt*dhdx < 0.0) dhdt = 0.0
          cff_new = max((dhdx*dhdx) + (dhdy*dhdy), eps)
          rx_new = min(dhdt*dhdx, cff_new*rx_max)
          ry_new = min(cff_new,max(dhdt*dhdy,-cff_new))
          if (gamma_u < 1.0) then
            rx_avg = (1.0-gamma_u)*segment%rx_norm_obl(I,j,k) + gamma_u*rx_new
            ry_avg = (1.0-gamma_u)*segment%ry_norm_obl(I,j,k) + gamma_u*ry_new
            cff_avg = (1.0-gamma_u)*segment%cff_normal(I,j,k) + gamma_u*cff_new
          else
            rx_avg = rx_new
            ry_avg = ry_new
            cff_avg = cff_new
          endif
          segment%rx_norm_obl(I,j,k) = rx_avg
          segment%ry_norm_obl(I,j,k) = ry_avg
          segment%cff_normal(I,j,k) = cff_avg
          segment%normal_vel(I,j,k) = ((cff_avg*u_new(I,j,k) + rx_avg*u_new(I-1,j,k)) - &
                             (max(ry_avg,0.0)*segment%grad_normal(J-1,2,k) + &
                              min(ry_avg,0.0)*segment%grad_normal(J,2,k))) / &
                           (cff_avg + rx_avg)
          if (gamma_u < 1.0) then
            ! Copy restart fields into 3-d arrays. This is an inefficient and temporary
            ! implementation as a work-around to limitations in restart capability
            OBC%rx_oblique_u(I,j,k) = segment%rx_norm_obl(I,j,k)
            OBC%ry_oblique_u(I,j,k) = segment%ry_norm_obl(I,j,k)
            OBC%cff_normal_u(I,j,k) = segment%cff_normal(I,j,k)
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
              dhdx = v_new(i,J,k)-v_new(i-1,J,k) !in new time backward sashay for I-1
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
!           if (G%mask2dCu(I-1,j) > 0.0 .and. G%mask2dCu(I-1,j+1) > 0.0) then
!             rx_avg = 0.5*(u_new(I-1,j,k) + u_new(I-1,j+1,k)) * dt * G%IdxBu(I-1,J)
!           elseif (G%mask2dCu(I-1,j) > 0.0) then
!             rx_avg = u_new(I-1,j,k) * dt * G%IdxBu(I-1,J)
!           elseif (G%mask2dCu(I-1,j+1) > 0.0) then
!             rx_avg = u_new(I-1,j+1,k) * dt * G%IdxBu(I-1,J)
!           else
!             rx_avg = 0.0
!           endif
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
              dhdx = v_new(i,J,k)-v_new(i-1,J,k) !in new time backward sashay for I-1
              if (dhdt*(segment%grad_tan(j,1,k) + segment%grad_tan(j+1,1,k)) > 0.0) then
                dhdy = segment%grad_tan(j,1,k)
              elseif (dhdt*(segment%grad_tan(j,1,k) + segment%grad_tan(j+1,1,k)) == 0.0) then
                dhdy = 0.0
              else
                dhdy = segment%grad_tan(j+1,1,k)
              endif
              if (dhdt*dhdx < 0.0) dhdt = 0.0
              cff_new = max((dhdx*dhdx) + (dhdy*dhdy), eps)
              rx_new = min(dhdt*dhdx, cff_new*rx_max)
              ry_new = min(cff_new,max(dhdt*dhdy,-cff_new))
              rx_tang_obl(I,J,k) = rx_new
              ry_tang_obl(I,J,k) = ry_new
              cff_tangential(I,J,k) = cff_new
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
          dhdx = (u_new(I+1,j,k) - u_new(I+2,j,k)) !in new time forward sashay for I+1
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
          dhdx = (u_new(I+1,j,k) - u_new(I+2,j,k)) !in new time forward sashay for I+1
          if (dhdt*(segment%grad_normal(J,1,k) + segment%grad_normal(J-1,1,k)) > 0.0) then
            dhdy = segment%grad_normal(J-1,1,k)
          elseif (dhdt*(segment%grad_normal(J,1,k) + segment%grad_normal(J-1,1,k)) == 0.0) then
            dhdy = 0.0
          else
            dhdy = segment%grad_normal(J,1,k)
          endif
          if (dhdt*dhdx < 0.0) dhdt = 0.0

          cff_new = max((dhdx*dhdx) + (dhdy*dhdy), eps)
          rx_new = min(dhdt*dhdx, cff_new*rx_max)
          ry_new = min(cff_new,max(dhdt*dhdy,-cff_new))
          if (gamma_u < 1.0) then
            rx_avg = (1.0-gamma_u)*segment%rx_norm_obl(I,j,k) + gamma_u*rx_new
            ry_avg = (1.0-gamma_u)*segment%ry_norm_obl(I,j,k) + gamma_u*ry_new
            cff_avg = (1.0-gamma_u)*segment%cff_normal(I,j,k) + gamma_u*cff_new
          else
            rx_avg = rx_new
            ry_avg = ry_new
            cff_avg = cff_new
          endif
          segment%rx_norm_obl(I,j,k) = rx_avg
          segment%ry_norm_obl(I,j,k) = ry_avg
          segment%cff_normal(I,j,k) = cff_avg
          segment%normal_vel(I,j,k) = ((cff_avg*u_new(I,j,k) + rx_avg*u_new(I+1,j,k)) - &
                                       (max(ry_avg,0.0)*segment%grad_normal(J-1,2,k) + &
                                        min(ry_avg,0.0)*segment%grad_normal(J,2,k))) / &
                                      (cff_avg + rx_avg)
          if (gamma_u < 1.0) then
            ! Copy restart fields into 3-d arrays. This is an inefficient and temporary issues
            ! implemented as a work-around to limitations in restart capability
            OBC%rx_oblique_u(I,j,k) = segment%rx_norm_obl(I,j,k)
            OBC%ry_oblique_u(I,j,k) = segment%ry_norm_obl(I,j,k)
            OBC%cff_normal_u(I,j,k) = segment%cff_normal(I,j,k)
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
              dhdx = v_new(i+1,J,k)-v_new(i+2,J,k) !in new time backward sashay for I-1
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
!           if (G%mask2dCu(I+1,j) > 0.0 .and. G%mask2dCu(I+1,j+1) > 0.0) then
!             rx_avg = 0.5*(u_new(I+1,j,k) + u_new(I+1,j+1,k)) * dt * G%IdxBu(I+1,J)
!           elseif (G%mask2dCu(I+1,j) > 0.0) then
!             rx_avg = u_new(I+1,j,k) * dt * G%IdxBu(I+1,J)
!           elseif (G%mask2dCu(I+1,j+1) > 0.0) then
!             rx_avg = u_new(I+1,j+1,k) * dt * G%IdxBu(I+1,J)
!           else
!             rx_avg = 0.0
!           endif
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
              dhdx = v_new(i+1,J,k)-v_new(i+2,J,k) !in new time backward sashay for I-1
              if (dhdt*(segment%grad_tan(j,1,k) + segment%grad_tan(j+1,1,k)) > 0.0) then
                dhdy = segment%grad_tan(j,1,k)
              elseif (dhdt*(segment%grad_tan(j,1,k) + segment%grad_tan(j+1,1,k)) == 0.0) then
                dhdy = 0.0
              else
                dhdy = segment%grad_tan(j+1,1,k)
              endif
              if (dhdt*dhdx < 0.0) dhdt = 0.0
              cff_new = max((dhdx*dhdx) + (dhdy*dhdy), eps)
              rx_new = min(dhdt*dhdx, cff_new*rx_max)
              ry_new = min(cff_new,max(dhdt*dhdy,-cff_new))
              rx_tang_obl(I,J,k) = rx_new
              ry_tang_obl(I,J,k) = ry_new
              cff_tangential(I,J,k) = cff_new
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
          dhdy = (v_new(i,J-1,k) - v_new(i,J-2,k)) !in new time backward sashay for J-1
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
          dhdy = (v_new(i,J-1,k) - v_new(i,J-2,k)) !in new time backward sashay for J-1
          if (dhdt*(segment%grad_normal(I,1,k) + segment%grad_normal(I-1,1,k)) > 0.0) then
            dhdx = segment%grad_normal(I-1,1,k)
          elseif (dhdt*(segment%grad_normal(I,1,k) + segment%grad_normal(I-1,1,k)) == 0.0) then
            dhdx = 0.0
          else
            dhdx = segment%grad_normal(I,1,k)
          endif
          if (dhdt*dhdy < 0.0) dhdt = 0.0
          cff_new = max((dhdx*dhdx) + (dhdy*dhdy), eps)
          ry_new = min(dhdt*dhdy, cff_new*ry_max)
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
          segment%rx_norm_obl(i,J,k) = rx_avg
          segment%ry_norm_obl(i,J,k) = ry_avg
          segment%cff_normal(i,J,k) = cff_avg
          segment%normal_vel(i,J,k) = ((cff_avg*v_new(i,J,k) + ry_avg*v_new(i,J-1,k)) - &
                                       (max(rx_avg,0.0)*segment%grad_normal(I-1,2,k) +&
                                        min(rx_avg,0.0)*segment%grad_normal(I,2,k))) / &
                                      (cff_avg + ry_avg)
          if (gamma_u < 1.0) then
            ! Copy restart fields into 3-d arrays. This is an inefficient and temporary issues
            ! implemented as a work-around to limitations in restart capability
            OBC%rx_oblique_v(i,J,k) = segment%rx_norm_obl(i,J,k)
            OBC%ry_oblique_v(i,J,k) = segment%ry_norm_obl(i,J,k)
            OBC%cff_normal_v(i,J,k) = segment%cff_normal(i,J,k)
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
              dhdy = u_new(I,j-1,k)-u_new(I,j-2,k) !in new time backward sashay for I-1
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
!           if (G%mask2dCv(i,J-1) > 0.0 .and. G%mask2dCv(i+1,J-1) > 0.0) then
!             ry_avg = 0.5*(v_new(i,J-1,k) + v_new(i+1,J-1,k) * dt * G%IdyBu(I,J-1))
!           elseif (G%mask2dCv(i,J-1) > 0.0) then
!             ry_avg = v_new(i,J-1,k) * dt *G%IdyBu(I,J-1)
!           elseif (G%mask2dCv(i+1,J-1) > 0.0) then
!             ry_avg = v_new(i+1,J-1,k) * dt *G%IdyBu(I,J-1)
!           else
!             ry_avg = 0.0
!           endif
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
              dhdy = u_new(I,j,k)-u_new(I,j-1,k) !in new time backward sashay for I-1
              if (dhdt*(segment%grad_tan(i,1,k) + segment%grad_tan(i+1,1,k)) > 0.0) then
                dhdx = segment%grad_tan(i,1,k)
              elseif (dhdt*(segment%grad_tan(i,1,k) + segment%grad_tan(i+1,1,k)) == 0.0) then
                dhdx = 0.0
              else
                dhdx = segment%grad_tan(i+1,1,k)
              endif
              if (dhdt*dhdy < 0.0) dhdt = 0.0
              cff_new = max((dhdx*dhdx) + (dhdy*dhdy), eps)
              ry_new = min(dhdt*dhdy, cff_new*ry_max)
              rx_new = min(cff_new,max(dhdt*dhdx,-cff_new))
              rx_tang_obl(I,J,k) = rx_new
              ry_tang_obl(I,J,k) = ry_new
              cff_tangential(I,J,k) = cff_new
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
          dhdy = (v_new(i,J+1,k) - v_new(i,J+2,k)) !in new time backward sashay for J-1
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
          dhdy = (v_new(i,J+1,k) - v_new(i,J+2,k)) !in new time backward sashay for J-1
          if (dhdt*(segment%grad_normal(I,1,k) + segment%grad_normal(I-1,1,k)) > 0.0) then
            dhdx = segment%grad_normal(I-1,1,k)
          elseif (dhdt*(segment%grad_normal(I,1,k) + segment%grad_normal(I-1,1,k)) == 0.0) then
            dhdx = 0.0
          else
            dhdx = segment%grad_normal(I,1,k)
          endif
          if (dhdt*dhdy < 0.0) dhdt = 0.0

          cff_new = max((dhdx*dhdx) + (dhdy*dhdy), eps)
          ry_new = min(dhdt*dhdy, cff_new*ry_max)
          rx_new = min(cff_new,max(dhdt*dhdx,-cff_new))
          if (gamma_u < 1.0) then
            rx_avg = (1.0-gamma_u)*segment%rx_norm_obl(i,J,k) + gamma_u*rx_new
            ry_avg = (1.0-gamma_u)*segment%ry_norm_obl(i,J,k) + gamma_u*ry_new
            cff_avg = (1.0-gamma_u)*segment%cff_normal(i,J,k) + gamma_u*cff_new
          else
            rx_avg = rx_new
            ry_avg = ry_new
            cff_avg = cff_new
          endif
          segment%rx_norm_obl(i,J,k) = rx_avg
          segment%ry_norm_obl(i,J,k) = ry_avg
          segment%cff_normal(i,J,k) = cff_avg
          segment%normal_vel(i,J,k) = ((cff_avg*v_new(i,J,k) + ry_avg*v_new(i,J+1,k)) - &
                                       (max(rx_avg,0.0)*segment%grad_normal(I-1,2,k) + &
                                        min(rx_avg,0.0)*segment%grad_normal(I,2,k))) / &
                                      (cff_avg + ry_avg)
          if (gamma_u < 1.0) then
            ! Copy restart fields into 3-d arrays. This is an inefficient and temporary issues
            ! implemented as a work-around to limitations in restart capability
            OBC%rx_oblique_v(i,J,k) = segment%rx_norm_obl(i,J,k)
            OBC%ry_oblique_v(i,J,k) = segment%ry_norm_obl(i,J,k)
            OBC%cff_normal_v(i,J,k) = segment%cff_normal(i,J,k)
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
              dhdy = u_new(I,j+1,k)-u_new(I,j+2,k) !in new time backward sashay for I-1
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
!           if (G%mask2dCv(i,J+1) > 0.0 .and. G%mask2dCv(i+1,J+1) > 0.0) then
!             ry_avg = 0.5*(v_new(i,J+1,k) + v_new(i+1,J+1,k)) * dt * G%IdyBu(I,J+1)
!           elseif (G%mask2dCv(i,J+1) > 0.0) then
!             ry_avg = v_new(i,J+1,k) * dt * G%IdyBu(I,J+1)
!           elseif (G%mask2dCv(i+1,J+1) > 0.0) then
!             ry_avg = v_new(i+1,J+1,k) * dt * G%IdyBu(I,J+1)
!           else
!             ry_avg = 0.0
!           endif
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
              dhdy = u_new(I,j+1,k)-u_new(I,j+2,k) !in new time backward sashay for I-1
              if (dhdt*(segment%grad_tan(i,1,k) + segment%grad_tan(i+1,1,k)) > 0.0) then
                dhdx = segment%grad_tan(i,1,k)
              elseif (dhdt*(segment%grad_tan(i,1,k) + segment%grad_tan(i+1,1,k)) == 0.0) then
                dhdx = 0.0
              else
                dhdx = segment%grad_tan(i+1,1,k)
              endif
              if (dhdt*dhdy < 0.0) dhdt = 0.0
              cff_new = max((dhdx*dhdx) + (dhdy*dhdy), eps)
              ry_new = min(dhdt*dhdy, cff_new*ry_max)
              rx_new = min(cff_new,max(dhdt*dhdx,-cff_new))
              rx_tang_obl(I,J,k) = rx_new
              ry_tang_obl(I,J,k) = ry_new
              cff_tangential(I,J,k) = cff_new
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
  call open_boundary_apply_normal_flow(OBC, G, GV, u_new, v_new)

  call pass_vector(u_new, v_new, G%Domain, clock=id_clock_pass)

  if (OBC%debug) then
    sym = G%Domain%symmetric
    if (OBC%radiation_BCs_exist_globally) then
      call uvchksum("radiation_OBCs: OBC%r[xy]_normal", OBC%rx_normal, OBC%ry_normal, G%HI, &
                  haloshift=0, symmetric=sym, scalar_pair=.true., unscale=1.0)
    endif
    if (OBC%oblique_BCs_exist_globally) then
      call uvchksum("radiation_OBCs: OBC%r[xy]_oblique_[uv]", OBC%rx_oblique_u, OBC%ry_oblique_v, G%HI, &
                  haloshift=0, symmetric=sym, scalar_pair=.true., unscale=1.0/US%L_T_to_m_s**2)
      call uvchksum("radiation_OBCs: OBC%r[yx]_oblique_[uv]", OBC%ry_oblique_u, OBC%rx_oblique_v, G%HI, &
                  haloshift=0, symmetric=sym, scalar_pair=.true., unscale=1.0/US%L_T_to_m_s**2)
      call uvchksum("radiation_OBCs: OBC%cff_normal_[uv]", OBC%cff_normal_u, OBC%cff_normal_v, G%HI, &
                  haloshift=0, symmetric=sym, scalar_pair=.true., unscale=1.0/US%L_T_to_m_s**2)
    endif
    if ((OBC%ntr > 0) .and. allocated(OBC%tres_x) .and. allocated(OBC%tres_y)) then
      do m=1,OBC%ntr
        write(var_num,'(I3.3)') m
        call uvchksum("radiation_OBCs: OBC%tres_[xy]_"//var_num, OBC%tres_x(:,:,:,m), OBC%tres_y(:,:,:,m), G%HI, &
                      haloshift=0, symmetric=sym, scalar_pair=.true., unscale=1.0)
      enddo
    endif
  endif

end subroutine radiation_open_bdry_conds

!> Applies OBC values stored in segments to 3d u,v fields
subroutine open_boundary_apply_normal_flow(OBC, G, GV, u, v)
  ! Arguments
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary control structure
  type(ocean_grid_type),                     intent(inout) :: G   !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV  !< The ocean's vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u   !< u field to update on open
                                                                  !! boundaries [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v   !< v field to update on open
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
        do k=1,GV%ke ;  do j=segment%HI%jsd,segment%HI%jed
          u(I,j,k) = segment%normal_vel(I,j,k)
        enddo ; enddo
      elseif (segment%is_N_or_S) then
        J=segment%HI%JsdB
        do k=1,GV%ke ;  do i=segment%HI%isd,segment%HI%ied
          v(i,J,k) = segment%normal_vel(i,J,k)
        enddo ; enddo
      endif
    endif
  enddo

end subroutine open_boundary_apply_normal_flow

!> Applies zero values to 3d u,v fields on OBC segments
subroutine open_boundary_zero_normal_flow(OBC, G, GV, u, v)
  ! Arguments
  type(ocean_OBC_type),                       pointer       :: OBC !< Open boundary control structure
  type(ocean_grid_type),                      intent(inout) :: G   !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV  !< The ocean's vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u   !< u field to update on open boundaries [arbitrary]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v   !< v field to update on open boundaries [arbitrary]
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
      do k=1,GV%ke ;  do j=segment%HI%jsd,segment%HI%jed
        u(I,j,k) = 0.
      enddo ; enddo
    elseif (segment%is_N_or_S) then
      J=segment%HI%JsdB
      do k=1,GV%ke ;  do i=segment%HI%isd,segment%HI%ied
        v(i,J,k) = 0.
      enddo ; enddo
    endif
  enddo

end subroutine open_boundary_zero_normal_flow

!> Calculate the tangential gradient of the normal flow at the boundary q-points.
subroutine gradient_at_q_points(G, GV, segment, uvel, vvel)
  type(ocean_grid_type),   intent(in) :: G !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV  !< The ocean's vertical grid structure
  type(OBC_segment_type), intent(inout) :: segment !< OBC segment structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)    :: uvel !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)    :: vvel !< meridional velocity [L T-1 ~> m s-1]
  integer :: i,j,k

  if (.not. segment%on_pe) return

  if (segment%is_E_or_W) then
    if (segment%direction == OBC_DIRECTION_E) then
      I=segment%HI%isdB
      do k=1,GV%ke
        do J=max(segment%HI%JsdB, G%HI%JsdB+1),min(segment%HI%JedB, G%HI%JedB-1)
          segment%grad_normal(J,1,k) = (uvel(I-1,j+1,k)-uvel(I-1,j,k)) * G%mask2dBu(I-1,J)
          segment%grad_normal(J,2,k) = (uvel(I,j+1,k)-uvel(I,j,k)) * G%mask2dBu(I,J)
        enddo
      enddo
      if (segment%oblique_tan) then
        do k=1,GV%ke
          do J=max(segment%HI%jsd-1, G%HI%jsd),min(segment%HI%jed+1, G%HI%jed)
            segment%grad_tan(j,1,k) = (vvel(i-1,J,k)-vvel(i-1,J-1,k)) * G%mask2dT(i-1,j)
            segment%grad_tan(j,2,k) = (vvel(i,J,k)-vvel(i,J-1,k)) * G%mask2dT(i,j)
          enddo
        enddo
      endif
      if (segment%oblique_grad) then
        do k=1,GV%ke
          do J=max(segment%HI%jsd, G%HI%jsd+1),min(segment%HI%jed, G%HI%jed-1)
            segment%grad_gradient(j,1,k) = (((vvel(i-1,J,k) - vvel(i-2,J,k))*G%IdxBu(I-2,J)) - &
                 ((vvel(i-1,J-1,k) - vvel(i-2,J-1,k))*G%IdxBu(I-2,J-1))) * G%mask2dCu(I-2,j)
            segment%grad_gradient(j,2,k) = (((vvel(i,J,k) - vvel(i-1,J,k))*G%IdxBu(I-1,J)) - &
                 ((vvel(i,J-1,k) - vvel(i-1,J-1,k))*G%IdxBu(I-1,J-1))) * G%mask2dCu(I-1,j)
          enddo
        enddo
      endif
    else ! western segment
      I=segment%HI%isdB
      do k=1,GV%ke
        do J=max(segment%HI%JsdB, G%HI%JsdB+1),min(segment%HI%JedB, G%HI%JedB-1)
          segment%grad_normal(J,1,k) = (uvel(I+1,j+1,k)-uvel(I+1,j,k)) * G%mask2dBu(I+1,J)
          segment%grad_normal(J,2,k) = (uvel(I,j+1,k)-uvel(I,j,k)) * G%mask2dBu(I,J)
        enddo
      enddo
      if (segment%oblique_tan) then
        do k=1,GV%ke
          do J=max(segment%HI%jsd-1, G%HI%jsd),min(segment%HI%jed+1, G%HI%jed)
            segment%grad_tan(j,1,k) = (vvel(i+2,J,k)-vvel(i+2,J-1,k)) * G%mask2dT(i+2,j)
            segment%grad_tan(j,2,k) = (vvel(i+1,J,k)-vvel(i+1,J-1,k)) * G%mask2dT(i+1,j)
          enddo
        enddo
      endif
      if (segment%oblique_grad) then
        do k=1,GV%ke
          do J=max(segment%HI%jsd, G%HI%jsd+1),min(segment%HI%jed, G%HI%jed-1)
            segment%grad_gradient(j,1,k) = (((vvel(i+3,J,k) - vvel(i+2,J,k))*G%IdxBu(I+2,J)) - &
                 ((vvel(i+3,J-1,k) - vvel(i+2,J-1,k))*G%IdxBu(I+2,J-1))) * G%mask2dCu(I+2,j)
            segment%grad_gradient(j,2,k) = (((vvel(i+2,J,k) - vvel(i+1,J,k))*G%IdxBu(I+1,J)) - &
                 ((vvel(i+2,J-1,k) - vvel(i+1,J-1,k))*G%IdxBu(I+1,J-1))) * G%mask2dCu(I+1,j)
          enddo
        enddo
      endif
    endif
  elseif (segment%is_N_or_S) then
    if (segment%direction == OBC_DIRECTION_N) then
      J=segment%HI%jsdB
      do k=1,GV%ke
        do I=max(segment%HI%IsdB, G%HI%IsdB+1),min(segment%HI%IedB, G%HI%IedB-1)
          segment%grad_normal(I,1,k) = (vvel(i+1,J-1,k)-vvel(i,J-1,k)) * G%mask2dBu(I,J-1)
          segment%grad_normal(I,2,k) = (vvel(i+1,J,k)-vvel(i,J,k)) * G%mask2dBu(I,J)
        enddo
      enddo
      if (segment%oblique_tan) then
        do k=1,GV%ke
          do I=max(segment%HI%isd-1, G%HI%isd),min(segment%HI%ied+1, G%HI%ied)
            segment%grad_tan(i,1,k) = (uvel(I,j-1,k)-uvel(I-1,j-1,k)) * G%mask2dT(i,j-1)
            segment%grad_tan(i,2,k) = (uvel(I,j,k)-uvel(I-1,j,k)) * G%mask2dT(i,j)
          enddo
        enddo
      endif
      if (segment%oblique_grad) then
        do k=1,GV%ke
          do I=max(segment%HI%isd, G%HI%isd+1),min(segment%HI%ied, G%HI%ied-1)
            segment%grad_gradient(i,1,k) = (((uvel(I,j-1,k) - uvel(I,j-2,k))*G%IdyBu(I,J-2)) - &
                 ((uvel(I-1,j-1,k) - uvel(I-1,j-2,k))*G%IdyBu(I-1,J-2))) * G%mask2dCv(i,J-2)
            segment%grad_gradient(i,2,k) = (((uvel(I,j,k) - uvel(I,j-1,k))*G%IdyBu(I,J-1)) - &
                 ((uvel(I-1,j,k) - uvel(I-1,j-1,k))*G%IdyBu(I-1,J-1))) * G%mask2dCv(i,J-1)
          enddo
        enddo
      endif
    else ! south segment
      J=segment%HI%jsdB
      do k=1,GV%ke
        do I=max(segment%HI%IsdB, G%HI%IsdB+1),min(segment%HI%IedB, G%HI%IedB-1)
          segment%grad_normal(I,1,k) = (vvel(i+1,J+1,k)-vvel(i,J+1,k)) * G%mask2dBu(I,J+1)
          segment%grad_normal(I,2,k) = (vvel(i+1,J,k)-vvel(i,J,k)) * G%mask2dBu(I,J)
        enddo
      enddo
      if (segment%oblique_tan) then
        do k=1,GV%ke
          do I=max(segment%HI%isd-1, G%HI%isd),min(segment%HI%ied+1, G%HI%ied)
            segment%grad_tan(i,1,k) = (uvel(I,j+2,k)-uvel(I-1,j+2,k)) * G%mask2dT(i,j+2)
            segment%grad_tan(i,2,k) = (uvel(I,j+1,k)-uvel(I-1,j+1,k)) * G%mask2dT(i,j+1)
          enddo
        enddo
      endif
      if (segment%oblique_grad) then
        do k=1,GV%ke
          do I=max(segment%HI%isd, G%HI%isd+1),min(segment%HI%ied, G%HI%ied-1)
            segment%grad_gradient(i,1,k) = (((uvel(I,j+3,k) - uvel(I,j+2,k))*G%IdyBu(I,J+2)) - &
                 ((uvel(I-1,j+3,k) - uvel(I-1,j+2,k))*G%IdyBu(I-1,J+2))) * G%mask2dCv(i,J+2)
            segment%grad_gradient(i,2,k) = (((uvel(I,j+2,k) - uvel(I,j+1,k))*G%IdyBu(I,J+1)) - &
                 ((uvel(I-1,j+2,k) - uvel(I-1,j+1,k))*G%IdyBu(I-1,J+1))) * G%mask2dCv(i,J+1)
          enddo
        enddo
      endif
    endif
  endif

end subroutine gradient_at_q_points


!> Return the field number on the segment for the named field, or -1 if there is no field with that name.
function lookup_seg_field(OBC_seg, field)
  type(OBC_segment_type), intent(in) :: OBC_seg !< OBC segment
  character(len=32), intent(in) :: field !< The field name
  integer :: lookup_seg_field
  ! Local variables
  integer :: n

  lookup_seg_field = -1
  do n=1,OBC_seg%num_fields
    if (trim(field) == OBC_seg%field(n)%name) then
      lookup_seg_field = n
      return
    endif
  enddo

end function lookup_seg_field

!> Return the tracer index from its name
function get_tracer_index(OBC_seg,tr_name)
  type(OBC_segment_type), pointer :: OBC_seg !< OBC segment
  character(len=*), intent(in) :: tr_name   !< The field name
  integer :: get_tracer_index, it
  get_tracer_index = -1
  it = 1
  do while(allocated(OBC_seg%tr_Reg%Tr(it)%t))
    if (trim(OBC_seg%tr_Reg%Tr(it)%name) == trim(tr_name)) then
      get_tracer_index = it
      exit
    endif
    it = it + 1
  enddo
end function get_tracer_index

!> Allocate segment data fields
subroutine allocate_OBC_segment_data(OBC, segment)
  type(ocean_OBC_type),   intent(in)    :: OBC     !< Open boundary structure
  type(OBC_segment_type), intent(inout) :: segment !< Open boundary segment
  ! Local variables
  integer :: isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB
  integer :: IscB, IecB, JscB, JecB

  isd = segment%HI%isd ; ied = segment%HI%ied
  jsd = segment%HI%jsd ; jed = segment%HI%jed
  IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
  JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB
  IscB = segment%HI%IscB ; IecB = segment%HI%IecB
  JscB = segment%HI%JscB ; JecB = segment%HI%JecB

  if (.not. segment%on_pe) return

  if (segment%is_E_or_W) then
    ! If these are just Flather, change update_OBC_segment_data accordingly
    allocate(segment%Cg(IsdB:IedB,jsd:jed), source=0.0)
    allocate(segment%Htot(IsdB:IedB,jsd:jed), source=0.0)
    ! Allocate dZtot with extra values at the end to avoid segmentation faults in cases where
    ! it is interpolated to OBC vorticity points.
    allocate(segment%dZtot(IsdB:IedB,jsd-1:jed+1), source=0.0)
    allocate(segment%h(IsdB:IedB,jsd:jed,OBC%ke), source=0.0)
    allocate(segment%SSH(IsdB:IedB,jsd:jed), source=0.0)
    if (segment%radiation) &
      allocate(segment%rx_norm_rad(IsdB:IedB,jsd:jed,OBC%ke), source=0.0)
    allocate(segment%normal_vel(IsdB:IedB,jsd:jed,OBC%ke), source=0.0)
    allocate(segment%normal_vel_bt(IsdB:IedB,jsd:jed), source=0.0)
    allocate(segment%normal_trans(IsdB:IedB,jsd:jed,OBC%ke), source=0.0)
    if (segment%nudged) &
      allocate(segment%nudged_normal_vel(IsdB:IedB,jsd:jed,OBC%ke), source=0.0)
    if (segment%radiation_tan .or. segment%nudged_tan .or. segment%specified_tan .or. &
        segment%oblique_tan .or. OBC%computed_vorticity .or. OBC%computed_strain) &
      allocate(segment%tangential_vel(IsdB:IedB,JsdB:JedB,OBC%ke), source=0.0)
    if (segment%nudged_tan) &
      allocate(segment%nudged_tangential_vel(IsdB:IedB,JsdB:JedB,OBC%ke), source=0.0)
    if (segment%nudged_grad) &
      allocate(segment%nudged_tangential_grad(IsdB:IedB,JsdB:JedB,OBC%ke), source=0.0)
    if (OBC%specified_vorticity .or. OBC%specified_strain .or. segment%radiation_grad .or. &
              segment%oblique_grad .or. segment%specified_grad) &
      allocate(segment%tangential_grad(IsdB:IedB,JsdB:JedB,OBC%ke), source=0.0)
    if (segment%oblique) then
      allocate(segment%grad_normal(JsdB:JedB,2,OBC%ke), source=0.0)
      allocate(segment%rx_norm_obl(IsdB:IedB,jsd:jed,OBC%ke), source=0.0)
      allocate(segment%ry_norm_obl(IsdB:IedB,jsd:jed,OBC%ke), source=0.0)
      allocate(segment%cff_normal(IsdB:IedB,jsd:jed,OBC%ke), source=0.0)
    endif
    if (segment%oblique_tan) &
      allocate(segment%grad_tan(jsd-1:jed+1,2,OBC%ke), source=0.0)
    if (segment%oblique_grad) &
      allocate(segment%grad_gradient(jsd:jed,2,OBC%ke), source=0.0)
  endif

  if (segment%is_N_or_S) then
    ! If these are just Flather, change update_OBC_segment_data accordingly
    allocate(segment%Cg(isd:ied,JsdB:JedB), source=0.0)
    allocate(segment%Htot(isd:ied,JsdB:JedB), source=0.0)
    ! Allocate dZtot with extra values at the end to avoid segmentation faults in cases where
    ! it is interpolated to OBC vorticity points.
    allocate(segment%dZtot(isd-1:ied+1,JsdB:JedB), source=0.0)
    allocate(segment%h(isd:ied,JsdB:JedB,OBC%ke), source=0.0)
    allocate(segment%SSH(isd:ied,JsdB:JedB), source=0.0)
    if (segment%radiation) &
      allocate(segment%ry_norm_rad(isd:ied,JsdB:JedB,OBC%ke), source=0.0)
    allocate(segment%normal_vel(isd:ied,JsdB:JedB,OBC%ke), source=0.0)
    allocate(segment%normal_vel_bt(isd:ied,JsdB:JedB), source=0.0)
    allocate(segment%normal_trans(isd:ied,JsdB:JedB,OBC%ke), source=0.0)
    if (segment%nudged) &
      allocate(segment%nudged_normal_vel(isd:ied,JsdB:JedB,OBC%ke), source=0.0)
    if (segment%radiation_tan .or. segment%nudged_tan .or. segment%specified_tan .or. &
        segment%oblique_tan .or. OBC%computed_vorticity .or. OBC%computed_strain) &
      allocate(segment%tangential_vel(IsdB:IedB,JsdB:JedB,OBC%ke), source=0.0)
    if (segment%nudged_tan) &
      allocate(segment%nudged_tangential_vel(IsdB:IedB,JsdB:JedB,OBC%ke), source=0.0)
    if (segment%nudged_grad) &
      allocate(segment%nudged_tangential_grad(IsdB:IedB,JsdB:JedB,OBC%ke), source=0.0)
    if (OBC%specified_vorticity .or. OBC%specified_strain .or. segment%radiation_grad .or. &
              segment%oblique_grad .or. segment%specified_grad) &
      allocate(segment%tangential_grad(IsdB:IedB,JsdB:JedB,OBC%ke), source=0.0)
    if (segment%oblique) then
      allocate(segment%grad_normal(IsdB:IedB,2,OBC%ke), source=0.0)
      allocate(segment%rx_norm_obl(isd:ied,JsdB:JedB,OBC%ke), source=0.0)
      allocate(segment%ry_norm_obl(isd:ied,JsdB:JedB,OBC%ke), source=0.0)
      allocate(segment%cff_normal(isd:ied,JsdB:JedB,OBC%ke), source=0.0)
    endif
    if (segment%oblique_tan) &
      allocate(segment%grad_tan(isd-1:ied+1,2,OBC%ke), source=0.0)
    if (segment%oblique_grad) &
      allocate(segment%grad_gradient(isd:ied,2,OBC%ke), source=0.0)
  endif

end subroutine allocate_OBC_segment_data

!> Deallocate segment data fields
subroutine deallocate_OBC_segment_data(segment)
  type(OBC_segment_type), intent(inout) :: segment !< Open boundary segment

  if (.not. segment%on_pe) return

  if (allocated(segment%Cg)) deallocate(segment%Cg)
  if (allocated(segment%Htot)) deallocate(segment%Htot)
  if (allocated(segment%dZtot)) deallocate(segment%dZtot)
  if (allocated(segment%h)) deallocate(segment%h)
  if (allocated(segment%SSH)) deallocate(segment%SSH)
  if (allocated(segment%rx_norm_rad)) deallocate(segment%rx_norm_rad)
  if (allocated(segment%ry_norm_rad)) deallocate(segment%ry_norm_rad)
  if (allocated(segment%rx_norm_obl)) deallocate(segment%rx_norm_obl)
  if (allocated(segment%ry_norm_obl)) deallocate(segment%ry_norm_obl)
  if (allocated(segment%cff_normal)) deallocate(segment%cff_normal)
  if (allocated(segment%grad_normal)) deallocate(segment%grad_normal)
  if (allocated(segment%grad_tan)) deallocate(segment%grad_tan)
  if (allocated(segment%grad_gradient)) deallocate(segment%grad_gradient)
  if (allocated(segment%normal_vel)) deallocate(segment%normal_vel)
  if (allocated(segment%normal_vel_bt)) deallocate(segment%normal_vel_bt)
  if (allocated(segment%normal_trans)) deallocate(segment%normal_trans)
  if (allocated(segment%nudged_normal_vel)) deallocate(segment%nudged_normal_vel)
  if (allocated(segment%tangential_vel)) deallocate(segment%tangential_vel)
  if (allocated(segment%nudged_tangential_vel)) deallocate(segment%nudged_tangential_vel)
  if (allocated(segment%nudged_tangential_grad)) deallocate(segment%nudged_tangential_grad)
  if (allocated(segment%tangential_grad)) deallocate(segment%tangential_grad)

  if (associated(segment%tr_Reg)) call segment_tracer_registry_end(segment%tr_Reg)


end subroutine deallocate_OBC_segment_data

!> Set tangential velocities outside of open boundaries to silly values
!! (used for checking the interior state is independent of values outside
!! of the domain).
subroutine open_boundary_test_extern_uv(G, GV, OBC, u, v)
  type(ocean_grid_type),                     intent(in)    :: G !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV  !< The ocean's vertical grid structure
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),intent(inout) :: u !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),intent(inout) :: v !< Meridional velocity [L T-1 ~> m s-1]
  ! Local variables
  integer :: i, j, k, n

  if (.not. associated(OBC)) return

  do n = 1, OBC%number_of_segments
    do k = 1, GV%ke
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

  silly_h = GV%Z_to_H * OBC%silly_h  ! This rescaling is here because GV was initialized after OBC.

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
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h    !< Thickness [H ~> m or kg m-2]
  type(time_type),                           intent(in)    :: Time !< Model time

  ! Local variables
  integer :: c, i, j, k, is, ie, js, je, isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB, n, m, nz, nt, nk_dst
  type(OBC_segment_type), pointer :: segment => NULL()
  integer, dimension(4) :: siz
  real, dimension(:,:,:), pointer :: tmp_buffer_in => NULL()  ! Unrotated input [various units]
  integer :: ni_seg, nj_seg  ! number of src gridpoints along the segments
  integer :: ni_buf, nj_buf  ! Number of filled values in tmp_buffer
  integer :: is_obc, ie_obc, js_obc, je_obc  ! segment indices within local domain
  integer :: ishift, jshift  ! offsets for staggered locations
  real    :: dz(SZI_(G),SZJ_(G),SZK_(GV)) ! Distance between the interfaces around a layer [Z ~> m]
  real, dimension(:,:,:), allocatable, target :: tmp_buffer ! A buffer for input data [various units]
  real, dimension(:), allocatable :: dz_stack  ! Distance between the interfaces at corner points [Z ~> m]
  integer :: is_obc2, js_obc2
  integer :: i_seg_offset, j_seg_offset, bug_offset
  real :: net_dz_src  ! Total vertical extent of the incoming flow in the source field [Z ~> m]
  real :: net_dz_int  ! Total vertical extent of the incoming flow in the model [Z ~> m]
  real :: scl_fac     ! A scaling factor to compensate for differences in total thicknesses [nondim]
  real :: tidal_vel   ! Interpolated tidal velocity at the OBC points [L T-1 ~> m s-1]
  real :: tidal_elev  ! Interpolated tidal elevation at the OBC points [Z ~> m]
  real :: ramp_value  ! If OBC%ramp is True, where we are on the ramp from 0 to 1, or 1 otherwise [nondim].
  real, allocatable :: normal_trans_bt(:,:) ! barotropic transport [H L2 T-1 ~> m3 s-1]
  integer :: turns    ! Number of index quarter turns
  real :: time_delta  ! Time since tidal reference date [T ~> s]
  logical :: flip_buffer ! If true, the input buffer needs to be transposed

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  nz=GV%ke

  turns = modulo(G%HI%turns, 4)

  if (.not. associated(OBC)) return

  if (OBC%add_tide_constituents) time_delta = US%s_to_T * time_type_to_real(Time - OBC%time_ref)

  if (OBC%number_of_segments >= 1) then
    dz(:,:,:) = 0.0
    call thickness_to_dz(h, tv, dz, G, GV, US)
    call pass_var(dz, G%Domain)
  endif

  do n = 1, OBC%number_of_segments
    segment => OBC%segment(n)

    if (.not. segment%on_pe) cycle ! continue to next segment if not in computational domain

    ! NOTE: segment%is_obc and segment%ie_obc are range of indices for the full segment.
    !  The other data set here are in segment%HI, but here they defined slightly differently.
    ni_seg = segment%ie_obc-segment%is_obc+1
    nj_seg = segment%je_obc-segment%js_obc+1
    is_obc = max(segment%is_obc,isd-1)
    ie_obc = min(segment%ie_obc,ied)
    js_obc = max(segment%js_obc,jsd-1)
    je_obc = min(segment%je_obc,jed)
    i_seg_offset = G%idg_offset - segment%HI%Isgb
    j_seg_offset = G%jdg_offset - segment%HI%Jsgb

! Calculate auxiliary fields at staggered locations.
! Segment indices are on q points:
!
!       |-----------|------------|-----------|-----------|  J_obc
!     Is_obc                                          Ie_obc
!
! i2 has to start at Is_obc+1 and end at Ie_obc.
! j2 is J_obc and jshift has to be +1 at both the north and south.

    ! calculate auxiliary fields at staggered locations
    ishift = 0 ; jshift = 0
    segment%Htot(:,:) = 0.0
    segment%dZtot(:,:) = 0.0
    if (segment%is_E_or_W) then
      allocate(normal_trans_bt(segment%HI%IsdB:segment%HI%IedB,segment%HI%jsd:segment%HI%jed), source=0.0)
      if (segment%direction == OBC_DIRECTION_W) ishift=1
      I=segment%HI%IsdB
      ! dZtot may extend one point past the end of the segment on the current PE for use at vorticity points
      do k=1,GV%ke ; do j = max(segment%HI%jsd-1,G%jsd), min(segment%HI%jed+1,G%jed)
        segment%dZtot(I,j) = segment%dZtot(I,j) + dz(i+ishift,j,k)
      enddo ; enddo
      do k=1,GV%ke ; do j=segment%HI%jsd,segment%HI%jed
        segment%h(I,j,k) = h(i+ishift,j,k)
        segment%Htot(I,j) = segment%Htot(I,j) + segment%h(I,j,k)
      enddo ; enddo
      do j=segment%HI%jsd,segment%HI%jed
        segment%Cg(I,j) = sqrt(GV%g_prime(1) * max(0.0, segment%dZtot(I,j)))
      enddo
    else ! (segment%direction == OBC_DIRECTION_N .or. segment%direction == OBC_DIRECTION_S)
      allocate(normal_trans_bt(segment%HI%isd:segment%HI%ied,segment%HI%JsdB:segment%HI%JedB), source=0.0)
      if (segment%direction == OBC_DIRECTION_S) jshift=1
      J=segment%HI%JsdB
      ! dZtot may extend one point past the end of the segment on the current PE for use at vorticity points
      do k=1,GV%ke ; do i = max(segment%HI%isd-1,G%isd), min(segment%HI%ied+1,G%ied)
        segment%dZtot(i,J) = segment%dZtot(i,J) + dz(i,j+jshift,k)
      enddo ; enddo
      do k=1,GV%ke ; do i=segment%HI%isd,segment%HI%ied
        segment%h(i,J,k) = h(i,j+jshift,k)
        segment%Htot(i,J) = segment%Htot(i,J) + segment%h(i,J,k)
      enddo ; enddo
      do i=segment%HI%isd,segment%HI%ied
        segment%Cg(i,J) = sqrt(GV%g_prime(1) * max(0.0, segment%dZtot(i,J)))
      enddo
    endif

    allocate(dz_stack(GV%ke), source=0.0)
    do m = 1,segment%num_fields
      !This field may not require a high frequency OBC segment update and might be allowed
      !a less frequent update as set by the parameter update_OBC_period_max in MOM.F90.
      !Cycle if it is not the time to update OBC segment data for this field.
      if (trim(segment%field(m)%genre) == 'obgc' .and. (.not. OBC%update_OBC_seg_data)) cycle
      if (segment%field(m)%use_IO) then
        siz(1) = size(segment%field(m)%buffer_src,1)
        siz(2) = size(segment%field(m)%buffer_src,2)
        siz(3) = size(segment%field(m)%buffer_src,3)
        if (.not.allocated(segment%field(m)%buffer_dst)) then
          if (siz(3) /= segment%field(m)%nk_src) call MOM_error(FATAL,'nk_src inconsistency')

          nk_dst = GV%ke
          if (field_is_tidal(segment%field(m)%name)) nk_dst = siz(3)
          if (segment%field(m)%nk_src <= 1) nk_dst = 1
          if (.not.segment%field(m)%on_face) then
            allocate(segment%field(m)%buffer_dst(is_obc:ie_obc, js_obc:je_obc, nk_dst), source=0.0)
          elseif (segment%is_E_or_W) then
            allocate(segment%field(m)%buffer_dst(is_obc:ie_obc, js_obc+1:je_obc, nk_dst), source=0.0)
          else
            allocate(segment%field(m)%buffer_dst(is_obc+1:ie_obc, js_obc:je_obc, nk_dst), source=0.0)
          endif
        endif
        ! read source data interpolated to the current model time
        ! NOTE: buffer is sized for vertex points, but may be used for faces
        if (siz(1)==1) then
          if (OBC%brushcutter_mode) then
            allocate(tmp_buffer(1,nj_seg*2-1,segment%field(m)%nk_src))  ! segment data is currently on supergrid
          else
            allocate(tmp_buffer(1,nj_seg,segment%field(m)%nk_src))  ! segment data is currently on native grid
          endif
        else
          if (OBC%brushcutter_mode) then
            allocate(tmp_buffer(ni_seg*2-1,1,segment%field(m)%nk_src))  ! segment data is currently on supergrid
          else
            allocate(tmp_buffer(ni_seg,1,segment%field(m)%nk_src))  ! segment data is currently on native grid
          endif
        endif

        ! TODO: Since we conditionally rotate a subset of tmp_buffer_in after
        !   reading the value, it is currently not possible to use the rotated
        !   implementation of time_interp_extern.
        !   For now, we must explicitly allocate and rotate this array.
        if (turns /= 0) then
          if (modulo(turns, 2) /= 0) then
            allocate(tmp_buffer_in(size(tmp_buffer, 2), size(tmp_buffer, 1), size(tmp_buffer, 3)))
          else
            allocate(tmp_buffer_in(size(tmp_buffer, 1), size(tmp_buffer, 2), size(tmp_buffer, 3)))
          endif
        else
          tmp_buffer_in => tmp_buffer
        endif

        ! This is where the data values are actually read in.
        call time_interp_external(segment%field(m)%handle, Time, tmp_buffer_in, scale=segment%field(m)%scale)

        ! NOTE: Rotation of face-points require that we skip the final value when not in brushcutter mode.
        if (turns /= 0) then
          flip_buffer = ((turns==1) .or. (turns==3))
          if (OBC%brushcutter_mode .or. (.not.flip_buffer)) then
            call rotate_array(tmp_buffer_in, turns, tmp_buffer)
          elseif (flip_buffer .and. segment%is_E_or_W .and. segment%field(m)%on_face) then
            nj_buf = size(tmp_buffer, 2) - 1
            call rotate_array(tmp_buffer_in(:nj_buf,:,:), turns, tmp_buffer(:,:nj_buf,:))
          elseif (flip_buffer .and. segment%is_N_or_S .and. segment%field(m)%on_face) then
            ni_buf = size(tmp_buffer, 1) - 1
            call rotate_array(tmp_buffer_in(:,:ni_buf,:), turns, tmp_buffer(:ni_buf,:,:))
          else
            call rotate_array(tmp_buffer_in, turns, tmp_buffer)
          endif

          if (((segment%field(m)%name == 'U') .and. ((turns==1).or.(turns==2))) .or. &
              ((segment%field(m)%name == 'V') .and. ((turns==2).or.(turns==3))) .or. &
              ((segment%field(m)%name == 'Vamp') .and. ((turns==2).or.(turns==3))) .or. &
              ((segment%field(m)%name == 'Uamp') .and. ((turns==1).or.(turns==2))) .or. &
              ((segment%field(m)%name == 'DVDX') .and. ((turns==1).or.(turns==3))) .or. &
              ((segment%field(m)%name == 'DUDY') .and. ((turns==1).or.(turns==3))) ) then
            tmp_buffer(:,:,:) = -tmp_buffer(:,:,:)
          endif
        endif

        if (OBC%brushcutter_mode) then
          ! In brushcutter mode, the input data includes vales at both the vorticity point nodes and
          ! the velocity point faces of the OBC segments.  The vorticity node values are at the odd
          ! positions in tmp_buffer, while the faces are at the even points.  The bug that is being
          ! corrected here is the use of the odd indexed points for both the corners and the faces.
          bug_offset = 0 ; if (OBC%hor_index_bug) bug_offset = -1
          if (segment%is_E_or_W) then
            if (.not.segment%field(m)%on_face) then
              segment%field(m)%buffer_src(is_obc,:,:) = &
                  tmp_buffer(1, 2*(js_obc+j_seg_offset+1)-1:2*(je_obc+j_seg_offset)+1:2, :)
            else
              segment%field(m)%buffer_src(is_obc,:,:) = &
                  tmp_buffer(1, 2*(js_obc+j_seg_offset+1)+bug_offset:2*(je_obc+j_seg_offset):2, :)
            endif
          else
            if (.not.segment%field(m)%on_face) then
              segment%field(m)%buffer_src(:,js_obc,:) = &
                  tmp_buffer(2*(is_obc+i_seg_offset+1)-1:2*(ie_obc+i_seg_offset)+1:2, 1, :)
            else
              segment%field(m)%buffer_src(:,js_obc,:) = &
                  tmp_buffer(2*(is_obc+i_seg_offset+1)+bug_offset:2*(ie_obc+i_seg_offset):2, 1, :)
            endif
          endif
        else  ! Not brushcutter_mode.
          if (segment%is_E_or_W) then
            if (.not.segment%field(m)%on_face) then
              segment%field(m)%buffer_src(is_obc,:,:) = &
                   tmp_buffer(1,js_obc+j_seg_offset+1:je_obc+j_seg_offset+1,:)
            else
              segment%field(m)%buffer_src(is_obc,:,:) = &
                   tmp_buffer(1,js_obc+j_seg_offset+1:je_obc+j_seg_offset,:)
            endif
          else
            if (.not.segment%field(m)%on_face) then
              segment%field(m)%buffer_src(:,js_obc,:) = &
                   tmp_buffer(is_obc+i_seg_offset+1:ie_obc+i_seg_offset+1,1,:)
            else
              segment%field(m)%buffer_src(:,js_obc,:) = &
                   tmp_buffer(is_obc+i_seg_offset+1:ie_obc+i_seg_offset,1,:)
            endif
          endif
        endif

        ! no dz for tidal variables
        if (segment%field(m)%nk_src <= 1) then  ! This is 2-d data with no remapping.
          segment%field(m)%buffer_dst(:,:,1) = segment%field(m)%buffer_src(:,:,1)
        elseif (field_is_tidal(segment%field(m)%name)) then
          ! The 3rd axis for tidal variables is the tidal constituent, so there is no remapping.
          segment%field(m)%buffer_dst(:,:,:) = segment%field(m)%buffer_src(:,:,:)
        else
          ! Read in 3-d data that may need to be remapped onto the new grid
          ! This is also where the 2-d tidal data values (apart from phase and amp) are actually read in.
          call time_interp_external(segment%field(m)%dz_handle, Time, tmp_buffer_in, scale=US%m_to_Z)

          if (turns /= 0) then
            flip_buffer = ((turns==1) .or. (turns==3))
            if (flip_buffer .and. segment%is_E_or_W .and. segment%field(m)%on_face) then
              nj_buf = size(tmp_buffer, 2) - 1
              call rotate_array(tmp_buffer_in(:nj_buf,:,:), turns, tmp_buffer(:,:nj_buf,:))
            elseif (flip_buffer .and. segment%is_N_or_S .and. segment%field(m)%on_face) then
              ni_buf = size(tmp_buffer, 1) - 1
              call rotate_array(tmp_buffer_in(:,:ni_buf,:), turns, tmp_buffer(:ni_buf,:,:))
            else
              call rotate_array(tmp_buffer_in, turns, tmp_buffer)
            endif
          endif ! End of rotation

          if (OBC%brushcutter_mode) then
            bug_offset = 0 ; if (OBC%hor_index_bug) bug_offset = -1
            if (segment%is_E_or_W) then
              if (.not.segment%field(m)%on_face) then
                segment%field(m)%dz_src(is_obc,:,:) = &
                    tmp_buffer(1, 2*(js_obc+j_seg_offset+1)-1:2*(je_obc+j_seg_offset)+1:2, :)
              else
                segment%field(m)%dz_src(is_obc,:,:) = &
                    tmp_buffer(1, 2*(js_obc+j_seg_offset+1)+bug_offset:2*(je_obc+j_seg_offset):2, :)
              endif
            else
              if (.not.segment%field(m)%on_face) then
                segment%field(m)%dz_src(:,js_obc,:) = &
                    tmp_buffer(2*(is_obc+i_seg_offset+1)-1:2*(ie_obc+i_seg_offset)+1:2, 1, :)
              else
                segment%field(m)%dz_src(:,js_obc,:) = &
                    tmp_buffer(2*(is_obc+i_seg_offset+1)+bug_offset:2*(ie_obc+i_seg_offset):2, 1, :)
              endif
            endif
          else  ! Not brushcutter_mode.
            if (segment%is_E_or_W) then
              if (.not.segment%field(m)%on_face) then
                segment%field(m)%dz_src(is_obc,:,:) = &
                    tmp_buffer(1,js_obc+j_seg_offset+1:je_obc+j_seg_offset+1,:)
              else
                segment%field(m)%dz_src(is_obc,:,:) = &
                    tmp_buffer(1,js_obc+j_seg_offset+1:je_obc+j_seg_offset,:)
              endif
            else
              if (.not.segment%field(m)%on_face) then
                segment%field(m)%dz_src(:,js_obc,:) = &
                    tmp_buffer(is_obc+i_seg_offset+1:ie_obc+i_seg_offset+1,1,:)
              else
                segment%field(m)%dz_src(:,js_obc,:) = &
                    tmp_buffer(is_obc+i_seg_offset+1:ie_obc+i_seg_offset,1,:)
              endif
            endif
          endif

          if ((.not.segment%field(m)%on_face) .and. (.not.OBC%hor_index_bug)) then
            ! This point is at the OBC vorticity point nodes, rather than the OBC velocity point faces.
            call adjustSegmentEtaToFitBathymetry(G, GV, US, segment, m, at_node=.true.)
          else
            call adjustSegmentEtaToFitBathymetry(G, GV, US, segment, m, at_node=.false.)
          endif

          if (segment%is_E_or_W) then
            ishift=1
            if (segment%direction == OBC_DIRECTION_E) ishift=0
            I=is_obc
            if (.not.segment%field(m)%on_face) then
              ! Do q points for the whole segment
              do J=max(js_obc,jsd),min(je_obc,jed-1)
                ! Using the h remapping approach
                ! Pretty sure we need to check for source/target grid consistency here
                !### For a concave corner between OBC segments, there are 3 thicknesses we might
                ! consider using.
                segment%field(m)%buffer_dst(I,J,:) = 0.0  ! initialize remap destination buffer
                if (G%mask2dCu(I,j)>0. .and. G%mask2dCu(I,j+1)>0.) then
                  dz_stack(:) = 0.5*(dz(i+ishift,j,:) + dz(i+ishift,j+1,:))
                  call remapping_core_h(OBC%remap_z_CS, &
                       segment%field(m)%nk_src, segment%field(m)%dz_src(I,J,:), &
                       segment%field(m)%buffer_src(I,J,:), &
                       GV%ke, dz_stack, segment%field(m)%buffer_dst(I,J,:))
                elseif (G%mask2dCu(I,j)>0.) then
                  dz_stack(:) = dz(i+ishift,j,:)
                  call remapping_core_h(OBC%remap_z_CS, &
                       segment%field(m)%nk_src, segment%field(m)%dz_src(I,J,:), &
                       segment%field(m)%buffer_src(I,J,:), &
                       GV%ke, dz_stack, segment%field(m)%buffer_dst(I,J,:))
                elseif (G%mask2dCu(I,j+1)>0.) then
                  dz_stack(:) = dz(i+ishift,j+1,:)
                  call remapping_core_h(OBC%remap_z_CS, &
                       segment%field(m)%nk_src, segment%field(m)%dz_src(I,j,:), &
                       segment%field(m)%buffer_src(I,J,:), &
                       GV%ke, dz_stack, segment%field(m)%buffer_dst(I,J,:))
                endif
              enddo
            else
              do j=js_obc+1,je_obc
                ! Using the h remapping approach
                ! Pretty sure we need to check for source/target grid consistency here
                segment%field(m)%buffer_dst(I,j,:) = 0.0  ! initialize remap destination buffer
                if (G%mask2dCu(I,j)>0.) then
                  net_dz_src = sum( segment%field(m)%dz_src(I,j,:) )
                  net_dz_int = sum( dz(i+ishift,j,:) )
                  scl_fac = net_dz_int / net_dz_src
                  call remapping_core_h(OBC%remap_z_CS, &
                       segment%field(m)%nk_src,  scl_fac*segment%field(m)%dz_src(I,j,:), &
                       segment%field(m)%buffer_src(I,j,:), &
                       GV%ke, dz(i+ishift,j,:), segment%field(m)%buffer_dst(I,j,:))
                endif
              enddo
            endif
          else
            jshift=1
            if (segment%direction == OBC_DIRECTION_N) jshift=0
            J=js_obc
            if (.not.segment%field(m)%on_face) then
              ! Do q points for the whole segment
              do I=max(is_obc,isd),min(ie_obc,ied-1)
                segment%field(m)%buffer_dst(I,J,:) = 0.0  ! initialize remap destination buffer
                if (G%mask2dCv(i,J)>0. .and. G%mask2dCv(i+1,J)>0.) then
              ! Using the h remapping approach
              ! Pretty sure we need to check for source/target grid consistency here
                  dz_stack(:) = 0.5*(dz(i,j+jshift,:) + dz(i+1,j+jshift,:))
                  call remapping_core_h(OBC%remap_z_CS, &
                       segment%field(m)%nk_src, segment%field(m)%dz_src(I,J,:), &
                       segment%field(m)%buffer_src(I,J,:), &
                       GV%ke, dz_stack, segment%field(m)%buffer_dst(I,J,:))
                elseif (G%mask2dCv(i,J)>0.) then
                  dz_stack(:) = dz(i,j+jshift,:)
                  call remapping_core_h(OBC%remap_z_CS, &
                       segment%field(m)%nk_src, segment%field(m)%dz_src(I,J,:), &
                       segment%field(m)%buffer_src(I,J,:), &
                       GV%ke, dz_stack, segment%field(m)%buffer_dst(I,J,:))
                elseif (G%mask2dCv(i+1,J)>0.) then
                  dz_stack(:) = dz(i+1,j+jshift,:)
                  call remapping_core_h(OBC%remap_z_CS, &
                       segment%field(m)%nk_src, segment%field(m)%dz_src(I,J,:), &
                       segment%field(m)%buffer_src(I,J,:), &
                       GV%ke, dz_stack, segment%field(m)%buffer_dst(I,J,:))
                endif
              enddo
            else
              do i=is_obc+1,ie_obc
              ! Using the h remapping approach
              ! Pretty sure we need to check for source/target grid consistency here
                segment%field(m)%buffer_dst(i,J,:) = 0.0  ! initialize remap destination buffer
                if (G%mask2dCv(i,J)>0.) then
                  net_dz_src = sum( segment%field(m)%dz_src(i,J,:) )
                  net_dz_int = sum( dz(i,j+jshift,:) )
                  scl_fac = net_dz_int / net_dz_src
                  call remapping_core_h(OBC%remap_z_CS, &
                       segment%field(m)%nk_src, scl_fac* segment%field(m)%dz_src(i,J,:), &
                       segment%field(m)%buffer_src(i,J,:), &
                       GV%ke, dz(i,j+jshift,:), segment%field(m)%buffer_dst(i,J,:))
                endif
              enddo
            endif
          endif
        endif
        deallocate(tmp_buffer)
        if (turns /= 0) deallocate(tmp_buffer_in)
      else ! use_IO = .false. (Uniform value)
        if (.not. allocated(segment%field(m)%buffer_dst)) then
          nk_dst = GV%ke
          if (field_is_tidal(segment%field(m)%name)) nk_dst = 1
          if (segment%field(m)%name == 'SSH') nk_dst = 1
          if (.not.segment%field(m)%on_face) then
            allocate(segment%field(m)%buffer_dst(is_obc:ie_obc, js_obc:je_obc, nk_dst), &
                     source=segment%field(m)%value)
          elseif (segment%is_E_or_W) then
            allocate(segment%field(m)%buffer_dst(is_obc:ie_obc, js_obc+1:je_obc, nk_dst), &
                     source=segment%field(m)%value)
          else
            allocate(segment%field(m)%buffer_dst(is_obc+1:ie_obc, js_obc:je_obc, nk_dst), &
                     source=segment%field(m)%value)
          endif
        endif
      endif
    enddo ! end field loop

    ! Start second loop to update all fields now that data for all fields are available.
    ! (split because tides depend on multiple variables).
    do m = 1,segment%num_fields
      !cycle if it is not the time to update OBGC tracers from source
      if (trim(segment%field(m)%genre) == 'obgc' .and. (.not. OBC%update_OBC_seg_data)) cycle
      ! if (segment%field(m)%use_IO) then
      ! calculate external BT velocity and transport if needed
      if (trim(segment%field(m)%name) == 'U' .or. trim(segment%field(m)%name) == 'V') then
        if (trim(segment%field(m)%name) == 'U' .and. segment%is_E_or_W) then
          I=is_obc
          do j=js_obc+1,je_obc
            normal_trans_bt(I,j) = 0.0
            tidal_vel = 0.0
            if (OBC%add_tide_constituents) then
              do c=1,OBC%n_tide_constituents
                tidal_vel = tidal_vel + (OBC%tide_fn(c) * segment%field(segment%uamp_index)%buffer_dst(I,j,c)) * &
                  cos((time_delta*OBC%tide_frequencies(c) - segment%field(segment%uphase_index)%buffer_dst(I,j,c)) &
                      + (OBC%tide_eq_phases(c) + OBC%tide_un(c)))
              enddo
            endif
            do k=1,GV%ke
              segment%normal_vel(I,j,k) = segment%field(m)%buffer_dst(I,j,k) + tidal_vel
              segment%normal_trans(I,j,k) = segment%normal_vel(I,j,k)*segment%h(I,j,k) * G%dyCu(I,j)
              normal_trans_bt(I,j) = normal_trans_bt(I,j) + segment%normal_trans(I,j,k)
            enddo
            segment%normal_vel_bt(I,j) = normal_trans_bt(I,j) &
                / (max(segment%Htot(I,j), 1.e-12 * GV%m_to_H) * G%dyCu(I,j))
            if (allocated(segment%nudged_normal_vel)) segment%nudged_normal_vel(I,j,:) = segment%normal_vel(I,j,:)
          enddo
        elseif (trim(segment%field(m)%name) == 'V' .and. segment%is_N_or_S) then
          J=js_obc
          do i=is_obc+1,ie_obc
            normal_trans_bt(i,J) = 0.0
            tidal_vel = 0.0
            if (OBC%add_tide_constituents) then
              do c=1,OBC%n_tide_constituents
                tidal_vel = tidal_vel + (OBC%tide_fn(c) * segment%field(segment%vamp_index)%buffer_dst(I,j,c)) * &
                  cos((time_delta*OBC%tide_frequencies(c) - segment%field(segment%vphase_index)%buffer_dst(I,j,c)) &
                      + (OBC%tide_eq_phases(c) + OBC%tide_un(c)))
              enddo
            endif
            do k=1,GV%ke
              segment%normal_vel(i,J,k) = segment%field(m)%buffer_dst(i,J,k) + tidal_vel
              segment%normal_trans(i,J,k) = segment%normal_vel(i,J,k)*segment%h(i,J,k) * &
                        G%dxCv(i,J)
              normal_trans_bt(i,J) = normal_trans_bt(i,J) + segment%normal_trans(i,J,k)
            enddo
            segment%normal_vel_bt(i,J) = normal_trans_bt(i,J) &
                / (max(segment%Htot(i,J), 1.e-12 * GV%m_to_H) * G%dxCv(i,J))
            if (allocated(segment%nudged_normal_vel)) segment%nudged_normal_vel(i,J,:) = segment%normal_vel(i,J,:)
          enddo
        elseif (trim(segment%field(m)%name) == 'V' .and. segment%is_E_or_W .and. &
                allocated(segment%tangential_vel)) then
          I=is_obc
          do J=js_obc,je_obc
            tidal_vel = 0.0
            if (OBC%add_tide_constituents) then
              do c=1,OBC%n_tide_constituents
                tidal_vel = tidal_vel + (OBC%tide_fn(c) * segment%field(segment%vamp_index)%buffer_dst(I,j,c)) * &
                  cos((time_delta*OBC%tide_frequencies(c) - segment%field(segment%vphase_index)%buffer_dst(I,j,c)) &
                      + (OBC%tide_eq_phases(c) + OBC%tide_un(c)))
              enddo
            endif
            do k=1,GV%ke
              segment%tangential_vel(I,J,k) = segment%field(m)%buffer_dst(I,J,k) + tidal_vel
            enddo
            if (allocated(segment%nudged_tangential_vel)) &
              segment%nudged_tangential_vel(I,J,:) = segment%tangential_vel(I,J,:)
          enddo
        elseif (trim(segment%field(m)%name) == 'U' .and. segment%is_N_or_S .and. &
                allocated(segment%tangential_vel)) then
          J=js_obc
          do I=is_obc,ie_obc
            tidal_vel = 0.0
            if (OBC%add_tide_constituents) then
              do c=1,OBC%n_tide_constituents
                tidal_vel = tidal_vel + (OBC%tide_fn(c) * segment%field(segment%uamp_index)%buffer_dst(I,j,c)) * &
                    cos((time_delta*OBC%tide_frequencies(c) - segment%field(segment%uphase_index)%buffer_dst(I,j,c)) &
                        + (OBC%tide_eq_phases(c) + OBC%tide_un(c)))
              enddo
            endif
            do k=1,GV%ke
              segment%tangential_vel(I,J,k) = segment%field(m)%buffer_dst(I,J,k) + tidal_vel
            enddo
            if (allocated(segment%nudged_tangential_vel)) &
              segment%nudged_tangential_vel(I,J,:) = segment%tangential_vel(I,J,:)
          enddo
        endif
      elseif (trim(segment%field(m)%name) == 'DVDX' .and. segment%is_E_or_W .and. &
              allocated(segment%tangential_grad)) then
        I=is_obc
        do J=js_obc,je_obc
          do k=1,GV%ke
            segment%tangential_grad(I,J,k) = segment%field(m)%buffer_dst(I,J,k)
            if (allocated(segment%nudged_tangential_grad)) &
              segment%nudged_tangential_grad(I,J,:) = segment%tangential_grad(I,J,:)
          enddo
        enddo
      elseif (trim(segment%field(m)%name) == 'DUDY' .and. segment%is_N_or_S .and. &
              allocated(segment%tangential_grad)) then
        J=js_obc
        do I=is_obc,ie_obc
          do k=1,GV%ke
            segment%tangential_grad(I,J,k) = segment%field(m)%buffer_dst(I,J,k)
            if (allocated(segment%nudged_tangential_grad)) &
              segment%nudged_tangential_grad(I,J,:) = segment%tangential_grad(I,J,:)
          enddo
        enddo
      endif

      ! endif

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
        ramp_value = 1.0
        if (OBC%ramp) ramp_value = OBC%ramp_value
        do j=js_obc2,je_obc ; do i=is_obc2,ie_obc
          tidal_elev = 0.0
          if (OBC%add_tide_constituents) then
            do c=1,OBC%n_tide_constituents
              tidal_elev = tidal_elev + (OBC%tide_fn(c) * segment%field(segment%zamp_index)%buffer_dst(i,j,c)) * &
                  cos((time_delta*OBC%tide_frequencies(c) - segment%field(segment%zphase_index)%buffer_dst(i,j,c)) &
                      + (OBC%tide_eq_phases(c) + OBC%tide_un(c)))
            enddo
          endif
          segment%SSH(i,j) = ramp_value * (segment%field(m)%buffer_dst(i,j,1) + tidal_elev)
        enddo ; enddo
      endif

      ! Set the inflow and reservoir data for tracers.
      if ((trim(segment%field(m)%name) == 'TEMP') .or. (trim(segment%field(m)%name) == 'SALT') .or. &
         (trim(segment%field(m)%genre) == 'obgc')) then
        if (trim(segment%field(m)%name) == 'TEMP') then
          nt = 1
        elseif (trim(segment%field(m)%name) == 'SALT') then
          nt = 2
        elseif (trim(segment%field(m)%genre) == 'obgc') then
          nt = get_tracer_index(segment,trim(segment%field(m)%name))
          if (nt < 0) call MOM_error(FATAL,"update_OBC_segment_data: Did not find tracer "//trim(segment%field(m)%name))
        endif
        if (allocated(segment%field(m)%buffer_dst)) then
          do k=1,nz; do j=js_obc2, je_obc; do i=is_obc2,ie_obc
            segment%tr_Reg%Tr(nt)%t(i,j,k) = segment%field(m)%buffer_dst(i,j,k)
          enddo ; enddo ; enddo
          if (.not. segment%tr_Reg%Tr(nt)%is_initialized) then
            ! If the tracer reservoir has not yet been initialized, then set to external value.
            do k=1,nz; do j=js_obc2, je_obc; do i=is_obc2,ie_obc
              segment%tr_Reg%Tr(nt)%tres(i,j,k) = segment%tr_Reg%Tr(nt)%t(i,j,k)
            enddo ; enddo ; enddo
            segment%tr_Reg%Tr(nt)%is_initialized=.true.
          endif
        else
          segment%tr_Reg%Tr(nt)%OBC_inflow_conc = segment%field(m)%value
        endif
      endif

    enddo ! end field loop
    deallocate(dz_stack)
    deallocate(normal_trans_bt)

  enddo ! end segment loop

end subroutine update_OBC_segment_data

!> Update the OBC ramp value as a function of time.
!! If called with the optional argument activate=.true., record the
!! value of Time as the beginning of the ramp period.
subroutine update_OBC_ramp(Time, OBC, US, activate)
  type(time_type), target, intent(in)    :: Time     !< Current model time
  type(ocean_OBC_type),    intent(inout) :: OBC      !< Open boundary structure
  type(unit_scale_type),   intent(in)    :: US       !< A dimensional unit scaling type
  logical, optional,       intent(in)    :: activate !< Specify whether to record the value of
                                                     !! Time as the beginning of the ramp period

  ! Local variables
  real :: deltaTime ! The time since start of ramping [T ~> s]
  real :: wghtA     ! A temporary variable used to set OBC%ramp_value [nondim]
  character(len=12) :: msg

  if (.not. OBC%ramp) return ! This indicates the ramping is turned off

  ! We use the optional argument to indicate this Time should be recorded as the
  ! beginning of the ramp-up period.
  if (present(activate)) then
    if (activate) then
      OBC%ramp_start_time = Time ! Record the current time
      OBC%ramping_is_activated = .true.
      OBC%trunc_ramp_time = OBC%ramp_timescale ! times 3.0 for tanh
    endif
  endif
  if (.not.OBC%ramping_is_activated) return
  deltaTime = max( 0., US%s_to_T*time_type_to_real( Time - OBC%ramp_start_time ) )
  if (deltaTime >= OBC%trunc_ramp_time) then
    OBC%ramp_value = 1.0
    OBC%ramp = .false. ! This turns off ramping after this call
  else
    wghtA = min( 1., deltaTime / OBC%ramp_timescale ) ! Linear profile in time
    !wghtA = wghtA*wghtA ! Convert linear profile to parabolic profile in time
    !wghtA = wghtA*wghtA*(3. - 2.*wghtA) ! Convert linear profile to cosine profile
    !wghtA = 1. - ( (1. - wghtA)**2 ) ! Convert linear profile to inverted parabolic profile
    !wghtA = tanh(wghtA)              ! Convert linear profile to tanh
    OBC%ramp_value = wghtA
  endif
  write(msg(1:12),'(es12.3)') OBC%ramp_value
  call MOM_error(NOTE, "MOM_open_boundary: update_OBC_ramp set OBC"// &
                       " ramp to "//trim(msg))
end subroutine update_OBC_ramp

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
    call MOM_error(FATAL,"MOM register_OBC: "//mesg)
  endif
  Reg%nobc = Reg%nobc + 1
  nobc     = Reg%nobc

  Reg%OB(nobc)%name = name

  if (Reg%locked) call MOM_error(FATAL, &
      "MOM register_OBC was called for OBC "//trim(Reg%OB(nobc)%name)//&
      " with a locked OBC registry.")

end subroutine register_OBC

!> This routine include declares and sets the variable "version".
subroutine OBC_registry_init(param_file, Reg)
  type(param_file_type),   intent(in) :: param_file !< open file to parse for model parameters
  type(OBC_registry_type), pointer    :: Reg        !< pointer to OBC registry

  integer, save :: init_calls = 0

# include "version_variable.h"
  character(len=256) :: mesg    ! Message for error messages.

  if (.not.associated(Reg)) then ; allocate(Reg)
  else ; return ; endif

  ! Read all relevant parameters and write them to the model log.
! call log_version(param_file, mdl, version, "")

  init_calls = init_calls + 1
  if (init_calls > 1) then
    write(mesg,'("OBC_registry_init called ",I3, &
      &" times with different registry pointers.")') init_calls
    if (is_root_pe()) call MOM_error(WARNING,"MOM_open_boundary"//mesg)
  endif

end subroutine OBC_registry_init

!> Add file to OBC registry.
function register_file_OBC(param_file, CS, US, OBC_Reg)
  type(param_file_type),    intent(in) :: param_file !< parameter file.
  type(file_OBC_CS),        pointer    :: CS         !< file control structure.
  type(unit_scale_type),    intent(in) :: US         !< A dimensional unit scaling type
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
# include "version_variable.h"
  character(len=40)  :: mdl = "segment_tracer_registry_init" ! This routine's name.
  !character(len=256) :: mesg    ! Message for error messages.

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

!> Register a tracer array that is active on an OBC segment, potentially also specifying how the
!! tracer inflow values are specified.
subroutine register_segment_tracer(tr_ptr, ntr_index, param_file, GV, segment, &
                                   OBC_scalar, OBC_array, scale, fd_index)
  type(verticalGrid_type), intent(in)   :: GV         !< ocean vertical grid structure
  type(tracer_type), target             :: tr_ptr     !< A target that can be used to set a pointer to the
                                                      !! stored value of tr. This target must be
                                                      !! an enduring part of the control structure,
                                                      !! because the tracer registry will use this memory,
                                                      !! but it also means that any updates to this
                                                      !! structure in the calling module will be
                                                      !! available subsequently to the tracer registry.
  integer, intent(in)                   :: ntr_index  !< index of segment tracer in the global tracer registry
  type(param_file_type),  intent(in)    :: param_file !< file to parse for model parameter values
  type(OBC_segment_type), intent(inout) :: segment    !< current segment data structure
  real,         optional, intent(in)    :: OBC_scalar !< If present, use scalar value for segment tracer
                                                      !! inflow concentration, including any rescaling to
                                                      !! put the tracer concentration into its internal units,
                                                      !! like [S ~> ppt] for salinity.
  logical,      optional, intent(in)    :: OBC_array  !< If true, use array values for segment tracer
                                                      !! inflow concentration.
  real,         optional, intent(in)    :: scale      !< A scaling factor that should be used with any
                                                      !! data that is read in to convert it to the internal
                                                      !! units of this tracer, in units like [S ppt-1 ~> 1]
                                                      !! for salinity.
  integer,      optional, intent(in)    :: fd_index   !< index of segment tracer in the input field

! Local variables
  real :: rescale ! A multiplicatively corrected scaling factor, in units like [S ppt-1 ~> 1] for
                  ! salinity, or other various units depending on what rescaling has occurred previously.
  integer :: ntseg, m, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  character(len=256) :: mesg    ! Message for error messages.

  call segment_tracer_registry_init(param_file, segment)

  if (segment%tr_Reg%ntseg>=MAX_FIELDS_) then
    write(mesg,'("Increase MAX_FIELDS_ in MOM_memory.h to at least ",I3," to allow for &
        &all the tracers being registered via register_segment_tracer.")') segment%tr_Reg%ntseg+1
    call MOM_error(FATAL,"MOM register_segment_tracer: "//mesg)
  endif
  segment%tr_Reg%ntseg = segment%tr_Reg%ntseg + 1
  ntseg     = segment%tr_Reg%ntseg

  isd = segment%HI%isd ; ied = segment%HI%ied
  jsd = segment%HI%jsd ; jed = segment%HI%jed
  IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
  JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB

  segment%tr_Reg%Tr(ntseg)%Tr => tr_ptr
  segment%tr_Reg%Tr(ntseg)%name = tr_ptr%name
  segment%tr_Reg%Tr(ntseg)%ntr_index = ntr_index
  if (present(fd_index)) segment%tr_Reg%Tr(ntseg)%fd_index = fd_index

  segment%tr_Reg%Tr(ntseg)%scale = 1.0
  if (present(scale)) then
    segment%tr_Reg%Tr(ntseg)%scale = scale
    do m=1,segment%num_fields
      ! Store the scaling factor for fields with exactly matching names, and possibly
      ! rescale the previously stored input values.  Note that calls to register_segment_tracer
      ! can come before or after calls to initialize_segment_data.
      if (uppercase(segment%field(m)%name) == uppercase(segment%tr_Reg%Tr(ntseg)%name)) then
        if (.not. segment%field(m)%use_IO) then
          rescale = scale
          if ((segment%field(m)%scale /= 0.0) .and. (segment%field(m)%scale /= 1.0)) &
            rescale = scale / segment%field(m)%scale
          segment%field(m)%value = rescale * segment%field(m)%value
        endif
        segment%field(m)%scale = scale
      endif
    enddo
  endif

  if (segment%tr_Reg%locked) call MOM_error(FATAL, &
      "MOM register_segment_tracer was called for variable "//trim(segment%tr_Reg%Tr(ntseg)%name)//&
      " with a locked tracer registry.")

  if (present(OBC_scalar)) segment%tr_Reg%Tr(ntseg)%OBC_inflow_conc = OBC_scalar ! initialize tracer value later
  if (present(OBC_array)) then
    if (segment%is_E_or_W) then
      allocate(segment%tr_Reg%Tr(ntseg)%t(IsdB:IedB,jsd:jed,1:GV%ke), source=0.0)
      allocate(segment%tr_Reg%Tr(ntseg)%tres(IsdB:IedB,jsd:jed,1:GV%ke), source=0.0)
      segment%tr_Reg%Tr(ntseg)%is_initialized = .false.
    elseif (segment%is_N_or_S) then
      allocate(segment%tr_Reg%Tr(ntseg)%t(isd:ied,JsdB:JedB,1:GV%ke), source=0.0)
      allocate(segment%tr_Reg%Tr(ntseg)%tres(isd:ied,JsdB:JedB,1:GV%ke), source=0.0)
      segment%tr_Reg%Tr(ntseg)%is_initialized = .false.
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
      if (allocated(Reg%Tr(n)%t)) deallocate(Reg%Tr(n)%t)
    enddo
    deallocate(Reg)
  endif
end subroutine segment_tracer_registry_end

!> Registers the temperature and salinity in the segment tracer registry.
subroutine register_temp_salt_segments(GV, US, OBC, tr_Reg, param_file)
  type(verticalGrid_type),    intent(in)    :: GV         !< ocean vertical grid structure
  type(unit_scale_type),      intent(in)    :: US         !< Unit scaling type
  type(ocean_OBC_type),       pointer       :: OBC        !< Open boundary structure
  type(tracer_registry_type), pointer       :: tr_Reg     !< Tracer registry
  type(param_file_type),      intent(in)    :: param_file !< file to parse for  model parameter values

  ! Local variables
  integer :: n, ntr_id
  character(len=32) :: name
  type(OBC_segment_type), pointer :: segment => NULL() ! pointer to segment type list
  type(tracer_type), pointer :: tr_ptr => NULL()

  if (.not. associated(OBC)) return

  do n=1, OBC%number_of_segments
    segment=>OBC%segment(n)
    if (.not. segment%on_pe) cycle

    if (associated(segment%tr_Reg)) &
         call MOM_error(FATAL,"register_temp_salt_segments: tracer array was previously allocated")

    name = 'temp'
    call tracer_name_lookup(tr_Reg, ntr_id, tr_ptr, name)
    call register_segment_tracer(tr_ptr, ntr_id, param_file, GV, segment, &
                                 OBC_array=segment%temp_segment_data_exists, scale=US%degC_to_C)
    name = 'salt'
    call tracer_name_lookup(tr_Reg, ntr_id, tr_ptr, name)
    call register_segment_tracer(tr_ptr, ntr_id, param_file, GV, segment, &
                                 OBC_array=segment%salt_segment_data_exists, scale=US%ppt_to_S)
  enddo

end subroutine register_temp_salt_segments

!> Sets the OBC properties of external obgc tracers, such as their source file and field name
subroutine set_obgc_segments_props(OBC,tr_name,obc_src_file_name,obc_src_field_name,lfac_in,lfac_out)
  type(ocean_OBC_type),pointer  :: OBC                !< Open boundary structure
  character(len=*),  intent(in) :: tr_name            !< Tracer name
  character(len=*),  intent(in) :: obc_src_file_name  !< OBC source file name
  character(len=*),  intent(in) :: obc_src_field_name !< name of the field in the source file
  real,              intent(in) :: lfac_in            !< factors for tracer reservoir inbound length scales [nondim]
  real,              intent(in) :: lfac_out           !< factors for tracer reservoir outbound length scales [nondim]

  type(external_tracers_segments_props),pointer :: node_ptr => NULL() !pointer to type that keeps
                                                                    ! the tracer segment properties
  allocate(node_ptr)
  node_ptr%tracer_name = trim(tr_name)
  node_ptr%tracer_src_file = trim(obc_src_file_name)
  node_ptr%tracer_src_field = trim(obc_src_field_name)
  node_ptr%lfac_in  = lfac_in
  node_ptr%lfac_out = lfac_out
  ! Reversed Linked List implementation! Make this new node to be the head of the list.
  node_ptr%next => OBC%obgc_segments_props
  OBC%obgc_segments_props => node_ptr
  OBC%num_obgc_tracers = OBC%num_obgc_tracers+1
end subroutine set_obgc_segments_props

!> Get the OBC properties of external obgc tracers, such as their source file, field name,
!! reservoir length scale factors
subroutine get_obgc_segments_props(node, tr_name,obc_src_file_name,obc_src_field_name,lfac_in,lfac_out)
  type(external_tracers_segments_props),pointer :: node !< pointer to tracer segment properties
  character(len=*), intent(out) :: tr_name            !< Tracer name
  character(len=*), intent(out) :: obc_src_file_name  !< OBC source file name
  character(len=*), intent(out) :: obc_src_field_name !< name of the field in the source file
  real,             intent(out) :: lfac_in   !< multiplicative factor for inbound  reservoir length scale [nondim]
  real,             intent(out) :: lfac_out  !< multiplicative factor for outbound reservoir length scale [nondim]
  tr_name = trim(node%tracer_name)
  obc_src_file_name = trim(node%tracer_src_file)
  obc_src_field_name = trim(node%tracer_src_field)
  lfac_in = node%lfac_in
  lfac_out = node%lfac_out
  node => node%next
end subroutine get_obgc_segments_props

!> Registers a named tracer in the segment tracer registries for the OBC segments on which it is active.
subroutine register_obgc_segments(GV, OBC, tr_Reg, param_file, tr_name)
  type(verticalGrid_type),    intent(in)    :: GV         !< ocean vertical grid structure
  type(ocean_OBC_type),       pointer       :: OBC        !< Open boundary structure
  type(tracer_registry_type), pointer       :: tr_Reg     !< Tracer registry
  type(param_file_type),      intent(in)    :: param_file !< file to parse for  model parameter values
  character(len=*),           intent(in)    :: tr_name    !< Tracer name
! Local variables
  integer :: isd, ied, IsdB, IedB, jsd, jed, JsdB, JedB, nz, nf, ntr_id, fd_id
  integer :: i, j, k, n, m
  type(OBC_segment_type), pointer :: segment => NULL() ! pointer to segment type list
  type(tracer_type), pointer      :: tr_ptr => NULL()

  if (.not. associated(OBC)) return

  do n=1, OBC%number_of_segments
    segment=>OBC%segment(n)
    if (.not. segment%on_pe) cycle
    call tracer_name_lookup(tr_Reg, ntr_id, tr_ptr, tr_name)
    ! get the obgc field index
    fd_id = -1
    do m=1,segment%num_fields
      if (lowercase(segment%field(m)%name) == lowercase(tr_name)) fd_id = m
    enddo
    call register_segment_tracer(tr_ptr, ntr_id, param_file, GV, segment, OBC_array=.True., fd_index=fd_id)
  enddo

end subroutine register_obgc_segments

!> Stores the interior tracer values on the segment, and in some cases also sets the tracer reservoir values.
subroutine fill_obgc_segments(G, GV, OBC, tr_ptr, tr_name)
  type(ocean_grid_type),      intent(inout) :: G      !< Ocean grid structure
  type(verticalGrid_type),    intent(in)    :: GV     !< ocean vertical grid structure
  type(ocean_OBC_type),       pointer       :: OBC    !< Open boundary structure
  real, dimension(:,:,:),     pointer       :: tr_ptr !< Pointer to tracer field in scaled concentration
                                                      !! units, like [S ~> ppt] for salinity.
  character(len=*),           intent(in)    :: tr_name !< Tracer name
! Local variables
  integer :: isd, ied, IsdB, IedB, jsd, jed, JsdB, JedB, n, nz, nt
  integer :: i, j, k
  type(OBC_segment_type), pointer :: segment => NULL() ! pointer to segment type list
  real :: I_scale  ! A factor that unscales the internal units of a tracer, like [ppt S-1 ~> 1] for salinity

  if (.not. associated(OBC)) return
  call pass_var(tr_ptr, G%Domain)
  nz = G%ke
  do n=1, OBC%number_of_segments
    segment => OBC%segment(n)
    if (.not. segment%on_pe) cycle
    nt = get_tracer_index(segment, tr_name)
    if (nt < 0) then
      call MOM_error(FATAL,"fill_obgc_segments: Did not find tracer "// tr_name)
    endif
    isd = segment%HI%isd ; ied = segment%HI%ied
    jsd = segment%HI%jsd ; jed = segment%HI%jed
    IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
    JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB

    ! Fill segments with Tracer values
    if (segment%direction == OBC_DIRECTION_W) then
      I = segment%HI%IsdB
      do k=1,nz ; do j=segment%HI%jsd,segment%HI%jed
        segment%tr_Reg%Tr(nt)%t(I,j,k) = tr_ptr(i+1,j,k)
      enddo ; enddo
    elseif (segment%direction == OBC_DIRECTION_E) then
      I = segment%HI%IsdB
      do k=1,nz ; do j=segment%HI%jsd,segment%HI%jed
        segment%tr_Reg%Tr(nt)%t(I,j,k) = tr_ptr(i,j,k)
      enddo ; enddo
    elseif (segment%direction == OBC_DIRECTION_S) then
      J = segment%HI%JsdB
      do k=1,nz ; do i=segment%HI%isd,segment%HI%ied
        segment%tr_Reg%Tr(nt)%t(i,J,k) = tr_ptr(i,j+1,k)
      enddo ; enddo
    elseif (segment%direction == OBC_DIRECTION_N) then
      J = segment%HI%JsdB
      do k=1,nz ; do i=segment%HI%isd,segment%HI%ied
        segment%tr_Reg%Tr(nt)%t(i,J,k) = tr_ptr(i,j,k)
      enddo ; enddo
    endif

    if (.not.segment%tr_Reg%Tr(nt)%is_initialized) &
      segment%tr_Reg%Tr(nt)%tres(:,:,:) = segment%tr_Reg%Tr(nt)%t(:,:,:)

    if (OBC%reservoir_init_bug) then
      ! OBC%tres_x and OBC%tres_y should not be set here, but in a subsequent call to setup_OBC_tracer_reservoirs.
      ! Note that fill_obgc_segments is not called for runs that start from a restart file.
      I_scale = 1.0
      if (segment%tr_Reg%Tr(nt)%scale /= 0.0) I_scale = 1.0 / segment%tr_Reg%Tr(nt)%scale
      if (segment%is_E_or_W) then
        if (allocated(OBC%tres_x)) then
          I = segment%HI%IsdB
          do k=1,nz ; do j=segment%HI%jsd,segment%HI%jed
            OBC%tres_x(I,j,k,nt) = I_scale * segment%tr_Reg%Tr(nt)%tres(I,j,k)
          enddo ; enddo
        endif
      else  ! segment%is_N_or_S
        if (allocated(OBC%tres_y)) then
          J = segment%HI%JsdB
          do k=1,nz ; do i=segment%HI%isd,segment%HI%ied
            OBC%tres_y(i,J,k,nt) = I_scale * segment%tr_Reg%Tr(nt)%tres(i,J,k)
          enddo ; enddo
        endif
      endif
    endif

  enddo ! End of loop over segments.

end subroutine fill_obgc_segments

!> Set the value of temperatures and salinities on OBC segments
subroutine fill_temp_salt_segments(G, GV, US, OBC, tv)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US  !< Unit scaling
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundary structure
  type(thermo_var_ptrs),   intent(in)    :: tv  !< Thermodynamics structure

  integer :: isd, ied, IsdB, IedB, jsd, jed, JsdB, JedB, n, nz
  integer :: i, j, k
  type(OBC_segment_type), pointer :: segment => NULL() ! pointer to segment type list

  if (.not. associated(OBC)) return
  if (.not. associated(tv%T) .and. associated(tv%S)) return
  ! Both temperature and salinity fields

  nz = GV%ke

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
    if (.not.segment%tr_Reg%Tr(1)%is_initialized) &
      segment%tr_Reg%Tr(1)%tres(:,:,:) = segment%tr_Reg%Tr(1)%t(:,:,:)
    if (.not.segment%tr_Reg%Tr(2)%is_initialized) &
      segment%tr_Reg%Tr(2)%tres(:,:,:) = segment%tr_Reg%Tr(2)%t(:,:,:)
  enddo

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
  integer :: i, j
  integer :: l_seg
  logical :: fatal_error = .False.
  real    :: min_depth  ! The minimum depth for ocean points [Z ~> m]
  real    :: mask_depth ! The masking depth for ocean points [Z ~> m]
  real    :: Dmask      ! The depth for masking in the same units as G%bathyT [Z ~> m].
  integer, parameter :: cin = 3, cout = 4, cland = -1, cedge = -2
  character(len=256) :: mesg    ! Message for error messages.
  real, allocatable, dimension(:,:) :: color, color2  ! For sorting inside from outside,
                                                      ! two different ways [nondim]

  if (.not. associated(OBC)) return

  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 units="m", default=0.0, scale=US%m_to_Z, do_not_log=.true.)
  call get_param(param_file, mdl, "MASKING_DEPTH", mask_depth, &
                 units="m", default=-9999.0, scale=US%m_to_Z, do_not_log=.true.)

  Dmask = mask_depth
  if (mask_depth == -9999.0*US%m_to_Z) Dmask = min_depth

  ! The reference depth on a dyn_horgrid is 0, otherwise would need:  min_depth = min_depth - G%Z_ref

  allocate(color(G%isd:G%ied, G%jsd:G%jed), source=0.0)
  allocate(color2(G%isd:G%ied, G%jsd:G%jed), source=0.0)

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
    if (OBC%segnum_u(I,j) < 0) then      !  OBC_DIRECTION_W
      if (color(i,j) == 0.0) color(i,j) = cout
      if (color(i+1,j) == 0.0) color(i+1,j) = cin
    elseif (OBC%segnum_u(I,j) > 0) then  !  OBC_DIRECTION_E
      if (color(i,j) == 0.0) color(i,j) = cin
      if (color(i+1,j) == 0.0) color(i+1,j) = cout
    endif
  enddo ; enddo
  do J=G%JsdB+1,G%JedB-1 ; do i=G%isd,G%ied
    if (OBC%segnum_v(i,J) < 0) then      ! OBC_DIRECTION_S
      if (color(i,j) == 0.0) color(i,j) = cout
      if (color(i,j+1) == 0.0) color(i,j+1) = cin
    elseif (OBC%segnum_v(i,J) > 0) then  ! OBC_DIRECTION_N
      if (color(i,j) == 0.0) color(i,j) = cin
      if (color(i,j+1) == 0.0) color(i,j+1) = cout
    endif
  enddo ; enddo

  do J=G%JsdB+1,G%JedB-1 ; do i=G%isd,G%ied
    if (OBC%segnum_v(i,J) < 0) then      ! OBC_DIRECTION_S
      if (color2(i,j) == 0.0) color2(i,j) = cout
      if (color2(i,j+1) == 0.0) color2(i,j+1) = cin
    elseif (OBC%segnum_v(i,J) > 0) then  ! OBC_DIRECTION_N
      if (color2(i,j) == 0.0) color2(i,j) = cin
      if (color2(i,j+1) == 0.0) color2(i,j+1) = cout
    endif
  enddo ; enddo
  do j=G%jsd,G%jed ; do i=G%IsdB+1,G%IedB-1
    if (OBC%segnum_u(I,j) < 0) then      !  OBC_DIRECTION_W
      if (color2(i,j) == 0.0) color2(i,j) = cout
      if (color2(i+1,j) == 0.0) color2(i+1,j) = cin
    elseif (OBC%segnum_u(I,j) > 0) then  !  OBC_DIRECTION_E
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
          &"the masking of the outside grid points.")') i, j
      call MOM_error(WARNING,"MOM mask_outside_OBCs: "//mesg, all_print=.true.)
    endif
    if (color(i,j) == cout) G%bathyT(i,j) = Dmask
  enddo ; enddo
  if (fatal_error) call MOM_error(FATAL, &
      "MOM_open_boundary: inconsistent OBC segments.")

  deallocate(color)
  deallocate(color2)
end subroutine mask_outside_OBCs

!> flood the cin, cout values
subroutine flood_fill(G, color, cin, cout, cland)
  type(dyn_horgrid_type), intent(inout) :: G      !< Ocean grid structure
  real, dimension(:,:),   intent(inout) :: color  !< For sorting inside from outside [nondim]
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
  real, dimension(:,:),   intent(inout) :: color   !< For sorting inside from outside [nondim]
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
subroutine open_boundary_register_restarts(HI, GV, US, OBC, Reg, param_file, restart_CS, &
                                           use_temperature)
  type(hor_index_type),    intent(in) :: HI !< Horizontal indices
  type(verticalGrid_type), pointer    :: GV !< Container for vertical grid information
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  type(ocean_OBC_type),    pointer    :: OBC !< OBC data structure, data intent(inout)
  type(tracer_registry_type), pointer :: Reg !< pointer to tracer registry
  type(param_file_type),   intent(in) :: param_file !< Parameter file handle
  type(MOM_restart_CS),    intent(inout) :: restart_CS !< MOM restart control structure
  logical,                 intent(in) :: use_temperature !< If true, T and S are used
  ! Local variables
  type(vardesc) :: vd(2)
  integer       :: m
  character(len=100) :: mesg, var_name

  if (.not. associated(OBC)) &
    call MOM_error(FATAL, "open_boundary_register_restarts: Called with "//&
                      "uninitialized OBC control structure")

  ! ### This is a temporary work around for restarts with OBC segments.
  ! This implementation uses 3D arrays solely for restarts. We need
  ! to be able to add 2D ( x,z or y,z ) data to restarts to avoid using
  ! so much memory and disk space.
  if (OBC%radiation_BCs_exist_globally) then
    allocate(OBC%rx_normal(HI%isdB:HI%iedB,HI%jsd:HI%jed,GV%ke), source=0.0)
    allocate(OBC%ry_normal(HI%isd:HI%ied,HI%jsdB:HI%jedB,GV%ke), source=0.0)

    vd(1) = var_desc("rx_normal", "gridpoint timestep-1", "Normal Phase Speed for EW radiation OBCs", 'u', 'L')
    vd(2) = var_desc("ry_normal", "gridpoint timestep-1", "Normal Phase Speed for NS radiation OBCs", 'v', 'L')
    call register_restart_pair(OBC%rx_normal, OBC%ry_normal, vd(1), vd(2), .false., restart_CS, scalar_pair=.true.)
    ! The rx_normal and ry_normal arrays used with radiation OBCs are currently in units of grid
    ! points per timestep, but if this were to be corrected to [L T-1 ~> m s-1] or [T-1 ~> s-1] to
    ! permit timesteps to change between calls to the OBC code, the following would be needed instead:
    ! vd(1) = var_desc("rx_normal", "m s-1", "Normal Phase Speed for EW radiation OBCs", 'u', 'L')
    ! vd(2) = var_desc("ry_normal", "m s-1", "Normal Phase Speed for NS radiation OBCs", 'v', 'L')
    ! call register_restart_pair(OBC%rx_normal, OBC%ry_normal, vd(1), vd(2), .false., restart_CS, &
    !                            conversion=US%L_T_to_m_s, scalar_pair=.true.)
  endif

  if (OBC%oblique_BCs_exist_globally) then
    allocate(OBC%rx_oblique_u(HI%isdB:HI%iedB,HI%jsd:HI%jed,GV%ke), source=0.0)
    allocate(OBC%ry_oblique_u(HI%isdB:HI%iedB,HI%jsd:HI%jed,GV%ke), source=0.0)
    allocate(OBC%cff_normal_u(HI%IsdB:HI%IedB,HI%jsd:HI%jed,GV%ke), source=0.0)
    allocate(OBC%rx_oblique_v(HI%isd:HI%ied,HI%jsdB:HI%jedB,GV%ke), source=0.0)
    allocate(OBC%ry_oblique_v(HI%isd:HI%ied,HI%jsdB:HI%jedB,GV%ke), source=0.0)
    allocate(OBC%cff_normal_v(HI%isd:HI%ied,HI%jsdB:HI%jedB,GV%ke), source=0.0)

    vd(1) = var_desc("rx_oblique_u", "m2 s-2", "X-Direction Radiation Speed Squared for EW oblique OBCs", 'u', 'L')
    vd(2) = var_desc("ry_oblique_v", "m2 s-2", "Y-Direction Radiation Speed Squared for NS oblique OBCs", 'v', 'L')
    call register_restart_pair(OBC%rx_oblique_u, OBC%ry_oblique_v, vd(1), vd(2), .false., &
                               restart_CS, conversion=US%L_T_to_m_s**2)
    vd(1) = var_desc("ry_oblique_u", "m2 s-2", "Y-Direction Radiation Speed Squared for EW oblique OBCs", 'u', 'L')
    vd(2) = var_desc("rx_oblique_v", "m2 s-2", "X-Direction Radiation Speed Squared for NS oblique OBCs", 'v', 'L')
    call register_restart_pair(OBC%ry_oblique_u, OBC%rx_oblique_v, vd(1), vd(2), .false., &
                               restart_CS, conversion=US%L_T_to_m_s**2)

    vd(1) = var_desc("norm_oblique_u", "m2 s-2", "Denominator for normalizing EW oblique OBC radiation rates", &
                     'u', 'L')
    vd(2) = var_desc("norm_oblique_v", "m2 s-2", "Denominator for normalizing NS oblique OBC radiation rates", &
                     'v', 'L')
    call register_restart_pair(OBC%cff_normal_u, OBC%cff_normal_v, vd(1), vd(2), .false., &
                               restart_CS, conversion=US%L_T_to_m_s**2)
  endif

  if (Reg%ntr == 0) return
  if (.not. allocated(OBC%tracer_x_reservoirs_used)) then
    OBC%ntr = Reg%ntr
    allocate(OBC%tracer_x_reservoirs_used(Reg%ntr), source=.false.)
    allocate(OBC%tracer_y_reservoirs_used(Reg%ntr), source=.false.)
    call parse_for_tracer_reservoirs(OBC, param_file, use_temperature)
  else
    ! This would be coming from user code such as DOME.
    if (OBC%ntr /= Reg%ntr) then
!        call MOM_error(FATAL, "open_boundary_register_restarts: Inconsistent value for ntr")
      write(mesg,'("Inconsistent values for ntr ", I8," and ",I8,".")') OBC%ntr, Reg%ntr
      call MOM_error(WARNING, 'open_boundary_register_restarts: '//mesg)
    endif
  endif

  ! Still painfully inefficient, now in four dimensions.
  if (any(OBC%tracer_x_reservoirs_used)) then
    allocate(OBC%tres_x(HI%isdB:HI%iedB,HI%jsd:HI%jed,GV%ke,OBC%ntr), source=0.0)
    do m=1,OBC%ntr
      if (OBC%tracer_x_reservoirs_used(m)) then
        if (modulo(HI%turns, 2) /= 0) then
          write(var_name,'("tres_y_",I3.3)') m
          call register_restart_field(OBC%tres_x(:,:,:,m), var_name, .false., restart_CS, &
                   longname="Tracer concentration for NS OBCs", units="Conc", hor_grid='v')
        else
          write(var_name,'("tres_x_",I3.3)') m
          call register_restart_field(OBC%tres_x(:,:,:,m), var_name, .false., restart_CS, &
                   longname="Tracer concentration for EW OBCs", units="Conc", hor_grid='u')
        endif
      endif
    enddo
  endif
  if (any(OBC%tracer_y_reservoirs_used)) then
    allocate(OBC%tres_y(HI%isd:HI%ied,HI%jsdB:HI%jedB,GV%ke,OBC%ntr), source=0.0)
    do m=1,OBC%ntr
      if (OBC%tracer_y_reservoirs_used(m)) then
        if (modulo(HI%turns, 2) /= 0) then
          write(var_name,'("tres_x_",I3.3)') m
          call register_restart_field(OBC%tres_y(:,:,:,m), var_name, .false., restart_CS, &
                   longname="Tracer concentration for EW OBCs", units="Conc", hor_grid='u')
        else
          write(var_name,'("tres_y_",I3.3)') m
          call register_restart_field(OBC%tres_y(:,:,:,m), var_name, .false., restart_CS, &
                   longname="Tracer concentration for NS OBCs", units="Conc", hor_grid='v')
        endif
      endif
    enddo
  endif

end subroutine open_boundary_register_restarts

!> Update the OBC tracer reservoirs after the tracers have been updated.
subroutine update_segment_tracer_reservoirs(G, GV, uhr, vhr, h, OBC, dt, Reg)
  type(ocean_grid_type),                      intent(in) :: G   !< The ocean's grid structure
  type(verticalGrid_type),                    intent(in) :: GV  !<  Ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: uhr !< accumulated volume/mass flux through
                                                                !! the zonal face [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: vhr !< accumulated volume/mass flux through
                                                                !! the meridional face [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in) :: h   !< layer thickness after advection
                                                                !! [H ~> m or kg m-2]
  type(ocean_OBC_type),                       pointer    :: OBC !< Open boundary structure
  real,                                       intent(in) :: dt  !< time increment [T ~> s]
  type(tracer_registry_type),                 pointer    :: Reg !< pointer to tracer registry

  ! Local variable
  type(OBC_segment_type), pointer :: segment=>NULL()
  real :: u_L_in, u_L_out ! The zonal distance moved in or out of a cell, normalized by the reservoir
                          ! length scale [nondim]
  real :: v_L_in, v_L_out ! The meridional distance moved in or out of a cell, normalized by the reservoir
                          ! length scale [nondim]
  real :: fac1            ! The denominator of the expression for tracer updates [nondim]
  real :: I_scale         ! The inverse of the scaling factor for the tracers.
                          ! For salinity the units would be [ppt S-1 ~> 1]
  integer :: i, j, k, m, n, ntr, nz, ntr_id, fd_id
  integer :: ishift, idir, jshift, jdir
  real :: resrv_lfac_out  ! The reservoir inverse length scale scaling factor for the outward
                          ! direction per field [nondim]
  real :: resrv_lfac_in   ! The reservoir inverse length scale scaling factor for the inward
                          ! direction per field [nondim]
  real :: b_in, b_out     ! The 0 and 1 switch for tracer reservoirs
                          ! 1 if the length scale of reservoir is zero [nondim]
  real :: a_in, a_out     ! The 0 and 1(-1) switch for reservoir source weights
                          ! e.g. a_in is -1 only if b_in ==1 and uhr or vhr is inward
                          ! e.g. a_out is 1 only if b_out==1 and uhr or vhr is outward
                          ! It's clear that a_in and a_out cannot be both non-zero [nondim]
  nz = GV%ke
  ntr = Reg%ntr

  if (associated(OBC)) then ; if (OBC%OBC_pe) then ; do n=1,OBC%number_of_segments
    segment=>OBC%segment(n)
    if (.not. associated(segment%tr_Reg)) cycle
    b_in  = 0.0; if (segment%Tr_InvLscale_in  == 0.0) b_in  = 1.0
    b_out = 0.0; if (segment%Tr_InvLscale_out == 0.0) b_out = 1.0
    if (segment%is_E_or_W) then
      I = segment%HI%IsdB
      do j=segment%HI%jsd,segment%HI%jed
        ! ishift+I corresponds to the nearest interior tracer cell index
        ! idir switches the sign of the flow so that positive is into the reservoir
        if (segment%direction == OBC_DIRECTION_W) then
          ishift = 1 ; idir = -1
        else
          ishift = 0 ; idir = 1
        endif
        ! Can keep this or take it out, either way
        if (G%mask2dT(I+ishift,j) == 0.0) cycle
        ! Update the reservoir tracer concentration implicitly using a Backward-Euler timestep
        do m=1,segment%tr_Reg%ntseg
          ntr_id = segment%tr_Reg%Tr(m)%ntr_index
          fd_id = segment%tr_Reg%Tr(m)%fd_index
          if (fd_id == -1) then
            resrv_lfac_out = 1.0
            resrv_lfac_in  = 1.0
          else
            resrv_lfac_out = segment%field(fd_id)%resrv_lfac_out
            resrv_lfac_in  = segment%field(fd_id)%resrv_lfac_in
          endif
          I_scale = 1.0 ; if (segment%tr_Reg%Tr(m)%scale /= 0.0) I_scale = 1.0 / segment%tr_Reg%Tr(m)%scale
          if (allocated(segment%tr_Reg%Tr(m)%tres)) then ; do k=1,nz
            ! Calculate weights. Both a and u_L are nondim. Adding them together has no meaning.
            ! However, since they cannot be both non-zero, adding them works like a switch.
            ! When InvLscale_out is 0 and outflow, only interior data is applied to reservoirs
            ! When InvLscale_in is 0 and inflow, only nudged data is applied to reservoirs
            a_out = b_out * max(0.0, sign(1.0, idir*uhr(I,j,k)))
            a_in  = b_in  * min(0.0, sign(1.0, idir*uhr(I,j,k)))
            u_L_out = max(0.0, (idir*uhr(I,j,k))*segment%Tr_InvLscale_out*resrv_lfac_out / &
                      ((h(i+ishift,j,k) + GV%H_subroundoff)*G%dyCu(I,j)))
            u_L_in  = min(0.0, (idir*uhr(I,j,k))*segment%Tr_InvLscale_in*resrv_lfac_in  / &
                      ((h(i+ishift,j,k) + GV%H_subroundoff)*G%dyCu(I,j)))
            fac1 = (1.0 - (a_out - a_in)) + ((u_L_out + a_out) - (u_L_in + a_in))
            segment%tr_Reg%Tr(m)%tres(I,j,k) = (1.0/fac1) * &
                              ((1.0-a_out+a_in)*segment%tr_Reg%Tr(m)%tres(I,j,k)+ &
                              ((u_L_out+a_out)*Reg%Tr(ntr_id)%t(I+ishift,j,k) - &
                               (u_L_in+a_in)*segment%tr_Reg%Tr(m)%t(I,j,k)))
            if (allocated(OBC%tres_x)) OBC%tres_x(I,j,k,m) = I_scale * segment%tr_Reg%Tr(m)%tres(I,j,k)
          enddo ; endif
        enddo
      enddo
    elseif (segment%is_N_or_S) then
      J = segment%HI%JsdB
      do i=segment%HI%isd,segment%HI%ied
        ! jshift+J corresponds to the nearest interior tracer cell index
        ! jdir switches the sign of the flow so that positive is into the reservoir
        if (segment%direction == OBC_DIRECTION_S) then
          jshift = 1 ; jdir = -1
        else
          jshift = 0 ; jdir = 1
        endif
        ! Can keep this or take it out, either way
        if (G%mask2dT(i,j+jshift) == 0.0) cycle
        ! Update the reservoir tracer concentration implicitly using a Backward-Euler timestep
        do m=1,segment%tr_Reg%ntseg
          ntr_id = segment%tr_Reg%Tr(m)%ntr_index
          fd_id = segment%tr_Reg%Tr(m)%fd_index
          if (fd_id == -1) then
            resrv_lfac_out = 1.0
            resrv_lfac_in  = 1.0
          else
            resrv_lfac_out = segment%field(fd_id)%resrv_lfac_out
            resrv_lfac_in  = segment%field(fd_id)%resrv_lfac_in
          endif
          I_scale = 1.0 ; if (segment%tr_Reg%Tr(m)%scale /= 0.0) I_scale = 1.0 / segment%tr_Reg%Tr(m)%scale
          if (allocated(segment%tr_Reg%Tr(m)%tres)) then ; do k=1,nz
            a_out = b_out * max(0.0, sign(1.0, jdir*vhr(i,J,k)))
            a_in  = b_in  * min(0.0, sign(1.0, jdir*vhr(i,J,k)))
            v_L_out = max(0.0, (jdir*vhr(i,J,k))*segment%Tr_InvLscale_out*resrv_lfac_out / &
                      ((h(i,j+jshift,k) + GV%H_subroundoff)*G%dxCv(i,J)))
            v_L_in  = min(0.0, (jdir*vhr(i,J,k))*segment%Tr_InvLscale_in*resrv_lfac_in  / &
                      ((h(i,j+jshift,k) + GV%H_subroundoff)*G%dxCv(i,J)))
            fac1 = (1.0 - (a_out - a_in)) + ((v_L_out + a_out) - (v_L_in + a_in))
            segment%tr_Reg%Tr(m)%tres(i,J,k) = (1.0/fac1) * &
                              ((1.0-a_out+a_in)*segment%tr_Reg%Tr(m)%tres(i,J,k) + &
                              ((v_L_out+a_out)*Reg%Tr(ntr_id)%t(i,J+jshift,k) - &
                               (v_L_in+a_in)*segment%tr_Reg%Tr(m)%t(i,J,k)))
            if (allocated(OBC%tres_y)) OBC%tres_y(i,J,k,m) = I_scale * segment%tr_Reg%Tr(m)%tres(i,J,k)
          enddo ; endif
        enddo
      enddo
    endif
  enddo ; endif ; endif

end subroutine update_segment_tracer_reservoirs

!> Vertically remap the OBC tracer reservoirs and radiation rates that are filtered in time.
subroutine remap_OBC_fields(G, GV, h_old, h_new, OBC, PCM_cell)
  type(ocean_grid_type),                     intent(in) :: G     !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in) :: GV    !<  Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h_old !< Thickness of source grid [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h_new !< Thickness of destination grid [H ~> m or kg m-2]
  type(ocean_OBC_type),                      pointer    :: OBC   !< Open boundary structure
  logical, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                   optional, intent(in) :: PCM_cell !< Use PCM remapping in cells where true

  ! Local variables
  type(OBC_segment_type), pointer :: segment => NULL() ! A pointer to the various segments, used just for shorthand.

  real :: tr_column(GV%ke)  ! A column of updated tracer concentrations in internally scaled units.
                        ! For salinity the units would be [S ~> ppt].
  real :: r_norm_col(GV%ke) ! A column of updated radiation rates, in grid points per timestep [nondim]
  real :: rxy_col(GV%ke) ! A column of updated radiation rates for oblique OBCs [L2 T-2 ~> m2 s-2]
  real :: h1(GV%ke)     ! A column of source grid layer thicknesses [H ~> m or kg m-2]
  real :: h2(GV%ke)     ! A column of target grid layer thicknesses [H ~> m or kg m-2]
  real :: I_scale       ! The inverse of the scaling factor for the tracers.
                        ! For salinity the units would be [ppt S-1 ~> 1].
  logical :: PCM(GV%ke) ! If true, do PCM remapping from a cell.
  integer :: i, j, k, m, n, ntr, nz

  if (.not.associated(OBC)) return

  nz = GV%ke
  ntr = OBC%ntr

  if (.not.present(PCM_cell)) PCM(:) = .false.

  if (associated(OBC)) then ; if (OBC%OBC_pe) then ; do n=1,OBC%number_of_segments
    segment => OBC%segment(n)
    if (.not.associated(segment%tr_Reg)) cycle

    if (segment%is_E_or_W) then
      I = segment%HI%IsdB
      do j=segment%HI%jsd,segment%HI%jed

        ! Store a column of the start and final grids
        if (segment%direction == OBC_DIRECTION_W) then
          if (G%mask2dT(i+1,j) == 0.0) cycle
          h1(:) = h_old(i+1,j,:)
          h2(:) = h_new(i+1,j,:)
          if (present(PCM_cell)) then ; PCM(:) = PCM_cell(i+1,j,:) ; endif
        else
          if (G%mask2dT(i,j) == 0.0) cycle
          h1(:) = h_old(i,j,:)
          h2(:) = h_new(i,j,:)
          if (present(PCM_cell)) then ; PCM(:) = PCM_cell(i,j,:) ; endif
        endif

        ! Vertically remap the reservoir tracer concentrations
        do m=1,ntr ; if (allocated(segment%tr_Reg%Tr(m)%tres)) then
          I_scale = 1.0 ; if (segment%tr_Reg%Tr(m)%scale /= 0.0) I_scale = 1.0 / segment%tr_Reg%Tr(m)%scale

          if (present(PCM_cell)) then
            call remapping_core_h(OBC%remap_h_CS, nz, h1, segment%tr_Reg%Tr(m)%tres(I,j,:), nz, h2, tr_column, &
                                  PCM_cell=PCM)
          else
            call remapping_core_h(OBC%remap_h_CS, nz, h1, segment%tr_Reg%Tr(m)%tres(I,j,:), nz, h2, tr_column)
          endif

          ! Possibly underflow any very tiny tracer concentrations to 0?

          ! Update tracer concentrations
          segment%tr_Reg%Tr(m)%tres(I,j,:) = tr_column(:)
          if (allocated(OBC%tres_x)) then ; do k=1,nz
            OBC%tres_x(I,j,k,m) = I_scale * segment%tr_Reg%Tr(m)%tres(I,j,k)
          enddo ; endif

        endif ; enddo

        if (segment%radiation .and. (OBC%gamma_uv < 1.0)) then
          call remapping_core_h(OBC%remap_h_CS, nz, h1, segment%rx_norm_rad(I,j,:), nz, h2, r_norm_col, &
                                PCM_cell=PCM)

          do k=1,nz
            segment%rx_norm_rad(I,j,k) = r_norm_col(k)
            OBC%rx_normal(I,j,k) = segment%rx_norm_rad(I,j,k)
          enddo
        endif

        if (segment%oblique .and. (OBC%gamma_uv < 1.0)) then
          call remapping_core_h(OBC%remap_h_CS, nz, h1, segment%rx_norm_obl(I,j,:), nz, h2, rxy_col, &
                                PCM_cell=PCM)
          segment%rx_norm_obl(I,j,:) = rxy_col(:)
          call remapping_core_h(OBC%remap_h_CS, nz, h1, segment%ry_norm_obl(I,j,:), nz, h2, rxy_col, &
                                PCM_cell=PCM)
          segment%ry_norm_obl(I,j,:) = rxy_col(:)
          call remapping_core_h(OBC%remap_h_CS, nz, h1, segment%cff_normal(I,j,:), nz, h2, rxy_col, &
                                PCM_cell=PCM)
          segment%cff_normal(I,j,:) = rxy_col(:)

          do k=1,nz
            OBC%rx_oblique_u(I,j,k) = segment%rx_norm_obl(I,j,k)
            OBC%ry_oblique_u(I,j,k) = segment%ry_norm_obl(I,j,k)
            OBC%cff_normal_u(I,j,k) = segment%cff_normal(I,j,k)
          enddo
        endif

      enddo
    elseif (segment%is_N_or_S) then
      J = segment%HI%JsdB
      do i=segment%HI%isd,segment%HI%ied

        ! Store a column of the start and final grids
        if (segment%direction == OBC_DIRECTION_S) then
          if (G%mask2dT(i,j+1) == 0.0) cycle
          h1(:) = h_old(i,j+1,:)
          h2(:) = h_new(i,j+1,:)
          if (present(PCM_cell)) then ; PCM(:) = PCM_cell(i,j+1,:) ; endif
        else
          if (G%mask2dT(i,j) == 0.0) cycle
          h1(:) = h_old(i,j,:)
          h2(:) = h_new(i,j,:)
          if (present(PCM_cell)) then ; PCM(:) = PCM_cell(i,j,:) ; endif
        endif

        ! Vertically remap the reservoir tracer concentrations
        do m=1,ntr ; if (allocated(segment%tr_Reg%Tr(m)%tres)) then
          I_scale = 1.0 ; if (segment%tr_Reg%Tr(m)%scale /= 0.0) I_scale = 1.0 / segment%tr_Reg%Tr(m)%scale

          if (present(PCM_cell)) then
            call remapping_core_h(OBC%remap_h_CS, nz, h1, segment%tr_Reg%Tr(m)%tres(i,J,:), nz, h2, tr_column, &
                                  PCM_cell=PCM)
          else
            call remapping_core_h(OBC%remap_h_CS, nz, h1, segment%tr_Reg%Tr(m)%tres(i,J,:), nz, h2, tr_column)
          endif

          ! Possibly underflow any very tiny tracer concentrations to 0?

          ! Update tracer concentrations
          segment%tr_Reg%Tr(m)%tres(i,J,:) = tr_column(:)
          if (allocated(OBC%tres_y)) then ; do k=1,nz
            OBC%tres_y(i,J,k,m) = I_scale * segment%tr_Reg%Tr(m)%tres(i,J,k)
          enddo ; endif

        endif ; enddo

        if (segment%radiation .and. (OBC%gamma_uv < 1.0)) then
          call remapping_core_h(OBC%remap_h_CS, nz, h1, segment%ry_norm_rad(i,J,:), nz, h2, r_norm_col, &
                                PCM_cell=PCM)

          do k=1,nz
            segment%ry_norm_rad(i,J,k) = r_norm_col(k)
            OBC%ry_normal(i,J,k) = segment%ry_norm_rad(i,J,k)
          enddo
        endif

        if (segment%oblique .and. (OBC%gamma_uv < 1.0)) then
          call remapping_core_h(OBC%remap_h_CS, nz, h1, segment%rx_norm_obl(i,J,:), nz, h2, rxy_col, &
                                PCM_cell=PCM)
          segment%rx_norm_obl(i,J,:) = rxy_col(:)
          call remapping_core_h(OBC%remap_h_CS, nz, h1, segment%ry_norm_obl(i,J,:), nz, h2, rxy_col, &
                                PCM_cell=PCM)
          segment%ry_norm_obl(i,J,:) = rxy_col(:)
          call remapping_core_h(OBC%remap_h_CS, nz, h1, segment%cff_normal(i,J,:), nz, h2, rxy_col, &
                                PCM_cell=PCM)
          segment%cff_normal(i,J,:) = rxy_col(:)

          do k=1,nz
            OBC%rx_oblique_v(i,J,k) = segment%rx_norm_obl(i,J,k)
            OBC%ry_oblique_v(i,J,k) = segment%ry_norm_obl(i,J,k)
            OBC%cff_normal_v(i,J,k) = segment%cff_normal(i,J,k)
          enddo
        endif

      enddo
    endif
  enddo ; endif ; endif
  if (OBC%radiation_BCs_exist_globally) call pass_vector(OBC%rx_normal, OBC%ry_normal, G%Domain, &
                     To_All+Scalar_Pair)
  if (OBC%oblique_BCs_exist_globally) then
    call do_group_pass(OBC%pass_oblique, G%Domain)
  endif

end subroutine remap_OBC_fields


!> Adjust interface heights to fit the bathymetry and diagnose layer thickness.
!!
!! If the bottom most interface is below the topography then the bottom-most
!! layers are contracted to GV%Angstrom_Z.
!! If the bottom most interface is above the topography then the entire column
!! is dilated (expanded) to fill the void.
!!   @remark{There is a (hard-wired) "tolerance" parameter such that the
!! criteria for adjustment must equal or exceed 10cm.}
subroutine adjustSegmentEtaToFitBathymetry(G, GV, US, segment, fld, at_node)
  type(ocean_grid_type),   intent(in)    :: G   !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(OBC_segment_type),  intent(inout) :: segment !< OBC segment
  integer,                 intent(in)    :: fld  !< field index to adjust thickness
  logical,                 intent(in)    :: at_node !< True this point is at the OBC nodes rather than the faces

  integer :: i, j, k, is, ie, js, je, nz, contractions, dilations
  real, allocatable, dimension(:,:,:) :: eta ! Segment source data interface heights [Z ~> m]
  real, allocatable, dimension(:,:)   :: dz_tot ! Segment total thicknesses [Z ~> m]
  real :: hTolerance = 0.1 !<  Tolerance to exceed adjustment criteria [Z ~> m]
  ! real :: dilate      ! A factor by which to dilate the water column [nondim]
  !character(len=100) :: mesg

  hTolerance = 0.1*US%m_to_Z

  nz = size(segment%field(fld)%dz_src,3)

  if (segment%is_E_or_W) then
    is = segment%HI%IsdB ; ie = segment%HI%IedB
    if (at_node) then   ! This point is at the OBC nodes, rather than the cell face centers.
      Js = max(segment%Js_obc, G%jsd)
      Je = min(segment%Je_obc, G%jed-1)
    else   ! Segment thicknesses are defined at cell face centers.
      js = segment%HI%jsd ; je = segment%HI%jed
    endif
  else ! segment%is_N_or_S
    js = segment%HI%jsdB ; je = segment%HI%jedB
    if (at_node) then  ! This point is at the OBC nodes, rather than the cell face centers.
      is = max(segment%HI%IsdB, G%isd)
      ie = min(segment%HI%IedB, G%ied-1)
    else   ! Segment thicknesses are defined at cell face centers.
      is = segment%HI%isd ; ie = segment%HI%ied
    endif
  endif
  allocate(eta(is:ie,js:je,nz+1))
  allocate(dz_tot(is:ie,js:je), source=0.0)

  if (at_node) then
    if (segment%is_E_or_W) then
      I = Is
      do J=Js,Je
        dz_tot(I,J) = 0.5*(segment%dZtot(I,j) + segment%dZtot(I,j+1))
      enddo
      ! Do not extrapolate past the end of a global segment.
      ! ### For a concave corner between segments, perhaps we should do something more sophisticated.
      if (Js == segment%Js_obc) dz_tot(I,Js) = segment%dZtot(I,js+1)
      if (Je == segment%Js_obc) dz_tot(I,Je) = segment%dZtot(I,je)
    else
      J = Js
      do I=Is,Ie
        dz_tot(I,J) = 0.5*(segment%dZtot(i,J) + segment%dZtot(i+1,J))
      enddo
      ! Do not extrapolate past the end of a global segment.
      if (Is == segment%Is_obc) dz_tot(Is,J) = segment%dZtot(is+1,J)
      if (Ie == segment%Is_obc) dz_tot(Ie,J) = segment%dZtot(ie,J)
    endif
  else
    do j=js,je ; do i=is,ie
      dz_tot(i,j) = segment%dZtot(i,j)
    enddo ; enddo
  endif

  contractions = 0 ; dilations = 0
  do j=js,je ; do i=is,ie
    eta(i,j,1) = 0.0  ! segment data are assumed to be located on a static grid
    ! For remapping calls, the entire column will be dilated
    ! by a factor equal to the ratio of the sum of the geopotential referenced
    ! source data thicknesses, and the current model thicknesses. This could be
    ! an issue to be addressed, for instance if we are placing open boundaries
    ! under ice shelf cavities.
    do k=2,nz+1
      eta(i,j,k) = eta(i,j,k-1) - segment%field(fld)%dz_src(i,j,k-1)
    enddo
    ! The normal slope at the boundary is zero by a
    ! previous call to open_boundary_impose_normal_slope
    do k=nz+1,1,-1
      if (-eta(i,j,k) > dz_tot(i,j) + hTolerance) then
        eta(i,j,k) = -dz_tot(i,j)
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
    if (-eta(i,j,nz+1) < dz_tot(i,j) - hTolerance) then
      dilations = dilations + 1
      ! expand bottom-most cell only
      eta(i,j,nz+1) = -dz_tot(i,j)
      segment%field(fld)%dz_src(i,j,nz) = eta(i,j,nz) - eta(i,j,nz+1)
      ! if (eta(i,j,1) <= eta(i,j,nz+1)) then
      !   do k=1,nz ; segment%field(fld)%dz_src(i,j,k) = (eta(i,j,1) + G%bathyT(i,j)) / real(nz) ; enddo
      ! else
      !   dilate = (eta(i,j,1) + G%bathyT(i,j)) / (eta(i,j,1) - eta(i,j,nz+1))
      !   do k=1,nz ; segment%field(fld)%dz_src(i,j,k) = segment%field(fld)%dz_src(i,j,k) * dilate ; enddo
      ! endif
      !do k=nz,2,-1 ; eta(i,j,K) = eta(i,j,K+1) + segment%field(fld)%dz_src(i,j,k) ; enddo
    endif
  enddo ; enddo

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

  deallocate(eta, dz_tot)

end subroutine adjustSegmentEtaToFitBathymetry

!> This is more of a rotate initialization than an actual rotate
subroutine rotate_OBC_config(OBC_in, G_in, OBC, G, turns)
  type(ocean_OBC_type), pointer, intent(in)    :: OBC_in !< Input OBC
  type(dyn_horgrid_type),        intent(in)    :: G_in  !< Input grid
  type(ocean_OBC_type), pointer, intent(inout) :: OBC   !< Rotated OBC
  type(dyn_horgrid_type),        intent(in)    :: G     !< Rotated grid
  integer,                       intent(in)    :: turns !< Number of quarter turns

  integer :: c, n, l_seg

  if (OBC_in%number_of_segments==0) return

  ! Scalar and logical transfer
  OBC%number_of_segments = OBC_in%number_of_segments
  OBC%ke = OBC_in%ke
  OBC%user_BCs_set_globally = OBC_in%user_BCs_set_globally

  ! These are conditionally read and set if number_of_segments > 0
  OBC%zero_vorticity = OBC_in%zero_vorticity
  OBC%freeslip_vorticity = OBC_in%freeslip_vorticity
  OBC%computed_vorticity = OBC_in%computed_vorticity
  OBC%specified_vorticity = OBC_in%specified_vorticity
  OBC%zero_strain = OBC_in%zero_strain
  OBC%freeslip_strain = OBC_in%freeslip_strain
  OBC%computed_strain = OBC_in%computed_strain
  OBC%specified_strain = OBC_in%specified_strain
  OBC%zero_biharmonic = OBC_in%zero_biharmonic
  OBC%silly_h = OBC_in%silly_h
  OBC%silly_u = OBC_in%silly_u

  ! Segment rotation
  allocate(OBC%segment(0:OBC%number_of_segments))
  do l_seg = 1, OBC%number_of_segments
    call rotate_OBC_segment_config(OBC_in%segment(l_seg), G_in, OBC%segment(l_seg), G, turns)
    ! Data up to setup_[uv]_point_obc is needed for allocate_obc_segment_data!
    call allocate_OBC_segment_data(OBC, OBC%segment(l_seg))
  enddo

  ! The horizontal segment map
  allocate(OBC%segnum_u(G%IsdB:G%IedB,G%jsd:G%jed), source=0)
  allocate(OBC%segnum_v(G%isd:G%ied,G%JsdB:G%JedB), source=0)
  call rotate_array_pair(OBC_in%segnum_u, OBC_in%segnum_v, turns, OBC%segnum_u, OBC%segnum_v)
  call set_segnum_signs(OBC, G)

  ! These are conditionally enabled during segment configuration
  if (modulo(turns,2) == 0) then
    OBC%open_u_BCs_exist_globally = OBC_in%open_u_BCs_exist_globally
    OBC%open_v_BCs_exist_globally = OBC_in%open_v_BCs_exist_globally
    OBC%Flather_u_BCs_exist_globally = OBC_in%Flather_u_BCs_exist_globally
    OBC%Flather_v_BCs_exist_globally = OBC_in%Flather_v_BCs_exist_globally
    OBC%nudged_u_BCs_exist_globally = OBC_in%nudged_u_BCs_exist_globally
    OBC%nudged_v_BCs_exist_globally = OBC_in%nudged_v_BCs_exist_globally
    OBC%specified_u_BCs_exist_globally= OBC_in%specified_u_BCs_exist_globally
    OBC%specified_v_BCs_exist_globally= OBC_in%specified_v_BCs_exist_globally
  else  ! Swap information for u- and v- OBCs
    OBC%open_u_BCs_exist_globally = OBC_in%open_v_BCs_exist_globally
    OBC%open_v_BCs_exist_globally = OBC_in%open_u_BCs_exist_globally
    OBC%Flather_u_BCs_exist_globally = OBC_in%Flather_v_BCs_exist_globally
    OBC%Flather_v_BCs_exist_globally = OBC_in%Flather_u_BCs_exist_globally
    OBC%nudged_u_BCs_exist_globally = OBC_in%nudged_v_BCs_exist_globally
    OBC%nudged_v_BCs_exist_globally = OBC_in%nudged_u_BCs_exist_globally
    OBC%specified_u_BCs_exist_globally= OBC_in%specified_v_BCs_exist_globally
    OBC%specified_v_BCs_exist_globally= OBC_in%specified_u_BCs_exist_globally
  endif
  OBC%oblique_BCs_exist_globally = OBC_in%oblique_BCs_exist_globally
  OBC%radiation_BCs_exist_globally = OBC_in%radiation_BCs_exist_globally

  ! These are set by initialize_segment_data
  OBC%brushcutter_mode = OBC_in%brushcutter_mode
  OBC%update_OBC = OBC_in%update_OBC
  OBC%needs_IO_for_data = OBC_in%needs_IO_for_data
  OBC%any_needs_IO_for_data = OBC_in%any_needs_IO_for_data

  OBC%update_OBC_seg_data = OBC_in%update_OBC_seg_data
  OBC%ntr = OBC_in%ntr
  if (OBC%ntr > 0) then
    allocate(OBC%tracer_x_reservoirs_used(OBC%ntr), source=.false.)
    allocate(OBC%tracer_y_reservoirs_used(OBC%ntr), source=.false.)
    if (modulo(turns,2) == 0) then
      do n=1,OBC%ntr
        OBC%tracer_x_reservoirs_used(n) = OBC_in%tracer_x_reservoirs_used(n)
        OBC%tracer_y_reservoirs_used(n) = OBC_in%tracer_y_reservoirs_used(n)
      enddo
    else  ! Swap information for u- and v- OBCs
      do n=1,OBC%ntr
        OBC%tracer_x_reservoirs_used(n) = OBC_in%tracer_y_reservoirs_used(n)
        OBC%tracer_y_reservoirs_used(n) = OBC_in%tracer_x_reservoirs_used(n)
      enddo
    endif
  endif

  OBC%gamma_uv = OBC_in%gamma_uv
  OBC%rx_max = OBC_in%rx_max
  OBC%OBC_pe = OBC_in%OBC_pe

  ! These are run-time parameters that are read in via open_boundary_config
  OBC%debug = OBC_in%debug
  OBC%ramp = OBC_in%ramp
  OBC%ramping_is_activated = OBC_in%ramping_is_activated
  OBC%ramp_timescale = OBC_in%ramp_timescale
  OBC%trunc_ramp_time = OBC_in%trunc_ramp_time
  OBC%ramp_value = OBC_in%ramp_value
  OBC%ramp_start_time = OBC_in%ramp_start_time
  OBC%remap_answer_date = OBC_in%remap_answer_date
  OBC%check_reconstruction  = OBC_in%check_reconstruction
  OBC%check_remapping = OBC_in%check_remapping
  OBC%force_bounds_in_subcell = OBC_in%force_bounds_in_subcell
  OBC%om4_remap_via_sub_cells = OBC_in%om4_remap_via_sub_cells
  OBC%remappingScheme = OBC_in%remappingScheme
  OBC%exterior_OBC_bug = OBC_in%exterior_OBC_bug
  OBC%hor_index_bug = OBC_in%hor_index_bug
  OBC%n_tide_constituents = OBC_in%n_tide_constituents
  OBC%add_tide_constituents = OBC_in%add_tide_constituents

  ! These are read in via initialize_obc_tides when n_tide_constituents > 0
  if (OBC%add_tide_constituents .and. (OBC%n_tide_constituents>0)) then
    OBC%add_eq_phase = OBC_in%add_eq_phase
    OBC%add_nodal_terms = OBC_in%add_nodal_terms
    OBC%time_ref = OBC_in%time_ref

    allocate(OBC%tide_names(OBC%n_tide_constituents))
    allocate(OBC%tide_frequencies(OBC%n_tide_constituents))
    allocate(OBC%tide_eq_phases(OBC%n_tide_constituents))
    allocate(OBC%tide_fn(OBC%n_tide_constituents))
    allocate(OBC%tide_un(OBC%n_tide_constituents))
    do c=1,OBC%n_tide_constituents
      OBC%tide_names(c) = OBC_in%tide_names(c)
      OBC%tide_frequencies(c) = OBC_in%tide_frequencies(c)
      OBC%tide_eq_phases(c) = OBC_in%tide_eq_phases(c)
      OBC%tide_fn(c) = OBC_in%tide_fn(c)
      OBC%tide_un(c) = OBC_in%tide_un(c)
    enddo

    if (OBC%add_eq_phase .or.  OBC%add_nodal_terms)  &
      OBC%tidal_longitudes = OBC_in%tidal_longitudes
  endif

  ! remap_z_CS and remap_h_CS are set up by initialize_segment_data, so we copy the fields here.
  if (ASSOCIATED(OBC_in%remap_z_CS)) then
    allocate(OBC%remap_z_CS)
    OBC%remap_z_CS = OBC_in%remap_z_CS
  endif
  if (ASSOCIATED(OBC_in%remap_h_CS)) then
    allocate(OBC%remap_h_CS)
    OBC%remap_h_CS = OBC_in%remap_h_CS
  endif

  ! TODO: The OBC registry is a list of "registered" OBC types that is set up
  !   on the rotated grid after rotate_OBC_config is called.
  !   It does not appear to be used, so for now we skip this record.
  ! OBC%OBC_Reg => OBC_in%OBC_Reg
  ! OBC%num_obgc_tracers = OBC_in%num_obgc_tracers
  ! if (associated(OBC_in%OBC_Reg) .and. (.not.associated(OBC%OBC_Reg))) allocate(OBC%OBC_Reg)
end subroutine rotate_OBC_config

!> Rotate the OBC segment configuration data from the input to model index map.
subroutine rotate_OBC_segment_config(segment_in, G_in, segment, G, turns)
  type(OBC_segment_type), intent(in) :: segment_in  !< Input OBC segment
  type(dyn_horgrid_type),  intent(in) :: G_in       !< Input grid metric
  type(OBC_segment_type), intent(inout) :: segment  !< Rotated OBC segment
  type(dyn_horgrid_type),  intent(in) :: G          !< Rotated grid metric
  integer, intent(in) :: turns                      !< Number of quarter turns

  ! Global segment indices
  integer :: Is_obc_in, Ie_obc_in, Js_obc_in, Je_obc_in ! Input domain global indices
  integer :: Is_obc, Ie_obc, Js_obc, Je_obc             ! Rotated domain global indices
  integer :: qturns ! The number of quarter turns in the range of 0 to 3

  ! NOTE: A "rotation" of the OBC segment string would allow us to use
  !   setup_[uv]_point_obc to set up most of this.  For now, we just copy/swap
  !   flags and manually rotate the indices.

  ! This is set if the segment is in the local grid
  segment%on_pe = segment_in%on_pe

  qturns = modulo(turns, 4)

  ! Transfer configuration flags
  segment%Flather = segment_in%Flather
  segment%radiation = segment_in%radiation
  segment%radiation_tan = segment_in%radiation_tan
  segment%radiation_grad = segment_in%radiation_grad
  segment%oblique = segment_in%oblique
  segment%oblique_tan = segment_in%oblique_tan
  segment%oblique_grad = segment_in%oblique_grad
  segment%nudged = segment_in%nudged
  segment%nudged_tan = segment_in%nudged_tan
  segment%nudged_grad = segment_in%nudged_grad
  segment%specified = segment_in%specified
  segment%specified_tan = segment_in%specified_tan
  segment%specified_grad = segment_in%specified_grad
  segment%open = segment_in%open
  segment%gradient = segment_in%gradient

  if ((qturns == 0) .or. (qturns == 2)) then
    segment%u_values_needed = segment_in%u_values_needed
    segment%v_values_needed = segment_in%v_values_needed
    segment%uamp_values_needed = segment_in%uamp_values_needed
    segment%vamp_values_needed = segment_in%vamp_values_needed
    segment%uphase_values_needed = segment_in%uphase_values_needed
    segment%vphase_values_needed = segment_in%vphase_values_needed
    segment%uamp_index = segment_in%uamp_index  ! ### Perhaps this should not be set here.
    segment%vamp_index = segment_in%vamp_index
    segment%uphase_index = segment_in%uphase_index
    segment%vphase_index = segment_in%vphase_index
  else ! NOTE: [uv]_values_needed are swapped
    segment%u_values_needed = segment_in%v_values_needed
    segment%v_values_needed = segment_in%u_values_needed
    segment%uamp_values_needed = segment_in%vamp_values_needed
    segment%vamp_values_needed = segment_in%uamp_values_needed
    segment%uphase_values_needed = segment_in%vphase_values_needed
    segment%vphase_values_needed = segment_in%uphase_values_needed
    segment%uamp_index = segment_in%vamp_index  ! ### Perhaps this should not be set here.
    segment%vamp_index = segment_in%uamp_index
    segment%uphase_index = segment_in%vphase_index
    segment%vphase_index = segment_in%uphase_index
  endif
  segment%z_values_needed = segment_in%z_values_needed
  segment%g_values_needed = segment_in%g_values_needed
  segment%t_values_needed = segment_in%t_values_needed
  segment%s_values_needed = segment_in%s_values_needed
  segment%zamp_values_needed = segment_in%zamp_values_needed
  segment%zphase_values_needed = segment_in%zphase_values_needed
  segment%zamp_index = segment_in%zamp_index  ! ### Perhaps this should not be set here.
  segment%zphase_index = segment_in%zphase_index

  segment%values_needed = segment_in%values_needed

  ! These are conditionally set if nudged
  segment%Velocity_nudging_timescale_in = segment_in%Velocity_nudging_timescale_in
  segment%Velocity_nudging_timescale_out= segment_in%Velocity_nudging_timescale_out

  ! Rotate segment indices

  ! Reverse engineer the input [IJ][se]_obc segment indices
  ! NOTE: The values stored in the segment are always saved in ascending order,
  !   e.g. (is < ie).  In order to use setup_segment_indices, we reorder the
  !   indices here to indicate face direction.
  !   Segment indices are also indexed locally, so here we convert to global indices
  if (segment_in%direction == OBC_DIRECTION_N) then
    Is_obc_in = segment_in%Ie_obc + G_in%idg_offset
    Ie_obc_in = segment_in%Is_obc + G_in%idg_offset
  else
    Is_obc_in = segment_in%Is_obc + G_in%idg_offset
    Ie_obc_in = segment_in%Ie_obc + G_in%idg_offset
  endif

  if (segment_in%direction == OBC_DIRECTION_W) then
    Js_obc_in = segment_in%Je_obc + G_in%jdg_offset
    Je_obc_in = segment_in%Js_obc + G_in%jdg_offset
  else
    Js_obc_in = segment_in%Js_obc + G_in%jdg_offset
    Je_obc_in = segment_in%Je_obc + G_in%jdg_offset
  endif

  ! Rotate the global indices of the segment according to the number of turns.
  if (qturns == 0) then
    Is_obc = Is_obc_in ; Ie_obc = Ie_obc_in
    Js_obc = Js_obc_in ; Je_obc = Je_obc_in
  elseif (qturns == 1) then
    Is_obc = G_in%JegB - Js_obc_in ; Ie_obc = G_in%JegB - Je_obc_in
    Js_obc = Is_obc_in ; Je_obc = Ie_obc_in
  elseif (qturns == 2) then
    Is_obc = G_in%IegB - Is_obc_in ; Ie_obc = G_in%IegB - Ie_obc_in
    Js_obc = G_in%JegB - Js_obc_in ; Je_obc = G_in%JegB - Je_obc_in
  elseif (qturns == 3) then
    Is_obc = Js_obc_in ; Ie_obc = Je_obc_in
    Js_obc = G_in%IegB - Is_obc_in ; Je_obc = G_in%IegB - Ie_obc_in
  endif

  ! Orientation is based on the index ordering, and setup_segment_indices
  ! is based on the the original order in the intput files.
  call setup_segment_indices(G, segment, Is_obc, Ie_obc, Js_obc, Je_obc)

  ! Re-order [IJ][se]_obc back to ascending, and remove the global indexing offset.
  if (Is_obc > Ie_obc) then
    segment%Is_obc = Ie_obc - G%idg_offset
    segment%Ie_obc = Is_obc - G%idg_offset
  else
    segment%Is_obc = Is_obc - G%idg_offset
    segment%Ie_obc = Ie_obc - G%idg_offset
  endif

  if (Js_obc > Je_obc) then
    segment%Js_obc = Je_obc - G%jdg_offset
    segment%Je_obc = Js_obc - G%jdg_offset
  else
    segment%Js_obc = Js_obc - G%jdg_offset
    segment%Je_obc = Je_obc - G%jdg_offset
  endif

  ! Reconfigure the directional flags
  segment%direction = rotate_OBC_segment_direction(segment_in%direction, turns)

  segment%is_E_or_W_2 = ((segment%direction == OBC_DIRECTION_E) .or. &
                         (segment%direction == OBC_DIRECTION_W))
  segment%is_E_or_W = segment_in%on_PE .and. segment%is_E_or_W_2
  segment%is_N_or_S = segment_in%on_PE .and. &
                        ((segment%direction == OBC_DIRECTION_N) .or. &
                         (segment%direction == OBC_DIRECTION_S))

  ! These are conditionally set if Lscale_{in,out} are present
  segment%Tr_InvLscale_in = segment_in%Tr_InvLscale_in
  segment%Tr_InvLscale_out = segment_in%Tr_InvLscale_out
end subroutine rotate_OBC_segment_config


!> Return the direction of an OBC segment on after rotation to the new grid.  Note that
!! rotate_OBC_seg_direction(rotate_OBC_seg_direction(direction, turns), -turns) = direction.
function rotate_OBC_segment_direction(direction, turns) result(rotated_dir)
  integer, intent(in) :: direction  !< The orientation of an OBC segment on the original grid
  integer, intent(in) :: turns      !< Number of quarter turns
  integer :: rotated_dir  !< An integer encoding the new rotated segment direction

  integer :: qturns ! The number of quarter turns in the range of 0 to 3

  qturns = modulo(turns, 4)

  if ((qturns == 0) .or. (direction == OBC_NONE)) then
    rotated_dir = direction
  else  ! Determine the segment direction on a rotated grid
    select case (direction)
      case (OBC_DIRECTION_N)
        if (qturns == 0) rotated_dir = OBC_DIRECTION_N
        if (qturns == 1) rotated_dir = OBC_DIRECTION_W
        if (qturns == 2) rotated_dir = OBC_DIRECTION_S
        if (qturns == 3) rotated_dir = OBC_DIRECTION_E
      case (OBC_DIRECTION_W)
        if (qturns == 0) rotated_dir = OBC_DIRECTION_W
        if (qturns == 1) rotated_dir = OBC_DIRECTION_S
        if (qturns == 2) rotated_dir = OBC_DIRECTION_E
        if (qturns == 3) rotated_dir = OBC_DIRECTION_N
      case (OBC_DIRECTION_S)
        if (qturns == 0) rotated_dir = OBC_DIRECTION_S
        if (qturns == 1) rotated_dir = OBC_DIRECTION_E
        if (qturns == 2) rotated_dir = OBC_DIRECTION_N
        if (qturns == 3) rotated_dir = OBC_DIRECTION_W
      case (OBC_DIRECTION_E)
        if (qturns == 0) rotated_dir = OBC_DIRECTION_E
        if (qturns == 1) rotated_dir = OBC_DIRECTION_N
        if (qturns == 2) rotated_dir = OBC_DIRECTION_W
        if (qturns == 3) rotated_dir = OBC_DIRECTION_S
      case (OBC_NONE)
        rotated_dir = OBC_NONE
      case default ! This should never happen.
        rotated_dir = direction
    end select
  endif

end function rotate_OBC_segment_direction


!> Initialize the segments and field-related data of a rotated OBC.
subroutine rotate_OBC_segment_fields(OBC_in, G, OBC)
  type(ocean_OBC_type), intent(in) :: OBC_in            !< OBC on input map
  type(ocean_grid_type), intent(in) :: G                !< Rotated grid metric
  type(ocean_OBC_type), pointer, intent(inout) :: OBC   !< Rotated OBC

  integer :: l_seg

  ! update_OBC may have been updated during initialization.
  OBC%update_OBC = OBC_in%update_OBC

  do l_seg = 1, OBC%number_of_segments
    call rotate_OBC_segment_data(OBC_in%segment(l_seg), OBC%segment(l_seg), G%HI%turns)
  enddo

end subroutine rotate_OBC_segment_fields


!> Initialize the segment tracer registry of a rotated OBC.
subroutine rotate_OBC_segment_tracer_registry(OBC_in, G, OBC)
  type(ocean_OBC_type), intent(in) :: OBC_in            !< OBC on input map
  type(ocean_grid_type), intent(in) :: G                !< Rotated grid metric
  type(ocean_OBC_type), pointer, intent(inout) :: OBC   !< Rotated OBC

  integer :: l_seg

  do l_seg = 1, OBC%number_of_segments
    call rotate_OBC_segment_tracer_data(OBC_in%segment(l_seg), OBC%segment(l_seg), G%HI%turns)
  enddo

end subroutine rotate_OBC_segment_tracer_registry


!> Rotate an OBC segment's fields from the input to the model index map.
subroutine rotate_OBC_segment_data(segment_in, segment, turns)
  type(OBC_segment_type), intent(in) :: segment_in  !< The unrotated segment to use as a source
  type(OBC_segment_type), intent(inout) :: segment  !< The rotated segment to initialize
  integer, intent(in) :: turns  !< The number of quarter turns of the grid to apply

  ! Local variables
  logical :: flip_normal_vel_sign, flip_tang_vel_sign
  integer :: n
  integer :: num_fields
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, ke

  num_fields = segment_in%num_fields
  allocate(segment%field(num_fields))

  isd = segment%HI%isd ; ied = segment%HI%ied
  jsd = segment%HI%jsd ; jed = segment%HI%jed
  IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
  JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB

  if ((turns == 0) .or. (turns == 2)) then
    segment%uamp_index = segment_in%uamp_index
    segment%vamp_index = segment_in%vamp_index
    segment%uphase_index = segment_in%uphase_index
    segment%vphase_index = segment_in%vphase_index
  else ! NOTE: [uv]_values_needed are swapped
    segment%uamp_index = segment_in%vamp_index
    segment%vamp_index = segment_in%uamp_index
    segment%uphase_index = segment_in%vphase_index
    segment%vphase_index = segment_in%uphase_index
  endif
  segment%zamp_index = segment_in%zamp_index
  segment%zphase_index = segment_in%zphase_index

  segment%num_fields = segment_in%num_fields
  do n = 1, num_fields
    segment%field(n)%handle = segment_in%field(n)%handle
    segment%field(n)%dz_handle = segment_in%field(n)%dz_handle
    segment%field(n)%use_IO = segment_in%field(n)%use_IO
    segment%field(n)%genre = segment_in%field(n)%genre
    segment%field(n)%scale = segment_in%field(n)%scale
    segment%field(n)%resrv_lfac_in = segment_in%field(n)%resrv_lfac_in
    segment%field(n)%resrv_lfac_out = segment_in%field(n)%resrv_lfac_out
    segment%field(n)%on_face = segment_in%field(n)%on_face

    if (allocated(segment_in%field(n)%buffer_dst)) then
      call allocate_rotated_seg_data(segment_in%field(n)%buffer_dst, segment_in%HI, &
                                     segment%field(n)%buffer_dst, segment)
      call rotate_array(segment_in%field(n)%buffer_dst, turns, segment%field(n)%buffer_dst)
    endif

    if (modulo(turns, 2) /= 0) then
      select case (segment_in%field(n)%name)
        case ('U')
          segment%field(n)%name = 'V'
        case ('Uamp')
          segment%field(n)%name = 'Vamp'
        case ('Uphase')
          segment%field(n)%name = 'Vphase'
        case ('V')
          segment%field(n)%name = 'U'
        case ('Vamp')
          segment%field(n)%name = 'Uamp'
        case ('Vphase')
          segment%field(n)%name = 'Uphase'
        case ('DVDX')
          segment%field(n)%name = 'DUDY'
        case ('DUDY')
          segment%field(n)%name = 'DVDX'
        case default
          segment%field(n)%name = segment_in%field(n)%name
      end select
    else
      segment%field(n)%name = segment_in%field(n)%name
    endif

    if (allocated(segment_in%field(n)%buffer_src)) then
      call allocate_rotated_seg_data(segment_in%field(n)%buffer_src, segment_in%HI, &
                                     segment%field(n)%buffer_src, segment)
      call rotate_array(segment_in%field(n)%buffer_src, turns, segment%field(n)%buffer_src)
    endif

    segment%field(n)%nk_src = segment_in%field(n)%nk_src

    if (allocated(segment_in%field(n)%dz_src)) then
      call allocate_rotated_seg_data(segment_in%field(n)%dz_src, segment_in%HI, segment%field(n)%dz_src, segment)
      call rotate_array(segment_in%field(n)%dz_src, turns, segment%field(n)%dz_src)
    endif

    segment%field(n)%value = segment_in%field(n)%value
  enddo

  if (allocated(segment_in%SSH)) &
      call rotate_array(segment_in%SSH, turns, segment%SSH)
  if (allocated(segment_in%cg)) &
      call rotate_array(segment_in%cg, turns, segment%cg)
  if (allocated(segment_in%htot)) &
      call rotate_array(segment_in%htot, turns, segment%htot)
  if (allocated(segment_in%dztot)) &
      call rotate_array(segment_in%dztot, turns, segment%dztot)
  if (allocated(segment_in%h)) &
      call rotate_array(segment_in%h, turns, segment%h)
  if (allocated(segment_in%normal_vel)) &
      call rotate_array(segment_in%normal_vel, turns, segment%normal_vel)
  if (allocated(segment_in%normal_trans)) &
      call rotate_array(segment_in%normal_trans, turns, segment%normal_trans)
  if (allocated(segment_in%normal_vel_bt)) &
      call rotate_array(segment_in%normal_vel_bt, turns, segment%normal_vel_bt)
  if (allocated(segment_in%tangential_vel)) &
      call rotate_array(segment_in%tangential_vel, turns, segment%tangential_vel)
  if (allocated(segment_in%tangential_grad)) &
      call rotate_array(segment_in%tangential_grad, turns, segment%tangential_grad)
  if (allocated(segment_in%grad_normal)) &
      call rotate_array(segment_in%grad_normal, turns, segment%grad_normal)
  if (allocated(segment_in%grad_tan)) &
      call rotate_array(segment_in%grad_tan, turns, segment%grad_tan)
  if (allocated(segment_in%grad_gradient)) &
      call rotate_array(segment_in%grad_gradient, turns, segment%grad_gradient)
  if (modulo(turns, 2) /= 0) then
    if (allocated(segment_in%rx_norm_rad)) &
        call rotate_array(segment_in%rx_norm_rad, turns, segment%ry_norm_rad)
    if (allocated(segment_in%ry_norm_rad)) &
        call rotate_array(segment_in%ry_norm_rad, turns, segment%rx_norm_rad)
    if (allocated(segment_in%rx_norm_obl)) &
        call rotate_array(segment_in%rx_norm_obl, turns, segment%ry_norm_obl)
    if (allocated(segment_in%ry_norm_obl)) &
        call rotate_array(segment_in%ry_norm_obl, turns, segment%rx_norm_obl)
  else
    if (allocated(segment_in%rx_norm_rad)) &
        call rotate_array(segment_in%rx_norm_rad, turns, segment%rx_norm_rad)
    if (allocated(segment_in%ry_norm_rad)) &
        call rotate_array(segment_in%ry_norm_rad, turns, segment%ry_norm_rad)
    if (allocated(segment_in%rx_norm_obl)) &
        call rotate_array(segment_in%rx_norm_obl, turns, segment%rx_norm_obl)
    if (allocated(segment_in%ry_norm_obl)) &
        call rotate_array(segment_in%ry_norm_obl, turns, segment%ry_norm_obl)
  endif
  if (allocated(segment_in%cff_normal)) &
      call rotate_array(segment_in%cff_normal, turns, segment%cff_normal)
  if (allocated(segment_in%nudged_normal_vel)) &
      call rotate_array(segment_in%nudged_normal_vel, turns, segment%nudged_normal_vel)
  if (allocated(segment_in%nudged_tangential_vel)) &
      call rotate_array(segment_in%nudged_tangential_vel, turns, segment%nudged_tangential_vel)
  if (allocated(segment_in%nudged_tangential_grad)) &
      call rotate_array(segment_in%nudged_tangential_grad, turns, segment%nudged_tangential_grad)

  ! Change the sign of the normal or tangential velocities or transports that have been read in from
  ! a file, depending on the orientation of the face and the number of quarter turns of the grid.
  flip_normal_vel_sign = .false. ; flip_tang_vel_sign = .false.
  do n = 1, num_fields
    if (((segment%field(n)%name == 'U') .or. (segment%field(n)%name == 'Uamp')) .and. &
        ((modulo(turns, 4) == 1) .or. (modulo(turns, 4) == 2)) ) then
      if (allocated(segment%field(n)%buffer_dst)) &
        segment%field(n)%buffer_dst(:,:,:) = -segment%field(n)%buffer_dst(:,:,:)
      segment%field(n)%value = -segment%field(n)%value
      if (segment%is_E_or_W) flip_normal_vel_sign = .true.
      if (segment%is_N_or_S) flip_tang_vel_sign = .true.
    elseif (((segment%field(n)%name == 'V') .or. (segment%field(n)%name == 'Vamp')) .and. &
            ((modulo(turns, 4) == 3) .or. (modulo(turns, 4) == 2)) ) then
      if (allocated(segment%field(n)%buffer_dst)) &
        segment%field(n)%buffer_dst(:,:,:) = -segment%field(n)%buffer_dst(:,:,:)
      segment%field(n)%value = -segment%field(n)%value
      if (segment%is_N_or_S) flip_normal_vel_sign = .true.
      if (segment%is_E_or_W) flip_tang_vel_sign = .true.
    endif
  enddo

  if (flip_normal_vel_sign) then
    segment%normal_trans(:,:,:) = -segment%normal_trans(:,:,:)
    segment%normal_vel(:,:,:) = -segment%normal_vel(:,:,:)
    segment%normal_vel_bt(:,:) = -segment%normal_vel_bt(:,:)
    if (allocated(segment%nudged_normal_vel)) &
      segment%nudged_normal_vel(:,:,:) = -segment%nudged_normal_vel(:,:,:)
  endif

  if (flip_tang_vel_sign) then
    if (allocated(segment%tangential_vel)) &
      segment%tangential_vel(:,:,:) = -segment%tangential_vel(:,:,:)
    if (allocated(segment%nudged_tangential_vel)) &
      segment%nudged_tangential_vel(:,:,:) = -segment%nudged_tangential_vel(:,:,:)
  endif

  segment%temp_segment_data_exists = segment_in%temp_segment_data_exists
  segment%salt_segment_data_exists = segment_in%salt_segment_data_exists
end subroutine rotate_OBC_segment_data

!> Rotate an OBC segment's tracer registry fields fields from the input to the model index map.
subroutine rotate_OBC_segment_tracer_data(segment_in, segment, turns)
  type(OBC_segment_type), intent(in) :: segment_in  !< The unrotated segment to use as a source
  type(OBC_segment_type), intent(inout) :: segment  !< The rotated segment to initialize
  integer, intent(in) :: turns  !< The number of quarter turns of the grid to apply

  ! Local variables
  integer :: n
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, ke

  if (associated(segment_in%tr_Reg)) then
    isd = segment%HI%isd ; ied = segment%HI%ied ; IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
    jsd = segment%HI%jsd ; jed = segment%HI%jed ; JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB

    if (.not.associated(segment%tr_Reg)) allocate(segment%tr_Reg)
    segment%tr_Reg%ntseg = segment_in%tr_Reg%ntseg

    do n = 1, segment_in%tr_Reg%ntseg
      ! segment_in already points to the rotated tracer fields in the registry.
      if (associated(segment_in%tr_Reg%Tr(n)%Tr)) &
        segment%tr_Reg%Tr(n)%Tr => segment_in%tr_Reg%Tr(n)%Tr
      segment%tr_Reg%Tr(n)%name = segment_in%tr_Reg%Tr(n)%name
      segment%tr_Reg%Tr(n)%ntr_index = segment_in%tr_Reg%Tr(n)%ntr_index
      segment%tr_Reg%Tr(n)%fd_index = segment_in%tr_Reg%Tr(n)%fd_index
      segment%tr_Reg%Tr(n)%scale = segment_in%tr_Reg%Tr(n)%scale
      segment%tr_Reg%Tr(n)%OBC_inflow_conc = segment_in%tr_Reg%Tr(n)%OBC_inflow_conc

      if (allocated(segment_in%tr_Reg%tr(n)%t)) then
        if (.not.allocated(segment%tr_Reg%tr(n)%t)) then
          ke = size(segment_in%tr_Reg%tr(n)%t, 3)
          if (segment%is_E_or_W) then
            allocate(segment%tr_Reg%Tr(n)%t(IsdB:IedB,jsd:jed,1:ke), source=0.0)
          elseif (segment%is_N_or_S) then
            allocate(segment%tr_Reg%Tr(n)%t(isd:ied,JsdB:JedB,1:ke), source=0.0)
          endif
        endif
        call rotate_array(segment_in%tr_Reg%tr(n)%t, turns, segment%tr_Reg%tr(n)%t)
      endif

      if (allocated(segment_in%tr_Reg%tr(n)%tres)) then
        if (.not.allocated(segment%tr_Reg%tr(n)%tres)) then
          ke = size(segment_in%tr_Reg%tr(n)%tres, 3)
          if (segment%is_E_or_W) then
            allocate(segment%tr_Reg%Tr(n)%tres(IsdB:IedB,jsd:jed,1:ke), source=0.0)
          elseif (segment%is_N_or_S) then
            allocate(segment%tr_Reg%Tr(n)%tres(isd:ied,JsdB:JedB,1:ke), source=0.0)
          endif
        endif
        call rotate_array(segment_in%tr_Reg%tr(n)%tres, turns, segment%tr_Reg%tr(n)%tres)
      endif
      segment%tr_Reg%Tr(n)%is_initialized = segment_in%tr_Reg%Tr(n)%is_initialized
    enddo
  endif

end subroutine rotate_OBC_segment_tracer_data


!> Allocate an array of data for a field on a segment based on the size of a potentially rotated source array
subroutine allocate_rotated_seg_data(src_array, HI_in, tgt_array, segment)
  real, dimension(:,:,:), intent(in) :: src_array !< The segment data on the unrotated source grid
  type(hor_index_type),   intent(in) :: HI_in !< Horizontal indices on the source grid
  real, dimension(:,:,:), allocatable, intent(inout) :: tgt_array !< The segment data that is being allocated
  type(OBC_segment_type), intent(inout) :: segment !< OBC segment on the target grid

  ! Local variables
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nk
  logical :: corner  ! True if this field is discretized at the OBC segment nodes rather than the faces.

  isd = segment%HI%isd ; ied = segment%HI%ied ; IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
  jsd = segment%HI%jsd ; jed = segment%HI%jed ; JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB
  nk = size(src_array, 3)

  ! Determine whether the source array is allocated at a segment face or at the corners.
  corner = (size(src_array, 1) == abs(HI_in%IedB - HI_in%IsdB) + 1 ) .and. &
           (size(src_array, 2) == abs(HI_in%JedB - HI_in%JsdB) + 1 )

  if (corner) then
    allocate(tgt_array(IsdB:IedB,JsdB:JedB,nk), source=0.0)
  elseif (segment%is_E_or_W) then
    allocate(tgt_array(IsdB:IedB,jsd:jed,nk), source=0.0)
  elseif (segment%is_N_or_S) then
    allocate(tgt_array(isd:ied,JsdB:JedB,nk), source=0.0)
  endif
end subroutine allocate_rotated_seg_data


!> Write out information about the contents of the OBC control structure
subroutine write_OBC_info(OBC, G, GV, US)
  type(ocean_OBC_type),    pointer    :: OBC   !< An open boundary condition control structure
  type(ocean_grid_type),   intent(in) :: G     !< Rotated grid metric
  type(verticalGrid_type), intent(in) :: GV    !< Vertical grid
  type(unit_scale_type),   intent(in) :: US    !< Unit scaling

  ! Local variables
  type(OBC_segment_type), pointer :: segment => NULL() ! pointer to segment type list
  integer :: turns    ! Number of index quarter turns
  integer :: c, n, dir, unrot_dir
  character(len=1024) :: mesg

  turns = modulo(G%HI%turns, 4)

  write(mesg, '("OBC has ", I3, " segments.")') OBC%number_of_segments
  call MOM_mesg(mesg, verb=1)
  !  call MOM_error(WARNING, mesg)

  if (modulo(turns, 2) == 0) then
    if (OBC%open_u_BCs_exist_globally) call MOM_mesg("open_u_BCs_exist_globally", verb=1)
    if (OBC%open_v_BCs_exist_globally) call MOM_mesg("open_v_BCs_exist_globally", verb=1)
    if (OBC%Flather_u_BCs_exist_globally) call MOM_mesg("Flather_u_BCs_exist_globally", verb=1)
    if (OBC%Flather_v_BCs_exist_globally) call MOM_mesg("Flather_v_BCs_exist_globally", verb=1)
    if (OBC%nudged_u_BCs_exist_globally) call MOM_mesg("nudged_u_BCs_exist_globally", verb=1)
    if (OBC%nudged_v_BCs_exist_globally) call MOM_mesg("nudged_v_BCs_exist_globally", verb=1)
    if (OBC%specified_u_BCs_exist_globally) call MOM_mesg("specified_u_BCs_exist_globally", verb=1)
    if (OBC%specified_v_BCs_exist_globally) call MOM_mesg("specified_v_BCs_exist_globally", verb=1)
  else  ! The u- and v-directions are swapped.
    if (OBC%open_v_BCs_exist_globally) call MOM_mesg("open_u_BCs_exist_globally", verb=1)
    if (OBC%open_u_BCs_exist_globally) call MOM_mesg("open_v_BCs_exist_globally", verb=1)
    if (OBC%Flather_v_BCs_exist_globally) call MOM_mesg("Flather_u_BCs_exist_globally", verb=1)
    if (OBC%Flather_u_BCs_exist_globally) call MOM_mesg("Flather_v_BCs_exist_globally", verb=1)
    if (OBC%nudged_v_BCs_exist_globally) call MOM_mesg("nudged_u_BCs_exist_globally", verb=1)
    if (OBC%nudged_u_BCs_exist_globally) call MOM_mesg("nudged_v_BCs_exist_globally", verb=1)
    if (OBC%specified_v_BCs_exist_globally) call MOM_mesg("specified_u_BCs_exist_globally", verb=1)
    if (OBC%specified_u_BCs_exist_globally) call MOM_mesg("specified_v_BCs_exist_globally", verb=1)
  endif

  if (OBC%oblique_BCs_exist_globally) call MOM_mesg("oblique_BCs_exist_globally", verb=1)
  if (OBC%radiation_BCs_exist_globally) call MOM_mesg("radiation_BCs_exist_globally", verb=1)
  if (OBC%user_BCs_set_globally) call MOM_mesg("user_BCs_set_globally", verb=1)
  if (OBC%update_OBC) call MOM_mesg("update_OBC", verb=1)
  if (OBC%update_OBC_seg_data) call MOM_mesg("update_OBC_seg_data", verb=1)
  if (OBC%needs_IO_for_data) call MOM_mesg("needs_IO_for_data", verb=1)
  if (OBC%any_needs_IO_for_data) call MOM_mesg("any_needs_IO_for_data", verb=1)
  if (OBC%zero_vorticity) call MOM_mesg("zero_vorticity", verb=1)
  if (OBC%freeslip_vorticity) call MOM_mesg("freeslip_vorticity", verb=1)
  if (OBC%computed_vorticity) call MOM_mesg("computed_vorticity", verb=1)
  if (OBC%specified_vorticity) call MOM_mesg("specified_vorticity", verb=1)
  if (OBC%zero_strain) call MOM_mesg("zero_strain", verb=1)
  if (OBC%freeslip_strain) call MOM_mesg("freeslip_strain", verb=1)
  if (OBC%computed_strain) call MOM_mesg("computed_strain", verb=1)
  if (OBC%specified_strain) call MOM_mesg("specified_strain", verb=1)
  if (OBC%zero_biharmonic) call MOM_mesg("zero_biharmonic", verb=1)
  if (OBC%brushcutter_mode) call MOM_mesg("brushcutter_mode", verb=1)
  if (OBC%check_reconstruction) call MOM_mesg("check_reconstruction", verb=1)
  if (OBC%check_remapping) call MOM_mesg("check_remapping", verb=1)
  if (OBC%force_bounds_in_subcell) call MOM_mesg("force_bounds_in_subcell", verb=1)
  if (OBC%om4_remap_via_sub_cells) call MOM_mesg("om4_remap_via_sub_cells", verb=1)
  if (OBC%exterior_OBC_bug) call MOM_mesg("exterior_OBC_bug", verb=1)
  if (OBC%hor_index_bug) call MOM_mesg("hor_index_bug", verb=1)
  if (OBC%debug) call MOM_mesg("debug", verb=1)
  if (OBC%ramp) call MOM_mesg("ramp", verb=1)
  if (OBC%ramping_is_activated) call MOM_mesg("ramping_is_activated", verb=1)
  write(mesg, '("n_tide_constituents ", I3)') OBC%n_tide_constituents
  call MOM_mesg(mesg, verb=1)
  if (OBC%n_tide_constituents > 0) then
    do c=1,OBC%n_tide_constituents
      write(mesg, '(" properties ", 4ES16.6)') &
            US%s_to_T*OBC%tide_frequencies(c), OBC%tide_eq_phases(c), OBC%tide_fn(c), OBC%tide_un(c)
      call MOM_mesg(trim(OBC%tide_names(c))//mesg, verb=1)
    enddo
  endif
  if (OBC%ramp) then
    write(mesg, '("ramp_values ", 3ES16.6)') OBC%ramp_timescale, OBC%trunc_ramp_time, OBC%ramp_value
    call MOM_mesg(mesg, verb=1)
  endif
  write(mesg, '("gamma_uv ", ES16.6)') OBC%gamma_uv
  call MOM_mesg(mesg, verb=1)
  write(mesg, '("rx_max ", ES16.6)') OBC%rx_max
  call MOM_mesg(mesg, verb=1)

  call MOM_mesg("remappingScheme = "//trim(OBC%remappingScheme), verb=1)

  do n=1,OBC%number_of_segments
    segment => OBC%segment(n)
    dir = segment%direction

    unrot_dir = rotate_OBC_segment_direction(dir, -turns)
    write(mesg, '(" Segment ", I3, " has direction ", I3)') n, unrot_dir
    if (unrot_dir == OBC_DIRECTION_N)  write(mesg, '(" Segment ", I3, " is Northern")') n
    if (unrot_dir == OBC_DIRECTION_S)  write(mesg, '(" Segment ", I3, " is Southern")') n
    if (unrot_dir == OBC_DIRECTION_E)  write(mesg, '(" Segment ", I3, " is Eastern")') n
    if (unrot_dir == OBC_DIRECTION_W)  write(mesg, '(" Segment ", I3, " is Western")') n
    call MOM_mesg(mesg, verb=1)

    ! write(mesg, '("  range: ", 4I3)') segment%Is_obc, segment%Ie_obc, segment%Js_obc, segment%Je_obc
    if (modulo(turns, 2) == 0) then
      write(mesg, '("  size: ", 4I3)') 1+abs(segment%Ie_obc-segment%Is_obc), 1+abs(segment%Je_obc-segment%Js_obc)
    else
      write(mesg, '("  size: ", 4I3)') 1+abs(segment%Je_obc-segment%Js_obc), 1+abs(segment%Ie_obc-segment%Is_obc)
    endif
    call MOM_mesg(mesg, verb=1)

    if (segment%on_pe)          call MOM_mesg("  Segment is on PE.", verb=1)

    if (segment%Flather)        call MOM_mesg("  Flather", verb=1)
    if (segment%radiation)      call MOM_mesg("  radiation", verb=1)
    if (segment%radiation_tan)  call MOM_mesg("  radiation_tan", verb=1)
    if (segment%radiation_grad) call MOM_mesg("  radiation_grad", verb=1)
    if (segment%oblique)        call MOM_mesg("  oblique", verb=1)
    if (segment%oblique_tan)    call MOM_mesg("  oblique_tan", verb=1)
    if (segment%oblique_grad)   call MOM_mesg("  oblique_grad", verb=1)
    if (segment%nudged)         call MOM_mesg("  nudged", verb=1)
    if (segment%nudged_tan)     call MOM_mesg("  nudged_tan", verb=1)
    if (segment%nudged_grad)    call MOM_mesg("  nudged_grad", verb=1)
    if (segment%specified)      call MOM_mesg("  specified", verb=1)
    if (segment%specified_tan)  call MOM_mesg("  specified_tan", verb=1)
    if (segment%specified_grad) call MOM_mesg("  specified_grad", verb=1)
    if (segment%open)           call MOM_mesg("  open", verb=1)
    if (segment%gradient)       call MOM_mesg("  gradient", verb=1)
    if (segment%values_needed)  call MOM_mesg("  values_needed", verb=1)
    if (modulo(turns, 2) == 0) then
      if (segment%is_N_or_S)      call MOM_mesg("  is_N_or_S", verb=1)
      if (segment%is_E_or_W)      call MOM_mesg("  is_E_or_W", verb=1)
      if (segment%u_values_needed) call MOM_mesg("  u_values_needed", verb=1)
      if (segment%uamp_values_needed) call MOM_mesg("  uamp_values_needed", verb=1)
      if (segment%uphase_values_needed) call MOM_mesg("  uphase_values_needed", verb=1)
      if (segment%v_values_needed) call MOM_mesg("  v_values_needed", verb=1)
      if (segment%vamp_values_needed) call MOM_mesg("  vamp_values_needed", verb=1)
      if (segment%vphase_values_needed) call MOM_mesg("  vphase_values_needed", verb=1)
    else  ! The x- and y-directions are swapped.
      if (segment%is_E_or_W)      call MOM_mesg("  is_N_or_S", verb=1)
      if (segment%is_N_or_S)      call MOM_mesg("  is_E_or_W", verb=1)
      if (segment%v_values_needed) call MOM_mesg("  u_values_needed", verb=1)
      if (segment%vamp_values_needed) call MOM_mesg("  uamp_values_needed", verb=1)
      if (segment%vphase_values_needed) call MOM_mesg("  uphase_values_needed", verb=1)
      if (segment%u_values_needed) call MOM_mesg("  v_values_needed", verb=1)
      if (segment%uamp_values_needed) call MOM_mesg("  vamp_values_needed", verb=1)
      if (segment%uphase_values_needed) call MOM_mesg("  vphase_values_needed", verb=1)
    endif
    if (segment%t_values_needed) call MOM_mesg("  t_values_needed", verb=1)
    if (segment%s_values_needed) call MOM_mesg("  s_values_needed", verb=1)
    if (segment%z_values_needed) call MOM_mesg("  z_values_needed", verb=1)
    if (segment%zamp_values_needed) call MOM_mesg("  zamp_values_needed", verb=1)
    if (segment%zphase_values_needed) call MOM_mesg("  zphase_values_needed", verb=1)
    if (segment%g_values_needed) call MOM_mesg("  g_values_needed", verb=1)
!    if (segment%is_E_or_W_2)    call MOM_mesg("  is_E_or_W_2", verb=1)
    if (segment%temp_segment_data_exists) call MOM_mesg("  temp_segment_data_exists", verb=1)
    if (segment%salt_segment_data_exists) call MOM_mesg("  salt_segment_data_exists", verb=1)

    write(mesg, '("  Tr_InvLscale_out ", ES16.6)') segment%Tr_InvLscale_out*US%m_to_L
    call MOM_mesg(mesg, verb=1)
    write(mesg, '("  Tr_InvLscale_in ", ES16.6)') segment%Tr_InvLscale_in*US%m_to_L
    call MOM_mesg(mesg, verb=1)

  enddo

  call chksum_OBC_segments(OBC, G, GV, US, 0)

end subroutine write_OBC_info

!> Write checksums and perhaps the values of all the allocated arrays on an OBC segments.
subroutine chksum_OBC_segments(OBC, G, GV, US, nk)
  type(ocean_OBC_type),    pointer    :: OBC   !< OBC on input map
  type(ocean_grid_type),   intent(in) :: G     !< Rotated grid metric
  type(verticalGrid_type), intent(in) :: GV    !< Vertical grid
  type(unit_scale_type),   intent(in) :: US    !< Unit scaling
  integer,                 intent(in) :: nk    !< The number of layers to print

  ! Local variables
  type(OBC_segment_type), pointer :: segment => NULL() ! pointer to segment type list
  real :: norm ! A sign change used when rotating a normal component [nondim]
  real :: tang ! A sign change used when rotating a tangential component [nondim]
  character(len=8) :: sn, segno
  character(len=1024) :: mesg
  integer :: c, n, dir

  do n=1,OBC%number_of_segments
    segment => OBC%segment(n)
    dir = segment%direction

    write(segno, '(I3)') n
    sn = '('//trim(adjustl(segno))//')'

    ! Turn each segment and write it as though it is an eastern face.
    norm = 0.0 ; tang = 0.0
    if (dir == OBC_DIRECTION_E) then
      norm = 1.0 ; tang = 1.0
    elseif (dir == OBC_DIRECTION_N) then
      norm = 1.0 ; tang = -1.0
    elseif (dir == OBC_DIRECTION_W) then
      norm = -1.0 ; tang = -1.0
    elseif (dir == OBC_DIRECTION_S) then
      norm = -1.0 ; tang = 1.0
    endif

    if (allocated(segment%Cg)) call write_2d_array_vals("Cg"//trim(sn), segment%Cg, dir, nk, unscale=US%L_T_to_m_s)
    if (allocated(segment%Htot)) call write_2d_array_vals("Htot"//trim(sn), segment%Htot, dir, nk, unscale=GV%H_to_mks)
    if (allocated(segment%dZtot)) call write_2d_array_vals("dZtot"//trim(sn), segment%dZtot, dir, nk, unscale=US%Z_to_m)
    if (allocated(segment%SSH)) call write_2d_array_vals("SSH"//trim(sn), segment%SSH, dir, nk, unscale=US%Z_to_m)
    if (allocated(segment%h)) call write_3d_array_vals("h"//trim(sn), segment%h, dir, nk, unscale=GV%H_to_mks)
    if (allocated(segment%normal_vel)) &
      call write_3d_array_vals("normal_vel"//trim(sn), segment%normal_vel, dir, nk, unscale=norm*US%L_T_to_m_s)
    if (allocated(segment%normal_vel_bt)) &
      call write_2d_array_vals("normal_vel_bt"//trim(sn), segment%normal_vel_bt, dir, nk, unscale=norm*US%L_T_to_m_s)
    if (allocated(segment%tangential_vel)) &
      call write_3d_array_vals("tangential_vel"//trim(sn), segment%tangential_vel, dir, nk, unscale=tang*US%L_T_to_m_s)
    if (allocated(segment%tangential_grad)) &
      call write_3d_array_vals("tangential_grad"//trim(sn), segment%tangential_grad, dir, nk, &
                    unscale=tang*norm*US%s_to_T)
    if (allocated(segment%normal_trans)) &
      call write_3d_array_vals("normal_trans"//trim(sn), segment%normal_trans, dir, nk, &
                    unscale=norm*GV%H_to_mks*US%L_T_to_m_s*US%L_to_m)
    if (allocated(segment%grad_normal)) &
      call write_3d_array_vals("grad_normal"//trim(sn), segment%grad_normal, dir, nk, unscale=norm*tang*US%L_T_to_m_s)
    if (allocated(segment%grad_tan)) &
      call write_3d_array_vals("grad_tan"//trim(sn), segment%grad_tan, dir, nk, unscale=1.0*US%L_T_to_m_s)
    if (allocated(segment%grad_gradient)) &
      call write_3d_array_vals("grad_gradient"//trim(sn), segment%grad_gradient, dir, nk, unscale=norm*US%s_to_T)

    if (allocated(segment%rx_norm_rad)) &
      call write_3d_array_vals("rxy_norm_rad"//trim(sn), segment%rx_norm_rad, dir, nk, unscale=1.0)
    if (allocated(segment%ry_norm_rad)) &
      call write_3d_array_vals("rxy_norm_rad"//trim(sn), segment%ry_norm_rad, dir, nk, unscale=1.0)
    if (segment%is_E_or_W) then
      if (allocated(segment%rx_norm_obl)) &
        call write_3d_array_vals("rx_norm_obl"//trim(sn), segment%rx_norm_obl, dir, nk, unscale=US%L_T_to_m_s**2)
      if (allocated(segment%ry_norm_obl)) &
        call write_3d_array_vals("ry_norm_obl"//trim(sn), segment%ry_norm_obl, dir, nk, unscale=US%L_T_to_m_s**2)
    else ! The x- and y- directions are swapped.
      if (allocated(segment%ry_norm_obl)) &
        call write_3d_array_vals("rx_norm_obl"//trim(sn), segment%ry_norm_obl, dir, nk, unscale=US%L_T_to_m_s**2)
      if (allocated(segment%rx_norm_obl)) &
        call write_3d_array_vals("ry_norm_obl"//trim(sn), segment%rx_norm_obl, dir, nk, unscale=US%L_T_to_m_s**2)
    endif

    if (allocated(segment%cff_normal)) &
      call write_3d_array_vals("cff_normal"//trim(sn), segment%cff_normal, dir, nk, unscale=US%L_T_to_m_s**2)
    if (allocated(segment%nudged_normal_vel)) &
      call write_3d_array_vals("nudged_normal_vel"//trim(sn), segment%nudged_normal_vel, dir, nk, &
                    unscale=norm*US%L_T_to_m_s)
    if (allocated(segment%nudged_tangential_vel)) &
      call write_3d_array_vals("nudged_tangential_vel"//trim(sn), segment%nudged_tangential_vel, dir, nk, &
                    unscale=tang*US%L_T_to_m_s)
    if (allocated(segment%nudged_tangential_grad)) &
      call write_3d_array_vals("nudged_tangential_grad"//trim(sn), segment%nudged_tangential_grad, dir, nk, &
                    unscale=tang*norm*US%s_to_T)
  enddo

  contains

  !> Write out the values in a named 2-d segment data array
  subroutine write_2d_array_vals(name, Array, seg_dir, nkp, unscale)
    character(len=*),     intent(in) :: name    !< The name of the variable
    real, dimension(:,:), intent(in) :: Array   !< The 2-d array to write [A ~> a]
    integer,              intent(in) :: seg_dir !< The direction of the segment
    integer,              intent(in) :: nkp     !< Print all the values if this is greater than 0
    real,       optional, intent(in) :: unscale !< A factor that undoes the scaling of the array [a A-1 ~> 1]
    ! Local variables
    real :: scale  !  A factor that undoes the scaling of the array [a A-1 ~> 1]
    character(len=1024) :: mesg
    character(len=24) :: val
    integer :: i, j, n, iounit

    scale = 1.0 ; if (present(unscale)) scale = unscale
    iounit = stderr

    if (nkp > 0) then
      write(iounit, '(2X,A,":")') trim(name)
      mesg = "" ; n = 0
      if ((seg_dir == OBC_DIRECTION_N) .or. (seg_dir == OBC_DIRECTION_W)) then
        do j=size(Array,2),1,-1 ; do i=size(Array,1),1,-1
          write(val, '(ES16.6)') scale*Array(i,j)
          mesg = trim(mesg)//" "//trim(val) ;  n = n + 1
          if (n >= 12) then
            write(iounit, '(2X,A)') trim(mesg)
            mesg = "" ; n = 0
          endif
        enddo ; enddo
      else
        do j=1,size(Array,2) ; do i=1,size(Array,1)
          write(val, '(ES16.6)') scale*Array(i,j)
          mesg = trim(mesg)//" "//trim(val) ;  n = n + 1
          if (n >= 12) then
            write(iounit, '(2X,A)') trim(mesg)
            mesg = "" ; n = 0
          endif
        enddo ; enddo
      endif
      if (n > 0) write(iounit, '(2X,A)') trim(mesg)
    endif

    if (scale == 1.0) then
      call chksum(Array, name)
    else
      call chksum(scale*Array(:,:), name)
    endif
  end subroutine write_2d_array_vals

  !> Write out the values in a 3-d segment data array
  subroutine write_3d_array_vals(name, Array, seg_dir, nkp, unscale)
    character(len=*),       intent(in) :: name    !< The name of the variable
    real, dimension(:,:,:), intent(in) :: Array   !< The 3-d array to write
    integer,                intent(in) :: seg_dir !< The direction of the segment
    integer,                intent(in) :: nkp     !< The number of layers to print
    real,         optional, intent(in) :: unscale !< A factor that undoes the scaling of the array [a A-1 ~> 1]
    ! Local variables
    real :: scale  !  A factor that undoes the scaling of the array [a A-1 ~> 1]
    logical :: reverse
    character(len=1024) :: mesg
    character(len=24) :: val
    integer :: i, j, k, n, nk, iounit

    scale = 1.0 ; if (present(unscale)) scale = unscale
    iounit = stderr

    if (nkp > 0) then
      nk = min(nkp, size(Array,3))
      write(iounit, '(2X,A,":")') trim(name)
      do k=1,nk
        mesg = "" ; n = 0
        if ((seg_dir == OBC_DIRECTION_N) .or. (seg_dir == OBC_DIRECTION_W)) then
          do j=size(Array,2),1,-1 ; do i=size(Array,1),1,-1
            write(val, '(ES16.6)') scale*Array(i,j,k)
            mesg = trim(mesg)//" "//trim(val) ; n = n + 1
            if (n >= 12) then
              write(iounit, '(2X,A)') trim(mesg)
              mesg = "" ; n = 0
            endif
          enddo ; enddo
        else
          do j=1,size(Array,2) ; do i=1,size(Array,1)
            write(val, '(ES16.6)') scale*Array(i,j,k)
            mesg = trim(mesg)//" "//trim(val) ; n = n + 1
            if (n >= 12) then
              write(iounit, '(2X,A)') trim(mesg)
              mesg = "" ; n = 0
            endif
          enddo ; enddo
        endif
        if (n > 0) write(iounit, '(2X,A)') trim(mesg)
      enddo
    endif

    if (scale == 1.0) then
      call chksum(Array, name)
    else
      call chksum(scale*Array(:,:,:), name)
    endif

  end subroutine write_3d_array_vals

end subroutine chksum_OBC_segments

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
