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
use MOM_restart,              only : register_restart_field, MOM_restart_CS
use MOM_obsolete_params,      only : obsolete_logical, obsolete_int, obsolete_real, obsolete_char
use MOM_string_functions,     only : extract_word, remove_spaces
use MOM_tracer_registry,      only : tracer_registry_type
use MOM_variables,            only : thermo_var_ptrs
use time_interp_external_mod, only : init_external_field, time_interp_external
use MOM_remapping,            only : remappingSchemesDoc, remappingDefaultScheme, remapping_CS
use MOM_remapping,            only : initialize_remapping, remapping_core_h, end_remapping
use MOM_regridding,           only : regridding_CS
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
public update_obc_segment_data
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

integer, parameter, public :: OBC_NONE = 0, OBC_SIMPLE = 1, OBC_WALL = 2
integer, parameter, public :: OBC_FLATHER = 3
integer, parameter, public :: OBC_RADIATION = 4
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
  real, dimension(:,:,:), pointer :: dz_src=>NULL()     !< vertical grid cell spacing of the incoming segment data (m)
  real, dimension(:,:,:), pointer :: buffer_dst=>NULL() !< buffer src data remapped to the target vertical grid
  real, dimension(:,:), pointer   :: bt_vel=>NULL()     !< barotropic velocity (m s-1)
  real                            :: value              !< constant value if fid is equal to -1
end type OBC_segment_data_type

!> Tracer segment data structure, for putting into an array of objects, not all the same shape.
!type, public :: segment_tracer_type
!  real, dimension(:,:,:), pointer :: tr         => NULL()  !< tracer concentration array
!end type segment_tracer_type

!> Tracer on OBC segment data structure, for putting into a segment tracer registry.
type, public :: OBC_segment_tracer_type
  real, dimension(:,:,:), pointer :: t          => NULL()  !< tracer concentration array
  real                            :: OBC_inflow_conc = 0.0 !< tracer concentration for generic inflows
  character(len=32)               :: name                  !< tracer name used for error messages
  type(vardesc), pointer          :: vd         => NULL()  !< metadata describing the tracer
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
  logical :: oblique        !< Oblique waves supported at radiation boundary.
  logical :: nudged         !< Optional supplement to radiation boundary.
  logical :: specified      !< Boundary fixed to external value.
  logical :: open           !< Boundary is open for continuity solver.
  logical :: gradient       !< Zero gradient at boundary.
  logical :: values_needed  !< Whether or not external OBC fields are needed.
  integer :: direction      !< Boundary faces one of the four directions.
  logical :: is_N_or_S      !< True is the OB is facing North or South and exists on this PE.
  logical :: is_E_or_W      !< True is the OB is facing East or West and exists on this PE.
  type(OBC_segment_data_type), pointer, dimension(:) :: field=>NULL()   !<  OBC data
  integer :: num_fields     !< number of OBC data fields (e.g. u_normal,u_parallel and eta for Flather)
  character(len=32), pointer, dimension(:) :: field_names=>NULL() !< field names for this segment
  integer :: Is_obc         !< i-indices of boundary segment.
  integer :: Ie_obc         !< i-indices of boundary segment.
  integer :: Js_obc         !< j-indices of boundary segment.
  integer :: Je_obc         !< j-indices of boundary segment.
  real :: Tnudge_in         !< Inverse nudging timescale on inflow (1/s).
  real :: Tnudge_out        !< Inverse nudging timescale on outflow (1/s).
  logical :: on_pe          !< true if segment is located in the computational domain
  logical :: temp_segment_data_exists !< true if temperature data arrays are present
  logical :: salt_segment_data_exists !< true if salinity data arrays are present
  real, pointer, dimension(:,:)   :: Cg=>NULL()     !< The external gravity
                                                    !! wave speed (m -s) at OBC-points.
  real, pointer, dimension(:,:)   :: Htot=>NULL()   !< The total column thickness (m) at OBC-points.
  real, pointer, dimension(:,:,:) :: h=>NULL()      !< The cell thickness (m) at OBC-points.
  real, pointer, dimension(:,:,:) :: normal_vel=>NULL()     !< The layer velocity normal to the OB
                                                            !! segment (m s-1).
  real, pointer, dimension(:,:,:) :: normal_trans=>NULL()   !< The layer transport normal to the OB
                                                            !! segment (m3 s-1).
  real, pointer, dimension(:,:)   :: normal_vel_bt=>NULL()  !< The barotropic velocity normal to
                                                            !! the OB segment (m s-1).
  real, pointer, dimension(:,:)   :: normal_trans_bt=>NULL()!< The barotropic transport normal to
                                                            !! the OB segment (m3 s-1).
  real, pointer, dimension(:,:)   :: eta=>NULL()            !< The sea-surface elevation along the segment (m).
  real, pointer, dimension(:,:,:) :: grad_normal=>NULL()    !< The gradient of the normal flow along the
                                                            !! segment (m s-1)
  real, pointer, dimension(:,:,:) :: rx_normal=>NULL()      !< The rx_old_u value for radiation coeff
                                                            !! for normal velocity
  real, pointer, dimension(:,:,:) :: nudged_normal_vel=>NULL() !< The layer velocity normal to the OB segment
                                                            !! that values should be nudged towards (m s-1).
  type(segment_tracer_registry_type), pointer  :: tr_Reg=> NULL()!< A pointer to the tracer registry for the segment.
  type(hor_index_type) :: HI !< Horizontal index ranges
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
  logical :: user_BCs_set_globally = .false.          !< True if any OBC_USER_CONFIG is set
                                                      !! for input from user directory.
  logical :: update_OBC = .false.                     !< Is OBC data time-dependent
  logical :: needs_IO_for_data = .false.              !< Is any i/o needed for OBCs
  logical :: zero_vorticity = .false.                 !< If True, sets relative vorticity to zero on open boundaries.
  logical :: freeslip_vorticity = .false.             !< If True, sets normal gradient of tangential velocity to zero
                                                      !! in the relative vorticity on open boundaries.
  logical :: zero_strain = .false.                    !< If True, sets strain to zero on open boundaries.
  logical :: freeslip_strain = .false.                !< If True, sets normal gradient of tangential velocity to zero
                                                      !! in the strain on open boundaries.
  logical :: zero_biharmonic = .false.                !< If True, zeros the Laplacian of flow on open boundaries for
                                                      !! use in the biharmonic viscosity term.
  logical :: extend_segments = .false.                !< If True, extend OBC segments (for testing)
  real :: g_Earth
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
  real :: gamma_h  !< The relative weighting for the baroclinic radiation
                   !! velocities (or speed of characteristics) at the
                   !! new time level (1) or the running mean (0) for thicknesses.
                   !! Valid values range from 0 to 1, with a default of 0.2.
  real :: rx_max   !< The maximum magnitude of the baroclinic radiation
                   !! velocity (or speed of characteristics), in m s-1.  The
                   !! default value is 10 m s-1.
  logical :: OBC_pe !< Is there an open boundary on this tile?
  type(remapping_CS),      pointer :: remap_CS   !< ALE remapping control structure for segments only
  type(OBC_registry_type), pointer :: OBC_Reg => NULL()  !< Registry type for boundaries
  real, pointer, dimension(:,:,:) :: rx_normal => NULL()  !< Array storage for restarts
  real, pointer, dimension(:,:,:) :: ry_normal => NULL()  !< Array storage for restarts
  real :: silly_h  !< A silly value of thickness outside of the domain that
                   !! can be used to test the independence of the OBCs to
                   !! this external data, in m.
  real :: silly_u  !< A silly value of velocity outside of the domain that
                   !! can be used to test the independence of the OBCs to
                   !! this external data, in m/s.
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

integer :: id_clock_pass

character(len=40)  :: mdl = "MOM_open_boundary" ! This module's name.
! This include declares and sets the variable "version".
#include "version_variable.h"

contains

!> Enables OBC module and reads configuration parameters
!> This routine is called from MOM_initialize_fixed which
!> occurs before the initialization of the vertical coordinate
!> and ALE_init.  Therefore segment data are not fully initialized
!> here. The remainder of the segment data are initialized in a
!> later call to update_open_boundary_data

subroutine open_boundary_config(G, param_file, OBC)
  type(dyn_horgrid_type),  intent(inout) :: G !< Ocean grid structure
  type(param_file_type),   intent(in)    :: param_file !< Parameter file handle
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundary control structure
  ! Local variables
  integer :: l ! For looping over segments
  logical :: debug_OBC, debug, mask_outside
  character(len=15) :: segment_param_str ! The run-time parameter name for each segment
  character(len=100) :: segment_str      ! The contents (rhs) for parameter "segment_param_str"
  character(len=200) :: config1          ! String for OBC_USER_CONFIG

  allocate(OBC)

  call log_version(param_file, mdl, version, "Controls where open boundaries are located, what "//&
                 "kind of boundary condition to impose, and what data to apply, if any.")
  call get_param(param_file, mdl, "OBC_NUMBER_OF_SEGMENTS", OBC%number_of_segments, &
                 "The number of open boundary segments.", &
                 default=0)
  call get_param(param_file, mdl, "G_EARTH", OBC%g_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80)
  call get_param(param_file, mdl, "OBC_USER_CONFIG", config1, &
                 "A string that sets how the open boundary conditions are \n"//&
                 " configured: \n", default="none", do_not_log=.true.)
  call get_param(param_file, mdl, "NK", OBC%ke, &
                 "The number of model layers", default=0, do_not_log=.true.)

  if (config1 .ne. "none") OBC%user_BCs_set_globally = .true.

  call get_param(param_file, mdl, "EXTEND_OBC_SEGMENTS", OBC%extend_segments, &
                   "If true, extend OBC segments. This option is used to recover\n"//&
                   "legacy solutions dependent on an incomplete implementaion of OBCs.\n"//&
                   "This option will be obsoleted in the future.", default=.false.)

  if (OBC%number_of_segments > 0) then
    call get_param(param_file, mdl, "OBC_ZERO_VORTICITY", OBC%zero_vorticity, &
         "If true, sets relative vorticity to zero on open boundaries.", &
         default=.false.)
    call get_param(param_file, mdl, "OBC_FREESLIP_VORTICITY", OBC%freeslip_vorticity, &
         "If true, sets the normal gradient of tangential velocity to\n"// &
         "zero in the relative vorticity on open boundaries. This cannot\n"// &
         "be true if OBC_ZERO_VORTICITY is True.", default=.false.)
    if (OBC%zero_vorticity .and. OBC%freeslip_vorticity) call MOM_error(FATAL, &
         "MOM_open_boundary.F90, open_boundary_config: "//&
         "Only one of OBC_ZERO_VORTICITY and OBC_FREESLIP_VORTICITY can be True at once.")
    call get_param(param_file, mdl, "OBC_ZERO_STRAIN", OBC%zero_strain, &
         "If true, sets the strain used in the stress tensor to zero on open boundaries.", &
         default=.false.)
    call get_param(param_file, mdl, "OBC_FREESLIP_STRAIN", OBC%freeslip_strain, &
         "If true, sets the normal gradient of tangential velocity to\n"// &
         "zero in the strain use in the stress tensor on open boundaries. This cannot\n"// &
         "be true if OBC_ZERO_STRAIN is True.", default=.false.)
    if (OBC%zero_strain .and. OBC%freeslip_strain) call MOM_error(FATAL, &
         "MOM_open_boundary.F90, open_boundary_config: "//&
         "Only one of OBC_ZERO_STRAIN and OBC_FREESLIP_STRAIN can be True at once.")
    call get_param(param_file, mdl, "OBC_ZERO_BIHARMONIC", OBC%zero_biharmonic, &
         "If true, zeros the Laplacian of flow on open boundaries in the biharmonic\n"//&
         "viscosity term.", default=.false.)
    call get_param(param_file, mdl, "MASK_OUTSIDE_OBCS", mask_outside, &
         "If true, set the areas outside open boundaries to be land.", &
         default=.false.)

    call get_param(param_file, mdl, "DEBUG", debug, default=.false.)
    call get_param(param_file, mdl, "DEBUG_OBC", debug_OBC, default=.false.)
    if (debug_OBC .or. debug) &
      call log_param(param_file, mdl, "DEBUG_OBC", debug_OBC, &
                 "If true, do additional calls to help debug the performance \n"//&
                 "of the open boundary condition code.", default=.false.)

    call get_param(param_file, mdl, "OBC_SILLY_THICK", OBC%silly_h, &
                 "A silly value of thicknesses used outside of open boundary \n"//&
                 "conditions for debugging.", units="m", default=0.0, &
                 do_not_log=.not.debug_OBC)
    call get_param(param_file, mdl, "OBC_SILLY_VEL", OBC%silly_u, &
                 "A silly value of velocities used outside of open boundary \n"//&
                 "conditions for debugging.", units="m/s", default=0.0, &
                 do_not_log=.not.debug_OBC)

    ! Allocate everything
    ! Note the 0-segment is needed when %segnum_u/v(:,:) = 0
    allocate(OBC%segment(0:OBC%number_of_segments))
    do l=0,OBC%number_of_segments
      OBC%segment(l)%Flather = .false.
      OBC%segment(l)%radiation = .false.
      OBC%segment(l)%oblique = .false.
      OBC%segment(l)%nudged = .false.
      OBC%segment(l)%specified = .false.
      OBC%segment(l)%open = .false.
      OBC%segment(l)%gradient = .false.
      OBC%segment(l)%values_needed = .false.
      OBC%segment(l)%direction = OBC_NONE
      OBC%segment(l)%is_N_or_S = .false.
      OBC%segment(l)%is_E_or_W = .false.
      OBC%segment(l)%Tnudge_in = 0.0
      OBC%segment(l)%Tnudge_out = 0.0
      OBC%segment(l)%num_fields = 0.0
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
        call setup_u_point_obc(OBC, G, segment_str, l)
      elseif (segment_str(1:2) == 'J=') then
        call setup_v_point_obc(OBC, G, segment_str, l)
      else
        call MOM_error(FATAL, "MOM_open_boundary.F90, open_boundary_config: "//&
             "Unable to interpret "//segment_param_str//" = "//trim(segment_str))
      endif
    enddo

    !    if (open_boundary_query(OBC, needs_ext_seg_data=.true.)) &
    call initialize_segment_data(G, OBC, param_file)

    if ( OBC%Flather_u_BCs_exist_globally .or. OBC%Flather_v_BCs_exist_globally ) then
      call get_param(param_file, mdl, "OBC_RADIATION_MAX", OBC%rx_max, &
                   "The maximum magnitude of the baroclinic radiation \n"//&
                   "velocity (or speed of characteristics).  This is only \n"//&
                   "used if one of the open boundary segments is using Orlanski.", &
                   units="m s-1", default=10.0)
      call get_param(param_file, mdl, "OBC_RAD_VEL_WT", OBC%gamma_uv, &
                   "The relative weighting for the baroclinic radiation \n"//&
                   "velocities (or speed of characteristics) at the new \n"//&
                   "time level (1) or the running mean (0) for velocities. \n"//&
                   "Valid values range from 0 to 1. This is only used if \n"//&
                   "one of the open boundary segments is using Orlanski.", &
                   units="nondim",  default=0.3)
      call get_param(param_file, mdl, "OBC_RAD_THICK_WT", OBC%gamma_h, &
                   "The relative weighting for the baroclinic radiation \n"//&
                   "velocities (or speed of characteristics) at the new \n"//&
                   "time level (1) or the running mean (0) for thicknesses. \n"//&
                   "Valid values range from 0 to 1. This is only used if \n"//&
                   "one of the open boundary segments is using Orlanski.", &
                   units="nondim",  default=0.2)
    endif
    if (mask_outside) call mask_outside_OBCs(G, param_file, OBC)
  endif

    ! Safety check
  if ((OBC%open_u_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally) .and. &
       .not.G%symmetric ) call MOM_error(FATAL, &
       "MOM_open_boundary, open_boundary_config: "//&
       "Symmetric memory must be used when using Flather OBCs.")

  if (.not.(OBC%specified_u_BCs_exist_globally .or. OBC%specified_v_BCs_exist_globally .or. &
              OBC%open_u_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally)) then
    ! No open boundaries have been requested
    call open_boundary_dealloc(OBC)
  endif

end subroutine open_boundary_config

subroutine initialize_segment_data(G, OBC, PF)
  use mpp_mod, only : mpp_pe, mpp_set_current_pelist, mpp_get_current_pelist,mpp_npes

  type(dyn_horgrid_type),  intent(in) :: G   !< Ocean grid structure
  type(ocean_OBC_type), intent(inout) :: OBC !< Open boundary control structure
  type(param_file_type), intent(in)   :: PF  !< Parameter file handle

  integer :: n,m,num_fields
  character(len=256) :: segstr, filename
  character(len=20)  :: segnam, suffix
  character(len=32)  :: varnam, fieldname
  real               :: value
  integer            :: orient
  character(len=32), dimension(MAX_OBC_FIELDS) :: fields  ! segment field names
  character(len=128) :: inputdir
  type(OBC_segment_type), pointer :: segment ! pointer to segment type list
  character(len=32)  :: remappingScheme
  logical :: check_reconstruction, check_remapping, force_bounds_in_subcell
  integer, dimension(4) :: siz,siz2
  integer :: is, ie, js, je
  integer :: isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB
  integer, dimension(:), allocatable :: saved_pelist
  integer :: current_pe
  integer, dimension(1) :: single_pelist

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  ! There is a problem with the order of the OBC initialization
  ! with respect to ALE_init. Currently handling this by copying the
  ! param file so that I can use it later in step_MOM in order to finish
  ! initializing segments on the first step.

  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)

  call get_param(PF, mdl, "REMAPPING_SCHEME", remappingScheme, &
          "This sets the reconstruction scheme used\n"//&
          "for vertical remapping for all variables.\n"//&
          "It can be one of the following schemes:\n"//&
          trim(remappingSchemesDoc), default=remappingDefaultScheme,do_not_log=.true.)
  call get_param(PF, mdl, "FATAL_CHECK_RECONSTRUCTIONS", check_reconstruction, &
          "If true, cell-by-cell reconstructions are checked for\n"//&
          "consistency and if non-monotonicity or an inconsistency is\n"//&
          "detected then a FATAL error is issued.", default=.false.,do_not_log=.true.)
  call get_param(PF, mdl, "FATAL_CHECK_REMAPPING", check_remapping, &
          "If true, the results of remapping are checked for\n"//&
          "conservation and new extrema and if an inconsistency is\n"//&
          "detected then a FATAL error is issued.", default=.false.,do_not_log=.true.)
  call get_param(PF, mdl, "REMAP_BOUND_INTERMEDIATE_VALUES", force_bounds_in_subcell, &
          "If true, the values on the intermediate grid used for remapping\n"//&
          "are forced to be bounded, which might not be the case due to\n"//&
          "round off.", default=.false.,do_not_log=.true.)

  allocate(OBC%remap_CS)
  call initialize_remapping(OBC%remap_CS, remappingScheme, boundary_extrapolation = .false., &
       check_reconstruction=check_reconstruction, &
       check_remapping=check_remapping, force_bounds_in_subcell=force_bounds_in_subcell)

  if (OBC%user_BCs_set_globally) return


  !< temporarily disable communication in order to read segment data independently

  allocate(saved_pelist(0:mpp_npes()-1))
  call mpp_get_current_pelist(saved_pelist)
  current_pe = mpp_pe()
  single_pelist(1) = current_pe
  call mpp_set_current_pelist(single_pelist)

  do n=1, OBC%number_of_segments
    segment => OBC%segment(n)

    write(segnam,"('OBC_SEGMENT_',i3.3,'_DATA')") n
    write(suffix,"('_segment_',i3.3)") n
    call get_param(PF, mdl, segnam, segstr)

    call parse_segment_data_str(trim(segstr), fields=fields, num_fields=num_fields)
    if (num_fields == 0) then
        print *,'num_fields = 0';cycle ! cycle to next segment
    endif

    allocate(segment%field(num_fields))

    if (segment%Flather) then
      if (num_fields < 3) call MOM_error(FATAL, &
           "MOM_open_boundary, initialize_segment_data: "//&
           "Need at least three inputs for Flather")
      segment%num_fields = num_fields  ! these are at least three input fields required for the Flather option
    endif

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
        segment%values_needed = .true. ! Indicates that i/o will be needed for this segment
        segment%field(m)%name = trim(fields(m))
        if (segment%field(m)%name == 'TEMP') &
           segment%temp_segment_data_exists=.true.
        if (segment%field(m)%name == 'SALT') &
           segment%salt_segment_data_exists=.true.
        filename = trim(inputdir)//trim(filename)
        fieldname = trim(fieldname)//trim(suffix)
        call field_size(filename,fieldname,siz,no_domain=.true.)
        if (siz(4) == 1) segment%values_needed = .false.
        if (segment%on_pe) then
          if (modulo(siz(1),2) == 0 .or. modulo(siz(2),2) == 0) then
            call MOM_error(FATAL,'segment data are not on the supergrid')
          endif
          siz2(1)=1
          if (siz(1)>1) then
            siz2(1)=(siz(1)-1)/2
          endif
          siz2(2)=1
          if (siz(2)>1) then
            siz2(2)=(siz(2)-1)/2
          endif
          siz2(3)=siz(3)

          if (segment%is_E_or_W) then
            allocate(segment%field(m)%buffer_src(IsdB:IedB,jsd:jed,siz2(3)))
          else
            allocate(segment%field(m)%buffer_src(isd:ied,JsdB:JedB,siz2(3)))
          endif
          segment%field(m)%buffer_src(:,:,:)=0.0
          segment%field(m)%fid = init_external_field(trim(filename),&
               trim(fieldname),ignore_axis_atts=.true.,threading=SINGLE_FILE)
          if (siz(3) > 1) then
            fieldname = 'dz_'//trim(fieldname)
            call field_size(filename,fieldname,siz,no_domain=.true.)
            if (segment%is_E_or_W) then
              allocate(segment%field(m)%dz_src(IsdB:IedB,jsd:jed,siz(3)))
            else
              allocate(segment%field(m)%dz_src(isd:ied,JsdB:JedB,siz(3)))
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
      endif
    enddo
  enddo

  call mpp_set_current_pelist(saved_pelist)

end subroutine initialize_segment_data

!< Define indices for segment and store in hor_index_type
!< using global segment bounds corresponding to q-points
subroutine setup_segment_indices(G, seg, Is_obc, Ie_obc, Js_obc, Je_obc)
  type(dyn_horgrid_type), intent(in) :: G !< grid type
  type(OBC_segment_type), intent(inout) :: seg  !< Open boundary segment
  integer, intent(in) :: Is_obc !< Q-point global i-index of start of segment
  integer, intent(in) :: Ie_obc !< Q-point global i-index of end of segment
  integer, intent(in) :: Js_obc !< Q-point global j-index of start of segment
  integer, intent(in) :: Je_obc !< Q-point global j-index of end of segment
  ! Local variables
  integer :: Isg,Ieg,Jsg,Jeg

!  if (.not. G%Domain%symmetric) call MOM_error(FATAL, "MOM_open_boundary.F90, setup_segment_indices: "//&
!                       "Need to compile in symmetric mode")

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
subroutine setup_u_point_obc(OBC, G, segment_str, l_seg)
  type(ocean_OBC_type),    pointer    :: OBC !< Open boundary control structure
  type(dyn_horgrid_type),  intent(in) :: G !< Ocean grid structure
  character(len=*),        intent(in) :: segment_str !< A string in form of "I=%,J=%:%,string"
  integer,                 intent(in) :: l_seg !< which segment is this?
  ! Local variables
  integer :: I_obc, Js_obc, Je_obc ! Position of segment in global index space
  integer :: j, a_loop
  character(len=32) :: action_str(5)

  ! This returns the global indices for the segment
  call parse_segment_str(G%ieg, G%jeg, segment_str, I_obc, Js_obc, Je_obc, action_str )

  call setup_segment_indices(G, OBC%segment(l_seg),I_obc,I_obc,Js_obc,Je_obc)

  I_obc = I_obc - G%idg_offset ! Convert to local tile indices on this tile
  Js_obc = Js_obc - G%jdg_offset ! Convert to local tile indices on this tile
  Je_obc = Je_obc - G%jdg_offset ! Convert to local tile indices on this tile

  ! Hack to extend segment by one point
  if (OBC%extend_segments) then
    if (Js_obc<Je_obc) then
      Js_obc = Js_obc - 1 ; Je_obc = Je_obc + 1
    else
      Js_obc = Js_obc + 1 ; Je_obc = Je_obc - 1
    endif
  endif

  if (Je_obc>Js_obc) then
    OBC%segment(l_seg)%direction = OBC_DIRECTION_E
  else if (Je_obc<Js_obc) then
    OBC%segment(l_seg)%direction = OBC_DIRECTION_W
    j=js_obc;js_obc=je_obc;je_obc=j
  endif

  OBC%segment(l_seg)%on_pe = .false.

  do a_loop = 1,5 ! up to 5 options available
    if (len_trim(action_str(a_loop)) == 0) then
      cycle
    elseif (trim(action_str(a_loop)) == 'FLATHER') then
      OBC%segment(l_seg)%Flather = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%Flather_u_BCs_exist_globally = .true.
      OBC%open_u_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'ORLANSKI') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%open_u_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'OBLIQUE') then
      OBC%segment(l_seg)%oblique = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%oblique_BCs_exist_globally = .true.
      OBC%open_u_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED') then
      OBC%segment(l_seg)%nudged = .true.
      OBC%segment(l_seg)%Tnudge_in = 1.0/(3*86400)
      OBC%segment(l_seg)%Tnudge_out = 1.0/(360*86400)
      OBC%nudged_u_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'GRADIENT') then
      OBC%segment(l_seg)%gradient = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%open_u_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'LEGACY') then
      OBC%segment(l_seg)%Flather = .true.
      OBC%segment(l_seg)%radiation = .true.
      OBC%Flather_u_BCs_exist_globally = .true.
      OBC%open_u_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE') then
      OBC%segment(l_seg)%specified = .true.
      OBC%specified_u_BCs_exist_globally = .true. ! This avoids deallocation
      ! Hack to undo the hack above for SIMPLE BCs
      if (OBC%extend_segments) then
        Js_obc = Js_obc + 1
        Je_obc = Je_obc - 1
      endif
    else
      call MOM_error(FATAL, "MOM_open_boundary.F90, setup_u_point_obc: "//&
                     "String '"//trim(action_str(a_loop))//"' not understood.")
    endif

    if (I_obc<=G%HI%IsdB .or. I_obc>=G%HI%IedB) return ! Boundary is not on tile
    if (Je_obc<=G%HI%JsdB .or. Js_obc>=G%HI%JedB) return ! Segment is not on tile
  enddo ! a_loop

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

end subroutine setup_u_point_obc

!> Parse an OBC_SEGMENT_%%% string starting with "J=" and configure placement and type of OBC accordingly
subroutine setup_v_point_obc(OBC, G, segment_str, l_seg)
  type(ocean_OBC_type),    pointer    :: OBC !< Open boundary control structure
  type(dyn_horgrid_type),  intent(in) :: G !< Ocean grid structure
  character(len=*),        intent(in) :: segment_str !< A string in form of "J=%,I=%:%,string"
  integer,                 intent(in) :: l_seg !< which segment is this?
  ! Local variables
  integer :: J_obc, Is_obc, Ie_obc ! Position of segment in global index space
  integer :: i, a_loop
  character(len=32) :: action_str(5)

  ! This returns the global indices for the segment
  call parse_segment_str(G%ieg, G%jeg, segment_str, J_obc, Is_obc, Ie_obc, action_str )

  call setup_segment_indices(G, OBC%segment(l_seg),Is_obc,Ie_obc,J_obc,J_obc)

  J_obc = J_obc - G%jdg_offset ! Convert to local tile indices on this tile
  Is_obc = Is_obc - G%idg_offset ! Convert to local tile indices on this tile
  Ie_obc = Ie_obc - G%idg_offset ! Convert to local tile indices on this tile

  ! Hack to extend segment by one point
  if (OBC%extend_segments) then
    if (Is_obc<Ie_obc) then
      Is_obc = Is_obc - 1 ; Ie_obc = Ie_obc + 1
    else
      Is_obc = Is_obc + 1 ; Ie_obc = Ie_obc - 1
    endif
  endif

  if (Ie_obc>Is_obc) then
     OBC%segment(l_seg)%direction = OBC_DIRECTION_S
  else if (Ie_obc<Is_obc) then
     OBC%segment(l_seg)%direction = OBC_DIRECTION_N
     i=Is_obc;Is_obc=Ie_obc;Ie_obc=i
  endif

  OBC%segment(l_seg)%on_pe = .false.

  do a_loop = 1,5
    if (len_trim(action_str(a_loop)) == 0) then
      cycle
    elseif (trim(action_str(a_loop)) == 'FLATHER') then
      OBC%segment(l_seg)%Flather = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%Flather_v_BCs_exist_globally = .true.
      OBC%open_v_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'ORLANSKI') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%open_v_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'OBLIQUE') then
      OBC%segment(l_seg)%oblique = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%oblique_BCs_exist_globally = .true.
      OBC%open_v_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED') then
      OBC%segment(l_seg)%nudged = .true.
      OBC%segment(l_seg)%Tnudge_in = 1.0/(3*86400)
      OBC%segment(l_seg)%Tnudge_out = 1.0/(360*86400)
      OBC%nudged_v_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'GRADIENT') then
      OBC%segment(l_seg)%gradient = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%open_v_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'LEGACY') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%segment(l_seg)%Flather = .true.
      OBC%Flather_v_BCs_exist_globally = .true.
      OBC%open_v_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE') then
      OBC%segment(l_seg)%specified = .true.
      OBC%specified_v_BCs_exist_globally = .true. ! This avoids deallocation
      ! Hack to undo the hack above for SIMPLE BCs
      if (OBC%extend_segments) then
        Is_obc = Is_obc + 1
        Ie_obc = Ie_obc - 1
      endif
    else
      call MOM_error(FATAL, "MOM_open_boundary.F90, setup_v_point_obc: "//&
                     "String '"//trim(action_str(a_loop))//"' not understood.")
    endif

    if (J_obc<=G%HI%JsdB .or. J_obc>=G%HI%JedB) return ! Boundary is not on tile
    if (Ie_obc<=G%HI%IsdB .or. Is_obc>=G%HI%IedB) return ! Segment is not on tile
  enddo ! a_loop

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

end subroutine setup_v_point_obc

!> Parse an OBC_SEGMENT_%%% string
subroutine parse_segment_str(ni_global, nj_global, segment_str, l, m, n, action_str )
  integer,          intent(in)  :: ni_global !< Number of h-points in zonal direction
  integer,          intent(in)  :: nj_global !< Number of h-points in meridional direction
  character(len=*), intent(in)  :: segment_str !< A string in form of "I=l,J=m:n,string" or "J=l,I=m,n,string"
  integer,          intent(out) :: l !< The value of I=l, if segment_str begins with I=l, or the value of J=l
  integer,          intent(out) :: m !< The value of J=m, if segment_str begins with I=, or the value of I=m
  integer,          intent(out) :: n !< The value of J=n, if segment_str begins with I=, or the value of I=n
  character(len=*), intent(out) :: action_str(:) !< The "string" part of segment_str
  ! Local variables
  character(len=24) :: word1, word2, m_word, n_word !< Words delineated by commas in a string in form of "I=%,J=%:%,string"
  integer :: l_max !< Either ni_global or nj_global, depending on whether segment_str begins with "I=" or "J="
  integer :: mn_max !< Either nj_global or ni_global, depending on whether segment_str begins with "I=" or "J="
  integer :: j

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
  if (m<-1 .or. m>mn_max+1) then
    call MOM_error(FATAL, "MOM_open_boundary.F90, parse_segment_str: "//&
                   "Beginning of range in string '"//trim(segment_str)//"' is outside of the physical domain.")
  endif

  ! Read m
  n_word = extract_word(word2(3:24),':',2)
  n = interpret_int_expr( n_word, mn_max )
  if (n<-1 .or. n>mn_max+1) then
    call MOM_error(FATAL, "MOM_open_boundary.F90, parse_segment_str: "//&
                   "End of range in string '"//trim(segment_str)//"' is outside of the physical domain.")
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

  ! Returns integer value interpreted from string in form of %I, N or N-%I
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
      read(string(2:slen),*,err=911) interpret_int_expr
      interpret_int_expr = imax - interpret_int_expr
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
   character(len=*), intent(in)             :: segment_str !< A string in form of "VAR1=file:foo1.nc(varnam1),VAR2=file:foo2.nc(varnam2),..."
   character(len=*), intent(in),  optional  :: var         !< The name of the variable for which parameters are needed
   character(len=*), intent(out), optional  :: filenam     !< The name of the input file if using "file" method
   character(len=*), intent(out), optional  :: fieldnam    !< The name of the variable in the input file if using "file" method
   real,             intent(out), optional  :: value       !< A constant value if using the "value" method
   character(len=*), dimension(MAX_OBC_FIELDS), intent(out), optional :: fields   !< List of fieldnames for each segment
   integer, intent(out), optional           :: num_fields
   logical, intent(in), optional            :: debug
   ! Local variables
   character(len=128) :: word1, word2, word3, method
   integer :: lword, nfields, n, m, orient
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

!> Initialize open boundary control structure
subroutine open_boundary_init(G, param_file, OBC)
  type(ocean_grid_type), intent(in)    :: G !< Ocean grid structure
  type(param_file_type), intent(in)    :: param_file !< Parameter file handle
  type(ocean_OBC_type),  pointer       :: OBC !< Open boundary control structure
  ! Local variables

  if (.not.associated(OBC)) return

  id_clock_pass = cpu_clock_id('(Ocean OBC halo updates)', grain=CLOCK_ROUTINE)

end subroutine open_boundary_init

logical function open_boundary_query(OBC, apply_open_OBC, apply_specified_OBC, apply_Flather_OBC, apply_nudged_OBC, needs_ext_seg_data)
  type(ocean_OBC_type), pointer     :: OBC !< Open boundary control structure
  logical, optional,    intent(in)  :: apply_open_OBC      !< If present, returns True if specified_*_BCs_exist_globally is true
  logical, optional,    intent(in)  :: apply_specified_OBC !< If present, returns True if specified_*_BCs_exist_globally is true
  logical, optional,    intent(in)  :: apply_Flather_OBC   !< If present, returns True if Flather_*_BCs_exist_globally is true
  logical, optional,    intent(in)  :: apply_nudged_OBC    !< If present, returns True if nudged_*_BCs_exist_globally is true
  logical, optional,    intent(in)  :: needs_ext_seg_data  !< If present, returns True if external segment data needed
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
  type(OBC_segment_type), pointer :: segment
  integer :: n

  if (.not. associated(OBC)) return

  do n=1, OBC%number_of_segments
    segment => OBC%segment(n)
    call deallocate_OBC_segment_data(OBC, segment)
  enddo
  if (associated(OBC%segment)) deallocate(OBC%segment)
  if (associated(OBC%segnum_u)) deallocate(OBC%segnum_u)
  if (associated(OBC%segnum_v)) deallocate(OBC%segnum_v)
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
  type(OBC_segment_type), pointer :: segment

  if (.not.associated(OBC)) return

  if (.not.(OBC%open_u_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally)) &
    return

  do n=1,OBC%number_of_segments
    segment=>OBC%segment(n)
    if (.not. segment%on_pe .or. segment%specified) cycle
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
!! Also adjust u- and v-point cell area on specified open boundaries.
subroutine open_boundary_impose_land_mask(OBC, G, areaCu, areaCv)
  type(ocean_OBC_type),              pointer       :: OBC !< Open boundary control structure
  type(dyn_horgrid_type),            intent(inout) :: G !< Ocean grid structure
  real, dimension(SZIB_(G),SZJ_(G)), intent(inout) :: areaCu !< Area of a u-cell (m2)
  real, dimension(SZI_(G),SZJB_(G)), intent(inout) :: areaCv !< Area of a u-cell (m2)
  ! Local variables
  integer :: i, j, n
  type(OBC_segment_type), pointer :: segment
  logical :: any_U, any_V

  if (.not.associated(OBC)) return

  do n=1,OBC%number_of_segments
    segment=>OBC%segment(n)
    if (.not. segment%on_pe .or. segment%specified) cycle
    if (segment%is_E_or_W) then
      ! Sweep along u-segments and delete the OBC for blocked points.
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
          areaCu(I,j) = G%areaT(i,j)
        else   ! West
          areaCu(I,j) = G%areaT(i+1,j)
        endif
      enddo
    else
      ! Sweep along v-segments and for %specified BC points reset the v-point area which was masked out
      J=segment%HI%JsdB
      do i=segment%HI%isd,segment%HI%ied
        if (segment%direction == OBC_DIRECTION_S) then
          areaCv(i,J) = G%areaT(i,j+1)
        else      ! North
          areaCu(i,J) = G%areaT(i,j)
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

!> Apply radiation conditions to 3D  u,v at open boundaries
subroutine radiation_open_bdry_conds(OBC, u_new, u_old, v_new, v_old, G, dt)
  type(ocean_grid_type),                     intent(inout) :: G !< Ocean grid structure
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary control structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: u_new !< On exit, new u values on open boundaries
                                                                    !! On entry, the old time-level v but
                                                                    !! including barotropic accelerations.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u_old !< Original unadjusted u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: v_new !< On exit, new v values on open boundaries.
                                                                    !! On entry, the old time-level v but
                                                                    !! including barotropic accelerations.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v_old !< Original unadjusted v
  real,                                      intent(in)    :: dt    !< Appropriate timestep
  ! Local variables
  real :: dhdt, dhdx, dhdy, gamma_u, gamma_h, gamma_v
  real :: cff, Cx, Cy, tau
  real :: rx_max, ry_max ! coefficients for radiation
  real :: rx_new, rx_avg ! coefficients for radiation
  real :: ry_new, ry_avg ! coefficients for radiation
  real, parameter :: eps = 1.0e-20
  type(OBC_segment_type), pointer :: segment
  integer :: i, j, k, is, ie, js, je, nz, n
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not.associated(OBC)) return

  if (.not.(OBC%open_u_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally)) &
    return

  !! Copy previously calculated phase velocity from global arrays into segments
  !! This is terribly inefficient and temporary solution for continuity across restarts
  !! and needs to be revisited in the future.
  do n=1,OBC%number_of_segments
    segment=>OBC%segment(n)
    if (.not. segment%on_pe) cycle
    if (segment%is_E_or_W .and. segment%radiation) then
      do k=1,G%ke
        I=segment%HI%IsdB
        do j=segment%HI%jsd,segment%HI%jed
          segment%rx_normal(I,j,k) = OBC%rx_normal(I,j,k)
        enddo
      enddo
    elseif (segment%is_N_or_S .and. segment%radiation) then
      do k=1,G%ke
        J=segment%HI%JsdB
        do i=segment%HI%isd,segment%HI%ied
          segment%rx_normal(i,J,k) = OBC%ry_normal(i,J,k)
        enddo
      enddo
    endif
  enddo

  gamma_u = OBC%gamma_uv ; gamma_v = OBC%gamma_uv ; gamma_h = OBC%gamma_h
  rx_max = OBC%rx_max ; ry_max = OBC%rx_max
  do n=1,OBC%number_of_segments
     segment=>OBC%segment(n)
     if (.not. segment%on_pe) cycle
     if (segment%oblique) call gradient_at_q_points(G,segment,u_new,v_new)
     if (segment%direction == OBC_DIRECTION_E) then
       I=segment%HI%IsdB
       if (I<G%HI%IscB) cycle
       do k=1,nz ;  do j=segment%HI%jsd,segment%HI%jed
         if (segment%radiation) then
           dhdt = u_old(I-1,j,k)-u_new(I-1,j,k) !old-new
           dhdx = u_new(I-1,j,k)-u_new(I-2,j,k) !in new time backward sasha for I-1
           rx_new = 0.0
           if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max) ! outward phase speed
           rx_avg = (1.0-gamma_u)*segment%rx_normal(I,j,k) + gamma_u*rx_new
           segment%rx_normal(I,j,k) = rx_avg
           ! The new boundary value is interpolated between future interior
           ! value, u_new(I-1) and past boundary value but with barotropic
           ! accelerations, u_new(I).
           segment%normal_vel(I,j,k) = (u_new(I,j,k) + rx_avg*u_new(I-1,j,k)) / (1.0+rx_avg)
           ! Copy restart fields into 3-d arrays. This is an inefficient and temporary issues
           ! implemented as a work-around to limitations in restart capability
           OBC%rx_normal(I,j,k) = segment%rx_normal(I,j,k)
         elseif (segment%oblique) then
           dhdt = u_old(I-1,j,k)-u_new(I-1,j,k) !old-new
           dhdx = u_new(I-1,j,k)-u_new(I-2,j,k) !in new time backward sasha for I-1
!          if (segment%oblique) then
             if (dhdt*(segment%grad_normal(J,1,k) + segment%grad_normal(J-1,1,k)) > 0.0) then
               dhdy = segment%grad_normal(J-1,1,k)
             elseif (dhdt*(segment%grad_normal(J,1,k) + segment%grad_normal(J-1,1,k)) == 0.0) then
               dhdy = 0.0
             else
               dhdy = segment%grad_normal(J,1,k)
             endif
!          endif
           if (dhdt*dhdx < 0.0) dhdt = 0.0
           Cx = min(dhdt*dhdx,rx_max) ! default to normal radiation
!          Cy = 0.0
           cff = max(dhdx*dhdx,eps)
!          if (segment%oblique) then
             cff = max(dhdx*dhdx + dhdy*dhdy, eps)
             Cy = min(cff,max(dhdt*dhdy,-cff))
!          endif
           segment%normal_vel(I,j,k) = ((cff*u_new(I,j,k) + Cx*u_new(I-1,j,k)) - &
              (max(Cy,0.0)*segment%grad_normal(J-1,2,k) + min(Cy,0.0)*segment%grad_normal(J,2,k))) / (cff + Cx)
         elseif (segment%gradient) then
           segment%normal_vel(I,j,k) = u_new(I-1,j,k)
         endif
         if ((segment%radiation .or. segment%oblique) .and. segment%nudged) then
           if (dhdt*dhdx < 0.0) then
             tau = segment%Tnudge_in
           else
             tau = segment%Tnudge_out
           endif
           segment%normal_vel(I,j,k) = u_new(I,j,k) + dt*tau*(segment%nudged_normal_vel(I,j,k) - u_new(I,j,k))
         endif
       enddo; enddo
     endif

     if (segment%direction == OBC_DIRECTION_W) then
       I=segment%HI%IsdB
       if (I>G%HI%IecB) cycle
       do k=1,nz ;  do j=segment%HI%jsd,segment%HI%jed
         if (segment%radiation) then
           dhdt = u_old(I+1,j,k)-u_new(I+1,j,k) !old-new
           dhdx = u_new(I+1,j,k)-u_new(I+2,j,k) !in new time forward sasha for I+1
           rx_new = 0.0
           if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
           rx_avg = (1.0-gamma_u)*segment%rx_normal(I,j,k) + gamma_u*rx_new
           segment%rx_normal(I,j,k) = rx_avg
           ! The new boundary value is interpolated between future interior
           ! value, u_new(I+1) and past boundary value but with barotropic
           ! accelerations, u_new(I).
           segment%normal_vel(I,j,k) = (u_new(I,j,k) + rx_avg*u_new(I+1,j,k)) / (1.0+rx_avg)
           ! Copy restart fields into 3-d arrays. This is an inefficient and temporary issues
           ! implemented as a work-around to limitations in restart capability
           OBC%rx_normal(I,j,k) = segment%rx_normal(I,j,k)
         elseif (segment%oblique) then
           dhdt = u_old(I+1,j,k)-u_new(I+1,j,k) !old-new
           dhdx = u_new(I+1,j,k)-u_new(I+2,j,k) !in new time forward sasha for I+1
!          if (segment%oblique) then
             if (dhdt*(segment%grad_normal(J,1,k) + segment%grad_normal(J-1,1,k)) > 0.0) then
               dhdy = segment%grad_normal(J-1,1,k)
             elseif (dhdt*(segment%grad_normal(J,1,k) + segment%grad_normal(J-1,1,k)) == 0.0) then
               dhdy = 0.0
             else
               dhdy = segment%grad_normal(J,1,k)
             endif
!          endif
           if (dhdt*dhdx < 0.0) dhdt = 0.0
           Cx = min(dhdt*dhdx,rx_max) ! default to normal flow only
!          Cy = 0.
           cff = max(dhdx*dhdx, eps)
!          if (segment%oblique) then
             cff = max(dhdx*dhdx + dhdy*dhdy, eps)
             Cy = min(cff,max(dhdt*dhdy,-cff))
!          endif
           segment%normal_vel(I,j,k) = ((cff*u_new(I,j,k) + Cx*u_new(I+1,j,k)) - &
             (max(Cy,0.0)*segment%grad_normal(J-1,2,k) + min(Cy,0.0)*segment%grad_normal(J,2,k))) / (cff + Cx)
         elseif (segment%gradient) then
           segment%normal_vel(I,j,k) = u_new(I+1,j,k)
         endif
         if ((segment%radiation .or. segment%oblique) .and. segment%nudged) then
           if (dhdt*dhdx < 0.0) then
             tau = segment%Tnudge_in
           else
             tau = segment%Tnudge_out
           endif
           segment%normal_vel(I,j,k) = u_new(I,j,k) + dt*tau*(segment%nudged_normal_vel(I,j,k) - u_new(I,j,k))
         endif
       enddo; enddo
     endif

     if (segment%direction == OBC_DIRECTION_N) then
       J=segment%HI%JsdB
       if (J<G%HI%JscB) cycle
       do k=1,nz ;  do i=segment%HI%isd,segment%HI%ied
         if (segment%radiation) then
           dhdt = v_old(i,J-1,k)-v_new(i,J-1,k) !old-new
           dhdy = v_new(i,J-1,k)-v_new(i,J-2,k) !in new time backward sasha for J-1
           ry_new = 0.0
           if (dhdt*dhdy > 0.0) ry_new = min( (dhdt/dhdy), ry_max)
           ry_avg = (1.0-gamma_v)*segment%rx_normal(I,j,k) + gamma_v*ry_new
           segment%rx_normal(i,J,k) = ry_avg
           ! The new boundary value is interpolated between future interior
           ! value, v_new(J-1) and past boundary value but with barotropic
           ! accelerations, v_new(J).
           segment%normal_vel(i,J,k) = (v_new(i,J,k) + ry_avg*v_new(i,J-1,k)) / (1.0+ry_avg)
           ! Copy restart fields into 3-d arrays. This is an inefficient and temporary issues
           ! implemented as a work-around to limitations in restart capability
           OBC%ry_normal(i,J,k) = segment%rx_normal(i,J,k)
         elseif (segment%oblique) then
           dhdt = v_old(i,J-1,k)-v_new(i,J-1,k) !old-new
           dhdy = v_new(i,J-1,k)-v_new(i,J-2,k) !in new time backward sasha for J-1
!          if (segment%oblique) then
             if (dhdt*(segment%grad_normal(I,1,k) + segment%grad_normal(I-1,1,k)) > 0.0) then
               dhdx = segment%grad_normal(I-1,1,k)
             elseif (dhdt*(segment%grad_normal(I,1,k) + segment%grad_normal(I-1,1,k)) == 0.0) then
               dhdx = 0.0
             else
               dhdx = segment%grad_normal(I,1,k)
             endif
!          endif
           if (dhdt*dhdy < 0.0) dhdt = 0.0
           Cy = min(dhdt*dhdy,rx_max) ! default to normal flow only
!          Cx = 0
           cff = max(dhdy*dhdy, eps)
!          if (segment%oblique) then
             cff = max(dhdx*dhdx + dhdy*dhdy, eps)
             Cx = min(cff,max(dhdt*dhdx,-cff))
!          endif
           segment%normal_vel(i,J,k) = ((cff*v_new(i,J,k) + Cy*v_new(i,J-1,k)) - &
              (max(Cx,0.0)*segment%grad_normal(I-1,2,k) + min(Cx,0.0)*segment%grad_normal(I,2,k))) / (cff + Cy)
         elseif (segment%gradient) then
           segment%normal_vel(i,J,k) = v_new(i,J-1,k)
         endif
         if ((segment%radiation .or. segment%oblique) .and. segment%nudged) then
           if (dhdt*dhdy < 0.0) then
             tau = segment%Tnudge_in
           else
             tau = segment%Tnudge_out
           endif
           segment%normal_vel(i,J,k) = v_new(i,J,k) + dt*tau*(segment%nudged_normal_vel(i,J,k) - v_new(i,J,k))
         endif
       enddo; enddo
     endif


     if (segment%direction == OBC_DIRECTION_S) then
       J=segment%HI%JsdB
       if (J>G%HI%JecB) cycle
       do k=1,nz ;  do i=segment%HI%isd,segment%HI%ied
         if (segment%radiation) then
           dhdt = v_old(i,J+1,k)-v_new(i,J+1,k) !old-new
           dhdy = v_new(i,J+1,k)-v_new(i,J+2,k) !in new time backward sasha for J-1
           ry_new = 0.0
           if (dhdt*dhdy > 0.0) ry_new = min( (dhdt/dhdy), ry_max)
           ry_avg = (1.0-gamma_v)*segment%rx_normal(I,j,k) + gamma_v*ry_new
           segment%rx_normal(i,J,k) = ry_avg
           ! The new boundary value is interpolated between future interior
           ! value, v_new(J+1) and past boundary value but with barotropic
           ! accelerations, v_new(J).
           segment%normal_vel(i,J,k) = (v_new(i,J,k) + ry_avg*v_new(i,J+1,k)) / (1.0+ry_avg)
           ! Copy restart fields into 3-d arrays. This is an inefficient and temporary issues
           ! implemented as a work-around to limitations in restart capability
           OBC%ry_normal(i,J,k) = segment%rx_normal(i,J,k)
         elseif (segment%oblique) then
           dhdt = v_old(i,J+1,k)-v_new(i,J+1,k) !old-new
           dhdy = v_new(i,J+1,k)-v_new(i,J+2,k) !in new time backward sasha for J-1
!          if (segment%oblique) then
             if (dhdt*(segment%grad_normal(I,1,k) + segment%grad_normal(I-1,1,k)) > 0.0) then
               dhdx = segment%grad_normal(I-1,1,k)
             elseif (dhdt*(segment%grad_normal(I,1,k) + segment%grad_normal(I-1,1,k)) == 0.0) then
               dhdx = 0.0
             else
               dhdx = segment%grad_normal(I,1,k)
             endif
!          endif
           if (dhdt*dhdy < 0.0) dhdt = 0.0
           Cy = min(dhdt*dhdy,rx_max) ! default to normal flow only
!          Cx = 0
           cff = max(dhdy*dhdy, eps)
!          if (segment%oblique) then
             cff = max(dhdx*dhdx + dhdy*dhdy, eps)
             Cx = min(cff,max(dhdt*dhdx,-cff))
!          endif
           segment%normal_vel(i,J,k) = ((cff*v_new(i,J,k) + Cy*v_new(i,J+1,k)) - &
              (max(Cx,0.0)*segment%grad_normal(I-1,2,k) + min(Cx,0.0)*segment%grad_normal(I,2,k))) / (cff + Cy)
         elseif (segment%gradient) then
           segment%normal_vel(i,J,k) = v_new(i,J+1,k)
         endif
         if ((segment%radiation .or. segment%oblique) .and. segment%nudged) then
           if (dhdt*dhdy < 0.0) then
             tau = segment%Tnudge_in
           else
             tau = segment%Tnudge_out
           endif
           segment%normal_vel(i,J,k) = v_new(i,J,k) + dt*tau*(segment%nudged_normal_vel(i,J,k) - v_new(i,J,k))
         endif
       enddo; enddo
     end if
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
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: u   !< u field to update on open boundaries
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: v   !< v field to update on open boundaries
  ! Local variables
  integer :: i, j, k, n
  type(OBC_segment_type), pointer :: segment

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
        enddo; enddo
      elseif (segment%is_N_or_S) then
        J=segment%HI%JsdB
        do k=1,G%ke ;  do i=segment%HI%isd,segment%HI%ied
          v(i,J,k) = segment%normal_vel(i,J,k)
        enddo; enddo
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
  type(OBC_segment_type), pointer :: segment

  if (.not.associated(OBC)) return ! Bail out if OBC is not available

  do n=1,OBC%number_of_segments
    segment => OBC%segment(n)
    if (.not. segment%on_pe) then
      cycle
    elseif (segment%is_E_or_W) then
      I=segment%HI%IsdB
      do k=1,G%ke ;  do j=segment%HI%jsd,segment%HI%jed
        u(I,j,k) = 0.
      enddo; enddo
    elseif (segment%is_N_or_S) then
      J=segment%HI%JsdB
      do k=1,G%ke ;  do i=segment%HI%isd,segment%HI%ied
        v(i,J,k) = 0.
      enddo; enddo
    endif
  enddo

end subroutine open_boundary_zero_normal_flow

!> Calculate the tangential gradient of the normal flow at the boundary q-points.
subroutine gradient_at_q_points(G,segment,uvel,vvel)
  type(ocean_grid_type), intent(in) :: G !< Ocean grid structure
  type(OBC_segment_type), pointer :: segment !< OBC segment structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: uvel !< zonal velocity
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: vvel !< meridional velocity
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
    else ! western segment
      I=segment%HI%isdB
      do k=1,G%ke
        do J=max(segment%HI%JsdB, G%HI%JsdB+1),min(segment%HI%JedB, G%HI%JedB-1)
          segment%grad_normal(J,1,k) = (uvel(I+1,j+1,k)-uvel(I+1,j,k)) * G%mask2dBu(I+1,J)
          segment%grad_normal(J,2,k) = (uvel(I,j+1,k)-uvel(I,j,k)) * G%mask2dBu(I,J)
        enddo
      enddo
    endif
  else if (segment%is_N_or_S) then
    if (segment%direction == OBC_DIRECTION_N) then
      J=segment%HI%jsdB
      do k=1,G%ke
        do I=max(segment%HI%IsdB, G%HI%IsdB+1),min(segment%HI%IedB, G%HI%IedB-1)
          segment%grad_normal(I,1,k) = (vvel(i+1,J-1,k)-vvel(i,J-1,k)) * G%mask2dBu(I,J-1)
          segment%grad_normal(I,2,k) = (vvel(i+1,J,k)-vvel(i,J,k)) * G%mask2dBu(I,J)
        enddo
      enddo
    else ! south segment
      J=segment%HI%jsdB
      do k=1,G%ke
        do I=max(segment%HI%IsdB, G%HI%IsdB+1),min(segment%HI%IedB, G%HI%IedB-1)
          segment%grad_normal(I,1,k) = (vvel(i+1,J+1,k)-vvel(i,J+1,k)) * G%mask2dBu(I,J+1)
          segment%grad_normal(I,2,k) = (vvel(i+1,J,k)-vvel(i,J,k)) * G%mask2dBu(I,J)
        enddo
      enddo
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
  type(OBC_segment_type), pointer :: segment ! pointer to segment type list
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
        enddo; enddo
      elseif (segment%direction == OBC_DIRECTION_W) then
        I=segment%HI%IsdB
        do k=1,G%ke ;  do j=segment%HI%jsd,segment%HI%jed
          tv%T(i,j,k) = tv%T(i+1,j,k) ; tv%S(i,j,k) = tv%S(i+1,j,k)
        enddo; enddo
      elseif (segment%direction == OBC_DIRECTION_N) then
        J=segment%HI%JsdB
        do k=1,G%ke ;  do i=segment%HI%isd,segment%HI%ied
          tv%T(i,j+1,k) = tv%T(i,j,k) ; tv%S(i,j+1,k) = tv%S(i,j,k)
        enddo; enddo
      elseif (segment%direction == OBC_DIRECTION_S) then
        J=segment%HI%JsdB
        do k=1,G%ke ;  do i=segment%HI%isd,segment%HI%ied
          tv%T(i,j,k) = tv%T(i,j+1,k) ; tv%S(i,j,k) = tv%S(i,j+1,k)
        enddo; enddo
      endif
    enddo
  endif

! do n=1,OBC%number_of_segments
!   segment => OBC%segment(n)
!   if (.not. segment%on_pe) cycle

!   if (segment%direction == OBC_DIRECTION_E) then
!     I=segment%HI%IsdB
!     do k=1,G%ke ;  do j=segment%HI%jsd,segment%HI%jed
!       h(i+1,j,k) = h(i,j,k)
!     enddo; enddo
!   elseif (segment%direction == OBC_DIRECTION_W) then
!     I=segment%HI%IsdB
!     do k=1,G%ke ;  do j=segment%HI%jsd,segment%HI%jed
!       h(i,j,k) = h(i+1,j,k)
!     enddo; enddo
!   elseif (segment%direction == OBC_DIRECTION_N) then
!     J=segment%HI%JsdB
!     do k=1,G%ke ;  do i=segment%HI%isd,segment%HI%ied
!       h(i,j+1,k) = h(i,j,k)
!     enddo; enddo
!   elseif (segment%direction == OBC_DIRECTION_S) then
!     J=segment%HI%JsdB
!     do k=1,G%ke ;  do i=segment%HI%isd,segment%HI%ied
!       h(i,j,k) = h(i,j+1,k)
!     enddo; enddo
!   endif
! enddo

end subroutine set_tracer_data

function lookup_seg_field(OBC_seg,field)
  type(OBC_segment_type), pointer :: OBC_seg
  character(len=32), intent(in) :: field ! The field name
  integer :: lookup_seg_field

  integer :: n,m

  lookup_seg_field=-1
  do n=1,OBC_seg%num_fields
   if (trim(field) == OBC_seg%field(m)%name) then
     lookup_seg_field=n
     return
   endif
  enddo

  return

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
    allocate(segment%normal_trans_bt(IsdB:IedB,jsd:jed));     segment%normal_trans_bt(:,:)=0.0
    if (segment%radiation) then
      allocate(segment%rx_normal(IsdB:IedB,jsd:jed,OBC%ke));    segment%rx_normal(:,:,:)=0.0
    endif
    allocate(segment%normal_vel(IsdB:IedB,jsd:jed,OBC%ke));   segment%normal_vel(:,:,:)=0.0
    allocate(segment%normal_vel_bt(IsdB:IedB,jsd:jed));       segment%normal_vel_bt(:,:)=0.0
    allocate(segment%normal_trans(IsdB:IedB,jsd:jed,OBC%ke)); segment%normal_trans(:,:,:)=0.0
    if (segment%nudged) then
      allocate(segment%nudged_normal_vel(IsdB:IedB,jsd:jed,OBC%ke)); segment%nudged_normal_vel(:,:,:)=0.0
    endif
    if (segment%oblique) then
      allocate(segment%grad_normal(JsdB:JedB,2,OBC%ke)); segment%grad_normal(:,:,:) = 0.0
    endif
  endif

  if (segment%is_N_or_S) then
    ! If these are just Flather, change update_OBC_segment_data accordingly
    allocate(segment%Cg(isd:ied,JsdB:JedB));                  segment%Cg(:,:)=0.
    allocate(segment%Htot(isd:ied,JsdB:JedB));                segment%Htot(:,:)=0.0
    allocate(segment%h(isd:ied,JsdB:JedB,OBC%ke));            segment%h(:,:,:)=0.0
    allocate(segment%eta(isd:ied,JsdB:JedB));                 segment%eta(:,:)=0.0
    allocate(segment%normal_trans_bt(isd:ied,JsdB:JedB));     segment%normal_trans_bt(:,:)=0.0
    if (segment%radiation) then
      allocate(segment%rx_normal(isd:ied,JsdB:JedB,OBC%ke));    segment%rx_normal(:,:,:)=0.0
    endif
    allocate(segment%normal_vel(isd:ied,JsdB:JedB,OBC%ke));   segment%normal_vel(:,:,:)=0.0
    allocate(segment%normal_vel_bt(isd:ied,JsdB:JedB));       segment%normal_vel_bt(:,:)=0.0
    allocate(segment%normal_trans(isd:ied,JsdB:JedB,OBC%ke)); segment%normal_trans(:,:,:)=0.0
    if (segment%nudged) then
      allocate(segment%nudged_normal_vel(isd:ied,JsdB:JedB,OBC%ke)); segment%nudged_normal_vel(:,:,:)=0.0
    endif
    if (segment%oblique) then
      allocate(segment%grad_normal(IsdB:IedB,2,OBC%ke)); segment%grad_normal(:,:,:) = 0.0
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
  if (associated (segment%normal_trans_bt)) deallocate(segment%normal_trans_bt)
  if (associated (segment%rx_normal)) deallocate(segment%rx_normal)
  if (associated (segment%normal_vel)) deallocate(segment%normal_vel)
  if (associated (segment%normal_vel_bt)) deallocate(segment%normal_vel_bt)
  if (associated (segment%normal_trans)) deallocate(segment%normal_trans)
  if (associated (segment%nudged_normal_vel)) deallocate(segment%nudged_normal_vel)
  if (associated (segment%tr_Reg)) call segment_tracer_registry_end(segment%tr_Reg)


end subroutine deallocate_OBC_segment_data

!> Set tangential velocities outside of open boundaries to silly values
!! (used for checking the interior state is independent of values outside
!! of the domain).
subroutine open_boundary_test_extern_uv(G, OBC, u, v)
  type(ocean_grid_type),                     intent(in)    :: G !< Ocean grid structure
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary structure
  real, dimension(SZIB_(G),SZJ_(G), SZK_(G)),intent(inout) :: u !< Zonal velocity (m/s)
  real, dimension(SZI_(G),SZJB_(G), SZK_(G)),intent(inout) :: v !< Meridional velocity (m/s)
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
subroutine open_boundary_test_extern_h(G, OBC, h)
  type(ocean_grid_type),                    intent(in)    :: G   !< Ocean grid structure
  type(ocean_OBC_type),                     pointer       :: OBC !< Open boundary structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)),intent(inout) :: h   !< Layer thickness (m or kg/m2)
  ! Local variables
  integer :: i, j, k, n

  if (.not. associated(OBC)) return

  do n = 1, OBC%number_of_segments
    do k = 1, G%ke
      if (OBC%segment(n)%is_N_or_S) then
        J = OBC%segment(n)%HI%JsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
          do i = OBC%segment(n)%HI%isd, OBC%segment(n)%HI%ied
            h(i,j+1,k) = OBC%silly_h
          enddo
        else
          do i = OBC%segment(n)%HI%isd, OBC%segment(n)%HI%ied
            h(i,j,k) = OBC%silly_h
          enddo
        endif
      elseif (OBC%segment(n)%is_E_or_W) then
        I = OBC%segment(n)%HI%IsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
          do j = OBC%segment(n)%HI%jsd, OBC%segment(n)%HI%jed
            h(i+1,j,k) = OBC%silly_h
          enddo
        else
          do j = OBC%segment(n)%HI%jsd, OBC%segment(n)%HI%jed
            h(i,j,k) = OBC%silly_h
          enddo
        endif
      endif
    enddo
  enddo

end subroutine open_boundary_test_extern_h

!> Update the OBC values on the segments.
subroutine update_OBC_segment_data(G, GV, OBC, tv, h, Time)
  type(ocean_grid_type),                     intent(in)    :: G    !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV   !<  Ocean vertical grid structure
  type(ocean_OBC_type),                      pointer       :: OBC  !< Open boundary structure
  type(thermo_var_ptrs),                     intent(in)    :: tv   !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(inout) :: h    !< Thickness
! real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(inout) :: e    !< Layer interface height
! real, dimension(SZI_(G),SZJ_(G))         , intent(inout) :: eta  !< Thickness
  type(time_type),                           intent(in)    :: Time !< Time
  ! Local variables

  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB, n, m, nz
  character(len=40)  :: mdl = "set_OBC_segment_data" ! This subroutine's name.
  character(len=200) :: filename, OBC_file, inputdir ! Strings for file/path
  type(OBC_segment_type), pointer :: segment
  integer, dimension(4) :: siz,siz2
  real :: sumh ! column sum of thicknesses (m)
  integer :: ni_seg, nj_seg  ! number of src gridpoints along the segments
  integer :: i2, j2          ! indices for referencing local domain array
  integer :: is_obc, ie_obc, js_obc, je_obc  ! segment indices within local domain
  integer :: ishift, jshift  ! offsets for staggered locations
  real, dimension(:,:), pointer :: seg_vel => NULL()  ! pointer to segment velocity array
  real, dimension(:,:), pointer :: seg_trans => NULL()  ! pointer to segment transport array
  real, dimension(:,:,:), allocatable :: tmp_buffer

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


    if (segment%is_E_or_W) then
      nj_seg=nj_seg-1
      js_obc=js_obc+1
    else
      ni_seg=ni_seg-1
      is_obc=is_obc+1
    endif

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
      if (segment%direction == OBC_DIRECTION_W) ishift=1
      I=segment%HI%IsdB
      do j=segment%HI%jsd,segment%HI%jed
        segment%Cg(I,j) = sqrt(GV%g_prime(1)*G%bathyT(i+ishift,j))
 !       if (GV%Boussinesq) then
        segment%Htot(I,j) = G%bathyT(i+ishift,j)*GV%m_to_H! + eta(i+ishift,j)
 !       else
 !         segment%Htot(I,j) =  eta(i+ishift,j)
 !       endif
        do k=1,G%ke
          segment%h(I,j,k) = h(i+ishift,j,k)
        enddo
      enddo


    else! (segment%direction == OBC_DIRECTION_N .or. segment%direction == OBC_DIRECTION_S)
      if (segment%direction == OBC_DIRECTION_S) jshift=1
      J=segment%HI%JsdB
      do i=segment%HI%isd,segment%HI%ied
        segment%Cg(i,J) = sqrt(GV%g_prime(1)*G%bathyT(i,j+jshift))
!       if (GV%Boussinesq) then
        segment%Htot(i,J) = G%bathyT(i,j+jshift)*GV%m_to_H! + eta(i,j+jshift)
!          else
!            segment%Htot(i,J) = eta(i,j+jshift)
!          endif
        do k=1,G%ke
          segment%h(i,J,k) = h(i,j+jshift,k)
!           segment%e(i,J,k) = e(i,j+jshift,k)
        enddo
      enddo
    endif

    do m = 1,segment%num_fields
      if (segment%field(m)%fid > 0) then
        siz(1)=size(segment%field(m)%buffer_src,1)
        siz(2)=size(segment%field(m)%buffer_src,2)
        siz(3)=size(segment%field(m)%buffer_src,3)
        if (.not.associated(segment%field(m)%buffer_dst)) then
          if (siz(3) /= segment%field(m)%nk_src) call MOM_error(FATAL,'nk_src inconsistency')
          if (segment%field(m)%nk_src > 1) then
            allocate(segment%field(m)%buffer_dst(is_obc:ie_obc,js_obc:je_obc,G%ke))
          else
            allocate(segment%field(m)%buffer_dst(is_obc:ie_obc,js_obc:je_obc,1))
          endif
          segment%field(m)%buffer_dst(:,:,:)=0.0
          if (trim(segment%field(m)%name) == 'U' .or. trim(segment%field(m)%name) == 'V') then
            allocate(segment%field(m)%bt_vel(is_obc:ie_obc,js_obc:je_obc))
            segment%field(m)%bt_vel(:,:)=0.0
          endif
        endif
        ! read source data interpolated to the current model time
        if (siz(1)==1) then
          allocate(tmp_buffer(1,nj_seg*2+1,segment%field(m)%nk_src))  ! segment data is currrently on supergrid
        else
          allocate(tmp_buffer(ni_seg*2+1,1,segment%field(m)%nk_src))  ! segment data is currrently on supergrid
        endif

        call time_interp_external(segment%field(m)%fid,Time, tmp_buffer)
        if (siz(1)==1) then
          segment%field(m)%buffer_src(is_obc,:,:)=tmp_buffer(1,2*(js_obc+G%jdg_offset)-1:2*(je_obc+G%jdg_offset)-1:2,:)
        else
          segment%field(m)%buffer_src(:,js_obc,:)=tmp_buffer(2*(is_obc+G%idg_offset)-1:2*(ie_obc+G%idg_offset)-1:2,1,:)
        endif
        if (segment%field(m)%nk_src > 1) then
          call time_interp_external(segment%field(m)%fid_dz,Time, tmp_buffer)
          if (siz(1)==1) then
            segment%field(m)%dz_src(is_obc,:,:)=tmp_buffer(1,2*(js_obc+G%jdg_offset)-1:2*(je_obc+G%jdg_offset)-1:2,:)
          else
            segment%field(m)%dz_src(:,js_obc,:)=tmp_buffer(2*(is_obc+G%idg_offset)-1:2*(ie_obc+G%idg_offset)-1:2,1,:)
          endif
          do j=js_obc,je_obc
            do i=is_obc,ie_obc

              ! Using the h remapping approach
              ! Pretty sure we need to check for source/target grid consistency here
              segment%field(m)%buffer_dst(i,j,:)=0.0  ! initialize remap destination buffer
              if (G%mask2dT(i,j)>0.) then
                call remapping_core_h(OBC%remap_CS, &
                     segment%field(m)%nk_src,segment%field(m)%dz_src(i,j,:), &
                     segment%field(m)%buffer_src(i,j,:), &
                     G%ke, h(i,j,:), segment%field(m)%buffer_dst(i,j,:))
              endif
            enddo
          enddo
        else  ! 2d data
          segment%field(m)%buffer_dst(:,:,1)=segment%field(m)%buffer_src(:,:,1)  ! initialize remap destination buffer
        endif
        deallocate(tmp_buffer)
      else ! fid <= 0
        if (.not. ASSOCIATED(segment%field(m)%buffer_dst)) then
          allocate(segment%field(m)%buffer_dst(is_obc:ie_obc,js_obc:je_obc,G%ke))
          segment%field(m)%buffer_dst(:,:,:)=segment%field(m)%value
          if (trim(segment%field(m)%name) == 'U' .or. trim(segment%field(m)%name) == 'V') then
            allocate(segment%field(m)%bt_vel(is_obc:ie_obc,js_obc:je_obc))
            segment%field(m)%bt_vel(:,:)=segment%field(m)%value
          endif
        endif
      endif

      if (trim(segment%field(m)%name) == 'U' .or. trim(segment%field(m)%name) == 'V') then
        if (segment%field(m)%fid>0) then ! calculate external BT velocity and transport if needed
          if((trim(segment%field(m)%name) == 'U' .and. segment%is_E_or_W) .or.  &
             (trim(segment%field(m)%name) == 'V' .and. segment%is_N_or_S)) then
            do j=js_obc,je_obc
              do i=is_obc,ie_obc
                segment%normal_trans_bt(i,j) = 0.0
                do k=1,G%ke
                  segment%normal_vel(i,j,k) = segment%field(m)%buffer_dst(i,j,k)
                  segment%normal_trans(i,j,k) = segment%field(m)%buffer_dst(i,j,k)*segment%h(i,j,k)
                  segment%normal_trans_bt(i,j)= segment%normal_trans_bt(i,j)+segment%normal_trans(i,j,k)
                enddo
                segment%normal_vel_bt(i,j) = segment%normal_trans_bt(i,j)/max(segment%Htot(i,j),1.e-12)
                if (associated(segment%nudged_normal_vel)) segment%nudged_normal_vel(i,j,:) = segment%normal_vel(i,j,:)
              enddo
            enddo
          endif
        endif
      endif

      if (trim(segment%field(m)%name) == 'SSH') then
        do j=js_obc,je_obc
          do i=is_obc,ie_obc
            segment%eta(i,j) = segment%field(m)%buffer_dst(i,j,1)
          enddo
        enddo
      endif

      if (trim(segment%field(m)%name) == 'TEMP') then
        if (associated(segment%field(m)%buffer_dst)) then
          do k=1,nz; do j=js_obc, je_obc;do i=is_obc,ie_obc
            segment%tr_Reg%Tr(1)%t(i,j,k) = segment%field(m)%buffer_dst(i,j,k)
          enddo; enddo; enddo
        else
          segment%tr_Reg%Tr(1)%OBC_inflow_conc = segment%field(m)%value
        endif
      elseif (trim(segment%field(m)%name) == 'SALT') then
        if (associated(segment%field(m)%buffer_dst)) then
          do k=1,nz; do j=js_obc, je_obc;do i=is_obc,ie_obc
            segment%tr_Reg%Tr(2)%t(i,j,k) = segment%field(m)%buffer_dst(i,j,k)
          enddo; enddo; enddo
        else
          segment%tr_Reg%Tr(2)%OBC_inflow_conc = segment%field(m)%value
        endif
      endif

    enddo

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

subroutine register_segment_tracer(tr_desc, param_file, GV, segment, tr_desc_ptr, &
                                   OBC_scalar, OBC_array)
  type(verticalGrid_type),        intent(in)    :: GV           !< ocean vertical grid structure
  type(vardesc),         intent(in)             :: tr_desc      !< metadata about the tracer
  type(param_file_type), intent(in)             :: param_file   !< file to parse for  model parameter values
  type(OBC_segment_type), intent(inout)            :: segment      !< current segment data structure
  type(vardesc), target, optional               :: tr_desc_ptr  !< A target that can be used to set a pointer to the
                                                                !! stored value of tr%tr_desc. This target must be
                                                                !! an enduring part of the control structure,
                                                                !! because the tracer registry will use this memory,
                                                                !! but it also means that any updates to this
                                                                !! structure in the calling module will be
                                                                !! available subsequently to the tracer registry.
  real, optional                                :: OBC_scalar   !< If present, use scalar value for segment tracer inflow concentration.
  logical, optional                             :: OBC_array    !< If true, use array values for segment tracer inflow concentration.


! Local variables
  integer :: ntseg
  integer :: isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB
  character(len=256) :: mesg    ! Message for error messages.

!  if (.not. associated(segment%tr_Reg)) call segment_tracer_registry_init(param_file, segment)
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

  if (present(tr_desc_ptr)) then
    segment%tr_Reg%Tr(ntseg)%vd => tr_desc_ptr
  else
    allocate(segment%tr_Reg%Tr(ntseg)%vd) ; segment%tr_Reg%Tr(ntseg)%vd = tr_desc
  endif

  call query_vardesc(segment%tr_Reg%Tr(ntseg)%vd, name=segment%tr_Reg%Tr(ntseg)%name)

  if (segment%tr_Reg%locked) call MOM_error(FATAL, &
      "MOM register_tracer was called for variable "//trim(segment%tr_Reg%Tr(ntseg)%name)//&
      " with a locked tracer registry.")

  if (present(OBC_scalar)) segment%tr_Reg%Tr(ntseg)%OBC_inflow_conc = OBC_scalar ! initialize tracer value later
  if (present(OBC_array)) then
    if (segment%is_E_or_W) then
      allocate(segment%tr_Reg%Tr(ntseg)%t(IsdB:IedB,jsd:jed,1:GV%ke))
    elseif (segment%is_N_or_S) then
      allocate(segment%tr_Reg%Tr(ntseg)%t(isd:ied,JsdB:JedB,1:GV%ke))
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

subroutine register_temp_salt_segments(GV, OBC, tv, vd_T, vd_S, param_file)
  type(verticalGrid_type),    intent(in)    :: GV         !< ocean vertical grid structure
  type(ocean_OBC_type),       pointer       :: OBC        !< Open boundary structure
  type(thermo_var_ptrs),      intent(in)    :: tv         !< Thermodynamics structure
  type(vardesc),              intent(in)    :: vd_T       !< Temperature descriptor
  type(vardesc),              intent(in)    :: vd_S       !< Salinity descriptor
  type(param_file_type),      intent(in)    :: param_file !< file to parse for  model parameter values

! Local variables
  integer :: isd, ied, IsdB, IedB, jsd, jed, JsdB, JedB, nz, nf
  integer :: i, j, k, n
  type(OBC_segment_type), pointer :: segment => NULL() ! pointer to segment type list

  if (.not. associated(OBC)) return

  do n=1, OBC%number_of_segments
    segment=>OBC%segment(n)
    if (.not. segment%on_pe) cycle

    if (associated(segment%tr_Reg)) &
         call MOM_error(FATAL,"register_temp_salt_segments: tracer array was previously allocated")

    call register_segment_tracer(vd_T, param_file, GV, segment, &
                                 OBC_array=segment%temp_segment_data_exists)
    call register_segment_tracer(vd_S, param_file, GV, segment, &
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
  type(OBC_segment_type), pointer :: segment ! pointer to segment type list

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
  enddo
end subroutine fill_temp_salt_segments

!> Find the region outside of all open boundary segments and
!! make sure it is set to land mask. Gonna need to know global land
!! mask as well to get it right...
subroutine mask_outside_OBCs(G, param_file, OBC)
  type(dyn_horgrid_type),  intent(inout) :: G          !< Ocean grid structure
  type(param_file_type),   intent(in)    :: param_file !< Parameter file handle
  type(ocean_OBC_type),    pointer       :: OBC        !< Open boundary structure

! Local variables
  integer :: isd, ied, IsdB, IedB, jsd, jed, JsdB, JedB, n
  integer :: i, j
  real    :: min_depth
  integer, parameter  :: cin = 3, cout = 4, cland = -1, cedge = -2
  type(OBC_segment_type), pointer :: segment ! pointer to segment type list
  real, allocatable, dimension(:,:) :: color, color2  ! For sorting inside from outside,
                                                      ! two different ways

  if (.not. associated(OBC)) return

  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 default=0.0, do_not_log=.true.)

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
    if (color(i,j) /= color2(i,j)) call MOM_error(FATAL, &
        "MOM_open_boundary: inconsistent OBC segments.")
    if (color(i,j) == cout) G%bathyT(i,j) = min_depth
  enddo ; enddo

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

subroutine open_boundary_register_restarts(HI, GV, OBC_CS,restart_CSp)
  type(hor_index_type), intent(in) :: HI
  type(verticalGrid_type), pointer, intent(in) :: GV
  type(ocean_OBC_type), pointer, intent(inout) :: OBC_CS !< OBC data structure
  type(MOM_restart_CS), pointer, intent(inout) :: restart_CSp
  type(vardesc) :: vd
  logical :: rx_normal_associated
  integer :: n

  if (.not. associated(OBC_CS)) &
       call MOM_error(FATAL, "open_boundary_register_restarts: Called with "//&
                      "uninitialized OBC control structure")

  if (associated(OBC_CS%rx_normal) .or. associated(OBC_CS%ry_normal)) &
       call MOM_error(FATAL, "open_boundary_register_restarts: Restart "//&
                      "arrays were previously allocated")


  rx_normal_associated = .false.

  do n=1,OBC_CS%number_of_segments
    if (associated(OBC_CS%segment(n)%rx_normal)) rx_normal_associated = .true.
  enddo

  if (rx_normal_associated) then
    allocate(OBC_CS%rx_normal(HI%isdB:HI%iedB,HI%jsd:HI%jed,GV%ke))
    OBC_CS%rx_normal(:,:,:) = 0.0
    vd = var_desc("rx_normal","m s-1", "Normal Phase Speed for EW OBCs",'u','L')
    call register_restart_field(OBC_CS%rx_normal, vd, .true., restart_CSp)
    allocate(OBC_CS%ry_normal(HI%isd:HI%ied,HI%jsdB:HI%jedB,GV%ke))
    OBC_CS%ry_normal(:,:,:) = 0.0
    vd = var_desc("ry_normal","m s-1", "Normal Phase Speed for NS OBCs",'v','L')
    call register_restart_field(OBC_CS%ry_normal, vd, .true., restart_CSp)
  endif

end subroutine open_boundary_register_restarts

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
