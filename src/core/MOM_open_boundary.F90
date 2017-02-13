!> Controls where open boundary conditions are applied
module MOM_open_boundary

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_domains, only : pass_var, pass_vector
use MOM_domains, only : To_All, SCALAR_PAIR, CGRID_NE
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type, log_param
use MOM_grid, only : ocean_grid_type, hor_index_type
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_io, only : EAST_FACE, NORTH_FACE
use MOM_io, only : slasher, read_data, field_size
use MOM_obsolete_params, only : obsolete_logical, obsolete_int, obsolete_real, obsolete_char
use MOM_string_functions, only : extract_word, remove_spaces
use MOM_tracer_registry, only : add_tracer_OBC_values, tracer_registry_type
use MOM_variables, only : thermo_var_ptrs
use time_interp_external_mod, only : init_external_field, time_interp_external
use MOM_remapping, only : remappingSchemesDoc, remappingDefaultScheme, remapping_CS, initialize_remapping
use MOM_remapping, only : remapping_core_h, end_remapping
use MOM_regridding, only : regridding_CS
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public open_boundary_config
public open_boundary_init
public open_boundary_query
public open_boundary_end
public open_boundary_impose_normal_slope
public open_boundary_impose_land_mask
public radiation_open_bdry_conds
public set_Flather_data
public update_obc_segment_data
public fill_OBC_halos
public open_boundary_test_extern_uv

integer, parameter, public :: OBC_NONE = 0, OBC_SIMPLE = 1, OBC_WALL = 2
integer, parameter, public :: OBC_FLATHER = 3
integer, parameter, public :: OBC_RADIATION = 4
integer, parameter, public :: OBC_DIRECTION_N = 100 !< Indicates the boundary is an effective northern boundary
integer, parameter, public :: OBC_DIRECTION_S = 200 !< Indicates the boundary is an effective southern boundary
integer, parameter, public :: OBC_DIRECTION_E = 300 !< Indicates the boundary is an effective eastern boundary
integer, parameter, public :: OBC_DIRECTION_W = 400 !< Indicates the boundary is an effective western boundary
integer, parameter         :: MAX_OBC_FIELDS = 100  !< Maximum number of data fields needed for OBC segments

type, public :: OBC_segment_data_type
  integer :: fid                                ! handle from FMS associated with segment data on disk
  integer :: fid_dz                             ! handle from FMS associated with segment thicknesses on disk
  character(len=8)                :: name       ! a name identifier for the segment data
  real, pointer, dimension(:,:,:) :: buffer_src=>NULL() ! buffer for segment data located at cell faces
                                                ! and on the original vertical grid
  integer                         :: nk_src     ! Number of vertical levels in the source data
  real, dimension(:,:,:), pointer :: dz_src=>NULL()     ! vertical grid cell spacing of the incoming segment data (m)
  real, dimension(:,:,:), pointer :: buffer_dst=>NULL() ! buffer src data remapped to the target vertical grid
  real, dimension(:,:), pointer   :: bt_vel=>NULL() ! barotropic velocity (m s-1)
  real                            :: value              ! constant value if fid is equal to -1
end type OBC_segment_data_type

type, public :: OBC_segment_type
  logical :: Flather        !< If true, applies Flather + Chapman radiation of barotropic gravity waves.
  logical :: radiation      !< If true, 1D Orlanksi radiation boundary conditions are applied.
                            !! If False, a gradient condition is applied.
  logical :: oblique        !< Oblique waves supported at radiation boundary.
  logical :: nudged         !< Optional supplement to radiation boundary.
  logical :: specified      !< Boundary fixed to external value.
  logical :: gradient       !< Zero gradient at boundary.
  logical :: values_needed  !< Whether or not external OBC fields are needed.
  logical :: legacy         !< Old code for tangential BT velocities.
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
  real, pointer, dimension(:,:)   :: Cg=>NULL()     !< The external gravity
                                                    !< wave speed (m -s) at OBC-points.
  real, pointer, dimension(:,:)   :: Htot=>NULL()   !< The total column thickness (m) at OBC-points.
  real, pointer, dimension(:,:,:) :: h=>NULL()      !< The cell thickness (m) at OBC-points.
  real, pointer, dimension(:,:,:) :: e=>NULL()      !< The interface height (m?) at OBC-points.
  real, pointer, dimension(:,:,:) :: tangent_vel=>NULL()    !< The layer velocity tangential to the
                                                            !! OB segment (m s-1).
  real, pointer, dimension(:,:)   :: tangent_vel_bt=>NULL() !< The barotropic velocity tangential
                                                            !! to the OB segment (m s-1).
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
  real, pointer, dimension(:,:,:) :: rx_normal => NULL()    !< The rx_old_u value for radiation coeff
                                                            !! for normal velocity
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
  real :: g_Earth
  ! Properties of the segments used.
  type(OBC_segment_type), pointer, dimension(:) :: &
    segment => NULL()   !< List of segment objects.
  ! Which segment object describes the current point.
  integer, pointer, dimension(:,:) :: &
    OBC_segment_u => NULL(), &   !< Segment number of u-points.
    OBC_segment_v => NULL()      !< Segment number of v-points.

  !   The following can be used to specify the outer-domain values of the
  ! surface height and barotropic velocity.  If these are not allocated, the
  ! default with Flather boundary conditions is the same as if they were
  ! filled with zeros.  With simple OBCs, these should not be allocated.
  real, pointer, dimension(:,:) :: &
    ubt_outer => NULL(), &    !< The u-velocity in the outer domain, in m s-1.
    vbt_outer => NULL(), &    !< The v-velocity in the outer domain, in m s-1.
    eta_outer_u => NULL(), &  !< The SSH anomaly in the outer domain, in m or kg m-2.
    eta_outer_v => NULL()     !< The SSH anomaly in the outer domain, in m or kg m-2.

  ! The following apply at points with OBC_kind_[uv] = OBC_SIMPLE and/or
  ! if nudging is turned on.
  real, pointer, dimension(:,:,:) :: &
    u => NULL(), &  !< The prescribed values of the zonal velocity (u) at OBC points.
    v => NULL(), &  !< The prescribed values of the meridional velocity (v) at OBC points.
    uh => NULL(), & !< The prescribed values of the zonal volume transport (uh) at OBC points.
    vh => NULL()    !< The prescribed values of the meridional volume transport (vh) at OBC points.

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
  character(len=200) :: OBC_user_config
  type(remapping_CS), pointer         :: remap_CS   ! ALE remapping control structure for segments only
end type ocean_OBC_type

integer :: id_clock_pass

character(len=40)  :: mod = "MOM_open_boundary" ! This module's name.
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
  type(dyn_horgrid_type),  intent(in)    :: G !< Ocean grid structure
  type(param_file_type),   intent(in)    :: param_file !< Parameter file handle
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundary control structure
  ! Local variables
  integer :: l ! For looping over segments
  character(len=15) :: segment_param_str ! The run-time parameter name for each segment
  character(len=100) :: segment_str      ! The contents (rhs) for parameter "segment_param_str"
  character(len=200) :: config1          ! String for OBC_USER_CONFIG

  allocate(OBC)

  call log_version(param_file, mod, version, "Controls where open boundaries are located, what "//&
                 "kind of boundary condition to impose, and what data to apply, if any.")
  call get_param(param_file, mod, "OBC_NUMBER_OF_SEGMENTS", OBC%number_of_segments, &
                 "The number of open boundary segments.", &
                 default=0)
  call get_param(param_file, mod, "G_EARTH", OBC%g_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80)
  call get_param(param_file, mod, "OBC_USER_CONFIG", config1, &
                 "A string that sets how the open boundary conditions are \n"//&
                 " configured: \n", default="none", do_not_log=.true.)
  call get_param(param_file, mod, "NK", OBC%ke, &
                 "The number of model layers", default=0, do_not_log=.true.)

  if (config1 .ne. "none") OBC%user_BCs_set_globally = .true.

  if (OBC%number_of_segments > 0) then
    call get_param(param_file, mod, "OBC_ZERO_VORTICITY", OBC%zero_vorticity, &
                   "If true, sets relative vorticity to zero on open boundaries.", &
                   default=.false.)
    call get_param(param_file, mod, "OBC_FREESLIP_VORTICITY", OBC%freeslip_vorticity, &
                   "If true, sets the normal gradient of tangential velocity to\n"// &
                   "zero in the relative vorticity on open boundaries. This cannot\n"// &
                   "be true if OBC_ZERO_VORTICITY is True.", default=.false.)
    if (OBC%zero_vorticity .and. OBC%freeslip_vorticity) call MOM_error(FATAL, &
                   "MOM_open_boundary.F90, open_boundary_config: "//&
                   "Only one of OBC_ZERO_VORTICITY and OBC_FREESLIP_VORTICITY can be True at once.")
    call get_param(param_file, mod, "OBC_ZERO_STRAIN", OBC%zero_strain, &
                   "If true, sets the strain used in the stress tensor to zero on open boundaries.", &
                   default=.false.)
    call get_param(param_file, mod, "OBC_FREESLIP_STRAIN", OBC%freeslip_strain, &
                   "If true, sets the normal gradient of tangential velocity to\n"// &
                   "zero in the strain use in the stress tensor on open boundaries. This cannot\n"// &
                   "be true if OBC_ZERO_STRAIN is True.", default=.false.)
    if (OBC%zero_strain .and. OBC%freeslip_strain) call MOM_error(FATAL, &
                   "MOM_open_boundary.F90, open_boundary_config: "//&
                   "Only one of OBC_ZERO_STRAIN and OBC_FREESLIP_STRAIN can be True at once.")
    call get_param(param_file, mod, "OBC_ZERO_BIHARMONIC", OBC%zero_biharmonic, &
                   "If true, zeros the Laplacian of flow on open boundaries in the biharmonic\n"//&
                   "viscosity term.", default=.false.)
    ! Allocate everything
    ! Note the 0-segment is needed when %OBC_segment_u/v(:,:) = 0
    allocate(OBC%segment(0:OBC%number_of_segments))
    do l=0,OBC%number_of_segments
      OBC%segment(l)%Flather = .false.
      OBC%segment(l)%radiation = .false.
      OBC%segment(l)%oblique = .false.
      OBC%segment(l)%nudged = .false.
      OBC%segment(l)%specified = .false.
      OBC%segment(l)%gradient = .false.
      OBC%segment(l)%values_needed = .false.
      OBC%segment(l)%legacy = .false.
      OBC%segment(l)%direction = OBC_NONE
      OBC%segment(l)%is_N_or_S = .false.
      OBC%segment(l)%is_E_or_W = .false.
      OBC%segment(l)%Tnudge_in = 0.0
      OBC%segment(l)%Tnudge_out = 0.0
      OBC%segment(l)%num_fields = 0.0
    enddo
    allocate(OBC%OBC_segment_u(G%IsdB:G%IedB,G%jsd:G%jed)) ; OBC%OBC_segment_u(:,:) = OBC_NONE
    allocate(OBC%OBC_segment_v(G%isd:G%ied,G%JsdB:G%JedB)) ; OBC%OBC_segment_v(:,:) = OBC_NONE

    do l = 1, OBC%number_of_segments
      write(segment_param_str(1:15),"('OBC_SEGMENT_',i3.3)") l
      call get_param(param_file, mod, segment_param_str, segment_str, &
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

    if (open_boundary_query(OBC, needs_ext_seg_data=.true.)) &
      call initialize_segment_data(G, OBC, param_file)
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
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  ! There is a problem with the order of the OBC initialization
  ! with respect to ALE_init. Currently handling this by copying the
  ! param file so that I can use it later in step_MOM in order to finish
  ! initializing segments on the first step.

  call get_param(PF, mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)

  call get_param(PF, mod, "REMAPPING_SCHEME", remappingScheme, &
          "This sets the reconstruction scheme used\n"//&
          "for vertical remapping for all variables.\n"//&
          "It can be one of the following schemes:\n"//&
          trim(remappingSchemesDoc), default=remappingDefaultScheme)
  call get_param(PF, mod, "FATAL_CHECK_RECONSTRUCTIONS", check_reconstruction, &
          "If true, cell-by-cell reconstructions are checked for\n"//&
          "consistency and if non-monotonicity or an inconsistency is\n"//&
          "detected then a FATAL error is issued.", default=.false.)
  call get_param(PF, mod, "FATAL_CHECK_REMAPPING", check_remapping, &
          "If true, the results of remapping are checked for\n"//&
          "conservation and new extrema and if an inconsistency is\n"//&
          "detected then a FATAL error is issued.", default=.false.)
  call get_param(PF, mod, "REMAP_BOUND_INTERMEDIATE_VALUES", force_bounds_in_subcell, &
          "If true, the values on the intermediate grid used for remapping\n"//&
          "are forced to be bounded, which might not be the case due to\n"//&
          "round off.", default=.false.)

  allocate(OBC%remap_CS)
  call initialize_remapping(OBC%remap_CS, remappingScheme, boundary_extrapolation = .false., &
       check_reconstruction=check_reconstruction, &
       check_remapping=check_remapping, force_bounds_in_subcell=force_bounds_in_subcell)

  if (OBC%user_BCs_set_globally) return

  do n=1, OBC%number_of_segments
    segment => OBC%segment(n)

    write(segnam,"('OBC_SEGMENT_',i3.3,'_DATA')") n
    write(suffix,"('_segment_',i3.3)") n
    call get_param(PF, mod, segnam, segstr)

    call parse_segment_data_str(trim(segstr), fields=fields, num_fields=num_fields)
    if (num_fields == 0) cycle ! cycle to next segment

    if (segment%on_pe) then
      allocate(segment%field(num_fields))

      if (segment%Flather) then
        if (num_fields /= 3) call MOM_error(FATAL, &
                   "MOM_open_boundary, initialize_segment_data: "//&
                   "Need three inputs for Flather")
        segment%num_fields = 3 ! these are the input fields required for the Flather option
                                       ! note that this is assuming that the inputs are coming in this order
                                       ! and independent of the input param string . Needs cleanup - mjh
        allocate(segment%field_names(segment%num_fields))
        segment%field_names(:)='None'
        segment%field_names(1)='UO'
        segment%field_names(2)='VO'
        segment%field_names(3)='ZOS'
      endif
    endif

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
        if (segment%on_pe) then
         segment%values_needed = .true. ! Indicates that i/o will be needed for this segment
         segment%field(m)%name = trim(fields(m))
         filename = trim(inputdir)//trim(filename)
         fieldname = trim(fieldname)//trim(suffix)
         call field_size(filename,fieldname,siz,no_domain=.true.)
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
         segment%field(m)%fid = init_external_field(trim(filename),trim(fieldname))
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
            segment%field(m)%fid_dz = init_external_field(trim(filename),trim(fieldname))
          else
            segment%field(m)%nk_src=1
          endif
       else
          segment%field(m)%fid = -1
          segment%field(m)%value = value
       endif
      endif
    enddo
  enddo

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
  ! The values are nonsence if the segment is not on this PE.
  seg%HI%IsdB = min( max(Isg, G%HI%IsdB), G%HI%IedB)
  seg%HI%IedB = min( max(Ieg, G%HI%IsdB), G%HI%IedB)
  seg%HI%isd = min( max(Isg+1, G%HI%isd), G%HI%ied)
  seg%HI%ied = min( max(Ieg, G%HI%isd), G%HI%ied)
  seg%HI%IscB = min( max(Isg, G%HI%IscB), G%HI%IecB)
  seg%HI%IecB = min( max(Ieg, G%HI%IscB), G%HI%IecB)
  seg%HI%isc = min( max(Isg+1, G%HI%isc), G%HI%iec)
  seg%HI%iec = min( max(Ieg, G%HI%isc), G%HI%iec)

  ! This is the j-extent of the segment on this PE.
  ! The values are nonsence if the segment is not on this PE.
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
  integer :: j, this_kind, a_loop
  character(len=32) :: action_str(5)

  ! This returns the global indices for the segment
  call parse_segment_str(G%ieg, G%jeg, segment_str, I_obc, Js_obc, Je_obc, action_str )

  call setup_segment_indices(G, OBC%segment(l_seg),I_obc,I_obc,Js_obc,Je_obc)

  I_obc = I_obc - G%idg_offset ! Convert to local tile indices on this tile
  Js_obc = Js_obc - G%jdg_offset ! Convert to local tile indices on this tile
  Je_obc = Je_obc - G%jdg_offset ! Convert to local tile indices on this tile

  this_kind = OBC_NONE

  ! Hack to extend segment by one point
  if (Js_obc<Je_obc) then
    Js_obc = Js_obc - 1 ; Je_obc = Je_obc + 1
  else
    Js_obc = Js_obc + 1 ; Je_obc = Je_obc - 1
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
      this_kind = OBC_FLATHER
      OBC%segment(l_seg)%Flather = .true.
      OBC%Flather_u_BCs_exist_globally = .true.
      OBC%open_u_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'ORLANSKI') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%Flather_u_BCs_exist_globally = .true.
      OBC%open_u_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'OBLIQUE') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%segment(l_seg)%oblique = .true.
      OBC%Flather_u_BCs_exist_globally = .true.
      OBC%open_u_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED') then
      OBC%segment(l_seg)%nudged = .true.
      OBC%segment(l_seg)%Tnudge_in = 1.0/(3*86400)
      OBC%segment(l_seg)%Tnudge_out = 1.0/(360*86400)
      OBC%nudged_u_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'GRADIENT') then
      OBC%segment(l_seg)%gradient = .true.
      OBC%open_u_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'LEGACY') then
      this_kind = OBC_FLATHER
      OBC%segment(l_seg)%legacy = .true.
      OBC%segment(l_seg)%Flather = .true.
      OBC%segment(l_seg)%radiation = .true.
      OBC%Flather_u_BCs_exist_globally = .true.
      OBC%open_u_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE') then
      OBC%segment(l_seg)%specified = .true.
      OBC%specified_u_BCs_exist_globally = .true. ! This avoids deallocation
      ! Hack to undo the hack above for SIMPLE BCs
      Js_obc = Js_obc + 1 ; Je_obc = Je_obc - 1
    else
      call MOM_error(FATAL, "MOM_open_boundary.F90, setup_u_point_obc: "//&
                     "String '"//trim(action_str(a_loop))//"' not understood.")
    endif

    if (I_obc<G%HI%IsdB .or. I_obc>G%HI%IedB) return ! Boundary is not on tile
    if (Js_obc<G%HI%JsdB .and. Je_obc<G%HI%JsdB) return ! Segment is not on tile
    if (Js_obc>G%HI%JedB) return ! Segment is not on tile
  enddo ! a_loop

  OBC%segment(l_seg)%on_pe = .true.
  OBC%segment(l_seg)%is_E_or_W = .true.

  do j=G%HI%jsd, G%HI%jed
    if (j>Js_obc .and. j<=Je_obc) then
      OBC%OBC_segment_u(I_obc,j) = l_seg
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
  integer :: i, this_kind, a_loop
  character(len=32) :: action_str(5)

  ! This returns the global indices for the segment
  call parse_segment_str(G%ieg, G%jeg, segment_str, J_obc, Is_obc, Ie_obc, action_str )

  call setup_segment_indices(G, OBC%segment(l_seg),Is_obc,Ie_obc,J_obc,J_obc)

  J_obc = J_obc - G%jdg_offset ! Convert to local tile indices on this tile
  Is_obc = Is_obc - G%idg_offset ! Convert to local tile indices on this tile
  Ie_obc = Ie_obc - G%idg_offset ! Convert to local tile indices on this tile
  this_kind = OBC_NONE

  ! Hack to extend segment by one point
  if (Is_obc<Ie_obc) then
    Is_obc = Is_obc - 1 ; Ie_obc = Ie_obc + 1
  else
    Is_obc = Is_obc + 1 ; Ie_obc = Ie_obc - 1
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
      this_kind = OBC_FLATHER
      OBC%segment(l_seg)%Flather = .true.
      OBC%Flather_v_BCs_exist_globally = .true.
      OBC%open_v_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'ORLANSKI') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%Flather_v_BCs_exist_globally = .true.
      OBC%open_v_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'OBLIQUE') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%segment(l_seg)%oblique = .true.
      OBC%Flather_v_BCs_exist_globally = .true.
      OBC%open_v_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED') then
      OBC%segment(l_seg)%nudged = .true.
      OBC%segment(l_seg)%Tnudge_in = 1.0/(3*86400)
      OBC%segment(l_seg)%Tnudge_out = 1.0/(360*86400)
      OBC%nudged_v_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'GRADIENT') then
      OBC%segment(l_seg)%gradient = .true.
      OBC%open_v_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'LEGACY') then
      this_kind = OBC_FLATHER
      OBC%segment(l_seg)%legacy = .true.
      OBC%segment(l_seg)%radiation = .true.
      OBC%segment(l_seg)%Flather = .true.
      OBC%Flather_v_BCs_exist_globally = .true.
      OBC%open_v_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE') then
      OBC%segment(l_seg)%specified = .true.
      OBC%specified_v_BCs_exist_globally = .true. ! This avoids deallocation
      ! Hack to undo the hack above for SIMPLE BCs
      Is_obc = Is_obc + 1 ; Ie_obc = Ie_obc - 1
    else
      call MOM_error(FATAL, "MOM_open_boundary.F90, setup_v_point_obc: "//&
                     "String '"//trim(action_str(a_loop))//"' not understood.")
    endif

    if (J_obc<G%HI%JsdB .or. J_obc>G%HI%JedB) return ! Boundary is not on tile
    if (Is_obc<G%HI%IsdB .and. Ie_obc<G%HI%IsdB) return ! Segment is not on tile
    if (Is_obc>G%HI%IedB) return ! Segment is not on tile
  enddo ! a_loop

  OBC%segment(l_seg)%on_pe = .true.
  OBC%segment(l_seg)%is_N_or_S = .true.

  do i=G%HI%isd, G%HI%ied
    if (i>Is_obc .and. i<=Ie_obc) then
      OBC%OBC_segment_v(i,J_obc) = l_seg
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

  if ( OBC%Flather_u_BCs_exist_globally .or. OBC%Flather_v_BCs_exist_globally ) then
    call get_param(param_file, mod, "OBC_RADIATION_MAX", OBC%rx_max, &
                   "The maximum magnitude of the baroclinic radiation \n"//&
                   "velocity (or speed of characteristics).  This is only \n"//&
                   "used if one of the open boundary segments is using Orlanski.", &
                   units="m s-1", default=10.0)
    call get_param(param_file, mod, "OBC_RAD_VEL_WT", OBC%gamma_uv, &
                   "The relative weighting for the baroclinic radiation \n"//&
                   "velocities (or speed of characteristics) at the new \n"//&
                   "time level (1) or the running mean (0) for velocities. \n"//&
                   "Valid values range from 0 to 1. This is only used if \n"//&
                   "one of the open boundary segments is using Orlanski.", &
                   units="nondim",  default=0.3)
    call get_param(param_file, mod, "OBC_RAD_THICK_WT", OBC%gamma_h, &
                   "The relative weighting for the baroclinic radiation \n"//&
                   "velocities (or speed of characteristics) at the new \n"//&
                   "time level (1) or the running mean (0) for thicknesses. \n"//&
                   "Valid values range from 0 to 1. This is only used if \n"//&
                   "one of the open boundary segments is using Orlanski.", &
                   units="nondim",  default=0.2)
  endif

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
  if (.not. associated(OBC)) return
  if (associated(OBC%segment)) deallocate(OBC%segment)
  if (associated(OBC%OBC_segment_u)) deallocate(OBC%OBC_segment_u)
  if (associated(OBC%OBC_segment_v)) deallocate(OBC%OBC_segment_v)
  if (associated(OBC%ubt_outer)) deallocate(OBC%ubt_outer)
  if (associated(OBC%vbt_outer)) deallocate(OBC%vbt_outer)
  if (associated(OBC%eta_outer_u)) deallocate(OBC%eta_outer_u)
  if (associated(OBC%eta_outer_v)) deallocate(OBC%eta_outer_v)
  if (associated(OBC%u)) deallocate(OBC%u)
  if (associated(OBC%v)) deallocate(OBC%v)
  if (associated(OBC%uh)) deallocate(OBC%uh)
  if (associated(OBC%vh)) deallocate(OBC%vh)
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
  integer :: i, j
  logical :: bc_north, bc_south, bc_east, bc_west

  if (.not.associated(OBC)) return

  do J=G%jsd+1,G%jed-1 ; do i=G%isd+1,G%ied-1
    bc_north = .false. ; bc_south = .false. ; bc_east = .false. ; bc_west = .false.
    if (associated(OBC%OBC_segment_u)) then
      if (OBC%segment(OBC%OBC_segment_u(I,j))%direction == OBC_DIRECTION_E &
          .and. .not. OBC%segment(OBC%OBC_segment_u(I,j))%specified) bc_east = .true.
      if (OBC%segment(OBC%OBC_segment_u(I-1,j))%direction == OBC_DIRECTION_W &
          .and. .not. OBC%segment(OBC%OBC_segment_u(I-1,j))%specified) bc_west = .true.
    endif
    if (associated(OBC%OBC_segment_v)) then
      if (OBC%segment(OBC%OBC_segment_v(i,J))%direction == OBC_DIRECTION_N &
          .and. .not. OBC%segment(OBC%OBC_segment_v(i,J))%specified) bc_north = .true.
      if (OBC%segment(OBC%OBC_segment_v(i,J-1))%direction == OBC_DIRECTION_S &
          .and. .not. OBC%segment(OBC%OBC_segment_v(i,J-1))%specified) bc_south = .true.
    endif
    if (bc_north) depth(i,j+1) = depth(i,j)
    if (bc_south) depth(i,j-1) = depth(i,j)
    if (bc_east) depth(i+1,j) = depth(i,j)
    if (bc_west) depth(i-1,j) = depth(i,j)
    ! Convex corner cases
    if (bc_north.and.bc_east) depth(i+1,j+1) = depth(i,j)
    if (bc_north.and.bc_west) depth(i-1,j+1) = depth(i,j)
    if (bc_south.and.bc_east) depth(i+1,j-1) = depth(i,j)
    if (bc_south.and.bc_west) depth(i-1,j-1) = depth(i,j)
  enddo ; enddo

end subroutine open_boundary_impose_normal_slope

!> Reconcile masks and open boundaries, deallocate OBC on PEs where it is not needed.
!! Also adjust u- and v-point cell area on specified open boundaries.
subroutine open_boundary_impose_land_mask(OBC, G, areaCu, areaCv)
  type(ocean_OBC_type),              pointer       :: OBC !< Open boundary control structure
  type(dyn_horgrid_type),            intent(in)    :: G !< Ocean grid structure
  real, dimension(SZIB_(G),SZJ_(G)), intent(inout) :: areaCu !< Area of a u-cell (m2)
  real, dimension(SZI_(G),SZJB_(G)), intent(inout) :: areaCv !< Area of a u-cell (m2)
  ! Local variables
  integer :: i, j
  logical :: any_U, any_V

  if (.not.associated(OBC)) return

  ! Sweep along u-segments and delete the OBC for blocked points.
  if (associated(OBC%OBC_segment_u)) then
    do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB
      if (G%mask2dCu(I,j) == 0 .and. (OBC%OBC_segment_u(I,j) /= OBC_NONE)) then
        if (.not. OBC%segment(OBC%OBC_segment_u(I,j))%specified) then
          OBC%OBC_segment_u(I,j) = OBC_NONE
        endif
      endif
    enddo ; enddo
  endif

  ! Sweep along v-segments and delete the OBC for blocked points.
  if (associated(OBC%OBC_segment_v)) then
    do J=G%JsdB,G%JedB ; do i=G%isd,G%ied
      if (G%mask2dCv(i,J) == 0 .and. (OBC%OBC_segment_v(i,J) /= OBC_NONE)) then
        if (.not. OBC%segment(OBC%OBC_segment_v(i,J))%specified) then
          OBC%OBC_segment_v(I,j) = OBC_NONE
        endif
      endif
    enddo ; enddo
  endif

  ! Sweep along u-segments and for %specified BC points reset the u-point area which was masked out
  if (associated(OBC%OBC_segment_u)) then
    do j=G%jsd,G%jed ; do I=G%isd,G%ied-1
      if (OBC%OBC_segment_u(I,j) /= OBC_NONE) then
        if (OBC%segment(OBC%OBC_segment_u(I,j))%specified) then
          if (OBC%segment(OBC%OBC_segment_u(I,j))%direction == OBC_DIRECTION_W) then
            areaCu(I,j) = G%areaT(i+1,j)
           !G%IareaCu(I,j) = G%IareaT(i+1,j) ?
          elseif (OBC%segment(OBC%OBC_segment_u(I,j))%direction == OBC_DIRECTION_E) then
            areaCu(I,j) = G%areaT(i,j)
           !G%IareaCu(I,j) = G%IareaT(i,j) ?
          endif
        endif
      endif
    enddo ; enddo
  endif

  ! Sweep along v-segments and for %specified BC points reset the v-point area which was masked out
  if (associated(OBC%OBC_segment_v)) then
    do J=G%jsd,G%jed-1 ; do i=G%isd,G%ied
      if (OBC%OBC_segment_v(i,J) /= OBC_NONE) then
        if (OBC%segment(OBC%OBC_segment_v(i,J))%specified) then
          if (OBC%segment(OBC%OBC_segment_v(i,J))%direction == OBC_DIRECTION_S) then
            areaCv(i,J) = G%areaT(i,j+1)
           !G%IareaCv(i,J) = G%IareaT(i,j+1) ?
          elseif (OBC%segment(OBC%OBC_segment_v(i,J))%direction == OBC_DIRECTION_N) then
            areaCu(i,J) = G%areaT(i,j)
           !G%IareaCu(i,J) = G%IareaT(i,j) ?
          endif
        endif
      endif
    enddo ; enddo
  endif

  any_U = .false.
  if (associated(OBC%OBC_segment_u)) then
    do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB
      ! G%mask2du will be open wherever bathymetry allows it.
      ! Bathymetry outside of the open boundary was adjusted to match
      ! the bathymetry inside so these points will be open unless the
      ! bathymetry inside the boundary was do shallow and flagged as land.
      if (OBC%OBC_segment_u(I,j) /= OBC_NONE) any_U = .true.
    enddo ; enddo
  endif

  any_V = .false.
  if (associated(OBC%OBC_segment_v)) then
    do J=G%JsdB,G%JedB ; do i=G%isd,G%ied
      if (OBC%OBC_segment_v(i,J) /= OBC_NONE) any_V = .true.
    enddo ; enddo
  endif

  OBC%OBC_pe = .true.
  if (.not.(any_U .or. any_V)) OBC%OBC_pe = .false.

end subroutine open_boundary_impose_land_mask

!> Apply radiation conditions to 3D  u,v (,h) at open boundaries
subroutine radiation_open_bdry_conds(OBC, u_new, u_old, v_new, v_old, &
                                     h_new, h_old, G, dt)
  type(ocean_grid_type),                     intent(inout) :: G !< Ocean grid structure
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary control structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: u_new !< New u values on open boundaries
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u_old !< Original unadjusted u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: v_new !< New v values on open boundaries
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v_old !< Original unadjusted v
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: h_new !< New h values on open boundaries
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h_old !< Original h values
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

  gamma_u = OBC%gamma_uv ; gamma_v = OBC%gamma_uv ; gamma_h = OBC%gamma_h
  rx_max = OBC%rx_max ; ry_max = OBC%rx_max
  do n=1,OBC%number_of_segments
     segment=>OBC%segment(n)
     if (.not. segment%on_pe) cycle
     if (segment%radiation) call gradient_at_q_points(G,segment,u_old,v_old)
     if (segment%direction == OBC_DIRECTION_E) then
       I=segment%HI%IscB
       do k=1,nz ;  do j=segment%HI%jsc,segment%HI%jec
         if (segment%legacy) then
           dhdt = u_old(I-1,j,k)-u_new(I-1,j,k) !old-new
           dhdx = u_new(I-1,j,k)-u_new(I-2,j,k) !in new time backward sasha for I-1
           rx_new = 0.0
           if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max) ! outward phase speed
           rx_avg = (1.0-gamma_u)*segment%rx_normal(I,j,k) + gamma_u*rx_new
           segment%rx_normal(I,j,k) = rx_avg
           u_new(I,j,k) = (u_old(I,j,k) + rx_avg*u_new(I-1,j,k)) / (1.0+rx_avg)
         elseif (segment%radiation) then
           dhdt = u_old(I-1,j,k)-u_new(I-1,j,k) !old-new
           dhdx = u_new(I-1,j,k)-u_new(I-2,j,k) !in new time backward sasha for I-1
           if (segment%oblique) then
             if (dhdt*(segment%grad_normal(J,1,k) + segment%grad_normal(J-1,1,k)) > 0.0) then
               dhdy = segment%grad_normal(J-1,1,k)
             elseif (dhdt*(segment%grad_normal(J,1,k) + segment%grad_normal(J-1,1,k)) == 0.0) then
               dhdy = 0.0
             else
               dhdy = segment%grad_normal(J,1,k)
             endif
           endif
           if (dhdt*dhdx < 0.0) dhdt = 0.0
           if (dhdx == 0.0) dhdx=eps  ! avoid segv
           Cx = min(dhdt/dhdx,rx_max) ! default to normal radiation
           Cy = 0.0
           cff = max(dhdx*dhdx,eps)
           if (segment%oblique) then
             cff = max(dhdx*dhdx + dhdy*dhdy, eps)
             if (dhdy==0.) dhdy=eps ! avoid segv
             Cy = min(cff,max(dhdt/dhdy,-cff))
           endif
           u_new(I,j,k) = ((cff*u_old(I,j,k) + Cx*u_new(I-1,j,k)) - &
              (max(Cy,0.0)*segment%grad_normal(J-1,2,k) + min(Cy,0.0)*segment%grad_normal(J,2,k))) / (cff + Cx)
         endif
         if ((segment%radiation .or. segment%legacy) .and. segment%nudged) then
           if (dhdt*dhdx < 0.0) then
             tau = segment%Tnudge_in
           else
             tau = segment%Tnudge_out
           endif
           u_new(I,j,k) = u_new(I,j,k) + dt*tau*(segment%normal_vel(I,j,k) - u_old(I,j,k))
!          u_new(I,j,k) = u_new(I,j,k) + dt*tau*(OBC%u(I,j,k) - u_old(I,j,k))
         endif
       enddo; enddo
     endif

     if (segment%direction == OBC_DIRECTION_W) then
       I=segment%HI%IscB
       do k=1,nz ;  do j=segment%HI%jsc,segment%HI%jec
         if (segment%legacy) then
           dhdt = u_old(I+1,j,k)-u_new(I+1,j,k) !old-new
           dhdx = u_new(I+1,j,k)-u_new(I+2,j,k) !in new time forward sasha for I+1
           rx_new = 0.0
           if (dhdt*dhdx > 0.0) rx_new = min( (dhdt/dhdx), rx_max)
           rx_avg = (1.0-gamma_u)*segment%rx_normal(I,j,k) + gamma_u*rx_new
           segment%rx_normal(I,j,k) = rx_avg
           u_new(I,j,k) = (u_old(I,j,k) + rx_avg*u_new(I+1,j,k)) / (1.0+rx_avg)
         elseif (segment%radiation) then
           dhdt = u_old(I+1,j,k)-u_new(I+1,j,k) !old-new
           dhdx = u_new(I+1,j,k)-u_new(I+2,j,k) !in new time forward sasha for I+1
           if (segment%oblique) then
             if (dhdt*(segment%grad_normal(J,1,k) + segment%grad_normal(J-1,1,k)) > 0.0) then
               dhdy = segment%grad_normal(J-1,1,k)
             elseif (dhdt*(segment%grad_normal(J,1,k) + segment%grad_normal(J-1,1,k)) == 0.0) then
               dhdy = 0.0
             else
               dhdy = segment%grad_normal(J,1,k)
             endif
           endif
           if (dhdt*dhdx < 0.0) dhdt = 0.0
           if (dhdx == 0.0) dhdx=eps  ! avoid segv
           Cx = min(dhdt/dhdx,rx_max) ! default to normal flow only
           Cy = 0.
           cff = max(dhdx*dhdx, eps)
           if (segment%oblique) then
             cff = max(dhdx*dhdx + dhdy*dhdy, eps)
             if (dhdy==0.) dhdy=eps ! avoid segv
             Cy = min(cff,max(dhdt/dhdy,-cff))
           endif
           u_new(I,j,k) = ((cff*u_old(I,j,k) + Cx*u_new(I+1,j,k)) - &
             (max(Cy,0.0)*segment%grad_normal(J-1,2,k) + min(Cy,0.0)*segment%grad_normal(J,2,k))) / (cff + Cx)
         endif
         if ((segment%radiation .or. segment%legacy) .and. segment%nudged) then
           if (dhdt*dhdx < 0.0) then
             tau = segment%Tnudge_in
           else
             tau = segment%Tnudge_out
           endif
           u_new(I,j,k) = u_new(I,j,k) + dt*tau*(segment%normal_vel(I,j,k) - u_old(I,j,k))
!          u_new(I,j,k) = u_new(I,j,k) + dt*tau*(OBC%u(I,j,k) - u_old(I,j,k))
         endif
       enddo; enddo
     endif

     if (segment%direction == OBC_DIRECTION_N) then
       J=segment%HI%JscB
       do k=1,nz ;  do i=segment%HI%isc,segment%HI%iec
         if (segment%legacy) then
           dhdt = v_old(i,J-1,k)-v_new(i,J-1,k) !old-new
           dhdy = v_new(i,J-1,k)-v_new(i,J-2,k) !in new time backward sasha for J-1
           ry_new = 0.0
           if (dhdt*dhdy > 0.0) ry_new = min( (dhdt/dhdy), ry_max)
           ry_avg = (1.0-gamma_v)*segment%rx_normal(I,j,k) + gamma_v*ry_new
           segment%rx_normal(i,J,k) = ry_avg
           v_new(i,J,k) = (v_old(i,J,k) + ry_avg*v_new(i,J-1,k)) / (1.0+ry_avg)
         elseif (segment%radiation) then
           dhdt = v_old(i,J-1,k)-v_new(i,J-1,k) !old-new
           dhdy = v_new(i,J-1,k)-v_new(i,J-2,k) !in new time backward sasha for J-1
           if (segment%oblique) then
             if (dhdt*(segment%grad_normal(I,1,k) + segment%grad_normal(I-1,1,k)) > 0.0) then
               dhdx = segment%grad_normal(I-1,1,k)
             elseif (dhdt*(segment%grad_normal(I,1,k) + segment%grad_normal(I-1,1,k)) == 0.0) then
               dhdx = 0.0
             else
               dhdx = segment%grad_normal(I,1,k)
             endif
           endif
           if (dhdt*dhdy < 0.0) dhdt = 0.0
           if (dhdy == 0.0) dhdy=eps  ! avoid segv
           Cy = min(dhdt/dhdy,rx_max) ! default to normal flow only
           Cx = 0
           cff = max(dhdy*dhdy, eps)
           if (segment%oblique) then
             cff = max(dhdx*dhdx + dhdy*dhdy, eps)
             if (dhdx==0.) dhdx=eps ! avoid segv
             Cx = min(cff,max(dhdt/dhdx,-cff))
           endif
           v_new(i,J,k) = ((cff*v_old(i,J,k) + Cy*v_new(i,J-1,k)) - &
              (max(Cx,0.0)*segment%grad_normal(I-1,2,k) + min(Cx,0.0)*segment%grad_normal(I,2,k))) / (cff + Cy)
         endif
         if ((segment%radiation .or. segment%legacy) .and. segment%nudged) then
           if (dhdt*dhdy < 0.0) then
             tau = segment%Tnudge_in
           else
             tau = segment%Tnudge_out
           endif
           v_new(i,J,k) = v_new(i,J,k) + dt*tau*(segment%normal_vel(i,J,k) - v_old(i,J,k))
!          v_new(i,J,k) = v_new(i,J,k) + dt*tau*(OBC%v(i,J,k) - v_old(i,J,k))
         endif
       enddo; enddo
     endif


     if (segment%direction == OBC_DIRECTION_S) then
       J=segment%HI%JscB
       do k=1,nz ;  do i=segment%HI%isc,segment%HI%iec
         if (segment%legacy) then
           dhdt = v_old(i,J+1,k)-v_new(i,J+1,k) !old-new
           dhdy = v_new(i,J+1,k)-v_new(i,J+2,k) !in new time backward sasha for J-1
           ry_new = 0.0
           if (dhdt*dhdy > 0.0) ry_new = min( (dhdt/dhdy), ry_max)
           ry_avg = (1.0-gamma_v)*segment%rx_normal(I,j,k) + gamma_v*ry_new
           segment%rx_normal(i,J,k) = ry_avg
           v_new(i,J,k) = (v_old(i,J,k) + ry_avg*v_new(i,J+1,k)) / (1.0+ry_avg)
         elseif (segment%radiation) then
           dhdt = v_old(i,J+1,k)-v_new(i,J+1,k) !old-new
           dhdy = v_new(i,J+1,k)-v_new(i,J+2,k) !in new time backward sasha for J-1
           if (segment%oblique) then
             if (dhdt*(segment%grad_normal(I,1,k) + segment%grad_normal(I-1,1,k)) > 0.0) then
               dhdx = segment%grad_normal(I-1,1,k)
             elseif (dhdt*(segment%grad_normal(I,1,k) + segment%grad_normal(I-1,1,k)) == 0.0) then
               dhdx = 0.0
             else
               dhdx = segment%grad_normal(I,1,k)
             endif
           endif
           if (dhdt*dhdy < 0.0) dhdt = 0.0
           if (dhdy == 0.0) dhdy=eps  ! avoid segv
           Cy = min(dhdt/dhdy,rx_max) ! default to normal flow only
           Cx = 0
           cff = max(dhdy*dhdy, eps)
           if (segment%oblique) then
             cff = max(dhdx*dhdx + dhdy*dhdy, eps)
             if (dhdx==0.) dhdx=eps ! avoid segv
             Cx = min(cff,max(dhdt/dhdx,-cff))
           endif
           v_new(i,J,k) = ((cff*v_old(i,J,k) + Cy*v_new(i,J+1,k)) - &
              (max(Cx,0.0)*segment%grad_normal(I-1,2,k) + min(Cx,0.0)*segment%grad_normal(I,2,k))) / (cff + Cy)
         endif
         if ((segment%radiation .or. segment%legacy) .and. segment%nudged) then
           if (dhdt*dhdy < 0.0) then
             tau = segment%Tnudge_in
           else
             tau = segment%Tnudge_out
           endif
           v_new(i,J,k) = v_new(i,J,k) + dt*tau*(segment%normal_vel(i,J,k) - v_old(i,J,k))
!          v_new(i,J,k) = v_new(i,J,k) + dt*tau*(OBC%v(i,J,k) - v_old(i,J,k))
         endif
       enddo; enddo
     end if
  enddo

  call cpu_clock_begin(id_clock_pass)
  call pass_vector(u_new, v_new, G%Domain)
!  call pass_var(h_new, G%Domain)
  call cpu_clock_end(id_clock_pass)

end subroutine radiation_open_bdry_conds


!< Calculate the tangential gradient of the normal flow at the boundary q-points.
subroutine gradient_at_q_points(G,segment,uvel,vvel)
  type(ocean_grid_type), intent(in) :: G !< Ocean grid structure
  type(OBC_segment_type), pointer :: segment !< OBC segment structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: uvel !< zonal velocity
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: vvel !< meridional velocity
  integer :: i,j,k

  if (.not. segment%on_pe) return

  if (segment%is_E_or_W) then

    if (.not.ASSOCIATED(segment%grad_normal)) then
      allocate(segment%grad_normal(segment%HI%JscB:segment%HI%JecB,2,G%ke))
    endif

    if (segment%direction == OBC_DIRECTION_E) then
      I=segment%HI%iscB
      do k=1,G%ke
        do J=segment%HI%JscB,segment%HI%JecB
          segment%grad_normal(J,1,k) = uvel(I-1,j+1,k)-uvel(I-1,j,k)
          segment%grad_normal(J,2,k) = uvel(I,j+1,k)-uvel(I,j,k)
        enddo
      enddo
    else ! western segment
      I=segment%HI%iscB
      do k=1,G%ke
        do J=segment%HI%JscB,segment%HI%JecB
          segment%grad_normal(J,1,k) = uvel(I+1,j+1,k)-uvel(I+1,j,k)
          segment%grad_normal(J,2,k) = uvel(I,j+1,k)-uvel(I,j,k)
        enddo
      enddo
    endif
  else if (segment%is_N_or_S) then

    if (.not.ASSOCIATED(segment%grad_normal)) then
      allocate(segment%grad_normal(segment%HI%IscB:segment%HI%IecB,2,G%ke))
    endif

    if (segment%direction == OBC_DIRECTION_N) then
      J=segment%HI%jscB
      do k=1,G%ke
        do I=segment%HI%IscB,segment%HI%IecB
          segment%grad_normal(I,1,k) = vvel(i+1,J-1,k)-vvel(i,J-1,k)
          segment%grad_normal(I,2,k) = vvel(i+1,J,k)-vvel(i,J,k)
        enddo
      enddo
    else ! south segment
      J=segment%HI%jscB
      do k=1,G%ke
        do I=segment%HI%IscB,segment%HI%IecB
          segment%grad_normal(I,1,k) = vvel(i+1,J+1,k)-vvel(i,J+1,k)
          segment%grad_normal(I,2,k) = vvel(i+1,J,k)-vvel(i,J,k)
        enddo
      enddo
    endif
  endif

end subroutine gradient_at_q_points


!> Sets the initial definitions of the characteristic open boundary conditions.
!! \author Mehmet Ilicak
subroutine set_Flather_data(OBC, tv, h, G, PF, tracer_Reg)
  type(ocean_grid_type),                     intent(inout) :: G !< Ocean grid structure
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary structure
  type(thermo_var_ptrs),                     intent(inout) :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(inout) :: h !< Thickness
  type(param_file_type),                     intent(in)    :: PF !< Parameter file handle
  type(tracer_registry_type),                pointer       :: tracer_Reg !< Tracer registry
  ! Local variables
  logical :: read_OBC_eta = .false.
  logical :: read_OBC_uv = .false.
  logical :: read_OBC_TS = .false.
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: isd_off, jsd_off
  integer :: IsdB, IedB, JsdB, JedB
  character(len=40)  :: mod = "set_Flather_data" ! This subroutine's name.
  character(len=200) :: filename, OBC_file, inputdir ! Strings for file/path

  real :: temp_u(G%domain%niglobal+1,G%domain%njglobal)
  real :: temp_v(G%domain%niglobal,G%domain%njglobal+1)

  real, pointer, dimension(:,:,:) :: &
    OBC_T_u => NULL(), &    ! These arrays should be allocated and set to
    OBC_T_v => NULL(), &    ! specify the values of T and S that should come
    OBC_S_u => NULL(), &
    OBC_S_v => NULL()

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call get_param(PF, mod, "READ_OBC_UV", read_OBC_uv, &
                 "If true, read the values for the velocity open boundary \n"//&
                 "conditions from the file specified by OBC_FILE.", &
                 default=.false.)
  call get_param(PF, mod, "READ_OBC_ETA", read_OBC_eta, &
                 "If true, read the values for the sea surface height \n"//&
                 "open boundary conditions from the file specified by \n"//&
                 "OBC_FILE.", default=.false.)
  call get_param(PF, mod, "READ_OBC_TS", read_OBC_TS, &
                 "If true, read the values for the temperature and \n"//&
                 "salinity open boundary conditions from the file \n"//&
                 "specified by OBC_FILE.", default=.false.)
  if (read_OBC_uv .or. read_OBC_eta .or. read_OBC_TS) then
    call get_param(PF, mod, "OBC_FILE", OBC_file, &
                 "The file from which the appropriate open boundary \n"//&
                 "condition values are read.", default="MOM_OBC_FILE.nc")
    call get_param(PF, mod, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    filename = trim(inputdir)//trim(OBC_file)
    call log_param(PF, mod, "INPUTDIR/OBC_FILE", filename)
  endif

  if (open_boundary_query(OBC, apply_Flather_OBC=.true.)) then
    if (.not.associated(OBC%vbt_outer)) then
      allocate(OBC%vbt_outer(isd:ied,JsdB:JedB)) ; OBC%vbt_outer(:,:) = 0.0
    endif

    if (.not.associated(OBC%ubt_outer)) then
      allocate(OBC%ubt_outer(IsdB:IedB,jsd:jed)) ; OBC%ubt_outer(:,:) = 0.0
    endif

    if (.not.associated(OBC%eta_outer_u)) then
      allocate(OBC%eta_outer_u(IsdB:IedB,jsd:jed)) ; OBC%eta_outer_u(:,:) = 0.0
    endif

    if (.not.associated(OBC%eta_outer_v)) then
      allocate(OBC%eta_outer_v(isd:ied,JsdB:JedB)) ; OBC%eta_outer_v(:,:) = 0.0
    endif

    if (read_OBC_uv) then
      call read_data(filename, 'ubt', OBC%ubt_outer, &
                     domain=G%Domain%mpp_domain, position=EAST_FACE)
      call read_data(filename, 'vbt', OBC%vbt_outer, &
                     domain=G%Domain%mpp_domain, position=NORTH_FACE)
    endif

    if (read_OBC_eta) then
      call read_data(filename, 'eta_outer_u', OBC%eta_outer_u, &
                     domain=G%Domain%mpp_domain, position=EAST_FACE)
      call read_data(filename, 'eta_outer_v', OBC%eta_outer_v, &
                     domain=G%Domain%mpp_domain, position=NORTH_FACE)
    endif

    call pass_vector(OBC%eta_outer_u,OBC%eta_outer_v,G%Domain, To_All+SCALAR_PAIR, CGRID_NE)
    call pass_vector(OBC%ubt_outer,OBC%vbt_outer,G%Domain)
  endif

  ! Define radiation coefficients r[xy]_old_[uvh] as needed.  For now, there are
  ! no radiation conditions applied to the thicknesses, since the thicknesses
  ! might not be physically motivated.  Instead, sponges should be used to
  ! enforce the near-boundary layer structure.
  if (OBC%nudged_u_BCs_exist_globally) then
    allocate(OBC%u(IsdB:IedB,jsd:jed,nz)) ; OBC%u(:,:,:) = 0.0
    allocate(OBC%uh(IsdB:IedB,jsd:jed,nz)) ; OBC%uh(:,:,:) = 0.0
  endif
  if (OBC%nudged_v_BCs_exist_globally) then
    allocate(OBC%v(isd:ied,JsdB:JedB,nz)) ; OBC%v(:,:,:) = 0.0
    allocate(OBC%vh(isd:ied,JsdB:JedB,nz)) ; OBC%vh(:,:,:) = 0.0
  endif

  if (associated(tv%T)) then
    allocate(OBC_T_u(IsdB:IedB,jsd:jed,nz)) ; OBC_T_u(:,:,:) = 0.0
    allocate(OBC_S_u(IsdB:IedB,jsd:jed,nz)) ; OBC_S_u(:,:,:) = 0.0
    allocate(OBC_T_v(isd:ied,JsdB:JedB,nz)) ; OBC_T_v(:,:,:) = 0.0
    allocate(OBC_S_v(isd:ied,JsdB:JedB,nz)) ; OBC_S_v(:,:,:) = 0.0

    if (read_OBC_TS) then
      call read_data(filename, 'OBC_T_u', OBC_T_u, &
                     domain=G%Domain%mpp_domain, position=EAST_FACE)
      call read_data(filename, 'OBC_S_u', OBC_S_u, &
                     domain=G%Domain%mpp_domain, position=EAST_FACE)

      call read_data(filename, 'OBC_T_v', OBC_T_v, &
                     domain=G%Domain%mpp_domain, position=NORTH_FACE)
      call read_data(filename, 'OBC_S_v', OBC_S_v, &
                     domain=G%Domain%mpp_domain, position=NORTH_FACE)
    else
      call pass_var(tv%T, G%Domain)
      call pass_var(tv%S, G%Domain)
      do k=1,nz ; do j=js,je ; do I=is-1,ie
        if (OBC%OBC_segment_u(I,j) /= OBC_NONE) then
          if (OBC%segment(OBC%OBC_segment_u(I,j))%direction == OBC_DIRECTION_E) then
            OBC_T_u(I,j,k) = tv%T(i,j,k)
            OBC_S_u(I,j,k) = tv%S(i,j,k)
          elseif (OBC%segment(OBC%OBC_segment_u(I,j))%direction == OBC_DIRECTION_W) then
            OBC_T_u(I,j,k) = tv%T(i+1,j,k)
            OBC_S_u(I,j,k) = tv%S(i+1,j,k)
          elseif (G%mask2dT(i,j) + G%mask2dT(i+1,j) > 0) then
            OBC_T_u(I,j,k) = (G%mask2dT(i,j)*tv%T(i,j,k) + G%mask2dT(i+1,j)*tv%T(i+1,j,k)) / &
                             (G%mask2dT(i,j) + G%mask2dT(i+1,j))
            OBC_S_u(I,j,k) = (G%mask2dT(i,j)*tv%S(i,j,k) + G%mask2dT(i+1,j)*tv%S(i+1,j,k)) / &
                             (G%mask2dT(i,j) + G%mask2dT(i+1,j))
          else ! This probably shouldn't happen or maybe it doesn't matter?
            OBC_T_u(I,j,k) = 0.5*(tv%T(i,j,k)+tv%T(i+1,j,k))
            OBC_S_u(I,j,k) = 0.5*(tv%S(i,j,k)+tv%S(i+1,j,k))
          endif
        else
          OBC_T_u(I,j,k) = 0.5*(tv%T(i,j,k)+tv%T(i+1,j,k))
          OBC_S_u(I,j,k) = 0.5*(tv%S(i,j,k)+tv%S(i+1,j,k))
        endif
      enddo; enddo ; enddo

      do k=1,nz ; do J=js-1,je ; do i=is,ie
        if (OBC%OBC_segment_v(i,J) /= OBC_NONE) then
          if (OBC%segment(OBC%OBC_segment_v(i,J))%direction == OBC_DIRECTION_N) then
            OBC_T_v(i,J,k) = tv%T(i,j,k)
            OBC_S_v(i,J,k) = tv%S(i,j,k)
          elseif (OBC%segment(OBC%OBC_segment_v(i,J))%direction == OBC_DIRECTION_S) then
            OBC_T_v(i,J,k) = tv%T(i,j+1,k)
            OBC_S_v(i,J,k) = tv%S(i,j+1,k)
          elseif (G%mask2dT(i,j) + G%mask2dT(i,j+1) > 0) then
            OBC_T_v(i,J,k) = (G%mask2dT(i,j)*tv%T(i,j,k) + G%mask2dT(i,j+1)*tv%T(i,j+1,k)) / &
                             (G%mask2dT(i,j) + G%mask2dT(i,j+1))
            OBC_S_v(i,J,k) = (G%mask2dT(i,j)*tv%S(i,j,k) + G%mask2dT(i,j+1)*tv%S(i,j+1,k)) / &
                             (G%mask2dT(i,j) + G%mask2dT(i,j+1))
          else ! This probably shouldn't happen or maybe it doesn't matter?
            OBC_T_v(i,J,k) = 0.5*(tv%T(i,j,k)+tv%T(i,j+1,k))
            OBC_S_v(i,J,k) = 0.5*(tv%S(i,j,k)+tv%S(i,j+1,k))
          endif
        else
          OBC_T_v(i,J,k) = 0.5*(tv%T(i,j,k)+tv%T(i,j+1,k))
          OBC_S_v(i,J,k) = 0.5*(tv%S(i,j,k)+tv%S(i,j+1,k))
        endif
      enddo; enddo ; enddo
    endif

    call pass_vector(OBC_T_u, OBC_T_v, G%Domain, To_All+SCALAR_PAIR, CGRID_NE)
    call pass_vector(OBC_S_u, OBC_S_v, G%Domain, To_All+SCALAR_PAIR, CGRID_NE)

    call add_tracer_OBC_values("T", tracer_Reg, OBC_in_u=OBC_T_u, &
                                                OBC_in_v=OBC_T_v)
    call add_tracer_OBC_values("S", tracer_Reg, OBC_in_u=OBC_S_u, &
                                                OBC_in_v=OBC_S_v)
    do k=1,nz ; do j=jsd,jed ; do I=isd,ied-1
      if (OBC%segment(OBC%OBC_segment_u(I,j))%direction == OBC_DIRECTION_E) then
        tv%T(i+1,j,k) = tv%T(i,j,k) ; tv%S(i+1,j,k) = tv%S(i,j,k)
      elseif (OBC%segment(OBC%OBC_segment_u(I,j))%direction == OBC_DIRECTION_W) then
        tv%T(i,j,k) = tv%T(i+1,j,k) ; tv%S(i,j,k) = tv%S(i+1,j,k)
      endif
    enddo ; enddo ; enddo
    do k=1,nz ; do J=jsd,jed-1 ; do i=isd,ied
      if (OBC%segment(OBC%OBC_segment_v(i,J))%direction == OBC_DIRECTION_N) then
        tv%T(i,j+1,k) = tv%T(i,j,k) ; tv%S(i,j+1,k) = tv%S(i,j,k)
      elseif (OBC%segment(OBC%OBC_segment_v(i,J))%direction == OBC_DIRECTION_S) then
        tv%T(i,j,k) = tv%T(i,j+1,k) ; tv%S(i,j,k) = tv%S(i,j+1,k)
      endif
    enddo ; enddo ; enddo
  endif

  do k=1,nz ; do j=jsd,jed ; do I=isd,ied-1
    if (OBC%segment(OBC%OBC_segment_u(I,j))%direction == OBC_DIRECTION_E) &
                        h(i+1,j,k) = h(i,j,k)
    if (OBC%segment(OBC%OBC_segment_u(I,j))%direction == OBC_DIRECTION_W) &
                        h(i,j,k) = h(i+1,j,k)
  enddo ; enddo ; enddo
  do k=1,nz ; do J=jsd,jed-1 ; do i=isd,ied
    if (OBC%segment(OBC%OBC_segment_v(i,J))%direction == OBC_DIRECTION_N) &
                       h(i,j+1,k) = h(i,j,k)
    if (OBC%segment(OBC%OBC_segment_v(i,J))%direction == OBC_DIRECTION_S) &
                       h(i,j,k) = h(i,j+1,k)
  enddo ; enddo ; enddo
! When we do not extend segments, this commented block was needed to
! get the same'ish h's.
! do k=1,nz ; do j=jsd,jed-1 ; do i=isd,ied
!   if (OBC%OBC_direction_v(i,J) == OBC_DIRECTION_N) h(i,j+1,k) = h(i,j,k)
! enddo ; enddo ; enddo
! do k=1,nz ; do j=jsd+1,jed ; do i=isd,ied
!   if (OBC%OBC_direction_v(i,J-1) == OBC_DIRECTION_S) h(i,j-1,k) = h(i,j,k)
! enddo ; enddo ; enddo
! do k=1,nz ; do j=jsd,jed ; do i=isd,ied-1
!   if (OBC%OBC_direction_u(I,j) == OBC_DIRECTION_E) h(i+1,j,k) = h(i,j,k)
! enddo ; enddo ; enddo
! do k=1,nz ; do j=jsd,jed ; do i=isd+1,ied
!   if (OBC%OBC_direction_u(I-1,j) == OBC_DIRECTION_W) h(i-1,j,k) = h(i,j,k)
! enddo ; enddo ; enddo
! do k=1,nz ; do j=jsd,jed-1 ; do i=isd,ied-1
!   if (OBC%OBC_direction_v(i,J) == OBC_DIRECTION_N .and. &
!       OBC%OBC_direction_u(I,j) == OBC_DIRECTION_E) h(i+1,j+1,k) = h(i,j,k)
! enddo ; enddo ; enddo
! do k=1,nz ; do j=jsd,jed-1 ; do i=isd+1,ied
!   if (OBC%OBC_direction_v(i,J) == OBC_DIRECTION_N .and. &
!       OBC%OBC_direction_u(I-1,j) == OBC_DIRECTION_W) h(i-1,j+1,k) = h(i,j,k)
! enddo ; enddo ; enddo
! do k=1,nz ; do j=jsd+1,jed ; do i=isd,ied-1
!   if (OBC%OBC_direction_v(i,J-1) == OBC_DIRECTION_S .and. &
!       OBC%OBC_direction_u(I,j) == OBC_DIRECTION_E) h(i+1,j-1,k) = h(i,j,k)
! enddo ; enddo ; enddo
! do k=1,nz ; do j=jsd+1,jed ; do i=isd+1,ied
!   if (OBC%OBC_direction_v(i,J-1) == OBC_DIRECTION_S .and. &
!       OBC%OBC_direction_u(I-1,j) == OBC_DIRECTION_W) h(i-1,j-1,k) = h(i,j,k)
! enddo ; enddo ; enddo

end subroutine set_Flather_data

function lookup_seg_field(OBC_seg,field)
  type(OBC_segment_type), pointer :: OBC_seg
  character(len=32), intent(in) :: field ! The field name
  integer :: lookup_seg_field

  integer :: n,m

  lookup_seg_field=-1
  do n=1,OBC_seg%num_fields
   if (trim(field) == OBC_seg%field_names(n)) then
     lookup_seg_field=n
     return
   endif
  enddo

  return

end function lookup_seg_field


subroutine fill_OBC_halos(G, GV, OBC, tv, h, tracer_Reg)

  type(ocean_grid_type),                     intent(inout) :: G !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !<  Ocean vertical grid structure
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary structure
  type(thermo_var_ptrs),                     intent(inout) :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(inout)    :: h !< Thickness
  type(tracer_registry_type),                pointer       :: tracer_Reg !< Tracer registry
  ! Local variables

  integer :: i, j, k, isd, ied, jsd, jed, is, ie, js, je
  integer :: n, nz
  character(len=40)  :: mod = "fill_OBC_halos" ! This subroutine's name.
  type(OBC_segment_type), pointer :: segment
  real, pointer, dimension(:,:,:) :: &
    OBC_T_u => NULL(), &    ! These arrays should be allocated and set to
    OBC_T_v => NULL(), &    ! specify the values of T and S that should come
    OBC_S_u => NULL(), &
    OBC_S_v => NULL()

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  nz=G%ke

  if (.not. associated(OBC) .or. .not. associated(tv%T)) return

  call pass_var(tv%T,G%Domain)
  call pass_var(tv%S,G%Domain)

  do n = 1, OBC%number_of_segments
    segment => OBC%segment(n)

    if (.not. segment%on_pe) cycle ! continue to next segment if not in computational domain

    do j=G%jsdB+1,G%jedB ; do I=G%isdB,G%iedB-1
      if (segment%direction == OBC_DIRECTION_E .and. OBC%OBC_segment_u(I,j) /= OBC_NONE ) then
        do k=1,nz
          tv%T(i+1,j,k) = tv%T(i,j,k) ; tv%S(i+1,j,k) = tv%S(i,j,k); h(i+1,j,k) = h(i,j,k)
        enddo
      endif
    enddo ; enddo

    do j=G%jsdB+1,G%jedB ; do I=G%isdB+1,G%iedB
      if (segment%direction == OBC_DIRECTION_W .and. OBC%OBC_segment_u(I,j) /= OBC_NONE ) then
        do k=1,nz
          tv%T(i-1,j,k) = tv%T(i,j,k) ; tv%S(i-1,j,k) = tv%S(i,j,k); h(i-1,j,k) = h(i,j,k)
        enddo
      endif
    enddo ; enddo

    do j=G%jsdB,G%jedB-1 ; do I=G%isdB+1,G%iedB
      if (segment%direction == OBC_DIRECTION_N .and. OBC%OBC_segment_v(I,j) /= OBC_NONE ) then
        do k=1,nz
          tv%T(i,j+1,k) = tv%T(i,j,k) ; tv%S(i,j+1,k) = tv%S(i,j,k); h(i,j+1,k) = h(i,j,k)
        enddo
      endif
    enddo ; enddo

    do j=G%jsdB+1,G%jedB ; do I=G%isdB+1,G%iedB
      if (segment%direction == OBC_DIRECTION_S .and. OBC%OBC_segment_v(I,j) /= OBC_NONE ) then
        do k=1,nz
          tv%T(i,j-1,k) = tv%T(i,j,k) ; tv%S(i,j-1,k) = tv%S(i,j,k); h(i,j-1,k) = h(i,j,k)
        enddo
      endif
    enddo ; enddo

  enddo

  if (associated(tv%T)) then
     if (.not. associated(OBC_T_u)) then
       allocate(OBC_T_u(G%IsdB:G%IedB,G%jsd:G%jed,nz)) ; OBC_T_u(:,:,:) = 0.0
       allocate(OBC_S_u(G%IsdB:G%IedB,G%jsd:G%jed,nz)) ; OBC_S_u(:,:,:) = 0.0
       allocate(OBC_T_v(G%Isd:G%Ied,G%jsdB:G%jedB,nz)) ; OBC_T_v(:,:,:) = 0.0
       allocate(OBC_S_v(G%Isd:G%Ied,G%jsdB:G%jedB,nz)) ; OBC_S_v(:,:,:) = 0.0
     endif
     do k=1,nz ; do j=G%jsd,G%jed ; do I=G%isd,G%ied-1
        if (OBC%OBC_segment_u(I,j) /= OBC_NONE) then
          if (OBC%segment(OBC%OBC_segment_u(I,j))%direction == OBC_DIRECTION_E) then
            OBC_T_u(I,j,k) = tv%T(i,j,k)
            OBC_S_u(I,j,k) = tv%S(i,j,k)
          elseif (OBC%segment(OBC%OBC_segment_u(I,j))%direction == OBC_DIRECTION_W) then
            OBC_T_u(I,j,k) = tv%T(i+1,j,k)
            OBC_S_u(I,j,k) = tv%S(i+1,j,k)
          elseif (G%mask2dT(i,j) + G%mask2dT(i+1,j) > 0) then
            OBC_T_u(I,j,k) = (G%mask2dT(i,j)*tv%T(i,j,k) + G%mask2dT(i+1,j)*tv%T(i+1,j,k)) / &
                             (G%mask2dT(i,j) + G%mask2dT(i+1,j))
            OBC_S_u(I,j,k) = (G%mask2dT(i,j)*tv%S(i,j,k) + G%mask2dT(i+1,j)*tv%S(i+1,j,k)) / &
                             (G%mask2dT(i,j) + G%mask2dT(i+1,j))
          else ! This probably shouldn't happen or maybe it doesn't matter?
            OBC_T_u(I,j,k) = 0.5*(tv%T(i,j,k)+tv%T(i+1,j,k))
            OBC_S_u(I,j,k) = 0.5*(tv%S(i,j,k)+tv%S(i+1,j,k))
          endif
        else
          OBC_T_u(I,j,k) = 0.5*(tv%T(i,j,k)+tv%T(i+1,j,k))
          OBC_S_u(I,j,k) = 0.5*(tv%S(i,j,k)+tv%S(i+1,j,k))
        endif
     enddo; enddo ; enddo

     do k=1,nz ; do J=G%jsd,G%jed-1 ; do i=G%isd,G%ied
        if (OBC%OBC_segment_v(i,J) /= OBC_NONE) then
          if (OBC%segment(OBC%OBC_segment_v(i,J))%direction == OBC_DIRECTION_N) then
            OBC_T_v(i,J,k) = tv%T(i,j,k)
            OBC_S_v(i,J,k) = tv%S(i,j,k)
          elseif (OBC%segment(OBC%OBC_segment_v(i,J))%direction == OBC_DIRECTION_S) then
            OBC_T_v(i,J,k) = tv%T(i,j+1,k)
            OBC_S_v(i,J,k) = tv%S(i,j+1,k)
          elseif (G%mask2dT(i,j) + G%mask2dT(i,j+1) > 0) then
            OBC_T_v(i,J,k) = (G%mask2dT(i,j)*tv%T(i,j,k) + G%mask2dT(i,j+1)*tv%T(i,j+1,k)) / &
                             (G%mask2dT(i,j) + G%mask2dT(i,j+1))
            OBC_S_v(i,J,k) = (G%mask2dT(i,j)*tv%S(i,j,k) + G%mask2dT(i,j+1)*tv%S(i,j+1,k)) / &
                             (G%mask2dT(i,j) + G%mask2dT(i,j+1))
          else ! This probably shouldn't happen or maybe it doesn't matter?
            OBC_T_v(i,J,k) = 0.5*(tv%T(i,j,k)+tv%T(i,j+1,k))
            OBC_S_v(i,J,k) = 0.5*(tv%S(i,j,k)+tv%S(i,j+1,k))
          endif
        else
          OBC_T_v(i,J,k) = 0.5*(tv%T(i,j,k)+tv%T(i,j+1,k))
          OBC_S_v(i,J,k) = 0.5*(tv%S(i,j,k)+tv%S(i,j+1,k))
        endif
     enddo; enddo ; enddo

    call pass_vector(OBC_T_u, OBC_T_v, G%Domain, To_All+SCALAR_PAIR, CGRID_NE)
    call pass_vector(OBC_S_u, OBC_S_v, G%Domain, To_All+SCALAR_PAIR, CGRID_NE)

    call add_tracer_OBC_values("T",tracer_Reg, OBC_in_u=OBC_T_u,OBC_in_v=OBC_T_v)
    call add_tracer_OBC_values("S",tracer_Reg, OBC_in_u=OBC_S_u,OBC_in_v=OBC_S_v)


  endif
end subroutine fill_OBC_halos

!> Allocate segment data fields
subroutine allocate_OBC_segment_data(OBC, segment)
  type(ocean_OBC_type),   pointer       :: OBC     !< Open boundary structure
  type(OBC_segment_type), intent(inout) :: segment !< Open boundary segment
  ! Local variables
  integer :: isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB
  character(len=40)  :: mod = "allocate_OBC_segment_data" ! This subroutine's name.

  isd = segment%HI%isd ; ied = segment%HI%ied
  jsd = segment%HI%jsd ; jed = segment%HI%jed
  IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
  JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB

  if (.not. segment%on_pe) return

  if (segment%is_E_or_W) then
    ! If these are just Flather, change update_OBC_segment_data accordingly
    allocate(segment%Cg(IsdB:IedB,jsd:jed));                    segment%Cg(:,:)=0.
    allocate(segment%Htot(IsdB:IedB,jsd:jed));                  segment%Htot(:,:)=0.0
    allocate(segment%h(IsdB:IedB,jsd:jed,OBC%ke));              segment%h(:,:,:)=0.0
    if (segment%Flather) then
!     allocate(segment%e(IsdB:IedB,jsd:jed,OBC%ke));              segment%e(:,:,:)=0.0
      allocate(segment%eta(IsdB:IedB,jsd:jed));                   segment%eta(:,:)=0.0
      allocate(segment%rx_normal(IsdB:IedB,jsd:jed,OBC%ke));      segment%rx_normal(:,:,:)=0.0
!     allocate(segment%normal_trans_bt(IsdB:IedB,jsd:jed));       segment%normal_trans_bt(:,:)=0.0
    endif
    if (segment%Flather .or. segment%specified) then
      allocate(segment%normal_vel_bt(IsdB:IedB,jsd:jed));         segment%normal_vel_bt(:,:)=0.0
    endif
    if (segment%nudged .or. segment%specified) then
      allocate(segment%normal_vel(IsdB:IedB,jsd:jed,OBC%ke));     segment%normal_vel(:,:,:)=0.0
    endif
    if (segment%specified) then
      allocate(segment%normal_trans(IsdB:IedB,jsd:jed,OBC%ke));   segment%normal_trans(:,:,:)=0.0
    endif
!   allocate(segment%tangent_vel(IsdB:IedB,JsdB:JedB,OBC%ke));  segment%tangent_vel(:,:,:)=0.0
!   allocate(segment%tangent_vel_bt(IsdB:IedB,JsdB:JedB));      segment%tangent_vel_bt(:,:)=0.0
  else
    allocate(segment%Cg(isd:ied,JsdB:JedB));                    segment%Cg(:,:)=0.
    allocate(segment%Htot(isd:ied,JsdB:JedB));                  segment%Htot(:,:)=0.0
    allocate(segment%h(isd:ied,JsdB:JedB,OBC%ke));              segment%h(:,:,:)=0.0
    if (segment%Flather) then
!     allocate(segment%e(isd:ied,JsdB:JedB,OBC%ke));              segment%e(:,:,:)=0.0
      allocate(segment%eta(isd:ied,JsdB:JedB));                   segment%eta(:,:)=0.0
      allocate(segment%rx_normal(isd:ied,JsdB:JedB,OBC%ke));      segment%rx_normal(:,:,:)=0.0
!     allocate(segment%normal_trans_bt(isd:ied,JsdB:JedB));       segment%normal_trans_bt(:,:)=0.0
    endif
    if (segment%Flather .or. segment%specified) then
      allocate(segment%normal_vel_bt(isd:ied,JsdB:JedB));         segment%normal_vel_bt(:,:)=0.0
    endif
    if (segment%nudged .or. segment%specified) then
      allocate(segment%normal_vel(isd:ied,JsdB:JedB,OBC%ke));     segment%normal_vel(:,:,:)=0.0
    endif
    if (segment%specified) then
      allocate(segment%normal_trans(isd:ied,JsdB:JedB,OBC%ke));   segment%normal_trans(:,:,:)=0.0
    endif
!   allocate(segment%tangent_vel(IsdB:IedB,JsdB:JedB,OBC%ke));  segment%tangent_vel(:,:,:)=0.0
!   allocate(segment%tangent_vel_bt(IsdB:IedB,JsdB:JedB));      segment%tangent_vel_bt(:,:)=0.0
  endif
end subroutine allocate_OBC_segment_data

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
  real, parameter :: silly_value = 1.E40

  if (.not. associated(OBC)) return

  do n = 1, OBC%number_of_segments
    do k = 1, G%ke
      if (OBC%segment(n)%is_N_or_S) then
        J = OBC%segment(n)%HI%JsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
          do I = OBC%segment(n)%HI%IsdB, OBC%segment(n)%HI%IedB
            u(I,j+1,k) = silly_value
          enddo
        else
          do I = OBC%segment(n)%HI%IsdB, OBC%segment(n)%HI%IedB
            u(I,j,k) = silly_value
          enddo
        endif
      elseif (OBC%segment(n)%is_E_or_W) then
        I = OBC%segment(n)%HI%IsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
          do J = OBC%segment(n)%HI%JsdB, OBC%segment(n)%HI%JedB
            v(i+1,J,k) = silly_value
          enddo
        else
          do J = OBC%segment(n)%HI%JsdB, OBC%segment(n)%HI%JedB
            v(i,J,k) = silly_value
          enddo
        endif
      endif
    enddo
  enddo

end subroutine open_boundary_test_extern_uv

subroutine update_OBC_segment_data(G, GV, OBC, tv, h, Time)
  type(ocean_grid_type),                     intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV  !<  Ocean vertical grid structure
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundary structure
  type(thermo_var_ptrs),                     intent(in)    :: tv  !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(inout) :: h   !< Thickness
! real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(inout) :: e   !< Layer interface height
! real, dimension(SZI_(G),SZJ_(G))         , intent(inout) :: eta !< Thickness
  type(time_type),                           intent(in)    :: Time
  ! Local variables

  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB, n, m, nz
  character(len=40)  :: mod = "set_OBC_segment_data" ! This subroutine's name.
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
    is_obc = max(segment%is_obc,is-1)
    ie_obc = min(segment%ie_obc,ie)
    js_obc = max(segment%js_obc,js-1)
    je_obc = min(segment%je_obc,je)


    if (segment%is_E_or_W) then
      nj_seg=nj_seg-1
      js_obc=js_obc+1
    else
      ni_seg=ni_seg-1
      is_obc=is_obc+1
    endif

!    do j=jsd,jed ; do I=isd,ied-1
!      if (segment%direction == OBC_DIRECTION_E .and. OBC%OBC_segment_u(I,j) /= OBC_NONE ) then
!        do k=1,nz
!          tv%T(i+1,j,k) = tv%T(i,j,k) ; tv%S(i+1,j,k) = tv%S(i,j,k); h(i+1,j,k) = h(i,j,k)
!        enddo
!      else if (segment%direction == OBC_DIRECTION_W .and. OBC%OBC_segment_u(I,j) /= OBC_NONE ) then
!        tv%T(i,j,k) = tv%T(i+1,j,k) ; tv%S(i,j,k) = tv%S(i+1,j,k); h(i,j,k) = h(i+1,j,k)
!      endif
!    enddo ; enddo

!    do j=jsd,jed-1 ; do I=isd,ied-1
!      if (segment%direction == OBC_DIRECTION_N .and. OBC%OBC_segment_v(I,j) /= OBC_NONE ) then
!        do k=1,nz
!          tv%T(i,j+1,k) = tv%T(i,j,k) ; tv%S(i,j+1,k) = tv%S(i,j,k); h(i,j+1,k) = h(i,j,k)
!        enddo
!      else if (segment%direction == OBC_DIRECTION_S .and. OBC%OBC_segment_v(I,j) /= OBC_NONE ) then
!        tv%T(i,j,k) = tv%T(i,j+1,k) ; tv%S(i,j,k) = tv%S(i,j+1,k); h(i,j,k) = h(i,j+1,k)
!      endif
!    enddo ; enddo

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
    if (segment%direction == OBC_DIRECTION_W .or. segment%direction == OBC_DIRECTION_E) then
      if (segment%direction == OBC_DIRECTION_W) ishift=1
      do j=js_obc,je_obc
        do I=Is_obc,Ie_obc
!        i2 =  segment%Is_obc + i - 1
!        j2 =  segment%Js_obc + j - 1
!        if ((i2 .gt. ied .or. i2 .lt. isd) .or. (j2 .gt. jed .or. j2 .lt. jsd)) cycle
!        if (OBC%OBC_segment_u(i2,j2) /= n) cycle
          segment%Cg(I,j) = sqrt(GV%g_prime(1)*G%bathyT(i+ishift,j))
          if (GV%Boussinesq) then
            segment%Htot(I,j) = G%bathyT(i+ishift,j)*GV%m_to_H! + eta(i+ishift,j)
!          else
!            segment%Htot(I,j) =  eta(i+ishift,j)
          endif
          do k=1,G%ke
            segment%h(I,j,k) = h(i+ishift,j,k)
!           segment%e(I,j,k) = e(i+ishift,j,k)
          enddo
        enddo
      enddo

    else! (segment%direction == OBC_DIRECTION_N .or. segment%direction == OBC_DIRECTION_S)
      if (segment%direction == OBC_DIRECTION_S) jshift=1

      do J=Js_obc,Je_obc
        do i=is_obc,ie_obc
          segment%Cg(i,J) = sqrt(GV%g_prime(1)*G%bathyT(i,j+jshift))
          if (GV%Boussinesq) then
            segment%Htot(i,J) = G%bathyT(i,j+jshift)*GV%m_to_H! + eta(i,j+jshift)
!          else
!            segment%Htot(i,J) = eta(i,j+jshift)
          endif
          do k=1,G%ke
            segment%h(i,J,k) = h(i,j+jshift,k)
!           segment%e(i,J,k) = e(i,j+jshift,k)
          enddo
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
          if (trim(segment%field(m)%name) == 'UO' .or. trim(segment%field(m)%name) == 'VO') then
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
                   call remapping_core_h(segment%field(m)%nk_src,segment%field(m)%dz_src(i,j,:),&
                        segment%field(m)%buffer_src(i,j,:),G%ke, h(i,j,:),&
                        segment%field(m)%buffer_dst(i,j,:),OBC%remap_CS)
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
           if (trim(segment%field(m)%name) == 'UO' .or. trim(segment%field(m)%name) == 'VO') then
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
    enddo


  enddo ! end segment loop

end subroutine update_OBC_segment_data

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
