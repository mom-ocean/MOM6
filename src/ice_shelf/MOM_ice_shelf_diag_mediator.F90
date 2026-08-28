! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> The subroutines here provide convenient wrappers to the FMS diag_manager
!! interfaces with additional diagnostic capabilities.
module MOM_IS_diag_mediator

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_checksums,          only : chksum0, hchksum, uchksum, vchksum, Bchksum
use MOM_coms,               only : PE_here
use MOM_cpu_clock,          only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,          only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_manager_infra, only : MOM_diag_manager_init
use MOM_diag_manager_infra, only : MOM_diag_axis_init, get_MOM_diag_axis_name
use MOM_diag_manager_infra, only : send_data_infra, MOM_diag_field_add_attribute, EAST, NORTH
use MOM_diag_manager_infra, only : register_diag_field_infra, register_static_field_infra
use MOM_diag_manager_infra, only : get_MOM_diag_field_id, DIAG_FIELD_NOT_FOUND
use MOM_diag_manager_infra, only : diag_send_complete_infra
use MOM_error_handler,      only : MOM_error, FATAL, is_root_pe, assert, callTree_showQuery
use MOM_error_handler,      only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,        only : get_param, log_version, param_file_type
use MOM_grid,               only : ocean_grid_type
use MOM_io,                 only : get_filename_appendix
use MOM_safe_alloc,         only : safe_alloc_ptr, safe_alloc_alloc
use MOM_string_functions,   only : lowercase, uppercase, slasher, ints_to_string, trim_trailing_commas
use MOM_time_manager,       only : time_type, get_time
use MOM_unit_scaling,       only : unit_scale_type

implicit none ; private

public MOM_IS_diag_mediator_infrastructure_init
public MOM_IS_diag_mediator_init, MOM_IS_diag_mediator_end, set_IS_diag_mediator_grid
public set_IS_axes_info, MOM_diag_axis_init
public register_MOM_IS_diag_field, register_MOM_IS_static_field, register_MOM_IS_scalar_field
public post_IS_data, post_IS_data_0d, MOM_IS_diag_send_complete
public safe_alloc_ptr, safe_alloc_alloc, time_type
public enable_averaging, enable_averages, disable_averaging, query_averaging_enabled
public MOM_IS_diag_mediator_close_registration, get_diag_time_end
public define_axes_group, diag_masks_set
public diag_register_area_ids, found_in_diagtable

!> Make a diagnostic available for averaging or output.
interface post_IS_data
  module procedure post_IS_data_2d, post_IS_data_0d
end interface post_IS_data

!> Registers a non-array scalar diagnostic, returning an integer handle
interface register_MOM_IS_scalar_field
  module procedure register_scalar_field_CS, register_scalar_field_axes
end interface register_MOM_IS_scalar_field

!> A group of 1D axes that comprise a 1D/2D/3D mesh
type, public :: axes_grp
  character(len=15) :: id   !< The id string for this particular combination of handles.
  integer           :: rank !< Number of dimensions in the list of axes.
  integer, dimension(:), allocatable :: handles !< Handles to 1D axes.
  type(diag_ctrl), pointer :: diag_cs => null() !< Circular link back to the main diagnostics control structure
                                                !! (Used to avoid passing said structure into every possible call).
  ! ID's for cell_methods
  character(len=9) :: x_cell_method = '' !< Default nature of data representation, if axes group
                                         !! includes x-direction.
  character(len=9) :: y_cell_method = '' !< Default nature of data representation, if axes group
                                         !! includes y-direction.
  ! For detecting position on the grid
  logical :: is_h_point = .false. !< If true, indicates that this axes group is for an h-point located field.
  logical :: is_q_point = .false. !< If true, indicates that this axes group is for a q-point located field.
  logical :: is_u_point = .false. !< If true, indicates that this axes group is for a u-point located field.
  logical :: is_v_point = .false. !< If true, indicates that this axes group is for a v-point located field.

  ! ID's for cell_measures
  integer :: id_area = -1 !< The diag_manager id for area to be used for cell_measure of variables with this axes_grp.
  ! For masking
  real, pointer, dimension(:,:)   :: mask2d => null() !< Mask for 2d (x-y) axes [nondim]
  real, pointer, dimension(:,:)   :: mask2d_comp => null() !< Mask for 2-d axes on the computational
                                                      !! domain for this diagnostic [nondim]
end type axes_grp

!> This type is used to represent a diagnostic at the diag_mediator level.
!!
!! There can be both 'primary' and 'secondary' diagnostics. The primaries
!! reside in the diag_cs%diags array. They have an id which is an index
!! into this array. The secondaries are 'variations' on the primary diagnostic.
!! For example the CMOR diagnostics are secondary. The secondary diagnostics
!! are kept in a list with the primary diagnostic as the head.
type, private :: diag_type
  logical :: in_use              !< True if this entry is being used.
  integer :: fms_diag_id         !< Underlying FMS diag_manager id.
  character(len=64) :: debug_str = '' !< The diagnostic name and module for FATAL errors and debugging.
  type(axes_grp), pointer :: axes => null() !< The axis group for this diagnostic
  type(diag_type), pointer :: next => null() !< Pointer to the next diagnostic
  real :: conversion_factor = 0. !< If non-zero, a factor to multiply data by before posting to FMS,
                                 !! often including factors to undo internal scaling in units of [a A-1 ~> 1]
end type diag_type

!>   The diag_ctrl data type contains times to regulate diagnostics along with masks and
!! axes to use with diagnostics, and a list of structures with data about each diagnostic.
type, public :: diag_ctrl
  integer :: available_diag_doc_unit = -1 !< The unit number of a diagnostic documentation file.
                                          !! This file is open if available_diag_doc_unit is > 0.
  integer :: chksum_iounit = -1           !< The unit number of a diagnostic documentation file.
                                          !! This file is open if available_diag_doc_unit is > 0.
  logical :: diag_as_chksum  !< If true, log chksums in a text file instead of posting diagnostics
  logical :: show_call_tree  !< Display the call tree while running. Set by VERBOSITY level.
  logical :: index_space_axes !< If true, diagnostic horizontal coordinates axes are in index space.

  ! The following fields are used for the output of the data.
  ! These give the computational-domain sizes, and are relative to a start value
  ! of 1 in memory for the tracer-point arrays.
  integer :: is  !< The start i-index of cell centers within the computational domain
  integer :: ie  !< The end i-index of cell centers within the computational domain
  integer :: js  !< The start j-index of cell centers within the computational domain
  integer :: je  !< The end j-index of cell centers within the computational domain
  ! These give the memory-domain sizes, and can be start at any value on each PE.
  integer :: isd !< The start i-index of cell centers within the data domain
  integer :: ied !< The end i-index of cell centers within the data domain
  integer :: jsd !< The start j-index of cell centers within the data domain
  integer :: jed !< The end j-index of cell centers within the data domain
  real :: time_int              !< The time interval for any fields
                                !! that are offered for averaging [s].
  type(time_type) :: time_end   !< The end time of the valid interval for any offered field.
  logical :: ave_enabled = .false. !< True if averaging is enabled.

  !>@{ The following are 3D and 2D axis groups defined for output.  The names indicate
  !! the horizontal locations (B, T, Cu, or Cv) and vertical locations (here just 1).
  type(axes_grp) :: axesB1, axesT1, axesCu1, axesCv1
  !>@}
  type(axes_grp) :: axesNull !< An axis group for scalars

  ! Mask arrays for 2D diagnostics
  real, dimension(:,:),   pointer :: mask2dT   => null() !< 2D mask array for cell-center points [nondim]
  real, dimension(:,:),   pointer :: mask2dBu  => null() !< 2D mask array for cell-corner points [nondim]
  real, dimension(:,:),   pointer :: mask2dCu  => null() !< 2D mask array for east-face points [nondim]
  real, dimension(:,:),   pointer :: mask2dCv  => null() !< 2D mask array for north-face points [nondim]
  real, dimension(:,:),   pointer :: mask2dT_comp => null() !< 2D cell-center mask on the computational domain [nondim]

! Space for diagnostics is dynamically allocated as it is needed.
! The chunk size is how much the array should grow on each new allocation.
#define DIAG_ALLOC_CHUNK_SIZE 15
  type(diag_type), dimension(:), allocatable :: diags !< The list of diagnostics
  integer :: next_free_diag_id !< The next unused diagnostic ID

  !> default missing value to be sent to ALL diagnostics registrations [various]
  real :: missing_value = -1.0e34

  type(ocean_grid_type), pointer :: G => null()  !< The ocean grid type
  type(unit_scale_type), pointer :: US => null() !< A dimensional unit scaling type

  !> Number of checksum-only diagnostics
  integer :: num_chksum_diags

end type diag_ctrl

!>@{ CPU clocks
integer :: id_clock_diag_mediator
!>@}

contains

!> Set up the grid and axis information for use by the ice shelf model.
subroutine set_IS_axes_info(G, diag_cs, axes_set_name)
  type(ocean_grid_type), intent(in)    :: G   !< The horizontal grid type
  type(diag_ctrl),       intent(inout) :: diag_cs !< A structure that is used to regulate diagnostic output
  character(len=*), optional, intent(in) :: axes_set_name !<  A name to use for this set of axes.
                                                !! The default is "ice".
!   This subroutine sets up the grid and axis information for use by the ice shelf model.

  ! Local variables
  integer :: id_xq, id_yq, id_xh, id_yh, id_null
  integer :: i, j
  character(len=80) :: set_name
  real, allocatable, dimension(:) :: IaxB, iax ! Index-based integer and half-integer i-axis labels [nondim]
  real, allocatable, dimension(:) :: JaxB, jax ! Index-based integer and half-integer j-axis labels [nondim]

  set_name = "ice_shelf" ; if (present(axes_set_name)) set_name = trim(axes_set_name)

  if (diag_cs%index_space_axes) then
    allocate(IaxB(G%IsgB:G%IegB))
    do I=G%IsgB,G%IegB
      Iaxb(I) = real(I)
    enddo
    allocate(iax(G%isg:G%ieg))
    do i=G%isg,G%ieg
      iax(i) = real(i)-0.5
    enddo
    allocate(JaxB(G%JsgB:G%JegB))
    do J=G%JsgB,G%JegB
      JaxB(J) = real(J)
    enddo
    allocate(jax(G%jsg:G%jeg))
    do j=G%jsg,G%jeg
      jax(j) = real(j)-0.5
    enddo
  endif

  ! Horizontal axes for the native grids.
  if (diag_cs%index_space_axes) then
    if (G%symmetric) then
      id_xq = MOM_diag_axis_init('Iq', IaxB(G%IsgB:G%IegB), 'none', 'x', &
          'Boundary (q) point grid-space longitude', G%Domain, position=EAST, set_name=set_name)
      id_yq = MOM_diag_axis_init('Jq', JaxB(G%JsgB:G%JegB), 'none', 'y', &
          'Boundary (q) point grid-space latitude', G%Domain, position=NORTH, set_name=set_name)
    else
      id_xq = MOM_diag_axis_init('Iq', IaxB(G%isg:G%ieg), 'none', 'x', &
          'Boundary (q) point grid-space longitude', G%Domain, position=EAST, set_name=set_name)
      id_yq = MOM_diag_axis_init('Jq', JaxB(G%jsg:G%jeg), 'none', 'y', &
          'Boundary (q) point grid-space latitude', G%Domain, position=NORTH, set_name=set_name)
    endif

    id_xh = MOM_diag_axis_init('ih', iax, 'none', 'x', &
        'Tracer (h) point grid-space longitude', G%Domain, set_name=set_name)
    id_yh = MOM_diag_axis_init('jh', jax, 'none', 'y', &
        'Tracer (h) point grid-space latitude', G%Domain, set_name=set_name)
  else
    if (G%symmetric) then
      id_xq = MOM_diag_axis_init('xB', G%gridLonB(G%isgB:G%iegB), G%x_axis_units, 'x', &
          'Boundary point nominal longitude', G%Domain, position=EAST, set_name=set_name)
      id_yq = MOM_diag_axis_init('yB', G%gridLatB(G%jsgB:G%jegB), G%y_axis_units, 'y', &
          'Boundary point nominal latitude', G%Domain, position=NORTH, set_name=set_name)

    else
      id_xq = MOM_diag_axis_init('xB', G%gridLonB(G%isg:G%ieg), G%x_axis_units, 'x', &
          'Boundary point nominal longitude', G%Domain, position=EAST, set_name=set_name)
      id_yq = MOM_diag_axis_init('yB', G%gridLatB(G%jsg:G%jeg), G%y_axis_units, 'y', &
          'Boundary point nominal latitude', G%Domain, position=NORTH, set_name=set_name)

    endif
    id_xh = MOM_diag_axis_init('xT', G%gridLonT(G%isg:G%ieg), G%x_axis_units, 'x', &
        'Tracer point nominal longitude', G%Domain, set_name=set_name)
    id_yh = MOM_diag_axis_init('yT', G%gridLatT(G%jsg:G%jeg), G%y_axis_units, 'y', &
        'Tracer point nominal latitude', G%Domain, set_name=set_name)
  endif

  ! Axis groupings for 2-D arrays
  call define_axes_group(diag_cs, (/id_xh, id_yh/), diag_cs%axesT1, &
       x_cell_method='mean', y_cell_method='mean', is_h_point=.true.)
  call define_axes_group(diag_cs, (/id_xq, id_yq/), diag_cs%axesB1, &
       x_cell_method='point', y_cell_method='point', is_q_point=.true.)
  call define_axes_group(diag_cs, (/id_xq, id_yh/), diag_cs%axesCu1, &
       x_cell_method='point', y_cell_method='mean', is_u_point=.true.)
  call define_axes_group(diag_cs, (/id_xh, id_yq/), diag_cs%axesCv1, &
       x_cell_method='mean', y_cell_method='point', is_v_point=.true.)

  ! Axis group for special null axis for scalars from diag manager.
  id_null = MOM_diag_axis_init('scalar_axis', (/0./), 'none', 'N', 'none', null_axis=.true.)
  call define_axes_group(diag_cs, (/ id_null /), diag_cs%axesNull)

  if (diag_cs%index_space_axes) then
    deallocate(IaxB, iax, JaxB, jax)
  endif

end subroutine set_IS_axes_info

!> Attaches the id of cell areas to axes groups for use with cell_measures
subroutine diag_register_area_ids(diag_cs, id_area_t, id_area_q)
  type(diag_ctrl), intent(inout) :: diag_cs   !< Diagnostics control structure
  integer,   optional, intent(in)    :: id_area_t !< Diag_mediator id for area of h-cells
  integer,   optional, intent(in)    :: id_area_q !< Diag_mediator id for area of q-cells
  ! Local variables
  integer :: fms_id, i
  if (present(id_area_t)) then
    fms_id = diag_cs%diags(id_area_t)%fms_diag_id
    diag_cs%axesT1%id_area = fms_id
  endif
  if (present(id_area_q)) then
    fms_id = diag_cs%diags(id_area_q)%fms_diag_id
    diag_cs%axesB1%id_area = fms_id
  endif
end subroutine diag_register_area_ids

!> Define a group of "axes" from a list of handles and associate a mask with it
subroutine define_axes_group(diag_cs, handles, axes, &
                             x_cell_method, y_cell_method, &
                             is_h_point, is_q_point, is_u_point, is_v_point)
  type(diag_ctrl), target,    intent(in)  :: diag_cs !< Structure used to regulate diagnostic output
  integer, dimension(:),      intent(in)  :: handles !< A list of 1D axis handles that define the axis group
  type(axes_grp),             intent(out) :: axes    !< The group of axes that is set up here
  character(len=*), optional, intent(in)  :: x_cell_method !< A x-direction cell method used to construct the
                                                           !! "cell_methods" attribute in CF convention
  character(len=*), optional, intent(in)  :: y_cell_method !< A y-direction cell method used to construct the
                                                           !! "cell_methods" attribute in CF convention
  logical,          optional, intent(in)  :: is_h_point !< If true, indicates this axes group for h-point
                                                        !! located fields
  logical,          optional, intent(in)  :: is_q_point !< If true, indicates this axes group for q-point
                                                        !! located fields
  logical,          optional, intent(in)  :: is_u_point !< If true, indicates this axes group for
                                                        !! u-point located fields
  logical,          optional, intent(in)  :: is_v_point !< If true, indicates this axes group for
                                                        !! v-point located fields

  ! Local variables
  integer :: n

  n = size(handles)
  if (n<1 .or. n>2) call MOM_error(FATAL, "define_axes_group: wrong size for list of handles!")
  allocate( axes%handles(n) )
  axes%id = ints_to_string(handles, max(n,2)) ! Identifying string
  axes%rank = n
  axes%handles(:) = handles(:)
  axes%diag_cs => diag_cs ! A (circular) link back to the diag_ctrl structure

  if ((axes%rank<2) .and. (present(x_cell_method) .or. present(x_cell_method))) &
    call MOM_error(FATAL, 'define_axes_group: Can not set x_cell_method or y_cell_method for rank<2.')
  axes%x_cell_method = '' ; if (present(x_cell_method)) axes%x_cell_method = trim(x_cell_method)
  axes%y_cell_method = '' ; if (present(y_cell_method)) axes%y_cell_method = trim(y_cell_method)

  if (present(is_h_point)) axes%is_h_point = is_h_point
  if (present(is_q_point)) axes%is_q_point = is_q_point
  if (present(is_u_point)) axes%is_u_point = is_u_point
  if (present(is_v_point)) axes%is_v_point = is_v_point

  ! Setup masks for this axes group
  axes%mask2d => null()
  if (axes%rank==2) then
    if (axes%is_h_point) axes%mask2d => diag_cs%mask2dT
    if (axes%is_h_point) axes%mask2d_comp => diag_cs%mask2dT_comp
    if (axes%is_u_point) axes%mask2d => diag_cs%mask2dCu
    if (axes%is_v_point) axes%mask2d => diag_cs%mask2dCv
    if (axes%is_q_point) axes%mask2d => diag_cs%mask2dBu
  endif

end subroutine define_axes_group

!> Set up the array extents for doing diagnostics
subroutine set_IS_diag_mediator_grid(G, diag_cs)
  type(ocean_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(diag_ctrl),     intent(inout) :: diag_cs !< Structure used to regulate diagnostic output

  diag_cs%is = G%isc - (G%isd-1) ; diag_cs%ie = G%iec - (G%isd-1)
  diag_cs%js = G%jsc - (G%jsd-1) ; diag_cs%je = G%jec - (G%jsd-1)
  diag_cs%isd = G%isd ; diag_cs%ied = G%ied
  diag_cs%jsd = G%jsd ; diag_cs%jed = G%jed

end subroutine set_IS_diag_mediator_grid

!> Make a real ice shelf scalar diagnostic available for averaging or output
subroutine post_IS_data_0d(diag_field_id, field, diag_cs, is_static)
  integer,           intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                 !! previous call to register_MOM_IS_diag_field.
  real,              intent(in) :: field         !< real value being offered for output or averaging
                                                 !! in internally scaled arbitrary units [A ~> a]
  type(diag_ctrl), target, intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.

  ! Local variables
  real :: locfield ! The field being offered in arbitrary unscaled units [a]
  logical :: used, is_stat
  type(diag_type), pointer :: diag => null()

  integer :: time_days
  integer :: time_seconds
  character(len=300) :: debug_mesg

  if (id_clock_diag_mediator>0) call cpu_clock_begin(id_clock_diag_mediator)
  is_stat = .false. ; if (present(is_static)) is_stat = is_static

  ! Iterate over list of diag 'variants', e.g. CMOR aliases, call send_data
  ! for each one.
  call assert(diag_field_id < diag_cs%next_free_diag_id, &
              'post_IS_data_0d: Unregistered diagnostic id')
  diag => diag_cs%diags(diag_field_id)

  do while (associated(diag))
    locfield = field
    if (diag%conversion_factor /= 0.) &
      locfield = locfield * diag%conversion_factor

    if (diag_cs%diag_as_chksum) then
      ! Append timestep to mesg
      call get_time(diag_cs%time_end, time_seconds, days=time_days)
      write(debug_mesg, '(a, 1x, i0, 1x, i0)') &
          trim(diag%debug_str), time_days, time_seconds

      call chksum0(locfield, debug_mesg, logunit=diag_cs%chksum_iounit)
    elseif (is_stat) then
      used = send_data_infra(diag%fms_diag_id, locfield)
    elseif (diag_cs%ave_enabled) then
      used = send_data_infra(diag%fms_diag_id, locfield, diag_cs%time_end)
    endif

    diag => diag%next
  enddo

  if (id_clock_diag_mediator>0) call cpu_clock_end(id_clock_diag_mediator)
end subroutine post_IS_data_0d


!> Make a real 2-d array diagnostic available for averaging or output
subroutine post_IS_data_2d(diag_field_id, field, diag_cs, is_static, mask)
  integer,           intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                 !! previous call to register_MOM_IS_diag_field.
  real,      target, intent(in) :: field(:,:)    !< 2-d array being offered for output or averaging
                                                 !! in internally scaled arbitrary units [A ~> a]
  type(diag_ctrl), target, intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.
  real,    optional, intent(in) :: mask(:,:) !< If present, use this real array as the data mask [nondim]

  ! Local variables
  type(diag_type), pointer :: diag => NULL()

  if (id_clock_diag_mediator>0) call cpu_clock_begin(id_clock_diag_mediator)

  ! Iterate over list of diag 'variants' (e.g. CMOR aliases) and post each.
  call assert(diag_field_id < diag_cs%next_free_diag_id, &
              'post_IS_data_2d: Unregistered diagnostic id')
  diag => diag_cs%diags(diag_field_id)
  do while (associated(diag))
    call post_data_2d_low(diag, field, diag_cs, is_static, mask)
    diag => diag%next
  enddo

  if (id_clock_diag_mediator>0) call cpu_clock_end(id_clock_diag_mediator)
end subroutine post_IS_data_2d

!> Make a real 2-d array diagnostic available for averaging or output
!! using a diag_type instead of an integer id.
subroutine post_data_2d_low(diag, field, diag_cs, is_static, mask)
  type(diag_type),   intent(in) :: diag       !< A structure describing the diagnostic to post
  real,    target,   intent(in) :: field(:,:) !< 2-d array being offered for output or averaging
                                              !! in internally scaled arbitrary units [A ~> a]
  type(diag_ctrl),   intent(in) :: diag_CS   !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.
  real, optional, target, intent(in) :: mask(:,:) !< If present, use this real array as the data mask [nondim]

  ! Local variables
  real, dimension(:,:), pointer :: locfield ! The field being offered in arbitrary unscaled units [a]
  real, dimension(:,:), pointer :: locmask  ! A pointer to the data mask to use [nondim]
  logical :: used  ! The return value of send_data is not used for anything.
  logical :: is_stat
  logical :: i_data, j_data ! True if the field is on the data domain in the i or j directions.
  integer :: cszi, cszj, dszi, dszj
  integer :: isv, iev, jsv, jev, i, j
  integer :: time_days, time_seconds
  character(len=300) :: mesg
  character(len=300) :: debug_mesg

  locfield => NULL()
  locmask => NULL()
  is_stat = .false. ; if (present(is_static)) is_stat = is_static

  ! Determine the proper array indices, noting that because of the (:,:)
  ! declaration of field, symmetric arrays are using a SW-grid indexing,
  ! but non-symmetric arrays are using a NE-grid indexing.  Send_data
  ! actually only uses the difference between ie and is to determine
  ! the output data size and assumes that halos are symmetric.
  isv = diag_cs%is ; iev = diag_cs%ie ; jsv = diag_cs%js ; jev = diag_cs%je

  cszi = (diag_cs%ie-diag_cs%is) +1 ; dszi = (diag_cs%ied-diag_cs%isd) +1
  cszj = (diag_cs%je-diag_cs%js) +1 ; dszj = (diag_cs%jed-diag_cs%jsd) +1
  if ( size(field,1) == dszi ) then
    isv = diag_cs%is ; iev = diag_cs%ie ; i_data = .true.   ! Data domain
  elseif ( size(field,1) == dszi + 1 ) then
    isv = diag_cs%is ; iev = diag_cs%ie+1 ; i_data = .true. ! Symmetric data domain
  elseif ( size(field,1) == cszi ) then
    isv = 1 ; iev = cszi ; i_data = .false. ! Computational domain
  elseif ( size(field,1) == cszi + 1 ) then
    isv = 1 ; iev = cszi+1 ; i_data = .false. ! Symmetric computational domain
  else
    write (mesg,*) " peculiar size ",size(field,1)," in i-direction\n"//&
       "does not match one of ", cszi, cszi+1, dszi, dszi+1
    call MOM_error(FATAL,"post_IS_data_2d_low: "//trim(diag%debug_str)//trim(mesg))
  endif

  if ( size(field,2) == dszj ) then
    jsv = diag_cs%js ; jev = diag_cs%je ; j_data = .true.   ! Data domain
  elseif ( size(field,2) == dszj + 1 ) then
    jsv = diag_cs%js ; jev = diag_cs%je+1 ; j_data = .true. ! Symmetric data domain
  elseif ( size(field,2) == cszj ) then
    jsv = 1 ; jev = cszj ; j_data = .false. ! Computational domain
  ! This was: elseif ( size(field,1) == cszj + 1 ) then
  elseif ( size(field,2) == cszj + 1 ) then
    jsv = 1 ; jev = cszj+1 ; j_data = .false. ! Symmetric computational domain
  else
    write (mesg,*) " peculiar size ",size(field,2)," in j-direction\n"//&
       "does not match one of ", cszj, cszj+1, dszj, dszj+1
    call MOM_error(FATAL,"post_IS_data_2d_low: "//trim(diag%debug_str)//trim(mesg))
  endif

  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) then
    allocate( locfield( lbound(field,1):ubound(field,1), lbound(field,2):ubound(field,2) ) )
    do j=jsv,jev ; do i=isv,iev
      if (field(i,j) == diag_cs%missing_value) then
        locfield(i,j) = diag_cs%missing_value
      else
        locfield(i,j) = field(i,j) * diag%conversion_factor
      endif
    enddo ; enddo
    locfield(isv:iev,jsv:jev) = field(isv:iev,jsv:jev) * diag%conversion_factor
  else
    locfield => field
  endif

  ! Handle cases where the data and computational domain are the same size.
  if (diag_cs%ied-diag_cs%isd == diag_cs%ie-diag_cs%is) i_data = j_data
  if (diag_cs%jed-diag_cs%jsd == diag_cs%je-diag_cs%js) j_data = i_data
  if ( i_data .NEQV. j_data ) then
    call MOM_error(FATAL, "post_IS_data_2d: post_IS_data called for "//&
                   trim(diag%debug_str)//" with mixed computational and data domain array sizes.")
  endif

  if (present(mask)) then
    locmask => mask
  elseif (.not.is_stat) then  ! Static fields do not have assigned axes.
    if (i_data .and. associated(diag%axes%mask2d)) then
      locmask => diag%axes%mask2d
    elseif ((.not.i_data) .and. associated(diag%axes%mask2d_comp)) then
      locmask => diag%axes%mask2d_comp
    endif
  endif
  if (associated(locmask)) call assert(size(locfield) == size(locmask), &
        'post_data_2d_low: mask size mismatch: '//trim(diag%debug_str))

  if (diag_cs%diag_as_chksum) then
    ! Append timestep to mesg
    call get_time(diag_cs%time_end, time_seconds, days=time_days)
    write(debug_mesg, '(a, 1x, i0, 1x, i0)') &
        trim(diag%debug_str), time_days, time_seconds

    if (diag%axes%is_h_point) then
      call hchksum(locfield, debug_mesg, diag_cs%G%HI, &
                   logunit=diag_cs%chksum_iounit)
    elseif (diag%axes%is_u_point) then
      call uchksum(locfield, debug_mesg, diag_cs%G%HI, &
                   logunit=diag_cs%chksum_iounit)
    elseif (diag%axes%is_v_point) then
      call vchksum(locfield, debug_mesg, diag_cs%G%HI, &
                   logunit=diag_cs%chksum_iounit)
    elseif (diag%axes%is_q_point) then
      call Bchksum(locfield, debug_mesg, diag_cs%G%HI, &
                   logunit=diag_cs%chksum_iounit)
    else
      call MOM_error(FATAL, "post_data_2d_low: unknown axis type.")
    endif
  else
    if (is_stat) then
      if (associated(locmask)) then
        used = send_data_infra(diag%fms_diag_id, locfield, &
                         is_in=isv, ie_in=iev, js_in=jsv, je_in=jev, rmask=locmask)
      else
        used = send_data_infra(diag%fms_diag_id, locfield, &
                         is_in=isv, ie_in=iev, js_in=jsv, je_in=jev)
      endif
    elseif (diag_cs%ave_enabled) then
      if (associated(locmask)) then
        used = send_data_infra(diag%fms_diag_id, locfield, &
                         is_in=isv, ie_in=iev, js_in=jsv, je_in=jev, &
                         time=diag_cs%time_end, weight=diag_cs%time_int, rmask=locmask)
      else
        used = send_data_infra(diag%fms_diag_id, locfield, &
                         is_in=isv, ie_in=iev, js_in=jsv, je_in=jev, &
                         time=diag_cs%time_end, weight=diag_cs%time_int)
      endif
    endif
  endif

  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) deallocate( locfield )

end subroutine post_data_2d_low

!> Enable the accumulation of time averages over the specified time interval.
subroutine enable_averaging(time_int_in, time_end_in, diag_cs)
  real,            intent(in)    :: time_int_in !< The time interval [s] over which any
                                                !! values that are offered are valid.
  type(time_type), intent(in)    :: time_end_in !< The end time of the valid interval
  type(diag_ctrl), intent(inout) :: diag_cs     !< Structure used to regulate diagnostic output
  ! This subroutine enables the accumulation of time averages over the specified time interval.

!  if (num_file==0) return
  diag_cs%time_int = time_int_in
  diag_cs%time_end = time_end_in
  diag_cs%ave_enabled = .true.
end subroutine enable_averaging

!> Enable the accumulation of time averages over the specified time interval in time units.
subroutine enable_averages(time_int, time_end, diag_CS, T_to_s)
  real,            intent(in)    :: time_int !< The time interval over which any values
                                             !! that are offered are valid [T ~> s].
  type(time_type), intent(in)    :: time_end !< The end time of the valid interval.
  type(diag_ctrl), intent(inout) :: diag_CS  !< A structure that is used to regulate diagnostic output
  real,  optional, intent(in)    :: T_to_s   !< A conversion factor for time_int to seconds [s T-1 ~> 1].
  ! This subroutine enables the accumulation of time averages over the specified time interval.

  if (present(T_to_s)) then
    diag_cs%time_int = time_int*T_to_s
  elseif (associated(diag_CS%US)) then
    diag_cs%time_int = time_int*diag_CS%US%T_to_s
  else
    diag_cs%time_int = time_int
  endif
  diag_cs%time_end = time_end
  diag_cs%ave_enabled = .true.
end subroutine enable_averages

!> Call this subroutine to avoid averaging any offered fields.
subroutine disable_averaging(diag_cs)
  type(diag_ctrl), intent(inout) :: diag_cs !< Structure used to regulate diagnostic output

  diag_cs%time_int = 0.0
  diag_cs%ave_enabled = .false.
end subroutine disable_averaging

!> Indicate whether averaging diagnostics is currently enabled
logical function query_averaging_enabled(diag_cs, time_int, time_end)
  type(diag_ctrl),           intent(in)  :: diag_cs  !< Structure used to regulate diagnostic output
  real,            optional, intent(out) :: time_int !< Current setting of diag_cs%time_int [s]
  type(time_type), optional, intent(out) :: time_end !< Current setting of diag_cs%time_end

  if (present(time_int)) time_int = diag_cs%time_int
  if (present(time_end)) time_end = diag_cs%time_end
  query_averaging_enabled = diag_cs%ave_enabled
end function query_averaging_enabled

!> This subroutine initializes the diag_manager via the MOM6 infrastructure
subroutine MOM_IS_diag_mediator_infrastructure_init(err_msg)
  character(len=*), optional, intent(out)   :: err_msg !< An error message

  call MOM_diag_manager_init(err_msg=err_msg)
end subroutine MOM_IS_diag_mediator_infrastructure_init

!> This function returns the valid end time for use with diagnostics that are
!! handled outside of the MOM6 diagnostics infrastructure.
function get_diag_time_end(diag_cs)
  type(diag_ctrl), intent(in)  :: diag_cs !< Structure used to regulate diagnostic output
  type(time_type) :: get_diag_time_end
  !   This function returns the valid end time for diagnostics that are handled
  ! outside of the MOM6 infrastructure, such as via the generic tracer code.

  get_diag_time_end = diag_cs%time_end
end function get_diag_time_end

!> Returns the "diag_mediator" handle for a group (native, CMOR, ...) of diagnostics
!! derived from one field.
function register_MOM_IS_diag_field(module_name, field_name, axes_in, init_time, &
            long_name, units, missing_value, range, mask_variant, standard_name, &
            verbose, do_not_log, err_msg, interp_method, tile_count, cmor_field_name, &
            cmor_long_name, cmor_units, cmor_standard_name, cell_methods, &
            x_cell_method, y_cell_method, conversion) result (register_diag_field)
  integer :: register_diag_field  !< The returned diagnostic handle
  character(len=*),           intent(in) :: module_name !< Name of this module, usually "ice_model"
  character(len=*),           intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp),     target, intent(in) :: axes_in   !< Container with up to 3 integer handles that
                                                      !! indicates axes for this field
  type(time_type),            intent(in) :: init_time !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values in
                                                          !! output files, in unscaled arbitrary units [a]
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
                                                     !! in arbitrary units [a]
  logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided with
                                                         !! post_IS_data calls (not used in MOM?)
  logical,          optional, intent(in) :: verbose !< If true, FMS is verbose (not used in MOM?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be
                                                    !! placed (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not
                                                          !! be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count    !< no clue (not used in MOM?)
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
  character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
  character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
  character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
  character(len=*), optional, intent(in) :: cell_methods !< String to append as cell_methods attribute. Use '' to
                                                         !! have no attribute.  If present, this overrides the
                                                         !! default constructed from the default for
                                                         !! each individual axis direction.
  character(len=*), optional, intent(in) :: x_cell_method !< Specifies the cell method for the x-direction.
                                                         !! Use '' have no method.
  character(len=*), optional, intent(in) :: y_cell_method !< Specifies the cell method for the y-direction.
                                                         !! Use '' have no method.
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to files,
                                                       !! often including factors to undo internal scaling and
                                                       !! in units of [a A-1 ~> 1]

  ! Local variables
  real :: MOM_missing_value ! A value used to indicate missing values in output files, in arbitrary units [a]
  type(diag_ctrl), pointer :: diag_cs => NULL() ! A structure that is used to regulate diagnostic output
  type(axes_grp), pointer :: axes
  integer :: dm_id
  character(len=256) :: msg
  character(len=256) :: cm_string ! A string describing the cell methods returned from attach_cell_methods.
  character(len=256) :: new_module_name
  character(len=480) :: module_list, var_list
  character(len=24)  :: dimensions
  integer :: num_modnm, num_varnm
  logical :: active

  diag_cs => axes_in%diag_cs

  ! Check if the axes match a standard grid axis.
  ! If not, allocate the new axis and copy the contents.
  if (axes_in%id == diag_cs%axesT1%id) then
    axes => diag_cs%axesT1
  elseif (axes_in%id == diag_cs%axesB1%id) then
    axes => diag_cs%axesB1
  elseif (axes_in%id == diag_cs%axesCu1%id) then
    axes => diag_cs%axesCu1
  elseif (axes_in%id == diag_cs%axesCv1%id) then
    axes => diag_cs%axesCv1
  else
    allocate(axes)
    axes = axes_in
  endif

  MOM_missing_value = axes%diag_cs%missing_value
  if (present(missing_value)) MOM_missing_value = missing_value

  diag_cs => axes%diag_cs
  dm_id = -1

  module_list = "{"//trim(module_name)
  num_modnm = 1

  ! Register the native diagnostic
  active = register_diag_field_expand_cmor(dm_id, module_name, field_name, axes, &
             init_time, long_name=long_name, units=units, missing_value=MOM_missing_value, &
             range=range, mask_variant=mask_variant, standard_name=standard_name, &
             verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
             interp_method=interp_method, tile_count=tile_count, &
             cmor_field_name=cmor_field_name, cmor_long_name=cmor_long_name, &
             cmor_units=cmor_units, cmor_standard_name=cmor_standard_name, &
             cell_methods=cell_methods, x_cell_method=x_cell_method, y_cell_method=y_cell_method, &
             conversion=conversion)
  num_varnm = 1 ; var_list = "{"//trim(field_name)
  if (present(cmor_field_name)) then
    num_varnm = num_varnm + 1
    var_list = trim(var_list)//","//trim(cmor_field_name)
  endif
  var_list = trim(var_list)//"}"

  dimensions = ""
  if (axes_in%is_h_point)   dimensions = trim(dimensions)//" xh, yh,"
  if (axes_in%is_q_point)   dimensions = trim(dimensions)//" xq, yq,"
  if (axes_in%is_u_point)   dimensions = trim(dimensions)//" xq, yh,"
  if (axes_in%is_v_point)   dimensions = trim(dimensions)//" xh, yq,"
  if (len_trim(dimensions) > 0) dimensions = trim_trailing_commas(dimensions)

  if (is_root_pe() .and. (diag_CS%available_diag_doc_unit > 0)) then
    msg = ''
    if (present(cmor_field_name)) msg = 'CMOR equivalent is "'//trim(cmor_field_name)//'"'
    call attach_cell_methods(-1, axes, cm_string, cell_methods, x_cell_method, y_cell_method)
    module_list = trim(module_list)//"}"
    if (num_modnm <= 1) module_list = module_name
    if (num_varnm <= 1) var_list = ''

    call log_available_diag(dm_id>0, module_list, field_name, cm_string, msg, diag_CS, &
                            long_name, units, standard_name, variants=var_list, dimensions=dimensions)
  endif

  register_diag_field = dm_id

end function register_MOM_IS_diag_field

!> Returns True if either the native or CMOR version of the diagnostic were registered. Updates 'dm_id'
!! after calling register_diag_field_expand_axes() for both native and CMOR variants of the field.
logical function register_diag_field_expand_cmor(dm_id, module_name, field_name, axes, init_time, &
            long_name, units, missing_value, range, mask_variant, standard_name,      &
            verbose, do_not_log, err_msg, interp_method, tile_count, cmor_field_name, &
            cmor_long_name, cmor_units, cmor_standard_name, cell_methods, &
            x_cell_method, y_cell_method, conversion)
  integer,          intent(inout) :: dm_id !< The diag_mediator ID for this diagnostic group
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ice_model" or "ice_model_fast"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp),   intent(in) :: axes !< Container with up to 3 integer handles that indicates axes
                                             !! for this field
  type(time_type),  intent(in) :: init_time !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values in
                                                          !! output files, in unscaled arbitrary units [a]
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
                                                     !! in arbitrary units [a]
  logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided
                                                         !! with post_data calls (not used in MOM?)
  logical,          optional, intent(in) :: verbose !< If true, FMS is verbose (not used in MOM?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be
                                                         !! placed (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should
                                                         !! not be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count !< no clue (not used in MOM?)
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
  character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
  character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
  character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
  character(len=*), optional, intent(in) :: cell_methods !< String to append as cell_methods attribute.
                                                         !! Use '' to have no attribute. If present, this
                                                         !! overrides the default constructed from the default
                                                         !! for each individual axis direction.
  character(len=*), optional, intent(in) :: x_cell_method !< Specifies the cell method for the x-direction.
                                                         !! Use '' have no method.
  character(len=*), optional, intent(in) :: y_cell_method !< Specifies the cell method for the y-direction.
                                                         !! Use '' have no method.
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to files,
                                                       !! often including factors to undo internal scaling and
                                                       !! in units of [a A-1 ~> 1]
  ! Local variables
  real :: MOM_missing_value ! A value used to indicate missing values in output files, in arbitrary units [a]
  type(diag_ctrl), pointer :: diag_cs => null()
  type(diag_type), pointer :: this_diag => null()
  integer :: fms_id
  character(len=256) :: posted_cmor_units, posted_cmor_standard_name, posted_cmor_long_name
  character(len=256) :: cm_string ! A string describing the cell methods returned from attach_cell_methods.

  MOM_missing_value = axes%diag_cs%missing_value
  if (present(missing_value)) MOM_missing_value = missing_value

  register_diag_field_expand_cmor = .false.
  diag_cs => axes%diag_cs

  ! Set up the 'primary' diagnostic, first get an underlying FMS id
  fms_id = register_diag_field_expand_axes(module_name, field_name, axes, init_time, &
             long_name=long_name, units=units, missing_value=MOM_missing_value, &
             range=range, mask_variant=mask_variant, standard_name=standard_name, &
             verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
             interp_method=interp_method, tile_count=tile_count)
  if (.not. diag_cs%diag_as_chksum) &
    call attach_cell_methods(fms_id, axes, cm_string, cell_methods, x_cell_method, y_cell_method)

  this_diag => null()
  if (fms_id /= DIAG_FIELD_NOT_FOUND) then
    call add_diag_to_list(diag_cs, dm_id, fms_id, this_diag, axes, module_name, field_name)
    if (present(conversion)) this_diag%conversion_factor = conversion
    register_diag_field_expand_cmor = .true.
  endif

  ! For the CMOR variation of the above diagnostic
  if (present(cmor_field_name) .and. .not. diag_cs%diag_as_chksum) then
    ! Fallback values for strings set to "NULL"
    posted_cmor_units = "not provided"         !
    posted_cmor_standard_name = "not provided" ! Values might be able to be replaced with a CS%missing field?
    posted_cmor_long_name = "not provided"     !

    ! If attributes are present for MOM variable names, use them first for the register_MOM_IS_diag_field
    ! call for CMOR verison of the variable
    if (present(units)) posted_cmor_units = units
    if (present(standard_name)) posted_cmor_standard_name = standard_name
    if (present(long_name)) posted_cmor_long_name = long_name

    ! If specified in the call to register_MOM_IS_diag_field, override attributes with the CMOR versions
    if (present(cmor_units)) posted_cmor_units = cmor_units
    if (present(cmor_standard_name)) posted_cmor_standard_name = cmor_standard_name
    if (present(cmor_long_name)) posted_cmor_long_name = cmor_long_name

    fms_id = register_diag_field_expand_axes(module_name, cmor_field_name, axes, init_time,    &
               long_name=trim(posted_cmor_long_name), units=trim(posted_cmor_units),                  &
               missing_value=MOM_missing_value, range=range, mask_variant=mask_variant,               &
               standard_name=trim(posted_cmor_standard_name), verbose=verbose, do_not_log=do_not_log, &
               err_msg=err_msg, interp_method=interp_method, tile_count=tile_count)
    call attach_cell_methods(fms_id, axes, cm_string, cell_methods, x_cell_method, y_cell_method)

    this_diag => null()
    if (fms_id /= DIAG_FIELD_NOT_FOUND) then
      call add_diag_to_list(diag_cs, dm_id, fms_id, this_diag, axes, module_name, field_name)
      if (present(conversion)) this_diag%conversion_factor = conversion
      register_diag_field_expand_cmor = .true.
    endif
  endif

end function register_diag_field_expand_cmor

!> Returns an FMS id from register_diag_field_fms (the diag_manager routine) after expanding axes
!! (axes-group) into handles and conditionally adding an FMS area_id for cell_measures.
integer function register_diag_field_expand_axes(module_name, field_name, axes, init_time, &
            long_name, units, missing_value, range, mask_variant, standard_name,  &
            verbose, do_not_log, err_msg, interp_method, tile_count)
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ocean_model"
                                              !! or "ice_shelf_model"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp), target, intent(in) :: axes !< Container with up to 3 integer handles that indicates
                                             !! axes for this field
  type(time_type),  intent(in) :: init_time !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values in
                                                          !! output files, in unscaled arbitrary units [a]
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
                                                     !! in arbitrary units [a]
  logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided
                                                         !! with post_data calls (not used in MOM?)
  logical,          optional, intent(in) :: verbose !< If true, FMS is verbose (not used in MOM?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something
                                                       !! (not used in MOM?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be
                                                         !! placed (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should
                                                         !! not be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count !< no clue (not used in MOM?)
  ! Local variables
  integer :: fms_id, area_id

  ! This gets the cell area associated with the grid location of this variable
  area_id = axes%id_area

  ! Get the FMS diagnostic id
  if (axes%diag_cs%diag_as_chksum) then
    fms_id = axes%diag_cs%num_chksum_diags + 1
    axes%diag_cs%num_chksum_diags = fms_id
  elseif (present(interp_method) .or. axes%is_h_point) then
    ! If interp_method is provided we must use it
    if (area_id>0) then
      fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
               init_time, long_name=long_name, units=units, missing_value=missing_value, &
               range=range, mask_variant=mask_variant, standard_name=standard_name, &
               verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
               interp_method=interp_method, tile_count=tile_count, area=area_id)
    else
      fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
               init_time, long_name=long_name, units=units, missing_value=missing_value, &
               range=range, mask_variant=mask_variant, standard_name=standard_name, &
               verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
               interp_method=interp_method, tile_count=tile_count)
    endif
  else
    ! If interp_method is not provided and the field is not at an h-point then interp_method='none'
    if (area_id>0) then
      fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
               init_time, long_name=long_name, units=units, missing_value=missing_value, &
               range=range, mask_variant=mask_variant, standard_name=standard_name, &
               verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
               interp_method='none', tile_count=tile_count, area=area_id)
    else
      fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
               init_time, long_name=long_name, units=units, missing_value=missing_value, &
               range=range, mask_variant=mask_variant, standard_name=standard_name, &
               verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
               interp_method='none', tile_count=tile_count)
    endif
  endif

  register_diag_field_expand_axes = fms_id

end function register_diag_field_expand_axes

!> Create a diagnostic type and attached to list
subroutine add_diag_to_list(diag_cs, dm_id, fms_id, this_diag, axes, module_name, field_name)
  type(diag_ctrl),        pointer       :: diag_cs !< Diagnostics mediator control structure
  integer,                intent(inout) :: dm_id !< The diag_mediator ID for this diagnostic group
  integer,                intent(in)    :: fms_id !< The FMS diag_manager ID for this diagnostic
  type(diag_type),        pointer       :: this_diag !< This diagnostic
  type(axes_grp), target, intent(in)    :: axes !< Container with up to 3 integer handles that
                                                !! indicates axes for this field
  character(len=*),       intent(in)    :: module_name !< Name of this module, usually
                                                       !! "ocean_model" or "ice_shelf_model"
  character(len=*),       intent(in)    :: field_name !< Name of diagnostic

  ! If the diagnostic is needed obtain a diag_mediator ID (if needed)
  if (dm_id == -1) dm_id = get_new_diag_id(diag_cs)
  ! Create a new diag_type to store links in
  call alloc_diag_with_id(dm_id, diag_cs, this_diag)
  call assert(associated(this_diag), 'add_diag_to_list: allocation failed for '//trim(field_name))
  ! Record FMS id, masks and conversion factor, in diag_type
  this_diag%fms_diag_id = fms_id
  this_diag%debug_str = trim(module_name)//"-"//trim(field_name)
  this_diag%axes => axes

end subroutine add_diag_to_list


!> Attaches "cell_methods" attribute to a variable based on defaults for axes_grp or optional arguments.
subroutine attach_cell_methods(id, axes, ostring, cell_methods, x_cell_method, y_cell_method)
  integer,                    intent(in)  :: id !< Handle to diagnostic
  type(axes_grp),             intent(in)  :: axes !< Container with up to 3 integer handles that indicates
                                                  !! axes for this field
  character(len=*),           intent(out) :: ostring !< The cell_methods strings that would appear in the file
  character(len=*), optional, intent(in)  :: cell_methods !< String to append as cell_methods attribute.
                                                         !! Use '' to have no attribute. If present, this
                                                         !! overrides the default constructed from the default
                                                         !! for each individual axis direction.
  character(len=*), optional, intent(in)  :: x_cell_method !< Specifies the cell method for the x-direction.
                                                         !! Use '' have no method.
  character(len=*), optional, intent(in)  :: y_cell_method !< Specifies the cell method for the y-direction.
                                                         !! Use '' have no method.
  ! Local variables
  character(len=9) :: axis_name
  logical :: x_mean, y_mean, x_sum, y_sum

  x_mean = .false.
  y_mean = .false.
  x_sum = .false.
  y_sum = .false.

  ostring = ''
  if (present(cell_methods)) then
    if (present(x_cell_method) .or. present(y_cell_method)) then
      call MOM_error(FATAL, "attach_cell_methods: " // &
           'Individual direction cell method was specified along with a "cell_methods" string.')
    endif
    if (len(trim(cell_methods))>0) then
      call MOM_diag_field_add_attribute(id, 'cell_methods', trim(cell_methods))
      ostring = trim(cell_methods)
    endif
  else
    if (present(x_cell_method)) then
      if (len(trim(x_cell_method))>0) then
        call get_MOM_diag_axis_name(axes%handles(1), axis_name)
        call MOM_diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(x_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(x_cell_method)
        if (trim(x_cell_method)=='mean') x_mean=.true.
        if (trim(x_cell_method)=='sum') x_sum=.true.
      endif
    else
      if (len(trim(axes%x_cell_method))>0) then
        call get_MOM_diag_axis_name(axes%handles(1), axis_name)
        call MOM_diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(axes%x_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(axes%x_cell_method)
        if (trim(axes%x_cell_method)=='mean') x_mean=.true.
        if (trim(axes%x_cell_method)=='sum') x_sum=.true.
      endif
    endif
    if (present(y_cell_method)) then
      if (len(trim(y_cell_method))>0) then
        call get_MOM_diag_axis_name(axes%handles(2), axis_name)
        call MOM_diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(y_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(y_cell_method)
        if (trim(y_cell_method)=='mean') y_mean=.true.
        if (trim(y_cell_method)=='sum') y_sum=.true.
      endif
    else
      if (len(trim(axes%y_cell_method))>0) then
        call get_MOM_diag_axis_name(axes%handles(2), axis_name)
        call MOM_diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(axes%y_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(axes%y_cell_method)
        if (trim(axes%y_cell_method)=='mean') y_mean=.true.
        if (trim(axes%y_cell_method)=='sum') y_sum=.true.
      endif
    endif
    if (x_mean .and. y_mean) then
      call MOM_diag_field_add_attribute(id, 'cell_methods', 'area:mean')
      ostring = trim(adjustl(ostring))//' area:mean'
    elseif (x_sum .and. y_sum) then
      call MOM_diag_field_add_attribute(id, 'cell_methods', 'area:sum')
      ostring = trim(adjustl(ostring))//' area:sum'
    endif
  endif
  ostring = adjustl(ostring)
end subroutine attach_cell_methods

!> Registers a non-array scalar diagnostic, returning an integer handle
function register_scalar_field_axes(module_name, field_name, axes, init_time, &
            long_name, units, missing_value, range, standard_name, &
            do_not_log, err_msg, interp_method, cmor_field_name, &
            cmor_long_name, cmor_units, cmor_standard_name, conversion) result (register_scalar_field)
  integer :: register_scalar_field !< An integer handle for a diagnostic array.
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ocean_model"
                                              !! or "ice_shelf_model"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp), target, intent(in) :: axes !< Container with up to 3 integer handles that
                                             !! indicates axes for this field
  type(time_type),  intent(in) :: init_time !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values in
                                                          !! output files, in unscaled arbitrary units [a]
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
                                                     !! in arbitrary units [a]
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be
                                                         !! placed (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not
                                                         !! be interpolated as a scalar
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
  character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
  character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
  character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to files,
                                                       !! often including factors to undo internal scaling and
                                                       !! in units of [a A-1 ~> 1]

  register_scalar_field = register_scalar_field_CS(module_name, field_name, init_time, axes%diag_cs, &
            long_name, units, missing_value, range, standard_name, &
            do_not_log, err_msg, interp_method, cmor_field_name, &
            cmor_long_name, cmor_units, cmor_standard_name, conversion)

end function register_scalar_field_axes

!> Registers a non-array scalar diagnostic, returning an integer handle
function register_scalar_field_CS(module_name, field_name, init_time, diag_cs, &
            long_name, units, missing_value, range, standard_name, &
            do_not_log, err_msg, interp_method, cmor_field_name, &
            cmor_long_name, cmor_units, cmor_standard_name, conversion) result (register_scalar_field)
  integer :: register_scalar_field !< An integer handle for a diagnostic array.
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ocean_model"
                                              !! or "ice_shelf_model"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(time_type),  intent(in) :: init_time !< Time at which a field is first available?
  type(diag_ctrl),  intent(inout) :: diag_CS !< Structure used to regulate diagnostic output
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values in
                                                          !! output files, in unscaled arbitrary units [a]
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
                                                     !! in arbitrary units [a]
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be
                                                         !! placed (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not
                                                         !! be interpolated as a scalar
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
  character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
  character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
  character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to files,
                                                       !! often including factors to undo internal scaling and
                                                       !! in units of [a A-1 ~> 1]

  ! Local variables
  real :: MOM_missing_value ! A value used to indicate missing values in output files, in arbitrary units [a]
  integer :: dm_id, fms_id
  type(diag_type), pointer :: diag => null(), cmor_diag => null()
  character(len=256) :: posted_cmor_units, posted_cmor_standard_name, posted_cmor_long_name
  character(len=16)  :: dimensions

  MOM_missing_value = diag_cs%missing_value
  if (present(missing_value)) MOM_missing_value = missing_value

  dm_id = -1
  diag => null()
  cmor_diag => null()

  if (diag_cs%diag_as_chksum) then
    fms_id = diag_cs%num_chksum_diags + 1
    diag_cs%num_chksum_diags = fms_id
  else
    fms_id = register_diag_field_infra(module_name, field_name, init_time, &
                long_name=long_name, units=units, missing_value=MOM_missing_value, &
                range=range, standard_name=standard_name, do_not_log=do_not_log, &
                err_msg=err_msg)
  endif

  if (fms_id /= DIAG_FIELD_NOT_FOUND) then
    dm_id = get_new_diag_id(diag_cs)
    call alloc_diag_with_id(dm_id, diag_cs, diag)
    call assert(associated(diag), 'register_scalar_field: diag allocation failed')
    diag%fms_diag_id = fms_id
    diag%debug_str = trim(module_name)//"-"//trim(field_name)
    if (present(conversion)) diag%conversion_factor = conversion
  endif

  if (present(cmor_field_name)) then
    ! Fallback values for strings set to "not provided"
    posted_cmor_units = "not provided"
    posted_cmor_standard_name = "not provided"
    posted_cmor_long_name = "not provided"

    ! If attributes are present for MOM variable names, use them as defaults for the
    ! register_diag_field_infra call for CMOR verison of the variable
    if (present(units)) posted_cmor_units = units
    if (present(standard_name)) posted_cmor_standard_name = standard_name
    if (present(long_name)) posted_cmor_long_name = long_name

    ! If specified in the call to register_MOM_IS_scalar_field, override attributes with the CMOR versions
    if (present(cmor_units)) posted_cmor_units = cmor_units
    if (present(cmor_standard_name)) posted_cmor_standard_name = cmor_standard_name
    if (present(cmor_long_name)) posted_cmor_long_name = cmor_long_name

    fms_id = register_diag_field_infra(module_name, cmor_field_name, init_time, &
           long_name=trim(posted_cmor_long_name), units=trim(posted_cmor_units), &
           missing_value=MOM_missing_value, range=range, &
           standard_name=trim(posted_cmor_standard_name), do_not_log=do_not_log, err_msg=err_msg)
    if (fms_id /= DIAG_FIELD_NOT_FOUND) then
      if (dm_id == -1) then
        dm_id = get_new_diag_id(diag_cs)
      endif
      call alloc_diag_with_id(dm_id, diag_cs, cmor_diag)
      cmor_diag%fms_diag_id = fms_id
      cmor_diag%debug_str = trim(module_name)//"-"//trim(cmor_field_name)
      if (present(conversion)) cmor_diag%conversion_factor = conversion
    endif
  endif

  dimensions = "scalar"

  ! Document diagnostics in list of available diagnostics
  if (is_root_pe() .and. diag_CS%available_diag_doc_unit > 0) then
    if (present(cmor_field_name)) then
      call log_available_diag(associated(diag), module_name, field_name, '', '', diag_CS, &
                              long_name, units, standard_name, &
                              variants="{"//trim(field_name)//","//trim(cmor_field_name)//"}", &
                              dimensions=dimensions)
    else
      call log_available_diag(associated(diag), module_name, field_name, '', '', diag_CS, &
                              long_name, units, standard_name, dimensions=dimensions)
    endif
  endif

  register_scalar_field = dm_id

end function register_scalar_field_CS

!> Registers a static diagnostic, returning an integer handle
function register_MOM_IS_static_field(module_name, field_name, axes, &
            long_name, units, missing_value, range, mask_variant, standard_name, &
            do_not_log, interp_method, tile_count, &
            cmor_field_name, cmor_long_name, cmor_units, cmor_standard_name, area, &
            x_cell_method, y_cell_method, area_cell_method, conversion) result(register_static_field)
  integer :: register_static_field !< An integer handle for a diagnostic array.
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ocean_model"
                                              !! or "ice_shelf_model"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp), target, intent(in) :: axes !< Container with up to 3 integer handles that
                                             !! indicates axes for this field
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values in
                                                          !! output files, in unscaled arbitrary units [a]
  real,             optional, intent(in) :: range(2) !< Valid range of a variable in arbitrary units [a]
  logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided with
                                                         !! post_IS_data calls (not used in MOM?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not
                                                         !! be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count   !< no clue (not used in MOM?)
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
  character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
  character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
  character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
  integer,          optional, intent(in) :: area !< fms_id for area_t
  character(len=*), optional, intent(in) :: x_cell_method !< Specifies the cell method for the x-direction.
  character(len=*), optional, intent(in) :: y_cell_method !< Specifies the cell method for the y-direction.
  character(len=*), optional, intent(in) :: area_cell_method !< Specifies the cell method for area
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to files,
                                                       !! often including factors to undo internal scaling and
                                                       !! in units of [a A-1 ~> 1]

  ! Local variables
  real :: MOM_missing_value ! A value used to indicate missing values in output files, in arbitrary units [a]
  type(diag_ctrl), pointer :: diag_cs => null() !< A structure that is used to regulate diagnostic output
  type(diag_type), pointer :: diag => null(), cmor_diag => null()
  integer :: dm_id, fms_id
  character(len=256) :: posted_cmor_units, posted_cmor_standard_name, posted_cmor_long_name
  character(len=9) :: axis_name
  character(len=24) :: dimensions

  MOM_missing_value = axes%diag_cs%missing_value
  if (present(missing_value)) MOM_missing_value = missing_value

  diag_cs => axes%diag_cs
  dm_id = -1
  diag => null()
  cmor_diag => null()

  if (diag_cs%diag_as_chksum) then
    fms_id = diag_cs%num_chksum_diags + 1
    diag_cs%num_chksum_diags = fms_id
  else
    fms_id = register_static_field_infra(module_name, field_name, axes%handles, &
           long_name=long_name, units=units, missing_value=MOM_missing_value, &
           range=range, mask_variant=mask_variant, standard_name=standard_name, &
           do_not_log=do_not_log, &
           interp_method=interp_method, tile_count=tile_count, area=area)
  endif

  if (fms_id /= DIAG_FIELD_NOT_FOUND) then
    dm_id = get_new_diag_id(diag_cs)
    call alloc_diag_with_id(dm_id, diag_cs, diag)
    call assert(associated(diag), 'register_static_field: diag allocation failed')
    diag%fms_diag_id = fms_id
    diag%debug_str = trim(module_name)//"-"//trim(field_name)
    if (present(conversion)) diag%conversion_factor = conversion

    if (diag_cs%diag_as_chksum) then
      diag%axes => axes
    else
      if (present(x_cell_method)) then
        call get_MOM_diag_axis_name(axes%handles(1), axis_name)
        call MOM_diag_field_add_attribute(fms_id, 'cell_methods', &
            trim(axis_name)//':'//trim(x_cell_method))
      endif
      if (present(y_cell_method)) then
        call get_MOM_diag_axis_name(axes%handles(2), axis_name)
        call MOM_diag_field_add_attribute(fms_id, 'cell_methods', &
            trim(axis_name)//':'//trim(y_cell_method))
      endif
      if (present(area_cell_method)) then
        call MOM_diag_field_add_attribute(fms_id, 'cell_methods', &
            'area:'//trim(area_cell_method))
      endif
    endif
  endif

  if (present(cmor_field_name) .and. .not. diag_cs%diag_as_chksum) then
    ! Fallback values for strings set to "not provided"
    posted_cmor_units = "not provided"
    posted_cmor_standard_name = "not provided"
    posted_cmor_long_name = "not provided"

    ! If attributes are present for MOM variable names, use them first for the register_static_field
    ! call for CMOR verison of the variable
    if (present(units)) posted_cmor_units = units
    if (present(standard_name)) posted_cmor_standard_name = standard_name
    if (present(long_name)) posted_cmor_long_name = long_name

    ! If specified in the call to register_static_field, override attributes with the CMOR versions
    if (present(cmor_units)) posted_cmor_units = cmor_units
    if (present(cmor_standard_name)) posted_cmor_standard_name = cmor_standard_name
    if (present(cmor_long_name)) posted_cmor_long_name = cmor_long_name

    fms_id = register_static_field_infra(module_name, cmor_field_name, axes%handles, &
                long_name=trim(posted_cmor_long_name), units=trim(posted_cmor_units), &
                missing_value=MOM_missing_value, range=range, mask_variant=mask_variant, &
                standard_name=trim(posted_cmor_standard_name), do_not_log=do_not_log, &
                interp_method=interp_method, tile_count=tile_count, area=area)
    if (fms_id /= DIAG_FIELD_NOT_FOUND) then
      if (dm_id == -1) then
        dm_id = get_new_diag_id(diag_cs)
      endif
      call alloc_diag_with_id(dm_id, diag_cs, cmor_diag)
      cmor_diag%fms_diag_id = fms_id
      cmor_diag%debug_str = trim(module_name)//"-"//trim(cmor_field_name)
      if (present(conversion)) cmor_diag%conversion_factor = conversion
      if (present(x_cell_method)) then
        call get_MOM_diag_axis_name(axes%handles(1), axis_name)
        call MOM_diag_field_add_attribute(fms_id, 'cell_methods', trim(axis_name)//':'//trim(x_cell_method))
      endif
      if (present(y_cell_method)) then
        call get_MOM_diag_axis_name(axes%handles(2), axis_name)
        call MOM_diag_field_add_attribute(fms_id, 'cell_methods', trim(axis_name)//':'//trim(y_cell_method))
      endif
      if (present(area_cell_method)) then
        call MOM_diag_field_add_attribute(fms_id, 'cell_methods', 'area:'//trim(area_cell_method))
      endif
    endif
  endif

  dimensions = ""
  if (axes%is_h_point)   dimensions = trim(dimensions)//" xh, yh,"
  if (axes%is_q_point)   dimensions = trim(dimensions)//" xq, yq,"
  if (axes%is_u_point)   dimensions = trim(dimensions)//" xq, yh,"
  if (axes%is_v_point)   dimensions = trim(dimensions)//" xh, yq,"
  if (len_trim(dimensions) > 0) dimensions = trim_trailing_commas(dimensions)

  ! Document diagnostics in list of available diagnostics
  if (is_root_pe() .and. diag_CS%available_diag_doc_unit > 0) then
    if (present(cmor_field_name)) then
      call log_available_diag(associated(diag), module_name, field_name, '', '', diag_CS, &
                              long_name, units, standard_name, &
                              variants="{"//trim(field_name)//","//trim(cmor_field_name)//"}", &
                              dimensions=dimensions)
    else
      call log_available_diag(associated(diag), module_name, field_name, '', '', diag_CS, &
                              long_name, units, standard_name, dimensions=dimensions)
    endif
  endif

  register_static_field = dm_id

end function register_MOM_IS_static_field

!> Add a description of an option to the documentation file
subroutine describe_option(opt_name, value, diag_CS)
  character(len=*),    intent(in) :: opt_name !< The name of the option
  character(len=*),    intent(in) :: value    !< The value of the option
  type(diag_ctrl), intent(in) :: diag_CS  !< Structure used to regulate diagnostic output

  ! Local variables
  character(len=480) :: mesg
  integer :: len_ind

  len_ind = len_trim(value)

  mesg = "    ! "//trim(opt_name)//": "//trim(value)
  write(diag_CS%available_diag_doc_unit, '(a)') trim(mesg)
end subroutine describe_option

!> Initialize the MOM_IS diag_mediator and opens the available diagnostics file, if appropriate.
subroutine MOM_IS_diag_mediator_init(G, US, param_file, diag_cs, component, err_msg, &
                                  doc_file_dir)
  type(ocean_grid_type), target, intent(inout) :: G  !< The horizontal grid type
  type(unit_scale_type), target, intent(in) :: US !< A dimensional unit scaling type
  type(param_file_type),      intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl),            intent(inout) :: diag_cs !< A structure that is used to regulate diagnostic output
  character(len=*), optional, intent(in)    :: component !< An optional component name
  character(len=*), optional, intent(out)   :: err_msg !< A string for a returned error message
  character(len=*), optional, intent(in)    :: doc_file_dir !< A directory in which to create the file

  ! This subroutine initializes the diag_mediator and the diag_manager.
  ! The grid type should have its dimensions set by this point, but it
  ! is not necessary that the metrics and axis labels be set up yet.

  ! Local variables
  integer :: ios, i, new_unit
  logical :: opened, new_file
  character(len=8)   :: this_pe
  character(len=240) :: doc_file, doc_file_dflt, doc_path
  character(len=40)  :: doc_file_param
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40) :: mdl = "MOM_IS_diag_mediator" ! This module's name.
  character(len=32) :: filename_appendix = '' !fms appendix to filename for ensemble runs

  call MOM_diag_manager_init(err_msg=err_msg)

  id_clock_diag_mediator = cpu_clock_id('(Ice shelf diagnostics framework)', grain=CLOCK_MODULE)

  ! Allocate and initialize list of all diagnostics (and variants)
  allocate(diag_cs%diags(DIAG_ALLOC_CHUNK_SIZE))
  diag_cs%next_free_diag_id = 1
  do i=1, DIAG_ALLOC_CHUNK_SIZE
    call initialize_diag_type(diag_cs%diags(i))
  enddo

  diag_cs%show_call_tree = callTree_showQuery()

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, 'USE_INDEX_DIAGNOSTIC_AXES', diag_cs%index_space_axes, &
                 'If true, use a grid index coordinate convention for diagnostic axes. ',&
                 default=.false.)

  call get_param(param_file, mdl, 'DIAG_MISVAL', diag_cs%missing_value, &
                 'Set the default missing value to use for diagnostics.', &
                 units="various", default=-1.e34)
  call get_param(param_file, mdl, 'DIAG_AS_CHKSUM', diag_cs%diag_as_chksum, &
                 'Instead of writing diagnostics to the diag manager, write '//&
                 'a text file containing the checksum (bitcount) of the array.',  &
                 default=.false.)

  if (diag_cs%diag_as_chksum) &
    diag_cs%num_chksum_diags = 0

  ! Keep pointers to the grid for diagnostic checksums
  diag_cs%G => G
  diag_cs%US => US

  diag_cs%is = G%isc - (G%isd-1) ; diag_cs%ie = G%iec - (G%isd-1)
  diag_cs%js = G%jsc - (G%jsd-1) ; diag_cs%je = G%jec - (G%jsd-1)
  diag_cs%isd = G%isd ; diag_cs%ied = G%ied
  diag_cs%jsd = G%jsd ; diag_cs%jed = G%jed

  ! Initialize available diagnostic log file
  if (is_root_pe() .and. (diag_CS%available_diag_doc_unit < 0)) then
    if (present(component)) then
      doc_file_dflt = trim(component)//".available_diags"
      doc_file_param = trim(uppercase(component))//"_AVAILABLE_DIAGS_FILE"
    else
      write(this_pe,'(i6.6)') PE_here()
      doc_file_dflt = "MOM_IS.available_diags."//this_pe
      doc_file_param = "AVAILABLE_MOM_IS_DIAGS_FILE"
    endif
    call get_param(param_file, mdl, trim(doc_file_param), doc_file, &
                 "A file into which to write a list of all available "//&
                 "ice shelf diagnostics that can be included in a diag_table.", &
                 default=doc_file_dflt, do_not_log=(diag_CS%available_diag_doc_unit/=-1))
    if (len_trim(doc_file) > 0) then
      new_file = .true. ; if (diag_CS%available_diag_doc_unit /= -1) new_file = .false.
    ! Find an unused unit number.
      do new_unit=512,42,-1
        inquire( new_unit, opened=opened)
        if (.not.opened) exit
      enddo
      if (opened) call MOM_error(FATAL, &
          "diag_mediator_init failed to find an unused unit number.")

      doc_path = doc_file
      if (present(doc_file_dir)) then ; if (len_trim(doc_file_dir) > 0) then
        doc_path = trim(slasher(doc_file_dir))//trim(doc_file)
      endif ; endif

      diag_CS%available_diag_doc_unit = new_unit

      if (new_file) then
        open(diag_CS%available_diag_doc_unit, file=trim(doc_path), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='REPLACE', iostat=ios)
      else ! This file is being reopened, and should be appended.
        open(diag_CS%available_diag_doc_unit, file=trim(doc_path), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='OLD', position='APPEND', iostat=ios)
      endif
      inquire(diag_CS%available_diag_doc_unit, opened=opened)
      if ((.not.opened) .or. (ios /= 0)) then
        call MOM_error(FATAL, "Failed to open available diags file "//trim(doc_path)//".")
      endif
    endif
  endif

  if (is_root_pe() .and. (diag_CS%chksum_iounit < 0) .and. diag_CS%diag_as_chksum) then
    !write(this_pe,'(i6.6)') PE_here()
    !doc_file_dflt = "chksum_diag."//this_pe
    doc_file_dflt = "chksum_diag"
    call get_param(param_file, mdl, "CHKSUM_DIAG_FILE", doc_file, &
                 "A file into which to write all checksums of the "//&
                 "diagnostics listed in the diag_table.", &
                 default=doc_file_dflt, do_not_log=(diag_CS%chksum_iounit/=-1))

    call get_filename_appendix(filename_appendix)
    if (len_trim(filename_appendix) > 0) then
      doc_file = trim(doc_file) //'.'//trim(filename_appendix)
    endif
#ifdef STATSLABEL
    doc_file = trim(doc_file)//"."//trim(adjustl(STATSLABEL))
#endif

    if (len_trim(doc_file) > 0) then
      new_file = .true. ; if (diag_CS%chksum_iounit /= -1) new_file = .false.
    ! Find an unused unit number.
      do new_unit=512,42,-1
        inquire( new_unit, opened=opened)
        if (.not.opened) exit
      enddo
      if (opened) call MOM_error(FATAL, &
          "diag_mediator_init failed to find an unused unit number.")

      doc_path = doc_file
      if (present(doc_file_dir)) then ; if (len_trim(doc_file_dir) > 0) then
        doc_path = trim(slasher(doc_file_dir))//trim(doc_file)
      endif ; endif

      diag_CS%chksum_iounit = new_unit

      if (new_file) then
        open(diag_CS%chksum_iounit, file=trim(doc_path), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='REPLACE', iostat=ios)
      else ! This file is being reopened, and should be appended.
        open(diag_CS%chksum_iounit, file=trim(doc_path), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='OLD', position='APPEND', iostat=ios)
      endif
      inquire(diag_CS%chksum_iounit, opened=opened)
      if ((.not.opened) .or. (ios /= 0)) then
        call MOM_error(FATAL, "Failed to open checksum diags file "//trim(doc_path)//".")
      endif
    endif
  endif

  call diag_masks_set(G, diag_cs%missing_value, diag_cs)

end subroutine MOM_IS_diag_mediator_init

!> Sets up the 2d masks for native diagnostics
subroutine diag_masks_set(G, missing_value, diag_cs)
  type(ocean_grid_type), target, intent(in)    :: G   !< The horizontal grid type
  real,                          intent(in)    :: missing_value !< A fill value for missing points
  type(diag_ctrl),               intent(inout) :: diag_cs !< Structure used to regulate diagnostic output

  ! Local variables
  integer :: i, j

  ! 2d masks point to the model masks since they are identical
  diag_cs%mask2dT  => G%mask2dT
  diag_cs%mask2dBu => G%mask2dBu
  diag_cs%mask2dCu => G%mask2dCu
  diag_cs%mask2dCv => G%mask2dCv

  allocate(diag_cs%mask2dT_comp(G%isc:G%iec,G%jsc:G%jec))
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    diag_cs%mask2dT_comp(i,j) = diag_cs%mask2dT(i,j)
  enddo ; enddo

  diag_cs%missing_value = missing_value

end subroutine diag_masks_set

!> Prevent the registration of additional diagnostics, so that the creation of files can occur
subroutine MOM_IS_diag_mediator_close_registration(diag_CS)
  type(diag_ctrl), intent(inout) :: diag_CS !< Structure used to regulate diagnostic output

  if (diag_CS%available_diag_doc_unit > -1) then
    close(diag_CS%available_diag_doc_unit) ; diag_CS%available_diag_doc_unit = -2
  endif

end subroutine MOM_IS_diag_mediator_close_registration

!> Deallocate memory associated with the MOM_IS diag mediator
subroutine MOM_IS_diag_mediator_end(diag_CS)
  type(diag_ctrl), intent(inout) :: diag_CS !< Structure used to regulate diagnostic output

  ! Local variables
  type(diag_type), pointer :: diag, next_diag
  integer :: i

  if (diag_CS%available_diag_doc_unit > -1) then
    close(diag_CS%available_diag_doc_unit) ; diag_CS%available_diag_doc_unit = -3
  endif
  if (diag_CS%chksum_iounit > -1) then
    close(diag_CS%chksum_iounit) ; diag_CS%chksum_iounit = -3
  endif

  do i=1, diag_cs%next_free_diag_id - 1
    if (associated(diag_cs%diags(i)%next)) then
      next_diag => diag_cs%diags(i)%next
      do while (associated(next_diag))
        diag => next_diag
        next_diag => diag%next
        deallocate(diag)
      enddo
    endif
  enddo

  deallocate(diag_cs%diags)

  ! These points to arrays in the grid type, so they can not be deallocated here.
  if (associated(diag_cs%mask2dT))  diag_cs%mask2dT => NULL()
  if (associated(diag_cs%mask2dBu)) diag_cs%mask2dBu => NULL()
  if (associated(diag_cs%mask2dCu)) diag_cs%mask2dCu => NULL()
  if (associated(diag_cs%mask2dCv)) diag_cs%mask2dCv => NULL()
  if (associated(diag_cs%mask2dT_comp)) deallocate(diag_cs%mask2dT_comp)

end subroutine MOM_IS_diag_mediator_end

!> Returns a new diagnostic id, it may be necessary to expand the diagnostics array.
integer function get_new_diag_id(diag_cs)
  type(diag_ctrl), intent(inout) :: diag_cs !< Diagnostics control structure
  ! Local variables
  type(diag_type), dimension(:), allocatable :: tmp
  integer :: i

  if (diag_cs%next_free_diag_id > size(diag_cs%diags)) then
    call assert(diag_cs%next_free_diag_id - size(diag_cs%diags) == 1, &
                'get_new_diag_id: inconsistent diag id')

    ! Increase the size of diag_cs%diags and copy data over.
    ! Do not use move_alloc() because it is not supported by Fortran 90
    allocate(tmp(size(diag_cs%diags)))
    tmp(:) = diag_cs%diags(:)
    deallocate(diag_cs%diags)
    allocate(diag_cs%diags(size(tmp) + DIAG_ALLOC_CHUNK_SIZE))
    diag_cs%diags(1:size(tmp)) = tmp(:)
    deallocate(tmp)

    ! Initialize new part of the diag array.
    do i=diag_cs%next_free_diag_id, size(diag_cs%diags)
      call initialize_diag_type(diag_cs%diags(i))
    enddo
  endif

  get_new_diag_id = diag_cs%next_free_diag_id
  diag_cs%next_free_diag_id = diag_cs%next_free_diag_id + 1

end function get_new_diag_id

!> Initializes a diag_type (used after allocating new memory)
subroutine initialize_diag_type(diag)
  type(diag_type), intent(inout) :: diag !< diag_type to be initialized

  diag%in_use = .false.
  diag%fms_diag_id = -1
  diag%axes => null()
  diag%next => null()
  diag%conversion_factor = 0.

end subroutine initialize_diag_type

!> Make a new diagnostic. Either use memory which is in the array of 'primary'
!! diagnostics, or if that is in use, insert it to the list of secondary diags.
subroutine alloc_diag_with_id(diag_id, diag_cs, diag)
  integer,                 intent(in   ) :: diag_id !< id for the diagnostic
  type(diag_ctrl), target, intent(inout) :: diag_cs !< structure used to regulate diagnostic output
  type(diag_type),         pointer       :: diag    !< structure representing a diagnostic (inout)

  type(diag_type), pointer :: tmp => NULL()

  if (.not. diag_cs%diags(diag_id)%in_use) then
    diag => diag_cs%diags(diag_id)
  else
    allocate(diag)
    tmp => diag_cs%diags(diag_id)%next
    diag_cs%diags(diag_id)%next => diag
    diag%next => tmp
  endif
  diag%in_use = .true.

end subroutine alloc_diag_with_id

!> Log a diagnostic to the available diagnostics file.
subroutine log_available_diag(used, module_name, field_name, cell_methods_string, comment, &
                              diag_CS, long_name, units, standard_name, variants, dimensions)
  logical,          intent(in) :: used !< Whether this diagnostic was in the diag_table or not
  character(len=*), intent(in) :: module_name !< Name of the diagnostic module
  character(len=*), intent(in) :: field_name !< Name of this diagnostic field
  character(len=*), intent(in) :: cell_methods_string !< The spatial component of the CF cell_methods attribute
  character(len=*), intent(in) :: comment !< A comment to append after [Used|Unused]
  type(diag_ctrl), intent(in) :: diag_CS  !< The diagnotics control structure
  character(len=*), optional, intent(in) :: dimensions !< Descriptor of the horizontal and vertical dimensions
  character(len=*), optional, intent(in) :: long_name !< CF long name of diagnostic
  character(len=*), optional, intent(in) :: units !< Units for diagnostic
  character(len=*), optional, intent(in) :: standard_name !< CF standardized name of diagnostic
  character(len=*), optional, intent(in) :: variants !< Alternate modules and variable names for
                                                     !! this diagnostic and derived diagnostics
  ! Local variables
  character(len=240) :: mesg

  if (used) then
    mesg = '"'//trim(field_name)//'"  [Used]'
  else
    mesg = '"'//trim(field_name)//'"  [Unused]'
  endif
  if (len(trim((comment)))>0) then
    write(diag_CS%available_diag_doc_unit, '(a,1x,"(",a,")")') trim(mesg),trim(comment)
  else
    write(diag_CS%available_diag_doc_unit, '(a)') trim(mesg)
  endif
  call describe_option("modules", module_name, diag_CS)
  if (present(dimensions)) then ; if (len(trim(dimensions)) > 0) then
    call describe_option("dimensions", dimensions, diag_CS)
  endif ; endif
  if (present(long_name)) call describe_option("long_name", long_name, diag_CS)
  if (present(units)) call describe_option("units", units, diag_CS)
  if (present(standard_name)) &
    call describe_option("standard_name", standard_name, diag_CS)
  if (len(trim((cell_methods_string)))>0) &
    call describe_option("cell_methods", trim(cell_methods_string), diag_CS)
  if (present(variants)) then ; if (len(trim(variants)) > 0) then
    call describe_option("variants", variants, diag_CS)
  endif ; endif
end subroutine log_available_diag

!> Log the diagnostic chksum to the chksum diag file
subroutine log_chksum_diag(docunit, description, chksum)
  integer,          intent(in) :: docunit     !< Handle of the log file
  character(len=*), intent(in) :: description !< Name of the diagnostic module
  integer,          intent(in) :: chksum      !< chksum of the diagnostic

  write(docunit, '(a,1x,i9.8)') description, chksum
  flush(docunit)

end subroutine log_chksum_diag

!> Fakes a register of a diagnostic to find out if an obsolete
!! parameter appears in the diag_table.
logical function found_in_diagtable(diag, varName)
  type(diag_ctrl), intent(in) :: diag     !< A structure used to control diagnostics.
  character(len=*),    intent(in) :: varName  !< The obsolete diagnostic name
  ! Local
  integer :: handle ! Integer handle returned from diag_manager

  ! We use register_static_field_fms() instead of register_static_field() so
  ! that the diagnostic does not appear in the available diagnostics list.
  handle = register_static_field_infra('ice_shelf_model', varName, diag%axesT1%handles)

  found_in_diagtable = (handle>0)

end function found_in_diagtable

!> Finishes the diag manager reduction methods as needed for the time_step
subroutine MOM_IS_diag_send_complete()
  call diag_send_complete_infra()
end subroutine MOM_IS_diag_send_complete

end module MOM_IS_diag_mediator
