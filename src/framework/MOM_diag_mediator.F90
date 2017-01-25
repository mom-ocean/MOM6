module MOM_diag_mediator

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

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*    The subroutines here provide convenient wrappers to the fms      *
!*  diag_manager interfaces with additional diagnostic capabilies.     *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_coms,             only : PE_here
use MOM_cpu_clock,        only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,        only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_error_handler,    only : MOM_error, FATAL, is_root_pe, assert
use MOM_file_parser,      only : get_param, log_version, param_file_type
use MOM_grid,             only : ocean_grid_type
use MOM_io,               only : slasher, vardesc, query_vardesc, mom_read_data
use MOM_safe_alloc,       only : safe_alloc_ptr, safe_alloc_alloc
use MOM_string_functions, only : lowercase
use MOM_time_manager,     only : time_type
use MOM_verticalGrid,     only : verticalGrid_type
use MOM_EOS,              only : EOS_type
use MOM_diag_remap,       only : diag_remap_ctrl
use MOM_diag_remap,       only : diag_remap_update
use MOM_diag_remap,       only : diag_remap_init, diag_remap_end, diag_remap_do_remap
use MOM_diag_remap,       only : vertically_reintegrate_diag_field, vertically_interpolate_diag_field
use MOM_diag_remap,       only : diag_remap_configure_axes, diag_remap_axes_configured
use MOM_diag_remap,       only : diag_remap_get_axes_info, diag_remap_set_active
use MOM_diag_remap,       only : diag_remap_diag_registration_closed
use MOM_diag_remap,       only : horizontally_average_diag_field

use diag_axis_mod, only : get_diag_axis_name
use diag_manager_mod, only : diag_manager_init, diag_manager_end
use diag_manager_mod, only : send_data, diag_axis_init, diag_field_add_attribute
! The following module is needed for PGI since the following line does not compile with PGI 6.5.0
! was: use diag_manager_mod, only : register_diag_field_fms=>register_diag_field
use MOM_diag_manager_wrapper, only : register_diag_field_fms
use diag_manager_mod, only : register_static_field_fms=>register_static_field
use diag_manager_mod, only : get_diag_field_id_fms=>get_diag_field_id
use diag_manager_mod, only : DIAG_FIELD_NOT_FOUND

implicit none ; private

#define __DO_SAFETY_CHECKS__
#define IMPLIES(A, B) ((.not. (A)) .or. (B))

public set_axes_info, post_data, register_diag_field, time_type
public post_data_1d_k
public safe_alloc_ptr, safe_alloc_alloc
public enable_averaging, disable_averaging, query_averaging_enabled
public diag_mediator_init, diag_mediator_end, set_diag_mediator_grid
public diag_mediator_infrastructure_init
public diag_mediator_close_registration, get_diag_time_end
public diag_axis_init, ocean_register_diag, register_static_field
public register_scalar_field
public define_axes_group, diag_masks_set
public diag_register_area_ids
public diag_register_volume_ids
public diag_set_state_ptrs, diag_update_remap_grids

interface post_data
  module procedure post_data_3d, post_data_2d, post_data_0d
end interface post_data

!> A group of 1D axes that comprise a 1D/2D/3D mesh
type, public :: axes_grp
  character(len=15) :: id   !< The id string for this particular combination of handles.
  integer           :: rank !< Number of dimensions in the list of axes.
  integer, dimension(:), allocatable :: handles !< Handles to 1D axes.
  type(diag_ctrl), pointer :: diag_cs => null() !< Circular link back to the main diagnostics control structure
                                                !! (Used to avoid passing said structure into every possible call).
  ! ID's for cell_methods
  character(len=9) :: x_cell_method = '' !< Default nature of data representation, if axes group includes x-direction.
  character(len=9) :: y_cell_method = '' !< Default nature of data representation, if axes group includes y-direction.
  character(len=9) :: v_cell_method = '' !< Default nature of data representation, if axes group includes vertical direction.
  ! For remapping
  integer :: nz = 0 !< Vertical dimension of diagnostic
  integer :: vertical_coordinate_number = 0 !< Index of the corresponding diag_remap_ctrl for this axis group
  ! For detecting position on the grid
  logical :: is_h_point = .false. !< If true, indicates that this axes group is for an h-point located field.
  logical :: is_q_point = .false. !< If true, indicates that this axes group is for a q-point located field.
  logical :: is_u_point = .false. !< If true, indicates that this axes group is for a u-point located field.
  logical :: is_v_point = .false. !< If true, indicates that this axes group is for a v-point located field.
  logical :: is_layer = .false. !< If true, indicates that this axes group is for a layer vertically-located field.
  logical :: is_interface = .false. !< If true, indicates that this axes group is for an interface vertically-located field.
  logical :: is_native = .true. !< If true, indicates that this axes group is for a native model grid. False for any other
                                !! grid. Used for rank>2.
  logical :: needs_remapping = .false. !< If true, indicates that this axes group is for a intensive layer-located field
                                       !! that must be remapped to these axes. Used for rank>2.
  logical :: needs_interpolating = .false. !< If true, indicates that this axes group is for a sampled interface-located field
                                           !! that must be interpolated to these axes. Used for rank>2.
  ! For horizontally averaged diagnositcs (applies to 2d and 3d fields only)
  type(axes_grp), pointer :: xyave_axes => null() !< The associated 1d axes for horizontall area-averaged diagnostics
  ! ID's for cell_measures
  integer :: id_area = -1 !< The diag_manager id for area to be used for cell_measure of variables with this axes_grp.
  integer :: id_volume = -1 !< The diag_manager id for volume to be used for cell_measure of variables with this axes_grp.
end type axes_grp

!> This type is used to represent a diagnostic at the diag_mediator level.
!! There can be both 'primary' and 'seconday' diagnostics. The primaries
!! reside in the diag_cs%diags array. They have an id which is an index
!! into this array. The secondaries are 'variations' on the primary diagnostic.
!! For example the CMOR diagnostics are secondary. The secondary diagnostics
!! are kept in a list with the primary diagnostic as the head.
type, private :: diag_type
  logical :: in_use !< True if this entry is being used.
  integer :: fms_diag_id !< Underlying FMS diag_manager id.
  integer :: fms_xyave_diag_id = -1 !< For a horizontally area-averaged diagnostic.
  character(32) :: debug_str = '' !< For FATAL errors and debugging.
  type(axes_grp), pointer :: axes => null()
  real, pointer, dimension(:,:)   :: mask2d => null()
  real, pointer, dimension(:,:,:) :: mask3d => null()
  type(diag_type), pointer :: next => null() !< Pointer to the next diag.
  real :: conversion_factor = 0. !< A factor to multiply data by before posting to FMS, if non-zero.
  logical :: v_extensive = .false. !< True for vertically extensive fields (vertically integrated). False for intensive (concentrations).
end type diag_type

!> The following data type a list of diagnostic fields an their variants,
!! as well as variables that control the handling of model output.
type, public :: diag_ctrl
  integer :: doc_unit = -1 !< The unit number of a diagnostic documentation file.
                           !! This file is open if doc_unit is > 0.

! The following fields are used for the output of the data.
  integer :: is, ie, js, je
  integer :: isd, ied, jsd, jed
  real :: time_int              !< The time interval in s for any fields
                                !! that are offered for averaging.
  type(time_type) :: time_end   !< The end time of the valid
                                !! interval for any offered field.
  logical :: ave_enabled = .false. !< True if averaging is enabled.

  ! The following are axis types defined for output.
  type(axes_grp) :: axesBL, axesTL, axesCuL, axesCvL
  type(axes_grp) :: axesBi, axesTi, axesCui, axesCvi
  type(axes_grp) :: axesB1, axesT1, axesCu1, axesCv1
  type(axes_grp) :: axesZi, axesZL

  ! Mask arrays for diagnostics
  real, dimension(:,:),   pointer :: mask2dT   => null()
  real, dimension(:,:),   pointer :: mask2dBu  => null()
  real, dimension(:,:),   pointer :: mask2dCu  => null()
  real, dimension(:,:),   pointer :: mask2dCv  => null()
  real, dimension(:,:,:), pointer :: mask3dTL  => null()
  real, dimension(:,:,:), pointer :: mask3dBL => null()
  real, dimension(:,:,:), pointer :: mask3dCuL => null()
  real, dimension(:,:,:), pointer :: mask3dCvL => null()
  real, dimension(:,:,:), pointer :: mask3dTi  => null()
  real, dimension(:,:,:), pointer :: mask3dBi => null()
  real, dimension(:,:,:), pointer :: mask3dCui => null()
  real, dimension(:,:,:), pointer :: mask3dCvi => null()

! Space for diagnostics is dynamically allocated as it is needed.
! The chunk size is how much the array should grow on each new allocation.
#define DIAG_ALLOC_CHUNK_SIZE 100
  type(diag_type), dimension(:), allocatable :: diags
  integer :: next_free_diag_id

  !default missing value to be sent to ALL diagnostics registrations
  real :: missing_value = -1.0e+34

  !> Number of diagnostic vertical coordinates (remapped)
  integer :: num_diag_coords
  !> Control structure for each possible coordinate
  type(diag_remap_ctrl), dimension(:), allocatable :: diag_remap_cs

  !> Axes groups for each possible coordinate (these will all be 3D groups)
  type(axes_grp), dimension(:), allocatable :: remap_axesZL, remap_axesZi
  type(axes_grp), dimension(:), allocatable :: remap_axesTL, remap_axesBL, remap_axesCuL, remap_axesCvL
  type(axes_grp), dimension(:), allocatable :: remap_axesTi, remap_axesBi, remap_axesCui, remap_axesCvi

  ! Pointer to H, G and T&S needed for remapping
  real, dimension(:,:,:), pointer :: h => null()
  real, dimension(:,:,:), pointer :: T => null()
  real, dimension(:,:,:), pointer :: S => null()
  type(EOS_type),  pointer :: eqn_of_state => null()
  type(ocean_grid_type), pointer :: G => null()

#if defined(DEBUG) || defined(__DO_SAFETY_CHECKS__)
  ! Keep a copy of h so that we know whether it has changed. If it has then
  ! need the target grid for vertical remapping needs to have been updated.
  real, dimension(:,:,:), allocatable :: h_old
#endif

end type diag_ctrl

! CPU clocks
integer :: id_clock_diag_mediator, id_clock_diag_remap, id_clock_diag_grid_updates

contains

!> Sets up diagnostics axes
subroutine set_axes_info(G, GV, param_file, diag_cs, set_vertical)
  type(ocean_grid_type), intent(inout) :: G !< Ocean grid structure
  type(verticalGrid_type), intent(in)  :: GV !< ocean vertical grid structure
  type(param_file_type), intent(in)    :: param_file !< Parameter file structure
  type(diag_ctrl),       intent(inout) :: diag_cs !< Diagnostics control structure
  logical, optional,     intent(in)    :: set_vertical !< If true or missing, set up
                                                       !! vertical axes
  ! Local variables
  integer :: id_xq, id_yq, id_zl, id_zi, id_xh, id_yh
  integer :: i, k, nz
  real :: zlev(GV%ke), zinter(GV%ke+1)
  logical :: set_vert

  set_vert = .true. ; if (present(set_vertical)) set_vert = set_vertical

  if(G%symmetric) then
    id_xq = diag_axis_init('xq', G%gridLonB(G%isgB:G%iegB), G%x_axis_units, 'x', &
              'q point nominal longitude', Domain2=G%Domain%mpp_domain)
    id_yq = diag_axis_init('yq', G%gridLatB(G%jsgB:G%jegB), G%y_axis_units, 'y', &
              'q point nominal latitude', Domain2=G%Domain%mpp_domain)
  else
    id_xq = diag_axis_init('xq', G%gridLonB(G%isg:G%ieg), G%x_axis_units, 'x', &
              'q point nominal longitude', Domain2=G%Domain%mpp_domain)
    id_yq = diag_axis_init('yq', G%gridLatB(G%jsg:G%jeg), G%y_axis_units, 'y', &
              'q point nominal latitude', Domain2=G%Domain%mpp_domain)
  endif
  id_xh = diag_axis_init('xh', G%gridLonT(G%isg:G%ieg), G%x_axis_units, 'x', &
              'h point nominal longitude', Domain2=G%Domain%mpp_domain)
  id_yh = diag_axis_init('yh', G%gridLatT(G%jsg:G%jeg), G%y_axis_units, 'y', &
              'h point nominal latitude', Domain2=G%Domain%mpp_domain)

  if (set_vert) then
    nz = GV%ke
    zinter(1:nz+1) = GV%sInterface(1:nz+1)
    zlev(1:nz) = GV%sLayer(1:nz)
    id_zl = diag_axis_init('zl', zlev, trim(GV%zAxisUnits), 'z', &
                           'Layer '//trim(GV%zAxisLongName),     &
                           direction=GV%direction)
    id_zi = diag_axis_init('zi', zinter, trim(GV%zAxisUnits), 'z', &
                           'Interface '//trim(GV%zAxisLongName),   &
                           direction=GV%direction)
  else
    id_zl = -1 ; id_zi = -1
  endif

  ! Vertical axes for the interfaces and layers
  call define_axes_group(diag_cs, (/ id_zi /), diag_cs%axesZi, &
       v_cell_method='point', is_interface=.true.)
  call define_axes_group(diag_cs, (/ id_zL /), diag_cs%axesZL, &
       v_cell_method='mean', is_layer=.true.)

  ! Axis groupings for the model layers
  call define_axes_group(diag_cs, (/ id_xh, id_yh, id_zL /), diag_cs%axesTL, &
       x_cell_method='mean', y_cell_method='mean', v_cell_method='mean', &
       is_h_point=.true., is_layer=.true., xyave_axes=diag_cs%axesZL)
  call define_axes_group(diag_cs, (/ id_xq, id_yq, id_zL /), diag_cs%axesBL, &
       x_cell_method='point', y_cell_method='point', v_cell_method='mean', &
       is_q_point=.true., is_layer=.true.)
  call define_axes_group(diag_cs, (/ id_xq, id_yh, id_zL /), diag_cs%axesCuL, &
       x_cell_method='point', y_cell_method='mean', v_cell_method='mean', &
       is_u_point=.true., is_layer=.true., xyave_axes=diag_cs%axesZL)
  call define_axes_group(diag_cs, (/ id_xh, id_yq, id_zL /), diag_cs%axesCvL, &
       x_cell_method='mean', y_cell_method='point', v_cell_method='mean', &
       is_v_point=.true., is_layer=.true., xyave_axes=diag_cs%axesZL)

  ! Axis groupings for the model interfaces
  call define_axes_group(diag_cs, (/ id_xh, id_yh, id_zi /), diag_cs%axesTi, &
       x_cell_method='mean', y_cell_method='mean', v_cell_method='point', &
       is_h_point=.true., is_interface=.true., xyave_axes=diag_cs%axesZi)
  call define_axes_group(diag_cs, (/ id_xq, id_yq, id_zi /), diag_cs%axesBi, &
       x_cell_method='point', y_cell_method='point', v_cell_method='point', &
       is_q_point=.true., is_interface=.true.)
  call define_axes_group(diag_cs, (/ id_xq, id_yh, id_zi /), diag_cs%axesCui, &
       x_cell_method='point', y_cell_method='mean', v_cell_method='point', &
       is_u_point=.true., is_interface=.true., xyave_axes=diag_cs%axesZi)
  call define_axes_group(diag_cs, (/ id_xh, id_yq, id_zi /), diag_cs%axesCvi, &
       x_cell_method='mean', y_cell_method='point', v_cell_method='point', &
       is_v_point=.true., is_interface=.true., xyave_axes=diag_cs%axesZi)

  ! Axis groupings for 2-D arrays
  call define_axes_group(diag_cs, (/ id_xh, id_yh /), diag_cs%axesT1, &
       x_cell_method='mean', y_cell_method='mean', is_h_point=.true.)
  call define_axes_group(diag_cs, (/ id_xq, id_yq /), diag_cs%axesB1, &
       x_cell_method='point', y_cell_method='point', is_q_point=.true.)
  call define_axes_group(diag_cs, (/ id_xq, id_yh /), diag_cs%axesCu1, &
       x_cell_method='point', y_cell_method='mean', is_u_point=.true.)
  call define_axes_group(diag_cs, (/ id_xh, id_yq /), diag_cs%axesCv1, &
       x_cell_method='mean', y_cell_method='point', is_v_point=.true.)

  if (diag_cs%num_diag_coords>0) then
    allocate(diag_cs%remap_axesZL(diag_cs%num_diag_coords))
    allocate(diag_cs%remap_axesTL(diag_cs%num_diag_coords))
    allocate(diag_cs%remap_axesBL(diag_cs%num_diag_coords))
    allocate(diag_cs%remap_axesCuL(diag_cs%num_diag_coords))
    allocate(diag_cs%remap_axesCvL(diag_cs%num_diag_coords))
    allocate(diag_cs%remap_axesZi(diag_cs%num_diag_coords))
    allocate(diag_cs%remap_axesTi(diag_cs%num_diag_coords))
    allocate(diag_cs%remap_axesBi(diag_cs%num_diag_coords))
    allocate(diag_cs%remap_axesCui(diag_cs%num_diag_coords))
    allocate(diag_cs%remap_axesCvi(diag_cs%num_diag_coords))
  endif

  do i=1, diag_cs%num_diag_coords
    ! For each possible diagnostic coordinate
    call diag_remap_configure_axes(diag_cs%diag_remap_cs(i), GV, param_file)

    ! This vertical coordinate has been configured so can be used.
    if (diag_remap_axes_configured(diag_cs%diag_remap_cs(i))) then

      ! This fetches the 1D-axis id for layers and interfaces and overwrite
      ! id_zl and id_zi from above. It also returns the number of layers.
      call diag_remap_get_axes_info(diag_cs%diag_remap_cs(i), nz, id_zL, id_zi)

      ! Axes for z layers
      call define_axes_group(diag_cs, (/ id_zL /), diag_cs%remap_axesZL(i), &
           nz=nz, vertical_coordinate_number=i, &
           v_cell_method='mean', &
           is_h_point=.true., is_layer=.true., is_native=.false., needs_remapping=.true.)
      call define_axes_group(diag_cs, (/ id_xh, id_yh, id_zL /), diag_cs%remap_axesTL(i), &
           nz=nz, vertical_coordinate_number=i, &
           x_cell_method='mean', y_cell_method='mean', v_cell_method='mean', &
           is_h_point=.true., is_layer=.true., is_native=.false., needs_remapping=.true., &
           xyave_axes=diag_cs%remap_axesZL(i))

       !! \note Remapping for B points is not yet implemented so needs_remapping is not provided for remap_axesBL
      call define_axes_group(diag_cs, (/ id_xq, id_yq, id_zL /), diag_cs%remap_axesBL(i), &
           nz=nz, vertical_coordinate_number=i, &
           x_cell_method='point', y_cell_method='point', v_cell_method='mean', &
           is_q_point=.true., is_layer=.true., is_native=.false.)

      call define_axes_group(diag_cs, (/ id_xq, id_yh, id_zL /), diag_cs%remap_axesCuL(i), &
           nz=nz, vertical_coordinate_number=i, &
           x_cell_method='point', y_cell_method='mean', v_cell_method='mean', &
           is_u_point=.true., is_layer=.true., is_native=.false., needs_remapping=.true., &
           xyave_axes=diag_cs%remap_axesZL(i))

      call define_axes_group(diag_cs, (/ id_xh, id_yq, id_zL /), diag_cs%remap_axesCvL(i), &
           nz=nz, vertical_coordinate_number=i, &
           x_cell_method='mean', y_cell_method='point', v_cell_method='mean', &
           is_v_point=.true., is_layer=.true., is_native=.false., needs_remapping=.true., &
           xyave_axes=diag_cs%remap_axesZL(i))

      ! Axes for z interfaces
      call define_axes_group(diag_cs, (/ id_zi /), diag_cs%remap_axesZi(i), &
           nz=nz, vertical_coordinate_number=i, &
           v_cell_method='point', &
           is_h_point=.true., is_interface=.true., is_native=.false., needs_interpolating=.true.)
      call define_axes_group(diag_cs, (/ id_xh, id_yh, id_zi /), diag_cs%remap_axesTi(i), &
           nz=nz, vertical_coordinate_number=i, &
           x_cell_method='mean', y_cell_method='mean', v_cell_method='point', &
           is_h_point=.true., is_interface=.true., is_native=.false., needs_interpolating=.true., &
           xyave_axes=diag_cs%remap_axesZi(i))

      !! \note Remapping for B points is not yet implemented so needs_remapping is not provided for remap_axesBi
      call define_axes_group(diag_cs, (/ id_xq, id_yq, id_zi /), diag_cs%remap_axesBi(i), &
           nz=nz, vertical_coordinate_number=i, &
           x_cell_method='point', y_cell_method='point', v_cell_method='point', &
           is_q_point=.true., is_interface=.true., is_native=.false.)

      call define_axes_group(diag_cs, (/ id_xq, id_yh, id_zi /), diag_cs%remap_axesCui(i), &
           nz=nz, vertical_coordinate_number=i, &
           x_cell_method='point', y_cell_method='mean', v_cell_method='point', &
           is_u_point=.true., is_interface=.true., is_native=.false., &
           needs_interpolating=.true., xyave_axes=diag_cs%remap_axesZi(i))

      call define_axes_group(diag_cs, (/ id_xh, id_yq, id_zi /), diag_cs%remap_axesCvi(i), &
           nz=nz, vertical_coordinate_number=i, &
           x_cell_method='mean', y_cell_method='point', v_cell_method='point', &
           is_v_point=.true., is_interface=.true., is_native=.false., &
           needs_interpolating=.true., xyave_axes=diag_cs%remap_axesZi(i))
    endif
  enddo

end subroutine set_axes_info

!> Attaches the id of cell areas to axes groups for use with cell_measures
subroutine diag_register_area_ids(diag_cs, id_area_t, id_area_q)
  type(diag_ctrl),   intent(inout) :: diag_cs   !< Diagnostics control structure
  integer, optional, intent(in)    :: id_area_t !< Diag_mediator id for area of h-cells
  integer, optional, intent(in)    :: id_area_q !< Diag_mediator id for area of q-cells
  ! Local variables
  integer :: fms_id, i
  if (present(id_area_t)) then
    fms_id = diag_cs%diags(id_area_t)%fms_diag_id
    diag_cs%axesT1%id_area = fms_id
    diag_cs%axesTi%id_area = fms_id
    diag_cs%axesTL%id_area = fms_id
    do i=1, diag_cs%num_diag_coords
      diag_cs%remap_axesTL(i)%id_area = fms_id
      ! Note to AJA: why am I not doing TZi too?
    enddo
  endif
  if (present(id_area_q)) then
    fms_id = diag_cs%diags(id_area_q)%fms_diag_id
    diag_cs%axesB1%id_area = fms_id
    diag_cs%axesBi%id_area = fms_id
    diag_cs%axesBL%id_area = fms_id
    do i=1, diag_cs%num_diag_coords
      diag_cs%remap_axesBL(i)%id_area = fms_id
    enddo
  endif
end subroutine diag_register_area_ids

!> Attaches the id of cell volumes to axes groups for use with cell_measures
subroutine diag_register_volume_ids(diag_cs, id_vol_t)
  type(diag_ctrl),   intent(inout) :: diag_cs   !< Diagnostics control structure
  integer, optional, intent(in)    :: id_vol_t !< Diag_manager id for volume of h-cells
  ! Local variables
  integer :: fms_id
  if (present(id_vol_t)) then
    fms_id = diag_cs%diags(id_vol_t)%fms_diag_id
    call MOM_error(FATAL,"diag_register_volume_ids: not implemented yet!")
  endif
end subroutine diag_register_volume_ids

!> Defines a group of "axes" from list of handles
subroutine define_axes_group(diag_cs, handles, axes, nz, vertical_coordinate_number, &
                             x_cell_method, y_cell_method, v_cell_method, &
                             is_h_point, is_q_point, is_u_point, is_v_point, &
                             is_layer, is_interface, &
                             is_native, needs_remapping, needs_interpolating, &
                             xyave_axes)
  type(diag_ctrl), target,    intent(in)  :: diag_cs !< Diagnostics control structure
  integer, dimension(:),      intent(in)  :: handles !< A list of 1D axis handles
  type(axes_grp),             intent(out) :: axes    !< The group of 1D axes
  integer,          optional, intent(in)  :: nz      !< Number of layers in this diagnostic grid
  integer,          optional, intent(in)  :: vertical_coordinate_number !< Index number for vertical coordinate
  character(len=*), optional, intent(in)  :: x_cell_method !< A x-direction cell method used to construct the "cell_methods" attribute in CF convention
  character(len=*), optional, intent(in)  :: y_cell_method !< A y-direction cell method used to construct the "cell_methods" attribute in CF convention
  character(len=*), optional, intent(in)  :: v_cell_method !< A vertical direction cell method used to construct the "cell_methods" attribute in CF convention
  logical,          optional, intent(in)  :: is_h_point !< If true, indicates this axes group for h-point located fields
  logical,          optional, intent(in)  :: is_q_point !< If true, indicates this axes group for q-point located fields
  logical,          optional, intent(in)  :: is_u_point !< If true, indicates this axes group for u-point located fields
  logical,          optional, intent(in)  :: is_v_point !< If true, indicates this axes group for v-point located fields
  logical,          optional, intent(in)  :: is_layer !< If true, indicates that this axes group is for a layer vertically-located field.
  logical,          optional, intent(in)  :: is_interface !< If true, indicates that this axes group is for an interface vertically-located field.
  logical,          optional, intent(in)  :: is_native !< If true, indicates that this axes group is for a native model grid. False for any other grid.
  logical,          optional, intent(in)  :: needs_remapping !< If true, indicates that this axes group is for a intensive layer-located field
                                                             !! that must be remapped to these axes. Used for rank>2.
  logical,          optional, intent(in)  :: needs_interpolating !< If true, indicates that this axes group is for a sampled interface-located field
                                                                 !! that must be interpolated to these axes. Used for rank>2.
  type(axes_grp),   optional, target      :: xyave_axes !< The corresponding axes group for horizontally area-average diagnostics
  ! Local variables
  integer :: n

  n = size(handles)
  if (n<1 .or. n>3) call MOM_error(FATAL, "define_axes_group: wrong size for list of handles!")
  allocate( axes%handles(n) )
  axes%id = i2s(handles, n) ! Identifying string
  axes%rank = n
  axes%handles(:) = handles(:)
  axes%diag_cs => diag_cs ! A [circular] link back to the diag_cs structure
  if (present(x_cell_method)) then
    if (axes%rank<2) call MOM_error(FATAL, 'define_axes_group: ' // &
                                           'Can not set x_cell_method for rank<2.')
    axes%x_cell_method = trim(x_cell_method)
  else
    axes%x_cell_method = ''
  endif
  if (present(y_cell_method)) then
    if (axes%rank<2) call MOM_error(FATAL, 'define_axes_group: ' // &
                                           'Can not set y_cell_method for rank<2.')
    axes%y_cell_method = trim(y_cell_method)
  else
    axes%y_cell_method = ''
  endif
  if (present(v_cell_method)) then
    if (axes%rank/=1 .and. axes%rank/=3) call MOM_error(FATAL, 'define_axes_group: ' // &
                                           'Can not set v_cell_method for rank<>1 or 3.')
    axes%v_cell_method = trim(v_cell_method)
  else
    axes%v_cell_method = ''
  endif
  if (present(nz)) axes%nz = nz
  if (present(vertical_coordinate_number)) axes%vertical_coordinate_number = vertical_coordinate_number
  if (present(is_h_point)) axes%is_h_point = is_h_point
  if (present(is_q_point)) axes%is_q_point = is_q_point
  if (present(is_u_point)) axes%is_u_point = is_u_point
  if (present(is_v_point)) axes%is_v_point = is_v_point
  if (present(is_layer)) axes%is_layer = is_layer
  if (present(is_interface)) axes%is_interface = is_interface
  if (present(is_native)) axes%is_native = is_native
  if (present(needs_remapping)) axes%needs_remapping = needs_remapping
  if (present(needs_interpolating)) axes%needs_interpolating = needs_interpolating
  if (present(xyave_axes)) axes%xyave_axes => xyave_axes

end subroutine define_axes_group

subroutine set_diag_mediator_grid(G, diag_cs)
  type(ocean_grid_type), intent(inout) :: G
  type(diag_ctrl),  intent(inout) :: diag_cs

! Arguments:
!  (inout)    G   - ocean grid structure
!  (inout)   diag - structure used to regulate diagnostic output

  diag_cs%is = G%isc - (G%isd-1) ; diag_cs%ie = G%iec - (G%isd-1)
  diag_cs%js = G%jsc - (G%jsd-1) ; diag_cs%je = G%jec - (G%jsd-1)
  diag_cs%isd = G%isd ; diag_cs%ied = G%ied
  diag_cs%jsd = G%jsd ; diag_cs%jed = G%jed

end subroutine set_diag_mediator_grid

subroutine post_data_0d(diag_field_id, field, diag_cs, is_static)
  integer,           intent(in) :: diag_field_id
  real,              intent(in) :: field
  type(diag_ctrl), target, intent(in) :: diag_cs
  logical, optional, intent(in) :: is_static

! Arguments:
!  (in) diag_field_id  - the id for an output variable returned by a
!                            previous call to register_diag_field.
!  (in)      field     - 0-d array being offered for output or averaging.
!  (inout)   diag_cs - structure used to regulate diagnostic output.
!  (in,opt)  is_static - If true, this is a static field that is always offered.
!  (in,opt)  mask      - If present, use this real array as the data mask.

  logical :: used, is_stat
  type(diag_type), pointer :: diag => null()

  if (id_clock_diag_mediator>0) call cpu_clock_begin(id_clock_diag_mediator)
  is_stat = .false. ; if (present(is_static)) is_stat = is_static

  ! Iterate over list of diag 'variants', e.g. CMOR aliases, call send_data
  ! for each one.
  call assert(diag_field_id < diag_cs%next_free_diag_id, &
              'post_data_0d: Unregistered diagnostic id')
  diag => diag_cs%diags(diag_field_id)
  do while (associated(diag))
    if (is_stat) then
      used = send_data(diag%fms_diag_id, field)
    elseif (diag_cs%ave_enabled) then
      used = send_data(diag%fms_diag_id, field, diag_cs%time_end)
    endif
    diag => diag%next
  enddo

  if (id_clock_diag_mediator>0) call cpu_clock_end(id_clock_diag_mediator)
end subroutine post_data_0d

subroutine post_data_1d_k(diag_field_id, field, diag_cs, is_static)
  integer,           intent(in) :: diag_field_id
  real,              intent(in) :: field(:)
  type(diag_ctrl), target, intent(in) :: diag_cs
  logical, optional, intent(in) :: is_static

! Arguments:
!  (in) diag_field_id - id for an output variable returned by a
!                       previous call to register_diag_field.
!  (in)         field - 3-d array being offered for output or averaging
!  (inout)    diag_cs - structure used to regulate diagnostic output
!  (in)        static - If true, this is a static field that is always offered.

  logical :: used  ! The return value of send_data is not used for anything.
  logical :: is_stat
  integer :: isv, iev, jsv, jev
  type(diag_type), pointer :: diag => null()

  if (id_clock_diag_mediator>0) call cpu_clock_begin(id_clock_diag_mediator)
  is_stat = .false. ; if (present(is_static)) is_stat = is_static

  ! Iterate over list of diag 'variants', e.g. CMOR aliases.
  call assert(diag_field_id < diag_cs%next_free_diag_id, &
              'post_data_1d_k: Unregistered diagnostic id')
  diag => diag_cs%diags(diag_field_id)
  do while (associated(diag))
    if (is_stat) then
      used = send_data(diag%fms_diag_id, field)
    elseif (diag_cs%ave_enabled) then
      used = send_data(diag%fms_diag_id, field, diag_cs%time_end, weight=diag_cs%time_int)
    endif
    diag => diag%next
  enddo

  if (id_clock_diag_mediator>0) call cpu_clock_end(id_clock_diag_mediator)
end subroutine post_data_1d_k

subroutine post_data_2d(diag_field_id, field, diag_cs, is_static, mask)
  integer,           intent(in) :: diag_field_id
  real,              intent(in) :: field(:,:)
  type(diag_ctrl), target, intent(in) :: diag_cs
  logical, optional, intent(in) :: is_static
  real,    optional, intent(in) :: mask(:,:)

! Arguments:
!  (in) diag_field_id  - id for an output variable returned by a
!                        previous call to register_diag_field.
!  (in)         field  - 2-d array being offered for output or averaging.
!  (inout)    diag_cs  - structure used to regulate diagnostic output.
!  (in,opt)  is_static - If true, this is a static field that is always offered.
!  (in,opt)   mask     - If present, use this real array as the data mask.

  type(diag_type), pointer :: diag => null()

  if (id_clock_diag_mediator>0) call cpu_clock_begin(id_clock_diag_mediator)

  ! Iterate over list of diag 'variants' (e.g. CMOR aliases) and post each.
  call assert(diag_field_id < diag_cs%next_free_diag_id, &
              'post_data_2d: Unregistered diagnostic id')
  diag => diag_cs%diags(diag_field_id)
  do while (associated(diag))
     call post_data_2d_low(diag, field, diag_cs, is_static, mask)
    diag => diag%next
  enddo

  if (id_clock_diag_mediator>0) call cpu_clock_end(id_clock_diag_mediator)
end subroutine post_data_2d

subroutine post_data_2d_low(diag, field, diag_cs, is_static, mask)
  type(diag_type),   intent(in) :: diag
  real,    target,   intent(in) :: field(:,:)
  type(diag_ctrl), intent(in) :: diag_cs
  logical, optional, intent(in) :: is_static
  real,    optional, intent(in) :: mask(:,:)

! Arguments:
!  (in) diag          - structure representing the diagnostic to post
!  (in)        field  - 2-d array being offered for output or averaging
!  (inout) diag_cs    - structure used to regulate diagnostic output
!  (in,opt) is_static - If true, this is a static field that is always offered.
!  (in,opt)  mask     - If present, use this real array as the data mask.

  real, dimension(:,:), pointer :: locfield => NULL()
  logical :: used, is_stat
  integer :: isv, iev, jsv, jev

  is_stat = .false. ; if (present(is_static)) is_stat = is_static

  ! Determine the propery array indices, noting that because of the (:,:)
  ! declaration of field, symmetric arrays are using a SW-grid indexing,
  ! but non-symmetric arrays are using a NE-grid indexing.  Send_data
  ! actually only uses the difference between ie and is to determine
  ! the output data size and assumes that halos are symmetric.
  isv = diag_cs%is ; iev = diag_cs%ie ; jsv = diag_cs%js ; jev = diag_cs%je

  if ( size(field,1) == diag_cs%ied-diag_cs%isd +1 ) then
    isv = diag_cs%is ; iev = diag_cs%ie        ! Data domain
  elseif ( size(field,1) == diag_cs%ied-diag_cs%isd +2 ) then
    isv = diag_cs%is ; iev = diag_cs%ie+1      ! Symmetric data domain
  elseif ( size(field,1) == diag_cs%ie-diag_cs%is +1 ) then
    isv = 1 ; iev = diag_cs%ie + 1-diag_cs%is  ! Computational domain
  elseif ( size(field,1) == diag_cs%ie-diag_cs%is +2 ) then
    isv = 1 ; iev = diag_cs%ie + 2-diag_cs%is  ! Symmetric computational domain
  else
    call MOM_error(FATAL,"post_data_2d_low: peculiar size in i-direction")
  endif
  if ( size(field,2) == diag_cs%jed-diag_cs%jsd +1 ) then
    jsv = diag_cs%js ; jev = diag_cs%je        ! Data domain
  elseif ( size(field,2) == diag_cs%jed-diag_cs%jsd +2 ) then
    jsv = diag_cs%js ; jev = diag_cs%je+1      ! Symmetric data domain
  elseif ( size(field,2) == diag_cs%je-diag_cs%js +1 ) then
    jsv = 1 ; jev = diag_cs%je + 1-diag_cs%js  ! Computational domain
  elseif ( size(field,1) == diag_cs%je-diag_cs%js +2 ) then
    jsv = 1 ; jev = diag_cs%je + 2-diag_cs%js  ! Symmetric computational domain
  else
    call MOM_error(FATAL,"post_data_2d_low: peculiar size in j-direction")
  endif

  if (diag%conversion_factor/=0.) then
    allocate( locfield( lbound(field,1):ubound(field,1), lbound(field,2):ubound(field,2) ) )
    locfield(isv:iev,jsv:jev) = field(isv:iev,jsv:jev) * diag%conversion_factor
  else
    locfield => field
  endif

  if (is_stat) then
    if (present(mask)) then
      call assert(size(locfield) == size(mask), &
          'post_data_2d_low is_stat: mask size mismatch: '//diag%debug_str)
      used = send_data(diag%fms_diag_id, locfield, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, rmask=mask)
   !elseif(associated(diag%mask2d)) then
   !  used = send_data(diag%fms_diag_id, locfield, &
   !                   is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, rmask=diag%mask2d)
    else
      used = send_data(diag%fms_diag_id, locfield, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev)
    endif
  elseif (diag_cs%ave_enabled) then
    if (present(mask)) then
      call assert(size(locfield) == size(mask), &
          'post_data_2d_low: mask size mismatch: '//diag%debug_str)
      used = send_data(diag%fms_diag_id, locfield, diag_cs%time_end, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                       weight=diag_cs%time_int, rmask=mask)
    elseif(associated(diag%mask2d)) then
      used = send_data(diag%fms_diag_id, locfield, diag_cs%time_end, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                       weight=diag_cs%time_int, rmask=diag%mask2d)
    else
      used = send_data(diag%fms_diag_id, locfield, diag_cs%time_end, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                       weight=diag_cs%time_int)
    endif
  endif
  if (diag%conversion_factor/=0.) deallocate( locfield )

end subroutine post_data_2d_low

subroutine post_data_3d(diag_field_id, field, diag_cs, is_static, mask)

  integer,           intent(in) :: diag_field_id
  real,              intent(in) :: field(:,:,:)
  type(diag_ctrl), target, intent(in) :: diag_cs
  logical, optional, intent(in) :: is_static
  real,    optional, intent(in) :: mask(:,:,:)

! Arguments:
!  (in) diag_field_id - id for an output variable returned by a
!                       previous call to register_diag_field.
!  (in)         field - 3-d array being offered for output or averaging
!  (inout)       diag - structure used to regulate diagnostic output
!  (in)        static - If true, this is a static field that is always offered.
!  (in,opt)      mask - If present, use this real array as the data mask.

  type(diag_type), pointer :: diag => null()
  integer :: nz, i, j, k
  real, dimension(:,:,:), allocatable :: remapped_field
  logical :: staggered_in_x, staggered_in_y

  if (id_clock_diag_mediator>0) call cpu_clock_begin(id_clock_diag_mediator)

  ! Iterate over list of diag 'variants', e.g. CMOR aliases, different vertical
  ! grids, and post each.
  call assert(diag_field_id < diag_cs%next_free_diag_id, &
              'post_data_3d: Unregistered diagnostic id')
  diag => diag_cs%diags(diag_field_id)
  do while (associated(diag))
    call assert(associated(diag%axes), 'post_data_3d: axes is not associated')

    staggered_in_x = diag%axes%is_u_point .or. diag%axes%is_q_point
    staggered_in_y = diag%axes%is_v_point .or. diag%axes%is_q_point

    if (diag%v_extensive .and. .not.diag%axes%is_native) then
      ! The field is vertically integrated and needs to be re-gridded
      if (present(mask)) then
        call MOM_error(FATAL,"post_data_3d: no mask for regridded field.")
      endif

      if (id_clock_diag_remap>0) call cpu_clock_begin(id_clock_diag_remap)
      allocate(remapped_field(size(field,1), size(field,2), diag%axes%nz))
      call vertically_reintegrate_diag_field( &
              diag_cs%diag_remap_cs(diag%axes%vertical_coordinate_number), &
              diag_cs%G, diag_cs%h, staggered_in_x, staggered_in_y, &
              diag%mask3d, diag_cs%missing_value, field, remapped_field)
      if (id_clock_diag_remap>0) call cpu_clock_end(id_clock_diag_remap)
      if (associated(diag%mask3d)) then
        ! Since 3d masks do not vary in the vertical, just use as much as is
        ! needed.
        call post_data_3d_low(diag, remapped_field, diag_cs, is_static, &
                              mask=diag%mask3d(:,:,:diag%axes%nz))
      else
        call post_data_3d_low(diag, remapped_field, diag_cs, is_static)
      endif
      if (id_clock_diag_remap>0) call cpu_clock_begin(id_clock_diag_remap)
      deallocate(remapped_field)
      if (id_clock_diag_remap>0) call cpu_clock_end(id_clock_diag_remap)
    elseif (diag%axes%needs_remapping) then
      ! Remap this field to another vertical coordinate.
      if (present(mask)) then
        call MOM_error(FATAL,"post_data_3d: no mask for regridded field.")
      endif

      if (id_clock_diag_remap>0) call cpu_clock_begin(id_clock_diag_remap)
      allocate(remapped_field(size(field,1), size(field,2), diag%axes%nz))
      call diag_remap_do_remap(diag_cs%diag_remap_cs( &
              diag%axes%vertical_coordinate_number), &
              diag_cs%G, diag_cs%h, staggered_in_x, staggered_in_y, &
              diag%mask3d, diag_cs%missing_value, field, remapped_field)
      if (id_clock_diag_remap>0) call cpu_clock_end(id_clock_diag_remap)
      if (associated(diag%mask3d)) then
        ! Since 3d masks do not vary in the vertical, just use as much as is
        ! needed.
        call post_data_3d_low(diag, remapped_field, diag_cs, is_static, &
                              mask=diag%mask3d(:,:,:diag%axes%nz))
      else
        call post_data_3d_low(diag, remapped_field, diag_cs, is_static)
      endif
      if (id_clock_diag_remap>0) call cpu_clock_begin(id_clock_diag_remap)
      deallocate(remapped_field)
      if (id_clock_diag_remap>0) call cpu_clock_end(id_clock_diag_remap)
    elseif (diag%axes%needs_interpolating) then
      ! Interpolate this field to another vertical coordinate.
      if (present(mask)) then
        call MOM_error(FATAL,"post_data_3d: no mask for regridded field.")
      endif

      if (id_clock_diag_remap>0) call cpu_clock_begin(id_clock_diag_remap)
      allocate(remapped_field(size(field,1), size(field,2), diag%axes%nz+1))
      call vertically_interpolate_diag_field(diag_cs%diag_remap_cs( &
              diag%axes%vertical_coordinate_number), &
              diag_cs%G, diag_cs%h, staggered_in_x, staggered_in_y, &
              diag%mask3d, diag_cs%missing_value, field, remapped_field)
      if (id_clock_diag_remap>0) call cpu_clock_end(id_clock_diag_remap)
      if (associated(diag%mask3d)) then
        ! Since 3d masks do not vary in the vertical, just use as much as is
        ! needed.
        call post_data_3d_low(diag, remapped_field, diag_cs, is_static, &
                              mask=diag%mask3d(:,:,:diag%axes%nz+1))
      else
        call post_data_3d_low(diag, remapped_field, diag_cs, is_static)
      endif
      if (id_clock_diag_remap>0) call cpu_clock_begin(id_clock_diag_remap)
      deallocate(remapped_field)
      if (id_clock_diag_remap>0) call cpu_clock_end(id_clock_diag_remap)
    else
      call post_data_3d_low(diag, field, diag_cs, is_static, mask)
    endif
    diag => diag%next
  enddo
  if (id_clock_diag_mediator>0) call cpu_clock_end(id_clock_diag_mediator)

end subroutine post_data_3d

subroutine post_data_3d_low(diag, field, diag_cs, is_static, mask)
  type(diag_type),   intent(in) :: diag
  real,    target,   intent(in) :: field(:,:,:)
  type(diag_ctrl),   intent(in) :: diag_cs
  logical, optional, intent(in) :: is_static
  real,    optional, intent(in) :: mask(:,:,:)

! Arguments:
!  (in) diag          - the diagnostic to post.
!  (in)         field - 3-d array being offered for output or averaging
!  (inout)    diag_cs - structure used to regulate diagnostic output
!  (in)        static - If true, this is a static field that is always offered.
!  (in,opt)      mask - If present, use this real array as the data mask.

  real, dimension(:,:,:), pointer :: locfield => NULL()
  logical :: used  ! The return value of send_data is not used for anything.
  logical :: is_stat
  integer :: isv, iev, jsv, jev

  is_stat = .false. ; if (present(is_static)) is_stat = is_static

  ! Determine the proper array indices, noting that because of the (:,:)
  ! declaration of field, symmetric arrays are using a SW-grid indexing,
  ! but non-symmetric arrays are using a NE-grid indexing.  Send_data
  ! actually only uses the difference between ie and is to determine
  ! the output data size and assumes that halos are symmetric.
  isv = diag_cs%is ; iev = diag_cs%ie ; jsv = diag_cs%js ; jev = diag_cs%je

  if ( size(field,1) == diag_cs%ied-diag_cs%isd +1 ) then
    isv = diag_cs%is ; iev = diag_cs%ie        ! Data domain
  elseif ( size(field,1) == diag_cs%ied-diag_cs%isd +2 ) then
    isv = diag_cs%is ; iev = diag_cs%ie+1      ! Symmetric data domain
  elseif ( size(field,1) == diag_cs%ie-diag_cs%is +1 ) then
    isv = 1 ; iev = diag_cs%ie + 1-diag_cs%is  ! Computational domain
  elseif ( size(field,1) == diag_cs%ie-diag_cs%is +2 ) then
    isv = 1 ; iev = diag_cs%ie + 2-diag_cs%is  ! Symmetric computational domain
  else
    call MOM_error(FATAL,"post_data_3d_low: peculiar size in i-direction")
  endif
  if ( size(field,2) == diag_cs%jed-diag_cs%jsd +1 ) then
    jsv = diag_cs%js ; jev = diag_cs%je        ! Data domain
  elseif ( size(field,2) == diag_cs%jed-diag_cs%jsd +2 ) then
    jsv = diag_cs%js ; jev = diag_cs%je+1      ! Symmetric data domain
  elseif ( size(field,2) == diag_cs%je-diag_cs%js +1 ) then
    jsv = 1 ; jev = diag_cs%je + 1-diag_cs%js  ! Computational domain
  elseif ( size(field,1) == diag_cs%je-diag_cs%js +2 ) then
    jsv = 1 ; jev = diag_cs%je + 2-diag_cs%js  ! Symmetric computational domain
  else
    call MOM_error(FATAL,"post_data_3d_low: peculiar size in j-direction")
  endif

  if (diag%conversion_factor/=0.) then
    allocate( locfield( lbound(field,1):ubound(field,1), lbound(field,2):ubound(field,2), &
                        lbound(field,3):ubound(field,3) ) )
    locfield(isv:iev,jsv:jev,:) = field(isv:iev,jsv:jev,:) * diag%conversion_factor
  else
    locfield => field
  endif

  if (diag%fms_diag_id>0) then
    if (is_stat) then
      if (present(mask)) then
        call assert(size(locfield) == size(mask), &
            'post_data_3d_low is_stat: mask size mismatch: '//diag%debug_str)
        used = send_data(diag%fms_diag_id, locfield, &
                         is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, rmask=mask)
     !elseif(associated(diag%mask3d)) then
     !  used = send_data(diag_field_id, locfield, &
     !                   is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, rmask=diag%mask3d)
      else
        used = send_data(diag%fms_diag_id, locfield, &
                         is_in=isv, js_in=jsv, ie_in=iev, je_in=jev)
      endif
    elseif (diag_cs%ave_enabled) then
      if (present(mask)) then
        call assert(size(locfield) == size(mask), &
            'post_data_3d_low: mask size mismatch: '//diag%debug_str)
        used = send_data(diag%fms_diag_id, locfield, diag_cs%time_end, &
                         is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                         weight=diag_cs%time_int, rmask=mask)
      elseif(associated(diag%mask3d)) then
        call assert(size(locfield) == size(diag%mask3d), &
            'post_data_3d_low: mask3d size mismatch: '//diag%debug_str)
        used = send_data(diag%fms_diag_id, locfield, diag_cs%time_end, &
                         is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                         weight=diag_cs%time_int, rmask=diag%mask3d)
      else
        used = send_data(diag%fms_diag_id, locfield, diag_cs%time_end, &
                         is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                         weight=diag_cs%time_int)
      endif
    endif
  endif
  if (diag%fms_xyave_diag_id>0) then
    call post_xy_average(diag_cs, diag, locfield)
  endif
  if (diag%conversion_factor/=0.) deallocate( locfield )

end subroutine post_data_3d_low

!> Post the horizontally area-averaged diagnostic
subroutine post_xy_average(diag_cs, diag, field)
  type(diag_type),   intent(in) :: diag !< This diagnostic
  real,    target,   intent(in) :: field(:,:,:) !< Diagnostic field
  type(diag_ctrl),   intent(in) :: diag_cs !< Diagnostics mediator control structure
  ! Local variable
  real, dimension(size(field,3)) :: averaged_field
  logical :: staggered_in_x, staggered_in_y, used
  integer :: nz, remap_nz, coord

  if (.not. diag_cs%ave_enabled) then
    return
  endif

  staggered_in_x = diag%axes%is_u_point .or. diag%axes%is_q_point
  staggered_in_y = diag%axes%is_v_point .or. diag%axes%is_q_point

  if (diag%axes%is_native) then
    call horizontally_average_diag_field(diag_cs%G, diag_cs%h, &
                                         staggered_in_x, staggered_in_y, &
                                         diag%axes%is_layer, diag%v_extensive, &
                                         diag_cs%missing_value, field, averaged_field)
  else
    nz = size(field, 3)
    coord = diag%axes%vertical_coordinate_number
    remap_nz = diag_cs%diag_remap_cs(coord)%nz

    call assert(diag_cs%diag_remap_cs(coord)%initialized, &
                'post_xy_average: remap_cs not initialized.')

    call assert(IMPLIES(diag%axes%is_layer, nz == remap_nz), &
              'post_xy_average: layer field dimension mismatch.')
    call assert(IMPLIES(.not. diag%axes%is_layer, nz == remap_nz+1), &
              'post_xy_average: interface field dimension mismatch.')

    call horizontally_average_diag_field(diag_cs%G, diag_cs%diag_remap_cs(coord)%h, &
                                         staggered_in_x, staggered_in_y, &
                                         diag%axes%is_layer, diag%v_extensive, &
                                         diag_cs%missing_value, field, averaged_field)
  endif

  used = send_data(diag%fms_xyave_diag_id, averaged_field, diag_cs%time_end, &
                   weight=diag_cs%time_int)
end subroutine post_xy_average

subroutine enable_averaging(time_int_in, time_end_in, diag_cs)
  real, intent(in) :: time_int_in
  type(time_type), intent(in) :: time_end_in
  type(diag_ctrl), intent(inout) :: diag_cs

! This subroutine enables the accumulation of time averages over the
! specified time interval.

! Arguments:
!  (in)      time_int_in - time interval in s over which any
!                          values that are offered are valid.
!  (in)      time_end_in - end time in s of the valid interval
!  (inout)   diag        - structure used to regulate diagnostic output

!  if (num_file==0) return
  diag_cs%time_int = time_int_in
  diag_cs%time_end = time_end_in
  diag_cs%ave_enabled = .true.
end subroutine enable_averaging

! Call this subroutine to avoid averaging any offered fields.
subroutine disable_averaging(diag_cs)
  type(diag_ctrl), intent(inout) :: diag_cs

! Argument:
! diag - structure used to regulate diagnostic output

  diag_cs%time_int = 0.0
  diag_cs%ave_enabled = .false.

end subroutine disable_averaging

! Call this subroutine to determine whether the averaging is
! currently enabled.  .true. is returned if it is.
function query_averaging_enabled(diag_cs, time_int, time_end)
  type(diag_ctrl),           intent(in) :: diag_cs
  real,            optional, intent(out) :: time_int
  type(time_type), optional, intent(out) :: time_end
  logical :: query_averaging_enabled

! Arguments:
!  (in)          diag - structure used to regulate diagnostic output
!  (out,opt) time_int - current setting of diag%time_int, in s
!  (out,opt) time_end - current setting of diag%time_end

  if (present(time_int)) time_int = diag_cs%time_int
  if (present(time_end)) time_end = diag_cs%time_end
  query_averaging_enabled = diag_cs%ave_enabled
end function query_averaging_enabled

function get_diag_time_end(diag_cs)
  type(diag_ctrl), intent(in)  :: diag_cs
  type(time_type) :: get_diag_time_end

! Argument:
! (in) diag - structure used to regulate diagnostic output

!   This function returns the valid end time for diagnostics that are handled
! outside of the MOM6 infrastructure, such as via the generic tracer code.

  get_diag_time_end = diag_cs%time_end
end function get_diag_time_end

!> Returns the "diag_mediator" handle for a group (native, CMOR, z-coord, ...) of diagnostics derived from one field.
integer function register_diag_field(module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, mask_variant, standard_name,      &
     verbose, do_not_log, err_msg, interp_method, tile_count, cmor_field_name, &
     cmor_long_name, cmor_units, cmor_standard_name, cell_methods, &
     x_cell_method, y_cell_method, v_cell_method, conversion, v_extensive)
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ocean_model" or "ice_shelf_model"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp), target, intent(in) :: axes !< Container w/ up to 3 integer handles that indicates axes for this field
  type(time_type),  intent(in) :: init_time !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values.
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
  logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided with post_data calls (not used in MOM?)
  logical,          optional, intent(in) :: verbose !< If true, FMS is verbose (not used in MOM?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be placed (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count !< no clue (not used in MOM?)
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
  character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
  character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
  character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
  character(len=*), optional, intent(in) :: cell_methods !< String to append as cell_methods attribute. Use '' to have no attribute.
                                                         !! If present, this overrides the default constructed from the default for
                                                         !! each individual axis direction.
  character(len=*), optional, intent(in) :: x_cell_method !< Specifies the cell method for the x-direction. Use '' have no method.
  character(len=*), optional, intent(in) :: y_cell_method !< Specifies the cell method for the y-direction. Use '' have no method.
  character(len=*), optional, intent(in) :: v_cell_method !< Specifies the cell method for the vertical direction. Use '' have no method.
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to file
  logical,          optional, intent(in) :: v_extensive !< True for vertically extensive fields (vertically integrated). Default/absent for intensive.
  ! Local variables
  real :: MOM_missing_value
  type(diag_ctrl), pointer :: diag_cs
  type(axes_grp), pointer :: remap_axes => null()
  integer :: dm_id, i
  character(len=256) :: new_module_name
  logical :: active

  MOM_missing_value = axes%diag_cs%missing_value
  if(present(missing_value)) MOM_missing_value = missing_value

  diag_cs => axes%diag_cs
  dm_id = -1

  ! Register the native diagnostic
  active = register_diag_field_expand_cmor(dm_id, module_name, field_name, axes, &
             init_time, long_name=long_name, units=units, missing_value=MOM_missing_value, &
             range=range, mask_variant=mask_variant, standard_name=standard_name, &
             verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
             interp_method=interp_method, tile_count=tile_count, &
             cmor_field_name=cmor_field_name, cmor_long_name=cmor_long_name, &
             cmor_units=cmor_units, cmor_standard_name=cmor_standard_name, &
             cell_methods=cell_methods, x_cell_method=x_cell_method, &
             y_cell_method=y_cell_method, v_cell_method=v_cell_method, &
             conversion=conversion, v_extensive=v_extensive)

  ! For each diagnostic coordinate register the diagnostic again under a different module name
  do i=1,diag_cs%num_diag_coords
    new_module_name = trim(module_name)//'_'//trim(diag_cs%diag_remap_cs(i)%diag_module_suffix)

    ! Register diagnostics remapped to z vertical coordinate
    if (axes%rank == 3) then
      remap_axes => null()
      if ((axes%id .eq. diag_cs%axesTL%id)) then
          remap_axes => diag_cs%remap_axesTL(i)
      elseif(axes%id .eq. diag_cs%axesBL%id) then
          remap_axes => diag_cs%remap_axesBL(i)
      elseif(axes%id .eq. diag_cs%axesCuL%id ) then
          remap_axes => diag_cs%remap_axesCuL(i)
      elseif(axes%id .eq. diag_cs%axesCvL%id) then
          remap_axes => diag_cs%remap_axesCvL(i)
      elseif(axes%id .eq. diag_cs%axesTi%id) then
          remap_axes => diag_cs%remap_axesTi(i)
      elseif(axes%id .eq. diag_cs%axesBi%id) then
          remap_axes => diag_cs%remap_axesBi(i)
      elseif(axes%id .eq. diag_cs%axesCui%id ) then
          remap_axes => diag_cs%remap_axesCui(i)
      elseif(axes%id .eq. diag_cs%axesCvi%id) then
          remap_axes => diag_cs%remap_axesCvi(i)
      endif
      ! When the MOM_diag_to_Z module has been obsoleted we can assume remap_axes will
      ! always exist but in the mean-time we have to do this check:
      ! call assert(associated(remap_axes), 'register_diag_field: remap_axes not set')
      if (associated(remap_axes)) then
        if (remap_axes%needs_remapping .or. remap_axes%needs_interpolating) then
          active = register_diag_field_expand_cmor(dm_id, new_module_name, field_name, remap_axes, &
                     init_time, long_name=long_name, units=units, missing_value=MOM_missing_value, &
                     range=range, mask_variant=mask_variant, standard_name=standard_name, &
                     verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                     interp_method=interp_method, tile_count=tile_count, &
                     cmor_field_name=cmor_field_name, cmor_long_name=cmor_long_name, &
                     cmor_units=cmor_units, cmor_standard_name=cmor_standard_name, &
                     cell_methods=cell_methods, x_cell_method=x_cell_method, &
                     y_cell_method=y_cell_method, v_cell_method=v_cell_method, &
                     conversion=conversion, v_extensive=v_extensive)
          if (active) then
            call diag_remap_set_active(diag_cs%diag_remap_cs(i))
          endif
        endif ! remap_axes%needs_remapping
      endif ! associated(remap_axes)
    endif ! axes%rank == 3
  enddo ! i

  register_diag_field = dm_id

end function register_diag_field

!> Returns True if either the native of CMOr version of the diagnostic were registered. Updates 'dm_id'
!! after calling register_diag_field_expand_axes() for both native and CMOR variants of the field.
logical function register_diag_field_expand_cmor(dm_id, module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, mask_variant, standard_name,      &
     verbose, do_not_log, err_msg, interp_method, tile_count, cmor_field_name, &
     cmor_long_name, cmor_units, cmor_standard_name, cell_methods, &
     x_cell_method, y_cell_method, v_cell_method, conversion, v_extensive)
  integer,          intent(inout) :: dm_id !< The diag_mediator ID for this diagnostic group
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ocean_model" or "ice_shelf_model"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp), target, intent(in) :: axes !< Container w/ up to 3 integer handles that indicates axes for this field
  type(time_type),  intent(in) :: init_time !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values.
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
  logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided with post_data calls (not used in MOM?)
  logical,          optional, intent(in) :: verbose !< If true, FMS is verbose (not used in MOM?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be placed (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count !< no clue (not used in MOM?)
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
  character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
  character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
  character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
  character(len=*), optional, intent(in) :: cell_methods !< String to append as cell_methods attribute. Use '' to have no attribute.
                                                         !! If present, this overrides the default constructed from the default for
                                                         !! each individual axis direction.
  character(len=*), optional, intent(in) :: x_cell_method !< Specifies the cell method for the x-direction. Use '' have no method.
  character(len=*), optional, intent(in) :: y_cell_method !< Specifies the cell method for the y-direction. Use '' have no method.
  character(len=*), optional, intent(in) :: v_cell_method !< Specifies the cell method for the vertical direction. Use '' have no method.
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to file
  logical,          optional, intent(in) :: v_extensive !< True for vertically extensive fields (vertically integrated). Default/absent for intensive.
  ! Local variables
  real :: MOM_missing_value
  type(diag_ctrl), pointer :: diag_cs
  type(diag_type), pointer :: this_diag => null()
  integer :: fms_id, fms_xyave_id
  character(len=256) :: posted_cmor_units, posted_cmor_standard_name, posted_cmor_long_name, cm_string, msg

  MOM_missing_value = axes%diag_cs%missing_value
  if(present(missing_value)) MOM_missing_value = missing_value

  register_diag_field_expand_cmor = .false.
  diag_cs => axes%diag_cs

  ! Set up the 'primary' diagnostic, first get an underlying FMS id
  fms_id = register_diag_field_expand_axes(module_name, field_name, axes, init_time, &
             long_name=long_name, units=units, missing_value=MOM_missing_value, &
             range=range, mask_variant=mask_variant, standard_name=standard_name, &
             verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
             interp_method=interp_method, tile_count=tile_count)
  call attach_cell_methods(fms_id, axes, cm_string, &
                           cell_methods, x_cell_method, y_cell_method, v_cell_method, &
                           v_extensive=v_extensive)
  if (is_root_pe() .and. diag_CS%doc_unit > 0) then
    msg = ''
    if (present(cmor_field_name)) msg = 'CMOR equivalent is "'//trim(cmor_field_name)//'"'
    call log_available_diag(fms_id>0, module_name, field_name, cm_string, &
                            msg, diag_CS, long_name, units, standard_name)
  endif
  ! Associated horizontally area-averaged diagnostic
  fms_xyave_id = DIAG_FIELD_NOT_FOUND
  if (associated(axes%xyave_axes)) then
    fms_xyave_id = register_diag_field_expand_axes(module_name, trim(field_name)//'_xyave', &
             axes%xyave_axes, init_time, &
             long_name=long_name, units=units, missing_value=MOM_missing_value, &
             range=range, mask_variant=mask_variant, standard_name=standard_name, &
             verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
             interp_method=interp_method, tile_count=tile_count)
    call attach_cell_methods(fms_xyave_id, axes%xyave_axes, cm_string, &
                             cell_methods, v_cell_method, v_extensive=v_extensive)
    if (is_root_pe() .and. diag_CS%doc_unit > 0) then
      msg = ''
      if (present(cmor_field_name)) msg = 'CMOR equivalent is "'//trim(cmor_field_name)//'_xyave"'
      call log_available_diag(fms_xyave_id>0, module_name, trim(field_name)//'_xyave', cm_string, &
                              msg, diag_CS, long_name, units, standard_name)
    endif
  endif
  this_diag => null()
  if (fms_id /= DIAG_FIELD_NOT_FOUND .or. fms_xyave_id /= DIAG_FIELD_NOT_FOUND) then
    call add_diag_to_list(diag_cs, dm_id, fms_id, this_diag, axes, module_name, field_name, msg)
    this_diag%fms_xyave_diag_id = fms_xyave_id

    if (present(v_extensive)) this_diag%v_extensive = v_extensive
    if (present(conversion)) this_diag%conversion_factor = conversion
    register_diag_field_expand_cmor = .true.
  endif

  ! For the CMOR variation of the above diagnostic
  if (present(cmor_field_name)) then
    ! Fallback values for strings set to "NULL"
    posted_cmor_units = "not provided"           !
    posted_cmor_standard_name = "not provided"   ! Values might be able to be replaced with a CS%missing field?
    posted_cmor_long_name = "not provided"       !

    ! If attributes are present for MOM variable names, use them first for the register_diag_field
    ! call for CMOR verison of the variable
    if (present(units)) posted_cmor_units = units
    if (present(standard_name)) posted_cmor_standard_name = standard_name
    if (present(long_name)) posted_cmor_long_name = long_name

    ! If specified in the call to register_diag_field, override attributes with the CMOR versions
    if (present(cmor_units)) posted_cmor_units = cmor_units
    if (present(cmor_standard_name)) posted_cmor_standard_name = cmor_standard_name
    if (present(cmor_long_name)) posted_cmor_long_name = cmor_long_name

    fms_id = register_diag_field_expand_axes(module_name, cmor_field_name, axes, init_time,    &
               long_name=trim(posted_cmor_long_name), units=trim(posted_cmor_units),                    &
               missing_value=MOM_missing_value, range=range, mask_variant=mask_variant,                 &
               standard_name=trim(posted_cmor_standard_name), verbose=verbose, do_not_log=do_not_log,   &
               err_msg=err_msg, interp_method=interp_method, tile_count=tile_count)
    call attach_cell_methods(fms_id, axes, cm_string, &
                             cell_methods, x_cell_method, y_cell_method, v_cell_method, &
                             v_extensive=v_extensive)
    if (is_root_pe() .and. diag_CS%doc_unit > 0) then
      msg = 'native name is "'//trim(field_name)//'"'
      call log_available_diag(fms_id>0, module_name, cmor_field_name, cm_string, &
                              msg, diag_CS, posted_cmor_long_name, posted_cmor_units, &
                              posted_cmor_standard_name)
    endif
    ! Associated horizontally area-averaged diagnostic
    fms_xyave_id = DIAG_FIELD_NOT_FOUND
    if (associated(axes%xyave_axes)) then
      fms_xyave_id = register_diag_field_expand_axes(module_name, trim(cmor_field_name)//'_xyave', &
               axes%xyave_axes, init_time, &
               long_name=trim(posted_cmor_long_name), units=trim(posted_cmor_units),                    &
               missing_value=MOM_missing_value, range=range, mask_variant=mask_variant,                 &
               standard_name=trim(posted_cmor_standard_name), verbose=verbose, do_not_log=do_not_log,   &
               err_msg=err_msg, interp_method=interp_method, tile_count=tile_count)
      call attach_cell_methods(fms_xyave_id, axes%xyave_axes, cm_string, &
                               cell_methods, v_cell_method, v_extensive=v_extensive)
      if (is_root_pe() .and. diag_CS%doc_unit > 0) then
        msg = 'native name is "'//trim(field_name)//'_xyave"'
        call log_available_diag(fms_xyave_id>0, module_name, trim(cmor_field_name)//'_xyave', cm_string, &
                                msg, diag_CS, posted_cmor_long_name, posted_cmor_units, &
                                posted_cmor_standard_name)
      endif
    endif
    this_diag => null()
    if (fms_id /= DIAG_FIELD_NOT_FOUND .or. fms_xyave_id /= DIAG_FIELD_NOT_FOUND) then
      call add_diag_to_list(diag_cs, dm_id, fms_id, this_diag, axes, module_name, field_name, msg)
      this_diag%fms_xyave_diag_id = fms_xyave_id

      if (present(v_extensive)) this_diag%v_extensive = v_extensive
      if (present(conversion)) this_diag%conversion_factor = conversion
      register_diag_field_expand_cmor = .true.
    endif
  endif

end function register_diag_field_expand_cmor

!> Returns an FMS id from register_diag_field_fms (the diag_manager routine) after expanding axes (axes-group)
!! into handles and conditionally adding an FMS area_id for cell_measures.
integer function register_diag_field_expand_axes(module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, mask_variant, standard_name,  &
     verbose, do_not_log, err_msg, interp_method, tile_count)
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ocean_model" or "ice_shelf_model"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp), target, intent(in) :: axes !< Container w/ up to 3 integer handles that indicates axes for this field
  type(time_type),  intent(in) :: init_time !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values.
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
  logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided with post_data calls (not used in MOM?)
  logical,          optional, intent(in) :: verbose !< If true, FMS is verbose (not used in MOM?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be placed (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count !< no clue (not used in MOM?)
  ! Local variables
  integer :: fms_id, area_id

  ! This gets the cell area associated with the grid location of this variable
  area_id = axes%id_area

  ! Get the FMS diagnostic id
  if (present(interp_method) .or. axes%is_h_point) then
    ! If interp_method is provided we must use it
    if (area_id>0) then
      fms_id = register_diag_field_fms(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method=interp_method, tile_count=tile_count, area=area_id)
    else
      fms_id = register_diag_field_fms(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method=interp_method, tile_count=tile_count)
    endif
  else
    ! If interp_method is not provided and the field is not at an h-point then interp_method='none'
    if (area_id>0) then
      fms_id = register_diag_field_fms(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method='none', tile_count=tile_count, area=area_id)
    else
      fms_id = register_diag_field_fms(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method='none', tile_count=tile_count)
    endif
  endif

  register_diag_field_expand_axes = fms_id

end function register_diag_field_expand_axes

!> Create a diagnostic type and attached to list
subroutine add_diag_to_list(diag_cs, dm_id, fms_id, this_diag, axes, module_name, field_name, msg)
  type(diag_ctrl),        pointer       :: diag_cs !< Diagnostics mediator control structure
  integer,                intent(inout) :: dm_id !< The diag_mediator ID for this diagnostic group
  integer,                intent(in)    :: fms_id !< The FMS diag_manager ID for this diagnostic
  type(diag_type),        pointer       :: this_diag !< This diagnostic
  type(axes_grp), target, intent(in)    :: axes !< Container w/ up to 3 integer handles that indicates axes for this field
  character(len=*),       intent(in)    :: module_name !< Name of this module, usually "ocean_model" or "ice_shelf_model"
  character(len=*),       intent(in)    :: field_name !< Name of diagnostic
  character(len=*),       intent(in)    :: msg !< Message for errors

  ! If the diagnostic is needed obtain a diag_mediator ID (if needed)
  if (dm_id == -1) dm_id = get_new_diag_id(diag_cs)
  ! Create a new diag_type to store links in
  call alloc_diag_with_id(dm_id, diag_cs, this_diag)
  call assert(associated(this_diag), trim(msg)//': diag_type allocation failed')
  ! Record FMS id, masks and conversion factor, in diag_type
  this_diag%fms_diag_id = fms_id
  this_diag%debug_str = trim(module_name)//"-"//trim(field_name)
  call set_diag_mask(this_diag, diag_cs, axes)
  this_diag%axes => axes

end subroutine add_diag_to_list

!> Attaches "cell_methods" attribute to a variable based on defaults for axes_grp or optional arguments.
subroutine attach_cell_methods(id, axes, ostring, cell_methods, &
                               x_cell_method, y_cell_method, v_cell_method, v_extensive)
  integer,                    intent(in)  :: id !< Handle to diagnostic
  type(axes_grp),             intent(in)  :: axes !< Container w/ up to 3 integer handles that indicates axes for this field
  character(len=*),           intent(out) :: ostring !< The cell_methods strings that would appear in the file
  character(len=*), optional, intent(in)  :: cell_methods !< String to append as cell_methods attribute. Use '' to have no attribute.
                                                         !! If present, this overrides the default constructed from the default for
                                                         !! each individual axis direction.
  character(len=*), optional, intent(in)  :: x_cell_method !< Specifies the cell method for the x-direction. Use '' have no method.
  character(len=*), optional, intent(in)  :: y_cell_method !< Specifies the cell method for the y-direction. Use '' have no method.
  character(len=*), optional, intent(in)  :: v_cell_method !< Specifies the cell method for the vertical direction. Use '' have no method.
  logical,          optional, intent(in)  :: v_extensive !< True for vertically extensive fields (vertically integrated). Default/absent for intensive.
  ! Local variables
  character(len=9) :: axis_name

  ostring = ''
  if (present(cell_methods)) then
    if (present(x_cell_method) .or. present(y_cell_method) .or. present(v_cell_method) &
        .or. present(v_extensive)) then
      call MOM_error(FATAL, "attach_cell_methods: " // &
           'Individual direction cell method was specified along with a "cell_methods" string.')
    endif
    if (len(trim(cell_methods))>0) then
      call diag_field_add_attribute(id, 'cell_methods', trim(cell_methods))
      ostring = trim(cell_methods)
    endif
  else
    if (present(x_cell_method)) then
      if (len(trim(x_cell_method))>0) then
        call get_diag_axis_name(axes%handles(1), axis_name)
        call diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(x_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(x_cell_method)
      endif
    else
      if (len(trim(axes%x_cell_method))>0) then
        call get_diag_axis_name(axes%handles(1), axis_name)
        call diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(axes%x_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(axes%x_cell_method)
      endif
    endif
    if (present(y_cell_method)) then
      if (len(trim(y_cell_method))>0) then
        call get_diag_axis_name(axes%handles(2), axis_name)
        call diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(y_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(y_cell_method)
      endif
    else
      if (len(trim(axes%y_cell_method))>0) then
        call get_diag_axis_name(axes%handles(2), axis_name)
        call diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(axes%y_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(axes%y_cell_method)
      endif
    endif
    if (present(v_cell_method)) then
      if (present(v_extensive)) call MOM_error(FATAL, "attach_cell_methods: " // &
           'Vertical cell method was specified along with the vertically extensive flag.')
      if (len(trim(v_cell_method))>0) then
        if (axes%rank==1) then
          call get_diag_axis_name(axes%handles(1), axis_name)
        elseif (axes%rank==3) then
          call get_diag_axis_name(axes%handles(3), axis_name)
        endif
        call diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(v_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(v_cell_method)
      endif
    elseif (present(v_extensive)) then
      if (axes%rank==1) then
        call get_diag_axis_name(axes%handles(1), axis_name)
      elseif (axes%rank==3) then
        call get_diag_axis_name(axes%handles(3), axis_name)
      endif
      call diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':sum')
      ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':sum'
    else
      if (len(trim(axes%v_cell_method))>0) then
        if (axes%rank==1) then
          call get_diag_axis_name(axes%handles(1), axis_name)
        elseif (axes%rank==3) then
          call get_diag_axis_name(axes%handles(3), axis_name)
        endif
        call diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(axes%v_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(axes%v_cell_method)
      endif
    endif
  endif
  ostring = adjustl(ostring)
end subroutine attach_cell_methods

function register_scalar_field(module_name, field_name, init_time, diag_cs, &
     long_name, units, missing_value, range, standard_name, &
     do_not_log, err_msg, interp_method, cmor_field_name, &
     cmor_long_name, cmor_units, cmor_standard_name)
  integer :: register_scalar_field
  character(len=*), intent(in) :: module_name, field_name
  type(time_type),  intent(in) :: init_time
  type(diag_ctrl),  intent(inout) :: diag_cs
  character(len=*), optional, intent(in) :: long_name, units, standard_name
  real,             optional, intent(in) :: missing_value, range(2)
  logical,          optional, intent(in) :: do_not_log
  character(len=*), optional, intent(out):: err_msg
  character(len=*), optional, intent(in) :: interp_method
  character(len=*), optional, intent(in) :: cmor_field_name, cmor_long_name
  character(len=*), optional, intent(in) :: cmor_units, cmor_standard_name

  ! Output:    An integer handle for a diagnostic array.
  ! Arguments:
  !  (in)      module_name   - name of this module, usually "ocean_model" or "ice_shelf_model".
  !  (in)      field_name    - name of the diagnostic field.
  !  (in)      init_time     - time at which a field is first available?
  !  (inout)   diag_cs     - structure used to regulate diagnostic output
  !  (in,opt)  long_name     - long name of a field
  !  (in,opt)  units         - units of a field
  !  (in,opt)  missing_value - indicates missing values
  !  (in,opt)  standard_name - standardized name associated with a field

  ! Following params have yet to be used in MOM.
  !  (in,opt)  range         - valid range of a variable
  !  (in,opt)  verbose       - If true, FMS is verbosed
  !  (in,opt)  do_not_log    - If true, do not log something
  !  (out,opt) err_msg       - character string into which an error message might be placed
  !  (in,opt)  interp_method - If 'none' indicates the field should not be interpolated as a scalar
  !  (in,opt)  tile_count    - no clue

  real :: MOM_missing_value
  integer :: dm_id, fms_id
  type(diag_type), pointer :: diag => null(), cmor_diag => null()
  character(len=256) :: posted_cmor_units, posted_cmor_standard_name, posted_cmor_long_name

  MOM_missing_value = diag_cs%missing_value
  if(present(missing_value)) MOM_missing_value = missing_value

  dm_id = -1
  diag => null()
  cmor_diag => null()

  fms_id = register_diag_field_fms(module_name, field_name, init_time, &
      long_name=long_name, units=units, missing_value=MOM_missing_value, &
      range=range, standard_name=standard_name, do_not_log=do_not_log, err_msg=err_msg)
  if (fms_id /= DIAG_FIELD_NOT_FOUND) then
    dm_id = get_new_diag_id(diag_cs)
    call alloc_diag_with_id(dm_id, diag_cs, diag)
    call assert(associated(diag), 'register_scalar_field: diag allocation failed')
    diag%fms_diag_id = fms_id
    diag%debug_str = trim(module_name)//"-"//trim(field_name)
  endif

  if (present(cmor_field_name)) then
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

    fms_id = register_diag_field_fms(module_name, cmor_field_name, init_time, &
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
    endif
  endif

  ! Document diagnostics in list of available diagnostics
  if (is_root_pe() .and. diag_CS%doc_unit > 0) then
    call log_available_diag(associated(diag), module_name, field_name, '', '', diag_CS, &
                            long_name, units, standard_name)
    if (present(cmor_field_name)) then
      call log_available_diag(associated(cmor_diag), module_name, cmor_field_name, &
                              '', '', diag_CS, posted_cmor_long_name, posted_cmor_units, &
                              posted_cmor_standard_name)
    endif
  endif

  register_scalar_field = dm_id

end function register_scalar_field

!> Registers a static diagnostic, returning an integer handle
function register_static_field(module_name, field_name, axes, &
     long_name, units, missing_value, range, mask_variant, standard_name, &
     do_not_log, interp_method, tile_count, &
     cmor_field_name, cmor_long_name, cmor_units, cmor_standard_name, area)
  integer :: register_static_field
  character(len=*), intent(in) :: module_name, field_name
  type(axes_grp),   intent(in) :: axes
  character(len=*), optional, intent(in) :: long_name, units, standard_name
  real,             optional, intent(in) :: missing_value, range(2)
  logical,          optional, intent(in) :: mask_variant, do_not_log
  character(len=*), optional, intent(in) :: interp_method
  integer,          optional, intent(in) :: tile_count
  character(len=*), optional, intent(in) :: cmor_field_name, cmor_long_name
  character(len=*), optional, intent(in) :: cmor_units, cmor_standard_name
  integer,          optional, intent(in) :: area !< fms_id for area_t

  ! Output:    An integer handle for a diagnostic array.
  ! Arguments:
  !  (in)      module_name   - name of this module, usually "ocean_model" or "ice_shelf_model".
  !  (in)      field_name    - name of the diagnostic field
  !  (in)      axes          - container with up to 3 integer handles that indicates axes for this field
  !  (in,opt)  long_name     - long name of a field
  !  (in,opt)  units         - units of a field
  !  (in,opt)  missing_value - A value that indicates missing values.
  !  (in,opt)  standard_name - standardized name associated with a field

  ! Following params have yet to be used in MOM.
  !  (in,opt)  range          - valid range of a variable
  !  (in,opt)  mask_variant   - If true a logical mask must be provided with post_data calls
  !  (in,opt)  do_not_log     - If true, do not log something
  !  (in,opt)  interp_method  - If 'none' indicates the field should not be interpolated as a scalar
  !  (in,opt)  tile_count     - no clue

  real :: MOM_missing_value
  type(diag_ctrl), pointer :: diag_cs
  type(diag_type), pointer :: diag => null(), cmor_diag => null()
  integer :: dm_id, fms_id, cmor_id
  character(len=256) :: posted_cmor_units, posted_cmor_standard_name, posted_cmor_long_name

  MOM_missing_value = axes%diag_cs%missing_value
  if(present(missing_value)) MOM_missing_value = missing_value

  diag_cs => axes%diag_cs
  dm_id = -1
  diag => null()
  cmor_diag => null()

  fms_id = register_static_field_fms(module_name, field_name, axes%handles, &
         long_name=long_name, units=units, missing_value=MOM_missing_value, &
         range=range, mask_variant=mask_variant, standard_name=standard_name, &
         do_not_log=do_not_log, &
         interp_method=interp_method, tile_count=tile_count, area=area)
  if (fms_id /= DIAG_FIELD_NOT_FOUND) then
    dm_id = get_new_diag_id(diag_cs)
    call alloc_diag_with_id(dm_id, diag_cs, diag)
    call assert(associated(diag), 'register_static_field: diag allocation failed')
    diag%fms_diag_id = fms_id
    diag%debug_str = trim(module_name)//"-"//trim(field_name)
  endif

  if (present(cmor_field_name)) then
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

    fms_id = register_static_field_fms(module_name, cmor_field_name, &
      axes%handles, long_name=trim(posted_cmor_long_name), units=trim(posted_cmor_units), &
      missing_value=MOM_missing_value, range=range, mask_variant=mask_variant,            &
      standard_name=trim(posted_cmor_standard_name), do_not_log=do_not_log,               &
      interp_method=interp_method, tile_count=tile_count, area=area)
    if (fms_id /= DIAG_FIELD_NOT_FOUND) then
      if (dm_id == -1) then
        dm_id = get_new_diag_id(diag_cs)
      endif
      call alloc_diag_with_id(dm_id, diag_cs, cmor_diag)
      cmor_diag%fms_diag_id = fms_id
      cmor_diag%debug_str = trim(module_name)//"-"//trim(cmor_field_name)
    endif
  endif

  ! Document diagnostics in list of available diagnostics
  if (is_root_pe() .and. diag_CS%doc_unit > 0) then
    call log_available_diag(associated(diag), module_name, field_name, '', '', diag_CS, &
                            long_name, units, standard_name)
    if (present(cmor_field_name)) then
      call log_available_diag(associated(cmor_diag), module_name, cmor_field_name, &
                              '', '', diag_CS, posted_cmor_long_name, posted_cmor_units, &
                              posted_cmor_standard_name)
    endif
  endif

  register_static_field = dm_id

end function register_static_field

subroutine describe_option(opt_name, value, diag_CS)
  character(len=*), intent(in) :: opt_name, value
  type(diag_ctrl), intent(in) :: diag_CS

  character(len=240) :: mesg
  integer :: len_ind

  len_ind = len_trim(value)  ! Add error handling for long values?

  mesg = "    ! "//trim(opt_name)//": "//trim(value)
  write(diag_CS%doc_unit, '(a)') trim(mesg)
end subroutine describe_option

!> Registers a diagnostic using the information encapsulated in the vardesc
!! type argument and returns an integer handle to this diagostic.  That
!! integer handle is negative if the diagnostic is unused.
function ocean_register_diag(var_desc, G, diag_CS, day)
  integer :: ocean_register_diag  !< An integer handle to this diagnostic.
  type(vardesc),         intent(in) :: var_desc !< The vardesc type describing the diagnostic
  type(ocean_grid_type), intent(in) :: G        !< The ocean's grid type
  type(diag_ctrl), intent(in), target :: diag_CS  !< The diagnotic control structure
  type(time_type),       intent(in) :: day      !< The current model time

  character(len=64) :: var_name         ! A variable's name.
  character(len=48) :: units            ! A variable's units.
  character(len=240) :: longname        ! A variable's longname.
  character(len=8) :: hor_grid, z_grid  ! Variable grid info.
  type(axes_grp), pointer :: axes

  call query_vardesc(var_desc, units=units, longname=longname, hor_grid=hor_grid, &
                     z_grid=z_grid, caller="ocean_register_diag")

  ! Use the hor_grid and z_grid components of vardesc to determine the
  ! desired axes to register the diagnostic field for.
  select case (z_grid)

    case ("L")
      select case (hor_grid)
        case ("q")
          axes => diag_cs%axesBL
        case ("h")
          axes => diag_cs%axesTL
        case ("u")
          axes => diag_cs%axesCuL
        case ("v")
          axes => diag_cs%axesCvL
        case ("Bu")
          axes => diag_cs%axesBL
        case ("T")
          axes => diag_cs%axesTL
        case ("Cu")
          axes => diag_cs%axesCuL
        case ("Cv")
          axes => diag_cs%axesCvL
        case ("z")
          axes => diag_cs%axeszL
        case default
          call MOM_error(FATAL, "ocean_register_diag: " // &
              "unknown hor_grid component "//trim(hor_grid))
      end select

    case ("i")
      select case (hor_grid)
        case ("q")
          axes => diag_cs%axesBi
        case ("h")
          axes => diag_cs%axesTi
        case ("u")
          axes => diag_cs%axesCui
        case ("v")
          axes => diag_cs%axesCvi
        case ("Bu")
          axes => diag_cs%axesBi
        case ("T")
          axes => diag_cs%axesTi
        case ("Cu")
          axes => diag_cs%axesCui
        case ("Cv")
          axes => diag_cs%axesCvi
        case ("z")
          axes => diag_cs%axeszi
        case default
          call MOM_error(FATAL, "ocean_register_diag: " // &
            "unknown hor_grid component "//trim(hor_grid))
      end select

    case ("1")
      select case (hor_grid)
        case ("q")
          axes => diag_cs%axesB1
        case ("h")
          axes => diag_cs%axesT1
        case ("u")
          axes => diag_cs%axesCu1
        case ("v")
          axes => diag_cs%axesCv1
        case ("Bu")
          axes => diag_cs%axesB1
        case ("T")
          axes => diag_cs%axesT1
        case ("Cu")
          axes => diag_cs%axesCu1
        case ("Cv")
          axes => diag_cs%axesCv1
        case default
          call MOM_error(FATAL, "ocean_register_diag: " // &
            "unknown hor_grid component "//trim(hor_grid))
      end select

    case default
      call MOM_error(FATAL,&
        "ocean_register_diag: unknown z_grid component "//trim(z_grid))
  end select

  ocean_register_diag = register_diag_field("ocean_model", trim(var_name), &
          axes, day, trim(longname), trim(units),  missing_value = -1.0e+34)

end function ocean_register_diag

subroutine diag_mediator_infrastructure_init(err_msg)
  ! This subroutine initializes the FMS diag_manager.
  character(len=*), optional, intent(out)   :: err_msg

  call diag_manager_init(err_msg=err_msg)
end subroutine diag_mediator_infrastructure_init

!> diag_mediator_init initializes the MOM diag_mediator and opens the available
!! diagnostics file, if appropriate.
subroutine diag_mediator_init(G, nz, param_file, diag_cs, doc_file_dir)
  type(ocean_grid_type), target, intent(inout) :: G  !< The ocean grid type.
  integer,                    intent(in)    :: nz    !< The number of layers in the model's native grid.
  type(param_file_type),      intent(in)    :: param_file !< Parameter file structure
  type(diag_ctrl),            intent(inout) :: diag_cs !< A pointer to a type with many variables
                                                     !! used for diagnostics
  character(len=*), optional, intent(in)    :: doc_file_dir !< A directory in which to create the
                                                     !! file

  ! This subroutine initializes the diag_mediator and the diag_manager.
  ! The grid type should have its dimensions set by this point, but it
  ! is not necessary that the metrics and axis labels be set up yet.
  integer :: ios, i, new_unit
  logical :: opened, new_file
  character(len=8)   :: this_pe
  character(len=240) :: doc_file, doc_file_dflt, doc_path
  character(len=240), allocatable :: diag_coords(:)
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod  = "MOM_diag_mediator" ! This module's name.

  id_clock_diag_mediator = cpu_clock_id('(Ocean diagnostics framework)', grain=CLOCK_MODULE)
  id_clock_diag_remap = cpu_clock_id('(Ocean diagnostics remapping)', grain=CLOCK_ROUTINE)
  id_clock_diag_grid_updates = cpu_clock_id('(Ocean diagnostics grid updates)', grain=CLOCK_ROUTINE)

  ! Allocate and initialize list of all diagnostics (and variants)
  allocate(diag_cs%diags(DIAG_ALLOC_CHUNK_SIZE))
  diag_cs%next_free_diag_id = 1
  do i=1, DIAG_ALLOC_CHUNK_SIZE
    call initialize_diag_type(diag_cs%diags(i))
  enddo

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")

  call get_param(param_file, mod, 'NUM_DIAG_COORDS', diag_cs%num_diag_coords, &
                 'The number of diagnostic vertical coordinates to use.\n'//&
                 'For each coordinate, an entry in DIAG_COORDS must be provided.', &
                 default=1)
  if (diag_cs%num_diag_coords>0) then
    allocate(diag_coords(diag_cs%num_diag_coords))
    if (diag_cs%num_diag_coords==1) then ! The default is to provide just one instance of Z*
      call get_param(param_file, mod, 'DIAG_COORDS', diag_coords, &
                 'A list of string tuples associating diag_table modules to\n'//&
                 'a coordinate definition used for diagnostics. Each string\n'//&
                 'is of the form "MODULE_SUFFIX PARAMETER_SUFFIX COORDINATE_NAME".', &
                 default='z Z ZSTAR')
    else ! If using more than 1 diagnostic coordinate, all must be explicitly defined
      call get_param(param_file, mod, 'DIAG_COORDS', diag_coords, &
                 'A list of string tuples associating diag_table modules to\n'//&
                 'a coordinate definition used for diagnostics. Each string\n'//&
                 'is of the form "MODULE_SUFFIX,PARAMETER_SUFFIX,COORDINATE_NAME".', &
                 fail_if_missing=.true.)
    endif
    allocate(diag_cs%diag_remap_cs(diag_cs%num_diag_coords))
    ! Initialize each diagnostic vertical coordinate
    do i=1, diag_cs%num_diag_coords
      call diag_remap_init(diag_cs%diag_remap_cs(i), diag_coords(i))
    enddo
    deallocate(diag_coords)
  endif

  ! Keep pointers grid, h, T, S needed diagnostic remapping
  diag_cs%G => G
  diag_cs%h => null()
  diag_cs%T => null()
  diag_cs%S => null()
  diag_cs%eqn_of_state => null()

#if defined(DEBUG) || defined(__DO_SAFETY_CHECKS__)
  allocate(diag_cs%h_old(G%isd:G%ied,G%jsd:G%jed,nz))
  diag_cs%h_old(:,:,:) = 0.0
#endif

  diag_cs%is = G%isc - (G%isd-1) ; diag_cs%ie = G%iec - (G%isd-1)
  diag_cs%js = G%jsc - (G%jsd-1) ; diag_cs%je = G%jec - (G%jsd-1)
  diag_cs%isd = G%isd ; diag_cs%ied = G%ied
  diag_cs%jsd = G%jsd ; diag_cs%jed = G%jed

  if (is_root_pe() .and. (diag_CS%doc_unit < 0)) then
    write(this_pe,'(i6.6)') PE_here()
    doc_file_dflt = "available_diags."//this_pe
    call get_param(param_file, mod, "AVAILABLE_DIAGS_FILE", doc_file, &
                 "A file into which to write a list of all available \n"//&
                 "ocean diagnostics that can be included in a diag_table.", &
                 default=doc_file_dflt, do_not_log=(diag_CS%doc_unit/=-1))
    if (len_trim(doc_file) > 0) then
      new_file = .true. ; if (diag_CS%doc_unit /= -1) new_file = .false.
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

      diag_CS%doc_unit = new_unit

      if (new_file) then
        open(diag_CS%doc_unit, file=trim(doc_path), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='REPLACE', iostat=ios)
      else ! This file is being reopened, and should be appended.
        open(diag_CS%doc_unit, file=trim(doc_path), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='OLD', position='APPEND', iostat=ios)
      endif
      inquire(diag_CS%doc_unit, opened=opened)
      if ((.not.opened) .or. (ios /= 0)) then
        call MOM_error(FATAL, "Failed to open available diags file "//trim(doc_path)//".")
      endif
    endif
  endif

end subroutine diag_mediator_init

subroutine diag_set_state_ptrs(h, T, S, eqn_of_state, diag_cs)

  real, dimension(:,:,:), target, intent(in) :: h, T, S
  type(EOS_type),  pointer, intent(in)       :: eqn_of_state !< Equation of state structure
  type(diag_ctrl), intent(inout) :: diag_cs

  !  (inout) diag_cs - diag mediator control structure
  !  (in) h - a pointer to model thickness
  !  (in) T - a pointer to model temperature
  !  (in) S - a pointer to model salinity

  ! Keep pointers to h, T, S needed for the diagnostic remapping
  diag_cs%h => h
  diag_cs%T => T
  diag_cs%S => S
  diag_cs%eqn_of_state => eqn_of_state

end subroutine

!> Build/update vertical grids for diagnostic remapping.
!! \note The target grids need to be updated whenever sea surface
!! height changes.
subroutine diag_update_remap_grids(diag_cs)
  type(diag_ctrl), intent(inout) :: diag_cs !< Diagnostics control structure
  ! Local variables
  integer :: i

  if (id_clock_diag_grid_updates>0) call cpu_clock_begin(id_clock_diag_grid_updates)

  do i=1, diag_cs%num_diag_coords
    call diag_remap_update(diag_cs%diag_remap_cs(i), &
                           diag_cs%G, diag_cs%h, diag_cs%T, diag_cs%S, &
                           diag_cs%eqn_of_state)
  enddo

#if defined(DEBUG) || defined(__DO_SAFETY_CHECKS__)
  ! Keep a copy of H - used to check whether grids are up-to-date
  ! when doing remapping.
  diag_cs%h_old(:,:,:) = diag_cs%h(:,:,:)
#endif

  if (id_clock_diag_grid_updates>0) call cpu_clock_end(id_clock_diag_grid_updates)

end subroutine diag_update_remap_grids

!> diag_masks_set sets up the 2d and 3d masks for diagnostics
subroutine diag_masks_set(G, nz, missing_value, diag_cs)
  type(ocean_grid_type), target, intent(in) :: G  !< The ocean grid type.
  integer,                       intent(in) :: nz !< The number of layers in the model's native grid.
  real,                          intent(in) :: missing_value !< A value to use for masked points.
  type(diag_ctrl),               pointer    :: diag_cs !< A pointer to a type with many variables
                                                  !! used for diagnostics
  ! Local variables
  integer :: k

  diag_cs%mask2dT => G%mask2dT
  diag_cs%mask2dBu=> G%mask2dBu
  diag_cs%mask2dCu=> G%mask2dCu
  diag_cs%mask2dCv=> G%mask2dCv
  allocate(diag_cs%mask3dTL(G%isd:G%ied,G%jsd:G%jed,1:nz))
  allocate(diag_cs%mask3dBL(G%IsdB:G%IedB,G%JsdB:G%JedB,1:nz))
  allocate(diag_cs%mask3dCuL(G%IsdB:G%IedB,G%jsd:G%jed,1:nz))
  allocate(diag_cs%mask3dCvL(G%isd:G%ied,G%JsdB:G%JedB,1:nz))
  do k=1,nz
    diag_cs%mask3dTL(:,:,k)  = diag_cs%mask2dT (:,:)
    diag_cs%mask3dBL(:,:,k) = diag_cs%mask2dBu(:,:)
    diag_cs%mask3dCuL(:,:,k) = diag_cs%mask2dCu(:,:)
    diag_cs%mask3dCvL(:,:,k) = diag_cs%mask2dCv(:,:)
  enddo
  allocate(diag_cs%mask3dTi(G%isd:G%ied,G%jsd:G%jed,1:nz+1))
  allocate(diag_cs%mask3dBi(G%IsdB:G%IedB,G%JsdB:G%JedB,1:nz+1))
  allocate(diag_cs%mask3dCui(G%IsdB:G%IedB,G%jsd:G%jed,1:nz+1))
  allocate(diag_cs%mask3dCvi(G%isd:G%ied,G%JsdB:G%JedB,1:nz+1))
  do k=1,nz+1
    diag_cs%mask3dTi(:,:,k)  = diag_cs%mask2dT (:,:)
    diag_cs%mask3dBi(:,:,k) = diag_cs%mask2dBu(:,:)
    diag_cs%mask3dCui(:,:,k) = diag_cs%mask2dCu(:,:)
    diag_cs%mask3dCvi(:,:,k) = diag_cs%mask2dCv(:,:)
  enddo

  diag_cs%missing_value = missing_value

end subroutine diag_masks_set

subroutine diag_mediator_close_registration(diag_CS)
  type(diag_ctrl), intent(inout) :: diag_CS

  integer :: i

  if (diag_CS%doc_unit > -1) then
    close(diag_CS%doc_unit) ; diag_CS%doc_unit = -2
  endif

  do i=1, diag_cs%num_diag_coords
    call diag_remap_diag_registration_closed(diag_cs%diag_remap_cs(i))
  enddo

end subroutine diag_mediator_close_registration

subroutine diag_mediator_end(time, diag_CS, end_diag_manager)
  type(time_type),   intent(in)  :: time
  type(diag_ctrl), intent(inout) :: diag_cs
  logical, optional, intent(in)  :: end_diag_manager !< If true, call diag_manager_end()

  ! Local variables
  integer :: i

  if (diag_CS%doc_unit > -1) then
    close(diag_CS%doc_unit) ; diag_CS%doc_unit = -3
  endif

  deallocate(diag_cs%diags)

  do i=1, diag_cs%num_diag_coords
    call diag_remap_end(diag_cs%diag_remap_cs(i))
  enddo

  deallocate(diag_cs%mask3dTL)
  deallocate(diag_cs%mask3dBL)
  deallocate(diag_cs%mask3dCuL)
  deallocate(diag_cs%mask3dCvL)
  deallocate(diag_cs%mask3dTi)
  deallocate(diag_cs%mask3dBi)
  deallocate(diag_cs%mask3dCui)
  deallocate(diag_cs%mask3dCvi)

#if defined(DEBUG) || defined(__DO_SAFETY_CHECKS__)
  deallocate(diag_cs%h_old)
#endif

  if (present(end_diag_manager)) then
    if (end_diag_manager) call diag_manager_end(time)
  endif

end subroutine diag_mediator_end

function i2s(a,n_in)
!   "Convert the first n elements of an integer array to a string."
    integer, dimension(:), intent(in) :: a
    integer, optional    , intent(in) :: n_in
    character(len=15) :: i2s

    character(len=15) :: i2s_temp
    integer :: i,n

    n=size(a)
    if(present(n_in)) n = n_in

    i2s = ''
    do i=1,n
       write (i2s_temp, '(I4.4)') a(i)
       i2s = trim(i2s) //'_'// trim(i2s_temp)
    enddo
    i2s = adjustl(i2s)
end function i2s

!> Associates the mask pointers within diag with the appropriate mask based on the axes group.
subroutine set_diag_mask(diag, diag_cs, axes)
  type(diag_ctrl), target, intent(in) :: diag_cs !< Diag_mediator control structure
  type(diag_type), pointer, intent(inout) :: diag !< This diag type
  type(axes_grp),  intent(in) :: axes !< Axes group

  diag%mask2d => null()
  diag%mask3d => null()

  if (axes%rank .eq. 3) then
    if (axes%is_layer) then
      if (axes%is_h_point) then
        diag%mask3d => diag_cs%mask3dTL
      elseif (axes%is_q_point) then
        diag%mask3d => diag_cs%mask3dBL
      elseif (axes%is_u_point) then
        diag%mask3d => diag_cs%mask3dCuL
      elseif (axes%is_v_point) then
        diag%mask3d => diag_cs%mask3dCvL
      endif
    elseif (axes%is_interface) then
      if (axes%is_h_point) then
        diag%mask3d => diag_cs%mask3dTi
      elseif (axes%is_q_point) then
        diag%mask3d => diag_cs%mask3dBi
      elseif (axes%is_u_point) then
        diag%mask3d => diag_cs%mask3dCui
      elseif (axes%is_v_point) then
        diag%mask3d => diag_cs%mask3dCvi
      endif
    endif

    !call assert(associated(diag%mask3d), "set_diag_mask: Invalid 3d axes id."// &
    !                                     " diag:"//diag%debug_str)
  elseif(axes%rank .eq. 2) then

    if (axes%is_h_point) then
      diag%mask2d =>  diag_cs%mask2dT
    elseif (axes%is_q_point) then
      diag%mask2d =>  diag_cs%mask2dBu
    elseif (axes%is_u_point) then
        diag%mask2d =>  diag_cs%mask2dCu
    elseif (axes%is_v_point) then
        diag%mask2d =>  diag_cs%mask2dCv
    endif

    !call assert(associated(diag%mask2d), "set_diag_mask.F90: Invalid 2d axes id."// &
    !                                     " diag:"//diag%debug_str)
  endif

end subroutine set_diag_mask

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
  diag%mask2d => null()
  diag%mask3d => null()
  diag%next => null()
  diag%conversion_factor = 0.

end subroutine initialize_diag_type

! Make a new diagnostic. Either use memory which is in the array of 'primary'
! diagnostics, or if that is in use, insert it to the list of secondary diags.
subroutine alloc_diag_with_id(diag_id, diag_cs, diag)
  integer, intent(in) :: diag_id
  type(diag_ctrl), target, intent(inout) :: diag_cs
  type(diag_type), pointer, intent(out) :: diag

  ! Arguments:
  !  (in)      diag_id  - new id for the diag.
  !  (inout)   diag_cs  - structure used to regulate diagnostic output
  !  (inout)   diag     - structure representing a diagnostic

  type(diag_type), pointer :: tmp

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
                              diag_CS, long_name, units, standard_name)
  logical,          intent(in) :: used !< Whether this diagnostic was in the diag_table or not
  character(len=*), intent(in) :: module_name !< Name of the diagnostic module
  character(len=*), intent(in) :: field_name !< Name of this diagnostic field
  character(len=*), intent(in) :: cell_methods_string !< The spatial component of the CF cell_methods attribute
  character(len=*), intent(in) :: comment !< A comment to append after [Used|Unused]
  type(diag_ctrl),  intent(in) :: diag_CS  !< The diagnotics control structure
  character(len=*), optional, intent(in) :: long_name !< CF long name of diagnostic
  character(len=*), optional, intent(in) :: units !< Units for diagnostic
  character(len=*), optional, intent(in) :: standard_name !< CF standardized name of diagnostic
  ! Local variables
  character(len=240) :: mesg

  if (used) then
    mesg = '"'//trim(module_name)//'", "'//trim(field_name)//'"  [Used]'
  else
    mesg = '"'//trim(module_name)//'", "'//trim(field_name)//'"  [Unused]'
  endif
  if (len(trim((comment)))>0) then
    write(diag_CS%doc_unit, '(a,x,"(",a,")")') trim(mesg),trim(comment)
  else
    write(diag_CS%doc_unit, '(a)') trim(mesg)
  endif
  if (present(long_name)) call describe_option("long_name", long_name, diag_CS)
  if (present(units)) call describe_option("units", units, diag_CS)
  if (present(standard_name)) &
    call describe_option("standard_name", standard_name, diag_CS)
  if (len(trim((cell_methods_string)))>0) &
    call describe_option("cell_methods", trim(cell_methods_string), diag_CS)

end subroutine log_available_diag

end module MOM_diag_mediator
