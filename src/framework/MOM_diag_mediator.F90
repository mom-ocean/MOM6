!> The subroutines here provide convenient wrappers to the fms diag_manager
!! interfaces with additional diagnostic capabilies.
module MOM_diag_mediator

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_checksums,        only : chksum_general
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
use MOM_diag_remap,       only : diag_remap_calc_hmask
use MOM_diag_remap,       only : diag_remap_init, diag_remap_end, diag_remap_do_remap
use MOM_diag_remap,       only : vertically_reintegrate_diag_field, vertically_interpolate_diag_field
use MOM_diag_remap,       only : diag_remap_configure_axes, diag_remap_axes_configured
use MOM_diag_remap,       only : diag_remap_get_axes_info, diag_remap_set_active
use MOM_diag_remap,       only : diag_remap_diag_registration_closed
use MOM_diag_remap,       only : horizontally_average_diag_field

use diag_axis_mod, only : get_diag_axis_name
use diag_data_mod, only : null_axis_id
use diag_manager_mod, only : diag_manager_init, diag_manager_end
use diag_manager_mod, only : send_data, diag_axis_init, diag_field_add_attribute
! The following module is needed for PGI since the following line does not compile with PGI 6.5.0
! was: use diag_manager_mod, only : register_diag_field_fms=>register_diag_field
use MOM_diag_manager_wrapper, only : register_diag_field_fms
use diag_manager_mod, only : register_static_field_fms=>register_static_field
use diag_manager_mod, only : get_diag_field_id_fms=>get_diag_field_id
use diag_manager_mod, only : DIAG_FIELD_NOT_FOUND

implicit none ; private

#undef __DO_SAFETY_CHECKS__
#define IMPLIES(A, B) ((.not. (A)) .or. (B))

public set_axes_info, post_data, register_diag_field, time_type
public set_masks_for_axes
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
public register_cell_measure, diag_associate_volume_cell_measure
public diag_get_volume_cell_measure_dm_id
public diag_set_state_ptrs, diag_update_remap_grids
public diag_grid_storage_init, diag_grid_storage_end
public diag_copy_diag_to_storage, diag_copy_storage_to_diag
public diag_save_grids, diag_restore_grids

!> Make a diagnostic available for averaging or output.
interface post_data
  module procedure post_data_3d, post_data_2d, post_data_1d_k, post_data_0d
end interface post_data

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
  character(len=9) :: v_cell_method = '' !< Default nature of data representation, if axes group
                                         !! includes vertical direction.
  ! For remapping
  integer :: nz = 0 !< Vertical dimension of diagnostic
  integer :: vertical_coordinate_number = 0 !< Index of the corresponding diag_remap_ctrl for this axis group
  ! For detecting position on the grid
  logical :: is_h_point = .false. !< If true, indicates that this axes group is for an h-point located field.
  logical :: is_q_point = .false. !< If true, indicates that this axes group is for a q-point located field.
  logical :: is_u_point = .false. !< If true, indicates that this axes group is for a u-point located field.
  logical :: is_v_point = .false. !< If true, indicates that this axes group is for a v-point located field.
  logical :: is_layer = .false. !< If true, indicates that this axes group is for a layer vertically-located field.
  logical :: is_interface = .false. !< If true, indicates that this axes group is for an interface
                                    !! vertically-located field.
  logical :: is_native = .true. !< If true, indicates that this axes group is for a native model grid.
                                !! False for any other grid. Used for rank>2.
  logical :: needs_remapping = .false. !< If true, indicates that this axes group is for a intensive layer-located
                                       !! field that must be remapped to these axes. Used for rank>2.
  logical :: needs_interpolating = .false. !< If true, indicates that this axes group is for a sampled
                                         !! interface-located field that must be interpolated to
                                         !! these axes. Used for rank>2.
  ! For horizontally averaged diagnositcs (applies to 2d and 3d fields only)
  type(axes_grp), pointer :: xyave_axes => null() !< The associated 1d axes for horizontall area-averaged diagnostics
  ! ID's for cell_measures
  integer :: id_area = -1 !< The diag_manager id for area to be used for cell_measure of variables with this axes_grp.
  integer :: id_volume = -1 !< The diag_manager id for volume to be used for cell_measure of variables
                            !! with this axes_grp.
  ! For masking
  real, pointer, dimension(:,:)   :: mask2d => null() !< Mask for 2d (x-y) axes
  real, pointer, dimension(:,:,:) :: mask3d => null() !< Mask for 3d axes
end type axes_grp

!> Contains an array to store a diagnostic target grid
type, private :: diag_grids_type
  real, dimension(:,:,:), allocatable :: h  !< Target grid for remapped coordinate
end type diag_grids_type

!> Stores all the remapping grids and the model's native space thicknesses
type, public :: diag_grid_storage
  integer                                          :: num_diag_coords !< Number of target coordinates
  real, dimension(:,:,:), allocatable              :: h_state         !< Layer thicknesses in native space
  type(diag_grids_type), dimension(:), allocatable :: diag_grids      !< Primarily empty, except h field
end type diag_grid_storage

!> This type is used to represent a diagnostic at the diag_mediator level.
!!
!! There can be both 'primary' and 'seconday' diagnostics. The primaries
!! reside in the diag_cs%diags array. They have an id which is an index
!! into this array. The secondaries are 'variations' on the primary diagnostic.
!! For example the CMOR diagnostics are secondary. The secondary diagnostics
!! are kept in a list with the primary diagnostic as the head.
type, private :: diag_type
  logical :: in_use !< True if this entry is being used.
  integer :: fms_diag_id !< Underlying FMS diag_manager id.
  integer :: fms_xyave_diag_id = -1 !< For a horizontally area-averaged diagnostic.
  character(64) :: debug_str = '' !< For FATAL errors and debugging.
  type(axes_grp), pointer :: axes => null() !< The axis group for this diagnostic
  type(diag_type), pointer :: next => null() !< Pointer to the next diagnostic
  real :: conversion_factor = 0. !< A factor to multiply data by before posting to FMS, if non-zero.
  logical :: v_extensive = .false. !< True for vertically extensive fields (vertically integrated).
                                   !! False for intensive (concentrations).
end type diag_type

!> The following data type a list of diagnostic fields an their variants,
!! as well as variables that control the handling of model output.
type, public :: diag_ctrl
  integer :: available_diag_doc_unit = -1 !< The unit number of a diagnostic documentation file.
                                          !! This file is open if available_diag_doc_unit is > 0.
  integer :: chksum_diag_doc_unit = -1 !< The unit number of a diagnostic documentation file.
                                          !! This file is open if available_diag_doc_unit is > 0.
  logical :: diag_as_chksum !< If true, log chksums in a text file instead of posting diagnostics

! The following fields are used for the output of the data.
  integer :: is  !< The start i-index of cell centers within the computational domain
  integer :: ie  !< The end i-index of cell centers within the computational domain
  integer :: js  !< The start j-index of cell centers within the computational domain
  integer :: je  !< The end j-index of cell centers within the computational domain

  integer :: isd !< The start i-index of cell centers within the data domain
  integer :: ied !< The end i-index of cell centers within the data domain
  integer :: jsd !< The start j-index of cell centers within the data domain
  integer :: jed !< The end j-index of cell centers within the data domain
  real :: time_int              !< The time interval in s for any fields
                                !! that are offered for averaging.
  type(time_type) :: time_end   !< The end time of the valid
                                !! interval for any offered field.
  logical :: ave_enabled = .false. !< True if averaging is enabled.

  !>@{ The following are 3D and 2D axis groups defined for output.  The names
  !! indicate the horizontal (B, T, Cu, or Cv) and vertical (L, i, or 1) locations.
  type(axes_grp) :: axesBL, axesTL, axesCuL, axesCvL
  type(axes_grp) :: axesBi, axesTi, axesCui, axesCvi
  type(axes_grp) :: axesB1, axesT1, axesCu1, axesCv1
  !!@}
  type(axes_grp) :: axesZi !< A 1-D z-space axis at interfaces
  type(axes_grp) :: axesZL !< A 1-D z-space axis at layer centers
  type(axes_grp) :: axesNull !< An axis group for scalars

  real, dimension(:,:),   pointer :: mask2dT   => null() !< 2D mask array for cell-center points
  real, dimension(:,:),   pointer :: mask2dBu  => null() !< 2D mask array for cell-corner points
  real, dimension(:,:),   pointer :: mask2dCu  => null() !< 2D mask array for east-face points
  real, dimension(:,:),   pointer :: mask2dCv  => null() !< 2D mask array for north-face points
  !>@{ 3D mask arrays for diagnostics at layers (mask...L) and interfaces (mask...i)
  real, dimension(:,:,:), pointer :: mask3dTL  => null()
  real, dimension(:,:,:), pointer :: mask3dBL  => null()
  real, dimension(:,:,:), pointer :: mask3dCuL => null()
  real, dimension(:,:,:), pointer :: mask3dCvL => null()
  real, dimension(:,:,:), pointer :: mask3dTi  => null()
  real, dimension(:,:,:), pointer :: mask3dBi  => null()
  real, dimension(:,:,:), pointer :: mask3dCui => null()
  real, dimension(:,:,:), pointer :: mask3dCvi => null()
  !!@}

! Space for diagnostics is dynamically allocated as it is needed.
! The chunk size is how much the array should grow on each new allocation.
#define DIAG_ALLOC_CHUNK_SIZE 100
  type(diag_type), dimension(:), allocatable :: diags !< The list of diagnostics
  integer :: next_free_diag_id !< The next unused diagnostic ID

  !> default missing value to be sent to ALL diagnostics registrations
  real :: missing_value = -1.0e+34

  !> Number of diagnostic vertical coordinates (remapped)
  integer :: num_diag_coords
  !> Control structure for each possible coordinate
  type(diag_remap_ctrl), dimension(:), allocatable :: diag_remap_cs
  type(diag_grid_storage) :: diag_grid_temp !< Stores the remapped diagnostic grid
  logical :: diag_grid_overridden = .false. !< True if the diagnostic grids have been overriden

  type(axes_grp), dimension(:), allocatable :: &
    remap_axesZL, &  !< The 1-D z-space cell-centered axis for remapping
    remap_axesZi     !< The 1-D z-space interface axis for remapping
  !!@{
  type(axes_grp), dimension(:), allocatable :: remap_axesTL, remap_axesBL, remap_axesCuL, remap_axesCvL
  type(axes_grp), dimension(:), allocatable :: remap_axesTi, remap_axesBi, remap_axesCui, remap_axesCvi
  !!@}

  ! Pointer to H, G and T&S needed for remapping
  real, dimension(:,:,:), pointer :: h => null() !< The thicknesses needed for remapping
  real, dimension(:,:,:), pointer :: T => null() !< The temperatures needed for remapping
  real, dimension(:,:,:), pointer :: S => null() !< The salinities needed for remapping
  type(EOS_type),  pointer :: eqn_of_state => null() !< The equation of state type
  type(ocean_grid_type), pointer :: G => null()  !< The ocean grid type
  type(verticalGrid_type), pointer :: GV => null()  !< The model's vertical ocean grid

  !> The volume cell measure (special diagnostic) manager id
  integer :: volume_cell_measure_dm_id = -1

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

  if (G%symmetric) then
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

  ! Axis group for special null axis from diag manager
  call define_axes_group(diag_cs, (/ null_axis_id /), diag_cs%axesNull)

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

       !! \note Remapping for B points is not yet implemented so needs_remapping is not
       !! provided for remap_axesBL
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

  call diag_grid_storage_init(diag_CS%diag_grid_temp, G, diag_CS)

end subroutine set_axes_info

!> set_masks_for_axes sets up the 2d and 3d masks for diagnostics using the current grid
!! recorded after calling diag_update_remap_grids()
subroutine set_masks_for_axes(G, diag_cs)
  type(ocean_grid_type), target, intent(in) :: G !< The ocean grid type.
  type(diag_ctrl),               pointer    :: diag_cs !< A pointer to a type with many variables
                                                       !! used for diagnostics
  ! Local variables
  integer :: c, nk, i, j, k
  type(axes_grp), pointer :: axes => NULL(), h_axes => NULL() ! Current axes, for convenience

  do c=1, diag_cs%num_diag_coords
    ! This vertical coordinate has been configured so can be used.
    if (diag_remap_axes_configured(diag_cs%diag_remap_cs(c))) then

      ! Level/layer h-points in diagnostic coordinate
      axes => diag_cs%remap_axesTL(c)
      nk = axes%nz
      allocate( axes%mask3d(G%isd:G%ied,G%jsd:G%jed,nk) ) ; axes%mask3d(:,:,:) = 0.
      call diag_remap_calc_hmask(diag_cs%diag_remap_cs(c), G, axes%mask3d)

      h_axes => diag_cs%remap_axesTL(c) ! Use the h-point masks to generate the u-, v- and q- masks

      ! Level/layer u-points in diagnostic coordinate
      axes => diag_cs%remap_axesCuL(c)
      call assert(axes%nz == nk, 'set_masks_for_axes: vertical size mismatch at u-layers')
      call assert(.not. associated(axes%mask3d), 'set_masks_for_axes: already associated')
      allocate( axes%mask3d(G%IsdB:G%IedB,G%jsd:G%jed,nk) ) ; axes%mask3d(:,:,:) = 0.
      do k = 1, nk ; do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
        if (h_axes%mask3d(i,j,k) + h_axes%mask3d(i+1,j,k) > 0.) axes%mask3d(I,j,k) = 1.
      enddo ; enddo ; enddo

      ! Level/layer v-points in diagnostic coordinate
      axes => diag_cs%remap_axesCvL(c)
      call assert(axes%nz == nk, 'set_masks_for_axes: vertical size mismatch at v-layers')
      call assert(.not. associated(axes%mask3d), 'set_masks_for_axes: already associated')
      allocate( axes%mask3d(G%isd:G%ied,G%JsdB:G%JedB,nk) ) ; axes%mask3d(:,:,:) = 0.
      do k = 1, nk ; do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
        if (h_axes%mask3d(i,j,k) + h_axes%mask3d(i,j+1,k) > 0.) axes%mask3d(i,J,k) = 1.
      enddo ; enddo ; enddo

      ! Level/layer q-points in diagnostic coordinate
      axes => diag_cs%remap_axesBL(c)
      call assert(axes%nz == nk, 'set_masks_for_axes: vertical size mismatch at q-layers')
      call assert(.not. associated(axes%mask3d), 'set_masks_for_axes: already associated')
      allocate( axes%mask3d(G%IsdB:G%IedB,G%JsdB:G%JedB,nk) ) ; axes%mask3d(:,:,:) = 0.
      do k = 1, nk ; do J=G%jsc-1,G%jec ; do I=G%isc-1,G%iec
        if (h_axes%mask3d(i,j,k) + h_axes%mask3d(i+1,j+1,k) + &
            h_axes%mask3d(i+1,j,k) + h_axes%mask3d(i,j+1,k) > 0.) axes%mask3d(I,J,k) = 1.
      enddo ; enddo ; enddo

      ! Interface h-points in diagnostic coordinate (w-point)
      axes => diag_cs%remap_axesTi(c)
      call assert(axes%nz == nk, 'set_masks_for_axes: vertical size mismatch at h-interfaces')
      call assert(.not. associated(axes%mask3d), 'set_masks_for_axes: already associated')
      allocate( axes%mask3d(G%isd:G%ied,G%jsd:G%jed,nk+1) ) ; axes%mask3d(:,:,:) = 0.
      do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
        if (h_axes%mask3d(i,j,1) > 0.) axes%mask3d(i,J,1) = 1.
        do K = 2, nk
          if (h_axes%mask3d(i,j,k-1) + h_axes%mask3d(i,j,k) > 0.) axes%mask3d(i,J,k) = 1.
        enddo
        if (h_axes%mask3d(i,j,nk) > 0.) axes%mask3d(i,J,nk+1) = 1.
      enddo ; enddo

      h_axes => diag_cs%remap_axesTi(c) ! Use the w-point masks to generate the u-, v- and q- masks

      ! Interface u-points in diagnostic coordinate
      axes => diag_cs%remap_axesCui(c)
      call assert(axes%nz == nk, 'set_masks_for_axes: vertical size mismatch at u-interfaces')
      call assert(.not. associated(axes%mask3d), 'set_masks_for_axes: already associated')
      allocate( axes%mask3d(G%IsdB:G%IedB,G%jsd:G%jed,nk+1) ) ; axes%mask3d(:,:,:) = 0.
      do k = 1, nk+1 ; do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
        if (h_axes%mask3d(i,j,k) + h_axes%mask3d(i+1,j,k) > 0.) axes%mask3d(I,j,k) = 1.
      enddo ; enddo ; enddo

      ! Interface v-points in diagnostic coordinate
      axes => diag_cs%remap_axesCvi(c)
      call assert(axes%nz == nk, 'set_masks_for_axes: vertical size mismatch at v-interfaces')
      call assert(.not. associated(axes%mask3d), 'set_masks_for_axes: already associated')
      allocate( axes%mask3d(G%isd:G%ied,G%JsdB:G%JedB,nk+1) ) ; axes%mask3d(:,:,:) = 0.
      do k = 1, nk+1 ; do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
        if (h_axes%mask3d(i,j,k) + h_axes%mask3d(i,j+1,k) > 0.) axes%mask3d(i,J,k) = 1.
      enddo ; enddo ; enddo

      ! Interface q-points in diagnostic coordinate
      axes => diag_cs%remap_axesBi(c)
      call assert(axes%nz == nk, 'set_masks_for_axes: vertical size mismatch at q-interfaces')
      call assert(.not. associated(axes%mask3d), 'set_masks_for_axes: already associated')
      allocate( axes%mask3d(G%IsdB:G%IedB,G%JsdB:G%JedB,nk+1) ) ; axes%mask3d(:,:,:) = 0.
      do k = 1, nk ; do J=G%jsc-1,G%jec ; do I=G%isc-1,G%iec
        if (h_axes%mask3d(i,j,k) + h_axes%mask3d(i+1,j+1,k) + &
            h_axes%mask3d(i+1,j,k) + h_axes%mask3d(i,j+1,k) > 0.) axes%mask3d(I,J,k) = 1.
      enddo ; enddo ; enddo

    endif
  enddo

end subroutine set_masks_for_axes

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
      diag_cs%remap_axesTi(i)%id_area = fms_id
    enddo
  endif
  if (present(id_area_q)) then
    fms_id = diag_cs%diags(id_area_q)%fms_diag_id
    diag_cs%axesB1%id_area = fms_id
    diag_cs%axesBi%id_area = fms_id
    diag_cs%axesBL%id_area = fms_id
    do i=1, diag_cs%num_diag_coords
      diag_cs%remap_axesBL(i)%id_area = fms_id
      diag_cs%remap_axesBi(i)%id_area = fms_id
    enddo
  endif
end subroutine diag_register_area_ids

!> Sets a handle inside diagnostics mediator to associate 3d cell measures
subroutine register_cell_measure(G, diag, Time)
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean grid structure
  type(diag_ctrl), target, intent(inout) :: diag !< Regulates diagnostic output
  type(time_type),         intent(in)    :: Time !< Model time
  ! Local variables
  integer :: id
  id = register_diag_field('ocean_model', 'volcello', diag%axesTL, &
                           Time, 'Ocean grid-cell volume', 'm3', &
                           standard_name='ocean_volume', v_extensive=.true., &
                           x_cell_method='sum', y_cell_method='sum')
  call diag_associate_volume_cell_measure(diag, id)

end subroutine register_cell_measure

!> Attaches the id of cell volumes to axes groups for use with cell_measures
subroutine diag_associate_volume_cell_measure(diag_cs, id_h_volume)
  type(diag_ctrl),   intent(inout) :: diag_cs     !< Diagnostics control structure
  integer,           intent(in)    :: id_h_volume !< Diag_manager id for volume of h-cells
  ! Local variables
  type(diag_type), pointer :: tmp => NULL()

  if (id_h_volume<=0) return ! Do nothing
  diag_cs%volume_cell_measure_dm_id = id_h_volume ! Record for diag_get_volume_cell_measure_dm_id()

  ! Set the cell measure for this axes group to the FMS id in this coordinate system
  diag_cs%diags(id_h_volume)%axes%id_volume = diag_cs%diags(id_h_volume)%fms_diag_id

  tmp => diag_cs%diags(id_h_volume)%next ! First item in the list, if any
  do while (associated(tmp))
    ! Set the cell measure for this axes group to the FMS id in this coordinate system
    tmp%axes%id_volume = tmp%fms_diag_id
    tmp => tmp%next ! Move to next axes group for this field
  enddo

end subroutine diag_associate_volume_cell_measure

!> Returns diag_manager id for cell measure of h-cells
integer function diag_get_volume_cell_measure_dm_id(diag_cs)
  type(diag_ctrl),   intent(in) :: diag_cs   !< Diagnostics control structure

  diag_get_volume_cell_measure_dm_id = diag_cs%volume_cell_measure_dm_id

end function diag_get_volume_cell_measure_dm_id

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
  character(len=*), optional, intent(in)  :: x_cell_method !< A x-direction cell method used to construct the
                                                           !! "cell_methods" attribute in CF convention
  character(len=*), optional, intent(in)  :: y_cell_method !< A y-direction cell method used to construct the
                                                           !! "cell_methods" attribute in CF convention
  character(len=*), optional, intent(in)  :: v_cell_method !< A vertical direction cell method used to construct
                                                        !! the "cell_methods" attribute in CF convention
  logical,          optional, intent(in)  :: is_h_point !< If true, indicates this axes group for h-point
                                                        !! located fields
  logical,          optional, intent(in)  :: is_q_point !< If true, indicates this axes group for q-point
                                                        !! located fields
  logical,          optional, intent(in)  :: is_u_point !< If true, indicates this axes group for
                                                        !! u-point located fields
  logical,          optional, intent(in)  :: is_v_point !< If true, indicates this axes group for
                                                        !! v-point located fields
  logical,          optional, intent(in)  :: is_layer   !< If true, indicates that this axes group is
                                                        !! for a layer vertically-located field.
  logical,          optional, intent(in)  :: is_interface !< If true, indicates that this axes group
                                                        !! is for an interface vertically-located field.
  logical,          optional, intent(in)  :: is_native  !< If true, indicates that this axes group is
                                                        !! for a native model grid. False for any other grid.
  logical,          optional, intent(in)  :: needs_remapping !< If true, indicates that this axes group is
                                                        !! for a intensive layer-located field that must
                                                        !! be remapped to these axes. Used for rank>2.
  logical,          optional, intent(in)  :: needs_interpolating !< If true, indicates that this axes group
                                                        !! is for a sampled interface-located field that must
                                                        !! be interpolated to these axes. Used for rank>2.
  type(axes_grp),   optional, target      :: xyave_axes !< The corresponding axes group for horizontally
                                                        !! area-average diagnostics
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

  ! Setup masks for this axes group
  axes%mask2d => null()
  if (axes%rank==2) then
    if (axes%is_h_point) axes%mask2d => diag_cs%mask2dT
    if (axes%is_u_point) axes%mask2d => diag_cs%mask2dCu
    if (axes%is_v_point) axes%mask2d => diag_cs%mask2dCv
    if (axes%is_q_point) axes%mask2d => diag_cs%mask2dBu
  endif
  ! A static 3d mask for non-native coordinates can only be setup when a grid is available
  axes%mask3d => null()
  if (axes%rank==3 .and. axes%is_native) then
    ! Native variables can/should use the native masks copied into diag_cs
    if (axes%is_layer) then
      if (axes%is_h_point) axes%mask3d => diag_cs%mask3dTL
      if (axes%is_u_point) axes%mask3d => diag_cs%mask3dCuL
      if (axes%is_v_point) axes%mask3d => diag_cs%mask3dCvL
      if (axes%is_q_point) axes%mask3d => diag_cs%mask3dBL
    elseif (axes%is_interface) then
      if (axes%is_h_point) axes%mask3d => diag_cs%mask3dTi
      if (axes%is_u_point) axes%mask3d => diag_cs%mask3dCui
      if (axes%is_v_point) axes%mask3d => diag_cs%mask3dCvi
      if (axes%is_q_point) axes%mask3d => diag_cs%mask3dBi
    endif
  endif

end subroutine define_axes_group

!> Set up the array extents for doing diagnostics
subroutine set_diag_mediator_grid(G, diag_cs)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
  type(diag_ctrl),  intent(inout) :: diag_CS !< Structure used to regulate diagnostic output

  diag_cs%is = G%isc - (G%isd-1) ; diag_cs%ie = G%iec - (G%isd-1)
  diag_cs%js = G%jsc - (G%jsd-1) ; diag_cs%je = G%jec - (G%jsd-1)
  diag_cs%isd = G%isd ; diag_cs%ied = G%ied
  diag_cs%jsd = G%jsd ; diag_cs%jed = G%jed

end subroutine set_diag_mediator_grid

!> Make a real scalar diagnostic available for averaging or output
subroutine post_data_0d(diag_field_id, field, diag_cs, is_static)
  integer,           intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                 !! previous call to register_diag_field.
  real,              intent(in) :: field         !< real value being offered for output or averaging
  type(diag_ctrl), target, intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.

  ! Local variables
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

!> Make a real 1-d array diagnostic available for averaging or output
subroutine post_data_1d_k(diag_field_id, field, diag_cs, is_static)
  integer,           intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                 !! previous call to register_diag_field.
  real, target,      intent(in) :: field(:)      !< 1-d array being offered for output or averaging
  type(diag_ctrl), target, intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.

  ! Local variables
  logical :: used  ! The return value of send_data is not used for anything.
  real, dimension(:), pointer :: locfield => NULL()
  logical :: is_stat
  integer :: k, ks, ke
  type(diag_type), pointer :: diag => null()

  if (id_clock_diag_mediator>0) call cpu_clock_begin(id_clock_diag_mediator)
  is_stat = .false. ; if (present(is_static)) is_stat = is_static

  ! Iterate over list of diag 'variants', e.g. CMOR aliases.
  call assert(diag_field_id < diag_cs%next_free_diag_id, &
              'post_data_1d_k: Unregistered diagnostic id')
  diag => diag_cs%diags(diag_field_id)
  do while (associated(diag))

    if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) then
      ks = lbound(field,1) ; ke = ubound(field,1)
      allocate( locfield( ks:ke ) )

      do k=ks,ke
        if (field(k) == diag_cs%missing_value) then
          locfield(k) = diag_cs%missing_value
        else
          locfield(k) = field(k) * diag%conversion_factor
        endif
      enddo
    else
      locfield => field
    endif

    if (is_stat) then
      used = send_data(diag%fms_diag_id, locfield)
    elseif (diag_cs%ave_enabled) then
      used = send_data(diag%fms_diag_id, locfield, diag_cs%time_end, weight=diag_cs%time_int)
    endif
    if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) deallocate( locfield )

    diag => diag%next
  enddo

  if (id_clock_diag_mediator>0) call cpu_clock_end(id_clock_diag_mediator)
end subroutine post_data_1d_k

!> Make a real 2-d array diagnostic available for averaging or output
subroutine post_data_2d(diag_field_id, field, diag_cs, is_static, mask)
  integer,           intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                 !! previous call to register_diag_field.
  real,              intent(in) :: field(:,:)    !< 2-d array being offered for output or averaging
  type(diag_ctrl), target, intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.
  real,    optional, intent(in) :: mask(:,:) !< If present, use this real array as the data mask.

  ! Local variables
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

!> Make a real 2-d array diagnostic available for averaging or output
!! using a diag_type instead of an integer id.
subroutine post_data_2d_low(diag, field, diag_cs, is_static, mask)
  type(diag_type),   intent(in) :: diag       !< A structure describing the diagnostic to post
  real,    target,   intent(in) :: field(:,:) !< 2-d array being offered for output or averaging
  type(diag_ctrl),   intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.
  real,    optional, intent(in) :: mask(:,:) !< If present, use this real array as the data mask.

  ! Local variables
  real, dimension(:,:), pointer :: locfield => NULL()
  character(len=300) :: mesg
  logical :: used, is_stat
  integer :: cszi, cszj, dszi, dszj
  integer :: isv, iev, jsv, jev, i, j, chksum

  is_stat = .false. ; if (present(is_static)) is_stat = is_static

  ! Determine the propery array indices, noting that because of the (:,:)
  ! declaration of field, symmetric arrays are using a SW-grid indexing,
  ! but non-symmetric arrays are using a NE-grid indexing.  Send_data
  ! actually only uses the difference between ie and is to determine
  ! the output data size and assumes that halos are symmetric.
  isv = diag_cs%is ; iev = diag_cs%ie ; jsv = diag_cs%js ; jev = diag_cs%je

  cszi = diag_cs%ie-diag_cs%is +1 ; dszi = diag_cs%ied-diag_cs%isd +1
  cszj = diag_cs%je-diag_cs%js +1 ; dszj = diag_cs%jed-diag_cs%jsd +1
  if ( size(field,1) == dszi ) then
    isv = diag_cs%is ; iev = diag_cs%ie     ! Data domain
  elseif ( size(field,1) == dszi + 1 ) then
    isv = diag_cs%is ; iev = diag_cs%ie+1   ! Symmetric data domain
  elseif ( size(field,1) == cszi) then
    isv = 1 ; iev = cszi                    ! Computational domain
  elseif ( size(field,1) == cszi + 1 ) then
    isv = 1 ; iev = cszi+1                  ! Symmetric computational domain
  else
    write (mesg,*) " peculiar size ",size(field,1)," in i-direction\n"//&
       "does not match one of ", cszi, cszi+1, dszi, dszi+1
    call MOM_error(FATAL,"post_data_2d_low: "//trim(diag%debug_str)//trim(mesg))
  endif

  if ( size(field,2) == dszj ) then
    jsv = diag_cs%js ; jev = diag_cs%je     ! Data domain
  elseif ( size(field,2) == dszj + 1 ) then
    jsv = diag_cs%js ; jev = diag_cs%je+1   ! Symmetric data domain
  elseif ( size(field,2) == cszj ) then
    jsv = 1 ; jev = cszj                    ! Computational domain
  elseif ( size(field,2) == cszj+1 ) then
    jsv = 1 ; jev = cszj+1                  ! Symmetric computational domain
  else
    write (mesg,*) " peculiar size ",size(field,2)," in j-direction\n"//&
       "does not match one of ", cszj, cszj+1, dszj, dszj+1
    call MOM_error(FATAL,"post_data_2d_low: "//trim(diag%debug_str)//trim(mesg))
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
  if (diag_cs%diag_as_chksum) then
    chksum = chksum_general(locfield)
    if (is_root_pe()) then
      call log_chksum_diag(diag_cs%chksum_diag_doc_unit, diag%debug_str, chksum)
    endif
  else
    if (is_stat) then
      if (present(mask)) then
        call assert(size(locfield) == size(mask), &
            'post_data_2d_low is_stat: mask size mismatch: '//diag%debug_str)
        used = send_data(diag%fms_diag_id, locfield, &
                         is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, rmask=mask)
     !elseif (associated(diag%axes%mask2d)) then
     !  used = send_data(diag%fms_diag_id, locfield, &
     !                   is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, rmask=diag%axes%mask2d)
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
      elseif (associated(diag%axes%mask2d)) then
        used = send_data(diag%fms_diag_id, locfield, diag_cs%time_end, &
                         is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                         weight=diag_cs%time_int, rmask=diag%axes%mask2d)
      else
        used = send_data(diag%fms_diag_id, locfield, diag_cs%time_end, &
                         is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                         weight=diag_cs%time_int)
      endif
    endif
  endif
  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) &
    deallocate( locfield )

end subroutine post_data_2d_low

!> Make a real 3-d array diagnostic available for averaging or output.
subroutine post_data_3d(diag_field_id, field, diag_cs, is_static, mask, alt_h)

  integer,           intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                 !! previous call to register_diag_field.
  real,              intent(in) :: field(:,:,:)  !< 3-d array being offered for output or averaging
  type(diag_ctrl), target, intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.
  real,    optional, intent(in) :: mask(:,:,:) !< If present, use this real array as the data mask.
  real, dimension(:,:,:), &
         target, optional, intent(in) :: alt_h  !< An alternate thickness to use for vertically
                                                !! remapping this diagnostic, in H.

  ! Local variables
  type(diag_type), pointer :: diag => null()
  integer :: nz, i, j, k
  real, dimension(:,:,:), allocatable :: remapped_field
  logical :: staggered_in_x, staggered_in_y
  real, dimension(:,:,:), pointer :: h_diag => NULL()

  if (present(alt_h)) then
    h_diag => alt_h
  else
    h_diag => diag_cs%h
  endif

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
              diag_cs%G, h_diag, staggered_in_x, staggered_in_y, &
              diag%axes%mask3d, diag_cs%missing_value, field, remapped_field)
      if (id_clock_diag_remap>0) call cpu_clock_end(id_clock_diag_remap)
      if (associated(diag%axes%mask3d)) then
        ! Since 3d masks do not vary in the vertical, just use as much as is
        ! needed.
        call post_data_3d_low(diag, remapped_field, diag_cs, is_static, &
                              mask=diag%axes%mask3d)
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
              diag_cs%G, diag_cs%GV, h_diag, staggered_in_x, staggered_in_y, &
              diag%axes%mask3d, diag_cs%missing_value, field, remapped_field)
      if (id_clock_diag_remap>0) call cpu_clock_end(id_clock_diag_remap)
      if (associated(diag%axes%mask3d)) then
        ! Since 3d masks do not vary in the vertical, just use as much as is
        ! needed.
        call post_data_3d_low(diag, remapped_field, diag_cs, is_static, &
                              mask=diag%axes%mask3d)
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
              diag_cs%G, h_diag, staggered_in_x, staggered_in_y, &
              diag%axes%mask3d, diag_cs%missing_value, field, remapped_field)
      if (id_clock_diag_remap>0) call cpu_clock_end(id_clock_diag_remap)
      if (associated(diag%axes%mask3d)) then
        ! Since 3d masks do not vary in the vertical, just use as much as is
        ! needed.
        call post_data_3d_low(diag, remapped_field, diag_cs, is_static, &
                              mask=diag%axes%mask3d)
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

!> Make a real 3-d array diagnostic available for averaging or output
!! using a diag_type instead of an integer id.
subroutine post_data_3d_low(diag, field, diag_cs, is_static, mask)
  type(diag_type),   intent(in) :: diag       !< A structure describing the diagnostic to post
  real,    target,   intent(in) :: field(:,:,:) !< 3-d array being offered for output or averaging
  type(diag_ctrl),   intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.
  real,    optional, intent(in) :: mask(:,:,:) !< If present, use this real array as the data mask.

  ! Local variables
  real, dimension(:,:,:), pointer :: locfield => NULL()
  character(len=300) :: mesg
  logical :: used  ! The return value of send_data is not used for anything.
  logical :: staggered_in_x, staggered_in_y
  logical :: is_stat
  integer :: cszi, cszj, dszi, dszj
  integer :: isv, iev, jsv, jev, ks, ke, i, j, k, isv_c, jsv_c
  integer :: chksum

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
    isv = diag_cs%is ; iev = diag_cs%ie     ! Data domain
  elseif ( size(field,1) == dszi + 1 ) then
    isv = diag_cs%is ; iev = diag_cs%ie+1   ! Symmetric data domain
  elseif ( size(field,1) == cszi) then
    isv = 1 ; iev = cszi                    ! Computational domain
  elseif ( size(field,1) == cszi + 1 ) then
    isv = 1 ; iev = cszi+1                  ! Symmetric computational domain
  else
    write (mesg,*) " peculiar size ",size(field,1)," in i-direction\n"//&
       "does not match one of ", cszi, cszi+1, dszi, dszi+1
    call MOM_error(FATAL,"post_data_3d_low: "//trim(diag%debug_str)//trim(mesg))
  endif

  if ( size(field,2) == dszj ) then
    jsv = diag_cs%js ; jev = diag_cs%je     ! Data domain
  elseif ( size(field,2) == dszj + 1 ) then
    jsv = diag_cs%js ; jev = diag_cs%je+1   ! Symmetric data domain
  elseif ( size(field,2) == cszj ) then
    jsv = 1 ; jev = cszj                    ! Computational domain
  elseif ( size(field,2) == cszj+1 ) then
    jsv = 1 ; jev = cszj+1                  ! Symmetric computational domain
  else
    write (mesg,*) " peculiar size ",size(field,2)," in j-direction\n"//&
       "does not match one of ", cszj, cszj+1, dszj, dszj+1
    call MOM_error(FATAL,"post_data_3d_low: "//trim(diag%debug_str)//trim(mesg))
  endif

  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) then
    ks = lbound(field,3) ; ke = ubound(field,3)
    allocate( locfield( lbound(field,1):ubound(field,1), lbound(field,2):ubound(field,2), ks:ke ) )
    ! locfield(:,:,:) = 0.0  ! Zeroing out this array would be a good idea, but it appears
                             ! not to be necessary.
    isv_c = isv ; jsv_c = jsv
    if (diag%fms_xyave_diag_id>0) then
      staggered_in_x = diag%axes%is_u_point .or. diag%axes%is_q_point
      staggered_in_y = diag%axes%is_v_point .or. diag%axes%is_q_point
      ! When averaging a staggered field, edge points are always required.
      if (staggered_in_x) isv_c = iev - (diag_cs%ie - diag_cs%is) - 1
      if (staggered_in_y) jsv_c = jev - (diag_cs%je - diag_cs%js) - 1
      if (isv_c < lbound(locfield,1)) call MOM_error(FATAL, &
        "It is an error to average a staggered diagnostic field that does not "//&
        "have i-direction space to represent the symmetric computational domain.")
      if (jsv_c < lbound(locfield,2)) call MOM_error(FATAL, &
        "It is an error to average a staggered diagnostic field that does not "//&
        "have j-direction space to represent the symmetric computational domain.")
    endif

    do k=ks,ke ; do j=jsv_c,jev ; do i=isv_c,iev
      if (field(i,j,k) == diag_cs%missing_value) then
        locfield(i,j,k) = diag_cs%missing_value
      else
        locfield(i,j,k) = field(i,j,k) * diag%conversion_factor
      endif
    enddo ; enddo ; enddo
  else
    locfield => field
  endif

  if (diag%fms_diag_id>0) then
    if (diag_cs%diag_as_chksum) then
      chksum = chksum_general(locfield)
      if (is_root_pe()) then
        call log_chksum_diag(diag_cs%chksum_diag_doc_unit, diag%debug_str, chksum)
      endif
    else
      if (is_stat) then
        if (present(mask)) then
          call assert(size(locfield) == size(mask), &
              'post_data_3d_low is_stat: mask size mismatch: '//diag%debug_str)
          used = send_data(diag%fms_diag_id, locfield, &
                           is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, rmask=mask)
       !elseif (associated(diag%axes%mask3d)) then
       !  used = send_data(diag_field_id, locfield, &
       !                   is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, rmask=diag%axes%mask3d)
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
        elseif (associated(diag%axes%mask3d)) then
          call assert(size(locfield) == size(diag%axes%mask3d), &
              'post_data_3d_low: mask3d size mismatch: '//diag%debug_str)
          used = send_data(diag%fms_diag_id, locfield, diag_cs%time_end, &
                           is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                           weight=diag_cs%time_int, rmask=diag%axes%mask3d)
        else
          used = send_data(diag%fms_diag_id, locfield, diag_cs%time_end, &
                           is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                           weight=diag_cs%time_int)
        endif
      endif
    endif
  endif
  if (diag%fms_xyave_diag_id>0) then
    call post_xy_average(diag_cs, diag, locfield)
  endif
  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) &
    deallocate( locfield )

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

!> This subroutine enables the accumulation of time averages over the specified time interval.
subroutine enable_averaging(time_int_in, time_end_in, diag_cs)
  real,            intent(in)    :: time_int_in !< The time interval in s over which any
                                                !!  values that are offered are valid.
  type(time_type), intent(in)    :: time_end_in !< The end time of the valid interval
  type(diag_ctrl), intent(inout) :: diag_CS !< Structure used to regulate diagnostic output

! This subroutine enables the accumulation of time averages over the
! specified time interval.

!  if (num_file==0) return
  diag_cs%time_int = time_int_in
  diag_cs%time_end = time_end_in
  diag_cs%ave_enabled = .true.
end subroutine enable_averaging

!> Call this subroutine to avoid averaging any offered fields.
subroutine disable_averaging(diag_cs)
  type(diag_ctrl), intent(inout) :: diag_CS !< Structure used to regulate diagnostic output

  diag_cs%time_int = 0.0
  diag_cs%ave_enabled = .false.

end subroutine disable_averaging

!> Call this subroutine to determine whether the averaging is
!! currently enabled.  .true. is returned if it is.
function query_averaging_enabled(diag_cs, time_int, time_end)
  type(diag_ctrl),           intent(in)  :: diag_CS  !< Structure used to regulate diagnostic output
  real,            optional, intent(out) :: time_int !< Current setting of diag%time_int, in s
  type(time_type), optional, intent(out) :: time_end !< Current setting of diag%time_end
  logical :: query_averaging_enabled

  if (present(time_int)) time_int = diag_cs%time_int
  if (present(time_end)) time_end = diag_cs%time_end
  query_averaging_enabled = diag_cs%ave_enabled
end function query_averaging_enabled

!> This function returns the valid end time for use with diagnostics that are
!! handled outside of the MOM6 diagnostics infrastructure.
function get_diag_time_end(diag_cs)
  type(diag_ctrl), intent(in)  :: diag_CS !< Structure used to regulate diagnostic output
  type(time_type) :: get_diag_time_end
  !   This function returns the valid end time for diagnostics that are handled
  ! outside of the MOM6 infrastructure, such as via the generic tracer code.

  get_diag_time_end = diag_cs%time_end
end function get_diag_time_end

!> Returns the "diag_mediator" handle for a group (native, CMOR, z-coord, ...) of diagnostics
!! derived from one field.
integer function register_diag_field(module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, mask_variant, standard_name,      &
     verbose, do_not_log, err_msg, interp_method, tile_count, cmor_field_name, &
     cmor_long_name, cmor_units, cmor_standard_name, cell_methods, &
     x_cell_method, y_cell_method, v_cell_method, conversion, v_extensive)
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ocean_model"
                                              !! or "ice_shelf_model"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp), target, intent(in) :: axes !< Container w/ up to 3 integer handles that
                                             !! indicates axes for this field
  type(time_type),  intent(in) :: init_time !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values.
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
  logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided with
                                                         !! post_data calls (not used in MOM?)
  logical,          optional, intent(in) :: verbose !< If true, FMS is verbose (not used in MOM?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be
                                                         !! placed (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not
                                                         !! be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count !< no clue (not used in MOM?)
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
  character(len=*), optional, intent(in) :: v_cell_method !< Specifies the cell method for the vertical direction.
                                                         !! Use '' have no method.
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to file
  logical,          optional, intent(in) :: v_extensive !< True for vertically extensive fields (vertically
                                                         !! integrated). Default/absent for intensive.
  ! Local variables
  real :: MOM_missing_value
  type(diag_ctrl), pointer :: diag_cs => NULL()
  type(axes_grp), pointer :: remap_axes => null()
  integer :: dm_id, i
  character(len=256) :: new_module_name
  logical :: active

  MOM_missing_value = axes%diag_cs%missing_value
  if (present(missing_value)) MOM_missing_value = missing_value

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
      if ((axes%id == diag_cs%axesTL%id)) then
          remap_axes => diag_cs%remap_axesTL(i)
      elseif (axes%id == diag_cs%axesBL%id) then
          remap_axes => diag_cs%remap_axesBL(i)
      elseif (axes%id == diag_cs%axesCuL%id ) then
          remap_axes => diag_cs%remap_axesCuL(i)
      elseif (axes%id == diag_cs%axesCvL%id) then
          remap_axes => diag_cs%remap_axesCvL(i)
      elseif (axes%id == diag_cs%axesTi%id) then
          remap_axes => diag_cs%remap_axesTi(i)
      elseif (axes%id == diag_cs%axesBi%id) then
          remap_axes => diag_cs%remap_axesBi(i)
      elseif (axes%id == diag_cs%axesCui%id ) then
          remap_axes => diag_cs%remap_axesCui(i)
      elseif (axes%id == diag_cs%axesCvi%id) then
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
  type(axes_grp), target, intent(in) :: axes !< Container w/ up to 3 integer handles that indicates axes
                                             !! for this field
  type(time_type),  intent(in) :: init_time !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values.
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
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
  character(len=*), optional, intent(in) :: v_cell_method !< Specifies the cell method for the vertical direction.
                                                         !! Use '' have no method.
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to file
  logical,          optional, intent(in) :: v_extensive !< True for vertically extensive fields (vertically
                                                         !! integrated). Default/absent for intensive.
  ! Local variables
  real :: MOM_missing_value
  type(diag_ctrl), pointer :: diag_cs => null()
  type(diag_type), pointer :: this_diag => null()
  integer :: fms_id, fms_xyave_id
  character(len=256) :: posted_cmor_units, posted_cmor_standard_name, posted_cmor_long_name, cm_string, msg

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
  call attach_cell_methods(fms_id, axes, cm_string, &
                           cell_methods, x_cell_method, y_cell_method, v_cell_method, &
                           v_extensive=v_extensive)
  if (is_root_pe() .and. diag_CS%available_diag_doc_unit > 0) then
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
    if (is_root_pe() .and. diag_CS%available_diag_doc_unit > 0) then
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
    posted_cmor_units = "not provided"         !
    posted_cmor_standard_name = "not provided" ! Values might be able to be replaced with a CS%missing field?
    posted_cmor_long_name = "not provided"     !

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
               long_name=trim(posted_cmor_long_name), units=trim(posted_cmor_units),                  &
               missing_value=MOM_missing_value, range=range, mask_variant=mask_variant,               &
               standard_name=trim(posted_cmor_standard_name), verbose=verbose, do_not_log=do_not_log, &
               err_msg=err_msg, interp_method=interp_method, tile_count=tile_count)
    call attach_cell_methods(fms_id, axes, cm_string, &
                             cell_methods, x_cell_method, y_cell_method, v_cell_method, &
                             v_extensive=v_extensive)
    if (is_root_pe() .and. diag_CS%available_diag_doc_unit > 0) then
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
               long_name=trim(posted_cmor_long_name), units=trim(posted_cmor_units),                  &
               missing_value=MOM_missing_value, range=range, mask_variant=mask_variant,               &
               standard_name=trim(posted_cmor_standard_name), verbose=verbose, do_not_log=do_not_log, &
               err_msg=err_msg, interp_method=interp_method, tile_count=tile_count)
      call attach_cell_methods(fms_xyave_id, axes%xyave_axes, cm_string, &
                               cell_methods, v_cell_method, v_extensive=v_extensive)
      if (is_root_pe() .and. diag_CS%available_diag_doc_unit > 0) then
        msg = 'native name is "'//trim(field_name)//'_xyave"'
        call log_available_diag(fms_xyave_id>0, module_name, trim(cmor_field_name)//'_xyave', &
                                cm_string, msg, diag_CS, posted_cmor_long_name, posted_cmor_units, &
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

!> Returns an FMS id from register_diag_field_fms (the diag_manager routine) after expanding axes
!! (axes-group) into handles and conditionally adding an FMS area_id for cell_measures.
integer function register_diag_field_expand_axes(module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, mask_variant, standard_name,  &
     verbose, do_not_log, err_msg, interp_method, tile_count)
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ocean_model"
                                              !! or "ice_shelf_model"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp), target, intent(in) :: axes !< Container w/ up to 3 integer handles that indicates
                                             !! axes for this field
  type(time_type),  intent(in) :: init_time !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values.
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
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
  integer :: fms_id, area_id, volume_id

  ! This gets the cell area associated with the grid location of this variable
  area_id = axes%id_area
  volume_id = axes%id_volume

  ! Get the FMS diagnostic id
  if (present(interp_method) .or. axes%is_h_point) then
    ! If interp_method is provided we must use it
    if (area_id>0) then
      if (volume_id>0) then
        fms_id = register_diag_field_fms(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method=interp_method, tile_count=tile_count, area=area_id, volume=volume_id)
      else
        fms_id = register_diag_field_fms(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method=interp_method, tile_count=tile_count, area=area_id)
      endif
    else
      if (volume_id>0) then
        fms_id = register_diag_field_fms(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method=interp_method, tile_count=tile_count, volume=volume_id)
      else
        fms_id = register_diag_field_fms(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method=interp_method, tile_count=tile_count)
      endif
    endif
  else
    ! If interp_method is not provided and the field is not at an h-point then interp_method='none'
    if (area_id>0) then
      if (volume_id>0) then
        fms_id = register_diag_field_fms(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method='none', tile_count=tile_count, area=area_id, volume=volume_id)
      else
        fms_id = register_diag_field_fms(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method='none', tile_count=tile_count, area=area_id)
      endif
    else
      if (volume_id>0) then
        fms_id = register_diag_field_fms(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method='none', tile_count=tile_count, volume=volume_id)
      else
        fms_id = register_diag_field_fms(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method='none', tile_count=tile_count)
      endif
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
  type(axes_grp), target, intent(in)    :: axes !< Container w/ up to 3 integer handles that
                                                !! indicates axes for this field
  character(len=*),       intent(in)    :: module_name !< Name of this module, usually
                                                       !! "ocean_model" or "ice_shelf_model"
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
  this_diag%axes => axes

end subroutine add_diag_to_list

!> Attaches "cell_methods" attribute to a variable based on defaults for axes_grp or optional arguments.
subroutine attach_cell_methods(id, axes, ostring, cell_methods, &
                               x_cell_method, y_cell_method, v_cell_method, v_extensive)
  integer,                    intent(in)  :: id !< Handle to diagnostic
  type(axes_grp),             intent(in)  :: axes !< Container w/ up to 3 integer handles that indicates
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
  character(len=*), optional, intent(in)  :: v_cell_method !< Specifies the cell method for the vertical direction.
                                                         !! Use '' have no method.
  logical,          optional, intent(in)  :: v_extensive !< True for vertically extensive fields
                                                         !! (vertically integrated). Default/absent for intensive.
  ! Local variables
  character(len=9) :: axis_name
  logical :: x_mean, y_mean, x_sum, y_sum

  x_mean = .false.
  y_mean = .false.
  x_sum = .false.
  y_sum = .false.

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
        if (trim(x_cell_method)=='mean') x_mean=.true.
        if (trim(x_cell_method)=='sum') x_sum=.true.
      endif
    else
      if (len(trim(axes%x_cell_method))>0) then
        call get_diag_axis_name(axes%handles(1), axis_name)
        call diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(axes%x_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(axes%x_cell_method)
        if (trim(axes%x_cell_method)=='mean') x_mean=.true.
        if (trim(axes%x_cell_method)=='sum') x_sum=.true.
      endif
    endif
    if (present(y_cell_method)) then
      if (len(trim(y_cell_method))>0) then
        call get_diag_axis_name(axes%handles(2), axis_name)
        call diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(y_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(y_cell_method)
        if (trim(y_cell_method)=='mean') y_mean=.true.
        if (trim(y_cell_method)=='sum') y_sum=.true.
      endif
    else
      if (len(trim(axes%y_cell_method))>0) then
        call get_diag_axis_name(axes%handles(2), axis_name)
        call diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(axes%y_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(axes%y_cell_method)
        if (trim(axes%y_cell_method)=='mean') y_mean=.true.
        if (trim(axes%y_cell_method)=='sum') y_sum=.true.
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
    if (x_mean .and. y_mean) then
      call diag_field_add_attribute(id, 'cell_methods', 'area:mean')
      ostring = trim(adjustl(ostring))//' area:mean'
    elseif (x_sum .and. y_sum) then
      call diag_field_add_attribute(id, 'cell_methods', 'area:sum')
      ostring = trim(adjustl(ostring))//' area:sum'
    endif
  endif
  ostring = adjustl(ostring)
end subroutine attach_cell_methods

function register_scalar_field(module_name, field_name, init_time, diag_cs, &
     long_name, units, missing_value, range, standard_name, &
     do_not_log, err_msg, interp_method, cmor_field_name, &
     cmor_long_name, cmor_units, cmor_standard_name)
  integer :: register_scalar_field !< An integer handle for a diagnostic array.
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ocean_model"
                                              !! or "ice_shelf_model"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(time_type),  intent(in) :: init_time !< Time at which a field is first available?
  type(diag_ctrl),  intent(inout) :: diag_CS !< Structure used to regulate diagnostic output
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values.
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be
                                                         !! placed (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not
                                                         !! be interpolated as a scalar
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
  character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
  character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
  character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field

  ! Local variables
  real :: MOM_missing_value
  integer :: dm_id, fms_id
  type(diag_type), pointer :: diag => null(), cmor_diag => null()
  character(len=256) :: posted_cmor_units, posted_cmor_standard_name, posted_cmor_long_name

  MOM_missing_value = diag_cs%missing_value
  if (present(missing_value)) MOM_missing_value = missing_value

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
  if (is_root_pe() .and. diag_CS%available_diag_doc_unit > 0) then
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
     cmor_field_name, cmor_long_name, cmor_units, cmor_standard_name, area, &
     x_cell_method, y_cell_method, area_cell_method)
  integer :: register_static_field !< An integer handle for a diagnostic array.
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ocean_model"
                                              !! or "ice_shelf_model"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp), target, intent(in) :: axes !< Container w/ up to 3 integer handles that
                                             !! indicates axes for this field
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values.
  real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
  logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided with
                                                         !! post_data calls (not used in MOM?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not
                                                         !! be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count !< no clue (not used in MOM?)
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
  character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
  character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
  character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
  integer,          optional, intent(in) :: area !< fms_id for area_t
  character(len=*), optional, intent(in) :: x_cell_method !< Specifies the cell method for the x-direction.
  character(len=*), optional, intent(in) :: y_cell_method !< Specifies the cell method for the y-direction.
  character(len=*), optional, intent(in) :: area_cell_method !< Specifies the cell method for area

  ! Local variables
  real :: MOM_missing_value
  type(diag_ctrl), pointer :: diag_cs => null()
  type(diag_type), pointer :: diag => null(), cmor_diag => null()
  integer :: dm_id, fms_id, cmor_id
  character(len=256) :: posted_cmor_units, posted_cmor_standard_name, posted_cmor_long_name
  character(len=9) :: axis_name

  MOM_missing_value = axes%diag_cs%missing_value
  if (present(missing_value)) MOM_missing_value = missing_value

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
    if (present(x_cell_method)) then
      call get_diag_axis_name(axes%handles(1), axis_name)
      call diag_field_add_attribute(fms_id, 'cell_methods', trim(axis_name)//':'//trim(x_cell_method))
    endif
    if (present(y_cell_method)) then
      call get_diag_axis_name(axes%handles(2), axis_name)
      call diag_field_add_attribute(fms_id, 'cell_methods', trim(axis_name)//':'//trim(y_cell_method))
    endif
    if (present(area_cell_method)) then
      call diag_field_add_attribute(fms_id, 'cell_methods', 'area:'//trim(area_cell_method))
    endif
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
      if (present(x_cell_method)) then
        call get_diag_axis_name(axes%handles(1), axis_name)
        call diag_field_add_attribute(fms_id, 'cell_methods', trim(axis_name)//':'//trim(x_cell_method))
      endif
      if (present(y_cell_method)) then
        call get_diag_axis_name(axes%handles(2), axis_name)
        call diag_field_add_attribute(fms_id, 'cell_methods', trim(axis_name)//':'//trim(y_cell_method))
      endif
      if (present(area_cell_method)) then
        call diag_field_add_attribute(fms_id, 'cell_methods', 'area:'//trim(area_cell_method))
      endif
    endif
  endif

  ! Document diagnostics in list of available diagnostics
  if (is_root_pe() .and. diag_CS%available_diag_doc_unit > 0) then
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

!> Describe an option setting in the diagnostic files.
subroutine describe_option(opt_name, value, diag_CS)
  character(len=*), intent(in) :: opt_name !< The name of the option
  character(len=*), intent(in) :: value   !< A character string with the setting of the option.
  type(diag_ctrl),  intent(in) :: diag_CS !< Structure used to regulate diagnostic output

  character(len=240) :: mesg
  integer :: len_ind

  len_ind = len_trim(value)  ! Add error handling for long values?

  mesg = "    ! "//trim(opt_name)//": "//trim(value)
  write(diag_CS%available_diag_doc_unit, '(a)') trim(mesg)
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
  type(axes_grp), pointer :: axes => NULL()

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
          axes, day, trim(longname), trim(units), missing_value=-1.0e+34)

end function ocean_register_diag

subroutine diag_mediator_infrastructure_init(err_msg)
  ! This subroutine initializes the FMS diag_manager.
  character(len=*), optional, intent(out)   :: err_msg !< An error message

  call diag_manager_init(err_msg=err_msg)
end subroutine diag_mediator_infrastructure_init

!> diag_mediator_init initializes the MOM diag_mediator and opens the available
!! diagnostics file, if appropriate.
subroutine diag_mediator_init(G, GV, nz, param_file, diag_cs, doc_file_dir)
  type(ocean_grid_type), target, intent(inout) :: G  !< The ocean grid type.
  type(verticalGrid_type), target, intent(in)  :: GV !< The ocean vertical grid structure
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
  character(len=40) :: mdl = "MOM_diag_mediator" ! This module's name.

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
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, 'NUM_DIAG_COORDS', diag_cs%num_diag_coords, &
                 'The number of diagnostic vertical coordinates to use.\n'//&
                 'For each coordinate, an entry in DIAG_COORDS must be provided.', &
                 default=1)
  if (diag_cs%num_diag_coords>0) then
    allocate(diag_coords(diag_cs%num_diag_coords))
    if (diag_cs%num_diag_coords==1) then ! The default is to provide just one instance of Z*
      call get_param(param_file, mdl, 'DIAG_COORDS', diag_coords, &
                 'A list of string tuples associating diag_table modules to\n'//&
                 'a coordinate definition used for diagnostics. Each string\n'//&
                 'is of the form "MODULE_SUFFIX PARAMETER_SUFFIX COORDINATE_NAME".', &
                 default='z Z ZSTAR')
    else ! If using more than 1 diagnostic coordinate, all must be explicitly defined
      call get_param(param_file, mdl, 'DIAG_COORDS', diag_coords, &
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

  call get_param(param_file, mdl, 'DIAG_MISVAL', diag_cs%missing_value, &
                 'Set the default missing value to use for diagnostics.', &
                 default=1.e20)
  call get_param(param_file, mdl, 'DIAG_AS_CHKSUM', diag_cs%diag_as_chksum, &
                 'Instead of writing diagnostics to the diag manager, write\n' //&
                 'a textfile containing the checksum (bitcount) of the array.',  &
                 default=.false.)

  ! Keep pointers grid, h, T, S needed diagnostic remapping
  diag_cs%G => G
  diag_cs%GV => GV
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

  ! Initialze available diagnostic log file
  if (is_root_pe() .and. (diag_CS%available_diag_doc_unit < 0)) then
    write(this_pe,'(i6.6)') PE_here()
    doc_file_dflt = "available_diags."//this_pe
    call get_param(param_file, mdl, "AVAILABLE_DIAGS_FILE", doc_file, &
                 "A file into which to write a list of all available \n"//&
                 "ocean diagnostics that can be included in a diag_table.", &
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

  if (is_root_pe() .and. (diag_CS%chksum_diag_doc_unit < 0) .and. diag_CS%diag_as_chksum) then
    write(this_pe,'(i6.6)') PE_here()
    doc_file_dflt = "chksum_diag."//this_pe
    call get_param(param_file, mdl, "CHKSUM_DIAG_FILE", doc_file, &
                 "A file into which to write all checksums of the \n"//&
                 "diagnostics listed in the diag_table.", &
                 default=doc_file_dflt, do_not_log=(diag_CS%chksum_diag_doc_unit/=-1))
    if (len_trim(doc_file) > 0) then
      new_file = .true. ; if (diag_CS%chksum_diag_doc_unit /= -1) new_file = .false.
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

      diag_CS%chksum_diag_doc_unit = new_unit

      if (new_file) then
        open(diag_CS%chksum_diag_doc_unit, file=trim(doc_path), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='REPLACE', iostat=ios)
      else ! This file is being reopened, and should be appended.
        open(diag_CS%chksum_diag_doc_unit, file=trim(doc_path), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='OLD', position='APPEND', iostat=ios)
      endif
      inquire(diag_CS%chksum_diag_doc_unit, opened=opened)
      if ((.not.opened) .or. (ios /= 0)) then
        call MOM_error(FATAL, "Failed to open checksum diags file "//trim(doc_path)//".")
      endif
    endif
  endif

end subroutine diag_mediator_init

!> Set pointers to the default state fields used to remap diagnostics.
subroutine diag_set_state_ptrs(h, T, S, eqn_of_state, diag_cs)
  real, dimension(:,:,:), target, intent(in   ) :: h !< the model thickness array
  real, dimension(:,:,:), target, intent(in   ) :: T !< the model temperature array
  real, dimension(:,:,:), target, intent(in   ) :: S !< the model salinity array
  type(EOS_type),         target, intent(in   ) :: eqn_of_state !< Equation of state structure
  type(diag_ctrl),                intent(inout) :: diag_cs !< diag mediator control structure

  ! Keep pointers to h, T, S needed for the diagnostic remapping
  diag_cs%h => h
  diag_cs%T => T
  diag_cs%S => S
  diag_cs%eqn_of_state => eqn_of_state

end subroutine

!> Build/update vertical grids for diagnostic remapping.
!! \note The target grids need to be updated whenever sea surface
!! height changes.
subroutine diag_update_remap_grids(diag_cs, alt_h, alt_T, alt_S)
  type(diag_ctrl),        intent(inout) :: diag_cs      !< Diagnostics control structure
  real, target, optional, intent(in   ) :: alt_h(:,:,:) !< Used if remapped grids should be something other than
                                                        !! the current thicknesses
  real, target, optional, intent(in   ) :: alt_T(:,:,:) !< Used if remapped grids should be something other than
                                                        !! the current temperatures
  real, target, optional, intent(in   ) :: alt_S(:,:,:) !< Used if remapped grids should be something other than
                                                        !! the current salinity
  ! Local variables
  integer :: i
  real, dimension(:,:,:), pointer :: h_diag => NULL()
  real, dimension(:,:,:), pointer :: T_diag => NULL(), S_diag => NULL()

  if (present(alt_h)) then
    h_diag => alt_h
  else
    h_diag => diag_cs%h
  endif

  if (present(alt_T)) then
    T_diag => alt_T
  else
    T_diag => diag_CS%T
  endif

  if (present(alt_S)) then
    S_diag => alt_S
  else
    S_diag => diag_CS%S
  endif

  if (id_clock_diag_grid_updates>0) call cpu_clock_begin(id_clock_diag_grid_updates)

  if (diag_cs%diag_grid_overridden) then
     call MOM_error(FATAL, "diag_update_remap_grids was called, but current grids in "// &
                           "diagnostic structure have been overridden")
  endif

  do i=1, diag_cs%num_diag_coords
    call diag_remap_update(diag_cs%diag_remap_cs(i), &
                           diag_cs%G, diag_cs%GV, h_diag, T_diag, S_diag, &
                           diag_cs%eqn_of_state)
  enddo

#if defined(DEBUG) || defined(__DO_SAFETY_CHECKS__)
  ! Keep a copy of H - used to check whether grids are up-to-date
  ! when doing remapping.
  diag_cs%h_old(:,:,:) = diag_cs%h(:,:,:)
#endif

  if (id_clock_diag_grid_updates>0) call cpu_clock_end(id_clock_diag_grid_updates)

end subroutine diag_update_remap_grids

!> Sets up the 2d and 3d masks for native diagnostics
subroutine diag_masks_set(G, nz, diag_cs)
  type(ocean_grid_type), target, intent(in) :: G  !< The ocean grid type.
  integer,                       intent(in) :: nz !< The number of layers in the model's native grid.
  type(diag_ctrl),               pointer    :: diag_cs !< A pointer to a type with many variables
                                                       !! used for diagnostics
  ! Local variables
  integer :: k

  ! 2d masks point to the model masks since they are identical
  diag_cs%mask2dT  => G%mask2dT
  diag_cs%mask2dBu => G%mask2dBu
  diag_cs%mask2dCu => G%mask2dCu
  diag_cs%mask2dCv => G%mask2dCv

  ! 3d native masks are needed by diag_manager but the native variables
  ! can only be masked 2d - for ocean points, all layers exists.
  allocate(diag_cs%mask3dTL(G%isd:G%ied,G%jsd:G%jed,1:nz))
  allocate(diag_cs%mask3dBL(G%IsdB:G%IedB,G%JsdB:G%JedB,1:nz))
  allocate(diag_cs%mask3dCuL(G%IsdB:G%IedB,G%jsd:G%jed,1:nz))
  allocate(diag_cs%mask3dCvL(G%isd:G%ied,G%JsdB:G%JedB,1:nz))
  do k=1,nz
    diag_cs%mask3dTL(:,:,k) = diag_cs%mask2dT(:,:)
    diag_cs%mask3dBL(:,:,k) = diag_cs%mask2dBu(:,:)
    diag_cs%mask3dCuL(:,:,k) = diag_cs%mask2dCu(:,:)
    diag_cs%mask3dCvL(:,:,k) = diag_cs%mask2dCv(:,:)
  enddo
  allocate(diag_cs%mask3dTi(G%isd:G%ied,G%jsd:G%jed,1:nz+1))
  allocate(diag_cs%mask3dBi(G%IsdB:G%IedB,G%JsdB:G%JedB,1:nz+1))
  allocate(diag_cs%mask3dCui(G%IsdB:G%IedB,G%jsd:G%jed,1:nz+1))
  allocate(diag_cs%mask3dCvi(G%isd:G%ied,G%JsdB:G%JedB,1:nz+1))
  do k=1,nz+1
    diag_cs%mask3dTi(:,:,k) = diag_cs%mask2dT(:,:)
    diag_cs%mask3dBi(:,:,k) = diag_cs%mask2dBu(:,:)
    diag_cs%mask3dCui(:,:,k) = diag_cs%mask2dCu(:,:)
    diag_cs%mask3dCvi(:,:,k) = diag_cs%mask2dCv(:,:)
  enddo

end subroutine diag_masks_set

subroutine diag_mediator_close_registration(diag_CS)
  type(diag_ctrl), intent(inout) :: diag_CS !< Structure used to regulate diagnostic output

  integer :: i

  if (diag_CS%available_diag_doc_unit > -1) then
    close(diag_CS%available_diag_doc_unit) ; diag_CS%available_diag_doc_unit = -2
  endif

  do i=1, diag_cs%num_diag_coords
    call diag_remap_diag_registration_closed(diag_cs%diag_remap_cs(i))
  enddo

end subroutine diag_mediator_close_registration

subroutine diag_mediator_end(time, diag_CS, end_diag_manager)
  type(time_type),   intent(in)  :: time !< The current model time
  type(diag_ctrl), intent(inout) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in)  :: end_diag_manager !< If true, call diag_manager_end()

  ! Local variables
  integer :: i

  if (diag_CS%available_diag_doc_unit > -1) then
    close(diag_CS%available_diag_doc_unit) ; diag_CS%available_diag_doc_unit = -3
  endif
  if (diag_CS%chksum_diag_doc_unit > -1) then
    close(diag_CS%chksum_diag_doc_unit) ; diag_CS%chksum_diag_doc_unit = -3
  endif

  deallocate(diag_cs%diags)

  do i=1, diag_cs%num_diag_coords
    call diag_remap_end(diag_cs%diag_remap_cs(i))
  enddo

  call diag_grid_storage_end(diag_cs%diag_grid_temp)
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

!> Convert the first n elements (up to 3) of an integer array to an underscore delimited string.
function i2s(a,n_in)
  ! "Convert the first n elements of an integer array to a string."
  ! Perhaps this belongs elsewhere in the MOM6 code?
  integer, dimension(:), intent(in) :: a    !< The array of integers to translate
  integer, optional    , intent(in) :: n_in !< The number of elements to translate, by default all
  character(len=15) :: i2s !< The returned string

  character(len=15) :: i2s_temp
  integer :: i,n

  n=size(a)
  if (present(n_in)) n = n_in

  i2s = ''
  do i=1,min(n,3)
    write (i2s_temp, '(I4.4)') a(i)
    i2s = trim(i2s) //'_'// trim(i2s_temp)
  enddo
  i2s = adjustl(i2s)
end function i2s

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
    write(diag_CS%available_diag_doc_unit, '(a,x,"(",a,")")') trim(mesg),trim(comment)
  else
    write(diag_CS%available_diag_doc_unit, '(a)') trim(mesg)
  endif
  if (present(long_name)) call describe_option("long_name", long_name, diag_CS)
  if (present(units)) call describe_option("units", units, diag_CS)
  if (present(standard_name)) &
    call describe_option("standard_name", standard_name, diag_CS)
  if (len(trim((cell_methods_string)))>0) &
    call describe_option("cell_methods", trim(cell_methods_string), diag_CS)

end subroutine log_available_diag

!> Log the diagnostic chksum to the chksum diag file
subroutine log_chksum_diag(docunit, description, chksum)
  integer,          intent(in) :: docunit     !< Handle of the log file
  character(len=*), intent(in) :: description !< Name of the diagnostic module
  integer,          intent(in) :: chksum      !< chksum of the diagnostic

  write(docunit, '(a,x,i9.8)') description, chksum
  flush(docunit)

end subroutine log_chksum_diag

!> Allocates fields necessary to store diagnostic remapping fields
subroutine diag_grid_storage_init(grid_storage, G, diag)
  type(diag_grid_storage), intent(inout) :: grid_storage !< Structure containing a snapshot of the target grids
  type(ocean_grid_type),   intent(in)    :: G           !< Horizontal grid
  type(diag_ctrl),         intent(in)    :: diag        !< Diagnostic control structure used as the contructor
                                                        !! template for this routine

  integer :: m, nz
  grid_storage%num_diag_coords = diag%num_diag_coords

  ! Don't do anything else if there are no remapped coordinates
  if (grid_storage%num_diag_coords < 1) return

  ! Allocate memory for the native space
  allocate(grid_storage%h_state(G%isd:G%ied,G%jsd:G%jed, G%ke))
  ! Allocate diagnostic remapping structures
  allocate(grid_storage%diag_grids(diag%num_diag_coords))
  ! Loop through and allocate memory for the grid on each target coordinate
  do m = 1, diag%num_diag_coords
    nz = diag%diag_remap_cs(m)%nz
    allocate(grid_storage%diag_grids(m)%h(G%isd:G%ied,G%jsd:G%jed, nz))
  enddo

end subroutine diag_grid_storage_init

!> Copy from the main diagnostic arrays to the grid storage as well as the native thicknesses
subroutine diag_copy_diag_to_storage(grid_storage, h_state, diag)
  type(diag_grid_storage), intent(inout) :: grid_storage !< Structure containing a snapshot of the target grids
  real, dimension(:,:,:),  intent(in)    :: h_state     !< Current model thicknesses
  type(diag_ctrl),         intent(in)    :: diag     !< Diagnostic control structure used as the contructor

  integer :: m

  ! Don't do anything else if there are no remapped coordinates
  if (grid_storage%num_diag_coords < 1) return

  grid_storage%h_state(:,:,:) = h_state(:,:,:)
  do m = 1,grid_storage%num_diag_coords
    if (diag%diag_remap_cs(m)%nz > 0) &
      grid_storage%diag_grids(m)%h(:,:,:) = diag%diag_remap_cs(m)%h(:,:,:)
  enddo

end subroutine diag_copy_diag_to_storage

!> Copy from the stored diagnostic arrays to the main diagnostic grids
subroutine diag_copy_storage_to_diag(diag, grid_storage)
  type(diag_ctrl),         intent(inout) :: diag     !< Diagnostic control structure used as the contructor
  type(diag_grid_storage), intent(in)    :: grid_storage !< Structure containing a snapshot of the target grids

  integer :: m

  ! Don't do anything else if there are no remapped coordinates
  if (grid_storage%num_diag_coords < 1) return

  diag%diag_grid_overridden = .true.
  do m = 1,grid_storage%num_diag_coords
    if (diag%diag_remap_cs(m)%nz > 0) &
      diag%diag_remap_cs(m)%h(:,:,:) = grid_storage%diag_grids(m)%h(:,:,:)
  enddo

end subroutine diag_copy_storage_to_diag

!> Save the current diagnostic grids in the temporary structure within diag
subroutine diag_save_grids(diag)
  type(diag_ctrl),         intent(inout) :: diag     !< Diagnostic control structure used as the contructor

  integer :: m

  ! Don't do anything else if there are no remapped coordinates
  if (diag%num_diag_coords < 1) return

  do m = 1,diag%num_diag_coords
    if (diag%diag_remap_cs(m)%nz > 0) &
      diag%diag_grid_temp%diag_grids(m)%h(:,:,:) = diag%diag_remap_cs(m)%h(:,:,:)
  enddo

end subroutine diag_save_grids

!> Restore the diagnostic grids from the temporary structure within diag
subroutine diag_restore_grids(diag)
  type(diag_ctrl),         intent(inout) :: diag     !< Diagnostic control structure used as the contructor

  integer :: m

  ! Don't do anything else if there are no remapped coordinates
  if (diag%num_diag_coords < 1) return

  diag%diag_grid_overridden = .false.
  do m = 1,diag%num_diag_coords
    if (diag%diag_remap_cs(m)%nz > 0) &
      diag%diag_remap_cs(m)%h(:,:,:) = diag%diag_grid_temp%diag_grids(m)%h(:,:,:)
  enddo

end subroutine diag_restore_grids

!> Deallocates the fields in the remapping fields container
subroutine diag_grid_storage_end(grid_storage)
  type(diag_grid_storage), intent(inout) :: grid_storage !< Structure containing a snapshot of the target grids
  ! Local variables
  integer :: m, nz

  ! Don't do anything else if there are no remapped coordinates
  if (grid_storage%num_diag_coords < 1) return

  ! Deallocate memory for the native space
  deallocate(grid_storage%h_state)
  ! Loop through and deallocate memory for the grid on each target coordinate
  do m = 1, grid_storage%num_diag_coords
    deallocate(grid_storage%diag_grids(m)%h)
  enddo
  ! Deallocate diagnostic remapping structures
  deallocate(grid_storage%diag_grids)
end subroutine diag_grid_storage_end

end module MOM_diag_mediator
