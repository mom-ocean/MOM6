!> The subroutines here provide convenient wrappers to the fms diag_manager
!! interfaces with additional diagnostic capabilies.
module MOM_diag_mediator

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_checksums,        only : chksum0, zchksum
use MOM_checksums,        only : hchksum, uchksum, vchksum, Bchksum
use MOM_coms,             only : PE_here
use MOM_cpu_clock,        only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,        only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_manager_infra, only : MOM_diag_manager_init, MOM_diag_manager_end
use MOM_diag_manager_infra, only : diag_axis_init=>MOM_diag_axis_init, get_MOM_diag_axis_name
use MOM_diag_manager_infra, only : send_data_infra, MOM_diag_field_add_attribute, EAST, NORTH
use MOM_diag_manager_infra, only : register_diag_field_infra, register_static_field_infra
use MOM_diag_manager_infra, only : get_MOM_diag_field_id, DIAG_FIELD_NOT_FOUND
use MOM_diag_remap,       only : diag_remap_ctrl, diag_remap_update, diag_remap_calc_hmask
use MOM_diag_remap,       only : diag_remap_init, diag_remap_end, diag_remap_do_remap
use MOM_diag_remap,       only : vertically_reintegrate_diag_field, vertically_interpolate_diag_field
use MOM_diag_remap,       only : horizontally_average_diag_field, diag_remap_get_axes_info
use MOM_diag_remap,       only : diag_remap_configure_axes, diag_remap_axes_configured
use MOM_diag_remap,       only : diag_remap_diag_registration_closed, diag_remap_set_active
use MOM_EOS,              only : EOS_type
use MOM_error_handler,    only : MOM_error, FATAL, WARNING, is_root_pe, assert
use MOM_file_parser,      only : get_param, log_version, param_file_type
use MOM_grid,             only : ocean_grid_type
use MOM_io,               only : slasher, vardesc, query_vardesc, MOM_read_data
use MOM_io,               only : get_filename_appendix
use MOM_safe_alloc,       only : safe_alloc_ptr, safe_alloc_alloc
use MOM_string_functions, only : lowercase
use MOM_time_manager,     only : time_type
use MOM_unit_scaling,     only : unit_scale_type
use MOM_verticalGrid,     only : verticalGrid_type

implicit none ; private

#undef __DO_SAFETY_CHECKS__
#define IMPLIES(A, B) ((.not. (A)) .or. (B))
#define MAX_DSAMP_LEV 2

public set_axes_info, post_data, register_diag_field, time_type
public set_masks_for_axes
public post_data_1d_k
public safe_alloc_ptr, safe_alloc_alloc
public enable_averaging, enable_averages, disable_averaging, query_averaging_enabled
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
public found_in_diagtable

!> Make a diagnostic available for averaging or output.
interface post_data
  module procedure post_data_3d, post_data_2d, post_data_1d_k, post_data_0d
end interface post_data

!> Down sample a field
interface downsample_field
  module procedure downsample_field_2d, downsample_field_3d
end interface downsample_field

!> Down sample the mask of a field
interface downsample_mask
  module procedure downsample_mask_2d, downsample_mask_3d
end interface downsample_mask

!> Down sample a diagnostic field
interface downsample_diag_field
  module procedure downsample_diag_field_2d, downsample_diag_field_3d
end interface downsample_diag_field

!> Contained for down sampled masks
type, private :: diag_dsamp
  real, pointer, dimension(:,:)   :: mask2d => null() !< Mask for 2d (x-y) axes
  real, pointer, dimension(:,:,:) :: mask3d => null() !< Mask for 3d axes
end type diag_dsamp

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
  integer :: downsample_level = 1 !< If greater than 1, the factor by which this diagnostic/axes/masks be downsampled
  ! For horizontally averaged diagnositcs (applies to 2d and 3d fields only)
  type(axes_grp), pointer :: xyave_axes => null() !< The associated 1d axes for horizontall area-averaged diagnostics
  ! ID's for cell_measures
  integer :: id_area = -1 !< The diag_manager id for area to be used for cell_measure of variables with this axes_grp.
  integer :: id_volume = -1 !< The diag_manager id for volume to be used for cell_measure of variables
                            !! with this axes_grp.
  ! For masking
  real, pointer, dimension(:,:)   :: mask2d => null() !< Mask for 2d (x-y) axes
  real, pointer, dimension(:,:,:) :: mask3d => null() !< Mask for 3d axes
  type(diag_dsamp), dimension(2:MAX_DSAMP_LEV) :: dsamp !< Downsample container
end type axes_grp

!> Contains an array to store a diagnostic target grid
type, private :: diag_grids_type
  real, dimension(:,:,:), allocatable :: h  !< Target grid for remapped coordinate
end type diag_grids_type

!> Stores all the remapping grids and the model's native space thicknesses
type, public :: diag_grid_storage
  integer                                          :: num_diag_coords !< Number of target coordinates
  real, dimension(:,:,:), allocatable              :: h_state         !< Layer thicknesses in native
                                                                      !! space [H ~> m or kg m-2]
  type(diag_grids_type), dimension(:), allocatable :: diag_grids      !< Primarily empty, except h field
end type diag_grid_storage

! Integers to encode the total cell methods
!integer :: PPP=111  ! x:point,y:point,z:point, this kind of diagnostic is not currently present in diag_table.MOM6
!integer :: PPS=112  ! x:point,y:point,z:sum  , this kind of diagnostic is not currently present in diag_table.MOM6
!integer :: PPM=113  ! x:point,y:point,z:mean , this kind of diagnostic is not currently present in diag_table.MOM6
integer :: PSP=121  !< x:point,y:sum,z:point
integer :: PSS=122  !< x:point,y:sum,z:point
integer :: PSM=123  !< x:point,y:sum,z:mean
integer :: PMP=131  !< x:point,y:mean,z:point
integer :: PMM=133  !< x:point,y:mean,z:mean
integer :: SPP=211  !< x:sum,y:point,z:point
integer :: SPS=212  !< x:sum,y:point,z:sum
integer :: SSP=221  !< x:sum;y:sum,z:point
integer :: MPP=311  !< x:mean,y:point,z:point
integer :: MPM=313  !< x:mean,y:point,z:mean
integer :: MMP=331  !< x:mean,y:mean,z:point
integer :: MMS=332  !< x:mean,y:mean,z:sum
integer :: SSS=222  !< x:sum,y:sum,z:sum
integer :: MMM=333  !< x:mean,y:mean,z:mean
integer :: MSK=-1   !< Use the downsample method of a mask

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
  integer :: downsample_diag_id = -1 !< For a horizontally area-downsampled diagnostic.
  character(64) :: debug_str = '' !< For FATAL errors and debugging.
  type(axes_grp), pointer :: axes => null() !< The axis group for this diagnostic
  type(diag_type), pointer :: next => null() !< Pointer to the next diagnostic
  real :: conversion_factor = 0. !< A factor to multiply data by before posting to FMS, if non-zero.
  logical :: v_extensive = .false. !< True for vertically extensive fields (vertically integrated).
                                   !! False for intensive (concentrations).
  integer :: xyz_method = 0 !< A 3 digit integer encoding the diagnostics cell method
                            !! It can be used to determine the downsample algorithm
end type diag_type

!> Container for down sampling information
type diagcs_dsamp
  integer :: isc !< The start i-index of cell centers within the computational domain
  integer :: iec !< The end i-index of cell centers within the computational domain
  integer :: jsc !< The start j-index of cell centers within the computational domain
  integer :: jec !< The end j-index of cell centers within the computational domain
  integer :: isd !< The start i-index of cell centers within the data domain
  integer :: ied !< The end i-index of cell centers within the data domain
  integer :: jsd !< The start j-index of cell centers within the data domain
  integer :: jed !< The end j-index of cell centers within the data domain
  integer :: isg !< The start i-index of cell centers within the global domain
  integer :: ieg !< The end i-index of cell centers within the global domain
  integer :: jsg !< The start j-index of cell centers within the global domain
  integer :: jeg !< The end j-index of cell centers within the global domain
  integer :: isgB !< The start i-index of cell corners within the global domain
  integer :: iegB !< The end i-index of cell corners within the global domain
  integer :: jsgB !< The start j-index of cell corners within the global domain
  integer :: jegB !< The end j-index of cell corners within the global domain

  !>@{ Axes for each location on a diagnostic grid
  type(axes_grp)  :: axesBL, axesTL, axesCuL, axesCvL
  type(axes_grp)  :: axesBi, axesTi, axesCui, axesCvi
  type(axes_grp)  :: axesB1, axesT1, axesCu1, axesCv1
  type(axes_grp), dimension(:), allocatable :: remap_axesTL, remap_axesBL, remap_axesCuL, remap_axesCvL
  type(axes_grp), dimension(:), allocatable :: remap_axesTi, remap_axesBi, remap_axesCui, remap_axesCvi
  !>@}

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
  !>@}
end type diagcs_dsamp

!> The following data type a list of diagnostic fields an their variants,
!! as well as variables that control the handling of model output.
type, public :: diag_ctrl
  integer :: available_diag_doc_unit = -1 !< The unit number of a diagnostic documentation file.
                                          !! This file is open if available_diag_doc_unit is > 0.
  integer :: chksum_iounit = -1           !< The unit number of a diagnostic documentation file.
                                          !! This file is open if available_diag_doc_unit is > 0.
  logical :: diag_as_chksum  !< If true, log chksums in a text file instead of posting diagnostics
  logical :: grid_space_axes !< If true, diagnostic horizontal coordinates axes are in grid space.
! The following fields are used for the output of the data.
  integer :: is  !< The start i-index of cell centers within the computational domain
  integer :: ie  !< The end i-index of cell centers within the computational domain
  integer :: js  !< The start j-index of cell centers within the computational domain
  integer :: je  !< The end j-index of cell centers within the computational domain

  integer :: isd !< The start i-index of cell centers within the data domain
  integer :: ied !< The end i-index of cell centers within the data domain
  integer :: jsd !< The start j-index of cell centers within the data domain
  integer :: jed !< The end j-index of cell centers within the data domain
  real :: time_int              !< The time interval for any fields
                                !! that are offered for averaging [s].
  type(time_type) :: time_end   !< The end time of the valid
                                !! interval for any offered field.
  logical :: ave_enabled = .false. !< True if averaging is enabled.

  !>@{ The following are 3D and 2D axis groups defined for output.  The names
  !! indicate the horizontal (B, T, Cu, or Cv) and vertical (L, i, or 1) locations.
  type(axes_grp) :: axesBL, axesTL, axesCuL, axesCvL
  type(axes_grp) :: axesBi, axesTi, axesCui, axesCvi
  type(axes_grp) :: axesB1, axesT1, axesCu1, axesCv1
  !>@}
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

  type(diagcs_dsamp), dimension(2:MAX_DSAMP_LEV) :: dsamp !< Downsample control container

  !>@}

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
  !>@{ Axes used for remapping
  type(axes_grp), dimension(:), allocatable :: remap_axesTL, remap_axesBL, remap_axesCuL, remap_axesCvL
  type(axes_grp), dimension(:), allocatable :: remap_axesTi, remap_axesBi, remap_axesCui, remap_axesCvi
  !>@}

  ! Pointer to H, G and T&S needed for remapping
  real, dimension(:,:,:), pointer :: h => null() !< The thicknesses needed for remapping [H ~> m or kg m-2]
  real, dimension(:,:,:), pointer :: T => null() !< The temperatures needed for remapping [degC]
  real, dimension(:,:,:), pointer :: S => null() !< The salinities needed for remapping [ppt]
  type(EOS_type),  pointer :: eqn_of_state => null() !< The equation of state type
  type(ocean_grid_type), pointer :: G => null()  !< The ocean grid type
  type(verticalGrid_type), pointer :: GV => null()  !< The model's vertical ocean grid
  type(unit_scale_type), pointer :: US => null() !< A dimensional unit scaling type

  !> The volume cell measure (special diagnostic) manager id
  integer :: volume_cell_measure_dm_id = -1

#if defined(DEBUG) || defined(__DO_SAFETY_CHECKS__)
  ! Keep a copy of h so that we know whether it has changed [H ~> m or kg m-2]. If it has then
  ! need the target grid for vertical remapping needs to have been updated.
  real, dimension(:,:,:), allocatable :: h_old
#endif

  !> Number of checksum-only diagnostics
  integer :: num_chksum_diags

  real, dimension(:,:,:), allocatable :: h_begin !< Layer thicknesses at the beginning of the timestep used
                                                 !! for remapping of extensive variables

end type diag_ctrl

!>@{ CPU clocks
integer :: id_clock_diag_mediator, id_clock_diag_remap, id_clock_diag_grid_updates
!>@}

contains

!> Sets up diagnostics axes
subroutine set_axes_info(G, GV, US, param_file, diag_cs, set_vertical)
  type(ocean_grid_type),   intent(inout) :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file structure
  type(diag_ctrl),         intent(inout) :: diag_cs !< Diagnostics control structure
  logical,       optional, intent(in)    :: set_vertical !< If true or missing, set up
                                                       !! vertical axes
  ! Local variables
  integer :: id_xq, id_yq, id_zl, id_zi, id_xh, id_yh, id_null
  integer :: id_zl_native, id_zi_native
  integer :: i, j, k, nz
  real :: zlev(GV%ke), zinter(GV%ke+1)
  logical :: set_vert
  real, allocatable, dimension(:) :: IaxB,iax
  real, allocatable, dimension(:) :: JaxB,jax


  set_vert = .true. ; if (present(set_vertical)) set_vert = set_vertical


  if (diag_cs%grid_space_axes) then
    allocate(IaxB(G%IsgB:G%IegB))
    do i=G%IsgB, G%IegB
      Iaxb(i)=real(i)
    enddo
    allocate(iax(G%isg:G%ieg))
    do i=G%isg, G%ieg
      iax(i)=real(i)-0.5
    enddo
    allocate(JaxB(G%JsgB:G%JegB))
    do j=G%JsgB, G%JegB
      JaxB(j)=real(j)
    enddo
    allocate(jax(G%jsg:G%jeg))
    do j=G%jsg, G%jeg
      jax(j)=real(j)-0.5
    enddo
  endif

  ! Horizontal axes for the native grids
  if (G%symmetric) then
    if (diag_cs%grid_space_axes) then
      id_xq = diag_axis_init('iq', IaxB(G%isgB:G%iegB), 'none', 'x', &
          'q point grid-space longitude', G%Domain, position=EAST)
      id_yq = diag_axis_init('jq', JaxB(G%jsgB:G%jegB), 'none', 'y', &
          'q point grid space latitude', G%Domain, position=NORTH)
    else
      id_xq = diag_axis_init('xq', G%gridLonB(G%isgB:G%iegB), G%x_axis_units, 'x', &
          'q point nominal longitude', G%Domain, position=EAST)
      id_yq = diag_axis_init('yq', G%gridLatB(G%jsgB:G%jegB), G%y_axis_units, 'y', &
          'q point nominal latitude', G%Domain, position=NORTH)
    endif
  else
    if (diag_cs%grid_space_axes) then
      id_xq = diag_axis_init('Iq', IaxB(G%isg:G%ieg), 'none', 'x', &
          'q point grid-space longitude', G%Domain, position=EAST)
      id_yq = diag_axis_init('Jq', JaxB(G%jsg:G%jeg), 'none', 'y', &
          'q point grid space latitude', G%Domain, position=NORTH)
    else
      id_xq = diag_axis_init('xq', G%gridLonB(G%isg:G%ieg), G%x_axis_units, 'x', &
          'q point nominal longitude', G%Domain, position=EAST)
      id_yq = diag_axis_init('yq', G%gridLatB(G%jsg:G%jeg), G%y_axis_units, 'y', &
          'q point nominal latitude', G%Domain, position=NORTH)
    endif
  endif

  if (diag_cs%grid_space_axes) then
    id_xh = diag_axis_init('ih', iax(G%isg:G%ieg), 'none', 'x', &
        'h point grid-space longitude', G%Domain, position=EAST)
    id_yh = diag_axis_init('jh', jax(G%jsg:G%jeg), 'none', 'y', &
        'h point grid space latitude', G%Domain, position=NORTH)
  else
    id_xh = diag_axis_init('xh', G%gridLonT(G%isg:G%ieg), G%x_axis_units, 'x', &
        'h point nominal longitude', G%Domain)
    id_yh = diag_axis_init('yh', G%gridLatT(G%jsg:G%jeg), G%y_axis_units, 'y', &
        'h point nominal latitude', G%Domain)
  endif

  if (set_vert) then
    nz = GV%ke
    zinter(1:nz+1) = GV%sInterface(1:nz+1)
    zlev(1:nz) = GV%sLayer(1:nz)
    id_zl = diag_axis_init('zl', zlev, trim(GV%zAxisUnits), 'z', &
                           'Layer '//trim(GV%zAxisLongName), direction=GV%direction)
    id_zi = diag_axis_init('zi', zinter, trim(GV%zAxisUnits), 'z', &
                           'Interface '//trim(GV%zAxisLongName), direction=GV%direction)
  else
    id_zl = -1 ; id_zi = -1
  endif
  id_zl_native = id_zl ; id_zi_native = id_zi
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

  ! Axis group for special null axis from diag manager.
  id_null = diag_axis_init('scalar_axis', (/0./), 'none', 'N', 'none', null_axis=.true.)
  call define_axes_group(diag_cs, (/ id_null /), diag_cs%axesNull)

  !Non-native Non-downsampled
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
    call diag_remap_configure_axes(diag_cs%diag_remap_cs(i), GV, US, param_file)

    ! Allocate these arrays since the size of the diagnostic array is now known
    allocate(diag_cs%diag_remap_cs(i)%h(G%isd:G%ied,G%jsd:G%jed, diag_cs%diag_remap_cs(i)%nz))
    allocate(diag_cs%diag_remap_cs(i)%h_extensive(G%isd:G%ied,G%jsd:G%jed, diag_cs%diag_remap_cs(i)%nz))

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

  if (diag_cs%grid_space_axes) then
    deallocate(IaxB, iax, JaxB, jax)
  endif
  !Define the downsampled axes
  call set_axes_info_dsamp(G, GV, param_file, diag_cs, id_zl_native, id_zi_native)

  call diag_grid_storage_init(diag_CS%diag_grid_temp, G, GV, diag_CS)

end subroutine set_axes_info

subroutine set_axes_info_dsamp(G, GV, param_file, diag_cs, id_zl_native, id_zi_native)
  type(ocean_grid_type), intent(in) :: G !< Ocean grid structure
  type(verticalGrid_type), intent(in)  :: GV !< ocean vertical grid structure
  type(param_file_type), intent(in)    :: param_file !< Parameter file structure
  type(diag_ctrl),       intent(inout) :: diag_cs !< Diagnostics control structure
  integer,               intent(in)    :: id_zl_native !< ID of native layers
  integer,               intent(in)    :: id_zi_native !< ID of native interfaces

  ! Local variables
  integer :: id_xq, id_yq, id_zl, id_zi, id_xh, id_yh
  integer :: i, j, k, nz, dl
  real, dimension(:), pointer :: gridLonT_dsamp =>NULL()
  real, dimension(:), pointer :: gridLatT_dsamp =>NULL()
  real, dimension(:), pointer :: gridLonB_dsamp =>NULL()
  real, dimension(:), pointer :: gridLatB_dsamp =>NULL()

  id_zl = id_zl_native ; id_zi = id_zi_native
  !Axes group for native downsampled diagnostics
  do dl=2,MAX_DSAMP_LEV
    if (dl /= 2) call MOM_error(FATAL, "set_axes_info_dsamp: Downsample level other than 2 is not supported yet!")
    if (G%symmetric) then
      allocate(gridLonB_dsamp(diag_cs%dsamp(dl)%isgB:diag_cs%dsamp(dl)%iegB))
      allocate(gridLatB_dsamp(diag_cs%dsamp(dl)%jsgB:diag_cs%dsamp(dl)%jegB))
      do i=diag_cs%dsamp(dl)%isgB,diag_cs%dsamp(dl)%iegB;  gridLonB_dsamp(i) = G%gridLonB(G%isgB+dl*i); enddo
      do j=diag_cs%dsamp(dl)%jsgB,diag_cs%dsamp(dl)%jegB;  gridLatB_dsamp(j) = G%gridLatB(G%jsgB+dl*j); enddo
      id_xq = diag_axis_init('xq', gridLonB_dsamp, G%x_axis_units, 'x', &
            'q point nominal longitude', G%Domain, coarsen=2)
      id_yq = diag_axis_init('yq', gridLatB_dsamp, G%y_axis_units, 'y', &
            'q point nominal latitude', G%Domain, coarsen=2)
      deallocate(gridLonB_dsamp,gridLatB_dsamp)
    else
      allocate(gridLonB_dsamp(diag_cs%dsamp(dl)%isg:diag_cs%dsamp(dl)%ieg))
      allocate(gridLatB_dsamp(diag_cs%dsamp(dl)%jsg:diag_cs%dsamp(dl)%jeg))
      do i=diag_cs%dsamp(dl)%isg,diag_cs%dsamp(dl)%ieg;  gridLonB_dsamp(i) = G%gridLonB(G%isg+dl*i-2); enddo
      do j=diag_cs%dsamp(dl)%jsg,diag_cs%dsamp(dl)%jeg;  gridLatB_dsamp(j) = G%gridLatB(G%jsg+dl*j-2); enddo
      id_xq = diag_axis_init('xq', gridLonB_dsamp, G%x_axis_units, 'x', &
            'q point nominal longitude', G%Domain, coarsen=2)
      id_yq = diag_axis_init('yq', gridLatB_dsamp, G%y_axis_units, 'y', &
            'q point nominal latitude', G%Domain, coarsen=2)
      deallocate(gridLonB_dsamp,gridLatB_dsamp)
    endif

    allocate(gridLonT_dsamp(diag_cs%dsamp(dl)%isg:diag_cs%dsamp(dl)%ieg))
    allocate(gridLatT_dsamp(diag_cs%dsamp(dl)%jsg:diag_cs%dsamp(dl)%jeg))
    do i=diag_cs%dsamp(dl)%isg,diag_cs%dsamp(dl)%ieg;  gridLonT_dsamp(i) = G%gridLonT(G%isg+dl*i-2); enddo
    do j=diag_cs%dsamp(dl)%jsg,diag_cs%dsamp(dl)%jeg;  gridLatT_dsamp(j) = G%gridLatT(G%jsg+dl*j-2); enddo
    id_xh = diag_axis_init('xh', gridLonT_dsamp, G%x_axis_units, 'x', &
          'h point nominal longitude', G%Domain, coarsen=2)
    id_yh = diag_axis_init('yh', gridLatT_dsamp, G%y_axis_units, 'y', &
          'h point nominal latitude', G%Domain, coarsen=2)

    deallocate(gridLonT_dsamp,gridLatT_dsamp)

    ! Axis groupings for the model layers
    call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yh, id_zL /), diag_cs%dsamp(dl)%axesTL, dl, &
            x_cell_method='mean', y_cell_method='mean', v_cell_method='mean', &
            is_h_point=.true., is_layer=.true., xyave_axes=diag_cs%axesZL)
    call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yq, id_zL /), diag_cs%dsamp(dl)%axesBL, dl, &
            x_cell_method='point', y_cell_method='point', v_cell_method='mean', &
            is_q_point=.true., is_layer=.true.)
    call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yh, id_zL /), diag_cs%dsamp(dl)%axesCuL, dl, &
            x_cell_method='point', y_cell_method='mean', v_cell_method='mean', &
            is_u_point=.true., is_layer=.true., xyave_axes=diag_cs%axesZL)
    call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yq, id_zL /), diag_cs%dsamp(dl)%axesCvL, dl, &
            x_cell_method='mean', y_cell_method='point', v_cell_method='mean', &
            is_v_point=.true., is_layer=.true., xyave_axes=diag_cs%axesZL)

    ! Axis groupings for the model interfaces
    call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yh, id_zi /), diag_cs%dsamp(dl)%axesTi, dl, &
            x_cell_method='mean', y_cell_method='mean', v_cell_method='point', &
            is_h_point=.true., is_interface=.true., xyave_axes=diag_cs%axesZi)
    call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yq, id_zi /), diag_cs%dsamp(dl)%axesBi, dl, &
            x_cell_method='point', y_cell_method='point', v_cell_method='point', &
            is_q_point=.true., is_interface=.true.)
    call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yh, id_zi /), diag_cs%dsamp(dl)%axesCui, dl, &
            x_cell_method='point', y_cell_method='mean', v_cell_method='point', &
            is_u_point=.true., is_interface=.true., xyave_axes=diag_cs%axesZi)
    call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yq, id_zi /), diag_cs%dsamp(dl)%axesCvi, dl, &
            x_cell_method='mean', y_cell_method='point', v_cell_method='point', &
            is_v_point=.true., is_interface=.true., xyave_axes=diag_cs%axesZi)

    ! Axis groupings for 2-D arrays
    call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yh /), diag_cs%dsamp(dl)%axesT1, dl, &
            x_cell_method='mean', y_cell_method='mean', is_h_point=.true.)
    call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yq /), diag_cs%dsamp(dl)%axesB1, dl, &
            x_cell_method='point', y_cell_method='point', is_q_point=.true.)
    call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yh /), diag_cs%dsamp(dl)%axesCu1, dl, &
            x_cell_method='point', y_cell_method='mean', is_u_point=.true.)
    call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yq /), diag_cs%dsamp(dl)%axesCv1, dl, &
            x_cell_method='mean', y_cell_method='point', is_v_point=.true.)

    !Non-native axes
    if (diag_cs%num_diag_coords>0) then
      allocate(diag_cs%dsamp(dl)%remap_axesTL(diag_cs%num_diag_coords))
      allocate(diag_cs%dsamp(dl)%remap_axesBL(diag_cs%num_diag_coords))
      allocate(diag_cs%dsamp(dl)%remap_axesCuL(diag_cs%num_diag_coords))
      allocate(diag_cs%dsamp(dl)%remap_axesCvL(diag_cs%num_diag_coords))
      allocate(diag_cs%dsamp(dl)%remap_axesTi(diag_cs%num_diag_coords))
      allocate(diag_cs%dsamp(dl)%remap_axesBi(diag_cs%num_diag_coords))
      allocate(diag_cs%dsamp(dl)%remap_axesCui(diag_cs%num_diag_coords))
      allocate(diag_cs%dsamp(dl)%remap_axesCvi(diag_cs%num_diag_coords))
    endif

    do i=1, diag_cs%num_diag_coords
      ! For each possible diagnostic coordinate
      !call diag_remap_configure_axes(diag_cs%diag_remap_cs(i), GV, param_file)

      ! This vertical coordinate has been configured so can be used.
      if (diag_remap_axes_configured(diag_cs%diag_remap_cs(i))) then

        ! This fetches the 1D-axis id for layers and interfaces and overwrite
        ! id_zl and id_zi from above. It also returns the number of layers.
        call diag_remap_get_axes_info(diag_cs%diag_remap_cs(i), nz, id_zL, id_zi)

        ! Axes for z layers
        call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yh, id_zL /), diag_cs%dsamp(dl)%remap_axesTL(i), dl, &
                nz=nz, vertical_coordinate_number=i, &
                x_cell_method='mean', y_cell_method='mean', v_cell_method='mean', &
                is_h_point=.true., is_layer=.true., is_native=.false., needs_remapping=.true., &
                xyave_axes=diag_cs%remap_axesZL(i))

        !! \note Remapping for B points is not yet implemented so needs_remapping is not
        !! provided for remap_axesBL
        call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yq, id_zL /), diag_cs%dsamp(dl)%remap_axesBL(i), dl, &
                nz=nz, vertical_coordinate_number=i, &
                x_cell_method='point', y_cell_method='point', v_cell_method='mean', &
                is_q_point=.true., is_layer=.true., is_native=.false.)

        call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yh, id_zL /), diag_cs%dsamp(dl)%remap_axesCuL(i), dl, &
                nz=nz, vertical_coordinate_number=i, &
                x_cell_method='point', y_cell_method='mean', v_cell_method='mean', &
                is_u_point=.true., is_layer=.true., is_native=.false., needs_remapping=.true., &
                xyave_axes=diag_cs%remap_axesZL(i))

        call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yq, id_zL /), diag_cs%dsamp(dl)%remap_axesCvL(i), dl, &
                nz=nz, vertical_coordinate_number=i, &
                x_cell_method='mean', y_cell_method='point', v_cell_method='mean', &
                is_v_point=.true., is_layer=.true., is_native=.false., needs_remapping=.true., &
                xyave_axes=diag_cs%remap_axesZL(i))

        ! Axes for z interfaces
        call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yh, id_zi /), diag_cs%dsamp(dl)%remap_axesTi(i), dl, &
                nz=nz, vertical_coordinate_number=i, &
                x_cell_method='mean', y_cell_method='mean', v_cell_method='point', &
                is_h_point=.true., is_interface=.true., is_native=.false., needs_interpolating=.true., &
                xyave_axes=diag_cs%remap_axesZi(i))

        !! \note Remapping for B points is not yet implemented so needs_remapping is not provided for remap_axesBi
        call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yq, id_zi /), diag_cs%dsamp(dl)%remap_axesBi(i), dl, &
                nz=nz, vertical_coordinate_number=i, &
                x_cell_method='point', y_cell_method='point', v_cell_method='point', &
                is_q_point=.true., is_interface=.true., is_native=.false.)

        call define_axes_group_dsamp(diag_cs, (/ id_xq, id_yh, id_zi /), diag_cs%dsamp(dl)%remap_axesCui(i), dl, &
                nz=nz, vertical_coordinate_number=i, &
                x_cell_method='point', y_cell_method='mean', v_cell_method='point', &
                is_u_point=.true., is_interface=.true., is_native=.false., &
                needs_interpolating=.true., xyave_axes=diag_cs%remap_axesZi(i))

        call define_axes_group_dsamp(diag_cs, (/ id_xh, id_yq, id_zi /), diag_cs%dsamp(dl)%remap_axesCvi(i), dl, &
                nz=nz, vertical_coordinate_number=i, &
                x_cell_method='mean', y_cell_method='point', v_cell_method='point', &
                is_v_point=.true., is_interface=.true., is_native=.false., &
                needs_interpolating=.true., xyave_axes=diag_cs%remap_axesZi(i))
      endif
    enddo
  enddo

end subroutine set_axes_info_dsamp


!> set_masks_for_axes sets up the 2d and 3d masks for diagnostics using the current grid
!! recorded after calling diag_update_remap_grids()
subroutine set_masks_for_axes(G, diag_cs)
  type(ocean_grid_type), target, intent(in) :: G !< The ocean grid type.
  type(diag_ctrl),               pointer    :: diag_cs !< A pointer to a type with many variables
                                                       !! used for diagnostics
  ! Local variables
  integer :: c, nk, i, j, k, ii, jj
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
      do J=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
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

  !Allocate and initialize the downsampled masks for the axes
  call set_masks_for_axes_dsamp(G, diag_cs)

end subroutine set_masks_for_axes

subroutine set_masks_for_axes_dsamp(G, diag_cs)
  type(ocean_grid_type), target, intent(in) :: G !< The ocean grid type.
  type(diag_ctrl),               pointer    :: diag_cs !< A pointer to a type with many variables
                                                       !! used for diagnostics
  ! Local variables
  integer :: c, nk, i, j, k, ii, jj
  integer :: dl
  type(axes_grp), pointer :: axes => NULL(), h_axes => NULL() ! Current axes, for convenience

  !Each downsampled axis needs both downsampled and non-downsampled mask
  !The downsampled mask is needed for sending out the diagnostics output via diag_manager
  !The non-downsampled mask is needed for downsampling the diagnostics field
  do dl=2,MAX_DSAMP_LEV
    if (dl /= 2) call MOM_error(FATAL, "set_masks_for_axes_dsamp: Downsample level other than 2 is not supported!")
    do c=1, diag_cs%num_diag_coords
      ! Level/layer h-points in diagnostic coordinate
      axes => diag_cs%remap_axesTL(c)
      call downsample_mask(axes%mask3d, diag_cs%dsamp(dl)%remap_axesTL(c)%dsamp(dl)%mask3d, dl,G%isc, G%jsc,  &
              G%HId2%isc, G%HId2%iec, G%HId2%jsc, G%HId2%jec, G%HId2%isd, G%HId2%ied, G%HId2%jsd, G%HId2%jed)
      diag_cs%dsamp(dl)%remap_axesTL(c)%mask3d => axes%mask3d !set non-downsampled mask
      ! Level/layer u-points in diagnostic coordinate
      axes => diag_cs%remap_axesCuL(c)
      call downsample_mask(axes%mask3d, diag_cs%dsamp(dl)%remap_axesCuL(c)%dsamp(dl)%mask3d, dl,G%IscB,G%JscB, &
               G%HId2%IscB,G%HId2%IecB,G%HId2%jsc, G%HId2%jec,G%HId2%IsdB,G%HId2%IedB,G%HId2%jsd, G%HId2%jed)
      diag_cs%dsamp(dl)%remap_axesCul(c)%mask3d => axes%mask3d !set non-downsampled mask
      ! Level/layer v-points in diagnostic coordinate
      axes => diag_cs%remap_axesCvL(c)
      call downsample_mask(axes%mask3d, diag_cs%dsamp(dl)%remap_axesCvL(c)%dsamp(dl)%mask3d, dl,G%isc ,G%JscB, &
              G%HId2%isc ,G%HId2%iec, G%HId2%JscB,G%HId2%JecB,G%HId2%isd ,G%HId2%ied, G%HId2%JsdB,G%HId2%JedB)
      diag_cs%dsamp(dl)%remap_axesCvL(c)%mask3d => axes%mask3d !set non-downsampled mask
      ! Level/layer q-points in diagnostic coordinate
      axes => diag_cs%remap_axesBL(c)
      call downsample_mask(axes%mask3d, diag_cs%dsamp(dl)%remap_axesBL(c)%dsamp(dl)%mask3d, dl,G%IscB,G%JscB, &
              G%HId2%IscB,G%HId2%IecB,G%HId2%JscB,G%HId2%JecB,G%HId2%IsdB,G%HId2%IedB,G%HId2%JsdB,G%HId2%JedB)
      diag_cs%dsamp(dl)%remap_axesBL(c)%mask3d => axes%mask3d !set non-downsampled mask
      ! Interface h-points in diagnostic coordinate (w-point)
      axes => diag_cs%remap_axesTi(c)
      call downsample_mask(axes%mask3d, diag_cs%dsamp(dl)%remap_axesTi(c)%dsamp(dl)%mask3d, dl,G%isc, G%jsc,  &
              G%HId2%isc, G%HId2%iec, G%HId2%jsc, G%HId2%jec, G%HId2%isd, G%HId2%ied, G%HId2%jsd, G%HId2%jed)
      diag_cs%dsamp(dl)%remap_axesTi(c)%mask3d => axes%mask3d !set non-downsampled mask
      ! Interface u-points in diagnostic coordinate
      axes => diag_cs%remap_axesCui(c)
      call downsample_mask(axes%mask3d, diag_cs%dsamp(dl)%remap_axesCui(c)%dsamp(dl)%mask3d, dl,G%IscB,G%JscB, &
              G%HId2%IscB,G%HId2%IecB,G%HId2%jsc, G%HId2%jec,G%HId2%IsdB,G%HId2%IedB,G%HId2%jsd, G%HId2%jed)
      diag_cs%dsamp(dl)%remap_axesCui(c)%mask3d => axes%mask3d !set non-downsampled mask
      ! Interface v-points in diagnostic coordinate
      axes => diag_cs%remap_axesCvi(c)
      call downsample_mask(axes%mask3d, diag_cs%dsamp(dl)%remap_axesCvi(c)%dsamp(dl)%mask3d, dl,G%isc ,G%JscB, &
              G%HId2%isc ,G%HId2%iec, G%HId2%JscB,G%HId2%JecB,G%HId2%isd ,G%HId2%ied, G%HId2%JsdB,G%HId2%JedB)
      diag_cs%dsamp(dl)%remap_axesCvi(c)%mask3d => axes%mask3d !set non-downsampled mask
      ! Interface q-points in diagnostic coordinate
      axes => diag_cs%remap_axesBi(c)
      call downsample_mask(axes%mask3d, diag_cs%dsamp(dl)%remap_axesBi(c)%dsamp(dl)%mask3d, dl,G%IscB,G%JscB, &
              G%HId2%IscB,G%HId2%IecB,G%HId2%JscB,G%HId2%JecB,G%HId2%IsdB,G%HId2%IedB,G%HId2%JsdB,G%HId2%JedB)
      diag_cs%dsamp(dl)%remap_axesBi(c)%mask3d => axes%mask3d !set non-downsampled mask
    enddo
  enddo
end subroutine set_masks_for_axes_dsamp

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

!> Defines a group of downsampled "axes" from list of handles
subroutine define_axes_group_dsamp(diag_cs, handles, axes, dl, nz, vertical_coordinate_number, &
                             x_cell_method, y_cell_method, v_cell_method, &
                             is_h_point, is_q_point, is_u_point, is_v_point, &
                             is_layer, is_interface, &
                             is_native, needs_remapping, needs_interpolating, &
                             xyave_axes)
  type(diag_ctrl), target,    intent(in)  :: diag_cs !< Diagnostics control structure
  integer, dimension(:),      intent(in)  :: handles !< A list of 1D axis handles
  type(axes_grp),             intent(out) :: axes    !< The group of 1D axes
  integer,                    intent(in)  :: dl      !< Downsample level
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
  axes%downsample_level = dl
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

  axes%dsamp(dl)%mask2d => null()
  if (axes%rank==2) then
    if (axes%is_h_point) axes%dsamp(dl)%mask2d => diag_cs%dsamp(dl)%mask2dT
    if (axes%is_u_point) axes%dsamp(dl)%mask2d => diag_cs%dsamp(dl)%mask2dCu
    if (axes%is_v_point) axes%dsamp(dl)%mask2d => diag_cs%dsamp(dl)%mask2dCv
    if (axes%is_q_point) axes%dsamp(dl)%mask2d => diag_cs%dsamp(dl)%mask2dBu
  endif
  ! A static 3d mask for non-native coordinates can only be setup when a grid is available
  axes%dsamp(dl)%mask3d => null()
  if (axes%rank==3 .and. axes%is_native) then
    ! Native variables can/should use the native masks copied into diag_cs
    if (axes%is_layer) then
      if (axes%is_h_point) axes%dsamp(dl)%mask3d => diag_cs%dsamp(dl)%mask3dTL
      if (axes%is_u_point) axes%dsamp(dl)%mask3d => diag_cs%dsamp(dl)%mask3dCuL
      if (axes%is_v_point) axes%dsamp(dl)%mask3d => diag_cs%dsamp(dl)%mask3dCvL
      if (axes%is_q_point) axes%dsamp(dl)%mask3d => diag_cs%dsamp(dl)%mask3dBL
    elseif (axes%is_interface) then
      if (axes%is_h_point) axes%dsamp(dl)%mask3d => diag_cs%dsamp(dl)%mask3dTi
      if (axes%is_u_point) axes%dsamp(dl)%mask3d => diag_cs%dsamp(dl)%mask3dCui
      if (axes%is_v_point) axes%dsamp(dl)%mask3d => diag_cs%dsamp(dl)%mask3dCvi
      if (axes%is_q_point) axes%dsamp(dl)%mask3d => diag_cs%dsamp(dl)%mask3dBi
    endif
  endif

end subroutine define_axes_group_dsamp

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
  real :: locfield
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
    locfield = field
    if (diag%conversion_factor /= 0.) &
      locfield = locfield * diag%conversion_factor

    if (diag_cs%diag_as_chksum) then
      call chksum0(locfield, diag%debug_str, logunit=diag_cs%chksum_iounit)
    elseif (is_stat) then
      used = send_data_infra(diag%fms_diag_id, locfield)
    elseif (diag_cs%ave_enabled) then
      used = send_data_infra(diag%fms_diag_id, locfield, diag_cs%time_end)
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

    if (diag_cs%diag_as_chksum) then
      call zchksum(locfield, diag%debug_str, logunit=diag_cs%chksum_iounit)
    elseif (is_stat) then
      used = send_data_infra(diag%fms_diag_id, locfield)
    elseif (diag_cs%ave_enabled) then
      used = send_data_infra(diag%fms_diag_id, locfield, time=diag_cs%time_end, weight=diag_cs%time_int)
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
  real,    optional,target, intent(in) :: mask(:,:) !< If present, use this real array as the data mask.

  ! Local variables
  real, dimension(:,:), pointer :: locfield
  real, dimension(:,:), pointer :: locmask
  character(len=300) :: mesg
  logical :: used, is_stat
  integer :: cszi, cszj, dszi, dszj
  integer :: isv, iev, jsv, jev, i, j, chksum, isv_o,jsv_o
  real, dimension(:,:), allocatable, target :: locfield_dsamp
  real, dimension(:,:), allocatable, target :: locmask_dsamp
  integer :: dl

  locfield => NULL()
  locmask => NULL()
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

  if (present(mask)) then
    locmask => mask
  elseif (.NOT. is_stat) then
    if (associated(diag%axes%mask2d)) locmask => diag%axes%mask2d
  endif

  dl=1
  if (.NOT. is_stat) dl = diag%axes%downsample_level !static field downsample i not supported yet
  !Downsample the diag field and mask (if present)
  if (dl > 1) then
    isv_o = isv ; jsv_o = jsv
    call downsample_diag_field(locfield, locfield_dsamp, dl, diag_cs, diag,isv,iev,jsv,jev, mask)
    if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) deallocate( locfield )
    locfield => locfield_dsamp
    if (present(mask)) then
      call downsample_field_2d(locmask, locmask_dsamp, dl, MSK, locmask, diag_cs,diag,isv_o,jsv_o,isv,iev,jsv,jev)
      locmask => locmask_dsamp
    elseif (associated(diag%axes%dsamp(dl)%mask2d)) then
      locmask => diag%axes%dsamp(dl)%mask2d
    endif
  endif

  if (diag_cs%diag_as_chksum) then
    if (diag%axes%is_h_point) then
      call hchksum(locfield, diag%debug_str, diag_cs%G%HI, &
                   logunit=diag_cs%chksum_iounit)
    elseif (diag%axes%is_u_point) then
      call uchksum(locfield, diag%debug_str, diag_cs%G%HI, &
                   logunit=diag_cs%chksum_iounit)
    elseif (diag%axes%is_v_point) then
      call vchksum(locfield, diag%debug_str, diag_cs%G%HI, &
                   logunit=diag_cs%chksum_iounit)
    elseif (diag%axes%is_q_point) then
      call Bchksum(locfield, diag%debug_str, diag_cs%G%HI, &
                   logunit=diag_cs%chksum_iounit)
    else
      call MOM_error(FATAL, "post_data_2d_low: unknown axis type.")
    endif
  else
    if (is_stat) then
      if (present(mask)) then
        call assert(size(locfield) == size(locmask), &
            'post_data_2d_low is_stat: mask size mismatch: '//diag%debug_str)
        used = send_data_infra(diag%fms_diag_id, locfield, &
                         is_in=isv, ie_in=iev, js_in=jsv, je_in=jev, rmask=locmask)
     !elseif (associated(diag%axes%mask2d)) then
     !  used = send_data(diag%fms_diag_id, locfield, &
     !                   is_in=isv, ie_in=iev, js_in=jsv, je_in=jev, rmask=diag%axes%mask2d)
      else
        used = send_data_infra(diag%fms_diag_id, locfield, &
                         is_in=isv, ie_in=iev, js_in=jsv, je_in=jev)
      endif
    elseif (diag_cs%ave_enabled) then
      if (associated(locmask)) then
        call assert(size(locfield) == size(locmask), &
            'post_data_2d_low: mask size mismatch: '//diag%debug_str)
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
  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.) .and. dl<2) &
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
                                                !! remapping this diagnostic [H ~> m or kg m-2].

  ! Local variables
  type(diag_type), pointer :: diag => null()
  integer :: nz, i, j, k
  real, dimension(:,:,:), allocatable :: remapped_field
  logical :: staggered_in_x, staggered_in_y
  real, dimension(:,:,:), pointer :: h_diag => NULL()

  if (id_clock_diag_mediator>0) call cpu_clock_begin(id_clock_diag_mediator)

  ! For intensive variables only, we can choose to use a different diagnostic grid
  ! to map to
  if (present(alt_h)) then
    h_diag => alt_h
  else
    h_diag => diag_cs%h
  endif

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
      call vertically_reintegrate_diag_field(                                    &
        diag_cs%diag_remap_cs(diag%axes%vertical_coordinate_number), diag_cs%G,  &
        diag_cs%h_begin,                                                         &
        diag_cs%diag_remap_cs(diag%axes%vertical_coordinate_number)%h_extensive, &
        staggered_in_x, staggered_in_y, diag%axes%mask3d, field, remapped_field)
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
      call diag_remap_do_remap(diag_cs%diag_remap_cs(diag%axes%vertical_coordinate_number), &
              diag_cs%G, diag_cs%GV, h_diag, staggered_in_x, staggered_in_y, &
              diag%axes%mask3d, field, remapped_field)
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
              diag%axes%mask3d, field, remapped_field)
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
  real,    optional,target, intent(in) :: mask(:,:,:) !< If present, use this real array as the data mask.

  ! Local variables
  real, dimension(:,:,:), pointer :: locfield
  real, dimension(:,:,:), pointer :: locmask
  character(len=300) :: mesg
  logical :: used  ! The return value of send_data is not used for anything.
  logical :: staggered_in_x, staggered_in_y
  logical :: is_stat
  integer :: cszi, cszj, dszi, dszj
  integer :: isv, iev, jsv, jev, ks, ke, i, j, k, isv_c, jsv_c, isv_o,jsv_o
  integer :: chksum
  real, dimension(:,:,:), allocatable, target :: locfield_dsamp
  real, dimension(:,:,:), allocatable, target :: locmask_dsamp
  integer :: dl

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

  ks = lbound(field,3) ; ke = ubound(field,3)
  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) then
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

    do k=ks,ke ; do j=jsv,jev ; do i=isv,iev
      if (field(i,j,k) == diag_cs%missing_value) then
        locfield(i,j,k) = diag_cs%missing_value
      else
        locfield(i,j,k) = field(i,j,k) * diag%conversion_factor
      endif
    enddo ; enddo ; enddo
  else
    locfield => field
  endif

  if (present(mask)) then
    locmask => mask
  elseif (associated(diag%axes%mask3d)) then
    locmask => diag%axes%mask3d
  endif

  dl=1
  if (.NOT. is_stat) dl = diag%axes%downsample_level !static field downsample i not supported yet
  !Downsample the diag field and mask (if present)
  if (dl > 1) then
    isv_o = isv ; jsv_o = jsv
    call downsample_diag_field(locfield, locfield_dsamp, dl, diag_cs, diag,isv,iev,jsv,jev, mask)
    if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) deallocate( locfield )
    locfield => locfield_dsamp
    if (present(mask)) then
      call downsample_field_3d(locmask, locmask_dsamp, dl, MSK, locmask, diag_cs,diag,isv_o,jsv_o,isv,iev,jsv,jev)
      locmask => locmask_dsamp
    elseif (associated(diag%axes%dsamp(dl)%mask3d)) then
      locmask => diag%axes%dsamp(dl)%mask3d
    endif
  endif

  if (diag%fms_diag_id>0) then
    if (diag_cs%diag_as_chksum) then
      if (diag%axes%is_h_point) then
        call hchksum(locfield, diag%debug_str, diag_cs%G%HI, &
                     logunit=diag_cs%chksum_iounit)
      elseif (diag%axes%is_u_point) then
        call uchksum(locfield, diag%debug_str, diag_cs%G%HI, &
                     logunit=diag_cs%chksum_iounit)
      elseif (diag%axes%is_v_point) then
        call vchksum(locfield, diag%debug_str, diag_cs%G%HI, &
                     logunit=diag_cs%chksum_iounit)
      elseif (diag%axes%is_q_point) then
        call Bchksum(locfield, diag%debug_str, diag_cs%G%HI, &
                     logunit=diag_cs%chksum_iounit)
      else
        call MOM_error(FATAL, "post_data_3d_low: unknown axis type.")
      endif
    else
      if (is_stat) then
        if (present(mask)) then
          call assert(size(locfield) == size(locmask), &
              'post_data_3d_low is_stat: mask size mismatch: '//diag%debug_str)
          used = send_data_infra(diag%fms_diag_id, locfield, &
                         is_in=isv, ie_in=iev, js_in=jsv, je_in=jev, rmask=locmask)
       !elseif (associated(diag%axes%mask2d)) then
       !  used = send_data(diag%fms_diag_id, locfield, &
       !                   is_in=isv, ie_in=iev, js_in=jsv, je_in=jev, rmask=diag%axes%mask2d)
        else
          used = send_data_infra(diag%fms_diag_id, locfield, &
                           is_in=isv, ie_in=iev, js_in=jsv, je_in=jev)
        endif
      elseif (diag_cs%ave_enabled) then
        if (associated(locmask)) then
          call assert(size(locfield) == size(locmask), &
              'post_data_3d_low: mask size mismatch: '//diag%debug_str)
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
  endif

  if (diag%fms_xyave_diag_id>0) then
    call post_xy_average(diag_cs, diag, locfield)
  endif

  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.) .and. dl<2) &
    deallocate( locfield )

end subroutine post_data_3d_low

!> Post the horizontally area-averaged diagnostic
subroutine post_xy_average(diag_cs, diag, field)
  type(diag_type),   intent(in) :: diag !< This diagnostic
  real,    target,   intent(in) :: field(:,:,:) !< Diagnostic field
  type(diag_ctrl),   intent(in) :: diag_cs !< Diagnostics mediator control structure
  ! Local variable
  real, dimension(size(field,3)) :: averaged_field
  logical, dimension(size(field,3)) :: averaged_mask
  logical :: staggered_in_x, staggered_in_y, used
  integer :: nz, remap_nz, coord

  if (.not. diag_cs%ave_enabled) then
    return
  endif

  staggered_in_x = diag%axes%is_u_point .or. diag%axes%is_q_point
  staggered_in_y = diag%axes%is_v_point .or. diag%axes%is_q_point

  if (diag%axes%is_native) then
    call horizontally_average_diag_field(diag_cs%G, diag_cs%GV, diag_cs%h, &
                                         staggered_in_x, staggered_in_y, &
                                         diag%axes%is_layer, diag%v_extensive, &
                                         field, &
                                         averaged_field, averaged_mask)
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

    call horizontally_average_diag_field(diag_cs%G, diag_cs%GV, &
                                         diag_cs%diag_remap_cs(coord)%h, &
                                         staggered_in_x, staggered_in_y, &
                                         diag%axes%is_layer, diag%v_extensive, &
                                         field, averaged_field, averaged_mask)
  endif

  if (diag_cs%diag_as_chksum) then
    call zchksum(averaged_field, trim(diag%debug_str)//'_xyave', &
                 logunit=diag_CS%chksum_iounit)
  else
    used = send_data_infra(diag%fms_xyave_diag_id, averaged_field, &
                           time=diag_cs%time_end, weight=diag_cs%time_int, mask=averaged_mask)
  endif
end subroutine post_xy_average

!> This subroutine enables the accumulation of time averages over the specified time interval.
subroutine enable_averaging(time_int_in, time_end_in, diag_cs)
  real,            intent(in)    :: time_int_in !< The time interval [s] over which any
                                                !!  values that are offered are valid.
  type(time_type), intent(in)    :: time_end_in !< The end time of the valid interval
  type(diag_ctrl), intent(inout) :: diag_CS !< Structure used to regulate diagnostic output

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
  real,  optional, intent(in)    :: T_to_s   !< A conversion factor for time_int to [s].
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
  type(diag_ctrl), intent(inout) :: diag_CS !< Structure used to regulate diagnostic output

  diag_cs%time_int = 0.0
  diag_cs%ave_enabled = .false.

end subroutine disable_averaging

!> Call this subroutine to determine whether the averaging is
!! currently enabled.  .true. is returned if it is.
function query_averaging_enabled(diag_cs, time_int, time_end)
  type(diag_ctrl),           intent(in)  :: diag_CS  !< Structure used to regulate diagnostic output
  real,            optional, intent(out) :: time_int !< Current setting of diag%time_int [s]
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
integer function register_diag_field(module_name, field_name, axes_in, init_time, &
            long_name, units, missing_value, range, mask_variant, standard_name,      &
            verbose, do_not_log, err_msg, interp_method, tile_count, cmor_field_name, &
            cmor_long_name, cmor_units, cmor_standard_name, cell_methods, &
            x_cell_method, y_cell_method, v_cell_method, conversion, v_extensive)
  character(len=*),           intent(in) :: module_name !< Name of this module, usually "ocean_model"
                                                        !! or "ice_shelf_model"
  character(len=*),           intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp),     target, intent(in) :: axes_in   !< Container w/ up to 3 integer handles that
                                                      !! indicates axes for this field
  type(time_type),            intent(in) :: init_time !< Time at which a field is first available?
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
  type(axes_grp), pointer :: axes => null()
  integer :: dm_id, i, dl
  character(len=256) :: msg, cm_string
  character(len=256) :: new_module_name
  character(len=480) :: module_list, var_list
  integer :: num_modnm, num_varnm
  logical :: active

  axes => axes_in
  MOM_missing_value = axes%diag_cs%missing_value
  if (present(missing_value)) MOM_missing_value = missing_value

  diag_cs => axes%diag_cs
  dm_id = -1

  if (axes_in%id == diag_cs%axesTL%id) then
    axes => diag_cs%axesTL
  elseif (axes_in%id == diag_cs%axesBL%id) then
    axes => diag_cs%axesBL
  elseif (axes_in%id == diag_cs%axesCuL%id ) then
    axes => diag_cs%axesCuL
  elseif (axes_in%id == diag_cs%axesCvL%id) then
    axes => diag_cs%axesCvL
  elseif (axes_in%id == diag_cs%axesTi%id) then
    axes => diag_cs%axesTi
  elseif (axes_in%id == diag_cs%axesBi%id) then
    axes => diag_cs%axesBi
  elseif (axes_in%id == diag_cs%axesCui%id ) then
    axes => diag_cs%axesCui
  elseif (axes_in%id == diag_cs%axesCvi%id) then
    axes => diag_cs%axesCvi
  endif

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
             cell_methods=cell_methods, x_cell_method=x_cell_method, &
             y_cell_method=y_cell_method, v_cell_method=v_cell_method, &
             conversion=conversion, v_extensive=v_extensive)
  if (associated(axes%xyave_axes)) then
    num_varnm = 2 ; var_list = "{"//trim(field_name)//","//trim(field_name)//"_xyave"
  else
    num_varnm = 1 ; var_list = "{"//trim(field_name)
  endif
  if (present(cmor_field_name)) then
    if (associated(axes%xyave_axes)) then
      num_varnm = num_varnm + 2
      var_list = trim(var_list)//","//trim(cmor_field_name)//","//trim(cmor_field_name)//"_xyave"
    else
      num_varnm = num_varnm + 1
      var_list = trim(var_list)//","//trim(cmor_field_name)
    endif
  endif
  var_list = trim(var_list)//"}"

  ! For each diagnostic coordinate register the diagnostic again under a different module name
  do i=1,diag_cs%num_diag_coords
    new_module_name = trim(module_name)//'_'//trim(diag_cs%diag_remap_cs(i)%diag_module_suffix)

    ! Register diagnostics remapped to z vertical coordinate
    if (axes_in%rank == 3) then
      remap_axes => null()
      if ((axes_in%id == diag_cs%axesTL%id)) then
        remap_axes => diag_cs%remap_axesTL(i)
      elseif (axes_in%id == diag_cs%axesBL%id) then
        remap_axes => diag_cs%remap_axesBL(i)
      elseif (axes_in%id == diag_cs%axesCuL%id ) then
        remap_axes => diag_cs%remap_axesCuL(i)
      elseif (axes_in%id == diag_cs%axesCvL%id) then
        remap_axes => diag_cs%remap_axesCvL(i)
      elseif (axes_in%id == diag_cs%axesTi%id) then
        remap_axes => diag_cs%remap_axesTi(i)
      elseif (axes_in%id == diag_cs%axesBi%id) then
        remap_axes => diag_cs%remap_axesBi(i)
      elseif (axes_in%id == diag_cs%axesCui%id ) then
        remap_axes => diag_cs%remap_axesCui(i)
      elseif (axes_in%id == diag_cs%axesCvi%id) then
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
          module_list = trim(module_list)//","//trim(new_module_name)
          num_modnm = num_modnm + 1
        endif ! remap_axes%needs_remapping
      endif ! associated(remap_axes)
    endif ! axes%rank == 3
  enddo ! i

  !Register downsampled diagnostics
  do dl=2,MAX_DSAMP_LEV
    ! Do not attempt to checksum the downsampled diagnostics
    if (diag_cs%diag_as_chksum) cycle

    new_module_name = trim(module_name)//'_d2'

    if (axes_in%rank == 3 .or. axes_in%rank == 2 ) then
      axes => null()
      if (axes_in%id == diag_cs%axesTL%id) then
        axes => diag_cs%dsamp(dl)%axesTL
      elseif (axes_in%id == diag_cs%axesBL%id) then
        axes => diag_cs%dsamp(dl)%axesBL
      elseif (axes_in%id == diag_cs%axesCuL%id ) then
        axes => diag_cs%dsamp(dl)%axesCuL
      elseif (axes_in%id == diag_cs%axesCvL%id) then
        axes => diag_cs%dsamp(dl)%axesCvL
      elseif (axes_in%id == diag_cs%axesTi%id) then
        axes => diag_cs%dsamp(dl)%axesTi
      elseif (axes_in%id == diag_cs%axesBi%id) then
        axes => diag_cs%dsamp(dl)%axesBi
      elseif (axes_in%id == diag_cs%axesCui%id ) then
        axes => diag_cs%dsamp(dl)%axesCui
      elseif (axes_in%id == diag_cs%axesCvi%id) then
        axes => diag_cs%dsamp(dl)%axesCvi
      elseif (axes_in%id == diag_cs%axesT1%id) then
        axes => diag_cs%dsamp(dl)%axesT1
      elseif (axes_in%id == diag_cs%axesB1%id) then
        axes => diag_cs%dsamp(dl)%axesB1
      elseif (axes_in%id == diag_cs%axesCu1%id ) then
        axes => diag_cs%dsamp(dl)%axesCu1
      elseif (axes_in%id == diag_cs%axesCv1%id) then
        axes => diag_cs%dsamp(dl)%axesCv1
      else
        !Niki: Should we worry about these, e.g., diag_to_Z_CS?
        call MOM_error(WARNING,"register_diag_field: Could not find a proper axes for " &
              //trim(new_module_name)//"-"//trim(field_name))
      endif
    endif
    ! Register the native diagnostic
    if (associated(axes)) then
       active = register_diag_field_expand_cmor(dm_id, new_module_name, field_name, axes, &
                init_time, long_name=long_name, units=units, missing_value=MOM_missing_value, &
                range=range, mask_variant=mask_variant, standard_name=standard_name, &
                verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                interp_method=interp_method, tile_count=tile_count, &
                cmor_field_name=cmor_field_name, cmor_long_name=cmor_long_name, &
                cmor_units=cmor_units, cmor_standard_name=cmor_standard_name, &
                cell_methods=cell_methods, x_cell_method=x_cell_method, &
                y_cell_method=y_cell_method, v_cell_method=v_cell_method, &
                conversion=conversion, v_extensive=v_extensive)
      module_list = trim(module_list)//","//trim(new_module_name)
      num_modnm = num_modnm + 1
    endif

    ! For each diagnostic coordinate register the diagnostic again under a different module name
    do i=1,diag_cs%num_diag_coords
      new_module_name = trim(module_name)//'_'//trim(diag_cs%diag_remap_cs(i)%diag_module_suffix)//'_d2'

      ! Register diagnostics remapped to z vertical coordinate
      if (axes_in%rank == 3) then
        remap_axes => null()
        if ((axes_in%id == diag_cs%axesTL%id)) then
          remap_axes => diag_cs%dsamp(dl)%remap_axesTL(i)
        elseif (axes_in%id == diag_cs%axesBL%id) then
          remap_axes => diag_cs%dsamp(dl)%remap_axesBL(i)
        elseif (axes_in%id == diag_cs%axesCuL%id ) then
          remap_axes => diag_cs%dsamp(dl)%remap_axesCuL(i)
        elseif (axes_in%id == diag_cs%axesCvL%id) then
          remap_axes => diag_cs%dsamp(dl)%remap_axesCvL(i)
        elseif (axes_in%id == diag_cs%axesTi%id) then
          remap_axes => diag_cs%dsamp(dl)%remap_axesTi(i)
        elseif (axes_in%id == diag_cs%axesBi%id) then
          remap_axes => diag_cs%dsamp(dl)%remap_axesBi(i)
        elseif (axes_in%id == diag_cs%axesCui%id ) then
          remap_axes => diag_cs%dsamp(dl)%remap_axesCui(i)
        elseif (axes_in%id == diag_cs%axesCvi%id) then
          remap_axes => diag_cs%dsamp(dl)%remap_axesCvi(i)
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
            module_list = trim(module_list)//","//trim(new_module_name)
            num_modnm = num_modnm + 1
          endif ! remap_axes%needs_remapping
        endif ! associated(remap_axes)
      endif ! axes%rank == 3
    enddo ! i
  enddo

  if (is_root_pe() .and. (diag_CS%available_diag_doc_unit > 0)) then
    msg = ''
    if (present(cmor_field_name)) msg = 'CMOR equivalent is "'//trim(cmor_field_name)//'"'
    call attach_cell_methods(-1, axes, cm_string, cell_methods, &
                             x_cell_method, y_cell_method, v_cell_method, &
                             v_extensive=v_extensive)
    module_list = trim(module_list)//"}"
    if (num_modnm <= 1) module_list = module_name
    if (num_varnm <= 1) var_list = ""

    call log_available_diag(dm_id>0, module_list, field_name, cm_string, msg, diag_CS, &
                            long_name, units, standard_name, variants=var_list)
  endif

  register_diag_field = dm_id

end function register_diag_field

!> Returns True if either the native or CMOr version of the diagnostic were registered. Updates 'dm_id'
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
  if (.not. diag_cs%diag_as_chksum) &
    call attach_cell_methods(fms_id, axes, cm_string, cell_methods, &
                             x_cell_method, y_cell_method, v_cell_method, &
                             v_extensive=v_extensive)
  ! Associated horizontally area-averaged diagnostic
  fms_xyave_id = DIAG_FIELD_NOT_FOUND
  if (associated(axes%xyave_axes)) then
    fms_xyave_id = register_diag_field_expand_axes(module_name, trim(field_name)//'_xyave', &
             axes%xyave_axes, init_time, &
             long_name=long_name, units=units, missing_value=MOM_missing_value, &
             range=range, mask_variant=mask_variant, standard_name=standard_name, &
             verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
             interp_method=interp_method, tile_count=tile_count)
    if (.not. diag_cs%diag_as_chksum) &
      call attach_cell_methods(fms_xyave_id, axes%xyave_axes, cm_string, &
                               cell_methods, v_cell_method, v_extensive=v_extensive)
  endif
  this_diag => null()
  if (fms_id /= DIAG_FIELD_NOT_FOUND .or. fms_xyave_id /= DIAG_FIELD_NOT_FOUND) then
    call add_diag_to_list(diag_cs, dm_id, fms_id, this_diag, axes, module_name, field_name, msg)
    this_diag%fms_xyave_diag_id = fms_xyave_id
    !Encode and save the cell methods for this diag
    call add_xyz_method(this_diag, axes, x_cell_method, y_cell_method, v_cell_method, v_extensive)
    if (present(v_extensive)) this_diag%v_extensive = v_extensive
    if (present(conversion)) this_diag%conversion_factor = conversion
    register_diag_field_expand_cmor = .true.
  endif

  ! For the CMOR variation of the above diagnostic
  if (present(cmor_field_name) .and. .not. diag_cs%diag_as_chksum) then
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
    endif
    this_diag => null()
    if (fms_id /= DIAG_FIELD_NOT_FOUND .or. fms_xyave_id /= DIAG_FIELD_NOT_FOUND) then
      call add_diag_to_list(diag_cs, dm_id, fms_id, this_diag, axes, module_name, field_name, msg)
      this_diag%fms_xyave_diag_id = fms_xyave_id
      !Encode and save the cell methods for this diag
      call add_xyz_method(this_diag, axes, x_cell_method, y_cell_method, v_cell_method, v_extensive)
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
  if (axes%diag_cs%diag_as_chksum) then
    fms_id = axes%diag_cs%num_chksum_diags + 1
    axes%diag_cs%num_chksum_diags = fms_id
  elseif (present(interp_method) .or. axes%is_h_point) then
    ! If interp_method is provided we must use it
    if (area_id>0) then
      if (volume_id>0) then
        fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method=interp_method, tile_count=tile_count, area=area_id, volume=volume_id)
      else
        fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method=interp_method, tile_count=tile_count, area=area_id)
      endif
    else
      if (volume_id>0) then
        fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method=interp_method, tile_count=tile_count, volume=volume_id)
      else
        fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
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
        fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method='none', tile_count=tile_count, area=area_id, volume=volume_id)
      else
        fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method='none', tile_count=tile_count, area=area_id)
      endif
    else
      if (volume_id>0) then
        fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
                 init_time, long_name=long_name, units=units, missing_value=missing_value, &
                 range=range, mask_variant=mask_variant, standard_name=standard_name, &
                 verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
                 interp_method='none', tile_count=tile_count, volume=volume_id)
      else
        fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
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

!> Adds the encoded "cell_methods" for a diagnostics as a diag% property
!! This allows access to the cell_method for a given diagnostics at the time of sending
subroutine add_xyz_method(diag, axes, x_cell_method, y_cell_method, v_cell_method, v_extensive)
  type(diag_type),          pointer       :: diag !< This diagnostic
  type(axes_grp),             intent(in)  :: axes !< Container w/ up to 3 integer handles that indicates
                                                  !! axes for this field
  character(len=*), optional, intent(in)  :: x_cell_method !< Specifies the cell method for the x-direction.
                                                         !! Use '' have no method.
  character(len=*), optional, intent(in)  :: y_cell_method !< Specifies the cell method for the y-direction.
                                                         !! Use '' have no method.
  character(len=*), optional, intent(in)  :: v_cell_method !< Specifies the cell method for the vertical direction.
                                                         !! Use '' have no method.
  logical,          optional, intent(in)  :: v_extensive !< True for vertically extensive fields
                                                         !! (vertically integrated). Default/absent for intensive.
  integer :: xyz_method
  character(len=9) :: mstr

  !This is a simple way to encode the cell method information made from 3 strings
  !(x_cell_method,y_cell_method,v_cell_method) in a 3 digit integer xyz
  !x_cell_method,y_cell_method,v_cell_method can each be 'point' or 'sum' or 'mean'
  !We can encode these with setting  1 for 'point', 2 for 'sum, 3 for 'mean' in
  !the 100s position for x, 10s position for y, 1s position for z
  !E.g., x:sum,y:point,z:mean is 213

  xyz_method = 111

  mstr = diag%axes%v_cell_method
  if (present(v_extensive)) then
    if (present(v_cell_method)) call MOM_error(FATAL, "attach_cell_methods: " // &
       'Vertical cell method was specified along with the vertically extensive flag.')
    if (v_extensive) then
      mstr='sum'
    else
      mstr='mean'
    endif
  elseif (present(v_cell_method)) then
    mstr = v_cell_method
  endif
  if (trim(mstr)=='sum') then
    xyz_method = xyz_method + 1
  elseif (trim(mstr)=='mean') then
    xyz_method = xyz_method + 2
  endif

  mstr = diag%axes%y_cell_method
  if (present(y_cell_method)) mstr = y_cell_method
  if (trim(mstr)=='sum') then
    xyz_method = xyz_method + 10
  elseif (trim(mstr)=='mean') then
    xyz_method = xyz_method + 20
  endif

  mstr = diag%axes%x_cell_method
  if (present(x_cell_method)) mstr = x_cell_method
  if (trim(mstr)=='sum') then
    xyz_method = xyz_method + 100
  elseif (trim(mstr)=='mean') then
    xyz_method = xyz_method + 200
  endif

  diag%xyz_method = xyz_method
end subroutine add_xyz_method

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
    if (present(v_cell_method)) then
      if (present(v_extensive)) call MOM_error(FATAL, "attach_cell_methods: " // &
           'Vertical cell method was specified along with the vertically extensive flag.')
      if (len(trim(v_cell_method))>0) then
        if (axes%rank==1) then
          call get_MOM_diag_axis_name(axes%handles(1), axis_name)
        elseif (axes%rank==3) then
          call get_MOM_diag_axis_name(axes%handles(3), axis_name)
        endif
        call MOM_diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(v_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(v_cell_method)
      endif
    elseif (present(v_extensive)) then
      if (v_extensive) then
        if (axes%rank==1) then
          call get_MOM_diag_axis_name(axes%handles(1), axis_name)
        elseif (axes%rank==3) then
          call get_MOM_diag_axis_name(axes%handles(3), axis_name)
        endif
        call MOM_diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':sum')
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':sum'
      endif
    else
      if (len(trim(axes%v_cell_method))>0) then
        if (axes%rank==1) then
          call get_MOM_diag_axis_name(axes%handles(1), axis_name)
        elseif (axes%rank==3) then
          call get_MOM_diag_axis_name(axes%handles(3), axis_name)
        endif
        call MOM_diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(axes%v_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(axes%v_cell_method)
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
    endif
  endif

  ! Document diagnostics in list of available diagnostics
  if (is_root_pe() .and. diag_CS%available_diag_doc_unit > 0) then
    if (present(cmor_field_name)) then
      call log_available_diag(associated(diag), module_name, field_name, '', '', diag_CS, &
                              long_name, units, standard_name, &
                              variants="{"//trim(field_name)//","//trim(cmor_field_name)//"}")
    else
      call log_available_diag(associated(diag), module_name, field_name, '', '', diag_CS, &
                              long_name, units, standard_name)
    endif
  endif

  register_scalar_field = dm_id

end function register_scalar_field

!> Registers a static diagnostic, returning an integer handle
function register_static_field(module_name, field_name, axes, &
            long_name, units, missing_value, range, mask_variant, standard_name, &
            do_not_log, interp_method, tile_count, &
            cmor_field_name, cmor_long_name, cmor_units, cmor_standard_name, area, &
            x_cell_method, y_cell_method, area_cell_method, conversion)
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
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to file

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

  ! Document diagnostics in list of available diagnostics
  if (is_root_pe() .and. diag_CS%available_diag_doc_unit > 0) then
    if (present(cmor_field_name)) then
      call log_available_diag(associated(diag), module_name, field_name, '', '', diag_CS, &
                              long_name, units, standard_name, &
                              variants="{"//trim(field_name)//","//trim(cmor_field_name)//"}")
    else
      call log_available_diag(associated(diag), module_name, field_name, '', '', diag_CS, &
                              long_name, units, standard_name)
    endif
  endif

  register_static_field = dm_id

end function register_static_field

!> Describe an option setting in the diagnostic files.
subroutine describe_option(opt_name, value, diag_CS)
  character(len=*), intent(in) :: opt_name !< The name of the option
  character(len=*), intent(in) :: value   !< A character string with the setting of the option.
  type(diag_ctrl),  intent(in) :: diag_CS !< Structure used to regulate diagnostic output

  character(len=480) :: mesg
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
        case ("q")  ; axes => diag_cs%axesBL
        case ("h")  ; axes => diag_cs%axesTL
        case ("u")  ; axes => diag_cs%axesCuL
        case ("v")  ; axes => diag_cs%axesCvL
        case ("Bu") ; axes => diag_cs%axesBL
        case ("T")  ; axes => diag_cs%axesTL
        case ("Cu") ; axes => diag_cs%axesCuL
        case ("Cv") ; axes => diag_cs%axesCvL
        case ("z")  ; axes => diag_cs%axeszL
        case default ; call MOM_error(FATAL, "ocean_register_diag: " // &
                                      "unknown hor_grid component "//trim(hor_grid))
      end select

    case ("i")
      select case (hor_grid)
        case ("q")  ; axes => diag_cs%axesBi
        case ("h")  ; axes => diag_cs%axesTi
        case ("u")  ; axes => diag_cs%axesCui
        case ("v")  ; axes => diag_cs%axesCvi
        case ("Bu") ; axes => diag_cs%axesBi
        case ("T")  ; axes => diag_cs%axesTi
        case ("Cu") ; axes => diag_cs%axesCui
        case ("Cv") ; axes => diag_cs%axesCvi
        case ("z")  ; axes => diag_cs%axeszi
        case default ; call MOM_error(FATAL, "ocean_register_diag: " // &
                                      "unknown hor_grid component "//trim(hor_grid))
      end select

    case ("1")
      select case (hor_grid)
        case ("q")  ; axes => diag_cs%axesB1
        case ("h")  ; axes => diag_cs%axesT1
        case ("u")  ; axes => diag_cs%axesCu1
        case ("v")  ; axes => diag_cs%axesCv1
        case ("Bu") ; axes => diag_cs%axesB1
        case ("T")  ; axes => diag_cs%axesT1
        case ("Cu") ; axes => diag_cs%axesCu1
        case ("Cv") ; axes => diag_cs%axesCv1
        case default ; call MOM_error(FATAL, "ocean_register_diag: " // &
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

  call MOM_diag_manager_init(err_msg=err_msg)
end subroutine diag_mediator_infrastructure_init

!> diag_mediator_init initializes the MOM diag_mediator and opens the available
!! diagnostics file, if appropriate.
subroutine diag_mediator_init(G, GV, US, nz, param_file, diag_cs, doc_file_dir)
  type(ocean_grid_type), target, intent(inout) :: G  !< The ocean grid type.
  type(verticalGrid_type), target, intent(in)  :: GV !< The ocean vertical grid structure
  type(unit_scale_type),   target, intent(in)  :: US !< A dimensional unit scaling type
  integer,                    intent(in)    :: nz    !< The number of layers in the model's native grid.
  type(param_file_type),      intent(in)    :: param_file !< Parameter file structure
  type(diag_ctrl),            intent(inout) :: diag_cs !< A pointer to a type with many variables
                                                     !! used for diagnostics
  character(len=*), optional, intent(in)    :: doc_file_dir !< A directory in which to create the
                                                     !! file

  ! This subroutine initializes the diag_mediator and the diag_manager.
  ! The grid type should have its dimensions set by this point, but it
  ! is not necessary that the metrics and axis labels be set up yet.

  ! Local variables
  integer :: ios, i, new_unit
  logical :: opened, new_file
  logical :: answers_2018, default_2018_answers
  character(len=8)   :: this_pe
  character(len=240) :: doc_file, doc_file_dflt, doc_path
  character(len=240), allocatable :: diag_coords(:)
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40) :: mdl = "MOM_diag_mediator" ! This module's name.
  character(len=32) :: filename_appendix = '' !fms appendix to filename for ensemble runs

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
                 'The number of diagnostic vertical coordinates to use. '//&
                 'For each coordinate, an entry in DIAG_COORDS must be provided.', &
                 default=1)
  call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false.)
  call get_param(param_file, mdl, "REMAPPING_2018_ANSWERS", answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)
  call get_param(param_file, mdl, 'USE_GRID_SPACE_DIAGNOSTIC_AXES', diag_cs%grid_space_axes, &
                 'If true, use a grid index coordinate convention for diagnostic axes. ',&
                 default=.false.)

  if (diag_cs%num_diag_coords>0) then
    allocate(diag_coords(diag_cs%num_diag_coords))
    if (diag_cs%num_diag_coords==1) then ! The default is to provide just one instance of Z*
      call get_param(param_file, mdl, 'DIAG_COORDS', diag_coords, &
                 'A list of string tuples associating diag_table modules to '//&
                 'a coordinate definition used for diagnostics. Each string '//&
                 'is of the form "MODULE_SUFFIX PARAMETER_SUFFIX COORDINATE_NAME".', &
                 default='z Z ZSTAR')
    else ! If using more than 1 diagnostic coordinate, all must be explicitly defined
      call get_param(param_file, mdl, 'DIAG_COORDS', diag_coords, &
                 'A list of string tuples associating diag_table modules to '//&
                 'a coordinate definition used for diagnostics. Each string '//&
                 'is of the form "MODULE_SUFFIX,PARAMETER_SUFFIX,COORDINATE_NAME".', &
                 fail_if_missing=.true.)
    endif
    allocate(diag_cs%diag_remap_cs(diag_cs%num_diag_coords))
    ! Initialize each diagnostic vertical coordinate
    do i=1, diag_cs%num_diag_coords
      call diag_remap_init(diag_cs%diag_remap_cs(i), diag_coords(i), answers_2018=answers_2018)
    enddo
    deallocate(diag_coords)
  endif

  call get_param(param_file, mdl, 'DIAG_MISVAL', diag_cs%missing_value, &
                 'Set the default missing value to use for diagnostics.', &
                 default=1.e20)
  call get_param(param_file, mdl, 'DIAG_AS_CHKSUM', diag_cs%diag_as_chksum, &
                 'Instead of writing diagnostics to the diag manager, write '//&
                 'a text file containing the checksum (bitcount) of the array.',  &
                 default=.false.)

  if (diag_cs%diag_as_chksum) &
    diag_cs%num_chksum_diags = 0

  ! Keep pointers grid, h, T, S needed diagnostic remapping
  diag_cs%G => G
  diag_cs%GV => GV
  diag_cs%US => US
  diag_cs%h => null()
  diag_cs%T => null()
  diag_cs%S => null()
  diag_cs%eqn_of_state => null()

  allocate(diag_cs%h_begin(G%isd:G%ied,G%jsd:G%jed,nz))
#if defined(DEBUG) || defined(__DO_SAFETY_CHECKS__)
  allocate(diag_cs%h_old(G%isd:G%ied,G%jsd:G%jed,nz))
  diag_cs%h_old(:,:,:) = 0.0
#endif

  diag_cs%is = G%isc - (G%isd-1) ; diag_cs%ie = G%iec - (G%isd-1)
  diag_cs%js = G%jsc - (G%jsd-1) ; diag_cs%je = G%jec - (G%jsd-1)
  diag_cs%isd = G%isd ; diag_cs%ied = G%ied
  diag_cs%jsd = G%jsd ; diag_cs%jed = G%jed

  !Downsample indices for dl=2 (should be generalized to arbitrary dl, perhaps via a G array)
  diag_cs%dsamp(2)%isc = G%HId2%isc - (G%HId2%isd-1) ; diag_cs%dsamp(2)%iec = G%HId2%iec - (G%HId2%isd-1)
  diag_cs%dsamp(2)%jsc = G%HId2%jsc - (G%HId2%jsd-1) ; diag_cs%dsamp(2)%jec = G%HId2%jec - (G%HId2%jsd-1)
  diag_cs%dsamp(2)%isd = G%HId2%isd ; diag_cs%dsamp(2)%ied = G%HId2%ied
  diag_cs%dsamp(2)%jsd = G%HId2%jsd ; diag_cs%dsamp(2)%jed = G%HId2%jed
  diag_cs%dsamp(2)%isg = G%HId2%isg ; diag_cs%dsamp(2)%ieg = G%HId2%ieg
  diag_cs%dsamp(2)%jsg = G%HId2%jsg ; diag_cs%dsamp(2)%jeg = G%HId2%jeg
  diag_cs%dsamp(2)%isgB = G%HId2%isgB ; diag_cs%dsamp(2)%iegB = G%HId2%iegB
  diag_cs%dsamp(2)%jsgB = G%HId2%jsgB ; diag_cs%dsamp(2)%jegB = G%HId2%jegB

  ! Initialze available diagnostic log file
  if (is_root_pe() .and. (diag_CS%available_diag_doc_unit < 0)) then
    write(this_pe,'(i6.6)') PE_here()
    doc_file_dflt = "available_diags."//this_pe
    call get_param(param_file, mdl, "AVAILABLE_DIAGS_FILE", doc_file, &
                 "A file into which to write a list of all available "//&
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

end subroutine diag_mediator_init

!> Set pointers to the default state fields used to remap diagnostics.
subroutine diag_set_state_ptrs(h, T, S, eqn_of_state, diag_cs)
  real, dimension(:,:,:), target, intent(in   ) :: h !< the model thickness array [H ~> m or kg m-2]
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
subroutine diag_update_remap_grids(diag_cs, alt_h, alt_T, alt_S, update_intensive, update_extensive )
  type(diag_ctrl),        intent(inout) :: diag_cs      !< Diagnostics control structure
  real, target, optional, intent(in   ) :: alt_h(:,:,:) !< Used if remapped grids should be something other than
                                                        !! the current thicknesses [H ~> m or kg m-2]
  real, target, optional, intent(in   ) :: alt_T(:,:,:) !< Used if remapped grids should be something other than
                                                        !! the current temperatures
  real, target, optional, intent(in   ) :: alt_S(:,:,:) !< Used if remapped grids should be something other than
                                                        !! the current salinity
  logical, optional,      intent(in   ) :: update_intensive !< If true (default), update the grids used for
                                                            !! intensive diagnostics
  logical, optional,      intent(in   ) :: update_extensive !< If true (not default), update the grids used for
                                                            !! intensive diagnostics
  ! Local variables
  integer :: i
  real, dimension(:,:,:), pointer :: h_diag => NULL() ! The layer thickneses for diagnostics [H ~> m or kg m-2]
  real, dimension(:,:,:), pointer :: T_diag => NULL(), S_diag => NULL()
  logical :: update_intensive_local, update_extensive_local

  ! Set values based on optional input arguments
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

  ! Defaults here are based on wanting to update intensive quantities frequently as soon as the model state changes.
  ! Conversely, for extensive quantities, in an effort to close budgets and to be consistent with the total time
  ! tendency, we construct the diagnostic grid at the beginning of the baroclinic timestep and remap all extensive
  ! quantities to the same grid
  update_intensive_local = .true.
  if (present(update_intensive)) update_intensive_local = update_intensive
  update_extensive_local = .false.
  if (present(update_extensive)) update_extensive_local = update_extensive

  if (id_clock_diag_grid_updates>0) call cpu_clock_begin(id_clock_diag_grid_updates)

  if (diag_cs%diag_grid_overridden) then
    call MOM_error(FATAL, "diag_update_remap_grids was called, but current grids in "// &
                          "diagnostic structure have been overridden")
  endif

  if (update_intensive_local) then
    do i=1, diag_cs%num_diag_coords
      call diag_remap_update(diag_cs%diag_remap_cs(i), diag_cs%G, diag_cs%GV, diag_cs%US, h_diag, T_diag, S_diag, &
                             diag_cs%eqn_of_state, diag_cs%diag_remap_cs(i)%h)
    enddo
  endif
  if (update_extensive_local) then
    diag_cs%h_begin(:,:,:) = diag_cs%h(:,:,:)
    do i=1, diag_cs%num_diag_coords
      call diag_remap_update(diag_cs%diag_remap_cs(i), diag_cs%G, diag_cs%GV, diag_cs%US, h_diag, T_diag, S_diag, &
                             diag_cs%eqn_of_state, diag_cs%diag_remap_cs(i)%h_extensive)
    enddo
  endif

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

  !Allocate and initialize the downsampled masks
  call downsample_diag_masks_set(G, nz, diag_cs)

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
  if (diag_CS%chksum_iounit > -1) then
    close(diag_CS%chksum_iounit) ; diag_CS%chksum_iounit = -3
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
  do i=2,MAX_DSAMP_LEV
    deallocate(diag_cs%dsamp(i)%mask2dT)
    deallocate(diag_cs%dsamp(i)%mask2dBu)
    deallocate(diag_cs%dsamp(i)%mask2dCu)
    deallocate(diag_cs%dsamp(i)%mask2dCv)
    deallocate(diag_cs%dsamp(i)%mask3dTL)
    deallocate(diag_cs%dsamp(i)%mask3dBL)
    deallocate(diag_cs%dsamp(i)%mask3dCuL)
    deallocate(diag_cs%dsamp(i)%mask3dCvL)
    deallocate(diag_cs%dsamp(i)%mask3dTi)
    deallocate(diag_cs%dsamp(i)%mask3dBi)
    deallocate(diag_cs%dsamp(i)%mask3dCui)
    deallocate(diag_cs%dsamp(i)%mask3dCvi)
  enddo

#if defined(DEBUG) || defined(__DO_SAFETY_CHECKS__)
  deallocate(diag_cs%h_old)
#endif

  if (present(end_diag_manager)) then
    if (end_diag_manager) call MOM_diag_manager_end(time)
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
                              diag_CS, long_name, units, standard_name, variants)
  logical,          intent(in) :: used !< Whether this diagnostic was in the diag_table or not
  character(len=*), intent(in) :: module_name !< Name of the diagnostic module
  character(len=*), intent(in) :: field_name !< Name of this diagnostic field
  character(len=*), intent(in) :: cell_methods_string !< The spatial component of the CF cell_methods attribute
  character(len=*), intent(in) :: comment !< A comment to append after [Used|Unused]
  type(diag_ctrl),  intent(in) :: diag_CS  !< The diagnotics control structure
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
    write(diag_CS%available_diag_doc_unit, '(a,x,"(",a,")")') trim(mesg),trim(comment)
  else
    write(diag_CS%available_diag_doc_unit, '(a)') trim(mesg)
  endif
  call describe_option("modules", module_name, diag_CS)
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

  write(docunit, '(a,x,i9.8)') description, chksum
  flush(docunit)

end subroutine log_chksum_diag

!> Allocates fields necessary to store diagnostic remapping fields
subroutine diag_grid_storage_init(grid_storage, G, GV, diag)
  type(diag_grid_storage), intent(inout) :: grid_storage !< Structure containing a snapshot of the target grids
  type(ocean_grid_type),   intent(in)    :: G           !< Horizontal grid
  type(verticalGrid_type), intent(in)    :: GV          !< ocean vertical grid structure
  type(diag_ctrl),         intent(in)    :: diag        !< Diagnostic control structure used as the contructor
                                                        !! template for this routine

  integer :: m, nz
  grid_storage%num_diag_coords = diag%num_diag_coords

  ! Don't do anything else if there are no remapped coordinates
  if (grid_storage%num_diag_coords < 1) return

  ! Allocate memory for the native space
  allocate( grid_storage%h_state(G%isd:G%ied, G%jsd:G%jed, GV%ke))
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
  real, dimension(:,:,:),  intent(in)    :: h_state     !< Current model thicknesses [H ~> m or kg m-2]
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

!< Allocate and initialize the masks for downsampled diagostics in diag_cs
!! The downsampled masks in the axes would later "point" to these.
subroutine downsample_diag_masks_set(G, nz, diag_cs)
  type(ocean_grid_type), target, intent(in) :: G  !< The ocean grid type.
  integer,                       intent(in) :: nz !< The number of layers in the model's native grid.
  type(diag_ctrl),               pointer    :: diag_cs !< A pointer to a type with many variables
                                                       !! used for diagnostics
  ! Local variables
  integer :: i,j,k,ii,jj,dl

!print*,'original c extents ',G%isc,G%iec,G%jsc,G%jec
!print*,'original c extents ',G%iscb,G%iecb,G%jscb,G%jecb
!print*,'coarse   c extents ',G%HId2%isc,G%HId2%iec,G%HId2%jsc,G%HId2%jec
!print*,'original d extents ',G%isd,G%ied,G%jsd,G%jed
!print*,'coarse   d extents ',G%HId2%isd,G%HId2%ied,G%HId2%jsd,G%HId2%jed
! original c  extents           5          52           5          52
! original cB-nonsym extents    5          52           5          52
! original cB-sym    extents    4          52           4          52
! coarse   c extents            3          26           3          26
! original d extents            1          56           1          56
! original dB-nonsym extents    1          56           1          56
! original dB-sym extents       0          56           0          56
! coarse   d extents            1          28           1          28

  do dl=2,MAX_DSAMP_LEV
    ! 2d mask
    call downsample_mask(G%mask2dT, diag_cs%dsamp(dl)%mask2dT,  dl,G%isc, G%jsc,  &
            G%HId2%isc, G%HId2%iec, G%HId2%jsc, G%HId2%jec, G%HId2%isd, G%HId2%ied, G%HId2%jsd, G%HId2%jed)
    call downsample_mask(G%mask2dBu,diag_cs%dsamp(dl)%mask2dBu, dl,G%IscB,G%JscB, &
            G%HId2%IscB,G%HId2%IecB,G%HId2%JscB,G%HId2%JecB,G%HId2%IsdB,G%HId2%IedB,G%HId2%JsdB,G%HId2%JedB)
    call downsample_mask(G%mask2dCu,diag_cs%dsamp(dl)%mask2dCu, dl,G%IscB,G%JscB, &
            G%HId2%IscB,G%HId2%IecB,G%HId2%jsc, G%HId2%jec,G%HId2%IsdB,G%HId2%IedB,G%HId2%jsd, G%HId2%jed)
    call downsample_mask(G%mask2dCv,diag_cs%dsamp(dl)%mask2dCv, dl,G%isc ,G%JscB, &
            G%HId2%isc ,G%HId2%iec, G%HId2%JscB,G%HId2%JecB,G%HId2%isd ,G%HId2%ied, G%HId2%JsdB,G%HId2%JedB)
    ! 3d native masks are needed by diag_manager but the native variables
    ! can only be masked 2d - for ocean points, all layers exists.
    allocate(diag_cs%dsamp(dl)%mask3dTL(G%HId2%isd:G%HId2%ied,G%HId2%jsd:G%HId2%jed,1:nz))
    allocate(diag_cs%dsamp(dl)%mask3dBL(G%HId2%IsdB:G%HId2%IedB,G%HId2%JsdB:G%HId2%JedB,1:nz))
    allocate(diag_cs%dsamp(dl)%mask3dCuL(G%HId2%IsdB:G%HId2%IedB,G%HId2%jsd:G%HId2%jed,1:nz))
    allocate(diag_cs%dsamp(dl)%mask3dCvL(G%HId2%isd:G%HId2%ied,G%HId2%JsdB:G%HId2%JedB,1:nz))
    do k=1,nz
      diag_cs%dsamp(dl)%mask3dTL(:,:,k) = diag_cs%dsamp(dl)%mask2dT(:,:)
      diag_cs%dsamp(dl)%mask3dBL(:,:,k) = diag_cs%dsamp(dl)%mask2dBu(:,:)
      diag_cs%dsamp(dl)%mask3dCuL(:,:,k) = diag_cs%dsamp(dl)%mask2dCu(:,:)
      diag_cs%dsamp(dl)%mask3dCvL(:,:,k) = diag_cs%dsamp(dl)%mask2dCv(:,:)
    enddo
    allocate(diag_cs%dsamp(dl)%mask3dTi(G%HId2%isd:G%HId2%ied,G%HId2%jsd:G%HId2%jed,1:nz+1))
    allocate(diag_cs%dsamp(dl)%mask3dBi(G%HId2%IsdB:G%HId2%IedB,G%HId2%JsdB:G%HId2%JedB,1:nz+1))
    allocate(diag_cs%dsamp(dl)%mask3dCui(G%HId2%IsdB:G%HId2%IedB,G%HId2%jsd:G%HId2%jed,1:nz+1))
    allocate(diag_cs%dsamp(dl)%mask3dCvi(G%HId2%isd:G%HId2%ied,G%HId2%JsdB:G%HId2%JedB,1:nz+1))
    do k=1,nz+1
      diag_cs%dsamp(dl)%mask3dTi(:,:,k) = diag_cs%dsamp(dl)%mask2dT(:,:)
      diag_cs%dsamp(dl)%mask3dBi(:,:,k) = diag_cs%dsamp(dl)%mask2dBu(:,:)
      diag_cs%dsamp(dl)%mask3dCui(:,:,k) = diag_cs%dsamp(dl)%mask2dCu(:,:)
      diag_cs%dsamp(dl)%mask3dCvi(:,:,k) = diag_cs%dsamp(dl)%mask2dCv(:,:)
    enddo
  enddo
end subroutine downsample_diag_masks_set

!> Get the diagnostics-compute indices (to be passed to send_data) based on the shape of
!! the diag field (the same way they are deduced for non-downsampled fields)
subroutine downsample_diag_indices_get(fo1, fo2, dl, diag_cs, isv, iev, jsv, jev)
  integer,           intent(in)  :: fo1     !< The size of the diag field in x
  integer,           intent(in)  :: fo2     !< The size of the diag field in y
  integer,           intent(in)  :: dl      !< Integer downsample level
  type(diag_ctrl),   intent(in)  :: diag_CS !< Structure used to regulate diagnostic output
  integer,           intent(out) :: isv     !< i-start index for diagnostics
  integer,           intent(out) :: iev     !< i-end index for diagnostics
  integer,           intent(out) :: jsv     !< j-start index for diagnostics
  integer,           intent(out) :: jev     !< j-end index for diagnostics
  ! Local variables
  integer :: dszi,cszi,dszj,cszj,f1,f2
  character(len=500) :: mesg
  logical, save :: first_check = .true.

  !Check ONCE that the downsampled diag-compute domain is commensurate with the original
  !non-downsampled diag-compute domain.
  !This is a major limitation of the current implementation of the downsampled diagnostics.
  !We assume that the compute domain can be subdivided to dl*dl cells, hence avoiding the need of halo updates.
  !We want this check to error out only if there was a downsampled diagnostics requested and about to post that is
  !why the check is here and not in the init routines. This check need to be done only once, hence the outer if.
  if (first_check) then
    if (mod(diag_cs%ie-diag_cs%is+1, dl) /= 0 .OR. mod(diag_cs%je-diag_cs%js+1, dl) /= 0) then
      write (mesg,*) "Non-commensurate downsampled domain is not supported. "//&
             "Please choose a layout such that NIGLOBAL/Layout_X and NJGLOBAL/Layout_Y are both divisible by dl=",dl,&
             " Current domain extents: ", diag_cs%is,diag_cs%ie, diag_cs%js,diag_cs%je
      call MOM_error(FATAL,"downsample_diag_indices_get: "//trim(mesg))
    endif
    first_check = .false.
  endif

  cszi = diag_cs%dsamp(dl)%iec-diag_cs%dsamp(dl)%isc +1 ; dszi = diag_cs%dsamp(dl)%ied-diag_cs%dsamp(dl)%isd +1
  cszj = diag_cs%dsamp(dl)%jec-diag_cs%dsamp(dl)%jsc +1 ; dszj = diag_cs%dsamp(dl)%jed-diag_cs%dsamp(dl)%jsd +1
  isv = diag_cs%dsamp(dl)%isc ; iev = diag_cs%dsamp(dl)%iec
  jsv = diag_cs%dsamp(dl)%jsc ; jev = diag_cs%dsamp(dl)%jec
  f1 = fo1/dl
  f2 = fo2/dl
  !Correction for the symmetric case
  if (diag_cs%G%symmetric) then
    f1 = f1 + mod(fo1,dl)
    f2 = f2 + mod(fo2,dl)
  endif
  if ( f1 == dszi ) then
    isv = diag_cs%dsamp(dl)%isc ; iev = diag_cs%dsamp(dl)%iec   ! field on Data domain, take compute domain indcies
  !The rest is not taken with the full MOM6 diag_table
  elseif ( f1 == dszi + 1 ) then
    isv = diag_cs%dsamp(dl)%isc ; iev = diag_cs%dsamp(dl)%iec+1   ! Symmetric data domain
  elseif ( f1 == cszi) then
    isv = 1 ; iev = (diag_cs%dsamp(dl)%iec-diag_cs%dsamp(dl)%isc) +1  ! Computational domain
  elseif ( f1 == cszi + 1 ) then
    isv = 1 ; iev = (diag_cs%dsamp(dl)%iec-diag_cs%dsamp(dl)%isc) +2  ! Symmetric computational domain
  else
    write (mesg,*) " peculiar size ",f1," in i-direction\n"//&
          "does not match one of ", cszi, cszi+1, dszi, dszi+1
    call MOM_error(FATAL,"downsample_diag_indices_get: "//trim(mesg))
  endif
  if ( f2 == dszj ) then
    jsv = diag_cs%dsamp(dl)%jsc ; jev = diag_cs%dsamp(dl)%jec     ! Data domain
  elseif ( f2 == dszj + 1 ) then
    jsv = diag_cs%dsamp(dl)%jsc ; jev = diag_cs%dsamp(dl)%jec+1   ! Symmetric data domain
  elseif ( f2 == cszj) then
    jsv = 1 ; jev = (diag_cs%dsamp(dl)%jec-diag_cs%dsamp(dl)%jsc) +1  ! Computational domain
  elseif ( f2 == cszj + 1 ) then
    jsv = 1 ; jev = (diag_cs%dsamp(dl)%jec-diag_cs%dsamp(dl)%jsc) +2  ! Symmetric computational domain
  else
    write (mesg,*) " peculiar size ",f2," in j-direction\n"//&
          "does not match one of ", cszj, cszj+1, dszj, dszj+1
    call MOM_error(FATAL,"downsample_diag_indices_get: "//trim(mesg))
  endif
end subroutine downsample_diag_indices_get

!> This subroutine allocates and computes a downsampled array from an input array
!! It also determines the diagnostics-compurte indices for the downsampled array
!! 3d interface
subroutine downsample_diag_field_3d(locfield, locfield_dsamp, dl, diag_cs, diag, isv, iev, jsv, jev, mask)
  real, dimension(:,:,:), pointer :: locfield  !< Input array pointer
  real, dimension(:,:,:), allocatable, intent(inout) :: locfield_dsamp !< Output (downsampled) array
  type(diag_ctrl),   intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  type(diag_type),   intent(in) :: diag    !< A structure describing the diagnostic to post
  integer, intent(in) :: dl                !< Level of down sampling
  integer, intent(inout) :: isv            !< i-start index for diagnostics
  integer, intent(inout) :: iev            !< i-end index for diagnostics
  integer, intent(inout) :: jsv            !< j-start index for diagnostics
  integer, intent(inout) :: jev            !< j-end index for diagnostics
  real,    optional,target, intent(in) :: mask(:,:,:) !< If present, use this real array as the data mask.
  ! Locals
  real, dimension(:,:,:), pointer :: locmask
  integer :: f1,f2,isv_o,jsv_o

  locmask => NULL()
  !Get the correct indices corresponding to input field
  !Shape of the input diag field
  f1 = size(locfield, 1)
  f2 = size(locfield, 2)
  !Save the extents of the original (fine) domain
  isv_o = isv ; jsv_o = jsv
  !Get the shape of the downsampled field and overwrite isv,iev,jsv,jev with them
  call downsample_diag_indices_get(f1, f2, dl, diag_cs, isv, iev, jsv, jev)
  !Set the non-downsampled mask, it must be associated and initialized
  if (present(mask)) then
    locmask => mask
  elseif (associated(diag%axes%mask3d)) then
    locmask => diag%axes%mask3d
  else
    call MOM_error(FATAL, "downsample_diag_field_3d: Cannot downsample without a mask!!! ")
  endif

  call downsample_field(locfield, locfield_dsamp, dl, diag%xyz_method, locmask, diag_cs, diag, &
                        isv_o, jsv_o, isv, iev, jsv, jev)

end subroutine downsample_diag_field_3d

!> This subroutine allocates and computes a downsampled array from an input array
!! It also determines the diagnostics-compurte indices for the downsampled array
!! 2d interface
subroutine downsample_diag_field_2d(locfield, locfield_dsamp, dl, diag_cs, diag, isv, iev, jsv, jev, mask)
  real, dimension(:,:), pointer :: locfield !< Input array pointer
  real, dimension(:,:), allocatable, intent(inout) :: locfield_dsamp !< Output (downsampled) array
  type(diag_ctrl),   intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  type(diag_type),   intent(in) :: diag    !< A structure describing the diagnostic to post
  integer, intent(in) :: dl                !< Level of down sampling
  integer, intent(inout) :: isv            !< i-start index for diagnostics
  integer, intent(inout) :: iev            !< i-end index for diagnostics
  integer, intent(inout) :: jsv            !< j-start index for diagnostics
  integer, intent(inout) :: jev            !< j-end index for diagnostics
  real,    optional,target, intent(in) :: mask(:,:) !< If present, use this real array as the data mask.
  ! Locals
  real, dimension(:,:), pointer :: locmask
  integer :: f1,f2,isv_o,jsv_o

  locmask => NULL()
  !Get the correct indices corresponding to input field
  !Shape of the input diag field
  f1 = size(locfield,1)
  f2 = size(locfield,2)
  !Save the extents of the original (fine) domain
  isv_o = isv ; jsv_o = jsv
  !Get the shape of the downsampled field and overwrite isv,iev,jsv,jev with them
  call downsample_diag_indices_get(f1,f2, dl, diag_cs,isv,iev,jsv,jev)
  !Set the non-downsampled mask, it must be associated and initialized
  if (present(mask)) then
    locmask => mask
  elseif (associated(diag%axes%mask2d)) then
    locmask => diag%axes%mask2d
  else
    call MOM_error(FATAL, "downsample_diag_field_2d: Cannot downsample without a mask!!! ")
  endif

  call downsample_field(locfield, locfield_dsamp, dl, diag%xyz_method, locmask, diag_cs,diag, &
                      isv_o,jsv_o,isv,iev,jsv,jev)

end subroutine downsample_diag_field_2d

!> \section downsampling The down sample algorithm
!!
!! The down sample method could be deduced (before send_data call)
!!  from the diag%x_cell_method, diag%y_cell_method and diag%v_cell_method
!!
!! This is the summary of the down sample algoritm for a diagnostic field f:
!!  \f[
!!     f(Id,Jd) = \sum_{i,j} f(Id+i,Jd+j) * weight(Id+i,Jd+j) / [ \sum_{i,j} weight(Id+i,Jd+j)]
!!  \f]
!!  Here, i and j run from 0 to dl-1 (dl being the down sample level).
!!  Id,Jd are the down sampled (coarse grid) indices run over the coarsened compute grid,
!!  if and jf are the original (fine grid) indices.
!!
!! \verbatim
!! Example   x_cell y_cell v_cell algorithm_id    implemented weight(if,jf)
!! ---------------------------------------------------------------------------------------
!! theta     mean   mean   mean   MMM =222        G%areaT(if,jf)*h(if,jf)
!! u         point  mean   mean   PMM =022        dyCu(if,jf)*h(if,jf)*delta(if,Id)
!! v         mean   point  mean   MPM =202        dxCv(if,jf)*h(if,jf)*delta(jf,Jd)
!! ?         point  sum    mean   PSM =012        h(if,jf)*delta(if,Id)
!! volcello  sum    sum    sum    SSS =111        1
!! T_dfxy_co sum    sum    point  SSP =110        1
!! umo       point  sum    sum    PSS =011        1*delta(if,Id)
!! vmo       sum    point  sum    SPS =101        1*delta(jf,Jd)
!! umo_2d    point  sum    point  PSP =010        1*delta(if,Id)
!! vmo_2d    sum    point  point  SPP =100        1*delta(jf,Jd)
!! ?         point  mean   point  PMP =020        dyCu(if,jf)*delta(if,Id)
!! ?         mean   point  point  MPP =200        dxCv(if,jf)*delta(jf,Jd)
!! w         mean   mean   point  MMP =220        G%areaT(if,jf)
!! h*theta   mean   mean   sum    MMS =221        G%areaT(if,jf)
!!
!! delta is the Kronecker delta
!! \endverbatim

!> This subroutine allocates and computes a down sampled 3d array given an input array
!! The down sample method is based on the "cell_methods" for the diagnostics as explained
!! in the above table
subroutine downsample_field_3d(field_in, field_out, dl, method, mask, diag_cs, diag,isv_o,jsv_o,isv_d,iev_d,jsv_d,jev_d)
  real, dimension(:,:,:), pointer :: field_in      !< Original field to be down sampled
  real, dimension(:,:,:), allocatable :: field_out !< down sampled field
  integer, intent(in) :: dl                !< Level of down sampling
  integer,  intent(in) :: method           !< Sampling method
  real,  dimension(:,:,:), pointer :: mask !< Mask for field
  type(diag_ctrl), intent(in) :: diag_CS   !< Structure used to regulate diagnostic output
  type(diag_type), intent(in) :: diag      !< A structure describing the diagnostic to post
  integer, intent(in) :: isv_o             !< Original i-start index
  integer, intent(in) :: jsv_o             !< Original j-start index
  integer, intent(in) :: isv_d             !< i-start index of down sampled data
  integer, intent(in) :: iev_d             !< i-end index of down sampled data
  integer, intent(in) :: jsv_d             !< j-start index of down sampled data
  integer, intent(in) :: jev_d             !< j-end index of down sampled data
  ! Locals
  character(len=240) :: mesg
  integer :: i,j,ii,jj,i0,j0,f1,f2,f_in1,f_in2
  integer :: k,ks,ke
  real :: ave,total_weight,weight
  real :: eps_vol   ! A negligibly small volume or mass [H L2 ~> m3 or kg]
  real :: eps_area  ! A negligibly small area [L2 ~> m2]
  real :: eps_face  ! A negligibly small face area [H L ~> m2 or kg m-1]

  ks = 1 ; ke = size(field_in,3)
  eps_face = 1.0e-20 * diag_cs%G%US%m_to_L * diag_cs%GV%m_to_H
  eps_area = 1.0e-20 * diag_cs%G%US%m_to_L**2
  eps_vol = 1.0e-20 * diag_cs%G%US%m_to_L**2 * diag_cs%GV%m_to_H

  ! Allocate the down sampled field on the down sampled data domain
!  allocate(field_out(diag_cs%dsamp(dl)%isd:diag_cs%dsamp(dl)%ied,diag_cs%dsamp(dl)%jsd:diag_cs%dsamp(dl)%jed,ks:ke))
!  allocate(field_out(1:size(field_in,1)/dl,1:size(field_in,2)/dl,ks:ke))
  f_in1 = size(field_in,1)
  f_in2 = size(field_in,2)
  f1 = f_in1/dl
  f2 = f_in2/dl
  !Correction for the symmetric case
  if (diag_cs%G%symmetric) then
    f1 = f1 + mod(f_in1,dl)
    f2 = f2 + mod(f_in2,dl)
  endif
  allocate(field_out(1:f1,1:f2,ks:ke))

  ! Fill the down sampled field on the down sampled diagnostics (almost always compuate) domain
  if (method == MMM) then
    do k=ks,ke ; do j=jsv_d,jev_d ; do i=isv_d,iev_d
      i0 = isv_o+dl*(i-isv_d)
      j0 = jsv_o+dl*(j-jsv_d)
      ave = 0.0
      total_weight = 0.0
      do jj=j0,j0+dl-1 ; do ii=i0,i0+dl-1
!     do ii=i0,i0+dl-1 ; do jj=j0,j0+dl-1 !This seems to be faster!!!!
        weight = mask(ii,jj,k) * diag_cs%G%areaT(ii,jj) * diag_cs%h(ii,jj,k)
        total_weight = total_weight + weight
        ave = ave+field_in(ii,jj,k) * weight
      enddo ; enddo
      field_out(i,j,k)  = ave/(total_weight + eps_vol)  !Avoid zero mask at all aggregating cells where ave=0.0
    enddo ; enddo ; enddo
  elseif (method == SSS) then   !e.g., volcello
    do k=ks,ke ; do j=jsv_d,jev_d ; do i=isv_d,iev_d
      i0 = isv_o+dl*(i-isv_d)
      j0 = jsv_o+dl*(j-jsv_d)
      ave = 0.0
      do jj=j0,j0+dl-1 ; do ii=i0,i0+dl-1
        weight = mask(ii,jj,k)
        ave = ave+field_in(ii,jj,k)*weight
      enddo ; enddo
      field_out(i,j,k)  = ave !Masked Sum (total_weight=1)
    enddo ; enddo ; enddo
  elseif (method == MMP .or. method == MMS) then   !e.g., T_advection_xy
    do k=ks,ke ; do j=jsv_d,jev_d ; do i=isv_d,iev_d
      i0 = isv_o+dl*(i-isv_d)
      j0 = jsv_o+dl*(j-jsv_d)
      ave = 0.0
      total_weight = 0.0
      do jj=j0,j0+dl-1 ; do ii=i0,i0+dl-1
!     do ii=i0,i0+dl-1 ; do jj=j0,j0+dl-1
        weight = mask(ii,jj,k) * diag_cs%G%areaT(ii,jj)
        total_weight = total_weight + weight
        ave = ave+field_in(ii,jj,k)*weight
      enddo ; enddo
      field_out(i,j,k)  = ave / (total_weight+eps_area)  !Avoid zero mask at all aggregating cells where ave=0.0
    enddo ; enddo ; enddo
  elseif (method == PMM) then
    do k=ks,ke ; do j=jsv_d,jev_d ; do i=isv_d,iev_d
      i0 = isv_o+dl*(i-isv_d)
      j0 = jsv_o+dl*(j-jsv_d)
      ave = 0.0
      total_weight = 0.0
      ii=i0
      do jj=j0,j0+dl-1
        weight = mask(ii,jj,k) * diag_cs%G%dyCu(ii,jj) * diag_cs%h(ii,jj,k)
        total_weight = total_weight +weight
        ave = ave+field_in(ii,jj,k)*weight
      enddo
      field_out(i,j,k)  = ave/(total_weight+eps_face)  !Avoid zero mask at all aggregating cells where ave=0.0
    enddo ; enddo ; enddo
  elseif (method == PSS) then    !e.g. umo
    do k=ks,ke ; do j=jsv_d,jev_d ; do i=isv_d,iev_d
      i0 = isv_o+dl*(i-isv_d)
      j0 = jsv_o+dl*(j-jsv_d)
      ave = 0.0
      ii=i0
      do jj=j0,j0+dl-1
        weight = mask(ii,jj,k)
        ave = ave+field_in(ii,jj,k)*weight
      enddo
      field_out(i,j,k)  = ave  !Masked Sum (total_weight=1)
    enddo ; enddo ; enddo
  elseif (method == SPS) then   !e.g. vmo
    do k=ks,ke ; do j=jsv_d,jev_d ; do i=isv_d,iev_d
      i0 = isv_o+dl*(i-isv_d)
      j0 = jsv_o+dl*(j-jsv_d)
      ave = 0.0
      jj=j0
      do ii=i0,i0+dl-1
        weight = mask(ii,jj,k)
        ave = ave+field_in(ii,jj,k)*weight
      enddo
      field_out(i,j,k)  = ave  !Masked Sum (total_weight=1)
    enddo ; enddo ; enddo
  elseif (method == MPM) then
    do k=ks,ke ; do j=jsv_d,jev_d ; do i=isv_d,iev_d
      i0 = isv_o+dl*(i-isv_d)
      j0 = jsv_o+dl*(j-jsv_d)
      ave = 0.0
      total_weight = 0.0
      jj=j0
      do ii=i0,i0+dl-1
        weight = mask(ii,jj,k) * diag_cs%G%dxCv(ii,jj) * diag_cs%h(ii,jj,k)
        total_weight = total_weight + weight
        ave = ave+field_in(ii,jj,k)*weight
      enddo
      field_out(i,j,k)  = ave/(total_weight+eps_face)  !Avoid zero mask at all aggregating cells where ave=0.0
    enddo ; enddo ; enddo
  elseif (method == MSK) then !The input field is a mask, subsample
    field_out(:,:,:) = 0.0
    do k=ks,ke ; do j=jsv_d,jev_d ; do i=isv_d,iev_d
      i0 = isv_o+dl*(i-isv_d)
      j0 = jsv_o+dl*(j-jsv_d)
      ave = 0.0
      do jj=j0,j0+dl-1 ; do ii=i0,i0+dl-1
        ave = ave+field_in(ii,jj,k)
      enddo ; enddo
      if (ave > 0.0) field_out(i,j,k)=1.0
    enddo ; enddo ; enddo
  else
    write (mesg,*) " unknown sampling method: ",method
    call MOM_error(FATAL, "downsample_field_3d: "//trim(mesg)//" "//trim(diag%debug_str))
  endif

end subroutine downsample_field_3d

!> This subroutine allocates and computes a down sampled 2d array given an input array
!! The down sample method is based on the "cell_methods" for the diagnostics as explained
!! in the above table
subroutine downsample_field_2d(field_in, field_out, dl, method, mask, diag_cs, diag, &
                               isv_o, jsv_o, isv_d, iev_d, jsv_d, jev_d)
  real, dimension(:,:), pointer :: field_in      !< Original field to be down sampled
  real, dimension(:,:), allocatable :: field_out !< Down sampled field
  integer, intent(in) :: dl                !< Level of down sampling
  integer,  intent(in) :: method           !< Sampling method
  real, dimension(:,:), pointer :: mask    !< Mask for field
  type(diag_ctrl),   intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  type(diag_type),   intent(in) :: diag    !< A structure describing the diagnostic to post
  integer, intent(in) :: isv_o             !< Original i-start index
  integer, intent(in) :: jsv_o             !< Original j-start index
  integer, intent(in) :: isv_d             !< i-start index of down sampled data
  integer, intent(in) :: iev_d             !< i-end index of down sampled data
  integer, intent(in) :: jsv_d             !< j-start index of down sampled data
  integer, intent(in) :: jev_d             !< j-end index of down sampled data
  ! Locals
  character(len=240) :: mesg
  integer :: i,j,ii,jj,i0,j0,f1,f2,f_in1,f_in2
  real :: ave, total_weight, weight
  real :: epsilon = 1.0e-20  ! A negligibly small count of weights [nondim]
  real :: eps_area  ! A negligibly small area [L2 ~> m2]
  real :: eps_len   ! A negligibly small horizontal length [L ~> m]

  eps_len = 1.0e-20 * diag_cs%G%US%m_to_L
  eps_area = 1.0e-20 * diag_cs%G%US%m_to_L**2

  ! Allocate the down sampled field on the down sampled data domain
!  allocate(field_out(diag_cs%dsamp(dl)%isd:diag_cs%dsamp(dl)%ied,diag_cs%dsamp(dl)%jsd:diag_cs%dsamp(dl)%jed))
!  allocate(field_out(1:size(field_in,1)/dl,1:size(field_in,2)/dl))
  ! Fill the down sampled field on the down sampled diagnostics (almost always compuate) domain
  f_in1 = size(field_in,1)
  f_in2 = size(field_in,2)
  f1 = f_in1/dl
  f2 = f_in2/dl
  ! Correction for the symmetric case
  if (diag_cs%G%symmetric) then
    f1 = f1 + mod(f_in1,dl)
    f2 = f2 + mod(f_in2,dl)
  endif
  allocate(field_out(1:f1,1:f2))

  if (method == MMP) then
    do j=jsv_d,jev_d ; do i=isv_d,iev_d
      i0 = isv_o+dl*(i-isv_d)
      j0 = jsv_o+dl*(j-jsv_d)
      ave = 0.0
      total_weight = 0.0
      do jj=j0,j0+dl-1 ; do ii=i0,i0+dl-1
!      do ii=i0,i0+dl-1 ; do jj=j0,j0+dl-1
        weight = mask(ii,jj)*diag_cs%G%areaT(ii,jj)
        total_weight = total_weight + weight
        ave = ave+field_in(ii,jj)*weight
      enddo ; enddo
      field_out(i,j) = ave/(total_weight + eps_area)  !Avoid zero mask at all aggregating cells where ave=0.0
    enddo ; enddo
  elseif (method == SSP) then    ! e.g., T_dfxy_cont_tendency_2d
    do j=jsv_d,jev_d ; do i=isv_d,iev_d
      i0 = isv_o+dl*(i-isv_d)
      j0 = jsv_o+dl*(j-jsv_d)
      ave = 0.0
      total_weight = 0.0
      do jj=j0,j0+dl-1 ; do ii=i0,i0+dl-1
!      do ii=i0,i0+dl-1 ; do jj=j0,j0+dl-1
        weight = mask(ii,jj)
        total_weight = total_weight + weight
        ave = ave+field_in(ii,jj)*weight
      enddo ; enddo
      field_out(i,j) = ave/(total_weight+epsilon)  !Avoid zero mask at all aggregating cells where ave=0.0
    enddo ; enddo
  elseif (method == PSP) then   ! e.g., umo_2d
    do j=jsv_d,jev_d ; do i=isv_d,iev_d
      i0 = isv_o+dl*(i-isv_d)
      j0 = jsv_o+dl*(j-jsv_d)
      ave = 0.0
      total_weight = 0.0
      ii=i0
      do jj=j0,j0+dl-1
        weight = mask(ii,jj)
        total_weight = total_weight +weight
        ave = ave+field_in(ii,jj)*weight
      enddo
      field_out(i,j) = ave/(total_weight+epsilon)  !Avoid zero mask at all aggregating cells where ave=0.0
    enddo ; enddo
  elseif (method == SPP) then   ! e.g., vmo_2d
    do j=jsv_d,jev_d ; do i=isv_d,iev_d
      i0 = isv_o+dl*(i-isv_d)
      j0 = jsv_o+dl*(j-jsv_d)
      ave = 0.0
      total_weight = 0.0
      jj=j0
      do ii=i0,i0+dl-1
        weight = mask(ii,jj)
        total_weight = total_weight +weight
        ave = ave+field_in(ii,jj)*weight
      enddo
      field_out(i,j) = ave/(total_weight+epsilon)  !Avoid zero mask at all aggregating cells where ave=0.0
    enddo ; enddo
  elseif (method == PMP) then
    do j=jsv_d,jev_d ; do i=isv_d,iev_d
      i0 = isv_o+dl*(i-isv_d)
      j0 = jsv_o+dl*(j-jsv_d)
      ave = 0.0
      total_weight = 0.0
      ii=i0
      do jj=j0,j0+dl-1
        weight = mask(ii,jj) * diag_cs%G%dyCu(ii,jj)!*diag_cs%h(ii,jj,1) !Niki?
        total_weight = total_weight +weight
        ave = ave+field_in(ii,jj)*weight
      enddo
      field_out(i,j) = ave/(total_weight+eps_len)  !Avoid zero mask at all aggregating cells where ave=0.0
    enddo ; enddo
  elseif (method == MPP) then
    do j=jsv_d,jev_d ; do i=isv_d,iev_d
      i0 = isv_o+dl*(i-isv_d)
      j0 = jsv_o+dl*(j-jsv_d)
      ave = 0.0
      total_weight = 0.0
      jj=j0
      do ii=i0,i0+dl-1
        weight = mask(ii,jj)* diag_cs%G%dxCv(ii,jj)!*diag_cs%h(ii,jj,1) !Niki?
        total_weight = total_weight +weight
        ave = ave+field_in(ii,jj)*weight
      enddo
      field_out(i,j) = ave/(total_weight+eps_len)  !Avoid zero mask at all aggregating cells where ave=0.0
    enddo ; enddo
  elseif (method == MSK) then !The input field is a mask, subsample
    field_out(:,:) = 0.0
    do j=jsv_d,jev_d ; do i=isv_d,iev_d
      i0 = isv_o+dl*(i-isv_d)
      j0 = jsv_o+dl*(j-jsv_d)
      ave = 0.0
      do jj=j0,j0+dl-1 ; do ii=i0,i0+dl-1
        ave = ave+field_in(ii,jj)
      enddo ; enddo
      if (ave > 0.0) field_out(i,j)=1.0
    enddo ; enddo
  else
    write (mesg,*) " unknown sampling method: ",method
    call MOM_error(FATAL, "downsample_field_2d: "//trim(mesg)//" "//trim(diag%debug_str))
  endif

end subroutine downsample_field_2d

!> Allocate and compute the 2d down sampled mask
!! The masks are down sampled based on a minority rule, i.e., a coarse cell is open (1)
!! if at least one of the sub-cells are open, otherwise it's closed (0)
subroutine downsample_mask_2d(field_in, field_out, dl, isc_o, jsc_o, isc_d, iec_d, jsc_d, jec_d, &
                              isd_d, ied_d, jsd_d, jed_d)
  real, dimension(:,:), intent(in) :: field_in !< Original field to be down sampled
  real, dimension(:,:), pointer :: field_out   !< Down sampled field
  integer, intent(in) :: dl    !< Level of down sampling
  integer, intent(in) :: isc_o !< Original i-start index
  integer, intent(in) :: jsc_o !< Original j-start index
  integer, intent(in) :: isc_d !< Computational i-start index of down sampled data
  integer, intent(in) :: iec_d !< Computational i-end index of down sampled data
  integer, intent(in) :: jsc_d !< Computational j-start index of down sampled data
  integer, intent(in) :: jec_d !< Computational j-end index of down sampled data
  integer, intent(in) :: isd_d !< Computational i-start index of down sampled data
  integer, intent(in) :: ied_d !< Computational i-end index of down sampled data
  integer, intent(in) :: jsd_d !< Computational j-start index of down sampled data
  integer, intent(in) :: jed_d !< Computational j-end index of down sampled data
  ! Locals
  integer :: i,j,ii,jj,i0,j0
  real    :: tot_non_zero
  ! down sampled mask = 0 unless the mask value of one of the down sampling cells is 1
  allocate(field_out(isd_d:ied_d,jsd_d:jed_d))
  field_out(:,:) = 0.0
  do j=jsc_d,jec_d ; do i=isc_d,iec_d
    i0 = isc_o+dl*(i-isc_d)
    j0 = jsc_o+dl*(j-jsc_d)
    tot_non_zero = 0.0
    do jj=j0,j0+dl-1 ; do ii=i0,i0+dl-1
      tot_non_zero = tot_non_zero + field_in(ii,jj)
    enddo ; enddo
    if (tot_non_zero > 0.0) field_out(i,j)=1.0
  enddo ; enddo
end subroutine downsample_mask_2d

!> Allocate and compute the 3d down sampled mask
!! The masks are down sampled based on a minority rule, i.e., a coarse cell is open (1)
!! if at least one of the sub-cells are open, otherwise it's closed (0)
subroutine downsample_mask_3d(field_in, field_out, dl, isc_o, jsc_o, isc_d, iec_d, jsc_d, jec_d, &
                              isd_d, ied_d, jsd_d, jed_d)
  real, dimension(:,:,:), intent(in) :: field_in !< Original field to be down sampled
  real, dimension(:,:,:), pointer :: field_out   !< down sampled field
  integer, intent(in) :: dl    !< Level of down sampling
  integer, intent(in) :: isc_o !< Original i-start index
  integer, intent(in) :: jsc_o !< Original j-start index
  integer, intent(in) :: isc_d !< Computational i-start index of down sampled data
  integer, intent(in) :: iec_d !< Computational i-end index of down sampled data
  integer, intent(in) :: jsc_d !< Computational j-start index of down sampled data
  integer, intent(in) :: jec_d !< Computational j-end index of down sampled data
  integer, intent(in) :: isd_d !< Computational i-start index of down sampled data
  integer, intent(in) :: ied_d !< Computational i-end index of down sampled data
  integer, intent(in) :: jsd_d !< Computational j-start index of down sampled data
  integer, intent(in) :: jed_d !< Computational j-end index of down sampled data
  ! Locals
  integer :: i,j,ii,jj,i0,j0,k,ks,ke
  real    :: tot_non_zero
  ! down sampled mask = 0 unless the mask value of one of the down sampling cells is 1
  ks = lbound(field_in,3) ; ke = ubound(field_in,3)
  allocate(field_out(isd_d:ied_d,jsd_d:jed_d,ks:ke))
  field_out(:,:,:) = 0.0
  do k=ks,ke ; do j=jsc_d,jec_d ; do i=isc_d,iec_d
    i0 = isc_o+dl*(i-isc_d)
    j0 = jsc_o+dl*(j-jsc_d)
    tot_non_zero = 0.0
    do jj=j0,j0+dl-1 ; do ii=i0,i0+dl-1
      tot_non_zero = tot_non_zero + field_in(ii,jj,k)
    enddo ; enddo
    if (tot_non_zero > 0.0) field_out(i,j,k)=1.0
  enddo ; enddo ; enddo
end subroutine downsample_mask_3d

!> Fakes a register of a diagnostic to find out if an obsolete
!! parameter appears in the diag_table.
logical function found_in_diagtable(diag, varName)
  type(diag_ctrl),            intent(in) :: diag       !< A structure used to control diagnostics.
  character(len=*),           intent(in) :: varName    !< The obsolete diagnostic name
  ! Local
  integer :: handle ! Integer handle returned from diag_manager

  ! We use register_static_field_fms() instead of register_static_field() so
  ! that the diagnostic does not appear in the available diagnostics list.
  handle = register_static_field_infra('ocean_model', varName, diag%axesT1%handles)

  found_in_diagtable = (handle>0)

end function found_in_diagtable

end module MOM_diag_mediator
