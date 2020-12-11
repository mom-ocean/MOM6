!> Convenient wrappers to the FMS diag_manager interfaces with additional diagnostic capabilies.
module MOM_IS_diag_mediator

! This file is a part of SIS2. See LICENSE.md for the license.

use MOM_grid, only : ocean_grid_type

use MOM_coms, only : PE_here
use MOM_error_handler, only : MOM_error, FATAL, is_root_pe, assert
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_safe_alloc, only : safe_alloc_ptr, safe_alloc_alloc
use MOM_string_functions, only : lowercase, uppercase, slasher
use MOM_time_manager, only : time_type

use diag_manager_mod, only : diag_manager_init
use diag_manager_mod, only : send_data, diag_axis_init,EAST,NORTH
use diag_manager_mod, only : register_diag_field_fms=>register_diag_field
use diag_manager_mod, only : register_static_field_fms=>register_static_field

implicit none ; private

public diag_mediator_infrastructure_init
public set_axes_info, post_data, register_MOM_IS_diag_field, time_type
public safe_alloc_ptr, safe_alloc_alloc
public enable_averaging, disable_averaging, query_averaging_enabled
public enable_averages
public diag_mediator_init, diag_mediator_end, set_diag_mediator_grid
public diag_mediator_close_registration, get_diag_time_end
public diag_axis_init, register_static_field

!> Make a diagnostic available for averaging or output.
!interface post_data
!   module procedure post_data_2d
!end interface post_data

!> 2D/3D axes type to contain 1D axes handles and pointers to masks
type, public :: axesType
  character(len=15) :: id   !< The id string for this particular combination of handles.
  integer           :: rank !< Number of dimensions in the list of axes.
  integer, dimension(:), allocatable :: handles !< Handles to 1D axes.
  type(diag_ctrl), pointer :: diag_cs => null() !< A structure that is used to regulate diagnostic output
end type axesType

!> This type is used to represent a diagnostic at the diag_mediator level.
type, private :: diag_type
  logical :: in_use              !< This diagnostic is in use
  integer :: fms_diag_id         !< underlying FMS diag id
  character(len=24) :: name      !< The diagnostic name
  real :: conversion_factor = 0. !< A factor to multiply data by before posting to FMS, if non-zero.
  real, pointer, dimension(:,:)   :: mask2d => null()      !< A 2-d mask on the data domain for this diagnostic
  real, pointer, dimension(:,:)   :: mask2d_comp => null() !< A 2-d mask on the computational domain for this diagnostic
end type diag_type

!>   The SIS_diag_ctrl data type contains times to regulate diagnostics along with masks and
!! axes to use with diagnostics, and a list of structures with data about each diagnostic.
type, public :: diag_ctrl
  integer :: doc_unit = -1 !< The unit number of a diagnostic documentation file.
                           !! This file is open if doc_unit is > 0.

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
  real :: time_int              !< The time interval in s for any fields that are offered for averaging.
  type(time_type) :: time_end   !< The end time of the valid interval for any offered field.
  logical :: ave_enabled = .false. !< .true. if averaging is enabled.

  !>@{ The following are 3D and 2D axis groups defined for output.  The names indicate
  !! the horizontal locations (B, T, Cu, or Cv), vertical locations (L, i, or 1) and
  !! thickness categories (c, c0, or 1).
  type(axesType) :: axesBL, axesTL, axesCuL, axesCvL
  type(axesType) :: axesBi, axesTi, axesCui, axesCvi
  type(axesType) :: axesBc, axesTc, axesCuc, axesCvc
  type(axesType) :: axesBc0, axesTc0, axesCuc0, axesCvc0
  type(axesType) :: axesB1, axesT1, axesCu1, axesCv1
  !!@}

  ! Mask arrays for diagnostics
  real, dimension(:,:),   pointer :: mask2dT   => null() !< 2D mask array for cell-center points
  real, dimension(:,:),   pointer :: mask2dBu  => null() !< 2D mask array for cell-corners
  real, dimension(:,:),   pointer :: mask2dCu  => null() !< 2D mask array for east-faces
  real, dimension(:,:),   pointer :: mask2dCv  => null() !< 2D mask array for north-faces
  !> Computational domain mask arrays for diagnostics.
  real, dimension(:,:),   pointer :: mask2dT_comp => null()

#define DIAG_ALLOC_CHUNK_SIZE 15
  type(diag_type), dimension(:), allocatable :: diags !< The array of diagnostics
  integer :: next_free_diag_id !< The next unused diagnostic ID
  !> default missing value to be sent to ALL diagnostics registerations
  real :: missing_value = -1.0e34

end type diag_ctrl

contains

!> Set up the grid and axis information for use by the ice shelf model.
subroutine set_axes_info(G, param_file, diag_cs, axes_set_name)
  type(ocean_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl),     intent(inout) :: diag_cs !< A structure that is used to regulate diagnostic output
  character(len=*), optional, intent(in) :: axes_set_name !<  A name to use for this set of axes.
                                                !! The default is "ice".
!   This subroutine sets up the grid and axis information for use by the ice shelf model.

  ! Local variables
  integer :: id_xq, id_yq, id_zl, id_zi, id_xh, id_yh, id_ct, id_ct0
  integer :: k
  logical :: Cartesian_grid
  character(len=80) :: grid_config, units_temp, set_name
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_IS_diag_mediator" ! This module's name.

  set_name = "ice_shelf" ; if (present(axes_set_name)) set_name = trim(axes_set_name)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version)
  call get_param(param_file, mdl, "GRID_CONFIG", grid_config, &
                 "The method for defining the horizontal grid.  Valid "//&
                 "entries include:\n"//&
                 "\t file - read the grid from GRID_FILE \n"//&
                 "\t mosaic - read the grid from a mosaic grid file \n"//&
                 "\t cartesian - a Cartesian grid \n"//&
                 "\t spherical - a spherical grid \n"//&
                 "\t mercator  - a Mercator grid", fail_if_missing=.true.)

  G%x_axis_units = "degrees_E"
  G%y_axis_units = "degrees_N"
  if (index(lowercase(trim(grid_config)),"cartesian") > 0) then
    ! This is a cartesian grid, and may have different axis units.
    Cartesian_grid = .true.
    call get_param(param_file, mdl, "AXIS_UNITS", units_temp, &
                 "The units for the x- and y- axis labels.  AXIS_UNITS "//&
                 "should be defined as 'k' for km, 'm' for m, or 'd' "//&
                 "for degrees of latitude and longitude (the default). "//&
                 "Except on a Cartesian grid, only degrees are currently "//&
                 "implemented.", default='degrees')
    if (units_temp(1:1) == 'k') then
      G%x_axis_units = "kilometers" ; G%y_axis_units = "kilometers"
    elseif (units_temp(1:1) == 'm') then
      G%x_axis_units = "meters" ; G%y_axis_units = "meters"
    endif
    call log_param(param_file, mdl, "explicit AXIS_UNITS", G%x_axis_units)
  else
    Cartesian_grid = .false.
  endif

  if (G%symmetric) then
     id_xq = diag_axis_init('xB', G%gridLonB(G%isgB:G%iegB), G%x_axis_units, 'x', &
          'Boundary point nominal longitude',set_name=set_name, &
          Domain2=G%Domain%mpp_domain, domain_position=EAST)
     id_yq = diag_axis_init('yB', G%gridLatB(G%jsgB:G%jegB), G%y_axis_units, 'y', &
          'Boundary point nominal latitude', set_name=set_name, &
          Domain2=G%Domain%mpp_domain, domain_position=NORTH)

  else
     id_xq = diag_axis_init('xB', G%gridLonB(G%isg:G%ieg), G%x_axis_units, 'x', &
          'Boundary point nominal longitude',set_name=set_name, &
          Domain2=G%Domain%mpp_domain, domain_position=EAST)
     id_yq = diag_axis_init('yB', G%gridLatB(G%jsg:G%jeg), G%y_axis_units, 'y', &
          'Boundary point nominal latitude', set_name=set_name, &
          Domain2=G%Domain%mpp_domain, domain_position=NORTH)

  endif
  id_xh = diag_axis_init('xT', G%gridLonT(G%isg:G%ieg), G%x_axis_units, 'x', &
          'T point nominal longitude', set_name=set_name, &
          Domain2=G%Domain%mpp_domain)
  id_yh = diag_axis_init('yT', G%gridLatT(G%jsg:G%jeg), G%y_axis_units, 'y', &
          'T point nominal latitude', set_name=set_name, &
          Domain2=G%Domain%mpp_domain)

  ! Axis groupings for 2-D arrays.
  call defineAxes(diag_cs, [id_xh, id_yh], diag_cs%axesT1)
  call defineAxes(diag_cs, [id_xq, id_yq], diag_cs%axesB1)
  call defineAxes(diag_cs, [id_xq, id_yh], diag_cs%axesCu1)
  call defineAxes(diag_cs, [id_xh, id_yq], diag_cs%axesCv1)

end subroutine set_axes_info

!> Define an a group of axes from a list of handles
subroutine defineAxes(diag_cs, handles, axes)
  ! Defines "axes" from list of handle and associates mask
  type(diag_ctrl), target, intent(in)  :: diag_cs !< A structure that is used to regulate diagnostic output
  integer, dimension(:),       intent(in)  :: handles !< A set of axis handles that define the axis group
  type(axesType),              intent(out) :: axes    !< A group of axes that is set up here

  ! Local variables
  integer :: n
  n = size(handles)
  if (n<1 .or. n>3) call MOM_error(FATAL,"defineAxes: wrong size for list of handles!")
  allocate( axes%handles(n) )
  axes%id = i2s(handles, n) ! Identifying string
  axes%rank = n
  axes%handles(:) = handles(:)
  axes%diag_cs => diag_cs ! A (circular) link back to the MOM_IS_diag_ctrl structure
end subroutine defineAxes

!> Set up the current grid for the diag mediator
subroutine set_diag_mediator_grid(G, diag_cs)
  type(ocean_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(diag_ctrl),     intent(inout) :: diag_cs !< A structure that is used to regulate diagnostic output

  diag_cs%is = G%isc - (G%isd-1) ; diag_cs%ie = G%iec - (G%isd-1)
  diag_cs%js = G%jsc - (G%jsd-1) ; diag_cs%je = G%jec - (G%jsd-1)
  diag_cs%isd = G%isd ; diag_cs%ied = G%ied ; diag_cs%jsd = G%jsd ; diag_cs%jed = G%jed
end subroutine set_diag_mediator_grid

!> Offer a 2d diagnostic field for output or averaging
subroutine post_data(diag_field_id, field, diag_cs, is_static, mask)
  integer,           intent(in) :: diag_field_id !< the id for an output variable returned by a
                                              !! previous call to register_diag_field.
  real,    target,   intent(in) :: field(:,:) !< The 2-d array being offered for output or averaging.
  type(diag_ctrl), target, &
                     intent(in) :: diag_cs !< A structure that is used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.
  logical, optional, intent(in) :: mask(:,:) !< If present, use this logical array as the data mask.

  ! Local variables
  real, dimension(:,:), pointer :: locfield
  logical :: used, is_stat
  logical :: i_data, j_data
  integer :: isv, iev, jsv, jev, i, j
  integer :: fms_diag_id
  type(diag_type), pointer :: diag => NULL()

  locfield => NULL()
  is_stat = .false. ; if (present(is_static)) is_stat = is_static

  ! Get a pointer to the diag type for this id, and the FMS-level diag id.
  call assert(diag_field_id < diag_cs%next_free_diag_id, &
              'post_data: Unregistered diagnostic id')
  diag => diag_cs%diags(diag_field_id)
  fms_diag_id = diag%fms_diag_id

  ! Determine the proper array indices, noting that because of the (:,:)
  ! declaration of field, symmetric arrays are using a SW-grid indexing,
  ! but non-symmetric arrays are using a NE-grid indexing.  Send_data
  ! actually only uses the difference between ie and is to determine
  ! the output data size and assumes that halos are symmetric.
  isv = diag_cs%is ; iev = diag_cs%ie ; jsv = diag_cs%js ; jev = diag_cs%je

  if ( size(field,1) == diag_cs%ied-diag_cs%isd +1 ) then
    isv = diag_cs%is ; iev = diag_cs%ie ; i_data = .true.   ! Data domain
  elseif ( size(field,1) == diag_cs%ied-diag_cs%isd +2 ) then
    isv = diag_cs%is ; iev = diag_cs%ie+1 ; i_data = .true. ! Symmetric data domain
  elseif ( size(field,1) == diag_cs%ie-diag_cs%is +1 ) then
    isv = 1 ; iev = diag_cs%ie + 1-diag_cs%is ; i_data = .false. ! Computational domain
  elseif ( size(field,1) == diag_cs%ie-diag_cs%is +2 ) then
    isv = 1 ; iev = diag_cs%ie + 2-diag_cs%is ; i_data = .false. ! Symmetric computational domain
  else
    call MOM_error(FATAL,"post_MOM_IS_data_2d: peculiar size in i-direction of "//trim(diag%name))
  endif
  if ( size(field,2) == diag_cs%jed-diag_cs%jsd +1 ) then
    jsv = diag_cs%js ; jev = diag_cs%je ; j_data = .true.   ! Data domain
  elseif ( size(field,2) == diag_cs%jed-diag_cs%jsd +2 ) then
    jsv = diag_cs%js ; jev = diag_cs%je+1 ; j_data = .true. ! Symmetric data domain
  elseif ( size(field,2) == diag_cs%je-diag_cs%js +1 ) then
    jsv = 1 ; jev = diag_cs%je + 1-diag_cs%js ; j_data = .false. ! Computational domain
  elseif ( size(field,1) == diag_cs%je-diag_cs%js +2 ) then
    jsv = 1 ; jev = diag_cs%je + 2-diag_cs%js ; j_data = .false. ! Symmetric computational domain
  else
    call MOM_error(FATAL,"post_MOM_IS_data_2d: peculiar size in j-direction "//trim(diag%name))
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

  if (present(mask)) then
    if ((size(field,1) /= size(mask,1)) .or. &
        (size(field,2) /= size(mask,2))) then
      call MOM_error(FATAL, "post_MOM_IS_data_2d: post_MOM_IS_data called with a mask "//&
                            "that does not match the size of field "//trim(diag%name))
    endif
  elseif ( i_data .NEQV. j_data ) then
    call MOM_error(FATAL, "post_MOM_IS_data_2d: post_MOM_IS_data called for "//&
                   trim(diag%name)//" with mixed computational and data domain array sizes.")
  endif

  if (is_stat) then
    if (present(mask)) then
      used = send_data(fms_diag_id, locfield, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, mask=mask)
    elseif(i_data .and. associated(diag%mask2d)) then
      used = send_data(fms_diag_id, locfield, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, rmask=diag%mask2d)
    elseif((.not.i_data) .and. associated(diag%mask2d_comp)) then
      used = send_data(fms_diag_id, locfield, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, rmask=diag%mask2d_comp)
    else
      used = send_data(fms_diag_id, locfield, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev)
    endif
  elseif (diag_cs%ave_enabled) then
    if (present(mask)) then
      used = send_data(fms_diag_id, locfield, diag_cs%time_end, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                       weight=diag_cs%time_int, mask=mask)
    elseif(i_data .and. associated(diag%mask2d)) then
      used = send_data(fms_diag_id, locfield, diag_cs%time_end, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                       weight=diag_cs%time_int, rmask=diag%mask2d)
    elseif((.not.i_data) .and. associated(diag%mask2d_comp)) then
      used = send_data(fms_diag_id, locfield, diag_cs%time_end, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                       weight=diag_cs%time_int, rmask=diag%mask2d_comp)
    else
      used = send_data(fms_diag_id, locfield, diag_cs%time_end, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                       weight=diag_cs%time_int)
    endif
  endif

  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.) ) deallocate( locfield )

end subroutine post_data


!> Enable the accumulation of time averages over the specified time interval.
subroutine enable_averaging(time_int_in, time_end_in, diag_cs)
  real,                intent(in)    :: time_int_in !< The time interval over which any values
!                                                   !! that are offered are valid [s].
  type(time_type),     intent(in)    :: time_end_in !< The end time of the valid interval.
  type(diag_ctrl), intent(inout) :: diag_cs !< A structure that is used to regulate diagnostic output
! This subroutine enables the accumulation of time averages over the
! specified time interval.

!  if (num_file==0) return
  diag_cs%time_int = time_int_in
  diag_cs%time_end = time_end_in
  diag_cs%ave_enabled = .true.
end subroutine enable_averaging

! Put a block on averaging any offered fields.
subroutine disable_averaging(diag_cs)
  type(diag_ctrl), intent(inout) :: diag_cs !< A structure that is used to regulate diagnostic output

  diag_cs%time_int = 0.0
  diag_cs%ave_enabled = .false.

end subroutine disable_averaging

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
!  elseif (associated(diag_CS%US)) then
!    diag_cs%time_int = time_int*diag_CS%US%T_to_s
  else
    diag_cs%time_int = time_int
  endif
  diag_cs%time_end = time_end
  diag_cs%ave_enabled = .true.
end subroutine enable_averages

!> Indicate whether averaging diagnostics is currently enabled
logical function query_averaging_enabled(diag_cs, time_int, time_end)
  type(diag_ctrl),           intent(in)  :: diag_cs !< A structure that is used to regulate diagnostic output
  real,            optional, intent(out) :: time_int !< The current setting of diag_cs%time_int [s].
  type(time_type), optional, intent(out) :: time_end !< The current setting of diag_cs%time_end.

  if (present(time_int)) time_int = diag_cs%time_int
  if (present(time_end)) time_end = diag_cs%time_end
  query_averaging_enabled = diag_cs%ave_enabled
end function query_averaging_enabled

subroutine diag_mediator_infrastructure_init(err_msg)
  ! This subroutine initializes the FMS diag_manager.
  character(len=*), optional, intent(out)   :: err_msg !< An error message

  call diag_manager_init(err_msg=err_msg)
end subroutine diag_mediator_infrastructure_init

!> diag_mediator_init initializes the MOM diag_mediator and opens the available

!> Return the currently specified valid end time for diagnostics
function get_diag_time_end(diag_cs)
  type(diag_ctrl),           intent(in)  :: diag_cs !< A structure that is used to regulate diagnostic output
  type(time_type) :: get_diag_time_end

!   This function returns the valid end time for diagnostics that are handled
! outside of the MOM6 infrastructure, such as via the generic tracer code.

  get_diag_time_end = diag_cs%time_end
end function get_diag_time_end

!> Returns the "MOM_IS_diag_mediator" handle for a group of diagnostics derived from one field.
function register_MOM_IS_diag_field(module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, mask_variant, standard_name, &
     verbose, do_not_log, err_msg, interp_method, tile_count, conversion) result (register_diag_field)
  integer :: register_diag_field  !< The returned diagnostic handle
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ice_model"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axesType),   intent(in) :: axes       !< The axis group for this field
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
  integer,          optional, intent(in) :: tile_count   !< no clue (not used in MOM_IS?)
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to file

  ! Local variables
  character(len=240) :: mesg
  real :: MOM_missing_value
  integer :: primary_id, fms_id
  type(diag_ctrl), pointer :: diag_cs => NULL() ! A structure that is used
                                               ! to regulate diagnostic output
  type(diag_type), pointer :: diag => NULL()

  MOM_missing_value = axes%diag_cs%missing_value
  if(present(missing_value)) MOM_missing_value = missing_value

  diag_cs => axes%diag_cs
  primary_id = -1

  fms_id = register_diag_field_fms(module_name, field_name, axes%handles, &
         init_time, long_name=long_name, units=units, missing_value=MOM_missing_value, &
         range=range, mask_variant=mask_variant, standard_name=standard_name, &
         verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
         interp_method=interp_method, tile_count=tile_count)
  if (fms_id > 0) then
    primary_id = get_new_diag_id(diag_cs)
    diag => diag_cs%diags(primary_id)
    diag%fms_diag_id = fms_id
    if (len(field_name) > len(diag%name)) then
      diag%name = field_name(1:len(diag%name))
    else ; diag%name = field_name ; endif

    if (present(conversion)) diag%conversion_factor = conversion
  endif

  if (is_root_pe() .and. diag_CS%doc_unit > 0) then
    if (primary_id > 0) then
       mesg = '"'//trim(module_name)//'", "'//trim(field_name)//'"  [Used]'
    else
       mesg = '"'//trim(module_name)//'", "'//trim(field_name)//'"  [Unused]'
    endif
    write(diag_CS%doc_unit, '(a)') trim(mesg)
    if (present(long_name)) call describe_option("long_name", long_name, diag_CS)
    if (present(units)) call describe_option("units", units, diag_CS)
    if (present(standard_name)) &
      call describe_option("standard_name", standard_name, diag_CS)
  endif

  !Decide what mask to use based on the axes info
  if (primary_id > 0) then
  !3d masks
    !2d masks
    if (axes%rank == 2) then
      diag%mask2d => null() ; diag%mask2d_comp => null()
      if (axes%id == diag_cs%axesT1%id) then
        diag%mask2d =>  diag_cs%mask2dT
        diag%mask2d_comp => diag_cs%mask2dT_comp
      elseif (axes%id == diag_cs%axesB1%id) then
        diag%mask2d =>  diag_cs%mask2dBu
      elseif (axes%id == diag_cs%axesCu1%id) then
        diag%mask2d =>  diag_cs%mask2dCu
      elseif (axes%id == diag_cs%axesCv1%id) then
        diag%mask2d =>  diag_cs%mask2dCv
  !   else
  !       call SIS_error(FATAL, "SIS_diag_mediator:register_diag_field: " // &
  !            "unknown axes for diagnostic variable "//trim(field_name))
      endif
    else
      call MOM_error(FATAL, "MOM_IS_diag_mediator:register_diag_field: " // &
           "unknown axes for diagnostic variable "//trim(field_name))
    endif
  endif ! if (primary_id>-1)

  register_diag_field = primary_id

end function register_MOM_IS_diag_field

!> Registers a static diagnostic, returning an integer handle
function register_static_field(module_name, field_name, axes, &
     long_name, units, missing_value, range, mask_variant, standard_name, &
     do_not_log, interp_method, tile_count)
  integer :: register_static_field !< The returned diagnostic handle
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ice_model"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axesType),   intent(in) :: axes       !< The axis group for this field
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
  integer,          optional, intent(in) :: tile_count   !< no clue (not used in MOM_IS?)

  ! Local variables
  character(len=240) :: mesg
  real :: MOM_missing_value
  integer :: primary_id, fms_id
  type(diag_ctrl), pointer :: diag_cs !< A structure that is used to regulate diagnostic output

  MOM_missing_value = axes%diag_cs%missing_value
  if(present(missing_value)) MOM_missing_value = missing_value

  diag_cs => axes%diag_cs
  primary_id = -1

  fms_id = register_static_field_fms(module_name, field_name, axes%handles, &
       long_name=long_name, units=units, missing_value=MOM_missing_value, &
       range=range, mask_variant=mask_variant, standard_name=standard_name, &
       do_not_log=do_not_log, &
       interp_method=interp_method, tile_count=tile_count)
  if (fms_id > 0) then
    primary_id = get_new_diag_id(diag_cs)
    diag_cs%diags(primary_id)%fms_diag_id = fms_id
  endif

  register_static_field = primary_id

end function register_static_field

!> Add a description of an option to the documentation file
subroutine describe_option(opt_name, value, diag_CS)
  character(len=*),    intent(in) :: opt_name !< The name of the option
  character(len=*),    intent(in) :: value    !< The value of the option
  type(diag_ctrl), intent(in) :: diag_CS  !< Diagnostic being documented

  ! Local variables
  character(len=240) :: mesg
  integer :: start_ind = 1, end_ind, len_ind

  len_ind = len_trim(value)

  mesg = "    ! "//trim(opt_name)//": "//trim(value)
  write(diag_CS%doc_unit, '(a)') trim(mesg)
end subroutine describe_option

!> Convert the first n elements (up to 3) of an integer array to an underscore delimited string.
function i2s(a, n_in)
  integer, dimension(:), intent(in) :: a    !< The array of integers to translate
  integer, optional    , intent(in) :: n_in !< The number of elements to translate, by default all
  character(len=15) :: i2s !< The returned string

  ! Local variables
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

!> Initialize the MOM_IS diag_mediator and opens the available diagnostics file.
subroutine diag_mediator_init(G, param_file, diag_cs, component, err_msg, &
                                  doc_file_dir)
  type(ocean_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(param_file_type),      intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl),        intent(inout) :: diag_cs !< A structure that is used to regulate diagnostic output
  character(len=*), optional, intent(in)    :: component !< An opitonal component name
  character(len=*), optional, intent(out)   :: err_msg !< A string for a returned error message
  character(len=*), optional, intent(in)    :: doc_file_dir !< A directory in which to create the file

  ! This subroutine initializes the diag_mediator and the diag_manager.
  ! The grid type should have its dimensions set by this point, but it
  ! is not necessary that the metrics and axis labels be set up yet.

  ! Local variables
  integer :: ios, new_unit
  logical :: opened, new_file
  character(len=8)   :: this_pe
  character(len=240) :: doc_file, doc_file_dflt, doc_path
  character(len=40)  :: doc_file_param
  character(len=40)  :: mdl = "MOM_IS_diag_mediator" ! This module's name.

  call diag_manager_init(err_msg=err_msg)

  ! Allocate list of all diagnostics
  allocate(diag_cs%diags(DIAG_ALLOC_CHUNK_SIZE))
  diag_cs%next_free_diag_id = 1
  diag_cs%diags(:)%in_use = .false.

  diag_cs%is = G%isc - (G%isd-1) ; diag_cs%ie = G%iec - (G%isd-1)
  diag_cs%js = G%jsc - (G%jsd-1) ; diag_cs%je = G%jec - (G%jsd-1)
  diag_cs%isd = G%isd ; diag_cs%ied = G%ied ; diag_cs%jsd = G%jsd ; diag_cs%jed = G%jed

  if (is_root_pe() .and. (diag_CS%doc_unit < 0)) then
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
                 default=doc_file_dflt)
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

  call diag_masks_set(G, -1.0e34, diag_cs)

end subroutine diag_mediator_init

subroutine diag_masks_set(G, missing_value, diag_cs)
! Setup the 2d masks for diagnostics
  type(ocean_grid_type), target, intent(in)    :: G   !< The horizontal grid type
  real,                            intent(in)    :: missing_value !< A fill value for missing points
  type(diag_ctrl),             intent(inout) :: diag_cs !< A structure that is used to regulate diagnostic output

  ! Local variables
  integer :: i, j, k, NkIce, CatIce


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
subroutine diag_mediator_close_registration(diag_CS)
  type(diag_ctrl), intent(inout) :: diag_CS !< A structure that is used to regulate diagnostic output

  if (diag_CS%doc_unit > -1) then
    close(diag_CS%doc_unit) ; diag_CS%doc_unit = -2
  endif

end subroutine diag_mediator_close_registration

!> Deallocate memory associated with the MOM_IS diag mediator
subroutine diag_mediator_end(time, diag_CS)
  type(time_type), intent(in) :: time !< The current model time
  type(diag_ctrl), intent(inout) :: diag_CS !< A structure that is used to regulate diagnostic output

  if (diag_CS%doc_unit > -1) then
    close(diag_CS%doc_unit) ; diag_CS%doc_unit = -3
  endif

end subroutine diag_mediator_end

!> Allocate a new diagnostic id, noting that it may be necessary to expand the diagnostics array.
function get_new_diag_id(diag_cs)

  integer :: get_new_diag_id !< The returned ID for the new diagnostic
  type(diag_ctrl), intent(inout) :: diag_cs !< A structure that is used to regulate diagnostic output

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

    ! Initialise new part of the diag array.
    do i=diag_cs%next_free_diag_id, size(diag_cs%diags)
      diag_cs%diags(i)%in_use = .false.
    enddo
  endif

  get_new_diag_id = diag_cs%next_free_diag_id
  diag_cs%next_free_diag_id = diag_cs%next_free_diag_id + 1

end function get_new_diag_id

end module MOM_IS_diag_mediator
