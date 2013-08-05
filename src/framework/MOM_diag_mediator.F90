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

use MOM_coms, only : PE_here
use MOM_error_handler, only : MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : vardesc
use MOM_safe_alloc, only : safe_alloc_ptr, safe_alloc_alloc
use MOM_string_functions, only : lowercase
use MOM_time_manager, only : time_type

use diag_manager_mod, only : diag_manager_init, diag_manager_end
use diag_manager_mod, only : send_data, diag_axis_init
use diag_manager_mod, only : register_diag_field_fms=>register_diag_field
use diag_manager_mod, only : register_static_field

implicit none ; private

public set_axes_info, post_data, register_diag_field, time_type
public safe_alloc_ptr, safe_alloc_alloc
public enable_averaging, disable_averaging, query_averaging_enabled
public diag_mediator_init, diag_mediator_end, set_diag_mediator_grid
public diag_mediator_close_registration, get_diag_time_end
public diag_axis_init, ocean_register_diag, register_static_field

interface post_data
  module procedure post_data_3d, post_data_2d
end interface post_data

!   The following data type contains pointers to diagnostic fields that might
! be shared between modules, and also to the variables that control the handling
! of model output.
type, public :: diag_ctrl

! The following fields are used for the output of the data.
  integer :: is, ie, js, je
  integer :: isd, ied, jsd, jed
  real :: time_int              ! The time interval in s for any fields
                                ! that are offered for averaging.
  type(time_type) :: time_end   ! The end time of the valid
                                ! interval for any offered field.
  logical :: ave_enabled = .false. ! .true. if averaging is enabled.

  ! The following are axis types defined for output.
  integer, dimension(3) :: axesBL, axesTL, axesCuL, axesCvL
  integer, dimension(3) :: axesBi, axesTi, axesCui, axesCvi
  integer, dimension(2) :: axesB1, axesT1, axesCu1, axesCv1
  integer, dimension(1) :: axeszi, axeszL

end type diag_ctrl

integer :: doc_unit = -1

contains

subroutine set_axes_info(G, param_file, diag, set_vertical)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
  type(diag_ctrl),       intent(inout) :: diag
  logical, optional,     intent(in)    :: set_vertical
! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (inout)   diag - A structure that is used to regulate diagnostic output.
!  (in,opt)  set_vertical - If true (or missing), set up the vertical axes.
  integer :: id_xq, id_yq, id_zl, id_zi, id_xh, id_yh, k, nz
  real :: zlev(G%ks:G%ke), zinter(G%ks:G%ke+1)
  logical :: set_vert, Cartesian_grid
  character(len=80) :: grid_config, units_temp
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod  = "MOM_diag_mediator" ! This module's name.
  nz = G%ke

  set_vert = .true. ; if (present(set_vertical)) set_vert = set_vertical

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version)
  call get_param(param_file, mod, "GRID_CONFIG", grid_config, &
                 "The method for defining the horizontal grid.  Valid \n"//&
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
    call get_param(param_file, mod, "AXIS_UNITS", units_temp, &
                 "The units for the x- and y- axis labels.  AXIS_UNITS \n"//&
                 "should be defined as 'k' for km, 'm' for m, or 'd' \n"//&
                 "for degrees of latitude and longitude (the default). \n"//&
                 "Except on a Cartesian grid, only degrees are currently \n"//&
                 "implemented.", default='degrees')
    if (units_temp(1:1) == 'k') then
      G%x_axis_units = "kilometers" ; G%y_axis_units = "kilometers"
    elseif (units_temp(1:1) == 'm') then
      G%x_axis_units = "meters" ; G%y_axis_units = "meters"
    endif
    call log_param(param_file, mod, "explicit AXIS_UNITS", G%x_axis_units)
  else
    Cartesian_grid = .false.
  endif
  
  do k=1,nz ; zlev(k) = G%Rlay(k) ; enddo
  zinter(1) = 1.5*G%Rlay(1) - 0.5*G%Rlay(2)
  do k=2,nz ; zinter(k) = 0.5*(G%Rlay(k) + G%Rlay(k-1)) ; enddo
  zinter(nz+1) = 1.5*G%Rlay(nz) - 0.5*G%Rlay(nz-1)

!  do i=1,nz ; zlev(i) = real(i) ; enddo
!  do i=1,nz+1 ; zinter(i) = real(i) - 0.5 ; enddo
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
    id_zl = diag_axis_init('zl', zlev, 'layer', 'z', 'cell depth')
    id_zi = diag_axis_init('zi', zinter, 'interface', 'z', 'cell interface depth')
  else
    id_zl = -1 ; id_zi = -1
  endif

  ! Vertical axes for the interfaces and layers.
  diag%axeszi(1) = id_zi ; diag%axeszL(1) = id_zL

  ! Axis groupings for the model layers.
  diag%axesTL(:) = (/ id_xh, id_yh, id_zL /)
  diag%axesBL(:) = (/ id_xq, id_yq, id_zL /)
  diag%axesCuL(:) = (/ id_xq, id_yh, id_zL /)
  diag%axesCvL(:) = (/ id_xh, id_yq, id_zL /)

  ! Axis groupings for the model interfaces.
  diag%axesTi(:) = (/ id_xh, id_yh, id_zi /)
  diag%axesCui(:) = (/ id_xq, id_yh, id_zi /)
  diag%axesCvi(:) = (/ id_xh, id_yq, id_zi /)
  diag%axesBi(:) = (/ id_xq, id_yq, id_zi /)

  ! Axis groupings for 2-D arrays.
  diag%axesT1(:) = (/ id_xh, id_yh /)
  diag%axesB1(:) = (/ id_xq, id_yq /)
  diag%axesCu1(:) = (/ id_xq, id_yh /)
  diag%axesCv1(:) = (/ id_xh, id_yq /)
 
end subroutine set_axes_info

subroutine set_diag_mediator_grid(G, diag)
  type(ocean_grid_type), intent(inout) :: G
  type(diag_ctrl),       intent(inout) :: diag
! Arguments: G - The ocean's grid structure.
!  (inout)   diag - A structure that is used to regulate diagnostic output.
  diag%is = G%isc ; diag%ie = G%iec ; diag%js = G%jsc ; diag%je = G%jec
  diag%isd = G%isd ; diag%ied = G%ied ; diag%jsd = G%jsd ; diag%jed = G%jed
end subroutine set_diag_mediator_grid

subroutine post_data_2d(diag_field_id, field, diag, is_static, mask)
  integer,           intent(in) :: diag_field_id
  real,              intent(in) :: field(:,:)
  type(diag_ctrl),   intent(in) :: diag
  logical, optional, intent(in) :: is_static
  real,    optional, intent(in) :: mask(:,:)
! Arguments: diag_field_id - the id for an output variable returned by a
!                            previous call to register_diag_field.
!  (in)      field - The 2-d array being offered for output or averaging.
!  (inout)   diag - A structure that is used to regulate diagnostic output.
!  (in,opt)  is_static - If true, this is a static field that is always offered.
!  (in,opt)  mask - If present, use this real array as the data mask.
  logical :: used, is_stat
  integer :: ishift, jshift

  is_stat = .false. ; if (present(is_static)) is_stat = is_static 

  ishift = 0 ; jshift = 0
  if ( size(field,1) == diag%ied-diag%isd +1 ) then
    ishift = 0
  elseif ( size(field,1) == diag%ied-diag%isd +2 ) then
    ishift = 1
!  else
!    call MOM_error(FATAL,"post_data_2d: peculiar size in i-direction")
  endif
  if ( size(field,2) == diag%jed-diag%jsd +1 ) then
    jshift = 0
  elseif ( size(field,2) == diag%jed-diag%jsd +2 ) then
    jshift = 1
!  else
!    call MOM_error(FATAL,"post_data_2d: peculiar size in j-direction")
  endif

  if (present(mask)) then
    if ((size(field,1) /= size(mask,1)) .or. &
        (size(field,2) /= size(mask,2))) then
      call MOM_error(FATAL, "post_data_2d: post_data called with a mask "//&
                             "that does not match the size of field.")
    endif
  endif

  if (is_stat) then
    if (present(mask)) then
      used = send_data(diag_field_id, field, &
                       is_in = diag%is-ishift, js_in = diag%js-jshift, &
                       ie_in = diag%ie, je_in = diag%je, rmask=mask)
    else
      used = send_data(diag_field_id, field, &
                       is_in = diag%is-ishift, js_in = diag%js-jshift, &
                       ie_in = diag%ie, je_in = diag%je)
    endif
  elseif (diag%ave_enabled) then
    if (present(mask)) then
      used = send_data(diag_field_id, field, diag%time_end, &
                       is_in = diag%is-ishift, js_in = diag%js-jshift, &
                       ie_in = diag%ie, je_in = diag%je, &
                       weight=diag%time_int, rmask=mask)
    else
      used = send_data(diag_field_id, field, diag%time_end, &
                       is_in = diag%is-ishift, js_in = diag%js-jshift, &
                       ie_in = diag%ie, je_in = diag%je, &
                       weight=diag%time_int)
    endif
  endif

end subroutine post_data_2d

subroutine post_data_3d(diag_field_id, field, diag, is_static, mask)
  integer,           intent(in) :: diag_field_id
  real,              intent(in) :: field(:,:,:)
  type(diag_ctrl),   intent(in) :: diag
  logical, optional, intent(in) :: is_static
  real,    optional, intent(in) :: mask(:,:,:)
! Arguments: diag_field_id - the id for an output variable returned by a
!                            previous call to register_diag_field.
!  (in)      field - The 3-d array being offered for output or averaging.
!  (inout)   diag - A structure that is used to regulate diagnostic output.
!  (in)      static - If true, this is a static field that is always offered.
!  (in,opt)  mask - If present, use this real array as the data mask.
  logical :: used  ! The return value of send_data is not used for anything.
  logical :: is_stat
  integer :: ishift, jshift
  is_stat = .false. ; if (present(is_static)) is_stat = is_static 
  
  ishift = 0 ; jshift = 0
  if ( size(field,1) == diag%ied-diag%isd +1 ) then
    ishift = 0
  elseif ( size(field,1) == diag%ied-diag%isd +2 ) then
    ishift = 1
!  else
!    call MOM_error(FATAL,"post_data_3d: peculiar size in i-direction")
  endif
  if ( size(field,2) == diag%jed-diag%jsd +1 ) then
    jshift = 0
  elseif ( size(field,2) == diag%jed-diag%jsd +2 ) then
    jshift = 1  
!  else
!    call MOM_error(FATAL,"post_data_3d: peculiar size in j-direction")
  endif

  if (present(mask)) then
    if ((size(field,1) /= size(mask,1)) .or. &
        (size(field,2) /= size(mask,2)) .or. &
        (size(field,3) /= size(mask,3))) then
      call MOM_error(FATAL, "post_data_3d: post_data called with a mask "//&
                             "that does not match the size of field.")
    endif
  endif

  if (is_stat) then
    if (present(mask)) then
      used = send_data(diag_field_id, field, &
                       is_in = diag%is-ishift, js_in = diag%js-jshift, &
                       ie_in = diag%ie, je_in = diag%je, rmask=mask)
    else
      used = send_data(diag_field_id, field, &
                       is_in = diag%is-ishift, js_in = diag%js-jshift, &
                       ie_in = diag%ie, je_in = diag%je)
    endif
  elseif (diag%ave_enabled) then
    if (present(mask)) then
      used = send_data(diag_field_id, field, diag%time_end, &
                       is_in = diag%is-ishift, js_in = diag%js-jshift, &
                       ie_in = diag%ie, je_in = diag%je, &
                       weight=diag%time_int, rmask=mask)
    else
      used = send_data(diag_field_id, field, diag%time_end, &
                       is_in = diag%is-ishift, js_in = diag%js-jshift, &
                       ie_in = diag%ie, je_in = diag%je, &
                       weight=diag%time_int)
    endif
  endif

end subroutine post_data_3d


subroutine enable_averaging(time_int_in, time_end_in, diag)
  real, intent(in) :: time_int_in
  type(time_type), intent(in) :: time_end_in
  type(diag_ctrl), intent(inout) :: diag
! This subroutine enables the accumulation of time averages over the
! specified time interval.

! Arguments: time_int_in - the time interval in s over which any
!                          values that are offered are valid.
!  (in)      time_end_in - the end time in s of the valid interval.
!  (inout)   diag - A structure that is used to regulate diagnostic output.

!  if (num_file==0) return
  diag%time_int = time_int_in
  diag%time_end = time_end_in
  diag%ave_enabled = .true.
end subroutine enable_averaging

! Call this subroutine to avoid averaging any offered fields.
subroutine disable_averaging(diag)
  type(diag_ctrl), intent(inout) :: diag
! Argument: diag - A structure that is used to regulate diagnostic output.

  diag%time_int = 0.0
  diag%ave_enabled = .false.

end subroutine disable_averaging

! Call this subroutine to determine whether the averaging is
! currently enabled.  .true. is returned if it is.
function query_averaging_enabled(diag, time_int, time_end)
  type(diag_ctrl),           intent(in)  :: diag
  real,            optional, intent(out) :: time_int
  type(time_type), optional, intent(out) :: time_end
  logical :: query_averaging_enabled
! Arguments: diag - A structure that is used to regulate diagnostic output.
!  (out,opt) time_int - The current setting of diag%time_int, in s.
!  (out,opt) time_end - The current setting of diag%time_end.

  if (present(time_int)) time_int = diag%time_int
  if (present(time_end)) time_end = diag%time_end
  query_averaging_enabled = diag%ave_enabled
end function query_averaging_enabled

function get_diag_time_end(diag)
  type(diag_ctrl),           intent(in)  :: diag
  type(time_type) :: get_diag_time_end
! Argument: diag - A structure that is used to regulate diagnostic output.

!   This function returns the valid end time for diagnostics that are handled
! outside of the MOM6 infrastructure, such as via the generic tracer code.

  get_diag_time_end = diag%time_end
end function get_diag_time_end

function register_diag_field(module_name, field_name, axes, init_time, &
             long_name, units, missing_value, range, mask_variant, standard_name, &
             verbose, do_not_log, err_msg, interp_method, tile_count)
  integer :: register_diag_field
  character(len=*), intent(in) :: module_name, field_name
  integer,          intent(in) :: axes(:)
  type(time_type),  intent(in) :: init_time
  character(len=*), optional, intent(in) :: long_name, units, standard_name
  real,             optional, intent(in) :: missing_value, range(2)
  logical,          optional, intent(in) :: mask_variant, verbose, do_not_log
  character(len=*), optional, intent(out):: err_msg
  character(len=*), optional, intent(in) :: interp_method
  integer,          optional, intent(in) :: tile_count
! Output:    An integer handle for a diagnostic array.
! Arguments: module_name - The name of this module, usually "ocean_model" or "ice_shelf_model".
!  (in)      field_name - The name of the diagnostic field.
!  (in)      axes - A set of up to 3 integers that indicates the axes for this field.
!  (in)      init_time - The time at which a field is first available?
!  (in,opt)  long_name - The long name of a field.
!  (in,opt)  units - The units of a field.
!  (in,opt)  standard_name - The standardized name associated with a field. (Not yet used in MOM.)
!  (in,opt)  missing_value - A value that indicates missing values.
!  (in,opt)  range - The valid range of a variable. (Not used in MOM.)
!  (in,opt)  mask_variant - If true a logical mask must be provided with post_data calls.  (Not used in MOM.)
!  (in,opt)  verbose - If true, FMS is verbosed. (Not used in MOM.)
!  (in,opt)  do_not_log - If true, do not log something. (Not used in MOM.)
!  (out,opt) err_msg - An character string into which an error message might be placed. (Not used in MOM.)
!  (in,opt)  interp_method - No clue. (Not used in MOM.)
!  (in,opt)  tile_count - No clue. (Not used in MOM.)
  character(len=240) :: mesg

  register_diag_field = register_diag_field_fms(module_name, field_name, axes, &
         init_time, long_name=long_name, units=units, missing_value=missing_value, &
         range=range, mask_variant=mask_variant, standard_name=standard_name, &
         verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
         interp_method=interp_method, tile_count=tile_count)

  if (is_root_pe() .and. doc_unit > 0) then
    if (register_diag_field > 0) then
      mesg = '"'//trim(module_name)//'", "'//trim(field_name)//'"  [Used]'
    else
      mesg = '"'//trim(module_name)//'", "'//trim(field_name)//'"  [Unused]'
    endif
    write(doc_unit, '(a)') trim(mesg)
    if (present(long_name)) call describe_option("long_name", long_name)
    if (present(units)) call describe_option("units", units)
    if (present(standard_name)) call describe_option("standard_name", standard_name)
  endif

end function register_diag_field

subroutine describe_option(opt_name, value)
  character(len=*), intent(in) :: opt_name, value

  character(len=240) :: mesg
  integer :: len_ind

  len_ind = len_trim(value)  ! Add error handling for long values?

  mesg = "    ! "//trim(opt_name)//": "//trim(value)
  write(doc_unit, '(a)') trim(mesg)
end subroutine describe_option

function ocean_register_diag(var_desc, G, diag, day)
  integer :: ocean_register_diag
  type(vardesc),         intent(in) :: var_desc
  type(ocean_grid_type), intent(in) :: G
  type(diag_ctrl),       intent(in) :: diag
  type(time_type),       intent(in) :: day

  integer, dimension(:), allocatable :: axes

  ! Use the hor_grid and z_grid components of vardesc to determine the 
  ! desired axes to register the diagnostic field for.
  select case (var_desc%z_grid)

    case ("L")
      select case (var_desc%hor_grid)
        case ("q")
          allocate(axes(3)) ; axes(:) = diag%axesBL(:)
        case ("h")
          allocate(axes(3)) ; axes(:) = diag%axesTL(:)
        case ("u")
          allocate(axes(3)) ; axes(:) = diag%axesCuL(:)
        case ("v")
          allocate(axes(3)) ; axes(:) = diag%axesCvL(:)
        case ("Bu")
          allocate(axes(3)) ; axes(:) = diag%axesBL(:)
        case ("T")
          allocate(axes(3)) ; axes(:) = diag%axesTL(:)
        case ("Cu")
          allocate(axes(3)) ; axes(:) = diag%axesCuL(:)
        case ("Cv")
          allocate(axes(3)) ; axes(:) = diag%axesCvL(:)
        case ("z")
          allocate(axes(1)) ; axes(:) = diag%axeszL(:)
        case default
          call MOM_error(FATAL, "ocean_register_diag: " // &
              "unknown hor_grid component "//trim(var_desc%hor_grid))
      end select

    case ("i")
      select case (var_desc%hor_grid)
        case ("q")
          allocate(axes(3)) ; axes(:) = diag%axesBi(:)
        case ("h")
          allocate(axes(3)) ; axes(:) = diag%axesTi(:)
        case ("u")
          allocate(axes(3)) ; axes(:) = diag%axesCui(:)
        case ("v")
          allocate(axes(3)) ; axes(:) = diag%axesCvi(:)
        case ("Bu")
          allocate(axes(3)) ; axes(:) = diag%axesBi(:)
        case ("T")
          allocate(axes(3)) ; axes(:) = diag%axesTi(:)
        case ("Cu")
          allocate(axes(3)) ; axes(:) = diag%axesCui(:)
        case ("Cv")
          allocate(axes(3)) ; axes(:) = diag%axesCvi(:)
        case ("z")
          allocate(axes(1)) ; axes(:) = diag%axeszi(:)
        case default
          call MOM_error(FATAL, "ocean_register_diag: " // &
            "unknown hor_grid component "//trim(var_desc%hor_grid))
      end select

    case ("1")
      allocate(axes(2))
      select case (var_desc%hor_grid)
        case ("q")
          axes(:) = diag%axesB1(:)
        case ("h")
          axes(:) = diag%axesT1(:)
        case ("u")
          axes(:) = diag%axesCu1(:)
        case ("v")
          axes(:) = diag%axesCv1(:)
        case ("Bu")
          axes(:) = diag%axesB1(:)
        case ("T")
          axes(:) = diag%axesT1(:)
        case ("Cu")
          axes(:) = diag%axesCu1(:)
        case ("Cv")
          axes(:) = diag%axesCv1(:)
        case default
          call MOM_error(FATAL, "ocean_register_diag: " // &
            "unknown hor_grid component "//trim(var_desc%hor_grid))
      end select

    case default
      call MOM_error(FATAL,&
        "ocean_register_diag: unknown z_grid component "//trim(var_desc%z_grid))
  end select

  ocean_register_diag = register_diag_field("ocean_model", trim(var_desc%name), axes, &
        day, trim(var_desc%longname), trim(var_desc%units), missing_value = -1.0e+34)

  if (allocated(axes)) deallocate(axes)

end function ocean_register_diag

subroutine diag_mediator_init(G, param_file, diag, err_msg)
  type(ocean_grid_type),      intent(inout) :: G
  type(param_file_type),      intent(in)    :: param_file
  type(diag_ctrl),            intent(inout) :: diag
  character(len=*), optional, intent(out)   :: err_msg

  ! This subroutine initializes the diag_mediator and the diag_manager.
  ! The grid type should have its dimensions set by this point, but it
  ! is not necessary that the metrics and axis labels be set up yet.
  integer :: ios
  logical :: opened, new_file
  character(len=8)   :: this_pe
  character(len=240) :: doc_file, doc_file_dflt
  character(len=40)  :: mod  = "MOM_diag_mediator" ! This module's name.

  call diag_manager_init(err_msg=err_msg)

  diag%is = G%isc ; diag%ie = G%iec ; diag%js = G%jsc ; diag%je = G%jec
  diag%isd = G%isd ; diag%ied = G%ied ; diag%jsd = G%jsd ; diag%jed = G%jed

  if (is_root_pe()) then
    write(this_pe,'(i6.6)') PE_here()
    doc_file_dflt = "available_diags."//this_pe
    call get_param(param_file, mod, "AVAILABLE_DIAGS_FILE", doc_file, &
                 "A file into which to write a list of all available \n"//&
                 "ocean diagnostics that can be included in a diag_table.", &
                 default=doc_file_dflt)
    if (len_trim(doc_file) > 0) then
      new_file = .true. ; if (doc_unit /= -1) new_file = .false.
    ! Find an unused unit number.
      do doc_unit=512,42,-1
        inquire( doc_unit, opened=opened)
        if (.not.opened) exit
      enddo

      if (opened) call MOM_error(FATAL, &
          "diag_mediator_init failed to find an unused unit number.")

      if (new_file) then
        open(doc_unit, file=trim(doc_file), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='REPLACE', iostat=ios)
      else ! This file is being reopened, and should be appended.
        open(doc_unit, file=trim(doc_file), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='OLD', position='APPEND', iostat=ios)
      endif
      inquire(doc_unit, opened=opened)
      if ((.not.opened) .or. (ios /= 0)) then
        call MOM_error(FATAL, "Failed to open available diags file "//trim(doc_file)//".")
      endif
    endif
  endif

end subroutine diag_mediator_init

subroutine diag_mediator_close_registration( )

  if (doc_unit > -1) then
    close(doc_unit) ; doc_unit = -2
  endif

end subroutine diag_mediator_close_registration

subroutine diag_mediator_end(time)
  type(time_type), intent(in) :: time

  call diag_manager_end(time)

  if (doc_unit > -1) then
    close(doc_unit) ; doc_unit = -3
  endif

end subroutine diag_mediator_end

end module MOM_diag_mediator
