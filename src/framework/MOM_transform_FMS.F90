!> Support functions and interfaces to permit transformed model domains to
!! interact with FMS operations registered on the non-transformed domains.

module MOM_transform_FMS

use MOM_array_transform, only : allocate_rotated_array, rotate_array
use MOM_error_handler, only : MOM_error, FATAL
use horiz_interp_mod, only : horiz_interp_type
use time_manager_mod, only : time_type
use time_interp_external_mod, only : time_interp_external

implicit none ; private

public rotated_time_interp_external

!> Read a field based on model time, and rotate to the model domain
interface rotated_time_interp_external
  module procedure rotated_time_interp_external_0d
  module procedure rotated_time_interp_external_2d
  module procedure rotated_time_interp_external_3d
end interface rotated_time_interp_external

contains

! NOTE: No transformations are applied to the 0d and 1d field implementations,
!   but are provided to maintain compatibility with the FMS interfaces.

!> Read a scalar field based on model time
!! This function is provided to support the full FMS time_interp_external
!! interface.
subroutine rotated_time_interp_external_0d(fms_id, time, data_in, verbose, &
    turns)
  integer, intent(in) :: fms_id                   !< FMS field ID
  type(time_type), intent(in) :: time             !< Model time
  real, intent(inout) :: data_in  !< field to write data
  logical, intent(in), optional :: verbose        !< Verbose output
  integer, intent(in), optional :: turns          !< Number of quarter turns

  if (present(turns)) &
    call MOM_error(FATAL, "Rotation not supported for 0d fields.")

  call time_interp_external(fms_id, time, data_in, verbose=verbose)
end subroutine rotated_time_interp_external_0d

!> Read a 2d field based on model time, and rotate to the model grid
subroutine rotated_time_interp_external_2d(fms_id, time, data_in, interp, &
    verbose, horz_interp, mask_out, is_in, ie_in, js_in, je_in, window_id, &
    turns)
  integer, intent(in) :: fms_id
  type(time_type), intent(in) :: time
  real, dimension(:,:), intent(inout) :: data_in
  integer, intent(in), optional :: interp
  logical, intent(in), optional :: verbose
  type(horiz_interp_type),intent(in), optional :: horz_interp
  logical, dimension(:,:), intent(out), optional :: mask_out
  integer, intent(in), optional :: is_in, ie_in, js_in, je_in
  integer, intent(in), optional :: window_id
  integer, intent(in), optional :: turns

  real, allocatable :: data_pre(:,:)
  integer :: qturns

  ! TODO: Mask rotation requires logical array rotation support
  if (present(mask_out)) &
    call MOM_error(FATAL, "Rotation of masked output not yet support")

  qturns = 0
  if (present(turns)) &
    qturns = modulo(turns, 4)


  if (qturns == 0) then
    call time_interp_external(fms_id, time, data_in, interp=interp, &
        verbose=verbose, horz_interp=horz_interp, mask_out=mask_out, &
        is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in, &
        window_id=window_id)
  else
    call allocate_rotated_array(data_in, [1,1], -qturns, data_pre)
    call time_interp_external(fms_id, time, data_pre, interp=interp, &
        verbose=verbose, horz_interp=horz_interp, mask_out=mask_out, &
        is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in, &
        window_id=window_id)
    call rotate_array(data_pre, turns, data_in)
    deallocate(data_pre)
  endif
end subroutine rotated_time_interp_external_2d


!> Read a 3d field based on model time, and rotate to the model grid
subroutine rotated_time_interp_external_3d(fms_id, time, data_in, interp, &
    verbose, horz_interp, mask_out, is_in, ie_in, js_in, je_in, window_id, &
    turns)
  integer, intent(in) :: fms_id
  type(time_type), intent(in) :: time
  real, dimension(:,:,:), intent(inout) :: data_in
  integer, intent(in), optional :: interp
  logical, intent(in), optional :: verbose
  type(horiz_interp_type),intent(in), optional :: horz_interp
  logical, dimension(:,:,:), intent(out), optional :: mask_out
  integer, intent(in), optional :: is_in, ie_in, js_in, je_in
  integer, intent(in), optional :: window_id
  integer, intent(in), optional :: turns

  real, allocatable :: data_pre(:,:,:)
  integer :: qturns

  ! TODO: Mask rotation requires logical array rotation support
  if (present(mask_out)) &
    call MOM_error(FATAL, "Rotation of masked output not yet support")

  qturns = 0
  if (present(turns)) &
    qturns = modulo(turns, 4)

  if (qturns == 0) then
    call time_interp_external(fms_id, time, data_in, interp=interp, &
        verbose=verbose, horz_interp=horz_interp, mask_out=mask_out, &
        is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in, &
        window_id=window_id)
  else
    call allocate_rotated_array(data_in, [1,1,1], -qturns, data_pre)
    call time_interp_external(fms_id, time, data_pre, interp=interp, &
        verbose=verbose, horz_interp=horz_interp, mask_out=mask_out, &
        is_in=is_in, ie_in=ie_in, js_in=js_in, je_in=je_in, &
        window_id=window_id)
    call rotate_array(data_pre, turns, data_in)
    deallocate(data_pre)
  endif
end subroutine rotated_time_interp_external_3d

end module MOM_transform_FMS
