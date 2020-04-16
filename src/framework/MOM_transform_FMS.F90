!> Support functions and interfaces to permit transformed model domains to
!! interact with FMS operations registered on the non-transformed domains.

module MOM_transform_FMS

use horiz_interp_mod, only : horiz_interp_type
use MOM_error_handler, only : MOM_error, FATAL
use MOM_io, only : fieldtype, write_field
use mpp_domains_mod, only : domain2D
use fms_mod, only : mpp_chksum
use time_manager_mod, only : time_type
use time_interp_external_mod, only : time_interp_external

use MOM_array_transform, only : allocate_rotated_array, rotate_array

implicit none

private
public rotated_mpp_chksum
public rotated_write_field
public rotated_time_interp_external

!> Rotate and compute the FMS (mpp) checksum of a field
interface rotated_mpp_chksum
  module procedure rotated_mpp_chksum_real_0d
  module procedure rotated_mpp_chksum_real_1d
  module procedure rotated_mpp_chksum_real_2d
  module procedure rotated_mpp_chksum_real_3d
  module procedure rotated_mpp_chksum_real_4d
end interface rotated_mpp_chksum

!> Rotate and write a registered field to an FMS output file
interface rotated_write_field
  module procedure rotated_write_field_real_0d
  module procedure rotated_write_field_real_1d
  module procedure rotated_write_field_real_2d
  module procedure rotated_write_field_real_3d
  module procedure rotated_write_field_real_4d
end interface rotated_write_field

!> Read a field based on model time, and rotate to the model domain
interface rotated_time_interp_external
  module procedure rotated_time_interp_external_0d
  module procedure rotated_time_interp_external_2d
  module procedure rotated_time_interp_external_3d
end interface rotated_time_interp_external

contains

! NOTE: No transformations are applied to the 0d and 1d field implementations,
!   but are provided to maintain compatibility with the FMS interfaces.


!> Compute the FMS (mpp) checksum of a scalar.
!! This function is provided to support the full FMS mpp_chksum interface.
function rotated_mpp_chksum_real_0d(field, pelist, mask_val, turns) &
    result(chksum)
  real, intent(in) :: field                   !> Input scalar
  integer, optional, intent(in) :: pelist(:)  !> PE list of ranks to checksum
  real, optional, intent(in) :: mask_val      !> FMS mask value
  integer, optional, intent(in) :: turns      !> Number of quarter turns
  integer :: chksum                           !> FMS checksum of scalar

  if (present(turns)) &
    call MOM_error(FATAL, "Rotation not supported for 0d fields.")

  chksum = mpp_chksum(field, pelist=pelist, mask_val=mask_val)
end function rotated_mpp_chksum_real_0d


!> Compute the FMS (mpp) checksum of a 1d field.
!! This function is provided to support the full FMS mpp_chksum interface.
function rotated_mpp_chksum_real_1d(field, pelist, mask_val, turns) &
    result(chksum)
  real, intent(in) :: field(:)                !> Input field
  integer, optional, intent(in) :: pelist(:)  !> PE list of ranks to checksum
  real, optional, intent(in) :: mask_val      !> FMS mask value
  integer, optional, intent(in) :: turns      !> Number of quarter-turns
  integer :: chksum                           !> FMS checksum of field

  if (present(turns)) &
    call MOM_error(FATAL, "Rotation not supported for 1d fields.")

  chksum = mpp_chksum(field, pelist=pelist, mask_val=mask_val)
end function rotated_mpp_chksum_real_1d


!> Compute the FMS (mpp) checksum of a rotated 2d field.
function rotated_mpp_chksum_real_2d(field, pelist, mask_val, turns) &
    result(chksum)
  real, intent(in) :: field(:,:)              !> Unrotated input field
  integer, optional, intent(in) :: pelist(:)  !> PE list of ranks to checksum
  real, optional, intent(in) :: mask_val      !> FMS mask value
  integer, optional, intent(in) :: turns      !> Number of quarter-turns
  integer :: chksum                           !> FMS checksum of field

  real, allocatable :: field_rot(:,:)
  integer :: qturns

  qturns = 0
  if (present(turns)) &
    qturns = modulo(turns, 4)

  if (qturns == 0) then
    chksum = mpp_chksum(field, pelist=pelist, mask_val=mask_val)
  else
    call allocate_rotated_array(field, [1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    chksum = mpp_chksum(field_rot, pelist=pelist, mask_val=mask_val)
    deallocate(field_rot)
  endif
end function rotated_mpp_chksum_real_2d


!> Compute the FMS (mpp) checksum of a rotated 3d field.
function rotated_mpp_chksum_real_3d(field, pelist, mask_val, turns) &
    result(chksum)
  real, intent(in) :: field(:,:,:)            !> Unrotated input field
  integer, optional, intent(in) :: pelist(:)  !> PE list of ranks to checksum
  real, optional, intent(in) :: mask_val      !> FMS mask value
  integer, optional, intent(in) :: turns      !> Number of quarter-turns
  integer :: chksum                           !> FMS checksum of field

  real, allocatable :: field_rot(:,:,:)
  integer :: qturns

  qturns = 0
  if (present(turns)) &
    qturns = modulo(turns, 4)

  if (qturns == 0) then
    chksum = mpp_chksum(field, pelist=pelist, mask_val=mask_val)
  else
    call allocate_rotated_array(field, [1,1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    chksum = mpp_chksum(field_rot, pelist=pelist, mask_val=mask_val)
    deallocate(field_rot)
  endif
end function rotated_mpp_chksum_real_3d


!> Compute the FMS (mpp) checksum of a rotated 4d field.
function rotated_mpp_chksum_real_4d(field, pelist, mask_val, turns) &
    result(chksum)
  real, intent(in) :: field(:,:,:,:)          !> Unrotated input field
  integer, optional, intent(in) :: pelist(:)  !> PE list of ranks to checksum
  real, optional, intent(in) :: mask_val      !> FMS mask value
  integer, optional, intent(in) :: turns      !> Number of quarter-turns
  integer :: chksum                           !> FMS checksum of field

  real, allocatable :: field_rot(:,:,:,:)
  integer :: qturns

  qturns = 0
  if (present(turns)) &
    qturns = modulo(turns, 4)

  if (qturns == 0) then
    chksum = mpp_chksum(field, pelist=pelist, mask_val=mask_val)
  else
    call allocate_rotated_array(field, [1,1,1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    chksum = mpp_chksum(field_rot, pelist=pelist, mask_val=mask_val)
    deallocate(field_rot)
  endif
end function rotated_mpp_chksum_real_4d


! NOTE: In MOM_io, write_field points to mpp_write, which supports a very broad
! range of interfaces.  Here, we only support the much more narrow family of
! mpp_write_2ddecomp functions used to write tiled data.


!> Write the rotation of a 1d field to an FMS output file
!! This function is provided to support the full FMS write_field interface.
subroutine rotated_write_field_real_0d(io_unit, field_md, field, tstamp, turns)
  integer, intent(in) :: io_unit              !> File I/O unit handle
  type(fieldtype), intent(in) :: field_md     !> FMS field metadata
  real, intent(inout) :: field                !> Unrotated field array
  real, optional, intent(in) :: tstamp        !> Model timestamp
  integer, optional, intent(in) :: turns      !> Number of quarter-turns

  if (present(turns)) &
    call MOM_error(FATAL, "Rotation not supported for 0d fields.")

  call write_field(io_unit, field_md, field, tstamp=tstamp)
end subroutine rotated_write_field_real_0d


!> Write the rotation of a 1d field to an FMS output file
!! This function is provided to support the full FMS write_field interface.
subroutine rotated_write_field_real_1d(io_unit, field_md, field, tstamp, turns)
  integer, intent(in) :: io_unit              !> File I/O unit handle
  type(fieldtype), intent(in) :: field_md     !> FMS field metadata
  real, intent(inout) :: field(:)             !> Unrotated field array
  real, optional, intent(in) :: tstamp        !> Model timestamp
  integer, optional, intent(in) :: turns      !> Number of quarter-turns

  if (present(turns)) &
    call MOM_error(FATAL, "Rotation not supported for 0d fields.")

  call write_field(io_unit, field_md, field, tstamp=tstamp)
end subroutine rotated_write_field_real_1d


!> Write the rotation of a 2d field to an FMS output file
subroutine rotated_write_field_real_2d(io_unit, field_md, domain, field, &
    tstamp, tile_count, default_data, turns)
  integer, intent(in) :: io_unit              !> File I/O unit handle
  type(fieldtype), intent(in) :: field_md     !> FMS field metadata
  type(domain2D), intent(inout) :: domain     !> FMS MPP domain
  real, intent(inout) :: field(:,:)           !> Unrotated field array
  real, optional, intent(in) :: tstamp        !> Model timestamp
  integer, optional, intent(in) :: tile_count !> PEs per tile (default: 1)
  real, optional, intent(in) :: default_data  !> Default fill value
  integer, optional, intent(in) :: turns      !> Number of quarter-turns

  real, allocatable :: field_rot(:,:)
  integer :: qturns

  qturns = 0
  if (present(turns)) &
    qturns = modulo(turns, 4)

  if (qturns == 0) then
    call write_field(io_unit, field_md, domain, field, tstamp=tstamp, &
        tile_count=tile_count, default_data=default_data)
  else
    call allocate_rotated_array(field, [1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    call write_field(io_unit, field_md, domain, field_rot, tstamp=tstamp, &
        tile_count=tile_count, default_data=default_data)
    deallocate(field_rot)
  endif
end subroutine rotated_write_field_real_2d


!> Write the rotation of a 3d field to an FMS output file
subroutine rotated_write_field_real_3d(io_unit, field_md, domain, field, &
    tstamp, tile_count, default_data, turns)
  integer, intent(in) :: io_unit              !> File I/O unit handle
  type(fieldtype), intent(in) :: field_md     !> FMS field metadata
  type(domain2D), intent(inout) :: domain     !> FMS MPP domain
  real, intent(inout) :: field(:,:,:)         !> Unrotated field array
  real, optional, intent(in) :: tstamp        !> Model timestamp
  integer, optional, intent(in) :: tile_count !> PEs per tile (default: 1)
  real, optional, intent(in) :: default_data  !> Default fill value
  integer, optional, intent(in) :: turns      !> Number of quarter-turns

  real, allocatable :: field_rot(:,:,:)
  integer :: qturns

  qturns = 0
  if (present(turns)) &
    qturns = modulo(turns, 4)

  if (qturns == 0) then
    call write_field(io_unit, field_md, domain, field, tstamp=tstamp, &
        tile_count=tile_count, default_data=default_data)
  else
    call allocate_rotated_array(field, [1,1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    call write_field(io_unit, field_md, domain, field_rot, tstamp=tstamp, &
        tile_count=tile_count, default_data=default_data)
    deallocate(field_rot)
  endif
end subroutine rotated_write_field_real_3d


!> Write the rotation of a 4d field to an FMS output file
subroutine rotated_write_field_real_4d(io_unit, field_md, domain, field, &
    tstamp, tile_count, default_data, turns)
  integer, intent(in) :: io_unit              !> File I/O unit handle
  type(fieldtype), intent(in) :: field_md     !> FMS field metadata
  type(domain2D), intent(inout) :: domain     !> FMS MPP domain
  real, intent(inout) :: field(:,:,:,:)       !> Unrotated field array
  real, optional, intent(in) :: tstamp        !> Model timestamp
  integer, optional, intent(in) :: tile_count !> PEs per tile (default: 1)
  real, optional, intent(in) :: default_data  !> Default fill value
  integer, optional, intent(in) :: turns      !> Number of quarter-turns

  real, allocatable :: field_rot(:,:,:,:)
  integer :: qturns

  qturns = 0
  if (present(turns)) &
    qturns = modulo(turns, 4)

  if (qturns == 0) then
    call write_field(io_unit, field_md, domain, field, tstamp=tstamp, &
        tile_count=tile_count, default_data=default_data)
  else
    call allocate_rotated_array(field, [1,1,1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    call write_field(io_unit, field_md, domain, field_rot, tstamp=tstamp, &
        tile_count=tile_count, default_data=default_data)
    deallocate(field_rot)
  endif
end subroutine rotated_write_field_real_4d


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
