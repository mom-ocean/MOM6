!> This module wraps the FMS temporal and spatial interpolation routines
module MOM_interpolate

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_array_transform, only : allocate_rotated_array, rotate_array
use MOM_error_handler, only : MOM_error, FATAL
use MOM_io_wrapper, only : axistype
use horiz_interp_mod,  only : horiz_interp_new, horiz_interp, horiz_interp_init, horiz_interp_type
use time_interp_external_mod, only : time_interp_external_fms=>time_interp_external
use time_interp_external_mod, only : init_external_field, time_interp_external_init
use time_interp_external_mod, only : get_external_field_size
use time_interp_external_mod, only : get_external_field_axes, get_external_field_missing
use time_manager_mod, only : time_type

implicit none ; private

public :: time_interp_extern, init_external_field, time_interp_external_init
public :: get_external_field_info
public :: horiz_interp_type, horiz_interp_init, horiz_interp, horiz_interp_new

!> Read a field based on model time, and rotate to the model domain.
! This inerface does not share the name time_interp_external with the module it primarily
! wraps because of errors (perhaps a bug) that arise with the PGI 19.10.0 compiler.
interface time_interp_extern
  module procedure time_interp_external_0d
  module procedure time_interp_external_2d
  module procedure time_interp_external_3d
end interface time_interp_extern

contains

!> Get information about the external fields.
subroutine get_external_field_info(field_id, size, axes, missing)
  integer,                                intent(in)    :: field_id !< The integer index of the external
                                                                    !! field returned from a previous
                                                                    !! call to init_external_field()
  integer,        dimension(4), optional, intent(inout) :: size    !< Dimension sizes for the input data
  type(axistype), dimension(4), optional, intent(inout) :: axes    !< Axis types for the input data
  real,                         optional, intent(inout) :: missing !< Missing value for the input data

  if (present(size)) then
    size(1:4) = get_external_field_size(field_id)
  endif

  if (present(axes)) then
    axes(1:4) = get_external_field_axes(field_id)
  endif

  if (present(missing)) then
    missing = get_external_field_missing(field_id)
  endif

end subroutine get_external_field_info


!> Read a scalar field based on model time.
subroutine time_interp_external_0d(field_id, time, data_in, verbose)
  integer,           intent(in)    :: field_id !< The integer index of the external field returned
                                               !! from a previous call to init_external_field()
  type(time_type),   intent(in)    :: time     !< The target time for the data
  real,              intent(inout) :: data_in  !< The interpolated value
  logical, optional, intent(in)    :: verbose  !< If true, write verbose output for debugging

  call time_interp_external_fms(field_id, time, data_in, verbose=verbose)
end subroutine time_interp_external_0d

!> Read a 2d field from an external based on model time, potentially including horizontal
!! interpolation and rotation of the data
subroutine time_interp_external_2d(field_id, time, data_in, interp, verbose, horz_interp, mask_out, turns)
  integer,              intent(in)    :: field_id !< The integer index of the external field returned
                                                  !! from a previous call to init_external_field()
  type(time_type),      intent(in)    :: time     !< The target time for the data
  real, dimension(:,:), intent(inout) :: data_in  !< The array in which to store the interpolated values
  integer,    optional, intent(in)    :: interp   !< A flag indicating the temporal interpolation method
  logical,    optional, intent(in)    :: verbose  !< If true, write verbose output for debugging
  type(horiz_interp_type), &
              optional, intent(in)    :: horz_interp !< A structure to control horizontal interpolation
  logical, dimension(:,:), &
              optional, intent(out)   :: mask_out !< An array that is true where there is valid data
  integer,    optional, intent(in)    :: turns    !< Number of quarter turns to rotate the data

  real, allocatable :: data_pre_rot(:,:) ! The input data before rotation
  integer :: qturns ! The number of quarter turns to rotate the data

  ! TODO: Mask rotation requires logical array rotation support
  if (present(mask_out)) &
    call MOM_error(FATAL, "Rotation of masked output not yet support")

  qturns = 0
  if (present(turns)) qturns = modulo(turns, 4)

  if (qturns == 0) then
    call time_interp_external_fms(field_id, time, data_in, interp=interp, &
                                  verbose=verbose, horz_interp=horz_interp)
  else
    call allocate_rotated_array(data_in, [1,1], -qturns, data_pre_rot)
    call time_interp_external_fms(field_id, time, data_pre_rot, interp=interp, &
                                  verbose=verbose, horz_interp=horz_interp)
    call rotate_array(data_pre_rot, turns, data_in)
    deallocate(data_pre_rot)
  endif
end subroutine time_interp_external_2d


!> Read a 3d field based on model time, and rotate to the model grid
subroutine time_interp_external_3d(field_id, time, data_in, interp, &
    verbose, horz_interp, mask_out, turns)
  integer,                intent(in)    :: field_id !< The integer index of the external field returned
                                                    !! from a previous call to init_external_field()
  type(time_type),        intent(in)    :: time     !< The target time for the data
  real, dimension(:,:,:), intent(inout) :: data_in  !< The array in which to store the interpolated values
  integer,      optional, intent(in)    :: interp   !< A flag indicating the temporal interpolation method
  logical,      optional, intent(in)    :: verbose  !< If true, write verbose output for debugging
  type(horiz_interp_type), &
                optional, intent(in)    :: horz_interp !< A structure to control horizontal interpolation
  logical, dimension(:,:,:), &
                optional, intent(out)   :: mask_out !< An array that is true where there is valid data
  integer,      optional, intent(in)    :: turns    !< Number of quarter turns to rotate the data

  real, allocatable :: data_pre_rot(:,:,:) ! The input data before rotation
  integer :: qturns  ! The number of quarter turns to rotate the data

  ! TODO: Mask rotation requires logical array rotation support
  if (present(mask_out)) &
    call MOM_error(FATAL, "Rotation of masked output not yet support")

  qturns = 0
  if (present(turns)) qturns = modulo(turns, 4)

  if (qturns == 0) then
    call time_interp_external_fms(field_id, time, data_in, interp=interp, &
                                  verbose=verbose, horz_interp=horz_interp)
  else
    call allocate_rotated_array(data_in, [1,1,1], -qturns, data_pre_rot)
    call time_interp_external_fms(field_id, time, data_pre_rot, interp=interp, &
                                  verbose=verbose, horz_interp=horz_interp)
    call rotate_array(data_pre_rot, turns, data_in)
    deallocate(data_pre_rot)
  endif
end subroutine time_interp_external_3d

end module MOM_interpolate
