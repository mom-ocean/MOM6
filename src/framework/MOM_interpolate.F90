!> This module provides added functionality to the FMS temporal and spatial interpolation routines
module MOM_interpolate

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_array_transform, only : allocate_rotated_array, rotate_array
use MOM_error_handler,   only : MOM_error, FATAL
use MOM_interp_infra,    only : time_interp_extern, init_external_field=>init_extern_field
use MOM_interp_infra,    only : time_interp_external_init=>time_interp_extern_init
use MOM_interp_infra,    only : horiz_interp_type, get_external_field_info
use MOM_interp_infra,    only : run_horiz_interp, build_horiz_interp_weights
use MOM_interp_infra,    only : external_field
use MOM_time_manager, only : time_type, set_date, operator(+), operator(<), operator(>)

implicit none ; private

!> Data type used to store information about forcing datasets that are time series
!! E.g. how do we align the data in the model with the time axis in the file?
type, public :: forcing_timeseries_dataset
    character(len=200) :: file_name  !< name of file containing river flux forcing
    logical :: l_time_varying        !< .true. => forcing is dependent on model time, .false. => static forcing
    ! logical :: l_FMS_modulo        !< .true. => let FMS handle determining time level to read (e.g. for climatologies)
    type(time_type) :: data_forcing  !< convert data_forcing_year to time type
    type(time_type) :: data_start    !< convert data_start_year to time type
    type(time_type) :: data_end      !< convert data_end_year to time type
    type(time_type) :: m2d_offset    !< add to model time to get data time
end type forcing_timeseries_dataset

public :: time_interp_external, init_external_field, time_interp_external_init, get_external_field_info
public :: horiz_interp_type, run_horiz_interp, build_horiz_interp_weights
public :: external_field
public :: forcing_timeseries_set_time_type_vars
public :: map_model_time_to_forcing_time

!> Read a field based on model time, and rotate to the model domain.
interface time_interp_external
  module procedure time_interp_external_0d
  module procedure time_interp_external_2d
  module procedure time_interp_external_3d
end interface time_interp_external

contains

!> Read a scalar field based on model time.
subroutine time_interp_external_0d(field, time, data_in, verbose, scale)
  type(external_field), intent(in) :: field    !< Handle for time interpolated field
  type(time_type),   intent(in)    :: time     !< The target time for the data
  real,              intent(inout) :: data_in  !< The interpolated value in arbitrary units [A ~> a]
  logical, optional, intent(in)    :: verbose  !< If true, write verbose output for debugging
  real,    optional, intent(in)    :: scale    !< A scaling factor that new values of data_in are
                                               !! multiplied by before it is returned [A a-1 ~> 1]
  real :: data_in_pre_scale ! The input data before rescaling [a]
  real :: I_scale ! The inverse of scale [a A-1 ~> 1]

  ! Store the input value in case the scaling factor is perfectly invertable.
  data_in_pre_scale = data_in
  I_scale = 1.0
  if (present(scale)) then ; if ((scale /= 1.0) .and. (scale /= 0.0)) then
    ! Because time_interp_extern has the ability to only set some values, but no clear
    ! mechanism to determine which values have been set, the input data has to
    ! be unscaled so that it will have the right values when it is returned.
    I_scale = 1.0 / scale
    data_in = data_in * I_scale
  endif ; endif

  call time_interp_extern(field, time, data_in, verbose=verbose)

  if (present(scale)) then ; if (scale /= 1.0) then
    ! Rescale data that has been newly set and restore the scaling of unset data.
    if (data_in == I_scale * data_in_pre_scale) then
      data_in = data_in_pre_scale
    else
      data_in = scale * data_in
    endif
  endif ; endif

end subroutine time_interp_external_0d

!> Read a 2d field from an external based on model time, potentially including horizontal
!! interpolation and rotation of the data
subroutine time_interp_external_2d(field, time, data_in, interp, &
                                   verbose, horz_interp, mask_out, turns, scale)
  type(external_field), intent(in)    :: field    !< Handle for time interpolated field
  type(time_type),      intent(in)    :: time     !< The target time for the data
  real, dimension(:,:), intent(inout) :: data_in  !< The array in which to store the interpolated
                                                  !! values in arbitrary units [A ~> a]
  integer,    optional, intent(in)    :: interp   !< A flag indicating the temporal interpolation method
  logical,    optional, intent(in)    :: verbose  !< If true, write verbose output for debugging
  type(horiz_interp_type), &
              optional, intent(in)    :: horz_interp !< A structure to control horizontal interpolation
  logical, dimension(:,:), &
              optional, intent(out)   :: mask_out !< An array that is true where there is valid data
  integer,    optional, intent(in)    :: turns    !< Number of quarter turns to rotate the data
  real,       optional, intent(in)    :: scale    !< A scaling factor that new values of data_in are
                                                  !! multiplied by before it is returned [A a-1 ~> 1]

  real, allocatable :: data_in_pre_scale(:,:) ! The input data before rescaling [a]
  real, allocatable :: data_pre_rot(:,:)      ! The unscaled input data before rotation [a]
  real    :: I_scale ! The inverse of scale [a A-1 ~> 1]
  integer :: qturns ! The number of quarter turns to rotate the data
  integer :: i, j

  ! TODO: Mask rotation requires logical array rotation support
  if (present(mask_out)) &
    call MOM_error(FATAL, "Rotation of masked output not yet support")

  if (present(scale)) then ; if ((scale /= 1.0) .and. (scale /= 0.0)) then
    ! Because time_interp_extern has the ability to only set some values, but no clear mechanism
    ! to determine which values have been set, the input data has to be unscaled so that it will
    ! have the right values when it is returned.  It may be a problem for some compiler settings
    ! if there are NaNs in data_in, but they will not spread.
    if (abs(fraction(scale)) /= 1.0) then
      ! This scaling factor may not be perfectly invertable, so store the input value
      allocate(data_in_pre_scale, source=data_in)
    endif
    I_scale = 1.0 / scale
    data_in(:,:) = I_scale * data_in(:,:)
  endif ; endif

  qturns = 0 ; if (present(turns)) qturns = modulo(turns, 4)

  if (qturns == 0) then
    call time_interp_extern(field, time, data_in, interp=interp, &
                            verbose=verbose, horz_interp=horz_interp)
  else
    call allocate_rotated_array(data_in, [1,1], -qturns, data_pre_rot)
    call time_interp_extern(field, time, data_pre_rot, interp=interp, &
                            verbose=verbose, horz_interp=horz_interp)
    call rotate_array(data_pre_rot, turns, data_in)
    deallocate(data_pre_rot)
  endif

  if (present(scale)) then ; if (scale /= 1.0) then
    ! Rescale data that has been newly set and restore the scaling of unset data.
    if ((abs(fraction(scale)) /= 1.0) .and. (scale /= 0.0)) then
      do j=LBOUND(data_in,2),UBOUND(data_in,2) ; do i=LBOUND(data_in,1),UBOUND(data_in,1)
        ! This handles the case where scale is not exactly invertable for data
        ! values that have not been modified by time_interp_extern.
        if (data_in(i,j) == I_scale * data_in_pre_scale(i,j)) then
          data_in(i,j) = data_in_pre_scale(i,j)
        else
          data_in(i,j) = scale * data_in(i,j)
        endif
      enddo ; enddo
    else
      data_in(:,:) = scale * data_in(:,:)
    endif
  endif ; endif

end subroutine time_interp_external_2d


!> Read a 3d field based on model time, and rotate to the model grid
subroutine time_interp_external_3d(field, time, data_in, interp, &
                                   verbose, horz_interp, mask_out, turns, scale)
  type(external_field), intent(in)      :: field    !< Handle for time interpolated field
  type(time_type),        intent(in)    :: time     !< The target time for the data
  real, dimension(:,:,:), intent(inout) :: data_in  !< The array in which to store the interpolated
                                                    !! values in arbitrary units [A ~> a]
  integer,      optional, intent(in)    :: interp   !< A flag indicating the temporal interpolation method
  logical,      optional, intent(in)    :: verbose  !< If true, write verbose output for debugging
  type(horiz_interp_type), &
                optional, intent(in)    :: horz_interp !< A structure to control horizontal interpolation
  logical, dimension(:,:,:), &
                optional, intent(out)   :: mask_out !< An array that is true where there is valid data
  integer,      optional, intent(in)    :: turns    !< Number of quarter turns to rotate the data
  real,         optional, intent(in)    :: scale    !< A scaling factor that new values of data_in are
                                                    !! multiplied by before it is returned [A a-1 ~> 1]

  real, allocatable :: data_in_pre_scale(:,:,:) ! The input data before rescaling [a]
  real, allocatable :: data_pre_rot(:,:,:)      ! The unscaled input data before rotation [a]
  real    :: I_scale ! The inverse of scale [a A-1 ~> 1]
  integer :: qturns  ! The number of quarter turns to rotate the data
  integer :: i, j, k

  ! TODO: Mask rotation requires logical array rotation support
  if (present(mask_out)) &
    call MOM_error(FATAL, "Rotation of masked output not yet support")

  if (present(scale)) then ; if ((scale /= 1.0) .and. (scale /= 0.0)) then
    ! Because time_interp_extern has the ability to only set some values, but no clear mechanism
    ! to determine which values have been set, the input data has to be unscaled so that it will
    ! have the right values when it is returned.  It may be a problem for some compiler settings
    ! if there are NaNs in data_in, but they will not spread.
    if (abs(fraction(scale)) /= 1.0) then
      ! This scaling factor may not be perfectly invertable, so store the input value
      allocate(data_in_pre_scale, source=data_in)
    endif
    I_scale = 1.0 / scale
    data_in(:,:,:) = I_scale * data_in(:,:,:)
  endif ; endif

  qturns = 0 ; if (present(turns)) qturns = modulo(turns, 4)

  if (qturns == 0) then
    call time_interp_extern(field, time, data_in, interp=interp, &
                            verbose=verbose, horz_interp=horz_interp)
  else
    call allocate_rotated_array(data_in, [1,1,1], -qturns, data_pre_rot)
    call time_interp_extern(field, time, data_pre_rot, interp=interp, &
                            verbose=verbose, horz_interp=horz_interp)
    call rotate_array(data_pre_rot, turns, data_in)
    deallocate(data_pre_rot)
  endif

  if (present(scale)) then ; if (scale /= 1.0) then
    ! Rescale data that has been newly set and restore the scaling of unset data.
    if ((abs(fraction(scale)) /= 1.0) .and. (scale /= 0.0)) then
      do k=LBOUND(data_in,3),UBOUND(data_in,3)
        do j=LBOUND(data_in,2),UBOUND(data_in,2)
          do i=LBOUND(data_in,1),UBOUND(data_in,1)
            ! This handles the case where scale is not exactly invertable for data
            ! values that have not been modified by time_interp_extern.
            if (data_in(i,j,k) == I_scale * data_in_pre_scale(i,j,k)) then
              data_in(i,j,k) = data_in_pre_scale(i,j,k)
            else
              data_in(i,j,k) = scale * data_in(i,j,k)
            endif
          enddo
        enddo
      enddo
    else
      data_in(:,:,:) = scale * data_in(:,:,:)
    endif
  endif ; endif

end subroutine time_interp_external_3d

!> Set time_type variables in forcing_timeseries_dataset type based on integer input
!! TODO: make this part of forcing_timeseries_dataset class if OO is okay in MOM6?
subroutine forcing_timeseries_set_time_type_vars(data_start_year, data_end_year, data_ref_year, &
  model_ref_year, data_forcing_year, forcing_dataset)

  integer,                          intent(in)    :: data_start_year    !< first year of data to read
                                                                        !! (this is ignored for static forcing)
  integer,                          intent(in)    :: data_end_year      !< last year of data to read
                                                                        !! (this is ignored for static forcing)
  integer,                          intent(in)    :: data_ref_year      !< for time-varying forcing, align
                                                                        !! data_ref_year in file with
                                                                        !! model_ref_year in model
  integer,                          intent(in)    :: model_ref_year     !< for time-varying forcing, align
                                                                        !! data_ref_year in file with
                                                                        !! model_ref_year in model
  integer,                          intent(in)    :: data_forcing_year  !< for static forcing, read file at this
                                                                        !! date (this is ignored for time-varying
                                                                        !! forcing)
  type(forcing_timeseries_dataset), intent(inout) :: forcing_dataset    !< information about forcing file

  if (forcing_dataset%l_time_varying) then
    forcing_dataset%data_start = set_date(data_start_year, 1, 1)
    forcing_dataset%data_end = set_date(data_end_year, 1, 1)
    forcing_dataset%m2d_offset = set_date(data_ref_year - model_ref_year, 1, 1)
  else
    forcing_dataset%data_forcing = set_date(data_forcing_year, 1, 1)
  endif

end subroutine forcing_timeseries_set_time_type_vars

!> If necessary, apply an offset to convert from model time to forcing time and then
!! ensure result is within acceptable bounds
function map_model_time_to_forcing_time(Time, forcing_dataset)

  type(time_type),                  intent(in)  :: Time             !< Model time
  type(forcing_timeseries_dataset), intent(in)  :: forcing_dataset  !< information about forcing file
  type(time_type) :: map_model_time_to_forcing_time                 !< time to read forcing file

  if (forcing_dataset%l_time_varying) then
    map_model_time_to_forcing_time = Time + forcing_dataset%m2d_offset
    ! If Time + offset is not between data_start and data_end, use whichever of those values is closer
    if (map_model_time_to_forcing_time < forcing_dataset%data_start) &
      map_model_time_to_forcing_time = forcing_dataset%data_start
    if (map_model_time_to_forcing_time > forcing_dataset%data_end) &
      map_model_time_to_forcing_time = forcing_dataset%data_end
  else
    map_model_time_to_forcing_time = forcing_dataset%data_forcing
  endif

end function map_model_time_to_forcing_time

end module MOM_interpolate
