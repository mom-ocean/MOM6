!> This module wraps the FMS temporal and spatial interpolation routines
module MOM_interp_infra

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domain_infra,    only : MOM_domain_type, domain2d
use MOM_io,              only : axis_info
use MOM_io,              only : get_var_axes_info
use MOM_time_manager,    only : time_type
use MOM_error_handler, only : MOM_error, FATAL
use MOM_string_functions, only : lowercase
use horiz_interp_mod, only : horiz_interp_new, horiz_interp, horiz_interp_init, horiz_interp_type
use netcdf_io_mod, only : FmsNetcdfFile_t, netcdf_file_open, netcdf_file_close
use netcdf_io_mod, only : get_num_variables, get_variable_names
use time_interp_external2_mod, only : time_interp_external
use time_interp_external2_mod, only : init_external_field, time_interp_external_init
use time_interp_external2_mod, only : get_external_field_size
use time_interp_external2_mod, only : get_external_field_missing

implicit none ; private

public :: horiz_interp_type, horizontal_interp_init
public :: time_interp_extern, init_extern_field, time_interp_extern_init
public :: get_external_field_info
public :: run_horiz_interp, build_horiz_interp_weights
public :: external_field

!< Handle of an external field for interpolation
type :: external_field
  private
  integer :: id
    !< FMS ID for the interpolated field
  character(len=:), allocatable :: filename
    !< Filename containing the field values
  character(len=:), allocatable :: label
    !< Field name in the file
end type external_field

!> Read a field based on model time, and rotate to the model domain.
interface time_interp_extern
  module procedure time_interp_extern_0d
  module procedure time_interp_extern_2d
  module procedure time_interp_extern_3d
end interface time_interp_extern

!> perform horizontal interpolation of field
interface run_horiz_interp
  module procedure horiz_interp_from_weights_field2d
  module procedure horiz_interp_from_weights_field3d
end interface

!> build weights for horizontal interpolation of field
interface build_horiz_interp_weights
  module procedure build_horiz_interp_weights_2d_to_2d
end interface build_horiz_interp_weights

contains

!> Do any initialization for the horizontal interpolation
subroutine horizontal_interp_init()
  call horiz_interp_init()
end subroutine horizontal_interp_init

!> Do any initialization for the time and space interpolation infrastructure
subroutine time_interp_extern_init()
  call time_interp_external_init()
end subroutine time_interp_extern_init

!> perform horizontal interpolation of a 2d field using pre-computed weights
!! source and destination coordinates are 2d
subroutine horiz_interp_from_weights_field2d(Interp, data_in, data_out, verbose, mask_in, mask_out, &
                                             missing_value, missing_permit, err_msg)

  type(horiz_interp_type),        intent(in)  :: Interp   !< type containing interpolation options and weights
  real, dimension(:,:),           intent(in)  :: data_in  !< input data
  real, dimension(:,:),           intent(out) :: data_out !< output data
  integer,              optional, intent(in)  :: verbose  !< verbosity level
  real, dimension(:,:), optional, intent(in)  :: mask_in  !< mask for input data
  real, dimension(:,:), optional, intent(out) :: mask_out !< mask for output data
  real,                 optional, intent(in)  :: missing_value  !< A value indicating missing data
  integer,              optional, intent(in)  :: missing_permit !< number of allowed points with
                                                          !! missing value for interpolation (0-3)
  character(len=*),     optional, intent(out) :: err_msg  !< error message

  call horiz_interp(Interp, data_in, data_out, verbose, &
                    mask_in, mask_out, missing_value, missing_permit, &
                    err_msg, new_missing_handle=.true. )

end subroutine horiz_interp_from_weights_field2d

!> perform horizontal interpolation of a 3d field using pre-computed weights
!! source and destination coordinates are 2d
subroutine horiz_interp_from_weights_field3d(Interp, data_in, data_out, verbose, mask_in, mask_out, &
                                             missing_value, missing_permit, err_msg)

  type(horiz_interp_type),          intent(in)  :: Interp   !< type containing interpolation options and weights
  real, dimension(:,:,:),           intent(in)  :: data_in  !< input data
  real, dimension(:,:,:),           intent(out) :: data_out !< output data
  integer,                optional, intent(in)  :: verbose  !< verbosity level
  real, dimension(:,:,:), optional, intent(in)  :: mask_in  !< mask for input data
  real, dimension(:,:,:), optional, intent(out) :: mask_out !< mask for output data
  real,                   optional, intent(in)  :: missing_value !< A value indicating missing data
  integer,                optional, intent(in)  :: missing_permit !< number of allowed points with
                                                            !! missing value for interpolation (0-3)
  character(len=*),       optional, intent(out) :: err_msg  !< error message

  call horiz_interp(Interp, data_in, data_out, verbose, mask_in, mask_out, &
                    missing_value, missing_permit, err_msg)

end subroutine horiz_interp_from_weights_field3d


!> build horizontal interpolation weights from source grid defined by 2d lon/lat to destination grid
!! defined by 2d lon/lat
subroutine build_horiz_interp_weights_2d_to_2d(Interp, lon_in, lat_in, lon_out, lat_out, &
                                               verbose, interp_method, num_nbrs, max_dist, &
                                               src_modulo, mask_in, mask_out, &
                                               is_latlon_in, is_latlon_out)

  type(horiz_interp_type), intent(inout) :: Interp         !< type containing interpolation options and weights
  real, dimension(:,:),       intent(in) :: lon_in         !< input longitude 2d
  real, dimension(:,:),       intent(in) :: lat_in         !< input latitude 2d
  real, dimension(:,:),       intent(in) :: lon_out        !< output longitude 2d
  real, dimension(:,:),       intent(in) :: lat_out        !< output latitude 2d
  integer,          optional, intent(in) :: verbose        !< verbosity level
  character(len=*), optional, intent(in) :: interp_method  !< interpolation method
  integer,          optional, intent(in) :: num_nbrs       !< number of nearest neighbors
  real,             optional, intent(in) :: max_dist       !< maximum region of influence
  logical,          optional, intent(in) :: src_modulo     !< periodicity of E-W boundary
  real, dimension(:,:), optional, intent(in) :: mask_in    !< mask for input data
  real, dimension(:,:), optional, intent(inout) :: mask_out !< mask for output data
  logical,          optional, intent(in) :: is_latlon_in   !< input grid is regular lat/lon grid
  logical,          optional, intent(in) :: is_latlon_out  !< output grid is regular lat/lon grid

  call horiz_interp_new(Interp, lon_in, lat_in, lon_out, lat_out, &
                        verbose, interp_method, num_nbrs, max_dist, &
                        src_modulo, mask_in, mask_out, &
                        is_latlon_in, is_latlon_out)

end subroutine build_horiz_interp_weights_2d_to_2d


!> get size of an external field from field index
function get_extern_field_size(index)

  integer, intent(in) :: index         !< field index
  integer :: get_extern_field_size(4)  !< field size

  get_extern_field_size = get_external_field_size(index)

end function get_extern_field_size


!> get axes of an external field from field index
function get_extern_field_axes(field) result(axes)
  type(external_field), intent(in) :: field   !< Field handle
  type(axis_info), dimension(4) :: axes        !< Field axes

  call get_var_axes_info(field%filename, field%label, axes)
end function get_extern_field_axes


!> get missing value of an external field from field index
function get_extern_field_missing(index)

  integer, intent(in) :: index     !< field index
  real :: get_extern_field_missing !< field missing value

  get_extern_field_missing = get_external_field_missing(index)

end function get_extern_field_missing


!> Get information about the external fields.
subroutine get_external_field_info(field, size, axes, missing)
  type(external_field),                   intent(in)    :: field    !< Handle for time interpolated external
                                                                    !! field returned from a previous
                                                                    !! call to init_external_field()
  integer,         dimension(4), optional, intent(inout) :: size    !< Dimension sizes for the input data
  type(axis_info), dimension(4), optional, intent(inout) :: axes    !< Axis types for the input data
  real,                          optional, intent(inout) :: missing !< Missing value for the input data

  if (present(size)) then
    size(1:4) = get_extern_field_size(field%id)
  endif

  if (present(axes)) then
    axes(1:4) = get_extern_field_axes(field)
  endif

  if (present(missing)) then
    missing = get_extern_field_missing(field%id)
  endif

end subroutine get_external_field_info


!> Read a scalar field based on model time.
subroutine time_interp_extern_0d(field, time, data_in, verbose)
  type(external_field), intent(in) :: field    !< Handle for time interpolated field
  type(time_type),   intent(in)    :: time     !< The target time for the data
  real,              intent(inout) :: data_in  !< The interpolated value
  logical, optional, intent(in)    :: verbose  !< If true, write verbose output for debugging

  call time_interp_external(field%id, time, data_in, verbose=verbose)
end subroutine time_interp_extern_0d


!> Read a 2d field from an external based on model time, potentially including horizontal
!! interpolation and rotation of the data
subroutine time_interp_extern_2d(field, time, data_in, interp, verbose, horz_interp, mask_out)
  type(external_field), intent(in)    :: field    !< Handle for time interpolated field
  type(time_type),      intent(in)    :: time     !< The target time for the data
  real, dimension(:,:), intent(inout) :: data_in  !< The array in which to store the interpolated values
  integer,    optional, intent(in)    :: interp   !< A flag indicating the temporal interpolation method
  logical,    optional, intent(in)    :: verbose  !< If true, write verbose output for debugging
  type(horiz_interp_type), &
              optional, intent(in)    :: horz_interp !< A structure to control horizontal interpolation
  logical, dimension(:,:), &
              optional, intent(out)   :: mask_out !< An array that is true where there is valid data

  call time_interp_external(field%id, time, data_in, interp=interp, verbose=verbose, &
                            horz_interp=horz_interp, mask_out=mask_out)
end subroutine time_interp_extern_2d


!> Read a 3d field based on model time, and rotate to the model grid
subroutine time_interp_extern_3d(field, time, data_in, interp, verbose, horz_interp, mask_out)
  type(external_field),   intent(in)    :: field    !< Handle for time interpolated field
  type(time_type),        intent(in)    :: time     !< The target time for the data
  real, dimension(:,:,:), intent(inout) :: data_in  !< The array in which to store the interpolated values
  integer,      optional, intent(in)    :: interp   !< A flag indicating the temporal interpolation method
  logical,      optional, intent(in)    :: verbose  !< If true, write verbose output for debugging
  type(horiz_interp_type), &
                optional, intent(in)    :: horz_interp !< A structure to control horizontal interpolation
  logical, dimension(:,:,:), &
                optional, intent(out)   :: mask_out !< An array that is true where there is valid data

  call time_interp_external(field%id, time, data_in, interp=interp, verbose=verbose, &
                            horz_interp=horz_interp, mask_out=mask_out)
end subroutine time_interp_extern_3d


!> initialize an external field
function init_extern_field(file, fieldname, MOM_domain, domain, verbose, &
    threading, ierr, ignore_axis_atts, correct_leap_year_inconsistency) &
    result(field)

  character(len=*),         intent(in)  :: file  !< The name of the file to read
  character(len=*),         intent(in)  :: fieldname !< The name of the field in the file
  integer,        optional, intent(in)  :: threading !< A flag specifying whether the root PE reads
                                                 !! the data and broadcasts it (SINGLE_FILE) or all
                                                 !! processors read (MULTIPLE, the default).
  logical,        optional, intent(in)  :: verbose !< If true, write verbose output for debugging
  type(domain2d), optional, intent(in)  :: domain !< A domain2d type that describes the decomposition
  type(MOM_domain_type), &
                  optional, intent(in)  :: MOM_Domain !< A MOM_Domain that describes the decomposition
  integer,        optional, intent(out) :: ierr  !< Returns a non-zero error code in case of failure
  logical,        optional, intent(in)  :: ignore_axis_atts !< If present and true, do not issue a
                                                 !! fatal error if the axis Cartesian attribute is
                                                 !! not set to a recognized value.
  logical,        optional, intent(in)  :: correct_leap_year_inconsistency !< If present and true,
                                                 !! then if, (1) a calendar containing leap years
                                                 !! is in use, and (2) the modulo time period of the
                                                 !! data is an integer number of years, then map
                                                 !! a model date of Feb 29. onto a common year on Feb. 28.
  type(external_field) :: field                  !< Handle to external field

  type(FmsNetcdfFile_t) :: extern_file
    ! Local instance of netCDF file used to locate case-insensitive field name
  integer :: num_fields
    ! Number of fields in external file
  character(len=256), allocatable :: extern_fieldnames(:)
    ! List of field names in file
    ! NOTE: length should NF90_MAX_NAME, but I don't know how to read it
  character(len=:), allocatable :: label
    ! Case-insensitive match to fieldname in file
  logical :: rc
    ! Return status
  integer :: i
    ! Loop index

  field%filename = file

  ! FMS2's init_external_field is case sensitive, so we must replicate the
  !   case-insensitivity of FMS1.  This requires opening the file twice.

  rc = netcdf_file_open(extern_file, file, 'read')
  if (.not. rc) then
    call MOM_error(FATAL, 'init_extern_file: file ' // trim(file) &
        // ' could not be opened.')
  endif

  ! TODO: broadcast = .false.?
  num_fields = get_num_variables(extern_file)

  allocate(extern_fieldnames(num_fields))
  call get_variable_names(extern_file, extern_fieldnames)

  do i = 1, num_fields
    if (lowercase(extern_fieldnames(i)) == lowercase(fieldname)) then
      field%label = extern_fieldnames(i)
      exit
    endif
  enddo

  call netcdf_file_close(extern_file)

  if (.not. allocated(field%label)) then
    call MOM_error(FATAL, 'init_extern_field: field ' // trim(fieldname) &
        // ' not found in ' // trim(file) // '.')
  endif

  ! Pass to FMS2 implementation of init_external_field

  ! NOTE: external fields are currently assumed to be on-grid, which holds
  ! across the current codebase.  In the future, we may need to either enforce
  ! this or somehow relax this requirement.

  if (present(MOM_Domain)) then
    field%id = init_external_field(file, field%label, domain=MOM_domain%mpp_domain, &
             verbose=verbose, ierr=ierr, ignore_axis_atts=ignore_axis_atts, &
             correct_leap_year_inconsistency=correct_leap_year_inconsistency, &
             ongrid=.true.)
  else
    field%id = init_external_field(file, field%label, domain=domain, &
             verbose=verbose, ierr=ierr, ignore_axis_atts=ignore_axis_atts, &
             correct_leap_year_inconsistency=correct_leap_year_inconsistency, &
             ongrid=.true.)
  endif
end function init_extern_field

end module MOM_interp_infra
