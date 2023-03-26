!> MOM6 interface to netCDF operations
module MOM_netcdf

! This file is part of MOM6. See LICENSE.md for the license.

use, intrinsic :: iso_fortran_env, only : real32, real64

use netcdf, only : nf90_create, nf90_open, nf90_close
use netcdf, only : nf90_sync
use netcdf, only : NF90_CLOBBER, NF90_NOCLOBBER, NF90_WRITE, NF90_NOWRITE
use netcdf, only : nf90_enddef
use netcdf, only : nf90_def_dim, nf90_def_var
use netcdf, only : NF90_UNLIMITED
use netcdf, only : nf90_get_var
use netcdf, only : nf90_put_var, nf90_put_att
use netcdf, only : NF90_FLOAT, NF90_DOUBLE
use netcdf, only : nf90_strerror, NF90_NOERR
use netcdf, only : NF90_GLOBAL
use netcdf, only : nf90_inquire, nf90_inquire_dimension, nf90_inquire_variable
use netcdf, only : nf90_inq_dimids, nf90_inq_varids
use netcdf, only : NF90_MAX_NAME

use MOM_error_handler, only : MOM_error, FATAL
use MOM_io_infra, only : READONLY_FILE, WRITEONLY_FILE
use MOM_io_infra, only : APPEND_FILE, OVERWRITE_FILE

implicit none ; private

public :: netcdf_file_type
public :: netcdf_axis
public :: netcdf_field
public :: open_netcdf_file
public :: close_netcdf_file
public :: flush_netcdf_file
public :: register_netcdf_axis
public :: register_netcdf_field
public :: write_netcdf_field
public :: write_netcdf_axis
public :: write_netcdf_attribute
public :: get_netcdf_size
public :: get_netcdf_fields
public :: read_netcdf_field


!> Internal time value used to indicate an uninitialized time
real, parameter :: NULLTIME = -1
! NOTE: For now, we use the FMS-compatible value, but may change in the future.


!> netCDF file abstraction
type :: netcdf_file_type
  private
  integer :: ncid
    !< netCDF file ID
  character(len=:), allocatable :: filename
    !< netCDF filename
  logical :: define_mode
    !< True if file is in define mode.
  integer :: time_id
    !< Time axis variable ID
  real :: time
    !< Current model time
  integer :: time_level
    !< Current time level for output
end type netcdf_file_type


!> Dimension axis for a netCDF file
type :: netcdf_axis
  private
  character(len=:), allocatable, public :: label
    !< Axis label name
  real, allocatable :: points(:)
    !< Grid points along the axis
  integer :: dimid
    !< netCDF dimension ID associated with axis
  integer :: varid
    !< netCDF variable ID associated with axis
end type netcdf_axis


!> Field variable for a netCDF file
type netcdf_field
  private
  character(len=:), allocatable, public :: label
    !< Variable name
  integer :: varid
    !< netCDF variable ID for field
end type netcdf_field


!> Write values to a field of a netCDF file
interface write_netcdf_field
  module procedure write_netcdf_field_4d
  module procedure write_netcdf_field_3d
  module procedure write_netcdf_field_2d
  module procedure write_netcdf_field_1d
  module procedure write_netcdf_field_0d
end interface write_netcdf_field

contains

subroutine open_netcdf_file(handle, filename, mode)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  character(len=*), intent(in) :: filename
    !< netCDF filename
  integer, intent(in), optional :: mode
    !< Input MOM I/O mode

  integer :: io_mode
    ! MOM I/O mode
  integer :: cmode
    ! netCDF creation mode
  integer :: rc
    ! nf90_create return code
  character(len=:), allocatable :: msg
    ! netCDF error message buffer

  ! I/O configuration
  io_mode = WRITEONLY_FILE
  if (present(mode)) io_mode = mode

  ! Translate the MOM I/O config to the netCDF mode
  select case(io_mode)
    case (WRITEONLY_FILE)
      rc = nf90_create(filename, nf90_noclobber, handle%ncid)
      handle%define_mode = .true.
    case (OVERWRITE_FILE)
      rc = nf90_create(filename, nf90_clobber, handle%ncid)
      handle%define_mode = .true.
    case (APPEND_FILE)
      rc = nf90_open(filename, nf90_write, handle%ncid)
      handle%define_mode = .false.
    case (READONLY_FILE)
      rc = nf90_open(filename, nf90_nowrite, handle%ncid)
      handle%define_mode = .false.
    case default
      call MOM_error(FATAL, &
          'open_netcdf_file: File ' // filename // ': Unknown mode.')
  end select
  call check_netcdf_call(rc, 'open_netcdf_file', 'File ' // filename)

  handle%filename = filename

  ! FMS writes the filename as an attribute
  if (any(io_mode == [WRITEONLY_FILE, OVERWRITE_FILE])) &
    call write_netcdf_attribute(handle, 'filename', filename)
end subroutine open_netcdf_file


!> Close an opened netCDF file.
subroutine close_netcdf_file(handle)
  type(netcdf_file_type), intent(in) :: handle

  integer :: rc

  rc = nf90_close(handle%ncid)
  call check_netcdf_call(rc, 'close_netcdf_file', &
      'File "' // handle%filename // '"')
end subroutine close_netcdf_file


!> Flush buffered output to the netCDF file
subroutine flush_netcdf_file(handle)
  type(netcdf_file_type), intent(in) :: handle

  integer :: rc

  rc = nf90_sync(handle%ncid)
  call check_netcdf_call(rc, 'flush_netcdf_file', &
    'File "' // handle%filename // '"')
end subroutine flush_netcdf_file


!> Change netCDF mode of handle from 'define' to 'write'.
subroutine enable_netcdf_write(handle)
  type(netcdf_file_type), intent(inout) :: handle

  integer :: rc

  if (handle%define_mode) then
    rc = nf90_enddef(handle%ncid)
    call check_netcdf_call(rc, 'enable_netcdf_write', &
        'File "' // handle%filename // '"')
    handle%define_mode = .false.
  endif
end subroutine enable_netcdf_write


!> Register a netCDF variable
function register_netcdf_field(handle, label, axes, longname, units) &
    result(field)
  type(netcdf_file_type), intent(in) :: handle
    !< netCDF file handle
  character(len=*), intent(in) :: label
    !< netCDF field name in the file
  type(netcdf_axis), intent(in) :: axes(:)
    !< Axes along which field is defined
  character(len=*), intent(in) :: longname
    !< Long name of the netCDF field
  character(len=*), intent(in) :: units
    !< Field units of measurement
  type(netcdf_field) :: field
    !< netCDF field

  integer :: rc
    ! netCDF function return code
  integer :: i
    ! Loop index
  integer, allocatable :: dimids(:)
    ! netCDF dimension IDs of axes
  integer :: xtype
    ! netCDF data type

  ! Gather the axis netCDF dimension IDs
  allocate(dimids(size(axes)))
  dimids(:) = [(axes(i)%dimid, i = 1, size(axes))]

  ! Determine the corresponding netCDF data type
  ! TODO: Support a `pack`-like argument
  select case (kind(1.0))
    case (real32)
      xtype = NF90_FLOAT
    case (real64)
      xtype = NF90_DOUBLE
    case default
      call MOM_error(FATAL, "register_netcdf_axis: Unknown kind(real).")
  end select

  ! Register the field variable
  rc = nf90_def_var(handle%ncid, label, xtype, dimids, field%varid)
  call check_netcdf_call(rc, 'register_netcdf_field', &
      'File "' // handle%filename // '", Field "' // label // '"')

  ! Assign attributes

  rc = nf90_put_att(handle%ncid, field%varid, 'long_name', longname)
  call check_netcdf_call(rc, 'register_netcdf_field', &
    'Attribute "long_name" of variable "' // label // '" in file "' &
    // handle%filename // '"')

  rc = nf90_put_att(handle%ncid, field%varid, 'units', units)
  call check_netcdf_call(rc, 'register_netcdf_field', &
    'Attribute "units" of variable "' // label // '" in file "' &
    // handle%filename // '"')
end function register_netcdf_field


!> Create an axis and associated dimension in a netCDF file
function register_netcdf_axis(handle, label, units, longname, points, &
    cartesian, sense) result(axis)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  character(len=*), intent(in) :: label
    !< netCDF axis name in the file
  character(len=*), intent(in), optional :: units
    !< Axis units of measurement
  character(len=*), intent(in), optional :: longname
    !< Long name of the axis
  real, intent(in), optional :: points(:)
    !< Values of axis points (for fixed axes)
  character(len=*), intent(in), optional :: cartesian
    !< Character denoting axis direction: X, Y, Z, T, or N for none
  integer, intent(in), optional :: sense
    !< Axis direction; +1 if axis increases upward or -1 if downward

  type(netcdf_axis) :: axis
    !< netCDF coordinate axis

  integer :: xtype
    ! netCDF external data type
  integer :: rc
    ! netCDF function return code
  logical :: unlimited
    ! True if the axis is unlimited in size (e.g. time)
  integer :: axis_size
    ! Either the number of points in the axis, or unlimited flag
  integer :: axis_sense
    ! Axis direction; +1 if axis increases upward or -1 if downward
  character(len=:), allocatable :: sense_attr
    ! CF-compiant value of sense attribute (as 'positive')

  ! Create the axis dimension
  unlimited = .false.
  if (present(cartesian)) then
    if (cartesian == 'T') unlimited = .true.
  endif

  ! Either the axis is explicitly set with data or is declared as unlimited
  if (present(points) .eqv. unlimited) then
    call MOM_error(FATAL, &
        "Axis must either have explicit points or be a time axis ('T').")
  endif

  if (present(points)) then
    axis_size = size(points)
    allocate(axis%points(axis_size))
    axis%points(:) = points(:)
  else
    axis_size = NF90_UNLIMITED
  endif

  rc = nf90_def_dim(handle%ncid, label, axis_size, axis%dimid)
  call check_netcdf_call(rc, 'register_netcdf_axis', &
      'Dimension "' // label // '" in file "' // handle%filename // '"')

  ! Determine the corresponding netCDF data type
  ! TODO: Support a `pack`-like argument
  select case (kind(1.0))
    case (real32)
      xtype = NF90_FLOAT
    case (real64)
      xtype = NF90_DOUBLE
    case default
      call MOM_error(FATAL, "register_netcdf_axis: Unknown kind(real).")
  end select

  ! Create a variable corresponding to the axis
  rc = nf90_def_var(handle%ncid, label, xtype, axis%dimid, axis%varid)
  call check_netcdf_call(rc, 'register_netcdf_axis', &
      'Variable ' // label // ' in file ' // handle%filename)

  ! Define the time axis, if available
  if (unlimited) then
    handle%time_id = axis%varid
    handle%time_level = 0
    handle%time = NULLTIME
  endif

  ! Assign attributes if present
  if (present(longname)) then
    rc = nf90_put_att(handle%ncid, axis%varid, 'long_name', longname)
    call check_netcdf_call(rc, 'register_netcdf_axis', &
      'Attribute ''long_name'' of variable ' // label // ' in file ' &
      // handle%filename)
  endif

  if (present(units)) then
    rc = nf90_put_att(handle%ncid, axis%varid, 'units', units)
    call check_netcdf_call(rc, 'register_netcdf_axis', &
      'Attribute ''units'' of variable ' // label // ' in file ' &
      // handle%filename)
  endif

  if (present(cartesian)) then
    rc = nf90_put_att(handle%ncid, axis%varid, 'cartesian_axis', cartesian)
    call check_netcdf_call(rc, 'register_netcdf_axis', &
      'Attribute ''cartesian_axis'' of variable ' // label // ' in file ' &
      // handle%filename)
  endif

  axis_sense = 0
  if (present(sense)) axis_sense = sense

  if (axis_sense /= 0) then
    select case (axis_sense)
      case (1)
        sense_attr = 'up'
      case (-1)
        sense_attr = 'down'
      case default
        call MOM_error(FATAL, 'register_netcdf_axis: sense must be either ' &
          // '0, 1, or -1.')
    end select
    rc = nf90_put_att(handle%ncid, axis%varid, 'positive', sense_attr)
    call check_netcdf_call(rc, 'register_netcdf_axis', &
      'Attribute "positive" of variable "' // label // '" in file "' &
      // handle%filename // '"')
  endif
end function register_netcdf_axis


!> Write a 4D array to a compatible netCDF field
subroutine write_netcdf_field_4d(handle, field, values, time)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  type(netcdf_field), intent(in) :: field
    !< Field metadata
  real, intent(in) :: values(:,:,:,:)
    !< Field values
  real, intent(in), optional :: time
    !< Timestep index to write data

  integer :: rc
    ! netCDF return code
  integer :: start(5)
    ! Start indices, if timestep is included

  ! Verify write mode
  if (handle%define_mode) &
    call enable_netcdf_write(handle)

  if (present(time)) then
    call update_netcdf_timestep(handle, time)
    start(:4) = 1
    start(5) = handle%time_level
    rc = nf90_put_var(handle%ncid, field%varid, values, start)
  else
    rc = nf90_put_var(handle%ncid, field%varid, values)
  endif
  call check_netcdf_call(rc, 'write_netcdf_file', &
      'File "' // handle%filename // '", Field "' // field%label // '"')
end subroutine write_netcdf_field_4d


!> Write a 3D array to a compatible netCDF field
subroutine write_netcdf_field_3d(handle, field, values, time)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  type(netcdf_field), intent(in) :: field
    !< Field metadata
  real, intent(in) :: values(:,:,:)
    !< Field values
  real, intent(in), optional :: time
    !< Timestep index to write data

  integer :: rc
    ! netCDF return code
  integer :: start(4)
    ! Start indices, if timestep is included

  ! Verify write mode
  if (handle%define_mode) &
    call enable_netcdf_write(handle)

  if (present(time)) then
    call update_netcdf_timestep(handle, time)
    start(:3) = 1
    start(4) = handle%time_level
    rc = nf90_put_var(handle%ncid, field%varid, values, start)
  else
    rc = nf90_put_var(handle%ncid, field%varid, values)
  endif
  call check_netcdf_call(rc, 'write_netcdf_file', &
      'File "' // handle%filename // '", Field "' // field%label // '"')
end subroutine write_netcdf_field_3d


!> Write a 2D array to a compatible netCDF field
subroutine write_netcdf_field_2d(handle, field, values, time)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  type(netcdf_field), intent(in) :: field
    !< Field metadata
  real, intent(in) :: values(:,:)
    !< Field values
  real, intent(in), optional :: time
    !< Timestep index to write data

  integer :: rc
    ! netCDF return code
  integer :: start(3)
    ! Start indices, if timestep is included

  ! Verify write mode
  if (handle%define_mode) &
    call enable_netcdf_write(handle)

  if (present(time)) then
    call update_netcdf_timestep(handle, time)
    start(:2) = 1
    start(3) = handle%time_level
    rc = nf90_put_var(handle%ncid, field%varid, values, start)
  else
    rc = nf90_put_var(handle%ncid, field%varid, values)
  endif
  call check_netcdf_call(rc, 'write_netcdf_file', &
      'File "' // handle%filename // '", Field "' // field%label // '"')
end subroutine write_netcdf_field_2d


!> Write a 1D array to a compatible netCDF field
subroutine write_netcdf_field_1d(handle, field, values, time)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  type(netcdf_field), intent(in) :: field
    !< Field metadata
  real, intent(in) :: values(:)
    !< Field values
  real, intent(in), optional :: time
    !< Timestep index to write data

  integer :: rc
    ! netCDF return code
  integer :: start(2)
    ! Start indices, if timestep is included

  ! Verify write mode
  if (handle%define_mode) &
    call enable_netcdf_write(handle)

  if (present(time)) then
    call update_netcdf_timestep(handle, time)
    start(1) = 1
    start(2) = handle%time_level
    rc = nf90_put_var(handle%ncid, field%varid, values, start)
  else
    rc = nf90_put_var(handle%ncid, field%varid, values)
  endif
  call check_netcdf_call(rc, 'write_netcdf_file', &
      'File "' // handle%filename // '", Field "' // field%label // '"')
end subroutine write_netcdf_field_1d


!> Write a scalar to a compatible netCDF field
subroutine write_netcdf_field_0d(handle, field, scalar, time)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  type(netcdf_field), intent(in) :: field
    !< Field metadata
  real, intent(in) :: scalar
    !< Field values
  real, intent(in), optional :: time
    !< Timestep index to write data

  integer :: rc
    ! netCDF return code
  integer :: start(1)
    ! Start indices, if timestep is included

  ! Verify write mode
  if (handle%define_mode) &
    call enable_netcdf_write(handle)

  if (present(time)) then
    call update_netcdf_timestep(handle, time)
    start(1) = handle%time_level
    rc = nf90_put_var(handle%ncid, field%varid, scalar, start)
  else
    rc = nf90_put_var(handle%ncid, field%varid, scalar)
  endif
  call check_netcdf_call(rc, 'write_netcdf_file', &
      'File "' // handle%filename // '", Field "' // field%label // '"')
end subroutine write_netcdf_field_0d


!> Write axis points to associated netCDF variable
subroutine write_netcdf_axis(handle, axis)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  type(netcdf_axis), intent(in) :: axis
    !< field variable

  integer :: rc
    ! netCDF return code

  ! Verify write mode
  if (handle%define_mode) &
    call enable_netcdf_write(handle)

  rc = nf90_put_var(handle%ncid, axis%varid, axis%points)
  call check_netcdf_call(rc, 'write_netcdf_axis', &
      'File "' // handle%filename // '", Axis "' // axis%label // '"')
end subroutine write_netcdf_axis


!> Write a global attribute to a netCDF file
subroutine write_netcdf_attribute(handle, label, attribute)
  type(netcdf_file_type), intent(in) :: handle
    !< netCDF file handle
  character(len=*), intent(in) :: label
    !< File attribute
  character(len=*), intent(in) :: attribute
    !< File attribute value

  integer :: rc
    ! netCDF return code

  rc = nf90_put_att(handle%ncid, NF90_GLOBAL, label, attribute)
  call check_netcdf_call(rc, 'write_netcdf_attribute', &
      'File "' // handle%filename // '", Attribute "' // label // '"')
end subroutine write_netcdf_attribute


! This is a thin interface to nf90_inquire, designed to mirror the existing
! I/O API.  A more axis-aware system might not need this, but for now it's here
!> Get the number of dimensions, variables, and timesteps in a netCDF file
subroutine get_netcdf_size(handle, ndims, nvars, nsteps)
  type(netcdf_file_type), intent(in) :: handle
    !< netCDF input file
  integer, intent(out), optional :: ndims
    !< number of dimensions in the file
  integer, intent(out), optional :: nvars
    !< number of variables in the file
  integer, intent(out), optional :: nsteps
    !< number of values in the file's unlimited axis

  integer :: rc
    ! netCDF return code
  integer :: unlimited_dimid
    ! netCDF dimension ID for unlimited time axis

  rc = nf90_inquire(handle%ncid, &
      nDimensions=ndims, &
      nVariables=nvars, &
      unlimitedDimId=unlimited_dimid &
  )
  call check_netcdf_call(rc, 'get_netcdf_size', &
      'File "' // handle%filename // '"')

  rc = nf90_inquire_dimension(handle%ncid, unlimited_dimid, len=nsteps)
  call check_netcdf_call(rc, 'get_netcdf_size', &
      'File "' // handle%filename // '"')
end subroutine get_netcdf_size


!> Get the metadata of the registered fields in a netCDF file
subroutine get_netcdf_fields(handle, axes, fields)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  type(netcdf_axis), intent(inout), allocatable :: axes(:)
    !< netCDF file axes
  type(netcdf_field), intent(inout), allocatable :: fields(:)
    !< netCDF file fields

  integer :: ndims
    ! Number of netCDF dimensions
  integer :: nvars
    ! Number of netCDF dimensions
  type(netcdf_field), allocatable :: vars(:)
    ! netCDF variables in handle
  integer :: nfields
    ! Number of fields in the file (i.e. non-axis variables)
  integer, allocatable :: dimids(:)
    ! netCDF dimension IDs of file
  integer, allocatable :: varids(:)
    ! netCDF variable IDs of file
  integer :: unlim_dimid
    ! netCDF dimension ID for the unlimited axis variable, if present
  integer :: unlim_index
    ! Index of the unlimited axis in axes(:), if present
  character(len=NF90_MAX_NAME) :: label
    ! Current dimension or variable label
  integer :: len
    ! Current dimension length
  integer :: rc
    ! netCDF return code
  integer :: grp_ndims, grp_nvars
    ! Group-based counts for nf90_inq_* (unused)
  logical :: is_axis
    ! True if the current variable is an axis
  integer :: i, j, n

  integer, save :: no_parent_groups = 0
    ! Flag indicating exclusion of parent groups in netCDF file
    ! NOTE: This must be passed as a variable, and cannot be declared as a
    !   parameter.

  rc = nf90_inquire(handle%ncid, &
      nDimensions=ndims, &
      nVariables=nvars, &
      unlimitedDimId=unlim_dimid &
  )
  call check_netcdf_call(rc, 'get_netcdf_fields', &
      'File "' // handle%filename // '"')

  allocate(dimids(ndims))
  rc = nf90_inq_dimids(handle%ncid, grp_ndims, dimids, no_parent_groups)
  call check_netcdf_call(rc, 'get_netcdf_fields', &
      'File "' // handle%filename // '"')

  allocate(varids(nvars))
  rc = nf90_inq_varids(handle%ncid, grp_nvars, varids)
  call check_netcdf_call(rc, 'get_netcdf_fields', &
      'File "' // trim(handle%filename) // '"')

  ! Initialize unlim_index with an unreachable value (outside [1,ndims])
  unlim_index = -1

  allocate(axes(ndims))
  do i = 1, ndims
    rc = nf90_inquire_dimension(handle%ncid, dimids(i), name=label, len=len)
    call check_netcdf_call(rc, 'get_netcdf_fields', &
        'File "' // trim(handle%filename) // '"')

    ! Check for the unlimited axis
    if (dimids(i) == unlim_dimid) unlim_index = i

    axes(i)%dimid = dimids(i)
    axes(i)%label = trim(label)
    allocate(axes(i)%points(len))
  enddo

  ! We cannot know if every axis also has a variable representation, so we
  ! over-allocate vars(:) and fill as fields are identified.
  allocate(vars(nvars))

  nfields = 0
  do i = 1, nvars
    rc = nf90_inquire_variable(handle%ncid, varids(i), name=label)
    call check_netcdf_call(rc, 'get_netcdf_fields', &
        'File "' // trim(handle%filename) // '"')

    ! Check if variable is an axis
    is_axis = .false.
    do j = 1, ndims
      if (label == axes(j)%label) then
        rc = nf90_get_var(handle%ncid, varids(i), axes(j)%points)
        call check_netcdf_call(rc, 'get_netcdf_fields', &
            'File "' // trim(handle%filename) // '"')
        axes(j)%varid = varids(i)

        if (j == unlim_index) then
          handle%time_id = varids(i)
          handle%time_level = size(axes(j)%points)
          handle%time = NULLTIME
        endif

        is_axis = .true.
        exit
      endif
    enddo
    if (is_axis) cycle

    nfields = nfields + 1
    vars(nfields)%label = trim(label)
    vars(nfields)%varid = varids(i)
  enddo

  allocate(fields(nfields))
  fields(:) = vars(:nfields)
end subroutine get_netcdf_fields


!> Read the values of a field from a netCDF file
subroutine read_netcdf_field(handle, field, values, bounds)
  type(netcdf_file_type), intent(in) :: handle
  type(netcdf_field), intent(in) :: field
  real, intent(out) :: values(:,:)
  integer, optional, intent(in) :: bounds(2,2)

  integer :: rc
    ! netCDF return code
  integer :: istart(2)
    ! Axis start index
  integer :: icount(2)
    ! Axis index count

  if (present(bounds)) then
    istart(:) = bounds(1,:)
    icount(:) = bounds(2,:) - bounds(1,:) + 1
    rc = nf90_get_var(handle%ncid, field%varid, values, start=istart, count=icount)
  else
    rc = nf90_get_var(handle%ncid, field%varid, values)
  endif
  call check_netcdf_call(rc, 'read_netcdf_field', &
      'File "' // trim(handle%filename) // '", Field "' // trim(field%label) // '"')
end subroutine read_netcdf_field


!> Set the current timestep of an open netCDF file
subroutine update_netcdf_timestep(handle, time)
  type(netcdf_file_type), intent(inout) :: handle
    !< netCDF file handle
  real, intent(in) :: time
    !< New model time

  integer :: start(1)
    !< Time axis start index array
  integer :: rc
    !< netCDF return code

  if (time > handle%time + epsilon(time)) then
    handle%time = time
    handle%time_level = handle%time_level + 1

    ! Write new value to time axis
    start = [handle%time_level]
    rc = nf90_put_var(handle%ncid, handle%time_id, time, start=start)
    call check_netcdf_call(rc, 'update_netcdf_timestep', &
        'File "' // handle%filename // '"')
  endif
end subroutine update_netcdf_timestep


!> Check netCDF function return codes, report the error log, and abort the run.
subroutine check_netcdf_call(ncerr, header, message)
  integer, intent(in) :: ncerr
    !< netCDF error code
  character(len=*), intent(in) :: header
    !< Message header (usually calling subroutine)
  character(len=*), intent(in) :: message
    !< Error message (usually action which instigated the error)

  character(len=:), allocatable :: errmsg
    ! Full error message, including netCDF message

  if (ncerr /= nf90_noerr) then
    errmsg = trim(header) // ": " // trim(message) // new_line('/') &
      // trim(nf90_strerror(ncerr))
    call MOM_error(FATAL, errmsg)
  endif
end subroutine check_netcdf_call

end module MOM_netcdf
