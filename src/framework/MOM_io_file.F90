!> This module contains the MOM file handler types
module MOM_io_file

! This file is part of MOM6. See LICENSE.md for the license.

use, intrinsic :: iso_fortran_env, only : int64

use MOM_domains, only : MOM_domain_type, domain1D
use MOM_domains, only : clone_MOM_domain
use MOM_domains, only : deallocate_MOM_domain
use MOM_io_infra, only : file_type, get_file_info, get_file_fields
use MOM_io_infra, only : open_file, close_file, flush_file
use MOM_io_infra, only : fms2_file_is_open => file_is_open
use MOM_io_infra, only : fieldtype
use MOM_io_infra, only : get_file_times, axistype
use MOM_io_infra, only : write_field, write_metadata
use MOM_io_infra, only : get_field_atts
use MOM_io_infra, only : read_field_chksum
use MOM_io_infra, only : SINGLE_FILE

use MOM_hor_index, only : hor_index_type
use MOM_hor_index, only : hor_index_init

use MOM_netcdf, only : netcdf_file_type
use MOM_netcdf, only : netcdf_axis
use MOM_netcdf, only : netcdf_field
use MOM_netcdf, only : open_netcdf_file
use MOM_netcdf, only : close_netcdf_file
use MOM_netcdf, only : flush_netcdf_file
use MOM_netcdf, only : register_netcdf_axis
use MOM_netcdf, only : register_netcdf_field
use MOM_netcdf, only : write_netcdf_field
use MOM_netcdf, only : write_netcdf_axis
use MOM_netcdf, only : write_netcdf_attribute
use MOM_netcdf, only : get_netcdf_size
use MOM_netcdf, only : get_netcdf_fields
use MOM_netcdf, only : get_netcdf_filename
use MOM_netcdf, only : read_netcdf_field

use MOM_error_handler, only : MOM_error, FATAL
use MOM_error_handler, only : is_root_PE

implicit none ; private

public :: MOM_file
public :: MOM_infra_file
public :: MOM_netcdf_file
public :: MOM_axis
public :: MOM_field


! Internal types

! NOTE: MOM_axis and MOM_field do not contain the actual axes and fields stored
! in the file.  They are very thin wrappers to the keys (as strings) used to
! reference the associated object inside of the MOM_file.

!> Handle for axis in MOM file
type :: MOM_axis
  character(len=:), allocatable :: label
    !< Identifier for the axis in handle's list
end type MOM_axis


!> Linked list of framework axes
type :: axis_list_infra
  private
  type(axis_node_infra), pointer :: head => null()
    !< Head of axis linked list
  type(axis_node_infra), pointer :: tail => null()
    !< Tail of axis linked list
contains
  !> Initialize the framework axis list
  procedure :: init => initialize_axis_list_infra
  !> Append a new axis to the framework axis list
  procedure :: append => append_axis_list_infra
  !> Get an axis from the framework axis list
  procedure :: get => get_axis_list_infra
  !> Deallocate the framework axis list
  procedure :: finalize => finalize_axis_list_infra
end type axis_list_infra


!> Framework axis linked list node
type :: axis_node_infra
  private
  character(len=:), allocatable :: label
    !< Axis identifier
  type(axis_node_infra), pointer :: next => null()
    !< Pointer to next axis node
  type(axistype) :: axis
    !< Axis node contents
end type axis_node_infra


!> Linked list of framework axes
type :: axis_list_nc
  private
  type(axis_node_nc), pointer :: head => null()
    !< Head of axis linked list
  type(axis_node_nc), pointer :: tail => null()
    !< Tail of axis linked list
contains
  !> Initialize the netCDF axis list
  procedure :: init => initialize_axis_list_nc
  !> Append a new axis to the netCDF axis list
  procedure :: append => append_axis_list_nc
  !> Get an axis from the netCDF axis list
  procedure :: get => get_axis_list_nc
  !> Deallocate the netCDF axis list
  procedure :: finalize => finalize_axis_list_nc
end type axis_list_nc


!> Framework axis linked list node
type :: axis_node_nc
  private
  character(len=:), allocatable :: label
    !< Axis identifier
  type(axis_node_nc), pointer :: next => null()
    !< Pointer to next axis node
  type(netcdf_axis) :: axis
    !< Axis node contents
end type axis_node_nc


!> Handle for field in MOM file
type :: MOM_field
  character(len=:), allocatable :: label
    !< Identifier for the field in the handle's list
  real :: conversion
    !< A factor to use to rescale the field before output [a A-1 ~> 1]
end type MOM_field


!> Linked list of framework fields
type :: field_list_infra
  private
  type(field_node_infra), pointer :: head => null()
    !< Head of field linked list
  type(field_node_infra), pointer :: tail => null()
    !< Tail of field linked list
contains
  !> Initialize the framework field list
  procedure :: init => initialize_field_list_infra
  !> Append a new axis to the framework field list
  procedure :: append => append_field_list_infra
  !> Get an axis from the framework field list
  procedure :: get => get_field_list_infra
  !> Deallocate the framework field list
  procedure :: finalize => finalize_field_list_infra
end type field_list_infra


!> Framework field linked list node
type :: field_node_infra
  private
  character(len=:), allocatable :: label
    !< Field identifier
  type(fieldtype) :: field
    !< Field node contents
  type(field_node_infra), pointer :: next => null()
    !< Pointer to next field node
end type field_node_infra


!> Linked list of framework fields
type :: field_list_nc
  private
  type(field_node_nc), pointer :: head => null()
    !< Head of field linked list
  type(field_node_nc), pointer :: tail => null()
    !< Tail of field linked list
contains
  !> Initialize the netCDF field list
  procedure :: init => initialize_field_list_nc
  !> Append a new axis to the netCDF field list
  procedure :: append => append_field_list_nc
  !> Get an axis from the netCDF field list
  procedure :: get => get_field_list_nc
  !> Deallocate the netCDF field list
  procedure :: finalize => finalize_field_list_nc
end type field_list_nc


!> Framework field linked list node
type :: field_node_nc
  private
  character(len=:), allocatable :: label
    !< Field identifier
  type(netcdf_field) :: field
    !< Field node contents
  type(field_node_nc), pointer :: next => null()
    !< Pointer to next field node
end type field_node_nc


!> Generic MOM file abstraction for common operations
type, abstract :: MOM_file
  private

  contains

  !> Open a file and connect to the MOM_file object
  procedure(i_open_file), deferred :: open
  !> Close the MOM file
  procedure(i_close_file), deferred :: close
  !> Flush buffered output to the MOM file
  procedure(i_flush_file), deferred :: flush

  !> Register an axis to the MOM file
  procedure(i_register_axis), deferred :: register_axis
  !> Register a field to the MOM file
  procedure(i_register_field), deferred :: register_field
  !> Write metadata to the MOM file
  procedure(i_write_attribute), deferred :: write_attribute

  !> Write field to a MOM file
  generic :: write_field => &
      write_field_4d, &
      write_field_3d, &
      write_field_2d, &
      write_field_1d, &
      write_field_0d, &
      write_field_axis

  !> Write a 4D field to the MOM file
  procedure(i_write_field_4d), deferred :: write_field_4d
  !> Write a 3D field to the MOM file
  procedure(i_write_field_3d), deferred :: write_field_3d
  !> Write a 2D field to the MOM file
  procedure(i_write_field_2d), deferred :: write_field_2d
  !> Write a 1D field to the MOM file
  procedure(i_write_field_1d), deferred :: write_field_1d
  !> Write a 0D field to the MOM file
  procedure(i_write_field_0d), deferred :: write_field_0d
  !> Write an axis field to the MOM file
  procedure(i_write_field_axis), deferred :: write_field_axis

  !> Return true if MOM file has been opened
  procedure(i_file_is_open), deferred :: file_is_open
  !> Return number of dimensions, variables, or time levels in a MOM file
  procedure(i_get_file_info), deferred :: get_file_info
  !> Get field objects from a MOM file
  procedure(i_get_file_fields), deferred :: get_file_fields
  !> Get attributes from a field
  procedure(i_get_field_atts), deferred :: get_field_atts
  !> Get checksum from a field
  procedure(i_read_field_chksum), deferred :: read_field_chksum
end type MOM_file


!> MOM file from the supporting framework ("infra") layer
type, extends(MOM_file) :: MOM_infra_file
  private

  type(MOM_domain_type), public, pointer :: domain => null()
    !< Internal domain used for single-file IO

  ! NOTE: This will be made private after the API transition
  type(file_type), public :: handle_infra
    !< Framework-specific file handler content
  type(axis_list_infra) :: axes
    !< List of axes in file
  type(field_list_infra) :: fields
    !< List of fields in file

  contains

  !> Open a framework file and connect to the MOM_file object
  procedure :: open => open_file_infra
  !> Close the MOM framework file
  procedure :: close => close_file_infra
  !> Flush buffered output to the MOM framework file
  procedure :: flush => flush_file_infra

  !> Register an axis to the MOM framework file
  procedure :: register_axis => register_axis_infra
  !> Register a field to the MOM framework file
  procedure :: register_field => register_field_infra
  !> Write global metadata to the MOM framework file
  procedure :: write_attribute => write_attribute_infra

  !> Write a 4D field to the MOM framework file
  procedure :: write_field_4d => write_field_4d_infra
  !> Write a 3D field to the MOM framework file
  procedure :: write_field_3d => write_field_3d_infra
  !> Write a 2D field to the MOM framework file
  procedure :: write_field_2d => write_field_2d_infra
  !> Write a 1D field to the MOM framework file
  procedure :: write_field_1d => write_field_1d_infra
  !> Write a 0D field to the MOM framework file
  procedure :: write_field_0d => write_field_0d_infra
  !> Write an axis field to the MOM framework file
  procedure :: write_field_axis => write_field_axis_infra

  !> Return true if MOM infra file has been opened
  procedure :: file_is_open => file_is_open_infra
  !> Return number of dimensions, variables, or time levels in a MOM infra file
  procedure :: get_file_info => get_file_info_infra
  !> Get field metadata from a MOM infra file
  procedure :: get_file_fields => get_file_fields_infra
  !> Get attributes from a field
  procedure :: get_field_atts => get_field_atts_infra
  !> Get checksum from a field
  procedure :: read_field_chksum => read_field_chksum_infra

  ! MOM_infra_file methods
  ! NOTE: These could naturally reside in MOM_file but is currently not needed.

  !> Get time levels of a MOM framework file
  procedure :: get_file_times => get_file_times_infra

  !> Get the fields as fieldtypes from a file
  procedure :: get_file_fieldtypes
  ! NOTE: This is provided to support the legacy API and may be removed.
end type MOM_infra_file


!> MOM file using netCDF backend
type, extends(MOM_file) :: MOM_netcdf_file
  private

  !> Framework-specific file handler content
  type(netcdf_file_type) :: handle_nc
  !> List of netCDF axes
  type(axis_list_nc) :: axes
  !> List of netCDF fields
  type(field_list_nc) :: fields
  !> True if the file has been opened
  logical :: is_open = .false.
  !> True if I/O content is domain-decomposed
  logical :: domain_decomposed = .false.
  !> True if I/O content is domain-decomposed
  type(hor_index_type) :: HI

  contains

  !> Open a framework file and connect to the MOM_netcdf_file object
  procedure :: open => open_file_nc
  !> Close the MOM netcdf file
  procedure :: close => close_file_nc
  !> Flush buffered output to the MOM netcdf file
  procedure :: flush => flush_file_nc

  !> Register an axis to the MOM netcdf file
  procedure :: register_axis => register_axis_nc
  !> Register a field to the MOM netcdf file
  procedure :: register_field => register_field_nc
  !> Write global metadata to the MOM netcdf file
  procedure :: write_attribute => write_attribute_nc

  !> Write a 4D field to the MOM netcdf file
  procedure :: write_field_4d => write_field_4d_nc
  !> Write a 3D field to the MOM netcdf file
  procedure :: write_field_3d => write_field_3d_nc
  !> Write a 2D field to the MOM netcdf file
  procedure :: write_field_2d => write_field_2d_nc
  !> Write a 1D field to the MOM netcdf file
  procedure :: write_field_1d => write_field_1d_nc
  !> Write a 0D field to the MOM netcdf file
  procedure :: write_field_0d => write_field_0d_nc
  !> Write an axis field to the MOM netcdf file
  procedure :: write_field_axis => write_field_axis_nc

  !> Return true if MOM netcdf file has been opened
  procedure :: file_is_open => file_is_open_nc
  !> Return number of dimensions, variables, or time levels in a MOM netcdf file
  procedure :: get_file_info => get_file_info_nc
  !> Get field metadata from a MOM netcdf file
  procedure :: get_file_fields => get_file_fields_nc
  !> Get attributes from a netCDF field
  procedure :: get_field_atts => get_field_atts_nc
  !> Get checksum from a netCDF field
  procedure :: read_field_chksum => read_field_chksum_nc

  ! NOTE: These are currently exclusive to netCDF I/O but could be generalized
  !> Read the values of a netCDF field
  procedure :: read => get_field_nc
  !> Update the axes and fields descriptors of a MOM netCDF file
  procedure :: update => update_file_contents_nc
end type MOM_netcdf_file


interface
  !> Interface for opening a MOM file
  subroutine i_open_file(handle, filename, action, MOM_domain, threading, fileset)
    import :: MOM_file, MOM_domain_type

    class(MOM_file), intent(inout) :: handle
      !< The handle for the opened file
    character(len=*), intent(in) :: filename
      !< The path name of the file being opened
    integer, optional, intent(in) :: action
      !< A flag indicating whether the file can be read or written to and how
      !! to handle existing files.  The default is WRITE_ONLY.
    type(MOM_domain_type), optional, intent(in) :: MOM_Domain
      !< A MOM_Domain that describes the decomposition
    integer, optional, intent(in) :: threading
      !< A flag indicating whether one (SINGLE_FILE) or multiple PEs (MULTIPLE)
      !! participate in I/O.  With the default, the root PE does I/O.
    integer, optional, intent(in) :: fileset
      !< A flag indicating whether multiple PEs doing I/O due to
      !! threading=MULTIPLE write to the same file (SINGLE_FILE) or to one file
      !! per PE (MULTIPLE, the default).
  end subroutine i_open_file


  !> Interface for closing a MOM file
  subroutine i_close_file(handle)
    import :: MOM_file
    class(MOM_file), intent(inout) :: handle
      !< The MOM file to be closed
  end subroutine i_close_file


  !> Interface for flushing I/O in a MOM file
  subroutine i_flush_file(handle)
    import :: MOM_file
    class(MOM_file), intent(in) :: handle
      !< The MOM file to be flushed
  end subroutine i_flush_file


  !> Interface to register an axis to a MOM file
  function i_register_axis(handle, label, units, longname, cartesian, sense, &
      domain, data, edge_axis, calendar) result(axis)
    import :: MOM_file, MOM_axis, domain1D

    class(MOM_file), intent(inout) :: handle
      !< Handle for a file that is open for writing
    character(len=*), intent(in) :: label
      !< The name in the file of this axis
    character(len=*), intent(in) :: units
      !< The units of this axis
    character(len=*), intent(in) :: longname
      !< The long description of this axis
    character(len=*), optional, intent(in) :: cartesian
      !< A variable indicating which direction this axis corresponds with.
      !! Valid values include 'X', 'Y', 'Z', 'T', and 'N' for none.
    integer, optional, intent(in) :: sense
      !< This is 1 for axes whose values increase upward, or -1 if they
      !! increase downward.
    type(domain1D), optional, intent(in) :: domain
      !< The domain decomposion for this axis
    real, dimension(:), optional, intent(in) :: data
      !< The coordinate values of the points on this axis
    logical, optional, intent(in) :: edge_axis
      !< If true, this axis marks an edge of the tracer cells
    character(len=*), optional, intent(in) :: calendar
      !< The name of the calendar used with a time axis
    type(MOM_axis) :: axis
      !< IO handle for axis in MOM_file
  end function i_register_axis


  !> Interface to register a field to a netCDF file
  function i_register_field(handle, axes, label, units, longname, &
      pack, standard_name, checksum, conversion) result(field)
    import :: MOM_file, MOM_axis, MOM_field, int64
    class(MOM_file), intent(inout) :: handle
        !< Handle for a file that is open for writing
    type(MOM_axis), intent(in) :: axes(:)
      !< Handles for the axis used for this variable
    character(len=*), intent(in) :: label
      !< The name in the file of this variable
    character(len=*), intent(in) :: units
      !< The units of this variable
    character(len=*), intent(in) :: longname
      !< The long description of this variable
    integer, optional, intent(in) :: pack
      !< A precision reduction factor with which the variable.  The default, 1,
      !! has no reduction, but 2 is not uncommon.
    character(len=*), optional, intent(in) :: standard_name
      !< The standard (e.g., CMOR) name for this variable
    integer(kind=int64), dimension(:), optional, intent(in) :: checksum
      !< Checksum values that can be used to verify reads.
    real, optional, intent(in) :: conversion
      !< A factor to use to rescale the field before output [a A-1 ~> 1]
    type(MOM_field) :: field
      !< IO handle for field in MOM_file
  end function i_register_field


  !> Interface for writing global metata to a MOM file
  subroutine i_write_attribute(handle, name, attribute)
    import :: MOM_file
    class(MOM_file), intent(in) :: handle
      !< Handle for a file that is open for writing
    character(len=*), intent(in) :: name
      !< The name in the file of this global attribute
    character(len=*), intent(in) :: attribute
      !< The value of this attribute
  end subroutine i_write_attribute


  !> Interface to write_field_4d()
  subroutine i_write_field_4d(handle, field_md, MOM_domain, field, tstamp, &
                              tile_count, fill_value)
    import :: MOM_file, MOM_field, MOM_domain_type
    class(MOM_file), intent(inout) :: handle
      !< Handle for a file that is open for writing
    type(MOM_field), intent(in) :: field_md
      !< Field type with metadata
    type(MOM_domain_type), intent(in) :: MOM_domain
      !< The MOM_Domain that describes the decomposition
    real, intent(inout) :: field(:,:,:,:)
      !< Field to write
    real, optional, intent(in) :: tstamp
      !< Model time of this field
    integer, optional, intent(in) :: tile_count
      !< PEs per tile (default: 1)
    real, optional, intent(in) :: fill_value
      !< Missing data fill value
  end subroutine i_write_field_4d


  !> Interface to write_field_3d()
  subroutine i_write_field_3d(handle, field_md, MOM_domain, field, tstamp, &
                              tile_count, fill_value)
    import :: MOM_file, MOM_field, MOM_domain_type
    class(MOM_file), intent(inout) :: handle
      !< Handle for a file that is open for writing
    type(MOM_field), intent(in) :: field_md
      !< Field type with metadata
    type(MOM_domain_type), intent(in) :: MOM_domain
      !< The MOM_Domain that describes the decomposition
    real, intent(inout) :: field(:,:,:)
      !< Field to write
    real, optional, intent(in) :: tstamp
      !< Model time of this field
    integer, optional, intent(in) :: tile_count
      !< PEs per tile (default: 1)
    real, optional, intent(in) :: fill_value
      !< Missing data fill value
  end subroutine i_write_field_3d


  !> Interface to write_field_2d()
  subroutine i_write_field_2d(handle, field_md, MOM_domain, field, tstamp, &
                              tile_count, fill_value)
    import :: MOM_file, MOM_field, MOM_domain_type
    class(MOM_file), intent(inout) :: handle
      !< Handle for a file that is open for writing
    type(MOM_field), intent(in) :: field_md
      !< Field type with metadata
    type(MOM_domain_type), intent(in) :: MOM_domain
      !< The MOM_Domain that describes the decomposition
    real, dimension(:,:), intent(inout) :: field
      !< Field to write
    real, optional, intent(in) :: tstamp
      !< Model time of this field
    integer, optional, intent(in) :: tile_count
      !< PEs per tile (default: 1)
    real, optional, intent(in) :: fill_value
      !< Missing data fill value
  end subroutine i_write_field_2d


  !> Interface to write_field_1d()
  subroutine i_write_field_1d(handle, field_md, field, tstamp)
    import :: MOM_file, MOM_field
    class(MOM_file), intent(inout) :: handle
      !< Handle for a file that is open for writing
    type(MOM_field), intent(in) :: field_md
      !< Field type with metadata
    real, dimension(:), intent(in) :: field
      !< Field to write
    real, optional, intent(in) :: tstamp
      !< Model time of this field
  end subroutine i_write_field_1d


  !> Interface to write_field_0d()
  subroutine i_write_field_0d(handle, field_md, field, tstamp)
    import :: MOM_file, MOM_field
    class(MOM_file), intent(inout) :: handle
      !< Handle for a file that is open for writing
    type(MOM_field), intent(in) :: field_md
      !< Field type with metadata
    real, intent(in) :: field
      !< Field to write
    real, optional, intent(in) :: tstamp
      !< Model time of this field
  end subroutine i_write_field_0d


  !> Interface to write_field_axis()
  subroutine i_write_field_axis(handle, axis)
    import :: MOM_file, MOM_axis
    class(MOM_file), intent(inout) :: handle
      !< Handle for a file that is open for writing
    type(MOM_axis), intent(in) :: axis
      !< An axis type variable with information to write
  end subroutine i_write_field_axis


  !> Interface to file_is_open()
  logical function i_file_is_open(handle)
    import :: MOM_file
    class(MOM_file), intent(in) :: handle
      !< Handle to a file to inquire about
  end function i_file_is_open


  !> Interface to get_file_info()
  subroutine i_get_file_info(handle, ndim, nvar, ntime)
    import :: MOM_file
    class(MOM_file), intent(in) :: handle
      !< Handle for a file that is open for I/O
    integer, optional, intent(out) :: ndim
      !< The number of dimensions in the file
    integer, optional, intent(out) :: nvar
      !< The number of variables in the file
    integer, optional, intent(out) :: ntime
      !< The number of time levels in the file
  end subroutine i_get_file_info


  !> Interface to get_file_fields()
  subroutine i_get_file_fields(handle, fields)
    import :: MOM_file, MOM_field
    class(MOM_file), intent(inout) :: handle
      !< Handle for a file that is open for I/O
    type(MOM_field), dimension(:), intent(inout) :: fields
      !< Field-type descriptions of all of the variables in a file.
  end subroutine i_get_file_fields


  !> Interface to get_field_atts()
  subroutine i_get_field_atts(handle, field, name, units, longname, checksum)
    import :: MOM_file, MOM_field, int64
    class(MOM_file), intent(in) :: handle
      !< File where field is stored
    type(MOM_field), intent(in) :: field
      !< The field to extract information from
    character(len=*), optional, intent(out) :: name
      !< The variable name
    character(len=*), optional, intent(out) :: units
      !< The units of the variable
    character(len=*), optional, intent(out) :: longname
      !< The long name of the variable
    integer(kind=int64), optional, intent(out) :: checksum(:)
      !< The checksums of the variable in a file
  end subroutine i_get_field_atts


  !> Interface to read_field_chksum
  subroutine i_read_field_chksum(handle, field, chksum, valid_chksum)
    import :: MOM_file, MOM_field, int64
    class(MOM_file), intent(in) :: handle
      !< File where field is stored
    type(MOM_field), intent(in) :: field
      !< The field whose checksum attribute is to be read
    integer(kind=int64), intent(out) :: chksum
      !< The checksum for the field.
    logical, intent(out) :: valid_chksum
      !< If true, chksum has been successfully read
  end subroutine i_read_field_chksum
end interface

contains

!> Initialize the linked list of framework axes
subroutine initialize_axis_list_infra(list)
  class(axis_list_infra), intent(inout) :: list

  ! Pre-allocate the first node and set the tail to this empty node
  allocate(list%head)
  list%tail => list%head
end subroutine initialize_axis_list_infra


!> Append a new axis to the list
subroutine append_axis_list_infra(list, axis, label)
  class(axis_list_infra), intent(inout) :: list
  type(axistype), intent(in) :: axis
  character(len=*), intent(in) :: label

  type(axis_node_infra), pointer :: empty_node

  ! Transfer value to tail
  list%tail%label = label
  list%tail%axis = axis

  ! Extend list to next empty node
  allocate(empty_node)
  list%tail%next => empty_node
  list%tail => empty_node
end subroutine append_axis_list_infra


!> Get axis based on label
function get_axis_list_infra(list, label) result(axis)
  class(axis_list_infra), intent(in) :: list
  character(len=*), intent(in) :: label
  type(axistype) :: axis

  type(axis_node_infra), pointer :: node

  ! NOTE: The tail is a pre-allocated empty node, so we check node%next
  node => list%head
  do while(associated(node%next))
    if (node%label == label) exit
    node => node%next
  enddo
  if (.not. associated(node)) &
    call MOM_error(FATAL, "axis associated with " // label // " not found.")

  axis = node%axis
end function get_axis_list_infra


!> Deallocate axes of list
subroutine finalize_axis_list_infra(list)
  class(axis_list_infra), intent(inout) :: list

  type(axis_node_infra), pointer :: node, next_node

  node => list%head
  do while(associated(node))
    next_node => node
    node => node%next
    deallocate(next_node)
  enddo
end subroutine finalize_axis_list_infra


!> Initialize the linked list of framework axes
subroutine initialize_axis_list_nc(list)
  class(axis_list_nc), intent(inout) :: list

  ! Pre-allocate the first node and set the tail to this empty node
  allocate(list%head)
  list%tail => list%head
end subroutine initialize_axis_list_nc


!> Append a new axis to the list
subroutine append_axis_list_nc(list, axis, label)
  class(axis_list_nc), intent(inout) :: list
  type(netcdf_axis), intent(in) :: axis
  character(len=*), intent(in) :: label

  type(axis_node_nc), pointer :: empty_node

  ! Transfer value to tail
  list%tail%label = label
  list%tail%axis = axis

  ! Extend list to next empty node
  allocate(empty_node)
  list%tail%next => empty_node
  list%tail => empty_node
end subroutine append_axis_list_nc


!> Get axis based on label
function get_axis_list_nc(list, label) result(axis)
  class(axis_list_nc), intent(in) :: list
  character(len=*), intent(in) :: label
  type(netcdf_axis) :: axis

  type(axis_node_nc), pointer :: node

  ! NOTE: The tail is a pre-allocated empty node, so we check node%next
  node => list%head
  do while(associated(node%next))
    if (node%label == label) exit
    node => node%next
  enddo
  if (.not. associated(node)) &
    call MOM_error(FATAL, "axis associated with " // label // " not found.")

  axis = node%axis
end function get_axis_list_nc


!> Deallocate axes of list
subroutine finalize_axis_list_nc(list)
  class(axis_list_nc), intent(inout) :: list

  type(axis_node_nc), pointer :: node, next_node

  node => list%head
  do while(associated(node))
    next_node => node
    node => node%next
    deallocate(next_node)
  enddo
end subroutine finalize_axis_list_nc


!> Initialize the linked list of framework axes
subroutine initialize_field_list_infra(list)
  class(field_list_infra), intent(inout) :: list

  ! Pre-allocate the first node and set the tail to this empty node
  allocate(list%head)
  list%tail => list%head
end subroutine initialize_field_list_infra


!> Append a new field to the list
subroutine append_field_list_infra(list, field, label)
  class(field_list_infra), intent(inout) :: list
  type(fieldtype), intent(in) :: field
  character(len=*), intent(in) :: label

  type(field_node_infra), pointer :: empty_node

  ! Transfer value to tail
  list%tail%label = label
  list%tail%field = field

  ! Extend list to next empty node
  allocate(empty_node)
  list%tail%next => empty_node
  list%tail => empty_node
end subroutine append_field_list_infra


!> Get axis based on label
function get_field_list_infra(list, label) result(field)
  class(field_list_infra), intent(in) :: list
  character(len=*), intent(in) :: label
  type(fieldtype) :: field

  type(field_node_infra), pointer :: node

  ! NOTE: The tail is a pre-allocated empty node, so we check node%next
  node => list%head
  do while(associated(node%next))
    if (node%label == label) exit
    node => node%next
  enddo
  if (.not. associated(node)) &
    call MOM_error(FATAL, "field associated with " // label // " not found.")

  field = node%field
end function get_field_list_infra


!> Deallocate fields of list
subroutine finalize_field_list_infra(list)
  class(field_list_infra), intent(inout) :: list

  type(field_node_infra), pointer :: node, next_node

  node => list%head
  do while(associated(node))
    next_node => node
    node => node%next
    deallocate(next_node)
  enddo
end subroutine finalize_field_list_infra


!> Initialize the linked list of framework axes
subroutine initialize_field_list_nc(list)
  class(field_list_nc), intent(inout) :: list

  ! Pre-allocate the first node and set the tail to this empty node
  allocate(list%head)
  list%tail => list%head
end subroutine initialize_field_list_nc


!> Append a new field to the list
subroutine append_field_list_nc(list, field, label)
  class(field_list_nc), intent(inout) :: list
  type(netcdf_field), intent(in) :: field
  character(len=*), intent(in) :: label

  type(field_node_nc), pointer :: empty_node

  ! Transfer value to tail
  list%tail%label = label
  list%tail%field = field

  ! Extend list to next empty node
  allocate(empty_node)
  list%tail%next => empty_node
  list%tail => empty_node
end subroutine append_field_list_nc


!> Get axis based on label
function get_field_list_nc(list, label) result(field)
  class(field_list_nc), intent(in) :: list
  character(len=*), intent(in) :: label
  type(netcdf_field) :: field

  type(field_node_nc), pointer :: node

  ! NOTE: The tail is a pre-allocated empty node, so we check node%next
  node => list%head
  do while(associated(node%next))
    if (node%label == label) exit
    node => node%next
  enddo
  if (.not. associated(node)) &
    call MOM_error(FATAL, "field associated with " // label // " not found.")

  field = node%field
end function get_field_list_nc


!> Deallocate fields of list
subroutine finalize_field_list_nc(list)
  class(field_list_nc), intent(inout) :: list

  type(field_node_nc), pointer :: node, next_node

  node => list%head
  do while(associated(node))
    next_node => node
    node => node%next
    deallocate(next_node)
  enddo
end subroutine finalize_field_list_nc


!> Open a MOM framework file
subroutine open_file_infra(handle, filename, action, MOM_domain, threading, fileset)
  class(MOM_infra_file), intent(inout) :: handle
  character(len=*), intent(in) :: filename
  integer, intent(in), optional :: action
  type(MOM_domain_type), optional, intent(in) :: MOM_domain
  integer, intent(in), optional :: threading
  integer, intent(in), optional :: fileset

  logical :: use_single_file_domain
    ! True if the domain is replaced with a single-file IO layout.

  use_single_file_domain = .false.
  if (present(MOM_domain) .and. present(fileset)) then
    if (fileset == SINGLE_FILE) &
      use_single_file_domain = .true.
  endif

  if (use_single_file_domain) then
    call clone_MOM_domain(MOM_domain, handle%domain, io_layout=[1,1])
    call open_file(handle%handle_infra, filename, action=action, &
        MOM_domain=handle%domain, threading=threading, fileset=fileset)
  else
    call open_file(handle%handle_infra, filename, action=action, &
        MOM_domain=MOM_domain, threading=threading, fileset=fileset)
  endif

  call handle%axes%init()
  call handle%fields%init()
end subroutine open_file_infra

!> Close a MOM framework file
subroutine close_file_infra(handle)
  class(MOM_infra_file), intent(inout) :: handle

  if (associated(handle%domain)) &
    call deallocate_MOM_domain(handle%domain)

  call close_file(handle%handle_infra)
  call handle%axes%finalize()
  call handle%fields%finalize()
end subroutine close_file_infra

!> Flush the buffer of a MOM framework file
subroutine flush_file_infra(handle)
  class(MOM_infra_file), intent(in) :: handle

  call flush_file(handle%handle_infra)
end subroutine flush_file_infra


!> Register an axis to the MOM framework file
function register_axis_infra(handle, label, units, longname, &
    cartesian, sense, domain, data, edge_axis, calendar) result(axis)

  class(MOM_infra_file), intent(inout) :: handle
    !< Handle for a file that is open for writing
  character(len=*), intent(in) :: label
    !< The name in the file of this axis
  character(len=*), intent(in) :: units
    !< The units of this axis
  character(len=*), intent(in) :: longname
    !< The long description of this axis
  character(len=*), optional, intent(in) :: cartesian
    !< A variable indicating which direction this axis corresponds with.
    !! Valid values include 'X', 'Y', 'Z', 'T', and 'N' for none.
  integer, optional, intent(in) :: sense
    !< This is 1 for axes whose values increase upward, or -1 if they increase
    !! downward.
  type(domain1D), optional, intent(in) :: domain
    !< The domain decomposion for this axis
  real, dimension(:), optional, intent(in) :: data
    !< The coordinate values of the points on this axis
  logical, optional, intent(in) :: edge_axis
    !< If true, this axis marks an edge of the tracer cells
  character(len=*), optional, intent(in) :: calendar
    !< The name of the calendar used with a time axis
  type(MOM_axis) :: axis
    !< The axis type where this information is stored

  type(axistype) :: ax_infra

  ! Create new infra axis and assign to pre-allocated tail of axes
  call write_metadata(handle%handle_infra, ax_infra, label, units, longname, &
      cartesian=cartesian, sense=sense, domain=domain, data=data, &
      edge_axis=edge_axis, calendar=calendar)

  call handle%axes%append(ax_infra, label)
  axis%label = label
end function register_axis_infra


!> Register a field to the MOM framework file
function register_field_infra(handle, axes, label, units, longname, pack, &
    standard_name, checksum, conversion) result(field)
  class(MOM_infra_file), intent(inout) :: handle
    !< Handle for a file that is open for writing
  type(MOM_axis), dimension(:), intent(in) :: axes
    !< Handles for the axis used for this variable
  character(len=*), intent(in) :: label
    !< The name in the file of this variable
  character(len=*), intent(in) :: units
    !< The units of this variable
  character(len=*), intent(in) :: longname
    !< The long description of this variable
  integer, optional, intent(in) :: pack
    !< A precision reduction factor with which the variable.  The default, 1,
    !! has no reduction, but 2 is not uncommon.
  character(len=*), optional, intent(in) :: standard_name
    !< The standard (e.g., CMOR) name for this variable
  integer(kind=int64), dimension(:), optional, intent(in) :: checksum
    !< Checksum values that can be used to verify reads.
  real, optional, intent(in) :: conversion
    !< A factor to use to rescale the field before output [a A-1 ~> 1]
  type(MOM_field) :: field
    !< The field type where this information is stored

  type(fieldtype) :: field_infra
  type(axistype), allocatable :: field_axes(:)
  integer :: i

  ! Construct array of framework axes
  allocate(field_axes(size(axes)))
  do i = 1, size(axes)
    field_axes(i) = handle%axes%get(axes(i)%label)
  enddo

  call write_metadata(handle%handle_infra, field_infra, field_axes, label, &
      units, longname, pack=pack, standard_name=standard_name, checksum=checksum)

  call handle%fields%append(field_infra, label)
  field%label = label
  field%conversion = 1.0 ; if (present(conversion)) field%conversion = conversion
end function register_field_infra


!> Write a 4D field to the MOM framework file
subroutine write_field_4d_infra(handle, field_md, MOM_domain, field, tstamp, &
                                tile_count, fill_value)
  class(MOM_infra_file), intent(inout) :: handle
    !< Handle for a file that is open for writing
  type(MOM_field), intent(in) :: field_md
    !< Field type with metadata
  type(MOM_domain_type), intent(in) :: MOM_domain
    !< The MOM_Domain that describes the decomposition
  real, intent(inout) :: field(:,:,:,:)
    !< Field to write
  real, optional, intent(in) :: tstamp
    !< Model time of this field
  integer, optional, intent(in) :: tile_count
    !< PEs per tile (default: 1)
  real, optional, intent(in) :: fill_value
    !< Missing data fill value

  type(fieldtype) :: field_infra
  real, allocatable :: unscaled_field(:,:,:,:) ! An unscaled version of field for output [a]

  field_infra = handle%fields%get(field_md%label)
  if (field_md%conversion == 1.0) then
    call write_field(handle%handle_infra, field_infra, MOM_domain, field, &
        tstamp=tstamp, tile_count=tile_count, fill_value=fill_value)
  else
    allocate(unscaled_field, source=field)
    unscaled_field(:,:,:,:) = field_md%conversion * field(:,:,:,:)
    call write_field(handle%handle_infra, field_infra, MOM_domain, unscaled_field, &
        tstamp=tstamp, tile_count=tile_count, fill_value=fill_value)
    deallocate(unscaled_field)
  endif
end subroutine write_field_4d_infra


!> Write a 3D field to the MOM framework file
subroutine write_field_3d_infra(handle, field_md, MOM_domain, field, tstamp, &
                                tile_count, fill_value)
  class(MOM_infra_file), intent(inout) :: handle
    !< Handle for a file that is open for writing
  type(MOM_field), intent(in) :: field_md
    !< Field type with metadata
  type(MOM_domain_type), intent(in) :: MOM_domain
    !< The MOM_Domain that describes the decomposition
  real, intent(inout) :: field(:,:,:)
    !< Field to write, perhaps in arbitrary rescaled units [A ~> a]
  real, optional, intent(in) :: tstamp
    !< Model time of this field
  integer, optional, intent(in) :: tile_count
    !< PEs per tile (default: 1)
  real, optional, intent(in) :: fill_value
    !< Missing data fill value

  type(fieldtype) :: field_infra
  real, allocatable :: unscaled_field(:,:,:) ! An unscaled version of field for output [a]

  field_infra = handle%fields%get(field_md%label)
  if (field_md%conversion == 1.0) then
    call write_field(handle%handle_infra, field_infra, MOM_domain, field, &
        tstamp=tstamp, tile_count=tile_count, fill_value=fill_value)
  else
    allocate(unscaled_field, source=field)
    unscaled_field(:,:,:) = field_md%conversion * field(:,:,:)
    call write_field(handle%handle_infra, field_infra, MOM_domain, unscaled_field, &
        tstamp=tstamp, tile_count=tile_count, fill_value=fill_value)
    deallocate(unscaled_field)
  endif

end subroutine write_field_3d_infra


!> Write a 2D field to the MOM framework file
subroutine write_field_2d_infra(handle, field_md, MOM_domain, field, tstamp, &
                                tile_count, fill_value)
  class(MOM_infra_file), intent(inout) :: handle
    !< Handle for a file that is open for writing
  type(MOM_field), intent(in) :: field_md
    !< Field type with metadata
  type(MOM_domain_type), intent(in) :: MOM_domain
    !< The MOM_Domain that describes the decomposition
  real, dimension(:,:), intent(inout) :: field
    !< Field to write
  real, optional, intent(in) :: tstamp
    !< Model time of this field
  integer, optional, intent(in) :: tile_count
    !< PEs per tile (default: 1)
  real, optional, intent(in) :: fill_value
    !< Missing data fill value

  type(fieldtype) :: field_infra
  real, allocatable :: unscaled_field(:,:) ! An unscaled version of field for output [a]

  field_infra = handle%fields%get(field_md%label)
  if (field_md%conversion == 1.0) then
    call write_field(handle%handle_infra, field_infra, MOM_domain, field, &
        tstamp=tstamp, tile_count=tile_count, fill_value=fill_value)
  else
    allocate(unscaled_field, source=field)
    unscaled_field(:,:) = field_md%conversion * field(:,:)
    call write_field(handle%handle_infra, field_infra, MOM_domain, unscaled_field, &
        tstamp=tstamp, tile_count=tile_count, fill_value=fill_value)
    deallocate(unscaled_field)
  endif
end subroutine write_field_2d_infra


!> Write a 1D field to the MOM framework file
subroutine write_field_1d_infra(handle, field_md, field, tstamp)
  class(MOM_infra_file), intent(inout) :: handle
    !< Handle for a file that is open for writing
  type(MOM_field), intent(in) :: field_md
    !< Field type with metadata
  real, dimension(:), intent(in) :: field
    !< Field to write
  real, optional, intent(in) :: tstamp
    !< Model time of this field

  type(fieldtype) :: field_infra
  real, allocatable :: unscaled_field(:) ! An unscaled version of field for output [a]

  field_infra = handle%fields%get(field_md%label)
  if (field_md%conversion == 1.0) then
    call write_field(handle%handle_infra, field_infra, field, tstamp=tstamp)
  else
    allocate(unscaled_field, source=field)
    unscaled_field(:) = field_md%conversion * field(:)
    call write_field(handle%handle_infra, field_infra, unscaled_field, tstamp=tstamp)
    deallocate(unscaled_field)
  endif
end subroutine write_field_1d_infra


!> Write a 0D field to the MOM framework file
subroutine write_field_0d_infra(handle, field_md, field, tstamp)
  class(MOM_infra_file), intent(inout) :: handle
    !< Handle for a file that is open for writing
  type(MOM_field), intent(in) :: field_md
    !< Field type with metadata
  real, intent(in) :: field
    !< Field to write
  real, optional, intent(in) :: tstamp
    !< Model time of this field

  type(fieldtype) :: field_infra
  real :: unscaled_field ! An unscaled version of field for output [a]

  field_infra = handle%fields%get(field_md%label)
  unscaled_field = field_md%conversion*field
  call write_field(handle%handle_infra, field_infra, unscaled_field, tstamp=tstamp)
end subroutine write_field_0d_infra


!> Write an axis field to the MOM framework file
subroutine write_field_axis_infra(handle, axis)
  class(MOM_infra_file), intent(inout) :: handle
    !< Handle for a file that is open for writing
  type(MOM_axis), intent(in) :: axis
    !< An axis type variable with information to write

  type(axistype) :: axis_infra
    !< An axis type variable with information to write

  axis_infra = handle%axes%get(axis%label)
  call write_field(handle%handle_infra, axis_infra)
end subroutine write_field_axis_infra


!> Write global metadata to the MOM framework file
subroutine write_attribute_infra(handle, name, attribute)
  class(MOM_infra_file), intent(in) :: handle
    !< Handle for a file that is open for writing
  character(len=*), intent(in) :: name
    !< The name in the file of this global attribute
  character(len=*), intent(in) :: attribute
    !< The value of this attribute

  call write_metadata(handle%handle_infra, name, attribute)
end subroutine write_attribute_infra


!> True if the framework file has been opened
logical function file_is_open_infra(handle)
  class(MOM_infra_file), intent(in) :: handle
    !< Handle to a file to inquire about

  file_is_open_infra = fms2_file_is_open(handle%handle_infra)
end function file_is_open_infra


!> Return number of dimensions, variables, or time levels in a MOM infra file
subroutine get_file_info_infra(handle, ndim, nvar, ntime)
  class(MOM_infra_file), intent(in) :: handle
    !< Handle for a file that is open for I/O
  integer, optional, intent(out) :: ndim
    !< The number of dimensions in the file
  integer, optional, intent(out) :: nvar
    !< The number of variables in the file
  integer,  optional, intent(out) :: ntime
    !< The number of time levels in the file

  call get_file_info(handle%handle_infra, ndim, nvar, ntime)
end subroutine get_file_info_infra


!> Return the field metadata associated with a MOM framework file
subroutine get_file_fields_infra(handle, fields)
  class(MOM_infra_file), intent(inout) :: handle
    !< Handle for a file that is open for I/O
  type(MOM_field), intent(inout) :: fields(:)
    !< Field-type descriptions of all of the variables in a file.

  type(fieldtype), allocatable :: fields_infra(:)
  integer :: i
  character(len=64) :: label

  allocate(fields_infra(size(fields)))
  call get_file_fields(handle%handle_infra, fields_infra)

  do i = 1, size(fields)
    call get_field_atts(fields_infra(i), name=label)
    call handle%fields%append(fields_infra(i), trim(label))
    fields(i)%label = trim(label)
  enddo
end subroutine get_file_fields_infra


!> Get time levels of a MOM framework file
subroutine get_file_times_infra(handle, time_values, ntime)
  class(MOM_infra_file), intent(in) :: handle
    !< Handle for a file that is open for I/O
  real, allocatable, dimension(:), intent(inout) :: time_values
    !< The real times for the records in file.
  integer, optional, intent(out) :: ntime
    !< The number of time levels in the file

  call get_file_times(handle%handle_infra, time_values, ntime=ntime)
end subroutine get_file_times_infra


!> Get attributes from a field
subroutine get_field_atts_infra(handle, field, name, units, longname, checksum)
  class(MOM_infra_file), intent(in) :: handle
    !< File where field is stored
  type(MOM_field), intent(in) :: field
    !< The field to extract information from
  character(len=*), optional, intent(out) :: name
    !< The variable name
  character(len=*), optional, intent(out) :: units
    !< The units of the variable
  character(len=*), optional, intent(out) :: longname
    !< The long name of the variable
  integer(kind=int64), optional, intent(out) :: checksum(:)
    !< The checksums of the variable in a file

  type(fieldtype) :: field_infra

  field_infra = handle%fields%get(field%label)
  call get_field_atts(field_infra, name, units, longname, checksum)
end subroutine get_field_atts_infra


!> Interface to read_field_chksum
subroutine read_field_chksum_infra(handle, field, chksum, valid_chksum)
  class(MOM_infra_file), intent(in) :: handle
    !< File where field is stored
  type(MOM_field), intent(in) :: field
    !< The field whose checksum attribute is to be read
  integer(kind=int64), intent(out) :: chksum
    !< The checksum for the field.
  logical, intent(out) :: valid_chksum
    !< If true, chksum has been successfully read

  type(fieldtype) :: field_infra

  field_infra = handle%fields%get(field%label)
  call read_field_chksum(field_infra, chksum, valid_chksum)
end subroutine read_field_chksum_infra

!> Get the native (fieldtype) fields of a MOM framework file
subroutine get_file_fieldtypes(handle, fields)
  class(MOM_infra_file), intent(in) :: handle
  type(fieldtype), intent(out) :: fields(:)

  type(field_node_infra), pointer :: node
  integer :: i

  ! NOTE: The tail is a pre-allocated empty node, so we check node%next
  node => handle%fields%head
  do i = 1, size(fields)
    if (.not. associated(node%next)) &
      call MOM_error(FATAL, 'fields(:) size exceeds number of registered fields.')
    fields(i) = node%field
    node => node%next
  enddo
end subroutine get_file_fieldtypes


! MOM_netcdf_file methods

!> Open a MOM netCDF file
subroutine open_file_nc(handle, filename, action, MOM_domain, threading, fileset)
  class(MOM_netcdf_file), intent(inout) :: handle
  character(len=*), intent(in) :: filename
  integer, intent(in), optional :: action
  type(MOM_domain_type), optional, intent(in) :: MOM_domain
  integer, intent(in), optional :: threading
  integer, intent(in), optional :: fileset

  if (.not. present(MOM_domain) .and. .not. is_root_PE()) return

  call open_netcdf_file(handle%handle_nc, filename, action)
  handle%is_open = .true.

  if (present(MOM_domain)) then
    handle%domain_decomposed = .true.

    ! Input files use unrotated indexing.
    if (associated(MOM_domain%domain_in)) then
      call hor_index_init(MOM_domain%domain_in, handle%HI)
    else
      call hor_index_init(MOM_domain, handle%HI)
    endif
  endif

  call handle%axes%init()
  call handle%fields%init()
end subroutine open_file_nc


!> Close a MOM netCDF file
subroutine close_file_nc(handle)
  class(MOM_netcdf_file), intent(inout) :: handle

  if (.not. handle%domain_decomposed .and. .not. is_root_PE()) return

  handle%is_open = .false.
  call close_netcdf_file(handle%handle_nc)
end subroutine close_file_nc


!> Flush the buffer of a MOM netCDF file
subroutine flush_file_nc(handle)
  class(MOM_netcdf_file), intent(in) :: handle

  if (.not. is_root_PE()) return

  call flush_netcdf_file(handle%handle_nc)
end subroutine flush_file_nc


!> Register an axis to the MOM netcdf file
function register_axis_nc(handle, label, units, longname, cartesian, sense, &
    domain, data, edge_axis, calendar) result(axis)
  class(MOM_netcdf_file), intent(inout) :: handle
    !< Handle for a netCDF file that is open for writing
  character(len=*), intent(in) :: label
    !< The name in the file of this axis
  character(len=*), intent(in) :: units
    !< The units of this axis
  character(len=*), intent(in) :: longname
    !< The long description of this axis
  character(len=*), optional, intent(in) :: cartesian
    !< A variable indicating which direction this axis corresponds with.
    !! Valid values include 'X', 'Y', 'Z', 'T', and 'N' for none.
  integer, optional, intent(in) :: sense
    !< This is 1 for axes whose values increase upward, or -1 if they increase
    !! downward.
  type(domain1D), optional, intent(in) :: domain
    !< The domain decomposion for this axis
  real, dimension(:), optional, intent(in) :: data
    !< The coordinate values of the points on this axis
  logical, optional, intent(in) :: edge_axis
    !< If true, this axis marks an edge of the tracer cells
  character(len=*), optional, intent(in) :: calendar
    !< The name of the calendar used with a time axis
  type(MOM_axis) :: axis

  type(netcdf_axis) :: axis_nc

  if (is_root_PE()) then
    axis_nc = register_netcdf_axis(handle%handle_nc, label, units, longname, &
        data, cartesian, sense)

    call handle%axes%append(axis_nc, label)
  endif
  axis%label = label
end function register_axis_nc


!> Register a field to the MOM netcdf file
function register_field_nc(handle, axes, label, units, longname, pack, &
    standard_name, checksum, conversion) result(field)
  class(MOM_netcdf_file), intent(inout) :: handle
    !< Handle for a file that is open for writing
  type(MOM_axis), intent(in) :: axes(:)
    !< Handles for the axis used for this variable
  character(len=*), intent(in) :: label
    !< The name in the file of this variable
  character(len=*), intent(in) :: units
    !< The units of this variable
  character(len=*), intent(in) :: longname
    !< The long description of this variable
  integer, optional, intent(in) :: pack
    !< A precision reduction factor with which the variable.  The default, 1,
    !! has no reduction, but 2 is not uncommon.
  character(len=*), optional, intent(in) :: standard_name
    !< The standard (e.g., CMOR) name for this variable
  integer(kind=int64), dimension(:), optional, intent(in) :: checksum
    !< Checksum values that can be used to verify reads.
  real, optional, intent(in) :: conversion
    !< A factor to use to rescale the field before output [a A-1 ~> 1]
  type(MOM_field) :: field

  type(netcdf_field) :: field_nc
  type(netcdf_axis), allocatable :: axes_nc(:)
  integer :: i

  if (is_root_PE()) then
    allocate(axes_nc(size(axes)))
    do i = 1, size(axes)
      axes_nc(i) = handle%axes%get(axes(i)%label)
    enddo

    field_nc = register_netcdf_field(handle%handle_nc, label, axes_nc, longname, units)

    call handle%fields%append(field_nc, label)
  endif
  field%label = label
  field%conversion = 1.0 ; if (present(conversion)) field%conversion = conversion
end function register_field_nc


!> Write global metadata to the MOM netcdf file
subroutine write_attribute_nc(handle, name, attribute)
  class(MOM_netcdf_file), intent(in) :: handle
    !< Handle for a file that is open for writing
  character(len=*), intent(in) :: name
    !< The name in the file of this global attribute
  character(len=*), intent(in) :: attribute
    !< The value of this attribute

  if (.not. is_root_PE()) return

  call write_netcdf_attribute(handle%handle_nc, name, attribute)
end subroutine write_attribute_nc


!> Write a 4D field to the MOM netcdf file
subroutine write_field_4d_nc(handle, field_md, MOM_domain, field, tstamp, &
    tile_count, fill_value)
  class(MOM_netcdf_file), intent(inout) :: handle
    !< Handle for a file that is open for writing
  type(MOM_field), intent(in) :: field_md
    !< Field type with metadata
  type(MOM_domain_type), intent(in) :: MOM_domain
    !< The MOM_Domain that describes the decomposition
  real, intent(inout) :: field(:,:,:,:)
    !< Field to write
  real, optional, intent(in) :: tstamp
    !< Model time of this field
  integer, optional, intent(in) :: tile_count
    !< PEs per tile (default: 1)
  real, optional, intent(in) :: fill_value
    !< Missing data fill value

  type(netcdf_field) :: field_nc
  real, allocatable :: unscaled_field(:,:,:,:) ! An unscaled version of field for output [a]

  if (.not. is_root_PE()) return

  field_nc = handle%fields%get(field_md%label)
  if (field_md%conversion == 1.0) then
    call write_netcdf_field(handle%handle_nc, field_nc, field, time=tstamp)
  else
    allocate(unscaled_field, source=field)
    unscaled_field(:,:,:,:) = field_md%conversion * field(:,:,:,:)
    call write_netcdf_field(handle%handle_nc, field_nc, unscaled_field, time=tstamp)
    deallocate(unscaled_field)
  endif
end subroutine write_field_4d_nc


!> Write a 3D field to the MOM netcdf file
subroutine write_field_3d_nc(handle, field_md, MOM_domain, field, tstamp, &
    tile_count, fill_value)
  class(MOM_netcdf_file), intent(inout) :: handle
    !< Handle for a file that is open for writing
  type(MOM_field), intent(in) :: field_md
    !< Field type with metadata
  type(MOM_domain_type), intent(in) :: MOM_domain
    !< The MOM_Domain that describes the decomposition
  real, intent(inout) :: field(:,:,:)
    !< Field to write
  real, optional, intent(in) :: tstamp
    !< Model time of this field
  integer, optional, intent(in) :: tile_count
    !< PEs per tile (default: 1)
  real, optional, intent(in) :: fill_value
    !< Missing data fill value

  type(netcdf_field) :: field_nc
  real, allocatable :: unscaled_field(:,:,:) ! An unscaled version of field for output [a]

  if (.not. is_root_PE()) return

  field_nc = handle%fields%get(field_md%label)
  if (field_md%conversion == 1.0) then
    call write_netcdf_field(handle%handle_nc, field_nc, field, time=tstamp)
  else
    allocate(unscaled_field, source=field)
    unscaled_field(:,:,:) = field_md%conversion * field(:,:,:)
    call write_netcdf_field(handle%handle_nc, field_nc, unscaled_field, time=tstamp)
    deallocate(unscaled_field)
  endif
end subroutine write_field_3d_nc


!> Write a 2D field to the MOM netcdf file
subroutine write_field_2d_nc(handle, field_md, MOM_domain, field, tstamp, &
    tile_count, fill_value)
  class(MOM_netcdf_file), intent(inout) :: handle
    !< Handle for a file that is open for writing
  type(MOM_field), intent(in) :: field_md
    !< Field type with metadata
  type(MOM_domain_type), intent(in) :: MOM_domain
    !< The MOM_Domain that describes the decomposition
  real, dimension(:,:), intent(inout) :: field
    !< Field to write
  real, optional, intent(in) :: tstamp
    !< Model time of this field
  integer, optional, intent(in) :: tile_count
    !< PEs per tile (default: 1)
  real, optional, intent(in) :: fill_value
    !< Missing data fill value

  type(netcdf_field) :: field_nc
  real, allocatable :: unscaled_field(:,:) ! An unscaled version of field for output [a]

  if (.not. is_root_PE()) return

  field_nc = handle%fields%get(field_md%label)
  if (field_md%conversion == 1.0) then
    call write_netcdf_field(handle%handle_nc, field_nc, field, time=tstamp)
  else
    allocate(unscaled_field, source=field)
    unscaled_field(:,:) = field_md%conversion * field(:,:)
    call write_netcdf_field(handle%handle_nc, field_nc, unscaled_field, time=tstamp)
    deallocate(unscaled_field)
  endif
end subroutine write_field_2d_nc


!> Write a 1D field to the MOM netcdf file
subroutine write_field_1d_nc(handle, field_md, field, tstamp)
  class(MOM_netcdf_file), intent(inout) :: handle
    !< Handle for a file that is open for writing
  type(MOM_field), intent(in) :: field_md
    !< Field type with metadata
  real, dimension(:), intent(in) :: field
    !< Field to write
  real, optional, intent(in) :: tstamp
    !< Model time of this field

  type(netcdf_field) :: field_nc
  real, allocatable :: unscaled_field(:) ! An unscaled version of field for output [a]

  if (.not. is_root_PE()) return

  field_nc = handle%fields%get(field_md%label)
  if (field_md%conversion == 1.0) then
    call write_netcdf_field(handle%handle_nc, field_nc, field, time=tstamp)
  else
    allocate(unscaled_field, source=field)
    unscaled_field(:) = field_md%conversion * field(:)
    call write_netcdf_field(handle%handle_nc, field_nc, unscaled_field, time=tstamp)
    deallocate(unscaled_field)
  endif
end subroutine write_field_1d_nc


!> Write a 0D field to the MOM netcdf file
subroutine write_field_0d_nc(handle, field_md, field, tstamp)
  class(MOM_netcdf_file), intent(inout) :: handle
    !< Handle for a file that is open for writing
  type(MOM_field), intent(in) :: field_md
    !< Field type with metadata
  real, intent(in) :: field
    !< Field to write
  real, optional, intent(in) :: tstamp
    !< Model time of this field

  type(netcdf_field) :: field_nc
  real :: unscaled_field ! An unscaled version of field for output [a]

  if (.not. is_root_PE()) return

  field_nc = handle%fields%get(field_md%label)
  unscaled_field = field_md%conversion * field
  call write_netcdf_field(handle%handle_nc, field_nc, unscaled_field, time=tstamp)
end subroutine write_field_0d_nc


!> Write an axis field to the MOM netcdf file
subroutine write_field_axis_nc(handle, axis)
  class(MOM_netcdf_file), intent(inout) :: handle
    !< Handle for a file that is open for writing
  type(MOM_axis), intent(in) :: axis
    !< An axis type variable with information to write

  type(netcdf_axis) :: axis_nc

  if (.not. is_root_PE()) return

  axis_nc = handle%axes%get(axis%label)
  call write_netcdf_axis(handle%handle_nc, axis_nc)
end subroutine write_field_axis_nc


!> True if the framework file has been opened
logical function file_is_open_nc(handle)
  class(MOM_netcdf_file), intent(in) :: handle
    !< Handle to a file to inquire about

  file_is_open_nc = handle%is_open
end function file_is_open_nc


!> Return number of dimensions, variables, or time levels in a MOM netcdf file
subroutine get_file_info_nc(handle, ndim, nvar, ntime)
  class(MOM_netcdf_file), intent(in) :: handle
    !< Handle for a file that is open for I/O
  integer, optional, intent(out) :: ndim
    !< The number of dimensions in the file
  integer, optional, intent(out) :: nvar
    !< The number of variables in the file
  integer,  optional, intent(out) :: ntime
    !< The number of time levels in the file

  integer :: ndim_nc, nvar_nc

  if (.not. is_root_PE()) return

  call get_netcdf_size(handle%handle_nc, ndims=ndim_nc, nvars=nvar_nc, nsteps=ntime)

  ! MOM I/O follows legacy FMS behavior and excludes axes from field count
  if (present(ndim)) ndim = ndim_nc
  if (present(nvar)) nvar = nvar_nc - ndim_nc
end subroutine get_file_info_nc


!> Update the axes and fields descriptors of a MOM netCDF file
subroutine update_file_contents_nc(handle)
  class(MOM_netcdf_file), intent(inout) :: handle
    !< Handle for a file that is open for I/O

  type(netcdf_axis), allocatable :: axes_nc(:)
    ! netCDF axis descriptors
  type(netcdf_field), allocatable :: fields_nc(:)
    ! netCDF field descriptors
  integer :: i
    ! Index counter

  if (.not. handle%domain_decomposed .and. .not. is_root_PE()) return

  call get_netcdf_fields(handle%handle_nc, axes_nc, fields_nc)

  do i = 1, size(axes_nc)
    call handle%axes%append(axes_nc(i), axes_nc(i)%label)
  enddo

  do i = 1, size(fields_nc)
    call handle%fields%append(fields_nc(i), fields_nc(i)%label)
  enddo
end subroutine update_file_contents_nc


!> Return the field descriptors of a MOM netCDF file
subroutine get_file_fields_nc(handle, fields)
  class(MOM_netcdf_file), intent(inout) :: handle
    !< Handle for a file that is open for I/O
  type(MOM_field), intent(inout) :: fields(:)
    !< Field-type descriptions of all of the variables in a file.

  type(field_node_nc), pointer :: node => null()
    ! Current field list node
  integer :: n
    ! Field counter

  if (.not. is_root_PE()) return

  ! Generate the manifest of axes and fields
  call handle%update()

  n = 0
  node => handle%fields%head
  do while (associated(node%next))
    n = n + 1
    fields(n)%label = trim(node%label)
    node => node%next
  enddo
end subroutine get_file_fields_nc


!> Get attributes from a netCDF field
subroutine get_field_atts_nc(handle, field, name, units, longname, checksum)
  class(MOM_netcdf_file), intent(in) :: handle
    !< File where field is stored
  type(MOM_field), intent(in) :: field
    !< The field to extract information from
  character(len=*), optional, intent(out) :: name
    !< The variable name
  character(len=*), optional, intent(out) :: units
    !< The units of the variable
  character(len=*), optional, intent(out) :: longname
    !< The long name of the variable
  integer(kind=int64), optional, intent(out) :: checksum(:)
    !< The checksums of the variable in a file

  call MOM_error(FATAL, 'get_field_atts over netCDF is not yet implemented.')
end subroutine get_field_atts_nc


!> Interface to read_field_chksum
subroutine read_field_chksum_nc(handle, field, chksum, valid_chksum)
  class(MOM_netcdf_file), intent(in) :: handle
    !< File where field is stored
  type(MOM_field), intent(in) :: field
    !< The field whose checksum attribute is to be read
  integer(kind=int64), intent(out) :: chksum
    !< The checksum for the field.
  logical, intent(out) :: valid_chksum
    !< If true, chksum has been successfully read

  call MOM_error(FATAL, 'read_field_chksum over netCDF is not yet implemented.')
  chksum = -1_int64
  valid_chksum = .false.
end subroutine read_field_chksum_nc


!> Read the values of a netCDF field into an array that might have halos
subroutine get_field_nc(handle, label, values, rescale)
  class(MOM_netcdf_file), intent(in) :: handle
    !< Handle of netCDF file to be read
  character(len=*), intent(in) :: label
    !< Field variable name
  real, intent(inout) :: values(:,:)
    !< Field values read from the file.  It would be intent(out) but for the
    !! need to preserve any initialized values in the halo regions.
  real, optional, intent(in) :: rescale
    !< A multiplicative rescaling factor for the values that are read.
    !! Omitting this is the same as setting it to 1.

  logical :: data_domain
    ! True if values matches the data domain size
  logical :: compute_domain
    ! True if values matches the compute domain size
  type(netcdf_field) :: field_nc
    ! netCDF field associated with label
  integer :: isc, iec, jsc, jec
    ! Index bounds of compute domain
  integer :: isd, ied, jsd, jed
    ! Index bounds of data domain
  integer :: iscl, iecl, jscl, jecl
    ! Local 1-based index bounds of compute domain
  integer :: bounds(2,2)
    ! Index bounds of domain
  real, allocatable :: values_c(:,:)
    ! Field values on the compute domain, used for copying to a data domain

  isc = handle%HI%isc
  iec = handle%HI%iec
  jsc = handle%HI%jsc
  jec = handle%HI%jec

  isd = handle%HI%isd
  ied = handle%HI%ied
  jsd = handle%HI%jsd
  jed = handle%HI%jed

  data_domain = all(shape(values) == [ied-isd+1, jed-jsd+1])
  compute_domain = all(shape(values) == [iec-isc+1, jec-jsc+1])

  ! NOTE: Data on face and vertex points is not yet supported.  This is a
  ! temporary check to detect such cases, but may be removed in the future.
  if (.not. (compute_domain .or. data_domain)) &
    call MOM_error(FATAL, 'get_field_nc trying to read '//trim(label)//' from '//&
                   trim(get_netcdf_filename(handle%handle_nc))//&
                   ': Only compute and data domains are currently supported.')

  field_nc = handle%fields%get(label)

  if (data_domain) &
    allocate(values_c(1:iec-isc+1,1:jec-jsc+1))

  if (handle%domain_decomposed) then
    bounds(1,:) = [isc, jsc] + [handle%HI%idg_offset, handle%HI%jdg_offset]
    bounds(2,:) = [iec, jec] + [handle%HI%idg_offset, handle%HI%jdg_offset]
    if (data_domain) then
      call read_netcdf_field(handle%handle_nc, field_nc, values_c, bounds=bounds)
    else
      call read_netcdf_field(handle%handle_nc, field_nc, values, bounds=bounds)
    endif
  else
    if (data_domain) then
      call read_netcdf_field(handle%handle_nc, field_nc, values_c)
    else
      call read_netcdf_field(handle%handle_nc, field_nc, values)
    endif
  endif

  if (data_domain) then
    iscl = isc - isd + 1
    iecl = iec - isd + 1
    jscl = jsc - jsd + 1
    jecl = jec - jsd + 1

    values(iscl:iecl,jscl:jecl) = values_c(:,:)
  else
    iscl = 1
    iecl = iec - isc + 1
    jscl = 1
    jecl = jec - jsc + 1
  endif

  ! NOTE: It is more efficient to do the rescale in-place while copying
  ! values_c(:,:) to values(:,:).  But since rescale is only present for
  ! debugging, we can probably disregard this impact on performance.
  if (present(rescale)) then
    if (rescale /= 1.0) then
      values(iscl:iecl,jscl:jecl) = rescale * values(iscl:iecl,jscl:jecl)
    endif
  endif
end subroutine get_field_nc


!> \namespace MOM_IO_file
!!
!! This file defines the MOM_file classes used to inferface with the internal
!! IO handlers, such as the configured "infra" layer (FMS) or native netCDF.
!!
!! `MOM_file`: The generic class used to reference any file type
!!    Cannot be used in a variable declaration.
!!
!! `MOM_infra_file`: A file handler for use by the infra layer.  Currently this
!!    means an FMS file, such a restart or diagnostic output.
!!
!! `MOM_netcdf_file`: A netCDF file handler for MOM-specific I/O.  This may
!!    include operations outside the scope of FMS or other infra frameworks.

end module MOM_io_file
