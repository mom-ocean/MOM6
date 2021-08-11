!> This module contains a thin inteface to mpp and fms I/O code
module MOM_io_infra

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domain_infra,     only : MOM_domain_type, rescale_comp_data, AGRID, BGRID_NE, CGRID_NE
use MOM_domain_infra,     only : domain2d, domain1d, CENTER, CORNER, NORTH_FACE, EAST_FACE
use MOM_error_infra,      only : MOM_error=>MOM_err, NOTE, FATAL, WARNING

use fms_mod,              only : write_version_number, open_namelist_file, check_nml_error
use fms_io_mod,           only : file_exist, field_exist, field_size, read_data
use fms_io_mod,           only : fms_io_exit, get_filename_appendix
use mpp_io_mod,           only : mpp_open, mpp_close, mpp_flush
use mpp_io_mod,           only : mpp_write_meta, mpp_write, mpp_read
use mpp_io_mod,           only : mpp_get_atts, mpp_attribute_exist
use mpp_io_mod,           only : mpp_get_axes, axistype, mpp_get_axis_data
use mpp_io_mod,           only : mpp_get_fields, fieldtype
use mpp_io_mod,           only : mpp_get_info, mpp_get_times
use mpp_io_mod,           only : mpp_io_init
use mpp_mod,              only : stdout_if_root=>stdout
! These are encoding constants.
use mpp_io_mod,           only : APPEND_FILE=>MPP_APPEND, WRITEONLY_FILE=>MPP_WRONLY
use mpp_io_mod,           only : OVERWRITE_FILE=>MPP_OVERWR, READONLY_FILE=>MPP_RDONLY
use mpp_io_mod,           only : NETCDF_FILE=>MPP_NETCDF, ASCII_FILE=>MPP_ASCII
use mpp_io_mod,           only : MULTIPLE=>MPP_MULTI, SINGLE_FILE=>MPP_SINGLE
use mpp_mod,              only : lowercase
use iso_fortran_env,      only : int64

implicit none ; private

! These interfaces are actually implemented or have explicit interfaces in this file.
public :: open_file, open_ASCII_file, file_is_open, close_file, flush_file, file_exists
public :: get_file_info, get_file_fields, get_file_times, get_filename_suffix
public :: MOM_read_data, MOM_read_vector, write_metadata, write_field
public :: field_exists, get_field_atts, get_field_size, get_axis_data, read_field_chksum
public :: io_infra_init, io_infra_end, MOM_namelist_file, check_namelist_error, write_version
public :: stdout_if_root
! These types are inherited from underlying infrastructure code, to act as containers for
! information about fields and axes, respectively, and are opaque to this module.
public :: fieldtype, axistype
! These are encoding constant parmeters.
public :: ASCII_FILE, NETCDF_FILE, SINGLE_FILE, MULTIPLE
public :: APPEND_FILE, READONLY_FILE, OVERWRITE_FILE, WRITEONLY_FILE
public :: CENTER, CORNER, NORTH_FACE, EAST_FACE

!> Indicate whether a file exists, perhaps with domain decomposition
interface file_exists
  module procedure FMS_file_exists
  module procedure MOM_file_exists
end interface

!> Open a file (or fileset) for parallel or single-file I/).
interface open_file
  module procedure open_file_type, open_file_unit
end interface open_file

!> Read a data field from a file
interface MOM_read_data
  module procedure MOM_read_data_4d
  module procedure MOM_read_data_3d
  module procedure MOM_read_data_2d, MOM_read_data_2d_region
  module procedure MOM_read_data_1d, MOM_read_data_1d_int
  module procedure MOM_read_data_0d, MOM_read_data_0d_int
end interface

!> Write a registered field to an output file
interface write_field
  module procedure write_field_4d
  module procedure write_field_3d
  module procedure write_field_2d
  module procedure write_field_1d
  module procedure write_field_0d
  module procedure MOM_write_axis
end interface write_field

!> Read a pair of data fields representing the two components of a vector from a file
interface MOM_read_vector
  module procedure MOM_read_vector_3d
  module procedure MOM_read_vector_2d
end interface MOM_read_vector

!> Write metadata about a variable or axis to a file and store it for later reuse
interface write_metadata
  module procedure write_metadata_axis, write_metadata_field, write_metadata_global
end interface write_metadata

!> Close a file (or fileset).  If the file handle does not point to an open file,
!! close_file simply returns without doing anything.
interface close_file
  module procedure close_file_type, close_file_unit
end interface close_file

!> Ensure that the output stream associated with a file handle is fully sent to disk
interface flush_file
  module procedure flush_file_type, flush_file_unit
end interface flush_file

!> Type for holding a handle to an open file and related information
type, public :: file_type ; private
  integer :: unit = -1 !< The framework identfier or netCDF unit number of an output file
  character(len=:), allocatable :: filename !< The path to this file, if it is open
  logical :: open_to_read  = .false. !< If true, this file or fileset can be read
  logical :: open_to_write = .false. !< If true, this file or fileset can be written to
end type file_type

contains

!> Reads the checksum value for a field that was recorded in a file, along with a flag indicating
!! whether the file contained a valid checksum for this field.
subroutine read_field_chksum(field, chksum, valid_chksum)
  type(fieldtype),     intent(in)  :: field !< The field whose checksum attribute is to be read.
  integer(kind=int64), intent(out) :: chksum !< The checksum for the field.
  logical,             intent(out) :: valid_chksum  !< If true, chksum has been successfully read.
  ! Local variables
  integer(kind=int64), dimension(3) :: checksum_file

  checksum_file(:) = -1
  valid_chksum = mpp_attribute_exist(field, "checksum")
  if (valid_chksum) then
    call get_field_atts(field, checksum=checksum_file)
    chksum = checksum_file(1)
  else
    chksum = -1
  endif
end subroutine read_field_chksum

!> Returns true if the named file or its domain-decomposed variant exists.
logical function MOM_file_exists(filename, MOM_Domain)
  character(len=*),       intent(in) :: filename   !< The name of the file being inquired about
  type(MOM_domain_type),  intent(in) :: MOM_Domain !< The MOM_Domain that describes the decomposition

! This function uses the fms_io function file_exist to determine whether
! a named file (or its decomposed variant) exists.

  MOM_file_exists = file_exist(filename, MOM_Domain%mpp_domain)

end function MOM_file_exists

!> Returns true if the named file or its domain-decomposed variant exists.
logical function FMS_file_exists(filename, domain, no_domain)
  character(len=*),         intent(in) :: filename  !< The name of the file being inquired about
  type(domain2d), optional, intent(in) :: domain    !< The mpp domain2d that describes the decomposition
  logical,        optional, intent(in) :: no_domain !< This file does not use domain decomposition
! This function uses the fms_io function file_exist to determine whether
! a named file (or its decomposed variant) exists.

  FMS_file_exists = file_exist(filename, domain, no_domain)

end function FMS_file_exists

!> indicates whether an I/O handle is attached to an open file
logical function file_is_open(IO_handle)
  type(file_type), intent(in) :: IO_handle !< Handle to a file to inquire about

  file_is_open = (IO_handle%unit >= 0)
end function file_is_open

!> closes a file (or fileset).  If the file handle does not point to an open file,
!! close_file_type simply returns without doing anything.
subroutine close_file_type(IO_handle)
  type(file_type), intent(inout) :: IO_handle   !< The I/O handle for the file to be closed

  call mpp_close(IO_handle%unit)
  if (allocated(IO_handle%filename)) deallocate(IO_handle%filename)
  IO_handle%open_to_read = .false. ; IO_handle%open_to_write = .false.
end subroutine close_file_type

!> closes a file.  If the unit does not point to an open file,
!! close_file_unit simply returns without doing anything.
subroutine close_file_unit(unit)
  integer, intent(inout) :: unit   !< The I/O unit for the file to be closed

  call mpp_close(unit)
end subroutine close_file_unit

!> Ensure that the output stream associated with a file handle is fully sent to disk.
subroutine flush_file_type(file)
  type(file_type), intent(in) :: file    !< The I/O handle for the file to flush

  call mpp_flush(file%unit)
end subroutine flush_file_type

!> Ensure that the output stream associated with a unit is fully sent to disk.
subroutine flush_file_unit(unit)
  integer, intent(in) :: unit    !< The I/O unit for the file to flush

  call mpp_flush(unit)
end subroutine flush_file_unit

!> Initialize the underlying I/O infrastructure
subroutine io_infra_init(maxunits)
  integer,   optional, intent(in) :: maxunits !< An optional maximum number of file
                                              !! unit numbers that can be used.
  call mpp_io_init(maxunit=maxunits)
end subroutine io_infra_init

!> Gracefully close out and terminate the underlying I/O infrastructure
subroutine io_infra_end()
  call fms_io_exit()
end subroutine io_infra_end

!> Open a single namelist file that is potentially readable by all PEs.
function MOM_namelist_file(file) result(unit)
  character(len=*), optional, intent(in) :: file !< The file to open, by default "input.nml".
  integer                                :: unit !< The opened unit number of the namelist file
  unit = open_namelist_file(file)
end function MOM_namelist_file

!> Checks the iostat argument that is returned after reading a namelist variable and writes a
!! message if there is an error.
subroutine check_namelist_error(IOstat, nml_name)
  integer,          intent(in) :: IOstat   !< An I/O status field from a namelist read call
  character(len=*), intent(in) :: nml_name !< The name of the namelist
  integer :: ierr
  ierr = check_nml_error(IOstat, nml_name)
end subroutine check_namelist_error

!> Write a file version number to the log file or other output file
subroutine write_version(version, tag, unit)
  character(len=*),           intent(in) :: version !< A string that contains the routine name and version
  character(len=*), optional, intent(in) :: tag  !< A tag name to add to the message
  integer,          optional, intent(in) :: unit !< An alternate unit number for output

  call write_version_number(version, tag, unit)
end subroutine write_version

!> open_file opens a file for parallel or single-file I/O.
subroutine open_file_unit(unit, filename, action, form, threading, fileset, nohdrs, domain, MOM_domain)
  integer,                  intent(out) :: unit   !< The I/O unit for the opened file
  character(len=*),         intent(in)  :: filename !< The name of the file being opened
  integer,        optional, intent(in)  :: action !< A flag indicating whether the file can be read
                                                  !! or written to and how to handle existing files.
  integer,        optional, intent(in)  :: form   !< A flag indicating the format of a new file.  The
                                                  !! default is ASCII_FILE, but NETCDF_FILE is also common.
  integer,        optional, intent(in)  :: threading !< A flag indicating whether one (SINGLE_FILE)
                                                  !! or multiple PEs (MULTIPLE) participate in I/O.
                                                  !! With the default, the root PE does I/O.
  integer,        optional, intent(in)  :: fileset !< A flag indicating whether multiple PEs doing I/O due
                                                  !! to threading=MULTIPLE write to the same file (SINGLE_FILE)
                                                  !! or to one file per PE (MULTIPLE, the default).
  logical,        optional, intent(in)  :: nohdrs !< If nohdrs is .TRUE., headers are not written to
                                                  !! ASCII files.  The default is .false.
  type(domain2d), optional, intent(in)  :: domain !< A domain2d type that describes the decomposition
  type(MOM_domain_type), optional, intent(in) :: MOM_Domain !< A MOM_Domain that describes the decomposition

  if (present(MOM_Domain)) then
    call mpp_open(unit, filename, action=action, form=form, threading=threading, fileset=fileset, &
                  nohdrs=nohdrs, domain=MOM_Domain%mpp_domain)
  else
    call mpp_open(unit, filename, action=action, form=form, threading=threading, fileset=fileset, &
                  nohdrs=nohdrs, domain=domain)
  endif
end subroutine open_file_unit

!> open_file opens a file for parallel or single-file I/O.
subroutine open_file_type(IO_handle, filename, action, MOM_domain, threading, fileset)
  type(file_type),          intent(inout) :: IO_handle !< The handle for the opened file
  character(len=*),         intent(in)    :: filename !< The path name of the file being opened
  integer,        optional, intent(in)    :: action !< A flag indicating whether the file can be read
                                                    !! or written to and how to handle existing files.
                                                    !! The default is WRITE_ONLY.
  type(MOM_domain_type), &
                  optional, intent(in)    :: MOM_Domain !< A MOM_Domain that describes the decomposition
  integer,        optional, intent(in)    :: threading !< A flag indicating whether one (SINGLE_FILE)
                                                    !! or multiple PEs (MULTIPLE) participate in I/O.
                                                    !! With the default, the root PE does I/O.
  integer,        optional, intent(in)    :: fileset !< A flag indicating whether multiple PEs doing I/O due
                                                    !! to threading=MULTIPLE write to the same file (SINGLE_FILE)
                                                    !! or to one file per PE (MULTIPLE, the default).

  if (present(MOM_Domain)) then
    call mpp_open(IO_handle%unit, filename, action=action, form=NETCDF_FILE, threading=threading, &
                  fileset=fileset, domain=MOM_Domain%mpp_domain)
  else
    call mpp_open(IO_handle%unit, filename, action=action, form=NETCDF_FILE, threading=threading, &
                  fileset=fileset)
  endif
  IO_handle%filename = trim(filename)
  if (present(action)) then
    if (action == READONLY_FILE) then
      IO_handle%open_to_read = .true. ; IO_handle%open_to_write = .false.
    else
      IO_handle%open_to_read = .false. ; IO_handle%open_to_write = .true.
    endif
  else
    IO_handle%open_to_read = .false. ; IO_handle%open_to_write = .true.
  endif

end subroutine open_file_type

!> open_file opens an ascii file for parallel or single-file I/O using Fortran read and write calls.
subroutine open_ASCII_file(unit, file, action, threading, fileset)
  integer,                  intent(out) :: unit   !< The I/O unit for the opened file
  character(len=*),         intent(in)  :: file   !< The name of the file being opened
  integer,        optional, intent(in)  :: action !< A flag indicating whether the file can be read
                                                  !! or written to and how to handle existing files.
  integer,        optional, intent(in)  :: threading !< A flag indicating whether one (SINGLE_FILE)
                                                  !! or multiple PEs (MULTIPLE) participate in I/O.
                                                  !! With the default, the root PE does I/O.
  integer,        optional, intent(in)  :: fileset !< A flag indicating whether multiple PEs doing I/O due
                                                  !! to threading=MULTIPLE write to the same file (SINGLE_FILE)
                                                  !! or to one file per PE (MULTIPLE, the default).

  call mpp_open(unit, file, action=action, form=ASCII_FILE, threading=threading, fileset=fileset, &
                  nohdrs=.true.)

end subroutine open_ASCII_file


!> Provide a string to append to filenames, to differentiate ensemble members, for example.
subroutine get_filename_suffix(suffix)
  character(len=*), intent(out) :: suffix !< A string to append to filenames

  call get_filename_appendix(suffix)
end subroutine get_filename_suffix


!> Get information about the number of dimensions, variables and time levels
!! in the file associated with an open file unit
subroutine get_file_info(IO_handle, ndim, nvar, ntime)
  type(file_type),    intent(in)  :: IO_handle !< Handle for a file that is open for I/O
  integer,  optional, intent(out) :: ndim  !< The number of dimensions in the file
  integer,  optional, intent(out) :: nvar  !< The number of variables in the file
  integer,  optional, intent(out) :: ntime !< The number of time levels in the file

  ! Local variables
  integer :: ndims, nvars, natts, ntimes

  call mpp_get_info(IO_handle%unit, ndims, nvars, natts, ntimes )

  if (present(ndim)) ndim = ndims
  if (present(nvar)) nvar = nvars
  if (present(ntime)) ntime = ntimes

end subroutine get_file_info


!> Get the times of records from a file
 !### Modify this to also convert to time_type, using information about the dimensions?
subroutine get_file_times(IO_handle, time_values, ntime)
  type(file_type),                 intent(in)    :: IO_handle !< Handle for a file that is open for I/O
  real, allocatable, dimension(:), intent(inout) :: time_values !< The real times for the records in file.
  integer,               optional, intent(out)   :: ntime !< The number of time levels in the file

  integer :: ntimes

  if (allocated(time_values)) deallocate(time_values)
  call get_file_info(IO_handle, ntime=ntimes)
  if (present(ntime)) ntime = ntimes
  if (ntimes > 0) then
    allocate(time_values(ntimes))
    call mpp_get_times(IO_handle%unit, time_values)
  endif
end subroutine get_file_times

!> Set up the field information (e.g., names and metadata) for all of the variables in a file.  The
!! argument fields must be allocated with a size that matches the number of variables in a file.
subroutine get_file_fields(IO_handle, fields)
  type(file_type),               intent(in)    :: IO_handle !< Handle for a file that is open for I/O
  type(fieldtype), dimension(:), intent(inout) :: fields !< Field-type descriptions of all of
                                                         !! the variables in a file.
  call mpp_get_fields(IO_handle%unit, fields)
end subroutine get_file_fields

!> Extract information from a field type, as stored or as found in a file
subroutine get_field_atts(field, name, units, longname, checksum)
  type(fieldtype),            intent(in)  :: field !< The field to extract information from
  character(len=*), optional, intent(out) :: name  !< The variable name
  character(len=*), optional, intent(out) :: units !< The units of the variable
  character(len=*), optional, intent(out) :: longname  !< The long name of the variable
  integer(kind=int64),  dimension(:), &
                    optional, intent(out) :: checksum !< The checksums of the variable in a file
  call mpp_get_atts(field, name=name, units=units, longname=longname, checksum=checksum)
end subroutine get_field_atts

!> Field_exists returns true if the field indicated by field_name is present in the
!! file file_name.  If file_name does not exist, it returns false.
function field_exists(filename, field_name, domain, no_domain, MOM_domain)
  character(len=*),                 intent(in) :: filename   !< The name of the file being inquired about
  character(len=*),                 intent(in) :: field_name !< The name of the field being sought
  type(domain2d), target, optional, intent(in) :: domain     !< A domain2d type that describes the decomposition
  logical,                optional, intent(in) :: no_domain  !< This file does not use domain decomposition
  type(MOM_domain_type),  optional, intent(in) :: MOM_Domain !< A MOM_Domain that describes the decomposition
  logical                                      :: field_exists !< True if filename exists and field_name is in filename

  if (present(MOM_domain)) then
    field_exists = field_exist(filename, field_name, domain=MOM_domain%mpp_domain, no_domain=no_domain)
  else
    field_exists = field_exist(filename, field_name, domain=domain, no_domain=no_domain)
  endif

end function field_exists

!> Given filename and fieldname, this subroutine returns the size of the field in the file
subroutine get_field_size(filename, fieldname, sizes, field_found, no_domain)
  character(len=*),      intent(in)    :: filename  !< The name of the file to read
  character(len=*),      intent(in)    :: fieldname !< The name of the variable whose sizes are returned
  integer, dimension(:), intent(inout) :: sizes     !< The sizes of the variable in each dimension
  logical,     optional, intent(out)   :: field_found !< This indicates whether the field was found in
                                                    !! the input file.  Without this argument, there
                                                    !! is a fatal error if the field is not found.
  logical,     optional, intent(in)    :: no_domain !< If present and true, do not check for file
                                                    !! names with an appended tile number

  call field_size(filename, fieldname, sizes, field_found=field_found, no_domain=no_domain)

end subroutine get_field_size

!> Extracts and returns the axis data stored in an axistype.
subroutine get_axis_data( axis, dat )
  type(axistype),     intent(in)  :: axis !< An axis type
  real, dimension(:), intent(out) :: dat  !< The data in the axis variable

  call mpp_get_axis_data( axis, dat )
end subroutine get_axis_data

!> This routine uses the fms_io subroutine read_data to read a scalar named
!! "fieldname" from a single or domain-decomposed file "filename".
subroutine MOM_read_data_0d(filename, fieldname, data, timelevel, scale, MOM_Domain, &
                            global_file, file_may_be_4d)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real,                   intent(inout) :: data      !< The 1-dimensional array into which the data
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.
  type(MOM_domain_type), &
                optional, intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  logical,      optional, intent(in)    :: global_file !< If true, read from a single global file
  logical,      optional, intent(in)    :: file_may_be_4d !< If true, this file may have 4-d arrays,
                                                     !! in which case a more elaborate set of calls
                                                     !! is needed to read it due to FMS limitations.

  ! Local variables
  character(len=80)  :: varname             ! The name of a variable in the file
  type(fieldtype), allocatable :: fields(:) ! An array of types describing all the variables in the file
  logical :: use_fms_read_data, file_is_global
  integer :: n, unit, ndim, nvar, natt, ntime

  use_fms_read_data = .true. ; if (present(file_may_be_4d)) use_fms_read_data = .not.file_may_be_4d
  file_is_global = .true. ; if (present(global_file)) file_is_global = global_file

  if (.not.use_fms_read_data) then
    if (file_is_global) then
      call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                    threading=MULTIPLE, fileset=SINGLE_FILE) !, domain=MOM_Domain%mpp_domain )
    else
      call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                    threading=MULTIPLE, fileset=MULTIPLE, domain=MOM_Domain%mpp_domain )
    endif
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(fields(nvar))
    call mpp_get_fields(unit, fields(1:nvar))
    do n=1, nvar
      call mpp_get_atts(fields(n), name=varname)
      if (lowercase(trim(varname)) == lowercase(trim(fieldname))) then
        ! Maybe something should be done depending on the value of ntime.
        call mpp_read(unit, fields(n), data, timelevel)
        exit
      endif
    enddo

    deallocate(fields)
    call mpp_close(unit)
  elseif (present(MOM_Domain)) then
    call read_data(filename, fieldname, data, MOM_Domain%mpp_domain, timelevel=timelevel)
  else
    call read_data(filename, fieldname, data, timelevel=timelevel, no_domain=.true.)
  endif

  if (present(scale)) then ; if (scale /= 1.0) then
    data = scale*data
  endif ; endif

end subroutine MOM_read_data_0d

!> This routine uses the fms_io subroutine read_data to read a 1-D data field named
!! "fieldname" from a single or domain-decomposed file "filename".
subroutine MOM_read_data_1d(filename, fieldname, data, timelevel, scale, MOM_Domain, &
                            global_file, file_may_be_4d)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:),     intent(inout) :: data      !< The 1-dimensional array into which the data
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before they are returned.
  type(MOM_domain_type), &
                optional, intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  logical,      optional, intent(in)    :: global_file !< If true, read from a single global file
  logical,      optional, intent(in)    :: file_may_be_4d !< If true, this file may have 4-d arrays,
                                                     !! in which case a more elaborate set of calls
                                                     !! is needed to read it due to FMS limitations.

  ! Local variables
  character(len=80)  :: varname             ! The name of a variable in the file
  type(fieldtype), allocatable :: fields(:) ! An array of types describing all the variables in the file
  logical :: use_fms_read_data, file_is_global
  integer :: n, unit, ndim, nvar, natt, ntime

  use_fms_read_data = .true. ; if (present(file_may_be_4d)) use_fms_read_data = .not.file_may_be_4d
  file_is_global = .true. ; if (present(global_file)) file_is_global = global_file

  if (.not.use_fms_read_data) then
    if (file_is_global) then
      call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                    threading=MULTIPLE, fileset=SINGLE_FILE) !, domain=MOM_Domain%mpp_domain )
    else
      call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                    threading=MULTIPLE, fileset=MULTIPLE, domain=MOM_Domain%mpp_domain )
    endif
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(fields(nvar))
    call mpp_get_fields(unit, fields(1:nvar))
    do n=1, nvar
      call mpp_get_atts(fields(n), name=varname)
      if (lowercase(trim(varname)) == lowercase(trim(fieldname))) then
        call MOM_error(NOTE, "Reading 1-d variable "//trim(fieldname)//" from file "//trim(filename))
        ! Maybe something should be done depending on the value of ntime.
        call mpp_read(unit, fields(n), data, timelevel)
        exit
      endif
    enddo
    if ((n == nvar+1) .or. (nvar < 1)) call MOM_error(WARNING, &
      "MOM_read_data apparently did not find 1-d variable "//trim(fieldname)//" in file "//trim(filename))

    deallocate(fields)
    call mpp_close(unit)
  elseif (present(MOM_Domain)) then
    call read_data(filename, fieldname, data, MOM_Domain%mpp_domain, timelevel=timelevel)
  else
    call read_data(filename, fieldname, data, timelevel=timelevel, no_domain=.true.)
  endif

  if (present(scale)) then ; if (scale /= 1.0) then
    data(:) = scale*data(:)
  endif ; endif

end subroutine MOM_read_data_1d

!> This routine uses the fms_io subroutine read_data to read a distributed
!! 2-D data field named "fieldname" from file "filename".  Valid values for
!! "position" include CORNER, CENTER, EAST_FACE and NORTH_FACE.
subroutine MOM_read_data_2d(filename, fieldname, data, MOM_Domain, &
                            timelevel, position, scale, global_file, file_may_be_4d)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:),   intent(inout) :: data      !< The 2-dimensional array into which the data
                                                     !! should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: position  !< A flag indicating where this data is located
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.
  logical,      optional, intent(in)    :: global_file !< If true, read from a single global file
  logical,      optional, intent(in)    :: file_may_be_4d !< If true, this file may have 4-d arrays,
                                                     !! in which case a more elaborate set of calls
                                                     !! is needed to read it due to FMS limitations.

  ! Local variables
  character(len=80)  :: varname             ! The name of a variable in the file
  type(fieldtype), allocatable :: fields(:) ! An array of types describing all the variables in the file
  logical :: use_fms_read_data, file_is_global
  integer :: n, unit, ndim, nvar, natt, ntime

  use_fms_read_data = .true. ; if (present(file_may_be_4d)) use_fms_read_data = .not.file_may_be_4d
  file_is_global = .true. ; if (present(global_file)) file_is_global = global_file

  if (use_fms_read_data) then
    call read_data(filename, fieldname, data, MOM_Domain%mpp_domain, &
                   timelevel=timelevel, position=position)
  else
    if (file_is_global) then
      call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                    threading=MULTIPLE, fileset=SINGLE_FILE) !, domain=MOM_Domain%mpp_domain )
    else
      call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                    threading=MULTIPLE, fileset=MULTIPLE, domain=MOM_Domain%mpp_domain )
    endif
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(fields(nvar))
    call mpp_get_fields(unit, fields(1:nvar))
    do n=1, nvar
      call mpp_get_atts(fields(n), name=varname)
      if (lowercase(trim(varname)) == lowercase(trim(fieldname))) then
        call MOM_error(NOTE, "Reading 2-d variable "//trim(fieldname)//" from file "//trim(filename))
        ! Maybe something should be done depending on the value of ntime.
        call mpp_read(unit, fields(n), MOM_Domain%mpp_domain, data, timelevel)
        exit
      endif
    enddo
    if ((n == nvar+1) .or. (nvar < 1)) call MOM_error(WARNING, &
      "MOM_read_data apparently did not find 2-d variable "//trim(fieldname)//" in file "//trim(filename))

    deallocate(fields)
    call mpp_close(unit)
  endif

  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(MOM_Domain, data, scale)
  endif ; endif

end subroutine MOM_read_data_2d

!> This routine uses the fms_io subroutine read_data to read a region from a distributed or
!! global 2-D data field named "fieldname" from file "filename".
subroutine MOM_read_data_2d_region(filename, fieldname, data, start, nread, MOM_domain, &
                                   no_domain, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:),   intent(inout) :: data      !< The 2-dimensional array into which the data
                                                     !! should be read
  integer, dimension(:),  intent(in)    :: start     !< The starting index to read in each of 4
                                                     !! dimensions.  For this 2-d read, the 3rd
                                                     !! and 4th values are always 1.
  integer, dimension(:),  intent(in)    :: nread     !< The number of points to read in each of 4
                                                     !! dimensions.  For this 2-d read, the 3rd
                                                     !! and 4th values are always 1.
  type(MOM_domain_type), &
                optional, intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  logical,      optional, intent(in)    :: no_domain !< If present and true, this variable does not
                                                     !! use domain decomposion.
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.

  if (present(MOM_Domain)) then
    call read_data(filename, fieldname, data, start, nread, domain=MOM_Domain%mpp_domain, &
                   no_domain=no_domain)
  else
    call read_data(filename, fieldname, data, start, nread, no_domain=no_domain)
  endif

  if (present(scale)) then ; if (scale /= 1.0) then
    if (present(MOM_Domain)) then
      call rescale_comp_data(MOM_Domain, data, scale)
    else
      ! Dangerously rescale the whole array
      data(:,:) = scale*data(:,:)
    endif
  endif ; endif

end subroutine MOM_read_data_2d_region

!> This routine uses the fms_io subroutine read_data to read a distributed
!! 3-D data field named "fieldname" from file "filename".  Valid values for
!! "position" include CORNER, CENTER, EAST_FACE and NORTH_FACE.
subroutine MOM_read_data_3d(filename, fieldname, data, MOM_Domain, &
                            timelevel, position, scale, global_file, file_may_be_4d)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:), intent(inout) :: data      !< The 3-dimensional array into which the data
                                                     !! should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: position  !< A flag indicating where this data is located
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.
  logical,      optional, intent(in)    :: global_file !< If true, read from a single global file
  logical,      optional, intent(in)    :: file_may_be_4d !< If true, this file may have 4-d arrays,
                                                     !! in which case a more elaborate set of calls
                                                     !! is needed to read it due to FMS limitations.

  ! Local variables
  character(len=80)  :: varname             ! The name of a variable in the file
  type(fieldtype), allocatable :: fields(:) ! An array of types describing all the variables in the file
  logical :: use_fms_read_data, file_is_global
  integer :: n, unit, ndim, nvar, natt, ntime

  use_fms_read_data = .true. ; if (present(file_may_be_4d)) use_fms_read_data = .not.file_may_be_4d
  file_is_global = .true. ; if (present(global_file)) file_is_global = global_file

  if (use_fms_read_data) then
    call read_data(filename, fieldname, data, MOM_Domain%mpp_domain, &
                   timelevel=timelevel, position=position)
  else
    if (file_is_global) then
      call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                    threading=MULTIPLE, fileset=SINGLE_FILE) !, domain=MOM_Domain%mpp_domain )
    else
      call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                    threading=MULTIPLE, fileset=MULTIPLE, domain=MOM_Domain%mpp_domain )
    endif
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(fields(nvar))
    call mpp_get_fields(unit, fields(1:nvar))
    do n=1, nvar
      call mpp_get_atts(fields(n), name=varname)
      if (lowercase(trim(varname)) == lowercase(trim(fieldname))) then
        call MOM_error(NOTE, "Reading 3-d variable "//trim(fieldname)//" from file "//trim(filename))
        ! Maybe something should be done depending on the value of ntime.
        call mpp_read(unit, fields(n), MOM_Domain%mpp_domain, data, timelevel)
        exit
      endif
    enddo
    if ((n == nvar+1) .or. (nvar < 1)) call MOM_error(WARNING, &
      "MOM_read_data apparently did not find 3-d variable "//trim(fieldname)//" in file "//trim(filename))

    deallocate(fields)
    call mpp_close(unit)
  endif

  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(MOM_Domain, data, scale)
  endif ; endif

end subroutine MOM_read_data_3d

!> This routine uses the fms_io subroutine read_data to read a distributed
!! 4-D data field named "fieldname" from file "filename".  Valid values for
!! "position" include CORNER, CENTER, EAST_FACE and NORTH_FACE.
subroutine MOM_read_data_4d(filename, fieldname, data, MOM_Domain, &
                            timelevel, position, scale, global_file)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:,:), intent(inout) :: data    !< The 4-dimensional array into which the data
                                                     !! should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: position  !< A flag indicating where this data is located
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.
  logical,      optional, intent(in)    :: global_file !< If true, read from a single global file

  ! Local variables
  character(len=80)  :: varname             ! The name of a variable in the file
  type(fieldtype), allocatable :: fields(:) ! An array of types describing all the variables in the file
  logical :: use_fms_read_data, file_is_global
  integer :: n, unit, ndim, nvar, natt, ntime
  integer :: is, ie, js, je

  ! This single call does not work for a 4-d array due to FMS limitations, so multiple calls are
  ! needed.
  ! call read_data(filename, fieldname, data, MOM_Domain%mpp_domain, &
  !                timelevel=timelevel, position=position)

  file_is_global = .true. ; if (present(global_file)) file_is_global = global_file

  if (file_is_global) then
    call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                  threading=MULTIPLE, fileset=SINGLE_FILE) !, domain=MOM_Domain%mpp_domain )
  else
    call mpp_open(unit, trim(filename), form=NETCDF_FILE, action=READONLY_FILE, &
                  threading=MULTIPLE, fileset=MULTIPLE, domain=MOM_Domain%mpp_domain )
  endif
  call mpp_get_info(unit, ndim, nvar, natt, ntime)
  allocate(fields(nvar))
  call mpp_get_fields(unit, fields(1:nvar))
  do n=1, nvar
    call mpp_get_atts(fields(n), name=varname)
    if (lowercase(trim(varname)) == lowercase(trim(fieldname))) then
        call MOM_error(NOTE, "Reading 4-d variable "//trim(fieldname)//" from file "//trim(filename))
      ! Maybe something should be done depending on the value of ntime.
      call mpp_read(unit, fields(n), MOM_Domain%mpp_domain, data, timelevel)
      exit
    endif
  enddo
  if ((n == nvar+1) .or. (nvar < 1)) call MOM_error(WARNING, &
    "MOM_read_data apparently did not find 4-d variable "//trim(fieldname)//" in file "//trim(filename))

  deallocate(fields)
  call mpp_close(unit)

  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(MOM_Domain, data, scale)
  endif ; endif

end subroutine MOM_read_data_4d

!> This routine uses the fms_io subroutine read_data to read a scalar integer
!! data field named "fieldname" from file "filename".
subroutine MOM_read_data_0d_int(filename, fieldname, data, timelevel)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  integer,                intent(inout) :: data      !< The 1-dimensional array into which the data
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read

  call read_data(filename, fieldname, data, timelevel=timelevel, no_domain=.true.)

end subroutine MOM_read_data_0d_int

!> This routine uses the fms_io subroutine read_data to read a 1-D integer
!! data field named "fieldname" from file "filename".
subroutine MOM_read_data_1d_int(filename, fieldname, data, timelevel)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  integer, dimension(:),  intent(inout) :: data      !< The 1-dimensional array into which the data
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read

  call read_data(filename, fieldname, data, timelevel=timelevel, no_domain=.true.)

end subroutine MOM_read_data_1d_int


!> This routine uses the fms_io subroutine read_data to read a pair of distributed
!! 2-D data fields with names given by "[uv]_fieldname" from file "filename".  Valid values for
!! "stagger" include CGRID_NE, BGRID_NE, and AGRID.
subroutine MOM_read_vector_2d(filename, u_fieldname, v_fieldname, u_data, v_data, MOM_Domain, &
                              timelevel, stagger, scalar_pair, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: u_fieldname !< The variable name of the u data in the file
  character(len=*),       intent(in)    :: v_fieldname !< The variable name of the v data in the file
  real, dimension(:,:),   intent(inout) :: u_data    !< The 2 dimensional array into which the
                                                     !! u-component of the data should be read
  real, dimension(:,:),   intent(inout) :: v_data    !< The 2 dimensional array into which the
                                                     !! v-component of the data should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: stagger   !< A flag indicating where this vector is discretized
  logical,      optional, intent(in)    :: scalar_pair !< If true, a pair of scalars are to be read
  real,         optional, intent(in)    :: scale     !< A scaling factor that the fields are multiplied
                                                     !! by before they are returned.
  integer :: u_pos, v_pos

  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == BGRID_NE) then ; u_pos = CORNER ; v_pos = CORNER
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif

  call read_data(filename, u_fieldname, u_data, MOM_Domain%mpp_domain, &
                 timelevel=timelevel, position=u_pos)
  call read_data(filename, v_fieldname, v_data, MOM_Domain%mpp_domain, &
                 timelevel=timelevel, position=v_pos)

  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(MOM_Domain, u_data, scale)
    call rescale_comp_data(MOM_Domain, v_data, scale)
  endif ; endif

end subroutine MOM_read_vector_2d

!> This routine uses the fms_io subroutine read_data to read a pair of distributed
!! 3-D data fields with names given by "[uv]_fieldname" from file "filename".  Valid values for
!! "stagger" include CGRID_NE, BGRID_NE, and AGRID.
subroutine MOM_read_vector_3d(filename, u_fieldname, v_fieldname, u_data, v_data, MOM_Domain, &
                              timelevel, stagger, scalar_pair, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: u_fieldname !< The variable name of the u data in the file
  character(len=*),       intent(in)    :: v_fieldname !< The variable name of the v data in the file
  real, dimension(:,:,:), intent(inout) :: u_data    !< The 3 dimensional array into which the
                                                     !! u-component of the data should be read
  real, dimension(:,:,:), intent(inout) :: v_data    !< The 3 dimensional array into which the
                                                     !! v-component of the data should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: stagger   !< A flag indicating where this vector is discretized
  logical,      optional, intent(in)    :: scalar_pair !< If true, a pair of scalars are to be read.cretized
  real,         optional, intent(in)    :: scale     !< A scaling factor that the fields are multiplied
                                                     !! by before they are returned.

  integer :: u_pos, v_pos

  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == BGRID_NE) then ; u_pos = CORNER ; v_pos = CORNER
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif

  call read_data(filename, u_fieldname, u_data, MOM_Domain%mpp_domain, &
                 timelevel=timelevel, position=u_pos)
  call read_data(filename, v_fieldname, v_data, MOM_Domain%mpp_domain, &
                 timelevel=timelevel, position=v_pos)

  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(MOM_Domain, u_data, scale)
    call rescale_comp_data(MOM_Domain, v_data, scale)
  endif ; endif

end subroutine MOM_read_vector_3d


!> Write a 4d field to an output file.
subroutine write_field_4d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, fill_value)
  type(file_type),          intent(in)    :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),          intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),    intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:,:,:), intent(inout) :: field      !< Field to write
  real,           optional, intent(in)    :: tstamp     !< Model time of this field
  integer,        optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,           optional, intent(in)    :: fill_value !< Missing data fill value

  call mpp_write(IO_handle%unit, field_md, MOM_domain%mpp_domain, field, tstamp=tstamp, &
                 tile_count=tile_count, default_data=fill_value)
end subroutine write_field_4d

!> Write a 3d field to an output file.
subroutine write_field_3d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, fill_value)
  type(file_type),        intent(in)    :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),  intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:,:), intent(inout) :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model time of this field
  integer,      optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,         optional, intent(in)    :: fill_value !< Missing data fill value

  call mpp_write(IO_handle%unit, field_md, MOM_domain%mpp_domain, field, tstamp=tstamp, &
                   tile_count=tile_count, default_data=fill_value)
end subroutine write_field_3d

!> Write a 2d field to an output file.
subroutine write_field_2d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, fill_value)
  type(file_type),        intent(in)    :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),  intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:),   intent(inout) :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model time of this field
  integer,      optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,         optional, intent(in)    :: fill_value !< Missing data fill value

  call mpp_write(IO_handle%unit, field_md, MOM_domain%mpp_domain, field, tstamp=tstamp, &
                   tile_count=tile_count, default_data=fill_value)
end subroutine write_field_2d

!> Write a 1d field to an output file.
subroutine write_field_1d(IO_handle, field_md, field, tstamp)
  type(file_type),        intent(in)    :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  real, dimension(:),     intent(in)    :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model time of this field

  call mpp_write(IO_handle%unit, field_md, field, tstamp=tstamp)
end subroutine write_field_1d

!> Write a 0d field to an output file.
subroutine write_field_0d(IO_handle, field_md, field, tstamp)
  type(file_type),        intent(in)    :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  real,                   intent(in)    :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model time of this field

  call mpp_write(IO_handle%unit, field_md, field, tstamp=tstamp)
end subroutine write_field_0d

!> Write the data for an axis
subroutine MOM_write_axis(IO_handle, axis)
  type(file_type), intent(in) :: IO_handle  !< Handle for a file that is open for writing
  type(axistype),  intent(in) :: axis       !< An axis type variable with information to write

  call mpp_write(IO_handle%unit, axis)

end subroutine MOM_write_axis

!> Store information about an axis in a previously defined axistype and write this
!! information to the file indicated by unit.
subroutine write_metadata_axis(IO_handle, axis, name, units, longname, cartesian, sense, domain, &
                               data, edge_axis, calendar)
  type(file_type),            intent(in)    :: IO_handle  !< Handle for a file that is open for writing
  type(axistype),             intent(inout) :: axis  !< The axistype where this information is stored.
  character(len=*),           intent(in)    :: name  !< The name in the file of this axis
  character(len=*),           intent(in)    :: units !< The units of this axis
  character(len=*),           intent(in)    :: longname !< The long description of this axis
  character(len=*), optional, intent(in)    :: cartesian !< A variable indicating which direction
                                                     !! this axis corresponds with. Valid values
                                                     !! include 'X', 'Y', 'Z', 'T', and 'N' for none.
  integer,          optional, intent(in)    :: sense !< This is 1 for axes whose values increase upward, or
                                                     !! -1 if they increase downward.
  type(domain1D),   optional, intent(in)    :: domain !< The domain decomposion for this axis
  real, dimension(:), optional, intent(in)  :: data   !< The coordinate values of the points on this axis
  logical,          optional, intent(in)    :: edge_axis !< If true, this axis marks an edge of the tracer cells
  character(len=*), optional, intent(in)    :: calendar !< The name of the calendar used with a time axis

  call mpp_write_meta(IO_handle%unit, axis, name, units, longname, cartesian=cartesian, sense=sense, &
                      domain=domain, data=data, calendar=calendar)
end subroutine write_metadata_axis

!> Store information about an output variable in a previously defined fieldtype and write this
!! information to the file indicated by unit.
subroutine write_metadata_field(IO_handle, field, axes, name, units, longname, &
                                pack, standard_name, checksum)
  type(file_type),            intent(in)    :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),            intent(inout) :: field !< The fieldtype where this information is stored
  type(axistype), dimension(:), intent(in)  :: axes  !< Handles for the axis used for this variable
  character(len=*),           intent(in)    :: name  !< The name in the file of this variable
  character(len=*),           intent(in)    :: units !< The units of this variable
  character(len=*),           intent(in)    :: longname !< The long description of this variable
  integer,          optional, intent(in)    :: pack  !< A precision reduction factor with which the
                                                     !! variable.  The default, 1, has no reduction,
                                                     !! but 2 is not uncommon.
  character(len=*), optional, intent(in)    :: standard_name !< The standard (e.g., CMOR) name for this variable
  integer(kind=int64), dimension(:), &
                    optional, intent(in)    :: checksum !< Checksum values that can be used to verify reads.


  call mpp_write_meta(IO_handle%unit, field, axes, name, units, longname, &
                      pack=pack, standard_name=standard_name, checksum=checksum)
  ! unused opt. args: min=min, max=max, fill=fill, scale=scale, add=add, &

end subroutine write_metadata_field

!> Write a global text attribute to a file.
subroutine write_metadata_global(IO_handle, name, attribute)
  type(file_type),            intent(in)    :: IO_handle !< Handle for a file that is open for writing
  character(len=*),           intent(in)    :: name      !< The name in the file of this global attribute
  character(len=*),           intent(in)    :: attribute !< The value of this attribute

  call mpp_write_meta(IO_handle%unit, name, cval=attribute)
end subroutine write_metadata_global

end module MOM_io_infra
