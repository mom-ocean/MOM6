!> This module contains a thin inteface to mpp and fms I/O code
module MOM_io_infra

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domain_infra,     only : MOM_domain_type, rescale_comp_data, AGRID, BGRID_NE, CGRID_NE
use MOM_domain_infra,     only : domain2d, domain1d, CENTER, CORNER, NORTH_FACE, EAST_FACE
use MOM_error_infra,      only : MOM_error=>MOM_err, NOTE, FATAL, WARNING, is_root_PE
use MOM_string_functions, only : lowercase

use fms2_io_mod,          only : fms2_open_file => open_file, check_if_open, fms2_close_file => close_file
use fms2_io_mod,          only : FmsNetcdfDomainFile_t, FmsNetcdfFile_t, fms2_read_data => read_data
use fms2_io_mod,          only : get_unlimited_dimension_name, get_num_dimensions, get_num_variables
use fms2_io_mod,          only : get_variable_names, variable_exists, get_variable_size, get_variable_units
use fms2_io_mod,          only : register_field, write_data, register_variable_attribute, register_global_attribute
use fms2_io_mod,          only : variable_att_exists, get_variable_attribute, get_variable_num_dimensions
use fms2_io_mod,          only : get_variable_dimension_names, is_dimension_registered, get_dimension_size
use fms2_io_mod,          only : is_dimension_unlimited, register_axis, unlimited
use fms2_io_mod,          only : get_global_io_domain_indices

use fms_mod,              only : write_version_number, open_namelist_file, check_nml_error
use fms_io_mod,           only : file_exist, field_exist, field_size, read_data
use fms_io_mod,           only : fms_io_exit, get_filename_appendix
use mpp_domains_mod,      only : mpp_get_compute_domain, mpp_get_global_domain
use mpp_io_mod,           only : mpp_open, mpp_close, mpp_flush
use mpp_io_mod,           only : mpp_write_meta, mpp_write
use mpp_io_mod,           only : mpp_get_atts, mpp_attribute_exist
use mpp_io_mod,           only : mpp_get_axes, mpp_axistype=>axistype, mpp_get_axis_data
use mpp_io_mod,           only : mpp_get_fields, mpp_fieldtype=>fieldtype
use mpp_io_mod,           only : mpp_get_info, mpp_get_times
use mpp_io_mod,           only : mpp_io_init
use mpp_mod,              only : stdout_if_root=>stdout
! These are encoding constants.
use mpp_io_mod,           only : APPEND_FILE=>MPP_APPEND, WRITEONLY_FILE=>MPP_WRONLY
use mpp_io_mod,           only : OVERWRITE_FILE=>MPP_OVERWR, READONLY_FILE=>MPP_RDONLY
use mpp_io_mod,           only : NETCDF_FILE=>MPP_NETCDF, ASCII_FILE=>MPP_ASCII
use mpp_io_mod,           only : MULTIPLE=>MPP_MULTI, SINGLE_FILE=>MPP_SINGLE
use iso_fortran_env,      only : int64

implicit none ; private

! These interfaces are actually implemented or have explicit interfaces in this file.
public :: open_file, open_ASCII_file, file_is_open, close_file, flush_file, file_exists
public :: get_file_info, get_file_fields, get_file_times, get_filename_suffix
public :: MOM_read_data, MOM_read_vector, write_metadata, write_field
public :: field_exists, get_field_atts, get_field_size, get_axis_data, read_field_chksum
public :: io_infra_init, io_infra_end, MOM_namelist_file, check_namelist_error, write_version
public :: stdout_if_root
! These types act as containers for information about files, fields and axes, respectively,
! and may also wrap opaque types from the underlying infrastructure.
public :: file_type, fieldtype, axistype
! These are encoding constant parmeters.
public :: ASCII_FILE, NETCDF_FILE, SINGLE_FILE, MULTIPLE
public :: APPEND_FILE, READONLY_FILE, OVERWRITE_FILE, WRITEONLY_FILE
public :: CENTER, CORNER, NORTH_FACE, EAST_FACE

!> Indicate whether a file exists, perhaps with domain decomposition
interface file_exists
  module procedure FMS_file_exists
  module procedure MOM_file_exists
end interface

!> Open a file (or fileset) for parallel or single-file I/O.
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
type :: file_type ; private
  integer :: unit = -1 !< The framework identfier or netCDF unit number of an output file
  type(FmsNetcdfDomainFile_t), pointer :: fileobj => NULL() !< A domain-decomposed
                       !! file object that is open for writing
  character(len=:), allocatable :: filename !< The path to this file, if it is open
  logical :: open_to_read  = .false. !< If true, this file or fileset can be read
  logical :: open_to_write = .false. !< If true, this file or fileset can be written to
  integer :: num_times !< The number of time levels in this file
  real    :: file_time !< The time of the latest entry in the file.
  logical :: FMS2_file !< If true, this file-type is to be used with FMS2 interfaces.
end type file_type

!> This type is a container for information about a variable in a file.
type :: fieldtype ; private
  character(len=256)  :: name !< The name of this field in the files.
  type(mpp_fieldtype) :: FT !< The FMS1 field-type that this type wraps
  character(len=:), allocatable :: longname !< The long name for this field
  character(len=:), allocatable :: units    !< The units for this field
  integer(kind=int64) :: chksum_read !< A checksum that has been read from a file
  logical :: valid_chksum !< If true, this field has a valid checksum value.
  logical :: FMS2_field  !< If true, this field-type should be used with FMS2 interfaces.
end type fieldtype

!> This type is a container for information about an axis in a file.
type :: axistype ; private
  character(len=256) :: name !< The name of this axis in the files.
  type(mpp_axistype) :: AT   !< The FMS1 axis-type that this type wraps
  real, allocatable, dimension(:) :: ax_data !< The values of the data on the axis.
  logical :: domain_decomposed = .false.  !< True if axis is domain-decomposed
end type axistype

!> For now, these module-variables are hard-coded to exercise the new FMS2 interfaces.
logical :: FMS2_reads  = .true.
logical :: FMS2_writes = .true.

contains

!> Reads the checksum value for a field that was recorded in a file, along with a flag indicating
!! whether the file contained a valid checksum for this field.
subroutine read_field_chksum(field, chksum, valid_chksum)
  type(fieldtype),     intent(in)  :: field !< The field whose checksum attribute is to be read.
  integer(kind=int64), intent(out) :: chksum !< The checksum for the field.
  logical,             intent(out) :: valid_chksum  !< If true, chksum has been successfully read.

  chksum = -1
  valid_chksum = field%valid_chksum
  if (valid_chksum) chksum = field%chksum_read

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

  file_is_open = ((IO_handle%unit >= 0) .or. associated(IO_handle%fileobj))
end function file_is_open

!> closes a file (or fileset).  If the file handle does not point to an open file,
!! close_file_type simply returns without doing anything.
subroutine close_file_type(IO_handle)
  type(file_type), intent(inout) :: IO_handle   !< The I/O handle for the file to be closed

  if (associated(IO_handle%fileobj)) then
    call fms2_close_file(IO_handle%fileobj)
    deallocate(IO_handle%fileobj)
  else
    call mpp_close(IO_handle%unit)
  endif
  if (allocated(IO_handle%filename)) deallocate(IO_handle%filename)
  IO_handle%open_to_read = .false. ; IO_handle%open_to_write = .false.
  IO_handle%num_times = 0 ; IO_handle%file_time = 0.0
  IO_handle%FMS2_file = .false.
end subroutine close_file_type

!> closes a file.  If the unit does not point to an open file,
!! close_file_unit simply returns without doing anything.
subroutine close_file_unit(unit)
  integer, intent(inout) :: unit   !< The I/O unit for the file to be closed

  call mpp_close(unit)
end subroutine close_file_unit

!> Ensure that the output stream associated with a file handle is fully sent to disk.
subroutine flush_file_type(IO_handle)
  type(file_type), intent(in) :: IO_handle    !< The I/O handle for the file to flush

  if (associated(IO_handle%fileobj)) then
    ! There does not appear to be an fms2 flush call.
  else
    call mpp_flush(IO_handle%unit)
  endif
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

  ! Local variables
  type(FmsNetcdfDomainFile_t) :: fileobj_read ! A handle to a domain-decomposed file for obtaining information
                                              ! about the exiting time axis entries in append mode.
  logical :: success         ! If true, the file was opened successfully
  integer :: file_mode       ! An integer that encodes whether the file is to be opened for
                             ! reading, writing or appending
  character(len=40)  :: mode ! A character string that encodes whether the file is to be opened for
                             ! reading, writing or appending
  character(len=:), allocatable :: filename_tmp  ! A copy of filename with .nc appended if necessary.
  character(len=256) :: dim_unlim_name ! name of the unlimited dimension in the file
  integer :: index_nc

  if (IO_handle%open_to_write) then
    call MOM_error(WARNING, "open_file_type called for file "//trim(filename)//&
        " with an IO_handle that is already open to to write.")
    return
  endif
  if (IO_handle%open_to_read) then
    call MOM_error(FATAL, "open_file_type called for file "//trim(filename)//&
        " with an IO_handle that is already open to to read.")
  endif

  file_mode = WRITEONLY_FILE ; if (present(action)) file_mode = action

  if (FMS2_writes .and. present(MOM_Domain)) then
    if (.not.associated(IO_handle%fileobj)) allocate (IO_handle%fileobj)

    ! The FMS1 interface automatically appends .nc if necessary, but FMS2 interface does not.
    index_nc = index(trim(filename), ".nc")
    if (index_nc > 0) then
      filename_tmp = trim(filename)
    else
      filename_tmp = trim(filename)//".nc"
      if (is_root_PE()) call MOM_error(WARNING, "Open_file is appending .nc to the filename "//trim(filename))
    endif

    if (file_mode == WRITEONLY_FILE) then ; mode = "write"
    elseif (file_mode == APPEND_FILE) then ; mode = "append"
    elseif (file_mode == OVERWRITE_FILE) then ; mode = "overwrite"
    elseif (file_mode == READONLY_FILE) then ; mode = "read"
    else
      call MOM_error(FATAL, "open_file_type called with unrecognized action.")
    endif

    IO_handle%num_times = 0
    IO_handle%file_time = 0.0
    if ((file_mode == APPEND_FILE) .and. file_exists(filename_tmp, MOM_Domain)) then
      ! Determine the latest file time and number of records so far.
      success = fms2_open_file(fileObj_read, trim(filename_tmp), "read", MOM_domain%mpp_domain)
      call get_unlimited_dimension_name(fileObj_read, dim_unlim_name)
      if (len_trim(dim_unlim_name) > 0) &
        call get_dimension_size(fileObj_read, trim(dim_unlim_name), IO_handle%num_times)
      if (IO_handle%num_times > 0) &
        call fms2_read_data(fileObj_read, trim(dim_unlim_name), IO_handle%file_time, &
                            unlim_dim_level=IO_handle%num_times)
      call fms2_close_file(fileObj_read)
    endif

    success = fms2_open_file(IO_handle%fileobj, trim(filename_tmp), trim(mode), MOM_domain%mpp_domain)
    if (.not.success) call MOM_error(FATAL, "Unable to open file "//trim(filename_tmp))
    IO_handle%FMS2_file = .true.
  elseif (present(MOM_Domain)) then
    call mpp_open(IO_handle%unit, filename, action=file_mode, form=NETCDF_FILE, threading=threading, &
                  fileset=fileset, domain=MOM_Domain%mpp_domain)
    IO_handle%FMS2_file = .false.
  else
    call mpp_open(IO_handle%unit, filename, action=file_mode, form=NETCDF_FILE, threading=threading, &
                  fileset=fileset)
    IO_handle%FMS2_file = .false.
  endif
  IO_handle%filename = trim(filename)

  if (file_mode == READONLY_FILE) then
    IO_handle%open_to_read = .true. ; IO_handle%open_to_write = .false.
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
  character(len=256) :: dim_unlim_name ! name of the unlimited dimension in the file
  integer :: ndims, nvars, natts, ntimes

  if (IO_handle%FMS2_file) then
    if (present(ndim)) ndim = get_num_dimensions(IO_handle%fileobj)
    if (present(nvar)) nvar = get_num_variables(IO_handle%fileobj)
    if (present(ntime)) then
      ntime = 0
      call get_unlimited_dimension_name(IO_handle%fileobj, dim_unlim_name)
      if (len_trim(dim_unlim_name) > 0) &
        call get_dimension_size(IO_handle%fileobj, trim(dim_unlim_name), ntime)
    endif
  else
    call mpp_get_info(IO_handle%unit, ndims, nvars, natts, ntimes )

    if (present(ndim)) ndim = ndims
    if (present(nvar)) nvar = nvars
    if (present(ntime)) ntime = ntimes
  endif

end subroutine get_file_info


!> Get the times of records from a file
subroutine get_file_times(IO_handle, time_values, ntime)
  type(file_type),                 intent(in)    :: IO_handle !< Handle for a file that is open for I/O
  real, allocatable, dimension(:), intent(inout) :: time_values !< The real times for the records in file.
  integer,               optional, intent(out)   :: ntime !< The number of time levels in the file

  character(len=256) :: dim_unlim_name ! name of the unlimited dimension in the file
  integer :: ntimes  ! The number of time levels in the file

  !### Modify this routine to optionally convert to time_type, using information about the dimensions?

  if (allocated(time_values)) deallocate(time_values)
  call get_file_info(IO_handle, ntime=ntimes)
  if (present(ntime)) ntime = ntimes
  if (ntimes > 0) then
    allocate(time_values(ntimes))
    if (IO_handle%FMS2_file) then
      call get_unlimited_dimension_name(IO_handle%fileobj, dim_unlim_name)
      call fms2_read_data(IO_handle%fileobj, trim(dim_unlim_name), time_values)
    else
      call mpp_get_times(IO_handle%unit, time_values)
    endif
  endif
end subroutine get_file_times

!> Set up the field information (e.g., names and metadata) for all of the variables in a file.  The
!! argument fields must be allocated with a size that matches the number of variables in a file.
subroutine get_file_fields(IO_handle, fields)
  type(file_type),               intent(in)    :: IO_handle !< Handle for a file that is open for I/O
  type(fieldtype), dimension(:), intent(inout) :: fields !< Field-type descriptions of all of
                                                         !! the variables in a file.
  type(mpp_fieldtype), dimension(size(fields)) :: mpp_fields ! Fieldtype structures for the variables
  character(len=256),  dimension(size(fields)) :: var_names ! The names of all variables
  character(len=256)  :: units    ! The units of a variable as recorded in the file
  character(len=2048) :: longname ! The long-name of a variable as recorded in the file
  character(len=64)   :: checksum_char ! The hexadecimal checksum read from the file
  integer(kind=int64), dimension(3) :: checksum_file ! The checksums for a variable in the file
  integer :: nvar  ! The number of variables in the file
  integer :: i

  nvar = size(fields)
  ! Local variables
  if (IO_handle%FMS2_file) then
    call get_variable_names(IO_handle%fileobj, var_names)
    do i=1,nvar
      fields(i)%name = trim(var_names(i))
      longname = ""
      if (variable_att_exists(IO_handle%fileobj, var_names(i), "long_name")) &
        call get_variable_attribute(IO_handle%fileobj, var_names(i), "long_name", longname)
      fields(i)%longname = trim(longname)
      units = ""
      if (variable_att_exists(IO_handle%fileobj, var_names(i), "units")) &
        call get_variable_attribute(IO_handle%fileobj, var_names(i), "units", units)
      fields(i)%units = trim(units)

      fields(i)%valid_chksum = variable_att_exists(IO_handle%fileobj, var_names(i), "checksum")
      if (fields(i)%valid_chksum) then
        call get_variable_attribute(IO_handle%fileobj, var_names(i), 'checksum', checksum_char)
        ! If there are problems, there might need to be code added to handle commas.
        read (checksum_char(1:16), '(Z16)') fields(i)%chksum_read
      endif
    enddo
  else
    call mpp_get_fields(IO_handle%unit, mpp_fields)
    do i=1,nvar
      fields(i)%FT = mpp_fields(i)
      call mpp_get_atts(fields(i)%FT, name=fields(i)%name, units=units, longname=longname, &
                        checksum=checksum_file)
      fields(i)%longname = trim(longname)
      fields(i)%units = trim(units)
      fields(i)%valid_chksum = mpp_attribute_exist(fields(i)%FT, "checksum")
      if (fields(i)%valid_chksum) fields(i)%chksum_read = checksum_file(1)
    enddo
  endif

end subroutine get_file_fields

!> Extract information from a field type, as stored or as found in a file
subroutine get_field_atts(field, name, units, longname, checksum)
  type(fieldtype),            intent(in)  :: field !< The field to extract information from
  character(len=*), optional, intent(out) :: name  !< The variable name
  character(len=*), optional, intent(out) :: units !< The units of the variable
  character(len=*), optional, intent(out) :: longname  !< The long name of the variable
  integer(kind=int64),  dimension(:), &
                    optional, intent(out) :: checksum !< The checksums of the variable in a file

  if (present(name)) name = trim(field%name)
  if (present(units)) units = trim(field%units)
  if (present(longname)) longname = trim(field%longname)
  if (present(checksum)) checksum = field%chksum_read

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

  ! Local variables
  type(FmsNetcdfDomainFile_t) :: fileObj_dd ! A handle to a domain-decomposed file for obtaining information
                                            ! about the exiting time axis entries in append mode.
  type(FmsNetcdfFile_t) :: fileObj_simple   ! A handle to a non-domain-decomposed file for obtaining information
                                            ! about the exiting time axis entries in append mode.
  logical :: success         ! If true, the file was opened successfully
  logical :: domainless      ! If true, this file does not use a domain-decomposed file.

  domainless = .not.(present(MOM_domain) .or. present(domain))
  if (present(no_domain)) then
    if (domainless .and. .not.no_domain) call MOM_error(FATAL, &
        "field_exists: When no_domain is present and false, a domain must be supplied in query about "//&
        trim(field_name)//" in file "//trim(filename))
    domainless = no_domain
  endif

  if (FMS2_reads) then
    field_exists = .false.
    if (file_exists(filename)) then
      if (domainless) then
        success = fms2_open_file(fileObj_simple, trim(filename), "read")
        if (success) then
          field_exists = variable_exists(fileObj_simple, field_name)
          call fms2_close_file(fileObj_simple)
        endif
      else
        if (present(MOM_domain)) then
          success = fms2_open_file(fileObj_dd, trim(filename), "read", MOM_domain%mpp_domain)
        else
          success = fms2_open_file(fileObj_dd, trim(filename), "read", domain)
        endif
        if (success) then
          field_exists = variable_exists(fileobj_dd, field_name)
          call fms2_close_file(fileObj_dd)
        endif
      endif
    endif
  elseif (present(MOM_domain)) then
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

  ! Local variables
  type(FmsNetcdfFile_t) :: fileobj_read ! A handle to a non-domain-decomposed file for obtaining information
                                              ! about the exiting time axis entries in append mode.
  logical :: success         ! If true, the file was opened successfully
  logical :: field_exists    ! True if filename exists and field_name is in filename
  integer :: i, ndims

  if (FMS2_reads) then
    field_exists = .false.
    if (file_exists(filename)) then
      success = fms2_open_file(fileObj_read, trim(filename), "read")
      if (success) then
        field_exists = variable_exists(fileobj_read, fieldname)
        if (field_exists) then
          ndims = get_variable_num_dimensions(fileobj_read, fieldname)
          if (ndims > size(sizes)) call MOM_error(FATAL, &
            "get_field_size called with too few sizes for "//trim(fieldname)//" in "//trim(filename))
          call get_variable_size(fileobj_read, fieldname, sizes(1:ndims))
          do i=ndims+1,size(sizes) ; sizes(i) = 0 ; enddo
        endif
      endif
    endif
    if (present(field_found)) field_found = field_exists
  else
    call field_size(filename, fieldname, sizes, field_found=field_found, no_domain=no_domain)
  endif

end subroutine get_field_size

!> Extracts and returns the axis data stored in an axistype.
subroutine get_axis_data( axis, dat )
  type(axistype),     intent(in)  :: axis !< An axis type
  real, dimension(:), intent(out) :: dat  !< The data in the axis variable

  integer :: i

  ! This routine might not be needed for MOM6.
  if (allocated(axis%ax_data)) then
    if (size(axis%ax_data) > size(dat)) call MOM_error(FATAL, &
      "get_axis_data called with too small of an output data array for "//trim(axis%name))
    do i=1,size(axis%ax_data) ; dat(i) = axis%ax_data(i) ; enddo
  elseif (.not.FMS2_writes) then
    call mpp_get_axis_data( axis%AT, dat )
  endif

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
  logical,      optional, intent(in)    :: file_may_be_4d !< If true, this file may have 4-d arrays, but
                                                     !! with the FMS2 I/O interfaces this does not matter.

  ! Local variables
  type(FmsNetcdfFile_t)       :: fileObj ! A handle to a non-domain-decomposed file
  type(FmsNetcdfDomainFile_t) :: fileobj_DD ! A handle to a domain-decomposed file object
  character(len=96) :: var_to_read ! Name of variable to read from the netcdf file
  logical :: has_time_dim          ! True if the variable has an unlimited time axis.
  logical :: success               ! True if the file was successfully opened

  if (present(MOM_Domain) .and. FMS2_reads) then
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileobj_DD, filename, "read", MOM_domain%mpp_domain)
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file and prepare to read it.
    call prepare_to_read_var(fileobj_DD, fieldname, "MOM_read_data_0d: ", filename, &
                             var_to_read, has_time_dim, timelevel)

    ! Read the data.
    if (present(timelevel) .and. has_time_dim) then
      call fms2_read_data(fileobj_DD, var_to_read, data, unlim_dim_level=timelevel)
    else
      call fms2_read_data(fileobj_DD, var_to_read, data)
    endif

    ! Close the file-set.
    if (check_if_open(fileobj_DD)) call fms2_close_file(fileobj_DD)
  elseif (FMS2_reads) then
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileObj, trim(filename), "read")
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file, and determine whether it
    ! has a time dimension.
    call find_varname_in_file(fileObj, fieldname, "MOM_read_data_0d: ", filename, &
                              var_to_read, has_time_dim, timelevel)

    ! Read the data.
    if (present(timelevel) .and. has_time_dim) then
      call fms2_read_data(fileobj, var_to_read, data, unlim_dim_level=timelevel)
    else
      call fms2_read_data(fileobj, var_to_read, data)
    endif

    ! Close the file-set.
    if (check_if_open(fileobj)) call fms2_close_file(fileobj)
  elseif (present(MOM_Domain)) then ! Read the variable using the FMS-1 interface.
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
  logical,      optional, intent(in)    :: file_may_be_4d !< If true, this file may have 4-d arrays, but
                                                     !! with the FMS2 I/O interfaces this does not matter.

  ! Local variables
  type(FmsNetcdfFile_t)       :: fileObj ! A handle to a non-domain-decomposed file
  type(FmsNetcdfDomainFile_t) :: fileobj_DD ! A handle to a domain-decomposed file object
  character(len=96) :: var_to_read ! Name of variable to read from the netcdf file
  logical :: has_time_dim          ! True if the variable has an unlimited time axis.
  logical :: success               ! True if the file was successfully opened

  if (present(MOM_Domain) .and. FMS2_reads) then
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileobj_DD, filename, "read", MOM_domain%mpp_domain)
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file and prepare to read it.
    call prepare_to_read_var(fileobj_DD, fieldname, "MOM_read_data_1d: ", filename, &
                             var_to_read, has_time_dim, timelevel)

    ! Read the data.
    if (present(timelevel) .and. has_time_dim) then
      call fms2_read_data(fileobj_DD, var_to_read, data, unlim_dim_level=timelevel)
    else
      call fms2_read_data(fileobj_DD, var_to_read, data)
    endif

    ! Close the file-set.
    if (check_if_open(fileobj_DD)) call fms2_close_file(fileobj_DD)
  elseif (FMS2_reads) then
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileObj, trim(filename), "read")
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file, and determine whether it
    ! has a time dimension.
    call find_varname_in_file(fileObj, fieldname, "MOM_read_data_1d: ", filename, &
                              var_to_read, has_time_dim, timelevel)

    ! Read the data.
    if (present(timelevel) .and. has_time_dim) then
      call fms2_read_data(fileobj, var_to_read, data, unlim_dim_level=timelevel)
    else
      call fms2_read_data(fileobj, var_to_read, data)
    endif

    ! Close the file-set.
    if (check_if_open(fileobj)) call fms2_close_file(fileobj)
  elseif (present(MOM_Domain)) then ! Read the variable using the FMS-1 interface.
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
  logical,      optional, intent(in)    :: file_may_be_4d !< If true, this file may have 4-d arrays, but
                                                     !! with the FMS2 I/O interfaces this does not matter.

  ! Local variables
  type(FmsNetcdfDomainFile_t) :: fileobj ! A handle to a domain-decomposed file object
  character(len=96) :: var_to_read ! Name of variable to read from the netcdf file
  logical :: has_time_dim          ! True if the variable has an unlimited time axis.
  logical :: success               ! True if the file was successfully opened

  if (FMS2_reads) then
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileobj, filename, "read", MOM_domain%mpp_domain)
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file and prepare to read it.
    call prepare_to_read_var(fileobj, fieldname, "MOM_read_data_2d: ", filename, &
                             var_to_read, has_time_dim, timelevel, position)

    ! Read the data.
    if (present(timelevel) .and. has_time_dim) then
      call fms2_read_data(fileobj, var_to_read, data, unlim_dim_level=timelevel)
    else
      call fms2_read_data(fileobj, var_to_read, data)
    endif

    ! Close the file-set.
    if (check_if_open(fileobj)) call fms2_close_file(fileobj)
  else ! Read the variable using the FMS-1 interface.
    call read_data(filename, fieldname, data, MOM_Domain%mpp_domain, &
                   timelevel=timelevel, position=position)
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

  ! Local variables
  type(FmsNetcdfFile_t)       :: fileObj ! A handle to a non-domain-decomposed file
  type(FmsNetcdfDomainFile_t) :: fileobj_DD ! A handle to a domain-decomposed file object
  character(len=96) :: var_to_read ! Name of variable to read from the netcdf file
  logical :: success               ! True if the file was successfully opened

  if (present(MOM_Domain) .and. FMS2_reads) then
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileobj_DD, filename, "read", MOM_domain%mpp_domain)
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file and prepare to read it.
    call prepare_to_read_var(fileobj_DD, fieldname, "MOM_read_data_2d_region: ", &
                             filename, var_to_read)

    ! Read the data.
    call fms2_read_data(fileobj_DD, var_to_read, data, corner=start(1:2), edge_lengths=nread(1:2))

    ! Close the file-set.
    if (check_if_open(fileobj_DD)) call fms2_close_file(fileobj_DD)
  elseif (FMS2_reads) then
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileObj, trim(filename), "read")
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file, and determine whether it
    ! has a time dimension.
    call find_varname_in_file(fileObj, fieldname, "MOM_read_data_2d_region: ", filename, var_to_read)

    ! Read the data.
    call fms2_read_data(fileobj, var_to_read, data, corner=start(1:2), edge_lengths=nread(1:2))

    ! Close the file-set.
    if (check_if_open(fileobj)) call fms2_close_file(fileobj)
  elseif (present(MOM_Domain)) then ! Read the variable using the FMS-1 interface.
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
  logical,      optional, intent(in)    :: file_may_be_4d !< If true, this file may have 4-d arrays, but
                                                     !! with the FMS2 I/O interfaces this does not matter.

  ! Local variables
  type(FmsNetcdfDomainFile_t) :: fileobj ! A handle to a domain-decomposed file object
  character(len=96) :: var_to_read ! Name of variable to read from the netcdf file
  logical :: has_time_dim          ! True if the variable has an unlimited time axis.
  logical :: success               ! True if the file was successfully opened

  if (FMS2_reads) then
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileobj, filename, "read", MOM_domain%mpp_domain)
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file and prepare to read it.
    call prepare_to_read_var(fileobj, fieldname, "MOM_read_data_3d: ", filename, &
                             var_to_read, has_time_dim, timelevel, position)

    ! Read the data.
    if (present(timelevel) .and. has_time_dim) then
      call fms2_read_data(fileobj, var_to_read, data, unlim_dim_level=timelevel)
    else
      call fms2_read_data(fileobj, var_to_read, data)
    endif

    ! Close the file-set.
    if (check_if_open(fileobj)) call fms2_close_file(fileobj)
  else ! Read the variable using the FMS-1 interface.
    call read_data(filename, fieldname, data, MOM_Domain%mpp_domain, &
                   timelevel=timelevel, position=position)
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
  type(FmsNetcdfDomainFile_t) :: fileobj ! A handle to a domain-decomposed file object
  logical :: has_time_dim          ! True if the variable has an unlimited time axis.
  character(len=96) :: var_to_read ! Name of variable to read from the netcdf file
  logical :: success  ! True if the file was successfully opened

  if (FMS2_reads) then
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileobj, filename, "read", MOM_domain%mpp_domain)
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file and prepare to read it.
    call prepare_to_read_var(fileobj, fieldname, "MOM_read_data_4d: ", filename, &
                             var_to_read, has_time_dim, timelevel, position)

    ! Read the data.
    if (present(timelevel) .and. has_time_dim) then
      call fms2_read_data(fileobj, var_to_read, data, unlim_dim_level=timelevel)
    else
      call fms2_read_data(fileobj, var_to_read, data)
    endif

    ! Close the file-set.
    if (check_if_open(fileobj)) call fms2_close_file(fileobj)
  else ! Read the variable using the FMS-1 interface.
    call read_data(filename, fieldname, data, MOM_Domain%mpp_domain, &
                   timelevel=timelevel, position=position)
  endif

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

  ! Local variables
  type(FmsNetcdfFile_t) :: fileObj ! A handle to a non-domain-decomposed file
  logical :: has_time_dim          ! True if the variable has an unlimited time axis.
  character(len=96) :: var_to_read ! Name of variable to read from the netcdf file
  logical :: success               ! If true, the file was opened successfully

  ! This routine might not be needed for MOM6.
  if (FMS2_reads) then
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileObj, trim(filename), "read")
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file, and determine whether it
    ! has a time dimension.
    call find_varname_in_file(fileObj, fieldname, "MOM_read_data_0d_int: ", filename, &
                              var_to_read, has_time_dim, timelevel)

    ! Read the data.
    if (present(timelevel) .and. has_time_dim) then
      call fms2_read_data(fileobj, var_to_read, data, unlim_dim_level=timelevel)
    else
      call fms2_read_data(fileobj, var_to_read, data)
    endif

    ! Close the file-set.
    if (check_if_open(fileobj)) call fms2_close_file(fileobj)
  else
    call read_data(filename, fieldname, data, timelevel=timelevel, no_domain=.true.)
  endif

end subroutine MOM_read_data_0d_int

!> This routine uses the fms_io subroutine read_data to read a 1-D integer
!! data field named "fieldname" from file "filename".
subroutine MOM_read_data_1d_int(filename, fieldname, data, timelevel)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  integer, dimension(:),  intent(inout) :: data      !< The 1-dimensional array into which the data
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read

  ! Local variables
  type(FmsNetcdfFile_t) :: fileObj ! A handle to a non-domain-decomposed file for obtaining information
                                   ! about the exiting time axis entries in append mode.
  logical :: has_time_dim          ! True if the variable has an unlimited time axis.
  character(len=96) :: var_to_read ! Name of variable to read from the netcdf file
  logical :: success               ! If true, the file was opened successfully

  ! This routine might not be needed for MOM6.
  if (FMS2_reads) then
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileObj, trim(filename), "read")
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive variable name in the file, and determine whether it
    ! has a time dimension.
    call find_varname_in_file(fileObj, fieldname, "MOM_read_data_1d_int: ", filename, &
                              var_to_read, has_time_dim, timelevel)

    ! Read the data.
    if (present(timelevel) .and. has_time_dim) then
      call fms2_read_data(fileobj, var_to_read, data, unlim_dim_level=timelevel)
    else
      call fms2_read_data(fileobj, var_to_read, data)
    endif

    ! Close the file-set.
    if (check_if_open(fileobj)) call fms2_close_file(fileobj)
  else
    call read_data(filename, fieldname, data, timelevel=timelevel, no_domain=.true.)
  endif

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
  ! Local variables
  type(FmsNetcdfDomainFile_t) :: fileobj ! A handle to a domain-decomposed file object
  logical :: has_time_dim           ! True if the variables have an unlimited time axis.
  character(len=96) :: u_var, v_var ! Name of u and v variables to read from the netcdf file
  logical :: success                ! True if the file was successfully opened
  integer :: u_pos, v_pos           ! Flags indicating the positions of the u- and v- components.

  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == BGRID_NE) then ; u_pos = CORNER ; v_pos = CORNER
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif

  if (FMS2_reads) then
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileobj, filename, "read", MOM_domain%mpp_domain)
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive u- and v-variable names in the file and prepare to read them.
    call prepare_to_read_var(fileobj, u_fieldname, "MOM_read_vector_2d: ", filename, &
                             u_var, has_time_dim, timelevel, position=u_pos)
    call prepare_to_read_var(fileobj, v_fieldname, "MOM_read_vector_2d: ", filename, &
                             v_var, has_time_dim, timelevel, position=v_pos)

    ! Read the u-data and v-data. There would already been an error message for one
    ! of the variables if they are inconsistent in having an unlimited dimension.
    if (present(timelevel) .and. has_time_dim) then
      call fms2_read_data(fileobj, u_var, u_data, unlim_dim_level=timelevel)
      call fms2_read_data(fileobj, v_var, v_data, unlim_dim_level=timelevel)
    else
      call fms2_read_data(fileobj, u_var, u_data)
      call fms2_read_data(fileobj, v_var, v_data)
    endif

    ! Close the file-set.
    if (check_if_open(fileobj)) call fms2_close_file(fileobj)
  else ! Read the variable using the FMS-1 interface.
    call read_data(filename, u_fieldname, u_data, MOM_Domain%mpp_domain, &
                   timelevel=timelevel, position=u_pos)
    call read_data(filename, v_fieldname, v_data, MOM_Domain%mpp_domain, &
                   timelevel=timelevel, position=v_pos)
  endif

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
  logical,      optional, intent(in)    :: scalar_pair !< If true, a pair of scalars are to be read.
  real,         optional, intent(in)    :: scale     !< A scaling factor that the fields are multiplied
                                                     !! by before they are returned.

  ! Local variables
  type(FmsNetcdfDomainFile_t) :: fileobj ! A handle to a domain-decomposed file object
  logical :: has_time_dim           ! True if the variables have an unlimited time axis.
  character(len=96) :: u_var, v_var ! Name of u and v variables to read from the netcdf file
  logical :: success                ! True if the file was successfully opened
  integer :: u_pos, v_pos           ! Flags indicating the positions of the u- and v- components.

  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == BGRID_NE) then ; u_pos = CORNER ; v_pos = CORNER
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif

  if (FMS2_reads) then
    ! Open the FMS2 file-set.
    success = fms2_open_file(fileobj, filename, "read", MOM_domain%mpp_domain)
    if (.not.success) call MOM_error(FATAL, "Failed to open "//trim(filename))

    ! Find the matching case-insensitive u- and v-variable names in the file and prepare to read them.
    call prepare_to_read_var(fileobj, u_fieldname, "MOM_read_vector_3d: ", filename, &
                             u_var, has_time_dim, timelevel, position=u_pos)
    call prepare_to_read_var(fileobj, v_fieldname, "MOM_read_vector_3d: ", filename, &
                             v_var, has_time_dim, timelevel, position=v_pos)

    ! Read the u-data and v-data, dangerously assuming either both or neither have time dimensions.
    ! There would already been an error message for one of the variables if they are inconsistent.
    if (present(timelevel) .and. has_time_dim) then
      call fms2_read_data(fileobj, u_var, u_data, unlim_dim_level=timelevel)
      call fms2_read_data(fileobj, v_var, v_data, unlim_dim_level=timelevel)
    else
      call fms2_read_data(fileobj, u_var, u_data)
      call fms2_read_data(fileobj, v_var, v_data)
    endif

    ! Close the file-set.
    if (check_if_open(fileobj)) call fms2_close_file(fileobj)
  else ! Read the variable using the FMS-1 interface.
    call read_data(filename, u_fieldname, u_data, MOM_Domain%mpp_domain, &
                   timelevel=timelevel, position=u_pos)
    call read_data(filename, v_fieldname, v_data, MOM_Domain%mpp_domain, &
                   timelevel=timelevel, position=v_pos)
  endif

  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(MOM_Domain, u_data, scale)
    call rescale_comp_data(MOM_Domain, v_data, scale)
  endif ; endif

end subroutine MOM_read_vector_3d


!> Find the case-sensitive name of the variable in a netCDF file with a case-insensitive name match.
!! Optionally also determine whether this variable has an unlimited time dimension.
subroutine find_varname_in_file(fileobj, fieldname, err_header, filename, var_to_read, has_time_dim, timelevel)
  type(FmsNetcdfFile_t),       intent(inout) :: fileobj     !< An FMS2 handle to an open NetCDF file
  character(len=*),            intent(in)    :: fieldname   !< The variable name to seek in the file
  character(len=*),            intent(in)    :: err_header  !< A descriptive prefix for error messages
  character(len=*),            intent(in)    :: filename    !< The name of the file to read
  character(len=*),            intent(out)   :: var_to_read !< The variable name to read from the file
  logical,           optional, intent(out)   :: has_time_dim !< Indicates whether fieldname has a time dimension
  integer,           optional, intent(in)    :: timelevel   !< A time level to read

  ! Local variables
  logical :: variable_found ! Is a case-insensitive version of the variable found in the netCDF file?
  character(len=256), allocatable, dimension(:) :: var_names ! The names of all the variables in the netCDF file
  character(len=256), allocatable :: dim_names(:) ! The names of a variable's dimensions
  integer :: nvars          ! The number of variables in the file
  integer :: dim_unlim_size ! The current size of the unlimited (time) dimension in the file.
  integer :: num_var_dims   ! The number of dimensions a variable has in the file.
  integer :: time_dim       ! The position of the unlimited (time) dimension for a variable, or -1
                            ! if it has no unlimited dimension.
  integer :: i

  ! Open the file if necessary
  if (.not.check_if_open(fileobj))  &
    call MOM_error(FATAL, trim(err_header)//trim(filename)//" was not open in call to find_varname_in_file.")

  ! Search for the variable in the file, looking for the case-sensitive name first.
  if (variable_exists(fileobj, trim(fieldname))) then
    var_to_read = trim(fieldname)
  else ! Look for case-insensitive variable name matches.
    nvars = get_num_variables(fileobj)
    if (nvars < 1) call MOM_error(FATAL, "nvars is less than 1 for file "//trim(filename))
    allocate(var_names(nvars))
    call get_variable_names(fileobj, var_names)

    ! search for the variable in the file
    variable_found = .false.
    do i=1,nvars
      if (lowercase(trim(var_names(i))) == lowercase(trim(fieldname))) then
        variable_found = .true.
        var_to_read = trim(var_names(i))
        exit
      endif
    enddo
    if (.not.(variable_found)) &
      call MOM_error(FATAL, trim(err_header)//trim(fieldname)//" not found in "//trim(filename))
    deallocate(var_names)
  endif

  ! FMS2 can not handle a timelevel argument if the variable does not have one in the file,
  ! so some error checking and logic are required.
  if (present(has_time_dim) .or. present(timelevel)) then
    time_dim = -1

    num_var_dims = get_variable_num_dimensions(fileobj, trim(var_to_read))
    allocate(dim_names(num_var_dims)) ; dim_names(:) = ""
    call get_variable_dimension_names(fileobj, trim(var_to_read), dim_names)

    do i=1,num_var_dims
      if (is_dimension_unlimited(fileobj, dim_names(i))) then
        time_dim = i
        if (present(timelevel)) then
          call get_dimension_size(fileobj, dim_names(i), dim_unlim_size)
          if ((timelevel > dim_unlim_size) .and. is_root_PE()) call MOM_error(FATAL, &
                trim(err_header)//"Attempting to read a time level of "//trim(var_to_read)//&
                " that exceeds the size of the time dimension in "//trim(filename))
        endif
        exit
      endif
    enddo
    deallocate(dim_names)

    if (present(timelevel) .and. (time_dim < 0) .and. is_root_PE()) &
      call MOM_error(WARNING, trim(err_header)//"time level specified, but the variable "//&
                   trim(var_to_read)//" does not have an unlimited dimension in "//trim(filename))
    if ((.not.present(timelevel)) .and. (time_dim > 0) .and. is_root_PE()) &
      call MOM_error(WARNING, trim(err_header)//"The variable "//trim(var_to_read)//&
                    " has an unlimited dimension in "//trim(filename)//" but no time level is specified.")
    if (present(has_time_dim)) has_time_dim = (time_dim > 0)
  endif

end subroutine find_varname_in_file


!> Find the case-insensitive name match with a variable in an open domain-decomposed file-set,
!! prepare FMS2 to read this variable, and return some information needed to call fms2_read_data
!! correctly for this variable and file.
subroutine prepare_to_read_var(fileobj, fieldname, err_header, filename, var_to_read, &
                               has_time_dim, timelevel, position)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj     !< An FMS2 handle to an open domain-decomposed file
  character(len=*),            intent(in)    :: fieldname   !< The variable name to seek in the file
  character(len=*),            intent(in)    :: err_header  !< A descriptive prefix for error messages
  character(len=*),            intent(in)    :: filename    !< The name of the file to read
  character(len=*),            intent(out)   :: var_to_read !< The variable name to read from the file
  logical,           optional, intent(out)   :: has_time_dim !< Indicates whether fieldname has a time dimension
  integer,           optional, intent(in)    :: timelevel   !< A time level to read
  integer,           optional, intent(in)    :: position    !< A flag indicating where this variable is discretized

  ! Local variables
  logical :: variable_found ! Is a case-insensitive version of the variable found in the netCDF file?
  character(len=256), allocatable, dimension(:) :: var_names ! The names of all the variables in the netCDF file
  character(len=256), allocatable :: dim_names(:) ! The names of a variable's dimensions
  integer :: nvars          ! The number of variables in the file.
  integer :: dim_unlim_size ! The current size of the unlimited (time) dimension in the file.
  integer :: num_var_dims   ! The number of dimensions a variable has in the file.
  integer :: time_dim       ! The position of the unlimited (time) dimension for a variable, or -1
                            ! if it has no unlimited dimension.
  integer :: i

  ! Open the file if necessary
  if (.not.check_if_open(fileobj))  &
    call MOM_error(FATAL, trim(err_header)//trim(filename)//" was not open in call to prepare_to_read_var.")

  ! Search for the variable in the file, looking for the case-sensitive name first.
  if (variable_exists(fileobj, trim(fieldname))) then
    var_to_read = trim(fieldname)
  else  ! Look for case-insensitive variable name matches.
    nvars = get_num_variables(fileobj)
    if (nvars < 1) call MOM_error(FATAL, "nvars is less than 1 for file "//trim(filename))
    allocate(var_names(nvars))
    call get_variable_names(fileobj, var_names)

    variable_found = .false.
    do i=1,nvars
      if (lowercase(trim(var_names(i))) == lowercase(trim(fieldname))) then
        variable_found = .true.
        var_to_read = trim(var_names(i))
        exit
      endif
    enddo
    if (.not.(variable_found)) &
      call MOM_error(FATAL, trim(err_header)//trim(fieldname)//" not found in "//trim(filename))
    deallocate(var_names)
  endif

  ! FMS2 can not handle a timelevel argument if the variable does not have one in the file,
  ! so some error checking and logic are required.
  if (present(has_time_dim) .or. present(timelevel)) then
    time_dim = -1

    num_var_dims = get_variable_num_dimensions(fileobj, trim(var_to_read))
    allocate(dim_names(num_var_dims)) ; dim_names(:) = ""
    call get_variable_dimension_names(fileobj, trim(var_to_read), dim_names)

    do i=1,num_var_dims
      if (is_dimension_unlimited(fileobj, dim_names(i))) then
        time_dim = i
        if (present(timelevel)) then
          call get_dimension_size(fileobj, dim_names(i), dim_unlim_size)
          if ((timelevel > dim_unlim_size) .and. is_root_PE()) call MOM_error(FATAL, &
                trim(err_header)//"Attempting to read a time level of "//trim(var_to_read)//&
                " that exceeds the size of the time dimension in "//trim(filename))
        endif
        exit
      endif
    enddo
    deallocate(dim_names)

    if (present(timelevel) .and. (time_dim < 0) .and. is_root_PE()) &
      call MOM_error(WARNING, trim(err_header)//"time level specified, but the variable "//&
                   trim(var_to_read)//" does not have an unlimited dimension in "//trim(filename))
    if ((.not.present(timelevel)) .and. (time_dim > 0) .and. is_root_PE()) &
      call MOM_error(WARNING, trim(err_header)//"The variable "//trim(var_to_read)//&
                    " has an unlimited dimension in "//trim(filename)//" but no time level is specified.")
    if (present(has_time_dim)) has_time_dim = (time_dim > 0)
  endif

  ! Registering the variable axes essentially just specifies the discrete position of this variable.
  call MOM_register_variable_axes(fileobj, var_to_read, filename, position)

end subroutine prepare_to_read_var

!> register axes associated with a variable from a domain-decomposed netCDF file
subroutine MOM_register_variable_axes(fileObj, variableName, filename, position)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObj !< Handle to an open FMS2 netCDF file object
  character(len=*),  intent(in) :: variableName !< name of the variable
  character(len=*),  intent(in) :: filename     !< The name of the file to read
  integer, optional, intent(in) :: position     !< A flag indicating where this data is discretized

  ! Local variables
  character(len=256), allocatable, dimension(:) :: dim_names ! variable dimension names
  integer, allocatable, dimension(:) :: dimSizes ! variable dimension sizes
  logical, allocatable, dimension(:) :: is_x ! Is this a (likely domain-decomposed) x-axis
  logical, allocatable, dimension(:) :: is_y ! Is this a (likely domain-decomposed) y-axis
  logical, allocatable, dimension(:) :: is_t ! Is this a time axis or another unlimited axis
  integer :: ndims ! number of dimensions
  integer :: xPos, yPos ! Discrete positions for x and y axes. Default is CENTER
  integer :: i

  xPos = CENTER ; yPos = CENTER
  if (present(position)) then
    if ((position == CORNER) .or. (position == EAST_FACE)) xPos = EAST_FACE
    if ((position == CORNER) .or. (position == NORTH_FACE)) yPos = NORTH_FACE
  endif

  ! get variable dimension names and lengths
  ndims = get_variable_num_dimensions(fileObj, trim(variableName))
  allocate(dimSizes(ndims))
  allocate(dim_names(ndims))
  allocate(is_x(ndims)) ; is_x(:) = .false.
  allocate(is_y(ndims)) ; is_y(:) = .false.
  allocate(is_t(ndims)) ; is_t(:) = .false.
  call get_variable_size(fileObj, trim(variableName), dimSizes)
  call get_variable_dimension_names(fileObj, trim(variableName), dim_names)
  call categorize_axes(fileObj, filename, ndims, dim_names, is_x, is_y, is_t)

  ! register the axes
  do i=1,ndims
    if ( .not.is_dimension_registered(fileobj, trim(dim_names(i))) ) then
      if (is_x(i)) then
        call register_axis(fileObj, trim(dim_names(i)), "x", domain_position=xPos)
      elseif (is_y(i)) then
        call register_axis(fileObj, trim(dim_names(i)), "y", domain_position=yPos)
      else
        call register_axis(fileObj, trim(dim_names(i)), dimSizes(i))
      endif
    endif
  enddo

  deallocate(dimSizes, dim_names, is_x, is_y, is_t)
end subroutine MOM_register_variable_axes

!> Determine whether a variable's axes are associated with x-, y- or time-dimensions.  Other
!! unlimited dimensions are also labeled as time axes for these purposes.
subroutine categorize_axes(fileObj, filename, ndims, dim_names, is_x, is_y, is_t)
  type(FmsNetcdfDomainFile_t), intent(in)  :: fileObj  !< Handle to an open FMS2 netCDF file object
  character(len=*),            intent(in)  :: filename !< The name of the file to read
  integer,                     intent(in)  :: ndims    !< The number of dimensions associated with a variable
  character(len=*), dimension(ndims), intent(in) :: dim_names !< Names of the dimensions associated with a variable
  logical, dimension(ndims),   intent(out) :: is_x !< Indicates if each dimension a (likely decomposed) x-axis
  logical, dimension(ndims),   intent(out) :: is_y !< Indicates if each dimension a (likely decomposed) y-axis
  logical, dimension(ndims),   intent(out) :: is_t !< Indicates if each dimension unlimited (usually time) axis

  ! Local variables
  character(len=128) :: cartesian ! A flag indicating a Cartesian direction - usually a single character.
  character(len=512) :: dim_list  ! A concatenated list of dimension names.
  character(len=128) :: units ! units corresponding to a specific variable dimension
  logical :: x_found, y_found ! Indicate whether an x- or y- dimension have been found.
  integer :: i

  x_found = .false. ; y_found = .false.
  is_x(:) = .false. ; is_y(:) = .false.
  do i=1,ndims
    is_t(i) = is_dimension_unlimited(fileObj, trim(dim_names(i)))
    ! First look for indicative variable attributes
    if (.not.is_t(i)) then
      if (variable_exists(fileobj, trim(dim_names(i)))) then
        if (variable_att_exists(fileobj, trim(dim_names(i)), "cartesian_axis")) then
          call get_variable_attribute(fileobj, trim(dim_names(i)), "cartesian_axis", cartesian)
          cartesian = adjustl(cartesian)
          if ((index(cartesian, "X") == 1) .or. (index(cartesian, "x") == 1)) is_x(i) = .true.
          if ((index(cartesian, "Y") == 1) .or. (index(cartesian, "y") == 1)) is_y(i) = .true.
          if ((index(cartesian, "T") == 1) .or. (index(cartesian, "t") == 1)) is_t(i) = .true.
        endif
      endif
    endif
    if (is_x(i)) x_found = .true.
    if (is_y(i)) y_found = .true.
  enddo

  if (.not.(x_found .and. y_found)) then
    ! Next look for hints from axis names for uncharacterized axes
    do i=1,ndims ; if (.not.(is_x(i) .or. is_y(i) .or. is_t(i))) then
      call categorize_axis_from_name(dim_names(i), is_x(i), is_y(i))
      if (is_x(i)) x_found = .true.
      if (is_y(i)) y_found = .true.
    endif ; enddo
  endif

  if (.not.(x_found .and. y_found)) then
    ! Look for hints from CF-compliant axis units for uncharacterized axes
    do i=1,ndims ; if (.not.(is_x(i) .or. is_y(i) .or. is_t(i))) then
      call get_variable_units(fileobj, trim(dim_names(i)), units)
      call categorize_axis_from_units(units, is_x(i), is_y(i))
      if (is_x(i)) x_found = .true.
      if (is_y(i)) y_found = .true.
    endif ; enddo
  endif

  if (.not.(x_found .and. y_found) .and. ((ndims>2) .or. ((ndims==2) .and. .not.is_t(ndims)))) then
    ! This is a case where one would expect to find x-and y-dimensions, but none have been found.
    if (is_root_pe()) then
      dim_list = trim(dim_names(1))//", "//trim(dim_names(2))
      do i=3,ndims ; dim_list = trim(dim_list)//", "//trim(dim_names(i)) ; enddo
      call MOM_error(WARNING, "categorize_axes: Failed to identify x- and y- axes in the axis list ("//&
                     trim(dim_list)//") of a variable being read from "//trim(filename))
    endif
  endif

end subroutine categorize_axes

!> Determine whether an axis is associated with the x- or y-directions based on a comparison of
!! its units with CF-compliant variants of latitude or longitude units.
subroutine categorize_axis_from_units(unit_string, is_x, is_y)
  character(len=*), intent(in) :: unit_string !< string of units
  logical, intent(out) :: is_x !< Indicates if the axis units are associated with an x-direction axis
  logical, intent(out) :: is_y !< Indicates if the axis units are associated with an y-direction axis

  is_x = .false. ; is_y = .false.
  select case (lowercase(trim(unit_string)))
    case ("degrees_north"); is_y = .true.
    case ("degree_north") ; is_y = .true.
    case ("degrees_n")    ; is_y = .true.
    case ("degree_n")     ; is_y = .true.
    case ("degreen")      ; is_y = .true.
    case ("degreesn")     ; is_y = .true.
    case ("degrees_east") ; is_x = .true.
    case ("degree_east")  ; is_x = .true.
    case ("degreese")     ; is_x = .true.
    case ("degreee")      ; is_x = .true.
    case ("degree_e")     ; is_x = .true.
    case ("degrees_e")    ; is_x = .true.
    case default ; is_x = .false. ; is_y = .false.
  end select

end subroutine categorize_axis_from_units

!> Tries to determine whether the axis name is commonly associated with an x- or y- axis.  This
!! approach is fragile and unreliable, but it a backup to reading a CARTESIAN file attribute.
subroutine categorize_axis_from_name(dimname, is_x, is_y)
  character(len=*), intent(in) :: dimname !< A dimension name
  logical, intent(out) :: is_x !< Indicates if the axis name is associated with an x-direction axis
  logical, intent(out) :: is_y !< Indicates if the axis name is associated with an y-direction axis

  is_x = .false. ; is_y = .false.
  select case(trim(lowercase(dimname)))
    case ("grid_x_t")  ; is_x = .true.
    case ("nx")        ; is_x = .true.
    case ("nxp")       ; is_x = .true.
    case ("longitude") ; is_x = .true.
    case ("long")      ; is_x = .true.
    case ("lon")       ; is_x = .true.
    case ("lonh")      ; is_x = .true.
    case ("lonq")      ; is_x = .true.
    case ("xh")        ; is_x = .true.
    case ("xq")        ; is_x = .true.
    case ("i")         ; is_x = .true.

    case ("grid_y_t")  ; is_y = .true.
    case ("ny")        ; is_y = .true.
    case ("nyp")       ; is_y = .true.
    case ("latitude")  ; is_y = .true.
    case ("lat")       ; is_y = .true.
    case ("lath")      ; is_y = .true.
    case ("latq")      ; is_y = .true.
    case ("yh")        ; is_y = .true.
    case ("yq")        ; is_y = .true.
    case ("j")         ; is_y = .true.

    case default ; is_x = .false. ; is_y = .false.
  end select

end subroutine categorize_axis_from_name


!> Write a 4d field to an output file.
subroutine write_field_4d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, fill_value)
  type(file_type),          intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),          intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),    intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:,:,:), intent(inout) :: field      !< Field to write
  real,           optional, intent(in)    :: tstamp     !< Model time of this field
  integer,        optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,           optional, intent(in)    :: fill_value !< Missing data fill value


  ! Local variables
  integer :: time_index

  if (IO_handle%FMS2_file .and. present(tstamp)) then
    time_index = write_time_if_later(IO_handle, tstamp)
    call write_data(IO_handle%fileobj, trim(field_md%name), field, unlim_dim_level=time_index)
  elseif (IO_handle%FMS2_file) then
    call write_data(IO_handle%fileobj, trim(field_md%name), field)
  else
    call mpp_write(IO_handle%unit, field_md%FT, MOM_domain%mpp_domain, field, tstamp=tstamp, &
                   tile_count=tile_count, default_data=fill_value)
  endif
end subroutine write_field_4d

!> Write a 3d field to an output file.
subroutine write_field_3d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, fill_value)
  type(file_type),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),  intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:,:), intent(inout) :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model time of this field
  integer,      optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,         optional, intent(in)    :: fill_value !< Missing data fill value

  ! Local variables
  integer :: time_index

  if (IO_handle%FMS2_file .and. present(tstamp)) then
    time_index = write_time_if_later(IO_handle, tstamp)
    call write_data(IO_handle%fileobj, trim(field_md%name), field, unlim_dim_level=time_index)
  elseif (IO_handle%FMS2_file) then
    call write_data(IO_handle%fileobj, trim(field_md%name), field)
  else
    call mpp_write(IO_handle%unit, field_md%FT, MOM_domain%mpp_domain, field, tstamp=tstamp, &
                   tile_count=tile_count, default_data=fill_value)
  endif
end subroutine write_field_3d

!> Write a 2d field to an output file.
subroutine write_field_2d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, fill_value)
  type(file_type),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),  intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:),   intent(inout) :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model time of this field
  integer,      optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,         optional, intent(in)    :: fill_value !< Missing data fill value

  ! Local variables
  integer :: time_index

  if (IO_handle%FMS2_file .and. present(tstamp)) then
    time_index = write_time_if_later(IO_handle, tstamp)
    call write_data(IO_handle%fileobj, trim(field_md%name), field, unlim_dim_level=time_index)
  elseif (IO_handle%FMS2_file) then
    call write_data(IO_handle%fileobj, trim(field_md%name), field)
  else
    call mpp_write(IO_handle%unit, field_md%FT, MOM_domain%mpp_domain, field, tstamp=tstamp, &
                   tile_count=tile_count, default_data=fill_value)
  endif
end subroutine write_field_2d

!> Write a 1d field to an output file.
subroutine write_field_1d(IO_handle, field_md, field, tstamp)
  type(file_type),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  real, dimension(:),     intent(in)    :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model time of this field

  ! Local variables
  integer :: time_index

  if (IO_handle%FMS2_file .and. present(tstamp)) then
    time_index = write_time_if_later(IO_handle, tstamp)
    call write_data(IO_handle%fileobj, trim(field_md%name), field, unlim_dim_level=time_index)
  elseif (IO_handle%FMS2_file) then
    call write_data(IO_handle%fileobj, trim(field_md%name), field)
  else
    call mpp_write(IO_handle%unit, field_md%FT, field, tstamp=tstamp)
  endif
end subroutine write_field_1d

!> Write a 0d field to an output file.
subroutine write_field_0d(IO_handle, field_md, field, tstamp)
  type(file_type),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  real,                   intent(in)    :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model time of this field

  ! Local variables
  integer :: time_index

  if (IO_handle%FMS2_file .and. present(tstamp)) then
    time_index = write_time_if_later(IO_handle, tstamp)
    call write_data(IO_handle%fileobj, trim(field_md%name), field, unlim_dim_level=time_index)
  elseif (IO_handle%FMS2_file) then
    call write_data(IO_handle%fileobj, trim(field_md%name), field)
  else
    call mpp_write(IO_handle%unit, field_md%FT, field, tstamp=tstamp)
  endif
end subroutine write_field_0d

!> Returns the integer time index for a write in this file, also writing the time variable to
!! the file if this time is later than what is already in the file.
integer function write_time_if_later(IO_handle, field_time)
  type(file_type), intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  real,            intent(in)    :: field_time !< Model time of this field

  ! Local variables
  character(len=256) :: dim_unlim_name ! name of the unlimited dimension in the file

  if ((field_time > IO_handle%file_time) .or. (IO_handle%num_times == 0)) then
    IO_handle%file_time = field_time
    IO_handle%num_times = IO_handle%num_times + 1
    if (IO_handle%FMS2_file) then
      call get_unlimited_dimension_name(IO_handle%fileobj, dim_unlim_name)
      call write_data(IO_handle%fileobj, trim(dim_unlim_name), (/field_time/), &
                      corner=(/IO_handle%num_times/), edge_lengths=(/1/))
    endif
  endif

  write_time_if_later = IO_handle%num_times
end function write_time_if_later

!> Write the data for an axis
subroutine MOM_write_axis(IO_handle, axis)
  type(file_type), intent(in) :: IO_handle  !< Handle for a file that is open for writing
  type(axistype),  intent(in) :: axis       !< An axis type variable with information to write

  integer :: is, ie

  if (IO_handle%FMS2_file) then
    if (axis%domain_decomposed) then
      ! FMS2 does not domain-decompose 1d arrays, so we explicitly slice it
      call get_global_io_domain_indices(IO_handle%fileobj, trim(axis%name), is, ie)
      call write_data(IO_handle%fileobj, trim(axis%name), axis%ax_data(is:ie))
    else
      call write_data(IO_handle%fileobj, trim(axis%name), axis%ax_data)
    endif
  else
    call mpp_write(IO_handle%unit, axis%AT)
  endif

end subroutine MOM_write_axis

!> Store information about an axis in a previously defined axistype and write this
!! information to the file indicated by unit.
subroutine write_metadata_axis(IO_handle, axis, name, units, longname, cartesian, sense, domain, data, &
                               edge_axis, calendar)
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

  character(len=:), allocatable :: cart ! A left-adjusted and trimmed copy of cartesian
  logical :: is_x, is_y, is_t  ! If true, this is a domain-decomposed axis in one of the directions.
  integer :: position    ! A flag indicating the axis staggering position.
  integer :: i, isc, iec, global_size

  if (IO_handle%FMS2_file) then
    if (is_dimension_registered(IO_handle%fileobj, trim(name))) then
      call MOM_error(FATAL, "write_metadata_axis was called more than once for axis "//trim(name)//&
                            " in file "//trim(IO_handle%filename))
      return
    endif
  endif

  axis%name = trim(name)
  if (present(data) .and. allocated(axis%ax_data)) call MOM_error(FATAL, &
        "Data is already allocated in a call to write_metadata_axis for axis "//&
        trim(name)//" in file "//trim(IO_handle%filename))

  if (IO_handle%FMS2_file) then
    is_x = .false. ; is_y = .false. ; is_t = .false.
    position = CENTER
    if (present(cartesian)) then
      cart = trim(adjustl(cartesian))
      if ((index(cart, "X") == 1) .or. (index(cart, "x") == 1)) is_x = .true.
      if ((index(cart, "Y") == 1) .or. (index(cart, "y") == 1)) is_y = .true.
      if ((index(cart, "T") == 1) .or. (index(cart, "t") == 1)) is_t = .true.
    endif

    ! For now, we assume that all horizontal axes are domain-decomposed.
    if (is_x .or. is_y) &
      axis%domain_decomposed = .true.

    if (is_x) then
      if (present(edge_axis)) then ; if (edge_axis) position = EAST_FACE ; endif
      call register_axis(IO_handle%fileobj, trim(name), 'x', domain_position=position)
    elseif (is_y) then
      if (present(edge_axis)) then ; if (edge_axis) position = NORTH_FACE ; endif
      call register_axis(IO_handle%fileobj, trim(name), 'y', domain_position=position)
    elseif (is_t .and. .not.present(data)) then
      ! This is the unlimited (time) dimension.
      call register_axis(IO_handle%fileobj, trim(name), unlimited)
    else
      if (.not.present(data)) call MOM_error(FATAL,"MOM_io:register_diagnostic_axis: "//&
                        "An axis_length argument is required to register the axis "//trim(name))
      call register_axis(IO_handle%fileobj, trim(name), size(data))
    endif

    if (present(data)) then
      ! With FMS2, the data for the axis labels has to match the computational domain on this PE.
      if (present(domain)) then
        ! The commented-out code on the next ~11 lines runs but there is missing data in the output file
        ! call mpp_get_compute_domain(domain, isc, iec)
        ! call mpp_get_global_domain(domain, size=global_size)
        ! if (size(data) == global_size) then
        !   allocate(axis%ax_data(iec+1-isc)) ; axis%ax_data(:) = data(isc:iec)
        !   ! A simpler set of labels: do i=1,iec-isc ; axis%ax_data(i) = real(isc + i) - 1.0 ; enddo
        ! elseif (size(data) == global_size+1) then
        !   ! This is an edge axis.  Note the effective SW indexing convention here.
        !   allocate(axis%ax_data(iec+2-isc)) ; axis%ax_data(:) = data(isc:iec+1)
        !   ! A simpler set of labels: do i=1,iec+1-isc ; axis%ax_data(i) = real(isc + i) - 1.5 ; enddo
        ! else
        !   call MOM_error(FATAL, "Unexpected size of data for "//trim(name)//" in write_metadata_axis.")
        ! endif

        ! This works for a simple 1x1 IO layout, but gives errors for nontrivial IO layouts
        allocate(axis%ax_data(size(data))) ; axis%ax_data(:) = data(:)

      else  ! Store the entire array of axis labels.
        allocate(axis%ax_data(size(data))) ; axis%ax_data(:) = data(:)
      endif
    endif


    ! Now create the variable that describes this axis.
    call register_field(IO_handle%fileobj, trim(name), "double", dimensions=(/name/))
    if (len_trim(longname) > 0) &
      call register_variable_attribute(IO_handle%fileobj, trim(name), 'long_name', &
                                       trim(longname), len_trim(longname))
    if (len_trim(units) > 0) &
      call register_variable_attribute(IO_handle%fileobj, trim(name), 'units', &
                                       trim(units), len_trim(units))
    if (present(cartesian)) &
      call register_variable_attribute(IO_handle%fileobj, trim(name), 'cartesian_axis', &
                                       trim(cartesian), len_trim(cartesian))
    if (present(sense)) &
      call register_variable_attribute(IO_handle%fileobj, trim(name), 'sense', sense)
  else
    if (present(data)) then
      allocate(axis%ax_data(size(data))) ; axis%ax_data(:) = data(:)
    endif

    call mpp_write_meta(IO_handle%unit, axis%AT, name, units, longname, cartesian=cartesian, sense=sense, &
                        domain=domain, data=data, calendar=calendar)
  endif
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

  ! Local variables
  character(len=256), dimension(size(axes)) :: dim_names ! The names of the dimensions
  type(mpp_axistype), dimension(size(axes)) :: mpp_axes  ! The array of mpp_axistypes for this variable
  character(len=16) :: prec_string     ! A string specifying the precision with which to save this variable
  character(len=64) :: checksum_string ! checksum character array created from checksum argument
  integer :: i, ndims

  ndims = size(axes)
  if (IO_handle%FMS2_file) then
    do i=1,ndims ; dim_names(i) = trim(axes(i)%name) ; enddo
    prec_string = "double" ; if (present(pack)) then ; if (pack > 1) prec_string = "float" ; endif
    call register_field(IO_handle%fileobj, trim(name), trim(prec_string), dimensions=dim_names)
    if (len_trim(longname) > 0) &
      call register_variable_attribute(IO_handle%fileobj, trim(name), 'long_name', &
                                       trim(longname), len_trim(longname))
    if (len_trim(units) > 0) &
      call register_variable_attribute(IO_handle%fileobj, trim(name), 'units', &
                                       trim(units), len_trim(units))
    if (present(standard_name)) &
      call register_variable_attribute(IO_handle%fileobj, trim(name), 'standard_name', &
                                       trim(standard_name), len_trim(standard_name))
    if (present(checksum)) then
      write (checksum_string,'(Z16)') checksum(1) ! Z16 is the hexadecimal format code
      call register_variable_attribute(IO_handle%fileobj, trim(name), "checksum", &
                                       trim(checksum_string), len_trim(checksum_string))
    endif
  else
    do i=1,ndims ; mpp_axes(i) = axes(i)%AT ; enddo
    call mpp_write_meta(IO_handle%unit, field%FT, mpp_axes, name, units, longname, &
                        pack=pack, standard_name=standard_name, checksum=checksum)
    ! unused opt. args: min=min, max=max, fill=fill, scale=scale, add=add, &
  endif

  ! Store information in the field-type, regardless of which interfaces are used.
  field%name = trim(name)
  field%longname = trim(longname)
  field%units = trim(units)
  field%chksum_read = -1
  field%valid_chksum = .false.

end subroutine write_metadata_field

!> Write a global text attribute to a file.
subroutine write_metadata_global(IO_handle, name, attribute)
  type(file_type),            intent(in)    :: IO_handle !< Handle for a file that is open for writing
  character(len=*),           intent(in)    :: name      !< The name in the file of this global attribute
  character(len=*),           intent(in)    :: attribute !< The value of this attribute

  if (IO_handle%FMS2_file) then
    call register_global_attribute(IO_handle%fileobj, name, attribute, len_trim(attribute))
  else
    call mpp_write_meta(IO_handle%unit, name, cval=attribute)
  endif

end subroutine write_metadata_global

end module MOM_io_infra
