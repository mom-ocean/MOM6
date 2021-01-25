!> This module contains a thin inteface to mpp and fms I/O code
module MOM_io_infra

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domain_infra,     only : MOM_domain_type, AGRID, BGRID_NE, CGRID_NE
use MOM_domain_infra,     only : get_simple_array_i_ind, get_simple_array_j_ind
use MOM_domain_infra,     only : domain2d, CENTER, CORNER, NORTH_FACE, EAST_FACE
use MOM_error_infra,      only : MOM_error=>MOM_err, NOTE, FATAL, WARNING

use ensemble_manager_mod, only : get_ensemble_id
use fms_mod,              only : write_version_number, open_namelist_file, check_nml_error
use fms_io_mod,           only : file_exist, field_exist, field_size, read_data
use fms_io_mod,           only : fms_io_exit, get_filename_appendix
use mpp_io_mod,           only : mpp_open, mpp_close, mpp_flush
use mpp_io_mod,           only : write_metadata=>mpp_write_meta, mpp_write
use mpp_io_mod,           only : mpp_get_atts, mpp_attribute_exist
use mpp_io_mod,           only : mpp_get_axes, axistype, get_axis_data=>mpp_get_axis_data
use mpp_io_mod,           only : mpp_get_fields, fieldtype
use mpp_io_mod,           only : mpp_get_info
use mpp_io_mod,           only : get_file_times=>mpp_get_times
use mpp_io_mod,           only : mpp_io_init
! These are encoding constants.
use mpp_io_mod,           only : APPEND_FILE=>MPP_APPEND, ASCII_FILE=>MPP_ASCII
use mpp_io_mod,           only : MULTIPLE=>MPP_MULTI, NETCDF_FILE=>MPP_NETCDF
use mpp_io_mod,           only : OVERWRITE_FILE=>MPP_OVERWR, READONLY_FILE=>MPP_RDONLY
use mpp_io_mod,           only : SINGLE_FILE=>MPP_SINGLE, WRITEONLY_FILE=>MPP_WRONLY

implicit none ; private

! These interfaces are actually implemented or have explicit interfaces in this file.
public :: MOM_read_data, MOM_read_vector, write_field, open_file, close_file, flush_file
public :: file_exists, field_exists, read_field_chksum
public :: get_file_info, get_file_fields, get_field_atts, io_infra_init, io_infra_end
! The following are simple pass throughs of routines from other modules.  They need
! to have explicit interfaces added to this file.
public :: fieldtype, axistype, field_size, get_filename_appendix
public :: get_file_times, read_data, get_axis_data
public :: write_metadata, write_version_number, get_ensemble_id
public :: open_namelist_file, check_nml_error
! These are encoding constants.
public :: APPEND_FILE, ASCII_FILE, MULTIPLE, NETCDF_FILE, OVERWRITE_FILE
public :: READONLY_FILE, SINGLE_FILE, WRITEONLY_FILE
public :: CENTER, CORNER, NORTH_FACE, EAST_FACE

!> Indicate whether a file exists, perhaps with domain decomposition
interface file_exists
  module procedure FMS_file_exists
  module procedure MOM_file_exists
end interface

!> Read a data field from a file
interface MOM_read_data
  module procedure MOM_read_data_4d
  module procedure MOM_read_data_3d
  module procedure MOM_read_data_2d
  module procedure MOM_read_data_1d
  module procedure MOM_read_data_0d
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
end interface

contains

!> Reads the checksum value for a field that was recorded in a file, along with a flag indicating
!! whether the file contained a valid checksum for this field.
subroutine read_field_chksum(field, chksum, valid_chksum)
  type(fieldtype), intent(in)  :: field !< The field whose checksum attribute is to be read.
  integer(kind=8), intent(out) :: chksum !< The checksum for the field.
  logical,         intent(out) :: valid_chksum  !< If true, chksum has been successfully read.
  ! Local variables
  integer(kind=8), dimension(3) :: checksum_file

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
function MOM_file_exists(filename, MOM_Domain)
  character(len=*),       intent(in) :: filename   !< The name of the file being inquired about
  type(MOM_domain_type),  intent(in) :: MOM_Domain !< The MOM_Domain that describes the decomposition

! This function uses the fms_io function file_exist to determine whether
! a named file (or its decomposed variant) exists.

  logical :: MOM_file_exists

  MOM_file_exists = file_exist(filename, MOM_Domain%mpp_domain)

end function MOM_file_exists

!> Returns true if the named file or its domain-decomposed variant exists.
function FMS_file_exists(filename, domain, no_domain)
  character(len=*),         intent(in) :: filename  !< The name of the file being inquired about
  type(domain2d), optional, intent(in) :: domain    !< The mpp domain2d that describes the decomposition
  logical,        optional, intent(in) :: no_domain !< This file does not use domain decomposition
! This function uses the fms_io function file_exist to determine whether
! a named file (or its decomposed variant) exists.

  logical :: FMS_file_exists

  FMS_file_exists = file_exist(filename, domain, no_domain)

end function FMS_file_exists

!> close_file closes a file (or fileset).  If the file handle does not point to an open file,
!! close_file simply returns without doing anything.
subroutine close_file(unit)
  integer,                  intent(out) :: unit   !< The I/O unit for the file to be closed

  call mpp_close(unit)
end subroutine close_file

!> Ensure that the output stream associated with a unit is fully sent to dis.
subroutine flush_file(unit)
  integer,                  intent(out) :: unit   !< The I/O unit for the file to flush

  call mpp_flush(unit)
end subroutine flush_file

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


!> open_file opens a file for parallel or single-file I/O.
subroutine open_file(unit, file, action, form, threading, fileset, nohdrs, domain, MOM_domain)
  integer,                  intent(out) :: unit   !< The I/O unit for the opened file
  character(len=*),         intent(in)  :: file   !< The name of the file being opened
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
    call mpp_open(unit, file, action=action, form=form, threading=threading, fileset=fileset, &
                  nohdrs=nohdrs, domain=MOM_Domain%mpp_domain)
  else
    call mpp_open(unit, file, action=action, form=form, threading=threading, fileset=fileset, &
                  nohdrs=nohdrs, domain=domain)
  endif
end subroutine open_file

!> Get information about the number of dimensions, variables, global attributes and time levels
!! in the file associated with an open file unit
subroutine get_file_info(unit, ndim, nvar, natt, ntime)
  integer,            intent(in)  :: unit  !< The I/O unit for the open file
  integer,  optional, intent(out) :: ndim  !< The number of dimensions in the file
  integer,  optional, intent(out) :: nvar  !< The number of variables in the file
  integer,  optional, intent(out) :: natt  !< The number of global attributes in the file
  integer,  optional, intent(out) :: ntime !< The number of time levels in the file

  ! Local variables
  integer :: ndims, nvars, natts, ntimes

  call mpp_get_info( unit, ndims, nvars, natts, ntimes )

  if (present(ndim)) ndim = ndims
  if (present(nvar)) nvar = nvars
  if (present(natt)) natt = natts
  if (present(ntime)) ntime = ntimes

end subroutine get_file_info

!> Set up the field information (e.g., names and metadata) for all of the variables in a file.  The
!! argument fields must be allocated with a size that matches the number of variables in a file.
subroutine get_file_fields(unit, fields)
  integer,                       intent(in)    :: unit   !< The I/O unit for the open file
  type(fieldtype), dimension(:), intent(inout) :: fields !< Field-type descriptions of all of
                                                         !! the variables in a file.
  call mpp_get_fields(unit, fields)
end subroutine get_file_fields

!> Extract information from a field type, as stored or as found in a file
subroutine get_field_atts(field, name, units, longname, checksum)
  type(fieldtype),            intent(in)  :: field !< The field to extract information from
  character(len=*), optional, intent(out) :: name  !< The variable name
  character(len=*), optional, intent(out) :: units !< The units of the variable
  character(len=*), optional, intent(out) :: longname  !< The long name of the variable
  integer(kind=8),  dimension(:), &
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

!> This function uses the fms_io function read_data to read a scalar
!! data field named "fieldname" from file "filename".
subroutine MOM_read_data_0d(filename, fieldname, data, timelevel, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real,                   intent(inout) :: data      !< The 1-dimensional array into which the data
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.

  call read_data(filename, fieldname, data, timelevel=timelevel, no_domain=.true.)

  if (present(scale)) then ; if (scale /= 1.0) then
    data = scale*data
  endif ; endif

end subroutine MOM_read_data_0d

!> This function uses the fms_io function read_data to read a 1-D
!! data field named "fieldname" from file "filename".
subroutine MOM_read_data_1d(filename, fieldname, data, timelevel, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:),     intent(inout) :: data      !< The 1-dimensional array into which the data
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before they are returned.

  call read_data(filename, fieldname, data, timelevel=timelevel, no_domain=.true.)

  if (present(scale)) then ; if (scale /= 1.0) then
    data(:) = scale*data(:)
  endif ; endif

end subroutine MOM_read_data_1d

!> This function uses the fms_io function read_data to read a distributed
!! 2-D data field named "fieldname" from file "filename".  Valid values for
!! "position" include CORNER, CENTER, EAST_FACE and NORTH_FACE.
subroutine MOM_read_data_2d(filename, fieldname, data, MOM_Domain, &
                            timelevel, position, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:),   intent(inout) :: data      !< The 2-dimensional array into which the data
                                                     !! should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: position  !< A flag indicating where this data is located
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.

  integer :: is, ie, js, je

  call read_data(filename, fieldname, data, MOM_Domain%mpp_domain, &
                 timelevel=timelevel, position=position)

  if (present(scale)) then ; if (scale /= 1.0) then
    call get_simple_array_i_ind(MOM_Domain, size(data,1), is, ie)
    call get_simple_array_j_ind(MOM_Domain, size(data,2), js, je)
    data(is:ie,js:je) = scale*data(is:ie,js:je)
  endif ; endif

end subroutine MOM_read_data_2d

!> This function uses the fms_io function read_data to read a distributed
!! 3-D data field named "fieldname" from file "filename".  Valid values for
!! "position" include CORNER, CENTER, EAST_FACE and NORTH_FACE.
subroutine MOM_read_data_3d(filename, fieldname, data, MOM_Domain, &
                            timelevel, position, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:), intent(inout) :: data      !< The 3-dimensional array into which the data
                                                     !! should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: position  !< A flag indicating where this data is located
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.

  integer :: is, ie, js, je

  call read_data(filename, fieldname, data, MOM_Domain%mpp_domain, &
                 timelevel=timelevel, position=position)

  if (present(scale)) then ; if (scale /= 1.0) then
    call get_simple_array_i_ind(MOM_Domain, size(data,1), is, ie)
    call get_simple_array_j_ind(MOM_Domain, size(data,2), js, je)
    data(is:ie,js:je,:) = scale*data(is:ie,js:je,:)
  endif ; endif

end subroutine MOM_read_data_3d

!> This function uses the fms_io function read_data to read a distributed
!! 4-D data field named "fieldname" from file "filename".  Valid values for
!! "position" include CORNER, CENTER, EAST_FACE and NORTH_FACE.
subroutine MOM_read_data_4d(filename, fieldname, data, MOM_Domain, &
                            timelevel, position, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:,:), intent(inout) :: data    !< The 4-dimensional array into which the data
                                                     !! should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: position  !< A flag indicating where this data is located
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before it is returned.

  integer :: is, ie, js, je

  call read_data(filename, fieldname, data, MOM_Domain%mpp_domain, &
                 timelevel=timelevel, position=position)

  if (present(scale)) then ; if (scale /= 1.0) then
    call get_simple_array_i_ind(MOM_Domain, size(data,1), is, ie)
    call get_simple_array_j_ind(MOM_Domain, size(data,2), js, je)
    data(is:ie,js:je,:,:) = scale*data(is:ie,js:je,:,:)
  endif ; endif

end subroutine MOM_read_data_4d


!> This function uses the fms_io function read_data to read a pair of distributed
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
  integer :: is, ie, js, je
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
    call get_simple_array_i_ind(MOM_Domain, size(u_data,1), is, ie)
    call get_simple_array_j_ind(MOM_Domain, size(u_data,2), js, je)
    u_data(is:ie,js:je) = scale*u_data(is:ie,js:je)
    call get_simple_array_i_ind(MOM_Domain, size(v_data,1), is, ie)
    call get_simple_array_j_ind(MOM_Domain, size(v_data,2), js, je)
    v_data(is:ie,js:je) = scale*v_data(is:ie,js:je)
  endif ; endif

end subroutine MOM_read_vector_2d

!> This function uses the fms_io function read_data to read a pair of distributed
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

  integer :: is, ie, js, je
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
    call get_simple_array_i_ind(MOM_Domain, size(u_data,1), is, ie)
    call get_simple_array_j_ind(MOM_Domain, size(u_data,2), js, je)
    u_data(is:ie,js:je,:) = scale*u_data(is:ie,js:je,:)
    call get_simple_array_i_ind(MOM_Domain, size(v_data,1), is, ie)
    call get_simple_array_j_ind(MOM_Domain, size(v_data,2), js, je)
    v_data(is:ie,js:je,:) = scale*v_data(is:ie,js:je,:)
  endif ; endif

end subroutine MOM_read_vector_3d


!> Write a 4d field to an output file.
subroutine write_field_4d(io_unit, field_md, MOM_domain, field, tstamp, tile_count, fill_value)
  integer,                  intent(in)    :: io_unit    !< File I/O unit handle
  type(fieldtype),          intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),    intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:,:,:), intent(inout) :: field      !< Field to write
  real,           optional, intent(in)    :: tstamp     !< Model timestamp
  integer,        optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,           optional, intent(in)    :: fill_value !< Missing data fill value

  call mpp_write(io_unit, field_md, MOM_domain%mpp_domain, field, tstamp=tstamp, &
                 tile_count=tile_count, default_data=fill_value)
end subroutine write_field_4d

!> Write a 3d field to an output file.
subroutine write_field_3d(io_unit, field_md, MOM_domain, field, tstamp, tile_count, fill_value)
  integer,                intent(in)    :: io_unit    !< File I/O unit handle
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),  intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:,:), intent(inout) :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model timestamp
  integer,      optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,         optional, intent(in)    :: fill_value !< Missing data fill value

  call mpp_write(io_unit, field_md, MOM_domain%mpp_domain, field, tstamp=tstamp, &
                   tile_count=tile_count, default_data=fill_value)
end subroutine write_field_3d

!> Write a 2d field to an output file.
subroutine write_field_2d(io_unit, field_md, MOM_domain, field, tstamp, tile_count, fill_value)
  integer,                intent(in)    :: io_unit    !< File I/O unit handle
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),  intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:),   intent(inout) :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model timestamp
  integer,      optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,         optional, intent(in)    :: fill_value !< Missing data fill value

  call mpp_write(io_unit, field_md, MOM_domain%mpp_domain, field, tstamp=tstamp, &
                   tile_count=tile_count, default_data=fill_value)
end subroutine write_field_2d

!> Write a 1d field to an output file.
subroutine write_field_1d(io_unit, field_md, field, tstamp)
  integer,                intent(in)    :: io_unit    !< File I/O unit handle
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  real, dimension(:),     intent(in)    :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model timestamp

  call mpp_write(io_unit, field_md, field, tstamp=tstamp)
end subroutine write_field_1d

!> Write a 0d field to an output file.
subroutine write_field_0d(io_unit, field_md, field, tstamp)
  integer,                intent(in)    :: io_unit    !< File I/O unit handle
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  real,                   intent(in)    :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model timestamp

  call mpp_write(io_unit, field_md, field, tstamp=tstamp)
end subroutine write_field_0d

subroutine MOM_write_axis(io_unit, axis)
  integer,        intent(in)  :: io_unit    !< File I/O unit handle
  type(axistype), intent(in)  :: axis       !< An axis type variable with information to write

  call mpp_write(io_unit, axis)

end subroutine MOM_write_axis

end module MOM_io_infra
