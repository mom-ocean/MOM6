!> This module contains I/O framework code
module MOM_io

! This file is part of MOM6. See LICENSE.md for the license.


use MOM_error_handler,    only : MOM_error, NOTE, FATAL, WARNING
use MOM_domains,          only : MOM_domain_type, AGRID, BGRID_NE, CGRID_NE
use MOM_domains,          only : get_simple_array_i_ind, get_simple_array_j_ind
use MOM_file_parser,      only : log_version, param_file_type
use MOM_grid,             only : ocean_grid_type
use MOM_dyn_horgrid,      only : dyn_horgrid_type
use MOM_string_functions, only : lowercase, slasher
use MOM_string_functions, only : append_substring
use MOM_time_manager,     only : time_type, time_type_to_real
use MOM_verticalGrid,     only : verticalGrid_type
use ensemble_manager_mod, only : get_ensemble_id
use fms_mod,              only : write_version_number, open_namelist_file, check_nml_error
use MOM_string_functions,  only : extract_word
use mpp_mod,              only : mpp_max 
use mpp_domains_mod,      only : domain1d, domain2d, domainug, mpp_get_domain_components
use mpp_domains_mod,      only : CENTER, CORNER, NORTH_FACE=>NORTH, EAST_FACE=>EAST
use mpp_io_mod,           only : mpp_open_file => mpp_open, mpp_close_file => mpp_close
use mpp_io_mod,           only : mpp_write_meta
use mpp_io_mod,           only : axistype
use mpp_io_mod,           only : fieldtype, axistype, flush_file => mpp_flush
use mpp_io_mod,           only : APPEND_FILE=>MPP_APPEND, ASCII_FILE=>MPP_ASCII
use mpp_io_mod,           only : MULTIPLE=>MPP_MULTI, NETCDF_FILE=>MPP_NETCDF
use mpp_io_mod,           only : OVERWRITE_FILE=>MPP_OVERWR, READONLY_FILE=>MPP_RDONLY
use mpp_io_mod,           only : SINGLE_FILE=>MPP_SINGLE, WRITEONLY_FILE=>MPP_WRONLY
use mpp_io_mod,           only : MPP_APPEND, MPP_MULTI, MPP_OVERWR, MPP_NETCDF, MPP_RDONLY
use mpp_io_mod,           only : io_infra_init=>mpp_io_init

use fms2_io_mod,          only: check_if_open, &
                                get_dimension_names, &           
                                get_dimension_size, &                      
                                get_compute_domain_dimension_indices, &
                                get_global_attribute, &
                                get_global_io_domain_indices, &
                                get_num_dimensions, &
                                get_num_variables, &
                                get_variable_dimension_names, &
                                get_variable_num_dimensions, &
                                get_variable_size, &
                                get_variable_units, &
                                get_variable_unlimited_dimension_index, &
                                global_att_exists, &
                                is_dimension_unlimited, &
                                read_data, &
                                read_restart, &
                                register_restart_field, &
                                register_axis, &
                                register_field, &
                                register_variable_attribute, &
                                open_file, &
                                close_file, &
                                write_data, &
                                write_restart, &
                                attribute_exists => variable_att_exists, &
                                variable_exists, &
                                dimension_exists, &
                                file_exists, &
                                FmsNetcdfDomainFile_t, &
                                FmsNetcdfFile_t, & 
                                FmsNetcdfUnstructuredDomainFile_t, &
                                unlimited

use netcdf

implicit none ; private

public :: mpp_close_file, mpp_open_file, fieldtype, flush_file
public :: MOM_read_vector, ensembler, num_timelevels
public :: slasher, MOM_io_init
public :: open_namelist_file, check_nml_error, io_infra_init
public :: APPEND_FILE, ASCII_FILE, MULTIPLE, NETCDF_FILE, OVERWRITE_FILE
public :: READONLY_FILE, SINGLE_FILE, WRITEONLY_FILE
public :: CENTER, CORNER, NORTH_FACE, EAST_FACE
public :: var_desc, modify_vardesc, query_vardesc, cmor_long_std
public :: scale_data
! new FMS-IO routines and wrappers
public :: attribute_exists
public :: check_if_open
public :: close_file
public :: dimension_exists
public :: file_exists
public :: FmsNetcdfFile_t
public :: FmsNetcdfDomainFile_t
public :: FmsNetcdfUnstructuredDomainFile_t
public :: get_compute_domain_dimension_indices
public :: get_dimension_names
public :: get_dimension_size
public :: get_global_io_domain_indices
public :: get_global_attribute
public :: get_horizontal_grid_logic
public :: get_num_dimensions
public :: get_num_variables
public :: get_time_units
public :: get_var_dimension_features
public :: get_variable_dimension_names
public :: get_variable_byte_size
public :: get_variable_num_dimensions
public :: get_variable_size
public :: get_variable_units
public :: get_variable_unlimited_dimension_index
public :: global_att_exists
public :: is_dimension_unlimited
public :: MOM_get_diagnostic_axis_data
public :: MOM_get_nc_corner_edgelengths
public :: MOM_open_file
public :: MOM_read_data
public :: MOM_register_diagnostic_axis
public :: MOM_register_variable_axes
public :: open_file
public :: read_data
public :: read_restart
public :: register_axis
public :: register_field
public :: register_restart_field
public :: register_variable_attribute
public :: variable_exists
public :: write_data
public :: write_restart
public :: unlimited

!> Type for describing a variable, typically a tracer
type, public :: vardesc
  character(len=64)  :: name               !< Variable name in a NetCDF file
  character(len=48)  :: units              !< Physical dimensions of the variable
  character(len=240) :: longname           !< Long name of the variable
  character(len=8)   :: hor_grid           !< Horizontal grid:  u, v, h, q, Cu, Cv, T, Bu, or 1
  character(len=8)   :: z_grid             !< Vertical grid:  L, i, or 1
  character(len=8)   :: t_grid             !< Time description: s, p, or 1
  character(len=64)  :: cmor_field_name    !< CMOR name
  character(len=64)  :: cmor_units         !< CMOR physical dimensions of the variable
  character(len=240) :: cmor_longname      !< CMOR long name of the variable
  real               :: conversion         !< for unit conversions, such as needed to
                                           !! convert from intensive to extensive
end type vardesc

!> A type for making arrays of pointers to 1-d arrays
type p1d
  real, dimension(:), pointer :: p => NULL() !< A pointer to a 1d array
end type p1d

!> A structure with information about a single axis variable
type axis_atts
  character(len=64)  :: name                    !< Names of the axis
  character(len=48)  :: units                   !< Physical dimensions of the axis
  character(len=240) :: longname                !< Long name of the axis
  character(len=8)   :: positive                !< Positive-definite direction: up, down, east, west, north, south
  integer            :: horgrid_position        !< Horizontal grid position
  logical            :: is_domain_decomposed    !< if .true. the axis data are domain-decomposed
                                                !! and need to be indexed by the compute domain
                                                !! before passing to write_data
end type axis_atts

!> Type for describing an axis variable (e.g., lath, lonh, Time)
type, public :: axis_data_type
  !> An array of descriptions of the registered axes
  type(axis_atts), pointer :: axis(:) => NULL()  !< structure with axis attributes
  type(p1d), pointer       :: data(:) => NULL()  !< pointer to the axis data
end type axis_data_type

!> interface for setting start and nread arrays to read in data from a netCDF file
interface MOM_get_nc_corner_edgelengths
  module procedure MOM_get_nc_corner_edgelengths_DD
  module procedure MOM_get_nc_corner_edgelengths_noDD

end interface 
!> Open a netCDF file 
interface MOM_open_file
  module procedure MOM_open_file_DD_ocean_grid
  module procedure MOM_open_file_DD_supergrid
  module procedure MOM_open_file_DD_dyn_horgrid
  module procedure MOM_open_file_noDD
end interface

! interface to read data from a netcdf file
interface MOM_read_data
  module procedure MOM_read_data_4d_DD
  module procedure MOM_read_data_3d_DD
  module procedure MOM_read_data_2d_DD
  module procedure MOM_read_data_1d_DD
  module procedure MOM_read_data_scalar
  module procedure MOM_read_data_4d_noDD
  module procedure MOM_read_data_3d_noDD
  module procedure MOM_read_data_2d_noDD
  module procedure MOM_read_data_1d_noDD
end interface 

!> Read a pair of data fields representing the two components of a vector from a netcdf file
interface MOM_read_vector
  module procedure MOM_read_vector_3d
  module procedure MOM_read_vector_2d
end interface

! interface to scale data after reading in a field
interface scale_data
  module procedure scale_data_4d
  module procedure scale_data_3d
  module procedure scale_data_2d
  module procedure scale_data_1d
end interface 


contains

!> Open domain-decomposed file(s) with the base file name 'filename' to read, write/append, or overwrite data.  
!! The domain comes from the ocean_grid_type structure G.
function MOM_open_file_DD_ocean_grid(MOMfileObj, filename, mode, G, is_restart) result(fileOpenSuccess)
  type(FmsNetcdfDomainFile_t), intent(inout) :: MOMfileObj !< netCDF file object 
  character(len=*),       intent(in) :: filename !< The base filename of the file(s) to search for
  character(len=*),       intent(in) :: mode !< read or write(checks if file exists to append)
  type(ocean_grid_type),      intent(in) :: G !< The ocean's grid structure
  logical, intent(in) :: is_restart !< indicates whether to check for restart file(s)
  ! local
  logical :: fileOpenSuccess ! returns .true. if the file(s) is(are) opened
  character(len=512) :: mesg ! A message for warnings.
   
  select case (trim(mode))
     case("read")
        fileOpenSuccess = open_file(MOMfileObj, filename, "read", & 
                          G%Domain%mpp_domain, is_restart = is_restart)
     case("write")
        ! check if file(s) already exists and can be appended
        fileOpenSuccess = open_file(MOMfileObj, filename, "append", & 
                                   G%Domain%mpp_domain, is_restart = is_restart)
        if (.not.(fileOpenSuccess)) then
           ! create and open new file(s) for domain-decomposed write
           fileOpenSuccess = open_file(MOMfileObj, filename, "write", & 
                                   G%Domain%mpp_domain, is_restart = is_restart)
        endif
     case("overwrite")
        ! overwrite existing file
        fileOpenSuccess = open_file(MOMfileObj, filename, "overwrite", &
                                  G%Domain%mpp_domain, is_restart = is_restart) 
     case default
        write(mesg,'( "ERROR, file mode must be read or write to open ",A)') trim(filename)
        call MOM_error(FATAL,"MOM_io::MOM_open_file_DD_ocean_grid: "//mesg)
  end select     
end function MOM_open_file_DD_ocean_grid

!> Open domain-decomposed file with the base file name 'filename' to read from, overwrite, or write/append to. 
!! The domain comes from the MOM_domain_type structure G.
function MOM_open_file_DD_supergrid(MOMfileObj, filename, mode, G, is_restart) result(fileOpenSuccess)
  type(FmsNetcdfDomainFile_t), intent(inout) :: MOMfileObj !< netCDF file object 
  character(len=*),       intent(in) :: filename !< The base filename of the file(s) to search for
  character(len=*),       intent(in) :: mode !< read or write(checks if file exists to append)
  type(MOM_domain_type),  intent(in)  :: G ! Supergrid domain defined in MOM_grid_initialize.F90
  logical, intent(in) :: is_restart !< indicates whether to check for restart file(s)
  ! local
  logical :: fileOpenSuccess ! returns .true. if the file(s) is(are) opened
  character(len=512) :: mesg ! A message for warnings.
   
  select case (trim(mode))
     case("read")
        fileOpenSuccess = open_file(MOMfileObj, filename, "read", & 
                          G%mpp_domain, is_restart = is_restart)
     case("write")
        ! check if file(s) already exists and can be appended
        fileOpenSuccess = open_file(MOMfileObj, filename, "append", & 
                                   G%mpp_domain, is_restart = is_restart)
        if (.not.(fileOpenSuccess)) then
           ! create and open new file(s) for domain-decomposed write
           fileOpenSuccess = open_file(MOMfileObj, filename, "write", & 
                                   G%mpp_domain, is_restart = is_restart)
        endif
     case("overwrite")
        ! overwrite existing file
        fileOpenSuccess = open_file(MOMfileObj, filename, "overwrite", & 
                                   G%mpp_domain, is_restart = is_restart)
     case default
        write(mesg,'( "ERROR, file mode must be read or write to open ",A)') trim(filename)
        call MOM_error(FATAL,"MOM_io::MOM_open_file_DD_supergrid: "//mesg)
  end select     
end function MOM_open_file_DD_supergrid

!> Open domain-decomposed file with the base file name 'filename' to read, overwrite, or write/append data. 
!! The domain comes from the dyn_horgrid_type structure G.
function MOM_open_file_DD_dyn_horgrid(MOMfileObj, filename, mode, G, is_restart) result(fileOpenSuccess)
  type(FmsNetcdfDomainFile_t), intent(inout) :: MOMfileObj !< netCDF file object 
  character(len=*),       intent(in) :: filename !< The base filename of the file(s) to search for
  character(len=*),       intent(in) :: mode !< read or write(checks if file exists to append)
  type(dyn_horgrid_type),  intent(in)  :: G !< Supergrid domain defined in MOM_grid_initialize.F90
  logical, intent(in) :: is_restart !< indicates whether to check for restart file(s)
  ! local
  logical :: fileOpenSuccess ! returns .true. if the file(s) is(are) opened
  character(len=512) :: mesg ! A message for warnings.
   
  select case (trim(mode))
     case("read")
        fileOpenSuccess = open_file(MOMfileObj, filename, "read", & 
                          G%Domain%mpp_domain, is_restart = is_restart)
     case("write")
        ! check if file(s) already exists and can be appended
        fileOpenSuccess = open_file(MOMfileObj, filename, "append", & 
                                   G%Domain%mpp_domain, is_restart = is_restart)
        if (.not.(fileOpenSuccess)) then
           ! create and open new file(s) for domain-decomposed write
           fileOpenSuccess = open_file(MOMfileObj, filename, "write", & 
                                   G%Domain%mpp_domain, is_restart = is_restart)
        endif
     case("overwrite")
        ! overwrite existing file
        fileOpenSuccess = open_file(MOMfileObj, filename, "overwrite", & 
                                   G%Domain%mpp_domain, is_restart = is_restart)
     case default
        write(mesg,'( "ERROR, file mode must be read or write to open ",A)') trim(filename)
        call MOM_error(FATAL,"MOM_io::MOM_open_file_DD_supergrid: "//mesg)
  end select     
end function MOM_open_file_DD_dyn_horgrid

!> Open non-domain-decomposed file(s) with the base file name 'filename' to read, overwrite, or write/append data.
function MOM_open_file_noDD(MOMfileObj, filename, mode, is_restart) result(fileOpenSuccess)
  type(FmsNetcdfFile_t), intent(inout) :: MOMfileObj !< netCDF file object 
  character(len=*),       intent(in) :: filename !< The base filename of the file(s) to search for
  character(len=*),       intent(in) :: mode !< read or write(checks if file exists to append)
  logical, intent(in) :: is_restart !< indicates whether to check for restart file(s)
  ! local
  logical :: fileOpenSuccess ! returns .true. if the file(s) is(are) opened
  character(len=512) :: mesg  ! A message for warnings.
   
  fileOpenSuccess = .false.
  select case (trim(mode))
     case("read")
        fileOpenSuccess = open_file(MOMfileObj, filename, "read", & 
                          is_restart = is_restart)
     case("write")
        ! check if file(s) already exists and can be appended
        fileOpenSuccess = open_file(MOMfileObj, filename, "append", & 
                                   is_restart = is_restart)
        if (.not.(fileOpenSuccess)) then
           ! create and open new file(s) for non-domain-decomposed write
           fileOpenSuccess = open_file(MOMfileObj, filename, "write", & 
                                   is_restart = is_restart)
        endif
     case("overwrite")
        ! overwirte existing file
        fileOpenSuccess = open_file(MOMfileObj, filename, "overwrite", & 
                                   is_restart = is_restart)
     case default
        write(mesg,'( "ERROR, file mode must be read or write to open ",A)') trim(filename)
        call MOM_error(FATAL,"MOM_io::MOM_open_file_DD_ocean_grid: "//mesg)
  end select     
end function MOM_open_file_noDD


!> This function uses the fms_io function read_data to read 1-D domain-decomposed data field named "fieldname" 
!! from file "filename".
subroutine MOM_read_data_1d_DD(filename, fieldname, data, domain, corner, edgeLengths, timeLevel, scale)
  character(len=*),       intent(in) :: filename !< The name of the file to read
  character(len=*),       intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:),     intent(inout) :: data !< The 1-dimensional data array to pass to read_data
  type(MOM_domain_type),  intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  integer,      optional, intent(in) :: corner !< starting index of data buffer. Default is 1
  integer,      optional, intent(in) :: edgeLengths !< number of data values to read in; default is the variable size
  integer,      optional, intent(in) :: timeLevel !< time level to read
  real,         optional, intent(in) :: scale !< A scaling factor that the field is multiplied by
  ! local
  type(FmsNetcdfDomainFile_t) :: fileObjRead ! netCDF file object returned by call to MOM_open_file
  logical :: fileOpenSuccess !.true. if call to MOM_open_file is successful
  integer :: i
  integer, dimension(1) :: start, nread ! indices for first data value and number of values to read
  character(len=40), dimension(1) :: dimNames ! variable dimension names
  
  ! open the file
  if (.not.(check_if_open(fileObjRead))) &
    fileOpenSuccess = open_file(fileObjRead, filename, "read", domain%mpp_domain, is_restart=.false.)
  ! register the variableaxes
  !> @note: the user will need to change the xUnits and yUnits if they expect different values for the
  !! x/longitude and/or y/latitude axes units
  call MOM_register_variable_axes(fileObjRead, trim(fieldname), xUnits="degrees_east", yUnits="degrees_north")
  
  if (present(corner) .or. present(edgeLengths) .or. present(timeLevel)) then
    call get_variable_dimension_names(fileObjRead, trim(fieldname), dimNames)
  endif
  
  start(1) = 1
  if (present(corner)) start(1) = corner
  if (present(edgeLengths)) then
    nread(1) = edgeLengths
  else
    call get_dimension_size(fileObjRead, trim(dimNames(1)), nread(1))
  endif
    
   ! read the data
  call read_data(fileObjRead, trim(fieldname), data, corner=start, edge_Lengths=nread)
  ! close the file
  if (check_if_open(fileObjRead)) call close_file(fileObjRead)

  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif

end subroutine MOM_read_data_1d_DD

!> This function uses the fms_io function read_data to read 2-D domain-decomposed data field named "fieldname" 
!! from file "filename".
subroutine MOM_read_data_2d_DD(filename, fieldname, data, domain, corner, edgeLengths, timeLevel, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:),   intent(inout) :: data !< The 2-dimensional data array to pass to read_data
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  integer, dimension(2),  optional, intent(in) :: corner !< starting indices of data buffer. Default is 1
  integer, dimension(2),  optional, intent(in) :: edgeLengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timeLevel !< time level to read
  real,         optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  ! local
  type(FmsNetcdfDomainFile_t) :: fileObjRead !netCDF file object returned by call to MOM_open_file
  logical :: fileOpenSuccess !.true. if call to MOM_open_file is successful
  integer :: i, dimUnlimIndex
  integer, dimension(2) :: start, nread ! indices for first data value and number of values to read
  character(len=40), dimension(2) :: dimNames ! variable dimension names
  
  ! open the file
  if (.not.(check_if_open(fileObjRead))) &
    fileOpenSuccess = open_file(fileObjRead, filename, "read", domain%mpp_domain, is_restart=.false.)
  ! register the variable axes
  !> @note: the user will need to change the xUnits and yUnits if they expect different values for the
  !! x/longitude and/or y/latitude axes units
  call MOM_register_variable_axes(fileObjRead, trim(fieldname), xUnits="degrees_east", yUnits="degrees_north")
  
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  if (present(corner) .or. present(edgeLengths) .or. present(timeLevel)) then
    call get_variable_dimension_names(fileObjRead, trim(fieldname), dimNames)
  endif
  
  start(:) = 1
  if (present(corner)) start = corner
  
  if (present(edgeLengths)) then
    nread = edgeLengths
  else
    do i=1,2
      call get_dimension_size(fileObjRead, trim(dimNames(i)), nread(i))
    enddo
  endif

  if (present(timeLevel)) then
    dimUnlimIndex=0
    do i=1,2
      if (is_dimension_unlimited(fileObjRead,dimNames(i))) then 
        dimUnlimIndex=i
        start(i)=timeLevel
        nread(i)=1
      endif
    enddo
    if (dimUnlimIndex .LE. 0) &
      call MOM_error(FATAL, "MOM_io::MOM_read_data_2d_DD: time level specified, but variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
  endif
 
  ! read the data
  call read_data(fileObjRead, trim(fieldname), data, corner=start, edge_Lengths=nread)
  ! close the file
  if (check_if_open(fileObjRead)) call close_file(fileObjRead)

  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale, domain)
  endif ; endif

end subroutine MOM_read_data_2d_DD

!> This function uses the fms_io function read_data to read 3-D domain-decomposed data field named "fieldname" 
!! from file "filename".
subroutine MOM_read_data_3d_DD(filename, fieldname, data, domain, corner, edgeLengths, timeLevel, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:),   intent(inout) :: data !< The 3-dimensional data array to pass to read_data
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  integer, dimension(3),  optional, intent(in) :: corner !< starting indices of data buffer. Default is 1
  integer, dimension(3),  optional, intent(in) :: edgeLengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timeLevel !< time level to read
  real,         optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  ! local
  type(FmsNetcdfDomainFile_t) :: fileObjRead !netCDF file object returned by call to MOM_open_file
  logical :: fileOpenSuccess !.true. if call to MOM_open_file is successful
  integer :: i, dimUnlimIndex
  integer, dimension(3) :: start, nread ! indices for first data value and number of values to read
  character(len=40), dimension(3) :: dimNames ! variable dimension names
  
  ! open the file
  if (.not.(check_if_open(fileObjRead))) &
    fileOpenSuccess = open_file(fileObjRead, filename, "read", domain%mpp_domain, is_restart=.false.)
  ! register the variable axes
  !> @note: the user will need to change the xUnits and yUnits if they expect different values for the
  !! x/longitude and/or y/latitude axes units
  call MOM_register_variable_axes(fileObjRead, trim(fieldname), xUnits="degrees_east", yUnits="degrees_north")
  
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  if (present(corner) .or. present(edgeLengths) .or. present(timeLevel)) then
    call get_variable_dimension_names(fileObjRead, trim(fieldname), dimNames)
  endif
  
  start(:) = 1
  if (present(corner)) start = corner
  
  if (present(edgeLengths)) then
    nread = edgeLengths
  else
    do i=1,3
      call get_dimension_size(fileObjRead, trim(dimNames(i)), nread(i))
    enddo
  endif

  if (present(timeLevel)) then
    dimUnlimIndex=0
    do i=1,3
      if (is_dimension_unlimited(fileObjRead,dimNames(i))) then 
        dimUnlimIndex=i
        start(i)=timeLevel
        nread(i)=1
      endif
    enddo
    if (dimUnlimIndex .LE. 0) &
      call MOM_error(FATAL, "MOM_io::MOM_read_data_3d_DD: time level specified, but variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
  endif
 
  ! read the data
  call read_data(fileObjRead, trim(fieldname), data, corner=start, edge_Lengths=nread)
  ! close the file
  if (check_if_open(fileObjRead)) call close_file(fileObjRead)

  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale, domain)
  endif ; endif

end subroutine MOM_read_data_3d_DD

!> This function uses the fms_io function read_data to read 4-D domain-decomposed data field named "fieldname" 
!! from file "filename".
subroutine MOM_read_data_4d_DD(filename, fieldname, data, domain, corner, edgeLengths, timeLevel, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:,:),   intent(inout) :: data !< The 4-dimensional data array to pass to read_data
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  integer, dimension(4),  optional, intent(in) :: corner !< starting indices of data buffer. Default is 1
  integer, dimension(4),  optional, intent(in) :: edgeLengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timeLevel !< time level to read
  real,         optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  ! local
  type(FmsNetcdfDomainFile_t) :: fileObjRead !netCDF file object returned by call to MOM_open_file
  logical :: fileOpenSuccess !.true. if call to MOM_open_file is successful
  integer :: i, dimUnlimIndex
  integer, dimension(4) :: start, nread ! indices for first data value and number of values to read
  character(len=40), dimension(4) :: dimNames ! variable dimension names
  
  ! open the file
  if (.not.(check_if_open(fileObjRead))) &
    fileOpenSuccess = open_file(fileObjRead, filename, "read", domain%mpp_domain, is_restart=.false.)
  ! register the variable axes
  !> @note: the user will need to change the xUnits and yUnits if they expect different values for the
  !! x/longitude and/or y/latitude axes units
  call MOM_register_variable_axes(fileObjRead, trim(fieldname), xUnits="degrees_east", yUnits="degrees_north")
  
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  if (present(corner) .or. present(edgeLengths) .or. present(timeLevel)) then
    call get_variable_dimension_names(fileObjRead, trim(fieldname), dimNames)
  endif
  
  start(:) = 1
  if (present(corner)) start = corner
  
  if (present(edgeLengths)) then
    nread = edgeLengths
  else
    do i=1,4
      call get_dimension_size(fileObjRead, trim(dimNames(i)), nread(i))
    enddo
  endif

  if (present(timeLevel)) then
    dimUnlimIndex=0
    do i=1,4
      if (is_dimension_unlimited(fileObjRead,dimNames(i))) then 
        dimUnlimIndex=i
        start(i)=timeLevel
        nread(i)=1
      endif
    enddo
    if (dimUnlimIndex .LE. 0) &
      call MOM_error(FATAL, "MOM_io::MOM_read_data_4d_DD: time level specified, but variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
  endif
 
  ! read the data
  call read_data(fileObjRead, trim(fieldname), data, corner=start, edge_Lengths=nread)
  ! close the file
  if (check_if_open(fileObjRead)) call close_file(fileObjRead)

  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale, domain)
  endif ; endif

end subroutine MOM_read_data_4d_DD

!> This function uses the fms_io function read_data to read a scalar data value named "fieldname" 
!! from file "filename".
subroutine MOM_read_data_scalar(filename, fieldname, data)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, intent(inout) :: data !< data buffer to pass to read_data
  
  ! local
  type(FmsNetcdfFile_t) :: fileObjRead ! netCDF file object returned by call to MOM_open_file
  logical :: fileOpenSuccess !.true. if call to MOM_open_file is successful
  ! open the file
  if (.not.(check_if_open(fileObjRead))) &
    fileOpenSuccess = open_file(fileObjRead, filename, "read", is_restart=.false.)
  ! read the data
  call read_data(fileObjRead, trim(fieldname), data)
  ! close the file
  if (check_if_open(fileObjRead)) call close_file(fileObjRead)

end subroutine MOM_read_data_scalar

!> This function uses the fms_io function read_data to read 1-D data field named "fieldname" 
!! from file "filename".
subroutine MOM_read_data_1d_noDD(filename, fieldname, data, corner, edgeLengths, timeLevel, scale)
  character(len=*),       intent(in) :: filename !< The name of the file to read
  character(len=*),       intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:),     intent(inout) :: data !< The 1-dimensional data array to pass to read_data
  integer,      optional, intent(in) :: corner !< starting index of data buffer. Default is 1
  integer,      optional, intent(in) :: edgeLengths !< number of data values to read in; default is the variable size
  integer,      optional, intent(in) :: timeLevel !< time level to read
  real,         optional, intent(in) :: scale !< A scaling factor that the field is multiplied by
  ! local
  type(FmsNetcdfFile_t) :: fileObjRead ! netCDF file object returned by call to MOM_open_file
  logical :: fileOpenSuccess !.true. if call to MOM_open_file is successful
  integer :: i, ndims
  integer, allocatable, dimension(:) :: start, nread ! indices for first data value and number of values to read
  character(len=40), allocatable, dimension(:) :: dimNames ! variable dimension names
  
  ! open the file
  if (.not.(check_if_open(fileObjRead))) &
    fileOpenSuccess = open_file(fileObjRead, filename, "read", is_restart=.false.)

  ndims = get_variable_num_dimensions(fileObjRead,trim(fieldname))
  allocate(start(ndims))
  allocate(nread(ndims))
  allocate(dimNames(ndims))
  
  if (present(corner)) start(1) = corner
  if (present(edgeLengths)) then
    nread(1) = edgeLengths
  else
    call get_variable_dimension_names(fileObjRead, trim(fieldname), dimNames)
    call get_dimension_size(fileObjRead, dimNames(1), nread(1))
  endif
    
  ! read the data
  call read_data(fileObjRead, trim(fieldname), data, corner=start, edge_Lengths=nread)
  ! close the file
  if (check_if_open(fileObjRead)) call close_file(fileObjRead)

  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif

  deallocate(start)
  deallocate(nread)
  deallocate(dimNames)
end subroutine MOM_read_data_1d_noDD

!> This function uses the fms_io function read_data to read 2-D data field named "fieldname" 
!! from file "filename".
subroutine MOM_read_data_2d_noDD(filename, fieldname, data, corner, edgeLengths, timeLevel, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:),   intent(inout) :: data !< The 2-dimensional data array to pass to read_data
  integer, dimension(2),  optional, intent(in) :: corner !< starting indices of data buffer. Default is 1
  integer, dimension(2),  optional, intent(in) :: edgeLengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timeLevel !< time level to read
  real,         optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  ! local
  type(FmsNetcdfFile_t) :: fileObjRead !netCDF file object returned by call to MOM_open_file
  logical :: fileOpenSuccess !.true. if call to MOM_open_file is successful
  integer :: i, dimUnlimIndex
  integer, dimension(2) :: start, nread ! indices for first data value and number of values to read
  character(len=40), dimension(2) :: dimNames ! variable dimension names
  
  ! open the file
  if (.not.(check_if_open(fileObjRead))) &
    fileOpenSuccess = open_file(fileObjRead, filename, "read", is_restart=.false.)

  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  if (present(corner) .or. present(edgeLengths) .or. present(timeLevel)) then
    call get_variable_dimension_names(fileObjRead, trim(fieldname), dimNames)
  endif
  
  start(:) = 1
  if (present(corner)) start = corner
  
  if (present(edgeLengths)) then
    nread = edgeLengths
  else
    do i=1,2
      call get_dimension_size(fileObjRead, trim(dimNames(i)), nread(i))
    enddo
  endif

  if (present(timeLevel)) then
    dimUnlimIndex=0
    do i=1,2
      if (is_dimension_unlimited(fileObjRead,dimNames(i))) then 
        dimUnlimIndex=i
        start(i)=timeLevel
        nread(i)=1
      endif
    enddo
    if (dimUnlimIndex .LE. 0) &
      call MOM_error(FATAL, "MOM_io::MOM_read_data_2d_noDD: time level specified, but variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
  endif
 
  ! read the data
  call read_data(fileObjRead, trim(fieldname), data, corner=start, edge_Lengths=nread)
  ! close the file
  if (check_if_open(fileObjRead)) call close_file(fileObjRead)

  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif

end subroutine MOM_read_data_2d_noDD

!> This function uses the fms_io function read_data to read 3-D data field named "fieldname" 
!! from file "filename".
subroutine MOM_read_data_3d_noDD(filename, fieldname, data, corner, edgeLengths, timeLevel, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:),   intent(inout) :: data !< The 3-dimensional data array to pass to read_data
  integer, dimension(3),  optional, intent(in) :: corner !< starting indices of data buffer. Default is 1
  integer, dimension(3),  optional, intent(in) :: edgeLengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timeLevel !< time level to read
  real,         optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  ! local
  type(FmsNetcdfFile_t) :: fileObjRead !netCDF file object returned by call to MOM_open_file
  logical :: fileOpenSuccess !.true. if call to MOM_open_file is successful
  integer :: i, dimUnlimIndex
  integer, dimension(3) :: start, nread ! indices for first data value and number of values to read
  character(len=40), dimension(3) :: dimNames ! variable dimension names
  
  ! open the file
  if (.not.(check_if_open(fileObjRead))) &
    fileOpenSuccess = open_file(fileObjRead, filename, "read", is_restart=.false.)

  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  if (present(corner) .or. present(edgeLengths) .or. present(timeLevel)) then
    call get_variable_dimension_names(fileObjRead, trim(fieldname), dimNames)
  endif
  
  start(:) = 1
  if (present(corner)) start = corner
  
  if (present(edgeLengths)) then
    nread = edgeLengths
  else
    do i=1,3
      call get_dimension_size(fileObjRead, trim(dimNames(i)), nread(i))
    enddo
  endif

  if (present(timeLevel)) then
    dimUnlimIndex=0
    do i=1,3
      if (is_dimension_unlimited(fileObjRead,dimNames(i))) then 
        dimUnlimIndex=i
        start(i)=timeLevel
        nread(i)=1
      endif
    enddo
    if (dimUnlimIndex .LE. 0) &
      call MOM_error(FATAL, "MOM_io::MOM_read_data_3d_noDD: time level specified, but variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
  endif
 
  ! read the data
  call read_data(fileObjRead, trim(fieldname), data, corner=start, edge_Lengths=nread)
  ! close the file
  if (check_if_open(fileObjRead)) call close_file(fileObjRead)

  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif

end subroutine MOM_read_data_3d_noDD

!> This function uses the fms_io function read_data to read 3-D data field named "fieldname" 
!! from file "filename".
subroutine MOM_read_data_4d_noDD(filename, fieldname, data, corner, edgeLengths, timeLevel, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:,:),   intent(inout) :: data !< The 4-dimensional array to pass to read_data
  integer, dimension(4),  optional, intent(in) :: corner !< starting indices of data buffer. Default is 1
  integer, dimension(4),  optional, intent(in) :: edgeLengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timeLevel !< time level to read
  real,         optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  ! local
  type(FmsNetcdfFile_t) :: fileObjRead !netCDF file object returned by call to MOM_open_file
  logical :: fileOpenSuccess !.true. if call to MOM_open_file is successful
  integer :: i, dimUnlimIndex
  integer, dimension(4) :: start, nread ! indices for first data value and number of values to read
  character(len=40), dimension(4) :: dimNames ! variable dimension names
  
  ! open the file
  if (.not.(check_if_open(fileObjRead))) &
    fileOpenSuccess = open_file(fileObjRead, filename, "read", is_restart=.false.)

  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  if (present(corner) .or. present(edgeLengths) .or. present(timeLevel)) then
    call get_variable_dimension_names(fileObjRead, trim(fieldname), dimNames)
  endif
  
  start(:) = 1
  if (present(corner)) start = corner
  
  if (present(edgeLengths)) then
    nread = edgeLengths
  else
    do i=1,4
      call get_dimension_size(fileObjRead, trim(dimNames(i)), nread(i))
    enddo
  endif

  if (present(timeLevel)) then
    dimUnlimIndex=0
    do i=1,4
      if (is_dimension_unlimited(fileObjRead,dimNames(i))) then 
        dimUnlimIndex=i
        start(i)=timeLevel
        nread(i)=1
      endif
    enddo
    if (dimUnlimIndex .LE. 0) &
      call MOM_error(FATAL, "MOM_io::MOM_read_data_4d_noDD: time level specified, but variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
  endif
 
  ! read the data
  call read_data(fileObjRead, trim(fieldname), data, corner=start, edge_Lengths=nread)
  ! close the file
  if (check_if_open(fileObjRead)) call close_file(fileObjRead)

  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif

end subroutine MOM_read_data_4d_noDD

!> register a MOM diagnostic axis to a domain-decomposed file 
subroutine MOM_register_diagnostic_axis(fileObj, axisName, axisLength)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObj !< netCDF file object returned by call to MOM_open_file
  character(len=*), intent(in) :: axisName !< name of the axis to register to file
  integer, intent(in), optional :: axisLength !< length of axis/dimension ;only needed for Layer, Interface, Time,
                                              !! Period
  select case (trim(axisName))      
    case ('latq'); call register_axis(fileObj,'latq','y', domain_position=NORTH_FACE)
    case ('lath'); call register_axis(fileObj,'lath','y', domain_position=CENTER) 
    case ('lonq'); call register_axis(fileObj,'lonq','x', domain_position=EAST_FACE) 
    case ('lonh'); call register_axis(fileObj,'lonh','x', domain_position=CENTER)
    case default
      if (.not. present(axisLength)) call MOM_error(FATAL,"MOM_io:register_diagnostic_axis: "//&
                        "An axis_length argument is required to register the axis "//trim(axisName))
      call register_axis(fileObj, axisName, axisLength) 
  end select
end subroutine MOM_register_diagnostic_axis

!> register axes associated with a variable from a domain-decomposed netCDF file
!> @note The user must specify units for variables with longitude/x-axis and/or latitude/y-axis axes to obtain 
!! the correct domain decomposition for the data buffer. 
subroutine MOM_register_variable_axes(fileObj, variableName, xUnits, yUnits, xPosition, yPosition)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObj !< netCDF file object returned by call to MOM_open_file
  character(len=*), intent(in) :: variableName !< name of the variable
  character(len=*), intent(in), optional :: xUnits !< x-axis (longitude) units to search for
  character(len=*), intent(in), optional :: yUnits !< y-axis (latitude) units to search for
  integer, intent(in), optional :: xPosition !< domain position of the x-axis
  integer, intent(in), optional :: yPosition !< domain position of the y-axis
  ! local
  character(len=40) :: units ! units corresponding to a specific variable dimension
  character(len=40), allocatable, dimension(:) :: dimNames ! variable dimension names
  integer :: i
  integer :: ndims ! number of dimensions
  integer :: xPos, yPos ! domain positions for x and y axes. Default is CENTER
  integer, allocatable, dimension(:) :: dimSizes ! variable dimension sizes
   
  if (.not. check_if_open(fileObj)) call MOM_error(FATAL,"MOM_io:register_variable_axes: The fileObj has "// &
                                                  "not been opened. Call MOM_open_file(fileObj,...) before "// &
                                                  "passing the fileObj argument to this function.")
  xPos=CENTER
  yPos=CENTER
  if (present(xPosition)) xPos=xPosition
  if (present(yPosition)) yPos=yPosition
 ! get variable dimension names and lengths
  ndims = get_variable_num_dimensions(fileObj, trim(variableName))
  allocate(dimSizes(ndims))
  allocate(dimNames(ndims))
  call get_variable_size(fileObj, trim(variableName), dimSizes, broadcast=.true.)
  call get_variable_dimension_names(fileObj, trim(variableName), dimNames)
  ! register the axes
  do i=1,ndims
    units=""
    if (present(xUnits)) then
      call get_variable_units(fileObj, dimNames(i), units)
      if (trim(lowercase(units)) .eq. trim(lowercase(xUnits))) then
        call register_axis(fileObj, dimNames(i),"x", domain_position=xPos)
      endif
    elseif (present(yUnits)) then
      call get_variable_units(fileObj, dimNames(i), units)
      if (trim(lowercase(units)) .eq. trim(lowercase(yUnits))) then
        call register_axis(fileObj, dimNames(i),"y", domain_position=yPos)
      endif
    else
      call register_axis(fileObj, dimNames(i), dimSizes(i))     
    endif
  enddo
  
  deallocate(dimSizes)
  deallocate(dimNames)
end subroutine MOM_register_variable_axes

!> set the "start" (corner) and "nread" (edge_lengths) arrays for domain-decomposed netCDF input data buffers 
subroutine MOM_get_nc_corner_edgelengths_DD(fileObj, variableName, corner, edgeLengths, myCorner, myCornerIndices, &
                                            myEdgeLengths, myEdgeLengthIndices)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObj !< netCDF file object returned by call to MOM_open_file
  character(len=*), intent(in), optional :: variableName !< name of the Varibl
  integer, allocatable, dimension(:), intent(out) :: corner !< array of corner indices to pass to read_data
  integer, allocatable, dimension(:), intent(out) :: edgeLengths !< array of edge_lengths indices to pass to read_data
  integer, dimension(:), intent(in), optional :: myCorner !< array of user-specified corner indices
  integer, dimension(:), intent(in), optional :: myCornerIndices !< positional indices for userCorner values
  integer, dimension(:), intent(in), optional :: myEdgeLengths !< array of user-specified edge_lengths indices
  integer, dimension(:), intent(in), optional :: myEdgeLengthIndices !< positional indices for myEdgelengths
  ! local
  integer :: i, idx, ndims ! counter, index, number of dimensions
   
  if (.not. check_if_open(fileObj)) call MOM_error(FATAL,"MOM_io:MOM_get_nc_corner_edgelengths_DD: "// &
                                                  "The fileObj has not been opened. Call MOM_open_file(fileObj,...)"// &
                                                  "before passing the fileObj argument to this function.")

 if (allocated(corner)) deallocate(corner)
 if (allocated(edgeLengths)) deallocate(edgeLengths)
 ! get variable dimension sizes and lengths
  ndims = get_variable_num_dimensions(fileObj, trim(variableName))
  allocate(corner(ndims))
  corner(:)=1
  allocate(edgeLengths(ndims))
  call get_variable_size(fileObj, trim(variableName), edgeLengths, broadcast=.true.)
  ! set user-specified corner values
  if (present(myCorner)) then
     if (.not.(present(myCornerIndices))) then 
       call MOM_error(FATAL,"MOM_io:MOM_get_nc_corner_edgelengths_DD: You passed the myCorner argument, but did "// &
                      "not pass myCornerIndices that defines the indices corresponding to the myCorner values.")
     endif
     do i=1,size(myCorner)
       idx = myCornerIndices(i)
       corner(idx)=myCorner(i)
     enddo
  endif
  
  ! set user-specified EdgeLengths values
  if (present(myEdgeLengths)) then
     if (.not.(present(myCornerIndices))) then
       call MOM_error(FATAL,"MOM_io:MOM_get_nc_corner_edgelengths_DD:You passed the myEdgeLengths argument, but "// &
                      "did not pass the myEdgeLengthIndices that defines the indices corresponding to the "// &
                      "myEdgeLengths values.")
     endif
     do i=1,size(myEdgeLengths)
       idx = myEdgeLengthIndices(i)
       edgeLengths(idx)=myEdgeLengths(i)
     enddo
  endif     

end subroutine MOM_get_nc_corner_edgelengths_DD

!> set the corner (start) and edgeLengths (count) arrays for non-domain-decomposed netCDF input data buffers
subroutine MOM_get_nc_corner_edgelengths_noDD(fileObj, variableName, corner, edgeLengths, myCorner, myCornerIndices, &
                                            myEdgeLengths, myEdgeLengthIndices)
  type(FmsNetcdfFile_t), intent(inout) :: fileObj !< netCDF file object returned by call to MOM_open_file
  character(len=*), intent(in), optional :: variableName !< name of the Varibl
  integer, allocatable, dimension(:), intent(out) :: corner !< array of corner indices to pass to read_data
  integer, allocatable, dimension(:), intent(out) :: edgeLengths !< array of edge_lengths indices to pass to read_data
  integer, dimension(:), intent(in), optional :: myCorner !< array of user-specified corner indices
  integer, dimension(:), intent(in), optional :: myCornerIndices !< positional indices for userCorner values
  integer, dimension(:), intent(in), optional :: myEdgeLengths !< array of user-specified edge_lengths indices
  integer, dimension(:), intent(in), optional :: myEdgeLengthIndices !< positional indices for myEdgelengths
  ! local
  integer :: i, idx, ndims ! counter, index, number of dimensions
   
  if (.not. check_if_open(fileObj)) call MOM_error(FATAL,"MOM_io:MOM_get_nc_corner_edgelengths_DD: "// &
                                                  "The fileObj has not been opened. Call "// &
                                                  "MOM_open_file(fileObj,...) before passing the fileObj argument "// &
                                                  "to this function.")

 if (allocated(corner)) deallocate(corner)
 if (allocated(edgeLengths)) deallocate(edgeLengths)
 ! get variable dimension sizes and lengths
  ndims = get_variable_num_dimensions(fileObj, trim(variableName))
  allocate(corner(ndims))
  corner(:)=1
  allocate(edgeLengths(ndims))
  call get_variable_size(fileObj, trim(variableName), edgeLengths, broadcast=.true.)
  ! set user-specified corner values
  if (present(myCorner)) then
     if (.not.(present(myCornerIndices))) then 
       call MOM_error(FATAL,"MOM_io:MOM_get_nc_corner_edgelengths_noDD: You passed the myCorner argument, but did "// &
                            "not pass myCornerIndices that defines the indices corresponding to the myCorner values")
     endif
     do i=1,size(myCorner)
       idx = myCornerIndices(i)
       corner(idx)=myCorner(i)
     enddo
  endif
  
  ! set user-specified EdgeLengths values
  if (present(myEdgeLengths)) then
     if (.not.(present(myCornerIndices))) then 
       call MOM_error(FATAL,"MOM_io:MOM_get_nc_corner_edgelengths_noDD: You passed the myEdgeLengths argument, "// &
                            "but did not pass myEdgeLengthIndices that defines the indices corresponding to the "// &
                            "myEdgeLengths values.")
     endif
     do i=1,size(myEdgeLengths)
       idx = myEdgeLengthIndices(i)
       edgeLengths(idx)=myEdgeLengths(i)
     enddo
  endif     
  
end subroutine MOM_get_nc_corner_edgelengths_noDD

!> Get the horizontal grid, vertical grid, and/or time dimension names and lengths
!! for a single variable from the grid ids returned by a prior call to query_vardesc
subroutine get_var_dimension_features(hor_grid, z_grid, t_grid_in, &
                                  dim_names, dim_lengths, num_dims, G, dG, GV)

  character(len=*), intent(in) :: hor_grid !< horizontal grid
  character(len=*), intent(in) :: z_grid !< vertical grid
  character(len=*), intent(in) :: t_grid_in !< time grid
  character(len=*), dimension(:), intent(out) :: dim_names !< array of dimension names
  integer, dimension(:), intent(out) :: dim_lengths !< array of dimension sizes
  integer, intent(out) ::  num_dims !< number of axes to register in the restart file
  type(ocean_grid_type),   optional,   intent(in)  :: G !< The ocean's grid structure
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                      !! is required if the new file uses any
                                                      !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure
  
  ! local
  logical :: use_lath
  logical :: use_lonh
  logical :: use_latq
  logical :: use_lonq
 
  character(len=8) :: t_grid
  character(len=8) :: t_grid_read
  integer :: isg, ieg, jsg, jeg, IsgB, IegB, JsgB, JegB
  real, pointer, dimension(:) :: gridLatT => NULL(), & ! The latitude or longitude of T or B points for
     gridLatB => NULL(), & ! the purpose of labeling the output axes.
     gridLonT => NULL(), &
     gridLonB => NULL()
  
  ! set the ocean grid coordinates
  if (present(G)) then
     gridLatT => G%gridLatT ; gridLatB => G%gridLatB
     gridLonT => G%gridLonT ; gridLonB => G%gridLonB
     isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg
     IsgB = G%IsgB ; IegB = G%IegB ; JsgB = G%JsgB ; JegB = G%JegB
  elseif (present(dG)) then
     gridLatT => dG%gridLatT ; gridLatB => dG%gridLatB
     gridLonT => dG%gridLonT ; gridLonB => dG%gridLonB
     isg = dG%isg ; ieg = dG%ieg ; jsg = dG%jsg ; jeg = dG%jeg
     IsgB = dG%IsgB ; IegB = dG%IegB ; JsgB = dG%JsgB ; JegB = dG%JegB
  endif

  use_lath = .false.
  use_lonh = .false.
  use_latq = .false.
  use_lonq = .false.
  
  call get_horizontal_grid_logic(hor_grid, use_lath, use_lonh, use_latq, use_lonq)
  ! add longitude name to dimension name array
  if (use_lonh) then
     num_dims = num_dims+1 
     dim_names(num_dims) = ''
     dim_names(num_dims)(1:len_trim('lonh')) = 'lonh'
     dim_lengths(num_dims) = size(gridLonT(isg:ieg))
  elseif (use_lonq) then
     num_dims = num_dims+1
     dim_names(num_dims) = ''
     dim_names(num_dims)(1:len_trim('lonq')) ='lonq'
     dim_lengths(num_dims) = size(gridLonB(IsgB:IegB)) 
  endif
  ! add latitude name to dimension name array
  if (use_lath) then
     num_dims = num_dims+1 
     dim_names(num_dims) = ''
     dim_names(num_dims)(1:len_trim('lath')) = 'lath'
     dim_lengths(num_dims) = size(gridLatT(jsg:jeg))
  elseif (use_latq) then
     num_dims = num_dims+1 
     dim_names(num_dims) = ''
     dim_names(num_dims)(1:len_trim('latq')) = 'latq'
     dim_lengths(num_dims) = size(gridLatB(JsgB:JegB))
  endif
  ! vertical grid
  select case (trim(z_grid))
     case ('L')
        num_dims = num_dims+1
        dim_names(num_dims) = ''
        dim_names(num_dims)(1:len_trim('Layer')) = 'Layer'
        dim_lengths(num_dims) = size(GV%sLayer(1:GV%ke))  
     case ('i')
        num_dims = num_dims+1
        dim_names(num_dims) = ''
        dim_names(num_dims)(1:len_trim('Interface')) = 'Interface'
        dim_lengths(num_dims) = size(GV%sInterface(1:GV%ke+1))
     case ('1') ! Do nothing.
     case default
        call MOM_error(FATAL, "MOM_io: get_var_dimension_features: "//&
                      " has an unrecognized z_grid argument"//trim(z_grid))
  end select
  ! time
  t_grid = adjustl(t_grid_in)
  select case (t_grid(1:1))
     case ('s', 'a', 'm')
        num_dims = num_dims+1
        dim_names(num_dims) = ''
        dim_names(num_dims)(1:len_trim('Time')) = 'Time'
        dim_lengths(num_dims) = unlimited
     case ('p')
        if (len_trim(t_grid(2:8)) <= 0) then
            call MOM_error(FATAL,"MOM_io:get_var_dimension_features: "//&
                           "No periodic axis length was specified in "//trim(t_grid))
        endif
        num_dims = num_dims+1
        dim_names(num_dims) = ''
        dim_names(num_dims)(1:len_trim('Period')) = 'Period'
        dim_lengths(num_dims) = unlimited
     case ('1') ! Do nothing.
     case default
           call MOM_error(WARNING, "MOM_io: get_var_dimension_features: "//&
                       "Unrecognized t_grid "//trim(t_grid))
  end select
end subroutine get_var_dimension_features

!> Populate the axis_data structure with axis data and attributes for diagnostic and restart files
subroutine MOM_get_diagnostic_axis_data(axis_data_CS, axis_name, axis_number, G, dG, GV, time_val, time_units)

  type(axis_data_type), intent(inout) :: axis_data_CS !< structure containing the axis data and metadata
  character(len=*), intent(in) :: axis_name !< name of the axis
  integer, intent(in) :: axis_number !< positional value (wrt to file) of the axis to register  
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                      !! is required if the file uses any
                                                      !! horizontal grid axes.
  type(verticalGrid_type), target, optional, intent(in) :: GV !< ocean vertical grid structure
  real,dimension(:), target, optional, intent(in) :: time_val !< time value
  character(len=*), optional,intent(in) :: time_units!< units for non-periodic time axis
  ! local
  character(len=40) :: x_axis_units='', y_axis_units=''
  integer :: isg, ieg, jsg, jeg, IsgB, IegB, JsgB, JegB
  real, pointer, dimension(:) :: gridLatT => NULL(), & ! The latitude or longitude of T or B points for
     gridLatB => NULL(), & ! the purpose of labeling the output axes.
     gridLonT => NULL(), &
     gridLonB => NULL()

  ! set the ocean grid coordinates
  if (present(G)) then
     gridLatT => G%gridLatT ; gridLatB => G%gridLatB
     gridLonT => G%gridLonT ; gridLonB => G%gridLonB
     x_axis_units = G%x_axis_units ; y_axis_units = G%y_axis_units
     isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg
     IsgB = G%IsgB ; IegB = G%IegB ; JsgB = G%JsgB ; JegB = G%JegB
  elseif (present(dG)) then
     gridLatT => dG%gridLatT ; gridLatB => dG%gridLatB
     gridLonT => dG%gridLonT ; gridLonB => dG%gridLonB
     x_axis_units = dG%x_axis_units ; y_axis_units = dG%y_axis_units
     isg = dG%isg ; ieg = dG%ieg ; jsg = dG%jsg ; jeg = dG%jeg
     IsgB = dG%IsgB ; IegB = dG%IegB ; JsgB = dG%JsgB ; JegB = dG%JegB
  endif
  
  ! initialize axis_data_CS elements
  axis_data_CS%axis(axis_number)%name = ''
  axis_data_CS%axis(axis_number)%longname = ''
  axis_data_CS%axis(axis_number)%units = ''
  axis_data_CS%axis(axis_number)%horgrid_position = 0
  axis_data_CS%axis(axis_number)%is_domain_decomposed = .false.
  axis_data_CS%axis(axis_number)%positive = ''
  axis_data_CS%data(axis_number)%p => NULL()
  
  select case(trim(axis_name))
     case('lath')
        axis_data_CS%data(axis_number)%p(jsg:jeg)=>gridLatT(jsg:jeg)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Latitude'
        axis_data_CS%axis(axis_number)%units = y_axis_units
        axis_data_CS%axis(axis_number)%horgrid_position  = CENTER
        axis_data_CS%axis(axis_number)%is_domain_decomposed  = .true.
     case('lonh')
        axis_data_CS%data(axis_number)%p(isg:ieg)=>gridLonT(isg:ieg)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%horgrid_position  = CENTER
        axis_data_CS%axis(axis_number)%longname = 'Longitude'
        axis_data_CS%axis(axis_number)%units = x_axis_units
        axis_data_CS%axis(axis_number)%is_domain_decomposed  = .true.
     case('latq')
        axis_data_CS%data(axis_number)%p(JsgB:JegB)=>gridLatB(JsgB:JegB)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Latitude'
        axis_data_CS%axis(axis_number)%units = y_axis_units
        axis_data_CS%axis(axis_number)%horgrid_position = NORTH_FACE
        axis_data_CS%axis(axis_number)%is_domain_decomposed = .true.
     case('lonq')
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%data(axis_number)%p(IsgB:IegB)=>gridLonB(IsgB:IegB)
        axis_data_CS%axis(axis_number)%longname  = 'Longitude'
        axis_data_CS%axis(axis_number)%units = x_axis_units
        axis_data_CS%axis(axis_number)%horgrid_position = EAST_FACE
        axis_data_CS%axis(axis_number)%is_domain_decomposed = .true.
     case('Layer')
        axis_data_CS%data(axis_number)%p=>GV%sLayer(1:GV%ke)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Layer pseudo-depth, -z*'
        axis_data_CS%axis(axis_number)%units = GV%zAxisUnits
        axis_data_CS%axis(axis_number)%positive  = 'up'
     case('Interface')
        axis_data_CS%data(axis_number)%p=>GV%sInterface(1:GV%ke+1)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Interface pseudo-depth, -z*'
        axis_data_CS%axis(axis_number)%units = GV%zAxisUnits
        axis_data_CS%axis(axis_number)%positive = 'up'
     case('Time')
        if (.not.(present(time_val))) then
           call MOM_error(FATAL, "MOM_io::get_diagnostic_axis_data: requires time_val"//&
                          " and time_units arguments for "//trim(axis_name))
        endif

        axis_data_CS%data(axis_number)%p=>time_val
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Time'

        if (present(time_units)) then
           axis_data_CS%axis(axis_number)%units = time_units
        else
           axis_data_CS%axis(axis_number)%units = 'days'
        endif
     case('Period')
        if (.not.(present(time_val))) then
           call MOM_error(FATAL, "MOM_io::get_diagnostic_axis_data: requires a time_val argument "// &
                          "for "//trim(axis_name))
        endif
        axis_data_CS%data(axis_number)%p=>time_val
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Periods for cyclical variables'
     case default
        call MOM_error(WARNING, "MOM_io::get_diagnostic_axis_data:"//trim(axis_name)//"is an unrecognized axis")
  end select

end subroutine MOM_get_diagnostic_axis_data

!> return the logic 
subroutine get_horizontal_grid_logic(grid_string_id, use_lath,use_lonh,use_latq,use_lonq)
  character(len=*), intent(in) :: grid_string_id !< horizontal grid string
  logical, intent(out) :: use_lath, use_lonh, use_latq, use_lonq !< logic indicating the x,y grid coordinates
                                                                 !! corresponding to hor_grid value
  use_lath = .false.
  use_lonh = .false.
  use_latq = .false.
  use_lonq = .false.
  select case (trim(grid_string_id))
     case ('h') ; use_lath = .true. ; use_lonh = .true. ! x=CENTER, y=CENTER
     case ('q') ; use_latq = .true. ; use_lonq = .true. ! x=EAST_FACE, y=NORTH_FACE
     case ('u') ; use_lath = .true. ; use_lonq = .true. ! x=EAST_FACE, y=CENTER
     case ('v') ; use_latq = .true. ; use_lonh = .true. ! x=CENTER, y=NORTH_FACE
     case ('T')  ; use_lath = .true. ; use_lonh = .true. ! x=CENTER, y=CENTER
     case ('Bu') ; use_latq = .true. ; use_lonq = .true. ! x=EAST_FACE, y=NORTH_FACE
     case ('Cu') ; use_lath = .true. ; use_lonq = .true. ! x=EAST_FACE, y=CENTER
     case ('Cv') ; use_latq = .true. ; use_lonh = .true. ! x=CENTER, y=NORTH_FACE
     case ('1') ; ! x=0, y=0
     case default
        call MOM_error(FATAL, "MOM_io:get_var_dimension_features "//&
                        "Unrecognized hor_grid argument "//trim(grid_string_id))
  end select
end subroutine get_horizontal_grid_logic
!> get the size of a variable in bytes
function get_variable_byte_size(hor_grid, z_grid, t_grid, G, num_zlevels) result(var_sz)
  character(len=*), intent(in) :: hor_grid !< horizontal grid string
  character(len=*), intent(in) :: z_grid !< vertical grid string
  character(len=*), intent(in) :: t_grid !< time string
  type(ocean_grid_type), intent(in) :: G !< The ocean's grid structure;
  integer, intent(in) :: num_zlevels     !< number of vertical levels
  ! local
  integer(kind=8) :: var_sz !< The size in bytes of each variable
  integer :: var_periods
  character(len=8) :: t_grid_read=''
  
  var_periods = 0
  
  if (trim(hor_grid) == '1') then
     var_sz = 8
  else
     var_sz = 8*(G%Domain%niglobal+1)*(G%Domain%njglobal+1)
  endif
  
  select case (trim(z_grid))
     case ('L') ; var_sz = var_sz * num_zlevels
     case ('i') ; var_sz = var_sz * (num_zlevels+1)
  end select

  if (adjustl(t_grid(1:1)) == 'p') then
     if (len_trim(t_grid(2:8)) > 0) then
        var_periods = -1
        t_grid_read = adjustl(t_grid(2:8))
        read(t_grid_read,*) var_periods
        if (var_periods > 1) var_sz = var_sz * var_periods
     endif
  endif

end function get_variable_byte_size

!> Define the time units for the input time value
function get_time_units(time_value) result(time_units_out)
   real, intent(in) :: time_value !< numerical time value in seconds
                                  !! i.e., before dividing by 86400.
   ! local
   character(len=10) :: time_units ! time units
   character(len=10) :: time_units_out ! time units trimmed
   time_units = ''
   time_units_out = ''
   if (time_value < 0.0) then
      time_units = "days" ! The default value.
   elseif (mod(time_value,86400.0)==0.0) then
      time_units = "days"
   elseif ((time_value >= 0.99) .and. (time_value < 1.01)) then
              time_units = "seconds"
   elseif ((time_value >= 3599.0) .and. (time_value < 3601.0)) then
              time_units = "hours"
   elseif ((time_value >= 86399.0) .and. (time_value < 86401.0)) then
              time_units = "days"
   elseif ((time_value >= 3.0e7) .and. (time_value < 3.2e7)) then
              time_units = "years"
   else
       write(time_units,'(es8.2," s")') time_value
   endif
   time_units_out = trim(time_units)
end function get_time_units

!> This function determines how many time levels a variable has.
function num_timelevels(filename, varname, min_dims) result(n_time)
  character(len=*),  intent(in) :: filename   !< name of the file to read
  character(len=*),  intent(in) :: varname    !< variable whose number of time levels
                                              !! are to be returned
  integer, optional, intent(in) :: min_dims   !< The minimum number of dimensions a variable must have
                                              !! if it has a time dimension.  If the variable has 1 less
                                              !! dimension than this, then 0 is returned.
  integer :: n_time                           !< number of time levels varname has in filename

  ! local
  logical :: fileOpenSuccess ! .true. if call to open_file is successful
  logical :: variableExists ! .true. if variable is found in file
  character(len=200) :: msg
  character(len=40), allocatable, dimension(:) :: dimNames ! variable dimension names
  integer :: i, ndims
  type(fmsNetcdfFile_t) :: fileObjread !netcdf file object returned by open_file

  n_time = -1

  ! open the file
  if (.not.(check_if_open(fileObjRead))) &
    fileOpenSuccess = open_file(fileObjRead, filename, "read", is_restart=.false.)

  ! check that variable is in the file
  if (.not.(variable_exists(fileObjRead, trim(varname)))) call MOM_error(FATAL, "num_timelevels: variable"//&
    trim(varname)//" not found in "//trim(filename))
  
  ! get the number of variable dimensions
  ndims = get_variable_num_dimensions(fileObjRead, trim(varname))

  if (present(min_dims)) then
    if (ndims .LT. min_dims-1) then
      write(msg, '(I3)') min_dims
      call MOM_error(WARNING, "num_timelevels: variable "//trim(varname)//&
        " in file "//trim(filename)//" has fewer than min_dims = "//trim(msg)//&
        " dimensions.")
    elseif (ndims .EQ. min_dims - 1) then
      n_time = 0 ; return
    endif
  endif

  ! check for the unlimited dimension and set n_time to the length of the unlimited dimension
  allocate(dimNames(ndims))

  call get_variable_dimension_names(fileObjRead, trim(varname), dimNames)
  do i=1,ndims
    if (is_dimension_unlimited(fileObjRead, trim(dimNames(i)))) &
      call get_dimension_size(fileObjRead, trim(dimNames(i)), n_time)
  enddo

  deallocate(dimNames)

end function num_timelevels

!> Returns a vardesc type whose elements have been filled with the provided
!! fields.  The argument name is required, while the others are optional and
!! have default values that are empty strings or are appropriate for a 3-d
!! tracer field at the tracer cell centers.
function var_desc(name, units, longname, hor_grid, z_grid, t_grid, &
                  cmor_field_name, cmor_units, cmor_longname, conversion, caller) result(vd)
  character(len=*),           intent(in) :: name               !< variable name
  character(len=*), optional, intent(in) :: units              !< variable units
  character(len=*), optional, intent(in) :: longname           !< variable long name
  character(len=*), optional, intent(in) :: hor_grid           !< variable horizonal staggering
  character(len=*), optional, intent(in) :: z_grid             !< variable vertical staggering
  character(len=*), optional, intent(in) :: t_grid             !< time description: s, p, or 1
  character(len=*), optional, intent(in) :: cmor_field_name    !< CMOR name
  character(len=*), optional, intent(in) :: cmor_units         !< CMOR physical dimensions of variable
  character(len=*), optional, intent(in) :: cmor_longname      !< CMOR long name
  real            , optional, intent(in) :: conversion         !< for unit conversions, such as needed to
                                                               !! convert from intensive to extensive
  character(len=*), optional, intent(in) :: caller             !< calling routine?
  type(vardesc)                          :: vd                 !< vardesc type that is created

  character(len=120) :: cllr
  cllr = "var_desc"
  if (present(caller)) cllr = trim(caller)

  call safe_string_copy(name, vd%name, "vd%name", cllr)

  vd%longname = "" ; vd%units = ""
  vd%hor_grid = 'h' ; vd%z_grid = 'L' ; vd%t_grid = 's'

  vd%cmor_field_name  =  ""
  vd%cmor_units       =  ""
  vd%cmor_longname    =  ""
  vd%conversion       =  1.0

  call modify_vardesc(vd, units=units, longname=longname, hor_grid=hor_grid, &
                      z_grid=z_grid, t_grid=t_grid,                          &
                      cmor_field_name=cmor_field_name,cmor_units=cmor_units, &
                      cmor_longname=cmor_longname, conversion=conversion, caller=cllr)

end function var_desc


!> This routine modifies the named elements of a vardesc type.
!! All arguments are optional, except the vardesc type to be modified.
subroutine modify_vardesc(vd, name, units, longname, hor_grid, z_grid, t_grid, &
                 cmor_field_name, cmor_units, cmor_longname, conversion, caller)
  type(vardesc),              intent(inout) :: vd              !< vardesc type that is modified
  character(len=*), optional, intent(in)    :: name            !< name of variable
  character(len=*), optional, intent(in)    :: units           !< units of variable
  character(len=*), optional, intent(in)    :: longname        !< long name of variable
  character(len=*), optional, intent(in)    :: hor_grid        !< horizonal staggering of variable
  character(len=*), optional, intent(in)    :: z_grid          !< vertical staggering of variable
  character(len=*), optional, intent(in)    :: t_grid          !< time description: s, p, or 1
  character(len=*), optional, intent(in)    :: cmor_field_name !< CMOR name
  character(len=*), optional, intent(in)    :: cmor_units      !< CMOR physical dimensions of variable
  character(len=*), optional, intent(in)    :: cmor_longname   !< CMOR long name
  real            , optional, intent(in)    :: conversion      !< for unit conversions, such as needed
                                                               !! to convert from intensive to extensive
  character(len=*), optional, intent(in)    :: caller          !< calling routine?

  character(len=120) :: cllr
  cllr = "mod_vardesc"
  if (present(caller)) cllr = trim(caller)

  if (present(name))      call safe_string_copy(name, vd%name, "vd%name", cllr)

  if (present(longname))  call safe_string_copy(longname, vd%longname, &
                               "vd%longname of "//trim(vd%name), cllr)
  if (present(units))     call safe_string_copy(units, vd%units,       &
                               "vd%units of "//trim(vd%name), cllr)
  if (present(hor_grid))  call safe_string_copy(hor_grid, vd%hor_grid, &
                               "vd%hor_grid of "//trim(vd%name), cllr)
  if (present(z_grid))    call safe_string_copy(z_grid, vd%z_grid,     &
                               "vd%z_grid of "//trim(vd%name), cllr)
  if (present(t_grid))    call safe_string_copy(t_grid, vd%t_grid,     &
                               "vd%t_grid of "//trim(vd%name), cllr)

  if (present(cmor_field_name)) call safe_string_copy(cmor_field_name, vd%cmor_field_name, &
                                     "vd%cmor_field_name of "//trim(vd%name), cllr)
  if (present(cmor_units))      call safe_string_copy(cmor_units, vd%cmor_units, &
                                     "vd%cmor_units of "//trim(vd%name), cllr)
  if (present(cmor_longname))   call safe_string_copy(cmor_longname, vd%cmor_longname, &
                                     "vd%cmor_longname of "//trim(vd%name), cllr)

end subroutine modify_vardesc

!> This function returns the CMOR standard name given a CMOR longname, based on
!! the standard pattern of character conversions.
function cmor_long_std(longname) result(std_name)
  character(len=*), intent(in) :: longname  !< The CMOR longname being converted
  character(len=len(longname)) :: std_name  !< The CMOR standard name generated from longname

  integer :: k

  std_name = lowercase(longname)

  do k=1, len_trim(std_name)
    if (std_name(k:k) == ' ') std_name(k:k) = '_'
  enddo

end function cmor_long_std

!> This routine queries vardesc
subroutine query_vardesc(vd, name, units, longname, hor_grid, z_grid, t_grid, &
                         cmor_field_name, cmor_units, cmor_longname, conversion, caller)
  type(vardesc),              intent(in)  :: vd                 !< vardesc type that is queried
  character(len=*), optional, intent(out) :: name               !< name of variable
  character(len=*), optional, intent(out) :: units              !< units of variable
  character(len=*), optional, intent(out) :: longname           !< long name of variable
  character(len=*), optional, intent(out) :: hor_grid           !< horiz staggering of variable
  character(len=*), optional, intent(out) :: z_grid             !< vert staggering of variable
  character(len=*), optional, intent(out) :: t_grid             !< time description: s, p, or 1
  character(len=*), optional, intent(out) :: cmor_field_name    !< CMOR name
  character(len=*), optional, intent(out) :: cmor_units         !< CMOR physical dimensions of variable
  character(len=*), optional, intent(out) :: cmor_longname      !< CMOR long name
  real            , optional, intent(out) :: conversion         !< for unit conversions, such as needed to
                                                                !! convert from intensive to extensive
  character(len=*), optional, intent(in)  :: caller             !< calling routine?


  character(len=120) :: cllr
  cllr = "mod_vardesc"
  if (present(caller)) cllr = trim(caller)

  if (present(name))      call safe_string_copy(vd%name, name,         &
                               "vd%name of "//trim(vd%name), cllr)
  if (present(longname))  call safe_string_copy(vd%longname, longname, &
                               "vd%longname of "//trim(vd%name), cllr)
  if (present(units))     call safe_string_copy(vd%units, units,       &
                               "vd%units of "//trim(vd%name), cllr)
  if (present(hor_grid))  call safe_string_copy(vd%hor_grid, hor_grid, &
                               "vd%hor_grid of "//trim(vd%name), cllr)
  if (present(z_grid))    call safe_string_copy(vd%z_grid, z_grid,     &
                               "vd%z_grid of "//trim(vd%name), cllr)
  if (present(t_grid))    call safe_string_copy(vd%t_grid, t_grid,     &
                               "vd%t_grid of "//trim(vd%name), cllr)

  if (present(cmor_field_name)) call safe_string_copy(vd%cmor_field_name, cmor_field_name, &
                                     "vd%cmor_field_name of "//trim(vd%name), cllr)
  if (present(cmor_units))      call safe_string_copy(vd%cmor_units, cmor_units,          &
                                     "vd%cmor_units of "//trim(vd%name), cllr)
  if (present(cmor_longname))   call safe_string_copy(vd%cmor_longname, cmor_longname, &
                                     "vd%cmor_longname of "//trim(vd%name), cllr)

end subroutine query_vardesc


!> Copies a string
subroutine safe_string_copy(str1, str2, fieldnm, caller)
  character(len=*),           intent(in)  :: str1    !< The string being copied
  character(len=*),           intent(out) :: str2    !< The string being copied into
  character(len=*), optional, intent(in)  :: fieldnm !< The name of the field for error messages
  character(len=*), optional, intent(in)  :: caller  !< The calling routine for error messages

  if (len(trim(str1)) > len(str2)) then
    if (present(fieldnm) .and. present(caller)) then
      call MOM_error(FATAL, trim(caller)//" attempted to copy the overly long"//&
        " string "//trim(str1)//" into "//trim(fieldnm))
    else
      call MOM_error(FATAL, "safe_string_copy: The string "//trim(str1)//&
                     " is longer than its intended target.")
    endif
  endif
  str2 = trim(str1)
end subroutine safe_string_copy


!> Returns a name with "%#E" or "%E" replaced with the ensemble member number.
function ensembler(name, ens_no_in) result(en_nm)
  character(len=*),  intent(in) :: name       !< The name to be modified
  integer, optional, intent(in) :: ens_no_in  !< The number of the current ensemble member
  character(len=len(name)) :: en_nm  !< The name encoded with the ensemble number

  ! This function replaces "%#E" or "%E" with the ensemble number anywhere it
  ! occurs in name, with %E using 4 or 6 digits (depending on the ensemble size)
  ! and %#E using # digits, where # is a number from 1 to 9.

  character(len=len(name)) :: tmp
  character(10) :: ens_num_char
  character(3)  :: code_str
  integer :: ens_no
  integer :: n, is, ie

  en_nm = trim(name)
  if (index(name,"%") == 0) return

  if (present(ens_no_in)) then
    ens_no = ens_no_in
  else
    ens_no = get_ensemble_id()
  endif

  write(ens_num_char, '(I10)') ens_no ; ens_num_char = adjustl(ens_num_char)
  do
    is = index(en_nm,"%E")
    if (is == 0) exit
    if (len(en_nm) < len(trim(en_nm)) + len(trim(ens_num_char)) - 2) &
      call MOM_error(FATAL, "MOM_io ensembler: name "//trim(name)// &
      " is not long enough for %E expansion for ens_no "//trim(ens_num_char))
    tmp = en_nm(1:is-1)//trim(ens_num_char)//trim(en_nm(is+2:))
    en_nm = tmp
  enddo

  if (index(name,"%") == 0) return

  write(ens_num_char, '(I10.10)') ens_no
  do n=1,9 ; do
    write(code_str, '("%",I1,"E")') n

    is = index(en_nm,code_str)
    if (is == 0) exit
    if (ens_no < 10**n) then
      if (len(en_nm) < len(trim(en_nm)) + n-3) call MOM_error(FATAL, &
        "MOM_io ensembler: name "//trim(name)//" is not long enough for %E expansion.")
      tmp = en_nm(1:is-1)//trim(ens_num_char(11-n:10))//trim(en_nm(is+3:))
    else
      call MOM_error(FATAL, "MOM_io ensembler: Ensemble number is too large "//&
          "to be encoded with "//code_str//" in "//trim(name))
    endif
    en_nm = tmp
  enddo ; enddo

end function ensembler

!> This routine uses the fms_io function read_data to read a pair of distributed
!! 2-D data fields with names given by "[uv]_fieldname" from file "filename".  Valid values for
!! "stagger" include CGRID_NE, BGRID_NE, and AGRID.
subroutine MOM_read_vector_2d(filename, u_fieldname, v_fieldname, u_data, v_data, MOM_Domain, &
                              timelevel, stagger, scale)
  character(len=*),       intent(in)    :: filename !< name of the netcdf file to read
  character(len=*),       intent(in)    :: u_fieldname !< The variable name of the u data in the file
  character(len=*),       intent(in)    :: v_fieldname !< The variable name of the v data in the file
  real, dimension(:,:),   intent(inout) :: u_data    !< The 2 dimensional array into which the
                                                     !! u-component of the data should be read
  real, dimension(:,:),   intent(inout) :: v_data    !< The 2 dimensional array into which the
                                                     !! v-component of the data should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: stagger   !< A flag indicating where this vector is discretized
  real,         optional, intent(in)    :: scale     !< A scaling factor that the fields are multiplied
                                                     !! by before they are returned.
  ! local
  type(FmsNetcdfDomainFile_t) :: fileObjRead !< netcdf file object returned by call to MOM_open_file
  integer :: is, ie, js, je, i
  integer :: u_pos, v_pos
  integer, dimension(2) :: start, dim_sizes_u, dim_sizes_v
  character(len=32), dimension(2) :: dim_names_u, dim_names_v, units_u, units_v
  logical :: fileOpenSuccess ! .true. if open file is successful
  
 
  !open the file
  if (.not.(check_if_open(fileObjRead))) &
    fileOpenSuccess = open_file(fileObjRead, filename, "read", MOM_domain%mpp_domain, is_restart=.false.)
  if (.not. fileOpenSuccess) call MOM_error(FATAL, "MOM_read_vector_2d: netcdf file "//trim(filename)//" not opened.")

  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE .or. stagger == BGRID_NE ) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif

  !call old_fms_read_data(filename, u_fieldname, u_data, MOM_Domain%mpp_domain, &
   !              timelevel=timelevel, position=u_pos)
  !call old_fms_read_data(filename, v_fieldname, v_data, MOM_Domain%mpp_domain, &
  !               timelevel=timelevel, position=v_pos)
  start(:) = 1
  call get_variable_size(fileObjRead, u_fieldname, dim_sizes_u, broadcast=.true.)
  call get_variable_size(fileObjRead, v_fieldname, dim_sizes_v, broadcast=.true.)
  call get_variable_dimension_names(fileObjRead, u_fieldname, dim_names_u, broadcast=.true.)
  call get_variable_dimension_names(fileObjRead, v_fieldname, dim_names_v, broadcast=.true.)
 
  do i=1,2
    ! register the u axes
    call get_variable_units(fileObjRead, dim_names_u(i), units_u(i))
    if (trim(lowercase(units_u(i))) .eq. "degrees_east") then
      call register_axis(fileObjRead, dim_names_u(i), "x", domain_position=u_pos)
    elseif (trim(lowercase(units_u(i))) .eq. "degrees_north") then 
      call register_axis(fileObjRead, dim_names_u(i), "y", domain_position=u_pos)  
    else
      if (is_dimension_unlimited(fileObjRead, dim_names_u(i))) then
        if (present(timelevel)) then
          start(i)=timelevel
          dim_sizes_u(i) = 1
          dim_sizes_v(i) = 1
        endif
        call register_axis(fileObjRead, dim_names_u(i),dim_sizes_u(i))
      endif
    endif
    ! Register the v axes if they differ from the u axes
    if (trim(lowercase(dim_names_v(i))) .ne. trim(lowercase(dim_names_u(i)))) then 
      call get_variable_units(fileObjRead, dim_names_v(i), units_v(i))
      if (trim(lowercase(units_v(i))) .eq. "degrees_east") then
        call register_axis(fileObjRead, dim_names_v(i), "x", domain_position=v_pos)
      elseif (trim(lowercase(units_v(i))) .eq. "degrees_north") then
        call register_axis(fileObjRead, dim_names_v(i), "y", domain_position=v_pos)  
      endif 
    endif
  enddo 

  call read_data(fileObjRead,u_fieldname, u_data, corner=start, edge_lengths=dim_sizes_u)
  call read_data(fileObjRead,v_fieldname, v_data, corner=start, edge_lengths=dim_sizes_v)

  ! close the file 
  if (check_if_open(fileObjRead)) call close_file(fileObjRead)

  if (present(scale)) then ; if (scale /= 1.0) then
    call get_simple_array_i_ind(MOM_Domain, size(u_data,1), is, ie)
    call get_simple_array_j_ind(MOM_Domain, size(u_data,2), js, je)
    u_data(is:ie,js:je) = scale*u_data(is:ie,js:je)
    call get_simple_array_i_ind(MOM_Domain, size(v_data,1), is, ie)
    call get_simple_array_j_ind(MOM_Domain, size(v_data,2), js, je)
    v_data(is:ie,js:je) = scale*v_data(is:ie,js:je)
  endif ; endif

end subroutine MOM_read_vector_2d

!> This routine uses the fms_io function read_data to read a pair of distributed
!! 3-D data fields with names given by "[uv]_fieldname" from file "filename".  Valid values for
!! "stagger" include CGRID_NE, BGRID_NE, and AGRID.
subroutine MOM_read_vector_3d(filename, u_fieldname, v_fieldname, u_data, v_data, MOM_Domain, &
                              timelevel, stagger, scale)
  character(len=*),       intent(in)    :: filename !< name of the netcdf file to read
  character(len=*),       intent(in)    :: u_fieldname !< The variable name of the u data in the file
  character(len=*),       intent(in)    :: v_fieldname !< The variable name of the v data in the file
  real, dimension(:,:,:), intent(inout) :: u_data    !< The 3 dimensional array into which the
                                                     !! u-component of the data should be read
  real, dimension(:,:,:), intent(inout) :: v_data    !< The 3 dimensional array into which the
                                                     !! v-component of the data should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: stagger   !< A flag indicating where this vector is discretized
  real,         optional, intent(in)    :: scale     !< A scaling factor that the fields are multiplied
                                                     !! by before they are returned.
  ! local
  type(FmsNetcdfDomainFile_t) :: fileObjRead !< netcdf file object returned by call to MOM_open_file
  integer :: is, ie, js, je, i
  integer :: u_pos, v_pos
  integer, dimension(3) :: start, dim_sizes_u, dim_sizes_v
  character(len=32), dimension(3) :: dim_names_u, dim_names_v, units_u, units_v
  logical :: fileOpenSuccess ! .true. if open file is successful
  ! open the file
  if (.not.(check_if_open(fileObjRead))) &
    fileOpenSuccess = open_file(fileObjRead, filename, "read", MOM_domain%mpp_domain, is_restart=.false.)
  if (.not. fileOpenSuccess) call MOM_error(FATAL, "MOM_read_vector_3d: netcdf file "//trim(filename)//" not opened.")

  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == BGRID_NE) then ; u_pos = CORNER ; v_pos = CORNER
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif

  !call old_fms_read_data(filename, u_fieldname, u_data, MOM_Domain%mpp_domain, &
  !               timelevel=timelevel, position=u_pos)
  !call old_fms_read_data(filename, v_fieldname, v_data, MOM_Domain%mpp_domain, &
  !               timelevel=timelevel, position=v_pos)
  start(:) = 1
  call get_variable_size(fileObjRead, u_fieldname, dim_sizes_u, broadcast=.true.)
  call get_variable_size(fileObjRead, v_fieldname, dim_sizes_v, broadcast=.true.)
  call get_variable_dimension_names(fileObjRead, u_fieldname, dim_names_u, broadcast=.true.)
  call get_variable_dimension_names(fileObjRead, v_fieldname, dim_names_v, broadcast=.true.)
 
  do i=1,3
    ! register the u axes
    call get_variable_units(fileObjRead, dim_names_u(i), units_u(i))
    if (trim(lowercase(units_u(i))) .eq. "degrees_east") then
      call register_axis(fileObjRead, dim_names_u(i), "x", domain_position=u_pos)
    elseif (trim(lowercase(units_u(i))) .eq. "degrees_north") then 
      call register_axis(fileObjRead, dim_names_u(i), "y", domain_position=u_pos)  
    else
      if (is_dimension_unlimited(fileObjRead, dim_names_u(i))) then
        if (present(timelevel)) then
          start(i)=timelevel
          dim_sizes_u(i) = 1
          dim_sizes_v(i) = 1
        endif
        call register_axis(fileObjRead, dim_names_u(i), dim_sizes_u(i))
      endif
    endif
  ! register the v axes if the differ from the u axes
    if (trim(lowercase(dim_names_v(i))) .ne. trim(lowercase(dim_names_u(i)))) then 
      call get_variable_units(fileObjRead, dim_names_v(i), units_v(i))
      if (trim(lowercase(units_v(i))) .eq. "degrees_east") then
        call register_axis(fileObjRead, dim_names_v(i), "x", domain_position=v_pos)
      elseif (trim(lowercase(units_v(i))) .eq. "degrees_north") then
        call register_axis(fileObjRead, dim_names_v(i), "y", domain_position=v_pos)  
      endif 
    endif
  enddo 

  call read_data(fileObjRead,u_fieldname, u_data, corner=start, edge_lengths=dim_sizes_u)
  call read_data(fileObjRead,v_fieldname, v_data, corner=start, edge_lengths=dim_sizes_v)

  ! close the file 
  if (check_if_open(fileObjRead)) call close_file(fileObjRead)

  if (present(scale)) then ; if (scale /= 1.0) then
    call get_simple_array_i_ind(MOM_Domain, size(u_data,1), is, ie)
    call get_simple_array_j_ind(MOM_Domain, size(u_data,2), js, je)
    u_data(is:ie,js:je,:) = scale*u_data(is:ie,js:je,:)
    call get_simple_array_i_ind(MOM_Domain, size(v_data,1), is, ie)
    call get_simple_array_j_ind(MOM_Domain, size(v_data,2), js, je)
    v_data(is:ie,js:je,:) = scale*v_data(is:ie,js:je,:)
  endif ; endif

end subroutine MOM_read_vector_3d

subroutine scale_data_1d(data, scale_factor)
  real, dimension(:), intent(inout) :: data !< The 1-dimensional data array
  real, intent(in) :: scale_factor !< Scale factor

  if (scale_factor /= 1.0) then
    data(:) = scale_factor*data(:)
  endif
end subroutine scale_data_1d

subroutine scale_data_2d(data, scale_factor, MOM_domain)
  real, dimension(:,:), intent(inout) :: data !< The 2-dimensional data array
  real, intent(in) :: scale_factor !< Scale factor
  type(MOM_domain_type), optional, intent(in) :: MOM_Domain !< The domain that describes the decomposition
  ! local
  integer :: is, ie, js, je

  if (scale_factor /= 1.0) then
    if (present(MOM_domain)) then
      call get_simple_array_i_ind(MOM_Domain, size(data,1), is, ie)
      call get_simple_array_j_ind(MOM_Domain, size(data,2), js, je)
      data(is:ie,js:je) = scale_factor*data(is:ie,js:je)
    else
      data(:,:) = scale_factor*data(:,:)
    endif
  endif
end subroutine scale_data_2d

subroutine scale_data_3d(data, scale_factor, MOM_domain)
  real, dimension(:,:,:), intent(inout) :: data !< The 3-dimensional data array
  real, intent(in) :: scale_factor !< Scale factor
  type(MOM_domain_type), optional, intent(in) :: MOM_Domain !< The domain that describes the decomposition
  ! local
  integer :: is, ie, js, je

  if (scale_factor /= 1.0) then
    if (present(MOM_domain)) then
      call get_simple_array_i_ind(MOM_Domain, size(data,1), is, ie)
      call get_simple_array_j_ind(MOM_Domain, size(data,2), js, je)
      data(is:ie,js:je,:) = scale_factor*data(is:ie,js:je,:)
    else
      data(:,:,:) = scale_factor*data(:,:,:)
    endif
  endif
end subroutine scale_data_3d

subroutine scale_data_4d(data, scale_factor, MOM_domain)
  real, dimension(:,:,:,:), intent(inout) :: data !< The 4-dimensional data array
  real, intent(in) :: scale_factor !< Scale factor
  type(MOM_domain_type), optional, intent(in) :: MOM_Domain !< The domain that describes the decomposition
  ! local
  integer :: is, ie, js, je

  if (scale_factor /= 1.0) then
    if (present(MOM_domain)) then
      call get_simple_array_i_ind(MOM_Domain, size(data,1), is, ie)
      call get_simple_array_j_ind(MOM_Domain, size(data,2), js, je)
      data(is:ie,js:je,:,:) = scale_factor*data(is:ie,js:je,:,:)
    else
      data(:,:,:,:) = scale_factor*data(:,:,:,:)
    endif
  endif
end subroutine scale_data_4d

!> Initialize the MOM_io module
subroutine MOM_io_init(param_file)
  type(param_file_type), intent(in) :: param_file  !< structure indicating the open file to
                                                   !! parse for model parameter values.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_io" ! This module's name.

  call log_version(param_file, mdl, version)

end subroutine MOM_io_init


!> \namespace mom_io
!!
!!   This file contains a number of subroutines that manipulate
!!  NetCDF files and handle input and output of fields.  These
!!  subroutines, along with their purpose, are:
!!
!!   * MOM_open_file: open a netcdf file in read, write, overwrite, or append mode
!!   * close_file: close an open netcdf file.
!!   * write_data: write a field to an open file.
!!   * MOM_read_data: read a field from a netcdf file and apply a scaling factor if specified.
!!   * MOM_read_vector : read in the components (u,v) of a vector field and apply a scaling factor to the data
!!    if specified
!!   * scale_data: apply a scaling factor to a data field

end module MOM_io
