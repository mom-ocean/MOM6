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
use MOM_string_functions, only : extract_word
use mpp_mod,              only : mpp_max
use mpp_domains_mod,      only : domain1d, domain2d, domainug, mpp_get_domain_components
use mpp_domains_mod,      only : mpp_get_domain_npes, mpp_define_io_domain, mpp_get_io_domain
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
use fms2_io_mod,          only : check_if_open, &
                                get_dimension_names, &
                                get_dimension_size, &
                                get_compute_domain_dimension_indices, &
                                get_global_attribute, &
                                get_global_io_domain_indices, &
                                get_num_dimensions, &
                                get_num_variables, &
                                get_unlimited_dimension_name, &
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
public :: num_timelevels, MOM_read_vector, ensembler, create_file
public :: slasher, MOM_io_init, field_exists, field_size, read_axis_data
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
public :: get_unlimited_dimension_name
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
public :: MOM_get_nc_corner_edge_lengths
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
public :: write_field
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
interface MOM_get_nc_corner_edge_lengths
  module procedure MOM_get_nc_corner_edge_lengths_DD
  module procedure MOM_get_nc_corner_edge_lengths_noDD
end interface

!> interface to read data from a netcdf file
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
  module procedure MOM_read_data_2d_noDD_diag_axes
end interface

!> Read a pair of data fields representing the two components of a vector from a netcdf file
interface MOM_read_vector
  module procedure MOM_read_vector_3d
  module procedure MOM_read_vector_2d
end interface

!> interface to write data to a netcdf file generated by create_file
interface write_field
  module procedure write_field_4d_DD
  module procedure write_field_3d_DD
  module procedure write_field_2d_DD
  module procedure write_field_1d_DD
  module procedure write_scalar
  module procedure write_field_4d_noDD
  module procedure write_field_3d_noDD
  module procedure write_field_2d_noDD
  module procedure write_field_1d_noDD
end interface

!> interface to scale data after reading in a field
interface scale_data
  module procedure scale_data_4d
  module procedure scale_data_3d
  module procedure scale_data_2d
  module procedure scale_data_1d
end interface

!> interface to read the most recent time from a netCDF file
interface read_most_recent_time
  module procedure read_most_recent_time_DD
  module procedure read_most_recent_time_noDD
end interface

contains

!> This routine opens a netcdf file in "write" or "overwrite" mode, registers the global diagnostic axes, and writes
!! the axis data and metadata to the file
subroutine create_file(filename, vars, numVariables, threading, timeUnit, register_time, G, DG, GV, checksums)
  character(len=*),      intent(in)               :: filename !< full path to the netcdf file
  type(vardesc), dimension(:), intent(in)         :: vars !< structures describing the output
  integer,               intent(in)               :: numVariables !< number of variables to write to the file
  integer, optional,     intent(in)               :: threading !< SINGLE_FILE or MULTIPLE
  real, optional,        intent(in)               :: timeUnit !< length of the units for time [s]. The
                                                             !! default value is 86400.0, for 1 day.
  logical, optional, intent(in) :: register_time !< if .true., register a time dimension to the file
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums(:,:)  !< checksums of the variables

  ! local
  type(FmsNetcdfFile_t) :: fileObjNoDD ! non-domain-decomposed netcdf file object returned by open_file
  type(FmsNetcdfDomainFile_t) :: fileObjDD ! domain-decomposed netcdf file object returned by open_file
  type(axis_data_type) :: axis_data_CS ! structure for coordinate variable metadata
  type(MOM_domain_type), pointer :: Domain => NULL()
  logical :: file_open_successDD, file_open_successNoDD ! true if netcdf file is opened
  logical :: one_file, domain_set ! indicates whether the file will be domain-decomposed or not
  logical :: reg_time ! register the time if .true.
  character(len=10) :: timeUnits, nc_mode
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=1024) :: filename_temp
  character(len=48), allocatable, dimension(:,:) :: dim_names ! variable dimension names
  integer :: i, j, total_axes, substring_index
  integer :: num_dims !< number of dimensions
  integer :: thread ! indicates whether threading is used
  integer, dimension(4) :: dim_lengths !< variable dimension lengths
  real :: time

  ! determine whether the file will be domain-decomposed or not
  thread = SINGLE_FILE
  if (PRESENT(threading)) thread = threading

  domain_set=.false.
  if (present(G)) then
    domain_set = .true. ; Domain => G%Domain
  elseif (present(dG)) then
    domain_set = .true. ; Domain => dG%Domain
  endif

  one_file = .true.
  if (domain_set) one_file = (thread == SINGLE_FILE)

  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index < 1) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif

  nc_mode = ""
  if (file_exists(trim(filename_temp))) then
    nc_mode = "overwrite"
  else
    nc_mode = "write"
  endif

  reg_time = .false.
  if (present(register_time)) reg_time = .true.

  ! set the time units
  timeUnits=""
  if (present(timeUnit)) then
    timeUnits = get_time_units(timeUnit)
  else
    timeUnits ="days"
  endif

  ! open the file
  file_open_successNoDD=.false.
  file_open_successDD=.false.

  if (domain_set) then
    ! define the io domain if on one pe and the io domain is not set
    if (mpp_get_domain_npes(domain%mpp_domain) .eq. 1 ) then
      if (.not. associated(mpp_get_io_domain(domain%mpp_domain))) &
        call mpp_define_io_domain(domain%mpp_domain, (/1,1/))
    endif

    file_open_successDD=open_file(fileObjDD, filename_temp, trim(nc_mode), Domain%mpp_domain, is_restart=.false.)
  else
    file_open_successNoDD=open_file(fileObjNoDD, filename_temp, trim(nc_mode), is_restart=.false.)
  endif

  ! allocate the output data variable dimension attributes
  allocate(dim_names(numVariables,4))

  ! allocate the axis data and attribute types for the file
  !> \note The user should increase the sizes of the axis and data attributes to accommodate more axes if necessary.
  allocate(axis_data_CS%axis(7))
  allocate(axis_data_CS%data(7))

  ! axis registration procedure for the domain-decomposed case

  if (file_open_successDD) then
    do i=1,numVariables
      num_dims=0
      dim_lengths(:) = 0
      dim_names(:,:) = ""

      if (present(G)) then
        call get_var_dimension_features(vars(i)%hor_grid, vars(i)%z_grid, vars(i)%t_grid, dim_names(i,:), &
                                          dim_lengths, num_dims, G=G)
      elseif(present(dG)) then
        call get_var_dimension_features(vars(i)%hor_grid, vars(i)%z_grid, vars(i)%t_grid, dim_names(i,:), &
                                          dim_lengths, num_dims, dG=dG)
      endif
      if(present(GV)) &
        call get_var_dimension_features(vars(i)%hor_grid, vars(i)%z_grid, vars(i)%t_grid, dim_names(i,:), &
                                        dim_lengths, num_dims, GV=GV)

      if (num_dims .le. 0) call MOM_error(FATAL, "MOM_io:create_file: num_dims is an invalid value.")
      ! register the global axes to the file
      do j=1,num_dims
        if (dim_lengths(j) .gt. 0) then
          if (.not.(dimension_exists(fileObjDD, dim_names(i,j)))) then
            if (present(G)) then
              if (present(timeUnit)) then
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, G=G, &
                                                  time_val=(/timeUnit/), time_units=timeUnits)
              else
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, G=G, &
                                                  time_val=(/1.0/), time_units=timeUnits)
              endif

            elseif (present(dG)) then
              if (present(timeUnit)) then
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, dG=dG, &
                                                    time_val=(/timeUnit/), time_units=timeUnits)
              else
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, dG=dG, &
                                                    time_val=(/1.0/), time_units=timeUnits)
              endif
            endif

            if (present(GV)) then
              if (present(timeUnit)) then
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, GV=GV, &
                                                    time_val=(/timeUnit/), time_units=timeUnits)
              else
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, GV=GV, &
                                                    time_val=(/1.0/), time_units=timeUnits)
              endif
            endif

            call MOM_register_diagnostic_axis(fileObjDD, trim(dim_names(i,j)), dim_lengths(j))
          endif
          ! register the axis attributes and write the axis data to the file
          if (.not.(variable_exists(fileObjDD, trim(axis_data_CS%axis(j)%name)))) then
            if (fileObjDD%is_root) then
              if (associated(axis_data_CS%data(j)%p)) then
                call register_field(fileObjDD, trim(axis_data_CS%axis(j)%name), &
                                  "double", dimensions=(/trim(axis_data_CS%axis(j)%name)/))

                call register_variable_attribute(fileObjDD, trim(axis_data_CS%axis(j)%name), &
                                               'long_name', axis_data_CS%axis(j)%longname)

                call register_variable_attribute(fileObjDD, trim(axis_data_CS%axis(j)%name), &
                                               'units', trim(axis_data_CS%axis(j)%units))

                !> \note: create_file does not write the initial time value
                if (lowercase(trim(axis_data_CS%axis(i)%name)) .ne. 'time') then
                  call write_data(fileObjDD, trim(axis_data_CS%axis(j)%name), axis_data_CS%data(j)%p)
                endif
              endif
            endif
          endif
        endif
      enddo
    enddo

    if (reg_time) then
     if (.not.(dimension_exists(fileObjDD,"Time"))) &
       call register_axis(fileObjDD, "Time", unlimited)
    endif

    if (check_if_open(fileObjDD)) call close_file(fileObjDD)
  ! axis registration procedure for the non-domain-decomposed case
  elseif (file_open_successNoDD) then
    do i=1,numVariables
      num_dims=0
      dim_lengths(:) = 0
      dim_names(:,:) = ""

      if (present(G)) then
        call get_var_dimension_features(vars(i)%hor_grid, vars(i)%z_grid, vars(i)%t_grid, dim_names(i,:), &
                                        dim_lengths, num_dims, G=G)
      elseif(present(dG)) then
        call get_var_dimension_features(vars(i)%hor_grid, vars(i)%z_grid, vars(i)%t_grid, dim_names(i,:), &
                                        dim_lengths, num_dims, dG=dG)
      endif
      if (present(GV)) &
        call get_var_dimension_features(vars(i)%hor_grid, vars(i)%z_grid, vars(i)%t_grid, dim_names(i,:), &
                                        dim_lengths, num_dims, GV=GV)

      if (num_dims .le. 0) call MOM_error(FATAL, "MOM_io:create_file: num_dims is an invalid value.")
      ! register the global axes to the file
      do j=1,num_dims
        if (dim_lengths(j) .gt. 0) then
          if (.not.(dimension_exists(fileObjNoDD, dim_names(i,j)))) then
            !total_axes=total_axes+1
            if (present(G)) then
              if (present(timeUnit)) then
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, G=G, &
                                                  time_val=(/timeUnit/), time_units=timeUnits)
              else
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, G=G, &
                                                  time_val=(/1.0/), time_units=timeUnits)
              endif

            elseif (present(dG)) then
              if (present(timeUnit)) then
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, dG=dG, &
                                                  time_val=(/timeUnit/), time_units=timeUnits)
              else
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, dG=dG, &
                                                  time_val=(/1.0/), time_units=timeUnits)
              endif
            endif

            if (present(GV)) then
              if (present(timeUnit)) then
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, GV=GV, &
                                                  time_val=(/timeUnit/), time_units=timeUnits)
              else
                call MOM_get_diagnostic_axis_data(axis_data_CS, dim_names(i,j), j, GV=GV, &
                                                  time_val=(/1.0/), time_units=timeUnits)
              endif
            endif

            call register_axis(fileObjNoDD, trim(dim_names(i,j)), dim_lengths(j))
          endif
          ! register the axis attributes and write the axis data to the file
          if (.not.(variable_exists(fileObjNoDD, trim(axis_data_CS%axis(j)%name)))) then
            if (fileObjNoDD%is_root) then
              if (associated(axis_data_CS%data(j)%p)) then
                call register_field(fileObjNoDD, trim(axis_data_CS%axis(j)%name), &
                                  "double", dimensions=(/trim(axis_data_CS%axis(j)%name)/))

                call register_variable_attribute(fileObjNoDD, trim(axis_data_CS%axis(j)%name), &
                                               'long_name', axis_data_CS%axis(j)%longname)

                call register_variable_attribute(fileObjNoDD, trim(axis_data_CS%axis(j)%name), &
                                               'units', trim(axis_data_CS%axis(j)%units))

                !> \note: create_file does not write the initial time value
                if (lowercase(trim(axis_data_CS%axis(j)%name)) .ne. 'time') then
                  call write_data(fileObjNoDD, trim(axis_data_CS%axis(j)%name), axis_data_CS%data(j)%p)
                endif
              endif
            endif
          endif
        endif
      enddo
    enddo

    if (reg_time) then
     if (.not.(dimension_exists(fileObjNoDD,"Time"))) &
       call register_axis(fileObjNoDD, "Time" , unlimited)
    endif

    if (check_if_open(fileObjNoDD)) call close_file(fileObjNoDD)
  endif

  deallocate(dim_names)
  deallocate(axis_data_CS%axis)
  deallocate(axis_data_CS%data)
  nullify(Domain)

end subroutine create_file

!> This function uses the fms_io function write_data to write a 1-D domain-decomposed data field named "fieldname"
!! to the file "filename" in "write", "overwrite", or "append" mode. It should be called after create_file in the MOM
!! file write procedure.
subroutine write_field_1d_DD(filename, fieldname, data, mode, domain, var_desc, &
                             corner, edge_lengths, time_level, scale, checksums, G, dG, GV)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:), intent(in) :: data !< The 1-dimensional data array to pass to read_data
  character(len=*), intent(in) :: mode !< "write", "overwrite", or "append"
  type(MOM_domain_type),intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  type(vardesc), optional, intent(in) :: var_desc !< structure describing variable output
  integer, dimension(1), optional, intent(in) :: corner !< starting index of data buffer. Default is 1
  integer, dimension(1), optional, intent(in) :: edge_lengths !< number of data values to read in; default is the
                                                          !! variable size
  real, optional, intent(in) :: time_level !< time value to write
  real, optional, intent(in) :: scale !< A scaling factor that the fields are multiplied by before they are written.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums  !< variable checksum
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  ! local
  type(FmsNetcdfDomainFile_t) :: fileobj ! netCDF domain-decomposed file object returned by call to open_file

  logical :: file_open_success !.true. if call to open_file is successful
  real :: file_time ! most recent time currently written to file
  real, allocatable, dimension(:) :: data_tmp
  integer :: num_dims, time_index, substring_index
  integer :: dim_unlim_size, dim_unlim_index ! size of the unlimited dimension
  integer :: is, ie
  integer, dimension(1) :: start, nwrite ! indices for first data value and number of values to write
  character(len=40) :: dim_unlim_name ! name of the unlimited dimension in the file
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=1024) :: filename_temp
  character(len=48), dimension(2) :: dim_names !< variable dimension names (or name, in the 1-D case)
  integer, dimension(2) :: dim_lengths !< variable dimension lengths (or length, in the 1-D case)

  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index <= 0) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif

  ! define the io domain for 1-pe jobs because it is required to write domain-decomposed files
  if (mpp_get_domain_npes(domain%mpp_domain) .eq. 1 ) then
    if (.not. associated(mpp_get_io_domain(domain%mpp_domain))) &
      call mpp_define_io_domain(domain%mpp_domain, (/1,1/))
  endif
  ! open the file for a domain-decomposed write
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, filename_temp, lowercase(trim(mode)), domain%mpp_domain, &
                                      is_restart=.false.)
  ! register the field if it is not in the file
  dim_unlim_index=0
  dim_names(:) = ""
  dim_lengths(:) = 0
  num_dims=0
  if (variable_exists(fileobj, trim(fieldname))) then
      dim_unlim_index = get_variable_unlimited_dimension_index(fileobj, trim(fieldname))
  else
    ! get the dimension names and lengths
    if (present(G)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, G=G)
    elseif(present(dG)) then  
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, dG=dG)
    endif

    if present(GV)) &
    call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                     dim_lengths, num_dims, GV=GV)

    call register_field(fileObj, trim(fieldname), "double", dimensions=dim_names(1:num_dims))
    call register_variable_attribute(fileObj, trim(fieldname), 'units', trim(var_desc%units))
    call register_variable_attribute(fileObj, trim(fieldname), 'long_name', trim(var_desc%longname))
    ! write the checksum attribute
    if (present(checksums)) then
      ! convert the checksum to a string
      checksum_char = ''
      checksum_char = convert_checksum_to_string(checksums(1,1))
      call register_variable_attribute(fileObj, trim(fieldname), "checksum", checksum_char)
    endif
  endif

  start(:) = 1
  nwrite(:) = size(data)
  if (present(corner)) then
    if (corner(1) .le. nwrite(1)) then
      if (corner(1) .gt. 0) then
          start(1) = corner(1)
      endif
    endif
  endif

  if (present(edge_lengths)) then
    if (edge_lengths(1) .le. nwrite(1)) then
      if (edge_lengths(1) .gt. 0) then
        nwrite(1) = edge_lengths(1)
      endif
    endif
  endif

  allocate(data_tmp(size(data,1)))
  data_tmp = data
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data_tmp,scale)
  endif ; endif

  time_index=1
  dim_unlim_size=0
  dim_unlim_name=""
  ! write the data
  call get_global_io_domain_indices(fileObjWrite, trim(axis_data_CS%axis(i)%name), is, ie)
  if (lowercase(trim(mode)) .eq. "append") then
    if (present(time_level)) then
      ! write the time value if it is not already written to the file
      file_time = read_most_recent_time(fileobj)
      if (time_level .gt. file_time+EPSILON(time_level)) then
        call get_unlimited_dimension_name(fileobj,dim_unlim_name)
        call get_dimension_size(fileobj, trim(dim_unlim_name), dim_unlim_size)
        if (dim_unlim_size .lt. unlimited) time_index=dim_unlim_size+1
        call write_data(fileobj, trim(dim_unlim_name), (/time_level/), corner=(/time_index/),edge_lengths=(/1/))
      endif
      if (dim_unlim_index .gt. 0) &
        call write_data(fileobj, trim(fieldname), data_tmp(is:ie), corner=start, edge_lengths=nwrite, &
                        unlim_dim_level=time_index)
    endif
  else
    call write_data(fileobj, trim(fieldname), data_tmp(is:ie), corner=start, edge_lengths=nwrite)
  endif
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)
  deallocate(data_tmp)

end subroutine write_field_1d_DD

!> This function uses the fms_io function write_data to write a 2-D domain-decomposed data field named "fieldname"
!! to the file "filename" in "write", "overwrite", or "append" mode. It should be called after create_file in the MOM
!! file write procedure.
subroutine write_field_2d_DD(filename, fieldname, data, mode, domain, var_desc, &
                             corner, edge_lengths, time_level, scale, checksums, G, dG, GV)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:), intent(in) :: data !< The 2-dimensional data array to pass to read_data
  character(len=*), intent(in) :: mode !< "write", "overwrite", or "append"
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  type(vardesc), optional, intent(in) :: var_desc !< structure describing variable output
  integer, dimension(2), optional, intent(in) :: corner !< starting index of data buffer. Default is 1
  integer, dimension(2), optional, intent(in) :: edge_lengths !< number of data values to read in; default is the
                                                          !! variable size
  real, optional, intent(in) :: time_level !< time value to write
  real, optional, intent(in) :: scale !< A scaling factor that the fields are multiplied by before they are written.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums  !< variable checksum
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  ! local
  type(FmsNetcdfDomainFile_t) :: fileobj ! netCDF file object returned by call to open_file
  type(FmsNetcdfFile_t) :: fileobj_nodd ! netCDF file object returned by call to open_file
  logical :: file_dd_success !.true. if call to open_file is successful
  logical :: var_exists ! .true. if variable exists in file
  real :: file_time ! most recent time currently written to file
  real, allocatable, dimension(:,:) :: data_tmp
  integer :: i, j, is, ie, js, je, ndims, num_dims, time_index, substring_index
  integer :: dim_unlim_size, dim_unlim_index ! size and dimension index of the unlimited dimension
  integer, dimension(2) :: start, nwrite ! indices for starting points and number of values to write
  character(len=40) :: dim_unlim_name ! name of the unlimited dimension in the file
  character(len=1024) :: filename_temp
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=48), dimension(3) :: dim_names !< variable dimension names
  integer, dimension(3) :: dim_lengths !< variable dimension lengths

  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index <= 0) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif
  ! if the job is running on more than 1 pe, and the io_layout is defined, use the domain-decomposed IO interfaces
  dim_names(:) = ""
  dim_lengths(:) = 0
  num_dims = 0
  dim_unlim_index = 0
  file_dd_success = .false.
  var_exists = .false.
  
  ! define the io domain for 1-pe jobs because it is required to write domain-decomposed files
  if (mpp_get_domain_npes(domain%mpp_domain) .eq. 1 ) then
    if (.not. associated(mpp_get_io_domain(domain%mpp_domain))) call mpp_define_io_domain(domain%mpp_domain, (/1,1/))
  endif

  if (.not.(check_if_open(fileobj))) &
    file_dd_success = open_file(fileobj, trim(filename_temp), lowercase(trim(mode)), domain%mpp_domain, &
                                is_restart=.false.)
  if (file_dd_success) then
    if (variable_exists(fileobj, trim(fieldname))) then
      dim_unlim_index = get_variable_unlimited_dimension_index(fileobj, trim(fieldname))
      var_exists = .true.
    endif
  else
    call MOM_error(FATAL, "MOM_io:write_field_2d_dd: unable to open file "//trim(filename_temp))
  endif
  if (.not.(var_exists)) then
    ! get the dimension names and lengths
    if (present(G)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, G=G)
    elseif(present(dG)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, dG=dG)
    endif
    if (present(GV)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, GV=GV)
    endif
    ! register the variable and its attributes
    call register_field(fileObj, trim(fieldname), "double", dimensions=dim_names(1:num_dims))
    call register_variable_attribute(fileObj, trim(fieldname), 'units', trim(var_desc%units))
    call register_variable_attribute(fileObj, trim(fieldname), 'long_name', trim(var_desc%longname))
    ! write the checksum attribute
    if (present(checksums)) then
      ! convert the checksum to a string
      checksum_char = ''
      checksum_char = convert_checksum_to_string(checksums(1,1))
      call register_variable_attribute(fileObj_nodd, trim(fieldname), "checksum", checksum_char)
    endif
  endif

  ndims=2
  start(:) = 1
  nwrite(:) = (/size(data,1), size(data,2)/)
  if (present(corner)) then
    do i=1,ndims
      if (corner(i) .le. nwrite(i)) then
        if (corner(i) .gt. 0) then
          if (i .ne. dim_unlim_index) start(i) = corner(i)
        endif
      endif
    enddo
  endif

  if (present(edge_lengths)) then
    do i=1,ndims
      if (edge_lengths(i) .le. nwrite(i)) then
        if (edge_lengths(i) .gt. 0) then
          if (i .ne. dim_unlim_index) nwrite(i) = edge_lengths(i)
        endif
      endif
    enddo
  endif

  allocate(data_tmp(size(data,1), size(data,2)))
  data_tmp = data

  !> /note remove this DEBUGGING LOOP
  do i=1,size(data_tmp,1)
    do j=1,size(data_tmp,2)
      write(*, "(A,A,F5.2)") trim(fieldname), " value is ", data_tmp(I,J)
    enddo
  enddo
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data_tmp,scale)
  endif ; endif

  time_index=1
  dim_unlim_size=0
  dim_unlim_name=""
  ! write the data
  call get_global_io_domain_indices(fileObjWrite, trim(axis_data_CS%axis(i)%name), is, ie, js, je)
  if (lowercase(trim(mode)) .eq. "append") then
    if (present(time_level)) then
      ! write the time value if it is not already written to the file
      file_time = read_most_recent_time(fileobj)
      if (time_level .gt. file_time+EPSILON(time_level)) then
        call get_unlimited_dimension_name(fileobj,dim_unlim_name)
        call get_dimension_size(fileobj, trim(dim_unlim_name), dim_unlim_size)
        if (dim_unlim_size .lt. unlimited) time_index=dim_unlim_size+1

        call write_data(fileobj, trim(dim_unlim_name), (/time_level/), corner=(/time_index/), &
                        edge_lengths=(/1/))
      endif
      if (dim_unlim_index .gt. 0) &
        call write_data(fileobj, trim(fieldname), data_tmp(is:ie,js:je), corner=start, edge_lengths=nwrite, &
                        unlim_dim_level=time_index)
    endif
  else
    call write_data(fileobj, trim(fieldname), data_tmp(is:ie,js:je), corner=start, edge_lengths=nwrite)
  endif
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)

  deallocate(data_tmp)
end subroutine write_field_2d_DD

!> This function uses the fms_io function write_data to write a 3-D domain-decomposed data field named "fieldname"
!! to the file "filename" in "write", "overwrite", or "append" mode. It should be called after create_file in the MOM
!! file write procedure.
subroutine write_field_3d_DD(filename, fieldname, data, mode, domain, var_desc, &
                             corner, edge_lengths, time_level, scale, checksums, G, dG, GV)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:), intent(in) :: data !< The 3-dimensional data array to pass to read_data
  character(len=*), intent(in) :: mode !< "write", "overwrite", or "append"
  type(MOM_domain_type),  intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  type(vardesc), optional, intent(in) :: var_desc !< structure describing variable output
  integer, dimension(3), optional, intent(in) :: corner !< starting index of data buffer. Default is 1
  integer, dimension(3), optional, intent(in) :: edge_lengths !< number of data values to read in; default is the
                                                          !! variable size
  real, optional, intent(in) :: time_level !< time value to write
  real, optional, intent(in) :: scale !< A scaling factor that the fields are multiplied by before they are written.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums  !< variable checksum
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  ! local
  type(FmsNetcdfDomainFile_t) :: fileobj ! netCDF file object returned by call to open_file
  logical :: file_open_success !.true. if call to open_file is successful
  real :: file_time ! most recent time currently written to file
  real, allocatable, dimension(:,:,:) :: data_tmp
  integer :: i, ndims, num_dims, time_index, substring_index
  integer :: dim_unlim_size, dim_unlim_index ! size and dimension index of the unlimited dimension
  integer, dimension(3) :: start, nwrite ! indices for first data value and number of values to write
  character(len=40) :: dim_unlim_name ! name of the unlimited dimension in the file
  character(len=1024) :: filename_temp
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=48), dimension(4) :: dim_names !< variable dimension names
  integer, dimension(4) :: dim_lengths !< variable dimension lengths

  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index <= 0) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif
  ! define the io domain for 1-pe jobs because it is required to write domain-decomposed files
  if (mpp_get_domain_npes(domain%mpp_domain) .eq. 1 ) then
    if (.not. associated(mpp_get_io_domain(domain%mpp_domain))) &
      call mpp_define_io_domain(domain%mpp_domain, (/1,1/))
  endif
  ! open the file for a domain-decomposed write
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, trim(filename_temp), lowercase(trim(mode)), domain%mpp_domain, &
                                  is_restart=.false.)
  dim_names(:) = ""
  dim_lengths(:) = 0
  num_dims = 0
  dim_unlim_index = 0
  if (variable_exists(fileobj, trim(fieldname))) then
      dim_unlim_index = get_variable_unlimited_dimension_index(fileobj, trim(fieldname))
  else
    ! get the dimension names and lengths
    if (present(G)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, G=G)
    elseif(present(dG)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, dG=dG)
    endif

    if (present(GV)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, GV=GV)

    call register_field(fileObj, trim(fieldname), "double", dimensions=dim_names(1:num_dims))
    call register_variable_attribute(fileObj, trim(fieldname), 'units', trim(var_desc%units))
    call register_variable_attribute(fileObj, trim(fieldname), 'long_name', trim(var_desc%longname))
    ! write the checksum attribute
    if (present(checksums)) then
      ! convert the checksum to a string
      checksum_char = ''
      checksum_char = convert_checksum_to_string(checksums(1,1))
      call register_variable_attribute(fileObj, trim(fieldname), "checksum", checksum_char)
    endif
  endif

  ndims = 3
  start(:) = 1
  nwrite(:) = (/size(data,1), size(data,2), size(data,3)/)
  if (present(corner)) then
    do i=1,ndims
      if (corner(i) .le. nwrite(i)) then
        if (corner(i) .gt. 0) then
          if (i .ne. dim_unlim_index) start(i) = corner(i)
        endif
      endif
    enddo
  endif

  if (present(edge_lengths)) then
    do i=1,ndims
      if (edge_lengths(i) .le. nwrite(i)) then
        if (edge_lengths(i) .gt. 0) then
          if (i .ne. dim_unlim_index) nwrite(i) = edge_lengths(i)
        endif
      endif
    enddo
  endif

  allocate(data_tmp(size(data,1), size(data,2), size(data,3)))
  data_tmp = data
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data_tmp,scale)
  endif ; endif

  time_index=1
  dim_unlim_size=0
  dim_unlim_name=""
  ! write the data
  if (lowercase(trim(mode)) .eq. "append") then
    if (present(time_level)) then
      ! write the time value if it is not already written to the file
      file_time = read_most_recent_time(fileobj)
      if (time_level .gt. file_time+EPSILON(time_level)) then
        call get_unlimited_dimension_name(fileobj,dim_unlim_name)
        call get_dimension_size(fileobj, trim(dim_unlim_name), dim_unlim_size)
        if (dim_unlim_size .lt. unlimited) time_index=dim_unlim_size+1
        call write_data(fileobj, trim(dim_unlim_name), (/time_level/), corner=(/time_index/),edge_lengths=(/1/))
      endif

    if (dim_unlim_index .gt. 0) &
      call write_data(fileobj, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite, &
                      unlim_dim_level=time_index)
    endif
  else
    call write_data(fileobj, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite)
  endif
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)
  deallocate(data_tmp)

end subroutine write_field_3d_DD

!> This function uses the fms_io function write_data to write a 4-D domain-decomposed data field named "fieldname"
!! to the file "filename" in "write", "overwrite", or "append" mode. It should be called after create_file in the MOM
!! file write procedure.
subroutine write_field_4d_DD(filename, fieldname, data, mode, domain, var_desc, &
                             corner, edge_lengths, time_level, scale, checksums, G, dG, GV)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:,:), intent(in) :: data !< The 4-dimensional data array to pass to read_data
  character(len=*), intent(in) :: mode !< "write", "overwrite", or "append"
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  type(vardesc), optional, intent(in) :: var_desc !< structure describing variable output
  integer, dimension(4), optional, intent(in) :: corner !< starting index of data buffer. Default is 1
  integer, dimension(4), optional, intent(in) :: edge_lengths !< number of data values to read in; default is the
                                                          !! variable size
  real, optional, intent(in) :: time_level !< time value to write
  real, optional, intent(in) :: scale !< A scaling factor that the fields are multiplied by before they are written.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums  !< variable checksum
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  ! local
  type(FmsNetcdfDomainFile_t) :: fileobj ! netCDF file object returned by call to open_file
  logical :: file_open_success !.true. if call to open_file is successful
  real, allocatable, dimension(:,:,:,:) :: data_tmp
  real :: file_time ! most recent time currently written to file
  integer :: i, ndims, num_dims, time_index, substring_index
  integer :: dim_unlim_size, dim_unlim_index ! size and dimension index of the unlimited dimension
  integer, dimension(4) :: start, nwrite ! indices for first data value and number of values to write
  character(len=40) :: dim_unlim_name ! name of the unlimited dimension in the file
  character(len=1024) :: filename_temp
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=48), dimension(4) :: dim_names !< variable dimension names
  integer, dimension(4) :: dim_lengths !< variable dimension lengths

  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index <= 0) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif
  ! define the io domain for 1-pe jobs because it is required to write domain-decomposed files
  if (mpp_get_domain_npes(domain%mpp_domain) .eq. 1 ) then
    if (.not. associated(mpp_get_io_domain(domain%mpp_domain))) &
      call mpp_define_io_domain(domain%mpp_domain, (/1,1/))
  endif
  ! open the file for a domain-decomposed write
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, trim(filename_temp), lowercase(trim(mode)), domain%mpp_domain, &
                                  is_restart=.false.)
  dim_names(:) = ""
  dim_lengths(:) = 0
  num_dims = 0
  dim_unlim_index = 0
  if (variable_exists(fileobj, trim(fieldname))) then
      dim_unlim_index = get_variable_unlimited_dimension_index(fileobj, trim(fieldname))
  else
    ! get the dimension names and lengths
    if (present(G)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, G=G)
    elseif(present(dG)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                       dim_lengths, num_dims, dG=dG)
    endif

    if (present(GV)) &
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, GV=GV)

    ! register the field to the file
    call register_field(fileObj, trim(fieldname), "double", dimensions=dim_names(1:num_dims))

    call register_variable_attribute(fileObj, trim(fieldname), 'units', trim(var_desc%units))
    call register_variable_attribute(fileObj, trim(fieldname), 'long_name', trim(var_desc%longname))
    ! write the checksum attribute
    if (present(checksums)) then
      ! convert the checksum to a string
      checksum_char = ''
      checksum_char = convert_checksum_to_string(checksums(1,1))
      call register_variable_attribute(fileObj, trim(fieldname), "checksum", checksum_char)
    endif
  endif

  ndims = 4
  start(:) = 1
  nwrite(:) = (/size(data,1), size(data,2), size(data,3), size(data,4)/)
  if (present(corner)) then
    do i=1,ndims
      if (corner(i) .le. nwrite(i)) then
        if (corner(i) .gt. 0) then
          if (i .ne. dim_unlim_index) start(i) = corner(i)
        endif
      endif
    enddo
  endif

  if (present(edge_lengths)) then
    do i=1,ndims
      if (edge_lengths(i) .le. nwrite(i)) then
        if (edge_lengths(i) .gt. 0) then
          if (i .ne. dim_unlim_index) nwrite(i) = edge_lengths(i)
        endif
      endif
    enddo
  endif

  allocate(data_tmp(size(data,1), size(data,2), size(data,3), size(data,4)))
  data_tmp = data
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data_tmp,scale)
  endif ; endif

  time_index=1
  dim_unlim_size=0
  dim_unlim_name=""
  ! write the data
  if (lowercase(trim(mode)) .eq. "append") then
    if (present(time_level)) then
      ! write the time value if it is not already written to the file
      file_time = read_most_recent_time(fileobj)
      if (time_level .gt. file_time+EPSILON(time_level)) then
        call get_unlimited_dimension_name(fileobj,dim_unlim_name)
        call get_dimension_size(fileobj, trim(dim_unlim_name), dim_unlim_size)
        if (dim_unlim_size .lt. unlimited) time_index=dim_unlim_size+1
        call write_data(fileobj, trim(dim_unlim_name), (/time_level/), corner=(/time_index/),edge_lengths=(/1/))
      endif
      if (dim_unlim_index .gt. 0) &
        call write_data(fileobj, trim(fieldname), data_tmp, corner=start, edge_Lengths=nwrite, &
                        unlim_dim_level=time_index)
    endif
  else
    call write_data(fileobj, trim(fieldname), data_tmp, corner=start, edge_Lengths=nwrite)
  endif
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)
  deallocate(data_tmp)

end subroutine write_field_4d_DD

!> This routine uses the fms_io function write_data to write a scalar variable named "fieldname"
!! to the file "filename" in "write", "overwrite", or "append" mode. It should be called after create_file in the MOM
!! file write procedure.
subroutine write_scalar(filename, fieldname, data, mode, time_level, var_desc, G, dG, GV)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, intent(in) :: data !< The 1-dimensional data array to pass to read_data
  character(len=*), intent(in) :: mode !< "write", "overwrite", or "append"
  real, optional, intent(in) :: time_level !< time value to write
  type(vardesc), optional, intent(in) :: var_desc !< structure describing variable output
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  ! local
  type(FmsNetcdfFile_t) :: fileobj ! netCDF file object returned by call to open_file
  logical :: file_open_success !.true. if call to open_file is successful
  integer :: num_dims, time_index, substring_index
  integer :: dim_unlim_size, dim_unlim_index ! size and dimension index of the unlimited dimension
  real :: file_time ! most recent time currently written to file
  character(len=40) :: dim_unlim_name ! name of the unlimited dimension in the file
  character(len=1024) :: filename_temp
  character(len=48), dimension(2) :: dim_names !< variable dimension names
  integer, dimension(2) :: dim_lengths !< variable dimension lengths

  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index <= 0) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif
  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, trim(filename_temp), lowercase(trim(mode)), is_restart=.false.)

  dim_unlim_index=0
  num_dims=0
  dim_names=""
  dim_lengths(:)=0
  if (variable_exists(fileobj, trim(fieldname))) then
    dim_unlim_index = get_variable_unlimited_dimension_index(fileobj, trim(fieldname))
  else
    ! get the dimension names and lengths
    if (present(G)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, G=G)
    elseif(present(dG)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, dG=dG)
    endif
    if (present(GV)) &
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, GV=GV)
    endif
    call register_field(fileObj, trim(fieldname), "double", dimensions=dim_names(1:num_dims))
    call register_variable_attribute(fileObj, trim(fieldname), 'units', trim(var_desc%units))
    call register_variable_attribute(fileObj, trim(fieldname), 'long_name', trim(var_desc%longname))
  endif
  time_index=1
  dim_unlim_size=0
  dim_unlim_name=""
  ! write the data
  if (lowercase(trim(mode)) .eq. "append") then
    if (present(time_level)) then
      ! write the time value if it is not already written to the file
      file_time = read_most_recent_time(fileobj)
      if (time_level .gt. file_time+EPSILON(time_level)) then
        call get_unlimited_dimension_name(fileobj,dim_unlim_name)
        call get_dimension_size(fileobj, trim(dim_unlim_name), dim_unlim_size)
        if (dim_unlim_size .lt. unlimited) time_index=dim_unlim_size+1

        call write_data(fileobj, trim(dim_unlim_name), (/time_level/), corner=(/time_index/))
      endif
      if (dim_unlim_index .gt. 0) &
          call write_data(fileobj, trim(fieldname), data, unlim_dim_level=time_index)
    endif
  else
    call write_data(fileobj, trim(fieldname), data)
  endif
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)
end subroutine write_scalar

!> This function uses the fms_io function write_data to write a 1-D non-domain-decomposed data field named "fieldname"
!! to the file "filename" in "write", "overwrite", or "append" mode. It should be called after create_file in the MOM
!! file write procedure.
subroutine write_field_1d_noDD(filename, fieldname, data, mode, var_desc, &
                               corner, edge_lengths, time_level, scale, checksums, G, dG, GV)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:), intent(in) :: data !< The 1-dimensional data array to pass to read_data
  character(len=*), intent(in) :: mode !< "write", "overwrite", or "append"
  type(vardesc), optional, intent(in) :: var_desc !< structure describing variable output
  integer, dimension(1), optional, intent(in) :: corner !< starting index of data buffer. Default is 1
  integer, dimension(1), optional, intent(in) :: edge_lengths !< number of data values to read in; default is the
                                                          !! variable size
  real, optional, intent(in) :: time_level !< time value to write
  real, optional, intent(in) :: scale !< A scaling factor that the fields are multiplied by before they are written.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums  !< variable checksum
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  ! local
  type(FmsNetcdfFile_t) :: fileobj ! netCDF file object returned by call to open_file
  logical :: file_open_success !.true. if call to open_file is successful
  real :: file_time ! most recent time currently written to file
  real, allocatable, dimension(:) :: data_tmp
  integer :: i, ndims, num_dims, time_index, substring_index
  integer :: dim_unlim_size, dim_unlim_index ! size of the unlimited dimension
  integer, dimension(1) :: start, nwrite ! indices for first data value and number of values to write
  character(len=40) :: dim_unlim_name ! name of the unlimited dimension in the file
  character(len=1024) :: filename_temp
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=48), dimension(2) :: dim_names !< variable dimension names (up to 2 if appended at time level)
  integer, dimension(2) :: dim_lengths !< variable dimension lengths

  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index <= 0) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif
  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, trim(filename_temp), lowercase(trim(mode)), is_restart=.false.)

  dim_names(:) = ""
  dim_lengths(:) = 0
  num_dims=0
  dim_unlim_index=0
  if (variable_exists(fileobj, trim(fieldname))) then
    dim_unlim_index = get_variable_unlimited_dimension_index(fileobj, trim(fieldname))
  else
    ! get the dimension names and lengths
    if (present(G)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, G=G)
    elseif(present(dG)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, dG=dG)
    endif
  
    if (present(GV)) &
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, GV=GV)

    call register_field(fileObj, trim(fieldname), "double", dimensions=dim_names(1:num_dims))
    call register_variable_attribute(fileObj, trim(fieldname), 'units', trim(var_desc%units))
    call register_variable_attribute(fileObj, trim(fieldname), 'long_name', trim(var_desc%longname))
    ! write the checksum attribute
    if (present(checksums)) then
      ! convert the checksum to a string
      checksum_char = ''
      checksum_char = convert_checksum_to_string(checksums(1,1))
      call register_variable_attribute(fileObj, trim(fieldname), "checksum", checksum_char)
    endif
  endif

  start(:) = 1
  nwrite(:) = size(data)
  if (present(corner)) then
    if (corner(1) .le. nwrite(1)) then
      if (corner(1) .gt. 0) then
        start(1) = corner(1)
      endif
    endif
  endif

  if (present(edge_lengths)) then
    if (edge_lengths(1) .le. nwrite(1)) then
      if (edge_lengths(1) .gt. 0) then
        nwrite(1) = edge_lengths(1)
      endif
    endif
  endif

  allocate(data_tmp(size(data,1)))
  data_tmp = data
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data_tmp,scale)
  endif ; endif

  time_index=1
  dim_unlim_size=0
  dim_unlim_name=""
  ! write the data
  if (lowercase(trim(mode)) .eq. "append") then
    if (present(time_level)) then
      ! write the time value if it is not already written to the file
      file_time = read_most_recent_time(fileobj)
      if (time_level .gt. file_time+EPSILON(time_level)) then
        call get_unlimited_dimension_name(fileobj,dim_unlim_name)
        call get_dimension_size(fileobj, trim(dim_unlim_name), dim_unlim_size)
        if (dim_unlim_size .lt. unlimited) time_index=dim_unlim_size+1

        call write_data(fileobj, trim(dim_unlim_name), (/time_level/), corner=(/time_index/),edge_lengths=(/1/))
      endif
      if (dim_unlim_index .gt. 0) &
          call write_data(fileobj, trim(fieldname), data_tmp, corner=start, edge_Lengths=nwrite, &
                          unlim_dim_level=time_index)
    endif
  else
    call write_data(fileobj, trim(fieldname), data_tmp, corner=start, edge_Lengths=nwrite)
  endif
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)
  deallocate(data_tmp)
end subroutine write_field_1d_noDD

!> This function uses the fms_io function write_data to write a scalar variable named "fieldname"
!! to the file "filename" in "write", "overwrite", or "append" mode. It should be called after create_file in the MOM
!! file write procedure.
subroutine write_field_2d_noDD(filename, fieldname, data, mode, var_desc, &
                               corner, edge_lengths, time_level, scale, checksums, G, dG, GV)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:), intent(in) :: data !< The 2-dimensional data array to pass to read_data
  character(len=*), intent(in) :: mode !< "write", "overwrite", or "append"
  type(vardesc), optional, intent(in) :: var_desc !< structure describing variable output
  integer, dimension(2), optional, intent(in) :: corner !< starting index of data buffer. Default is 1
  integer, dimension(2), optional, intent(in) :: edge_lengths !< number of data values to read in; default is the
                                                          !! variable size
  real, optional, intent(in) :: time_level !< time value to write
  real, optional, intent(in) :: scale !< A scaling factor that the fields are multiplied by before they are written.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums  !< variable checksum
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  ! local
  type(FmsNetcdfFile_t) :: fileobj ! netCDF file object returned by call to open_file
  logical :: file_open_success !.true. if call to open_file is successful
  real :: file_time ! most recent time currently written to file
  real, allocatable, dimension (:,:) :: data_tmp
  integer :: i, ndims, num_dims, time_index, substring_index
  integer :: dim_unlim_size, dim_unlim_index ! size and dimension index of the unlimited dimension
  integer, dimension(2) :: start, nwrite ! indices for starting points and number of values to write
  character(len=40) :: dim_unlim_name ! name of the unlimited dimension in the file
  character(len=1024) :: filename_temp
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=48), dimension(3) :: dim_names !< variable dimension names
  integer, dimension(3) :: dim_lengths !< variable dimension lengths

  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index <= 0) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif
  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, trim(filename_temp), lowercase(trim(mode)), is_restart=.false.)

  dim_names(:) = ""
  dim_lengths(:) = 0
  num_dims=0
  dim_unlim_index=0
  if (variable_exists(fileobj, trim(fieldname))) then
    dim_unlim_index = get_variable_unlimited_dimension_index(fileobj, trim(fieldname))
  else
    ! get the dimension names and lengths
    if (present(G)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, G=G)
    elseif(present(dG)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, dG=dG)
    endif

    if (present(GV)) &
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, GV=GV)

    call register_field(fileObj, trim(fieldname), "double", dimensions=dim_names(1:num_dims))
    call register_variable_attribute(fileObj, trim(fieldname), 'units', trim(var_desc%units))
    call register_variable_attribute(fileObj, trim(fieldname), 'long_name', trim(var_desc%longname))
    ! write the checksum attribute
    if (present(checksums)) then
      ! convert the checksum to a string
      checksum_char = ''
      checksum_char = convert_checksum_to_string(checksums(1,1))
      call register_variable_attribute(fileObj, trim(fieldname), "checksum", checksum_char)
    endif
  endif

  ndims=2
  start(:) = 1
  nwrite(:) = (/size(data,1), size(data,2)/)
  if (present(corner)) then
    do i=1,ndims
      if (corner(i) .le. nwrite(i)) then
        if (corner(i) .gt. 0) then
          if (i .ne. dim_unlim_index) start(i) = corner(i)
        endif
      endif
    enddo
  endif

  if (present(edge_lengths)) then
    do i=1,ndims
      if (edge_lengths(i) .le. nwrite(i)) then
        if (edge_lengths(i) .gt. 0) then
          if (i .ne. dim_unlim_index) nwrite(i) = edge_lengths(i)
        endif
      endif
    enddo
  endif

  allocate(data_tmp(size(data,1), size(data,2)))
  data_tmp = data
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data_tmp,scale)
  endif ; endif


  time_index=1
  dim_unlim_size=0
  dim_unlim_name=""
  ! write the data
  if (lowercase(trim(mode)) .eq. "append") then
    if (present(time_level)) then
      ! write the time value if it is not already written to the file
      file_time = read_most_recent_time(fileobj)
      if (time_level .gt. file_time+EPSILON(time_level)) then
        call get_unlimited_dimension_name(fileobj,dim_unlim_name)
        call get_dimension_size(fileobj, trim(dim_unlim_name), dim_unlim_size)
        if (dim_unlim_size .lt. unlimited) time_index=dim_unlim_size+1
        call write_data(fileobj, trim(dim_unlim_name), (/time_level/), corner=(/time_index/),edge_lengths=(/1/))
      endif
      if (dim_unlim_index .gt. 0) &
          call write_data(fileobj, trim(fieldname), data_tmp, corner=start, edge_Lengths=nwrite, &
                          unlim_dim_level=time_index)
    endif
  else
    call write_data(fileobj, trim(fieldname), data_tmp, corner=start, edge_Lengths=nwrite)
  endif
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)
  deallocate(data_tmp)
end subroutine write_field_2d_noDD

!> This function uses the fms_io function write_data to write a 3-D non-domain-decomposed data field named "fieldname"
!! to the file "filename" in "write", "overwrite", or "append" mode. It should be called after create_file in the MOM
!! file write procedure.
subroutine write_field_3d_noDD(filename, fieldname, data, mode, var_desc, &
                               corner, edge_lengths, time_level, scale, checksums, G, dG, GV)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:), intent(in) :: data !< The 3-dimensional data array to pass to read_data
  character(len=*), intent(in) :: mode !< "write", "overwrite", or "append"
  type(vardesc), optional, intent(in) :: var_desc !< structure describing variable output
  integer, dimension(3), optional, intent(in) :: corner !< starting index of data buffer. Default is 1
  integer, dimension(3), optional, intent(in) :: edge_lengths !< number of data values to read in; default is the
                                                          !! variable size
  real, optional, intent(in) :: time_level !< time value to write
  real, optional, intent(in) :: scale !< A scaling factor that the fields are multiplied by before they are written.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums  !< variable checksum
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  ! local
  type(FmsNetcdfFile_t) :: fileobj ! netCDF file object returned by call to open_file
  logical :: file_open_success !.true. if call to open_file is successful
  real :: file_time ! most recent time currently written to file
  real, allocatable, dimension(:,:,:) :: data_tmp
  integer :: i, ndims, num_dims, time_index, substring_index
  integer :: dim_unlim_size, dim_unlim_index ! size and dimension index of the unlimited dimension
  integer, dimension(3) :: start, nwrite ! indices for first data value and number of values to write
  character(len=40) :: dim_unlim_name ! name of the unlimited dimension in the file
  character(len=1024) :: filename_temp
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=48), dimension(4) :: dim_names !< variable dimension names
  integer, dimension(4) :: dim_lengths !< variable dimension lengths

  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index <= 0) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif
  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, trim(filename_temp), lowercase(trim(mode)), is_restart=.false.)

  dim_names(:) = ""
  dim_lengths(:) = 0
  num_dims=0
  dim_unlim_index=0
  if (variable_exists(fileobj, trim(fieldname))) then
    dim_unlim_index = get_variable_unlimited_dimension_index(fileobj, trim(fieldname))
  else
    ! get the dimension names and lengths
    if (present(G)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, G=G)
    elseif(present(dG)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, dG=dG)
    endif
    
    if (present(GV)) &
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, GV=GV)
    call register_field(fileObj, trim(fieldname), "double", dimensions=dim_names(1:num_dims))
    call register_variable_attribute(fileObj, trim(fieldname), 'units', trim(var_desc%units))
    call register_variable_attribute(fileObj, trim(fieldname), 'long_name', trim(var_desc%longname))
    ! write the checksum attribute
    if (present(checksums)) then
      ! convert the checksum to a string
      checksum_char = ""
      checksum_char = convert_checksum_to_string(checksums(1,1))
      call register_variable_attribute(fileObj, trim(fieldname), "checksum", checksum_char)
    endif
  endif

  ndims = 3
  start(:) = 1
  nwrite(:) = (/size(data,1), size(data,2), size(data,3)/)
  if (present(corner)) then
    do i=1,ndims
      if (corner(i) .le. nwrite(i)) then
        if (corner(i) .gt. 0) then
          if (i .ne. dim_unlim_index) start(i) = corner(i)
        endif
      endif
    enddo
  endif

  if (present(edge_lengths)) then
    do i=1,ndims
      if (edge_lengths(i) .le. nwrite(i)) then
        if (edge_lengths(i) .gt. 0) then
          if (i .ne. dim_unlim_index) nwrite(i) = edge_lengths(i)
        endif
      endif
    enddo
  endif

  allocate(data_tmp(size(data,1), size(data,2), size(data,3)))
  data_tmp = data
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data_tmp,scale)
  endif ; endif

  time_index = 1
  dim_unlim_size = 0
  dim_unlim_name = ""
  ! write the data
  if (lowercase(trim(mode)) .eq. "append") then
    if (present(time_level)) then
      ! write the time value if it is not already written to the file
      file_time = read_most_recent_time(fileobj)
      if (time_level .gt. file_time+EPSILON(time_level)) then
        call get_unlimited_dimension_name(fileobj,dim_unlim_name)
        call get_dimension_size(fileobj, trim(dim_unlim_name), dim_unlim_size)
        if (dim_unlim_size .lt. unlimited) time_index=dim_unlim_size+1

        call write_data(fileobj, trim(dim_unlim_name), (/time_level/), corner=(/time_index/),edge_lengths=(/1/))
      endif
      if (dim_unlim_index .gt. 0) &
          call write_data(fileobj, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite, &
                          unlim_dim_level=time_index)
    endif
  else
    call write_data(fileobj, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite)
  endif
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)
  deallocate(data_tmp)
end subroutine write_field_3d_noDD

!> This function uses the fms_io function write_data to write a 4-D non-domain-decomposed data field named "fieldname"
!! to the file "filename" in "write", "overwrite", or "append" mode. It should be called after create_file in the MOM
!! file write procedure.
subroutine write_field_4d_noDD(filename, fieldname, data, mode, var_desc, &
                               corner, edge_lengths, time_level, scale, checksums, G, dG, GV)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:,:), intent(in) :: data !< The 4-dimensional data array to pass to read_data
  character(len=*), intent(in) :: mode !< "write", "overwrite", or "append"
  type(vardesc), optional, intent(in) :: var_desc !< structure describing variable output
  integer, dimension(4), optional, intent(in) :: corner !< starting index of data buffer. Default is 1
  integer, dimension(4), optional, intent(in) :: edge_lengths !< number of data values to read in; default is the
                                                          !! variable size
  real, optional, intent(in) :: time_level !< time value to write
  real, optional, intent(in) :: scale !< A scaling factor that the fields are multiplied by before they are written.
  integer(kind=8), dimension(:,:), optional, intent(in) :: checksums  !< variable checksum
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  ! local
  type(FmsNetcdfFile_t) :: fileobj ! netCDF file object returned by call to open_file
  logical :: file_open_success !.true. if call to open_file is successful
  real :: file_time ! most recent time currently written to file
  real, allocatable, dimension(:,:,:,:) :: data_tmp
  integer :: i, ndims, num_dims, time_index, substring_index
  integer :: dim_unlim_size, dim_unlim_index ! size and dimension index of the unlimited dimension
  integer, dimension(4) :: start, nwrite ! indices for first data value and number of values to write
  character(len=40) :: dim_unlim_name ! name of the unlimited dimension in the file
  character(len=1024) :: filename_temp
  character(len=64) :: checksum_char ! checksum character array created from checksum argument
  character(len=48), dimension(4) :: dim_names !< variable dimension names
  integer, dimension(4) :: dim_lengths !< variable dimension lengths

  ! append '.nc' to the file name if it is missing
  filename_temp = ""
  substring_index = 0
  substring_index = index(trim(filename), ".nc")
  if (substring_index <= 0) then
    filename_temp = append_substring(filename,".nc")
  else
    filename_temp = filename
  endif
  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, trim(filename_temp), lowercase(trim(mode)), is_restart=.false.)

  dim_names(:) = ""
  dim_lengths(:) = 0
  num_dims=0
  dim_unlim_index=0
  if (variable_exists(fileobj, trim(fieldname))) then
    dim_unlim_index = get_variable_unlimited_dimension_index(fileobj, trim(fieldname))
  else
    ! get the dimension names and lengths
    if (present(G)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, G=G)
    elseif(present(dG)) then
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, dG=dG)
    endif
    if (present(GV)) &
      call get_var_dimension_features(var_desc%hor_grid, var_desc%z_grid, var_desc%t_grid, dim_names, &
                                      dim_lengths, num_dims, GV=GV)
    
    call register_field(fileObj, trim(fieldname), "double", dimensions=dim_names(1:num_dims))
    call register_variable_attribute(fileObj, trim(fieldname), 'units', trim(var_desc%units))
    call register_variable_attribute(fileObj, trim(fieldname), 'long_name', trim(var_desc%longname))
    ! write the checksum attribute
    if (present(checksums)) then
      ! convert the checksum to a string
      checksum_char = ''
      checksum_char = convert_checksum_to_string(checksums(1,1))
      call register_variable_attribute(fileObj, trim(fieldname), "checksum", checksum_char)
    endif
  endif

  ndims = 4
  start(:) = 1
  nwrite(:) = (/size(data,1), size(data,2), size(data,3), size(data,4)/)
  if (present(corner)) then
    do i=1,ndims
      if (corner(i) .le. nwrite(i)) then
        if (corner(i) .gt. 0) then
          if (i .ne. dim_unlim_index) start(i) = corner(i)
        endif
      endif
    enddo
  endif

  if (present(edge_lengths)) then
    do i=1,ndims
      if (edge_lengths(i) .le. nwrite(i)) then
        if (edge_lengths(i) .gt. 0) then
          if (i .ne. dim_unlim_index) nwrite(i) = edge_lengths(i)
        endif
      endif
    enddo
  endif

  allocate(data_tmp(size(data,1), size(data,2), size(data,3), size(data,4)))
  data_tmp = data
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data_tmp,scale)
  endif ; endif

  time_index=1
  dim_unlim_size=0
  dim_unlim_name=""

  ! write the data
  if (lowercase(trim(mode)) .eq. "append") then
    if (present(time_level)) then
      ! write the time value if it is not already written to the file
      file_time = read_most_recent_time(fileobj)
      if (time_level .gt. file_time+EPSILON(time_level)) then
        call get_unlimited_dimension_name(fileobj,dim_unlim_name)
        call get_dimension_size(fileobj, trim(dim_unlim_name), dim_unlim_size)
        if (dim_unlim_size .lt. unlimited) time_index=dim_unlim_size+1

        call write_data(fileobj, trim(dim_unlim_name), (/time_level/), corner=(/time_index/),edge_lengths=(/1/))
      endif

      if (dim_unlim_index .gt. 0) &
          call write_data(fileobj, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite, &
                          unlim_dim_level=time_index)
    endif
  else
    call write_data(fileobj, trim(fieldname), data_tmp, corner=start, edge_lengths=nwrite)
  endif
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)
  deallocate(data_tmp)
end subroutine write_field_4d_noDD

!> This function uses the fms_io function read_data to read 1-D domain-decomposed data field named "fieldname"
!! from file "filename".
subroutine MOM_read_data_1d_DD(filename, fieldname, data, domain, corner, edge_lengths, time_level, scale, &
                               x_position, y_position, x_units, y_units)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:), intent(inout) :: data !< The 1-dimensional data array to pass to read_data
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  integer, dimension(1), optional, intent(in) :: corner !< starting index of data buffer. Default is 1
  integer, dimension(1), optional, intent(in) :: edge_lengths !< number of data values to read in; default is the
                                                             !! variable size
  integer, optional, intent(in) :: time_level !< time level to read
  real, optional, intent(in) :: scale !< A scaling factor that the field is multiplied by
  integer, intent(in), optional :: x_position !< domain position of x-dimension; CENTER (default) or EAST_FACE
  integer, intent(in), optional :: y_position !< domain position of y-dimension; CENTER (default) or NORTH_FACE
  character(len=*), intent(in), optional :: x_units !< x-dimension units; default is "degrees_east"
  character(len=*), intent(in), optional :: y_units !< y-dimension units; default is "degrees_north"
  ! local
  type(FmsNetcdfDomainFile_t) :: fileobj ! netCDF file object returned by call to open_file
  logical :: file_open_success !.true. if call to open_file is successful
  integer :: i
  integer, dimension(1) :: start, nread ! indices for first data value and number of values to read
  character(len=64) :: xunits, yunits ! x- and y-dimension units
  character(len=40), dimension(1) :: dim_names ! variable dimension names
  integer :: xpos, ypos ! x and y domain positions
  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, filename, "read", domain%mpp_domain, is_restart=.false.)
  ! register the variable axes
  !> \note: the user will need to change the xUnits and yUnits if they expect different values for the
  !! x/longitude and/or y/latitude axes units
  xpos = CENTER
  ypos = CENTER
  if (present(x_position)) xpos = x_position
  if (present(y_position)) ypos = y_position

  xunits=""
  yunits=""
  if (present(x_units)) then
    xunits(1:len_trim(x_units)) = x_units
  else
    xunits = "degrees_east"
  endif
  if (present(y_units)) then
    yunits(1:len_trim(y_units)) = y_units
  else
    yunits = "degrees_north"
  endif

  call MOM_register_variable_axes(fileobj, trim(fieldname), xUnits=trim(xunits), yUnits=trim(yunits), &
                                   xPosition=xpos, yPosition=ypos)
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  call get_variable_dimension_names(fileobj, trim(fieldname), dim_names)

  start(1) = 1
  if (present(corner)) start(1) = corner(1)
  if (present(edge_lengths)) then
    nread(1) = edge_lengths(1)
  else
    call get_dimension_size(fileobj, trim(dim_names(1)), nread(1))
  endif
  ! read the data
  call read_data(fileobj, trim(fieldname), data, corner=start, edge_Lengths=nread)
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)

  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif

end subroutine MOM_read_data_1d_DD

!> This function uses the fms_io function read_data to read 2-D domain-decomposed data field named "fieldname"
!! from file "filename".
subroutine MOM_read_data_2d_DD(filename, fieldname, data, domain, corner, edge_lengths, time_level, scale, &
                               x_position, y_position, x_units, y_units)
  character(len=*), intent(in) :: filename  !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:), intent(inout) :: data !< The 2-dimensional data array to pass to read_data
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  integer, dimension(2), optional, intent(in) :: corner !< starting indices of data buffer. Default is 1
  integer, dimension(2), optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: time_level !< time level to read
  real, optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  integer, intent(in), optional :: x_position !< domain position of x-dimension; CENTER (default) or EAST_FACE
  integer, intent(in), optional :: y_position !< domain position of y-dimension; CENTER (default) or NORTH_FACE
  character(len=*), intent(in), optional :: x_units !< x-dimension units; default is "degrees_east"
  character(len=*), intent(in), optional :: y_units !< y-dimension units; default is "degrees_north"
  ! local
  type(FmsNetcdfDomainFile_t) :: fileobj !netCDF file object returned by call to open_file
  logical :: file_open_success !.true. if call to open_file is successful
  integer :: i, dim_unlim_index
  integer, dimension(2) :: start, nread ! indices for first data value and number of values to read
  character(len=64) :: xunits, yunits ! x- and y-dimension units
  character(len=40), dimension(2) :: dim_names ! variable dimension names
  integer :: xpos, ypos ! x and y domain positions

  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, filename, "read", domain%mpp_domain, is_restart=.false.)
  ! register the variable axes
  !> \note: the user will need to change the xUnits and yUnits if they expect different values for the
  !! x/longitude and/or y/latitude axes units
  xpos = CENTER
  ypos = CENTER
  if (present(x_position)) xpos = x_position
  if (present(y_position)) ypos = y_position
  xunits=""
  yunits=""
  if (present(x_units)) then
    xunits(1:len_trim(x_units)) = x_units
  else
    xunits = "degrees_east"
  endif
  if (present(y_units)) then
    yunits(1:len_trim(y_units)) = y_units
  else
    yunits = "degrees_north"
  endif
  call MOM_register_variable_axes(fileobj, trim(fieldname), xUnits=trim(xunits), yUnits=trim(yunits), &
                                  xPosition=xpos, yPosition=ypos)

  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  if (present(corner) .or. present(edge_lengths) .or. present(time_level)) then
    call get_variable_dimension_names(fileobj, trim(fieldname), dim_names)
  endif

  start(:) = 1
  if (present(corner)) start = corner

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    do i=1,2
      call get_dimension_size(fileobj, trim(dim_names(i)), nread(i))
    enddo
  endif

  if (present(time_level)) then
    dim_unlim_index=0
    do i=1,2
      if (is_dimension_unlimited(fileobj,dim_names(i))) then
        dim_unlim_index=i
        start(i)=time_level
        nread(i)=1
      endif
    enddo
    if (dim_unlim_index .LE. 0) &
      call MOM_error(FATAL, "MOM_io::MOM_read_data_2d_DD: time level specified, but variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
  endif
  ! read the data
  call read_data(fileobj, trim(fieldname), data, corner=start, edge_Lengths=nread)
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale, domain)
  endif ; endif

end subroutine MOM_read_data_2d_DD

!> This function uses the fms_io function read_data to read 3-D domain-decomposed data field named "fieldname"
!! from file "filename".
subroutine MOM_read_data_3d_DD(filename, fieldname, data, domain, corner, edge_lengths, time_level, scale, &
                               x_position, y_position, x_units, y_units)
  character(len=*), intent(in) :: filename  !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:), intent(inout) :: data !< The 3-dimensional data array to pass to read_data
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  integer, dimension(3), optional, intent(in) :: corner !< starting indices of data buffer. Default is 1
  integer, dimension(3), optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                             !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: time_level !< time level to read
  real,         optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  integer, intent(in), optional :: x_position !< domain position of x-dimension; CENTER (default) or EAST_FACE
  integer, intent(in), optional :: y_position !< domain position of y-dimension; CENTER (default) or NORTH_FACE
  character(len=*), intent(in), optional :: x_units !< x-dimension units; default is "degrees_east"
  character(len=*), intent(in), optional :: y_units !< y-dimension units; default is "degrees_north"
  ! local
  type(FmsNetcdfDomainFile_t) :: fileobj !netCDF file object returned by call to open_file
  logical :: file_open_success !.true. if call to open_file is successful
  integer :: i, dim_unlim_index
  integer, dimension(3) :: start, nread ! indices for first data value and number of values to read
  character(len=64) :: xunits, yunits ! x- and y-dimension units
  character(len=40), dimension(3) :: dim_names ! variable dimension names
  integer :: xpos, ypos ! x and y domain positions

  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, filename, "read", domain%mpp_domain, is_restart=.false.)
  ! register the variable axes
  !> \note: the user will need to change the xUnits and yUnits if they expect different values for the
  !! x/longitude and/or y/latitude axes units
  xpos = CENTER
  ypos = CENTER
  if (present(x_position)) xpos = x_position
  if (present(y_position)) ypos = y_position
  xunits=""
  yunits=""
  if (present(x_units)) then
    xunits(1:len_trim(x_units)) = x_units
  else
    xunits = "degrees_east"
  endif
  if (present(y_units)) then
    yunits(1:len_trim(y_units)) = y_units
  else
    yunits = "degrees_north"
  endif

  call MOM_register_variable_axes(fileobj, trim(fieldname), xUnits=trim(xunits), yUnits=trim(yunits), &
                                  xPosition=xpos, yPosition=ypos)
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  if (present(corner) .or. present(edge_lengths) .or. present(time_level)) then
    call get_variable_dimension_names(fileobj, trim(fieldname), dim_names)
  endif

  start(:) = 1
  if (present(corner)) start = corner

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    do i=1,3
      call get_dimension_size(fileobj, trim(dim_names(i)), nread(i))
    enddo
  endif

  if (present(time_level)) then
    dim_unlim_index=0
    do i=1,3
      if (is_dimension_unlimited(fileobj,dim_names(i))) then
        dim_unlim_index=i
        start(i)=time_level
        nread(i)=1
      endif
    enddo
    if (dim_unlim_index .LE. 0) &
      call MOM_error(FATAL, "MOM_io::MOM_read_data_3d_DD: time level specified, but variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
  endif
  ! read the data
  call read_data(fileobj, trim(fieldname), data, corner=start, edge_Lengths=nread)
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale, domain)
  endif ; endif

end subroutine MOM_read_data_3d_DD

!> This function uses the fms_io function read_data to read 4-D domain-decomposed data field named "fieldname"
!! from file "filename".
subroutine MOM_read_data_4d_DD(filename, fieldname, data, domain, corner, edge_lengths, time_level, scale, &
                               x_position, y_position, x_units, y_units)
  character(len=*),       intent(in) :: filename  !< The name of the file to read
  character(len=*),       intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:,:), intent(inout) :: data !< The 4-dimensional data array to pass to read_data
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  integer, dimension(4),  optional, intent(in) :: corner !< starting indices of data buffer. Default is 1
  integer, dimension(4),  optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: time_level !< time level to read
  real,         optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  integer, intent(in), optional :: x_position !< domain position of x-dimension; CENTER (default) or EAST_FACE
  integer, intent(in), optional :: y_position !< domain position of y-dimension; CENTER (default) or NORTH_FACE
  character(len=*), intent(in), optional :: x_units !< x-dimension units; default is "degrees_east"
  character(len=*), intent(in), optional :: y_units !< y-dimension units; default is "degrees_north"

  ! local
  type(FmsNetcdfDomainFile_t) :: fileobj !netCDF file object returned by call to open_file
  logical :: file_open_success !.true. if call to open_file is successful
  integer :: i, dim_unlim_index
  integer, dimension(4) :: start, nread ! indices for first data value and number of values to read
  character(len=40), dimension(4) :: dim_names ! variable dimension names
  character(len=64) :: xunits, yunits ! x- and y-dimension units
  integer :: xpos, ypos ! x and y domain positions

  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, filename, "read", domain%mpp_domain, is_restart=.false.)
  ! register the variable axes
  !> \note: the user will need to change the xUnits and yUnits if they expect different values for the
  !! x/longitude and/or y/latitude axes units
  xpos = CENTER
  ypos = CENTER
  if (present(x_position)) xpos = x_position
  if (present(y_position)) ypos = y_position
  xunits=""
  yunits=""
  if (present(x_units)) then
    xunits(1:len_trim(x_units)) = x_units
  else
    xunits = "degrees_east"
  endif
  if (present(y_units)) then
    yunits(1:len_trim(y_units)) = y_units
  else
    yunits = "degrees_north"
  endif

  call MOM_register_variable_axes(fileobj, trim(fieldname), xUnits=trim(xunits), yUnits=trim(yunits), &
                                  xPosition=xpos, yPosition=ypos)
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  if (present(corner) .or. present(edge_lengths) .or. present(time_level)) then
    call get_variable_dimension_names(fileobj, trim(fieldname), dim_names)
  endif

  start(:) = 1
  if (present(corner)) start = corner

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    do i=1,4
      call get_dimension_size(fileobj, trim(dim_names(i)), nread(i))
    enddo
  endif

  if (present(time_level)) then
    dim_unlim_index=0
    do i=1,4
      if (is_dimension_unlimited(fileobj,dim_names(i))) then
        dim_unlim_index=i
        start(i)=time_level
        nread(i)=1
      endif
    enddo
    if (dim_unlim_index .LE. 0) &
      call MOM_error(FATAL, "MOM_io::MOM_read_data_4d_DD: time level specified, but variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
  endif
  ! read the data
  call read_data(fileobj, trim(fieldname), data, corner=start, edge_Lengths=nread)
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)
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
  type(FmsNetcdfFile_t) :: fileobj ! netCDF file object returned by call to open_file
  logical :: file_open_success !.true. if call to open_file is successful
  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, filename, "read", is_restart=.false.)
  ! read the data
  call read_data(fileobj, trim(fieldname), data)
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)

end subroutine MOM_read_data_scalar

!> This function uses the fms_io function read_data to read 1-D data field named "fieldname"
!! from file "filename".
subroutine MOM_read_data_1d_noDD(filename, fieldname, data, corner, edge_lengths, time_level, scale)
  character(len=*),       intent(in) :: filename !< The name of the file to read
  character(len=*),       intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:),     intent(inout) :: data !< The 1-dimensional data array to pass to read_data
  integer, dimension(1), optional, intent(in) :: corner !< starting index of data buffer. Default is 1
  integer, dimension(1), optional, intent(in) :: edge_lengths !< number of data values to read in; default is
                                                             !! the variable size
  integer,      optional, intent(in) :: time_level !< time level to read
  real,         optional, intent(in) :: scale !< A scaling factor that the field is multiplied by
  ! local
  type(FmsNetcdfFile_t) :: fileobj ! netCDF file object returned by call to open_file
  logical :: file_open_success !.true. if call to open_file is successful
  integer :: i, ndims
  integer, allocatable, dimension(:) :: start, nread ! indices for first data value and number of values to read
  character(len=40), allocatable, dimension(:) :: dim_names ! variable dimension names

  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, filename, "read", is_restart=.false.)

  ndims = get_variable_num_dimensions(fileobj,trim(fieldname))
  allocate(start(ndims))
  allocate(nread(ndims))
  allocate(dim_names(ndims))

  if (present(corner)) start(1) = corner(1)
  if (present(edge_lengths)) then
    nread(1) = edge_lengths(1)
  else
    call get_variable_dimension_names(fileobj, trim(fieldname), dim_names)
    call get_dimension_size(fileobj, dim_names(1), nread(1))
  endif
  ! read the data
  call read_data(fileobj, trim(fieldname), data, corner=start, edge_Lengths=nread)
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif

  deallocate(start)
  deallocate(nread)
  deallocate(dim_names)
end subroutine MOM_read_data_1d_noDD

!> This function uses the fms_io function read_data to read 2-D data field named "fieldname"
!! from file "filename".
subroutine MOM_read_data_2d_noDD(filename, fieldname, data, corner, edge_lengths, time_level, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:),   intent(inout) :: data !< The 2-dimensional data array to pass to read_data
  integer, dimension(2),  optional, intent(in) :: corner !< starting indices of data buffer. Default is 1
  integer, dimension(2),  optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: time_level !< time level to read
  real, optional, intent(in):: scale !< A scaling factor that the field is multiplied by

  ! local
  type(FmsNetcdfFile_t) :: fileobj !netCDF file object returned by call to open_file
  logical :: file_open_success !.true. if call to open_file is successful
  integer :: i, dim_unlim_index
  integer, dimension(2) :: start, nread ! indices for first data value and number of values to read
  character(len=40), dimension(2) :: dim_names ! variable dimension names

  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, filename, "read", is_restart=.false.)
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  if (present(corner) .or. present(edge_lengths) .or. present(time_level)) then
    call get_variable_dimension_names(fileobj, trim(fieldname), dim_names)
  endif

  start(:) = 1
  if (present(corner)) start = corner

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    do i=1,2
      call get_dimension_size(fileobj, trim(dim_names(i)), nread(i))
    enddo
  endif

  if (present(time_level)) then
    dim_unlim_index=0
    do i=1,2
      if (is_dimension_unlimited(fileobj,dim_names(i))) then
        dim_unlim_index=i
        start(i)=time_level
        nread(i)=1
      endif
    enddo
    if (dim_unlim_index .LE. 0) &
      call MOM_error(FATAL, "MOM_io::MOM_read_data_2d_noDD: time level specified, but variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
  endif

  ! read the data
  call read_data(fileobj, trim(fieldname), data, corner=start, edge_Lengths=nread)
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif

end subroutine MOM_read_data_2d_noDD

!> This function uses the fms_io function read_data to read 3-D data field named "fieldname"
!! from file "filename".
subroutine MOM_read_data_3d_noDD(filename, fieldname, data, corner, edge_lengths, time_level, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:), intent(inout) :: data !< The 3-dimensional data array to pass to read_data
  integer, dimension(3),  optional, intent(in) :: corner !< starting indices of data buffer. Default is 1
  integer, dimension(3),  optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: time_level !< time level to read
  real,         optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  ! local
  type(FmsNetcdfFile_t) :: fileobj !netCDF file object returned by call to open_file
  logical :: file_open_success !.true. if call to open_file is successful
  integer :: i, dim_unlim_index
  integer, dimension(3) :: start, nread ! indices for first data value and number of values to read
  character(len=40), dimension(3) :: dim_names ! variable dimension names

  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, filename, "read", is_restart=.false.)
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  if (present(corner) .or. present(edge_lengths) .or. present(time_level)) then
    call get_variable_dimension_names(fileobj, trim(fieldname), dim_names)
  endif

  start(:) = 1
  if (present(corner)) start = corner

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    do i=1,3
      call get_dimension_size(fileobj, trim(dim_names(i)), nread(i))
    enddo
  endif

  if (present(time_level)) then
    dim_unlim_index=0
    do i=1,3
      if (is_dimension_unlimited(fileobj,dim_names(i))) then
        dim_unlim_index=i
        start(i)=time_level
        nread(i)=1
      endif
    enddo
    if (dim_unlim_index .LE. 0) &
      call MOM_error(FATAL, "MOM_io::MOM_read_data_3d_noDD: time level specified, but variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
  endif
  ! read the data
  call read_data(fileobj, trim(fieldname), data, corner=start, edge_Lengths=nread)
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif

end subroutine MOM_read_data_3d_noDD

!> This function uses the fms_io function read_data to read 3-D data field named "fieldname"
!! from file "filename".
subroutine MOM_read_data_4d_noDD(filename, fieldname, data, corner, edge_lengths, time_level, scale)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:,:),   intent(inout) :: data !< The 4-dimensional array to pass to read_data
  integer, dimension(4),  optional, intent(in) :: corner !< starting indices of data buffer. Default is 1
  integer, dimension(4),  optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: time_level !< time level to read
  real,         optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  ! local
  type(FmsNetcdfFile_t) :: fileobj !netCDF file object returned by call to open_file
  logical :: file_open_success !.true. if call to open_file is successful
  integer :: i, dim_unlim_index
  integer, dimension(4) :: start, nread ! indices for first data value and number of values to read
  character(len=40), dimension(4) :: dim_names ! variable dimension names

  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, filename, "read", is_restart=.false.)
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  if (present(corner) .or. present(edge_lengths) .or. present(time_level)) then
    call get_variable_dimension_names(fileobj, trim(fieldname), dim_names)
  endif

  start(:) = 1
  if (present(corner)) start = corner

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    do i=1,4
      call get_dimension_size(fileobj, trim(dim_names(i)), nread(i))
    enddo
  endif

  if (present(time_level)) then
    dim_unlim_index=0
    do i=1,4
      if (is_dimension_unlimited(fileobj,dim_names(i))) then
        dim_unlim_index=i
        start(i)=time_level
        nread(i)=1
      endif
    enddo
    if (dim_unlim_index .LE. 0) &
      call MOM_error(FATAL, "MOM_io::MOM_read_data_4d_noDD: time level specified, but variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
  endif
  ! read the data
  call read_data(fileobj, trim(fieldname), data, corner=start, edge_Lengths=nread)
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif

end subroutine MOM_read_data_4d_noDD

!> This is a kluge interface to obtain the correct x and y values used to define the diagnostic axes
!! lath, latq, lonh, and lonq in MOM_grid_intialize. This routine allocates a buffer of the same size as the field to
!! read. If reading in 'y' (latitude), it searches the first dimension for the point that contains the full range
!! of values used to define "lath" and "latq".
subroutine MOM_read_data_2d_noDD_diag_axes(filename, fieldname, data, define_diagnostic_axes, G, corner, edge_lengths, &
                                           time_level, scale, grid_type)
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:),   intent(inout) :: data !< The 2-dimensional data array to pass to read_data
  logical, intent(in) :: define_diagnostic_axes !< if .true., read in the full data array, search for the
                                                   !! full ranges of x/lon and y/lat values needed to define
                                                   !! the diagnostic axes
  type(dyn_horgrid_type), optional, intent(in) :: G !< The dynamic horizontal grid type; required if
                                                    !! define_diagnostic_axes=.true.
  integer, dimension(2),  optional, intent(in) :: corner !< starting indices of data buffer. Default is 1
  integer, dimension(2),  optional, intent(in) :: edge_lengths !< number of data values to read in.
  !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: time_level !< time level to read
  real, optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  character(len=*), optional, intent(in) :: grid_type !< grid type; either 'b' or 't'
  ! local
  type(FmsNetcdfFile_t) :: fileobj !netCDF file object returned by call to open_file
  logical :: file_open_success !.true. if call to open_file is successful
  logical, dimension(:), allocatable :: yuse! Mask for geoLatT and geoLatB values
  integer :: i, j, dim_unlim_index, xidx, yidx1, yidx2
  integer, dimension(2) :: start, nread, dimSizes ! indices for first data value and number of values to read,
                                                  ! variable dimension sizes
  character(len=40), dimension(2) :: dim_names ! variable dimension names
  character(len=40) :: units ! variable units
  character(len=1) :: gtype ! grid type
  real :: ymax1, ymax2 ! max values from latitude array
  real, allocatable, dimension(:,:) :: tmpGlbl ! temporary array to hold data

  ! register the global axes
  !do i=1,ndims

  !  call get_dimension_size(fileObjRead, dim_names(i), globalDimSize)

  !  if (globalDimSize .eq. size(tmpT,1)) then
  !    call register_axis(fileObjRead, trim(dim_names(i)),'x', domain_position=CENTER)
  !  elseif (globalDimSize .eq. size(tmpT,2)) then
  !    call register_axis(fileObjRead, trim(dim_names(i)),'y', domain_position=CENTER)
   ! elseif (globalDimSize .eq. size(tmpV,1)) then
  !    call register_axis(fileObjRead, trim(dim_names(i)),'x', domain_position=EAST_FACE)
   ! elseif (globalDimSize .eq. size(tmpU,2)) then
  !    call register_axis(fileObjRead, trim(dim_names(i)),'y', domain_position=NORTH_FACE)
   ! endif
 ! enddo

  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, filename, "read", is_restart=.false.)

  call get_variable_dimension_names(fileobj, trim(fieldname), dim_names)
  start(:) = 1
  if (present(corner)) start = corner

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    do i=1,2
      call get_dimension_size(fileobj, trim(dim_names(i)), nread(i))
    enddo
  endif

  if (present(time_level)) then
    dim_unlim_index=0
    do i=1,2
      if (is_dimension_unlimited(fileobj,dim_names(i))) then
        dim_unlim_index=i
        start(i)=time_level
        nread(i)=1
      endif
    enddo
    if (dim_unlim_index .LE. 0) &
      call MOM_error(FATAL, "MOM_io::MOM_read_data_2d_noDD_diag_axes: time level specified, but variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
 endif

 !> \note: the following indexing procedure searches for global x-dimension indices that correspond to the full ranges
 !! of y (geoLatT and geoLatB) values needed for the diagnostic indices (lath and latq).
 !! Defining the gridLatT and gridLatB values using the previous values in tmpGlbl(1,:) did not result in the correct
 !! diagnostic index values due to different indexing conventions for non-domain-decomposed IO in the new and previous
 !! procedures.

 if (define_diagnostic_axes) then
   if (.not.(present(grid_type))) call MOM_error(FATAL, "MOM_io::MOM_read_data_2d_noDD_diag_axes: grid_type"// &
      " argument required if define_diag_axes=.true.")
   gtype=""
   allocate(tmpGlbl(nread(1),nread(2)))
   tmpGlbl(:,:) = 999.0

   call get_variable_units(fileobj, fieldname, units)
   if (lowercase(trim(units)) .eq. "degrees_north") then
      ! create a mask for the T-grid latitude values in the tmpGlbl array
      allocate(yuse(nread(2)))
      yuse(:) = .FALSE.
      gtype = lowercase(trim(grid_type))
      select case(gtype)
        case ('t')
          do j=G%jsg,G%jeg
            yuse(2*(j-G%jsg)+2) = .TRUE.
          enddo
        case ('b')
          do j=G%jsg-1,G%jeg
            yuse(2*(j-G%jsg)+3) = .TRUE.
          enddo
        case default
          call MOM_error(FATAL, "MOM_io::MOM_read_data_2d_noDD_diag_axes: grid_type must be either t or b")
      end select

      ! get the index for the x-dimension with the maximum T-grid latitude (geoLatT) value in the tmpGlbl array
      xidx = 0
      yidx1 = 0
      ymax1 = 999.0
      yidx2 = 0
      ymax2 = 999.0

      do i=1,nread(1)
        ! find index of the maximum T-grid latitude value for the ith x (longitude) dimension
        yidx1 = MAXLOC(tmpGlbl(i,:), 1, MASK=yuse)
        ymax1 = tmpGLbl(i,yidx1)
        ! if the new max is greater than the current max, set the current max to the new max value
        if ( MAX(ymax1,ymax2) .EQ. ymax1) then
          yidx2 = yidx1
          ymax2 = ymax1
          xidx = i
        endif
      enddo

      if (xidx .LT. 1) call MOM_error(FATAL, "MOM_io::MOM_read_data_2d_nodDD_diag_axes: xidx is less than 1")
      if (ymax2 .GT. 90.0) &
      call MOM_error(FATAL, "MOM_io::MOM_read_data_2d_nodDD_diag_axes: ymax2 is greater than 90.0")
      ! read the data
      call read_data(fileobj, trim(fieldname), tmpGlbl, corner=start, edge_lengths=nread)
      data(1,1:nread(2)) = tmpGlbl(xidx,:)

      deallocate(yuse)

    elseif (lowercase(trim(units)) .eq. "degrees_north") then
        ! read the data
        call read_data(fileobj, trim(fieldname), tmpGlbl, corner=start, edge_Lengths=nread)
        ! all latitude indices contain the full range of x/longitude values required for the diagnostic axes
        data(1:nread(1),:) = tmpGlbl(:,1:size(data,2))
    endif
    deallocate(tmpGlbl)
  else ! just read the data into the user-specified data buffer
     call read_data(fileobj, trim(fieldname), data, corner=start, edge_Lengths=nread)
  endif

  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)
  ! scale the data
  if (present(scale)) then ; if (scale /= 1.0) then
    call scale_data(data, scale)
  endif ; endif

end subroutine MOM_read_data_2d_noDD_diag_axes

!> register a MOM diagnostic axis to a domain-decomposed file
subroutine MOM_register_diagnostic_axis(fileObj, axisName, axisLength)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObj !< netCDF file object returned by call to open_file
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
!> \note The user must specify units for variables with longitude/x-axis and/or latitude/y-axis axes to obtain
!! the correct domain decomposition for the data buffer.
subroutine MOM_register_variable_axes(fileObj, variableName, xUnits, yUnits, xPosition, yPosition)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObj !< netCDF file object returned by call to open_file
  character(len=*), intent(in) :: variableName !< name of the variable
  character(len=*), intent(in), optional :: xUnits !< x-axis (longitude) units to search for
  character(len=*), intent(in), optional :: yUnits !< y-axis (latitude) units to search for
  integer, intent(in), optional :: xPosition !< domain position of the x-axis
  integer, intent(in), optional :: yPosition !< domain position of the y-axis
  ! local
  character(len=40) :: units ! units corresponding to a specific variable dimension
  character(len=40), allocatable, dimension(:) :: dim_names ! variable dimension names
  integer :: i
  integer :: ndims ! number of dimensions
  integer :: xPos, yPos ! domain positions for x and y axes. Default is CENTER
  integer, allocatable, dimension(:) :: dimSizes ! variable dimension sizes

  if (.not. check_if_open(fileObj)) call MOM_error(FATAL,"MOM_io:register_variable_axes: The fileObj has "// &
                                                  "not been opened. Call open_file(fileObj,...) before "// &
                                                  "passing the fileObj argument to this function.")
  xPos=CENTER
  yPos=CENTER
  if (present(xPosition)) xPos=xPosition
  if (present(yPosition)) yPos=yPosition
 ! get variable dimension names and lengths
  ndims = get_variable_num_dimensions(fileObj, trim(variableName))
  allocate(dimSizes(ndims))
  allocate(dim_names(ndims))
  call get_variable_size(fileObj, trim(variableName), dimSizes, broadcast=.true.)
  call get_variable_dimension_names(fileObj, trim(variableName), dim_names)
  ! register the axes
  do i=1,ndims
    units=""
    if (present(xUnits)) then
      call get_variable_units(fileObj, dim_names(i), units)
      if (trim(lowercase(units)) .eq. trim(lowercase(xUnits))) then
        call register_axis(fileObj, dim_names(i),"x", domain_position=xPos)
      endif
    elseif (present(yUnits)) then
      call get_variable_units(fileObj, dim_names(i), units)
      if (trim(lowercase(units)) .eq. trim(lowercase(yUnits))) then
        call register_axis(fileObj, dim_names(i),"y", domain_position=yPos)
      endif
    else
      call register_axis(fileObj, dim_names(i), dimSizes(i))
    endif
  enddo

  deallocate(dimSizes)
  deallocate(dim_names)
end subroutine MOM_register_variable_axes

!> set the "start" (corner) and "nread" (edge_lengths) arrays for domain-decomposed netCDF input data buffers
subroutine MOM_get_nc_corner_edge_lengths_DD(fileObj, variableName, corner, edge_lengths, myCorner, myCornerIndices, &
                                            myedge_lengths, myEdgeLengthIndices)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObj !< netCDF file object returned by call to open_file
  character(len=*), intent(in), optional :: variableName !< name of the Varibl
  integer, allocatable, dimension(:), intent(out) :: corner !< array of corner indices to pass to read_data
  integer, allocatable, dimension(:), intent(out) :: edge_lengths !< array of edge_lengths indices to pass to read_data
  integer, dimension(:), intent(in), optional :: myCorner !< array of user-specified corner indices
  integer, dimension(:), intent(in), optional :: myCornerIndices !< positional indices for userCorner values
  integer, dimension(:), intent(in), optional :: myedge_lengths !< array of user-specified edge_lengths indices
  integer, dimension(:), intent(in), optional :: myEdgeLengthIndices !< positional indices for myedge_lengths
  ! local
  integer :: i, idx, ndims ! counter, index, number of dimensions

  if (.not. check_if_open(fileObj)) call MOM_error(FATAL,"MOM_io:MOM_get_nc_corner_edge_lengths_DD: "// &
                                                  "The fileObj has not been opened. Call open_file(fileObj,...)"// &
                                                  "before passing the fileObj argument to this function.")

 if (allocated(corner)) deallocate(corner)
 if (allocated(edge_lengths)) deallocate(edge_lengths)
 ! get variable dimension sizes and lengths
  ndims = get_variable_num_dimensions(fileObj, trim(variableName))
  allocate(corner(ndims))
  corner(:)=1
  allocate(edge_lengths(ndims))
  call get_variable_size(fileObj, trim(variableName), edge_lengths, broadcast=.true.)
  ! set user-specified corner values
  if (present(myCorner)) then
     if (.not.(present(myCornerIndices))) then
       call MOM_error(FATAL,"MOM_io:MOM_get_nc_corner_edge_lengths_DD: You passed the myCorner argument, but did "// &
                      "not pass myCornerIndices that defines the indices corresponding to the myCorner values.")
     endif
     do i=1,size(myCorner)
       idx = myCornerIndices(i)
       corner(idx)=myCorner(i)
     enddo
  endif

  ! set user-specified edge_lengths values
  if (present(myedge_lengths)) then
     if (.not.(present(myCornerIndices))) then
       call MOM_error(FATAL,"MOM_io:MOM_get_nc_corner_edge_lengths_DD:You passed the myedge_lengths argument, but "// &
                      "did not pass the myEdgeLengthIndices that defines the indices corresponding to the "// &
                      "myedge_lengths values.")
     endif
     do i=1,size(myedge_lengths)
       idx = myEdgeLengthIndices(i)
       edge_lengths(idx)=myedge_lengths(i)
     enddo
  endif

end subroutine MOM_get_nc_corner_edge_lengths_DD

!> set the corner (start) and edge_lengths (count) arrays for non-domain-decomposed netCDF input data buffers
subroutine MOM_get_nc_corner_edge_lengths_noDD(fileObj, variableName, corner, edge_lengths, myCorner, myCornerIndices, &
                                            myedge_lengths, myEdgeLengthIndices)
  type(FmsNetcdfFile_t), intent(inout) :: fileObj !< netCDF file object returned by call to open_file
  character(len=*), intent(in), optional :: variableName !< name of the Varibl
  integer, allocatable, dimension(:), intent(out) :: corner !< array of corner indices to pass to read_data
  integer, allocatable, dimension(:), intent(out) :: edge_lengths !< array of edge_lengths indices to pass to read_data
  integer, dimension(:), intent(in), optional :: myCorner !< array of user-specified corner indices
  integer, dimension(:), intent(in), optional :: myCornerIndices !< positional indices for userCorner values
  integer, dimension(:), intent(in), optional :: myedge_lengths !< array of user-specified edge_lengths indices
  integer, dimension(:), intent(in), optional :: myEdgeLengthIndices !< positional indices for myedge_lengths
  ! local
  integer :: i, idx, ndims ! counter, index, number of dimensions

  if (.not. check_if_open(fileObj)) call MOM_error(FATAL,"MOM_io:MOM_get_nc_corner_edge_lengths_DD: "// &
                                                  "The fileObj has not been opened. Call "// &
                                                  "open_file(fileObj,...) before passing the fileObj argument "// &
                                                  "to this function.")

 if (allocated(corner)) deallocate(corner)
 if (allocated(edge_lengths)) deallocate(edge_lengths)
 ! get variable dimension sizes and lengths
  ndims = get_variable_num_dimensions(fileObj, trim(variableName))
  allocate(corner(ndims))
  corner(:)=1
  allocate(edge_lengths(ndims))
  call get_variable_size(fileObj, trim(variableName), edge_lengths, broadcast=.true.)
  ! set user-specified corner values
  if (present(myCorner)) then
     if (.not.(present(myCornerIndices))) then
       call MOM_error(FATAL,"MOM_io:MOM_get_nc_corner_edge_lengths_noDD: You passed the myCorner argument, but did "// &
                            "not pass myCornerIndices that defines the indices corresponding to the myCorner values")
     endif
     do i=1,size(myCorner)
       idx = myCornerIndices(i)
       corner(idx)=myCorner(i)
     enddo
  endif

  ! set user-specified edge_lengths values
  if (present(myedge_lengths)) then
     if (.not.(present(myCornerIndices))) then
       call MOM_error(FATAL,"MOM_io:MOM_get_nc_corner_edge_lengths_noDD: You passed the myedge_lengths argument, "// &
                            "but did not pass myEdgeLengthIndices that defines the indices corresponding to the "// &
                            "myedge_lengths values.")
     endif
     do i=1,size(myedge_lengths)
       idx = myEdgeLengthIndices(i)
       edge_lengths(idx)=myedge_lengths(i)
     enddo
  endif

end subroutine MOM_get_nc_corner_edge_lengths_noDD

!> Get the horizontal grid, vertical grid, and/or time dimension names and lengths
!! for a single variable from the grid ids returned by a prior call to query_vardesc
subroutine get_var_dimension_features(hor_grid, z_grid, t_grid_in, &
                                      dim_names, dim_lengths, num_dims, G, dG, GV)

  character(len=*), intent(in) :: hor_grid !< horizontal grid
  character(len=*), intent(in) :: z_grid !< vertical grid
  character(len=*), intent(in) :: t_grid_in !< time grid
  character(len=*), dimension(:), intent(inout) :: dim_names !< array of dimension names
  integer, dimension(:), intent(inout) :: dim_lengths !< array of dimension sizes
  integer, intent(inout) :: num_dims !< number of axes to register in the restart file
  type(ocean_grid_type),  optional, intent(in) :: G !< The ocean's grid structure
  type(dyn_horgrid_type), optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
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
  !integer :: npes
  real, pointer, dimension(:) :: gridLatT => NULL(), & ! The latitude or longitude of T or B points for
     gridLatB => NULL(), & ! the purpose of labeling the output axes.
     gridLonT => NULL(), &
     gridLonB => NULL()
  type(MOM_domain_type), pointer :: domain => NULL() ! Domain used to get the pe count
  ! set the ocean grid coordinates

  !npes=0
  if (present(G)) then
    gridLatT => G%gridLatT ; gridLatB => G%gridLatB
    gridLonT => G%gridLonT ; gridLonB => G%gridLonB
    isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg
    IsgB = G%IsgB ; IegB = G%IegB ; JsgB = G%JsgB ; JegB = G%JegB
    domain => G%domain
    !npes = mpp_get_domain_npes(domain%mpp_domain)
  elseif (present(dG)) then
    gridLatT => dG%gridLatT ; gridLatB => dG%gridLatB
    gridLonT => dG%gridLonT ; gridLonB => dG%gridLonB
    isg = dG%isg ; ieg = dG%ieg ; jsg = dG%jsg ; jeg = dG%jeg
    IsgB = dG%IsgB ; IegB = dG%IegB ; JsgB = dG%JsgB ; JegB = dG%JegB
    domain => dG%domain
    !npes = mpp_get_domain_npes(domain%mpp_domain)
  endif

  use_lath = .false.
  use_lonh = .false.
  use_latq = .false.
  use_lonq = .false.

  call get_horizontal_grid_logic(hor_grid, use_lath, use_lonh, use_latq, use_lonq)
  ! add longitude name to dimension name array
  if (use_lonh) then
    num_dims = num_dims+1
    dim_names(num_dims)(1:len_trim('lonh')) = 'lonh'
    !if (npes .gt. 1) then
      dim_lengths(num_dims) = size(gridLonT(isg:ieg))
    !else
    !  dim_lengths(num_dims) = size(gridLonT)
    !endif

  elseif (use_lonq) then
    num_dims = num_dims+1
    dim_names(num_dims)(1:len_trim('lonq')) = 'lonq'
    !if (npes .gt. 1) then
      dim_lengths(num_dims) = size(gridLonB(IsgB:IegB))
    !else
    !  dim_lengths(num_dims) = size(gridLonB)
    !endif

  endif
  ! add latitude name to dimension name array
  if (use_lath) then
    num_dims = num_dims+1
    dim_names(num_dims)(1:len_trim('lath')) = 'lath'
    !if (npes .gt. 1) then
      dim_lengths(num_dims) = size(gridLatT(jsg:jeg))
    !else
    !  dim_lengths(num_dims) = size(gridLatT)
    !endif
  elseif (use_latq) then
    num_dims = num_dims+1
    dim_names(num_dims)(1:len_trim('latq')) = 'latq'
    !if (npes .gt. 1) then
      dim_lengths(num_dims) = size(gridLatB(JsgB:JegB))
    !else
    !  dim_lengths(num_dims) = size(gridLatB)
   ! endif
  endif

  if (associated(domain)) nullify(domain)
  ! vertical grid
  select case (trim(z_grid))
    case ('L')
      num_dims = num_dims+1
      dim_names(num_dims)(1:len_trim('Layer')) = 'Layer'
      dim_lengths(num_dims) = size(GV%sLayer(1:GV%ke))
    case ('i')
      num_dims = num_dims+1
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
      dim_names(num_dims)(1:len_trim('Time')) = 'Time'
      dim_lengths(num_dims) = unlimited
    case ('p')
      if (len_trim(t_grid(2:8)) <= 0) then
          call MOM_error(FATAL,"MOM_io:get_var_dimension_features: "//&
                           "No periodic axis length was specified in "//trim(t_grid))
      endif
      num_dims = num_dims+1
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
  integer :: npes
  real, pointer, dimension(:) :: gridLatT => NULL(), & ! The latitude or longitude of T or B points for
     gridLatB => NULL(), & ! the purpose of labeling the output axes.
     gridLonT => NULL(), &
     gridLonB => NULL()

  type(MOM_domain_type), pointer :: domain => NULL() ! Domain used to get the pe count 
  ! set the ocean grid coordinates

  npes=0
  if (present(G)) then
    gridLatT => G%gridLatT ; gridLatB => G%gridLatB
    gridLonT => G%gridLonT ; gridLonB => G%gridLonB
    x_axis_units = G%x_axis_units ; y_axis_units = G%y_axis_units
    isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg
    IsgB = G%IsgB ; IegB = G%IegB ; JsgB = G%JsgB ; JegB = G%JegB
    domain => G%domain

    npes = mpp_get_domain_npes(domain%mpp_domain)
  elseif (present(dG)) then
    gridLatT => dG%gridLatT ; gridLatB => dG%gridLatB
    gridLonT => dG%gridLonT ; gridLonB => dG%gridLonB
    x_axis_units = dG%x_axis_units ; y_axis_units = dG%y_axis_units
    isg = dG%isg ; ieg = dG%ieg ; jsg = dG%jsg ; jeg = dG%jeg
    IsgB = dG%IsgB ; IegB = dG%IegB ; JsgB = dG%JsgB ; JegB = dG%JegB
    domain => dG%domain

    npes = mpp_get_domain_npes(domain%mpp_domain)
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
      if (associated(gridLatT)) then 
        axis_data_CS%data(axis_number)%p=>gridLatT(jsg:jeg)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Latitude'
        axis_data_CS%axis(axis_number)%units = y_axis_units
        axis_data_CS%axis(axis_number)%horgrid_position = CENTER
        axis_data_CS%axis(axis_number)%is_domain_decomposed = .true.
      endif
    case('lonh')
      if (associated(gridLonT)) then
        axis_data_CS%data(axis_number)%p=>gridLonT(isg:ieg)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%horgrid_position  = CENTER
        axis_data_CS%axis(axis_number)%longname = 'Longitude'
        axis_data_CS%axis(axis_number)%units = x_axis_units
        axis_data_CS%axis(axis_number)%is_domain_decomposed  = .true.
      endif
    case('latq')
      if (associated(gridLatB)) then
        axis_data_CS%data(axis_number)%p=>gridLatB(JsgB:JegB)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Latitude'
        axis_data_CS%axis(axis_number)%units = y_axis_units
        axis_data_CS%axis(axis_number)%horgrid_position = NORTH_FACE
        axis_data_CS%axis(axis_number)%is_domain_decomposed = .true.
      endif
    case('lonq')
      if (associated(gridLonB)) then
        axis_data_CS%data(axis_number)%p=>gridLonB(IsgB:IegB)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Longitude'
        axis_data_CS%axis(axis_number)%units = x_axis_units
        axis_data_CS%axis(axis_number)%horgrid_position = EAST_FACE
        axis_data_CS%axis(axis_number)%is_domain_decomposed = .true.
      endif
    case('Layer')
      if (present(GV)) then
        axis_data_CS%data(axis_number)%p=>GV%sLayer(1:GV%ke)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Layer pseudo-depth, -z*'
        axis_data_CS%axis(axis_number)%units = GV%zAxisUnits
        axis_data_CS%axis(axis_number)%positive = 'up'
      endif
    case('Interface')
      if (present(GV)) then
        axis_data_CS%data(axis_number)%p=>GV%sInterface(1:GV%ke+1)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Interface pseudo-depth, -z*'
        axis_data_CS%axis(axis_number)%units = GV%zAxisUnits
        axis_data_CS%axis(axis_number)%positive = 'up'
      endif
    case('Time')
      if (.not.(present(time_val))) &
           call MOM_error(FATAL, "MOM_io::get_diagnostic_axis_data: requires time_val"//&
                          " and time_units arguments for "//trim(axis_name))

      axis_data_CS%data(axis_number)%p=>time_val
      axis_data_CS%axis(axis_number)%name = trim(axis_name)
      axis_data_CS%axis(axis_number)%longname = 'Time'

      if (present(time_units)) then
        axis_data_CS%axis(axis_number)%units = time_units
      else
        axis_data_CS%axis(axis_number)%units = 'days'
      endif
    case('Period')
      if (.not.(present(time_val))) &
        call MOM_error(FATAL, "MOM_io::get_diagnostic_axis_data: requires a time_val argument "// &
                       "for "//trim(axis_name))
      axis_data_CS%data(axis_number)%p=>time_val
      axis_data_CS%axis(axis_number)%name = trim(axis_name)
      axis_data_CS%axis(axis_number)%longname = 'Periods for cyclical variables'
    case default
      call MOM_error(WARNING, "MOM_io::get_diagnostic_axis_data:"//trim(axis_name)//"is an unrecognized axis")
  end select

  if (associated(domain)) nullify(domain)

end subroutine MOM_get_diagnostic_axis_data

!> check a netcdf file to see whether it contains the specified field
function field_exists(file_name, variable_name) result (var_exists)
  character(len=*),  intent(in) :: file_name   !< name of the file to check
  character(len=*),  intent(in) :: variable_name  !< name of the variable to check for

  ! local
  logical :: var_exists ! .true. if variable is found in file
  logical :: file_open_success ! .true. if call to open_file is successful
  type(fmsNetcdfFile_t) :: fileobj ! netcdf file object returned by open_file

  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, file_name, "read", is_restart=.false.)
  var_exists = variable_exists(fileobj, trim(variable_name))
  if (check_if_open(fileobj)) call close_file(fileobj)
end function field_exists

!> get the dimesion sizes of a variable
subroutine field_size(file_name, variable_name, dim_sizes)
  character(len=*),  intent(in) :: file_name   !< name of the file to check
  character(len=*),  intent(in) :: variable_name  !< name of the variable to check for
  integer, dimension(:), intent(inout) :: dim_sizes !< variable dimension sizes
  ! local
  logical :: file_open_success ! .true. if call to open_file is successful
  type(fmsNetcdfFile_t) :: fileobj ! netcdf file object returned by open_file
  integer :: i, ndims

  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, file_name, "read", is_restart=.false.)
  if (.not. (variable_exists(fileobj, trim(variable_name)))) &
    call MOM_error(FATAL, "MOM_io::field_size: variable "//trim(variable_name)// &
                   " not found in file "//trim(file_name))

  ndims = get_variable_num_dimensions(fileObj, trim(variable_name))
  if (size(dim_sizes) .ne. ndims) &
    call MOM_error(FATAL, "MOM_io::field_size: The number of dimensions of variable "//trim(variable_name)// &
    " does not match the size of the dim_sizes argument passed to the routine")
  call get_variable_size(fileObj, trim(variable_name), dim_sizes)

  if (check_if_open(fileobj)) call close_file(fileobj)
end subroutine field_size

!> set the logical variables that determine which diagnositic axes to use
subroutine get_horizontal_grid_logic(grid_string_id, use_lath, use_lonh, use_latq, use_lonq)
  character(len=*), intent(in) :: grid_string_id !< horizontal grid string
  logical, intent(out) :: use_lath !< if .true., y-axis is oriented in CENTER position
  logical, intent(out) :: use_lonh !< if .true., x-axis is oriented in CENTER position
  logical, intent(out) :: use_latq !< if .true., y-axis is oriented in NORTH_FACE position
  logical, intent(out) :: use_lonq !< if .true., x-axis is oriented in EAST_FACE position

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

!> Read the data associated with a named axis in a file
!>\todo: this routine is redundant; replace calls to this routine with MOM_read_data
subroutine read_axis_data(filename, axis_name, var)
  character(len=*),   intent(in)  :: filename  !< Name of the file to read
  character(len=*),   intent(in)  :: axis_name !< Name of the axis to read
  real, dimension(:), intent(out) :: var       !< The axis location data

  ! local
  logical :: file_open_success ! .true. if call to open_file is successful
  type(fmsNetcdfFile_t) :: fileobj ! netcdf file object returned by open_file
  integer :: i, ndims

  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, filename, "read", is_restart=.false.)
  if (.not. (variable_exists(fileobj, trim(axis_name)))) &
    call MOM_error(FATAL, "MOM_io::read_axis_data: variable "//trim(axis_name)// &
                   " not found in file "//trim(filename))

  ndims = get_variable_num_dimensions(fileObj, trim(axis_name))
  if (size(var) .ne. ndims) &
    call MOM_error(FATAL, "MOM_io::read_axis_data: The number of dimensions of variable "//trim(axis_name)// &
    " does not match the size of the var argument passed to the routine")
  call read_data(fileObj, trim(axis_name), var)

  if (check_if_open(fileobj)) call close_file(fileobj)

end subroutine read_axis_data

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
  logical :: file_open_success ! .true. if call to open_file is successful
  logical :: variableExists ! .true. if variable is found in file
  character(len=200) :: msg
  character(len=40), allocatable, dimension(:) :: dim_names ! variable dimension names
  integer :: i, ndims
  type(fmsNetcdfFile_t) :: fileobj !netcdf file object returned by open_file

  n_time = -1

  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, filename, "read", is_restart=.false.)

  ! check that variable is in the file
  if (.not.(variable_exists(fileobj, trim(varname)))) call MOM_error(FATAL, "num_time_levels: variable"//&
    trim(varname)//" not found in "//trim(filename))
  ! get the number of variable dimensions
  ndims = get_variable_num_dimensions(fileobj, trim(varname))

  if (present(min_dims)) then
    if (ndims .LT. min_dims-1) then
      write(msg, '(I3)') min_dims
      call MOM_error(WARNING, "num_time_levels: variable "//trim(varname)//&
        " in file "//trim(filename)//" has fewer than min_dims = "//trim(msg)//&
        " dimensions.")
    elseif (ndims .EQ. min_dims - 1) then
      n_time = 0 ; return
    endif
  endif
  ! check for the unlimited dimension and set n_time to the length of the unlimited dimension
  allocate(dim_names(ndims))

  call get_variable_dimension_names(fileobj, trim(varname), dim_names)
  do i=1,ndims
    if (is_dimension_unlimited(fileobj, trim(dim_names(i)))) &
      call get_dimension_size(fileobj, trim(dim_names(i)), n_time)
  enddo

  deallocate(dim_names)

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
                              time_level, stagger, scale)
  character(len=*),       intent(in)    :: filename !< name of the netcdf file to read
  character(len=*),       intent(in)    :: u_fieldname !< The variable name of the u data in the file
  character(len=*),       intent(in)    :: v_fieldname !< The variable name of the v data in the file
  real, dimension(:,:),   intent(inout) :: u_data    !< The 2 dimensional array into which the
                                                     !! u-component of the data should be read
  real, dimension(:,:),   intent(inout) :: v_data    !< The 2 dimensional array into which the
                                                     !! v-component of the data should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: time_level !< The time level in the file to read
  integer,      optional, intent(in)    :: stagger   !< A flag indicating where this vector is discretized
  real,         optional, intent(in)    :: scale     !< A scaling factor that the fields are multiplied
                                                     !! by before they are returned.
  ! local
  type(FmsNetcdfDomainFile_t) :: fileobj !< netcdf file object returned by call to open_file
  integer :: is, ie, js, je, i
  integer :: u_pos, v_pos
  integer, dimension(2) :: start, dim_sizes_u, dim_sizes_v
  character(len=32), dimension(2) :: dim_names_u, dim_names_v, units_u, units_v
  logical :: file_open_success ! .true. if open file is successful

  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, filename, "read", MOM_domain%mpp_domain, is_restart=.false.)
  if (.not. file_open_success) call MOM_error(FATAL, "MOM_read_vector_2d: netcdf file "//trim(filename)//" not opened.")

  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE .or. stagger == BGRID_NE ) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif

  start(:) = 1
  call get_variable_size(fileobj, u_fieldname, dim_sizes_u, broadcast=.true.)
  call get_variable_size(fileobj, v_fieldname, dim_sizes_v, broadcast=.true.)
  call get_variable_dimension_names(fileobj, u_fieldname, dim_names_u, broadcast=.true.)
  call get_variable_dimension_names(fileobj, v_fieldname, dim_names_v, broadcast=.true.)

  do i=1,2
    ! register the u axes
    call get_variable_units(fileobj, dim_names_u(i), units_u(i))
    if (trim(lowercase(units_u(i))) .eq. "degrees_east") then
      call register_axis(fileobj, dim_names_u(i), "x", domain_position=u_pos)
    elseif (trim(lowercase(units_u(i))) .eq. "degrees_north") then
      call register_axis(fileobj, dim_names_u(i), "y", domain_position=u_pos)
    else
      if (is_dimension_unlimited(fileobj, dim_names_u(i))) then
        if (present(time_level)) then
          start(i)=time_level
          dim_sizes_u(i) = 1
          dim_sizes_v(i) = 1
        endif
        call register_axis(fileobj, dim_names_u(i),dim_sizes_u(i))
      endif
    endif
    ! Register the v axes if they differ from the u axes
    if (trim(lowercase(dim_names_v(i))) .ne. trim(lowercase(dim_names_u(i)))) then
      call get_variable_units(fileobj, dim_names_v(i), units_v(i))
      if (trim(lowercase(units_v(i))) .eq. "degrees_east") then
        call register_axis(fileobj, dim_names_v(i), "x", domain_position=v_pos)
      elseif (trim(lowercase(units_v(i))) .eq. "degrees_north") then
        call register_axis(fileobj, dim_names_v(i), "y", domain_position=v_pos)
      endif
    endif
  enddo
  ! read the data
  call read_data(fileobj,u_fieldname, u_data, corner=start, edge_lengths=dim_sizes_u)
  call read_data(fileobj,v_fieldname, v_data, corner=start, edge_lengths=dim_sizes_v)
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)
  ! scale the data
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
                              time_level, stagger, scale)
  character(len=*),       intent(in)    :: filename !< name of the netcdf file to read
  character(len=*),       intent(in)    :: u_fieldname !< The variable name of the u data in the file
  character(len=*),       intent(in)    :: v_fieldname !< The variable name of the v data in the file
  real, dimension(:,:,:), intent(inout) :: u_data    !< The 3 dimensional array into which the
                                                     !! u-component of the data should be read
  real, dimension(:,:,:), intent(inout) :: v_data    !< The 3 dimensional array into which the
                                                     !! v-component of the data should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: time_level !< The time level in the file to read
  integer,      optional, intent(in)    :: stagger   !< A flag indicating where this vector is discretized
  real,         optional, intent(in)    :: scale     !< A scaling factor that the fields are multiplied
                                                     !! by before they are returned.
  ! local
  type(FmsNetcdfDomainFile_t) :: fileobj !< netcdf file object returned by call to open_file
  integer :: is, ie, js, je, i
  integer :: u_pos, v_pos
  integer, dimension(3) :: start, dim_sizes_u, dim_sizes_v
  character(len=32), dimension(3) :: dim_names_u, dim_names_v, units_u, units_v
  logical :: file_open_success ! .true. if open file is successful
  ! open the file
  if (.not.(check_if_open(fileobj))) &
    file_open_success = open_file(fileobj, filename, "read", MOM_domain%mpp_domain, is_restart=.false.)
  if (.not. file_open_success) call MOM_error(FATAL, "MOM_read_vector_3d: netcdf file "//trim(filename)//" not opened.")

  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == BGRID_NE) then ; u_pos = CORNER ; v_pos = CORNER
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif

  start(:) = 1
  call get_variable_size(fileobj, u_fieldname, dim_sizes_u, broadcast=.true.)
  call get_variable_size(fileobj, v_fieldname, dim_sizes_v, broadcast=.true.)
  call get_variable_dimension_names(fileobj, u_fieldname, dim_names_u, broadcast=.true.)
  call get_variable_dimension_names(fileobj, v_fieldname, dim_names_v, broadcast=.true.)

  do i=1,3
    ! register the u axes
    call get_variable_units(fileobj, dim_names_u(i), units_u(i))
    if (trim(lowercase(units_u(i))) .eq. "degrees_east") then
      call register_axis(fileobj, dim_names_u(i), "x", domain_position=u_pos)
    elseif (trim(lowercase(units_u(i))) .eq. "degrees_north") then
      call register_axis(fileobj, dim_names_u(i), "y", domain_position=u_pos)
    else
      if (is_dimension_unlimited(fileobj, dim_names_u(i))) then
        if (present(time_level)) then
          start(i)=time_level
          dim_sizes_u(i) = 1
          dim_sizes_v(i) = 1
        endif
        call register_axis(fileobj, dim_names_u(i), dim_sizes_u(i))
      endif
    endif
  ! register the v axes if the differ from the u axes
    if (trim(lowercase(dim_names_v(i))) .ne. trim(lowercase(dim_names_u(i)))) then
      call get_variable_units(fileobj, dim_names_v(i), units_v(i))
      if (trim(lowercase(units_v(i))) .eq. "degrees_east") then
        call register_axis(fileobj, dim_names_v(i), "x", domain_position=v_pos)
      elseif (trim(lowercase(units_v(i))) .eq. "degrees_north") then
        call register_axis(fileobj, dim_names_v(i), "y", domain_position=v_pos)
      endif
    endif
  enddo
  ! read the data
  call read_data(fileobj,u_fieldname, u_data, corner=start, edge_lengths=dim_sizes_u)
  call read_data(fileobj,v_fieldname, v_data, corner=start, edge_lengths=dim_sizes_v)
  ! close the file
  if (check_if_open(fileobj)) call close_file(fileobj)
  ! scale the data
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

!> convert the variable checksum integer(s) to a single string
!! If there is more than 1 checksum, commas are inserted between
!! each checksum value in the output string
function convert_checksum_to_string(checksum_int) result (checksum_string)
  integer(kind=8), intent(in) :: checksum_int !< checksum integer values
! local
  character(len=64) :: checksum_string
  integer :: i

  checksum_string = ''

  write (checksum_string,'(Z16)') checksum_int ! Z16 is the hexadecimal format code

end function convert_checksum_to_string

!> convert the variable checksum string read in from a restart file
!! to integer value(s)
function convert_checksum_string_to_int(checksum_char) result(checksum_file)
  character(len=*), intent(in) :: checksum_char !< checksum character array
  ! local
  integer(kind=8),dimension(3) :: checksum_file !< checksum string corresponds to
                                                   !< values from up to 3 times
  integer :: last
  integer :: start
  integer(kind=8) :: checksumh
  integer :: num_checksumh
  integer :: k

  start =0
  last = 0
  checksumh = 0
  num_checksumh = 1
  last = len_trim(checksum_char)
  start = index(trim(checksum_char),",") ! A start value of 0 implies only 1 checksum value
  ! Scan checksum character array for the ',' delimiter, which indicates that the corresponding variable
  ! has multiple time levels.
  do while ((start > 0) .and. (start < (last-15)))
     start = start + scan(checksum_char(start:last), "," ) ! move starting pointer after ","
     num_checksumh = num_checksumh + 1
  enddo

  start = 1

  do k = 1, num_checksumh
     read(checksum_char(start:start+15),'(Z16)') checksumh ! Z=hexadecimal integer: Z16 is for 64-bit data types
     checksum_file(k) = checksumh
     start = start+ 17 ! Move start index past the ',' in checksum_char
  enddo

end function convert_checksum_string_to_int

!> function to obtain the most recent time in a domain-decomposed netCDF file
function read_most_recent_time_DD(fileObj) result (file_time)
  type(fmsNetcdfDomainFile_t) fileObj ! netCDF file object returned by open_file
  ! local
  real :: file_time ! most recent time read from netcdf file. If there is no time value, time = 0
  integer :: dim_unlim_size
  character(len=40) :: dim_unlim_name

  if (.not. check_if_open(fileObj)) &
   call MOM_error(FATAL, "MOM_io:read_most_recent_time_DD : fileobj must be opened before calling this function")

  call get_unlimited_dimension_name(fileObj,dim_unlim_name)
  call get_dimension_size(fileObj, trim(dim_unlim_name), dim_unlim_size)
  if (dim_unlim_size .le. 0) then
    file_time = 0
  else
    call read_data(fileobj,trim(dim_unlim_name), file_time, unlim_dim_level=dim_unlim_size)
  endif
end function read_most_recent_time_DD

!> function to obtain the most recent time in a non-domain-decomposed netCDF file
function read_most_recent_time_noDD(fileObj) result (file_time)
  type(fmsNetcdfFile_t) fileObj ! netCDF file object returned by open_file
  ! local
  real :: file_time ! most recent time read from netcdf file. If there is no time value, file_time = 0
  integer :: dim_unlim_size
  character(len=40) :: dim_unlim_name

  if (.not. check_if_open(fileObj)) &
   call MOM_error(FATAL, "MOM_io:read_most_recent_time_DD : fileobj must be opened before calling this function")
  call get_unlimited_dimension_name(fileObj,dim_unlim_name)
  call get_dimension_size(fileObj, trim(dim_unlim_name), dim_unlim_size)
  if (dim_unlim_size .le. 0) then
    file_time = 0
  else
    call read_data(fileobj,trim(dim_unlim_name), file_time, unlim_dim_level=dim_unlim_size)
  endif
end function read_most_recent_time_noDD

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
!!  This file contains a number of subroutines that manipulate
!!  NetCDF files and handle input and output of fields.  These
!!  subroutines, along with their purpose, are:
!!
!!   * create_file: create a netCDF file for and register the global axes and variables to the file
!!   * write_field: write a field to an open file generated by a previous call to create_file
!!   * MOM_read_data: read a field from a netcdf file and apply a scaling factor if specified.
!!   * MOM_read_vector : read in the components (u,v) of a vector field and apply a scaling factor to the data
!!    if specified
!!   *

end module MOM_io
