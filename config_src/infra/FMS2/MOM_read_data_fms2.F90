!> This module contains routines that wrap the fms2 read_data calls
module MOM_read_data_fms2

! This file is part of MOM6. See LICENSE.md for the license.
use MOM_axis,             only : MOM_register_variable_axes
use MOM_error_infra,      only : MOM_error=>MOM_err, NOTE, FATAL, WARNING
use MOM_domain_infra,     only : MOM_domain_type, AGRID, BGRID_NE, CGRID_NE
use MOM_domain_infra,     only : domain2d, CENTER, CORNER, NORTH_FACE, EAST_FACE
use MOM_domain_infra,     only : rescale_comp_data
use MOM_string_functions, only : lowercase
use fms2_io_mod,          only : read_data, attribute_exists => variable_att_exists, variable_exists
use fms2_io_mod,          only : register_field, register_variable_attribute, fms2_open_file => open_file
use fms2_io_mod,          only : fms2_close_file => close_file, write_data, get_variable_dimension_names
use fms2_io_mod,          only : check_if_open, get_dimension_names, get_dimension_size
use fms2_io_mod,          only : is_dimension_registered, register_axis, get_variable_size
use fms2_io_mod,          only : FmsNetcdfDomainFile_t, FmsNetcdfFile_t, unlimited, get_variable_names
use fms2_io_mod,          only : get_variable_num_dimensions, get_variable_units, is_dimension_unlimited
use fms2_io_mod,          only : get_num_variables

implicit none ; private

public MOM_read_data_scalar, MOM_read_vector_2d_fms2, MOM_read_vector_3d_fms2
public MOM_read_data_4d_noDD, MOM_read_data_3d_noDD, MOM_read_data_2d_noDD, MOM_read_data_1d_noDD
public MOM_read_data_4d_DD, MOM_read_data_3d_DD, MOM_read_data_2d_DD, MOM_read_data_1d_DD

! CAUTION: The following variables are saved by default, and are only necessary for consecutive calls to
! MOM_read_data with the same file name. The user should ensure that fms2_close_file on
! the fileobj_read structures are called at every requisite time step at after the last
! variable is written to the file by omitting the optional leave_file_open argument, or setting it to .false.

!> netCDF domain-decomposed file object returned by call to
!! open_file in MOM_read_data_DD calls
type(FmsNetcdfDomainFile_t), private :: fileobj_read_dd

!> netCDF domain-decomposed file object returned by call to
!! open_file in MOM_read_data_noDD calls
type(FmsNetcdfFile_t), private :: fileobj_read

!> Type with variable metadata for a netCDF file opened to read domain-decomposed data
type file_variable_meta_DD
  integer :: nvars = 0!< number of variables in a netCDF file opened to read domain-decomposed data
  character(len=96), allocatable, dimension(:) :: var_names !< array for names of variables in a netCDF
                                                            !! file opened to read domain-decomposed data
end type file_variable_meta_DD

!> Type with variable metadata for a netCDF file opened to read non-domain-decomposed data
type file_variable_meta_noDD
  integer :: nvars = 0 !< number of variables in a netCDF file opened to read non-domain-decomposed data
  character(len=96), allocatable, dimension(:) :: var_names !< array for names of variables in a netCDF
                                                            !! file opened to read non-domain-decomposed data
end type file_variable_meta_noDD
!> type to hold metadata for variables in a domain-decomposed file
type (file_variable_meta_DD), private :: file_var_meta_DD

!> type to hold metadata for variables in a non-domain-decomposed file
type (file_variable_meta_noDD), private :: file_var_meta_noDD

! !> index of the time_level value that is written to netCDF file by the write_field routines.
! integer, private :: write_field_time_index

contains

!> This routine calls the fms_io read_data subroutine to read 1-D domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call.
subroutine MOM_read_data_1d_DD(filename, fieldname, data, domain, start_index, edge_lengths, &
                               timelevel, scale, x_position, y_position, leave_file_open)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:), intent(inout) :: data !< The 1-dimensional data array to pass to read_data
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  integer, dimension(1), optional, intent(in) :: start_index !< starting index of data buffer. Default is 1
  integer, dimension(1), optional, intent(in) :: edge_lengths !< number of data values to read in
                                                              !! default is the variable size
  integer, optional, intent(in) :: timelevel !< time level to read
  real, optional, intent(in) :: scale !< A scaling factor that the field is multiplied by
  integer, optional, intent(in) :: x_position !< domain position of x-dimension; CENTER (default) or EAST_FACE
  integer, optional, intent(in) :: y_position !< domain position of y-dimension; CENTER (default) or NORTH_FACE
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: variable_found ! .true. if lowercase(fieldname) matches one of the lowercase file variable names
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  integer :: i, num_var_dims, dim_unlim_size
  integer, dimension(1) :: start, nread ! indices for first data value and number of values to read
  character(len=40), allocatable :: dim_names(:) ! variable dimension names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file
  integer :: xpos, ypos ! x and y domain positions

  xpos = CENTER
  ypos = CENTER
  if (present(x_position)) xpos = x_position
  if (present(y_position)) ypos = y_position

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  ! open the file
  if (.not.(check_if_open(fileobj_read_dd))) then
    file_open_success = fms2_open_file(fileobj_read_dd, filename, "read", domain%mpp_domain, is_restart=.false.)
    file_var_meta_DD%nvars = get_num_variables(fileobj_read_dd)
    if (file_var_meta_DD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_DD%var_names))) allocate(file_var_meta_DD%var_names(file_var_meta_DD%nvars))
    call get_variable_names(fileobj_read_dd, file_var_meta_DD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_DD%nvars
    if (lowercase(trim(file_var_meta_DD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_DD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) call MOM_error(FATAL, "MOM_read_data_fms2:MOM_read_data_1d_DD: "//&
    trim(fieldname)//" not found in"//trim(filename))
  ! register the variable axes
  call MOM_register_variable_axes(fileobj_read_dd, trim(variable_to_read), xPosition=xpos, yPosition=ypos)
  ! set the start and nread values that will be passed as the read_data start_index and edge_lengths arguments
  allocate(dim_names(num_var_dims))
  dim_names(:) = ""
  call get_variable_dimension_names(fileobj_read_dd, trim(variable_to_read), dim_names)

  start(1)=1
  if (present(timelevel)) then
    if (is_dimension_unlimited(fileobj_read_dd, dim_names(1))) start(1) = timelevel
  elseif (present(start_index)) then
    start(1) = start_index(1)
  endif

  if (present(edge_lengths)) then
    nread(1) = edge_lengths(1)
  else
    call get_dimension_size(fileobj_read_dd, trim(dim_names(1)), nread(1))
  endif
  ! read the data
  dim_unlim_size = 0
  if (present(timelevel)) then
    do i=1,num_var_dims
      if (is_dimension_unlimited(fileobj_read_dd, dim_names(i))) then
        call get_dimension_size(fileobj_read_dd, dim_names(i), dim_unlim_size)
        exit
      endif
    enddo
    if (dim_unlim_size .gt. 0) then
      call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread, &
                     unlim_dim_level=timelevel)
    else
      call MOM_error(WARNING, "MOM_read_data_fms2::MOM_read_data_1d_DD: time level specified, but the variable "//&
                      trim(fieldName)// " does not have an unlimited dimension.")
      call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread)
    endif
  else
    call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread)
  endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read_dd)) call fms2_close_file(fileobj_read_dd)
     if (allocated(file_var_meta_DD%var_names)) deallocate(file_var_meta_DD%var_names)
     file_var_meta_DD%nvars = 0
  endif
  if (allocated(dim_names)) deallocate(dim_names)

  ! Rescale the data that was read if necessary.
  if (present(scale)) then ; if (scale /= 1.0) then
    data(:) = scale*data(:)
  endif ; endif

end subroutine MOM_read_data_1d_DD

!> This routine calls the fms_io read_data subroutine to read 2-D domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call.
subroutine MOM_read_data_2d_DD(filename, fieldname, data, domain,start_index, edge_lengths, &
                               timelevel, scale, x_position, y_position, leave_file_open)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:), intent(inout) :: data !< The 2-dimensional data array to pass to read_data
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  integer, dimension(2), optional, intent(in) :: start_index !< starting indices of data buffer. Default is 1
  integer, dimension(2), optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timelevel !< time level to read
  real, optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  integer, intent(in), optional :: x_position !< domain position of x-dimension; CENTER (default) or EAST_FACE
  integer, intent(in), optional :: y_position !< domain position of y-dimension; CENTER (default) or NORTH_FACE
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: variable_found ! .true. if lowercase(fieldname) matches one of the lowercase file variable names
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  integer :: i, dim_unlim_size, num_var_dims, first(2), last(2)
  integer :: start(2), nread(2) ! indices for first data value and number of values to read
  character(len=40), allocatable :: dim_names(:) ! variable dimension names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file
  integer :: xpos, ypos, pos ! x and y domain positions
  integer :: isc, iec, jsc, jec, isg, ieg, jsg, jeg

  xpos = CENTER
  ypos = CENTER
  if (present(x_position)) xpos = x_position
  if (present(y_position)) ypos = y_position

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)

  ! open the file
  if (.not.(check_if_open(fileobj_read_dd))) then
    file_open_success = fms2_open_file(fileobj_read_dd, filename, "read", domain%mpp_domain, is_restart=.false.)
    file_var_meta_DD%nvars = get_num_variables(fileobj_read_dd)
    if (file_var_meta_DD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_DD%var_names))) allocate(file_var_meta_DD%var_names(file_var_meta_DD%nvars))
    call get_variable_names(fileobj_read_dd, file_var_meta_DD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_DD%nvars
    if (lowercase(trim(file_var_meta_DD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_DD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) call MOM_error(FATAL, "MOM_read_data_fms2:MOM_read_data_2d_DD: "//&
    trim(fieldname)//" not found in "//trim(filename))
  ! register the variable axes
  call MOM_register_variable_axes(fileobj_read_dd, trim(variable_to_read), xPosition=xpos, yPosition=ypos)

  pos = CENTER
  if (present(x_position)) then
    if (present(y_position)) then
      pos = CORNER
    else
      pos = xpos
    endif
  elseif (present(y_position)) then
    pos = ypos
  endif
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths argument
  num_var_dims = get_variable_num_dimensions(fileobj_read_dd, trim(variable_to_read))
  allocate(dim_names(num_var_dims))
  dim_names(:) = ""
  call get_variable_dimension_names(fileobj_read_dd, trim(variable_to_read), dim_names)
  !last(1) = iec - isg + 1 ! get array indices for the axis data
  !last(2) = jec - jsg + 1
  !first(1) = isc - isg + 1
  !first(2) = jsc - jsg + 1

  start(:) = 1
  if (present(start_index)) then
    start = start_index
  !else
  !  start(:) = first(:)
  endif

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    nread = shape(data)
  endif
  ! read the data
  dim_unlim_size=0
  if (present(timelevel)) then
    do i=1,num_var_dims
      if (is_dimension_unlimited(fileobj_read_dd, dim_names(i))) then
        call get_dimension_size(fileobj_read_dd, dim_names(i), dim_unlim_size)
      endif
    enddo
    if (dim_unlim_size .gt. 0) then
      call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread, &
                     unlim_dim_level=timelevel)
    else
      call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread)
    endif
  else
    call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread)
  endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read_dd)) call fms2_close_file(fileobj_read_dd)
    if (allocated(file_var_meta_DD%var_names)) deallocate(file_var_meta_DD%var_names)
    file_var_meta_DD%nvars = 0
  endif
  if (allocated(dim_names)) deallocate(dim_names)

  ! Rescale the data that was read if necessary.
  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(domain, data, scale)
  endif ; endif

end subroutine MOM_read_data_2d_DD

!> This routine calls the fms_io read_data subroutine to read 3-D domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call.
subroutine MOM_read_data_3d_DD(filename, fieldname, data, domain, start_index, edge_lengths, &
                               timelevel, scale, x_position, y_position, leave_file_open)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:), intent(inout) :: data !< The 3-dimensional data array to pass to read_data
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  integer, dimension(3), optional, intent(in) :: start_index !< starting indices of data buffer. Default is 1
  integer, dimension(3), optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                             !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timelevel !< time level to read
  real,         optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  integer, intent(in), optional :: x_position !< domain position of x-dimension; CENTER (default) or EAST_FACE
  integer, intent(in), optional :: y_position !< domain position of y-dimension; CENTER (default) or NORTH_FACE
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: variable_found ! if .true., the variable was found in the netCDF file
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  integer :: i, dim_unlim_size, num_var_dims
  integer, dimension(3) :: start, nread, first, last ! indices for first data value and number of values to read
  character(len=40), allocatable :: dim_names(:) ! variable dimension names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file
  integer :: xpos, ypos, pos ! x and y domain positions
  integer :: isc, iec, jsc, jec, isg, ieg, jsg, jeg

  xpos = CENTER
  ypos = CENTER
  if (present(x_position)) xpos = x_position
  if (present(y_position)) ypos = y_position

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  ! open the file
  if (.not.(check_if_open(fileobj_read_dd))) then
    file_open_success = fms2_open_file(fileobj_read_dd, filename, "read", domain%mpp_domain, is_restart=.false.)
    file_var_meta_DD%nvars = get_num_variables(fileobj_read_dd)
    if (file_var_meta_DD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_DD%var_names))) allocate(file_var_meta_DD%var_names(file_var_meta_DD%nvars))
    call get_variable_names(fileobj_read_dd, file_var_meta_DD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_DD%nvars
    if (lowercase(trim(file_var_meta_DD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_DD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) call MOM_error(FATAL, "MOM_read_data_fms2:MOM_read_data_3d_DD: "//&
    trim(fieldname)//" not found in"//trim(filename))
  ! register the variable axes
  call MOM_register_variable_axes(fileobj_read_dd, trim(variable_to_read), xPosition=xpos, yPosition=ypos)
  pos = CENTER
  if (present(x_position)) then
    if (present(y_position)) then
      pos = CORNER
    else
      pos = xpos
    endif
  elseif (present(y_position)) then
    pos = ypos
  endif
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths argument
  num_var_dims = get_variable_num_dimensions(fileobj_read_dd, trim(variable_to_read))
  allocate(dim_names(num_var_dims))
  dim_names(:) = ""
  call get_variable_dimension_names(fileobj_read_dd, trim(variable_to_read), dim_names)

  start(:) = 1
  if (present(start_index)) then
    start = start_index
  !else
  !  start(1:2) = first(1:2)
  endif

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    !nread(1) = last(1) - first(1) + 1
    !nread(2) = last(2) - first(2) + 1
    nread = shape(data)
  endif
 ! read the data
  dim_unlim_size=0
  if (present(timelevel)) then
    do i=1,num_var_dims
      if (is_dimension_unlimited(fileobj_read_dd, dim_names(i))) then
        call get_dimension_size(fileobj_read_dd, dim_names(i), dim_unlim_size)
      endif
    enddo
    if (dim_unlim_size .gt. 0) then
      call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread, &
                     unlim_dim_level=timelevel)
    else
      call MOM_error(WARNING, "MOM_read_data_fms2::MOM_read_data_3d_DD: time level specified, but the variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
      call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread)
    endif
  else
    call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread)
  endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read_dd)) call fms2_close_file(fileobj_read_dd)
     if (allocated(file_var_meta_DD%var_names)) deallocate(file_var_meta_DD%var_names)
    file_var_meta_DD%nvars = 0
  endif

  if (allocated(dim_names)) deallocate(dim_names)

  ! Rescale the data that was read if necessary.
  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(domain, data, scale)
  endif ; endif

end subroutine MOM_read_data_3d_DD

!> This routine calls the fms_io read_data subroutine to read 4-D domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call.
subroutine MOM_read_data_4d_DD(filename, fieldname, data, domain, start_index, edge_lengths, &
                               timelevel, scale, x_position, y_position, leave_file_open)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:,:), intent(inout) :: data !< The 1-dimensional data array to pass to read_data
  type(MOM_domain_type), intent(in) :: domain !< MOM domain attribute with the mpp_domain decomposition
  integer, dimension(4),  optional, intent(in) :: start_index !< starting indices of data buffer. Default is 1
  integer, dimension(4),  optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timelevel !< time level to read
  real,         optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  integer, intent(in), optional :: x_position !< domain position of x-dimension; CENTER (default) or EAST_FACE
  integer, intent(in), optional :: y_position !< domain position of y-dimension; CENTER (default) or NORTH_FACE
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: variable_found ! .true. if lowercase(fieldname) matches one of the lowercase file variable names
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  integer :: i, dim_unlim_size, num_var_dims
  integer, dimension(4) :: start, nread, first, last ! indices for first data value and number of values to read
  character(len=40), allocatable :: dim_names(:) ! variable dimension names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file
  integer :: xpos, ypos, pos ! x and y domain positions
  integer :: isc, iec, jsc, jec, isg, ieg, jsg, jeg

  xpos = CENTER
  ypos = CENTER
  if (present(x_position)) xpos = x_position
  if (present(y_position)) ypos = y_position

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  ! open the file
  if (.not.(check_if_open(fileobj_read_dd))) then
    file_open_success = fms2_open_file(fileobj_read_dd, filename, "read", domain%mpp_domain, is_restart=.false.)
    file_var_meta_DD%nvars = get_num_variables(fileobj_read_dd)
    if (file_var_meta_DD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_DD%var_names))) allocate(file_var_meta_DD%var_names(file_var_meta_DD%nvars))
    call get_variable_names(fileobj_read_dd, file_var_meta_DD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_DD%nvars
    if (lowercase(trim(file_var_meta_DD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_DD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) &
    call MOM_error(FATAL, "MOM_read_data_fms2:MOM_read_data_4d_DD: "//trim(fieldname)//" not found in"//&
                   trim(filename))
  ! register the variable axes
  call MOM_register_variable_axes(fileobj_read_dd, trim(variable_to_read), xPosition=xpos, yPosition=ypos)
  pos = CENTER
  if (present(x_position)) then
    if (present(y_position)) then
      pos = CORNER
    else
      pos = xpos
    endif
  elseif (present(y_position)) then
    pos = ypos
  endif
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths argument
  num_var_dims = get_variable_num_dimensions(fileobj_read_dd, trim(variable_to_read))
  allocate(dim_names(num_var_dims))
  dim_names(:) = ""
  call get_variable_dimension_names(fileobj_read_dd, trim(variable_to_read), dim_names)

  start(:) = 1
  if (present(start_index)) then
    start(:) = start_index(:)
  !else
    !start(1:2) = first(1:2)
  endif

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    !nread(1) = last(1) - first(1) + 1
    !nread(2) = last(2) - first(2) + 1
    nread = shape(data)
  endif
  ! read the data
  dim_unlim_size=0
  if (present(timelevel)) then
    do i=1, num_var_dims
      if (is_dimension_unlimited(fileobj_read_dd, dim_names(i))) then
        call get_dimension_size(fileobj_read_dd, dim_names(i), dim_unlim_size)
      endif
      if (i .eq. 4) then
        nread(i) = 1
        start(i) = timelevel
      endif
    enddo
    if (dim_unlim_size .gt. 0) then
      call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread, &
                     unlim_dim_level=timelevel)
    else
      call MOM_error(WARNING, "MOM_read_data_fms2::MOM_read_data_4d_DD: time level specified, but the variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
      call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread)
    endif
  else
    call read_data(fileobj_read_dd, trim(variable_to_read), data, corner=start, edge_lengths=nread)
  endif

  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read_dd)) call fms2_close_file(fileobj_read_dd)
     if (allocated(file_var_meta_DD%var_names)) deallocate(file_var_meta_DD%var_names)
    file_var_meta_DD%nvars = 0
  endif
  if (allocated(dim_names)) deallocate(dim_names)

  ! Rescale the data that was read if necessary.
  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(domain, data, scale)
  endif ; endif

end subroutine MOM_read_data_4d_DD

!!> This routine calls the fms_io read_data subroutine to read a scalar (0-D) field named "fieldname"
!! from file "filename".
subroutine MOM_read_data_scalar(filename, fieldname, data, scale, leave_file_open)
  character(len=*), intent(in) :: filename !< The name of the file to read
  character(len=*), intent(in) :: fieldname !< The variable name of the data in the file
  real, intent(inout) :: data !< data buffer to pass to read_data
  real,   optional, intent(in) :: scale     !< A scaling factor that the scalar is multiplied
                                            !! by before it is returned.
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  integer :: i
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  logical :: variable_found ! .true. if lowercase(fieldname) matches one of the lowercase file variable names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  ! open the file
  if (.not.(check_if_open(fileobj_read))) then
    file_open_success = fms2_open_file(fileobj_read, filename, "read", is_restart=.false.)
    file_var_meta_noDD%nvars = get_num_variables(fileobj_read)
    if (file_var_meta_noDD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_noDD%var_names))) &
      allocate(file_var_meta_noDD%var_names(file_var_meta_noDD%nvars))
    call get_variable_names(fileobj_read, file_var_meta_noDD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_noDD%nvars
    if (lowercase(trim(file_var_meta_noDD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_noDD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) call MOM_error(FATAL, "MOM_read_data_fms2:MOM_read_data_scalar: "//trim(fieldname)// &
                                            " not found in"//trim(filename))
  ! read the data
  call read_data(fileobj_read, trim(fieldname), data)
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read)) call fms2_close_file(fileobj_read)
    if (allocated(file_var_meta_noDD%var_names)) deallocate(file_var_meta_noDD%var_names)
    file_var_meta_noDD%nvars = 0
  endif

  ! Rescale the data that was read if necessary.
  if (present(scale)) then ; if (scale /= 1.0) then
    data = scale*data
  endif ; endif

end subroutine MOM_read_data_scalar

!> This routine calls the fms_io read_data subroutine to read 1-D non-domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call.
subroutine MOM_read_data_1d_noDD(filename, fieldname, data, start_index, &
                                 edge_lengths, timelevel, scale, leave_file_open)
  character(len=*),       intent(in) :: filename !< The name of the file to read
  character(len=*),       intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:),     intent(inout) :: data !< The 1-dimensional data array to pass to read_data
  integer, dimension(1), optional, intent(in) :: start_index !< starting index of data buffer. Default is 1
  integer, dimension(1), optional, intent(in) :: edge_lengths !< number of data values to read in; default is
                                                             !! the variable size
  integer,      optional, intent(in) :: timelevel !< time level to read
  real,         optional, intent(in) :: scale !< A scaling factor that the field is multiplied by
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  logical :: variable_found ! .true. if lowercase(fieldname) matches one of the lowercase file variable names
  integer :: i, num_var_dims, dim_unlim_size
  integer, dimension(1) :: start, nread ! indices for first data value and number of values to read
  character(len=40), allocatable:: dim_names(:) ! variable dimension names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  ! open the file
  if (.not.(check_if_open(fileobj_read))) then
    file_open_success = fms2_open_file(fileobj_read, filename, "read", is_restart=.false.)
    file_var_meta_noDD%nvars = get_num_variables(fileobj_read)
    if (file_var_meta_noDD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_noDD%var_names))) &
      allocate(file_var_meta_noDD%var_names(file_var_meta_noDD%nvars))
    call get_variable_names(fileobj_read, file_var_meta_noDD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_noDD%nvars
    if (lowercase(trim(file_var_meta_noDD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_noDD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) call MOM_error(FATAL, "MOM_io:MOM_read_data_1d_noDD: "//trim(fieldname)//&
                                            " not found in "//trim(filename))

  num_var_dims = get_variable_num_dimensions(fileobj_read, trim(fieldname))
  allocate(dim_names(num_var_dims))
  dim_names(:) = ""
  call get_variable_dimension_names(fileobj_read, trim(variable_to_read), dim_names)

  ! set the start and nread values that will be passed as the read_data start_index and edge_lengths arguments
  start(1)=1
  if (present(timelevel)) then
    if (is_dimension_unlimited(fileobj_read, dim_names(1))) start(1) = timelevel
  elseif (present(start_index)) then
    start(1) = start_index(1)
  endif

  if (present(edge_lengths)) then
    nread(1) = edge_lengths(1)
  else
    nread = shape(data)
  endif
  ! read the data
  dim_unlim_size = 0
  if (present(timelevel)) then
    do i=1,num_var_dims
      if (is_dimension_unlimited(fileobj_read, dim_names(i))) then
        call get_dimension_size(fileobj_read, dim_names(i), dim_unlim_size)
        exit
      endif
    enddo
    if (dim_unlim_size .gt. 0) then
      call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread, &
                     unlim_dim_level=timelevel)
    else
      call MOM_error(WARNING, "MOM_read_data_fms2::MOM_read_data_1d_noDD: time level specified, but the variable "//&
                      trim(fieldName)// " does not have an unlimited dimension.")
      call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread)
    endif
  else
    call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread)
  endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read)) call fms2_close_file(fileobj_read)
    if (allocated(file_var_meta_noDD%var_names)) deallocate(file_var_meta_noDD%var_names)
    file_var_meta_noDD%nvars = 0
  endif
  if (allocated(dim_names)) deallocate(dim_names)

  ! Rescale the data that was read if necessary.
  if (present(scale)) then ; if (scale /= 1.0) then
    data(:) = scale*data(:)
  endif ; endif

end subroutine MOM_read_data_1d_noDD

!> This routine calls the fms_io read_data subroutine to read 2-D non-domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call.
subroutine MOM_read_data_2d_noDD(filename, fieldname, data, start_index, &
                                 edge_lengths, timelevel, scale, leave_file_open)
  character(len=*),       intent(in) :: filename !< The name of the file to read
  character(len=*),       intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:),     intent(inout) :: data !< The 2-dimensional data array to pass to read_data
  integer, dimension(2),  optional, intent(in) :: start_index !< starting indices of data buffer. Default is 1
  integer, dimension(2),  optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timelevel !< time level to read
  real, optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: variable_found ! .true. if lowercase(fieldname) matches one of the lowercase file variable names
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  integer :: i, dim_unlim_size, num_var_dims
  integer, dimension(2) :: start, nread ! indices for first data value and number of values to read
  character(len=40), allocatable :: dim_names(:) ! variable dimension names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  ! open the file
  if (.not.(check_if_open(fileobj_read))) then
    file_open_success = fms2_open_file(fileobj_read, filename, "read", is_restart=.false.)
    file_var_meta_noDD%nvars = get_num_variables(fileobj_read)
    if (file_var_meta_noDD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_noDD%var_names))) &
      allocate(file_var_meta_noDD%var_names(file_var_meta_noDD%nvars))
    call get_variable_names(fileobj_read, file_var_meta_noDD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_noDD%nvars
    if (lowercase(trim(file_var_meta_noDD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_noDD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) &
    call MOM_error(FATAL, "MOM_io:MOM_read_data_2d_noDD: "//trim(fieldname)//&
                                            " not found in "//trim(filename))
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  start(:) = 1
  if (present(start_index)) start = start_index

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    nread = shape(data)
  endif
  ! read the data
  dim_unlim_size=0
  if (present(timelevel)) then
    num_var_dims = get_variable_num_dimensions(fileobj_read, trim(fieldname))
    allocate(dim_names(num_var_dims))
    call get_variable_dimension_names(fileobj_read, trim(variable_to_read), dim_names)
    dim_names(:) = ""
    do i=1,num_var_dims
      if (is_dimension_unlimited(fileobj_read, dim_names(i))) then
        call get_dimension_size(fileobj_read, dim_names(i), dim_unlim_size)
      endif
    enddo
    if (dim_unlim_size .LE. 0) then
      call MOM_error(WARNING, "MOM_io::MOM_read_data_2d_noDD: time level specified, but the variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
      call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread)
    else
      call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread, &
                     unlim_dim_level=timelevel)
    endif
  else
    call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread)
  endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read)) call fms2_close_file(fileobj_read)
    if (allocated(file_var_meta_noDD%var_names)) deallocate(file_var_meta_noDD%var_names)
    file_var_meta_noDD%nvars = 0
  endif
  if(allocated(dim_names)) deallocate(dim_names)

  ! Rescale the data that was read if necessary.
  if (present(scale)) then ; if (scale /= 1.0) then
    data(:,:) = scale*data(:,:)
  endif ; endif

end subroutine MOM_read_data_2d_noDD

!> This routine calls the fms_io read_data subroutine to read 3-D non-domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call.
subroutine MOM_read_data_3d_noDD(filename, fieldname, data, start_index, &
                                 edge_lengths, timelevel, scale, leave_file_open)
  character(len=*),       intent(in) :: filename !< The name of the file to read
  character(len=*),       intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:),     intent(inout) :: data !< The 3-dimensional data array to pass to read_data
  integer, dimension(3), optional, intent(in) :: start_index !< starting indices of data buffer. Default is 1
  integer, dimension(3), optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timelevel !< time level to read
  real, optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: variable_found ! .true. if lowercase(fieldname) matches one of the lowercase file variable names
  logical :: close_the_file ! indicates whether to close the file after write_data is called; default is .true.
  integer :: i, dim_unlim_size, num_var_dims
  integer, dimension(3) :: start, nread ! indices for first data value and number of values to read
  character(len=40), allocatable :: dim_names(:) ! variable dimension names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  ! open the file
  if (.not.(check_if_open(fileobj_read))) then
    file_open_success = fms2_open_file(fileobj_read, filename, "read", is_restart=.false.)
    file_var_meta_noDD%nvars = get_num_variables(fileobj_read)
    if (file_var_meta_noDD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_noDD%var_names))) &
      allocate(file_var_meta_noDD%var_names(file_var_meta_noDD%nvars))
    call get_variable_names(fileobj_read, file_var_meta_noDD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_noDD%nvars
    if (lowercase(trim(file_var_meta_noDD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_noDD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) call MOM_error(FATAL, "MOM_io:MOM_read_data_3d_noDD: "//trim(fieldname)//&
                                            " not found in "//trim(filename))
  ! get the variable dimensions
  num_var_dims = get_variable_num_dimensions(fileobj_read, trim(fieldname))
  allocate(dim_names(num_var_dims))
  dim_names(:) = ""
  call get_variable_dimension_names(fileobj_read, trim(variable_to_read), dim_names)
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  start(:) = 1
  if (present(start_index)) start = start_index

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    nread = shape(data)
  endif
  ! read the data
  dim_unlim_size=0
  if (present(timelevel)) then
    do i=1,num_var_dims
      if (is_dimension_unlimited(fileobj_read, dim_names(i))) then
        call get_dimension_size(fileobj_read, dim_names(i), dim_unlim_size)
      endif
    enddo
    if (dim_unlim_size .LE. 0) then
      call MOM_error(WARNING, "MOM_read_data_fms2::MOM_read_data_3d_noDD: time level specified, but the variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
      call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread)
    else
      call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread, &
                     unlim_dim_level=timelevel)
    endif
  else
    call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread)
  endif

  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read)) call fms2_close_file(fileobj_read)
    if (allocated(file_var_meta_noDD%var_names)) deallocate(file_var_meta_noDD%var_names)
    file_var_meta_noDD%nvars = 0
  endif
  if (allocated(dim_names)) deallocate(dim_names)

  ! Rescale the data that was read if necessary.
  if (present(scale)) then ; if (scale /= 1.0) then
    data(:,:,:) = scale*data(:,:,:)
  endif ; endif

end subroutine MOM_read_data_3d_noDD

!> This routine calls the fms_io read_data subroutine to read 4-D non-domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call.
subroutine MOM_read_data_4d_noDD(filename, fieldname, data, start_index, &
                                 edge_lengths, timelevel, scale, leave_file_open)
  character(len=*),       intent(in) :: filename !< The name of the file to read
  character(len=*),       intent(in) :: fieldname !< The variable name of the data in the file
  real, dimension(:,:,:,:),     intent(inout) :: data !< The 4-dimensional data array to pass to read_data
  integer, dimension(4),  optional, intent(in) :: start_index !< starting indices of data buffer. Default is 1
  integer, dimension(4),  optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer, optional, intent(in) :: timelevel !< time level to read
  real,         optional, intent(in):: scale !< A scaling factor that the field is multiplied by
  ! local
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: variable_found ! .true. if lowercase(fieldname) matches one of the lowercase file variable names
  logical :: close_the_file ! indicates whether to close the file after read_data is called; default is .true.
  integer :: i, dim_unlim_size, num_var_dims
  integer, dimension(4) :: start, nread ! indices for first data value and number of values to read
  character(len=40), allocatable :: dim_names(:) ! variable dimension names
  character(len=96) :: variable_to_read ! variable to read from the netcdf file
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  ! open the file
  if (.not.(check_if_open(fileobj_read))) then
    file_open_success = fms2_open_file(fileobj_read, filename, "read", is_restart=.false.)
    file_var_meta_noDD%nvars = get_num_variables(fileobj_read)
    if (file_var_meta_noDD%nvars .lt. 1) call MOM_error(FATAL, "nvars is less than 1 for file "// &
                                                           trim(filename))
    if (.not.(allocated(file_var_meta_noDD%var_names))) &
      allocate(file_var_meta_noDD%var_names(file_var_meta_noDD%nvars))
    call get_variable_names(fileobj_read, file_var_meta_noDD%var_names)
  endif
  ! search for the variable in the file
  variable_to_read = ""
  variable_found = .false.
  do i=1,file_var_meta_noDD%nvars
    if (lowercase(trim(file_var_meta_noDD%var_names(i))) .eq. lowercase(trim(fieldname))) then
      variable_found = .true.
      variable_to_read = trim(file_var_meta_noDD%var_names(i))
      exit
    endif
  enddo
  if (.not.(variable_found)) call MOM_error(FATAL, "MOM_read_data_fms2:MOM_read_data_4d_noDD: "//&
                                            trim(fieldname)//" not found in "//trim(filename))
  ! get the variable dimensions
  num_var_dims = get_variable_num_dimensions(fileobj_read, trim(fieldname))
  allocate(dim_names(num_var_dims))
  dim_names(:) = ""
  call get_variable_dimension_names(fileobj_read, trim(variable_to_read), dim_names)
  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  start(:) = 1
  if (present(start_index)) start = start_index

  if (present(edge_lengths)) then
    nread = edge_lengths
  else
    nread = shape(data)
  endif
  ! read the data
  dim_unlim_size=0
  if (present(timelevel)) then
    do i=1, num_var_dims
      if (is_dimension_unlimited(fileobj_read, dim_names(i))) then
        call get_dimension_size(fileobj_read, dim_names(i), dim_unlim_size)
      endif
      if (i .eq. 4) then
        nread(i) = 1
        start(i) = timelevel
      endif
    enddo
    if (dim_unlim_size .LE. 0) then
      call MOM_error(WARNING, "MOM_io::MOM_read_data_4d_noDD: time level specified, but the variable "//&
                     trim(fieldName)// " does not have an unlimited dimension.")
      call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread)
    else
      call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread, &
                     unlim_dim_level=timelevel)
    endif
  else
    call read_data(fileobj_read, trim(variable_to_read), data, corner=start, edge_lengths=nread)
  endif
 ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read)) call fms2_close_file(fileobj_read)
    if (allocated(file_var_meta_noDD%var_names)) deallocate(file_var_meta_noDD%var_names)
    file_var_meta_noDD%nvars = 0
  endif
  if (allocated(dim_names)) deallocate(dim_names)

  ! Rescale the data that was read if necessary.
  if (present(scale)) then ; if (scale /= 1.0) then
    data(:,:,:,:) = scale*data(:,:,:,:)
  endif ; endif

end subroutine MOM_read_data_4d_noDD


!> This routine uses the fms2_io read_data interface to read a pair of distributed
!! 2-D data fields with names given by "[uv]_fieldname" from file "filename".  Valid values for
!! "stagger" include CGRID_NE, BGRID_NE, and AGRID.
subroutine MOM_read_vector_2d_fms2(filename, u_fieldname, v_fieldname, u_data, v_data, MOM_Domain, &
                                   timelevel, stagger, scale, leave_file_open)
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
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  integer :: is, ie, js, je, i, ndims, dim_unlim_index
  integer :: u_pos, v_pos
  integer, allocatable :: dim_sizes_u(:), dim_sizes_v(:)
  character(len=32), allocatable :: dim_names_u(:), dim_names_v(:), units_u(:), units_v(:)
  character(len=1) :: x_or_y ! orientation of cartesian coordinate axis
  logical :: is_valid
  logical :: file_open_success ! .true. if open file is successful
  logical :: close_the_file ! indicates whether to close the file after MOM_read_vector is called; default is .true.

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)

  ! open the file
  if (.not.(check_if_open(fileobj_read_dd))) &
    file_open_success = fms2_open_file(fileobj_read_dd, filename, "read", MOM_domain%mpp_domain, is_restart=.false.)
  if (.not. file_open_success) call MOM_error(FATAL, "MOM_read_vector_2d_fms2: netcdf file "//&
                                              trim(filename)//" not opened.")

  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE .or. stagger == BGRID_NE ) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif

  ndims = get_variable_num_dimensions(fileobj_read_dd, u_fieldname)
  allocate(dim_sizes_u(ndims))
  allocate(dim_sizes_v(ndims))
  allocate(dim_names_u(ndims))
  allocate(dim_names_v(ndims))
  allocate(units_u(ndims))
  allocate(units_v(ndims))

  units_u(:) = ""
  units_v(:) = ""
  dim_names_u(:) = ""
  dim_names_v(:) = ""
  dim_sizes_u(:) = 0
  dim_sizes_v(:) = 0

  call get_variable_size(fileobj_read_dd, u_fieldname, dim_sizes_u)
  call get_variable_size(fileobj_read_dd, v_fieldname, dim_sizes_v)
  call get_variable_dimension_names(fileobj_read_dd, u_fieldname, dim_names_u)
  call get_variable_dimension_names(fileobj_read_dd, v_fieldname, dim_names_v)

  do i=1,ndims
    ! register the u axes
    if (.not.(is_dimension_registered(fileobj_read_dd, dim_names_u(i)))) then
      call get_variable_units(fileobj_read_dd, dim_names_u(i), units_u(i))
      call validate_lat_lon_units(units_u(i), x_or_y, is_valid)
      if (is_valid) then
        call register_axis(fileobj_read_dd, dim_names_u(i), x_or_y, domain_position=u_pos)
      else
        call register_axis(fileobj_read_dd, dim_names_u(i), dim_sizes_u(i))
      endif
    endif
    ! Register the v axes if they differ from the u axes
    if (trim(lowercase(dim_names_v(i))) .ne. trim(lowercase(dim_names_u(i)))) then
      if (.not.(is_dimension_registered(fileobj_read_dd, dim_names_v(i)))) then
        call get_variable_units(fileobj_read_dd, dim_names_v(i), units_v(i))
        call validate_lat_lon_units(units_v(i), x_or_y, is_valid)
        if (is_valid) then
          call register_axis(fileobj_read_dd, dim_names_v(i), x_or_y, domain_position=v_pos)
        else
          call register_axis(fileobj_read_dd, dim_names_v(i), dim_sizes_v(i))
        endif
      endif
    endif
  enddo
  ! read the data
  dim_unlim_index = 0
  if (present(timelevel)) then
    do i=1,ndims
      if (is_dimension_unlimited(fileobj_read_dd, dim_names_u(i))) then
        dim_unlim_index = i
        exit
      endif
    enddo
    if (dim_unlim_index .gt. 0) then
      call read_data(fileobj_read_dd, u_fieldname,u_data, unlim_dim_level=timelevel)
      call read_data(fileobj_read_dd, v_fieldname, v_data, unlim_dim_level=timelevel)
    else
      call read_data(fileobj_read_dd, u_fieldname, u_data)
      call read_data(fileobj_read_dd, v_fieldname, v_data)
    endif
  else
    call read_data(fileobj_read_dd, u_fieldname, u_data)
    call read_data(fileobj_read_dd, v_fieldname, v_data)
  endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read_dd)) call fms2_close_file(fileobj_read_dd)
  endif
  if (allocated(dim_names_u)) deallocate(dim_names_u)
  if (allocated(dim_names_v)) deallocate(dim_names_v)
  if (allocated(dim_sizes_u)) deallocate(dim_sizes_u)
  if (allocated(dim_sizes_v)) deallocate(dim_sizes_v)
  if (allocated(units_u)) deallocate(units_u)
  if (allocated(units_v)) deallocate(units_v)

  ! Rescale the data that was read if necessary.
  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(MOM_Domain, u_data, scale)
    call rescale_comp_data(MOM_Domain, v_data, scale)
  endif ; endif

end subroutine MOM_read_vector_2d_fms2

!> This routine uses the fms2_io read_data interface to read a pair of distributed
!! 3-D data fields with names given by "[uv]_fieldname" from file "filename".  Valid values for
!! "stagger" include CGRID_NE, BGRID_NE, and AGRID.
subroutine MOM_read_vector_3d_fms2(filename, u_fieldname, v_fieldname, u_data, v_data, MOM_Domain, &
                                   timelevel, stagger, scale, leave_file_open)
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
  logical, optional, intent(in) :: leave_file_open !< if .true., leave file open
  ! local
  integer :: is, ie, js, je, i, dim_unlim, ndims
  integer :: u_pos, v_pos
  integer, allocatable :: dim_sizes_u(:), dim_sizes_v(:)
  character(len=32), allocatable :: dim_names_u(:), dim_names_v(:), units_u(:), units_v(:)
  character(len=1) :: x_or_y
  logical :: is_valid
  logical :: file_open_success ! .true. if open file is successful
  logical :: close_the_file ! indicates whether to close the file after MOM_read_vector is called; default is .true.

  close_the_file = .true.
  if (present(leave_file_open)) close_the_file = .not.(leave_file_open)

  ! open the file
  if (.not.(check_if_open(fileobj_read_dd))) then
    file_open_success = fms2_open_file(fileobj_read_dd, filename, "read", MOM_domain%mpp_domain, is_restart=.false.)
    if (.not. file_open_success) &
      call MOM_error(FATAL, "MOM_read_vector_3d_fms2: netcdf file "//trim(filename)//" not opened.")
  endif

  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == BGRID_NE) then ; u_pos = CORNER ; v_pos = CORNER
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif

  ndims = get_variable_num_dimensions(fileobj_read_dd, u_fieldname)
  allocate(dim_sizes_u(ndims))
  allocate(dim_sizes_v(ndims))
  allocate(dim_names_u(ndims))
  allocate(dim_names_v(ndims))
  allocate(units_u(ndims))
  allocate(units_v(ndims))

  units_u(:) = ""
  units_v(:) = ""
  dim_names_u(:) = ""
  dim_names_v(:) = ""

  call get_variable_size(fileobj_read_dd, u_fieldname, dim_sizes_u, broadcast=.true.)
  call get_variable_size(fileobj_read_dd, v_fieldname, dim_sizes_v, broadcast=.true.)
  call get_variable_dimension_names(fileobj_read_dd, u_fieldname, dim_names_u, broadcast=.true.)
  call get_variable_dimension_names(fileobj_read_dd, v_fieldname, dim_names_v, broadcast=.true.)

  do i=1,ndims
    ! register the u axes
    if (.not.(is_dimension_registered(fileobj_read_dd, dim_names_u(i)))) then
      call get_variable_units(fileobj_read_dd, dim_names_u(i), units_u(i))
      call validate_lat_lon_units(units_u(i), x_or_y, is_valid)
      if (is_valid) then
        call register_axis(fileobj_read_dd, dim_names_u(i), x_or_y, domain_position=u_pos)
      else
        call register_axis(fileobj_read_dd, dim_names_u(i), dim_sizes_u(i))
      endif
    endif
    ! Register the v axes if they differ from the u axes
    if (trim(lowercase(dim_names_v(i))) .ne. trim(lowercase(dim_names_u(i)))) then
      if (.not.(is_dimension_registered(fileobj_read_dd, dim_names_v(i)))) then
        call get_variable_units(fileobj_read_dd, dim_names_v(i), units_v(i))
        call validate_lat_lon_units(units_v(i), x_or_y, is_valid)
        if (is_valid) then
          call register_axis(fileobj_read_dd, dim_names_v(i), x_or_y, domain_position=v_pos)
        else
          call register_axis(fileobj_read_dd, dim_names_v(i), dim_sizes_v(i))
        endif
      endif
    endif
  enddo
  ! read the data
  dim_unlim = 0
  if (present(timelevel)) then
    do i=1,ndims
      if (is_dimension_unlimited(fileobj_read_dd, dim_names_u(i))) then
        dim_unlim = i
        exit
      endif
    enddo
    if (dim_unlim .gt. 0) then
      call read_data(fileobj_read_dd, u_fieldname, u_data, unlim_dim_level=timelevel)
      call read_data(fileobj_read_dd, v_fieldname, v_data, unlim_dim_level=timelevel)
    else
      call read_data(fileobj_read_dd, u_fieldname, u_data, edge_lengths=dim_sizes_u)
      call read_data(fileobj_read_dd, v_fieldname, v_data, edge_lengths=dim_sizes_v)
    endif
  else
    call read_data(fileobj_read_dd, u_fieldname, u_data, edge_lengths=dim_sizes_u)
    call read_data(fileobj_read_dd, v_fieldname, v_data, edge_lengths=dim_sizes_v)
  endif
  ! close the file
  if (close_the_file) then
    if (check_if_open(fileobj_read_dd)) call fms2_close_file(fileobj_read_dd)
  endif
  if (allocated(dim_names_u)) deallocate(dim_names_u)
  if (allocated(dim_names_v)) deallocate(dim_names_v)
  if (allocated(dim_sizes_u)) deallocate(dim_sizes_u)
  if (allocated(dim_sizes_v)) deallocate(dim_sizes_v)
  if (allocated(units_u)) deallocate(units_u)
  if (allocated(units_v)) deallocate(units_v)

  ! Rescale the data that was read if necessary.
  if (present(scale)) then ; if (scale /= 1.0) then
    call rescale_comp_data(MOM_Domain, u_data, scale)
    call rescale_comp_data(MOM_Domain, v_data, scale)
  endif ; endif

end subroutine MOM_read_vector_3d_fms2


!> check that latitude or longitude units are valid CF-compliant values
!! return true or false and x_or_y character value corresponding to the axis direction
subroutine validate_lat_lon_units(unit_string, x_or_y, units_are_valid)
character(len=*), intent(in) :: unit_string !< string of units
character(len=1), intent(out) :: x_or_y !< "x" for longitude or "y" latitude
logical, intent(out) :: units_are_valid !< .true. if units match acceptable values; default is .false.

select case (lowercase(trim(unit_string)))
  case ("degrees_north"); units_are_valid = .true.; x_or_y = "y"
  case ("degree_north"); units_are_valid = .true.; x_or_y = "y"
  case ("degrees_n"); units_are_valid = .true.; x_or_y = "y"
  case ("degree_n"); units_are_valid = .true.; x_or_y = "y"
  case ("degreen"); units_are_valid = .true.; x_or_y = "y"
  case ("degreesn"); units_are_valid = .true.; x_or_y = "y"
  case ("degrees_east"); units_are_valid = .true.; x_or_y = "x"
  case ("degree_east"); units_are_valid = .true.;x_or_y = "x"
  case ("degreese"); units_are_valid = .true.; x_or_y =  "x"
  case ("degreee"); units_are_valid = .true.; x_or_y =  "x"
  case ("degree_e"); units_are_valid = .true.; x_or_y =  "x"
  case ("degrees_e"); units_are_valid = .true.; x_or_y = "x"
  case default; units_are_valid = .false.; x_or_y = ""
end select

end subroutine validate_lat_lon_units

end module MOM_read_data_fms2
