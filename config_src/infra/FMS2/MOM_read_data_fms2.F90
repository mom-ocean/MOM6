!> This module contains routines that wrap the fms2 read_data calls
module MOM_read_data_fms2

! This file is part of MOM6. See LICENSE.md for the license.
use MOM_error_infra,      only : MOM_error=>MOM_err, NOTE, FATAL, WARNING, is_root_PE
use MOM_domain_infra,     only : MOM_domain_type, AGRID, BGRID_NE, CGRID_NE
use MOM_domain_infra,     only : domain2d, CENTER, CORNER, NORTH_FACE, EAST_FACE
use MOM_string_functions, only : lowercase
use fms2_io_mod,          only : FmsNetcdfDomainFile_t, FmsNetcdfFile_t
use fms2_io_mod,          only : fms2_open_file => open_file, fms2_close_file => close_file
use fms2_io_mod,          only : get_num_variables, get_variable_names, check_if_open
use fms2_io_mod,          only : read_data, variable_exists, get_variable_size, get_variable_units
use fms2_io_mod,          only : get_variable_attribute, attribute_exists => variable_att_exists
use fms2_io_mod,          only : get_variable_num_dimensions, get_variable_dimension_names
use fms2_io_mod,          only : is_dimension_unlimited, get_dimension_size
use fms2_io_mod,          only : is_dimension_registered, register_axis

implicit none ; private

public prepare_to_read_var
! public MOM_read_data_scalar, MOM_read_data_2d_noDD, MOM_read_data_1d_noDD

contains

!> Find the case-insensitive name match with a variable in a domain-decomposed file-set
!! opening the file(s) as necessary, prepare FMS2 to read this variable, and return some
!! information needed to call read_data correctly for this variable and file.
subroutine prepare_to_read_var(fileobj, fieldname, domain, err_header, filename, var_to_read, &
                               has_time_dim, timelevel, position)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj     !< A handle to an FMS2 file object, that
                                                            !! will be opened if necessary
  character(len=*),            intent(in)    :: fieldname   !< The variable name to seek in the file
  type(MOM_domain_type),       intent(in)    :: domain      !< MOM domain attribute with the mpp_domain decomposition
  character(len=*),            intent(in)    :: err_header  !< A descriptive prefix for error messages
  character(len=*),            intent(in)    :: filename    !< The name of the file to read
  character(len=*),            intent(out)   :: var_to_read !< The variable name to read from the file
  logical,           optional, intent(out)   :: has_time_dim !< Indicates whether fieldname has a time dimension
  integer,           optional, intent(in)    :: timelevel   !< A time level to read
  integer,           optional, intent(in)    :: position    !< A flag indicating where this variable is discretized

  ! Local variables
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: variable_found ! Is a case-insensitive version of the variable found in the netCDF file?
  character(len=96), allocatable, dimension(:) :: var_names !< array for names of variables in a netCDF
                                                            !! file opened to read
  character(len=96), allocatable :: dim_names(:) ! variable dimension names
  integer :: nvars ! The number of variables in the file.
  integer :: i, dim_unlim_size, num_var_dims, time_dim

  ! Open the file if necessary
  if (.not.(check_if_open(fileobj))) then
    file_open_success = fms2_open_file(fileobj, filename, "read", domain%mpp_domain, is_restart=.false.)
    if (.not.file_open_success) call MOM_error(FATAL, trim(err_header)//" failed to open "//trim(filename))
  endif

  ! Search for the variable in the file, looking for the case-sensitive name first.
  if (variable_exists(fileobj, trim(fieldname))) then
    var_to_read = trim(fieldname)
    variable_found = .true.
  else  ! Look for case-insensitive variable name matches.
    var_to_read = ""
    variable_found = .false.

    nvars = get_num_variables(fileobj)
    if (nvars < 1) call MOM_error(FATAL, "nvars is less than 1 for file "//trim(filename))
    allocate(var_names(nvars))
    call get_variable_names(fileobj, var_names)

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
!> @note The user must specify units for variables with longitude/x-axis and/or latitude/y-axis axes
!!  to obtain the correct domain decomposition for the data buffer.
subroutine MOM_register_variable_axes(fileObj, variableName, filename, position)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObj !< Handle to an open FMS2 netCDF file object
  character(len=*),  intent(in) :: variableName !< name of the variable
  character(len=*),  intent(in) :: filename     !< The name of the file to read
  integer, optional, intent(in) :: position     !< A flag indicating where this data is discretized

  ! Local variables
  character(len=40) :: units ! units corresponding to a specific variable dimension
  character(len=40), allocatable, dimension(:) :: dim_names ! variable dimension names
  integer, allocatable, dimension(:) :: dimSizes ! variable dimension sizes
  logical, allocatable, dimension(:) :: is_x ! Is this a (likely domain-decomposed) x-axis
  logical, allocatable, dimension(:) :: is_y ! Is this a (likely domain-decomposed) y-axis
  logical, allocatable, dimension(:) :: is_t ! Is this a time axis or another unlimited axis
  integer :: ndims ! number of dimensions
  integer :: i
  integer :: xPos, yPos ! domain positions for x and y axes. Default is CENTER

  if (.not. check_if_open(fileObj)) call MOM_error(FATAL,"MOM_axis:register_variable_axes: The fileObj has "// &
                                                  "not been opened. Call fms2_open_file(fileObj,...) before "// &
                                                  "passing the fileObj argument to this function.")
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

  deallocate(dimSizes)
  deallocate(dim_names)
  deallocate(is_x, is_y, is_t)
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

  integer :: i
  character(len=256) :: cartesian ! A flag indicating a Cartesian direction - usually a single character.
  character(len=512) :: dim_list  ! A concatenated list of dimension names.
  character(len=40)  :: units ! units corresponding to a specific variable dimension
  logical :: x_found, y_found ! Indicate whether an x- or y- dimension have been found.

  x_found = .false. ; y_found = .false.
  is_x(:) = .false. ; is_y(:) = .false.
  do i=1,ndims
    is_t(i) = is_dimension_unlimited(fileObj, trim(dim_names(i)))
    ! First look for indicative variable attributes
    if (.not.is_t(i)) then
      if (variable_exists(fileobj, trim(dim_names(i)))) then
        if (attribute_exists(fileobj, trim(dim_names(i)), "cartesian_axis")) then
          call get_variable_attribute(fileobj, trim(dim_names(i)), "cartesian_axis", cartesian)
          cartesian = adjustl(cartesian)
          if ((index(cartesian, "X") == 1) .or. (index(cartesian, "x") == 1)) is_x(i) = .true.
          if ((index(cartesian, "Y") == 1) .or. (index(cartesian, "y") == 1)) is_y(i) = .true.
          if ((index(cartesian, "T") == 1) .or. (index(cartesian, "t") == 1)) is_t(i) = .true.
          ! if (is_root_pe() .and. is_x(i)) &
          !   call MOM_error(NOTE, "X-dimension determined from cartesian_axis for "//trim(dim_names(i)))
          ! if (is_root_pe() .and. is_y(i)) &
          !   call MOM_error(NOTE, "Y-dimension determined from cartesian_axis for "//trim(dim_names(i)))
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

  if (.not.(x_found .and. y_found) .and. (ndims>2) .or. ((ndims==2) .and. .not.is_t(ndims))) then
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

!===== Everything below this pertains to reading non-decomposed variables ===!
!===== using FMS2 interfaces will probably be discarded eventually. =========!

!!> This routine calls the fms_io read_data subroutine to read a scalar (0-D) field named "fieldname"
!! from file "filename".
subroutine MOM_read_data_scalar(filename, fieldname, data, timelevel, scale, leave_file_open)
  character(len=*),         intent(in)    :: filename  !< The name of the file to read
  character(len=*),         intent(in)    :: fieldname !< The variable name of the data in the file
  real,                     intent(inout) :: data      !< The variable to read from read_data
  integer,        optional, intent(in)    :: timelevel !< time level to read
  real,           optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied by
  logical,        optional, intent(in)    :: leave_file_open !< if .true., leave file open

  ! Local variables
  type(FmsNetcdfFile_t) :: fileobj       ! A handle to a simple netCDF file
  logical :: close_the_file ! indicates whether to close the file after read_data is called.
  character(len=96) :: var_to_read ! variable to read from the netcdf file
  character(len=48) :: err_header ! A preamble for error messages

  err_header = "MOM_read_data_fms2:MOM_read_data_scalar: "

  ! Find the matching variable name in the file, opening it and reading metadata if necessary.
  call find_varname_in_file(fileobj, fieldname, err_header, filename, var_to_read)

  ! read the data
  if (present(timelevel)) then
    call read_data(fileobj, trim(var_to_read), data, unlim_dim_level=timelevel)
  else
    call read_data(fileobj, trim(var_to_read), data)
  endif

  ! Close the file, if necessary
  close_the_file = .true. ; if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  if (close_the_file .and. check_if_open(fileobj)) call fms2_close_file(fileobj)

  ! Rescale the data that was read if necessary.
  if (present(scale)) then ; if (scale /= 1.0) then
    data = scale*data
  endif ; endif

end subroutine MOM_read_data_scalar

!> This routine calls the fms_io read_data subroutine to read 1-D non-domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call.
subroutine MOM_read_data_1d_noDD(filename, fieldname, data, start_index, &
                                 edge_lengths, timelevel, scale, leave_file_open)
  character(len=*),         intent(in)    :: filename  !< The name of the file to read
  character(len=*),         intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:),       intent(inout) :: data      !< The 1-dimensional data array to pass to read_data
  integer, dimension(1), optional, intent(in) :: start_index !< starting index of data buffer. Default is 1
  integer, dimension(1), optional, intent(in) :: edge_lengths !< number of data values to read in; default is
                                                              !! the variable size
  integer,        optional, intent(in)    :: timelevel !< time level to read
  real,           optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied by
  logical,        optional, intent(in)    :: leave_file_open !< if .true., leave file open

  ! Local variables
  type(FmsNetcdfFile_t) :: fileobj       ! A handle to a simple netCDF file
  logical :: close_the_file ! indicates whether to close the file after read_data is called.
  integer :: time_dim       ! The dimension position of a variables unlimited time axis, or -1 if it has none.
  integer, parameter :: ndim = 1 ! The dimensionality of the array being read
  integer, dimension(ndim) :: start, nread ! indices for first data value and number of values to read
  character(len=96) :: var_to_read ! variable to read from the netcdf file
  character(len=48) :: err_header ! A preamble for error messages

  err_header = "MOM_read_data_fms2:MOM_read_data_1d_noDD: "

  ! Find the matching case-insensitive variable name in the file, opening the file if necessary.
  call find_varname_in_file(fileobj, fieldname, err_header, filename, var_to_read)

  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  start(:) = 1 ; if (present(start_index)) start(:) = start_index(:)
  nread(:) = shape(data) ; if (present(edge_lengths)) nread(:) = edge_lengths(:)

  time_dim = -1
  if (present(timelevel)) then
    time_dim = get_time_dim(fileobj, var_to_read, err_header, filename, timelevel)
    if (time_dim == ndim) then ; nread(ndim) = 1 ; start(ndim) = timelevel ; endif
  endif

  ! read the data
  if (time_dim > 0) then
    call read_data(fileobj, trim(var_to_read), data, corner=start, edge_lengths=nread, &
                   unlim_dim_level=timelevel)
  else
    call read_data(fileobj, trim(var_to_read), data, corner=start, edge_lengths=nread)
  endif

  ! Close the file, if necessary
  close_the_file = .true. ; if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  if (close_the_file .and. check_if_open(fileobj)) call fms2_close_file(fileobj)

  ! Rescale the data that was read if necessary.
  if (present(scale)) then ; if (scale /= 1.0) then
    data(:) = scale*data(:)
  endif ; endif

end subroutine MOM_read_data_1d_noDD

!> This routine calls the fms_io read_data subroutine to read a 2-D non-domain-decomposed data field named "fieldname"
!! from file "filename". The routine multiplies the data by "scale" if the optional argument is included in the call.
subroutine MOM_read_data_2d_noDD(filename, fieldname, data, start_index, &
                                 edge_lengths, timelevel, position, scale, leave_file_open)
  character(len=*),         intent(in)    :: filename  !< The name of the file to read
  character(len=*),         intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:,:),     intent(inout) :: data      !< The 2-dimensional data array to pass to read_data
  integer, dimension(2), optional, intent(in) :: start_index !< starting indices of data buffer. Default is 1
  integer, dimension(2), optional, intent(in) :: edge_lengths !< number of data values to read in.
                                                              !! Default values are the variable dimension sizes
  integer,        optional, intent(in)    :: timelevel !< time level to read
  integer,        optional, intent(in)    :: position  !< A flag indicating where this data is located
  real,           optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied by
  logical,        optional, intent(in)    :: leave_file_open !< if .true., leave file open

  ! Local variables
  type(FmsNetcdfFile_t) :: fileobj       ! A handle to a simple netCDF file
  logical :: close_the_file ! indicates whether to close the file after read_data is called.
  integer :: time_dim       ! The dimension position of a variables unlimited time axis, or -1 if it has none.
  integer, parameter :: ndim = 2 ! The dimensionality of the array being read
  integer, dimension(ndim) :: start, nread ! indices for first data value and number of values to read
  character(len=96) :: var_to_read ! variable to read from the netcdf file
  character(len=48) :: err_header ! A preamble for error messages

  err_header = "MOM_read_data_fms2:MOM_read_data_2d_DD: "

  ! Find the matching case-insensitive variable name in the file, opening the file if necessary.
  call find_varname_in_file(fileobj, fieldname, err_header, filename, var_to_read)

  ! set the start and nread values that will be passed as the read_data corner and edge_lengths arguments
  start(:) = 1 ; if (present(start_index)) start(:) = start_index(:)
  nread(:) = shape(data) ; if (present(edge_lengths)) nread(:) = edge_lengths(:)

  time_dim = -1
  if (present(timelevel)) then
    time_dim = get_time_dim(fileobj, var_to_read, err_header, filename, timelevel)
    if (time_dim == ndim) then ; nread(ndim) = 1 ; start(ndim) = timelevel ; endif
  endif

  ! read the data
  if (time_dim > 0) then
    call read_data(fileobj, trim(var_to_read), data, corner=start, edge_lengths=nread, &
                   unlim_dim_level=timelevel)
  else
    call read_data(fileobj, trim(var_to_read), data, corner=start, edge_lengths=nread)
  endif

  ! Close the file, if necessary
  close_the_file = .true. ; if (present(leave_file_open)) close_the_file = .not.(leave_file_open)
  if (close_the_file .and. check_if_open(fileobj)) call fms2_close_file(fileobj)

  ! Rescale the data that was read if necessary.
  if (present(scale)) then ; if (scale /= 1.0) then
    data(:,:) = scale*data(:,:)
  endif ; endif

end subroutine MOM_read_data_2d_noDD


!> Find the case-sensitive name of the variable in a netCDF file with a case-insensitive name match.
subroutine find_varname_in_file(fileobj, fieldname, err_header, filename, var_to_read)
  type(FmsNetcdfFile_t),       intent(inout) :: fileobj     !< A handle to a file object, that
                                                            !! will be opened if necessary
  character(len=*),            intent(in)    :: fieldname   !< The variable name to seek in the file
  character(len=*),            intent(in)    :: err_header  !< A descriptive prefix for error messages
  character(len=*),            intent(in)    :: filename    !< The name of the file to read
  character(len=*),            intent(out)   :: var_to_read !< The variable name to read from the file

  ! Local variables
  logical :: file_open_success !.true. if call to open_file is successful
  logical :: variable_found ! Is a case-insensitive version of the variable found in the netCDF file?
  character(len=96), allocatable, dimension(:) :: var_names !< array for names of variables in a netCDF
                                                            !! file opened to read
  integer :: nvars ! The number of variables in the file.
  integer :: i

  var_to_read = ""

  ! Open the file if necessary
  if (.not.(check_if_open(fileobj))) then
    file_open_success = fms2_open_file(fileobj, filename, "read", is_restart=.false.)
    if (.not.file_open_success) call MOM_error(FATAL, trim(err_header)//" failed to open "//trim(filename))
  endif

  if (variable_exists(fileobj, fieldname)) then
    var_to_read = fieldname
  else
    variable_found = .false.
    nvars = get_num_variables(fileobj)
    if (nvars < 1) call MOM_error(FATAL, "nvars is less than 1 for file "//trim(filename))
    allocate(var_names(nvars))
    call get_variable_names(fileobj, var_names)

    ! search for the variable in the file
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

end subroutine find_varname_in_file

!> Return the number of the time dimension for a variable in an open non-domain-decomposed file,
!! or -1 if it has no time (or other unlimited) dimension.
integer function get_time_dim(fileobj, var_to_read, err_header, filename, timelevel)
  type(FmsNetcdfFile_t),       intent(in) :: fileobj     !< A handle to an open file object
  character(len=*),            intent(in) :: var_to_read !< The variable name to read from the file
  character(len=*),            intent(in) :: err_header  !< A descriptive prefix for error messages
  character(len=*),            intent(in) :: filename    !< The name of the file to read
  integer,           optional, intent(in) :: timelevel   !< A time level to read

  ! Local variables
  integer :: i, dim_unlim_size, num_var_dims
  character(len=96), allocatable :: dim_names(:) ! variable dimension names

  num_var_dims = get_variable_num_dimensions(fileobj, trim(var_to_read))
  allocate(dim_names(num_var_dims)) ; dim_names(:) = ""
  call get_variable_dimension_names(fileobj, trim(var_to_read), dim_names)

  get_time_dim = -1
  do i=1,num_var_dims
    if (is_dimension_unlimited(fileobj, dim_names(i))) then
      get_time_dim = i
      if (present(timelevel)) then
        call get_dimension_size(fileobj, dim_names(i), dim_unlim_size)
        if (timelevel > dim_unlim_size) call MOM_error(FATAL, trim(err_header)//&
              "Attempting to read a time level of "//trim(var_to_read)//&
              " that exceeds the size of "//trim(filename))
      endif
      exit
    endif
  enddo
  if (get_time_dim < 0) &
    call MOM_error(WARNING, trim(err_header)//"time level specified, but the variable "//&
                   trim(var_to_read)//" does not have an unlimited dimension in "//trim(filename))
  deallocate(dim_names)

end function get_time_dim

end module MOM_read_data_fms2
