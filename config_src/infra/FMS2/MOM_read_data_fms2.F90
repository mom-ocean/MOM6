!> This module contains routines that encapsulate common preparatory work for FMS2 read_data calls
module MOM_read_data_fms2

! This file is part of MOM6. See LICENSE.md for the license.
use MOM_error_infra,      only : MOM_error=>MOM_err, NOTE, FATAL, WARNING, is_root_PE
use MOM_domain_infra,     only : CENTER, CORNER, NORTH_FACE, EAST_FACE
use MOM_string_functions, only : lowercase
use fms2_io_mod,          only : FmsNetcdfDomainFile_t, FmsNetcdfFile_t, check_if_open
use fms2_io_mod,          only : variable_exists, get_num_variables, get_variable_names
use fms2_io_mod,          only : get_variable_size, get_variable_units
use fms2_io_mod,          only : get_variable_attribute, variable_att_exists
use fms2_io_mod,          only : get_variable_num_dimensions, get_variable_dimension_names
use fms2_io_mod,          only : is_dimension_unlimited, get_dimension_size
use fms2_io_mod,          only : is_dimension_registered, register_axis

implicit none ; private

public prepare_to_read_var, find_varname_in_file

contains

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

end module MOM_read_data_fms2
