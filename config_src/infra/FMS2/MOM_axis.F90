!> This module contains routines that define and register axes to files
module MOM_axis

! This file is part of MOM6. See LICENSE.md for the license.
use MOM_domains,          only : MOM_domain_type
use MOM_error_handler,    only : MOM_error, NOTE, FATAL, WARNING
use MOM_grid,             only : ocean_grid_type
use MOM_dyn_horgrid,      only : dyn_horgrid_type
use MOM_string_functions, only : lowercase
use MOM_verticalGrid,     only : verticalGrid_type
use fms2_io_mod,          only : is_dimension_registered, register_axis, is_dimension_unlimited
use fms2_io_mod,          only : FmsNetcdfDomainFile_t, FmsNetcdfFile_t, unlimited
use fms2_io_mod,          only : get_variable_size, get_variable_num_dimensions, check_if_open
use fms2_io_mod,          only : fms2_open_file=>open_file, fms2_close_file=>close_file
use fms2_io_mod,          only : get_variable_dimension_names, read_data, get_unlimited_dimension_name
use fms2_io_mod,          only : get_dimension_size
use mpp_domains_mod,      only : domain2d, CENTER, CORNER, NORTH_FACE=>NORTH, EAST_FACE=>EAST
use mpp_domains_mod,      only : mpp_get_compute_domain
use netcdf
implicit none ; private

public MOM_register_diagnostic_axis, get_var_dimension_metadata, get_time_units
public MOM_get_diagnostic_axis_data, MOM_register_variable_axes, get_time_index
public convert_checksum_to_string
!> A type for making arrays of pointers to real 1-d arrays
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

!> interface for registering axes associated with a variable to a netCDF file object
interface MOM_register_variable_axes
  module procedure MOM_register_variable_axes_subdomain
  module procedure MOM_register_variable_axes_full
end interface MOM_register_variable_axes

contains

!> register a MOM diagnostic axis to a domain-decomposed file
subroutine MOM_register_diagnostic_axis(fileObj, axisName, axisLength)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObj !< netCDF file object returned by call to open_file
  character(len=*), intent(in) :: axisName !< name of the axis to register to file
  integer, intent(in), optional :: axisLength !< length of axis/dimension ;only needed for Layer, Interface, Time,
                                              !! Period
  select case (trim(lowercase(axisName)))
    case ('latq'); call register_axis(fileObj,'latq','y', domain_position=NORTH_FACE)
    case ('lath'); call register_axis(fileObj,'lath','y', domain_position=CENTER)
    case ('lonq'); call register_axis(fileObj,'lonq','x', domain_position=EAST_FACE)
    case ('lonh'); call register_axis(fileObj,'lonh','x', domain_position=CENTER)
    case default
      if (.not. present(axisLength)) call MOM_error(FATAL,"MOM_io:register_diagnostic_axis: "//&
                        "An axis_length argument is required to register the axis "//trim(axisName))
      call register_axis(fileObj, trim(axisName), axisLength)
  end select
end subroutine MOM_register_diagnostic_axis


!> Get the horizontal grid, vertical grid, and/or time dimension names and lengths
!! for a single variable from the hor_grid, t_grid, and z_grid values returned by a prior call to query_vardesc
subroutine get_var_dimension_metadata(hor_grid, z_grid, t_grid_in, &
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

  use_lath = .false.
  use_lonh = .false.
  use_latq = .false.
  use_lonq = .false.

  ! set the ocean grid coordinates

  if (present(G)) then
    gridLatT => G%gridLatT ; gridLatB => G%gridLatB
    gridLonT => G%gridLonT ; gridLonB => G%gridLonB
    isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg
    IsgB = G%IsgB ; IegB = G%IegB ; JsgB = G%JsgB ; JegB = G%JegB

    call get_horizontal_grid_logic(hor_grid, use_lath, use_lonh, use_latq, use_lonq)
  elseif (present(dG)) then
    gridLatT => dG%gridLatT ; gridLatB => dG%gridLatB
    gridLonT => dG%gridLonT ; gridLonB => dG%gridLonB
    isg = dG%isg ; ieg = dG%ieg ; jsg = dG%jsg ; jeg = dG%jeg
    IsgB = dG%IsgB ; IegB = dG%IegB ; JsgB = dG%JsgB ; JegB = dG%JegB

    call get_horizontal_grid_logic(hor_grid, use_lath, use_lonh, use_latq, use_lonq)
  endif

  ! add longitude name to dimension name array
  if (use_lonh) then
    num_dims = num_dims+1
    dim_names(num_dims) = ""
    dim_names(num_dims)(1:len_trim("lonh")) = "lonh"
    dim_lengths(num_dims) = size(gridLonT(isg:ieg))
  elseif (use_lonq) then
    num_dims = num_dims+1
    dim_names(num_dims) = ""
    dim_names(num_dims)(1:len_trim("lonq")) = "lonq"
    dim_lengths(num_dims) = size(gridLonB(IsgB:IegB))
  endif
  ! add latitude name to dimension name array
  if (use_lath) then
    num_dims = num_dims+1
    dim_names(num_dims) = ""
    dim_names(num_dims)(1:len_trim("lath")) = "lath"
    dim_lengths(num_dims) = size(gridLatT(jsg:jeg))
  elseif (use_latq) then
    num_dims = num_dims+1
    dim_names(num_dims) = ""
    dim_names(num_dims)(1:len_trim("latq")) = "latq"
    dim_lengths(num_dims) = size(gridLatB(JsgB:JegB))
  endif

  if (present(GV)) then
    ! vertical grid
    select case (trim(z_grid))
      case ('L')
        num_dims = num_dims+1
        dim_names(num_dims) = ""
        dim_names(num_dims)(1:len_trim("Layer")) = "Layer"
        dim_lengths(num_dims) = GV%ke
      case ('i')
        num_dims = num_dims+1
        dim_names(num_dims) = ""
        dim_names(num_dims)(1:len_trim("Interface")) = "Interface"
        dim_lengths(num_dims) = GV%ke+1
      case ('1') ! Do nothing.
      case default
      call MOM_error(FATAL, "MOM_io: get_var_dimension_features: "//&
                     " has an unrecognized z_grid argument"//trim(z_grid))
    end select
  endif
  ! time
  t_grid = adjustl(t_grid_in)
  select case (t_grid(1:1))
    case ('s', 'a', 'm')
      num_dims = num_dims+1
      dim_names(num_dims) = ""
      dim_names(num_dims)(1:len_trim("Time")) = "Time"
      dim_lengths(num_dims) = unlimited
    case ('p')
      if (len_trim(t_grid(2:8)) <= 0) then
          call MOM_error(FATAL,"MOM_io:get_var_dimension_features: "//&
                           "No periodic axis length was specified in "//trim(t_grid))
      endif
      num_dims = num_dims+1
      dim_names(num_dims) = ""
      dim_names(num_dims)(1:len_trim("Period")) = "Period"
      dim_lengths(num_dims) = unlimited
    case ('1') ! Do nothing.
    case default
      call MOM_error(WARNING, "MOM_io: get_var_dimension_metadata: "//&
                     "Unrecognized t_grid "//trim(t_grid))
  end select
end subroutine get_var_dimension_metadata


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

  ! initialize axis_data_CS elements
  axis_data_CS%axis(axis_number)%name = ''
  axis_data_CS%axis(axis_number)%longname = ''
  axis_data_CS%axis(axis_number)%units = ''
  axis_data_CS%axis(axis_number)%horgrid_position = 0
  axis_data_CS%axis(axis_number)%is_domain_decomposed = .false.
  axis_data_CS%axis(axis_number)%positive = ''
  axis_data_CS%data(axis_number)%p => NULL()

  ! set the ocean grid coordinates and metadata
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

  select case(trim(lowercase(axis_name)))
    case('lath')
      if (associated(gridLatT)) &
        axis_data_CS%data(axis_number)%p=>gridLatT(jsg:jeg)

      axis_data_CS%axis(axis_number)%name = trim(axis_name)
      axis_data_CS%axis(axis_number)%longname = 'Latitude'
      axis_data_CS%axis(axis_number)%units = y_axis_units
      axis_data_CS%axis(axis_number)%horgrid_position = CENTER
      axis_data_CS%axis(axis_number)%is_domain_decomposed = .true.
    case('lonh')
      if (associated(gridLonT)) &
        axis_data_CS%data(axis_number)%p=>gridLonT(isg:ieg)

      axis_data_CS%axis(axis_number)%name = trim(axis_name)
      axis_data_CS%axis(axis_number)%horgrid_position  = CENTER
      axis_data_CS%axis(axis_number)%longname = 'Longitude'
      axis_data_CS%axis(axis_number)%units = x_axis_units
      axis_data_CS%axis(axis_number)%is_domain_decomposed  = .true.
    case('latq')
      if (associated(gridLatB)) &
        axis_data_CS%data(axis_number)%p=>gridLatB(JsgB:JegB)

      axis_data_CS%axis(axis_number)%name = trim(axis_name)
      axis_data_CS%axis(axis_number)%longname = 'Latitude'
      axis_data_CS%axis(axis_number)%units = y_axis_units
      axis_data_CS%axis(axis_number)%horgrid_position = NORTH_FACE
      axis_data_CS%axis(axis_number)%is_domain_decomposed = .true.
    case('lonq')
      if (associated(gridLonB)) &
        axis_data_CS%data(axis_number)%p=>gridLonB(IsgB:IegB)

        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Longitude'
        axis_data_CS%axis(axis_number)%units = x_axis_units
        axis_data_CS%axis(axis_number)%horgrid_position = EAST_FACE
        axis_data_CS%axis(axis_number)%is_domain_decomposed = .true.
    case('layer')
      if (present(GV)) then
        axis_data_CS%data(axis_number)%p=>GV%sLayer(1:GV%ke)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Layer pseudo-depth, -z*'
        axis_data_CS%axis(axis_number)%units = GV%zAxisUnits
        axis_data_CS%axis(axis_number)%positive = 'up'
      endif
    case('interface')
      if (present(GV)) then
        axis_data_CS%data(axis_number)%p=>GV%sInterface(1:GV%ke+1)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Interface pseudo-depth, -z*'
        axis_data_CS%axis(axis_number)%units = GV%zAxisUnits
        axis_data_CS%axis(axis_number)%positive = 'up'
      endif
    case('time')
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
    case('period')
      if (.not.(present(time_val))) &
        call MOM_error(FATAL, "MOM_axis::get_diagnostic_axis_data: requires a time_val argument "// &
                       "for "//trim(axis_name))
      axis_data_CS%data(axis_number)%p=>time_val
      axis_data_CS%axis(axis_number)%name = trim(axis_name)
      axis_data_CS%axis(axis_number)%longname = 'Periods for cyclical variables'
    case default
      call MOM_error(WARNING, "MOM_axis::get_diagnostic_axis_data:"//trim(axis_name)//" is an unrecognized axis")
  end select

end subroutine MOM_get_diagnostic_axis_data


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
        call MOM_error(FATAL, "MOM_axis:get_var_dimension_features "//&
                        "Unrecognized hor_grid argument "//trim(grid_string_id))
  end select
end subroutine get_horizontal_grid_logic

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

!> function to get the index of a time_value from a netCDF file
function get_time_index(filename, time_to_find) result (time_index)
  character(len=*) :: filename ! name of the file to read in
  real, intent(in) :: time_to_find ! time value to search for in file
  ! local
  type(fmsNetcdfFile_t) :: fileobj ! netCDF file object returned by open_file
  real, allocatable, dimension(:) :: file_times ! array of time values read from file
  integer :: dim_unlim_size, i, time_index
  character(len=nf90_max_name) :: dim_unlim_name ! name of the unlimited dimension in the file
  logical :: file_open_success

  time_index = 1
  dim_unlim_size = 0
  dim_unlim_name = ""
  file_open_success = .false.

  if (.not. check_if_open(fileobj)) &
    !call MOM_error(FATAL, "get_time_index_nodd: netcdf file object must be open.")
    file_open_success=fms2_open_file(fileobj, trim(filename), "read", is_restart=.false.)

  call get_unlimited_dimension_name(fileobj, dim_unlim_name)
  call get_dimension_size(fileObj, trim(dim_unlim_name), dim_unlim_size)
  ! time index will be one more than the unlimited dimension size if the time_to_find is not in the file
  if (dim_unlim_size .gt. 0) then
    time_index = dim_unlim_size+1
    allocate(file_times(dim_unlim_size))
    call read_data(fileobj,trim(dim_unlim_name), file_times)

    do i=1,dim_unlim_size
      if (ABS(file_times(i)-time_to_find) .gt. TINY(time_to_find)) then
        continue
      else
        time_index = i
        exit
      endif
    enddo
    deallocate(file_times)
  endif
  if (check_if_open(fileobj)) call fms2_close_file(fileobj)
end function get_time_index

!> register axes associated with a variable from a domain-decomposed netCDF file that are mapped to
!! a sub-domain (e.g., a supergrid).
!> \note The user must specify units for variables with longitude/x-axis and/or latitude/y-axis axes to obtain
!! the correct domain decomposition for the data buffer.
subroutine MOM_register_variable_axes_subdomain(fileObj, variableName, io_domain, position)
  type(FmsNetcdfFile_t), intent(inout) :: fileObj !< netCDF file object returned by call to open_file
  character(len=*),  intent(in) :: variableName !< name of the variable
  type(domain2d),    intent(in) :: io_domain !< type that contains the mpp io domain
  integer, optional, intent(in) :: position  !< A flag indicating where this data is discretized

  ! Local variables
  character(len=40) :: units ! units corresponding to a specific variable dimension
  character(len=40), allocatable, dimension(:) :: dim_names ! variable dimension names
  integer :: i, isg, ieg, isc, iec, jsg, jeg, jsc, jec, xlen, ylen
  integer :: ndims ! number of dimensions
  integer :: pos   ! Discrete variable position. Default is CENTER
  integer, allocatable, dimension(:) :: dimSizes ! variable dimension sizes

  if (.not. check_if_open(fileObj)) call MOM_error(FATAL,"MOM_axis:register_variable_axes_subdomain: The fileObj "// &
                                                  " has not been opened. Call fms2_open_file(fileObj,...) "// &
                                                  "before passing the fileObj argument to this function.")

  ! get variable dimension names and lengths
  ndims = get_variable_num_dimensions(fileObj, trim(variableName))
  allocate(dimSizes(ndims))
  allocate(dim_names(ndims))
  call get_variable_size(fileObj, trim(variableName), dimSizes, broadcast=.true.)
  call get_variable_dimension_names(fileObj, trim(variableName), dim_names)

  ! Get the lengths of the global indicies, using the discrete position of this variable
  pos = CORNER ; if (present(position)) pos = position
  call mpp_get_compute_domain(io_domain, xsize=xlen, ysize=ylen, position=pos)
  ! register the axes
  !>\note: This is not a comprehensive check for all possible supported horizontal axes associated with variables
  !! read from netCDF files. Developers should add/remove cases as needed.
  do i=1,ndims
    !if (.not.(is_dimension_registered(fileObj, trim(dim_names(i))))) then
      select case(trim(lowercase(dim_names(i))))
        case ("grid_x_t")
          call register_axis(fileObj, trim(dim_names(i)), xlen)
        case ("nx")
          call register_axis(fileObj, trim(dim_names(i)), xlen)
        case("nxp")
          call register_axis(fileObj, trim(dim_names(i)), xlen)
        case("longitude")
          call register_axis(fileObj, trim(dim_names(i)), xlen)
        case("long")
          call register_axis(fileObj, trim(dim_names(i)), xlen)
        case("lon")
          call register_axis(fileObj, trim(dim_names(i)), xlen)
        case("lonh")
          call register_axis(fileObj, trim(dim_names(i)), xlen)
        case("lonq")
          call register_axis(fileObj, trim(dim_names(i)), xlen)
        case("xh")
          call register_axis(fileObj, trim(dim_names(i)), xlen)
        case ("grid_y_t")
          call register_axis(fileObj, trim(dim_names(i)), ylen)
        case ("ny")
          call register_axis(fileObj, trim(dim_names(i)), ylen)
        case("nyp")
          call register_axis(fileObj, trim(dim_names(i)), ylen)
        case("latitude")
          call register_axis(fileObj, trim(dim_names(i)), ylen)
        case("lat")
          call register_axis(fileObj, trim(dim_names(i)), ylen)
        case("lath")
          call register_axis(fileObj, trim(dim_names(i)), ylen)
        case("latq")
          call register_axis(fileObj, trim(dim_names(i)), ylen)
        case("yh")
          call register_axis(fileObj, trim(dim_names(i)), ylen)
        case default ! assumes that the axis is not domain-decomposed
          if (.not. is_dimension_unlimited(fileObj, trim(dim_names(i)))) &
            call MOM_error(WARNING,"MOM_register_variable_axes_subdomain: the axis "//trim(dim_names(i))//&
              "is not included in the valid x and y dimension cases. If the code hangs, check the whether "//&
              "an x or y axis is being registered as a non-domain-decomposed variable, "//&
              "and add it to the accepted cases if necessary.")
          call register_axis(fileObj, trim(dim_names(i)), dimSizes(i))
      end select
   ! endif
  enddo

  if (allocated(dimSizes)) deallocate(dimSizes)
  if (allocated(dim_names)) deallocate(dim_names)
end subroutine MOM_register_variable_axes_subdomain

!> register axes associated with a variable from a domain-decomposed netCDF file
!> @note The user must specify units for variables with longitude/x-axis and/or latitude/y-axis axes
!!  to obtain the correct domain decomposition for the data buffer.
subroutine MOM_register_variable_axes_full(fileObj, variableName, position)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObj !< netCDF file object returned by call to open_file
  character(len=*),  intent(in) :: variableName !< name of the variable
  integer, optional, intent(in) :: position  !< A flag indicating where this data is discretized

  ! Local variables
  character(len=40) :: units ! units corresponding to a specific variable dimension
  character(len=40), allocatable, dimension(:) :: dim_names ! variable dimension names
  integer :: i
  integer :: ndims ! number of dimensions
  integer :: xPos, yPos ! domain positions for x and y axes. Default is CENTER
  integer, allocatable, dimension(:) :: dimSizes ! variable dimension sizes

  if (.not. check_if_open(fileObj)) call MOM_error(FATAL,"MOM_axis:register_variable_axes: The fileObj has "// &
                                                  "not been opened. Call fms2_open_file(fileObj,...) before "// &
                                                  "passing the fileObj argument to this function.")
  xpos = CENTER ; ypos = CENTER
  if (present(position)) then
    if ((position == CORNER) .or. (position == EAST_FACE)) xpos = EAST_FACE
    if ((position == CORNER) .or. (position == NORTH_FACE)) ypos = NORTH_FACE
  endif

  ! get variable dimension names and lengths
  ndims = get_variable_num_dimensions(fileObj, trim(variableName))
  allocate(dimSizes(ndims))
  allocate(dim_names(ndims))
  call get_variable_size(fileObj, trim(variableName), dimSizes)
  call get_variable_dimension_names(fileObj, trim(variableName), dim_names)
  ! register the axes
  !>@note: This is not a comprehensive check for all possible supported horizontal axes associated with variables
  !! read from netCDF files. Developers should add/remove cases as needed.
  do i=1,ndims
    if (.not.(is_dimension_registered(fileobj, trim(dim_names(i))))) then
      select case(trim(lowercase(dim_names(i))))
        case ("grid_x_t")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case ("nx")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case("nxp")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case("longitude")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case("long")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case("lon")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case("lonh")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case("lonq")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case("xh")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case("i")
          call register_axis(fileObj, trim(dim_names(i)),"x", domain_position=xPos)
        case ("grid_y_t")
          call register_axis(fileObj, trim(dim_names(i)),"y", domain_position=yPos)
        case ("ny")
          call register_axis(fileObj, trim(dim_names(i)),"y", domain_position=yPos)
        case("nyp")
          call register_axis(fileObj, trim(dim_names(i)),"y", domain_position=yPos)
        case("latitude")
          call register_axis(fileObj, trim(dim_names(i)),"y", domain_position=yPos)
        case("lat")
          call register_axis(fileObj, trim(dim_names(i)),"y", domain_position=yPos)
        case("lath")
          call register_axis(fileObj, trim(dim_names(i)),"y", domain_position=yPos)
        case("latq")
          call register_axis(fileObj, trim(dim_names(i)),"y", domain_position=yPos)
        case("yh")
          call register_axis(fileObj, trim(dim_names(i)),"y", domain_position=yPos)
        case("j")
          call register_axis(fileObj, trim(dim_names(i)),"y", domain_position=yPos)
        case default ! assumes that the axis is not domain-decomposed
          if (.not. is_dimension_unlimited(fileObj, trim(dim_names(i)))) &
            call MOM_error(WARNING,"MOM_register_variable_axes_full: the axis "//trim(dim_names(i))//" is not "//&
              "included in the valid x and y dimension cases. If the code hangs, check the whether "//&
              "an x or y axis is being registered as a non-domain-decomposed variable, "//&
              "and add it to the accepted cases if necessary.")
          call register_axis(fileObj, trim(dim_names(i)), dimSizes(i))
      end select
    endif
  enddo

  deallocate(dimSizes)
  deallocate(dim_names)
end subroutine MOM_register_variable_axes_full


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


end module MOM_axis
