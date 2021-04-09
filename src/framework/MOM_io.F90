!> This module contains I/O framework code
module MOM_io

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_array_transform,  only : allocate_rotated_array, rotate_array
use MOM_domains,          only : MOM_domain_type, domain1D, broadcast, get_domain_components
use MOM_domains,          only : rescale_comp_data, AGRID, BGRID_NE, CGRID_NE
use MOM_dyn_horgrid,      only : dyn_horgrid_type
use MOM_ensemble_manager, only : get_ensemble_id
use MOM_error_handler,    only : MOM_error, NOTE, FATAL, WARNING, is_root_PE
use MOM_file_parser,      only : log_version, param_file_type
use MOM_grid,             only : ocean_grid_type
use MOM_io_infra,         only : MOM_read_data, MOM_read_vector, read_field_chksum
use MOM_io_infra,         only : read_data=>MOM_read_data ! read_data will be removed soon.
use MOM_io_infra,         only : file_type, file_exists, get_file_info, get_file_fields
use MOM_io_infra,         only : open_file, open_ASCII_file, close_file, flush_file, file_is_open
use MOM_io_infra,         only : get_field_size, fieldtype, field_exists, get_field_atts
use MOM_io_infra,         only : get_file_times, axistype, get_axis_data, get_filename_suffix
use MOM_io_infra,         only : write_field, write_metadata, write_version
use MOM_io_infra,         only : MOM_namelist_file, check_namelist_error, io_infra_init, io_infra_end
use MOM_io_infra,         only : APPEND_FILE, ASCII_FILE, MULTIPLE, NETCDF_FILE, OVERWRITE_FILE
use MOM_io_infra,         only : READONLY_FILE, SINGLE_FILE, WRITEONLY_FILE
use MOM_io_infra,         only : CENTER, CORNER, NORTH_FACE, EAST_FACE
use MOM_string_functions, only : lowercase, slasher
use MOM_verticalGrid,     only : verticalGrid_type

use iso_fortran_env,      only : int32, int64, stdout_iso=>output_unit, stderr_iso=>error_unit
use netcdf,               only : NF90_open, NF90_inq_varid, NF90_inq_varids, NF90_inquire, NF90_close
use netcdf,               only : NF90_inquire_variable, NF90_get_var, NF90_get_att, NF90_inquire_attribute
use netcdf,               only : NF90_strerror, NF90_inquire_dimension
use netcdf,               only : NF90_NOWRITE, NF90_NOERR, NF90_GLOBAL, NF90_ENOTATT, NF90_CHAR

implicit none ; private

! These interfaces are actually implemented in this file.
public :: create_file, reopen_file, cmor_long_std, ensembler, MOM_io_init
public :: MOM_write_field, var_desc, modify_vardesc, query_vardesc, position_from_horgrid
public :: open_namelist_file, check_namelist_error, check_nml_error
public :: get_var_sizes, verify_variable_units, num_timelevels, read_variable, read_attribute
public :: open_file_to_read, close_file_to_read
! The following are simple pass throughs of routines from MOM_io_infra or other modules.
public :: file_exists, open_file, open_ASCII_file, close_file, flush_file, file_type
public :: get_file_info, field_exists, get_file_fields, get_file_times, get_filename_appendix
public :: fieldtype, field_size, get_field_atts
public :: axistype, get_axis_data
public :: MOM_read_data, MOM_read_vector, read_field_chksum
public :: slasher, write_field, write_version_number
public :: io_infra_init, io_infra_end
! This is used to set up information descibing non-domain-decomposed axes.
public :: axis_info, set_axis_info, delete_axis_info
! This is used to set up global file attributes
public :: attribute_info, set_attribute_info, delete_attribute_info
! This API is here just to support potential use by non-FMS drivers, and should not persist.
public :: read_data
!> These encoding constants are used to indicate the file format
public :: ASCII_FILE, NETCDF_FILE
!> These encoding constants are used to indicate whether the file is domain decomposed
public :: MULTIPLE, SINGLE_FILE
!> These encoding constants are used to indicate the access mode for a file
public :: APPEND_FILE, OVERWRITE_FILE, READONLY_FILE, WRITEONLY_FILE
!> These encoding constants are used to indicate the discretization position of a variable
public :: CENTER, CORNER, NORTH_FACE, EAST_FACE

!> Write a registered field to an output file, potentially with rotation
interface MOM_write_field
  module procedure MOM_write_field_4d
  module procedure MOM_write_field_3d
  module procedure MOM_write_field_2d
  module procedure MOM_write_field_1d
  module procedure MOM_write_field_0d
end interface MOM_write_field

!> Read an entire named variable from a named netCDF file using netCDF calls directly, rather
!! than any infrastructure routines and broadcast it from the root PE to the other PEs.
interface read_variable
  module procedure read_variable_0d, read_variable_0d_int
  module procedure read_variable_1d, read_variable_1d_int
end interface read_variable

!> Read a global or variable attribute from a named netCDF file using netCDF calls
!! directly, in some cases reading from the root PE before broadcasting to the other PEs.
interface read_attribute
  module procedure read_attribute_str, read_attribute_real
  module procedure read_attribute_int32, read_attribute_int64
end interface read_attribute

!> Type for describing a 3-d variable for output
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
  character(len=32)  :: dim_names(5)       !< The names in the file of the axes for this variable
  integer            :: position = -1      !< An integer encoding the horizontal position, it may
                                           !! CENTER, CORNER, EAST_FACE, NORTH_FACE, or 0.
end type vardesc

!> Type that stores information that can be used to create a non-decomposed axis.
type :: axis_info ; private
  character(len=32)  :: name = ""       !< The name of this axis for use in files
  character(len=256) :: longname = ""   !< A longer name describing this axis
  character(len=48)  :: units = ""      !< The units of the axis labels
  character(len=8)   :: cartesian = "N" !< A variable indicating which direction
                                        !! this axis corresponds with. Valid values
                                        !! include 'X', 'Y', 'Z', 'T', and 'N' for none.
  integer            :: sense = 0       !< This is 1 for axes whose values increase upward, or -1
                                        !! if they increase downward.  The default, 0, is ignored.
  integer            :: ax_size = 0     !< The number of elements in this axis
  real, allocatable, dimension(:) :: ax_data !< The values of the data on the axis.
end type axis_info

!> Type that stores for a global file attribute
type :: attribute_info ; private
  character(len=:), allocatable :: name    !< The name of this attribute
  character(len=:), allocatable :: att_val !< The values of this attribute
end type attribute_info


integer, public :: stdout = stdout_iso  !< standard output unit
integer, public :: stderr = stderr_iso  !< standard output unit

contains

!> Routine creates a new NetCDF file.  It also sets up fieldtype
!! structures that describe this file and variables that will
!! later be written to this file.
subroutine create_file(IO_handle, filename, vars, novars, fields, threading, timeunit, &
                       G, dG, GV, checksums, extra_axes, global_atts)
  type(file_type),       intent(inout) :: IO_handle  !< Handle for a files or fileset that is to be
                                                     !! opened or reopened for writing
  character(len=*),      intent(in)    :: filename   !< full path to the file to create
  type(vardesc),         intent(in)    :: vars(:)    !< structures describing fields written to filename
  integer,               intent(in)    :: novars     !< number of fields written to filename
  type(fieldtype),       intent(inout) :: fields(:)  !< array of fieldtypes for each variable
  integer, optional,     intent(in)    :: threading  !< SINGLE_FILE or MULTIPLE
  real, optional,        intent(in)    :: timeunit   !< length of the units for time [s]. The
                                                     !! default value is 86400.0, for 1 day.
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if the new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  integer(kind=int64),     optional, intent(in) :: checksums(:,:) !< checksums of vars
  type(axis_info),         optional, intent(in) :: extra_axes(:)  !< Types with information about
                                                     !! some axes that might be used in this file
  type(attribute_info),    optional, intent(in) :: global_atts(:) !< Global attributes to
                                                     !! write to this file

  logical        :: use_lath, use_lonh, use_latq, use_lonq, use_time
  logical        :: use_layer, use_int, use_periodic
  logical        :: one_file, domain_set, dim_found
  logical, dimension(:), allocatable :: use_extra_axis
  type(axistype) :: axis_lath, axis_latq, axis_lonh, axis_lonq
  type(axistype) :: axis_layer, axis_int, axis_time, axis_periodic
  type(axistype),  dimension(:), allocatable :: more_axes ! Axes generated from extra_axes
  type(axistype) :: axes(5)     ! The axes of a variable
  type(MOM_domain_type), pointer :: Domain => NULL()
  type(domain1d) :: x_domain, y_domain
  integer        :: position, numaxes, pack, thread, k, n, m
  integer        :: num_extra_dims ! The number of extra possible dimensions from extra_axes
  integer        :: isg, ieg, jsg, jeg, IsgB, IegB, JsgB, JegB
  integer        :: var_periods, num_periods=0
  real, dimension(:), allocatable :: axis_val
  real, pointer, dimension(:) :: &
    gridLatT => NULL(), & ! The latitude or longitude of T or B points for
    gridLatB => NULL(), & ! the purpose of labeling the output axes.
    gridLonT => NULL(), gridLonB => NULL()
  character(len=40) :: time_units, x_axis_units, y_axis_units
  character(len=8)  :: t_grid, t_grid_read
  character(len=64) :: ax_name(5)  ! The axis names of a variable

  use_lath  = .false. ; use_lonh     = .false.
  use_latq  = .false. ; use_lonq     = .false.
  use_time  = .false. ; use_periodic = .false.
  use_layer = .false. ; use_int      = .false.
  num_extra_dims = 0
  if (present(extra_axes)) then
    num_extra_dims = size(extra_axes)
    if (num_extra_dims > 0) then
      allocate(use_extra_axis(num_extra_dims)) ; use_extra_axis = .false.
      allocate(more_axes(num_extra_dims))
    endif
  endif

  thread = SINGLE_FILE
  if (PRESENT(threading)) thread = threading

  domain_set = .false.
  if (present(G)) then
    domain_set = .true. ; Domain => G%Domain
    gridLatT => G%gridLatT ; gridLatB => G%gridLatB
    gridLonT => G%gridLonT ; gridLonB => G%gridLonB
    x_axis_units = G%x_axis_units ; y_axis_units = G%y_axis_units
    isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg
    IsgB = G%IsgB ; IegB = G%IegB ; JsgB = G%JsgB ; JegB = G%JegB
  elseif (present(dG)) then
    domain_set = .true. ; Domain => dG%Domain
    gridLatT => dG%gridLatT ; gridLatB => dG%gridLatB
    gridLonT => dG%gridLonT ; gridLonB => dG%gridLonB
    x_axis_units = dG%x_axis_units ; y_axis_units = dG%y_axis_units
    isg = dG%isg ; ieg = dG%ieg ; jsg = dG%jsg ; jeg = dG%jeg
    IsgB = dG%IsgB ; IegB = dG%IegB ; JsgB = dG%JsgB ; JegB = dG%JegB
  endif

  one_file = .true.
  if (domain_set) one_file = (thread == SINGLE_FILE)

  if (one_file) then
    call open_file(IO_handle, filename, OVERWRITE_FILE, threading=thread)
  else
    call open_file(IO_handle, filename, OVERWRITE_FILE, MOM_domain=Domain)
  endif

! Define the coordinates.
  do k=1,novars
    position = vars(k)%position
    if (position == -1) position = position_from_horgrid(vars(k)%hor_grid)
    select case (position)
      case (CENTER) ; use_lath = .true. ; use_lonh = .true.
      case (CORNER) ; use_latq = .true. ; use_lonq = .true.
      case (EAST_FACE) ; use_lath = .true. ; use_lonq = .true.
      case (NORTH_FACE) ; use_latq = .true. ; use_lonh = .true.
      case (0) ! Do nothing.
      case default
        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//" has an unrecognized value of postion")
    end select
    select case (vars(k)%z_grid)
      case ('L') ; use_layer = .true.
      case ('i') ; use_int = .true.
      case ('1') ! Do nothing.
      case default
        call MOM_error(FATAL, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized z_grid "//trim(vars(k)%z_grid))
    end select
    t_grid = adjustl(vars(k)%t_grid)
    select case (t_grid(1:1))
      case ('s', 'a', 'm') ; use_time = .true.
      case ('p') ; use_periodic = .true.
        if (len_trim(t_grid(2:8)) <= 0) call MOM_error(FATAL, &
          "MOM_io create_file: No periodic axis length was specified in "//&
          trim(vars(k)%t_grid) // " in the periodic axes of variable "//&
          trim(vars(k)%name)//" in file "//trim(filename))
        var_periods = -9999999
        t_grid_read = adjustl(t_grid(2:8))
        read(t_grid_read,*) var_periods
        if (var_periods == -9999999) call MOM_error(FATAL, &
          "MOM_io create_file: Failed to read the number of periods from "//&
          trim(vars(k)%t_grid) // " in the periodic axes of variable "//&
          trim(vars(k)%name)//" in file "//trim(filename))
        if (var_periods < 1) call MOM_error(FATAL, "MOM_io create_file: "//&
           "variable "//trim(vars(k)%name)//" in file "//trim(filename)//&
           " uses a periodic time axis, and must have a positive "//&
           "value for the number of periods in "//vars(k)%t_grid )
        if ((num_periods > 0) .and. (var_periods /= num_periods)) &
          call MOM_error(FATAL, "MOM_io create_file: "//&
            "Only one value of the number of periods can be used in the "//&
            "create_file call for file "//trim(filename)//".  The second is "//&
            "variable "//trim(vars(k)%name)//" with t_grid "//vars(k)%t_grid )

        num_periods = var_periods
      case ('1') ! Do nothing.
      case default
        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized t_grid "//trim(vars(k)%t_grid))
    end select

    do n=1,5 ; if (len_trim(vars(k)%dim_names(n)) > 0) then
      dim_found = .false.
      do m=1,num_extra_dims
        if (lowercase(trim(vars(k)%dim_names(n))) == lowercase(trim(extra_axes(m)%name))) then
          use_extra_axis(m) = .true.
          dim_found = .true.
          exit
        endif
      enddo
      if (.not.dim_found) call MOM_error(FATAL, "Unable to find a match for dimension "//&
          trim(vars(k)%dim_names(n))//" for variable "//trim(vars(k)%name)//" in file "//trim(filename))
    endif ; enddo
  enddo

  if ((use_lath .or. use_lonh .or. use_latq .or. use_lonq)) then
    if (.not.domain_set) call MOM_error(FATAL, "create_file: "//&
      "An ocean_grid_type or dyn_horgrid_type is required to create a file with a horizontal coordinate.")

    call get_domain_components(Domain, x_domain, y_domain)
  endif
  if ((use_layer .or. use_int) .and. .not.present(GV)) call MOM_error(FATAL, &
    "create_file: A vertical grid type is required to create a file with a vertical coordinate.")

  if (use_lath) &
    call write_metadata(IO_handle, axis_lath, name="lath", units=y_axis_units, longname="Latitude", &
                        cartesian='Y', domain=y_domain, data=gridLatT(jsg:jeg))

  if (use_lonh) &
    call write_metadata(IO_handle, axis_lonh, name="lonh", units=x_axis_units, longname="Longitude", &
                        cartesian='X', domain=x_domain, data=gridLonT(isg:ieg))

  if (use_latq) &
    call write_metadata(IO_handle, axis_latq, name="latq", units=y_axis_units, longname="Latitude", &
                        cartesian='Y', domain=y_domain, data=gridLatB(JsgB:JegB), edge_axis=.true.)

  if (use_lonq) &
    call write_metadata(IO_handle, axis_lonq, name="lonq", units=x_axis_units, longname="Longitude", &
                        cartesian='X', domain=x_domain, data=gridLonB(IsgB:IegB), edge_axis=.true.)

  if (use_layer) &
    call write_metadata(IO_handle, axis_layer, name="Layer", units=trim(GV%zAxisUnits), &
                        longname="Layer "//trim(GV%zAxisLongName), cartesian='Z', &
                        sense=1, data=GV%sLayer(1:GV%ke))

  if (use_int) &
    call write_metadata(IO_handle, axis_int, name="Interface", units=trim(GV%zAxisUnits), &
                        longname="Interface "//trim(GV%zAxisLongName), cartesian='Z', &
                        sense=1, data=GV%sInterface(1:GV%ke+1))

  if (use_time) then ; if (present(timeunit)) then
    ! Set appropriate units, depending on the value.
    if (timeunit < 0.0) then
      time_units = "days" ! The default value.
    elseif ((timeunit >= 0.99) .and. (timeunit < 1.01)) then
      time_units = "seconds"
    elseif ((timeunit >= 3599.0) .and. (timeunit < 3601.0)) then
      time_units = "hours"
    elseif ((timeunit >= 86399.0) .and. (timeunit < 86401.0)) then
      time_units = "days"
    elseif ((timeunit >= 3.0e7) .and. (timeunit < 3.2e7)) then
      time_units = "years"
    else
      write(time_units,'(es8.2," s")') timeunit
    endif

    call write_metadata(IO_handle, axis_time, name="Time", units=time_units, longname="Time", cartesian='T')
  else
    call write_metadata(IO_handle, axis_time, name="Time", units="days", longname="Time", cartesian='T')
  endif ; endif

  if (use_periodic) then
    if (num_periods <= 1) call MOM_error(FATAL, "MOM_io create_file: "//&
      "num_periods for file "//trim(filename)//" must be at least 1.")
    ! Define a periodic axis with unit labels.
    allocate(axis_val(num_periods))
    do k=1,num_periods ; axis_val(k) = real(k) ; enddo
    call write_metadata(IO_handle, axis_periodic, name="Period", units="nondimensional", &
                        longname="Periods for cyclical variables", cartesian='T', data=axis_val)
    deallocate(axis_val)
  endif

  do m=1,num_extra_dims ; if (use_extra_axis(m)) then
    if (allocated(extra_axes(m)%ax_data)) then
      call write_metadata(IO_handle, more_axes(m), name=extra_axes(m)%name, units=extra_axes(m)%units, &
                          longname=extra_axes(m)%longname, cartesian=extra_axes(m)%cartesian, &
                          sense=extra_axes(m)%sense, data=extra_axes(m)%ax_data)
    elseif (trim(extra_axes(m)%cartesian) == "T") then
      call write_metadata(IO_handle, more_axes(m), name=extra_axes(m)%name, units=extra_axes(m)%units, &
                          longname=extra_axes(m)%longname, cartesian=extra_axes(m)%cartesian)
    else
      ! FMS requires that non-time axes have variables that label their values, even if they are trivial.
      allocate (axis_val(extra_axes(m)%ax_size))
      do k=1,extra_axes(m)%ax_size ; axis_val(k) = real(k) ; enddo
      call write_metadata(IO_handle, more_axes(m), name=extra_axes(m)%name, units=extra_axes(m)%units, &
                          longname=extra_axes(m)%longname, cartesian=extra_axes(m)%cartesian, &
                          sense=extra_axes(m)%sense, data=axis_val)
      deallocate(axis_val)
    endif
  endif ; enddo

  do k=1,novars
    numaxes = 0
    position = vars(k)%position
    if (position == -1) position = position_from_horgrid(vars(k)%hor_grid)
    select case (position)
      case (CENTER)
        numaxes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_lath ; ax_name(1) = "lonh" ; ax_name(2) = "lath"
      case (CORNER)
        numaxes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_latq ; ax_name(1) = "lonq" ; ax_name(2) = "latq"
      case (EAST_FACE)
        numaxes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_lath ; ax_name(1) = "lonq" ; ax_name(2) = "lath"
      case (NORTH_FACE)
        numaxes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_latq ; ax_name(1) = "lonh" ; ax_name(2) = "latq"
      case (0) ! Do nothing.
      case default
        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized position, hor_grid = "//trim(vars(k)%hor_grid))
    end select
    select case (vars(k)%z_grid)
      case ('L') ; numaxes = numaxes+1 ; axes(numaxes) = axis_layer ; ax_name(numaxes) = "Layer"
      case ('i') ; numaxes = numaxes+1 ; axes(numaxes) = axis_int ; ax_name(numaxes) = "Interface"
      case ('1') ! Do nothing.
      case default
        call MOM_error(FATAL, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized z_grid "//trim(vars(k)%z_grid))
    end select

    do n=1,numaxes
      if ( (len_trim(vars(k)%dim_names(n)) > 0) .and. (trim(ax_name(n)) /= trim(vars(k)%dim_names(n))) ) &
        call MOM_error(WARNING, "MOM_io create_file: dimension "//trim(ax_name(n))//&
               " of variable "//trim(vars(k)%name)//" in "//trim(filename)//&
               " is being set inconsistently as "//trim(vars(k)%dim_names(n)))
    enddo
    do n=numaxes+1,5 ; if (len_trim(vars(k)%dim_names(n)) > 0) then
      dim_found = .false.
      do m=1,num_extra_dims
        if (lowercase(trim(vars(k)%dim_names(n))) == lowercase(trim(extra_axes(m)%name))) then
          numaxes = numaxes+1 ; axes(numaxes) = more_axes(m)
          exit
        endif
      enddo
    endif ; enddo

    t_grid = adjustl(vars(k)%t_grid)
    select case (t_grid(1:1))
      case ('s', 'a', 'm') ; numaxes = numaxes+1 ; axes(numaxes) = axis_time
      case ('p')           ; numaxes = numaxes+1 ; axes(numaxes) = axis_periodic
      case ('1') ! Do nothing.
      case default
        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized t_grid "//trim(vars(k)%t_grid))
    end select

    pack = 1
    if (present(checksums)) then
      call write_metadata(IO_handle, fields(k), axes(1:numaxes), vars(k)%name, vars(k)%units, &
                          vars(k)%longname, pack=pack, checksum=checksums(k,:))
    else
      call write_metadata(IO_handle, fields(k), axes(1:numaxes), vars(k)%name, vars(k)%units, &
                          vars(k)%longname, pack=pack)
    endif
  enddo

  if (present(global_atts)) then
    do n=1,size(global_atts)
      if (allocated(global_atts(n)%name) .and. allocated(global_atts(n)%att_val)) &
        call write_metadata(IO_handle, global_atts(n)%name, global_atts(n)%att_val)
    enddo
  endif

  ! Now actualy write the variables with the axis label values
  if (use_lath) call write_field(IO_handle, axis_lath)
  if (use_latq) call write_field(IO_handle, axis_latq)
  if (use_lonh) call write_field(IO_handle, axis_lonh)
  if (use_lonq) call write_field(IO_handle, axis_lonq)
  if (use_layer) call write_field(IO_handle, axis_layer)
  if (use_int) call write_field(IO_handle, axis_int)
  if (use_periodic) call write_field(IO_handle, axis_periodic)
  do m=1,num_extra_dims ; if (use_extra_axis(m)) then
    call write_field(IO_handle, more_axes(m))
  endif ; enddo

  if (num_extra_dims > 0) then
    deallocate(use_extra_axis, more_axes)
  endif

end subroutine create_file


!> This routine opens an existing NetCDF file for output.  If it
!! does not find the file, a new file is created.  It also sets up
!! structures that describe this file and the variables that will
!! later be written to this file.
subroutine reopen_file(IO_handle, filename, vars, novars, fields, threading, timeunit, &
                       G, dG, GV, extra_axes, global_atts)
  type(file_type),       intent(inout) :: IO_handle  !< Handle for a file or fileset that is to be
                                                     !! opened or reopened for writing
  character(len=*),      intent(in)    :: filename   !< full path to the file to create
  type(vardesc),         intent(in)    :: vars(:)    !< structures describing fields written to filename
  integer,               intent(in)    :: novars     !< number of fields written to filename
  type(fieldtype),       intent(inout) :: fields(:)  !< array of fieldtypes for each variable
  integer, optional,     intent(in)    :: threading  !< SINGLE_FILE or MULTIPLE
  real, optional,        intent(in)    :: timeunit   !< length of the units for time [s]. The
                                                     !! default value is 86400.0, for 1 day.
  type(ocean_grid_type),   optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                     !! is required if a new file uses any
                                                     !! horizontal grid axes.
  type(dyn_horgrid_type),  optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                     !! is required if a new file uses any
                                                     !! horizontal grid axes.
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if a new file uses any
                                                     !! vertical grid axes.
  type(axis_info),         optional, intent(in) :: extra_axes(:)  !< Types with information about
                                                     !! some axes that might be used in this file
  type(attribute_info),    optional, intent(in) :: global_atts(:) !< Global attributes to
                                                     !! write to this file

  type(MOM_domain_type), pointer :: Domain => NULL()
  character(len=200) :: check_name, mesg
  integer :: length, nvar, thread
  logical :: exists, one_file, domain_set

  thread = SINGLE_FILE
  if (PRESENT(threading)) thread = threading

  check_name = filename
  length = len(trim(check_name))
  if (check_name(length-2:length) /= ".nc") check_name = trim(check_name)//".nc"
  if (thread /= SINGLE_FILE) check_name = trim(check_name)//".0000"

  inquire(file=check_name,EXIST=exists)

  if (.not.exists) then
    call create_file(IO_handle, filename, vars, novars, fields, threading, timeunit, &
                     G=G, dG=dG, GV=GV, extra_axes=extra_axes, global_atts=global_atts)
  else

    domain_set = .false.
    if (present(G)) then
      domain_set = .true. ; Domain => G%Domain
    elseif (present(dG)) then
      domain_set = .true. ; Domain => dG%Domain
    endif

    one_file = .true.
    if (domain_set) one_file = (thread == SINGLE_FILE)

    if (one_file) then
      call open_file(IO_handle, filename, APPEND_FILE, threading=thread)
    else
      call open_file(IO_handle, filename, APPEND_FILE, MOM_domain=Domain)
    endif
    if (.not.file_is_open(IO_handle)) return

    call get_file_info(IO_handle, nvar=nvar)

    if (nvar == -1) then
      write (mesg,*) "Reopening file ",trim(filename)," apparently had ",nvar,&
                     " variables. Clobbering and creating file with ",novars," instead."
      call MOM_error(WARNING,"MOM_io: "//mesg)
      call create_file(IO_handle, filename, vars, novars, fields, threading, timeunit, &
                       G=G, dG=dG, GV=GV, extra_axes=extra_axes, global_atts=global_atts)
    elseif (nvar /= novars) then
      write (mesg,*) "Reopening file ",trim(filename)," with ",novars,&
                     " variables instead of ",nvar,"."
      call MOM_error(FATAL,"MOM_io: "//mesg)
    endif

    if (nvar > 0) call get_file_fields(IO_handle, fields(1:nvar))

    ! Check for inconsistent field names...
!    do i=1,nvar
!      call get_field_atts(fields(i), name)
!      !if (trim(name) /= trim(vars%name)) then
!      !  write (mesg, '("Reopening file ",a," variable ",a," is called ",a,".")',&
!      !         trim(filename), trim(vars%name), trim(name))
!      !  call MOM_error(NOTE, "MOM_io: "//trim(mesg))
!      !endif
!    enddo
  endif

end subroutine reopen_file


!> This function determines how many time levels a variable has in a file.
function num_timelevels(filename, varname, min_dims) result(n_time)
  character(len=*),  intent(in) :: filename   !< name of the file to read
  character(len=*),  intent(in) :: varname    !< variable whose number of time levels
                                              !! are to be returned
  integer, optional, intent(in) :: min_dims   !< The minimum number of dimensions a variable must have
                                              !! if it has a time dimension.  If the variable has 1 less
                                              !! dimension than this, then 0 is returned.
  integer :: n_time                           !< number of time levels varname has in filename

  character(len=256) :: msg
  integer :: ncid, status, varid, ndims
  integer :: sizes(8)

  n_time = -1

  ! To do almost the same via MOM_io_infra calls, we could do the following:
  !   found = field_exists(filename, varname)
  !   if (found) then
  !     call open_file(ncid, filename, action=READONLY_FILE, form=NETCDF_FILE, threading=MULTIPLE)
  !     call get_file_info(ncid, ntime=n_time)
  !   endif
  ! However, this does not handle the case where the time axis for the variable is not the record
  ! axis and min_dims is not used.

  call get_var_sizes(filename, varname, ndims, sizes, match_case=.false., caller="num_timelevels")

  n_time = sizes(ndims)

  if (present(min_dims)) then
    if (ndims < min_dims-1) then
      write(msg, '(I3)') min_dims
      call MOM_error(WARNING, "num_timelevels: variable "//trim(varname)//" in file "//&
        trim(filename)//" has fewer than min_dims = "//trim(msg)//" dimensions.")
      n_time = -1
    elseif (ndims == min_dims - 1) then
      n_time = 0
    endif
  endif

end function num_timelevels


!> get_var_sizes returns the number and size of dimensions associate with a variable in a file.
!! Usually only the root PE does the read, and then the information is broadcast
subroutine get_var_sizes(filename, varname, ndims, sizes, match_case, caller, all_read, dim_names, ncid_in)
  character(len=*),      intent(in)  :: filename   !< Name of the file to read, used here in messages
  character(len=*),      intent(in)  :: varname    !< The variable name, used here for messages
  integer,               intent(out) :: ndims      !< The number of dimensions to the variable
  integer, dimension(:), intent(out) :: sizes      !< The dimension sizes, or 0 for extra values
  logical,     optional, intent(in)  :: match_case !< If false, allow for variables name matches to be
                                                   !! case insensitive, but take a perfect match if
                                                   !! found.  The default is true.
  character(len=*), optional, intent(in) :: caller !< The name of a calling routine for use in error messages
  logical,     optional, intent(in)  :: all_read   !< If present and true, all PEs that call this
                                                   !! routine actually do the read, otherwise only
                                                   !! root PE reads and then it broadcasts the results.
  character(len=*), dimension(:), &
               optional, intent(out) :: dim_names  !< The names of the dimensions for this variable
  integer,     optional, intent(in)  :: ncid_in    !< The netCDF ID of an open file.  If absent, the
                                                   !! file is opened and closed within this routine.

  logical :: do_read, do_broadcast
  integer, allocatable :: size_msg(:)  ! An array combining the number of dimensions and the sizes.
  integer :: n, nval

  do_read = is_root_pe()
  if (present(all_read)) do_read = all_read .or. do_read
  do_broadcast = .true. ; if (present(all_read)) do_broadcast = .not.all_read

  if (do_read) call read_var_sizes(filename, varname, ndims, sizes, match_case, caller, dim_names, ncid_in)

  if (do_broadcast) then
    ! Distribute the sizes from the root PE.
    nval = size(sizes) + 1

    allocate(size_msg(nval))
    size_msg(1) = ndims
    do n=2,nval ; size_msg(n) = sizes(n-1) ; enddo

    call broadcast(size_msg, nval, blocking=.true.)

    ndims = size_msg(1)
    do n=2,nval ;  sizes(n-1) = size_msg(n) ; enddo
    deallocate(size_msg)

    if (present(dim_names)) then
      nval = min(ndims, size(dim_names))
      call broadcast(dim_names(1:nval), len(dim_names(1)), blocking=.true.)
    endif
  endif

end subroutine get_var_sizes

!> read_var_sizes returns the number and size of dimensions associate with a variable in a file.
!! Every processor for which this is called does the reading.
subroutine read_var_sizes(filename, varname, ndims, sizes, match_case, caller, dim_names, ncid_in)
  character(len=*),      intent(in)  :: filename   !< Name of the file to read, used here in messages
  character(len=*),      intent(in)  :: varname    !< The variable name, used here for messages
  integer,               intent(out) :: ndims      !< The number of dimensions to the variable
  integer, dimension(:), intent(out) :: sizes      !< The dimension sizes, or 0 for extra values
  logical,     optional, intent(in)  :: match_case !< If false, allow for variables name matches to be
                                                   !! case insensitive, but take a perfect match if
                                                   !! found.  The default is true.
  character(len=*), &
               optional, intent(in)  :: caller     !< The name of a calling routine for use in error messages
  character(len=*), dimension(:), &
               optional, intent(out) :: dim_names  !< The names of the dimensions for this variable
  integer,     optional, intent(in)  :: ncid_in    !< The netCDF ID of an open file.  If absent, the
                                                   !! file is opened and closed within this routine.

  character(len=256) :: hdr, dimname
  integer, allocatable :: dimids(:)
  integer :: varid, ncid, n, status
  logical :: success
  hdr = "get_var_size: " ; if (present(caller)) hdr = trim(hdr)//": "
  sizes(:) = 0 ; ndims = -1

  if (present(ncid_in)) then
    ncid = ncid_in
  else
    call open_file_to_read(filename, ncid, success=success)
    if (.not.success) return
  endif

  ! Get the dimension sizes of the variable varname.
  call get_varid(varname, ncid, filename, varid, match_case=match_case)
  if (varid < 0) return

  status = NF90_inquire_variable(ncid, varid, ndims=ndims)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING, trim(hdr) // trim(NF90_STRERROR(status)) //&
      " Getting number of dimensions of "//trim(varname)//" in "//trim(filename))
    return
  endif
  if (ndims < 1) return

  allocate(dimids(ndims))
  status = NF90_inquire_variable(ncid, varid, dimids=dimids(1:ndims))
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING, trim(hdr) // trim(NF90_STRERROR(status)) //&
      " Getting dimension IDs for "//trim(varname)//" in "//trim(filename))
    deallocate(dimids) ; return
  endif

  do n = 1, min(ndims,size(sizes))
    status = NF90_Inquire_Dimension(ncid, dimids(n), name=dimname, len=sizes(n))
    if (status /= NF90_NOERR) call MOM_error(WARNING, trim(hdr) // trim(NF90_STRERROR(status)) //&
        " Getting dimension length for "//trim(varname)//" in "//trim(filename))
    if (present(dim_names)) then
      if (n <= size(dim_names)) dim_names = trim(dimname)
    endif
  enddo
  deallocate(dimids)

  if (.not.present(ncid_in)) call close_file_to_read(ncid, filename)

end subroutine read_var_sizes

!> Read a real scalar variable from a netCDF file with the root PE, and broadcast the
!! results to all the other PEs.
subroutine read_variable_0d(filename, varname, var, ncid_in, scale)
  character(len=*),  intent(in)    :: filename !< The name of the file to read
  character(len=*),  intent(in)    :: varname  !< The variable name of the data in the file
  real,              intent(inout) :: var      !< The scalar into which to read the data
  integer, optional, intent(in)    :: ncid_in  !< The netCDF ID of an open file.  If absent, the
                                               !! file is opened and closed within this routine
  real,    optional, intent(in)    :: scale    !< A scaling factor that the variable is
                                               !! multiplied by before it is returned

  integer :: varid, ncid, rc
  character(len=256) :: hdr
  hdr = "read_variable_0d"

  if (is_root_pe()) then
    if (present(ncid_in)) then
      ncid = ncid_in
    else
      call open_file_to_read(filename, ncid)
    endif

    call get_varid(varname, ncid, filename, varid, match_case=.false.)
    if (varid < 0) call MOM_error(FATAL, "Unable to get netCDF varid for "//trim(varname)//&
                                         " in "//trim(filename))
    rc = NF90_get_var(ncid, varid, var)
    if (rc /= NF90_NOERR) call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //&
          " Difficulties reading "//trim(varname)//" from "//trim(filename))

    if (.not.present(ncid_in)) call close_file_to_read(ncid, filename)

    if (present(scale)) var = scale * var
  endif

  call broadcast(var, blocking=.true.)
end subroutine read_variable_0d

!> Read a 1-d real variable from a netCDF file with the root PE, and broadcast the
!! results to all the other PEs.
subroutine read_variable_1d(filename, varname, var, ncid_in, scale)
  character(len=*),   intent(in)    :: filename !< The name of the file to read
  character(len=*),   intent(in)    :: varname  !< The variable name of the data in the file
  real, dimension(:), intent(inout) :: var      !< The 1-d array into which to read the data
  integer,  optional, intent(in)    :: ncid_in  !< The netCDF ID of an open file.  If absent, the
                                                !! file is opened and closed within this routine
  real,     optional, intent(in)    :: scale    !< A scaling factor that the variable is
                                                !! multiplied by before it is returned

  integer :: varid, ncid, rc
  character(len=256) :: hdr
  hdr = "read_variable_1d"

  if (is_root_pe()) then
    if (present(ncid_in)) then
      ncid = ncid_in
    else
      call open_file_to_read(filename, ncid)
    endif

    call get_varid(varname, ncid, filename, varid, match_case=.false.)
    if (varid < 0) call MOM_error(FATAL, "Unable to get netCDF varid for "//trim(varname)//&
                                         " in "//trim(filename))
    rc = NF90_get_var(ncid, varid, var)
    if (rc /= NF90_NOERR) call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //&
          " Difficulties reading "//trim(varname)//" from "//trim(filename))

    if (.not.present(ncid_in)) call close_file_to_read(ncid, filename)

    if (present(scale)) then ; if (scale /= 1.0) then
      var(:) = scale * var(:)
    endif ; endif
  endif

  call broadcast(var, size(var), blocking=.true.)
end subroutine read_variable_1d

!> Read a integer scalar variable from a netCDF file with the root PE, and broadcast the
!! results to all the other PEs.
subroutine read_variable_0d_int(filename, varname, var, ncid_in)
  character(len=*),  intent(in)    :: filename !< The name of the file to read
  character(len=*),  intent(in)    :: varname  !< The variable name of the data in the file
  integer,           intent(inout) :: var      !< The scalar into which to read the data
  integer, optional, intent(in)    :: ncid_in  !< The netCDF ID of an open file.  If absent, the
                                               !! file is opened and closed within this routine.

  integer :: varid, ncid, rc
  character(len=256) :: hdr
  hdr = "read_variable_0d_int"

  if (is_root_pe()) then
    if (present(ncid_in)) then
      ncid = ncid_in
    else
      call open_file_to_read(filename, ncid)
    endif

    call get_varid(varname, ncid, filename, varid, match_case=.false.)
    if (varid < 0) call MOM_error(FATAL, "Unable to get netCDF varid for "//trim(varname)//&
                                         " in "//trim(filename))
    rc = NF90_get_var(ncid, varid, var)
    if (rc /= NF90_NOERR) call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //&
          " Difficulties reading "//trim(varname)//" from "//trim(filename))

    if (.not.present(ncid_in)) call close_file_to_read(ncid, filename)
  endif

  call broadcast(var, blocking=.true.)
end subroutine read_variable_0d_int

!> Read a 1-d integer variable from a netCDF file with the root PE, and broadcast the
!! results to all the other PEs.
subroutine read_variable_1d_int(filename, varname, var, ncid_in)
  character(len=*),      intent(in)    :: filename !< The name of the file to read
  character(len=*),      intent(in)    :: varname  !< The variable name of the data in the file
  integer, dimension(:), intent(inout) :: var      !< The 1-d array into which to read the data
  integer, optional,     intent(in)    :: ncid_in  !< The netCDF ID of an open file.  If absent, the
                                                   !! file is opened and closed within this routine.

  integer :: varid, ncid, rc
  character(len=256) :: hdr
  hdr = "read_variable_1d_int"

  if (is_root_pe()) then
    if (present(ncid_in)) then
      ncid = ncid_in
    else
      call open_file_to_read(filename, ncid)
    endif

    call get_varid(varname, ncid, filename, varid, match_case=.false.)
    if (varid < 0) call MOM_error(FATAL, "Unable to get netCDF varid for "//trim(varname)//&
                                         " in "//trim(filename))
    rc = NF90_get_var(ncid, varid, var)
    if (rc /= NF90_NOERR) call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //&
          " Difficulties reading "//trim(varname)//" from "//trim(filename))

    if (.not.present(ncid_in)) call close_file_to_read(ncid, filename)
  endif

  call broadcast(var, size(var), blocking=.true.)
end subroutine read_variable_1d_int

!> Read a character-string global or variable attribute
subroutine read_attribute_str(filename, attname, att_val, varname, found, all_read, ncid_in)
  character(len=*),           intent(in)  :: filename !< Name of the file to read
  character(len=*),           intent(in)  :: attname  !< Name of the attribute to read
  character(:), allocatable,  intent(out) :: att_val  !< The value of the attribute
  character(len=*), optional, intent(in)  :: varname  !< The name of the variable whose attribute will
                                                      !! be read. If missing, read a global attribute.
  logical,          optional, intent(out) :: found    !< Returns true if the attribute is found
  logical,          optional, intent(in)  :: all_read !< If present and true, all PEs that call this
                                                      !! routine actually do the read, otherwise only
                                                      !! root PE reads and then broadcasts the results.
  integer,          optional, intent(in)  :: ncid_in  !< The netCDF ID of an open file.  If absent, the
                                                      !! file is opened and closed within this routine.

  logical :: do_read, do_broadcast
  integer :: rc, ncid, varid, att_type, att_len, info(2)
  character(len=256) :: hdr, att_str
  character(len=:), dimension(:), allocatable :: tmp_str
  hdr = "read_attribute_str"
  att_len = 0

  do_read = is_root_pe() ; if (present(all_read)) do_read = all_read .or. do_read
  do_broadcast = .true. ; if (present(all_read)) do_broadcast = .not.all_read

  if (do_read) then
    if (present(ncid_in)) then
      ncid = ncid_in
    else
      call open_file_to_read(filename, ncid, success=found)
      if (present(found)) then ; if (.not.found) do_read = .false. ; endif
    endif
  endif

  if (do_read) then
    rc = NF90_ENOTATT ; att_len = 0
    if (present(varname)) then  ! Read a variable attribute
      call get_varid(varname, ncid, filename, varid, match_case=.false., found=found)
      att_str = "att "//trim(attname)//" for "//trim(varname)//" from "//trim(filename)
    else   ! Read a global attribute
      varid = NF90_GLOBAL
      att_str = "global att "//trim(attname)//" from "//trim(filename)
    endif
    if ((varid > 0) .or. (varid == NF90_GLOBAL)) then ! The named variable does exist, and found would be true.
      rc = NF90_inquire_attribute(ncid, varid, attname, xtype=att_type, len=att_len)
      if ((rc /= NF90_NOERR) .and. (rc /= NF90_ENOTATT)) &
        call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //" Error getting info for "//trim(att_str))
      if (att_type /= NF90_CHAR) &
        call MOM_error(FATAL, trim(hdr)//": Attribute data type is not a char for "//trim(att_str))
!      if (att_len > len(att_val)) &
!        call MOM_error(FATAL, trim(hdr)//": Insufficiently long string passed in to read "//trim(att_str))
      allocate(character(att_len) :: att_val)

      if (rc == NF90_NOERR) then
        rc = NF90_get_att(ncid, varid, attname, att_val)
        if ((rc /= NF90_NOERR) .and. (rc /= NF90_ENOTATT)) &
          call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //" Difficulties reading "//trim(att_str))
      endif
    endif
    if (present(found)) found = (rc == NF90_NOERR)

    if (.not.present(ncid_in)) call close_file_to_read(ncid, filename)
  endif

  if (do_broadcast) then
    ! Communicate the string length
    info(1) = att_len ; info(2) = 0 ; if (do_read .and. found) info(2) = 1
    call broadcast(info, 2, blocking=.true.)
    att_len = info(1)

    if (att_len > 0) then
      ! These extra copies are here because broadcast only supports arrays of strings.
      allocate(character(att_len) :: tmp_str(1))
      if (.not.do_read) allocate(character(att_len) :: att_val)
      if (do_read) tmp_str(1) = att_val
      call broadcast(tmp_str, att_len, blocking=.true.)
      att_val = tmp_str(1)
      if (present(found)) found = (info(2) /= 0)
    elseif (.not.allocated(att_val)) then
      allocate(character(4) :: att_val) ; att_val = ''
    endif
  elseif (.not.allocated(att_val)) then
    allocate(character(4) :: att_val) ; att_val = ''
  endif
end subroutine read_attribute_str


!> Read a 32-bit integer global or variable attribute
subroutine read_attribute_int32(filename, attname, att_val, varname, found, all_read, ncid_in)
  character(len=*),           intent(in)  :: filename !< Name of the file to read
  character(len=*),           intent(in)  :: attname  !< Name of the attribute to read
  integer(kind=int32),        intent(out) :: att_val  !< The value of the attribute
  character(len=*), optional, intent(in)  :: varname  !< The name of the variable whose attribute will
                                                      !! be read. If missing, read a global attribute.
  logical,          optional, intent(out) :: found    !< Returns true if the attribute is found
  logical,          optional, intent(in)  :: all_read !< If present and true, all PEs that call this
                                                      !! routine actually do the read, otherwise only
                                                      !! root PE reads and then broadcasts the results.
  integer,        optional, intent(in)    :: ncid_in  !< The netCDF ID of an open file.  If absent, the
                                                      !! file is opened and closed within this routine.

  logical :: do_read, do_broadcast
  integer :: rc, ncid, varid, is_found
  character(len=256) :: hdr
  hdr = "read_attribute_int32"
  att_val = 0

  do_read = is_root_pe() ; if (present(all_read)) do_read = all_read .or. do_read
  do_broadcast = .true. ; if (present(all_read)) do_broadcast = .not.all_read

  if (do_read) then
    if (present(ncid_in)) then
      ncid = ncid_in
    else
      call open_file_to_read(filename, ncid, success=found)
      if (present(found)) then ; if (.not.found) do_read = .false. ; endif
    endif
  endif

  if (do_read) then
    rc = NF90_ENOTATT
    if (present(varname)) then  ! Read a variable attribute
      call get_varid(varname, ncid, filename, varid, match_case=.false., found=found)
      if (varid >= 0) then ! The named variable does exist, and found would be true.
        rc = NF90_get_att(ncid, varid, attname, att_val)
        if ((rc /= NF90_NOERR) .and. (rc /= NF90_ENOTATT)) &
          call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //" Difficulties reading att "//&
                trim(attname)//" for "//trim(varname)//" from "//trim(filename))
      endif
    else  ! Read a global attribute
      rc = NF90_get_att(ncid, NF90_GLOBAL, attname, att_val)
      if ((rc /= NF90_NOERR) .and. (rc /= NF90_ENOTATT)) &
        call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //&
                " Difficulties reading global att "//trim(attname)//" from "//trim(filename))
    endif
    if (present(found)) found = (rc == NF90_NOERR)

    if (.not.present(ncid_in)) call close_file_to_read(ncid, filename)
  endif

  if (do_broadcast) then
    if (present(found)) then
      is_found = 0 ; if (is_root_pe() .and. found) is_found = 1
      call broadcast(is_found, blocking=.false.)
    endif
    call broadcast(att_val, blocking=.true.)
    if (present(found)) found = (is_found /= 0)
  endif

end subroutine read_attribute_int32


!> Read a 64-bit integer global or variable attribute
subroutine read_attribute_int64(filename, attname, att_val, varname, found, all_read, ncid_in)
  character(len=*),           intent(in)  :: filename !< Name of the file to read
  character(len=*),           intent(in)  :: attname  !< Name of the attribute to read
  integer(kind=int64),        intent(out) :: att_val  !< The value of the attribute
  character(len=*), optional, intent(in)  :: varname  !< The name of the variable whose attribute will
                                                      !! be read. If missing, read a global attribute.
  logical,          optional, intent(out) :: found    !< Returns true if the attribute is found
  logical,          optional, intent(in)  :: all_read !< If present and true, all PEs that call this
                                                      !! routine actually do the read, otherwise only
                                                      !! root PE reads and then broadcasts the results.
  integer,        optional, intent(in)    :: ncid_in  !< The netCDF ID of an open file.  If absent, the
                                                      !! file is opened and closed within this routine.

  logical :: do_read, do_broadcast
  integer :: rc, ncid, varid, is_found
  character(len=256) :: hdr
  hdr = "read_attribute_int64"
  att_val = 0

  do_read = is_root_pe() ; if (present(all_read)) do_read = all_read .or. do_read
  do_broadcast = .true. ; if (present(all_read)) do_broadcast = .not.all_read

  if (do_read) then
    if (present(ncid_in)) then
      ncid = ncid_in
    else
      call open_file_to_read(filename, ncid, success=found)
      if (present(found)) then ; if (.not.found) do_read = .false. ; endif
    endif
  endif

  if (do_read) then
    rc = NF90_ENOTATT
    if (present(varname)) then  ! Read a variable attribute
      call get_varid(varname, ncid, filename, varid, match_case=.false., found=found)
      if (varid >= 0) then ! The named variable does exist, and found would be true.
        rc = NF90_get_att(ncid, varid, attname, att_val)
        if ((rc /= NF90_NOERR) .and. (rc /= NF90_ENOTATT)) &
          call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //" Difficulties reading att "//&
                trim(attname)//" for "//trim(varname)//" from "//trim(filename))
      endif
    else  ! Read a global attribute
      rc = NF90_get_att(ncid, NF90_GLOBAL, attname, att_val)
      if ((rc /= NF90_NOERR) .and. (rc /= NF90_ENOTATT)) &
        call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //&
                " Difficulties reading global att "//trim(attname)//" from "//trim(filename))
    endif
    if (present(found)) found = (rc == NF90_NOERR)

    rc = NF90_close(ncid)
  endif

  if (do_broadcast) then
    if (present(found)) then
      is_found = 0 ; if (is_root_pe() .and. found) is_found = 1
      call broadcast(is_found, blocking=.false.)
    endif
    call broadcast(att_val, blocking=.true.)
    if (present(found)) found = (is_found /= 0)
  endif

end subroutine read_attribute_int64

!> Read a real global or variable attribute
subroutine read_attribute_real(filename, attname, att_val, varname, found, all_read, ncid_in)
  character(len=*),           intent(in)  :: filename !< Name of the file to read
  character(len=*),           intent(in)  :: attname  !< Name of the attribute to read
  real,                       intent(out) :: att_val  !< The value of the attribute
  character(len=*), optional, intent(in)  :: varname  !< The name of the variable whose attribute will
                                                      !! be read. If missing, read a global attribute.
  logical,          optional, intent(out) :: found    !< Returns true if the attribute is found
  logical,          optional, intent(in)  :: all_read !< If present and true, all PEs that call this
                                                      !! routine actually do the read, otherwise only
                                                      !! root PE reads and then broadcasts the results.
  integer,          optional, intent(in)  :: ncid_in  !< The netCDF ID of an open file.  If absent, the
                                                      !! file is opened and closed within this routine.

  logical :: do_read, do_broadcast
  integer :: rc, ncid, varid, is_found
  character(len=256) :: hdr
  hdr = "read_attribute_real"
  att_val = 0.0

  do_read = is_root_pe() ; if (present(all_read)) do_read = all_read .or. do_read
  do_broadcast = .true. ; if (present(all_read)) do_broadcast = .not.all_read

  if (do_read) then
    if (present(ncid_in)) then
      ncid = ncid_in
    else
      call open_file_to_read(filename, ncid, success=found)
      if (present(found)) then ; if (.not.found) do_read = .false. ; endif
    endif
  endif

  if (do_read) then
    rc = NF90_ENOTATT
    if (present(varname)) then  ! Read a variable attribute
      call get_varid(varname, ncid, filename, varid, match_case=.false., found=found)
      if (varid >= 0) then ! The named variable does exist, and found would be true.
        rc = NF90_get_att(ncid, varid, attname, att_val)
        if ((rc /= NF90_NOERR) .and. (rc /= NF90_ENOTATT)) &
          call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //" Difficulties reading att "//&
                trim(attname)//" for "//trim(varname)//" from "//trim(filename))
      endif
    else  ! Read a global attribute
      rc = NF90_get_att(ncid, NF90_GLOBAL, attname, att_val)
      if ((rc /= NF90_NOERR) .and. (rc /= NF90_ENOTATT)) &
        call MOM_error(FATAL, trim(hdr) // trim(NF90_STRERROR(rc)) //&
                " Difficulties reading global att "//trim(attname)//" from "//trim(filename))
    endif
    if (present(found)) found = (rc == NF90_NOERR)

    if (.not.present(ncid_in)) call close_file_to_read(ncid, filename)
  endif

  if (do_broadcast) then
    if (present(found)) then
      is_found = 0 ; if (is_root_pe() .and. found) is_found = 1
      call broadcast(is_found, blocking=.false.)
    endif
    call broadcast(att_val, blocking=.true.)
    if (present(found)) found = (is_found /= 0)
  endif

end subroutine read_attribute_real

!> Open a netcdf file for reading, with error handling
subroutine open_file_to_read(filename, ncid, success)
  character(len=*),  intent(in)  :: filename   !< path and name of the file to open for reading
  integer,           intent(out) :: ncid       !< The netcdf handle for the file
  logical, optional, intent(out) :: success    !< Returns true if the file was opened, or if this
                                               !! argument is not present, failure is fatal error.
  ! Local variables
  integer rc

  rc = NF90_open(trim(filename), NF90_NOWRITE, ncid)
  if (present(success)) then
    success = (rc == NF90_NOERR)
  elseif (rc /= NF90_NOERR) then
    call MOM_error(FATAL, "Difficulties opening "//trim(filename)//" - "//trim(NF90_STRERROR(rc)) )
  endif

end subroutine open_file_to_read

!> Close a netcdf file that had been opened for reading, with error handling
subroutine close_file_to_read(ncid, filename)
  integer,                    intent(inout) :: ncid       !< The netcdf handle for the file to close
  character(len=*), optional, intent(in)    :: filename   !< path and name of the file to close
  integer :: rc
  if (ncid >= 0) then
    rc = NF90_close(ncid)
    if (present(filename) .and. (rc /= NF90_NOERR)) then
      call MOM_error(WARNING, "Difficulties closing "//trim(filename)//": "//trim(NF90_STRERROR(rc)))
    elseif (rc /= NF90_NOERR) then
      call MOM_error(WARNING, "Difficulties closing file: "//trim(NF90_STRERROR(rc)))
    endif
  endif
  ncid = -1
end subroutine close_file_to_read

!> get_varid finds the netcdf handle for the potentially case-insensitive variable name in a file
subroutine get_varid(varname, ncid, filename, varid, match_case, found)
  character(len=*),  intent(in)  :: varname    !< The name of the variable that is being sought
  integer,           intent(in)  :: ncid       !< The open netcdf handle for the file
  character(len=*),  intent(in)  :: filename   !< name of the file to read, used here in messages
  integer,           intent(out) :: varid      !< The netcdf handle for the variable
  logical, optional, intent(in)  :: match_case !< If false, allow for variables name matches to be
                                               !! case insensitive, but take a perfect match if
                                               !! found.  The default is true.
  logical, optional, intent(out) :: found      !< Returns true if the attribute is found

  logical :: var_found, insensitive
  character(len=256) :: name
  integer, allocatable :: varids(:)
  integer :: nvars, status, n

  varid = -1
  var_found = .false.
  insensitive = .false. ; if (present(match_case)) insensitive = .not.match_case

  if (insensitive) then
    ! This code ounddoes a case-insensitive search for a variable in the file.
    status = NF90_inquire(ncid, nVariables=nvars)
    if (present(found) .and. ((status /= NF90_NOERR) .or. (nvars < 1))) then
      found = .false. ; return
    elseif (status /= NF90_NOERR) then
      call MOM_error(FATAL, "get_varid:  Difficulties getting the number of variables in file "//&
          trim(filename)//" - "//trim(NF90_STRERROR(status)))
    elseif (nvars < 1) then
      call MOM_error(FATAL, "get_varid: There appear not to be any variables in "//trim(filename))
    endif

    allocate(varids(nvars))

    status = nf90_inq_varids(ncid, nvars, varids)
    if (status /= NF90_NOERR) then
      call MOM_error(WARNING, "get_varid: Difficulties getting the variable IDs in file "//&
          trim(filename)//" - "//trim(NF90_STRERROR(status)))
      nvars = -1  ! Full error handling will occur after the do-loop.
    endif

    do n = 1,nvars
      status = nf90_inquire_variable(ncid, varids(n), name=name)
      if (status /= NF90_NOERR) then
        call MOM_error(WARNING, "get_varid:  Difficulties getting a variable name in file "//&
            trim(filename)//" - "//trim(NF90_STRERROR(status)))
      endif

      if (trim(lowercase(name)) == trim(lowercase(varname))) then
        if (var_found) then
          call MOM_error(WARNING, "get_varid: Two variables match the case-insensitive name "//&
                  trim(varname)//" in file "//trim(filename))
          ! Replace the first variable if the second one is a case-sensitive match
          if (trim(name) == trim(varname)) varid = varids(n)
        else
          varid = varids(n) ; var_found = .true.
        endif
      endif
    enddo
    if (present(found)) found = var_found
    if ((.not.var_found) .and. .not.present(found)) call MOM_error(FATAL, &
        "get_varid: variable "//trim(varname)//" was not found in file "//trim(filename))

    deallocate(varids)
  else
    status = NF90_INQ_VARID(ncid, trim(varname), varid)
    if (present(found)) found = (status == NF90_NOERR)
    if ((status /= NF90_NOERR) .and. .not.present(found)) then
      call MOM_error(FATAL, "get_varid: Difficulties getting a variable id for "//&
          trim(varname)//" in file "//trim(filename)//" - "//trim(NF90_STRERROR(status)))
    endif
  endif

end subroutine get_varid

!> Verify that a file contains a named variable with the expected units.
subroutine verify_variable_units(filename, varname, expected_units, msg, ierr, alt_units)
  character(len=*),           intent(in)    :: filename  !< File name
  character(len=*),           intent(in)    :: varname   !< Variable name
  character(len=*),           intent(in)    :: expected_units !< Expected units of variable
  character(len=*),           intent(inout) :: msg       !< Message to use for errors
  logical,                    intent(out)   :: ierr      !< True if an error occurs
  character(len=*), optional, intent(in)    :: alt_units !< Alterate acceptable units of variable

  ! Local variables
  character (len=200) :: units
  logical :: units_correct, success
  integer :: i, ncid, status, vid

  if (.not.is_root_pe()) then ! Only the root PE should do the verification.
    ierr = .false. ; msg = '' ; return
  endif

  ierr = .true.
  call open_file_to_read(filename, ncid, success)
  if (.not.success) then
    msg = 'File not found: '//trim(filename)
    return
  endif

  status = NF90_INQ_VARID(ncid, trim(varname), vid)
  if (status /= NF90_NOERR) then
    msg = 'Var not found: '//trim(varname)
  else
    status = NF90_GET_ATT(ncid, vid, "units", units)
    if (status /= NF90_NOERR) then
      msg = 'Attribute not found: units'
    else
      ! NF90_GET_ATT can return attributes with null characters, which TRIM will not truncate.
      ! This loop replaces any null characters with a space so that the subsequent check
      ! between the read units and the expected units will pass
      do i=1,LEN_TRIM(units)
        if (units(i:i) == CHAR(0)) units(i:i) = " "
      enddo

      units_correct = (trim(units) == trim(expected_units))
      if (present(alt_units)) then
        units_correct = units_correct .or. (trim(units) == trim(alt_units))
      endif
      if (units_correct) then
        ierr = .false.
        msg = ''
      else
        msg = 'Units incorrect: '//trim(units)//' /= '//trim(expected_units)
      endif
    endif
  endif

  status = NF90_close(ncid)

end subroutine verify_variable_units

!> Returns a vardesc type whose elements have been filled with the provided
!! fields.  The argument name is required, while the others are optional and
!! have default values that are empty strings or are appropriate for a 3-d
!! tracer field at the tracer cell centers.
function var_desc(name, units, longname, hor_grid, z_grid, t_grid, cmor_field_name, &
                  cmor_units, cmor_longname, conversion, caller, position, dim_names, fixed) result(vd)
  character(len=*),           intent(in) :: name            !< variable name
  character(len=*), optional, intent(in) :: units           !< variable units
  character(len=*), optional, intent(in) :: longname        !< variable long name
  character(len=*), optional, intent(in) :: hor_grid        !< A character string indicating the horizontal
                                                            !! position of this variable
  character(len=*), optional, intent(in) :: z_grid          !< variable vertical staggering
  character(len=*), optional, intent(in) :: t_grid          !< time description: s, p, or 1
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name
  character(len=*), optional, intent(in) :: cmor_units      !< CMOR physical dimensions of variable
  character(len=*), optional, intent(in) :: cmor_longname   !< CMOR long name
  real            , optional, intent(in) :: conversion      !< for unit conversions, such as needed to
                                                            !! convert from intensive to extensive
  character(len=*), optional, intent(in) :: caller          !< The calling routine for error messages
  integer,          optional, intent(in) :: position        !< A coded integer indicating the horizontal position
                                                            !! of this variable if it has such dimensions.
                                                            !! Valid values include CORNER, CENTER, EAST_FACE
                                                            !! NORTH_FACE, and 0 for no horizontal dimensions.
  character(len=*), dimension(:), &
                    optional, intent(in) :: dim_names       !< The names of the dimensions of this variable
  logical,          optional, intent(in) :: fixed           !< If true, this does not evolve with time
  type(vardesc)                          :: vd              !< vardesc type that is created

  character(len=120) :: cllr
  cllr = "var_desc"
  if (present(caller)) cllr = trim(caller)

  call safe_string_copy(name, vd%name, "vd%name", cllr)

  vd%longname = "" ; vd%units = ""
  vd%hor_grid = 'h' ; vd%position = CENTER ; vd%z_grid = 'L' ; vd%t_grid = 's'
  if (present(dim_names)) vd%z_grid = '1' ! In this case the names are used to set the non-horizontal axes
  if (present(fixed)) then ; if (fixed) vd%t_grid = '1' ; endif

  vd%cmor_field_name  =  ""
  vd%cmor_units       =  ""
  vd%cmor_longname    =  ""
  vd%conversion       =  1.0
  vd%dim_names(:)     =  ""

  call modify_vardesc(vd, units=units, longname=longname, hor_grid=hor_grid, &
                      z_grid=z_grid, t_grid=t_grid, position=position, dim_names=dim_names, &
                      cmor_field_name=cmor_field_name, cmor_units=cmor_units, &
                      cmor_longname=cmor_longname, conversion=conversion, caller=cllr)

end function var_desc


!> This routine modifies the named elements of a vardesc type.
!! All arguments are optional, except the vardesc type to be modified.
subroutine modify_vardesc(vd, name, units, longname, hor_grid, z_grid, t_grid, &
                 cmor_field_name, cmor_units, cmor_longname, conversion, caller, position, dim_names)
  type(vardesc),              intent(inout) :: vd              !< vardesc type that is modified
  character(len=*), optional, intent(in)    :: name            !< name of variable
  character(len=*), optional, intent(in)    :: units           !< units of variable
  character(len=*), optional, intent(in)    :: longname        !< long name of variable
  character(len=*), optional, intent(in)    :: hor_grid        !< horizontal staggering of variable
  character(len=*), optional, intent(in)    :: z_grid          !< vertical staggering of variable
  character(len=*), optional, intent(in)    :: t_grid          !< time description: s, p, or 1
  character(len=*), optional, intent(in)    :: cmor_field_name !< CMOR name
  character(len=*), optional, intent(in)    :: cmor_units      !< CMOR physical dimensions of variable
  character(len=*), optional, intent(in)    :: cmor_longname   !< CMOR long name
  real            , optional, intent(in)    :: conversion      !< A multiplicative factor for unit conversions,
                                                               !! such as needed to convert from intensive to
                                                               !! extensive or dimensional consistency testing
  character(len=*), optional, intent(in)    :: caller          !< The calling routine for error messages
  integer,          optional, intent(in)    :: position        !< A coded integer indicating the horizontal position
                                                               !! of this variable if it has such dimensions.
                                                               !! Valid values include CORNER, CENTER, EAST_FACE
                                                               !! NORTH_FACE, and 0 for no horizontal dimensions.
  character(len=*), dimension(:), &
                    optional, intent(in)    :: dim_names       !< The names of the dimensions of this variable

  character(len=120) :: cllr
  integer :: n

  cllr = "mod_vardesc" ; if (present(caller)) cllr = trim(caller)

  if (present(name))      call safe_string_copy(name, vd%name, "vd%name", cllr)

  if (present(longname))  call safe_string_copy(longname, vd%longname, &
                               "vd%longname of "//trim(vd%name), cllr)
  if (present(units))     call safe_string_copy(units, vd%units,       &
                               "vd%units of "//trim(vd%name), cllr)
  if (present(position)) then
    vd%position = position
    select case (position)
      case (CENTER)     ; vd%hor_grid = 'T'
      case (CORNER)     ; vd%hor_grid = 'Bu'
      case (EAST_FACE)  ; vd%hor_grid = 'Cu'
      case (NORTH_FACE) ; vd%hor_grid = 'Cv'
      case (0)          ; vd%hor_grid = '1'
      case default
        call MOM_error(FATAL, "modify_vardesc: "//trim(vd%name)//" has unrecognized position argument")
    end select
  endif
  if (present(hor_grid)) then
    call safe_string_copy(hor_grid, vd%hor_grid, "vd%hor_grid of "//trim(vd%name), cllr)
    vd%position = position_from_horgrid(vd%hor_grid)
    if (present(caller) .and. (vd%position == -1)) then
      call MOM_error(FATAL, "modify_vardesc called by "//trim(caller)//": "//trim(vd%name)//&
                   " has an unrecognized hor_grid argument "//trim(vd%hor_grid))
    elseif (vd%position == -1) then
      call MOM_error(FATAL, "modify_vardesc called with bad hor_grid argument "//trim(vd%hor_grid))
    endif
  endif
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

  if (present(dim_names)) then
    do n=1,min(5,size(dim_names)) ; if (len_trim(dim_names(n)) > 0) then
      call safe_string_copy(dim_names(n), vd%dim_names(n), "vd%dim_names of "//trim(vd%name), cllr)
    endif ; enddo
  endif

end subroutine modify_vardesc

integer function position_from_horgrid(hor_grid)
  character(len=*), intent(in)    :: hor_grid        !< horizontal staggering of variable

  select case (trim(hor_grid))
    case ('h')   ; position_from_horgrid = CENTER
    case ('q')   ; position_from_horgrid = CORNER
    case ('u')   ; position_from_horgrid = EAST_FACE
    case ('v')   ; position_from_horgrid = NORTH_FACE
    case ('T')   ; position_from_horgrid = CENTER
    case ('Bu')  ; position_from_horgrid = CORNER
    case ('Cu')  ; position_from_horgrid = EAST_FACE
    case ('Cv')  ; position_from_horgrid = NORTH_FACE
    case ('1')   ; position_from_horgrid = 0
    case default ; position_from_horgrid = -1 ! This is a bad-value flag.
  end select
end function position_from_horgrid

!> Store information that can be used to create an axis in a subsequent call to create_file.
subroutine set_axis_info(axis, name, units, longname, ax_size, ax_data, cartesian, sense)
  type(axis_info),              intent(inout) :: axis  !< A type with information about a named axis
  character(len=*),             intent(in)    :: name  !< The name of this axis for use in files
  character(len=*),   optional, intent(in)    :: units !< The units of the axis labels
  character(len=*),   optional, intent(in)    :: longname  !< Long name of the axis variable
  integer,            optional, intent(in)    :: ax_size !< The number of elements in this axis
  real, dimension(:), optional, intent(in)    :: ax_data !< The values of the data on the axis
  character(len=*),   optional, intent(in)    :: cartesian !< A variable indicating which direction this axis
                                                       !! axis corresponds with. Valid values
                                                       !! include 'X', 'Y', 'Z', 'T', and 'N' (the default) for none.
  integer,            optional, intent(in)    :: sense !< This is 1 for axes whose values increase upward, or -1
                                                       !! if they increase downward.  The default, 0, is ignored.

  call safe_string_copy(name, axis%name, "axis%name of "//trim(name), "set_axis_info")
  ! Set the default values.
  axis%longname = trim(axis%name) ; axis%units = "" ; axis%cartesian = "N" ; axis%sense = 0

  if (present(longname))  call safe_string_copy(longname, axis%longname, &
                                                "axis%longname of "//trim(name), "set_axis_info")
  if (present(units))     call safe_string_copy(units, axis%units, &
                                                "axis%units of "//trim(name), "set_axis_info")
  if (present(cartesian)) call safe_string_copy(cartesian, axis%cartesian, &
                                                "axis%cartesian of "//trim(name), "set_axis_info")
  if (present(sense)) axis%sense = sense

  if (.not.(present(ax_size) .or. present(ax_data)) ) then
    call MOM_error(FATAL, "set_axis_info called for "//trim(name)//&
      "without either an ax_size or an ax_data argument.")
  elseif (present(ax_size) .and. present(ax_data)) then
    if (size(ax_data) /= ax_size) call MOM_error(FATAL, "set_axis_info called for "//trim(name)//&
      "with an inconsistent value of ax_size and size of ax_data.")
  endif

  if (present(ax_size)) then
    axis%ax_size = ax_size
  else
    axis%ax_size = size(ax_data)
  endif
  if (present(ax_data)) then
    allocate(axis%ax_data(axis%ax_size)) ; axis%ax_data(:) = ax_data(:)
  endif

end subroutine set_axis_info

!> Delete the information in an array of axis_info types and deallocate memory in them.
subroutine delete_axis_info(axes)
  type(axis_info), dimension(:), intent(inout) :: axes  !< An array with information about named axes

  integer :: n
  do n=1,size(axes)
    axes(n)%name = "" ; axes(n)%longname = "" ; axes(n)%units = "" ; axes(n)%cartesian = "N"
    axes(n)%sense = 0 ; axes(n)%ax_size = 0
    if (allocated(axes(n)%ax_data)) deallocate(axes(n)%ax_data)
  enddo
end subroutine delete_axis_info

!> Store information that can be used to create an attribute in a subsequent call to create_file.
subroutine set_attribute_info(attribute, name, str_value)
  type(attribute_info), intent(inout) :: attribute !< A type with information about a named attribute
  character(len=*),     intent(in)    :: name      !< The name of this attribute for use in files
  character(len=*),     intent(in)    :: str_value !< The value of this attribute

  attribute%name = trim(name)
  attribute%att_val = trim(str_value)
end subroutine set_attribute_info

!> Delete the information in an array of attribute_info types and deallocate memory in them.
subroutine delete_attribute_info(atts)
  type(attribute_info), dimension(:), intent(inout) :: atts  !< An array of global attributes

  integer :: n
  do n=1,size(atts)
    if (allocated(atts(n)%name)) deallocate(atts(n)%name)
    if (allocated(atts(n)%att_val)) deallocate(atts(n)%att_val)
  enddo
end subroutine delete_attribute_info


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
                         cmor_field_name, cmor_units, cmor_longname, conversion, caller, &
                         position, dim_names)
  type(vardesc),              intent(in)  :: vd                 !< vardesc type that is queried
  character(len=*), optional, intent(out) :: name               !< name of variable
  character(len=*), optional, intent(out) :: units              !< units of variable
  character(len=*), optional, intent(out) :: longname           !< long name of variable
  character(len=*), optional, intent(out) :: hor_grid           !< horizontal staggering of variable
  character(len=*), optional, intent(out) :: z_grid             !< verticle staggering of variable
  character(len=*), optional, intent(out) :: t_grid             !< time description: s, p, or 1
  character(len=*), optional, intent(out) :: cmor_field_name    !< CMOR name
  character(len=*), optional, intent(out) :: cmor_units         !< CMOR physical dimensions of variable
  character(len=*), optional, intent(out) :: cmor_longname      !< CMOR long name
  real            , optional, intent(out) :: conversion         !< for unit conversions, such as needed to
                                                                !! convert from intensive to extensive
  character(len=*), optional, intent(in)  :: caller             !< calling routine?
  integer,          optional, intent(out) :: position        !< A coded integer indicating the horizontal position
                                                            !! of this variable if it has such dimensions.
                                                            !! Valid values include CORNER, CENTER, EAST_FACE
                                                            !! NORTH_FACE, and 0 for no horizontal dimensions.
  character(len=*), dimension(:), &
                    optional, intent(out) :: dim_names       !< The names of the dimensions of this variable

  integer :: n
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
  if (present(position)) then
    position = vd%position
    if (position == -1) position = position_from_horgrid(vd%hor_grid)
  endif
  if (present(dim_names)) then
    do n=1,min(5,size(dim_names))
      call safe_string_copy(vd%dim_names(n), dim_names(n), "vd%dim_names of "//trim(vd%name), cllr)
    enddo
  endif

end subroutine query_vardesc


!> Write a 4d field to an output file, potentially with rotation
subroutine MOM_write_field_4d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, &
                              fill_value, turns, scale)
  type(file_type),          intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),          intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),    intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:,:,:), intent(inout) :: field      !< Unrotated field to write
  real,           optional, intent(in)    :: tstamp     !< Model timestamp
  integer,        optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,           optional, intent(in)    :: fill_value !< Missing data fill value
  integer,        optional, intent(in)    :: turns      !< Number of quarter-turns to rotate the data
  real,           optional, intent(in)    :: scale      !< A scaling factor that the field is
                                                        !! multiplied by before it is written

  real, allocatable :: field_rot(:,:,:,:)  ! A rotated version of field, with the same units or rescaled
  real :: scale_fac ! A scaling factor to use before writing the array
  integer :: qturns ! The number of quarter turns through which to rotate field
  integer :: is, ie, js, je ! The extent of the computational domain

  qturns = 0 ; if (present(turns)) qturns = modulo(turns, 4)
  scale_fac = 1.0 ; if (present(scale)) scale_fac = scale

  if ((qturns == 0) .and. (scale_fac == 1.0)) then
    call write_field(IO_handle, field_md, MOM_domain, field, tstamp=tstamp, &
                         tile_count=tile_count, fill_value=fill_value)
  else
    call allocate_rotated_array(field, [1,1,1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    if (scale_fac /= 1.0) call rescale_comp_data(MOM_Domain, field_rot, scale_fac)
    call write_field(IO_handle, field_md, MOM_domain, field_rot, tstamp=tstamp, &
                         tile_count=tile_count, fill_value=fill_value)
    deallocate(field_rot)
  endif
end subroutine MOM_write_field_4d

!> Write a 3d field to an output file, potentially with rotation
subroutine MOM_write_field_3d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, &
                              fill_value, turns, scale)
  type(file_type),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),  intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:,:), intent(inout) :: field      !< Unrotated field to write
  real,         optional, intent(in)    :: tstamp     !< Model timestamp
  integer,      optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,         optional, intent(in)    :: fill_value !< Missing data fill value
  integer,      optional, intent(in)    :: turns      !< Number of quarter-turns to rotate the data
  real,         optional, intent(in)    :: scale      !< A scaling factor that the field is
                                                      !! multiplied by before it is written

  real, allocatable :: field_rot(:,:,:)  ! A rotated version of field, with the same units or rescaled
  real :: scale_fac ! A scaling factor to use before writing the array
  integer :: qturns ! The number of quarter turns through which to rotate field
  integer :: is, ie, js, je ! The extent of the computational domain

  qturns = 0 ; if (present(turns)) qturns = modulo(turns, 4)
  scale_fac = 1.0 ; if (present(scale)) scale_fac = scale

  if ((qturns == 0) .and. (scale_fac == 1.0)) then
    call write_field(IO_handle, field_md, MOM_domain, field, tstamp=tstamp, &
                         tile_count=tile_count, fill_value=fill_value)
  else
    call allocate_rotated_array(field, [1,1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    if (scale_fac /= 1.0) call rescale_comp_data(MOM_Domain, field_rot, scale_fac)
    call write_field(IO_handle, field_md, MOM_domain, field_rot, tstamp=tstamp, &
                         tile_count=tile_count, fill_value=fill_value)
    deallocate(field_rot)
  endif
end subroutine MOM_write_field_3d

!> Write a 2d field to an output file, potentially with rotation
subroutine MOM_write_field_2d(IO_handle, field_md, MOM_domain, field, tstamp, tile_count, &
                              fill_value, turns, scale)
  type(file_type),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  type(MOM_domain_type),  intent(in)    :: MOM_domain !< The MOM_Domain that describes the decomposition
  real, dimension(:,:),   intent(inout) :: field      !< Unrotated field to write
  real,         optional, intent(in)    :: tstamp     !< Model timestamp
  integer,      optional, intent(in)    :: tile_count !< PEs per tile (default: 1)
  real,         optional, intent(in)    :: fill_value !< Missing data fill value
  integer,      optional, intent(in)    :: turns      !< Number of quarter-turns to rotate the data
  real,         optional, intent(in)    :: scale      !< A scaling factor that the field is
                                                      !! multiplied by before it is written

  real, allocatable :: field_rot(:,:)  ! A rotated version of field, with the same units
  real :: scale_fac ! A scaling factor to use before writing the array
  integer :: qturns ! The number of quarter turns through which to rotate field
  integer :: is, ie, js, je ! The extent of the computational domain

  qturns = 0 ; if (present(turns)) qturns = modulo(turns, 4)
  scale_fac = 1.0 ; if (present(scale)) scale_fac = scale

  if ((qturns == 0) .and. (scale_fac == 1.0)) then
    call write_field(IO_handle, field_md, MOM_domain, field, tstamp=tstamp, &
                         tile_count=tile_count, fill_value=fill_value)
  else
    call allocate_rotated_array(field, [1,1], qturns, field_rot)
    call rotate_array(field, qturns, field_rot)
    if (scale_fac /= 1.0) call rescale_comp_data(MOM_Domain, field_rot, scale_fac)
    call write_field(IO_handle, field_md, MOM_domain, field_rot, tstamp=tstamp, &
                         tile_count=tile_count, fill_value=fill_value)
    deallocate(field_rot)
  endif
end subroutine MOM_write_field_2d

!> Write a 1d field to an output file
subroutine MOM_write_field_1d(IO_handle, field_md, field, tstamp, fill_value, scale)
  type(file_type),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  real, dimension(:),     intent(in)    :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model timestamp
  real,         optional, intent(in)    :: fill_value !< Missing data fill value
  real,         optional, intent(in)    :: scale      !< A scaling factor that the field is
                                                      !! multiplied by before it is written

  real, dimension(:), allocatable :: array ! A rescaled copy of field
  real :: scale_fac ! A scaling factor to use before writing the array
  integer :: i

  scale_fac = 1.0 ; if (present(scale)) scale_fac = scale

  if (scale_fac == 1.0) then
    call write_field(IO_handle, field_md, field, tstamp=tstamp)
  else
    allocate(array(size(field)))
    array(:) = scale_fac * field(:)
    if (present(fill_value)) then
      do i=1,size(field) ; if (field(i) == fill_value) array(i) = fill_value ; enddo
    endif
    call write_field(IO_handle, field_md, array, tstamp=tstamp)
    deallocate(array)
  endif
end subroutine MOM_write_field_1d

!> Write a 0d field to an output file
subroutine MOM_write_field_0d(IO_handle, field_md, field, tstamp, fill_value, scale)
  type(file_type),        intent(inout) :: IO_handle  !< Handle for a file that is open for writing
  type(fieldtype),        intent(in)    :: field_md   !< Field type with metadata
  real,                   intent(in)    :: field      !< Field to write
  real,         optional, intent(in)    :: tstamp     !< Model timestamp
  real,         optional, intent(in)    :: fill_value !< Missing data fill value
  real,         optional, intent(in)    :: scale      !< A scaling factor that the field is
                                                      !! multiplied by before it is written
  real :: scaled_val ! A rescaled copy of field

  scaled_val = field
  if (present(scale)) scaled_val = scale*field
  if (present(fill_value)) then ; if (field == fill_value) scaled_val = fill_value ; endif

  call write_field(IO_handle, field_md, scaled_val, tstamp=tstamp)
end subroutine MOM_write_field_0d

!> Given filename and fieldname, this subroutine returns the size of the field in the file
subroutine field_size(filename, fieldname, sizes, field_found, no_domain, ndims, ncid_in)
  character(len=*),      intent(in)    :: filename  !< The name of the file to read
  character(len=*),      intent(in)    :: fieldname !< The name of the variable whose sizes are returned
  integer, dimension(:), intent(inout) :: sizes     !< The sizes of the variable in each dimension
  logical,     optional, intent(out)   :: field_found !< This indicates whether the field was found in
                                                    !! the input file.  Without this argument, there
                                                    !! is a fatal error if the field is not found.
  logical,     optional, intent(in)    :: no_domain !< If present and true, do not check for file
                                                    !! names with an appended tile number.  If
                                                    !! ndims is present, the default changes to true.
  integer,     optional, intent(out)   :: ndims     !< The number of dimensions to the variable
  integer,     optional, intent(in)    :: ncid_in   !< The netCDF ID of an open file.  If absent, the
                                                    !! file is opened and closed within this routine.

  if (present(ndims)) then
    if (present(no_domain)) then ; if (.not.no_domain) call MOM_error(FATAL, &
          "field_size does not support the ndims argument when no_domain is present and false.")
    endif
    call get_var_sizes(filename, fieldname, ndims, sizes, match_case=.false., ncid_in=ncid_in)
    if (present(field_found)) field_found = (ndims >= 0)
    if ((ndims < 0) .and. .not.present(field_found)) then
      call MOM_error(FATAL, "Variable "//trim(fieldname)//" not found in "//trim(filename) )
    endif
  else
    call get_field_size(filename, fieldname, sizes, field_found=field_found, no_domain=no_domain)
  endif

end subroutine field_size


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

!> Provide a string to append to filenames, to differentiate ensemble members, for example.
subroutine get_filename_appendix(suffix)
  character(len=*), intent(out) :: suffix !< A string to append to filenames

  call get_filename_suffix(suffix)
end subroutine get_filename_appendix

!> Write a file version number to the log file or other output file
subroutine write_version_number(version, tag, unit)
  character(len=*),           intent(in) :: version !< A string that contains the routine name and version
  character(len=*), optional, intent(in) :: tag  !< A tag name to add to the message
  integer,          optional, intent(in) :: unit !< An alternate unit number for output

  call write_version(version, tag, unit)
end subroutine write_version_number


!> Open a single namelist file that is potentially readable by all PEs.
function open_namelist_file(file) result(unit)
  character(len=*), optional, intent(in) :: file !< The file to open, by default "input.nml"
  integer                                :: unit !< The opened unit number of the namelist file
  unit = MOM_namelist_file(file)
end function open_namelist_file

!> Checks the iostat argument that is returned after reading a namelist variable and writes a
!! message if there is an error.
function check_nml_error(IOstat, nml_name) result(ierr)
  integer,          intent(in) :: IOstat   !< An I/O status field from a namelist read call
  character(len=*), intent(in) :: nml_name !< The name of the namelist
  integer :: ierr    !< A copy of IOstat that is returned to preserve legacy function behavior
  call check_namelist_error(IOstat, nml_name)
  ierr = IOstat
end function check_nml_error

!> Initialize the MOM_io module
subroutine MOM_io_init(param_file)
  type(param_file_type), intent(in) :: param_file  !< structure indicating the open file to
                                                   !! parse for model parameter values.

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_io" ! This module's name.

  call log_version(param_file, mdl, version)

end subroutine MOM_io_init

!> \namespace mom_io
!!
!!   This file contains a number of subroutines that manipulate
!!  NetCDF files and handle input and output of fields.  These
!!  subroutines, along with their purpose, are:
!!
!!   * create_file: create a new file and set up structures that are
!!       needed for subsequent output and write out the coordinates.
!!   * reopen_file: reopen an existing file for writing and set up
!!       structures that are needed for subsequent output.
!!   * open_input_file: open the indicated file for reading only.
!!   * close_file: close an open file.
!!   * synch_file: flush the buffers, completing all pending output.
!!
!!   * write_field: write a field to an open file.
!!   * write_time: write a value of the time axis to an open file.
!!   * read_data: read a variable from an open file.
!!   * read_time: read a time from an open file.
!!
!!   * name_output_file: provide a name for an output file based on a
!!       name root and the time of the output.
!!   * find_input_file: find a file that has been previously written by
!!       MOM and named by name_output_file and open it for reading.
!!
!!   * handle_error: write an error code and quit.

end module MOM_io
