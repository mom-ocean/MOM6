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
use fms_io_mod,           only : file_exist, field_size, old_fms_read_data => read_data
use fms_io_mod,           only : field_exists => field_exist, io_infra_end=>fms_io_exit
use fms_io_mod,           only : get_filename_appendix => get_filename_appendix ! FYI: this function only trims strings if used without calling set_filename_appendix
use MOM_string_functions  only : extract_word
use mpp_mod,              only : mpp_max 
use mpp_domains_mod,      only : domain1d, domain2d, domainug, mpp_get_domain_components
use mpp_domains_mod,      only : CENTER, CORNER, NORTH_FACE=>NORTH, EAST_FACE=>EAST
use mpp_io_mod,           only : mpp_open_file => mpp_open, mpp_close_file => mpp_close
use mpp_io_mod,           only : mpp_write_meta, write_field => mpp_write, mpp_get_info
use mpp_io_mod,           only : mpp_get_atts, mpp_get_axes, get_axis_data=>mpp_get_axis_data, axistype
use mpp_io_mod,           only : mpp_get_fields, fieldtype, axistype, flush_file => mpp_flush
use mpp_io_mod,           only : APPEND_FILE=>MPP_APPEND, ASCII_FILE=>MPP_ASCII
use mpp_io_mod,           only : MULTIPLE=>MPP_MULTI, NETCDF_FILE=>MPP_NETCDF
use mpp_io_mod,           only : OVERWRITE_FILE=>MPP_OVERWR, READONLY_FILE=>MPP_RDONLY
use mpp_io_mod,           only : SINGLE_FILE=>MPP_SINGLE, WRITEONLY_FILE=>MPP_WRONLY
use mpp_io_mod,           only : MPP_APPEND, MPP_MULTI, MPP_OVERWR, MPP_NETCDF, MPP_RDONLY
use mpp_io_mod,           only : get_file_info=>mpp_get_info, get_file_atts=>mpp_get_atts
use mpp_io_mod,           only : get_file_fields=>mpp_get_fields, get_file_times=>mpp_get_times
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

public :: mpp_close_file, mpp_open_file, create_file, field_exists, field_size, fieldtype, get_filename_appendix
public :: flush_file, get_file_info, get_file_atts, get_file_fields
public :: get_file_times, read_axis_data
public :: num_timelevels, MOM_read_data, MOM_read_vector, ensembler
public :: reopen_file, slasher, write_field, write_version_number, MOM_io_init
public :: open_namelist_file, check_nml_error, io_infra_init, io_infra_end
public :: APPEND_FILE, ASCII_FILE, MULTIPLE, NETCDF_FILE, OVERWRITE_FILE
public :: READONLY_FILE, SINGLE_FILE, WRITEONLY_FILE
public :: CENTER, CORNER, NORTH_FACE, EAST_FACE
public :: var_desc, modify_vardesc, query_vardesc, cmor_long_std
public :; scale_data
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
public :: get_horizontal_grid_position
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
public :: MOM_get_axis_data
public :: MOM_open_file
public :: MOM_register_diagnostic_axis
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
  character(len=8)   :: axis                    !< Name of the cartesian axis: X,Y,Z,T
  character(len=8)   :: positive                !< Positive-definite direction: 
                                                !! up, down, east, west, north, south
  integer            :: horgrid_position        !< Horizontal grid position
  integer            :: x_position              !< x-direction grid position
  integer            :: y_position              !< y-direction grid position
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

!> Open a netCDF file 
interface MOM_open_file
  module procedure MOM_open_file_DD_ocean_grid
  module procedure MOM_open_file_DD_supergrid
  module procedure MOM_open_file_DD_dyn_horgrid
  module procedure MOM_open_file_noDD
end interface
!> Register axes to a netCDF file
interface MOM_register_axis
  module procedure MOM_register_axis_DD
  module procedure MOM_register_axis_noDD
end interface 

!> Read a data field from a file
interface MOM_read_data
  module procedure MOM_read_data_4d
  module procedure MOM_read_data_3d
  module procedure MOM_read_data_2d
  module procedure MOM_read_data_1d
end interface

!> Read a pair of data fields representing the two components of a vector from a file
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

!> Routine creates a new NetCDF file.  It also sets up
!! structures that describe this file and variables that will
!! later be written to this file. Type for describing a variable, typically a tracer
subroutine create_file(unit, filename, vars, novars, fields, threading, timeunit, G, dG, GV, checksums)
  integer,               intent(out)   :: unit       !< unit id of an open file or -1 on a
                                                     !! nonwriting PE with single file output
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
  integer(kind=8), optional,      intent(in)    :: checksums(:,:)  !< checksums of vars

  logical        :: use_lath, use_lonh, use_latq, use_lonq, use_time
  logical        :: use_layer, use_int, use_periodic
  logical        :: one_file, domain_set
  type(axistype) :: axis_lath, axis_latq, axis_lonh, axis_lonq
  type(axistype) :: axis_layer, axis_int, axis_time, axis_periodic
  type(axistype) :: axes(4)
  type(MOM_domain_type), pointer :: Domain => NULL()
  type(domain1d) :: x_domain, y_domain
  integer        :: numaxes, pack, thread, k
  integer        :: isg, ieg, jsg, jeg, IsgB, IegB, JsgB, JegB
  integer        :: var_periods, num_periods=0
  real, dimension(:), allocatable :: period_val
  real, pointer, dimension(:) :: &
    gridLatT => NULL(), & ! The latitude or longitude of T or B points for
    gridLatB => NULL(), & ! the purpose of labeling the output axes.
    gridLonT => NULL(), gridLonB => NULL()
  character(len=40) :: time_units, x_axis_units, y_axis_units
  character(len=8)  :: t_grid, t_grid_read

  use_lath  = .false. ; use_lonh     = .false.
  use_latq  = .false. ; use_lonq     = .false.
  use_time  = .false. ; use_periodic = .false.
  use_layer = .false. ; use_int      = .false.

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
    call mpp_open_file(unit, filename, MPP_OVERWR, MPP_NETCDF, threading=thread)
  else
    call mpp_open_file(unit, filename, MPP_OVERWR, MPP_NETCDF, domain=Domain%mpp_domain)
  endif

! Define the coordinates.
  do k=1,novars
    select case (vars(k)%hor_grid)
      case ('h') ; use_lath = .true. ; use_lonh = .true.
      case ('q') ; use_latq = .true. ; use_lonq = .true.
      case ('u') ; use_lath = .true. ; use_lonq = .true.
      case ('v') ; use_latq = .true. ; use_lonh = .true.
      case ('T')  ; use_lath = .true. ; use_lonh = .true.
      case ('Bu') ; use_latq = .true. ; use_lonq = .true.
      case ('Cu') ; use_lath = .true. ; use_lonq = .true.
      case ('Cv') ; use_latq = .true. ; use_lonh = .true.
      case ('1') ! Do nothing.
      case default
        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized hor_grid "//trim(vars(k)%hor_grid))
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
  enddo

  if ((use_lath .or. use_lonh .or. use_latq .or. use_lonq)) then
    if (.not.domain_set) call MOM_error(FATAL, "create_file: "//&
      "An ocean_grid_type or dyn_horgrid_type is required to create a file with a horizontal coordinate.")

    call mpp_get_domain_components(Domain%mpp_domain, x_domain, y_domain)
  endif
  if ((use_layer .or. use_int) .and. .not.present(GV)) call MOM_error(FATAL, &
    "create_file: A vertical grid type is required to create a file with a vertical coordinate.")

! Specify all optional arguments to mpp_write_meta: name, units, longname, cartesian, calendar, sense,
! domain, data, min). Otherwise if optional arguments are added to mpp_write_meta the compiler may
! (and in case of GNU does) get confused and crash.
  if (use_lath) &
    call mpp_write_meta(unit, axis_lath, name="lath", units=y_axis_units, longname="Latitude", &
                   cartesian='Y', domain = y_domain, data=gridLatT(jsg:jeg))

  if (use_lonh) &
    call mpp_write_meta(unit, axis_lonh, name="lonh", units=x_axis_units, longname="Longitude", &
                   cartesian='X', domain = x_domain, data=gridLonT(isg:ieg))

  if (use_latq) &
    call mpp_write_meta(unit, axis_latq, name="latq", units=y_axis_units, longname="Latitude", &
                   cartesian='Y', domain = y_domain, data=gridLatB(JsgB:JegB))

  if (use_lonq) &
    call mpp_write_meta(unit, axis_lonq, name="lonq", units=x_axis_units, longname="Longitude", &
                   cartesian='X', domain = x_domain, data=gridLonB(IsgB:IegB))

  if (use_layer) &
    call mpp_write_meta(unit, axis_layer, name="Layer", units=trim(GV%zAxisUnits), &
          longname="Layer "//trim(GV%zAxisLongName), cartesian='Z', &
          sense=1, data=GV%sLayer(1:GV%ke))

  if (use_int) &
    call mpp_write_meta(unit, axis_int, name="Interface", units=trim(GV%zAxisUnits), &
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

    call mpp_write_meta(unit, axis_time, name="Time", units=time_units, longname="Time", cartesian='T')
  else
    call mpp_write_meta(unit, axis_time, name="Time", units="days", longname="Time",cartesian= 'T')
  endif ; endif

  if (use_periodic) then
    if (num_periods <= 1) call MOM_error(FATAL, "MOM_io create_file: "//&
      "num_periods for file "//trim(filename)//" must be at least 1.")
    ! Define a periodic axis with unit labels.
    allocate(period_val(num_periods))
    do k=1,num_periods ; period_val(k) = real(k) ; enddo
    call mpp_write_meta(unit, axis_periodic, name="Period", units="nondimensional", &
          longname="Periods for cyclical varaiables", cartesian= 't', data=period_val)
    deallocate(period_val)
  endif

  do k=1,novars
    numaxes = 0
    select case (vars(k)%hor_grid)
      case ('h')  ; numaxes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_lath
      case ('q')  ; numaxes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_latq
      case ('u')  ; numaxes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_lath
      case ('v')  ; numaxes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_latq
      case ('T')  ; numaxes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_lath
      case ('Bu') ; numaxes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_latq
      case ('Cu') ; numaxes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_lath
      case ('Cv') ; numaxes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_latq
      case ('1') ! Do nothing.
      case default
        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized hor_grid "//trim(vars(k)%hor_grid))
    end select
    select case (vars(k)%z_grid)
      case ('L') ; numaxes = numaxes+1 ; axes(numaxes) = axis_layer
      case ('i') ; numaxes = numaxes+1 ; axes(numaxes) = axis_int
      case ('1') ! Do nothing.
      case default
        call MOM_error(FATAL, "MOM_io create_file: "//trim(vars(k)%name)//&
                        " has unrecognized z_grid "//trim(vars(k)%z_grid))
    end select
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
       call mpp_write_meta(unit, fields(k), axes(1:numaxes), vars(k)%name, vars(k)%units, &
           vars(k)%longname, pack = pack, checksum=checksums(k,:))
    else
       call mpp_write_meta(unit, fields(k), axes(1:numaxes), vars(k)%name, vars(k)%units, &
           vars(k)%longname, pack = pack)
    endif
  enddo

  if (use_lath) call write_field(unit, axis_lath)
  if (use_latq) call write_field(unit, axis_latq)
  if (use_lonh) call write_field(unit, axis_lonh)
  if (use_lonq) call write_field(unit, axis_lonq)
  if (use_layer) call write_field(unit, axis_layer)
  if (use_int) call write_field(unit, axis_int)
  if (use_periodic) call write_field(unit, axis_periodic)

end subroutine create_file


!> This routine opens an existing NetCDF file for output.  If it
!! does not find the file, a new file is created.  It also sets up
!! structures that describe this file and the variables that will
!! later be written to this file.
subroutine reopen_file(unit, filename, vars, novars, fields, threading, timeunit, G, dG, GV)
  integer,               intent(out)   :: unit       !< unit id of an open file or -1 on a
                                                     !! nonwriting PE with single file output
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

  type(MOM_domain_type), pointer :: Domain => NULL()
  character(len=200) :: check_name, mesg
  integer :: length, ndim, nvar, natt, ntime, thread
  logical :: exists, one_file, domain_set

  thread = SINGLE_FILE
  if (PRESENT(threading)) thread = threading

  check_name = filename
  length = len(trim(check_name))
  if (check_name(length-2:length) /= ".nc") check_name = trim(check_name)//".nc"
  if (thread /= SINGLE_FILE) check_name = trim(check_name)//".0000"

  inquire(file=check_name,EXIST=exists)

  if (.not.exists) then
    call create_file(unit, filename, vars, novars, fields, threading, timeunit, &
                     G=G, dG=dG, GV=GV)
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
      call mpp_open_file(unit, filename, MPP_APPEND, MPP_NETCDF, threading=thread)
    else
      call mpp_open_file(unit, filename, MPP_APPEND, MPP_NETCDF, domain=Domain%mpp_domain)
    endif
    if (unit < 0) return

    call mpp_get_info(unit, ndim, nvar, natt, ntime)

    if (nvar == -1) then
      write (mesg,*) "Reopening file ",trim(filename)," apparently had ",nvar,&
                     " variables. Clobbering and creating file with ",novars," instead."
      call MOM_error(WARNING,"MOM_io: "//mesg)
      call create_file(unit, filename, vars, novars, fields, threading, timeunit, G=G, GV=GV)
    elseif (nvar /= novars) then
      write (mesg,*) "Reopening file ",trim(filename)," with ",novars,&
                     " variables instead of ",nvar,"."
      call MOM_error(FATAL,"MOM_io: "//mesg)
    endif

    if (nvar>0) call mpp_get_fields(unit,fields(1:nvar))

    ! Check the field names...
!    do i=1,nvar
!      call mpp_get_field_atts(fields(i),name)
!      !if (trim(name) /= trim(vars%name) then
!      !write (mesg,'("Reopening file ",a," variable ",a," is called ",a,".")',&
!      !    filename,vars%name,name)
!      !call MOM_error(NOTE,"MOM_io: "//mesg)
!    enddo
  endif

end subroutine reopen_file

!> register an axis to a domain-decomposed file
!! This routine specifies 
subroutine MOM_register_diagnostic_axis(fileObj, axis_name, axis_length)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObj !< file object returned by prior call to open_file
  character(len=*), intent(in) :: axis_name !< name of the restart file axis to register to file
  integer, optional, intent(in) :: axis_length !< length of axis/dimension ;only needed for Layer, Interface, Time,
                                                !! Period
  select case (trim(axis_name))      
    case ('latq'); call register_axis(fileObj,'latq','y', domain_position=NORTH_FACE)
    case ('lath'); call register_axis(fileObj,'lath','y', domain_position=CENTER) 
    case ('lonq'); call register_axis(fileObj,'lonq','x', domain_position=EAST_FACE) 
    case ('lonh'); call register_axis(fileObj,'lonh','x', domain_position=CENTER)
    case default
      if (.not. present(axis_length)) call MOM_error(FATAL,"MOM_io:register_axis_DD: "//&
                        "An axis_length argument is required to register the axis "//trim(axis_name))
      call register_axis(fileObj, axis_name, axis_length) 
  end select
end subroutine MOM_register_diagnostic_axis

!> register an axis to a non domain-decomposed file
subroutine MOM_register_axis_noDD(fileObj, axis_name, axis_length)
  type(FmsNetcdfFile_t), intent(inout) :: fileObj !< file object returned by prior call to open_file
  character(len=*), intent(in) :: axis_name !< name of the restart file axis to register to file
  integer, intent(in) :: axis_length !< length of axis/dimension
  call register_axis(fileObj, axis_name, axis_length)
end subroutine MOM_register_axis_noDD

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
  integer :: x_pos, y_pos
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
  
  call get_horizontal_grid_logic(hor_grid, use_lath,use_lonh,use_latq,use_lonq)
  
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
        call MOM_error(FATAL, "MOM_io: get_dimension_var_features: "//&
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

!> get the axis data from the name and return the 
!! structure with data and meta data
subroutine MOM_get_axis_data(axis_data_CS, axis_name, axis_number, & 
                             G, dG, GV, time_val, time_units)

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
  axis_data_CS%axis(axis_number)%x_position = CENTER
  axis_data_CS%axis(axis_number)%y_position = CENTER
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
        axis_data_CS%axis(axis_number)%axis  = 'Y'
        axis_data_CS%axis(axis_number)%is_domain_decomposed  = .true.
     case('lonh')
        axis_data_CS%data(axis_number)%p(isg:ieg)=>gridLonT(isg:ieg)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%horgrid_position  = CENTER
        axis_data_CS%axis(axis_number)%longname = 'Longitude'
        axis_data_CS%axis(axis_number)%units = x_axis_units
        axis_data_CS%axis(axis_number)%axis = 'X'
        axis_data_CS%axis(axis_number)%is_domain_decomposed  = .true.
     case('latq')
        axis_data_CS%data(axis_number)%p(JsgB:JegB)=>gridLatB(JsgB:JegB)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Latitude'
        axis_data_CS%axis(axis_number)%units = y_axis_units
        axis_data_CS%axis(axis_number)%axis = 'Y'
        axis_data_CS%axis(axis_number)%horgrid_position = NORTH_FACE
        axis_data_CS%axis(axis_number)%is_domain_decomposed = .true.
     case('lonq')
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%data(axis_number)%p(IsgB:IegB)=>gridLonB(IsgB:IegB)
        axis_data_CS%axis(axis_number)%longname  = 'Longitude'
        axis_data_CS%axis(axis_number)%units = x_axis_units
        axis_data_CS%axis(axis_number)%axis = 'X'
        axis_data_CS%axis(axis_number)%horgrid_position = EAST_FACE
        axis_data_CS%axis(axis_number)%is_domain_decomposed = .true.
     case('Layer')
        axis_data_CS%data(axis_number)%p=>GV%sLayer(1:GV%ke)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Layer pseudo-depth, -z*'
        axis_data_CS%axis(axis_number)%units = GV%zAxisUnits
        axis_data_CS%axis(axis_number)%axis = 'Z'
        axis_data_CS%axis(axis_number)%positive  = 'up'
        axis_data_CS%axis(axis_number)%horgrid_position  = CENTER ! dummy value for the domain-decomposed write
     case('Interface')
        axis_data_CS%data(axis_number)%p=>GV%sInterface(1:GV%ke+1)
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Interface pseudo-depth, -z*'
        axis_data_CS%axis(axis_number)%units = GV%zAxisUnits
        axis_data_CS%axis(axis_number)%axis = 'Z'
        axis_data_CS%axis(axis_number)%positive = 'up'
        axis_data_CS%axis(axis_number)%horgrid_position = CENTER ! dummy value for the domain-decomposed write
     case('Time')
        
        if (.not.(present(time_val))) then
           call MOM_error(FATAL, "MOM_io::get_axis_data: requires time_val"//&
                          " and time_units arguments for "//trim(axis_name))
        endif

        axis_data_CS%data(axis_number)%p=>time_val
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname  = 'Time'

        if (present(time_units)) then
           axis_data_CS%axis(axis_number)%units  = time_units
        else
           axis_data_CS%axis(axis_number)%units  = 'days'
        endif
        axis_data_CS%axis(axis_number)%axis = 'T'
        axis_data_CS%axis(axis_number)%horgrid_position = CENTER ! dummy value for the domain-decomposed write
     case('Period')
        if (.not.(present(time_val))) then
           call MOM_error(FATAL, "MOM_io::get_axis_data: requires a time_val argument"//&
                          " for "//trim(axis_name))
        endif

        axis_data_CS%data(axis_number)%p=>time_val
        axis_data_CS%axis(axis_number)%name = trim(axis_name)
        axis_data_CS%axis(axis_number)%longname = 'Periods for cyclical variables'
        axis_data_CS%axis(axis_number)%horgrid_position = CENTER ! dummy value for the domain-decomposed write
        axis_data_CS%axis(axis_number)%axis = 'T'
     case default
        call MOM_error(WARNING, "MOM_io::get_axis_data:"//trim(axis_name)//&
                       " is an unrecognized axis")
  end select

end subroutine MOM_get_axis_data

!> get the position parameter value from the horizontal grid (hor_grid) string id
subroutine get_horizontal_grid_position(grid_string_id, x_pos,y_pos)
  character(len=*), intent(in) :: grid_string_id !< horizontal grid string
  integer, intent(out) :: x_pos, y_pos !< integers corresponding to the x and y grid positions

  select case (trim(grid_string_id))
     case ('h') ; x_pos = CENTER; y_pos = CENTER
     case ('q') ; x_pos = EAST_FACE; y_pos=NORTH_FACE
     case ('u') ; x_pos = EAST_FACE; y_pos=CENTER
     case ('v') ; x_pos = CENTER; y_pos=NORTH_FACE
     case ('T')  ; x_pos = CENTER; y_pos=CENTER
     case ('Bu') ; x_pos = EAST_FACE; y_pos=NORTH_FACE
     case ('Cu') ; x_pos = EAST_FACE; y_pos=CENTER
     case ('Cv') ; x_pos = CENTER; y_pos=NORTH_FACE
     case ('1') ; x_pos = 0; y_pos=0 
     case default
        call MOM_error(FATAL, "MOM_io:get_horizontal_grid_position "//&
                        "Unrecognized grid_string_id argument "//trim(grid_string_id))
  end select

end subroutine get_horizontal_grid_position

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
     case ('h') ; use_lath = .true. ; use_lonh = .true.
     case ('q') ; use_latq = .true. ; use_lonq = .true.
     case ('u') ; use_lath = .true. ; use_lonq = .true.
     case ('v') ; use_latq = .true. ; use_lonh = .true.
     case ('T')  ; use_lath = .true. ; use_lonh = .true.
     case ('Bu') ; use_latq = .true. ; use_lonq = .true.
     case ('Cu') ; use_lath = .true. ; use_lonq = .true.
     case ('Cv') ; use_latq = .true. ; use_lonh = .true.
     case ('1') ; 
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
   character(len=10) :: time_units !< time units
   character(len=10) :: time_units_out !< time units trimmed
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

!> Read the data associated with a named axis in a file
subroutine read_axis_data(filename, axis_name, var)
  character(len=*),   intent(in)  :: filename  !< Name of the file to read
  character(len=*),   intent(in)  :: axis_name !< Name of the axis to read
  real, dimension(:), intent(out) :: var       !< The axis location data

  integer :: i,len,unit, ndim, nvar, natt, ntime
  logical :: axis_found
  type(axistype), allocatable :: axes(:)
  type(axistype) :: time_axis
  character(len=32) :: name, units

  call mpp_open_file(unit, trim(filename), action=MPP_RDONLY, form=MPP_NETCDF, &
                 threading=MPP_MULTI, fileset=SINGLE_FILE)

!Find the number of variables (nvar) in this file
  call mpp_get_info(unit, ndim, nvar, natt, ntime)
! -------------------------------------------------------------------
! Allocate space for the number of axes in the data file.
! -------------------------------------------------------------------
  allocate(axes(ndim))
  call mpp_get_axes(unit, axes, time_axis)

  axis_found = .false.
  do i = 1, ndim
    call mpp_get_atts(axes(i), name=name,len=len,units=units)
    if (name == axis_name) then
      axis_found = .true.
      call get_axis_data(axes(i),var)
      exit
    endif
  enddo

  if (.not.axis_found) call MOM_error(FATAL, "MOM_io read_axis_data: "//&
    "Unable to find axis "//trim(axis_name)//" in file "//trim(filename))

  deallocate(axes)

end subroutine read_axis_data

!> This function determines how many time levels a variable has.
function num_timelevels(filename, varname, min_dims) result(n_time)
  character(len=*),  intent(in) :: filename   !< name of the file to read
  character(len=*),  intent(in) :: varname    !< variable whose number of time levels
                                              !! are to be returned
  integer, optional, intent(in) :: min_dims   !< The minimum number of dimensions a variable must have
                                              !! if it has a time dimension.  If the variable has 1 less
                                              !! dimension than this, then 0 is returned.
  integer :: n_time                           !< number of time levels varname has in filename

  logical :: found
  character(len=200) :: msg
  character(len=nf90_max_name) :: name
  integer :: ncid, nvars, status, varid, ndims, n
  integer, allocatable :: varids(:)
  integer, dimension(nf90_max_var_dims) :: dimids

  n_time = -1
  found = .false.

  status = NF90_OPEN(filename, NF90_NOWRITE, ncid)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,"num_timelevels: "//&
        " Difficulties opening "//trim(filename)//" - "//&
        trim(NF90_STRERROR(status)))
    return
  endif

  status = NF90_INQUIRE(ncid, nVariables=nvars)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,"num_timelevels: "//&
        " Difficulties getting the number of variables in file "//&
        trim(filename)//" - "//trim(NF90_STRERROR(status)))
    return
  endif

  if (nvars < 1) then
    call MOM_error(WARNING,"num_timelevels: "//&
        " There appear not to be any variables in "//trim(filename))
    return
  endif


  allocate(varids(nvars))

  status = nf90_inq_varids(ncid, nvars, varids)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,"num_timelevels: "//&
        " Difficulties getting the variable IDs in file "//&
        trim(filename)//" - "//trim(NF90_STRERROR(status)))
    deallocate(varids) ; return
  endif

  do n = 1,nvars
    status = nf90_inquire_variable(ncid, varids(n), name=name)
    if (status /= NF90_NOERR) then
      call MOM_error(WARNING,"num_timelevels: "//&
          " Difficulties getting a variable name in file "//&
          trim(filename)//" - "//trim(NF90_STRERROR(status)))
    endif

    if (trim(lowercase(name)) == trim(lowercase(varname))) then
      if (found) then
        call MOM_error(WARNING,"num_timelevels: "//&
          " Two variables match the case-insensitive name "//trim(varname)//&
          " in file "//trim(filename)//" - "//trim(NF90_STRERROR(status)))
      else
        varid = varids(n) ; found = .true.
      endif
    endif
  enddo

  deallocate(varids)

  if (.not.found) then
    call MOM_error(WARNING,"num_timelevels: "//&
        " variable "//trim(varname)//" was not found in file "//&
        trim(filename))
    return
  endif

  status = nf90_inquire_variable(ncid, varid, ndims = ndims)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,"num_timelevels: "//&
      trim(NF90_STRERROR(status))//" Getting number of dimensions of "//&
      trim(varname)//" in "//trim(filename))
    return
  endif

  if (present(min_dims)) then
    if (ndims < min_dims-1) then
      write(msg, '(I3)') min_dims
      call MOM_error(WARNING, "num_timelevels: variable "//trim(varname)//&
        " in file "//trim(filename)//" has fewer than min_dims = "//trim(msg)//&
        " dimensions.")
    elseif (ndims == min_dims - 1) then
      n_time = 0 ; return
    endif
  endif

  status = nf90_inquire_variable(ncid, varid, dimids = dimids(1:ndims))
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,"num_timelevels: "//&
      trim(NF90_STRERROR(status))//" Getting last dimension ID for "//&
      trim(varname)//" in "//trim(filename))
    return
  endif

  status = nf90_Inquire_Dimension(ncid, dimids(ndims), len=n_time)
  if (status /= NF90_NOERR) call MOM_error(WARNING,"num_timelevels: "//&
      trim(NF90_STRERROR(status))//" Getting number of time levels of "//&
      trim(varname)//" in "//trim(filename))

  return

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

!> Open domain-decomposed file(s) with the base file name
!! 'filename' to read or write/append data.  
!! The domain comes from the ocean_grid_type structure G.
function MOM_open_file_DD_ocean_grid(MOMfileObj, filename, mode, G, is_restart) result(file_open_success)
  type(FmsNetcdfDomainFile_t), intent(inout) :: MOMfileObj !< netCDF file object 
  character(len=*),       intent(in) :: filename !< The base filename of the file(s) to search for
  character(len=*),       intent(in) :: mode !< read or write(checks if file exists to append)
  type(ocean_grid_type),      intent(in) :: G !< The ocean's grid structure
  logical, intent(in) :: is_restart !< indicates whether to check for restart file(s)

  logical :: file_open_success !< returns .true. if the file(s) is(are) opened
  character(len=512) :: mesg      ! A message for warnings.
   
  select case (trim(mode))
     case("read")
        file_open_success = open_file(MOMfileObj, filename, "read", & 
                          G%Domain%mpp_domain, is_restart = is_restart)
     case("write")
        ! check if file(s) already exists and can be appended
        file_open_success = open_file(MOMfileObj, filename, "append", & 
                                   G%Domain%mpp_domain, is_restart = is_restart)
        if (.not.(file_open_success)) then
           ! create and open new file(s) for domain-decomposed write
           file_open_success = open_file(MOMfileObj, filename, "write", & 
                                   G%Domain%mpp_domain, is_restart = is_restart)
        endif
     case("overwrite")
        ! check if file(s) already exists and can be overwritten
        file_open_success = open_file(MOMfileObj, filename, "overwrite", &
                                  G%Domain%mpp_domain, is_restart = is_restart) 
     case default
        write(mesg,'( "ERROR, file mode must be read or write to open ",A)') trim(filename)
        call MOM_error(FATAL,"MOM_io::MOM_open_file_DD_ocean_grid: "//mesg)
  end select     
end function MOM_open_file_DD_ocean_grid

!> Open domain-decomposed file with the base file name
!! 'filename' to read from or write/append to. The domain comes from the MOM_domain_type structure G.
function MOM_open_file_DD_supergrid(MOMfileObj, filename, mode, G, is_restart) result(file_open_success)
  type(FmsNetcdfDomainFile_t), intent(inout) :: MOMfileObj !< netCDF file object 
  character(len=*),       intent(in) :: filename !< The base filename of the file(s) to search for
  character(len=*),       intent(in) :: mode !< read or write(checks if file exists to append)
  type(MOM_domain_type),  intent(in)  :: G ! Supergrid domain defined in MOM_grid_initialize.F90
  logical, intent(in) :: is_restart !< indicates whether to check for restart file(s)

  logical :: file_open_success !< returns .true. if the file(s) is(are) opened
  character(len=512) :: mesg      ! A message for warnings.
   
  select case (trim(mode))
     case("read")
        file_open_success = open_file(MOMfileObj, filename, "read", & 
                          G%mpp_domain, is_restart = is_restart)
     case("write")
        ! check if file(s) already exists and can be appended
        file_open_success = open_file(MOMfileObj, filename, "append", & 
                                   G%mpp_domain, is_restart = is_restart)
        if (.not.(file_open_success)) then
           ! create and open new file(s) for domain-decomposed write
           file_open_success = open_file(MOMfileObj, filename, "write", & 
                                   G%mpp_domain, is_restart = is_restart)
        endif
     case("overwrite")
        ! check if file(s) already exists and can be overwritten
        file_open_success = open_file(MOMfileObj, filename, "overwrite", & 
                                   G%mpp_domain, is_restart = is_restart)
     case default
        write(mesg,'( "ERROR, file mode must be read or write to open ",A)') trim(filename)
        call MOM_error(FATAL,"MOM_io::MOM_open_file_DD_supergrid: "//mesg)
  end select     
end function MOM_open_file_DD_supergrid

!> Open domain-decomposed file with the base file name
!! 'filename' to read, or write/append data. 
!! The domain comes from the dyn_horgrid_type structure G.
function MOM_open_file_DD_dyn_horgrid(MOMfileObj, filename, mode, G, is_restart) result(file_open_success)
  type(FmsNetcdfDomainFile_t), intent(inout) :: MOMfileObj !< netCDF file object 
  character(len=*),       intent(in) :: filename !< The base filename of the file(s) to search for
  character(len=*),       intent(in) :: mode !< read or write(checks if file exists to append)
  type(dyn_horgrid_type),  intent(in)  :: G ! Supergrid domain defined in MOM_grid_initialize.F90
  logical, intent(in) :: is_restart !< indicates whether to check for restart file(s)

  logical :: file_open_success !< returns .true. if the file(s) is(are) opened
  character(len=512) :: mesg      ! A message for warnings.
   
  select case (trim(mode))
     case("read")
        file_open_success = open_file(MOMfileObj, filename, "read", & 
                          G%Domain%mpp_domain, is_restart = is_restart)
     case("write")
        ! check if file(s) already exists and can be appended
        file_open_success = open_file(MOMfileObj, filename, "append", & 
                                   G%Domain%mpp_domain, is_restart = is_restart)
        if (.not.(file_open_success)) then
           ! create and open new file(s) for domain-decomposed write
           file_open_success = open_file(MOMfileObj, filename, "write", & 
                                   G%Domain%mpp_domain, is_restart = is_restart)
        endif
     case("overwrite")
        ! check if file(s) already exists and can be overwritten
        file_open_success = open_file(MOMfileObj, filename, "overwrite", & 
                                   G%Domain%mpp_domain, is_restart = is_restart)
     case default
        write(mesg,'( "ERROR, file mode must be read or write to open ",A)') trim(filename)
        call MOM_error(FATAL,"MOM_io::MOM_open_file_DD_supergrid: "//mesg)
  end select     
end function MOM_open_file_DD_dyn_horgrid
!> Open non-domain-decomposed file(s) with the base file name
!! 'filename' to read, write/append data.
function MOM_open_file_noDD(MOMfileObj, filename, mode, is_restart) result(file_open_success)
  type(FmsNetcdfFile_t), intent(inout) :: MOMfileObj !< netCDF file object 
  character(len=*),       intent(in) :: filename !< The base filename of the file(s) to search for
  character(len=*),       intent(in) :: mode !< read or write(checks if file exists to append)
  logical, intent(in) :: is_restart !< indicates whether to check for restart file(s)

  logical :: file_open_success !< returns .true. if the file(s) is(are) opened
  character(len=512) :: mesg      ! A message for warnings.
   
  file_open_success = .false.
  select case (trim(mode))
     case("read")
        file_open_success = open_file(MOMfileObj, filename, "read", & 
                          is_restart = is_restart)
     case("write")
        ! check if file(s) already exists and can be appended
        file_open_success = open_file(MOMfileObj, filename, "append", & 
                                   is_restart = is_restart)
        if (.not.(file_open_success)) then
           ! create and open new file(s) for non-domain-decomposed write
           file_open_success = open_file(MOMfileObj, filename, "write", & 
                                   is_restart = is_restart)
        endif
     case("overwrite")
        ! check if file(s) already exists and can be overwritten
        file_open_success = open_file(MOMfileObj, filename, "overwrite", & 
                                   is_restart = is_restart)
     case default
        write(mesg,'( "ERROR, file mode must be read or write to open ",A)') trim(filename)
        call MOM_error(FATAL,"MOM_io::MOM_open_file_DD_ocean_grid: "//mesg)
  end select     
end function MOM_open_file_noDD

!> This function uses the fms_io function read_data to read 1-D
!! data field named "fieldname" from file "filename".
subroutine MOM_read_data_1d(filename, fieldname, data, timelevel, scale)
!  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObj !< netCDF file object 
  character(len=*),       intent(in)    :: filename  !< The name of the file to read
  character(len=*),       intent(in)    :: fieldname !< The variable name of the data in the file
  real, dimension(:),     intent(inout) :: data      !< The 1-dimensional array into which the data
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  real,         optional, intent(in)    :: scale     !< A scaling factor that the field is multiplied
                                                     !! by before they are returned.

  call old_fms_read_data(filename, fieldname, data, timelevel=timelevel, no_domain=.true.)
!  call read_data(fileObj, fieldname, data)

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
                                                     !! by before they are returned.

  integer :: is, ie, js, je

  call old_fms_read_data(filename, fieldname, data, MOM_Domain%mpp_domain, &
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
                                                     !! by before they are returned.

  integer :: is, ie, js, je

  call old_fms_read_data(filename, fieldname, data, MOM_Domain%mpp_domain, &
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
                                                     !! by before they are returned.

  integer :: is, ie, js, je

  call old_fms_read_data(filename, fieldname, data, MOM_Domain%mpp_domain, &
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
subroutine MOM_read_vector_2d(fileObj, u_fieldname, v_fieldname, u_data, v_data, MOM_Domain, &
                              timelevel, stagger, scalar_pair, scale)
  type(FmsNetcdfDomainFile_t) :: fileObj !< netcdf file object returned by call to MOM_open_file
  character(len=*),       intent(in)    :: u_fieldname !< The variable name of the u data in the file
  character(len=*),       intent(in)    :: v_fieldname !< The variable name of the v data in the file
  real, dimension(:,:),   intent(inout) :: u_data    !< The 2 dimensional array into which the
                                                     !! u-component of the data should be read
  real, dimension(:,:),   intent(inout) :: v_data    !< The 2 dimensional array into which the
                                                     !! v-component of the data should be read
  type(MOM_domain_type),  intent(in)    :: MOM_Domain !< The MOM_Domain that describes the decomposition
  integer,      optional, intent(in)    :: timelevel !< The time level in the file to read
  integer,      optional, intent(in)    :: stagger   !< A flag indicating where this vector is discretized
  logical,      optional, intent(in)    :: scalar_pair !< If true, a pair of scalars are to be read.cretized
  real,         optional, intent(in)    :: scale     !< A scaling factor that the fields are multiplied
                                                     !! by before they are returned.
  integer :: is, ie, js, je
  integer :: u_pos, v_pos
  integer :: start(2), nread(2), dim_sizes_u(2), dim_sizes_v(2)
  character(len=32), dimension(2) :: dim_names_u, dim_names_v, units_u, units_v
  
  if (.not. check_if_open(fileObj)) call MOM_error(FATAL, "MOM_read_vector_2d: netcdf fileObj not open.")

  u_pos = EAST_FACE ; v_pos = NORTH_FACE
  if (present(stagger)) then
    if (stagger == CGRID_NE .or. stagger == BGRID_NE ) then ; u_pos = EAST_FACE ; v_pos = NORTH_FACE
    elseif (stagger == AGRID) then ; u_pos = CENTER ; v_pos = CENTER ; endif
  endif
  
  start(:) = 1
  call get_variable_size(fileObj, u_fieldname, dim_sizes_u, broadcast=.true.)
  call get_variable_size(fileObj, v_fieldname, dim_sizes_v, broadcast=.true.)

  if present(timelevel) then
    start(1) = timelevel 
    dim_sizes_u(1) = timelevel
    dim_sizes_v(1) = timelevel
  endif
  !call old_fms_read_data(filename, u_fieldname, u_data, MOM_Domain%mpp_domain, &
   !              timelevel=timelevel, position=u_pos)
  !call old_fms_read_data(filename, v_fieldname, v_data, MOM_Domain%mpp_domain, &
  !               timelevel=timelevel, position=v_pos)
  
  call get_variable_dimension_names(fileObj, u_fieldname, dim_names_u, broadcast=.true.)
  call get_variable_dimension_names(fileObj, v_fieldname, dim_names_v, broadcast=.true.)
 
  do i=1,2
    call get_variable_units(u_fieldname, dim_names_u(i), units_u(i))
    call get_variable_units(u_fieldname, dim_names_v(i), units_v(i))
    select case (trim(lowercase(units_u(i)))
      case ("degrees_east"); call register_axis(fileObj, dim_names_u(i), "x", domain_position=u_pos)
      case ("degrees_north"); call register_axis(fileObj, dim_names_u(i), "y", domain_position=u_pos)  
      case default
        call register_axis(fileObj, dim_names_u(i),dim_sizes_u(i))   
  
  do i=1,2
    do j=1,size(axis_names)
      call MOM_register_diagnostic_axis(fileObj, dim_names_u(i), )
      else
     
  enddo

  call read_data(fileObj,u_fieldname, u_data, corner=start, nread)
  call read_data(fileObj,v_fieldname, v_data, corner=start)

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
subroutine MOM_read_vector_3d(fileObj, filename, u_fieldname, v_fieldname, u_data, v_data, MOM_Domain, &
                              timelevel, stagger, scalar_pair, scale)
  type(FmsNetcdfDomainFile_t) :: fileObj !< netcdf file object returned by call to MOM_open_file
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

  call old_fms_read_data(filename, u_fieldname, u_data, MOM_Domain%mpp_domain, &
                 timelevel=timelevel, position=u_pos)
  call old_fms_read_data(filename, v_fieldname, v_data, MOM_Domain%mpp_domain, &
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
  type(MOM_domain_type),  intent(in) :: MOM_Domain !< The domain that describes the decomposition
  ! local
  integer :: is, ie, js, je

  if (scale_factor /= 1.0) then
    call get_simple_array_i_ind(MOM_Domain, size(data,1), is, ie)
    call get_simple_array_j_ind(MOM_Domain, size(data,2), js, je)
    data(is:ie,js:je) = scale_factor*data(is:ie,js:je)
  endif
end subroutine scale_data_2d

subroutine scale_data_3d(data, scale_factor, MOM_domain)
  real, dimension(:,:,:), intent(inout) :: data !< The 3-dimensional data array
  real, intent(in) :: scale_factor !< Scale factor
  type(MOM_domain_type),  intent(in) :: MOM_Domain !< The domain that describes the decomposition
  ! local
  integer :: is, ie, js, je

  if (scale_factor /= 1.0) then
    call get_simple_array_i_ind(MOM_Domain, size(data,1), is, ie)
    call get_simple_array_j_ind(MOM_Domain, size(data,2), js, je)
    data(is:ie,js:je,:) = scale_factor*data(is:ie,js:je,:)
  endif
end subroutine scale_data_3d

subroutine scale_data_4d(data, scale_factor, MOM_domain)
  real, dimension(:,:,:,:), intent(inout) :: data !< The 4-dimensional data array
  real, intent(in) :: scale_factor !< Scale factor
  type(MOM_domain_type),  intent(in) :: MOM_Domain !< The domain that describes the decomposition
  ! local
  integer :: is, ie, js, je

  if (scale_factor /= 1.0) then
    call get_simple_array_i_ind(MOM_Domain, size(data,1), is, ie)
    call get_simple_array_j_ind(MOM_Domain, size(data,2), js, je)
    data(is:ie,js:je,:,:) = scale_factor*data(is:ie,js:je,:,:)
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
!!   * create_file: create a new file and set up structures that are
!!       needed for subsequent output and write out the coordinates.
!!   * reopen_file: reopen an existing file for writing and set up
!!       structures that are needed for subsequent output.
!!   * open_input_file: open the indicated file for reading only.
!!   * mpp_close_file: close an open file.
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
