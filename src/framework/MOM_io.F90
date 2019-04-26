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
use fms_io_mod,           only : file_exist, field_size, read_data
use fms_io_mod,           only : field_exists => field_exist, io_infra_end=>fms_io_exit
use fms_io_mod,           only : get_filename_appendix => get_filename_appendix ! FYI: this function only trims strings if used without calling set_filename_appendix
use mpp_domains_mod,      only : domain1d, domain2d, mpp_get_domain_components
use mpp_domains_mod,      only : CENTER, CORNER, NORTH_FACE=>NORTH, EAST_FACE=>EAST
use mpp_io_mod,           only : open_file => mpp_open, close_file => mpp_close
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

use fms2_io_mod,          only: fms2_get_dimension_size => get_dimension_size, &
                                fms2_get_global_io_domain_indices => get_global_io_domain_indices, &
                                fms2_get_num_variables => get_num_variables, &
                                fms2_register_restart_field => register_restart_field, &
                                fms2_register_axis => register_axis, &
                                fms2_register_field => register_field, &
                                fms2_register_variable_attribute => register_variable_attribute, &
                                fms2_open_file => open_file, &
                                fms2_close_file => close_file, &
                                fms2_write_data => write_data, &
                                fms2_attribute_exists => variable_att_exists, &
                                fms2_variable_exists => variable_exists, &
                                fms2_dimension_exists => dimension_exists, &
                                fms2_file_exists => file_exists, &
                                FmsNetcdfDomainFile_t, unlimited

use netcdf

implicit none ; private

public :: close_file, create_file, field_exists, field_size, fieldtype, get_filename_appendix
public :: file_exists, flush_file, get_file_info, get_file_atts, get_file_fields
public :: get_file_times, open_file, read_axis_data, read_data
public :: num_timelevels, MOM_read_data, MOM_read_vector, ensembler
public :: reopen_file, slasher, write_field, write_version_number, MOM_io_init
public :: open_namelist_file, check_nml_error, io_infra_init, io_infra_end
public :: APPEND_FILE, ASCII_FILE, MULTIPLE, NETCDF_FILE, OVERWRITE_FILE
public :: READONLY_FILE, SINGLE_FILE, WRITEONLY_FILE
public :: CENTER, CORNER, NORTH_FACE, EAST_FACE
public :: var_desc, modify_vardesc, query_vardesc, cmor_long_std

public :: get_dimension_features
public :: get_variable_byte_size
public :: get_horizontal_grid_position
public :: get_time_values
public :: get_time_units
public :: MOM_get_axis_data
public :: MOM_open_file
public :: MOM_close_file
public :: MOM_register_axis
public :: MOM_register_variable_attribute
public :: MOM_write_data
public :: MOM_write_IC

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

!> Type for describing a variable axis
type, public :: axis_data_type
  character(len=64)  :: name = ''               !< Name of the axis
  character(len=48)  :: units = ''              !< Physical dimensions of the axis
  character(len=240) :: longname = ''           !< Long name of the axis
  integer   :: horgrid_position = 0             !< Horizontal grid position
  logical :: is_domain_decomposed = .false.     !< if .true. the axis data are domain-decomposed
                                                !! and need to be indexed by the compute domain
                                                !! before passing to fms2_write_data
  real, pointer,dimension(:) :: data => NULL()  !< pointer to the axis data
end type axis_data_type

!> Indicate whether a file exists, perhaps with domain decomposition
interface file_exists
  module procedure MOM_file_exists
end interface
!> Open a netCDF file 
interface MOM_open_file
  module procedure MOM_open_file_DD
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

!> Write a data field to a netCDF file
interface MOM_write_data
  module procedure MOM_write_data_4d
  module procedure MOM_write_data_3d
  module procedure MOM_write_data_2d
  module procedure MOM_write_data_1d
  module procedure MOM_write_data_0d
end interface

!> Write initial conditions to a netCDF file
interface MOM_write_IC
  module procedure MOM_write_IC_4d
  module procedure MOM_write_IC_3d
  module procedure MOM_write_IC_2d
  module procedure MOM_write_IC_1d
end interface

interface MOM_register_variable_attribute
  module procedure register_variable_attribute_string
  module procedure register_variable_attribute_integer
  module procedure register_variable_attribute_real  
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
    call open_file(unit, filename, MPP_OVERWR, MPP_NETCDF, threading=thread)
  else
    call open_file(unit, filename, MPP_OVERWR, MPP_NETCDF, domain=Domain%mpp_domain)
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
      call open_file(unit, filename, MPP_APPEND, MPP_NETCDF, threading=thread)
    else
      call open_file(unit, filename, MPP_APPEND, MPP_NETCDF, domain=Domain%mpp_domain)
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

!> register an axis to a restart file
subroutine MOM_register_axis(fileObj, axis_name, axis_length)
   type(FmsNetcdfDomainFile_t), intent(inout) :: fileObj !< file object returned by prior call to fms2_open_file
   character(len=*), intent(in) :: axis_name !< name of the restart file axis to register to file
   integer, optional, intent(in) :: axis_length !< length of axis/dimension
                                                !! (only needed for z-grid and Time)
   select case (trim(axis_name))
         case ('latq'); call fms2_register_axis(fileObj,'latq','y')
         case ('lath'); call fms2_register_axis(fileObj,'lath','y') 
         case ('lonq'); call fms2_register_axis(fileObj,'lonq','x') 
         case ('lonh'); call fms2_register_axis(fileObj,'lonh','x')
         case ('Layer')
            if (.not.(present(axis_length))) then
                call MOM_error(FATAL,"MOM_restart::register_restart_axis: "//&
                     "axis_length argument required to register the Layer axis")
            endif 
            call fms2_register_axis(fileObj,'Layer',axis_length)
         case ('Interface')
            if (.not.(present(axis_length))) then
                call MOM_error(FATAL,"MOM_restart::register_restart_axis: "//&
                     "axis_length argument required to register the Interface axis")
            endif 
            call fms2_register_axis(fileObj,'Interface',axis_length)
         case ('Time')
            if (.not.(present(axis_length))) then
                call MOM_error(FATAL,"MOM_restart::register_restart_axis: "//&
                     "axis_length argument required to register the Time axis")
            endif 
            call fms2_register_axis(fileObj,'Time', axis_length)
         case ('Period')
            if (.not.(present(axis_length))) then
                call MOM_error(FATAL,"MOM_restart::register_restart_axis: "//&
                     "axis_length argument required to register the Period axis")
            endif 
            call fms2_register_axis(fileObj,'Period',axis_length)
   end select
end subroutine MOM_register_axis

!> register a string variable attribute to a netCDF file
subroutine register_variable_attribute_string(fileObjWrite, var_name, att_name, att_value)
   type(FmsNetcdfDomainFile_t), intent(inout) :: fileObjWrite !< file object returned by prior call to fms2_open_file
   character(len=*), intent(in) :: var_name  !< Name of the variable
   character(len=*), intent(in) :: att_name  !< Name of the variable attribute to register to the file
   character(len=*), intent(in) :: att_value !< The variable attribute value

   call fms2_register_variable_attribute(fileObjWrite, trim(var_name), trim(att_name), trim(att_value))
end subroutine register_variable_attribute_string

!> register an integer variable attribute to a netCDF file
subroutine register_variable_attribute_integer(fileObjWrite, var_name, att_name, att_value)
   type(FmsNetcdfDomainFile_t), intent(inout) :: fileObjWrite !< file object returned by prior call to fms2_open_file
   character(len=*), intent(in) :: var_name  !< Name of the variable
   character(len=*), intent(in) :: att_name  !< Name of the variable attribute to register to the file
   integer, intent(in) :: att_value !< The variable attribute value

   call fms2_register_variable_attribute(fileObjWrite, trim(var_name), trim(att_name), att_value)
end subroutine register_variable_attribute_integer

!> register an integer variable attribute to a netCDF file
subroutine register_variable_attribute_real(fileObjWrite, var_name, att_name, att_value)
   type(FmsNetcdfDomainFile_t), intent(inout) :: fileObjWrite !< file object returned by prior call to fms2_open_file
   character(len=*), intent(in) :: var_name  !< Name of the variable
   character(len=*), intent(in) :: att_name  !< Name of the variable attribute to register to the file
   real, intent(in) :: att_value !< The variable attribute value

   call fms2_register_variable_attribute(fileObjWrite, trim(var_name), trim(att_name), att_value)
end subroutine register_variable_attribute_real

!> write 4d data to a netcdf file
subroutine MOM_write_data_4d(fileObjWrite, field_name, field_data)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObjWrite  !< file object returned by prior call to fms2_open_file
  character(len=*),           intent(in) :: field_name        !< Name of the field
  real, dimension(:,:,:,:),   intent(in) :: field_data        !< data to write to the file
 
  call fms2_write_data(fileObjWrite, field_name, field_data)

end subroutine MOM_write_data_4d

!> write 3d data to a netcdf file
subroutine MOM_write_data_3d(fileObjWrite, field_name, field_data)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObjWrite  !< file object returned by prior call to fms2_open_file
  character(len=*),           intent(in) :: field_name        !< Name of the field
  real, dimension(:,:,:),     intent(in) :: field_data        !< data to write to the file

  call fms2_write_data(fileObjWrite, field_name, field_data)

end subroutine MOM_write_data_3d

!> write 2d data to a netcdf file
subroutine MOM_write_data_2d(fileObjWrite, field_name, field_data)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObjWrite  !< file object returned by prior call to fms2_open_file
  character(len=*),           intent(in) :: field_name        !< Name of the field
  real, dimension(:,:),       intent(in) :: field_data        !< data to write to the file
 
  call fms2_write_data(fileObjWrite, field_name, field_data)

end subroutine MOM_write_data_2d

!> write 1d data to a netcdf file
subroutine MOM_write_data_1d(fileObjWrite, field_name, field_data)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObjWrite  !< file object returned by prior call to fms2_open_file
  character(len=*),           intent(in) :: field_name        !< Name of the field
  real, dimension(:),         intent(in) :: field_data        !< data to write to the file
 
  call fms2_write_data(fileObjWrite, field_name, field_data)
end subroutine MOM_write_data_1d

!> write 0d data to a netcdf file
subroutine MOM_write_data_0d(fileObjWrite, field_name, field_data)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObjWrite  !< file object returned by prior call to fms2_open_file
  character(len=*),           intent(in) :: field_name        !< Name of the field
  real,                        intent(in) :: field_data        !< data to write to the file
  call fms2_write_data(fileObjWrite, field_name, field_data)

end subroutine MOM_write_data_0d

!> Get the horizontal grid, vertical grid, and/or time dimension names and lengths
!! from the grid ids returned by a prior call to query_vardesc
subroutine get_dimension_features(hor_grid, z_grid, t_grid_in, G, GV, &
                                  dim_names, dim_length, num_axes)
  character(len=*), intent(in) :: hor_grid !< horizontal grid
  character(len=*), intent(in) :: z_grid !< vertical grid
  character(len=*), intent(in) :: t_grid_in !< time grid
  type(ocean_grid_type), intent(in)  :: G !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV !< ocean vertical grid structure
  character(len=*), dimension(:), intent(out) :: dim_names !< array of dimension names
  integer, dimension(:), intent(out) :: dim_length !< array of dimension sizes
  integer, intent(out) ::  num_axes !< number of axes to register in the restart file
  
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
  
  num_axes = 0
 
  ! set the ocean grid coordinates
  gridLatT => G%gridLatT
  gridLatB => G%gridLatB
  gridLonT => G%gridLonT
  gridLonB => G%gridLonB
  isg = G%isg
  ieg = G%ieg 
  jsg = G%jsg
  jeg = G%jeg
  IsgB = G%IsgB
  IegB = G%IegB
  JsgB = G%JsgB
  JegB = G%JegB

  use_lath = .false.
  use_lonh = .false.
  use_latq = .false.
  use_lonq = .false.
  
  select case (trim(hor_grid))
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
        call MOM_error(FATAL, "MOM_io:get_dimension_features "//&
                        "Unrecognized hor_grid argument "//trim(hor_grid))
  end select
  
  ! add longitude name to dimension name array
  if (use_lonh) then
     num_axes = num_axes+1 
     dim_names(num_axes) = ''
     dim_names(num_axes)(1:len_trim('lonh')) = 'lonh'
     dim_length(num_axes) = size(gridLonT(isg:ieg))
  elseif (use_lonq) then
     num_axes = num_axes+1
     dim_names(num_axes) = ''
     dim_names(num_axes)(1:len_trim('lonq')) ='lonq'
     dim_length(num_axes) = size(gridLonB(IsgB:IegB)) 
  endif
 
  ! add latitude name to dimension name array
  if (use_lath) then
     num_axes = num_axes+1 
     dim_names(num_axes) = ''
     dim_names(num_axes)(1:len_trim('lath')) = 'lath'
     dim_length(num_axes) = size(gridLatT(jsg:jeg))
  elseif (use_latq) then
     num_axes = num_axes+1 
     dim_names(num_axes) = ''
     dim_names(num_axes)(1:len_trim('latq')) = 'latq'
     dim_length(num_axes) = size(gridLatB(JsgB:JegB))
  endif

  ! vertical grid
 
  select case (trim(z_grid))
     case ('L')
        num_axes = num_axes+1
        dim_names(num_axes) = ''
        dim_names(num_axes)(1:len_trim('Layer')) = 'Layer'
        dim_length(num_axes) = size(GV%sLayer(1:GV%ke))  
     case ('i')
        num_axes = num_axes+1
        dim_names(num_axes) = ''
        dim_names(num_axes)(1:len_trim('Interface')) = 'Interface'
        dim_length(num_axes) = size(GV%sInterface(1:GV%ke+1))
     case ('1') ! Do nothing.
     case default
        call MOM_error(FATAL, "MOM_io: get_dimension_features: "//&
                      " has an unrecognized z_grid argument"//trim(z_grid))
  end select
  
  ! time
  t_grid = adjustl(t_grid_in)
  select case (t_grid(1:1))
     case ('s', 'a', 'm')
        num_axes = num_axes+1
        dim_names(num_axes) = ''
        dim_names(num_axes)(1:len_trim('Time')) = 'Time'
        dim_length(num_axes) = unlimited
     case ('p')
        if (len_trim(t_grid(2:8)) <= 0) then
            call MOM_error(FATAL,"MOM_io:get_dimension_features: "//&
                           "No periodic axis length was specified in "//trim(t_grid))
        endif
        num_axes = num_axes+1
        dim_names(num_axes) = ''
        dim_names(num_axes)(1:len_trim('Period')) = 'Period'
        dim_length(num_axes) = unlimited
     case ('1') ! Do nothing.
     case default
           call MOM_error(WARNING, "MOM_io: get_dimension_features: "//&
                       "Unrecognized t_grid "//trim(t_grid))
  end select

end subroutine get_dimension_features

!> get the position parameter value from the horizontal grid (hor_grid) string id
function get_horizontal_grid_position(grid_string_id) result(grid_position)
  character(len=*), intent(in) :: grid_string_id !< horizontal grid string
  integer :: grid_position !< integer corresponding to the grid position

  select case (grid_string_id)
     case ('h') ; grid_position = CENTER
     case ('q') ; grid_position = CORNER
     case ('u') ; grid_position = EAST_FACE
     case ('v') ; grid_position = NORTH_FACE
     case ('T')  ; grid_position = CENTER
     case ('Bu') ; grid_position = CORNER
     case ('Cu') ; grid_position = EAST_FACE
     case ('Cv') ; grid_position = NORTH_FACE
     case ('1') ; grid_position = 0 
     case default
        call MOM_error(FATAL, "MOM_io:get_horizontal_grid_position "//&
                        "Unrecognized grid_string_id argument "//trim(grid_string_id))
  end select

end function get_horizontal_grid_position

!> get the axis data from the name and return the 
!! structure with data and meta data
subroutine MOM_get_axis_data(axis_data_CS, axis_name, G, GV, time_val, time_units)
  type(axis_data_type), intent(inout) :: axis_data_CS !< structure containing the axis data and metadata
  character(len=*), intent(in) :: axis_name !< name of the axis
  type(ocean_grid_type), intent(in) :: G !< ocean horizontal grid structure; G or dG
  type(verticalGrid_type), target, intent(in) :: GV !< ocean vertical grid structure
  real,dimension(:), target, optional, intent(in) :: time_val !< time value
  character(len=*), optional,intent(in) :: time_units!< units for non-periodic time axis
  ! local
  integer :: isg, ieg, jsg, jeg, IsgB, IegB, JsgB, JegB
  real, pointer, dimension(:) :: gridLatT => NULL(), & ! The latitude or longitude of T or B points for
     gridLatB => NULL(), & ! the purpose of labeling the output axes.
     gridLonT => NULL(), &
     gridLonB => NULL()

  ! set the ocean grid coordinates
  gridLatT => G%gridLatT
  gridLatB => G%gridLatB
  gridLonT => G%gridLonT
  gridLonB => G%gridLonB
  isg = G%isg
  ieg = G%ieg 
  jsg = G%jsg
  jeg = G%jeg
  IsgB = G%IsgB
  IegB = G%IegB
  JsgB = G%JsgB
  JegB = G%JegB

  axis_data_CS%name = ''
  axis_data_CS%name = trim(axis_name)
  axis_data_CS%data => NULL()
  axis_data_CS%longname = ''
  axis_data_CS%units = ''
  axis_data_CS%horgrid_position = 0
  axis_data_CS%is_domain_decomposed = .false.
  
  select case(trim(axis_name))
     case('lath')
        axis_data_CS%data=>gridLatT(jsg:jeg)
        axis_data_CS%longname = 'Latitude'
        axis_data_CS%units = G%y_axis_units
        axis_data_CS%horgrid_position = CENTER
        axis_data_CS%is_domain_decomposed = .true.
     case('lonh')
        axis_data_CS%data=>gridLonT(isg:ieg)
        axis_data_CS%horgrid_position = CENTER
        axis_data_CS%longname = 'Longitude'
        axis_data_CS%units = G%x_axis_units
        axis_data_CS%is_domain_decomposed = .true.
     case('latq')
        axis_data_CS%data=>gridLatB(JsgB:JegB)
        axis_data_CS%longname = 'Latitude'
        axis_data_CS%units = G%y_axis_units
        axis_data_CS%horgrid_position = CORNER
        axis_data_CS%is_domain_decomposed = .true.
     case('lonq')
        axis_data_CS%data=>gridLonB(IsgB:IegB)
        axis_data_CS%longname = 'Longitude'
        axis_data_CS%units = G%x_axis_units
        axis_data_CS%horgrid_position = CORNER
        axis_data_CS%is_domain_decomposed = .true.
     case('Layer')
        axis_data_CS%data=>GV%sLayer(1:GV%ke)
        axis_data_CS%longname = 'Layer'
        axis_data_CS%units = GV%zAxisUnits
        axis_data_CS%horgrid_position = CENTER ! dummy value for the domain-decomposed write
     case('Interface')
        axis_data_CS%data=>GV%sInterface(1:GV%ke+1)
        axis_data_CS%longname = 'Interface'
        axis_data_CS%units = GV%zAxisUnits
        axis_data_CS%horgrid_position = CENTER ! dummy value for the domain-decomposed write
     case('Time')
        if (.not.(present(time_val))) then
           call MOM_error(FATAL, "MOM_io::get_axis_data: requires time_val"//&
                          " and time_units arguments for "//trim(axis_name))
        endif
        axis_data_CS%data=>time_val
        axis_data_CS%longname = 'Time'
        if (present(time_units)) then
           axis_data_CS%units = time_units
        else
           axis_data_CS%units = 'days'
        endif
        axis_data_CS%horgrid_position = CENTER ! dummy value for the domain-decomposed write
     case('Period')
        if (.not.(present(time_val))) then
           call MOM_error(FATAL, "MOM_io::get_axis_data: requires a time_val argument"//&
                          " for "//trim(axis_name))
        endif
        axis_data_CS%data=>time_val
        axis_data_CS%longname = 'Periods for cyclical variables'
        axis_data_CS%horgrid_position = CENTER ! dummy value for the domain-decomposed write
         
     case default
        call MOM_error(WARNING, "MOM_io::get_axis_data:"//trim(axis_name)//&
                       " is an unrecognized axis")
  end select

end subroutine MOM_get_axis_data

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

!> get the time values from the t_grid string
!>@note: Be sure to deallocate the time_values array
!! at the end of the routine that calls this function to
!! avoid memory leaks.
function get_time_values(t_grid_in, array_size) result(time_values)
  character(len=*), intent(in) :: t_grid_in !< string for the time grid
  integer, optional, intent(in) :: array_size !< array size to allocate
  ! local
  character(len=12) :: t_grid
  character(len=12) :: t_grid_read
  integer :: var_periods
  integer :: k
  real, dimension(:), allocatable :: time_values
  
  t_grid = ''
  t_grid_read = ''
  var_periods = -9999999
  t_grid(1:len_trim(adjustl(t_grid_in))) = trim(adjustl(t_grid_in))

  select case (t_grid(1:1))
     case ('s', 'a', 'm') ! allocate an empty array that will be populated after the function call
        allocate(time_values(array_size))
     case ('p')
        if (len_trim(t_grid(2:8)) <= 0) then
             call MOM_error(FATAL, &
             "MOM_io::get_time_values: No periodic axis length was specified in "//&
          trim(t_grid_in) // " in the periodic axis argument")
        endif
        var_periods = -9999999
        t_grid_read = adjustl(t_grid(2:8))
        read(t_grid_read,*) var_periods

        if (var_periods == -9999999) then
             call MOM_error(FATAL, &
             "MOM_io:: get_period_value: Failed to read the number of periods from "//&
              trim(t_grid_in))
        endif

        if (var_periods < 1) then 
            call MOM_error(FATAL, "MOM_io::get_time_values: "//&
           "Period value must be positive.")
        endif
        ! Define a periodic axis array
        allocate(time_values(var_periods))
        do k=1,var_periods
           time_values(k) = real(k)
        enddo
     case ('1') ! Do nothing.
     case default
        call MOM_error(WARNING, "MOM_io::get_time_values:"//trim(t_grid_in)//&
                       " is an unrecognized t_grid value.")
  end select

end function get_time_values

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

  call open_file(unit, trim(filename), action=MPP_RDONLY, form=MPP_NETCDF, &
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

!> write a 4d initial condition field to a netCDF file
subroutine MOM_write_IC_4d(directory, filename,variable_name, field_data, variable_required, &
                       G, GV, time, hor_grid, z_grid, t_grid, longname, units)
  character(len=*), intent(in) :: directory !< location of the IC file
  character(len=*), intent(in) :: filename !< name of the IC file
  character(len=*),         intent(in) :: variable_name      !< name of variable to writie to the IC file
  real, dimension(:,:,:,:), intent(in) :: field_data     !< Field to write
  logical,                  intent(in) :: variable_required !< If true, the run will abort if this field is not
                                                      !! successfully read from the IC file.
  type(ocean_grid_type),    intent(in) :: G         !< The ocean's grid structure
  type(verticalGrid_type),  intent(in) :: GV !< ocean vertical grid structure
  type(time_type),          intent(in) :: time        !< model time                        
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizonal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
 
  ! local
  type(vardesc) :: vd
  type(FmsNetcdfDomainFile_t) :: fileObjWrite
  type(axis_data_type) :: axis_data_CS
  integer :: horgrid_position = 1
  integer :: substring_index = 0
  integer :: name_length = 0
  integer :: num_axes = 0
  integer :: i, is, ie
  integer, dimension(4) :: dim_lengths
  logical :: file_open_success = .false.
  logical :: axis_exists = .false.
  logical :: variable_exists =.false.
  character(len=200) :: base_file_name = ''
  character(len=200) :: dim_names(4)
  character(len=20) :: time_units = ''
  character(len=20) :: t_grid_read = ''
  real :: ic_time
  real, dimension(:), allocatable :: time_vals
  real, dimension(:), allocatable :: data_temp

  ! append '.nc' to the restart file name if it is missing
  substring_index = index('.nc', trim(filename))
  if (substring_index <= 0) then
      base_file_name = append_substring(trim(directory)//trim(filename),'.nc')
  else
      name_length = len(trim(directory)//trim(filename))
      base_file_name(1:name_length) = trim(directory)//trim(filename)
  endif

  vd = var_desc(variable_name, units=units, longname=longname, hor_grid=hor_grid, &
                z_grid=z_grid, t_grid=t_grid)

  horgrid_position = get_horizontal_grid_position(vd%hor_grid)

  num_axes=0
  call get_dimension_features(vd%hor_grid, vd%z_grid, vd%t_grid, G, GV, &
                              dim_names, dim_lengths, num_axes)

  if (num_axes <= 0) then
     call MOM_error(FATAL,"MOM_io::write_IC_data_4d: num_axes is an invalid value.")
  endif
  
  ! get the time units
  ic_time = time_type_to_real(time) / 86400.0
  time_units = get_time_units(ic_time)
  ! get an array of time values
  if (.not.(allocated(time_vals))) then
     t_grid_read = adjustl(vd%t_grid)
     if (t_grid_read(1:1) /= 'p') then
        time_vals = get_time_values(vd%t_grid, 1)
        time_vals(1) = ic_time
     else
        time_vals = get_time_values(vd%t_grid)
     endif
  endif
 
  file_open_success = MOM_open_file(fileObjWrite, base_file_name, "write", G, .false.)
  
  if (.not. (file_open_success)) then
     call MOM_error(FATAL,"MOM_io::write_IC_data_4d: Failed to open file "//trim(base_file_name))
  endif
  ! register the axes, and write the axis variables to the file if they do not exist
  do i=1,num_axes
     axis_exists = fms2_dimension_exists(fileObjWrite, dim_names(i))
     if (.not.(axis_exists)) then    
        call MOM_get_axis_data(axis_data_CS, dim_names(i), G, GV, &
                                        time_vals, time_units)
        call MOM_register_axis(fileObjWrite, axis_data_CS%name, dim_lengths(i))
     endif
  enddo
  ! write the axis data to the file
  do i=1,num_axes
     variable_exists = fms2_variable_exists(fileObjWrite, dim_names(i))
     if (.not.(variable_exists)) then
        call MOM_get_axis_data(axis_data_CS, dim_names(i), G, GV, &
                                        time_vals, time_units)
        if (associated(axis_data_CS%data)) then
           call fms2_register_field(fileObjWrite, trim(axis_data_CS%name), "double", &
                                   dimensions=(/trim(axis_data_CS%name)/))
         
           allocate(data_temp(size(axis_data_CS%data)))
           data_temp = axis_data_CS%data

           if (axis_data_CS%is_domain_decomposed) then
               call fms2_get_global_io_domain_indices(fileObjWrite, trim(axis_data_CS%name), is, ie)
               call fms2_write_data(fileObjWrite, axis_data_CS%name, data_temp(is:ie))
           else
               call fms2_write_data(fileObjWrite, axis_data_CS%name, data_temp) 
           endif

           deallocate(data_temp)

           call MOM_register_variable_attribute(fileObjWrite, axis_data_CS%name, &
                                                      'long_name',axis_data_CS%longname)
           call MOM_register_variable_attribute(fileObjWrite, axis_data_CS%name, &
                                                      'units',axis_data_CS%units)
        endif
     endif
  enddo
   
  call fms2_register_field(fileObjWrite, variable_name, "double", &
                          dimensions=dim_names(1:num_axes), domain_position=horgrid_position) 
  call MOM_write_data(fileObjWrite, variable_name, field_data)
  call MOM_register_variable_attribute(fileObjWrite, variable_name, 'units', units)
  call MOM_register_variable_attribute(fileObjWrite, variable_name, 'long_name', longname) 
        
  call MOM_close_file(fileObjWrite)

  if(allocated(time_vals)) deallocate(time_vals)

end subroutine MOM_write_IC_4d

!> Write a 3d initial condition field to a netCDF file
subroutine MOM_write_IC_3d(directory, filename,variable_name, field_data, variable_required, &
                       G, GV, time, hor_grid, z_grid, t_grid, longname, units)
  character(len=*), intent(in) :: directory !< location of the IC file
  character(len=*), intent(in) :: filename !< name of the IC file
  character(len=*),         intent(in) :: variable_name      !< name of variable to writie to the IC file
  real, dimension(:,:,:),   intent(in) :: field_data     !< Field to write
  logical,                  intent(in) :: variable_required !< If true, the run will abort if this field is not
                                                      !! successfully read from the IC file.
  type(ocean_grid_type),    intent(in) :: G         !< The ocean's grid structure
  type(verticalGrid_type),  intent(in) :: GV !< ocean vertical grid structure
  type(time_type),          intent(in) :: time        !< model time                        
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizonal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
 
  ! local
  type(vardesc) :: vd
  type(FmsNetcdfDomainFile_t) :: fileObjWrite
  type(axis_data_type) :: axis_data_CS
  integer :: horgrid_position = 1
  integer :: substring_index = 0
  integer :: name_length = 0
  integer :: num_axes = 0
  integer :: i, is, ie
  integer, dimension(4) :: dim_lengths
  logical :: file_open_success = .false.
  logical :: axis_exists = .false.
  logical :: variable_exists =.false.
  character(len=200) :: base_file_name = ''
  character(len=200) :: dim_names(4)
  character(len=20) :: time_units = ''
  character(len=20) :: t_grid_read = ''
  real :: ic_time
  real, dimension(:), allocatable :: time_vals
  real, dimension(:), allocatable :: data_temp

  ! append '.nc' to the restart file name if it is missing
  substring_index = index('.nc', trim(filename))
  if (substring_index <= 0) then
      base_file_name = append_substring(trim(directory)//trim(filename),'.nc')
  else
      name_length = len(trim(directory)//trim(filename))
      base_file_name(1:name_length) = trim(directory)//trim(filename)
  endif

  vd = var_desc(variable_name, units=units, longname=longname, hor_grid=hor_grid, &
                z_grid=z_grid, t_grid=t_grid)

  horgrid_position = get_horizontal_grid_position(vd%hor_grid)

  num_axes=0
  call get_dimension_features(vd%hor_grid, vd%z_grid, vd%t_grid, G, GV, &
                              dim_names, dim_lengths, num_axes)
  if (num_axes <= 0) then
     call MOM_error(FATAL,"MOM_io::write_IC_data_3d: num_axes is an invalid value.")
  endif

  ! get the time units
  ic_time = time_type_to_real(time) / 86400.0
  time_units = get_time_units(ic_time)
  ! get an array of time values
  if (.not.(allocated(time_vals))) then
     t_grid_read = adjustl(vd%t_grid)
     if (t_grid_read(1:1) /= 'p') then
        time_vals = get_time_values(vd%t_grid, 1)
        time_vals(1) = ic_time
     else
        time_vals = get_time_values(vd%t_grid)
     endif
  endif

  file_open_success = MOM_open_file(fileObjWrite, base_file_name, "write", G, .false.)
  
  if (.not. (file_open_success)) then
     call MOM_error(FATAL,"MOM_io::write_IC_data_3d: Failed to open file "//trim(base_file_name))
  endif
  ! register the axes, and write the axis variables to the file if they do not exist
  do i=1,num_axes
     axis_exists = fms2_dimension_exists(fileObjWrite, dim_names(i))
     if (.not.(axis_exists)) then
        call MOM_get_axis_data(axis_data_CS, dim_names(i), G, GV, &
                                        time_vals, time_units)
        call MOM_register_axis(fileObjWrite, axis_data_CS%name, dim_lengths(i))
     endif
  enddo
  ! write the axis data to the file
  do i=1,num_axes
     variable_exists = fms2_variable_exists(fileObjWrite, dim_names(i))
     if (.not.(variable_exists)) then
        call MOM_get_axis_data(axis_data_CS, dim_names(i), G, GV, &
                                        time_vals, time_units)
      
        if (associated(axis_data_CS%data)) then

           call fms2_register_field(fileObjWrite, axis_data_CS%name, "double", &
                                   dimensions=(/trim(axis_data_CS%name)/))

           allocate(data_temp(size(axis_data_CS%data)))
           data_temp = axis_data_CS%data

           if (axis_data_CS%is_domain_decomposed) then
               call fms2_get_global_io_domain_indices(fileObjWrite, trim(axis_data_CS%name), is, ie)
               call fms2_write_data(fileObjWrite, axis_data_CS%name, data_temp(is:ie))
           else
               call fms2_write_data(fileObjWrite, axis_data_CS%name, data_temp) 
           endif

           deallocate(data_temp)

           call MOM_register_variable_attribute(fileObjWrite, axis_data_CS%name, &
                                                      'long_name',axis_data_CS%longname)
           call MOM_register_variable_attribute(fileObjWrite, axis_data_CS%name, &
                                                      'units',axis_data_CS%units)
        endif
     endif
  enddo

  call fms2_register_field(fileObjWrite, variable_name, "double", &
                          dimensions=dim_names(1:num_axes), domain_position=horgrid_position) 
  call MOM_write_data(fileObjWrite, variable_name, field_data)
  call MOM_register_variable_attribute(fileObjWrite, variable_name, 'units', units)
  call MOM_register_variable_attribute(fileObjWrite, variable_name, 'long_name', longname) 
        
  call MOM_close_file(fileObjWrite)

  if(allocated(time_vals)) deallocate(time_vals)

end subroutine MOM_write_IC_3d

!> Write a 2d initial condition field to a netCDF file
subroutine MOM_write_IC_2d(directory, filename,variable_name, field_data, variable_required, &
                       G, GV, time, hor_grid, z_grid, t_grid, longname, units)
  character(len=*), intent(in) :: directory !< location of the IC file
  character(len=*), intent(in) :: filename !< name of the IC file
  character(len=*),         intent(in) :: variable_name      !< name of variable to writie to the IC file
  real, dimension(:,:), intent(in) :: field_data     !< Field to write
  logical,                  intent(in) :: variable_required !< If true, the run will abort if this field is not
                                                      !! successfully read from the IC file.
  type(ocean_grid_type),    intent(in) :: G         !< The ocean's grid structure
  type(verticalGrid_type),  intent(in) :: GV !< ocean vertical grid structure
  type(time_type),          intent(in) :: time        !< model time                        
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizonal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
 
  ! local
  type(vardesc) :: vd
  type(FmsNetcdfDomainFile_t) :: fileObjWrite
  type(axis_data_type) :: axis_data_CS
  integer :: horgrid_position = 1
  integer :: substring_index = 0
  integer :: name_length = 0
  integer :: num_axes
  integer :: i, is, ie
  integer, dimension(4) :: dim_lengths
  logical :: file_open_success = .false.
  logical :: axis_exists = .false.
  logical :: variable_exists = .false.
  character(len=200) :: base_file_name = ''
  character(len=200) :: dim_names(3)
  character(len=20) :: time_units = ''
  character(len=20) :: t_grid_read =''
  real :: ic_time
  real, dimension(:), allocatable :: time_vals
  real, dimension(:), allocatable :: data_temp

  ! append '.nc' to the restart file name if it is missing
  substring_index = index('.nc', trim(filename))
  if (substring_index <= 0) then
      base_file_name = append_substring(trim(directory)//trim(filename),'.nc')
  else
      name_length = len(trim(directory)//trim(filename))
      base_file_name(1:name_length) = trim(directory)//trim(filename)
  endif

  vd = var_desc(variable_name, units=units, longname=longname, hor_grid=hor_grid, &
                z_grid=z_grid, t_grid=t_grid)

  horgrid_position = get_horizontal_grid_position(vd%hor_grid)

  num_axes=0
  call get_dimension_features(vd%hor_grid, vd%z_grid, vd%t_grid, G, GV, &
                              dim_names, dim_lengths, num_axes)
  if (num_axes <= 0) then
     call MOM_error(FATAL,"MOM_io::write_IC_data_2d: num_axes is an invalid value.")
  endif
  ! get the time units
  ic_time = time_type_to_real(time) / 86400.0
  time_units = get_time_units(ic_time)
  ! get an array of time values
  if (.not.(allocated(time_vals))) then
     t_grid_read = adjustl(vd%t_grid)
     if (t_grid_read(1:1) /= 'p') then
        time_vals = get_time_values(vd%t_grid, 1)
        time_vals(1) = ic_time
     else
        time_vals = get_time_values(vd%t_grid)
     endif
  endif

  file_open_success = MOM_open_file(fileObjWrite, base_file_name, "write", G, .false.)
  
  if (.not. (file_open_success)) then
     call MOM_error(FATAL,"MOM_io::write_IC_data_2d: Failed to open file "//trim(base_file_name))
  endif
  ! register the axes, and write the axis variables to the file if they do not exist
  do i=1,num_axes
     axis_exists = fms2_dimension_exists(fileObjWrite, dim_names(i))
     if (.not.(axis_exists)) then
        call MOM_get_axis_data(axis_data_CS, dim_names(i), G, GV, &
                                        time_vals, time_units)
        call MOM_register_axis(fileObjWrite, axis_data_CS%name, dim_lengths(i))
     endif
  enddo
  ! write the axis data to the file
  do i=1,num_axes
     variable_exists = fms2_variable_exists(fileObjWrite, dim_names(i))
     if (.not.(variable_exists)) then
        call MOM_get_axis_data(axis_data_CS, dim_names(i), G, GV, &
                                        time_vals, time_units)
        if (associated(axis_data_CS%data)) then
           call fms2_register_field(fileObjWrite, axis_data_CS%name, "double", &
                                   dimensions=(/trim(axis_data_CS%name)/))

          if (axis_data_CS%is_domain_decomposed) then
               call fms2_get_global_io_domain_indices(fileObjWrite, trim(axis_data_CS%name), is, ie)
               call fms2_write_data(fileObjWrite, axis_data_CS%name, data_temp(is:ie))
           else
               call fms2_write_data(fileObjWrite, axis_data_CS%name, data_temp) 
           endif

           allocate(data_temp(size(axis_data_CS%data)))
           data_temp = axis_data_CS%data

           call MOM_write_data(fileObjWrite, axis_data_CS%name, data_temp(is:ie))
           deallocate(data_temp)

           call MOM_register_variable_attribute(fileObjWrite, axis_data_CS%name, &
                                                      'long_name',axis_data_CS%longname)
           call MOM_register_variable_attribute(fileObjWrite, axis_data_CS%name, &
                                                      'units',axis_data_CS%units)
        endif              
     endif
  enddo
   
  call fms2_register_field(fileObjWrite, variable_name, "double", &
                          dimensions=dim_names(1:num_axes), domain_position=horgrid_position) 
  call MOM_write_data(fileObjWrite, variable_name, field_data)
  call MOM_register_variable_attribute(fileObjWrite, variable_name, 'units', units)
  call MOM_register_variable_attribute(fileObjWrite, variable_name, 'long_name', longname) 
        
  call MOM_close_file(fileObjWrite)

  if(allocated(time_vals)) deallocate(time_vals)

end subroutine MOM_write_IC_2d

!> Write a 1d initial condition field to a netCDF file
subroutine MOM_write_IC_1d(directory, filename,variable_name, field_data, variable_required, &
                       G, GV, time, hor_grid, z_grid, t_grid, longname, units)
  character(len=*), intent(in) :: directory !< location of the IC file
  character(len=*), intent(in) :: filename !< name of the IC file
  character(len=*),         intent(in) :: variable_name      !< name of variable to writie to the IC file
  real, dimension(:), intent(in) :: field_data     !< Field to write
  logical,                  intent(in) :: variable_required !< If true, the run will abort if this field is not
                                                      !! successfully read from the IC file.
  type(ocean_grid_type),    intent(in) :: G         !< The ocean's grid structure
  type(verticalGrid_type),  intent(in) :: GV !< ocean vertical grid structure
  type(time_type),          intent(in) :: time        !< model time                        
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizonal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
 
  ! local
  type(vardesc) :: vd
  type(FmsNetcdfDomainFile_t) :: fileObjWrite
  type(axis_data_type) :: axis_data_CS
  integer :: horgrid_position = 1
  integer :: substring_index = 0
  integer :: name_length = 0
  integer :: num_axes = 0
  integer :: i, is, ie
  integer, dimension(4) :: dim_lengths
  logical :: file_open_success = .false.
  logical :: axis_exists = .false.
  logical :: variable_exists =.false.
  character(len=200) :: base_file_name = ''
  character(len=200) :: dim_names(2)
  character(len=10) :: time_units = ''
  character(len=20) :: t_grid_read = ''
  real :: ic_time
  real, dimension(:), allocatable :: time_vals
  real, dimension(:), allocatable :: data_temp

  ! append '.nc' to the restart file name if it is missing
  substring_index = index('.nc', trim(filename))
  if (substring_index <= 0) then
      base_file_name = append_substring(trim(directory)//trim(filename),'.nc')
  else
      name_length = len(trim(directory)//trim(filename))
      base_file_name(1:name_length) = trim(directory)//trim(filename)
  endif

  vd = var_desc(variable_name, units=units, longname=longname, hor_grid=hor_grid, &
                z_grid=z_grid, t_grid=t_grid)

  horgrid_position = get_horizontal_grid_position(vd%hor_grid)

  call get_dimension_features(vd%hor_grid, vd%z_grid, vd%t_grid, G, GV, &
                              dim_names, dim_lengths, num_axes)
  if (num_axes <= 0) then
     call MOM_error(FATAL,"MOM_io::write_IC_data_1d: num_axes is an invalid value.")
  endif

  ! get the time units
  ic_time = time_type_to_real(time) / 86400.0
  time_units = get_time_units(ic_time)
  ! get an array of time values
  if (.not.(allocated(time_vals))) then
     t_grid_read = adjustl(vd%t_grid)
     if (t_grid_read(1:1) /= 'p') then
        time_vals = get_time_values(vd%t_grid, 1)
        time_vals(1) = ic_time
     else
        time_vals = get_time_values(vd%t_grid)
     endif
  endif

  file_open_success = MOM_open_file(fileObjWrite, base_file_name, "write", G, .false.)
  
  if (.not. (file_open_success)) then
     call MOM_error(FATAL,"MOM_io::write_IC_data_1d: Failed to open file "//trim(base_file_name))
  endif
  ! register the axes, and write the axis variables to the file if they do not exist
  do i=1,num_axes
     axis_exists = fms2_dimension_exists(fileObjWrite, dim_names(i))
     if (.not.(axis_exists)) then
        call MOM_get_axis_data(axis_data_CS, dim_names(i), G, GV, &
                                        time_vals, time_units)
        call MOM_register_axis(fileObjWrite, axis_data_CS%name, dim_lengths(i))
     endif
  enddo
  ! write the axis data to the file
  do i=1,num_axes
     variable_exists = fms2_variable_exists(fileObjWrite, dim_names(i))
     if (.not.(variable_exists)) then
        call MOM_get_axis_data(axis_data_CS, dim_names(i), G, GV, &
                                        time_vals, time_units)
        if (associated(axis_data_CS%data)) then
           call fms2_register_field(fileObjWrite, axis_data_CS%name, "double", &
                                   dimensions=(/trim(axis_data_CS%name)/))

           allocate(data_temp(size(axis_data_CS%data)))
           data_temp = axis_data_CS%data

           if (axis_data_CS%is_domain_decomposed) then
               call fms2_get_global_io_domain_indices(fileObjWrite, trim(axis_data_CS%name), is, ie)
               call fms2_write_data(fileObjWrite, axis_data_CS%name, data_temp(is:ie))
           else
               call fms2_write_data(fileObjWrite, axis_data_CS%name, data_temp) 
           endif
           deallocate(data_temp)

           call MOM_register_variable_attribute(fileObjWrite, axis_data_CS%name, &
                                                      'long_name',axis_data_CS%longname)
           call MOM_register_variable_attribute(fileObjWrite, axis_data_CS%name, &
                                                      'units',axis_data_CS%units)
        endif
     endif
  enddo
   
  call fms2_register_field(fileObjWrite, variable_name, "double", &
                          dimensions=dim_names(1:num_axes), domain_position=horgrid_position) 
  call MOM_write_data(fileObjWrite, variable_name, field_data)
  call MOM_register_variable_attribute(fileObjWrite, variable_name, 'units', units)
  call MOM_register_variable_attribute(fileObjWrite, variable_name, 'long_name', longname) 
        
  call MOM_close_file(fileObjWrite)

  if(allocated(time_vals)) deallocate(time_vals)

end subroutine MOM_write_IC_1d

!> Returns true if the named file exists.
function MOM_file_exists(filename)
  character(len=*),       intent(in) :: filename   !< The name of the file being inquired about

! This function uses the fms_io function file_exist to determine whether
! a named file (or its decomposed variant) exists.
  logical :: MOM_file_exists
  MOM_file_exists = fms2_file_exists(filename)
end function MOM_file_exists

!> Open domain-decomposed file(s) with the base file name
!! 'filename' to read from or write/append to
function MOM_open_file_DD(fileObj, filename, mode, G, is_restart) result(file_open_success)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObj !< netCDF file object 
  character(len=*),       intent(in) :: filename !< The base filename of the file(s) to search for
  character(len=*),       intent(in) :: mode !< read or write(checks if file exists to append)
  type(ocean_grid_type),      intent(in) :: G !< The ocean's grid structure
  logical, intent(in) :: is_restart !< indicates whether to check for restart file(s)

  logical :: file_open_success !< returns .true. if the file(s) is(are) opened
  character(len=512) :: mesg      ! A message for warnings.
   
  select case (trim(mode))
     case("read")
        file_open_success=fms2_open_file(fileObj, filename, "read", & 
                          G%Domain%mpp_domain, is_restart = is_restart)
     case("write")
        ! check if file(s) already exists and can be appended
        file_open_success=fms2_open_file(fileObj, filename, "append", & 
                                   G%Domain%mpp_domain, is_restart = is_restart)
        if (.not.(file_open_success)) then
           ! create and open new file(s) for domain-decomposed write
           file_open_success=fms2_open_file(fileObj, filename, "write", & 
                                   G%Domain%mpp_domain, is_restart = is_restart)
        endif
     case default
        write(mesg,'( "ERROR, file mode must be read or write to open ",A)') trim(filename)
        call MOM_error(FATAL,"MOM_io::MOM_open_file_DD: "//mesg)
  end select     
end function MOM_open_file_DD

!> wrapper for fms2_close_file that closes a netcdf file
subroutine MOM_close_file(fileObj)
   type(FmsNetcdfDomainFile_t), intent(inout) :: fileObj !< netCDF file object 
   call fms2_close_file(fileObj)
end subroutine MOM_close_file

!> This function uses the fms_io function read_data to read 1-D
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
                                                     !! by before they are returned.

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
                                                     !! by before they are returned.

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
                                                     !! by before they are returned.

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
