!> The MOM6 facility for reading and writing restart files, and querying what has been read.
module MOM_restart

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domains, only : pe_here, num_PEs, MOM_domain_type
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_string_functions, only : lowercase, append_substring
use MOM_grid, only : ocean_grid_type
use MOM_io, only : create_file, fieldtype,open_file, close_file
use MOM_io, only : write_field, MOM_read_data, read_data, get_filename_appendix
use MOM_io, only : get_file_info, get_file_atts, get_file_fields, get_file_times
use MOM_io, only : vardesc, var_desc, query_vardesc, modify_vardesc
use MOM_io, only : MULTIPLE, NETCDF_FILE, READONLY_FILE, SINGLE_FILE
use MOM_io, only : CENTER, CORNER, NORTH_FACE, EAST_FACE
use MOM_io, only : get_horizontal_grid_coordinates, get_horizontal_grid_position
use MOM_io, only : get_vertical_grid_coordinates, get_time_coordinates
use MOM_io, only : write_axis_data
use MOM_time_manager, only : time_type, time_type_to_real, real_to_time
use MOM_time_manager, only : days_in_month, get_date, set_date
use MOM_verticalGrid, only : verticalGrid_type
use mpp_mod,         only:  mpp_chksum, mpp_pe
!use mpp_domains_mod, only : domain1d, domain2d
!use mpp_io_mod,      only: mpp_attribute_exist, mpp_get_atts
!use mpp_io_mod,      only: axistype
use fms2_io_mod,     only: fms2_register_restart_field => register_restart_field, &
                           fms2_register_axis => register_axis, &
                           fms2_register_field => register_field, &
                           fms2_write_data => write_data, &
                           fms2_read_data => read_data, &
                           fms2_open_file => open_file, &
                           fms2_close_file => close_file, &
                           fms2_get_variable_attribute => get_variable_attribute, &
                           fms2_attribute_exists => variable_att_exists, &
                           fms2_get_variable_names => get_variable_names, &
                           fms2_register_variable_attribute => register_variable_attribute, &
                           fms2_get_dimension_size => get_dimension_size, &
                           fms2_get_num_variables => get_num_variables, &
                           fms2_variable_exists => variable_exists, &
                           fms2_dimension_exists => dimension_exists, &
                           fms2_file_exists => file_exists, &
                           FmsNetcdfDomainFile_t, unlimited
#include <fms_platform.h>
!!                           
implicit none ; private

public restart_init, restart_end, restore_state, register_restart_field
public save_restart, query_initialized, restart_init_end, vardesc
public restart_files_exist, determine_is_new_run, is_new_run

!> A type for making arrays of pointers to 4-d arrays
type p4d
  real, dimension(:,:,:,:), pointer :: p => NULL() !< A pointer to a 4d array
end type p4d

!> A type for making arrays of pointers to 3-d arrays
type p3d
  real, dimension(:,:,:), pointer :: p => NULL() !< A pointer to a 3d array
end type p3d

!> A type for making arrays of pointers to 2-d arrays
type p2d
  real, dimension(:,:), pointer :: p => NULL() !< A pointer to a 2d array
end type p2d

!> A type for making arrays of pointers to 1-d arrays
type p1d
  real, dimension(:), pointer :: p => NULL() !< A pointer to a 1d array
end type p1d

!> A type for making arrays of pointers to scalars
type p0d
  real, pointer :: p => NULL() !< A pointer to a scalar
end type p0d

!> A structure with information about a single restart field
type field_restart
  type(vardesc) :: vars         !< Description of a field that is to be read from or written
                                !! to the restart file.
  logical :: mand_var           !< If .true. the run will abort if this field is not successfully
                                !! read from the restart file.
  logical :: initialized        !< .true. if this field has been read from the restart file.
  character(len=32) :: var_name !< A name by which a variable may be queried.
end type field_restart

!> A restart registry and the control structure for restarts
type, public :: MOM_restart_CS ; private
  logical :: restart    !< restart is set to .true. if the run has been started from a full restart
                        !! file.  Otherwise some fields must be initialized approximately.
  integer :: novars = 0 !< The number of restart fields that have been registered.
  logical :: parallel_restartfiles  !< If true, each PE writes its own restart file,
                                    !! otherwise they are combined internally.
  logical :: large_file_support     !< If true, NetCDF 3.6 or later is being used
                                    !! and large-file-support is enabled.
  logical :: new_run                !< If true, the input filenames and restart file existence will
                                    !! result in a new run that is not initialized from restart files.
  logical :: new_run_set = .false.  !< If true, new_run has been determined for this restart_CS.
  logical :: checksum_required      !< If true, require the restart checksums to match and error out otherwise.
                                    !! Users may want to avoid this comparison if for example the restarts are
                                    !! made from a run with a different mask_table than the current run,
                                    !! in which case the checksums will not match and cause crash.
  character(len=240) :: restartfile !< The name or name root for MOM restart files.
  logical :: restart_file_created = .false. !< If true, one or more restart files with the restartfile name root have
                                    !! were created (nc mode = 'write') during a previous call 
                                    !! to register_restart_field

  !> An array of descriptions of the registered fields
  type(field_restart), pointer :: restart_field(:) => NULL()

  !>@{ Pointers to the fields that have been registered for restarts
  type(p0d), pointer :: var_ptr0d(:) => NULL()
  type(p1d), pointer :: var_ptr1d(:) => NULL()
  type(p2d), pointer :: var_ptr2d(:) => NULL()
  type(p3d), pointer :: var_ptr3d(:) => NULL()
  type(p4d), pointer :: var_ptr4d(:) => NULL()
  !!@}
  integer :: max_fields !< The maximum number of restart fields
  type(FmsNetcdfDomainFile_t) :: fileObjRead
  type(FmsNetcdfDomainFile_t) :: fileObjWrite
  
end type MOM_restart_CS

!> Register fields for restarts
interface register_restart_field
  module procedure register_restart_field_ptr4d, register_restart_field_4d
  module procedure register_restart_field_ptr3d, register_restart_field_3d
  module procedure register_restart_field_ptr2d, register_restart_field_2d
  module procedure register_restart_field_ptr1d, register_restart_field_1d
  module procedure register_restart_field_ptr0d, register_restart_field_0d
end interface

!> Indicate whether a field has been read from a restart file
interface query_initialized
  module procedure query_initialized_name
  module procedure query_initialized_0d, query_initialized_0d_name
  module procedure query_initialized_1d, query_initialized_1d_name
  module procedure query_initialized_2d, query_initialized_2d_name
  module procedure query_initialized_3d, query_initialized_3d_name
  module procedure query_initialized_4d, query_initialized_4d_name
end interface

contains

!> Register a 3-d field for restarts, providing the metadata in a structure
subroutine register_restart_field_ptr3d(f_ptr, var_desc, mandatory, CS)
  real, dimension(:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  type(vardesc),              intent(in) :: var_desc  !< A structure with metadata about this variable
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "register_restart_field: Module must be initialized before it is used.")

  CS%novars = CS%novars+1! remove .res. from the file name since fms read automatically appends it to
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                     ! once the total number of fields is known.

  CS%restart_field(CS%novars)%vars = var_desc
  CS%restart_field(CS%novars)%mand_var = mandatory
  CS%restart_field(CS%novars)%initialized = .false.
  call query_vardesc(CS%restart_field(CS%novars)%vars, &
                     name=CS%restart_field(CS%novars)%var_name, &
                     caller="register_restart_field_ptr3d")

  CS%var_ptr3d(CS%novars)%p => f_ptr
  CS%var_ptr4d(CS%novars)%p => NULL()
  CS%var_ptr2d(CS%novars)%p => NULL()
  CS%var_ptr1d(CS%novars)%p => NULL()
  CS%var_ptr0d(CS%novars)%p => NULL()

end subroutine register_restart_field_ptr3d

!> Register a 4-d field for restarts, providing the metadata in a structure
subroutine register_restart_field_ptr4d(f_ptr, var_desc, mandatory, CS)
  real, dimension(:,:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  type(vardesc),              intent(in) :: var_desc  !< A structure with metadata about this variable
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "register_restart_field: Module must be initialized before it is used.")

  CS%novars = CS%novars+1
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                     ! once the total number of fields is known.

  CS%restart_field(CS%novars)%vars = var_desc
  CS%restart_field(CS%novars)%mand_var = mandatory
  CS%restart_field(CS%novars)%initialized = .false.
  call query_vardesc(CS%restart_field(CS%novars)%vars, &
                     name=CS%restart_field(CS%novars)%var_name, &
                     caller="register_restart_field_ptr4d")

  CS%var_ptr4d(CS%novars)%p => f_ptr
  CS%var_ptr3d(CS%novars)%p => NULL()
  CS%var_ptr2d(CS%novars)%p => NULL()
  CS%var_ptr1d(CS%novars)%p => NULL()
  CS%var_ptr0d(CS%novars)%p => NULL()

end subroutine register_restart_field_ptr4d

!> Register a 2-d field for restarts, providing the metadata in a structure
subroutine register_restart_field_ptr2d(f_ptr, var_desc, mandatory, CS)
  real, dimension(:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  type(vardesc),              intent(in) :: var_desc  !< A structure with metadata about this variable
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "register_restart_field: Module must be initialized before it is used.")

  CS%novars = CS%novars+1
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                     ! once the total number of fields is known.

  CS%restart_field(CS%novars)%vars = var_desc
  CS%restart_field(CS%novars)%mand_var = mandatory
  CS%restart_field(CS%novars)%initialized = .false.
  call query_vardesc(CS%restart_field(CS%novars)%vars, &
                     name=CS%restart_field(CS%novars)%var_name, &
                     caller="register_restart_field_ptr2d")

  CS%var_ptr2d(CS%novars)%p => f_ptr
  CS%var_ptr4d(CS%novars)%p => NULL()
  CS%var_ptr3d(CS%novars)%p => NULL()
  CS%var_ptr1d(CS%novars)%p => NULL()
  CS%var_ptr0d(CS%novars)%p => NULL()

end subroutine register_restart_field_ptr2d

!> Register a 1-d field for restarts, providing the metadata in a structure
subroutine register_restart_field_ptr1d(f_ptr, var_desc, mandatory, CS)
  real, dimension(:), target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  type(vardesc),              intent(in) :: var_desc  !< A structure with metadata about this variable
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "register_restart_field: Module must be initialized before it is used.")

  CS%novars = CS%novars+1
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                     ! once the total number of fields is known.

  CS%restart_field(CS%novars)%vars = var_desc
  CS%restart_field(CS%novars)%mand_var = mandatory
  CS%restart_field(CS%novars)%initialized = .false.
  call query_vardesc(CS%restart_field(CS%novars)%vars, &
                     name=CS%restart_field(CS%novars)%var_name, &
                     caller="register_restart_field_ptr1d")

  CS%var_ptr1d(CS%novars)%p => f_ptr
  CS%var_ptr4d(CS%novars)%p => NULL()
  CS%var_ptr3d(CS%novars)%p => NULL()
  CS%var_ptr2d(CS%novars)%p => NULL()
  CS%var_ptr0d(CS%novars)%p => NULL()

end subroutine register_restart_field_ptr1d

!> Register a 0-d field for restarts, providing the metadata in a structure
subroutine register_restart_field_ptr0d(f_ptr, var_desc, mandatory, CS)
  real,               target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  type(vardesc),              intent(in) :: var_desc  !< A structure with metadata about this variable
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "register_restart_field: Module must be initialized before it is used.")

  CS%novars = CS%novars+1
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                     ! once the total number of fields is known.

  CS%restart_field(CS%novars)%vars = var_desc
  CS%restart_field(CS%novars)%mand_var = mandatory
  CS%restart_field(CS%novars)%initialized = .false.
  call query_vardesc(CS%restart_field(CS%novars)%vars, &
                     name=CS%restart_field(CS%novars)%var_name, &
                     caller="register_restart_field_ptr0d")

  CS%var_ptr0d(CS%novars)%p => f_ptr
  CS%var_ptr4d(CS%novars)%p => NULL()
  CS%var_ptr3d(CS%novars)%p => NULL()
  CS%var_ptr2d(CS%novars)%p => NULL()
  CS%var_ptr1d(CS%novars)%p => NULL()

end subroutine register_restart_field_ptr0d

! The following provide alternate interfaces to register restarts.

!> Register a 4-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_4d(f_ptr, name, mandatory, CS, G, GV, filename, & 
                                     longname, units, hor_grid, z_grid, t_grid)
  real, dimension(:,:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  type(ocean_grid_type),      intent(in) :: G         !< The ocean's grid structure
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                      !! required if the new file uses any vertical grid axes.
  character(len=*), optional, intent(in) :: filename  !< user-specified file name if different from cs%restartfile
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizonal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent
 
  ! local
  type(vardesc) :: vd
  type(MOM_restart_CS) :: fileObjWrite
  logical :: file_open_success = .false.
  character(len=200) :: base_file_name
  character(len=200) :: restart_file_name
  character(len=200) :: dim_names(4)
  character(len=16) :: nc_action
  character(len=200) :: mesg
  integer :: name_length = 0
  integer :: num_axes, i
  integer :: substring_index = 0
  integer :: horgrid_position = 1
 
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: "//&
      "register_restart_field_4d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  WRITE(mpp_pe()+2000,*) "register_restart_field_4d: registering restart variable ", trim(name)
  call flush(mpp_pe()+2000)

  vd = var_desc(name, units=units, longname=longname, hor_grid=hor_grid, &
                z_grid=z_grid, t_grid=t_grid)

  base_file_name = ''
  restart_file_name = ''
  nc_action = ''

  if (present(filename)) then
     restart_file_name(1:len_trim(filename))=trim(filename)
  else
     restart_file_name(1:len_trim(CS%restartfile))=trim(CS%restartfile)
  endif
  ! append '.nc' to the restart file name if it is missing
  substring_index = index('.nc', trim(restart_file_name))
  if (substring_index <= 0) then
     base_file_name = append_substring(trim(restart_file_name),'.nc')
  else
     name_length = len_trim(restart_file_name)
     base_file_name(1:name_length) = trim(restart_file_name)
  endif

  ! check if file is an existing restart or initial conditions file
  CS%restart_file_created = fms2_file_exists(base_file_name)

  if (CS%restart_file_created) then
     nc_action = "append"
  else
     nc_action = "write" 
  endif                 

  ! open the restart file for domain-decomposed write
  file_open_success=fms2_open_file(CS%fileObjWrite, base_file_name, nc_action, & 
                                   G%Domain%mpp_domain, is_restart = .true.)
  if (.not. file_open_success) then 
     write(mesg,'( "ERROR, unable to open restart file ",A) ') trim(base_file_name)
     call MOM_error(FATAL,"MOM_restart:register_restart_field_0d: "//mesg)
  endif

  ! get the axis coordinate values, 
  ! check if axes are registered in file, 
  ! and register them if they are not
                                   
  ! 4d variables are lon x lat x vertical level x time
  ! horizontal grid (hor_grid)
  num_axes = 0

  call get_horizontal_grid_coordinates(CS%fileObjWrite, dim_names, num_axes, vd%hor_grid, &
                                       G, horgrid_position)
  ! Vertical (z) grid 
  call get_vertical_grid_coordinates(CS%fileObjWrite, dim_names, num_axes, GV, vd%z_grid)

  ! time (t) grid  
  call get_time_coordinates(CS%fileObjWrite, dim_names, num_axes, vd%t_grid)
  
  do i=1,num_axes
     WRITE(mpp_pe()+2000,*) "register_restart_field_4d: dim name ", trim(dim_names(i))
     call flush(mpp_pe()+2000)
  enddo

  ! register the restart field
  call register_restart_field_ptr4d(f_ptr, vd, mandatory, CS)
  ! Need to get the dimension names to register the domain-decomposed variables
  ! The NUMBER of dimensions is defined in MOM_io::create_file
  ! 1. axis metadata are written to file, returns a an axis structure with the metadaa
  ! 2. axes for a given variable are determined
  ! 3. axes metadata are associated with variable by referencing the axis structure
  
  call fms2_register_restart_field(CS%fileObjWrite, name, f_ptr, & 
       dimensions=dim_names(1:num_axes), domain_position=horgrid_position)
  
  ! register variable attributes
  if (present(units)) call fms2_register_variable_attribute(CS%fileObjWrite,name,'units',vd%units)
  if (present(longname)) call fms2_register_variable_attribute(CS%fileObjWrite,name,'long_name',vd%longname)

  call fms2_close_file(CS%fileObjWrite)
 
end subroutine register_restart_field_4d

!> Register a 3-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_3d(f_ptr, name, mandatory, CS, G, GV, filename, & 
                                     longname, units, hor_grid, z_grid, t_grid)
  real, dimension(:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  type(ocean_grid_type),      intent(in) :: G         !< The ocean's grid structure
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                      !! required if the new file uses any vertical grid axes.
  character(len=*), optional, intent(in) :: filename  !< user-specified file name if different from cs%restartfile
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizonal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent

  ! local
  type(vardesc) :: vd
  type(MOM_restart_CS) :: fileObjWrite
  logical :: file_open_success = .false.
  character(len=200) :: base_file_name 
  character(len=200) :: restart_file_name
  character(len=200) :: dim_names(4)
  character(len=16) :: nc_action
  character(len=200) :: mesg
  integer :: horgrid_position = 1
  integer :: num_axes, i
  integer :: substring_index = 0
  integer :: name_length
          
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: "//&
      "register_restart_field_3d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  WRITE(mpp_pe()+2000,*) "register_restart_field_3d: registering restart variable ", trim(name)
  call flush(mpp_pe()+2000)

  vd = var_desc(name, units=units, longname=longname, hor_grid=hor_grid, &
                z_grid=z_grid, t_grid=t_grid)
  
  base_file_name = ''
  restart_file_name = ''
  nc_action = ''

  if (present(filename)) then
     restart_file_name(1:len_trim(filename))=trim(filename)
  else
     restart_file_name(1:len_trim(CS%restartfile))=trim(CS%restartfile)
  endif
  ! append '.nc' to the restart file name if it is missing
  substring_index = index('.nc', trim(restart_file_name))
  if (substring_index <= 0) then
     base_file_name = append_substring(trim(restart_file_name),'.nc')
  else
     name_length = len_trim(restart_file_name)
     base_file_name(1:name_length) = trim(restart_file_name)
  endif

  ! check if file is an existing restart or initial conditions file
  CS%restart_file_created = fms2_file_exists(base_file_name)

  if (CS%restart_file_created) then
     nc_action = "append"
  else
     nc_action = "write" 
  endif                          

  ! open the restart file for domain-decomposed write
  file_open_success = fms2_open_file(CS%fileObjWrite, base_file_name, nc_action, &
                                      G%Domain%mpp_domain, is_restart = .true.)

  if (.not. file_open_success) then 
     write(mesg,'( "ERROR, unable to open restart file ",A) ') trim(base_file_name)
     call MOM_error(FATAL,"MOM_restart:register_restart_field_3d: "//mesg)
  endif

  ! get the axis coordinate values, 
  ! check if axes are registered in file, 
  ! and register them if they are not

  ! horizontal grid (hor_grid)
  num_axes = 0

  call get_horizontal_grid_coordinates(CS%fileObjWrite, dim_names, num_axes, vd%hor_grid, &
                                       G, horgrid_position)
  ! Vertical (z) grid 
  call get_vertical_grid_coordinates(CS%fileObjWrite, dim_names, num_axes, GV, vd%z_grid)
 
  ! time (t) grid  
  call get_time_coordinates(CS%fileObjWrite, dim_names, num_axes, vd%t_grid)

  do i=1,num_axes
     WRITE(mpp_pe()+2000,*) "register_restart_field_3d: dim name ", trim(dim_names(i))
     call flush(mpp_pe()+2000)
  enddo
 
  ! register the restart field
  call register_restart_field_ptr3d(f_ptr, vd, mandatory, CS)

  call fms2_register_restart_field(CS%fileObjWrite, name, f_ptr, & 
       dimensions=dim_names(1:num_axes), domain_position=horgrid_position)
 
  ! register variable attributes
  if (present(units)) call fms2_register_variable_attribute(CS%fileObjWrite, name,'units',vd%units)
  if (present(longname)) call fms2_register_variable_attribute(CS%fileObjWrite ,name,'long_name',vd%longname)

  call fms2_close_file(CS%fileObjWrite)

end subroutine register_restart_field_3d

!> Register a 2-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_2d(f_ptr, name, mandatory, CS, G, GV, filename, & 
                                     longname, units, hor_grid, z_grid, t_grid)
  real, dimension(:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  type(ocean_grid_type),      intent(in) :: G         !< The ocean's grid structure
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                      !! required if the new file uses any vertical grid axes.
  character(len=*), optional, intent(in) :: filename  !< user-specified file name if different from cs%restartfile
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizonal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent
  ! local
  type(vardesc) :: vd
  character(len=8) :: Zgrid
  type(MOM_restart_CS) :: fileObjWrite
  logical :: file_open_success = .false.
  character(len=200) :: base_file_name
  character(len=200) :: restart_file_name
  character(len=200) :: dim_names(3)
  character(len=16) ::  nc_action
  character(len=200) :: mesg
  integer :: horgrid_position = 1
  integer :: num_axes, i
  integer :: substring_index = 0
  integer :: name_length=0
          
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: "//&
      "register_restart_field_2d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  WRITE(mpp_pe()+2000,*) "register_restart_field_2d: registering restart variable ", trim(Name)
  call flush(mpp_pe()+2000)

  
  if (present(z_grid)) then 
     Zgrid = z_grid
  else
     Zgrid = '1' ; 
  endif

  vd = var_desc(name, units=units, longname=longname, hor_grid=hor_grid, &
                z_grid=Zgrid, t_grid=t_grid)

  base_file_name = ''
  restart_file_name = ''
  nc_action = ''

  if (present(filename)) then
     restart_file_name(1:len_trim(filename))=trim(filename)
  else
     restart_file_name(1:len_trim(CS%restartfile))=trim(CS%restartfile)
  endif
  ! append '.nc' to the restart file name if it is missing
  substring_index = index('.nc', trim(restart_file_name))
  if (substring_index <= 0) then
     base_file_name = append_substring(trim(restart_file_name),'.nc')
  else
     name_length = len_trim(restart_file_name)
     base_file_name(1:name_length) = trim(restart_file_name)
  endif

  CS%restart_file_created = fms2_file_exists(base_file_name)

  if (CS%restart_file_created) then
     nc_action = "append"
  else
     nc_action = "write" 
  endif                 

  ! open the restart file for domain-decomposed write
  file_open_success=fms2_open_file(CS%fileObjWrite, base_file_name, nc_action, &
                                   G%Domain%mpp_domain, is_restart = .true.)

  if (.not. file_open_success) then 
     write(mesg,'( "ERROR, unable to open restart file ",A) ') trim(base_file_name)
     call MOM_error(FATAL,"MOM_restart:register_restart_field_2d: "//mesg)
  endif
                     
  ! get the axis coordinate values, 
  ! check if axes are registered in file, 
  ! and register them if they are not

  ! horizontal grid (hor_grid)
  num_axes = 0

  call get_horizontal_grid_coordinates(CS%fileObjWrite, dim_names, num_axes, vd%hor_grid, &
                                       G, horgrid_position)
  if (num_axes < 2) then
     ! Vertical (z) grid 
     call get_vertical_grid_coordinates(CS%fileObjWrite, dim_names, num_axes, GV, vd%z_grid)
  endif

  ! time (t) grid   
  call get_time_coordinates(CS%fileObjWrite, dim_names, num_axes, vd%t_grid)

  do i=1,num_axes
     WRITE(mpp_pe()+2000,*) "register_restart_field_2d: dim name ", trim(dim_names(i))
     call flush(mpp_pe()+2000)
  enddo

  ! register the restart field
  call register_restart_field_ptr2d(f_ptr, vd, mandatory, CS)

  call fms2_register_restart_field(CS%fileObjWrite, name, f_ptr, & 
       dimensions=dim_names(1:num_axes), domain_position=horgrid_position)
  
  ! register variable attributes
  if (present(units)) call fms2_register_variable_attribute(CS%fileObjWrite,name,'units',vd%units)
  if (present(longname)) call fms2_register_variable_attribute(CS%fileObjWrite,name,'long_name',vd%longname)

  call fms2_close_file(CS%fileObjWrite)

end subroutine register_restart_field_2d

!> Register a 1-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_1d(f_ptr, name, mandatory, CS, G, GV, filename, & 
                                     longname, units, hor_grid, z_grid, t_grid)
  real, dimension(:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  type(ocean_grid_type),      intent(in) :: G         !< The ocean's grid structure
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                      !! required if the new file uses any vertical grid axes.
  character(len=*), optional, intent(in) :: filename  !< user-specified file name if different from cs%restartfile
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizonal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent

  ! local 
  type(vardesc) :: vd
  character(len=8) :: hgrid
  type(MOM_restart_CS) :: fileObjWrite
  logical :: file_open_success = .false.
  character(len=200) :: dim_names(2)
  character(len=200) :: base_file_name
  character(len=200) :: restart_file_name
  character(len=16) :: nc_action
  character(len=200) :: mesg
  integer :: horgrid_position = 1
  integer :: num_axes, i
  integer :: substring_index = 0
  integer :: name_length = 0

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: " // &
      "register_restart_field_3d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  WRITE(mpp_pe()+2000,*) "register_restart_field_1d: registering restart variable ", trim(Name)
  call flush(mpp_pe()+2000)

  if (present(hor_grid)) then 
     hgrid = hor_grid
  else 
     hgrid = '1'
  endif 

  vd = var_desc(name, units=units, longname=longname, hor_grid=hgrid, &
                z_grid=z_grid, t_grid=t_grid)
  base_file_name = ''
  restart_file_name = ''
  nc_action = ''

  if (present(filename)) then
     restart_file_name(1:len_trim(filename))=trim(filename)
  else
     restart_file_name(1:len_trim(CS%restartfile))=trim(CS%restartfile)
  endif
  ! append '.nc' to the restart file name if it is missing
  substring_index = index('.nc', trim(restart_file_name))
  if (substring_index <= 0) then
     base_file_name = append_substring(trim(restart_file_name),'.nc')
  else
     name_length = len_trim(restart_file_name)
     base_file_name(1:name_length) = trim(restart_file_name)
  endif

  CS%restart_file_created = fms2_file_exists(base_file_name)

  if (CS%restart_file_created) then
     nc_action = "append"
  else
     nc_action = "write" 
  endif           

  ! open the restart file for domain-decomposed write
  file_open_success=fms2_open_file(CS%fileObjWrite, base_file_name, nc_action, &
                                   G%Domain%mpp_domain, is_restart = .true.)

  if (.not. file_open_success) then 
     write(mesg,'( "ERROR, unable to open restart file ",A) ') trim(base_file_name)
     call MOM_error(FATAL,"MOM_restart:register_restart_field_1d: "//mesg)
  endif

  ! get the axis coordinate values, 
  ! check if axes are registered in file, 
  ! and register them if they are not
                              
  ! Vertical (z) grid
  num_axes = 0
  call get_vertical_grid_coordinates(CS%fileObjWrite, dim_names, num_axes, GV, vd%z_grid)

  ! time (t) grid  
  call get_time_coordinates(CS%fileObjWrite, dim_names, num_axes, vd%t_grid)
  
  do i=1,num_axes
     WRITE(mpp_pe()+2000,*) "register_restart_field_1d: dim name ", trim(dim_names(i))
     call flush(mpp_pe()+2000)
  enddo

  ! register the restart field
  call register_restart_field_ptr1d(f_ptr, vd, mandatory, CS)

  if (is_restart_file) then
     call fms2_register_restart_field(CS%fileObjWrite, name, f_ptr, & 
       dimensions=dim_names(1:num_axes), domain_position=horgrid_position)
  else
     call fms2_register_field(CS%fileObjWrite, name, "double", & 
       dimensions=dim_names(1:num_axes), domain_position=horgrid_position)
  endif

  ! register variable attributes
  if (present(units)) call fms2_register_variable_attribute(CS%fileObjWrite, name,'units',vd%units)
  if (present(longname)) call fms2_register_variable_attribute(CS%fileObjWrite, name,'long_name',vd%longname)

  call fms2_close_file(CS%fileObjWrite)

end subroutine register_restart_field_1d

!> Register a 0-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_0d(f_ptr, name, mandatory, CS, G, GV, filename, & 
                                     longname, units, t_grid)
  real,               target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  type(ocean_grid_type),      intent(in) :: G         !< The ocean's grid structure
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                      !! required if the new file uses any vertical grid axes.
  character(len=*), optional, intent(in) :: filename  !< user-specified file name if different from cs%restartfile
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent

  !local
  type(vardesc) :: vd
  type(MOM_restart_CS) :: fileObjWrite
  logical :: file_open_success = .false.
  character(len=200) :: dim_names(1)
  character(len=16) :: nc_action
  character(len=200) :: base_file_name
  character(len=200) :: restart_file_name
  character(len=200) :: mesg
  integer :: horgrid_position = 1
  integer :: num_axes, i
  integer :: substring_index = 0
  integer :: name_length
          
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: " // &
      "register_restart_field_0d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  WRITE(mpp_pe()+2000,*) "register_restart_field_0d: registering restart variable ", trim(Name)
  call flush(mpp_pe()+2000)
  
  vd = var_desc(name, units=units, longname=longname, hor_grid='1', &
                z_grid='1', t_grid=t_grid)

  base_file_name = ''
  restart_file_name = ''
  nc_action = ''

  if (present(filename)) then
     restart_file_name(1:len_trim(filename))=trim(filename)
  else
     restart_file_name(1:len_trim(CS%restartfile))=trim(CS%restartfile)
  endif
  ! append '.nc' to the restart file name if it is missing
  substring_index = index('.nc', trim(restart_file_name))
  if (substring_index <= 0) then
     base_file_name = append_substring(trim(restart_file_name),'.nc')
  else
     name_length = len_trim(restart_file_name)
     base_file_name(1:name_length) = trim(restart_file_name)
  endif

  CS%restart_file_created = fms2_file_exists(base_file_name)

  if (CS%restart_file_created) then
     nc_action = "append"
  else
     nc_action = "write" 
  endif                            

  ! open the restart file for domain-decomposed write
  file_open_success = fms2_open_file(CS%fileObjWrite, base_file_name, nc_action, &
                                      G%Domain%mpp_domain, is_restart = .true.)

  if (.not. file_open_success) then 
     write(mesg,'( "ERROR, unable to open restart file ",A) ') trim(base_file_name)
     call MOM_error(FATAL,"MOM_restart:register_restart_field_0d: "//mesg)
  endif
                 
  ! get the axis coordinate values, 
  ! check if axes are registered in file, 
  ! and register them if they are not
 
  num_axes = 0

  ! time (t) grid  
  call get_time_coordinates(CS%fileObjWrite, dim_names, num_axes, vd%t_grid)

  do i=1,num_axes
     WRITE(mpp_pe()+2000,*) "register_restart_field_0d: dim name ", trim(dim_names(i))
     call flush(mpp_pe()+2000)
  enddo

  ! register the restart field
  call register_restart_field_ptr0d(f_ptr, vd, mandatory, CS)
       
  call fms2_register_restart_field(CS%fileObjWrite, name, f_ptr, & 
                                   dimensions=dim_names, domain_position=horgrid_position)
  ! register variable attributes
  if (present(units)) call fms2_register_variable_attribute(CS%fileObjWrite, name,'units',vd%units)
  if (present(longname)) call fms2_register_variable_attribute(CS%fileObjWrite, name,'long_name',vd%longname)

  call fms2_close_file(CS%fileObjWrite)

end subroutine register_restart_field_0d

!> query_initialized_name determines whether a named field has been successfully
!! read from a restart file yet.
function query_initialized_name(name, CS) result(query_initialized)
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine returns .true. if the field referred to by name has
! initialized from a restart file, and .false. otherwise.

  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (trim(name) == CS%restart_field(m)%var_name) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.
  if ((n==CS%novars+1) .and. (is_root_pe())) &
    call MOM_error(NOTE,"MOM_restart: Unknown restart variable "//name// &
                        " queried for initialization.")

  if ((is_root_pe()) .and. query_initialized) &
    call MOM_error(NOTE,"MOM_restart: "//name// &
                         " initialization confirmed by name.")

end function query_initialized_name

!> Indicate whether the field pointed to by f_ptr has been initialized from a restart file.
function query_initialized_0d(f_ptr, CS) result(query_initialized)
  real,         target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr has
! been initialized from a restart file.

  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr0d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_0d

!> Indicate whether the field pointed to by f_ptr has been initialized from a restart file.
function query_initialized_1d(f_ptr, CS) result(query_initialized)
  real, dimension(:), target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  type(MOM_restart_CS),       pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr has
! been initialized from a restart file.

  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr1d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_1d

!> Indicate whether the field pointed to by f_ptr has been initialized from a restart file.
function query_initialized_2d(f_ptr, CS) result(query_initialized)
  real, dimension(:,:), &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr has
! been initialized from a restart file.

  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr2d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_2d

!> Indicate whether the field pointed to by f_ptr has been initialized from a restart file.
function query_initialized_3d(f_ptr, CS) result(query_initialized)
  real, dimension(:,:,:), &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr has
! been initialized from a restart file.

  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr3d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_3d

!> Indicate whether the field pointed to by f_ptr has been initialized from a restart file.
function query_initialized_4d(f_ptr, CS) result(query_initialized)
  real, dimension(:,:,:,:),  &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr has
! been initialized from a restart file.

  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr4d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_4d

!> Indicate whether the field pointed to by f_ptr or with the specified variable
!! name has been initialized from a restart file.
function query_initialized_0d_name(f_ptr, name, CS) result(query_initialized)
  real,         target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr or with the
! specified variable name has been initialized from a restart file.

  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr0d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.
  if (n==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE,"MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "probably because of the suspect comparison of pointers by ASSOCIATED.")
    query_initialized = query_initialized_name(name, CS)
  endif

end function query_initialized_0d_name

!> Indicate whether the field pointed to by f_ptr or with the specified variable
!! name has been initialized from a restart file.
function query_initialized_1d_name(f_ptr, name, CS) result(query_initialized)
  real, dimension(:),  &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr or with the
! specified variable name has been initialized from a restart file.

  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr1d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.
  if (n==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE,"MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "probably because of the suspect comparison of pointers by ASSOCIATED.")
    query_initialized = query_initialized_name(name, CS)
  endif

end function query_initialized_1d_name

!> Indicate whether the field pointed to by f_ptr or with the specified variable
!! name has been initialized from a restart file.
function query_initialized_2d_name(f_ptr, name, CS) result(query_initialized)
  real, dimension(:,:),  &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr or with the
! specified variable name has been initialized from a restart file.

  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr2d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.
  if (n==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE,"MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "probably because of the suspect comparison of pointers by ASSOCIATED.")
    query_initialized = query_initialized_name(name, CS)
  endif

end function query_initialized_2d_name

!> Indicate whether the field pointed to by f_ptr or with the specified variable
!! name has been initialized from a restart file.
function query_initialized_3d_name(f_ptr, name, CS) result(query_initialized)
  real, dimension(:,:,:),  &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr or with the
! specified variable name has been initialized from a restart file.

  integer :: m, n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr3d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.
  if (n==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE, "MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "possibly because of the suspect comparison of pointers by ASSOCIATED.")
    query_initialized = query_initialized_name(name, CS)
  endif

end function query_initialized_3d_name

!> Indicate whether the field pointed to by f_ptr or with the specified variable
!! name has been initialized from a restart file.
function query_initialized_4d_name(f_ptr, name, CS) result(query_initialized)
  real, dimension(:,:,:,:),  &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr or with the
! specified variable name has been initialized from a restart file.

  integer :: m, n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr4d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.
  if (n==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE, "MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "possibly because of the suspect comparison of pointers by ASSOCIATED.")
    query_initialized = query_initialized_name(name, CS)
  endif

end function query_initialized_4d_name

!> save_restart saves all registered variables to restart files.
subroutine save_restart(directory, time, G, CS, time_stamped, filename, GV)
  character(len=*),        intent(in)    :: directory !< The directory where the restart files
                                                  !! are to be written
  type(time_type),         intent(in)    :: time  !< The current model time
  type(ocean_grid_type),   intent(inout) :: G     !< The ocean's grid structure
  type(MOM_restart_CS),    pointer       :: CS    !< The control structure returned by a previous
                                                  !! call to restart_init.
  logical,          optional, intent(in) :: time_stamped !< If present and true, add time-stamp
                                                  !! to the restart file names.
  character(len=*), optional, intent(in) :: filename !< A filename that overrides the name in CS%restartfile.
  type(verticalGrid_type), optional, intent(in) :: GV   !< The ocean's vertical grid structure

  ! Local variables
  type(vardesc) :: vars(CS%max_fields)  ! Descriptions of the fields that
                                        ! are to be read from the restart file.
  type(fieldtype) :: fields(CS%max_fields) !
  character(len=512) :: restartpath     ! The restart file path (dir/file).
  character(len=512) :: restartpath2     ! The restart file path (dir/file).
  character(len=256) :: restartname     ! The restart file name (no dir).
  character(len=256) :: base_file_name  ! Temporary location for restart file name (no dir)
  character(len=8)   :: suffix          ! A suffix (like _2) that is appended
                                        ! to the name of files after the first.
  integer(kind=8) :: var_sz, size_in_file ! The size in bytes of each variable
                                        ! and the variables already in a file.
  integer(kind=8) :: max_file_size = 2147483647_8 ! The maximum size in bytes
                                        ! for any one file.  With NetCDF3,
                                        ! this should be 2 Gb or less.
  integer :: start_var, next_var        ! The starting variables of the
                                        ! current and next files.
  integer :: unit                       ! The mpp unit of the open file.
  integer :: m, nz, num_files, var_periods
  integer :: seconds, days, year, month, hour, minute
  character(len=8) :: hor_grid, z_grid, t_grid ! Variable grid info.
  character(len=8) :: t_grid_read
  character(len=16) :: nc_action
  character(len=64) :: var_name         ! A variable's name.
  character(len=256) :: restartnameapp
  real :: restart_time
  character(len=32) :: filename_appendix !fms appendix to filename for ensemble runs
  integer :: length
  integer(kind=8) :: check_val(CS%max_fields,1)
  integer :: isL, ieL, jsL, jeL, pos
  integer :: substring_index = 0
  type(MOM_restart_CS) :: fileObjWrite
  logical :: file_open_success = .false.
  logical :: axis_exists = .false.
  logical :: variable_exists = .false.

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "save_restart: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  ! With parallel read & write, it is possible to disable the following...

  ! The maximum file size is 4294967292, according to the NetCDF documentation.
  if (CS%large_file_support) max_file_size = 4294967292_8

  num_files = 0
  next_var = 0
  nz = 1 ; if (present(GV)) nz = GV%ke

  restart_time = time_type_to_real(time) / 86400.0

  pos = 1
  restartname = ''
  base_file_name = ''
  nc_action = ''

  if (present(filename)) then 
     base_file_name = trim(filename) 
  else
     base_file_name=trim(CS%restartfile)
  endif

  if (PRESENT(time_stamped)) then ; if (time_stamped) then
    call get_date(time,year,month,days,hour,minute,seconds)
    ! Compute the year-day, because I don't like months. - RWH
    do m=1,month-1
      days = days + days_in_month(set_date(year,m,2,0,0,0))
    enddo
    seconds = seconds + 60*minute + 3600*hour
    if (year <= 9999) then
      write(restartname,'("_Y",I4.4,"_D",I3.3,"_S",I5.5)') year, days, seconds
    elseif (year <= 99999) then
      write(restartname,'("_Y",I5.5,"_D",I3.3,"_S",I5.5)') year, days, seconds
    else
      write(restartname,'("_Y",I10.10,"_D",I3.3,"_S",I5.5)') year, days, seconds
    endif
    restartname = trim(base_file_name)//trim(restartname)
  endif ; 
 else
    restartname = trim(base_file_name)
 endif
 
  next_var = 1
  do while (next_var <= CS%novars )
     start_var = next_var
     size_in_file = 8*(2*G%Domain%niglobal+2*G%Domain%njglobal+2*nz+1000)
     
     do m=start_var,CS%novars
        call query_vardesc(CS%restart_field(m)%vars, hor_grid=hor_grid, &
                           z_grid=z_grid, t_grid=t_grid, caller="save_restart")
        if (hor_grid == '1') then
           var_sz = 8
        else
           var_sz = 8*(G%Domain%niglobal+1)*(G%Domain%njglobal+1)
        endif
        select case (z_grid)
           case ('L') ; var_sz = var_sz * nz
           case ('i') ; var_sz = var_sz * (nz+1)
        end select
        t_grid = adjustl(t_grid)
        if (t_grid(1:1) == 'p') then
           if (len_trim(t_grid(2:8)) > 0) then
              var_periods = -1
              t_grid_read = adjustl(t_grid(2:8))
              read(t_grid_read,*) var_periods
              if (var_periods > 1) var_sz = var_sz * var_periods
           endif
        endif

        if ((m==start_var) .OR. (size_in_file < max_file_size-var_sz)) then
           size_in_file = size_in_file + var_sz
        else ; exit
        endif

     enddo
     next_var = m

     restartpath = ''
     restartpath2 = ''
     filename_appendix = ''
     restartnameapp = '' 

     !query fms_io if there is a filename_appendix (for ensemble runs)
     call get_filename_appendix(filename_appendix)
     if (len_trim(filename_appendix) > 0 .and. (trim(filename_appendix) /= '')) then
        length = len_trim(restartname)
        if (restartname(length-2:length) == '.nc') then
           restartnameapp = restartname(1:length-3)//'.'//trim(filename_appendix)//'.nc'
        else
           restartnameapp = restartname(1:length)  //'.'//trim(filename_appendix)
        endif
        restartpath = trim(directory)//trim(restartnameapp)
     else
        restartpath = trim(directory)//trim(restartname)
     endif

     ! append '.nc' to the restart file name if it is missing
     substring_index = index('.nc', trim(restartpath))
     if (substring_index <= 0) then
        restartpath2 = append_substring(restartpath,'.nc')
     else
        restartpath2(1:len_trim(restartpath)) = trim(restartpath)
     endif

     if (num_files < 10) then
        write(suffix,'("_",I1)') num_files
     else
        write(suffix,'("_",I2)') num_files
     endif

     if (num_files > 0) restartpath2 = trim(restartpath2) // trim(suffix)
       
     !Prepare the checksum of the restart fields to be written to restart files
     call get_checksum_loop_ranges(G, pos, isL, ieL, jsL, jeL)
     do m=start_var,next_var-1
        if (associated(CS%var_ptr3d(m)%p)) then
           check_val(m-start_var+1,1) = mpp_chksum(CS%var_ptr3d(m)%p(isL:ieL,jsL:jeL,:))
        elseif (associated(CS%var_ptr2d(m)%p)) then
           check_val(m-start_var+1,1) = mpp_chksum(CS%var_ptr2d(m)%p(isL:ieL,jsL:jeL))
        elseif (associated(CS%var_ptr4d(m)%p)) then
           check_val(m-start_var+1,1) = mpp_chksum(CS%var_ptr4d(m)%p(isL:ieL,jsL:jeL,:,:))
        elseif (associated(CS%var_ptr1d(m)%p)) then
           check_val(m-start_var+1,1) = mpp_chksum(CS%var_ptr1d(m)%p)
        elseif (associated(CS%var_ptr0d(m)%p)) then
           check_val(m-start_var+1,1) = mpp_chksum(CS%var_ptr0d(m)%p,pelist=(/mpp_pe()/))
        endif
     enddo
     
     file_open_success = .false.
     file_open_success = fms2_open_file(CS%fileObjWrite, restartpath2,"append", &
                                        G%Domain%mpp_domain, is_restart=.true.)
     if (.not. (file_open_success)) then
        call MOM_error(FATAL,"MOM_restart::save_restart: Failed to open file "//trim(restartpath))
     endif

     call query_vardesc(CS%restart_field(m)%vars, hor_grid=hor_grid, &
                           z_grid=z_grid, t_grid=t_grid, caller="save_restart")
     ! write the axis (dimension) data to the restart file
     axis_exists = fms2_dimension_exists(CS%fileObjWrite,'lath')
     variable_exists = fms2_variable_exists(CS%fileobjWrite, 'lath')
     if (axis_exists .and. .not.(variable_exists)) call write_axis_data(CS%fileObjWrite,'lath', &
                                                                   G=G, is_restart_file = .true.) 
     axis_exists = fms2_dimension_exists(CS%fileObjWrite,'lonh')
     variable_exists = fms2_dimension_exists(CS%fileObjWrite,'lonh')
     if (axis_exists .and. .not.(variable_exists)) call write_axis_data(CS%fileObjWrite,'lonh', &
                                                                   G=G, is_restart_file = .true.)
     axis_exists = fms2_dimension_exists(CS%fileObjWrite,'latq')
     variable_exists = fms2_dimension_exists(CS%fileObjWrite,'latq')
     if (axis_exists .and. .not.(variable_exists)) call write_axis_data(CS%fileObjWrite,'latq', &
                                                                   G=G, is_restart_file = .true.)
     axis_exists = fms2_dimension_exists(CS%fileObjWrite,'lonq')
     variable_exists = fms2_dimension_exists(CS%fileObjWrite,'lonq')
     if (axis_exists .and. .not.(variable_exists)) call write_axis_data(CS%fileObjWrite,'lonq', &
                                                                   G=G, is_restart_file = .true.)
     axis_exists = fms2_dimension_exists(CS%fileObjWrite,'Layer')
     variable_exists = fms2_dimension_exists(CS%fileObjWrite,'Layer')
     if (axis_exists .and. .not.(variable_exists)) call write_axis_data(CS%fileObjWrite,'Layer', &
                                                                   GV=GV, is_restart_file = .true.)
     axis_exists = fms2_dimension_exists(CS%fileObjWrite,'Interface')
     variable_exists = fms2_dimension_exists(CS%fileObjWrite,'Interface')
     if (axis_exists .and. .not.(variable_exists)) call write_axis_data(CS%fileObjWrite,'Interface', &
                                                                   GV=GV, is_restart_file = .true.)
     axis_exists = fms2_dimension_exists(CS%fileObjWrite,'Time')
     variable_exists = fms2_dimension_exists(CS%fileObjWrite,'Time')
     if (axis_exists .and. .not.(variable_exists)) then
        call write_axis_data(CS%fileObjWrite,'Time', &
                             restart_time_in_days=restart_time, is_restart_file = .true.)
     else
        axis_exists = fms2_dimension_exists(CS%fileObjWrite,'Period')
        variable_exists = fms2_dimension_exists(CS%fileObjWrite,'Period')
        if (axis_exists .and. .not.(variable_exists)) call write_axis_data(CS%fileObjWrite, 'Period', & 
                                                                        t_grid_in=t_grid, is_restart_file = .true.)
     endif
     
     do m=start_var,next_var-1     
        if (associated(CS%var_ptr3d(m)%p)) then
           !call write_field(unit,fields(m-start_var+1), G%Domain%mpp_domain, &
           !                 CS%var_ptr3d(m)%p, restart_time)
           call fms2_write_data(CS%fileObjWrite,vars(m)%name, CS%var_ptr3d(m)%p)
        elseif (associated(CS%var_ptr2d(m)%p)) then
           !call write_field(unit,fields(m-start_var+1), G%Domain%mpp_domain, &
           !                 CS%var_ptr2d(m)%p, restart_time)
           call fms2_write_data(CS%fileObjWrite,vars(m)%name, CS%var_ptr2d(m)%p)
        elseif (associated(CS%var_ptr4d(m)%p)) then
           !call write_field(unit,fields(m-start_var+1), G%Domain%mpp_domain, &
           !                 CS%var_ptr4d(m)%p, restart_time)
           call fms2_write_data(CS%fileObjWrite,vars(m)%name, CS%var_ptr4d(m)%p)
        elseif (associated(CS%var_ptr1d(m)%p)) then
           !call write_field(unit, fields(m-start_var+1), CS%var_ptr1d(m)%p, &
           !                 restart_time)
           call fms2_write_data(CS%fileObjWrite,vars(m)%name, CS%var_ptr1d(m)%p)
        elseif (associated(CS%var_ptr0d(m)%p)) then
           !call write_field(unit, fields(m-start_var+1), CS%var_ptr0d(m)%p, &
           !                 restart_time)
           call fms2_write_data(CS%fileObjWrite,vars(m)%name, CS%var_ptr0d(m)%p)
        endif
        ! write the checksum
        call fms2_register_variable_attribute(CS%fileObjWrite,vars(m)%name,'checksum',check_val(m,1))      
     enddo

     !call close_file(unit)
     call fms2_close_file(CS%fileObjWrite)
     num_files = num_files+1
  enddo

end subroutine save_restart

!> restore_state reads the model state from previously generated files.  All
!! restart variables are read from the first file in the input filename list
!! in which they are found.
subroutine restore_state(filename, directory, day, G, CS)
  character(len=*),      intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files.
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(time_type),       intent(out) :: day       !< The time of the restarted run
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  pointer     :: CS        !< The control structure returned by a previous
                                                  !! call to restart_init.

!  This subroutine reads the model state from previously
!  generated files.  All restart variables are read from the first
!  file in the input filename list in which they are found.

  ! Local variables
  character(len=200) :: filepath  ! The path (dir/file) to the file being opened.
  character(len=80) :: fname      ! The name of the current file.
  character(len=8)  :: suffix     ! A suffix (like "_2") that is added to any
                                  ! additional restart files.
  character(len=512) :: mesg      ! A message for warnings.
  character(len=80) :: varname    ! A variable's name.
  integer :: num_file        ! The number of files (restart files and others
                             ! explicitly in filename) that are open.
  integer :: i, n, m, missing_fields, ic
  integer :: isL, ieL, jsL, jeL, is0, js0
  integer :: sizes(7)
  integer :: ndim, nvar, natt, ntime, pos

  integer :: unit(CS%max_fields) ! The mpp unit of all open files.
  character(len=200) :: unit_path(CS%max_fields) ! The file names.
  logical :: unit_is_global(CS%max_fields) ! True if the file is global.

  character(len=8)   :: hor_grid ! Variable grid info.
  real    :: t1, t2 ! Two times.
  real, allocatable :: time_vals(:)
  type(fieldtype), allocatable :: fields(:)
  logical                          :: check_exist, is_there_a_checksum
  integer(LONG_KIND),dimension(3)  :: checksum_file
  integer(kind=8)                  :: checksum_data
  logical :: file_open_success = .false. ! returned by call to fms2_open_file 
  type(MOM_restart_CS) :: fileObjRead  ! fms2 data structure
  integer :: str_split_index = 1
  integer :: str_end_index = 1
  character(len=200) :: file_path_1,file_path_2
  character(len=80), dimension(:), allocatable :: variable_names(:) ! File variable names
  character(len=64) :: checksum_char
  integer(LONG_KIND) :: checksumh
  integer :: num_checksumh, last, is, k
  logical :: var_exists = .false.

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "restore_state: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)
    
! Get number of restart files and full paths to files 
  if ((LEN_TRIM(filename) == 1) .and. (filename(1:1) == 'F')) then
     num_file = open_restart_units('r', directory, G, CS, &
                     file_paths=unit_path, global_files=unit_is_global)
  else
    num_file = open_restart_units(filename, directory, G, CS, &
                     file_paths=unit_path, global_files=unit_is_global)
  endif

  if (num_file == 0) then
     write(mesg,'("Unable to find any restart files specified by  ",A,"  in directory ",A,".")') &
                  trim(filename), trim(directory)
     call MOM_error(FATAL,"MOM_restart: "//mesg)
  endif

! Get the time from the first file in the list that has one.
  do n=1,num_file
     file_path_1(1:len_trim(unit_path(n))) = trim(unit_path(n))
     str_split_index = INDEX(file_path_1,'.res')
     str_end_index = INDEX(file_path_1,'.nc') 
     if (str_split_index > 1 .and. str_end_index > 1 )then ! remove .res. from the file name since fms read automatically appends it to the file name
         file_path_2 = trim(file_path_1(1:str_split_index-1)// &
                         file_path_1(str_split_index+4:str_end_index-1))//'.nc'                         
     else 
         file_path_2 = trim(file_path_1)
     endif
     file_open_success=fms2_open_file(CS%fileObjRead, trim(file_path_2),"read",  G%Domain%mpp_domain, is_restart = .true.)
                                                
     if(.not. file_open_success) then 
        write(mesg,'( "ERROR, unable to open restart file  ",A) ') trim(unit_path(n))
        call MOM_error(FATAL,"MOM_restart: "//mesg)
     endif
     call fms2_get_dimension_size(CS%fileObjRead, "Time", ntime)

     if (ntime < 1) cycle

     allocate(time_vals(ntime))
     call fms2_register_restart_field(CS%fileObjRead, "Time","double")
     call fms2_read_data(CS%fileObjRead,"Time", time_vals)
     t1 = time_vals(1)
     deallocate(time_vals)
     day = real_to_time(t1*86400.0)

     call fms2_close_file(CS%fileObjRead)
     exit
  enddo

  if (n>num_file) call MOM_error(WARNING,"MOM_restart: " // &
                                 "No times found in restart files.")

! Check the remaining files, if any, for different times and issue a warning
! if they differ from the first time.
  if (is_root_pe()) then
     do m = n+1,num_file
        file_path_1(1:len_trim(unit_path(n))) = trim(unit_path(n))
        str_split_index = INDEX(file_path_1,'.res')
        str_end_index = INDEX(file_path_1,'.nc')
        if (str_split_index > 1 .and. str_end_index > 1 )then
           file_path_2 = trim(file_path_1(1:str_split_index-1)// &
                         file_path_1(str_split_index+4:str_end_index-1))//'.nc'                         
        else 
          file_path_2 = trim(file_path_1)
        endif
        file_open_success=fms2_open_file(CS%fileObjRead, trim(file_path_2),"read", G%Domain%mpp_domain, is_restart = .true.)
                                                
        if(.not. file_open_success) then 
           write(mesg,'( "ERROR, unable to open restart file  ",A) ') trim(unit_path(n))
           call MOM_error(FATAL,"MOM_restart: "//mesg)
        endif
        
        call fms2_get_dimension_size(CS%fileObjRead, "Time", ntime)

        if (ntime < 1) cycle

        allocate(time_vals(ntime))
        call fms2_register_restart_field(CS%fileObjRead, "Time","double")
        call fms2_read_data(CS%fileObjRead,"Time", time_vals)
        t2 = time_vals(1)
        deallocate(time_vals)

        if (t1 /= t2) then
           write(mesg,'("WARNING: Restart file ",I2," has time ",F10.4,"whereas simulation is restarted at ",F10.4," (differing by ",F10.4,").")')&
              m,t1,t2,t1-t2
           call MOM_error(WARNING, "MOM_restart: "//mesg)
        endif
        call fms2_close_file(CS%fileObjRead)
     enddo
  endif

  if (n>num_file) call MOM_error(WARNING,"MOM_restart: " // &
                            "No times found in restart files.")

! Read each variable from the first file in which it is found.
  do n=1,num_file
     !call get_file_info(unit(n), ndim, nvar, natt, ntime)
     file_path_1(1:len_trim(unit_path(n))) = trim(unit_path(n))
     str_split_index = INDEX(file_path_1,'.res')
     str_end_index = INDEX(file_path_1,'.nc')
     if (str_split_index > 1 .and. str_end_index > 1 )then
        file_path_2 = trim(file_path_1(1:str_split_index-1)// &
                         file_path_1(str_split_index+4:str_end_index-1))//'.nc'                         
     else 
        file_path_2 = trim(file_path_1)
     endif
     file_open_success=fms2_open_file(CS%fileObjRead, trim(file_path_2),"read", G%Domain%mpp_domain, is_restart = .true.)
                                                
     if (.not. file_open_success) then 
        write(mesg,'( "ERROR, unable to open restart file  ",A )') trim(unit_path(n))
        call MOM_error(FATAL,"MOM_restart: "//mesg)
     endif
     ! register the horizontal axes
     call fms2_register_axis(CS%fileObjRead,'latq','y')
     call fms2_register_axis(CS%fileObjRead,'lath','y')
     call fms2_register_axis(CS%fileObjRead,'lonq','x')
     call fms2_register_axis(CS%fileObjRead,'lonh','x')

     ! get number of variables in the file
     nvar=fms2_get_num_variables(CS%fileObjRead)
     ! get the names of the variables in the file
     allocate(variable_names(nvar))
     call fms2_get_variable_names(CS%fileObjRead, variable_names)
     ! allocate the fields to hold the variable data
     allocate(fields(nvar))
    
     missing_fields = 0

     do m=1,CS%novars
        if (CS%restart_field(m)%initialized) cycle
        call query_vardesc(CS%restart_field(m)%vars, hor_grid=hor_grid, &
                             caller="restore_state")
        select case (hor_grid)
           case ('q') ; pos = CORNER
           case ('h') ; pos = CENTER
           case ('u') ; pos = EAST_FACE
           case ('v') ; pos = NORTH_FACE
           case ('Bu') ; pos = CORNER
           case ('T')  ; pos = CENTER
           case ('Cu') ; pos = EAST_FACE
           case ('Cv') ; pos = NORTH_FACE
           case ('1') ; pos = 0
           case default ; pos = 0
        end select

        call get_checksum_loop_ranges(G, pos, isL, ieL, jsL, jeL)
        do i=1, nvar
           varname = " "
           checksum_char = " "
           var_exists = .false.
           varname(1:len_trim(CS%restart_field(m)%var_name)) = trim(CS%restart_field(m)%var_name)
           ! check if variable is in the restart file
           var_exists = fms2_variable_exists(CS%fileObjRead, varname)

           if (var_exists) then
              if (.NOT. CS%checksum_required) then
                 is_there_a_checksum = .false. ! Do not need to do data checksumming.
              else
                 check_exist = fms2_attribute_exists(CS%fileObjRead, varname, "checksum")
                 checksum_file(:) = -1
                 checksum_data = -1
                 is_there_a_checksum = .false.
                 if ( check_exist ) then
                    ! The following checksum conversion proceure is adapted from mpp_get_atts 
                    call fms2_get_variable_attribute(CS%fileObjRead,varname,"checksum",checksum_char)

                    last = len_trim(checksum_char)
                    is = index(trim(checksum_char),",") ! A value of 0 implies only 1 checksum value
                    ! Scan checksum character array for the ',' delimiter, which indicates that the corresponding variable
                    ! has multiple time levels.
                    checksumh = 0
                    num_checksumh = 1
                    do while ((is > 0) .and. (is < (last-15)))
                       is = is + scan(checksum_char(is:last), "," ) ! move starting pointer after ","
                       num_checksumh = num_checksumh + 1
                    enddo
           
                    is = 1
  
                    do k = 1, num_checksumh
                       read(checksum_char(is:is+15),'(Z16)') checksumh ! Z=hexadecimal integer: Z16 is for 64-bit data types
                       checksum_file(k) = checksumh 
                       is = is + 17 ! Move index past the ',' in checksum_char
                    enddo

                    is_there_a_checksum = .true.
                 endif
              endif
              ! register the restart variable
              call fms2_register_restart_field(CS%fileObjRead,varname,'real')

              if (associated(CS%var_ptr1d(m)%p))  then
                 call fms2_read_data(CS%fileObjRead, varname, CS%var_ptr1d(m)%p)

                 if (is_there_a_checksum) checksum_data = mpp_chksum(CS%var_ptr1d(m)%p)
              elseif (associated(CS%var_ptr0d(m)%p)) then ! Read a scalar
                 call fms2_read_data(CS%fileObjRead,varname,  CS%var_ptr0d(m)%p)
                 if (is_there_a_checksum) checksum_data = mpp_chksum(CS%var_ptr0d(m)%p,pelist=(/mpp_pe()/))
                   
              elseif (associated(CS%var_ptr2d(m)%p)) then  ! Read a 2d array
                 if (pos /= 0) then
                    call fms2_read_data(CS%fileObjRead, varname, CS%var_ptr2d(m)%p) ! domain-decomposed read
                 !else ! This array is not domain-decomposed.  This variant may be under-tested.
                 !   call read_data(unit_path(n), varname, CS%var_ptr2d(m)%p, &
                 !          no_domain=.true., timelevel=1)
                 endif
                 if (is_there_a_checksum) checksum_data = mpp_chksum(CS%var_ptr2d(m)%p(isL:ieL,jsL:jeL))
              elseif (associated(CS%var_ptr3d(m)%p)) then  ! Read a 3d array.
                 if (pos /= 0) then
                    call fms2_read_data(CS%fileObjRead, varname,  CS%var_ptr3d(m)%p) ! domain-decomposed read

                 !else ! This array is not domain-decomposed.  This variant may be under-tested.
                 !   call read_data(unit_path(n), varname, CS%var_ptr3d(m)%p, &
                 !                  no_domain=.true., timelevel=1)
                 endif
                 if (is_there_a_checksum) checksum_data = mpp_chksum(CS%var_ptr3d(m)%p(isL:ieL,jsL:jeL,:))
              elseif (associated(CS%var_ptr4d(m)%p)) then  ! Read a 4d array.
                 if (pos /= 0) then
                     call fms2_read_data(CS%fileObjRead, varname, CS%var_ptr4d(m)%p) ! domain-decomposed read
                 !else ! This array is not domain-decomposed.  This variant may be under-tested.
                 !    call read_data(unit_path(n), varname, CS%var_ptr4d(m)%p, &
                 !           no_domain=.true., timelevel=1)
                 endif

                 if (is_there_a_checksum) then 
                    checksum_data = mpp_chksum(CS%var_ptr4d(m)%p(isL:ieL,jsL:jeL,:,:))
                 endif                 
              else
                 call MOM_error(FATAL, "MOM_restart restore_state: No pointers set for "//trim(varname))
              endif

              if (is_root_pe() .and. is_there_a_checksum .and. (checksum_file(1) /= checksum_data)) then
                 write (mesg,'(a,Z16,a,Z16,a)') "Checksum of input field "// trim(varname)//" ",checksum_data,&
                        " does not match value ", checksum_file(1), &
                        " stored in "//trim(unit_path(n)//"." )
                 call MOM_error(FATAL, "MOM_restart(restore_state): "//trim(mesg) )
              endif

              CS%restart_field(m)%initialized = .true.
               
              exit ! Start search for next restart variable.
           endif
        enddo
        if (i>nvar) missing_fields = missing_fields+1
     enddo

     deallocate(fields)
     deallocate(variable_names)

     call fms2_close_file(CS%fileObjRead)

     if (missing_fields == 0) exit   
  enddo

! Check whether any mandatory fields have not been found.
  CS%restart = .true.
  do m=1,CS%novars
     if (.not.(CS%restart_field(m)%initialized)) then
        CS%restart = .false.
        if (CS%restart_field(m)%mand_var) then
           call MOM_error(FATAL,"MOM_restart: Unable to find mandatory variable " &
                       //trim(CS%restart_field(m)%var_name)//" in restart files.")
        endif
     endif
  enddo

end subroutine restore_state

!> restart_files_exist determines whether any restart files exist.
function restart_files_exist(filename, directory, G, CS)
  character(len=*),      intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files.
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  pointer     :: CS        !< The control structure returned by a previous
                                                  !! call to restart_init.
  logical :: restart_files_exist                  !< The function result, which indicates whether
                                                  !! any of the explicitly or automatically named
                                                  !! restart files exist in directory.
  integer :: num_files

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "restart_files_exist: Module must be initialized before it is used.")

  if ((LEN_TRIM(filename) == 1) .and. (filename(1:1) == 'F')) then
    num_files = open_restart_units('r', directory, G, CS)
  else
    num_files = open_restart_units(filename, directory, G, CS)
  endif
  restart_files_exist = (num_files > 0)

end function restart_files_exist

!> determine_is_new_run determines from the value of filename and the existence
!! automatically named restart files in directory whether this would be a new,
!! and as a side effect stores this information in CS.
function determine_is_new_run(filename, directory, G, CS) result(is_new_run)
  character(len=*),      intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files.
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  pointer     :: CS        !< The control structure returned by a previous
                                                  !! call to restart_init.
  logical :: is_new_run                           !< The function result, which indicates whether
                                                  !! this is a new run, based on the value of
                                                  !! filename and whether restart files exist.

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "determine_is_new_run: Module must be initialized before it is used.")
  if (LEN_TRIM(filename) > 1) then
    CS%new_run = .false.
  elseif (LEN_TRIM(filename) == 0) then
    CS%new_run = .true.
  elseif (filename(1:1) == 'n') then
    CS%new_run = .true.
  elseif (filename(1:1) == 'F') then
    CS%new_run = (open_restart_units('r', directory, G, CS) == 0)
  else
    CS%new_run = .false.
  endif

  CS%new_run_set = .true.
  is_new_run = CS%new_run
end function determine_is_new_run

!> is_new_run returns whether this is going to be a new run based on the
!! information stored in CS by a previous call to determine_is_new_run.
function is_new_run(CS)
  type(MOM_restart_CS),  pointer :: CS !< The control structure returned by a previous
                                       !! call to restart_init.
  logical :: is_new_run                !< The function result, which indicates whether
                                       !! this is a new run, based on the value of
                                       !! filename and whether restart files exist.

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "is_new_run: Module must be initialized before it is used.")
  if (.not.CS%new_run_set) call MOM_error(FATAL, "MOM_restart " // &
      "determine_is_new_run must be called for a restart file before is_new_run.")

  is_new_run = CS%new_run
end function is_new_run

!> open_restart_units determines the number of existing restart files and optionally opens
!! them and returns unit ids, paths and whether the files are global or spatially decomposed.
function open_restart_units(filename, directory, G, CS, units, file_paths, &
                            global_files) result(num_files)
  character(len=*),      intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files.
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  pointer     :: CS        !< The control structure returned by a previous
                                                  !! call to restart_init.
  integer, dimension(:), &
               optional, intent(out) :: units     !< The mpp units of all opened files.
  character(len=*), dimension(:), &
               optional, intent(out) :: file_paths   !< The full paths to open files.
  logical, dimension(:), &
               optional, intent(out) :: global_files !< True if a file is global.

  integer :: num_files  !< The number of files (both automatically named restart
                        !! files and others explicitly in filename) that have been opened.

!  This subroutine reads the model state from previously
!  generated files.  All restart variables are read from the first
!  file in the input filename list in which they are found.

  ! Local variables
  character(len=256) :: filepath  ! The path (dir/file) to the file being opened.
  character(len=256) :: fname     ! The name of the current file.
  character(len=8)   :: suffix    ! A suffix (like "_2") that is added to any
                                  ! additional restart files.
! character(len=256) :: mesg      ! A message for warnings.
  integer :: num_restart     ! The number of restart files that have already
                             ! been opened.
  integer :: start_char      ! The location of the starting character in the
                             ! current file name.
  integer :: n, m, err, length


  logical :: fexists
  character(len=32) :: filename_appendix = '' !fms appendix to filename for ensemble runs
  character(len=80) :: restartname

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "open_restart_units: Module must be initialized before it is used.")

! Get NetCDF ids for all of the restart files.
  num_restart = 0 ; n = 1 ; start_char = 1
  do while (start_char <= len_trim(filename) )
    do m=start_char,len_trim(filename)
      if (filename(m:m) == ' ') exit
    enddo
    fname = filename(start_char:m-1)
    start_char = m
    do while (start_char <= len_trim(filename))
      if (filename(start_char:start_char) == ' ') then
        start_char = start_char + 1
      else
        exit
      endif
    enddo

    if ((fname(1:1)=='r') .and. ( len_trim(fname) == 1)) then
      err = 0
      if (num_restart > 0) err = 1 ! Avoid going through the file list twice.
      do while (err == 0)
        restartname = trim(CS%restartfile)

       !query fms_io if there is a filename_appendix (for ensemble runs)
       call get_filename_appendix(filename_appendix)
       if (len_trim(filename_appendix) > 0) then
         length = len_trim(restartname)
         if (restartname(length-2:length) == '.nc') then
         restartname = restartname(1:length-3)//'.'//trim(filename_appendix)//'.nc'
         else
           restartname = restartname(1:length)  //'.'//trim(filename_appendix)
         endif
        endif
        filepath = trim(directory) // trim(restartname)

        if (num_restart < 10) then
          write(suffix,'("_",I1)') num_restart
        else
          write(suffix,'("_",I2)') num_restart
        endif
        if (num_restart > 0) filepath = trim(filepath) // suffix

        ! if (.not.file_exists(filepath)) &
          filepath = trim(filepath)//".nc"

        num_restart = num_restart + 1
        inquire(file=filepath, exist=fexists)
        if (fexists) then
          if (present(units)) &
            call open_file(units(n), trim(filepath), READONLY_FILE, NETCDF_FILE, &
                           threading = MULTIPLE, fileset = SINGLE_FILE)
          if (present(global_files)) global_files(n) = .true.
        elseif (CS%parallel_restartfiles) then
          ! Look for decomposed files using the I/O Layout.
          !fexists = file_exists(filepath, G%Domain)
          fexists = fms2_file_exists(filepath)
          if (fexists .and. (present(units))) &
            call open_file(units(n), trim(filepath), READONLY_FILE, NETCDF_FILE, &
                           domain=G%Domain%mpp_domain)
          if (fexists .and. present(global_files)) global_files(n) = .false.
        endif

        if (fexists) then
          if (present(file_paths)) file_paths(n) = filepath
          n = n + 1
          if (is_root_pe() .and. (present(units))) &
            call MOM_error(NOTE, "MOM_restart: MOM run restarted using : "//trim(filepath))
        else
          err = 1 ; exit
        endif
      enddo ! while (err == 0) loop
    else
      filepath = trim(directory)//trim(fname)
      inquire(file=filepath, exist=fexists)
      if (.not. fexists) filepath = trim(filepath)//".nc"

      inquire(file=filepath, exist=fexists)
      if (fexists) then
        if (present(units)) &
          call open_file(units(n), trim(filepath), READONLY_FILE, NETCDF_FILE, &
                       threading = MULTIPLE, fileset = SINGLE_FILE)
        if (present(global_files)) global_files(n) = .true.
        if (present(file_paths)) file_paths(n) = filepath
        n = n + 1
        if (is_root_pe() .and. (present(units))) &
          call MOM_error(NOTE,"MOM_restart: MOM run restarted using : "//trim(filepath))
      else
        if (present(units)) &
          call MOM_error(WARNING,"MOM_restart: Unable to find restart file : "//trim(filepath))
      endif

    endif
  enddo ! while (start_char < strlen(filename)) loop
  num_files = n-1

end function open_restart_units

!> Initialize this module and set up a restart control structure.
subroutine restart_init(param_file, CS, restart_root)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
  type(MOM_restart_CS),  pointer    :: CS !< A pointer to a MOM_restart_CS object that is allocated here
  character(len=*), optional, &
                         intent(in) :: restart_root !< A filename root that overrides the value
                                          !! set by RESTARTFILE to enable the use of this module by
                                          !! other components than MOM.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_restart"   ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "restart_init called with an associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "PARALLEL_RESTARTFILES", &
                                CS%parallel_restartfiles, &
                 "If true, each processor writes its own restart file, \n"//&
                 "otherwise a single restart file is generated", &
                 default=.false.)

  if (present(restart_root)) then
    CS%restartfile = restart_root
    call log_param(param_file, mdl, "RESTARTFILE from argument", CS%restartfile)
  else
    call get_param(param_file, mdl, "RESTARTFILE", CS%restartfile, &
                 "The name-root of the restart file.", default="MOM.res")
  endif
  call get_param(param_file, mdl, "LARGE_FILE_SUPPORT", CS%large_file_support, &
                 "If true, use the file-size limits with NetCDF large \n"//&
                 "file support (4Gb), otherwise the limit is 2Gb.", &
                 default=.true.)
  call get_param(param_file, mdl, "MAX_FIELDS", CS%max_fields, &
                 "The maximum number of restart fields that can be used.", &
                 default=100)
  call get_param(param_file, mdl, "RESTART_CHECKSUMS_REQUIRED", CS%checksum_required, &
                 "If true, require the restart checksums to match and error out otherwise. \n"//&
                 "Users may want to avoid this comparison if for example the restarts are  \n"//&
                 "made from a run with a different mask_table than the current run,  \n"//&
                 "in which case the checksums will not match and cause crash.",&
                 default=.true.)

  allocate(CS%restart_field(CS%max_fields))
  allocate(CS%var_ptr0d(CS%max_fields))
  allocate(CS%var_ptr1d(CS%max_fields))
  allocate(CS%var_ptr2d(CS%max_fields))
  allocate(CS%var_ptr3d(CS%max_fields))
  allocate(CS%var_ptr4d(CS%max_fields))

end subroutine restart_init

!> Indicate that all variables have now been registered.
subroutine restart_init_end(CS)
  type(MOM_restart_CS),  pointer    :: CS !< A pointer to a MOM_restart_CS object

  if (associated(CS)) then
    if (CS%novars == 0) call restart_end(CS)
  endif

end subroutine restart_init_end

!> Deallocate memory associated with a MOM_restart_CS variable.
subroutine restart_end(CS)
  type(MOM_restart_CS),  pointer    :: CS !< A pointer to a MOM_restart_CS object

  if (associated(CS%restart_field)) deallocate(CS%restart_field)
  if (associated(CS%var_ptr0d)) deallocate(CS%var_ptr0d)
  if (associated(CS%var_ptr1d)) deallocate(CS%var_ptr1d)
  if (associated(CS%var_ptr2d)) deallocate(CS%var_ptr2d)
  if (associated(CS%var_ptr3d)) deallocate(CS%var_ptr3d)
  if (associated(CS%var_ptr4d)) deallocate(CS%var_ptr4d)
  deallocate(CS)

end subroutine restart_end

subroutine restart_error(CS)
  type(MOM_restart_CS),  pointer    :: CS !< A pointer to a MOM_restart_CS object

  character(len=16)  :: num  ! String for error messages

  if (CS%novars > CS%max_fields) then
    write(num,'(I0)') CS%novars
    call MOM_error(FATAL,"MOM_restart: Too many fields registered for " // &
           "restart.  Set MAX_FIELDS to be at least " // &
           trim(adjustl(num)) // " in the MOM input file.")
  else
    call MOM_error(FATAL,"MOM_restart: Unspecified fatal error.")
  endif
end subroutine restart_error

!> Return bounds for computing checksums to store in restart files
subroutine get_checksum_loop_ranges(G, pos, isL, ieL, jsL, jeL)
  type(ocean_grid_type), intent(in)  :: G !< The ocean's grid structure
  integer,               intent(in)  :: pos !< An integer indicating staggering of variable
  integer,               intent(out) :: isL !< i-start for checksum
  integer,               intent(out) :: ieL !< i-end for checksum
  integer,               intent(out) :: jsL !< j-start for checksum
  integer,               intent(out) :: jeL !< j-end for checksum

  ! Regular non-symmetric compute domain
  isL = G%isc-G%isd+1
  ieL = G%iec-G%isd+1
  jsL = G%jsc-G%jsd+1
  jeL = G%jec-G%jsd+1

  ! Expand range east or south for symmetric arrays
  if (G%symmetric) then
    if ((pos == EAST_FACE) .or. (pos == CORNER)) then ! For u-, q-points only
      if (G%idg_offset == 0) isL = isL - 1 ! include western edge in checksums only for western PEs
    endif
    if ((pos == NORTH_FACE) .or. (pos == CORNER)) then ! For v-, q-points only
      if (G%jdg_offset == 0) jsL = jsL - 1 ! include western edge in checksums only for southern PEs
    endif
  endif

end subroutine get_checksum_loop_ranges

!> check restart file for an axis, and register it if it is unregistered

subroutine check_for_restart_axis(fileObjWrite, axis_name, axis_length)
   type(FmsNetcdfDomainFile_t), intent(inout) :: fileObjWrite !< file object returned by prior call to fms2_open_file
   character(len=*), intent(in) :: axis_name ! name of the restart file axis to register to file
   integer, intent(in) :: axis_length ! length of axis/dimension (only needed for Z and Time)
   ! local  
   logical :: axis_exists 
   
   axis_exists = .false.
 
   axis_exists = fms2_dimension_exists(fileObjWrite, axis_name)
   if (.not. (axis_exists)) then
      select case (trim(axis_name))
         case ('latq'); call fms2_register_axis(fileObjWrite,'latq','y')
         case ('lath'); call fms2_register_axis(fileObjWrite,'lath','y') 
         case ('lonq'); call fms2_register_axis(fileObjWrite,'lonq','x') 
         case ('lonh'); call fms2_register_axis(fileObjWrite,'lonh','x')
         case ('Layer'); call fms2_register_axis(fileObjWrite,'Layer',axis_length)
         case ('Interface'); call fms2_register_axis(fileObjWrite,'Interface',axis_length)
         case ('Time'); call fms2_register_axis(fileObjWrite,'Time', axis_length)
         case ('Period'); call fms2_register_axis(fileObjWrite,'Period',axis_length)
      end select
   endif
  
end subroutine check_for_restart_axis

end module MOM_restart
