!> The MOM6 facility for reading and writing restart files, and querying what has been read.
module MOM_restart

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domains, only : pe_here, num_PEs, MOM_domain_type
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_string_functions, only : lowercase, append_substring
use MOM_grid, only : ocean_grid_type
use MOM_io, only : fieldtype,open_file, close_file
use MOM_io, only : MOM_read_data
use MOM_io, only : get_file_info, get_file_atts, get_file_fields, get_file_times
use MOM_io, only : vardesc, var_desc, query_vardesc, modify_vardesc
use MOM_io, only : MULTIPLE, NETCDF_FILE, READONLY_FILE, SINGLE_FILE
use MOM_io, only : CENTER, CORNER, NORTH_FACE, EAST_FACE
use MOM_io, only : file_exists
use MOM_io, only : get_var_dimension_features
use MOM_io, only : get_horizontal_grid_position
use MOM_io, only : get_time_units
use MOM_io, only : get_variable_byte_size
use MOM_io, only : MOM_register_axis
use MOM_io, only : MOM_register_variable_attribute
use MOM_io, only : MOM_open_file
use MOM_io, only : MOM_write_data
use MOM_io, only : axis_data_type
use MOM_io, only : MOM_get_axis_data

use MOM_time_manager, only : time_type, time_type_to_real, real_to_time
use MOM_time_manager, only : days_in_month, get_date, set_date
use MOM_verticalGrid, only : verticalGrid_type
use mpp_mod,         only:  mpp_chksum, mpp_pe, mpp_max
use fms2_io_mod,     only: fms2_register_restart_field => register_restart_field, &
                           fms2_register_axis => register_axis, &
                           fms2_read_data => read_data, &
                           fms2_read_restart => read_restart, &
                           fms2_write_restart => write_restart,&
                           fms2_open_file => open_file, &
                           fms2_close_file => close_file, &
                           fms2_global_att_exists => global_att_exists, &
                           fms2_attribute_exists => variable_att_exists, &
                           fms2_get_global_attribute => get_global_attribute, &
                           fms2_get_compute_domain_dimension_indices => get_compute_domain_dimension_indices, &
                           fms2_get_global_io_domain_indices => get_global_io_domain_indices, &
                           fms2_get_compute_domain_dimension_indices => get_compute_domain_dimension_indices, &
                           fms2_get_variable_attribute => get_variable_attribute, &
                           fms2_get_variable_names => get_variable_names, &
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
public register_restart_field_as_obsolete

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

!> A structure to store information about restart fields that are no longer used
type obsolete_restart
   character(len=32) :: field_name       !< Name of restart field that is no longer in use
   character(len=32) :: replacement_name !< Name of replacement restart field, if applicable
end type obsolete_restart

!> A restart registry and the control structure for restarts
type, public :: MOM_restart_CS ; private
  logical :: restart    !< restart is set to .true. if the run has been started from a full restart
                        !! file.  Otherwise some fields must be initialized approximately.
  integer :: novars = 0 !< The number of restart fields that have been registered.
  integer :: num_obsolete_vars = 0  !< The number of obsolete restart fields that have been registered.
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

  !> An array of descriptions of the registered fields
  type(field_restart), pointer :: restart_field(:) => NULL()

  !> An array of obsolete restart fields
  type(obsolete_restart), pointer :: restart_obsolete(:) => NULL()

  !>@{ Pointers to the fields that have been registered for restarts
  type(p0d), pointer :: var_ptr0d(:) => NULL()
  type(p1d), pointer :: var_ptr1d(:) => NULL()
  type(p2d), pointer :: var_ptr2d(:) => NULL()
  type(p3d), pointer :: var_ptr3d(:) => NULL()
  type(p4d), pointer :: var_ptr4d(:) => NULL()
  !!@}
  integer :: max_fields !< The maximum number of restart fields
  
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
!!> Register a restart field as obsolete
subroutine register_restart_field_as_obsolete(field_name, replacement_name, CS)
  character(*), intent(in) :: field_name       !< Name of restart field that is no longer in use
  character(*), intent(in) :: replacement_name !< Name of replacement restart field, if applicable
  type(MOM_restart_CS), pointer :: CS          !< A pointer to a MOM_restart_CS object (intent in/out)

  CS%num_obsolete_vars = CS%num_obsolete_vars+1
  CS%restart_obsolete(CS%num_obsolete_vars)%field_name = field_name
  CS%restart_obsolete(CS%num_obsolete_vars)%replacement_name = replacement_name
end subroutine register_restart_field_as_obsolete

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
subroutine register_restart_field_4d(f_ptr, name, mandatory, CS, G, & 
                                     longname, units, hor_grid, z_grid, t_grid)
  real, dimension(:,:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  type(ocean_grid_type), optional,     intent(in) :: G         !< The ocean's grid structure
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizonal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent
 
  ! local
  type(vardesc) :: vd
  character(len=200) :: mesg
 
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: "//&
      "register_restart_field_4d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  vd = var_desc(name, units=units, longname=longname, hor_grid=hor_grid, &
                z_grid=z_grid, t_grid=t_grid)

  call register_restart_field_ptr4d(f_ptr, vd, mandatory, CS)

end subroutine register_restart_field_4d

!> Register a 3-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_3d(f_ptr, name, mandatory, CS, G, & 
                                     longname, units, hor_grid, z_grid, t_grid)
  real, dimension(:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  type(ocean_grid_type),      intent(in) :: G         !< The ocean's grid structure
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizonal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent

  ! local
  type(vardesc) :: vd
  character(len=200) :: mesg
        
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: "//&
      "register_restart_field_3d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  vd = var_desc(name, units=units, longname=longname, hor_grid=hor_grid, &
                z_grid=z_grid, t_grid=t_grid)

  ! register the restart field
  call register_restart_field_ptr3d(f_ptr, vd, mandatory, CS)

end subroutine register_restart_field_3d

!> Register a 2-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_2d(f_ptr, name, mandatory, CS, G, & 
                                     longname, units, hor_grid, z_grid, t_grid)
  real, dimension(:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  type(ocean_grid_type),      intent(in) :: G         !< The ocean's grid structure
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizonal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent
  ! local
  type(vardesc) :: vd
  character(len=8) :: Zgrid
  character(len=200) :: mesg
          
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: "//&
      "register_restart_field_2d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  if (present(z_grid)) then 
     Zgrid = z_grid
  else
     Zgrid = '1' ; 
  endif

  vd = var_desc(name, units=units, longname=longname, hor_grid=hor_grid, &
                z_grid=Zgrid, t_grid=t_grid)

  call register_restart_field_ptr2d(f_ptr, vd, mandatory, CS)

end subroutine register_restart_field_2d

!> Register a 1-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_1d(f_ptr, name, mandatory, CS, G, & 
                                     longname, units, hor_grid, z_grid, t_grid)
  real, dimension(:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  type(ocean_grid_type),      intent(in) :: G         !< The ocean's grid structure
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizonal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent

  ! local 
  type(vardesc) :: vd
  character(len=8) :: hgrid
  character(len=200) :: mesg
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: " // &
      "register_restart_field_3d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  if (present(hor_grid)) then 
     hgrid = hor_grid
  else 
     hgrid = '1'
  endif 

  vd = var_desc(name, units=units, longname=longname, hor_grid=hgrid, &
                z_grid=z_grid, t_grid=t_grid)
 
  call register_restart_field_ptr1d(f_ptr, vd, mandatory, CS)

end subroutine register_restart_field_1d

!> Register a 0-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_0d(f_ptr, name, mandatory, CS, G, & 
                                     longname, units, t_grid)
  real,               target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  type(ocean_grid_type),      intent(in) :: G         !< The ocean's grid structure
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent

  !local
  type(vardesc) :: vd
  character(len=200) :: mesg

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: " // &
      "register_restart_field_0d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  vd = var_desc(name, units=units, longname=longname, hor_grid='1', &
                z_grid='1', t_grid=t_grid)

  call register_restart_field_ptr0d(f_ptr, vd, mandatory, CS)
       
end subroutine register_restart_field_0d

!> query_initialized_name determines whether a named field has been successfully
!! read from a restart file yet.
function query_initialized_name(name, CS) result(query_initialized)
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!  This subroutine returns .true. if the field referred to by name has
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
  type(FmsNetcdfDomainFile_t) :: fileObjWrite  ! file object returned by a call to fms2_open_file
  character(len=1024) :: restartpath     ! The restart file path (dir/file).
  character(len=1024) :: restartname     ! The restart file name (no dir).
  character(len=1024) :: restartpath_temp ! temporary location for the restart file path (dir/file).
  character(len=1024) :: restartname_temp ! temporary location for restart name
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
  integer :: m, nz, i, k, num_files, pos
  integer :: seconds, days, year, month, hour, minute
  character(len=8) :: hor_grid, z_grid, t_grid ! Variable grid info.
  character(len=64) :: var_name         ! A variable's name.
  character(len=256) :: date_appendix   ! date string to append to a file name if desired
  character(len=64) :: axis_names(4)    ! Array to hold up to 4 strings for the variable axis names 
  integer, dimension(4) :: axis_lengths ! Array of integer lengths corresponding to the name(s) in axis_names
  integer :: name_length
  integer(kind=8) :: check_val(CS%max_fields,1)
  integer :: isL, ieL, jsL, jeL
  integer :: is, ie
  integer :: substring_index
  integer :: horgrid_position = 1
  integer :: num_axes, total_axes
  integer :: var_periods
  logical :: file_open_success = .false.
  logical :: axis_exists = .false.
  logical :: variable_exists = .false.
  real :: restart_time
  character(len=16) :: restart_time_units
  character(len=64) :: checksum_char
  character(len=64) :: units
  character(len=256) :: longname
  character(len=16) :: t_grid_read, t_grid_str
  real, dimension(:), allocatable :: time_vals
  real, dimension(:), allocatable :: data_temp
  type(axis_data_type) :: axis_data_CS

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "save_restart: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  ! With parallel read & write, it is possible to disable the following...

  ! The maximum file size is 4294967292, according to the NetCDF documentation.
  if (CS%large_file_support) max_file_size = 4294967292_8

  name_length = 0
  num_files = 0
  restartname = ''
  base_file_name = ''
  restartname_temp = ''
  date_appendix = ''
  restart_time_units = ''

  ! get the number of vertical levels
  nz = 1 ; if (present(GV)) nz = GV%ke
  
  if (present(filename)) then 
     base_file_name = trim(filename) 
  else
     base_file_name=trim(CS%restartfile)
  endif

  ! append a time stamp to the file name if time_stamp is specified
  if (PRESENT(time_stamped)) then
     if (time_stamped) then
        call get_date(time,year,month,days,hour,minute,seconds)
        ! Compute the year-day, because I don't like months. - RWH
        do m=1,month-1
          days = days + days_in_month(set_date(year,m,2,0,0,0))
        enddo
        seconds = seconds + 60*minute + 3600*hour
        if (year <= 9999) then
           write(date_appendix,'("_Y",I4.4,"_D",I3.3,"_S",I5.5)') year, days, seconds
        elseif (year <= 99999) then
           write(date_appendix,'("_Y",I5.5,"_D",I3.3,"_S",I5.5)') year, days, seconds
        else
           write(date_appendix,'("_Y",I10.10,"_D",I3.3,"_S",I5.5)') year, days, seconds
        endif
        restartname_temp = trim(base_file_name)//trim(date_appendix)
     endif 
  else
     restartname_temp = trim(base_file_name)
  endif

  ! append '.nc' to the restart file name if it is missing
  substring_index = index('.nc', trim(restartname_temp))
  if (substring_index <= 0) then
     restartname = append_substring(restartname_temp,'.nc')
  else 
     restartname = restartname_temp
  endif

 ! get the restart time units

  restart_time = time_type_to_real(time) / 86400.0
  restart_time_units = get_time_units(restart_time*86400.0)

  !WRITE(mpp_pe()+2000, '(A,F8.2)') "save_restart: the restart_time is ", restart_time
  !call flush(mpp_pe()+2000)

  next_var = 1
  do while (next_var <= CS%novars )
     start_var = next_var
     restartpath = ''
     restartpath_temp = ''
     suffix = ''

     name_length = len_trim(trim(directory)//trim(restartname))
     restartpath_temp(1:name_length) = trim(directory)//trim(restartname)
     if (num_files < 10) then
        write(suffix,'("_",I1)') num_files
     else
        write(suffix,'("_",I2)') num_files
     endif
  
     if (num_files > 0) then
        name_length = len_trim(trim(restartpath_temp) // trim(suffix))
        restartpath(1:name_length) = trim(restartpath_temp) // trim(suffix)
     else
        name_length = len_trim(restartpath_temp)
        restartpath(1:name_length) = trim(restartpath_temp)
     endif 

     file_open_success = MOM_open_file(fileObjWrite, trim(restartpath),"write", &
                    G, is_restart = .true.)
   
     if (.not. (file_open_success)) then
        call MOM_error(FATAL,"MOM_restart::save_restart: Failed to open file "//trim(restartpath))
     endif

     ! get variable sizes in bytes
     size_in_file = 8*(2*G%Domain%niglobal+2*G%Domain%njglobal+2*nz+1000)

     total_axes=0 ! total number of axes registered and written to restart file
     do m=start_var,CS%novars
        !WRITE(mpp_pe()+2000,*) "save_restart: getting axis stuff for ", trim(CS%restart_field(m)%var_name)
        !call flush(mpp_pe()+2000)
   
        call query_vardesc(CS%restart_field(m)%vars, hor_grid=hor_grid, &
                           z_grid=z_grid, t_grid=t_grid, caller="save_restart")

        ! get the restart time value(s) and assign to time_vals array
        if (.not.(allocated(time_vals))) then
           t_grid_str = ''
           t_grid_str = adjustl(t_grid)
           select case (t_grid_str(1:1))
              case ('s', 'a', 'm') ! allocate an empty array that will be populated after the function call
                 allocate(time_vals(1))
                 time_vals(1) = restart_time
              case ('p')
                 if (len_trim(t_grid(2:8)) > 0) then
                    var_periods = -1
                    t_grid_read = adjustl(t_grid(2:8))
                    read(t_grid_read,*) var_periods
                    if (var_periods < 1) then 
                       call MOM_error(FATAL, "MOM_restart::save_restart: "//&
                                       "Period value must be positive.")
                    endif
                    ! Define a periodic axis array
                    allocate(time_vals(var_periods))
                    do k=1,var_periods
                       time_vals(k) = real(k)
                    enddo
                 endif
           end select
        endif
        
        var_sz = get_variable_byte_size(hor_grid, z_grid, t_grid, G, nz)
     
        if ((m==start_var) .OR. (size_in_file < max_file_size-var_sz)) then
           size_in_file = size_in_file + var_sz
        !else ; exit
        endif
        
        ! get the axis (dimension) names and lengths for variable 'm'                                
        ! note: 4d variables are lon x lat x vertical level x time
        num_axes = 0

        call get_var_dimension_features(hor_grid, z_grid, t_grid, G, &
                                     axis_names, axis_lengths, num_axes, GV)
        
        ! register all of the restart variable axes to the file if they do not exist
        if (num_axes <= 0) then
           call MOM_error(FATAL,"MOM_restart::save_restart: num_axes is an invalid value.")
        endif
        
        do i=1,num_axes
           axis_exists = fms2_dimension_exists(fileObjWrite, axis_names(i))
           if (.not.(axis_exists)) then
              total_axes=total_axes+1
              call MOM_get_axis_data(axis_data_CS, axis_names(i), total_axes, G, GV, &
                                     time_vals, restart_time_units)
              call MOM_register_axis(fileObjWrite, axis_data_CS%name(total_axes), axis_lengths(i))
           endif
        enddo
    
     enddo

     ! register the axis variables
     do i=1,total_axes
        variable_exists = fms2_variable_exists(fileObjWrite, trim(axis_data_CS%name(i)))
        if (.not.(variable_exists)) then 
           if (associated(axis_data_CS%data(i)%p)) then
                
              !allocate(data_temp(size(axis_data_CS%data(i)%p)))
              !data_temp = axis_data_CS%data(i)%p
                 
              if (axis_data_CS%is_domain_decomposed(i)) then
                 call fms2_get_global_io_domain_indices(fileObjWrite, trim(axis_data_CS%name(i)), is, ie)
                 call fms2_register_restart_field(fileObjWrite, axis_data_CS%name(i), axis_data_CS%data(i)%p(is:ie), &
                                                  dimensions=(/trim(axis_data_CS%name(i))/), &
                                                  domain_position=axis_data_CS%horgrid_position(i))

                 !call MOM_write_data(fileObjWrite, trim(axis_data_CS%name), data_temp(is:ie))
              else
                 call fms2_register_restart_field(fileObjWrite, axis_data_CS%name(i), axis_data_CS%data(i)%p, &
                                                  dimensions=(/trim(axis_data_CS%name(i))/), &
                                                  domain_position=axis_data_CS%horgrid_position(i))

                    !call MOM_write_data(fileObjWrite, trim(axis_data_CS%name), data_temp) 
              endif
 
              !deallocate(data_temp)

              call MOM_register_variable_attribute(fileObjWrite, trim(axis_data_CS%name(i)), &
                                                      'long_name',axis_data_CS%longname(i))
              call MOM_register_variable_attribute(fileObjWrite, trim(axis_data_CS%name(i)), &
                                                      'units',axis_data_CS%units(i))
              call MOM_register_variable_attribute(fileObjWrite, trim(axis_data_CS%name(i)), &
                                                      'cartesian_axis',axis_data_CS%cartesian_axis(i))
              if (len_trim(axis_data_CS%positive(i))>1) then
                  call MOM_register_variable_attribute(fileObjWrite, trim(axis_data_CS%name(i)), &
                                                       'positive',axis_data_CS%positive(i))
              endif
                 
           endif
        endif
     enddo   

     next_var = m
     
     do m=start_var,next_var-1
        units=''
        longname=''
        call query_vardesc(CS%restart_field(m)%vars, hor_grid=hor_grid, &
                           z_grid=z_grid, t_grid=t_grid, longname=longname, &
                           units=units, caller="save_restart")
        horgrid_position = get_horizontal_grid_position(hor_grid)  
        
        call get_checksum_loop_ranges(G, horgrid_position, isL, ieL, jsL, jeL)
        
        num_axes = 0

        call get_var_dimension_features(hor_grid, z_grid, t_grid, G, &
                                    axis_names, axis_lengths, num_axes, GV)
        
        ! register and write the restart variables to the file
        if (associated(CS%var_ptr3d(m)%p)) then
           call fms2_register_restart_field(fileObjWrite, CS%restart_field(m)%var_name, CS%var_ptr3d(m)%p, & 
               dimensions=axis_names(1:num_axes), domain_position=horgrid_position)

           ! prepare the restart field checksum
           !check_val(m-start_var+1,1) = mpp_chksum(CS%var_ptr3d(m)%p(isL:ieL,jsL:jeL,:))
        elseif (associated(CS%var_ptr2d(m)%p)) then

           call fms2_register_restart_field(fileObjWrite, CS%restart_field(m)%var_name, CS%var_ptr2d(m)%p, & 
               dimensions=axis_names(1:num_axes), domain_position=horgrid_position)

           ! prepare the restart field checksum
           !check_val(m-start_var+1,1) = mpp_chksum(CS%var_ptr2d(m)%p(isL:ieL,jsL:jeL))
        elseif (associated(CS%var_ptr4d(m)%p)) then

           call fms2_register_restart_field(fileObjWrite, CS%restart_field(m)%var_name, CS%var_ptr4d(m)%p, & 
               dimensions=axis_names(1:num_axes), domain_position=horgrid_position)

           ! prepare the restart field checksum
           !check_val(m-start_var+1,1) = mpp_chksum(CS%var_ptr4d(m)%p(isL:ieL,jsL:jeL,:,:))
        elseif (associated(CS%var_ptr1d(m)%p)) then
           ! need to explicitly define axis_names array for 1-D variable
           call fms2_register_restart_field(fileObjWrite, CS%restart_field(m)%var_name, CS%var_ptr1d(m)%p, & 
               dimensions=(/axis_names(1:num_axes)/), domain_position=horgrid_position)

           ! prepare the restart field checksum
           !check_val(m-start_var+1,1) = mpp_chksum(CS%var_ptr1d(m)%p)
        elseif (associated(CS%var_ptr0d(m)%p)) then
           ! need to explicitly define axis_names array for scalar variable
           call fms2_register_restart_field(fileObjWrite, CS%restart_field(m)%var_name, CS%var_ptr0d(m)%p, & 
               dimensions=(/axis_names(1:num_axes)/), domain_position=horgrid_position)

           ! prepare the restart field checksum
           !check_val(m-start_var+1,1) = mpp_chksum(CS%var_ptr0d(m)%p,pelist=(/mpp_pe()/))
        endif
        ! convert the checksum to a string
        !checksum_char = ''
        !checksum_char = convert_checksum_to_string(check_val(m,1))
        ! register the variable attributes

        call MOM_register_variable_attribute(fileObjWrite, CS%restart_field(m)%var_name, 'units', units)
        call MOM_register_variable_attribute(fileObjWrite, CS%restart_field(m)%var_name, 'long_name', longname) 
        !call MOM_register_variable_attribute(fileObjWrite, CS%restart_field(m)%var_name, 'checksum', trim(checksum_char))
        
        if(allocated(time_vals)) deallocate(time_vals)     
     enddo
     
     call fms2_write_restart(fileObjWrite)
     call fms2_close_file(fileObjWrite)

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
  integer :: i, m
  integer :: isL, ieL, jsL, jeL, is0, js0
  integer :: ntime, pos
  integer :: unit(CS%max_fields) ! The mpp unit of all open files..
  logical :: unit_is_global(CS%max_fields) ! True if the file is global.
  character(len=200) :: base_file_name
  character(len=1024) :: temp_file_name
  character(len=8)   :: hor_grid ! Variable grid info.
  real    :: t1, t2 ! Two times.
  real, allocatable :: time_vals(:)
  logical                          :: check_exist, is_there_a_checksum
  integer(LONG_KIND),dimension(3)  :: checksum_file
  integer(kind=8)                  :: checksum_data
  logical :: file_open_success = .false. ! returned by call to fms2_open_file 
  type(FmsNetcdfDomainFile_t) :: fileObjRead  ! fms2 data structure
  integer :: str_index
  character(len=96), dimension(:), allocatable :: variable_names(:) ! File variable names
  character(len=64) :: checksum_char
  character(len=4) :: axis_names(4)

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "restore_state: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)


  str_index = 0
  ! get the base restart file name
  temp_file_name=''
  if ((LEN_TRIM(filename) == 1 .and. filename(1:1) == 'F') & 
     .or. trim(filename)=='r') then
     temp_file_name = trim(CS%restartfile)
  else
     temp_file_name = trim(filename)
  endif

  ! append '.nc.' to the file name if it is missing
  base_file_name = ''

  str_index = INDEX(temp_file_name,'.nc') 
  if (str_index <=0) then
     base_file_name = trim(append_substring(temp_file_name, '.nc'))                         
  else 
     base_file_name = trim(temp_file_name)
  endif

  !Open the restart file.
  file_open_success = MOM_open_file(fileObjRead, trim(directory)//trim(base_file_name), "read", G, is_restart=.true.)

  call fms2_get_dimension_size(fileObjRead, "Time", ntime)
  if (ntime .lt. 1) then
    call MOM_error(FATAL, "MOM_restart: less than one time level in restart file.")
  endif
  allocate(time_vals(ntime))
  call fms2_read_data(fileObjRead, "Time", time_vals)
  t1 = time_vals(1)
  deallocate(time_vals)
  t2 = t1
  call mpp_max(t2)
  if (t1 .ne. t2) then
    call MOM_error(FATAL, "times are different in different restart files.")
  endif
  
  day = real_to_time(t1*86400.0)

 ! Register the horizontal axes that correspond to x and y of the domain.
  axis_names = (/'lath', 'lonh', 'latq', 'lonq'/)
  do i = 1,size(axis_names)
    if (fms2_dimension_exists(fileObjRead, trim(axis_names(i)))) then
      call MOM_register_axis(fileObjRead, trim(axis_names(i)))
    endif
  enddo

 !Read in each variable from the restart files.
  do m = 1, CS%novars
    varname = ''
    varname = trim(CS%restart_field(m)%var_name)

    !Check for obsolete fields
    do i = 1,CS%num_obsolete_vars
      if (adjustl(lowercase(trim(varname))) .eq. adjustl(lowercase(trim(CS%restart_obsolete(i)%field_name)))) then
        call MOM_error(FATAL, "MOM_restart restore_state: Attempting to use obsolete restart field "//&
                               trim(varname)//" - the new corresponding restart field is "//&
                               trim(CS%restart_obsolete(i)%replacement_name))
      endif
    enddo

    !Skip fields that have already been initialized
    if (CS%restart_field(m)%initialized) then
      cycle
    endif

    !Check if the variable is mandatory and present in the restart file(s)
    if (.not. fms2_variable_exists(fileObjRead, trim(varname))) then
      if (CS%restart_field(m)%mand_var) then
        call MOM_error(FATAL, "MOM_restart: Unable to find mandatory variable " &
                              //trim(varname)//" in restart files.")
      else
        CS%restart_field(m)%initialized = .false.
        cycle
      endif
    endif

    !Get the variable's "domain position."
    call query_vardesc(CS%restart_field(m)%vars, hor_grid=hor_grid, caller="restore_state")
    pos = get_horizontal_grid_position(hor_grid)

    !Register the restart fields and compute the checksums.
    if (associated(CS%var_ptr1d(m)%p)) then
      call fms2_register_restart_field(fileObjRead, trim(varname), CS%var_ptr1d(m)%p)
    elseif (associated(CS%var_ptr0d(m)%p)) then
      call fms2_register_restart_field(fileObjRead, trim(varname), CS%var_ptr0d(m)%p)
    elseif (associated(CS%var_ptr2d(m)%p)) then
      call fms2_register_restart_field(fileObjRead, trim(varname), CS%var_ptr2d(m)%p, domain_position=pos)
    elseif (associated(CS%var_ptr3d(m)%p)) then
      call fms2_register_restart_field(fileObjRead, trim(varname), CS%var_ptr3d(m)%p, domain_position=pos)
    elseif (associated(CS%var_ptr4d(m)%p)) then
      call fms2_register_restart_field(fileObjRead, trim(varname), CS%var_ptr4d(m)%p, domain_position=pos)
    else
      call MOM_error(FATAL, "MOM_restart restore_state: No pointers set for "//trim(varname))
    endif
    CS%restart_field(m)%initialized = .true.
  enddo

  !Read in restart data and then close the file.
  call fms2_read_restart(fileObjRead)

  call fms2_close_file(fileObjRead)
          
end subroutine restore_state

!> restart_files_exist determines whether any restart files exist in the restart directory
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
    num_files = get_num_restart_files('r', directory, G, CS)
  else
    num_files = get_num_restart_files(filename, directory, G, CS)
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
    CS%new_run = (get_num_restart_files('r', directory, G, CS) == 0)
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

!> get_num_restart_files determines the number of existing restart files.
!> @note This function replaces open_restart_units
function get_num_restart_files(filename, directory, G, CS) result(num_files)
  character(len=*),      intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files.
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  pointer     :: CS        !< The control structure returned by a previous
                                                  !! call to restart_init.
 
  integer :: num_files  !< The number of files (both automatically named restart
                        !! files and others explicitly in filename) that have been opened

  ! Local variables
  character(len=256) :: filepath  ! The path (dir/file) to the file being opened.
  character(len=256) :: fname     ! The name of the current file.
  character(len=8)   :: suffix    ! A suffix (like "_2") that is added to any
                                  ! additional restart files.
  integer :: num_restart     ! The number of restart files that have already
                             ! been opened.
  integer :: start_char      ! The location of the starting character in the
                             ! current file name.
  integer :: m, err
  logical :: fexists
  character(len=80) :: restartname
  character(len=1024) :: restartname_temp ! temporary location for restart name
  integer :: substring_index
  type(FmsNetcdfDomainFile_t) :: fileObjRead ! fms2 data structure

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "get_num_restart_files: Module must be initialized before it is used.")

! Check whether restart file(s) in 'filename' or, if filename is 'r', restart file(s) with the base name
! in CS%restart file, exist
  num_restart = 0
  start_char = 1

  do while (start_char <= len_trim(filename) )
     ! parse filename
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
    
     restartname = ''
     restartname_temp=''

     if ((fname(1:1)=='r') .and. ( len_trim(fname) == 1)) then
        restartname_temp = trim(CS%restartfile)
     else
        restartname_temp = trim(fname)
     endif
      
     err = 0
     if (num_restart > 0) err = 1 ! Avoid going through the file list twice.
     do while (err == 0)
        ! append '.nc' to the restart file name if it is missing
        substring_index = index('.nc', trim(restartname_temp))
        if (substring_index <= 0) then
           restartname = append_substring(restartname_temp,'.nc')
        else 
           restartname = restartname_temp
        endif
        filepath = trim(directory) // trim(restartname)

        fexists = MOM_open_file(fileObjRead, trim(filepath), "read", G, is_restart=.true.)
        if (fexists) then
           if (fms2_global_att_exists(fileObjRead,'NumFilesInSet')) then
              call fms2_get_global_attribute(fileObjRead, 'NumFilesInSet', num_restart)
           else
              num_restart = num_restart + 1
           endif
   
           if (is_root_pe()) then
              call MOM_error(NOTE, "MOM_restart: MOM run restarted using : "//trim(filepath))
           else
              err = 1 ; exit
           endif
        else
            call MOM_error(WARNING,"MOM_restart: Unable to find restart file(s) with base name : "//trim(filepath))
        endif

        call fms2_close_file(fileObjRead)
 
     enddo ! while (err == 0) loop
      
  enddo ! while (start_char < strlen(filename)) loop

  num_files = num_restart
  
end function get_num_restart_files

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
  allocate(CS%restart_obsolete(CS%max_fields))
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
  if (associated(CS%restart_obsolete)) deallocate(CS%restart_obsolete)
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
  integer(LONG_KIND),dimension(3) :: checksum_file !< checksum string corresponds to
                                                   !< values from up to 3 times 
  integer :: last
  integer :: start
  integer(LONG_KIND) :: checksumh
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

end module MOM_restart
