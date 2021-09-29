!> The MOM6 facility for reading and writing restart files, and querying what has been read.
module MOM_restart

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_checksums, only : chksum => rotated_field_chksum
use MOM_domains, only : PE_here, num_PEs
use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : create_file, file_type, fieldtype, file_exists, open_file, close_file
use MOM_io, only : MOM_read_data, read_data, MOM_write_field, read_field_chksum
use MOM_io, only : get_file_info, get_file_fields, get_field_atts, get_file_times
use MOM_io, only : vardesc, var_desc, query_vardesc, modify_vardesc, get_filename_appendix
use MOM_io, only : MULTIPLE, READONLY_FILE, SINGLE_FILE
use MOM_io, only : CENTER, CORNER, NORTH_FACE, EAST_FACE
use MOM_string_functions, only : lowercase
use MOM_time_manager,  only : time_type, time_type_to_real, real_to_time
use MOM_time_manager,  only : days_in_month, get_date, set_date
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

public restart_init, restart_end, restore_state, register_restart_field
public save_restart, query_initialized, restart_registry_lock, restart_init_end, vardesc
public restart_files_exist, determine_is_new_run, is_new_run
public register_restart_field_as_obsolete, register_restart_pair

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
  integer :: turns                  !< Number of quarter turns from input to model domain
  logical :: locked = .false.       !< If true this registry has been locked and no further restart
                                    !! fields can be added without explicitly unlocking the registry.

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
  !>@}
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

!> Register a pair of restart fields whose rotations map onto each other
interface register_restart_pair
  module procedure register_restart_pair_ptr2d
  module procedure register_restart_pair_ptr3d
  module procedure register_restart_pair_ptr4d
end interface register_restart_pair

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

!> Register a restart field as obsolete
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

  call lock_check(CS, var_desc)

  CS%novars = CS%novars+1
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

  call lock_check(CS, var_desc)

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

  call lock_check(CS, var_desc)

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

  call lock_check(CS, var_desc)

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

  call lock_check(CS, var_desc)

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


!> Register a pair of rotationally equivalent 2d restart fields
subroutine register_restart_pair_ptr2d(a_ptr, b_ptr, a_desc, b_desc, &
    mandatory, CS)
  real, dimension(:,:), target, intent(in) :: a_ptr   !< First field pointer
  real, dimension(:,:), target, intent(in) :: b_ptr   !< Second field pointer
  type(vardesc), intent(in) :: a_desc   !< First field descriptor
  type(vardesc), intent(in) :: b_desc   !< Second field descriptor
  logical, intent(in) :: mandatory      !< If true, abort if field is missing
  type(MOM_restart_CS), pointer :: CS   !< MOM restart control structure

  call lock_check(CS, a_desc)

  if (modulo(CS%turns, 2) /= 0) then
    call register_restart_field(b_ptr, a_desc, mandatory, CS)
    call register_restart_field(a_ptr, b_desc, mandatory, CS)
  else
    call register_restart_field(a_ptr, a_desc, mandatory, CS)
    call register_restart_field(b_ptr, b_desc, mandatory, CS)
  endif
end subroutine register_restart_pair_ptr2d


!> Register a pair of rotationally equivalent 3d restart fields
subroutine register_restart_pair_ptr3d(a_ptr, b_ptr, a_desc, b_desc, &
    mandatory, CS)
  real, dimension(:,:,:), target, intent(in) :: a_ptr   !< First field pointer
  real, dimension(:,:,:), target, intent(in) :: b_ptr   !< Second field pointer
  type(vardesc), intent(in) :: a_desc   !< First field descriptor
  type(vardesc), intent(in) :: b_desc   !< Second field descriptor
  logical, intent(in) :: mandatory      !< If true, abort if field is missing
  type(MOM_restart_CS), pointer :: CS   !< MOM restart control structure

  call lock_check(CS, a_desc)

  if (modulo(CS%turns, 2) /= 0) then
    call register_restart_field(b_ptr, a_desc, mandatory, CS)
    call register_restart_field(a_ptr, b_desc, mandatory, CS)
  else
    call register_restart_field(a_ptr, a_desc, mandatory, CS)
    call register_restart_field(b_ptr, b_desc, mandatory, CS)
  endif
end subroutine register_restart_pair_ptr3d


!> Register a pair of rotationally equivalent 2d restart fields
subroutine register_restart_pair_ptr4d(a_ptr, b_ptr, a_desc, b_desc, &
    mandatory, CS)
  real, dimension(:,:,:,:), target, intent(in) :: a_ptr !< First field pointer
  real, dimension(:,:,:,:), target, intent(in) :: b_ptr !< Second field pointer
  type(vardesc), intent(in) :: a_desc   !< First field descriptor
  type(vardesc), intent(in) :: b_desc   !< Second field descriptor
  logical, intent(in) :: mandatory      !< If true, abort if field is missing
  type(MOM_restart_CS), pointer :: CS   !< MOM restart control structure

  call lock_check(CS, a_desc)

  if (modulo(CS%turns, 2) /= 0) then
    call register_restart_field(b_ptr, a_desc, mandatory, CS)
    call register_restart_field(a_ptr, b_desc, mandatory, CS)
  else
    call register_restart_field(a_ptr, a_desc, mandatory, CS)
    call register_restart_field(b_ptr, b_desc, mandatory, CS)
  endif
end subroutine register_restart_pair_ptr4d


! The following provide alternate interfaces to register restarts.

!> Register a 4-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_4d(f_ptr, name, mandatory, CS, longname, units, &
                                     hor_grid, z_grid, t_grid)
  real, dimension(:,:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizontal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent

  type(vardesc) :: vd

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: " // &
      "register_restart_field_4d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  call lock_check(CS, name=name)

  vd = var_desc(name, units=units, longname=longname, hor_grid=hor_grid, &
                z_grid=z_grid, t_grid=t_grid)

  call register_restart_field_ptr4d(f_ptr, vd, mandatory, CS)

end subroutine register_restart_field_4d

!> Register a 3-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_3d(f_ptr, name, mandatory, CS, longname, units, &
                                     hor_grid, z_grid, t_grid)
  real, dimension(:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizontal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent

  type(vardesc) :: vd

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: " // &
      "register_restart_field_3d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  call lock_check(CS, name=name)

  vd = var_desc(name, units=units, longname=longname, hor_grid=hor_grid, &
                z_grid=z_grid, t_grid=t_grid)

  call register_restart_field_ptr3d(f_ptr, vd, mandatory, CS)

end subroutine register_restart_field_3d

!> Register a 2-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_2d(f_ptr, name, mandatory, CS, longname, units, &
                                     hor_grid, z_grid, t_grid)
  real, dimension(:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizontal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, '1' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent

  type(vardesc) :: vd
  character(len=8) :: Zgrid

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: " // &
      "register_restart_field_2d: Module must be initialized before "//&
      "it is used to register "//trim(name))
  zgrid = '1' ; if (present(z_grid)) zgrid = z_grid

  call lock_check(CS, name=name)

  vd = var_desc(name, units=units, longname=longname, hor_grid=hor_grid, &
                z_grid=zgrid, t_grid=t_grid)

  call register_restart_field_ptr2d(f_ptr, vd, mandatory, CS)

end subroutine register_restart_field_2d

!> Register a 1-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_1d(f_ptr, name, mandatory, CS, longname, units, &
                                     hor_grid, z_grid, t_grid)
  real, dimension(:), target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizontal staggering, '1' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent

  type(vardesc) :: vd
  character(len=8) :: hgrid

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: " // &
      "register_restart_field_3d: Module must be initialized before "//&
      "it is used to register "//trim(name))
  hgrid = '1' ; if (present(hor_grid)) hgrid = hor_grid

  call lock_check(CS, name=name)

  vd = var_desc(name, units=units, longname=longname, hor_grid=hgrid, &
                z_grid=z_grid, t_grid=t_grid)

  call register_restart_field_ptr1d(f_ptr, vd, mandatory, CS)

end subroutine register_restart_field_1d

!> Register a 0-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_0d(f_ptr, name, mandatory, CS, longname, units, &
                                     t_grid)
  real,               target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent

  type(vardesc) :: vd

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: " // &
      "register_restart_field_0d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  call lock_check(CS, name=name)

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

  integer :: m, n

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
  ! Assume that you are going to initialize it now, so set flag to initialized if queried again.
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

  integer :: m, n

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
  ! Assume that you are going to initialize it now, so set flag to initialized if queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_0d

!> Indicate whether the field pointed to by f_ptr has been initialized from a restart file.
function query_initialized_1d(f_ptr, CS) result(query_initialized)
  real, dimension(:), target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  type(MOM_restart_CS),       pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized

  integer :: m, n

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
  ! Assume that you are going to initialize it now, so set flag to initialized if queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_1d

!> Indicate whether the field pointed to by f_ptr has been initialized from a restart file.
function query_initialized_2d(f_ptr, CS) result(query_initialized)
  real, dimension(:,:), &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized

  integer :: m, n

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
  ! Assume that you are going to initialize it now, so set flag to initialized if queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_2d

!> Indicate whether the field pointed to by f_ptr has been initialized from a restart file.
function query_initialized_3d(f_ptr, CS) result(query_initialized)
  real, dimension(:,:,:), &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized

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
  ! Assume that you are going to initialize it now, so set flag to initialized if queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_3d

!> Indicate whether the field pointed to by f_ptr has been initialized from a restart file.
function query_initialized_4d(f_ptr, CS) result(query_initialized)
  real, dimension(:,:,:,:),  &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized

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
  ! Assume that you are going to initialize it now, so set flag to initialized if queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_4d

!> Indicate whether the field pointed to by f_ptr or with the specified variable
!! name has been initialized from a restart file.
function query_initialized_0d_name(f_ptr, name, CS) result(query_initialized)
  real,         target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized

  integer :: m, n

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
  ! Assume that you are going to initialize it now, so set flag to initialized if queried again.
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

  integer :: m, n

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
  ! Assume that you are going to initialize it now, so set flag to initialized if queried again.
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

  integer :: m, n

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
  ! Assume that you are going to initialize it now, so set flag to initialized if queried again.
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
  ! Assume that you are going to initialize it now, so set flag to initialized if queried again.
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
  ! Assume that you are going to initialize it now, so set flag to initialized if queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.
  if (n==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE, "MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "possibly because of the suspect comparison of pointers by ASSOCIATED.")
    query_initialized = query_initialized_name(name, CS)
  endif

end function query_initialized_4d_name

!> save_restart saves all registered variables to restart files.
subroutine save_restart(directory, time, G, CS, time_stamped, filename, GV, num_rest_files, write_IC)
  character(len=*),        intent(in)    :: directory !< The directory where the restart files
                                                  !! are to be written
  type(time_type),         intent(in)    :: time  !< The current model time
  type(ocean_grid_type),   intent(inout) :: G     !< The ocean's grid structure
  type(MOM_restart_CS),    pointer       :: CS    !< The control structure returned by a previous
                                                  !! call to restart_init
  logical,       optional, intent(in)    :: time_stamped !< If present and true, add time-stamp
                                                  !! to the restart file names
  character(len=*), optional, intent(in) :: filename !< A filename that overrides the name in CS%restartfile
  type(verticalGrid_type), &
                 optional, intent(in)    :: GV    !< The ocean's vertical grid structure
  integer,       optional, intent(out)   :: num_rest_files !< number of restart files written
  logical,       optional, intent(in)    :: write_IC !< If present and true, initial conditions
                                                  !! are being written

  ! Local variables
  type(vardesc) :: vars(CS%max_fields)  ! Descriptions of the fields that
                                        ! are to be read from the restart file.
  type(fieldtype) :: fields(CS%max_fields) ! Opaque types containing metadata describing
                                        ! each variable that will be written.
  character(len=512) :: restartpath     ! The restart file path (dir/file).
  character(len=256) :: restartname     ! The restart file name (no dir).
  character(len=8)   :: suffix          ! A suffix (like _2) that is appended
                                        ! to the name of files after the first.
  integer(kind=8) :: var_sz, size_in_file ! The size in bytes of each variable
                                        ! and the variables already in a file.
  integer(kind=8) :: max_file_size = 2147483647_8 ! The maximum size in bytes
                                        ! for any one file.  With NetCDF3,
                                        ! this should be 2 Gb or less.
  integer :: start_var, next_var        ! The starting variables of the
                                        ! current and next files.
  type(file_type) :: IO_handle          ! The I/O handle of the open fileset
  integer :: m, nz
  integer :: num_files                  ! The number of restart files that will be used.
  integer :: seconds, days, year, month, hour, minute
  character(len=8) :: hor_grid, z_grid, t_grid ! Variable grid info.
  character(len=64) :: var_name         ! A variable's name.
  real :: restart_time
  character(len=32) :: filename_appendix = '' ! Appendix to filename for ensemble runs
  integer :: length                     ! The length of a text string.
  integer(kind=8) :: check_val(CS%max_fields,1)
  integer :: isL, ieL, jsL, jeL, pos
  integer :: turns

  turns = CS%turns

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

  restartname = trim(CS%restartfile)
  if (present(filename)) restartname = trim(filename)
  if (PRESENT(time_stamped)) then ; if (time_stamped) then
    call get_date(time, year, month, days, hour, minute, seconds)
    ! Compute the year-day, because I don't like months. - RWH
    do m=1,month-1
      days = days + days_in_month(set_date(year, m, 2, 0, 0, 0))
    enddo
    seconds = seconds + 60*minute + 3600*hour
    if (year <= 9999) then
      write(restartname,'("_Y",I4.4,"_D",I3.3,"_S",I5.5)') year, days, seconds
    elseif (year <= 99999) then
      write(restartname,'("_Y",I5.5,"_D",I3.3,"_S",I5.5)') year, days, seconds
    else
      write(restartname,'("_Y",I10.10,"_D",I3.3,"_S",I5.5)') year, days, seconds
    endif
    restartname = trim(CS%restartfile)//trim(restartname)
  endif ; endif

  next_var = 1
  do while (next_var <= CS%novars )
    start_var = next_var
    size_in_file = 8*(2*G%Domain%niglobal+2*G%Domain%njglobal+2*nz+1000)

    do m=start_var,CS%novars
      call query_vardesc(CS%restart_field(m)%vars, hor_grid=hor_grid, &
                         z_grid=z_grid, t_grid=t_grid, caller="save_restart")
      var_sz = get_variable_byte_size(hor_grid, z_grid, t_grid, G, nz)

      if ((m==start_var) .OR. (size_in_file < max_file_size-var_sz)) then
        size_in_file = size_in_file + var_sz
      else ; exit
      endif

    enddo
    next_var = m

    ! Determine if there is a filename_appendix (used for ensemble runs).
    call get_filename_appendix(filename_appendix)
    if (len_trim(filename_appendix) > 0) then
      length = len_trim(restartname)
      if (restartname(length-2:length) == '.nc') then
        restartname = restartname(1:length-3)//'.'//trim(filename_appendix)//'.nc'
      else
        restartname = restartname(1:length)  //'.'//trim(filename_appendix)
      endif
    endif

    restartpath = trim(directory) // trim(restartname)

    if (num_files < 10) then
      write(suffix,'("_",I1)') num_files
    else
      write(suffix,'("_",I2)') num_files
    endif

    length = len_trim(restartpath)
    if (length < 3) then  ! This case is very uncommon but this test avoids segmentation-faults.
      if (num_files > 0) restartpath = trim(restartpath) // suffix
      restartpath = trim(restartpath)//".nc"
    elseif (restartpath(length-2:length) == ".nc") then
      if (num_files > 0) restartpath = restartpath(1:length-3)//trim(suffix)//".nc"
    else
      if (num_files > 0) restartpath = trim(restartpath) // suffix
      restartpath = trim(restartpath)//".nc"
    endif

    do m=start_var,next_var-1
      vars(m-start_var+1) = CS%restart_field(m)%vars
    enddo
    call query_vardesc(vars(1), t_grid=t_grid, hor_grid=hor_grid, caller="save_restart")
    t_grid = adjustl(t_grid)
    if (t_grid(1:1) /= 'p') &
      call modify_vardesc(vars(1), t_grid='s', caller="save_restart")
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

    !Prepare the checksum of the restart fields to be written to restart files
    if (modulo(turns, 2) /= 0) then
      call get_checksum_loop_ranges(G, pos, jsL, jeL, isL, ieL)
    else
      call get_checksum_loop_ranges(G, pos, isL, ieL, jsL, jeL)
    endif
    do m=start_var,next_var-1
      if (associated(CS%var_ptr3d(m)%p)) then
        check_val(m-start_var+1,1) = chksum(CS%var_ptr3d(m)%p(isL:ieL,jsL:jeL,:), turns=-turns)
      elseif (associated(CS%var_ptr2d(m)%p)) then
        check_val(m-start_var+1,1) = chksum(CS%var_ptr2d(m)%p(isL:ieL,jsL:jeL), turns=-turns)
      elseif (associated(CS%var_ptr4d(m)%p)) then
        check_val(m-start_var+1,1) = chksum(CS%var_ptr4d(m)%p(isL:ieL,jsL:jeL,:,:), turns=-turns)
      elseif (associated(CS%var_ptr1d(m)%p)) then
        check_val(m-start_var+1,1) = chksum(CS%var_ptr1d(m)%p)
      elseif (associated(CS%var_ptr0d(m)%p)) then
        check_val(m-start_var+1,1) = chksum(CS%var_ptr0d(m)%p, pelist=(/PE_here()/))
      endif
    enddo

    if (CS%parallel_restartfiles) then
      call create_file(IO_handle, trim(restartpath), vars, (next_var-start_var), &
                       fields, MULTIPLE, G=G, GV=GV, checksums=check_val)
    else
      call create_file(IO_handle, trim(restartpath), vars, (next_var-start_var), &
                       fields, SINGLE_FILE, G=G, GV=GV, checksums=check_val)
    endif

    do m=start_var,next_var-1
      if (associated(CS%var_ptr3d(m)%p)) then
        call MOM_write_field(IO_handle, fields(m-start_var+1), G%Domain, &
                         CS%var_ptr3d(m)%p, restart_time, turns=-turns)
      elseif (associated(CS%var_ptr2d(m)%p)) then
        call MOM_write_field(IO_handle, fields(m-start_var+1), G%Domain, &
                         CS%var_ptr2d(m)%p, restart_time, turns=-turns)
      elseif (associated(CS%var_ptr4d(m)%p)) then
        call MOM_write_field(IO_handle, fields(m-start_var+1), G%Domain, &
                         CS%var_ptr4d(m)%p, restart_time, turns=-turns)
      elseif (associated(CS%var_ptr1d(m)%p)) then
        call MOM_write_field(IO_handle, fields(m-start_var+1), CS%var_ptr1d(m)%p, restart_time)
      elseif (associated(CS%var_ptr0d(m)%p)) then
        call MOM_write_field(IO_handle, fields(m-start_var+1), CS%var_ptr0d(m)%p, restart_time)
      endif
    enddo

    call close_file(IO_handle)

    num_files = num_files+1

  enddo

  if (present(num_rest_files)) num_rest_files = num_files

end subroutine save_restart

!> restore_state reads the model state from previously generated files.  All
!! restart variables are read from the first file in the input filename list
!! in which they are found.
subroutine restore_state(filename, directory, day, G, CS)
  character(len=*),      intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(time_type),       intent(out) :: day       !< The time of the restarted run
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  pointer     :: CS        !< The control structure returned by a previous
                                                  !! call to restart_init

  ! Local variables
  character(len=200) :: filepath  ! The path (dir/file) to the file being opened.
  character(len=80) :: fname      ! The name of the current file.
  character(len=8)  :: suffix     ! A suffix (like "_2") that is added to any
                                  ! additional restart files.
  character(len=512) :: mesg      ! A message for warnings.
  character(len=80) :: varname    ! A variable's name.
  integer :: num_file        ! The number of files (restart files and others
                             ! explicitly in filename) that are open.
  integer :: i, n, m, missing_fields
  integer :: isL, ieL, jsL, jeL, is0, js0
  integer :: sizes(7)
  integer :: nvar, ntime, pos

  type(file_type) :: IO_handles(CS%max_fields) ! The I/O units of all open files.
  character(len=200) :: unit_path(CS%max_fields) ! The file names.
  logical :: unit_is_global(CS%max_fields) ! True if the file is global.

  character(len=8)   :: hor_grid ! Variable grid info.
  real    :: t1, t2 ! Two times.
  real, allocatable :: time_vals(:)
  type(fieldtype), allocatable :: fields(:)
  logical            :: is_there_a_checksum ! Is there a valid checksum that should be checked.
  integer(kind=8)    :: checksum_file  ! The checksum value recorded in the input file.
  integer(kind=8)    :: checksum_data  ! The checksum value for the data that was read in.

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "restore_state: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  ! Get NetCDF ids for all of the restart files.
  if ((LEN_TRIM(filename) == 1) .and. (filename(1:1) == 'F')) then
    num_file = open_restart_units('r', directory, G, CS, IO_handles=IO_handles, &
                     file_paths=unit_path, global_files=unit_is_global)
  else
    num_file = open_restart_units(filename, directory, G, CS, IO_handles=IO_handles, &
                     file_paths=unit_path, global_files=unit_is_global)
  endif

  if (num_file == 0) then
    write(mesg,'("Unable to find any restart files specified by  ",A,"  in directory ",A,".")') &
                  trim(filename), trim(directory)
    call MOM_error(FATAL,"MOM_restart: "//mesg)
  endif

  ! Get the time from the first file in the list that has one.
  do n=1,num_file
    call get_file_times(IO_handles(n), time_vals, ntime)
    if (ntime < 1) cycle

    t1 = time_vals(1)
    deallocate(time_vals)

    day = real_to_time(t1*86400.0)
    exit
  enddo

  if (n>num_file) call MOM_error(WARNING,"MOM_restart: " // &
                                 "No times found in restart files.")

  ! Check the remaining files for different times and issue a warning
  ! if they differ from the first time.
    do m = n+1,num_file
      call get_file_times(IO_handles(n), time_vals, ntime)
      if (ntime < 1) cycle

      t2 = time_vals(1)
      deallocate(time_vals)

      if (t1 /= t2 .and. is_root_PE()) then
        write(mesg,'("WARNING: Restart file ",I2," has time ",F10.4,"whereas &
         &simulation is restarted at ",F10.4," (differing by ",F10.4,").")')&
               m,t1,t2,t1-t2
        call MOM_error(WARNING, "MOM_restart: "//mesg)
      endif
    enddo

  ! Read each variable from the first file in which it is found.
  do n=1,num_file
    call get_file_info(IO_handles(n), nvar=nvar)

    allocate(fields(nvar))
    call get_file_fields(IO_handles(n), fields(1:nvar))

    do m=1, nvar
      call get_field_atts(fields(m), name=varname)
      do i=1,CS%num_obsolete_vars
        if (adjustl(lowercase(trim(varname))) == adjustl(lowercase(trim(CS%restart_obsolete(i)%field_name)))) then
            call MOM_error(FATAL, "MOM_restart restore_state: Attempting to use obsolete restart field "//&
                           trim(varname)//" - the new corresponding restart field is "//&
                           trim(CS%restart_obsolete(i)%replacement_name))
        endif
      enddo
    enddo

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
        call get_field_atts(fields(i), name=varname)
        if (lowercase(trim(varname)) == lowercase(trim(CS%restart_field(m)%var_name))) then
          checksum_data = -1
          if (CS%checksum_required) then
            call read_field_chksum(fields(i), checksum_file, is_there_a_checksum)
          else
            checksum_file = -1
            is_there_a_checksum = .false. ! Do not need to do data checksumming.
          endif

          if (associated(CS%var_ptr1d(m)%p))  then
            ! Read a 1d array, which should be invariant to domain decomposition.
            call MOM_read_data(unit_path(n), varname, CS%var_ptr1d(m)%p, &
                               timelevel=1, MOM_Domain=G%Domain)
            if (is_there_a_checksum) checksum_data = chksum(CS%var_ptr1d(m)%p)
          elseif (associated(CS%var_ptr0d(m)%p)) then ! Read a scalar...
            call MOM_read_data(unit_path(n), varname, CS%var_ptr0d(m)%p, &
                               timelevel=1, MOM_Domain=G%Domain)
            if (is_there_a_checksum) checksum_data = chksum(CS%var_ptr0d(m)%p, pelist=(/PE_here()/))
          elseif (associated(CS%var_ptr2d(m)%p)) then  ! Read a 2d array.
            if (pos /= 0) then
              call MOM_read_data(unit_path(n), varname, CS%var_ptr2d(m)%p, &
                                 G%Domain, timelevel=1, position=pos)
            else ! This array is not domain-decomposed.  This variant may be under-tested.
              call MOM_error(FATAL, &
                        "MOM_restart does not support 2-d arrays without domain decomposition.")
              ! call read_data(unit_path(n), varname, CS%var_ptr2d(m)%p,no_domain=.true., timelevel=1)
            endif
            if (is_there_a_checksum) checksum_data = chksum(CS%var_ptr2d(m)%p(isL:ieL,jsL:jeL))
          elseif (associated(CS%var_ptr3d(m)%p)) then  ! Read a 3d array.
            if (pos /= 0) then
              call MOM_read_data(unit_path(n), varname, CS%var_ptr3d(m)%p, &
                                 G%Domain, timelevel=1, position=pos)
            else ! This array is not domain-decomposed.  This variant may be under-tested.
              call MOM_error(FATAL, &
                        "MOM_restart does not support 3-d arrays without domain decomposition.")
              ! call read_data(unit_path(n), varname, CS%var_ptr3d(m)%p, no_domain=.true., timelevel=1)
            endif
            if (is_there_a_checksum) checksum_data = chksum(CS%var_ptr3d(m)%p(isL:ieL,jsL:jeL,:))
          elseif (associated(CS%var_ptr4d(m)%p)) then  ! Read a 4d array.
            if (pos /= 0) then
              call MOM_read_data(unit_path(n), varname, CS%var_ptr4d(m)%p, &
                                 G%Domain, timelevel=1, position=pos)
            else ! This array is not domain-decomposed.  This variant may be under-tested.
              call MOM_error(FATAL, &
                        "MOM_restart does not support 4-d arrays without domain decomposition.")
              ! call read_data(unit_path(n), varname, CS%var_ptr4d(m)%p, no_domain=.true., timelevel=1)
            endif
            if (is_there_a_checksum) checksum_data = chksum(CS%var_ptr4d(m)%p(isL:ieL,jsL:jeL,:,:))
          else
            call MOM_error(FATAL, "MOM_restart restore_state: No pointers set for "//trim(varname))
          endif

          if (is_root_pe() .and. is_there_a_checksum .and. (checksum_file /= checksum_data)) then
             write (mesg,'(a,Z16,a,Z16,a)') "Checksum of input field "// trim(varname)//" ",checksum_data,&
                                          " does not match value ", checksum_file, &
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
    if (missing_fields == 0) exit
  enddo

  do n=1,num_file
    call close_file(IO_handles(n))
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

  ! Lock the restart registry so that no further variables can be registered.
  CS%locked = .true.

end subroutine restore_state

!> restart_files_exist determines whether any restart files exist.
function restart_files_exist(filename, directory, G, CS)
  character(len=*),      intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  pointer     :: CS        !< The control structure returned by a previous
                                                  !! call to restart_init
  logical :: restart_files_exist                  !< The function result, which indicates whether
                                                  !! any of the explicitly or automatically named
                                                  !! restart files exist in directory
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
                                                  !! character 'r' to read automatically named files
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  pointer     :: CS        !< The control structure returned by a previous
                                                  !! call to restart_init
  logical :: is_new_run                           !< The function result, which indicates whether
                                                  !! this is a new run, based on the value of
                                                  !! filename and whether restart files exist

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
                                       !! call to restart_init
  logical :: is_new_run                !< The function result, which had been stored in CS during
                                       !! a previous call to determine_is_new_run

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "is_new_run: Module must be initialized before it is used.")
  if (.not.CS%new_run_set) call MOM_error(FATAL, "MOM_restart " // &
      "determine_is_new_run must be called for a restart file before is_new_run.")

  is_new_run = CS%new_run
end function is_new_run

!> open_restart_units determines the number of existing restart files and optionally opens
!! them and returns unit ids, paths and whether the files are global or spatially decomposed.
function open_restart_units(filename, directory, G, CS, IO_handles, file_paths, &
                            global_files) result(num_files)
  character(len=*),      intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  pointer     :: CS        !< The control structure returned by a previous
                                                  !! call to restart_init
  type(file_type), dimension(:), &
               optional, intent(out) :: IO_handles !< The I/O handles of all opened files
  character(len=*), dimension(:), &
               optional, intent(out) :: file_paths   !< The full paths to open files
  logical, dimension(:), &
               optional, intent(out) :: global_files !< True if a file is global

  integer :: num_files  !< The number of files (both automatically named restart
                        !! files and others explicitly in filename) that have been opened.

  ! Local variables
  character(len=256) :: filepath  ! The path (dir/file) to the file being opened.
  character(len=256) :: fname     ! The name of the current file.
  character(len=8)   :: suffix    ! A suffix (like "_2") that is added to any
                                  ! additional restart files.
  integer :: num_restart     ! The number of restart files that have already
                             ! been opened using their numbered suffix.
  integer :: start_char      ! The location of the starting character in the
                             ! current file name.
  integer :: nf              ! The number of files that have been found so far
  integer :: m, length
  logical :: still_looking   ! If true, the code is still looking for automatically named files
  logical :: fexists         ! True if a file has been found
  character(len=32) :: filename_appendix = '' ! Filename appendix for ensemble runs
  character(len=80) :: restartname

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "open_restart_units: Module must be initialized before it is used.")

  ! Get NetCDF ids for all of the restart files.
  num_restart = 0 ; nf = 0 ; start_char = 1
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
      still_looking = (num_restart <= 0) ! Avoid going through the file list twice.
      do while (still_looking)
        restartname = trim(CS%restartfile)

        ! Determine if there is a filename_appendix (used for ensemble runs).
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

        filepath = trim(filepath)//".nc"

        num_restart = num_restart + 1
        ! Look for a global netCDF file.
        inquire(file=filepath, exist=fexists)
        if (fexists) then
          nf = nf + 1
          if (present(IO_handles)) &
            call open_file(IO_handles(nf), trim(filepath), READONLY_FILE, &
                           threading=MULTIPLE, fileset=SINGLE_FILE)
          if (present(global_files)) global_files(nf) = .true.
          if (present(file_paths)) file_paths(nf) = filepath
        elseif (CS%parallel_restartfiles) then
          ! Look for decomposed files using the I/O Layout.
          fexists = file_exists(filepath, G%Domain)
          if (fexists) then
            nf = nf + 1
            if (present(IO_handles)) &
              call open_file(IO_handles(nf), trim(filepath), READONLY_FILE, MOM_domain=G%Domain)
            if (present(global_files)) global_files(nf) = .false.
            if (present(file_paths)) file_paths(nf) = filepath
          endif
        endif

        if (fexists) then
          if (is_root_pe() .and. (present(IO_handles))) &
            call MOM_error(NOTE, "MOM_restart: MOM run restarted using : "//trim(filepath))
        else
          still_looking = .false. ; exit
        endif
      enddo ! while (still_looking) loop
    else
      filepath = trim(directory)//trim(fname)
      inquire(file=filepath, exist=fexists)
      if (.not. fexists) filepath = trim(filepath)//".nc"

      inquire(file=filepath, exist=fexists)
      if (fexists) then
        nf = nf + 1
        if (present(IO_handles)) &
          call open_file(IO_handles(nf), trim(filepath), READONLY_FILE, &
                       threading=MULTIPLE, fileset=SINGLE_FILE)
        if (present(global_files)) global_files(nf) = .true.
        if (present(file_paths)) file_paths(nf) = filepath
        if (is_root_pe() .and. (present(IO_handles))) &
          call MOM_error(NOTE,"MOM_restart: MOM run restarted using : "//trim(filepath))
      else
        if (present(IO_handles)) &
          call MOM_error(WARNING,"MOM_restart: Unable to find restart file : "//trim(filepath))
      endif

    endif
  enddo ! while (start_char < len_trim(filename)) loop
  num_files = nf

end function open_restart_units

!> get_num_restart_files returns the number of existing restart files that match the provided
!! directory structure and other information stored in the control structure and optionally
!! also provides the full paths to these files.
function get_num_restart_files(filenames, directory, G, CS, file_paths) result(num_files)
  character(len=*),      intent(in)  :: filenames !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  pointer     :: CS        !< The control structure returned by a previous
                                                  !! call to restart_init
  character(len=*), dimension(:), &
               optional, intent(out) :: file_paths !< The full paths to the restart files.
  integer :: num_files  !< The function result, the number of files (both automatically named
                        !! restart files and others explicitly in filename) that have been opened

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "get_num_restart_files: Module must be initialized before it is used.")

  ! This call uses open_restart_units without the optional arguments needed to actually
  ! open the files to determine the number of restart files.
  num_files = open_restart_units(filenames, directory, G, CS, file_paths=file_paths)

end function get_num_restart_files


!> Initialize this module and set up a restart control structure.
subroutine restart_init(param_file, CS, restart_root)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
  type(MOM_restart_CS),  pointer    :: CS !< A pointer to a MOM_restart_CS object that is allocated here
  character(len=*), optional, &
                         intent(in) :: restart_root !< A filename root that overrides the value
                                          !! set by RESTARTFILE to enable the use of this module by
                                          !! other components than MOM.

  logical :: rotate_index

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_restart"   ! This module's name.
  logical :: all_default   ! If true, all parameters are using their default values.

  if (associated(CS)) then
    call MOM_error(WARNING, "restart_init called with an associated control structure.")
    return
  endif
  allocate(CS)

  ! Determine whether all paramters are set to their default values.
  call get_param(param_file, mdl, "PARALLEL_RESTARTFILES", CS%parallel_restartfiles, &
                 default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "LARGE_FILE_SUPPORT", CS%large_file_support, &
                 default=.true., do_not_log=.true.)
  call get_param(param_file, mdl, "MAX_FIELDS", CS%max_fields, default=100, do_not_log=.true.)
  call get_param(param_file, mdl, "RESTART_CHECKSUMS_REQUIRED", CS%checksum_required, &
                 default=.true., do_not_log=.true.)
  all_default = ((.not.CS%parallel_restartfiles) .and. (CS%large_file_support) .and. &
                 (CS%max_fields == 100) .and. (CS%checksum_required))
  if (.not.present(restart_root)) then
    call get_param(param_file, mdl, "RESTARTFILE", CS%restartfile, &
                   default="MOM.res", do_not_log=.true.)
    all_default = (all_default .and. (trim(CS%restartfile) == trim("MOM.res")))
  endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "", all_default=all_default)
  call get_param(param_file, mdl, "PARALLEL_RESTARTFILES", CS%parallel_restartfiles, &
                 "If true, each processor writes its own restart file, "//&
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
                 "If true, use the file-size limits with NetCDF large "//&
                 "file support (4Gb), otherwise the limit is 2Gb.", &
                 default=.true.)
  call get_param(param_file, mdl, "MAX_FIELDS", CS%max_fields, &
                 "The maximum number of restart fields that can be used.", &
                 default=100)
  call get_param(param_file, mdl, "RESTART_CHECKSUMS_REQUIRED", CS%checksum_required, &
                 "If true, require the restart checksums to match and error out otherwise. "//&
                 "Users may want to avoid this comparison if for example the restarts are "//&
                 "made from a run with a different mask_table than the current run, "//&
                 "in which case the checksums will not match and cause crash.",&
                 default=.true.)

  ! Maybe not the best place to do this?
  call get_param(param_file, mdl, "ROTATE_INDEX", rotate_index, &
      default=.false., do_not_log=.true.)

  CS%turns = 0
  if (rotate_index) then
    call get_param(param_file, mdl, "INDEX_TURNS", CS%turns, &
        default=1, do_not_log=.true.)
  endif

  allocate(CS%restart_field(CS%max_fields))
  allocate(CS%restart_obsolete(CS%max_fields))
  allocate(CS%var_ptr0d(CS%max_fields))
  allocate(CS%var_ptr1d(CS%max_fields))
  allocate(CS%var_ptr2d(CS%max_fields))
  allocate(CS%var_ptr3d(CS%max_fields))
  allocate(CS%var_ptr4d(CS%max_fields))

  CS%locked = .false.

end subroutine restart_init

!> Issue an error message if the restart_registry is locked.
subroutine lock_check(CS, var_desc, name)
  type(MOM_restart_CS),       intent(in) :: CS        !< A MOM_restart_CS object (intent in)
  type(vardesc),    optional, intent(in) :: var_desc  !< A structure with metadata about this variable
  character(len=*), optional, intent(in) :: name      !< variable name to be used in the restart file

  character(len=256) :: var_name  ! A variable name.

  if (CS%locked) then
    if (present(var_desc)) then
      call query_vardesc(var_desc, name=var_name)
      call MOM_error(FATAL, "Attempted to register "//trim(var_name)//" but the restart registry is locked.")
    elseif (present(name)) then
      call MOM_error(FATAL, "Attempted to register "//trim(name)//" but the restart registry is locked.")
    else
      call MOM_error(FATAL, "Attempted to register a variable but the restart registry is locked.")
    endif
  endif

end subroutine lock_check

!> Lock the restart registry so that an error is issued if any further restart variables are registered.
subroutine restart_registry_lock(CS, unlocked)
  type(MOM_restart_CS), intent(inout) :: CS        !< A MOM_restart_CS object (intent inout)
  logical, optional,    intent(in)    :: unlocked  !< If present and true, unlock the registry

  CS%locked = .true.
  if (present(unlocked)) CS%locked = .not.unlocked
end subroutine restart_registry_lock

!> Indicate that all variables have now been registered and lock the registry.
subroutine restart_init_end(CS)
  type(MOM_restart_CS),  pointer    :: CS !< A pointer to a MOM_restart_CS object

  if (associated(CS)) then
    CS%locked = .true.

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

!> get the size of a variable in bytes
function get_variable_byte_size(hor_grid, z_grid, t_grid, G, num_z) result(var_sz)
  character(len=8),      intent(in) :: hor_grid !< The horizontal grid string to interpret
  character(len=8),      intent(in) :: z_grid   !< The vertical grid string to interpret
  character(len=8),      intent(in) :: t_grid   !< A time string to interpret
  type(ocean_grid_type), intent(in) :: G        !< The ocean's grid structure
  integer,               intent(in) :: num_z    !< The number of vertical layers in the grid
  integer(kind=8) :: var_sz !< The function result, the size in bytes of a variable

  ! Local variables
  integer :: var_periods  ! The number of entries in a time-periodic axis
  character(len=8) :: t_grid_read, t_grid_tmp ! Modified versions of t_grid

  if (trim(hor_grid) == '1') then
    var_sz = 8
  else ! This may be an overestimate, as it is based on symmetric-memory corner points.
    var_sz = 8*(G%Domain%niglobal+1)*(G%Domain%njglobal+1)
  endif

  select case (trim(z_grid))
    case ('L') ; var_sz = var_sz * num_z
    case ('i') ; var_sz = var_sz * (num_z+1)
  end select

  t_grid_tmp = adjustl(t_grid)
  if (t_grid_tmp(1:1) == 'p') then
    if (len_trim(t_grid_tmp(2:8)) > 0) then
      var_periods = -1
      t_grid_read = adjustl(t_grid_tmp(2:8))
      read(t_grid_read,*) var_periods
      if (var_periods > 1) var_sz = var_sz * var_periods
    endif
  endif

end function get_variable_byte_size

end module MOM_restart
