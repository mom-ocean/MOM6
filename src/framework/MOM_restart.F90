module MOM_restart
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, April 1994 - June 2002                         *
!*                                                                     *
!*    This file contains four subroutines associated with saving       *
!*  restart files or restoring the model state from files.             *
!*                                                                     *
!*    register_restart_field is used to specify the fields that will   *
!*  be written to restart files.                                       *
!*                                                                     *
!*    Save_restart saves a restart file from which a simulation can    *
!*  be restarted with results that are identical to those which would  *
!*  have been attained if there had been no interruption.  If this     *
!*  file would be larger than 2 Gbytes, it is broken up into a number  *
!*  of smaller files.                                                  *
!*                                                                     *
!*    The subroutine restore_state initializes the fields for the      *
!*  simulations from a number of restart files or other NetCDF files.  *
!*  Each restart field is initialized from the first file in the       *
!*  list in which it is found.  The files are separated by spaces,     *
!*  and all must be in the specified directory.  If 'r' is included    *
!*  in the list, it is expanded to include all of the restart files    *
!*  that are found in the directory.                                   *
!*                                                                     *
!*    query_initialized returns true if a field (or the entire restart *
!*  file) has been initialized from a restart file and false otherwise.*
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, CoriolisBu                            *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, bathyT, tr                            *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1                                                  *
!*           i  i+1                                                    *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_domains, only : pe_here, num_PEs
use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_string_functions, only : lowercase
use MOM_grid, only : ocean_grid_type
use MOM_io, only : create_file, fieldtype, file_exists, open_file, close_file
use MOM_io, only : read_field, write_field, read_data, get_filename_appendix
use MOM_io, only : get_file_info, get_file_atts, get_file_fields, get_file_times
use MOM_io, only : vardesc, query_vardesc, modify_vardesc
use MOM_io, only : MULTIPLE, NETCDF_FILE, READONLY_FILE, SINGLE_FILE
use MOM_io, only : CENTER, CORNER, NORTH_FACE, EAST_FACE
use MOM_time_manager, only : time_type, get_time, get_date, set_date, set_time
use MOM_time_manager, only : days_in_month
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

public restart_init, restart_end, restore_state, register_restart_field
public save_restart, query_initialized, restart_init_end, vardesc

type p4d
  real, dimension(:,:,:,:), pointer :: p => NULL()
end type p4d

type p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d

type p2d
  real, dimension(:,:), pointer :: p => NULL()
end type p2d

type p1d
  real, dimension(:), pointer :: p => NULL()
end type p1d

type p0d
  real, pointer :: p => NULL()
end type p0d

type field_restart
  type(vardesc) :: vars         ! Descriptions of the fields that
                                ! are to be read from or  written
                                ! to the restart file.
  logical :: mand_var           ! If .true. the run will abort if this
                                ! field is not successfully read
                                ! from the restart file.
  logical :: initialized        ! .true. if this field has been read
                                ! from the restart file.
  character(len=32) :: var_name ! A name by which a variable may be queried.
end type field_restart

type, public :: MOM_restart_CS ; private
  logical :: restart    ! restart is set to .true. if the run has been started
                        ! from a full restart file.  Otherwise some fields must
                        ! be initialized approximately.
  integer :: novars = 0 ! The number of restart fields that have been registered.
  logical :: parallel_restartfiles  ! If true, each PE writes its own restart file,
                                    ! otherwise they are combined internally.
  logical :: large_file_support     ! If true, NetCDF 3.6 or later is being used
                                    ! and large-file-support is enabled.
  character(len=240) :: restartfile ! The name or name root for MOM restart files.

  type(field_restart), pointer :: restart_field(:) => NULL()
  type(p0d), pointer :: var_ptr0d(:) => NULL()
  type(p1d), pointer :: var_ptr1d(:) => NULL()
  type(p2d), pointer :: var_ptr2d(:) => NULL()
  type(p3d), pointer :: var_ptr3d(:) => NULL()
  type(p4d), pointer :: var_ptr4d(:) => NULL()
  integer :: max_fields
end type MOM_restart_CS

interface register_restart_field
  module procedure register_restart_field_ptr4d
  module procedure register_restart_field_ptr3d
  module procedure register_restart_field_ptr2d
  module procedure register_restart_field_ptr1d
  module procedure register_restart_field_ptr0d
end interface

interface query_initialized
  module procedure query_initialized_name
  module procedure query_initialized_0d, query_initialized_0d_name
  module procedure query_initialized_1d, query_initialized_1d_name
  module procedure query_initialized_2d, query_initialized_2d_name
  module procedure query_initialized_3d, query_initialized_3d_name
  module procedure query_initialized_4d, query_initialized_4d_name
end interface

contains

subroutine register_restart_field_ptr3d(f_ptr, var_desc, mandatory, CS)
  real, dimension(:,:,:), target :: f_ptr
  type(vardesc),      intent(in) :: var_desc
  logical,            intent(in) :: mandatory
  type(MOM_restart_CS),  pointer :: CS
!  Set up a field that will be written to and read from restart
!  files.
!
! Arguments: f_ptr - A pointer to the field to be read or written.
!  (in)      var_desc - The descriptive structure for the field.
!  (in)      mandatory - If .true. the run will abort if this field is not
!                        successfully read from the restart file.  If .false.,
!                        alternate techniques are provided to initialize this
!                        field if it is cannot be read from the file.
!  (in/out)  CS - The control structure returned by a previous call to
!                 restart_init.
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
                     caller="register_restart_field_ptr3d")

  CS%var_ptr3d(CS%novars)%p => f_ptr
  CS%var_ptr4d(CS%novars)%p => NULL()
  CS%var_ptr2d(CS%novars)%p => NULL()
  CS%var_ptr1d(CS%novars)%p => NULL()
  CS%var_ptr0d(CS%novars)%p => NULL()

end subroutine register_restart_field_ptr3d

subroutine register_restart_field_ptr4d(f_ptr, var_desc, mandatory, CS)
  real, dimension(:,:,:,:), target :: f_ptr
  type(vardesc),      intent(in) :: var_desc
  logical,            intent(in) :: mandatory
  type(MOM_restart_CS),  pointer :: CS
!  Set up a field that will be written to and read from restart
!  files.
!
! Arguments: f_ptr - A pointer to the field to be read or written.
!  (in)      var_desc - The descriptive structure for the field.
!  (in)      mandatory - If .true. the run will abort if this field is not
!                        successfully read from the restart file.  If .false.,
!                        alternate techniques are provided to initialize this
!                        field if it is cannot be read from the file.
!  (in/out)  CS - The control structure returned by a previous call to
!                 restart_init.
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

subroutine register_restart_field_ptr2d(f_ptr, var_desc, mandatory, CS)
  real, dimension(:,:), target :: f_ptr
  type(vardesc), intent(in) :: var_desc
  logical, intent(in)       :: mandatory
  type(MOM_restart_CS),  pointer :: CS
!  Set up a field that will be written to and read from restart
!  files.
!
! Arguments: f_ptr - A pointer to the field to be read or written.
!  (in)      var_desc - The descriptive structure for the field.
!  (in)      mandatory - If .true. the run will abort if this field is not
!                        successfully read from the restart file.  If .false.,
!                        alternate techniques are provided to initialize this
!                        field if it is cannot be read from the file.
!  (in/out)  CS - The control structure returned by a previous call to
!                 restart_init.
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

subroutine register_restart_field_ptr1d(f_ptr, var_desc, mandatory, CS)
  real, dimension(:), target :: f_ptr
  type(vardesc), intent(in) :: var_desc
  logical, intent(in)       :: mandatory
  type(MOM_restart_CS),  pointer :: CS
!  Set up a field that will be written to and read from restart
!  files.
!
! Arguments: f_ptr - A pointer to the field to be read or written.
!  (in)      var_desc - The descriptive structure for the field.
!  (in)      mandatory - If .true. the run will abort if this field is not
!                        successfully read from the restart file.  If .false.,
!                        alternate techniques are provided to initialize this
!                        field if it is cannot be read from the file.
!  (in/out)  CS - The control structure returned by a previous call to
!                 restart_init.
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

subroutine register_restart_field_ptr0d(f_ptr, var_desc, mandatory, CS)
  real, target :: f_ptr
  type(vardesc), intent(in) :: var_desc
  logical, intent(in)       :: mandatory
  type(MOM_restart_CS),  pointer :: CS
!  Set up a field that will be written to and read from restart
!  files.
!
! Arguments: f_ptr - A pointer to the field to be read or written.
!  (in)      var_desc - The descriptive structure for the field.
!  (in)      mandatory - If .true. the run will abort if this field is not
!                        successfully read from the restart file.  If .false.,
!                        alternate techniques are provided to initialize this
!                        field if it is cannot be read from the file.
!  (in/out)  CS - The control structure returned by a previous call to
!                 restart_init.
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

function query_initialized_name(name, CS) result(query_initialized)
  character(len=*) :: name
  type(MOM_restart_CS), pointer :: CS
  logical :: query_initialized
!   This subroutine returns .true. if the field referred to by name has
! initialized from a restart file, and .false. otherwise.
!
! Arguments: name - A pointer to the field that is being queried.
!  (in)      CS - The control structure returned by a previous call to
!                 restart_init.
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

function query_initialized_0d(f_ptr, CS) result(query_initialized)
  real,                 target  :: f_ptr
  type(MOM_restart_CS), pointer :: CS
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr has
! been initialized from a restart file.
!
! Arguments: f_ptr - A pointer to the field that is being queried.
!  (in)      CS - The control structure returned by a previous call to
!                 restart_init.
  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (ASSOCIATED(CS%var_ptr0d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_0d

function query_initialized_1d(f_ptr, CS) result(query_initialized)
  real, dimension(:),   target  :: f_ptr
  type(MOM_restart_CS), pointer :: CS
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr has
! been initialized from a restart file.
!
! Arguments: f_ptr - A pointer to the field that is being queried.
!  (in)      CS - The control structure returned by a previous call to
!                 restart_init.
  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (ASSOCIATED(CS%var_ptr1d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_1d

function query_initialized_2d(f_ptr, CS) result(query_initialized)
  real, dimension(:,:), target  :: f_ptr
  type(MOM_restart_CS), pointer :: CS
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr has
! been initialized from a restart file.
!
! Arguments: f_ptr - A pointer to the field that is being queried.
!  (in)      CS - The control structure returned by a previous call to
!                 restart_init.
  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (ASSOCIATED(CS%var_ptr2d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_2d

function query_initialized_3d(f_ptr, CS) result(query_initialized)
  real, dimension(:,:,:), target  :: f_ptr
  type(MOM_restart_CS),   pointer :: CS
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr has
! been initialized from a restart file.
!
! Arguments: f_ptr - A pointer to the field that is being queried.
!  (in)      CS - The control structure returned by a previous call to
!                 restart_init.
  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (ASSOCIATED(CS%var_ptr3d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_3d

function query_initialized_4d(f_ptr, CS) result(query_initialized)
  real, dimension(:,:,:,:), target  :: f_ptr
  type(MOM_restart_CS),     pointer :: CS
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr has
! been initialized from a restart file.
!
! Arguments: f_ptr - A pointer to the field that is being queried.
!  (in)      CS - The control structure returned by a previous call to
!                 restart_init.
  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (ASSOCIATED(CS%var_ptr4d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_4d

function query_initialized_0d_name(f_ptr, name, CS) result(query_initialized)
  real,                 target  :: f_ptr
  character(len=*)              :: name
  type(MOM_restart_CS), pointer :: CS
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr or with the
! specified variable name has been initialized from a restart file.
!
! Arguments: f_ptr - A pointer to the field that is being queried.
!  (in)      name - The name of the field that is being queried.
!  (in)      CS - The control structure returned by a previous call to
!                 restart_init.
  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (ASSOCIATED(CS%var_ptr0d(m)%p,f_ptr)) then
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

function query_initialized_1d_name(f_ptr, name, CS) result(query_initialized)
  real, dimension(:),   target  :: f_ptr
  character(len=*)              :: name
  type(MOM_restart_CS), pointer :: CS
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr or with the
! specified variable name has been initialized from a restart file.
!
! Arguments: f_ptr - A pointer to the field that is being queried.
!  (in)      name - The name of the field that is being queried.
!  (in)      CS - The control structure returned by a previous call to
!                 restart_init.
  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (ASSOCIATED(CS%var_ptr1d(m)%p,f_ptr)) then
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

function query_initialized_2d_name(f_ptr, name, CS) result(query_initialized)
  real, dimension(:,:), target  :: f_ptr
  character(len=*)              :: name
  type(MOM_restart_CS), pointer :: CS
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr or with the
! specified variable name has been initialized from a restart file.
!
! Arguments: f_ptr - A pointer to the field that is being queried.
!  (in)      name - The name of the field that is being queried.
!  (in)      CS - The control structure returned by a previous call to
!                 restart_init.
  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (ASSOCIATED(CS%var_ptr2d(m)%p,f_ptr)) then
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

function query_initialized_3d_name(f_ptr, name, CS) result(query_initialized)
  real, dimension(:,:,:), target  :: f_ptr
  character(len=*)                :: name
  type(MOM_restart_CS),   pointer :: CS
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr or with the
! specified variable name has been initialized from a restart file.
!
! Arguments: f_ptr - A pointer to the field that is being queried.
!  (in)      name - The name of the field that is being queried.
!  (in)      CS - The control structure returned by a previous call to
!                 restart_init.
  integer :: m, n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (ASSOCIATED(CS%var_ptr3d(m)%p,f_ptr)) then
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

function query_initialized_4d_name(f_ptr, name, CS) result(query_initialized)
  real, dimension(:,:,:,:), target  :: f_ptr
  character(len=*)                  :: name
  type(MOM_restart_CS),     pointer :: CS
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr or with the
! specified variable name has been initialized from a restart file.
!
! Arguments: f_ptr - A pointer to the field that is being queried.
!  (in)      name - The name of the field that is being queried.
!  (in)      CS - The control structure returned by a previous call to
!                 restart_init.
  integer :: m, n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (ASSOCIATED(CS%var_ptr4d(m)%p,f_ptr)) then
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

subroutine save_restart(directory, time, G, CS, time_stamped, filename, GV)
!  save_restart saves all registered variables to restart files.
  character(len=*),        intent(in)    :: directory
  type(time_type),         intent(in)    :: time
  type(ocean_grid_type),   intent(inout) :: G
  type(MOM_restart_CS),    pointer       :: CS
  logical,          optional, intent(in) :: time_stamped
  character(len=*), optional, intent(in) :: filename
  type(verticalGrid_type), optional, intent(in) :: GV
! Arguments: directory - The directory where the restart file goes.
!  (in)      time - The time of this restart file.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 restart_init.
!  (in, opt) time_stamped - If true, the restart file names include
!                           a unique time stamp.  The default is false.
!  (in, opt) filename - A filename that overrides the name in CS%restartfile.
!
!  (in, opt) GV - The ocean's vertical grid structure.
  type(vardesc) :: vars(CS%max_fields)  ! Descriptions of the fields that
                                        ! are to be read from the restart file.
  type(fieldtype) :: fields(CS%max_fields) !
  character(len=200) :: restartpath     ! The restart file path (dir/file).
  character(len=80)  :: restartname     ! The restart file name (no dir).
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
  character(len=64) :: var_name         ! A variable's name.
  real :: restart_time
  character(len=32) :: filename_appendix = '' !fms appendix to filename for ensemble runs
  integer :: length

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "save_restart: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  ! With parallel read & write, it is possible to disable the following...

! jgj: this was set to 4294967292, changed to 4294967295 (see mpp_parameter.F90)
  if (CS%large_file_support) max_file_size = 4294967295_8

  num_files = 0
  next_var = 0
  nz = 1 ; if (present(GV)) nz = GV%ke

  call get_time(time,seconds,days)
  restart_time = real(days) + real(seconds)/86400.0

  restartname = trim(CS%restartfile)
  if (present(filename)) restartname = trim(filename)
  if (PRESENT(time_stamped)) then ; if (time_stamped) then
    call get_date(time,year,month,days,hour,minute,seconds)
    ! Compute the year-day, because I don't like months. - RWH
    do m=1,month-1
      days = days + days_in_month(set_date(year,m,2,0,0,0))
    enddo
    seconds = seconds + 60*minute + 3600*hour
    if (year <= 9999) then
      write(restartname,'("_Y",I4.4,"_D",I3.3,"_S",I5.5)') year, days, seconds
    else if (year <= 99999) then
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

    !query fms_io if there is a filename_appendix (for ensemble runs)
    call get_filename_appendix(filename_appendix)
    if(len_trim(filename_appendix) > 0) then
      length = len_trim(restartname)
      if(restartname(length-2:length) == '.nc') then
        restartname = restartname(1:length-3)//'.'//trim(filename_appendix)//'.nc'
      else
        restartname = restartname(1:length)  //'.'//trim(filename_appendix)
      end if
    end if

    restartpath = trim(directory)// trim(restartname)

    if (num_files < 10) then
      write(suffix,'("_",I1)') num_files
    else
      write(suffix,'("_",I2)') num_files
    endif

    if (num_files > 0) restartpath = trim(restartpath) // trim(suffix)

    do m=start_var,next_var-1
      vars(m-start_var+1) = CS%restart_field(m)%vars
    enddo
    call query_vardesc(vars(1), t_grid=t_grid, caller="save_restart")
    t_grid = adjustl(t_grid)
    if (t_grid(1:1) /= 'p') &
      call modify_vardesc(vars(1), t_grid='s', caller="save_restart")

    if (CS%parallel_restartfiles) then
      call create_file(unit, trim(restartpath), vars, (next_var-start_var), &
                       fields, MULTIPLE, G=G, GV=GV)
    else
      call create_file(unit, trim(restartpath), vars, (next_var-start_var), &
                       fields, SINGLE_FILE, G=G, GV=GV)
    endif

    do m=start_var,next_var-1

      if (ASSOCIATED(CS%var_ptr3d(m)%p)) then
        call write_field(unit,fields(m-start_var+1), G%Domain%mpp_domain, &
                         CS%var_ptr3d(m)%p, restart_time)
      elseif (ASSOCIATED(CS%var_ptr2d(m)%p)) then
        call write_field(unit,fields(m-start_var+1), G%Domain%mpp_domain, &
                         CS%var_ptr2d(m)%p, restart_time)
      elseif (ASSOCIATED(CS%var_ptr4d(m)%p)) then
        call write_field(unit,fields(m-start_var+1), G%Domain%mpp_domain, &
                         CS%var_ptr4d(m)%p, restart_time)
      elseif (ASSOCIATED(CS%var_ptr1d(m)%p)) then
        call write_field(unit, fields(m-start_var+1), CS%var_ptr1d(m)%p, &
                         restart_time)
      elseif (ASSOCIATED(CS%var_ptr0d(m)%p)) then
        call write_field(unit, fields(m-start_var+1), CS%var_ptr0d(m)%p, &
                         restart_time)
      endif
    enddo

    call close_file(unit)

    num_files = num_files+1

  enddo
end subroutine save_restart


subroutine restore_state(filename, directory, day, G, CS)
  character(len=*),      intent(in)  :: filename
  character(len=*),      intent(in)  :: directory
  type(time_type),       intent(out) :: day
  type(ocean_grid_type), intent(in)  :: G
  type(MOM_restart_CS),  pointer     :: CS
!    This subroutine reads the model state from previously
!  generated files.  All restart variables are read from the first
!  file in the input filename list in which they are found.

! Arguments: filename - A series of space delimited strings, each of
!                       which is either "r" or the name of a file
!                       from which the run is to be restarted.
!  (in)      directory - The directory where the restart or save
!                        files should be found.
!  (out)     day - The time of the restarted run.
!  (in)      G - The ocean's grid structure.
!  (in/out)  CS - The control structure returned by a previous call to
!                 restart_init.

  character(len=200) :: filepath  ! The path (dir/file) to the file being opened.
  character(len=80) :: fname      ! The name of the current file.
  character(len=8)  :: suffix     ! A suffix (like "_2") that is added to any
                                  ! additional restart files.
  character(len=256) :: mesg      ! A message for warnings.
  character(len=80) :: varname    ! A variable's name.
  integer :: num_restart     ! The number of restart files that have already
                             ! been opened.
  integer :: num_file        ! The number of files (restart files and others
                             ! explicitly in filename) that are open.
  integer :: start_char      ! The location of the starting character in the
                             ! current file name.
  integer :: n, m, start_of_day, num_days
  integer :: isL, ieL, jsL, jeL, is0, js0
  integer :: sizes(7)
  integer :: ndim, nvar, natt, ntime, pos
  integer :: unit(CS%max_fields) ! The mpp unit of all open files.
  logical :: unit_is_global(CS%max_fields) ! True if the file is global.
  character(len=8)   :: hor_grid ! Variable grid info.
  character(len=200) :: unit_path(CS%max_fields) ! The file names.
  logical :: fexists
  real, allocatable :: time_vals(:)
  type(fieldtype), allocatable :: fields(:)
  integer :: i, missing_fields
  real    :: t1, t2
  integer :: err
  character(len=32) :: filename_appendix = '' !fms appendix to filename for ensemble runs
  character(len=80) :: restartname
  integer :: length

  num_restart = 0 ; n = 1 ; start_char = 1
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "restore_state: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

! Get NetCDF ids for all of the restart files.
  do while (start_char <= len_trim(filename) )
    do m=start_char,len_trim(filename)
      if (filename(m:m) == ' ') exit
    enddo
    fname = filename(start_char:m-1)
    start_char = m
    do while ((start_char <= len_trim(filename)) .and. (filename(start_char:start_char) == ' '))
      start_char = start_char + 1
    enddo

    if ((fname(1:1)=='r') .and. ( len_trim(fname) == 1)) then
      err = 0
      if (num_restart > 0) err = 1 ! Avoid going through the file list twice.
      do while (err == 0)
        restartname = trim(CS%restartfile)

       !query fms_io if there is a filename_appendix (for ensemble runs)
       call get_filename_appendix(filename_appendix)
       if(len_trim(filename_appendix) > 0) then
         length = len_trim(restartname)
         if(restartname(length-2:length) == '.nc') then
           restartname = restartname(1:length-3)//'.'//trim(filename_appendix)//'.nc'
         else
           restartname = restartname(1:length)  //'.'//trim(filename_appendix)
         end if
        end if
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
          call open_file(unit(n), trim(filepath), READONLY_FILE, NETCDF_FILE, &
                         threading = MULTIPLE, fileset = SINGLE_FILE)
          unit_is_global(n) = .true.
        elseif (CS%parallel_restartfiles) then
          if (G%Domain%use_io_layout) then
            ! Look for decomposed files using the I/O Layout.
            fexists = file_exists(filepath, G%Domain)
            if (fexists) &
              call open_file(unit(n), trim(filepath), READONLY_FILE, NETCDF_FILE, &
                             domain=G%Domain%mpp_domain)
          else
            ! Look for any PE-specific files of the form NAME.nc.####.
            if (num_PEs()>10000) then
              write(filepath, '(a,i6.6)' ) trim(filepath)//'.', pe_here()
            else
              write(filepath, '(a,i4.4)' ) trim(filepath)//'.', pe_here()
            endif
            inquire(file=filepath, exist=fexists)
            if (fexists) &
              call open_file(unit(n), trim(filepath), READONLY_FILE, NETCDF_FILE, &
                             threading = MULTIPLE, fileset = SINGLE_FILE)
          endif
          if (fexists) unit_is_global(n) = .false.
        endif

        if (fexists) then
          unit_path(n) = filepath
          n = n + 1
          if (is_root_pe()) &
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
        call open_file(unit(n), trim(filepath), READONLY_FILE, NETCDF_FILE, &
                       threading = MULTIPLE, fileset = SINGLE_FILE)
        unit_is_global(n) = .true.
        unit_path(n) = filepath
        n = n + 1
        if (is_root_pe()) &
          call MOM_error(NOTE,"MOM_restart: MOM run restarted using : "//trim(filepath))
      else
        call MOM_error(WARNING,"MOM_restart: Unable to find restart file : "//trim(filepath))
      endif

    endif
  enddo ! while (start_char < strlen(filename)) loop
  num_file = n-1

  if (num_file == 0) then
    write(mesg,'("Unable to find any restart files specified by  ",A,"  in directory ",A,".")') &
                  trim(filename), trim(directory)
    call MOM_error(FATAL,"MOM_restart: "//mesg)
  endif

! Get the time from the first file in the list that has one.
  do n=1,num_file
    call get_file_info(unit(n), ndim, nvar, natt, ntime)
    if (ntime < 1) cycle

    allocate(time_vals(ntime))
    call get_file_times(unit(n), time_vals)
    t1 = time_vals(1)
    deallocate(time_vals)

    start_of_day = INT((t1 - INT(t1)) *86400) ! Number of seconds.
    num_days = INT(t1)
    day = set_time(start_of_day, num_days)
    exit
  enddo

  if (n>num_file) call MOM_error(WARNING,"MOM_restart: " // &
                                 "No times found in restart files.")

! Check the remaining files for different times and issue a warning
! if they differ from the first time.
  if (is_root_pe()) then
    do m = n+1,num_file
      call get_file_info(unit(n), ndim, nvar, natt, ntime)
      if (ntime < 1) cycle

      allocate(time_vals(ntime))
      call get_file_times(unit(n), time_vals)
      t2 = time_vals(1)
      deallocate(time_vals)

      if (t1 /= t2) then
        write(mesg,'("WARNING: Restart file ",I2," has time ",F10.4,"whereas &
         &simulation is restarted at ",F10.4," (differing by ",F10.4,").")')&
               m,t1,t2,t1-t2
        call MOM_error(WARNING, "MOM_restart: "//mesg)
      endif
    enddo
  endif

! Read each variable from the first file in which it is found.
  do n=1,num_file
    call get_file_info(unit(n), ndim, nvar, natt, ntime)

    allocate(fields(nvar))
    call get_file_fields(unit(n),fields(1:nvar))

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

      do i=1, nvar
        call get_file_atts(fields(i),name=varname)
        if (lowercase(trim(varname)) == lowercase(trim(CS%restart_field(m)%var_name))) then
          if (ASSOCIATED(CS%var_ptr1d(m)%p))  then
            ! Read a 1d array, which should be invariant to domain decomposition.
            call read_data(unit_path(n), varname, CS%var_ptr1d(m)%p, &
                           no_domain=.true., timelevel=1)
          elseif (ASSOCIATED(CS%var_ptr0d(m)%p)) then ! Read a scalar...
            call read_data(unit_path(n), varname, CS%var_ptr0d(m)%p, &
                           no_domain=.true., timelevel=1)
          elseif ((pos == 0) .and. ASSOCIATED(CS%var_ptr2d(m)%p)) then  ! Read a non-decomposed 2d array.
            ! Probably should query the field type to make sure that the sizes are right.
            call read_data(unit_path(n), varname, CS%var_ptr2d(m)%p, &
                           no_domain=.true., timelevel=1)
          elseif ((pos == 0) .and. ASSOCIATED(CS%var_ptr3d(m)%p)) then  ! Read a non-decomposed 3d array.
            ! Probably should query the field type to make sure that the sizes are right.
            call read_data(unit_path(n), varname, CS%var_ptr3d(m)%p, &
                           no_domain=.true., timelevel=1)
          elseif ((pos == 0) .and. ASSOCIATED(CS%var_ptr4d(m)%p)) then  ! Read a non-decomposed 4d array.
            ! Probably should query the field type to make sure that the sizes are right.
            call read_data(unit_path(n), varname, CS%var_ptr4d(m)%p, &
                           no_domain=.true., timelevel=1)
          elseif (unit_is_global(n) .or. G%Domain%use_io_layout) then
            if (ASSOCIATED(CS%var_ptr3d(m)%p)) then
              ! Read 3d array...  Time level 1 is always used.
              call read_data(unit_path(n), varname, CS%var_ptr3d(m)%p, &
                             G%Domain%mpp_domain, 1, position=pos)
            elseif (ASSOCIATED(CS%var_ptr2d(m)%p)) then ! Read 2d array...
              call read_data(unit_path(n), varname, CS%var_ptr2d(m)%p, &
                             G%Domain%mpp_domain, 1, position=pos)
            elseif (ASSOCIATED(CS%var_ptr4d(m)%p)) then ! Read 4d array...
              call read_data(unit_path(n), varname, CS%var_ptr4d(m)%p, &
                             G%Domain%mpp_domain, 1, position=pos)
            else
              call MOM_error(FATAL, "MOM_restart restore_state: "//&
                              "No pointers set for "//trim(varname))
            endif
          else ! Do not use an io_layout.  !### GET RID OF THIS BRANCH ONCE read_data_4d_new IS AVAILABLE.
            ! This file is decomposed onto the current processors.  We need
            ! to check whether the sizes look right, and abort if not.
            call get_file_atts(fields(i),ndim=ndim,siz=sizes)

            !   NOTE: The index ranges f var_ptrs always start with 1, so with
            ! symmetric memory the staggering is swapped from NE to SW!
            is0 = 1-G%isd
            if ((pos == EAST_FACE) .or. (pos == CORNER)) is0 = 1-G%IsdB
            if (sizes(1) == G%iec-G%isc+1) then
              isL = G%isc+is0 ; ieL = G%iec+is0
            elseif (sizes(1) == G%IecB-G%IscB+1) then
              isL = G%IscB+is0 ; ieL = G%IecB+is0
            elseif (((pos == EAST_FACE) .or. (pos == CORNER)) .and. &
                    (G%IscB == G%isc) .and. (sizes(1) == G%iec-G%isc+2)) then
              ! This is reading a symmetric file in a non-symmetric model.
              isL = G%isc-1+is0 ; ieL = G%iec+is0
            else
              call MOM_error(WARNING, "MOM_restart restore_state, "//trim(varname)//&
                    " has the wrong i-size in "//trim(filepath))
              exit
            endif

            js0 = 1-G%jsd
            if ((pos == NORTH_FACE) .or. (pos == CORNER)) js0 = 1-G%JsdB
            if (sizes(2) == G%jec-G%jsc+1) then
              jsL = G%jsc+js0 ; jeL = G%jec+js0
            elseif (sizes(2) == G%jecB-G%jscB+1) then
              jsL = G%jscB+js0 ; jeL = G%jecB+js0
            elseif (((pos == NORTH_FACE) .or. (pos == CORNER)) .and. &
                    (G%JscB == G%jsc) .and. (sizes(2) == G%jec-G%jsc+2)) then
              ! This is reading a symmetric file in a non-symmetric model.
              jsL = G%jsc-1+js0 ; jeL = G%jec+js0
            else
              call MOM_error(WARNING, "MOM_restart restore_state, "//trim(varname)//&
                    " has the wrong j-size in "//trim(filepath))
              exit
            endif

            if (ASSOCIATED(CS%var_ptr3d(m)%p)) then
              if (ntime == 0) then
                call read_field(unit(n), fields(i), &
                                CS%var_ptr3d(m)%p(isL:ieL,jsL:jeL,:))
              else
                call read_field(unit(n), fields(i), &
                                CS%var_ptr3d(m)%p(isL:ieL,jsL:jeL,:), 1)
              endif
            elseif (ASSOCIATED(CS%var_ptr2d(m)%p)) then
              if (ntime == 0) then
                call read_field(unit(n), fields(i), &
                                CS%var_ptr2d(m)%p(isL:ieL,jsL:jeL))
              else
                call read_field(unit(n), fields(i), &
                                CS%var_ptr2d(m)%p(isL:ieL,jsL:jeL), 1)
              endif
            elseif (ASSOCIATED(CS%var_ptr4d(m)%p)) then
              if (ntime == 0) then
                call read_field(unit(n), fields(i), &
                                CS%var_ptr4d(m)%p(isL:ieL,jsL:jeL,:,:))
              else
                call read_field(unit(n), fields(i), &
                                CS%var_ptr4d(m)%p(isL:ieL,jsL:jeL,:,:), 1)
              endif
            else
              call MOM_error(FATAL, "MOM_restart restore_state: "//&
                              "No pointers set for "//trim(varname))
            endif
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
    call close_file(unit(n))
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

subroutine restart_init(param_file, CS, restart_root)
  type(param_file_type), intent(in) :: param_file
  type(MOM_restart_CS),  pointer    :: CS
  character(len=*), optional, intent(in) :: restart_root
! Arguments: param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module.
!  (in,opt)  restart_root - A filename root that overrides the value in
!                           RESTARTFILE.  This will enable the use of this
!                           module by other components.
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_restart"   ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "restart_init called with an associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "PARALLEL_RESTARTFILES", &
                                CS%parallel_restartfiles, &
                 "If true, each processor writes its own restart file, \n"//&
                 "otherwise a single restart file is generated", &
                 default=.false.)

  if (present(restart_root)) then
    CS%restartfile = restart_root
    call log_param(param_file, mod, "RESTARTFILE from argument", CS%restartfile)
  else
    call get_param(param_file, mod, "RESTARTFILE", CS%restartfile, &
                 "The name-root of the restart file.", default="MOM.res")
  endif
  call get_param(param_file, mod, "LARGE_FILE_SUPPORT", CS%large_file_support, &
                 "If true, use the file-size limits with NetCDF large \n"//&
                 "file support (4Gb), otherwise the limit is 2Gb.", &
                 default=.true.)
  call get_param(param_file, mod, "MAX_FIELDS", CS%max_fields, &
                 "The maximum number of restart fields that can be used.", &
                 default=100)

  allocate(CS%restart_field(CS%max_fields))
  allocate(CS%var_ptr0d(CS%max_fields))
  allocate(CS%var_ptr1d(CS%max_fields))
  allocate(CS%var_ptr2d(CS%max_fields))
  allocate(CS%var_ptr3d(CS%max_fields))
  allocate(CS%var_ptr4d(CS%max_fields))

end subroutine restart_init

subroutine restart_init_end(CS)
  type(MOM_restart_CS),  pointer    :: CS

  if (associated(CS)) then
    if (CS%novars == 0) call restart_end(CS)
  endif

end subroutine restart_init_end

subroutine restart_end(CS)
  type(MOM_restart_CS),  pointer    :: CS

  if (associated(CS%restart_field)) deallocate(CS%restart_field)
  if (associated(CS%var_ptr0d)) deallocate(CS%var_ptr0d)
  if (associated(CS%var_ptr1d)) deallocate(CS%var_ptr1d)
  if (associated(CS%var_ptr2d)) deallocate(CS%var_ptr2d)
  if (associated(CS%var_ptr3d)) deallocate(CS%var_ptr3d)
  if (associated(CS%var_ptr4d)) deallocate(CS%var_ptr4d)
  deallocate(CS)

end subroutine restart_end

subroutine restart_error(CS)
  type(MOM_restart_CS),  pointer    :: CS
! Arguments: CS - A pointer that is set to point to the control structure
!                 for this module.  (Intent in.)
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

end module MOM_restart
