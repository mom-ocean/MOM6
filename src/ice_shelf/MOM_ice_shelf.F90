module MOM_ice_shelf
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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                              !
!   This is a null version of the ice shelf code with the same interfaces as   !
! the active version.                                                          !
!                                                                              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : read_param, param_file_type, log_param, log_version
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_initialization, only : Get_MOM_Input
use MOM_io, only : write_version_number
use MOM_time_manager, only : time_type, set_time, time_type_to_real
use MOM_variables, only : directories, surface

implicit none ; private

public shelf_calc_flux, add_shelf_flux, initialize_ice_shelf, ice_shelf_end
public ice_shelf_save_restart, solo_time_step

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
type, public :: ice_shelf_CS ; private
  logical :: isshelf=.false.  ! True if a shelf model is to be used.
end type ice_shelf_CS

contains

subroutine shelf_calc_flux(state, fluxes, Time, time_step, CS)
  type(surface),         intent(inout) :: state
  type(forcing),         intent(inout) :: fluxes
  type(time_type),       intent(in)    :: Time
  real,                  intent(in)    :: time_step
  type(ice_shelf_CS),    pointer       :: CS

! This sets up the fluxes between the ocean and an ice-shelf.
!
! Arguments: state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (inout)   fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      Time - Start time of the fluxes.
!  (in)      time_step - Length of time over which these fluxes
!                        will be applied, in s.
!  (in)      CS - A pointer to the control structure returned by a previous
!                 call to initialize_ice_shelf.

  if (.not. associated(CS)) call MOM_error(FATAL, "shelf_calc_flux: "// &
       "initialize_ice_shelf must be called before shelf_calc_flux.")

  ! In the null version of the code, this subroutine does nothing.

end subroutine shelf_calc_flux

subroutine add_shelf_flux(G, CS, state, fluxes)
  type(ocean_grid_type),              intent(in)    :: G
  type(ice_shelf_CS),                 intent(in)    :: CS
  type(surface),                      intent(in)    :: state
  type(forcing),                      intent(inout) :: fluxes
! Arguments:
!  (in)      fluxes - A structure of surface fluxes that may be used.
!  (in)      visc - A structure containing vertical viscosities, bottom boundary
!                   layer properies, and related fields.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - This module's control structure.

  ! In the null version of the code, this subroutine does nothing.

end subroutine add_shelf_flux

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! shelf_model_init - initializes shelf model data, parameters and diagnostics  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine initialize_ice_shelf(Time, CS, fluxes, Time_in, solo_mode_in)
  type(time_type),           intent(inout) :: Time
  type(ice_shelf_CS),        pointer       :: CS
  type(forcing),   optional, intent(inout) :: fluxes
  type(time_type), optional, intent(in)    :: Time_in
  logical,         optional, intent(in)    :: solo_mode_in

  ! In the null version of the code, this subroutine does nothing.

  type(param_file_type) :: param_file
  type(directories)  :: dirs
  logical :: isshelf ! True if a shelf model is to be used.
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'
  character(len=40)  :: mod = "MOM_ice_shelf"  ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "shelf_model_init called with an associated "// &
                             "control structure.")
    return
  endif
    
  !   Go through all of the infrastructure initialization calls, since this is
  ! being treated as an independent component that just happens to use the
  ! MOM's grid and infrastructure.
  call Get_MOM_Input(param_file, dirs)
  isshelf = .false. ; call read_param(param_file,"ICE_SHELF",isshelf)
  if (.not.isshelf) return

  call log_version(param_file, mod, version, tagname, "")

  call MOM_error(FATAL, "This version of the code does not have an active "//&
                         "ice shelf available.")

end subroutine initialize_ice_shelf

subroutine ice_shelf_save_restart(CS, Time, directory, time_stamped, filename_suffix)
  type(ice_shelf_CS),         pointer    :: CS
  type(time_type),            intent(in) :: Time
  character(len=*), optional, intent(in) :: directory
  logical,          optional, intent(in) :: time_stamped
  character(len=*), optional, intent(in) :: filename_suffix

! Arguments: CS - A structure containing the internal ocean state (in).
!  (in)      Time - The model time at this call.  This is needed for mpp_write calls.
!  (in, opt) directory - An optional directory into which to write these restart files.
!  (in, opt) time_stamped - If true, the restart file names include
!                           a unique time stamp.  The default is false.
!  (in, opt) filename_suffix - An optional suffix (e.g., a time-stamp) to append
!                              to the restart file names.

  ! In the null version of the code, this subroutine does nothing.

end subroutine ice_shelf_save_restart

subroutine ice_shelf_end(CS)
  type(ice_shelf_CS), pointer   :: CS
  ! This subroutine deallocates all memory associated with this module.

  if (.not.associated(CS)) return

  ! In the null version of the code, this subroutine does nothing.

  deallocate(CS)

end subroutine ice_shelf_end

subroutine solo_time_step(CS, time_step, n, Time, min_time_step_in)
  type(ice_shelf_CS), pointer       :: CS
  real,               intent(in)    :: time_step
  integer,            intent(inout) :: n
  type(time_type),    intent(inout) :: Time
  real, optional,     intent(in)    :: min_time_step_in

  ! In the null version of the code, this subroutine does nothing.

end subroutine solo_time_step

end module MOM_ice_shelf
