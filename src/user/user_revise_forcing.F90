module user_revise_forcing
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
!
!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  This module provides a method for updating the forcing fluxes      *
!*  using user-written code without the need to duplicate the          *
!*  extensive code used to create or obtain the fluxes.                *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**
use MOM_domains, only : pass_var, pass_vector, AGRID
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, read_data
use MOM_restart, only : register_restart_field, MOM_restart_CS
use MOM_time_manager, only : time_type, operator(+), operator(/), get_time
use MOM_tracer_flow_control, only : call_tracer_set_forcing
use MOM_tracer_flow_control, only : tracer_flow_control_CS
use MOM_variables, only : surface

implicit none ; private

public user_alter_forcing, user_revise_forcing_init

type, public :: user_revise_forcing_CS ; private
  real    :: cdrag               ! The quadratic bottom drag coefficient.
end type user_revise_forcing_CS

contains

!> This subroutine sets the surface wind stresses.
subroutine user_alter_forcing(state, fluxes, day, G, CS)
  type(surface),            intent(in)    :: state  !< A structure containing fields that
                                                    !! describe the surface state of the ocean.
  type(forcing),            intent(inout) :: fluxes !< A structure containing pointers to any
                                                    !! possible forcing fields. Unused fields
                                                    !! have NULL ptrs.
  type(time_type),          intent(in)    :: day    !< Time of the fluxes.
  type(ocean_grid_type),    intent(in)    :: G      !< The ocean's grid structure.
  type(user_revise_forcing_CS), pointer   :: CS     !< A pointer to the control structure
                                                    !! returned by a previous call to
                                                    !! surface_forcing_init.

end subroutine user_alter_forcing

subroutine user_revise_forcing_init(param_file,CS)
  type(param_file_type), intent(in) :: param_file   !< !< A structure indicating the open file to
                                                    !! parse for model parameter values.
  type(user_revise_forcing_CS), pointer   :: CS     !< A pointer to the control structure
                                                    !! returned by a previous call to
                                                    !! surface_forcing_init.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "user_revise_forcing" ! This module's name.

  call log_version(param_file, mod, version)

end subroutine user_revise_forcing_init

end module user_revise_forcing
