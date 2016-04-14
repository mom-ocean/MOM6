module user_initialization
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
!*    This subroutine initializes the fields for the simulations.      *
!*  The one argument passed to initialize, Time, is set to the         *
!*  current time of the simulation.  The fields which are initialized  *
!*  here are:                                                          *
!*    u - Zonal velocity in m s-1.                                     *
!*    v - Meridional velocity in m s-1.                                *
!*    h - Layer thickness in m.  (Must be positive.)                   *
!*    G%bathyT - Basin depth in m.  (Must be positive.)                *
!*    G%CoriolisBu - The Coriolis parameter, in s-1.                   *
!*    GV%g_prime - The reduced gravity at each interface, in m s-2.    *
!*    GV%Rlay - Layer potential density (coordinate variable), kg m-3. *
!*  If ENABLE_THERMODYNAMICS is defined:                               *
!*    T - Temperature in C.                                            *
!*    S - Salinity in psu.                                             *
!*  If BULKMIXEDLAYER is defined:                                      *
!*    Rml - Mixed layer and buffer layer potential densities in        *
!*          units of kg m-3.                                           *
!*  If SPONGE is defined:                                              *
!*    A series of subroutine calls are made to set up the damping      *
!*    rates and reference profiles for all variables that are damped   *
!*    in the sponge.                                                   *
!*  Any user provided tracer code is also first linked through this    *
!*  subroutine.                                                        *
!*                                                                     *
!*    Forcing-related fields (taux, tauy, buoy, ustar, etc.) are set   *
!*  in MOM_surface_forcing.F90.                                        *
!*                                                                     *
!*    These variables are all set in the set of subroutines (in this   *
!*  file) USER_initialize_bottom_depth, USER_initialize_thickness,     *
!*  USER_initialize_velocity,  USER_initialize_temperature_salinity,   *
!*  USER_initialize_mixed_layer_density, USER_initialize_sponges,      *
!*  USER_set_coord, and USER_set_ref_profile.                          *
!*                                                                     *
!*    The names of these subroutines should be self-explanatory. They  *
!*  start with "USER_" to indicate that they will likely have to be    *
!*  modified for each simulation to set the initial conditions and     *
!*  boundary conditions.  Most of these take two arguments: an integer *
!*  argument specifying whether the fields are to be calculated        *
!*  internally or read from a NetCDF file; and a string giving the     *
!*  path to that file.  If the field is initialized internally, the    *
!*  path is ignored.                                                   *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, CoriolisBu                            *
!*    j+1  > o > o >   At ^:  v, tauy                                  *
!*    j    x ^ x ^ x   At >:  u, taux                                  *
!*    j    > o > o >   At o:  h, bathyT, buoy, tr, T, S, Rml, ustar    *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : close_file, create_file, fieldtype, file_exists
use MOM_io, only : open_file, read_data, read_axis_data, SINGLE_FILE
use MOM_io, only : write_field, slasher
use MOM_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use MOM_tracer_registry, only : tracer_registry_type, add_tracer_OBC_values
use MOM_variables, only : thermo_var_ptrs, ocean_OBC_type, OBC_NONE, OBC_SIMPLE
use MOM_variables, only : OBC_FLATHER_E, OBC_FLATHER_W, OBC_FLATHER_N, OBC_FLATHER_S
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
implicit none ; private

#include <MOM_memory.h>

public USER_set_coord, USER_initialize_topography, USER_initialize_thickness
public USER_initialize_velocity, USER_init_temperature_salinity
public USER_init_mixed_layer_density, USER_initialize_sponges
public USER_set_Open_Bdry_Conds, USER_set_rotation

logical :: first_call = .true.

contains

subroutine USER_set_coord(Rlay, g_prime, G, param_file, eqn_of_state)
  type(ocean_grid_type), intent(in) :: G
  real, dimension(:), intent(out) :: Rlay, g_prime
  type(param_file_type), intent(in) :: param_file
  type(EOS_type),        pointer    :: eqn_of_state

  call MOM_error(FATAL, &
   "USER_initialization.F90, USER_set_coord: " // &
   "Unmodified user routine called - you must edit the routine to use it")
  Rlay(:) = 0.0
  g_prime(:) = 0.0
  
  if (first_call) call write_user_log(param_file)

end subroutine USER_set_coord

subroutine USER_initialize_topography(D, G, param_file)
  type(ocean_grid_type), intent(in) :: G
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: D
  type(param_file_type), intent(in) :: param_file
  call MOM_error(FATAL, &
   "USER_initialization.F90, USER_initialize_topography: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  D(:,:) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_topography

subroutine USER_initialize_thickness(h, G, param_file, T)
  type(ocean_grid_type), intent(in) :: G
  real, intent(out), dimension(SZI_(G),SZJ_(G), SZK_(G)) :: h
  type(param_file_type), intent(in) :: param_file
  real, intent(in), dimension(SZI_(G),SZJ_(G), SZK_(G))  :: T
  call MOM_error(FATAL, &
   "USER_initialization.F90, USER_initialize_thickness: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  h(:,:,1) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_thickness

subroutine USER_initialize_velocity(u, v, G, param_file)
  type(ocean_grid_type),                       intent(in)  :: G
  real, dimension(SZIB_(G), SZJ_(G), SZK_(G)), intent(out) :: u
  real, dimension(SZI_(G), SZJB_(G), SZK_(G)), intent(out) :: v
  type(param_file_type),                       intent(in)  :: param_file
  call MOM_error(FATAL, &
   "USER_initialization.F90, USER_initialize_velocity: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  u(:,:,1) = 0.0
  v(:,:,1) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_velocity

subroutine USER_init_temperature_salinity(T, S, G, param_file, eqn_of_state)
  type(ocean_grid_type),                      intent(in)  :: G
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)), intent(out) :: T, S
  type(param_file_type),                      intent(in)  :: param_file
  type(EOS_type),                             pointer     :: eqn_of_state
  call MOM_error(FATAL, &
   "USER_initialization.F90, USER_init_temperature_salinity: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  T(:,:,1) = 0.0
  S(:,:,1) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_init_temperature_salinity

subroutine USER_init_mixed_layer_density(Rml, G, param_file, use_temperature, &
                                         eqn_of_state, T, S, P_Ref)
  type(ocean_grid_type),                         intent(in)  :: G
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)),    intent(out) :: Rml
  type(param_file_type),                         intent(in)  :: param_file
  logical,                                       intent(in)  :: use_temperature
  type(EOS_type),                      optional, pointer     :: eqn_of_state
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)), optional, intent(in) :: T, S
  real,                                optional, intent(in)  :: P_Ref
  call MOM_error(FATAL, &
   "USER_initialization.F90, USER_init_mixed_layer_density: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  Rml(:,:,1) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_init_mixed_layer_density

subroutine USER_initialize_sponges(G, use_temperature, tv, param_file, CSp, h)
  type(ocean_grid_type), intent(in) :: G
  logical,               intent(in) :: use_temperature
  type(thermo_var_ptrs), intent(in) :: tv
  type(param_file_type), intent(in) :: param_file
  type(sponge_CS),       pointer    :: CSp
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)), intent(in) :: h
  call MOM_error(FATAL, &
   "USER_initialization.F90, USER_initialize_sponges: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_sponges

subroutine USER_set_Open_Bdry_Conds(OBC, tv, G, param_file, tr_Reg)
  type(ocean_OBC_type),       pointer    :: OBC
  type(thermo_var_ptrs),      intent(in) :: tv
  type(ocean_grid_type),      intent(in) :: G
  type(param_file_type),      intent(in) :: param_file
  type(tracer_registry_type), pointer    :: tr_Reg
  call MOM_error(FATAL, &
   "USER_initialization.F90, USER_set_Open_Bdry_Conds: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  if (first_call) call write_user_log(param_file)

end subroutine USER_set_Open_Bdry_Conds

subroutine USER_set_rotation(G, param_file)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
  call MOM_error(FATAL, &
   "USER_initialization.F90, USER_set_rotation: " // &
   "Unmodified user routine called - you must edit the routine to use it")

  if (first_call) call write_user_log(param_file)

end subroutine USER_set_rotation

subroutine write_user_log(param_file)
  type(param_file_type), intent(in) :: param_file

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "user_initialization" ! This module's name.

  call log_version(param_file, mod, version)
  first_call = .false.

end subroutine write_user_log

end module user_initialization
