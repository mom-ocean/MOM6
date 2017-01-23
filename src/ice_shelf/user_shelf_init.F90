module user_shelf_init
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                        *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT *
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
!*    D - Basin depth in m.  (Must be positive.)                       *
!*    f - The Coriolis parameter, in s-1.                              *
!*    g - The reduced gravity at each interface, in m s-2.             *
!*    Rlay - Layer potential density (coordinate variable) in kg m-3.  *
!*  If TEMPERATURE is defined:                                         *
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
!*  in MOM_surface_forcing.F90.                                       *
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
!*  Macros written all in capital letters are defined in MOM_memory.h.*
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, f                                     *
!*    j+1  > o > o >   At ^:  v, tauy                                  *
!*    j    x ^ x ^ x   At >:  u, taux                                  *
!*    j    > o > o >   At o:  h, D, buoy, tr, T, S, Rml, ustar         *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

! use MOM_domains, only : sum_across_PEs
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_time_manager, only : time_type, set_time, time_type_to_real

use mpp_mod, only : mpp_pe, mpp_sync
! use MOM_io, only : close_file, fieldtype, file_exists
! use MOM_io, only : open_file, read_data, read_axis_data, SINGLE_FILE
! use MOM_io, only : write_field, slasher, vardesc
implicit none ; private

#include <MOM_memory.h>

public USER_initialize_shelf_mass, USER_update_shelf_mass
public USER_init_ice_thickness
logical :: first_call = .true.

type, public :: user_ice_shelf_CS ; private
  real :: Rho_ocean  ! The ocean's typical density, in kg m-3.
  real :: max_draft  ! The maximum ocean draft of the ice shelf, in m.
  real :: min_draft  ! The minimum ocean draft of the ice shelf, in m.
  real :: flat_shelf_width ! The range over which the shelf is min_draft thick.
  real :: shelf_slope_scale ! The range over which the shelf slopes.
  real :: pos_shelf_edge_0
  real :: shelf_speed
end type user_ice_shelf_CS

contains

subroutine USER_initialize_shelf_mass(mass_shelf, area_shelf_h, h_shelf, hmask, G, CS, param_file, new_sim)

  type(ocean_grid_type),            intent(in)  :: G
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: mass_shelf, area_shelf_h, hmask, h_shelf
  type(user_ice_shelf_CS),          pointer     :: CS
  type(param_file_type),            intent(in)  :: param_file
  logical                                       :: new_sim

! Arguments: mass_shelf - The mass per unit area averaged over the full ocean
!                         cell, in kg m-2. (Intent out)
!  (out)     area_shelf_h - The area of the ocean cell that is covered by the
!                           rigid ice shelf, in m2.
!  (in)      G          - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.


! just check for cvs
! This subroutine sets up the initial mass and area covered by the ice shelf.
  real :: Rho_ocean  ! The ocean's typical density, in kg m-3.
  real :: max_draft  ! The maximum ocean draft of the ice shelf, in m.
  real :: min_draft  ! The minimum ocean draft of the ice shelf, in m.
  real :: flat_shelf_width ! The range over which the shelf is min_draft thick.
  real :: c1 ! The maximum depths in m.
  character(len=40) :: mod = "USER_initialize_shelf_mass" ! This subroutine's name.
  integer :: i, j

  ! call MOM_error(FATAL, "USER_shelf_init.F90, USER_set_shelf_mass: " // &
  !  "Unmodified user routine called - you must edit the routine to use it")

  if (.not.associated(CS)) allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  if (first_call) call write_user_log(param_file)
  call get_param(param_file, mod, "RHO_0", CS%Rho_ocean, &
                 "The mean ocean density used with BOUSSINESQ true to \n"//&
                 "calculate accelerations and the mass for conservation \n"//&
                 "properties, or with BOUSSINSEQ false to convert some \n"//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0)
  call get_param(param_file, mod, "SHELF_MAX_DRAFT", CS%max_draft, &
                 units="m", default=1.0)
  call get_param(param_file, mod, "SHELF_MIN_DRAFT", CS%min_draft, &
                 units="m", default=1.0)
  call get_param(param_file, mod, "FLAT_SHELF_WIDTH", CS%flat_shelf_width, &
                 units="axis_units", default=0.0)
  call get_param(param_file, mod, "SHELF_SLOPE_SCALE", CS%shelf_slope_scale, &
                 units="axis_units", default=0.0)
  call get_param(param_file, mod, "SHELF_EDGE_POS_0", CS%pos_shelf_edge_0, &
                 units="axis_units", default=0.0)
  call get_param(param_file, mod, "SHELF_SPEED", CS%shelf_speed, &
                 units="axis_units day-1", default=0.0)

  call USER_update_shelf_mass(mass_shelf, area_shelf_h, h_shelf, hmask, G, CS, set_time(0,0), new_sim)


end subroutine USER_initialize_shelf_mass

subroutine USER_init_ice_thickness(h_shelf, area_shelf_h, hmask, G, param_file)
  type(ocean_grid_type),            intent(in)  :: G
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: area_shelf_h, hmask, h_shelf
  type(param_file_type),            intent(in)  :: param_file

  ! This subroutine initializes the ice shelf thickness.  Currently it does so
  ! calling USER_initialize_shelf_mass, but this can be revised as needed.
  real, dimension(SZI_(G),SZJ_(G)) :: mass_shelf
  type(user_ice_shelf_CS), pointer :: CS => NULL()

  call USER_initialize_shelf_mass(mass_shelf, area_shelf_h, h_shelf, hmask, G, CS, param_file, .true.)

end subroutine USER_init_ice_thickness

subroutine USER_update_shelf_mass(mass_shelf, area_shelf_h, h_shelf, hmask, G, CS, Time, new_sim)
  type(ocean_grid_type),            intent(in)    :: G
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: mass_shelf, area_shelf_h, hmask, h_shelf
  type(user_ice_shelf_CS),          pointer       :: CS
  type(time_type),                  intent(in)    :: Time
  logical,                          intent(in)    :: new_sim

! Arguments: mass_shelf - The mass per unit area averaged over the full ocean
!                         cell, in kg m-2. (Intent out)
!  (out)     area_shelf_h - The area of the ocean cell that is covered by the
!                           rigid ice shelf, in m2.
!  (in)      G          - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

  real :: c1, edge_pos, slope_pos
  integer :: i, j

  edge_pos = CS%pos_shelf_edge_0 + CS%shelf_speed*(time_type_to_real(Time) / 86400.0)

  slope_pos = edge_pos - CS%flat_shelf_width
  c1 = 0.0 ; if (CS%shelf_slope_scale > 0.0) c1 = 1.0 / CS%shelf_slope_scale


  do j=G%jsd,G%jed ;

   if (((j+G%jdg_offset) .le. G%domain%njglobal+G%domain%njhalo) .AND. &
       ((j+G%jdg_offset) .ge. G%domain%njhalo+1)) then

    do i=G%isc,G%iec

!    if (((i+G%idg_offset) <= G%domain%niglobal+G%domain%nihalo) .AND. &
!           ((i+G%idg_offset) >= G%domain%nihalo+1)) then

    if ((j.ge.G%jsc) .and. (j.le.G%jec)) then

      if (new_sim) then ; if (G%geoLonCu(i-1,j) >= edge_pos) then
        ! Everything past the edge is open ocean.
        mass_shelf(i,j) = 0.0
        area_shelf_h(i,j) = 0.0
        hmask (i,j) = 0.0
        h_shelf (i,j) = 0.0
      else
        if (G%geoLonCu(i,j) > edge_pos) then
          area_shelf_h(i,j) = G%areaT(i,j) * (edge_pos - G%geoLonCu(i-1,j)) / &
                              (G%geoLonCu(i,j) - G%geoLonCu(i-1,j))
          hmask (i,j) = 2.0
        else
          area_shelf_h(i,j) = G%areaT(i,j)
          hmask (i,j) = 1.0
        endif

        if (G%geoLonT(i,j) > slope_pos) then
          h_shelf (i,j) = CS%min_draft
          mass_shelf(i,j) = CS%Rho_ocean * CS%min_draft
        else
          mass_shelf(i,j) = CS%Rho_ocean * (CS%min_draft + &
                 (CS%max_draft - CS%min_draft) * &
                 min(1.0, (c1*(slope_pos - G%geoLonT(i,j)))**2) )
          h_shelf(i,j) = (CS%min_draft + &
                 (CS%max_draft - CS%min_draft) * &
                 min(1.0, (c1*(slope_pos - G%geoLonT(i,j)))**2) )
        endif

    endif ; endif ; endif

    if ((i+G%idg_offset) .eq. G%domain%nihalo+1) then
      hmask(i-1,j) = 3.0
    endif

  enddo ; endif ; enddo

end subroutine USER_update_shelf_mass

subroutine write_user_log(param_file)
  type(param_file_type), intent(in) :: param_file

  character(len=128) :: version = '$Id: user_shelf_init.F90,v 1.1.2.7 2012/06/19 22:15:52 Robert.Hallberg Exp $'
  character(len=128) :: tagname = '$Name: MOM_ogrp $'
  character(len=40)  :: mod = "user_shelf_init" ! This module's name.

  call log_version(param_file, mod, version, tagname)
  first_call = .false.

end subroutine write_user_log

end module user_shelf_init
