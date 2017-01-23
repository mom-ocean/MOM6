module BFB_initialization
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
!*    G%g_prime - The reduced gravity at each interface, in m s-2.     *
!*    G%Rlay - Layer potential density (coordinate variable) in kg m-3.*
!*  If SPONGE is defined:                                              *
!*    A series of subroutine calls are made to set up the damping      *
!*    rates and reference profiles for all variables that are damped   *
!*    in the sponge.                                                   *
!*                                                                     *
!*    These variables are all set in the set of subroutines (in this   *
!*  file) BFB_initialize_sponges_southonly and BFB_set_coord.          *
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
use MOM_variables, only : thermo_var_ptrs
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use MOM_verticalGrid, only : verticalGrid_type
implicit none ; private

#include <MOM_memory.h>

public BFB_set_coord
public BFB_initialize_sponges_southonly

logical :: first_call = .true.

contains

subroutine BFB_set_coord(Rlay, g_prime, GV, param_file, eqn_of_state)
! This subroutine specifies the vertical coordinate in terms of temperature at the surface and at the bottom. This case is set up in
! such a way that the temperature of the topmost layer is equal to the SST at the southern edge of the domain. The temperatures are
! then converted to densities of the top and bottom layers and linearly interpolated for the intermediate layers.
  real, dimension(NKMEM_), intent(out) :: Rlay, g_prime
  type(verticalGrid_type), intent(in)  :: GV
  type(param_file_type),   intent(in)  :: param_file
  type(EOS_type),          pointer     :: eqn_of_state
  real                                 :: drho_dt, SST_s, T_bot, rho_top, rho_bot
  integer                              :: k, nz
  character(len=40)  :: mod = "BFB_set_coord" ! This subroutine's name.

  call get_param(param_file, mod, "DRHO_DT", drho_dt, &
          "Rate of change of density with temperature.", &
           units="kg m-3 K-1", default=-0.2)
  call get_param(param_file, mod, "SST_S", SST_s, &
          "SST at the suothern edge of the domain.", units="C", default=20.0)
  call get_param(param_file, mod, "T_BOT", T_bot, &
                 "Bottom Temp", units="C", default=5.0)
  rho_top = GV%rho0 + drho_dt*SST_s
  rho_bot = GV%rho0 + drho_dt*T_bot
  nz = GV%ke

  !call MOM_error(FATAL, &
  ! "BFB_initialization.F90, BFB_set_coord: " // &
  ! "Unmodified user routine called - you must edit the routine to use it")
  do k = 1,nz
    Rlay(k) = (rho_bot - rho_top)/(nz-1)*real(k-1) + rho_top
    if (k >1) then
      g_prime(k) = (Rlay(k) - Rlay(k-1))*GV%g_earth/GV%rho0
    else
      g_prime(k) = GV%g_earth
    end if
    !Rlay(:) = 0.0
    !g_prime(:) = 0.0
  end do

  if (first_call) call write_BFB_log(param_file)

end subroutine BFB_set_coord

subroutine BFB_initialize_sponges_southonly(G, use_temperature, tv, param_file, CSp, h)
! This subroutine sets up the sponges for the southern bouundary of the domain. Maximum damping occurs within 2 degrees lat of the
! boundary. The damping linearly decreases northward over the next 2 degrees.
  type(ocean_grid_type), intent(in)                   :: G
  logical,               intent(in)                   :: use_temperature
  type(thermo_var_ptrs), intent(in)                   :: tv
  type(param_file_type), intent(in)                   :: param_file
  type(sponge_CS),       pointer                      :: CSp
  real, dimension(NIMEM_, NJMEM_, NKMEM_), intent(in) :: h
  !call MOM_error(FATAL, &
  ! "BFB_initialization.F90, BFB_initialize_sponges: " // &
  ! "Unmodified user routine called - you must edit the routine to use it")

  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1) ! A temporary array for eta.
  real :: Idamp(SZI_(G),SZJ_(G))    ! The inverse damping rate, in s-1.

  real :: H0(SZK_(G))
  real :: min_depth
  real :: damp, e_dense, damp_new, slat, wlon, lenlat, lenlon, nlat
  character(len=40)  :: mod = "BFB_initialize_sponges_southonly" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  eta(:,:,:) = 0.0 ; Idamp(:,:) = 0.0

!  Here the inverse damping time, in s-1, is set. Set Idamp to 0     !
!  wherever there is no sponge, and the subroutines that are called  !
!  will automatically set up the sponges only where Idamp is positive!
!  and mask2dT is 1.                                                   !

!   Set up sponges for DOME configuration
  call get_param(param_file, mod, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)

  call get_param(param_file, mod, "SOUTHLAT", slat, &
                 "The southern latitude of the domain.", units="degrees")
  call get_param(param_file, mod, "LENLAT", lenlat, &
                 "The latitudinal length of the domain.", units="degrees")
  call get_param(param_file, mod, "WESTLON", wlon, &
                 "The western longitude of the domain.", units="degrees", default=0.0)
  call get_param(param_file, mod, "LENLON", lenlon, &
                 "The longitudinal length of the domain.", units="degrees")
  nlat = slat + lenlat
  do k=1,nz ; H0(k) = -G%max_depth * real(k-1) / real(nz) ; enddo
!  do k=1,nz ; H0(k) = -G%max_depth * real(k-1) / real(nz-1) ; enddo ! Use for meridional thickness profile initialization
  do i=is,ie; do j=js,je
    if (G%geoLatT(i,j) < slat+2.0) then ; damp = 1.0
    elseif (G%geoLatT(i,j) < slat+4.0) then
       damp_new = 1.0*(slat+4.0-G%geoLatT(i,j))/2.0
    else ; damp = 0.0
    endif

    ! These will be streched inside of apply_sponge, so they can be in
    ! depth space for Boussinesq or non-Boussinesq models.

    ! This section is used for uniform thickness initialization
    do k = 1,nz; eta(i,j,k) = H0(k); enddo

    ! The below section is used for meridional temperature profile thickness initiation
    ! do k = 1,nz; eta(i,j,k) = H0(k); enddo
    ! if (G%geoLatT(i,j) > 40.0) then
    !   do k = 1,nz
    !     eta(i,j,k) = -G%Angstrom_z*(k-1)
    !   enddo
    ! elseif (G%geoLatT(i,j) > 20.0) then
    !   do k = 1,nz
    !     eta(i,j,k) = min(H0(k) + (G%geoLatT(i,j) - 20.0)*(G%max_depth - nz*G%Angstrom_z)/20.0, -(k-1)*G%angstrom_z)
    !   enddo
    ! endif
    eta(i,j,nz+1) = -G%max_depth

    if (G%bathyT(i,j) > min_depth) then
      Idamp(i,j) = damp/86400.0
    else ; Idamp(i,j) = 0.0 ; endif
  enddo ; enddo

!  This call sets up the damping rates and interface heights.
!  This sets the inverse damping timescale fields in the sponges.    !
  call initialize_sponge(Idamp, eta, G, param_file, CSp)

!   Now register all of the fields which are damped in the sponge.   !
! By default, momentum is advected vertically within the sponge, but !
! momentum is typically not damped within the sponge.                !

  if (first_call) call write_BFB_log(param_file)

end subroutine BFB_initialize_sponges_southonly

!> Write output about the parameter values being used.
subroutine write_BFB_log(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure indicating the
                                                  !! open file to parse for model
                                                  !! parameter values.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "BFB_initialization" ! This module's name.

  call log_version(param_file, mod, version)
  first_call = .false.

end subroutine write_BFB_log

end module BFB_initialization
