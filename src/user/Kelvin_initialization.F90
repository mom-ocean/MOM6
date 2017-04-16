module Kelvin_initialization
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

use MOM_dyn_horgrid,    only : dyn_horgrid_type
use MOM_error_handler,  only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_grid,           only : ocean_grid_type
use MOM_open_boundary,  only : ocean_OBC_type, OBC_NONE
use MOM_open_boundary,  only : OBC_segment_type
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_time_manager,   only : time_type, set_time, time_type_to_real

implicit none ; private

#include <MOM_memory.h>

public Kelvin_set_OBC_data
public Kelvin_initialize

contains

!> This subroutine sets the properties of flow at open boundary conditions.
subroutine Kelvin_initialize(OBC, G, h, param_file)
  type(ocean_OBC_type),   pointer    :: OBC  !< This open boundary condition type specifies
                                             !! whether, where, and what open boundary
                                             !! conditions are used.
  type(ocean_grid_type),  intent(in) :: G    !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in) :: h !< layer thickness.
  type(param_file_type),  intent(in) :: param_file !< Parameter file structure

  ! The following variables are used to set up the transport in the Kelvin example.
  real :: time_sec, cff
  real :: PI
  character(len=40)  :: mod = "Kelvin_initialize" ! This subroutine's name.
  integer :: i, j, k, n, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB
  integer :: mode
  real    :: fac, omega, x, y

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) call MOM_error(FATAL, 'Kelvin_initialization.F90: '// &
        'Kelvin_set_OBC_data() was called but OBC type was not initialized!')

  call get_param(param_file, mod, "KELVIN_WAVE_MODE", mode, &
                 "Vertical Kelvin wave mode imposed at upstream open boundary.", &
                 default=0)

end subroutine Kelvin_initialize

!> This subroutine sets the properties of flow at open boundary conditions.
subroutine Kelvin_set_OBC_data(OBC, G, h, Time)
  type(ocean_OBC_type),   pointer    :: OBC  !< This open boundary condition type specifies
                                             !! whether, where, and what open boundary
                                             !! conditions are used.
  type(ocean_grid_type),  intent(in) :: G    !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in) :: h !< layer thickness.
  type(time_type),        intent(in) :: Time !< model time.

  ! The following variables are used to set up the transport in the Kelvin example.
  real :: time_sec, cff
  real :: PI
  character(len=40)  :: mod = "Kelvin_set_OBC_data" ! This subroutine's name.
  integer :: i, j, k, n, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB
  integer :: mode = 0
  real    :: fac, omega, x, y
  real    :: val1, val2
  type(OBC_segment_type), pointer :: segment

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) call MOM_error(FATAL, 'Kelvin_initialization.F90: '// &
        'Kelvin_set_OBC_data() was called but OBC type was not initialized!')

  time_sec = time_type_to_real(Time)
  PI = 4.0*atan(1.0)

  if (mode == 0) then
    fac = 1.0
    omega = 2.0 * PI / (12.42 * 3600.0)      ! M2 Tide period
    val1 = sin(omega * time_sec)
  else
    omega = 2.0 * PI / (12.42 * 3600.0)      ! M2 Tide period
  endif

  do n=1,OBC%number_of_segments
    segment => OBC%segment(n)
    if (.not. segment%on_pe) cycle
    ! Apply values to the inflow end only.
    if (segment%Is_obc /= 4) cycle

    isdB = segment%HI%isdB ; iedB = segment%HI%iedB
    Jsd = segment%HI%Jsd ; Jed = segment%HI%Jed
    do J=Jsd,Jed ; do i=isdB,iedB
      x = 0.5 * 1000. * (G%geoLonT(i+1,j) + G%geoLonT(i,j))
      y = 0.5 * 1000. * (G%geoLatT(i+1,j) + G%geoLatT(i,j))
      if (mode == 0) then
        cff = sqrt(G%g_Earth * 0.5 * (G%bathyT(i+1,j) + G%bathyT(i,j)))
        val2 = fac * exp(- 0.5 * (G%CoriolisBu(I,J) + G%CoriolisBu(I,J-1)) * y / cff)
        OBC%eta_outer_u(I,j) = val2 * cos(omega * time_sec)
        OBC%ubt_outer(I,j) = val1 * cff / (0.5 *              &
               (G%bathyT(i+1,j) + G%bathyT(i,j))) * val2
! New way, not yet
!       segment%eta(I,j) = val2 * cos(omega * time_sec)
!       segment%normal_vel_bt(I,j) = val1 * cff / (0.5 *      &
!                (G%bathyT(i+1,j) + G%bathyT(i,j))) * val2
      else
      endif
    enddo ; enddo
  enddo

end subroutine Kelvin_set_OBC_data

!> \class Kelvin_Initialization
!!
!! The module configures the model for the Kelvin wave experiment.
!! Kelvin = coastally-trapped Kelvin waves from the ROMS examples.
!! Initialize with level surfaces and drive the wave in at the west,
!! radiate out at the east.
end module Kelvin_initialization
