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
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_time_manager,   only : time_type, set_time, time_type_to_real

implicit none ; private

#include <MOM_memory.h>

public Kelvin_initialize_thickness
public Kelvin_initialize_velocity 
public Kelvin_set_OBC_data

contains

!> Initialization of thicknesses in Kelvin wave test
subroutine Kelvin_initialize_thickness(h, G)
  type(ocean_grid_type),   intent(in) :: G                 !< Grid structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: h !< Thickness

  integer :: i, j, k, is, ie, js, je, nz
  real    :: x, y
  real    :: val1, val2, PI
  character(len=40) :: verticalCoordinate

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call MOM_mesg("Kelvin_initialization.F90, Kelvin_initialize_thickness: setting thickness")

  PI = 4.0*atan(1.0) ;
  val1 = 1.0                             ! zeta0
  val2 = 2.0 * PI / (12.42 * 3600.0)       ! M2 Tide period

! do j = G%jsc,G%jec ; do i = G%isc,G%iec
!   do k = 1, nz
!     x = G%geoLonT(i,j)
!     y = G%geoLatT(i,j)
!     h(i,j,k)=val1 * exp(-G%CoriolisBu(i,j) * y /        &   
!                  sqrt(G%G_earth * G%depthT(i,j))) *     &   
!                cos(val2 * x / sqrt(g * G%depthT(i,j)))
!   enddo
! end do ; end do
  h(:,:,:) = 0.0

end subroutine Kelvin_initialize_thickness


!> Initialization of u and v in the Kelvin wave test
subroutine Kelvin_initialize_velocity(u, v, h, G)
  type(ocean_grid_type),                  intent(in)     :: G  !< Grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(out) :: u  !< i-component of velocity [m/s]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(out) :: v  !< j-component of velocity [m/s]
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(in)  :: h  !< Thickness [H]

  real    :: x, y, PI
  real    :: val1, val2
  integer :: i, j, k, is, ie, js, je, nz
  
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  PI = 4.0*atan(1.0) ;
  val1 = 1.0                             ! zeta0
  val2 = 2.0 * PI / (12.42*3600.0)       ! M2 Tide period

  v(:,:,:) = 0.0
  u(:,:,:) = 0.0
  
! do j = G%jsc,G%jec ; do I = G%isc-1,G%iec+1
!   do k = 1, nz
!     x = 0.5*(G%geoLonT(i+1,j)+G%geoLonT(i,j))
!     y = 0.5*(G%geoLatT(i+1,j)+G%geoLatT(i,j))
!     u(I,j,k) = ...
!   enddo
! enddo ; enddo

end subroutine Kelvin_initialize_velocity

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
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB
  real    :: fac, omega, x, y
  real    :: val1, val2

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) return

  time_sec = time_type_to_real(Time)
  PI = 4.0*atan(1.0)
  fac = 1.0
  omega = 2.0 * PI / (12.42 * 3600.0)      ! M2 Tide period
  val1 = fac * sin(omega * time_sec)

  do j=jsd,jed ; do I=IsdB,IedB
    if (OBC%OBC_segment_u(I,j) /= OBC_NONE) then
      x = 0.5 * (G%geoLonT(i-1,j) + G%geoLonT(i,j))
      y = 0.5 * (G%geoLatT(i-1,j) + G%geoLatT(i,j))
      cff = sqrt(G%g_Earth * 0.5 * (G%bathyT(i-1,j) + G%bathyT(i,j)))
      val2 = fac * exp(-G%CoriolisBu(i,j) * y / cff)
      OBC%eta_outer_u(I,j) = val2 * cos(omega * time_sec)
      OBC%ubt_outer(I,j) = val1 * cff / 0.5 *            &
             (G%bathyT(i-1,j) + G%bathyT(i,j)) * val2
    endif
  enddo ; enddo

  do J=JsdB,JedB ; do i=isd,ied
    if (OBC%OBC_segment_v(i,J) /= OBC_NONE) then
      x = 0.5 * (G%geoLonT(i-1,j+1) + G%geoLonT(i-1,j))
      y = 0.5 * (G%geoLatT(i-1,j+1) + G%geoLatT(i-1,j))
      cff = sqrt(G%g_Earth * 0.5 * (G%bathyT(i-1,j+1) + G%bathyT(i-1,j)))
      val2 = fac * exp(-G%CoriolisBu(i,j) * y / cff)
      OBC%eta_outer_v(i,J) = cff
      OBC%vbt_outer(i,J) = 0.0
    endif
  enddo ; enddo

end subroutine Kelvin_set_OBC_data

!> \class Kelvin_Initialization
!!
!! The module configures the model for the Kelvin wave experiment.
!! Kelvin = coastally-trapped Kelvin waves from the ROMS examples.
end module Kelvin_initialization
