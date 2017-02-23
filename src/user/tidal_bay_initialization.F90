module tidal_bay_initialization
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

public tidal_bay_set_OBC_data

contains

!> This subroutine sets the properties of flow at open boundary conditions.
subroutine tidal_bay_set_OBC_data(OBC, G, h, Time)
  type(ocean_OBC_type),   pointer    :: OBC  !< This open boundary condition type specifies
                                             !! whether, where, and what open boundary
                                             !! conditions are used.
  type(ocean_grid_type),  intent(in) :: G    !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in) :: h !< layer thickness.
  type(time_type),        intent(in) :: Time !< model time.

  ! The following variables are used to set up the transport in the tidal_bay example.
  real :: time_sec, cff, cff2, tide_flow
  real :: my_area, my_flux
  real :: PI
  character(len=40)  :: mod = "tidal_bay_set_OBC_data" ! This subroutine's name.
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, nz, n
  integer :: IsdB, IedB, JsdB, JedB
  type(OBC_segment_type), pointer :: segment

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  PI = 4.0*atan(1.0) ;

  if (.not.associated(OBC)) return

  time_sec = time_type_to_real(Time)
  cff = 0.1*sin(2.0*PI*time_sec/(12.0*3600.0))
  tide_flow = 3.0e6
  my_area=0.0
  my_flux=0.0
  do j=jsd,jed ; do I=IsdB,IedB
    if (OBC%OBC_segment_u(I,j) /= OBC_NONE) then
      do k=1,nz
        cff2 = h(I,j,k)*G%dyCu(I,j)
        my_area = my_area + cff2
      enddo
    endif
  enddo ; enddo
  my_flux = -tide_flow*SIN(2.0*PI*time_sec/(12.0*3600.0))

  ! Old way
  segment => OBC%segment(1)
  do j=jsd,jed ; do I=IsdB,IedB
    if (OBC%OBC_segment_u(I,j) /= OBC_NONE) then
      OBC%eta_outer_u(I,j) = cff
      OBC%ubt_outer(I,j) = my_flux/my_area
      if (segment%nudged) then
        do k=1,nz
          OBC%u(I,j,k) = my_flux/my_area
          OBC%uh(I,j,k) = 0.5*OBC%u(I,j,k)*(h(i,j,k) + h(i+1,j,k))
        enddo
      endif
    endif
  enddo ; enddo
  do J=JsdB,JedB ; do i=isd,ied
    if (OBC%OBC_segment_v(i,J) /= OBC_NONE) then
      OBC%eta_outer_v(i,J) = cff
      OBC%vbt_outer(i,J) = 0.0
    endif
  enddo ; enddo

  ! New way
  do n = 1, OBC%number_of_segments
    segment => OBC%segment(n)

    if (.not. segment%on_pe) cycle

    segment%normal_vel_bt(:,:) = my_flux/my_area
    segment%eta(:,:) = cff

  enddo ! end segment loop

end subroutine tidal_bay_set_OBC_data

!> \class tidal_bay_Initialization
!!
!! The module configures the model for the "tidal_bay" experiment.
!! tidal_bay = Tidally resonant bay from Zygmunt Kowalik's class on tides.
end module tidal_bay_initialization
