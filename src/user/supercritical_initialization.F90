module supercritical_initialization
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
use MOM_open_boundary,  only : ocean_OBC_type, OBC_NONE, OBC_SIMPLE
use MOM_open_boundary,  only : open_boundary_query
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_time_manager,   only : time_type, set_time, time_type_to_real

implicit none ; private

#include <MOM_memory.h>

public supercritical_initialize_topography
public supercritical_set_OBC_data

contains

! -----------------------------------------------------------------------------
!> This subroutine sets up the supercritical topography
subroutine supercritical_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),             intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),              intent(in)  :: param_file !< Parameter file structure
  real,                               intent(in)  :: max_depth  !< Maximum depth of model in m

  real :: min_depth ! The minimum and maximum depths in m.
  real :: PI
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "supercritical_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  PI = 4.0*atan(1.0) ;

  call MOM_mesg("  supercritical_initialization.F90, supercritical_initialize_topography: setting topography", 5)

  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)

  do j=js,je ; do i=is,ie
    D(i,j)=max_depth
    if ((G%geoLonT(i,j) > 10.0).AND. &
            (atan2(G%geoLatT(i,j),G%geoLonT(i,j)-10.0) < 8.95*PI/180.)) then
      D(i,j)=0.5*min_depth
    endif

    if (D(i,j) > max_depth) D(i,j) = max_depth
    if (D(i,j) < min_depth) D(i,j) = 0.5*min_depth
  enddo ; enddo

end subroutine supercritical_initialize_topography
! -----------------------------------------------------------------------------
!> This subroutine sets the properties of flow at open boundary conditions.
subroutine supercritical_set_OBC_data(OBC, G)
  type(ocean_OBC_type),   pointer    :: OBC  !< This open boundary condition type specifies
                                             !! whether, where, and what open boundary
                                             !! conditions are used.
  type(ocean_grid_type),  intent(in) :: G    !< The ocean's grid structure.

  ! The following variables are used to set up the transport in the TIDAL_BAY example.
  character(len=40)  :: mod = "supercritical_set_OBC_data" ! This subroutine's name.
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) return

  do j=jsd,jed ; do I=IsdB,IedB
    if (OBC%OBC_mask_u(I,j)) then
      OBC%eta_outer_u(I,j) = 0.0
      OBC%ubt_outer(I,j) = 8.57
    endif
  enddo ; enddo
  do J=JsdB,JedB ; do i=isd,ied
    if (OBC%OBC_mask_v(i,J)) then
      OBC%eta_outer_v(i,J) = 0.0
      OBC%vbt_outer(i,J) = 0.0
    endif
  enddo ; enddo

end subroutine supercritical_set_OBC_data

!> \class supercritical_initialization
!!
!! The module configures the model for the "supercritical" experiment.
!! https://marine.rutgers.edu/po/index.php?model=test-problems&title=supercritical
end module supercritical_initialization
