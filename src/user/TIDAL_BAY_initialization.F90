module TIDAL_BAY_initialization
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
use MOM_open_boundary,  only : open_boundary_query, set_Flather_positions
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_time_manager,   only : time_type, set_time, time_type_to_real

implicit none ; private

#include <MOM_memory.h>

public TIDAL_BAY_set_OBC_positions
public TIDAL_BAY_alloc_OBC_data
public TIDAL_BAY_set_OBC_data

contains

!> Set the positions of the open boundary needed for the TIDAL_BAY experiment.
subroutine TIDAL_BAY_set_OBC_positions(G, param_file, OBC)
  type(dyn_horgrid_type),     intent(inout) :: G   !< Grid structure.
  type(param_file_type),      intent(in)    :: param_file !< Parameter file handle.
  type(ocean_OBC_type),       pointer       :: OBC !< Open boundary control structure.
  ! Local variables
  character(len=40)  :: mod = "TIDAL_BAY_set_OBC_positions" ! This subroutine's name.
  integer :: i, j

  if (.not.associated(OBC)) call MOM_error(FATAL, &
           "TIDAL_BAY_initialization, TIDAL_BAY_set_OBC_positions: OBC type was not allocated!")

  ! This isn't called when APPLY_OBC_U is requested.
  if (open_boundary_query(OBC, apply_orig_Flather=.true.)) then
    call set_Flather_positions(G, OBC)
    call TIDAL_BAY_alloc_OBC_data(OBC, G)
  endif
  OBC%update_OBC = .true.
  if (OBC%apply_OBC_v) then
    ! Set where v points are determined by OBCs.
    !allocate(OBC_mask_v(IsdB:IedB,jsd:jed)) ; OBC_mask_v(:,:) = .false.
    call MOM_error(FATAL,"TIDAL_BAY_initialization, TIDAL_BAY_set_OBC_positions: "//&
                   "APPLY_OBC_V=True is not coded for the TIDAL_BAY experiment")
  endif

end subroutine TIDAL_BAY_set_OBC_positions

!> This subroutine allocates the arrays for open boundary conditions.
subroutine TIDAL_BAY_alloc_OBC_data(OBC, G)
  type(ocean_OBC_type),   pointer    :: OBC  !< This open boundary condition type specifies
                                             !! whether, where, and what open boundary
                                             !! conditions are used.
  type(dyn_horgrid_type),  intent(in) :: G    !< The ocean's grid structure.

  logical :: apply_OBC_u, apply_OBC_v
  character(len=40)  :: mod = "TIDAL_BAY_set_OBC_data" ! This subroutine's name.
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) return
  if (.not.(OBC%apply_OBC_u .or. OBC%apply_OBC_v)) return

  if (.not.associated(OBC%vbt_outer)) then
    allocate(OBC%vbt_outer(isd:ied,JsdB:JedB)) ; OBC%vbt_outer(:,:) = 0.0 
  endif

  if (.not.associated(OBC%ubt_outer)) then
    allocate(OBC%ubt_outer(IsdB:IedB,jsd:jed)) ; OBC%ubt_outer(:,:) = 0.0 
  endif

  if (.not.associated(OBC%eta_outer_u)) then
    allocate(OBC%eta_outer_u(IsdB:IedB,jsd:jed)) ; OBC%eta_outer_u(:,:) = 0.0 
  endif

  if (.not.associated(OBC%eta_outer_v)) then
    allocate(OBC%eta_outer_v(isd:ied,JsdB:JedB)) ; OBC%eta_outer_v(:,:) = 0.0
  endif

!  call pass_vector(OBC%eta_outer_u,OBC%eta_outer_v,G%Domain, To_All+SCALAR_PAIR, CGRID_NE)
!  call pass_vector(OBC%ubt_outer,OBC%vbt_outer,G%Domain)

end subroutine TIDAL_BAY_alloc_OBC_data

!> This subroutine sets the properties of flow at open boundary conditions.
subroutine TIDAL_BAY_set_OBC_data(OBC, G, Time)
  type(ocean_OBC_type),   pointer    :: OBC  !< This open boundary condition type specifies
                                             !! whether, where, and what open boundary
                                             !! conditions are used.
  type(ocean_grid_type),  intent(in) :: G    !< The ocean's grid structure.
  type(time_type),        intent(in) :: Time !< model time.

  logical :: apply_OBC_u, apply_OBC_v
  ! The following variables are used to set up the transport in the TIDAL_BAY example.
  real :: time_sec, cff, cff2, tide_flow
  real :: my_area, my_flux
  real, parameter :: pi = 3.1415926535
  character(len=40)  :: mod = "TIDAL_BAY_set_OBC_data" ! This subroutine's name.
  integer :: i, j, itt, is, ie, js, je, isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) return
  if (.not.(OBC%apply_OBC_u .or. OBC%apply_OBC_v)) return

  if (OBC%apply_OBC_u) then
    time_sec = time_type_to_real(Time)
    cff = 0.1*sin(2.0*pi*time_sec/(12.0*3600.0))
    tide_flow = 3.0e6
    my_area=0.0
    my_flux=0.0
    do J=JsdB,JedB ; do i=isd,ied
! HACK to fix
!            cff2 = 0.5*(zeta(Iend  ,j,knew)+h(Iend  ,j)+                &
!     &                  zeta(Iend+1,j,knew)+h(Iend+1,j))/pn(Iend,j)
!            my_area = my_area+cff2
      if (OBC%OBC_mask_u(I,j)) then
        cff2 = 35*2000.
        my_area = my_area+cff2
      endif
    enddo ; enddo
    my_flux = -tide_flow*SIN(2.0*pi*time_sec/(12.0*3600.0))

    do J=JsdB,JedB ; do i=isd,ied
      if (OBC%OBC_mask_u(I,j)) then
        OBC%eta_outer_u(I,j) = cff
        OBC%ubt_outer(I,j) = my_flux/my_area
      endif
      if (OBC%OBC_mask_v(i,J)) then
        OBC%eta_outer_v(i,J) = cff
        OBC%vbt_outer(i,J) = 0.0
      endif
    enddo ; enddo
  endif

end subroutine TIDAL_BAY_set_OBC_data

!> \class TIDAL_BAY_initialization
!!
!! The module configures the model for the "TIDAL_BAY" experiment.
!! TIDAL_BAY = Tidally resonant bay from Zygmunt Kowalik's class on tides.
end module TIDAL_BAY_initialization
