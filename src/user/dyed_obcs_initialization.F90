module dyed_obcs_initialization
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

use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE, OBC_SIMPLE
use MOM_open_boundary, only : OBC_segment_type
use MOM_tracer_registry, only : tracer_registry_type, add_tracer_OBC_values
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public dyed_obcs_set_OBC_data

contains

!> This subroutine sets the dye properties at open boundary conditions.
subroutine dyed_obcs_set_OBC_data(OBC, G, GV, param_file, tr_Reg)
  type(ocean_OBC_type),       pointer    :: OBC !< This open boundary condition type specifies
                                                !! whether, where, and what open boundary
                                                !! conditions are used.
  type(ocean_grid_type),      intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in) :: GV  !< The ocean's vertical grid structure.
  type(param_file_type),      intent(in) :: param_file !< A structure indicating the open file
                              !! to parse for model parameter values.
  type(tracer_registry_type), pointer    :: tr_Reg !< Tracer registry.

  real, pointer, dimension(:,:,:) :: &
    OBC_dye_1_v => NULL(), &    ! Specify the dye concentrations at the boundaries,
    OBC_dye_2_v => NULL(), &    ! at both u and v points.
    OBC_dye_3_u => NULL(), &
    OBC_dye_4_u => NULL()

  character(len=40)  :: mdl = "dyed_obcs_set_OBC_data" ! This subroutine's name.
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB
  type(OBC_segment_type), pointer :: segment

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) return

  if (OBC%number_of_segments .ne. 4) then
    print *, 'Error in dyed_obcs segment setup'
    return   !!! Need a better error message here
  endif

  allocate(OBC_dye_1_v(isd:ied,JsdB:JedB,nz)) ; OBC_dye_1_v(:,:,:) = 0.0
  allocate(OBC_dye_2_v(isd:ied,JsdB:JedB,nz)) ; OBC_dye_2_v(:,:,:) = 0.0
  allocate(OBC_dye_3_u(IsdB:IedB,jsd:jed,nz)) ; OBC_dye_3_u(:,:,:) = 0.0
  allocate(OBC_dye_4_u(IsdB:IedB,jsd:jed,nz)) ; OBC_dye_4_u(:,:,:) = 0.0

! do n=1,OBC%number_of_segments
!   segment => OBC%segment(n)
!   if (.not. segment%on_pe) return
!   ! New way (not yet)
!   isd = segment%HI%isd ; ied = segment%HI%ied
!   JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB

!   do k=1,nz
!     do J=JsdB,JedB ; do i=isd,ied
!       segment%
!       segment%
!     enddo ; enddo
!   enddo
! enddo

  ! Set the inflow values of the dyes, one per segment.
  ! We know the order: north, south, east, west
  segment => OBC%segment(1)
  isd = segment%HI%isd ; ied = segment%HI%ied
  JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB
  do k=1,nz
    do J=JsdB,JedB ; do i=isd,ied
      OBC_dye_1_v(i,J,k) = 1.0
    enddo ; enddo
  enddo
  segment => OBC%segment(2)
  isd = segment%HI%isd ; ied = segment%HI%ied
  JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB
  do k=1,nz
    do J=JsdB,JedB ; do i=isd,ied
      OBC_dye_2_v(i,J,k) = 1.0
    enddo ; enddo
  enddo
  segment => OBC%segment(3)
  IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
  jsd = segment%HI%jsd ; jed = segment%HI%jed
  do k=1,nz
    do j=jsd,jed ; do I=IsdB,IedB
      OBC_dye_3_u(I,j,k) = 1.0
    enddo ; enddo
  enddo
  segment => OBC%segment(4)
  IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
  jsd = segment%HI%jsd ; jed = segment%HI%jed
  do k=1,nz
    do j=jsd,jed ; do I=IsdB,IedB
      OBC_dye_4_u(I,j,k) = 1.0
    enddo ; enddo
  enddo

  call add_tracer_OBC_values("dye_1", tr_Reg, OBC_in_v=OBC_dye_1_v)
  call add_tracer_OBC_values("dye_2", tr_Reg, OBC_in_v=OBC_dye_2_v)
  call add_tracer_OBC_values("dye_3", tr_Reg, OBC_in_u=OBC_dye_3_u)
  call add_tracer_OBC_values("dye_4", tr_Reg, OBC_in_u=OBC_dye_4_u)

end subroutine dyed_obcs_set_OBC_data

!> \namespace dyed_obcs_initialization
!! Setting dyes, one for painting the inflow on each side.
end module dyed_obcs_initialization
