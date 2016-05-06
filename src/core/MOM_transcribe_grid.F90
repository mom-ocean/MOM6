module MOM_transcribe_grid

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

use MOM_hor_index, only : hor_index_type, hor_index_init
use MOM_domains, only : MOM_domain_type, get_domain_extent, compute_block_extent
use MOM_domains, only : pass_var, pass_vector, sum_across_PEs, broadcast
use MOM_domains, only : root_PE, To_All, SCALAR_PAIR, CGRID_NE, AGRID, BGRID_NE, CORNER
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type, set_derived_metrics
! use MOM_grid, only :  ! , set_first_direction

implicit none ; private

#include <MOM_memory.h>

contains

!> Copies information from a dynamic (shared) horizontal grid type into an
!! ocean_grid_type.
subroutine copy_dyngrid_to_MOM_grid(dG, G)
  type(dyn_horgrid_type), intent(in)    :: dG  !< Common horizontal grid type
  type(ocean_grid_type),  intent(inout) :: G   !< Ocean grid type

  integer :: isd, ied, jsd, jed, nk  ! Common data domains.
  integer :: IsdB, IedB, JsdB, JedB  ! Common data domains.
  integer :: ido, jdo  ! Indexing offsets between the grids.
  integer :: i, j, k

  ! MOM_grid_init and create_dyn_horgrid are called outside of this routine.
  ! This routine copies over the fields that were set by MOM_initilized_fixed.

  ! Determine the indexing offsets between the grids.
  ido = dG%idg_offset - G%idg_offset
  jdo = dG%jdg_offset - G%jdg_offset

  isd = max(G%isd, dG%isd+ido) ; jsd = max(G%jsd, dG%jsd+jdo)
  ied = min(G%ied, dG%ied+ido) ; jed = min(G%jed, dG%jed+jdo)
  IsdB = max(G%IsdB, dG%IsdB+ido) ; JsdB = max(G%JsdB, dG%JsdB+jdo)
  IedB = min(G%IedB, dG%IedB+ido) ; JedB = min(G%JedB, dG%JedB+jdo)

  ! Check that the grids conform.
  if ((isd > G%isc) .or. (ied < G%ied) .or. (jsd > G%jsc) .or. (jed > G%jed)) &
    call MOM_error(FATAL, "copy_dyngrid_to_MOM_grid called with incompatible grids.")

  do i=isd,ied ; do j=jsd,jed
    G%geoLonT(i,j) = dG%geoLonT(i+ido,j+jdo)
    G%geoLatT(i,j) = dG%geoLatT(i+ido,j+jdo)
    G%dxT(i,j) = dG%dxT(i+ido,j+jdo)
    G%dyT(i,j) = dG%dyT(i+ido,j+jdo)
    G%areaT(i,j) = dG%areaT(i+ido,j+jdo)
    G%bathyT(i,j) = dG%bathyT(i+ido,j+jdo)
  enddo ; enddo

  do I=IsdB,IedB ; do j=jsd,jed
    G%geoLonCu(I,j) = dG%geoLonCu(I+ido,j+jdo)
    G%geoLatCu(I,j) = dG%geoLatCu(I+ido,j+jdo)
    G%dxCu(I,j) = dG%dxCu(I+ido,j+jdo)
    G%dyCu(I,j) = dG%dyCu(I+ido,j+jdo)
  enddo ; enddo

  do i=isd,ied ; do J=JsdB,JedB
    G%geoLonCv(i,J) = dG%geoLonCv(i+ido,J+jdo)
    G%geoLatCv(i,J) = dG%geoLatCv(i+ido,J+jdo)
    G%dxCv(i,J) = dG%dxCv(i+ido,J+jdo)
    G%dyCv(i,J) = dG%dyCv(i+ido,J+jdo)
  enddo ; enddo

  do I=IsdB,IedB ; do J=JsdB,JedB
    G%geoLonBu(I,J) = dG%geoLonBu(I+ido,J+jdo)
    G%geoLatBu(I,J) = dG%geoLatBu(I+ido,J+jdo)
    G%dxBu(I,J) = dG%dxBu(I+ido,J+jdo)
    G%dyBu(I,J) = dG%dyBu(I+ido,J+jdo)
    G%areaBu(i,j) = dG%areaBu(I+ido,J+jdo)
  enddo ; enddo

! Update the halos in case the dynamic grid has smaller halos than the ocean grid.
  call pass_vector(G%dyCu, G%dxCv, G%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(G%dxCu, G%dyCv, G%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(G%dxBu, G%dyBu, G%Domain, To_All+Scalar_Pair, BGRID_NE)
  call pass_var(G%areaT, G%Domain)
  call pass_var(G%areaBu, G%Domain, position=CORNER)

  call set_derived_metrics(G)

end subroutine copy_dyngrid_to_MOM_grid

! subroutine set_first_direction(G, y_first)
!   type(ocean_grid_type), intent(inout) :: G
!   integer,               intent(in) :: y_first

!   G%first_direction = y_first
! end subroutine set_first_direction

!---------------------------------------------------------------------
!> Allocate memory used by the ocean_grid_type and related structures.
subroutine allocate_metrics(G)
  type(ocean_grid_type), intent(inout) :: G !< The horizontal grid type
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, isg, ieg, jsg, jeg

  ! This subroutine allocates the lateral elements of the ocean_grid_type that
  ! are always used and zeros them out.

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg

  ALLOC_(G%dxT(isd:ied,jsd:jed))       ; G%dxT(:,:) = 0.0
  ALLOC_(G%dxCu(IsdB:IedB,jsd:jed))    ; G%dxCu(:,:) = 0.0
  ALLOC_(G%dxCv(isd:ied,JsdB:JedB))    ; G%dxCv(:,:) = 0.0
  ALLOC_(G%dxBu(IsdB:IedB,JsdB:JedB))  ; G%dxBu(:,:) = 0.0
  ALLOC_(G%IdxT(isd:ied,jsd:jed))      ; G%IdxT(:,:) = 0.0
  ALLOC_(G%IdxCu(IsdB:IedB,jsd:jed))   ; G%IdxCu(:,:) = 0.0
  ALLOC_(G%IdxCv(isd:ied,JsdB:JedB))   ; G%IdxCv(:,:) = 0.0
  ALLOC_(G%IdxBu(IsdB:IedB,JsdB:JedB)) ; G%IdxBu(:,:) = 0.0

  ALLOC_(G%dyT(isd:ied,jsd:jed))       ; G%dyT(:,:) = 0.0
  ALLOC_(G%dyCu(IsdB:IedB,jsd:jed))    ; G%dyCu(:,:) = 0.0
  ALLOC_(G%dyCv(isd:ied,JsdB:JedB))    ; G%dyCv(:,:) = 0.0
  ALLOC_(G%dyBu(IsdB:IedB,JsdB:JedB))  ; G%dyBu(:,:) = 0.0
  ALLOC_(G%IdyT(isd:ied,jsd:jed))      ; G%IdyT(:,:) = 0.0
  ALLOC_(G%IdyCu(IsdB:IedB,jsd:jed))   ; G%IdyCu(:,:) = 0.0
  ALLOC_(G%IdyCv(isd:ied,JsdB:JedB))   ; G%IdyCv(:,:) = 0.0
  ALLOC_(G%IdyBu(IsdB:IedB,JsdB:JedB)) ; G%IdyBu(:,:) = 0.0

  ALLOC_(G%areaT(isd:ied,jsd:jed))       ; G%areaT(:,:) = 0.0
  ALLOC_(G%IareaT(isd:ied,jsd:jed))      ; G%IareaT(:,:) = 0.0
  ALLOC_(G%areaBu(IsdB:IedB,JsdB:JedB))  ; G%areaBu(:,:) = 0.0
  ALLOC_(G%IareaBu(IsdB:IedB,JsdB:JedB)) ; G%IareaBu(:,:) = 0.0

  ALLOC_(G%mask2dT(isd:ied,jsd:jed))      ; G%mask2dT(:,:) = 0.0
  ALLOC_(G%mask2dCu(IsdB:IedB,jsd:jed))   ; G%mask2dCu(:,:) = 0.0
  ALLOC_(G%mask2dCv(isd:ied,JsdB:JedB))   ; G%mask2dCv(:,:) = 0.0
  ALLOC_(G%mask2dBu(IsdB:IedB,JsdB:JedB)) ; G%mask2dBu(:,:) = 0.0
  ALLOC_(G%geoLatT(isd:ied,jsd:jed))      ; G%geoLatT(:,:) = 0.0
  ALLOC_(G%geoLatCu(IsdB:IedB,jsd:jed))   ; G%geoLatCu(:,:) = 0.0
  ALLOC_(G%geoLatCv(isd:ied,JsdB:JedB))   ; G%geoLatCv(:,:) = 0.0
  ALLOC_(G%geoLatBu(IsdB:IedB,JsdB:JedB)) ; G%geoLatBu(:,:) = 0.0
  ALLOC_(G%geoLonT(isd:ied,jsd:jed))      ; G%geoLonT(:,:) = 0.0
  ALLOC_(G%geoLonCu(IsdB:IedB,jsd:jed))   ; G%geoLonCu(:,:) = 0.0
  ALLOC_(G%geoLonCv(isd:ied,JsdB:JedB))   ; G%geoLonCv(:,:) = 0.0
  ALLOC_(G%geoLonBu(IsdB:IedB,JsdB:JedB)) ; G%geoLonBu(:,:) = 0.0

  ALLOC_(G%dx_Cv(isd:ied,JsdB:JedB))     ; G%dx_Cv(:,:) = 0.0
  ALLOC_(G%dy_Cu(IsdB:IedB,jsd:jed))     ; G%dy_Cu(:,:) = 0.0
  ALLOC_(G%dx_Cv_obc(isd:ied,JsdB:JedB)) ; G%dx_Cv_obc(:,:) = 0.0
  ALLOC_(G%dy_Cu_obc(IsdB:IedB,jsd:jed)) ; G%dy_Cu_obc(:,:) = 0.0

  ALLOC_(G%areaCu(IsdB:IedB,jsd:jed))  ; G%areaCu(:,:) = 0.0
  ALLOC_(G%areaCv(isd:ied,JsdB:JedB))  ; G%areaCv(:,:) = 0.0
  ALLOC_(G%IareaCu(IsdB:IedB,jsd:jed)) ; G%IareaCu(:,:) = 0.0
  ALLOC_(G%IareaCv(isd:ied,JsdB:JedB)) ; G%IareaCv(:,:) = 0.0

  ALLOC_(G%bathyT(isd:ied, jsd:jed)) ; G%bathyT(:,:) = 0.0
  ALLOC_(G%CoriolisBu(IsdB:IedB, JsdB:JedB)) ; G%CoriolisBu(:,:) = 0.0
  ALLOC_(G%dF_dx(isd:ied, jsd:jed)) ; G%dF_dx(:,:) = 0.0
  ALLOC_(G%dF_dy(isd:ied, jsd:jed)) ; G%dF_dy(:,:) = 0.0

  ALLOC_(G%sin_rot(isd:ied,jsd:jed)) ; G%sin_rot(:,:) = 0.0
  ALLOC_(G%cos_rot(isd:ied,jsd:jed)) ; G%cos_rot(:,:) = 1.0

  allocate(G%gridLonT(isg:ieg))   ; G%gridLonT(:) = 0.0
  allocate(G%gridLonB(G%IsgB:G%IegB)) ; G%gridLonB(:) = 0.0
  allocate(G%gridLatT(jsg:jeg))   ; G%gridLatT(:) = 0.0
  allocate(G%gridLatB(G%JsgB:G%JegB)) ; G%gridLatB(:) = 0.0

end subroutine allocate_metrics


end module MOM_transcribe_grid
