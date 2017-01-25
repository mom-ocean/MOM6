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

use MOM_domains, only : pass_var, pass_vector
use MOM_domains, only : To_All, SCALAR_PAIR, CGRID_NE, AGRID, BGRID_NE, CORNER
use MOM_dyn_horgrid, only : dyn_horgrid_type, set_derived_dyn_horgrid
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_grid, only : ocean_grid_type, set_derived_metrics

implicit none ; private

public copy_dyngrid_to_MOM_grid, copy_MOM_grid_to_dyngrid

contains

!> Copies information from a dynamic (shared) horizontal grid type into an
!! ocean_grid_type.
subroutine copy_dyngrid_to_MOM_grid(dG, oG)
  type(dyn_horgrid_type), intent(in)    :: dG  !< Common horizontal grid type
  type(ocean_grid_type),  intent(inout) :: oG  !< Ocean grid type

  integer :: isd, ied, jsd, jed      ! Common data domains.
  integer :: IsdB, IedB, JsdB, JedB  ! Common data domains.
  integer :: ido, jdo, Ido2, Jdo2    ! Indexing offsets between the grids.
  integer :: Igst, Jgst              ! Global starting indices.
  integer :: i, j

  ! MOM_grid_init and create_dyn_horgrid are called outside of this routine.
  ! This routine copies over the fields that were set by MOM_initialized_fixed.

  ! Determine the indexing offsets between the grids.
  ido = dG%idg_offset - oG%idg_offset
  jdo = dG%jdg_offset - oG%jdg_offset

  isd = max(oG%isd, dG%isd+ido) ; jsd = max(oG%jsd, dG%jsd+jdo)
  ied = min(oG%ied, dG%ied+ido) ; jed = min(oG%jed, dG%jed+jdo)
  IsdB = max(oG%IsdB, dG%IsdB+ido) ; JsdB = max(oG%JsdB, dG%JsdB+jdo)
  IedB = min(oG%IedB, dG%IedB+ido) ; JedB = min(oG%JedB, dG%JedB+jdo)

  ! Check that the grids conform.
  if ((isd > oG%isc) .or. (ied < oG%ied) .or. (jsd > oG%jsc) .or. (jed > oG%jed)) &
    call MOM_error(FATAL, "copy_dyngrid_to_MOM_grid called with incompatible grids.")

  do i=isd,ied ; do j=jsd,jed
    oG%geoLonT(i,j) = dG%geoLonT(i+ido,j+jdo)
    oG%geoLatT(i,j) = dG%geoLatT(i+ido,j+jdo)
    oG%dxT(i,j) = dG%dxT(i+ido,j+jdo)
    oG%dyT(i,j) = dG%dyT(i+ido,j+jdo)
    oG%areaT(i,j) = dG%areaT(i+ido,j+jdo)
    oG%bathyT(i,j) = dG%bathyT(i+ido,j+jdo)

    oG%dF_dx(i,j) = dG%dF_dx(i+ido,j+jdo)
    oG%dF_dy(i,j) = dG%dF_dy(i+ido,j+jdo)
    oG%sin_rot(i,j) = dG%sin_rot(i+ido,j+jdo)
    oG%cos_rot(i,j) = dG%cos_rot(i+ido,j+jdo)
    oG%mask2dT(i,j) = dG%mask2dT(i+ido,j+jdo)
  enddo ; enddo

  do I=IsdB,IedB ; do j=jsd,jed
    oG%geoLonCu(I,j) = dG%geoLonCu(I+ido,j+jdo)
    oG%geoLatCu(I,j) = dG%geoLatCu(I+ido,j+jdo)
    oG%dxCu(I,j) = dG%dxCu(I+ido,j+jdo)
    oG%dyCu(I,j) = dG%dyCu(I+ido,j+jdo)
    oG%dy_Cu(I,j) = dG%dy_Cu(I+ido,j+jdo)
    oG%dy_Cu_obc(I,j) = dG%dy_Cu_obc(I+ido,j+jdo)

    oG%mask2dCu(I,j) = dG%mask2dCu(I+ido,j+jdo)
    oG%areaCu(I,j) = dG%areaCu(I+ido,j+jdo)
    oG%IareaCu(I,j) = dG%IareaCu(I+ido,j+jdo)
  enddo ; enddo

  do i=isd,ied ; do J=JsdB,JedB
    oG%geoLonCv(i,J) = dG%geoLonCv(i+ido,J+jdo)
    oG%geoLatCv(i,J) = dG%geoLatCv(i+ido,J+jdo)
    oG%dxCv(i,J) = dG%dxCv(i+ido,J+jdo)
    oG%dyCv(i,J) = dG%dyCv(i+ido,J+jdo)
    oG%dx_Cv(i,J) = dG%dx_Cv(i+ido,J+jdo)
    oG%dx_Cv_obc(i,J) = dG%dx_Cv_obc(i+ido,J+jdo)

    oG%mask2dCv(i,J) = dG%mask2dCv(i+ido,J+jdo)
    oG%areaCv(i,J) = dG%areaCv(i+ido,J+jdo)
    oG%IareaCv(i,J) = dG%IareaCv(i+ido,J+jdo)
  enddo ; enddo

  do I=IsdB,IedB ; do J=JsdB,JedB
    oG%geoLonBu(I,J) = dG%geoLonBu(I+ido,J+jdo)
    oG%geoLatBu(I,J) = dG%geoLatBu(I+ido,J+jdo)
    oG%dxBu(I,J) = dG%dxBu(I+ido,J+jdo)
    oG%dyBu(I,J) = dG%dyBu(I+ido,J+jdo)
    oG%areaBu(I,J) = dG%areaBu(I+ido,J+jdo)
    oG%CoriolisBu(I,J) = dG%CoriolisBu(I+ido,J+jdo)
    oG%mask2dBu(I,J) = dG%mask2dBu(I+ido,J+jdo)
  enddo ; enddo

  oG%bathymetry_at_vel = dG%bathymetry_at_vel
  if (oG%bathymetry_at_vel) then
    do I=IsdB,IedB ; do j=jsd,jed
      oG%Dblock_u(I,j) = dG%Dblock_u(I+ido,j+jdo)
      oG%Dopen_u(I,j) = dG%Dopen_u(I+ido,j+jdo)
    enddo ; enddo
    do i=isd,ied ; do J=JsdB,JedB
      oG%Dblock_v(i,J) = dG%Dblock_v(i+ido,J+jdo)
      oG%Dopen_v(i,J) = dG%Dopen_v(i+ido,J+jdo)
    enddo ; enddo
  endif

  oG%gridLonT(oG%isg:oG%ieg) = dG%gridLonT(dG%isg:dG%ieg)
  oG%gridLatT(oG%jsg:oG%jeg) = dG%gridLatT(dG%jsg:dG%jeg)
  ! The more complicated logic here avoids segmentation faults if one grid uses
  ! global symmetric memory while the other does not.  Because a northeast grid
  ! convention is being used, the upper bounds for each array correspond.
  !   Note that the dynamic grid always uses symmetric memory.
  Ido2 = dG%IegB-oG%IegB ; Igst = max(oG%IsgB, (dG%isg-1)-Ido2)
  Jdo2 = dG%JegB-oG%JegB ; Jgst = max(oG%JsgB, (dG%jsg-1)-Jdo2)
  do I=Igst,oG%IegB ; oG%gridLonB(I) = dG%gridLonB(I+Ido2) ; enddo
  do J=Jgst,oG%JegB ; oG%gridLatB(J) = dG%gridLatB(J+Jdo2) ; enddo

  ! Copy various scalar variables and strings.
  oG%x_axis_units = dG%x_axis_units ; oG%y_axis_units = dG%y_axis_units
  oG%areaT_global = dG%areaT_global ; oG%IareaT_global = dG%IareaT_global
  oG%south_lat = dG%south_lat ; oG%west_lon  = dG%west_lon
  oG%len_lat = dG%len_lat ; oG%len_lon = dG%len_lon
  oG%Rad_Earth = dG%Rad_Earth ; oG%max_depth = dG%max_depth

! Update the halos in case the dynamic grid has smaller halos than the ocean grid.
  call pass_var(oG%areaT, oG%Domain)
  call pass_var(oG%bathyT, oG%Domain)
  call pass_var(oG%geoLonT, oG%Domain)
  call pass_var(oG%geoLatT, oG%Domain)
  call pass_vector(oG%dxT, oG%dyT, oG%Domain, To_All+Scalar_Pair, AGRID)
  call pass_vector(oG%dF_dx, oG%dF_dy, oG%Domain, To_All, AGRID)
  call pass_vector(oG%cos_rot, oG%sin_rot, oG%Domain, To_All, AGRID)
  call pass_var(oG%mask2dT, oG%Domain)

  call pass_vector(oG%areaCu, oG%areaCv, oG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(oG%dyCu, oG%dxCv, oG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(oG%dxCu, oG%dyCv, oG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(oG%dy_Cu, oG%dx_Cv, oG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(oG%dy_Cu_obc, oG%dx_Cv_obc, oG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(oG%mask2dCu, oG%mask2dCv, oG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(oG%IareaCu, oG%IareaCv, oG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(oG%IareaCu, oG%IareaCv, oG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(oG%geoLatCu, oG%geoLatCv, oG%Domain, To_All+Scalar_Pair, CGRID_NE)

  call pass_var(oG%areaBu, oG%Domain, position=CORNER)
  call pass_var(oG%geoLonBu, oG%Domain, position=CORNER)
  call pass_var(oG%geoLatBu, oG%Domain, position=CORNER)
  call pass_vector(oG%dxBu, oG%dyBu, oG%Domain, To_All+Scalar_Pair, BGRID_NE)
  call pass_var(oG%CoriolisBu, oG%Domain, position=CORNER)
  call pass_var(oG%mask2dBu, oG%Domain, position=CORNER)

  if (oG%bathymetry_at_vel) then
    call pass_vector(oG%Dblock_u, oG%Dblock_v, oG%Domain, To_All+Scalar_Pair, CGRID_NE)
    call pass_vector(oG%Dopen_u, oG%Dopen_v, oG%Domain, To_All+Scalar_Pair, CGRID_NE)
  endif

  call set_derived_metrics(oG)

end subroutine copy_dyngrid_to_MOM_grid


!> Copies information from an ocean_grid_type into a dynamic (shared)
!! horizontal grid type.
subroutine copy_MOM_grid_to_dyngrid(oG, dG)
  type(ocean_grid_type),  intent(in)    :: oG  !< Ocean grid type
  type(dyn_horgrid_type), intent(inout) :: dG  !< Common horizontal grid type

  integer :: isd, ied, jsd, jed      ! Common data domains.
  integer :: IsdB, IedB, JsdB, JedB  ! Common data domains.
  integer :: ido, jdo, Ido2, Jdo2    ! Indexing offsets between the grids.
  integer :: Igst, Jgst              ! Global starting indices.
  integer :: i, j

  ! MOM_grid_init and create_dyn_horgrid are called outside of this routine.
  ! This routine copies over the fields that were set by MOM_initialized_fixed.

  ! Determine the indexing offsets between the grids.
  ido = oG%idG_offset - dG%idG_offset
  jdo = oG%jdG_offset - dG%jdG_offset

  isd = max(dG%isd, oG%isd+ido) ; jsd = max(dG%jsd, oG%jsd+jdo)
  ied = min(dG%ied, oG%ied+ido) ; jed = min(dG%jed, oG%jed+jdo)
  IsdB = max(dG%IsdB, oG%IsdB+ido) ; JsdB = max(dG%JsdB, oG%JsdB+jdo)
  IedB = min(dG%IedB, oG%IedB+ido) ; JedB = min(dG%JedB, oG%JedB+jdo)

  ! Check that the grids conform.
  if ((isd > dG%isc) .or. (ied < dG%ied) .or. (jsd > dG%jsc) .or. (jed > dG%jed)) &
    call MOM_error(FATAL, "copy_dyngrid_to_MOM_grid called with incompatible grids.")

  do i=isd,ied ; do j=jsd,jed
    dG%geoLonT(i,j) = oG%geoLonT(i+ido,j+jdo)
    dG%geoLatT(i,j) = oG%geoLatT(i+ido,j+jdo)
    dG%dxT(i,j) = oG%dxT(i+ido,j+jdo)
    dG%dyT(i,j) = oG%dyT(i+ido,j+jdo)
    dG%areaT(i,j) = oG%areaT(i+ido,j+jdo)
    dG%bathyT(i,j) = oG%bathyT(i+ido,j+jdo)

    dG%dF_dx(i,j) = oG%dF_dx(i+ido,j+jdo)
    dG%dF_dy(i,j) = oG%dF_dy(i+ido,j+jdo)
    dG%sin_rot(i,j) = oG%sin_rot(i+ido,j+jdo)
    dG%cos_rot(i,j) = oG%cos_rot(i+ido,j+jdo)
    dG%mask2dT(i,j) = oG%mask2dT(i+ido,j+jdo)
  enddo ; enddo

  do I=IsdB,IedB ; do j=jsd,jed
    dG%geoLonCu(I,j) = oG%geoLonCu(I+ido,j+jdo)
    dG%geoLatCu(I,j) = oG%geoLatCu(I+ido,j+jdo)
    dG%dxCu(I,j) = oG%dxCu(I+ido,j+jdo)
    dG%dyCu(I,j) = oG%dyCu(I+ido,j+jdo)
    dG%dy_Cu(I,j) = oG%dy_Cu(I+ido,j+jdo)
    dG%dy_Cu_obc(I,j) = oG%dy_Cu_obc(I+ido,j+jdo)

    dG%mask2dCu(I,j) = oG%mask2dCu(I+ido,j+jdo)
    dG%areaCu(I,j) = oG%areaCu(I+ido,j+jdo)
    dG%IareaCu(I,j) = oG%IareaCu(I+ido,j+jdo)
  enddo ; enddo

  do i=isd,ied ; do J=JsdB,JedB
    dG%geoLonCv(i,J) = oG%geoLonCv(i+ido,J+jdo)
    dG%geoLatCv(i,J) = oG%geoLatCv(i+ido,J+jdo)
    dG%dxCv(i,J) = oG%dxCv(i+ido,J+jdo)
    dG%dyCv(i,J) = oG%dyCv(i+ido,J+jdo)
    dG%dx_Cv(i,J) = oG%dx_Cv(i+ido,J+jdo)
    dG%dx_Cv_obc(i,J) = oG%dx_Cv_obc(i+ido,J+jdo)

    dG%mask2dCv(i,J) = oG%mask2dCv(i+ido,J+jdo)
    dG%areaCv(i,J) = oG%areaCv(i+ido,J+jdo)
    dG%IareaCv(i,J) = oG%IareaCv(i+ido,J+jdo)
  enddo ; enddo

  do I=IsdB,IedB ; do J=JsdB,JedB
    dG%geoLonBu(I,J) = oG%geoLonBu(I+ido,J+jdo)
    dG%geoLatBu(I,J) = oG%geoLatBu(I+ido,J+jdo)
    dG%dxBu(I,J) = oG%dxBu(I+ido,J+jdo)
    dG%dyBu(I,J) = oG%dyBu(I+ido,J+jdo)
    dG%areaBu(I,J) = oG%areaBu(I+ido,J+jdo)
    dG%CoriolisBu(I,J) = oG%CoriolisBu(I+ido,J+jdo)
    dG%mask2dBu(I,J) = oG%mask2dBu(I+ido,J+jdo)
  enddo ; enddo

  dG%bathymetry_at_vel = oG%bathymetry_at_vel
  if (dG%bathymetry_at_vel) then
    do I=IsdB,IedB ; do j=jsd,jed
      dG%Dblock_u(I,j) = oG%Dblock_u(I+ido,j+jdo)
      dG%Dopen_u(I,j) = oG%Dopen_u(I+ido,j+jdo)
    enddo ; enddo
    do i=isd,ied ; do J=JsdB,JedB
      dG%Dblock_v(i,J) = oG%Dblock_v(i+ido,J+jdo)
      dG%Dopen_v(i,J) = oG%Dopen_v(i+ido,J+jdo)
    enddo ; enddo
  endif

  dG%gridLonT(dG%isg:dG%ieg) = oG%gridLonT(oG%isg:oG%ieg)
  dG%gridLatT(dG%jsg:dG%jeg) = oG%gridLatT(oG%jsg:oG%jeg)

  ! The more complicated logic here avoids segmentation faults if one grid uses
  ! global symmetric memory while the other does not.  Because a northeast grid
  ! convention is being used, the upper bounds for each array correspond.
  !   Note that the dynamic grid always uses symmetric memory.
  Ido2 = oG%IegB-dG%IegB ; Igst = max(dG%isg-1, oG%IsgB-Ido2)
  Jdo2 = oG%JegB-dG%JegB ; Jgst = max(dG%jsg-1, oG%JsgB-Jdo2)
  do I=Igst,dG%IegB ; dG%gridLonB(I) = oG%gridLonB(I+Ido2) ; enddo
  do J=Jgst,dG%JegB ; dG%gridLatB(J) = oG%gridLatB(J+Jdo2) ; enddo

  ! Copy various scalar variables and strings.
  dG%x_axis_units = oG%x_axis_units ; dG%y_axis_units = oG%y_axis_units
  dG%areaT_global = oG%areaT_global ; dG%IareaT_global = oG%IareaT_global
  dG%south_lat = oG%south_lat ; dG%west_lon  = oG%west_lon
  dG%len_lat = oG%len_lat ; dG%len_lon = oG%len_lon
  dG%Rad_Earth = oG%Rad_Earth ; dG%max_depth = oG%max_depth

! Update the halos in case the dynamic grid has smaller halos than the ocean grid.
  call pass_var(dG%areaT, dG%Domain)
  call pass_var(dG%bathyT, dG%Domain)
  call pass_var(dG%geoLonT, dG%Domain)
  call pass_var(dG%geoLatT, dG%Domain)
  call pass_vector(dG%dxT, dG%dyT, dG%Domain, To_All+Scalar_Pair, AGRID)
  call pass_vector(dG%dF_dx, dG%dF_dy, dG%Domain, To_All, AGRID)
  call pass_vector(dG%cos_rot, dG%sin_rot, dG%Domain, To_All, AGRID)
  call pass_var(dG%mask2dT, dG%Domain)

  call pass_vector(dG%areaCu, dG%areaCv, dG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dG%dyCu, dG%dxCv, dG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dG%dxCu, dG%dyCv, dG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dG%dy_Cu, dG%dx_Cv, dG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dG%dy_Cu_obc, dG%dx_Cv_obc, dG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dG%mask2dCu, dG%mask2dCv, dG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dG%IareaCu, dG%IareaCv, dG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dG%IareaCu, dG%IareaCv, dG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dG%geoLatCu, dG%geoLatCv, dG%Domain, To_All+Scalar_Pair, CGRID_NE)

  call pass_var(dG%areaBu, dG%Domain, position=CORNER)
  call pass_var(dG%geoLonBu, dG%Domain, position=CORNER)
  call pass_var(dG%geoLatBu, dG%Domain, position=CORNER)
  call pass_vector(dG%dxBu, dG%dyBu, dG%Domain, To_All+Scalar_Pair, BGRID_NE)
  call pass_var(dG%CoriolisBu, dG%Domain, position=CORNER)
  call pass_var(dG%mask2dBu, dG%Domain, position=CORNER)

  if (dG%bathymetry_at_vel) then
    call pass_vector(dG%Dblock_u, dG%Dblock_v, dG%Domain, To_All+Scalar_Pair, CGRID_NE)
    call pass_vector(dG%Dopen_u, dG%Dopen_v, dG%Domain, To_All+Scalar_Pair, CGRID_NE)
  endif

  call  set_derived_dyn_horgrid(dG)

end subroutine copy_MOM_grid_to_dyngrid

end module MOM_transcribe_grid
