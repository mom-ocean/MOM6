!> The module calculates interface heights, including free surface height.
module MOM_interface_heights

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
!*  By Robert Hallberg, February 2001                                  *
!*                                                                     *
!*    The subroutines here calculate the heights of interfaces in a    *
!*  way that is consistent with the calculations in the pressure       *
!*  gradient acceleration calculation.  In a Boussinseq model this is  *
!*  pretty simple vertical sum, but in a non-Boussinesq model it uses  *
!*  integrals of the equation of state.                                *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, CoriolisBu                            *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, bathyT                                *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser, only : log_version
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : int_specific_vol_dp

implicit none ; private

#include <MOM_memory.h>

public find_eta

interface find_eta
  module procedure find_eta_2d, find_eta_3d
end interface find_eta

contains

subroutine find_eta_3d(h, tv, G_Earth, G, GV, eta, eta_bt, halo_size)
  type(ocean_grid_type),                      intent(in)  :: G
  type(verticalGrid_type),                    intent(in)  :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)  :: h
  type(thermo_var_ptrs),                      intent(in)  :: tv
  real,                                       intent(in)  :: G_Earth
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(out) :: eta
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)  :: eta_bt
  integer,                          optional, intent(in)  :: halo_size

!   This subroutine determines the heights of all interfaces between layers,
! using the appropriate form for consistency with the calculation of the
! pressure gradient forces.  Additionally, these height may be dilated for
! consistency with the corresponding time-average quantity from the barotropic
! calculation.

! Arguments:
!  (in)      h       - layer thickness H (meter or kg/m2)
!  (in)      tv      - structure pointing to thermodynamic variables
!  (in)      G_Earth - Earth gravitational acceleration (m/s2)
!  (in)      G       - ocean grid structure
!  (in)      GV      - The ocean's vertical grid structure.
!  (out)     eta     - layer interface heights (meter)
!  (in,opt)  eta_bt  - optional barotropic variable that gives the "correct"
!                      free surface height (Boussinesq) or total water column
!                      mass per unit aread (non-Boussinesq).  This is used to
!                      dilate the layer thicknesses when calculating interface
!                      heights, in H (m or kg m-2).
! (in,opt) halo_size - width of halo points on which to calculate eta

  real :: p(SZI_(G),SZJ_(G),SZK_(G)+1)
  real :: dz_geo(SZI_(G),SZJ_(G),SZK_(G)) ! The change in geopotential height
                                          ! across a layer, in m2 s-2.
  real :: dilate(SZI_(G))                 ! non-dimensional dilation factor
  real :: htot(SZI_(G))                   ! total thickness H
  real :: I_gEarth
  integer i, j, k, isv, iev, jsv, jev, nz, halo

  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)

  isv = G%isc-halo ; iev = G%iec+halo ; jsv = G%jsc-halo ; jev = G%jec+halo
  nz = G%ke

  if ((isv<G%isd) .or. (iev>G%ied) .or. (jsv<G%jsd) .or. (jev>G%jed)) &
    call MOM_error(FATAL,"find_eta called with an overly large halo_size.")

  I_gEarth = 1.0 / G_Earth

!$OMP parallel default(none) shared(isv,iev,jsv,jev,nz,eta,G,GV,h,eta_bt,tv,p, &
!$OMP                               G_Earth,dz_geo,halo,I_gEarth) &
!$OMP                       private(dilate,htot)
!$OMP do
  do j=jsv,jev ; do i=isv,iev ; eta(i,j,nz+1) = -G%bathyT(i,j) ; enddo ; enddo

  if (GV%Boussinesq) then
!$OMP do
    do j=jsv,jev ; do k=nz,1,-1; do i=isv,iev
      eta(i,j,K) = eta(i,j,K+1) + h(i,j,k)*GV%H_to_m
    enddo ; enddo ; enddo
    if (present(eta_bt)) then
      ! Dilate the water column to agree with the free surface height
      ! that is used for the dynamics.
!$OMP do
      do j=jsv,jev
        do i=isv,iev
          dilate(i) = (eta_bt(i,j)*GV%H_to_m + G%bathyT(i,j)) / &
                      (eta(i,j,1) + G%bathyT(i,j))
        enddo
        do k=1,nz ; do i=isv,iev
          eta(i,j,K) = dilate(i) * (eta(i,j,K) + G%bathyT(i,j)) - G%bathyT(i,j)
        enddo ; enddo
      enddo
    endif
  else
    if (associated(tv%eqn_of_state)) then
      ! ### THIS SHOULD BE P_SURF, IF AVAILABLE.
!$OMP do
      do j=jsv,jev
        do i=isv,iev ; p(i,j,1) = 0.0 ; enddo
        do k=1,nz ; do i=isv,iev
          p(i,j,K+1) = p(i,j,K) + G_Earth*GV%H_to_kg_m2*h(i,j,k)
        enddo ; enddo
      enddo
!$OMP do
      do k=1,nz
        call int_specific_vol_dp(tv%T(:,:,k), tv%S(:,:,k), p(:,:,K), p(:,:,K+1), &
                                 0.0, G%HI, tv%eqn_of_state, dz_geo(:,:,k), halo_size=halo)
      enddo
!$OMP do
      do j=jsv,jev
        do k=nz,1,-1 ; do i=isv,iev
          eta(i,j,K) = eta(i,j,K+1) + I_gEarth * dz_geo(i,j,k)
        enddo ; enddo
      enddo
    else
!$OMP do
      do j=jsv,jev ;  do k=nz,1,-1; do i=isv,iev
        eta(i,j,K) = eta(i,j,K+1) + GV%H_to_kg_m2*h(i,j,k)/GV%Rlay(k)
      enddo ; enddo ; enddo
    endif
    if (present(eta_bt)) then
      ! Dilate the water column to agree with the free surface height
      ! from the time-averaged barotropic solution.
!$OMP do
      do j=jsv,jev
        do i=isv,iev ; htot(i) = GV%H_subroundoff ; enddo
        do k=1,nz ; do i=isv,iev ; htot(i) = htot(i) + h(i,j,k) ; enddo ; enddo
        do i=isv,iev ; dilate(i) = eta_bt(i,j) / htot(i) ; enddo
        do k=1,nz ; do i=isv,iev
          eta(i,j,K) = dilate(i) * (eta(i,j,K) + G%bathyT(i,j)) - G%bathyT(i,j)
        enddo ; enddo
      enddo
    endif
  endif
!$OMP end parallel

end subroutine find_eta_3d


subroutine find_eta_2d(h, tv, G_Earth, G, GV, eta, eta_bt, halo_size)
  type(ocean_grid_type),                      intent(in)  :: G
  type(verticalGrid_type),                    intent(in)  :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)  :: h
  type(thermo_var_ptrs),                      intent(in)  :: tv
  real,                                       intent(in)  :: G_Earth
  real, dimension(SZI_(G),SZJ_(G)),           intent(out) :: eta
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)  :: eta_bt
  integer,                          optional, intent(in)  :: halo_size

!   This subroutine determines the free surface height, using the appropriate
! form for consistency with the calculation of the pressure gradient forces.
! Additionally, the sea surface height may be adjusted for consistency with the
! corresponding time-average quantity from the barotropic calculation.

! Arguments:
!  (in)       h      - layer thickness (meter or kg/m2)
!  (in)      tv      - structure pointing to various thermodynamic variables
!  (in)      G_Earth - Earth gravitational acceleration (m/s2)
!  (in)      G       - ocean grid structure
!  (in)      GV      - The ocean's vertical grid structure.
!  (out)     eta     - free surface height relative to mean sea level (z=0) (m)
!  (in,opt)  eta_bt  - optional barotropic variable that gives the "correct"
!                      free surface height (Boussinesq) or total water column
!                      mass per unit aread (non-Boussinesq), in H (m or kg m-2)
!  (in,opt)  halo_size - width of halo points on which to calculate eta

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: &
    p     ! The pressure in Pa.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    dz_geo     ! The change in geopotential height across a layer, in m2 s-2.
  real :: htot(SZI_(G))  ! The sum of all layers' thicknesses, in kg m-2 or m.
  real :: I_gEarth
  integer i, j, k, is, ie, js, je, nz, halo

  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)
  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  nz = G%ke

  I_gEarth = 1.0 / G_Earth

!$OMP parallel default(none) shared(is,ie,js,je,nz,eta,G,GV,eta_bt,h,tv,p, &
!$OMP                               G_Earth,dz_geo,halo,I_gEarth) &
!$OMP                       private(htot)
!$OMP do
  do j=js,je ; do i=is,ie ; eta(i,j) = -G%bathyT(i,j) ; enddo ; enddo

  if (GV%Boussinesq) then
    if (present(eta_bt)) then
!$OMP do
      do j=js,je ; do i=is,ie
        eta(i,j) = eta_bt(i,j)
      enddo ; enddo
    else
!$OMP do
      do j=js,je ; do k=1,nz ; do i=is,ie
        eta(i,j) = eta(i,j) + h(i,j,k)*GV%H_to_m
      enddo ; enddo ; enddo
    endif
  else
    if (associated(tv%eqn_of_state)) then
!$OMP do
      do j=js,je
        do i=is,ie ; p(i,j,1) = 0.0 ; enddo

        do k=1,nz ; do i=is,ie
          p(i,j,k+1) = p(i,j,k) + G_Earth*GV%H_to_kg_m2*h(i,j,k)
        enddo ; enddo
      enddo
!$OMP do
      do k = 1, nz
        call int_specific_vol_dp(tv%T(:,:,k), tv%S(:,:,k), p(:,:,k), p(:,:,k+1), 0.0, &
                                 G%HI, tv%eqn_of_state, dz_geo(:,:,k), halo_size=halo)
      enddo
!$OMP do
      do j=js,je ; do k=1,nz ; do i=is,ie
          eta(i,j) = eta(i,j) + I_gEarth * dz_geo(i,j,k)
      enddo ; enddo ; enddo
    else
!$OMP do
      do j=js,je ; do k=1,nz ; do i=is,ie
        eta(i,j) = eta(i,j) + GV%H_to_kg_m2*h(i,j,k)/GV%Rlay(k)
      enddo ; enddo ; enddo
    endif
    if (present(eta_bt)) then
      !   Dilate the water column to agree with the the time-averaged column
      ! mass from the barotropic solution.
!$OMP do
      do j=js,je
        do i=is,ie ; htot(i) = GV%H_subroundoff ; enddo
        do k=1,nz ; do i=is,ie ; htot(i) = htot(i) + h(i,j,k) ; enddo ; enddo
        do i=is,ie
          eta(i,j) = (eta_bt(i,j) / htot(i)) * (eta(i,j) + G%bathyT(i,j)) - &
                     G%bathyT(i,j)
        enddo
      enddo
    endif
  endif
!$OMP end parallel

end subroutine find_eta_2d

end module MOM_interface_heights
