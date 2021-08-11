!> Functions for calculating interface heights, including free surface height.
module MOM_interface_heights

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser, only : log_version
use MOM_grid, only : ocean_grid_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_density_integrals, only : int_specific_vol_dp

implicit none ; private

#include <MOM_memory.h>

public find_eta

!> Calculates the heights of sruface or all interfaces from layer thicknesses.
interface find_eta
  module procedure find_eta_2d, find_eta_3d
end interface find_eta

contains

!> Calculates the heights of all interfaces between layers, using the appropriate
!! form for consistency with the calculation of the pressure gradient forces.
!! Additionally, these height may be dilated for consistency with the
!! corresponding time-average quantity from the barotropic calculation.
subroutine find_eta_3d(h, tv, G, GV, US, eta, eta_bt, halo_size, eta_to_m)
  type(ocean_grid_type),                      intent(in)  :: G   !< The ocean's grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< A structure pointing to various
                                                                 !! thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(out) :: eta !< layer interface heights
                                                                 !! [Z ~> m] or [1/eta_to_m m].
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)  :: eta_bt !< optional barotropic
             !! variable that gives the "correct" free surface height (Boussinesq) or total water
             !! column mass per unit area (non-Boussinesq).  This is used to dilate the layer.
             !! thicknesses when calculating interfaceheights [H ~> m or kg m-2].
  integer,                          optional, intent(in)  :: halo_size !< width of halo points on
                                                                 !! which to calculate eta.
  real,                             optional, intent(in)  :: eta_to_m  !< The conversion factor from
             !! the units of eta to m; by default this is US%Z_to_m.

  ! Local variables
  real :: p(SZI_(G),SZJ_(G),SZK_(GV)+1)   ! Hydrostatic pressure at each interface [R L2 T-2 ~> Pa]
  real :: dz_geo(SZI_(G),SZJ_(G),SZK_(GV)) ! The change in geopotential height
                                           ! across a layer [L2 T-2 ~> m2 s-2].
  real :: dilate(SZI_(G))                 ! non-dimensional dilation factor
  real :: htot(SZI_(G))                   ! total thickness [H ~> m or kg m-2]
  real :: I_gEarth          ! The inverse of the gravitational acceleration times the
                            ! rescaling factor derived from eta_to_m [T2 Z L-2 ~> s2 m-1]
  real :: Z_to_eta, H_to_eta, H_to_rho_eta ! Unit conversion factors with obvious names.
  integer i, j, k, isv, iev, jsv, jev, nz, halo

  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)

  isv = G%isc-halo ; iev = G%iec+halo ; jsv = G%jsc-halo ; jev = G%jec+halo
  nz = GV%ke

  if ((isv<G%isd) .or. (iev>G%ied) .or. (jsv<G%jsd) .or. (jev>G%jed)) &
    call MOM_error(FATAL,"find_eta called with an overly large halo_size.")

  Z_to_eta = 1.0 ; if (present(eta_to_m)) Z_to_eta = US%Z_to_m / eta_to_m
  H_to_eta = GV%H_to_Z * Z_to_eta
  H_to_rho_eta =  GV%H_to_RZ * Z_to_eta
  I_gEarth = Z_to_eta / GV%g_Earth

!$OMP parallel default(shared) private(dilate,htot)
!$OMP do
  do j=jsv,jev ; do i=isv,iev ; eta(i,j,nz+1) = -Z_to_eta*G%bathyT(i,j) ; enddo ; enddo

  if (GV%Boussinesq) then
!$OMP do
    do j=jsv,jev ; do k=nz,1,-1 ; do i=isv,iev
      eta(i,j,K) = eta(i,j,K+1) + h(i,j,k)*H_to_eta
    enddo ; enddo ; enddo
    if (present(eta_bt)) then
      ! Dilate the water column to agree with the free surface height
      ! that is used for the dynamics.
!$OMP do
      do j=jsv,jev
        do i=isv,iev
          dilate(i) = (eta_bt(i,j)*H_to_eta + Z_to_eta*G%bathyT(i,j)) / &
                      (eta(i,j,1) + Z_to_eta*G%bathyT(i,j))
        enddo
        do k=1,nz ; do i=isv,iev
          eta(i,j,K) = dilate(i) * (eta(i,j,K) + Z_to_eta*G%bathyT(i,j)) - Z_to_eta*G%bathyT(i,j)
        enddo ; enddo
      enddo
    endif
  else
    if (associated(tv%eqn_of_state)) then
!$OMP do
      do j=jsv,jev
        if (associated(tv%p_surf)) then
          do i=isv,iev ; p(i,j,1) = tv%p_surf(i,j) ; enddo
        else
          do i=isv,iev ; p(i,j,1) = 0.0 ; enddo
        endif
        do k=1,nz ; do i=isv,iev
          p(i,j,K+1) = p(i,j,K) + GV%g_Earth*GV%H_to_RZ*h(i,j,k)
        enddo ; enddo
      enddo
!$OMP do
      do k=1,nz
        call int_specific_vol_dp(tv%T(:,:,k), tv%S(:,:,k), p(:,:,K), p(:,:,K+1), &
                                 0.0, G%HI, tv%eqn_of_state, US, dz_geo(:,:,k), halo_size=halo)
      enddo
!$OMP do
      do j=jsv,jev
        do k=nz,1,-1 ; do i=isv,iev
          eta(i,j,K) = eta(i,j,K+1) + I_gEarth * dz_geo(i,j,k)
        enddo ; enddo
      enddo
    else
!$OMP do
      do j=jsv,jev ;  do k=nz,1,-1 ; do i=isv,iev
        eta(i,j,K) = eta(i,j,K+1) + H_to_rho_eta*h(i,j,k) / GV%Rlay(k)
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
          eta(i,j,K) = dilate(i) * (eta(i,j,K) + Z_to_eta*G%bathyT(i,j)) - Z_to_eta*G%bathyT(i,j)
        enddo ; enddo
      enddo
    endif
  endif
!$OMP end parallel

end subroutine find_eta_3d

!> Calculates the free surface height, using the appropriate form for consistency
!! with the calculation of the pressure gradient forces.  Additionally, the sea
!! surface height may be adjusted for consistency with the corresponding
!! time-average quantity from the barotropic calculation.
subroutine find_eta_2d(h, tv, G, GV, US, eta, eta_bt, halo_size, eta_to_m)
  type(ocean_grid_type),                      intent(in)  :: G   !< The ocean's grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< A structure pointing to various
                                                                 !! thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G)),           intent(out) :: eta !< free surface height relative to
                                                                 !! mean sea level (z=0) often [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)  :: eta_bt !< optional barotropic
                   !! variable that gives the "correct" free surface height (Boussinesq) or total
                   !! water column mass per unit area (non-Boussinesq) [H ~> m or kg m-2].
  integer,                          optional, intent(in)  :: halo_size !< width of halo points on
                                                                 !! which to calculate eta.
  real,                             optional, intent(in)  :: eta_to_m  !< The conversion factor from
             !! the units of eta to m; by default this is US%Z_to_m.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: &
    p          ! Hydrostatic pressure at each interface [R L2 T-2 ~> Pa]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    dz_geo     ! The change in geopotential height across a layer [L2 T-2 ~> m2 s-2].
  real :: htot(SZI_(G))  ! The sum of all layers' thicknesses [H ~> m or kg m-2].
  real :: I_gEarth          ! The inverse of the gravitational acceleration times the
                            ! rescaling factor derived from eta_to_m [T2 Z L-2 ~> s2 m-1]
  real :: Z_to_eta, H_to_eta, H_to_rho_eta ! Unit conversion factors with obvious names.
  integer i, j, k, is, ie, js, je, nz, halo

  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)
  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  nz = GV%ke

  Z_to_eta = 1.0 ; if (present(eta_to_m)) Z_to_eta = US%Z_to_m / eta_to_m
  H_to_eta = GV%H_to_Z * Z_to_eta
  H_to_rho_eta =  GV%H_to_RZ * Z_to_eta
  I_gEarth = Z_to_eta / GV%g_Earth

!$OMP parallel default(shared) private(htot)
!$OMP do
  do j=js,je ; do i=is,ie ; eta(i,j) = -Z_to_eta*G%bathyT(i,j) ; enddo ; enddo

  if (GV%Boussinesq) then
    if (present(eta_bt)) then
!$OMP do
      do j=js,je ; do i=is,ie
        eta(i,j) = H_to_eta*eta_bt(i,j)
      enddo ; enddo
    else
!$OMP do
      do j=js,je ; do k=1,nz ; do i=is,ie
        eta(i,j) = eta(i,j) + h(i,j,k)*H_to_eta
      enddo ; enddo ; enddo
    endif
  else
    if (associated(tv%eqn_of_state)) then
!$OMP do
      do j=js,je
        if (associated(tv%p_surf)) then
          do i=is,ie ; p(i,j,1) = tv%p_surf(i,j) ; enddo
        else
          do i=is,ie ; p(i,j,1) = 0.0 ; enddo
        endif

        do k=1,nz ; do i=is,ie
          p(i,j,k+1) = p(i,j,k) + GV%g_Earth*GV%H_to_RZ*h(i,j,k)
        enddo ; enddo
      enddo
!$OMP do
      do k = 1, nz
        call int_specific_vol_dp(tv%T(:,:,k), tv%S(:,:,k), p(:,:,k), p(:,:,k+1), 0.0, &
                                 G%HI, tv%eqn_of_state, US, dz_geo(:,:,k), halo_size=halo)
      enddo
!$OMP do
      do j=js,je ; do k=1,nz ; do i=is,ie
          eta(i,j) = eta(i,j) + I_gEarth * dz_geo(i,j,k)
      enddo ; enddo ; enddo
    else
!$OMP do
      do j=js,je ; do k=1,nz ; do i=is,ie
        eta(i,j) = eta(i,j) + H_to_rho_eta*h(i,j,k) / GV%Rlay(k)
      enddo ; enddo ; enddo
    endif
    if (present(eta_bt)) then
      !   Dilate the water column to agree with the time-averaged column
      ! mass from the barotropic solution.
!$OMP do
      do j=js,je
        do i=is,ie ; htot(i) = GV%H_subroundoff ; enddo
        do k=1,nz ; do i=is,ie ; htot(i) = htot(i) + h(i,j,k) ; enddo ; enddo
        do i=is,ie
          eta(i,j) = (eta_bt(i,j) / htot(i)) * (eta(i,j) + Z_to_eta*G%bathyT(i,j)) - &
                     Z_to_eta*G%bathyT(i,j)
        enddo
      enddo
    endif
  endif
!$OMP end parallel

end subroutine find_eta_2d

end module MOM_interface_heights
