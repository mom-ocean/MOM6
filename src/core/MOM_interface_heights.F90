!> Functions for calculating interface heights, including free surface height.
module MOM_interface_heights

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_density_integrals, only : int_specific_vol_dp, avg_specific_vol
use MOM_debugging,     only : hchksum
use MOM_error_handler, only : MOM_error, FATAL
use MOM_EOS,           only : calculate_density, average_specific_vol, EOS_type, EOS_domain
use MOM_file_parser,   only : log_version
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public find_eta, dz_to_thickness, thickness_to_dz, dz_to_thickness_simple
public calc_derived_thermo
public find_rho_bottom, find_col_avg_SpV

!> Calculates the heights of the free surface or all interfaces from layer thicknesses.
interface find_eta
  module procedure find_eta_2d, find_eta_3d
end interface find_eta

!> Calculates layer thickness in thickness units from geometric distance between the
!! interfaces around that layer in height units.
interface dz_to_thickness
  module procedure dz_to_thickness_tv, dz_to_thickness_EoS
end interface dz_to_thickness

!> Converts layer thickness in thickness units into the vertical distance between the
!! interfaces around a layer in height units.
interface thickness_to_dz
  module procedure thickness_to_dz_3d, thickness_to_dz_jslice
end interface thickness_to_dz

contains

!> Calculates the heights of all interfaces between layers, using the appropriate
!! form for consistency with the calculation of the pressure gradient forces.
!! Additionally, these height may be dilated for consistency with the
!! corresponding time-average quantity from the barotropic calculation.
subroutine find_eta_3d(h, tv, G, GV, US, eta, eta_bt, halo_size, dZref)
  type(ocean_grid_type),                      intent(in)  :: G   !< The ocean's grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< A structure pointing to various
                                                                 !! thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(out) :: eta !< layer interface heights [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)  :: eta_bt !< optional barotropic variable
                    !! that gives the "correct" free surface height (Boussinesq) or total water
                    !! column mass per unit area (non-Boussinesq).  This is used to dilate the layer
                    !! thicknesses when calculating interface heights [H ~> m or kg m-2].
                    !! In Boussinesq mode, eta_bt and G%bathyT use the same reference height.
  integer,                          optional, intent(in)  :: halo_size !< width of halo points on
                                                                 !! which to calculate eta.
  real,                             optional, intent(in)  :: dZref !< The difference in the
                    !! reference height between G%bathyT and eta [Z ~> m]. The default is 0.

  ! Local variables
  real :: p(SZI_(G),SZJ_(G),SZK_(GV)+1)   ! Hydrostatic pressure at each interface [R L2 T-2 ~> Pa]
  real :: dz_geo(SZI_(G),SZJ_(G),SZK_(GV)) ! The change in geopotential height
                                           ! across a layer [L2 T-2 ~> m2 s-2].
  real :: dilate(SZI_(G))                 ! non-dimensional dilation factor
  real :: htot(SZI_(G))                   ! total thickness [H ~> m or kg m-2]
  real :: I_gEarth          ! The inverse of the gravitational acceleration times the
                            ! rescaling factor derived from eta_to_m [T2 Z L-2 ~> s2 m-1]
  real :: dZ_ref    ! The difference in the reference height between G%bathyT and eta [Z ~> m].
                    ! dZ_ref is 0 unless the optional argument dZref is present.
  integer i, j, k, isv, iev, jsv, jev, nz, halo

  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)

  isv = G%isc-halo ; iev = G%iec+halo ; jsv = G%jsc-halo ; jev = G%jec+halo
  nz = GV%ke

  if ((isv<G%isd) .or. (iev>G%ied) .or. (jsv<G%jsd) .or. (jev>G%jed)) &
    call MOM_error(FATAL,"find_eta called with an overly large halo_size.")

  I_gEarth = 1.0 / GV%g_Earth
  dZ_ref = 0.0 ; if (present(dZref)) dZ_ref = dZref

  !$OMP parallel default(shared) private(dilate,htot)
  !$OMP do
  do j=jsv,jev ; do i=isv,iev ; eta(i,j,nz+1) = -(G%bathyT(i,j) + dZ_ref) ; enddo ; enddo

  if (GV%Boussinesq) then
    !$OMP do
    do j=jsv,jev ; do k=nz,1,-1 ; do i=isv,iev
      eta(i,j,K) = eta(i,j,K+1) + h(i,j,k)*GV%H_to_Z
    enddo ; enddo ; enddo
    if (present(eta_bt)) then
      ! Dilate the water column to agree with the free surface height
      ! that is used for the dynamics.
      !$OMP do
      do j=jsv,jev
        do i=isv,iev
          dilate(i) = (eta_bt(i,j)*GV%H_to_Z + G%bathyT(i,j)) / &
                      (eta(i,j,1) + (G%bathyT(i,j) + dZ_ref))
        enddo
        do k=1,nz ; do i=isv,iev
          eta(i,j,K) = dilate(i) * (eta(i,j,K) + (G%bathyT(i,j) + dZ_ref)) - &
                       (G%bathyT(i,j) + dZ_ref)
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
        eta(i,j,K) = eta(i,j,K+1) + GV%H_to_RZ*h(i,j,k) / GV%Rlay(k)
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
          eta(i,j,K) = dilate(i) * (eta(i,j,K) + (G%bathyT(i,j) + dZ_ref)) - &
                       (G%bathyT(i,j) + dZ_ref)
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
subroutine find_eta_2d(h, tv, G, GV, US, eta, eta_bt, halo_size, dZref)
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
                    !! In Boussinesq mode, eta_bt and G%bathyT use the same reference height.
  integer,                          optional, intent(in)  :: halo_size !< width of halo points on
                                                                 !! which to calculate eta.
  real,                             optional, intent(in)  :: dZref !< The difference in the
                    !! reference height between G%bathyT and eta [Z ~> m]. The default is 0.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: &
    p          ! Hydrostatic pressure at each interface [R L2 T-2 ~> Pa]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    dz_geo     ! The change in geopotential height across a layer [L2 T-2 ~> m2 s-2].
  real :: htot(SZI_(G))  ! The sum of all layers' thicknesses [H ~> m or kg m-2].
  real :: I_gEarth          ! The inverse of the gravitational acceleration times the
                            ! rescaling factor derived from eta_to_m [T2 Z L-2 ~> s2 m-1]
  real :: dZ_ref    ! The difference in the reference height between G%bathyT and eta [Z ~> m].
                    ! dZ_ref is 0 unless the optional argument dZref is present.
  integer i, j, k, is, ie, js, je, nz, halo

  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)
  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  nz = GV%ke

  I_gEarth = 1.0 / GV%g_Earth
  dZ_ref = 0.0 ; if (present(dZref)) dZ_ref = dZref

  !$OMP parallel default(shared) private(htot)
  !$OMP do
  do j=js,je ; do i=is,ie ; eta(i,j) = -(G%bathyT(i,j) + dZ_ref) ; enddo ; enddo

  if (GV%Boussinesq) then
    if (present(eta_bt)) then
      !$OMP do
      do j=js,je ; do i=is,ie
        eta(i,j) = GV%H_to_Z*eta_bt(i,j) - dZ_ref
      enddo ; enddo
    else
      !$OMP do
      do j=js,je ; do k=1,nz ; do i=is,ie
        eta(i,j) = eta(i,j) + h(i,j,k)*GV%H_to_Z
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
        eta(i,j) = eta(i,j) + GV%H_to_RZ*h(i,j,k) / GV%Rlay(k)
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
          eta(i,j) = (eta_bt(i,j) / htot(i)) * (eta(i,j) + (G%bathyT(i,j) + dZ_ref)) - &
                     (G%bathyT(i,j) + dZ_ref)
        enddo
      enddo
    endif
  endif
  !$OMP end parallel

end subroutine find_eta_2d


!> Calculate derived thermodynamic quantities for re-use later.
subroutine calc_derived_thermo(tv, h, G, GV, US, halo, debug)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(thermo_var_ptrs),   intent(inout) :: tv !< A structure pointing to various
                                               !! thermodynamic variables, some of
                                               !! which will be set here.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2].
  integer,       optional, intent(in)    :: halo !< Width of halo within which to
                                               !! calculate thicknesses
  logical,       optional, intent(in)    :: debug !< If present and true, write debugging checksums
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: p_t  ! Hydrostatic pressure atop a layer [R L2 T-2 ~> Pa]
  real, dimension(SZI_(G),SZJ_(G)) :: dp   ! Pressure change across a layer [R L2 T-2 ~> Pa]
  real, dimension(SZK_(GV)) :: SpV_lay     ! The specific volume of each layer when no equation of
                                           ! state is used [R-1 ~> m3 kg-1]
  logical :: do_debug  ! If true, write checksums for debugging.
  integer :: i, j, k, is, ie, js, je, halos, nz

  do_debug = .false. ; if (present(debug)) do_debug = debug
  halos = 0 ; if (present(halo)) halos = max(0,halo)
  is = G%isc-halos ; ie = G%iec+halos ; js = G%jsc-halos ; je = G%jec+halos ; nz = GV%ke

  if (allocated(tv%Spv_avg) .and. associated(tv%eqn_of_state)) then
    if (associated(tv%p_surf)) then
      do j=js,je ; do i=is,ie ; p_t(i,j) = tv%p_surf(i,j) ; enddo ; enddo
    else
      do j=js,je ; do i=is,ie ; p_t(i,j) = 0.0 ; enddo ; enddo
    endif
    do k=1,nz
      do j=js,je ; do i=is,ie
        dp(i,j) = GV%g_Earth*GV%H_to_RZ*h(i,j,k)
      enddo ; enddo
      call avg_specific_vol(tv%T(:,:,k), tv%S(:,:,k), p_t, dp, G%HI, tv%eqn_of_state, tv%SpV_avg(:,:,k), halo)
      if (k<nz) then ; do j=js,je ; do i=is,ie
        p_t(i,j) = p_t(i,j) + dp(i,j)
      enddo ; enddo ; endif
    enddo
    tv%valid_SpV_halo = halos

    if (do_debug) then
      call hchksum(h, "derived_thermo h", G%HI, haloshift=halos, scale=GV%H_to_MKS)
      if (associated(tv%p_surf)) call hchksum(tv%p_surf, "derived_thermo p_surf", G%HI, &
                                              haloshift=halos, scale=US%RL2_T2_to_Pa)
      call hchksum(tv%T, "derived_thermo T", G%HI, haloshift=halos, scale=US%C_to_degC)
      call hchksum(tv%S, "derived_thermo S", G%HI, haloshift=halos, scale=US%S_to_ppt)
    endif
  elseif (allocated(tv%Spv_avg)) then
    do k=1,nz ; SpV_lay(k) = 1.0 / GV%Rlay(k) ; enddo
    do k=1,nz ; do j=js,je ; do i=is,ie
      tv%SpV_avg(i,j,k) = SpV_lay(k)
    enddo ; enddo ; enddo
    tv%valid_SpV_halo = halos
  endif

end subroutine calc_derived_thermo


!> Determine the column average specific volumes.
subroutine find_col_avg_SpV(h, SpV_avg, tv, G, GV, US, halo_size)
  type(ocean_grid_type),    intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),  intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), &
                            intent(inout) :: SpV_avg !< Column average specific volume [R-1 ~> m3 kg-1]
                                                  ! SpV_avg is intent inout to retain excess halo values.
  type(thermo_var_ptrs),    intent(in)    :: tv   !< Structure containing pointers to any available
                                                  !! thermodynamic fields.
  integer,        optional, intent(in)    :: halo_size !< width of halo points on which to work

  ! Local variables
  real :: h_tot(SZI_(G))        ! Sum of the layer thicknesses [H ~> m or kg m-3]
  real :: SpV_x_h_tot(SZI_(G))  ! Vertical sum of the layer average specific volume times
                                ! the layer thicknesses [H R-1 ~> m4 kg-1 or m]
  real :: I_rho                 ! The inverse of the Boussiensq reference density [R-1 ~> m3 kg-1]
  real :: SpV_lay(SZK_(GV))     ! The inverse of the layer target potential densities [R-1 ~> m3 kg-1]
  character(len=128) :: mesg    ! A string for error messages
  integer i, j, k, is, ie, js, je, nz, halo

  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)

  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  nz = GV%ke

  if (GV%Boussinesq) then
    I_rho = 1.0 / GV%Rho0
    do j=js,je ; do i=is,ie
      SpV_avg(i,j) = I_rho
    enddo ; enddo
  elseif (.not.allocated(tv%SpV_avg)) then
    do k=1,nz ; Spv_lay(k) = 1.0 / GV%Rlay(k) ; enddo
    do j=js,je
      do i=is,ie ; SpV_x_h_tot(i) = 0.0 ; h_tot(i) = 0.0 ; enddo
      do k=1,nz ; do i=is,ie
        h_tot(i) = h_tot(i) + max(h(i,j,k), GV%H_subroundoff)
        SpV_x_h_tot(i) = SpV_x_h_tot(i) + Spv_lay(k)*max(h(i,j,k), GV%H_subroundoff)
      enddo ; enddo
      do i=is,ie ; SpV_avg(i,j) = SpV_x_h_tot(i) / h_tot(i) ; enddo
    enddo
  else
    ! Check that SpV_avg has been set.
    if ((allocated(tv%SpV_avg)) .and. (tv%valid_SpV_halo < halo)) then
      if (tv%valid_SpV_halo < 0) then
        mesg = "invalid values of SpV_avg."
      else
        write(mesg, '("insufficiently large SpV_avg halos of width ", i2, " but ", i2," is needed.")') &
                     tv%valid_SpV_halo, halo
      endif
      call MOM_error(FATAL, "find_col_avg_SpV called in fully non-Boussinesq mode with "//trim(mesg))
    endif

    do j=js,je
      do i=is,ie ; SpV_x_h_tot(i) = 0.0 ; h_tot(i) = 0.0 ; enddo
      do k=1,nz ; do i=is,ie
        h_tot(i) = h_tot(i) + max(h(i,j,k), GV%H_subroundoff)
        SpV_x_h_tot(i) = SpV_x_h_tot(i) + tv%SpV_avg(i,j,k)*max(h(i,j,k), GV%H_subroundoff)
      enddo ; enddo
      do i=is,ie ; SpV_avg(i,j) = SpV_x_h_tot(i) / h_tot(i) ; enddo
    enddo
  endif

end subroutine find_col_avg_SpV


!> Determine the in situ density averaged over a specified distance from the bottom,
!! calculating it as the inverse of the mass-weighted average specific volume.
subroutine find_rho_bottom(h, dz, pres_int, dz_avg, tv, j, G, GV, US, Rho_bot)
  type(ocean_grid_type),    intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type),  intent(in)  :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),    intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZK_(GV)), &
                            intent(in)  :: dz   !< Height change across layers [Z ~> m]
  real, dimension(SZI_(G),SZK_(GV)+1), &
                            intent(in)  :: pres_int !< Pressure at each interface [R L2 T-2 ~> Pa]
  real, dimension(SZI_(G)), intent(in)  :: dz_avg !< The vertical distance over which to average [Z ~> m]
  type(thermo_var_ptrs),    intent(in)  :: tv   !< Structure containing pointers to any available
                                                !! thermodynamic fields.
  integer,                  intent(in)  :: j    !< j-index of row to work on
  real, dimension(SZI_(G)), intent(out) :: Rho_bot  !< Near-bottom density [R ~> kg m-3].

  ! Local variables
  real :: hb(SZI_(G))         ! Running sum of the thickness in the bottom boundary layer [H ~> m or kg m-2]
  real :: SpV_h_bot(SZI_(G))  ! Running sum of the specific volume times thickness in the bottom
                              ! boundary layer [R-1 H ~> m4 kg-1 or m]
  real :: dz_bbl_rem(SZI_(G)) ! Vertical extent of the boundary layer that has yet to be accounted
                              ! for [Z ~> m]
  real :: h_bbl_frac(SZI_(G)) ! Thickness of the fractional layer that makes up the top of the
                              ! boundary layer [H ~> m or kg m-2]
  real :: T_bbl(SZI_(G))      ! Temperature of the fractional layer that makes up the top of the
                              ! boundary layer [C ~> degC]
  real :: S_bbl(SZI_(G))      ! Salinity of the fractional layer that makes up the top of the
                              ! boundary layer [S ~> ppt]
  real :: P_bbl(SZI_(G))      ! Pressure the top of the boundary layer [R L2 T-2 ~> Pa]
  real :: dp(SZI_(G))         ! Pressure change across the fractional layer that makes up the top
                              ! of the boundary layer [R L2 T-2 ~> Pa]
  real :: SpV_bbl(SZI_(G))    ! In situ specific volume of the fractional layer that makes up the
                              ! top of the boundary layer [R-1 ~> m3 kg-1]
  real :: frac_in             ! The fraction of a layer that is within the bottom boundary layer [nondim]
  logical :: do_i(SZI_(G)), do_any
  logical :: use_EOS
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, k, is, ie, nz

  is = G%isc ; ie = G%iec ; nz = GV%ke

  use_EOS = associated(tv%T) .and. associated(tv%S) .and. associated(tv%eqn_of_state)

  if (GV%Boussinesq .or. GV%semi_Boussinesq .or. .not.allocated(tv%SpV_avg)) then
    do i=is,ie
      rho_bot(i) = GV%Rho0
    enddo
  else
    ! Check that SpV_avg has been set.
    if (tv%valid_SpV_halo < 0) call MOM_error(FATAL, &
        "find_rho_bottom called in fully non-Boussinesq mode with invalid values of SpV_avg.")

    ! Set the bottom density to the inverse of the in situ specific volume averaged over the
    ! specified distance, with care taken to avoid having compressibility lead to an imprint
    ! of the layer thicknesses on this density.
    do i=is,ie
      hb(i) = 0.0 ; SpV_h_bot(i) = 0.0
      dz_bbl_rem(i) = G%mask2dT(i,j) * max(0.0, dz_avg(i))
      do_i(i) = .true.
      if (G%mask2dT(i,j) <= 0.0) then
        ! Set acceptable values for calling the equation of state over land.
        T_bbl(i) = 0.0 ; S_bbl(i) = 0.0 ; dp(i) = 0.0 ; P_bbl(i) = 0.0
        SpV_bbl(i) = 1.0 ! This value is arbitrary, provided it is non-zero.
        h_bbl_frac(i) = 0.0
        do_i(i) = .false.
      endif
    enddo

    do k=nz,1,-1
      do_any = .false.
      do i=is,ie ; if (do_i(i)) then
        if (dz(i,k) < dz_bbl_rem(i)) then
          ! This layer is fully within the averaging depth.
          SpV_h_bot(i) = SpV_h_bot(i) + h(i,j,k) * tv%SpV_avg(i,j,k)
          dz_bbl_rem(i) = dz_bbl_rem(i) - dz(i,k)
          hb(i) = hb(i) + h(i,j,k)
          do_any = .true.
        else
          if (dz(i,k) > 0.0) then
            frac_in = dz_bbl_rem(i) / dz(i,k)
          else
            frac_in = 0.0
          endif
          if (use_EOS) then
            ! Store the properties of this layer to determine the average
            ! specific volume of the portion that is within the BBL.
            T_bbl(i) = tv%T(i,j,k) ; S_bbl(i) = tv%S(i,j,k)
            dp(i) = frac_in * (GV%g_Earth*GV%H_to_RZ * h(i,j,k))
            P_bbl(i) = pres_int(i,K) + (1.0-frac_in) * (GV%g_Earth*GV%H_to_RZ * h(i,j,k))
          else
            SpV_bbl(i) = tv%SpV_avg(i,j,k)
          endif
          h_bbl_frac(i) = frac_in * h(i,j,k)
          dz_bbl_rem(i) = 0.0
          do_i(i) = .false.
        endif
      endif ; enddo
      if (.not.do_any) exit
    enddo
    do i=is,ie ; if (do_i(i)) then
      ! The nominal bottom boundary layer is thicker than the water column, but layer 1 is
      ! already included in the averages.  These values are set so that the call to find
      ! the layer-average specific volume will behave sensibly.
      if (use_EOS) then
        T_bbl(i) = tv%T(i,j,1) ; S_bbl(i) = tv%S(i,j,1)
        dp(i) = 0.0
        P_bbl(i) = pres_int(i,1)
      else
        SpV_bbl(i) = tv%SpV_avg(i,j,1)
      endif
      h_bbl_frac(i) = 0.0
    endif ; enddo

    if (use_EOS) then
      ! Find the average specific volume of the fractional layer atop the BBL.
      EOSdom(:) = EOS_domain(G%HI)
      call average_specific_vol(T_bbl, S_bbl, P_bbl, dp, SpV_bbl, tv%eqn_of_state, EOSdom)
    endif

    do i=is,ie
      if (hb(i) + h_bbl_frac(i) < GV%H_subroundoff) h_bbl_frac(i) = GV%H_subroundoff
      rho_bot(i) = G%mask2dT(i,j) * (hb(i) + h_bbl_frac(i)) / (SpV_h_bot(i) + h_bbl_frac(i)*SpV_bbl(i))
    enddo
  endif

end subroutine find_rho_bottom


!> Converts thickness from geometric height units to thickness units, perhaps via an
!! inversion of the integral of the density in pressure using variables stored in
!! the thermo_var_ptrs type when in non-Boussinesq mode.
subroutine dz_to_thickness_tv(dz, tv, h, G, GV, US, halo_size)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dz !< Geometric layer thicknesses in height units [Z ~> m]
  type(thermo_var_ptrs),   intent(in)    :: tv !< A structure pointing to various
                                               !! thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h  !< Output thicknesses in thickness units [H ~> m or kg m-2].
                                               !! This is essentially intent out, but declared as intent
                                               !! inout to preserve any initialized values in halo points.
  integer,         optional, intent(in)  :: halo_size !< Width of halo within which to
                                               !! calculate thicknesses
  ! Local variables
  integer :: i, j, k, is, ie, js, je, halo, nz

  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)
  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo ; nz = GV%ke

  if (GV%Boussinesq) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      h(i,j,k) = GV%Z_to_H * dz(i,j,k)
    enddo ; enddo ; enddo
  else
    if (associated(tv%eqn_of_state)) then
      if (associated(tv%p_surf)) then
        call dz_to_thickness_EOS(dz, tv%T, tv%S, tv%eqn_of_state, h, G, GV, US, halo, tv%p_surf)
      else
        call dz_to_thickness_EOS(dz, tv%T, tv%S, tv%eqn_of_state, h, G, GV, US, halo)
      endif
    else
      do k=1,nz ; do j=js,je ; do i=is,ie
        h(i,j,k) = (GV%RZ_to_H * GV%Rlay(k)) * dz(i,j,k)
      enddo ; enddo ; enddo
    endif
  endif

end subroutine dz_to_thickness_tv

!> Converts thickness from geometric height units to thickness units, working via an
!! inversion of the integral of the density in pressure when in non-Boussinesq mode.
subroutine dz_to_thickness_EOS(dz, Temp, Saln, EoS, h, G, GV, US, halo_size, p_surf)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dz !< Geometric layer thicknesses in height units [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: Temp !< Input layer temperatures [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: Saln !< Input layer salinities [S ~> ppt]
  type(EOS_type),          intent(in)    :: EoS  !< Equation of state structure
    real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h  !< Output thicknesses in thickness units [H ~> m or kg m-2].
                                               !! This is essentially intent out, but declared as intent
                                               !! inout to preserve any initialized values in halo points.
  integer,         optional, intent(in)  :: halo_size !< Width of halo within which to
                                               !! calculate thicknesses
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)  :: p_surf !< Surface pressures [R L2 T-2 ~> Pa]
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    p_top, p_bot                  ! Pressure at the interfaces above and below a layer [R L2 T-2 ~> Pa]
  real :: dp(SZI_(G),SZJ_(G))     ! Pressure change across a layer [R L2 T-2 ~> Pa]
  real :: dz_geo(SZI_(G),SZJ_(G)) ! The change in geopotential height across a layer [L2 T-2 ~> m2 s-2]
  real :: rho(SZI_(G))            ! The in situ density [R ~> kg m-3]
  real :: dp_adj                  ! The amount by which to change the bottom pressure in an
                                  ! iteration [R L2 T-2 ~> Pa]
  real :: I_gEarth                ! Unit conversion factors divided by the gravitational
                                  ! acceleration [H T2 R-1 L-2 ~> s2 m2 kg-1 or s2 m-1]
  logical :: do_more(SZI_(G),SZJ_(G)) ! If true, additional iterations would be beneficial.
  logical :: do_any               ! True if there are points in this layer that need more itertions.
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, k, is, ie, js, je, halo, nz
  integer :: itt, max_itt

  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)
  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo ; nz = GV%ke
  max_itt = 10

  if (GV%Boussinesq) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      h(i,j,k) = GV%Z_to_H * dz(i,j,k)
    enddo ; enddo ; enddo
  else
    I_gEarth = GV%RZ_to_H / GV%g_Earth

    if (present(p_surf)) then
      do j=js,je ; do i=is,ie
        p_bot(i,j) = 0.0 ; p_top(i,j) = p_surf(i,j)
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        p_bot(i,j) = 0.0 ; p_top(i,j) = 0.0
      enddo ; enddo
    endif
    EOSdom(:) = EOS_domain(G%HI)

    ! The iterative approach here is inherited from very old code that was in the
    ! MOM_state_initialization module.  It does converge, but it is very inefficient and
    ! should be revised, although doing so would change answers in non-Boussinesq mode.
    do k=1,nz
      do j=js,je
        do i=is,ie ; p_top(i,j) = p_bot(i,j) ; enddo
        call calculate_density(Temp(:,j,k), Saln(:,j,k), p_top(:,j), rho, &
                               EoS, EOSdom)
        ! The following two expressions are mathematically equivalent.
        if (GV%semi_Boussinesq) then
          do i=is,ie
            p_bot(i,j) = p_top(i,j) + (GV%g_Earth*GV%H_to_Z) * ((GV%Z_to_H*dz(i,j,k)) * rho(i))
            dp(i,j) = (GV%g_Earth*GV%H_to_Z) * ((GV%Z_to_H*dz(i,j,k)) * rho(i))
          enddo
        else
          do i=is,ie
            p_bot(i,j) = p_top(i,j) + rho(i) * (GV%g_Earth * dz(i,j,k))
            dp(i,j) = rho(i) * (GV%g_Earth * dz(i,j,k))
          enddo
        endif
      enddo

      do_more(:,:) = .true.
      do itt=1,max_itt
        do_any = .false.
        call int_specific_vol_dp(Temp(:,:,k), Saln(:,:,k), p_top, p_bot, 0.0, G%HI, EoS, US, dz_geo)
        if (itt < max_itt) then ; do j=js,je
          call calculate_density(Temp(:,j,k), Saln(:,j,k), p_bot(:,j), rho, EoS, EOSdom)
          ! Use Newton's method to correct the bottom value.
          ! The hydrostatic equation is sufficiently linear that no bounds-checking is needed.
          if (GV%semi_Boussinesq) then
            do i=is,ie
              dp_adj = rho(i) * ((GV%g_Earth*GV%H_to_Z)*(GV%Z_to_H*dz(i,j,k)) - dz_geo(i,j))
              p_bot(i,j) = p_bot(i,j) + dp_adj
              dp(i,j) = dp(i,j) + dp_adj
            enddo
            do_any = .true. ! To avoid changing answers, always use the maximum number of itertions.
          else
            do i=is,ie ; if (do_more(i,j)) then
              dp_adj = rho(i) * (GV%g_Earth*dz(i,j,k) - dz_geo(i,j))
              p_bot(i,j) = p_bot(i,j) + dp_adj
              dp(i,j) = dp(i,j) + dp_adj
              ! Check for convergence to roundoff.
              do_more(i,j) = (abs(dp_adj) > 1.0e-15*dp(i,j))
              if (do_more(i,j)) do_any = .true.
            endif ; enddo
          endif
        enddo ; endif
        if (.not.do_any) exit
      enddo

      if (GV%semi_Boussinesq) then
        do j=js,je ; do i=is,ie
          h(i,j,k) = (p_bot(i,j) - p_top(i,j)) * I_gEarth
        enddo ; enddo
      else
        do j=js,je ; do i=is,ie
          h(i,j,k) = dp(i,j) * I_gEarth
        enddo ; enddo
      endif
    enddo
  endif

end subroutine dz_to_thickness_EOS

!> Converts thickness from geometric height units to thickness units, perhaps using
!! a simple conversion factor that may be problematic in non-Boussinesq mode.
subroutine dz_to_thickness_simple(dz, h, G, GV, US, halo_size, layer_mode)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dz !< Geometric layer thicknesses in height units [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h  !< Output thicknesses in thickness units [H ~> m or kg m-2].
                                               !! This is essentially intent out, but declared as intent
                                               !! inout to preserve any initialized values in halo points.
  integer,         optional, intent(in)  :: halo_size !< Width of halo within which to
                                               !! calculate thicknesses
  logical,         optional, intent(in)  :: layer_mode !< If present and true, do the conversion that
                                               !! is appropriate in pure isopycnal layer mode with
                                               !! no state variables or equation of state.  Otherwise
                                               !! use a simple constant rescaling factor and avoid the
                                               !! use of GV%Rlay.
  ! Local variables
  logical :: layered  ! If true and the model is non-Boussinesq, do calculations appropriate for use
                      ! in pure isopycnal layered mode with no state variables or equation of state.
  integer :: i, j, k, is, ie, js, je, halo, nz

  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)
  layered = .false. ; if (present(layer_mode)) layered = layer_mode
  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo ; nz = GV%ke

  if (GV%Boussinesq) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      h(i,j,k) = GV%Z_to_H * dz(i,j,k)
    enddo ; enddo ; enddo
  elseif (layered) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      h(i,j,k) = (GV%RZ_to_H * GV%Rlay(k)) * dz(i,j,k)
    enddo ; enddo ; enddo
  else
    do k=1,nz ; do j=js,je ; do i=is,ie
      h(i,j,k) = (US%Z_to_m * GV%m_to_H) * dz(i,j,k)
    enddo ; enddo ; enddo
  endif

end subroutine dz_to_thickness_simple

!> Converts layer thicknesses in thickness units to the vertical distance between edges in height
!! units, perhaps by multiplication by the precomputed layer-mean specific volume stored in an
!! array in the thermo_var_ptrs type when in non-Boussinesq mode.
subroutine thickness_to_dz_3d(h, tv, dz, G, GV, US, halo_size)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Input thicknesses in thickness units [H ~> m or kg m-2].
  type(thermo_var_ptrs),   intent(in)    :: tv !< A structure pointing to various
                                               !! thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: dz !< Geometric layer thicknesses in height units [Z ~> m]
                                               !! This is essentially intent out, but declared as intent
                                               !! inout to preserve any initialized values in halo points.
  integer,       optional, intent(in)    :: halo_size !< Width of halo within which to
                                               !! calculate thicknesses
  ! Local variables
  character(len=128) :: mesg    ! A string for error messages
  integer :: i, j, k, is, ie, js, je, halo, nz

  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)
  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo ; nz = GV%ke

  if ((.not.GV%Boussinesq) .and. allocated(tv%SpV_avg))  then
    if ((allocated(tv%SpV_avg)) .and. (tv%valid_SpV_halo < halo)) then
      if (tv%valid_SpV_halo < 0) then
        mesg = "invalid values of SpV_avg."
      else
        write(mesg, '("insufficiently large SpV_avg halos of width ", i2, " but ", i2," is needed.")') &
                     tv%valid_SpV_halo, halo
      endif
      call MOM_error(FATAL, "thickness_to_dz called in fully non-Boussinesq mode with "//trim(mesg))
    endif

    do k=1,nz ; do j=js,je ; do i=is,ie
      dz(i,j,k) = GV%H_to_RZ * h(i,j,k) * tv%SpV_avg(i,j,k)
    enddo ; enddo ; enddo
  else
    do k=1,nz ; do j=js,je ; do i=is,ie
      dz(i,j,k) = GV%H_to_Z * h(i,j,k)
    enddo ; enddo ; enddo
  endif

end subroutine thickness_to_dz_3d


!> Converts a vertical i- / k- slice of layer thicknesses in thickness units to the vertical
!! distance between edges in height units, perhaps by multiplication by the precomputed layer-mean
!! specific volume stored in an array in the thermo_var_ptrs type when in non-Boussinesq mode.
subroutine thickness_to_dz_jslice(h, tv, dz, j, G, GV, halo_size)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
   real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Input thicknesses in thickness units [H ~> m or kg m-2].
  type(thermo_var_ptrs),   intent(in)    :: tv !< A structure pointing to various
                                               !! thermodynamic variables
  real, dimension(SZI_(G),SZK_(GV)), &
                           intent(inout) :: dz !< Geometric layer thicknesses in height units [Z ~> m]
                                               !! This is essentially intent out, but declared as intent
                                               !! inout to preserve any initialized values in halo points.
  integer,                 intent(in)    :: j  !< The second (j-) index of the input thicknesses to work with
  integer,       optional, intent(in)    :: halo_size !< Width of halo within which to
                                               !! calculate thicknesses
  ! Local variables
  character(len=128) :: mesg    ! A string for error messages
  integer :: i, k, is, ie, halo, nz

  halo = 0 ; if (present(halo_size)) halo = max(0,halo_size)
  is = G%isc-halo ; ie = G%iec+halo ; nz = GV%ke

  if ((.not.GV%Boussinesq) .and. allocated(tv%SpV_avg))  then
    if ((allocated(tv%SpV_avg)) .and. (tv%valid_SpV_halo < halo)) then
      if (tv%valid_SpV_halo < 0) then
        mesg = "invalid values of SpV_avg."
      else
        write(mesg, '("insufficiently large SpV_avg halos of width ", i2, " but ", i2," is needed.")') &
                     tv%valid_SpV_halo, halo
      endif
      call MOM_error(FATAL, "thickness_to_dz called in fully non-Boussinesq mode with "//trim(mesg))
    endif

    do k=1,nz ; do i=is,ie
      dz(i,k) = GV%H_to_RZ * h(i,j,k) * tv%SpV_avg(i,j,k)
    enddo ; enddo
  else
    do k=1,nz ; do i=is,ie
      dz(i,k) = GV%H_to_Z * h(i,j,k)
    enddo ; enddo
  endif

end subroutine thickness_to_dz_jslice

end module MOM_interface_heights
