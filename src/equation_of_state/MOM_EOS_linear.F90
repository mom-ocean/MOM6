module MOM_EOS_linear

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

!***********************************************************************
!*  The subroutines in this file implement a simple linear equation of *
!*  state for sea water with constant coefficients set as parameters.  *
!***********************************************************************

use MOM_grid, only : ocean_grid_type

implicit none ; private

#include <MOM_memory.h>

public calculate_compress_linear, calculate_density_linear
public calculate_density_derivs_linear, calculate_2_densities_linear
public calculate_density_scalar_linear, calculate_density_array_linear
public int_density_dz_linear, int_spec_vol_dp_linear

interface calculate_density_linear
  module procedure calculate_density_scalar_linear, calculate_density_array_linear
end interface calculate_density_linear

contains

subroutine calculate_density_scalar_linear(T, S, pressure, rho, start, npts, &
                                           Rho_T0_S0, dRho_dT, dRho_dS)
  real,    intent(in)  :: T, S, pressure
  real,    intent(out) :: rho
  integer, intent(in)  :: start, npts
  real,    intent(in)  :: Rho_T0_S0, dRho_dT, dRho_dS
! *  This subroutine computes the density of sea water with a trivial  *
! *  linear equation of state (in kg/m^3) from salinity (sal in psu),  *
! *  potential temperature (T in deg C), and pressure in Pa.           *
! *                                                                    *
! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     rho - in situ density in kg m-3.                        *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
! *  (in)      Rho_T0_S0 - The density at T=0, S=0, in kg m-3.         *
! *  (in)      dRho_dT - The derivatives of density with temperature   *
! *  (in)      dRho_dS - and salinity, in kg m-3 C-1 and kg m-3 psu-1. *

  rho = Rho_T0_S0 + dRho_dT*T + dRho_dS*S

end subroutine calculate_density_scalar_linear

subroutine calculate_density_array_linear(T, S, pressure, rho, start, npts, &
                                          Rho_T0_S0, dRho_dT, dRho_dS)
  real,    intent(in),  dimension(:) :: T, S, pressure
  real,    intent(out), dimension(:) :: rho
  integer, intent(in)                :: start, npts
  real,    intent(in)  :: Rho_T0_S0, dRho_dT, dRho_dS
! *  This subroutine computes the density of sea water with a trivial  *
! *  linear equation of state (in kg/m^3) from salinity (sal in psu),  *
! *  potential temperature (T in deg C), and pressure in Pa.           *
! *                                                                    *
! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     rho - in situ density in kg m-3.                        *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
! *  (in)      Rho_T0_S0 - The density at T=0, S=0, in kg m-3.         *
! *  (in)      dRho_dT - The derivatives of density with temperature   *
! *  (in)      dRho_dS - and salinity, in kg m-3 C-1 and kg m-3 psu-1. *
  real :: al0, p0, lambda
  integer :: j

  do j=start,start+npts-1
    rho(j) = Rho_T0_S0 + dRho_dT*T(j) + dRho_dS*S(j)
  enddo
end subroutine calculate_density_array_linear

subroutine calculate_density_derivs_linear(T, S, pressure, drho_dT_out, &
                       drho_dS_out, start, npts, Rho_T0_S0, dRho_dT, dRho_dS)
  real,    intent(in),  dimension(:) :: T, S, pressure
  real,    intent(out), dimension(:) :: drho_dT_out, drho_dS_out
  integer, intent(in)                :: start, npts
  real,    intent(in)                :: Rho_T0_S0, dRho_dT, dRho_dS
! *   This subroutine calculates the partial derivatives of density    *
! * with potential temperature and salinity.                           *
! *                                                                    *
! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     drho_dT_out - the partial derivative of density with    *
! *                      potential temperature, in kg m-3 K-1.         *
! *  (out)     drho_dS_out - the partial derivative of density with    *
! *                      salinity, in kg m-3 psu-1.                    *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
! *  (in)      Rho_T0_S0 - The density at T=0, S=0, in kg m-3.         *
! *  (in)      dRho_dT - The derivatives of density with temperature   *
! *  (in)      dRho_dS - and salinity, in kg m-3 C-1 and kg m-3 psu-1. *
  integer :: j

  do j=start,start+npts-1
    drho_dT_out(j) = dRho_dT
    drho_dS_out(j) = dRho_dS
  enddo

end subroutine calculate_density_derivs_linear

subroutine calculate_compress_linear(T, S, pressure, rho, drho_dp, start, npts,&
                                     Rho_T0_S0, dRho_dT, dRho_dS)
  real,    intent(in),  dimension(:) :: T, S, pressure
  real,    intent(out), dimension(:) :: rho, drho_dp
  integer, intent(in)                :: start, npts
  real,    intent(in)                :: Rho_T0_S0, dRho_dT, dRho_dS
! *  This subroutine computes the in situ density of sea water (rho)   *
! *  and the compressibility (drho/dp == C_sound^-2) at the given      *
! *  salinity, potential temperature, and pressure.                    *
! *                                                                    *
! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     rho - in situ density in kg m-3.                        *
! *  (out)     drho_dp - the partial derivative of density with        *
! *                      pressure (also the inverse of the square of   *
! *                      sound speed) in s2 m-2.                       *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
! *  (in)      Rho_T0_S0 - The density at T=0, S=0, in kg m-3.         *
! *  (in)      dRho_dT - The derivatives of density with temperature   *
! *  (in)      dRho_dS - and salinity, in kg m-3 C-1 and kg m-3 psu-1. *

  integer :: j

  do j=start,start+npts-1
    rho(j) = Rho_T0_S0 + dRho_dT*T(j) + dRho_dS*S(j)
    drho_dp(j) = 0.0
  enddo
end subroutine calculate_compress_linear

subroutine calculate_2_densities_linear(T, S, pressure1, pressure2, rho1, rho2,&
                                        start, npts, Rho_T0_S0, dRho_dT, dRho_dS)
  real,    intent(in),  dimension(:) :: T, S
  real,    intent(in)                :: pressure1, pressure2
  real,    intent(out), dimension(:) :: rho1, rho2
  integer, intent(in)                :: start, npts
  real,    intent(in)                :: Rho_T0_S0, dRho_dT, dRho_dS
! *  This subroutine computes the densities of sea water (rho1 and     *
! *  rho2) at two reference pressures (pressure1 and pressure2) from   *
! *  salinity and potential temperature.                               *
! *                                                                    *
! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
! *  (in)      pressure1 - the first pressure in Pa.                   *
! *  (in)      pressure2 -  the second pressure in Pa.                 *
! *  (out)     rho1 - density at pressure1 in kg m-3.                  *
! *  (out)     rho2 - density at pressure2 in kg m-3.                  *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
! *  (in)      Rho_T0_S0 - The density at T=0, S=0, in kg m-3.         *
! *  (in)      dRho_dT - The derivatives of density with temperature   *
! *  (in)      dRho_dS - and salinity, in kg m-3 C-1 and kg m-3 psu-1. *
  integer :: j

  do j=start, start+npts-1
    rho1(j) = Rho_T0_S0 +  dRho_dT*T(j) + dRho_dS*S(j);
    rho2(j) = Rho_T0_S0 +  dRho_dT*T(j) + dRho_dS*S(j);
  enddo
end subroutine calculate_2_densities_linear

subroutine int_density_dz_linear(T, S, z_t, z_b, rho_ref, rho_0_pres, G_e, G, &
                 Rho_T0_S0, dRho_dT, dRho_dS, dpa, intz_dpa, intx_dpa, inty_dpa)
  real, dimension(NIMEM_,NJMEM_),  intent(in)  :: T, S, z_t, z_b
  real,                            intent(in)  :: rho_ref, rho_0_pres, G_e
  type(ocean_grid_type),           intent(in)  :: G
  real,                            intent(in)  :: Rho_T0_S0, dRho_dT, dRho_dS
  real, dimension(NIMEM_,NJMEM_),  intent(out) :: dpa
  real, dimension(NIMEM_,NJMEM_),  optional, intent(out) :: intz_dpa
  real, dimension(NIMEMB_,NJMEM_), optional, intent(out) :: intx_dpa
  real, dimension(NIMEM_,NJMEMB_), optional, intent(out) :: inty_dpa
!   This subroutine calculates analytical and nearly-analytical integrals of
! pressure anomalies across layers, which are required for calculating the
! finite-volume form pressure accelerations in a Boussinesq model.
!
! Arguments: T - potential temperature relative to the surface in C. 
!  (in)      S - salinity in PSU.                                    
!  (in)      z_t - height at the top of the layer in m.           
!  (in)      z_b - height at the top of the layer in m.           
!  (in)      rho_ref - A mean density, in kg m-3, that is subtracted out to reduce
!                    the magnitude of each of the integrals.
!  (in)      rho_0_pres - A density, in kg m-3, that is used to calculate the
!                    pressure (as p~=-z*rho_0_pres*G_e) used in the equation of
!                    state. rho_0_pres is not used here.
!  (in)      G_e - The Earth's gravitational acceleration, in m s-2. 
!  (in)      G - The ocean's grid structure.
!  (in)      Rho_T0_S0 - The density at T=0, S=0, in kg m-3.
!  (in)      dRho_dT - The derivative of density with temperature in kg m-3 C-1.
!  (in)      dRho_dS - The derivative of density with salinity, in kg m-3 psu-1.
!  (out)     dpa - The change in the pressure anomaly across the layer, in Pa.                                  
!  (out,opt) intz_dpa - The integral through the thickness of the layer of the
!                       pressure anomaly relative to the anomaly at the top of
!                       the layer, in Pa m. 
!  (out,opt) intx_dpa - The integral in x of the difference between the
!                       pressure anomaly at the top and bottom of the layer
!                       divided by the x grid spacing, in Pa.
!  (out,opt) inty_dpa - The integral in y of the difference between the
!                       pressure anomaly at the top and bottom of the layer
!                       divided by the y grid spacing, in Pa.
  real :: rho_anom      ! The density anomaly from rho_ref, in kg m-3.
  real :: raL, raR      ! rho_anom to the left and right, in kg m-3.
  real :: dz, dzL, dzR  ! Layer thicknesses in m.
  real :: C1_6
  integer :: Isq, Ieq, Jsq, Jeq, i, j

  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq
  C1_6 = 1.0 / 6.0

  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    dz = z_t(i,j) - z_b(i,j)
    rho_anom = (Rho_T0_S0 - rho_ref) + dRho_dT*T(i,j) + dRho_dS*S(i,j)
    dpa(i,j) = G_e*rho_anom*dz
    if (present(intz_dpa)) intz_dpa(i,j) = 0.5*G_e*rho_anom*dz**2
  enddo ; enddo

  if (present(intx_dpa)) then ; do j=G%jsc,G%jec ; do I=Isq,Ieq
    dzL = z_t(i,j) - z_b(i,j) ; dzR = z_t(i+1,j) - z_b(i+1,j)
    raL = (Rho_T0_S0 - rho_ref) + (dRho_dT*T(i,j) + dRho_dS*S(i,j))
    raR = (Rho_T0_S0 - rho_ref) + (dRho_dT*T(i+1,j) + dRho_dS*S(i+1,j))

    intx_dpa(i,j) = G_e*C1_6 * (dzL*(2.0*raL + raR) + dzR*(2.0*raR + raL))
  enddo ; enddo ; endif

  if (present(inty_dpa)) then ; do J=Jsq,Jeq ; do i=G%isc,G%iec
    dzL = z_t(i,j) - z_b(i,j) ; dzR = z_t(i,j+1) - z_b(i,j+1)
    raL = (Rho_T0_S0 - rho_ref) + (dRho_dT*T(i,j) + dRho_dS*S(i,j))
    raR = (Rho_T0_S0 - rho_ref) + (dRho_dT*T(i,j+1) + dRho_dS*S(i,j+1))

    inty_dpa(i,j) = G_e*C1_6 * (dzL*(2.0*raL + raR) + dzR*(2.0*raR + raL))
  enddo ; enddo ; endif
end subroutine int_density_dz_linear

subroutine int_spec_vol_dp_linear(T, S, p_t, p_b, alpha_ref, G, Rho_T0_S0, &
               dRho_dT, dRho_dS, dza, intp_dza, intx_dza, inty_dza, halo_size)
  real, dimension(NIMEM_,NJMEM_),  intent(in)  :: T, S, p_t, p_b
  real,                            intent(in)  :: alpha_ref
  type(ocean_grid_type),           intent(in)  :: G
  real,    intent(in)  :: Rho_T0_S0, dRho_dT, dRho_dS
  real, dimension(NIMEM_,NJMEM_),  intent(out) :: dza
  real, dimension(NIMEM_,NJMEM_),  optional, intent(out) :: intp_dza
  real, dimension(NIMEMB_,NJMEM_), optional, intent(out) :: intx_dza
  real, dimension(NIMEM_,NJMEMB_), optional, intent(out) :: inty_dza
  integer,                         optional, intent(in)  :: halo_size
!   This subroutine calculates analytical and nearly-analytical integrals in
! pressure across layers of geopotential anomalies, which are required for
! calculating the finite-volume form pressure accelerations in a non-Boussinesq
! model.  Specific volume is assumed to vary linearly between adjacent points.
!
! Arguments: T - potential temperature relative to the surface in C. 
!  (in)      S - salinity in PSU.                                    
!  (in)      p_t - pressure at the top of the layer in Pa.           
!  (in)      p_b - pressure at the top of the layer in Pa.           
!  (in)      alpha_ref - A mean specific volume that is subtracted out to reduce
!                        the magnitude of each of the integrals, m3 kg-1.
!                        The calculation is mathematically identical with
!                        different values of alpha_ref, but this reduces the
!                        effects of roundoff.       
!  (in)      G - The ocean's grid structure.
!  (in)      Rho_T0_S0 - The density at T=0, S=0, in kg m-3.
!  (in)      dRho_dT - The derivative of density with temperature in kg m-3 C-1.
!  (in)      dRho_dS - The derivative of density with salinity, in kg m-3 psu-1.
!  (out)     dza - The change in the geopotential anomaly across the layer, 
!                  in m2 s-2.                                  
!  (out,opt) intp_dza - The integral in pressure through the layer of the 
!                       geopotential anomaly relative to the anomaly at the
!                       bottom of the layer, in Pa m2 s-2. 
!  (out,opt) intx_dza - The integral in x of the difference between the
!                       geopotential anomaly at the top and bottom of the layer
!                       divided by the x grid spacing, in m2 s-2.
!  (out,opt) inty_dza - The integral in y of the difference between the
!                       geopotential anomaly at the top and bottom of the layer
!                       divided by the y grid spacing, in m2 s-2.
  real :: dRho_TS       ! The density anomaly due to T and S, in kg m-3.
  real :: alpha_anom    ! The specific volume anomaly from 1/rho_ref, in m3 kg-1.
  real :: aaL, aaR      ! rho_anom to the left and right, in kg m-3.
  real :: dp, dpL, dpR  ! Layer pressure thicknesses in Pa.
  real :: C1_6
  integer :: Isq, Ieq, Jsq, Jeq, ish, ieh, jsh, jeh, i, j, halo

  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq
  halo = 0 ; if (present(halo_size)) halo = MAX(halo_size,0)
  ish = G%isc-halo ; ieh = G%iec+halo ; jsh = G%jsc-halo ; jeh = G%jec+halo
  if (present(intx_dza)) then ; ish = MIN(Isq,ish) ; ieh = MAX(Ieq+1,ieh); endif
  if (present(inty_dza)) then ; jsh = MIN(Jsq,jsh) ; jeh = MAX(Jeq+1,jeh); endif
  C1_6 = 1.0 / 6.0

  do j=jsh,jeh ; do i=ish,ieh
    dp = p_b(i,j) - p_t(i,j)
    dRho_TS = dRho_dT*T(i,j) + dRho_dS*S(i,j)
    ! alpha_anom = 1.0/(Rho_T0_S0  + dRho_TS)) - alpha_ref
    alpha_anom = ((1.0-Rho_T0_S0*alpha_ref) - dRho_TS*alpha_ref) / (Rho_T0_S0 + dRho_TS)
    dza(i,j) = alpha_anom*dp
    if (present(intp_dza)) intp_dza(i,j) = 0.5*alpha_anom*dp**2
  enddo ; enddo

  if (present(intx_dza)) then ; do j=G%jsc,G%jec ; do I=Isq,Ieq
    dpL = p_b(i,j) - p_t(i,j) ; dpR = p_b(i+1,j) - p_t(i+1,j)
    dRho_TS = dRho_dT*T(i,j) + dRho_dS*S(i,j)
    aaL = ((1.0 - Rho_T0_S0*alpha_ref) - dRho_TS*alpha_ref) / (Rho_T0_S0 + dRho_TS)
    dRho_TS = dRho_dT*T(i+1,j) + dRho_dS*S(i+1,j)
    aaR = ((1.0 - Rho_T0_S0*alpha_ref) - dRho_TS*alpha_ref) / (Rho_T0_S0 + dRho_TS)

    intx_dza(i,j) = C1_6 * (2.0*(dpL*aaL + dpR*aaR) + (dpL*aaR + dpR*aaL))
  enddo ; enddo ; endif

  if (present(inty_dza)) then ; do J=Jsq,Jeq ; do i=G%isc,G%iec
    dpL = p_b(i,j) - p_t(i,j) ; dpR = p_b(i,j+1) - p_t(i,j+1)
    dRho_TS = dRho_dT*T(i,j) + dRho_dS*S(i,j)
    aaL = ((1.0 - Rho_T0_S0*alpha_ref) - dRho_TS*alpha_ref) / (Rho_T0_S0 + dRho_TS)
    dRho_TS = dRho_dT*T(i,j+1) + dRho_dS*S(i,j+1)
    aaR = ((1.0 - Rho_T0_S0*alpha_ref) - dRho_TS*alpha_ref) / (Rho_T0_S0 + dRho_TS)

    inty_dza(i,j) = C1_6 * (2.0*(dpL*aaL + dpR*aaR) + (dpL*aaR + dpR*aaL))
  enddo ; enddo ; endif
end subroutine int_spec_vol_dp_linear

end module MOM_EOS_linear
