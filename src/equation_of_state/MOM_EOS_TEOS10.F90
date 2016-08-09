module MOM_EOS_TEOS10
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
!*  The subroutines in this file implement the equation of state for   *
!*  sea water using the formulae given by  Wright, 1997, J. Atmos.     *
!*  Ocean. Tech., 14, 735-740.  Coded by R. Hallberg, 7/00.            *
!***********************************************************************

use MOM_hor_index, only : hor_index_type
use gsw_mod_toolbox, only : gsw_rho, gsw_rho_first_derivatives, gsw_specvol_first_derivatives

implicit none ; private

#include <MOM_memory.h>

public calculate_compress_teos10, calculate_density_teos10
public calculate_density_derivs_teos10, calculate_2_densities_teos10
public calculate_specvol_derivs_teos10
public calculate_density_scalar_teos10, calculate_density_array_teos10

interface calculate_density_teos10
  module procedure calculate_density_scalar_teos10, calculate_density_array_teos10
end interface calculate_density_teos10

!real :: a0, a1, a2, b0, b1, b2, b3, b4, b5, c0, c1, c2, c3, c4, c5
!    One of the two following blocks of values should be commented out.
!  Following are the values for the full range formula.
!
!real, parameter :: a0 = 7.133718e-4, a1 = 2.724670e-7, a2 = -1.646582e-7
!real, parameter :: b0 = 5.613770e8,  b1 = 3.600337e6,  b2 = -3.727194e4
!real, parameter :: b3 = 1.660557e2,  b4 = 6.844158e5,  b5 = -8.389457e3
!real, parameter :: c0 = 1.609893e5,  c1 = 8.427815e2,  c2 = -6.931554
!real, parameter :: c3 = 3.869318e-2, c4 = -1.664201e2, c5 = -2.765195


! Following are the values for the reduced range formula.
real, parameter :: a0 = 7.057924e-4, a1 = 3.480336e-7, a2 = -1.112733e-7
real, parameter :: b0 = 5.790749e8,  b1 = 3.516535e6,  b2 = -4.002714e4
real, parameter :: b3 = 2.084372e2,  b4 = 5.944068e5,  b5 = -9.643486e3
real, parameter :: c0 = 1.704853e5,  c1 = 7.904722e2,  c2 = -7.984422
real, parameter :: c3 = 5.140652e-2, c4 = -2.302158e2, c5 = -3.079464

contains

subroutine calculate_density_scalar_teos10(T, S, pressure, rho)
real,    intent(in)  :: T, S, pressure
real,    intent(out) :: rho
! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     rho - in situ density in kg m-3.                        *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *

! *====================================================================*
! *  This subroutine computes the in situ density of sea water (rho in *
! *  units of kg/m^3) from salinity (S in psu), potential temperature  *
! *  (T in deg C), and pressure in Pa.  It uses the expression from    *
! *  Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.                *
! *  Coded by R. Hallberg, 7/00                                        *
! *====================================================================*

  real :: al0, p0, lambda
  integer :: j
  real, dimension(1) :: T0, S0, pressure0
  real, dimension(1) :: rho0

  T0(1) = T
  S0(1) = S
  pressure0(1) = pressure

  call calculate_density_array_teos10(T0, S0, pressure0, rho0, 1, 1)
  rho = rho0(1)

end subroutine calculate_density_scalar_teos10

subroutine calculate_density_array_teos10(T, S, pressure, rho, start, npts)
  real,    intent(in),  dimension(:) :: T, S, pressure
  real,    intent(out), dimension(:) :: rho
  integer, intent(in)                :: start, npts
! * Arguments: T - conservative temperature in C.                      *
! *  (in)      S - absolute salinity in g/kg.                          *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     rho - in situ density in kg m-3.                        *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *

! *====================================================================*
! *  This subroutine computes the in situ density of sea water (rho in *
! *  units of kg/m^3) from absolute salinity (S in g/Kg),              *
! *  conservative temperature (T in deg C), and pressure in Pa.        *
! *  It uses the functions from TEOS10 website                         *
! *====================================================================*
  real :: gsw_rho !external TEOS10 function
  real :: p_dbar 
  integer :: j

  do j=start,start+npts-1
    if(S(j).lt.-1.0e-10) cycle !Can we assume safely that this is a missing value?
    p_dbar = pressure(j)*1.0e-4 !convert pressure to dbar    
    rho(j) = gsw_rho(S(j),T(j),p_dbar)
 enddo
end subroutine calculate_density_array_teos10

subroutine calculate_density_derivs_teos10(T, S, pressure, drho_dT, drho_dS, start, npts)
  real,    intent(in),  dimension(:) ::  T, S, pressure
  real,    intent(out), dimension(:) :: drho_dT, drho_dS
  integer, intent(in)                :: start, npts
! * Arguments: T - conservative temperature in C.                      *
! *  (in)      S - absolute salinity in g/kg.                          *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     drho_dT - the partial derivative of density with        *
! *                      potential temperature, in kg m-3 K-1.         *
! *  (out)     drho_dS - the partial derivative of density with        *
! *                      salinity, in kg m-3 psu-1.                    *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
  real :: p_dbar
  integer :: j

  do j=start,start+npts-1
    if(S(j).lt.-1.0e-10) cycle !Can we assume safely that this is a missing value?
    p_dbar = pressure(j)*1.0e-4 !convert pressure to dbar    
    call gsw_rho_first_derivatives(S(j),T(j), p_dbar, drho_dsa=drho_dS(j), drho_dct=drho_dT(j))
  enddo

end subroutine calculate_density_derivs_teos10

subroutine calculate_specvol_derivs_teos10(T, S, pressure, dSV_dT, dSV_dS, start, npts)
  real,    intent(in),  dimension(:) ::  T, S, pressure
  real,    intent(out), dimension(:) :: dSV_dT, dSV_dS
  integer, intent(in)                :: start, npts
! * Arguments: T - conservative temperature in C.                      *
! *  (in)      S - absolute salinity in g/kg.                          *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     dSV_dT - the partial derivative of specific volume with *
! *                     potential temperature, in m3 kg-1 K-1.         *
! *  (out)     dSV_dS - the partial derivative of specific volume with *
! *                      salinity, in m3 kg-1 / (g/kg).                *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
  real :: p_dbar
  integer :: j

  do j=start,start+npts-1
    if(S(j).lt.-1.0e-10) cycle !Can we assume safely that this is a missing value?
    p_dbar = pressure(j)*1.0e-4 !convert pressure to dbar    
    call gsw_specvol_first_derivatives(S(j),T(j), p_dbar, v_sa=dSV_dS(j), v_ct=dSV_dT(j))
  enddo

end subroutine calculate_specvol_derivs_teos10

subroutine calculate_compress_teos10(T, S, pressure, rho, drho_dp, start, npts)
  real,    intent(in),  dimension(:) :: T, S, pressure
  real,    intent(out), dimension(:) :: rho, drho_dp
  integer, intent(in)                :: start, npts
! * Arguments: T - conservative temperature in C.                      *
! *  (in)      S - absolute salinity in g/kg.                          *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     rho - in situ density in kg m-3.                        *
! *  (out)     drho_dp - the partial derivative of density with        *
! *                      pressure (also the inverse of the square of   *
! *                      sound speed) in s2 m-2.                       *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
! *====================================================================*
! *  This subroutine computes the in situ density of sea water (rho in *
! *  units of kg/m^3) and the compressibility (drho/dp = C_sound^-2)   *
! *  (drho_dp in units of s2 m-2) from salinity (sal in psu), potential*
! *  temperature (T in deg C), and pressure in Pa.  It uses the        *
! *  subroutines from TEOS10 website                                   *
! *====================================================================*
  real :: gsw_rho !external TEOS10 function
  real :: p_dbar 
  integer :: j

  do j=start,start+npts-1
    if(S(j).lt.-1.0e-10) cycle !Can we assume safely that this is a missing value?
    p_dbar = pressure(j)*1.0e-4 !convert pressure to dbar    
    rho(j) = gsw_rho(S(j),T(j),p_dbar)
    call gsw_rho_first_derivatives(S(j),T(j), p_dbar, drho_dp=drho_dp(j))
 enddo
end subroutine calculate_compress_teos10

subroutine calculate_2_densities_teos10( T, S, pressure1, pressure2, rho1, rho2, start, npts)
  real,    intent(in),  dimension(:) :: T, S
  real,    intent(in)                :: pressure1, pressure2
  real,    intent(out), dimension(:) :: rho1, rho2
  integer, intent(in)                :: start, npts
! * Arguments: T - conservative temperature in C.                      *
! *  (in)      S - absolute salinity in g/kg.                          *
! *  (in)      pressure1 - the first pressure in Pa.                   *
! *  (in)      pressure2 -  the second pressure in Pa.                 *
! *  (out)     rho1 - density at pressure1 in kg m-3.                  *
! *  (out)     rho2 - density at pressure2 in kg m-3.                  *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *

  real :: gsw_rho !external TEOS10 function
  real :: p1_dbar,p2_dbar 
  integer :: j

  p1_dbar = pressure1*1.0e-4 !convert pressure to dbar    
  p2_dbar = pressure2*1.0e-4 !convert pressure to dbar    

  do j=start,start+npts-1
    if(S(j).lt.-1.0e-10) cycle !Can we assume safely that this is a missing value?
    rho1(j) = gsw_rho(S(j),T(j),p1_dbar)
    rho2(j) = gsw_rho(S(j),T(j),p2_dbar)
 enddo
end subroutine calculate_2_densities_teos10


end module MOM_EOS_TEOS10
