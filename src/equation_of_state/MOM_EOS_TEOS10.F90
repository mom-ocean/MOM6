
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
!*  sea water using the TEOS10 functions                               *
!***********************************************************************

use gsw_mod_toolbox, only : gsw_sp_from_sr, gsw_pt_from_ct
use gsw_mod_toolbox, only : gsw_rho, gsw_rho_first_derivatives, gsw_specvol_first_derivatives
!use gsw_mod_toolbox, only : gsw_sr_from_sp, gsw_ct_from_pt

implicit none ; private

public calculate_compress_teos10, calculate_density_teos10
public calculate_density_derivs_teos10, calculate_specvol_derivs_teos10
public calculate_density_scalar_teos10, calculate_density_array_teos10
public gsw_sp_from_sr, gsw_pt_from_ct

interface calculate_density_teos10
  module procedure calculate_density_scalar_teos10, calculate_density_array_teos10
end interface calculate_density_teos10

real, parameter :: Pa2db  = 1.e-4  ! The conversion factor from Pa to dbar.

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
! *  TEOS10 website.                                                   *
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
  real :: zs,zt,zp
  integer :: j

  do j=start,start+npts-1
    !Conversions
    zs = S(j) !gsw_sr_from_sp(S(j))       !Convert practical salinity to absolute salinity
    zt = T(j) !gsw_ct_from_pt(S(j),T(j))  !Convert potantial temp to conservative temp
    zp = pressure(j)* Pa2db         !Convert pressure from Pascal to decibar

    if(S(j).lt.-1.0e-10) cycle !Can we assume safely that this is a missing value?
    rho(j) = gsw_rho(zs,zt,zp)
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
  real :: zs,zt,zp
  integer :: j

  do j=start,start+npts-1
    !Conversions
    zs = S(j) !gsw_sr_from_sp(S(j))       !Convert practical salinity to absolute salinity
    zt = T(j) !gsw_ct_from_pt(S(j),T(j))  !Convert potantial temp to conservative temp
    zp = pressure(j)* Pa2db         !Convert pressure from Pascal to decibar
    if(S(j).lt.-1.0e-10) cycle !Can we assume safely that this is a missing value?
    call gsw_rho_first_derivatives(zs, zt, zp, drho_dsa=drho_dS(j), drho_dct=drho_dT(j))
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
  real :: zs, zt, zp
  integer :: j

  do j=start,start+npts-1
    !Conversions
    zs = S(j) !gsw_sr_from_sp(S(j))       !Convert practical salinity to absolute salinity
    zt = T(j) !gsw_ct_from_pt(S(j),T(j))  !Convert potantial temp to conservative temp
    zp = pressure(j)* Pa2db         !Convert pressure from Pascal to decibar
    if(S(j).lt.-1.0e-10) cycle !Can we assume safely that this is a missing value?
    call gsw_specvol_first_derivatives(zs,zt,zp, v_sa=dSV_dS(j), v_ct=dSV_dT(j))
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
  real :: zs,zt,zp
  integer :: j

  do j=start,start+npts-1
    !Conversions
    zs = S(j) !gsw_sr_from_sp(S(j))       !Convert practical salinity to absolute salinity
    zt = T(j) !gsw_ct_from_pt(S(j),T(j))  !Convert potantial temp to conservative temp
    zp = pressure(j)* Pa2db         !Convert pressure from Pascal to decibar
    if(S(j).lt.-1.0e-10) cycle !Can we assume safely that this is a missing value?
    rho(j) = gsw_rho(zs,zt,zp)
    call gsw_rho_first_derivatives(zs,zt,zp, drho_dp=drho_dp(j))
 enddo
end subroutine calculate_compress_teos10

end module MOM_EOS_TEOS10
