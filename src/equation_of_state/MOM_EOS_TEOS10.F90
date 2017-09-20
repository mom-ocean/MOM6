module MOM_EOS_TEOS10

! This file is part of MOM6. See LICENSE.md for the license.

!***********************************************************************
!*  The subroutines in this file implement the equation of state for   *
!*  sea water using the TEOS10 functions                               *
!***********************************************************************

use gsw_mod_toolbox, only : gsw_sp_from_sr, gsw_pt_from_ct
use gsw_mod_toolbox, only : gsw_rho, gsw_rho_first_derivatives, gsw_specvol_first_derivatives
use gsw_mod_toolbox, only : gsw_rho_second_derivatives
!use gsw_mod_toolbox, only : gsw_sr_from_sp, gsw_ct_from_pt

implicit none ; private

public calculate_compress_teos10, calculate_density_teos10
public calculate_density_derivs_teos10
public calculate_specvol_derivs_teos10
public calculate_density_second_derivs_teos10
public gsw_sp_from_sr, gsw_pt_from_ct

interface calculate_density_teos10
  module procedure calculate_density_scalar_teos10, calculate_density_array_teos10
end interface calculate_density_teos10

interface calculate_density_derivs_teos10
  module procedure calculate_density_derivs_scalar_teos10, calculate_density_derivs_array_teos10
end interface calculate_density_derivs_teos10

interface calculate_density_second_derivs_teos10
  module procedure calculate_density_second_derivs_scalar_teos10, calculate_density_second_derivs_array_teos10
end interface calculate_density_second_derivs_teos10

real, parameter :: Pa2db  = 1.e-4  ! The conversion factor from Pa to dbar.

contains

!> This subroutine computes the in situ density of sea water (rho in
!! units of kg/m^3) from salinity (S in psu), potential temperature
!! (T in deg C), and pressure in Pa.  It uses the expression from
!! TEOS10 website.
subroutine calculate_density_scalar_teos10(T, S, pressure, rho)
real,    intent(in)  :: T        !< Conservative temperature in C.
real,    intent(in)  :: S        !< Absolute salinity in g/kg.
real,    intent(in)  :: pressure !< Pressure in Pa.
real,    intent(out) :: rho      !< In situ density in kg m-3.
! * Arguments: T - conservative temperature in C.                      *
! *  (in)      S - absolute salinity in g/kg.                          *
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

subroutine calculate_density_derivs_array_teos10(T, S, pressure, drho_dT, drho_dS, start, npts)
  real,    intent(in),  dimension(:) :: T        !< Conservative temperature in C.
  real,    intent(in),  dimension(:) :: S        !< Absolute salinity in g/kg.
  real,    intent(in),  dimension(:) :: pressure !< Pressure in Pa.
  real,    intent(out), dimension(:) :: drho_dT  !< The partial derivative of density with potential
                                                 !! temperature, in kg m-3 K-1.
  real,    intent(out), dimension(:) :: drho_dS  !< The partial derivative of density with salinity,
                                                 !! in kg m-3 psu-1.
  integer, intent(in)                :: start    !< The starting point in the arrays.
  integer, intent(in)                :: npts     !< The number of values to calculate.
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

end subroutine calculate_density_derivs_array_teos10

subroutine calculate_density_derivs_scalar_teos10(T, S, pressure, drho_dT, drho_dS)
  real,    intent(in)  ::  T, S, pressure
  real,    intent(out) :: drho_dT, drho_dS
  ! Local variables
  real :: zs,zt,zp
  !Conversions
  zs = S !gsw_sr_from_sp(S)       !Convert practical salinity to absolute salinity
  zt = T !gsw_ct_from_pt(S,T)  !Convert potantial temp to conservative temp
  zp = pressure* Pa2db         !Convert pressure from Pascal to decibar
  if(S.lt.-1.0e-10) return !Can we assume safely that this is a missing value?
  call gsw_rho_first_derivatives(zs, zt, zp, drho_dsa=drho_dS, drho_dct=drho_dT)
end subroutine calculate_density_derivs_scalar_teos10

subroutine calculate_specvol_derivs_teos10(T, S, pressure, dSV_dT, dSV_dS, start, npts)
  real,    intent(in),  dimension(:) :: T        !< Conservative temperature in C.
  real,    intent(in),  dimension(:) :: S        !< Absolute salinity in g/kg.
  real,    intent(in),  dimension(:) :: pressure !< Pressure in Pa.
  real,    intent(out), dimension(:) :: dSV_dT   !< The partial derivative of specific volume with
                                                 !! potential temperature, in m3 kg-1 K-1.
  real,    intent(out), dimension(:) :: dSV_dS   !< The partial derivative of specific volume with
                                                 !! salinity, in m3 kg-1 / (g/kg).
  integer, intent(in)                :: start    !< The starting point in the arrays.
  integer, intent(in)                :: npts     !< The number of values to calculate.
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

!> Calculate the 5 second derivatives of the equation of state for scalar inputs
subroutine calculate_density_second_derivs_scalar_teos10(T, S, pressure, drho_dS_dS, drho_dS_dT, drho_dT_dT, drho_dS_dP, &
                                                         drho_dT_dP)
  real, intent(in)     :: T, S, pressure
  real, intent(out)    :: drho_dS_dS !< Partial derivative of beta with respect to S
  real, intent(out)    :: drho_dS_dT !< Partial derivative of beta with resepct to T
  real, intent(out)    :: drho_dT_dT !< Partial derivative of alpha with respect to T
  real, intent(out)    :: drho_dS_dP !< Partial derivative of beta with respect to pressure
  real, intent(out)    :: drho_dT_dP !< Partial derivative of alpha with respect to pressure
! * Arguments: T - conservative temperature in C.                      *
! *  (in)      S - absolute salinity in g/kg.                          *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     drho_dT - the partial derivative of density with        *
! *                      potential temperature, in kg m-3 K-1.         *
! *  (out)     drho_dS - the partial derivative of density with        *
! *                      salinity, in kg m-3 psu-1.                    *
  real :: zs,zt,zp

  !Conversions
  zs = S !gsw_sr_from_sp(S)       !Convert practical salinity to absolute salinity
  zt = T !gsw_ct_from_pt(S,T)  !Convert potantial temp to conservative temp
  zp = pressure* Pa2db         !Convert pressure from Pascal to decibar
  if(S.lt.-1.0e-10) return !Can we assume safely that this is a missing value?
  call gsw_rho_second_derivatives(zs, zt, zp, rho_sa_sa=drho_dS_dS, rho_sa_ct=drho_dS_dT, &
                                     rho_ct_ct=drho_dT_dT, rho_sa_p=drho_dS_dP, rho_ct_p=drho_dT_dP)

end subroutine calculate_density_second_derivs_scalar_teos10

!> Calculate the 5 second derivatives of the equation of state for scalar inputs
subroutine calculate_density_second_derivs_array_teos10(T, S, pressure, drho_dS_dS, drho_dS_dT, drho_dT_dT, drho_dS_dP, &
                                                        drho_dT_dP, start, npts)
  real, dimension(:), intent(in)     :: T, S, pressure
  real, dimension(:), intent(out)    :: drho_dS_dS !< Partial derivative of beta with respect to S
  real, dimension(:), intent(out)    :: drho_dS_dT !< Partial derivative of beta with resepct to T
  real, dimension(:), intent(out)    :: drho_dT_dT !< Partial derivative of alpha with respect to T
  real, dimension(:), intent(out)    :: drho_dS_dP !< Partial derivative of beta with respect to pressure
  real, dimension(:), intent(out)    :: drho_dT_dP !< Partial derivative of alpha with respect to pressure
  integer, intent(in)  :: start    !< The starting point in the arrays.
  integer, intent(in)  :: npts     !< The number of values to calculate.
! * Arguments: T - conservative temperature in C.                      *
! *  (in)      S - absolute salinity in g/kg.                          *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     drho_dT - the partial derivative of density with        *
! *                      potential temperature, in kg m-3 K-1.         *
! *  (out)     drho_dS - the partial derivative of density with        *
! *                      salinity, in kg m-3 psu-1.                    *
  real :: zs,zt,zp
  integer :: j
  do j=start,start+npts-1
    !Conversions
    zs = S(j) !gsw_sr_from_sp(S)       !Convert practical salinity to absolute salinity
    zt = T(j) !gsw_ct_from_pt(S,T)  !Convert potantial temp to conservative temp
    zp = pressure(j)* Pa2db         !Convert pressure from Pascal to decibar
    if(zs .lt. -1.0e-10) return !Can we assume safely that this is a missing value?
    call gsw_rho_second_derivatives(zs, zt, zp, rho_sa_sa=drho_dS_dS(j), rho_sa_ct=drho_dS_dT(j), &
                                    rho_ct_ct=drho_dT_dT(j), rho_sa_p=drho_dS_dP(j), rho_ct_p=drho_dT_dP(j))
  enddo

end subroutine calculate_density_second_derivs_array_teos10

!> This subroutine computes the in situ density of sea water (rho in *
!! units of kg/m^3) and the compressibility (drho/dp = C_sound^-2)   *
!! (drho_dp in units of s2 m-2) from salinity (sal in psu), potential*
!! temperature (T in deg C), and pressure in Pa.  It uses the        *
!! subroutines from TEOS10 website
subroutine calculate_compress_teos10(T, S, pressure, rho, drho_dp, start, npts)
  real,    intent(in),  dimension(:) :: T        !< Conservative temperature in C.
  real,    intent(in),  dimension(:) :: S        !< Absolute salinity in g/kg.
  real,    intent(in),  dimension(:) :: pressure !< Pressure in Pa.
  real,    intent(out), dimension(:) :: rho      !< In situ density in kg m-3.
  real,    intent(out), dimension(:) :: drho_dp  !< The partial derivative of density with pressure
                                                 !! (also the inverse of the square of sound speed)
                                                 !! in s2 m-2.
  integer, intent(in)                :: start    !< The starting point in the arrays.
  integer, intent(in)                :: npts     !< The number of values to calculate.
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
