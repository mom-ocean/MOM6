module MOM_EOS_Wright

! This file is part of MOM6. See LICENSE.md for the license.

!***********************************************************************
!*  The subroutines in this file implement the equation of state for   *
!*  sea water using the formulae given by  Wright, 1997, J. Atmos.     *
!*  Ocean. Tech., 14, 735-740.  Coded by R. Hallberg, 7/00.            *
!***********************************************************************

use MOM_hor_index, only : hor_index_type

implicit none ; private

#include <MOM_memory.h>

public calculate_compress_wright, calculate_density_wright
public calculate_density_derivs_wright, calculate_specvol_derivs_wright
public calculate_density_second_derivs_wright
public int_density_dz_wright, int_spec_vol_dp_wright

interface calculate_density_wright
  module procedure calculate_density_scalar_wright, calculate_density_array_wright
end interface calculate_density_wright

interface calculate_density_derivs_wright
  module procedure calculate_density_derivs_scalar_wright, calculate_density_derivs_array_wright
end interface

interface calculate_density_second_derivs_wright
  module procedure calculate_density_second_derivs_scalar_wright, calculate_density_second_derivs_array_wright
end interface

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

!> This subroutine computes the in situ density of sea water (rho in
!! units of kg/m^3) from salinity (S in psu), potential temperature
!! (T in deg C), and pressure in Pa.  It uses the expression from
!! Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.
!! Coded by R. Hallberg, 7/00
subroutine calculate_density_scalar_wright(T, S, pressure, rho)
real,    intent(in)  :: T        !< Potential temperature relative to the surface in C.
real,    intent(in)  :: S        !< Salinity in PSU.
real,    intent(in)  :: pressure !< Pressure in Pa.
real,    intent(out) :: rho      !< In situ density in kg m-3.

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

  call calculate_density_array_wright(T0, S0, pressure0, rho0, 1, 1)
  rho = rho0(1)

end subroutine calculate_density_scalar_wright

!> This subroutine computes the in situ density of sea water (rho in
!! units of kg/m^3) from salinity (S in psu), potential temperature
!! (T in deg C), and pressure in Pa.  It uses the expression from
!! Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.
!! Coded by R. Hallberg, 7/00
subroutine calculate_density_array_wright(T, S, pressure, rho, start, npts)
  real,    intent(in),  dimension(:) :: T        !< potential temperature relative to the surface
                                                 !! in C.
  real,    intent(in),  dimension(:) :: S        !< salinity in PSU.
  real,    intent(in),  dimension(:) :: pressure !< pressure in Pa.
  real,    intent(out), dimension(:) :: rho      !< in situ density in kg m-3.
  integer, intent(in)                :: start    !< the starting point in the arrays.
  integer, intent(in)                :: npts     !< the number of values to calculate.

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

  do j=start,start+npts-1
    al0 = (a0 + a1*T(j)) +a2*S(j)
    p0 = (b0 + b4*S(j)) + T(j) * (b1 + T(j)*((b2 + b3*T(j))) + b5*S(j))
    lambda = (c0 +c4*S(j)) + T(j) * (c1 + T(j)*((c2 + c3*T(j))) + c5*S(j))

    rho(j) = (pressure(j) + p0) / (lambda + al0*(pressure(j) + p0))
  enddo
end subroutine calculate_density_array_wright

!> For a given thermodynamic state, return the thermal/haline expansion coefficients
subroutine calculate_density_derivs_array_wright(T, S, pressure, drho_dT, drho_dS, start, npts)
  real,    intent(in),  dimension(:) :: T        !< Potential temperature relative to the surface
                                                 !! in C.
  real,    intent(in),  dimension(:) :: S        !< Salinity in PSU.
  real,    intent(in),  dimension(:) :: pressure !< Pressure in Pa.
  real,    intent(out), dimension(:) :: drho_dT  !< The partial derivative of density with potential
                                                 !! temperature, in kg m-3 K-1.
  real,    intent(out), dimension(:) :: drho_dS  !< The partial derivative of density with salinity,
                                                 !! in kg m-3 psu-1.
  integer, intent(in)                :: start    !< The starting point in the arrays.
  integer, intent(in)                :: npts     !< The number of values to calculate.

! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     drho_dT - the partial derivative of density with        *
! *                      potential temperature, in kg m-3 K-1.         *
! *  (out)     drho_dS - the partial derivative of density with        *
! *                      salinity, in kg m-3 psu-1.                    *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
  real :: al0, p0, lambda, I_denom2
  integer :: j

  do j=start,start+npts-1
    al0 = (a0 + a1*T(j)) + a2*S(j)
    p0 = (b0 + b4*S(j)) + T(j) * (b1 + T(j)*((b2 + b3*T(j))) + b5*S(j))
    lambda = (c0 +c4*S(j)) + T(j) * (c1 + T(j)*((c2 + c3*T(j))) + c5*S(j))

    I_denom2 = 1.0 / (lambda + al0*(pressure(j) + p0))
    I_denom2 = I_denom2 *I_denom2
    drho_dT(j) = I_denom2 * &
      (lambda* (b1 + T(j)*(2.0*b2 + 3.0*b3*T(j)) + b5*S(j)) - &
       (pressure(j)+p0) * ( (pressure(j)+p0)*a1 + &
        (c1 + T(j)*(c2*2.0 + c3*3.0*T(j)) + c5*S(j)) ))
    drho_dS(j) = I_denom2 * (lambda* (b4 + b5*T(j)) - &
      (pressure(j)+p0) * ( (pressure(j)+p0)*a2 + (c4 + c5*T(j)) ))
  enddo

end subroutine calculate_density_derivs_array_wright

!> The scalar version of calculate_density_derivs which promotes scalar inputs to a 1-element array and then
!! demotes the output back to a scalar
subroutine calculate_density_derivs_scalar_wright(T, S, pressure, drho_dT, drho_dS)
  real,    intent(in) :: T        !< Potential temperature relative to the surface
                                  !! in C.
  real,    intent(in) :: S        !< Salinity in PSU.
  real,    intent(in) :: pressure !< Pressure in Pa.
  real,    intent(out) :: drho_dT  !< The partial derivative of density with potential
                                   !! temperature, in kg m-3 K-1.
  real,    intent(out) :: drho_dS  !< The partial derivative of density with salinity,
                                   !! in kg m-3 psu-1.

  ! Local variables needed to promote the input/output scalars to 1-element arrays
  real, dimension(1) :: T0, S0, P0
  real, dimension(1) :: drdt0, drds0

  T0(1) = T
  S0(1) = S
  P0(1) = pressure
  call calculate_density_derivs_array_wright(T0, S0, P0, drdt0, drds0, 1, 1)
  drho_dT = drdt0(1)
  drho_dS = drds0(1)

end subroutine calculate_density_derivs_scalar_wright

!> Second derivatives of density with respect to temperature, salinity, and pressure
subroutine calculate_density_second_derivs_array_wright(T, S, P, drho_ds_ds, drho_ds_dt, drho_dt_dt, &
                                                         drho_ds_dp, drho_dt_dp, start, npts)
  real, dimension(:), intent(in   ) :: T !< Potential temperature referenced to 0 dbar
  real, dimension(:), intent(in   ) :: S !< Salinity in PSU
  real, dimension(:), intent(in   ) :: P !< Pressure in Pa
  real, dimension(:), intent(  out) :: drho_ds_ds !< Partial derivative of beta with respect to S
  real, dimension(:), intent(  out) :: drho_ds_dt !< Partial derivative of beta with resepct to T
  real, dimension(:), intent(  out) :: drho_dt_dt !< Partial derivative of alpha with respect to T
  real, dimension(:), intent(  out) :: drho_ds_dp !< Partial derivative of beta with respect to pressure
  real, dimension(:), intent(  out) :: drho_dt_dp !< Partial derivative of alpha with respect to pressure
  integer,            intent(in   ) :: start !< Starting index in T,S,P
  integer,            intent(in   ) :: npts  !< Number of points to loop over

  integer :: j
  ! Based on the above expression with common terms factored, there probably exists a more numerically stable
  ! and/or efficient expression
  real :: z0, z1, z2, z3, z4, z5, z6 ,z7, z8, z9, z10, z11, z2_2, z2_3

  do j = start,start+npts-1
    z0 = T(j)*(b1 + b5*S(j) + T(j)*(b2 + b3*T(j)))
    z1 = (b0 + P(j) + b4*S(j) + z0)
    z3 = (b1 + b5*S(j) + T(j)*(2.*b2 + 2.*b3*T(j)))
    z4 = (c0 + c4*S(j) + T(j)*(c1 + c5*S(j) + T(j)*(c2 + c3*T(j))))
    z5 = (b1 + b5*S(j) + T(j)*(b2 + b3*T(j)) + T(j)*(b2 + 2.*b3*T(j)))
    z6 = c1 + c5*S(j) + T(j)*(c2 + c3*T(j)) + T(j)*(c2 + 2.*c3*T(j))
    z7 = (c4 + c5*T(j) + a2*z1)
    z8 = (c1 + c5*S(j) + T(j)*(2.*c2 + 3.*c3*T(j)) + a1*z1)
    z9 = (a0 + a2*S(j) + a1*T(j))
    z10 = (b4 + b5*T(j))
    z11 = (z10*z4 - z1*z7)
    z2 = (c0 + c4*S(j) + T(j)*(c1 + c5*S(j) + T(j)*(c2 + c3*T(j))) + z9*z1)
    z2_2 = z2*z2
    z2_3 = z2_2*z2

    drho_ds_ds(j) = (z10*(c4 + c5*T(j)) - a2*z10*z1 - z10*z7)/z2_2 - (2.*(c4 + c5*T(j) + z9*z10 + a2*z1)*z11)/z2_3
    drho_ds_dt(j) = (z10*z6 - z1*(c5 + a2*z5) + b5*z4 - z5*z7)/z2_2 - (2.*(z6 + z9*z5 + a1*z1)*z11)/z2_3
    drho_dt_dt(j) = (z3*z6 - z1*(2.*c2 + 6.*c3*T(j) + a1*z5) + (2.*b2 + 4.*b3*T(j))*z4 - z5*z8)/z2_2 - &
                    (2.*(z6 + z9*z5 + a1*z1)*(z3*z4 - z1*z8))/z2_3
    drho_ds_dp(j) = (-c4 - c5*T(j) - 2.*a2*z1)/z2_2 - (2.*z9*z11)/z2_3
    drho_dt_dp(j) = (-c1 - c5*S(j) - T(j)*(2.*c2 + 3.*c3*T(j)) - 2.*a1*z1)/z2_2 - (2.*z9*(z3*z4 - z1*z8))/z2_3
  enddo

end subroutine calculate_density_second_derivs_array_wright

!> Second derivatives of density with respect to temperature, salinity, and pressure for scalar inputs. Inputs
!! promoted to 1-element array and output demoted to scalar
subroutine calculate_density_second_derivs_scalar_wright(T, S, P, drho_ds_ds, drho_ds_dt, drho_dt_dt, &
                                                         drho_ds_dp, drho_dt_dp)
  real, intent(in   ) :: T          !< Potential temperature referenced to 0 dbar
  real, intent(in   ) :: S          !< Salinity in PSU
  real, intent(in   ) :: P          !< Pressure in Pa
  real, intent(  out) :: drho_ds_ds !< Partial derivative of beta with respect to S
  real, intent(  out) :: drho_ds_dt !< Partial derivative of beta with resepct to T
  real, intent(  out) :: drho_dt_dt !< Partial derivative of alpha with respect to T
  real, intent(  out) :: drho_ds_dp !< Partial derivative of beta with respect to pressure
  real, intent(  out) :: drho_dt_dp !< Partial derivative of alpha with respect to pressure
  ! Local variables
  real, dimension(1) :: T0, S0, P0
  real, dimension(1) :: drdsds, drdsdt, drdtdt, drdsdp, drdtdp

  T0(1) = T
  S0(1) = S
  P0(1) = P
  call calculate_density_second_derivs_array_wright(T0, S0, P0, drdsds, drdsdt, drdtdt, drdsdp, drdtdp, 1, 1)
  drho_ds_ds = drdsds(1)
  drho_ds_dt = drdsdt(1)
  drho_dt_dt = drdtdt(1)
  drho_ds_dp = drdsdp(1)
  drho_dt_dp = drdtdp(1)

end subroutine calculate_density_second_derivs_scalar_wright

subroutine calculate_specvol_derivs_wright(T, S, pressure, dSV_dT, dSV_dS, start, npts)
  real,    intent(in),  dimension(:) :: T        !< Potential temperature relative to the surface
                                                 !! in C.
  real,    intent(in),  dimension(:) :: S        !< Salinity in g/kg.
  real,    intent(in),  dimension(:) :: pressure !< Pressure in Pa.
  real,    intent(out), dimension(:) :: dSV_dT   !< The partial derivative of specific volume with
                                                 !! potential temperature, in m3 kg-1 K-1.
  real,    intent(out), dimension(:) :: dSV_dS   !< The partial derivative of specific volume with
                                                 !! salinity, in m3 kg-1 / (g/kg).
  integer, intent(in)                :: start    !< The starting point in the arrays.
  integer, intent(in)                :: npts     !< The number of values to calculate.

! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in g/kg.                                   *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     dSV_dT - the partial derivative of specific volume with *
! *                     potential temperature, in m3 kg-1 K-1.         *
! *  (out)     dSV_dS - the partial derivative of specific volume with *
! *                      salinity, in m3 kg-1 / (g/kg).                *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
  real :: al0, p0, lambda, I_denom
  integer :: j

  do j=start,start+npts-1
!    al0 = (a0 + a1*T(j)) + a2*S(j)
    p0 = (b0 + b4*S(j)) + T(j) * (b1 + T(j)*((b2 + b3*T(j))) + b5*S(j))
    lambda = (c0 +c4*S(j)) + T(j) * (c1 + T(j)*((c2 + c3*T(j))) + c5*S(j))

    ! SV = al0 + lambda / (pressure(j) + p0)

    I_denom = 1.0 / (pressure(j) + p0)
    dSV_dT(j) = (a1 + I_denom * (c1 + T(j)*((2.0*c2 + 3.0*c3*T(j))) + c5*S(j))) - &
                (I_denom**2 * lambda) *  (b1 + T(j)*((2.0*b2 + 3.0*b3*T(j))) + b5*S(j))
    dSV_dS(j) = (a2 + I_denom * (c4 + c5*T(j))) - &
                (I_denom**2 * lambda) *  (b4 + b5*T(j))
  enddo

end subroutine calculate_specvol_derivs_wright

!> This subroutine computes the in situ density of sea water (rho in
!! units of kg/m^3) and the compressibility (drho/dp = C_sound^-2)
!! (drho_dp in units of s2 m-2) from salinity (sal in psu), potential
!! temperature (T in deg C), and pressure in Pa.  It uses the expressions
!! from Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.
!! Coded by R. Hallberg, 1/01
subroutine calculate_compress_wright(T, S, pressure, rho, drho_dp, start, npts)
  real,    intent(in),  dimension(:) :: T        !< Potential temperature relative to the surface
                                                 !! in C.
  real,    intent(in),  dimension(:) :: S        !< Salinity in PSU.
  real,    intent(in),  dimension(:) :: pressure !< Pressure in Pa.
  real,    intent(out), dimension(:) :: rho      !< In situ density in kg m-3.
  real,    intent(out), dimension(:) :: drho_dp  !< The partial derivative of density with pressure
                                                 !! (also the inverse of the square of sound speed)
                                                 !! in s2 m-2.
  integer, intent(in)                :: start    !< The starting point in the arrays.
  integer, intent(in)                :: npts     !< The number of values to calculate.

! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
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
! *  temperature (T in deg C), and pressure in Pa.  It uses the expres-*
! *  sions from Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.     *
! *  Coded by R. Hallberg, 1/01                                        *
! *====================================================================*
  real :: al0, p0, lambda, I_denom
  integer :: j

  do j=start,start+npts-1
    al0 = (a0 + a1*T(j)) +a2*S(j)
    p0 = (b0 + b4*S(j)) + T(j) * (b1 + T(j)*((b2 + b3*T(j))) + b5*S(j))
    lambda = (c0 +c4*S(j)) + T(j) * (c1 + T(j)*((c2 + c3*T(j))) + c5*S(j))

    I_denom = 1.0 / (lambda + al0*(pressure(j) + p0))
    rho(j) = (pressure(j) + p0) * I_denom
    drho_dp(j) = lambda * I_denom * I_denom
  enddo
end subroutine calculate_compress_wright

!> This subroutine calculates analytical and nearly-analytical integrals of
!! pressure anomalies across layers, which are required for calculating the
!! finite-volume form pressure accelerations in a Boussinesq model.
subroutine int_density_dz_wright(T, S, z_t, z_b, rho_ref, rho_0, G_e, HII, HIO, &
                                 dpa, intz_dpa, intx_dpa, inty_dpa)
  type(hor_index_type), intent(in)  :: HII, HIO
  real, dimension(HII%isd:HII%ied,HII%jsd:HII%jed), &
                        intent(in)  :: T        !< Potential temperature relative to the surface
                                                !! in C.
  real, dimension(HII%isd:HII%ied,HII%jsd:HII%jed), &
                        intent(in)  :: S        !< Salinity in PSU.
  real, dimension(HII%isd:HII%ied,HII%jsd:HII%jed), &
                        intent(in)  :: z_t      !< Height at the top of the layer in m.
  real, dimension(HII%isd:HII%ied,HII%jsd:HII%jed), &
                        intent(in)  :: z_b      !< Height at the top of the layer in m.
  real,                 intent(in)  :: rho_ref  !< A mean density, in kg m-3, that is subtracted out
                                                !! to reduce the magnitude of each of the integrals.
                                                !! (The pressure is calucated as p~=-z*rho_0*G_e.)
  real,                 intent(in)  :: rho_0    !< Density, in kg m-3, that is used to calculate the
                                                !! pressure (as p~=-z*rho_0*G_e) used in the
                                                !! equation of state.
  real,                 intent(in)  :: G_e      !< The Earth's gravitational acceleration, in m s-2.
  real, dimension(HIO%isd:HIO%ied,HIO%jsd:HIO%jed), &
                        intent(out) :: dpa      !< The change in the pressure anomaly across the
                                                !! layer, in Pa.
  real, dimension(HIO%isd:HIO%ied,HIO%jsd:HIO%jed), &
              optional, intent(out) :: intz_dpa !< The integral through the thickness of the layer
                                                !! of the pressure anomaly relative to the anomaly
                                                !! at the top of the layer, in Pa m.
  real, dimension(HIO%IsdB:HIO%IedB,HIO%jsd:HIO%jed), &
              optional, intent(out) :: intx_dpa !< The integral in x of the difference between the
                                                !! pressure anomaly at the top and bottom of the
                                                !! layer divided by the x grid spacing, in Pa.
  real, dimension(HIO%isd:HIO%ied,HIO%JsdB:HIO%JedB), &
              optional, intent(out) :: inty_dpa !< The integral in y of the difference between the
                                                !! pressure anomaly at the top and bottom of the
                                                !! layer divided by the y grid spacing, in Pa.

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
!                    (The pressure is calucated as p~=-z*rho_0*G_e.)
!  (in)      rho_0 - A density, in kg m-3, that is used to calculate the pressure
!                    (as p~=-z*rho_0*G_e) used in the equation of state.
!  (in)      G_e - The Earth's gravitational acceleration, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (out)     dpa - The change in the pressure anomaly across the layer,
!                  in Pa.
!  (out,opt) intz_dpa - The integral through the thickness of the layer of the
!                       pressure anomaly relative to the anomaly at the top of
!                       the layer, in Pa m.
!  (out,opt) intx_dpa - The integral in x of the difference between the
!                       pressure anomaly at the top and bottom of the layer
!                       divided by the x grid spacing, in Pa.
!  (out,opt) inty_dpa - The integral in y of the difference between the
!                       pressure anomaly at the top and bottom of the layer
!                       divided by the y grid spacing, in Pa.

  real, dimension(HII%isd:HII%ied,HII%jsd:HII%jed) :: al0_2d, p0_2d, lambda_2d
  real :: al0, p0, lambda
  real :: eps, eps2, rho_anom, rem
  real :: w_left, w_right, intz(5)
  real, parameter :: C1_3 = 1.0/3.0, C1_7 = 1.0/7.0    ! Rational constants.
  real, parameter :: C1_9 = 1.0/9.0, C1_90 = 1.0/90.0  ! Rational constants.
  real :: GxRho, I_Rho
  real :: dz, p_ave, I_al0, I_Lzz
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, i, j, ioff, joff, m

  ioff = HIO%idg_offset - HII%idg_offset
  joff = HIO%jdg_offset - HII%jdg_offset

  ! These array bounds work for the indexing convention of the input arrays, but
  ! on the computational domain defined for the output arrays.
  Isq = HIO%IscB + ioff ; Ieq = HIO%IecB + ioff
  Jsq = HIO%JscB + joff ; Jeq = HIO%JecB + joff
  is = HIO%isc + ioff ; ie = HIO%iec + ioff
  js = HIO%jsc + joff ; je = HIO%jec + joff

  GxRho = G_e * rho_0
  I_Rho = 1.0 / rho_0

  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    al0_2d(i,j) = (a0 + a1*T(i,j)) + a2*S(i,j)
    p0_2d(i,j) = (b0 + b4*S(i,j)) + T(i,j) * (b1 + T(i,j)*((b2 + b3*T(i,j))) + b5*S(i,j))
    lambda_2d(i,j) = (c0 +c4*S(i,j)) + T(i,j) * (c1 + T(i,j)*((c2 + c3*T(i,j))) + c5*S(i,j))

    al0 = al0_2d(i,j) ; p0 = p0_2d(i,j) ; lambda = lambda_2d(i,j)

    dz = z_t(i,j) - z_b(i,j)
    p_ave = -0.5*GxRho*(z_t(i,j)+z_b(i,j))

    I_al0 = 1.0 / al0
    I_Lzz = 1.0 / (p0 + (lambda * I_al0) + p_ave)
    eps = 0.5*GxRho*dz*I_Lzz ; eps2 = eps*eps

!     rho(j) = (pressure(j) + p0) / (lambda + al0*(pressure(j) + p0))

    rho_anom = (p0 + p_ave)*(I_Lzz*I_al0) - rho_ref
    rem = I_Rho * (lambda * I_al0**2) * eps2 * &
          (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2)))
    dpa(i-ioff,j-joff) = G_e*rho_anom*dz - 2.0*eps*rem
    if (present(intz_dpa)) &
      intz_dpa(i-ioff,j-joff) = 0.5*G_e*rho_anom*dz**2 - dz*(1.0+eps)*rem
  enddo ; enddo

  if (present(intx_dpa)) then ; do j=js,je ; do I=Isq,Ieq
    intz(1) = dpa(i-ioff,j-joff) ; intz(5) = dpa(i+1-ioff,j-joff)
    do m=2,4
      w_left = 0.25*real(m-1) ; w_right = 1.0-w_left
      al0 = w_left*al0_2d(i,j) + w_right*al0_2d(i+1,j)
      p0 = w_left*p0_2d(i,j) + w_right*p0_2d(i+1,j)
      lambda = w_left*lambda_2d(i,j) + w_right*lambda_2d(i+1,j)

      dz = w_left*(z_t(i,j) - z_b(i,j)) + w_right*(z_t(i+1,j) - z_b(i+1,j))
      p_ave = -0.5*GxRho*(w_left*(z_t(i,j)+z_b(i,j)) + &
                          w_right*(z_t(i+1,j)+z_b(i+1,j)))

      I_al0 = 1.0 / al0
      I_Lzz = 1.0 / (p0 + (lambda * I_al0) + p_ave)
      eps = 0.5*GxRho*dz*I_Lzz ; eps2 = eps*eps

      intz(m) = G_e*dz*((p0 + p_ave)*(I_Lzz*I_al0) - rho_ref) - 2.0*eps * &
               I_Rho * (lambda * I_al0**2) * eps2 * &
               (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2)))
    enddo
    ! Use Bode's rule to integrate the values.
    intx_dpa(i-ioff,j-joff) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                           12.0*intz(3))
  enddo ; enddo ; endif

  if (present(inty_dpa)) then ; do J=Jsq,Jeq ; do i=is,ie
    intz(1) = dpa(i-ioff,j-joff) ; intz(5) = dpa(i-ioff,j+1-joff)
    do m=2,4
      w_left = 0.25*real(m-1) ; w_right = 1.0-w_left
      al0 = w_left*al0_2d(i,j) + w_right*al0_2d(i,j+1)
      p0 = w_left*p0_2d(i,j) + w_right*p0_2d(i,j+1)
      lambda = w_left*lambda_2d(i,j) + w_right*lambda_2d(i,j+1)

      dz = w_left*(z_t(i,j) - z_b(i,j)) + w_right*(z_t(i,j+1) - z_b(i,j+1))
      p_ave = -0.5*GxRho*(w_left*(z_t(i,j)+z_b(i,j)) + &
                          w_right*(z_t(i,j+1)+z_b(i,j+1)))

      I_al0 = 1.0 / al0
      I_Lzz = 1.0 / (p0 + (lambda * I_al0) + p_ave)
      eps = 0.5*GxRho*dz*I_Lzz ; eps2 = eps*eps

      intz(m) = G_e*dz*((p0 + p_ave)*(I_Lzz*I_al0) - rho_ref) - 2.0*eps * &
               I_Rho * (lambda * I_al0**2) * eps2 * &
               (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2)))
    enddo
    ! Use Bode's rule to integrate the values.
    inty_dpa(i-ioff,j-joff) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                           12.0*intz(3))
  enddo ; enddo ; endif
end subroutine int_density_dz_wright

!>   This subroutine calculates analytical and nearly-analytical integrals in
!! pressure across layers of geopotential anomalies, which are required for
!! calculating the finite-volume form pressure accelerations in a non-Boussinesq
!! model.  There are essentially no free assumptions, apart from the use of
!! Bode's rule to do the horizontal integrals, and from a truncation in the
!! series for log(1-eps/1+eps) that assumes that |eps| < 0.34.
subroutine int_spec_vol_dp_wright(T, S, p_t, p_b, alpha_ref, HI, dza, &
                                  intp_dza, intx_dza, inty_dza, halo_size)
  type(hor_index_type), intent(in)  :: HI
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T         !< Potential temperature relative to the surface
                                                 !! in C.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S         !< Salinity in PSU.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: p_t       !< Pressure at the top of the layer in Pa.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: p_b       !< Pressure at the top of the layer in Pa.
  real,                 intent(in)  :: alpha_ref !< A mean specific volume that is subtracted out
           !! to reduce the magnitude of each of the integrals, m3 kg-1.The calculation is
           !! mathematically identical with different values of alpha_ref, but this reduces the
           !! effects of roundoff.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(out) :: dza       !< The change in the geopotential anomaly across
                                                 !! the layer, in m2 s-2.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(out) :: intp_dza  !< The integral in pressure through the layer of
                                                 !! the geopotential anomaly relative to the anomaly
                                                 !! at the bottom of the layer, in Pa m2 s-2.
  real, dimension(HI%IsdB:HI%IedB,HI%jsd:HI%jed), &
              optional, intent(out) :: intx_dza  !< The integral in x of the difference between the
                                                 !! geopotential anomaly at the top and bottom of
                                                 !! the layer divided by the x grid spacing,
                                                 !! in m2 s-2.
  real, dimension(HI%isd:HI%ied,HI%JsdB:HI%JedB), &
              optional, intent(out) :: inty_dza  !< The integral in y of the difference between the
                                                 !! geopotential anomaly at the top and bottom of
                                                 !! the layer divided by the y grid spacing,
                                                 !! in m2 s-2.
  integer,    optional, intent(in)  :: halo_size !< The width of halo points on which to calculate
                                                 !! dza.

!   This subroutine calculates analytical and nearly-analytical integrals in
! pressure across layers of geopotential anomalies, which are required for
! calculating the finite-volume form pressure accelerations in a non-Boussinesq
! model.  There are essentially no free assumptions, apart from the use of
! Bode's rule to do the horizontal integrals, and from a truncation in the
! series for log(1-eps/1+eps) that assumes that |eps| < 0.34.
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
!  (in)      HI - The ocean's horizontal index structure.
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
!  (in,opt)  halo_size - The width of halo points on which to calculate dza.

  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed) :: al0_2d, p0_2d, lambda_2d
  real :: al0, p0, lambda
  real :: alpha_anom, dp, p_ave
  real :: rem, eps, eps2
  real :: w_left, w_right, intp(5)
  real, parameter :: C1_3 = 1.0/3.0, C1_7 = 1.0/7.0    ! Rational constants.
  real, parameter :: C1_9 = 1.0/9.0, C1_90 = 1.0/90.0  ! Rational constants.
  integer :: Isq, Ieq, Jsq, Jeq, ish, ieh, jsh, jeh, i, j, m, halo

  Isq = HI%IscB ; Ieq = HI%IecB ; Jsq = HI%JscB ; Jeq = HI%JecB
  halo = 0 ; if (present(halo_size)) halo = MAX(halo_size,0)
  ish = HI%isc-halo ; ieh = HI%iec+halo ; jsh = HI%jsc-halo ; jeh = HI%jec+halo
  if (present(intx_dza)) then ; ish = MIN(Isq,ish) ; ieh = MAX(Ieq+1,ieh); endif
  if (present(inty_dza)) then ; jsh = MIN(Jsq,jsh) ; jeh = MAX(Jeq+1,jeh); endif

  if (present(intp_dza)) then
    do j=jsh,jeh ; do i=ish,ieh
      al0_2d(i,j) = (a0 + a1*T(i,j)) + a2*S(i,j)
      p0_2d(i,j) = (b0 + b4*S(i,j)) + T(i,j) * (b1 + T(i,j)*((b2 + b3*T(i,j))) + b5*S(i,j))
      lambda_2d(i,j) = (c0 +c4*S(i,j)) + T(i,j) * (c1 + T(i,j)*((c2 + c3*T(i,j))) + c5*S(i,j))

      al0 = al0_2d(i,j) ; p0 = p0_2d(i,j) ; lambda = lambda_2d(i,j)
      dp = p_b(i,j) - p_t(i,j)
      p_ave = 0.5*(p_t(i,j)+p_b(i,j))

      eps = 0.5 * dp / (p0 + p_ave) ; eps2 = eps*eps
      alpha_anom = al0 + lambda / (p0 + p_ave) - alpha_ref
      rem = lambda * eps2 * (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2)))
      dza(i,j) = alpha_anom*dp + 2.0*eps*rem
      intp_dza(i,j) = 0.5*alpha_anom*dp**2 - dp*(1.0-eps)*rem
    enddo ; enddo
  else
    do j=jsh,jeh ; do i=ish,ieh
      al0_2d(i,j) = (a0 + a1*T(i,j)) + a2*S(i,j)
      p0_2d(i,j) = (b0 + b4*S(i,j)) + T(i,j) * (b1 + T(i,j)*((b2 + b3*T(i,j))) + b5*S(i,j))
      lambda_2d(i,j) = (c0 +c4*S(i,j)) + T(i,j) * (c1 + T(i,j)*((c2 + c3*T(i,j))) + c5*S(i,j))

      al0 = al0_2d(i,j) ; p0 = p0_2d(i,j) ; lambda = lambda_2d(i,j)
      dp = p_b(i,j) - p_t(i,j)
      p_ave = 0.5*(p_t(i,j)+p_b(i,j))

      eps = 0.5 * dp / (p0 + p_ave) ; eps2 = eps*eps
      alpha_anom = al0 + lambda / (p0 + p_ave) - alpha_ref
      rem = lambda * eps2 * (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2)))
      dza(i,j) = alpha_anom*dp + 2.0*eps*rem
    enddo ; enddo
  endif

  if (present(intx_dza)) then ; do j=HI%jsc,HI%jec ; do I=Isq,Ieq
    intp(1) = dza(i,j) ; intp(5) = dza(i+1,j)
    do m=2,4
      w_left = 0.25*real(5-m) ; w_right = 1.0-w_left
      al0 = w_left*al0_2d(i,j) + w_right*al0_2d(i+1,j)
      p0 = w_left*p0_2d(i,j) + w_right*p0_2d(i+1,j)
      lambda = w_left*lambda_2d(i,j) + w_right*lambda_2d(i+1,j)

      dp = w_left*(p_b(i,j) - p_t(i,j)) + w_right*(p_b(i+1,j) - p_t(i+1,j))
      p_ave = 0.5*(w_left*(p_t(i,j)+p_b(i,j)) + w_right*(p_t(i+1,j)+p_b(i+1,j)))

      eps = 0.5 * dp / (p0 + p_ave) ; eps2 = eps*eps
      intp(m) = (al0 + lambda / (p0 + p_ave) - alpha_ref)*dp + 2.0*eps* &
               lambda * eps2 * (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2)))
    enddo
    ! Use Bode's rule to integrate the values.
    intx_dza(i,j) = C1_90*(7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4)) + &
                           12.0*intp(3))
  enddo ; enddo ; endif

  if (present(inty_dza)) then ; do J=Jsq,Jeq ; do i=HI%isc,HI%iec
    intp(1) = dza(i,j) ; intp(5) = dza(i,j+1)
    do m=2,4
      w_left = 0.25*real(5-m) ; w_right = 1.0-w_left
      al0 = w_left*al0_2d(i,j) + w_right*al0_2d(i,j+1)
      p0 = w_left*p0_2d(i,j) + w_right*p0_2d(i,j+1)
      lambda = w_left*lambda_2d(i,j) + w_right*lambda_2d(i,j+1)

      dp = w_left*(p_b(i,j) - p_t(i,j)) + w_right*(p_b(i,j+1) - p_t(i,j+1))
      p_ave = 0.5*(w_left*(p_t(i,j)+p_b(i,j)) + w_right*(p_t(i,j+1)+p_b(i,j+1)))

      eps = 0.5 * dp / (p0 + p_ave) ; eps2 = eps*eps
      intp(m) = (al0 + lambda / (p0 + p_ave) - alpha_ref)*dp + 2.0*eps* &
               lambda * eps2 * (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2)))
    enddo
    ! Use Bode's rule to integrate the values.
    inty_dza(i,j) = C1_90*(7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4)) + &
                           12.0*intp(3))
  enddo ; enddo ; endif
end subroutine int_spec_vol_dp_wright

end module MOM_EOS_Wright
