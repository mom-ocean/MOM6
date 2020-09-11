!> The equation of state using the Wright 1997 expressions
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

public calculate_compress_wright, calculate_density_wright, calculate_spec_vol_wright
public calculate_density_derivs_wright, calculate_specvol_derivs_wright
public calculate_density_second_derivs_wright
public int_density_dz_wright, int_spec_vol_dp_wright

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.


!> Compute the in situ density of sea water (in [kg m-3]), or its anomaly with respect to
!! a reference density, from salinity (in psu), potential temperature (in deg C), and pressure [Pa],
!! using the expressions from Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.
interface calculate_density_wright
  module procedure calculate_density_scalar_wright, calculate_density_array_wright
end interface calculate_density_wright

!> Compute the in situ specific volume of sea water (in [m3 kg-1]), or an anomaly with respect
!! to a reference specific volume, from salinity (in psu), potential temperature (in deg C), and
!! pressure [Pa], using the expressions from Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.
interface calculate_spec_vol_wright
  module procedure calculate_spec_vol_scalar_wright, calculate_spec_vol_array_wright
end interface calculate_spec_vol_wright

!> For a given thermodynamic state, return the derivatives of density with temperature and salinity
interface calculate_density_derivs_wright
  module procedure calculate_density_derivs_scalar_wright, calculate_density_derivs_array_wright
end interface

!> For a given thermodynamic state, return the second derivatives of density with various combinations
!! of temperature, salinity, and pressure
interface calculate_density_second_derivs_wright
  module procedure calculate_density_second_derivs_scalar_wright, calculate_density_second_derivs_array_wright
end interface

!>@{ Parameters in the Wright equation of state
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
real, parameter :: a0 = 7.057924e-4, a1 = 3.480336e-7, a2 = -1.112733e-7 ! a0/a1 ~= 2028 ; a0/a2 ~= -6343
real, parameter :: b0 = 5.790749e8,  b1 = 3.516535e6,  b2 = -4.002714e4  ! b0/b1 ~= 165  ; b0/b4 ~= 974
real, parameter :: b3 = 2.084372e2,  b4 = 5.944068e5,  b5 = -9.643486e3
real, parameter :: c0 = 1.704853e5,  c1 = 7.904722e2,  c2 = -7.984422    ! c0/c1 ~= 216  ; c0/c4 ~= -740
real, parameter :: c3 = 5.140652e-2, c4 = -2.302158e2, c5 = -3.079464
!>@}

contains

!> This subroutine computes the in situ density of sea water (rho in
!! [kg m-3]) from salinity (S [PSU]), potential temperature
!! (T [degC]), and pressure [Pa].  It uses the expression from
!! Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.
subroutine calculate_density_scalar_wright(T, S, pressure, rho, rho_ref)
  real,           intent(in)  :: T        !< Potential temperature relative to the surface [degC].
  real,           intent(in)  :: S        !< Salinity [PSU].
  real,           intent(in)  :: pressure !< pressure [Pa].
  real,           intent(out) :: rho      !< In situ density [kg m-3].
  real, optional, intent(in)  :: rho_ref  !< A reference density [kg m-3].

! *====================================================================*
! *  This subroutine computes the in situ density of sea water (rho in *
! *  [kg m-3]) from salinity (S [PSU]), potential temperature  *
! *  (T [degC]), and pressure [Pa].  It uses the expression from    *
! *  Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.                *
! *  Coded by R. Hallberg, 7/00                                        *
! *====================================================================*

  real, dimension(1) :: T0, S0, pressure0, rho0

  T0(1) = T
  S0(1) = S
  pressure0(1) = pressure

  call calculate_density_array_wright(T0, S0, pressure0, rho0, 1, 1, rho_ref)
  rho = rho0(1)

end subroutine calculate_density_scalar_wright

!> This subroutine computes the in situ density of sea water (rho in
!! [kg m-3]) from salinity (S [PSU]), potential temperature
!! (T [degC]), and pressure [Pa].  It uses the expression from
!! Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.
subroutine calculate_density_array_wright(T, S, pressure, rho, start, npts, rho_ref)
  real, dimension(:), intent(in)    :: T        !< potential temperature relative to the surface [degC].
  real, dimension(:), intent(in)    :: S        !< salinity [PSU].
  real, dimension(:), intent(in)    :: pressure !< pressure [Pa].
  real, dimension(:), intent(inout) :: rho      !< in situ density [kg m-3].
  integer,            intent(in)    :: start    !< the starting point in the arrays.
  integer,            intent(in)    :: npts     !< the number of values to calculate.
  real,     optional, intent(in)    :: rho_ref  !< A reference density [kg m-3].

  ! Original coded by R. Hallberg, 7/00, anomaly coded in 3/18.
  ! Local variables
  real :: al0, p0, lambda
  real :: al_TS, p_TSp, lam_TS, pa_000
  integer :: j

  if (present(rho_ref)) pa_000 = (b0*(1.0 - a0*rho_ref) - rho_ref*c0)
  if (present(rho_ref)) then ; do j=start,start+npts-1
    al_TS = a1*T(j) +a2*S(j)
    al0 = a0 + al_TS
    p_TSp = pressure(j) + (b4*S(j) + T(j) * (b1 + (T(j)*(b2 + b3*T(j)) + b5*S(j))))
    lam_TS = c4*S(j) + T(j) * (c1 + (T(j)*(c2 + c3*T(j)) + c5*S(j)))

    ! The following two expressions are mathematically equivalent.
    ! rho(j) = (b0 + p0_TSp) / ((c0 + lam_TS) + al0*(b0 + p0_TSp)) - rho_ref
    rho(j) = (pa_000 + (p_TSp - rho_ref*(p_TSp*al0 + (b0*al_TS + lam_TS)))) / &
             ( (c0 + lam_TS) + al0*(b0 + p_TSp) )
  enddo ; else ; do j=start,start+npts-1
    al0 = (a0 + a1*T(j)) +a2*S(j)
    p0 = (b0 + b4*S(j)) + T(j) * (b1 + T(j)*(b2 + b3*T(j)) + b5*S(j))
    lambda = (c0 +c4*S(j)) + T(j) * (c1 + T(j)*(c2 + c3*T(j)) + c5*S(j))
    rho(j) = (pressure(j) + p0) / (lambda + al0*(pressure(j) + p0))
  enddo ; endif

end subroutine calculate_density_array_wright

!> This subroutine computes the in situ specific volume of sea water (specvol in
!! [m3 kg-1]) from salinity (S [PSU]), potential temperature (T [degC])
!! and pressure [Pa].  It uses the expression from
!! Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.
!! If spv_ref is present, specvol is an anomaly from spv_ref.
subroutine calculate_spec_vol_scalar_wright(T, S, pressure, specvol, spv_ref)
  real,           intent(in)  :: T        !< potential temperature relative to the surface [degC].
  real,           intent(in)  :: S        !< salinity [PSU].
  real,           intent(in)  :: pressure !< pressure [Pa].
  real,           intent(out) :: specvol  !< in situ specific volume [m3 kg-1].
  real, optional, intent(in)  :: spv_ref  !< A reference specific volume [m3 kg-1].

  ! Local variables
  real, dimension(1) :: T0, S0, pressure0, spv0

  T0(1) = T ; S0(1) = S ; pressure0(1) = pressure

  call calculate_spec_vol_array_wright(T0, S0, pressure0, spv0, 1, 1, spv_ref)
  specvol = spv0(1)
end subroutine calculate_spec_vol_scalar_wright

!> This subroutine computes the in situ specific volume of sea water (specvol in
!! [m3 kg-1]) from salinity (S [PSU]), potential temperature (T [degC])
!! and pressure [Pa].  It uses the expression from
!! Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.
!! If spv_ref is present, specvol is an anomaly from spv_ref.
subroutine calculate_spec_vol_array_wright(T, S, pressure, specvol, start, npts, spv_ref)
  real, dimension(:), intent(in)    :: T        !< potential temperature relative to the
                                              !! surface [degC].
  real, dimension(:), intent(in)    :: S        !< salinity [PSU].
  real, dimension(:), intent(in)    :: pressure !< pressure [Pa].
  real, dimension(:), intent(inout) :: specvol  !< in situ specific volume [m3 kg-1].
  integer,            intent(in)    :: start    !< the starting point in the arrays.
  integer,            intent(in)    :: npts     !< the number of values to calculate.
  real,     optional, intent(in)    :: spv_ref  !< A reference specific volume [m3 kg-1].

  ! Local variables
  real :: al0, p0, lambda
  integer :: j

  do j=start,start+npts-1
    al0 = (a0 + a1*T(j)) +a2*S(j)
    p0 = (b0 + b4*S(j)) + T(j) * (b1 + T(j)*((b2 + b3*T(j))) + b5*S(j))
    lambda = (c0 +c4*S(j)) + T(j) * (c1 + T(j)*((c2 + c3*T(j))) + c5*S(j))

    if (present(spv_ref)) then
      specvol(j) = (lambda + (al0 - spv_ref)*(pressure(j) + p0)) / (pressure(j) + p0)
    else
      specvol(j) = (lambda + al0*(pressure(j) + p0)) / (pressure(j) + p0)
    endif
  enddo
end subroutine calculate_spec_vol_array_wright

!> For a given thermodynamic state, return the thermal/haline expansion coefficients
subroutine calculate_density_derivs_array_wright(T, S, pressure, drho_dT, drho_dS, start, npts)
  real,    intent(in),    dimension(:) :: T        !< Potential temperature relative to the
                                                   !! surface [degC].
  real,    intent(in),    dimension(:) :: S        !< Salinity [PSU].
  real,    intent(in),    dimension(:) :: pressure !< pressure [Pa].
  real,    intent(inout), dimension(:) :: drho_dT  !< The partial derivative of density with potential
                                                   !! temperature [kg m-3 degC-1].
  real,    intent(inout), dimension(:) :: drho_dS  !< The partial derivative of density with salinity,
                                                   !! in [kg m-3 PSU-1].
  integer, intent(in)                  :: start    !< The starting point in the arrays.
  integer, intent(in)                  :: npts     !< The number of values to calculate.

  ! Local variables
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
  real,    intent(in)  :: T        !< Potential temperature relative to the surface [degC].
  real,    intent(in)  :: S        !< Salinity [PSU].
  real,    intent(in)  :: pressure !< pressure [Pa].
  real,    intent(out) :: drho_dT  !< The partial derivative of density with potential
                                   !! temperature [kg m-3 degC-1].
  real,    intent(out) :: drho_dS  !< The partial derivative of density with salinity,
                                   !! in [kg m-3 PSU-1].

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
  real, dimension(:), intent(in   ) :: T !< Potential temperature referenced to 0 dbar [degC]
  real, dimension(:), intent(in   ) :: S !< Salinity [PSU]
  real, dimension(:), intent(in   ) :: P !< Pressure [Pa]
  real, dimension(:), intent(inout) :: drho_ds_ds !< Partial derivative of beta with respect
                                                  !! to S [kg m-3 PSU-2]
  real, dimension(:), intent(inout) :: drho_ds_dt !< Partial derivative of beta with respcct
                                                  !! to T [kg m-3 PSU-1 degC-1]
  real, dimension(:), intent(inout) :: drho_dt_dt !< Partial derivative of alpha with respect
                                                  !! to T [kg m-3 degC-2]
  real, dimension(:), intent(inout) :: drho_ds_dp !< Partial derivative of beta with respect
                                                  !! to pressure [kg m-3 PSU-1 Pa-1]
  real, dimension(:), intent(inout) :: drho_dt_dp !< Partial derivative of alpha with respect
                                                  !! to pressure [kg m-3 degC-1 Pa-1]
  integer,            intent(in   ) :: start !< Starting index in T,S,P
  integer,            intent(in   ) :: npts  !< Number of points to loop over

  ! Local variables
  real :: z0, z1, z2, z3, z4, z5, z6 ,z7, z8, z9, z10, z11, z2_2, z2_3
  integer :: j
  ! Based on the above expression with common terms factored, there probably exists a more numerically stable
  ! and/or efficient expression

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
  real, intent(in   ) :: S          !< Salinity [PSU]
  real, intent(in   ) :: P          !< pressure [Pa]
  real, intent(  out) :: drho_ds_ds !< Partial derivative of beta with respect
                                    !! to S [kg m-3 PSU-2]
  real, intent(  out) :: drho_ds_dt !< Partial derivative of beta with respcct
                                    !! to T [kg m-3 PSU-1 degC-1]
  real, intent(  out) :: drho_dt_dt !< Partial derivative of alpha with respect
                                    !! to T [kg m-3 degC-2]
  real, intent(  out) :: drho_ds_dp !< Partial derivative of beta with respect
                                    !! to pressure [kg m-3 PSU-1 Pa-1]
  real, intent(  out) :: drho_dt_dp !< Partial derivative of alpha with respect
                                    !! to pressure [kg m-3 degC-1 Pa-1]
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

!> For a given thermodynamic state, return the partial derivatives of specific volume
!! with temperature and salinity
subroutine calculate_specvol_derivs_wright(T, S, pressure, dSV_dT, dSV_dS, start, npts)
  real,    intent(in),    dimension(:) :: T        !< Potential temperature relative to the surface [degC].
  real,    intent(in),    dimension(:) :: S        !< Salinity [PSU].
  real,    intent(in),    dimension(:) :: pressure !< pressure [Pa].
  real,    intent(inout), dimension(:) :: dSV_dT   !< The partial derivative of specific volume with
                                                   !! potential temperature [m3 kg-1 degC-1].
  real,    intent(inout), dimension(:) :: dSV_dS   !< The partial derivative of specific volume with
                                                   !! salinity [m3 kg-1 / Pa].
  integer, intent(in)                  :: start    !< The starting point in the arrays.
  integer, intent(in)                  :: npts     !< The number of values to calculate.

  ! Local variables
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
!! [kg m-3]) and the compressibility (drho/dp = C_sound^-2)
!! (drho_dp [s2 m-2]) from salinity (sal in psu), potential
!! temperature (T [degC]), and pressure [Pa].  It uses the expressions
!! from Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.
!! Coded by R. Hallberg, 1/01
subroutine calculate_compress_wright(T, S, pressure, rho, drho_dp, start, npts)
  real,    intent(in),    dimension(:) :: T        !< Potential temperature relative to the surface [degC].
  real,    intent(in),    dimension(:) :: S        !< Salinity [PSU].
  real,    intent(in),    dimension(:) :: pressure !< pressure [Pa].
  real,    intent(inout), dimension(:) :: rho      !< In situ density [kg m-3].
  real,    intent(inout), dimension(:) :: drho_dp  !< The partial derivative of density with pressure
                                                   !! (also the inverse of the square of sound speed)
                                                   !! [s2 m-2].
  integer, intent(in)                  :: start    !< The starting point in the arrays.
  integer, intent(in)                  :: npts     !< The number of values to calculate.

  ! Coded by R. Hallberg, 1/01
  ! Local variables
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
subroutine int_density_dz_wright(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                                 dpa, intz_dpa, intx_dpa, inty_dpa, &
                                 bathyT, dz_neglect, useMassWghtInterp, rho_scale, pres_scale)
  type(hor_index_type), intent(in)  :: HI       !< The horizontal index type for the arrays.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T        !< Potential temperature relative to the surface
                                                !! [degC].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S        !< Salinity [PSU].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_t      !< Height at the top of the layer in depth units [Z ~> m].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_b      !< Height at the top of the layer [Z ~> m].
  real,                 intent(in)  :: rho_ref  !< A mean density [R ~> kg m-3] or [kg m-3], that is subtracted
                                                !! out to reduce the magnitude of each of the integrals.
                                                !! (The pressure is calucated as p~=-z*rho_0*G_e.)
  real,                 intent(in)  :: rho_0    !< Density [R ~> kg m-3] or [kg m-3], that is used
                                                !! to calculate the pressure (as p~=-z*rho_0*G_e)
                                                !! used in the equation of state.
  real,                 intent(in)  :: G_e      !< The Earth's gravitational acceleration
                                                !! [L2 Z-1 T-2 ~> m s-2] or [m2 Z-1 s-2 ~> m s-2].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(inout) :: dpa    !< The change in the pressure anomaly across the
                                                !! layer [R L2 T-2 ~> Pa] or [Pa].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(inout) :: intz_dpa !< The integral through the thickness of the layer
                                                !! of the pressure anomaly relative to the anomaly
                                                !! at the top of the layer [R Z L2 T-2 ~> Pa m].
  real, dimension(HI%IsdB:HI%IedB,HI%jsd:HI%jed), &
              optional, intent(inout) :: intx_dpa !< The integral in x of the difference between the
                                                !! pressure anomaly at the top and bottom of the
                                                !! layer divided by the x grid spacing [R L2 T-2 ~> Pa].
  real, dimension(HI%isd:HI%ied,HI%JsdB:HI%JedB), &
              optional, intent(inout) :: inty_dpa !< The integral in y of the difference between the
                                                !! pressure anomaly at the top and bottom of the
                                                !! layer divided by the y grid spacing [R L2 T-2 ~> Pa].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: bathyT   !< The depth of the bathymetry [Z ~> m].
  real,       optional, intent(in)  :: dz_neglect !< A miniscule thickness change [Z ~> m].
  logical,    optional, intent(in)  :: useMassWghtInterp !< If true, uses mass weighting to
                                                !! interpolate T/S for top and bottom integrals.
  real,       optional, intent(in)  :: rho_scale !< A multiplicative factor by which to scale density
                                                 !! from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real,       optional, intent(in)  :: pres_scale !< A multiplicative factor to convert pressure
                                                 !! into Pa [Pa T2 R-1 L-2 ~> 1].

  ! Local variables
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed) :: al0_2d, p0_2d, lambda_2d
  real :: al0, p0, lambda
  real :: rho_anom   ! The density anomaly from rho_ref [kg m-3].
  real :: eps, eps2, rem
  real :: GxRho      ! The gravitational acceleration times density and unit conversion factors [Pa Z-1 ~> kg m-2 s-2]
  real :: g_Earth    ! The gravitational acceleration [m2 Z-1 s-2 ~> m s-2]
  real :: I_Rho      ! The inverse of the Boussinesq density [m3 kg-1]
  real :: rho_ref_mks ! The reference density in MKS units, never rescaled from kg m-3 [kg m-3]
  real :: p_ave, I_al0, I_Lzz
  real :: dz         ! The layer thickness [Z ~> m].
  real :: hWght      ! A pressure-thickness below topography [Z ~> m].
  real :: hL, hR     ! Pressure-thicknesses of the columns to the left and right [Z ~> m].
  real :: iDenom     ! The inverse of the denominator in the weights [Z-2 ~> m-2].
  real :: hWt_LL, hWt_LR ! hWt_LA is the weighted influence of A on the left column [nondim].
  real :: hWt_RL, hWt_RR ! hWt_RA is the weighted influence of A on the right column [nondim].
  real :: wt_L, wt_R ! The linear weights of the left and right columns [nondim].
  real :: wtT_L, wtT_R ! The weights for tracers from the left and right columns [nondim].
  real :: intz(5)    ! The integrals of density with height at the
                     ! 5 sub-column locations [R L2 T-2 ~> Pa].
  real :: Pa_to_RL2_T2 ! A conversion factor of pressures from Pa to the output units indicated by
                       ! pres_scale [R L2 T-2 Pa-1 ~> 1] or [1].
  logical :: do_massWeight ! Indicates whether to do mass weighting.
  real, parameter :: C1_3 = 1.0/3.0, C1_7 = 1.0/7.0    ! Rational constants.
  real, parameter :: C1_9 = 1.0/9.0, C1_90 = 1.0/90.0  ! Rational constants.
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, i, j, m

  ! These array bounds work for the indexing convention of the input arrays, but
  ! on the computational domain defined for the output arrays.
  Isq = HI%IscB ; Ieq = HI%IecB
  Jsq = HI%JscB ; Jeq = HI%JecB
  is = HI%isc ; ie = HI%iec
  js = HI%jsc ; je = HI%jec

  if (present(pres_scale)) then
    GxRho = pres_scale * G_e * rho_0 ; g_Earth = pres_scale * G_e
    Pa_to_RL2_T2 = 1.0 / pres_scale
  else
    GxRho = G_e * rho_0 ; g_Earth = G_e
    Pa_to_RL2_T2 = 1.0
  endif
  if (present(rho_scale)) then
    g_Earth = g_Earth * rho_scale
    rho_ref_mks = rho_ref / rho_scale ; I_Rho = rho_scale / rho_0
  else
    rho_ref_mks = rho_ref ; I_Rho = 1.0 / rho_0
  endif

  do_massWeight = .false.
  if (present(useMassWghtInterp)) then ; if (useMassWghtInterp) then
    do_massWeight = .true.
  ! if (.not.present(bathyT)) call MOM_error(FATAL, "int_density_dz_generic: "//&
  !     "bathyT must be present if useMassWghtInterp is present and true.")
  ! if (.not.present(dz_neglect)) call MOM_error(FATAL, "int_density_dz_generic: "//&
  !     "dz_neglect must be present if useMassWghtInterp is present and true.")
  endif ; endif

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

    rho_anom = (p0 + p_ave)*(I_Lzz*I_al0) - rho_ref_mks
    rem = I_Rho * (lambda * I_al0**2) * eps2 * &
          (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2)))
    dpa(i,j) = Pa_to_RL2_T2 * (g_Earth*rho_anom*dz - 2.0*eps*rem)
    if (present(intz_dpa)) &
      intz_dpa(i,j) = Pa_to_RL2_T2 * (0.5*g_Earth*rho_anom*dz**2 - dz*(1.0+eps)*rem)
  enddo ; enddo

  if (present(intx_dpa)) then ; do j=js,je ; do I=Isq,Ieq
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., -bathyT(i,j)-z_t(i+1,j), -bathyT(i+1,j)-z_t(i,j))
    if (hWght > 0.) then
      hL = (z_t(i,j) - z_b(i,j)) + dz_neglect
      hR = (z_t(i+1,j) - z_b(i+1,j)) + dz_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    intz(1) = dpa(i,j) ; intz(5) = dpa(i+1,j)
    do m=2,4
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR

      al0 = wtT_L*al0_2d(i,j) + wtT_R*al0_2d(i+1,j)
      p0 = wtT_L*p0_2d(i,j) + wtT_R*p0_2d(i+1,j)
      lambda = wtT_L*lambda_2d(i,j) + wtT_R*lambda_2d(i+1,j)

      dz = wt_L*(z_t(i,j) - z_b(i,j)) + wt_R*(z_t(i+1,j) - z_b(i+1,j))
      p_ave = -0.5*GxRho*(wt_L*(z_t(i,j)+z_b(i,j)) + &
                          wt_R*(z_t(i+1,j)+z_b(i+1,j)))

      I_al0 = 1.0 / al0
      I_Lzz = 1.0 / (p0 + (lambda * I_al0) + p_ave)
      eps = 0.5*GxRho*dz*I_Lzz ; eps2 = eps*eps

      intz(m) = Pa_to_RL2_T2 * ( g_Earth*dz*((p0 + p_ave)*(I_Lzz*I_al0) - rho_ref_mks) - 2.0*eps * &
                I_Rho * (lambda * I_al0**2) * eps2 * (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2))) )
    enddo
    ! Use Bode's rule to integrate the values.
    intx_dpa(i,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + 12.0*intz(3))
  enddo ; enddo ; endif

  if (present(inty_dpa)) then ; do J=Jsq,Jeq ; do i=is,ie
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., -bathyT(i,j)-z_t(i,j+1), -bathyT(i,j+1)-z_t(i,j))
    if (hWght > 0.) then
      hL = (z_t(i,j) - z_b(i,j)) + dz_neglect
      hR = (z_t(i,j+1) - z_b(i,j+1)) + dz_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    intz(1) = dpa(i,j) ; intz(5) = dpa(i,j+1)
    do m=2,4
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR

      al0 = wtT_L*al0_2d(i,j) + wtT_R*al0_2d(i,j+1)
      p0 = wtT_L*p0_2d(i,j) + wtT_R*p0_2d(i,j+1)
      lambda = wtT_L*lambda_2d(i,j) + wtT_R*lambda_2d(i,j+1)

      dz = wt_L*(z_t(i,j) - z_b(i,j)) + wt_R*(z_t(i,j+1) - z_b(i,j+1))
      p_ave = -0.5*GxRho*(wt_L*(z_t(i,j)+z_b(i,j)) + &
                          wt_R*(z_t(i,j+1)+z_b(i,j+1)))

      I_al0 = 1.0 / al0
      I_Lzz = 1.0 / (p0 + (lambda * I_al0) + p_ave)
      eps = 0.5*GxRho*dz*I_Lzz ; eps2 = eps*eps

      intz(m) = Pa_to_RL2_T2 * ( g_Earth*dz*((p0 + p_ave)*(I_Lzz*I_al0) - rho_ref_mks) - 2.0*eps * &
                I_Rho * (lambda * I_al0**2) * eps2 * (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2))) )
    enddo
    ! Use Bode's rule to integrate the values.
    inty_dpa(i,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + 12.0*intz(3))
  enddo ; enddo ; endif

end subroutine int_density_dz_wright

!>   This subroutine calculates analytical and nearly-analytical integrals in
!! pressure across layers of geopotential anomalies, which are required for
!! calculating the finite-volume form pressure accelerations in a non-Boussinesq
!! model.  There are essentially no free assumptions, apart from the use of
!! Bode's rule to do the horizontal integrals, and from a truncation in the
!! series for log(1-eps/1+eps) that assumes that |eps| < 0.34.
subroutine int_spec_vol_dp_wright(T, S, p_t, p_b, spv_ref, HI, dza, &
                                  intp_dza, intx_dza, inty_dza, halo_size, &
                                  bathyP, dP_neglect, useMassWghtInterp, SV_scale, pres_scale)
  type(hor_index_type), intent(in)  :: HI        !< The ocean's horizontal index type.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T         !< Potential temperature relative to the surface
                                                 !! [degC].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S         !< Salinity [PSU].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: p_t       !< Pressure at the top of the layer [R L2 T-2 ~> Pa] or [Pa].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: p_b       !< Pressure at the top of the layer [R L2 T-2 ~> Pa] or [Pa].
  real,                 intent(in)  :: spv_ref   !< A mean specific volume that is subtracted out
                            !! to reduce the magnitude of each of the integrals [R-1 ~> m3 kg-1].
                            !! The calculation is mathematically identical with different values of
                            !! spv_ref, but this reduces the effects of roundoff.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(inout) :: dza     !< The change in the geopotential anomaly across
                                                 !! the layer [T-2 ~> m2 s-2] or [m2 s-2].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(inout) :: intp_dza !< The integral in pressure through the layer of
                                                 !! the geopotential anomaly relative to the anomaly
                                                 !! at the bottom of the layer [R L4 T-4 ~> Pa m2 s-2]
                                                 !! or [Pa m2 s-2].
  real, dimension(HI%IsdB:HI%IedB,HI%jsd:HI%jed), &
              optional, intent(inout) :: intx_dza !< The integral in x of the difference between the
                                                 !! geopotential anomaly at the top and bottom of
                                                 !! the layer divided by the x grid spacing
                                                 !! [L2 T-2 ~> m2 s-2] or [m2 s-2].
  real, dimension(HI%isd:HI%ied,HI%JsdB:HI%JedB), &
              optional, intent(inout) :: inty_dza !< The integral in y of the difference between the
                                                 !! geopotential anomaly at the top and bottom of
                                                 !! the layer divided by the y grid spacing
                                                 !! [L2 T-2 ~> m2 s-2] or [m2 s-2].
  integer,    optional, intent(in)  :: halo_size !< The width of halo points on which to calculate
                                                 !! dza.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: bathyP    !< The pressure at the bathymetry [R L2 T-2 ~> Pa] or [Pa]
  real,       optional, intent(in)  :: dP_neglect !< A miniscule pressure change with
                                                 !! the same units as p_t [R L2 T-2 ~> Pa] or [Pa]
  logical,    optional, intent(in)  :: useMassWghtInterp !< If true, uses mass weighting
                            !! to interpolate T/S for top and bottom integrals.
  real,       optional, intent(in)  :: SV_scale  !< A multiplicative factor by which to scale specific
                            !! volume from m3 kg-1 to the desired units [kg m-3 R-1 ~> 1]
  real,       optional, intent(in)  :: pres_scale !< A multiplicative factor to convert pressure
                            !! into Pa [Pa T2 R-1 L-2 ~> 1].

  ! Local variables
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed) :: al0_2d, p0_2d, lambda_2d
  real :: al0        ! A term in the Wright EOS [R-1 ~> m3 kg-1]
  real :: p0         ! A term in the Wright EOS [R L2 T-2 ~> Pa]
  real :: lambda     ! A term in the Wright EOS [L2 T-2 ~> m2 s-2]
  real :: al0_scale  ! Scaling factor to convert al0 from MKS units [R-1 kg m-3 ~> 1]
  real :: p0_scale   ! Scaling factor to convert p0 from MKS units [R L2 T-2 Pa-1 ~> 1]
  real :: lam_scale  ! Scaling factor to convert lambda from MKS units [L2 s2 T-2 m-2 ~> 1]
  real :: p_ave      ! The layer average pressure [R L2 T-2 ~> Pa]
  real :: rem        ! [L2 T-2 ~> m2 s-2]
  real :: eps, eps2  ! A nondimensional ratio and its square [nondim]
  real :: alpha_anom ! The depth averaged specific volume anomaly [R-1 ~> m3 kg-1].
  real :: dp         ! The pressure change through a layer [R L2 T-2 ~> Pa].
  real :: hWght      ! A pressure-thickness below topography [R L2 T-2 ~> Pa].
  real :: hL, hR     ! Pressure-thicknesses of the columns to the left and right [R L2 T-2 ~> Pa].
  real :: iDenom     ! The inverse of the denominator in the weights [T4 R-2 L-4 ~> Pa-2].
  real :: hWt_LL, hWt_LR ! hWt_LA is the weighted influence of A on the left column [nondim].
  real :: hWt_RL, hWt_RR ! hWt_RA is the weighted influence of A on the right column [nondim].
  real :: wt_L, wt_R ! The linear weights of the left and right columns [nondim].
  real :: wtT_L, wtT_R ! The weights for tracers from the left and right columns [nondim].
  real :: intp(5)    ! The integrals of specific volume with pressure at the
                     ! 5 sub-column locations [L2 T-2 ~> m2 s-2].
  logical :: do_massWeight ! Indicates whether to do mass weighting.
  real, parameter :: C1_3 = 1.0/3.0, C1_7 = 1.0/7.0    ! Rational constants.
  real, parameter :: C1_9 = 1.0/9.0, C1_90 = 1.0/90.0  ! Rational constants.
  integer :: Isq, Ieq, Jsq, Jeq, ish, ieh, jsh, jeh, i, j, m, halo

  Isq = HI%IscB ; Ieq = HI%IecB ; Jsq = HI%JscB ; Jeq = HI%JecB
  halo = 0 ; if (present(halo_size)) halo = MAX(halo_size,0)
  ish = HI%isc-halo ; ieh = HI%iec+halo ; jsh = HI%jsc-halo ; jeh = HI%jec+halo
  if (present(intx_dza)) then ; ish = MIN(Isq,ish) ; ieh = MAX(Ieq+1,ieh); endif
  if (present(inty_dza)) then ; jsh = MIN(Jsq,jsh) ; jeh = MAX(Jeq+1,jeh); endif


  al0_scale = 1.0 ; if (present(SV_scale)) al0_scale = SV_scale
  p0_scale = 1.0
  if (present(pres_scale)) then ; if (pres_scale /= 1.0) then
    p0_scale = 1.0 / pres_scale
  endif ; endif
  lam_scale = al0_scale * p0_scale

  do_massWeight = .false.
  if (present(useMassWghtInterp)) then ; if (useMassWghtInterp) then
    do_massWeight = .true.
!    if (.not.present(bathyP)) call MOM_error(FATAL, "int_spec_vol_dp_generic: "//&
!        "bathyP must be present if useMassWghtInterp is present and true.")
!    if (.not.present(dP_neglect)) call MOM_error(FATAL, "int_spec_vol_dp_generic: "//&
!        "dP_neglect must be present if useMassWghtInterp is present and true.")
  endif ; endif

  !  alpha(j) = (lambda + al0*(pressure(j) + p0)) / (pressure(j) + p0)
  do j=jsh,jeh ; do i=ish,ieh
    al0_2d(i,j) = al0_scale * ( (a0 + a1*T(i,j)) + a2*S(i,j) )
    p0_2d(i,j) = p0_scale * ( (b0 + b4*S(i,j)) + T(i,j) * (b1 + T(i,j)*((b2 + b3*T(i,j))) + b5*S(i,j)) )
    lambda_2d(i,j) = lam_scale * ( (c0 + c4*S(i,j)) + T(i,j) * (c1 + T(i,j)*((c2 + c3*T(i,j))) + c5*S(i,j)) )

    al0 = al0_2d(i,j) ; p0 = p0_2d(i,j) ; lambda = lambda_2d(i,j)
    dp = p_b(i,j) - p_t(i,j)
    p_ave = 0.5*(p_t(i,j)+p_b(i,j))

    eps = 0.5 * dp / (p0 + p_ave) ; eps2 = eps*eps
    alpha_anom = al0 + lambda / (p0 + p_ave) - spv_ref
    rem = lambda * eps2 * (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2)))
    dza(i,j) = alpha_anom*dp + 2.0*eps*rem
    if (present(intp_dza)) &
      intp_dza(i,j) = 0.5*alpha_anom*dp**2 - dp*(1.0-eps)*rem
  enddo ; enddo

  if (present(intx_dza)) then ; do j=HI%jsc,HI%jec ; do I=Isq,Ieq
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., bathyP(i,j)-p_t(i+1,j), bathyP(i+1,j)-p_t(i,j))
    if (hWght > 0.) then
      hL = (p_b(i,j) - p_t(i,j)) + dP_neglect
      hR = (p_b(i+1,j) - p_t(i+1,j)) + dP_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    intp(1) = dza(i,j) ; intp(5) = dza(i+1,j)
    do m=2,4
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR

      ! T, S, and p are interpolated in the horizontal.  The p interpolation
      ! is linear, but for T and S it may be thickness wekghted.
      al0 = wtT_L*al0_2d(i,j) + wtT_R*al0_2d(i+1,j)
      p0 = wtT_L*p0_2d(i,j) + wtT_R*p0_2d(i+1,j)
      lambda = wtT_L*lambda_2d(i,j) + wtT_R*lambda_2d(i+1,j)

      dp = wt_L*(p_b(i,j) - p_t(i,j)) + wt_R*(p_b(i+1,j) - p_t(i+1,j))
      p_ave = 0.5*(wt_L*(p_t(i,j)+p_b(i,j)) + wt_R*(p_t(i+1,j)+p_b(i+1,j)))

      eps = 0.5 * dp / (p0 + p_ave) ; eps2 = eps*eps
      intp(m) = (al0 + lambda / (p0 + p_ave) - spv_ref)*dp + 2.0*eps* &
               lambda * eps2 * (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2)))
    enddo
    ! Use Bode's rule to integrate the values.
    intx_dza(i,j) = C1_90*(7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4)) + &
                           12.0*intp(3))
  enddo ; enddo ; endif

  if (present(inty_dza)) then ; do J=Jsq,Jeq ; do i=HI%isc,HI%iec
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., bathyP(i,j)-p_t(i,j+1), bathyP(i,j+1)-p_t(i,j))
    if (hWght > 0.) then
      hL = (p_b(i,j) - p_t(i,j)) + dP_neglect
      hR = (p_b(i,j+1) - p_t(i,j+1)) + dP_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    intp(1) = dza(i,j) ; intp(5) = dza(i,j+1)
    do m=2,4
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR

      ! T, S, and p are interpolated in the horizontal.  The p interpolation
      ! is linear, but for T and S it may be thickness wekghted.
      al0 = wt_L*al0_2d(i,j) + wt_R*al0_2d(i,j+1)
      p0 = wt_L*p0_2d(i,j) + wt_R*p0_2d(i,j+1)
      lambda = wt_L*lambda_2d(i,j) + wt_R*lambda_2d(i,j+1)

      dp = wt_L*(p_b(i,j) - p_t(i,j)) + wt_R*(p_b(i,j+1) - p_t(i,j+1))
      p_ave = 0.5*(wt_L*(p_t(i,j)+p_b(i,j)) + wt_R*(p_t(i,j+1)+p_b(i,j+1)))

      eps = 0.5 * dp / (p0 + p_ave) ; eps2 = eps*eps
      intp(m) = (al0 + lambda / (p0 + p_ave) - spv_ref)*dp + 2.0*eps* &
               lambda * eps2 * (C1_3 + eps2*(0.2 + eps2*(C1_7 + C1_9*eps2)))
    enddo
    ! Use Bode's rule to integrate the values.
    inty_dza(i,j) = C1_90*(7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4)) + &
                           12.0*intp(3))
  enddo ; enddo ; endif
end subroutine int_spec_vol_dp_wright

end module MOM_EOS_Wright
