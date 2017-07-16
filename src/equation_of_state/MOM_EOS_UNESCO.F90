module MOM_EOS_UNESCO
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
!*  sea water using the fit to the UNESCO equation of state given by   *
!*  the expressions from Jackett and McDougall, 1995, J. Atmos.        *
!*  Ocean. Tech., 12, 381-389.  Coded by J. Stephens, 9/99.            *
!***********************************************************************

implicit none ; private

public calculate_compress_UNESCO, calculate_density_UNESCO
public calculate_density_derivs_UNESCO
public calculate_density_scalar_UNESCO, calculate_density_array_UNESCO

interface calculate_density_UNESCO
  module procedure calculate_density_scalar_UNESCO, calculate_density_array_UNESCO
end interface calculate_density_UNESCO

! The following constants are used to calculate rho0.  The notation
! is Rab for the contribution to rho0 from T^aS^b.
real, parameter ::  R00 = 999.842594, R10 = 6.793952e-2, R20 = -9.095290e-3, &
  R30 = 1.001685e-4, R40 = -1.120083e-6, R50 = 6.536332e-9, R01 = 0.824493, &
  R11 = -4.0899e-3, R21 = 7.6438e-5, R31 = -8.2467e-7, R41 = 5.3875e-9, &
  R032 = -5.72466e-3, R132 = 1.0227e-4, R232 = -1.6546e-6, R02 = 4.8314e-4;

! The following constants are used to calculate the secant bulk mod-
! ulus. The notation here is Sab for terms proportional to T^a*S^b,
! Spab for terms proportional to p*T^a*S^b, and SPab for terms
! proportional to p^2*T^a*S^b.
real, parameter ::  S00 = 1.965933e4, S10 = 1.444304e2, S20 = -1.706103, &
  S30 = 9.648704e-3, S40 = -4.190253e-5, S01 = 52.84855, S11 = -3.101089e-1, &
  S21 = 6.283263e-3, S31 = -5.084188e-5, S032 = 3.886640e-1, S132 = 9.085835e-3, &
  S232 = -4.619924e-4, Sp00 = 3.186519, Sp10 = 2.212276e-2, Sp20 = -2.984642e-4, &
  Sp30 = 1.956415e-6, Sp01 = 6.704388e-3, Sp11 = -1.847318e-4, Sp21 = 2.059331e-7, &
  Sp032 = 1.480266e-4, SP000 = 2.102898e-4, SP010 = -1.202016e-5, SP020 = 1.394680e-7, &
  SP001 = -2.040237e-6, SP011 = 6.128773e-8, SP021 = 6.207323e-10;


contains

subroutine calculate_density_scalar_UNESCO(T, S, pressure, rho)
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

  real, dimension(1) :: T0, S0, pressure0
  real, dimension(1) :: rho0

  T0(1) = T
  S0(1) = S
  pressure0(1) = pressure

  call calculate_density_array_UNESCO(T0, S0, pressure0, rho0, 1, 1)
  rho = rho0(1)

end subroutine calculate_density_scalar_UNESCO

subroutine calculate_density_array_UNESCO(T, S, pressure, rho, start, npts)
  real,    intent(in),  dimension(:) :: T, S, pressure
  real,    intent(out), dimension(:) :: rho
  integer, intent(in)                :: start, npts

! *  This subroutine computes the in situ density of sea water (rho in *
! *  units of kg/m^3) from salinity (S in psu), potential temperature  *
! *  (T in deg C), and pressure in Pa.                                 *

! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     rho - in situ density in kg m-3.                        *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *

  real :: t_local, t2, t3, t4, t5; ! Temperature to the 1st - 5th power.
  real :: s_local, s32, s2;        ! Salinity to the 1st, 3/2, & 2nd power.
  real :: p1, p2;      ! Pressure (in bars) to the 1st and 2nd power.
  real :: rho0;        ! Density at 1 bar pressure, in kg m-3.
  real :: ks;          ! The secant bulk modulus in bar.
  integer :: j

  do j=start,start+npts-1
    p1 = pressure(j)*1.0e-5; p2 = p1*p1;
    t_local = T(j); t2 = t_local*t_local; t3 = t_local*t2; t4 = t2*t2; t5 = t3*t2;
    s_local = S(j); s2 = s_local*s_local; s32 = s_local*sqrt(s_local);

!  Compute rho(s,theta,p=0) - (same as rho(s,t_insitu,p=0) ).

    rho0 = R00 + R10*t_local + R20*t2 + R30*t3 + R40*t4 + R50*t5 + &
           s_local*(R01 + R11*t_local + R21*t2 + R31*t3 + R41*t4) + &
           s32*(R032 + R132*t_local + R232*t2) + R02*s2;

!  Compute rho(s,theta,p), first calculating the secant bulk modulus.

    ks = S00 + S10*t_local + S20*t2 + S30*t3 + S40*t4 + s_local*(S01 + S11*t_local + S21*t2 + S31*t3) + &
         s32*(S032 + S132*t_local + S232*t2) + &
         p1*(Sp00 + Sp10*t_local + Sp20*t2 + Sp30*t3 + &
             s_local*(Sp01 + Sp11*t_local + Sp21*t2) + Sp032*s32) + &
         p2*(SP000 + SP010*t_local + SP020*t2 + s_local*(SP001 + SP011*t_local + SP021*t2));

    rho(j) = rho0*ks / (ks - p1);
  enddo
end subroutine calculate_density_array_UNESCO

subroutine calculate_density_derivs_UNESCO(T, S, pressure, drho_dT, drho_dS, start, npts)
  real,    intent(in),  dimension(:) ::  T, S, pressure
  real,    intent(out), dimension(:) :: drho_dT, drho_dS
  integer, intent(in)                :: start, npts
! *   This subroutine calculates the partial derivatives of density    *
! * with potential temperature and salinity.                           *
! *                                                                    *
! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     drho_dT - the partial derivative of density with        *
! *                      potential temperature, in kg m-3 K-1.         *
! *  (out)     drho_dS - the partial derivative of density with        *
! *                      salinity, in kg m-3 psu-1.                    *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
  real :: t_local, t2, t3, t4, t5; ! Temperature to the 1st - 5th power.
  real :: s12, s_local, s32, s2;   ! Salinity to the 1/2 - 2nd powers.
  real :: p1, p2;         ! Pressure (in bars) to the 1st & 2nd power.
  real :: rho0;           ! Density at 1 bar pressure, in kg m-3.
  real :: ks;             ! The secant bulk modulus, in bar.
  real :: drho0_dT;       ! Derivative of rho0 with T, in kg m-3 K-1.
  real :: drho0_dS;       ! Derivative of rho0 with S, kg m-3 psu-1.
  real :: dks_dT;         ! Derivative of ks with T, in bar K-1.
  real :: dks_dS;         ! Derivative of ks with S, in bar psu-1.
  real :: denom;          ! 1.0 / (ks - p1) in bar-1.
  integer :: j

  do j=start,start+npts-1
    p1 = pressure(j)*1.0e-5; p2 = p1*p1;
    t_local = T(j); t2 = t_local*t_local; t3 = t_local*t2; t4 = t2*t2; t5 = t3*t2;
    s_local = S(j); s2 = s_local*s_local; s12 = sqrt(s_local); s32 = s_local*s12;

!       compute rho(s,theta,p=0) - (same as rho(s,t_insitu,p=0) )

    rho0 = R00 + R10*t_local + R20*t2 + R30*t3 + R40*t4 + R50*t5 + &
           s_local*(R01 + R11*t_local + R21*t2 + R31*t3 + R41*t4) + &
           s32*(R032 + R132*t_local + R232*t2) + R02*s2;
    drho0_dT = R10 + 2.0*R20*t_local + 3.0*R30*t2 + 4.0*R40*t3 + 5.0*R50*t4 + &
               s_local*(R11 + 2.0*R21*t_local + 3.0*R31*t2 + 4.0*R41*t3) + &
               s32*(R132 + 2.0*R232*t_local);
    drho0_dS = (R01 + R11*t_local + R21*t2 + R31*t3 + R41*t4) + &
               1.5*s12*(R032 + R132*t_local + R232*t2) + 2.0*R02*s_local;

!       compute rho(s,theta,p)

    ks = S00 + S10*t_local + S20*t2 + S30*t3 + S40*t4 + s_local*(S01 + S11*t_local + S21*t2 + S31*t3) + &
         s32*(S032 + S132*t_local + S232*t2) + &
         p1*(Sp00 + Sp10*t_local + Sp20*t2 + Sp30*t3 + &
             s_local*(Sp01 + Sp11*t_local + Sp21*t2) + Sp032*s32) + &
         p2*(SP000 + SP010*t_local + SP020*t2 + s_local*(SP001 + SP011*t_local + SP021*t2));
    dks_dT = S10 + 2.0*S20*t_local + 3.0*S30*t2 + 4.0*S40*t3 + &
             s_local*(S11 + 2.0*S21*t_local + 3.0*S31*t2) + s32*(S132 + 2.0*S232*t_local) + &
             p1*(Sp10 + 2.0*Sp20*t_local + 3.0*Sp30*t2 + s_local*(Sp11 + 2.0*Sp21*t_local)) + &
             p2*(SP010 + 2.0*SP020*t_local + s_local*(SP011 + 2.0*SP021*t_local));
    dks_dS = (S01 + S11*t_local + S21*t2 + S31*t3) + 1.5*s12*(S032 + S132*t_local + S232*t2) + &
             p1*(Sp01 + Sp11*t_local + Sp21*t2 + 1.5*Sp032*s12) + &
             p2*(SP001 + SP011*t_local + SP021*t2);

    denom = 1.0 / (ks - p1);
    drho_dT(j) = denom*(ks*drho0_dT - rho0*p1*denom*dks_dT);
    drho_dS(j) = denom*(ks*drho0_dS - rho0*p1*denom*dks_dS);
  enddo

end subroutine calculate_density_derivs_UNESCO

subroutine calculate_compress_UNESCO(T, S, pressure, rho, drho_dp, start, npts)
  real,    intent(in),  dimension(:) :: T, S, pressure
  real,    intent(out), dimension(:) :: rho, drho_dp
  integer, intent(in)                :: start, npts
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

  real :: t_local, t2, t3, t4, t5; ! Temperature to the 1st - 5th power.
  real :: s_local, s32, s2;        ! Salinity to the 1st, 3/2, & 2nd power.
  real :: p1, p2;      ! Pressure (in bars) to the 1st and 2nd power.
  real :: rho0;        ! Density at 1 bar pressure, in kg m-3.
  real :: ks;          ! The secant bulk modulus in bar.
  real :: ks_0, ks_1, ks_2;
  real :: dks_dp;      ! The derivative of the secant bulk modulus
                       ! with pressure, nondimensional.
  integer :: j

  do j=start,start+npts-1
    p1 = pressure(j)*1.0e-5; p2 = p1*p1;
    t_local = T(j); t2 = t_local*t_local; t3 = t_local*t2; t4 = t2*t2; t5 = t3*t2;
    s_local = S(j); s2 = s_local*s_local; s32 = s_local*sqrt(s_local);

!  Compute rho(s,theta,p=0) - (same as rho(s,t_insitu,p=0) ).

    rho0 = R00 + R10*t_local + R20*t2 + R30*t3 + R40*t4 + R50*t5 + &
           s_local*(R01 + R11*t_local + R21*t2 + R31*t3 + R41*t4) + &
           s32*(R032 + R132*t_local + R232*t2) + R02*s2;

!  Compute rho(s,theta,p), first calculating the secant bulk modulus.
    ks_0 = S00 + S10*t_local + S20*t2 + S30*t3 + S40*t4 + &
           s_local*(S01 + S11*t_local + S21*t2 + S31*t3) + s32*(S032 + S132*t_local + S232*t2);
    ks_1 = Sp00 + Sp10*t_local + Sp20*t2 + Sp30*t3 + &
           s_local*(Sp01 + Sp11*t_local + Sp21*t2) + Sp032*s32;
    ks_2 = SP000 + SP010*t_local + SP020*t2 + s_local*(SP001 + SP011*t_local + SP021*t2);

    ks = ks_0 + p1*ks_1 + p2*ks_2;
    dks_dp = ks_1 + 2.0*p1*ks_2;

    rho(j) = rho0*ks / (ks - p1);
! The factor of 1.0e-5 is because pressure here is in bars, not Pa.
    drho_dp(j) = 1.0e-5 * (rho(j) / (ks - p1)) * (1.0 - dks_dp*p1/ks);
  enddo
end subroutine calculate_compress_UNESCO


end module MOM_EOS_UNESCO
