module MOM_ocmip2_co2calc_mod  !{
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
!
!
!<CONTACT EMAIL="Richard.Slater@noaa.gov"> Richard D. Slater
!</CONTACT>
!
!<REVIEWER EMAIL="John.Dunne@noaa.gov"> John P. Dunne
!</REVIEWER>
!
!<OVERVIEW>
! Surface fCO2 calculation
!</OVERVIEW>
!
!<DESCRIPTION>
! Calculate the fugacity of CO2 at the surface in thermodynamic
! equilibrium with the current alkalinity (Alk) and total dissolved
! inorganic carbon (DIC) at a particular temperature and salinity
! using an initial guess for the total hydrogen
! ion concentration (htotal)
!</DESCRIPTION>
!

!
!------------------------------------------------------------------
!
!       Global definitions
!
!------------------------------------------------------------------
!

implicit none

private

public  :: MOM_ocmip2_co2calc, CO2_dope_vector

! This include declares and sets the variable "version".
#include "version_variable.h"

type CO2_dope_vector
  integer  :: isc, iec, jsc, jec
  integer  :: isd, ied, jsd, jed
end type CO2_dope_vector
!
!-----------------------------------------------------------------------
!
!       Subroutine and function definitions
!
!-----------------------------------------------------------------------
!

contains


!#######################################################################
! <SUBROUTINE NAME="MOM_ocmip2_co2calc">
!
! <DESCRIPTION>
!       Calculate co2* from total alkalinity and total CO2 at
! temperature (t) and salinity (s).
! It is assumed that init_ocmip2_co2calc has already been called with
! the T and S to calculate the various coefficients.
!
! INPUT
!
!       dope_vec   = an array of indices corresponding to the compute
!                    and data domain boundaries.
!
!       mask       = land mask array (0.0 = land)
!
!       dic_in     = total inorganic carbon (mol/kg)
!                    where 1 T = 1 metric ton = 1000 kg
!
!       ta_in      = total alkalinity (eq/kg)
!
!       pt_in      = inorganic phosphate (mol/kg)
!
!       sit_in     = inorganic silicate (mol/kg)
!
!       htotallo   = lower limit of htotal range
!
!       htotalhi   = upper limit of htotal range
!
!       htotal     = H+ concentration (mol/kg)
!
! OUTPUT
!       co2star    = CO2*water, or H2CO3 concentration (mol/kg)
!       alpha      = Solubility of CO2 for air (mol/kg/atm)
!       pco2surf   = oceanic pCO2 (ppmv)
!       co3_ion    = Carbonate ion, or CO3-- concentration (mol/kg)
!
! FILES and PROGRAMS NEEDED: drtsafe, ta_iter_1
!
! IMPORTANT: co2star and alpha need to be multiplied by rho before being
! passed to the atmosphere.
!
! </DESCRIPTION>

subroutine MOM_ocmip2_co2calc(dope_vec, mask,                      &
                          t_in, s_in, dic_in, pt_in, sit_in, ta_in, htotallo, &
                          htotalhi, htotal, co2star, alpha, pCO2surf, co3_ion)  !{

implicit none

!
!       local parameters
!

real, parameter :: permeg = 1.e-6
real, parameter :: xacc = 1.0e-10

!
!       arguments
!
type(CO2_dope_vector), intent(in)           :: dope_vec
real, dimension(dope_vec%isd:dope_vec%ied,dope_vec%jsd:dope_vec%jed), &
      intent(in)::             mask, &
                               t_in, &
                               s_in, &
                               dic_in, &
                               pt_in, &
                               sit_in, &
                               ta_in, &
                               htotallo, &
                               htotalhi
real, dimension(dope_vec%isd:dope_vec%ied,dope_vec%jsd:dope_vec%jed), &
      intent(inout)         :: htotal
real, dimension(dope_vec%isd:dope_vec%ied,dope_vec%jsd:dope_vec%jed), &
      intent(out), optional :: alpha, &
                               pCO2surf, &
                               co2star, &
                               co3_ion
!
!       local variables
!
integer :: isc, iec, jsc, jec
integer :: i,j
real :: alpha_internal
real :: bt
real :: co2star_internal
real :: dlogtk
real :: ft
real :: htotal2
real :: invtk
real :: is
real :: is2
real :: k0
real :: k1
real :: k2
real :: k1p
real :: k2p
real :: k3p
real :: kb
real :: kf
real :: ks
real :: ksi
real :: kw
real :: log100
real :: s2
real :: scl
real :: sqrtis
real :: sqrts
real :: st
real :: tk
real :: tk100
real :: tk1002
real :: logf_of_s

! Set the loop indices.
  isc = dope_vec%isc ; iec = dope_vec%iec
  jsc = dope_vec%jsc ; jec = dope_vec%jec

!
!       Initialize the module
!
  log100 = log(100.0)

  do j = jsc, jec  !{
    do i = isc, iec  !{
!
!---------------------------------------------------------------------
!
!***********************************************************************
! Calculate all constants needed to convert between various measured
! carbon species. References for each equation are noted in the code.
! Once calculated, the constants are stored and passed in the common
! block "const". The original version of this code was based on
! the code by Dickson in Version 2 of "Handbook of Methods for the
! Analysis of the Various Parameters of the Carbon Dioxide System
! in Seawater", DOE, 1994 (SOP No. 3, p25-26).
        tk        = 273.15 + t_in(i,j)
        tk100     = tk / 100.0
        tk1002     = tk100**2
        invtk     = 1.0 / tk
        dlogtk    = log(tk)
        is        = 19.924 * s_in(i,j) /(1000.0 -1.005 * s_in(i,j))
        is2       = is * is
        sqrtis    = sqrt(is)
        s2        = s_in(i,j) * s_in(i,j)
        sqrts     = sqrt(s_in(i,j))
        scl       = s_in(i,j) / 1.80655
        logf_of_s = log(1.0 - 0.001005 * s_in(i,j))
!
! k0 from Weiss 1974
!

        k0 = exp(93.4517/tk100 - 60.2409 + 23.3585 * log(tk100) +       &
                 s_in(i,j) * (0.023517 - 0.023656 * tk100 +             &
                 0.0047036 * tk1002))
!
! k1 = [H][HCO3]/[H2CO3]
! k2 = [H][CO3]/[HCO3]
!
! Millero p.664 (1995) using Mehrbach et al. data on seawater scale
!

        k1 = 10.0**(-(3670.7 * invtk - 62.008 + 9.7944 * dlogtk -       &
                             0.0118 * s_in(i,j) + 0.000116 * s2))
        k2 = 10.0**(-(1394.7 * invtk + 4.777 -                          &
                             0.0184 * s_in(i,j) + 0.000118 * s2))
!
! kb = [H][BO2]/[HBO2]
!
! Millero p.669 (1995) using data from Dickson (1990)
!

        kb = exp((-8966.90 - 2890.53 * sqrts - 77.942 * s_in(i,j) +     &
                 1.728 * sqrts**3 - 0.0996 * s2) * invtk + (148.0248 +  &
                 137.1942 * sqrts + 1.62142 * s_in(i,j)) + (-24.4344 -  &
                 25.085 * sqrts - 0.2474 * s_in(i,j)) * dlogtk +        &
                 0.053105 * sqrts * tk)
!
! k1p = [H][H2PO4]/[H3PO4]
!
! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
!

        k1p = exp(-4576.752 * invtk + 115.525 - 18.453 * dlogtk +       &
                  (-106.736 * invtk + 0.69171) * sqrts + (-0.65643 *    &
                  invtk - 0.01844) * s_in(i,j))
!
! k2p = [H][HPO4]/[H2PO4]
!
! DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
!

        k2p = exp(-8814.715 * invtk + 172.0883 - 27.927 * (-160.340 *   &
                  invtk + 1.3566) * sqrts + (0.37335 * invtk -          &
                  0.05778) * s_in(i,j))
!
!-----------------------------------------------------------------------
! k3p = [H][PO4]/[HPO4]
!
! DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
!

        k3p = exp(-3070.75 * invtk - 18.141 +(17.27039 * invtk +        &
                  2.81197) * sqrts + (-44.99486 * invtk - 0.09984) *    &
                  s_in(i,j))
!
!-----------------------------------------------------------------------
! ksi = [H][SiO(OH)3]/[Si(OH)4]
!
! Millero p.671 (1995) using data from Yao and Millero (1995)
!
        ksi = exp(-8904.2 * invtk + 117.385 - 19.334 * dlogtk +         &
                  (-458.79 * invtk + 3.5913) * sqrtis + (188.74 *       &
                  invtk - 1.5998) * is + (-12.1652 * invtk + 0.07871) * &
                  is2 + logf_of_s)
!
!-----------------------------------------------------------------------
! kw = [H][OH]
!
! Millero p.670 (1995) using composite data
!

        kw = exp(-13847.26 * invtk + 148.9652 - 23.6521 *  dlogtk +     &
                 (118.67 * invtk - 5.977 + 1.0495 * dlogtk) * sqrts -   &
                 0.01615 * s_in(i,j))
!
!-----------------------------------------------------------------------
! ks = [H][SO4]/[HSO4]
!
! Dickson (1990, J. chem. Thermodynamics 22, 113)
!
        ks = exp(-4276.1 * invtk + 141.328 - 23.093 * dlogtk +          &
                 (-13856.0 * invtk + 324.57 - 47.986 * dlogtk) *        &
                 sqrtis + (35474.0 * invtk - 771.54 + 114.723 *         &
                 dlogtk) * is - 2698.0 * invtk * sqrtis**3 +            &
                 1776.0 * invtk * is2 + logf_of_s)
!
!-----------------------------------------------------------------------
! kf = [H][F]/[HF]
!
! Dickson and Riley (1979) -- change pH scale to total
!
        kf = exp(1590.2 * invtk - 12.641 + 1.525 * sqrtis + logf_of_s + &
                 log(1.0 + (0.1400 / 96.062) * scl / ks))
!
!-----------------------------------------------------------------------
! Calculate concentrations for borate, sulfate, and fluoride
!
! Uppstrom (1974)
!
        bt = 0.000232 / 10.811 * scl
!
! Morris & Riley (1966)
!
        st = 0.14 / 96.062 * scl
!
! Riley (1965)
!
        ft = 0.000067 / 18.9984 * scl
!
!***********************************************************************
!
! Calculate [H+] total when DIC and TA are known at T, S and 1 atm.
! The solution converges to err of xacc. The solution must be within
! the range x1 to x2.
!
! If DIC and TA are known then either a root finding or iterative method
! must be used to calculate htotal. In this case we use the
! Newton-Raphson "safe" method taken from "Numerical Recipes"
! (function "rtsafe.f" with error trapping removed).
!
! As currently set, this procedure iterates about 12 times. The x1
! and x2 values set below will accomodate ANY oceanographic values.
! If an initial guess of the pH is known, then the number of
! iterations can be reduced to about 5 by narrowing the gap between
! x1 and x2. It is recommended that the first few time steps be run
! with x1 and x2 set as below. After that, set x1 and x2 to the
! previous value of the pH +/- ~0.5. The current setting of xacc will
! result in co2star accurate to 3 significant figures (xx.y). Making
! xacc bigger will result in faster convergence also, but this is not
! recommended (xacc of 10**-9 drops precision to 2 significant
! figures).
!
      if (mask(i,j) .ne. 0.0) then  !{
        htotal(i,j) = drtsafe(k0, k1, k2, kb, k1p, k2p, k3p, ksi, kw,   &
                                ks, kf, bt, dic_in(i,j), ft, pt_in(i,j),&
                                sit_in(i,j), st, ta_in(i,j),            &
                                htotalhi(i,j), htotallo(i,j), xacc)
      endif
!
! Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2,
! ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
!
        htotal2 = htotal(i,j) * htotal(i,j)
        co2star_internal = dic_in(i,j) * htotal2 / (htotal2 +           &
                         k1 * htotal(i,j) + k1 * k2)
        if (present(co2star)) co2star(i,j) = co2star_internal
        if (present(co3_ion)) co3_ion(i,j) = co2star_internal * k1 * k2 / htotal2
!
! Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6
!                values)
!
        if (present(alpha) .or. present(pCO2surf)) then
          alpha_internal = exp(-162.8301 + 218.2968 / tk100 + 90.9241 *     &
                         (dlogtk -log100) - 1.47696 * tk1002 +          &
                         s_in(i,j) * (0.025695 - 0.025225 * tk100 +   &
                         0.0049867 * tk1002))
        endif
        if (present(alpha)) alpha(i,j) = alpha_internal
        if (present(pCO2surf)) then
          pCO2surf(i,j) = co2star_internal / (alpha_internal * permeg)
        endif
    enddo  !} i
  enddo  !} j

return

end subroutine  MOM_ocmip2_co2calc  !}
! </SUBROUTINE> NAME="MOM_ocmip2_co2calc"


!#######################################################################
! <FUNCTION NAME="drtsafe">
!
! <DESCRIPTION>
!       File taken from Numerical Recipes. Modified  R. M. Key 4/94
! </DESCRIPTION>

function drtsafe(k0, k1, k2, kb, k1p, k2p, k3p, ksi, kw, ks, kf, &
               bt, dic, ft, pt, sit, st, ta, x1, x2, xacc)  !{

implicit none

!
!       arguments
!

real    :: k0, k1, k2, kb, k1p, k2p, k3p, ksi, kw, ks, kf
real    :: bt, dic, ft, pt, sit, st, ta
real    :: drtsafe
real    :: x1, x2, xacc

!
!       local parameters
!

integer, parameter      :: maxit = 100

!
!       local variables
!

integer :: j
real    :: fl, df, fh, swap, xl, xh, dxold, dx, f, temp

call ta_iter_1(k0, k1, k2, kb, k1p, k2p, k3p, ksi, kw, ks, kf, &
               bt, dic, ft, pt, sit, st, ta, x1, fl, df)
call ta_iter_1(k0, k1, k2, kb, k1p, k2p, k3p, ksi, kw, ks, kf, &
               bt, dic, ft, pt, sit, st, ta, x2, fh, df)
if(fl .lt. 0.0) then
  xl=x1
  xh=x2
else
  xh=x1
  xl=x2
  swap=fl
  fl=fh
  fh=swap
end if
drtsafe=0.5*(x1+x2)
dxold=abs(x2-x1)
dx=dxold
call ta_iter_1(k0, k1, k2, kb, k1p, k2p, k3p, ksi, kw, ks, kf, &
               bt, dic, ft, pt, sit, st, ta, drtsafe, f, df)
do j=1,maxit  !{
  if (((drtsafe-xh)*df-f)*((drtsafe-xl)*df-f) .ge. 0.0 .or.     &
      abs(2.0*f) .gt. abs(dxold*df)) then
    dxold=dx
    dx=0.5*(xh-xl)
    drtsafe=xl+dx
    if (xl .eq. drtsafe) then
!     write (6,*) 'Exiting drtsafe at A on iteration  ', j, ', ph = ', -log10(drtsafe)
      return
    endif
  else
    dxold=dx
    dx=f/df
    temp=drtsafe
    drtsafe=drtsafe-dx
    if (temp .eq. drtsafe) then
!     write (6,*) 'Exiting drtsafe at B on iteration  ', j, ', ph = ', -log10(drtsafe)
      return
    endif
  end if
  if (abs(dx) .lt. xacc) then
!     write (6,*) 'Exiting drtsafe at C on iteration  ', j, ', ph = ', -log10(drtsafe)
    return
  endif
  call ta_iter_1(k0, k1, k2, kb, k1p, k2p, k3p, ksi, kw, ks, kf, &
                 bt, dic, ft, pt, sit, st, ta, drtsafe, f, df)
  if(f .lt. 0.0) then
    xl=drtsafe
    fl=f
  else
    xh=drtsafe
    fh=f
  end if
enddo  !} j

return

end  function  drtsafe  !}
! </FUNCTION> NAME="drtsafe"


!#######################################################################
! <SUBROUTINE NAME="ta_iter_1">
!
! <DESCRIPTION>
! This routine expresses TA as a function of DIC, htotal and constants.
! It also calculates the derivative of this function with respect to
! htotal. It is used in the iterative solution for htotal. In the call
! "x" is the input value for htotal, "fn" is the calculated value for TA
! and "df" is the value for dTA/dhtotal
! </DESCRIPTION>

subroutine ta_iter_1(k0, k1, k2, kb, k1p, k2p, k3p, ksi, kw, ks, kf, &
                     bt, dic, ft, pt, sit, st, ta, x, fn, df)  !{

implicit none

!
!       arguments
!

real    :: k0, k1, k2, kb, k1p, k2p, k3p, ksi, kw, ks, kf
real    :: bt, dic, ft, pt, sit, st, ta, x, fn, df

!
!       local variables
!

real    :: x2, x3, k12, k12p, k123p, c, a, a2, da, b, b2, db

x2 = x*x
x3 = x2*x
k12 = k1*k2
k12p = k1p*k2p
k123p = k12p*k3p
c = 1.0 + st/ks
a = x3 + k1p*x2 + k12p*x + k123p
a2 = a*a
da = 3.0*x2 + 2.0*k1p*x + k12p
b = x2 + k1*x + k12
b2 = b*b
db = 2.0*x + k1
!
!     fn = hco3+co3+borate+oh+hpo4+2*po4+silicate+hfree+hso4+hf+h3po4-ta
!
fn = k1*x*dic/b + 2.0*dic*k12/b + bt/ (1.0 + x/kb) + kw/x +   &
     pt*k12p*x/a + 2.0*pt*k123p/a + sit/(1.0 + x/ksi) -       &
     x/c - st/(1.0 + ks/x/c) - ft/(1.0 + kf/x) - pt*x3/a - ta
!
!     df = dfn/dx
!
df = ((k1*dic*b) - k1*x*dic*db)/b2 - 2.0*dic*k12*db/b2 -      &
     bt/kb/(1.0+x/kb)**2 - kw/x2 + (pt*k12p*(a - x*da))/a2 -  &
     2.0*pt*k123p*da/a2 - sit/ksi/ (1.0+x/ksi)**2 - 1.0/c +   &
     st*(1.0 + ks/x/c)**(-2)*(ks/c/x2) +                      &
     ft*(1.0 + kf/x)**(-2)*kf/x2 - pt*x2*(3.0*a-x*da)/a2

return

end subroutine  ta_iter_1  !}
! </SUBROUTINE> NAME="ta_iter_1"

end module  MOM_ocmip2_co2calc_mod  !}
