!> Functions to convert between conservative and potential temperature
module MOM_temperature_convert

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

public poTemp_to_consTemp, consTemp_to_poTemp

!>@{ Parameters in the temperature conversion code
real, parameter :: Sprac_Sref = (35.0/35.16504) ! The TEOS 10 conversion factor to go from
                                    ! reference salinity to practical salinity [nondim]
real, parameter :: I_S0 = 0.025*Sprac_Sref ! The inverse of a plausible range of oceanic salinities [kg g-1]
real, parameter :: I_Ts = 0.025          ! The inverse of a plausible range of oceanic temperatures [degC-1]
real, parameter :: I_cp0 = 1.0/3991.86795711963 ! The inverse of the "specific heat" for use
                      !  with Conservative Temperature, as defined with TEOS10 [degC kg J-1]

! The following are coefficients of contributions to conservative temperature as a function of the square root
! of normalized absolute salinity with an offset (zS) and potential temperature (T) with a contribution
! Hab * zS**a * T**b.  The numbers here are copied directly from the corresponding gsw module, but
! the expressions here do not use the same nondimensionalization for pressure or temperature as they do.

real, parameter :: H00 = 61.01362420681071*I_cp0             ! Tp to Tc fit constant [degC]
real, parameter :: H01 = 168776.46138048015*(I_cp0*I_Ts)     ! Tp to Tc fit T coef. [nondim]
real, parameter :: H02 = -2735.2785605119625*(I_cp0*I_Ts**2) ! Tp to Tc fit T**2 coef. [degC-1]
real, parameter :: H03 = 2574.2164453821433*(I_cp0*I_Ts**3)  ! Tp to Tc fit T**3 coef. [degC-2]
real, parameter :: H04 = -1536.6644434977543*(I_cp0*I_Ts**4) ! Tp to Tc fit T**4 coef. [degC-3]
real, parameter :: H05 = 545.7340497931629*(I_cp0*I_Ts**5)   ! Tp to Tc fit T**5 coef. [degC-4]
real, parameter :: H06 = -50.91091728474331*(I_cp0*I_Ts**6)  ! Tp to Tc fit T**6 coef. [degC-5]
real, parameter :: H07 = -18.30489878927802*(I_cp0*I_Ts**7)  ! Tp to Tc fit T**7 coef. [degC-6]
real, parameter :: H20 = 268.5520265845071*I_cp0             ! Tp to Tc fit zS**2 coef. [degC]
real, parameter :: H21 = -12019.028203559312*(I_cp0*I_Ts)    ! Tp to Tc fit zS**2 * T coef. [nondim]
real, parameter :: H22 = 3734.858026725145*(I_cp0*I_Ts**2)   ! Tp to Tc fit zS**2 * T**2 coef. [degC-1]
real, parameter :: H23 = -2046.7671145057618*(I_cp0*I_Ts**3) ! Tp to Tc fit zS**2 * T**3 coef. [degC-2]
real, parameter :: H24 = 465.28655623826234*(I_cp0*I_Ts**4)  ! Tp to Tc fit zS**2 * T**4 coef. [degC-3]
real, parameter :: H25 = -0.6370820302376359*(I_cp0*I_Ts**5) ! Tp to Tc fit zS**2 * T**5 coef. [degC-4]
real, parameter :: H26 = -10.650848542359153*(I_cp0*I_Ts**6) ! Tp to Tc fit zS**2 * T**6 coef. [degC-5]
real, parameter :: H30 = 937.2099110620707*I_cp0             ! Tp to Tc fit zS**3 coef. [degC]
real, parameter :: H31 = 588.1802812170108*(I_cp0*I_Ts)      ! Tp to Tc fit zS** 3* T coef. [nondim]
real, parameter :: H32 = 248.39476522971285*(I_cp0*I_Ts**2)  ! Tp to Tc fit zS**3 * T**2 coef. [degC-1]
real, parameter :: H33 = -3.871557904936333*(I_cp0*I_Ts**3)  ! Tp to Tc fit zS**3 * T**3 coef. [degC-2]
real, parameter :: H34 = -2.6268019854268356*(I_cp0*I_Ts**4) ! Tp to Tc fit zS**3 * T**4 coef. [degC-3]
real, parameter :: H40 = -1687.914374187449*I_cp0            ! Tp to Tc fit zS**4 coef. [degC]
real, parameter :: H41 = 936.3206544460336*(I_cp0*I_Ts)      ! Tp to Tc fit zS**4 * T coef. [nondim]
real, parameter :: H42 = -942.7827304544439*(I_cp0*I_Ts**2)  ! Tp to Tc fit zS**4 * T**2 coef. [degC-1]
real, parameter :: H43 = 369.4389437509002*(I_cp0*I_Ts**3)   ! Tp to Tc fit zS**4 * T**3 coef. [degC-2]
real, parameter :: H44 = -33.83664947895248*(I_cp0*I_Ts**4)  ! Tp to Tc fit zS**4 * T**4 coef. [degC-3]
real, parameter :: H45 = -9.987880382780322*(I_cp0*I_Ts**5)  ! Tp to Tc fit zS**4 * T**5 coef. [degC-4]
real, parameter :: H50 = 246.9598888781377*I_cp0             ! Tp to Tc fit zS**5 coef. [degC]
real, parameter :: H60 = 123.59576582457964*I_cp0            ! Tp to Tc fit zS**6 coef. [degC]
real, parameter :: H70 = -48.5891069025409*I_cp0             ! Tp to Tc fit zS**7 coef. [degC]

!>@}

contains

!> Convert input potential temperature [degC] and absolute salinity [g kg-1] to returned
!! conservative temperature [degC] using the polynomial expressions from TEOS-10.
elemental real function poTemp_to_consTemp(T, Sa) result(Tc)
  real, intent(in) :: T  !< Potential temperature [degC]
  real, intent(in) :: Sa !< Absolute salinity [g kg-1]

  ! Local variables
  real :: x2 ! Absolute salinity normalized by a plausible salinity range [nondim]
  real :: x  ! Square root of normalized absolute salinity [nondim]

  x2 = max(I_S0 * Sa, 0.0)
  x = sqrt(x2)

  Tc = H00 + (T*(H01 + T*(H02 +  T*(H03 +  T*(H04  + T*(H05 + T*(H06 + T* H07)))))) &
                    + x2*(H20 + (T*(H21 +  T*(H22  + T*(H23 + T*(H24 + T*(H25 + T*H26))))) &
                              +  x*(H30 + (T*(H31  + T*(H32 + T*(H33 + T* H34))) &
                                         + x*(H40 + (T*(H41 + T*(H42 + T*(H43 + T*(H44 + T*H45)))) &
                                                  +  x*(H50 + x*(H60 + x* H70)) )) )) )) )

end function poTemp_to_consTemp


!> Return the partial derivative of conservative temperature with potential temperature [nondim]
!! based on the polynomial expressions from TEOS-10.
elemental real function dTc_dTp(T, Sa)
  real, intent(in) :: T  !< Potential temperature [degC]
  real, intent(in) :: Sa !< Absolute salinity [g kg-1]

  ! Local variables
  real :: x2 ! Absolute salinity normalized by a plausible salinity range [nondim]
  real :: x  ! Square root of normalized absolute salinity [nondim]

  x2 = max(I_S0 * Sa, 0.0)
  x = sqrt(x2)

  dTc_dTp = (     H01 + T*(2.*H02 + T*(3.*H03 + T*(4.*H04 + T*(5.*H05 + T*(6.*H06 + T*(7.*H07)))))) ) &
      + x2*(     (H21 + T*(2.*H22 + T*(3.*H23 + T*(4.*H24 + T*(5.*H25 + T*(6.*H26)))))) &
         +  x*(  (H31 + T*(2.*H32 + T*(3.*H33 + T*(4.*H34)))) &
             + x*(H41 + T*(2.*H42 + T*(3.*H43 + T*(4.*H44 + T*(5.*H45))))) ) )

end function dTc_dTp



!> Convert input potential temperature [degC] and absolute salinity [g kg-1] to returned
!! conservative temperature [degC] by inverting the polynomial expressions from TEOS-10.
elemental real function consTemp_to_poTemp(Tc, Sa) result(Tp)
  real, intent(in) :: Tc !< Conservative temperature [degC]
  real, intent(in) :: Sa !< Absolute salinity [g kg-1]

  real :: Tp_num    ! The numerator  of a simple expression for potential temperature [degC]
  real :: I_Tp_den  ! The inverse of the denominator of a simple expression for potential temperature [nondim]
  real :: Tc_diff   ! The difference between an estimate of conservative temperature and its target [degC]
  real :: Tp_old    ! A previous estimate of the potential tempearture [degC]
  real :: dTp_dTc   ! The partial derivative of potential temperature with conservative temperature [nondim]
  ! The following are coefficients in the nominator (TPNxx) or denominator (TPDxx) of a simple rational
  ! expression that approximately converts conservative temperature to potential temperature.
  real, parameter :: TPN00 = -1.446013646344788e-2               ! Simple fit numerator constant [degC]
  real, parameter :: TPN10 = -3.305308995852924e-3*Sprac_Sref    ! Simple fit numerator Sa coef. [degC ppt-1]
  real, parameter :: TPN20 =  1.062415929128982e-4*Sprac_Sref**2 ! Simple fit numerator Sa**2 coef. [degC ppt-2]
  real, parameter :: TPN01 =  9.477566673794488e-1               ! Simple fit numerator Tc coef. [nondim]
  real, parameter :: TPN11 =  2.166591947736613e-3*Sprac_Sref    ! Simple fit numerator Sa * Tc coef. [ppt-1]
  real, parameter :: TPN02 =  3.828842955039902e-3               ! Simple fit numerator Tc**2 coef. [degC-1]
  real, parameter :: TPD10 =  6.506097115635800e-4*Sprac_Sref    ! Simple fit denominator Sa coef. [ppt-1]
  real, parameter :: TPD01 =  3.830289486850898e-3               ! Simple fit denominator Tc coef. [degC-1]
  real, parameter :: TPD02 =  1.247811760368034e-6               ! Simple fit denominator Tc**2 coef. [degC-2]

  ! Estimate the potential temperature and its derivative from an approximate rational function fit.
  Tp_num = TPN00 + (Sa*(TPN10 + TPN20*Sa) + Tc*(TPN01 + (TPN11*Sa + TPN02*Tc)))
  I_Tp_den = 1.0 / (1.0 + (TPD10*Sa + Tc*(TPD01 + TPD02*Tc)))
  Tp = Tp_num*I_Tp_den
  dTp_dTc = ((TPN01 + (TPN11*Sa + 2.*TPN02*Tc)) - (TPD01 + 2.*TPD02*Tc)*Tp)*I_Tp_den

  ! Start the 1.5 iterations through the modified Newton-Raphson iterative method, which is also known
  ! as the Newton-McDougall method.  In this case 1.5 iterations converge to 64-bit machine precision
  ! for oceanographically relevant temperatures and salinities.

  Tc_diff = poTemp_to_consTemp(Tp, Sa) - Tc
  Tp_old = Tp
  Tp = Tp_old - Tc_diff*dTp_dTc

  dTp_dTc = 1.0 / dTc_dTp(0.5*(Tp + Tp_old), Sa)

  Tp = Tp_old - Tc_diff*dTp_dTc
  Tc_diff = poTemp_to_consTemp(Tp, Sa) - Tc
  Tp_old = Tp

  Tp = Tp_old - Tc_diff*dTp_dTc

end function consTemp_to_poTemp

!> \namespace MOM_temperature_conv
!!
!! \section MOM_temperature_conv Temperature conversions
!!
!!   This module has functions that convert potential temperature to conservative temperature
!! and the reverse, as described in the TEOS-10 manual.  This code was originally derived
!! from their corresponding routines in the gsw code package, but has had some refactoring so that the
!! answers are more likely to reproduce across compilers and levels of optimization.  A complete
!! discussion of the thermodynamics of seawater and the definition of conservative temperature
!! can be found in IOC et al. (2010).
!!
!! \subsection section_temperature_conv_references References
!!
!! IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010:
!!   Calculation and use of thermodynamic properties. Intergovernmental Oceanographic Commission,
!!   Manuals and Guides No. 56, UNESCO (English), 196 pp.
!!   (Available from www.teos-10.org/pubs/TEOS-10_Manual.pdf)

end module MOM_temperature_convert
