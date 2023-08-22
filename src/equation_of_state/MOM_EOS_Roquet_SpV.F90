!> The equation of state for specific volume (SpV) using the expressions of Roquet et al. 2015
module MOM_EOS_Roquet_Spv

! This file is part of MOM6. See LICENSE.md for the license.

!use gsw_mod_toolbox, only : gsw_sr_from_sp, gsw_ct_from_pt

implicit none ; private

public calculate_compress_Roquet_SpV, calculate_density_Roquet_SpV, calculate_spec_vol_Roquet_SpV
public calculate_density_derivs_Roquet_SpV, calculate_specvol_derivs_Roquet_SpV
public calculate_density_scalar_Roquet_SpV, calculate_density_array_Roquet_SpV
public calculate_density_second_derivs_Roquet_SpV, EoS_fit_range_Roquet_SpV

!> Compute the in situ density of sea water [kg m-3], or its anomaly with respect to
!! a reference density, from absolute salinity [g kg-1], conservative temperature [degC],
!! and pressure [Pa], using the specific volume polynomial fit from Roquet et al. (2015)
interface calculate_density_Roquet_SpV
  module procedure calculate_density_scalar_Roquet_SpV, calculate_density_array_Roquet_SpV
end interface calculate_density_Roquet_SpV

!> Compute the in situ specific volume of sea water (in [m3 kg-1]), or an anomaly with respect
!! to a reference specific volume, from absolute salinity ([g kg-1]), conservative
!! temperature (in degrees Celsius [degC]), and pressure [Pa], using the specific volume
!! polynomial fit from Roquet et al. (2015)
interface calculate_spec_vol_Roquet_SpV
  module procedure calculate_spec_vol_scalar_Roquet_SpV, calculate_spec_vol_array_Roquet_SpV
end interface calculate_spec_vol_Roquet_SpV

!> For a given thermodynamic state, return the derivatives of density with conservative temperature
!! and absolute salinity, using the specific volume polynomial fit from Roquet et al. (2015)
interface calculate_density_derivs_Roquet_SpV
  module procedure calculate_density_derivs_scalar_Roquet_SpV, calculate_density_derivs_array_Roquet_SpV
end interface calculate_density_derivs_Roquet_SpV

!> Compute the second derivatives of density with various combinations of temperature, salinity
!! and pressure using the specific volume polynomial fit from Roquet et al. (2015)
interface calculate_density_second_derivs_Roquet_SpV
  module procedure calculate_density_second_derivs_scalar_Roquet_SpV
  module procedure calculate_density_second_derivs_array_Roquet_SpV
end interface calculate_density_second_derivs_Roquet_SpV

real, parameter :: Pa2kb  = 1.e-8 !< Conversion factor between Pa and kbar [kbar Pa-1]
!>@{ Parameters in the Roquet specific volume polynomial equation of state
real, parameter :: rdeltaS = 24.          ! An offset to salinity before taking its square root [g kg-1]
real, parameter :: r1_S0 = 0.875/35.16504 ! The inverse of a plausible range of oceanic salinities [kg g-1]
real, parameter :: I_Ts = 0.025           ! The inverse of a plausible range of oceanic temperatures [degC-1]
! The following are the coefficients of the fit to the reference density profile (rho00p) as a function of
! pressure (P), with a contribution R0c * P**(c+1).  The nomenclature follows Roquet.
real, parameter :: V00 = -4.4015007269e-05*Pa2kb    ! SpV00p P coef.    [m3 kg-1 Pa-1]
real, parameter :: V01 = 6.9232335784e-06*Pa2kb**2  ! SpV00p P**2 coef. [m3 kg-1 Pa-2]
real, parameter :: V02 = -7.5004675975e-07*Pa2kb**3 ! SpV00p P**3 coef. [m3 kg-1 Pa-3]
real, parameter :: V03 = 1.7009109288e-08*Pa2kb**4  ! SpV00p P**4 coef. [m3 kg-1 Pa-4]
real, parameter :: V04 = -1.6884162004e-08*Pa2kb**5 ! SpV00p P**5 coef. [m3 kg-1 Pa-5]
real, parameter :: V05 = 1.9613503930e-09*Pa2kb**6  ! SpV00p P**6 coef. [m3 kg-1 Pa-6]

! The following terms are contributions to specific volume (SpV) as a function of the square root of
! normalized absolute salinity with an offset (zs), temperature (T) and pressure (P), with a contribution
! SPVabc * zs**a * T**b * P**c.  The numbers here are copied directly from Roquet et al. (2015), but
! the expressions here do not use the same nondimensionalization for pressure or temperature as they do.
real, parameter :: SPV000 = 1.0772899069e-03                  ! Constant SpV contribution  [m3 kg-1]
real, parameter :: SPV100 = -3.1263658781e-04                 ! SpV zs coef.               [m3 kg-1]
real, parameter :: SPV200 = 6.7615860683e-04                  ! SpV zs**2 coef.            [m3 kg-1]
real, parameter :: SPV300 = -8.6127884515e-04                 ! SpV zs**3 coef.            [m3 kg-1]
real, parameter :: SPV400 = 5.9010812596e-04                  ! SpV zs**4 coef.            [m3 kg-1]
real, parameter :: SPV500 = -2.1503943538e-04                 ! SpV zs**5 coef.            [m3 kg-1]
real, parameter :: SPV600 = 3.2678954455e-05                  ! SpV zs**6 coef.            [m3 kg-1]
real, parameter :: SPV010 = -1.4949652640e-05*I_Ts            ! SpV T coef.         [m3 kg-1 degC-1]
real, parameter :: SPV110 = 3.1866349188e-05*I_Ts             ! SpV zs * T coef.    [m3 kg-1 degC-1]
real, parameter :: SPV210 = -3.8070687610e-05*I_Ts            ! SpV zs**2 * T coef. [m3 kg-1 degC-1]
real, parameter :: SPV310 = 2.9818473563e-05*I_Ts             ! SpV zs**3 * T coef. [m3 kg-1 degC-1]
real, parameter :: SPV410 = -1.0011321965e-05*I_Ts            ! SpV zs**4 * T coef. [m3 kg-1 degC-1]
real, parameter :: SPV510 = 1.0751931163e-06*I_Ts             ! SpV zs**5 * T coef. [m3 kg-1 degC-1]
real, parameter :: SPV020 = 2.7546851539e-05*I_Ts**2          ! SpV T**2 coef.      [m3 kg-1 degC-2]
real, parameter :: SPV120 = -3.6597334199e-05*I_Ts**2         ! SpV zs * T**2 coef. [m3 kg-1 degC-2]
real, parameter :: SPV220 = 3.4489154625e-05*I_Ts**2          ! SpV zs**2 * T**2 coef. [m3 kg-1 degC-2]
real, parameter :: SPV320 = -1.7663254122e-05*I_Ts**2         ! SpV zs**3 * T**2 coef. [m3 kg-1 degC-2]
real, parameter :: SPV420 = 3.5965131935e-06*I_Ts**2          ! SpV zs**4 * T**2 coef. [m3 kg-1 degC-2]
real, parameter :: SPV030 = -1.6506828994e-05*I_Ts**3         ! SpV T**3 coef.      [m3 kg-1 degC-3]
real, parameter :: SPV130 = 2.4412359055e-05*I_Ts**3          ! SpV zs * T**3 coef. [m3 kg-1 degC-3]
real, parameter :: SPV230 = -1.4606740723e-05*I_Ts**3         ! SpV zs**2 * T**3 coef. [m3 kg-1 degC-3]
real, parameter :: SPV330 = 2.3293406656e-06*I_Ts**3          ! SpV zs**3 * T**3 coef. [m3 kg-1 degC-3]
real, parameter :: SPV040 = 6.7896174634e-06*I_Ts**4          ! SpV T**4 coef.      [m3 kg-1 degC-4]
real, parameter :: SPV140 = -8.7951832993e-06*I_Ts**4         ! SpV zs * T**4 coef. [m3 kg-1 degC-4]
real, parameter :: SPV240 = 4.4249040774e-06*I_Ts**4          ! SpV zs**2 * T**4 coef. [m3 kg-1 degC-4]
real, parameter :: SPV050 = -7.2535743349e-07*I_Ts**5         ! SpV T**5 coef.      [m3 kg-1 degC-5]
real, parameter :: SPV150 = -3.4680559205e-07*I_Ts**5         ! SpV zs * T**5 coef. [m3 kg-1 degC-5]
real, parameter :: SPV060 = 1.9041365570e-07*I_Ts**6          ! SpV T**6 coef.      [m3 kg-1 degC-6]
real, parameter :: SPV001 = -1.6889436589e-05*Pa2kb           ! SpV P coef.           [m3 kg-1 Pa-1]
real, parameter :: SPV101 = 2.1106556158e-05*Pa2kb            ! SpV zs * P coef.      [m3 kg-1 Pa-1]
real, parameter :: SPV201 = -2.1322804368e-05*Pa2kb           ! SpV zs**2 * P coef.   [m3 kg-1 Pa-1]
real, parameter :: SPV301 = 1.7347655458e-05*Pa2kb            ! SpV zs**3 * P coef.   [m3 kg-1 Pa-1]
real, parameter :: SPV401 = -4.3209400767e-06*Pa2kb           ! SpV zs**4 * P coef.   [m3 kg-1 Pa-1]
real, parameter :: SPV011 = 1.5355844621e-05*(I_Ts*Pa2kb)     ! SpV T * P coef. [m3 kg-1 degC-1 Pa-1]
real, parameter :: SPV111 = 2.0914122241e-06*(I_Ts*Pa2kb)     ! SpV zs * T * P coef. [m3 kg-1 degC-1 Pa-1]
real, parameter :: SPV211 = -5.7751479725e-06*(I_Ts*Pa2kb)    ! SpV zs**2 * T * P coef. [m3 kg-1 degC-1 Pa-1]
real, parameter :: SPV311 = 1.0767234341e-06*(I_Ts*Pa2kb)     ! SpV zs**3 * T * P coef. [m3 kg-1 degC-1 Pa-1]
real, parameter :: SPV021 = -9.6659393016e-06*(I_Ts**2*Pa2kb) ! SpV T**2 * P coef. [m3 kg-1 degC-2 Pa-1]
real, parameter :: SPV121 = -7.0686982208e-07*(I_Ts**2*Pa2kb) ! SpV zs * T**2 * P coef. [m3 kg-1 degC-2 Pa-1]
real, parameter :: SPV221 = 1.4488066593e-06*(I_Ts**2*Pa2kb)  ! SpV zs**2 * T**2 * P coef. [m3 kg-1 degC-2 Pa-1]
real, parameter :: SPV031 = 3.1134283336e-06*(I_Ts**3*Pa2kb)  ! SpV T**3 * P coef. [m3 kg-1 degC-3 Pa-1]
real, parameter :: SPV131 = 7.9562529879e-08*(I_Ts**3*Pa2kb)  ! SpV zs * T**3 * P coef. [m3 kg-1 degC-3 Pa-1]
real, parameter :: SPV041 = -5.6590253863e-07*(I_Ts**4*Pa2kb) ! SpV T**4 * P coef. [m3 kg-1 degC-4 Pa-1]
real, parameter :: SPV002 = 1.0500241168e-06*Pa2kb**2         ! SpV P**2 coef.        [m3 kg-1 Pa-2]
real, parameter :: SPV102 = 1.9600661704e-06*Pa2kb**2         ! SpV zs * P**2 coef.   [m3 kg-1 Pa-2]
real, parameter :: SPV202 = -2.1666693382e-06*Pa2kb**2        ! SpV zs**2 * P**2 coef. [m3 kg-1 Pa-2]
real, parameter :: SPV012 = -3.8541359685e-06*(I_Ts*Pa2kb**2) ! SpV T * P**2 coef. [m3 kg-1 degC-1 Pa-2]
real, parameter :: SPV112 = 1.0157632247e-06*(I_Ts*Pa2kb**2)  ! SpV zs * T * P**2 coef. [m3 kg-1 degC-1 Pa-2]
real, parameter :: SPV022 = 1.7178343158e-06*(I_Ts**2*Pa2kb**2) ! SpV T**2 * P**2 coef. [m3 kg-1 degC-2 Pa-2]
real, parameter :: SPV003 = -4.1503454190e-07*Pa2kb**3        ! SpV P**3 coef.        [m3 kg-1 Pa-3]
real, parameter :: SPV103 = 3.5627020989e-07*Pa2kb**3         ! SpV zs * P**3 coef.   [m3 kg-1 Pa-3]
real, parameter :: SPV013 = -1.1293871415e-07*(I_Ts*Pa2kb**3) ! SpV T * P**3 coef. [m3 kg-1 degC-1 Pa-3]

real, parameter :: ALP000 =    SPV010   ! Constant in the dSpV_dT fit               [m3 kg-1 degC-1]
real, parameter :: ALP100 =    SPV110   ! dSpV_dT fit zs coef.                      [m3 kg-1 degC-1]
real, parameter :: ALP200 =    SPV210   ! dSpV_dT fit zs**2 coef.                   [m3 kg-1 degC-1]
real, parameter :: ALP300 =    SPV310   ! dSpV_dT fit zs**3 coef.                   [m3 kg-1 degC-1]
real, parameter :: ALP400 =    SPV410   ! dSpV_dT fit zs**4 coef.                   [m3 kg-1 degC-1]
real, parameter :: ALP500 =    SPV510   ! dSpV_dT fit zs**5 coef.                   [m3 kg-1 degC-1]
real, parameter :: ALP010 = 2.*SPV020   ! dSpV_dT fit T coef.                       [m3 kg-1 degC-2]
real, parameter :: ALP110 = 2.*SPV120   ! dSpV_dT fit zs * T coef.                  [m3 kg-1 degC-2]
real, parameter :: ALP210 = 2.*SPV220   ! dSpV_dT fit zs**2 * T coef.               [m3 kg-1 degC-2]
real, parameter :: ALP310 = 2.*SPV320   ! dSpV_dT fit zs**3 * T coef.               [m3 kg-1 degC-2]
real, parameter :: ALP410 = 2.*SPV420   ! dSpV_dT fit zs**4 * T coef.               [m3 kg-1 degC-2]
real, parameter :: ALP020 = 3.*SPV030   ! dSpV_dT fit T**2 coef.                    [m3 kg-1 degC-3]
real, parameter :: ALP120 = 3.*SPV130   ! dSpV_dT fit zs * T**2 coef.               [m3 kg-1 degC-3]
real, parameter :: ALP220 = 3.*SPV230   ! dSpV_dT fit zs**2 * T**2 coef.            [m3 kg-1 degC-3]
real, parameter :: ALP320 = 3.*SPV330   ! dSpV_dT fit zs**3 * T**2 coef.            [m3 kg-1 degC-3]
real, parameter :: ALP030 = 4.*SPV040   ! dSpV_dT fit T**3 coef.                    [m3 kg-1 degC-4]
real, parameter :: ALP130 = 4.*SPV140   ! dSpV_dT fit zs * T**3 coef.               [m3 kg-1 degC-4]
real, parameter :: ALP230 = 4.*SPV240   ! dSpV_dT fit zs**2 * T**3 coef.            [m3 kg-1 degC-4]
real, parameter :: ALP040 = 5.*SPV050   ! dSpV_dT fit T**4 coef.                    [m3 kg-1 degC-5]
real, parameter :: ALP140 = 5.*SPV150   ! dSpV_dT fit zs* * T**4 coef.              [m3 kg-1 degC-5]
real, parameter :: ALP050 = 6.*SPV060   ! dSpV_dT fit T**5 coef.                    [m3 kg-1 degC-6]
real, parameter :: ALP001 =    SPV011   ! dSpV_dT fit P coef.                  [m3 kg-1 degC-1 Pa-1]
real, parameter :: ALP101 =    SPV111   ! dSpV_dT fit zs * P coef.             [m3 kg-1 degC-1 Pa-1]
real, parameter :: ALP201 =    SPV211   ! dSpV_dT fit zs**2 * P coef.          [m3 kg-1 degC-1 Pa-1]
real, parameter :: ALP301 =    SPV311   ! dSpV_dT fit zs**3 * P coef.          [m3 kg-1 degC-1 Pa-1]
real, parameter :: ALP011 = 2.*SPV021   ! dSpV_dT fit T * P coef.              [m3 kg-1 degC-2 Pa-1]
real, parameter :: ALP111 = 2.*SPV121   ! dSpV_dT fit zs * T * P coef.         [m3 kg-1 degC-2 Pa-1]
real, parameter :: ALP211 = 2.*SPV221   ! dSpV_dT fit zs**2 * T * P coef.      [m3 kg-1 degC-2 Pa-1]
real, parameter :: ALP021 = 3.*SPV031   ! dSpV_dT fit T**2 * P coef.           [m3 kg-1 degC-3 Pa-1]
real, parameter :: ALP121 = 3.*SPV131   ! dSpV_dT fit zs * T**2 * P coef.      [m3 kg-1 degC-3 Pa-1]
real, parameter :: ALP031 = 4.*SPV041   ! dSpV_dT fit T**3 * P coef.           [m3 kg-1 degC-4 Pa-1]
real, parameter :: ALP002 =    SPV012   ! dSpV_dT fit P**2 coef.               [m3 kg-1 degC-1 Pa-2]
real, parameter :: ALP102 =    SPV112   ! dSpV_dT fit zs * P**2 coef.          [m3 kg-1 degC-1 Pa-2]
real, parameter :: ALP012 = 2.*SPV022   ! dSpV_dT fit T * P**2 coef.           [m3 kg-1 degC-2 Pa-2]
real, parameter :: ALP003 =    SPV013   ! dSpV_dT fit P**3 coef.               [m3 kg-1 degC-1 Pa-3]

real, parameter :: BET000 = 0.5*SPV100*r1_S0  ! Constant in the dSpV_dS fit          [m3 kg-1 ppt-1]
real, parameter :: BET100 =     SPV200*r1_S0  ! dSpV_dS fit zs coef.                 [m3 kg-1 ppt-1]
real, parameter :: BET200 = 1.5*SPV300*r1_S0  ! dSpV_dS fit zs**2 coef.              [m3 kg-1 ppt-1]
real, parameter :: BET300 = 2.0*SPV400*r1_S0  ! dSpV_dS fit zs**3 coef.              [m3 kg-1 ppt-1]
real, parameter :: BET400 = 2.5*SPV500*r1_S0  ! dSpV_dS fit zs**4 coef.              [m3 kg-1 ppt-1]
real, parameter :: BET500 = 3.0*SPV600*r1_S0  ! dSpV_dS fit zs**5 coef.              [m3 kg-1 ppt-1]
real, parameter :: BET010 = 0.5*SPV110*r1_S0  ! dSpV_dS fit T coef.           [m3 kg-1 ppt-1 degC-1]
real, parameter :: BET110 =     SPV210*r1_S0  ! dSpV_dS fit zs * T coef.      [m3 kg-1 ppt-1 degC-1]
real, parameter :: BET210 = 1.5*SPV310*r1_S0  ! dSpV_dS fit zs**2 * T coef.   [m3 kg-1 ppt-1 degC-1]
real, parameter :: BET310 = 2.0*SPV410*r1_S0  ! dSpV_dS fit zs**3 * T coef.   [m3 kg-1 ppt-1 degC-1]
real, parameter :: BET410 = 2.5*SPV510*r1_S0  ! dSpV_dS fit zs**4 * T coef.   [m3 kg-1 ppt-1 degC-1]
real, parameter :: BET020 = 0.5*SPV120*r1_S0  ! dSpV_dS fit T**2 coef.        [m3 kg-1 ppt-1 degC-2]
real, parameter :: BET120 =     SPV220*r1_S0  ! dSpV_dS fit zs * T**2 coef.   [m3 kg-1 ppt-1 degC-2]
real, parameter :: BET220 = 1.5*SPV320*r1_S0  ! dSpV_dS fit zs**2 * T**2 coef. [m3 kg-1 ppt-1 degC-2]
real, parameter :: BET320 = 2.0*SPV420*r1_S0  ! dSpV_dS fit zs**3 * T**2 coef. [m3 kg-1 ppt-1 degC-2]
real, parameter :: BET030 = 0.5*SPV130*r1_S0  ! dSpV_dS fit T**3 coef.        [m3 kg-1 ppt-1 degC-3]
real, parameter :: BET130 =     SPV230*r1_S0  ! dSpV_dS fit zs * T**3 coef.   [m3 kg-1 ppt-1 degC-3]
real, parameter :: BET230 = 1.5*SPV330*r1_S0  ! dSpV_dS fit zs**2 * T**3 coef. [m3 kg-1 ppt-1 degC-3]
real, parameter :: BET040 = 0.5*SPV140*r1_S0  ! dSpV_dS fit T**4 coef.        [m3 kg-1 ppt-1 degC-4]
real, parameter :: BET140 =     SPV240*r1_S0  ! dSpV_dS fit zs * T**4 coef.   [m3 kg-1 ppt-1 degC-4]
real, parameter :: BET050 = 0.5*SPV150*r1_S0  ! dSpV_dS fit T**5 coef.        [m3 kg-1 ppt-1 degC-5]
real, parameter :: BET001 = 0.5*SPV101*r1_S0  ! dSpV_dS fit P coef.             [m3 kg-1 ppt-1 Pa-1]
real, parameter :: BET101 =     SPV201*r1_S0  ! dSpV_dS fit zs * P coef.        [m3 kg-1 ppt-1 Pa-1]
real, parameter :: BET201 = 1.5*SPV301*r1_S0  ! dSpV_dS fit zs**2 * P coef.     [m3 kg-1 ppt-1 Pa-1]
real, parameter :: BET301 = 2.0*SPV401*r1_S0  ! dSpV_dS fit zs**3 * P coef.     [m3 kg-1 ppt-1 Pa-1]
real, parameter :: BET011 = 0.5*SPV111*r1_S0  ! dSpV_dS fit T * P coef.  [m3 kg-1 ppt-1 degC-1 Pa-1]
real, parameter :: BET111 =     SPV211*r1_S0  ! dSpV_dS fit zs * T * P coef. [m3 kg-1 ppt-1 degC-1 Pa-1]
real, parameter :: BET211 = 1.5*SPV311*r1_S0  ! dSpV_dS fit zs**2 * T * P coef. [m3 kg-1 ppt-1 degC-1 Pa-1]
real, parameter :: BET021 = 0.5*SPV121*r1_S0  ! dSpV_dS fit T**2 * P coef. [m3 kg-1 ppt-1 degC-2 Pa-1]
real, parameter :: BET121 =     SPV221*r1_S0  ! dSpV_dS fit zs * T**2 * P coef. [m3 kg-1 ppt-1 degC-2 Pa-1]
real, parameter :: BET031 = 0.5*SPV131*r1_S0  ! dSpV_dS fit T**3 * P coef. [m3 kg-1 ppt-1 degC-3 Pa-1]
real, parameter :: BET002 = 0.5*SPV102*r1_S0  ! dSpV_dS fit P**2 coef.          [m3 kg-1 ppt-1 Pa-2]
real, parameter :: BET102 =     SPV202*r1_S0  ! dSpV_dS fit zs * P**2 coef.     [m3 kg-1 ppt-1 Pa-2]
real, parameter :: BET012 = 0.5*SPV112*r1_S0  ! dSpV_dS fit T * P**2 coef. [m3 kg-1 ppt-1 degC-1 Pa-2]
real, parameter :: BET003 = 0.5*SPV103*r1_S0  ! dSpV_dS fit P**3 coef.          [m3 kg-1 ppt-1 Pa-3]
!>@}

contains

!> Computes the Roquet et al. in situ specific volume of sea water for scalar inputs and outputs.
!!
!! Returns the in situ specific volume of sea water (specvol in [m3 kg-1]) from absolute salinity (S [g kg-1]),
!! conservative temperature (T [degC]) and pressure [Pa].  It uses the specific volume polynomial
!! fit from Roquet et al. (2015).
!! If spv_ref is present, specvol is an anomaly from spv_ref.
subroutine calculate_spec_vol_scalar_Roquet_SpV(T, S, pressure, specvol, spv_ref)
  real,           intent(in)  :: T        !< Conservative temperature [degC]
  real,           intent(in)  :: S        !< Absolute salinity [g kg-1]
  real,           intent(in)  :: pressure !< Pressure [Pa]
  real,           intent(out) :: specvol  !< In situ specific volume [m3 kg-1]
  real, optional, intent(in)  :: spv_ref  !< A reference specific volume [m3 kg-1]

  ! Local variables
  real, dimension(1) :: T0    ! A 1-d array with a copy of the conservative temperature [degC]
  real, dimension(1) :: S0    ! A 1-d array with a copy of the absolutes salinity [g kg-1]
  real, dimension(1) :: pres0 ! A 1-d array with a copy of the pressure [Pa]
  real, dimension(1) :: spv0  ! A 1-d array with a copy of the specific volume [m3 kg-1]

  T0(1) = T ; S0(1) = S ; pres0(1) = pressure

  call calculate_spec_vol_array_Roquet_SpV(T0, S0, pres0, spv0, 1, 1, spv_ref)
  specvol = spv0(1)

end subroutine calculate_spec_vol_scalar_Roquet_SpV

!> Computes the Roquet et al. in situ specific volume of sea water for 1-d array inputs and outputs.
!!
!! Returns the in situ specific volume of sea water (specvol in [m3 kg-1]) from absolute salinity (S [g kg-1]),
!! conservative temperature (T [degC]) and pressure [Pa].  It uses the specific volume polynomial
!! fit from Roquet et al. (2015).
!! If spv_ref is present, specvol is an anomaly from spv_ref.
subroutine calculate_spec_vol_array_Roquet_SpV(T, S, pressure, specvol, start, npts, spv_ref)
  real, dimension(:), intent(in)    :: T        !< Conservative temperature [degC]
  real, dimension(:), intent(in)    :: S        !< Absolute salinity [g kg-1]
  real, dimension(:), intent(in)    :: pressure !< pressure [Pa]
  real, dimension(:), intent(inout) :: specvol  !< in situ specific volume [m3 kg-1]
  integer,            intent(in)    :: start    !< The starting index for calculations
  integer,            intent(in)    :: npts     !< the number of values to calculate
  real,     optional, intent(in)    :: spv_ref  !< A reference specific volume [m3 kg-1]

  ! Local variables
  real :: zp     ! Pressure [Pa]
  real :: zt     ! Conservative temperature [degC]
  real :: zs     ! The square root of absolute salinity with an offset normalized
                 ! by an assumed salinity range [nondim]
  real :: SV_00p ! A pressure-dependent but temperature and salinity independent contribution to
                 ! specific volume at the reference temperature and salinity [m3 kg-1]
  real :: SV_TS  ! Specific volume without a pressure-dependent contribution [m3 kg-1]
  real :: SV_TS0 ! A contribution to specific volume from temperature and salinity anomalies at
                 ! the surface pressure [m3 kg-1]
  real :: SV_TS1 ! A temperature and salinity dependent specific volume contribution that is
                 ! proportional to pressure [m3 kg-1 Pa-1]
  real :: SV_TS2 ! A temperature and salinity dependent specific volume contribution that is
                 ! proportional to pressure**2 [m3 kg-1 Pa-2]
  real :: SV_TS3 ! A temperature and salinity dependent specific volume contribution that is
                 ! proportional to pressure**3 [m3 kg-1 Pa-3]
  real :: SV_0S0 ! Salinity dependent specific volume at the surface pressure and zero temperature [m3 kg-1]
  integer :: j

  ! The following algorithm was published by Roquet et al. (2015), intended for use in non-Boussinesq ocean models.
  do j=start,start+npts-1
    ! Conversions to the units used here.
    zt = T(j)
    zs = SQRT( ABS( S(j) + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
    zp = pressure(j)

    ! The next two lines should be used if it is necessary to convert potential temperature and
    ! practical salinity to conservative temperature and absolute salinity.
    ! zt = gsw_ct_from_pt(S(j),T(j)) ! Convert potential temp to conservative temp [degC]
    ! zs = SQRT( ABS( gsw_sr_from_sp(S(j)) + rdeltaS ) * r1_S0 ) ! Convert S from practical to absolute salinity.

    SV_TS3 = SPV003 + (zs*SPV103 + zt*SPV013)
    SV_TS2 = SPV002 + (zs*(SPV102 +  zs*SPV202) &
                     + zt*(SPV012 + (zs*SPV112 + zt*SPV022)) )
    SV_TS1 = SPV001 + (zs*(SPV101 +  zs*(SPV201 +  zs*(SPV301 +  zs*SPV401))) &
                     + zt*(SPV011 + (zs*(SPV111 +  zs*(SPV211 +  zs*SPV311)) &
                                   + zt*(SPV021 + (zs*(SPV121 +  zs*SPV221) &
                                                 + zt*(SPV031 + (zs*SPV131 + zt*SPV041)) )) )) )
    SV_TS0 = zt*(SPV010 &
               + (zs*(SPV110 +  zs*(SPV210 +  zs*(SPV310 +  zs*(SPV410 +  zs*SPV510)))) &
                + zt*(SPV020 + (zs*(SPV120 +  zs*(SPV220 +  zs*(SPV320 +  zs*SPV420))) &
                              + zt*(SPV030 + (zs*(SPV130 +  zs*(SPV230 +  zs*SPV330)) &
                                            + zt*(SPV040 + (zs*(SPV140 +  zs*SPV240) &
                                                          + zt*(SPV050 + (zs*SPV150 + zt*SPV060)) )) )) )) ) )

    SV_0S0 = SPV000 + zs*(SPV100 + zs*(SPV200 + zs*(SPV300 + zs*(SPV400 + zs*(SPV500 + zs*SPV600)))))

    SV_00p = zp*(V00 + zp*(V01 + zp*(V02 + zp*(V03 + zp*(V04 + zp*V05)))))

    if (present(spv_ref)) SV_0S0 = SV_0S0 - spv_ref

    SV_TS  = (SV_TS0 + SV_0S0) + zp*(SV_TS1 + zp*(SV_TS2 +  zp*SV_TS3))
    specvol(j) = SV_TS + SV_00p  ! In situ specific volume [m3 kg-1]
  enddo

end subroutine calculate_spec_vol_array_Roquet_SpV


!> Compute the in situ density of sea water at a point (rho in [kg m-3]) from absolute
!! salinity (S [g kg-1]), conservative temperature (T [degC]) and pressure [Pa], using the
!! specific volume polynomial fit from Roquet et al. (2015).
subroutine calculate_density_scalar_Roquet_SpV(T, S, pressure, rho, rho_ref)
  real,           intent(in)  :: T        !< Conservative temperature [degC]
  real,           intent(in)  :: S        !< Absolute salinity [g kg-1]
  real,           intent(in)  :: pressure !< Pressure [Pa]
  real,           intent(out) :: rho      !< In situ density [kg m-3]
  real, optional, intent(in)  :: rho_ref  !< A reference density [kg m-3]

  real, dimension(1) :: T0    ! A 1-d array with a copy of the conservative temperature [degC]
  real, dimension(1) :: S0    ! A 1-d array with a copy of the absolute salinity [g kg-1]
  real, dimension(1) :: pres0 ! A 1-d array with a copy of the pressure [Pa]
  real, dimension(1) :: spv   ! A 1-d array with the specific volume [m3 kg-1]

  T0(1) = T
  S0(1) = S
  pres0(1) = pressure

  if (present(rho_ref)) then
    call calculate_spec_vol_array_Roquet_SpV(T0, S0, pres0, spv, 1, 1, spv_ref=1.0/rho_ref)
    rho = -rho_ref**2*spv(1) / (rho_ref*spv(1) + 1.0)  ! In situ density [kg m-3]
  else
    call calculate_spec_vol_array_Roquet_SpV(T0, S0, pres0, spv, 1, 1)
    rho = 1.0 / spv(1)
  endif

end subroutine calculate_density_scalar_Roquet_SpV

!> Compute an array of in situ densities of sea water (rho in [kg m-3]) from absolute
!! salinity (S [g kg-1]), conservative temperature (T [degC]) and pressure [Pa],
!! using the specific volume polynomial fit from Roquet et al. (2015).
subroutine calculate_density_array_Roquet_SpV(T, S, pressure, rho, start, npts, rho_ref)
  real, dimension(:), intent(in)  :: T        !< Conservative temperature [degC]
  real, dimension(:), intent(in)  :: S        !< Absolute salinity [g kg-1]
  real, dimension(:), intent(in)  :: pressure !< Pressure [Pa]
  real, dimension(:), intent(out) :: rho      !< In situ density [kg m-3]
  integer,            intent(in)  :: start    !< The starting index for calculations
  integer,            intent(in)  :: npts     !< The number of values to calculate
  real,     optional, intent(in)  :: rho_ref  !< A reference density [kg m-3]

  ! Local variables
  real, dimension(size(T)) :: spv   ! The specific volume [m3 kg-1]
  integer :: j

  if (present(rho_ref)) then
    call calculate_spec_vol_array_Roquet_SpV(T, S, pressure, spv, start, npts, spv_ref=1.0/rho_ref)
    do j=start,start+npts-1
      rho(j) = -rho_ref**2*spv(j) / (rho_ref*spv(j) + 1.0)  ! In situ density [kg m-3]
    enddo
  else
    call calculate_spec_vol_array_Roquet_SpV(T, S, pressure, spv, start, npts)
    do j=start,start+npts-1
      rho(j) = 1.0 / spv(j)  ! In situ density [kg m-3]
    enddo
  endif

end subroutine calculate_density_array_Roquet_SpV

!> Return the partial derivatives of specific volume with temperature and salinity for 1-d array
!! inputs and outputs, using the specific volume polynomial fit from Roquet et al. (2015).
subroutine calculate_specvol_derivs_Roquet_SpV(T, S, pressure, dSV_dT, dSV_dS, start, npts)
  real,    intent(in),    dimension(:) :: T        !< Conservative temperature [degC]
  real,    intent(in),    dimension(:) :: S        !< Absolute salinity [g kg-1]
  real,    intent(in),    dimension(:) :: pressure !< Pressure [Pa]
  real,    intent(inout), dimension(:) :: dSV_dT   !< The partial derivative of specific volume with
                                                   !! conservative temperature [m3 kg-1 degC-1]
  real,    intent(inout), dimension(:) :: dSV_dS   !< The partial derivative of specific volume with
                                                   !! absolute salinity [m3 kg-1 ppt-1]
  integer, intent(in)                  :: start    !< The starting index for calculations
  integer, intent(in)                  :: npts     !< The number of values to calculate

  real :: zp      ! Pressure [Pa]
  real :: zt      ! Conservative temperature [degC]
  real :: zs      ! The square root of absolute salinity with an offset normalized
                  ! by an assumed salinity range [nondim]
  real :: dSVdzt0 ! A contribution to the partial derivative of specific volume with temperature
                  ! from temperature anomalies at the surface pressure [m3 kg-1 degC-1]
  real :: dSVdzt1 ! A contribution to the partial derivative of specific volume with temperature
                  ! that is proportional to pressure [m3 kg-1 degC-1 Pa-1]
  real :: dSVdzt2 ! A contribution to the partial derivative of specific volume with temperature
                  ! that is proportional to pressure**2 [m3 kg-1 degC-1 Pa-2]
  real :: dSVdzt3 ! A contribution to the partial derivative of specific volume with temperature
                  ! that is proportional to pressure**3 [m3 kg-1 degC-1 Pa-3]
  real :: dSVdzs0 ! A contribution to the partial derivative of specific volume with
                  ! salinity [m3 kg-1 ppt-1] from temperature anomalies at the surface pressure
  real :: dSVdzs1 ! A contribution to the partial derivative of specific volume with
                  ! salinity [m3 kg-1 ppt-1 Pa-1] proportional to pressure
  real :: dSVdzs2 ! A contribution to the partial derivative of specific volume with
                  ! salinity [m3 kg-1 ppt-1 Pa-2] proportional to pressure**2
  real :: dSVdzs3 ! A contribution to the partial derivative of specific volume with
                  ! salinity [m3 kg-1 ppt-1 Pa-3] proportional to pressure**3
  integer :: j

  do j=start,start+npts-1
    ! Conversions to the units used here.
    zt = T(j)
    zs = SQRT( ABS( S(j) + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
    zp = pressure(j)

    ! The next two lines should be used if it is necessary to convert potential temperature and
    ! practical salinity to conservative temperature and absolute salinity.
    ! zt = gsw_ct_from_pt(S(j),T(j)) ! Convert potential temp to conservative temp [degC]
    ! zs = SQRT( ABS( gsw_sr_from_sp(S(j)) + rdeltaS ) * r1_S0 ) ! Convert S from practical to absolute salinity.

    ! Find the partial derivative of specific volume with temperature
    dSVdzt3 = ALP003
    dSVdzt2 = ALP002 + (zs*ALP102 + zt*ALP012)
    dSVdzt1 = ALP001 + (zs*(ALP101 + zs*(ALP201 + zs*ALP301)) &
                      + zt*(ALP011 + (zs*(ALP111 + zs*ALP211) &
                                    + zt*(ALP021 + (zs*ALP121 + zt*ALP031)) )) )
    dSVdzt0 = ALP000 + (zs*(ALP100 +  zs*(ALP200 +  zs*(ALP300 + zs*(ALP400 + zs*ALP500)))) &
                      + zt*(ALP010 + (zs*(ALP110 +  zs*(ALP210 + zs*(ALP310 + zs*ALP410))) &
                                    + zt*(ALP020 + (zs*(ALP120 + zs*(ALP220 + zs*ALP320)) &
                                                  + zt*(ALP030 + (zt*(ALP040 + (zs*ALP140 + zt*ALP050)) &
                                                                + zs*(ALP130 + zs*ALP230) )) )) )) )

    dSV_dT(j) = dSVdzt0 + zp*(dSVdzt1 + zp*(dSVdzt2 + zp*dSVdzt3))

    ! Find the partial derivative of specific volume with salinity
    dSVdzs3 = BET003
    dSVdzs2 = BET002 + (zs*BET102 + zt*BET012)
    dSVdzs1 = BET001 + (zs*(BET101 + zs*(BET201 + zs*BET301)) &
                      + zt*(BET011 + (zs*(BET111 + zs*BET211) &
                                    + zt*(BET021 + (zs*BET121 + zt*BET031)) )) )
    dSVdzs0 = BET000 + (zs*(BET100 + zs*(BET200 + zs*(BET300 + zs*(BET400 + zs*BET500)))) &
                      + zt*(BET010 + (zs*(BET110 + zs*(BET210 + zs*(BET310 + zs*BET410))) &
                                    + zt*(BET020 + (zs*(BET120 + zs*(BET220 + zs*BET320)) &
                                                  + zt*(BET030 + (zt*(BET040 + (zs*BET140 + zt*BET050)) &
                                                                + zs*(BET130 + zs*BET230) )) )) )) )

    ! The division by zs here is because zs = sqrt(S + S0), so dSV_dS = dzs_dS * dSV_dzs = (0.5 / zs) * dSV_dzs
    dSV_dS(j) = (dSVdzs0 + zp*(dSVdzs1 + zp*(dSVdzs2 + zp * dSVdzs3))) / zs
  enddo

end subroutine calculate_specvol_derivs_Roquet_SpV


!> Compute an array of derivatives of densities of sea water with temperature (drho_dT in [kg m-3 degC-1])
!! and salinity (drho_dS in [kg m-3 ppt-1]) from absolute salinity (S [g kg-1]), conservative temperature
!! (T [degC]) and pressure [Pa], using the specific volume polynomial fit from Roquet et al. (2015).
subroutine calculate_density_derivs_array_Roquet_SpV(T, S, pressure, drho_dT, drho_dS, start, npts)
  real,    intent(in),  dimension(:) :: T        !< Conservative temperature [degC]
  real,    intent(in),  dimension(:) :: S        !< Absolute salinity [g kg-1]
  real,    intent(in),  dimension(:) :: pressure !< pressure [Pa]
  real,    intent(out), dimension(:) :: drho_dT  !< The partial derivative of density with
                                                 !! conservative temperature [kg m-3 degC-1]
  real,    intent(out), dimension(:) :: drho_dS  !< The partial derivative of density with
                                                 !! absolute salinity [kg m-3 ppt-1]
  integer, intent(in)                :: start    !< The starting index for calculations
  integer, intent(in)                :: npts     !< The number of values to calculate

  ! Local variables
  real, dimension(size(T)) :: specvol  ! The specific volume [m3 kg-1]
  real, dimension(size(T)) :: dSV_dT   ! The partial derivative of specific volume with
                                       ! conservative temperature [m3 kg-1 degC-1]
  real, dimension(size(T)) :: dSV_dS   ! The partial derivative of specific volume with
                                       ! absolute salinity [m3 kg-1 ppt-1]
  real :: rho  ! The in situ density [kg m-3]
  integer :: j

  call calculate_spec_vol_array_Roquet_SpV(T, S, pressure, specvol, start, npts)
  call calculate_specvol_derivs_Roquet_SpV(T, S, pressure, dSV_dT, dSV_dS, start, npts)

  do j=start,start+npts-1
    rho = 1.0 / specvol(j)
    drho_dT(j) = -dSv_dT(j) * rho**2
    drho_dS(j) = -dSv_dS(j) * rho**2
  enddo

end subroutine calculate_density_derivs_array_Roquet_SpV

!> Wrapper to calculate_density_derivs_array_Roquet_SpV for scalar inputs
subroutine calculate_density_derivs_scalar_Roquet_SpV(T, S, pressure, drho_dt, drho_ds)
  real,    intent(in)  :: T        !< Conservative temperature [degC]
  real,    intent(in)  :: S        !< Absolute salinity [g kg-1]
  real,    intent(in)  :: pressure !< Pressure [Pa]
  real,    intent(out) :: drho_dT  !< The partial derivative of density with
                                   !! conservative temperature [kg m-3 degC-1]
  real,    intent(out) :: drho_dS  !< The partial derivative of density with
                                   !! absolute salinity [kg m-3 ppt-1]
  ! Local variables
  real, dimension(1) :: T0    ! A 1-d array with a copy of the conservative temperature [degC]
  real, dimension(1) :: S0    ! A 1-d array with a copy of the absolute salinity [g kg-1]
  real, dimension(1) :: pres0 ! A 1-d array with a copy of the pressure [Pa]
  real, dimension(1) :: drdt0 ! A 1-d array with a copy of the derivative of density
                              ! with conservative temperature [kg m-3 degC-1]
  real, dimension(1) :: drds0 ! A 1-d array with a copy of the derivative of density
                              ! with absolute salinity [kg m-3 ppt-1]

  T0(1) = T
  S0(1) = S
  pres0(1) = pressure

  call calculate_density_derivs_array_Roquet_SpV(T0, S0, pres0, drdt0, drds0, 1, 1)
  drho_dt = drdt0(1)
  drho_ds = drds0(1)
end subroutine calculate_density_derivs_scalar_Roquet_SpV

!> Compute the in situ density of sea water (rho in [kg m-3]) and the compressibility
!! (drho/dp = C_sound^-2, stored as drho_dp [s2 m-2]) from absolute salinity (sal [g kg-1]),
!! conservative temperature (T [degC]), and pressure [Pa], using the specific volume
!! polynomial fit from Roquet et al. (2015).
subroutine calculate_compress_Roquet_SpV(T, S, pressure, rho, drho_dp, start, npts)
  real,    intent(in),  dimension(:) :: T        !< Conservative temperature [degC]
  real,    intent(in),  dimension(:) :: S        !< Absolute salinity [g kg-1]
  real,    intent(in),  dimension(:) :: pressure !< pressure [Pa]
  real,    intent(out), dimension(:) :: rho      !< In situ density [kg m-3]
  real,    intent(out), dimension(:) :: drho_dp  !< The partial derivative of density with pressure
                                                 !! (also the inverse of the square of sound speed)
                                                 !! [s2 m-2]
  integer, intent(in)                :: start    !< The starting index for calculations
  integer, intent(in)                :: npts     !< The number of values to calculate

  ! Local variables
  real :: zp     ! Pressure [Pa]
  real :: zt     ! Conservative temperature [degC]
  real :: zs     ! The square root of absolute salinity with an offset normalized
                 ! by an assumed salinity range [nondim]
  real :: dSV_00p_dp ! Derivative of the pressure-dependent reference specific volume profile with
                 ! pressure [m3 kg-1 Pa-1]
  real :: dSV_TS_dp  ! Derivative of the specific volume anomaly from the reference profile with
                 ! pressure [m3 kg-1 Pa-1]
  real :: SV_00p ! A pressure-dependent but temperature and salinity independent contribution to
                 ! specific volume at the reference temperature and salinity [m3 kg-1]
  real :: SV_TS  ! specific volume without a pressure-dependent contribution [m3 kg-1]
  real :: SV_TS0 ! A contribution to specific volume from temperature and salinity anomalies at
                 ! the surface pressure [m3 kg-1]
  real :: SV_TS1 ! A temperature and salinity dependent specific volume contribution that is
                 ! proportional to pressure [m3 kg-1 Pa-1]
  real :: SV_TS2 ! A temperature and salinity dependent specific volume contribution that is
                 ! proportional to pressure**2 [m3 kg-1 Pa-2]
  real :: SV_TS3 ! A temperature and salinity dependent specific volume contribution that is
                 ! proportional to pressure**3 [m3 kg-1 Pa-3]
  real :: SV_0S0 ! Salinity dependent specific volume at the surface pressure and zero temperature [m3 kg-1]
  real :: dSpecVol_dp ! The partial derivative of specific volume with pressure [m3 kg-1 Pa-1]
  integer :: j

  ! The following algorithm was published by Roquet et al. (2015), intended for use
  ! with NEMO, but it is not necessarily the algorithm used in NEMO ocean model.
  do j=start,start+npts-1
    ! Conversions to the units used here.
    zt = T(j)
    zs = SQRT( ABS( S(j) + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
    zp = pressure(j)

    ! The next two lines should be used if it is necessary to convert potential temperature and
    ! practical salinity to conservative temperature and absolute salinity.
    ! zt = gsw_ct_from_pt(S(j),T(j)) ! Convert potential temp to conservative temp [degC]
    ! zs = SQRT( ABS( gsw_sr_from_sp(S(j)) + rdeltaS ) * r1_S0 ) ! Convert S from practical to absolute salinity.

    SV_TS3 = SPV003 + (zs*SPV103 + zt*SPV013)
    SV_TS2 = SPV002 + (zs*(SPV102 +  zs*SPV202) &
                     + zt*(SPV012 + (zs*SPV112 + zt*SPV022)) )
    SV_TS1 = SPV001 + (zs*(SPV101 +  zs*(SPV201 +  zs*(SPV301 +  zs*SPV401))) &
                     + zt*(SPV011 + (zs*(SPV111 +  zs*(SPV211 +  zs*SPV311)) &
                                   + zt*(SPV021 + (zs*(SPV121 +  zs*SPV221) &
                                                 + zt*(SPV031 + (zs*SPV131 + zt*SPV041)) )) )) )

    SV_TS0 = zt*(SPV010 &
               + (zs*(SPV110 +  zs*(SPV210 +  zs*(SPV310 +  zs*(SPV410 +  zs*SPV510)))) &
                + zt*(SPV020 + (zs*(SPV120 +  zs*(SPV220 +  zs*(SPV320 +  zs*SPV420))) &
                              + zt*(SPV030 + (zs*(SPV130 +  zs*(SPV230 +  zs*SPV330)) &
                                            + zt*(SPV040 + (zs*(SPV140 +  zs*SPV240) &
                                                          + zt*(SPV050 + (zs*SPV150 + zt*SPV060)) )) )) )) ) )

    SV_0S0 = SPV000 + zs*(SPV100 + zs*(SPV200 + zs*(SPV300 + zs*(SPV400 + zs*(SPV500 + zs*SPV600)))))

    SV_00p = zp*(V00 + zp*(V01 + zp*(V02 + zp*(V03 + zp*(V04 + zp*V05)))))

    SV_TS  = (SV_TS0 + SV_0S0) + zp*(SV_TS1 + zp*(SV_TS2 +  zp*SV_TS3))
    ! specvol = SV_TS + SV_00p ! In situ specific volume [m3 kg-1]
    rho(j) = 1.0 / (SV_TS + SV_00p) ! In situ density [kg m-3]

    dSV_00p_dp = V00 + zp*(2.*V01 + zp*(3.*V02 + zp*(4.*V03 + zp*(5.*V04 + zp*(6.*V05)))))
    dSV_TS_dp  = SV_TS1 + zp*(2.*SV_TS2 + zp*(3.*SV_TS3))
    dSpecVol_dp = dSV_TS_dp + dSV_00p_dp  !  [m3 kg-1 Pa-1]
    drho_dp(j) = -dSpecVol_dp * rho(j)**2 ! Compressibility [s2 m-2]

  enddo
end subroutine calculate_compress_Roquet_SpV


!> Second derivatives of specific volume with respect to temperature, salinity, and pressure for a
!! 1-d array inputs and outputs using the specific volume polynomial fit from Roquet et al. (2015).
subroutine calc_spec_vol_second_derivs_array_Roquet_SpV(T, S, P, dSV_ds_ds, dSV_ds_dt, dSV_dt_dt, &
                                                        dSV_ds_dp, dSV_dt_dp, start, npts)
  real, dimension(:), intent(in   ) :: T !< Conservative temperature [degC]
  real, dimension(:), intent(in   ) :: S !< Absolute salinity [g kg-1]
  real, dimension(:), intent(in   ) :: P !< Pressure [Pa]
  real, dimension(:), intent(inout) :: dSV_ds_ds  !< Second derivative of specific volume with respect
                                                  !! to salinity [m3 kg-1 ppt-2]
  real, dimension(:), intent(inout) :: dSV_ds_dt  !< Second derivative of specific volume with respect
                                                  !! to salinity and temperature [m3 kg-1 ppt-1 degC-1]
  real, dimension(:), intent(inout) :: dSV_dt_dt  !< Second derivative of specific volume with respect
                                                  !! to temperature [m3 kg-1 degC-2]
  real, dimension(:), intent(inout) :: dSV_ds_dp  !< Second derivative of specific volume with respect to pressure
                                                  !! and salinity [m3 kg-1 ppt-1 Pa-1]
  real, dimension(:), intent(inout) :: dSV_dt_dp  !< Second derivative of specific volume with respect to pressure
                                                  !! and temperature [m3 kg-1 degC-1 Pa-1]
  integer,            intent(in   ) :: start      !< The starting index for calculations
  integer,            intent(in   ) :: npts       !< The number of values to calculate

  ! Local variables
  real :: zp      ! Pressure [Pa]
  real :: zt      ! Conservative temperature [degC]
  real :: zs      ! The square root of absolute salinity with an offset normalized
                  ! by an assumed salinity range [nondim]
  real :: I_s     ! The inverse of zs [nondim]
  real :: d2SV_p0 ! A contribution to one of the second derivatives that is independent of pressure [various]
  real :: d2SV_p1 ! A contribution to one of the second derivatives that is proportional to pressure [various]
  real :: d2SV_p2 ! A contribution to one of the second derivatives that is proportional to pressure**2 [various]
  real :: d2SV_p3 ! A contribution to one of the second derivatives that is proportional to pressure**3 [various]
  integer :: j

  do j = start,start+npts-1
    ! Conversions to the units used here.
    zt = T(j)
    zs = SQRT( ABS( S(j) + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
    zp = P(j)

    ! The next two lines should be used if it is necessary to convert potential temperature and
    ! practical salinity to conservative temperature and absolute salinity.
    ! zt = gsw_ct_from_pt(S(j),T(j)) ! Convert potential temp to conservative temp [degC]
    ! zs = SQRT( ABS( gsw_sr_from_sp(S(j)) + rdeltaS ) * r1_S0 )  ! Convert S from practical to absolute salinity.

    I_s = 1.0 / zs

    ! Find dSV_ds_ds
    d2SV_p3 = -SPV103*I_s**2
    d2SV_p2 = -(SPV102 + zt*SPV112)*I_s**2
    d2SV_p1 = (3.*SPV301 + (zt*(3.*SPV311) + zs*(8.*SPV401))) &
              - ( SPV101 + zt*(SPV111 + zt*(SPV121 + zt*SPV131)) )*I_s**2
    d2SV_p0 = (3.*SPV300 + (zs*(8.*SPV400 + zs*(15.*SPV500 + zs*(24.*SPV600))) &
                          + zt*(3.*SPV310 + (zs*(8.*SPV410 + zs*(15.*SPV510)) &
                                           + zt*(3.*SPV320 + (zs*(8.*SPV420) + zt*(3.*SPV330))) )) )) &
              - (SPV100 + zt*(SPV110 + zt*(SPV120 + zt*(SPV130 + zt*(SPV140 + zt*SPV150)))) )*I_s**2
    dSV_dS_dS(j) = (0.5*r1_S0)**2 * ((d2SV_p0 + zp*(d2SV_p1 + zp*(d2SV_p2 + zp*d2SV_p3))) * I_s)

    ! Find dSV_ds_dt
    d2SV_p2 = SPV112
    d2SV_p1 = SPV111 + (zs*(2.*SPV211 +  zs*(3.*SPV311)) &
                      + zt*(2.*SPV121 + (zs*(4.*SPV221) + zt*(3.*SPV131))) )
    d2SV_p0 = SPV110 + (zs*(2.*SPV210 +  zs*(3.*SPV310 +  zs*(4.*SPV410 +  zs*(5.*SPV510)))) &
                      + zt*(2.*SPV120 + (zs*(4.*SPV220 +  zs*(6.*SPV320 +  zs*(8.*SPV420))) &
                                       + zt*(3.*SPV130 + (zs*(6.*SPV230 +  zs*(9.*SPV330)) &
                                                        + zt*(4.*SPV140 + (zs*(8.*SPV240) &
                                                                         + zt*(5.*SPV150))) )) )) )
    dSV_ds_dt(j) = (0.5*r1_S0) * ((d2SV_p0 + zp*(d2SV_p1 + zp*d2SV_p2)) * I_s)

    ! Find dSV_dt_dt
    d2SV_p2 = 2.*SPV022
    d2SV_p1 = 2.*SPV021 + (zs*(2.*SPV121 +  zs*(2.*SPV221)) &
                         + zt*(6.*SPV031 + (zs*(6.*SPV131) + zt*(12.*SPV041))) )
    d2SV_p0 = 2.*SPV020 + (zs*(2.*SPV120 +  zs*( 2.*SPV220 +  zs*( 2.*SPV320 + zs * (2.*SPV420)))) &
                         + zt*(6.*SPV030 + (zs*( 6.*SPV130 +  zs*( 6.*SPV230 + zs * (6.*SPV330))) &
                                          + zt*(12.*SPV040 + (zs*(12.*SPV140 + zs *(12.*SPV240)) &
                                                            + zt*(20.*SPV050 + (zs*(20.*SPV150) &
                                                                              + zt*(30.*SPV060) )) )) )) )
    dSV_dt_dt(j) = d2SV_p0 + zp*(d2SV_p1 + zp*d2SV_p2)

    ! Find dSV_ds_dp
    d2SV_p2 = 3.*SPV103
    d2SV_p1 = 2.*SPV102 + (zs*(4.*SPV202) + zt*(2.*SPV112))
    d2SV_p0 = SPV101 + (zs*(2.*SPV201 + zs*(3.*SPV301 +  zs*(4.*SPV401))) &
                      + zt*(SPV111 +   (zs*(2.*SPV211 +  zs*(3.*SPV311)) &
                                      + zt*(   SPV121 + (zs*(2.*SPV221) + zt*SPV131)) )) )
    dSV_ds_dp(j) =  ((d2SV_p0 + zp*(d2SV_p1 + zp*d2SV_p2)) * I_s) * (0.5*r1_S0)

    ! Find dSV_dt_dp
    d2SV_p2 = 3.*SPV013
    d2SV_p1 = 2.*SPV012 + (zs*(2.*SPV112) + zt*(4.*SPV022))
    d2SV_p0 = SPV011 + (zs*(SPV111     + zs*(   SPV211 +  zs*    SPV311)) &
                      + zt*(2.*SPV021 + (zs*(2.*SPV121 +  zs*(2.*SPV221)) &
                                       + zt*(3.*SPV031 + (zs*(3.*SPV131) + zt*(4.*SPV041))) )) )
    dSV_dt_dp(j) =  d2SV_p0 + zp*(d2SV_p1 + zp*d2SV_p2)
  enddo

end subroutine calc_spec_vol_second_derivs_array_Roquet_SpV


!> Second derivatives of density with respect to temperature, salinity, and pressure for a
!! 1-d array inputs and outputs using the specific volume polynomial fit from Roquet et al. (2015).
subroutine calculate_density_second_derivs_array_Roquet_SpV(T, S, P, drho_ds_ds, drho_ds_dt, drho_dt_dt, &
                                                            drho_ds_dp, drho_dt_dp, start, npts)
  real, dimension(:), intent(in   ) :: T !< Conservative temperature [degC]
  real, dimension(:), intent(in   ) :: S !< Absolute salinity [g kg-1]
  real, dimension(:), intent(in   ) :: P !< Pressure [Pa]
  real, dimension(:), intent(inout) :: drho_ds_ds !< Second derivative of density with respect
                                                  !! to salinity [kg m-3 ppt-2]
  real, dimension(:), intent(inout) :: drho_ds_dt !< Second derivative of density with respect
                                                  !! to salinity and temperature [kg m-3 ppt-1 degC-1]
  real, dimension(:), intent(inout) :: drho_dt_dt !< Second derivative of density with respect
                                                  !! to temperature [kg m-3 degC-2]
  real, dimension(:), intent(inout) :: drho_ds_dp !< Second derivative of density with respect to pressure
                                                  !! and salinity [kg m-3 ppt-1 Pa-1] = [s2 m-2 ppt-1]
  real, dimension(:), intent(inout) :: drho_dt_dp !< Second derivative of density with respect to pressure
                                                  !! and temperature [kg m-3 degC-1 Pa-1] = [s2 m-2 degC-1]
  integer,            intent(in   ) :: start      !< The starting index for calculations
  integer,            intent(in   ) :: npts       !< The number of values to calculate

  ! Local variables
  real, dimension(size(T)) :: rho       ! The in situ density [kg m-3]
  real, dimension(size(T)) :: drho_dp   ! The partial derivative of density with pressure
                                        ! (also the inverse of the square of sound speed)  [s2 m-2]
  real, dimension(size(T)) :: dSV_dT    ! The partial derivative of specific volume with
                                        ! conservative temperature [m3 kg-1 degC-1]
  real, dimension(size(T)) :: dSV_dS    ! The partial derivative of specific volume with
                                        ! absolute salinity [m3 kg-1 ppt-1]
  real, dimension(size(T)) :: dSV_ds_ds ! Second derivative of specific volume with respect
                                        ! to salinity [m3 kg-1 ppt-2]
  real, dimension(size(T)) :: dSV_ds_dt ! Second derivative of specific volume with respect
                                        ! to salinity and temperature [m3 kg-1 ppt-1 degC-1]
  real, dimension(size(T)) :: dSV_dt_dt ! Second derivative of specific volume with respect
                                        ! to temperature [m3 kg-1 degC-2]
  real, dimension(size(T)) :: dSV_ds_dp ! Second derivative of specific volume with respect to pressure
                                        ! and salinity [m3 kg-1 ppt-1 Pa-1]
  real, dimension(size(T)) :: dSV_dt_dp ! Second derivative of specific volume with respect to pressure
                                        ! and temperature [m3 kg-1 degC-1 Pa-1]
  integer :: j

  call calc_spec_vol_second_derivs_array_Roquet_SpV(T, S, P, dSV_ds_ds, dSV_ds_dt, dSV_dt_dt, &
                                                    dSV_ds_dp, dSV_dt_dp, start, npts)
  call calculate_specvol_derivs_Roquet_SpV(T, S, P, dSV_dT, dSV_dS, start, npts)
  call calculate_compress_Roquet_SpV(T, S, P, rho, drho_dp, start, npts)

  do j = start,start+npts-1
    ! Find drho_ds_ds
    drho_dS_dS(j) = rho(j)**2 * (2.0*rho(j)*dSV_dS(j)**2 - dSV_dS_dS(j))

    ! Find drho_ds_dt
    drho_ds_dt(j) = rho(j)**2 * (2.0*rho(j)*(dSV_dT(j)*dSV_dS(j)) - dSV_dS_dT(j))

    ! Find drho_dt_dt
    drho_dT_dT(j) = rho(j)**2 * (2.0*rho(j)*dSV_dT(j)**2 - dSV_dT_dT(j))

    ! Find drho_ds_dp
    drho_ds_dp(j) =  -rho(j) * (2.0*dSV_dS(j) * drho_dp(j) + rho(j) * dSV_dS_dp(j))

    ! Find drho_dt_dp
    drho_dt_dp(j) =  -rho(j) * (2.0*dSV_dT(j) * drho_dp(j) + rho(j) * dSV_dT_dp(j))
  enddo

end subroutine calculate_density_second_derivs_array_Roquet_SpV

!> Second derivatives of density with respect to temperature, salinity, and pressure for scalar inputs.
!!
!! The scalar version of calculate_density_second_derivs promotes scalar inputs to 1-element array
!! and then demotes the output back to a scalar
subroutine calculate_density_second_derivs_scalar_Roquet_SpV(T, S, P, drho_ds_ds, drho_ds_dt, drho_dt_dt, &
                                                             drho_ds_dp, drho_dt_dp)
  real, intent(in   ) :: T          !< Conservative temperature [degC]
  real, intent(in   ) :: S          !< Absolute salinity [g kg-1]
  real, intent(in   ) :: P          !< pressure [Pa]
  real, intent(  out) :: drho_ds_ds !< Second derivative of density with respect
                                    !! to salinity [kg m-3 ppt-2]
  real, intent(  out) :: drho_ds_dt !< Second derivative of density with respect
                                    !! to salinity and temperature [kg m-3 ppt-1 degC-1]
  real, intent(  out) :: drho_dt_dt !< Second derivative of density with respect
                                    !! to temperature [kg m-3 degC-2]
  real, intent(  out) :: drho_ds_dp !< Second derivative of density with respect to pressure
                                    !! and salinity [kg m-3 ppt-1 Pa-1] = [s2 m-2 ppt-1]
  real, intent(  out) :: drho_dt_dp !< Second derivative of density with respect to pressure
                                    !! and temperature [kg m-3 degC-1 Pa-1] = [s2 m-2 degC-1]
  ! Local variables
  real, dimension(1) :: T0     ! A 1-d array with a copy of the temperature [degC]
  real, dimension(1) :: S0     ! A 1-d array with a copy of the salinity [g kg-1]
  real, dimension(1) :: p0     ! A 1-d array with a copy of the pressure [Pa]
  real, dimension(1) :: drdsds ! The second derivative of density with salinity [kg m-3 ppt-2]
  real, dimension(1) :: drdsdt ! The second derivative of density with salinity and
                               ! temperature [kg m-3 ppt-1 degC-1]
  real, dimension(1) :: drdtdt ! The second derivative of density with temperature [kg m-3 degC-2]
  real, dimension(1) :: drdsdp ! The second derivative of density with salinity and
                               ! pressure [kg m-3 ppt-1 Pa-1] = [s2 m-2 ppt-1]
  real, dimension(1) :: drdtdp ! The second derivative of density with temperature and
                               ! pressure [kg m-3 degC-1 Pa-1] = [s2 m-2 degC-1]

  T0(1) = T
  S0(1) = S
  P0(1) = P
  call calculate_density_second_derivs_array_Roquet_SpV(T0, S0, P0, drdsds, drdsdt, drdtdt, drdsdp, drdtdp, 1, 1)
  drho_ds_ds = drdsds(1)
  drho_ds_dt = drdsdt(1)
  drho_dt_dt = drdtdt(1)
  drho_ds_dp = drdsdp(1)
  drho_dt_dp = drdtdp(1)

end subroutine calculate_density_second_derivs_scalar_Roquet_SpV

!> Return the range of temperatures, salinities and pressures for which the Roquet et al. (2015)
!! expression for specific volume has been fitted to observations.  Care should be taken when
!! applying this equation of state outside of its fit range.
subroutine EoS_fit_range_Roquet_SpV(T_min, T_max, S_min, S_max, p_min, p_max)
  real, optional, intent(out) :: T_min !< The minimum conservative temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: T_max !< The maximum conservative temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: S_min !< The minimum absolute salinity over which this EoS is fitted [g kg-1]
  real, optional, intent(out) :: S_max !< The maximum absolute salinity over which this EoS is fitted [g kg-1]
  real, optional, intent(out) :: p_min !< The minimum pressure over which this EoS is fitted [Pa]
  real, optional, intent(out) :: p_max !< The maximum pressure over which this EoS is fitted [Pa]

  if (present(T_min)) T_min = -6.0
  if (present(T_max)) T_max = 40.0
  if (present(S_min)) S_min =  0.0
  if (present(S_max)) S_max = 42.0
  if (present(p_min)) p_min = 0.0
  if (present(p_max)) p_max = 1.0e8

end subroutine EoS_fit_range_Roquet_SpV

!> \namespace mom_eos_Roquet_SpV
!!
!! \section section_EOS_Roquet_SpV NEMO equation of state
!!
!!  Fabien Roquet and colleagues developed this equation of state using a simple polynomial fit
!! to the TEOS-10 equation of state expressions for specific, for efficiency when used with a
!! non-Boussinesq ocean model.  This particular equation of state is a balance between an
!! accuracy that matches the TEOS-10 density to better than observational uncertainty with a
!! polynomial form that can be evaluated quickly despite having 55 terms.
!!
!! \subsection section_EOS_Roquet_Spv_references References
!!
!! Roquet, F., Madec, G., McDougall, T. J., and Barker, P. M., 2015:
!!  Accurate polynomial expressions for the density and specific volume
!!  of seawater using the TEOS-10 standard. Ocean Modelling, 90:29-43.

end module MOM_EOS_Roquet_Spv
