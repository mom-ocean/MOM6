!> The equation of state using the expressions of Roquet et al. (2015) that are used in NEMO
module MOM_EOS_Roquet_rho

! This file is part of MOM6. See LICENSE.md for the license.

!use gsw_mod_toolbox, only : gsw_sr_from_sp, gsw_ct_from_pt

implicit none ; private

public calculate_compress_Roquet_rho, calculate_density_Roquet_rho
public calculate_density_derivs_Roquet_rho
public calculate_density_scalar_Roquet_rho, calculate_density_array_Roquet_rho
public calculate_density_second_derivs_Roquet_rho, EoS_fit_range_Roquet_rho

!> Compute the in situ density of sea water [kg m-3], or its anomaly with respect to
!! a reference density, from absolute salinity [g kg-1], conservative temperature [degC],
!! and pressure [Pa], using the expressions for density from Roquet et al. (2015)
interface calculate_density_Roquet_rho
  module procedure calculate_density_scalar_Roquet_rho, calculate_density_array_Roquet_rho
end interface calculate_density_Roquet_rho

!> For a given thermodynamic state, return the derivatives of density with conservative temperature
!! and absolute salinity, using the expressions for density from Roquet et al. (2015)
interface calculate_density_derivs_Roquet_rho
  module procedure calculate_density_derivs_scalar_Roquet_rho, calculate_density_derivs_array_Roquet_rho
end interface calculate_density_derivs_Roquet_rho

!> Compute the second derivatives of density with various combinations of temperature,
!! salinity, and pressure using the expressions for density from Roquet et al. (2015)
interface calculate_density_second_derivs_Roquet_rho
  module procedure calculate_density_second_derivs_scalar_Roquet_rho
  module procedure calculate_density_second_derivs_array_Roquet_rho
end interface calculate_density_second_derivs_Roquet_rho

real, parameter :: Pa2kb  = 1.e-8 !< Conversion factor between Pa and kbar [kbar Pa-1]
!>@{ Parameters in the Roquet_rho (Roquet density) equation of state
real, parameter :: rdeltaS = 32.          ! An offset to salinity before taking its square root [g kg-1]
real, parameter :: r1_S0 = 0.875/35.16504 ! The inverse of a plausible range of oceanic salinities [kg g-1]
real, parameter :: I_Ts = 0.025           ! The inverse of a plausible range of oceanic temperatures [degC-1]

! The following are the coefficients of the fit to the reference density profile (rho00p) as a function of
! pressure (P), with a contribution R0c * P**(c+1).  The nomenclature follows Roquet.
real, parameter :: R00 = 4.6494977072e+01*Pa2kb     ! rho00p P coef.    [kg m-3 Pa-1]
real, parameter :: R01 = -5.2099962525*Pa2kb**2     ! rho00p P**2 coef. [kg m-3 Pa-2]
real, parameter :: R02 = 2.2601900708e-01*Pa2kb**3  ! rho00p P**3 coef. [kg m-3 Pa-3]
real, parameter :: R03 = 6.4326772569e-02*Pa2kb**4  ! rho00p P**4 coef. [kg m-3 Pa-4]
real, parameter :: R04 = 1.5616995503e-02*Pa2kb**5  ! rho00p P**5 coef. [kg m-3 Pa-5]
real, parameter :: R05 = -1.7243708991e-03*Pa2kb**6 ! rho00p P**6 coef. [kg m-3 Pa-6]

! The following are coefficients of contributions to density as a function of the square root
! of normalized salinity with an offset (zs), temperature (T) and pressure (P), with a contribution
! EOSabc * zs**a * T**b * P**c.  The numbers here are copied directly from Roquet et al. (2015), but
! the expressions here do not use the same nondimensionalization for pressure or temperature as they do.
real, parameter :: EOS000 = 8.0189615746e+02                  ! A constant density contribution [kg m-3]
real, parameter :: EOS100 = 8.6672408165e+02                  ! EoS zs coef.                [kg m-3]
real, parameter :: EOS200 = -1.7864682637e+03                 ! EoS zs**2 coef.             [kg m-3]
real, parameter :: EOS300 = 2.0375295546e+03                  ! EoS zs**3 coef.             [kg m-3]
real, parameter :: EOS400 = -1.2849161071e+03                 ! EoS zs**4 coef.             [kg m-3]
real, parameter :: EOS500 = 4.3227585684e+02                  ! EoS zs**5 coef.             [kg m-3]
real, parameter :: EOS600 = -6.0579916612e+01                 ! EoS zs**6 coef.             [kg m-3]
real, parameter :: EOS010 = 2.6010145068e+01*I_Ts             ! EoS T coef.          [kg m-3 degC-1]
real, parameter :: EOS110 = -6.5281885265e+01*I_Ts            ! EoS zs * T coef.     [kg m-3 degC-1]
real, parameter :: EOS210 = 8.1770425108e+01*I_Ts             ! EoS zs**2 * T coef.  [kg m-3 degC-1]
real, parameter :: EOS310 = -5.6888046321e+01*I_Ts            ! EoS zs**3 * T coef.  [kg m-3 degC-1]
real, parameter :: EOS410 = 1.7681814114e+01*I_Ts             ! EoS zs**2 * T coef.  [kg m-3 degC-1]
real, parameter :: EOS510 = -1.9193502195*I_Ts                ! EoS zs**5 * T coef.  [kg m-3 degC-1]
real, parameter :: EOS020 = -3.7074170417e+01*I_Ts**2         ! EoS T**2 coef.       [kg m-3 degC-2]
real, parameter :: EOS120 = 6.1548258127e+01*I_Ts**2          ! EoS zs * T**2 coef.  [kg m-3 degC-2]
real, parameter :: EOS220 = -6.0362551501e+01*I_Ts**2         ! EoS zs**2 * T**2 coef. [kg m-3 degC-2]
real, parameter :: EOS320 = 2.9130021253e+01*I_Ts**2          ! EoS zs**3 * T**2 coef. [kg m-3 degC-2]
real, parameter :: EOS420 = -5.4723692739*I_Ts**2             ! EoS zs**4 * T**2 coef. [kg m-3 degC-2]
real, parameter :: EOS030 = 2.1661789529e+01*I_Ts**3          ! EoS T**3 coef.       [kg m-3 degC-3]
real, parameter :: EOS130 = -3.3449108469e+01*I_Ts**3         ! EoS zs * T**3 coef.  [kg m-3 degC-3]
real, parameter :: EOS230 = 1.9717078466e+01*I_Ts**3          ! EoS zs**2 * T**3 coef. [kg m-3 degC-3]
real, parameter :: EOS330 = -3.1742946532*I_Ts**3             ! EoS zs**3 * T**3 coef. [kg m-3 degC-3]
real, parameter :: EOS040 = -8.3627885467*I_Ts**4             ! EoS T**4 coef.       [kg m-3 degC-4]
real, parameter :: EOS140 = 1.1311538584e+01*I_Ts**4          ! EoS zs * T**4 coef.  [kg m-3 degC-4]
real, parameter :: EOS240 = -5.3563304045*I_Ts**4             ! EoS zs**2 * T**4 coef. [kg m-3 degC-4]
real, parameter :: EOS050 = 5.4048723791e-01*I_Ts**5          ! EoS T**5 coef.       [kg m-3 degC-5]
real, parameter :: EOS150 = 4.8169980163e-01*I_Ts**5          ! EoS zs * T**5 coef.  [kg m-3 degC-5]
real, parameter :: EOS060 = -1.9083568888e-01*I_Ts**6         ! EoS T**6             [kg m-3 degC-6]
real, parameter :: EOS001 = 1.9681925209e+01*Pa2kb            ! EoS P coef.            [kg m-3 Pa-1]
real, parameter :: EOS101 = -4.2549998214e+01*Pa2kb           ! EoS zs * P coef.       [kg m-3 Pa-1]
real, parameter :: EOS201 = 5.0774768218e+01*Pa2kb            ! EoS zs**2 * P coef.    [kg m-3 Pa-1]
real, parameter :: EOS301 = -3.0938076334e+01*Pa2kb           ! EoS zs**3 * P coef.    [kg m-3 Pa-1]
real, parameter :: EOS401 = 6.6051753097*Pa2kb                ! EoS zs**4 * P coef.    [kg m-3 Pa-1]
real, parameter :: EOS011 = -1.3336301113e+01*(I_Ts*Pa2kb)    ! EoS T * P coef. [kg m-3 degC-1 Pa-1]
real, parameter :: EOS111 = -4.4870114575*(I_Ts*Pa2kb)        ! EoS zs * T * P coef. [kg m-3 degC-1 Pa-1]
real, parameter :: EOS211 = 5.0042598061*(I_Ts*Pa2kb)         ! EoS zs**2 * T * P coef. [kg m-3 degC-1 Pa-1]
real, parameter :: EOS311 = -6.5399043664e-01*(I_Ts*Pa2kb)    ! EoS zs**3 * T * P coef. [kg m-3 degC-1 Pa-1]
real, parameter :: EOS021 = 6.7080479603*(I_Ts**2*Pa2kb)      ! EoS T**2 * P coef. [kg m-3 degC-2 Pa-1]
real, parameter :: EOS121 = 3.5063081279*(I_Ts**2*Pa2kb)      ! EoS zs * T**2 * P coef. [kg m-3 degC-2 Pa-1]
real, parameter :: EOS221 = -1.8795372996*(I_Ts**2*Pa2kb)     ! EoS zs**2 * T**2 * P coef. [kg m-3 degC-2 Pa-1]
real, parameter :: EOS031 = -2.4649669534*(I_Ts**3*Pa2kb)     ! EoS T**3 * P coef. [kg m-3 degC-3 Pa-1]
real, parameter :: EOS131 = -5.5077101279e-01*(I_Ts**3*Pa2kb) ! EoS zs * T**3 * P coef. [kg m-3 degC-3 Pa-1]
real, parameter :: EOS041 = 5.5927935970e-01*(I_Ts**4*Pa2kb)  ! EoS T**4 * P coef. [kg m-3 degC-4 Pa-1]
real, parameter :: EOS002 = 2.0660924175*Pa2kb**2             ! EoS P**2 coef.         [kg m-3 Pa-2]
real, parameter :: EOS102 = -4.9527603989*Pa2kb**2            ! EoS zs * P**2 coef.    [kg m-3 Pa-2]
real, parameter :: EOS202 = 2.5019633244*Pa2kb**2             ! EoS zs**2 * P**2 coef. [kg m-3 Pa-2]
real, parameter :: EOS012 = 2.0564311499*(I_Ts*Pa2kb**2)      ! EoS T * P**2 coef. [kg m-3 degC-1 Pa-2]
real, parameter :: EOS112 = -2.1311365518e-01*(I_Ts*Pa2kb**2) ! EoS zs * T * P**2 coef. [kg m-3 degC-1 Pa-2]
real, parameter :: EOS022 = -1.2419983026*(I_Ts**2*Pa2kb**2)  ! EoS T**2 * P**2 coef. [kg m-3 degC-2 Pa-2]
real, parameter :: EOS003 = -2.3342758797e-02*Pa2kb**3        ! EoS P**3 coef.         [kg m-3 Pa-3]
real, parameter :: EOS103 = -1.8507636718e-02*Pa2kb**3        ! EoS zs * P**3 coef.    [kg m-3 Pa-3]
real, parameter :: EOS013 = 3.7969820455e-01*(I_Ts*Pa2kb**3)  ! EoS T * P**3 coef. [kg m-3 degC-1 Pa-3]

real, parameter :: ALP000 =    EOS010   ! Constant in the drho_dT fit                [kg m-3 degC-1]
real, parameter :: ALP100 =    EOS110   ! drho_dT fit zs coef.                       [kg m-3 degC-1]
real, parameter :: ALP200 =    EOS210   ! drho_dT fit zs**2 coef.                    [kg m-3 degC-1]
real, parameter :: ALP300 =    EOS310   ! drho_dT fit zs**3 coef.                    [kg m-3 degC-1]
real, parameter :: ALP400 =    EOS410   ! drho_dT fit zs**4 coef.                    [kg m-3 degC-1]
real, parameter :: ALP500 =    EOS510   ! drho_dT fit zs**5 coef.                    [kg m-3 degC-1]
real, parameter :: ALP010 = 2.*EOS020   ! drho_dT fit T coef.                        [kg m-3 degC-2]
real, parameter :: ALP110 = 2.*EOS120   ! drho_dT fit zs * T coef.                   [kg m-3 degC-2]
real, parameter :: ALP210 = 2.*EOS220   ! drho_dT fit zs**2 * T coef.                [kg m-3 degC-2]
real, parameter :: ALP310 = 2.*EOS320   ! drho_dT fit zs**3 * T coef.                [kg m-3 degC-2]
real, parameter :: ALP410 = 2.*EOS420   ! drho_dT fit zs**4 * T coef.                [kg m-3 degC-2]
real, parameter :: ALP020 = 3.*EOS030   ! drho_dT fit T**2 coef.                     [kg m-3 degC-3]
real, parameter :: ALP120 = 3.*EOS130   ! drho_dT fit zs * T**2 coef.                [kg m-3 degC-3]
real, parameter :: ALP220 = 3.*EOS230   ! drho_dT fit zs**2 * T**2 coef.             [kg m-3 degC-3]
real, parameter :: ALP320 = 3.*EOS330   ! drho_dT fit zs**3 * T**2 coef.             [kg m-3 degC-3]
real, parameter :: ALP030 = 4.*EOS040   ! drho_dT fit T**3 coef.                     [kg m-3 degC-4]
real, parameter :: ALP130 = 4.*EOS140   ! drho_dT fit zs * T**3 coef.                [kg m-3 degC-4]
real, parameter :: ALP230 = 4.*EOS240   ! drho_dT fit zs**2 * T**3 coef.             [kg m-3 degC-4]
real, parameter :: ALP040 = 5.*EOS050   ! drho_dT fit T**4 coef.                     [kg m-3 degC-5]
real, parameter :: ALP140 = 5.*EOS150   ! drho_dT fit zs* * T**4 coef.               [kg m-3 degC-5]
real, parameter :: ALP050 = 6.*EOS060   ! drho_dT fit T**5 coef.                     [kg m-3 degC-6]
real, parameter :: ALP001 =    EOS011   ! drho_dT fit P coef.                   [kg m-3 degC-1 Pa-1]
real, parameter :: ALP101 =    EOS111   ! drho_dT fit zs * P coef.              [kg m-3 degC-1 Pa-1]
real, parameter :: ALP201 =    EOS211   ! drho_dT fit zs**2 * P coef.           [kg m-3 degC-1 Pa-1]
real, parameter :: ALP301 =    EOS311   ! drho_dT fit zs**3 * P coef.           [kg m-3 degC-1 Pa-1]
real, parameter :: ALP011 = 2.*EOS021   ! drho_dT fit T * P coef.               [kg m-3 degC-2 Pa-1]
real, parameter :: ALP111 = 2.*EOS121   ! drho_dT fit zs * T * P coef.          [kg m-3 degC-2 Pa-1]
real, parameter :: ALP211 = 2.*EOS221   ! drho_dT fit zs**2 * T * P coef.       [kg m-3 degC-2 Pa-1]
real, parameter :: ALP021 = 3.*EOS031   ! drho_dT fit T**2 * P coef.            [kg m-3 degC-3 Pa-1]
real, parameter :: ALP121 = 3.*EOS131   ! drho_dT fit zs * T**2 * P coef.       [kg m-3 degC-3 Pa-1]
real, parameter :: ALP031 = 4.*EOS041   ! drho_dT fit T**3 * P coef.            [kg m-3 degC-4 Pa-1]
real, parameter :: ALP002 =    EOS012   ! drho_dT fit P**2 coef.                [kg m-3 degC-1 Pa-2]
real, parameter :: ALP102 =    EOS112   ! drho_dT fit zs * P**2 coef.           [kg m-3 degC-1 Pa-2]
real, parameter :: ALP012 = 2.*EOS022   ! drho_dT fit T * P**2 coef.            [kg m-3 degC-2 Pa-2]
real, parameter :: ALP003 =    EOS013   ! drho_dT fit P**3 coef.                [kg m-3 degC-1 Pa-3]

real, parameter :: BET000 = 0.5*EOS100*r1_S0  ! Constant in the drho_dS fit           [kg m-3 ppt-1]
real, parameter :: BET100 =     EOS200*r1_S0  ! drho_dS fit zs coef.                  [kg m-3 ppt-1]
real, parameter :: BET200 = 1.5*EOS300*r1_S0  ! drho_dS fit zs**2 coef.               [kg m-3 ppt-1]
real, parameter :: BET300 = 2.0*EOS400*r1_S0  ! drho_dS fit zs**3 coef.               [kg m-3 ppt-1]
real, parameter :: BET400 = 2.5*EOS500*r1_S0  ! drho_dS fit zs**4 coef.               [kg m-3 ppt-1]
real, parameter :: BET500 = 3.0*EOS600*r1_S0  ! drho_dS fit zs**5 coef.               [kg m-3 ppt-1]
real, parameter :: BET010 = 0.5*EOS110*r1_S0  ! drho_dS fit T coef.            [kg m-3 ppt-1 degC-1]
real, parameter :: BET110 =     EOS210*r1_S0  ! drho_dS fit zs * T coef.       [kg m-3 ppt-1 degC-1]
real, parameter :: BET210 = 1.5*EOS310*r1_S0  ! drho_dS fit zs**2 * T coef.    [kg m-3 ppt-1 degC-1]
real, parameter :: BET310 = 2.0*EOS410*r1_S0  ! drho_dS fit zs**3 * T coef.    [kg m-3 ppt-1 degC-1]
real, parameter :: BET410 = 2.5*EOS510*r1_S0  ! drho_dS fit zs**4 * T coef.    [kg m-3 ppt-1 degC-1]
real, parameter :: BET020 = 0.5*EOS120*r1_S0  ! drho_dS fit T**2 coef.         [kg m-3 ppt-1 degC-2]
real, parameter :: BET120 =     EOS220*r1_S0  ! drho_dS fit zs * T**2 coef.    [kg m-3 ppt-1 degC-2]
real, parameter :: BET220 = 1.5*EOS320*r1_S0  ! drho_dS fit zs**2 * T**2 coef. [kg m-3 ppt-1 degC-2]
real, parameter :: BET320 = 2.0*EOS420*r1_S0  ! drho_dS fit zs**3 * T**2 coef. [kg m-3 ppt-1 degC-2]
real, parameter :: BET030 = 0.5*EOS130*r1_S0  ! drho_dS fit T**3 coef.         [kg m-3 ppt-1 degC-3]
real, parameter :: BET130 =     EOS230*r1_S0  ! drho_dS fit zs * T**3 coef.    [kg m-3 ppt-1 degC-3]
real, parameter :: BET230 = 1.5*EOS330*r1_S0  ! drho_dS fit zs**2 * T**3 coef. [kg m-3 ppt-1 degC-3]
real, parameter :: BET040 = 0.5*EOS140*r1_S0  ! drho_dS fit T**4 coef.         [kg m-3 ppt-1 degC-4]
real, parameter :: BET140 =     EOS240*r1_S0  ! drho_dS fit zs * T**4 coef.    [kg m-3 ppt-1 degC-4]
real, parameter :: BET050 = 0.5*EOS150*r1_S0  ! drho_dS fit T**5 coef.         [kg m-3 ppt-1 degC-5]
real, parameter :: BET001 = 0.5*EOS101*r1_S0  ! drho_dS fit P coef.              [kg m-3 ppt-1 Pa-1]
real, parameter :: BET101 =     EOS201*r1_S0  ! drho_dS fit zs * P coef.         [kg m-3 ppt-1 Pa-1]
real, parameter :: BET201 = 1.5*EOS301*r1_S0  ! drho_dS fit zs**2 * P coef.      [kg m-3 ppt-1 Pa-1]
real, parameter :: BET301 = 2.0*EOS401*r1_S0  ! drho_dS fit zs**3 * P coef.      [kg m-3 ppt-1 Pa-1]
real, parameter :: BET011 = 0.5*EOS111*r1_S0  ! drho_dS fit T * P coef.   [kg m-3 ppt-1 degC-1 Pa-1]
real, parameter :: BET111 =     EOS211*r1_S0  ! drho_dS fit zs * T * P coef. [kg m-3 ppt-1 degC-1 Pa-1]
real, parameter :: BET211 = 1.5*EOS311*r1_S0  ! drho_dS fit zs**2 * T * P coef. [kg m-3 ppt-1 degC-1 Pa-1]
real, parameter :: BET021 = 0.5*EOS121*r1_S0  ! drho_dS fit T**2 * P coef. [kg m-3 ppt-1 degC-2 Pa-1]
real, parameter :: BET121 =     EOS221*r1_S0  ! drho_dS fit zs * T**2 * P coef. [kg m-3 ppt-1 degC-2 Pa-1]
real, parameter :: BET031 = 0.5*EOS131*r1_S0  ! drho_dS fit T**3 * P coef. [kg m-3 ppt-1 degC-3 Pa-1]
real, parameter :: BET002 = 0.5*EOS102*r1_S0  ! drho_dS fit P**2 coef.           [kg m-3 ppt-1 Pa-2]
real, parameter :: BET102 =     EOS202*r1_S0  ! drho_dS fit zs * P**2 coef.      [kg m-3 ppt-1 Pa-2]
real, parameter :: BET012 = 0.5*EOS112*r1_S0  ! drho_dS fit T * P**2 coef. [kg m-3 ppt-1 degC-1 Pa-2]
real, parameter :: BET003 = 0.5*EOS103*r1_S0  ! drho_dS fit P**3 coef.           [kg m-3 ppt-1 Pa-3]
!>@}

contains

!> This subroutine computes the in situ density of sea water (rho in [kg m-3])
!! from absolute salinity (S [g kg-1]), conservative temperature (T [degC])
!! and pressure [Pa], using the density polynomial fit EOS from Roquet et al. (2015).
subroutine calculate_density_scalar_Roquet_rho(T, S, pres, rho, rho_ref)
  real,           intent(in)  :: T        !< Conservative temperature [degC]
  real,           intent(in)  :: S        !< Absolute salinity [g kg-1]
  real,           intent(in)  :: pres     !< Pressure [Pa]
  real,           intent(out) :: rho      !< In situ density [kg m-3]
  real, optional, intent(in)  :: rho_ref  !< A reference density [kg m-3]

  real, dimension(1) :: T0    ! A 1-d array with a copy of the conservative temperature [degC]
  real, dimension(1) :: S0    ! A 1-d array with a copy of the absolute salinity [g kg-1]
  real, dimension(1) :: pres0 ! A 1-d array with a copy of the pressure [Pa]
  real, dimension(1) :: rho0  ! A 1-d array with a copy of the density [kg m-3]

  T0(1) = T
  S0(1) = S
  pres0(1) = pres

  call calculate_density_array_Roquet_rho(T0, S0, pres0, rho0, 1, 1, rho_ref)
  rho = rho0(1)

end subroutine calculate_density_scalar_Roquet_rho

!> This subroutine computes an array of in situ densities of sea water (rho in [kg m-3])
!! from absolute salinity (S [g kg-1]), conservative temperature (T [degC]), and pressure
!! [Pa], using the density polynomial fit EOS from Roquet et al. (2015).
subroutine calculate_density_array_Roquet_rho(T, S, pres, rho, start, npts, rho_ref)
  real, dimension(:), intent(in)  :: T        !< Conservative temperature [degC]
  real, dimension(:), intent(in)  :: S        !< Absolute salinity [g kg-1]
  real, dimension(:), intent(in)  :: pres     !< Pressure [Pa]
  real, dimension(:), intent(out) :: rho      !< In situ density [kg m-3]
  integer,            intent(in)  :: start    !< The starting index for calculations
  integer,            intent(in)  :: npts     !< The number of values to calculate
  real,     optional, intent(in)  :: rho_ref  !< A reference density [kg m-3]

  ! Local variables
  real :: zp     ! Pressure [Pa]
  real :: zt     ! Conservative temperature [degC]
  real :: zs     ! The square root of absolute salinity with an offset normalized
                 ! by an assumed salinity range [nondim]
  real :: rho00p ! A pressure-dependent but temperature and salinity independent contribution to
                 ! density at the reference temperature and salinity [kg m-3]
  real :: rhoTS  ! Density without a pressure-dependent contribution [kg m-3]
  real :: rhoTS0 ! A contribution to density from temperature and salinity anomalies at the
                 ! surface pressure [kg m-3]
  real :: rhoTS1 ! A density contribution proportional to pressure [kg m-3 Pa-1]
  real :: rhoTS2 ! A density contribution proportional to pressure**2 [kg m-3 Pa-2]
  real :: rhoTS3 ! A density contribution proportional to pressure**3 [kg m-3 Pa-3]
  real :: rho0S0 ! Salinity dependent density at the surface pressure and zero temperature [kg m-3]
  integer :: j

  ! The following algorithm was published by Roquet et al. (2015), intended for use with NEMO.
  do j=start,start+npts-1
    ! Conversions to the units used here.
    zt = T(j)
    zs = SQRT( ABS( S(j) + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
    zp = pres(j)

    ! The next two lines should be used if it is necessary to convert potential temperature and
    ! practical salinity to conservative temperature and absolute salinity.
    ! zt = gsw_ct_from_pt(S(j),T(j)) ! Convert potential temp to conservative temp [degC]
    ! zs = SQRT( ABS( gsw_sr_from_sp(S(j)) + rdeltaS ) * r1_S0 ) ! Convert S from practical to absolute salinity.

    rhoTS3 = EOS003 + (zs*EOS103 + zt*EOS013)
    rhoTS2 = EOS002 + (zs*(EOS102 +  zs*EOS202) &
                     + zt*(EOS012 + (zs*EOS112 + zt*EOS022)) )
    rhoTS1 = EOS001 + (zs*(EOS101 +  zs*(EOS201 +  zs*(EOS301 +  zs*EOS401))) &
                     + zt*(EOS011 + (zs*(EOS111 +  zs*(EOS211 +  zs*EOS311)) &
                                   + zt*(EOS021 + (zs*(EOS121 +  zs*EOS221) &
                                                 + zt*(EOS031 + (zs*EOS131 + zt*EOS041)) )) )) )
    rhoTS0 = zt*(EOS010 &
               + (zs*(EOS110 +  zs*(EOS210 +  zs*(EOS310 +  zs*(EOS410 +  zs*EOS510)))) &
                + zt*(EOS020 + (zs*(EOS120 +  zs*(EOS220 +  zs*(EOS320 +  zs*EOS420))) &
                              + zt*(EOS030 + (zs*(EOS130 +  zs*(EOS230 +  zs*EOS330)) &
                                            + zt*(EOS040 + (zs*(EOS140 +  zs*EOS240) &
                                                          + zt*(EOS050 + (zs*EOS150 + zt*EOS060)) )) )) )) ) )

    rho0S0 = EOS000 + zs*(EOS100 + zs*(EOS200 + zs*(EOS300 + zs*(EOS400 + zs*(EOS500 + zs*EOS600)))))

    rho00p = zp*(R00 + zp*(R01 + zp*(R02 + zp*(R03 + zp*(R04 + zp*R05)))))

    if (present(rho_ref)) rho0S0 = rho0S0 - rho_ref

    rhoTS  = (rhoTS0 + rho0S0) + zp*(rhoTS1 + zp*(rhoTS2 +  zp*rhoTS3))
    rho(j) = rhoTS + rho00p  ! In situ density [kg m-3]

  enddo
end subroutine calculate_density_array_Roquet_rho

!> For a given thermodynamic state, calculate the derivatives of density with conservative
!! temperature and absolute salinity, using the density polynomial fit EOS from Roquet et al. (2015).
subroutine calculate_density_derivs_array_Roquet_rho(T, S, pres, drho_dT, drho_dS, start, npts)
  real,    intent(in),  dimension(:) :: T        !< Conservative temperature [degC]
  real,    intent(in),  dimension(:) :: S        !< Absolute salinity [g kg-1]
  real,    intent(in),  dimension(:) :: pres     !< Pressure [Pa]
  real,    intent(out), dimension(:) :: drho_dT  !< The partial derivative of density with
                                                 !! conservative temperature [kg m-3 degC-1]
  real,    intent(out), dimension(:) :: drho_dS  !< The partial derivative of density with
                                                 !! absolute salinity [kg m-3 ppt-1]
  integer, intent(in)                :: start    !< The starting index for calculations
  integer, intent(in)                :: npts     !< The number of values to calculate

  ! Local variables
  real :: zp      ! Pressure [Pa]
  real :: zt      ! Conservative temperature [degC]
  real :: zs      ! The square root of absolute salinity with an offset normalized
                  ! by an assumed salinity range [nondim]
  real :: dRdzt0  ! A contribution to the partial derivative of density with temperature [kg m-3 degC-1]
                  ! from temperature anomalies at the surface pressure
  real :: dRdzt1  ! A contribution to the partial derivative of density with temperature [kg m-3 degC-1 Pa-1]
                  ! proportional to pressure
  real :: dRdzt2  ! A contribution to the partial derivative of density with temperature [kg m-3 degC-1 Pa-2]
                  ! proportional to pressure**2
  real :: dRdzt3  ! A contribution to the partial derivative of density with temperature [kg m-3 degC-1 Pa-3]
                  ! proportional to pressure**3
  real :: dRdzs0  ! A contribution to the partial derivative of density with
                  ! salinity [kg m-3 ppt-1] from temperature anomalies at the surface pressure
  real :: dRdzs1  ! A contribution to the partial derivative of density with
                  ! salinity [kg m-3 ppt-1 Pa-1] proportional to pressure
  real :: dRdzs2  ! A contribution to the partial derivative of density with
                  ! salinity [kg m-3 ppt-1 Pa-2] proportional to pressure**2
  real :: dRdzs3  ! A contribution to the partial derivative of density with
                  ! salinity [kg m-3 ppt-1 Pa-3] proportional to pressure**3
  integer :: j

  do j=start,start+npts-1
    ! Conversions to the units used here.
    zt = T(j)
    zs = SQRT( ABS( S(j) + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
    zp = pres(j)

    ! The next two lines should be used if it is necessary to convert potential temperature and
    ! practical salinity to conservative temperature and absolute salinity.
    ! zt = gsw_ct_from_pt(S(j),T(j)) ! Convert potential temp to conservative temp [degC]
    ! zs = SQRT( ABS( gsw_sr_from_sp(S(j)) + rdeltaS ) * r1_S0 ) ! Convert S from practical to absolute salinity.

    ! Find the partial derivative of density with temperature
    dRdzt3 = ALP003
    dRdzt2 = ALP002 + (zs*ALP102 + zt*ALP012)
    dRdzt1 = ALP001 + (zs*(ALP101 + zs*(ALP201 + zs*ALP301)) &
                     + zt*(ALP011 + (zs*(ALP111 + zs*ALP211) &
                                   + zt*(ALP021 + (zs*ALP121 + zt*ALP031)) )) )
    dRdzt0 = ALP000 + (zs*(ALP100 +  zs*(ALP200 +  zs*(ALP300 + zs*(ALP400 + zs*ALP500)))) &
                     + zt*(ALP010 + (zs*(ALP110 +  zs*(ALP210 + zs*(ALP310 + zs*ALP410))) &
                                   + zt*(ALP020 + (zs*(ALP120 + zs*(ALP220 + zs*ALP320)) &
                                                 + zt*(ALP030 + (zt*(ALP040 + (zs*ALP140 + zt*ALP050)) &
                                                               + zs*(ALP130 + zs*ALP230) )) )) )) )

    drho_dT(j) = dRdzt0 + zp*(dRdzt1 + zp*(dRdzt2 + zp*dRdzt3))

    ! Find the partial derivative of density with salinity
    dRdzs3 = BET003
    dRdzs2 = BET002 + (zs*BET102 + zt*BET012)
    dRdzs1 = BET001 + (zs*(BET101 + zs*(BET201 + zs*BET301)) &
                     + zt*(BET011 + (zs*(BET111 + zs*BET211) &
                                   + zt*(BET021 + (zs*BET121 + zt*BET031)) )) )
    dRdzs0 = BET000 + (zs*(BET100 + zs*(BET200 + zs*(BET300 + zs*(BET400 + zs*BET500)))) &
                     + zt*(BET010 + (zs*(BET110 + zs*(BET210 + zs*(BET310 + zs*BET410))) &
                                   + zt*(BET020 + (zs*(BET120 + zs*(BET220 + zs*BET320)) &
                                                 + zt*(BET030 + (zt*(BET040 + (zs*BET140 + zt*BET050)) &
                                                               + zs*(BET130 + zs*BET230) )) )) )) )

    ! The division by zs here is because zs = sqrt(S + S0), so drho_dS = dzs_dS * drho_dzs = (0.5 / zs) * drho_dzs
    drho_dS(j) = (dRdzs0 + zp*(dRdzs1 + zp*(dRdzs2 + zp * dRdzs3))) / zs
  enddo

end subroutine calculate_density_derivs_array_Roquet_rho

!> Wrapper to calculate_density_derivs_array for scalar inputs
subroutine calculate_density_derivs_scalar_Roquet_rho(T, S, pres, drho_dt, drho_ds)
  real,    intent(in)  :: T        !< Conservative temperature [degC]
  real,    intent(in)  :: S        !< Absolute salinity [g kg-1]
  real,    intent(in)  :: pres     !< Pressure [Pa]
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
  pres0(1) = pres

  call calculate_density_derivs_array_Roquet_rho(T0, S0, pres0, drdt0, drds0, 1, 1)
  drho_dt = drdt0(1)
  drho_ds = drds0(1)
end subroutine calculate_density_derivs_scalar_Roquet_rho

!> Compute the in situ density of sea water (rho in [kg m-3]) and the compressibility
!! (drho/dp = C_sound^-2, stored as drho_dp [s2 m-2]) from absolute salinity (sal [g kg-1]),
!! conservative temperature (T [degC]), and pressure [Pa], using the density polynomial
!! fit EOS from Roquet et al. (2015).
subroutine calculate_compress_Roquet_rho(T, S, pres, rho, drho_dp, start, npts)
  real,    intent(in),  dimension(:) :: T        !< Conservative temperature [degC]
  real,    intent(in),  dimension(:) :: S        !< Absolute salinity [g kg-1]
  real,    intent(in),  dimension(:) :: pres     !< Pressure [Pa]
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
  real :: drho00p_dp ! Derivative of the pressure-dependent reference density profile with pressure [kg m-3 Pa-1]
  real :: drhoTS_dp  ! Derivative of the density anomaly from the reference profile with pressure [kg m-3 Pa-1]
  real :: rho00p ! The pressure-dependent (but temperature and salinity independent) reference
                 ! density profile [kg m-3]
  real :: rhoTS  ! Density anomaly from the reference profile [kg m-3]
  real :: rhoTS0 ! A contribution to density from temperature and salinity anomalies at the
                 ! surface pressure [kg m-3]
  real :: rhoTS1 ! A density contribution proportional to pressure [kg m-3 Pa-1]
  real :: rhoTS2 ! A density contribution proportional to pressure**2 [kg m-3 Pa-2]
  real :: rhoTS3 ! A density contribution proportional to pressure**3 [kg m-3 Pa-3]
  real :: rho0S0 ! Salinity dependent density at the surface pressure and zero temperature [kg m-3]
  integer :: j

  ! The following algorithm was published by Roquet et al. (2015), intended for use with NEMO.
  do j=start,start+npts-1
    ! Conversions to the units used here.
    zt = T(j)
    zs = SQRT( ABS( S(j) + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
    zp = pres(j)

    ! The next two lines should be used if it is necessary to convert potential temperature and
    ! practical salinity to conservative temperature and absolute salinity.
    ! zt = gsw_ct_from_pt(S(j),T(j)) ! Convert potential temp to conservative temp [degC]
    ! zs = SQRT( ABS( gsw_sr_from_sp(S(j)) + rdeltaS ) * r1_S0 ) ! Convert S from practical to absolute salinity.

    rhoTS3 = EOS003 + (zs*EOS103 + zt*EOS013)
    rhoTS2 = EOS002 + (zs*(EOS102 +  zs*EOS202) &
                     + zt*(EOS012 + (zs*EOS112 + zt*EOS022)) )
    rhoTS1 = EOS001 + (zs*(EOS101 +  zs*(EOS201 +  zs*(EOS301 +  zs*EOS401))) &
                     + zt*(EOS011 + (zs*(EOS111 +  zs*(EOS211 +  zs*EOS311)) &
                                   + zt*(EOS021 + (zs*(EOS121 +  zs*EOS221) &
                                                 + zt*(EOS031 + (zs*EOS131 + zt*EOS041)) )) )) )

    rhoTS0 = zt*(EOS010 &
               + (zs*(EOS110 +  zs*(EOS210 +  zs*(EOS310 +  zs*(EOS410 +  zs*EOS510)))) &
                + zt*(EOS020 + (zs*(EOS120 +  zs*(EOS220 +  zs*(EOS320 +  zs*EOS420))) &
                              + zt*(EOS030 + (zs*(EOS130 +  zs*(EOS230 +  zs*EOS330)) &
                                            + zt*(EOS040 + (zs*(EOS140 +  zs*EOS240) &
                                                          + zt*(EOS050 + (zs*EOS150 + zt*EOS060)) )) )) )) ) )

    rho0S0 = EOS000 + zs*(EOS100 + zs*(EOS200 + zs*(EOS300 + zs*(EOS400 + zs*(EOS500 + zs*EOS600)))))

    rho00p = zp*(R00 + zp*(R01 + zp*(R02 + zp*(R03 + zp*(R04 + zp*R05)))))

    rhoTS  = (rhoTS0 + rho0S0) + zp*(rhoTS1 + zp*(rhoTS2 +  zp*rhoTS3))
    rho(j) = rhoTS + rho00p ! In situ density [kg m-3]

    drho00p_dp = R00 + zp*(2.*R01 + zp*(3.*R02 + zp*(4.*R03 + zp*(5.*R04 + zp*(6.*R05)))))
    drhoTS_dp  = rhoTS1 + zp*(2.*rhoTS2 + zp*(3.*rhoTS3))
    drho_dp(j) = drhoTS_dp + drho00p_dp ! Compressibility [s2 m-2]

  enddo
end subroutine calculate_compress_Roquet_rho


!> Second derivatives of density with respect to temperature, salinity, and pressure for 1-d array
!! inputs and outputs.
subroutine calculate_density_second_derivs_array_Roquet_rho(T, S, P, drho_ds_ds, drho_ds_dt, drho_dt_dt, &
                                                      drho_ds_dp, drho_dt_dp, start, npts)
  real, dimension(:), intent(in   ) :: T          !< Conservative temperature [degC]
  real, dimension(:), intent(in   ) :: S          !< Absolute salinity [g kg-1] = [ppt]
  real, dimension(:), intent(in   ) :: P          !< Pressure [Pa]
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
  real :: zp     ! Pressure [Pa]
  real :: zt     ! Conservative temperature [degC]
  real :: zs     ! The square root of absolute salinity with an offset normalized
                 ! by an assumed salinity range [nondim]
  real :: I_s    ! The inverse of zs [nondim]
  real :: d2R_p0 ! A contribution to one of the second derivatives that is independent of pressure [various]
  real :: d2R_p1 ! A contribution to one of the second derivatives that is proportional to pressure [various]
  real :: d2R_p2 ! A contribution to one of the second derivatives that is proportional to pressure**2 [various]
  real :: d2R_p3 ! A contribution to one of the second derivatives that is proportional to pressure**3 [various]
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

    ! Find drho_ds_ds
    d2R_p3 = -EOS103*I_s**2
    d2R_p2 = -(EOS102 + zt*EOS112)*I_s**2
    d2R_p1 = (3.*EOS301 + (zt*(3.*EOS311) + zs*(8.*EOS401))) &
             - ( EOS101 + zt*(EOS111 + zt*(EOS121 + zt*EOS131)) )*I_s**2
    d2R_p0 = (3.*EOS300 + (zs*(8.*EOS400 + zs*(15.*EOS500 + zs*(24.*EOS600))) &
                         + zt*(3.*EOS310 + (zs*(8.*EOS410 + zs*(15.*EOS510)) &
                                          + zt*(3.*EOS320 + (zs*(8.*EOS420) + zt*(3.*EOS330))) )) )) &
             - (EOS100 + zt*(EOS110 + zt*(EOS120 + zt*(EOS130 + zt*(EOS140 + zt*EOS150)))) )*I_s**2
    drho_dS_dS(j) = (0.5*r1_S0)**2 * ((d2R_p0 + zp*(d2R_p1 + zp*(d2R_p2 + zp*d2R_p3))) * I_s)

    ! Find drho_ds_dt
    d2R_p2 = EOS112
    d2R_p1 = EOS111 + (zs*(2.*EOS211 +  zs*(3.*EOS311)) &
                     + zt*(2.*EOS121 + (zs*(4.*EOS221) + zt*(3.*EOS131))) )
    d2R_p0 = EOS110 + (zs*(2.*EOS210 +  zs*(3.*EOS310 +  zs*(4.*EOS410 +  zs*(5.*EOS510)))) &
                     + zt*(2.*EOS120 + (zs*(4.*EOS220 +  zs*(6.*EOS320 +  zs*(8.*EOS420))) &
                                      + zt*(3.*EOS130 + (zs*(6.*EOS230 +  zs*(9.*EOS330)) &
                                                       + zt*(4.*EOS140 + (zs*(8.*EOS240) &
                                                                        + zt*(5.*EOS150))) )) )) )
    drho_ds_dt(j) = (0.5*r1_S0) * ((d2R_p0 + zp*(d2R_p1 + zp*d2R_p2)) * I_s)

    ! Find drho_dt_dt
    d2R_p2 = 2.*EOS022
    d2R_p1 = 2.*EOS021 + (zs*(2.*EOS121 +  zs*(2.*EOS221)) &
                        + zt*(6.*EOS031 + (zs*(6.*EOS131) + zt*(12.*EOS041))) )
    d2R_p0 = 2.*EOS020 + (zs*(2.*EOS120 +  zs*( 2.*EOS220 +  zs*( 2.*EOS320 + zs * (2.*EOS420)))) &
                        + zt*(6.*EOS030 + (zs*( 6.*EOS130 +  zs*( 6.*EOS230 + zs * (6.*EOS330))) &
                                         + zt*(12.*EOS040 + (zs*(12.*EOS140 + zs *(12.*EOS240)) &
                                                           + zt*(20.*EOS050 + (zs*(20.*EOS150) &
                                                                             + zt*(30.*EOS060) )) )) )) )
    drho_dt_dt(j) = (d2R_p0 + zp*(d2R_p1 + zp*d2R_p2))

    ! Find drho_ds_dp
    d2R_p2 = 3.*EOS103
    d2R_p1 = 2.*EOS102 + (zs*(4.*EOS202) + zt*(2.*EOS112))
    d2R_p0 = EOS101 + (zs*(2.*EOS201 + zs*(3.*EOS301 +  zs*(4.*EOS401))) &
                     + zt*(EOS111 +   (zs*(2.*EOS211 +  zs*(3.*EOS311)) &
                                     + zt*(   EOS121 + (zs*(2.*EOS221) + zt*EOS131)) )) )
    drho_ds_dp(j) =  ((d2R_p0 + zp*(d2R_p1 + zp*d2R_p2)) * I_s) * (0.5*r1_S0)

    ! Find drho_dt_dp
    d2R_p2 = 3.*EOS013
    d2R_p1 = 2.*EOS012 + (zs*(2.*EOS112) + zt*(4.*EOS022))
    d2R_p0 = EOS011 + (zs*(EOS111     + zs*(   EOS211 +  zs*    EOS311)) &
                     + zt*(2.*EOS021 + (zs*(2.*EOS121 +  zs*(2.*EOS221)) &
                                      + zt*(3.*EOS031 + (zs*(3.*EOS131) + zt*(4.*EOS041))) )) )
    drho_dt_dp(j) =  (d2R_p0 + zp*(d2R_p1 + zp*d2R_p2))
  enddo

end subroutine calculate_density_second_derivs_array_Roquet_rho

!> Second derivatives of density with respect to temperature, salinity, and pressure for scalar inputs.
!!
!! The scalar version of calculate_density_second_derivs promotes scalar inputs to 1-element array
!! and then demotes the output back to a scalar
subroutine calculate_density_second_derivs_scalar_Roquet_rho(T, S, P, drho_ds_ds, drho_ds_dt, drho_dt_dt, &
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
  real, dimension(1) :: S0     ! A 1-d array with a copy of the salinity [g kg-1] = [ppt]
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
  call calculate_density_second_derivs_array_Roquet_rho(T0, S0, P0, drdsds, drdsdt, drdtdt, drdsdp, drdtdp, 1, 1)
  drho_ds_ds = drdsds(1)
  drho_ds_dt = drdsdt(1)
  drho_dt_dt = drdtdt(1)
  drho_ds_dp = drdsdp(1)
  drho_dt_dp = drdtdp(1)

end subroutine calculate_density_second_derivs_scalar_Roquet_rho

!> Return the range of temperatures, salinities and pressures for which the Roquet et al. (2015)
!! expression for in situ density has been fitted to observations.  Care should be taken when
!! applying this equation of state outside of its fit range.
subroutine EoS_fit_range_Roquet_rho(T_min, T_max, S_min, S_max, p_min, p_max)
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

end subroutine EoS_fit_range_Roquet_rho

!> \namespace mom_eos_Roquet_rho
!!
!! \section section_EOS_Roquet_rho Roquet_rho equation of state
!!
!!  Fabien Roquet and colleagues developed this equation of state using a simple polynomial fit
!! to the TEOS-10 equation of state, for efficiency when used in the NEMO ocean model.  Fabien
!! Roquet also graciously provided the MOM6 team with the original code implementing this
!! equation of state, although it has since been modified and extended to have capabilities
!! mirroring those available with other equations of state in MOM6.  This particular equation
!! of state is a balance between an accuracy that matches the TEOS-10 density to better than
!! observational uncertainty with a polynomial form that can be evaluated quickly despite having
!! 52 terms.
!!
!! \subsection section_EOS_Roquet_rho_references References
!!
!! Roquet, F., Madec, G., McDougall, T. J., and Barker, P. M., 2015:
!!  Accurate polynomial expressions for the density and specific volume
!!  of seawater using the TEOS-10 standard. Ocean Modelling, 90:29-43.

end module MOM_EOS_Roquet_rho
