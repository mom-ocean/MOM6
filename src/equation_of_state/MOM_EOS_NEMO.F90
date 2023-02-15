!> The equation of state using the expressions of Roquet et al. that are used in NEMO
module MOM_EOS_NEMO

! This file is part of MOM6. See LICENSE.md for the license.

!use gsw_mod_toolbox, only : gsw_sr_from_sp, gsw_ct_from_pt

implicit none ; private

public calculate_compress_nemo, calculate_density_nemo
public calculate_density_derivs_nemo
public calculate_density_scalar_nemo, calculate_density_array_nemo
public calculate_density_second_derivs_nemo

!> Compute the in situ density of sea water [kg m-3], or its anomaly with respect to
!! a reference density, from absolute salinity [g kg-1], conservative temperature [degC],
!! and pressure [Pa], using the expressions derived for use with NEMO
interface calculate_density_nemo
  module procedure calculate_density_scalar_nemo, calculate_density_array_nemo
end interface calculate_density_nemo

!> For a given thermodynamic state, return the derivatives of density with conservative temperature
!! and absolute salinity, the expressions derived for use with NEMO
interface calculate_density_derivs_nemo
  module procedure calculate_density_derivs_scalar_nemo, calculate_density_derivs_array_nemo
end interface calculate_density_derivs_nemo

!> Compute the second derivatives of density with various combinations
!! of temperature, salinity, and pressure
interface calculate_density_second_derivs_nemo
  module procedure calculate_density_second_derivs_scalar_nemo, calculate_density_second_derivs_array_nemo
end interface calculate_density_second_derivs_nemo

real, parameter :: Pa2db  = 1.e-4 !< Conversion factor between Pa and dbar [Pa dbar-1]
!>@{ Parameters in the NEMO equation of state
real, parameter :: rdeltaS = 32.    ! An offset to salinity before taking its square root [g kg-1]
real, parameter :: r1_S0  = 0.875/35.16504  ! The inverse of a plausible range of oceanic salinities [kg g-1]
real, parameter :: r1_T0  = 1./40.  ! The inverse of a plausible range of oceanic temperatures [degC-1]
real, parameter :: r1_P0  = 1.e-4   ! The inverse of a plausible range of oceanic pressures [dbar-1]
real, parameter :: R00 = 4.6494977072e+01  ! Contribution to zr0 proportional to zp [kg m-3]
real, parameter :: R01 = -5.2099962525     ! Contribution to zr0 proportional to zp**2 [kg m-3]
real, parameter :: R02 = 2.2601900708e-01  ! Contribution to zr0 proportional to zp**3 [kg m-3]
real, parameter :: R03 = 6.4326772569e-02  ! Contribution to zr0 proportional to zp**4 [kg m-3]
real, parameter :: R04 = 1.5616995503e-02  ! Contribution to zr0 proportional to zp**5 [kg m-3]
real, parameter :: R05 = -1.7243708991e-03 ! Contribution to zr0 proportional to zp**6 [kg m-3]

! The following terms are contributions to density as a function of the normalized square root of salinity
! with an offset (zs),  temperature (zt) and pressure, with a contribution EOSabc * zs**a * zt**b * zp**c
real, parameter :: EOS000 = 8.0189615746e+02  ! A constant density contribution [kg m-3]
real, parameter :: EOS100 = 8.6672408165e+02  ! Coefficient of the EOS proportional to zs [kg m-3]
real, parameter :: EOS200 = -1.7864682637e+03 ! Coefficient of the EOS proportional to zs**2 [kg m-3]
real, parameter :: EOS300 = 2.0375295546e+03  ! Coefficient of the EOS proportional to zs**3 [kg m-3]
real, parameter :: EOS400 = -1.2849161071e+03 ! Coefficient of the EOS proportional to zs**4 [kg m-3]
real, parameter :: EOS500 = 4.3227585684e+02  ! Coefficient of the EOS proportional to zs**5 [kg m-3]
real, parameter :: EOS600 = -6.0579916612e+01 ! Coefficient of the EOS proportional to zs**6 [kg m-3]
real, parameter :: EOS010 = 2.6010145068e+01  ! Coefficient of the EOS proportional to zt [kg m-3]
real, parameter :: EOS110 = -6.5281885265e+01 ! Coefficient of the EOS proportional to zs * zt [kg m-3]
real, parameter :: EOS210 = 8.1770425108e+01  ! Coefficient of the EOS proportional to zs**2 * zt [kg m-3]
real, parameter :: EOS310 = -5.6888046321e+01 ! Coefficient of the EOS proportional to zs**3 * zt [kg m-3]
real, parameter :: EOS410 = 1.7681814114e+01  ! Coefficient of the EOS proportional to zs**2 * zt [kg m-3]
real, parameter :: EOS510 = -1.9193502195     ! Coefficient of the EOS proportional to zs**5 * zt [kg m-3]
real, parameter :: EOS020 = -3.7074170417e+01 ! Coefficient of the EOS proportional to zt**2 [kg m-3]
real, parameter :: EOS120 = 6.1548258127e+01  ! Coefficient of the EOS proportional to zs * zt**2 [kg m-3]
real, parameter :: EOS220 = -6.0362551501e+01 ! Coefficient of the EOS proportional to zs**2 * zt**2 [kg m-3]
real, parameter :: EOS320 = 2.9130021253e+01  ! Coefficient of the EOS proportional to s**3 * zt**2 [kg m-3]
real, parameter :: EOS420 = -5.4723692739     ! Coefficient of the EOS proportional to zs**4 * zt**2 [kg m-3]
real, parameter :: EOS030 = 2.1661789529e+01  ! Coefficient of the EOS proportional to zt**3 [kg m-3]
real, parameter :: EOS130 = -3.3449108469e+01 ! Coefficient of the EOS proportional to zs * zt**3 [kg m-3]
real, parameter :: EOS230 = 1.9717078466e+01  ! Coefficient of the EOS proportional to zs**2 * zt**3 [kg m-3]
real, parameter :: EOS330 = -3.1742946532     ! Coefficient of the EOS proportional to zs**3 * zt**3 [kg m-3]
real, parameter :: EOS040 = -8.3627885467     ! Coefficient of the EOS proportional to zt**4 [kg m-3]
real, parameter :: EOS140 = 1.1311538584e+01  ! Coefficient of the EOS proportional to zs * zt**4 [kg m-3]
real, parameter :: EOS240 = -5.3563304045     ! Coefficient of the EOS proportional to zs**2 * zt**4 [kg m-3]
real, parameter :: EOS050 = 5.4048723791e-01  ! Coefficient of the EOS proportional to zt**5 [kg m-3]
real, parameter :: EOS150 = 4.8169980163e-01  ! Coefficient of the EOS proportional to zs * zt**5 [kg m-3]
real, parameter :: EOS060 = -1.9083568888e-01 ! Coefficient of the EOS proportional to zt**6 [kg m-3]
real, parameter :: EOS001 = 1.9681925209e+01  ! Coefficient of the EOS proportional to zp [kg m-3]
real, parameter :: EOS101 = -4.2549998214e+01 ! Coefficient of the EOS proportional to zs * zp [kg m-3]
real, parameter :: EOS201 = 5.0774768218e+01  ! Coefficient of the EOS proportional to zs**2 * zp [kg m-3]
real, parameter :: EOS301 = -3.0938076334e+01 ! Coefficient of the EOS proportional to zs**3 * zp [kg m-3]
real, parameter :: EOS401 = 6.6051753097      ! Coefficient of the EOS proportional to zs**4 * zp [kg m-3]
real, parameter :: EOS011 = -1.3336301113e+01 ! Coefficient of the EOS proportional to zt * zp [kg m-3]
real, parameter :: EOS111 = -4.4870114575     ! Coefficient of the EOS proportional to zs * zt * zp [kg m-3]
real, parameter :: EOS211 = 5.0042598061      ! Coefficient of the EOS proportional to zs**2 * zt * zp [kg m-3]
real, parameter :: EOS311 = -6.5399043664e-01 ! Coefficient of the EOS proportional to zs**3 * zt * zp [kg m-3]
real, parameter :: EOS021 = 6.7080479603      ! Coefficient of the EOS proportional to zt**2 * zp [kg m-3]
real, parameter :: EOS121 = 3.5063081279      ! Coefficient of the EOS proportional to zs * zt**2 * zp [kg m-3]
real, parameter :: EOS221 = -1.8795372996     ! Coefficient of the EOS proportional to zs**2 * zt**2 * zp [kg m-3]
real, parameter :: EOS031 = -2.4649669534     ! Coefficient of the EOS proportional to zt**3 * zp [kg m-3]
real, parameter :: EOS131 = -5.5077101279e-01 ! Coefficient of the EOS proportional to zs * zt**3 * zp [kg m-3]
real, parameter :: EOS041 = 5.5927935970e-01  ! Coefficient of the EOS proportional to zt**4 * zp [kg m-3]
real, parameter :: EOS002 = 2.0660924175      ! Coefficient of the EOS proportional to zp**2 [kg m-3]
real, parameter :: EOS102 = -4.9527603989     ! Coefficient of the EOS proportional to zs * zp**2 [kg m-3]
real, parameter :: EOS202 = 2.5019633244      ! Coefficient of the EOS proportional to zs**2 * zp**2 [kg m-3]
real, parameter :: EOS012 = 2.0564311499      ! Coefficient of the EOS proportional to zt * zp**2 [kg m-3]
real, parameter :: EOS112 = -2.1311365518e-01 ! Coefficient of the EOS proportional to zs * zt * zp**2 [kg m-3]
real, parameter :: EOS022 = -1.2419983026     ! Coefficient of the EOS proportional to zt**2 * zp**2 [kg m-3]
real, parameter :: EOS003 = -2.3342758797e-02 ! Coefficient of the EOS proportional to zp**3 [kg m-3]
real, parameter :: EOS103 = -1.8507636718e-02 ! Coefficient of the EOS proportional to zs * zp**3 [kg m-3]
real, parameter :: EOS013 = 3.7969820455e-01  ! Coefficient of the EOS proportional to zt * zp**3 [kg m-3]

real, parameter :: ALP000 =    EOS010*r1_T0   ! Constant in the drho_dT fit [kg m-3 degC-1]
real, parameter :: ALP100 =    EOS110*r1_T0   ! Coefficient of the drho_dT fit zs term [kg m-3 degC-1]
real, parameter :: ALP200 =    EOS210*r1_T0   ! Coefficient of the drho_dT fit zs**2 term [kg m-3 degC-1]
real, parameter :: ALP300 =    EOS310*r1_T0   ! Coefficient of the drho_dT fit zs**3 term [kg m-3 degC-1]
real, parameter :: ALP400 =    EOS410*r1_T0   ! Coefficient of the drho_dT fit zs**4 term [kg m-3 degC-1]
real, parameter :: ALP500 =    EOS510*r1_T0   ! Coefficient of the drho_dT fit zs**5 term [kg m-3 degC-1]
real, parameter :: ALP010 = 2.*EOS020*r1_T0   ! Coefficient of the drho_dT fit zt term [kg m-3 degC-1]
real, parameter :: ALP110 = 2.*EOS120*r1_T0   ! Coefficient of the drho_dT fit zs * zt term [kg m-3 degC-1]
real, parameter :: ALP210 = 2.*EOS220*r1_T0   ! Coefficient of the drho_dT fit zs**2 * zt term [kg m-3 degC-1]
real, parameter :: ALP310 = 2.*EOS320*r1_T0   ! Coefficient of the drho_dT fit zs**3 * zt term [kg m-3 degC-1]
real, parameter :: ALP410 = 2.*EOS420*r1_T0   ! Coefficient of the drho_dT fit zs**4 * zt term [kg m-3 degC-1]
real, parameter :: ALP020 = 3.*EOS030*r1_T0   ! Coefficient of the drho_dT fit zt**2 term [kg m-3 degC-1]
real, parameter :: ALP120 = 3.*EOS130*r1_T0   ! Coefficient of the drho_dT fit zs * zt**2 term [kg m-3 degC-1]
real, parameter :: ALP220 = 3.*EOS230*r1_T0   ! Coefficient of the drho_dT fit zs**2 * zt**2 term [kg m-3 degC-1]
real, parameter :: ALP320 = 3.*EOS330*r1_T0   ! Coefficient of the drho_dT fit zs**3 * zt**2 term [kg m-3 degC-1]
real, parameter :: ALP030 = 4.*EOS040*r1_T0   ! Coefficient of the drho_dT fit zt**3 term [kg m-3 degC-1]
real, parameter :: ALP130 = 4.*EOS140*r1_T0   ! Coefficient of the drho_dT fit zs * zt**3 term [kg m-3 degC-1]
real, parameter :: ALP230 = 4.*EOS240*r1_T0   ! Coefficient of the drho_dT fit zs**2 * zt**3 term [kg m-3 degC-1]
real, parameter :: ALP040 = 5.*EOS050*r1_T0   ! Coefficient of the drho_dT fit zt**4 term [kg m-3 degC-1]
real, parameter :: ALP140 = 5.*EOS150*r1_T0   ! Coefficient of the drho_dT fit zs* * zt**4 term [kg m-3 degC-1]
real, parameter :: ALP050 = 6.*EOS060*r1_T0   ! Coefficient of the drho_dT fit zt**5 term [kg m-3 degC-1]
real, parameter :: ALP001 =    EOS011*r1_T0   ! Coefficient of the drho_dT fit zp term [kg m-3 degC-1]
real, parameter :: ALP101 =    EOS111*r1_T0   ! Coefficient of the drho_dT fit zs * zp term [kg m-3 degC-1]
real, parameter :: ALP201 =    EOS211*r1_T0   ! Coefficient of the drho_dT fit zs**2 * zp term [kg m-3 degC-1]
real, parameter :: ALP301 =    EOS311*r1_T0   ! Coefficient of the drho_dT fit zs**3 * zp term [kg m-3 degC-1]
real, parameter :: ALP011 = 2.*EOS021*r1_T0   ! Coefficient of the drho_dT fit zt * zp term [kg m-3 degC-1]
real, parameter :: ALP111 = 2.*EOS121*r1_T0   ! Coefficient of the drho_dT fit zs * zt * zp term [kg m-3 degC-1]
real, parameter :: ALP211 = 2.*EOS221*r1_T0   ! Coefficient of the drho_dT fit zs**2 * zt * zp term [kg m-3 degC-1]
real, parameter :: ALP021 = 3.*EOS031*r1_T0   ! Coefficient of the drho_dT fit zt**2 * zp term [kg m-3 degC-1]
real, parameter :: ALP121 = 3.*EOS131*r1_T0   ! Coefficient of the drho_dT fit zs * zt**2 * zp term [kg m-3 degC-1]
real, parameter :: ALP031 = 4.*EOS041*r1_T0   ! Coefficient of the drho_dT fit zt**3 * zp term [kg m-3 degC-1]
real, parameter :: ALP002 =    EOS012*r1_T0   ! Coefficient of the drho_dT fit zp**2 term [kg m-3 degC-1]
real, parameter :: ALP102 =    EOS112*r1_T0   ! Coefficient of the drho_dT fit zs * zp**2 term [kg m-3 degC-1]
real, parameter :: ALP012 = 2.*EOS022*r1_T0   ! Coefficient of the drho_dT fit zt * zp**2 term [kg m-3 degC-1]
real, parameter :: ALP003 =    EOS013*r1_T0   ! Coefficient of the drho_dT fit zp**3 term [kg m-3 degC-1]

real, parameter :: BET000 = 0.5*EOS100*r1_S0  ! Constant in the drho_dS fit [kg m-3 ppt-1]
real, parameter :: BET100 =     EOS200*r1_S0  ! Coefficient of the drho_dS fit zs term [kg m-3 ppt-1]
real, parameter :: BET200 = 1.5*EOS300*r1_S0  ! Coefficient of the drho_dS fit zs**2 term [kg m-3 ppt-1]
real, parameter :: BET300 = 2.0*EOS400*r1_S0  ! Coefficient of the drho_dS fit zs**3 term [kg m-3 ppt-1]
real, parameter :: BET400 = 2.5*EOS500*r1_S0  ! Coefficient of the drho_dS fit zs**4 term [kg m-3 ppt-1]
real, parameter :: BET500 = 3.0*EOS600*r1_S0  ! Coefficient of the drho_dS fit zs**5 term [kg m-3 ppt-1]
real, parameter :: BET010 = 0.5*EOS110*r1_S0  ! Coefficient of the drho_dS fit zt term [kg m-3 ppt-1]
real, parameter :: BET110 =     EOS210*r1_S0  ! Coefficient of the drho_dS fit zs * zt term [kg m-3 ppt-1]
real, parameter :: BET210 = 1.5*EOS310*r1_S0  ! Coefficient of the drho_dS fit zs**2 * zt term [kg m-3 ppt-1]
real, parameter :: BET310 = 2.0*EOS410*r1_S0  ! Coefficient of the drho_dS fit zs**3 * zt term [kg m-3 ppt-1]
real, parameter :: BET410 = 2.5*EOS510*r1_S0  ! Coefficient of the drho_dS fit zs**4 * zt term [kg m-3 ppt-1]
real, parameter :: BET020 = 0.5*EOS120*r1_S0  ! Coefficient of the drho_dS fit zt**2 term [kg m-3 ppt-1]
real, parameter :: BET120 =     EOS220*r1_S0  ! Coefficient of the drho_dS fit zs * zt**2 term [kg m-3 ppt-1]
real, parameter :: BET220 = 1.5*EOS320*r1_S0  ! Coefficient of the drho_dS fit zs**2 * zt**2 term [kg m-3 ppt-1]
real, parameter :: BET320 = 2.0*EOS420*r1_S0  ! Coefficient of the drho_dS fit zs**3 * zt**2 term [kg m-3 ppt-1]
real, parameter :: BET030 = 0.5*EOS130*r1_S0  ! Coefficient of the drho_dS fit zt**3 term [kg m-3 ppt-1]
real, parameter :: BET130 =     EOS230*r1_S0  ! Coefficient of the drho_dS fit zs * zt**3 term [kg m-3 ppt-1]
real, parameter :: BET230 = 1.5*EOS330*r1_S0  ! Coefficient of the drho_dS fit zs**2 * zt**3 term [kg m-3 ppt-1]
real, parameter :: BET040 = 0.5*EOS140*r1_S0  ! Coefficient of the drho_dS fit zt**4 term [kg m-3 ppt-1]
real, parameter :: BET140 =     EOS240*r1_S0  ! Coefficient of the drho_dS fit zs * zt**4 term [kg m-3 ppt-1]
real, parameter :: BET050 = 0.5*EOS150*r1_S0  ! Coefficient of the drho_dS fit zt**5 term [kg m-3 ppt-1]
real, parameter :: BET001 = 0.5*EOS101*r1_S0  ! Coefficient of the drho_dS fit zp term [kg m-3 ppt-1]
real, parameter :: BET101 =     EOS201*r1_S0  ! Coefficient of the drho_dS fit zs * zp term [kg m-3 ppt-1]
real, parameter :: BET201 = 1.5*EOS301*r1_S0  ! Coefficient of the drho_dS fit zs**2 * zp term [kg m-3 ppt-1]
real, parameter :: BET301 = 2.0*EOS401*r1_S0  ! Coefficient of the drho_dS fit zs**3 * zp term [kg m-3 ppt-1]
real, parameter :: BET011 = 0.5*EOS111*r1_S0  ! Coefficient of the drho_dS fit zt * zp term [kg m-3 ppt-1]
real, parameter :: BET111 =     EOS211*r1_S0  ! Coefficient of the drho_dS fit zs * zt * zp term [kg m-3 ppt-1]
real, parameter :: BET211 = 1.5*EOS311*r1_S0  ! Coefficient of the drho_dS fit zs**2 * zt * zp term [kg m-3 ppt-1]
real, parameter :: BET021 = 0.5*EOS121*r1_S0  ! Coefficient of the drho_dS fit zt**2 * zp term [kg m-3 ppt-1]
real, parameter :: BET121 =     EOS221*r1_S0  ! Coefficient of the drho_dS fit zs * zt**2 * zp term [kg m-3 ppt-1]
real, parameter :: BET031 = 0.5*EOS131*r1_S0  ! Coefficient of the drho_dS fit zt**3 * zp term [kg m-3 ppt-1]
real, parameter :: BET002 = 0.5*EOS102*r1_S0  ! Coefficient of the drho_dS fit zp**2 term [kg m-3 ppt-1]
real, parameter :: BET102 =     EOS202*r1_S0  ! Coefficient of the drho_dS fit zs * zp**2 term [kg m-3 ppt-1]
real, parameter :: BET012 = 0.5*EOS112*r1_S0  ! Coefficient of the drho_dS fit zt * zp**2 term [kg m-3 ppt-1]
real, parameter :: BET003 = 0.5*EOS103*r1_S0  ! Coefficient of the drho_dS fit zp**3 term [kg m-3 ppt-1]
!>@}

contains

!> This subroutine computes the in situ density of sea water (rho in
!! [kg m-3]) from absolute salinity (S [g kg-1]), conservative temperature
!! (T [degC]), and pressure [Pa].  It uses the expressions derived for use
!! with NEMO.
subroutine calculate_density_scalar_nemo(T, S, pressure, rho, rho_ref)
  real,           intent(in)  :: T        !< Conservative temperature [degC].
  real,           intent(in)  :: S        !< Absolute salinity [g kg-1].
  real,           intent(in)  :: pressure !< pressure [Pa].
  real,           intent(out) :: rho      !< In situ density [kg m-3].
  real, optional, intent(in)  :: rho_ref  !< A reference density [kg m-3].

  real, dimension(1) :: T0    ! A 1-d array with a copy of the conservative temperature [degC]
  real, dimension(1) :: S0    ! A 1-d array with a copy of the absolute salinity [g kg-1]
  real, dimension(1) :: pressure0 ! A 1-d array with a copy of the pressure [Pa]
  real, dimension(1) :: rho0  ! A 1-d array with a copy of the density [kg m-3]

  T0(1) = T
  S0(1) = S
  pressure0(1) = pressure

  call calculate_density_array_nemo(T0, S0, pressure0, rho0, 1, 1, rho_ref)
  rho = rho0(1)

end subroutine calculate_density_scalar_nemo

!> This subroutine computes the in situ density of sea water (rho in
!! [kg m-3]) from absolute salinity (S [g kg-1]), conservative temperature
!! (T [degC]), and pressure [Pa].  It uses the expressions derived for use
!! with NEMO.
subroutine calculate_density_array_nemo(T, S, pressure, rho, start, npts, rho_ref)
  real, dimension(:), intent(in)  :: T        !< Conservative temperature [degC].
  real, dimension(:), intent(in)  :: S        !< Absolute salinity [g kg-1].
  real, dimension(:), intent(in)  :: pressure !< pressure [Pa].
  real, dimension(:), intent(out) :: rho      !< in situ density [kg m-3].
  integer,            intent(in)  :: start    !< the starting point in the arrays.
  integer,            intent(in)  :: npts     !< the number of values to calculate.
  real,     optional, intent(in)  :: rho_ref  !< A reference density [kg m-3].

  ! Local variables
  real :: zp  ! Pressure, first in [dbar], then normalized by an assumed pressure range [nondim]
  real :: zt  ! Conservative temperature, first in [degC], then normalized by an assumed temperature range [nondim]
  real :: zs  ! Absolute salinity, first in [g kg-1], then the square root of salinity with an offset normalized
              ! by an assumed salnity range [nondim]
  real :: zr0 ! A pressure-dependent but temperature and salinity independent contribution to
              ! density at the reference temperature and salinity [kg m-3]
  real :: zn  ! Density without a pressure-dependent contribution [kg m-3]
  real :: zn0 ! A contribution to density from temperature and salinity anomalies at the surface pressure [kg m-3]
  real :: zn1 ! A temperature and salinity dependent density contribution proportional to pressure [kg m-3]
  real :: zn2 ! A temperature and salinity dependent density contribution proportional to pressure^2 [kg m-3]
  real :: zn3 ! A temperature and salinity dependent density contribution proportional to pressure^3 [kg m-3]
  real :: zs0 ! Salinity dependent density at the surface pressure and temperature [kg m-3]
  integer :: j

  ! The following algorithm was published by Roquet et al. (2015), intended for use
  ! with NEMO, but it is not necessarily the algorithm used in NEMO ocean model.
  do j=start,start+npts-1
    ! Conversions to the units used here.
    zt = T(j) * r1_T0  ! Conservative temperature normalized by a plausible oceanic range [nondim]
    zs = SQRT( ABS( S(j) + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
    zp = pressure(j) * (Pa2db*r1_P0)    ! Convert pressure from Pascals to kilobars to normalize it [nondim]

    ! The next two lines should be used if it is necessary to convert potential temperature and
    ! pratical salinity to conservative temperature and absolute salinity.
    ! zt = r1_T0 * gsw_ct_from_pt(S(j),T(j)) ! Convert potential temp to conservative temp [degC]
    ! zs = SQRT( ABS( gsw_sr_from_sp(S(j)) + rdeltaS ) * r1_S0 ) ! Convert S from practical to absolute salinity.

    zn3 = EOS013*zt   &
       &   + EOS103*zs+EOS003

    zn2 = (EOS022*zt   &
       &   + EOS112*zs+EOS012)*zt   &
       &   + (EOS202*zs+EOS102)*zs+EOS002

    zn1 = (((EOS041*zt   &
       &   + EOS131*zs+EOS031)*zt   &
       &   + (EOS221*zs+EOS121)*zs+EOS021)*zt   &
       &   + ((EOS311*zs+EOS211)*zs+EOS111)*zs+EOS011)*zt   &
       &   + (((EOS401*zs+EOS301)*zs+EOS201)*zs+EOS101)*zs+EOS001

    zn0 = (((((EOS060*zt   &
       &   + EOS150*zs+EOS050)*zt   &
       &   + (EOS240*zs+EOS140)*zs+EOS040)*zt   &
       &   + ((EOS330*zs+EOS230)*zs+EOS130)*zs+EOS030)*zt   &
       &   + (((EOS420*zs+EOS320)*zs+EOS220)*zs+EOS120)*zs+EOS020)*zt   &
       &   + ((((EOS510*zs+EOS410)*zs+EOS310)*zs+EOS210)*zs+EOS110)*zs+EOS010)*zt

    zs0 = (((((EOS600*zs+EOS500)*zs+EOS400)*zs+EOS300)*zs+EOS200)*zs+EOS100)*zs + EOS000

    zr0 = (((((R05 * zp+R04) * zp+R03 ) * zp+R02 ) * zp+R01) * zp+R00) * zp

    if (present(rho_ref)) then
      zn  = ( ( zn3 * zp + zn2 ) * zp + zn1 ) * zp + (zn0 + (zs0 - rho_ref))
      rho(j) =  ( zn + zr0 ) ! density
    else
      zn  = ( ( zn3 * zp + zn2 ) * zp + zn1 ) * zp + (zn0 + zs0)
      rho(j) =  ( zn + zr0 ) ! density
    endif

  enddo
end subroutine calculate_density_array_nemo

!> For a given thermodynamic state, calculate the derivatives of density with conservative
!! temperature and absolute salinity, using the expressions derived for use with NEMO.
subroutine calculate_density_derivs_array_nemo(T, S, pressure, drho_dT, drho_dS, start, npts)
  real,    intent(in),  dimension(:) :: T        !< Conservative temperature [degC].
  real,    intent(in),  dimension(:) :: S        !< Absolute salinity [g kg-1].
  real,    intent(in),  dimension(:) :: pressure !< pressure [Pa].
  real,    intent(out), dimension(:) :: drho_dT  !< The partial derivative of density with potential
                                                 !! temperature [kg m-3 degC-1].
  real,    intent(out), dimension(:) :: drho_dS  !< The partial derivative of density with salinity,
                                                 !! in [kg m-3 ppt-1].
  integer, intent(in)                :: start    !< The starting point in the arrays.
  integer, intent(in)                :: npts     !< The number of values to calculate.

  ! Local variables
  real :: zp  ! Pressure, first in [dbar], then normalized by an assumed pressure range [nondim]
  real :: zt  ! Conservative temperature, first in [degC], then normalized by an assumed temperature range [nondim]
  real :: zs  ! Absolute salinity, first in [g kg-1], then the square root of salinity with an offset normalized
              ! by an assumed salnity range [nondim]
  real :: zn  ! Partial derivative of density with temperature [kg m-3 degC-1] or salinity [kg m-3 ppt-1]
              ! without a pressure-dependent contribution
  real :: zn0 ! A contribution to the partial derivative of density with temperature [kg m-3 degC-1] or
              ! salinity [kg m-3 ppt-1] from temperature anomalies at the surface pressure
  real :: zn1 ! A contribution to the partial derivative of density with temperature [kg m-3 degC-1] or
              ! salinity [kg m-3 ppt-1] proportional to pressure
  real :: zn2 ! A contribution to the partial derivative of density with temperature [kg m-3 degC-1] or
              ! salinity [kg m-3 ppt-1] proportional to pressure^2
  real :: zn3 ! A contribution to the partial derivative of density with temperature [kg m-3 degC-1] or
              ! salinity [kg m-3 ppt-1] proportional to pressure^3
  integer :: j

  do j=start,start+npts-1
    ! Conversions to the units used here.
    zt = T(j) * r1_T0  ! Conservative temperature normalized by a plausible oceanic range [nondim]
    zs = SQRT( ABS( S(j) + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
    zp = pressure(j) * (Pa2db*r1_P0)    ! Convert pressure from Pascals to kilobars to normalize it [nondim]

    ! The next two lines should be used if it is necessary to convert potential temperature and
    ! pratical salinity to conservative temperature and absolute salinity.
    ! zt = r1_T0 * gsw_ct_from_pt(S(j),T(j)) ! Convert potential temp to conservative temp [degC]
    ! zs = SQRT( ABS( gsw_sr_from_sp(S(j)) + rdeltaS ) * r1_S0 ) ! Convert S from practical to absolute salinity.

    !
    ! alpha
    zn3 = ALP003
    !
    zn2 = ALP012*zt + ALP102*zs+ALP002
    !
    zn1 = ((ALP031*zt   &
       &   + ALP121*zs+ALP021)*zt   &
       &   + (ALP211*zs+ALP111)*zs+ALP011)*zt   &
       &   + ((ALP301*zs+ALP201)*zs+ALP101)*zs+ALP001
       !
    zn0 = ((((ALP050*zt   &
       &   + ALP140*zs+ALP040)*zt   &
       &   + (ALP230*zs+ALP130)*zs+ALP030)*zt   &
       &   + ((ALP320*zs+ALP220)*zs+ALP120)*zs+ALP020)*zt   &
       &   + (((ALP410*zs+ALP310)*zs+ALP210)*zs+ALP110)*zs+ALP010)*zt   &
       &   + ((((ALP500*zs+ALP400)*zs+ALP300)*zs+ALP200)*zs+ALP100)*zs+ALP000
       !
    zn  = ( ( zn3 * zp + zn2 ) * zp + zn1 ) * zp + zn0
    !
    drho_dT(j) = zn
    !
    ! beta
    !
    zn3 = BET003
    !
    zn2 = BET012*zt + BET102*zs+BET002
    !
    zn1 = ((BET031*zt   &
       &   + BET121*zs+BET021)*zt   &
       &   + (BET211*zs+BET111)*zs+BET011)*zt   &
       &   + ((BET301*zs+BET201)*zs+BET101)*zs+BET001
       !
    zn0 = ((((BET050*zt   &
       &   + BET140*zs+BET040)*zt   &
       &   + (BET230*zs+BET130)*zs+BET030)*zt   &
       &   + ((BET320*zs+BET220)*zs+BET120)*zs+BET020)*zt   &
       &   + (((BET410*zs+BET310)*zs+BET210)*zs+BET110)*zs+BET010)*zt   &
       &   + ((((BET500*zs+BET400)*zs+BET300)*zs+BET200)*zs+BET100)*zs+BET000
       !
    zn  = ( ( zn3 * zp + zn2 ) * zp + zn1 ) * zp + zn0

    ! The division by zs here is because zs = sqrt(S + S0), so drho_dS = dzs_dS * drho_dzs = (0.5 / zs) * drho_dzs
    drho_dS(j) = zn / zs
  enddo

end subroutine calculate_density_derivs_array_nemo

!> Wrapper to calculate_density_derivs_array for scalar inputs
subroutine calculate_density_derivs_scalar_nemo(T, S, pressure, drho_dt, drho_ds)
  real,    intent(in)  :: T        !< Potential temperature relative to the surface [degC].
  real,    intent(in)  :: S        !< Salinity [g kg-1].
  real,    intent(in)  :: pressure !< Pressure [Pa].
  real,    intent(out) :: drho_dT  !< The partial derivative of density with potential
                                   !! temperature [kg m-3 degC-1].
  real,    intent(out) :: drho_dS  !< The partial derivative of density with salinity,
                                   !! in [kg m-3 ppt-1].
  ! Local variables
  real, dimension(1) :: T0    ! A 1-d array with a copy of the conservative temperature [degC]
  real, dimension(1) :: S0    ! A 1-d array with a copy of the absolute salinity [g kg-1]
  real, dimension(1) :: pressure0 ! A 1-d array with a copy of the pressure [Pa]
  real, dimension(1) :: drdt0 ! A 1-d array with a copy of the derivative of density
                              ! with potential temperature [kg m-3 degC-1]
  real, dimension(1) :: drds0 ! A 1-d array with a copy of the derivative of density
                              ! with salinity [kg m-3 ppt-1]

  T0(1) = T
  S0(1) = S
  pressure0(1) = pressure

  call calculate_density_derivs_array_nemo(T0, S0, pressure0, drdt0, drds0, 1, 1)
  drho_dt = drdt0(1)
  drho_ds = drds0(1)
end subroutine calculate_density_derivs_scalar_nemo

!> Compute the in situ density of sea water (rho in [kg m-3]) and the compressibility
!! (drho/dp = C_sound^-2, stored as drho_dp [s2 m-2]) from absolute salinity (sal [g kg-1]),
!! conservative temperature (T [degC]), and pressure [Pa], using the expressions
!! derived for use with NEMO.
subroutine calculate_compress_nemo(T, S, pressure, rho, drho_dp, start, npts)
  real,    intent(in),  dimension(:) :: T        !< Conservative temperature [degC].
  real,    intent(in),  dimension(:) :: S        !< Absolute salinity [g kg-1].
  real,    intent(in),  dimension(:) :: pressure !< pressure [Pa].
  real,    intent(out), dimension(:) :: rho      !< In situ density [kg m-3].
  real,    intent(out), dimension(:) :: drho_dp  !< The partial derivative of density with pressure
                                                 !! (also the inverse of the square of sound speed)
                                                 !! [s2 m-2].
  integer, intent(in)                :: start    !< The starting point in the arrays.
  integer, intent(in)                :: npts     !< The number of values to calculate.

  ! Local variables
  real :: zp  ! Pressure normalized by an assumed pressure range [nondim]
  real :: zt  ! Conservative temperature normalized by an assumed temperature range [nondim]
  real :: zs  ! The square root of absolute salinity with an offset normalized
              ! by an assumed salnity range [nondim]
  real :: dzr0_dp ! Derivative of the pressure-dependent reference density profile with normalized pressure [kg m-3]
  real :: dzn_dp  ! Derivative of the density anomaly from the reference profile with normalized pressure [kg m-3]
  real :: zr0 ! The pressure-dependent (but temperature and salinity independent) reference density profile [kg m-3]
  real :: zn  ! Density anomaly from the reference profile [kg m-3]
  real :: zn0 ! A contribution to density from temperature and salinity anomalies at the surface pressure [kg m-3]
  real :: zn1 ! A temperature and salinity dependent density contribution proportional to pressure [kg m-3]
  real :: zn2 ! A temperature and salinity dependent density contribution proportional to pressure^2 [kg m-3]
  real :: zn3 ! A temperature and salinity dependent density contribution proportional to pressure^3 [kg m-3]
  real :: zs0 ! Salinity dependent density at the surface pressure and temperature [kg m-3]
  integer :: j

  ! The following algorithm was published by Roquet et al. (2015), intended for use
  ! with NEMO, but it is not necessarily the algorithm used in NEMO ocean model.
  do j=start,start+npts-1
    ! Conversions to the units used here.
    zt = T(j) * r1_T0  ! Conservative temperature normalized by a plausible oceanic range [nondim]
    zs = SQRT( ABS( S(j) + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
    zp = pressure(j) * (Pa2db*r1_P0)    ! Convert pressure from Pascals to kilobars to normalize it [nondim]

    ! The next two lines should be used if it is necessary to convert potential temperature and
    ! pratical salinity to conservative temperature and absolute salinity.
    ! zt = r1_T0 * gsw_ct_from_pt(S(j),T(j)) ! Convert potential temp to conservative temp [degC]
    ! zs = SQRT( ABS( gsw_sr_from_sp(S(j)) + rdeltaS ) * r1_S0 ) ! Convert S from practical to absolute salinity.

    zn3 = EOS013*zt + EOS103*zs + EOS003

    zn2 = (EOS022*zt   &
       &   + EOS112*zs + EOS012)*zt   &
       &   + (EOS202*zs + EOS102)*zs + EOS002

    zn1 = (((EOS041*zt   &
       &   + EOS131*zs + EOS031)*zt   &
       &   + (EOS221*zs + EOS121)*zs + EOS021)*zt   &
       &   + ((EOS311*zs + EOS211)*zs + EOS111)*zs + EOS011)*zt   &
       &   + (((EOS401*zs + EOS301)*zs + EOS201)*zs + EOS101)*zs + EOS001

    zn0 = (((((EOS060*zt   &
       &   + EOS150*zs + EOS050)*zt   &
       &   + (EOS240*zs + EOS140)*zs + EOS040)*zt   &
       &   + ((EOS330*zs + EOS230)*zs + EOS130)*zs + EOS030)*zt   &
       &   + (((EOS420*zs + EOS320)*zs + EOS220)*zs + EOS120)*zs + EOS020)*zt   &
       &   + ((((EOS510*zs + EOS410)*zs + EOS310)*zs + EOS210)*zs + EOS110)*zs + EOS010)*zt

    zs0 = (((((EOS600*zs + EOS500)*zs + EOS400)*zs + EOS300)*zs + EOS200)*zs + EOS100)*zs + EOS000

    zr0 = (((((R05*zp + R04)*zp + R03)*zp + R02)*zp + R01)*zp + R00)*zp

    zn  = ( ( zn3*zp + zn2 )*zp + zn1 )*zp + (zn0 + zs0)
    rho(j) =  ( zn + zr0 ) ! density

    dzr0_dp = ((((6.*R05*zp + 5.*R04)*zp + 4.*R03)*zp + 3.*R02)*zp + 2.*R01)*zp + R00
    dzn_dp  = ( 3.*zn3*zp + 2.*zn2 )*zp + zn1
    drho_dp(j) =  ( dzn_dp + dzr0_dp ) * (Pa2db*r1_P0) ! density

  enddo
end subroutine calculate_compress_nemo


!> Second derivatives of density with respect to temperature, salinity, and pressure for 1-d array inputs and outputs.
subroutine calculate_density_second_derivs_array_NEMO(T, S, P, drho_ds_ds, drho_ds_dt, drho_dt_dt, &
                                                         drho_ds_dp, drho_dt_dp, start, npts)
  real, dimension(:), intent(in   ) :: T !< Potential temperature referenced to 0 dbar [degC]
  real, dimension(:), intent(in   ) :: S !< Salinity [PSU]
  real, dimension(:), intent(in   ) :: P !< Pressure [Pa]
  real, dimension(:), intent(inout) :: drho_ds_ds !< Partial derivative of beta with respect
                                                  !! to S [kg m-3 PSU-2]
  real, dimension(:), intent(inout) :: drho_ds_dt !< Partial derivative of beta with respect
                                                  !! to T [kg m-3 PSU-1 degC-1]
  real, dimension(:), intent(inout) :: drho_dt_dt !< Partial derivative of alpha with respect
                                                  !! to T [kg m-3 degC-2]
  real, dimension(:), intent(inout) :: drho_ds_dp !< Partial derivative of beta with respect
                                                  !! to pressure [kg m-3 PSU-1 Pa-1] = [s2 m-2 PSU-1]
  real, dimension(:), intent(inout) :: drho_dt_dp !< Partial derivative of alpha with respect
                                                  !! to pressure [kg m-3 degC-1 Pa-1] = [s2 m-2 degC-1]
  integer,            intent(in   ) :: start !< Starting index in T,S,P
  integer,            intent(in   ) :: npts  !< Number of points to loop over

  ! Local variables
  real :: zp  ! Pressure normalized by an assumed pressure range [nondim]
  real :: zt  ! Conservative temperature normalized by an assumed temperature range [nondim]
  real :: zs  ! The square root of absolute salinity with an offset normalized
              ! by an assumed salnity range [nondim]
  real :: I_s ! The inverse of zs [nondim]
  real :: dzr0_dp ! Derivative of the pressure-dependent reference density profile with normalized pressure [kg m-3]
  real :: dzn_dp  ! Derivative of the density anomaly from the reference profile with normalized pressure [kg m-3]
  real :: dzn_ds  ! Derivative of the density anomaly from the reference profile with zs [kg m-3]
  real :: zr0 ! The pressure-dependent (but temperature and salinity independent) reference density profile [kg m-3]
  real :: zn  ! Density anomaly from the reference profile [kg m-3]
  real :: zn0 ! A contribution to one of the second derivatives that is independent of pressure [various]
  real :: zn1 ! A contribution to one of the second derivatives that is proportional to pressure [various]
  real :: zn2 ! A contribution to one of the second derivatives that is proportional to pressure^2 [various]
  real :: zn3 ! A temperature and salinity dependent density contribution proportional to pressure^3 [various]
  integer :: j

  do j = start,start+npts-1
    ! Conversions to the units used here.
    zt = T(j) * r1_T0  ! Conservative temperature normalized by a plausible oceanic range [nondim]
    zs = SQRT( ABS( S(j) + rdeltaS ) * r1_S0 )  ! square root of normalized salinity plus an offset [nondim]
    zp = P(j) * (Pa2db*r1_P0)     ! Convert pressure from Pascals to kilobars to normalize it [nondim]

    ! The next two lines should be used if it is necessary to convert potential temperature and
    ! pratical salinity to conservative temperature and absolute salinity.
    ! zt = r1_T0 * gsw_ct_from_pt(S(j),T(j)) ! Convert potential temp to conservative temp [degC]
    ! zs = SQRT( ABS( gsw_sr_from_sp(S(j)) + rdeltaS ) * r1_S0 )  ! Convert S from practical to absolute salinity.

    I_s = 1.0 / zs

    ! Find drho_ds_ds
    zn3 = -EOS103*I_s**2
    zn2 = -(EOS112*zt + EOS102)*I_s**2
    zn1 = (3.*EOS311*zt + (8.*EOS401*zs + 3.*EOS301) ) &
         - ( ((EOS131*zt + EOS121)*zt + EOS111)*zt + EOS101 )*I_s**2
    zn0 = ( (( 3.*EOS330*zt + (8.*EOS420*zs + 3.*EOS320))*zt + &
             ((15.*EOS510*zs + 8.*EOS410)*zs + 3.*EOS310))*zt + &
            (((24.*EOS600*zs + 15.*EOS500)*zs + 8.*EOS400)*zs + 3.*EOS300) ) &
          - ( ((((EOS150*zt + EOS140)*zt + EOS130)*zt + EOS120)*zt + EOS110)*zt + EOS100 )*I_s**2
    zn  = ( ( zn3 * zp + zn2) * zp + zn1 ) * zp + zn0
    drho_dS_dS(j) = (0.5*r1_S0)**2 * (zn * I_s)

    ! Find drho_ds_dt
    zn2 = EOS112
    zn1 = ((3.*EOS131)*zt  + (4.*EOS221*zs + 2.*EOS121))*zt + &
          ((3.*EOS311*zs + 2.*EOS211)*zs + EOS111)
    zn0 = (((5.*EOS150*zt + (8.*EOS240*zs + 4.*EOS140))*zt + &
            ((9.*EOS330*zs + 6.*EOS230)*zs + 3.*EOS130))*zt + &
           ((((8.*EOS420*zs + 6.*EOS320)*zs + 4.*EOS220)*zs + 2.*EOS120)))*zt +  &
          ((((5.*EOS510*zs + 4.*EOS410)*zs + 3.*EOS310)*zs + 2.*EOS210)*zs + EOS110)
    zn  = ( zn2 * zp + zn1 ) * zp + zn0
    drho_ds_dt(j) = (0.5*r1_S0*r1_T0) * (zn * I_s)

    ! Find drho_dt_dt
    zn2 = 2.*EOS022
    zn1 = (12.*EOS041*zt  + 6.*(EOS131*zs + EOS031))*zt +  &
          2.*((EOS221*zs + EOS121)*zs + EOS021)
    zn0 = (((30.*EOS060*zt + 20.*(EOS150*zs + EOS050))*zt + &
            12.*((EOS240*zs + EOS140)*zs + EOS040))*zt  + &
           6.*(((EOS330*zs + EOS230)*zs + EOS130)*zs + EOS030))*zt +  &
          2.*((((EOS420*zs + EOS320)*zs + EOS220)*zs + EOS120)*zs + EOS020)
    zn  = ( zn2 * zp + zn1 ) * zp + zn0
    drho_dt_dt(j) = zn * r1_T0**2

    ! Find drho_ds_dp
    zn3 = EOS103
    zn2 = EOS112*zt + (2.*EOS202*zs + EOS102)
    zn1 = ((EOS131*zt + (2.*EOS221*zs + EOS121))*zt  + ((3.*EOS311*zs + 2.*EOS211)*zs + EOS111))*zt + &
          (((4.*EOS401*zs + 3.*EOS301)*zs + 2.*EOS201)*zs + EOS101)
    dzn_dp  = ( ( 3.*zn3 * zp + 2.*zn2 ) * zp + zn1 )
    drho_ds_dp(j) =  ( dzn_dp * I_s ) * (0.5*r1_S0 * Pa2db*r1_P0) ! Second derivative of density


    ! Find drho_dt_dp
    zn3 = EOS013
    zn2 = 2.*EOS022*zt + (EOS112*zs + EOS012)
    zn1 = ((4.*EOS041*zt  + 3.*(EOS131*zs + EOS031))*zt + 2.*((EOS221*zs + EOS121)*zs + EOS021))*zt + &
          (((EOS311*zs + EOS211)*zs + EOS111)*zs + EOS011)
    dzn_dp  = ( ( 3.*zn3 * zp + 2.*zn2 ) * zp + zn1 )
    drho_dt_dp(j) =  ( dzn_dp ) * (Pa2db*r1_P0* r1_T0) ! Second derivative of density
  enddo

end subroutine calculate_density_second_derivs_array_NEMO

!> Second derivatives of density with respect to temperature, salinity, and pressure for scalar inputs.
!!
!! The scalar version of calculate_density_second_derivs promotes scalar inputs to 1-element array
!! and then demotes the output back to a scalar
subroutine calculate_density_second_derivs_scalar_NEMO(T, S, P, drho_ds_ds, drho_ds_dt, drho_dt_dt, &
                                                         drho_ds_dp, drho_dt_dp)
  real, intent(in   ) :: T          !< Potential temperature referenced to 0 dbar
  real, intent(in   ) :: S          !< Salinity [PSU]
  real, intent(in   ) :: P          !< pressure [Pa]
  real, intent(  out) :: drho_ds_ds !< Partial derivative of beta with respect
                                    !! to S [kg m-3 PSU-2]
  real, intent(  out) :: drho_ds_dt !< Partial derivative of beta with respect
                                    !! to T [kg m-3 PSU-1 degC-1]
  real, intent(  out) :: drho_dt_dt !< Partial derivative of alpha with respect
                                    !! to T [kg m-3 degC-2]
  real, intent(  out) :: drho_ds_dp !< Partial derivative of beta with respect
                                    !! to pressure [kg m-3 PSU-1 Pa-1] = [s2 m-2 PSU-1]
  real, intent(  out) :: drho_dt_dp !< Partial derivative of alpha with respect
                                    !! to pressure [kg m-3 degC-1 Pa-1] = [s2 m-2 degC-1]
  ! Local variables
  real, dimension(1) :: T0    ! A 1-d array with a copy of the temperature [degC]
  real, dimension(1) :: S0    ! A 1-d array with a copy of the salinity [PSU]
  real, dimension(1) :: p0    ! A 1-d array with a copy of the pressure [Pa]
  real, dimension(1) :: drdsds ! The second derivative of density with salinity [kg m-3 PSU-2]
  real, dimension(1) :: drdsdt ! The second derivative of density with salinity and
                               ! temperature [kg m-3 PSU-1 degC-1]
  real, dimension(1) :: drdtdt ! The second derivative of density with temperature [kg m-3 degC-2]
  real, dimension(1) :: drdsdp ! The second derivative of density with salinity and
                               ! pressure [kg m-3 PSU-1 Pa-1] = [s2 m-2 PSU-1]
  real, dimension(1) :: drdtdp ! The second derivative of density with temperature and
                               ! pressure [kg m-3 degC-1 Pa-1] = [s2 m-2 degC-1]

  T0(1) = T
  S0(1) = S
  P0(1) = P
  call calculate_density_second_derivs_array_NEMO(T0, S0, P0, drdsds, drdsdt, drdtdt, drdsdp, drdtdp, 1, 1)
  drho_ds_ds = drdsds(1)
  drho_ds_dt = drdsdt(1)
  drho_dt_dt = drdtdt(1)
  drho_ds_dp = drdsdp(1)
  drho_dt_dp = drdtdp(1)

end subroutine calculate_density_second_derivs_scalar_NEMO

!> \namespace mom_eos_NEMO
!!
!! \section section_EOS_NEMO NEMO equation of state
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
!! The NEMO label used to describe this equation of state reflects that it was used in the NEMO
!! ocean model before it was used in MOM6, but it probably should be described as the Roquet
!! equation of.   However, these algorithms, especially as modified here, are not from
!! the standard NEMO codebase.
!!
!! \subsection section_EOS_NEMO_references References
!!
!! Roquet, F., Madec, G., McDougall, T. J., and Barker, P. M., 2015:
!!  Accurate polynomial expressions for the density and specific volume
!!  of seawater using the TEOS-10 standard. Ocean Modelling, 90:29-43.

end module MOM_EOS_NEMO
