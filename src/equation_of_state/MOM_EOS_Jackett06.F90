!> The equation of state using the Jackett et al 2006 expressions that are often used in Hycom
module MOM_EOS_Jackett06

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_hor_index, only : hor_index_type

implicit none ; private

public calculate_compress_Jackett06, calculate_density_Jackett06, calculate_spec_vol_Jackett06
public calculate_density_derivs_Jackett06, calculate_specvol_derivs_Jackett06
public calculate_density_second_derivs_Jackett06, EoS_fit_range_Jackett06

!> Compute the in situ density of sea water (in [kg m-3]), or its anomaly with respect to
!! a reference density, from salinity in practical salinity units ([PSU]), potential
!! temperature (in degrees Celsius [degC]), and pressure [Pa], using the expressions from
!! Jackett et al., 2006, J. Atmos. Ocean. Tech., 32, 1709-1728.
interface calculate_density_Jackett06
  module procedure calculate_density_scalar_Jackett, calculate_density_array_Jackett
end interface calculate_density_Jackett06

!> Compute the in situ specific volume of sea water (in [m3 kg-1]), or an anomaly with respect
!! to a reference specific volume, from salinity in practical salinity units ([PSU]), potential
!! temperature (in degrees Celsius [degC]), and pressure [Pa], using the expressions from
!! Jackett et al., 2006, J. Atmos. Ocean. Tech., 32, 1709-1728.
interface calculate_spec_vol_Jackett06
  module procedure calculate_spec_vol_scalar_Jackett, calculate_spec_vol_array_Jackett
end interface calculate_spec_vol_Jackett06

!> Compute the derivatives of density with temperature and salinity
interface calculate_density_derivs_Jackett06
  module procedure calculate_density_derivs_scalar_Jackett, calculate_density_derivs_array_Jackett
end interface calculate_density_derivs_Jackett06

!> Compute the second derivatives of density with various combinations
!! of temperature, salinity, and pressure
interface calculate_density_second_derivs_Jackett06
  module procedure calculate_density_second_derivs_scalar_Jackett, calculate_density_second_derivs_array_Jackett
end interface calculate_density_second_derivs_Jackett06

!>@{ Parameters in the Jackett et al. equation of state, which is a fit to the Fiestel (2003)
!    equation of state for the range: -2 < theta < 40 [degC], 0 < S < 42 [PSU], 0 < p < 1e8 [Pa].
!    The notation here is for terms in the numerator of the expression for density of
!    RNabc for terms proportional to S**a * T**b * P**c, and terms in the denominator as RDabc.
!    For terms proportional to S**1.5, 6 is used in this notation.

! --- coefficients for 25-term rational function sigloc().
real, parameter :: &
  RN000 =  9.9984085444849347d+02, & ! Density numerator constant coefficient [kg m-3]
  RN001 =  1.1798263740430364d-06, & ! Density numerator P       coefficient [kg m-3 Pa-1]
  RN002 = -2.5862187075154352d-16, & ! Density numerator P^2     coefficient [kg m-3 Pa-2]
  RN010 =  7.3471625860981584d+00, & ! Density numerator T       coefficient [kg m-3 degC-1]
  RN020 = -5.3211231792841769d-02, & ! Density numerator T^2     coefficient [kg m-3 degC-2]
  RN021 =  9.8920219266399117d-12, & ! Density numerator T^2 P   coefficient [kg m-3 degC-2 Pa-1]
  RN022 = -3.2921414007960662d-20, & ! Density numerator T^2 P^2 coefficient [kg m-3 degC-2 Pa-2]
  RN030 =  3.6492439109814549d-04, & ! Density numerator T^3     coefficient [kg m-3 degC-3]
  RN100 =  2.5880571023991390d+00, & ! Density numerator S       coefficient [kg m-3 PSU-1]
  RN101 =  4.6996642771754730d-10, & ! Density numerator S P     coefficient [kg m-3 PSU-1 Pa-1]
  RN110 = -6.7168282786692355d-03, & ! Density numerator S T     coefficient [kg m-3 degC-1 PSU-1]
  RN200 =  1.9203202055760151d-03, & ! Density numerator S^2      coefficient [kg m-3]

  RD001 =  6.7103246285651894d-10, & ! Density denominator P       coefficient [Pa-1]
  RD010 =  7.2815210113327091d-03, & ! Density denominator T       coefficient [degC-1]
  RD013 = -9.1534417604289062d-30, & ! Density denominator T P^3   coefficient [degC-1 Pa-3]
  RD020 = -4.4787265461983921d-05, & ! Density denominator T^2     coefficient [degC-2]
  RD030 =  3.3851002965802430d-07, & ! Density denominator T^3     coefficient [degC-3]
  RD032 = -2.4461698007024582d-25, & ! Density denominator T^3 P^2 coefficient [degC-3  Pa-2]
  RD040 =  1.3651202389758572d-10, & ! Density denominator T^4     coefficient [degC-4]
  RD100 =  1.7632126669040377d-03, & ! Density denominator S       coefficient [PSU-1]
  RD110 = -8.8066583251206474d-06, & ! Density denominator S T     coefficient [degC-1 PSU-1]
  RD130 = -1.8832689434804897d-10, & ! Density denominator S T^3   coefficient [degC-3 PSU-1]
  RD600 =  5.7463776745432097d-06, & ! Density denominator S^1.5   coefficient [PSU-1.5]
  RD620 =  1.4716275472242334d-09    ! Density denominator S^1.5 T^2 coefficient [degC-2 PSU-1.5]
!>@}

contains

!> Computes the in situ density of sea water for 1-d array inputs and outputs.
!!
!! Returns the in situ density of sea water (rho in [kg m-3]) from salinity (S [PSU]),
!! potential temperature (T [degC]), and pressure [Pa].  It uses the expression from
!! Jackett et al., 2006, J. Atmos. Ocean. Tech., 32, 1709-1728.
subroutine calculate_density_array_Jackett(T, S, pres, rho, start, npts, rho_ref)
  real, dimension(:), intent(in)    :: T        !< Potential temperature relative to the surface [degC].
  real, dimension(:), intent(in)    :: S        !< Salinity [PSU].
  real, dimension(:), intent(in)    :: pres     !< Pressure [Pa].
  real, dimension(:), intent(inout) :: rho      !< In situ density [kg m-3].
  integer,            intent(in)    :: start    !< The starting point in the arrays.
  integer,            intent(in)    :: npts     !< The number of values to calculate.
  real,     optional, intent(in)    :: rho_ref  !< A reference density [kg m-3].

  ! Local variables
  real :: num_STP ! State dependent part of the numerator of the rational expresion
                  ! for density [kg m-3]
  real :: den     ! Denominator of the rational expresion for density [nondim]
  real :: den_STP ! State dependent part of the denominator of the rational expresion
                  ! for density [nondim]
  real :: I_den   ! The inverse of the denominator of the rational expresion for density [nondim]
  real :: T2      ! Temperature squared [degC2]
  real :: S1_2    ! Limited square root of salinity [PSU1/2]
  real :: rho0    ! The surface density of fresh water at 0 degC, perhaps less the refernce density [kg m-3]
  integer :: j

  do j=start,start+npts-1
    S1_2 = sqrt(max(0.0,s(j)))
    T2 = T(j)*T(j)

    num_STP = (T(j)*(RN010 + T(j)*(RN020 + T(j)*RN030)) + &
               S(j)*(RN100 + (T(j)*RN110 + S(j)*RN200)) ) + &
              pres(j)*(RN001 + ((T2*RN021 + S(j)*RN101) + pres(j)*(RN002 + T2*RN022)))
    den = 1.0 + ((T(j)*(RD010 + T(j)*(RD020 + T(j)*(RD030 + T(j)* RD040))) + &
                  S(j)*(RD100 + (T(j)*(RD110 + T2*RD130) + S1_2*(RD600 + T2*RD620))) ) + &
                 pres(j)*(RD001 + pres(j)*T(j)*(T2*RD032 + pres(j)*RD013)) )
    I_den = 1.0 / den

    rho0 = RN000
    if (present(rho_ref)) rho0 = RN000 - rho_ref*den

    rho(j) = (rho0 + num_STP)*I_den
  enddo

end subroutine calculate_density_array_Jackett

!> Computes the Jackett et al. in situ specific volume of sea water for 1-d array inputs and outputs.
!!
!! Returns the in situ specific volume of sea water (specvol in [m3 kg-1]) from salinity (S [PSU]),
!! potential temperature (T [degC]) and pressure [Pa].  It uses the expression from
!! Jackett et al., 2006, J. Atmos. Ocean. Tech., 32, 1709-1728.
!! If spv_ref is present, specvol is an anomaly from spv_ref.
subroutine calculate_spec_vol_array_Jackett(T, S, pres, specvol, start, npts, spv_ref)
  real, dimension(:), intent(in)    :: T        !< potential temperature relative to the
                                                !! surface [degC].
  real, dimension(:), intent(in)    :: S        !< salinity [PSU].
  real, dimension(:), intent(in)    :: pres     !< pressure [Pa].
  real, dimension(:), intent(inout) :: specvol  !< in situ specific volume [m3 kg-1].
  integer,            intent(in)    :: start    !< the starting point in the arrays.
  integer,            intent(in)    :: npts     !< the number of values to calculate.
  real,     optional, intent(in)    :: spv_ref  !< A reference specific volume [m3 kg-1].

  ! Local variables
  real :: num_STP ! State dependent part of the numerator of the rational expresion
                  ! for density (not specific volume) [kg m-3]
  real :: den_STP ! State dependent part of the denominator of the rational expresion
                  ! for density (not specific volume) [nondim]
  real :: I_num   ! The inverse of the numerator of the rational expresion for density [nondim]
  real :: T2      ! Temperature squared [degC2]
  real :: S1_2    ! Limited square root of salinity [PSU1/2]
  integer :: j

  do j=start,start+npts-1
    S1_2 = sqrt(max(0.0,s(j)))
    T2 = T(j)*T(j)

    num_STP = (T(j)*(RN010 + T(j)*(RN020 + T(j)*RN030)) + &
               S(j)*(RN100 + (T(j)*RN110 + S(j)*RN200)) ) + &
              pres(j)*(RN001 + ((T2*RN021 + S(j)*RN101) + pres(j)*(RN002 + T2*RN022)))
    den_STP = (T(j)*(RD010 + T(j)*(RD020 + T(j)*(RD030 + T(j)* RD040))) + &
               S(j)*(RD100 + (T(j)*(RD110 + T2*RD130) + S1_2*(RD600 + T2*RD620))) ) + &
              pres(j)*(RD001 + pres(j)*T(j)*(T2*RD032 + pres(j)*RD013))
    I_num = 1.0 / (RN000 + num_STP)
    if (present(spv_ref)) then
      ! This form is slightly more complicated, but it cancels the leading terms better.
      specvol(j) = ((1.0 - spv_ref*RN000) + (den_STP - spv_ref*num_STP)) * I_num
    else
      specvol(j) = (1.0 + den_STP) * I_num
    endif
  enddo

end subroutine calculate_spec_vol_array_Jackett

!> Return the thermal/haline expansion coefficients for 1-d array inputs and outputs
subroutine calculate_density_derivs_array_Jackett(T, S, pres, drho_dT, drho_dS, start, npts)
  real,    intent(in),    dimension(:) :: T        !< Potential temperature relative to the
                                                   !! surface [degC].
  real,    intent(in),    dimension(:) :: S        !< Salinity [PSU].
  real,    intent(in),    dimension(:) :: pres     !< pressure [Pa].
  real,    intent(inout), dimension(:) :: drho_dT  !< The partial derivative of density with potential
                                                   !! temperature [kg m-3 degC-1].
  real,    intent(inout), dimension(:) :: drho_dS  !< The partial derivative of density with salinity,
                                                   !! in [kg m-3 PSU-1].
  integer, intent(in)                  :: start    !< The starting point in the arrays.
  integer, intent(in)                  :: npts     !< The number of values to calculate.

  ! Local variables
  real :: num     ! Numerator of the rational expresion for density [kg m-3]
  real :: den     ! Denominator of the rational expresion for density [nondim]
  real :: I_denom2 ! The inverse of the square of the denominator of the rational expression
                  ! for density [nondim]
  real :: dnum_dT ! The derivative of num with potential temperature [kg m-3 degC-1]
  real :: dnum_dS ! The derivative of num with salinity [kg m-3 PSU-1]
  real :: dden_dT ! The derivative of den with potential temperature [degC-1]
  real :: dden_dS ! The derivative of den with salinity PSU-1]
  real :: T2      ! Temperature squared [degC2]
  real :: S1_2    ! Limited square root of salinity [PSU1/2]
  integer :: j

  do j=start,start+npts-1
    S1_2 = sqrt(max(0.0,s(j)))
    T2 = T(j)*T(j)

    num = RN000 + ((T(j)*(RN010 + T(j)*(RN020 + T(j)*RN030)) + &
                    S(j)*(RN100 + (T(j)*RN110 + S(j)*RN200)) ) + &
                   pres(j)*(RN001 + ((T2*RN021 + S(j)*RN101) + pres(j)*(RN002 + T2*RN022))) )
    den = 1.0 + ((T(j)*(RD010 + T(j)*(RD020 + T(j)*(RD030 + T(j)* RD040))) + &
                  S(j)*(RD100 + (T(j)*(RD110 + T2*RD130) + S1_2*(RD600 + T2*RD620))) ) + &
                 pres(j)*(RD001 + pres(j)*T(j)*(T2*RD032 + pres(j)*RD013)) )

    dnum_dT = ((RN010 + T(j)*(2.*RN020 + T(j)*(3.*RN030))) + S(j)*RN110) + &
              pres(j)*T(j)*(2.*RN021 + pres(j)*(2.*RN022))
    dnum_dS = (RN100 + (T(j)*RN110 + S(j)*(2.*RN200))) + pres(j)*RN101
    dden_dT = ((RD010 + T(j)*((2.*RD020) + T(j)*((3.*RD030) + T(j)*(4.*RD040)))) + &
               S(j)*((RD110 + T2*(3.*RD130)) + S1_2*T(j)*(2.*RD620)) ) + &
              pres(j)**2*(T2*3.*RD032 + pres(j)*RD013)
    dden_dS = RD100 + (T(j)*(RD110 + T2*RD130) + S1_2*(1.5*RD600 + T2*(1.5*RD620)))
    I_denom2 = 1.0 / den**2

    ! rho(j) = num / den
    drho_dT(j) = (dnum_dT * den - num * dden_dT) * I_denom2
    drho_dS(j) = (dnum_dS * den - num * dden_dS) * I_denom2
  enddo

end subroutine calculate_density_derivs_array_Jackett

!> Return the partial derivatives of specific volume with temperature and salinity
!! for 1-d array inputs and outputs
subroutine calculate_specvol_derivs_Jackett06(T, S, pres, dSV_dT, dSV_dS, start, npts)
  real,    intent(in),    dimension(:) :: T        !< Potential temperature relative to the surface [degC].
  real,    intent(in),    dimension(:) :: S        !< Salinity [PSU].
  real,    intent(in),    dimension(:) :: pres     !< Pressure [Pa].
  real,    intent(inout), dimension(:) :: dSV_dT   !< The partial derivative of specific volume with
                                                   !! potential temperature [m3 kg-1 degC-1].
  real,    intent(inout), dimension(:) :: dSV_dS   !< The partial derivative of specific volume with
                                                   !! salinity [m3 kg-1 PSU-1].
  integer, intent(in)                  :: start    !< The starting point in the arrays.
  integer, intent(in)                  :: npts     !< The number of values to calculate.

  ! Local variables
  real :: num     ! Numerator of the rational expresion for density (not specific volume) [kg m-3]
  real :: den     ! Denominator of the rational expresion for density (not specific volume) [nondim]
  real :: I_num2  ! The inverse of the square of the numerator of the rational expression
                  ! for density [nondim]
  real :: dnum_dT ! The derivative of num with potential temperature [kg m-3 degC-1]
  real :: dnum_dS ! The derivative of num with salinity [kg m-3 PSU-1]
  real :: dden_dT ! The derivative of den with potential temperature [degC-1]
  real :: dden_dS ! The derivative of den with salinity PSU-1]
  real :: T2      ! Temperature squared [degC2]
  real :: S1_2    ! Limited square root of salinity [PSU1/2]
  integer :: j

  do j=start,start+npts-1
    S1_2 = sqrt(max(0.0,s(j)))
    T2 = T(j)*T(j)

    num = RN000 + ((T(j)*(RN010 + T(j)*(RN020 + T(j)*RN030)) + &
                    S(j)*(RN100 + (T(j)*RN110 + S(j)*RN200)) ) + &
                   pres(j)*(RN001 + ((T2*RN021 + S(j)*RN101) + pres(j)*(RN002 + T2*RN022))) )
    den = 1.0 + ((T(j)*(RD010 + T(j)*(RD020 + T(j)*(RD030 + T(j)* RD040))) + &
                  S(j)*(RD100 + (T(j)*(RD110 + T2*RD130) + S1_2*(RD600 + T2*RD620))) ) + &
                 pres(j)*(RD001 + pres(j)*T(j)*(T2*RD032 + pres(j)*RD013)) )

    dnum_dT = ((RN010 + T(j)*(2.*RN020 + T(j)*(3.*RN030))) + S(j)*RN110) + &
              pres(j)*T(j)*(2.*RN021 + pres(j)*(2.*RN022))
    dnum_dS = (RN100 + (T(j)*RN110 + S(j)*(2.*RN200))) + pres(j)*RN101
    dden_dT = ((RD010 + T(j)*((2.*RD020) + T(j)*((3.*RD030) + T(j)*(4.*RD040)))) + &
               S(j)*((RD110 + T2*(3.*RD130)) + S1_2*T(j)*(2.*RD620)) ) + &
              pres(j)**2*(T2*3.*RD032 + pres(j)*RD013)
    dden_dS = RD100 + (T(j)*(RD110 + T2*RD130) + S1_2*(1.5*RD600 + T2*(1.5*RD620)))
    I_num2 = 1.0 / num**2

    ! SV(j) = den / num
    dSV_dT(j) = (num * dden_dT - dnum_dT * den) * I_num2
    dSV_dS(j) = (num * dden_dS - dnum_dS * den) * I_num2
  enddo

end subroutine calculate_specvol_derivs_Jackett06

!> Computes the compressibility of seawater for 1-d array inputs and outputs
subroutine calculate_compress_Jackett06(T, S, pres, rho, drho_dp, start, npts)
  real,    intent(in),    dimension(:) :: T        !< Potential temperature relative to the surface [degC].
  real,    intent(in),    dimension(:) :: S        !< Salinity [PSU].
  real,    intent(in),    dimension(:) :: pres     !< Pressure [Pa].
  real,    intent(inout), dimension(:) :: rho      !< In situ density [kg m-3].
  real,    intent(inout), dimension(:) :: drho_dp  !< The partial derivative of density with pressure
                                                   !! (also the inverse of the square of sound speed)
                                                   !! [s2 m-2].
  integer, intent(in)                  :: start    !< The starting point in the arrays.
  integer, intent(in)                  :: npts     !< The number of values to calculate.

  ! Local variables
  real :: num     ! Numerator of the rational expresion for density [kg m-3]
  real :: den     ! Denominator of the rational expresion for density [nondim]
  real :: I_den   ! The inverse of the denominator of the rational expression for density [nondim]
  real :: dnum_dp ! The derivative of num with pressure [kg m-3 dbar-1]
  real :: dden_dp ! The derivative of den with pressure [dbar-1]
  real :: T2      ! Temperature squared [degC2]
  real :: S1_2    ! Limited square root of salinity [PSU1/2]
  integer :: j

  do j=start,start+npts-1
    S1_2 = sqrt(max(0.0,s(j)))
    T2 = T(j)*T(j)

    num = RN000 + ((T(j)*(RN010 + T(j)*(RN020 + T(j)*RN030)) + &
                    S(j)*(RN100 + (T(j)*RN110 + S(j)*RN200)) ) + &
                   pres(j)*(RN001 + ((T2*RN021 + S(j)*RN101) + pres(j)*(RN002 + T2*RN022))) )
    den = 1.0 + ((T(j)*(RD010 + T(j)*(RD020 + T(j)*(RD030 + T(j)* RD040))) + &
                  S(j)*(RD100 + (T(j)*(RD110 + T2*RD130) + S1_2*(RD600 + T2*RD620))) ) + &
                 pres(j)*(RD001 + pres(j)*T(j)*(T2*RD032 + pres(j)*RD013)) )
    dnum_dp = RN001 + ((T2*RN021 + S(j)*RN101) + pres(j)*(2.*RN002 + T2*(2.*RN022)))
    dden_dp = RD001 + pres(j)*T(j)*(T2*(2.*RD032) + pres(j)*(3.*RD013))

    I_den  = 1.0 / den
    rho(j) = num * I_den
    drho_dp(j) = (dnum_dp * den - num * dden_dp) * I_den**2
  enddo
end subroutine calculate_compress_Jackett06

!> Second derivatives of density with respect to temperature, salinity, and pressure for 1-d array inputs and outputs.
subroutine calculate_density_second_derivs_array_Jackett(T, S, P, drho_ds_ds, drho_ds_dt, drho_dt_dt, &
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
  real :: num         ! Numerator of the rational expresion for density [kg m-3]
  real :: den         ! Denominator of the rational expresion for density [nondim]
  real :: I_num2      ! The inverse of the square of the numerator of the rational expression
                      ! for density [nondim]
  real :: dnum_dT     ! The derivative of num with potential temperature [kg m-3 degC-1]
  real :: dnum_dS     ! The derivative of num with salinity [kg m-3 PSU-1]
  real :: dden_dT     ! The derivative of den with potential temperature [degC-1]
  real :: dden_dS     ! The derivative of den with salinity PSU-1]
  real :: dnum_dp     ! The derivative of num with pressure [kg m-3 dbar-1]
  real :: dden_dp     ! The derivative of det with pressure [dbar-1]
  real :: d2num_dT2   ! The second derivative of num with potential temperature [kg m-3 degC-2]
  real :: d2num_dT_dS ! The second derivative of num with potential temperature and
                      ! salinity [kg m-3 degC-1 PSU-1]
  real :: d2num_dS2   ! The second derivative of num with salinity [kg m-3 PSU-2]
  real :: d2num_dT_dp ! The second derivative of num with potential temperature and
                      ! pressure [kg m-3 degC-1 dbar-1]
  real :: d2num_dS_dp ! The second derivative of num with salinity and
                      ! pressure [kg m-3 PSU-1 dbar-1]
  real :: d2den_dT2   ! The second derivative of den with potential temperature [degC-2]
  real :: d2den_dT_dS ! The second derivative of den with potential temperature and salinity [degC-1 PSU-1]
  real :: d2den_dS2   ! The second derivative of den with salinity [PSU-2]
  real :: d2den_dT_dp ! The second derivative of den with potential temperature and pressure [degC-1 dbar-1]
  real :: d2den_dS_dp ! The second derivative of den with salinity and pressure [PSU-1 dbar-1]
  real :: T2          ! Temperature squared [degC2]
  real :: S1_2        ! Limited square root of salinity [PSU1/2]
  real :: I_s12       ! The inverse of the square root of salinity [PSU-1/2]
  real :: I_denom2    ! The inverse of the square of the denominator of the rational expression
                      ! for density [nondim]
  real :: I_denom3    ! The inverse of the cube of the denominator of the rational expression
                      ! for density [nondim]
  integer :: j

  do j = start,start+npts-1
    S1_2 = sqrt(max(0.0,s(j)))
    T2 = T(j)*T(j)

    num = RN000 + ((T(j)*(RN010 + T(j)*(RN020 + T(j)*RN030)) + &
                    S(j)*(RN100 + (T(j)*RN110 + S(j)*RN200)) ) + &
                   P(j)*(RN001 + ((T2*RN021 + S(j)*RN101) + P(j)*(RN002 + T2*RN022))) )
    den = 1.0 + ((T(j)*(RD010 + T(j)*(RD020 + T(j)*(RD030 + T(j)* RD040))) + &
                  S(j)*(RD100 + (T(j)*(RD110 + T2*RD130) + S1_2*(RD600 + T2*RD620))) ) + &
                 P(j)*(RD001 + P(j)*T(j)*(T2*RD032 + P(j)*RD013)) )
    ! rho(j) = num*I_den

    dnum_dT = ((RN010 + T(j)*(2.*RN020 + T(j)*(3.*RN030))) + S(j)*RN110) + &
              P(j)*T(j)*(2.*RN021 + P(j)*(2.*RN022))
    dnum_dS = (RN100 + (T(j)*RN110 + S(j)*(2.*RN200))) + P(j)*RN101
    dnum_dp = RN001 + ((T2*RN021 + S(j)*RN101) + P(j)*(2.*RN002 + T2*(2.*RN022)))
    d2num_dT2 = 2.*RN020 + T(j)*(6.*RN030) + P(j)*(2.*RN021 + P(j)*(2.*RN022))
    d2num_dT_dS = RN110
    d2num_dS2 = 2.*RN200
    d2num_dT_dp = T(j)*(2.*RN021 + P(j)*(4.*RN022))
    d2num_dS_dp = RN101

    dden_dT = ((RD010 + T(j)*((2.*RD020) + T(j)*((3.*RD030) + T(j)*(4.*RD040)))) + &
               S(j)*((RD110 + T2*(3.*RD130)) + S1_2*T(j)*(2.*RD620)) ) + &
              P(j)**2*(T2*3.*RD032 + P(j)*RD013)
    dden_dS = RD100 + (T(j)*(RD110 + T2*RD130) + S1_2*(1.5*RD600 + T2*(1.5*RD620)))
    dden_dp = RD001 + P(j)*T(j)*(T2*(2.*RD032) + P(j)*(3.*RD013))

    d2den_dT2 = (((2.*RD020) + T(j)*((6.*RD030) + T(j)*(12.*RD040))) + &
                 S(j)*(T(j)*(6.*RD130) + S1_2*(2.*RD620)) ) + P(j)**2*(T(j)*(6.*RD032))
    d2den_dT_dS = (RD110 + T2*3.*RD130) + (T(j)*S1_2)*(3.0*RD620)
    d2den_dT_dp = P(j)*(T2*(6.*RD032) + P(j)*(3.*RD013))
    d2den_dS_dp = 0.0

    ! The Jackett et al. 2006 equation of state is a fit to density, but it chooses a form that
    ! exhibits a singularity in the second derivatives with salinity for fresh water.  To avoid
    ! this, the square root of salinity can be treated with a floor such that the contribution from
    ! the S**1.5 terms to both the surface density and the secant bulk modulus are lost to roundoff.
    ! This salinity is given by (~1e-16/RD600)**(2/3) ~= 7e-8 PSU, or S1_2 ~= 2.6e-4
    I_S12 = 1.0 / (max(S1_2, 1.0e-4))
    d2den_dS2 = (0.75*RD600 + T2*(0.75*RD620)) * I_S12

    I_denom3 = 1.0 / den**3

    ! In deriving the following, it is useful to note that:
    !   drho_dp(j) = (dnum_dp * den - num * dden_dp) / den**2
    !   drho_dT(j) = (dnum_dT * den - num * dden_dT) / den**2
    !   drho_dS(j) = (dnum_dS * den - num * dden_dS) / den**2
    drho_dS_dS(j) = (den*(den*d2num_dS2 - 2.*dnum_dS*dden_dS) + num*(2.*dden_dS**2 - den*d2den_dS2)) * I_denom3
    drho_dS_dt(j) = (den*(den*d2num_dT_dS - (dnum_dT*dden_dS + dnum_dS*dden_dT)) + &
                     num*(2.*dden_dT*dden_dS - den*d2den_dT_dS)) * I_denom3
    drho_dT_dT(j) = (den*(den*d2num_dT2 - 2.*dnum_dT*dden_dT) + num*(2.*dden_dT**2 - den*d2den_dT2)) * I_denom3

    drho_dS_dp(j) = (den*(den*d2num_dS_dp - (dnum_dp*dden_dS + dnum_dS*dden_dp)) + &
                     num*(2.*dden_dS*dden_dp - den*d2den_dS_dp)) * I_denom3
    drho_dT_dp(j) = (den*(den*d2num_dT_dp - (dnum_dp*dden_dT + dnum_dT*dden_dp)) + &
                     num*(2.*dden_dT*dden_dp - den*d2den_dT_dp)) * I_denom3
  enddo

end subroutine calculate_density_second_derivs_array_Jackett

!> Computes the in situ density of sea water for scalar inputs and outputs.
!!
!! Returns the in situ density of sea water (rho in [kg m-3]) from salinity (S [PSU]),
!! potential temperature (T [degC]), and pressure [Pa].  It uses the expression from
!! Jackett et al., 2006, J. Atmos. Ocean. Tech., 32, 1709-1728.
subroutine calculate_density_scalar_Jackett(T, S, pressure, rho, rho_ref)
  real,           intent(in)  :: T        !< Potential temperature relative to the surface [degC].
  real,           intent(in)  :: S        !< Salinity [PSU].
  real,           intent(in)  :: pressure !< pressure [Pa].
  real,           intent(out) :: rho      !< In situ density [kg m-3].
  real, optional, intent(in)  :: rho_ref  !< A reference density [kg m-3].

  ! Local variables
  real, dimension(1) :: T0    ! A 1-d array with a copy of the potential temperature [degC]
  real, dimension(1) :: S0    ! A 1-d array with a copy of the salinity [PSU]
  real, dimension(1) :: pressure0 ! A 1-d array with a copy of the pressure [Pa]
  real, dimension(1) :: rho0  ! A 1-d array with a copy of the density [kg m-3]

  T0(1) = T ; S0(1) = S ; pressure0(1) = pressure
  call calculate_density_array_Jackett(T0, S0, pressure0, rho0, 1, 1, rho_ref)
  rho = rho0(1)

end subroutine calculate_density_scalar_Jackett

!> Computes the Jackett et al. 2006 in situ specific volume of sea water for scalar inputs and outputs.
!!
!! Returns the in situ specific volume of sea water (specvol in [m3 kg-1]) from salinity (S [PSU]),
!! potential temperature (T [degC]) and pressure [Pa].  It uses the expression from
!! Jackett et al., 2006, J. Atmos. Ocean. Tech., 32, 1709-1728.
!! If spv_ref is present, specvol is an anomaly from spv_ref.
subroutine calculate_spec_vol_scalar_Jackett(T, S, pressure, specvol, spv_ref)
  real,           intent(in)  :: T        !< potential temperature relative to the surface [degC].
  real,           intent(in)  :: S        !< salinity [PSU].
  real,           intent(in)  :: pressure !< pressure [Pa].
  real,           intent(out) :: specvol  !< in situ specific volume [m3 kg-1].
  real, optional, intent(in)  :: spv_ref  !< A reference specific volume [m3 kg-1].

  ! Local variables
  real, dimension(1) :: T0    ! A 1-d array with a copy of the potential temperature [degC]
  real, dimension(1) :: S0    ! A 1-d array with a copy of the salinity [PSU]
  real, dimension(1) :: pressure0 ! A 1-d array with a copy of the pressure [Pa]
  real, dimension(1) :: spv0  ! A 1-d array with a copy of the specific volume [m3 kg-1]

  T0(1) = T ; S0(1) = S ; pressure0(1) = pressure
  call calculate_spec_vol_array_Jackett(T0, S0, pressure0, spv0, 1, 1, spv_ref)
  specvol = spv0(1)
end subroutine calculate_spec_vol_scalar_Jackett

!> Return the thermal/haline expansion coefficients for scalar inputs and outputs
!!
!! The scalar version of calculate_density_derivs promotes scalar inputs to 1-element array
!! and then demotes the output back to a scalar
subroutine calculate_density_derivs_scalar_Jackett(T, S, pressure, drho_dT, drho_dS)
  real,    intent(in)  :: T        !< Potential temperature relative to the surface [degC].
  real,    intent(in)  :: S        !< Salinity [PSU].
  real,    intent(in)  :: pressure !< pressure [Pa].
  real,    intent(out) :: drho_dT  !< The partial derivative of density with potential
                                   !! temperature [kg m-3 degC-1].
  real,    intent(out) :: drho_dS  !< The partial derivative of density with salinity,
                                   !! in [kg m-3 PSU-1].

  ! Local variables needed to promote the input/output scalars to 1-element arrays
  real, dimension(1) :: T0    ! A 1-d array with a copy of the temperature [degC]
  real, dimension(1) :: S0    ! A 1-d array with a copy of the salinity [PSU]
  real, dimension(1) :: p0    ! A 1-d array with a copy of the pressure [Pa]
  real, dimension(1) :: drdt0 ! The derivative of density with temperature [kg m-3 degC-1]
  real, dimension(1) :: drds0 ! The derivative of density with salinity [kg m-3 PSU-1]

  T0(1) = T ; S0(1) = S ; P0(1) = pressure
  call calculate_density_derivs_array_Jackett(T0, S0, P0, drdt0, drds0, 1, 1)
  drho_dT = drdt0(1) ; drho_dS = drds0(1)

end subroutine calculate_density_derivs_scalar_Jackett

!> Second derivatives of density with respect to temperature, salinity, and pressure for scalar inputs.
!!
!! The scalar version of calculate_density_second_derivs promotes scalar inputs to 1-element array
!! and then demotes the output back to a scalar
subroutine calculate_density_second_derivs_scalar_Jackett(T, S, P, drho_ds_ds, drho_ds_dt, drho_dt_dt, &
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

  T0(1) = T ; S0(1) = S ; P0(1) = P
  call calculate_density_second_derivs_array_Jackett(T0, S0, P0, drdsds, drdsdt, drdtdt, drdsdp, drdtdp, 1, 1)
  drho_ds_ds = drdsds(1) ; drho_ds_dt = drdsdt(1) ; drho_dt_dt = drdtdt(1)
  drho_ds_dp = drdsdp(1) ; drho_dt_dp = drdtdp(1)

end subroutine calculate_density_second_derivs_scalar_Jackett

!> Return the range of temperatures, salinities and pressures for which the Jackett et al. (2006)
!! equation of state has been fitted to observations.  Care should be taken when applying this
!! equation of state outside of its fit range.
subroutine EoS_fit_range_Jackett06(T_min, T_max, S_min, S_max, p_min, p_max)
  real, optional, intent(out) :: T_min !< The minimum potential temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: T_max !< The maximum potential temperature over which this EoS is fitted [degC]
  real, optional, intent(out) :: S_min !< The minimum practical salinity over which this EoS is fitted [PSU]
  real, optional, intent(out) :: S_max !< The maximum practical salinity over which this EoS is fitted [PSU]
  real, optional, intent(out) :: p_min !< The minimum pressure over which this EoS is fitted [Pa]
  real, optional, intent(out) :: p_max !< The maximum pressure over which this EoS is fitted [Pa]

  ! Note that the actual fit range is given for the surface range of temperatures and salinities,
  ! but Jackett et al. use a more limited range of properties at higher pressures.
  if (present(T_min)) T_min = -4.5
  if (present(T_max)) T_max = 40.0
  if (present(S_min)) S_min =  0.0
  if (present(S_max)) S_max = 42.0
  if (present(p_min)) p_min = 0.0
  if (present(p_max)) p_max = 8.5e7

end subroutine EoS_fit_range_Jackett06

!> \namespace mom_eos_Jackett06
!!
!! \section section_EOS_Jackett06 Jackett et al. 2006 (Hycom-25-term) equation of state
!!
!! Jackett et al. (2006) provide an approximation for the in situ density as a function of
!! potential temperature, salinity, and pressure.  This 25 term equation of state is
!! frequently used in Hycom for a potential density, at which point it only has 17 terms
!! and so is commonly called the "17-term equation of state" there.  Here the full expressions
!! for the in situ densities are used.
!!
!! The functional form of this equation of state includes terms proportional to salinity to the
!! 3/2 power.  This introduces a singularity in the second derivative of density with salinity
!! at a salinity of 0, but this has been addressed here by setting a floor of 1e-8 PSU on the
!! salinity that is used in the denominator of these second derivative expressions.  This value
!! was chosen to imply a contribution that is smaller than numerical roundoff in the expression for
!! density, which is the field for which the Jackett et al. equation of state was originally derived.
!!
!! \subsection section_EOS_Jackett06_references References
!!
!! Jackett, D., T. McDougall, R. Feistel, D. Wright and S. Griffies (2006),
!!   Algorithms for density, potential temperature, conservative
!!   temperature, and the freezing temperature of seawater, JAOT
!!   doi.org/10.1175/JTECH1946.1

end module MOM_EOS_Jackett06
