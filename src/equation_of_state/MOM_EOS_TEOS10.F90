!> The equation of state using the TEOS10 expressions
module MOM_EOS_TEOS10

! This file is part of MOM6. See LICENSE.md for the license.

!***********************************************************************
!*  The subroutines in this file implement the equation of state for   *
!*  sea water using the TEOS10 functions                               *
!***********************************************************************

use gsw_mod_toolbox, only : gsw_sp_from_sr, gsw_pt_from_ct
use gsw_mod_toolbox, only : gsw_rho, gsw_specvol
use gsw_mod_toolbox, only : gsw_rho_first_derivatives, gsw_specvol_first_derivatives
use gsw_mod_toolbox, only : gsw_rho_second_derivatives
!use gsw_mod_toolbox, only : gsw_sr_from_sp, gsw_ct_from_pt

implicit none ; private

public calculate_compress_teos10, calculate_density_teos10, calculate_spec_vol_teos10
public calculate_density_derivs_teos10, calculate_specvol_derivs_teos10
public calculate_density_second_derivs_teos10, EoS_fit_range_teos10
public gsw_sp_from_sr, gsw_pt_from_ct

!> Compute the in situ density of sea water ([kg m-3]), or its anomaly with respect to
!! a reference density, from absolute salinity (g/kg), conservative temperature (in deg C),
!! and pressure [Pa], using the TEOS10 expressions.
interface calculate_density_teos10
  module procedure calculate_density_scalar_teos10, calculate_density_array_teos10
end interface calculate_density_teos10

!> Compute the in situ specific volume of sea water (in [m3 kg-1]), or an anomaly with respect
!! to a reference specific volume, from absolute salinity (in g/kg), conservative temperature
!! (in deg C), and pressure [Pa], using the TEOS10 expressions.
interface calculate_spec_vol_teos10
  module procedure calculate_spec_vol_scalar_teos10, calculate_spec_vol_array_teos10
end interface calculate_spec_vol_teos10

!> For a given thermodynamic state, return the derivatives of density with conservative temperature
!! and absolute salinity, using the TEOS10 expressions.
interface calculate_density_derivs_teos10
  module procedure calculate_density_derivs_scalar_teos10, calculate_density_derivs_array_teos10
end interface calculate_density_derivs_teos10

!> For a given thermodynamic state, return the second derivatives of density with various combinations
!! of conservative temperature, absolute salinity, and pressure, using the TEOS10 expressions.
interface calculate_density_second_derivs_teos10
  module procedure calculate_density_second_derivs_scalar_teos10, calculate_density_second_derivs_array_teos10
end interface calculate_density_second_derivs_teos10

real, parameter :: Pa2db  = 1.e-4  !< The conversion factor from Pa to dbar [dbar Pa-1]

contains

!> This subroutine computes the in situ density of sea water (rho in [kg m-3])
!! from absolute salinity (S [g kg-1]), conservative temperature (T [degC]),
!! and pressure [Pa].  It uses the expression from the TEOS10 website.
subroutine calculate_density_scalar_teos10(T, S, pressure, rho, rho_ref)
  real,           intent(in)  :: T        !< Conservative temperature [degC].
  real,           intent(in)  :: S        !< Absolute salinity [g kg-1].
  real,           intent(in)  :: pressure !< pressure [Pa].
  real,           intent(out) :: rho      !< In situ density [kg m-3].
  real, optional, intent(in)  :: rho_ref  !< A reference density [kg m-3].

  ! Local variables
  real, dimension(1) :: T0    ! A 1-d array with a copy of the conservative temperature [degC]
  real, dimension(1) :: S0    ! A 1-d array with a copy of the absolute salinity [g kg-1]
  real, dimension(1) :: pressure0 ! A 1-d array with a copy of the pressure [Pa]
  real, dimension(1) :: rho0  ! A 1-d array with a copy of the density [kg m-3]

  T0(1) = T
  S0(1) = S
  pressure0(1) = pressure

  call calculate_density_array_teos10(T0, S0, pressure0, rho0, 1, 1, rho_ref)
  rho = rho0(1)

end subroutine calculate_density_scalar_teos10

!> This subroutine computes the in situ density of sea water (rho in [kg m-3])
!! from absolute salinity (S [g kg-1]), conservative temperature (T [degC]),
!! and pressure [Pa].  It uses the expression from the
!! TEOS10 website.
subroutine calculate_density_array_teos10(T, S, pressure, rho, start, npts, rho_ref)
  real, dimension(:), intent(in)  :: T        !< Conservative temperature [degC].
  real, dimension(:), intent(in)  :: S        !< Absolute salinity [g kg-1]
  real, dimension(:), intent(in)  :: pressure !< pressure [Pa].
  real, dimension(:), intent(out) :: rho      !< in situ density [kg m-3].
  integer,            intent(in)  :: start    !< the starting point in the arrays.
  integer,            intent(in)  :: npts     !< the number of values to calculate.
  real,     optional, intent(in)  :: rho_ref  !< A reference density [kg m-3].

  ! Local variables
  real :: zs  ! Absolute salinity [g kg-1]
  real :: zt  ! Conservative temperature [degC]
  real :: zp  ! Pressure converted to decibars [dbar]
  integer :: j

  do j=start,start+npts-1
    !Conversions
    zs = S(j) !gsw_sr_from_sp(S(j))       !Convert practical salinity to absolute salinity
    zt = T(j) !gsw_ct_from_pt(S(j),T(j))  !Convert potential temp to conservative temp
    zp = pressure(j)* Pa2db         !Convert pressure from Pascal to decibar

    if (S(j) < -1.0e-10) then !Can we assume safely that this is a missing value?
      rho(j) = 1000.0
    else
      rho(j) = gsw_rho(zs,zt,zp)
    endif
    if (present(rho_ref)) rho(j) = rho(j) - rho_ref
  enddo
end subroutine calculate_density_array_teos10

!> This subroutine computes the in situ specific volume of sea water (specvol in
!! [m3 kg-1]) from absolute salinity (S [g kg-1]), conservative temperature (T [degC])
!! and pressure [Pa], using the TEOS10 equation of state.
!! If spv_ref is present, specvol is an anomaly from spv_ref.
subroutine calculate_spec_vol_scalar_teos10(T, S, pressure, specvol, spv_ref)
  real,           intent(in)  :: T        !< Conservative temperature [degC].
  real,           intent(in)  :: S        !< Absolute salinity [g kg-1]
  real,           intent(in)  :: pressure !< pressure [Pa].
  real,           intent(out) :: specvol  !< in situ specific volume [m3 kg-1].
  real, optional, intent(in)  :: spv_ref  !< A reference specific volume [m3 kg-1].

  ! Local variables
  real, dimension(1) :: T0    ! A 1-d array with a copy of the conservative temperature [degC]
  real, dimension(1) :: S0    ! A 1-d array with a copy of the absolute salinity [g kg-1]
  real, dimension(1) :: pressure0 ! A 1-d array with a copy of the pressure [Pa]
  real, dimension(1) :: spv0  ! A 1-d array with a copy of the specific volume [m3 kg-1]

  T0(1) = T ; S0(1) = S ; pressure0(1) = pressure

  call calculate_spec_vol_array_teos10(T0, S0, pressure0, spv0, 1, 1, spv_ref)
  specvol = spv0(1)
end subroutine calculate_spec_vol_scalar_teos10


!> This subroutine computes the in situ specific volume of sea water (specvol in
!! [m3 kg-1]) from absolute salinity (S [g kg-1]), conservative temperature (T [degC])
!! and pressure [Pa], using the TEOS10 equation of state.
!! If spv_ref is present, specvol is an anomaly from spv_ref.
subroutine calculate_spec_vol_array_teos10(T, S, pressure, specvol, start, npts, spv_ref)
  real, dimension(:), intent(in)  :: T        !< Conservative temperature [degC].
  real, dimension(:), intent(in)  :: S        !< salinity [g kg-1].
  real, dimension(:), intent(in)  :: pressure !< pressure [Pa].
  real, dimension(:), intent(out) :: specvol  !< in situ specific volume [m3 kg-1].
  integer,            intent(in)  :: start    !< the starting point in the arrays.
  integer,            intent(in)  :: npts     !< the number of values to calculate.
  real,     optional, intent(in)  :: spv_ref  !< A reference specific volume [m3 kg-1].

  ! Local variables
  real :: zs  ! Absolute salinity [g kg-1]
  real :: zt  ! Conservative temperature [degC]
  real :: zp  ! Pressure converted to decibars [dbar]
  integer :: j

  do j=start,start+npts-1
    !Conversions
    zs = S(j) !gsw_sr_from_sp(S(j))       !Convert practical salinity to absolute salinity
    zt = T(j) !gsw_ct_from_pt(S(j),T(j))  !Convert potential temp to conservative temp
    zp = pressure(j)* Pa2db         !Convert pressure from Pascal to decibar

    if (S(j) < -1.0e-10) then
      specvol(j) = 0.001 !Can we assume safely that this is a missing value?
    else
      specvol(j) = gsw_specvol(zs,zt,zp)
    endif
    if (present(spv_ref)) specvol(j) = specvol(j) - spv_ref
  enddo

end subroutine calculate_spec_vol_array_teos10

!> For a given thermodynamic state, calculate the derivatives of density with conservative
!! temperature and absolute salinity, using the TEOS10 expressions.
subroutine calculate_density_derivs_array_teos10(T, S, pressure, drho_dT, drho_dS, start, npts)
  real,    intent(in),  dimension(:) :: T        !< Conservative temperature [degC].
  real,    intent(in),  dimension(:) :: S        !< Absolute salinity [g kg-1].
  real,    intent(in),  dimension(:) :: pressure !< pressure [Pa].
  real,    intent(out), dimension(:) :: drho_dT  !< The partial derivative of density with conservative
                                                 !! temperature [kg m-3 degC-1].
  real,    intent(out), dimension(:) :: drho_dS  !< The partial derivative of density with absolute salinity,
                                                 !! [kg m-3 (g/kg)-1].
  integer, intent(in)                :: start    !< The starting point in the arrays.
  integer, intent(in)                :: npts     !< The number of values to calculate.

  ! Local variables
  real :: zs  ! Absolute salinity [g kg-1]
  real :: zt  ! Conservative temperature [degC]
  real :: zp  ! Pressure converted to decibars [dbar]
  integer :: j

  do j=start,start+npts-1
    !Conversions
    zs = S(j) !gsw_sr_from_sp(S(j))       !Convert practical salinity to absolute salinity
    zt = T(j) !gsw_ct_from_pt(S(j),T(j))  !Convert potential temp to conservative temp
    zp = pressure(j)* Pa2db         !Convert pressure from Pascal to decibar
    if (S(j) < -1.0e-10) then   !Can we assume safely that this is a missing value?
      drho_dT(j) = 0.0 ; drho_dS(j) = 0.0
    else
      call gsw_rho_first_derivatives(zs, zt, zp, drho_dsa=drho_dS(j), drho_dct=drho_dT(j))
    endif
  enddo

end subroutine calculate_density_derivs_array_teos10

!> For a given thermodynamic state, calculate the derivatives of density with conservative
!! temperature and absolute salinity, using the TEOS10 expressions.
subroutine calculate_density_derivs_scalar_teos10(T, S, pressure, drho_dT, drho_dS)
  real,    intent(in)  :: T        !< Conservative temperature [degC]
  real,    intent(in)  :: S        !< Absolute Salinity [g kg-1]
  real,    intent(in)  :: pressure !< pressure [Pa].
  real,    intent(out) :: drho_dT  !< The partial derivative of density with conservative
                                   !! temperature [kg m-3 degC-1].
  real,    intent(out) :: drho_dS  !< The partial derivative of density with absolute salinity,
                                   !! [kg m-3 (g/kg)-1].

  ! Local variables
  real :: zs  ! Absolute salinity [g kg-1]
  real :: zt  ! Conservative temperature [degC]
  real :: zp  ! Pressure converted to decibars [dbar]

  !Conversions
  zs = S !gsw_sr_from_sp(S)       !Convert practical salinity to absolute salinity
  zt = T !gsw_ct_from_pt(S,T)  !Convert potential temp to conservative temp
  zp = pressure* Pa2db         !Convert pressure from Pascal to decibar
  if (S < -1.0e-10) return !Can we assume safely that this is a missing value?
  call gsw_rho_first_derivatives(zs, zt, zp, drho_dsa=drho_dS, drho_dct=drho_dT)
end subroutine calculate_density_derivs_scalar_teos10

!> For a given thermodynamic state, calculate the derivatives of specific volume with conservative
!! temperature and absolute salinity, using the TEOS10 expressions.
subroutine calculate_specvol_derivs_teos10(T, S, pressure, dSV_dT, dSV_dS, start, npts)
  real,    intent(in),  dimension(:) :: T        !< Conservative temperature [degC].
  real,    intent(in),  dimension(:) :: S        !< Absolute salinity [g kg-1].
  real,    intent(in),  dimension(:) :: pressure !< pressure [Pa].
  real,    intent(out), dimension(:) :: dSV_dT   !< The partial derivative of specific volume with
                                                 !! conservative temperature [m3 kg-1 degC-1].
  real,    intent(out), dimension(:) :: dSV_dS   !< The partial derivative of specific volume with
                                                 !! absolute salinity [m3 kg-1 (g/kg)-1].
  integer, intent(in)                :: start    !< The starting point in the arrays.
  integer, intent(in)                :: npts     !< The number of values to calculate.

  ! Local variables
  real :: zs  ! Absolute salinity [g kg-1]
  real :: zt  ! Conservative temperature [degC]
  real :: zp  ! Pressure converted to decibars [dbar]
  integer :: j

  do j=start,start+npts-1
    !Conversions
    zs = S(j) !gsw_sr_from_sp(S(j))       !Convert practical salinity to absolute salinity
    zt = T(j) !gsw_ct_from_pt(S(j),T(j))  !Convert potential temp to conservative temp
    zp = pressure(j)* Pa2db         !Convert pressure from Pascal to decibar
    if (S(j) < -1.0e-10) then   !Can we assume safely that this is a missing value?
      dSV_dT(j) = 0.0 ; dSV_dS(j) = 0.0
    else
      call gsw_specvol_first_derivatives(zs,zt,zp, v_sa=dSV_dS(j), v_ct=dSV_dT(j))
    endif
  enddo

end subroutine calculate_specvol_derivs_teos10

!> Calculate the 5 second derivatives of the equation of state for scalar inputs
subroutine calculate_density_second_derivs_scalar_teos10(T, S, pressure, drho_dS_dS, drho_dS_dT, &
                                                         drho_dT_dT, drho_dS_dP, drho_dT_dP)
  real, intent(in)     :: T          !< Conservative temperature [degC]
  real, intent(in)     :: S          !< Absolute Salinity [g kg-1]
  real, intent(in)     :: pressure   !< pressure [Pa].
  real, intent(out)    :: drho_dS_dS !< Partial derivative of beta with respect
                                     !! to S [kg m-3 (g/kg)-2]
  real, intent(out)    :: drho_dS_dT !< Partial derivative of beta with respect
                                     !! to T [kg m-3 (g/kg)-1 degC-1]
  real, intent(out)    :: drho_dT_dT !< Partial derivative of alpha with respect
                                     !! to T [kg m-3 degC-2]
  real, intent(out)    :: drho_dS_dP !< Partial derivative of beta with respect
                                     !! to pressure [kg m-3 (g/kg)-1 Pa-1] = [s2 m-2 (g/kg)-1]
  real, intent(out)    :: drho_dT_dP !< Partial derivative of alpha with respect
                                     !! to pressure [kg m-3 degC-1 Pa-1] = [s2 m-2 degC-1]

  ! Local variables
  real :: zs  ! Absolute salinity [g kg-1]
  real :: zt  ! Conservative temperature [degC]
  real :: zp  ! Pressure converted to decibars [dbar]

  !Conversions
  zs = S !gsw_sr_from_sp(S)       !Convert practical salinity to absolute salinity
  zt = T !gsw_ct_from_pt(S,T)  !Convert potential temp to conservative temp
  zp = pressure* Pa2db         !Convert pressure from Pascal to decibar
  if (S < -1.0e-10) return !Can we assume safely that this is a missing value?
  call gsw_rho_second_derivatives(zs, zt, zp, rho_sa_sa=drho_dS_dS, rho_sa_ct=drho_dS_dT, &
                                     rho_ct_ct=drho_dT_dT, rho_sa_p=drho_dS_dP, rho_ct_p=drho_dT_dP)

end subroutine calculate_density_second_derivs_scalar_teos10

!> Calculate the 5 second derivatives of the equation of state for scalar inputs
subroutine calculate_density_second_derivs_array_teos10(T, S, pressure, drho_dS_dS, drho_dS_dT, &
                                                        drho_dT_dT, drho_dS_dP, drho_dT_dP, start, npts)
  real, dimension(:), intent(in)     :: T          !< Conservative temperature [degC]
  real, dimension(:), intent(in)     :: S          !< Absolute Salinity [g kg-1]
  real, dimension(:), intent(in)     :: pressure   !< pressure [Pa].
  real, dimension(:), intent(out)    :: drho_dS_dS !< Partial derivative of beta with respect
                                                   !! to S [kg m-3 (g/kg)-2]
  real, dimension(:), intent(out)    :: drho_dS_dT !< Partial derivative of beta with respect
                                                   !! to T [kg m-3 (g/kg)-1 degC-1]
  real, dimension(:), intent(out)    :: drho_dT_dT !< Partial derivative of alpha with respect
                                                   !! to T [kg m-3 degC-2]
  real, dimension(:), intent(out)    :: drho_dS_dP !< Partial derivative of beta with respect
                                                   !! to pressure [kg m-3 (g/kg)-1 Pa-1] = [s2 m-2 (g/kg)-1]
  real, dimension(:), intent(out)    :: drho_dT_dP !< Partial derivative of alpha with respect
                                                   !! to pressure [kg m-3 degC-1 Pa-1] = [s2 m-2 degC-1]
  integer, intent(in)  :: start    !< The starting point in the arrays.
  integer, intent(in)  :: npts     !< The number of values to calculate.

  ! Local variables
  real :: zs  ! Absolute salinity [g kg-1]
  real :: zt  ! Conservative temperature [degC]
  real :: zp  ! Pressure converted to decibars [dbar]
  integer :: j

  do j=start,start+npts-1
    !Conversions
    zs = S(j) !gsw_sr_from_sp(S)       !Convert practical salinity to absolute salinity
    zt = T(j) !gsw_ct_from_pt(S,T)  !Convert potential temp to conservative temp
    zp = pressure(j)* Pa2db         !Convert pressure from Pascal to decibar
    if (S(j) < -1.0e-10) then   !Can we assume safely that this is a missing value?
      drho_dS_dS(j) = 0.0 ; drho_dS_dT(j) = 0.0 ; drho_dT_dT(j) = 0.0
      drho_dS_dP(j) = 0.0 ; drho_dT_dP(j) = 0.0
    else
      call gsw_rho_second_derivatives(zs, zt, zp, rho_sa_sa=drho_dS_dS(j), rho_sa_ct=drho_dS_dT(j), &
                                      rho_ct_ct=drho_dT_dT(j), rho_sa_p=drho_dS_dP(j), rho_ct_p=drho_dT_dP(j))
    endif
  enddo

end subroutine calculate_density_second_derivs_array_teos10

!> This subroutine computes the in situ density of sea water (rho in
!! [kg m-3]) and the compressibility (drho/dp = C_sound^-2)
!! (drho_dp [s2 m-2]) from absolute salinity (sal [g kg-1]),
!! conservative temperature (T [degC]), and pressure [Pa].  It uses the
!! subroutines from TEOS10 website
subroutine calculate_compress_teos10(T, S, pressure, rho, drho_dp, start, npts)
  real,    intent(in),  dimension(:) :: T        !< Conservative temperature [degC].
  real,    intent(in),  dimension(:) :: S        !< Absolute salinity [g kg-1].
  real,    intent(in),  dimension(:) :: pressure !< Pressure [Pa].
  real,    intent(out), dimension(:) :: rho      !< In situ density [kg m-3].
  real,    intent(out), dimension(:) :: drho_dp  !< The partial derivative of density with pressure
                                                 !! (also the inverse of the square of sound speed)
                                                 !! [s2 m-2].
  integer, intent(in)                :: start    !< The starting point in the arrays.
  integer, intent(in)                :: npts     !< The number of values to calculate.

  ! Local variables
  real :: zs  ! Absolute salinity [g kg-1]
  real :: zt  ! Conservative temperature [degC]
  real :: zp  ! Pressure converted to decibars [dbar]
  integer :: j

  do j=start,start+npts-1
    !Conversions
    zs = S(j) !gsw_sr_from_sp(S(j))       !Convert practical salinity to absolute salinity
    zt = T(j) !gsw_ct_from_pt(S(j),T(j))  !Convert potential temp to conservative temp
    zp = pressure(j)* Pa2db         !Convert pressure from Pascal to decibar
    if (S(j) < -1.0e-10) then   !Can we assume safely that this is a missing value?
      rho(j) = 1000.0 ; drho_dp(j) = 0.0
    else
      rho(j) = gsw_rho(zs,zt,zp)
      call gsw_rho_first_derivatives(zs,zt,zp, drho_dp=drho_dp(j))
    endif
  enddo
end subroutine calculate_compress_teos10


!> Return the range of temperatures, salinities and pressures for which the TEOS-10
!! equation of state has been fitted to observations.  Care should be taken when
!! applying this equation of state outside of its fit range.
subroutine EoS_fit_range_teos10(T_min, T_max, S_min, S_max, p_min, p_max)
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

end subroutine EoS_fit_range_teos10

end module MOM_EOS_TEOS10
