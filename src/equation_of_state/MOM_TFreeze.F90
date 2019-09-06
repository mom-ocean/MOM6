!> Freezing point expressions
module MOM_TFreeze

! This file is part of MOM6. See LICENSE.md for the license.

!********+*********+*********+*********+*********+*********+*********+**
!*  The subroutines in this file determine the potential temperature   *
!* at which sea-water freezes.                                         *
!********+*********+*********+*********+*********+*********+*********+**
use gsw_mod_toolbox, only : gsw_ct_freezing_exact

implicit none ; private

public calculate_TFreeze_linear, calculate_TFreeze_Millero, calculate_TFreeze_teos10

!> Compute the freezing point potential temperature [degC] from salinity [ppt] and
!! pressure [Pa] using a simple linear expression, with coefficients passed in as arguments.
interface calculate_TFreeze_linear
  module procedure calculate_TFreeze_linear_scalar, calculate_TFreeze_linear_array
end interface calculate_TFreeze_linear

!> Compute the freezing point potential temperature [degC] from salinity [PSU] and
!! pressure [Pa] using the expression from Millero (1978) (and in appendix A of Gill 1982),
!! but with the of the pressure dependence changed from 7.53e-8 to 7.75e-8 to make this an
!! expression for potential temperature (not in situ temperature), using a
!! value that is correct at the freezing point at 35 PSU and 5e6 Pa (500 dbar).
interface calculate_TFreeze_Millero
  module procedure calculate_TFreeze_Millero_scalar, calculate_TFreeze_Millero_array
end interface calculate_TFreeze_Millero

!> Compute the freezing point conservative temperature [degC] from absolute salinity [g/kg]
!! and pressure [Pa] using the TEOS10 package.
interface calculate_TFreeze_teos10
  module procedure calculate_TFreeze_teos10_scalar, calculate_TFreeze_teos10_array
end interface calculate_TFreeze_teos10

contains

!>  This subroutine computes the freezing point potential temperature
!!  [degC] from salinity [ppt], and pressure [Pa] using a simple
!!  linear expression, with coefficients passed in as arguments.
subroutine calculate_TFreeze_linear_scalar(S, pres, T_Fr, TFr_S0_P0, &
                                           dTFr_dS, dTFr_dp)
  real,  intent(in)  :: S         !< salinity [ppt].
  real,  intent(in)  :: pres      !< pressure [Pa].
  real,  intent(out) :: T_Fr      !< Freezing point potential temperature [degC].
  real,  intent(in)  :: TFr_S0_P0 !< The freezing point at S=0, p=0 [degC].
  real,  intent(in)  :: dTFr_dS   !< The derivative of freezing point with salinity,
                                  !! [degC ppt-1].
  real,  intent(in)  :: dTFr_dp   !< The derivative of freezing point with pressure,
                                  !! [degC Pa-1].

  T_Fr = (TFr_S0_P0 + dTFr_dS*S) + dTFr_dp*pres

end subroutine calculate_TFreeze_linear_scalar

!>  This subroutine computes an array of freezing point potential temperatures
!!  [degC] from salinity [ppt], and pressure [Pa] using a simple
!!  linear expression, with coefficients passed in as arguments.
subroutine calculate_TFreeze_linear_array(S, pres, T_Fr, start, npts, &
                                          TFr_S0_P0, dTFr_dS, dTFr_dp)
  real,  dimension(:), intent(in)  :: S         !< salinity [ppt].
  real,  dimension(:), intent(in)  :: pres      !< pressure [Pa].
  real,  dimension(:), intent(out) :: T_Fr      !< Freezing point potential temperature [degC].
  integer,             intent(in)  :: start     !< the starting point in the arrays.
  integer,             intent(in)  :: npts      !< the number of values to calculate.
  real,                intent(in)  :: TFr_S0_P0 !< The freezing point at S=0, p=0, [degC].
  real,                intent(in)  :: dTFr_dS   !< The derivative of freezing point with salinity,
                                                !! [degC PSU-1].
  real,                intent(in)  :: dTFr_dp   !< The derivative of freezing point with pressure,
                                                !! [degC Pa-1].
  integer :: j

  do j=start,start+npts-1
    T_Fr(j) = (TFr_S0_P0 + dTFr_dS*S(j)) + dTFr_dp*pres(j)
  enddo

end subroutine calculate_TFreeze_linear_array

!> This subroutine computes the freezing point potential temperature
!! [degC] from salinity [ppt], and pressure [Pa] using the expression
!! from Millero (1978) (and in appendix A of Gill 1982), but with the of the
!! pressure dependence changed from 7.53e-8 to 7.75e-8 to make this an
!! expression for potential temperature (not in situ temperature), using a
!! value that is correct at the freezing point at 35 PSU and 5e6 Pa (500 dbar).
subroutine calculate_TFreeze_Millero_scalar(S, pres, T_Fr)
  real,    intent(in)  :: S    !< Salinity in PSU.
  real,    intent(in)  :: pres !< Pressure [Pa].
  real,    intent(out) :: T_Fr !< Freezing point potential temperature [degC].

  ! Local variables
  real, parameter :: cS1 = -0.0575, cS3_2 = 1.710523e-3, cS2 = -2.154996e-4
  real, parameter :: dTFr_dp = -7.75e-8

  T_Fr = S*(cS1 + (cS3_2 * sqrt(max(S,0.0)) + cS2 * S)) + dTFr_dp*pres

end subroutine calculate_TFreeze_Millero_scalar

!> This subroutine computes the freezing point potential temperature
!! [degC] from salinity [ppt], and pressure [Pa] using the expression
!! from Millero (1978) (and in appendix A of Gill 1982), but with the of the
!! pressure dependence changed from 7.53e-8 to 7.75e-8 to make this an
!! expression for potential temperature (not in situ temperature), using a
!! value that is correct at the freezing point at 35 PSU and 5e6 Pa (500 dbar).
subroutine calculate_TFreeze_Millero_array(S, pres, T_Fr, start, npts)
  real,  dimension(:), intent(in)  :: S     !< Salinity [PSU].
  real,  dimension(:), intent(in)  :: pres  !< Pressure [Pa].
  real,  dimension(:), intent(out) :: T_Fr  !< Freezing point potential temperature [degC].
  integer,             intent(in)  :: start !< The starting point in the arrays.
  integer,             intent(in)  :: npts  !< The number of values to calculate.

  ! Local variables
  real, parameter :: cS1 = -0.0575, cS3_2 = 1.710523e-3, cS2 = -2.154996e-4
  real, parameter :: dTFr_dp = -7.75e-8
  integer :: j

  do j=start,start+npts-1
    T_Fr(j) = S(j)*(cS1 + (cS3_2 * sqrt(max(S(j),0.0)) + cS2 * S(j))) + &
              dTFr_dp*pres(j)
  enddo

end subroutine calculate_TFreeze_Millero_array

!> This subroutine computes the freezing point conservative temperature
!! [degC] from absolute salinity [g/kg], and pressure [Pa] using the
!! TEOS10 package.
subroutine calculate_TFreeze_teos10_scalar(S, pres, T_Fr)
  real,    intent(in)  :: S    !< Absolute salinity [g/kg].
  real,    intent(in)  :: pres !< Pressure [Pa].
  real,    intent(out) :: T_Fr !< Freezing point conservative temperature [degC].

  ! Local variables
  real, dimension(1) :: S0, pres0
  real, dimension(1) :: tfr0

  S0(1) = S
  pres0(1) = pres

  call calculate_TFreeze_teos10_array(S0, pres0, tfr0, 1, 1)
  T_Fr = tfr0(1)

end subroutine calculate_TFreeze_teos10_scalar

!> This subroutine computes the freezing point conservative temperature
!! [degC] from absolute salinity [g/kg], and pressure [Pa] using the
!! TEOS10 package.
subroutine calculate_TFreeze_teos10_array(S, pres, T_Fr, start, npts)
  real, dimension(:), intent(in)  :: S     !< absolute salinity [g/kg].
  real, dimension(:), intent(in)  :: pres  !< pressure [Pa].
  real, dimension(:), intent(out) :: T_Fr  !< Freezing point conservative temperature [degC].
  integer,            intent(in)  :: start !< the starting point in the arrays.
  integer,            intent(in)  :: npts  !< the number of values to calculate.

  ! Local variables
  real, parameter :: Pa2db  = 1.e-4  ! The conversion factor from Pa to dbar.
  real :: zs,zp
  integer :: j
  ! Assume sea-water contains no dissolved air.
  real, parameter :: saturation_fraction = 0.0

  do j=start,start+npts-1
    !Conversions
    zs = S(j)
    zp = pres(j)* Pa2db         !Convert pressure from Pascal to decibar

    if (S(j) < -1.0e-10) cycle !Can we assume safely that this is a missing value?
    T_Fr(j) = gsw_ct_freezing_exact(zs,zp,saturation_fraction)
 enddo

end subroutine calculate_TFreeze_teos10_array

end module MOM_TFreeze
