!> Freezing point expressions
module MOM_TFreeze

! This file is part of MOM6. See LICENSE.md for the license.

!********+*********+*********+*********+*********+*********+*********+**
!*  The subroutines in this file determine the potential temperature   *
!* or conservative temperature at which sea-water freezes.             *
!********+*********+*********+*********+*********+*********+*********+**
use gsw_mod_toolbox, only : gsw_ct_freezing_exact

implicit none ; private

public calculate_TFreeze_linear, calculate_TFreeze_Millero, calculate_TFreeze_teos10
public calculate_TFreeze_TEOS_poly

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

!> Compute the freezing point conservative temperature [degC] from absolute salinity [g kg-1]
!! and pressure [Pa] using the TEOS10 package.
interface calculate_TFreeze_teos10
  module procedure calculate_TFreeze_teos10_scalar, calculate_TFreeze_teos10_array
end interface calculate_TFreeze_teos10

!> Compute the freezing point conservative temperature [degC] from absolute salinity [g kg-1] and
!! pressure [Pa] using a rescaled and refactored version of the expressions from the TEOS10 package.
interface calculate_TFreeze_TEOS_poly
  module procedure calculate_TFreeze_TEOS_poly_scalar, calculate_TFreeze_TEOS_poly_array
end interface calculate_TFreeze_TEOS_poly

contains

!>  This subroutine computes the freezing point potential temperature [degC] from
!!  salinity [ppt], and pressure [Pa] using a simple linear expression,
!!  with coefficients passed in as arguments.
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
                                                !! [degC ppt-1].
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
  real,    intent(in)  :: S    !< Salinity [PSU]
  real,    intent(in)  :: pres !< Pressure [Pa]
  real,    intent(out) :: T_Fr !< Freezing point potential temperature [degC]

  ! Local variables
  real, parameter :: cS1 = -0.0575      ! A term in the freezing point fit [degC PSU-1]
  real, parameter :: cS3_2 = 1.710523e-3 ! A term in the freezing point fit [degC PSU-3/2]
  real, parameter :: cS2 = -2.154996e-4 ! A term in the freezing point fit [degC PSU-2]
  real, parameter :: dTFr_dp = -7.75e-8 ! Derivative of freezing point with pressure [degC Pa-1]

  T_Fr = S*(cS1 + (cS3_2 * sqrt(max(S, 0.0)) + cS2 * S)) + dTFr_dp*pres

end subroutine calculate_TFreeze_Millero_scalar

!> This subroutine computes the freezing point potential temperature
!! [degC] from salinity [ppt], and pressure [Pa] using the expression
!! from Millero (1978) (and in appendix A of Gill 1982), but with the
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
  real, parameter :: cS1 = -0.0575      ! A term in the freezing point fit [degC PSU-1]
  real, parameter :: cS3_2 = 1.710523e-3 ! A term in the freezing point fit [degC PSU-3/2]
  real, parameter :: cS2 = -2.154996e-4 ! A term in the freezing point fit [degC PSU-2]
  real, parameter :: dTFr_dp = -7.75e-8 ! Derivative of freezing point with pressure [degC Pa-1]
  integer :: j

  do j=start,start+npts-1
    T_Fr(j) = S(j)*(cS1 + (cS3_2 * sqrt(max(S(j), 0.0)) + cS2 * S(j))) + &
              dTFr_dp*pres(j)
  enddo

end subroutine calculate_TFreeze_Millero_array

!> This subroutine computes the freezing point conservative temperature [degC]
!! from absolute salinity [g kg-1], and pressure [Pa] using a rescaled and
!! refactored version of the polynomial expressions from the TEOS10 package.
subroutine calculate_TFreeze_TEOS_poly_scalar(S, pres, T_Fr)
  real,    intent(in)  :: S    !< Absolute salinity [g kg-1].
  real,    intent(in)  :: pres !< Pressure [Pa].
  real,    intent(out) :: T_Fr !< Freezing point conservative temperature [degC].

  ! Local variables
  real, dimension(1) :: S0    ! Salinity at a point [g kg-1]
  real, dimension(1) :: pres0 ! Pressure at a point [Pa]
  real, dimension(1) :: tfr0  ! The freezing temperature [degC]

  S0(1) = S
  pres0(1) = pres

  call calculate_TFreeze_TEOS_poly_array(S0, pres0, tfr0, 1, 1)
  T_Fr = tfr0(1)

end subroutine calculate_TFreeze_TEOS_poly_scalar

!> This subroutine computes the freezing point conservative temperature [degC]
!! from absolute salinity [g kg-1], and pressure [Pa] using a rescaled and
!! refactored version of the polynomial expressions from the TEOS10 package.
subroutine calculate_TFreeze_TEOS_poly_array(S, pres, T_Fr, start, npts)
  real, dimension(:), intent(in)  :: S     !< absolute salinity [g kg-1].
  real, dimension(:), intent(in)  :: pres  !< Pressure [Pa].
  real, dimension(:), intent(out) :: T_Fr  !< Freezing point conservative temperature [degC].
  integer,            intent(in)  :: start !< The starting point in the arrays
  integer,            intent(in)  :: npts  !< The number of values to calculate

  ! Local variables
  real :: Sa    ! Absolute salinity [g kg-1] = [ppt]
  real :: rS    ! Square root of salinity [ppt1/2]
  ! The coefficients here use the notation TFab for contributions proportional to S**a/2 * P**b.
  real, parameter :: TF00 =  0.017947064327968736  ! Freezing point coefficient [degC]
  real, parameter :: TF20 = -6.076099099929818e-2  ! Freezing point coefficient [degC ppt-1]
  real, parameter :: TF30 =  4.883198653547851e-3  ! Freezing point coefficient [degC ppt-3/2]
  real, parameter :: TF40 = -1.188081601230542e-3  ! Freezing point coefficient [degC ppt-2]
  real, parameter :: TF50 =  1.334658511480257e-4  ! Freezing point coefficient [degC ppt-5/2]
  real, parameter :: TF60 = -8.722761043208607e-6  ! Freezing point coefficient [degC ppt-3]
  real, parameter :: TF70 =  2.082038908808201e-7  ! Freezing point coefficient [degC ppt-7/2]
  real, parameter :: TF01 = -7.389420998107497e-8  ! Freezing point coefficient [degC Pa-1]
  real, parameter :: TF21 = -9.891538123307282e-11 ! Freezing point coefficient [degC ppt-1 Pa-1]
  real, parameter :: TF31 = -8.987150128406496e-13 ! Freezing point coefficient [degC ppt-3/2 Pa-1]
  real, parameter :: TF41 =  1.054318231187074e-12 ! Freezing point coefficient [degC ppt-2 Pa-1]
  real, parameter :: TF51 =  3.850133554097069e-14 ! Freezing point coefficient [degC ppt-5/2 Pa-1]
  real, parameter :: TF61 = -2.079022768390933e-14 ! Freezing point coefficient [degC ppt-3 Pa-1]
  real, parameter :: TF71 =  1.242891021876471e-15 ! Freezing point coefficient [degC ppt-7/2 Pa-1]
  real, parameter :: TF02 = -2.110913185058476e-16 ! Freezing point coefficient [degC Pa-2]
  real, parameter :: TF22 =  3.831132432071728e-19 ! Freezing point coefficient [degC ppt-1 Pa-2]
  real, parameter :: TF32 =  1.065556599652796e-19 ! Freezing point coefficient [degC ppt-3/2 Pa-2]
  real, parameter :: TF42 = -2.078616693017569e-20 ! Freezing point coefficient [degC ppt-2 Pa-2]
  real, parameter :: TF52 =  1.596435439942262e-21 ! Freezing point coefficient [degC ppt-5/2 Pa-2]
  real, parameter :: TF03 =  2.295491578006229e-25 ! Freezing point coefficient [degC Pa-3]
  real, parameter :: TF23 = -7.997496801694032e-27 ! Freezing point coefficient [degC ppt-1 Pa-3]
  real, parameter :: TF33 =  8.756340772729538e-28 ! Freezing point coefficient [degC ppt-3/2 Pa-3]
  real, parameter :: TF43 =  1.338002171109174e-29 ! Freezing point coefficient [degC ppt-2 Pa-3]
  integer :: j

  do j=start,start+npts-1
    rS = sqrt(max(S(j), 0.0))
    T_Fr(j) =       (TF00 + S(j)*(TF20 + rS*(TF30 + rS*(TF40 + rS*(TF50 + rS*(TF60 + rS*TF70)))))) &
        + pres(j)*( (TF01 + S(j)*(TF21 + rS*(TF31 + rS*(TF41 + rS*(TF51 + rS*(TF61 + rS*TF71)))))) &
         + pres(j)*((TF02 + S(j)*(TF22 + rS*(TF32 + rS*(TF42 + rS* TF52)))) &
          + pres(j)*(TF03 + S(j)*(TF23 + rS*(TF33 + rS* TF43))) ) )
  enddo

end subroutine calculate_TFreeze_TEOS_poly_array

!> This subroutine computes the freezing point conservative temperature [degC]
!! from absolute salinity [g kg-1], and pressure [Pa] using the
!! TEOS10 package.
subroutine calculate_TFreeze_teos10_scalar(S, pres, T_Fr)
  real,    intent(in)  :: S    !< Absolute salinity [g kg-1].
  real,    intent(in)  :: pres !< Pressure [Pa].
  real,    intent(out) :: T_Fr !< Freezing point conservative temperature [degC].

  ! Local variables
  real, dimension(1) :: S0    ! Salinity at a point [g kg-1]
  real, dimension(1) :: pres0 ! Pressure at a point [Pa]
  real, dimension(1) :: tfr0  ! The freezing temperature [degC]

  S0(1) = S
  pres0(1) = pres

  call calculate_TFreeze_teos10_array(S0, pres0, tfr0, 1, 1)
  T_Fr = tfr0(1)

end subroutine calculate_TFreeze_teos10_scalar

!> This subroutine computes the freezing point conservative temperature [degC]
!! from absolute salinity [g kg-1], and pressure [Pa] using the
!! TEOS10 package.
subroutine calculate_TFreeze_teos10_array(S, pres, T_Fr, start, npts)
  real, dimension(:), intent(in)  :: S     !< absolute salinity [g kg-1].
  real, dimension(:), intent(in)  :: pres  !< pressure [Pa].
  real, dimension(:), intent(out) :: T_Fr  !< Freezing point conservative temperature [degC].
  integer,            intent(in)  :: start !< the starting point in the arrays.
  integer,            intent(in)  :: npts  !< the number of values to calculate.

  ! Local variables
  real, parameter :: Pa2db  = 1.e-4  ! The conversion factor from Pa to dbar [dbar Pa-1]
  real :: zp    ! Pressures in [dbar]
  integer :: j
  ! Assume sea-water contains no dissolved air.
  real, parameter :: saturation_fraction = 0.0 ! Air saturation fraction in seawater [nondim]

  do j=start,start+npts-1
    !Conversions
    zp = pres(j)* Pa2db         !Convert pressure from Pascal to decibar

    if (S(j) < -1.0e-10) cycle !Can we assume safely that this is a missing value?
    T_Fr(j) = gsw_ct_freezing_exact(S(j), zp, saturation_fraction)
  enddo

end subroutine calculate_TFreeze_teos10_array

end module MOM_TFreeze
