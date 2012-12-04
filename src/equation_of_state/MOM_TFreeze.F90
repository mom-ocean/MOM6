module MOM_TFreeze

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

!********+*********+*********+*********+*********+*********+*********+**
!*  The subroutines in this file determine the potential temperature   *
!* at which sea-water freezes.                                         *
!********+*********+*********+*********+*********+*********+*********+**

implicit none ; private

public calculate_TFreeze_linear, calculate_TFreeze_Millero

interface calculate_TFreeze_linear
  module procedure calculate_TFreeze_linear_scalar, calculate_TFreeze_linear_array
end interface calculate_TFreeze_linear

interface calculate_TFreeze_Millero
  module procedure calculate_TFreeze_Millero_scalar, calculate_TFreeze_Millero_array
end interface calculate_TFreeze_Millero

contains

subroutine calculate_TFreeze_linear_scalar(S, pres, T_Fr, TFr_S0_P0, &
                                           dTFr_dS, dTFr_dp)
  real,    intent(in)  :: S, pres
  real,    intent(out) :: T_Fr
  real,    intent(in)  :: TFr_S0_P0, dTFr_dS, dTFr_dp
!    This subroutine computes the freezing point potential temparature
!  (in deg C) from salinity (in psu), and pressure (in Pa) using a simple
!  linear expression, with coefficients passed in as arguments.
!
! Arguments: S - salinity in PSU.
!  (in)      pres - pressure in Pa.
!  (out)     T_Fr - Freezing point potential temperature in deg C.
!  (in)      TFr_S0_P0 - The freezing point at S=0, p=0, in deg C.
!  (in)      dTFr_dS - The derivatives of freezing point with salinity, in
!                      deg C PSU-1.
!  (in)      dTFr_dp - The derivatives of freezing point with pressure, in
!                      deg C Pa-1.

  T_Fr = (TFr_S0_P0 + dTFr_dS*S) + dTFr_dp*pres

end subroutine calculate_TFreeze_linear_scalar

subroutine calculate_TFreeze_linear_array(S, pres, T_Fr, start, npts, &
                                          TFr_S0_P0, dTFr_dS, dTFr_dp)
  real,  dimension(:), intent(in)  :: S, pres
  real,  dimension(:), intent(out) :: T_Fr
  integer,             intent(in)  :: start, npts
  real,                intent(in)  :: TFr_S0_P0, dTFr_dS, dTFr_dp
!    This subroutine computes the freezing point potential temparature
!  (in deg C) from salinity (in psu), and pressure (in Pa) using a simple
!  linear expression, with coefficients passed in as arguments.
!
! Arguments: S - salinity in PSU.
!  (in)      pres - pressure in Pa.
!  (out)     T_Fr - Freezing point potential temperature in deg C.
!  (in)      start - the starting point in the arrays.
!  (in)      npts - the number of values to calculate.
!  (in)      TFr_S0_P0 - The freezing point at S=0, p=0, in deg C.
!  (in)      dTFr_dS - The derivative of freezing point with salinity, in
!                      deg C PSU-1.
!  (in)      dTFr_dp - The derivative of freezing point with pressure, in
!                      deg C Pa-1.
  integer :: j

  do j=start,start+npts-1
    T_Fr(j) = (TFr_S0_P0 + dTFr_dS*S(j)) + dTFr_dp*pres(j)
  enddo

end subroutine calculate_TFreeze_linear_array

subroutine calculate_TFreeze_Millero_scalar(S, pres, T_Fr)
  real,    intent(in)  :: S, pres
  real,    intent(out) :: T_Fr
!    This subroutine computes the freezing point potential temparature
!  (in deg C) from salinity (in psu), and pressure (in Pa) using the expression
!  from Millero (1978) (and in appendix A of Gill 1982), but with the of the
!  pressure dependence changed from 7.53e-8 to 7.75e-8 to make this an
!  expression for potential temperature (not in situ temperature), using a
!  value that is correct at the freezing point at 35 PSU and 5e6 Pa (500 dbar).
!
! Arguments: S - salinity in PSU.
!  (in)      pres - pressure in Pa.
!  (out)     T_Fr - Freezing point potential temperature in deg C.
  real, parameter :: cS1 = -0.0575, cS3_2 = 1.710523e-3, cS2 = -2.154996e-4
  real, parameter :: dTFr_dp = -7.75e-8

  T_Fr = S*(cS1 + (cS3_2 * sqrt(max(S,0.0)) + cS2 * S)) + dTFr_dp*pres

end subroutine calculate_TFreeze_Millero_scalar


subroutine calculate_TFreeze_Millero_array(S, pres, T_Fr, start, npts)
  real,  dimension(:), intent(in)  :: S, pres
  real,  dimension(:), intent(out) :: T_Fr
  integer,             intent(in)  :: start, npts
!    This subroutine computes the freezing point potential temparature
!  (in deg C) from salinity (in psu), and pressure (in Pa) using the expression
!  from Millero (1978) (and in appendix A of Gill 1982), but with the of the
!  pressure dependence changed from 7.53e-8 to 7.75e-8 to make this an
!  expression for potential temperature (not in situ temperature), using a
!  value that is correct at the freezing point at 35 PSU and 5e6 Pa (500 dbar).
!
! Arguments: S - salinity in PSU.
!  (in)      pres - pressure in Pa.
!  (out)     T_Fr - Freezing point potential temperature in deg C.
!  (in)      start - the starting point in the arrays.
!  (in)      npts - the number of values to calculate.
  real, parameter :: cS1 = -0.0575, cS3_2 = 1.710523e-3, cS2 = -2.154996e-4
  real, parameter :: dTFr_dp = -7.75e-8
  integer :: j

  do j=start,start+npts-1
    T_Fr(j) = S(j)*(cS1 + (cS3_2 * sqrt(max(S(j),0.0)) + cS2 * S(j))) + &
              dTFr_dp*pres(j)
  enddo

end subroutine calculate_TFreeze_Millero_array


end module MOM_TFreeze
