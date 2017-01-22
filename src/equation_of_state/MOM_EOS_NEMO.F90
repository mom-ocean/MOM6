module MOM_EOS_NEMO
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

!***********************************************************************
!*  The subroutines in this file implement the equation of state for   *
!*  sea water using the formulae provided by NEMO developer Roquet     *
!*  in a private communication , Roquet et al, Ocean Modelling (2015)  *
!*  Roquet, F., Madec, G., McDougall, T. J., and Barker, P. M., 2015.  *
!*  Accurate polynomial expressions for the density and specific volume*
!*  of seawater using the TEOS-10 standard. Ocean Modelling, 90:29-43. *
!*  These algorithms are NOT from NEMO package!!                       *
!***********************************************************************

!use gsw_mod_toolbox, only : gsw_sr_from_sp, gsw_ct_from_pt
use gsw_mod_toolbox, only : gsw_rho_first_derivatives

implicit none ; private

public calculate_compress_nemo, calculate_density_nemo
public calculate_density_derivs_nemo
public calculate_density_scalar_nemo, calculate_density_array_nemo

interface calculate_density_nemo
  module procedure calculate_density_scalar_nemo, calculate_density_array_nemo
end interface calculate_density_nemo

   real, parameter :: Pa2db  = 1.e-4
   real, parameter :: rdeltaS = 32.
   real, parameter :: r1_S0  = 0.875/35.16504
   real, parameter :: r1_T0  = 1./40.
   real, parameter :: r1_P0  = 1.e-4
   real, parameter :: R00 = 4.6494977072e+01
   real, parameter :: R01 = -5.2099962525
   real, parameter :: R02 = 2.2601900708e-01
   real, parameter :: R03 = 6.4326772569e-02
   real, parameter :: R04 = 1.5616995503e-02
   real, parameter :: R05 = -1.7243708991e-03
   real, parameter :: EOS000 = 8.0189615746e+02
   real, parameter :: EOS100 = 8.6672408165e+02
   real, parameter :: EOS200 = -1.7864682637e+03
   real, parameter :: EOS300 = 2.0375295546e+03
   real, parameter :: EOS400 = -1.2849161071e+03
   real, parameter :: EOS500 = 4.3227585684e+02
   real, parameter :: EOS600 = -6.0579916612e+01
   real, parameter :: EOS010 = 2.6010145068e+01
   real, parameter :: EOS110 = -6.5281885265e+01
   real, parameter :: EOS210 = 8.1770425108e+01
   real, parameter :: EOS310 = -5.6888046321e+01
   real, parameter :: EOS410 = 1.7681814114e+01
   real, parameter :: EOS510 = -1.9193502195
   real, parameter :: EOS020 = -3.7074170417e+01
   real, parameter :: EOS120 = 6.1548258127e+01
   real, parameter :: EOS220 = -6.0362551501e+01
   real, parameter :: EOS320 = 2.9130021253e+01
   real, parameter :: EOS420 = -5.4723692739
   real, parameter :: EOS030 = 2.1661789529e+01
   real, parameter :: EOS130 = -3.3449108469e+01
   real, parameter :: EOS230 = 1.9717078466e+01
   real, parameter :: EOS330 = -3.1742946532
   real, parameter :: EOS040 = -8.3627885467
   real, parameter :: EOS140 = 1.1311538584e+01
   real, parameter :: EOS240 = -5.3563304045
   real, parameter :: EOS050 = 5.4048723791e-01
   real, parameter :: EOS150 = 4.8169980163e-01
   real, parameter :: EOS060 = -1.9083568888e-01
   real, parameter :: EOS001 = 1.9681925209e+01
   real, parameter :: EOS101 = -4.2549998214e+01
   real, parameter :: EOS201 = 5.0774768218e+01
   real, parameter :: EOS301 = -3.0938076334e+01
   real, parameter :: EOS401 = 6.6051753097
   real, parameter :: EOS011 = -1.3336301113e+01
   real, parameter :: EOS111 = -4.4870114575
   real, parameter :: EOS211 = 5.0042598061
   real, parameter :: EOS311 = -6.5399043664e-01
   real, parameter :: EOS021 = 6.7080479603
   real, parameter :: EOS121 = 3.5063081279
   real, parameter :: EOS221 = -1.8795372996
   real, parameter :: EOS031 = -2.4649669534
   real, parameter :: EOS131 = -5.5077101279e-01
   real, parameter :: EOS041 = 5.5927935970e-01
   real, parameter :: EOS002 = 2.0660924175
   real, parameter :: EOS102 = -4.9527603989
   real, parameter :: EOS202 = 2.5019633244
   real, parameter :: EOS012 = 2.0564311499
   real, parameter :: EOS112 = -2.1311365518e-01
   real, parameter :: EOS022 = -1.2419983026
   real, parameter :: EOS003 = -2.3342758797e-02
   real, parameter :: EOS103 = -1.8507636718e-02
   real, parameter :: EOS013 = 3.7969820455e-01
   real, parameter :: ALP000 = -6.5025362670e-01
   real, parameter :: ALP100 = 1.6320471316
   real, parameter :: ALP200 = -2.0442606277
   real, parameter :: ALP300 = 1.4222011580
   real, parameter :: ALP400 = -4.4204535284e-01
   real, parameter :: ALP500 = 4.7983755487e-02
   real, parameter :: ALP010 = 1.8537085209
   real, parameter :: ALP110 = -3.0774129064
   real, parameter :: ALP210 = 3.0181275751
   real, parameter :: ALP310 = -1.4565010626
   real, parameter :: ALP410 = 2.7361846370e-01
   real, parameter :: ALP020 = -1.6246342147
   real, parameter :: ALP120 = 2.5086831352
   real, parameter :: ALP220 = -1.4787808849
   real, parameter :: ALP320 = 2.3807209899e-01
   real, parameter :: ALP030 = 8.3627885467e-01
   real, parameter :: ALP130 = -1.1311538584
   real, parameter :: ALP230 = 5.3563304045e-01
   real, parameter :: ALP040 = -6.7560904739e-02
   real, parameter :: ALP140 = -6.0212475204e-02
   real, parameter :: ALP050 = 2.8625353333e-02
   real, parameter :: ALP001 = 3.3340752782e-01
   real, parameter :: ALP101 = 1.1217528644e-01
   real, parameter :: ALP201 = -1.2510649515e-01
   real, parameter :: ALP301 = 1.6349760916e-02
   real, parameter :: ALP011 = -3.3540239802e-01
   real, parameter :: ALP111 = -1.7531540640e-01
   real, parameter :: ALP211 = 9.3976864981e-02
   real, parameter :: ALP021 = 1.8487252150e-01
   real, parameter :: ALP121 = 4.1307825959e-02
   real, parameter :: ALP031 = -5.5927935970e-02
   real, parameter :: ALP002 = -5.1410778748e-02
   real, parameter :: ALP102 = 5.3278413794e-03
   real, parameter :: ALP012 = 6.2099915132e-02
   real, parameter :: ALP003 = -9.4924551138e-03
   real, parameter :: BET000 = 1.0783203594e+01
   real, parameter :: BET100 = -4.4452095908e+01
   real, parameter :: BET200 = 7.6048755820e+01
   real, parameter :: BET300 = -6.3944280668e+01
   real, parameter :: BET400 = 2.6890441098e+01
   real, parameter :: BET500 = -4.5221697773
   real, parameter :: BET010 = -8.1219372432e-01
   real, parameter :: BET110 = 2.0346663041
   real, parameter :: BET210 = -2.1232895170
   real, parameter :: BET310 = 8.7994140485e-01
   real, parameter :: BET410 = -1.1939638360e-01
   real, parameter :: BET020 = 7.6574242289e-01
   real, parameter :: BET120 = -1.5019813020
   real, parameter :: BET220 = 1.0872489522
   real, parameter :: BET320 = -2.7233429080e-01
   real, parameter :: BET030 = -4.1615152308e-01
   real, parameter :: BET130 = 4.9061350869e-01
   real, parameter :: BET230 = -1.1847737788e-01
   real, parameter :: BET040 = 1.4073062708e-01
   real, parameter :: BET140 = -1.3327978879e-01
   real, parameter :: BET050 = 5.9929880134e-03
   real, parameter :: BET001 = -5.2937873009e-01
   real, parameter :: BET101 = 1.2634116779
   real, parameter :: BET201 = -1.1547328025
   real, parameter :: BET301 = 3.2870876279e-01
   real, parameter :: BET011 = -5.5824407214e-02
   real, parameter :: BET111 = 1.2451933313e-01
   real, parameter :: BET211 = -2.4409539932e-02
   real, parameter :: BET021 = 4.3623149752e-02
   real, parameter :: BET121 = -4.6767901790e-02
   real, parameter :: BET031 = -6.8523260060e-03
   real, parameter :: BET002 = -6.1618945251e-02
   real, parameter :: BET102 = 6.2255521644e-02
   real, parameter :: BET012 = -2.6514181169e-03
   real, parameter :: BET003 = -2.3025968587e-04



contains

subroutine calculate_density_scalar_nemo(T, S, pressure, rho)
real,    intent(in)  :: T, S, pressure
real,    intent(out) :: rho
! * Arguments: T - conservative temperature in C.                      *
! *  (in)      S - absoulte salinity in g/Kg.                          *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     rho - in situ density in kg m-3.                        *

! *====================================================================*
! *  This subroutine computes the in situ density of sea water (rho in *
! *  units of kg/m^3) from absolute salinity (S in g/Kg), conservative *
! *  temperature (T in deg C), and pressure in Pa.                     *
! *====================================================================*

  real :: al0, p0, lambda
  integer :: j
  real, dimension(1) :: T0, S0, pressure0
  real, dimension(1) :: rho0

  T0(1) = T
  S0(1) = S
  pressure0(1) = pressure

  call calculate_density_array_nemo(T0, S0, pressure0, rho0, 1, 1)
  rho = rho0(1)
end subroutine calculate_density_scalar_nemo

subroutine calculate_density_array_nemo(T, S, pressure, rho, start, npts)
  real,    intent(in),  dimension(:) :: T, S, pressure
  real,    intent(out), dimension(:) :: rho
  integer, intent(in)                :: start, npts
! * Arguments: T - conservative temperature in C.                      *
! *  (in)      S - absoulte salinity in g/Kg.                          *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     rho - in situ density in kg m-3.                        *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *

! *====================================================================*
! *  This subroutine computes the in situ density of sea water (rho in *
! *  units of kg/m^3) from absolute salinity (S in g/Kg),              *
! *  conservative temperature (T in deg C), and pressure in Pa.        *
! *====================================================================*
  real :: zp,zt , zh , zs , zr0, zn , zn0, zn1, zn2, zn3
  integer :: j

  do j=start,start+npts-1
    !Conversions
    zs = S(j) !gsw_sr_from_sp(S(j))       !Convert practical salinity to absolute salinity
    zt = T(j) !gsw_ct_from_pt(S(j),T(j))  !Convert potantial temp to conservative temp
    zp = pressure(j)* Pa2db         !Convert pressure from Pascal to decibar

    !The following algorithm was provided by Roquet in a private communication.
    !It is not necessarily the algorithm used in NEMO ocean!
    zp  = zp * r1_P0 !pressure
    zt  = zt * r1_T0 !temperature
    zs  = SQRT( ABS( zs + rdeltaS ) * r1_S0 )   ! square root salinity
    !
    zn3 = EOS013*zt   &
       &   + EOS103*zs+EOS003
       !
    zn2 = (EOS022*zt   &
       &   + EOS112*zs+EOS012)*zt   &
       &   + (EOS202*zs+EOS102)*zs+EOS002
       !
    zn1 = (((EOS041*zt   &
       &   + EOS131*zs+EOS031)*zt   &
       &   + (EOS221*zs+EOS121)*zs+EOS021)*zt   &
       &   + ((EOS311*zs+EOS211)*zs+EOS111)*zs+EOS011)*zt   &
       &   + (((EOS401*zs+EOS301)*zs+EOS201)*zs+EOS101)*zs+EOS001
       !
    zn0 = (((((EOS060*zt   &
       &   + EOS150*zs+EOS050)*zt   &
       &   + (EOS240*zs+EOS140)*zs+EOS040)*zt   &
       &   + ((EOS330*zs+EOS230)*zs+EOS130)*zs+EOS030)*zt   &
       &   + (((EOS420*zs+EOS320)*zs+EOS220)*zs+EOS120)*zs+EOS020)*zt   &
       &   + ((((EOS510*zs+EOS410)*zs+EOS310)*zs+EOS210)*zs+EOS110)*zs+EOS010)*zt   &
       &   + (((((EOS600*zs+EOS500)*zs+EOS400)*zs+EOS300)*zs+EOS200)*zs+EOS100)*zs+EOS000
       !
    zn  = ( ( zn3 * zp + zn2 ) * zp + zn1 ) * zp + zn0
    !
    zr0 = (((((R05 * zp+R04) * zp+R03 ) * zp+R02 ) * zp+R01) * zp+R00) * zp
    !
    rho(j) =  ( zn + zr0 ) ! density

 enddo
end subroutine calculate_density_array_nemo

subroutine calculate_density_derivs_nemo(T, S, pressure, drho_dT, drho_dS, start, npts)
  real,    intent(in),  dimension(:) ::  T, S, pressure
  real,    intent(out), dimension(:) :: drho_dT, drho_dS
  integer, intent(in)                :: start, npts
! * Arguments: T - conservative temperature in C.                      *
! *  (in)      S - absolute salinity in g/kg.                          *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     drho_dT - the partial derivative of density with        *
! *                      potential temperature, in kg m-3 K-1.         *
! *  (out)     drho_dS - the partial derivative of density with        *
! *                      salinity, in kg m-3 psu-1.                    *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
  real :: zp,zt , zh , zs , zr0, zn , zn0, zn1, zn2, zn3
  integer :: j

  do j=start,start+npts-1
    !Conversions
    zs = S(j) !gsw_sr_from_sp(S(j))       !Convert practical salinity to absolute salinity
    zt = T(j) !gsw_ct_from_pt(S(j),T(j))  !Convert potantial temp to conservative temp
    zp = pressure(j)* Pa2db         !Convert pressure from Pascal to decibar

    !The following algorithm was provided by Roquet in a private communication.
    !It is not necessarily the algorithm used in NEMO ocean!
    zp  = zp * r1_P0  ! pressure (first converted to decibar)
    zt  = zt * r1_T0                ! temperature
    zs  = SQRT( ABS( zs + rdeltaS ) * r1_S0 )   ! square root salinity
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
    drho_dT(j) = -zn
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
    !
    drho_dS(j) = zn / zs
  enddo

end subroutine calculate_density_derivs_nemo

subroutine calculate_compress_nemo(T, S, pressure, rho, drho_dp, start, npts)
  real,    intent(in),  dimension(:) :: T, S, pressure
  real,    intent(out), dimension(:) :: rho, drho_dp
  integer, intent(in)                :: start, npts
! * Arguments: T - conservative temperature in C.                      *
! *  (in)      S - absolute salinity in g/kg.                          *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     rho - in situ density in kg m-3.                        *
! *  (out)     drho_dp - the partial derivative of density with        *
! *                      pressure (also the inverse of the square of   *
! *                      sound speed) in s2 m-2.                       *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
! *====================================================================*
  real ::  zs,zt,zp
  integer :: j

  call calculate_density_array_nemo(T, S, pressure, rho, start, npts)
  !
  !NOTE: The following calculates the TEOS10 approximation to compressibility
  !      since the corresponding NEMO approximation is not available yet.
  !
  do j=start,start+npts-1
   !Conversions
    zs = S(j) !gsw_sr_from_sp(S(j))       !Convert practical salinity to absolute salinity
    zt = T(j) !gsw_ct_from_pt(S(j),T(j))  !Convert potantial temp to conservative temp
    zp = pressure(j)* Pa2db         !Convert pressure from Pascal to decibar
    call gsw_rho_first_derivatives(zs,zt,zp, drho_dp=drho_dp(j))
 enddo
end subroutine calculate_compress_nemo

end module MOM_EOS_NEMO
