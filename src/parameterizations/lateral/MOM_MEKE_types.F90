module MOM_MEKE_types
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
!*                                                                     *
!*    This program contains the subroutine that calculates the         *
!*  effects of horizontal viscosity, including parameterizations of    *
!*  the value of the viscosity itself. mesosclae_EKE calculates        *
!*  the evolution of sub-grid scale mesoscale EKE.                     *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

implicit none ; private

type, public :: MEKE_type
  ! Variables
  real, dimension(:,:), pointer :: &
    MEKE => NULL(), &   ! Vertically averaged eddy kinetic energy, in m2 s-2.
    GM_src => NULL(), & ! MEKE source due to thickness mixing (GM), in W m-2.
    mom_src => NULL(),& ! MEKE source from lateral friction in the momentum
                        ! equations, in W m-2.
    Kh => NULL(), &     ! The MEKE-derived lateral mixing coefficient in m2 s-1.
    Rd_dx_h => NULL(), &! The deformation radius compared with the grid
                        ! spacing, copied from VarMix_CS, nondim.
    Ku => NULL()        ! The MEKE-derived lateral viscosity coefficient in m2 s-1.
                        ! This viscosity can be negative when representing backscatter
                        ! from unresolved eddies (see Jansen and Held, 2014).
  ! Parameters
  real :: KhTh_fac = 1.0 ! Multiplier to map Kh(MEKE) to KhTh, nondim
  real :: KhTr_fac = 1.0 ! Multiplier to map Kh(MEKE) to KhTr, nondim.
  real :: backscatter_Ro_pow = 0.0 ! Power in Rossby number function for backscatter.
  real :: backscatter_Ro_c = 0.0 ! Coefficient in Rossby number function for backscatter.
end type MEKE_type

end module MOM_MEKE_types
