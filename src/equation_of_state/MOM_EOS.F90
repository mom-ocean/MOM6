module MOM_EOS
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
!*  The subroutines in this file are drivers for selecting different   *
!*  equations of state of sea water, provided by other modules.        *
!*                                                                     *
!*  By default, the Wright equation of state is used (EOS_DEFAULT)     *
!*                                                                     *
!***********************************************************************

use MOM_EOS_linear
use MOM_EOS_Wright
use MOM_EOS_UNESCO
use MOM_TFreeze, only : calculate_TFreeze_linear, calculate_TFreeze_Millero
use MOM_error_handler, only : MOM_error, FATAL, MOM_mesg
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_file_parser, only : uppercase
use MOM_grid, only : ocean_grid_type

implicit none ; private

#include <MOM_memory.h>

public calculate_compress, calculate_density, query_compressible
public calculate_density_derivs, calculate_2_densities
public select_eqn_of_state, deselect_eqn_of_state
public int_density_dz, int_specific_vol_dp
public calculate_TFreeze

interface calculate_density
  module procedure calculate_density_scalar, calculate_density_array
end interface calculate_density

interface calculate_TFreeze
  module procedure calculate_TFreeze_scalar, calculate_TFreeze_array
end interface calculate_TFreeze

type, public :: EOS_type ; private
  integer :: form_of_EOS = 0 ! The equation of state to use.
  integer :: form_of_TFreeze = 0 ! The expression for the potential temperature
                             ! of the freezing point.
  logical :: EOS_quadrature  ! If true, always use the generic (quadrature) 
                             ! code for the integrals of density.
  logical :: Compressible = .true. ! If true, in situ density is a function 
                             ! of pressure.
! The following parameters are used with the linear equation of state only.
  real :: Rho_T0_S0   ! The density at T=0, S=0, in kg m-3.
  real :: dRho_dT     ! The partial derivatives of density with temperature
  real :: dRho_dS     ! and salnity, in kg m-3 K-1 and kg m-3 psu-1.
! The following parameters are use with the linear expression for the freezing
! point only.
  real :: TFr_S0_P0   ! The freezing potential temperature at S=0, P=0 in deg C.
  real :: dTFr_dS ! The derivative of freezing point with salinity, in deg C PSU-1.  
  real :: dTFr_dp ! The derivative of freezing point with pressure, in deg C Pa-1.  
end type EOS_type

! The named integers that might be stored in eqn_of_state_type%form_of_EOS.
integer, parameter :: EOS_LINEAR = 1
integer, parameter :: EOS_UNESCO = 2
integer, parameter :: EOS_WRIGHT = 3

character*(10), parameter :: EOS_LINEAR_STRING = "LINEAR"
character*(10), parameter :: EOS_UNESCO_STRING = "UNESCO"
character*(10), parameter :: EOS_WRIGHT_STRING = "WRIGHT"
character*(10), parameter :: EOS_DEFAULT = EOS_WRIGHT_STRING

integer, parameter :: TFREEZE_LINEAR = 1
integer, parameter :: TFREEZE_MILLERO = 2
character*(10), parameter :: TFREEZE_LINEAR_STRING = "LINEAR"
character*(10), parameter :: TFREEZE_MILLERO_STRING = "MILLERO_78"
character*(10), parameter :: TFREEZE_DEFAULT = TFREEZE_LINEAR_STRING

contains

subroutine calculate_density_scalar(T, S, pressure, rho, start, npts, EOS)
  real,           intent(in)  :: T, S, pressure
  real,           intent(out) :: rho
  integer,        intent(in)  :: start, npts
  type(EOS_type), pointer     :: EOS
! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     rho - in situ density in kg m-3.                        *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
! *  (in)      EOS - the equation of state type.                       *
! *====================================================================*
! *  This subroutine calls the appropriate subroutine to calculate     *
! *  density for a scalar.                                             *
! *====================================================================*

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_density_scalar called with an unassociated EOS_type EOS.")

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_density_scalar_linear(T, S, pressure, rho, start, npts, &
                                      EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS)
    case (EOS_UNESCO)
      call calculate_density_scalar_unesco(T, S, pressure, rho, start, npts)
    case (EOS_WRIGHT)
      call calculate_density_scalar_wright(T, S, pressure, rho, start, npts)
    case default
      call MOM_error(FATAL, &
           "calculate_density_scalar: EOS is not valid.")
  end select

end subroutine calculate_density_scalar

subroutine calculate_density_array(T, S, pressure, rho, start, npts, EOS)
  real,    intent(in),  dimension(:) :: T, S, pressure
  real,    intent(out), dimension(:) :: rho
  integer, intent(in)                :: start, npts
  type(EOS_type),        pointer     :: EOS
! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     rho - in situ density in kg m-3.                        *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
! *  (in)      EOS - the equation of state type.                       *
! *====================================================================*
! *  This subroutine the appropriate subroutine to calculate density   *
! *  for an array.                                                     *
! *====================================================================*

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_density_array called with an unassociated EOS_type EOS.")

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_density_array_linear(T, S, pressure, rho, start, npts, &
                                      EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS)
    case (EOS_UNESCO)
      call calculate_density_array_unesco(T, S, pressure, rho, start, npts)
    case (EOS_WRIGHT)
      call calculate_density_array_wright(T, S, pressure, rho, start, npts)
    case default
      call MOM_error(FATAL, &
           "calculate_density_array: EOS%form_of_EOS is not valid.")
  end select
end subroutine calculate_density_array

subroutine calculate_TFreeze_scalar(S, pressure, T_fr, EOS)
  real,           intent(in)  :: S, pressure
  real,           intent(out) :: T_fr
  type(EOS_type), pointer     :: EOS
! * Arguments: S - salinity in PSU.                                    *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     Tfreeze - Freezing point potential temperature relative *
! *                      to 1 atmosphere pressure, in C.               *
! *  (in)      EOS - the equation of state type.                       *
! *====================================================================*
! *  This subroutine calls the appropriate subroutine to calculate the *
! *  freezing point for a scalar.                                      *
! *====================================================================*

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_TFreeze_scalar called with an unassociated EOS_type EOS.")

  select case (EOS%form_of_TFreeze)
    case (TFREEZE_LINEAR)
      call calculate_TFreeze_linear(S, pressure, T_fr, EOS%TFr_S0_P0, &
                                    EOS%dTFr_dS, EOS%dTFr_dp)
    case (TFREEZE_MILLERO)
      call calculate_TFreeze_Millero(S, pressure, T_fr)
    case default
      call MOM_error(FATAL, &
           "calculate_TFreeze_scalar: form_of_TFreeze is not valid.")
  end select

end subroutine calculate_TFreeze_scalar

subroutine calculate_TFreeze_array(S, pressure, T_fr, start, npts, EOS)
  real, dimension(:), intent(in)  :: S, pressure
  real, dimension(:), intent(out) :: T_fr
  integer,            intent(in)  :: start, npts
  type(EOS_type),     pointer     :: EOS
! * Arguments: S - salinity in PSU.                                    *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     Tfreeze - Freezing point potential temperature relative *
! *                      to 1 atmosphere pressure, in C.               *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
! *  (in)      EOS - the equation of state type.                       *
! *====================================================================*
! *  This subroutine calls the appropriate subroutine to calculate the *
! *  freezing point for a scalar.                                      *
! *====================================================================*

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_TFreeze_scalar called with an unassociated EOS_type EOS.")

  select case (EOS%form_of_TFreeze)
    case (TFREEZE_LINEAR)
      call calculate_TFreeze_linear(S, pressure, T_fr, start, npts, &
                                    EOS%TFr_S0_P0, EOS%dTFr_dS, EOS%dTFr_dp)
    case (TFREEZE_MILLERO)
      call calculate_TFreeze_Millero(S, pressure, T_fr, start, npts)
    case default
      call MOM_error(FATAL, &
           "calculate_TFreeze_scalar: form_of_TFreeze is not valid.")
  end select

end subroutine calculate_TFreeze_array

subroutine calculate_density_derivs(T, S, pressure, drho_dT, drho_dS, start, npts, EOS)
  real,    intent(in),  dimension(:) ::  T, S, pressure
  real,    intent(out), dimension(:) :: drho_dT, drho_dS
  integer, intent(in)                :: start, npts
  type(EOS_type),        pointer     :: EOS
! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     drho_dT - the partial derivative of density with        *
! *                      potential tempetature, in kg m-3 K-1.         *
! *  (out)     drho_dS - the partial derivative of density with        *
! *                      salinity, in kg m-3 psu-1.                    *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
! *  (in)      EOS - the equation of state type.                       *
! *====================================================================*
! *  This subroutine the appropriate subroutine to calculate density   *
! *  derivatives for an array.                                         *
! *====================================================================*

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_density_derivs called with an unassociated EOS_type EOS.")

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_density_derivs_linear(T, S, pressure, drho_dT, drho_dS, start, &
                                           npts, EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS)
    case (EOS_UNESCO)
      call calculate_density_derivs_unesco(T, S, pressure, drho_dT, drho_dS, start, npts)
    case (EOS_WRIGHT)
      call calculate_density_derivs_wright(T, S, pressure, drho_dT, drho_dS, start, npts)
    case default
      call MOM_error(FATAL, &
           "calculate_density_derivs: EOS%form_of_EOS is not valid.")
  end select

end subroutine calculate_density_derivs

subroutine calculate_compress(T, S, pressure, rho, drho_dp, start, npts, EOS)
  real,    intent(in),  dimension(:) :: T, S, pressure
  real,    intent(out), dimension(:) :: rho, drho_dp
  integer, intent(in)                :: start, npts
  type(EOS_type),       pointer      :: EOS
! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     rho - in situ density in kg m-3.                        *
! *  (out)     drho_dp - the partial derivative of density with        *
! *                      pressure (also the inverse of the square of   *
! *                      sound speed) in s2 m-2.                       *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
! *  (in)      EOS - the equation of state type.                       *
! *====================================================================*
! *  This subroutine computes the in situ density of sea water (rho in *
! *  units of kg/m^3) and the compressibility (drho/dp = C_sound^-2)   *
! *  (drho_dp in units of s2 m-2) from salinity (sal in psu), potential*
! *  temperature (T in deg C), and pressure in Pa.  It uses the expres-*
! *  sions from Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.     *
! *  Coded by R. Hallberg, 1/01                                        *
! *====================================================================*
! *  This subroutine the appropriate subroutine to calculate           *
! *  compressibility for an array.                                     *
! *====================================================================*

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_compress called with an unassociated EOS_type EOS.")

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_compress_linear(T, S, pressure, rho, drho_dp, start, npts, &
                                     EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS)
    case (EOS_UNESCO)
      call calculate_compress_unesco(T, S, pressure, rho, drho_dp, start, npts)
    case (EOS_WRIGHT)
      call calculate_compress_wright(T, S, pressure, rho, drho_dp, start, npts)
    case default
      call MOM_error(FATAL, &
           "calculate_compress: EOS%form_of_EOS is not valid.")
  end select

end subroutine calculate_compress

subroutine calculate_2_densities( T, S, pressure1, pressure2, rho1, rho2, start, npts, EOS)
  real,    intent(in),  dimension(:) :: T, S
  real,    intent(in)                :: pressure1, pressure2
  real,    intent(out), dimension(:) :: rho1, rho2
  integer, intent(in)                :: start, npts
  type(EOS_type),        pointer     :: EOS
! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
! *  (in)      pressure1 - the first pressure in Pa.                   *
! *  (in)      pressure2 -  the second pressure in Pa.                 *
! *  (out)     rho1 - density at pressure1 in kg m-3.                  *
! *  (out)     rho2 - density at pressure2 in kg m-3.                  *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
! *  (in)      EOS - the equation of state type.                       *
! *====================================================================*
! *  This subroutine the appropriate subroutine to calculate density   *
! *  for two arrays.                                                   *
! *====================================================================*

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_2_densities called with an unassociated EOS_type EOS.")

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_2_densities_linear(T, S, pressure1, pressure2, rho1, rho2, start, &
                                        npts, EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS)
    case (EOS_UNESCO)
      call calculate_2_densities_unesco(T, S, pressure1, pressure2, rho1, rho2, start, npts)
    case (EOS_WRIGHT)
      call calculate_2_densities_wright(T, S, pressure1, pressure2, rho1, rho2, start, npts)
    case default
      call MOM_error(FATAL, &
           "calculate_2_densities: EOS%form_of_EOS is not valid.")
  end select

end subroutine calculate_2_densities

subroutine int_specific_vol_dp(T, S, p_t, p_b, alpha_ref, G, EOS, &
                               dza, intp_dza, intx_dza, inty_dza, halo_size)
  real, dimension(NIMEM_,NJMEM_),  intent(in)  :: T, S, p_t, p_b
  real,                            intent(in)  :: alpha_ref
  type(ocean_grid_type),           intent(in)  :: G
  type(EOS_type),                  pointer     :: EOS
  real, dimension(NIMEM_,NJMEM_),  intent(out) :: dza
  real, dimension(NIMEM_,NJMEM_),  optional, intent(out) :: intp_dza
  real, dimension(NIMEMB_,NJMEM_), optional, intent(out) :: intx_dza
  real, dimension(NIMEM_,NJMEMB_), optional, intent(out) :: inty_dza
  integer,                         optional, intent(in)  :: halo_size
!   This subroutine calculates analytical and nearly-analytical integrals in
! pressure across layers of geopotential anomalies, which are required for
! calculating the finite-volume form pressure accelerations in a non-Boussinesq
! model.  There are essentially no free assumptions, apart from the use of
! Bode's rule to do the horizontal integrals, and from a truncation in the
! series for log(1-eps/1+eps) that assumes that |eps| <  .
!
! Arguments: T - potential temperature relative to the surface in C. 
!  (in)      S - salinity in PSU.                                    
!  (in)      p_t - pressure at the top of the layer in Pa.           
!  (in)      p_b - pressure at the top of the layer in Pa.           
!  (in)      alpha_ref - A mean specific volume that is subtracted out to reduce
!                        the magnitude of each of the integrals, m3 kg-1.
!                        The calculation is mathematically identical with
!                        different values of alpha_ref, but this reduces the
!                        effects of roundoff.       
!  (in)      G - The ocean's grid structure.
!  (in)      EOS - The equation of state type.
!  (out)     dza - The change in the geopotential anomaly across the layer, 
!                  in m2 s-2.                                  
!  (out,opt) intp_dza - The integral in pressure through the layer of the 
!                       geopotential anomaly relative to the anomaly at the
!                       bottom of the layer, in Pa m2 s-2. 
!  (out,opt) intx_dza - The integral in x of the difference between the
!                       geopotential anomaly at the top and bottom of the layer
!                       divided by the x grid spacing, in m2 s-2.
!  (out,opt) inty_dza - The integral in y of the difference between the
!                       geopotential anomaly at the top and bottom of the layer
!                       divided by the y grid spacing, in m2 s-2.
!  (in,opt)  halo_size - The width of halo points on which to calculate dza.

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "int_specific_vol_dp called with an unassociated EOS_type EOS.")

  if (EOS%EOS_quadrature) then
    call int_spec_vol_dp_generic(T, S, p_t, p_b, alpha_ref, G, EOS, &
                                 dza, intp_dza, intx_dza, inty_dza, halo_size)
  else ; select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call int_spec_vol_dp_linear(T, S, p_t, p_b, alpha_ref, G, EOS%Rho_T0_S0, &
                                  EOS%dRho_dT, EOS%dRho_dS, dza, intp_dza, &
                                  intx_dza, inty_dza, halo_size)
    case (EOS_WRIGHT)
      call int_spec_vol_dp_wright(T, S, p_t, p_b, alpha_ref, G, dza, &
                                  intp_dza, intx_dza, inty_dza, halo_size)
    case default
      call int_spec_vol_dp_generic(T, S, p_t, p_b, alpha_ref, G, EOS, &
                                   dza, intp_dza, intx_dza, inty_dza, halo_size)
  end select ; endif

end subroutine int_specific_vol_dp

subroutine int_density_dz(T, S, z_t, z_b, rho_ref, rho_0, G_e, G, EOS, &
                          dpa, intz_dpa, intx_dpa, inty_dpa)
  real, dimension(NIMEM_,NJMEM_),  intent(in)  :: T, S, z_t, z_b
  real,                            intent(in)  :: rho_ref, rho_0, G_e
  type(ocean_grid_type),           intent(in)  :: G
  type(EOS_type),                  pointer     :: EOS
  real, dimension(NIMEM_,NJMEM_),  intent(out) :: dpa
  real, dimension(NIMEM_,NJMEM_),  optional, intent(out) :: intz_dpa
  real, dimension(NIMEMB_,NJMEM_), optional, intent(out) :: intx_dpa
  real, dimension(NIMEM_,NJMEMB_), optional, intent(out) :: inty_dpa
!   This subroutine calculates analytical and nearly-analytical integrals of
! pressure anomalies across layers, which are required for calculating the
! finite-volume form pressure accelerations in a Boussinesq model.  The one
! potentially dodgy assumtion here is that rho_0 is used both in the denominator
! of the accelerations, and in the pressure used to calculated density (the
! latter being -z*rho_0*G_e).  These two uses could be separated if need be.
!
! Arguments: T - potential temperature relative to the surface in C. 
!  (in)      S - salinity in PSU.                                    
!  (in)      z_t - height at the top of the layer in m.           
!  (in)      z_b - height at the top of the layer in m.           
!  (in)      rho_ref - A mean density, in kg m-3, that is subtracted out to reduce
!                    the magnitude of each of the integrals.
!                    (The pressure is calucated as p~=-z*rho_0*G_e.)
!  (in)      rho_0 - A density, in kg m-3, that is used to calculate the pressure
!                    (as p~=-z*rho_0*G_e) used in the equation of state.
!  (in)      G_e - The Earth's gravitational acceleration, in m s-2. 
!  (in)      G - The ocean's grid structure.
!  (in)      EOS - The equation of state type.
!  (out)     dpa - The change in the pressure anomaly across the layer, 
!                  in Pa.                                  
!  (out,opt) intz_dpa - The integral through the thickness of the layer of the
!                       pressure anomaly relative to the anomaly at the top of
!                       the layer, in Pa m. 
!  (out,opt) intx_dpa - The integral in x of the difference between the
!                       pressure anomaly at the top and bottom of the layer
!                       divided by the x grid spacing, in Pa.
!  (out,opt) inty_dpa - The integral in y of the difference between the
!                       pressure anomaly at the top and bottom of the layer
!                       divided by the y grid spacing, in Pa.

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "int_density_dz called with an unassociated EOS_type EOS.")

  if (EOS%EOS_quadrature) then
    call int_density_dz_generic(T, S, z_t, z_b, rho_ref, rho_0, G_e, G, &
                                EOS, dpa, intz_dpa, intx_dpa, inty_dpa)
  else ; select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call int_density_dz_linear(T, S, z_t, z_b, rho_ref, rho_0, G_e, G, &
                                 EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS, &
                                 dpa, intz_dpa, intx_dpa, inty_dpa)
    case (EOS_WRIGHT)
      call int_density_dz_wright(T, S, z_t, z_b, rho_ref, rho_0, G_e, G, &
                                 dpa, intz_dpa, intx_dpa, inty_dpa)
    case default
      call int_density_dz_generic(T, S, z_t, z_b, rho_ref, rho_0, G_e, G, &
                                  EOS, dpa, intz_dpa, intx_dpa, inty_dpa)
  end select ; endif

end subroutine int_density_dz

function query_compressible(EOS)
  type(EOS_type),        pointer    :: EOS
  logical :: query_compressible
! Argument:  EOS - the equation of state type.
!
!  This function indicates whether an equation of state with nonzero
! compressibility (i.e., drho/dp) is being used.

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "query_compressible called with an unassociated EOS_type EOS.")

  query_compressible = EOS%compressible
end function query_compressible

subroutine select_eqn_of_state(param_file, EOS)
  type(param_file_type), intent(in) :: param_file
  type(EOS_type),        pointer    :: EOS
! *====================================================================*
! *  (in)      param_file  - parameter file                            *
! *  (out)     EOS - equation of state type                            *
! *====================================================================*
! *  This subroutine reads the EQN_OF_STATE parameter from param_file  *
! *  and sets EOS%form_of_EOS to the appropriate integer that selects  *
! *  the equation of state in the other subroutines in this module.    *
! *  in the case of a linear equation of state, it also sets the       *
! *  run-time parseable parameters of the equation of state.           *
! *====================================================================*
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'
  character(len=40)  :: mod = "MOM_EOS" ! This module's name.
  character(len=40)  :: tmpstr

  if (.not.associated(EOS)) allocate(EOS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, tagname)

  call get_param(param_file, mod, "EQN_OF_STATE", tmpstr, &
                 "EQN_OF_STATE determines which ocean equation of state \n"//&
                 "should be used.  Currently, the valid choices are \n"//&
                 '"LINEAR", "UNESCO", and "WRIGHT". \n'//&
                 "This is only used if USE_EOS is true.", default=EOS_DEFAULT)
  select case (uppercase(tmpstr))
    case (EOS_LINEAR_STRING)
      EOS%form_of_EOS = EOS_LINEAR
    case (EOS_UNESCO_STRING)
      EOS%form_of_EOS = EOS_UNESCO
    case (EOS_WRIGHT_STRING)
      EOS%form_of_EOS = EOS_WRIGHT
    case default
      call MOM_error(FATAL, "interpret_eos_selection: EQN_OF_STATE "//&
                              trim(tmpstr) // "in input file is invalid.")
  end select
  call MOM_mesg('interpret_eos_selection: equation of state set to "' // &
                trim(tmpstr)//'"', 5)

  if (EOS%form_of_EOS == EOS_LINEAR) then
    EOS%Compressible = .false.
    call get_param(param_file, mod, "RHO_T0_S0", EOS%Rho_T0_S0, &
                 "When EQN_OF_STATE="//trim(EOS_LINEAR_STRING)//", \n"//&
                 "this is the density at T=0, S=0.", units="kg m-3", &
                 default=1000.0)
    call get_param(param_file, mod, "DRHO_DT", EOS%dRho_dT, &
                 "When EQN_OF_STATE="//trim(EOS_LINEAR_STRING)//", \n"//&
                 "this is the partial derivative of density with \n"//&
                 "temperature.", units="kg m-3 K-1", default=-0.2)
    call get_param(param_file, mod, "DRHO_DS", EOS%dRho_dS, &
                 "When EQN_OF_STATE="//trim(EOS_LINEAR_STRING)//", \n"//&
                 "this is the partial derivative of density with \n"//&
                 "salinity.", units="kg m-3 PSU-1", default=0.8)
  endif

  call get_param(param_file, mod, "EOS_QUADRATURE", EOS%EOS_quadrature, &
                 "If true, always use the generic (quadrature) code \n"//&
                 "code for the integrals of density.", default=.false.) 

  call get_param(param_file, mod, "TFREEZE_FORM", tmpstr, &
                 "TFREEZE_FORM determines which expression should be \n"//&
                 "used for the freezing point.  Currently, the valid \n"//&
                 'choices are "LINEAR", "MILLERO_78".', &
                 default=TFREEZE_DEFAULT)
  select case (uppercase(tmpstr))
    case (TFREEZE_LINEAR_STRING)
      EOS%form_of_TFreeze = TFREEZE_LINEAR
    case (TFREEZE_MILLERO_STRING)
      EOS%form_of_TFreeze = TFREEZE_MILLERO
    case default
      call MOM_error(FATAL, "interpret_eos_selection:  TFREEZE_FORM "//&
                              trim(tmpstr) // "in input file is invalid.")
  end select
  if (EOS%form_of_TFreeze == TFREEZE_LINEAR) then
    call get_param(param_file, mod, "TFREEZE_S0_P0",EOS%TFr_S0_P0, &
                 "When TFREEZE_FORM="//trim(TFREEZE_LINEAR_STRING)//", \n"//&
                 "this is the freezing potential temperature at \n"//&
                 "S=0, P=0.", units="deg C", default=0.0)
    call get_param(param_file, mod, "DTFREEZE_DS",EOS%dTFr_dS, &
                 "When TFREEZE_FORM="//trim(TFREEZE_LINEAR_STRING)//", \n"//&
                 "this is the derivative of the freezing potential \n"//&
                 "temperature with salinity.", &
                 units="deg C PSU-1", default=-0.054)
    call get_param(param_file, mod, "DTFREEZE_DP",EOS%dTFr_dP, &
                 "When TFREEZE_FORM="//trim(TFREEZE_LINEAR_STRING)//", \n"//&
                 "this is the derivative of the freezing potential \n"//&
                 "temperature with pressure.", &
                 units="deg C Pa-1", default=0.0)
  endif

end subroutine select_eqn_of_state

subroutine deselect_eqn_of_state(EOS)
  type(EOS_type), pointer :: EOS
! *====================================================================*
! *  (in/out)     EOS - equation of state identifier                   *
! *====================================================================*
! *  This subroutine deallocates EOS.                                  *
! *====================================================================*
  
  if (associated(EOS)) deallocate(EOS)
end subroutine deselect_eqn_of_state


subroutine int_density_dz_generic(T, S, z_t, z_b, rho_ref, rho_0, G_e, G, &
                                  EOS, dpa, intz_dpa, intx_dpa, inty_dpa)
  real, dimension(NIMEM_,NJMEM_),  intent(in)  :: T, S, z_t, z_b
  real,                            intent(in)  :: rho_ref, rho_0, G_e
  type(ocean_grid_type),           intent(in)  :: G
  type(EOS_type),                  pointer     :: EOS
  real, dimension(NIMEM_,NJMEM_),  intent(out) :: dpa
  real, dimension(NIMEM_,NJMEM_),  optional, intent(out) :: intz_dpa
  real, dimension(NIMEMB_,NJMEM_), optional, intent(out) :: intx_dpa
  real, dimension(NIMEM_,NJMEMB_), optional, intent(out) :: inty_dpa
!   This subroutine calculates (by numerical quadrature) integrals of
! pressure anomalies across layers, which are required for calculating the
! finite-volume form pressure accelerations in a Boussinesq model.  The one
! potentially dodgy assumtion here is that rho_0 is used both in the denominator
! of the accelerations, and in the pressure used to calculated density (the
! latter being -z*rho_0*G_e).  These two uses could be separated if need be.
!
! Arguments: T - potential temperature relative to the surface in C. 
!  (in)      S - salinity in PSU.                                    
!  (in)      z_t - height at the top of the layer in m.           
!  (in)      z_b - height at the top of the layer in m.           
!  (in)      rho_ref - A mean density, in kg m-3, that is subtracted out to reduce
!                    the magnitude of each of the integrals.
!                    (The pressure is calucated as p~=-z*rho_0*G_e.)
!  (in)      rho_0 - A density, in kg m-3, that is used to calculate the pressure
!                    (as p~=-z*rho_0*G_e) used in the equation of state.
!  (in)      G_e - The Earth's gravitational acceleration, in m s-2. 
!  (in)      G - The ocean's grid structure.
!  (in)      EOS - type that selects the eqn of state.
!  (out)     dpa - The change in the pressure anomaly across the layer, 
!                  in Pa.                                  
!  (out,opt) intz_dpa - The integral through the thickness of the layer of the
!                       pressure anomaly relative to the anomaly at the top of
!                       the layer, in Pa m. 
!  (out,opt) intx_dpa - The integral in x of the difference between the
!                       pressure anomaly at the top and bottom of the layer
!                       divided by the x grid spacing, in Pa.
!  (out,opt) inty_dpa - The integral in y of the difference between the
!                       pressure anomaly at the top and bottom of the layer
!                       divided by the y grid spacing, in Pa.

  real :: T5(5), S5(5), p5(5), r5(5)
  real :: rho_anom
  real :: w_left, w_right, intz(5)
  real, parameter :: C1_90 = 1.0/90.0  ! Rational constants.
  real :: GxRho, I_Rho
  real :: dz
  integer :: Isq, Ieq, Jsq, Jeq, i, j, m, n

  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq

  GxRho = G_e * rho_0
  I_Rho = 1.0 / rho_0

  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    dz = z_t(i,j) - z_b(i,j)
    do n=1,5
      T5(n) = T(i,j) ; S5(n) = S(i,j)
      p5(n) = -GxRho*(z_t(i,j) - 0.25*real(n-1)*dz)
    enddo
    call calculate_density(T5, S5, p5, r5, 1, 5, EOS)
    
    ! Use Bode's rule to estimate the pressure anomaly change.
    rho_anom = C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + 12.0*r5(3)) - &
               rho_ref
    dpa(i,j) = G_e*dz*rho_anom
    ! Use a Bode's-rule-like fifth-order accurate estimate of the double integral of
    ! the pressure anomaly.
    if (present(intz_dpa)) intz_dpa(i,j) = 0.5*G_e*dz**2 * &
          (rho_anom - C1_90*(16.0*(r5(4)-r5(2)) + 7.0*(r5(5)-r5(1))) )
  enddo ; enddo

  if (present(intx_dpa)) then ; do j=G%jsc,G%jec ; do I=Isq,Ieq
    intz(1) = dpa(i,j) ; intz(5) = dpa(i+1,j)
    do m=2,4
      w_left = 0.25*real(5-m) ; w_right = 1.0-w_left
      dz = w_left*(z_t(i,j) - z_b(i,j)) + w_right*(z_t(i+1,j) - z_b(i+1,j))
      T5(1) = w_left*T(i,j) + w_right*T(i+1,j)
      S5(1) = w_left*S(i,j) + w_right*S(i+1,j)
      p5(1) = -GxRho*(w_left*z_t(i,j) + w_right*z_t(i+1,j))
      do n=2,5
        T5(n) = T5(1) ; S5(n) = S5(1)
        p5(n) = p5(n-1) + GxRho*0.25*dz
      enddo
      call calculate_density(T5, S5, p5, r5, 1, 5, EOS)
    
    ! Use Bode's rule to estimate the pressure anomaly change.
      intz(m) = G_e*dz*( C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + &
                                12.0*r5(3)) - rho_ref)
    enddo
    ! Use Bode's rule to integrate the bottom pressure anomaly values in x.
    intx_dpa(i,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                           12.0*intz(3))
  enddo ; enddo ; endif

  if (present(inty_dpa)) then ; do J=Jsq,Jeq ; do i=G%isc,G%iec
    intz(1) = dpa(i,j) ; intz(5) = dpa(i,j+1)
    do m=2,4
      w_left = 0.25*real(5-m) ; w_right = 1.0-w_left
      dz = w_left*(z_t(i,j) - z_b(i,j)) + w_right*(z_t(i,j+1) - z_b(i,j+1))
      T5(1) = w_left*T(i,j) + w_right*T(i,j+1)
      S5(1) = w_left*S(i,j) + w_right*S(i,j+1)
      p5(1) = -GxRho*(w_left*z_t(i,j) + w_right*z_t(i,j+1))
      do n=2,5
        T5(n) = T5(1) ; S5(n) = S5(1)
        p5(n) = p5(n-1) + GxRho*0.25*dz
      enddo
      call calculate_density(T5, S5, p5, r5, 1, 5, EOS)
    
    ! Use Bode's rule to estimate the pressure anomaly change.
      intz(m) = G_e*dz*( C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + &
                                12.0*r5(3)) - rho_ref)
    enddo
    ! Use Bode's rule to integrate the values.
    inty_dpa(i,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                           12.0*intz(3))
  enddo ; enddo ; endif
end subroutine int_density_dz_generic

subroutine int_spec_vol_dp_generic(T, S, p_t, p_b, alpha_ref, G, EOS, &
                                   dza, intp_dza, intx_dza, inty_dza, halo_size)
  real, dimension(NIMEM_,NJMEM_),  intent(in)  :: T, S, p_t, p_b
  real,                            intent(in)  :: alpha_ref
  type(ocean_grid_type),           intent(in)  :: G
  type(EOS_type),                  pointer     :: EOS
  real, dimension(NIMEM_,NJMEM_),  intent(out) :: dza
  real, dimension(NIMEM_,NJMEM_),  optional, intent(out) :: intp_dza
  real, dimension(NIMEMB_,NJMEM_), optional, intent(out) :: intx_dza
  real, dimension(NIMEM_,NJMEMB_), optional, intent(out) :: inty_dza
  integer,                         optional, intent(in)  :: halo_size
!   This subroutine calculates analytical and nearly-analytical integrals in
! pressure across layers of geopotential anomalies, which are required for
! calculating the finite-volume form pressure accelerations in a non-Boussinesq
! model.  There are essentially no free assumptions, apart from the use of
! Bode's rule to do the horizontal integrals, and from a truncation in the
! series for log(1-eps/1+eps) that assumes that |eps| < 0.34.
!
! Arguments: T - potential temperature relative to the surface in C. 
!  (in)      S - salinity in PSU.                                    
!  (in)      p_t - pressure at the top of the layer in Pa.           
!  (in)      p_b - pressure at the top of the layer in Pa.           
!  (in)      alpha_ref - A mean specific volume that is subtracted out to reduce
!                        the magnitude of each of the integrals, m3 kg-1.
!                        The calculation is mathematically identical with
!                        different values of alpha_ref, but this reduces the
!                        effects of roundoff.       
!  (in)      G - The ocean's grid structure.
!  (in)      EOS - type that selects the eqn of state.
!  (out)     dza - The change in the geopotential anomaly across the layer, 
!                  in m2 s-2.                                  
!  (out,opt) intp_dza - The integral in pressure through the layer of the 
!                       geopotential anomaly relative to the anomaly at the
!                       bottom of the layer, in Pa m2 s-2. 
!  (out,opt) intx_dza - The integral in x of the difference between the
!                       geopotential anomaly at the top and bottom of the layer
!                       divided by the x grid spacing, in m2 s-2.
!  (out,opt) inty_dza - The integral in y of the difference between the
!                       geopotential anomaly at the top and bottom of the layer
!                       divided by the y grid spacing, in m2 s-2.
!  (in,opt)  halo_size - The width of halo points on which to calculate dza.
  real :: T5(5), S5(5), p5(5), r5(5), a5(5)
  real :: alpha_anom
  real :: w_left, w_right, intp(5)
  real, parameter :: C1_90 = 1.0/90.0  ! Rational constants.
  real :: dp         ! The pressure change through each layer, in Pa.
  integer :: Isq, Ieq, Jsq, Jeq, ish, ieh, jsh, jeh, i, j, m, n, halo

  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq
  halo = 0 ; if (present(halo_size)) halo = MAX(halo_size,0)
  ish = G%isc-halo ; ieh = G%iec+halo ; jsh = G%jsc-halo ; jeh = G%jec+halo
  if (present(intx_dza)) then ; ish = MIN(Isq,ish) ; ieh = MAX(Ieq+1,ieh); endif
  if (present(inty_dza)) then ; jsh = MIN(Jsq,jsh) ; jeh = MAX(Jeq+1,jeh); endif

  do j=jsh,jeh ; do i=ish,ieh
    dp = p_b(i,j) - p_t(i,j)
    do n=1,5
      T5(n) = T(i,j) ; S5(n) = S(i,j)
      p5(n) = p_b(i,j) - 0.25*real(n-1)*dp
    enddo
    call calculate_density(T5, S5, p5, r5, 1, 5, EOS)
    do n=1,5 ; a5(n) = 1.0 / r5(n) ; enddo
    
    ! Use Bode's rule to estimate the pressure anomaly change.
    alpha_anom = C1_90*(7.0*(a5(1)+a5(5)) + 32.0*(a5(2)+a5(4)) + 12.0*a5(3)) - &
                 alpha_ref
    dza(i,j) = dp*alpha_anom
    ! Use a Bode's-rule-like fifth-order accurate estimate of the double integral of
    ! the pressure anomaly.
    if (present(intp_dza)) intp_dza(i,j) = 0.5*dp**2 * &
          (alpha_anom - C1_90*(16.0*(a5(4)-a5(2)) + 7.0*(a5(5)-a5(1))) )
  enddo ; enddo

  if (present(intx_dza)) then ; do j=G%jsc,G%jec ; do I=Isq,Ieq
    intp(1) = dza(i,j) ; intp(5) = dza(i+1,j)
    do m=2,4
      w_left = 0.25*real(5-m) ; w_right = 1.0-w_left
      dp = w_left*(p_b(i,j) - p_t(i,j)) + w_right*(p_b(i+1,j) - p_t(i+1,j))
      T5(1) = w_left*T(i,j) + w_right*T(i+1,j)
      S5(1) = w_left*S(i,j) + w_right*S(i+1,j)
      p5(1) = w_left*p_b(i,j) + w_right*p_b(i+1,j)
      do n=2,5
        T5(n) = T5(1) ; S5(n) = S5(1) ; p5(n) = p5(n-1) - 0.25*dp
      enddo
      call calculate_density(T5, S5, p5, r5, 1, 5, EOS)
      do n=1,5 ; a5(n) = 1.0 / r5(n) ; enddo
    
    ! Use Bode's rule to estimate the pressure anomaly change.
      intp(m) = dp*( C1_90*(7.0*(a5(1)+a5(5)) + 32.0*(a5(2)+a5(4)) + &
                                12.0*a5(3)) - alpha_ref)
    enddo
    ! Use Bode's rule to integrate the bottom pressure anomaly values in x.
    intx_dza(i,j) = C1_90*(7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4)) + &
                           12.0*intp(3))
  enddo ; enddo ; endif

  if (present(inty_dza)) then ; do J=Jsq,Jeq ; do i=G%isc,G%iec
    intp(1) = dza(i,j) ; intp(5) = dza(i,j+1)
    do m=2,4
      w_left = 0.25*real(5-m) ; w_right = 1.0-w_left
      dp = w_left*(p_b(i,j) - p_t(i,j)) + w_right*(p_b(i,j+1) - p_t(i,j+1))
      T5(1) = w_left*T(i,j) + w_right*T(i,j+1)
      S5(1) = w_left*S(i,j) + w_right*S(i,j+1)
      p5(1) = w_left*p_b(i,j) + w_right*p_b(i,j+1)
      do n=2,5
        T5(n) = T5(1) ; S5(n) = S5(1) ; p5(n) = p5(n-1) - 0.25*dp
      enddo
      call calculate_density(T5, S5, p5, r5, 1, 5, EOS)
      do n=1,5 ; a5(n) = 1.0 / r5(n) ; enddo
    
    ! Use Bode's rule to estimate the pressure anomaly change.
      intp(m) = dp*( C1_90*(7.0*(a5(1)+a5(5)) + 32.0*(a5(2)+a5(4)) + &
                                12.0*a5(3)) - alpha_ref)
    enddo
    ! Use Bode's rule to integrate the bottom pressure anomaly values in y.
    inty_dza(i,j) = C1_90*(7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4)) + &
                           12.0*intp(3))
  enddo ; enddo ; endif

end subroutine int_spec_vol_dp_generic

end module MOM_EOS
