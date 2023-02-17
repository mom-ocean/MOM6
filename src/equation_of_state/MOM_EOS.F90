!> Provides subroutines for quantities specific to the equation of state
module MOM_EOS

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_EOS_linear, only : calculate_density_linear, calculate_spec_vol_linear
use MOM_EOS_linear, only : calculate_density_derivs_linear
use MOM_EOS_linear, only : calculate_specvol_derivs_linear, int_density_dz_linear
use MOM_EOS_linear, only : calculate_density_second_derivs_linear
use MOM_EOS_linear, only : calculate_compress_linear, int_spec_vol_dp_linear
use MOM_EOS_Wright, only : calculate_density_wright, calculate_spec_vol_wright
use MOM_EOS_Wright, only : calculate_density_derivs_wright
use MOM_EOS_Wright, only : calculate_specvol_derivs_wright, int_density_dz_wright
use MOM_EOS_Wright, only : calculate_compress_wright, int_spec_vol_dp_wright
use MOM_EOS_Wright, only : calculate_density_second_derivs_wright
use MOM_EOS_UNESCO, only : calculate_density_unesco, calculate_spec_vol_unesco
use MOM_EOS_UNESCO, only : calculate_density_derivs_unesco, calculate_density_unesco
use MOM_EOS_UNESCO, only : calculate_compress_unesco
use MOM_EOS_NEMO,   only : calculate_density_nemo
use MOM_EOS_NEMO,   only : calculate_density_derivs_nemo, calculate_density_nemo
use MOM_EOS_NEMO,   only : calculate_compress_nemo
use MOM_EOS_TEOS10, only : calculate_density_teos10, calculate_spec_vol_teos10
use MOM_EOS_TEOS10, only : calculate_density_derivs_teos10
use MOM_EOS_TEOS10, only : calculate_specvol_derivs_teos10
use MOM_EOS_TEOS10, only : calculate_density_second_derivs_teos10
use MOM_EOS_TEOS10, only : calculate_compress_teos10
use MOM_EOS_TEOS10, only : gsw_sp_from_sr, gsw_pt_from_ct
use MOM_TFreeze,    only : calculate_TFreeze_linear, calculate_TFreeze_Millero
use MOM_TFreeze,    only : calculate_TFreeze_teos10
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_hor_index,   only : hor_index_type
use MOM_string_functions, only : uppercase
use MOM_unit_scaling, only : unit_scale_type

implicit none ; private

#include <MOM_memory.h>

public EOS_domain
public EOS_init
public EOS_manual_init
public EOS_quadrature
public EOS_use_linear
public analytic_int_density_dz
public analytic_int_specific_vol_dp
public calculate_compress
public calculate_density
public calculate_density_derivs
public calculate_density_second_derivs
public calculate_spec_vol
public calculate_specific_vol_derivs
public calculate_TFreeze
public convert_temp_salt_for_TEOS10
public extract_member_EOS
public cons_temp_to_pot_temp
public abs_saln_to_prac_saln
public gsw_sp_from_sr
public gsw_pt_from_ct
public query_compressible

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Calculates density of sea water from T, S and P
interface calculate_density
  module procedure calculate_density_scalar
  module procedure calculate_density_array
  module procedure calculate_density_1d
  module procedure calculate_stanley_density_scalar
  module procedure calculate_stanley_density_array
  module procedure calculate_stanley_density_1d
end interface calculate_density

!> Calculates specific volume of sea water from T, S and P
interface calculate_spec_vol
  module procedure calc_spec_vol_scalar, calculate_spec_vol_array
  module procedure calc_spec_vol_1d
end interface calculate_spec_vol

!> Calculate the derivatives of density with temperature and salinity from T, S, and P
interface calculate_density_derivs
  module procedure calculate_density_derivs_scalar, calculate_density_derivs_array
  module procedure calculate_density_derivs_1d
end interface calculate_density_derivs

!> Calculate the derivatives of specific volume with temperature and salinity from T, S, and P
interface calculate_specific_vol_derivs
  module procedure calc_spec_vol_derivs_1d, calculate_spec_vol_derivs_array
end interface calculate_specific_vol_derivs

!> Calculates the second derivatives of density with various combinations of temperature,
!! salinity, and pressure from T, S and P
interface calculate_density_second_derivs
  module procedure calculate_density_second_derivs_scalar, calculate_density_second_derivs_1d
end interface calculate_density_second_derivs

!> Calculates the freezing point of sea water from T, S and P
interface calculate_TFreeze
  module procedure calculate_TFreeze_scalar, calculate_TFreeze_1d, calculate_TFreeze_array
end interface calculate_TFreeze

!> Calculates the compressibility of water from T, S, and P
interface calculate_compress
  module procedure calculate_compress_scalar, calculate_compress_1d
end interface calculate_compress

!> A control structure for the equation of state
type, public :: EOS_type ; private
  integer :: form_of_EOS = 0 !< The equation of state to use.
  integer :: form_of_TFreeze = 0 !< The expression for the potential temperature
                             !! of the freezing point.
  logical :: EOS_quadrature  !< If true, always use the generic (quadrature)
                             !! code for the integrals of density.
  logical :: Compressible = .true. !< If true, in situ density is a function of pressure.
! The following parameters are used with the linear equation of state only.
  real :: Rho_T0_S0 !< The density at T=0, S=0 [kg m-3]
  real :: dRho_dT   !< The partial derivative of density with temperature [kg m-3 degC-1]
  real :: dRho_dS   !< The partial derivative of density with salinity [kg m-3 ppt-1]
! The following parameters are use with the linear expression for the freezing
! point only.
  real :: TFr_S0_P0 !< The freezing potential temperature at S=0, P=0 [degC]
  real :: dTFr_dS   !< The derivative of freezing point with salinity [degC ppt-1]
  real :: dTFr_dp   !< The derivative of freezing point with pressure [degC Pa-1]

! Unit conversion factors (normally used for dimensional testing but could also allow for
! change of units of arguments to functions
  real :: m_to_Z = 1.      !< A constant that translates distances in meters to the units of depth [Z m-1 ~> 1]
  real :: kg_m3_to_R = 1.  !< A constant that translates kilograms per meter cubed to the
                           !! units of density [R m3 kg-1 ~> 1]
  real :: R_to_kg_m3 = 1.  !< A constant that translates the units of density to
                           !! kilograms per meter cubed [kg m-3 R-1 ~> 1]
  real :: RL2_T2_to_Pa = 1.!< Convert pressures from R L2 T-2 to Pa [Pa T2 R-1 L-2 ~> 1]
  real :: L_T_to_m_s = 1.  !< Convert lateral velocities from L T-1 to m s-1 [m T s-1 L-1 ~> 1]
  real :: degC_to_C = 1.   !< A constant that translates degrees Celsius to the units of temperature [C degC-1 ~> 1]
  real :: C_to_degC = 1.   !< A constant that translates the units of temperature to degrees Celsius [degC C-1 ~> 1]
  real :: ppt_to_S = 1.    !< A constant that translates parts per thousand to the units of salinity [S ppt-1 ~> 1]
  real :: S_to_ppt = 1.    !< A constant that translates the units of salinity to parts per thousand [ppt S-1 ~> 1]

!  logical :: test_EOS = .true. ! If true, test the equation of state
end type EOS_type

! The named integers that might be stored in eqn_of_state_type%form_of_EOS.
integer, parameter, public :: EOS_LINEAR = 1 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_UNESCO = 2 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_WRIGHT = 3 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_TEOS10 = 4 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_NEMO   = 5 !< A named integer specifying an equation of state

character*(10), parameter :: EOS_LINEAR_STRING = "LINEAR" !< A string for specifying the equation of state
character*(10), parameter :: EOS_UNESCO_STRING = "UNESCO" !< A string for specifying the equation of state
character*(10), parameter :: EOS_WRIGHT_STRING = "WRIGHT" !< A string for specifying the equation of state
character*(10), parameter :: EOS_TEOS10_STRING = "TEOS10" !< A string for specifying the equation of state
character*(10), parameter :: EOS_NEMO_STRING   = "NEMO"   !< A string for specifying the equation of state
character*(10), parameter :: EOS_DEFAULT = EOS_WRIGHT_STRING !< The default equation of state

integer, parameter :: TFREEZE_LINEAR = 1  !< A named integer specifying a freezing point expression
integer, parameter :: TFREEZE_MILLERO = 2 !< A named integer specifying a freezing point expression
integer, parameter :: TFREEZE_TEOS10 = 3  !< A named integer specifying a freezing point expression
character*(10), parameter :: TFREEZE_LINEAR_STRING = "LINEAR" !< A string for specifying the freezing point expression
character*(10), parameter :: TFREEZE_MILLERO_STRING = "MILLERO_78" !< A string for specifying
                                                              !! freezing point expression
character*(10), parameter :: TFREEZE_TEOS10_STRING = "TEOS10" !< A string for specifying the freezing point expression
character*(10), parameter :: TFREEZE_DEFAULT = TFREEZE_LINEAR_STRING !< The default freezing point expression

contains

!> Calls the appropriate subroutine to calculate density of sea water for scalar inputs.
!! If rho_ref is present, the anomaly with respect to rho_ref is returned.  The pressure and
!! density can be rescaled with the values stored in EOS.  If the scale argument is present the density
!! scaling uses the product of the two scaling factors.
subroutine calculate_density_scalar(T, S, pressure, rho, EOS, rho_ref, scale)
  real,           intent(in)  :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real,           intent(in)  :: S        !< Salinity [S ~> ppt]
  real,           intent(in)  :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real,           intent(out) :: rho      !< Density (in-situ if pressure is local) [R ~> kg m-3]
  type(EOS_type), intent(in)  :: EOS      !< Equation of state structure
  real, optional, intent(in)  :: rho_ref  !< A reference density [R ~> kg m-3]
  real, optional, intent(in)  :: scale    !< A multiplicative factor by which to scale output density in
                                          !! combination with scaling stored in EOS [various]

  real :: Ta(1)      ! An array of temperatures [degC]
  real :: Sa(1)      ! An array of salinities [ppt]
  real :: pres(1)    ! An mks version of the pressure to use [Pa]
  real :: rho_mks(1) ! An mks version of the density to be returned [kg m-3]
  real :: rho_scale  ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]

  pres(1) = EOS%RL2_T2_to_Pa * pressure
  Ta(1) = EOS%C_to_degC * T ; Sa(1) = EOS%S_to_ppt * S
  if (present(rho_ref)) then
    call calculate_density_array(Ta, Sa, pres, rho_mks, 1, 1, EOS, EOS%R_to_kg_m3*rho_ref)
  else
    call calculate_density_array(Ta, Sa, pres, rho_mks, 1, 1, EOS)
  endif

  ! Rescale the output density to the desired units.
  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  rho = rho_scale * rho_mks(1)

end subroutine calculate_density_scalar

!> Calls the appropriate subroutine to calculate density of sea water for scalar inputs
!! including the variance of T, S and covariance of T-S.
!! The calculation uses only the second order correction in a series as discussed
!! in Stanley et al., 2020.
!! If rho_ref is present, the anomaly with respect to rho_ref is returned. The
!! density can be rescaled using rho_ref.
subroutine calculate_stanley_density_scalar(T, S, pressure, Tvar, TScov, Svar, rho, EOS, rho_ref, scale)
  real,           intent(in)  :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real,           intent(in)  :: S        !< Salinity [S ~> ppt]
  real,           intent(in)  :: Tvar     !< Variance of potential temperature referenced to the surface [C2 ~> degC2]
  real,           intent(in)  :: TScov    !< Covariance of potential temperature and salinity [C S ~> degC ppt]
  real,           intent(in)  :: Svar     !< Variance of salinity [S2 ~> ppt2]
  real,           intent(in)  :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real,           intent(out) :: rho      !< Density (in-situ if pressure is local) [R ~> kg m-3]
  type(EOS_type), intent(in)  :: EOS      !< Equation of state structure
  real, optional, intent(in)  :: rho_ref  !< A reference density [R ~> kg m-3].
  real, optional, intent(in)  :: scale    !< A multiplicative factor by which to scale output density in
                                          !! combination with scaling stored in EOS [various]
  ! Local variables
  real :: d2RdTT   ! Second derivative of density with temperature [kg m-3 degC-2]
  real :: d2RdST   ! Second derivative of density with temperature and salinity [kg m-3 degC-1 ppt-1]
  real :: d2RdSS   ! Second derivative of density with salinity [kg m-3 ppt-2]
  real :: d2RdSp   ! Second derivative of density with salinity and pressure [kg m-3 ppt-1 Pa-1]
  real :: d2RdTp   ! Second derivative of density with temperature and pressure [kg m-3 degC-1 Pa-1]
  real :: p_scale  ! A factor to convert pressure to units of Pa [Pa T2 R-1 L-2 ~> 1]
  real :: T_scale  ! A factor to convert temperature to units of degC [degC C-1 ~> 1]
  real :: S_scale  ! A factor to convert salinity to units of ppt [ppt S-1 ~> 1]

  call calculate_density_scalar(T, S, pressure, rho, EOS, rho_ref)

  p_scale = EOS%RL2_T2_to_Pa
  T_scale = EOS%C_to_degC
  S_scale = EOS%S_to_ppt
  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_density_second_derivs_linear(T_scale*T, S_scale*S, p_scale*pressure, &
                                                  d2RdSS, d2RdST, d2RdTT, d2RdSp, d2RdTP)
    case (EOS_WRIGHT)
      call calculate_density_second_derivs_wright(T_scale*T, S_scale*S, p_scale*pressure, &
                                                  d2RdSS, d2RdST, d2RdTT, d2RdSp, d2RdTP)
    case (EOS_TEOS10)
      call calculate_density_second_derivs_teos10(T_scale*T, S_scale*S, p_scale*pressure, &
                                                  d2RdSS, d2RdST, d2RdTT, d2RdSp, d2RdTP)
    case default
      call MOM_error(FATAL, "calculate_stanley_density_scalar: EOS is not valid.")
  end select

  ! Equation 25 of Stanley et al., 2020.
  rho = rho + EOS%kg_m3_to_R * ( 0.5 * (T_scale**2 * d2RdTT) * Tvar + &
                 ( (S_scale*T_scale * d2RdST) * TScov + 0.5 * (S_scale**2 * d2RdSS) * Svar ) )

  if (present(scale)) rho = rho * scale

end subroutine calculate_stanley_density_scalar

!> Calls the appropriate subroutine to calculate the density of sea water for 1-D array inputs.
!! If rho_ref is present, the anomaly with respect to rho_ref is returned.
subroutine calculate_density_array(T, S, pressure, rho, start, npts, EOS, rho_ref, scale)
  real, dimension(:), intent(in)    :: T        !< Potential temperature referenced to the surface [degC]
  real, dimension(:), intent(in)    :: S        !< Salinity [ppt]
  real, dimension(:), intent(in)    :: pressure !< Pressure [Pa]
  real, dimension(:), intent(inout) :: rho      !< Density (in-situ if pressure is local) [kg m-3] or other
                                                !! units if rescaled via a scale argument
  integer,            intent(in)    :: start    !< Start index for computation
  integer,            intent(in)    :: npts     !< Number of point to compute
  type(EOS_type),     intent(in)    :: EOS      !< Equation of state structure
  real,     optional, intent(in)    :: rho_ref  !< A reference density [kg m-3]
  real,     optional, intent(in)    :: scale    !< A multiplicative factor by which to scale the output
                                                !! density, perhaps to other units than kg m-3 [various]
  integer :: j

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_density_linear(T, S, pressure, rho, start, npts, &
                                    EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS, rho_ref)
    case (EOS_UNESCO)
      call calculate_density_unesco(T, S, pressure, rho, start, npts, rho_ref)
    case (EOS_WRIGHT)
      call calculate_density_wright(T, S, pressure, rho, start, npts, rho_ref)
    case (EOS_TEOS10)
      call calculate_density_teos10(T, S, pressure, rho, start, npts, rho_ref)
    case (EOS_NEMO)
    call calculate_density_nemo(T, S, pressure, rho, start, npts, rho_ref)
    case default
      call MOM_error(FATAL, "calculate_density_array: EOS%form_of_EOS is not valid.")
  end select

  if (present(scale)) then ; if (scale /= 1.0) then ; do j=start,start+npts-1
    rho(j) = scale * rho(j)
  enddo ; endif ; endif

end subroutine calculate_density_array

!> Calls the appropriate subroutine to calculate the density of sea water for 1-D array inputs
!! including the variance of T, S and covariance of T-S.
!! The calculation uses only the second order correction in a series as discussed
!! in Stanley et al., 2020.
!! If rho_ref is present, the anomaly with respect to rho_ref is returned.
subroutine calculate_stanley_density_array(T, S, pressure, Tvar, TScov, Svar, rho, start, npts, EOS, rho_ref, scale)
  real, dimension(:), intent(in)    :: T        !< Potential temperature referenced to the surface [degC]
  real, dimension(:), intent(in)    :: S        !< Salinity [ppt]
  real, dimension(:), intent(in)    :: pressure !< Pressure [Pa]
  real, dimension(:), intent(in)    :: Tvar     !< Variance of potential temperature referenced to the surface [degC2]
  real, dimension(:), intent(in)    :: TScov    !< Covariance of potential temperature and salinity [degC ppt]
  real, dimension(:), intent(in)    :: Svar     !< Variance of salinity [ppt2]
  real, dimension(:), intent(inout) :: rho      !< Density (in-situ if pressure is local) [kg m-3]
  integer,            intent(in)    :: start    !< Start index for computation
  integer,            intent(in)    :: npts     !< Number of point to compute
  type(EOS_type),     intent(in)    :: EOS      !< Equation of state structure
  real,     optional, intent(in)    :: rho_ref  !< A reference density [kg m-3].
  real,     optional, intent(in)    :: scale    !< A multiplicative factor by which to scale the output
                                                !! density, perhaps to other units than kg m-3 [various]
  ! Local variables
  real, dimension(size(T)) :: &
    d2RdTT, &   ! Second derivative of density with temperature [kg m-3 degC-2]
    d2RdST, &   ! Second derivative of density with temperature and salinity [kg m-3 degC-1 ppt-1]
    d2RdSS, &   ! Second derivative of density with salinity [kg m-3 ppt-2]
    d2RdSp, &   ! Second derivative of density with salinity and pressure [kg m-3 ppt-1 Pa-1]
    d2RdTp      ! Second derivative of density with temperature and pressure [kg m-3 degC-1 Pa-1]
  integer :: j

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_density_linear(T, S, pressure, rho, start, npts, &
                                    EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS, rho_ref)
      call calculate_density_second_derivs_linear(T, S, pressure, d2RdSS, d2RdST, &
                                                  d2RdTT, d2RdSp, d2RdTP, start, npts)
    case (EOS_WRIGHT)
      call calculate_density_wright(T, S, pressure, rho, start, npts, rho_ref)
      call calculate_density_second_derivs_wright(T, S, pressure, d2RdSS, d2RdST, &
                                                  d2RdTT, d2RdSp, d2RdTP, start, npts)
    case (EOS_TEOS10)
      call calculate_density_teos10(T, S, pressure, rho, start, npts, rho_ref)
      call calculate_density_second_derivs_teos10(T, S, pressure, d2RdSS, d2RdST, &
                                                  d2RdTT, d2RdSp, d2RdTP, start, npts)
    case default
      call MOM_error(FATAL, "calculate_stanley_density_array: EOS%form_of_EOS is not valid.")
  end select

  ! Equation 25 of Stanley et al., 2020.
  do j=start,start+npts-1
    rho(j) = rho(j) &
             + ( 0.5 * d2RdTT(j) * Tvar(j) + ( d2RdST(j) * TScov(j) + 0.5 * d2RdSS(j) * Svar(j) ) )
  enddo

  if (present(scale)) then ; if (scale /= 1.0) then ; do j=start,start+npts-1
    rho(j) = scale * rho(j)
  enddo ; endif ; endif

end subroutine calculate_stanley_density_array

!> Calls the appropriate subroutine to calculate the density of sea water for 1-D array inputs,
!! potentially limiting the domain of indices that are worked on.
!! If rho_ref is present, the anomaly with respect to rho_ref is returned.
subroutine calculate_density_1d(T, S, pressure, rho, EOS, dom, rho_ref, scale)
  real, dimension(:),    intent(in)    :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(:),    intent(in)    :: S        !< Salinity [S ~> ppt]
  real, dimension(:),    intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:),    intent(inout) :: rho      !< Density (in-situ if pressure is local) [R ~> kg m-3]
  type(EOS_type),        intent(in)    :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.
  real,                  optional, intent(in) :: rho_ref !< A reference density [R ~> kg m-3]
  real,                  optional, intent(in) :: scale !< A multiplicative factor by which to scale density
                                                   !! in combination with scaling stored in EOS [various]
  ! Local variables
  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real, dimension(size(rho)) :: pres  ! Pressure converted to [Pa]
  real, dimension(size(rho)) :: Ta    ! Temperature converted to [degC]
  real, dimension(size(rho)) :: Sa    ! Salinity converted to [ppt]
  integer :: i, is, ie, npts

  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(rho) ; npts = 1 + ie - is
  endif

  if ((EOS%RL2_T2_to_Pa == 1.0) .and. (EOS%R_to_kg_m3 == 1.0) .and. &
      (EOS%C_to_degC == 1.0) .and. (EOS%S_to_ppt == 1.0)) then
    call calculate_density_array(T, S, pressure, rho, is, npts, EOS, rho_ref=rho_ref)
  else ! This is the same as above, but with some extra work to rescale variables.
    do i=is,ie
      pres(i) = EOS%RL2_T2_to_Pa * pressure(i)
      Ta(i) = EOS%C_to_degC * T(i)
      Sa(i) = EOS%S_to_ppt * S(i)
    enddo
    if (present(rho_ref)) then
      call calculate_density_array(Ta, Sa, pres, rho, is, npts, EOS, rho_ref=EOS%R_to_kg_m3*rho_ref)
    else
      call calculate_density_array(Ta, Sa, pres, rho, is, npts, EOS)
    endif
  endif

  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  if (rho_scale /= 1.0) then ; do i=is,ie
    rho(i) = rho_scale * rho(i)
  enddo ; endif

end subroutine calculate_density_1d

!> Calls the appropriate subroutine to calculate the density of sea water for 1-D array inputs
!! including the variance of T, S and covariance of T-S,
!! potentially limiting the domain of indices that are worked on.
!! The calculation uses only the second order correction in a series as discussed
!! in Stanley et al., 2020.
!! If rho_ref is present, the anomaly with respect to rho_ref is returned.
subroutine calculate_stanley_density_1d(T, S, pressure, Tvar, TScov, Svar, rho, EOS, dom, rho_ref, scale)
  real, dimension(:),    intent(in)    :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(:),    intent(in)    :: S        !< Salinity [S ~> ppt]
  real, dimension(:),    intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:),    intent(in)    :: Tvar     !< Variance of potential temperature [C2 ~> degC2]
  real, dimension(:),    intent(in)    :: TScov    !< Covariance of potential temperature and salinity [C S ~> degC ppt]
  real, dimension(:),    intent(in)    :: Svar     !< Variance of salinity [S2 ~> ppt2]
  real, dimension(:),    intent(inout) :: rho      !< Density (in-situ if pressure is local) [R ~> kg m-3]
  type(EOS_type),        intent(in)    :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.
  real,                  optional, intent(in) :: rho_ref !< A reference density [R ~> kg m-3]
  real,                  optional, intent(in) :: scale !< A multiplicative factor by which to scale density
                                                   !! in combination with scaling stored in EOS [various]
  ! Local variables
  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real :: T2_scale  ! A factor to convert temperature variance to units of degC2 [degC2 C-2 ~> 1]
  real :: S2_scale  ! A factor to convert salinity variance to units of ppt2 [ppt2 S-2 ~> 1]
  real :: TS_scale  ! A factor to convert temperature-salinity covariance to units of
                    ! degC ppt [degC ppt C-1 S-1 ~> 1]
  real :: rho_reference ! rho_ref converted to [kg m-3]
  real, dimension(size(rho)) :: pres  ! Pressure converted to [Pa]
  real, dimension(size(rho)) :: Ta    ! Temperature converted to [degC]
  real, dimension(size(rho)) :: Sa    ! Salinity converted to [ppt]
  real, dimension(size(T)) :: &
    d2RdTT, &   ! Second derivative of density with temperature [kg m-3 degC-2]
    d2RdST, &   ! Second derivative of density with temperature and salinity [kg m-3 degC-1 ppt-1]
    d2RdSS, &   ! Second derivative of density with salinity [kg m-3 ppt-2]
    d2RdSp, &   ! Second derivative of density with salinity and pressure [kg m-3 ppt-1 Pa-1]
    d2RdTp      ! Second derivative of density with temperature and pressure [kg m-3 degC-1 Pa-1]
  integer :: i, is, ie, npts

  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(rho) ; npts = 1 + ie - is
  endif

  do i=is,ie
    pres(i) = EOS%RL2_T2_to_Pa * pressure(i)
    Ta(i) = EOS%C_to_degC * T(i)
    Sa(i) = EOS%S_to_ppt * S(i)
  enddo
  T2_scale = EOS%C_to_degC**2
  S2_scale = EOS%S_to_ppt**2
  TS_scale = EOS%C_to_degC*EOS%S_to_ppt

  ! Rho_ref is seems like it is always present when calculate_Stanley_density is called, so
  ! always set rho_reference, even though a 0 value can change answers at roundoff with
  ! some equations of state.
  rho_reference = 0.0 ; if (present(rho_ref)) rho_reference = EOS%R_to_kg_m3*rho_ref

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_density_linear(Ta, Sa, pres, rho, is, npts, &
                                    EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS, rho_reference)
      call calculate_density_second_derivs_linear(Ta, Sa, pres, d2RdSS, d2RdST, &
                                                  d2RdTT, d2RdSp, d2RdTP, is, npts)
    case (EOS_WRIGHT)
      call calculate_density_wright(Ta, Sa, pres, rho, is, npts, rho_reference)
      call calculate_density_second_derivs_wright(Ta, Sa, pres, d2RdSS, d2RdST, &
                                                  d2RdTT, d2RdSp, d2RdTP, is, npts)
    case (EOS_TEOS10)
      call calculate_density_teos10(Ta, Sa, pres, rho, is, npts, rho_reference)
      call calculate_density_second_derivs_teos10(Ta, Sa, pres, d2RdSS, d2RdST, &
                                                  d2RdTT, d2RdSp, d2RdTP, is, npts)
    case default
      call MOM_error(FATAL, "calculate_stanley_density_scalar: EOS is not valid.")
  end select

  ! Equation 25 of Stanley et al., 2020.
  do i=is,ie
    rho(i) = rho(i) + ( 0.5 * (T2_scale * d2RdTT(i)) * Tvar(i) + &
                        ( (TS_scale * d2RdST(i)) * TScov(i) + &
                          0.5 * (S2_scale * d2RdSS(i)) * Svar(i) ) )
  enddo

  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  if (rho_scale /= 1.0) then ; do i=is,ie
    rho(i) = rho_scale * rho(i)
  enddo ; endif

end subroutine calculate_stanley_density_1d

!> Calls the appropriate subroutine to calculate the specific volume of sea water
!! for 1-D array inputs.
subroutine calculate_spec_vol_array(T, S, pressure, specvol, start, npts, EOS, spv_ref, scale)
  real, dimension(:), intent(in)    :: T        !< potential temperature relative to the surface [degC]
  real, dimension(:), intent(in)    :: S        !< salinity [ppt]
  real, dimension(:), intent(in)    :: pressure !< pressure [Pa]
  real, dimension(:), intent(inout) :: specvol  !< in situ specific volume [kg m-3]
  integer,            intent(in)    :: start    !< the starting point in the arrays.
  integer,            intent(in)    :: npts     !< the number of values to calculate.
  type(EOS_type),     intent(in)    :: EOS      !< Equation of state structure
  real,     optional, intent(in)    :: spv_ref  !< A reference specific volume [m3 kg-1]
  real,     optional, intent(in)    :: scale    !< A multiplicative factor by which to scale specific
                                                !! volume in combination with scaling stored in EOS [various]

  real, dimension(size(specvol))  :: rho   ! Density [kg m-3]
  integer :: j

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_spec_vol_linear(T, S, pressure, specvol, start, npts, &
               EOS%rho_T0_S0, EOS%drho_dT, EOS%drho_dS, spv_ref)
    case (EOS_UNESCO)
      call calculate_spec_vol_unesco(T, S, pressure, specvol, start, npts, spv_ref)
    case (EOS_WRIGHT)
      call calculate_spec_vol_wright(T, S, pressure, specvol, start, npts, spv_ref)
    case (EOS_TEOS10)
      call calculate_spec_vol_teos10(T, S, pressure, specvol, start, npts, spv_ref)
    case (EOS_NEMO)
      call calculate_density_nemo(T, S, pressure, rho, start, npts)
      if (present(spv_ref)) then
        specvol(:) = 1.0 / rho(:) - spv_ref
      else
        specvol(:) = 1.0 / rho(:)
      endif
    case default
      call MOM_error(FATAL, "calculate_spec_vol_array: EOS%form_of_EOS is not valid.")
  end select

  if (present(scale)) then ; if (scale /= 1.0) then ; do j=start,start+npts-1
    specvol(j) = scale * specvol(j)
  enddo ; endif ; endif

end subroutine calculate_spec_vol_array

!> Calls the appropriate subroutine to calculate specific volume of sea water
!! for scalar inputs.
subroutine calc_spec_vol_scalar(T, S, pressure, specvol, EOS, spv_ref, scale)
  real,           intent(in)  :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real,           intent(in)  :: S        !< Salinity [S ~> ppt]
  real,           intent(in)  :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real,           intent(out) :: specvol  !< In situ or potential specific volume [R-1 ~> m3 kg-1]
                                          !! or other units determined by the scale argument
  type(EOS_type), intent(in)  :: EOS      !< Equation of state structure
  real, optional, intent(in)  :: spv_ref  !< A reference specific volume [R-1 ~> m3 kg-1]
  real, optional, intent(in)  :: scale    !< A multiplicative factor by which to scale specific
                                          !! volume in combination with scaling stored in EOS [various]

  real, dimension(1) :: Ta   ! Rescaled single element array version of temperature [degC]
  real, dimension(1) :: Sa   ! Rescaled single element array version of salinity [ppt]
  real, dimension(1) :: pres ! Rescaled single element array version of pressure [Pa]
  real, dimension(1) :: spv  ! Rescaled single element array version of specific volume [m3 kg-1]
  real :: spv_scale ! A factor to convert specific volume from m3 kg-1 to the desired units [kg R-1 m-3 ~> 1]

  pres(1) = EOS%RL2_T2_to_Pa * pressure
  Ta(1) = EOS%C_to_degC * T ; Sa(1) = EOS%S_to_ppt * S

  if (present(spv_ref)) then
    call calculate_spec_vol_array(Ta, Sa, pres, spv, 1, 1, EOS, EOS%kg_m3_to_R*spv_ref)
  else
    call calculate_spec_vol_array(Ta, Sa, pres, spv, 1, 1, EOS)
  endif
  specvol = spv(1)

  spv_scale = EOS%R_to_kg_m3
  if (present(scale)) spv_scale = spv_scale * scale
  if (spv_scale /= 1.0) then
    specvol = spv_scale * specvol
  endif

end subroutine calc_spec_vol_scalar

!> Calls the appropriate subroutine to calculate the specific volume of sea water for 1-D array
!! inputs, potentially limiting the domain of indices that are worked on.
subroutine calc_spec_vol_1d(T, S, pressure, specvol, EOS, dom, spv_ref, scale)
  real, dimension(:),    intent(in)    :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(:),    intent(in)    :: S        !< Salinity [S ~> ppt]
  real, dimension(:),    intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:),    intent(inout) :: specvol  !< In situ specific volume [R-1 ~> m3 kg-1]
  type(EOS_type),        intent(in)    :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.
  real,                  optional, intent(in) :: spv_ref !< A reference specific volume [R-1 ~> m3 kg-1]
  real,                  optional, intent(in) :: scale !< A multiplicative factor by which to scale
                                                       !! output specific volume in combination with
                                                       !! scaling stored in EOS [various]
  ! Local variables
  real, dimension(size(T)) :: pres  ! Pressure converted to [Pa]
  real, dimension(size(T)) :: Ta    ! Temperature converted to [degC]
  real, dimension(size(T)) :: Sa    ! Salinity converted to [ppt]
  real :: spv_scale ! A factor to convert specific volume from m3 kg-1 to the desired units [kg m-3 R-1 ~> 1]
  integer :: i, is, ie, npts

  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(specvol) ; npts = 1 + ie - is
  endif

  if ((EOS%RL2_T2_to_Pa == 1.0) .and. (EOS%kg_m3_to_R == 1.0) .and. &
      (EOS%C_to_degC == 1.0) .and. (EOS%S_to_ppt == 1.0)) then
    call calculate_spec_vol_array(T, S, pressure, specvol, is, npts, EOS, spv_ref)
  else ! This is the same as above, but with some extra work to rescale variables.
    do i=is,ie
      pres(i) = EOS%RL2_T2_to_Pa * pressure(i)
      Ta(i) = EOS%C_to_degC * T(i)
      Sa(i) = EOS%S_to_ppt * S(i)
    enddo
    if (present(spv_ref)) then
      call calculate_spec_vol_array(Ta, Sa, pres, specvol, is, npts, EOS, EOS%kg_m3_to_R*spv_ref)
    else
      ! There is rescaling of variables, but spv_ref is not present. Passing a 0 value of spv_ref
      ! changes answers at roundoff for some equations of state, like Wright and UNESCO.
      call calculate_spec_vol_array(Ta, Sa, pres, specvol, is, npts, EOS)
    endif
  endif

  spv_scale = EOS%R_to_kg_m3
  if (present(scale)) spv_scale = spv_scale * scale
  if (spv_scale /= 1.0) then ; do i=is,ie
    specvol(i) = spv_scale * specvol(i)
  enddo ; endif

end subroutine calc_spec_vol_1d


!> Calls the appropriate subroutine to calculate the freezing point for scalar inputs.
subroutine calculate_TFreeze_scalar(S, pressure, T_fr, EOS, pres_scale, scale_from_EOS)
  real,           intent(in)  :: S    !< Salinity, [ppt] or [S ~> ppt] depending on scale_from_EOS
  real,           intent(in)  :: pressure !< Pressure, in [Pa] or [R L2 T-2 ~> Pa] depending on
                                      !! pres_scale or scale_from_EOS
  real,           intent(out) :: T_fr !< Freezing point potential temperature referenced to the
                                      !! surface [degC] or [C ~> degC] depending on scale_from_EOS
  type(EOS_type), intent(in)  :: EOS  !< Equation of state structure
  real, optional, intent(in)  :: pres_scale  !< A multiplicative factor to convert pressure
                                      !! into Pa [Pa T2 R-1 L-2 ~> 1].
  logical, optional, intent(in)  :: scale_from_EOS !< If present true use the dimensional scaling
                                      !! factors stored in EOS.  Omission is the same .false.

  ! Local variables
  real :: p_scale ! A factor to convert pressure to units of Pa [Pa T2 R-1 L-2 ~> 1]
  real :: S_scale ! A factor to convert salinity to units of ppt [ppt S-1 ~> 1]

  p_scale = 1.0 ; S_scale = 1.0
  if (present(pres_scale)) p_scale = pres_scale
  if (present(scale_from_EOS)) then ; if (scale_from_EOS) then
    p_scale = EOS%RL2_T2_to_Pa
    S_scale = EOS%S_to_ppt
  endif ; endif

  select case (EOS%form_of_TFreeze)
    case (TFREEZE_LINEAR)
      call calculate_TFreeze_linear(S_scale*S, p_scale*pressure, T_fr, EOS%TFr_S0_P0, &
                                    EOS%dTFr_dS, EOS%dTFr_dp)
    case (TFREEZE_MILLERO)
      call calculate_TFreeze_Millero(S_scale*S, p_scale*pressure, T_fr)
    case (TFREEZE_TEOS10)
      call calculate_TFreeze_teos10(S_scale*S, p_scale*pressure, T_fr)
    case default
      call MOM_error(FATAL, "calculate_TFreeze_scalar: form_of_TFreeze is not valid.")
  end select

  if (present(scale_from_EOS)) then ; if (scale_from_EOS) then
    T_fr = EOS%degC_to_C * T_fr
  endif ; endif

end subroutine calculate_TFreeze_scalar

!> Calls the appropriate subroutine to calculate the freezing point for a 1-D array.
subroutine calculate_TFreeze_array(S, pressure, T_fr, start, npts, EOS, pres_scale)
  real, dimension(:), intent(in)    :: S        !< Salinity [ppt]
  real, dimension(:), intent(in)    :: pressure !< Pressure, in [Pa] or [R L2 T-2 ~> Pa] depending on pres_scale
  real, dimension(:), intent(inout) :: T_fr     !< Freezing point potential temperature referenced
                                                !! to the surface [degC]
  integer,            intent(in)    :: start    !< Starting index within the array
  integer,            intent(in)    :: npts     !< The number of values to calculate
  type(EOS_type),     intent(in)    :: EOS      !< Equation of state structure
  real,     optional, intent(in)    :: pres_scale !< A multiplicative factor to convert pressure
                                                !! into Pa [Pa T2 R-1 L-2 ~> 1].

  ! Local variables
  real, dimension(size(pressure)) :: pres  ! Pressure converted to [Pa]
  real :: p_scale  ! A factor to convert pressure to units of Pa [Pa T2 R-1 L-2 ~> 1]
  integer :: j

  p_scale = 1.0 ; if (present(pres_scale)) p_scale = pres_scale

  if (p_scale == 1.0) then
    select case (EOS%form_of_TFreeze)
      case (TFREEZE_LINEAR)
        call calculate_TFreeze_linear(S, pressure, T_fr, start, npts, &
                                      EOS%TFr_S0_P0, EOS%dTFr_dS, EOS%dTFr_dp)
      case (TFREEZE_MILLERO)
        call calculate_TFreeze_Millero(S, pressure, T_fr, start, npts)
      case (TFREEZE_TEOS10)
        call calculate_TFreeze_teos10(S, pressure, T_fr, start, npts)
      case default
        call MOM_error(FATAL, "calculate_TFreeze_scalar: form_of_TFreeze is not valid.")
    end select
  else
    do j=start,start+npts-1 ; pres(j) = p_scale * pressure(j) ; enddo
    select case (EOS%form_of_TFreeze)
      case (TFREEZE_LINEAR)
        call calculate_TFreeze_linear(S, pres, T_fr, start, npts, &
                                      EOS%TFr_S0_P0, EOS%dTFr_dS, EOS%dTFr_dp)
      case (TFREEZE_MILLERO)
        call calculate_TFreeze_Millero(S, pres, T_fr, start, npts)
      case (TFREEZE_TEOS10)
        call calculate_TFreeze_teos10(S, pres, T_fr, start, npts)
      case default
        call MOM_error(FATAL, "calculate_TFreeze_scalar: form_of_TFreeze is not valid.")
    end select
  endif

end subroutine calculate_TFreeze_array

!> Calls the appropriate subroutine to calculate the freezing point for a 1-D array, taking
!! dimensionally rescaled arguments with factors stored in EOS.
subroutine calculate_TFreeze_1d(S, pressure, T_fr, EOS, dom)
  real, dimension(:), intent(in)    :: S        !< Salinity [S ~> ppt]
  real, dimension(:), intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:), intent(inout) :: T_fr     !< Freezing point potential temperature referenced
                                                !! to the surface [C ~> degC]
  type(EOS_type),     intent(in)    :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.

  ! Local variables
  real, dimension(size(T_fr)) :: pres  ! Pressure converted to [Pa]
  real, dimension(size(T_fr)) :: Sa    ! Salinity converted to [ppt]
  integer :: i, is, ie, npts

  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(T_Fr) ; npts = 1 + ie - is
  endif

  if ((EOS%RL2_T2_to_Pa == 1.0) .and. (EOS%S_to_ppt == 1.0)) then
    select case (EOS%form_of_TFreeze)
      case (TFREEZE_LINEAR)
        call calculate_TFreeze_linear(S, pressure, T_fr, is, npts, &
                                      EOS%TFr_S0_P0, EOS%dTFr_dS, EOS%dTFr_dp)
      case (TFREEZE_MILLERO)
        call calculate_TFreeze_Millero(S, pressure, T_fr, is, npts)
      case (TFREEZE_TEOS10)
        call calculate_TFreeze_teos10(S, pressure, T_fr, is, npts)
      case default
        call MOM_error(FATAL, "calculate_TFreeze_scalar: form_of_TFreeze is not valid.")
    end select
  else
    do i=is,ie
      pres(i) = EOS%RL2_T2_to_Pa * pressure(i)
      Sa(i) = EOS%S_to_ppt * S(i)
    enddo
    select case (EOS%form_of_TFreeze)
      case (TFREEZE_LINEAR)
        call calculate_TFreeze_linear(Sa, pres, T_fr, is, npts, &
                                      EOS%TFr_S0_P0, EOS%dTFr_dS, EOS%dTFr_dp)
      case (TFREEZE_MILLERO)
        call calculate_TFreeze_Millero(Sa, pres, T_fr, is, npts)
      case (TFREEZE_TEOS10)
        call calculate_TFreeze_teos10(Sa, pres, T_fr, is, npts)
      case default
        call MOM_error(FATAL, "calculate_TFreeze_scalar: form_of_TFreeze is not valid.")
    end select
  endif

  if (EOS%degC_to_C /= 1.0) then
    do i=is,ie ; T_fr(i) = EOS%degC_to_C * T_fr(i) ; enddo
  endif

end subroutine calculate_TFreeze_1d


!> Calls the appropriate subroutine to calculate density derivatives for 1-D array inputs.
subroutine calculate_density_derivs_array(T, S, pressure, drho_dT, drho_dS, start, npts, EOS, scale)
  real, dimension(:), intent(in)    :: T        !< Potential temperature referenced to the surface [degC]
  real, dimension(:), intent(in)    :: S        !< Salinity [ppt]
  real, dimension(:), intent(in)    :: pressure !< Pressure [Pa]
  real, dimension(:), intent(inout) :: drho_dT  !< The partial derivative of density with potential
                                                !! temperature [kg m-3 degC-1] or other units determined
                                                !! by the optional scale argument
  real, dimension(:), intent(inout) :: drho_dS  !< The partial derivative of density with salinity,
                                                !! in [kg m-3 ppt-1] or other units determined
                                                !! by the optional scale argument
  integer,            intent(in)    :: start    !< Starting index within the array
  integer,            intent(in)    :: npts     !< The number of values to calculate
  type(EOS_type),     intent(in)    :: EOS      !< Equation of state structure
  real,     optional, intent(in)    :: scale    !< A multiplicative factor by which to scale density
                                                !! in combination with scaling stored in EOS [various]

  ! Local variables
  integer :: j

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_density_derivs_linear(T, S, pressure, drho_dT, drho_dS, EOS%Rho_T0_S0, &
                                           EOS%dRho_dT, EOS%dRho_dS, start, npts)
    case (EOS_UNESCO)
      call calculate_density_derivs_unesco(T, S, pressure, drho_dT, drho_dS, start, npts)
    case (EOS_WRIGHT)
      call calculate_density_derivs_wright(T, S, pressure, drho_dT, drho_dS, start, npts)
    case (EOS_TEOS10)
      call calculate_density_derivs_teos10(T, S, pressure, drho_dT, drho_dS, start, npts)
    case (EOS_NEMO)
      call calculate_density_derivs_nemo(T, S, pressure, drho_dT, drho_dS, start, npts)
    case default
      call MOM_error(FATAL, "calculate_density_derivs_array: EOS%form_of_EOS is not valid.")
  end select

  if (present(scale)) then ; if (scale /= 1.0) then ; do j=start,start+npts-1
    drho_dT(j) = scale * drho_dT(j)
    drho_dS(j) = scale * drho_dS(j)
  enddo ; endif ; endif

end subroutine calculate_density_derivs_array


!> Calls the appropriate subroutine to calculate density derivatives for 1-D array inputs.
subroutine calculate_density_derivs_1d(T, S, pressure, drho_dT, drho_dS, EOS, dom, scale)
  real, dimension(:),    intent(in)    :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(:),    intent(in)    :: S        !< Salinity [S ~> ppt]
  real, dimension(:),    intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:),    intent(inout) :: drho_dT  !< The partial derivative of density with potential
                                                   !! temperature [R C-1 ~> kg m-3 degC-1]
  real, dimension(:),    intent(inout) :: drho_dS  !< The partial derivative of density with salinity
                                                   !! [R S-1 ~> kg m-3 ppt-1]
  type(EOS_type),        intent(in)    :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.
  real,                  optional, intent(in) :: scale !< A multiplicative factor by which to scale density
                                                       !! in combination with scaling stored in EOS [various]
  ! Local variables
  real, dimension(size(drho_dT)) :: pres  ! Pressure converted to [Pa]
  real, dimension(size(drho_dT)) :: Ta    ! Temperature converted to [degC]
  real, dimension(size(drho_dT)) :: Sa    ! Salinity converted to [ppt]
  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real :: dRdT_scale ! A factor to convert drho_dT to the desired units [R degC m3 C-1 kg-1 ~> 1]
  real :: dRdS_scale ! A factor to convert drho_dS to the desired units [R ppt m3 S-1 kg-1 ~> 1]
  integer :: i, is, ie, npts

  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(drho_dT) ; npts = 1 + ie - is
  endif

  if ((EOS%RL2_T2_to_Pa == 1.0) .and. (EOS%C_to_degC == 1.0) .and. (EOS%S_to_ppt == 1.0)) then
    call calculate_density_derivs_array(T, S, pressure, drho_dT, drho_dS, is, npts, EOS)
  else
    do i=is,ie
      pres(i) = EOS%RL2_T2_to_Pa * pressure(i)
      Ta(i) = EOS%C_to_degC * T(i)
      Sa(i) = EOS%S_to_ppt * S(i)
    enddo
    call calculate_density_derivs_array(Ta, Sa, pres, drho_dT, drho_dS, is, npts, EOS)
  endif

  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  dRdT_scale = rho_scale * EOS%C_to_degC
  dRdS_scale = rho_scale * EOS%S_to_ppt
  if ((dRdT_scale /= 1.0) .or. (dRdS_scale /= 1.0)) then ; do i=is,ie
    drho_dT(i) = dRdT_scale * drho_dT(i)
    drho_dS(i) = dRdS_scale * drho_dS(i)
  enddo ; endif

end subroutine calculate_density_derivs_1d


!> Calls the appropriate subroutines to calculate density derivatives by promoting a scalar
!! to a one-element array
subroutine calculate_density_derivs_scalar(T, S, pressure, drho_dT, drho_dS, EOS, scale)
  real,           intent(in)  :: T !< Potential temperature referenced to the surface [C ~> degC]
  real,           intent(in)  :: S !< Salinity [S ~> ppt]
  real,           intent(in)  :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real,           intent(out) :: drho_dT !< The partial derivative of density with potential
                                         !! temperature [R C-1 ~> kg m-3 degC-1] or other
                                         !! units determined by the optional scale argument
  real,           intent(out) :: drho_dS !< The partial derivative of density with salinity,
                                         !! in [R S-1 ~> kg m-3 ppt-1] or other units
                                         !! determined by the optional scale argument
  type(EOS_type), intent(in)  :: EOS     !< Equation of state structure
  real, optional, intent(in)  :: scale   !< A multiplicative factor by which to scale density
                                         !! in combination with scaling stored in EOS [various]
  ! Local variables
  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real :: dRdT_scale ! A factor to convert drho_dT to the desired units [R degC m3 C-1 kg-1 ~> 1]
  real :: dRdS_scale ! A factor to convert drho_dS to the desired units [R ppt m3 S-1 kg-1 ~> 1]
  real :: pres  ! Pressure converted to [Pa]
  real :: Ta    ! Temperature converted to [degC]
  real :: Sa    ! Salinity converted to [ppt]

  pres = EOS%RL2_T2_to_Pa*pressure
  Ta = EOS%C_to_degC * T
  Sa = EOS%S_to_ppt * S

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_density_derivs_linear(Ta, Sa, pres, drho_dT, drho_dS, &
                                           EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS)
    case (EOS_WRIGHT)
      call calculate_density_derivs_wright(Ta, Sa, pres, drho_dT, drho_dS)
    case (EOS_TEOS10)
      call calculate_density_derivs_teos10(Ta, Sa, pres, drho_dT, drho_dS)
    case default
      call MOM_error(FATAL, "calculate_density_derivs_scalar: EOS%form_of_EOS is not valid.")
  end select

  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  dRdT_scale = rho_scale * EOS%C_to_degC
  dRdS_scale = rho_scale * EOS%S_to_ppt
  if ((dRdT_scale /= 1.0) .or. (dRdS_scale /= 1.0)) then
    drho_dT = dRdT_scale * drho_dT
    drho_dS = dRdS_scale * drho_dS
  endif

end subroutine calculate_density_derivs_scalar

!> Calls the appropriate subroutine to calculate density second derivatives for 1-D array inputs.
subroutine calculate_density_second_derivs_1d(T, S, pressure, drho_dS_dS, drho_dS_dT, drho_dT_dT, &
                                              drho_dS_dP, drho_dT_dP, EOS, dom, scale)
  real, dimension(:), intent(in)  :: T !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(:), intent(in)  :: S !< Salinity [S ~> ppt]
  real, dimension(:), intent(in)  :: pressure   !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:), intent(inout) :: drho_dS_dS !< Partial derivative of beta with respect to S
                                                  !! [R S-2 ~> kg m-3 ppt-2]
  real, dimension(:), intent(inout) :: drho_dS_dT !< Partial derivative of beta with respect to T
                                                  !! [R S-1 C-1 ~> kg m-3 ppt-1 degC-1]
  real, dimension(:), intent(inout) :: drho_dT_dT !< Partial derivative of alpha with respect to T
                                                  !! [R C-2 ~> kg m-3 degC-2]
  real, dimension(:), intent(inout) :: drho_dS_dP !< Partial derivative of beta with respect to pressure
                                                  !! [T2 S-1 L-2 ~> kg m-3 ppt-1 Pa-1]
  real, dimension(:), intent(inout) :: drho_dT_dP !< Partial derivative of alpha with respect to pressure
                                                  !! [T2 C-1 L-2 ~> kg m-3 degC-1 Pa-1]
  type(EOS_type),     intent(in)    :: EOS        !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                  !! into account that arrays start at 1.
  real,     optional, intent(in)    :: scale      !< A multiplicative factor by which to scale density
                                                  !! in combination with scaling stored in EOS [various]
  ! Local variables
  real, dimension(size(T)) :: pres  ! Pressure converted to [Pa]
  real, dimension(size(T)) :: Ta    ! Temperature converted to [degC]
  real, dimension(size(T)) :: Sa    ! Salinity converted to [ppt]
  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  integer :: i, is, ie, npts

  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(T) ; npts = 1 + ie - is
  endif

  if ((EOS%RL2_T2_to_Pa == 1.0) .and. (EOS%C_to_degC == 1.0) .and. (EOS%S_to_ppt == 1.0)) then
    select case (EOS%form_of_EOS)
      case (EOS_LINEAR)
        call calculate_density_second_derivs_linear(T, S, pressure, drho_dS_dS, drho_dS_dT, &
                                                    drho_dT_dT, drho_dS_dP, drho_dT_dP, is, npts)
      case (EOS_WRIGHT)
        call calculate_density_second_derivs_wright(T, S, pressure, drho_dS_dS, drho_dS_dT, &
                                                    drho_dT_dT, drho_dS_dP, drho_dT_dP, is, npts)
      case (EOS_TEOS10)
        call calculate_density_second_derivs_teos10(T, S, pressure, drho_dS_dS, drho_dS_dT, &
                                                    drho_dT_dT, drho_dS_dP, drho_dT_dP, is, npts)
      case default
        call MOM_error(FATAL, "calculate_density_derivs: EOS%form_of_EOS is not valid.")
    end select
  else
    do i=is,ie
      pres(i) = EOS%RL2_T2_to_Pa * pressure(i)
      Ta(i) = EOS%C_to_degC * T(i)
      Sa(i) = EOS%S_to_ppt * S(i)
    enddo
    select case (EOS%form_of_EOS)
      case (EOS_LINEAR)
        call calculate_density_second_derivs_linear(Ta, Sa, pres, drho_dS_dS, drho_dS_dT, &
                                                    drho_dT_dT, drho_dS_dP, drho_dT_dP, is, npts)
      case (EOS_WRIGHT)
        call calculate_density_second_derivs_wright(Ta, Sa, pres, drho_dS_dS, drho_dS_dT, &
                                                    drho_dT_dT, drho_dS_dP, drho_dT_dP, is, npts)
      case (EOS_TEOS10)
        call calculate_density_second_derivs_teos10(Ta, Sa, pres, drho_dS_dS, drho_dS_dT, &
                                                    drho_dT_dT, drho_dS_dP, drho_dT_dP, is, npts)
      case default
        call MOM_error(FATAL, "calculate_density_derivs: EOS%form_of_EOS is not valid.")
    end select
  endif

  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  if (rho_scale /= 1.0) then ; do i=is,ie
    drho_dS_dS(i) = rho_scale * drho_dS_dS(i)
    drho_dS_dT(i) = rho_scale * drho_dS_dT(i)
    drho_dT_dT(i) = rho_scale * drho_dT_dT(i)
    drho_dS_dP(i) = rho_scale * drho_dS_dP(i)
    drho_dT_dP(i) = rho_scale * drho_dT_dP(i)
  enddo ; endif

  if (EOS%RL2_T2_to_Pa /= 1.0) then ; do i=is,ie
    drho_dS_dP(i) = EOS%RL2_T2_to_Pa * drho_dS_dP(i)
    drho_dT_dP(i) = EOS%RL2_T2_to_Pa * drho_dT_dP(i)
  enddo ; endif

  if (EOS%C_to_degC /= 1.0) then ; do i=is,ie
    drho_dS_dT(i) = EOS%C_to_degC * drho_dS_dT(i)
    drho_dT_dT(i) = EOS%C_to_degC**2 * drho_dT_dT(i)
    drho_dT_dP(i) = EOS%C_to_degC * drho_dT_dP(i)
  enddo ; endif

  if (EOS%S_to_ppt /= 1.0) then ; do i=is,ie
    drho_dS_dS(i) = EOS%S_to_ppt**2 * drho_dS_dS(i)
    drho_dS_dT(i) = EOS%S_to_ppt * drho_dS_dT(i)
    drho_dS_dP(i) = EOS%S_to_ppt * drho_dS_dP(i)
  enddo ; endif

end subroutine calculate_density_second_derivs_1d

!> Calls the appropriate subroutine to calculate density second derivatives for scalar inputs.
subroutine calculate_density_second_derivs_scalar(T, S, pressure, drho_dS_dS, drho_dS_dT, drho_dT_dT, &
                                                  drho_dS_dP, drho_dT_dP, EOS, scale)
  real, intent(in)  :: T !< Potential temperature referenced to the surface [C ~> degC]
  real, intent(in)  :: S !< Salinity [S ~> ppt]
  real, intent(in)  :: pressure   !< Pressure [R L2 T-2 ~> Pa]
  real, intent(out) :: drho_dS_dS !< Partial derivative of beta with respect to S
                                  !! [R S-2 ~> kg m-3 ppt-2]
  real, intent(out) :: drho_dS_dT !< Partial derivative of beta with respect to T
                                  !! [R S-1 C-1 ~> kg m-3 ppt-1 degC-1]
  real, intent(out) :: drho_dT_dT !< Partial derivative of alpha with respect to T
                                  !! [R C-2 ~> kg m-3 degC-2]
  real, intent(out) :: drho_dS_dP !< Partial derivative of beta with respect to pressure
                                  !! [T2 S-1 L-2 ~> kg m-3 ppt-1 Pa-1]
  real, intent(out) :: drho_dT_dP !< Partial derivative of alpha with respect to pressure
                                  !! [T2 C-1 L-2 ~> kg m-3 degC-1 Pa-1]
  type(EOS_type), intent(in) :: EOS !< Equation of state structure
  real, optional, intent(in) :: scale !< A multiplicative factor by which to scale density
                                  !! in combination with scaling stored in EOS [various]
  ! Local variables
  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real :: pres  ! Pressure converted to [Pa]
  real :: Ta    ! Temperature converted to [degC]
  real :: Sa    ! Salinity converted to [ppt]

  pres = EOS%RL2_T2_to_Pa*pressure
  Ta = EOS%C_to_degC * T
  Sa = EOS%S_to_ppt * S

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_density_second_derivs_linear(Ta, Sa, pres, drho_dS_dS, drho_dS_dT, &
                                                  drho_dT_dT, drho_dS_dP, drho_dT_dP)
    case (EOS_WRIGHT)
      call calculate_density_second_derivs_wright(Ta, Sa, pres, drho_dS_dS, drho_dS_dT, &
                                                  drho_dT_dT, drho_dS_dP, drho_dT_dP)
    case (EOS_TEOS10)
      call calculate_density_second_derivs_teos10(Ta, Sa, pres, drho_dS_dS, drho_dS_dT, &
                                                  drho_dT_dT, drho_dS_dP, drho_dT_dP)
    case default
      call MOM_error(FATAL, "calculate_density_derivs: EOS%form_of_EOS is not valid.")
  end select

  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  if (rho_scale /= 1.0) then
    drho_dS_dS = rho_scale * drho_dS_dS
    drho_dS_dT = rho_scale * drho_dS_dT
    drho_dT_dT = rho_scale * drho_dT_dT
    drho_dS_dP = rho_scale * drho_dS_dP
    drho_dT_dP = rho_scale * drho_dT_dP
  endif

  if (EOS%RL2_T2_to_Pa /= 1.0) then
    drho_dS_dP = EOS%RL2_T2_to_Pa * drho_dS_dP
    drho_dT_dP = EOS%RL2_T2_to_Pa * drho_dT_dP
  endif

  if (EOS%C_to_degC /= 1.0) then
    drho_dS_dT = EOS%C_to_degC * drho_dS_dT
    drho_dT_dT = EOS%C_to_degC**2 * drho_dT_dT
    drho_dT_dP = EOS%C_to_degC * drho_dT_dP
  endif

  if (EOS%S_to_ppt /= 1.0) then
    drho_dS_dS = EOS%S_to_ppt**2 * drho_dS_dS
    drho_dS_dT = EOS%S_to_ppt * drho_dS_dT
    drho_dS_dP = EOS%S_to_ppt * drho_dS_dP
  endif

end subroutine calculate_density_second_derivs_scalar

!> Calls the appropriate subroutine to calculate specific volume derivatives for an array.
subroutine calculate_spec_vol_derivs_array(T, S, pressure, dSV_dT, dSV_dS, start, npts, EOS)
  real, dimension(:), intent(in)  :: T !< Potential temperature referenced to the surface [degC]
  real, dimension(:), intent(in)  :: S !< Salinity [ppt]
  real, dimension(:), intent(in)  :: pressure !< Pressure [Pa]
  real, dimension(:), intent(inout) :: dSV_dT !< The partial derivative of specific volume with potential
                                              !! temperature [m3 kg-1 degC-1]
  real, dimension(:), intent(inout) :: dSV_dS !< The partial derivative of specific volume with salinity
                                              !! [m3 kg-1 ppt-1]
  integer,            intent(in)  :: start  !< Starting index within the array
  integer,            intent(in)  :: npts   !< The number of values to calculate
  type(EOS_type),     intent(in)  :: EOS    !< Equation of state structure

  ! Local variables
  real, dimension(size(T)) :: rho     ! In situ density [kg m-3]
  real, dimension(size(T)) :: dRho_dT ! Derivative of density with temperature [kg m-3 degC-1]
  real, dimension(size(T)) :: dRho_dS ! Derivative of density with salinity [kg m-3 ppt-1]
  integer :: j

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_specvol_derivs_linear(T, S, pressure, dSV_dT, dSV_dS, start, &
                                           npts, EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS)
    case (EOS_UNESCO)
      call calculate_density_unesco(T, S, pressure, rho, start, npts)
      call calculate_density_derivs_unesco(T, S, pressure, drho_dT, drho_dS, start, npts)
      do j=start,start+npts-1
        dSV_dT(j) = -dRho_DT(j)/(rho(j)**2)
        dSV_dS(j) = -dRho_DS(j)/(rho(j)**2)
      enddo
    case (EOS_WRIGHT)
      call calculate_specvol_derivs_wright(T, S, pressure, dSV_dT, dSV_dS, start, npts)
    case (EOS_TEOS10)
      call calculate_specvol_derivs_teos10(T, S, pressure, dSV_dT, dSV_dS, start, npts)
    case (EOS_NEMO)
      call calculate_density_nemo(T, S, pressure, rho, start, npts)
      call calculate_density_derivs_nemo(T, S, pressure, drho_dT, drho_dS, start, npts)
      do j=start,start+npts-1
        dSV_dT(j) = -dRho_DT(j)/(rho(j)**2)
        dSV_dS(j) = -dRho_DS(j)/(rho(j)**2)
      enddo
    case default
      call MOM_error(FATAL, "calculate_spec_vol_derivs_array: EOS%form_of_EOS is not valid.")
  end select

end subroutine calculate_spec_vol_derivs_array

!> Calls the appropriate subroutine to calculate specific volume derivatives for 1-d array inputs,
!! potentially limiting the domain of indices that are worked on.
subroutine calc_spec_vol_derivs_1d(T, S, pressure, dSV_dT, dSV_dS, EOS, dom, scale)
  real, dimension(:), intent(in)    :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(:), intent(in)    :: S        !< Salinity [S ~> ppt]
  real, dimension(:), intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:), intent(inout) :: dSV_dT   !< The partial derivative of specific volume with potential
                                                !! temperature [R-1 C-1 ~> m3 kg-1 degC-1]
  real, dimension(:), intent(inout) :: dSV_dS   !< The partial derivative of specific volume with salinity
                                                !! [R-1 S-1 ~> m3 kg-1 ppt-1]
  type(EOS_type),     intent(in)    :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.
  real,                  optional, intent(in) :: scale !< A multiplicative factor by which to scale specific
                                                !! volume in combination with scaling stored in EOS [various]

  ! Local variables
  real, dimension(size(T)) :: pres  ! Pressure converted to [Pa]
  real, dimension(size(T)) :: Ta    ! Temperature converted to [degC]
  real, dimension(size(T)) :: Sa    ! Salinity converted to [ppt]
  real :: spv_scale ! A factor to convert specific volume from m3 kg-1 to the desired units [kg R-1 m-3 ~> 1]
  real :: dSVdT_scale ! A factor to convert dSV_dT to the desired units [kg degC R-1 C-1 m-3 ~> 1]
  real :: dSVdS_scale ! A factor to convert dSV_dS to the desired units [kg ppt R-1 S-1 m-3 ~> 1]
  integer :: i, is, ie, npts

  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(dSV_dT) ; npts = 1 + ie - is
  endif

  if ((EOS%RL2_T2_to_Pa == 1.0) .and. (EOS%C_to_degC == 1.0) .and. (EOS%S_to_ppt == 1.0)) then
    call calculate_spec_vol_derivs_array(T, S, pressure, dSV_dT, dSV_dS, is, npts, EOS)
  else
    do i=is,ie
      pres(i) = EOS%RL2_T2_to_Pa * pressure(i)
      Ta(i) = EOS%C_to_degC * T(i)
      Sa(i) = EOS%S_to_ppt * S(i)
    enddo
    call calculate_spec_vol_derivs_array(Ta, Sa, pres, dSV_dT, dSV_dS, is, npts, EOS)
  endif

  spv_scale = EOS%R_to_kg_m3
  if (present(scale)) spv_scale = spv_scale * scale
  dSVdT_scale = spv_scale * EOS%C_to_degC
  dSVdS_scale = spv_scale * EOS%S_to_ppt
  if ((dSVdT_scale /= 1.0) .or. (dSVdS_scale /= 1.0)) then ; do i=is,ie
    dSV_dT(i) = dSVdT_scale * dSV_dT(i)
    dSV_dS(i) = dSVdS_scale * dSV_dS(i)
  enddo ; endif

end subroutine calc_spec_vol_derivs_1d


!> Calls the appropriate subroutine to calculate the density and compressibility for 1-D array
!! inputs.  The inputs and outputs use dimensionally rescaled units.
subroutine calculate_compress_1d(T, S, pressure, rho, drho_dp, EOS, dom)
  real, dimension(:), intent(in)    :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(:), intent(in)    :: S        !< Salinity [S ~> ppt]
  real, dimension(:), intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:), intent(inout) :: rho      !< In situ density [R ~> kg m-3]
  real, dimension(:), intent(inout) :: drho_dp  !< The partial derivative of density with pressure
                                                !! (also the inverse of the square of sound speed)
                                                !! [T2 L-2 ~> s2 m-2]
  type(EOS_type),     intent(in)  :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.

  ! Local variables
  real, dimension(size(T)) :: pres  ! Pressure converted to [Pa]
  real, dimension(size(T)) :: Ta    ! Temperature converted to [degC]
  real, dimension(size(T)) :: Sa    ! Salinity converted to [ppt]
  integer :: i, is, ie, npts

  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(rho) ; npts = 1 + ie - is
  endif

  do i=is,ie
    pres(i) = EOS%RL2_T2_to_Pa * pressure(i)
    Ta(i) = EOS%C_to_degC * T(i)
    Sa(i) = EOS%S_to_ppt * S(i)
  enddo

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_compress_linear(Ta, Sa, pres, rho, drho_dp, is, npts, &
                                     EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS)
    case (EOS_UNESCO)
      call calculate_compress_unesco(Ta, Sa, pres, rho, drho_dp, is, npts)
    case (EOS_WRIGHT)
      call calculate_compress_wright(Ta, Sa, pres, rho, drho_dp, is, npts)
    case (EOS_TEOS10)
      call calculate_compress_teos10(Ta, Sa, pres, rho, drho_dp, is, npts)
    case (EOS_NEMO)
      call calculate_compress_nemo(Ta, Sa, pres, rho, drho_dp, is, npts)
    case default
      call MOM_error(FATAL, "calculate_compress: EOS%form_of_EOS is not valid.")
  end select

  if (EOS%kg_m3_to_R /= 1.0) then ; do i=is,ie
    rho(i) = EOS%kg_m3_to_R * rho(i)
  enddo ; endif
  if (EOS%L_T_to_m_s /= 1.0) then ; do i=is,ie
    drho_dp(i) = EOS%L_T_to_m_s**2 * drho_dp(i)
  enddo ; endif

end subroutine calculate_compress_1d

!> Calculate density and compressibility for a scalar. This just promotes the scalar to an array
!! with a singleton dimension and calls calculate_compress_1d.  The inputs and outputs use
!! dimensionally rescaled units.
subroutine calculate_compress_scalar(T, S, pressure, rho, drho_dp, EOS)
  real, intent(in)        :: T        !< Potential temperature referenced to the surface [C ~> degC]
  real, intent(in)        :: S        !< Salinity [S ~> ppt]
  real, intent(in)        :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, intent(out)       :: rho      !< In situ density [R ~> kg m-3]
  real, intent(out)       :: drho_dp  !< The partial derivative of density with pressure (also the
                                      !! inverse of the square of sound speed) [T2 L-2 ~> s2 m-2]
  type(EOS_type), intent(in) :: EOS   !< Equation of state structure

  ! Local variables
  ! These arrays use the same units as their counterparts in calculate_compress_1d.
  real, dimension(1) :: pa    ! Pressure in a size-1 1d array [R L2 T-2 ~> Pa]
  real, dimension(1) :: Ta    ! Temperature in a size-1 1d array [C ~> degC]
  real, dimension(1) :: Sa    ! Salinity in a size-1 1d array [S ~> ppt]
  real, dimension(1) :: rhoa  ! In situ density in a size-1 1d array [R ~> kg m-3]
  real, dimension(1) :: drho_dpa ! The partial derivative of density with pressure (also the
                              ! inverse of the square of sound speed) in a 1d array [T2 L-2 ~> s2 m-2]

  Ta(1) = T ; Sa(1) = S ; pa(1) = pressure

  call calculate_compress_1d(Ta, Sa, pa, rhoa, drho_dpa, EOS)
  rho = rhoa(1) ; drho_dp = drho_dpa(1)

end subroutine calculate_compress_scalar


!> This subroutine returns a two point integer array indicating the domain of i-indices
!! to work on in EOS calls based on information from a hor_index type
function EOS_domain(HI, halo) result(EOSdom)
  type(hor_index_type), intent(in)  :: HI    !< The horizontal index structure
  integer,    optional, intent(in)  :: halo  !< The halo size to work on; missing is equivalent to 0.
  integer, dimension(2) :: EOSdom   !< The index domain that the EOS will work on, taking into account
                                    !! that the arrays inside the EOS routines will start at 1.

  ! Local variables
  integer :: halo_sz

  halo_sz = 0 ; if (present(halo)) halo_sz = halo

  EOSdom(1) = HI%isc - (HI%isd-1) - halo_sz
  EOSdom(2) = HI%iec - (HI%isd-1) + halo_sz

end function EOS_domain

!> Calls the appropriate subroutine to calculate analytical and nearly-analytical
!! integrals in pressure across layers of geopotential anomalies, which are
!! required for calculating the finite-volume form pressure accelerations in a
!! non-Boussinesq model.  There are essentially no free assumptions, apart from the
!! use of Boole's rule to do the horizontal integrals, and from a truncation in the
!! series for log(1-eps/1+eps) that assumes that |eps| < 0.34.
subroutine analytic_int_specific_vol_dp(T, S, p_t, p_b, alpha_ref, HI, EOS, &
                               dza, intp_dza, intx_dza, inty_dza, halo_size, &
                               bathyP, dP_tiny, useMassWghtInterp)
  type(hor_index_type), intent(in)  :: HI  !< The horizontal index structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T   !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S   !< Salinity [S ~> ppt]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: p_t !< Pressure at the top of the layer [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: p_b !< Pressure at the bottom of the layer [R L2 T-2 ~> Pa]
  real,                 intent(in)  :: alpha_ref !< A mean specific volume that is subtracted out
                            !! to reduce the magnitude of each of the integrals [R-1 ~> m3 kg-1]
                            !! The calculation is mathematically identical with different values of
                            !! alpha_ref, but this reduces the effects of roundoff.
  type(EOS_type),       intent(in)  :: EOS !< Equation of state structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(inout) :: dza !< The change in the geopotential anomaly across
                            !! the layer [L2 T-2 ~> m2 s-2]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(inout) :: intp_dza !< The integral in pressure through the layer of the
                            !! geopotential anomaly relative to the anomaly at the bottom of the
                            !! layer [R L4 T-4 ~> Pa m2 s-2]
  real, dimension(HI%IsdB:HI%IedB,HI%jsd:HI%jed), &
              optional, intent(inout) :: intx_dza !< The integral in x of the difference between the
                            !! geopotential anomaly at the top and bottom of the layer divided by
                            !! the x grid spacing [L2 T-2 ~> m2 s-2]
  real, dimension(HI%isd:HI%ied,HI%JsdB:HI%JedB), &
              optional, intent(inout) :: inty_dza !< The integral in y of the difference between the
                            !! geopotential anomaly at the top and bottom of the layer divided by
                            !! the y grid spacing [L2 T-2 ~> m2 s-2]
  integer,    optional, intent(in)  :: halo_size !< The width of halo points on which to calculate dza.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: bathyP  !< The pressure at the bathymetry [R L2 T-2 ~> Pa]
  real,       optional, intent(in)  :: dP_tiny !< A miniscule pressure change with
                            !! the same units as p_t [R L2 T-2 ~> Pa]
  logical,    optional, intent(in)  :: useMassWghtInterp !< If true, uses mass weighting
                            !! to interpolate T/S for top and bottom integrals.

  ! Local variables
  real :: dRdT_scale ! A factor to convert drho_dT to the desired units [R degC m3 C-1 kg-1 ~> 1]
  real :: dRdS_scale ! A factor to convert drho_dS to the desired units [R ppt m3 S-1 kg-1 ~> 1]



  ! We should never reach this point with quadrature. EOS_quadrature indicates that numerical
  ! integration be used instead of analytic. This is a safety check.
  if (EOS%EOS_quadrature) call MOM_error(FATAL, "EOS_quadrature is set!")

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      dRdT_scale = EOS%kg_m3_to_R * EOS%C_to_degC
      dRdS_scale = EOS%kg_m3_to_R * EOS%S_to_ppt
      call int_spec_vol_dp_linear(T, S, p_t, p_b, alpha_ref, HI, EOS%kg_m3_to_R*EOS%Rho_T0_S0, &
                                dRdT_scale*EOS%dRho_dT, dRdS_scale*EOS%dRho_dS, dza, &
                                intp_dza, intx_dza, inty_dza, halo_size, &
                                bathyP, dP_tiny, useMassWghtInterp)
    case (EOS_WRIGHT)
      call int_spec_vol_dp_wright(T, S, p_t, p_b, alpha_ref, HI, dza, intp_dza, intx_dza, &
                                  inty_dza, halo_size, bathyP, dP_tiny, useMassWghtInterp, &
                                  SV_scale=EOS%R_to_kg_m3, pres_scale=EOS%RL2_T2_to_Pa, &
                                  temp_scale=EOS%C_to_degC, saln_scale=EOS%S_to_ppt)
    case default
      call MOM_error(FATAL, "No analytic integration option is available with this EOS!")
  end select

end subroutine analytic_int_specific_vol_dp

!> This subroutine calculates analytical and nearly-analytical integrals of
!! pressure anomalies across layers, which are required for calculating the
!! finite-volume form pressure accelerations in a Boussinesq model.
subroutine analytic_int_density_dz(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, EOS, dpa, &
                          intz_dpa, intx_dpa, inty_dpa, bathyT, dz_neglect, useMassWghtInterp, Z_0p)
  type(hor_index_type), intent(in)  :: HI !< Ocean horizontal index structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T   !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S   !< Salinity [S ~> ppt]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_t !< Height at the top of the layer in depth units [Z ~> m]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_b !< Height at the bottom of the layer [Z ~> m]
  real,                 intent(in)  :: rho_ref !< A mean density [R ~> kg m-3], that is
                                           !! subtracted out to reduce the magnitude of each of the
                                           !! integrals.
  real,                 intent(in)  :: rho_0 !< A density [R ~> kg m-3], that is used
                                           !! to calculate the pressure (as p~=-z*rho_0*G_e)
                                           !! used in the equation of state.
  real,                 intent(in)  :: G_e !< The Earth's gravitational acceleration
                                           !! [L2 Z-1 T-2 ~> m s-2]
  type(EOS_type),       intent(in)  :: EOS !< Equation of state structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                      intent(inout) :: dpa !< The change in the pressure anomaly
                                           !! across the layer [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
            optional, intent(inout) :: intz_dpa !< The integral through the thickness of the
                                           !! layer of the pressure anomaly relative to the
                                           !! anomaly at the top of the layer [R L2 Z T-2 ~> Pa m]
  real, dimension(HI%IsdB:HI%IedB,HI%jsd:HI%jed), &
            optional, intent(inout) :: intx_dpa !< The integral in x of the difference between
                                          !! the pressure anomaly at the top and bottom of the
                                          !! layer divided by the x grid spacing [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%JsdB:HI%JedB), &
            optional, intent(inout) :: inty_dpa !< The integral in y of the difference between
                                          !! the pressure anomaly at the top and bottom of the
                                          !! layer divided by the y grid spacing [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: bathyT !< The depth of the bathymetry [Z ~> m]
  real,       optional, intent(in)  :: dz_neglect !< A miniscule thickness change [Z ~> m]
  logical,    optional, intent(in)  :: useMassWghtInterp !< If true, uses mass weighting to
                                           !! interpolate T/S for top and bottom integrals.
  real,       optional, intent(in)  :: Z_0p !< The height at which the pressure is 0 [Z ~> m]

  ! Local variables
  real :: rho_scale  ! A multiplicative factor by which to scale density from kg m-3 to the
                     ! desired units [R m3 kg-1 ~> 1]
  real :: dRdT_scale ! A factor to convert drho_dT to the desired units [R degC m3 C-1 kg-1 ~> 1]
  real :: dRdS_scale ! A factor to convert drho_dS to the desired units [R ppt m3 S-1 kg-1 ~> 1]
  real :: pres_scale ! A multiplicative factor to convert pressure into Pa [Pa T2 R-1 L-2 ~> 1]

  ! We should never reach this point with quadrature. EOS_quadrature indicates that numerical
  ! integration be used instead of analytic. This is a safety check.
  if (EOS%EOS_quadrature) call MOM_error(FATAL, "EOS_quadrature is set!")

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      rho_scale = EOS%kg_m3_to_R
      dRdT_scale = EOS%kg_m3_to_R * EOS%C_to_degC
      dRdS_scale = EOS%kg_m3_to_R * EOS%S_to_ppt
      if ((rho_scale /= 1.0) .or. (dRdT_scale /= 1.0) .or. (dRdS_scale /= 1.0)) then
        call int_density_dz_linear(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                         rho_scale*EOS%Rho_T0_S0, dRdT_scale*EOS%dRho_dT, dRdS_scale*EOS%dRho_dS, &
                         dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, dz_neglect, useMassWghtInterp)
      else
        call int_density_dz_linear(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                         EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS, &
                         dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, dz_neglect, useMassWghtInterp)
      endif
    case (EOS_WRIGHT)
      rho_scale = EOS%kg_m3_to_R
      pres_scale = EOS%RL2_T2_to_Pa
      if ((rho_scale /= 1.0) .or. (pres_scale /= 1.0) .or. (EOS%C_to_degC /= 1.0) .or. (EOS%S_to_ppt /= 1.0)) then
        call int_density_dz_wright(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                                   dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, &
                                   dz_neglect, useMassWghtInterp, rho_scale, pres_scale, &
                                   temp_scale=EOS%C_to_degC, saln_scale=EOS%S_to_ppt, Z_0p=Z_0p)
      else
        call int_density_dz_wright(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                                   dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, &
                                   dz_neglect, useMassWghtInterp, Z_0p=Z_0p)
      endif
    case default
      call MOM_error(FATAL, "No analytic integration option is available with this EOS!")
  end select

end subroutine analytic_int_density_dz

!> Returns true if the equation of state is compressible (i.e. has pressure dependence)
logical function query_compressible(EOS)
  type(EOS_type), intent(in) :: EOS !< Equation of state structure

  query_compressible = EOS%compressible
end function query_compressible

!> Initializes EOS_type by allocating and reading parameters.  The scaling factors in
!! US are stored in EOS for later use.
subroutine EOS_init(param_file, EOS, US)
  type(param_file_type), intent(in) :: param_file !< Parameter file structure
  type(EOS_type), intent(inout)     :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US  !< A dimensional unit scaling type
  optional :: US
  ! Local variables
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_EOS" ! This module's name.
  character(len=40)  :: tmpstr

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "EQN_OF_STATE", tmpstr, &
                 "EQN_OF_STATE determines which ocean equation of state "//&
                 "should be used.  Currently, the valid choices are "//&
                 '"LINEAR", "UNESCO", "WRIGHT", "NEMO" and "TEOS10". '//&
                 "This is only used if USE_EOS is true.", default=EOS_DEFAULT)
  select case (uppercase(tmpstr))
    case (EOS_LINEAR_STRING)
      EOS%form_of_EOS = EOS_LINEAR
    case (EOS_UNESCO_STRING)
      EOS%form_of_EOS = EOS_UNESCO
    case (EOS_WRIGHT_STRING)
      EOS%form_of_EOS = EOS_WRIGHT
    case (EOS_TEOS10_STRING)
      EOS%form_of_EOS = EOS_TEOS10
    case (EOS_NEMO_STRING)
      EOS%form_of_EOS = EOS_NEMO
    case default
      call MOM_error(FATAL, "interpret_eos_selection: EQN_OF_STATE "//&
                              trim(tmpstr) // "in input file is invalid.")
  end select
  call MOM_mesg('interpret_eos_selection: equation of state set to "' // &
                trim(tmpstr)//'"', 5)

  if (EOS%form_of_EOS == EOS_LINEAR) then
    EOS%Compressible = .false.
    call get_param(param_file, mdl, "RHO_T0_S0", EOS%Rho_T0_S0, &
                 "When EQN_OF_STATE="//trim(EOS_LINEAR_STRING)//", "//&
                 "this is the density at T=0, S=0.", units="kg m-3", &
                 default=1000.0)
    call get_param(param_file, mdl, "DRHO_DT", EOS%dRho_dT, &
                 "When EQN_OF_STATE="//trim(EOS_LINEAR_STRING)//", "//&
                 "this is the partial derivative of density with "//&
                 "temperature.", units="kg m-3 K-1", default=-0.2)
    call get_param(param_file, mdl, "DRHO_DS", EOS%dRho_dS, &
                 "When EQN_OF_STATE="//trim(EOS_LINEAR_STRING)//", "//&
                 "this is the partial derivative of density with "//&
                 "salinity.", units="kg m-3 PSU-1", default=0.8)
  endif

  call get_param(param_file, mdl, "EOS_QUADRATURE", EOS%EOS_quadrature, &
                 "If true, always use the generic (quadrature) code "//&
                 "code for the integrals of density.", default=.false.)

  call get_param(param_file, mdl, "TFREEZE_FORM", tmpstr, &
                 "TFREEZE_FORM determines which expression should be "//&
                 "used for the freezing point.  Currently, the valid "//&
                 'choices are "LINEAR", "MILLERO_78", "TEOS10"', &
                 default=TFREEZE_DEFAULT)
  select case (uppercase(tmpstr))
    case (TFREEZE_LINEAR_STRING)
      EOS%form_of_TFreeze = TFREEZE_LINEAR
    case (TFREEZE_MILLERO_STRING)
      EOS%form_of_TFreeze = TFREEZE_MILLERO
    case (TFREEZE_TEOS10_STRING)
      EOS%form_of_TFreeze = TFREEZE_TEOS10
    case default
      call MOM_error(FATAL, "interpret_eos_selection:  TFREEZE_FORM "//&
                              trim(tmpstr) // "in input file is invalid.")
  end select

  if (EOS%form_of_TFreeze == TFREEZE_LINEAR) then
    call get_param(param_file, mdl, "TFREEZE_S0_P0",EOS%TFr_S0_P0, &
                 "When TFREEZE_FORM="//trim(TFREEZE_LINEAR_STRING)//", "//&
                 "this is the freezing potential temperature at "//&
                 "S=0, P=0.", units="deg C", default=0.0)
    call get_param(param_file, mdl, "DTFREEZE_DS",EOS%dTFr_dS, &
                 "When TFREEZE_FORM="//trim(TFREEZE_LINEAR_STRING)//", "//&
                 "this is the derivative of the freezing potential "//&
                 "temperature with salinity.", &
                 units="deg C PSU-1", default=-0.054)
    call get_param(param_file, mdl, "DTFREEZE_DP",EOS%dTFr_dP, &
                 "When TFREEZE_FORM="//trim(TFREEZE_LINEAR_STRING)//", "//&
                 "this is the derivative of the freezing potential "//&
                 "temperature with pressure.", &
                 units="deg C Pa-1", default=0.0)
  endif

  if ((EOS%form_of_EOS == EOS_TEOS10 .OR. EOS%form_of_EOS == EOS_NEMO) .AND. &
      EOS%form_of_TFreeze /= TFREEZE_TEOS10) then
      call MOM_error(FATAL, "interpret_eos_selection:  EOS_TEOS10 or EOS_NEMO \n" //&
      "should only be used along with TFREEZE_FORM = TFREEZE_TEOS10 .")
  endif

  ! Unit conversions
  EOS%m_to_Z = 1. ; if (present(US)) EOS%m_to_Z = US%m_to_Z
  EOS%kg_m3_to_R = 1. ; if (present(US)) EOS%kg_m3_to_R = US%kg_m3_to_R
  EOS%R_to_kg_m3 = 1. ; if (present(US)) EOS%R_to_kg_m3 = US%R_to_kg_m3
  EOS%RL2_T2_to_Pa = 1. ; if (present(US)) EOS%RL2_T2_to_Pa = US%RL2_T2_to_Pa
  EOS%L_T_to_m_s = 1. ; if (present(US)) EOS%L_T_to_m_s = US%L_T_to_m_s
  EOS%degC_to_C = 1. ; if (present(US)) EOS%degC_to_C = US%degC_to_C
  EOS%C_to_degC = 1. ; if (present(US)) EOS%C_to_degC = US%C_to_degC
  EOS%ppt_to_S = 1. ; if (present(US)) EOS%ppt_to_S = US%ppt_to_S
  EOS%S_to_ppt = 1. ; if (present(US)) EOS%S_to_ppt = US%S_to_ppt

end subroutine EOS_init

!> Manually initialized an EOS type (intended for unit testing of routines which need a specific EOS)
subroutine EOS_manual_init(EOS, form_of_EOS, form_of_TFreeze, EOS_quadrature, Compressible, &
                           Rho_T0_S0, drho_dT, dRho_dS, TFr_S0_P0, dTFr_dS, dTFr_dp)
  type(EOS_type),    intent(inout) :: EOS !< Equation of state structure
  integer, optional, intent(in) :: form_of_EOS !< A coded integer indicating the equation of state to use.
  integer, optional, intent(in) :: form_of_TFreeze !< A coded integer indicating the expression for
                                       !! the potential temperature of the freezing point.
  logical, optional, intent(in) :: EOS_quadrature !< If true, always use the generic (quadrature)
                                       !! code for the integrals of density.
  logical, optional, intent(in) :: Compressible  !< If true, in situ density is a function of pressure.
  real   , optional, intent(in) :: Rho_T0_S0 !< Density at T=0 degC and S=0 ppt [kg m-3]
  real   , optional, intent(in) :: drho_dT   !< Partial derivative of density with temperature
                                             !! in [kg m-3 degC-1]
  real   , optional, intent(in) :: dRho_dS   !< Partial derivative of density with salinity
                                             !! in [kg m-3 ppt-1]
  real   , optional, intent(in) :: TFr_S0_P0 !< The freezing potential temperature at S=0, P=0 [degC]
  real   , optional, intent(in) :: dTFr_dS   !< The derivative of freezing point with salinity
                                             !! in [degC ppt-1]
  real   , optional, intent(in) :: dTFr_dp   !< The derivative of freezing point with pressure
                                             !! in [degC Pa-1]

  if (present(form_of_EOS    ))  EOS%form_of_EOS     = form_of_EOS
  if (present(form_of_TFreeze))  EOS%form_of_TFreeze = form_of_TFreeze
  if (present(EOS_quadrature ))  EOS%EOS_quadrature  = EOS_quadrature
  if (present(Compressible   ))  EOS%Compressible    = Compressible
  if (present(Rho_T0_S0      ))  EOS%Rho_T0_S0       = Rho_T0_S0
  if (present(drho_dT        ))  EOS%drho_dT         = drho_dT
  if (present(dRho_dS        ))  EOS%dRho_dS         = dRho_dS
  if (present(TFr_S0_P0      ))  EOS%TFr_S0_P0       = TFr_S0_P0
  if (present(dTFr_dS        ))  EOS%dTFr_dS         = dTFr_dS
  if (present(dTFr_dp        ))  EOS%dTFr_dp         = dTFr_dp

end subroutine EOS_manual_init

!> Set equation of state structure (EOS) to linear with given coefficients
!!
!! \note This routine is primarily for testing and allows a local copy of the
!! EOS_type (EOS argument) to be set to use the linear equation of state
!! independent from the rest of the model.
subroutine EOS_use_linear(Rho_T0_S0, dRho_dT, dRho_dS, EOS, use_quadrature)
  real,              intent(in) :: Rho_T0_S0 !< Density at T=0 degC and S=0 ppt [kg m-3]
  real,              intent(in) :: dRho_dT   !< Partial derivative of density with temperature [kg m-3 degC-1]
  real,              intent(in) :: dRho_dS   !< Partial derivative of density with salinity [kg m-3 ppt-1]
  logical, optional, intent(in) :: use_quadrature !< If true, always use the generic (quadrature)
                                             !! code for the integrals of density.
  type(EOS_type),    intent(inout) :: EOS    !< Equation of state structure

  EOS%form_of_EOS = EOS_LINEAR
  EOS%Compressible = .false.
  EOS%Rho_T0_S0 = Rho_T0_S0
  EOS%dRho_dT = dRho_dT
  EOS%dRho_dS = dRho_dS
  EOS%EOS_quadrature = .false.
  if (present(use_quadrature)) EOS%EOS_quadrature = use_quadrature

end subroutine EOS_use_linear


!> Convert T&S to Absolute Salinity and Conservative Temperature if using TEOS10
subroutine convert_temp_salt_for_TEOS10(T, S, HI, kd, mask_z, EOS)
  integer,               intent(in)    :: kd  !< The number of layers to work on
  type(hor_index_type),  intent(in)    :: HI       !< The horizontal index structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed,kd), &
                         intent(inout) :: T   !< Potential temperature referenced to the surface [C ~> degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed,kd), &
                         intent(inout) :: S   !< Salinity [S ~> ppt]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed,kd), &
                         intent(in)    :: mask_z !< 3d mask regulating which points to convert [nondim]
  type(EOS_type),        intent(in)    :: EOS !< Equation of state structure

  real :: gsw_sr_from_sp ! Reference salinity after conversion from practical salinity [ppt]
  real :: gsw_ct_from_pt ! Conservative temperature after conversion from potential temperature [degC]
  integer :: i, j, k

  if ((EOS%form_of_EOS /= EOS_TEOS10) .and. (EOS%form_of_EOS /= EOS_NEMO)) return

  do k=1,kd ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
    if (mask_z(i,j,k) >= 1.0) then
      S(i,j,k) = EOS%ppt_to_S*gsw_sr_from_sp(EOS%S_to_ppt*S(i,j,k))
!     Get absolute salinity from practical salinity, converting pressures from Pascal to dbar.
!     If this option is activated, pressure will need to be added as an argument, and it should be
!     moved out into module that is not shared between components, where the ocean_grid can be used.
!     S(i,j,k) = gsw_sa_from_sp(S(i,j,k),pres(i,j,k)*1.0e-4,G%geoLonT(i,j),G%geoLatT(i,j))
      T(i,j,k) = EOS%degC_to_C*gsw_ct_from_pt(EOS%S_to_ppt*S(i,j,k), EOS%S_to_ppt*T(i,j,k))
    endif
  enddo ; enddo ; enddo
end subroutine convert_temp_salt_for_TEOS10


!> Converts an array of conservative temperatures to potential temperatures.  The input arguments
!! use the dimesionally rescaling as specified within the EOS type.  The output potential
!! temperature uses this same scaling, but this can be replaced by the factor given by scale.
subroutine cons_temp_to_pot_temp(T, S, poTemp, EOS, dom, scale)
  real, dimension(:), intent(in)    :: T        !< Conservative temperature [C ~> degC]
  real, dimension(:), intent(in)    :: S        !< Absolute salinity [S ~> ppt]
  real, dimension(:), intent(inout) :: poTemp   !< The potential temperature with a reference pressure
                                                !! of 0 Pa, [C ~> degC]
  type(EOS_type),     intent(in)    :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom  !< The domain of indices to work on, taking
                                                !! into account that arrays start at 1.
  real,     optional, intent(in)    :: scale    !< A multiplicative factor by which to scale the output
                                                !! potential temperature in place of with scaling stored
                                                !! in EOS.  A value of 1.0 returns temperatures in [degC],
                                                !! while the default is equivalent to EOS%degC_to_C.

  ! Local variables
  real, dimension(size(T)) :: Ta    ! Temperature converted to [degC]
  real, dimension(size(S)) :: Sa    ! Salinity converted to [ppt]
  real :: T_scale ! A factor to convert potential temperature from degC to the desired units [C degC-1 ~> 1]
  integer :: i, is, ie

  if (present(dom)) then
    is = dom(1) ; ie = dom(2)
  else
    is = 1 ; ie = size(T)
  endif

  if ((EOS%C_to_degC == 1.0) .and. (EOS%S_to_ppt == 1.0)) then
    poTemp(is:ie) = gsw_pt_from_ct(S(is:ie), T(is:ie))
  else
    do i=is,ie
      Ta(i) = EOS%C_to_degC * T(i)
      Sa(i) = EOS%S_to_ppt * S(i)
    enddo
    poTemp(is:ie) = gsw_pt_from_ct(Sa(is:ie), Ta(is:ie))
  endif

  T_scale = EOS%degC_to_C
  if (present(scale)) T_scale = scale
  if (T_scale /= 1.0) then ; do i=is,ie
    poTemp(i) = T_scale * poTemp(i)
  enddo ; endif

end subroutine cons_temp_to_pot_temp


!> Converts an array of absolute salinity to practical salinity.  The input arguments
!! use the dimesionally rescaling as specified within the EOS type.  The output potential
!! temperature uses this same scaling, but this can be replaced by the factor given by scale.
subroutine abs_saln_to_prac_saln(S, prSaln, EOS, dom, scale)
  real, dimension(:), intent(in)    :: S        !< Absolute salinity [S ~> ppt]
  real, dimension(:), intent(inout) :: prSaln   !< Practical salinity [S ~> ppt]
  type(EOS_type),     intent(in)    :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom  !< The domain of indices to work on, taking
                                                !! into account that arrays start at 1.
  real,     optional, intent(in)    :: scale    !< A multiplicative factor by which to scale the output
                                                !! practical in place of with scaling stored
                                                !! in EOS.  A value of 1.0 returns salinities in [PSU],
                                                !! while the default is equivalent to EOS%ppt_to_S.

  ! Local variables
  real, dimension(size(S)) :: Sa    ! Salinity converted to [ppt]
  real :: S_scale ! A factor to convert practical salinity from ppt to the desired units [S ppt-1 ~> 1]
  integer :: i, is, ie

  if (present(dom)) then
    is = dom(1) ; ie = dom(2)
  else
    is = 1 ; ie = size(S)
  endif

  if ((EOS%C_to_degC == 1.0) .and. (EOS%S_to_ppt == 1.0)) then
    prSaln(is:ie) = gsw_sp_from_sr(Sa(is:ie))
  else
    do i=is,ie ; Sa(i) = EOS%S_to_ppt * S(i) ; enddo
    prSaln(is:ie) = gsw_sp_from_sr(Sa(is:ie))
  endif

  S_scale = EOS%ppt_to_S
  if (present(scale)) S_scale = scale
  if (S_scale /= 1.0) then ; do i=is,ie
    prSaln(i) = S_scale * prSaln(i)
  enddo ; endif

end subroutine abs_saln_to_prac_saln


!> Return value of EOS_quadrature
logical function EOS_quadrature(EOS)
  type(EOS_type), intent(in) :: EOS   !< Equation of state structure

  EOS_quadrature  = EOS%EOS_quadrature

end function EOS_quadrature

!> Extractor routine for the EOS type if the members need to be accessed outside this module
subroutine extract_member_EOS(EOS, form_of_EOS, form_of_TFreeze, EOS_quadrature, Compressible, &
                              Rho_T0_S0, drho_dT, dRho_dS, TFr_S0_P0, dTFr_dS, dTFr_dp)
  type(EOS_type),    intent(in)  :: EOS !< Equation of state structure
  integer, optional, intent(out) :: form_of_EOS !< A coded integer indicating the equation of state to use.
  integer, optional, intent(out) :: form_of_TFreeze !< A coded integer indicating the expression for
                                       !! the potential temperature of the freezing point.
  logical, optional, intent(out) :: EOS_quadrature !< If true, always use the generic (quadrature)
                                       !! code for the integrals of density.
  logical, optional, intent(out) :: Compressible !< If true, in situ density is a function of pressure.
  real   , optional, intent(out) :: Rho_T0_S0 !< Density at T=0 degC and S=0 ppt [kg m-3]
  real   , optional, intent(out) :: drho_dT   !< Partial derivative of density with temperature
                                              !! in [kg m-3 degC-1]
  real   , optional, intent(out) :: dRho_dS   !< Partial derivative of density with salinity
                                              !! in [kg m-3 ppt-1]
  real   , optional, intent(out) :: TFr_S0_P0 !< The freezing potential temperature at S=0, P=0 [degC]
  real   , optional, intent(out) :: dTFr_dS   !< The derivative of freezing point with salinity
                                              !! [degC PSU-1]
  real   , optional, intent(out) :: dTFr_dp   !< The derivative of freezing point with pressure
                                              !! [degC Pa-1]

  if (present(form_of_EOS    ))  form_of_EOS     = EOS%form_of_EOS
  if (present(form_of_TFreeze))  form_of_TFreeze = EOS%form_of_TFreeze
  if (present(EOS_quadrature ))  EOS_quadrature  = EOS%EOS_quadrature
  if (present(Compressible   ))  Compressible    = EOS%Compressible
  if (present(Rho_T0_S0      ))  Rho_T0_S0       = EOS%Rho_T0_S0
  if (present(drho_dT        ))  drho_dT         = EOS%drho_dT
  if (present(dRho_dS        ))  dRho_dS         = EOS%dRho_dS
  if (present(TFr_S0_P0      ))  TFr_S0_P0       = EOS%TFr_S0_P0
  if (present(dTFr_dS        ))  dTFr_dS         = EOS%dTFr_dS
  if (present(dTFr_dp        ))  dTFr_dp         = EOS%dTFr_dp

end subroutine extract_member_EOS

end module MOM_EOS

!> \namespace mom_eos
!!
!! The MOM_EOS module is a wrapper for various equations of state (e.g. Linear,
!! Wright, UNESCO) and provides a uniform interface to the rest of the model
!! independent of which equation of state is being used.
