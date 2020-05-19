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

public calculate_compress, calculate_density, query_compressible
public calculate_density_derivs, calculate_specific_vol_derivs
public calculate_density_second_derivs
public EOS_init, EOS_manual_init, EOS_end, EOS_allocate, EOS_domain
public EOS_use_linear, calculate_spec_vol
public int_density_dz, int_specific_vol_dp
public int_density_dz_generic_plm, int_density_dz_generic_ppm
public int_spec_vol_dp_generic_plm !, int_spec_vol_dz_generic_ppm
public int_density_dz_generic, int_spec_vol_dp_generic
public find_depth_of_pressure_in_cell
public calculate_TFreeze
public convert_temp_salt_for_TEOS10
public gsw_sp_from_sr, gsw_pt_from_ct
public extract_member_EOS

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Calculates density of sea water from T, S and P
interface calculate_density
  module procedure calculate_density_scalar, calculate_density_array, calculate_density_1d
end interface calculate_density

!> Calculates specific volume of sea water from T, S and P
interface calculate_spec_vol
  module procedure calc_spec_vol_scalar, calculate_spec_vol_array, &
                   calc_spec_vol_1d
end interface calculate_spec_vol

!> Calculate the derivatives of density with temperature and salinity from T, S, and P
interface calculate_density_derivs
  module procedure calculate_density_derivs_scalar, calculate_density_derivs_array, &
                   calculate_density_derivs_1d
end interface calculate_density_derivs

!> Calculate the derivatives of specific volume with temperature and salinity from T, S, and P
interface calculate_specific_vol_derivs
  module procedure calculate_spec_vol_derivs_array, calc_spec_vol_derivs_1d
end interface calculate_specific_vol_derivs

!> Calculates the second derivatives of density with various combinations of temperature,
!! salinity, and pressure from T, S and P
interface calculate_density_second_derivs
  module procedure calculate_density_second_derivs_scalar, calculate_density_second_derivs_array
end interface calculate_density_second_derivs

!> Calculates the freezing point of sea water from T, S and P
interface calculate_TFreeze
  module procedure calculate_TFreeze_scalar, calculate_TFreeze_array
end interface calculate_TFreeze

!> Calculates the compressibility of water from T, S, and P
interface calculate_compress
  module procedure calculate_compress_scalar, calculate_compress_array
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
! change of units of arguments to functions)
  real :: m_to_Z = 1.      !< A constant that translates distances in meters to the units of depth.
  real :: kg_m3_to_R = 1.  !< A constant that translates kilograms per meter cubed to the units of density.
  real :: R_to_kg_m3 = 1.  !< A constant that translates the units of density to kilograms per meter cubed.
  real :: RL2_T2_to_Pa = 1.!< Convert pressures from R L2 T-2 to Pa.
  real :: L_T_to_m_s = 1.  !< Convert lateral velocities from L T-1 to m s-1.

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
!! density can be rescaled with the US.  If both the US and scale arguments are present the density
!! scaling uses the product of the two scaling factors.
subroutine calculate_density_scalar(T, S, pressure, rho, EOS, rho_ref, scale)
  real,           intent(in)  :: T        !< Potential temperature referenced to the surface [degC]
  real,           intent(in)  :: S        !< Salinity [ppt]
  real,           intent(in)  :: pressure !< Pressure [Pa] or [R L2 T-2 ~> Pa]
  real,           intent(out) :: rho      !< Density (in-situ if pressure is local) [kg m-3] or [R ~> kg m-3]
  type(EOS_type), pointer     :: EOS      !< Equation of state structure
  real, optional, intent(in)  :: rho_ref  !< A reference density [kg m-3]
  real, optional, intent(in)  :: scale    !< A multiplicative factor by which to scale density in
                                          !! combination with scaling given by US [various]

  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real :: p_scale   ! A factor to convert pressure to units of Pa [Pa T2 R-1 L-2 ~> 1]

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_density_scalar called with an unassociated EOS_type EOS.")

  p_scale = EOS%RL2_T2_to_Pa

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_density_linear(T, S, p_scale*pressure, rho, &
                                      EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS, rho_ref)
    case (EOS_UNESCO)
      call calculate_density_unesco(T, S, p_scale*pressure, rho, rho_ref)
    case (EOS_WRIGHT)
      call calculate_density_wright(T, S, p_scale*pressure, rho, rho_ref)
    case (EOS_TEOS10)
      call calculate_density_teos10(T, S, p_scale*pressure, rho, rho_ref)
    case (EOS_NEMO)
      call calculate_density_nemo(T, S, p_scale*pressure, rho, rho_ref)
    case default
      call MOM_error(FATAL, "calculate_density_scalar: EOS is not valid.")
  end select

  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  rho = rho_scale * rho

end subroutine calculate_density_scalar

!> Calls the appropriate subroutine to calculate the density of sea water for 1-D array inputs.
!! If rho_ref is present, the anomaly with respect to rho_ref is returned.
subroutine calculate_density_array(T, S, pressure, rho, start, npts, EOS, rho_ref, scale)
  real, dimension(:), intent(in)    :: T        !< Potential temperature referenced to the surface [degC]
  real, dimension(:), intent(in)    :: S        !< Salinity [ppt]
  real, dimension(:), intent(in)    :: pressure !< Pressure [Pa] or [R L2 T-2 ~> Pa]
  real, dimension(:), intent(inout) :: rho      !< Density (in-situ if pressure is local) [kg m-3] or [R ~> kg m-3]
  integer,            intent(in)    :: start    !< Start index for computation
  integer,            intent(in)    :: npts     !< Number of point to compute
  type(EOS_type),     pointer       :: EOS      !< Equation of state structure
  real,                  optional, intent(in) :: rho_ref  !< A reference density [kg m-3]
  real,                  optional, intent(in) :: scale    !< A multiplicative factor by which to scale density
                                                !! in combination with scaling given by US [various]

  integer :: j

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_density_array called with an unassociated EOS_type EOS.")

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

!> Calls the appropriate subroutine to calculate the density of sea water for 1-D array inputs,
!! potentially limiting the domain of indices that are worked on.
!! If rho_ref is present, the anomaly with respect to rho_ref is returned.
subroutine calculate_density_1d(T, S, pressure, rho, EOS, dom, rho_ref, scale)
  real, dimension(:),    intent(in)    :: T        !< Potential temperature referenced to the surface [degC]
  real, dimension(:),    intent(in)    :: S        !< Salinity [ppt]
  real, dimension(:),    intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:),    intent(inout) :: rho      !< Density (in-situ if pressure is local) [R ~> kg m-3]
  type(EOS_type),        pointer       :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.
  real,                  optional, intent(in) :: rho_ref !< A reference density [kg m-3]
  real,                  optional, intent(in) :: scale !< A multiplicative factor by which to scale density
                                                   !! in combination with scaling given by US [various]
  ! Local variables
  real :: p_scale   ! A factor to convert pressure to units of Pa [Pa T2 R-1 L-2 ~> 1]
  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real :: rho_unscale ! A factor to convert density from R to kg m-3 [kg m-3 R-1 ~> 1]
  real :: rho_reference ! rho_ref converted to [kg m-3]
  real, dimension(size(rho)) :: pres  ! Pressure converted to [Pa]
  integer :: i, is, ie, npts

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_density_1d called with an unassociated EOS_type EOS.")

  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(rho) ; npts = 1 + ie - is
  endif

  p_scale = EOS%RL2_T2_to_Pa
  rho_unscale = EOS%R_to_kg_m3

  if ((p_scale == 1.0) .and. (rho_unscale == 1.0)) then
    call calculate_density_array(T, S, pressure, rho, is, npts, EOS, rho_ref=rho_ref)
  elseif (present(rho_ref)) then ! This is the same as above, but with some extra work to rescale variables.
    do i=is,ie ; pres(i) = p_scale * pressure(i) ; enddo
    rho_reference = rho_unscale*rho_ref
    call calculate_density_array(T, S, pres, rho, is, npts, EOS, rho_ref=rho_reference)
  else  ! There is rescaling of variables, but rho_ref is not present. Passing a 0 value of rho_ref
        ! changes answers at roundoff for some equations of state, like Wright and UNESCO.
    do i=is,ie ; pres(i) = p_scale * pressure(i) ; enddo
    call calculate_density_array(T, S, pres, rho, is, npts, EOS)
  endif

  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  if (rho_scale /= 1.0) then ; do i=is,ie
    rho(i) = rho_scale * rho(i)
  enddo ; endif

end subroutine calculate_density_1d

!> Calls the appropriate subroutine to calculate the specific volume of sea water
!! for 1-D array inputs.
subroutine calculate_spec_vol_array(T, S, pressure, specvol, start, npts, EOS, spv_ref, scale)
  real, dimension(:), intent(in)    :: T        !< potential temperature relative to the surface [degC]
  real, dimension(:), intent(in)    :: S        !< salinity [ppt]
  real, dimension(:), intent(in)    :: pressure !< pressure [Pa]
  real, dimension(:), intent(inout) :: specvol  !< in situ specific volume [kg m-3]
  integer,            intent(in)    :: start    !< the starting point in the arrays.
  integer,            intent(in)    :: npts     !< the number of values to calculate.
  type(EOS_type),     pointer       :: EOS      !< Equation of state structure
  real,     optional, intent(in)    :: spv_ref  !< A reference specific volume [m3 kg-1]
  real,     optional, intent(in)    :: scale    !< A multiplicative factor by which to scale specific
                                                !! volume in combination with scaling given by US [various]

  real, dimension(size(specvol))  :: rho   ! Density [kg m-3]
  integer :: j

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_spec_vol_array called with an unassociated EOS_type EOS.")

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
  real,           intent(in)  :: T        !< Potential temperature referenced to the surface [degC]
  real,           intent(in)  :: S        !< Salinity [ppt]
  real,           intent(in)  :: pressure !< Pressure [Pa] or [R L2 T-2 ~> Pa]
  real,           intent(out) :: specvol  !< In situ? specific volume [m3 kg-1] or [R-1 ~> m3 kg-1]
  type(EOS_type), pointer     :: EOS      !< Equation of state structure
  real, optional, intent(in)  :: spv_ref  !< A reference specific volume [m3 kg-1] or [R-1 m3 kg-1]
  real, optional, intent(in)  :: scale    !< A multiplicative factor by which to scale specific
                                          !! volume in combination with scaling given by US [various]

  real, dimension(1) :: Ta, Sa, pres, spv  ! Rescaled single element array versions of the arguments.
  real :: spv_reference ! spv_ref converted to [m3 kg-1]
  real :: spv_scale ! A factor to convert specific volume from m3 kg-1 to the desired units [kg R-1 m-3 ~> 1]

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calc_spec_vol_scalar called with an unassociated EOS_type EOS.")

  pres(1) = EOS%RL2_T2_to_Pa*pressure
  Ta(1) = T ; Sa(1) = S

  if (present(spv_ref)) then
    spv_reference = EOS%kg_m3_to_R*spv_ref
    call calculate_spec_vol_array(Ta, Sa, pres, spv, 1, 1, EOS, spv_reference)
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
  real, dimension(:),    intent(in)    :: T        !< Potential temperature referenced to the surface [degC]
  real, dimension(:),    intent(in)    :: S        !< Salinity [ppt]
  real, dimension(:),    intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:),    intent(inout) :: specvol  !< In situ specific volume [R-1 ~> m3 kg-1]
  type(EOS_type),        pointer       :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.
  real,                  optional, intent(in) :: spv_ref !< A reference specific volume [R-1 ~> m3 kg-1]
  real,                  optional, intent(in) :: scale !< A multiplicative factor by which to scale
                                                       !! output specific volume in combination with
                                                       !! scaling given by US [various]
  ! Local variables
  real, dimension(size(specvol)) :: pres  ! Pressure converted to [Pa]
  real :: p_scale   ! A factor to convert pressure to units of Pa [Pa T2 R-1 L-2 ~> 1]
  real :: spv_unscale ! A factor to convert specific volume from R-1 to m3 kg-1 [m3 kg-1 R ~> 1]
  real :: spv_scale ! A factor to convert specific volume from m3 kg-1 to the desired units [kg m-3 R-1 ~> 1]
  real :: spv_reference ! spv_ref converted to [m3 kg-1]
  integer :: i, is, ie, npts

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calc_spec_vol_1d called with an unassociated EOS_type EOS.")

  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(specvol) ; npts = 1 + ie - is
  endif

  p_scale = EOS%RL2_T2_to_Pa
  spv_unscale = EOS%kg_m3_to_R

  if ((p_scale == 1.0) .and. (spv_unscale == 1.0)) then
    call calculate_spec_vol_array(T, S, pressure, specvol, is, npts, EOS, spv_ref)
  elseif (present(spv_ref)) then ! This is the same as above, but with some extra work to rescale variables.
    do i=is,ie ; pres(i) = p_scale * pressure(i) ; enddo
    spv_reference = spv_unscale*spv_ref
    call calculate_spec_vol_array(T, S, pres, specvol, is, npts, EOS, spv_reference)
  else  ! There is rescaling of variables, but spv_ref is not present. Passing a 0 value of spv_ref
        ! changes answers at roundoff for some equations of state, like Wright and UNESCO.
    do i=is,ie ; pres(i) = p_scale * pressure(i) ; enddo
    call calculate_spec_vol_array(T, S, pres, specvol, is, npts, EOS)
  endif

  spv_scale = EOS%R_to_kg_m3
  if (present(scale)) spv_scale = spv_scale * scale
  if (spv_scale /= 1.0) then ; do i=is,ie
    specvol(i) = spv_scale * specvol(i)
  enddo ; endif

end subroutine calc_spec_vol_1d


!> Calls the appropriate subroutine to calculate the freezing point for scalar inputs.
subroutine calculate_TFreeze_scalar(S, pressure, T_fr, EOS, pres_scale)
  real,           intent(in)  :: S !< Salinity [ppt]
  real,           intent(in)  :: pressure !< Pressure [Pa] or [other]
  real,           intent(out) :: T_fr !< Freezing point potential temperature referenced
                                      !! to the surface [degC]
  type(EOS_type), pointer     :: EOS !< Equation of state structure
  real, optional, intent(in)  :: pres_scale !< A multiplicative factor to convert pressure into Pa

  ! Local variables
  real :: p_scale ! A factor to convert pressure to units of Pa.

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_TFreeze_scalar called with an unassociated EOS_type EOS.")

  p_scale = 1.0 ; if (present(pres_scale)) p_scale = pres_scale

  select case (EOS%form_of_TFreeze)
    case (TFREEZE_LINEAR)
      call calculate_TFreeze_linear(S, p_scale*pressure, T_fr, EOS%TFr_S0_P0, &
                                    EOS%dTFr_dS, EOS%dTFr_dp)
    case (TFREEZE_MILLERO)
      call calculate_TFreeze_Millero(S, p_scale*pressure, T_fr)
    case (TFREEZE_TEOS10)
      call calculate_TFreeze_teos10(S, p_scale*pressure, T_fr)
    case default
      call MOM_error(FATAL, "calculate_TFreeze_scalar: form_of_TFreeze is not valid.")
  end select

end subroutine calculate_TFreeze_scalar

!> Calls the appropriate subroutine to calculate the freezing point for a 1-D array.
subroutine calculate_TFreeze_array(S, pressure, T_fr, start, npts, EOS, pres_scale)
  real, dimension(:), intent(in)    :: S        !< Salinity [ppt]
  real, dimension(:), intent(in)    :: pressure !< Pressure [Pa] or [other]
  real, dimension(:), intent(inout) :: T_fr     !< Freezing point potential temperature referenced
                                                !! to the surface [degC]
  integer,            intent(in)    :: start    !< Starting index within the array
  integer,            intent(in)    :: npts     !< The number of values to calculate
  type(EOS_type),     pointer       :: EOS      !< Equation of state structure
  real,     optional, intent(in)    :: pres_scale !< A multiplicative factor to convert pressure into Pa.

  ! Local variables
  real, dimension(size(pressure)) :: pres  ! Pressure converted to [Pa]
  real :: p_scale ! A factor to convert pressure to units of Pa.
  integer :: j

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_TFreeze_scalar called with an unassociated EOS_type EOS.")

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

!> Calls the appropriate subroutine to calculate density derivatives for 1-D array inputs.
subroutine calculate_density_derivs_array(T, S, pressure, drho_dT, drho_dS, start, npts, EOS, scale)
  real, dimension(:), intent(in)    :: T        !< Potential temperature referenced to the surface [degC]
  real, dimension(:), intent(in)    :: S        !< Salinity [ppt]
  real, dimension(:), intent(in)    :: pressure !< Pressure [Pa] or [R L2 T-2 ~> Pa]
  real, dimension(:), intent(inout) :: drho_dT  !< The partial derivative of density with potential
                                                !! temperature [kg m-3 degC-1] or [R degC-1 ~> kg m-3 degC-1]
  real, dimension(:), intent(inout) :: drho_dS  !< The partial derivative of density with salinity,
                                                !! in [kg m-3 ppt-1] or [R degC-1 ~> kg m-3 ppt-1]
  integer,            intent(in)    :: start    !< Starting index within the array
  integer,            intent(in)    :: npts     !< The number of values to calculate
  type(EOS_type),     pointer       :: EOS      !< Equation of state structure
  real,     optional, intent(in)    :: scale !< A multiplicative factor by which to scale density
                                                !! in combination with scaling given by US [various]

  ! Local variables
  integer :: j

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_density_derivs called with an unassociated EOS_type EOS.")

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
  real, dimension(:),    intent(in)    :: T        !< Potential temperature referenced to the surface [degC]
  real, dimension(:),    intent(in)    :: S        !< Salinity [ppt]
  real, dimension(:),    intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:),    intent(inout) :: drho_dT  !< The partial derivative of density with potential
                                                   !! temperature [R degC-1 ~> kg m-3 degC-1]
  real, dimension(:),    intent(inout) :: drho_dS  !< The partial derivative of density with salinity
                                                   !! [R degC-1 ~> kg m-3 ppt-1]
  type(EOS_type),        pointer       :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.
  real,                  optional, intent(in) :: scale !< A multiplicative factor by which to scale density
                                                       !! in combination with scaling given by US [various]
  ! Local variables
  real, dimension(size(drho_dT)) :: pres  ! Pressure converted to [Pa]
  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real :: p_scale   ! A factor to convert pressure to units of Pa [Pa T2 R-1 L-2 ~> 1]
  integer :: i, is, ie, npts

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_density_derivs called with an unassociated EOS_type EOS.")

  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(drho_dT) ; npts = 1 + ie - is
  endif

  p_scale = EOS%RL2_T2_to_Pa

  if (p_scale == 1.0) then
    call calculate_density_derivs_array(T, S, pressure, drho_dT, drho_dS, is, npts, EOS)
  else
    do i=is,ie ; pres(i) = p_scale * pressure(i) ; enddo
    call calculate_density_derivs_array(T, S, pres, drho_dT, drho_dS, is, npts, EOS)
  endif

  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  if (rho_scale /= 1.0) then ; do i=is,ie
    drho_dT(i) = rho_scale * drho_dT(i)
    drho_dS(i) = rho_scale * drho_dS(i)
  enddo ; endif

end subroutine calculate_density_derivs_1d


!> Calls the appropriate subroutines to calculate density derivatives by promoting a scalar
!! to a one-element array
subroutine calculate_density_derivs_scalar(T, S, pressure, drho_dT, drho_dS, EOS, scale)
  real,           intent(in)  :: T !< Potential temperature referenced to the surface [degC]
  real,           intent(in)  :: S !< Salinity [ppt]
  real,           intent(in)  :: pressure !< Pressure [Pa] or [R L2 T-2 ~> Pa]
  real,           intent(out) :: drho_dT !< The partial derivative of density with potential
                                         !! temperature [kg m-3 degC-1] or [R degC-1 ~> kg m-3 degC-1]
  real,           intent(out) :: drho_dS !< The partial derivative of density with salinity,
                                         !! in [kg m-3 ppt-1] or [R ppt-1 ~> kg m-3 ppt-1]
  type(EOS_type), pointer     :: EOS     !< Equation of state structure
  real, optional, intent(in)  :: scale   !< A multiplicative factor by which to scale density
                                         !! in combination with scaling given by US [various]
  ! Local variables
  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real :: p_scale   ! A factor to convert pressure to units of Pa [Pa T2 R-1 L-2 ~> 1]
  integer :: j

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_density_derivs called with an unassociated EOS_type EOS.")

  p_scale = EOS%RL2_T2_to_Pa

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_density_derivs_linear(T, S, p_scale*pressure, drho_dT, drho_dS, &
                                           EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS)
    case (EOS_WRIGHT)
      call calculate_density_derivs_wright(T, S, p_scale*pressure, drho_dT, drho_dS)
    case (EOS_TEOS10)
      call calculate_density_derivs_teos10(T, S, p_scale*pressure, drho_dT, drho_dS)
    case default
      call MOM_error(FATAL, "calculate_density_derivs_scalar: EOS%form_of_EOS is not valid.")
  end select

  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  if (rho_scale /= 1.0) then
    drho_dT = rho_scale * drho_dT
    drho_dS = rho_scale * drho_dS
  endif

end subroutine calculate_density_derivs_scalar

!> Calls the appropriate subroutine to calculate density second derivatives for 1-D array inputs.
subroutine calculate_density_second_derivs_array(T, S, pressure, drho_dS_dS, drho_dS_dT, drho_dT_dT, &
                                                 drho_dS_dP, drho_dT_dP, start, npts, EOS, scale)
  real, dimension(:), intent(in)  :: T !< Potential temperature referenced to the surface [degC]
  real, dimension(:), intent(in)  :: S !< Salinity [ppt]
  real, dimension(:), intent(in)  :: pressure   !< Pressure [Pa] or [R L2 T-2 ~> Pa]
  real, dimension(:), intent(inout) :: drho_dS_dS !< Partial derivative of beta with respect to S
                                                  !!  [kg m-3 ppt-2] or [R ppt-2 ~> kg m-3 ppt-2]
  real, dimension(:), intent(inout) :: drho_dS_dT !< Partial derivative of beta with respect to T
                                                  !! [kg m-3 ppt-1 degC-1] or [R ppt-1 degC-1 ~> kg m-3 ppt-1 degC-1]
  real, dimension(:), intent(inout) :: drho_dT_dT !< Partial derivative of alpha with respect to T
                                                  !! [kg m-3 degC-2] or [R degC-2 ~> kg m-3 degC-2]
  real, dimension(:), intent(inout) :: drho_dS_dP !< Partial derivative of beta with respect to pressure
                                                  !! [kg m-3 ppt-1 Pa-1] or [R ppt-1 Pa-1 ~> kg m-3 ppt-1 Pa-1]
  real, dimension(:), intent(inout) :: drho_dT_dP !< Partial derivative of alpha with respect to pressure
                                                  !! [kg m-3 degC-1 Pa-1] or [R degC-1 Pa-1 ~> kg m-3 degC-1 Pa-1]
  integer,            intent(in)  :: start !< Starting index within the array
  integer,            intent(in)  :: npts  !< The number of values to calculate
  type(EOS_type),     pointer     :: EOS   !< Equation of state structure
  real,                  optional, intent(in) :: scale !< A multiplicative factor by which to scale density
                                                  !! in combination with scaling given by US [various]
  ! Local variables
  real, dimension(size(pressure)) :: pres  ! Pressure converted to [Pa]
  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real :: p_scale   ! A factor to convert pressure to units of Pa [Pa T2 R-1 L-2 ~> 1]
  real :: I_p_scale ! The inverse of the factor to convert pressure to units of Pa [R L2 T-2 Pa-1 ~> 1]
  integer :: j

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_density_derivs called with an unassociated EOS_type EOS.")

  p_scale = EOS%RL2_T2_to_Pa

  if (p_scale == 1.0) then
    select case (EOS%form_of_EOS)
      case (EOS_LINEAR)
        call calculate_density_second_derivs_linear(T, S, pressure, drho_dS_dS, drho_dS_dT, &
                                                    drho_dT_dT, drho_dS_dP, drho_dT_dP, start, npts)
      case (EOS_WRIGHT)
        call calculate_density_second_derivs_wright(T, S, pressure, drho_dS_dS, drho_dS_dT, &
                                                    drho_dT_dT, drho_dS_dP, drho_dT_dP, start, npts)
      case (EOS_TEOS10)
        call calculate_density_second_derivs_teos10(T, S, pressure, drho_dS_dS, drho_dS_dT, &
                                                    drho_dT_dT, drho_dS_dP, drho_dT_dP, start, npts)
      case default
        call MOM_error(FATAL, "calculate_density_derivs: EOS%form_of_EOS is not valid.")
    end select
  else
    do j=start,start+npts-1 ; pres(j) = p_scale * pressure(j) ; enddo
    select case (EOS%form_of_EOS)
      case (EOS_LINEAR)
        call calculate_density_second_derivs_linear(T, S, pres, drho_dS_dS, drho_dS_dT, &
                                                    drho_dT_dT, drho_dS_dP, drho_dT_dP, start, npts)
      case (EOS_WRIGHT)
        call calculate_density_second_derivs_wright(T, S, pres, drho_dS_dS, drho_dS_dT, &
                                                    drho_dT_dT, drho_dS_dP, drho_dT_dP, start, npts)
      case (EOS_TEOS10)
        call calculate_density_second_derivs_teos10(T, S, pres, drho_dS_dS, drho_dS_dT, &
                                                    drho_dT_dT, drho_dS_dP, drho_dT_dP, start, npts)
      case default
        call MOM_error(FATAL, "calculate_density_derivs: EOS%form_of_EOS is not valid.")
    end select
  endif

  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  if (rho_scale /= 1.0) then ; do j=start,start+npts-1
    drho_dS_dS(j) = rho_scale * drho_dS_dS(j)
    drho_dS_dT(j) = rho_scale * drho_dS_dT(j)
    drho_dT_dT(j) = rho_scale * drho_dT_dT(j)
    drho_dS_dP(j) = rho_scale * drho_dS_dP(j)
    drho_dT_dP(j) = rho_scale * drho_dT_dP(j)
  enddo ; endif

  if (p_scale /= 1.0) then
    I_p_scale = 1.0 / p_scale
    do j=start,start+npts-1
      drho_dS_dP(j) = I_p_scale * drho_dS_dP(j)
      drho_dT_dP(j) = I_p_scale * drho_dT_dP(j)
    enddo
  endif

end subroutine calculate_density_second_derivs_array

!> Calls the appropriate subroutine to calculate density second derivatives for scalar nputs.
subroutine calculate_density_second_derivs_scalar(T, S, pressure, drho_dS_dS, drho_dS_dT, drho_dT_dT, &
                                                  drho_dS_dP, drho_dT_dP, EOS, scale)
  real, intent(in)  :: T !< Potential temperature referenced to the surface [degC]
  real, intent(in)  :: S !< Salinity [ppt]
  real, intent(in)  :: pressure   !< Pressure [Pa] or [R L2 T-2 ~> Pa]
  real, intent(out) :: drho_dS_dS !< Partial derivative of beta with respect to S
                                  !! [kg m-3 ppt-2] or [R ppt-2 ~> kg m-3 ppt-2]
  real, intent(out) :: drho_dS_dT !< Partial derivative of beta with respect to T
                                  !! [kg m-3 ppt-1 degC-1] or [R ppt-1 degC-1 ~> kg m-3 ppt-1 degC-1]
  real, intent(out) :: drho_dT_dT !< Partial derivative of alpha with respect to T
                                  !! [kg m-3 degC-2] or [R degC-2 ~> kg m-3 degC-2]
  real, intent(out) :: drho_dS_dP !< Partial derivative of beta with respect to pressure
                                  !! [kg m-3 ppt-1 Pa-1] or [R ppt-1 Pa-1 ~> kg m-3 ppt-1 Pa-1]
  real, intent(out) :: drho_dT_dP !< Partial derivative of alpha with respect to pressure
                                  !! [kg m-3 degC-1 Pa-1] or [R degC-1 Pa-1 ~> kg m-3 degC-1 Pa-1]
  type(EOS_type), pointer    :: EOS !< Equation of state structure
  real, optional, intent(in) :: scale !< A multiplicative factor by which to scale density
                                  !! in combination with scaling given by US [various]
  ! Local variables
  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real :: p_scale   ! A factor to convert pressure to units of Pa [Pa T2 R-1 L-2 ~> 1]
  real :: I_p_scale ! The inverse of the factor to convert pressure to units of Pa [R L2 T-2 Pa-1 ~> 1]

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_density_derivs called with an unassociated EOS_type EOS.")

  p_scale = EOS%RL2_T2_to_Pa

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_density_second_derivs_linear(T, S, p_scale*pressure, drho_dS_dS, drho_dS_dT, &
                                                  drho_dT_dT, drho_dS_dP, drho_dT_dP)
    case (EOS_WRIGHT)
      call calculate_density_second_derivs_wright(T, S, p_scale*pressure, drho_dS_dS, drho_dS_dT, &
                                                  drho_dT_dT, drho_dS_dP, drho_dT_dP)
    case (EOS_TEOS10)
      call calculate_density_second_derivs_teos10(T, S, p_scale*pressure, drho_dS_dS, drho_dS_dT, &
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

  if (p_scale /= 1.0) then
    I_p_scale = 1.0 / p_scale
    drho_dS_dP = I_p_scale * drho_dS_dP
    drho_dT_dP = I_p_scale * drho_dT_dP
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
  type(EOS_type),     pointer     :: EOS    !< Equation of state structure

  ! Local variables
  real, dimension(size(T)) :: press   ! Pressure converted to [Pa]
  real, dimension(size(T)) :: rho     ! In situ density [kg m-3]
  real, dimension(size(T)) :: dRho_dT ! Derivative of density with temperature [kg m-3 degC-1]
  real, dimension(size(T)) :: dRho_dS ! Derivative of density with salinity [kg m-3 ppt-1]
  integer :: j

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_spec_vol_derivs_array called with an unassociated EOS_type EOS.")

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
  real, dimension(:), intent(in)    :: T        !< Potential temperature referenced to the surface [degC]
  real, dimension(:), intent(in)    :: S        !< Salinity [ppt]
  real, dimension(:), intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:), intent(inout) :: dSV_dT   !< The partial derivative of specific volume with potential
                                                !! temperature [R-1 degC-1 ~> m3 kg-1 degC-1]
  real, dimension(:), intent(inout) :: dSV_dS   !< The partial derivative of specific volume with salinity
                                                !! [R-1 ppt-1 ~> m3 kg-1 ppt-1]
  type(EOS_type),     pointer       :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.
  real,                  optional, intent(in) :: scale !< A multiplicative factor by which to scale specific
                                                !! volume in combination with scaling given by US [various]

  ! Local variables
  real, dimension(size(dSV_dT)) :: press   ! Pressure converted to [Pa]
  real :: spv_scale ! A factor to convert specific volume from m3 kg-1 to the desired units [kg R-1 m-3 ~> 1]
  real :: p_scale   ! A factor to convert pressure to units of Pa [Pa T2 R-1 L-2 ~> 1]
  integer :: i, is, ie, npts

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_spec_vol_derivs_1d called with an unassociated EOS_type EOS.")

  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(dSV_dT) ; npts = 1 + ie - is
  endif
  p_scale = EOS%RL2_T2_to_Pa

  if (p_scale == 1.0) then
    call calculate_spec_vol_derivs_array(T, S, pressure, dSV_dT, dSV_dS, is, npts, EOS)
  else
    do i=is,ie ; press(i) = p_scale * pressure(i) ; enddo
    call calculate_spec_vol_derivs_array(T, S, press, dSV_dT, dSV_dS, is, npts, EOS)
  endif

  spv_scale = EOS%R_to_kg_m3
  if (present(scale)) spv_scale = spv_scale * scale
  if (spv_scale /= 1.0) then ; do i=is,ie
    dSV_dT(i) = spv_scale * dSV_dT(i)
    dSV_dS(i) = spv_scale * dSV_dS(i)
  enddo ; endif

end subroutine calc_spec_vol_derivs_1d


!> Calls the appropriate subroutine to calculate the density and compressibility for 1-D array
!! inputs.  If US is present, the units of the inputs and outputs are rescaled.
subroutine calculate_compress_array(T, S, press, rho, drho_dp, start, npts, EOS)
  real, dimension(:), intent(in)  :: T        !< Potential temperature referenced to the surface [degC]
  real, dimension(:), intent(in)  :: S        !< Salinity [PSU]
  real, dimension(:), intent(in)  :: press    !< Pressure [Pa] or [R L2 T-2 ~> Pa]
  real, dimension(:), intent(inout) :: rho      !< In situ density [kg m-3] or [R ~> kg m-3]
  real, dimension(:), intent(inout) :: drho_dp  !< The partial derivative of density with pressure
                                                !! (also the inverse of the square of sound speed)
                                                !! [s2 m-2] or [T2 L-2]
  integer,            intent(in)  :: start    !< Starting index within the array
  integer,            intent(in)  :: npts     !< The number of values to calculate
  type(EOS_type),     pointer     :: EOS      !< Equation of state structure

  ! Local variables
  real, dimension(size(press)) :: pressure  ! Pressure converted to [Pa]
  integer :: i, is, ie

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_compress called with an unassociated EOS_type EOS.")

  is = start ; ie = is + npts - 1
  do i=is,ie ; pressure(i) = EOS%RL2_T2_to_Pa * press(i) ; enddo

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_compress_linear(T, S, pressure, rho, drho_dp, start, npts, &
                                     EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS)
    case (EOS_UNESCO)
      call calculate_compress_unesco(T, S, pressure, rho, drho_dp, start, npts)
    case (EOS_WRIGHT)
      call calculate_compress_wright(T, S, pressure, rho, drho_dp, start, npts)
    case (EOS_TEOS10)
      call calculate_compress_teos10(T, S, pressure, rho, drho_dp, start, npts)
    case (EOS_NEMO)
      call calculate_compress_nemo(T, S, pressure, rho, drho_dp, start, npts)
    case default
      call MOM_error(FATAL, "calculate_compress: EOS%form_of_EOS is not valid.")
  end select

  if (EOS%kg_m3_to_R /= 1.0) then ; do i=is,ie
    rho(i) = EOS%kg_m3_to_R * rho(i)
  enddo ; endif
  if (EOS%L_T_to_m_s /= 1.0) then ; do i=is,ie
    drho_dp(i) = EOS%L_T_to_m_s**2 * drho_dp(i)
  enddo ; endif

end subroutine calculate_compress_array

!> Calculate density and compressibility for a scalar. This just promotes the scalar to an array
!! with a singleton dimension and calls calculate_compress_array.  If US is present, the units of
!! the inputs and outputs are rescaled.
subroutine calculate_compress_scalar(T, S, pressure, rho, drho_dp, EOS)
  real, intent(in)        :: T        !< Potential temperature referenced to the surface [degC]
  real, intent(in)        :: S        !< Salinity [ppt]
  real, intent(in)        :: pressure !< Pressure [Pa] or [R L2 T-2 ~> Pa]
  real, intent(out)       :: rho      !< In situ density [kg m-3] or [R ~> kg m-3]
  real, intent(out)       :: drho_dp  !< The partial derivative of density with pressure (also the
                                      !! inverse of the square of sound speed) [s2 m-2] or [T2 L-2]
  type(EOS_type), pointer :: EOS      !< Equation of state structure

  ! Local variables
  real, dimension(1) :: Ta, Sa, pa, rhoa, drho_dpa

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_compress called with an unassociated EOS_type EOS.")
  Ta(1) = T ; Sa(1) = S; pa(1) = pressure

  call calculate_compress_array(Ta, Sa, pa, rhoa, drho_dpa, 1, 1, EOS)
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
!! use of Bode's rule to do the horizontal integrals, and from a truncation in the
!! series for log(1-eps/1+eps) that assumes that |eps| < 0.34.
subroutine int_specific_vol_dp(T, S, p_t, p_b, alpha_ref, HI, EOS, &
                               dza, intp_dza, intx_dza, inty_dza, halo_size, &
                               bathyP, dP_tiny, useMassWghtInterp)
  type(hor_index_type), intent(in)  :: HI  !< The horizontal index structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T   !< Potential temperature referenced to the surface [degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S   !< Salinity [ppt]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: p_t !< Pressure at the top of the layer [R L2 T-2 ~> Pa] or [Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: p_b !< Pressure at the bottom of the layer [R L2 T-2 ~> Pa] or [Pa]
  real,                 intent(in)  :: alpha_ref !< A mean specific volume that is subtracted out
                            !! to reduce the magnitude of each of the integrals [R-1 ~> m3 kg-1]
                            !! The calculation is mathematically identical with different values of
                            !! alpha_ref, but this reduces the effects of roundoff.
  type(EOS_type),       pointer     :: EOS !< Equation of state structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(inout) :: dza !< The change in the geopotential anomaly across
                            !! the layer [L2 T-2 ~> m2 s-2] or [m2 s-2]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(inout) :: intp_dza !< The integral in pressure through the layer of the
                            !! geopotential anomaly relative to the anomaly at the bottom of the
                            !! layer [R L4 T-4 ~> Pa m2 s-2] or [Pa m2 s-2]
  real, dimension(HI%IsdB:HI%IedB,HI%jsd:HI%jed), &
              optional, intent(inout) :: intx_dza !< The integral in x of the difference between the
                            !! geopotential anomaly at the top and bottom of the layer divided by
                            !! the x grid spacing [L2 T-2 ~> m2 s-2] or [m2 s-2]
  real, dimension(HI%isd:HI%ied,HI%JsdB:HI%JedB), &
              optional, intent(inout) :: inty_dza !< The integral in y of the difference between the
                            !! geopotential anomaly at the top and bottom of the layer divided by
                            !! the y grid spacing [L2 T-2 ~> m2 s-2] or [m2 s-2]
  integer,    optional, intent(in)  :: halo_size !< The width of halo points on which to calculate dza.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: bathyP  !< The pressure at the bathymetry [R L2 T-2 ~> Pa] or [Pa]
  real,       optional, intent(in)  :: dP_tiny !< A miniscule pressure change with
                            !! the same units as p_t [R L2 T-2 ~> Pa] or [Pa]
  logical,    optional, intent(in)  :: useMassWghtInterp !< If true, uses mass weighting
                            !! to interpolate T/S for top and bottom integrals.
  ! Local variables
  real :: pres_scale    ! A unit conversion factor from the rescaled units of pressure to Pa [Pa T2 R-1 L-2 ~> 1]
  real :: SV_scale      ! A multiplicative factor by which to scale specific
                        ! volume from m3 kg-1 to the desired units [kg m-3 R-1 ~> 1]

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "int_specific_vol_dp called with an unassociated EOS_type EOS.")

  if (EOS%EOS_quadrature) then
    call int_spec_vol_dp_generic(T, S, p_t, p_b, alpha_ref, HI, EOS, &
                                 dza, intp_dza, intx_dza, inty_dza, halo_size, &
                                 bathyP, dP_tiny, useMassWghtInterp)
  else ; select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call int_spec_vol_dp_linear(T, S, p_t, p_b, alpha_ref, HI, EOS%kg_m3_to_R*EOS%Rho_T0_S0, &
                                EOS%kg_m3_to_R*EOS%dRho_dT, EOS%kg_m3_to_R*EOS%dRho_dS, dza, &
                                intp_dza, intx_dza, inty_dza, halo_size, &
                                bathyP, dP_tiny, useMassWghtInterp)
    case (EOS_WRIGHT)
      call int_spec_vol_dp_wright(T, S, p_t, p_b, alpha_ref, HI, dza, intp_dza, intx_dza, &
                                  inty_dza, halo_size, bathyP, dP_tiny, useMassWghtInterp, &
                                  SV_scale=EOS%R_to_kg_m3, pres_scale=EOS%RL2_T2_to_Pa)
    case default
      call int_spec_vol_dp_generic(T, S, p_t, p_b, alpha_ref, HI, EOS, &
                                   dza, intp_dza, intx_dza, inty_dza, halo_size, &
                                   bathyP, dP_tiny, useMassWghtInterp)
  end select ; endif

end subroutine int_specific_vol_dp

!> This subroutine calculates analytical and nearly-analytical integrals of
!! pressure anomalies across layers, which are required for calculating the
!! finite-volume form pressure accelerations in a Boussinesq model.
subroutine int_density_dz(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, EOS, dpa, &
                          intz_dpa, intx_dpa, inty_dpa, bathyT, dz_neglect, useMassWghtInterp)
  type(hor_index_type), intent(in)  :: HI !< Ocean horizontal index structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T   !< Potential temperature referenced to the surface [degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S   !< Salinity [ppt]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_t !< Height at the top of the layer in depth units [Z ~> m]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_b !< Height at the bottom of the layer [Z ~> m]
  real,                 intent(in)  :: rho_ref !< A mean density [R ~> kg m-3] or [kg m-3], that is
                                           !! subtracted out to reduce the magnitude of each of the
                                           !! integrals.
  real,                 intent(in)  :: rho_0 !< A density [R ~> kg m-3] or [kg m-3], that is used
                                           !! to calculate the pressure (as p~=-z*rho_0*G_e)
                                           !! used in the equation of state.
  real,                 intent(in)  :: G_e !< The Earth's gravitational acceleration
                                           !! [L2 Z-1 T-2 ~> m s-2] or [m2 Z-1 s-2 ~> m s-2]
  type(EOS_type),       pointer     :: EOS !< Equation of state structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                      intent(inout) :: dpa !< The change in the pressure anomaly
                                           !! across the layer [R L2 T-2 ~> Pa] or [Pa]
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
  ! Local variables
  real :: rho_scale  ! A multiplicative factor by which to scale density from kg m-3 to the
                     ! desired units [R m3 kg-1 ~> 1]
  real :: pres_scale ! A multiplicative factor to convert pressure into Pa [Pa T2 R-1 L-2 ~> 1]

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "int_density_dz called with an unassociated EOS_type EOS.")

  if (EOS%EOS_quadrature) then
    call int_density_dz_generic(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, EOS, dpa, &
                                intz_dpa, intx_dpa, inty_dpa, bathyT, dz_neglect, useMassWghtInterp)
  else ; select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      rho_scale = EOS%kg_m3_to_R
      if (rho_scale /= 1.0) then
        call int_density_dz_linear(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                         rho_scale*EOS%Rho_T0_S0, rho_scale*EOS%dRho_dT, rho_scale*EOS%dRho_dS, &
                         dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, dz_neglect, useMassWghtInterp)
      else
        call int_density_dz_linear(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                         EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS, &
                         dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, dz_neglect, useMassWghtInterp)
      endif
    case (EOS_WRIGHT)
      rho_scale = EOS%kg_m3_to_R
      pres_scale = EOS%RL2_T2_to_Pa
      if ((rho_scale /= 1.0) .or. (pres_scale /= 1.0)) then
        call int_density_dz_wright(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                                   dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, &
                                   dz_neglect, useMassWghtInterp, rho_scale, pres_scale)
      else
        call int_density_dz_wright(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                                   dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, &
                                   dz_neglect, useMassWghtInterp)
      endif
    case default
      call int_density_dz_generic(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, EOS, dpa, &
                                  intz_dpa, intx_dpa, inty_dpa, bathyT, dz_neglect, useMassWghtInterp)
  end select ; endif

end subroutine int_density_dz

!> Returns true if the equation of state is compressible (i.e. has pressure dependence)
logical function query_compressible(EOS)
  type(EOS_type), pointer :: EOS !< Equation of state structure

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "query_compressible called with an unassociated EOS_type EOS.")

  query_compressible = EOS%compressible
end function query_compressible

!> Initializes EOS_type by allocating and reading parameters
subroutine EOS_init(param_file, EOS, US)
  type(param_file_type), intent(in) :: param_file !< Parameter file structure
  type(EOS_type),        pointer    :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US  !< A dimensional unit scaling type
  optional :: US
  ! Local variables
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_EOS" ! This module's name.
  character(len=40)  :: tmpstr

  if (.not.associated(EOS)) call EOS_allocate(EOS)

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

end subroutine EOS_init

!> Manually initialized an EOS type (intended for unit testing of routines which need a specific EOS)
subroutine EOS_manual_init(EOS, form_of_EOS, form_of_TFreeze, EOS_quadrature, Compressible, &
                           Rho_T0_S0, drho_dT, dRho_dS, TFr_S0_P0, dTFr_dS, dTFr_dp)
  type(EOS_type),    pointer    :: EOS !< Equation of state structure
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

!> Allocates EOS_type
subroutine EOS_allocate(EOS)
  type(EOS_type), pointer :: EOS !< Equation of state structure

  if (.not.associated(EOS)) allocate(EOS)
end subroutine EOS_allocate

!> Deallocates EOS_type
subroutine EOS_end(EOS)
  type(EOS_type), pointer :: EOS !< Equation of state structure

  if (associated(EOS)) deallocate(EOS)
end subroutine EOS_end

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
  type(EOS_type),    pointer    :: EOS       !< Equation of state structure

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "MOM_EOS.F90: EOS_use_linear() called with an unassociated EOS_type EOS.")

  EOS%form_of_EOS = EOS_LINEAR
  EOS%Compressible = .false.
  EOS%Rho_T0_S0 = Rho_T0_S0
  EOS%dRho_dT = dRho_dT
  EOS%dRho_dS = dRho_dS
  EOS%EOS_quadrature = .false.
  if (present(use_quadrature)) EOS%EOS_quadrature = use_quadrature

end subroutine EOS_use_linear

!>   This subroutine calculates (by numerical quadrature) integrals of
!! pressure anomalies across layers, which are required for calculating the
!! finite-volume form pressure accelerations in a Boussinesq model.
subroutine int_density_dz_generic(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                                  EOS, dpa, intz_dpa, intx_dpa, inty_dpa, &
                                  bathyT, dz_neglect, useMassWghtInterp)
  type(hor_index_type), intent(in)  :: HI !< Horizontal index type for variables.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T  !< Potential temperature of the layer [degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S  !< Salinity of the layer [ppt]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_t !< Height at the top of the layer in depth units [Z ~> m]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_b !< Height at the bottom of the layer [Z ~> m]
  real,                 intent(in)  :: rho_ref !< A mean density [R ~> kg m-3] or [kg m-3], that is
                                          !! subtracted out to reduce the magnitude
                                          !! of each of the integrals.
  real,                 intent(in)  :: rho_0 !< A density [R ~> kg m-3] or [kg m-3], that is used
                                          !! to calculate the pressure (as p~=-z*rho_0*G_e)
                                          !! used in the equation of state.
  real,                 intent(in)  :: G_e !< The Earth's gravitational acceleration
                                          !! [L2 Z-1 T-2 ~> m s-2] or [m2 Z-1 s-2 ~> m s-2]
  type(EOS_type),       pointer     :: EOS !< Equation of state structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                      intent(inout) :: dpa !< The change in the pressure anomaly
                                          !! across the layer [R L2 T-2 ~> Pa] or [Pa]
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
  ! Local variables
  real :: T5(5), S5(5) ! Temperatures and salinities at five quadrature points [degC] and [ppt]
  real :: p5(5)      ! Pressures at five quadrature points, never rescaled from Pa [Pa]
  real :: r5(5)      ! Densities at five quadrature points [R ~> kg m-3] or [kg m-3]
  real :: rho_anom   ! The depth averaged density anomaly [R ~> kg m-3] or [kg m-3]
  real :: w_left, w_right ! Left and right weights [nondim]
  real, parameter :: C1_90 = 1.0/90.0  ! Rational constants.
  real :: GxRho      ! The gravitational acceleration times density and unit conversion factors [Pa Z-1 ~> kg m-2 s-2]
  real :: I_Rho      ! The inverse of the Boussinesq density [R-1 ~> m3 kg-1] or [m3 kg-1]
  real :: rho_scale  ! A scaling factor for densities from kg m-3 to R [R m3 kg-1 ~> 1]
  real :: rho_ref_mks ! The reference density in MKS units, never rescaled from kg m-3 [kg m-3]
  real :: dz         ! The layer thickness [Z ~> m]
  real :: hWght      ! A pressure-thickness below topography [Z ~> m]
  real :: hL, hR     ! Pressure-thicknesses of the columns to the left and right [Z ~> m]
  real :: iDenom     ! The inverse of the denominator in the weights [Z-2 ~> m-2]
  real :: hWt_LL, hWt_LR ! hWt_LA is the weighted influence of A on the left column [nondim]
  real :: hWt_RL, hWt_RR ! hWt_RA is the weighted influence of A on the right column [nondim]
  real :: wt_L, wt_R ! The linear weights of the left and right columns [nondim]
  real :: wtT_L, wtT_R ! The weights for tracers from the left and right columns [nondim]
  real :: intz(5)    ! The gravitational acceleration times the integrals of density
                     ! with height at the 5 sub-column locations [R L2 T-2 ~> Pa] or [Pa]
  logical :: do_massWeight ! Indicates whether to do mass weighting.
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, i, j, m, n

  ! These array bounds work for the indexing convention of the input arrays, but
  ! on the computational domain defined for the output arrays.
  Isq = HI%IscB ; Ieq = HI%IecB
  Jsq = HI%JscB ; Jeq = HI%JecB
  is = HI%isc ; ie = HI%iec
  js = HI%jsc ; je = HI%jec

  rho_scale = EOS%kg_m3_to_R
  GxRho = EOS%RL2_T2_to_Pa * G_e * rho_0
  rho_ref_mks = rho_ref * EOS%R_to_kg_m3
  I_Rho = 1.0 / rho_0

  do_massWeight = .false.
  if (present(useMassWghtInterp)) then ; if (useMassWghtInterp) then
    do_massWeight = .true.
    if (.not.present(bathyT)) call MOM_error(FATAL, "int_density_dz_generic: "//&
        "bathyT must be present if useMassWghtInterp is present and true.")
    if (.not.present(dz_neglect)) call MOM_error(FATAL, "int_density_dz_generic: "//&
        "dz_neglect must be present if useMassWghtInterp is present and true.")
  endif ; endif

  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    dz = z_t(i,j) - z_b(i,j)
    do n=1,5
      T5(n) = T(i,j) ; S5(n) = S(i,j)
      p5(n) = -GxRho*(z_t(i,j) - 0.25*real(n-1)*dz)
    enddo
    if (rho_scale /= 1.0) then
      call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks, scale=rho_scale)
    else
      call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks)
    endif

    ! Use Bode's rule to estimate the pressure anomaly change.
    rho_anom = C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + 12.0*r5(3))
    dpa(i,j) = G_e*dz*rho_anom
    ! Use a Bode's-rule-like fifth-order accurate estimate of the double integral of
    ! the pressure anomaly.
    if (present(intz_dpa)) intz_dpa(i,j) = 0.5*G_e*dz**2 * &
          (rho_anom - C1_90*(16.0*(r5(4)-r5(2)) + 7.0*(r5(5)-r5(1))) )
  enddo ; enddo

  if (present(intx_dpa)) then ; do j=js,je ; do I=Isq,Ieq
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., -bathyT(i,j)-z_t(i+1,j), -bathyT(i+1,j)-z_t(i,j))
    if (hWght > 0.) then
      hL = (z_t(i,j) - z_b(i,j)) + dz_neglect
      hR = (z_t(i+1,j) - z_b(i+1,j)) + dz_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    intz(1) = dpa(i,j) ; intz(5) = dpa(i+1,j)
    do m=2,4
      ! T, S, and z are interpolated in the horizontal.  The z interpolation
      ! is linear, but for T and S it may be thickness weighted.
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR
      dz = wt_L*(z_t(i,j) - z_b(i,j)) + wt_R*(z_t(i+1,j) - z_b(i+1,j))
      T5(1) = wtT_L*T(i,j) + wtT_R*T(i+1,j)
      S5(1) = wtT_L*S(i,j) + wtT_R*S(i+1,j)
      p5(1) = -GxRho*(wt_L*z_t(i,j) + wt_R*z_t(i+1,j))
      do n=2,5
        T5(n) = T5(1) ; S5(n) = S5(1) ; p5(n) = p5(n-1) + GxRho*0.25*dz
      enddo
      if (rho_scale /= 1.0) then
        call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks, scale=rho_scale)
      else
        call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks)
      endif

    ! Use Bode's rule to estimate the pressure anomaly change.
      intz(m) = G_e*dz*( C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + 12.0*r5(3)))
    enddo
    ! Use Bode's rule to integrate the bottom pressure anomaly values in x.
    intx_dpa(i,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                           12.0*intz(3))
  enddo ; enddo ; endif

  if (present(inty_dpa)) then ; do J=Jsq,Jeq ; do i=is,ie
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., -bathyT(i,j)-z_t(i,j+1), -bathyT(i,j+1)-z_t(i,j))
    if (hWght > 0.) then
      hL = (z_t(i,j) - z_b(i,j)) + dz_neglect
      hR = (z_t(i,j+1) - z_b(i,j+1)) + dz_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    intz(1) = dpa(i,j) ; intz(5) = dpa(i,j+1)
    do m=2,4
      ! T, S, and z are interpolated in the horizontal.  The z interpolation
      ! is linear, but for T and S it may be thickness weighted.
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR
      dz = wt_L*(z_t(i,j) - z_b(i,j)) + wt_R*(z_t(i,j+1) - z_b(i,j+1))
      T5(1) = wtT_L*T(i,j) + wtT_R*T(i,j+1)
      S5(1) = wtT_L*S(i,j) + wtT_R*S(i,j+1)
      p5(1) = -GxRho*(wt_L*z_t(i,j) + wt_R*z_t(i,j+1))
      do n=2,5
        T5(n) = T5(1) ; S5(n) = S5(1)
        p5(n) = p5(n-1) + GxRho*0.25*dz
      enddo
      if (rho_scale /= 1.0) then
        call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks, scale=rho_scale)
      else
        call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks)
      endif

    ! Use Bode's rule to estimate the pressure anomaly change.
      intz(m) = G_e*dz*( C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + 12.0*r5(3)))
    enddo
    ! Use Bode's rule to integrate the values.
    inty_dpa(i,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                                     12.0*intz(3))
  enddo ; enddo ; endif
end subroutine int_density_dz_generic


! ==========================================================================
!> Compute pressure gradient force integrals by quadrature for the case where
!! T and S are linear profiles.
subroutine int_density_dz_generic_plm (T_t, T_b, S_t, S_b, z_t, z_b, rho_ref, &
                                       rho_0, G_e, dz_subroundoff, bathyT, HI, EOS, dpa, &
                                       intz_dpa, intx_dpa, inty_dpa, useMassWghtInterp)
  type(hor_index_type), intent(in)  :: HI !< Ocean horizontal index structures for the arrays
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T_t !< Potential temperatue at the cell top [degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T_b !< Potential temperatue at the cell bottom [degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S_t !< Salinity at the cell top [ppt]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S_b !< Salinity at the cell bottom [ppt]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_t !< The geometric height at the top of the layer [Z ~> m]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_b !< The geometric height at the bottom of the layer [Z ~> m]
  real,                 intent(in)  :: rho_ref !< A mean density [R ~> kg m-3] or [kg m-3], that is subtracted
                                           !! out to reduce the magnitude of each of the integrals.
  real,                 intent(in)  :: rho_0 !< A density [R ~> kg m-3] or [kg m-3], that is used to calculate
                                           !! the pressure (as p~=-z*rho_0*G_e) used in the equation of state.
  real,                 intent(in)  :: G_e !< The Earth's gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real,                 intent(in)  :: dz_subroundoff !< A miniscule thickness change [Z ~> m]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: bathyT !< The depth of the bathymetry [Z ~> m]
  type(EOS_type),       pointer     :: EOS !< Equation of state structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(inout) :: dpa !< The change in the pressure anomaly across the layer [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(inout) :: intz_dpa !< The integral through the thickness of the layer of
                                           !! the pressure anomaly relative to the anomaly at the
                                           !! top of the layer [R L2 Z T-2 ~> Pa Z]
  real, dimension(HI%IsdB:HI%IedB,HI%jsd:HI%jed), &
              optional, intent(inout) :: intx_dpa !< The integral in x of the difference between the
                                           !! pressure anomaly at the top and bottom of the layer
                                           !! divided by the x grid spacing [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%JsdB:HI%JedB), &
              optional, intent(inout) :: inty_dpa !< The integral in y of the difference between the
                                           !! pressure anomaly at the top and bottom of the layer
                                           !! divided by the y grid spacing [R L2 T-2 ~> Pa]
  logical,    optional, intent(in)  :: useMassWghtInterp !< If true, uses mass weighting to
                                           !! interpolate T/S for top and bottom integrals.

! This subroutine calculates (by numerical quadrature) integrals of
! pressure anomalies across layers, which are required for calculating the
! finite-volume form pressure accelerations in a Boussinesq model.  The one
! potentially dodgy assumtion here is that rho_0 is used both in the denominator
! of the accelerations, and in the pressure used to calculated density (the
! latter being -z*rho_0*G_e).  These two uses could be separated if need be.
!
! It is assumed that the salinity and temperature profiles are linear in the
! vertical. The top and bottom values within each layer are provided and
! a linear interpolation is used to compute intermediate values.

  ! Local variables
  real :: T5((5*HI%iscB+1):(5*(HI%iecB+2)))  ! Temperatures along a line of subgrid locations [degC]
  real :: S5((5*HI%iscB+1):(5*(HI%iecB+2)))  ! Salinities along a line of subgrid locations [ppt]
  real :: p5((5*HI%iscB+1):(5*(HI%iecB+2)))  ! Pressures along a line of subgrid locations, never
                                             ! rescaled from Pa [Pa]
  real :: r5((5*HI%iscB+1):(5*(HI%iecB+2)))  ! Densities anomalies along a line of subgrid
                                             ! locations [R ~> kg m-3] or [kg m-3]
  real :: T15((15*HI%iscB+1):(15*(HI%iecB+1))) ! Temperatures at an array of subgrid locations [degC]
  real :: S15((15*HI%iscB+1):(15*(HI%iecB+1))) ! Salinities at an array of subgrid locations [ppt]
  real :: p15((15*HI%iscB+1):(15*(HI%iecB+1))) ! Pressures at an array of subgrid locations [Pa]
  real :: r15((15*HI%iscB+1):(15*(HI%iecB+1))) ! Densities at an array of subgrid locations
                                               ! [R ~> kg m-3] or [kg m-3]
  real :: wt_t(5), wt_b(5)          ! Top and bottom weights [nondim]
  real :: rho_anom                  ! A density anomaly [R ~> kg m-3] or [kg m-3]
  real :: w_left, w_right           ! Left and right weights [nondim]
  real :: intz(5)    ! The gravitational acceleration times the integrals of density
                     ! with height at the 5 sub-column locations [R L2 T-2 ~> Pa] or [Pa]
  real, parameter :: C1_90 = 1.0/90.0  ! A rational constant [nondim]
  real :: GxRho      ! The gravitational acceleration times density and unit conversion factors [Pa Z-1 ~> kg m-2 s-2]
  real :: I_Rho      ! The inverse of the Boussinesq density [R-1 ~> m3 kg-1] or [m3 kg-1]
  real :: rho_scale  ! A scaling factor for densities from kg m-3 to R [R m3 kg-1 ~> 1]
  real :: rho_ref_mks ! The reference density in MKS units, never rescaled from kg m-3 [kg m-3]
  real :: dz(HI%iscB:HI%iecB+1)   ! Layer thicknesses at tracer points [Z ~> m]
  real :: dz_x(5,HI%iscB:HI%iecB) ! Layer thicknesses along an x-line of subrid locations [Z ~> m]
  real :: dz_y(5,HI%isc:HI%iec)   ! Layer thicknesses along a y-line of subrid locations [Z ~> m]
  real :: weight_t, weight_b        ! Nondimensional weights of the top and bottom [nondim]
  real :: massWeightToggle          ! A nondimensional toggle factor (0 or 1) [nondim]
  real :: Ttl, Tbl, Ttr, Tbr        ! Temperatures at the velocity cell corners [degC]
  real :: Stl, Sbl, Str, Sbr        ! Salinities at the velocity cell corners [ppt]
  real :: hWght                     ! A topographically limited thicknes weight [Z ~> m]
  real :: hL, hR                    ! Thicknesses to the left and right [Z ~> m]
  real :: iDenom                    ! The denominator of the thickness weight expressions [Z-2 ~> m-2]
  integer :: Isq, Ieq, Jsq, Jeq, i, j, m, n
  integer :: pos

  Isq = HI%IscB ; Ieq = HI%IecB ; Jsq = HI%JscB ; Jeq = HI%JecB

  rho_scale = EOS%kg_m3_to_R
  GxRho = EOS%RL2_T2_to_Pa * G_e * rho_0
  rho_ref_mks = rho_ref * EOS%R_to_kg_m3
  I_Rho = 1.0 / rho_0
  massWeightToggle = 0.
  if (present(useMassWghtInterp)) then
    if (useMassWghtInterp) massWeightToggle = 1.
  endif

  do n = 1, 5
    wt_t(n) = 0.25 * real(5-n)
    wt_b(n) = 1.0 - wt_t(n)
  enddo

  ! =============================
  ! 1. Compute vertical integrals
  ! =============================
  do j=Jsq,Jeq+1
    do i = Isq,Ieq+1
      dz(i) = z_t(i,j) - z_b(i,j)
      do n=1,5
        p5(i*5+n) = -GxRho*(z_t(i,j) - 0.25*real(n-1)*dz(i))
        ! Salinity and temperature points are linearly interpolated
        S5(i*5+n) = wt_t(n) * S_t(i,j) + wt_b(n) * S_b(i,j)
        T5(i*5+n) = wt_t(n) * T_t(i,j) + wt_b(n) * T_b(i,j)
      enddo
    enddo
    if (rho_scale /= 1.0) then
      call calculate_density_array(T5, S5, p5, r5, 1, (ieq-isq+2)*5, EOS, rho_ref=rho_ref_mks, scale=rho_scale)
    else
      call calculate_density_array(T5, S5, p5, r5, 1, (ieq-isq+2)*5, EOS, rho_ref=rho_ref_mks)
    endif

    do i=isq,ieq+1
    ! Use Bode's rule to estimate the pressure anomaly change.
      rho_anom = C1_90*(7.0*(r5(i*5+1)+r5(i*5+5)) + 32.0*(r5(i*5+2)+r5(i*5+4)) + 12.0*r5(i*5+3))
      dpa(i,j) = G_e*dz(i)*rho_anom
      if (present(intz_dpa)) then
      ! Use a Bode's-rule-like fifth-order accurate estimate of
      ! the double integral of the pressure anomaly.
        intz_dpa(i,j) = 0.5*G_e*dz(i)**2 * &
                (rho_anom - C1_90*(16.0*(r5(i*5+4)-r5(i*5+2)) + 7.0*(r5(i*5+5)-r5(i*5+1))) )
      endif
    enddo
  enddo ! end loops on j


  ! ==================================================
  ! 2. Compute horizontal integrals in the x direction
  ! ==================================================
  if (present(intx_dpa)) then ; do j=HI%jsc,HI%jec
    do I=Isq,Ieq
      ! Corner values of T and S
      ! hWght is the distance measure by which the cell is violation of
      ! hydrostatic consistency. For large hWght we bias the interpolation
      ! of T,S along the top and bottom integrals, almost like thickness
      ! weighting.
      ! Note: To work in terrain following coordinates we could offset
      ! this distance by the layer thickness to replicate other models.
      hWght = massWeightToggle * &
              max(0., -bathyT(i,j)-z_t(i+1,j), -bathyT(i+1,j)-z_t(i,j))
      if (hWght > 0.) then
        hL = (z_t(i,j) - z_b(i,j)) + dz_subroundoff
        hR = (z_t(i+1,j) - z_b(i+1,j)) + dz_subroundoff
        hWght = hWght * ( (hL-hR)/(hL+hR) )**2
        iDenom = 1./( hWght*(hR + hL) + hL*hR )
        Ttl = ( (hWght*hR)*T_t(i+1,j) + (hWght*hL + hR*hL)*T_t(i,j) ) * iDenom
        Ttr = ( (hWght*hL)*T_t(i,j) + (hWght*hR + hR*hL)*T_t(i+1,j) ) * iDenom
        Tbl = ( (hWght*hR)*T_b(i+1,j) + (hWght*hL + hR*hL)*T_b(i,j) ) * iDenom
        Tbr = ( (hWght*hL)*T_b(i,j) + (hWght*hR + hR*hL)*T_b(i+1,j) ) * iDenom
        Stl = ( (hWght*hR)*S_t(i+1,j) + (hWght*hL + hR*hL)*S_t(i,j) ) * iDenom
        Str = ( (hWght*hL)*S_t(i,j) + (hWght*hR + hR*hL)*S_t(i+1,j) ) * iDenom
        Sbl = ( (hWght*hR)*S_b(i+1,j) + (hWght*hL + hR*hL)*S_b(i,j) ) * iDenom
        Sbr = ( (hWght*hL)*S_b(i,j) + (hWght*hR + hR*hL)*S_b(i+1,j) ) * iDenom
      else
        Ttl = T_t(i,j); Tbl = T_b(i,j); Ttr = T_t(i+1,j); Tbr = T_b(i+1,j)
        Stl = S_t(i,j); Sbl = S_b(i,j); Str = S_t(i+1,j); Sbr = S_b(i+1,j)
      endif

      do m=2,4
        w_left = 0.25*real(5-m) ; w_right = 1.0-w_left
        dz_x(m,i) = w_left*(z_t(i,j) - z_b(i,j)) + w_right*(z_t(i+1,j) - z_b(i+1,j))

        ! Salinity and temperature points are linearly interpolated in
        ! the horizontal. The subscript (1) refers to the top value in
        ! the vertical profile while subscript (5) refers to the bottom
        ! value in the vertical profile.
        pos = i*15+(m-2)*5
        T15(pos+1) = w_left*Ttl + w_right*Ttr
        T15(pos+5) = w_left*Tbl + w_right*Tbr

        S15(pos+1) = w_left*Stl + w_right*Str
        S15(pos+5) = w_left*Sbl + w_right*Sbr

        p15(pos+1) = -GxRho*(w_left*z_t(i,j) + w_right*z_t(i+1,j))

        ! Pressure
        do n=2,5
          p15(pos+n) = p15(pos+n-1) + GxRho*0.25*dz_x(m,i)
        enddo

        ! Salinity and temperature (linear interpolation in the vertical)
        do n=2,4
          weight_t = 0.25 * real(5-n)
          weight_b = 1.0 - weight_t
          S15(pos+n) = weight_t * S15(pos+1) + weight_b * S15(pos+5)
          T15(pos+n) = weight_t * T15(pos+1) + weight_b * T15(pos+5)
        enddo
      enddo
    enddo

    if (rho_scale /= 1.0) then
      call calculate_density(T15, S15, p15, r15, 1, 15*(ieq-isq+1), EOS, rho_ref=rho_ref_mks, scale=rho_scale)
    else
      call calculate_density(T15, S15, p15, r15, 1, 15*(ieq-isq+1), EOS, rho_ref=rho_ref_mks)
    endif

    do I=Isq,Ieq
      intz(1) = dpa(i,j) ; intz(5) = dpa(i+1,j)

      ! Use Bode's rule to estimate the pressure anomaly change.
      do m = 2,4
        pos = i*15+(m-2)*5
        intz(m) = G_e*dz_x(m,i)*( C1_90*(7.0*(r15(pos+1)+r15(pos+5)) + 32.0*(r15(pos+2)+r15(pos+4)) + &
                          12.0*r15(pos+3)))
      enddo
      ! Use Bode's rule to integrate the bottom pressure anomaly values in x.
      intx_dpa(i,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                             12.0*intz(3))
    enddo
  enddo ; endif

  ! ==================================================
  ! 3. Compute horizontal integrals in the y direction
  ! ==================================================
  if (present(inty_dpa)) then ; do J=Jsq,Jeq
    do i=HI%isc,HI%iec
    ! Corner values of T and S
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation
    ! of T,S along the top and bottom integrals, almost like thickness
    ! weighting.
    ! Note: To work in terrain following coordinates we could offset
    ! this distance by the layer thickness to replicate other models.
      hWght = massWeightToggle * &
              max(0., -bathyT(i,j)-z_t(i,j+1), -bathyT(i,j+1)-z_t(i,j))
      if (hWght > 0.) then
        hL = (z_t(i,j) - z_b(i,j)) + dz_subroundoff
        hR = (z_t(i,j+1) - z_b(i,j+1)) + dz_subroundoff
        hWght = hWght * ( (hL-hR)/(hL+hR) )**2
        iDenom = 1./( hWght*(hR + hL) + hL*hR )
        Ttl = ( (hWght*hR)*T_t(i,j+1) + (hWght*hL + hR*hL)*T_t(i,j) ) * iDenom
        Ttr = ( (hWght*hL)*T_t(i,j) + (hWght*hR + hR*hL)*T_t(i,j+1) ) * iDenom
        Tbl = ( (hWght*hR)*T_b(i,j+1) + (hWght*hL + hR*hL)*T_b(i,j) ) * iDenom
        Tbr = ( (hWght*hL)*T_b(i,j) + (hWght*hR + hR*hL)*T_b(i,j+1) ) * iDenom
        Stl = ( (hWght*hR)*S_t(i,j+1) + (hWght*hL + hR*hL)*S_t(i,j) ) * iDenom
        Str = ( (hWght*hL)*S_t(i,j) + (hWght*hR + hR*hL)*S_t(i,j+1) ) * iDenom
        Sbl = ( (hWght*hR)*S_b(i,j+1) + (hWght*hL + hR*hL)*S_b(i,j) ) * iDenom
        Sbr = ( (hWght*hL)*S_b(i,j) + (hWght*hR + hR*hL)*S_b(i,j+1) ) * iDenom
      else
        Ttl = T_t(i,j); Tbl = T_b(i,j); Ttr = T_t(i,j+1); Tbr = T_b(i,j+1)
        Stl = S_t(i,j); Sbl = S_b(i,j); Str = S_t(i,j+1); Sbr = S_b(i,j+1)
      endif

      do m=2,4
        w_left = 0.25*real(5-m) ; w_right = 1.0-w_left
        dz_y(m,i) = w_left*(z_t(i,j) - z_b(i,j)) + w_right*(z_t(i,j+1) - z_b(i,j+1))

        ! Salinity and temperature points are linearly interpolated in
        ! the horizontal. The subscript (1) refers to the top value in
        ! the vertical profile while subscript (5) refers to the bottom
        ! value in the vertical profile.
        pos = i*15+(m-2)*5
        T15(pos+1) = w_left*Ttl + w_right*Ttr
        T15(pos+5) = w_left*Tbl + w_right*Tbr

        S15(pos+1) = w_left*Stl + w_right*Str
        S15(pos+5) = w_left*Sbl + w_right*Sbr

        p15(pos+1) = -GxRho*(w_left*z_t(i,j) + w_right*z_t(i,j+1))

        ! Pressure
        do n=2,5 ; p15(pos+n) = p15(pos+n-1) + GxRho*0.25*dz_y(m,i) ; enddo

        ! Salinity and temperature (linear interpolation in the vertical)
        do n=2,4
          weight_t = 0.25 * real(5-n)
          weight_b = 1.0 - weight_t
          S15(pos+n) = weight_t * S15(pos+1) + weight_b * S15(pos+5)
          T15(pos+n) = weight_t * T15(pos+1) + weight_b * T15(pos+5)
        enddo
      enddo
    enddo

    if (rho_scale /= 1.0) then
      call calculate_density_array(T15(15*HI%isc+1:), S15(15*HI%isc+1:), p15(15*HI%isc+1:), &
                                   r15(15*HI%isc+1:), 1, 15*(HI%iec-HI%isc+1), EOS, &
                                   rho_ref=rho_ref_mks, scale=rho_scale)
    else
      call calculate_density_array(T15(15*HI%isc+1:), S15(15*HI%isc+1:), p15(15*HI%isc+1:), &
                                   r15(15*HI%isc+1:), 1, 15*(HI%iec-HI%isc+1), EOS, rho_ref=rho_ref_mks)
    endif
    do i=HI%isc,HI%iec
      intz(1) = dpa(i,j) ; intz(5) = dpa(i,j+1)

      ! Use Bode's rule to estimate the pressure anomaly change.
      do m = 2,4
        pos = i*15+(m-2)*5
        intz(m) = G_e*dz_y(m,i)*( C1_90*(7.0*(r15(pos+1)+r15(pos+5)) + &
                                         32.0*(r15(pos+2)+r15(pos+4)) + &
                                         12.0*r15(pos+3)))
      enddo
      ! Use Bode's rule to integrate the values.
      inty_dpa(i,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                             12.0*intz(3))
    enddo
  enddo ; endif

end subroutine int_density_dz_generic_plm
! ==========================================================================
! Above is the routine where only the S and T profiles are modified
! The real topography is still used
! ==========================================================================

!> Find the depth at which the reconstructed pressure matches P_tgt
subroutine find_depth_of_pressure_in_cell(T_t, T_b, S_t, S_b, z_t, z_b, P_t, P_tgt, &
                       rho_ref, G_e, EOS, P_b, z_out, z_tol)
  real,           intent(in)  :: T_t !< Potential temperatue at the cell top [degC]
  real,           intent(in)  :: T_b !< Potential temperatue at the cell bottom [degC]
  real,           intent(in)  :: S_t !< Salinity at the cell top [ppt]
  real,           intent(in)  :: S_b !< Salinity at the cell bottom [ppt]
  real,           intent(in)  :: z_t !< Absolute height of top of cell [Z ~> m]   (Boussinesq ????)
  real,           intent(in)  :: z_b !< Absolute height of bottom of cell [Z ~> m]
  real,           intent(in)  :: P_t !< Anomalous pressure of top of cell, relative to g*rho_ref*z_t [R L2 T-2 ~> Pa]
  real,           intent(in)  :: P_tgt !< Target pressure at height z_out, relative to g*rho_ref*z_out [R L2 T-2 ~> Pa]
  real,           intent(in)  :: rho_ref !< Reference density with which calculation are anomalous to [R ~> kg m-3]
  real,           intent(in)  :: G_e !< Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  type(EOS_type), pointer     :: EOS !< Equation of state structure
  real,           intent(out) :: P_b !< Pressure at the bottom of the cell [R L2 T-2 ~> Pa]
  real,           intent(out) :: z_out !< Absolute depth at which anomalous pressure = p_tgt [Z ~> m]
  real, optional, intent(in)  :: z_tol !< The tolerance in finding z_out [Z ~> m]

  ! Local variables
  real :: dp    ! Pressure thickness of the layer [R L2 T-2 ~> Pa]
  real :: F_guess, F_l, F_r  ! Fractional positions [nondim]
  real :: GxRho ! The product of the gravitational acceleration and reference density [R L2 Z-1 T-2 ~> Pa m-1]
  real :: Pa, Pa_left, Pa_right, Pa_tol ! Pressure anomalies, P = integral of g*(rho-rho_ref) dz [R L2 T-2 ~> Pa]
  character(len=240) :: msg

  GxRho = G_e * rho_ref

  ! Anomalous pressure difference across whole cell
  dp = frac_dp_at_pos(T_t, T_b, S_t, S_b, z_t, z_b, rho_ref, G_e, 1.0, EOS)

  P_b = P_t + dp ! Anomalous pressure at bottom of cell

  if (P_tgt <= P_t ) then
    z_out = z_t
    return
  endif

  if (P_tgt >= P_b) then
    z_out = z_b
    return
  endif

  F_l = 0.
  Pa_left = P_t - P_tgt ! Pa_left < 0
  F_r = 1.
  Pa_right = P_b - P_tgt ! Pa_right > 0
  Pa_tol = GxRho * 1.0e-5*EOS%m_to_Z
  if (present(z_tol)) Pa_tol = GxRho * z_tol

  F_guess = F_l - Pa_left / (Pa_right - Pa_left) * (F_r - F_l)
  Pa = Pa_right - Pa_left ! To get into iterative loop
  do while ( abs(Pa) > Pa_tol )

    z_out = z_t + ( z_b - z_t ) * F_guess
    Pa = frac_dp_at_pos(T_t, T_b, S_t, S_b, z_t, z_b, rho_ref, G_e, F_guess, EOS) - ( P_tgt - P_t )

    if (Pa<Pa_left) then
      write(msg,*) Pa_left,Pa,Pa_right,P_t-P_tgt,P_b-P_tgt
      call MOM_error(FATAL, 'find_depth_of_pressure_in_cell out of bounds negative: /n'//msg)
    elseif (Pa<0.) then
      Pa_left = Pa
      F_l = F_guess
    elseif (Pa>Pa_right) then
      write(msg,*) Pa_left,Pa,Pa_right,P_t-P_tgt,P_b-P_tgt
      call MOM_error(FATAL, 'find_depth_of_pressure_in_cell out of bounds positive: /n'//msg)
    elseif (Pa>0.) then
      Pa_right = Pa
      F_r = F_guess
    else ! Pa == 0
      return
    endif
    F_guess = F_l - Pa_left / (Pa_right - Pa_left) * (F_r - F_l)

  enddo

end subroutine find_depth_of_pressure_in_cell

!> Returns change in anomalous pressure change from top to non-dimensional
!! position pos between z_t and z_b
real function frac_dp_at_pos(T_t, T_b, S_t, S_b, z_t, z_b, rho_ref, G_e, pos, EOS)
  real,           intent(in)  :: T_t !< Potential temperatue at the cell top [degC]
  real,           intent(in)  :: T_b !< Potential temperatue at the cell bottom [degC]
  real,           intent(in)  :: S_t !< Salinity at the cell top [ppt]
  real,           intent(in)  :: S_b !< Salinity at the cell bottom [ppt]
  real,           intent(in)  :: z_t !< The geometric height at the top of the layer [Z ~> m]
  real,           intent(in)  :: z_b !< The geometric height at the bottom of the layer [Z ~> m]
  real,           intent(in)  :: rho_ref !< A mean density [R ~> kg m-3], that is subtracted out to
                                     !! reduce the magnitude of each of the integrals.
  real,           intent(in)  :: G_e !< The Earth's gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real,           intent(in)  :: pos !< The fractional vertical position, 0 to 1 [nondim]
  type(EOS_type), pointer     :: EOS !< Equation of state structure
  real                        :: fract_dp_at_pos !< The change in pressure from the layer top to
                                     !! fractional position pos [R L2 T-2 ~> Pa]
  ! Local variables
  real, parameter :: C1_90 = 1.0/90.0  ! A rational constant [nondim]
  real :: dz                 ! Distance from the layer top [Z ~> m]
  real :: top_weight, bottom_weight ! Fractional weights at quadrature points [nondim]
  real :: rho_ave            ! Average density [R ~> kg m-3]
  real, dimension(5) :: T5   ! Tempratures at quadrature points [degC]
  real, dimension(5) :: S5   ! Salinities at quadrature points [ppt]
  real, dimension(5) :: p5   ! Pressures at quadrature points [R L2 T-2 ~> Pa]
  real, dimension(5) :: rho5 ! Densities at quadrature points [R ~> kg m-3]
  integer :: n

  do n=1,5
    ! Evalute density at five quadrature points
    bottom_weight = 0.25*real(n-1) * pos
    top_weight = 1.0 - bottom_weight
    ! Salinity and temperature points are linearly interpolated
    S5(n) = top_weight * S_t + bottom_weight * S_b
    T5(n) = top_weight * T_t + bottom_weight * T_b
    p5(n) = ( top_weight * z_t + bottom_weight * z_b ) * ( G_e * rho_ref )
  enddo
  call calculate_density_1d(T5, S5, p5, rho5, EOS)
  rho5(:) = rho5(:) !- rho_ref ! Work with anomalies relative to rho_ref

  ! Use Bode's rule to estimate the average density
  rho_ave = C1_90*(7.0*(rho5(1)+rho5(5)) + 32.0*(rho5(2)+rho5(4)) + 12.0*rho5(3))

  dz = ( z_t - z_b ) * pos
  frac_dp_at_pos = G_e * dz * rho_ave
end function frac_dp_at_pos


! ==========================================================================
!> Compute pressure gradient force integrals for the case where T and S
!! are parabolic profiles
subroutine int_density_dz_generic_ppm(T, T_t, T_b, S, S_t, S_b, &
                                      z_t, z_b, rho_ref, rho_0, G_e, HI, &
                                      EOS, dpa, intz_dpa, intx_dpa, inty_dpa)

  type(hor_index_type), intent(in)  :: HI !< Ocean horizontal index structures for the arrays
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T   !< Potential temperature referenced to the surface [degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T_t !< Potential temperatue at the cell top [degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T_b !< Potential temperatue at the cell bottom [degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S   !< Salinity [ppt]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S_t !< Salinity at the cell top [ppt]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S_b !< Salinity at the cell bottom [ppt]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_t !< Height at the top of the layer [Z ~> m]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_b !< Height at the bottom of the layer [Z ~> m]
  real,                 intent(in)  :: rho_ref !< A mean density [R ~> kg m-3] or [kg m-3], that is
                                           !! subtracted out to reduce the magnitude of each of the integrals.
  real,                 intent(in)  :: rho_0 !< A density [R ~> kg m-3] or [kg m-3], that is used to calculate
                                           !! the pressure (as p~=-z*rho_0*G_e) used in the equation of state.
  real,                 intent(in)  :: G_e !< The Earth's gravitational acceleration [m s-2]
  type(EOS_type),       pointer     :: EOS !< Equation of state structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(inout) :: dpa !< The change in the pressure anomaly across the layer [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(inout) :: intz_dpa !< The integral through the thickness of the layer of
                                           !! the pressure anomaly relative to the anomaly at the
                                           !! top of the layer [R L2 Z T-2 ~> Pa m]
  real, dimension(HI%IsdB:HI%IedB,HI%jsd:HI%jed), &
              optional, intent(inout) :: intx_dpa !< The integral in x of the difference between the
                                           !! pressure anomaly at the top and bottom of the layer
                                           !! divided by the x grid spacing [R L2 T-2 ~> Pa]
  real, dimension(HI%isd:HI%ied,HI%JsdB:HI%JedB), &
              optional, intent(inout) :: inty_dpa !< The integral in y of the difference between the
                                           !! pressure anomaly at the top and bottom of the layer
                                           !! divided by the y grid spacing [R L2 T-2 ~> Pa]

! This subroutine calculates (by numerical quadrature) integrals of
! pressure anomalies across layers, which are required for calculating the
! finite-volume form pressure accelerations in a Boussinesq model.  The one
! potentially dodgy assumtion here is that rho_0 is used both in the denominator
! of the accelerations, and in the pressure used to calculated density (the
! latter being -z*rho_0*G_e).  These two uses could be separated if need be.
!
! It is assumed that the salinity and temperature profiles are linear in the
! vertical. The top and bottom values within each layer are provided and
! a linear interpolation is used to compute intermediate values.

!### Please note that this subroutine has not been verified to work properly!

  ! Local variables
  real :: T5(5), S5(5)
  real :: p5(5)      ! Pressures at five quadrature points, never rescaled from Pa [Pa]
  real :: r5(5)      ! Density anomalies from rho_ref at quadrature points [R ~> kg m-3] or [kg m-3]
  real :: rho_anom   ! The integrated density anomaly [R ~> kg m-3] or [kg m-3]
  real :: w_left, w_right  ! Left and right weights [nondim]
  real :: intz(5)    ! The gravitational acceleration times the integrals of density
                     ! with height at the 5 sub-column locations [R L2 T-2 ~> Pa] or [Pa]
  real, parameter :: C1_90 = 1.0/90.0  ! Rational constants.
  real :: GxRho      ! The gravitational acceleration times density and unit conversion factors [Pa Z-1 ~> kg m-2 s-2]
  real :: I_Rho      ! The inverse of the Boussinesq density [R-1 ~> m3 kg-1] or [m3 kg-1]
  real :: rho_scale  ! A scaling factor for densities from kg m-3 to R [R m3 kg-1 ~> 1]
  real :: rho_ref_mks ! The reference density in MKS units, never rescaled from kg m-3 [kg m-3]
  real :: dz
  real :: weight_t, weight_b
  real :: s0, s1, s2                   ! parabola coefficients for S [ppt]
  real :: t0, t1, t2                   ! parabola coefficients for T [degC]
  real :: xi                           ! normalized coordinate
  real :: T_top, T_mid, T_bot
  real :: S_top, S_mid, S_bot
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, i, j, m, n
  real, dimension(4) :: x, y
  real, dimension(9) :: S_node, T_node, p_node, r_node


  call MOM_error(FATAL, &
    "int_density_dz_generic_ppm: the implementation is not done yet, contact developer")

  ! These array bounds work for the indexing convention of the input arrays, but
  ! on the computational domain defined for the output arrays.
  Isq = HI%IscB ; Ieq = HI%IecB
  Jsq = HI%JscB ; Jeq = HI%JecB
  is = HI%isc ; ie = HI%iec
  js = HI%jsc ; je = HI%jec

  rho_scale = EOS%kg_m3_to_R
  GxRho = EOS%RL2_T2_to_Pa * G_e * rho_0
  rho_ref_mks = rho_ref * EOS%R_to_kg_m3
  I_Rho = 1.0 / rho_0

  ! =============================
  ! 1. Compute vertical integrals
  ! =============================
  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    dz = z_t(i,j) - z_b(i,j)

    ! Coefficients of the parabola for S
    s0 = S_t(i,j)
    s1 = 6.0 * S(i,j) - 4.0 * S_t(i,j) - 2.0 * S_b(i,j)
    s2 = 3.0 * ( S_t(i,j) + S_b(i,j) - 2.0*S(i,j) )

    ! Coefficients of the parabola for T
    t0 = T_t(i,j)
    t1 = 6.0 * T(i,j) - 4.0 * T_t(i,j) - 2.0 * T_b(i,j)
    t2 = 3.0 * ( T_t(i,j) + T_b(i,j) - 2.0*T(i,j) )

    do n=1,5
      p5(n) = -GxRho*(z_t(i,j) - 0.25*real(n-1)*dz)

      ! Parabolic reconstruction for T and S
      xi = 0.25 * ( n - 1 )
      S5(n) = s0 + s1 * xi + s2 * xi**2
      T5(n) = t0 + t1 * xi + t2 * xi**2
    enddo

    if (rho_scale /= 1.0) then
      call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks, scale=rho_scale)
    else
      call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks)
    endif

    ! Use Bode's rule to estimate the pressure anomaly change.
    rho_anom = C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + 12.0*r5(3))

    dpa(i,j) = G_e*dz*rho_anom

    ! Use a Bode's-rule-like fifth-order accurate estimate of
    ! the double integral of the pressure anomaly.
    if (present(intz_dpa)) intz_dpa(i,j) = 0.5*G_e*dz**2 * &
          (rho_anom - C1_90*(16.0*(r5(4)-r5(2)) + 7.0*(r5(5)-r5(1))) )

  enddo ; enddo ! end loops on j and i

  ! ==================================================
  ! 2. Compute horizontal integrals in the x direction
  ! ==================================================
  if (present(intx_dpa)) then ; do j=js,je ; do I=Isq,Ieq
    intz(1) = dpa(i,j) ; intz(5) = dpa(i+1,j)
    do m=2,4
      w_left = 0.25*real(5-m) ; w_right = 1.0-w_left
      dz = w_left*(z_t(i,j) - z_b(i,j)) + w_right*(z_t(i+1,j) - z_b(i+1,j))

      ! Salinity and temperature points are linearly interpolated in
      ! the horizontal. The subscript (1) refers to the top value in
      ! the vertical profile while subscript (5) refers to the bottom
      ! value in the vertical profile.
      T_top = w_left*T_t(i,j) + w_right*T_t(i+1,j)
      T_mid = w_left*T(i,j)   + w_right*T(i+1,j)
      T_bot = w_left*T_b(i,j) + w_right*T_b(i+1,j)

      S_top = w_left*S_t(i,j) + w_right*S_t(i+1,j)
      S_mid = w_left*S(i,j)   + w_right*S(i+1,j)
      S_bot = w_left*S_b(i,j) + w_right*S_b(i+1,j)

      p5(1) = -GxRho*(w_left*z_t(i,j) + w_right*z_t(i+1,j))

      ! Pressure
      do n=2,5
        p5(n) = p5(n-1) + GxRho*0.25*dz
      enddo

      ! Coefficients of the parabola for S
      s0 = S_top
      s1 = 6.0 * S_mid - 4.0 * S_top - 2.0 * S_bot
      s2 = 3.0 * ( S_top + S_bot - 2.0*S_mid )

      ! Coefficients of the parabola for T
      t0 = T_top
      t1 = 6.0 * T_mid - 4.0 * T_top - 2.0 * T_bot
      t2 = 3.0 * ( T_top + T_bot - 2.0*T_mid )

      do n=1,5
        ! Parabolic reconstruction for T and S
        xi = 0.25 * ( n - 1 )
        S5(n) = s0 + s1 * xi + s2 * xi**2
        T5(n) = t0 + t1 * xi + t2 * xi**2
      enddo

      if (rho_scale /= 1.0) then
        call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks, scale=rho_scale)
      else
        call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks)
      endif

    ! Use Bode's rule to estimate the pressure anomaly change.
      intz(m) = G_e*dz*( C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + &
                            12.0*r5(3)) )
    enddo
    intx_dpa(i,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                           12.0*intz(3))

    ! Use Gauss quadrature rule to compute integral

    ! The following coordinates define the quadrilateral on which the integral
    ! is computed
    x(1) = 1.0
    x(2) = 0.0
    x(3) = 0.0
    x(4) = 1.0
    y(1) = z_t(i+1,j)
    y(2) = z_t(i,j)
    y(3) = z_b(i,j)
    y(4) = z_b(i+1,j)

    T_node = 0.0
    p_node = 0.0

    ! Nodal values for S

    ! Parabolic reconstruction on the left
    s0 = S_t(i,j)
    s1 = 6.0 * S(i,j) - 4.0 * S_t(i,j) - 2.0 * S_b(i,j)
    s2 = 3.0 * ( S_t(i,j) + S_b(i,j) - 2.0 * S(i,j) )
    S_node(2) = s0
    S_node(6) = s0 + 0.5 * s1 + 0.25 * s2
    S_node(3) = s0 + s1 + s2

    ! Parabolic reconstruction on the left
    s0 = S_t(i+1,j)
    s1 = 6.0 * S(i+1,j) - 4.0 * S_t(i+1,j) - 2.0 * S_b(i+1,j)
    s2 = 3.0 * ( S_t(i+1,j) + S_b(i+1,j) - 2.0 * S(i+1,j) )
    S_node(1) = s0
    S_node(8) = s0 + 0.5 * s1 + 0.25 * s2
    S_node(4) = s0 + s1 + s2

    S_node(5) = 0.5 * ( S_node(2) + S_node(1) )
    S_node(9) = 0.5 * ( S_node(6) + S_node(8) )
    S_node(7) = 0.5 * ( S_node(3) + S_node(4) )

    if (rho_scale /= 1.0) then
      call calculate_density( T_node, S_node, p_node, r_node, 1, 9, EOS, rho_ref=rho_ref_mks, scale=rho_scale )
    else
      call calculate_density( T_node, S_node, p_node, r_node, 1, 9, EOS, rho_ref=rho_ref_mks)
    endif
    r_node = r_node - rho_ref

    call compute_integral_quadratic( x, y, r_node, intx_dpa(i,j) )

    intx_dpa(i,j) = intx_dpa(i,j) * G_e

  enddo ; enddo ; endif

  ! ==================================================
  ! 3. Compute horizontal integrals in the y direction
  ! ==================================================
  if (present(inty_dpa)) then
    call MOM_error(WARNING, "int_density_dz_generic_ppm still needs to be written for inty_dpa!")
    do J=Jsq,Jeq ; do i=is,ie

      inty_dpa(i,j) = 0.0

    enddo ; enddo
  endif

end subroutine int_density_dz_generic_ppm



! =============================================================================
!> Compute the integral of the quadratic function
subroutine compute_integral_quadratic( x, y, f, integral )
  real, dimension(4), intent(in)  :: x  !< The x-position of the corners
  real, dimension(4), intent(in)  :: y  !< The y-position of the corners
  real, dimension(9), intent(in)  :: f  !< The function at the quadrature points
  real,               intent(out) :: integral !< The returned integral

  ! Local variables
  integer               :: i, k
  real, dimension(9)    :: weight, xi, eta          ! integration points
  real                  :: f_k
  real                  :: dxdxi, dxdeta
  real                  :: dydxi, dydeta
  real, dimension(4)    :: phiiso, dphiisodxi, dphiisodeta
  real, dimension(9)    :: phi, dphidxi, dphideta
  real                  :: jacobian_k
  real                  :: t

  ! Quadrature rule (4 points)
  !weight(:) = 1.0
  !xi(1) = - sqrt(3.0) / 3.0
  !xi(2) = sqrt(3.0) / 3.0
  !xi(3) = sqrt(3.0) / 3.0
  !xi(4) = - sqrt(3.0) / 3.0
  !eta(1) = - sqrt(3.0) / 3.0
  !eta(2) = - sqrt(3.0) / 3.0
  !eta(3) = sqrt(3.0) / 3.0
  !eta(4) = sqrt(3.0) / 3.0

  ! Quadrature rule (9 points)
  t = sqrt(3.0/5.0)
  weight(1) = 25.0/81.0 ; xi(1) = -t ; eta(1) = t
  weight(2) = 40.0/81.0 ; xi(2) = .0 ; eta(2) = t
  weight(3) = 25.0/81.0 ; xi(3) =  t ; eta(3) = t
  weight(4) = 40.0/81.0 ; xi(4) = -t ; eta(4) = .0
  weight(5) = 64.0/81.0 ; xi(5) = .0 ; eta(5) = .0
  weight(6) = 40.0/81.0 ; xi(6) =  t ; eta(6) = .0
  weight(7) = 25.0/81.0 ; xi(7) = -t ; eta(7) = -t
  weight(8) = 40.0/81.0 ; xi(8) = .0 ; eta(8) = -t
  weight(9) = 25.0/81.0 ; xi(9) =  t ; eta(9) = -t

  integral = 0.0

  ! Integration loop
  do k = 1,9

    ! Evaluate shape functions and gradients for isomorphism
    call evaluate_shape_bilinear( xi(k), eta(k), phiiso, &
                                  dphiisodxi, dphiisodeta )

    ! Determine gradient of global coordinate at integration point
    dxdxi  = 0.0
    dxdeta = 0.0
    dydxi  = 0.0
    dydeta = 0.0

    do i = 1,4
      dxdxi  = dxdxi  + x(i) * dphiisodxi(i)
      dxdeta = dxdeta + x(i) * dphiisodeta(i)
      dydxi  = dydxi  + y(i) * dphiisodxi(i)
      dydeta = dydeta + y(i) * dphiisodeta(i)
    enddo

    ! Evaluate Jacobian at integration point
    jacobian_k = dxdxi*dydeta - dydxi*dxdeta

    ! Evaluate shape functions for interpolation
    call evaluate_shape_quadratic( xi(k), eta(k), phi, dphidxi, dphideta )

    ! Evaluate function at integration point
    f_k = 0.0
    do i = 1,9
      f_k = f_k + f(i) * phi(i)
    enddo

    integral = integral + weight(k) * f_k * jacobian_k

  enddo ! end integration loop

end subroutine compute_integral_quadratic


! =============================================================================
!> Evaluation of the four bilinear shape fn and their gradients at (xi,eta)
subroutine evaluate_shape_bilinear( xi, eta, phi, dphidxi, dphideta )
  real,               intent(in)  :: xi  !< The x position to evaluate
  real,               intent(in)  :: eta !< The z position to evaluate
  real, dimension(4), intent(inout) :: phi !< The weights of the four corners at this point
  real, dimension(4), intent(inout) :: dphidxi  !< The x-gradient of the weights of the four
                                         !! corners at this point
  real, dimension(4), intent(inout) :: dphideta !< The z-gradient of the weights of the four
                                         !! corners at this point

  ! The shape functions within the parent element are defined as shown here:
  !
  !    (-1,1) 2 o------------o 1 (1,1)
  !             |            |
  !             |            |
  !             |            |
  !             |            |
  !   (-1,-1) 3 o------------o 4 (1,-1)
  !

  phi(1) = 0.25 * ( 1 + xi ) * ( 1 + eta )
  phi(2) = 0.25 * ( 1 - xi ) * ( 1 + eta )
  phi(3) = 0.25 * ( 1 - xi ) * ( 1 - eta )
  phi(4) = 0.25 * ( 1 + xi ) * ( 1 - eta )

  dphidxi(1) = 0.25 * ( 1 + eta )
  dphidxi(2) = - 0.25 * ( 1 + eta )
  dphidxi(3) = - 0.25 * ( 1 - eta )
  dphidxi(4) = 0.25 * ( 1 - eta )

  dphideta(1) = 0.25 * ( 1 + xi )
  dphideta(2) = 0.25 * ( 1 - xi )
  dphideta(3) = - 0.25 * ( 1 - xi )
  dphideta(4) = - 0.25 * ( 1 + xi )

end subroutine evaluate_shape_bilinear


! =============================================================================
!> Evaluation of the nine quadratic shape fn weights and their gradients at (xi,eta)
subroutine evaluate_shape_quadratic ( xi, eta, phi, dphidxi, dphideta )

  ! Arguments
  real,               intent(in)  :: xi  !< The x position to evaluate
  real,               intent(in)  :: eta !< The z position to evaluate
  real, dimension(9), intent(inout) :: phi !< The weights of the 9 bilinear quadrature points
                                         !! at this point
  real, dimension(9), intent(inout) :: dphidxi  !< The x-gradient of the weights of the 9 bilinear
                                         !! quadrature points corners at this point
  real, dimension(9), intent(inout) :: dphideta !< The z-gradient of the weights of the 9 bilinear
                                         !! quadrature points corners at this point

  ! The quadratic shape functions within the parent element are defined as shown here:
  !
  !                 5 (0,1)
  !    (-1,1) 2 o------o------o 1 (1,1)
  !             |             |
  !             |   9 (0,0)   |
  !    (-1,0) 6 o      o      o 8 (1,0)
  !             |             |
  !             |             |
  !   (-1,-1) 3 o------o------o 4 (1,-1)
  !                 7 (0,-1)
  !

  phi(:)   = 0.0
  dphidxi(:)  = 0.0
  dphideta(:) = 0.0

  phi(1) = 0.25 * xi * ( 1 + xi ) * eta * ( 1 + eta )
  phi(2) = - 0.25 * xi * ( 1 - xi ) * eta * ( 1 + eta )
  phi(3) = 0.25 * xi * ( 1 - xi ) * eta * ( 1 - eta )
  phi(4) = - 0.25 * xi * ( 1 + xi ) * eta * ( 1 - eta )
  phi(5) = 0.5 * ( 1 + xi ) * ( 1 - xi ) * eta * ( 1 + eta )
  phi(6) = - 0.5 * xi * ( 1 - xi ) * ( 1 - eta ) * ( 1 + eta )
  phi(7) = - 0.5 * ( 1 - xi ) * ( 1 + xi ) * eta * ( 1 - eta )
  phi(8) = 0.5 * xi * ( 1 + xi ) * ( 1 - eta ) * ( 1 + eta )
  phi(9) = ( 1 - xi ) * ( 1 + xi ) * ( 1 - eta ) * ( 1 + eta )

  !dphidxi(1) = 0.25 * ( 1 + 2*xi ) * eta * ( 1 + eta )
  !dphidxi(2) = - 0.25 * ( 1 - 2*xi ) * eta * ( 1 + eta )
  !dphidxi(3) = 0.25 * ( 1 - 2*xi ) * eta * ( 1 - eta )
  !dphidxi(4) = - 0.25 * ( 1 + 2*xi ) * eta * ( 1 - eta )
  !dphidxi(5) = - xi * eta * ( 1 + eta )
  !dphidxi(6) = - 0.5 * ( 1 - 2*xi ) * ( 1 - eta ) * ( 1 + eta )
  !dphidxi(7) = xi * eta * ( 1 - eta )
  !dphidxi(8) = 0.5 * ( 1 + 2*xi ) * ( 1 - eta ) * ( 1 + eta )
  !dphidxi(9) = - 2 * xi * ( 1 - eta ) * ( 1 + eta )

  !dphideta(1) = 0.25 * xi * ( 1 + xi ) * ( 1 + 2*eta )
  !dphideta(2) = - 0.25 * xi * ( 1 - xi ) * ( 1 + 2*eta )
  !dphideta(3) = 0.25 * xi * ( 1 - xi ) * ( 1 - 2*eta )
  !dphideta(4) = - 0.25 * xi * ( 1 + xi ) * ( 1 - 2*eta )
  !dphideta(5) = 0.5 * ( 1 + xi ) * ( 1 - xi ) * ( 1 + 2*eta )
  !dphideta(6) = xi * ( 1 - xi ) * eta
  !dphideta(7) = - 0.5 * ( 1 - xi ) * ( 1 + xi ) * ( 1 - 2*eta )
  !dphideta(8) = - xi * ( 1 + xi ) * eta
  !dphideta(9) = - 2 * ( 1 - xi ) * ( 1 + xi ) * eta

end subroutine evaluate_shape_quadratic
! ==============================================================================

!>   This subroutine calculates integrals of specific volume anomalies in
!! pressure across layers, which are required for calculating the finite-volume
!! form pressure accelerations in a non-Boussinesq model.  There are essentially
!! no free assumptions, apart from the use of Bode's rule quadrature to do the integrals.
subroutine int_spec_vol_dp_generic(T, S, p_t, p_b, alpha_ref, HI, EOS, dza, &
                                   intp_dza, intx_dza, inty_dza, halo_size, &
                                   bathyP, dP_neglect, useMassWghtInterp)
  type(hor_index_type), intent(in)  :: HI !< A horizontal index type structure.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T  !< Potential temperature of the layer [degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S  !< Salinity of the layer [ppt]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: p_t !< Pressure atop the layer [R L2 T-2 ~> Pa] or [Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: p_b !< Pressure below the layer [R L2 T-2 ~> Pa] or [Pa]
  real,                 intent(in)  :: alpha_ref !< A mean specific volume that is subtracted out
                            !! to reduce the magnitude of each of the integrals [R-1 ~> m3 kg-1]
                            !! The calculation is mathematically identical with different values of
                            !! alpha_ref, but alpha_ref alters the effects of roundoff, and
                            !! answers do change.
  type(EOS_type),       pointer     :: EOS !< Equation of state structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(inout) :: dza !< The change in the geopotential anomaly
                            !! across the layer [L2 T-2 ~> m2 s-2] or [m2 s-2]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(inout) :: intp_dza !< The integral in pressure through the layer of
                            !! the geopotential anomaly relative to the anomaly at the bottom of the
                            !! layer [R L4 T-4 ~> Pa m2 s-2] or [Pa m2 s-2]
  real, dimension(HI%IsdB:HI%IedB,HI%jsd:HI%jed), &
              optional, intent(inout) :: intx_dza  !< The integral in x of the difference between
                            !! the geopotential anomaly at the top and bottom of the layer divided
                            !! by the x grid spacing [L2 T-2 ~> m2 s-2] or [m2 s-2]
  real, dimension(HI%isd:HI%ied,HI%JsdB:HI%JedB), &
              optional, intent(inout) :: inty_dza  !< The integral in y of the difference between
                            !! the geopotential anomaly at the top and bottom of the layer divided
                            !! by the y grid spacing [L2 T-2 ~> m2 s-2] or [m2 s-2]
  integer,    optional, intent(in)  :: halo_size !< The width of halo points on which to calculate dza.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: bathyP !< The pressure at the bathymetry [R L2 T-2 ~> Pa] or [Pa]
  real,       optional, intent(in)  :: dP_neglect !< A miniscule pressure change with
                                             !! the same units as p_t [R L2 T-2 ~> Pa] or [Pa]
  logical,    optional, intent(in)  :: useMassWghtInterp !< If true, uses mass weighting
                            !! to interpolate T/S for top and bottom integrals.

!   This subroutine calculates analytical and nearly-analytical integrals in
! pressure across layers of geopotential anomalies, which are required for
! calculating the finite-volume form pressure accelerations in a non-Boussinesq
! model.  There are essentially no free assumptions, apart from the use of
! Bode's rule to do the horizontal integrals, and from a truncation in the
! series for log(1-eps/1+eps) that assumes that |eps| < 0.34.

  ! Local variables
  real :: T5(5)      ! Temperatures at five quadrature points [degC]
  real :: S5(5)      ! Salinities at five quadrature points [ppt]
  real :: p5(5)      ! Pressures at five quadrature points, scaled back to Pa if necessary [Pa]
  real :: a5(5)      ! Specific volumes at five quadrature points [R-1 ~> m3 kg-1] or [m3 kg-1]
  real :: alpha_anom ! The depth averaged specific density anomaly [R-1 ~> m3 kg-1]
  real :: dp         ! The pressure change through a layer [R L2 T-2 ~> Pa]
  real :: hWght      ! A pressure-thickness below topography [R L2 T-2 ~> Pa]
  real :: hL, hR     ! Pressure-thicknesses of the columns to the left and right [R L2 T-2 ~> Pa]
  real :: alpha_ref_mks ! The reference specific volume in MKS units, never rescaled from m3 kg-1 [m3 kg-1]
  real :: iDenom     ! The inverse of the denominator in the weights [T4 R-2 L-4 ~> Pa-2]
  real :: hWt_LL, hWt_LR ! hWt_LA is the weighted influence of A on the left column [nondim]
  real :: hWt_RL, hWt_RR ! hWt_RA is the weighted influence of A on the right column [nondim]
  real :: wt_L, wt_R ! The linear weights of the left and right columns [nondim]
  real :: wtT_L, wtT_R ! The weights for tracers from the left and right columns [nondim]
  real :: intp(5)    ! The integrals of specific volume with pressure at the
                     ! 5 sub-column locations [L2 T-2 ~> m2 s-2]
  real :: RL2_T2_to_Pa  ! A unit conversion factor from the rescaled units of pressure to Pa [Pa T2 R-1 L-2 ~> 1]
  real :: SV_scale   ! A multiplicative factor by which to scale specific
                     ! volume from m3 kg-1 to the desired units [kg m-3 R-1 ~> 1]
  logical :: do_massWeight ! Indicates whether to do mass weighting.
  real, parameter :: C1_90 = 1.0/90.0  ! A rational constant.
  integer :: Isq, Ieq, Jsq, Jeq, ish, ieh, jsh, jeh, i, j, m, n, halo

  Isq = HI%IscB ; Ieq = HI%IecB ; Jsq = HI%JscB ; Jeq = HI%JecB
  halo = 0 ; if (present(halo_size)) halo = MAX(halo_size,0)
  ish = HI%isc-halo ; ieh = HI%iec+halo ; jsh = HI%jsc-halo ; jeh = HI%jec+halo
  if (present(intx_dza)) then ; ish = MIN(Isq,ish) ; ieh = MAX(Ieq+1,ieh); endif
  if (present(inty_dza)) then ; jsh = MIN(Jsq,jsh) ; jeh = MAX(Jeq+1,jeh); endif

  SV_scale = EOS%R_to_kg_m3
  RL2_T2_to_Pa = EOS%RL2_T2_to_Pa
  alpha_ref_mks = alpha_ref * EOS%kg_m3_to_R

  do_massWeight = .false.
  if (present(useMassWghtInterp)) then ; if (useMassWghtInterp) then
    do_massWeight = .true.
    if (.not.present(bathyP)) call MOM_error(FATAL, "int_spec_vol_dp_generic: "//&
        "bathyP must be present if useMassWghtInterp is present and true.")
    if (.not.present(dP_neglect)) call MOM_error(FATAL, "int_spec_vol_dp_generic: "//&
        "dP_neglect must be present if useMassWghtInterp is present and true.")
  endif ; endif

  do j=jsh,jeh ; do i=ish,ieh
    dp = p_b(i,j) - p_t(i,j)
    do n=1,5
      T5(n) = T(i,j) ; S5(n) = S(i,j)
      p5(n) = RL2_T2_to_Pa * (p_b(i,j) - 0.25*real(n-1)*dp)
    enddo

    if (SV_scale /= 1.0) then
      call calculate_spec_vol(T5, S5, p5, a5, 1, 5, EOS, alpha_ref_mks, scale=SV_scale)
    else
      call calculate_spec_vol(T5, S5, p5, a5, 1, 5, EOS, alpha_ref_mks)
    endif

    ! Use Bode's rule to estimate the interface height anomaly change.
    alpha_anom = C1_90*(7.0*(a5(1)+a5(5)) + 32.0*(a5(2)+a5(4)) + 12.0*a5(3))
    dza(i,j) = dp*alpha_anom
    ! Use a Bode's-rule-like fifth-order accurate estimate of the double integral of
    ! the interface height anomaly.
    if (present(intp_dza)) intp_dza(i,j) = 0.5*dp**2 * &
          (alpha_anom - C1_90*(16.0*(a5(4)-a5(2)) + 7.0*(a5(5)-a5(1))) )
  enddo ; enddo

  if (present(intx_dza)) then ; do j=HI%jsc,HI%jec ; do I=Isq,Ieq
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., bathyP(i,j)-p_t(i+1,j), bathyP(i+1,j)-p_t(i,j))
    if (hWght > 0.) then
      hL = (p_b(i,j) - p_t(i,j)) + dP_neglect
      hR = (p_b(i+1,j) - p_t(i+1,j)) + dP_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    intp(1) = dza(i,j) ; intp(5) = dza(i+1,j)
    do m=2,4
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR

      ! T, S, and p are interpolated in the horizontal.  The p interpolation
      ! is linear, but for T and S it may be thickness wekghted.
      p5(1) = RL2_T2_to_Pa * (wt_L*p_b(i,j) + wt_R*p_b(i+1,j))
      dp = wt_L*(p_b(i,j) - p_t(i,j)) + wt_R*(p_b(i+1,j) - p_t(i+1,j))
      T5(1) = wtT_L*T(i,j) + wtT_R*T(i+1,j)
      S5(1) = wtT_L*S(i,j) + wtT_R*S(i+1,j)

      do n=2,5
        T5(n) = T5(1) ; S5(n) = S5(1) ; p5(n) = p5(n-1) - RL2_T2_to_Pa * 0.25*dp
      enddo
      if (SV_scale /= 1.0) then
        call calculate_spec_vol(T5, S5, p5, a5, 1, 5, EOS, alpha_ref_mks, scale=SV_scale)
      else
        call calculate_spec_vol(T5, S5, p5, a5, 1, 5, EOS, alpha_ref_mks)
      endif

    ! Use Bode's rule to estimate the interface height anomaly change.
      intp(m) = dp*( C1_90*(7.0*(a5(1)+a5(5)) + 32.0*(a5(2)+a5(4)) + &
                                12.0*a5(3)))
    enddo
    ! Use Bode's rule to integrate the interface height anomaly values in x.
    intx_dza(i,j) = C1_90*(7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4)) + &
                           12.0*intp(3))
  enddo ; enddo ; endif

  if (present(inty_dza)) then ; do J=Jsq,Jeq ; do i=HI%isc,HI%iec
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., bathyP(i,j)-p_t(i,j+1), bathyP(i,j+1)-p_t(i,j))
    if (hWght > 0.) then
      hL = (p_b(i,j) - p_t(i,j)) + dP_neglect
      hR = (p_b(i,j+1) - p_t(i,j+1)) + dP_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    intp(1) = dza(i,j) ; intp(5) = dza(i,j+1)
    do m=2,4
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR

      ! T, S, and p are interpolated in the horizontal.  The p interpolation
      ! is linear, but for T and S it may be thickness wekghted.
      p5(1) = RL2_T2_to_Pa * (wt_L*p_b(i,j) + wt_R*p_b(i,j+1))
      dp = wt_L*(p_b(i,j) - p_t(i,j)) + wt_R*(p_b(i,j+1) - p_t(i,j+1))
      T5(1) = wtT_L*T(i,j) + wtT_R*T(i,j+1)
      S5(1) = wtT_L*S(i,j) + wtT_R*S(i,j+1)
      do n=2,5
        T5(n) = T5(1) ; S5(n) = S5(1) ; p5(n) = RL2_T2_to_Pa * (p5(n-1) - 0.25*dp)
      enddo
      if (SV_scale /= 1.0) then
        call calculate_spec_vol(T5, S5, p5, a5, 1, 5, EOS, alpha_ref_mks, scale=SV_scale)
      else
        call calculate_spec_vol(T5, S5, p5, a5, 1, 5, EOS, alpha_ref_mks)
      endif

    ! Use Bode's rule to estimate the interface height anomaly change.
      intp(m) = dp*( C1_90*(7.0*(a5(1)+a5(5)) + 32.0*(a5(2)+a5(4)) + &
                                12.0*a5(3)))
    enddo
    ! Use Bode's rule to integrate the interface height anomaly values in y.
    inty_dza(i,j) = C1_90*(7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4)) + &
                           12.0*intp(3))
  enddo ; enddo ; endif

end subroutine int_spec_vol_dp_generic

!>   This subroutine calculates integrals of specific volume anomalies in
!! pressure across layers, which are required for calculating the finite-volume
!! form pressure accelerations in a non-Boussinesq model.  There are essentially
!! no free assumptions, apart from the use of Bode's rule quadrature to do the integrals.
subroutine int_spec_vol_dp_generic_plm(T_t, T_b, S_t, S_b, p_t, p_b, alpha_ref, &
                             dP_neglect, bathyP, HI, EOS, dza, &
                             intp_dza, intx_dza, inty_dza, useMassWghtInterp)
  type(hor_index_type), intent(in)  :: HI !< A horizontal index type structure.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T_t  !< Potential temperature at the top of the layer [degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T_b  !< Potential temperature at the bottom of the layer [degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S_t  !< Salinity at the top the layer [ppt]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S_b  !< Salinity at the bottom the layer [ppt]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: p_t !< Pressure atop the layer [R L2 T-2 ~> Pa] or [Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: p_b !< Pressure below the layer [R L2 T-2 ~> Pa] or [Pa]
  real,                 intent(in)  :: alpha_ref !< A mean specific volume that is subtracted out
                            !! to reduce the magnitude of each of the integrals [R-1 ~> m3 kg-1]
                            !! The calculation is mathematically identical with different values of
                            !! alpha_ref, but alpha_ref alters the effects of roundoff, and
                            !! answers do change.
  real,                 intent(in)  :: dP_neglect !<!< A miniscule pressure change with
                                             !! the same units as p_t [R L2 T-2 ~> Pa] or [Pa]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: bathyP !< The pressure at the bathymetry [R L2 T-2 ~> Pa] or [Pa]
  type(EOS_type),       pointer     :: EOS !< Equation of state structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(inout) :: dza !< The change in the geopotential anomaly
                            !! across the layer [L2 T-2 ~> m2 s-2]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(inout) :: intp_dza !< The integral in pressure through the layer of
                            !! the geopotential anomaly relative to the anomaly at the bottom of the
                            !! layer [R L4 T-4 ~> Pa m2 s-2] or [Pa m2 s-2]
  real, dimension(HI%IsdB:HI%IedB,HI%jsd:HI%jed), &
              optional, intent(inout) :: intx_dza  !< The integral in x of the difference between
                            !! the geopotential anomaly at the top and bottom of the layer divided
                            !! by the x grid spacing [L2 T-2 ~> m2 s-2] or [m2 s-2]
  real, dimension(HI%isd:HI%ied,HI%JsdB:HI%JedB), &
              optional, intent(inout) :: inty_dza  !< The integral in y of the difference between
                            !! the geopotential anomaly at the top and bottom of the layer divided
                            !! by the y grid spacing [L2 T-2 ~> m2 s-2] or [m2 s-2]
  logical,    optional, intent(in)  :: useMassWghtInterp !< If true, uses mass weighting
                            !! to interpolate T/S for top and bottom integrals.

!   This subroutine calculates analytical and nearly-analytical integrals in
! pressure across layers of geopotential anomalies, which are required for
! calculating the finite-volume form pressure accelerations in a non-Boussinesq
! model.  There are essentially no free assumptions, apart from the use of
! Bode's rule to do the horizontal integrals, and from a truncation in the
! series for log(1-eps/1+eps) that assumes that |eps| < 0.34.

  real :: T5(5)      ! Temperatures at five quadrature points [degC]
  real :: S5(5)      ! Salinities at five quadrature points [ppt]
  real :: p5(5)      ! Pressures at five quadrature points, scaled back to Pa as necessary [Pa]
  real :: a5(5)      ! Specific volumes at five quadrature points [R-1 ~> m3 kg-1] or [m3 kg-1]
  real :: T15(15)    ! Temperatures at fifteen interior quadrature points [degC]
  real :: S15(15)    ! Salinities at fifteen interior quadrature points [ppt]
  real :: p15(15)    ! Pressures at fifteen quadrature points, scaled back to Pa as necessary [Pa]
  real :: a15(15)    ! Specific volumes at fifteen quadrature points [R-1 ~> m3 kg-1] or [m3 kg-1]
  real :: wt_t(5), wt_b(5) ! Weights of top and bottom values at quadrature points [nondim]
  real :: T_top, T_bot, S_top, S_bot, P_top, P_bot

  real :: alpha_anom ! The depth averaged specific density anomaly [m3 kg-1]
  real :: dp         ! The pressure change through a layer [R L2 T-2 ~> Pa]
  real :: dp_90(2:4) ! The pressure change through a layer divided by 90 [R L2 T-2 ~> Pa]
  real :: hWght      ! A pressure-thickness below topography [R L2 T-2 ~> Pa]
  real :: hL, hR     ! Pressure-thicknesses of the columns to the left and right [R L2 T-2 ~> Pa]
  real :: alpha_ref_mks ! The reference specific volume in MKS units, never rescaled from m3 kg-1 [m3 kg-1]
  real :: iDenom     ! The inverse of the denominator in the weights [T4 R-2 L-4 ~> Pa-2]
  real :: hWt_LL, hWt_LR ! hWt_LA is the weighted influence of A on the left column [nondim]
  real :: hWt_RL, hWt_RR ! hWt_RA is the weighted influence of A on the right column [nondim]
  real :: wt_L, wt_R ! The linear weights of the left and right columns [nondim]
  real :: wtT_L, wtT_R ! The weights for tracers from the left and right columns [nondim]
  real :: intp(5)    ! The integrals of specific volume with pressure at the
                     ! 5 sub-column locations [L2 T-2 ~> m2 s-2]
  real :: RL2_T2_to_Pa  ! A unit conversion factor from the rescaled units of pressure to Pa [Pa T2 R-1 L-2 ~> 1]
  real :: SV_scale   ! A multiplicative factor by which to scale specific
                     ! volume from m3 kg-1 to the desired units [kg m-3 R-1 ~> 1]
  real, parameter :: C1_90 = 1.0/90.0  ! A rational constant.
  logical :: do_massWeight ! Indicates whether to do mass weighting.
  integer :: Isq, Ieq, Jsq, Jeq, i, j, m, n, pos

  Isq = HI%IscB ; Ieq = HI%IecB ; Jsq = HI%JscB ; Jeq = HI%JecB

  do_massWeight = .false.
  if (present(useMassWghtInterp)) do_massWeight = useMassWghtInterp

  SV_scale = EOS%R_to_kg_m3
  RL2_T2_to_Pa = EOS%RL2_T2_to_Pa
  alpha_ref_mks = alpha_ref * EOS%kg_m3_to_R

  do n = 1, 5 ! Note that these are reversed from int_density_dz.
    wt_t(n) = 0.25 * real(n-1)
    wt_b(n) = 1.0 - wt_t(n)
  enddo

  ! =============================
  ! 1. Compute vertical integrals
  ! =============================
  do j=Jsq,Jeq+1; do i=Isq,Ieq+1
    dp = p_b(i,j) - p_t(i,j)
    do n=1,5 ! T, S and p are linearly interpolated in the vertical.
      p5(n) = RL2_T2_to_Pa * (wt_t(n) * p_t(i,j) + wt_b(n) * p_b(i,j))
      S5(n) = wt_t(n) * S_t(i,j) + wt_b(n) * S_b(i,j)
      T5(n) = wt_t(n) * T_t(i,j) + wt_b(n) * T_b(i,j)
    enddo
    if (SV_scale /= 1.0) then
      call calculate_spec_vol(T5, S5, p5, a5, 1, 5, EOS, alpha_ref_mks, scale=SV_scale)
    else
      call calculate_spec_vol(T5, S5, p5, a5, 1, 5, EOS, alpha_ref_mks)
    endif

    ! Use Bode's rule to estimate the interface height anomaly change.
    alpha_anom = C1_90*((7.0*(a5(1)+a5(5)) + 32.0*(a5(2)+a5(4))) + 12.0*a5(3))
    dza(i,j) = dp*alpha_anom
    ! Use a Bode's-rule-like fifth-order accurate estimate of the double integral of
    ! the interface height anomaly.
    if (present(intp_dza)) intp_dza(i,j) = 0.5*dp**2 * &
          (alpha_anom - C1_90*(16.0*(a5(4)-a5(2)) + 7.0*(a5(5)-a5(1))) )
  enddo ; enddo

  ! ==================================================
  ! 2. Compute horizontal integrals in the x direction
  ! ==================================================
  if (present(intx_dza)) then ; do j=HI%jsc,HI%jec ; do I=Isq,Ieq
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation
    ! of T,S along the top and bottom integrals, almost like thickness
    ! weighting. Note: To work in terrain following coordinates we could
    ! offset this distance by the layer thickness to replicate other models.
    hWght = 0.0
    if (do_massWeight) &
      hWght =  max(0., bathyP(i,j)-p_t(i+1,j), bathyP(i+1,j)-p_t(i,j))
    if (hWght > 0.) then
      hL = (p_b(i,j) - p_t(i,j)) + dP_neglect
      hR = (p_b(i+1,j) - p_t(i+1,j)) + dP_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    do m=2,4
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR

      ! T, S, and p are interpolated in the horizontal.  The p interpolation
      ! is linear, but for T and S it may be thickness wekghted.
      P_top = wt_L*p_t(i,j) + wt_R*p_t(i+1,j)
      P_bot = wt_L*p_b(i,j) + wt_R*p_b(i+1,j)
      T_top = wtT_L*T_t(i,j) + wtT_R*T_t(i+1,j)
      T_bot = wtT_L*T_b(i,j) + wtT_R*T_b(i+1,j)
      S_top = wtT_L*S_t(i,j) + wtT_R*S_t(i+1,j)
      S_bot = wtT_L*S_b(i,j) + wtT_R*S_b(i+1,j)
      dp_90(m) = C1_90*(P_bot - P_top)

      ! Salinity, temperature and pressure with linear interpolation in the vertical.
      pos = (m-2)*5
      do n=1,5
        p15(pos+n) = RL2_T2_to_Pa * (wt_t(n) * P_top + wt_b(n) * P_bot)
        S15(pos+n) = wt_t(n) * S_top + wt_b(n) * S_bot
        T15(pos+n) = wt_t(n) * T_top + wt_b(n) * T_bot
      enddo
    enddo

    if (SV_scale /= 1.0) then
      call calculate_spec_vol(T15, S15, p15, a15, 1, 15, EOS, alpha_ref_mks, scale=SV_scale)
    else
      call calculate_spec_vol(T15, S15, p15, a15, 1, 15, EOS, alpha_ref_mks)
    endif

    intp(1) = dza(i,j) ; intp(5) = dza(i+1,j)
    do m=2,4
      ! Use Bode's rule to estimate the interface height anomaly change.
      ! The integrals at the ends of the segment are already known.
      pos = (m-2)*5
      intp(m) = dp_90(m)*((7.0*(a15(pos+1)+a15(pos+5)) + &
                          32.0*(a15(pos+2)+a15(pos+4))) + 12.0*a15(pos+3))
    enddo
    ! Use Bode's rule to integrate the interface height anomaly values in x.
    intx_dza(I,j) = C1_90*((7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4))) + &
                           12.0*intp(3))
  enddo ; enddo ; endif

  ! ==================================================
  ! 3. Compute horizontal integrals in the y direction
  ! ==================================================
  if (present(inty_dza)) then ; do J=Jsq,Jeq ; do i=HI%isc,HI%iec
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation
    ! of T,S along the top and bottom integrals, like thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., bathyP(i,j)-p_t(i,j+1), bathyP(i,j+1)-p_t(i,j))
    if (hWght > 0.) then
      hL = (p_b(i,j) - p_t(i,j)) + dP_neglect
      hR = (p_b(i,j+1) - p_t(i,j+1)) + dP_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    do m=2,4
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR

      ! T, S, and p are interpolated in the horizontal.  The p interpolation
      ! is linear, but for T and S it may be thickness wekghted.
      P_top = wt_L*p_t(i,j) + wt_R*p_t(i,j+1)
      P_bot = wt_L*p_b(i,j) + wt_R*p_b(i,j+1)
      T_top = wtT_L*T_t(i,j) + wtT_R*T_t(i,j+1)
      T_bot = wtT_L*T_b(i,j) + wtT_R*T_b(i,j+1)
      S_top = wtT_L*S_t(i,j) + wtT_R*S_t(i,j+1)
      S_bot = wtT_L*S_b(i,j) + wtT_R*S_b(i,j+1)
      dp_90(m) = C1_90*(P_bot - P_top)

      ! Salinity, temperature and pressure with linear interpolation in the vertical.
      pos = (m-2)*5
      do n=1,5
        p15(pos+n) = RL2_T2_to_Pa * (wt_t(n) * P_top + wt_b(n) * P_bot)
        S15(pos+n) = wt_t(n) * S_top + wt_b(n) * S_bot
        T15(pos+n) = wt_t(n) * T_top + wt_b(n) * T_bot
      enddo
    enddo

    if (SV_scale /= 1.0) then
      call calculate_spec_vol(T15, S15, p15, a15, 1, 15, EOS, alpha_ref_mks, scale=SV_scale)
    else
      call calculate_spec_vol(T15, S15, p15, a15, 1, 15, EOS, alpha_ref_mks)
    endif

    intp(1) = dza(i,j) ; intp(5) = dza(i,j+1)
    do m=2,4
      ! Use Bode's rule to estimate the interface height anomaly change.
      ! The integrals at the ends of the segment are already known.
      pos = (m-2)*5
      intp(m) = dp_90(m) * ((7.0*(a15(pos+1)+a15(pos+5)) + &
                            32.0*(a15(pos+2)+a15(pos+4))) + 12.0*a15(pos+3))
    enddo
    ! Use Bode's rule to integrate the interface height anomaly values in x.
    inty_dza(i,J) = C1_90*((7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4))) + &
                           12.0*intp(3))
  enddo ; enddo ; endif

end subroutine int_spec_vol_dp_generic_plm

!> Convert T&S to Absolute Salinity and Conservative Temperature if using TEOS10
subroutine convert_temp_salt_for_TEOS10(T, S, HI, kd, mask_z, EOS)
  integer,               intent(in)    :: kd  !< The number of layers to work on
  type(hor_index_type),  intent(in)    :: HI       !< The horizontal index structure
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed,kd), &
                         intent(inout) :: T   !< Potential temperature referenced to the surface [degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed,kd), &
                         intent(inout) :: S   !< Salinity [ppt]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed,kd), &
                         intent(in)    :: mask_z !< 3d mask regulating which points to convert.
  type(EOS_type),        pointer       :: EOS !< Equation of state structure

  integer :: i, j, k
  real :: gsw_sr_from_sp, gsw_ct_from_pt, gsw_sa_from_sp
  real :: p

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "convert_temp_salt_to_TEOS10 called with an unassociated EOS_type EOS.")

  if ((EOS%form_of_EOS /= EOS_TEOS10) .and. (EOS%form_of_EOS /= EOS_NEMO)) return

  do k=1,kd ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
    if (mask_z(i,j,k) >= 1.0) then
     S(i,j,k) = gsw_sr_from_sp(S(i,j,k))
!     Get absolute salnity from practical salinity, converting pressures from Pascal to dbar.
!     If this option is activated, pressure will need to be added as an argument, and it should be
!     moved out into module that is not shared between components, where the ocean_grid can be used.
!     S(i,j,k) = gsw_sa_from_sp(S(i,j,k),pres(i,j,k)*1.0e-4,G%geoLonT(i,j),G%geoLatT(i,j))
     T(i,j,k) = gsw_ct_from_pt(S(i,j,k), T(i,j,k))
    endif
  enddo ; enddo ; enddo
end subroutine convert_temp_salt_for_TEOS10

!> Extractor routine for the EOS type if the members need to be accessed outside this module
subroutine extract_member_EOS(EOS, form_of_EOS, form_of_TFreeze, EOS_quadrature, Compressible, &
                              Rho_T0_S0, drho_dT, dRho_dS, TFr_S0_P0, dTFr_dS, dTFr_dp)
  type(EOS_type),    pointer     :: EOS !< Equation of state structure
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
