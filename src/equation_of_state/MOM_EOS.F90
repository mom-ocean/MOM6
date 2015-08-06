!> Provides subroutines for quantities specific to the equation of state
module MOM_EOS

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_EOS_linear
use MOM_EOS_Wright
use MOM_EOS_UNESCO
use MOM_TFreeze, only : calculate_TFreeze_linear, calculate_TFreeze_Millero
use MOM_error_handler, only : MOM_error, FATAL, MOM_mesg
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_string_functions, only : uppercase
use MOM_grid, only : ocean_grid_type, ocean_block_type

implicit none ; private

#include <MOM_memory.h>

public calculate_compress, calculate_density, query_compressible
public calculate_density_derivs, calculate_2_densities
public calculate_specific_vol_derivs
public EOS_init, EOS_end, EOS_allocate
public EOS_use_linear
public int_density_dz, int_specific_vol_dp
public int_density_dz_generic_plm, int_density_dz_generic_ppm
public int_density_dz_generic_plm_analytic
public calculate_TFreeze

!> Calculates density of sea water from T, S and P
interface calculate_density
  module procedure calculate_density_scalar, calculate_density_array
end interface calculate_density

!> Calculates the freezing point of sea water from T, S and P
interface calculate_TFreeze
  module procedure calculate_TFreeze_scalar, calculate_TFreeze_array
end interface calculate_TFreeze

!> A control structure for the equation of state
type, public :: EOS_type ; private
  integer :: form_of_EOS = 0 !< The equation of state to use.
  integer :: form_of_TFreeze = 0 !< The expression for the potential temperature
                             !! of the freezing point.
  logical :: EOS_quadrature  !< If true, always use the generic (quadrature)
                             !! code for the integrals of density.
  logical :: Compressible = .true. !< If true, in situ density is a function
                             !! of pressure.
! The following parameters are used with the linear equation of state only.
  real :: Rho_T0_S0   !< The density at T=0, S=0, in kg m-3.
  real :: dRho_dT     !< The partial derivatives of density with temperature
  real :: dRho_dS     !< and salinity, in kg m-3 K-1 and kg m-3 psu-1.
! The following parameters are use with the linear expression for the freezing
! point only.
  real :: TFr_S0_P0   !< The freezing potential temperature at S=0, P=0 in deg C.
  real :: dTFr_dS !< The derivative of freezing point with salinity, in deg C PSU-1.
  real :: dTFr_dp !< The derivative of freezing point with pressure, in deg C Pa-1.
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

!> Calls the appropriate subroutine to calculate density of sea water for scalar inputs.
subroutine calculate_density_scalar(T, S, pressure, rho, EOS)
  real,           intent(in)  :: T !< Potential temperature referenced to the surface (degC)
  real,           intent(in)  :: S !< Salinity (PSU)
  real,           intent(in)  :: pressure !< Pressure (Pa)
  real,           intent(out) :: rho !< Density (in-situ if pressure is local) (kg m-3)
  type(EOS_type), pointer     :: EOS !< Equation of state structure

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_density_scalar called with an unassociated EOS_type EOS.")

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_density_scalar_linear(T, S, pressure, rho, &
                                      EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS)
    case (EOS_UNESCO)
      call calculate_density_scalar_unesco(T, S, pressure, rho)
    case (EOS_WRIGHT)
      call calculate_density_scalar_wright(T, S, pressure, rho)
    case default
      call MOM_error(FATAL, &
           "calculate_density_scalar: EOS is not valid.")
  end select

end subroutine calculate_density_scalar

!> Calls the appropriate subroutine to calculate the density of sea water for 1-D array inputs.
subroutine calculate_density_array(T, S, pressure, rho, start, npts, EOS)
  real, dimension(:), intent(in)  :: T !< Potential temperature referenced to the surface (degC)
  real, dimension(:), intent(in)  :: S !< Salinity (PSU)
  real, dimension(:), intent(in)  :: pressure !< Pressure (Pa)
  real, dimension(:), intent(out) :: rho !< Density (in-situ if pressure is local) (kg m-3)
  integer,            intent(in)  :: start !< Start index for computation
  integer,            intent(in)  :: npts !< Number of point to compute
  type(EOS_type),     pointer     :: EOS !< Equation of state structure

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

!> Calls the appropriate subroutine to calculate the freezing point for scalar inputs.
subroutine calculate_TFreeze_scalar(S, pressure, T_fr, EOS)
  real,           intent(in)  :: S !< Salinity (PSU)
  real,           intent(in)  :: pressure !< Pressure (Pa)
  real,           intent(out) :: T_fr !< Freezing point potential temperature referenced to the surface (degC)
  type(EOS_type), pointer     :: EOS !< Equation of state structure

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

!> Calls the appropriate subroutine to calculate the freezing point for a 1-D array.
subroutine calculate_TFreeze_array(S, pressure, T_fr, start, npts, EOS)
  real, dimension(:), intent(in)  :: S !< Salinity (PSU)
  real, dimension(:), intent(in)  :: pressure !< Pressure (Pa)
  real, dimension(:), intent(out) :: T_fr !< Freezing point potential temperature referenced to the surface (degC)
  integer,            intent(in)  :: start !< Starting index within the array
  integer,            intent(in)  :: npts !< The number of values to calculate
  type(EOS_type),     pointer     :: EOS !< Equation of state structure

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

!> Calls the appropriate subroutine to calculate density derivatives for 1-D array inputs.
subroutine calculate_density_derivs(T, S, pressure, drho_dT, drho_dS, start, npts, EOS)
  real, dimension(:), intent(in)  :: T !< Potential temperature referenced to the surface (degC)
  real, dimension(:), intent(in)  :: S !< Salinity (PSU)
  real, dimension(:), intent(in)  :: pressure !< Pressure (Pa)
  real, dimension(:), intent(out) :: drho_dT !< The partial derivative of density with potential tempetature, in kg m-3 K-1.
  real, dimension(:), intent(out) :: drho_dS !< The partial derivative of density with salinity, in kg m-3 psu-1.
  integer,            intent(in)  :: start !< Starting index within the array
  integer,            intent(in)  :: npts !< The number of values to calculate
  type(EOS_type),     pointer     :: EOS !< Equation of state structure

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

!> Calls the appropriate subroutine to calculate specific volume derivatives for an array.
subroutine calculate_specific_vol_derivs(T, S, pressure, dSV_dT, dSV_dS, start, npts, EOS)
  real, dimension(:), intent(in)  :: T !< Potential temperature referenced to the surface (degC)
  real, dimension(:), intent(in)  :: S !< Salinity (PSU)
  real, dimension(:), intent(in)  :: pressure !< Pressure (Pa)
  real, dimension(:), intent(out) :: dSV_dT !< The partial derivative of specific volume with potential temperature, in m3 kg-1 K-1.
  real, dimension(:), intent(out) :: dSV_dS !< The partial derivative of specific volume with salinity, in m3 kg-1 / (g/kg).
  integer,            intent(in)  :: start !< Starting index within the array
  integer,            intent(in)  :: npts !< The number of values to calculate
  type(EOS_type),     pointer     :: EOS !< Equation of state structure
  ! Local variables
  real, dimension(size(T)) :: dRho_dT, dRho_dS, rho
  integer :: j

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "calculate_density_derivs called with an unassociated EOS_type EOS.")

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
    case default
      call MOM_error(FATAL, &
           "calculate_density_derivs: EOS%form_of_EOS is not valid.")
  end select

end subroutine calculate_specific_vol_derivs

!> Calls the appropriate subroutine to calculate the density and compressibility for 1-D array inputs.
subroutine calculate_compress(T, S, pressure, rho, drho_dp, start, npts, EOS)
  real, dimension(:), intent(in)  :: T !< Potential temperature referenced to the surface (degC)
  real, dimension(:), intent(in)  :: S !< Salinity (PSU)
  real, dimension(:), intent(in)  :: pressure !< Pressure (Pa)
  real, dimension(:), intent(out) :: rho !< In situ density in kg m-3.
  real, dimension(:), intent(out) :: drho_dp !< The partial derivative of density with pressure
                                     !! (also the inverse of the square of sound speed) in s2 m-2.
  integer,            intent(in)  :: start !< Starting index within the array
  integer,            intent(in)  :: npts !< The number of values to calculate
  type(EOS_type),     pointer     :: EOS !< Equation of state structure

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

!> Calls the appropriate subroutine to calculate density for two arrays.
subroutine calculate_2_densities( T, S, pressure1, pressure2, rho1, rho2, start, npts, EOS)
  real, dimension(:), intent(in)  :: T !< Potential temperature referenced to the surface (degC)
  real, dimension(:), intent(in)  :: S !< Salinity (PSU)
  real,               intent(in)  :: pressure1 !< Pressure (Pa)
  real,               intent(in)  :: pressure2 !< A second pressure (Pa)
  real, dimension(:), intent(out) :: rho1 !< Density at pressure1, in kg m-3.
  real, dimension(:), intent(out) :: rho2 !< Density at pressure2, in kg m-3.
  integer,            intent(in)  :: start !< Starting index within the array
  integer,            intent(in)  :: npts !< The number of values to calculate
  type(EOS_type),     pointer     :: EOS !< Equation of state structure

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

!> Calls the appropriate subroutine to alculate analytical and nearly-analytical
!! integrals in pressure across layers of geopotential anomalies, which are
!! required for calculating the finite-volume form pressure accelerations in a
!! non-Boussinesq model.  There are essentially no free assumptions, apart from the
!! use of Bode's rule to do the horizontal integrals, and from a truncation in the
!! series for log(1-eps/1+eps) that assumes that |eps| <  .
subroutine int_specific_vol_dp(T, S, p_t, p_b, alpha_ref, G, EOS, &
                               dza, intp_dza, intx_dza, inty_dza, halo_size)
  !> Potential temperature referenced to the surface (degC)
    real, dimension(NIMEM_,NJMEM_),  intent(in)  :: T
  !> Salinity (PSU)
    real, dimension(NIMEM_,NJMEM_),  intent(in)  :: S
  !> Pressure at the top of the layer in Pa.
    real, dimension(NIMEM_,NJMEM_),  intent(in)  :: p_t
  !> Pressure at the bottom of the layer in Pa.
    real, dimension(NIMEM_,NJMEM_),  intent(in)  :: p_b
  !> A mean specific volume that is subtracted out to reduce the magnitude of
  !! each of the integrals, m3 kg-1. The calculation is mathematically identical
  !! with different values of alpha_ref, but this reduces the effects of roundoff.
    real,                            intent(in)  :: alpha_ref
  !> Grid structure
    type(ocean_grid_type),           intent(in)  :: G
  !> Equation of state structure
    type(EOS_type),                  pointer     :: EOS
  !> The change in the geopotential anomaly across the layer, in m2 s-2.
    real, dimension(NIMEM_,NJMEM_),  intent(out) :: dza
  !> The integral in pressure through the layer of the geopotential anomaly
  !! relative to the anomaly at the bottom of the layer, in Pa m2 s-2.
    real, dimension(NIMEM_,NJMEM_),  optional, intent(out) :: intp_dza
  !> The integral in x of the difference between the geopotential anomaly at the
  !! top and bottom of the layer divided by the x grid spacing, in m2 s-2.
    real, dimension(NIMEMB_,NJMEM_), optional, intent(out) :: intx_dza
  !> The integral in y of the difference between the geopotential anomaly at the
  !! top and bottom of the layer divided by the y grid spacing, in m2 s-2.
    real, dimension(NIMEM_,NJMEMB_), optional, intent(out) :: inty_dza
  !> The width of halo points on which to calculate dza.
    integer,                         optional, intent(in)  :: halo_size

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

!> This subroutine calculates analytical and nearly-analytical integrals of
!! pressure anomalies across layers, which are required for calculating the
!! finite-volume form pressure accelerations in a Boussinesq model.  The one
!! potentially dodgy assumtion here is that rho_0 is used both in the denominator
!! of the accelerations, and in the pressure used to calculated density (the latter
!! being -z*rho_0*G_e).  These two uses could be separated if need be.
subroutine int_density_dz(T, S, z_t, z_b, rho_ref, rho_0, G_e, B, EOS, &
                                dpa, intz_dpa, intx_dpa, inty_dpa)
  !> Potential temperature referenced to the surface (degC)
    real, dimension(NIMEM_,NJMEM_),  intent(in)  :: T
  !> Salinity (PSU)
    real, dimension(NIMEM_,NJMEM_),  intent(in)  :: S
  !> Height at the top of the layer in m.
    real, dimension(NIMEM_,NJMEM_),  intent(in)  :: z_t
  !> Height at the bottom of the layer in m.
    real, dimension(NIMEM_,NJMEM_),  intent(in)  :: z_b
  !> A mean density, in kg m-3, that is subtracted out to reduce the magnitude
  !! of each of the integrals. (The pressure is calculated as p~=-z*rho_0*G_e.)
    real,                            intent(in)  :: rho_ref
  !> A density, in kg m-3, that is used to calculate the pressure
  !! (as p~=-z*rho_0*G_e) used in the equation of state.
    real,                            intent(in)  :: rho_0
  !> The Earth's gravitational acceleration, in m s-2.
    real,                            intent(in)  :: G_e
  !> Ocean block structure.
    type(ocean_block_type),          intent(in)  :: B
  !> Equation of state structure
    type(EOS_type),                  pointer     :: EOS
  !> The change in the pressure anomaly across the layer, in Pa.
    real, dimension(NIMEM_BK_,NJMEM_BK_),            intent(out) :: dpa
  !> The integral through the thickness of the layer of the pressure anomaly
  !! relative to the anomaly at the top of the layer, in Pa m.
    real, dimension(NIMEM_BK_,NJMEM_BK_),  optional, intent(out) :: intz_dpa
  !> The integral in x of the difference between the pressure anomaly at the
  !! top and bottom of the layer divided by the x grid spacing, in Pa.
    real, dimension(NIMEMB_BK_,NJMEM_BK_), optional, intent(out) :: intx_dpa
  !> The integral in y of the difference between the pressure anomaly at the
  !! top and bottom of the layer divided by the y grid spacing, in Pa.
    real, dimension(NIMEM_BK_,NJMEMB_BK_), optional, intent(out) :: inty_dpa

  if (.not.associated(EOS)) call MOM_error(FATAL, &
    "int_density_dz called with an unassociated EOS_type EOS.")

  if (EOS%EOS_quadrature) then
    call int_density_dz_generic(T, S, z_t, z_b, rho_ref, rho_0, G_e, B, &
                                EOS, dpa, intz_dpa, intx_dpa, inty_dpa)
  else ; select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call int_density_dz_linear(T, S, z_t, z_b, rho_ref, rho_0, G_e, B, &
                                 EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS, &
                                 dpa, intz_dpa, intx_dpa, inty_dpa)
    case (EOS_WRIGHT)
      call int_density_dz_wright(T, S, z_t, z_b, rho_ref, rho_0, G_e, B, &
                                       dpa, intz_dpa, intx_dpa, inty_dpa)
    case default
      call int_density_dz_generic(T, S, z_t, z_b, rho_ref, rho_0, G_e, B, &
                                  EOS, dpa, intz_dpa, intx_dpa, inty_dpa)
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
subroutine EOS_init(param_file, EOS)
  type(param_file_type), intent(in) :: param_file !< Parameter file structure
  type(EOS_type),        pointer    :: EOS !< Equation of state structure
  ! Local variables
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_EOS" ! This module's name.
  character(len=40)  :: tmpstr

  if (.not.associated(EOS)) call EOS_allocate(EOS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version)

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

end subroutine EOS_init

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
  real,              intent(in) :: Rho_T0_S0 !< Density at T=0 degC and S=0 ppt (kg m-3)
  real,              intent(in) :: dRho_dT   !< Partial derivative of density with temperature (kg m-3 degC-1)
  real,              intent(in) :: dRho_dS   !< Partial derivative of density with salinity (kg m-3 ppt-1)
  logical, optional, intent(in) :: use_quadrature !< Partial derivative of density with salinity (kg m-3 ppt-1)
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

subroutine int_density_dz_generic(T, S, z_t, z_b, rho_ref, rho_0, G_e, B, &
                                  EOS, dpa, intz_dpa, intx_dpa, inty_dpa)
  real, dimension(NIMEM_BK_,NJMEM_BK_),  intent(in)  :: T, S, z_t, z_b
  real,                                  intent(in)  :: rho_ref, rho_0, G_e
  type(ocean_block_type),                intent(in)  :: B
  type(EOS_type),                        pointer     :: EOS !< Equation of state structure
  real, dimension(NIMEM_BK_,NJMEM_BK_),  intent(out) :: dpa
  real, dimension(NIMEM_BK_,NJMEM_BK_),  optional, intent(out) :: intz_dpa
  real, dimension(NIMEMB_BK_,NJMEM_BK_), optional, intent(out) :: intx_dpa
  real, dimension(NIMEM_BK_,NJMEMB_BK_), optional, intent(out) :: inty_dpa
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

  Isq = B%IscB ; Ieq = B%IecB ; Jsq = B%JscB ; Jeq = B%JecB

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

  if (present(intx_dpa)) then ; do j=B%jsc,B%jec ; do I=Isq,Ieq
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

  if (present(inty_dpa)) then ; do J=Jsq,Jeq ; do i=B%isc,B%iec
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


! ==============================================================================
subroutine int_density_dz_generic_cell (T_t_arg, T_b_arg, S_t_arg, S_b_arg, &
                                        z_t_arg, z_b_arg, depth, rho_ref, &
                                        rho_0, G_e, G, EOS, dpa, &
                                        intz_dpa, intx_dpa, inty_dpa)

  ! Arguments
  real, dimension(2), intent(in)     :: T_t_arg, T_b_arg, S_t_arg, S_b_arg
  real, dimension(2), intent(inout)  :: z_t_arg, z_b_arg
  real, dimension(2), intent(in)     :: depth
  real, intent(in)                   :: rho_ref, rho_0, G_e
  type(ocean_grid_type), intent(in)  :: G
  type(EOS_type), pointer            :: EOS !< Equation of state structure
  real, dimension(2), intent(out)    :: dpa
  real, dimension(2), intent(out)    :: intz_dpa
  real, intent(out)                  :: intx_dpa
  real, intent(out)                  :: inty_dpa

  ! Local variables
  real    :: T_t(2), T_b(2)             ! top and bottom temperatures
  real    :: S_t(2), S_b(2)             ! top and bottom salinities
  real    :: z_t(2), z_b(2)             ! top and bottom heights
  real    :: h1, h2                     ! cell thicknesses
  real    :: T5(5), S5(5), p5(5), r5(5) ! temperature, salinity, pressure and
                                        ! density at quadrature points
  real    :: rho_anom
  real    :: w_left, w_right, intz(5)
  real    :: GxRho
  real    :: dz
  real    :: weight_t, weight_b
  real    :: Dmin ! minimum depth
  integer :: n, m
  real, parameter :: C1_90 = 1.0/90.0 ! Rational constants

  GxRho = G_e * rho_0

  ! -------------------------------------------------------
  ! 0. Modify cell geometry to take topography into account
  ! -------------------------------------------------------
  Dmin = min ( depth(1), depth(2) )

  z_b(1) = max ( z_b_arg(1), -Dmin )
  z_b(2) = max ( z_b_arg(2), -Dmin )

  if ( z_b(1) .GT. z_t_arg(1) ) z_b(1) = z_t_arg(1)
  if ( z_b(2) .GT. z_t_arg(2) ) z_b(2) = z_t_arg(2)

  ! We do not modify the heights at the top of the cell
  z_t = z_t_arg

  h1 = z_t(1) - z_b(1)
  h2 = z_t(2) - z_b(2)

  ! Compute new salinities and temperatures at the bottom
  S_b(1) = S_t_arg(1) + (S_b_arg(1)-S_t_arg(1)) * (h1 / (z_t_arg(1)-z_b_arg(1)))
  S_b(2) = S_t_arg(2) + (S_b_arg(2)-S_t_arg(2)) * (h2 / (z_t_arg(2)-z_b_arg(2)))

  T_b(1) = T_t_arg(1) + (T_b_arg(1)-T_t_arg(1)) * (h1 / (z_t_arg(1)-z_b_arg(1)))
  T_b(2) = T_t_arg(2) + (T_b_arg(2)-T_t_arg(2)) * (h2 / (z_t_arg(2)-z_b_arg(2)))
  ! Temperatures and salinities at the top remain the same
  S_t = S_t_arg
  T_t = T_t_arg

  ! Save layer bottom interface heights for use outside this routine
  z_b_arg = z_b

  ! ----------------------------------------
  ! 1. Compute left side (vertical) integral
  ! ----------------------------------------
  dz = z_t(1) - z_b(1)

  do n = 1,5
    weight_t = 0.25 * real(5-n)
    weight_b = 1.0 - weight_t
    S5(n) = weight_t * S_t(1) + weight_b * S_b(1)
    T5(n) = weight_t * T_t(1) + weight_b * T_b(1)
    p5(n) = -GxRho*(z_t(1) - 0.25*real(n-1)*dz)
  enddo

  call calculate_density(T5, S5, p5, r5, 1, 5, EOS)

  ! Use Boole's rule to estimate the pressure anomaly change.
  rho_anom = C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + 12.0*r5(3)) - rho_ref
  dpa(1) = G_e*dz*rho_anom

  ! Use a Bode's-rule-like fifth-order accurate estimate of the
  ! double integral of the pressure anomaly.
  r5 = r5 - rho_ref
  intz_dpa(1) = 0.5*G_e*dz**2 * &
  (rho_anom - C1_90*(16.0*(r5(4)-r5(2)) + 7.0*(r5(5)-r5(1))) )

  ! -----------------------------------------
  ! 2. Compute right side (vertical) integral
  ! -----------------------------------------
  dz = z_t(2) - z_b(2)

  do n = 1,5
    weight_t = 0.25 * real(5-n)
    weight_b = 1.0 - weight_t
    S5(n) = weight_t * S_t(2) + weight_b * S_b(2)
    T5(n) = weight_t * T_t(2) + weight_b * T_b(2)
    p5(n) = -GxRho*(z_t(2) - 0.25*real(n-1)*dz)
  enddo

  call calculate_density(T5, S5, p5, r5, 1, 5, EOS)

  ! Use Boole's rule to estimate the pressure anomaly change.
  rho_anom = C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + 12.0*r5(3)) - rho_ref
  dpa(2) = G_e*dz*rho_anom

  ! Use a Bode's-rule-like fifth-order accurate estimate of the
  ! double integral of the pressure anomaly.
  r5 = r5 - rho_ref
  intz_dpa(2) = 0.5*G_e*dz**2 * &
  (rho_anom - C1_90*(16.0*(r5(4)-r5(2)) + 7.0*(r5(5)-r5(1))) )

  ! ----------------------
  ! 3. Compute x-intergral
  ! ----------------------
  !if (present(intx_dpa)) then

    intz(1) = dpa(1)
    intz(5) = dpa(2)

    do m = 2,4
      w_left = 0.25*real(5-m)
      w_right = 1.0-w_left

      dz = w_left*(z_t(1) - z_b(1)) + w_right*(z_t(2) - z_b(2))

      T5(1) = w_left*T_t(1) + w_right*T_t(2)
      T5(5) = w_left*T_b(1) + w_right*T_b(2)

      S5(1) = w_left*S_t(1) + w_right*S_t(2)
      S5(5) = w_left*S_b(1) + w_right*S_b(2)

      p5(1) = -GxRho*(w_left*z_t(1) + w_right*z_t(2))
      do n=2,5
        p5(n) = p5(n-1) + GxRho*0.25*dz
      enddo

      do n = 1,5
        weight_t = 0.25 * real(5-n)
        weight_b = 1.0 - weight_t
        S5(n) = weight_t * S5(1) + weight_b * S5(5)
      enddo

      call calculate_density (T5, S5, p5, r5, 1, 5, EOS)

      ! Use Bode's rule to estimate the pressure anomaly change.
      intz(m) = G_e*dz*( C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + &
                                12.0*r5(3)) - rho_ref)
    enddo

    ! Use Bode's rule to integrate the bottom pressure anomaly values in x.
    intx_dpa = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                      12.0*intz(3))
  !end if ! check if intx_dpa is present

  ! ---------------------
  ! 4. Compute y-intergal
  ! ---------------------
  !if (present(inty_dpa)) then

    intz(1) = dpa(1)
    intz(5) = dpa(2)

    do m=2,4
      w_left = 0.25*real(5-m)
      w_right = 1.0-w_left

      dz = w_left*(z_t(1) - z_b(2)) + w_right*(z_t(1) - z_b(2))

      S5(1) = w_left*S_t(1) + w_right*S_t(2)
      S5(5) = w_left*S_b(1) + w_right*S_b(2)
      S5(1) = w_left*S_t(1) + w_right*S_t(2)
      S5(5) = w_left*S_b(1) + w_right*S_b(2)
      p5(1) = -GxRho*(w_left*z_t(1) + w_right*z_t(2))
      do n=2,5
        p5(n) = p5(n-1) + GxRho*0.25*dz
      enddo
      do n=1,5
        weight_t = 0.25 * real(5-n)
        weight_b = 1.0 - weight_t
        S5(n) = weight_t * S5(1) + weight_b * S5(5)
      enddo
      call calculate_density(T5, S5, p5, r5, 1, 5, EOS)

      ! Use Bode's rule to estimate the pressure anomaly change.
      intz(m) = G_e*dz*( C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + &
                         12.0*r5(3)) - rho_ref)
      enddo

      ! Use Bode's rule to integrate the values.
      inty_dpa = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                        12.0*intz(3))

  !end if ! check if inty_dpa is present

end subroutine int_density_dz_generic_cell
! ============================================================================


! ==========================================================================
! Compute pressure gradient force integrals for the case where T and S
! are linear profiles.
! ==========================================================================
subroutine int_density_dz_generic_plm (T_t, T_b, S_t, S_b, z_t, z_b, rho_ref, &
                                       rho_0, G_e, H_subroundoff, bathyT, B, EOS, dpa, &
                                       intz_dpa, intx_dpa, inty_dpa, &
                                       useMassWghtInterp)
  real, dimension(NIMEM_BK_,NJMEM_BK_),  intent(in)  :: T_t, T_b, S_t, S_b, z_t, z_b
  real,                                  intent(in)  :: rho_ref, rho_0, G_e, H_subroundoff
  real, dimension(NIMEM_BK_,NJMEM_BK_),  intent(in)  :: bathyT
  type(ocean_block_type),                intent(in)  :: B
  type(EOS_type), pointer                            :: EOS !< Equation of state structure
  real, dimension(NIMEM_BK_,NJMEM_BK_),  intent(out) :: dpa
  real, dimension(NIMEM_BK_,NJMEM_BK_),  optional, intent(out) :: intz_dpa
  real, dimension(NIMEMB_BK_,NJMEM_BK_), optional, intent(out) :: intx_dpa
  real, dimension(NIMEM_BK_,NJMEMB_BK_), optional, intent(out) :: inty_dpa
  logical,                               optional, intent(in)  :: useMassWghtInterp
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
!
! Arguments: T - potential temperature relative to the surface in C
!                (the 't' and 'b' subscripts refer to the values at
!                 the top and the bottom of each layer)
!  (in)      S - salinity in PSU.    
!                (the 't' and 'b' subscripts refer to the values at
!                 the top and the bottom of each layer)
!  (in)      z_t - height at the top of the layer in m.
!  (in)      z_b - height at the top of the layer in m.
!  (in)      rho_ref - A mean density, in kg m-3, that is subtracted out to reduce
!                    the magnitude of each of the integrals.
!                    (The pressure is calucated as p~=-z*rho_0*G_e.)
!  (in)      rho_0 - A density, in kg m-3, that is used to calculate the pressure
!                    (as p~=-z*rho_0*G_e) used in the equation of state.
!  (in)      G_e - The Earth's gravitational acceleration, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      form_of_eos - integer that selects the eqn of state.
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
!  (in,opt) useMassWghtInterp - If true, uses mass weighting to interpolate
!                       T/S for top and bottom integrals.

  real :: T5((5*B%iscB+1):(5*(B%iecB+2)))
  real :: S5((5*B%iscB+1):(5*(B%iecB+2)))
  real :: p5((5*B%iscB+1):(5*(B%iecB+2)))
  real :: r5((5*B%iscB+1):(5*(B%iecB+2)))
  real :: u5((5*B%iscB+1):(5*(B%iecB+2)))
  real :: T15((15*B%iscB+1):(15*(B%iecB+1)))
  real :: S15((15*B%iscB+1):(15*(B%iecB+1)))
  real :: p15((15*B%iscB+1):(15*(B%iecB+1)))
  real :: r15((15*B%iscB+1):(15*(B%iecB+1)))
  real :: wt_t(5), wt_b(5)
  real :: rho_anom
  real :: w_left, w_right, intz(5)
  real, parameter :: C1_90 = 1.0/90.0  ! Rational constants.
  real :: GxRho, I_Rho
  real :: dz(B%iscB:B%iecB+1), dz_x(5,B%iscB:B%iecB), dz_y(5,B%isc:B%iec)
  real :: weight_t, weight_b, hWght, massWeightingToggle
  real :: Ttl, Tbl, Ttr, Tbr, Stl, Sbl, Str, Sbr, hL, hR, iDenom
  integer :: Isq, Ieq, Jsq, Jeq, i, j, m, n
  integer :: pos

  Isq = B%IscB ; Ieq = B%IecB ; Jsq = B%JscB ; Jeq = B%JecB

  GxRho = G_e * rho_0
  I_Rho = 1.0 / rho_0
  massWeightingToggle = 0.
  if (present(useMassWghtInterp)) then
    if (useMassWghtInterp) massWeightingToggle = 1.
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
    call calculate_density_array(T5, S5, p5, r5, 1, (ieq-isq+2)*5, EOS )
    u5 = r5 - rho_ref

    do i = isq, ieq+1
    ! Use Bode's rule to estimate the pressure anomaly change.
      rho_anom = C1_90*(7.0*(r5(i*5+1)+r5(i*5+5)) + 32.0*(r5(i*5+2)+r5(i*5+4)) + 12.0*r5(i*5+3)) - &
           rho_ref
      dpa(i,j) = G_e*dz(i)*rho_anom
      if (present(intz_dpa)) then
      ! Use a Bode's-rule-like fifth-order accurate estimate of
      ! the double integral of the pressure anomaly.
        intz_dpa(i,j) = 0.5*G_e*dz(i)**2 * &
                (rho_anom - C1_90*(16.0*(u5(i*5+4)-u5(i*5+2)) + 7.0*(u5(i*5+5)-u5(i*5+1))) )
      endif
    enddo
  enddo ! end loops on j


  ! ==================================================
  ! 2. Compute horizontal integrals in the x direction
  ! ==================================================
  if (present(intx_dpa)) then ; do j=B%jsc,B%jec
     do I=Isq,Ieq

    ! Corner values of T and S
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation
    ! of T,S along the top and bottom integrals, almost like thickness
    ! weighting.
    ! Note: To work in terrain following coordinates we could offset
    ! this distance by the layer thickness to replicate other models.
    hWght = massWeightingToggle * &
            max(0., -bathyT(i,j)-z_t(i+1,j), -bathyT(i+1,j)-z_t(i,j))
    if (hWght > 0.) then
      hL = (z_t(i,j) - z_b(i,j)) + H_subroundoff
      hR = (z_t(i+1,j) - z_b(i+1,j)) + H_subroundoff
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

    call calculate_density(T15, S15, p15, r15, 1, 15*(ieq-isq+1), EOS)

    do I=Isq,Ieq
       intz(1) = dpa(i,j) ; intz(5) = dpa(i+1,j)

    ! Use Bode's rule to estimate the pressure anomaly change.
    do m = 2,4
          pos = i*15+(m-2)*5
          intz(m) = G_e*dz_x(m,i)*( C1_90*(7.0*(r15(pos+1)+r15(pos+5)) + 32.0*(r15(pos+2)+r15(pos+4)) + &
                            12.0*r15(pos+3)) - rho_ref)
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
    do i=B%isc,B%iec
    ! Corner values of T and S
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation
    ! of T,S along the top and bottom integrals, almost like thickness
    ! weighting.
    ! Note: To work in terrain following coordinates we could offset
    ! this distance by the layer thickness to replicate other models.
      hWght = massWeightingToggle * &
              max(0., -bathyT(i,j)-z_t(i,j+1), -bathyT(i,j+1)-z_t(i,j))
      if (hWght > 0.) then
        hL = (z_t(i,j) - z_b(i,j)) + H_subroundoff
        hR = (z_t(i,j+1) - z_b(i,j+1)) + H_subroundoff
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

    call calculate_density_array(T15(15*B%isc+1:), S15(15*B%isc+1:), p15(15*B%isc+1:), &
                                 r15(15*B%isc+1:), 1, 15*(B%iec-B%isc+1), EOS)
    do i=B%isc,B%iec
      intz(1) = dpa(i,j) ; intz(5) = dpa(i,j+1)

      ! Use Bode's rule to estimate the pressure anomaly change.
      do m = 2,4
        pos = i*15+(m-2)*5
        intz(m) = G_e*dz_y(m,i)*( C1_90*(7.0*(r15(pos+1)+r15(pos+5)) + &
                                         32.0*(r15(pos+2)+r15(pos+4)) + &
                                         12.0*r15(pos+3)) - rho_ref)
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

! ==========================================================================
! Compute pressure gradient force integrals for the case where T and S
! are parabolic profiles
! ==========================================================================
subroutine int_density_dz_generic_ppm (T, T_t, T_b, S, S_t, S_b, &
                                       z_t, z_b, rho_ref, rho_0, G_e, B, &
                                       EOS, dpa, intz_dpa, intx_dpa, inty_dpa)
  
  real, dimension(NIMEM_BK_,NJMEM_BK_),  intent(in)  :: T, T_t, T_b, S, S_t, S_b, &
                                                        z_t, z_b
  real,                                  intent(in)  :: rho_ref, rho_0, G_e
  type(ocean_block_type),                intent(in)  :: B
  type(EOS_type), pointer                            :: EOS !< Equation of state structure
  real, dimension(NIMEM_BK_,NJMEM_BK_),  intent(out) :: dpa
  real, dimension(NIMEM_BK_,NJMEM_BK_),  optional, intent(out) :: intz_dpa
  real, dimension(NIMEMB_BK_,NJMEM_BK_), optional, intent(out) :: intx_dpa
  real, dimension(NIMEM_BK_,NJMEMB_BK_), optional, intent(out) :: inty_dpa
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
!
! Arguments: T - potential temperature relative to the surface in C
!                (the 't' and 'b' subscripts refer to the values at
!                 the top and the bottom of each layer)
!  (in)      S - salinity in PSU.    
!                (the 't' and 'b' subscripts refer to the values at
!                 the top and the bottom of each layer)
!  (in)      z_t - height at the top of the layer in m.
!  (in)      z_b - height at the top of the layer in m.
!  (in)      rho_ref - A mean density, in kg m-3, that is subtracted out to reduce
!                    the magnitude of each of the integrals.
!                    (The pressure is calucated as p~=-z*rho_0*G_e.)
!  (in)      rho_0 - A density, in kg m-3, that is used to calculate the pressure
!                    (as p~=-z*rho_0*G_e) used in the equation of state.
!  (in)      G_e - The Earth's gravitational acceleration, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      form_of_eos - integer that selects the eqn of state.
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
  real :: weight_t, weight_b
  real :: s0, s1, s2                   ! parabola coefficients for S
  real :: t0, t1, t2                   ! parabola coefficients for T
  real :: xi                           ! normalized coordinate
  real :: T_top, T_mid, T_bot
  real :: S_top, S_mid, S_bot
  integer :: Isq, Ieq, Jsq, Jeq, i, j, m, n
  real, dimension(4) :: x, y
  real, dimension(9) :: S_node, T_node, p_node, r_node


  call MOM_error(FATAL, &
    "int_density_dz_generic_ppm: the implementation is not done yet, contact developer")

  Isq = B%IscB ; Ieq = B%IecB ; Jsq = B%JscB ; Jeq = B%JecB

  GxRho = G_e * rho_0
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

    call calculate_density(T5, S5, p5, r5, 1, 5, EOS)

    ! Use Bode's rule to estimate the pressure anomaly change.
    !rho_anom = C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + 12.0*r5(3)) - &
    !       rho_ref

    rho_anom = 1000.0 + S(i,j) - rho_ref;
    dpa(i,j) = G_e*dz*rho_anom

    ! Use a Bode's-rule-like fifth-order accurate estimate of
    ! the double integral of the pressure anomaly.
    !r5 = r5 - rho_ref
    !if (present(intz_dpa)) intz_dpa(i,j) = 0.5*G_e*dz**2 * &
    !      (rho_anom - C1_90*(16.0*(r5(4)-r5(2)) + 7.0*(r5(5)-r5(1))) )

    intz_dpa(i,j) = 0.5 * G_e * dz**2 * ( 1000.0 - rho_ref + s0 + s1/3.0 + &
                                    s2/6.0 )
  enddo ; enddo ! end loops on j and i

  ! ==================================================
  ! 2. Compute horizontal integrals in the x direction
  ! ==================================================
  if (present(intx_dpa)) then ; do j=B%jsc,B%jec ; do I=Isq,Ieq
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

      call calculate_density(T5, S5, p5, r5, 1, 5, EOS)

    ! Use Bode's rule to estimate the pressure anomaly change.
      intz(m) = G_e*dz*( C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + &
                            12.0*r5(3)) - rho_ref)
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

    call calculate_density ( T_node, S_node, p_node, r_node, 1, 9, EOS )
    r_node = r_node - rho_ref

    call compute_integral_quadratic ( x, y, r_node, intx_dpa(i,j) )

    intx_dpa(i,j) = intx_dpa(i,j) * G_e

  enddo ; enddo ; endif

  ! ==================================================
  ! 3. Compute horizontal integrals in the y direction
  ! ==================================================
  if (present(inty_dpa)) then ; do J=Jsq,Jeq ; do i=B%isc,B%iec

    inty_dpa(i,j) = 0.0

  enddo ; enddo ; endif

end subroutine int_density_dz_generic_ppm


! ==========================================================================
! Compute pressure gradient force integrals for the case where T and S
! are linear profiles (analytical !!)
! ==========================================================================
subroutine int_density_dz_generic_plm_analytic (T_t, T_b, S_t, S_b, z_t, &
            z_b, rho_ref, rho_0, G_e, G, EOS, dpa, intz_dpa, intx_dpa, inty_dpa)
  
  real, dimension(NIMEM_,NJMEM_),  intent(in)  :: T_t, T_b, S_t, S_b, z_t, z_b
  real,                            intent(in)  :: rho_ref, rho_0, G_e
  type(ocean_grid_type),           intent(in)  :: G
  type(EOS_type), pointer                      :: EOS !< Equation of state structure
  real, dimension(NIMEM_,NJMEM_),  intent(out) :: dpa
  real, dimension(NIMEM_,NJMEM_),  optional, intent(out) :: intz_dpa
  real, dimension(NIMEMB_,NJMEM_), optional, intent(out) :: intx_dpa
  real, dimension(NIMEM_,NJMEMB_), optional, intent(out) :: inty_dpa
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
!
! Arguments: T - potential temperature relative to the surface in C
!                (the 't' and 'b' subscripts refer to the values at
!                 the top and the bottom of each layer)
!  (in)      S - salinity in PSU.    
!                (the 't' and 'b' subscripts refer to the values at
!                 the top and the bottom of each layer)
!  (in)      z_t - height at the top of the layer in m.
!  (in)      z_b - height at the top of the layer in m.
!  (in)      rho_ref - A mean density, in kg m-3, that is subtracted out to reduce
!                    the magnitude of each of the integrals.
!                    (The pressure is calucated as p~=-z*rho_0*G_e.)
!  (in)      rho_0 - A density, in kg m-3, that is used to calculate the pressure
!                    (as p~=-z*rho_0*G_e) used in the equation of state.
!  (in)      G_e - The Earth's gravitational acceleration, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      form_of_eos - integer that selects the eqn of state.
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
  real :: dz, dS
  real :: weight_t, weight_b
  integer :: Isq, Ieq, Jsq, Jeq, i, j, m, n

  real, dimension(4) :: x, y, f

  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  GxRho = G_e * rho_0
  I_Rho = 1.0 / rho_0

  ! =============================
  ! 1. Compute vertical integrals
  ! =============================
  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    dz = z_t(i,j) - z_b(i,j)
    dS = S_t(i,j) - S_b(i,j)
    do n=1,5
      weight_t = 0.25 * real(5-n)
      weight_b = 1.0 - weight_t

      p5(n) = -GxRho*(z_t(i,j) - 0.25*real(n-1)*dz)

      ! Salinity and temperature points are linearly interpolated
      S5(n) = weight_t * S_t(i,j) + weight_b * S_b(i,j)
      T5(n) = weight_t * T_t(i,j) + weight_b * T_b(i,j)
    enddo

    call calculate_density(T5, S5, p5, r5, 1, 5, EOS)

    ! Use Bode's rule to estimate the pressure anomaly change.
    rho_anom = C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + 12.0*r5(3)) - &
           rho_ref
    dpa(i,j) = G_e*dz*rho_anom

    ! Pressure anomaly change (computed analytically based on linear EOS)
    rho_anom = 1000.0 + S_b(i,j) + 0.5 * dS - rho_ref
    dpa(i,j) = G_e * dz * rho_anom

    ! Use a Bode's-rule-like fifth-order accurate estimate of
    ! the double integral of the pressure anomaly.
    r5 = r5 - rho_ref
    if (present(intz_dpa)) intz_dpa(i,j) = 0.5*G_e*dz**2 * &
          (rho_anom - C1_90*(16.0*(r5(4)-r5(2)) + 7.0*(r5(5)-r5(1))) )

    intz_dpa(i,j) = ( 0.5 * (S_b(i,j)+1000.0-rho_ref) + &
                      (1.0/3.0) * dS ) * G_e * dz**2

  enddo ; enddo ! end loops on j and i

  ! ==================================================
  ! 2. Compute horizontal integrals in the x direction
  ! ==================================================
  if (present(intx_dpa)) then ; do j=G%jsc,G%jec ; do I=Isq,Ieq

    ! Use Gauss quadrature rule to compute integral
    x(1) = 1.0
    x(2) = 0.0
    x(3) = 0.0
    x(4) = 1.0
    y(1) = z_t(i+1,j)
    y(2) = z_t(i,j)
    y(3) = z_b(i,j)
    y(4) = z_b(i+1,j)
    f(1) = 1000.0 + S_t(i+1,j) - rho_ref
    f(2) = 1000.0 + S_t(i,j) - rho_ref
    f(3) = 1000.0 + S_b(i,j) - rho_ref
    f(4) = 1000.0 + S_b(i+1,j) - rho_ref

    call compute_integral_bilinear ( x, y, f, intx_dpa(i,j) )

    intx_dpa(i,j) = intx_dpa(i,j) * G_e

  enddo ; enddo ; endif

  ! ==================================================
  ! 3. Compute horizontal integrals in the y direction
  ! ==================================================
  if (present(inty_dpa)) then ; do J=Jsq,Jeq ; do i=G%isc,G%iec

    ! Use Gauss quadrature rule to compute integral
    x(1) = 1.0
    x(2) = 0.0
    x(3) = 0.0
    x(4) = 1.0
    y(1) = z_t(i,j+1)
    y(2) = z_t(i,j)
    y(3) = z_b(i,j)
    y(4) = z_b(i,j+1)
    f(1) = 1000.0 + S_t(i,j+1) - rho_ref
    f(2) = 1000.0 + S_t(i,j) - rho_ref
    f(3) = 1000.0 + S_b(i,j) - rho_ref
    f(4) = 1000.0 + S_b(i,j+1) - rho_ref

    call compute_integral_bilinear ( x, y, f, inty_dpa(i,j) )

    inty_dpa(i,j) = inty_dpa(i,j) * G_e

  enddo ; enddo ; endif

end subroutine int_density_dz_generic_plm_analytic


! =============================================================================
! Compute integral of bilinear function
! =============================================================================
subroutine compute_integral_bilinear ( x, y, f, integral )

  ! Arguments
  real, intent(in), dimension(4)    :: x, y, f
  real, intent(out)                 :: integral

  ! Local variables
  integer               :: i, k
  real, dimension(4)    :: weight, xi, eta          ! integration points
  real                  :: f_k
  real                  :: dxdxi, dxdeta
  real                  :: dydxi, dydeta
  real, dimension(4)    :: phi, dphidxi, dphideta
  real                  :: jacobian_k

  ! Quadrature rule
  weight(:) = 1.0
  xi(1) = - sqrt(3.0) / 3.0
  xi(2) = sqrt(3.0) / 3.0
  xi(3) = sqrt(3.0) / 3.0
  xi(4) = - sqrt(3.0) / 3.0
  eta(1) = - sqrt(3.0) / 3.0
  eta(2) = - sqrt(3.0) / 3.0
  eta(3) = sqrt(3.0) / 3.0
  eta(4) = sqrt(3.0) / 3.0

  integral = 0.0

  ! Integration loop
  do k = 1,4

    ! Evaluate shape functions and gradients
    call evaluate_shape_bilinear ( xi(k), eta(k), phi, dphidxi, dphideta )

    ! Determine gradient of global coordinate at integration point
    dxdxi  = 0.0
    dxdeta = 0.0
    dydxi  = 0.0
    dydeta = 0.0

    do i = 1,4
      dxdxi  = dxdxi  + x(i) * dphidxi(i)
      dxdeta = dxdeta + x(i) * dphideta(i)
      dydxi  = dydxi  + y(i) * dphidxi(i)
      dydeta = dydeta + y(i) * dphideta(i)
    enddo

    ! Evaluate Jacobian at integration point
    jacobian_k = dxdxi*dydeta - dydxi*dxdeta

    ! Evaluate function at integration point
    f_k = 0.0
    do i = 1,4
      f_k = f_k + f(i) * phi(i)
    enddo

    integral = integral + weight(k) * f_k * jacobian_k

  enddo ! end integration loop

end subroutine compute_integral_bilinear


! =============================================================================
! Compute integral of quadratic function
! =============================================================================
subroutine compute_integral_quadratic ( x, y, f, integral )

  ! Arguments
  real, intent(in), dimension(4)    :: x, y
  real, intent(in), dimension(9)    :: f
  real, intent(out)                 :: integral

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
    call evaluate_shape_bilinear ( xi(k), eta(k), phiiso, &
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
    call evaluate_shape_quadratic ( xi(k), eta(k), phi, dphidxi, dphideta )

    ! Evaluate function at integration point
    f_k = 0.0
    do i = 1,9
      f_k = f_k + f(i) * phi(i)
    enddo

    integral = integral + weight(k) * f_k * jacobian_k

  enddo ! end integration loop

end subroutine compute_integral_quadratic


! =============================================================================
! Evaluation of the four bilinear shape fn and their gradients at (xi,eta)
! =============================================================================
subroutine evaluate_shape_bilinear ( xi, eta, phi, dphidxi, dphideta )

  ! Arguments
  real, intent(in)                  :: xi, eta
  real, dimension(4), intent(out)   :: phi, dphidxi, dphideta

  ! The shape functions within the parent element are defined as shown
  ! here:
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
! Evaluation of the nine quadratic shape fn and their gradients at (xi,eta)
! =============================================================================
subroutine evaluate_shape_quadratic ( xi, eta, phi, dphidxi, dphideta )

  ! Arguments
  real, intent(in)                  :: xi, eta
  real, dimension(9), intent(out)   :: phi, dphidxi, dphideta

  ! The quadratic shape functions within the parent element are
  ! defined as shown here:
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

  phi      = 0.0
  dphidxi  = 0.0
  dphideta = 0.0

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

subroutine int_spec_vol_dp_generic(T, S, p_t, p_b, alpha_ref, G, EOS, &
                                   dza, intp_dza, intx_dza, inty_dza, halo_size)
  real, dimension(NIMEM_,NJMEM_),  intent(in)  :: T, S, p_t, p_b
  real,                            intent(in)  :: alpha_ref
  type(ocean_grid_type),           intent(in)  :: G
  type(EOS_type),                  pointer     :: EOS !< Equation of state structure
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

  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
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

!> \namespace mom_eos
!!
!! The MOM_EOS module is a wrapper for various equations of state (e.g. Linear,
!! Wright, UNESCO) and provides a uniform interface to the rest of the model
!! independent of which equation of state is being used.
