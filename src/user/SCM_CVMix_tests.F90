!> Initial conditions and forcing for the single column model (SCM) CVMix test set.
module SCM_CVMix_tests

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domains,       only : pass_var, pass_vector, TO_ALL
use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_forcing_type,  only : forcing, mech_forcing
use MOM_grid,          only : ocean_grid_type
use MOM_verticalgrid,  only : verticalGrid_type
use MOM_safe_alloc,    only : safe_alloc_ptr
use MOM_unit_scaling,  only : unit_scale_type
use MOM_time_manager,  only : time_type, operator(+), operator(/), time_type_to_real
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs, surface

implicit none ; private

#include <MOM_memory.h>

public SCM_CVMix_tests_TS_init
public SCM_CVMix_tests_surface_forcing_init
public SCM_CVMix_tests_wind_forcing
public SCM_CVMix_tests_buoyancy_forcing
public SCM_CVMix_tests_CS

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Container for surface forcing parameters
type SCM_CVMix_tests_CS ; private
  logical :: UseWindStress  !< True to use wind stress
  logical :: UseHeatFlux    !< True to use heat flux
  logical :: UseEvaporation !< True to use evaporation
  logical :: UseDiurnalSW   !< True to use diurnal sw radiation
  real :: tau_x !< (Constant) Wind stress, X [Pa]
  real :: tau_y !< (Constant) Wind stress, Y [Pa]
  real :: surf_HF !< (Constant) Heat flux [degC Z T-1 ~> m degC s-1]
  real :: surf_evap !< (Constant) Evaporation rate [Z T-1 ~> m s-1]
  real :: Max_sw !< maximum of diurnal sw radiation [degC Z T-1 ~> degC m s-1]
  real :: Rho0 !< reference density [R ~> kg m-3]
end type

! This include declares and sets the variable "version".
#include "version_variable.h"

character(len=40)  :: mdl = "SCM_CVMix_tests" !< This module's name.

contains

!> Initializes temperature and salinity for the SCM CVMix test example
subroutine SCM_CVMix_tests_TS_init(T, S, h, G, GV, US, param_file, just_read_params)
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(out) :: T  !< Potential temperature [degC]
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(out) :: S  !< Salinity [psu]
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(in)  :: h  !< Layer thickness [H ~> m or kg m-2]
  type(ocean_grid_type),                  intent(in)  :: G  !< Grid structure
  type(verticalGrid_type),                intent(in)  :: GV !< Vertical grid structure
  type(unit_scale_type),                  intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                  intent(in)  :: param_file !< Input parameter structure
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.
  ! Local variables
  real :: UpperLayerTempMLD !< Upper layer Temp MLD thickness [Z ~> m].
  real :: UpperLayerSaltMLD !< Upper layer Salt MLD thickness [Z ~> m].
  real :: UpperLayerTemp !< Upper layer temperature (SST if thickness 0) [degC]
  real :: UpperLayerSalt !< Upper layer salinity (SSS if thickness 0) [ppt]
  real :: LowerLayerTemp !< Temp at top of lower layer [degC]
  real :: LowerLayerSalt !< Salt at top of lower layer [ppt]
  real :: LowerLayerdTdz !< Temp gradient in lower layer [degC / Z ~> degC m-1].
  real :: LowerLayerdSdz !< Salt gradient in lower layer [ppt / Z ~> ppt m-1].
  real :: LowerLayerMinTemp !< Minimum temperature in lower layer [degC]
  real :: zC, DZ, top, bottom ! Depths and thicknesses [Z ~> m].
  logical :: just_read    ! If true, just read parameters but set nothing.
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed


  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) call log_version(param_file, mdl, version)
  call get_param(param_file, mdl, "SCM_TEMP_MLD", UpperLayerTempMLD, &
                 'Initial temp mixed layer depth', &
                 units='m', default=0.0, scale=US%m_to_Z, do_not_log=just_read)
  call get_param(param_file, mdl, "SCM_SALT_MLD", UpperLayerSaltMLD, &
                 'Initial salt mixed layer depth', &
                 units='m', default=0.0, scale=US%m_to_Z, do_not_log=just_read)
  call get_param(param_file, mdl, "SCM_L1_SALT", UpperLayerSalt, &
                 'Layer 2 surface salinity', units='1e-3', default=35.0, do_not_log=just_read)
  call get_param(param_file, mdl, "SCM_L1_TEMP", UpperLayerTemp, &
                 'Layer 1 surface temperature', units='C', default=20.0, do_not_log=just_read)
  call get_param(param_file, mdl, "SCM_L2_SALT", LowerLayerSalt, &
                 'Layer 2 surface salinity', units='1e-3', default=35.0, do_not_log=just_read)
  call get_param(param_file, mdl, "SCM_L2_TEMP", LowerLayerTemp, &
                 'Layer 2 surface temperature', units='C', default=20.0, do_not_log=just_read)
  call get_param(param_file, mdl, "SCM_L2_DTDZ", LowerLayerdTdZ,     &
                 'Initial temperature stratification in layer 2', &
                 units='C/m', default=0.0, scale=US%Z_to_m, do_not_log=just_read)
  call get_param(param_file, mdl, "SCM_L2_DSDZ", LowerLayerdSdZ,  &
                 'Initial salinity stratification in layer 2', &
                 units='PPT/m', default=0.0, scale=US%Z_to_m, do_not_log=just_read)
  call get_param(param_file, mdl, "SCM_L2_MINTEMP",LowerLayerMinTemp, &
                 'Layer 2 minimum temperature', units='C', default=4.0, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  do j=js,je ; do i=is,ie
    top = 0. ! Reference to surface
    bottom = 0.
    do k=1,nz
      bottom = bottom - h(i,j,k)*GV%H_to_Z ! Interface below layer [Z ~> m]
      zC = 0.5*( top + bottom )        ! Z of middle of layer [Z ~> m]
      DZ = min(0., zC + UpperLayerTempMLD)
      T(i,j,k) = max(LowerLayerMinTemp,LowerLayerTemp + LowerLayerdTdZ * DZ)
      DZ = min(0., zC + UpperLayerSaltMLD)
      S(i,j,k) = LowerLayerSalt + LowerLayerdSdZ * DZ
      top = bottom
    enddo ! k
  enddo ; enddo

end subroutine SCM_CVMix_tests_TS_init

!> Initializes surface forcing for the CVMix test case suite
subroutine SCM_CVMix_tests_surface_forcing_init(Time, G, param_file, CS)
  type(time_type),          intent(in) :: Time       !< Model time
  type(ocean_grid_type),    intent(in) :: G          !< Grid structure
  type(param_file_type),    intent(in) :: param_file !< Input parameter structure
  type(SCM_CVMix_tests_CS), pointer    :: CS         !< Parameter container


  ! This include declares and sets the variable "version".
# include "version_variable.h"
  type(unit_scale_type), pointer :: US => NULL() !< A dimensional unit scaling type

  US => G%US

  if (associated(CS)) then
    call MOM_error(FATAL, "SCM_CVMix_tests_surface_forcing_init called with an associated "// &
                          "control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "SCM_USE_WIND_STRESS",              &
                 CS%UseWindStress, "Wind Stress switch "//            &
                 "used in the SCM CVMix surface forcing.",            &
                 units='', default=.false.)
  call get_param(param_file, mdl, "SCM_USE_HEAT_FLUX",                &
                 CS%UseHeatFlux, "Heat flux switch "//                &
                 "used in the SCM CVMix test surface forcing.",       &
                 units='', default=.false.)
  call get_param(param_file, mdl, "SCM_USE_EVAPORATION",              &
                 CS%UseEvaporation, "Evaporation switch "//           &
                 "used in the SCM CVMix test surface forcing.",       &
                 units='', default=.false.)
  call get_param(param_file, mdl, "SCM_USE_DIURNAL_SW",               &
                 CS%UseDiurnalSW, "Diurnal sw radation switch "//     &
                 "used in the SCM CVMix test surface forcing.",       &
                 units='', default=.false.)
  if (CS%UseWindStress) then
    call get_param(param_file, mdl, "SCM_TAU_X",                      &
                 CS%tau_x, "Constant X-dir wind stress "//            &
                 "used in the SCM CVMix test surface forcing.",       &
                 units='N/m2', scale=US%kg_m2s_to_RZ_T*US%m_s_to_L_T, fail_if_missing=.true.)
    call get_param(param_file, mdl, "SCM_TAU_Y",                      &
                 CS%tau_y, "Constant y-dir wind stress "//            &
                 "used in the SCM CVMix test surface forcing.",       &
                 units='N/m2', scale=US%kg_m2s_to_RZ_T*US%m_s_to_L_T, fail_if_missing=.true.)
  endif
  if (CS%UseHeatFlux) then
    call get_param(param_file, mdl, "SCM_HEAT_FLUX",                  &
                 CS%surf_HF, "Constant surface heat flux "//          &
                 "used in the SCM CVMix test surface forcing.",       &
                 units='m K/s', scale=US%m_to_Z*US%T_to_s, fail_if_missing=.true.)
  endif
  if (CS%UseEvaporation) then
    call get_param(param_file, mdl, "SCM_EVAPORATION",                &
                 CS%surf_evap, "Constant surface evaporation "//      &
                 "used in the SCM CVMix test surface forcing.",       &
                 units='m/s', scale=US%m_to_Z*US%T_to_s, fail_if_missing=.true.)
  endif
  if (CS%UseDiurnalSW) then
    call get_param(param_file, mdl, "SCM_DIURNAL_SW_MAX",             &
                 CS%Max_sw, "Maximum diurnal sw radiation "//         &
                 "used in the SCM CVMix test surface forcing.",       &
                 units='m K/s', scale=US%m_to_Z*US%T_to_s, fail_if_missing=.true.)
  endif
  call get_param(param_file, mdl, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINSEQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%kg_m3_to_R)

end subroutine SCM_CVMix_tests_surface_forcing_init

subroutine SCM_CVMix_tests_wind_forcing(sfc_state, forces, day, G, US, CS)
  type(surface),            intent(in)    :: sfc_state  !< Surface state structure
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),          intent(in)    :: day    !< Time in days
  type(ocean_grid_type),    intent(inout) :: G      !< Grid structure
  type(unit_scale_type),    intent(in)    :: US     !< A dimensional unit scaling type
  type(SCM_CVMix_tests_CS), pointer       :: CS     !< Container for SCM parameters
  ! Local variables
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  real    :: mag_tau
  ! Bounds for loops and memory allocation
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  do j=js,je ; do I=Isq,Ieq
    forces%taux(I,j) = CS%tau_x
  enddo ; enddo
  do J=Jsq,Jeq ; do i=is,ie
    forces%tauy(i,J) = CS%tau_y
  enddo ; enddo
  call pass_vector(forces%taux, forces%tauy, G%Domain, To_All)

  mag_tau = sqrt(CS%tau_x*CS%tau_x + CS%tau_y*CS%tau_y)
  if (associated(forces%ustar)) then ; do j=js,je ; do i=is,ie
    forces%ustar(i,j) = sqrt( US%L_to_Z * mag_tau / (CS%Rho0) )
  enddo ; enddo ; endif

end subroutine SCM_CVMix_tests_wind_forcing


subroutine SCM_CVMix_tests_buoyancy_forcing(sfc_state, fluxes, day, G, US, CS)
  type(surface),            intent(in)    :: sfc_state  !< Surface state structure
  type(forcing),            intent(inout) :: fluxes !< Surface fluxes structure
  type(time_type),          intent(in)    :: day    !< Current model time
  type(ocean_grid_type),    intent(inout) :: G      !< Grid structure
  type(unit_scale_type),    intent(in)    :: US     !< A dimensional unit scaling type
  type(SCM_CVMix_tests_CS), pointer       :: CS     !< Container for SCM parameters

  ! Local variables
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  real :: PI

  PI = 4.0*atan(1.0)

  ! Bounds for loops and memory allocation
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (CS%UseHeatFlux) then
    ! Note CVMix test inputs give Heat flux in [m K/s]
    ! therefore must convert to W/m2 by multiplying
    ! by Rho0*Cp
    do J=Jsq,Jeq ; do i=is,ie
      fluxes%sens(i,J) = CS%surf_HF * CS%Rho0 * fluxes%C_p
    enddo ; enddo
  endif

  if (CS%UseEvaporation) then
    do J=Jsq,Jeq ; do i=is,ie
    ! Note CVMix test inputs give evaporation in [m s-1]
    ! This therefore must be converted to mass flux in [R Z T-1 ~> kg m-2 s-1]
    ! by multiplying by density and some unit conversion factors.
      fluxes%evap(i,J) = CS%surf_evap * CS%Rho0
    enddo ; enddo
  endif

  if (CS%UseDiurnalSW) then
    do J=Jsq,Jeq ; do i=is,ie
    ! Note CVMix test inputs give max sw rad in [m degC/s]
    ! therefore must convert to W/m2 by multiplying by Rho0*Cp
    ! Note diurnal cycle peaks at Noon.
      fluxes%sw(i,J) = CS%Max_sw *  max(0.0, cos(2*PI*(time_type_to_real(DAY)/86400.0 - 0.5))) * CS%RHO0 * fluxes%C_p
    enddo ; enddo
  endif

end subroutine SCM_CVMix_tests_buoyancy_forcing

end module SCM_CVMix_tests
