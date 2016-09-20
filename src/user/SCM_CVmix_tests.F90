!> Initial conditions and forcing for the single column model (SCM) CVmix
!! test set.
module SCM_CVmix_tests

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing, allocate_forcing_type
use MOM_grid, only : ocean_grid_type
use MOM_verticalgrid, only: verticalGrid_type
use MOM_safe_alloc, only : safe_alloc_ptr
use MOM_time_manager, only : time_type, operator(+), operator(/), get_time,&
                             time_type_to_real
use MOM_variables, only : thermo_var_ptrs, surface
implicit none ; private

#include <MOM_memory.h>

public SCM_CVmix_tests_TS_init
public SCM_CVmix_tests_surface_forcing_init
public SCM_CVmix_tests_wind_forcing
public SCM_CVmix_tests_buoyancy_forcing
public SCM_CVmix_tests_CS

!> Container for surface forcing parameters
type SCM_CVmix_tests_CS ;
private
  logical :: UseWindStress  !< True to use wind stress
  logical :: UseHeatFlux  !< True to use heat flux
  logical :: UseEvaporation !< True to use evaporation
  logical :: UseDiurnalSW   !< True to use diurnal sw radiation
  real :: tau_x !< (Constant) Wind stress, X (Pa)
  real :: tau_y !< (Constant) Wind stress, Y (Pa)
  real :: surf_HF !< (Constant) Heat flux (m K s^{-1})
  real :: surf_evap !< (Constant) Evaporation rate (m/s)
  real :: Max_sw !< maximum of diurnal sw radiation (m K s^{-1})
  real,public :: Rho0 !< reference density copied for easy passing
end type

! This include declares and sets the variable "version".
#include "version_variable.h"

character(len=40)  :: mod = "SCM_CVmix_tests" ! This module's name.

contains

!> Initializes temperature and salinity for the SCM CVmix test example
subroutine SCM_CVmix_tests_TS_init(T, S, h, G, GV, param_file)
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(out) :: T !< Potential tempera\ture (degC)
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(out) :: S !< Salinity (psu)
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(in)  :: h !< Layer thickness (m or Pa)
  type(ocean_grid_type),                  intent(in)  :: G !< Grid structure
  type(verticalGrid_type),                intent(in)  :: GV!< Vertical grid structure
  type(param_file_type),                  intent(in)  :: param_file !< Input parameter structure
  ! Local variables
  real :: eta(SZK_(G)+1) ! The 1-d nominal positions of the interfaces.
  real :: UpperLayerTempMLD !< Upper layer Temp MLD thickness (m)
  real :: UpperLayerSaltMLD !< Upper layer Salt MLD thickness (m)
  real :: UpperLayerTemp !< Upper layer temperature (SST if thickness 0) (deg C)
  real :: UpperLayerSalt !< Upper layer salinity (SSS if thickness 0) (PPT)
  real :: LowerLayerTemp !< Temp at top of lower layer (deg C)
  real :: LowerLayerSalt !< Salt at top of lower layer (PPT)
  real :: LowerLayerdTdz !< Temp gradient in lower layer (deg C m^{-1})
  real :: LowerLayerdSdz !< Salt gradient in lower layer (PPT m^{-1})
  real :: zC, DZ
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed


  call log_version(param_file, mod, version)
  call get_param(param_file,mod,"SCM_TEMP_MLD",UpperLayerTempMLD, &
                 'Initial temp mixed layer depth', units='m',default=0.0)
  call get_param(param_file,mod,"SCM_SALT_MLD",UpperLayerSaltMLD, &
                 'Initial salt mixed layer depth', units='m',default=0.0)
  call get_param(param_file,mod,"SCM_L1_SALT",UpperLayerSalt, &
                 'Layer 2 surface salinity', units='1e-3',default=35.0)
  call get_param(param_file,mod,"SCM_L1_TEMP",UpperLayerTemp, &
                 'Layer 1 surface temperature', units='C', default=20.0)
  call get_param(param_file,mod,"SCM_L2_SALT",LowerLayerSalt, &
                 'Layer 2 surface salinity', units='1e-3',default=35.0)
  call get_param(param_file,mod,"SCM_L2_TEMP",LowerLayerTemp, &
                 'Layer 2 surface temperature', units='C', default=20.0)
  call get_param(param_file,mod,"SCM_L2_DTDZ",LowerLayerdTdZ,     &
                 'Initial temperature stratification in layer 2', &
                 units='C/m', default=0.00)
  call get_param(param_file,mod,"SCM_L2_DSDZ",LowerLayerdSdZ,  &
                 'Initial salinity stratification in layer 2', &
                 units='PPT/m', default=0.00)

  do j=js,je ; do i=is,ie
    eta(1) = 0. ! Reference to surface
    do k=1,nz
      eta(K+1) = eta(K) - h(i,j,k)*GV%H_to_m ! Interface below layer (in m)
      zC = 0.5*( eta(K) + eta(K+1) )        ! Z of middle of layer (in m)
      DZ = min(0., zC+UpperLayerTempMLD)
      if (DZ.ge.0.0) then ! in Layer 1
        T(i,j,k) = UpperLayerTemp
      else ! in Layer 2
        T(i,j,k) = LowerLayerTemp + LowerLayerdTdZ * DZ
      endif
      DZ = min(0., zC+UpperLayerSaltMLD)
      if (DZ.ge.0.0) then ! in Layer 1
        S(i,j,k) = UpperLayerSalt
      else ! in Layer 2
        S(i,j,k) = LowerLayerSalt + LowerLayerdSdZ * DZ
      endif
    enddo ! k
  enddo ; enddo

end subroutine SCM_CVmix_tests_TS_init

!> Initializes surface forcing for the CVmix test case suite
subroutine SCM_CVmix_tests_surface_forcing_init(Time, G, param_file, CS)
  type(time_type),              intent(in) :: Time       !< Time
  type(ocean_grid_type),        intent(in) :: G          !< Grid structure
  type(param_file_type),        intent(in) :: param_file !< Input parameter structure
  type(SCM_CVmix_tests_CS), pointer    :: CS         !< Parameter container

! This include declares and sets the variable "version".
#include "version_variable.h"

  if (associated(CS)) then
    call MOM_error(FATAL, "SCM_CVmix_tests_surface_forcing_init called with an associated "// &
                          "control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "SCM_USE_WIND_STRESS",              &
                 CS%UseWindStress, "Wind Stress switch "//            &
                 "used in the SCM CVmix surface forcing.",            &
                 units='', default=.false.)
  call get_param(param_file, mod, "SCM_USE_HEAT_FLUX",                &
                 CS%UseHeatFlux, "Heat flux switch "//                &
                 "used in the SCM CVmix test surface forcing.",       &
                 units='', default=.false.)
  call get_param(param_file, mod, "SCM_USE_EVAPORATION",              &
                 CS%UseEvaporation, "Evaporation switch "//           &
                 "used in the SCM CVmix test surface forcing.",       &
                 units='', default=.false.)
  call get_param(param_file, mod, "SCM_USE_DIURNAL_SW",               &
                 CS%UseDiurnalSW, "Diurnal sw radation switch "//     &
                 "used in the SCM CVmix test surface forcing.",       &
                 units='', default=.false.)
  if (CS%UseWindStress) then
    call get_param(param_file, mod, "SCM_TAU_X",                      &
                 CS%tau_x, "Constant X-dir wind stress "//            &
                 "used in the SCM CVmix test surface forcing.",       &
                 units='N/m2', fail_if_missing=.true.)
    call get_param(param_file, mod, "SCM_TAU_Y",                      &
                 CS%tau_y, "Constant y-dir wind stress "//            &
                 "used in the SCM CVmix test surface forcing.",       &
                 units='N/m2', fail_if_missing=.true.)
  endif
  if (CS%UseHeatFlux) then
    call get_param(param_file, mod, "SCM_HEAT_FLUX",                  &
                 CS%surf_HF, "Constant surface heat flux "//          &
                 "used in the SCM CVmix test surface forcing.",       &
                 units='m K/s', fail_if_missing=.true.)
  endif
  if (CS%UseEvaporation) then
    call get_param(param_file, mod, "SCM_EVAPORATION",                &
                 CS%surf_evap, "Constant surface evaporation "//      &
                 "used in the SCM CVmix test surface forcing.",       &
                 units='m/s', fail_if_missing=.true.)
  endif
  if (CS%UseDiurnalSW) then
    call get_param(param_file, mod, "SCM_DIURNAL_SW_MAX",             &
                 CS%Max_sw, "Maximum diurnal sw radiation "//         &
                 "used in the SCM CVmix test surface forcing.",       &
                 units='m K/s', fail_if_missing=.true.)
  endif

end subroutine SCM_CVmix_tests_surface_forcing_init

subroutine SCM_CVmix_tests_wind_forcing(state, fluxes, day, G, CS)
  type(surface),                    intent(in)    :: state  !< Surface state structure
  type(forcing),                    intent(inout) :: fluxes !< Surface fluxes structure
  type(time_type),                  intent(in)    :: day    !< Time in days
  type(ocean_grid_type),            intent(inout) :: G      !< Grid structure
  type(SCM_CVmix_tests_CS), pointer       :: CS     !< Container for SCM parameters
  ! Local variables
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  real    :: mag_tau
  ! Bounds for loops and memory allocation
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ! Allocate the forcing arrays, if necessary.
  call allocate_forcing_type(G, fluxes, stress=.true., ustar=.true. )

  do j=js,je ; do I=Isq,Ieq
    fluxes%taux(I,j) = CS%tau_x
  enddo ; enddo
  do J=Jsq,Jeq ; do i=is,ie
    fluxes%tauy(i,J) = CS%tau_y
  enddo ; enddo
  mag_tau = sqrt(CS%tau_x*CS%tau_x + CS%tau_y*CS%tau_y)
  if (associated(fluxes%ustar)) then ; do j=js,je ; do i=is,ie
    fluxes%ustar(i,j) = sqrt(  mag_tau / CS%Rho0 )
  enddo ; enddo ; endif

end subroutine SCM_CVmix_tests_wind_forcing


subroutine SCM_CVmix_tests_buoyancy_forcing(state, fluxes, day, G, CS)
  type(surface),                    intent(in)    :: state  !< Surface state structure
  type(forcing),                    intent(inout) :: fluxes !< Surface fluxes structure
  type(time_type),                  intent(in)    :: day    !< Time in days (seconds?)
  type(ocean_grid_type),            intent(inout) :: G      !< Grid structure
  type(SCM_CVmix_tests_CS), pointer       :: CS     !< Container for SCM parameters

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

  ! Allocate the forcing arrays, if necessary.
  call allocate_forcing_type(G, fluxes, water=.true., heat=.true. )

  if (CS%UseHeatFlux) then
    ! Note CVmix test inputs give Heat flux in [m K/s]
    ! therefore must convert to W/m2 by multiplying
    ! by Rho0*Cp
    do J=Jsq,Jeq ; do i=is,ie
      fluxes%sens(i,J) = CS%surf_HF * CS%Rho0 * fluxes%C_p
    enddo; enddo
  endif

  if (CS%UseEvaporation) then
    do J=Jsq,Jeq ; do i=is,ie
    ! Note CVmix test inputs give evaporation in m/s
    ! This therefore must be converted to mass flux
    ! by multiplying by density
      fluxes%evap(i,J) = CS%surf_evap * CS%Rho0
    enddo; enddo
  endif

  if (CS%UseDiurnalSW) then
    do J=Jsq,Jeq ; do i=is,ie
    ! Note CVmix test inputs give max sw rad in [m K/s]
    ! therefore must convert to W/m2 by multiplying
    ! by Rho0*Cp
    ! Note diurnal cycle peaks at Noon.
      fluxes%sw(i,J) = CS%Max_sw * max(0.0,cos(2*PI*     &
           (time_type_to_real(DAY)/86400.-0.5))) * CS%RHO0 * fluxes%C_p
    enddo; enddo
  endif

end subroutine SCM_CVmix_tests_buoyancy_forcing

end module SCM_CVmix_tests
