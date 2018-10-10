!> Initial conditions and forcing for the idealized hurricane example.
module Idealized_hurricane
! Renamed from SCM_idealized_hurricane to idealizeD_hurricane
!  This module is no longer exclusively for use in SCM mode.

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing, mech_forcing
use MOM_forcing_type, only : allocate_forcing_type, allocate_mech_forcing
use MOM_grid, only : ocean_grid_type
use MOM_safe_alloc, only : safe_alloc_ptr
use MOM_time_manager, only : time_type, operator(+), operator(/), time_type_to_real
use MOM_variables, only : thermo_var_ptrs, surface
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public idealized_hurricane_TS_init !Public interface to initialize TS as vertically
                                   ! uniform with prescribed vertical for hurricane
                                   ! experiments.  Used for other idealized
                                   ! configurations.
public idealized_hurricane_wind_init !Public interface to intialize the idealized
                                     ! hurricane wind profile.
public idealized_hurricane_wind_forcing !Public interface to update the idealized
                                        ! hurricane wind profile.

!> Container for parameters describing idealized wind structure
type, public :: idealized_hurricane_CS ; private

  ! Parameters used to compute Holland radial wind profile
  real    :: rho_a                !< Mean air density [kg/m3]
  real    :: pressure_ambient     !< Pressure at surface of ambient air [Pa]
  real    :: pressure_central     !< Pressure at surface at hurricane center [Pa]
  real    :: rad_max_wind         !< Radius of maximum winds [m]
  real    :: max_windspeed        !< Maximum wind speeds [m/s]
  real    :: hurr_translation_spd !< Hurricane translation speed [m/s]
  real    :: hurr_translation_dir !< Hurricane translation speed [m/s]
  real    :: gustiness            !< Gustiness (optional, used in u*) [m/s]
  real    :: Rho0                 !< A reference ocean density [kg/m3]
  real    :: Hurr_cen_Y0          !< The initial y position of the hurricane
                                  !!  This experiment is conducted in a Cartesian
                                  !!  grid and this is assumed to be in meters [m]
  real    :: Hurr_cen_X0          !< The initial x position of the hurricane
                                  !!  This experiment is conducted in a Cartesian
                                  !!  grid and this is assumed to be in meters [m]
  real    :: Holland_A            !< Parameter 'A' from the Holland formula
  real    :: Holland_B            !< Parameter 'B' from the Holland formula
  real    :: Holland_AxBxDP       !< 'A' x 'B' x (Pressure Ambient-Pressure central)
                                  !!  for the Holland prorfile calculation
  logical :: relative_tau         !< A logical to take difference between wind
                                  !!  and surface currents to compute the stress


  ! Parameters used if in SCM (single column model) mode
  logical :: SCM_mode        !< Single Column Model Mode [nd]
  logical :: BR_BENCH        !< A "benchmark" configuration (which is meant to
                             !!  provide identical wind to reproduce a previous
                             !!  experiment, where that wind formula contained
                             !!  an error)
  real    :: DY_from_center  !< (Fixed) distance in y from storm center path [m]

  ! Par
  real :: PI
  real :: Deg2Rad

end type

! This include declares and sets the variable "version".
#include "version_variable.h"

character(len=40)  :: mdl = "idealized_hurricane" !< This module's name.

contains

!> Initializes temperature and salinity for the idealized hurricane example
subroutine idealized_hurricane_TS_init(T, S, h, G, GV, param_file, just_read_params)
  type(ocean_grid_type), &
       intent(in)  :: G                !< Grid structure
  type(verticalGrid_type), &
       intent(in)  :: GV               !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), &
       intent(out) :: T                !< Potential temperature (degC)
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), &
       intent(out) :: S                !< Salinity (psu)
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), &
       intent(in)  :: h                !< Layer thickness in H (m or Pa)
  type(param_file_type), &
       intent(in)  :: param_file       !< Input parameter structure
  logical, optional, &
       intent(in)  :: just_read_params !< If present and true, this call will
                                       !! only read parameters without changing h.
  ! Local variables
  real :: top ! The 1-d nominal positions of the upper interface.
  real :: bot ! The 1-d nominal positions of the lower interface.
  real :: S_ref, SST_ref, dTdZ, MLD
  real :: zC
  logical :: just_read    ! If true, just read parameters but set nothing.
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz
  real :: Tbot

  Tbot = 4.0

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) call log_version(param_file, mdl, version)
  call get_param(param_file, mdl,"SALT_REF",S_ref, &
                 'Reference salinity', units='1e-3',default=35.0, &
                 do_not_log=just_read)
  call get_param(param_file, mdl,"TEMP_REF",SST_ref, &
                 'Reference surface temperature', units='C', &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl,"INTERIOR_DTDZ",dTdZ, &
                 'Initial temperature stratification below mixed layer', &
                 units='C/m', fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl,"REF_LAYER_DEPTH",MLD, &
                 'Initial mixed layer depth', units='m', &
                 fail_if_missing=.not.just_read, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  T(:,:,:) = 0.0
  S(:,:,:) = 0.0

  do j=jsd,jed
    do i=isd,ied
      top = 0.
      bot = 0.
      do k=1,nz
        ! Compute next interface
        bot = bot - h(i,j,k)*GV%H_to_m
        ! Depth of middle of layer
        zC = 0.5*( top + bot )
        ! Compute Temperature and Salinity based on decay rates
        T(i,j,k) = max(Tbot,SST_ref + dTdz * min(0., zC + MLD))
        S(i,j,k) = S_ref
        top = bot
      enddo ! k
    enddo
  enddo

end subroutine idealized_hurricane_TS_init

!> Initializes wind profile for the SCM idealized hurricane example
subroutine idealized_hurricane_wind_init(Time, G, param_file, CS)
  type(time_type), &
       intent(in) :: Time       !< Model time
  type(ocean_grid_type), &
       intent(in) :: G          !< Grid structure
  type(param_file_type), &
       intent(in) :: param_file !< Input parameter structure
  type(idealized_hurricane_CS), &
       pointer    :: CS         !< Parameter container

  real :: DP, C

! This include declares and sets the variable "version".
#include "version_variable.h"

  if (associated(CS)) then
    call MOM_error(FATAL, "idealized_hurricane_wind_init called "// &
                          "with an associated control structure.")
    return
  endif

  allocate(CS)

  CS%pi = 4.0*atan(1.0)
  CS%Deg2Rad = CS%pi/180.

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  ! Parameters for computing a wind profile
  call get_param(param_file, mdl, "IDL_HURR_RHO_AIR", CS%rho_a, &
                 "Air density used to compute the idealized hurricane"// &
                 "wind profile.", units='kg/m3', default=1.2)
  call get_param(param_file, mdl, "IDL_HURR_AMBIENT_PRESSURE", &
                 CS%pressure_ambient, "Ambient pressure used in the "// &
                 "idealized hurricane wind profile.", units='Pa', &
                 default=101200.)
  call get_param(param_file, mdl, "IDL_HURR_CENTRAL_PRESSURE", &
                 CS%pressure_central, "Central pressure used in the "// &
                 "idealized hurricane wind profile.", units='Pa', &
                 default=96800.)
  call get_param(param_file, mdl, "IDL_HURR_RAD_MAX_WIND", &
                 CS%rad_max_wind, "Radius of maximum winds used in the"// &
                 "idealized hurricane wind profile.", units='m', &
                 default=50.e3)
  call get_param(param_file, mdl, "IDL_HURR_MAX_WIND", CS%max_windspeed, &
                 "Maximum wind speed used in the idealized hurricane"// &
                 "wind profile.", units='m/s', default=65.)
  call get_param(param_file, mdl, "IDL_HURR_TRAN_SPEED", CS%hurr_translation_spd, &
                 "Translation speed of hurricane used in the idealized"// &
                 "hurricane wind profile.", units='m/s', default=5.0)
  call get_param(param_file, mdl, "IDL_HURR_TRAN_DIR", CS%hurr_translation_dir, &
                 "Translation direction (towards) of hurricane used in the "//&
                 "idealized hurricane wind profile.", units='degrees', &
                 default=180.0)
  CS%hurr_translation_dir = CS%hurr_translation_dir * CS%Deg2Rad
  call get_param(param_file, mdl, "IDL_HURR_X0", CS%Hurr_cen_X0, &
                 "Idealized Hurricane initial X position", &
                 units='m', default=0.)
  call get_param(param_file, mdl, "IDL_HURR_Y0", CS%Hurr_cen_Y0, &
                 "Idealized Hurricane initial Y position", &
                 units='m', default=0.)
  call get_param(param_file, mdl, "IDL_HURR_TAU_CURR_REL", CS%relative_tau, &
                 "Current relative stress switch"//                         &
                 "used in the idealized hurricane wind profile.", &
                 units='', default=.false.)

  ! Parameters for SCM mode
  call get_param(param_file, mdl, "IDL_HURR_SCM_BR_BENCH", CS%BR_BENCH, &
                 "Single column mode benchmark case switch, which is "// &
                 "invoking a modification (bug) in the wind profile meant to "//&
                 "reproduce a previous implementation.", units='', default=.false.)
  call get_param(param_file, mdl, "IDL_HURR_SCM", CS%SCM_MODE, &
                 "Single Column mode switch"//                         &
                 "used in the SCM idealized hurricane wind profile.", &
                 units='', default=.false.)
  call get_param(param_file, mdl, "IDL_HURR_SCM_LOCY", CS%DY_from_center, &
                 "Y distance of station used in the SCM idealized hurricane "// &
                 "wind profile.", units='m', default=50.e3)

  ! The following parameters are model run-time parameters which are used
  ! and logged elsewhere and so should not be logged here. The default
  ! value should be consistent with the rest of the model.
  call get_param(param_file, mdl, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to \n"//&
                 "calculate accelerations and the mass for conservation \n"//&
                 "properties, or with BOUSSINSEQ false to convert some \n"//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, do_not_log=.true.)
  call get_param(param_file, mdl, "GUST_CONST", CS%gustiness, &
                 "The background gustiness in the winds.", units="Pa", &
                 default=0.00, do_not_log=.true.)


  if (CS%BR_BENCH) then
    CS%rho_a = 1.2
  endif
  DP = CS%pressure_ambient - CS%pressure_central
  C = CS%max_windspeed / sqrt( DP )
  CS%Holland_B = C**2 * CS%rho_a * exp(1.0)
  CS%Holland_A = (CS%rad_max_wind)**CS%Holland_B
  CS%Holland_AxBxDP = CS%Holland_A*CS%Holland_B*DP

  return
end subroutine idealized_hurricane_wind_init

!> Computes the surface wind for the idealized hurricane test cases
subroutine idealized_hurricane_wind_forcing(state, forces, day, G, CS)
  type(surface), &
       intent(in)    :: state  !< Surface state structure
  type(mech_forcing), &
       intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type), &
       intent(in)    :: day    !< Time in days
  type(ocean_grid_type), &
       intent(inout) :: G      !< Grid structure
  type(idealized_hurricane_CS), &
       pointer       :: CS     !< Container for idealized hurricane parameters

  ! Local variables
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  real :: TX,TY       !< wind stress
  real :: Uocn, Vocn  !< Surface ocean velocity components
  real :: LAT, LON    !< Grid location
  real :: YY, XX      !< storm relative position
  real :: XC, YC      !< Storm center location
  real :: f           !< Coriolis
  real :: fbench      !< The benchmark 'f' value
  real :: fbench_fac  !< A factor that is set to 0 to use the
                      !!  benchmark 'f' value
  real :: rel_tau_fac !< A factor that is set to 0 to disable
                      !!  current relative stress calculation

  ! Bounds for loops and memory allocation
  is = G%isc    ; ie = G%iec    ; js = G%jsc    ; je = G%jec
  Isq = G%IscB  ; Ieq = G%IecB  ; Jsq = G%JscB  ; Jeq = G%JecB
  isd = G%isd   ; ied = G%ied   ; jsd = G%jsd   ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ! Allocate the forcing arrays, if necessary.
  call allocate_mech_forcing(G, forces, stress=.true., ustar=.true.)

  if (CS%relative_tau) then
     REL_TAU_FAC = 1.
  else
     REL_TAU_FAC = 0. !Multiplied to 0 surface current
  endif

  !> Compute storm center location
  XC = CS%Hurr_cen_X0 + (time_type_to_real(day)*CS%hurr_translation_spd*&
       cos(CS%hurr_translation_dir))
  YC = CS%Hurr_cen_Y0 + (time_type_to_real(day)*CS%hurr_translation_spd*&
       sin(CS%hurr_translation_dir))

  if (CS%BR_Bench) then
     ! f reset to value used in generated wind for benchmark test
     fbench = 5.5659e-05
     fbench_fac = 0.0
  else
     fbench = 0.0
     fbench_fac = 1.0
  endif

  !> Computes taux
  do j=js,je
    do I=is-1,Ieq
      Uocn = state%u(I,j)*REL_TAU_FAC
      Vocn = 0.25*(state%v(i,J)+state%v(i+1,J-1)&
             +state%v(i+1,J)+state%v(i,J-1))*REL_TAU_FAC
      f = abs(0.5*(G%CoriolisBu(I,J)+G%CoriolisBu(I,J-1)))*fbench_fac &
           + fbench
      ! Calculate position as a function of time.
      if (CS%SCM_mode) then
        YY = YC
        XX = XC
      else
        LAT = G%geoLatCu(I,j)*1000. !KM_to_m
        LON = G%geoLonCu(I,j)*1000. !KM_to_m
        YY = LAT - YC
        XX = LON - XC
      endif
      call idealized_hurricane_wind_profile(&
         CS,f,YY,XX,Uocn,Vocn,TX,TY)
      forces%taux(I,j) = G%mask2dCu(I,j) * TX
    enddo
  enddo
  !> Computes tauy
  do J=js-1,Jeq
    do i=is,ie
      Uocn = 0.25*(state%u(I,j)+state%u(I-1,j+1)&
            +state%u(I-1,j)+state%u(I,j+1))*REL_TAU_FAC
      Vocn = state%v(i,J)*REL_TAU_FAC
      f = abs(0.5*(G%CoriolisBu(I-1,J)+G%CoriolisBu(I,J)))*fbench_fac &
           + fbench
      ! Calculate position as a function of time.
      if (CS%SCM_mode) then
        YY = YC
        XX = XC
      else
        LAT = G%geoLatCv(i,J)*1000. !KM_to_m
        LON = G%geoLonCv(i,J)*1000. !KM_to_m
        YY = LAT - YC
        XX = LON - XC
      endif
      call idealized_hurricane_wind_profile(&
           CS,f,YY,XX,Uocn,Vocn,TX,TY)
      forces%tauy(i,J) = G%mask2dCv(i,J) * TY
    enddo
  enddo

  !> Get Ustar
  do j=js,je
    do i=is,ie
      !  This expression can be changed if desired, but need not be.
      forces%ustar(i,j) = G%mask2dT(i,j) * sqrt(CS%gustiness/CS%Rho0 + &
         sqrt(0.5*(forces%taux(I-1,j)**2 + forces%taux(I,j)**2) + &
            0.5*(forces%tauy(i,J-1)**2 + forces%tauy(i,J)**2))/CS%Rho0)
    enddo
  enddo

  return
end subroutine idealized_hurricane_wind_forcing

!> Calculate the wind speed at a location as a function of time.
subroutine idealized_hurricane_wind_profile(CS,absf,YY,XX,UOCN,VOCN,Tx,Ty)
  ! Author: Brandon Reichl
  ! Date: Nov-20-2014
  !       Aug-14-2018 Generalized for non-SCM configuration

  ! Input parameters
  type(idealized_hurricane_CS), &
        pointer     :: CS   !< Container for SCM parameters
  real, intent(in)  :: absf !<Input coriolis magnitude
  real, intent(in)  :: YY   !< Location in m relative to center y
  real, intent(in)  :: XX   !< Location in m relative to center x
  real, intent(in)  :: UOCN !< X surface current
  real, intent(in)  :: VOCN !< Y surface current
  real, intent(out) :: Tx   !< X stress
  real, intent(out) :: Ty   !< Y stress

  ! Local variables

  ! Wind profile terms
  real :: U10
  real :: radius
  real :: radius10
  real :: radius_km
  real :: radiusB
  real :: fcor
  real :: du10
  real :: du
  real :: dv
  real :: CD

  !Wind angle variables
  real :: Alph !< The resulting inflow angle (positive outward)
  real :: Rstr
  real :: A0
  real :: A1
  real :: P1
  real :: Adir
  real :: V_TS
  real :: U_TS

  ! Implementing Holland (1980) parameteric wind profile

  Radius = SQRT(XX**2.+YY**2.)

  !/ BGR
  ! rkm - r converted to km for Holland prof.
  !       used in km due to error, correct implementation should
  !       not need rkm, but to match winds w/ experiment this must
  !       be maintained.  Causes winds far from storm center to be a
  !       couple of m/s higher than the correct Holland prof.
  if (CS%BR_Bench) then
    radius_km = radius/1000.
  else
    ! if not comparing to benchmark, then use correct Holland prof.
    radius_km = radius
  endif
  radiusB = (radius)**CS%Holland_B

  !/
  ! Calculate U10 in the interior (inside of 10x radius of maximum wind),
  ! while adjusting U10 to 0 outside of 12x radius of maximum wind.
  if ( (radius/CS%rad_max_wind .gt. 0.001) .and. &
       (radius/CS%rad_max_wind .lt. 10.) ) then
    U10 = sqrt(CS%Holland_AxBxDP*exp(-CS%Holland_A/radiusB)/(CS%rho_A*radiusB)&
                +0.25*(radius_km*absf)**2) - 0.5*radius_km*absf
  elseif ( (radius/CS%rad_max_wind .gt. 10.) .and. &
           (radius/CS%rad_max_wind .lt. 15.) ) then

    radius10 = CS%rad_max_wind*10.
    if (CS%BR_Bench) then
      radius_km = radius10/1000.
    else
      radius_km = radius10
    endif
    radiusB=radius10**CS%Holland_B

    U10 = (sqrt(CS%Holland_AxBxDp*exp(-CS%Holland_A/radiusB)/(CS%rho_A*radiusB)&
                  +0.25*(radius_km*absf)**2)-0.5*radius_km*absf) &
           * (15.-radius/CS%rad_max_wind)/5.
  else
    U10 = 0.
  endif
  Adir = atan2(YY,xx)
  !\

  ! Wind angle model following Zhang and Ulhorn (2012)
  ! ALPH is inflow angle positive outward.
  RSTR = min(10.,radius / CS%rad_max_wind)
  A0 = -0.9*RSTR - 0.09*CS%max_windspeed - 14.33
  A1 = -A0*(0.04*RSTR + 0.05*CS%Hurr_translation_spd + 0.14)
  P1 = (6.88*RSTR - 9.60*CS%Hurr_translation_spd + 85.31) * CS%Deg2Rad
  ALPH = A0 - A1*cos(CS%hurr_translation_dir-Adir-P1)
  if ( (radius/CS%rad_max_wind.gt.10.) .and.&
       (radius/CS%rad_max_wind.lt.15.) ) then
     ALPH = ALPH*(15.0-radius/CS%rad_max_wind)/5.
  elseif (radius/CS%rad_max_wind.gt.15.) then
     ALPH = 0.0
  endif
  ALPH = ALPH * CS%Deg2Rad

  ! Calculate translation speed components
  U_TS = CS%hurr_translation_spd/2.*cos(CS%hurr_translation_dir)
  V_TS = CS%hurr_translation_spd/2.*sin(CS%hurr_translation_dir)

  ! Set output (relative) winds
  dU = U10*sin(Adir-CS%Pi-Alph) - UOCN + U_TS
  dV = U10*cos(Adir-Alph) - VOCN + V_TS

  !  Use a simple drag coefficient as a function of U10 (from Sullivan et al., 2010)
  du10 = sqrt(du**2+dv**2)
  if (du10.lt.11.) then
     Cd = 1.2e-3
  elseif (du10.lt.20.0) then
     Cd = (0.49 + 0.065*U10)*1.e-3
  else
     Cd = 1.8e-3
  endif

  ! Compute stress vector
  TX = CS%rho_A * Cd * sqrt(du**2+dV**2) * dU
  TY = CS%rho_A * Cd * sqrt(du**2+dV**2) * dV

  return
end subroutine idealized_hurricane_wind_profile

end module idealized_hurricane
