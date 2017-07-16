!> Initial conditions and forcing for the single column model (SCM) idealized
!! hurricane example.
module SCM_idealized_hurricane

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing, allocate_forcing_type
use MOM_grid, only : ocean_grid_type
use MOM_safe_alloc, only : safe_alloc_ptr
use MOM_time_manager, only : time_type, operator(+), operator(/), get_time,&
                             time_type_to_real
use MOM_variables, only : thermo_var_ptrs, surface
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public SCM_idealized_hurricane_TS_init
public SCM_idealized_hurricane_wind_init
public SCM_idealized_hurricane_wind_forcing
public SCM_idealized_hurricane_CS

!> Container for parameters describing idealized wind structure
type SCM_idealized_hurricane_CS ; private
  real :: rho_a  !< Air density
  real :: p_n    !< Ambient pressure
  real :: p_c    !< Central pressure
  real :: r_max  !< Radius of maximum winds
  real :: U_max  !< Maximum wind speeds
  real :: YY     !< Distance (positive north) of storm center
  real :: tran_speed !< Hurricane translation speed
  real :: gust_const !< Gustiness (used in u*)
  real :: Rho0   !< A reference ocean density in kg/m3
end type

! This include declares and sets the variable "version".
#include "version_variable.h"

character(len=40)  :: mod = "SCM_idealized_hurricane" ! This module's name.

contains

!> Initializes temperature and salinity for the SCM idealized hurricane example
subroutine SCM_idealized_hurricane_TS_init(T, S, h, G, GV, param_file)
  type(ocean_grid_type),                     intent(in)  :: G !< Grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: T !< Potential temperature (degC)
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: S !< Salinity (psu)
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(in)  :: h !< Layer thickness (m or Pa)
  type(param_file_type),                     intent(in)  :: param_file !< Input parameter structure
  ! Local variables
  real :: eta(SZK_(G)+1) ! The 1-d nominal positions of the interfaces.
  real :: S_ref, SST_ref, dTdZ, MLD
  real :: zC
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call log_version(param_file, mod, version)
  call get_param(param_file,mod,"SCM_S_REF",S_ref, &
                 'Reference salinity', units='1e-3',default=35.0)
  call get_param(param_file,mod,"SCM_SST_REF",SST_ref, &
                 'Reference surface temperature', units='C', fail_if_missing=.true.)
  call get_param(param_file,mod,"SCM_DTDZ",dTdZ,                         &
                 'Initial temperature stratification below mixed layer', &
                 units='C/m', fail_if_missing=.true.)
  call get_param(param_file,mod,"SCM_MLD",MLD, &
                 'Initial mixed layer depth', units='m', fail_if_missing=.true.)

  do j=js,je ; do i=is,ie
    eta(1) = 0. ! Reference to surface
    do k=1,nz
      eta(K+1) = eta(K) - h(i,j,k)*GV%H_to_m ! Interface below layer (in m)
      zC = 0.5*( eta(K) + eta(K+1) )        ! Z of middle of layer (in m)
      T(i,j,k) = SST_ref + dTdz*min(0., zC+MLD)
      S(i,j,k) = S_ref
    enddo ! k
  enddo ; enddo

end subroutine SCM_idealized_hurricane_TS_init

!> Initializes wind profile for the SCM idealized hurricane example
subroutine SCM_idealized_hurricane_wind_init(Time, G, param_file, CS)
  type(time_type),              intent(in) :: Time       !< Time
  type(ocean_grid_type),        intent(in) :: G          !< Grid structure
  type(param_file_type),        intent(in) :: param_file !< Input parameter structure
  type(SCM_idealized_hurricane_CS), pointer    :: CS         !< Parameter container

! This include declares and sets the variable "version".
#include "version_variable.h"

  if (associated(CS)) then
    call MOM_error(FATAL, "SCM_idealized_hurricane_wind_init called with an associated "// &
                          "control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "SCM_RHO_AIR", CS%rho_a,            &
                 "Air density "//                                     &
                 "used in the SCM idealized hurricane wind profile.", &
                 units='kg/m3', default=1.2)
  call get_param(param_file, mod, "SCM_AMBIENT_PRESSURE", CS%p_n,     &
                 "Ambient pressure "//                                &
                 "used in the SCM idealized hurricane wind profile.", &
                 units='Pa', default=101200.)
  call get_param(param_file, mod, "SCM_CENTRAL_PRESSURE", CS%p_c,     &
                 "Central pressure "//                                &
                 "used in the SCM idealized hurricane wind profile.", &
                 units='Pa', default=96800.)
  call get_param(param_file, mod, "SCM_RADIUS_MAX_WINDS", CS%r_max,   &
                 "Radius of maximum winds "//                         &
                 "used in the SCM idealized hurricane wind profile.", &
                 units='m', default=50.e3)
  call get_param(param_file, mod, "SCM_MAX_WIND_SPEED", CS%U_max,     &
                 "Maximum wind speed "//                              &
                 "used in the SCM idealized hurricane wind profile.", &
                 units='m/s', default=65.)
  call get_param(param_file, mod, "SCM_YY", CS%YY,     &
                 "Y distance of station "//                           &
                 "used in the SCM idealized hurricane wind profile.", &
                 units='m', default=50.e3)
  call get_param(param_file, mod, "SCM_TRAN_SPEED", CS%TRAN_SPEED,     &
                 "Translation speed of hurricane"//                   &
                 "used in the SCM idealized hurricane wind profile.", &
                 units='m/s', default=5.0)
  call get_param(param_file, mod, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to \n"//&
                 "calculate accelerations and the mass for conservation \n"//&
                 "properties, or with BOUSSINSEQ false to convert some \n"//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0)
  ! The following parameter is a model run-time parameter which is used
  ! and logged elsewhere and so should not be logged here. The default
  ! value should be consistent with the rest of the model.
  call get_param(param_file, mod, "GUST_CONST", CS%gust_const, &
                 "The background gustiness in the winds.", units="Pa", &
                 default=0.00, do_not_log=.true.)


end subroutine SCM_idealized_hurricane_wind_init

subroutine SCM_idealized_hurricane_wind_forcing(state, fluxes, day, G, CS)
  type(surface),                    intent(in)    :: state  !< Surface state structure
  type(forcing),                    intent(inout) :: fluxes !< Surface fluxes structure
  type(time_type),                  intent(in)    :: day    !< Time in days
  type(ocean_grid_type),            intent(inout) :: G      !< Grid structure
  type(SCM_idealized_hurricane_CS), pointer       :: CS     !< Container for SCM parameters
  ! Local variables
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  real :: pie, Deg2Rad
  real :: U10, A, B, C, r, f,du10,rkm ! For wind profile expression
  real :: xx, t0 !for location
  real :: dp, rB
  real :: Cd ! Air-sea drag coefficient
  real :: Uocn, Vocn ! Surface ocean velocity components
  real :: dU, dV ! Air-sea differential motion
  !Wind angle variables
  real :: Alph,Rstr, A0, A1, P1, Adir, transdir, V_TS, U_TS
  logical :: BR_Bench
  ! Bounds for loops and memory allocation
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ! Allocate the forcing arrays, if necessary.
  call allocate_forcing_type(G, fluxes, stress=.true., ustar=.true.)

  pie = 4.0*atan(1.0) ; Deg2Rad = pie/180.

  !/ BR
  ! Implementing Holland (1980) parameteric wind profile
  !------------------------------------------------------|
  BR_Bench = .true.   !true if comparing to LES runs     |
  t0 = 129600.        !TC 'eye' crosses (0,0) at 36 hours|
  transdir = pie      !translation direction (-x)        |
  !------------------------------------------------------|
  dp = CS%p_n - CS%p_c
  C = CS%U_max / sqrt( dp )
  B = C**2 * CS%rho_a * exp(1.0)
  if (BR_Bench) then
     ! rho_a reset to value used in generated wind for benchmark test
     B = C**2 * 1.2 * exp(1.0)
  endif
  A = (CS%r_max/1000.)**B
  f =G%CoriolisBu(is,js) ! f=f(x,y) but in the SCM is constant
  if (BR_Bench) then
     ! f reset to value used in generated wind for benchmark test
     f = 5.5659e-05
  endif
  !/ BR
  ! Calculate x position as a function of time.
  xx = ( t0 - time_type_to_real(day)) * CS%tran_speed * cos(transdir)
  r = sqrt(xx**2.+CS%YY**2.)
  !/ BR
  ! rkm - r converted to km for Holland prof.
  !       used in km due to error, correct implementation should
  !       not need rkm, but to match winds w/ experiment this must
  !       be maintained.  Causes winds far from storm center to be a
  !       couple of m/s higher than the correct Holland prof.
  if (BR_Bench) then
     rkm = r/1000.
     rB = (rkm)**B
  else
     ! if not comparing to benchmark, then use correct Holland prof.
     rkm = r
     rB = r**B
  endif
  !/ BR
  ! Calculate U10 in the interior (inside of 10x radius of maximum wind),
  ! while adjusting U10 to 0 outside of 12x radius of maximum wind.
  ! Note that rho_a is set to 1.2 following generated wind for experiment
  if (r/CS%r_max.gt.0.001 .AND. r/CS%r_max.lt.10.) then
     U10 = sqrt( A*B*dp*exp(-A/rB)/(1.2*rB) + 0.25*(rkm*f)**2 ) - 0.5*rkm*f
  elseif (r/CS%r_max.gt.10. .AND. r/CS%r_max.lt.12.) then
     r=CS%r_max*10.
     if (BR_Bench) then
        rkm = r/1000.
        rB=rkm**B
     else
        rkm = r
        rB = r**B
     endif
     U10 = ( sqrt( A*B*dp*exp(-A/rB)/(1.2*rB) + 0.25*(rkm*f)**2 ) - 0.5*rkm*f) &
           * (12. - r/CS%r_max)/2.
  else
     U10 = 0.
  end if
  Adir = atan2(CS%YY,xx)

  !/ BR
  ! Wind angle model following Zhang and Ulhorn (2012)
  ! ALPH is inflow angle positive outward.
  RSTR = min(10.,r / CS%r_max)
  A0 = -0.9*RSTR -0.09*CS%U_max -14.33
  A1 = -A0 *(0.04*RSTR +0.05*CS%tran_speed+0.14)
  P1 = (6.88*RSTR -9.60*CS%tran_speed+85.31)*pie/180.
  ALPH = A0 - A1*cos( (TRANSDIR - ADIR ) - P1)
  if (r/CS%r_max.gt.10. .AND. r/CS%r_max.lt.12.) then
     ALPH = ALPH* (12. - r/CS%r_max)/2.
  elseif (r/CS%r_max.gt.12.) then
     ALPH = 0.0
  endif
  ALPH = ALPH * Deg2Rad

  !/BR
  ! Prepare for wind calculation
  ! X_TS is component of translation speed added to wind vector
  ! due to background steering wind.
  U_TS = CS%tran_speed/2.*cos(transdir)
  V_TS = CS%tran_speed/2.*sin(transdir)

  ! Set the surface wind stresses, in units of Pa. A positive taux
  ! accelerates the ocean to the (pseudo-)east.
  !   The i-loop extends to is-1 so that taux can be used later in the
  ! calculation of ustar - otherwise the lower bound would be Isq.
  do j=js,je ; do I=is-1,Ieq
    !/BR
    ! Turn off surface current for stress calculation to be
    ! consistent with test case.
    Uocn = 0.!state%u(I,j)
    Vocn = 0.!0.25*( (state%v(i,J) + state%v(i+1,J-1)) &
             !    +(state%v(i+1,J) + state%v(i,J-1)) )
    !/BR
    ! Wind vector calculated from location/direction (sin/cos flipped b/c
    ! cyclonic wind is 90 deg. phase shifted from position angle).
    dU = U10*sin(Adir-pie-Alph) - Uocn + U_TS
    dV = U10*cos(Adir-Alph) - Vocn + V_TS
    !/----------------------------------------------------|
    !BR
    !  Add a simple drag coefficient as a function of U10 |
    !/----------------------------------------------------|
    du10=sqrt(du**2+dv**2)
    if (du10.LT.11.) then
       Cd = 1.2e-3
    elseif (du10.LT.20.) then
       Cd = (0.49 + 0.065 * U10 )*0.001
    else
       Cd = 0.0018
    endif
    fluxes%taux(I,j) = CS%rho_a * G%mask2dCu(I,j) * Cd*sqrt(du**2+dV**2)*dU
  enddo ; enddo
  !/BR
  ! See notes above
  do J=js-1,Jeq ; do i=is,ie
    Uocn = 0.!0.25*( (state%u(I,j) + state%u(I-1,j+1)) &
             !    +(state%u(I-1,j) + state%u(I,j+1)) )
    Vocn = 0.!state%v(i,J)
    dU = U10*sin(Adir-pie-Alph) - Uocn + U_TS
    dV = U10*cos(Adir-Alph) - Vocn + V_TS
    du10=sqrt(du**2+dv**2)
    if (du10.LT.11.) then
       Cd = 1.2e-3
    elseif (du10.LT.20.) then
       Cd = (0.49 + 0.065 * U10 )*0.001
    else
       Cd = 0.0018
    endif
    fluxes%tauy(I,j) = CS%rho_a * G%mask2dCv(I,j) * Cd*du10*dV
  enddo ; enddo
  ! Set the surface friction velocity, in units of m s-1. ustar is always positive.
  do j=js,je ; do i=is,ie
    !  This expression can be changed if desired, but need not be.
    fluxes%ustar(i,j) = G%mask2dT(i,j) * sqrt(CS%gust_const/CS%Rho0 + &
       sqrt(0.5*(fluxes%taux(I-1,j)**2 + fluxes%taux(I,j)**2) + &
            0.5*(fluxes%tauy(i,J-1)**2 + fluxes%tauy(i,J)**2))/CS%Rho0)
  enddo ; enddo

end subroutine SCM_idealized_hurricane_wind_forcing

end module SCM_idealized_hurricane
