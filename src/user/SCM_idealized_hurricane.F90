!> Initial conditions and forcing for the single column model (SCM) idealized
!! hurricane example.
module SCM_idealized_hurricane

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_safe_alloc, only : safe_alloc_ptr
use MOM_time_manager, only : time_type, operator(+), operator(/), get_time
use MOM_variables, only : thermo_var_ptrs, surface
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
  real :: gust_const !< Gustiness (used in u*)
end type

! This include declares and sets the variable "version".
#include "version_variable.h"

character(len=40)  :: mod = "SCM_idealized_hurricane" ! This module's name.

contains

!> Initializes temperature and salinity for the SCM idealized hurricane example
subroutine SCM_idealized_hurricane_TS_init(T, S, h, G, param_file)
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(out) :: T !< Potential temperature (degC)
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(out) :: S !< Salinity (psu)
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(in)  :: h !< Layer thickness (m or Pa)
  type(ocean_grid_type),                  intent(in)  :: G !< Grid structure
  type(param_file_type),                  intent(in)  :: param_file !< Input parameter structure
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
      eta(K+1) = eta(K) - h(i,j,k)*G%H_to_m ! Interface below layer (in m)
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
                 units='kg/m3', default=1.25)
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
  ! The following parameter is a model run-time parameter which is used
  ! and logged elsewhere and so should not be logged here. The default
  ! value should be consistent with the rest of the model.
  call get_param(param_file, mod, "GUST_CONST", CS%gust_const, &
                 "The background gustiness in the winds.", units="Pa", &
                 default=0.02, do_not_log=.true.)


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
  real :: U10, A, B, C, r, f ! For wind profile expression
  real :: dp, rB
  real :: Cd ! Air-sea drag coefficient
  real :: Uocn, Vocn ! Surface ocean velocity components
  real :: dU, dV ! Air-sea differential motion

  ! Bounds for loops and memory allocation
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ! Allocate the forcing arrays, if necessary.
  call safe_alloc_ptr(fluxes%taux, IsdB, IedB, jsd, jed)
  call safe_alloc_ptr(fluxes%tauy, isd, ied, JsdB, JedB)
  call safe_alloc_ptr(fluxes%ustar, isd, ied, jsd, jed)

  dp = CS%p_n - CS%p_c
  C = CS%U_max / sqrt( dp )
  B = C**2 * CS%rho_a * exp(1.0)
  A = CS%r_max**B
  f = G%CoriolisBu(is,js) ! f=f(x,y) but in the SCM is constant

  r = 70.e3 ! Should be function of time

  rB = r**B
  U10 = sqrt( A*B*dp*exp(-A/rB)/(CS%rho_a*rB) + 0.25*(r*f)**2 ) - 0.5*r*f

  Cd = 1.e-3 ! Probably should be a function of sea-state, U10, etc.

  ! Set the surface wind stresses, in units of Pa. A positive taux
  ! accelerates the ocean to the (pseudo-)east.
  !   The i-loop extends to is-1 so that taux can be used later in the
  ! calculation of ustar - otherwise the lower bound would be Isq.
  do j=js,je ; do I=is-1,Ieq
    Uocn = state%u(I,j)
    Vocn = 0.25*( (state%v(i,J) + state%v(i+1,J-1)) &
                 +(state%v(i+1,J) + state%v(i,J-1)) )
    dU = U10 - Uocn
    dV = 0. - Uocn
    fluxes%taux(I,j) = G%mask2dCu(I,j) * Cd*abs(du**2+dV**2)*dU
  enddo ; enddo
  do J=js-1,Jeq ; do i=is,ie
    Uocn = 0.25*( (state%u(I,j) + state%u(I-1,j+1)) &
                 +(state%u(I-1,j) + state%u(I,j+1)) )
    Vocn = state%v(i,J)
    dU = U10 - Uocn
    dV = 0. - Uocn
    fluxes%tauy(i,J) = G%mask2dCv(i,J) * Cd*abs(du**2+dV**2)*dV
  enddo ; enddo

  ! Set the surface friction velocity, in units of m s-1. ustar is always positive.
  do j=js,je ; do i=is,ie
    !  This expression can be changed if desired, but need not be.
    fluxes%ustar(i,j) = G%mask2dT(i,j) * sqrt(CS%gust_const/G%Rho0 + &
       sqrt(0.5*(fluxes%taux(I-1,j)**2 + fluxes%taux(I,j)**2) + &
            0.5*(fluxes%tauy(i,J-1)**2 + fluxes%tauy(i,J)**2))/G%Rho0)
  enddo ; enddo

end subroutine SCM_idealized_hurricane_wind_forcing

end module SCM_idealized_hurricane
