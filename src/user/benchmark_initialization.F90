!> Initialization for the "bench mark" configuration
module benchmark_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type

implicit none ; private

#include <MOM_memory.h>

public benchmark_initialize_topography
public benchmark_initialize_thickness
public benchmark_init_temperature_salinity

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> This subroutine sets up the benchmark test case topography.
subroutine benchmark_initialize_topography(D, G, param_file, max_depth, US)
  type(dyn_horgrid_type),          intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                   intent(out) :: D !< Ocean bottom depth in m or [Z ~> m] if US is present
  type(param_file_type),           intent(in)  :: param_file !< Parameter file structure
  real,                            intent(in)  :: max_depth !< Maximum model depth in the units of D
  type(unit_scale_type), optional, intent(in)  :: US !< A dimensional unit scaling type

  ! Local variables
  real :: min_depth            ! The minimum and maximum depths [Z ~> m].
  real :: PI                   ! 3.1415926... calculated as 4*atan(1)
  real :: D0                   ! A constant to make the maximum     !
                               ! basin depth MAXIMUM_DEPTH.         !
  real :: m_to_Z  ! A dimensional rescaling factor.
  real :: x, y
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "benchmark_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call MOM_mesg("  benchmark_initialization.F90, benchmark_initialize_topography: setting topography", 5)

  m_to_Z = 1.0 ; if (present(US)) m_to_Z = US%m_to_Z

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0, scale=m_to_Z)

  PI = 4.0*atan(1.0)
  D0 = max_depth / 0.5

!  Calculate the depth of the bottom.
  do j=js,je ; do i=is,ie
    x = (G%geoLonT(i,j)-G%west_lon) / G%len_lon
    y = (G%geoLatT(i,j)-G%south_lat) / G%len_lat
!  This sets topography that has a reentrant channel to the south.
    D(i,j) = -D0 * ( y*(1.0 + 0.6*cos(4.0*PI*x)) &
                   + 0.75*exp(-6.0*y) &
                   + 0.05*cos(10.0*PI*x) - 0.7 )
    if (D(i,j) > max_depth) D(i,j) = max_depth
    if (D(i,j) < min_depth) D(i,j) = 0.
  enddo ; enddo

end subroutine benchmark_initialize_topography

!> Initializes layer thicknesses for the benchmark test case,
!! by finding the depths of interfaces in a specified latitude-dependent
!! temperature profile with an exponentially decaying thermocline on top of a
!! linear stratification.
subroutine benchmark_initialize_thickness(h, depth_tot, G, GV, US, param_file, eqn_of_state, &
                                          P_Ref, just_read_params)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  type(EOS_type),          pointer     :: eqn_of_state !< Equation of state structure
  real,                    intent(in)  :: P_Ref       !< The coordinate-density
                                                      !! reference pressure [R L2 T-2 ~> Pa].
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.
  ! Local variables
  real :: e0(SZK_(GV)+1)     ! The resting interface heights, in depth units [Z ~> m],
                             ! usually negative because it is positive upward.
  real :: e_pert(SZK_(GV)+1) ! Interface height perturbations, positive upward,
                             ! in depth units [Z ~> m].
  real :: eta1D(SZK_(GV)+1)  ! Interface height relative to the sea surface
                             ! positive upward, in depth units [Z ~> m].
  real :: SST       !  The initial sea surface temperature [degC].
  real :: T_int     !  The initial temperature of an interface [degC].
  real :: ML_depth  !  The specified initial mixed layer depth, in depth units [Z ~> m].
  real :: thermocline_scale ! The e-folding scale of the thermocline, in depth units [Z ~> m].
  real, dimension(SZK_(GV)) :: &
    T0, S0, &       ! Profiles of temperature [degC] and salinity [ppt]
    rho_guess, &    ! Potential density at T0 & S0 [R ~> kg m-3].
    drho_dT, &      ! Derivative of density with temperature [R degC-1 ~> kg m-3 degC-1].
    drho_dS         ! Derivative of density with salinity [R ppt-1 ~> kg m-3 ppt-1].
  real :: pres(SZK_(GV))  ! Reference pressure [R L2 T-2 ~> Pa].
  real :: a_exp     ! The fraction of the overall stratification that is exponential.
  real :: I_ts, I_md ! Inverse lengthscales [Z-1 ~> m-1].
  real :: T_frac    ! A ratio of the interface temperature to the range
                    ! between SST and the bottom temperature.
  real :: err, derr_dz  ! The error between the profile's temperature and the
                    ! interface temperature for a given z and its derivative.
  real :: pi, z
  logical :: just_read
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "benchmark_initialize_thickness" ! This subroutine's name.
  integer :: i, j, k, k1, is, ie, js, je, nz, itt

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params
  if (.not.just_read) call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "BENCHMARK_ML_DEPTH_IC", ML_depth, &
                 "Initial mixed layer depth in the benchmark test case.", &
                 units='m', default=50.0, scale=US%m_to_Z, do_not_log=just_read)
  call get_param(param_file, mdl, "BENCHMARK_THERMOCLINE_SCALE", thermocline_scale, &
                 "Initial thermocline depth scale in the benchmark test case.", &
                 default=500.0, units="m", scale=US%m_to_Z, do_not_log=just_read)

  if (just_read) return ! This subroutine has no run-time parameters.

  call MOM_mesg("  benchmark_initialization.F90, benchmark_initialize_thickness: setting thickness", 5)

  k1 = GV%nk_rho_varies + 1

  a_exp = 0.9

! This block calculates T0(k) for the purpose of diagnosing where the
! interfaces will be found.
  do k=1,nz
    pres(k) = P_Ref ; S0(k) = 35.0
  enddo
  T0(k1) = 29.0
  call calculate_density(T0(k1), S0(k1), pres(k1), rho_guess(k1), eqn_of_state)
  call calculate_density_derivs(T0(k1), S0(k1), pres(k1), drho_dT(k1), drho_dS(k1), eqn_of_state)

! A first guess of the layers' temperatures.
  do k=1,nz
    T0(k) = T0(k1) + (GV%Rlay(k) - rho_guess(k1)) / drho_dT(k1)
  enddo

! Refine the guesses for each layer.
  do itt=1,6
    call calculate_density(T0, S0, pres, rho_guess, eqn_of_state)
    call calculate_density_derivs(T0, S0, pres, drho_dT, drho_dS, eqn_of_state)
    do k=1,nz
      T0(k) = T0(k) + (GV%Rlay(k) - rho_guess(k)) / drho_dT(k)
    enddo
  enddo

  pi = 4.0*atan(1.0)
  I_ts = 1.0 / thermocline_scale
  I_md = 1.0 / G%max_depth
  do j=js,je ; do i=is,ie
    SST = 0.5*(T0(k1)+T0(nz)) - 0.9*0.5*(T0(k1)-T0(nz)) * &
                               cos(pi*(G%geoLatT(i,j)-G%south_lat)/(G%len_lat))

    do k=1,nz ; e_pert(K) = 0.0 ; enddo

    !   This sets the initial thickness (in [H ~> m or kg m-2]) of the layers.  The thicknesses
    ! are set to insure that: 1. each layer is at least  Gv%Angstrom_m thick, and
    ! 2. the interfaces are where they should be based on the resting depths and interface
    ! height perturbations, as long at this doesn't interfere with 1.
    eta1D(nz+1) = -depth_tot(i,j)

    do k=nz,2,-1
      T_int = 0.5*(T0(k) + T0(k-1))
      T_frac = (T_int - T0(nz)) / (SST - T0(nz))
      ! Find the z such that T_frac = a exp(z/thermocline_scale) + (1-a) (z+D)/D
      z = 0.0
      do itt=1,6
        err = a_exp * exp(z*I_ts) + (1.0 - a_exp) * (z*I_md + 1.0) - T_frac
        derr_dz = a_exp * I_ts * exp(z*I_ts) + (1.0 - a_exp) * I_md
        z = z - err / derr_dz
      enddo
      e0(K) = z
!       e0(K) = -ML_depth + thermocline_scale * log((T_int - T0(nz)) / (SST - T0(nz)))

      eta1D(K) = e0(K) + e_pert(K)

      if (eta1D(K) > -ML_depth) eta1D(K) = -ML_depth

      if (eta1D(K) < eta1D(K+1) + GV%Angstrom_Z) &
        eta1D(K) = eta1D(K+1) + GV%Angstrom_Z

      h(i,j,k) = max(GV%Z_to_H * (eta1D(K) - eta1D(K+1)), GV%Angstrom_H)
    enddo
    h(i,j,1) = max(GV%Z_to_H * (0.0 - eta1D(2)), GV%Angstrom_H)

  enddo ; enddo

end subroutine benchmark_initialize_thickness

!> Initializes layer temperatures and salinities for benchmark
subroutine benchmark_init_temperature_salinity(T, S, G, GV, US, param_file, &
               eqn_of_state, P_Ref, just_read_params)
  type(ocean_grid_type),               intent(in)  :: G            !< The ocean's grid structure.
  type(verticalGrid_type),             intent(in)  :: GV           !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T      !< The potential temperature
                                                                   !! that is being initialized [degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S      !< The salinity that is being
                                                                   !! initialized [ppt]
  type(unit_scale_type),               intent(in)  :: US           !< A dimensional unit scaling type
  type(param_file_type),               intent(in)  :: param_file   !< A structure indicating the
                                                                   !! open file to parse for
                                                                   !! model parameter values.
  type(EOS_type),                      pointer     :: eqn_of_state !< Equation of state structure
  real,                                intent(in)  :: P_Ref        !< The coordinate-density
                                                                   !! reference pressure [R L2 T-2 ~> Pa].
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.
  ! Local variables
  real :: T0(SZK_(GV))       ! A profile of temperatures [degC]
  real :: S0(SZK_(GV))       ! A profile of salinities [ppt]
  real :: pres(SZK_(GV))     ! Reference pressure [R L2 T-2 ~> Pa].
  real :: drho_dT(SZK_(GV))  ! Derivative of density with temperature [R degC-1 ~> kg m-3 degC-1].
  real :: drho_dS(SZK_(GV))  ! Derivative of density with salinity [R ppt-1 ~> kg m-3 ppt-1].
  real :: rho_guess(SZK_(GV)) ! Potential density at T0 & S0 [R ~> kg m-3].
  real :: PI        ! 3.1415926... calculated as 4*atan(1)
  real :: SST       !  The initial sea surface temperature [degC].
  real :: lat
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=40)  :: mdl = "benchmark_init_temperature_salinity" ! This subroutine's name.
  integer :: i, j, k, k1, is, ie, js, je, nz, itt

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (just_read) return ! All run-time parameters have been read, so return.

  k1 = GV%nk_rho_varies + 1

  do k=1,nz
    pres(k) = P_Ref ; S0(k) = 35.0
  enddo

  T0(k1) = 29.0
  call calculate_density(T0(k1), S0(k1), pres(k1), rho_guess(k1), eqn_of_state)
  call calculate_density_derivs(T0, S0, pres, drho_dT, drho_dS, eqn_of_state, (/k1,k1/) )

! A first guess of the layers' temperatures.                         !
  do k=1,nz
    T0(k) = T0(k1) + (GV%Rlay(k) - rho_guess(k1)) / drho_dT(k1)
  enddo

! Refine the guesses for each layer.                                 !
  do itt = 1,6
    call calculate_density(T0, S0, pres, rho_guess, eqn_of_state)
    call calculate_density_derivs(T0, S0, pres, drho_dT, drho_dS, eqn_of_state)
    do k=1,nz
      T0(k) = T0(k) + (GV%Rlay(k) - rho_guess(k)) / drho_dT(k)
    enddo
  enddo

  do k=1,nz ; do i=is,ie ; do j=js,je
    T(i,j,k) = T0(k)
    S(i,j,k) = S0(k)
  enddo ; enddo ; enddo
  PI = 4.0*atan(1.0)
  do i=is,ie ; do j=js,je
    SST = 0.5*(T0(k1)+T0(nz)) - 0.9*0.5*(T0(k1)-T0(nz)) * &
                               cos(PI*(G%geoLatT(i,j)-G%south_lat)/(G%len_lat))
    do k=1,k1-1
      T(i,j,k) = SST
    enddo
  enddo ; enddo

end subroutine benchmark_init_temperature_salinity

end module benchmark_initialization
