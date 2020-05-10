!> Initial conditions for the 2D Rossby front test
module Rossby_front_2d_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA

implicit none ; private

#include <MOM_memory.h>

! Private (module-wise) parameters
character(len=40) :: mdl = "Rossby_front_2d_initialization" !< This module's name.
! This include declares and sets the variable "version".
#include "version_variable.h"

public Rossby_front_initialize_thickness
public Rossby_front_initialize_temperature_salinity
public Rossby_front_initialize_velocity

! Parameters defining the initial conditions of this test case
real, parameter :: frontFractionalWidth = 0.5 !< Width of front as fraction of domain [nondim]
real, parameter :: HMLmin = 0.25 !< Shallowest ML as fractional depth of ocean [nondim]
real, parameter :: HMLmax = 0.75 !< Deepest ML as fractional depth of ocean [nondim]

contains

!> Initialization of thicknesses in 2D Rossby front test
subroutine Rossby_front_initialize_thickness(h, G, GV, US, param_file, just_read_params)
  type(ocean_grid_type),   intent(in)  :: G           !< Grid structure
  type(verticalGrid_type), intent(in)  :: GV          !< Vertical grid structure
  type(unit_scale_type),   intent(in)  :: US          !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [H ~> m or kg m-2]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.

  integer :: i, j, k, is, ie, js, je, nz
  real    :: Tz, Dml, eta, stretch, h0
  real    :: min_thickness, T_range
  real    :: dRho_dT      ! The partial derivative of density with temperature [R degC-1 ~> kg m-3 degC-1]
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=40) :: verticalCoordinate

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) &
    call MOM_mesg("Rossby_front_2d_initialization.F90, Rossby_front_initialize_thickness: setting thickness")

  if (.not.just_read) call log_version(param_file, mdl, version, "")
  ! Read parameters needed to set thickness
  call get_param(param_file, mdl, "MIN_THICKNESS", min_thickness, &
                 'Minimum layer thickness',units='m',default=1.e-3, do_not_log=just_read)
  call get_param(param_file, mdl, "REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=just_read)
  call get_param(param_file, mdl, "T_RANGE", T_range, 'Initial temperature range', &
                 units='C', default=0.0, do_not_log=just_read)
  call get_param(param_file, mdl, "DRHO_DT", dRho_dT, default=-0.2, scale=US%kg_m3_to_R, do_not_log=.true.)

  if (just_read) return ! All run-time parameters have been read, so return.

  Tz = T_range / G%max_depth

  select case ( coordinateMode(verticalCoordinate) )

    case (REGRIDDING_LAYER, REGRIDDING_RHO)
      do j = G%jsc,G%jec ; do i = G%isc,G%iec
        Dml = Hml( G, G%geoLatT(i,j) )
        eta = -( -dRho_DT / GV%Rho0 ) * Tz * 0.5 * ( Dml * Dml )
        stretch = ( ( G%max_depth + eta ) / G%max_depth )
        h0 = ( G%max_depth / real(nz) ) * stretch
        do k = 1, nz
          h(i,j,k) = h0 * GV%Z_to_H
        enddo
      enddo ; enddo

    case (REGRIDDING_ZSTAR, REGRIDDING_SIGMA)
      do j = G%jsc,G%jec ; do i = G%isc,G%iec
        Dml = Hml( G, G%geoLatT(i,j) )
        eta = -( -dRho_DT / GV%Rho0 ) * Tz * 0.5 * ( Dml * Dml )
        stretch = ( ( G%max_depth + eta ) / G%max_depth )
        h0 = ( G%max_depth / real(nz) ) * stretch
        do k = 1, nz
          h(i,j,k) = h0 * GV%Z_to_H
        enddo
      enddo ; enddo

    case default
      call MOM_error(FATAL,"Rossby_front_initialize: "// &
      "Unrecognized i.c. setup - set REGRIDDING_COORDINATE_MODE")

  end select

end subroutine Rossby_front_initialize_thickness


!> Initialization of temperature and salinity in the Rossby front test
subroutine Rossby_front_initialize_temperature_salinity(T, S, h, G, GV, &
                   param_file, eqn_of_state, just_read_params)
  type(ocean_grid_type),                     intent(in)  :: G  !< Grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: T  !< Potential temperature [degC]
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: S  !< Salinity [ppt]
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(in)  :: h  !< Thickness [H ~> m or kg m-2]
  type(param_file_type),                     intent(in)  :: param_file   !< Parameter file handle
  type(EOS_type),                            pointer     :: eqn_of_state !< Equation of state structure
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing T & S.

  integer   :: i, j, k, is, ie, js, je, nz
  real      :: T_ref, S_ref ! Reference salinity and temerature within surface layer
  real      :: T_range      ! Range of salinities and temperatures over the vertical
  real      :: y, zc, zi, dTdz
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=40) :: verticalCoordinate
  real      :: PI                   ! 3.1415926... calculated as 4*atan(1)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  call get_param(param_file, mdl,"REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
            default=DEFAULT_COORDINATE_MODE, do_not_log=just_read)
  call get_param(param_file, mdl,"S_REF",S_ref,'Reference salinity', units='1e-3', &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl,"T_REF",T_ref,'Reference temperature',units='C',&
                 fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl,"T_RANGE",T_range,'Initial temperature range',&
                 units='C', default=0.0, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  T(:,:,:) = 0.0
  S(:,:,:) = S_ref
  dTdz = T_range / G%max_depth

  do j = G%jsc,G%jec ; do i = G%isc,G%iec
    zi = 0.
    do k = 1, nz
      zi = zi - h(i,j,k)              ! Bottom interface position
      zc = GV%H_to_Z * (zi - 0.5*h(i,j,k))    ! Position of middle of cell
      zc = min( zc, -Hml(G, G%geoLatT(i,j)) ) ! Bound by depth of mixed layer
      T(i,j,k) = T_ref + dTdz * zc ! Linear temperature profile
    enddo
  enddo ; enddo

end subroutine Rossby_front_initialize_temperature_salinity


!> Initialization of u and v in the Rossby front test
subroutine Rossby_front_initialize_velocity(u, v, h, G, GV, US, param_file, just_read_params)
  type(ocean_grid_type),      intent(in)  :: G  !< Grid structure
  type(verticalGrid_type),    intent(in)  :: GV !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                              intent(out) :: u  !< i-component of velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                              intent(out) :: v  !< j-component of velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), &
                              intent(in)  :: h  !< Thickness [H ~> m or kg m-2]
  type(unit_scale_type),      intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),      intent(in)  :: param_file !< A structure indicating the open file
                                                !! to parse for model parameter values.
  logical,          optional, intent(in)  :: just_read_params !< If present and true, this call
                                                !! will only read parameters without setting u & v.

  real    :: y            ! Non-dimensional coordinate across channel, 0..pi
  real    :: T_range      ! Range of salinities and temperatures over the vertical
  real    :: dUdT         ! Factor to convert dT/dy into dU/dz, g*alpha/f [L2 Z-1 T-1 degC-1 ~> m s-1 degC-1]
  real    :: dRho_dT      ! The partial derivative of density with temperature [R degC-1 ~> kg m-3 degC-1]
  real    :: Dml, zi, zc, zm ! Depths [Z ~> m].
  real    :: f            ! The local Coriolis parameter [T-1 ~> s-1]
  real    :: Ty           ! The meridional temperature gradient [degC L-1 ~> degC m-1]
  real    :: hAtU         ! Interpolated layer thickness [Z ~> m].
  integer :: i, j, k, is, ie, js, je, nz
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=40) :: verticalCoordinate

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  call get_param(param_file, mdl, "REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=just_read)
  call get_param(param_file, mdl, "T_RANGE", T_range, 'Initial temperature range', &
                 units='C', default=0.0, do_not_log=just_read)
  call get_param(param_file, mdl, "DRHO_DT", dRho_dT, default=-0.2, scale=US%kg_m3_to_R, do_not_log=.true.)

  if (just_read) return ! All run-time parameters have been read, so return.

  v(:,:,:) = 0.0
  u(:,:,:) = 0.0

  do j = G%jsc,G%jec ; do I = G%isc-1,G%iec+1
    f = 0.5* (G%CoriolisBu(I,j) + G%CoriolisBu(I,j-1) )
    dUdT = 0.0 ; if (abs(f) > 0.0) &
      dUdT = ( GV%g_Earth*dRho_dT ) / ( f * GV%Rho0 )
    Dml = Hml( G, G%geoLatT(i,j) )
    Ty = US%L_to_m*dTdy( G, T_range, G%geoLatT(i,j) )
    zi = 0.
    do k = 1, nz
      hAtU = 0.5*(h(i,j,k)+h(i+1,j,k)) * GV%H_to_Z
      zi = zi - hAtU              ! Bottom interface position
      zc = zi - 0.5*hAtU          ! Position of middle of cell
      zm = max( zc + Dml, 0. )    ! Height above bottom of mixed layer
      u(I,j,k) = dUdT * Ty * zm   ! Thermal wind starting at base of ML
    enddo
  enddo ; enddo

end subroutine Rossby_front_initialize_velocity

!> Pseudo coordinate across domain used by Hml() and dTdy()
!! returns a coordinate from -PI/2 .. PI/2 squashed towards the
!! center of the domain.
real function yPseudo( G, lat )
  type(ocean_grid_type), intent(in) :: G   !< Grid structure
  real,                  intent(in) :: lat !< Latitude
  ! Local
  real :: y, PI

  PI = 4.0 * atan(1.0)
  yPseudo = ( ( lat - G%south_lat ) / G%len_lat ) - 0.5 ! -1/2 .. 1/.2
  yPseudo = PI * max(-0.5, min(0.5, yPseudo / frontFractionalWidth))
end function yPseudo


!> Analytic prescription of mixed layer depth in 2d Rossby front test,
!! in the same units as G%max_depth
real function Hml( G, lat )
  type(ocean_grid_type), intent(in) :: G   !< Grid structure
  real,                  intent(in) :: lat !< Latitude
  ! Local
  real :: dHML, HMLmean

  dHML = 0.5 * ( HMLmax - HMLmin ) * G%max_depth
  HMLmean = 0.5 * ( HMLmin + HMLmax ) * G%max_depth
  Hml = HMLmean + dHML * sin( yPseudo(G, lat) )
end function Hml


!> Analytic prescription of mixed layer temperature gradient in 2d Rossby front test
real function dTdy( G, dT, lat )
  type(ocean_grid_type), intent(in) :: G     !< Grid structure
  real,                  intent(in) :: dT    !< Top to bottom temperature difference
  real,                  intent(in) :: lat   !< Latitude
  ! Local
  real :: PI, dHML, dHdy
  real :: km = 1.e3 ! AXIS_UNITS = 'k' (1000 m)

  PI = 4.0 * atan(1.0)
  dHML = 0.5 * ( HMLmax - HMLmin ) * G%max_depth
  dHdy = dHML * ( PI / ( frontFractionalWidth * G%len_lat * km ) ) * cos( yPseudo(G, lat) )
  dTdy = -( dT / G%max_depth ) * dHdy

end function dTdy


!> \namespace rossby_front_2d_initialization
!!
!! \section section_Rossby_front_2d Description of the 2d Rossby front initial conditions
!!
!! Consistent with a linear equation of state, the system has a constant stratification
!! below the mixed layer, stratified in temperature only. Isotherms are flat below the
!! mixed layer and vertical within. Salinity is constant. The mixed layer has a half sine
!! form so that there are no mixed layer or temperature gradients at the side walls.
!!
!! Below the mixed layer the potential temperature, \f$\theta(z)\f$, is given by
!! \f[ \theta(z) = \theta_0 - \Delta \theta \left( z + h_{ML} \right) \f]
!! where \f$ \theta_0 \f$ and \f$ \Delta \theta \f$ are external model parameters.
!!
!! The depth of the mixed layer, \f$H_{ML}\f$ is
!! \f[ h_{ML}(y) = h_{min} + \left( h_{max} - h_{min} \right) \cos{\pi y/L} \f].
!! The temperature in mixed layer is given by the reference temperature at \f$z=h_{ML}\f$
!! so that
!! \f[ \theta(y,z) =
!!     \theta_0 - \Delta \theta \left( z + h_{ML} \right) & \forall & z < h_{ML}(y) T.B.D.
!! \f]

end module Rossby_front_2d_initialization
