!> Initial conditions for the 2D Rossby front test
module Rossby_front_2d_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : close_file, fieldtype, file_exists
use MOM_io, only : open_file, read_data, read_axis_data, SINGLE_FILE
use MOM_io, only : write_field, slasher, vardesc
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA

implicit none ; private

#include <MOM_memory.h>

! Private (module-wise) parameters
character(len=40) :: mod = "Rossby_front_2d_initialization" !< This module's name.

public Rossby_front_initialize_thickness
public Rossby_front_initialize_temperature_salinity
public Rossby_front_initialize_velocity

! Parameters defining the initial conditions of this test case
real, parameter :: frontFractionalWidth = 0.5 !< Width of front as fraction of domain
real, parameter :: HMLmin = 0.25 !< Shallowest ML as fractional depth of ocean
real, parameter :: HMLmax = 0.75 !< Deepest ML as fractional depth of ocean

contains

!> Initialization of thicknesses in 2D Rossby front test
subroutine Rossby_front_initialize_thickness(h, G, GV, param_file )
  type(ocean_grid_type),   intent(in) :: G                 !< Grid structure
  type(verticalGrid_type), intent(in) :: GV                !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: h !< Thickness
  type(param_file_type),   intent(in) :: param_file        !< Parameter file handle

  integer :: i, j, k, is, ie, js, je, nz
  real    :: Tz, Dml, eta, stretch, h0
  real    :: min_thickness, T_range, dRho_dT
  character(len=40) :: verticalCoordinate

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call MOM_mesg("Rossby_front_2d_initialization.F90, Rossby_front_initialize_thickness: setting thickness")

  ! Read parameters needed to set thickness
  call get_param(param_file, mod, "MIN_THICKNESS", min_thickness, &
                 'Minimum layer thickness',units='m',default=1.e-3)
  call get_param(param_file, mod, "REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE)
  call get_param(param_file, mod, "T_RANGE", T_range, 'Initial temperature range', units='C', default=0.0)
  call get_param(param_file, mod, "DRHO_DT", dRho_dT, default=-0.2, do_not_log=.true.)

  Tz = T_range / G%max_depth

  select case ( coordinateMode(verticalCoordinate) )

    case (REGRIDDING_LAYER, REGRIDDING_RHO)
      do j = G%jsc,G%jec ; do i = G%isc,G%iec
        Dml = Hml( G, G%geoLatT(i,j) )
        eta = -( -dRho_DT / GV%Rho0 ) * Tz * 0.5 * ( Dml * Dml )
        stretch = ( ( G%max_depth + eta ) / G%max_depth )
        h0 = ( G%max_depth / real(nz) ) * stretch
        do k = 1, nz
          h(i,j,k) = h0
        enddo
      end do ; end do

    case (REGRIDDING_ZSTAR, REGRIDDING_SIGMA)
      do j = G%jsc,G%jec ; do i = G%isc,G%iec
        Dml = Hml( G, G%geoLatT(i,j) )
        eta = -( -dRho_DT / GV%Rho0 ) * Tz * 0.5 * ( Dml * Dml )
        stretch = ( ( G%max_depth + eta ) / G%max_depth )
        h0 = ( G%max_depth / real(nz) ) * stretch
        do k = 1, nz
          h(i,j,k) = h0
        enddo
      end do ; end do

    case default
      call MOM_error(FATAL,"Rossby_front_initialize: "// &
      "Unrecognized i.c. setup - set REGRIDDING_COORDINATE_MODE")

  end select

end subroutine Rossby_front_initialize_thickness


!> Initialization of temperature and salinity in the Rossby front test
subroutine Rossby_front_initialize_temperature_salinity(T, S, h, G, param_file, eqn_of_state)
  type(ocean_grid_type),                     intent(in)  :: G  !< Grid structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: T  !< Potential temperature [deg C]
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: S  !< Salinity [ppt]
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(in)  :: h  !< Thickness
  type(param_file_type),                     intent(in)  :: param_file   !< Parameter file handle
  type(EOS_type),                            pointer     :: eqn_of_state !< Equation of state structure

  integer   :: i, j, k, is, ie, js, je, nz
  real      :: T_ref, S_ref ! Reference salinity and temerature within surface layer
  real      :: T_range      ! Range of salinities and temperatures over the vertical
  real      :: y, zc, zi, dTdz
  character(len=40) :: verticalCoordinate
  real      :: PI                   ! 3.1415926... calculated as 4*atan(1)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call get_param(param_file,mod,"REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
            default=DEFAULT_COORDINATE_MODE)
  call get_param(param_file,mod,"S_REF",S_ref,'Reference salinity',units='1e-3',fail_if_missing=.true.)
  call get_param(param_file,mod,"T_REF",T_ref,'Reference temperature',units='C',fail_if_missing=.true.)
  call get_param(param_file,mod,"T_RANGE",T_range,'Initial temperature range',units='C',default=0.0)

  T(:,:,:) = 0.0
  S(:,:,:) = S_ref
  dTdz = T_range / G%max_depth

  do j = G%jsc,G%jec ; do i = G%isc,G%iec
    zi = 0.
    do k = 1, nz
      zi = zi - h(i,j,k)              ! Bottom interface position
      zc = zi - 0.5*h(i,j,k)          ! Position of middle of cell
      zc = min( zc, -Hml(G, G%geoLatT(i,j)) ) ! Bound by depth of mixed layer
      T(i,j,k) = T_ref + dTdz * zc ! Linear temperature profile
    enddo
  end do ; end do

end subroutine Rossby_front_initialize_temperature_salinity


!> Initialization of u and v in the Rossby front test
subroutine Rossby_front_initialize_velocity(u, v, h, G, GV, param_file)
  type(ocean_grid_type),                  intent(in)     :: G  !< Grid structure
  type(verticalGrid_type),                intent(in)     :: GV !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(out) :: u  !< i-component of velocity [m/s]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(out) :: v  !< j-component of velocity [m/s]
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(in)  :: h  !< Thickness [H]
  type(param_file_type),                  intent(in)     :: param_file !< A structure indicating the
                                                               !! open file to parse for model
                                                               !! parameter values.

  real    :: y              ! Non-dimensional coordinate across channel, 0..pi
  real    :: T_range        ! Range of salinities and temperatures over the vertical
  real    :: dUdT           ! Factor to convert dT/dy into dU/dz, g*alpha/f
  real    :: dRho_dT, zi, zc, zm, f, Ty, Dml, hAtU
  integer :: i, j, k, is, ie, js, je, nz
  character(len=40) :: verticalCoordinate

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call get_param(param_file, mod, "REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
            default=DEFAULT_COORDINATE_MODE)
  call get_param(param_file, mod, "T_RANGE", T_range, 'Initial temperature range', units='C', default=0.0)
  call get_param(param_file, mod, "DRHO_DT", dRho_dT, default=-0.2, do_not_log=.true.)

  v(:,:,:) = 0.0
  u(:,:,:) = 0.0

  do j = G%jsc,G%jec ; do I = G%isc-1,G%iec+1
    f = 0.5*( G%CoriolisBu(I,j) + G%CoriolisBu(I,j-1) )
    dUdT = 0.0 ; if (abs(f) > 0.0) &
      dUdT = ( GV%g_Earth * dRho_dT ) / ( f * GV%Rho0 )
    Dml = Hml( G, G%geoLatT(i,j) )
    Ty = dTdy( G, T_range, G%geoLatT(i,j) )
    zi = 0.
    do k = 1, nz
      hAtU = 0.5*(h(i,j,k)+h(i+1,j,k))
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


!> Analytic prescription of mixed layer depth in 2d Rossby front test
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


!> \namespace Rossby_front_2d_initialization
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
