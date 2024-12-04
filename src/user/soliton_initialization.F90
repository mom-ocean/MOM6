!> Initial conditions for the Equatorial Rossby soliton test (Boyd).
module soliton_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA

implicit none ; private

#include <MOM_memory.h>

! Private (module-wise) parameters
character(len=40) :: mdl = "soliton_initialization" !< This module's name.

public soliton_initialize_thickness
public soliton_initialize_velocity

contains

!> Initialization of thicknesses in equatorial Rossby soliton test, as described in section
!! 6.1 of Haidvogel and Beckman (1990) and in Boyd (1980, JPO) and Boyd (1985, JPO).
subroutine soliton_initialize_thickness(h, depth_tot, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h    !< The thickness that is being initialized [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.

  ! Local variables
  real    :: max_depth  ! Maximum depth of the model bathymetry [Z ~> m]
  real    :: cg_max     ! The external wave speed based on max_depth [L T-1 ~> m s-1]
  real    :: beta       ! The meridional gradient of the Coriolis parameter [T-1 L-1 ~> s-1 m-1]
  real    :: L_eq       ! The equatorial deformation radius used in nondimensionalizing this problem [L ~> m]
  real    :: scale_pos  ! A conversion factor to nondimensionalize the axis units, usually [m-1]
  real    :: x0    ! Initial x-position of the soliton in the same units as geoLonT, often [m].
  real    :: y0    ! Initial y-position of the soliton in the same units as geoLatT, often [m].
  real    :: x, y  ! Nondimensionalized positions [nondim]
  real    :: I_nz  ! The inverse of the number of layers [nondim]
  real    :: val1  ! A nondimensionlized zonal decay scale [nondim]
  real    :: val2  ! An overall surface height anomaly amplitude [L T-1 ~> m s-1]
  real    :: val3  ! A decay factor [nondim]
  real    :: val4  ! The local velocity amplitude [L T-1 ~> m s-1]
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.just_read) &
    call MOM_mesg("soliton_initialization.F90, soliton_initialize_thickness: setting thickness")

  if (.not.just_read) call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "MAXIMUM_DEPTH", max_depth, &
                 units="m", default=-1.e9, scale=US%m_to_Z, do_not_log=.true.)
  call get_param(param_file, mdl, "BETA", beta, &
                 "The northward gradient of the Coriolis parameter with the betaplane option.", &
                 units="m-1 s-1", default=0.0, scale=US%T_to_s*US%L_to_m, do_not_log=.true.)

  if (just_read) return ! All run-time parameters have been read, so return.

  if (max_depth <= 0.0) call MOM_error(FATAL, &
      "soliton_initialization, soliton_initialize_thickness: "//&
      "This module requires a positive value of MAXIMUM_DEPTH.")
  if (abs(beta) <= 0.0) call MOM_error(FATAL, &
      "soliton_initialization, soliton_initialize_thickness: "//&
      "This module requires a non-zero value of BETA.")
  if (G%grid_unit_to_L <= 0.) call MOM_error(FATAL, "soliton_initialization.F90: "//&
          "soliton_initialize_thickness() is only set to work with Cartesian axis units.")

  cg_max = sqrt(GV%g_Earth * max_depth)
  L_eq = sqrt(cg_max / abs(beta))
  scale_pos = G%grid_unit_to_L / L_eq
  I_nz = 1.0 / real(nz)

  x0 = 2.0*G%len_lon/3.0
  y0 = 0.0
  val1 = 0.395
  val2 = max_depth * 0.771*(val1*val1)

  do j = G%jsc,G%jec ; do i = G%isc,G%iec
    do k = 1, nz
      x = (G%geoLonT(i,j)-x0) * scale_pos
      y = (G%geoLatT(i,j)-y0) * scale_pos
      val3 = exp(-val1*x)
      val4 = val2 * ( 2.0*val3 / (1.0 + (val3*val3)) )**2
      h(i,j,k) = (0.25*val4*(6.0*y*y + 3.0) * exp(-0.5*y*y) + depth_tot(i,j)) * I_nz
    enddo
  enddo ; enddo

end subroutine soliton_initialize_thickness


!> Initialization of u and v in the equatorial Rossby soliton test, as described in section
!! 6.1 of Haidvogel and Beckman (1990) and in Boyd (1980, JPO) and Boyd (1985, JPO).
subroutine soliton_initialize_velocity(u, v, G, GV, US, param_file, just_read)
  type(ocean_grid_type),                      intent(in)  :: G  !< Grid structure
  type(verticalGrid_type),                    intent(in)  :: GV !< The ocean's vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out) :: u  !< i-component of velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: v  !< j-component of velocity [L T-1 ~> m s-1]
  type(unit_scale_type),                      intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.

  ! Local variables
  real    :: max_depth  ! Maximum depth of the model bathymetry [Z ~> m]
  real    :: cg_max     ! The external wave speed based on max_depth [L T-1 ~> m s-1]
  real    :: beta       ! The meridional gradient of the Coriolis parameter [T-1 L-1 ~> s-1 m-1]
  real    :: L_eq       ! The equatorial deformation radius used in nondimensionalizing this problem [L ~> m]
  real    :: scale_pos  ! A conversion factor to nondimensionalize the axis units, usually [m-1]
  real    :: x0    ! Initial x-position of the soliton in the same units as geoLonT, often [m].
  real    :: y0    ! Initial y-position of the soliton in the same units as geoLatT, often [m].
  real    :: x, y  ! Nondimensionalized positions [nondim]
  real    :: val1  ! A nondimensionlized zonal decay scale [nondim]
  real    :: val2  ! An overall velocity amplitude [L T-1 ~> m s-1]
  real    :: val3  ! A decay factor [nondim]
  real    :: val4  ! The local velocity amplitude [L T-1 ~> m s-1]
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.just_read) &
    call MOM_mesg("soliton_initialization.F90, soliton_initialize_thickness: setting thickness")

  call get_param(param_file, mdl, "MAXIMUM_DEPTH", max_depth, &
                 units="m", default=-1.e9, scale=US%m_to_Z, do_not_log=.true.)
  call get_param(param_file, mdl, "BETA", beta, &
                 "The northward gradient of the Coriolis parameter with the betaplane option.", &
                 units="m-1 s-1", default=0.0, scale=US%T_to_s*US%L_to_m, do_not_log=.true.)

  if (just_read) return ! All run-time parameters have been read, so return.

  if (max_depth <= 0.0) call MOM_error(FATAL, &
      "soliton_initialization, soliton_initialize_velocity: "//&
      "This module requires a positive value of MAXIMUM_DEPTH.")
  if (abs(beta) <= 0.0) call MOM_error(FATAL, &
      "soliton_initialization, soliton_initialize_velocity: "//&
      "This module requires a non-zero value of BETA.")
  if (G%grid_unit_to_L <= 0.) call MOM_error(FATAL, "soliton_initialization.F90: "//&
          "soliton_initialize_velocity() is only set to work with Cartesian axis units.")

  cg_max = sqrt(GV%g_Earth * max_depth)
  L_eq = sqrt(cg_max / abs(beta))
  scale_pos = G%grid_unit_to_L / L_eq

  x0 = 2.0*G%len_lon/3.0
  y0 = 0.0
  val1 = 0.395
  val2 = cg_max * 0.771*(val1*val1)

  v(:,:,:) = 0.0
  u(:,:,:) = 0.0

  do j = G%jsc,G%jec ; do I = G%isc-1,G%iec+1
    do k = 1, nz
      x = (0.5*(G%geoLonT(i+1,j)+G%geoLonT(i,j))-x0) * scale_pos
      y = (0.5*(G%geoLatT(i+1,j)+G%geoLatT(i,j))-y0) * scale_pos
      val3 = exp(-val1*x)
      val4 = val2*((2.0*val3/(1.0+(val3*val3)))**2)
      u(I,j,k) = 0.25*val4*(6.0*y*y-9.0) * exp(-0.5*y*y)
    enddo
  enddo ; enddo
  do j = G%jsc-1,G%jec+1 ; do I = G%isc,G%iec
    do k = 1, nz
      x = 0.5*(G%geoLonT(i,j+1)+G%geoLonT(i,j))-x0
      y = 0.5*(G%geoLatT(i,j+1)+G%geoLatT(i,j))-y0
      val3 = exp(-val1*x)
      val4 = val2*((2.0*val3/(1.0+(val3*val3)))**2)
      v(i,J,k) = 2.0*val4*y*(-2.0*val1*tanh(val1*x)) * exp(-0.5*y*y)
    enddo
  enddo ; enddo

end subroutine soliton_initialize_velocity


!> \namespace soliton_initialization
!!
!! \section section_soliton Description of the equatorial Rossby soliton initial
!! conditions
!!

end module soliton_initialization
