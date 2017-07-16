!> Initial conditions for the Equatorial Rossby soliton test (Boyd).
module soliton_initialization

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
character(len=40) :: mod = "soliton_initialization" !< This module's name.

public soliton_initialize_thickness
public soliton_initialize_velocity

contains

!> Initialization of thicknesses in Equatorial Rossby soliton test
subroutine soliton_initialize_thickness(h, G)
  type(ocean_grid_type),   intent(in) :: G                 !< Grid structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: h !< Thickness

  integer :: i, j, k, is, ie, js, je, nz
  real    :: x, y, x0, y0
  real    :: val1, val2, val3, val4
  character(len=40) :: verticalCoordinate

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call MOM_mesg("soliton_initialization.F90, soliton_initialize_thickness: setting thickness")

  x0 = 2.0*G%len_lon/3.0
  y0 = 0.0
  val1 = 0.395
  val2 = 0.771*(val1*val1)

  do j = G%jsc,G%jec ; do i = G%isc,G%iec
    do k = 1, nz
      x = G%geoLonT(i,j)-x0
      y = G%geoLatT(i,j)-y0
      val3 = exp(-val1*x)
      val4 = val2*((2.0*val3/(1.0+(val3*val3)))**2)
      h(i,j,k) = 0.25*val4*(6.0*y*y+3.0)*                 &
                exp(-0.5*y*y)
    enddo
  end do ; end do

end subroutine soliton_initialize_thickness


!> Initialization of u and v in the equatorial Rossby soliton test
subroutine soliton_initialize_velocity(u, v, h, G)
  type(ocean_grid_type),                  intent(in)     :: G  !< Grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(out) :: u  !< i-component of velocity [m/s]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(out) :: v  !< j-component of velocity [m/s]
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(in)  :: h  !< Thickness [H]

  real    :: x, y, x0, y0
  real    :: val1, val2, val3, val4
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  x0 = 2.0*G%len_lon/3.0
  y0 = 0.0
  val1 = 0.395
  val2 = 0.771*(val1*val1)

  v(:,:,:) = 0.0
  u(:,:,:) = 0.0

  do j = G%jsc,G%jec ; do I = G%isc-1,G%iec+1
    do k = 1, nz
      x = 0.5*(G%geoLonT(i+1,j)+G%geoLonT(i,j))-x0
      y = 0.5*(G%geoLatT(i+1,j)+G%geoLatT(i,j))-y0
      val3 = exp(-val1*x)
      val4 = val2*((2.0*val3/(1.0+(val3*val3)))**2)
      u(I,j,k) = 0.25*val4*(6.0*y*y-9.0)*         &
                      exp(-0.5*y*y)
    enddo
  enddo ; enddo
  do j = G%jsc-1,G%jec+1 ; do I = G%isc,G%iec
    do k = 1, nz
      x = 0.5*(G%geoLonT(i,j+1)+G%geoLonT(i,j))-x0
      y = 0.5*(G%geoLatT(i,j+1)+G%geoLatT(i,j))-y0
      val3 = exp(-val1*x)
      val4 = val2*((2.0*val3/(1.0+(val3*val3)))**2)
      v(i,J,k) = 2.0*val4*y*(-2.0*val1*tanh(val1*x))*        &
                     exp(-0.5*y*y)
    enddo
  enddo ; enddo

end subroutine soliton_initialize_velocity


!> \namespace soliton_initialization
!!
!! \section section_soliton Description of the equatorial Rossby soliton initial
!! conditions
!!

end module soliton_initialization
