module sloshing_initialization
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

use MOM_domains, only : sum_across_PEs
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : close_file, fieldtype, file_exists
use MOM_io, only : open_file, read_data, read_axis_data, SINGLE_FILE
use MOM_io, only : write_field, slasher
use MOM_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use MOM_tracer_registry, only : tracer_registry_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type

implicit none ; private

#include <MOM_memory.h>

! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public sloshing_initialize_topography
public sloshing_initialize_thickness
public sloshing_initialize_temperature_salinity

character(len=40)  :: mod = "sloshing_initialization" ! This module's name.

! -----------------------------------------------------------------------------
! This module contains the following routines
! -----------------------------------------------------------------------------
contains

!> Initialization of topography.
subroutine sloshing_initialize_topography ( D, G, param_file, max_depth )
  type(dyn_horgrid_type),             intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),              intent(in)  :: param_file !< Parameter file structure
  real,                               intent(in)  :: max_depth  !< Maximum depth of model in m

  ! Local variables
  integer   :: i, j

  do i=G%isc,G%iec
    do j=G%jsc,G%jec

      D(i,j) = max_depth

    enddo
  enddo

end subroutine sloshing_initialize_topography


!> Initialization of thicknesses
!! This routine is called when THICKNESS_CONFIG is set to 'sloshing'
!!
!! This routine initializes layer positions to set off a sloshing motion in
!! the zonal direction in a rectangular basin. All layers have initially the
!! same thickness but all interfaces (except bottom and sea surface) are
!! displaced according to a half-period cosine, with maximum value on the
!! left and minimum value on the right. This sets off a regular sloshing motion.
subroutine sloshing_initialize_thickness ( h, G, GV, param_file )
  type(ocean_grid_type), intent(in)           :: G          !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)         :: GV         !< The ocean's vertical grid structure.
  real, intent(out), dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h !< The thicknesses being
                                                            !! initialized.
  type(param_file_type), intent(in)           :: param_file !< A structure indicating the
                                                            !! open file to parse for model
                                                            !! parameter values.

  real    :: displ(SZK_(G)+1)
  real    :: z_unif(SZK_(G)+1)
  real    :: z_inter(SZK_(G)+1)
  real    :: x
  real    :: a0
  real    :: deltah
  real    :: total_height
  real    :: weight_z
  real    :: x1, y1, x2, y2
  real    :: t
  integer :: n

  integer :: i, j, k, is, ie, js, je, nx, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  deltah = G%max_depth / nz

  ! Define thicknesses
  do j=G%jsc,G%jec ; do i=G%isc,G%iec

    ! Define uniform interfaces
    do k = 0,nz
      z_unif(k+1) = -real(k)/real(nz)
    end do

    ! 1. Define stratification
    n = 3
    do k = 1,nz+1

      ! Thin pycnocline in the middle
      !z_inter(k) = (2.0**(n-1)) * (z_unif(k) + 0.5)**n - 0.5

      ! Thin pycnocline in the middle (piecewise linear profile)
      x1 = 0.30; y1 = 0.48; x2 = 0.70; y2 = 0.52

      x = -z_unif(k)

      if ( x .le. x1 ) then
        t = y1*x/x1;
      else if ( (x .gt. x1 ) .and. ( x .lt. x2 )) then
        t = y1 + (y2-y1) * (x-x1) / (x2-x1)
      else
        t = y2 + (1.0-y2) * (x-x2) / (1.0-x2)
      end if

      t = - z_unif(k)

      z_inter(k) = -t * G%max_depth

    end do

    ! 2. Define displacement
    a0 = 75.0;      ! Displacement amplitude (meters)
    do k = 1,nz+1

      weight_z = - 4.0 * ( z_unif(k) + 0.5 )**2 + 1

      x = G%geoLonT(i,j) / G%len_lon
      displ(k) = a0 * cos(acos(-1.0)*x) + weight_z;

      if ( k .EQ. 1 ) then
        displ(k) = 0.0
      end if

      if ( k .EQ. nz+1 ) then
        displ(k) = 0.0
      end if

      z_inter(k) = z_inter(k) + displ(k)

    end do

    ! 3. The last interface must coincide with the seabed
    z_inter(nz+1) = -G%bathyT(i,j)

    ! Modify interface heights to make sure all thicknesses
    ! are strictly positive
    do k = nz,1,-1

      if ( z_inter(k) .LT. (z_inter(k+1) + GV%Angstrom) ) then
        z_inter(k) = z_inter(k+1) + GV%Angstrom
      end if

    end do

    ! 4. Define layers
    total_height = 0.0
    do k = 1,nz
      h(i,j,k) = z_inter(k) - z_inter(k+1)

      total_height = total_height + h(i,j,k)
    end do

  enddo ; enddo

end subroutine sloshing_initialize_thickness


!------------------------------------------------------------------------------
!> Initialization of temperature and salinity
!!
!! This subroutine initializes linear profiles for T and S according to
!! reference surface layer salinity and temperature and a specified range.
!! Note that the linear distribution is set up with respect to the layer
!! number, not the physical position).
subroutine sloshing_initialize_temperature_salinity ( T, S, h, G, param_file, &
                                                      eqn_of_state)
  type(ocean_grid_type),                     intent(in)  :: G !< Ocean grid structure.
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: T !< Potential temperature (degC).
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: S !< Salinity (ppt).
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(in)  :: h !< Layer thickness (m or Pa).
  type(param_file_type),                     intent(in)  :: param_file !< A structure indicating the
                                                            !! open file to parse for model
                                                            !! parameter values.
  type(EOS_type),                            pointer     :: eqn_of_state !< Equation of state structure.

  integer :: i, j, k, is, ie, js, je, nz
  real    :: delta_S, delta_T
  real    :: S_ref, T_ref;      ! Reference salinity and temerature within
                                ! surface layer
  real    :: S_range, T_range;  ! Range of salinities and temperatures over the
                                ! vertical
  integer :: kdelta
  real    :: deltah
  real    :: xi0, xi1
  character(len=40)  :: mod = "initialize_temp_salt_linear" ! This subroutine's
                                                            ! name.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  call get_param(param_file,mod,"S_REF",S_ref,'Reference value for salinity',units='1e-3',fail_if_missing=.true.)
  call get_param(param_file,mod,"T_REF",T_ref,'Refernce value for temperature',units='C',fail_if_missing=.true.)

  ! The default is to assume an increase by 2 for the salinity and a uniform
  ! temperature
  call get_param(param_file,mod,"S_RANGE",S_range,'Initial salinity range.',units='1e-3',default=2.0)
  call get_param(param_file,mod,"T_RANGE",T_range,'Initial temperature range',units='C',default=0.0)

  ! Prescribe salinity
  !delta_S = S_range / ( G%ke - 1.0 )

  !S(:,:,1) = S_ref
  !do k = 2,G%ke
  !  S(:,:,k) = S(:,:,k-1) + delta_S
  !end do

  deltah = G%max_depth / nz;
  do j=js,je ; do i=is,ie
    xi0 = 0.0
    do k = 1,nz
      xi1 = xi0 + deltah / G%max_depth
      S(i,j,k) = 34.0 + 0.5 * S_range * (xi0 + xi1)
      xi0 = xi1
    enddo
  enddo ; enddo

  ! Prescribe temperature
  delta_T = T_range / ( G%ke - 1.0 )

  T(:,:,1) = T_ref
  do k = 2,G%ke
    T(:,:,k) = T(:,:,k-1) + delta_T
  end do
  kdelta = 2
  T(:,:,G%ke/2 - (kdelta-1):G%ke/2 + kdelta) = 1.0

end subroutine sloshing_initialize_temperature_salinity

!> \class sloshing_initialization
!!
!! The module configures the model for the non-rotating sloshing
!! test case.
end module sloshing_initialization
