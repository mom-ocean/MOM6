module seamount_initialization
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

!***********************************************************************
!*                                                                     *
!*  The module configures the model for the idealized seamount         *
!* test case.                                                          *
!*                                                                     *
!***********************************************************************


use MOM_domains, only : sum_across_PEs
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : close_file, create_file, fieldtype, file_exists
use MOM_io, only : open_file, read_data, read_axis_data, SINGLE_FILE
use MOM_io, only : write_field, slasher, vardesc
use MOM_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use MOM_tracer_registry, only : tracer_registry_type, add_tracer_OBC_values
use MOM_variables, only : thermo_var_ptrs, ocean_OBC_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA

implicit none ; private

#include <MOM_memory.h>

character(len=40) :: mod = "seamount_initialization" ! This module's name.

! -----------------------------------------------------------------------------
! Private (module-wise) parameters
! -----------------------------------------------------------------------------
integer, parameter :: RHO_LINEAR = 0;
integer, parameter :: RHO_EXP = 1;
integer, parameter :: RHO_PARABOLIC = 2;

integer, parameter :: density_profile = RHO_LINEAR;

! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public seamount_initialize_topography
public seamount_initialize_thickness
public seamount_initialize_temperature_salinity 

! -----------------------------------------------------------------------------
! This module contains the following routines
! -----------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
! Initialization of topography
!------------------------------------------------------------------------------
subroutine seamount_initialize_topography ( D, G, param_file, max_depth )
  ! Arguments 
  real, dimension(NIMEM_,NJMEM_), intent(out) :: D
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  real,                  intent(in) :: max_depth
  
  ! Local variables 
  integer   :: i, j
  real      :: x, delta, L

  call get_param(param_file,mod,"SEAMOUNT_DELTA",delta, &
                 "Non-dimensional height of seamount.", &
                 units="non-dim", default=0.5)
  call get_param(param_file,mod,"SEAMOUNT_LENGTH_SCALE",L, &
                 "Length scale of seamount.", &
                 units="Same as x,y", default=20.)
  
  ! Domain extent in kilometers
  
  L = L / G%len_lon
  do i=G%isc,G%iec 
    do j=G%jsc,G%jec 
      ! Compute normalized zonal coordinate (x=0 at center of domain)
      x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon - 0.5
      D(i,j) = G%max_depth * ( 1.0 - delta * exp(-(x/L)**2) )
    enddo
  enddo

end subroutine seamount_initialize_topography

!------------------------------------------------------------------------------
! Initialization of thicknesses
!------------------------------------------------------------------------------
subroutine seamount_initialize_thickness ( h, G, param_file )

  real, intent(out), dimension(NIMEM_,NJMEM_, NKMEM_) :: h
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file

! Arguments: h - The thickness that is being initialized.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

!  This subroutine initializes the layer thicknesses to be uniform.
  real :: e0(SZK_(G))     ! The resting interface heights, in m, usually !
                          ! negative because it is positive upward.      !
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface !
                          ! positive upward, in m.                       !
  integer :: i, j, k, is, ie, js, je, nz
  real    :: x;
  real    :: delta_h;
  real    :: min_thickness;
  character(len=20) :: verticalCoordinate

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call MOM_mesg("MOM_initialization.F90, initialize_thickness_uniform: setting thickness")

  call get_param(param_file,mod,"MIN_THICKNESS",min_thickness,'Minimum thickness for layer',units='m',default=1.0e-3)
  call get_param(param_file,mod,"REGRIDDING_COORDINATE_MODE",verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE)

 
  ! WARNING: this routine specifies the interface heights so that the last layer
  !          is vanished, even at maximum depth. In order to have a uniform
  !          layer distribution, use this line of code within the loop:
  !          e0(k) = -G%max_depth * real(k-1) / real(nz)
  !          To obtain a thickness distribution where the last layer is 
  !          vanished and the other thicknesses uniformly distributed, use:
  !          e0(k) = -G%max_depth * real(k-1) / real(nz-1)
  do k=1,nz
    e0(k) = -G%max_depth * real(k-1) / real(nz)
  enddo

    
  select case ( coordinateMode(verticalCoordinate) )
    
  case ( REGRIDDING_LAYER, REGRIDDING_RHO ) ! Initial thicknesses for isopycnal coordinates
    do j=js,je ; do i=is,ie
      eta1D(nz+1) = -1.0*G%bathyT(i,j)
      do k=nz,1,-1
        eta1D(k) = e0(k)
        if (eta1D(k) < (eta1D(k+1) + G%GV%Angstrom_z)) then
          eta1D(k) = eta1D(k+1) + G%GV%Angstrom_z
          h(i,j,k) = G%GV%Angstrom_z
        else
          h(i,j,k) = eta1D(k) - eta1D(k+1)
        endif
      enddo
    enddo ; enddo

  case ( REGRIDDING_ZSTAR )                       ! Initial thicknesses for z coordinates
    do j=js,je ; do i=is,ie
      eta1D(nz+1) = -1.0*G%bathyT(i,j)
      do k=nz,1,-1
        eta1D(k) = e0(k)
        if (eta1D(k) < (eta1D(k+1) + min_thickness)) then
          eta1D(k) = eta1D(k+1) + min_thickness
          h(i,j,k) = min_thickness
        else
          h(i,j,k) = eta1D(k) - eta1D(k+1)
        endif
      enddo
   enddo ; enddo

  case ( REGRIDDING_SIGMA )             ! Initial thicknesses for sigma coordinates
    do j=js,je ; do i=is,ie
      delta_h = G%bathyT(i,j) / dfloat(nz);
      h(i,j,:) = delta_h;
    end do ; end do 
end select
    
end subroutine seamount_initialize_thickness

!------------------------------------------------------------------------------
! Initialization of temperature and salinity
!------------------------------------------------------------------------------
subroutine seamount_initialize_temperature_salinity ( T, S, h, G, param_file, &
                                                      eqn_of_state)
                                                      
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(out) :: T, S
  real, intent(in), dimension(NIMEM_,NJMEM_, NKMEM_)  :: h
  type(ocean_grid_type),               intent(in)  :: G
  type(param_file_type),               intent(in)  :: param_file
  type(EOS_type),                      pointer     :: eqn_of_state

  integer :: i, j, k, is, ie, js, je, nz
  real    :: xi0, xi1, dxi, r
  character(len=20) :: verticalCoordinate
  
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call get_param(param_file,mod,"REGRIDDING_COORDINATE_MODE",verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE)
  
  select case ( coordinateMode(verticalCoordinate) )
    case ( REGRIDDING_LAYER ) ! Initial thicknesses for layer isopycnal coordinates
      do k=1,nz ; do j=js,je ; do i=is,ie
        xi1 = dfloat(k)/dfloat(nz)
        S(i,j,k) = 34.0 + (xi1);
      enddo ; enddo ; enddo
    case ( REGRIDDING_SIGMA, REGRIDDING_ZSTAR, REGRIDDING_RHO ) ! All other coordinate use FV initialization
      do j=js,je ; do i=is,ie
        xi0 = 0.0;
        do k = 1,nz
          xi1 = xi0 + h(i,j,k) / G%max_depth;
          if ( density_profile .eq. RHO_LINEAR ) then 
            ! ---------------------------
            ! Linear density profile
            ! ---------------------------
            S(i,j,k) = 34.0 + (xi0 + xi1);
          else if ( density_profile .eq. RHO_PARABOLIC ) then
            ! ---------------------------
            ! Parabolic density profile
            ! ---------------------------
            S(i,j,k) = 34.0 + (2.0 / 3.0) * (xi1**3 - xi0**3) / (xi1 - xi0);
          else if ( density_profile .eq. RHO_EXP ) then
            ! ---------------------------
            ! Exponential density profile
            ! ---------------------------
            r = 0.8;    ! small values give sharp profiles
            S(i,j,k) = 34.0 + r * (exp(xi1/r)-exp(xi0/r)) / (xi1 - xi0);
          end if    
          xi0 = xi1;
        enddo
      enddo ; enddo
  end select
  
end subroutine seamount_initialize_temperature_salinity

end module seamount_initialization
