module DOME2d_initialization
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
!*  The module configures the model for a 2-d version of the "DOME"    *
!*  experiment.   DOME = Dynamics of Overflows and Mixing Experiment   *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : close_file, create_file, fieldtype, file_exists
use MOM_io, only : open_file, read_data, read_axis_data, SINGLE_FILE
use MOM_io, only : write_field, slasher, vardesc
use MOM_variables, only : thermo_var_ptrs, ocean_OBC_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA

implicit none ; private

#include <MOM_memory.h>

! -----------------------------------------------------------------------------
! Private (module-wise) parameters
! -----------------------------------------------------------------------------
real, parameter :: dome2d_width_bay = 0.1;      ! width of bay
real, parameter :: dome2d_width_bottom = 0.3;   ! width of deepest part
real, parameter :: dome2d_depth_bay = 0.2;      ! normalized bay depth

character(len=40) :: mod = "DOEM2D_initialization" ! This module's name.

! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public DOME2d_initialize_topography
public DOME2d_initialize_thickness
public DOME2d_initialize_temperature_salinity 

! -----------------------------------------------------------------------------
! This module contains the following routines
! -----------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
! Initialization of topography
!------------------------------------------------------------------------------
subroutine DOME2d_initialize_topography ( D, G, param_file, max_depth )
  ! Arguments 
  real, dimension(NIMEM_,NJMEM_), intent(out) :: D
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  real,                  intent(in) :: max_depth
  
  ! Local variables 
  integer   :: i, j
  real      :: x, bay_depth, l1, l2
  
  ! location where downslope starts
  l1 = dome2d_width_bay
  
  ! location where downslope reaches maximum depth 
  l2 = 1.0 - dome2d_width_bottom
  
  bay_depth = dome2d_depth_bay

  do i=G%isc,G%iec 
    do j=G%jsc,G%jec 
    
      ! Compute normalized zonal coordinate
      x = G%geoLonT(i,j) / G%len_lon;
    
      if ( x .le. l1 ) then
        D(i,j) = bay_depth * max_depth
      else if (( x .gt. l1 ) .and. ( x .lt. l2 )) then
        D(i,j) = bay_depth * max_depth + (1.0-bay_depth) * max_depth * &
                 ( x - l1 ) / (l2 - l1)
      else
        D(i,j) = max_depth
      end if         
      
    enddo
  enddo
end subroutine DOME2d_initialize_topography

!------------------------------------------------------------------------------
! Initialization of thicknesses
!------------------------------------------------------------------------------
subroutine DOME2d_initialize_thickness ( h, G, param_file )

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
  real    :: x
  real    :: delta_h
  real    :: min_thickness
  character(len=40) :: verticalCoordinate

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call MOM_mesg("MOM_initialization.F90, initialize_thickness_uniform: setting thickness")

  call get_param(param_file,mod,"MIN_THICKNESS",min_thickness,'Minimum layer thickness',units='m',default=1.e-3)
  call get_param(param_file,mod,"REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
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

    case ( REGRIDDING_LAYER, REGRIDDING_RHO )
  
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
     
         x = G%geoLonT(i,j) / G%len_lon;
         if ( x .le. dome2d_width_bay ) then
           h(i,j,1:nz-1) = G%GV%Angstrom;
           h(i,j,nz) = dome2d_depth_bay * G%max_depth - (nz-1) * G%GV%Angstrom;
         end if
      
      end do ; end do   
    
 !  case ( IC_RHO_C )
 !  
 !    do j=js,je ; do i=is,ie
 !        eta1D(nz+1) = -1.0*G%bathyT(i,j)
 !        do k=nz,1,-1
 !          eta1D(k) = e0(k)
 !          if (eta1D(k) < (eta1D(k+1) + min_thickness)) then
 !            eta1D(k) = eta1D(k+1) + min_thickness
 !            h(i,j,k) = min_thickness
 !          else
 !            h(i,j,k) = eta1D(k) - eta1D(k+1)
 !          endif
 !       enddo
 !   
 !       x = G%geoLonT(i,j) / G%len_lon;
 !       if ( x .le. dome2d_width_bay ) then
 !         h(i,j,1:nz-1) = min_thickness;
 !         h(i,j,nz) = dome2d_depth_bay * G%max_depth - (nz-1) * min_thickness;
 !       end if
 !   
 !    enddo ; enddo

    case ( REGRIDDING_ZSTAR )
    
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

    case ( REGRIDDING_SIGMA )
      do j=js,je ; do i=is,ie
        delta_h = G%bathyT(i,j) / nz;
        h(i,j,:) = delta_h;
      end do ; end do   

    case default
      call MOM_error(FATAL,"dome2d_initialize: "// &
      "Unrecognized i.c. setup - set REGRIDDING_COORDINATE_MODE")

  end select

end subroutine DOME2d_initialize_thickness


!------------------------------------------------------------------------------
! Initialization of temperature and salinity
!------------------------------------------------------------------------------
subroutine DOME2d_initialize_temperature_salinity ( T, S, h, G, param_file, &
                                                    eqn_of_state)
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(out) :: T, S
  real, intent(in), dimension(NIMEM_,NJMEM_, NKMEM_)  :: h
  type(ocean_grid_type),               intent(in)  :: G
  type(param_file_type),               intent(in)  :: param_file
  type(EOS_type),                      pointer     :: eqn_of_state

  integer   :: i, j, k, is, ie, js, je, nz
  real      :: x;
  integer   :: index_bay_z;
  real      :: delta_S, delta_T;
  real      :: S_ref, T_ref;        ! Reference salinity and temerature within
                                    ! surface layer
  real      :: S_range, T_range;    ! Range of salinities and temperatures over the
                                    ! vertical
  real      :: xi0, xi1;
  character(len=40) :: verticalCoordinate
  
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call get_param(param_file,mod,"REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
            default=DEFAULT_COORDINATE_MODE)
  call get_param(param_file,mod,"S_REF",S_ref,'Reference salinity',units='1e-3',fail_if_missing=.true.)
  call get_param(param_file,mod,"T_REF",T_ref,'Refernce temperature',units='C',fail_if_missing=.true.)
  call get_param(param_file,mod,"S_RANGE",S_range,'Initial salinity range',units='1e-3',default=2.0)
  call get_param(param_file,mod,"T_RANGE",T_range,'Initial temperature range',units='1e-3',default=0.0)
  
  T(:,:,:) = 0.0
  S(:,:,:) = 0.0
  
  ! Linear salinity profile

  select case ( coordinateMode(verticalCoordinate) )

    case ( REGRIDDING_ZSTAR, REGRIDDING_SIGMA )

      do j=js,je ; do i=is,ie
        xi0 = 0.0;
        do k = 1,nz
          xi1 = xi0 + h(i,j,k) / G%max_depth;
          S(i,j,k) = 34.0 + 0.5 * S_range * (xi0 + xi1);
          xi0 = xi1;
        enddo
      enddo ; enddo
    
    case ( REGRIDDING_RHO )

      do j=js,je ; do i=is,ie
        xi0 = 0.0;
        do k = 1,nz
          xi1 = xi0 + h(i,j,k) / G%max_depth;
          S(i,j,k) = 34.0 + 0.5 * S_range * (xi0 + xi1);
          xi0 = xi1;
        enddo
        x = G%geoLonT(i,j) / G%len_lon;
        if ( x .le. dome2d_width_bay ) then
          S(i,j,nz) = 34.0 + S_range;
        end if  
      enddo ; enddo

    case ( REGRIDDING_LAYER )

      delta_S = S_range / ( G%ke - 1.0 );
      S(:,:,1) = S_ref;
      do k = 2,G%ke
        S(:,:,k) = S(:,:,k-1) + delta_S;
      end do 
    
    case default
      call MOM_error(FATAL,"dome2d_initialize: "// &
      "Unrecognized i.c. setup - set REGRIDDING_COORDINATE_MODE")

  end select

  ! Modify salinity and temperature when z coordinates are used  
  if ( coordinateMode(verticalCoordinate) .eq. REGRIDDING_ZSTAR ) then
    index_bay_z = Nint ( dome2d_depth_bay * G%ke );
    do j = G%jsc,G%jec ; do i = G%isc,G%iec 
      x = G%geoLonT(i,j) / G%len_lon;
      if ( x .le. dome2d_width_bay ) then
        S(i,j,1:index_bay_z) = S_ref + S_range; ! Use for z coordinates
        T(i,j,1:index_bay_z) = 1.0;             ! Use for z coordinates
      end if    
    end do ; end do ! i and j loops
  end if ! Z initial conditions

  ! Modify salinity and temperature when sigma coordinates are used  
  if ( coordinateMode(verticalCoordinate) .eq. REGRIDDING_SIGMA ) then
    do i = G%isc,G%iec ; do j = G%jsc,G%jec
      x = G%geoLonT(i,j) / G%len_lon;
      if ( x .le. dome2d_width_bay ) then
        S(i,j,1:G%ke) = S_ref + S_range;    ! Use for sigma coordinates
        T(i,j,1:G%ke) = 1.0;                ! Use for sigma coordinates
      end if    
    end do ; end do
  end if
  
  ! Modify temperature when rho coordinates are used
  T(G%isc:G%iec,G%jsc:G%jec,1:G%ke) = 0.0;
  if (( coordinateMode(verticalCoordinate) .eq. REGRIDDING_RHO ) .or. ( coordinateMode(verticalCoordinate) .eq. REGRIDDING_LAYER )) then 
    do i = G%isc,G%iec ; do j = G%jsc,G%jec
      x = G%geoLonT(i,j) / G%len_lon;
      if ( x .le. dome2d_width_bay ) then
        T(i,j,G%ke) = 1.0;
      end if    
    end do ; end do
  end if
  
  
end subroutine DOME2d_initialize_temperature_salinity

end module DOME2d_initialization
