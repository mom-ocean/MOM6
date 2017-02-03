module adjustment_initialization
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

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : close_file, fieldtype, file_exists
use MOM_io, only : open_file, read_data, read_axis_data, SINGLE_FILE
use MOM_io, only : write_field, slasher
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA

implicit none ; private

character(len=40) :: mod = "adjustment_initialization" ! This module's name.

#include <MOM_memory.h>

! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public adjustment_initialize_thickness
public adjustment_initialize_temperature_salinity

! -----------------------------------------------------------------------------
! This module contains the following routines
! -----------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!> Initialization of thicknesses.
!! This subroutine initializes the layer thicknesses to be uniform.
!------------------------------------------------------------------------------
subroutine adjustment_initialize_thickness ( h, G, GV, param_file )

  type(ocean_grid_type),   intent(in) :: G                    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV                   !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: h !< The thickness that is being initialized.
  type(param_file_type),   intent(in) :: param_file           !< A structure indicating the
                                         !! open file to parse for model parameter values.

  real :: e0(SZK_(G)+1)   ! The resting interface heights, in m, usually !
                          ! negative because it is positive upward.      !
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface !
                          ! positive upward, in m.                       !
  integer :: i, j, k, is, ie, js, je, nz
  real    :: x, y, yy, delta_S_strat, dSdz, delta_S, S_ref
  real    :: min_thickness, adjustment_width, adjustment_delta, adjustment_deltaS
  real    :: front_wave_amp, front_wave_length, front_wave_asym
  real    :: target_values(SZK_(G)+1)
  character(len=20) :: verticalCoordinate

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call MOM_mesg("initialize_thickness_uniform: setting thickness")

  ! Parameters used by main model initialization
  call get_param(param_file,mod,"S_REF",S_ref,fail_if_missing=.true.,do_not_log=.true.)
  call get_param(param_file,mod,"MIN_THICKNESS",min_thickness,'Minimum layer thickness', &
         units='m',default=1.0e-3)

  ! Parameters specific to this experiment configuration
  call get_param(param_file,mod,"REGRIDDING_COORDINATE_MODE",verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE)
  call get_param(param_file,mod,"ADJUSTMENT_WIDTH",adjustment_width,     &
                 "Width of frontal zone",                                &
                 units="same as x,y",fail_if_missing=.true.)
  call get_param(param_file,mod,"DELTA_S_STRAT",delta_S_strat,           &
                 "Top-to-bottom salinity difference of stratification",  &
                 units="1e-3",fail_if_missing=.true.)
  call get_param(param_file,mod,"ADJUSTMENT_DELTAS",adjustment_deltaS,   &
                 "Salinity difference across front",                     &
                 units="1e-3",fail_if_missing=.true.)
  call get_param(param_file,mod,"FRONT_WAVE_AMP",front_wave_amp,         &
                 "Amplitude of trans-frontal wave perturbation",         &
                 units="same as x,y",default=0.)
  call get_param(param_file,mod,"FRONT_WAVE_LENGTH",front_wave_length,   &
                 "Wave-length of trans-frontal wave perturbation",       &
                 units="same as x,y",default=0.)
  call get_param(param_file,mod,"FRONT_WAVE_ASYM",front_wave_asym,       &
                 "Amplitude of frontal asymmetric perturbation",         &
                 default=0.)

  ! WARNING: this routine specifies the interface heights so that the last layer
  !          is vanished, even at maximum depth. In order to have a uniform
  !          layer distribution, use this line of code within the loop:
  !          e0(k) = -G%max_depth * real(k-1) / real(nz)
  !          To obtain a thickness distribution where the last layer is
  !          vanished and the other thicknesses uniformly distributed, use:
  !          e0(k) = -G%max_depth * real(k-1) / real(nz-1)

  dSdz = -delta_S_strat/G%max_depth

  select case ( coordinateMode(verticalCoordinate) )

    case ( REGRIDDING_LAYER, REGRIDDING_RHO )
      if (delta_S_strat.ne.0.) then
        adjustment_delta = adjustment_deltaS / delta_S_strat * G%max_depth
        do k=1,nz+1
          e0(k) = adjustment_delta-(G%max_depth+2*adjustment_delta) * (real(k-1) / real(nz))
        enddo
      else
        adjustment_delta = 2.*G%max_depth
        do k=1,nz+1
          e0(k) = -(G%max_depth) * (real(k-1) / real(nz))
        enddo
      endif
      target_values(1)    = GV%Rlay(1)+0.5*(GV%Rlay(1)-GV%Rlay(2))
      target_values(nz+1) = GV%Rlay(nz)+0.5*(GV%Rlay(nz)-GV%Rlay(nz-1))
      do k = 2,nz
        target_values(k) = target_values(k-1) + ( GV%Rlay(nz) - GV%Rlay(1) ) / (nz-1)
      end do
      target_values = target_values - 1000.
      do j=js,je ; do i=is,ie
          if (front_wave_length.ne.0.) then
            y = ( 0.125 + G%geoLatT(i,j) / front_wave_length ) * ( 4. * acos(0.) )
            yy = 2. * ( G%geoLatT(i,j) - 0.5 * G%len_lat ) / adjustment_width
            yy = min(1.0, yy); yy = max(-1.0, yy)
            yy = yy * 2. * acos( 0. )
            y = front_wave_amp*sin(y) + front_wave_asym*sin(yy)
          else
            y = 0.
          endif
          x = ( ( G%geoLonT(i,j) - 0.5 * G%len_lon ) + y ) / adjustment_width
          x = min(1.0, x); x = max(-1.0, x)
          x = x * acos( 0. )
          delta_S = adjustment_deltaS * 0.5 * (1. - sin( x ) )
          do k=2,nz
            if (dSdz.ne.0.) then
              eta1D(k) = ( target_values(k) - ( S_ref + delta_S ) ) / dSdz
            else
              eta1D(k) = e0(k) - (0.5*adjustment_delta) * sin( x )
            endif
            eta1D(k) = max( eta1D(k), -G%max_depth )
            eta1D(k) = min( eta1D(k), 0. )
          enddo
          eta1D(1)=0.; eta1D(nz+1)=-G%max_depth
          do k=nz,1,-1
            if (eta1D(k) > 0.) then
              eta1D(k) = max( eta1D(k+1) + min_thickness, 0. )
              h(i,j,k) = max( eta1D(k) - eta1D(k+1), min_thickness )
            elseif (eta1D(k) <= (eta1D(k+1) + min_thickness)) then
              eta1D(k) = eta1D(k+1) + min_thickness
              h(i,j,k) = min_thickness
            else
              h(i,j,k) = eta1D(k) - eta1D(k+1)
            endif
          enddo
      enddo ; enddo

    case ( REGRIDDING_ZSTAR, REGRIDDING_SIGMA )
      do k=1,nz+1
        eta1D(k) = -(G%max_depth) * (real(k-1) / real(nz))
        eta1D(k) = max(min(eta1D(k),0.),-G%max_depth)
      enddo
      do j=js,je ; do i=is,ie
        do k=nz,1,-1
            h(i,j,k) = eta1D(k) - eta1D(k+1)
        enddo
      enddo ; enddo

    case default
      call MOM_error(FATAL,"adjustment_initialize_thickness: "// &
      "Unrecognized i.c. setup - set ADJUSTMENT_IC")

  end select

end subroutine adjustment_initialize_thickness


!------------------------------------------------------------------------------
!> Initialization of temperature and salinity.
!------------------------------------------------------------------------------
subroutine adjustment_initialize_temperature_salinity ( T, S, h, G, param_file, &
                                                    eqn_of_state)
  type(ocean_grid_type),   intent(in) :: G                    !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: T !< The temperature that is being initialized.
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: S !< The salinity that is being initialized.
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(in)  :: h !< The model thickness.
  type(param_file_type),   intent(in) :: param_file           !< A structure indicating the
                                         !! open file to parse for model parameter values.
  type(EOS_type),                 pointer     :: eqn_of_state !< Equation of state.

  integer   :: i, j, k, is, ie, js, je, nz
  real      :: x, y, yy
  integer   :: index_bay_z
  real      :: S_ref, T_ref         ! Reference salinity and temerature within
                                    ! surface layer
  real      :: S_range, T_range     ! Range of salinities and temperatures over the
                                    ! vertical
  real      :: xi0, xi1, dSdz, delta_S, delta_S_strat
  real      :: adjustment_width, adjustment_deltaS
  real       :: front_wave_amp, front_wave_length, front_wave_asym
  real      :: eta1d(SZK_(G)+1)
  character(len=20) :: verticalCoordinate

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  ! Parameters used by main model initialization
  call get_param(param_file,mod,"S_REF",S_ref,'Reference salinity',units='1e-3',fail_if_missing=.true.)
  call get_param(param_file,mod,"T_REF",T_ref,'Reference temperature',units='C',fail_if_missing=.true.)
  call get_param(param_file,mod,"S_RANGE",S_range,'Initial salinity range',units='1e-3', &
                 default=2.0)
  call get_param(param_file,mod,"T_RANGE",T_range,'Initial temperature range',units='C', &
                 default=0.0)
  ! Parameters specific to this experiment configuration BUT logged in previous s/r
  call get_param(param_file,mod,"REGRIDDING_COORDINATE_MODE",verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE)
  call get_param(param_file,mod,"ADJUSTMENT_WIDTH",adjustment_width,fail_if_missing=.true.,do_not_log=.true.)
  call get_param(param_file,mod,"ADJUSTMENT_DELTAS",adjustment_deltaS,fail_if_missing=.true.,do_not_log=.true.)
  call get_param(param_file,mod,"DELTA_S_STRAT",delta_S_strat,fail_if_missing=.true.,do_not_log=.true.)
  call get_param(param_file,mod,"FRONT_WAVE_AMP",front_wave_amp,default=0.,do_not_log=.true.)
  call get_param(param_file,mod,"FRONT_WAVE_LENGTH",front_wave_length,default=0.,do_not_log=.true.)
  call get_param(param_file,mod,"FRONT_WAVE_ASYM",front_wave_asym,default=0.,do_not_log=.true.)

  T(:,:,:) = 0.0
  S(:,:,:) = 0.0

  ! Linear salinity profile
  select case ( coordinateMode(verticalCoordinate) )

    case ( REGRIDDING_ZSTAR, REGRIDDING_SIGMA )
      dSdz = -delta_S_strat/G%max_depth
      do j=js,je ; do i=is,ie
          eta1d(nz+1)=-G%bathyT(i,j)
          do k=nz,1,-1
            eta1d(k)=eta1d(k+1)+h(i,j,k)
          enddo
          if (front_wave_length.ne.0.) then
            y = ( 0.125 + G%geoLatT(i,j) / front_wave_length ) * ( 4. * acos(0.) )
            yy = 2. * ( G%geoLatT(i,j) - 0.5 * G%len_lat ) / front_wave_length
            yy = min(1.0, yy); yy = max(-1.0, yy)
            yy = yy * 2. * acos( 0. )
            y = front_wave_amp*sin(y) + front_wave_asym*sin(yy)
          else
            y = 0.
          endif
          x = ( ( G%geoLonT(i,j) - 0.5 * G%len_lon ) + y ) / adjustment_width
          x = min(1.0, x); x = max(-1.0, x)
          x = x * acos( 0. )
          delta_S = adjustment_deltaS * 0.5 * (1. - sin( x ) )
          do k=1,nz
            S(i,j,k) = S_ref + delta_S + 0.5 * ( eta1D(k)+eta1D(k+1) ) * dSdz
            x = abs(S(i,j,k) - 0.5*real(nz-1)/real(nz)*S_range)/S_range*real(2*nz)
            x = 1.-min(1., x)
            T(i,j,k) = x
         enddo
   !     x=sum(T(i,j,:)*h(i,j,:))
   !     T(i,j,:)=T(i,j,:)/x*(G%max_depth*1.5/real(nz))
      enddo ; enddo

    case ( REGRIDDING_LAYER, REGRIDDING_RHO )
      do k = 1,nz
        S(:,:,k) = S_ref + S_range * ( (real(k)-0.5) / real( nz ) )
   !    x = abs(S(1,1,k) - 0.5*real(nz-1)/real(nz)*S_range)/S_range*real(2*nz)
   !    x = 1.-min(1., x)
   !    T(:,:,k) = x
      end do

    case default
      call MOM_error(FATAL,"adjustment_initialize_temperature_salinity: "// &
      "Unrecognized i.c. setup - set ADJUSTMENT_IC")

  end select

end subroutine adjustment_initialize_temperature_salinity

!> \class adjustment_initialization
!!
!! The module configures the model for the geostrophic adjustment
!! test case.
end module adjustment_initialization
