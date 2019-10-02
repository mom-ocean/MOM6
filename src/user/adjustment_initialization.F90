!> Configures the model for the geostrophic adjustment test case.
module adjustment_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
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

character(len=40) :: mdl = "adjustment_initialization" !< This module's name.

#include <MOM_memory.h>

public adjustment_initialize_thickness
public adjustment_initialize_temperature_salinity

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> Initializes the layer thicknesses in the adjustment test case
subroutine adjustment_initialize_thickness ( h, G, GV, US, param_file, just_read_params)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US          !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [H ~> m or kg m-2].
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.
  ! Local variables
  real :: e0(SZK_(G)+1)   ! The resting interface heights, in depth units [Z ~> m], usually
                          ! negative because it is positive upward.
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface
                          ! positive upward, in depth units [Z ~> m].
  real    :: dRho_dS      ! The partial derivative of density with salinity [R ppt-1 ~> kg m-3 ppt-1].
                          ! In this subroutine it is hard coded at 1.0 kg m-3 ppt-1.
  real    :: x, y, yy
  real    :: delta_S_strat, dSdz, delta_S, S_ref
  real    :: min_thickness, adjustment_width, adjustment_delta
  real    :: adjustment_deltaS
  real    :: front_wave_amp, front_wave_length, front_wave_asym
  real    :: target_values(SZK_(G)+1)  ! Target densities or density anomalies [R ~> kg m-3]
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=20) :: verticalCoordinate
! This include declares and sets the variable "version".
#include "version_variable.h"
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) &
    call MOM_mesg("initialize_thickness_uniform: setting thickness")

  ! Parameters used by main model initialization
  if (.not.just_read) call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl,"S_REF",S_ref,fail_if_missing=.true.,do_not_log=.true.)
  call get_param(param_file, mdl,"MIN_THICKNESS",min_thickness,'Minimum layer thickness', &
         units='m', default=1.0e-3, do_not_log=just_read, scale=US%m_to_Z)

  ! Parameters specific to this experiment configuration
  call get_param(param_file, mdl,"REGRIDDING_COORDINATE_MODE",verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=just_read)
  call get_param(param_file, mdl,"ADJUSTMENT_WIDTH",adjustment_width,     &
                 "Width of frontal zone",                                &
                 units="same as x,y", fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl,"DELTA_S_STRAT",delta_S_strat,           &
                 "Top-to-bottom salinity difference of stratification",  &
                 units="1e-3", fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl,"ADJUSTMENT_DELTAS",adjustment_deltaS,   &
                 "Salinity difference across front",                     &
                 units="1e-3", fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl,"FRONT_WAVE_AMP",front_wave_amp,         &
                 "Amplitude of trans-frontal wave perturbation",         &
                 units="same as x,y",default=0., do_not_log=just_read)
  call get_param(param_file, mdl,"FRONT_WAVE_LENGTH",front_wave_length,   &
                 "Wave-length of trans-frontal wave perturbation",       &
                 units="same as x,y",default=0., do_not_log=just_read)
  call get_param(param_file, mdl,"FRONT_WAVE_ASYM",front_wave_asym,       &
                 "Amplitude of frontal asymmetric perturbation",         &
                 default=0., do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  ! WARNING: this routine specifies the interface heights so that the last layer
  !          is vanished, even at maximum depth. In order to have a uniform
  !          layer distribution, use this line of code within the loop:
  !          e0(k) = -G%max_depth * real(k-1) / real(nz)
  !          To obtain a thickness distribution where the last layer is
  !          vanished and the other thicknesses uniformly distributed, use:
  !          e0(k) = -G%max_depth * real(k-1) / real(nz-1)

  dSdz = -delta_S_strat / G%max_depth

  select case ( coordinateMode(verticalCoordinate) )

    case ( REGRIDDING_LAYER, REGRIDDING_RHO )
      dRho_dS = 1.0 * US%kg_m3_to_R
      if (delta_S_strat /= 0.) then
        ! This was previously coded ambiguously.
        adjustment_delta = (adjustment_deltaS / delta_S_strat) * G%max_depth
        do k=1,nz+1
          e0(k) = adjustment_delta - (G%max_depth + 2*adjustment_delta) * (real(k-1) / real(nz))
        enddo
      else
        adjustment_delta = 2.*G%max_depth
        do k=1,nz+1
          e0(k) = -G%max_depth * (real(k-1) / real(nz))
        enddo
      endif
      target_values(1)    = ( GV%Rlay(1) + 0.5*(GV%Rlay(1)-GV%Rlay(2)) )
      target_values(nz+1) = ( GV%Rlay(nz) + 0.5*(GV%Rlay(nz)-GV%Rlay(nz-1)) )
      do k = 2,nz
        target_values(k) = target_values(k-1) + ( GV%Rlay(nz) - GV%Rlay(1) ) / (nz-1)
      enddo
      target_values(:) = target_values(:) - 1000.*US%kg_m3_to_R
      do j=js,je ; do i=is,ie
        if (front_wave_length /= 0.) then
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
          if (dRho_dS*dSdz /= 0.) then
            eta1D(k) = ( target_values(k) - dRho_dS*( S_ref + delta_S ) ) / (dRho_dS*dSdz)
          else
            eta1D(k) = e0(k) - (0.5*adjustment_delta) * sin( x )
          endif
          eta1D(k) = max( eta1D(k), -G%max_depth )
          eta1D(k) = min( eta1D(k), 0. )
        enddo
        eta1D(1) = 0.; eta1D(nz+1) = -G%max_depth
        do k=nz,1,-1
          if (eta1D(k) > 0.) then
            eta1D(k) = max( eta1D(k+1) + min_thickness, 0. )
            h(i,j,k) = GV%Z_to_H * max( eta1D(k) - eta1D(k+1), min_thickness )
          elseif (eta1D(k) <= (eta1D(k+1) + min_thickness)) then
            eta1D(k) = eta1D(k+1) + min_thickness
            h(i,j,k) = GV%Z_to_H * min_thickness
          else
            h(i,j,k) = GV%Z_to_H * (eta1D(k) - eta1D(k+1))
          endif
        enddo
      enddo ; enddo

    case ( REGRIDDING_ZSTAR, REGRIDDING_SIGMA )
      do k=1,nz+1
        eta1D(k) = -G%max_depth * (real(k-1) / real(nz))
        eta1D(k) = max(min(eta1D(k), 0.), -G%max_depth)
      enddo
      do j=js,je ; do i=is,ie
        do k=nz,1,-1
          h(i,j,k) = GV%Z_to_H * (eta1D(k) - eta1D(k+1))
        enddo
      enddo ; enddo

    case default
      call MOM_error(FATAL,"adjustment_initialize_thickness: "// &
                     "Unrecognized i.c. setup - set ADJUSTMENT_IC")

  end select

end subroutine adjustment_initialize_thickness

!> Initialization of temperature and salinity in the adjustment test case
subroutine adjustment_initialize_temperature_salinity(T, S, h, G, GV, param_file, &
                                                      eqn_of_state, just_read_params)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: T !< The temperature that is being initialized.
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: S !< The salinity that is being initialized.
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(in)  :: h !< The model thicknesses [H ~> m or kg m-2].
  type(param_file_type),   intent(in) :: param_file   !< A structure indicating the open file to
                                                      !! parse for model parameter values.
  type(EOS_type),                 pointer     :: eqn_of_state !< Equation of state.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing T & S.

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
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=20) :: verticalCoordinate

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  ! Parameters used by main model initialization
  call get_param(param_file, mdl,"S_REF",S_ref,'Reference salinity', units='1e-3', &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl,"T_REF",T_ref,'Reference temperature', units='C', &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl,"S_RANGE",S_range,'Initial salinity range', units='1e-3', &
                 default=2.0, do_not_log=just_read)
  call get_param(param_file, mdl,"T_RANGE",T_range,'Initial temperature range', units='C', &
                 default=0.0, do_not_log=just_read)
  ! Parameters specific to this experiment configuration BUT logged in previous s/r
  call get_param(param_file, mdl,"REGRIDDING_COORDINATE_MODE",verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=just_read)
  call get_param(param_file, mdl,"ADJUSTMENT_WIDTH", adjustment_width, &
                 fail_if_missing=.not.just_read, do_not_log=.true.)
  call get_param(param_file, mdl,"ADJUSTMENT_DELTAS", adjustment_deltaS, &
                 fail_if_missing=.not.just_read, do_not_log=.true.)
  call get_param(param_file, mdl,"DELTA_S_STRAT", delta_S_strat, &
                 fail_if_missing=.not.just_read, do_not_log=.true.)
  call get_param(param_file, mdl,"FRONT_WAVE_AMP", front_wave_amp, default=0., &
                 do_not_log=.true.)
  call get_param(param_file, mdl,"FRONT_WAVE_LENGTH",front_wave_length, &
                 default=0., do_not_log=.true.)
  call get_param(param_file, mdl,"FRONT_WAVE_ASYM", front_wave_asym, default=0., &
                 do_not_log=.true.)

  if (just_read) return ! All run-time parameters have been read, so return.

  T(:,:,:) = 0.0
  S(:,:,:) = 0.0

  ! Linear salinity profile
  select case ( coordinateMode(verticalCoordinate) )

    case ( REGRIDDING_ZSTAR, REGRIDDING_SIGMA )
      dSdz = -delta_S_strat / G%max_depth
      do j=js,je ; do i=is,ie
        eta1d(nz+1) = -G%bathyT(i,j)
        do k=nz,1,-1
          eta1d(k) = eta1d(k+1) + h(i,j,k)*GV%H_to_Z
        enddo
        if (front_wave_length /= 0.) then
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
          x = 1. - min(1., x)
          T(i,j,k) = x
        enddo
   !    x = GV%H_to_Z*sum(T(i,j,:)*h(i,j,:))
   !    T(i,j,:) = (T(i,j,:) / x) * (G%max_depth*1.5/real(nz))
      enddo ; enddo

    case ( REGRIDDING_LAYER, REGRIDDING_RHO )
      do k = 1,nz
        S(:,:,k) = S_ref + S_range * ( (real(k)-0.5) / real( nz ) )
   !    x = abs(S(1,1,k) - 0.5*real(nz-1)/real(nz)*S_range)/S_range*real(2*nz)
   !    x = 1.-min(1., x)
   !    T(:,:,k) = x
      enddo

    case default
      call MOM_error(FATAL,"adjustment_initialize_temperature_salinity: "// &
      "Unrecognized i.c. setup - set ADJUSTMENT_IC")

  end select

end subroutine adjustment_initialize_temperature_salinity

end module adjustment_initialization
