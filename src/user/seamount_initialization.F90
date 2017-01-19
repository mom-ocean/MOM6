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

use MOM_domains, only : sum_across_PEs
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : close_file, fieldtype, file_exists
use MOM_io, only : open_file, read_data, read_axis_data, SINGLE_FILE
use MOM_io, only : write_field, slasher, vardesc
use MOM_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use MOM_tracer_registry, only : tracer_registry_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA

implicit none ; private

#include <MOM_memory.h>

character(len=40) :: mod = "seamount_initialization" ! This module's name.

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

!> Initialization of topography.
subroutine seamount_initialize_topography ( D, G, param_file, max_depth )
  ! Arguments
  type(dyn_horgrid_type),             intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),              intent(in)  :: param_file !< Parameter file structure
  real,                               intent(in)  :: max_depth  !< Maximum depth of model in m

  ! Local variables
  integer   :: i, j
  real      :: x, y, delta, Lx, rLx, Ly, rLy

  call get_param(param_file,mod,"SEAMOUNT_DELTA",delta, &
                 "Non-dimensional height of seamount.", &
                 units="non-dim", default=0.5)
  call get_param(param_file,mod,"SEAMOUNT_X_LENGTH_SCALE",Lx, &
                 "Length scale of seamount in x-direction.\n"//&
                 "Set to zero make topography uniform in the x-direction.", &
                 units="Same as x,y", default=20.)
  call get_param(param_file,mod,"SEAMOUNT_Y_LENGTH_SCALE",Ly, &
                 "Length scale of seamount in y-direction.\n"//&
                 "Set to zero make topography uniform in the y-direction.", &
                 units="Same as x,y", default=0.)

  Lx = Lx / G%len_lon
  Ly = Ly / G%len_lat
  rLx = 0. ; if (Lx>0.) rLx = 1. / Lx
  rLy = 0. ; if (Ly>0.) rLy = 1. / Ly
  do i=G%isc,G%iec
    do j=G%jsc,G%jec
      ! Compute normalized zonal coordinates (x,y=0 at center of domain)
      x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon - 0.5
      y = ( G%geoLatT(i,j) - G%south_lat ) / G%len_lat - 0.5
      D(i,j) = G%max_depth * ( 1.0 - delta * exp(-(rLx*x)**2 -(rLy*y)**2) )
    enddo
  enddo

end subroutine seamount_initialize_topography

!> Initialization of thicknesses.
!! This subroutine initializes the layer thicknesses to be uniform.
subroutine seamount_initialize_thickness ( h, G, GV, param_file )
  type(ocean_grid_type), intent(in)           :: G          !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)         :: GV         !< The ocean's vertical grid structure.
  real, intent(out), dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h !< The thicknesses being
                                                            !! initialized.
  type(param_file_type), intent(in)           :: param_file !< A structure indicating the
                                                            !! open file to parse for model
                                                            !! parameter values.

  real :: e0(SZK_(G)+1)   ! The resting interface heights, in m, usually !
                          ! negative because it is positive upward.      !
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface !
                          ! positive upward, in m.                       !
  integer :: i, j, k, is, ie, js, je, nz
  real    :: x
  real    :: delta_h
  real    :: min_thickness, S_surf, S_range, S_ref, S_light, S_dense
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
  !do k=1,nz+1
  !  e0(k) = -G%max_depth * real(k-1) / real(nz)
  !enddo

  select case ( coordinateMode(verticalCoordinate) )

  case ( REGRIDDING_LAYER, REGRIDDING_RHO ) ! Initial thicknesses for isopycnal coordinates
    call get_param(param_file,mod,"INITIAL_SSS", S_surf, default=34., do_not_log=.true.)
    call get_param(param_file,mod,"INITIAL_S_RANGE", S_range, default=2., do_not_log=.true.)
    call get_param(param_file, mod, "S_REF", S_ref, default=35.0, do_not_log=.true.)
    call get_param(param_file, mod, "TS_RANGE_S_LIGHT", S_light, default = S_Ref, do_not_log=.true.)
    call get_param(param_file, mod, "TS_RANGE_S_DENSE", S_dense, default = S_Ref, do_not_log=.true.)
    do K=1,nz+1
      ! Salinity of layer k is S_light + (k-1)/(nz-1) * (S_dense - S_light)
      ! Salinity of interface K is S_light + (K-3/2)/(nz-1) * (S_dense - S_light)
      ! Salinity at depth z should be S(z) = S_surf - S_range * z/max_depth
      ! Equating: S_surf - S_range * z/max_depth = S_light + (K-3/2)/(nz-1) * (S_dense - S_light)
      ! Equating: - S_range * z/max_depth = S_light - S_surf + (K-3/2)/(nz-1) * (S_dense - S_light)
      ! Equating: z/max_depth = - ( S_light - S_surf + (K-3/2)/(nz-1) * (S_dense - S_light) ) / S_range
      e0(K) = - G%max_depth * ( ( S_light  - S_surf ) + ( S_dense - S_light ) * ( (real(K)-1.5) / real(nz-1) ) ) / S_range
      e0(K) = nint(2048.*e0(K))/2048. ! Force round numbers ... the above expression has irrational factors ...
      e0(K) = min(real(1-K)*GV%Angstrom_z, e0(K)) ! Bound by surface
      e0(K) = max(-G%max_depth, e0(K)) ! Bound by bottom
    enddo
    do j=js,je ; do i=is,ie
      eta1D(nz+1) = -1.0*G%bathyT(i,j)
      do k=nz,1,-1
        eta1D(k) = e0(k)
        if (eta1D(k) < (eta1D(k+1) + GV%Angstrom_z)) then
          eta1D(k) = eta1D(k+1) + GV%Angstrom_z
          h(i,j,k) = GV%Angstrom_z
        else
          h(i,j,k) = eta1D(k) - eta1D(k+1)
        endif
      enddo
    enddo ; enddo

  case ( REGRIDDING_ZSTAR )                       ! Initial thicknesses for z coordinates
    do j=js,je ; do i=is,ie
      eta1D(nz+1) = -1.0*G%bathyT(i,j)
      do k=nz,1,-1
        eta1D(k) =  -G%max_depth * real(k-1) / real(nz)
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
      delta_h = G%bathyT(i,j) / dfloat(nz)
      h(i,j,:) = delta_h
    end do ; end do
end select

end subroutine seamount_initialize_thickness

!> Initial values for temperature and salinity
subroutine seamount_initialize_temperature_salinity ( T, S, h, G, GV, param_file, eqn_of_state)
  type(ocean_grid_type),                     intent(in)  :: G !< Ocean grid structure
  type(verticalGrid_type),                   intent(in) :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: T !< Potential temperature (degC)
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: S !< Salinity (ppt)
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(in)  :: h !< Layer thickness (m or Pa)
  type(param_file_type),                     intent(in)  :: param_file !< Parameter file structure
  type(EOS_type),                            pointer     :: eqn_of_state !< Equation of state structure
  ! Local variables
  integer :: i, j, k, is, ie, js, je, nz, k_light
  real    :: xi0, xi1, dxi, r, S_surf, T_surf, S_range, T_range
  real    :: T_ref, T_Light, T_Dense, S_ref, S_Light, S_Dense, a1, frac_dense, k_frac, res_rat
  character(len=20) :: verticalCoordinate, density_profile

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call get_param(param_file, mod, "REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE)
  call get_param(param_file,mod,"INITIAL_DENSITY_PROFILE", density_profile, &
                 'Initial profile shape. Valid values are "linear", "parabolic"\n'// &
                 'and "exponential".', default='linear')
  call get_param(param_file,mod,"INITIAL_SSS", S_surf, &
                 'Initial surface salinity', units='1e-3', default=34.)
  call get_param(param_file,mod,"INITIAL_SST", T_surf, &
                 'Initial surface temperature', units='C', default=0.)
  call get_param(param_file,mod,"INITIAL_S_RANGE", S_range, &
                 'Initial salinity range (bottom - surface)', units='1e-3', default=2.)
  call get_param(param_file,mod,"INITIAL_T_RANGE", T_range, &
                 'Initial temperature range (bottom - surface)', units='C', default=0.)

  select case ( coordinateMode(verticalCoordinate) )
    case ( REGRIDDING_LAYER ) ! Initial thicknesses for layer isopycnal coordinates
      ! These parameters are used in MOM_fixed_initialization.F90 when CONFIG_COORD="ts_range"
      call get_param(param_file, mod, "T_REF", T_ref, default=10.0, do_not_log=.true.)
      call get_param(param_file, mod, "TS_RANGE_T_LIGHT", T_light, default=T_Ref, do_not_log=.true.)
      call get_param(param_file, mod, "TS_RANGE_T_DENSE", T_dense, default=T_Ref, do_not_log=.true.)
      call get_param(param_file, mod, "S_REF", S_ref, default=35.0, do_not_log=.true.)
      call get_param(param_file, mod, "TS_RANGE_S_LIGHT", S_light, default = S_Ref, do_not_log=.true.)
      call get_param(param_file, mod, "TS_RANGE_S_DENSE", S_dense, default = S_Ref, do_not_log=.true.)
      call get_param(param_file, mod, "TS_RANGE_RESOLN_RATIO", res_rat, default=1.0, do_not_log=.true.)
      ! Emulate the T,S used in the "ts_range" coordinate configuration code
      k_light = GV%nk_rho_varies + 1
      do j=js,je ; do i=is,ie
        T(i,j,k_light) = T_light ; S(i,j,k_light) = S_light
      enddo ; enddo
      a1 = 2.0 * res_rat / (1.0 + res_rat)
      do k=k_light+1,nz
        k_frac = real(k-k_light)/real(nz-k_light)
        frac_dense = a1 * k_frac + (1.0 - a1) * k_frac**2
        do j=js,je ; do i=is,ie
          T(i,j,k) = frac_dense * (T_Dense - T_Light) + T_Light
          S(i,j,k) = frac_dense * (S_Dense - S_Light) + S_Light
        enddo ; enddo
      enddo
    case ( REGRIDDING_SIGMA, REGRIDDING_ZSTAR, REGRIDDING_RHO ) ! All other coordinate use FV initialization
      do j=js,je ; do i=is,ie
        xi0 = 0.0
        do k = 1,nz
          xi1 = xi0 + h(i,j,k) / G%max_depth
          select case ( trim(density_profile) )
            case ('linear')
             !S(i,j,k) = S_surf + S_range * 0.5 * (xi0 + xi1)
              S(i,j,k) = S_surf + ( 0.5 * S_range ) * (xi0 + xi1) ! Coded this way to reproduce old hard-coded answers
              T(i,j,k) = T_surf + T_range * 0.5 * (xi0 + xi1)
            case ('parabolic')
              S(i,j,k) = S_surf + S_range * (2.0 / 3.0) * (xi1**3 - xi0**3) / (xi1 - xi0)
              T(i,j,k) = T_surf + T_range * (2.0 / 3.0) * (xi1**3 - xi0**3) / (xi1 - xi0)
            case ('exponential')
              S(i,j,k) = S_surf + S_range * (exp(xi1/r)-exp(xi0/r)) / (xi1 - xi0)
              T(i,j,k) = T_surf + T_range * (exp(xi1/r)-exp(xi0/r)) / (xi1 - xi0)
            case default
              call MOM_error(FATAL, 'Unknown value for "INITIAL_DENSITY_PROFILE"')
          end select
          xi0 = xi1
        enddo
      enddo ; enddo
  end select

end subroutine seamount_initialize_temperature_salinity

!> \class seamount_initialization
!!
!! The module configures the model for the idealized seamount
!! test case.
end module seamount_initialization
