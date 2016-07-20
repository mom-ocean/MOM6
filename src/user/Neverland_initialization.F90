module Neverland_initialization
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

use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_tracer_registry, only : tracer_registry_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type

implicit none ; private

#include <MOM_memory.h>

public Neverland_initialize_topography
public Neverland_initialize_thickness

contains

! -----------------------------------------------------------------------------
!> This subroutine sets up the Neverland test case topography.
subroutine Neverland_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),             intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),              intent(in)  :: param_file !< Parameter file structure
  real,                               intent(in)  :: max_depth  !< Maximum depth of model in m

! This subroutine sets up the Neverland test case topography
  real :: min_depth            ! The minimum and maximum depths in m.
  real :: PI                   ! 3.1415926... calculated as 4*atan(1)
  real :: D0                   ! A constant to make the maximum     !
                               ! basin depth MAXIMUM_DEPTH.         !
  real :: x, y
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "Neverland_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call MOM_mesg("  Neverland_initialization.F90, Neverland_initialize_topography: setting topography", 5)

  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)

  PI = 4.0*atan(1.0)

!  Calculate the depth of the bottom.
  do i=is,ie  
  do j=js,je
    x=(G%geoLonT(i,j)-G%west_lon)/G%len_lon
    y=(G%geoLatT(i,j)-G%south_lat)/G%len_lat
!  This sets topography that has a reentrant channel to the south.

    D(i,j) = 1.0 - (1.2 * spike(x,0.2) + 1.2 * spike(x-1.0,0.2)) * spike(MIN(0.0,y-0.3),0.2) & !< South America
                -  1.2 * spike(x-0.5,0.2) * spike(MIN(0.0,y-0.55),0.2)                       & !< Africa  
                -  1.1 * spike(y-1,0.12) - 1.1 * spike(y,0.12)                               & !< The great northern wall and Antarctica
                -  1.2 * (spike(x,0.12)  + spike(x-1,0.12)) * spike(MAX(0.0,y-0.06),0.12)    & !< Antarctic Peninsula
                -  0.1 * (cosbell(x,0.1) + cosbell(x-1,0.1))                                   & !< Drake Passage ridge
                -  0.5 * cosbell(x-0.16,0.05) * (cosbell(y-0.18,0.13)**0.4)                  & !< Scotia Arc East 
                -  0.4 * (cosbell(x-0.09,0.08)**0.4) * cosbell(y-0.26,0.05)                  & !< Scotia Arc North
                -  0.4 * (cosbell(x-0.08,0.08)**0.4) * cosbell(y-0.1,0.05)                   & !< Scotia Arc South
                -  0.05 * cos(14*PI*x) * sin(14*PI*y)                                        & !< roughness
                -  0.05 * cos(20*PI*x) * cos(20*PI*y)                                         !< roughness 

    if (D(i,j) < 0.0) D(i,j) = 0.0
    D(i,j) = D(i,j) * max_depth
  enddo 
  enddo
  
end subroutine Neverland_initialize_topography
! -----------------------------------------------------------------------------

!> Returns the value of a cosine-bell function evaluated at x/L
 real function cosbell(x,L)

   real , intent(in) :: x       !< non-dimensional position
   real , intent(in) :: L       !< non-dimensional width
   real              :: PI      !< 3.1415926... calculated as 4*atan(1) 

   PI      = 4.0*atan(1.0)
   cosbell = 0.5 * (1 + cos(PI*MIN(ABS(x/L),1.0)))
 end function cosbell

!> Returns the value of a sin-spike function evaluated at x/L
 real function spike(x,L)

   real , intent(in) :: x       !< non-dimensional position
   real , intent(in) :: L       !< non-dimensional width
   real              :: PI      !< 3.1415926... calculated as 4*atan(1) 

   PI    = 4.0*atan(1.0)
   spike = (1 - sin(PI*MIN(ABS(x/L),0.5)))
 end function spike


! -----------------------------------------------------------------------------
!> This subroutine initializes layer thicknesses for the Neverland test case,
!! by finding the depths of interfaces in a specified latitude-dependent
!! temperature profile with an exponentially decaying thermocline on top of a
!! linear stratification.
subroutine Neverland_initialize_thickness(h, G, GV, param_file, eqn_of_state, P_ref)
  type(ocean_grid_type),   intent(in) :: G                    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV                   !< The ocean's vertical grid structure.
  real, intent(out), dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h !< The thickness that is being
                                                              !! initialized.
  type(param_file_type),   intent(in) :: param_file           !< A structure indicating the open
                                                              !! file to parse for model
                                                              !! parameter values.
  type(EOS_type),          pointer    :: eqn_of_state         !< integer that selects the
                                                              !! equation of state.
  real,                    intent(in) :: P_Ref                !< The coordinate-density
                                                              !! reference pressure in Pa.

  real :: e0(SZK_(G)+1)     ! The resting interface heights, in m, usually !
                            ! negative because it is positive upward.      !
  real :: e_pert(SZK_(G)+1) ! Interface height perturbations, positive     !
                            ! upward, in m.                                !
  real :: eta1D(SZK_(G)+1)  ! Interface height relative to the sea surface !
                            ! positive upward, in m.                       !
  real :: SST       !  The initial sea surface temperature, in deg C.
  real :: T_int     !  The initial temperature of an interface, in deg C.
  real :: ML_depth  !  The specified initial mixed layer depth, in m.
  real :: thermocline_scale ! The e-folding scale of the thermocline, in m.
  real, dimension(SZK_(G)) :: T0, pres, S0, rho_guess, drho, drho_dT, drho_dS
  real :: a_exp      ! The fraction of the overall stratification that is exponential.
  real :: I_ts, I_md ! Inverse lengthscales in m-1.
  real :: T_frac     ! A ratio of the interface temperature to the range
                     ! between SST and the bottom temperature.
  real :: err, derr_dz  ! The error between the profile's temperature and the
                     ! interface temperature for a given z and its derivative.
  real :: pi, z
  character(len=40)  :: mod = "Neverland_initialize_thickness" ! This subroutine's name.
  integer :: i, j, k, k1, is, ie, js, je, nz, itt

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call MOM_mesg("  Neverland_initialization.F90, Neverland_initialize_thickness: setting thickness", 5)

  k1 = GV%nk_rho_varies + 1

  ML_depth = 50.0
  thermocline_scale = 500.0
  a_exp = 0.9

! This block calculates T0(k) for the purpose of diagnosing where the 
! interfaces will be found.
  do k=1,nz
    pres(k) = P_Ref ; S0(k) = 35.0
  enddo
  T0(k1) = 29.0
  call calculate_density(T0(k1),S0(k1),pres(k1),rho_guess(k1),eqn_of_state)
  call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,k1,1,eqn_of_state)

! A first guess of the layers' temperatures.
  do k=1,nz
    T0(k) = T0(k1) + (GV%Rlay(k) - rho_guess(k1)) / drho_dT(k1)
  enddo

! Refine the guesses for each layer.
  do itt=1,6
    call calculate_density(T0,S0,pres,rho_guess,1,nz,eqn_of_state)
    call calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,1,nz,eqn_of_state)
    do k=1,nz
      T0(k) = T0(k) + (GV%Rlay(k) - rho_guess(k)) / drho_dT(k)
    enddo
  enddo

  pi = 4.0*atan(1.0)
  I_ts = 1.0 / thermocline_scale
  I_md = 1.0 / G%max_depth
  do j=js,je ; do i=is,ie
    SST = 0.5*(T0(k1)+T0(nz)) - 0.9*0.5*(T0(k1)-T0(nz)) * &
                               cos(pi*(G%geoLatT(i,j)-G%south_lat)/(G%len_lat))

    do k=1,nz ; e_pert(K) = 0.0 ; enddo

!  The remainder of this subroutine should not be changed.           !

!    This sets the initial thickness (in m) of the layers.  The      !
!  thicknesses are set to insure that: 1.  each layer is at least    !
!  Gv%Angstrom_z thick, and 2.  the interfaces are where they should be    !
!  based on the resting depths and interface height perturbations,   !
!  as long at this doesn't interfere with 1.                         !
    eta1D(nz+1) = -1.0*G%bathyT(i,j)

    do k=nz,2,-1
      T_int = 0.5*(T0(k) + T0(k-1))
      T_frac = (T_int - T0(nz)) / (SST - T0(nz))
      ! Find the z such that T_frac = a exp(z/thermocline_scale) + (1-a) (z+D)/D
      z = 0.0
      do itt=1,6
        err = a_exp * exp(z*I_ts) + (1.0 - a_exp) * (z*I_md + 1.0) - T_frac
        derr_dz = a_exp * I_ts * exp(z*I_ts) + (1.0 - a_exp) * I_md
        z = z - err / derr_dz
      enddo
      e0(K) = z
!       e0(K) = -ML_depth + thermocline_scale * log((T_int - T0(nz)) / (SST - T0(nz)))

      eta1D(K) = e0(K) + e_pert(K)

      if (eta1D(K) > -ML_depth) eta1D(K) = -ML_depth

      if (eta1D(K) < eta1D(K+1) + GV%Angstrom_z) &
        eta1D(K) = eta1D(K+1) + GV%Angstrom_z

      h(i,j,k) = max(eta1D(K) - eta1D(K+1), GV%Angstrom_z)
    enddo
    h(i,j,1) = max(0.0 - eta1D(2), GV%Angstrom_z)

  enddo ; enddo

end subroutine Neverland_initialize_thickness
! -----------------------------------------------------------------------------

!! \class Neverland_initialization
!!
!! The module configures the model for the Neverland experiment.
end module Neverland_initialization
