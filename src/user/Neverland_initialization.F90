!> Initialization for the "Neverland" configuration
module Neverland_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type

use random_numbers_mod, only: initializeRandomNumberStream, getRandomNumbers, randomNumberStream

implicit none ; private

#include <MOM_memory.h>

public Neverland_initialize_topography
public Neverland_initialize_thickness

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> This subroutine sets up the Neverland test case topography.
subroutine Neverland_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),  intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                           intent(out) :: D !< Ocean bottom depth in the units of depth_max
  type(param_file_type),   intent(in)  :: param_file !< Parameter file structure
  real,                    intent(in)  :: max_depth !< Maximum ocean depth in arbitrary units

  ! Local variables
  real :: PI                   ! 3.1415926... calculated as 4*atan(1)
  real :: D0                   ! A constant to make the maximum     !
                               ! basin depth MAXIMUM_DEPTH.         !
  real :: x, y
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "Neverland_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  real :: nl_roughness_amp, nl_top_amp
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call MOM_mesg("  Neverland_initialization.F90, Neverland_initialize_topography: setting topography", 5)

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "NL_ROUGHNESS_AMP", nl_roughness_amp, &
                 "Amplitude of wavy signal in bathymetry.", default=0.05)
  call get_param(param_file, mdl, "NL_CONTINENT_AMP", nl_top_amp, &
                 "Scale factor for topography - 0.0 for no continents.", default=1.0)

  PI = 4.0*atan(1.0)

!  Calculate the depth of the bottom.
  do j=js,je ; do i=is,ie
    x = (G%geoLonT(i,j)-G%west_lon) / G%len_lon
    y =( G%geoLatT(i,j)-G%south_lat) / G%len_lat
!  This sets topography that has a reentrant channel to the south.
    D(i,j) = 1.0 - 1.1 * spike(y-1,0.12) - 1.1 * spike(y,0.12) - & !< The great northern wall and Antarctica
              nl_top_amp*( &
                (1.2 * spike(x,0.2) + 1.2 * spike(x-1.0,0.2)) * spike(MIN(0.0,y-0.3),0.2) & !< South America
              +  1.2 * spike(x-0.5,0.2) * spike(MIN(0.0,y-0.55),0.2)       & !< Africa
              +  1.2 * (spike(x,0.12)  + spike(x-1,0.12)) * spike(MAX(0.0,y-0.06),0.12)    & !< Antarctic Peninsula
              +  0.1 * (cosbell(x,0.1) + cosbell(x-1,0.1))                 & !< Drake Passage ridge
              +  0.5 * cosbell(x-0.16,0.05) * (cosbell(y-0.18,0.13)**0.4)  & !< Scotia Arc East
              +  0.4 * (cosbell(x-0.09,0.08)**0.4) * cosbell(y-0.26,0.05)  & !< Scotia Arc North
              +  0.4 * (cosbell(x-0.08,0.08)**0.4) * cosbell(y-0.1,0.05))   & !< Scotia Arc South
              -  nl_roughness_amp * cos(14*PI*x) * sin(14*PI*y)            & !< roughness
              -  nl_roughness_amp * cos(20*PI*x) * cos(20*PI*y)              !< roughness
    if (D(i,j) < 0.0) D(i,j) = 0.0
    D(i,j) = D(i,j) * max_depth
  enddo ; enddo

end subroutine Neverland_initialize_topography
! -----------------------------------------------------------------------------

!> Returns the value of a cosine-bell function evaluated at x/L
real function cosbell(x, L)
  real , intent(in) :: x       !< non-dimensional position
  real , intent(in) :: L       !< non-dimensional width
  real              :: PI      !< 3.1415926... calculated as 4*atan(1)

  PI      = 4.0*atan(1.0)
  cosbell = 0.5 * (1 + cos(PI*MIN(ABS(x/L),1.0)))
end function cosbell

!> Returns the value of a sin-spike function evaluated at x/L
real function spike(x, L)

  real , intent(in) :: x       !< non-dimensional position
  real , intent(in) :: L       !< non-dimensional width
  real              :: PI      !< 3.1415926... calculated as 4*atan(1)

  PI    = 4.0*atan(1.0)
  spike = (1 - sin(PI*MIN(ABS(x/L),0.5)))
end function spike

!> This subroutine initializes layer thicknesses for the Neverland test case,
!! by finding the depths of interfaces in a specified latitude-dependent
!! temperature profile with an exponentially decaying thermocline on top of a
!! linear stratification.
subroutine Neverland_initialize_thickness(h, G, GV, US, param_file, eqn_of_state, P_ref)
  type(ocean_grid_type),   intent(in) :: G                    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV                   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in) :: US                   !< A dimensional unit scaling type
  real, intent(out), dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h !< The thickness that is being
                                                              !! initialized [H ~> m or kg m-2].
  type(param_file_type),   intent(in) :: param_file           !< A structure indicating the open
                                                              !! file to parse for model
                                                              !! parameter values.
  type(EOS_type),          pointer    :: eqn_of_state         !< integer that selects the
                                                              !! equation of state.
  real,                    intent(in) :: P_Ref                !< The coordinate-density
                                                              !! reference pressure [R L2 T-2 ~> Pa].
  ! Local variables
  real :: e0(SZK_(G)+1)     ! The resting interface heights, in depth units [Z ~> m],
                            ! usually negative because it is positive upward.
  real, dimension(SZK_(G)) :: h_profile ! Vector of initial thickness profile [Z ~> m].
  real :: e_interface ! Current interface position [Z ~> m].
  real :: x,y,r1,r2 ! x,y and radial coordinates for computation of initial pert.
  real :: pert_amp ! Amplitude of perturbations measured in Angstrom_H
  real :: h_noise ! Amplitude of noise to scale h by
  real :: noise ! Noise
  type(randomNumberStream) :: rns ! Random numbers for stochastic tidal parameterization
  character(len=40)  :: mdl = "Neverland_initialize_thickness" ! This subroutine's name.
  integer :: i, j, k, k1, is, ie, js, je, nz, itt

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call MOM_mesg("  Neverland_initialization.F90, Neverland_initialize_thickness: setting thickness", 5)
  call get_param(param_file, mdl, "INIT_THICKNESS_PROFILE", h_profile, &
                 "Profile of initial layer thicknesses.", units="m", scale=US%m_to_Z, &
                 fail_if_missing=.true.)
  call get_param(param_file, mdl, "NL_THICKNESS_PERT_AMP", pert_amp, &
                 "Amplitude of finite scale perturbations as fraction of depth.", &
                 units="nondim", default=0.)
  call get_param(param_file, mdl, "NL_THICKNESS_NOISE_AMP", h_noise, &
                 "Amplitude of noise to scale layer by.", units="nondim", default=0.)

  ! e0 is the notional position of interfaces
  e0(1) = 0. ! The surface
  do k=1,nz
    e0(k+1) = e0(k) - h_profile(k)
  enddo

  do j=js,je ; do i=is,ie
    e_interface = -G%bathyT(i,j)
    do k=nz,2,-1
      h(i,j,k) = GV%Z_to_H * (e0(k) - e_interface) ! Nominal thickness
      x=(G%geoLonT(i,j)-G%west_lon)/G%len_lon
      y=(G%geoLatT(i,j)-G%south_lat)/G%len_lat
      r1=sqrt((x-0.7)**2+(y-0.2)**2)
      r2=sqrt((x-0.3)**2+(y-0.25)**2)
      h(i,j,k) = h(i,j,k) + pert_amp * (e0(k) - e0(nz+1)) * GV%Z_to_H * &
                            (spike(r1,0.15)-spike(r2,0.15)) ! Prescribed perturbation
      if (h_noise /= 0.) then
        rns = initializeRandomNumberStream( int( 4096*(x + (y+1.)) ) )
        call getRandomNumbers(rns, noise) ! x will be in (0,1)
        noise = h_noise * 2. * ( noise - 0.5 ) ! range -h_noise to h_noise
        h(i,j,k) = ( 1. + noise ) * h(i,j,k)
      endif
      h(i,j,k) = max( GV%Angstrom_H, h(i,j,k) ) ! Limit to non-negative
      e_interface = e_interface + GV%H_to_Z * h(i,j,k) ! Actual position of upper interface
    enddo
    h(i,j,1) = GV%Z_to_H * (e0(1) - e_interface) ! Nominal thickness
    h(i,j,1) = max( GV%Angstrom_H, h(i,j,1) ) ! Limit to non-negative
  enddo ; enddo

end subroutine Neverland_initialize_thickness

end module Neverland_initialization
