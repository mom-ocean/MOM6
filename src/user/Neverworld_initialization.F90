!> Initialization for the "Neverworld" configuration
module Neverworld_initialization

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

public Neverworld_initialize_topography
public Neverworld_initialize_thickness

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> This subroutine sets up the Neverworld test case topography.
subroutine Neverworld_initialize_topography(D, G, param_file, max_depth)
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
  character(len=40)  :: mdl = "Neverworld_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  real :: nl_roughness_amp, nl_top_amp
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call MOM_mesg("  Neverworld_initialization.F90, Neverworld_initialize_topography: setting topography", 5)

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

end subroutine Neverworld_initialize_topography

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

!> Returns the value of a triangular function centered at x=x0 with value 1
!! and linearly decreasing to 0 at x=x0+/-L, and 0 otherwise.
!! If clip is present the top of the cone is cut off at "clip", which
!! effectively defaults to 1.
real function cone(x, x0, L, clip)
  real,           intent(in) :: x    !< non-dimensional coordinate [nondim]
  real,           intent(in) :: x0   !< position of peak [nondim]
  real,           intent(in) :: L    !< half-width of base of cone [nondim]
  real, optional, intent(in) :: clip !< clipping height of cone [nondim]

  cone = max( 0., 1. - abs(x - x0) / L )
  if (present(clip)) cone = min(clip, cone)
end function cone

!> Returns an s-curve s(x) s.t. s(x0)<=0, s(x0+L)>=1 and cubic in between.
real function scurve(x, x0, L)
  real, intent(in) :: x       !< non-dimensional coordinate [nondim]
  real, intent(in) :: x0      !< position of peak [nondim]
  real, intent(in) :: L       !< half-width of base of cone [nondim]
  real :: s

  s = max( 0., min( 1.,( x - x0 ) / L ) )
  scurve = ( 3. - 2.*s ) * ( s * s )
end function scurve

!> Returns a "coastal" profile.
real function cstprof(x, x0, L, lf, bf, sf, sh)
  real, intent(in) :: x       !< non-dimensional coordinate [nondim]
  real, intent(in) :: x0      !< position of peak [nondim]
  real, intent(in) :: L       !< width of profile [nondim]
  real, intent(in) :: lf      !< fraction of width that is "land" [nondim]
  real, intent(in) :: bf      !< fraction of width that is "beach" [nondim]
  real, intent(in) :: sf      !< fraction of width that is "continental slope" [nondim]
  real, intent(in) :: sh      !< depth of shelf as fraction of full depth [nondim]
  real :: s

  s = max( 0., min( 1.,( x - x0 ) / L ) )
  cstprof = sh * scurve(s-lf,0.,bf) + (1.-sh) * scurve(s - (1.-sf),0.,sf)
end function cstprof

!> Distance between points x,y and a line segment (x0,y0) and (x0,y1).
real function dist_line_fixed_x(x, y, x0, y0, y1)
  real, intent(in) :: x       !< non-dimensional x-coordinate [nondim]
  real, intent(in) :: y       !< non-dimensional y-coordinate [nondim]
  real, intent(in) :: x0      !< x-position of line segment [nondim]
  real, intent(in) :: y0      !< y-position of line segment end[nondim]
  real, intent(in) :: y1      !< y-position of line segment end[nondim]
  real :: dx, yr, dy

  dx = x - x0
  yr = min( max(y0,y1), max( min(y0,y1), y ) ) ! bound y by y0,y1
  dy = y - yr ! =0 within y0<y<y1, =y0-y for y<y0, =y-y1 for y>y1
  dist_line_fixed_x = sqrt( dx*dx + dy*dy )
end function dist_line_fixed_x

!> Distance between points x,y and a line segment (x0,y0) and (x1,y0).
real function dist_line_fixed_y(x, y, x0, x1, y0)
  real, intent(in) :: x       !< non-dimensional x-coordinate [nondim]
  real, intent(in) :: y       !< non-dimensional y-coordinate [nondim]
  real, intent(in) :: x0      !< x-position of line segment end[nondim]
  real, intent(in) :: x1      !< x-position of line segment end[nondim]
  real, intent(in) :: y0      !< y-position of line segment [nondim]
  real :: dx, yr, dy

  dist_line_fixed_y = dist_line_fixed_x(y, x, y0, x0, x1)
end function dist_line_fixed_y

!> A "coast profile" applied in an N-S line from lon0,lat0 to lon0,lat1.
real function NS_coast(lon, lat, lon0, lat0, lat1, dlon, sh)
  real, intent(in) :: lon     !< Longitude [degrees_E]
  real, intent(in) :: lat     !< Latitude [degrees_N]
  real, intent(in) :: lon0    !< Longitude of coast [degrees_E]
  real, intent(in) :: lat0    !< Latitude of coast end [degrees_N]
  real, intent(in) :: lat1    !< Latitude of coast end [degrees_N]
  real, intent(in) :: dlon    !< "Radius" of coast profile [degrees]
  real, intent(in) :: sh      !< depth of shelf as fraction of full depth [nondim]
  real :: r

  r = dist_line_fixed_x( lon, lat, lon0, lat0, lat1 )
  NS_coast = cstprof(r, 0., dlon, 0.125, 0.125, 0.5, sh)
end function NS_coast

!> A "coast profile" applied in an E-W line from lon0,lat0 to lon1,lat0.
real function EW_coast(lon, lat, lon0, lon1, lat0, dlat, sh)
  real, intent(in) :: lon     !< Longitude [degrees_E]
  real, intent(in) :: lat     !< Latitude [degrees_N]
  real, intent(in) :: lon0    !< Longitude of coast end [degrees_E]
  real, intent(in) :: lon1    !< Longitude of coast end [degrees_E]
  real, intent(in) :: lat0    !< Latitude of coast [degrees_N]
  real, intent(in) :: dlat    !< "Radius" of coast profile [degrees]
  real, intent(in) :: sh      !< depth of shelf as fraction of full depth [nondim]
  real :: r

  r = dist_line_fixed_y( lon, lat, lon0, lon1, lat0 )
  EW_coast = cstprof(r, 0., dlat, 0.125, 0.125, 0.5, sh)
end function EW_coast

!> A NS ridge
real function NS_ridge(lon, lat, lon0, lat0, lat1, dlon, rh)
  real, intent(in) :: lon     !< Longitude [degrees_E]
  real, intent(in) :: lat     !< Latitude [degrees_N]
  real, intent(in) :: lon0    !< Longitude of ridge center [degrees_E]
  real, intent(in) :: lat0    !< Latitude of ridge end [degrees_N]
  real, intent(in) :: lat1    !< Latitude of ridge end [degrees_N]
  real, intent(in) :: dlon    !< "Radius" of ridge profile [degrees]
  real, intent(in) :: rh      !< depth of ridge as fraction of full depth [nondim]
  real :: r

  r = dist_line_fixed_x( lon, lat, lon0, lat0, lat1 )
  NS_ridge = 1. - rh * cone(r, 0., dlon)
end function NS_ridge


!> A circular ridge
real function circ_ridge(lon, lat, lon0, lat0, ring_radius, ring_thickness, ridge_height)
  real, intent(in) :: lon            !< Longitude [degrees_E]
  real, intent(in) :: lat            !< Latitude [degrees_N]
  real, intent(in) :: lon0           !< Longitude of center of ring [degrees_E]
  real, intent(in) :: lat0           !< Latitude of center of ring [degrees_N]
  real, intent(in) :: ring_radius    !< Radius of ring [degrees]
  real, intent(in) :: ring_thickness !< Radial thickness of ring [degrees]
  real, intent(in) :: ridge_height   !< Ridge height as fraction of full depth [nondim]
  real :: r

  r = sqrt( (lon - lon0)**2 + (lat - lat0)**2 ) ! Pseudo-distance from a point
  r = abs( r - ring_radius) ! Pseudo-distance from a circle
  r = cone(r, 0., ring_thickness, ridge_height) ! 0 .. frac_ridge_height
  circ_ridge = 1. - r ! nondim depths (1-frac_ridge_height) .. 1
end function circ_ridge

!> This subroutine initializes layer thicknesses for the Neverworld test case,
!! by finding the depths of interfaces in a specified latitude-dependent
!! temperature profile with an exponentially decaying thermocline on top of a
!! linear stratification.
subroutine Neverworld_initialize_thickness(h, G, GV, US, param_file, eqn_of_state, P_ref)
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
  character(len=40)  :: mdl = "Neverworld_initialize_thickness" ! This subroutine's name.
  integer :: i, j, k, k1, is, ie, js, je, nz, itt

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call MOM_mesg("  Neverworld_initialization.F90, Neverworld_initialize_thickness: setting thickness", 5)
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

end subroutine Neverworld_initialize_thickness

end module Neverworld_initialization
