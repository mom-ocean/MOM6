!> An idealized topography building system
module basin_builder

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_string_functions, only : lowercase
use MOM_unit_scaling, only : unit_scale_type

implicit none ; private

#include <MOM_memory.h>

public basin_builder_topography

! This include declares and sets the variable "version".
# include "version_variable.h"
character(len=40) :: mdl = "basin_builder" !< This module's name.

contains

!> Constructs idealized topography from simple functions
subroutine basin_builder_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),  intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                           intent(out) :: D !< Ocean bottom depth in the units of depth_max
  type(param_file_type),   intent(in)  :: param_file !< Parameter file structure
  real,                    intent(in)  :: max_depth !< Maximum ocean depth in arbitrary units
  ! Local variables
  character(len=17) :: pname1, pname2 ! For construction of parameter names
  character(len=20) :: funcs ! Basin build function
  real, dimension(20) :: pars ! Parameters for each function
  real :: lon ! Longitude [degrees_E}
  real :: lat ! Latitude [degrees_N]
  integer :: i, j, n, n_funcs

  call MOM_mesg("  basin_builder.F90, basin_builder_topography: setting topography", 5)
  call log_version(param_file, mdl, version, "")

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    D(i,j) = 1.0
  enddo ; enddo

  call get_param(param_file, mdl, "BBUILDER_N", n_funcs, &
                 "Number of pieces of topography to use.", fail_if_missing=.true.)

  do n=1,n_funcs
    write( pname1, "('BBUILDER_',i3.3,'_FUNC')" ) n
    write( pname2, "('BBUILDER_',i3.3,'_PARS')" ) n
    call get_param(param_file, mdl, pname1, funcs, &
                   "The basin builder function to apply with parameters "//&
                   trim(pname2)//". Choices are: NS_COAST, EW_COAST, "//&
                   "CIRC_CONIC_RIDGE, NS_CONIC_RIDGE, CIRC_SCURVE_RIDGE, "//&
                   "NS_SCURVE_RIDGE.", &
                   fail_if_missing=.true.)
    pars(:) = 0.
    if (trim(lowercase(funcs)) == 'ns_coast') then
      call get_param(param_file, mdl, pname2, pars(1:5), &
                     "NS_COAST parameters: longitude, starting latitude, "//&
                     "ending latitude, footprint radius, shelf depth.", &
                     units="degrees_E,degrees_N,degrees_N,degrees,m", &
                     fail_if_missing=.true.)
      pars(5) = pars(5) / max_depth
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        lon = G%geoLonT(i,j)
        lat = G%geoLatT(i,j)
        D(i,j) = min( D(i,j), NS_coast(lon, lat, pars(1), pars(2), pars(3), pars(4), pars(5)) )
      enddo ; enddo
    elseif (trim(lowercase(funcs)) == 'ns_conic_ridge') then
      call get_param(param_file, mdl, pname2, pars(1:5), &
                     "NS_CONIC_RIDGE parameters: longitude, starting latitude, "//&
                     "ending latitude, footprint radius, ridge height.", &
                     units="degrees_E,degrees_N,degrees_N,degrees,m", &
                     fail_if_missing=.true.)
      pars(5) = pars(5) / max_depth
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        lon = G%geoLonT(i,j)
        lat = G%geoLatT(i,j)
        D(i,j) = min( D(i,j), NS_conic_ridge(lon, lat, pars(1), pars(2), pars(3), pars(4), pars(5)) )
      enddo ; enddo
    elseif (trim(lowercase(funcs)) == 'ns_scurve_ridge') then
      call get_param(param_file, mdl, pname2, pars(1:5), &
                     "NS_SCURVE_RIDGE parameters: longitude, starting latitude, "//&
                     "ending latitude, footprint radius, ridge height.", &
                     units="degrees_E,degrees_N,degrees_N,degrees,m", &
                     fail_if_missing=.true.)
      pars(5) = pars(5) / max_depth
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        lon = G%geoLonT(i,j)
        lat = G%geoLatT(i,j)
        D(i,j) = min( D(i,j), NS_scurve_ridge(lon, lat, pars(1), pars(2), pars(3), pars(4), pars(5)) )
      enddo ; enddo
    elseif (trim(lowercase(funcs)) == 'ew_coast') then
      call get_param(param_file, mdl, pname2, pars(1:5), &
                     "EW_COAST parameters: latitude, starting longitude, "//&
                     "ending longitude, footprint radius, shelf depth.", &
                     units="degrees_N,degrees_E,degrees_E,degrees,m", &
                     fail_if_missing=.true.)
      pars(5) = pars(5) / max_depth
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        lon = G%geoLonT(i,j)
        lat = G%geoLatT(i,j)
        D(i,j) = min( D(i,j), EW_coast(lon, lat, pars(1), pars(2), pars(3), pars(4), pars(5)) )
      enddo ; enddo
    elseif (trim(lowercase(funcs)) == 'circ_conic_ridge') then
      call get_param(param_file, mdl, pname2, pars(1:5), &
                     "CIRC_CONIC_RIDGE parameters: center longitude, center latitude, "//&
                     "ring radius, footprint radius, ridge height.", &
                     units="degrees_E,degrees_N,degrees,degrees,m", &
                     fail_if_missing=.true.)
      pars(5) = pars(5) / max_depth
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        lon = G%geoLonT(i,j)
        lat = G%geoLatT(i,j)
        D(i,j) = min( D(i,j), circ_conic_ridge(lon, lat, pars(1), pars(2), pars(3), pars(4), pars(5)) )
      enddo ; enddo
    elseif (trim(lowercase(funcs)) == 'circ_scurve_ridge') then
      call get_param(param_file, mdl, pname2, pars(1:5), &
                     "CIRC_SCURVe_RIDGE parameters: center longitude, center latitude, "//&
                     "ring radius, footprint radius, ridge height.", &
                     units="degrees_E,degrees_N,degrees,degrees,m", &
                     fail_if_missing=.true.)
      pars(5) = pars(5) / max_depth
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        lon = G%geoLonT(i,j)
        lat = G%geoLatT(i,j)
        D(i,j) = min( D(i,j), circ_scurve_ridge(lon, lat, pars(1), pars(2), pars(3), pars(4), pars(5)) )
      enddo ; enddo
    else
      call MOM_error(FATAL, "basin_builder.F90, basin_builer_topography:\n"//&
                     "Unrecognized function "//trim(funcs))
    endif

  enddo ! n

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    ! Dimensionalize by scaling 1 to max_depth
    D(i,j) = D(i,j) * max_depth
  enddo ; enddo

end subroutine basin_builder_topography

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

!> A "coast profile" applied in an N-S line from lonC,lat0 to lonC,lat1.
real function NS_coast(lon, lat, lonC, lat0, lat1, dlon, sh)
  real, intent(in) :: lon     !< Longitude [degrees_E]
  real, intent(in) :: lat     !< Latitude [degrees_N]
  real, intent(in) :: lonC    !< Longitude of coast [degrees_E]
  real, intent(in) :: lat0    !< Latitude of coast end [degrees_N]
  real, intent(in) :: lat1    !< Latitude of coast end [degrees_N]
  real, intent(in) :: dlon    !< "Radius" of coast profile [degrees]
  real, intent(in) :: sh      !< depth of shelf as fraction of full depth [nondim]
  real :: r

  r = dist_line_fixed_x( lon, lat, lonC, lat0, lat1 )
  NS_coast = cstprof(r, 0., dlon, 0.125, 0.125, 0.5, sh)
end function NS_coast

!> A "coast profile" applied in an E-W line from lon0,latC to lon1,latC.
real function EW_coast(lon, lat, latC, lon0, lon1, dlat, sh)
  real, intent(in) :: lon     !< Longitude [degrees_E]
  real, intent(in) :: lat     !< Latitude [degrees_N]
  real, intent(in) :: latC    !< Latitude of coast [degrees_N]
  real, intent(in) :: lon0    !< Longitude of coast end [degrees_E]
  real, intent(in) :: lon1    !< Longitude of coast end [degrees_E]
  real, intent(in) :: dlat    !< "Radius" of coast profile [degrees]
  real, intent(in) :: sh      !< depth of shelf as fraction of full depth [nondim]
  real :: r

  r = dist_line_fixed_y( lon, lat, lon0, lon1, latC )
  EW_coast = cstprof(r, 0., dlat, 0.125, 0.125, 0.5, sh)
end function EW_coast

!> A NS ridge with a cone profile
real function NS_conic_ridge(lon, lat, lonC, lat0, lat1, dlon, rh)
  real, intent(in) :: lon     !< Longitude [degrees_E]
  real, intent(in) :: lat     !< Latitude [degrees_N]
  real, intent(in) :: lonC    !< Longitude of ridge center [degrees_E]
  real, intent(in) :: lat0    !< Latitude of ridge end [degrees_N]
  real, intent(in) :: lat1    !< Latitude of ridge end [degrees_N]
  real, intent(in) :: dlon    !< "Radius" of ridge profile [degrees]
  real, intent(in) :: rh      !< depth of ridge as fraction of full depth [nondim]
  real :: r

  r = dist_line_fixed_x( lon, lat, lonC, lat0, lat1 )
  NS_conic_ridge = 1. - rh * cone(r, 0., dlon)
end function NS_conic_ridge

!> A NS ridge with an scurve profile
real function NS_scurve_ridge(lon, lat, lonC, lat0, lat1, dlon, rh)
  real, intent(in) :: lon     !< Longitude [degrees_E]
  real, intent(in) :: lat     !< Latitude [degrees_N]
  real, intent(in) :: lonC    !< Longitude of ridge center [degrees_E]
  real, intent(in) :: lat0    !< Latitude of ridge end [degrees_N]
  real, intent(in) :: lat1    !< Latitude of ridge end [degrees_N]
  real, intent(in) :: dlon    !< "Radius" of ridge profile [degrees]
  real, intent(in) :: rh      !< depth of ridge as fraction of full depth [nondim]
  real :: r

  r = dist_line_fixed_x( lon, lat, lonC, lat0, lat1 )
  NS_scurve_ridge = 1. - rh * (1. - scurve(r, 0., dlon) )
end function NS_scurve_ridge

!> A circular ridge with cutoff conic profile
real function circ_conic_ridge(lon, lat, lon0, lat0, ring_radius, ring_thickness, ridge_height)
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
  circ_conic_ridge = 1. - r ! nondim depths (1-frac_ridge_height) .. 1
end function circ_conic_ridge

!> A circular ridge with cutoff scurve profile
real function circ_scurve_ridge(lon, lat, lon0, lat0, ring_radius, ring_thickness, ridge_height)
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
  r = 1. - scurve(r, 0., ring_thickness) ! 0 .. 1
  r = r * ridge_height ! 0 .. frac_ridge_height
  circ_scurve_ridge = 1. - r ! nondim depths (1-frac_ridge_height) .. 1
end function circ_scurve_ridge

end module basin_builder
