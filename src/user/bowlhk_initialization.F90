module bowlhk_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_sponge, only : sponge_CS, initialize_sponge, set_up_sponge_field, Apply_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_tracer_registry, only : tracer_registry_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type

implicit none ; private

#include <MOM_memory.h>

public bowlhk_initialize_topography

! This include declares and sets the variable "version".
#include "version_variable.h"

contains


! -----------------------------------------------------------------------------
!> This subroutine sets up hk's bowl test case topography.
subroutine bowlhk_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),             intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),              intent(in)  :: param_file !< Parameter file structure
  real,                               intent(in)  :: max_depth  !< Maximum depth of model in m

  real :: PI = 4.0*atan(1.0)   ! 3.1415926... calculated as 4*atan(1)
  real :: latext, lonext       ! latitude extent of the model area
  real :: ep = epsilon(1.)     ! an infinitesimally small quantity
  real :: x, y, ll=5.0         ! the width of continental slope, in degree
  real :: lx, ly, dx           ! non-dimensional longitudinal grid scale
  character(len=40)  :: mod = "bowlhk_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  latext = G%len_lat; lonext = G%len_lon
  lx = ll/lonext; ly = ll/latext                ! non-dimensional continental shelf width
  D = 0.0
  dx = (G%geoLonT(is+1,js)-G%geoLonT(is,js))/lonext     ! non-dimensional zonal grid increment

  call MOM_mesg("  bowlhk_initialization.F90, bowlhk_initialize_topography: setting topography", 5)

  call log_version(param_file, mod, version, "")


  !  Calculate the depth of the bottom.
  do j=js,je                ! meridional grid points
  do i=is,ie                ! zonal grid points
    x=(G%geoLonT(i,j)-G%west_lon) / lonext      ! non-dimensional longitude
    y=(G%geoLatT(i,j)-G%south_lat) / latext     ! non-dimensional latitude

    ! this sets up a bowl shaped topography, with sloping continental shelf in all 4 sides
    ! named bowlhk to dinstinguish from the bowl topo inherint in MOM6 by GFDL

    if (x<=lx .and. y<=ly) then
      D(i,j) = max(-1.0+cos(x/lx*PI/2), -1.0+cos(y/ly*PI/2))
    elseif (x<=lx .and. y>=1-ly) then
      D(i,j) = max(-1.0+cos(x/lx*PI/2), -1.0+cos((y-1)/ly*PI/2))
    elseif (x>=1-lx .and. y<=ly) then
      D(i,j) = max(-1.0+cos((x-1)/lx*PI/2), -1.0+cos(y/ly*PI/2))
    elseif (x>=1-lx .and. y>=1-ly) then
      D(i,j) = max(-1.0+cos((x-1)/lx*PI/2), -1.0+cos((y-1)/ly*PI/2)) 
    elseif (x<=lx) then
      D(i,j) = -1.0+cos(x/lx*PI/2)
    elseif (x>=1-lx) then
      D(i,j) = -1.0+cos((x-1)/lx*PI/2)
    elseif (y<=ly) then
      D(i,j) = -1.0+cos(y/ly*PI/2)
    elseif (y>=1-ly) then
      D(i,j) = -1.0+cos((y-1)/ly*PI/2)
    else
      D(i,j) = -1.0
    endif

    if (D(i,j) < -1.0) then
      D(i,j) = -1.0
    else if (D(i,j) > 0.0) then
      D(i,j) = 0.0
    endif
    
    ! make sure the region is zonally blocked 
    if (x<=dx/1.5) then
      D(i,j) = 0.0
    elseif (x>=1-dx/1.5) then
      D(i,j) = 0.0
    endif
 
    D(i,j) = -D(i,j) * max_depth
  enddo
  enddo

end subroutine bowlhk_initialize_topography


end module bowlhk_initialization
