module box_initialization

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

public box_initialize_topography

! This include declares and sets the variable "version".
#include "version_variable.h"

contains


! -----------------------------------------------------------------------------
!> This subroutine sets up the box test case topography.
subroutine box_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),             intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),              intent(in)  :: param_file !< Parameter file structure
  real,                               intent(in)  :: max_depth  !< Maximum depth of model in m

! This subroutine sets up the box test case topography
  real :: PI = 4.0*atan(1.0)   ! 3.1415926... calculated as 4*atan(1)
  real :: latext, lonext       ! latitude extent of the model area
  real :: x, dx                ! non-dimensional longitudinal grid scale

  character(len=40)  :: mod = "box_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  latext = G%len_lat
  lonext = G%len_lon
  D = 1.0
  dx = (G%geoLonT(is+1,js)-G%geoLonT(is,js))/lonext

  call MOM_mesg("  box_initialization.F90, box_initialize_topography: setting topography", 5)

  call log_version(param_file, mod, version, "")
  

  !  Calculate the depth of the bottom.
  do j=js,je                ! meridional grid points
  do i=is,ie                ! zonal grid points
    x=(G%geoLonT(i,j)-G%west_lon) / lonext      ! non-dimensional longitude
    ! meridional walls at western boundary
    if (x < dx/1.5) then
      D(i,j) = 0.0
    endif

    D(i,j) = D(i,j) * max_depth
  enddo
  enddo

end subroutine box_initialize_topography


end module box_initialization
