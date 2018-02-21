module channel4t_initialization

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

public channel4t_initialize_sponges

! This include declares and sets the variable "version".
#include "version_variable.h"

contains


! -----------------------------------------------------------------------------
!> Sets up the the inverse restoration time (Idamp), and
! the values towards which the interface heights and an arbitrary
! number of tracers should be restored within each sponge.
!
! difference from channel 4: use stratification at 32S to test the rebustmness
! of the sponge
! e.g. only eta0 is different

subroutine channel4t_initialize_sponges(G, GV, use_temperature, tv, param_file, CSp, h)
  type(ocean_grid_type), intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV  !< The ocean's vertical grid structure
                                             ! so as to be used as input of set_coord_from_file
  logical, intent(in) :: use_temperature    !< Switch for temperature.
  type(thermo_var_ptrs), intent(in) :: tv   !< A structure containing pointers
                                            !! to any available thermodynamic
                                            !! fields, potential temperature and
                                            !! salinity or mixed layer density.
                                            !! Absent fields have NULL ptrs.
  type(param_file_type), intent(in) :: param_file !< A structure indicating the
                                            !! open file to parse for model
                                            !! parameter values.
  type(sponge_CS),   pointer    :: CSp      !< A pointer that is set to point to
                                            !! the control structure for the
                                            !! sponge module.
  real, intent(in), dimension(SZI_(G),SZJ_(G), SZK_(G)) :: h !< Thickness field. - may not be used!
  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1) ! A temporary array for eta, m.
                                         ! eta only varies in z
  real :: Idamp(SZI_(G),SZJ_(G))         ! The inverse damping rate, in s-1.
                                         ! Idamp only varies in y; /= 0 only within the sponge layer
  real, dimension(SZK_(G)+1) :: eta0    ! target interface heights for density class in the sponge layer
  real :: damp_rate, damp, spongelen = 2.0,  min_depth, nlat, dx
  ! spongelen: thickness of sponge layer in dimensional degree
  ! dx is non-dimensional longitudinal grid increment
  integer :: i, j, k, is, ie, js, je, nz
  logical, save :: first_call = .true.
  character(len=40)  :: mdl = "channel4t_initialize_sponges" ! This subroutine's name

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  eta(:,:,:) = 0.0 ; Idamp(:,:) = 0.0; eta0(:) = 0.0; 
  dx = (G%geoLonT(is+1,js)-G%geoLonT(is,js))/G%len_lon

  ! target interface heights: all negative values
  ! corresponds to rho 8
  
    eta0 = (/0., 0., 0., 0., 0., &
          0.0, -27.2, -74.5, -132.9, -242.4, &
          -352.9, -448.5, -532.5, -610.2, -684.6, &
          -758.3, -834.6, -916.9, -1009.7, -1124.0, &
          -1279.0, -1516.3, -1811.5, -2098.0, -2389.3, &
          -2680.5, -2958.3, -3239.6, -3531.8, -3844.1, -4000. /)

  if (first_call) call log_version(param_file, mdl, version)
  first_call = .false.

  call get_param(param_file, mdl, "SPONGE_RATE", damp_rate, &
                 "The rate at which the zonal-mean sponges damp.", units="s-1", &
                 default = 1.0/(10.0*86400.0))
  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)


  nlat = G%south_lat + G%len_lat        ! should be -30.0 degree

  
  ! initialize the damping rate so it is 0 outside of the sponge layer &
  ! and increases linearly with latitude within the sponge layer
  do i = is, ie; do j = js, je
    if (G%geoLatT(i,j) >= nlat-spongelen) then
      damp = damp_rate/spongelen * (G%geoLatT(i,j)-nlat+spongelen)
    else
      damp = 0.0   ! outside of the sponge
    endif
    
    do k = 1,nz+1; eta(i,j,k) = eta0(k); enddo    ! initialize target heights for each (lat,lon) grid

    if (G%bathyT(i,j) > min_depth) then ! bathT: Ocean bottom depth (positive) at tracer points, in m.
      Idamp(i,j) = damp                 ! no need to divide this by 86400!!!
    else
      Idamp(i,j) = 0.0
    endif     ! so that at the side walls, Idamp = 0
  enddo; enddo

  call initialize_sponge(Idamp, eta, G, param_file, CSp)

!From MOM_sponge.F90
!subroutine initialize_sponge(Idamp, eta, G, param_file, CSp, &
!                             Iresttime_i_mean, int_height_i_mean)
!
! Arguments: Idamp - The inverse of the restoring time, in s-1.
!  (in)      eta - The interface heights to damp back toward, in m.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CSp - A pointer that is set to point to the control structure
!                 for this module


end subroutine channel4t_initialize_sponges


! -----------------------------------------------------------------------------
! define functions used in the above subroutines


!> Returns the value of a sinusoidal bell function 
  real function spike(x,L)

    real, intent(in) :: x
    real, intent(in) :: L
    real             :: PI = 4.0*atan(1.0)

    spike = 1-sin(PI*min(abs(x)/L, 0.5))

  end function spike

!> Returns the value of a cosine-bell function evaluated at x/L
 real function cosbell(x,L)

   real , intent(in) :: x                         !< non-dimensional position
   real , intent(in) :: L                         !< non-dimensional width
   real              :: PI = 4.0 * atan(1.0)      !< 3.1415926... calculated as 4*atan(1)

   cosbell = 0.5 * (1 + cos(PI*MIN(ABS(x/L),1.0)))
 end function cosbell

!< Return the value of a half-cosine-bell function evaluated at x/L;
!< i.e. from peak to trough only on one side of the bell
 real function cosbellh(x, L, dir)

  real, intent(in) :: x       !< non-dimensional position
  real, intent(in) :: L       !< non-dimensional width
  real             :: PI, xx  !< 3.1415926... calculated as 4*atan(1)
  real, intent(in) :: dir     !< direction flag; 1 for east/north; -1 for west/south
  PI        = 4.0*atan(1.0)

    !< if the grid falls on the opposite side of the bell, override x to be so big that x/L > 1
    if (x*dir .lt. 0.0) then
      xx = L+1
    else
      xx = x
    endif

    cosbellh  = cos(PI/2.0*MIN(abs(xx)/L, 1.0))
  end function cosbellh

  !< similar to cosbellh but takes a different shape of bell
    real function cosbellhnew(x, L, dir)

      real, intent(in) :: x       !< non-dimensional position
      real, intent(in) :: L       !< non-dimensional width
      real             :: PI, xx  !< 3.1415926... calculated as 4*atan(1)
      real, intent(in) :: dir     !< direction flag; 1 for east/north; -1 for west/south
      PI        = 4.0*atan(1.0)

      !< if the grid falls on the opposite side of the bell, override x to be so big that x/L > 1
      if (x*dir .lt. 0.0) then
        xx = L+1
      else
        xx = x
      endif

     cosbellhnew  = 0.5*(1+cos(PI*MIN(xx/L, 1.0)))
    end function cosbellhnew


 !< make sure the depth within L is homogeneous
   real function homo(x, L)

     real, intent(in) :: x       !< non-dimensional position
     real, intent(in) :: L       !< non-dimensional width

     !< if x falls within -L ~ L, assign 1 to the non-dimensional depth
    if (abs(x) .le. L) then
      homo = 1.0
    else
      homo = 0.0
    endif
   end function homo

end module channel4t_initialization
