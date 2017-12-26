module channel2_initialization

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

public channel2_initialize_topography
public channel2_initialize_sponges

! This include declares and sets the variable "version".
#include "version_variable.h"

contains


! -----------------------------------------------------------------------------
!> This subroutine sets up the channel2 test case topography.
!> channel2 is similar to channel but:
!> 1) with sloped side walls to mimic continental slope, and to reduce numerical instability
!> 2) the sponge layer has no slope
!> 3) the slope on west/east boundaries decay abruptly to zero at the edge of sponge

subroutine channel2_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),             intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),              intent(in)  :: param_file !< Parameter file structure
  real,                               intent(in)  :: max_depth  !< Maximum depth of model in m

  real :: PI = 4.0*atan(1.0)   ! 3.1415926... calculated as 4*atan(1)
  real :: latext, lonext       ! latitude extent of the model area
  real :: ep = epsilon(1.)     ! an infinitesimally small quantity
  real :: x, y, sa_dim = 1500.0! dimensional height of Scotia Arc
  real :: sa                   ! the non-dimensional height of Scotia Arc top;
                               ! default value is 1500/4000=0.375
  real :: dx                   ! non-dimensional longitudinal grid scale
  real :: reentrants, reentrantn  ! the non-dimensional latitudes of the southern and northern
                                  ! boundary of the reentrant channel
  real :: sdp = 6.0, ll=6.0, ssp = 2.0    ! the width of the slope in Drake Passage, the half width of continental slope, and sponge width

  character(len=40)  :: mod = "channel2_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  sa = sa_dim / max_depth
  latext = G%len_lat
  lonext = G%len_lon
  reentrants = 6.0/latext          ! non-dimensional southern bound of the reentrant zone
  reentrantn = 10.0/latext         ! non-dimensional northern bound of the reentrant zone
  sdp = sdp/latext
  ssp = ssp/latext
  D = 0.0
  dx = (G%geoLonT(is+1,js)-G%geoLonT(is,js))/lonext

  call MOM_mesg("  channel2_initialization.F90, channel2_initialize_topography: setting topography", 5)

  call log_version(param_file, mod, version, "")


  !  Calculate the depth of the bottom.
  do j = js,je                ! meridional grid points
  do i = is,ie                ! zonal grid points
    x = (G%geoLonT(i,j)-G%west_lon) / lonext      ! non-dimensional longitude
    y=(G%geoLatT(i,j)-G%south_lat) / latext     ! non-dimensional latitude

    D(i,j) = 1.0 - spike(x-dx/2, ll/lonext)*spike(min(0.0, y-reentrantn-sdp/2.0), sdp) &                 ! Patagonia, west
                -spike(x-1.0+dx/2, ll/lonext)*spike(min(0.0, y-reentrantn-sdp/2.0), sdp) &               ! Patanogina, east
                -spike(x-dx/2, ll/lonext)*spike(max(0.0, y-reentrants+sdp/2.0), sdp) &                   ! Antarctic Peninsula, west
                -spike(x-1.0+dx/2, ll/lonext)*spike(max(0.0, y-reentrants+sdp/2.0), sdp) &               ! Antarctic Peninsula, east
                -spike(y, ll/latext) &                                                                   ! Antarctica
                - sa *cosbell(x-20.0/lonext, 2.5/lonext) * homo(y-8.0/latext, 2.0/latext) &              !Scotia Arc East, center
                - sa * cosbell(x-20.0/lonext, 2.5/lonext) * cosbellh(y-10.0/latext-ep, 3.0/latext, 1.) & !Scotia Arc East, north slope
                - sa * cosbell(x-20.0/lonext, 2.5/lonext) * cosbellh(y-6.0/latext+ep, 3.0/latext, -1.) & !Scotia Arc East, south slope
                - sa * cosbell(y-12.0/latext, 2.5/latext) * cosbellh(x-18.0/lonext-ep, 2.5/lonext, 1.) & !Scotia Arc North, east half (slope side)
                - sa * cosbell(y-12.0/latext, 2.5/latext) * homo(x-9.0/lonext, 9.0/lonext) &             !Scotia Arc North, west half
                - sa * cosbell(y-4.0/latext, 2.5/latext) * cosbellh(x-18.0/lonext-ep, 2.5/lonext, 1.) &  !Scotia Arc South, east half (slope side)
                - sa * cosbell(y-4.0/latext, 2.5/latext) * homo(x-9.0/lonext, 9.0/lonext)                !Scotia Arc South, west half

      ! make sure no deeper than max depth and no shallower than Scotia Arc top IN the ocean interior
      if (D(i,j)<1.0 - sa .and. x>=dx/1.5+ll/2/lonext .and. x<=1.0-ll/2/lonext-dx/1.5 .and. y>=ll/2/latext .and. y<=1.0-ll/2/latext) then
        D(i,j) = 1 - sa
      else if (D(i,j) > 1.0) then
        D(i,j) = 1.0
      endif
      
      ! no continental slope in the sponge layer
      if (y >= 1.0-ssp) then
        D(i,j)=1.0
      endif

      ! make sure the model is not zonally reentrant outside of Drake Passage
      if (((y>=reentrantn+sdp/2 .or. y<=reentrants-sdp/2) .and. x<dx/1.5) &
                .or. ((y>=reentrantn+sdp/2 .or. y<=reentrants-sdp/2) .and. x>1.0-dx/1.5)) then 
        D(i,j) = 0.0
      endif

    D(i,j) = D(i,j) * max_depth
  enddo
  enddo


end subroutine channel2_initialize_topography



! -----------------------------------------------------------------------------
!> Sets up the the inverse restoration time (Idamp), and
! the values towards which the interface heights and an arbitrary
! number of tracers should be restored within each sponge.
subroutine channel2_initialize_sponges(G, GV, use_temperature, tv, param_file, CSp, h)
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
  character(len=40)  :: mdl = "channel2_initialize_sponges" ! This subroutine's name

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  eta(:,:,:) = 0.0 ; Idamp(:,:) = 0.0; eta0(:) = 0.0; 
  dx = (G%geoLonT(is+1,js)-G%geoLonT(is,js))/G%len_lon

  ! target interface heights: all negative values
  eta0 = (/0.0,0.0,0.0,-49.0,-99.0,-152.0,-211.0,-281.0,-361.0,-450.0, &
                -550.0,-659.0,-779.0,-908.0,-1047.0,-1197.0,-1357.0, &
                -1527.0,-1707.0,-1898.0,-2095.0,-2306.0,-2529.0,-2760.0, &
                -2998.0,-3244.0,-3492.0,-3750.0,-4000.0,-4000.0,-4000.0 /)

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


end subroutine channel2_initialize_sponges


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

end module channel2_initialization
