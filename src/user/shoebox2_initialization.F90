module shoebox2_initialization

! This file is part of MOM6. See LICENSE.md for the license.

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

public shoebox2_initialize_topography

contains

! -----------------------------------------------------------------------------
!> This subroutine sets up the shoebox2 test case topography.
!> Differences from shoebox: added continental slopes along all 4 side walls

subroutine shoebox2_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),             intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),              intent(in)  :: param_file !< Parameter file structure
  real,                               intent(in)  :: max_depth  !< Maximum depth of model in m

! This subroutine sets up the shoebox2 test case topography
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

 
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "shoebox2_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  sa = sa_dim / max_depth
  latext = G%len_lat
  lonext = G%len_lon
  reentrants = 6.0/latext          ! non-dimensional southern bound of the reentrant zone
  reentrantn = 10.0/latext         ! non-dimensional northern bound of the reentrant zone
  sdp = sdp/latext
  ssp = ssp/latext
  D = 0.0
  dx = (G%geoLonT(is+1,js) - G%geoLonT(is, js)) / lonext

  call MOM_mesg("  shoebox2_initialization.F90, shoebox2_initialize_topography: setting topography", 5)

  call log_version(param_file, mod, version, "")

  !  Calculate the depth of the bottom.
  do j=js,je                ! meridional grid points
  do i=is,ie                ! zonal grid points
    x=(G%geoLonT(i,j)-G%west_lon) / lonext      ! non-dimensional longitude
    y=(G%geoLatT(i,j)-G%south_lat) / latext     ! non-dimensional latitude

    !  This sets topography that has a reentrant channel to the south and a basin in the north

    D(i,j) = 1.0 - spike(x-dx/2, ll/lonext)*spike(min(0.0, y-reentrantn-sdp/2.0), sdp) &                 ! Patagonia, west
                -spike(x-1.0+dx/2, ll/lonext)*spike(min(0.0, y-reentrantn-sdp/2.0), sdp) &               ! Patagonia, east, original
                -sa * spike(x-1.0+dx/2, ll/lonext)*cosbell(y-12.0/latext, 2.5/latext) &                  ! Patagonia, east, extra slope
                -spike(x-dx/2, ll/lonext)*spike(max(0.0, y-reentrants+sdp/2.0), sdp) &                   ! Antarctic Peninsula, west
                -spike(x-1.0+dx/2, ll/lonext)*spike(max(0.0, y-reentrants+sdp/2.0), sdp) &               ! Antarctic Peninsula, east, original
                -sa * spike(x-1.0+dx/2, ll/lonext)*cosbell(y-4.0/latext, 2.5/latext) &                   ! Antarctic Peninsula, east, extra slope
                -spike(y, ll/latext) &                                                                   ! Antarctica
                -spike(y-1.0, ll/latext) &                                                               ! Northern Wall
                - sa *cosbell(x-dx/2-20.0/lonext, 2.5/lonext) * homo(y-8.0/latext, 2.0/latext) &              !Scotia Arc East, center
                - sa * cosbell(x-dx/2-20.0/lonext, 2.5/lonext) * cosbellh(y-10.0/latext-ep, 3.0/latext, 1.) & !Scotia Arc East, north slope
                - sa * cosbell(x-dx/2-20.0/lonext, 2.5/lonext) * cosbellh(y-6.0/latext+ep, 3.0/latext, -1.) & !Scotia Arc East, south slope
                - sa * cosbell(y-12.0/latext, 2.5/latext) * cosbellh(x-dx/2-18.0/lonext, 2.5/lonext, 1.) & !Scotia Arc North, east half (slope side)
                - sa * cosbell(y-12.0/latext, 2.5/latext) * homo(x-dx/2-9.0/lonext, 9.0/lonext) &             !Scotia Arc North, west half
                - sa * cosbell(y-4.0/latext, 2.5/latext) * cosbellh(x-dx/2-18.0/lonext, 2.5/lonext, 1.) &  !Scotia Arc South, east half (slope side)
                - sa * cosbell(y-4.0/latext, 2.5/latext) * homo(x-dx/2-9.0/lonext, 9.0/lonext)                !Scotia Arc South, west half

    ! make sure no deeper than max depth and no shallower than Scotia Arc top IN the ocean interior
    if (D(i,j)<1.0 - sa .and. x>=dx/1.5+ll/2/lonext .and. x<=1.0-ll/2/lonext-dx/1.5 &
        .and. y>=ll/2/latext .and. y<=1.0-ll/2/latext) then
      D(i,j) = 1 - sa
    else if (D(i,j) > 1.0) then
      D(i,j) = 1.0
    endif

      ! make sure the model is not zonally reentrant outside of Drake Passage
      if (((y>=reentrantn+sdp/2 .or. y<=reentrants-sdp/2) .and. x<dx/1.5) &
                .or. ((y>=reentrantn+sdp/2 .or. y<=reentrants-sdp/2) .and. x>1.0-dx/1.5)) then
        D(i,j) = 0.0
      endif

    D(i,j) = D(i,j) * max_depth
  enddo
  enddo


end subroutine shoebox2_initialize_topography

! -----------------------------------------------------------------------------
! define functions used in the above subroutines

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

end module shoebox2_initialization
