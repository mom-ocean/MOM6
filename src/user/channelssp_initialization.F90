module channelssp_initialization

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

public channelssp_initialize_topography
public channelssp_initialize_sponges

! This include declares and sets the variable "version".
#include "version_variable.h"

contains


! -----------------------------------------------------------------------------
!> This subroutine sets up the channelssp test case topography.
!> This sets topography that has a reentrant channel between southern and
!> northern boundaries
!> difference from channls: no slope in the northern boundary, to allow the
!> sponge to work properly

subroutine channelssp_initialize_topography(D, G, param_file, max_depth)
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
  real :: sdp = 4.0, ll=6.0    ! the width of the slope in Drake Passage, and the half width of continental slope

  character(len=40)  :: mod = "channelssp_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  sa = sa_dim / max_depth
  latext = G%len_lat
  lonext = G%len_lon
  reentrants = 6.0/latext          ! non-dimensional southern bound of the reentrant zone
  reentrantn = 10.0/latext         ! non-dimensional northern bound of the reentrant zone
  sdp = sdp/latext
  D = 0.0
  dx = (G%geoLonT(is+1,js)-G%geoLonT(is,js))/lonext     ! longitudinal increment for each grid; non-dimensional

  call MOM_mesg("  channelssp_initialization.F90, channelssp_initialize_topography: setting topography", 5)

  call log_version(param_file, mod, version, "")


  !  Calculate the depth of the bottom.
  do j=js,je                ! meridional grid points
  do i=is,ie                ! zonal grid points
    x=(G%geoLonT(i,j)-G%west_lon) / lonext      ! non-dimensional longitude
    y=(G%geoLatT(i,j)-G%south_lat) / latext     ! non-dimensional latitude

    D(i,j) = 1.0 - spike(x-dx/2, ll/lonext)*spike(min(0.0, y-reentrantn-sdp/2.0), sdp) &                 ! Patagonia, west
                -spike(x-1.0+dx/2, ll/lonext)*spike(min(0.0, y-reentrantn-sdp/2.0), sdp) &               ! Patagonia, east
                -spike(x-dx/2, ll/lonext)*spike(max(0.0, y-reentrants+sdp/2.0), sdp) &                   ! Antarctic Peninsula, west
                -spike(x-1.0+dx/2, ll/lonext)*spike(max(0.0, y-reentrants+sdp/2.0), sdp) &               ! Antarctic Peninsula, east
                -spike(y, ll/latext) &                                                                   ! Antarctica
                - sa * cosbell(x-20.0/lonext, 2.5/lonext) * homo(y-8.0/latext, 2.0/latext) &              !Scotia Arc East, center
                - sa * cosbell(x-20.0/lonext, 2.5/lonext) * cosbellh(y-10.0/latext-ep, 3.0/latext, 1.) & !Scotia Arc East, north slope
                - sa * cosbell(x-20.0/lonext, 2.5/lonext) * cosbellh(y-6.0/latext+ep, 3.0/latext, -1.) & !Scotia Arc East, south slope
                - sa * cosbell(y-12.0/latext, 2.5/latext) * cosbellh(x-18.0/lonext-ep, 2.5/lonext, 1.) &    !Scotia Arc North, east half (slope side)
                - sa * cosbell(y-12.0/latext, 2.5/latext) * homo(x-9.0/lonext, 9.0/lonext) &             !Scotia Arc North, west half
                - sa * cosbell(y-4.0/latext, 2.5/latext) * cosbellh(x-18.0/lonext-ep, 2.5/lonext, 1.) &     !Scotia Arc South, east half (slope side)
                - sa * cosbell(y-4.0/latext, 2.5/latext) * homo(x-9.0/lonext, 9.0/lonext)                !Scotia Arc South, west half

      ! make sure no deeper than max depth and no shallower than Scotia Arc top IN the ocean interior
      if (D(i,j)<1.0 - sa .and. x>=dx+ll/2/lonext .and. x<=1.0-ll/2/lonext-dx .and. y>=ll/2/latext .and. y<=1.0-ll/2/latext) then
        D(i,j) = 1 - sa
      else if (D(i,j) > 1.0) then
        D(i,j) = 1.0
      endif

    D(i,j) = D(i,j) * max_depth
  enddo
  enddo


end subroutine channelssp_initialize_topography


! -----------------------------------------------------------------------------
!> set up the sponge used for channelssp topography
!> difference from the sponge designed in channel_initialization.F90: 
!> no sponge to restore interface heights where there is continental slope
!> specific for channelssp: the northermost several degrees of west/east boundaries

subroutine channelssp_initialize_sponges(G, GV, use_temperature, tv, param_file, CSp, h)
  type(ocean_grid_type), intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV  !< The ocean's vertical grid structure
                                             ! so as to be used as input of
                                             ! set_coord_from_file
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
  real, dimension(SZK_(G)) :: rho           ! layer densities
  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1) ! A temporary array for eta, m.
                                         ! eta only varies in z
  real :: Idamp(SZI_(G),SZJ_(G))         ! The inverse damping rate, in s-1.
                                         ! Idamp only varies in y; /= 0 only
                                         ! within the sponge layer
  real :: eta0(SZK_(G)+1)                ! target interface heights for density class in the sponge layer
  real :: damp_rate, damp, spongelen = 2.0, etab, c = 7.7, rhot, rhob, min_depth, intrho, nlat, ll, dx
  ! spongelen: thickness of sponge layer in dimensional degree
  ! etab: bottom ocean height; rhot & rhob:  biggest/smallest rho in the sponge
  ! layer; 
  ! c: constant in fitting eta = eta(rho); intrho is interface rho; 
  ! ll is the double width of continental slope, dx is non-dimensional
  ! longitudinal grid increment
  integer :: i, j, k, is, ie, js, je, nz
  logical, save :: first_call = .true.
  character(len=40)  :: mdl = "channelssp_initialize_sponges" ! This subroutine's name

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  eta(:,:,:) = 0.0 ; Idamp(:,:) = 0.0; eta0(:) = 0.0; etab = - G%max_depth;
  ll = 3.0/G%len_lon            ! non-dimensional longitudinal range where continental slope is located  
  dx = (G%geoLonT(is+1,js)-G%geoLonT(is,js))/G%len_lon

  rho = GV%Rlay         ! loaded target potential density for each layer
  rhot = 1025.212       ! target surface rho in the sponge
  rhob = 1027.412       ! target bottom interface rho in the sponge

  if (first_call) call log_version(param_file, mdl, version)
  first_call = .false.


  call get_param(param_file, mdl, "SPONGE_RATE", damp_rate, &
                 "The rate at which the zonal-mean sponges damp.", units="s-1", &
                 default = 1.0/(10.0*86400.0))
  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)


  nlat = G%south_lat + G%len_lat        ! should be -30.0 degree

  ! compute target interface heights, stored in eta0
  do k = 4, nz-2                        ! top 3 interface heights = 0.0, so start from k = 4
    intrho = (rho(k-1)+rho(k))/2.0      ! interface rho
    eta0(k) = etab*(1-1/c*log(1-(exp(c)-1)*(intrho-rhob)/(rhob-rhot)))  ! target height at each interface
  enddo
  do k = nz-1, nz+1
    eta0(k) = etab                      ! interface heights reaches ocean bottom
  enddo

  ! initialize the damping rate so it is 0 outside of the sponge layer &
  ! and increases linearly with latitude within the sponge layer
  do i = is, ie; do j = js, je
    if ((G%geoLatT(i,j) >= nlat-spongelen .and. G%geoLonT(i,j)>=ll+dx/2) .or. &
        (G%geoLatT(i,j) >= nlat-spongelen .and. G%geoLonT(i,j)<=1.0-ll-dx/2)) then ! the sponge has meridional AND zonal range
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


end subroutine channelssp_initialize_sponges



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

    cosbellh  = cos(PI/2.0*MIN(xx/L, 1.0))
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

end module channelssp_initialization
