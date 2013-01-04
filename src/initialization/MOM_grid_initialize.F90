module MOM_grid_initialize
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

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, November 1998 - June 2002                      *
!*                                                                     *
!*    This program contains 2 externally callable subroutines.         *
!*  set_grid_metrics calculates the various metric terms that are used *
!*  by MOM.  This routine is intended to be modified by the user to    *
!*  enable the use of any general orthogonal grid.  initialize_masks   *
!*  initializes the land masks; it is in this file because it a key    *
!*  part of the physical grid description.                             *
!*                                                                     *
!*    This subroutine is also used by MOM-related preprocessing and    *
!*  postprocessing codes.                                              *
!*                                                                     *
!*    The metric terms have the form Dzp, IDzp, or DXDYp, where z can  *
!*  be X or Y, and p can be q, u, v, or h.  z describes the direction  *
!*  of the metric, while p describes the location.  IDzp is the        *
!*  inverse of Dzp, while DXDYp is the product of DXp and DYp except   *
!*  that areaT is calculated analytically from the latitudes and       *
!*  longitudes of the surrounding q points.                            *
!*                                                                     *
!*    On a sphere, a variety of grids can be implemented by defining   *
!*  analytic expressions for dx_di, dy_dj (where x and y are latitude  *
!*  and longitude, and i and j are grid indices) and the expressions   *
!*  for the integrals of their inverses in the four subroutines        *
!*  dy_dj, Int_dj_dy, dx_di, and Int_di_dx.                            *
!*                                                                     *
!*    initialize_masks sets up land masks based on the depth field.    *
!*  The one argument is the minimum ocean depth.  Depths that are      *
!*  less than this are interpreted as land points.                     *
!*                                                                     *
!*    Macros written all in capital letters are from MOM_memory.h.     *
!*                                                                     *
!*     A small fragment of the C-grid is shown below:                  *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, DXq, IDXq, DYq, IDYq, etc.            *
!*    j+1  > o > o >   At ^:  v, DXv, IDXv, DYv, IDYv, etc.            *
!*    j    x ^ x ^ x   At >:  u, DXu, IDXu, DYu, IDYu, etc.            *
!*    j    > o > o >   At o:  h, DXh, IDXh, DYh, IDYh, areaT, etc.     *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_diag_mediator, only : set_axes_info
use MOM_domains, only : pass_var, pass_vector, pe_here, root_PE, broadcast
use MOM_checksums, only : hchksum, qchksum, uchksum, vchksum
use MOM_domains, only : AGRID, BGRID_NE, CGRID_NE, To_All, Scalar_Pair
use MOM_domains, only : To_North, To_South, To_East, To_West
use MOM_domains, only : MOM_define_domain, MOM_define_IO_domain
use MOM_domains, only : MOM_domain_type
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, NOTE, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : read_data, slasher, file_exists
use MOM_io, only : CORNER, NORTH_FACE, EAST_FACE
use mpp_domains_mod, only : mpp_get_compute_domain,mpp_get_compute_domains
use mpp_domains_mod, only : mpp_get_data_domain
use mpp_domains_mod, only : domain1D, mpp_get_domain_components

implicit none ; private

#include <MOM_memory.h>

public set_grid_metrics, initialize_masks

type, public :: GPS ; private
  real :: len_lon
  real :: len_lat
  real :: west_lon
  real :: south_lat
  real :: Rad_Earth
  real :: Lat_enhance_factor
  real :: Lat_eq_enhance
  logical :: isotropic
  logical :: equator_reference
  integer :: niglobal, njglobal         ! Duplicates of niglobal and njglobal from MOM_dom
end type GPS

real, parameter :: Epsln = 1.0e-10  !   A distance used to replace negative
                                    ! distances in the metric arrays.

contains

function ds_di(x, y, GP)
  real, intent(in) :: x, y
  type(GPS), intent(in) :: GP
  real :: ds_di
! This function returns the grid spacing in the logical x direction.
! Arguments: x - The latitude in question.
!  (in)      y - The longitude in question.
  ds_di = GP%Rad_Earth * cos(y) * dx_di(x,GP)
! In general, this might be...
! ds_di = GP%Rad_Earth * sqrt( cos(y)*cos(y) * dx_di(x,y,GP)*dx_di(x,y,GP) + &
!                           dy_di(x,y,GP)*dy_di(x,y,GP))
end function ds_di

function ds_dj(x, y, GP)
  real, intent(in) :: x, y
  type(GPS), intent(in) :: GP
  real :: ds_dj
! This function returns the grid spacing in the logical y direction.
! Arguments: x - The latitude in question.
!  (in)      y - The longitude in question.
  ds_dj = GP%Rad_Earth * dy_dj(y,GP)
! In general, this might be...
! ds_dj = GP%Rad_Earth * sqrt( cos(y)*cos(y) * dx_dj(x,y,GP)*dx_dj(x,y,GP) + &
!                           dy_dj(x,y,GP)*dy_dj(x,y,GP))
end function ds_dj


function  dL(x1, x2, y1, y2)
  real, intent(in) :: x1, x2, y1, y2
  real :: dL
!  This subroutine calculates the contribution from the line integral
! along one of the four sides of a cell face to the area of a cell,
! assuming that the sides follow a linear path in latitude and long-
! itude (i.e., on a Mercator grid).
! Argumnts: x1 - Segment starting longitude.
!  (in)     x2 - Segment ending longitude.
!  (in)     y1 - Segment ending latitude.
!  (in)     y2 - Segment ending latitude.
  real :: r, dy

  dy = y2 - y1

  if (ABS(dy) > 2.5e-8) then
    r = ((1.0 - cos(dy))*cos(y1) + sin(dy)*sin(y1)) / dy
  else
    r = (0.5*dy*cos(y1) + sin(y1))
  endif
  dL = r * (x2 - x1)

end function  dL

function find_root( fn, dy_df, GP, fnval, y1, ymin, ymax, ittmax)
  real :: find_root
  real, external :: fn, dy_df
  type(GPS), intent(in) :: GP
  real, intent(in) :: fnval, y1, ymin, ymax
  integer, intent(out) :: ittmax
  real :: y
! This subroutine finds and returns the value of y at which the
! monotonic function fn takes the value fnval, also returning
! in ittmax the number of iterations of Newton's method that were
! used to polish the root.
  real :: ybot, ytop, fnbot, fntop
  integer :: itt
  character(len =256) :: warnmesg

  real :: dy_dfn, dy, fny

! For Fortran we need to copy the input y1 to y.
! Otherwise the value of y gets changed globally.
! i.e. y_h = find_root(fn,dy_df,fnval, y_q, ...) modifies y_q in
! the code above

!  Bracket the root.
  y = y1 ; ybot = y1
  fnbot = fn(ybot,GP) - fnval
  itt = 0
  do while (fnbot > 0.0)
    if ((ybot - 2.0*dy_df(ybot,GP)) < (0.5*(ybot+ymin))) then
      ybot = ybot - 2.0*dy_df(ybot,GP)
    else
      ybot = 0.5*(ybot+ymin) ; itt = itt + 1
    endif
    fnbot = fn(ybot,GP) - fnval

    if ((itt > 50) .and. (fnbot > 0.0)) then
      write(warnmesg, '("PE ",I2," unable to find bottom bound for grid function. &
        &x = ",ES10.4,", xmax = ",ES10.4,", fn = ",ES10.4,", dfn_dx = ",ES10.4,&
        &", seeking fn = ",ES10.4," - fn = ",ES10.4,".")') &
          pe_here(),ybot,ymin,fn(ybot,GP),dy_df(ybot,GP),fnval, fnbot

      call MOM_error(FATAL,warnmesg)
    endif
  enddo

  if ((y + 2.0*dy_df(y,GP)) < (0.5*(y+ymax))) then
    ytop = y + 2.0*dy_df(y,GP)
  else
    ytop = 0.5*(y+ymax)
  endif
  fntop = fn(ytop,GP) - fnval ; itt = 0
  do while (fntop < 0.0)
    if ((ytop + 2.0*dy_df(ytop,GP)) < (0.5*(ytop+ymax))) then
      ytop = ytop + 2.0*dy_df(ytop,GP)
    else
      ytop = 0.5*(ytop+ymax) ; itt = itt + 1
    endif
    fntop = fn(ytop,GP) - fnval

    if ((itt > 50) .and. (fntop < 0.0)) then
      write(warnmesg, '("PE ",I2," unable to find top bound for grid function. &
        &x = ",ES10.4,", xmax = ",ES10.4,", fn = ",ES10.4,", dfn_dx = ",ES10.4, &
        &", seeking fn = ",ES10.4," - fn = ",ES10.4,".")') &
          pe_here(),ytop,ymax,fn(ytop,GP),dy_df(ytop,GP),fnval,fntop

      call MOM_error(FATAL,warnmesg)
    endif
  enddo
!  Bisect several times to insure that the root is within the radius
!  of convergence in the Newton's method polisher.
  do itt=1,10
    y = 0.5*(ybot + ytop)
    fny = fn(y,GP) - fnval
    if (fny < 0.0) then
      fnbot = fny ; ybot = y
    else
      fntop = fny ; ytop = y
    endif
  enddo

!    Polish the root using Newton's method.
  do itt=1,10
    dy_dfn = dy_df(y,GP)
    fny = fn(y,GP) - fnval

    dy = -1.0* fny * dy_dfn
    y = y + dy
    if (y > ytop) y = ytop
    if (y < ybot) y = ybot
    if (ABS(dy) < (8.0e-15*ABS(y)+1.e-20)) exit
  enddo
  if (ABS(y) < 1e-12) y = 0.0

  ittmax = itt
  find_root = y
end function find_root

function dx_di(x, GP)
  real, intent(in) :: x
  type(GPS), intent(in) :: GP
  real :: dx_di
! This subroutine calculates and returns the value of dx/di, where
! x is the longitude in Radians, and i is the integral north-south
! grid index.

  dx_di = (GP%len_lon * 4.0*atan(1.0)) / (180.0 * GP%niglobal)

end function dx_di

function Int_di_dx(x, GP)
  real, intent(in) :: x
  type(GPS), intent(in) :: GP
  real :: Int_di_dx
! This subroutine calculates and returns the integral of the inverse
! of dx/di to the point x, in radians.

  Int_di_dx = x * ((180.0 * GP%niglobal) / (GP%len_lon * 4.0*atan(1.0)))

end function Int_di_dx

function dy_dj(y, GP)
  real, intent(in) :: y
  type(GPS), intent(in) :: GP
  real :: dy_dj
! This subroutine calculates and returns the value of dy/dj, where
! y is the latitude in Radians, and j is the integral north-south
! grid index.
  real :: PI            ! 3.1415926... calculated as 4*atan(1)
  real :: C0            ! The constant that converts the nominal y-spacing in
                        ! gridpoints to the nominal spacing in Radians.
  real :: y_eq_enhance  ! The latitude in radians within which the resolution
                        ! is enhanced.
  PI = 4.0*atan(1.0)
  if (GP%isotropic) then
    C0 = (GP%len_lon * PI) / (180.0 * GP%niglobal)
    y_eq_enhance = PI*abs(GP%lat_eq_enhance)/180.0
    if (ABS(y) < y_eq_enhance) then
      dy_dj = C0 * (cos(y) / (1.0 + 0.5*cos(y) * (GP%lat_enhance_factor - 1.0) * &
                         (1.0+cos(PI*y/y_eq_enhance)) ))
    else
      dy_dj = C0 * cos(y)
    endif
  else
    C0 = (GP%len_lat * PI) / (180.0 * GP%njglobal)
    dy_dj = C0
  endif

end function dy_dj

function Int_dj_dy(y, GP)
  real, intent(in) :: y
  type(GPS), intent(in) :: GP
  real :: Int_dj_dy
! This subroutine calculates and returns the integral of the inverse
! of dy/dj to the point y, in radians.
  real :: I_C0 = 0.0       !   The inverse of the constant that converts the
                           ! nominal spacing in gridpoints to the nominal
                           ! spacing in Radians.
  real :: PI               ! 3.1415926... calculated as 4*atan(1)
  real :: y_eq_enhance     ! The latitude in radians from
                           ! from the equator within which the
                           ! meridional grid spacing is enhanced by
                           ! a factor of GP%lat_enhance_factor.
  real :: r

  PI = 4.0*atan(1.0)
  if (GP%isotropic) then
    I_C0 = (180.0 * GP%niglobal) / (GP%len_lon * PI)
    y_eq_enhance = PI*ABS(GP%lat_eq_enhance)/180.0

    if (y >= 0.0) then
      r = I_C0 * log((1.0 + sin(y))/cos(y))
    else
      r = -1.0 * I_C0 * log((1.0 - sin(y))/cos(y))
    endif

    if (y >= y_eq_enhance) then
      r = r + I_C0*0.5*(GP%lat_enhance_factor - 1.0)*y_eq_enhance
    else if (y <= -y_eq_enhance) then
      r = r - I_C0*0.5*(GP%lat_enhance_factor - 1.0)*y_eq_enhance
    else
      r = r + I_C0*0.5*(GP%lat_enhance_factor - 1.0) * &
              (y + (y_eq_enhance/PI)*sin(PI*y/y_eq_enhance))
    endif
  else
    I_C0 = (180.0 * GP%njglobal) / (GP%len_lat * PI)
    r = I_C0 * y
  endif

  Int_dj_dy = r
end function Int_dj_dy

! ------------------------------------------------------------------------------

subroutine set_grid_metrics(G, param_file, set_vertical)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
  logical, optional,     intent(in)    :: set_vertical
! Arguments:
!  (inout)   G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in,opt)  set_vertical - If true (or missing), set up the vertical axes.

!    Calculate the values of the metric terms that might be used
!  and save them in arrays.
!    Within this subroutine, the x- and y- grid spacings and their
!  inverses and the cell areas centered on h, q, u, and v points are
!  calculated, as are the geographic locations of each of these 4
!  sets of points.
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'
  logical :: debug
  character(len=256) :: config

  call MOM_mesg("  MOM_grid_init.F90, set_grid_metrics: allocating metrics", 5)
 
  call allocate_metrics(G)

  call log_version(param_file, "MOM_grid_init", version, tagname, "")
  call get_param(param_file, "MOM_grid_init", "GRID_CONFIG", config, &
                 "A character string that determines the method for \n"//&
                 "defining the horizontal grid.  Current options are: \n"//&
                 " \t mosaic - read the grid from a mosaic (supergrid) \n"//&
                 " \t          file set by GRID_FILE.\n"//&
                 " \t cartesian - use a (flat) Cartesian grid.\n"//&
                 " \t spherical - use a simple spherical grid.\n"//&
                 " \t mercator - use a Mercator spherical grid.", &
                 fail_if_missing=.true.)
  call get_param(param_file, "MOM_grid_init", "DEBUG", debug, &
                 "If true, write out verbose debugging data.", default=.false.)

  select case (trim(config))
    case ("mosaic");    call set_grid_metrics_from_mosaic(G, param_file)
    case ("cartesian"); call set_grid_metrics_cartesian(G, param_file)
    case ("spherical"); call set_grid_metrics_spherical(G, param_file)
    case ("mercator");  call set_grid_metrics_mercator(G, param_file)
    case ("file"); call MOM_error(FATAL, "MOM_grid_init: set_grid_metrics "//&
           'GRID_CONFIG "file" is no longer a supported option.  Use a '//&
           'mosaic file ("mosaic") or one of the analytic forms instead.')
    case default ; call MOM_error(FATAL, "MOM_grid_init: set_grid_metrics "//&
           "Unrecognized grid configuration "//trim(config))
  end select

! Calculate derived metrics (i.e. reciprocals and products)
  call set_grid_derived_metrics(G, param_file)

! This call sets up the diagnostic axes.
  call set_axes_info(G%gridlatq, G%gridlath, G%gridlonq, G%gridlonh, G, &
                     param_file, set_vertical=set_vertical)

  if (debug) call grid_metrics_chksum('MOM_grid_init/set_grid_metrics',G)

end subroutine set_grid_metrics

! ------------------------------------------------------------------------------

subroutine set_grid_derived_metrics(G, param_file)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
! Arguments:
!  (inout)   G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

!    Calculate the values of the metric terms that might be used
!  and save them in arrays.
!    Within this subroutine, the x- and y- grid spacings and their
!  inverses and the cell areas centered on h, q, u, and v points are
!  calculated, as are the geographic locations of each of these 4
!  sets of points.
  character( len = 128) :: warnmesg
  integer :: i,j, isd, ied, jsd, jed
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, Isdq, Iedq, Jsdq, Jedq

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq
  Isdq = G%Isdq ; Iedq = G%Iedq ; Jsdq = G%Jsdq ; Jedq = G%Jedq

  call MOM_mesg("  MOM_grid_init.F90, set_grid_derived_metrics: deriving metrics", 5)
 
  do j=jsd,jed ; do i=isd,ied
    if (G%DXh(i,j) <= 0.0) then
      write(warnmesg,68)  pe_here(),"DXh",i,j,G%DXh(i,j),Epsln
      call MOM_error(NOTE, warnmesg, all_print=.true.)
      G%DXh(i,j) = Epsln
    endif
    if (G%DYh(i,j) <= 0.0) then
      write(warnmesg,68)  pe_here(),"DYh",i,j,G%DYh(i,j),Epsln
      call MOM_error(NOTE, warnmesg, all_print=.true.)
      G%DYh(i,j) = Epsln
    endif
    G%IDXh(i,j) = 1.0 / G%DXh(i,j)
    G%IDYh(i,j) = 1.0 / G%DYh(i,j)
    G%IareaT(i,j) = 1.0 / G%areaT(i,j)
  enddo ; enddo

  do j=jsd,jed ; do I=Isdq,Iedq
    if (G%DXu(I,j) <= 0.0) then
      write(warnmesg,68)  pe_here(),"DXu",I,j,G%DXu(I,j),Epsln
      call MOM_error(NOTE, warnmesg, all_print=.true.)
      G%DXu(I,j) = Epsln
    endif
    if (G%DYu(I,j) <= 0.0) then
      write(warnmesg,68)  pe_here(),"DYu",I,j,G%DYu(I,j),Epsln
      call MOM_error(NOTE, warnmesg, all_print=.true.)
      G%DYu(I,j) = Epsln
    endif
    G%IDXu(i,j) = 1.0 / G%DXu(i,j)
    G%IDYu(i,j) = 1.0 / G%DYu(i,j)
  enddo ; enddo

  do J=Jsdq,Jedq ; do i=isd,ied
    if (G%DXv(i,j) <= 0.0) then
      write(warnmesg,68)  pe_here(),"DXv",i,j,G%DXv(i,j),Epsln
      call MOM_error(NOTE, warnmesg, all_print=.true.)
      G%DXv(i,j) = Epsln
    endif
    if (G%DYv(i,j) <= 0.0) then
      write(warnmesg,68)  pe_here(),"DYv",i,j,G%DYv(i,j),Epsln
      call MOM_error(NOTE, warnmesg, all_print=.true.)
      G%DYv(i,j) = Epsln
    endif
    G%IDXv(i,j) = 1.0 / G%DXv(i,j)
    G%IDYv(i,j) = 1.0 / G%DYv(i,j)
  enddo ; enddo

  do J=Jsdq,Jedq ; do I=Isdq,Iedq
    if (G%DXq(I,J) <= 0.0) then
      write(warnmesg,68)  pe_here(),"DXq",I,J,G%DXq(I,J),Epsln
      call MOM_error(NOTE, warnmesg, all_print=.true.)
      G%DXq(I,J) = Epsln
    endif
    if (G%DYq(I,J) <= 0.0) then
      write(warnmesg,68)  pe_here(),"DYq",I,J,G%DYq(I,J),Epsln
      call MOM_error(NOTE, warnmesg, all_print=.true.)
      G%DYq(I,J) = Epsln
    endif

    G%IDXq(I,J) = 1.0 / G%DXq(I,J)
    G%IDYq(I,J) = 1.0 / G%DYq(I,J)
    G%areaBu(I,J) = G%DXq(I,J) * G%DYq(I,J)
    G%IareaBu(I,J) = G%IDXq(I,J) * G%IDYq(I,J)
  enddo ; enddo

68 FORMAT ("WARNING: PE ",I4," ",a3,"(",I4,",",I4,") = ",ES10.4, &
           " is being changed to ",ES10.4,".")

end subroutine set_grid_derived_metrics

! ------------------------------------------------------------------------------

subroutine grid_metrics_chksum(parent, G)
  character(len=*),      intent(in) :: parent
  type(ocean_grid_type), intent(in) :: G
! Arguments:
!  (in)          parent - String indentifying caller
!  (in)               G - The ocean's grid structure.
  real, dimension(G%isd :G%ied ,G%jsd :G%jed ) :: tempH
  real, dimension(G%Isdq:G%Iedq,G%Jsdq:G%Jedq) :: tempQ
  real, dimension(G%Isdq:G%Iedq,G%jsd :G%jed ) :: tempE
  real, dimension(G%isd :G%ied ,G%Jsdq:G%Jedq) :: tempN
  integer :: i, j, isd, ied, jsd, jed, halo
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, Isdq, Iedq, Jsdq, Jedq

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq
  Isdq = G%Isdq ; Iedq = G%Iedq ; Jsdq = G%Jsdq ; Jedq = G%Jedq
  halo = min(ied-ie, jed-je)
halo=1 ! AJA

  do i=isd,ied ; do j=jsd,jed ; tempH(i,j) = G%DXh(i,j) ; enddo ; enddo
  call hchksum(tempH,trim(parent)//': DXh',G,haloshift=halo)

  do I=Isdq,Iedq ; do j=jsd,jed ; tempE(I,j) = G%DXu(I,j) ; enddo ; enddo
  call uchksum(tempE,trim(parent)//': DXu',G,haloshift=halo)

  do i=isd,ied ; do J=Jsdq,Jedq ; tempN(i,J) = G%DXv(i,J) ; enddo ; enddo
  call vchksum(tempN,trim(parent)//': DXv',G,haloshift=halo)
 
  do I=Isdq,Iedq ; do J=Jsdq,Jedq ; tempQ(I,J) = G%DXq(I,J) ; enddo ; enddo
  call qchksum(tempQ,trim(parent)//': DXq',G,haloshift=halo)
 
  do i=isd,ied ; do j=jsd,jed ; tempH(i,j) = G%DYh(i,j) ; enddo ; enddo
  call hchksum(tempH,trim(parent)//': DYh',G,haloshift=halo)
 
  do I=Isdq,Iedq ; do j=jsd,jed ; tempE(I,j) = G%DYu(I,j) ; enddo ; enddo
  call uchksum(tempE,trim(parent)//': DYu',G,haloshift=halo)
 
  do i=isd,ied ; do J=Jsdq,Jedq ; tempN(i,J) = G%DYv(i,J) ; enddo ; enddo
  call vchksum(tempN,trim(parent)//': DYv',G,haloshift=halo)
 
  do I=Isdq,Iedq ; do J=Jsdq,Jedq ; tempQ(I,J) = G%DYq(I,J) ; enddo ; enddo
  call qchksum(tempQ,trim(parent)//': DYq',G,haloshift=halo)

  do i=isd,ied ; do j=jsd,jed ; tempH(i,j) = G%IDXh(i,j) ; enddo ; enddo
  call hchksum(tempH,trim(parent)//': IDXh',G,haloshift=halo)

  do I=Isdq,Iedq ; do j=jsd,jed ; tempE(I,j) = G%IDXu(I,j) ; enddo ; enddo
  call uchksum(tempE,trim(parent)//': IDXu',G,haloshift=halo)
 
  do i=isd,ied ; do J=Jsdq,Jedq ; tempN(i,J) = G%IDXv(i,J) ; enddo ; enddo
  call vchksum(tempN,trim(parent)//': IDXv',G,haloshift=halo)
 
  do I=Isdq,Iedq ; do J=Jsdq,Jedq ; tempQ(I,J) = G%IDXq(I,J) ; enddo ; enddo
  call qchksum(tempQ,trim(parent)//': IDXq',G,haloshift=halo)

  do i=isd,ied ; do j=jsd,jed ; tempH(i,j) = G%IDYh(i,j) ; enddo ; enddo
  call hchksum(tempH,trim(parent)//': IDYh',G,haloshift=halo)
 
  do I=Isdq,Iedq ; do j=jsd,jed ; tempE(I,j) = G%IDYu(I,j) ; enddo ; enddo
  call uchksum(tempE,trim(parent)//': IDYu',G,haloshift=halo)
 
  do i=isd,ied ; do J=Jsdq,Jedq ; tempN(i,J) = G%IDYv(i,J) ; enddo ; enddo
  call vchksum(tempN,trim(parent)//': IDYv',G,haloshift=halo)
 
  do I=Isdq,Iedq ; do J=Jsdq,Jedq ; tempQ(I,J) = G%IDYq(I,J) ; enddo ; enddo
  call qchksum(tempQ,trim(parent)//': IDYq',G,haloshift=halo)

  do i=isd,ied ; do j=jsd,jed ; tempH(i,j) = G%areaT(i,j) ; enddo ; enddo
  call hchksum(tempH,trim(parent)//': areaT',G,haloshift=halo)
 
  do I=Isdq,Iedq ; do J=Jsdq,Jedq ; tempQ(I,J) = G%areaBu(I,J) ; enddo ; enddo
  call qchksum(tempQ,trim(parent)//': areaBu',G,haloshift=halo)
 
  do i=isd,ied ; do j=jsd,jed ; tempH(i,j) = G%IareaT(i,j) ; enddo ; enddo
  call hchksum(tempH,trim(parent)//': IareaT',G,haloshift=halo)
 
  do I=Isdq,Iedq ; do J=Jsdq,Jedq ; tempQ(I,J) = G%IareaBu(I,J) ; enddo ; enddo
  call qchksum(tempQ,trim(parent)//': IareaBu',G,haloshift=halo)

  call hchksum(G%geolonh,trim(parent)//': geolonh',G,haloshift=halo)

  call hchksum(G%geolath,trim(parent)//': geolath',G,haloshift=halo)

  do I=Isdq,Iedq ; do J=Jsdq,Jedq ; tempQ(I,J) = G%geolonq(I,J) ; enddo ; enddo
  call qchksum(tempQ,trim(parent)//': geolonq',G,haloshift=halo)

  do I=Isdq,Iedq ; do J=Jsdq,Jedq ; tempQ(I,J) = G%geolatq(I,J) ; enddo ; enddo
  call qchksum(tempQ,trim(parent)//': geolatq',G,haloshift=halo)

  do I=Isdq,Iedq ; do j=jsd,jed ; tempE(I,J) = G%geolonu(I,J) ; enddo ; enddo
  call uchksum(tempE,trim(parent)//': geolonu',G,haloshift=halo)

  do I=Isdq,Iedq ; do j=jsd,jed ; tempE(I,J) = G%geolatu(I,J) ; enddo ; enddo
  call uchksum(tempE,trim(parent)//': geolatu',G,haloshift=halo)

  do i=isd,ied ; do J=Jsdq,Jedq ; tempN(I,J) = G%geolonv(I,J) ; enddo ; enddo
  call vchksum(tempN,trim(parent)//': geolonv',G,haloshift=halo)

  do i=isd,ied ; do J=Jsdq,Jedq ; tempN(I,J) = G%geolatv(I,J) ; enddo ; enddo
  call vchksum(tempN,trim(parent)//': geolatv',G,haloshift=halo)

end subroutine grid_metrics_chksum

! ------------------------------------------------------------------------------

subroutine set_grid_metrics_from_mosaic(G,param_file)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
! Arguments:
!  (inout)   G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

  real, dimension(G%isd :G%ied ,G%jsd :G%jed ) :: tempH1, tempH2, tempH3, tempH4
  real, dimension(G%Isdq:G%Iedq,G%Jsdq:G%Jedq) :: tempQ1, tempQ2, tempQ3, tempQ4
  real, dimension(G%Isdq:G%Iedq,G%jsd :G%jed ) :: tempE1, tempE2
  real, dimension(G%isd :G%ied ,G%Jsdq:G%Jedq) :: tempN1, tempN2
  ! These arrays are a holdover from earlier code in which the arrays in G were
  ! macros and may have had reduced dimensions.
  real, dimension(G%isd :G%ied ,G%jsd :G%jed ) :: dxh, dyh, areaT
  real, dimension(G%Isdq:G%Iedq,G%jsd :G%jed ) :: dxu, dyu
  real, dimension(G%isd :G%ied ,G%Jsdq:G%Jedq) :: dxv, dyv
  real, dimension(G%Isdq:G%Iedq,G%Jsdq:G%Jedq) :: dxq, dyq, areaBu
  ! This are symmetric arrays, corresponding to the data in the mosaic file
  real, dimension(2*G%isd-1:2*G%ied,2*G%jsd-1:2*G%jed) :: tmpT
  real, dimension(2*G%isd-2:2*G%ied,2*G%jsd-1:2*G%jed) :: tmpU
  real, dimension(2*G%isd-1:2*G%ied,2*G%jsd-2:2*G%jed) :: tmpV
  real, dimension(2*G%isd-2:2*G%ied,2*G%jsd-2:2*G%jed) :: tmpZ
  real, dimension(:,:), allocatable :: tmpGlbl
  character(len=200) :: filename, grid_file, inputdir
  character(len=64)  :: mod="MOM_grid_init set_grid_metrics_from_mosaic"
  integer :: err=0, dv(2,5), ni, nj, global_indices(4)
  type(MOM_domain_type) :: SGdom ! Supergrid domain
  integer :: i, j
  integer :: npei,npej
  integer :: isc,jsc,iec,jec ! Computational domain extents on the supergrid.
  integer :: isd,jsd,ied,jed ! Memory extents on the supergrid.
  integer, dimension(:), allocatable :: exni,exnj
  type(domain1D) :: domx, domy
  integer        :: start(4), nread(4)
 
  call MOM_mesg("   MOM_grid_init.F90, set_grid_metrics_from_mosaic: reading grid", 5)

  call get_param(param_file, mod, "GRID_FILE", grid_file, &
                 "Name of the file from which to read horizontal grid data.", &
                 fail_if_missing=.true.)
  call get_param(param_file,  mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  filename = trim(adjustl(inputdir)) // trim(adjustl(grid_file))
  call log_param(param_file, mod, "INPUTDIR/GRID_FILE", filename)
  if (.not.file_exists(filename)) &
    call MOM_error(FATAL," set_grid_metrics_from_mosaic: Unable to open "//&
                           trim(filename))

! Initialize everything to a small number
  dxu(:,:)=Epsln; dyu(:,:)=Epsln
  dxv(:,:)=Epsln; dyv(:,:)=Epsln
  dxq(:,:)=Epsln; dyq(:,:)=Epsln; areaBu(:,:)=Epsln

!<MISSING CODE TO READ REFINEMENT LEVEL>
  ni=2*(G%iec-G%isc+1) ! i size of supergrid
  nj=2*(G%jec-G%jsc+1) ! j size of supergrid
  dv(1,1)=2*G%isd-1
  dv(1,2)=2*G%ied
  dv(1,3)=2*G%isc-1
  dv(1,4)=2*G%iec
  dv(1,5)=2*G%isd_global-1 ! location of tmpT(1,1) in data file
  dv(2,1)=2*G%jsd-1
  dv(2,2)=2*G%jed
  dv(2,3)=2*G%jsc-1
  dv(2,4)=2*G%jec
  dv(2,5)=2*G%jsd_global-1 ! location of tmpT(1,1) in data file

! Define a domain for the supergrid (SGdom)
  npei=G%domain%layout(1)
  npej=G%domain%layout(2)
  allocate(exni(npei))
  allocate(exnj(npej))
  call mpp_get_domain_components(G%domain%mpp_domain, domx, domy)
  call mpp_get_compute_domains(domx,size=exni)
  call mpp_get_compute_domains(domy,size=exnj)
  allocate(SGdom%mpp_domain)
  SGdom%nihalo = 2*G%domain%nihalo
  SGdom%njhalo = 2*G%domain%njhalo
  SGdom%niglobal = 2*G%domain%niglobal
  SGdom%njglobal = 2*G%domain%njglobal
  SGdom%layout(:) = G%domain%layout(:)
  SGdom%use_io_layout = G%domain%use_io_layout
  SGdom%io_layout(:) = G%domain%io_layout(:)
  global_indices(1) = 1+SGdom%nihalo
  global_indices(2) = SGdom%niglobal+SGdom%nihalo
  global_indices(3) = 1+SGdom%njhalo
  global_indices(4) = SGdom%njglobal+SGdom%njhalo
  exni(:)=2*exni(:); exnj(:)=2*exnj(:)
  call MOM_define_domain(global_indices, SGdom%layout, SGdom%mpp_domain, &
         xflags=G%domain%X_FLAGS, yflags=G%domain%Y_FLAGS, &
         xhalo=SGdom%nihalo, yhalo=SGdom%njhalo, &
         xextent=exni,yextent=exnj, &
         symmetry=.true., name="MOM_MOSAIC")
  if (SGdom%use_io_layout) &
    call MOM_define_IO_domain(SGdom%mpp_domain, SGdom%io_layout)
!  call mpp_get_compute_domain(G%domain%mpp_domain,isc,iec,jsc,jec)
!  call mpp_get_data_domain(G%domain%mpp_domain,isd,ied,jsd,jed)
!  call mpp_get_compute_domain(SGdom%mpp_domain,isc,iec,jsc,jec)
!  call mpp_get_data_domain(SGdom%mpp_domain,isd,ied,jsd,jed)
  deallocate(exni)
  deallocate(exnj)

  isc=2*G%isc-1 ; iec=2*G%iec ; jsc=2*G%jsc-1 ; jec=2*G%jec
  isd=2*G%isd-1 ; ied=2*G%ied ; jsd=2*G%jsd-1 ; jed=2*G%jed

! Read X from the supergrid
  tmpZ(:,:)=999.
  call read_data(filename,'x',tmpZ,domain=SGdom%mpp_domain,position=CORNER)
 !call read_LRG_supergrid(filename,dv,x=tmpZ,err=err)
 !if (err.ne.0) &
 !  call MOM_error(FATAL," set_grid_metrics_from_mosaic: read_LRG(x) failed!")

  call pass_var(tmpZ, SGdom, position=CORNER)
  call extrapolate_metric(tmpZ,jsc-jsd+1)
  G%geolonh(G%isd:G%ied,G%jsd:G%jed)=tmpZ(isd:ied:2,jsd:jed:2)
  G%geolonq(G%isd:G%ied,G%jsd:G%jed)=tmpZ(isd+1:ied:2,jsd+1:jed:2)
  G%geolonu(G%isd:G%ied,G%jsd:G%jed)=tmpZ(isd+1:ied:2,jsd:jed:2)
  G%geolonv(G%isd:G%ied,G%jsd:G%jed)=tmpZ(isd:ied:2,jsd+1:jed:2)
 ! call pass_var(G%geolonq, G%domain, position=CORNER)

! Read Y from the supergrid
  tmpZ(:,:)=999.
  call read_data(filename,'y',tmpZ,domain=SGdom%mpp_domain,position=CORNER)
 !call read_LRG_supergrid(filename,dv,y=tmpZ,err=err)
 !if (err.ne.0) &
 !  call MOM_error(FATAL," set_grid_metrics_from_mosaic: read_LRG(y) failed!")

  call pass_var(tmpZ, SGdom, position=CORNER)
  call extrapolate_metric(tmpZ,jsc-jsd+1)
  G%geolath(G%isd:G%ied,G%jsd:G%jed)=tmpZ(isd:ied:2,jsd:jed:2)
  G%geolatq(G%isd:G%ied,G%jsd:G%jed)=tmpZ(isd+1:ied:2,jsd+1:jed:2)
  G%geolatu(G%isd:G%ied,G%jsd:G%jed)=tmpZ(isd+1:ied:2,jsd:jed:2)
  G%geolatv(G%isd:G%ied,G%jsd:G%jed)=tmpZ(isd:ied:2,jsd+1:jed:2)

! Read DX,DY from the supergrid
  tmpU(:,:)=0.; tmpV(:,:)=0.
  call read_data(filename,'dx',tmpV,domain=SGdom%mpp_domain,position=NORTH_FACE)
  call read_data(filename,'dy',tmpU,domain=SGdom%mpp_domain,position=EAST_FACE)
 !call read_LRG_supergrid(filename,dv,dx=tmpV,err=err)
 !if (err.ne.0) &
 !  call MOM_error(FATAL," set_grid_metrics_from_mosaic: read_LRG(dx) failed!")
 !call read_LRG_supergrid(filename,dv,dy=tmpU,err=err)
 !if (err.ne.0) &
 !  call MOM_error(FATAL," set_grid_metrics_from_mosaic: read_LRG(dy) failed!")
  call pass_vector(tmpU,tmpV,SGdom,To_All+Scalar_Pair,CGRID_NE)
  call extrapolate_metric(tmpV,jsc-jsd+1)
  call extrapolate_metric(tmpU,jsc-jsd+1)
  dxh(G%isd:G%ied,G%jsd:G%jed)=tmpV(isd:ied:2,jsd:jed:2)+tmpV(isd+1:ied:2,jsd:jed:2)
  dyh(G%isd:G%ied,G%jsd:G%jed)=tmpU(isd:ied:2,jsd:jed:2)+tmpU(isd:ied:2,jsd+1:jed:2)

  dxv(G%isd:G%ied,G%jsd:G%jed)=tmpV(isd:ied:2,jsd+1:jed:2)+tmpV(isd+1:ied:2,jsd+1:jed:2)
  dyu(G%isd:G%ied,G%jsd:G%jed)=tmpU(isd+1:ied:2,jsd:jed:2)+tmpU(isd+1:ied:2,jsd+1:jed:2)

  dxu(G%isd:G%ied-1,G%jsd:G%jed)=tmpV(isd+1:ied-2:2,jsd:jed:2)+tmpV(isd+2:ied-1:2,jsd:jed:2)
  dyv(G%isd:G%ied,G%jsd:G%jed-1)=tmpU(isd:ied:2,jsd+1:jed-2:2)+tmpU(isd:ied:2,jsd+2:jed-1:2)

  dxq(G%isd:G%ied-1,G%jsd:G%jed-1)=tmpV(isd+1:ied-2:2,jsd+1:jed-2:2)+tmpV(isd+2:ied-1:2,jsd+1:jed-2:2)
  dyq(G%isd:G%ied-1,G%jsd:G%jed-1)=tmpU(isd+1:ied-2:2,jsd+1:jed-2:2)+tmpU(isd+1:ied-2:2,jsd+2:jed-1:2)

! Read AREA from the supergrid
  tmpT(:,:)=0.
  call read_data(filename,'area',tmpT,domain=SGdom%mpp_domain)
 !call read_LRG_supergrid(filename,dv,area=tmpT,err=err)
 !if (err.ne.0) &
 !  call MOM_error(FATAL," set_grid_metrics_from_mosaic: read_LRG(A) failed!")
  call pass_var(tmpT, SGdom)
  call extrapolate_metric(tmpT,jsc-jsd+1)

  areaT(G%isd:G%ied,G%jsd:G%jed)=( &
        (tmpT(isd:ied-1:2,jsd:jed-1:2)+tmpT(isd+1:ied:2,jsd+1:jed:2)) &
       +(tmpT(isd+1:ied:2,jsd:jed-1:2)+tmpT(isd:ied-1:2,jsd+1:jed:2)) )
  areaBu(G%isd:G%ied-1,G%jsd:G%jed-1)=0.*( &
        (tmpT(isd+1:ied-2:2,jsd+1:jed-2:2)+tmpT(isd+2:ied-1:2,jsd+2:jed-1:2)) &
       +(tmpT(isd+1:ied-2:2,jsd+2:jed-1:2)+tmpT(isd+2:ied-1:2,jsd+1:jed-2:2)) )

  ni=SGdom%niglobal
  nj=SGdom%njglobal
  deallocate(SGdom%mpp_domain)

  call pass_vector(dyu, dxv, G%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dxu, dyv, G%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dxq, dyq, G%Domain, To_All+Scalar_Pair, BGRID_NE)
  call pass_var(areaT, G%Domain)
  call pass_var(areaBu, G%Domain, position=CORNER)

  do i=G%isd,G%ied ; do j=G%jsd,G%jed
    G%DXh(i,j) = dxh(i,j) ; G%DYh(i,j) = dyh(i,j) ; G%areaT(i,j) = areaT(i,j)
  enddo ; enddo
  do I=G%Isdq,G%Iedq ; do j=G%jsd,G%jed
    G%DXu(I,j) = dxu(I,j) ; G%DYu(I,j) = dyu(I,j)
  enddo ; enddo
  do i=G%isd,G%ied ; do J=G%Jsdq,G%Jedq
    G%DXv(i,J) = dxv(i,J) ; G%DYv(i,J) = dyv(i,J)
  enddo ; enddo
  do I=G%Isdq,G%Iedq ; do J=G%Jsdq,G%Jedq
    G%DXq(I,J) = dxq(I,J) ; G%DYq(I,J) = dyq(I,J) ; G%areaBu(I,J) = areaBu(I,J)
  enddo ; enddo

  ! Construct axes for diagnostic output (only necessary because "ferret" uses
  ! broken convention for interpretting netCDF files).
  start(:) = 1 ; nread(:) = 1
  start(2) = 2 ; nread(1) = ni+1 ; nread(2) = 2
  allocate( tmpGlbl(ni+1,2) )
  if (is_root_PE()) &
    call read_data(filename, "x", tmpGlbl, start, nread, no_domain=.TRUE.)
  call broadcast(tmpGlbl, 2*(ni+1), root_PE())
  
  G%gridlonh(G%domain%nihalo+1:G%domain%nihalo+G%domain%niglobal) = &
    tmpGlbl(2:ni:2,2)
  G%gridlonq(G%domain%nihalo+0:G%domain%nihalo+G%domain%niglobal) = &
    tmpGlbl(1:ni+1:2,1)
  deallocate( tmpGlbl )

  allocate  ( tmpGlbl(1, nj+1) )  
  start(:) = 1 ; nread(:) = 1
  start(1) = int(ni/4)+1 ; nread(2) = nj+1  
  if (is_root_PE()) &
    call read_data(filename, "y", tmpGlbl, start, nread, no_domain=.TRUE.)
  call broadcast(tmpGlbl, nj+1, root_PE())

  G%gridlath(G%domain%njhalo+1:G%domain%njhalo+G%domain%njglobal) = &
    tmpGlbl(1,2:nj:2)
  G%gridlatq(G%domain%njhalo+0:G%domain%njhalo+G%domain%njglobal) = &
    tmpGlbl(1,1:nj+1:2)
  deallocate( tmpGlbl )


end subroutine set_grid_metrics_from_mosaic


! ------------------------------------------------------------------------------

subroutine set_grid_metrics_cartesian(G, param_file)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
! Arguments:
!  (inout)   G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

!    Calculate the values of the metric terms for a Cartesian grid that
!  might be used and save them in arrays.
!    Within this subroutine, the x- and y- grid spacings and their
!  inverses and the cell areas centered on h, q, u, and v points are
!  calculated, as are the geographic locations of each of these 4
!  sets of points.
  integer :: i, j, isd, ied, jsd, jed, Isdq, Iedq, Jsdq, Jedq, X1off, Y1off
  integer :: niglobal, njglobal, nihalo, njhalo
  real :: grid_latq(0:G%Domain%njglobal+2*G%Domain%njhalo)
  real :: grid_lonq(0:G%Domain%niglobal+2*G%Domain%nihalo)
  real :: len_lon, len_lat, west_lon, south_lat, Rad_Earth, PI
  character(len=60) :: axis_units
  character(len=48)  :: mod  = "MOM_grid_init set_grid_metrics_cartesian"

  niglobal = G%Domain%niglobal ; njglobal = G%Domain%njglobal
  nihalo = G%Domain%nihalo ; njhalo = G%Domain%njhalo
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Isdq = G%Isdq ; Iedq = G%Iedq ; Jsdq = G%Jsdq ; Jedq = G%Jedq
  X1off = G%isd_global - isd ; Y1off = G%jsd_global - jsd;

  call MOM_mesg("  MOM_grid_init.F90, set_grid_metrics_cartesian: setting metrics", 5)
 
  PI = 4.0*atan(1.0) ;

  call get_param(param_file, mod, "AXIS_UNITS", axis_units, &
                 "The units for the Cartesian axes. Valid entries are: \n"//&
                 " \t degrees - degrees of latitude and longitude \n"//&
                 " \t m - meters \n \t k - kilometers", default="degrees")
  call get_param(param_file, mod, "SOUTHLAT", south_lat, &
                 "The southern latitude of the domain or the equivalent \n"//&
                 "starting value for the y-axis.", units=axis_units, &
                 fail_if_missing=.true.)
  call get_param(param_file, mod, "LENLAT", len_lat, &
                 "The latitudinal or y-direction length of the domain.", &
                 units=axis_units, fail_if_missing=.true.)
  call get_param(param_file, mod, "WESTLON", west_lon, &
                 "The western longitude of the domain or the equivalent \n"//&
                 "starting value for the x-axis.", units=axis_units, &
                 default=0.0)
  call get_param(param_file, mod, "LENLON", len_lon, &
                 "The longitudinal or x-direction length of the domain.", &
                 units=axis_units, fail_if_missing=.true.)
  call get_param(param_file, mod, "RAD_EARTH", Rad_Earth, &
                 "The radius of the Earth.", units="m", default=6.378e6)

  ! These are larger in case symmetric memory is being used.
  do J=0,njglobal+2*njhalo
    grid_latq(J) = south_lat + len_lat* REAL(J-njhalo)/REAL(njglobal)
  enddo 
  do I=0,niglobal+2*nihalo
    grid_lonq(I) = west_lon + len_lon*REAL(I-nihalo)/REAL(niglobal)
  enddo
  do j=1,njglobal+2*njhalo
    G%gridlatq(j) = south_lat + len_lat* REAL(j-njhalo)/REAL(njglobal)
    G%gridlath(j) = south_lat + len_lat*(REAL(j-njhalo)-0.5)/REAL(njglobal)
  enddo
  do i=1,niglobal+2*nihalo
    G%gridlonq(i) = west_lon + len_lon*REAL(i-nihalo)/REAL(niglobal)
    G%gridlonh(i) = west_lon + len_lon*(REAL(i-nihalo)-0.5)/REAL(niglobal)
  enddo

  do J=Jsdq,Jedq ; do I=Isdq,Iedq
    G%geolonq(i,j) = grid_lonq(i+X1off) ; G%geolatq(i,j) = grid_latq(j+Y1off)

    G%DXq(I,J) = Rad_Earth * len_lon * PI / (180.0 * niglobal)
    G%DYq(I,J) = Rad_Earth * len_lat * PI / (180.0 * njglobal)

    if (axis_units(1:1) == 'k') then ! Axes are measured in km.
      G%DXq(I,J) = 1000.0 * len_lon / (REAL(niglobal))
      G%DYq(I,J) = 1000.0 * len_lat / (REAL(njglobal))
    else if (axis_units(1:1) == 'm') then ! Axes are measured in m.
      G%DXq(I,J) = len_lon / (REAL(niglobal))
      G%DYq(I,J) = len_lat / (REAL(njglobal))
    else ! Axes are measured in degrees of latitude and longitude.
      G%DXq(I,J) = Rad_Earth * len_lon * PI / (180.0 * niglobal)
      G%DYq(I,J) = Rad_Earth * len_lat * PI / (180.0 * njglobal)
    endif

    G%IDXq(I,J) = 1.0 / G%DXq(I,J)
    G%IDYq(I,J) = 1.0 / G%DYq(I,J)
    G%areaBu(I,J) = G%DXq(I,J) * G%DYq(I,J)
    G%IareaBu(I,J) = G%IDXq(I,J) * G%IDYq(I,J)
  enddo ; enddo

  do j=jsd,jed ; do i=isd,ied
    G%geolonh(i,j) = G%gridlonh(i+X1off) ; G%geolath(i,j) = G%gridlath(j+Y1off)
    G%DXh(i,j) = G%DXq(I,J) ; G%IDXh(i,j) = G%IDXq(I,J)
    G%DYh(i,j) = G%DYq(I,J) ; G%IDYh(i,j) = G%IDYq(I,J)
    G%areaT(i,j) = G%areaBu(I,J) ; G%IareaT(i,j) = G%IareaBu(I,J)
  enddo ; enddo

  do j=jsd,jed ; do I=Isdq,Iedq
    G%geolonu(i,j) = grid_lonq(i+X1off) ; G%geolatu(i,j) = G%gridlath(j+Y1off)

    G%DXu(I,j) = G%DXq(I,J) ; G%IDXu(I,j) = G%IDXq(I,J)
    G%DYu(I,j) = G%DYq(I,J) ; G%IDYu(I,j) = G%IDYq(I,J)
  enddo ; enddo

  do J=Jsdq,Jedq ; do i=isd,ied
    G%geolonv(i,j) = G%gridlonh(i+X1off) ; G%geolatv(i,j) = grid_latq(j+Y1off)

    G%DXv(i,J) = G%DXq(I,J) ; G%IDXv(i,J) = G%IDXq(I,J)
    G%DYv(i,J) = G%DYq(I,J) ; G%IDYv(i,J) = G%IDYq(I,J)
  enddo ; enddo

end subroutine set_grid_metrics_cartesian

! ------------------------------------------------------------------------------

subroutine set_grid_metrics_spherical(G, param_file)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
! Arguments:
!  (inout)   G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

!    Calculate the values of the metric terms that might be used
!  and save them in arrays.
!    Within this subroutine, the x- and y- grid spacings and their
!  inverses and the cell areas centered on h, q, u, and v points are
!  calculated, as are the geographic locations of each of these 4
!  sets of points.
  real :: PI, PI_180! PI = 3.1415926... as 4*atan(1)
  integer :: i,j, isd, ied, jsd, jed
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, Isdq, Iedq, Jsdq, Jedq
  character(len=200) :: axis_units
  integer :: i_offset, j_offset
  real :: grid_latq(0:G%Domain%njglobal+2*G%Domain%njhalo)
  real :: grid_lonq(0:G%Domain%niglobal+2*G%Domain%nihalo)
  real :: dLon,dLat,latitude,longitude,dL_di
  real :: south_lat,len_lat,west_lon,len_lon,Rad_Earth
  character(len=48)  :: mod  = "MOM_grid_init set_grid_metrics_spherical"

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq
  Isdq = G%Isdq ; Iedq = G%Iedq ; Jsdq = G%Jsdq ; Jedq = G%Jedq
  i_offset = G%isd_global - isd; j_offset = G%jsd_global - jsd

  call MOM_mesg("  MOM_grid_init.F90, set_grid_metrics_simple_spherical: "// &
                 "Setting metrics.", 5)
 
!    Calculate the values of the metric terms that might be used
!  and save them in arrays.
  PI = 4.0*atan(1.0); PI_180 = atan(1.0)/45.
 
  call get_param(param_file, mod, "AXIS_UNITS", axis_units, default="degrees")
  if (trim(axis_units) == "") axis_units = "degrees"
  if (trim(axis_units) .ne. "degrees") call MOM_error(FATAL, &
    "MOM_grid_init.F90, set_grid_metrics_simple_spherical: "// &
    "axis_units must be degrees")

  call get_param(param_file, mod, "SOUTHLAT", south_lat, &
                 "The southern latitude of the domain.", units="degrees", &
                 fail_if_missing=.true.)
  call get_param(param_file, mod, "LENLAT", len_lat, &
                 "The latitudinal length of the domain.", units="degrees", &
                 fail_if_missing=.true.)
  call get_param(param_file, mod, "WESTLON", west_lon, &
                 "The western longitude of the domain.", units="degrees", &
                 default=0.0)
  call get_param(param_file, mod, "LENLON", len_lon, &
                 "The longitudinal length of the domain.", units="degrees", &
                 fail_if_missing=.true.)
  call get_param(param_file, mod, "RAD_EARTH", Rad_Earth, &
                 "The radius of the Earth.", units="m", default=6.378e6)

  dLon = len_lon/G%Domain%niglobal
  dLat = len_lat/G%Domain%njglobal

  do J=0,G%Domain%njglobal+2*G%Domain%njhalo
    latitude = south_lat + dLat* REAL(J-G%Domain%njhalo)
    grid_latq(J) = MIN(MAX(latitude,-90.),90.)
  enddo
  do I=1,G%Domain%niglobal+2*G%Domain%nihalo
    grid_lonq(I) = west_lon + dLon*REAL(I-G%Domain%nihalo)
  enddo

  do j=1,G%Domain%njglobal+2*G%Domain%njhalo
    G%gridlatq(J) = grid_latq(J)
    latitude = south_lat + dLat*(REAL(j-G%Domain%njhalo)-0.5)
    G%gridlath(j) = MIN(MAX(latitude,-90.),90.)
  enddo
  do i=1,G%Domain%niglobal+2*G%Domain%nihalo
    G%gridlonq(I) = grid_lonq(I)
    G%gridlonh(i) = west_lon + dLon*(REAL(i-G%Domain%nihalo)-0.5)
  enddo

  dL_di = (len_lon * 4.0*atan(1.0)) / (180.0 * G%Domain%niglobal)
  do J=Jsdq,Jedq ; do I=Isdq,Iedq
    G%geolonq(I,J) = grid_lonq(I+I_offset)
    G%geolatq(I,J) = grid_latq(J+J_offset)

! The following line is needed to reproduce the solution from
! set_grid_metrics_mercator when used to generate a simple spherical grid.
    G%DXq(I,J) = Rad_Earth * COS( G%geolatq(I,J)*PI_180 ) * dL_di
!   G%DXq(I,J) = Rad_Earth * dLon*PI_180 * COS( G%geolatq(I,J)*PI_180 )
    G%DYq(I,J) = Rad_Earth * dLat*PI_180
    G%areaBu(I,J) = G%DXq(I,J) * G%DYq(I,J)
  enddo; enddo

  do J=Jsdq,Jedq ; do i=isd,ied
    G%geolonv(i,J) = G%gridlonh(i+i_offset)
    G%geolatv(i,J) = grid_latq(j+j_offset)

! The following line is needed to reproduce the solution from
! set_grid_metrics_mercator when used to generate a simple spherical grid.
    G%DXv(i,J) = Rad_Earth * COS( G%geolatv(i,J)*PI_180 ) * dL_di
!   G%DXv(i,J) = Rad_Earth * (dLon*PI_180) * COS( G%geolatv(i,J)*PI_180 )
    G%DYv(i,J) = Rad_Earth * dLat*PI_180
  enddo; enddo

  do j=jsd,jed ; do I=Isdq,Iedq
    G%geolonu(I,j) = grid_lonq(i+i_offset)
    G%geolatu(I,j) = G%gridlath(j+j_offset)

! The following line is needed to reproduce the solution from
! set_grid_metrics_mercator when used to generate a simple spherical grid.
    G%DXu(I,j) = Rad_Earth * COS( G%geolatu(I,j)*PI_180 ) * dL_di
!   G%DXu(I,j) = Rad_Earth * dLon*PI_180 * COS( latitude )
    G%DYu(I,j) = Rad_Earth * dLat*PI_180
  enddo; enddo

  do j=jsd,jed ; do i=isd,ied
    G%geolonh(i,j) = G%gridlonh(i+i_offset)
    G%geolath(i,j) = G%gridlath(j+j_offset)

! The following line is needed to reproduce the solution from
! set_grid_metrics_mercator when used to generate a simple spherical grid.
    G%DXh(i,j) = Rad_Earth * COS( G%geolath(i,j)*PI_180 ) * dL_di
!   G%DXh(i,j) = Rad_Earth * dLon*PI_180 * COS( latitude )
    G%DYh(i,j) = Rad_Earth * dLat*PI_180

!   latitude = G%geolatv(i,J)*PI_180             ! In radians
!   dL_di    = G%geolatv(i,max(jsd,J-1))*PI_180  ! In radians
!   G%areaT(i,j) = Rad_Earth**2*dLon*dLat*ABS(SIN(latitude)-SIN(dL_di))
    G%areaT(i,j) = G%DXh(i,j) * G%DYh(i,j)
  enddo; enddo

end subroutine set_grid_metrics_spherical

! ------------------------------------------------------------------------------

subroutine set_grid_metrics_mercator(G, param_file)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
! Arguments:
!  (inout)   G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

!    Calculate the values of the metric terms that might be used
!  and save them in arrays.
!    Within this subroutine, the x- and y- grid spacings and their
!  inverses and the cell areas centered on h, q, u, and v points are
!  calculated, as are the geographic locations of each of these 4
!  sets of points.
  real :: grid_latq(0:G%Domain%njglobal+2*G%Domain%njhalo)
  real :: grid_lonq(0:G%Domain%niglobal+2*G%Domain%nihalo)
  integer :: i,j, isd, ied, jsd, jed
  integer :: X1off, Y1off
  type(GPS) :: GP
  character(len=128) :: warnmesg
  character(len=48)  :: mod = "MOM_grid_init set_grid_metrics_mercator"
  real :: PI, PI_2! PI = 3.1415926... as 4*atan(1), PI_2 = (PI) /2.0


!   All of the metric terms should be defined over the domain from
! isd to ied.  Outside of the physical domain, both the metrics
! and their inverses may be set to zero.

!  The metric terms within the computational domain are set here.
  real :: y_q, y_h, jd, x_q, x_h, id
  real, dimension(G%isd:G%ied,G%jsd:G%jed) :: &
    xh, yh ! Latitude and longitude of h points in radians.
  real, dimension(G%Isdq:G%Iedq,G%jsd:G%jed) :: &
    xu, yu ! Latitude and longitude of u points in radians.
  real, dimension(G%isd:G%ied,G%Jsdq:G%Jedq) :: &
    xv, yv ! Latitude and longitude of v points in radians.
  real, dimension(G%Isdq:G%Iedq,G%Jsdq:G%Jedq) :: &
    xq, yq ! Latitude and longitude of q points in radians.
  real :: fnRef           ! fnRef is the value of Int_dj_dy or
                          ! Int_dj_dy at a latitude or longitude that is
  real :: jRef, iRef      ! being set to be at grid index jRef or iRef.
  real :: dx_q, x_q_west
  integer :: itt1, itt2
  integer :: err
  logical :: debug = .FALSE., simple_area = .true.
  real    :: temp(G%isdq:G%iedq,G%jsdq:G%jedq)
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, Isdq, Iedq, Jsdq, Jedq

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Isq = G%Iscq ; Ieq = G%Iecq ; Jsq = G%Jscq ; Jeq = G%Jecq
  Isdq = G%Isdq ; Iedq = G%Iedq ; Jsdq = G%Jsdq ; Jedq = G%Jedq
  X1off = G%isd_global - isd ; Y1off = G%jsd_global - jsd;

  GP%niglobal = G%Domain%niglobal
  GP%njglobal = G%Domain%njglobal

  call MOM_mesg("  MOM_grid_init.F90, set_grid_metrics_mercator: setting metrics", 5)
 
!    Calculate the values of the metric terms that might be used
!  and save them in arrays.
  PI = 4.0*atan(1.0) ; PI_2 = 0.5*PI

  call get_param(param_file, mod, "SOUTHLAT", GP%south_lat, &
                 "The southern latitude of the domain.", units="degrees", &
                 fail_if_missing=.true.)
  call get_param(param_file, mod, "LENLAT", GP%len_lat, &
                 "The latitudinal length of the domain.", units="degrees", &
                 fail_if_missing=.true.)
  call get_param(param_file, mod, "WESTLON", GP%west_lon, &
                 "The western longitude of the domain.", units="degrees", &
                 default=0.0)
  call get_param(param_file, mod, "LENLON", GP%len_lon, &
                 "The longitudinal length of the domain.", units="degrees", &
                 fail_if_missing=.true.)
  call get_param(param_file, mod, "RAD_EARTH", GP%Rad_Earth, &
                 "The radius of the Earth.", units="m", default=6.378e6)
  call get_param(param_file, mod, "ISOTROPIC", GP%isotropic, &
                 "If true, an isotropic grid on a sphere (also known as \n"//&
                 "a Mercator grid) is used. With an isotropic grid, the \n"//&
                 "meridional extent of the domain (LENLAT), the zonal \n"//&
                 "extent (LENLON), and the number of grid points in each \n"//&
                 "direction are _not_ independent. In MOM the meridional \n"//&
                 "extent is determined to fit the zonal extent and the \n"//&
                 "number of grid points, while grid is perfectly isotropic.", &
                 default=.false.)
  call get_param(param_file, mod, "EQUATOR_REFERENCE", GP%equator_reference, &
                 "If true, the grid is defined to have the equator at the \n"//&
                 "nearest q- or h- grid point to (-LOWLAT*NJGLOBAL/LENLAT).", &
                 default=.true.)
  call get_param(param_file, mod, "LAT_ENHANCE_FACTOR", GP%Lat_enhance_factor, &
                 "The amount by which the meridional resolution is \n"//&
                 "enhanced within LAT_EQ_ENHANCE of the equator.", &
                 units="nondim", default=1.0)
  call get_param(param_file, mod, "LAT_EQ_ENHANCE", GP%Lat_eq_enhance, &
                 "The latitude range to the north and south of the equator \n"//&
                 "over which the resolution is enhanced.", units="degrees", &
                 default=0.0)

!    With an isotropic grid, the north-south extent of the domain,
!  the east-west extent, and the number of grid points in each
!  direction are _not_ independent.  Here the north-south extent
!  will be determined to fit the east-west extent and the number of
!  grid points.  The grid is perfectly isotropic.
  if (GP%equator_reference) then
! With the following expression, the equator will always be placed
! on either h or q points, in a position consistent with the ratio
! GP%south_lat to GP%len_lat.
    jRef =  G%Domain%njhalo-1 + 0.5*FLOOR(GP%njglobal*((-1.0*GP%south_lat*2.0)/GP%len_lat)+0.5)
    fnRef = Int_dj_dy(0.0,GP)
  else
!  The following line sets the refererence latitude GP%south_lat at j=js-1.
    jRef = G%Domain%njhalo-1
    fnRef = Int_dj_dy((GP%south_lat*PI/180.0),GP)
  endif

  y_q = GP%south_lat*PI/180.0
  if ((0 >= Jsdq+Y1off) .and. (0 <= Jedq+Y1off)) then
    do I=Isdq,Iedq ; yq(I, -Y1off) = y_q ; enddo
    do i=isd,ied ; yv(i, -Y1off) = y_q ; enddo
  endif
  do j=1,GP%njglobal+2*G%Domain%njhalo
    jd = fnRef + (j -1 - jRef) - 0.5
    y_h = find_root(Int_dj_dy,dy_dj,GP,jd,y_q,-1.0*PI_2,PI_2,itt1)

    jd = fnRef + (j -1 - jRef)
    y_q = find_root(Int_dj_dy,dy_dj,GP,jd,y_h,-1.0*PI_2,PI_2,itt2)

    G%gridlatq(j) = y_q*180.0/PI
    G%gridlath(j) = y_h*180.0/PI

    if ((j >= jsd+Y1off) .and. (j <= jed+Y1off)) then
      do i=isd,ied ; yh(i,j-Y1off) = y_h ; enddo
      do I=Isdq,Iedq ; yu(I,j-Y1off) = y_h ; enddo
    endif
    if ((J >= Jsdq+Y1off) .and. (J <= Jedq+Y1off)) then
      do I=Isdq,Iedq ; yq(I,J-Y1off) = y_q ; enddo
      do i=isd,ied ; yv(i,J-Y1off) = y_q ; enddo
    endif
  enddo

! Determine the longitudes of the various points.

! These two lines place the western edge of the domain at GP%west_lon.
  iRef = G%Domain%nihalo -1 + GP%niglobal
  fnRef = Int_di_dx(((GP%west_lon+GP%len_lon)*PI/180.0),GP)

  x_q = GP%west_lon*PI/180.0
  x_q_west = x_q
! If the model is in parallel in the X-direction, do the same set of
! calculations which would occur on a single processor.
  if ((0 >= Isdq+X1off) .and. (0 <= Iedq+X1off)) then
    do J=Jsdq,Jedq ; xq(-X1off,j) = x_q ; enddo
    do j=jsd,jed ; xu(-X1off,j) = x_q ; enddo
  endif
  do i=1,GP%niglobal+2*G%Domain%nihalo
    id = fnRef + (i - 1 - iRef) - 0.5
    x_h = find_root(Int_di_dx,dx_di,GP,id,x_q,-4.0*PI,4.0*PI,itt1)

    id = fnRef + (i - 1 -iRef)
    x_q = find_root(Int_di_dx,dx_di,GP,id,x_h,-4.0*PI,4.0*PI,itt2)
    if(i == G%isc) dx_q = x_q - x_q_west

    G%gridlonq(i) = x_q*180.0/PI
    G%gridlonh(i) = x_h*180.0/PI

    if ((i >= isd+X1off) .and. (i <= ied+X1off)) then
      do j=jsd,jed ; xh(i-X1off,j) = x_h ; enddo
      do J=Jsdq,Jedq ; xv(i-X1off,J) = x_h ; enddo
    endif
    if ((I >= Isdq+X1off) .and. (I <= Iedq+X1off)) then
      do J=Jsdq,Jedq ; xq(I-X1off,J) = x_q ; enddo
      do j=jsd,jed ; xu(I-X1off,j) = x_q ; enddo
    endif
  enddo

  do J=Jsdq,Jedq ; do I=Isdq,Iedq
    G%geolonq(i,j) = xq(i,j)*180.0/PI
    G%geolatq(i,j) = yq(i,j)*180.0/PI
    G%DXq(i,j) = ds_di(xq(i,j), yq(i,j), GP)
    G%DYq(i,j) = ds_dj(xq(i,j), yq(i,j), GP)

    G%areaBu(i,j) = G%DXq(i,j) * G%DYq(i,j)
    G%IareaBu(i,j) = 1.0 / G%areaBu(i,j)
  enddo ; enddo

  do j=jsd,jed ; do i=isd,ied
    G%geolonh(i,j) = xh(i,j)*180.0/PI
    G%geolath(i,j) = yh(i,j)*180.0/PI
    G%DXh(i,j) = ds_di(xh(i,j), yh(i,j), GP)
    G%DYh(i,j) = ds_dj(xh(i,j), yh(i,j), GP)

    G%areaT(i,j) = G%DXh(i,j)*G%DYh(i,j)
    G%IareaT(i,j) = 1.0 / G%areaT(i,j)
  enddo ; enddo

  do j=jsd,jed ; do I=Isdq,Iedq
    G%geolonu(i,j) = xu(i,j)*180.0/PI
    G%geolatu(i,j) = yu(i,j)*180.0/PI
    G%DXu(i,j) = ds_di(xu(i,j), yu(i,j), GP)
    G%DYu(i,j) = ds_dj(xu(i,j), yu(i,j), GP)
  enddo ; enddo

  do J=Jsdq,Jedq ; do i=isd,ied
    G%geolonv(i,j) = xv(i,j)*180.0/PI
    G%geolatv(i,j) = yv(i,j)*180.0/PI
    G%DXv(i,j) = ds_di(xv(i,j), yv(i,j), GP)
    G%DYv(i,j) = ds_dj(xv(i,j), yv(i,j), GP)
  enddo ; enddo

  if (.not.simple_area) then
    do j=Jsdq+1,jed ; do i=Isdq+1,ied
  !         The following test is to ensure parallel reproducibility.
      if (size(G%areaT(:,:),1) == 1) then !{
        temp(i,j) = GP%Rad_Earth**2 * &
            (dL(0.0,0.0,yq(I-1,J-1),yq(I-1,J)) + &
            (dL(0.0,dx_q,yq(I-1,J),yq(I,J)) +   &
            (dL(0.0,0.0,yq(I,J),yq(I,J-1)) +      &
             dL(dx_q,0.0,yq(I,J-1),yq(I-1,J-1)))))
      else
        temp(I,J) = GP%Rad_Earth**2 * &
            (dL(xq(I-1,J-1),xq(I-1,J),yq(I-1,J-1),yq(I-1,J)) + &
            (dL(xq(I-1,J),xq(I,J),yq(I-1,J),yq(I,J)) +          &
            (dL(xq(I,J),xq(I,J-1),yq(I,J),yq(I,J-1)) +          &
             dL(xq(I,J-1),xq(I-1,J-1),yq(I,J-1),yq(I-1,J-1)))))
      endif
    enddo ;enddo
  ! Fill in row and column one.
    if ((Isdq == isd) .or. (Jsdq == jsq)) then
      temp(isd,:) = temp(isd+1,jsd+1) ; temp(:,jsd) = temp(isd+1,jsd+1)
      temp(isd,:) = temp(isd+1,:)     ; temp(:,jsd) = temp(:,jsd+1)
      call pass_var(temp,G%Domain)
    endif
    do j=jsd,jed ; do i=isd,ied
      G%areaT(i,j) = temp(i,j)
      G%IareaT(i,j) = 1.0 / G%areaT(i,j)
    enddo ; enddo
  endif

end subroutine set_grid_metrics_mercator

! ------------------------------------------------------------------------------

subroutine extrapolate_metric(var, jh)
  real, dimension(:,:), intent(inout) ::  var
  integer, intent(in) :: jh
  integer :: i,j

  ! Fill in southern halo by extrapolating from the computational domain
  do j=lbound(var,2)+jh,lbound(var,2),-1; do i=lbound(var,1),ubound(var,1)
    if (var(i,j)==0.) var(i,j)=2.0*var(i,j+1)-var(i,j+2)
  enddo; enddo

  ! Fill in northern halo by extrapolating from the computational domain
  do j=ubound(var,2)-jh,ubound(var,2); do i=lbound(var,1),ubound(var,1)
    if (var(i,j)==0.) var(i,j)=2.0*var(i,j-1)-var(i,j-2)
  enddo; enddo
  
  ! Fill in western halo by extrapolating from the computational domain
  do j=lbound(var,2),ubound(var,2); do i=lbound(var,1)+jh,lbound(var,1),-1
    if (var(i,j)==0.) var(i,j)=2.0*var(i+1,j)-var(i+2,j)
  enddo; enddo

  ! Fill in eastern halo by extrapolating from the computational domain
  do j=lbound(var,2),ubound(var,2); do i=ubound(var,1)-jh,ubound(var,1)
    if (var(i,j)==0.) var(i,j)=2.0*var(i-1,j)-var(i-2,j)
  enddo; enddo

end subroutine extrapolate_metric


subroutine allocate_metrics(G)
  type(ocean_grid_type), intent(inout) :: G
  integer :: isd, ied, jsd, jed, Isdq, Iedq, Jsdq, Jedq

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Isdq = G%Isdq ; Iedq = G%Iedq ; Jsdq = G%Jsdq ; Jedq = G%Jedq

  ALLOC_(G%dxh(isd:ied,jsd:jed)) ; G%dxh(:,:) = 0.0
  ALLOC_(G%dxu(Isdq:Iedq,jsd:jed)) ; G%dxu(:,:) = 0.0
  ALLOC_(G%dxv(isd:ied,Jsdq:Jedq)) ; G%dxv(:,:) = 0.0
  ALLOC_(G%dxq(Isdq:Iedq,Jsdq:Jedq)) ; G%dxq(:,:) = 0.0
  ALLOC_(G%Idxh(isd:ied,jsd:jed)) ; G%Idxh(:,:) = 0.0
  ALLOC_(G%Idxu(Isdq:Iedq,jsd:jed)) ; G%Idxu(:,:) = 0.0
  ALLOC_(G%Idxv(isd:ied,Jsdq:Jedq)) ; G%Idxv(:,:) = 0.0
  ALLOC_(G%Idxq(Isdq:Iedq,Jsdq:Jedq)) ; G%Idxq(:,:) = 0.0

  ALLOC_(G%dyh(isd:ied,jsd:jed)) ; G%dyh(:,:) = 0.0
  ALLOC_(G%dyu(Isdq:Iedq,jsd:jed)) ; G%dyu(:,:) = 0.0
  ALLOC_(G%dyv(isd:ied,Jsdq:Jedq)) ; G%dyv(:,:) = 0.0
  ALLOC_(G%dyq(Isdq:Iedq,Jsdq:Jedq)) ; G%dyq(:,:) = 0.0
  ALLOC_(G%Idyh(isd:ied,jsd:jed)) ; G%Idyh(:,:) = 0.0
  ALLOC_(G%Idyu(Isdq:Iedq,jsd:jed)) ; G%Idyu(:,:) = 0.0
  ALLOC_(G%Idyv(isd:ied,Jsdq:Jedq)) ; G%Idyv(:,:) = 0.0
  ALLOC_(G%Idyq(Isdq:Iedq,Jsdq:Jedq)) ; G%Idyq(:,:) = 0.0

  ALLOC_(G%areaT(isd:ied,jsd:jed)) ; G%areaT(:,:) = 0.0
  ALLOC_(G%IareaT(isd:ied,jsd:jed)) ; G%IareaT(:,:) = 0.0
  ALLOC_(G%areaBu(Isdq:Iedq,Jsdq:Jedq)) ; G%areaBu(:,:) = 0.0
  ALLOC_(G%IareaBu(Isdq:Iedq,Jsdq:Jedq)) ; G%IareaBu(:,:) = 0.0

  ALLOC_(G%hmask(isd:ied,jsd:jed)) ; G%hmask(:,:) = 0.0
  ALLOC_(G%umask(Isdq:Iedq,jsd:jed)) ; G%umask(:,:) = 0.0
  ALLOC_(G%vmask(isd:ied,Jsdq:Jedq)) ; G%vmask(:,:) = 0.0
  ALLOC_(G%qmask(Isdq:Iedq,Jsdq:Jedq)) ; G%qmask(:,:) = 0.0
  ALLOC_(G%geolath(isd:ied,jsd:jed)) ; G%geolath(:,:) = 0.0
  ALLOC_(G%geolatu(Isdq:Iedq,jsd:jed)) ; G%geolatu(:,:) = 0.0
  ALLOC_(G%geolatv(isd:ied,Jsdq:Jedq)) ; G%geolatv(:,:) = 0.0
  ALLOC_(G%geolatq(Isdq:Iedq,Jsdq:Jedq)) ; G%geolatq(:,:) = 0.0
  ALLOC_(G%geolonh(isd:ied,jsd:jed)) ; G%geolonh(:,:) = 0.0
  ALLOC_(G%geolonu(Isdq:Iedq,jsd:jed)) ; G%geolonu(:,:) = 0.0
  ALLOC_(G%geolonv(isd:ied,Jsdq:Jedq)) ; G%geolonv(:,:) = 0.0
  ALLOC_(G%geolonq(Isdq:Iedq,Jsdq:Jedq)) ; G%geolonq(:,:) = 0.0

  ALLOC_(G%dx_v(isd:ied,Jsdq:Jedq)) ; G%dx_v(:,:) = 0.0
  ALLOC_(G%dy_u(Isdq:Iedq,jsd:jed)) ; G%dy_u(:,:) = 0.0
  ALLOC_(G%dx_v_obc(isd:ied,Jsdq:Jedq)) ; G%dx_v_obc(:,:) = 0.0
  ALLOC_(G%dy_u_obc(Isdq:Iedq,jsd:jed)) ; G%dy_u_obc(:,:) = 0.0  
  ALLOC_(G%areaCu(Isdq:Iedq,jsd:jed)) ; G%areaCu(:,:) = 0.0
  ALLOC_(G%areaCv(isd:ied,Jsdq:Jedq)) ; G%areaCv(:,:) = 0.0
  ALLOC_(G%IareaCu(Isdq:Iedq,jsd:jed)) ; G%IareaCu(:,:) = 0.0
  ALLOC_(G%IareaCv(isd:ied,Jsdq:Jedq)) ; G%IareaCv(:,:) = 0.0

end subroutine allocate_metrics


subroutine initialize_masks(G, PF)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: PF
! Arguments:
!  (inout)   G - The ocean's grid structure.
!  (in)      PF - A structure indicating the open file to parse for
!                 model parameter values.

!    Initialize_masks sets hmask, umask, vmask, and qmask to mask out
! flow over any points which are shallower than Dmin and permit an
! appropriate treatment of the boundary conditions.  umask and vmask
! are 0.0 at any points adjacent to a land point.  qmask is 0.0 at
! any land or boundary point.  For points in the interior, umask,
! vmask, and qmask are all 1.0.

  real :: Dmin, min_depth
  integer :: i, j
  integer :: isd, isd_global, jsd, jsd_global 
  logical :: apply_OBC_u_flather_east, apply_OBC_u_flather_west
  logical :: apply_OBC_v_flather_north, apply_OBC_v_flather_south
  character(len=40)  :: mod = "MOM_grid_init initialize_masks"

  call get_param(PF, "MOM_grid_init initialize_masks", &
                 "MINIMUM_DEPTH", min_depth, &
                 "The minimum ocean depth, anything shallower than which \n"//&
                 "is assumed to be land and all fluxes are masked out.", &
                 units="m", default=0.0)
  call get_param(PF, mod, "APPLY_OBC_U_FLATHER_EAST", apply_OBC_u_flather_east,&
                 "Apply a Flather open boundary condition on the eastern \n"//&
                 "side of the global domain", default=.false.)
  call get_param(PF, mod, "APPLY_OBC_U_FLATHER_WEST", apply_OBC_u_flather_west,&
                 "Apply a Flather open boundary condition on the western \n"//&
                 "side of the global domain", default=.false.)
  call get_param(PF, mod, "APPLY_OBC_V_FLATHER_NORTH", apply_OBC_v_flather_north,&
                 "Apply a Flather open boundary condition on the northern \n"//&
                 "side of the global domain", default=.false.)
  call get_param(PF, mod, "APPLY_OBC_V_FLATHER_SOUTH", apply_OBC_v_flather_south,&
                 "Apply a Flather open boundary condition on the southern \n"//&
                 "side of the global domain", default=.false.)

  if ((apply_OBC_u_flather_west .or. apply_OBC_v_flather_south) .and. &
      .not.G%symmetric ) &
    call MOM_error(FATAL, "Symmetric memory must be used when "//&
      "APPLY_OBC_U_FLATHER_WEST or APPLY_OBC_V_FLATHER_SOUTH is true.")

  Dmin = MAX(min_depth,2.0*G%Angstrom_z)

  call pass_var(G%bathyT, G%Domain)
  G%umask(:,:) = 0.0 ; G%vmask(:,:) = 0.0 ; G%qmask(:,:) = 0.0

  ! Extrapolate the bottom depths at any points that are subject to Flather
  ! open boundary conditions.  This should be generalized for Flather OBCs
  ! that are not necessarily at the edges of the domain.
  if (apply_OBC_u_flather_west) then
    do j=G%jsd,G%jed ; do I=G%isd+1,G%ied
      if ((I+G%isd_global-G%isd) == G%domain%nihalo+1) then
        G%bathyT(i-1,j) = G%bathyT(i,j)
      endif
    enddo; enddo       
  endif

  if (apply_OBC_u_flather_east) then
    do j=G%jsd,G%jed ; do I=G%isd,G%ied-1
      if ((i+G%isd_global-G%isd) == G%domain%niglobal+G%domain%nihalo) then
        G%bathyT(i+1,j) = G%bathyT(i,j)
      endif
    enddo; enddo    
  endif

  if (apply_OBC_v_flather_north) then
    do J=G%jsd,G%jed-1 ; do i=G%isd,G%ied
      if ((j+G%jsd_global-G%jsd) == G%domain%njglobal+G%domain%njhalo) then
        G%bathyT(i,j+1) = G%bathyT(i,j)
      endif
    enddo; enddo    
  endif

  if (apply_OBC_v_flather_south) then
    do J=G%jsd+1,G%jed ; do i=G%isd,G%ied
      if ((J+G%jsd_global-G%jsd) == G%domain%njhalo+1) then
        G%bathyT(i,j-1) = G%bathyT(i,j)
      endif
    enddo; enddo
  endif

  do j=G%jsd,G%jed ; do i=G%isd,G%ied
    if (G%bathyT(i,j) <= Dmin) then
      G%hmask(i,j) = 0.0
    else
      G%hmask(i,j) = 1.0
    endif
  enddo ; enddo

  do j=G%jsd,G%jed ; do I=G%isd,G%ied-1
    if ((G%bathyT(i,j) <= Dmin) .or. (G%bathyT(i+1,j) <= Dmin)) then
      G%umask(I,j) = 0.0
    else
      G%umask(I,j) = 1.0
    endif
  enddo ; enddo

  do J=G%jsd,G%jed-1 ; do i=G%isd,G%ied
    if ((G%bathyT(i,j) <= Dmin) .or. (G%bathyT(i,j+1) <= Dmin)) then
      G%vmask(i,J) = 0.0
    else
      G%vmask(i,J) = 1.0
    endif
  enddo ; enddo

  do J=G%jsd,G%jed-1 ; do I=G%isd,G%ied-1
    if ((G%bathyT(i+1,j) <= Dmin) .or. (G%bathyT(i+1,j+1) <= Dmin) .or. &
        (G%bathyT(i,j) <= Dmin) .or. (G%bathyT(i,j+1) <= Dmin)) then
      G%qmask(I,J) = 0.0
    else
      G%qmask(I,J) = 1.0
    endif
  enddo ; enddo

  call pass_vector(G%umask, G%vmask, G%Domain, To_All+Scalar_Pair, CGRID_NE)

  do j=G%Jsdq,G%Jedq ; do i=G%isd,G%ied
    G%dx_v(i,j) = G%vmask(i,j)*G%DXv(i,j)
    G%dx_v_obc(i,j) = G%vmask(i,j)*G%DXv(i,j)
    G%areaCv(i,j) = G%DYv(i,j)*G%dx_v(i,j)
    G%IareaCv(i,j) = 0.0
    if (G%areaCv(i,j) > 0.0) G%IareaCv(i,j) = G%vmask(i,j) / G%areaCv(i,j)
  enddo ; enddo
  do j=G%jsd,G%jed ; do i=G%Isdq,G%Iedq  
    G%dy_u(i,j) = G%umask(i,j)*G%DYu(i,j)
    G%dy_u_obc(i,j) = G%umask(i,j)*G%DYu(i,j)  
    G%areaCu(i,j) = G%DXu(i,j)*G%dy_u(i,j)
    G%IareaCu(i,j) = 0.0
    if (G%areaCu(i,j) > 0.0) G%IareaCu(i,j) = G%umask(i,j) / G%areaCu(i,j)
  enddo ; enddo

end subroutine initialize_masks

end module MOM_grid_initialize
