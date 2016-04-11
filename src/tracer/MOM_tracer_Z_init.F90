module MOM_tracer_Z_init
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
!*  By Robert Hallberg, September 2009                                 *
!*                                                                     *
!*    This file contains a subroutine to initialize tracers into the   *
!*  MOM vertical grid from a depth- (or z*-) space file.               *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:                                           *
!*    j+1  > o > o >   At ^:                                           *
!*    j    x ^ x ^ x   At >:                                           *
!*    j    > o > o >   At o:  tr                                       *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_diag_to_Z, only : find_overlap, find_limited_slope
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
! use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : read_data

use netcdf

implicit none ; private

#include <MOM_memory.h>

public tracer_Z_init

contains

function tracer_Z_init(tr, h, filename, tr_name, G, missing_val, land_val)
  logical :: tracer_Z_init
  type(ocean_grid_type),                 intent(in)    :: G
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out)   :: tr
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h
  character(len=*),                      intent(in)    :: filename, tr_name
! type(param_file_type),                 intent(in)    :: param_file
  real,                        optional, intent(in)    :: missing_val
  real,                        optional, intent(in)    :: land_val
!   This function initializes a tracer by reading a Z-space file, returning
! .true. if this appears to have been successful, and false otherwise.
! Arguments: tr - The tracer to initialize.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      filename - The name of the file to read from.
!  (in)      tr_name - The name of the tracer in the file.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in,opt)  missing_val - The missing value for the tracer.
!  (in,opt)  land_val - The value to use to fill in land points.
  integer, save :: init_calls = 0
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_tracer_Z_init" ! This module's name.
  character(len=256) :: mesg    ! Message for error messages.

  real, allocatable, dimension(:,:,:) :: &
    tr_in   ! The z-space array of tracer concentrations that is read in.
  real, allocatable, dimension(:) :: &
    z_edges, &  ! The depths of the cell edges or cell centers (depending on
                ! the value of has_edges) in the input z* data.
    tr_1d, &    ! A copy of the input tracer concentrations in a column.
    wt, &   ! The fractional weight for each layer in the range between
            ! k_top and k_bot, nondim.
    z1, &   ! z1 and z2 are the depths of the top and bottom limits of the part
    z2      ! of a z-cell that contributes to a layer, relative to the cell
            ! center and normalized by the cell thickness, nondim.
            ! Note that -1/2 <= z1 <= z2 <= 1/2.
  real    :: e(SZK_(G)+1)  ! The z-star interface heights in m.
  real    :: landval    ! The tracer value to use in land points.
  real    :: sl_tr      ! The normalized slope of the tracer
                        ! within the cell, in tracer units.
  real    :: htot(SZI_(G)) ! The vertical sum of h, in m or kg m-2.
  real    :: dilate     ! The amount by which the thicknesses are dilated to
                        ! create a z-star coordinate, nondim or in m3 kg-1.
  real    :: missing  ! The missing value for the tracer.

  logical :: has_edges, use_missing, zero_surface
  character(len=80) :: loc_msg
  integer :: k_top, k_bot, k_bot_prev
  integer :: i, j, k, kz, is, ie, js, je, nz, nz_in
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  landval = 0.0 ; if (present(land_val)) landval = land_val

  zero_surface = .false. ! Make this false for errors to be fatal.

  use_missing = .false.
  if (present(missing_val)) then
    use_missing = .true. ; missing = missing_val
  endif

  ! Find out the number of input levels and read the depth of the edges,
  ! also modifying their sign convention to be monotonically decreasing.
  call read_Z_edges(filename, tr_name, z_edges, nz_in, has_edges, use_missing, missing)
  if (nz_in < 1) then
    tracer_Z_init = .false.
    return
  endif

  allocate(tr_in(G%isd:G%ied,G%jsd:G%jed,nz_in)) ; tr_in(:,:,:) = 0.0
  allocate(tr_1d(nz_in)) ; tr_1d(:) = 0.0
  call read_data(filename, tr_name, tr_in(:,:,:), domain=G%Domain%mpp_domain)

  ! Fill missing values from above?  Use a "close" test to avoid problems
  ! from type-conversion rounoff.
  if (present(missing_val)) then
    do j=js,je ; do i=is,ie
      if (G%mask2dT(i,j) == 0.0) then
        tr_in(i,j,1) = landval
      elseif (abs(tr_in(i,j,1) - missing_val) <= 1e-6*abs(missing_val)) then
        write(loc_msg,'(f7.2," N ",f7.2," E")') G%geoLatT(i,j), G%geoLonT(i,j)
        if (zero_surface) then
          call MOM_error(WARNING, "tracer_Z_init: Missing value of "// &
                trim(tr_name)//" found in an ocean point at "//trim(loc_msg)// &
                " in "//trim(filename) )
          tr_in(i,j,1) = 0.0
        else
          call MOM_error(FATAL, "tracer_Z_init: Missing value of "// &
                trim(tr_name)//" found in an ocean point at "//trim(loc_msg)// &
                " in "//trim(filename) )
        endif
      endif
    enddo ; enddo
    do k=2,nz_in ; do j=js,je ; do i=is,ie
      if (abs(tr_in(i,j,k) - missing_val) <= 1e-6*abs(missing_val)) &
        tr_in(i,j,k) = tr_in(i,j,k-1)
    enddo ; enddo ; enddo
  endif

  allocate(wt(nz_in+1)) ; allocate(z1(nz_in+1)) ; allocate(z2(nz_in+1))

  ! This is a placeholder, and will be replaced with our full vertical
  ! interpolation machinery when it is in place.
  if (has_edges) then
    do j=js,je
      do i=is,ie ; htot(i) = 0.0 ; enddo
      do k=1,nz ; do i=is,ie ; htot(i) = htot(i) + h(i,j,k) ; enddo ; enddo

      do i=is,ie ; if (G%mask2dT(i,j)*htot(i) > 0.0) then
        ! Determine the z* heights of the model interfaces.
        dilate = (G%bathyT(i,j) - 0.0) / htot(i)
        e(nz+1) = -G%bathyT(i,j)
        do k=nz,1,-1 ; e(K) = e(K+1) + dilate * h(i,j,k) ; enddo

        ! Create a single-column copy of tr_in.  ### CHANGE THIS LATER?
        do k=1,nz_in ; tr_1d(k) = tr_in(i,j,k) ; enddo
        k_bot = 1 ; k_bot_prev = -1
        do k=1,nz
          if (e(K+1) > z_edges(1)) then
            tr(i,j,k) = tr_1d(1)
          elseif (e(K) < z_edges(nz_in+1)) then
            tr(i,j,k) = tr_1d(nz_in)
          else
            call find_overlap(z_edges, e(K), e(K+1), nz_in, &
                              k_bot, k_top, k_bot, wt, z1, z2)
            kz = k_top
            if (kz /= k_bot_prev) then
              ! Calculate the intra-cell profile.
              sl_tr = 0.0 ! ; cur_tr = 0.0
              if ((kz < nz_in) .and. (kz > 1)) call &
                find_limited_slope(tr_1d, z_edges, sl_tr, kz)
            endif
            ! This is the piecewise linear form.
            tr(i,j,k) = wt(kz) * &
                (tr_1d(kz) + 0.5*sl_tr*(z2(kz) + z1(kz)))
            ! For the piecewise parabolic form add the following...
            !     + C1_3*cur_tr*(z2(kz)**2 + z2(kz)*z1(kz) + z1(kz)**2))
            do kz=k_top+1,k_bot-1
              tr(i,j,k) = tr(i,j,k) + wt(kz)*tr_1d(kz)
            enddo
            if (k_bot > k_top) then
              kz = k_bot
              ! Calculate the intra-cell profile.
              sl_tr = 0.0 ! ; cur_tr = 0.0
              if ((kz < nz_in) .and. (kz > 1)) call &
                find_limited_slope(tr_1d, z_edges, sl_tr, kz)
              ! This is the piecewise linear form.
              tr(i,j,k) = tr(i,j,k) + wt(kz) * &
                  (tr_1d(kz) + 0.5*sl_tr*(z2(kz) + z1(kz)))
              ! For the piecewise parabolic form add the following...
              !     + C1_3*cur_tr*(z2(kz)**2 + z2(kz)*z1(kz) + z1(kz)**2))
            endif
            k_bot_prev = k_bot

            !   Now handle the unlikely case where the layer partially extends
            ! past the valid range of the input data by extrapolating using
            ! the top or bottom value.
            if ((e(K) > z_edges(1)) .and. (z_edges(nz_in+1) > e(K+1))) then
              tr(i,j,k) = (((e(K) - z_edges(1)) * tr_1d(1) + &
                           (z_edges(1) - z_edges(nz_in)) * tr(i,j,k)) + &
                           (z_edges(nz_in+1) - e(K+1)) * tr_1d(nz_in)) / &
                          (e(K) - e(K+1))
            elseif (e(K) > z_edges(1)) then
              tr(i,j,k) = ((e(K) - z_edges(1)) * tr_1d(1) + &
                           (z_edges(1) - e(K+1)) * tr(i,j,k)) / &
                          (e(K) - e(K+1))
            elseif (z_edges(nz_in) > e(K+1)) then
              tr(i,j,k) = ((e(K) - z_edges(nz_in+1)) * tr(i,j,k) + &
                           (z_edges(nz_in+1) - e(K+1)) * tr_1d(nz_in)) / &
                          (e(K) - e(K+1))
            endif
          endif
        enddo ! k-loop
      else
        do k=1,nz ; tr(i,j,k) = landval ; enddo
      endif ; enddo ! i-loop
    enddo ! j-loop
  else
    ! Without edge values, integrate a linear interpolation between cell centers.
    do j=js,je
      do i=is,ie ; htot(i) = 0.0 ; enddo
      do k=1,nz ; do i=is,ie ; htot(i) = htot(i) + h(i,j,k) ; enddo ; enddo

      do i=is,ie ; if (G%mask2dT(i,j)*htot(i) > 0.0) then
        ! Determine the z* heights of the model interfaces.
        dilate = (G%bathyT(i,j) - 0.0) / htot(i)
        e(nz+1) = -G%bathyT(i,j)
        do k=nz,1,-1 ; e(K) = e(K+1) + dilate * h(i,j,k) ; enddo

        ! Create a single-column copy of tr_in.  ### CHANGE THIS LATER?
        do k=1,nz_in ; tr_1d(k) = tr_in(i,j,k) ; enddo
        k_bot = 1
        do k=1,nz
          if (e(K+1) > z_edges(1)) then
            tr(i,j,k) = tr_1d(1)
          elseif (z_edges(nz_in) > e(K)) then
            tr(i,j,k) = tr_1d(nz_in)
          else
            call find_overlap(z_edges, e(K), e(K+1), nz_in-1, &
                              k_bot, k_top, k_bot, wt, z1, z2)

            kz = k_top
            if (k_top < nz_in) then
              tr(i,j,k) = wt(kz)*0.5*((tr_1d(kz) + tr_1d(kz+1)) + &
                                      (tr_1d(kz+1) - tr_1d(kz))*(z2(kz)+z1(kz)))
            else
              tr(i,j,k) = wt(kz)*tr_1d(nz_in)
            endif
            do kz=k_top+1,k_bot-1
              tr(i,j,k) = tr(i,j,k) + wt(kz)*0.5*(tr_1d(kz) + tr_1d(kz+1))
            enddo
            if (k_bot > k_top) then
              kz = k_bot
              tr(i,j,k) = tr(i,j,k) + wt(kz)*0.5*((tr_1d(kz) + tr_1d(kz+1)) + &
                                        (tr_1d(kz+1) - tr_1d(kz))*(z2(kz)+z1(kz)))
            endif

            ! Now handle the case where the layer partially extends past
            ! the valid range of the input data.
            if ((e(K) > z_edges(1)) .and. (z_edges(nz_in) > e(K+1))) then
              tr(i,j,k) = (((e(K) - z_edges(1)) * tr_1d(1) + &
                           (z_edges(1) - z_edges(nz_in)) * tr(i,j,k)) + &
                           (z_edges(nz_in) - e(K+1)) * tr_1d(nz_in)) / &
                          (e(K) - e(K+1))
            elseif (e(K) > z_edges(1)) then
              tr(i,j,k) = ((e(K) - z_edges(1)) * tr_1d(1) + &
                           (z_edges(1) - e(K+1)) * tr(i,j,k)) / &
                          (e(K) - e(K+1))
            elseif (z_edges(nz_in) > e(K+1)) then
              tr(i,j,k) = ((e(K) - z_edges(nz_in)) * tr(i,j,k) + &
                           (z_edges(nz_in) - e(K+1)) * tr_1d(nz_in)) / &
                          (e(K) - e(K+1))
            endif
          endif
        enddo
      else
        do k=1,nz ; tr(i,j,k) = landval ; enddo
      endif ; enddo ! i-loop
    enddo  ! j-loop
  endif

  deallocate(tr_in) ; deallocate(tr_1d) ; deallocate(z_edges)
  deallocate(wt) ; deallocate(z1) ; deallocate(z2)

  tracer_Z_init = .true.

end function tracer_Z_init


subroutine read_Z_edges(filename, tr_name, z_edges, nz_out, has_edges, &
                        use_missing, missing)
  character(len=*),                intent(in)    :: filename, tr_name
  real, allocatable, dimension(:), intent(out)   :: z_edges
  integer,                         intent(out)   :: nz_out
  logical,                         intent(out)   :: has_edges
  logical,                         intent(inout) :: use_missing
  real,                            intent(inout) :: missing
!   This subroutine reads the vertical coordinate data for a field from a
! NetCDF file.  It also might read the missing value attribute for that
! same field.
! Arguments: filename - The file to be read from.
!  (in)      tr_name - The name of the tracer to be read.
!  (out)     z_edges - The depths of the vertical edges of the tracer array.
!  (out)     nz_out - The number of vertical layers in the tracer array.
!  (out)     has_edges - If true, the values in z_edges are the edges of the
!                        tracer cells, otherwise they are the cell centers.
!  (inout)   use_missing - If false on input, see whether the tracer has a
!                          missing value, and if so return true.
!  (inout)   missing - The missing value, if one has been found.

  character(len=32) :: mod
  character(len=120) :: dim_name, edge_name, tr_msg, dim_msg
  logical :: monotonic
  integer :: ncid, status, intid, tr_id, layid, k
  integer :: nz_edge, ndim, tr_dim_ids(NF90_MAX_VAR_DIMS)

  mod = "MOM_tracer_Z_init read_Z_edges: "
  tr_msg = trim(tr_name)//" in "//trim(filename)

  status = NF90_OPEN(filename, NF90_NOWRITE, ncid);
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,mod//" Difficulties opening "//trim(filename)//&
        " - "//trim(NF90_STRERROR(status)))
    nz_out = -1 ; return
  endif

  status = NF90_INQ_VARID(ncid, tr_name, tr_id)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,mod//" Difficulties finding variable "//&
        trim(tr_msg)//" - "//trim(NF90_STRERROR(status)))
    nz_out = -1 ; status = NF90_CLOSE(ncid) ; return
  endif
  status = NF90_INQUIRE_VARIABLE(ncid, tr_id, ndims=ndim, dimids=tr_dim_ids)
  if (status /= NF90_NOERR) then
    call MOM_ERROR(WARNING,mod//" cannot inquire about "//trim(tr_msg))
  elseif ((ndim < 3) .or. (ndim > 4)) then
    call MOM_ERROR(WARNING,mod//" "//trim(tr_msg)//&
         " has too many or too few dimensions.")
    nz_out = -1 ; status = NF90_CLOSE(ncid) ; return
  endif

  if (.not.use_missing) then
    ! Try to find the missing value from the dataset.
    status = NF90_GET_ATT(ncid, tr_id, "missing_value", missing)
    if (status /= NF90_NOERR) use_missing = .true.
  endif

  ! Get the axis name and length.
  status = NF90_INQUIRE_DIMENSION(ncid, tr_dim_ids(3), dim_name, len=nz_out)
  if (status /= NF90_NOERR) then
    call MOM_ERROR(WARNING,mod//" cannot inquire about dimension(3) of "//&
                    trim(tr_msg))
  endif

  dim_msg = trim(dim_name)//" in "//trim(filename)
  status = NF90_INQ_VARID(ncid, dim_name, layid)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,mod//" Difficulties finding variable "//&
        trim(dim_msg)//" - "//trim(NF90_STRERROR(status)))
    nz_out = -1 ; status = NF90_CLOSE(ncid) ; return
  endif
  ! Find out if the Z-axis has an edges attribute
  status = NF90_GET_ATT(ncid, layid, "edges", edge_name)
  if (status /= NF90_NOERR) then
    call MOM_mesg(mod//" "//trim(dim_msg)//&
         " has no readable edges attribute - "//trim(NF90_STRERROR(status)))
    has_edges = .false.
  else
    has_edges = .true.
    status = NF90_INQ_VARID(ncid, edge_name, intid)
    if (status /= NF90_NOERR) then
      call MOM_error(WARNING,mod//" Difficulties finding edge variable "//&
          trim(edge_name)//" in "//trim(filename)//" - "//trim(NF90_STRERROR(status)))
      has_edges = .false.
    endif
  endif

  nz_edge = nz_out ; if (has_edges) nz_edge = nz_out+1
  allocate(z_edges(nz_edge)) ; z_edges(:) = 0.0

  if (nz_out < 1) return

  ! Read the right variable.
  if (has_edges) then
    dim_msg = trim(edge_name)//" in "//trim(filename)
    status = NF90_GET_VAR(ncid, intid, z_edges)
    if (status /= NF90_NOERR) then
      call MOM_error(WARNING,mod//" Difficulties reading variable "//&
          trim(dim_msg)//" - "//trim(NF90_STRERROR(status)))
      nz_out = -1 ; status = NF90_CLOSE(ncid) ; return
    endif
  else
    status = NF90_GET_VAR(ncid, layid, z_edges)
    if (status /= NF90_NOERR) then
      call MOM_error(WARNING,mod//" Difficulties reading variable "//&
          trim(dim_msg)//" - "//trim(NF90_STRERROR(status)))
      nz_out = -1 ; status = NF90_CLOSE(ncid) ; return
    endif
  endif

  status = NF90_CLOSE(ncid)
  if (status /= NF90_NOERR) call MOM_error(WARNING, mod// &
    " Difficulties closing "//trim(filename)//" - "//trim(NF90_STRERROR(status)))

  ! z_edges should be montonically decreasing with our sign convention.
  ! Change the sign sign convention if it looks like z_edges is increasing.
  if (z_edges(1) < z_edges(2)) then
    do k=1,nz_edge ; z_edges(k) = -z_edges(k) ; enddo
  endif
  ! Check that z_edges is now monotonically decreasing.
  monotonic = .true.
  do k=2,nz_edge ; if (z_edges(k) >= z_edges(k-1)) monotonic = .false. ; enddo
  if (.not.monotonic) &
    call MOM_error(WARNING,mod//" "//trim(dim_msg)//" is not monotonic.")

end subroutine read_Z_edges


end module MOM_tracer_Z_init
