!> Used to initialize tracers from a depth- (or z*-) space file.
module MOM_tracer_Z_init

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
! use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : MOM_read_data
use MOM_EOS, only : EOS_type, calculate_density, calculate_density_derivs
use MOM_unit_scaling, only : unit_scale_type

use netcdf

implicit none ; private

#include <MOM_memory.h>

public tracer_Z_init, tracer_Z_init_array, find_interfaces, determine_temperature

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!>   This function initializes a tracer by reading a Z-space file, returning
!! .true. if this appears to have been successful, and false otherwise.
function tracer_Z_init(tr, h, filename, tr_name, G, US, missing_val, land_val)
  logical :: tracer_Z_init !< A return code indicating if the initialization has been successful
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type), intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                         intent(out)   :: tr   !< The tracer to initialize
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                         intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  character(len=*),      intent(in)    :: filename !< The name of the file to read from
  character(len=*),      intent(in)    :: tr_name !< The name of the tracer in the file
! type(param_file_type), intent(in)    :: param_file !< A structure to parse for run-time parameters
  real,        optional, intent(in)    :: missing_val !< The missing value for the tracer
  real,        optional, intent(in)    :: land_val !< A value to use to fill in land points

  !   This function initializes a tracer by reading a Z-space file, returning true if this
  ! appears to have been successful, and false otherwise.
!
  integer, save :: init_calls = 0
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_tracer_Z_init" ! This module's name.
  character(len=256) :: mesg    ! Message for error messages.

  real, allocatable, dimension(:,:,:) :: &
    tr_in   ! The z-space array of tracer concentrations that is read in.
  real, allocatable, dimension(:) :: &
    z_edges, &  ! The depths of the cell edges or cell centers (depending on
                ! the value of has_edges) in the input z* data [Z ~> m].
    tr_1d, &    ! A copy of the input tracer concentrations in a column.
    wt, &   ! The fractional weight for each layer in the range between
            ! k_top and k_bot, nondim.
    z1, &   ! z1 and z2 are the depths of the top and bottom limits of the part
    z2      ! of a z-cell that contributes to a layer, relative to the cell
            ! center and normalized by the cell thickness, nondim.
            ! Note that -1/2 <= z1 <= z2 <= 1/2.
  real    :: e(SZK_(G)+1)  ! The z-star interface heights [Z ~> m].
  real    :: landval    ! The tracer value to use in land points.
  real    :: sl_tr      ! The normalized slope of the tracer
                        ! within the cell, in tracer units.
  real    :: htot(SZI_(G)) ! The vertical sum of h [H ~> m or kg m-2].
  real    :: dilate     ! The amount by which the thicknesses are dilated to
                        ! create a z-star coordinate, nondim or in m3 kg-1.
  real    :: missing    ! The missing value for the tracer.

  logical :: has_edges, use_missing, zero_surface
  character(len=80) :: loc_msg
  integer :: k_top, k_bot, k_bot_prev, k_start
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
  call read_Z_edges(filename, tr_name, z_edges, nz_in, has_edges, use_missing, &
                    missing, scale=US%m_to_Z)
  if (nz_in < 1) then
    tracer_Z_init = .false.
    return
  endif

  allocate(tr_in(G%isd:G%ied,G%jsd:G%jed,nz_in)) ; tr_in(:,:,:) = 0.0
  allocate(tr_1d(nz_in)) ; tr_1d(:) = 0.0
  call MOM_read_data(filename, tr_name, tr_in(:,:,:), G%Domain)

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

        ! Create a single-column copy of tr_in.  Efficiency is not an issue here.
        do k=1,nz_in ; tr_1d(k) = tr_in(i,j,k) ; enddo
        k_bot = 1 ; k_bot_prev = -1
        do k=1,nz
          if (e(K+1) > z_edges(1)) then
            tr(i,j,k) = tr_1d(1)
          elseif (e(K) < z_edges(nz_in+1)) then
            tr(i,j,k) = tr_1d(nz_in)
          else
            k_start = k_bot ! The starting point for this search
            call find_overlap(z_edges, e(K), e(K+1), nz_in, &
                              k_start, k_top, k_bot, wt, z1, z2)
            kz = k_top
            if (kz /= k_bot_prev) then
              ! Calculate the intra-cell profile.
              sl_tr = 0.0 ! ; cur_tr = 0.0
              if ((kz < nz_in) .and. (kz > 1)) &
                sl_tr = find_limited_slope(tr_1d, z_edges, kz)
            endif
            ! This is the piecewise linear form.
            tr(i,j,k) = wt(kz) * (tr_1d(kz) + 0.5*sl_tr*(z2(kz) + z1(kz)))
            ! For the piecewise parabolic form add the following...
            !     + C1_3*cur_tr*(z2(kz)**2 + z2(kz)*z1(kz) + z1(kz)**2))
            do kz=k_top+1,k_bot-1
              tr(i,j,k) = tr(i,j,k) + wt(kz)*tr_1d(kz)
            enddo
            if (k_bot > k_top) then
              kz = k_bot
              ! Calculate the intra-cell profile.
              sl_tr = 0.0 ! ; cur_tr = 0.0
              if ((kz < nz_in) .and. (kz > 1)) &
                sl_tr = find_limited_slope(tr_1d, z_edges, kz)
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

        ! Create a single-column copy of tr_in.  Efficiency is not an issue here.
        do k=1,nz_in ; tr_1d(k) = tr_in(i,j,k) ; enddo
        k_bot = 1
        do k=1,nz
          if (e(K+1) > z_edges(1)) then
            tr(i,j,k) = tr_1d(1)
          elseif (z_edges(nz_in) > e(K)) then
            tr(i,j,k) = tr_1d(nz_in)
          else
            k_start = k_bot ! The starting point for this search
            call find_overlap(z_edges, e(K), e(K+1), nz_in-1, &
                              k_start, k_top, k_bot, wt, z1, z2)

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

!> Layer model routine for remapping tracers
!! from pseudo-z coordinates into layers defined
!! by target interface positions.
subroutine tracer_z_init_array(tr_in, z_edges, e, nkml, nkbl, land_fill, wet, nlay, nlevs, &
                       eps_z, tr)
  real, dimension(:,:,:),           intent(in) :: tr_in !< The z-space array of tracer concentrations that is read in.
  real, dimension(size(tr_in,3)+1), intent(in) :: z_edges !< The depths of the cell edges in the input z* data
                                                          !! [Z ~> m or m]
  integer,                          intent(in) :: nlay !< The number of vertical layers in the target grid
  real, dimension(size(tr_in,1),size(tr_in,2),nlay+1), &
                                    intent(in) :: e !< The depths of the target layer interfaces [Z ~> m or m]
  integer,                          intent(in) :: nkml !< The number of mixed layers
  integer,                          intent(in) :: nkbl !< The number of buffer layers
  real,                             intent(in) :: land_fill !< fill in data over land (1)
  real, dimension(size(tr_in,1),size(tr_in,2)), &
                                    intent(in) :: wet !< The wet mask for the source data (valid points)
  integer, dimension(size(tr_in,1),size(tr_in,2)), &
                                    intent(in) :: nlevs !< The number of input levels with valid data
  real,                             intent(in) :: eps_z !< A negligibly thin layer thickness [Z ~> m].
  real, dimension(size(tr_in,1),size(tr_in,2),nlay), intent(out) :: tr !< tracers in layer space

  ! Local variables
  real, dimension(size(tr_in,3)) :: tr_1d !< a copy of the input tracer concentrations in a column.
  real, dimension(nlay+1) :: e_1d ! A 1-d column of intreface heights, in the same units as e.
  real, dimension(nlay) :: tr_    ! A 1-d column of tracer concentrations
  integer :: n,i,j,k,l,nx,ny,nz,nt,kz
  integer :: k_top,k_bot,k_bot_prev,kk,kstart
  real    :: sl_tr    ! The tracer concentration slope times the layer thickness, in tracer units.
  real, dimension(size(tr_in,3)) :: wt !< The fractional weight for each layer in the range between z1 and z2
  real, dimension(size(tr_in,3)) :: z1, z2 ! z1 and z2 are the fractional depths of the top and bottom
                                  ! limits of the part of a z-cell that contributes to a layer, relative
                                  ! to the cell center and normalized by the cell thickness [nondim].
                                  ! Note that -1/2 <= z1 <= z2 <= 1/2.

  nx = size(tr_in,1); ny=size(tr_in,2); nz = size(tr_in,3)


  do j=1,ny
    i_loop: do i=1,nx
      if (nlevs(i,j) == 0 .or. wet(i,j) == 0.) then
        tr(i,j,:) = land_fill
        cycle i_loop
      endif

      do k=1,nz
        tr_1d(k) = tr_in(i,j,k)
      enddo

      do k=1,nlay+1
        e_1d(k) = e(i,j,k)
      enddo
      k_bot = 1 ; k_bot_prev = -1
      do k=1,nlay
        if (e_1d(k+1) > z_edges(1)) then
          tr(i,j,k) = tr_1d(1)
        elseif (e_1d(k) < z_edges(nlevs(i,j)+1)) then
          tr(i,j,k) = tr_1d(nlevs(i,j))

        else
          kstart=k_bot
          call find_overlap(z_edges, e_1d(k), e_1d(k+1), nlevs(i,j), &
                            kstart, k_top, k_bot, wt, z1, z2)
          kz = k_top
          sl_tr=0.0; ! cur_tr=0.0
          if (kz /= k_bot_prev) then
            ! Calculate the intra-cell profile.
            if ((kz < nlevs(i,j)) .and. (kz > 1)) then
              sl_tr = find_limited_slope(tr_1d, z_edges, kz)
            endif
          endif
          if (kz > nlevs(i,j)) kz = nlevs(i,j)
          ! This is the piecewise linear form.
          tr(i,j,k) = wt(kz) * (tr_1d(kz) + 0.5*sl_tr*(z2(kz) + z1(kz)))
          ! For the piecewise parabolic form add the following...
          !     + C1_3*wt(kz) * cur_tr*(z2(kz)**2 + z2(kz)*z1(kz) + z1(kz)**2))
          do kz=k_top+1,k_bot-1
            tr(i,j,k) = tr(i,j,k) + wt(kz)*tr_1d(kz)
          enddo

          if (k_bot > k_top) then
            kz = k_bot
            ! Calculate the intra-cell profile.
            sl_tr = 0.0 ! ; cur_tr = 0.0
            if ((kz < nlevs(i,j)) .and. (kz > 1)) then
              sl_tr = find_limited_slope(tr_1d, z_edges, kz)
            endif
            ! This is the piecewise linear form.
            tr(i,j,k) = tr(i,j,k) + wt(kz) * &
                 (tr_1d(kz) + 0.5*sl_tr*(z2(kz) + z1(kz)))
            ! For the piecewise parabolic form add the following...
            !     + C1_3*cur_tr*(z2(kz)**2 + z2(kz)*z1(kz) + z1(kz)**2))
          endif
          k_bot_prev = k_bot

        endif
      enddo ! k-loop

      do k=2,nlay  ! simply fill vanished layers with adjacent value
        if (e_1d(k)-e_1d(k+1) <= eps_z) tr(i,j,k)=tr(i,j,k-1)
      enddo

    enddo i_loop
  enddo

end subroutine tracer_z_init_array

!> This subroutine reads the vertical coordinate data for a field from a NetCDF file.
!! It also might read the missing value attribute for that same field.
subroutine read_Z_edges(filename, tr_name, z_edges, nz_out, has_edges, &
                        use_missing, missing, scale)
  character(len=*), intent(in)    :: filename !< The name of the file to read from.
  character(len=*), intent(in)    :: tr_name !< The name of the tracer in the file.
  real, dimension(:), allocatable, &
                    intent(out)   :: z_edges !< The depths of the vertical edges of the tracer array
  integer,          intent(out)   :: nz_out  !< The number of vertical layers in the tracer array
  logical,          intent(out)   :: has_edges !< If true the values in z_edges are the edges of the
                                             !! tracer cells, otherwise they are the cell centers
  logical,          intent(inout) :: use_missing !< If false on input, see whether the tracer has a
                                             !! missing value, and if so return true
  real,             intent(inout) :: missing !< The missing value, if one has been found
  real,             intent(in)    :: scale   !< A scaling factor for z_edges into new units.

  !   This subroutine reads the vertical coordinate data for a field from a
  ! NetCDF file.  It also might read the missing value attribute for that same field.
  character(len=32) :: mdl
  character(len=120) :: dim_name, edge_name, tr_msg, dim_msg
  logical :: monotonic
  integer :: ncid, status, intid, tr_id, layid, k
  integer :: nz_edge, ndim, tr_dim_ids(NF90_MAX_VAR_DIMS)

  mdl = "MOM_tracer_Z_init read_Z_edges: "
  tr_msg = trim(tr_name)//" in "//trim(filename)

  status = NF90_OPEN(filename, NF90_NOWRITE, ncid)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,mdl//" Difficulties opening "//trim(filename)//&
        " - "//trim(NF90_STRERROR(status)))
    nz_out = -1 ; return
  endif

  status = NF90_INQ_VARID(ncid, tr_name, tr_id)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,mdl//" Difficulties finding variable "//&
        trim(tr_msg)//" - "//trim(NF90_STRERROR(status)))
    nz_out = -1 ; status = NF90_CLOSE(ncid) ; return
  endif
  status = NF90_INQUIRE_VARIABLE(ncid, tr_id, ndims=ndim, dimids=tr_dim_ids)
  if (status /= NF90_NOERR) then
    call MOM_ERROR(WARNING,mdl//" cannot inquire about "//trim(tr_msg))
  elseif ((ndim < 3) .or. (ndim > 4)) then
    call MOM_ERROR(WARNING,mdl//" "//trim(tr_msg)//&
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
    call MOM_ERROR(WARNING,mdl//" cannot inquire about dimension(3) of "//&
                    trim(tr_msg))
  endif

  dim_msg = trim(dim_name)//" in "//trim(filename)
  status = NF90_INQ_VARID(ncid, dim_name, layid)
  if (status /= NF90_NOERR) then
    call MOM_error(WARNING,mdl//" Difficulties finding variable "//&
        trim(dim_msg)//" - "//trim(NF90_STRERROR(status)))
    nz_out = -1 ; status = NF90_CLOSE(ncid) ; return
  endif
  ! Find out if the Z-axis has an edges attribute
  status = NF90_GET_ATT(ncid, layid, "edges", edge_name)
  if (status /= NF90_NOERR) then
    call MOM_mesg(mdl//" "//trim(dim_msg)//&
         " has no readable edges attribute - "//trim(NF90_STRERROR(status)))
    has_edges = .false.
  else
    has_edges = .true.
    status = NF90_INQ_VARID(ncid, edge_name, intid)
    if (status /= NF90_NOERR) then
      call MOM_error(WARNING,mdl//" Difficulties finding edge variable "//&
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
      call MOM_error(WARNING,mdl//" Difficulties reading variable "//&
          trim(dim_msg)//" - "//trim(NF90_STRERROR(status)))
      nz_out = -1 ; status = NF90_CLOSE(ncid) ; return
    endif
  else
    status = NF90_GET_VAR(ncid, layid, z_edges)
    if (status /= NF90_NOERR) then
      call MOM_error(WARNING,mdl//" Difficulties reading variable "//&
          trim(dim_msg)//" - "//trim(NF90_STRERROR(status)))
      nz_out = -1 ; status = NF90_CLOSE(ncid) ; return
    endif
  endif

  status = NF90_CLOSE(ncid)
  if (status /= NF90_NOERR) call MOM_error(WARNING, mdl// &
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
    call MOM_error(WARNING,mdl//" "//trim(dim_msg)//" is not monotonic.")

  if (scale /= 1.0) then ; do k=1,nz_edge ; z_edges(k) = scale*z_edges(k) ; enddo ; endif

end subroutine read_Z_edges

!### `find_overlap` and `find_limited_slope` were previously part of
!    MOM_diag_to_Z.F90, and are nearly identical to `find_overlap` in
!    `midas_vertmap.F90` with some slight differences.  We keep it here for
!    reproducibility, but the two should be merged at some point

!> Determines the layers bounded by interfaces e that overlap
!! with the depth range between Z_top and Z_bot, and the fractional weights
!! of each layer. It also calculates the normalized relative depths of the range
!! of each layer that overlaps that depth range.
subroutine find_overlap(e, Z_top, Z_bot, k_max, k_start, k_top, k_bot, wt, z1, z2)
  real, dimension(:), intent(in)  :: e      !< Column interface heights, [Z ~> m] or other units.
  real,               intent(in)  :: Z_top  !< Top of range being mapped to, in the units of e [Z ~> m].
  real,               intent(in)  :: Z_bot  !< Bottom of range being mapped to, in the units of e [Z ~> m].
  integer,            intent(in)  :: k_max  !< Number of valid layers.
  integer,            intent(in)  :: k_start !< Layer at which to start searching.
  integer,            intent(out) :: k_top  !< Indices of top layers that overlap with the depth range.
  integer,            intent(out) :: k_bot  !< Indices of bottom layers that overlap with the depth range.
  real, dimension(:), intent(out) :: wt     !< Relative weights of each layer from k_top to k_bot [nondim].
  real, dimension(:), intent(out) :: z1     !< Depth of the top limits of the part of
       !! a layer that contributes to a depth level, relative to the cell center and normalized
       !! by the cell thickness [nondim].  Note that -1/2 <= z1 < z2 <= 1/2.
  real, dimension(:), intent(out)   :: z2     !< Depths of the bottom limit of the part of
       !! a layer that contributes to a depth level, relative to the cell center and normalized
       !! by the cell thickness [nondim].  Note that -1/2 <= z1 < z2 <= 1/2.
  ! Local variables
  real    :: Ih, e_c, tot_wt, I_totwt
  integer :: k

  wt(:) = 0.0 ; z1(:) = 0.0 ; z2(:) = 0.0 ; k_bot = k_max

  do k=k_start,k_max ; if (e(K+1) < Z_top) exit ; enddo
  k_top = k
  if (k_top > k_max) return

  ! Determine the fractional weights of each layer.
  ! Note that by convention, e and Z_int decrease with increasing k.
  if (e(K+1) <= Z_bot) then
    wt(k) = 1.0 ; k_bot = k
    Ih = 0.0 ; if (e(K) /= e(K+1)) Ih = 1.0 / (e(K)-e(K+1))
    e_c = 0.5*(e(K)+e(K+1))
    z1(k) = (e_c - MIN(e(K), Z_top)) * Ih
    z2(k) = (e_c - Z_bot) * Ih
  else
    wt(k) = MIN(e(K),Z_top) - e(K+1) ; tot_wt = wt(k) ! These are always > 0.
    if (e(K) /= e(K+1)) then
      z1(k) = (0.5*(e(K)+e(K+1)) - MIN(e(K), Z_top)) / (e(K)-e(K+1))
    else ; z1(k) = -0.5 ; endif
    z2(k) = 0.5
    k_bot = k_max
    do k=k_top+1,k_max
      if (e(K+1) <= Z_bot) then
        k_bot = k
        wt(k) = e(K) - Z_bot ; z1(k) = -0.5
        if (e(K) /= e(K+1)) then
          z2(k) = (0.5*(e(K)+e(K+1)) - Z_bot) / (e(K)-e(K+1))
        else ; z2(k) = 0.5 ; endif
      else
        wt(k) = e(K) - e(K+1) ; z1(k) = -0.5 ; z2(k) = 0.5
      endif
      tot_wt = tot_wt + wt(k) ! wt(k) is always > 0.
      if (k>=k_bot) exit
    enddo

    I_totwt = 0.0 ; if (tot_wt > 0.0) I_totwt = 1.0 / tot_wt
    do k=k_top,k_bot ; wt(k) = I_totwt*wt(k) ; enddo
  endif

end subroutine find_overlap

!> This subroutine determines a limited slope for val to be advected with
!! a piecewise limited scheme.
function find_limited_slope(val, e, k) result(slope)
  real, dimension(:), intent(in) :: val !< An column the values that are being interpolated.
  real, dimension(:), intent(in) :: e   !< A column's interface heights [Z ~> m] or other units.
  integer,            intent(in) :: k   !< The layer whose slope is being determined.
  real :: slope !< The normalized slope in the intracell distribution of val.
  ! Local variables
  real :: amn, cmn
  real :: d1, d2

  if ((val(k)-val(k-1)) * (val(k)-val(k+1)) >= 0.0) then
    slope = 0.0 ! ; curvature = 0.0
  else
    d1 = 0.5*(e(K-1)-e(K+1)) ; d2 = 0.5*(e(K)-e(K+2))
    if (d1*d2 > 0.0) then
      slope = ((d1**2)*(val(k+1) - val(k)) + (d2**2)*(val(k) - val(k-1))) * &
              (e(K) - e(K+1)) / (d1*d2*(d1+d2))
      ! slope = 0.5*(val(k+1) - val(k-1))
      ! This is S.J. Lin's form of the PLM limiter.
      amn = min(abs(slope), 2.0*(max(val(k-1), val(k), val(k+1)) - val(k)))
      cmn = 2.0*(val(k) - min(val(k-1), val(k), val(k+1)))
      slope = sign(1.0, slope) * min(amn, cmn)

      ! min(abs(slope), 2.0*(max(val(k-1),val(k),val(k+1)) - val(k)), &
      !                 2.0*(val(k) - min(val(k-1),val(k),val(k+1))))
      ! curvature = 0.0
    else
      slope = 0.0 ! ; curvature = 0.0
    endif
  endif

end function find_limited_slope

!> Find interface positions corresponding to density profile
function find_interfaces(rho, zin, Rb, depth, nlevs, nkml, nkbl, hml, debug, eps_z, eps_rho) result(zi)
  real, dimension(:,:,:), &
                      intent(in) :: rho   !< potential density in z-space [kg m-3 or R ~> kg m-3]
  real, dimension(size(rho,3)), &
                      intent(in) :: zin   !< Input data levels [m or Z ~> m].
  real, dimension(:), intent(in) :: Rb    !< target interface densities [kg m-3 or R ~> kg m-3]
  real, dimension(size(rho,1),size(rho,2)), &
                      intent(in) :: depth !< ocean depth [Z ~> m].
  integer, dimension(size(rho,1),size(rho,2)), &
            optional, intent(in) :: nlevs !< number of valid points in each column
  logical,  optional, intent(in) :: debug !< optional debug flag
  integer,  optional, intent(in) :: nkml  !< number of mixed layer pieces
  integer,  optional, intent(in) :: nkbl  !< number of buffer layer pieces
  real,     optional, intent(in) :: hml   !< mixed layer depth [Z ~> m].
  real,     optional, intent(in) :: eps_z !< A negligibly small layer thickness [m or Z ~> m].
  real,     optional, intent(in) :: eps_rho !< A negligibly small density difference [kg m-3 or R ~> kg m-3].
  real, dimension(size(rho,1),size(rho,2),size(Rb,1)) :: zi !< The returned interface, in the same units az zin.

  ! Local variables
  real, dimension(size(rho,1),size(rho,3)) :: rho_ ! A slice of densities [R ~> kg m-3]
  real, dimension(size(rho,1)) :: depth_
  logical :: unstable
  integer :: dir
  integer, dimension(size(rho,1),size(Rb,1)) :: ki_
  real, dimension(size(rho,1),size(Rb,1)) :: zi_
  integer, dimension(size(rho,1),size(rho,2)) :: nlevs_data
  integer, dimension(size(rho,1)) :: lo, hi
  real :: slope,rsm,drhodz,hml_
  integer :: n,i,j,k,l,nx,ny,nz,nt
  integer :: nlay,kk,nkml_,nkbl_
  logical :: debug_ = .false.
  real    :: epsln_Z    ! A negligibly thin layer thickness [m or Z ~> m].
  real    :: epsln_rho  ! A negligibly small density change [kg m-3 or R ~> kg m-3].
  real, parameter :: zoff=0.999

  nlay=size(Rb)-1

  zi(:,:,:) = 0.0

  if (PRESENT(debug)) debug_=debug

  nx = size(rho,1); ny=size(rho,2); nz = size(rho,3)
  nlevs_data(:,:) = size(rho,3)

  nkml_ = 0 ;  if (PRESENT(nkml)) nkml_ = max(0, nkml)
  nkbl_ = 0 ;  if (PRESENT(nkbl)) nkbl_ = max(0, nkbl)
  hml_ = 0.0 ; if (PRESENT(hml)) hml_ = hml
  epsln_Z = 1.0e-10 ; if (PRESENT(eps_z)) epsln_Z = eps_z
  epsln_rho = 1.0e-10 ; if (PRESENT(eps_rho)) epsln_rho = eps_rho

  if (PRESENT(nlevs)) then
    nlevs_data(:,:) = nlevs(:,:)
  endif

  do j=1,ny
    rho_(:,:) = rho(:,j,:)
    i_loop: do i=1,nx
      if (debug_) then
        print *,'looking for interfaces, i,j,nlevs= ',i,j,nlevs_data(i,j)
        print *,'initial density profile= ', rho_(i,:)
      endif
      unstable=.true.
      dir=1
      do while (unstable)
        unstable=.false.
        if (dir == 1) then
          do k=2,nlevs_data(i,j)-1
            if (rho_(i,k) - rho_(i,k-1) < 0.0 ) then
              if (k == 2) then
                rho_(i,k-1) = rho_(i,k)-epsln_rho
              else
                drhodz = (rho_(i,k+1)-rho_(i,k-1)) / (zin(k+1)-zin(k-1))
                if (drhodz < 0.0) unstable=.true.
                rho_(i,k) = rho_(i,k-1) + drhodz*zoff*(zin(k)-zin(k-1))
              endif
            endif
          enddo
          dir = -1*dir
        else
          do k=nlevs_data(i,j)-1,2,-1
            if (rho_(i,k+1) - rho_(i,k) < 0.0) then
              if (k == nlevs_data(i,j)-1) then
                rho_(i,k+1) = rho_(i,k-1)+epsln_rho
              else
                drhodz = (rho_(i,k+1)-rho_(i,k-1))/(zin(k+1)-zin(k-1))
                if (drhodz  < 0.0) unstable=.true.
                rho_(i,k) = rho_(i,k+1)-drhodz*(zin(k+1)-zin(k))
              endif
            endif
          enddo
          dir = -1*dir
        endif
      enddo
      if (debug_) then
        print *,'final density profile= ', rho_(i,:)
      endif
    enddo i_loop

    ki_(:,:) = 0
    zi_(:,:) = 0.0
    depth_(:) = -1.0*depth(:,j)
    lo(:) = 1
    hi(:) = nlevs_data(:,j)
    ki_ = bisect_fast(rho_, Rb, lo, hi)
    ki_(:,:) = max(1, ki_(:,:)-1)
    do i=1,nx
      do l=2,nlay
        slope = (zin(ki_(i,l)+1) - zin(ki_(i,l))) / max(rho_(i,ki_(i,l)+1) - rho_(i,ki_(i,l)),epsln_rho)
        zi_(i,l) = -1.0*(zin(ki_(i,l)) + slope*(Rb(l)-rho_(i,ki_(i,l))))
        zi_(i,l) = max(zi_(i,l), depth_(i))
        zi_(i,l) = min(zi_(i,l), -1.0*hml_)
      enddo
      zi_(i,nlay+1) = depth_(i)
      do l=2,nkml_+1
        zi_(i,l) = max(hml_*((1.0-real(l))/real(nkml_)), depth_(i))
      enddo
      do l=nlay,nkml_+2,-1
        if (zi_(i,l) < zi_(i,l+1) + epsln_Z) zi_(i,l) = zi_(i,l+1) + epsln_Z
        if (zi_(i,l) > -1.0*hml_)  zi_(i,l) = max(-1.0*hml_, depth_(i))
      enddo
    enddo
    zi(:,j,:) = zi_(:,:)
  enddo

end function find_interfaces

!> This subroutine determines the potential temperature and salinity that
!! is consistent with the target density using provided initial guess
subroutine determine_temperature(temp, salt, R, p_ref, niter, land_fill, h, k_start, eos)
  real, dimension(:,:,:),        intent(inout) :: temp !< potential temperature [degC]
  real, dimension(:,:,:),        intent(inout) :: salt !< salinity [PSU]
  real, dimension(size(temp,3)), intent(in)    :: R !< desired potential density [kg m-3].
  real,                          intent(in)    :: p_ref !< reference pressure [Pa].
  integer,                       intent(in)    :: niter !< maximum number of iterations
  integer,                       intent(in)    :: k_start !< starting index (i.e. below the buffer layer)
  real,                          intent(in)    :: land_fill !< land fill value
  real, dimension(:,:,:),        intent(in)    :: h   !< layer thickness, used only to avoid working on massless layers
  type(eos_type),                pointer       :: eos !< seawater equation of state control structure

  real, parameter :: T_max = 31.0, T_min = -2.0
  ! Local variables (All of which need documentation!)
  real(kind=8), dimension(size(temp,1),size(temp,3)) :: T, S, dT, dS, rho, hin
  real(kind=8), dimension(size(temp,1),size(temp,3)) :: drho_dT, drho_dS
  real(kind=8), dimension(size(temp,1)) :: press
  integer :: nx, ny, nz, nt, i, j, k, n, itt
  real    :: dT_dS_gauge  ! The relative penalizing of temperature to salinity changes when
                          ! minimizing property changes while correcting density [degC ppt-1].
  real    :: I_denom      ! The inverse of the magnitude squared of the density gradient in
                          ! T-S space streched with dT_dS_gauge [m6 kg-2 ppt-1]
  logical :: adjust_salt, old_fit
  real, parameter :: S_min = 0.5, S_max=65.0
  real, parameter :: tol_T=1.e-4, tol_S=1.e-4, tol_rho=1.e-4
  real, parameter :: max_t_adj=1.0, max_s_adj = 0.5

  old_fit = .true.   ! reproduces siena behavior

  ! ### The whole determine_temperature subroutine needs to be reexamined, both the algorithms
  !     and the extensive use of hard-coded dimensional parameters.

  ! We will switch to the newer method which simultaneously adjusts
  ! temp and salt based on the ratio of the thermal and haline coefficients, once it is tested.

  nx=size(temp,1) ; ny=size(temp,2) ; nz=size(temp,3)

  press(:) = p_ref

  do j=1,ny
    dS(:,:) = 0. ! Needs to be zero everywhere since there is a maxval(abs(dS)) later...
    T=temp(:,j,:)
    S=salt(:,j,:)
    hin=h(:,j,:)
    dT=0.0
    adjust_salt = .true.
    iter_loop: do itt = 1,niter
      do k=1, nz
        call calculate_density(T(:,k), S(:,k), press, rho(:,k), 1, nx, eos)
        call calculate_density_derivs(T(:,k), S(:,k), press, drho_dT(:,k), drho_dS(:,k), 1, nx, eos)
      enddo
      do k=k_start,nz ; do i=1,nx

!       if (abs(rho(i,k)-R(k))>tol .and. hin(i,k)>epsln .and. abs(T(i,k)-land_fill) < epsln) then
        if (abs(rho(i,k)-R(k))>tol_rho) then
          if (old_fit) then
            dT(i,k) = max(min((R(k)-rho(i,k)) / drho_dT(i,k), max_t_adj), -max_t_adj)
            T(i,k) = max(min(T(i,k)+dT(i,k), T_max), T_min)
          else
            dT_dS_gauge = 10.0  ! 10 degC is weighted equivalently to 1 ppt.
            I_denom = 1.0 / (drho_dS(i,k)**2 + dT_dS_gauge**2*drho_dT(i,k)**2)
            dS(i,k) = (R(k)-rho(i,k)) * drho_dS(i,k) * I_denom
            dT(i,k) = (R(k)-rho(i,k)) * dT_dS_gauge**2*drho_dT(i,k) * I_denom

            T(i,k) = max(min(T(i,k)+dT(i,k), T_max), T_min)
            S(i,k) = max(min(S(i,k)+dS(i,k), S_max), S_min)
          endif
        endif
      enddo ; enddo
      if (maxval(abs(dT)) < tol_T) then
         adjust_salt = .false.
         exit iter_loop
      endif
    enddo iter_loop

    if (adjust_salt .and. old_fit) then ; do itt = 1,niter
      do k=1, nz
        call calculate_density(T(:,k),S(:,k),press,rho(:,k),1,nx,eos)
        call calculate_density_derivs(T(:,k),S(:,k),press,drho_dT(:,k),drho_dS(:,k),1,nx,eos)
      enddo
      do k=k_start,nz ; do i=1,nx
!       if (abs(rho(i,k)-R(k))>tol .and. hin(i,k)>epsln .and. abs(T(i,k)-land_fill) < epsln ) then
        if (abs(rho(i,k)-R(k)) > tol_rho) then
          dS(i,k) = max(min((R(k)-rho(i,k)) / drho_dS(i,k), max_s_adj), -max_s_adj)
          S(i,k) = max(min(S(i,k)+dS(i,k), S_max), S_min)
        endif
      enddo ; enddo
      if (maxval(abs(dS)) < tol_S) exit
    enddo ; endif

    temp(:,j,:)=T(:,:)
    salt(:,j,:)=S(:,:)
  enddo

end subroutine determine_temperature

!> Return the index where to insert item x in list a, assuming a is sorted.
!! The return values [i] is such that all e in a[:i-1] have e <= x, and all e in
!! a[i:] have e > x. So if x already appears in the list, will
!! insert just after the rightmost x already there.
!! Optional args lo (default 1) and hi (default len(a)) bound the
!! slice of a to be searched.
function bisect_fast(a, x, lo, hi) result(bi_r)
  real, dimension(:,:), intent(in) :: a !< Sorted list
  real, dimension(:), intent(in) :: x !< Item to be inserted
  integer, dimension(size(a,1)), optional, intent(in) :: lo !< Lower bracket of optional range to search
  integer, dimension(size(a,1)), optional, intent(in) :: hi !< Upper bracket of optional range to search
  integer, dimension(size(a,1),size(x,1))  :: bi_r

  integer :: mid,num_x,num_a,i
  integer, dimension(size(a,1))  :: lo_,hi_,lo0,hi0
  integer :: nprofs,j

  lo_=1;hi_=size(a,2);num_x=size(x,1);bi_r=-1;nprofs=size(a,1)

  if (PRESENT(lo)) then
     where (lo>0) lo_=lo
  endif
  if (PRESENT(hi)) then
     where (hi>0) hi_=hi
  endif

  lo0=lo_;hi0=hi_

  do j=1,nprofs
    do i=1,num_x
      lo_=lo0;hi_=hi0
      do while (lo_(j) < hi_(j))
        mid = (lo_(j)+hi_(j))/2
        if (x(i) < a(j,mid)) then
           hi_(j) = mid
        else
           lo_(j) = mid+1
        endif
      enddo
      bi_r(j,i)=lo_(j)
    enddo
  enddo


  return

end function bisect_fast

end module MOM_tracer_Z_init
