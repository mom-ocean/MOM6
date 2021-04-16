!> Horizontal interpolation
module MOM_horizontal_regridding

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging,     only : hchksum
use MOM_coms,          only : max_across_PEs, min_across_PEs, sum_across_PEs, broadcast
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_LOOP
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_interpolate,   only : time_interp_external, horiz_interp_init
use MOM_interpolate,   only : build_horiz_interp_weights, run_horiz_interp, horiz_interp_type
use MOM_interp_infra,  only : axistype, get_external_field_info, get_axis_data
use MOM_time_manager,  only : time_type

use netcdf, only : NF90_OPEN, NF90_NOWRITE, NF90_GET_ATT, NF90_GET_VAR
use netcdf, only : NF90_INQ_VARID, NF90_INQUIRE_VARIABLE, NF90_INQUIRE_DIMENSION

implicit none ; private

#include <MOM_memory.h>

public :: horiz_interp_and_extrap_tracer, myStats

!> Extrapolate and interpolate data
interface horiz_interp_and_extrap_tracer
  module procedure horiz_interp_and_extrap_tracer_record
  module procedure horiz_interp_and_extrap_tracer_fms_id
end interface

contains

!> Write to the terminal some basic statistics about the k-th level of an array
subroutine myStats(array, missing, is, ie, js, je, k, mesg)
  real, dimension(:,:), intent(in) :: array !< input array (ND)
  real, intent(in) :: missing !< missing value (ND)
  integer :: is   !< Start index in i
  integer :: ie   !< End index in i
  integer :: js   !< Start index in j
  integer :: je   !< End index in j
  integer :: k !< Level to calculate statistics for
  character(len=*) :: mesg !< Label to use in message
  ! Local variables
  real :: minA, maxA
  integer :: i,j
  logical :: found
  character(len=120) :: lMesg
  minA = 9.E24 ; maxA = -9.E24 ; found = .false.

  do j=js,je ; do i=is,ie
    if (array(i,j) /= array(i,j)) stop 'Nan!'
    if (abs(array(i,j)-missing) > 1.e-6*abs(missing)) then
      if (found) then
        minA = min(minA, array(i,j))
        maxA = max(maxA, array(i,j))
      else
        found = .true.
        minA = array(i,j)
        maxA = array(i,j)
      endif
    endif
  enddo ; enddo
  call min_across_PEs(minA)
  call max_across_PEs(maxA)
  if (is_root_pe()) then
    write(lMesg(1:120),'(2(a,es12.4),a,i3,x,a)') &
         'init_from_Z: min=',minA,' max=',maxA,' Level=',k,trim(mesg)
    call MOM_mesg(lMesg,2)
  endif

end subroutine myStats

!> Use ICE-9 algorithm to populate points (fill=1) with valid data (good=1).  If no information
!! is available, use a previous guess (prev). Optionally (smooth) blend the filled points to
!! achieve a more desirable result.
subroutine fill_miss_2d(aout, good, fill, prev, G, smooth, num_pass, relc, crit, debug, answers_2018)
  use MOM_coms, only : sum_across_PEs

  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(inout) :: aout !< The array with missing values to fill
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(in)    :: good !< Valid data mask for incoming array
                                               !! (1==good data; 0==missing data).
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(in)    :: fill !< Same shape array of points which need
                                               !! filling (1==fill;0==dont fill)
  real, dimension(SZI_(G),SZJ_(G)), &
               optional, intent(in)    :: prev !< First guess where isolated holes exist.
  logical,     optional, intent(in)    :: smooth !< If present and true, apply a number of
                                                 !! Laplacian iterations to the interpolated data
  integer,     optional, intent(in)    :: num_pass !< The maximum number of iterations
  real,        optional, intent(in)    :: relc !< A relaxation coefficient for Laplacian (ND)
  real,        optional, intent(in)    :: crit !< A minimal value for deltas between iterations.
  logical,     optional, intent(in)    :: debug !< If true, write verbose debugging messages.
  logical,     optional, intent(in)    :: answers_2018 !< If true, use expressions that give the same
                                                !! answers as the code did in late 2018.  Otherwise
                                                !! add parentheses for rotational symmetry.

  real, dimension(SZI_(G),SZJ_(G)) :: a_filled ! The aout with missing values filled in
  real, dimension(SZI_(G),SZJ_(G)) :: a_chg    ! The change in aout due to an iteration of smoothing
  real, dimension(SZI_(G),SZJ_(G)) :: fill_pts ! 1 for points that still need to be filled
  real, dimension(SZI_(G),SZJ_(G)) :: good_    ! The values that are valid for the current iteration
  real, dimension(SZI_(G),SZJ_(G)) :: good_new ! The values of good_ to use for the next iteration

  real    :: east, west, north, south ! Valid neighboring values or 0 for invalid values
  real    :: ge, gw, gn, gs  ! Flags indicating which neighbors have valid values
  real    :: ngood     ! The number of valid values in neighboring points
  logical :: do_smooth ! Indicates whether to do smoothing of the array
  real    :: nfill     ! The remaining number of points to fill
  real    :: nfill_prev ! The previous value of nfill
  character(len=256) :: mesg  ! The text of an error message
  integer :: i, j, k
  integer, parameter :: num_pass_default = 10000
  real, parameter :: relc_default = 0.25, crit_default = 1.e-3

  integer :: npass
  integer :: is, ie, js, je
  real    :: relax_coeff, acrit, ares
  logical :: debug_it, ans_2018

  debug_it=.false.
  if (PRESENT(debug)) debug_it=debug

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  npass = num_pass_default
  if (PRESENT(num_pass)) npass = num_pass

  relax_coeff = relc_default
  if (PRESENT(relc)) relax_coeff = relc

  acrit = crit_default
  if (PRESENT(crit)) acrit = crit

  do_smooth=.false.
  if (PRESENT(smooth)) do_smooth=smooth

  ans_2018 = .true. ; if (PRESENT(answers_2018)) ans_2018 = answers_2018

  fill_pts(:,:) = fill(:,:)

  nfill = sum(fill(is:ie,js:je))
  call sum_across_PEs(nfill)

  nfill_prev = nfill
  good_(:,:) = good(:,:)
  a_chg(:,:) = 0.0

  do while (nfill > 0.0)

    call pass_var(good_,G%Domain)
    call pass_var(aout,G%Domain)

    a_filled(:,:) = aout(:,:)
    good_new(:,:) = good_(:,:)

    do j=js,je ; do i=is,ie

      if (good_(i,j) == 1.0 .or. fill(i,j) == 0.) cycle

      ge=good_(i+1,j) ; gw=good_(i-1,j)
      gn=good_(i,j+1) ; gs=good_(i,j-1)
      east=0.0 ; west=0.0 ; north=0.0 ; south=0.0
      if (ge == 1.0) east = aout(i+1,j)*ge
      if (gw == 1.0) west = aout(i-1,j)*gw
      if (gn == 1.0) north = aout(i,j+1)*gn
      if (gs == 1.0) south = aout(i,j-1)*gs

      if (ans_2018) then
        ngood = ge+gw+gn+gs
      else
        ngood = (ge+gw) + (gn+gs)
      endif
      if (ngood > 0.) then
        if (ans_2018) then
          a_filled(i,j) = (east+west+north+south)/ngood
        else
          a_filled(i,j) = ((east+west) + (north+south))/ngood
        endif
        fill_pts(i,j) = 0.0
        good_new(i,j) = 1.0
      endif
    enddo ; enddo

    aout(is:ie,js:je) = a_filled(is:ie,js:je)
    good_(is:ie,js:je) = good_new(is:ie,js:je)
    nfill_prev = nfill
    nfill = sum(fill_pts(is:ie,js:je))
    call sum_across_PEs(nfill)

    if (nfill == nfill_prev .and. PRESENT(prev)) then
      do j=js,je ; do i=is,ie ; if (fill_pts(i,j) == 1.0) then
        aout(i,j) = prev(i,j)
        fill_pts(i,j) = 0.0
      endif ; enddo ; enddo
    elseif (nfill == nfill_prev) then
      call MOM_error(WARNING, &
           'Unable to fill missing points using either data at the same vertical level from a connected basin'//&
           'or using a point from a previous vertical level.  Make sure that the original data has some valid'//&
           'data in all basins.', .true.)
      write(mesg,*) 'nfill=',nfill
      call MOM_error(WARNING, mesg, .true.)
    endif

    ! Determine the number of remaining points to fill globally.
    nfill = sum(fill_pts(is:ie,js:je))
    call sum_across_PEs(nfill)

  enddo ! while block for remaining points to fill.

  ! Do Laplacian smoothing for the points that have been filled in.
  if (do_smooth) then ; do k=1,npass
    call pass_var(aout,G%Domain)
    do j=js,je ; do i=is,ie
      if (fill(i,j) == 1) then
        east = max(good(i+1,j),fill(i+1,j)) ; west = max(good(i-1,j),fill(i-1,j))
        north = max(good(i,j+1),fill(i,j+1)) ; south = max(good(i,j-1),fill(i,j-1))
        if (ans_2018) then
          a_chg(i,j) = relax_coeff*(south*aout(i,j-1)+north*aout(i,j+1) + &
                                    west*aout(i-1,j)+east*aout(i+1,j) - &
                                   (south+north+west+east)*aout(i,j))
        else
          a_chg(i,j) = relax_coeff*( ((south*aout(i,j-1) + north*aout(i,j+1)) + &
                                  (west*aout(i-1,j)+east*aout(i+1,j))) - &
                                 ((south+north)+(west+east))*aout(i,j) )
        endif
      else
        a_chg(i,j) = 0.
      endif
    enddo ; enddo
    ares = 0.0
    do j=js,je ; do i=is,ie
      aout(i,j) = a_chg(i,j) + aout(i,j)
      ares = max(ares, abs(a_chg(i,j)))
    enddo ; enddo
    call max_across_PEs(ares)
    if (ares <= acrit) exit
  enddo ; endif

  do j=js,je ; do i=is,ie
    if (good_(i,j) == 0.0 .and. fill_pts(i,j) == 1.0) then
      write(mesg,*) 'In fill_miss, fill, good,i,j= ',fill_pts(i,j),good_(i,j),i,j
      call MOM_error(WARNING, mesg, .true.)
      call MOM_error(FATAL,"MOM_initialize: "// &
           "fill is true and good is false after fill_miss, how did this happen? ")
    endif
 enddo ; enddo

end subroutine fill_miss_2d

!> Extrapolate and interpolate from a file record
subroutine horiz_interp_and_extrap_tracer_record(filename, varnam,  conversion, recnum, G, tr_z, &
                                                 mask_z, z_in, z_edges_in, missing_value, reentrant_x, &
                                                 tripolar_n, homogenize, m_to_Z, answers_2018, ongrid)

  character(len=*),      intent(in)    :: filename   !< Path to file containing tracer to be
                                                     !! interpolated.
  character(len=*),      intent(in)    :: varnam     !< Name of tracer in filee.
  real,                  intent(in)    :: conversion !< Conversion factor for tracer.
  integer,               intent(in)    :: recnum     !< Record number of tracer to be read.
  type(ocean_grid_type), intent(inout) :: G          !< Grid object
  real, allocatable, dimension(:,:,:)  :: tr_z       !< pointer to allocatable tracer array on local
                                                     !! model grid and input-file vertical levels.
  real, allocatable, dimension(:,:,:)  :: mask_z     !< pointer to allocatable tracer mask array on
                                                     !! local model grid and input-file vertical levels.
  real, allocatable,     dimension(:)  :: z_in       !< Cell grid values for input data.
  real, allocatable,     dimension(:)  :: z_edges_in !< Cell grid edge values for input data.
  real,                  intent(out)   :: missing_value !< The missing value in the returned array.
  logical,               intent(in)    :: reentrant_x !< If true, this grid is reentrant in the x-direction
  logical,               intent(in)    :: tripolar_n !< If true, this is a northern tripolar grid
  logical,     optional, intent(in)    :: homogenize !< If present and true, horizontally homogenize data
                                                     !! to produce perfectly "flat" initial conditions
  real,        optional, intent(in)    :: m_to_Z     !< A conversion factor from meters to the units
                                                     !! of depth.  If missing, G%bathyT must be in m.
  logical,     optional, intent(in)    :: answers_2018 !< If true, use expressions that give the same
                                                     !! answers as the code did in late 2018.  Otherwise
                                                     !! add parentheses for rotational symmetry.
  logical,     optional, intent(in)    :: ongrid     !< If true, then data are assumed to have been interpolated
                                                     !! to the model horizontal grid. In this case, only
                                                     !! extrapolation is performed by this routine

  ! Local variables
  real, dimension(:,:),  allocatable   :: tr_in      !< A 2-d array for holding input data on its
                                                     !! native horizontal grid.
  real, dimension(:,:),  allocatable   :: tr_inp     !< Native horizontal grid data extended to the poles.
  real, dimension(:,:),  allocatable   :: mask_in    ! A 2-d mask for extended input grid.

  real :: PI_180
  integer :: rcode, ncid, varid, ndims, id, jd, kd, jdp
  integer :: i, j, k
  integer, dimension(4) :: start, count, dims, dim_id
  real, dimension(:,:), allocatable :: x_in, y_in
  real, dimension(:), allocatable  :: lon_in, lat_in ! The longitude and latitude in the input file
  real, dimension(:), allocatable  :: lat_inp ! The input file latitudes expanded to the pole
  real :: max_lat, min_lat, pole, max_depth, npole
  real :: roundoff  ! The magnitude of roundoff, usually ~2e-16.
  real :: add_offset, scale_factor
  logical :: add_np
  logical :: is_ongrid
  character(len=8)  :: laynum
  type(horiz_interp_type) :: Interp
  integer :: is, ie, js, je     ! compute domain indices
  integer :: isc, iec, jsc, jec ! global compute domain indices
  integer :: isg, ieg, jsg, jeg ! global extent
  integer :: isd, ied, jsd, jed ! data domain indices
  integer :: id_clock_read
  character(len=12)  :: dim_name(4)
  logical :: debug=.false.
  real :: npoints, varAvg
  real, dimension(SZI_(G),SZJ_(G)) :: lon_out, lat_out ! The longitude and latitude of points on the model grid
  real, dimension(SZI_(G),SZJ_(G)) :: tr_out, mask_out ! The tracer and mask on the model grid
  real, dimension(SZI_(G),SZJ_(G)) :: good    ! Where the data is valid, this is 1.
  real, dimension(SZI_(G),SZJ_(G)) :: fill    ! 1 where the data needs to be filled in
  real, dimension(SZI_(G),SZJ_(G)) :: tr_outf ! The tracer concentrations after Ice-9
  real, dimension(SZI_(G),SZJ_(G)) :: tr_prev ! The tracer concentrations in the layer above
  real, dimension(SZI_(G),SZJ_(G)) :: good2   ! 1 where the data is valid after Ice-9
  real, dimension(SZI_(G),SZJ_(G)) :: fill2   ! 1 for points that still need to be filled after Ice-9

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg

  id_clock_read = cpu_clock_id('(Initialize tracer from Z) read', grain=CLOCK_LOOP)

  is_ongrid = .false.
  if (present(ongrid)) is_ongrid = ongrid

  if (allocated(tr_z)) deallocate(tr_z)
  if (allocated(mask_z)) deallocate(mask_z)
  if (allocated(z_edges_in)) deallocate(z_edges_in)

  PI_180 = atan(1.0)/45.

  ! Open NetCDF file and if present, extract data and spatial coordinate information
  ! The convention adopted here requires that the data be written in (i,j,k) ordering.

  call cpu_clock_begin(id_clock_read)

  rcode = NF90_OPEN(filename, NF90_NOWRITE, ncid)
  if (rcode /= 0) call MOM_error(FATAL,"error opening file "//trim(filename)//&
                           " in hinterp_extrap")
  rcode = NF90_INQ_VARID(ncid, varnam, varid)
  if (rcode /= 0) call MOM_error(FATAL,"error finding variable "//trim(varnam)//&
                                 " in file "//trim(filename)//" in hinterp_extrap")

  rcode = NF90_INQUIRE_VARIABLE(ncid, varid, ndims=ndims, dimids=dims)
  if (rcode /= 0) call MOM_error(FATAL,'error inquiring dimensions hinterp_extrap')
  if (ndims < 3) call MOM_error(FATAL,"Variable "//trim(varnam)//" in file "// &
              trim(filename)//" has too few dimensions.")

  rcode = NF90_INQUIRE_DIMENSION(ncid, dims(1), dim_name(1), len=id)
  if (rcode /= 0) call MOM_error(FATAL,"error reading dimension 1 data for "// &
                trim(varnam)//" in file "// trim(filename)//" in hinterp_extrap")
  rcode = NF90_INQ_VARID(ncid, dim_name(1), dim_id(1))
  if (rcode /= 0) call MOM_error(FATAL,"error finding variable "//trim(dim_name(1))//&
                                 " in file "//trim(filename)//" in hinterp_extrap")
  rcode = NF90_INQUIRE_DIMENSION(ncid, dims(2), dim_name(2), len=jd)
  if (rcode /= 0) call MOM_error(FATAL,"error reading dimension 2 data for "// &
                trim(varnam)//" in file "// trim(filename)//" in hinterp_extrap")
  rcode = NF90_INQ_VARID(ncid, dim_name(2), dim_id(2))
  if (rcode /= 0) call MOM_error(FATAL,"error finding variable "//trim(dim_name(2))//&
                                 " in file "//trim(filename)//" in hinterp_extrap")
  rcode = NF90_INQUIRE_DIMENSION(ncid, dims(3), dim_name(3), len=kd)
  if (rcode /= 0) call MOM_error(FATAL,"error reading dimension 3 data for "// &
                trim(varnam)//" in file "// trim(filename)//" in hinterp_extrap")
  rcode = NF90_INQ_VARID(ncid, dim_name(3), dim_id(3))
  if (rcode /= 0) call MOM_error(FATAL,"error finding variable "//trim(dim_name(3))//&
                                 " in file "//trim(filename)//" in hinterp_extrap")

  missing_value=0.0
  rcode = NF90_GET_ATT(ncid, varid, "_FillValue", missing_value)
  if (rcode /= 0) call MOM_error(FATAL,"error finding missing value for "//&
       trim(varnam)//" in file "// trim(filename)//" in hinterp_extrap")

  rcode = NF90_GET_ATT(ncid, varid, "add_offset", add_offset)
  if (rcode /= 0) add_offset = 0.0

  rcode = NF90_GET_ATT(ncid, varid, "scale_factor", scale_factor)
  if (rcode /= 0) scale_factor = 1.0

  if (allocated(lon_in)) deallocate(lon_in)
  if (allocated(lat_in)) deallocate(lat_in)
  if (allocated(z_in)) deallocate(z_in)
  if (allocated(z_edges_in)) deallocate(z_edges_in)
  if (allocated(tr_z)) deallocate(tr_z)
  if (allocated(mask_z)) deallocate(mask_z)

  allocate(lon_in(id), lat_in(jd), z_in(kd), z_edges_in(kd+1))
  allocate(tr_z(isd:ied,jsd:jed,kd), mask_z(isd:ied,jsd:jed,kd))

  start = 1 ; count = 1 ; count(1) = id
  rcode = NF90_GET_VAR(ncid, dim_id(1), lon_in, start, count)
  if (rcode /= 0) call MOM_error(FATAL,"error reading dimension 1 values for var_name "// &
                trim(varnam)//",dim_name "//trim(dim_name(1))//" in file "// trim(filename)//" in hinterp_extrap")
  start = 1 ; count = 1 ; count(1) = jd
  rcode = NF90_GET_VAR(ncid, dim_id(2), lat_in, start, count)
  if (rcode /= 0) call MOM_error(FATAL,"error reading dimension 2 values for var_name "// &
                trim(varnam)//",dim_name "//trim(dim_name(2))//" in file "// trim(filename)//" in  hinterp_extrap")
  start = 1 ; count = 1 ; count(1) = kd
  rcode = NF90_GET_VAR(ncid, dim_id(3), z_in, start, count)
  if (rcode /= 0) call MOM_error(FATAL,"error reading dimension 3 values for var_name "// &
                trim(varnam//",dim_name "//trim(dim_name(3)))//" in file "// trim(filename)//" in  hinterp_extrap")

  call cpu_clock_end(id_clock_read)

  if (present(m_to_Z)) then ; do k=1,kd ; z_in(k) = m_to_Z * z_in(k) ; enddo ; endif

  ! extrapolate the input data to the north pole using the northern-most latitude
  add_np = .false.
  jdp = jd
  if (.not. is_ongrid) then
    max_lat = maxval(lat_in)
    if (max_lat < 90.0) then
      add_np = .true.
      jdp = jd+1
      allocate(lat_inp(jdp))
      lat_inp(1:jd) = lat_in(:)
      lat_inp(jd+1) = 90.0
      deallocate(lat_in)
      allocate(lat_in(1:jdp))
      lat_in(:) = lat_inp(:)
    endif
  endif
  ! construct level cell boundaries as the mid-point between adjacent centers

  z_edges_in(1) = 0.0
  do K=2,kd
    z_edges_in(K) = 0.5*(z_in(k-1)+z_in(k))
  enddo
  z_edges_in(kd+1) = 2.0*z_in(kd) - z_in(kd-1)

  if (is_ongrid) then
    allocate(tr_in(is:ie,js:je)) ; tr_in(:,:)=0.0
    allocate(mask_in(is:ie,js:je)) ; mask_in(:,:)=0.0
  else
    call horiz_interp_init()
    lon_in = lon_in*PI_180
    lat_in = lat_in*PI_180
    allocate(x_in(id,jdp), y_in(id,jdp))
    call meshgrid(lon_in, lat_in, x_in, y_in)
    lon_out(:,:) = G%geoLonT(:,:)*PI_180
    lat_out(:,:) = G%geoLatT(:,:)*PI_180
    allocate(tr_in(id,jd)) ; tr_in(:,:) = 0.0
    allocate(tr_inp(id,jdp)) ; tr_inp(:,:) = 0.0
    allocate(mask_in(id,jdp)) ; mask_in(:,:) = 0.0
  endif

  max_depth = maxval(G%bathyT)
  call max_across_PEs(max_depth)

  if (z_edges_in(kd+1)<max_depth) z_edges_in(kd+1)=max_depth
  roundoff = 3.0*EPSILON(missing_value)

  ! loop through each data level and interpolate to model grid.
  ! after interpolating, fill in points which will be needed
  ! to define the layers
  do k=1,kd
    write(laynum,'(I8)') k ; laynum = adjustl(laynum)
    mask_in = 0.0
    if (is_ongrid) then
      start(1) = is+G%HI%idg_offset ; start(2) = js+G%HI%jdg_offset ; start(3) = k
      count(1) = ie-is+1 ; count(2) = je-js+1; count(3) = 1
      rcode = NF90_GET_VAR(ncid,varid, tr_in, start, count)
      if (rcode /= 0) call MOM_error(FATAL,"hinterp_and_extract_from_Fie: "//&
           "error reading level "//trim(laynum)//" of variable "//&
           trim(varnam)//" in file "// trim(filename))

      do j=js,je
        do i=is,ie
          if (abs(tr_in(i,j)-missing_value) > abs(roundoff*missing_value)) then
            mask_in(i,j) = 1.0
            tr_in(i,j) = (tr_in(i,j)*scale_factor+add_offset) * conversion
          else
            tr_in(i,j) = missing_value
          endif
        enddo
      enddo

    else
      if (is_root_pe()) then
        start = 1 ; start(3) = k ; count(:) = 1 ; count(1) = id ; count(2) = jd
        rcode = NF90_GET_VAR(ncid,varid, tr_in, start, count)
        if (rcode /= 0) call MOM_error(FATAL,"hinterp_and_extract_from_Fie: "//&
             "error reading level "//trim(laynum)//" of variable "//&
             trim(varnam)//" in file "// trim(filename))

        if (add_np) then
          pole = 0.0 ; npole = 0.0
          do i=1,id
            if (abs(tr_in(i,jd)-missing_value) > abs(roundoff*missing_value)) then
              pole = pole + tr_in(i,jd)
              npole = npole + 1.0
            endif
          enddo
          if (npole > 0) then
            pole=pole/npole
          else
            pole=missing_value
          endif
          tr_inp(:,1:jd) = tr_in(:,:)
          tr_inp(:,jdp) = pole
        else
          tr_inp(:,:) = tr_in(:,:)
        endif
      endif

      call broadcast(tr_inp, id*jdp, blocking=.true.)

      do j=1,jdp
        do i=1,id
          if (abs(tr_inp(i,j)-missing_value) > abs(roundoff*missing_value)) then
            mask_in(i,j) = 1.0
            tr_inp(i,j) = (tr_inp(i,j)*scale_factor+add_offset) * conversion
          else
            tr_inp(i,j) = missing_value
          endif
        enddo
      enddo

    endif

!   call fms routine horiz_interp to interpolate input level data to model horizontal grid
    if (.not. is_ongrid) then
      if (k == 1) then
        call build_horiz_interp_weights(Interp, x_in, y_in, lon_out(is:ie,js:je), lat_out(is:ie,js:je), &
                                        interp_method='bilinear', src_modulo=.true.)
      endif

      if (debug) then
        call myStats(tr_inp,missing_value, is, ie, js, je, k,'Tracer from file')
      endif
    endif

    tr_out(:,:) = 0.0
    if (is_ongrid) then
      tr_out(is:ie,js:je)=tr_in(is:ie,js:je)
    else
      call run_horiz_interp(Interp, tr_inp, tr_out(is:ie,js:je), missing_value=missing_value)
    endif

    mask_out=1.0
    do j=js,je
      do i=is,ie
        if (abs(tr_out(i,j)-missing_value) < abs(roundoff*missing_value)) mask_out(i,j)=0.
      enddo
    enddo

    fill = 0.0; good = 0.0

    nPoints = 0 ; varAvg = 0.
    do j=js,je
      do i=is,ie
        if (mask_out(i,j) < 1.0) then
          tr_out(i,j)=missing_value
        else
          good(i,j)=1.0
          nPoints = nPoints + 1
          varAvg = varAvg + tr_out(i,j)
        endif
        if (G%mask2dT(i,j) == 1.0 .and. z_edges_in(k) <= G%bathyT(i,j) .and. mask_out(i,j) < 1.0) &
          fill(i,j)=1.0
      enddo
    enddo
    call pass_var(fill,G%Domain)
    call pass_var(good,G%Domain)

    if (debug) then
      call myStats(tr_out,missing_value, is,ie,js,je,k,'variable from horiz_interp()')
    endif

    ! Horizontally homogenize data to produce perfectly "flat" initial conditions
    if (PRESENT(homogenize)) then
      if (homogenize) then
        call sum_across_PEs(nPoints)
        call sum_across_PEs(varAvg)
        if (nPoints>0) then
          varAvg = varAvg/real(nPoints)
        endif
        tr_out(:,:) = varAvg
      endif
    endif

    ! tr_out contains input z-space data on the model grid with missing values
    ! now fill in missing values using "ICE-nine" algorithm.
    tr_outf(:,:) = tr_out(:,:)
    if (k==1) tr_prev(:,:) = tr_outf(:,:)
    good2(:,:) = good(:,:)
    fill2(:,:) = fill(:,:)

    call fill_miss_2d(tr_outf, good2, fill2, tr_prev, G, smooth=.true., answers_2018=answers_2018)
    if (debug) then
      call myStats(tr_outf, missing_value, is, ie, js, je, k, 'field from fill_miss_2d()')
    endif

    tr_z(:,:,k) = tr_outf(:,:) * G%mask2dT(:,:)
    mask_z(:,:,k) = good2(:,:) + fill2(:,:)

    tr_prev(:,:) = tr_z(:,:,k)

    if (debug) then
      call hchksum(tr_prev,'field after fill ',G%HI)
    endif

  enddo ! kd

end subroutine horiz_interp_and_extrap_tracer_record

!> Extrapolate and interpolate using a FMS time interpolation handle
subroutine horiz_interp_and_extrap_tracer_fms_id(fms_id,  Time, conversion, G, tr_z, mask_z, &
                                                 z_in, z_edges_in, missing_value, reentrant_x, &
                                                 tripolar_n, homogenize, spongeOngrid, m_to_Z, answers_2018)

  integer,               intent(in)    :: fms_id     !< A unique id used by the FMS time interpolator
  type(time_type),       intent(in)    :: Time       !< A FMS time type
  real,                  intent(in)    :: conversion !< Conversion factor for tracer.
  type(ocean_grid_type), intent(inout) :: G          !< Grid object
  real, allocatable, dimension(:,:,:)  :: tr_z       !< pointer to allocatable tracer array on local
                                                     !! model grid and native vertical levels.
  real, allocatable, dimension(:,:,:)  :: mask_z     !< pointer to allocatable tracer mask array on
                                                     !! local model grid and native vertical levels.
  real, allocatable,     dimension(:)  :: z_in       !< Cell grid values for input data.
  real, allocatable,     dimension(:)  :: z_edges_in !< Cell grid edge values for input data. (Intent out)
  real,                  intent(out)   :: missing_value !< The missing value in the returned array.
  logical,               intent(in)    :: reentrant_x !< If true, this grid is reentrant in the x-direction
  logical,               intent(in)    :: tripolar_n !< If true, this is a northern tripolar grid
  logical,     optional, intent(in)    :: homogenize !< If present and true, horizontally homogenize data
                                                     !! to produce perfectly "flat" initial conditions
  logical,     optional, intent(in)    :: spongeOngrid !< If present and true, the sponge data are on the model grid
  real,        optional, intent(in)    :: m_to_Z     !< A conversion factor from meters to the units
                                                     !! of depth.  If missing, G%bathyT must be in m.
  logical,     optional, intent(in)    :: answers_2018 !< If true, use expressions that give the same
                                                     !! answers as the code did in late 2018.  Otherwise
                                                     !! add parentheses for rotational symmetry.

  ! Local variables
  real, dimension(:,:),  allocatable   :: tr_in      !< A 2-d array for holding input data on its
                                                     !! native horizontal grid.
  real, dimension(:,:),  allocatable   :: tr_inp     !< Native horizontal grid data extended to the poles.
  real, dimension(:,:,:), allocatable  :: data_in    !< A buffer for storing the full 3-d time-interpolated array
                                                     !! on the original grid
  real, dimension(:,:),  allocatable   :: mask_in    !< A 2-d mask for extended input grid.

  real :: PI_180
  integer :: rcode, ncid, varid, ndims, id, jd, kd, jdp
  integer :: i,j,k
  integer, dimension(4) :: start, count, dims, dim_id
  real, dimension(:,:), allocatable :: x_in, y_in
  real, dimension(:), allocatable  :: lon_in, lat_in ! The longitude and latitude in the input file
  real, dimension(:), allocatable  :: lat_inp ! The input file latitudes expanded to the pole
  real :: max_lat, min_lat, pole, max_depth, npole
  real :: roundoff  ! The magnitude of roundoff, usually ~2e-16.
  logical :: add_np
  character(len=8)  :: laynum
  type(horiz_interp_type) :: Interp
  type(axistype), dimension(4) :: axes_data
  integer :: is, ie, js, je     ! compute domain indices
  integer :: isc,iec,jsc,jec    ! global compute domain indices
  integer :: isg, ieg, jsg, jeg ! global extent
  integer :: isd, ied, jsd, jed ! data domain indices
  integer :: id_clock_read
  integer, dimension(4) :: fld_sz
  character(len=12)  :: dim_name(4)
  logical :: debug=.false.
  logical :: spongeDataOngrid
  logical :: ans_2018
  real :: npoints, varAvg
  real, dimension(SZI_(G),SZJ_(G)) :: lon_out, lat_out ! The longitude and latitude of points on the model grid
  real, dimension(SZI_(G),SZJ_(G)) :: tr_out, mask_out ! The tracer and mask on the model grid
  real, dimension(SZI_(G),SZJ_(G)) :: good    ! Where the data is valid, this is 1.
  real, dimension(SZI_(G),SZJ_(G)) :: fill    ! 1 where the data needs to be filled in
  real, dimension(SZI_(G),SZJ_(G)) :: tr_outf ! The tracer concentrations after Ice-9
  real, dimension(SZI_(G),SZJ_(G)) :: tr_prev ! The tracer concentrations in the layer above
  real, dimension(SZI_(G),SZJ_(G)) :: good2   ! 1 where the data is valid after Ice-9
  real, dimension(SZI_(G),SZJ_(G)) :: fill2   ! 1 for points that still need to be filled after Ice-9
  integer :: turns

  turns = G%HI%turns

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg

  id_clock_read = cpu_clock_id('(Initialize tracer from Z) read', grain=CLOCK_LOOP)

  PI_180 = atan(1.0)/45.

  ans_2018 = .true.;if (present(answers_2018)) ans_2018 = answers_2018

  ! Open NetCDF file and if present, extract data and spatial coordinate information
  ! The convention adopted here requires that the data be written in (i,j,k) ordering.

  call cpu_clock_begin(id_clock_read)

  call get_external_field_info(fms_id, size=fld_sz, axes=axes_data, missing=missing_value)
  if (allocated(lon_in)) deallocate(lon_in)
  if (allocated(lat_in)) deallocate(lat_in)
  if (allocated(z_in)) deallocate(z_in)
  if (allocated(z_edges_in)) deallocate(z_edges_in)
  if (allocated(tr_z)) deallocate(tr_z)
  if (allocated(mask_z)) deallocate(mask_z)

  id = fld_sz(1) ; jd  = fld_sz(2) ; kd = fld_sz(3)

  spongeDataOngrid = .false.
  if (PRESENT(spongeOngrid)) spongeDataOngrid = spongeOngrid
  if (.not. spongeDataOngrid) then
    allocate(lon_in(id), lat_in(jd))
    call get_axis_data(axes_data(1), lon_in)
    call get_axis_data(axes_data(2), lat_in)
  endif

  allocate(z_in(kd), z_edges_in(kd+1))

  allocate(tr_z(isd:ied,jsd:jed,kd), mask_z(isd:ied,jsd:jed,kd))

  call get_axis_data(axes_data(3), z_in)

  if (present(m_to_Z)) then ; do k=1,kd ; z_in(k) = m_to_Z * z_in(k) ; enddo ; endif

  call cpu_clock_end(id_clock_read)

  if (.not. spongeDataOngrid) then
    ! Extrapolate the input data to the north pole using the northerm-most latitude.
    max_lat = maxval(lat_in)
    add_np = .false.
    if (max_lat < 90.0) then
      add_np = .true.
      jdp = jd+1
      allocate(lat_inp(jdp))
      lat_inp(1:jd) = lat_in(:)
      lat_inp(jd+1) = 90.0
      deallocate(lat_in)
      allocate(lat_in(1:jdp))
      lat_in(:) = lat_inp(:)
    else
      jdp = jd
    endif
    call horiz_interp_init()
    lon_in = lon_in*PI_180
    lat_in = lat_in*PI_180
    allocate(x_in(id,jdp), y_in(id,jdp))
    call meshgrid(lon_in, lat_in, x_in, y_in)
    lon_out(:,:) = G%geoLonT(:,:)*PI_180
    lat_out(:,:) = G%geoLatT(:,:)*PI_180
    allocate(data_in(id,jd,kd)) ; data_in(:,:,:)=0.0
    allocate(tr_in(id,jd)) ; tr_in(:,:)=0.0
    allocate(tr_inp(id,jdp)) ; tr_inp(:,:)=0.0
    allocate(mask_in(id,jdp)) ; mask_in(:,:)=0.0
  else
    allocate(data_in(isd:ied,jsd:jed,kd))
  endif
  ! construct level cell boundaries as the mid-point between adjacent centers
  z_edges_in(1) = 0.0
  do k=2,kd
    z_edges_in(k) = 0.5*(z_in(k-1)+z_in(k))
  enddo
  z_edges_in(kd+1) = 2.0*z_in(kd) - z_in(kd-1)


  max_depth = maxval(G%bathyT)
  call max_across_PEs(max_depth)

  if (z_edges_in(kd+1)<max_depth) z_edges_in(kd+1)=max_depth

  !  roundoff = 3.0*EPSILON(missing_value)
  roundoff = 1.e-4

  if (.not.spongeDataOngrid) then
    if (is_root_pe()) &
      call time_interp_external(fms_id, Time, data_in, verbose=.true., turns=turns)
    ! Loop through each data level and interpolate to model grid.
    ! After interpolating, fill in points which will be needed to define the layers.
    do k=1,kd
      write(laynum,'(I8)') k ; laynum = adjustl(laynum)
      if (is_root_pe()) then
        tr_in(1:id,1:jd) = data_in(1:id,1:jd,k)
        if (add_np) then
          pole = 0.0 ; npole = 0.0
          do i=1,id
            if (abs(tr_in(i,jd)-missing_value) > abs(roundoff*missing_value)) then
              pole = pole + tr_in(i,jd)
              npole = npole+1.0
            endif
          enddo
          if (npole > 0) then
            pole = pole / npole
          else
            pole = missing_value
          endif
          tr_inp(:,1:jd) = tr_in(:,:)
          tr_inp(:,jdp) = pole
        else
          tr_inp(:,:) = tr_in(:,:)
        endif
      endif

      call broadcast(tr_inp, id*jdp, blocking=.true.)

      mask_in(:,:) = 0.0

      do j=1,jdp ; do i=1,id
        if (abs(tr_inp(i,j)-missing_value) > abs(roundoff*missing_value)) then
          mask_in(i,j) = 1.0
          tr_inp(i,j) = tr_inp(i,j) * conversion
        else
          tr_inp(i,j) = missing_value
        endif
      enddo ; enddo

      ! call fms routine horiz_interp to interpolate input level data to model horizontal grid
      if (k == 1) then
        call build_horiz_interp_weights(Interp, x_in, y_in, lon_out(is:ie,js:je), lat_out(is:ie,js:je), &
                                        interp_method='bilinear', src_modulo=.true.)
      endif

      if (debug) then
        call myStats(tr_in, missing_value, 1, id, 1, jd, k, 'Tracer from file')
      endif

      tr_out(:,:) = 0.0

      call run_horiz_interp(Interp, tr_inp, tr_out(is:ie,js:je), missing_value=missing_value)

      mask_out(:,:) = 1.0
      do j=js,je ; do i=is,ie
        if (abs(tr_out(i,j)-missing_value) < abs(roundoff*missing_value)) mask_out(i,j) = 0.
      enddo ; enddo

      fill(:,:) = 0.0 ; good(:,:) = 0.0

      nPoints = 0 ; varAvg = 0.
      do j=js,je ; do i=is,ie
        if (mask_out(i,j) < 1.0) then
          tr_out(i,j) = missing_value
        else
          good(i,j) = 1.0
          nPoints = nPoints + 1
          varAvg = varAvg + tr_out(i,j)
        endif
        if ((G%mask2dT(i,j) == 1.0) .and. (z_edges_in(k) <= G%bathyT(i,j)) .and. &
            (mask_out(i,j) < 1.0)) &
          fill(i,j) = 1.0
      enddo ;  enddo
      call pass_var(fill, G%Domain)
      call pass_var(good, G%Domain)

      if (debug) then
        call myStats(tr_out, missing_value, is, ie, js, je, k, 'variable from horiz_interp()')
      endif

      ! Horizontally homogenize data to produce perfectly "flat" initial conditions
      if (PRESENT(homogenize)) then ; if (homogenize) then
        call sum_across_PEs(nPoints)
        call sum_across_PEs(varAvg)
        if (nPoints>0) then
          varAvg = varAvg/real(nPoints)
        endif
        tr_out(:,:) = varAvg
      endif ; endif

      ! tr_out contains input z-space data on the model grid with missing values
      ! now fill in missing values using "ICE-nine" algorithm.
      tr_outf(:,:) = tr_out(:,:)
      if (k==1) tr_prev(:,:) = tr_outf(:,:)
      good2(:,:) = good(:,:)
      fill2(:,:) = fill(:,:)

      call fill_miss_2d(tr_outf, good2, fill2, tr_prev, G, smooth=.true., answers_2018=answers_2018)

!     if (debug) then
!       call hchksum(tr_outf, 'field from fill_miss_2d ', G%HI)
!     endif

!     call myStats(tr_outf, missing_value, is, ie, js, je, k, 'field from fill_miss_2d()')

      tr_z(:,:,k) = tr_outf(:,:)*G%mask2dT(:,:)
      mask_z(:,:,k) = good2(:,:) + fill2(:,:)
      tr_prev(:,:) = tr_z(:,:,k)

      if (debug) then
        call hchksum(tr_prev,'field after fill ',G%HI)
      endif

    enddo ! kd
  else
      call time_interp_external(fms_id, Time, data_in, verbose=.true., turns=turns)
      do k=1,kd
        do j=js,je
          do i=is,ie
            tr_z(i,j,k)=data_in(i,j,k)
            if (.not. ans_2018) mask_z(i,j,k) = 1.
            if (abs(tr_z(i,j,k)-missing_value) < abs(roundoff*missing_value)) mask_z(i,j,k) = 0.
          enddo
        enddo
      enddo
  endif

end subroutine horiz_interp_and_extrap_tracer_fms_id

!> Create a 2d-mesh of grid coordinates from 1-d arrays.
subroutine meshgrid(x, y, x_T, y_T)
  real, dimension(:),                   intent(in)    :: x  !< input 1-dimensional vector
  real, dimension(:),                   intent(in)    :: y  !< input 1-dimensional vector
  real, dimension(size(x,1),size(y,1)), intent(inout) :: x_T !< output 2-dimensional array
  real, dimension(size(x,1),size(y,1)), intent(inout) :: y_T !< output 2-dimensional array

  integer :: ni, nj, i, j

  ni = size(x,1) ; nj = size(y,1)

  do j=1,nj ; do i=1,ni
    x_T(i,j) = x(i)
  enddo ; enddo

  do j=1,nj ; do i=1,ni
    y_T(i,j) = y(j)
  enddo ; enddo

end subroutine meshgrid

end module MOM_horizontal_regridding
