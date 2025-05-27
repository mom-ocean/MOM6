!> Horizontal interpolation
module MOM_horizontal_regridding

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging,     only : hchksum
use MOM_coms,          only : max_across_PEs, min_across_PEs, sum_across_PEs, broadcast
use MOM_coms,          only : reproducing_sum
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_LOOP
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_error_handler, only : MOM_get_verbosity
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_interpolate,   only : time_interp_external
use MOM_interp_infra,  only : run_horiz_interp, build_horiz_interp_weights
use MOM_interp_infra,  only : horiz_interp_type, horizontal_interp_init
use MOM_interp_infra,  only : get_external_field_info
use MOM_interp_infra,  only : external_field
use MOM_time_manager,  only : time_type
use MOM_io,            only : axis_info, get_axis_info, get_var_axes_info, MOM_read_data
use MOM_io,            only : read_attribute, read_variable

implicit none ; private

#include <MOM_memory.h>

public :: horiz_interp_and_extrap_tracer, myStats, homogenize_field

!> Extrapolate and interpolate data
interface horiz_interp_and_extrap_tracer
  module procedure horiz_interp_and_extrap_tracer_record
  module procedure horiz_interp_and_extrap_tracer_fms_id
end interface

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.
! The functions in this module work with variables with arbitrary units, in which case the
! arbitrary rescaled units are indicated with [A ~> a], while the unscaled units are just [a].

contains

!> Write to the terminal some basic statistics about the k-th level of an array
subroutine myStats(array, missing, G, k, mesg, unscale, full_halo)
  type(ocean_grid_type), intent(in) :: G     !< Ocean grid type
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(in) :: array !< input array in arbitrary units [A ~> a]
  real,                  intent(in) :: missing !< missing value in arbitrary units [A ~> a]
  integer,               intent(in) :: k     !< Level to calculate statistics for
  character(len=*),      intent(in) :: mesg  !< Label to use in message
  real,        optional, intent(in) :: unscale !< A scaling factor for output that countacts
                                             !! any internal dimesional scaling [a A-1 ~> 1]
  logical,     optional, intent(in) :: full_halo !< If present and true, test values on the whole
                                             !! array rather than just the computational domain.
  ! Local variables
  real :: minA ! Minimum value in the array in the arbitrary units of the input array [A ~> a]
  real :: maxA ! Maximum value in the array in the arbitrary units of the input array [A ~> a]
  real :: scl  ! A factor for undoing any scaling of the array statistics for output [a A-1 ~> 1]
  integer :: i, j, is, ie, js, je
  logical :: found
  character(len=120) :: lMesg

  scl = 1.0 ; if (present(unscale)) scl = unscale
  minA = 9.E24 / scl ; maxA = -9.E24 / scl ; found = .false.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  if (present(full_halo)) then ; if (full_halo) then
    is = G%isd ; ie = G%ied ; js = G%jsd ; je = G%jed
  endif ; endif

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
    write(lMesg(1:120),'(2(a,es12.4),a,i3,1x,a)') &
         'init_from_Z: min=',minA*scl,' max=',maxA*scl,' Level=',k,trim(mesg)
    call MOM_mesg(lMesg,2)
  endif

end subroutine myStats

!> Use ICE-9 algorithm to populate points (fill=1) with valid data (good=1).  If no information
!! is available, use a previous guess (prev). Optionally (smooth) blend the filled points to
!! achieve a more desirable result.
subroutine fill_miss_2d(aout, good, fill, prev, G, acrit, num_pass, relc, debug, answer_date)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(inout) :: aout !< The array with missing values to fill [arbitrary]
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(in)    :: good !< Valid data mask for incoming array
                                               !! (1==good data; 0==missing data) [nondim].
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(in)    :: fill !< Same shape array of points which need
                                               !! filling (1==fill;0==dont fill) [nondim]
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(in)    :: prev !< First guess where isolated holes exist [arbitrary]
  real,                  intent(in)    :: acrit !< A minimal value for deltas between iterations that
                                               !! determines when the smoothing has converged [arbitrary].
  integer,     optional, intent(in)    :: num_pass !< The maximum number of iterations
  real,        optional, intent(in)    :: relc !< A relaxation coefficient for Laplacian [nondim]
  logical,     optional, intent(in)    :: debug !< If true, write verbose debugging messages.
  integer,     optional, intent(in)    :: answer_date !< The vintage of the expressions in the code.
                                                !! Dates before 20190101 give the same  answers
                                                !! as the code did in late 2018, while later versions
                                                !! add parentheses for rotational symmetry.

  real, dimension(SZI_(G),SZJ_(G)) :: a_filled ! The aout with missing values filled in [arbitrary]
  real, dimension(SZI_(G),SZJ_(G)) :: a_chg    ! The change in aout due to an iteration of smoothing [arbitrary]
  real, dimension(SZI_(G),SZJ_(G)) :: fill_pts ! 1 for points that still need to be filled [nondim]
  real, dimension(SZI_(G),SZJ_(G)) :: good_    ! The values that are valid for the current iteration [nondim]
  real, dimension(SZI_(G),SZJ_(G)) :: good_new ! The values of good_ to use for the next iteration [nondim]

  real    :: east, west, north, south ! Valid neighboring values or 0 for invalid values [arbitrary]
  real    :: ge, gw, gn, gs  ! Flags set to 0 or 1 indicating which neighbors have valid values [nondim]
  real    :: ngood     ! The number of valid values in neighboring points [nondim]
  real    :: nfill     ! The remaining number of points to fill [nondim]
  real    :: nfill_prev ! The previous value of nfill [nondim]
  character(len=256) :: mesg  ! The text of an error message
  integer :: i, j, k
  integer, parameter :: num_pass_default = 10000
  real, parameter :: relc_default = 0.25  ! The default relaxation coefficient [nondim]

  integer :: npass  ! The maximum number of passes of the Laplacian smoother
  integer :: is, ie, js, je
  real    :: relax_coeff  ! The grid-scale Laplacian relaxation coefficient per timestep [nondim]
  real    :: ares   ! The maximum magnitude change in aout [A]
  logical :: debug_it, ans_2018

  debug_it=.false.
  if (PRESENT(debug)) debug_it=debug

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  npass = num_pass_default
  if (PRESENT(num_pass)) npass = num_pass

  relax_coeff = relc_default
  if (PRESENT(relc)) relax_coeff = relc

  ans_2018 = .true. ; if (PRESENT(answer_date)) ans_2018 = (answer_date < 20190101)

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

    if (nfill == nfill_prev) then
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
  do k=1,npass
    call pass_var(aout,G%Domain)

    a_chg(:,:) = 0.0
    if (ans_2018) then
      do j=js,je ; do i=is,ie
        if (fill(i,j) == 1) then
          east = max(good(i+1,j),fill(i+1,j)) ; west = max(good(i-1,j),fill(i-1,j))
          north = max(good(i,j+1),fill(i,j+1)) ; south = max(good(i,j-1),fill(i,j-1))
          a_chg(i,j) = relax_coeff*(south*aout(i,j-1)+north*aout(i,j+1) + &
                                    west*aout(i-1,j)+east*aout(i+1,j) - &
                                   (south+north+west+east)*aout(i,j))
        endif
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        if (fill(i,j) == 1) then
          ge = max(good(i+1,j),fill(i+1,j)) ; gw = max(good(i-1,j),fill(i-1,j))
          gn = max(good(i,j+1),fill(i,j+1)) ; gs = max(good(i,j-1),fill(i,j-1))
          a_chg(i,j) = relax_coeff*( ((gs*aout(i,j-1) + gn*aout(i,j+1)) + &
                                      (gw*aout(i-1,j) + ge*aout(i+1,j))) - &
                                     ((gs + gn) + (gw + ge))*aout(i,j) )
        endif
      enddo ; enddo
    endif

    ares = 0.0
    do j=js,je ; do i=is,ie
      aout(i,j) = a_chg(i,j) + aout(i,j)
      ares = max(ares, abs(a_chg(i,j)))
    enddo ; enddo
    call max_across_PEs(ares)
    if (ares <= acrit) exit
  enddo

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
subroutine horiz_interp_and_extrap_tracer_record(filename, varnam, recnum, G, tr_z, mask_z, &
                                                 z_in, z_edges_in, missing_value, scale, &
                                                 homogenize, m_to_Z, answers_2018, ongrid, tr_iter_tol, answer_date)

  character(len=*),      intent(in)    :: filename   !< Path to file containing tracer to be
                                                     !! interpolated.
  character(len=*),      intent(in)    :: varnam     !< Name of tracer in file.
  integer,               intent(in)    :: recnum     !< Record number of tracer to be read.
  type(ocean_grid_type), intent(inout) :: G          !< Grid object
  real, allocatable, dimension(:,:,:), intent(out) :: tr_z
                                                     !< Allocatable tracer array on the horizontal
                                                     !! model grid and input-file vertical levels
                                                     !! in arbitrary units [A ~> a]
  real, allocatable, dimension(:,:,:), intent(out) :: mask_z
                                                     !< Allocatable tracer mask array on the horizontal
                                                     !! model grid and input-file vertical levels [nondim]
  real, allocatable, dimension(:), intent(out) :: z_in
                                                     !< Cell grid values for input data [Z ~> m]
  real, allocatable, dimension(:), intent(out) :: z_edges_in
                                                     !< Cell grid edge values for input data [Z ~> m]
  real,                  intent(out)   :: missing_value !< The missing value in the returned array, scaled
                                                     !! to avoid accidentally having valid values match
                                                     !! missing values in the same units as tr_z [A ~> a]
  real,                  intent(in)    :: scale      !< Scaling factor for tracer into the internal
                                                     !! units of the model for the units in the file [A a-1 ~> 1]
  logical,     optional, intent(in)    :: homogenize !< If present and true, horizontally homogenize data
                                                     !! to produce perfectly "flat" initial conditions
  real,        optional, intent(in)    :: m_to_Z     !< A conversion factor from meters to the units
                                                     !! of depth [Z m-1 ~> 1].  If missing, G%bathyT must be in m.
  logical,     optional, intent(in)    :: answers_2018 !< If true, use expressions that give the same
                                                     !! answers as the code did in late 2018.  Otherwise
                                                     !! add parentheses for rotational symmetry.
  logical,     optional, intent(in)    :: ongrid     !< If true, then data are assumed to have been interpolated
                                                     !! to the model horizontal grid. In this case, only
                                                     !! extrapolation is performed by this routine
  real,        optional, intent(in)    :: tr_iter_tol !< The tolerance for changes in tracer concentrations
                                                     !! between smoothing iterations that determines when to
                                                     !! stop iterating in the same units as tr_z [A ~> a]
  integer,     optional, intent(in)    :: answer_date !< The vintage of the expressions in the code.
                                                     !! Dates before 20190101 give the same  answers
                                                     !! as the code did in late 2018, while later versions
                                                     !! add parentheses for rotational symmetry.

  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums
  real, dimension(:,:),  allocatable   :: tr_in      !< A 2-d array for holding input data on its
                                                     !! native horizontal grid, with units that change
                                                     !! as the input data is interpreted [a] then [A ~> a]
  real, dimension(:,:,:), allocatable  :: tr_in_full  !< A 3-d array for holding input data on the
                                                     !! model horizontal grid, with units that change
                                                     !! as the input data is interpreted [a] then [A ~> a]
  real, dimension(:,:),  allocatable   :: tr_inp     !< Native horizontal grid data extended to the poles
                                                     !! with units that change as the input data is
                                                     !! interpreted [a] then [A ~> a]
  real, dimension(:,:),  allocatable   :: mask_in    ! A 2-d mask for extended input grid [nondim]

  real :: PI_180  ! A conversion factor from degrees to radians [radians degree-1]
  integer :: id, jd, kd, jdp ! Input dataset data sizes
  integer :: i, j, k
  integer, dimension(4) :: start, count
  real, dimension(:,:), allocatable :: x_in ! Input file longitudes [radians]
  real, dimension(:,:), allocatable :: y_in ! Input file latitudes [radians]
  real, dimension(:), allocatable :: lon_in ! The longitudes in the input file [degreesE] then [radians]
  real, dimension(:), allocatable :: lat_in ! The latitudes in the input file [degreesN] then [radians]
  real, dimension(:), allocatable :: lat_inp ! The input file latitudes expanded to the pole [degreesN] then [radians]
  real :: max_lat   ! The maximum latitude on the input grid [degreesN]
  real :: pole      ! The sum of tracer values at the pole [a]
  real :: max_depth ! The maximum depth of the ocean [Z ~> m]
  real :: npole     ! The number of points contributing to the pole value [nondim]
  real :: missing_val_in ! The missing value in the input field [a]
  real :: roundoff  ! The magnitude of roundoff, usually ~2e-16 [nondim]
  real :: add_offset, scale_factor  ! File-specific conversion factors [a] or [nondim]
  integer :: ans_date           ! The vintage of the expressions and order of arithmetic to use
  logical :: found_attr
  logical :: add_np
  logical :: is_ongrid
  type(horiz_interp_type) :: Interp
  type(axis_info), dimension(4) :: axes_info ! Axis information used for regridding
  integer :: is, ie, js, je     ! compute domain indices
  integer :: isg, ieg, jsg, jeg ! global extent
  integer :: isd, ied, jsd, jed ! data domain indices
  integer :: id_clock_read
  logical :: debug=.false.
  real :: I_scale               ! The inverse of the scale factor for diagnostic output [a A-1 ~> 1]
  real :: dtr_iter_stop         ! The tolerance for changes in tracer concentrations between smoothing
                                ! iterations that determines when to stop iterating [A ~> a]
  real, dimension(SZI_(G),SZJ_(G)) :: lon_out ! The longitude of points on the model grid [radians]
  real, dimension(SZI_(G),SZJ_(G)) :: lat_out ! The latitude of points on the model grid [radians]
  real, dimension(SZI_(G),SZJ_(G)) :: tr_out  ! The tracer on the model grid [A ~> a]
  real, dimension(SZI_(G),SZJ_(G)) :: mask_out ! The mask on the model grid [nondim]
  real, dimension(SZI_(G),SZJ_(G)) :: good    ! Where the data is valid, this is 1 [nondim]
  real, dimension(SZI_(G),SZJ_(G)) :: fill    ! 1 where the data needs to be filled in [nondim]
  real, dimension(SZI_(G),SZJ_(G)) :: tr_outf ! The tracer concentrations after Ice-9 [A ~> a]
  real, dimension(SZI_(G),SZJ_(G)) :: tr_prev ! The tracer concentrations in the layer above [A ~> a]
  real, dimension(SZI_(G),SZJ_(G)) :: good2   ! 1 where the data is valid after Ice-9 [nondim]
  real, dimension(SZI_(G),SZJ_(G)) :: fill2   ! 1 for points that still need to be filled after Ice-9 [nondim]

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg

  id_clock_read = cpu_clock_id('(Initialize tracer from Z) read', grain=CLOCK_LOOP)

  is_ongrid = .false.
  if (present(ongrid)) is_ongrid = ongrid

  dtr_iter_stop = 1.0e-3*scale
  if (present(tr_iter_tol)) dtr_iter_stop = tr_iter_tol

  I_scale = 1.0 / scale

  PI_180 = atan(1.0)/45.

  ans_date = 20181231
  if (present(answers_2018)) then ; if (.not.answers_2018) ans_date = 20190101 ; endif
  if (present(answer_date)) ans_date = answer_date

  ! Open NetCDF file and if present, extract data and spatial coordinate information
  ! The convention adopted here requires that the data be written in (i,j,k) ordering.

  call cpu_clock_begin(id_clock_read)

  ! A note by MJH copied from elsewhere suggests that this code may be using the model connectivity
  ! (e.g., reentrant or tripolar) but should use the dataset's connectivity instead.

  call get_var_axes_info(trim(filename), trim(varnam), axes_info)

  if (allocated(z_in)) deallocate(z_in)
  if (allocated(z_edges_in)) deallocate(z_edges_in)
  if (allocated(tr_z)) deallocate(tr_z)
  if (allocated(mask_z)) deallocate(mask_z)

  call get_axis_info(axes_info(1),ax_size=id)
  call get_axis_info(axes_info(2),ax_size=jd)
  call get_axis_info(axes_info(3),ax_size=kd)

  allocate(lon_in(id), lat_in(jd), z_in(kd), z_edges_in(kd+1))
  allocate(tr_z(isd:ied,jsd:jed,kd), source=0.0)
  allocate(mask_z(isd:ied,jsd:jed,kd), source=0.0)

  call get_axis_info(axes_info(1),ax_data=lon_in)
  call get_axis_info(axes_info(2),ax_data=lat_in)
  call get_axis_info(axes_info(3),ax_data=z_in)

  call cpu_clock_end(id_clock_read)

  if (present(m_to_Z)) then ; do k=1,kd ; z_in(k) = m_to_Z * z_in(k) ; enddo ; endif

  add_np = .false.
  jdp = jd
  if (.not. is_ongrid) then
    max_lat = maxval(lat_in)
    if (max_lat < 90.0) then
      ! Extrapolate the input data to the north pole using the northern-most latitude.
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

  ! Set the I/O attributes
  call read_attribute(trim(filename), "_FillValue", missing_val_in, &
                      varname=trim(varnam), found=found_attr)
  if (.not. found_attr) call MOM_error(FATAL, &
    "error finding missing value for " // trim(varnam) // &
    " in file " // trim(filename) // " in hinterp_extrap")
  missing_value = scale * missing_val_in

  call read_attribute(trim(filename), "scale_factor", scale_factor, &
                      varname=trim(varnam), found=found_attr)
  if (.not. found_attr) scale_factor = 1.

  call read_attribute(trim(filename), "add_offset", add_offset, &
                      varname=trim(varnam), found=found_attr)
  if (.not. found_attr) add_offset = 0.

  z_edges_in(1) = 0.0
  do K=2,kd
    z_edges_in(K) = 0.5*(z_in(k-1)+z_in(k))
  enddo
  z_edges_in(kd+1) = 2.0*z_in(kd) - z_in(kd-1)

  if (is_ongrid) then
    allocate(tr_in(is:ie,js:je), source=0.0)
    allocate(tr_in_full(is:ie,js:je,kd), source=0.0)
    allocate(mask_in(is:ie,js:je), source=0.0)
  else
    call horizontal_interp_init()
    lon_in = lon_in*PI_180
    lat_in = lat_in*PI_180
    allocate(x_in(id,jdp), y_in(id,jdp))
    call meshgrid(lon_in, lat_in, x_in, y_in)
    lon_out(:,:) = G%geoLonT(:,:)*PI_180
    lat_out(:,:) = G%geoLatT(:,:)*PI_180
    allocate(tr_in(id,jd), source=0.0)
    allocate(tr_inp(id,jdp), source=0.0)
    allocate(mask_in(id,jdp), source=0.0)
  endif

  max_depth = maxval(G%bathyT(:,:)) + G%Z_ref
  call max_across_PEs(max_depth)

  if (z_edges_in(kd+1) < max_depth) z_edges_in(kd+1) = max_depth
  roundoff = 3.0*EPSILON(missing_val_in)

  ! Loop through each data level and interpolate to model grid.
  ! After interpolating, fill in points which will be needed to define the layers.

  if (is_ongrid) then
    start(1) = is+G%HI%idg_offset ; start(2) = js+G%HI%jdg_offset ; start(3) = 1
    count(1) = ie-is+1 ; count(2) = je-js+1 ; count(3) = kd ; start(4) = 1 ; count(4) = 1
    call MOM_read_data(trim(filename), trim(varnam), tr_in_full, start, count, G%Domain)
  endif

  do k=1,kd
    mask_in(:,:)  = 0.0
    tr_out(:,:) = 0.0

    if (is_ongrid) then
      tr_in(is:ie,js:je) = tr_in_full(is:ie,js:je,k)
      do j=js,je
        do i=is,ie
          if (abs(tr_in(i,j)-missing_val_in) > abs(roundoff*missing_val_in)) then
            mask_in(i,j) = 1.0
            tr_in(i,j) = (tr_in(i,j)*scale_factor+add_offset) * scale
          else
            tr_in(i,j) = missing_value
          endif
        enddo
      enddo

      tr_out(is:ie,js:je) = tr_in(is:ie,js:je)

    else  ! .not.is_ongrid

      start(:) = 1 ; start(3) = k
      count(:) = 1 ; count(1) = id ; count(2) = jd
      call read_variable(trim(filename), trim(varnam), tr_in, start=start, nread=count)

      if (is_root_pe()) then
        if (add_np) then
          pole = 0.0 ; npole = 0.0
          do i=1,id
            if (abs(tr_in(i,jd)-missing_val_in) > abs(roundoff*missing_val_in)) then
              pole = pole + tr_in(i,jd)
              npole = npole + 1.0
            endif
          enddo
          if (npole > 0) then
            pole = pole / npole
          else
            pole = missing_val_in
          endif
          tr_inp(:,1:jd) = tr_in(:,:)
          tr_inp(:,jdp) = pole
        else
          tr_inp(:,:) = tr_in(:,:)
        endif
      endif

      call broadcast(tr_inp, id*jdp, blocking=.true.)

      do j=1,jdp ; do i=1,id
        if (abs(tr_inp(i,j)-missing_val_in) > abs(roundoff*missing_val_in)) then
          mask_in(i,j) = 1.0
          tr_inp(i,j) = (tr_inp(i,j)*scale_factor+add_offset) * scale
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
        call myStats(tr_inp, missing_value, G, k, 'Tracer from file', unscale=I_scale, full_halo=.true.)
      endif

      call run_horiz_interp(Interp, tr_inp, tr_out(is:ie,js:je), missing_value=missing_value)
    endif  ! End of .not.is_ongrid

    mask_out(:,:) = 1.0
    do j=js,je ; do i=is,ie
      if (abs(tr_out(i,j)-missing_value) < abs(roundoff*missing_value)) mask_out(i,j) = 0.
    enddo ; enddo

    fill(:,:) = 0.0 ; good(:,:) = 0.0

    do j=js,je ; do i=is,ie
      if (mask_out(i,j) < 1.0) then
        tr_out(i,j) = missing_value
      else
        good(i,j) = 1.0
      endif
      if ((G%mask2dT(i,j) == 1.0) .and. (z_edges_in(k) <= G%bathyT(i,j) + G%Z_ref) .and. &
          (mask_out(i,j) < 1.0)) &
        fill(i,j) = 1.0
    enddo ; enddo

    call pass_var(fill, G%Domain)
    call pass_var(good, G%Domain)

    if (debug) then
      call myStats(tr_out, missing_value, G, k, 'variable from horiz_interp()', unscale=I_scale)
    endif

    ! Horizontally homogenize data to produce perfectly "flat" initial conditions
    if (PRESENT(homogenize)) then ; if (homogenize) then
      call homogenize_field(tr_out, G, tmp_scale=I_scale, weights=mask_out, answer_date=answer_date)
    endif ; endif

    ! tr_out contains input z-space data on the model grid with missing values
    ! now fill in missing values using "ICE-nine" algorithm.
    tr_outf(:,:) = tr_out(:,:)
    if (k==1) tr_prev(:,:) = tr_outf(:,:)
    good2(:,:) = good(:,:)
    fill2(:,:) = fill(:,:)

    call fill_miss_2d(tr_outf, good2, fill2, tr_prev, G, dtr_iter_stop, answer_date=ans_date)
    if (debug) then
      call myStats(tr_outf, missing_value, G, k, 'field from fill_miss_2d()', unscale=I_scale)
    endif

    tr_z(:,:,k) = tr_outf(:,:) * G%mask2dT(:,:)
    mask_z(:,:,k) = good2(:,:) + fill2(:,:)

    tr_prev(:,:) = tr_z(:,:,k)

    if (debug) then
      call hchksum(tr_prev, 'field after fill ', G%HI, unscale=I_scale)
    endif

  enddo ! kd

  if (allocated(lat_inp)) deallocate(lat_inp)
  deallocate(tr_in)
  if (allocated(tr_inp)) deallocate(tr_inp)
  if (allocated(tr_in_full)) deallocate(tr_in_full)

end subroutine horiz_interp_and_extrap_tracer_record

!> Extrapolate and interpolate using a FMS time interpolation handle
subroutine horiz_interp_and_extrap_tracer_fms_id(field, Time, G, tr_z, mask_z, &
                                                 z_in, z_edges_in, missing_value, scale, &
                                                 homogenize, spongeOngrid, m_to_Z, &
                                                 answers_2018, tr_iter_tol, answer_date, &
                                                 axes)

  type(external_field), intent(in)     :: field      !< Handle for the time interpolated field
  type(time_type),       intent(in)    :: Time       !< A FMS time type
  type(ocean_grid_type), intent(inout) :: G          !< Grid object
  real, allocatable, dimension(:,:,:), intent(out) :: tr_z
                                                     !< Allocatable tracer array on the horizontal
                                                     !! model grid and input-file vertical levels
                                                     !! in arbitrary units [A ~> a]
  real, allocatable, dimension(:,:,:), intent(out) :: mask_z
                                                     !< Allocatable tracer mask array on the horizontal
                                                     !! model grid and input-file vertical levels [nondim]
  real, allocatable, dimension(:), intent(out) :: z_in
                                                     !< Cell grid values for input data [Z ~> m]
  real, allocatable, dimension(:), intent(out) :: z_edges_in
                                                     !< Cell grid edge values for input data [Z ~> m]
  real,                  intent(out)   :: missing_value !< The missing value in the returned array, scaled
                                                     !! to avoid accidentally having valid values match
                                                     !! missing values, in the same arbitrary units as tr_z [A ~> a]
  real,                  intent(in)    :: scale      !< Scaling factor for tracer into the internal
                                                     !! units of the model [A a-1 ~> 1]
  logical,     optional, intent(in)    :: homogenize !< If present and true, horizontally homogenize data
                                                     !! to produce perfectly "flat" initial conditions
  logical,     optional, intent(in)    :: spongeOngrid !< If present and true, the sponge data are on the model grid
  real,        optional, intent(in)    :: m_to_Z     !< A conversion factor from meters to the units
                                                     !! of depth [Z m-1 ~> 1].  If missing, G%bathyT must be in m.
  logical,     optional, intent(in)    :: answers_2018 !< If true, use expressions that give the same
                                                     !! answers as the code did in late 2018.  Otherwise
                                                     !! add parentheses for rotational symmetry.
  real,        optional, intent(in)    :: tr_iter_tol !< The tolerance for changes in tracer concentrations
                                                     !! between smoothing iterations that determines when to
                                                     !! stop iterating, in the same arbitrary units as tr_z [A ~> a]
  integer,     optional, intent(in)    :: answer_date !< The vintage of the expressions in the code.
                                                     !! Dates before 20190101 give the same  answers
                                                     !! as the code did in late 2018, while later versions
                                                     !! add parentheses for rotational symmetry.
  type(axis_info), allocatable, dimension(:), optional, intent(inout) :: axes !< Axis types for the input data

  ! Local variables
  ! In the following comments, [A] is used to indicate the arbitrary, possibly rescaled units of the
  ! input array while [a] indicates the unscaled (e.g., mks) units that can be used with the reproducing sums
  real, dimension(:,:),  allocatable   :: tr_in      !< A 2-d array for holding input data on its
                                                     !! native horizontal grid, with units that change
                                                     !! as the input data is interpreted [a] then [A ~> a]
  real, dimension(:,:),  allocatable   :: tr_inp     !< Native horizontal grid data extended to the poles
                                                     !! with units that change as the input data is
                                                     !! interpreted [a] then [A ~> a]
  real, dimension(:,:,:), allocatable  :: data_in    !< A buffer for storing the full 3-d time-interpolated array
                                                     !! on the original grid [a]
  real, dimension(:,:),  allocatable   :: mask_in    !< A 2-d mask for extended input grid [nondim]

  real :: PI_180  ! A conversion factor from degrees to radians [radians degree-1]
  integer :: id, jd, kd, jdp ! Input dataset data sizes
  integer :: i, j, k
  real, dimension(:,:), allocatable :: x_in ! Input file longitudes [radians]
  real, dimension(:,:), allocatable :: y_in ! Input file latitudes [radians]
  real, dimension(:), allocatable :: lon_in ! The longitudes in the input file [degreesE] then [radians]
  real, dimension(:), allocatable :: lat_in ! The latitudes in the input file [degreesN] then [radians]
  real, dimension(:), allocatable :: lat_inp ! The input file latitudes expanded to the pole [degreesN] then [radians]
  real :: max_lat   ! The maximum latitude on the input grid [degreesN]
  real :: pole      ! The sum of tracer values at the pole [a]
  real :: max_depth ! The maximum depth of the ocean [Z ~> m]
  real :: npole     ! The number of points contributing to the pole value [nondim]
  real :: missing_val_in ! The missing value in the input field [a]
  real :: roundoff  ! The magnitude of roundoff, usually ~2e-16 [nondim]
  logical :: add_np
  type(horiz_interp_type) :: Interp
  type(axis_info), dimension(4) :: axes_data
  integer :: is, ie, js, je     ! compute domain indices
  integer :: isg, ieg, jsg, jeg ! global extent
  integer :: isd, ied, jsd, jed ! data domain indices
  integer :: id_clock_read
  integer, dimension(4) :: fld_sz
  logical :: debug=.false.
  logical :: is_ongrid
  integer :: ans_date           ! The vintage of the expressions and order of arithmetic to use
  real :: I_scale               ! The inverse of the scale factor for diagnostic output [a A-1 ~> 1]
  real :: dtr_iter_stop         ! The tolerance for changes in tracer concentrations between smoothing
                                ! iterations that determines when to stop iterating [A ~> a]
  real, dimension(SZI_(G),SZJ_(G)) :: lon_out ! The longitude of points on the model grid [radians]
  real, dimension(SZI_(G),SZJ_(G)) :: lat_out ! The latitude of points on the model grid [radians]
  real, dimension(SZI_(G),SZJ_(G)) :: tr_out  ! The tracer on the model grid [A ~> a]
  real, dimension(SZI_(G),SZJ_(G)) :: mask_out ! The mask on the model grid [nondim]
  real, dimension(SZI_(G),SZJ_(G)) :: good    ! Where the data is valid, this is 1 [nondim]
  real, dimension(SZI_(G),SZJ_(G)) :: fill    ! 1 where the data needs to be filled in [nondim]
  real, dimension(SZI_(G),SZJ_(G)) :: tr_outf ! The tracer concentrations after Ice-9 [A ~> a]
  real, dimension(SZI_(G),SZJ_(G)) :: tr_prev ! The tracer concentrations in the layer above [A ~> a]
  real, dimension(SZI_(G),SZJ_(G)) :: good2   ! 1 where the data is valid after Ice-9 [nondim]
  real, dimension(SZI_(G),SZJ_(G)) :: fill2   ! 1 for points that still need to be filled after Ice-9 [nondim]
  integer :: turns
  integer :: verbosity

  turns = G%HI%turns

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg

  id_clock_read = cpu_clock_id('(Initialize tracer from Z) read', grain=CLOCK_LOOP)

  dtr_iter_stop = 1.0e-3*scale
  if (present(tr_iter_tol)) dtr_iter_stop = tr_iter_tol

  I_scale = 1.0 / scale

  PI_180 = atan(1.0)/45.

  ans_date = 20181231
  if (present(answers_2018)) then ; if (.not.answers_2018) ans_date = 20190101 ; endif
  if (present(answer_date)) ans_date = answer_date

  ! Open NetCDF file and if present, extract data and spatial coordinate information
  ! The convention adopted here requires that the data be written in (i,j,k) ordering.

  call cpu_clock_begin(id_clock_read)

  if (present(axes) .and. allocated(axes)) then
    call get_external_field_info(field, size=fld_sz, missing=missing_val_in)
    axes_data = axes
  else
    call get_external_field_info(field, size=fld_sz, axes=axes_data, missing=missing_val_in)
    if (present(axes)) then
      allocate(axes(4))
      axes = axes_data
    endif
  endif
  missing_value = scale*missing_val_in

  verbosity = MOM_get_verbosity()

  id = fld_sz(1) ; jd  = fld_sz(2) ; kd = fld_sz(3)

  is_ongrid = .false.
  if (PRESENT(spongeOngrid)) is_ongrid = spongeOngrid
  if (.not. is_ongrid) then
    allocate(lon_in(id), lat_in(jd))
    call get_axis_info(axes_data(1), ax_data=lon_in)
    call get_axis_info(axes_data(2), ax_data=lat_in)
  endif

  allocate(z_in(kd), z_edges_in(kd+1))

  allocate(tr_z(isd:ied,jsd:jed,kd), source=0.0)
  allocate(mask_z(isd:ied,jsd:jed,kd), source=0.0)

  call get_axis_info(axes_data(3), ax_data=z_in)

  if (present(m_to_Z)) then ; do k=1,kd ; z_in(k) = m_to_Z * z_in(k) ; enddo ; endif

  call cpu_clock_end(id_clock_read)

  if (.not. is_ongrid) then
    max_lat = maxval(lat_in)
    add_np = .false.
    if (max_lat < 90.0) then
      ! Extrapolate the input data to the north pole using the northern-most latitude.
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
    call horizontal_interp_init()
    lon_in = lon_in*PI_180
    lat_in = lat_in*PI_180
    allocate(x_in(id,jdp), y_in(id,jdp))
    call meshgrid(lon_in, lat_in, x_in, y_in)
    lon_out(:,:) = G%geoLonT(:,:)*PI_180
    lat_out(:,:) = G%geoLatT(:,:)*PI_180
    allocate(data_in(id,jd,kd), source=0.0)
    allocate(tr_in(id,jd), source=0.0)
    allocate(tr_inp(id,jdp), source=0.0)
    allocate(mask_in(id,jdp), source=0.0)
  else
    allocate(data_in(isd:ied,jsd:jed,kd))
  endif

  ! Construct level cell boundaries as the mid-point between adjacent centers.
  z_edges_in(1) = 0.0
  do K=2,kd
    z_edges_in(K) = 0.5*(z_in(k-1)+z_in(k))
  enddo
  z_edges_in(kd+1) = 2.0*z_in(kd) - z_in(kd-1)

  max_depth = maxval(G%bathyT(:,:)) + G%Z_ref
  call max_across_PEs(max_depth)

  if (z_edges_in(kd+1) < max_depth) z_edges_in(kd+1) = max_depth

  !  roundoff = 3.0*EPSILON(missing_value)
  roundoff = 1.e-4

  if (.not.is_ongrid) then
    if (is_root_pe()) &
      call time_interp_external(field, Time, data_in, verbose=(verbosity>5), turns=turns)

    ! Loop through each data level and interpolate to model grid.
    ! After interpolating, fill in points which will be needed to define the layers.
    do k=1,kd
      if (is_root_pe()) then
        tr_in(1:id,1:jd) = data_in(1:id,1:jd,k)
        if (add_np) then
          pole = 0.0 ; npole = 0.0
          do i=1,id
            if (abs(tr_in(i,jd)-missing_val_in) > abs(roundoff*missing_val_in)) then
              pole = pole + tr_in(i,jd)
              npole = npole + 1.0
            endif
          enddo
          if (npole > 0) then
            pole = pole / npole
          else
            pole = missing_val_in
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
        if (abs(tr_inp(i,j)-missing_val_in) > abs(roundoff*missing_val_in)) then
          mask_in(i,j) = 1.0
          tr_inp(i,j) = tr_inp(i,j) * scale
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
        call myStats(tr_inp, missing_value, G, k, 'Tracer from file', unscale=I_scale, full_halo=.true.)
      endif

      tr_out(:,:) = 0.0

      call run_horiz_interp(Interp, tr_inp, tr_out(is:ie,js:je), missing_value=missing_value)

      mask_out(:,:) = 1.0
      do j=js,je ; do i=is,ie
        if (abs(tr_out(i,j)-missing_value) < abs(roundoff*missing_value)) mask_out(i,j) = 0.
      enddo ; enddo

      fill(:,:) = 0.0 ; good(:,:) = 0.0

      do j=js,je ; do i=is,ie
        if (mask_out(i,j) < 1.0) then
          tr_out(i,j) = missing_value
        else
          good(i,j) = 1.0
        endif
        if ((G%mask2dT(i,j) == 1.0) .and. (z_edges_in(k) <= G%bathyT(i,j) + G%Z_ref) .and. &
            (mask_out(i,j) < 1.0)) &
          fill(i,j) = 1.0
      enddo ;  enddo
      call pass_var(fill, G%Domain)
      call pass_var(good, G%Domain)

      if (debug) then
        call myStats(tr_out, missing_value, G, k, 'variable from horiz_interp()', unscale=I_scale)
      endif

      ! Horizontally homogenize data to produce perfectly "flat" initial conditions
      if (PRESENT(homogenize)) then ; if (homogenize) then
        call homogenize_field(tr_out, G, tmp_scale=I_scale, weights=mask_out, answer_date=answer_date)
      endif ; endif

      ! tr_out contains input z-space data on the model grid with missing values
      ! now fill in missing values using "ICE-nine" algorithm.
      tr_outf(:,:) = tr_out(:,:)
      if (k==1) tr_prev(:,:) = tr_outf(:,:)
      good2(:,:) = good(:,:)
      fill2(:,:) = fill(:,:)

      call fill_miss_2d(tr_outf, good2, fill2, tr_prev, G, dtr_iter_stop, answer_date=ans_date)

!     if (debug) then
!       call hchksum(tr_outf, 'field from fill_miss_2d ', G%HI, unscale=I_scale)
!       call myStats(tr_outf, missing_value, G, k, 'field from fill_miss_2d()', unscale=I_scale)
!     endif

      tr_z(:,:,k) = tr_outf(:,:) * G%mask2dT(:,:)
      mask_z(:,:,k) = good2(:,:) + fill2(:,:)
      tr_prev(:,:) = tr_z(:,:,k)

      if (debug) then
        call hchksum(tr_prev, 'field after fill ', G%HI, unscale=I_scale)
      endif

    enddo ! kd
  else
    call time_interp_external(field, Time, data_in, verbose=(verbosity>5), turns=turns)
    do k=1,kd
      do j=js,je
        do i=is,ie
          tr_z(i,j,k) = data_in(i,j,k) * scale
          if (ans_date >= 20190101) mask_z(i,j,k) = 1.
          if (abs(tr_z(i,j,k)-missing_value) < abs(roundoff*missing_value)) mask_z(i,j,k) = 0.
        enddo
      enddo
    enddo
  endif

end subroutine horiz_interp_and_extrap_tracer_fms_id

!> Replace all values of a 2-d field with the weighted average over the valid points.
subroutine homogenize_field(field, G, tmp_scale, weights, answer_date, wt_unscale)
  type(ocean_grid_type),            intent(inout) :: G      !< Ocean grid type
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: field  !< The tracer on the model grid in arbitrary units [A ~> a]
  real,                   optional, intent(in)    :: tmp_scale !< A temporary rescaling factor for the
                                                            !! variable that is reversed in the
                                                            !! return value [a A-1 ~> 1]
  real, dimension(SZI_(G),SZJ_(G)), &
                          optional, intent(in)    :: weights !< The weights for the tracer in arbitrary units that
                                                            !! typically differ from those used by field [B ~> b]
  integer,                optional, intent(in)    :: answer_date !< The vintage of the expressions in the code.
                                                            !! Dates before 20230101 use non-reproducing sums
                                                            !! in their averages, while later versions use
                                                            !! reproducing sums for rotational symmetry and
                                                            !! consistency across PE layouts.
  real,                   optional, intent(in)    :: wt_unscale !< A factor that undoes any dimensional scaling
                                                            !! of the weights so that they can be used with
                                                            !! reproducing sums [b B-1 ~> 1]

  ! Local variables
  ! In the following comments, [A] and [B] are used to indicate the arbitrary, possibly rescaled
  ! units of the input field and the weighting array, while [a] and [b] indicate the corresponding
  ! unscaled (e.g., mks) units that can be used with the reproducing sums
  real, dimension(G%isc:G%iec, G%jsc:G%jec) :: field_for_Sums  ! The field times the weights [A B ~> a b]
  real, dimension(G%isc:G%iec, G%jsc:G%jec) :: weight ! A copy of weights, if it is present, or the
                      ! tracer-point grid mask if it weights is absent [B ~> b]
  real :: var_unscale ! The reciprocal of the scaling factor for the field and weights [a b A-1 B-1 ~> 1]
  real :: wt_sum      ! The sum of the weights, in [B ~> b]
  real :: varsum      ! The weighted sum of field being averaged [A B ~> a b]
  real :: varAvg      ! The average of the field [A ~> a]
  logical :: use_repro_sums  ! If true, use reproducing sums.
  integer :: i, j, is, ie, js, je

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  varAvg = 0.0  ! This value will be used if wt_sum is 0.

  use_repro_sums = .false. ; if (present(answer_date)) use_repro_sums = (answer_date >= 20230101)

  if (present(weights)) then
    do j=js,je ; do i=is,ie
      weight(i,j) = weights(i,j)
    enddo ; enddo
  else
    do j=js,je ; do i=is,ie
      weight(i,j) = G%mask2dT(i,j)
    enddo ; enddo
  endif

  if (use_repro_sums) then
    var_unscale = 1.0 ; if (present(tmp_scale)) var_unscale = tmp_scale
    if (present(wt_unscale)) var_unscale = wt_unscale * var_unscale

    do j=js,je ; do i=is,ie
      field_for_Sums(i,j) = field(i,j) * weight(i,j)
    enddo ; enddo

    wt_sum = reproducing_sum(weight, unscale=wt_unscale)
    if (abs(wt_sum) > 0.0) &
      varAvg = reproducing_sum(field_for_Sums, unscale=var_unscale) * (1.0 / wt_sum)

  else  ! Do the averages with order-dependent sums to reproduce older answers.
    wt_sum = 0 ; varsum = 0.
    do j=js,je ; do i=is,ie
      if (weight(i,j) > 0.0) then
        wt_sum = wt_sum + weight(i,j)
        varsum = varsum + field(i,j) * weight(i,j)
      endif
    enddo ; enddo

    ! Note that these averages will not reproduce across PE layouts or grid rotation.
    call sum_across_PEs(wt_sum)
    if (wt_sum > 0.0) then
      call sum_across_PEs(varsum)
      varAvg = varsum / wt_sum
    endif

  endif

  ! This seems like an unlikely case to ever be used, but it is needed to recreate previous behavior.
  if (present(tmp_scale)) then ; if (tmp_scale == 0.0) varAvg = 0.0 ; endif

  field(:,:) = varAvg

end subroutine homogenize_field


!> Create a 2d-mesh of grid coordinates from 1-d arrays.
subroutine meshgrid(x, y, x_T, y_T)
  real, dimension(:),                   intent(in)    :: x  !< input 1-dimensional vector [arbitrary]
  real, dimension(:),                   intent(in)    :: y  !< input 1-dimensional vector [arbitrary]
  real, dimension(size(x,1),size(y,1)), intent(inout) :: x_T !< output 2-dimensional array [arbitrary]
  real, dimension(size(x,1),size(y,1)), intent(inout) :: y_T !< output 2-dimensional array [arbitrary]

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
