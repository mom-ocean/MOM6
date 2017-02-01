module MOM_tracer_initialization_from_Z

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging, only : hchksum
use MOM_coms, only : max_across_PEs, min_across_PEs
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only :  CLOCK_ROUTINE, CLOCK_LOOP
use MOM_domains, only : pass_var, pass_vector, sum_across_PEs, broadcast
use MOM_domains, only : root_PE, To_All, SCALAR_PAIR, CGRID_NE, AGRID
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, read_param, log_param, param_file_type
use MOM_file_parser, only : log_version
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type, isPointInCell
use MOM_io, only : close_file, fieldtype, file_exists
use MOM_io, only : open_file, read_data, read_axis_data, SINGLE_FILE, MULTIPLE
use MOM_io, only : slasher, vardesc, write_field
use MOM_string_functions, only : uppercase
use MOM_time_manager, only : time_type, set_time
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : setVerticalGridAxes
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use MOM_EOS, only : int_specific_vol_dp
use MOM_ALE, only : ALE_remap_scalar
use MOM_regridding, only : regridding_CS
use MOM_remapping, only : remapping_CS, initialize_remapping
use MOM_remapping, only : remapping_core_h
use MOM_verticalGrid,     only : verticalGrid_type

use mpp_domains_mod, only  : mpp_global_field, mpp_get_compute_domain
use mpp_mod, only          : mpp_broadcast,mpp_root_pe,mpp_sync,mpp_sync_self
use mpp_mod, only          : mpp_max
use horiz_interp_mod, only : horiz_interp_new, horiz_interp,horiz_interp_type
use horiz_interp_mod, only : horiz_interp_init, horiz_interp_del

use netcdf

implicit none ; private

#include <MOM_memory.h>

public :: MOM_initialize_tracer_from_Z, horiz_interp_and_extrap_tracer

character(len=40)  :: mod = "MOM_tracer_initialization_from_Z" ! This module's name.

interface fill_boundaries
  module procedure fill_boundaries_real
  module procedure fill_boundaries_int
end interface

real, parameter :: epsln=1.e-10

contains

subroutine MOM_initialize_tracer_from_Z(h, tr, G, GV, PF, src_file, src_var_nam, &
                                src_var_unit_conversion, src_var_record, &
                                homogenize, useALEremapping, remappingScheme, src_var_gridspec )

! Arguments:
!  (in)     h  - Layer thickness, in m.
!  (inout)  tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (in)      G  - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.

  type(ocean_grid_type),                 intent(inout) :: G   !< Ocean grid structure
  type(verticalGrid_type),               intent(in)    :: GV  !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h
  real, dimension(:,:,:),                pointer       :: tr
  type(param_file_type),                 intent(in)    :: PF
  character(len=*),                      intent(in)    :: src_file, src_var_nam
  real, optional,                        intent(in)    :: src_var_unit_conversion
  integer, optional,                     intent(in)    :: src_var_record
  logical, optional,                     intent(in)    :: homogenize, useALEremapping
  character(len=*),  optional,           intent(in)    :: remappingScheme
  character(len=*),  optional,           intent(in)    :: src_var_gridspec ! Not implemented yet.

  real :: land_fill = 0.0
  character(len=200) :: inputdir ! The directory where NetCDF input files are.
  character(len=200) :: mesg
  real               :: convert
  integer            :: recnum
  character(len=10)  :: remapScheme
  logical            :: homog,useALE

! This include declares and sets the variable "version".
#include "version_variable.h"

  character(len=40)  :: mod = "MOM_initialize_tracers_from_Z" ! This module's name.

  integer :: is, ie, js, je, nz ! compute domain indices
  integer :: isc,iec,jsc,jec    ! global compute domain indices
  integer :: isg, ieg, jsg, jeg ! global extent
  integer :: isd, ied, jsd, jed ! data domain indices

  integer :: i, j, k, kd

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: zi
  real, allocatable, dimension(:,:,:), target :: tr_z, mask_z
  real, allocatable, dimension(:), target :: z_edges_in, z_in

  ! Local variables for ALE remapping
  real, dimension(:), allocatable :: h1, h2, hTarget, deltaE, tmpT1d
  real, dimension(:), allocatable :: tmpT1dIn
  real :: zTopOfCell, zBottomOfCell
  type(remapping_CS) :: remapCS ! Remapping parameters and work arrays

  real, dimension(:,:,:), allocatable :: hSrc

  real :: tempAvg, missing_value
  integer :: nPoints, ans
  integer :: id_clock_routine, id_clock_ALE
  logical :: reentrant_x, tripolar_n

  id_clock_routine = cpu_clock_id('(Initialize tracer from Z)', grain=CLOCK_ROUTINE)
  id_clock_ALE = cpu_clock_id('(Initialize tracer from Z) ALE', grain=CLOCK_LOOP)

  call cpu_clock_begin(id_clock_routine)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg

  call callTree_enter(trim(mod)//"(), MOM_state_initialization.F90")

  call mpp_get_compute_domain(G%domain%mpp_domain,isc,iec,jsc,jec)

  call get_param(PF, mod, "Z_INIT_HOMOGENIZE", homog, &
                 "If True, then horizontally homogenize the interpolated \n"//&
                 "initial conditions.", default=.false.)
  call get_param(PF, mod, "Z_INIT_ALE_REMAPPING", useALE, &
                 "If True, then remap straight to model coordinate from file.",&
                 default=.true.)
  call get_param(PF, mod, "Z_INIT_REMAPPING_SCHEME", remapScheme, &
                 "The remapping scheme to use if using Z_INIT_ALE_REMAPPING\n"//&
                 "is True.", default="PLM")

  ! These are model grid properties, but being applied to the data grid for now.
  ! need to revisit this (mjh)
  reentrant_x = .false. ;  call get_param(PF, mod, "REENTRANT_X", reentrant_x,default=.true.)
  tripolar_n = .false. ;  call get_param(PF, mod, "TRIPOLAR_N", tripolar_n, default=.false.)


  if (PRESENT(homogenize)) homog=homogenize
  if (PRESENT(useALEremapping)) useALE=useALEremapping
  if (PRESENT(remappingScheme)) remapScheme=remappingScheme
  recnum=1
  if (PRESENT(src_var_record)) recnum = src_var_record
  convert=1.0
  if (PRESENT(src_var_unit_conversion)) convert = src_var_unit_conversion


  call horiz_interp_and_extrap_tracer(src_file, src_var_nam, convert, recnum, &
       G, tr_z, mask_z, z_in, z_edges_in, missing_value, reentrant_x, tripolar_n, homog)

  kd = size(z_edges_in,1)-1
  call pass_var(tr_z,G%Domain)
  call pass_var(mask_z,G%Domain)

! Done with horizontal interpolation.
! Now remap to model coordinates
  if (useALE) then
    call cpu_clock_begin(id_clock_ALE)
    ! First we reserve a work space for reconstructions of the source data
    allocate( h1(kd) )
    allocate( hSrc(isd:ied,jsd:jed,kd) )
    allocate( tmpT1dIn(kd) )
    call initialize_remapping( remapCS, remapScheme, boundary_extrapolation=.false. ) ! Data for reconstructions
    ! Next we initialize the regridding package so that it knows about the target grid
    allocate( hTarget(nz) )
    allocate( h2(nz) )
    allocate( tmpT1d(nz) )
    allocate( deltaE(nz+1) )

    do j = js, je ; do i = is, ie
      if (G%mask2dT(i,j)>0.) then
        ! Build the source grid
        zTopOfCell = 0. ; zBottomOfCell = 0. ; nPoints = 0
        do k = 1, kd
          if (mask_z(i,j,k) > 0.) then
            zBottomOfCell = -min( z_edges_in(k+1), G%bathyT(i,j) )
            tmpT1dIn(k) = tr_z(i,j,k)
          elseif (k>1) then
            zBottomOfCell = -G%bathyT(i,j)
            tmpT1dIn(k) = tmpT1dIn(k-1)
          else ! This next block should only ever be reached over land
            tmpT1dIn(k) = -99.9
          endif
          h1(k) = zTopOfCell - zBottomOfCell
          if (h1(k)>0.) nPoints = nPoints + 1
          zTopOfCell = zBottomOfCell ! Bottom becomes top for next value of k
        enddo
        h1(kd) = h1(kd) + ( zTopOfCell + G%bathyT(i,j) ) ! In case data is deeper than model
      else
        tr(i,j,:) = 0.
      endif ! mask2dT
      hSrc(i,j,:) = h1(:)
    enddo ; enddo

    call ALE_remap_scalar(remapCS, G, GV, kd, hSrc, tr_z, h, tr, all_cells=.true. )

    deallocate( hSrc )
    deallocate( h1 )
    deallocate( h2 )
    deallocate( hTarget )
    deallocate( tmpT1d )
    deallocate( tmpT1dIn )
    deallocate( deltaE )

    do k=1,nz
      call myStats(tr(:,:,k),missing_value,is,ie,js,je,k,'Tracer from ALE()')
    enddo
    call cpu_clock_end(id_clock_ALE)
  endif ! useALEremapping

! Fill land values
  do k=1,nz ; do j=js,je ; do i=is,ie
    if (tr(i,j,k) == missing_value) then
      tr(i,j,k)=land_fill
    endif
  enddo ; enddo ; enddo


  call callTree_leave(trim(mod)//'()')
  call cpu_clock_end(id_clock_routine)


end subroutine MOM_initialize_tracer_from_Z

subroutine myStats(array, missing, is, ie, js, je, k, mesg)
  real, dimension(:,:), intent(in) :: array
  real, intent(in) :: missing
  integer :: is,ie,js,je,k
  character(len=*) :: mesg
  ! Local variables
  real :: minA, maxA
  integer :: i,j
  logical :: found
  character(len=120) :: lMesg
  minA = 9.E24 ; maxA = -9.E24 ; found = .false.

  do j = js, je
     do i = is, ie
        if (array(i,j) /= array(i,j)) stop 'Nan!'
        if (abs(array(i,j)-missing)>1.e-6*abs(missing)) then
           if (found) then
              minA = min(minA, array(i,j))
              maxA = max(maxA, array(i,j))
           else
              found = .true.
              minA = array(i,j)
              maxA = array(i,j)
           endif
        endif
     enddo
  enddo
  call min_across_PEs(minA)
  call max_across_PEs(maxA)
  if (is_root_pe()) then
     write(lMesg(1:120),'(2(a,es12.4),a,i3,x,a)') &
          'init_from_Z: min=',minA,' max=',maxA,' Level=',k,trim(mesg)
     call MOM_mesg(lMesg,8)
  endif
end subroutine myStats

subroutine fill_miss_2d(aout,good,fill,prev,G,smooth,num_pass,relc,crit,keep_bug,debug)
  !
  !# Use ICE-9 algorithm to populate points (fill=1) with
  !# valid data (good=1). If no information is available,
  !# Then use a previous guess (prev). Optionally (smooth)
  !# blend the filled points to achieve a more desirable result.
  !
  !  (in)        a   : input 2-d array with missing values
  !  (in)     good   : valid data mask for incoming array (1==good data; 0==missing data)
  !  (in)     fill   : same shape array of points which need filling (1==please fill;0==leave it alone)
  !  (in)     prev   : first guess where isolated holes exist,
  !
  use MOM_coms, only : sum_across_PEs

  type(ocean_grid_type), intent(inout)  :: G
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: aout
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: good, fill
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: prev
  logical, intent(in), optional :: smooth
  integer, intent(in), optional :: num_pass
  real, intent(in), optional    :: relc,crit
  logical, intent(in), optional :: keep_bug, debug


  real, dimension(SZI_(G),SZJ_(G)) :: b,r
  real, dimension(SZI_(G),SZJ_(G)) :: fill_pts,good_,good_new

  integer :: i,j,k
  real    :: east,west,north,south,sor
  real    :: ge,gw,gn,gs,ngood
  logical :: do_smooth,siena_bug
  real    :: nfill, nfill_prev
  integer, parameter :: num_pass_default = 10000
  real, parameter :: relc_default = 0.25, crit_default = 1.e-3

  integer :: npass
  integer :: is, ie, js, je, nz
  real    :: relax_coeff, acrit, ares
  logical :: debug_it

  debug_it=.false.
  if (PRESENT(debug)) debug_it=debug

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  npass = num_pass_default
  if (PRESENT(num_pass)) npass = num_pass

  relax_coeff = relc_default
  if (PRESENT(relc)) relax_coeff = relc

  acrit = crit_default
  if (PRESENT(crit)) acrit = crit

  siena_bug=.false.
  if (PRESENT(keep_bug)) siena_bug = keep_bug

  do_smooth=.false.
  if (PRESENT(smooth)) do_smooth=smooth

  fill_pts(:,:)=fill(:,:)

  nfill = sum(fill(is:ie,js:je))
  call sum_across_PEs(nfill)

  nfill_prev = nfill
  good_(:,:)=good(:,:)
  r(:,:)=0.0

  do while (nfill > 0.0)

     call pass_var(good_,G%Domain)
     call pass_var(aout,G%Domain)

     b(:,:)=aout(:,:)
     good_new(:,:)=good_(:,:)

     do j=js,je
        i_loop: do i=is,ie

           if (good_(i,j) .eq. 1.0 .or. fill(i,j) .eq. 0.) cycle i_loop

           ge=good_(i+1,j);gw=good_(i-1,j)
           gn=good_(i,j+1);gs=good_(i,j-1)
           east=0.0;west=0.0;north=0.0;south=0.0
           if (ge.eq.1.0) east=aout(i+1,j)*ge
           if (gw.eq.1.0) west=aout(i-1,j)*gw
           if (gn.eq.1.0) north=aout(i,j+1)*gn
           if (gs.eq.1.0) south=aout(i,j-1)*gs

           ngood = ge+gw+gn+gs
           if (ngood > 0.) then
              b(i,j)=(east+west+north+south)/ngood
              fill_pts(i,j)=0.0
              good_new(i,j)=1.0
           endif
        enddo i_loop
     enddo

     aout(is:ie,js:je)=b(is:ie,js:je)
     good_(is:ie,js:je)=good_new(is:ie,js:je)
     nfill_prev = nfill
     nfill = sum(fill_pts(is:ie,js:je))
     call sum_across_PEs(nfill)

     if (nfill == nfill_prev .and. PRESENT(prev)) then
        do j=js,je
           do i=is,ie
              if (fill_pts(i,j).eq.1.0) then
                 aout(i,j)=prev(i,j)
                 fill_pts(i,j)=0.0
              endif
           enddo
        enddo
     else if (nfill .eq. nfill_prev) then
        print *,&
             'Unable to fill missing points using either data at the same vertical level from a connected basin'//&
             'or using a point from a previous vertical level.  Make sure that the original data has some valid'//&
             'data in all basins.'
        print *,'nfill=',nfill
     endif

     nfill = sum(fill_pts(is:ie,js:je))
     call sum_across_PEs(nfill)

  end do

  if (do_smooth) then
     do k=1,npass
        call pass_var(aout,G%Domain)
        do j=js,je
           do i=is,ie
              if (fill(i,j) .eq. 1) then
                 east=max(good(i+1,j),fill(i+1,j));west=max(good(i-1,j),fill(i-1,j))
                 north=max(good(i,j+1),fill(i,j+1));south=max(good(i,j-1),fill(i,j-1))
                 r(i,j) = relax_coeff*(south*aout(i,j-1)+north*aout(i,j+1)+west*aout(i-1,j)+east*aout(i+1,j) - (south+north+west+east)*aout(i,j))
              else
                 r(i,j) = 0.
              endif
           enddo
        enddo
        aout(is:ie,js:je)=r(is:ie,js:je)+aout(is:ie,js:je)
        ares = maxval(abs(r))
        call max_across_PEs(ares)
        if (ares <= acrit) exit
     enddo
  endif

  do j=js,je
     do i=is,ie
        if (good_(i,j).eq.0.0 .and. fill_pts(i,j) .eq. 1.0) then
           print *,'in fill_miss, fill, good,i,j= ',fill_pts(i,j),good_(i,j),i,j
           call MOM_error(FATAL,"MOM_initialize: "// &
                "fill is true and good is false after fill_miss, how did this happen? ")
        endif
     enddo
  enddo

  return

end subroutine fill_miss_2d

subroutine horiz_interp_and_extrap_tracer(filename, varnam,  conversion, recnum, G, tr_z, mask_z, z_in, z_edges_in, missing_value, reentrant_x, tripolar_n, homogenize )

  character(len=*), intent(in) :: filename   ! Path to file containing tracer to be interpolated
  character(len=*), intent(in) :: varnam     ! name of tracer in filee
  real,             intent(in) :: conversion ! conversion factor for tracer
  integer,          intent(in) :: recnum     ! record number of tracer to be read
  type(ocean_grid_type), intent(inout) :: G     ! Grid object
  real, allocatable, dimension(:,:,:) :: tr_z    ! pointer to allocatable tracer array on local model grid
                                              ! and native vertical levels
  real, allocatable, dimension(:,:,:) :: mask_z   ! pointer to allocatable tracer mask array on local model grid
                                              ! and native vertical levels
  real, allocatable, dimension(:)     :: z_in  ! Cell grid values for input data
  real, allocatable, dimension(:)     :: z_edges_in  ! Cell grid edge values for input data
  real,                intent(out)  :: missing_value
  logical,          intent(in) :: reentrant_x, tripolar_n
  logical, intent(in), optional     :: homogenize

  real, dimension(:,:), allocatable :: tr_in,tr_inp ! A 2-d array for holding input data on native horizontal
                                                    ! grid and extended grid with poles
  real, dimension(:,:), allocatable :: mask_in      ! A 2-d mask for extended input grid

  real :: PI_180
  integer :: rcode, ncid, varid, ndims, id, jd, kd, jdp
  integer :: i,j,k
  integer, dimension(4) :: start, count, dims, dim_id
  real, dimension(:,:), allocatable :: x_in, y_in
  real, dimension(:), allocatable  :: lon_in, lat_in
  real, dimension(:), allocatable  :: lat_inp, last_row
  real :: max_lat, min_lat, pole, max_depth, npole
  real :: roundoff  ! The magnitude of roundoff, usually ~2e-16.
  logical :: add_np
  character(len=8)  :: laynum
  type(horiz_interp_type) :: Interp
  integer :: is, ie, js, je     ! compute domain indices
  integer :: isc,iec,jsc,jec    ! global compute domain indices
  integer :: isg, ieg, jsg, jeg ! global extent
  integer :: isd, ied, jsd, jed ! data domain indices
  integer :: ni, nj, nz         ! global grid size
  integer :: id_clock_read
  character(len=12)  :: dim_name(4)
  logical :: debug=.false.
  real :: npoints,varAvg
  real, dimension(SZI_(G),SZJ_(G)) :: lon_out, lat_out, tr_out, mask_out
  real, dimension(SZI_(G),SZJ_(G)) :: good, fill
  real, dimension(SZI_(G),SZJ_(G)) :: tr_outf,tr_prev
  real, dimension(SZI_(G),SZJ_(G))  :: good2,fill2
  real, dimension(SZI_(G),SZJ_(G))  :: nlevs

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg

  id_clock_read = cpu_clock_id('(Initialize tracer from Z) read', grain=CLOCK_LOOP)


  if (ALLOCATED(tr_z)) deallocate(tr_z)
  if (ALLOCATED(mask_z)) deallocate(mask_z)
  if (ALLOCATED(z_edges_in)) deallocate(z_edges_in)

  PI_180=atan(1.0)/45.

  ! Open NetCDF file and if present, extract data and spatial coordinate information
  ! The convention adopted here requires that the data be written in (i,j,k) ordering.

  call cpu_clock_begin(id_clock_read)


  rcode = NF90_OPEN(filename, NF90_NOWRITE, ncid)
  if (rcode .ne. 0) call MOM_error(FATAL,"error opening file "//trim(filename)//&
                           " in hinterp_extrap")
  rcode = NF90_INQ_VARID(ncid, varnam, varid)
  if (rcode .ne. 0) call MOM_error(FATAL,"error finding variable "//trim(varnam)//&
                                 " in file "//trim(filename)//" in hinterp_extrap")

  rcode = NF90_INQUIRE_VARIABLE(ncid, varid, ndims=ndims, dimids=dims)
  if (rcode .ne. 0) call MOM_error(FATAL,'error inquiring dimensions hinterp_extrap')
  if (ndims < 3) call MOM_error(FATAL,"Variable "//trim(varnam)//" in file "// &
              trim(filename)//" has too few dimensions.")

  rcode = NF90_INQUIRE_DIMENSION(ncid, dims(1), dim_name(1), len=id)
  if (rcode .ne. 0) call MOM_error(FATAL,"error reading dimension 1 data for "// &
                trim(varnam)//" in file "// trim(filename)//" in hinterp_extrap")
  rcode = NF90_INQ_VARID(ncid, dim_name(1), dim_id(1))
  if (rcode .ne. 0) call MOM_error(FATAL,"error finding variable "//trim(dim_name(1))//&
                                 " in file "//trim(filename)//" in hinterp_extrap")
  rcode = NF90_INQUIRE_DIMENSION(ncid, dims(2), dim_name(2), len=jd)
  if (rcode .ne. 0) call MOM_error(FATAL,"error reading dimension 2 data for "// &
                trim(varnam)//" in file "// trim(filename)//" in hinterp_extrap")
  rcode = NF90_INQ_VARID(ncid, dim_name(2), dim_id(2))
  if (rcode .ne. 0) call MOM_error(FATAL,"error finding variable "//trim(dim_name(2))//&
                                 " in file "//trim(filename)//" in hinterp_extrap")
  rcode = NF90_INQUIRE_DIMENSION(ncid, dims(3), dim_name(3), len=kd)
  if (rcode .ne. 0) call MOM_error(FATAL,"error reading dimension 3 data for "// &
                trim(varnam)//" in file "// trim(filename)//" in hinterp_extrap")
  rcode = NF90_INQ_VARID(ncid, dim_name(3), dim_id(3))
  if (rcode .ne. 0) call MOM_error(FATAL,"error finding variable "//trim(dim_name(3))//&
                                 " in file "//trim(filename)//" in hinterp_extrap")


  missing_value=0.0
  rcode = NF90_GET_ATT(ncid, varid, "_FillValue", missing_value)
  if (rcode .ne. 0) call MOM_error(FATAL,"error finding missing value for "//&
       trim(varnam)//" in file "// trim(filename)//" in hinterp_extrap")

  if (allocated(lon_in)) deallocate(lon_in)
  if (allocated(lat_in)) deallocate(lat_in)
  if (allocated(z_in)) deallocate(z_in)
  if (allocated(z_edges_in)) deallocate(z_edges_in)
  if (allocated(tr_z)) deallocate(tr_z)
  if (allocated(mask_z)) deallocate(mask_z)

  allocate(lon_in(id),lat_in(jd),z_in(kd),z_edges_in(kd+1))
  allocate(tr_z(isd:ied,jsd:jed,kd), mask_z(isd:ied,jsd:jed,kd))

  start = 1; count = 1; count(1) = id
  rcode = NF90_GET_VAR(ncid, dim_id(1), lon_in, start, count)
  if (rcode .ne. 0) call MOM_error(FATAL,"error reading dimension 1 values for var_name "// &
                trim(varnam)//",dim_name "//trim(dim_name(1))//" in file "// trim(filename)//" in hinterp_extrap")
  start = 1; count = 1; count(1) = jd
  rcode = NF90_GET_VAR(ncid, dim_id(2), lat_in, start, count)
  if (rcode .ne. 0) call MOM_error(FATAL,"error reading dimension 2 values for var_name "// &
                trim(varnam)//",dim_name "//trim(dim_name(2))//" in file "// trim(filename)//" in  hinterp_extrap")
  start = 1; count = 1; count(1) = kd
  rcode = NF90_GET_VAR(ncid, dim_id(3), z_in, start, count)
  if (rcode .ne. 0) call MOM_error(FATAL,"error reading dimension 3 values for var_name "// &
                trim(varnam//",dim_name "//trim(dim_name(3)))//" in file "// trim(filename)//" in  hinterp_extrap")

  call cpu_clock_end(id_clock_read)

! extrapolate the input data to the north pole using the northerm-most latitude

  max_lat = maxval(lat_in)
  add_np=.false.
  if (max_lat < 90.0) then
    add_np=.true.
    jdp=jd+1
    allocate(lat_inp(jdp))
    lat_inp(1:jd)=lat_in(:)
    lat_inp(jd+1)=90.0
    deallocate(lat_in)
    allocate(lat_in(1:jdp))
    lat_in(:)=lat_inp(:)
  else
    jdp=jd
  endif

! construct level cell boundaries as the mid-point between adjacent centers

  z_edges_in(1) = 0.0
  do k=2,kd
   z_edges_in(k)=0.5*(z_in(k-1)+z_in(k))
  enddo
  z_edges_in(kd+1)=2.0*z_in(kd) - z_in(kd-1)

  call horiz_interp_init()

  lon_in = lon_in*PI_180
  lat_in = lat_in*PI_180
  allocate(x_in(id,jdp),y_in(id,jdp))
  call meshgrid(lon_in,lat_in, x_in, y_in)

  lon_out(:,:) = G%geoLonT(:,:)*PI_180
  lat_out(:,:) = G%geoLatT(:,:)*PI_180


  allocate(tr_in(id,jd)) ; tr_in(:,:)=0.0
  allocate(tr_inp(id,jdp)) ; tr_inp(:,:)=0.0
  allocate(mask_in(id,jdp)) ; mask_in(:,:)=0.0
  allocate(last_row(id))    ; last_row(:)=0.0

  max_depth = maxval(G%bathyT)
  call mpp_max(max_depth)

  if (z_edges_in(kd+1)<max_depth) z_edges_in(kd+1)=max_depth


! loop through each data level and interpolate to model grid.
! after interpolating, fill in points which will be needed
! to define the layers

  roundoff = 3.0*EPSILON(missing_value)

  do k=1,kd
    write(laynum,'(I8)') k ; laynum = adjustl(laynum)

    if (is_root_pe()) then
      start = 1; start(3) = k; count = 1; count(1) = id; count(2) = jd
      rcode = NF90_GET_VAR(ncid,varid, tr_in, start, count)
      if (rcode .ne. 0) call MOM_error(FATAL,"hinterp_and_extract_from_Fie: "//&
           "error reading level "//trim(laynum)//" of variable "//&
           trim(varnam)//" in file "// trim(filename))

      if (add_np) then
         last_row(:)=tr_in(:,jd); pole=0.0;npole=0.0
         do i=1,id
            if (abs(tr_in(i,jd)-missing_value) .gt. abs(roundoff*missing_value)) then
               pole = pole+last_row(i)
               npole = npole+1.0
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

    call mpp_sync()
    call mpp_broadcast(tr_inp,id*jdp,root_PE())
    call mpp_sync_self ()

    mask_in=0.0

    do j=1,jdp
      do i=1,id
         if (abs(tr_inp(i,j)-missing_value) .gt. abs(roundoff*missing_value)) then
           mask_in(i,j)=1.0
            tr_inp(i,j) = tr_inp(i,j) * conversion
         else
           tr_inp(i,j)=missing_value
         endif
      enddo
    enddo


! call fms routine horiz_interp to interpolate input level data to model horizontal grid


    if (k == 1) then
      call horiz_interp_new(Interp,x_in,y_in,lon_out(is:ie,js:je),lat_out(is:ie,js:je), &
               interp_method='bilinear',src_modulo=reentrant_x)
    endif

    if (debug) then
       call myStats(tr_inp,missing_value, is,ie,js,je,k,'Tracer from file')
    endif

    tr_out(:,:) = 0.0

    call horiz_interp(Interp,tr_inp,tr_out(is:ie,js:je), missing_value=missing_value, new_missing_handle=.true.)

    mask_out=1.0
    do j=js,je
      do i=is,ie
        if (abs(tr_out(i,j)-missing_value) .lt. abs(roundoff*missing_value)) mask_out(i,j)=0.
      enddo
    enddo

    fill = 0.0; good = 0.0

    nPoints = 0 ; varAvg = 0.
    do j=js,je
      do i=is,ie
        if (mask_out(i,j) .lt. 1.0) then
          tr_out(i,j)=missing_value
        else
          good(i,j)=1.0
          nPoints = nPoints + 1
          varAvg = varAvg + tr_out(i,j)
        endif
        if (G%mask2dT(i,j) == 1.0 .and. z_edges_in(k) <= G%bathyT(i,j) .and. mask_out(i,j) .lt. 1.0) fill(i,j)=1.0
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

    tr_outf(:,:)=tr_out(:,:)
    if (k==1) tr_prev(:,:)=tr_outf(:,:)
    good2(:,:)=good(:,:)
    fill2(:,:)=fill(:,:)

    call fill_miss_2d(tr_outf,good2,fill2,tr_prev,G,smooth=.true.)
    call myStats(tr_outf,missing_value,is,ie,js,je,k,'field from fill_miss_2d()')

    tr_z(:,:,k) = tr_outf(:,:)*G%mask2dT(:,:)
    mask_z(:,:,k) = good2(:,:)+fill2(:,:)

    tr_prev(:,:)=tr_z(:,:,k)

    if (debug) then
      call hchksum(tr_prev,'field after fill ',G%HI)
    endif

  enddo ! kd

end subroutine horiz_interp_and_extrap_tracer

function tracer_z_init(tr_in,z_edges,e,nkml,nkbl,land_fill,wet,nlay,nlevs,debug,i_debug,j_debug) result(tr)
!
! Adopted from R. Hallberg
! Arguments:
!  (in)     tr_in  - The z-space array of tracer concentrations that is read in.
!  (in)   z_edges  - The depths of the cell edges in the input z* data (m)
!  (in)         e  - The depths of the layer interfaces (m)
!  (in)       nkml - number of mixed layer pieces
!  (in)       nkbl - number of buffer layer pieces
!  (in)  land_fill - fill in data over land
!  (in)        wet - wet mask (1=ocean)
!  (in)       nlay - number of layers
!  (in)      nlevs - number of levels

!  (out)        tr - tracers on layers

!    tr_1d    ! A copy of the input tracer concentrations in a column.
!    wt       ! The fractional weight for each layer in the range between
              ! k_top and k_bot, nondim.
!    z1       ! z1 and z2 are the depths of the top and bottom limits of the part
!    z2       ! of a z-cell that contributes to a layer, relative to the cell
!               center and normalized by the cell thickness, nondim.
!               Note that -1/2 <= z1 <= z2 <= 1/2.
!
real, dimension(:,:,:), intent(in)           :: tr_in
real, dimension(size(tr_in,3)+1), intent(in) :: z_edges
integer, intent(in)                          :: nlay
real, dimension(size(tr_in,1),size(tr_in,2),nlay+1), intent(in) :: e
integer, intent(in)                          :: nkml,nkbl
real, intent(in)                             :: land_fill
real, dimension(size(tr_in,1),size(tr_in,2)), intent(in) :: wet
real, dimension(size(tr_in,1),size(tr_in,2)), optional, intent(in) ::nlevs
logical, intent(in), optional                :: debug
integer, intent(in), optional                :: i_debug, j_debug

real, dimension(size(tr_in,1),size(tr_in,2),nlay) :: tr
real, dimension(size(tr_in,3)) :: tr_1d
real, dimension(nlay+1) :: e_1d
real, dimension(nlay) :: tr_
integer, dimension(size(tr_in,1),size(tr_in,2)) :: nlevs_data

integer :: n,i,j,k,l,nx,ny,nz,nt,kz
integer :: k_top,k_bot,k_bot_prev,kk,kstart
real    :: sl_tr
real, dimension(size(tr_in,3)) :: wt,z1,z2
logical :: debug_msg, debug_

nx = size(tr_in,1); ny=size(tr_in,2); nz = size(tr_in,3)

nlevs_data = size(tr_in,3)
if (PRESENT(nlevs)) then
  nlevs_data  = anint(nlevs)
endif

debug_=.false.
if (PRESENT(debug)) then
  debug_=debug
endif

debug_msg = .false.
if (debug_) then
  debug_msg=.true.
endif


do j=1,ny
  i_loop: do i=1,nx
    if (nlevs_data(i,j) .eq. 0 .or. wet(i,j) .eq. 0.) then
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
      elseif (e_1d(k) < z_edges(nlevs_data(i,j)+1)) then
        if (debug_msg) then
          print *,'*** WARNING : Found interface below valid range of z data '
          print *,'(i,j,z_bottom,interface)= ',&
               i,j,z_edges(nlevs_data(i,j)+1),e_1d(k)
          print *,'z_edges= ',z_edges
          print *,'e=',e_1d
          print *,'*** I will extrapolate below using the bottom-most valid values'
          debug_msg = .false.
        endif
        tr(i,j,k) = tr_1d(nlevs_data(i,j))

      else
        kstart=k_bot
        call find_overlap(z_edges, e_1d(k), e_1d(k+1), nlevs_data(i,j), &
             kstart, k_top, k_bot, wt, z1, z2)

        if (debug_) then
           if (PRESENT(i_debug)) then
             if (i.eq.i_debug.and.j.eq.j_debug) then
                print *,'0001 k,k_top,k_bot,sum(wt),sum(z2-z1) = ',k,k_top,k_bot,sum(wt),sum(z2-z1)
             endif
           endif
        endif
        kz = k_top
        sl_tr=0.0; ! cur_tr=0.0
        if (kz /= k_bot_prev) then
! Calculate the intra-cell profile.
          if ((kz < nlevs_data(i,j)) .and. (kz > 1)) then
            sl_tr = find_limited_slope(tr_1d, z_edges, kz)
          endif
        endif
        if (kz > nlevs_data(i,j)) kz = nlevs_data(i,j)
! This is the piecewise linear form.
        tr(i,j,k) = wt(kz) * &
             (tr_1d(kz) + 0.5*sl_tr*(z2(kz) + z1(kz)))
! For the piecewise parabolic form add the following...
!     + C1_3*cur_tr*(z2(kz)**2 + z2(kz)*z1(kz) + z1(kz)**2))
!        if (debug_) then
!          print *,'k,k_top,k_bot= ',k,k_top,k_bot
!        endif
        if (debug_) then
           if (PRESENT(i_debug)) then
             if (i.eq.i_debug.and.j.eq.j_debug) then
                print *,'0002 k,k_top,k_bot,k_bot_prev,sl_tr = ',k,k_top,k_bot,k_bot_prev,sl_tr
             endif
           endif
        endif

        do kz=k_top+1,k_bot-1
          tr(i,j,k) = tr(i,j,k) + wt(kz)*tr_1d(kz)
        enddo

        if (debug_) then
           if (PRESENT(i_debug)) then
             if (i.eq.i_debug.and.j.eq.j_debug) then
                print *,'0003 k,tr = ',k,tr(i,j,k)
             endif
           endif
        endif

        if (k_bot > k_top) then
          kz = k_bot
! Calculate the intra-cell profile.
          sl_tr = 0.0 ! ; cur_tr = 0.0
          if ((kz < nlevs_data(i,j)) .and. (kz > 1)) then
            sl_tr = find_limited_slope(tr_1d, z_edges, kz)
!            if (debug_) then
!              print *,'002 sl_tr,k,kz,nlevs= ',sl_tr,k,kz,nlevs_data(i,j),nlevs(i,j)
!            endif
          endif
! This is the piecewise linear form.
          tr(i,j,k) = tr(i,j,k) + wt(kz) * &
               (tr_1d(kz) + 0.5*sl_tr*(z2(kz) + z1(kz)))
! For the piecewise parabolic form add the following...
!     + C1_3*cur_tr*(z2(kz)**2 + z2(kz)*z1(kz) + z1(kz)**2))

          if (debug_) then
             if (PRESENT(i_debug)) then
               if (i.eq.i_debug.and.j.eq.j_debug) then
                  print *,'0004 k,kz,nlevs,sl_tr,tr = ',k,kz,nlevs_data(i,j),sl_tr,tr(i,j,k)
                  print *,'0005 k,kz,tr(kz),tr(kz-1),tr(kz+1) = ',k,kz,tr_1d(kz),tr_1d(kz-1),tr_1d(kz+1),z_edges(kz+2)
               endif
             endif
          endif

        endif
        k_bot_prev = k_bot

      endif
    enddo ! k-loop

    do k=2,nlay  ! simply fill vanished layers with adjacent value
      if (e_1d(k)-e_1d(k+1) .le. epsln) tr(i,j,k)=tr(i,j,k-1)
    enddo

  enddo i_loop
enddo

return

end function tracer_z_init


function bisect_fast(a, x, lo, hi) result(bi_r)
!
!  Return the index where to insert item x in list a, assuming a is sorted.
!  The return values [i] is such that all e in a[:i-1] have e <= x, and all e in
!  a[i:] have e > x. So if x already appears in the list, will
!  insert just after the rightmost x already there.
!  Optional args lo (default 1) and hi (default len(a)) bound the
!  slice of a to be searched.
!
!  (in)  a - sorted list
!  (in)  x - item to be inserted
!  (in)  lo, hi - optional range to search

real, dimension(:,:), intent(in) :: a
real, dimension(:), intent(in) :: x
integer, dimension(size(a,1)), intent(in), optional  :: lo,hi
integer, dimension(size(a,1),size(x,1))  :: bi_r

integer :: mid,num_x,num_a,i
integer, dimension(size(a,1))  :: lo_,hi_,lo0,hi0
integer :: nprofs,j

lo_=1;hi_=size(a,2);num_x=size(x,1);bi_r=-1;nprofs=size(a,1)

if (PRESENT(lo)) then
  where (lo>0) lo_=lo
end if
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


#ifdef PY_SOLO
subroutine determine_temperature(temp,salt,R,p_ref,niter,land_fill,h,k_start)

! # This subroutine determines the potential temperature and
! # salinity that is consistent with the target density
! # using provided initial guess
! #   (inout)     temp - potential temperature (degC)
! #   (inout)     salt - salinity (PSU)
! #   (in)           R - Desired potential density, in kg m-3.
! #   (in)       p_ref - Reference pressure, in Pa.
! #   (in)       niter - maximum number of iterations
! #   (in)           h - layer thickness . Do not iterate for massless layers
! #   (in)     k_start - starting index (i.e. below the buffer layer)
! #   (in)   land_fill - land fill value

real(kind=8), dimension(:,:,:), intent(inout) :: temp,salt
real(kind=8), dimension(size(temp,3)), intent(in) :: R
real, intent(in) :: p_ref
integer, intent(in) :: niter
integer, intent(in) :: k_start
real, intent(in) :: land_fill
real(kind=8), dimension(:,:,:), intent(in) :: h

real(kind=8), dimension(size(temp,1),size(temp,3)) :: T,S,dT,dS,rho,hin
real(kind=8), dimension(size(temp,1),size(temp,3)) :: drho_dT,drho_dS
real(kind=8), dimension(size(temp,1)) :: press

integer :: nx,ny,nz,nt,i,j,k,n,itt
logical :: adjust_salt , old_fit
real    :: dT_dS
real, parameter :: T_max = 35.0, T_min = -2.0
real, parameter :: S_min = 0.5, S_max=65.0
real, parameter :: tol=1.e-4, max_t_adj=1.0, max_s_adj = 0.5

#else

subroutine determine_temperature(temp,salt,R,p_ref,niter,land_fill,h,k_start,eos)

! # This subroutine determines the potential temperature and
! # salinity that is consistent with the target density
! # using provided initial guess
! #   (inout)     temp - potential temperature (degC)
! #   (inout)     salt - salinity (PSU)
! #   (in)           R - Desired potential density, in kg m-3.
! #   (in)       p_ref - Reference pressure, in Pa.
! #   (in)       niter - maximum number of iterations
! #   (in)           h - layer thickness . Do not iterate for massless layers
! #   (in)     k_start - starting index (i.e. below the buffer layer)
! #   (in)   land_fill - land fill value
! #   (in)        eos  - seawater equation of state

real, dimension(:,:,:), intent(inout) :: temp,salt
real, dimension(size(temp,3)), intent(in) :: R
real, intent(in) :: p_ref
integer, intent(in) :: niter
integer, intent(in) :: k_start
real, intent(in) :: land_fill
real, dimension(:,:,:), intent(in) :: h
type(eos_type), pointer, intent(in) :: eos

real(kind=8), dimension(size(temp,1),size(temp,3)) :: T,S,dT,dS,rho,hin
real(kind=8), dimension(size(temp,1),size(temp,3)) :: drho_dT,drho_dS
real(kind=8), dimension(size(temp,1)) :: press

integer :: nx,ny,nz,nt,i,j,k,n,itt
real    :: dT_dS
logical :: adjust_salt , old_fit
real, parameter :: T_max = 31.0, T_min = -2.0
real, parameter :: S_min = 0.5, S_max=65.0
real, parameter :: tol=1.e-4, max_t_adj=1.0, max_s_adj = 0.5


#endif


old_fit = .true.   ! reproduces siena behavior
                   ! will switch to the newer
                   ! method which simultaneously adjusts
                   ! temp and salt based on the ratio
                   ! of the thermal and haline coefficients.

nx=size(temp,1);ny=size(temp,2); nz=size(temp,3)

press(:) = p_ref

do j=1,ny
  dS(:,:) = 0. ! Needs to be zero everywhere since there is a maxval(abs(dS)) later...
  T=temp(:,j,:)
  S=salt(:,j,:)
  hin=h(:,j,:)
  dT=0.0
  adjust_salt = .true.
  iter_loop: do itt = 1,niter
#ifdef PY_SOLO
    rho=wright_eos_2d(T,S,p_ref)
    drho_dT=alpha_wright_eos_2d(T,S,p_ref)
#else
    do k=1, nz
      call calculate_density(T(:,k),S(:,k),press,rho(:,k),1,nx,eos)
      call calculate_density_derivs(T(:,k),S(:,k),press,drho_dT(:,k),drho_dS(:,k),1,nx,eos)
    enddo
#endif
    do k=k_start,nz
      do i=1,nx

!               if (abs(rho(i,k)-R(k))>tol .and. hin(i,k)>epsln .and. abs(T(i,k)-land_fill) < epsln) then
        if (abs(rho(i,k)-R(k))>tol) then
           if (old_fit) then
              dT(i,k)=(R(k)-rho(i,k))/drho_dT(i,k)
              if (dT(i,k)>max_t_adj) dT(i,k)=max_t_adj
              if (dT(i,k)<-1.0*max_t_adj) dT(i,k)=-1.0*max_t_adj
              T(i,k)=max(min(T(i,k)+dT(i,k),T_max),T_min)
           else
              dT_dS = 10.0 - min(-drho_dT(i,k)/drho_dS(i,k),10.)
              dS(i,k) = (R(k)-rho(i,k))/(drho_dS(i,k) - drho_dT(i,k)*dT_dS )
              dT(i,k)= -dT_dS*dS(i,k)
              !          if (dT(i,k)>max_t_adj) dT(i,k)=max_t_adj
              !          if (dT(i,k)<-1.0*max_t_adj) dT(i,k)=-1.0*max_t_adj
              T(i,k)=max(min(T(i,k)+dT(i,k),T_max),T_min)
              S(i,k)=max(min(S(i,k)+dS(i,k),S_max),S_min)
           endif
        endif
      enddo
    enddo
    if (maxval(abs(dT)) < tol) then
       adjust_salt = .false.
      exit iter_loop
    endif
  enddo iter_loop

  if (adjust_salt .and. old_fit) then
    iter_loop2: do itt = 1,niter
#ifdef PY_SOLO
      rho=wright_eos_2d(T,S,p_ref)
      drho_dS=beta_wright_eos_2d(T,S,p_ref)
#else
      do k=1, nz
        call calculate_density(T(:,k),S(:,k),press,rho(:,k),1,nx,eos)
        call calculate_density_derivs(T(:,k),S(:,k),press,drho_dT(:,k),drho_dS(:,k),1,nx,eos)
      enddo
#endif
      do k=k_start,nz
        do i=1,nx
!                   if (abs(rho(i,k)-R(k))>tol .and. hin(i,k)>epsln .and. abs(T(i,k)-land_fill) < epsln ) then
          if (abs(rho(i,k)-R(k))>tol ) then
             dS(i,k)=(R(k)-rho(i,k))/drho_dS(i,k)
             if (dS(i,k)>max_s_adj) dS(i,k)=max_s_adj
             if (dS(i,k)<-1.0*max_s_adj) dS(i,k)=-1.0*max_s_adj
             S(i,k)=max(min(S(i,k)+dS(i,k),S_max),S_min)
          endif
        enddo
      enddo
      if (maxval(abs(dS)) < tol) then
        exit iter_loop2
      endif
    enddo iter_loop2
  endif

  temp(:,j,:)=T(:,:)
  salt(:,j,:)=S(:,:)
enddo


return

end subroutine determine_temperature


subroutine find_overlap(e, Z_top, Z_bot, k_max, k_start, k_top, k_bot, wt, z1, z2)

!   This subroutine determines the layers bounded by interfaces e that overlap
! with the depth range between Z_top and Z_bot, and also the fractional weights
! of each layer. It also calculates the normalized relative depths of the range
! of each layer that overlaps that depth range.
!   Note that by convention, e decreases with increasing k and Z_top > Z_bot.
!
! Arguments: e - A column's interface heights, in m.
!  (in)      Z_top - The top of the range being mapped to, in m.
!  (in)      Z_bot - The bottom of the range being mapped to, in m.
!  (in)      k_max - The number of valid layers.
!  (in)      k_start - The layer at which to start searching.
!  (out)     k_top, k_bot - The indices of the top and bottom layers that
!                           overlap with the depth range.
!  (out)     wt - The relative weights of each layer from k_top to k_bot.
!  (out)     z1, z2 - z1 and z2 are the depths of the top and bottom limits of
!                     the part of a layer that contributes to a depth level,
!                     relative to the cell center and normalized by the cell
!                     thickness, nondim.  Note that -1/2 <= z1 < z2 <= 1/2.

real, dimension(:), intent(in) :: e
real, intent(in)   :: Z_top, Z_bot
integer, intent(in) :: k_max, k_start
integer, intent(out) :: k_top, k_bot
real, dimension(:), intent(out) :: wt, z1, z2

real :: Ih, e_c, tot_wt, I_totwt
integer :: k

wt(:)=0.0; z1(:)=0.0; z2(:)=0.0
k_top = k_start; k_bot= k_start; wt(1) = 1.0; z1(1)=-0.5; z2(1) = 0.5

do k=k_start,k_max ;if (e(k+1)<Z_top) exit ; enddo
k_top = k
if (k>k_max) return

! Determine the fractional weights of each layer.
! Note that by convention, e and Z_int decrease with increasing k.
if (e(k+1)<=Z_bot) then
  wt(k) = 1.0 ; k_bot = k
  Ih = 1.0 / (e(k)-e(k+1))
  e_c = 0.5*(e(k)+e(k+1))
  z1(k) = (e_c - MIN(e(k),Z_top)) * Ih
  z2(k) = (e_c - Z_bot) * Ih
else
  wt(k) = MIN(e(k),Z_top) - e(k+1) ; tot_wt = wt(k) ! These are always > 0.
  z1(k) = (0.5*(e(k)+e(k+1)) - MIN(e(k),Z_top)) / (e(k)-e(k+1))
  z2(k) = 0.5
  k_bot = k_max
  do k=k_top+1,k_max
    if (e(k+1)<=Z_bot) then
      k_bot = k
      wt(k) = e(k) - Z_bot ; z1(k) = -0.5
      z2(k) = (0.5*(e(k)+e(k+1)) - Z_bot) / (e(k)-e(k+1))
    else
      wt(k) = e(k) - e(k+1) ; z1(k) = -0.5 ; z2(k) = 0.5
    endif
    tot_wt = tot_wt + wt(k) ! wt(k) is always > 0.
    if (k>=k_bot) exit
  enddo

  I_totwt = 1.0 / tot_wt
  do k=k_top,k_bot ; wt(k) = I_totwt*wt(k) ; enddo
endif

return

end subroutine find_overlap


function find_limited_slope(val, e, k) result(slope)

!   This subroutine determines a limited slope for val to be advected with
! a piecewise limited scheme.

! Arguments: val - An column the values that are being interpolated.
!  (in)      e - A column's interface heights, in m.
!  (in)      slope - The normalized slope in the intracell distribution of val.
!  (in)      k - The layer whose slope is being determined.


real, dimension(:), intent(in) :: val
real, dimension(:), intent(in) :: e
integer, intent(in) :: k
real :: slope,amx,bmx,amn,bmn,cmn,dmn

real :: d1, d2

if ((val(k)-val(k-1)) * (val(k)-val(k+1)) >= 0.0) then
  slope = 0.0 ! ; curvature = 0.0
else
  d1 = 0.5*(e(k-1)-e(k+1)) ; d2 = 0.5*(e(k)-e(k+2))
  slope = ((d1**2)*(val(k+1) - val(k)) + (d2**2)*(val(k) - val(k-1))) * &
       (e(k) - e(k+1)) / (d1*d2*(d1+d2))
! slope = 0.5*(val(k+1) - val(k-1))
! This is S.J. Lin's form of the PLM limiter.
  amx=max(val(k-1),val(k))
  bmx = max(amx,val(k+1))
  amn = min(abs(slope),2.0*(bmx-val(k)))
  bmn = min(val(k-1),val(k))
  cmn = 2.0*(val(k)-min(bmn,val(k+1)))
  dmn = min(amn,cmn)
  slope = sign(1.0,slope) * dmn

! min(abs(slope), &
!             2.0*(max(val(k-1),val(k),val(k+1)) - val(k)), &
!             2.0*(val(k) - min(val(k-1),val(k),val(k+1))))
! curvature = 0.0
endif

return

end function find_limited_slope



function find_interfaces(rho,zin,Rb,depth,nlevs,nkml,nkbl,hml,debug) result(zi)
!  (in)      rho : potential density in z-space (kg m-3)
!  (in)      zin : levels (m)
!  (in)      Rb  : target interface densities (kg m-3)
!  (in)     depth: ocean depth (m)
!  (in)     nlevs: number of valid points in each column
!  (in)     nkml : number of mixed layer pieces
!  (in)     nkbl : number of buffer layer pieces
!  (in)      hml : mixed layer depth

real, dimension(:,:,:), intent(in) :: rho
real, dimension(size(rho,3)), intent(in) :: zin
real, dimension(:), intent(in) :: Rb
real, dimension(size(rho,1),size(rho,2)), intent(in) :: depth
real, dimension(size(rho,1),size(rho,2)), optional, intent(in) ::nlevs
logical, optional, intent(in) :: debug
real, dimension(size(rho,1),size(rho,2),size(Rb,1)) :: zi
integer, intent(in), optional :: nkml, nkbl
real, intent(in), optional    :: hml

real, dimension(size(rho,1),size(rho,3)) :: rho_
real, dimension(size(rho,1)) :: depth_
logical :: unstable
integer :: dir
integer, dimension(size(rho,1),size(Rb,1)) :: ki_
real, dimension(size(rho,1),size(Rb,1)) :: zi_
integer, dimension(size(rho,1),size(rho,2)) :: nlevs_data
integer, dimension(size(rho,1)) :: lo,hi
real :: slope,rsm,drhodz,hml_
integer :: n,i,j,k,l,nx,ny,nz,nt
integer :: nlay,kk,nkml_,nkbl_
logical :: debug_ = .false.

real, parameter :: zoff=0.999

nlay=size(Rb)-1

zi=0.0


if (PRESENT(debug)) debug_=debug

nx = size(rho,1); ny=size(rho,2); nz = size(rho,3)
nlevs_data(:,:) = size(rho,3)

nkml_=0;nkbl_=0;hml_=0.0
if (PRESENT(nkml)) nkml_=max(0,nkml)
if (PRESENT(nkbl)) nkbl_=max(0,nkbl)
if (PRESENT(hml)) hml_=hml

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
            if (k.eq.2) then
              rho_(i,k-1)=rho_(i,k)-epsln
            else
              drhodz = (rho_(i,k+1)-rho_(i,k-1))/(zin(k+1)-zin(k-1))
              if (drhodz  < 0.0) then
                unstable=.true.
              endif
              rho_(i,k) = rho_(i,k-1)+drhodz*zoff*(zin(k)-zin(k-1))
            endif
          endif
        enddo
        dir=-1*dir
      else
        do k=nlevs_data(i,j)-1,2,-1
          if (rho_(i,k+1) - rho_(i,k) < 0.0) then
            if (k .eq. nlevs_data(i,j)-1) then
              rho_(i,k+1)=rho_(i,k-1)+epsln
            else
              drhodz = (rho_(i,k+1)-rho_(i,k-1))/(zin(k+1)-zin(k-1))
              if (drhodz  < 0.0) then
                unstable=.true.
              endif
              rho_(i,k) = rho_(i,k+1)-drhodz*(zin(k+1)-zin(k))
            endif
          endif
        enddo
        dir=-1*dir
      endif
    enddo
    if (debug_) then
      print *,'final density profile= ', rho_(i,:)
    endif
  enddo i_loop

  ki_(:,:) = 0
  zi_(:,:) = 0.0
  depth_(:)=-1.0*depth(:,j)
  lo(:)=1
  hi(:)=nlevs_data(:,j)
  ki_ = bisect_fast(rho_,Rb,lo,hi)
  ki_(:,:) = max(1,ki_(:,:)-1)
  do i=1,nx
    do l=2,nlay
      slope = (zin(ki_(i,l)+1) - zin(ki_(i,l)))/max(rho_(i,ki_(i,l)+1) - rho_(i,ki_(i,l)),epsln)
      zi_(i,l) = -1.0*(zin(ki_(i,l)) + slope*(Rb(l)-rho_(i,ki_(i,l))))
      zi_(i,l) = max(zi_(i,l),depth_(i))
      zi_(i,l) = min(zi_(i,l),-1.0*hml_)
    enddo
    zi_(i,nlay+1)=depth_(i)
    do l=2,nkml_+1
      zi_(i,l)=max(((1.0-real(l))/real(nkml_))*hml_,depth_(i))
    enddo
    do l=nlay,nkml_+2,-1
      if (zi_(i,l) < zi_(i,l+1)+epsln) then
        zi_(i,l)=zi_(i,l+1)+epsln
      endif
      if (zi_(i,l)>-1.0*hml_) then
        zi_(i,l)=max(-1.0*hml_,depth_(i))
      endif
    enddo
  enddo
  zi(:,j,:)=zi_(:,:)
enddo

return


end function find_interfaces

subroutine meshgrid(x,y,x_T,y_T)

!  create a 2d-mesh of grid coordinates
!  from 1-d arrays.

real, dimension(:), intent(in) :: x,y
real, dimension(size(x,1),size(y,1)), intent(inout) :: x_T,y_T

integer :: ni,nj,i,j

ni=size(x,1);nj=size(y,1)

do j=1,nj
  x_T(:,j)=x(:)
enddo

do i=1,ni
  y_T(i,:)=y(:)
enddo

return

end subroutine meshgrid

subroutine smooth_heights(zi,fill,bad,sor,niter,cyclic_x, tripolar_n)
!
! Solve del2 (zi) = 0 using successive iterations
! with a 5 point stencil. Only points fill==1 are
! modified. Except where bad==1, information propagates
! isotropically in index space.  The resulting solution
! in each region is an approximation to del2(zi)=0 subject to
! boundary conditions along the valid points curve bounding this region.


real, dimension(:,:), intent(inout) :: zi
integer, dimension(size(zi,1),size(zi,2)), intent(in) :: fill
integer, dimension(size(zi,1),size(zi,2)), intent(in) :: bad
real, intent(in)  :: sor
integer, intent(in) :: niter
logical, intent(in) :: cyclic_x, tripolar_n

integer :: i,j,k,n
integer :: ni,nj

real, dimension(size(zi,1),size(zi,2)) :: res, m
integer, dimension(size(zi,1),size(zi,2),4) :: B
real, dimension(0:size(zi,1)+1,0:size(zi,2)+1) :: mp
integer, dimension(0:size(zi,1)+1,0:size(zi,2)+1) :: nm

real :: Isum, bsum

ni=size(zi,1); nj=size(zi,2)


mp=fill_boundaries(zi,cyclic_x,tripolar_n)

B(:,:,:)=0.0
nm=fill_boundaries(bad,cyclic_x,tripolar_n)

do j=1,nj
  do i=1,ni
    if (fill(i,j) .eq. 1) then
      B(i,j,1)=1-nm(i+1,j);B(i,j,2)=1-nm(i-1,j)
      B(i,j,3)=1-nm(i,j+1);B(i,j,4)=1-nm(i,j-1)
    endif
  enddo
enddo

do n=1,niter
  do j=1,nj
    do i=1,ni
      if (fill(i,j) .eq. 1) then
        bsum = real(B(i,j,1)+B(i,j,2)+B(i,j,3)+B(i,j,4))
        Isum = 1.0/bsum
        res(i,j)=Isum*(B(i,j,1)*mp(i+1,j)+B(i,j,2)*mp(i-1,j)+&
             B(i,j,3)*mp(i,j+1)+B(i,j,4)*mp(i,j-1)) - mp(i,j)
      endif
    enddo
  enddo
  res(:,:)=res(:,:)*sor

  do j=1,nj
    do i=1,ni
      mp(i,j)=mp(i,j)+res(i,j)
    enddo
  enddo

  zi(:,:)=mp(1:ni,1:nj)
  mp = fill_boundaries(zi,cyclic_x,tripolar_n)
end do



return

end subroutine smooth_heights

function fill_boundaries_int(m,cyclic_x,tripolar_n) result(mp)
!
! fill grid edges
!
integer, dimension(:,:), intent(in) :: m
logical, intent(in) :: cyclic_x, tripolar_n
real, dimension(size(m,1),size(m,2)) :: m_real
real, dimension(0:size(m,1)+1,0:size(m,2)+1) :: mp_real
integer, dimension(0:size(m,1)+1,0:size(m,2)+1) :: mp

m_real = real(m)

mp_real = fill_boundaries_real(m_real,cyclic_x,tripolar_n)

mp = int(mp_real)

return

end function fill_boundaries_int

function fill_boundaries_real(m,cyclic_x,tripolar_n) result(mp)
!
! fill grid edges
!
real, dimension(:,:), intent(in) :: m
logical, intent(in) :: cyclic_x, tripolar_n
real, dimension(0:size(m,1)+1,0:size(m,2)+1) :: mp

integer :: ni,nj,i,j

ni=size(m,1); nj=size(m,2)

mp(1:ni,1:nj)=m(:,:)

if (cyclic_x) then
  mp(0,1:nj)=m(ni,1:nj)
  mp(ni+1,1:nj)=m(1,1:nj)
else
  mp(0,1:nj)=m(1,1:nj)
  mp(ni+1,1:nj)=m(ni,1:nj)
endif

mp(1:ni,0)=m(1:ni,1)
if (tripolar_n) then
  do i=1,ni
    mp(i,nj+1)=m(ni-i+1,nj)
  enddo
else
  mp(1:ni,nj+1)=m(1:ni,nj)
endif

return

end function fill_boundaries_real



end module MOM_tracer_initialization_from_Z
