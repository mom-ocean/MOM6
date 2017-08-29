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
use MOM_horizontal_regridding, only : myStats, horiz_interp_and_extrap_tracer
implicit none ; private

#include <MOM_memory.h>

public :: MOM_initialize_tracer_from_Z

character(len=40)  :: mdl = "MOM_tracer_initialization_from_Z" ! This module's name.

real, parameter :: epsln=1.e-10

contains

subroutine MOM_initialize_tracer_from_Z(h, tr, G, GV, PF, src_file, src_var_nam, &
                                src_var_unit_conversion, src_var_record, &
                                homogenize, useALEremapping, remappingScheme, src_var_gridspec )

! Arguments:
!  (in)     h  - Layer thickness, in m.
!  (inout)  tr - pointer to array containing field to be initialized.
!  (in)      G  - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.

  type(ocean_grid_type),       intent(inout) :: G   !< Ocean grid structure.
  type(verticalGrid_type),     intent(in)    :: GV  !< Ocean vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                               intent(in)    :: h   !< Layer thickness, in m.
  real, dimension(:,:,:), pointer, intent(inout) :: tr  !< Pointer to array to be initialized
  type(param_file_type),       intent(in)    :: PF  !< parameter file
  character(len=*),            intent(in)    :: src_file, src_var_nam !< source filename and variable name on disk
  real, optional,              intent(in)    :: src_var_unit_conversion !< optional multiplicative unit conversion
  integer, optional,           intent(in)    :: src_var_record  !< record to read for multiple time-level files
  logical, optional,           intent(in)    :: homogenize !< optionally homogenize to mean value
  logical, optional,           intent(in)    :: useALEremapping !< to remap or not (optional)
  character(len=*),  optional, intent(in)    :: remappingScheme !< remapping scheme to use.
  character(len=*),  optional, intent(in)    :: src_var_gridspec ! Not implemented yet.

  real :: land_fill = 0.0
  character(len=200) :: inputdir ! The directory where NetCDF input files are.
  character(len=200) :: mesg
  real               :: convert
  integer            :: recnum
  character(len=10)  :: remapScheme
  logical            :: homog,useALE

! This include declares and sets the variable "version".
#include "version_variable.h"

  character(len=40)  :: mdl = "MOM_initialize_tracers_from_Z" ! This module's name.

  integer :: is, ie, js, je, nz ! compute domain indices
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

  call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")

  call get_param(PF, mdl, "Z_INIT_HOMOGENIZE", homog, &
                 "If True, then horizontally homogenize the interpolated \n"//&
                 "initial conditions.", default=.false.)
  call get_param(PF, mdl, "Z_INIT_ALE_REMAPPING", useALE, &
                 "If True, then remap straight to model coordinate from file.",&
                 default=.true.)
  call get_param(PF, mdl, "Z_INIT_REMAPPING_SCHEME", remapScheme, &
                 "The remapping scheme to use if using Z_INIT_ALE_REMAPPING\n"//&
                 "is True.", default="PLM")

  ! These are model grid properties, but being applied to the data grid for now.
  ! need to revisit this (mjh)
  reentrant_x = .false. ;  call get_param(PF, mdl, "REENTRANT_X", reentrant_x,default=.true.)
  tripolar_n = .false. ;  call get_param(PF, mdl, "TRIPOLAR_N", tripolar_n, default=.false.)


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

    call ALE_remap_scalar(remapCS, G, GV, kd, hSrc, tr_z, h, tr, all_cells=.false. )

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


  call callTree_leave(trim(mdl)//'()')
  call cpu_clock_end(id_clock_routine)


end subroutine MOM_initialize_tracer_from_Z






end module MOM_tracer_initialization_from_Z
