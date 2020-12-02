!> Initializes hydrography from z-coordinate climatology files
module MOM_tracer_initialization_from_Z

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging, only : hchksum
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_ROUTINE, CLOCK_LOOP
use MOM_domains, only : pass_var
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, param_file_type, log_version
use MOM_grid, only : ocean_grid_type
use MOM_horizontal_regridding, only : myStats, horiz_interp_and_extrap_tracer
use MOM_remapping, only : remapping_CS, initialize_remapping
use MOM_unit_scaling, only : unit_scale_type
use MOM_verticalGrid, only : verticalGrid_type
use MOM_ALE, only : ALE_remap_scalar

implicit none ; private

#include <MOM_memory.h>

public :: MOM_initialize_tracer_from_Z

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

character(len=40)  :: mdl = "MOM_tracer_initialization_from_Z" !< This module's name.

contains

!> Initializes a tracer from a z-space data file, including any lateral regridding that is needed.
subroutine MOM_initialize_tracer_from_Z(h, tr, G, GV, US, PF, src_file, src_var_nam, &
                          src_var_unit_conversion, src_var_record, homogenize, &
                          useALEremapping, remappingScheme, src_var_gridspec )
  type(ocean_grid_type),      intent(inout) :: G   !< Ocean grid structure.
  type(verticalGrid_type),    intent(in)    :: GV  !< Ocean vertical grid structure.
  type(unit_scale_type),      intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in)    :: h   !< Layer thickness [H ~> m or kg m-2].
  real, dimension(:,:,:),     pointer       :: tr  !< Pointer to array to be initialized
  type(param_file_type),      intent(in)    :: PF  !< parameter file
  character(len=*),           intent(in)    :: src_file !< source filename
  character(len=*),           intent(in)    :: src_var_nam !< variable name in file
  real,             optional, intent(in)    :: src_var_unit_conversion !< optional multiplicative unit conversion
  integer,          optional, intent(in)    :: src_var_record  !< record to read for multiple time-level files
  logical,          optional, intent(in)    :: homogenize !< optionally homogenize to mean value
  logical,          optional, intent(in)    :: useALEremapping !< to remap or not (optional)
  character(len=*), optional, intent(in)    :: remappingScheme !< remapping scheme to use.
  character(len=*), optional, intent(in)    :: src_var_gridspec !< Source variable name in a gridspec file.
                                                                !! This is not implemented yet.
  ! Local variables
  real :: land_fill = 0.0
  character(len=200) :: inputdir ! The directory where NetCDF input files are.
  character(len=200) :: mesg
  real               :: convert
  integer            :: recnum
  character(len=10)  :: remapScheme
  logical            :: homog,useALE

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_initialize_tracers_from_Z" ! This module's name.

  integer :: is, ie, js, je, nz ! compute domain indices
  integer :: isd, ied, jsd, jed ! data domain indices
  integer :: i, j, k, kd
  real, allocatable, dimension(:,:,:), target :: tr_z, mask_z
  real, allocatable, dimension(:), target :: z_edges_in, z_in

  ! Local variables for ALE remapping
  real, dimension(:,:,:), allocatable :: hSrc ! Source thicknesses [H ~> m or kg m-2].
  real, dimension(:), allocatable :: h1 ! A 1-d column of source thicknesses [Z ~> m].
  real :: zTopOfCell, zBottomOfCell, z_bathy  ! Heights [Z ~> m].
  type(remapping_CS) :: remapCS ! Remapping parameters and work arrays

  real :: missing_value
  integer :: nPoints
  integer :: id_clock_routine, id_clock_ALE
  logical :: answers_2018, default_2018_answers, hor_regrid_answers_2018
  logical :: reentrant_x, tripolar_n

  id_clock_routine = cpu_clock_id('(Initialize tracer from Z)', grain=CLOCK_ROUTINE)
  id_clock_ALE = cpu_clock_id('(Initialize tracer from Z) ALE', grain=CLOCK_LOOP)

  call cpu_clock_begin(id_clock_routine)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call callTree_enter(trim(mdl)//"(), MOM_state_initialization.F90")

  call get_param(PF, mdl, "Z_INIT_HOMOGENIZE", homog, &
                 "If True, then horizontally homogenize the interpolated "//&
                 "initial conditions.", default=.false.)
  call get_param(PF, mdl, "Z_INIT_ALE_REMAPPING", useALE, &
                 "If True, then remap straight to model coordinate from file.",&
                 default=.true.)
  call get_param(PF, mdl, "Z_INIT_REMAPPING_SCHEME", remapScheme, &
                 "The remapping scheme to use if using Z_INIT_ALE_REMAPPING is True.", &
                 default="PLM")
  call get_param(PF, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false.)
  if (useALE) then
    call get_param(PF, mdl, "REMAPPING_2018_ANSWERS", answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)
  endif
  call get_param(PF, mdl, "HOR_REGRID_2018_ANSWERS", hor_regrid_answers_2018, &
                 "If true, use the order of arithmetic for horizonal regridding that recovers "//&
                 "the answers from the end of 2018.  Otherwise, use rotationally symmetric "//&
                 "forms of the same expressions.", default=default_2018_answers)

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
       G, tr_z, mask_z, z_in, z_edges_in, missing_value, reentrant_x, tripolar_n, &
       homog, m_to_Z=US%m_to_Z, answers_2018=hor_regrid_answers_2018)

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
    ! Set parameters for reconstructions
    call initialize_remapping( remapCS, remapScheme, boundary_extrapolation=.false., answers_2018=answers_2018 )
    ! Next we initialize the regridding package so that it knows about the target grid

    do j = js, je ; do i = is, ie
      if (G%mask2dT(i,j)>0.) then
        ! Build the source grid
        zTopOfCell = 0. ; zBottomOfCell = 0. ; nPoints = 0
        z_bathy = G%bathyT(i,j)
        do k = 1, kd
          if (mask_z(i,j,k) > 0.) then
            zBottomOfCell = -min( z_edges_in(k+1), z_bathy )
          elseif (k>1) then
            zBottomOfCell = -z_bathy
          endif
          h1(k) = zTopOfCell - zBottomOfCell
          if (h1(k)>0.) nPoints = nPoints + 1
          zTopOfCell = zBottomOfCell ! Bottom becomes top for next value of k
        enddo
        h1(kd) = h1(kd) + ( zTopOfCell + z_bathy ) ! In case data is deeper than model
      else
        tr(i,j,:) = 0.
      endif ! mask2dT
      hSrc(i,j,:) = GV%Z_to_H * h1(:)
    enddo ; enddo

    call ALE_remap_scalar(remapCS, G, GV, kd, hSrc, tr_z, h, tr, all_cells=.false., answers_2018=answers_2018 )

    deallocate( hSrc )
    deallocate( h1 )

    do k=1,nz
      call myStats(tr(:,:,k), missing_value, is, ie, js, je, k, 'Tracer from ALE()')
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
