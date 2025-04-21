!> Initializes hydrography from z-coordinate climatology files
module MOM_tracer_initialization_from_Z

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging,     only : hchksum
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,     only : CLOCK_ROUTINE, CLOCK_LOOP
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,   only : get_param, param_file_type, log_version
use MOM_grid,          only : ocean_grid_type
use MOM_horizontal_regridding, only : myStats, horiz_interp_and_extrap_tracer
use MOM_interface_heights, only : dz_to_thickness_simple
use MOM_regridding,    only : set_dz_neglect, set_h_neglect
use MOM_remapping,     only : remapping_CS, initialize_remapping
use MOM_unit_scaling,  only : unit_scale_type
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_ALE,           only : ALE_remap_scalar

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
                          useALEremapping, remappingScheme, src_var_gridspec, h_in_Z_units, &
                          ongrid)
  type(ocean_grid_type),      intent(inout) :: G   !< Ocean grid structure.
  type(verticalGrid_type),    intent(in)    :: GV  !< Ocean vertical grid structure.
  type(unit_scale_type),      intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in)    :: h   !< Layer thicknesses, in [H ~> m or kg m-2] or
                                                   !! [Z ~> m] depending on the value of h_in_Z_units.
  real, dimension(:,:,:),     pointer       :: tr  !< Pointer to array to be initialized [CU ~> conc]
  type(param_file_type),      intent(in)    :: PF  !< parameter file
  character(len=*),           intent(in)    :: src_file !< source filename
  character(len=*),           intent(in)    :: src_var_nam !< variable name in file
  real,             optional, intent(in)    :: src_var_unit_conversion !< optional multiplicative unit conversion,
                                                   !! often used for rescaling into model units [CU conc-1 ~> 1]
  integer,          optional, intent(in)    :: src_var_record  !< record to read for multiple time-level files
  logical,          optional, intent(in)    :: homogenize !< optionally homogenize to mean value
  logical,          optional, intent(in)    :: useALEremapping !< to remap or not (optional)
  character(len=*), optional, intent(in)    :: remappingScheme !< remapping scheme to use.
  character(len=*), optional, intent(in)    :: src_var_gridspec !< Source variable name in a gridspec file.
                                                                !! This is not implemented yet.
  logical,          optional, intent(in)    :: h_in_Z_units !< If present and true, the input grid
                                                            !! thicknesses are in the units of height
                                                            !! ([Z ~> m]) instead of the usual units of
                                                            !! thicknesses ([H ~> m or kg m-2])
  logical,          optional, intent(in)    :: ongrid       !< If true, then data are assumed to have been
                                                            !! interpolated to the model horizontal grid. In this case,
                                                            !! only extrapolation is performed by
                                                            !! horiz_interp_and_extrap_tracer()
  ! Local variables
  real :: land_fill = 0.0  ! A value to use to replace missing values [CU ~> conc]
  real :: convert ! A conversion factor into the model's internal units [CU conc-1 ~> 1]
  integer            :: recnum
  character(len=64)  :: remapScheme
  logical            :: homog, useALE
  logical            :: h_is_in_Z_units

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_initialize_tracers_from_Z" ! This module's name.

  integer :: is, ie, js, je, nz ! compute domain indices
  integer :: isd, ied, jsd, jed ! data domain indices
  integer :: i, j, k, kd
  real, allocatable, dimension(:,:,:), target :: tr_z   ! Tracer array on the horizontal model grid
                                                        ! and input-file vertical levels [CU ~> conc]
  real, allocatable, dimension(:,:,:), target :: mask_z ! Missing value mask on the horizontal model grid
                                                        ! and input-file vertical levels [nondim]
  real, allocatable, dimension(:), target :: z_edges_in ! Cell edge depths for input data [Z ~> m]
  real, allocatable, dimension(:), target :: z_in       ! Cell center depths for input data [Z ~> m]

  ! Local variables for ALE remapping
  real, dimension(:,:,:), allocatable :: dzSrc ! Source thicknesses in height units [Z ~> m]
  real, dimension(:,:,:), allocatable :: hSrc  ! Source thicknesses [H ~> m or kg m-2]
  real, dimension(:), allocatable :: h1 ! A 1-d column of source thicknesses [Z ~> m].
  real :: zTopOfCell, zBottomOfCell, z_bathy  ! Heights [Z ~> m].
  type(remapping_CS) :: remapCS ! Remapping parameters and work arrays
  type(verticalGrid_type) :: GV_loc ! A temporary vertical grid structure

  real :: missing_value ! A value indicating that there is no valid input data at this point [CU ~> conc]
  real :: dz_neglect              ! A negligibly small vertical layer extent used in
                                  ! remapping cell reconstructions [Z ~> m] or [H ~> m or kg m-2]
  real :: dz_neglect_edge         ! A negligibly small vertical layer extent used in
                                  ! remapping edge value calculations [Z ~> m] or [H ~> m or kg m-2]
  logical :: om4_remap_via_sub_cells ! If true, use the OM4 remapping algorithm
  integer :: nPoints    ! The number of valid input data points in a column
  integer :: id_clock_routine, id_clock_ALE
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  integer :: remap_answer_date    ! The vintage of the order of arithmetic and expressions to use
                                  ! for remapping.  Values below 20190101 recover the remapping
                                  ! answers from 2018, while higher values use more robust
                                  ! forms of the same remapping expressions.
  integer :: hor_regrid_answer_date  ! The vintage of the order of arithmetic and expressions to use
                                  ! for horizontal regridding.  Values below 20190101 recover the
                                  ! answers from 2018, while higher values use expressions that have
                                  ! been rearranged for rotational invariance.

  id_clock_routine = cpu_clock_id('(Initialize tracer from Z)', grain=CLOCK_ROUTINE)
  id_clock_ALE = cpu_clock_id('(Initialize tracer from Z) ALE', grain=CLOCK_LOOP)

  call cpu_clock_begin(id_clock_routine)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call callTree_enter(trim(mdl)//"(), MOM_tracer_initialization_from_Z.F90")

  call get_param(PF, mdl, "Z_INIT_HOMOGENIZE", homog, &
                 "If True, then horizontally homogenize the interpolated "//&
                 "initial conditions.", default=.false.)
  call get_param(PF, mdl, "Z_INIT_ALE_REMAPPING", useALE, &
                 "If True, then remap straight to model coordinate from file.",&
                 default=.false.)
  call get_param(PF, mdl, "Z_INIT_REMAPPING_SCHEME", remapScheme, &
                 "The remapping scheme to use if using Z_INIT_ALE_REMAPPING is True.", &
                 default="PPM_IH4")
  call get_param(PF, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
  if (useALE) then
    call get_param(PF, mdl, "REMAPPING_ANSWER_DATE", remap_answer_date, &
                 "The vintage of the expressions and order of arithmetic to use for remapping.  "//&
                 "Values below 20190101 result in the use of older, less accurate expressions "//&
                 "that were in use at the end of 2018.  Higher values result in the use of more "//&
                 "robust and accurate forms of mathematically equivalent expressions.", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
    call get_param(PF, mdl, "REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                   do_not_log=.true., default=.true.)
    call get_param(PF, mdl, "Z_INIT_REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                 "If true, use the OM4 remapping-via-subcells algorithm for initialization. "//&
                 "See REMAPPING_USE_OM4_SUBCELLS for more details. "//&
                 "We recommend setting this option to false.", default=om4_remap_via_sub_cells)
    if (.not.GV%Boussinesq) remap_answer_date = max(remap_answer_date, 20230701)
  endif
  call get_param(PF, mdl, "HOR_REGRID_ANSWER_DATE", hor_regrid_answer_date, &
                 "The vintage of the order of arithmetic for horizontal regridding.  "//&
                 "Dates before 20190101 give the same answers as the code did in late 2018, "//&
                 "while later versions add parentheses for rotational symmetry.  "//&
                 "Dates after 20230101 use reproducing sums for global averages.", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
  if (.not.GV%Boussinesq) hor_regrid_answer_date = max(hor_regrid_answer_date, 20230701)

  if (PRESENT(homogenize)) homog=homogenize
  if (PRESENT(useALEremapping)) useALE=useALEremapping
  if (PRESENT(remappingScheme)) remapScheme=remappingScheme
  recnum = 1
  if (PRESENT(src_var_record)) recnum = src_var_record
  convert = 1.0
  if (PRESENT(src_var_unit_conversion)) convert = src_var_unit_conversion

  h_is_in_Z_units = .false. ; if (present(h_in_Z_units)) h_is_in_Z_units = h_in_Z_units

  call horiz_interp_and_extrap_tracer(src_file, src_var_nam, recnum, &
            G, tr_z, mask_z, z_in, z_edges_in, missing_value, &
            scale=convert, homogenize=homog, m_to_Z=US%m_to_Z, &
            answer_date=hor_regrid_answer_date, ongrid=ongrid)

  kd = size(z_edges_in,1)-1
  call pass_var(tr_z,G%Domain)
  call pass_var(mask_z,G%Domain)

! Done with horizontal interpolation.
! Now remap to model coordinates
  if (useALE) then
    call cpu_clock_begin(id_clock_ALE)
    ! First we reserve a work space for reconstructions of the source data
    allocate( h1(kd) )
    allocate( dzSrc(isd:ied,jsd:jed,kd) )
    allocate( hSrc(isd:ied,jsd:jed,kd) )
    ! Set parameters for reconstructions in the right units
    if (h_is_in_Z_units) then
      dz_neglect = set_dz_neglect(GV, US, remap_answer_date, dz_neglect_edge)
    else
      dz_neglect = set_h_neglect(GV, remap_answer_date, dz_neglect_edge)
    endif
    call initialize_remapping( remapCS, remapScheme, boundary_extrapolation=.false., &
                               om4_remap_via_sub_cells=om4_remap_via_sub_cells, answer_date=remap_answer_date, &
                               H_neglect=dz_neglect, H_neglect_edge=dz_neglect_edge )
    ! Next we initialize the regridding package so that it knows about the target grid

    do j = js, je ; do i = is, ie
      if (G%mask2dT(i,j)>0.) then
        ! Build the source grid
        zTopOfCell = 0. ; zBottomOfCell = 0. ; nPoints = 0
        z_bathy = G%bathyT(i,j) + G%Z_ref
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
      dzSrc(i,j,:) = h1(:)
    enddo ; enddo

    if (h_is_in_Z_units) then
      ! Because h is in units of [Z ~> m], dzSrc is already in the right units.
      call ALE_remap_scalar(remapCS, G, GV, kd, dzSrc, tr_z, h, tr, all_cells=.false.)
    else
      ! Equation of state data is not available, so a simpler rescaling will have to suffice,
      ! but it might be problematic in non-Boussinesq mode.
      GV_loc = GV ; GV_loc%ke = kd
      call dz_to_thickness_simple(dzSrc, hSrc, G, GV_loc, US)

      call ALE_remap_scalar(remapCS, G, GV, kd, hSrc, tr_z, h, tr, all_cells=.false.)
    endif

    deallocate( hSrc )
    deallocate( dzSrc )
    deallocate( h1 )

    do k=1,nz
      call myStats(tr(:,:,k), missing_value, G, k, 'Tracer from ALE()')
    enddo
    call cpu_clock_end(id_clock_ALE)
  endif ! useALEremapping

! Fill land values
  do k=1,nz ; do j=js,je ; do i=is,ie
    if (tr(i,j,k) == missing_value) then
      tr(i,j,k) = land_fill
    endif
  enddo ; enddo ; enddo

  call callTree_leave(trim(mdl)//'()')
  call cpu_clock_end(id_clock_routine)

end subroutine MOM_initialize_tracer_from_Z

end module MOM_tracer_initialization_from_Z
