module regrid_defs
!==============================================================================
!
! This file is part of MOM.
!
! Date of creation: 2008.12.15
! L. White
!
! This module contains the parameters and types used for the
! regridding/remapping.
!
!==============================================================================
implicit none ; public

! -----------------------------------------------------------------------------
! Definition of type 'regridding_opts' that contains regridding parameters
! -----------------------------------------------------------------------------
type, public :: regridding_opts_t

  ! Indicates whether regridding/remapping capabilities should be used at all
  ! (it should be true to use regridding and false to NOT use regridding).
  logical   :: use_regridding;

  ! Indicates which grid to use in the vertical (z*, sigma, target interface
  ! densities)
  integer   :: regridding_scheme;

  ! The following parameter is only relevant when used with the target
  ! interface densities regridding scheme. It indicates which interpolation
  ! to use to determine the grid.
  integer   :: interpolation_scheme;

  ! Indicates which remapping scheme to use to remap variables between grids
  integer   :: remapping_scheme;

  ! Indicate whether high-order boundary extrapolation should be used within
  ! boundary cells
  logical   :: boundary_extrapolation;

  ! Minimum thickness allowed when building the new grid through regridding
  real      :: min_thickness;

  ! Indicates whether integrals for FV pressure gradient calculation will
  ! use reconstruction of T/S.
  ! By default, it is true when use_regridding=True
  logical   :: reconstructForPressure

  ! The form of the reconstruction of T/S for FV pressure gradient calculation.
  ! By default, it is =1 (PLM)
  integer   :: pressureReconstructionScheme

end type regridding_opts_t

! -----------------------------------------------------------------------------
! Reconstruction schemes for integrals in FV pressure gradient calculation 
! -----------------------------------------------------------------------------

! List of reconstruction schemes
integer, parameter  :: PRESSURE_RECONSTRUCTION_PLM   = 1
integer, parameter  :: PRESSURE_RECONSTRUCTION_PPM   = 2

end module regrid_defs
