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

  ! Maximum domain depth
  real      :: max_depth;

  ! Minimum thickness allowed when building the new grid through regridding
  real      :: min_thickness;

end type regridding_opts_t

! -----------------------------------------------------------------------------
! Regridding parameters 
! -----------------------------------------------------------------------------
! Maximum number of regridding iterations
integer, parameter  :: NB_REGRIDDING_ITERATIONS = 1;
! Deviation tolerance between succesive grids in regridding iterations
real, parameter     :: DEVIATION_TOLERANCE = 1e-10; 
! Maximum number of Newton-Raphson iterations. Newton-Raphson iterations are
! used to build the new grid by finding the coordinates associated with 
! target densities and interpolations of degree larger than 1.
integer, parameter  :: NR_ITERATIONS = 8;
! Tolerance for Newton-Raphson iterations (stop when increment falls below this)
real, parameter     :: NR_TOLERANCE = 1e-12;
! When the N-R algorithm produces an estimate that lies outside [0,1], the
! estimate is set to be equal to the boundary location, 0 or 1, plus or minus
! an offset, respectively, when the derivative is zero at the boundary.
real, parameter     :: NR_OFFSET = 1e-6;

! List of regridding types
integer, parameter :: REGRIDDING_Z         = 0;     ! z* coordinates
integer, parameter :: REGRIDDING_RHO       = 1;     ! target interface densities
integer, parameter :: REGRIDDING_SIGMA     = 2;     ! sigma coordinates
integer, parameter :: REGRIDDING_ARBITRARY = 3;     ! arbitrary coordinates

! List of interpolation schemes
integer, parameter :: INTERPOLATION_P1M_H2     = 0; ! O(h^2)
integer, parameter :: INTERPOLATION_P1M_H4     = 1; ! O(h^2)
integer, parameter :: INTERPOLATION_P1M_IH4    = 2; ! O(h^2)
integer, parameter :: INTERPOLATION_PLM        = 3; ! O(h^2)
integer, parameter :: INTERPOLATION_PPM_H4     = 4; ! O(h^3)
integer, parameter :: INTERPOLATION_PPM_IH4    = 5; ! O(h^3)
integer, parameter :: INTERPOLATION_P3M_IH4IH3 = 6; ! O(h^4)
integer, parameter :: INTERPOLATION_P3M_IH6IH5 = 7; ! O(h^4)
integer, parameter :: INTERPOLATION_PQM_IH4IH3 = 8; ! O(h^4)
integer, parameter :: INTERPOLATION_PQM_IH6IH5 = 9; ! O(h^5)

! List of interpolant degrees
integer, parameter :: DEGREE_1 = 1;
integer, parameter :: DEGREE_2 = 2;
integer, parameter :: DEGREE_3 = 3;
integer, parameter :: DEGREE_4 = 4;

! -----------------------------------------------------------------------------
! Remapping parameters 
! -----------------------------------------------------------------------------

! List of remapping schemes
integer, parameter  :: REMAPPING_PCM        = 0; ! O(h^1)
integer, parameter  :: REMAPPING_PLM        = 1; ! O(h^2)
integer, parameter  :: REMAPPING_PPM_H4     = 2; ! O(h^3)
integer, parameter  :: REMAPPING_PPM_IH4    = 3; ! O(h^3)
integer, parameter  :: REMAPPING_PQM_IH4IH3 = 4; ! O(h^4)
integer, parameter  :: REMAPPING_PQM_IH6IH5 = 5; ! O(h^5)

! These control what routine to use for the remapping integration
integer, parameter  :: INTEGRATION_PCM = 0  ! scope: global
integer, parameter  :: INTEGRATION_PLM = 1  ! scope: global
integer, parameter  :: INTEGRATION_PPM = 3  ! scope: global
integer, parameter  :: INTEGRATION_PQM = 5  ! scope: global

end module regrid_defs
