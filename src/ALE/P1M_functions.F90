module P1M_functions
!==============================================================================
!
! This file is part of MOM.
!
! Date of creation: 2008.06.09
! L. White
!
! This module contains p1m (linear) interpolation routines.
!
! p1m interpolation is performed by estimating the edge values and
! linearly interpolating between them.

! Once the edge values are estimated, the limiting process takes care of
! ensuring that (1) edge values are bounded by neighoring cell averages
! and (2) discontinuous edge values are averaged in order to provide a
! fully continuous interpolant throughout the domain. This last step is
! essential for the regridding problem to yield a unique solution.
! Also, a routine is provided that takes care of linear extrapolation
! within the boundary cells.
!
! The module contains the following routines:
!
! P1M_interpolation (public)
! P1M_boundary_extrapolation (public)
!
!==============================================================================
use regrid_edge_values, only : bound_edge_values, average_discontinuous_edge_values

implicit none ; private

! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public P1M_interpolation, P1M_boundary_extrapolation

contains


!------------------------------------------------------------------------------
! p1m interpolation
!------------------------------------------------------------------------------
subroutine P1M_interpolation( N, h, u, ppoly_E, ppoly_coefficients )
! ------------------------------------------------------------------------------
! Linearly interpolate between edge values.
! The resulting piecewise interpolant is stored in 'ppoly'.
! See 'ppoly.F90' for a definition of this structure.
!
! The edge values MUST have been estimated prior to calling this routine.
!
! The estimated edge values must be limited to ensure monotonicity of the
! interpolant. We also make sure that edge values are NOT discontinuous.
!
! It is assumed that the size of the array 'u' is equal to the number of cells
! defining 'grid' and 'ppoly'. No consistency check is performed here.
! ------------------------------------------------------------------------------

  ! Arguments
  integer,            intent(in)    :: N ! Number of cells
  real, dimension(:), intent(in)    :: h ! cell widths (size N)
  real, dimension(:), intent(in)    :: u ! cell averages (size N)
  real, dimension(:,:), intent(inout) :: ppoly_E
  real, dimension(:,:), intent(inout) :: ppoly_coefficients

  ! Local variables
  integer   :: k            ! loop index
  real      :: u0_l, u0_r   ! edge values (left and right)

  ! Bound edge values (routine found in 'edge_values.F90')
  call bound_edge_values( N, h, u, ppoly_E )

  ! Systematically average discontinuous edge values (routine found in
  ! 'edge_values.F90')
  call average_discontinuous_edge_values( N, ppoly_E )

  ! Loop on interior cells to build interpolants
  do k = 1,N

    u0_l = ppoly_E(k,1)
    u0_r = ppoly_E(k,2)

    ppoly_coefficients(k,1) = u0_l
    ppoly_coefficients(k,2) = u0_r - u0_l

  end do ! end loop on interior cells

end subroutine P1M_interpolation


!------------------------------------------------------------------------------
! p1m boundary extrapolation
! -----------------------------------------------------------------------------
subroutine P1M_boundary_extrapolation( N, h, u, ppoly_E, ppoly_coefficients )
!------------------------------------------------------------------------------
! Interpolation by linear polynomials within boundary cells.
! The left and right edge values in the left and right boundary cells,
! respectively, are estimated using a linear extrapolation within the cells.
!
! N:     number of cells in grid
! h:     thicknesses of grid cells
! u:     cell averages to use in constructing piecewise polynomials
! ppoly_E : edge values of piecewise polynomials
! ppoly_coefficients : coefficients of piecewise polynomials
!
! It is assumed that the size of the array 'u' is equal to the number of cells
! defining 'grid' and 'ppoly'. No consistency check is performed here.
!------------------------------------------------------------------------------

  ! Arguments
  integer,            intent(in)    :: N ! Number of cells
  real, dimension(:), intent(in)    :: h ! cell widths (size N)
  real, dimension(:), intent(in)    :: u ! cell averages (size N)
  real, dimension(:,:), intent(inout) :: ppoly_E
  real, dimension(:,:), intent(inout) :: ppoly_coefficients

  ! Local variables
  real          :: u0, u1               ! cell averages
  real          :: h0, h1               ! corresponding cell widths
  real          :: slope                ! retained PLM slope
  real          :: u0_l, u0_r           ! edge values

  ! -----------------------------------------
  ! Left edge value in the left boundary cell
  ! -----------------------------------------
  h0 = h(1)
  h1 = h(2)

  u0 = u(1)
  u1 = u(2)

  ! The standard PLM slope is computed as a first estimate for the
  ! interpolation within the cell
  slope = 2.0 * ( u1 - u0 )

  ! The right edge value is then computed and we check whether this
  ! right edge value is consistent: it cannot be larger than the edge
  ! value in the neighboring cell if the data set is increasing.
  ! If the right value is found to too large, the slope is further limited
  ! by using the edge value in the neighboring cell.
  u0_r = u0 + 0.5 * slope

  if ( (u1 - u0) * (ppoly_E(2,1) - u0_r) .LT. 0.0 ) then
    slope = 2.0 * ( ppoly_E(2,1) - u0 )
  end if

  ! Using the limited slope, the left edge value is reevaluated and
  ! the interpolant coefficients recomputed
  if ( h0 .NE. 0.0 ) then
    ppoly_E(1,1) = u0 - 0.5 * slope
  else
    ppoly_E(1,1) = u0
  end if

  ppoly_coefficients(1,1) = ppoly_E(1,1)
  ppoly_coefficients(1,2) = ppoly_E(1,2) - ppoly_E(1,1)

  ! ------------------------------------------
  ! Right edge value in the left boundary cell
  ! ------------------------------------------
  h0 = h(N-1)
  h1 = h(N)

  u0 = u(N-1)
  u1 = u(N)

  slope = 2.0 * ( u1 - u0 )

  u0_l = u1 - 0.5 * slope

  if ( (u1 - u0) * (u0_l - ppoly_E(N-1,2)) .LT. 0.0 ) then
    slope = 2.0 * ( u1 - ppoly_E(N-1,2) )
  end if

  if ( h1 .NE. 0.0 ) then
    ppoly_E(N,2) = u1 + 0.5 * slope
  else
    ppoly_E(N,2) = u1
  end if

  ppoly_coefficients(N,1) = ppoly_E(N,1)
  ppoly_coefficients(N,2) = ppoly_E(N,2) - ppoly_E(N,1)

end subroutine P1M_boundary_extrapolation

end module P1M_functions
