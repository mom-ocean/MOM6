!> Linear interpolation functions
module P1M_functions

! This file is part of MOM6. See LICENSE.md for the license.

use regrid_edge_values, only : bound_edge_values, average_discontinuous_edge_values

implicit none ; private

! The following routines are visible to the outside world
public P1M_interpolation, P1M_boundary_extrapolation

contains

!> Linearly interpolate between edge values
!!
!! The resulting piecewise interpolant is stored in 'ppoly'.
!! See 'ppoly.F90' for a definition of this structure.
!!
!! The edge values MUST have been estimated prior to calling this routine.
!!
!! The estimated edge values must be limited to ensure monotonicity of the
!! interpolant. We also make sure that edge values are NOT discontinuous.
!!
!! It is assumed that the size of the array 'u' is equal to the number of cells
!! defining 'grid' and 'ppoly'. No consistency check is performed here.
subroutine P1M_interpolation( N, h, u, edge_values, ppoly_coef, h_neglect, answers_2018 )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell average properties (size N) [A]
  real, dimension(:,:), intent(inout) :: edge_values !< Potentially modified edge values [A]
  real, dimension(:,:), intent(inout) :: ppoly_coef !< Potentially modified
                                           !! piecewise polynomial coefficients, mainly [A]
  real,       optional, intent(in)    :: h_neglect !< A negligibly small width [H]
  logical,    optional, intent(in)    :: answers_2018 !< If true use older, less acccurate expressions.

  ! Local variables
  integer   :: k            ! loop index
  real      :: u0_l, u0_r   ! edge values (left and right)

  ! Bound edge values (routine found in 'edge_values.F90')
  call bound_edge_values( N, h, u, edge_values, h_neglect, answers_2018 )

  ! Systematically average discontinuous edge values (routine found in
  ! 'edge_values.F90')
  call average_discontinuous_edge_values( N, edge_values )

  ! Loop on interior cells to build interpolants
  do k = 1,N

    u0_l = edge_values(k,1)
    u0_r = edge_values(k,2)

    ppoly_coef(k,1) = u0_l
    ppoly_coef(k,2) = u0_r - u0_l

  enddo ! end loop on interior cells

end subroutine P1M_interpolation

!> Interpolation by linear polynomials within boundary cells
!!
!!  The left and right edge values in the left and right boundary cells,
!! respectively, are estimated using a linear extrapolation within the cells.
!!
!! It is assumed that the size of the array 'u' is equal to the number of cells
!! defining 'grid' and 'ppoly'. No consistency check is performed here.
subroutine P1M_boundary_extrapolation( N, h, u, edge_values, ppoly_coef )
  ! Arguments
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell averages (size N) [A]
  real, dimension(:,:), intent(inout) :: edge_values !< edge values of piecewise polynomials [A]
  real, dimension(:,:), intent(inout) :: ppoly_coef !< coefficients of piecewise polynomials, mainly [A]

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

  if ( (u1 - u0) * (edge_values(2,1) - u0_r) < 0.0 ) then
    slope = 2.0 * ( edge_values(2,1) - u0 )
  endif

  ! Using the limited slope, the left edge value is reevaluated and
  ! the interpolant coefficients recomputed
  if ( h0 /= 0.0 ) then
    edge_values(1,1) = u0 - 0.5 * slope
  else
    edge_values(1,1) = u0
  endif

  ppoly_coef(1,1) = edge_values(1,1)
  ppoly_coef(1,2) = edge_values(1,2) - edge_values(1,1)

  ! ------------------------------------------
  ! Right edge value in the left boundary cell
  ! ------------------------------------------
  h0 = h(N-1)
  h1 = h(N)

  u0 = u(N-1)
  u1 = u(N)

  slope = 2.0 * ( u1 - u0 )

  u0_l = u1 - 0.5 * slope

  if ( (u1 - u0) * (u0_l - edge_values(N-1,2)) < 0.0 ) then
    slope = 2.0 * ( u1 - edge_values(N-1,2) )
  endif

  if ( h1 /= 0.0 ) then
    edge_values(N,2) = u1 + 0.5 * slope
  else
    edge_values(N,2) = u1
  endif

  ppoly_coef(N,1) = edge_values(N,1)
  ppoly_coef(N,2) = edge_values(N,2) - edge_values(N,1)

end subroutine P1M_boundary_extrapolation

!> \namespace p1m_functions
!!
!! Date of creation: 2008.06.09
!! L. White
!!
!! This module contains p1m (linear) interpolation routines.
!!
!! p1m interpolation is performed by estimating the edge values and
!! linearly interpolating between them.
!
!! Once the edge values are estimated, the limiting process takes care of
!! ensuring that (1) edge values are bounded by neighoring cell averages
!! and (2) discontinuous edge values are averaged in order to provide a
!! fully continuous interpolant throughout the domain. This last step is
!! essential for the regridding problem to yield a unique solution.
!! Also, a routine is provided that takes care of linear extrapolation
!! within the boundary cells.

end module P1M_functions
