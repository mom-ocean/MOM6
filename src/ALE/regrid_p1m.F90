module regrid_p1m
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
! p1m_interpolation (public)
! p1m_boundary_extrapolation (public)
!
!==============================================================================
use regrid_grid1d_class   ! see 'regrid_grid1d_class.F90'
use regrid_ppoly_class    ! see 'regrid_ppoly.F90'
use regrid_edge_values    ! see 'regrid_edge_values.F90'

implicit none ; private

! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public p1m_interpolation 
public p1m_boundary_extrapolation

contains


!------------------------------------------------------------------------------
! p1m interpolation
!------------------------------------------------------------------------------
subroutine  p1m_interpolation ( grid, ppoly, u )
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
  type(grid1d_t), intent(in)      :: grid
  type(ppoly_t), intent(inout)    :: ppoly
  real, dimension(:), intent(in)  :: u

  ! Local variables
  integer   :: k;           ! loop index
  integer   :: N;           ! number of cells
  real      :: u0_l, u0_r;  ! edge values (left and right)

  N = grid%nb_cells

  ! Bound edge values (routine found in 'edge_values.F90')
  call bound_edge_values ( grid, u, ppoly%E )
  
  ! Systematically average discontinuous edge values (routine found in
  ! 'edge_values.F90')
  call average_discontinuous_edge_values ( grid, u, ppoly%E )
  
  ! Loop on interior cells to build interpolants
  do k = 1,N
  
    u0_l = ppoly%E(k,1)
    u0_r = ppoly%E(k,2)
    
    ppoly%coefficients(k,1) = u0_l
    ppoly%coefficients(k,2) = u0_r - u0_l
  
  end do ! end loop on interior cells

end subroutine p1m_interpolation


!------------------------------------------------------------------------------
! p1m boundary extrapolation
! -----------------------------------------------------------------------------
subroutine p1m_boundary_extrapolation ( grid, ppoly, u )
!------------------------------------------------------------------------------
! Interpolation by linear polynomials within boundary cells.
! The left and right edge values in the left and right boundary cells,
! respectively, are estimated using a linear extrapolation within the cells.
!
! grid:  one-dimensional grid (properly initialized)
! ppoly: piecewise linear polynomial to be reconstructed (properly initialized)
! u:     cell averages
!
! It is assumed that the size of the array 'u' is equal to the number of cells
! defining 'grid' and 'ppoly'. No consistency check is performed here.
!------------------------------------------------------------------------------

  ! Arguments
  type(grid1d_t), intent(in)      :: grid
  type(ppoly_t), intent(inout)    :: ppoly
  real, dimension(:), intent(in)  :: u

  ! Local variables
  integer       :: k;                   ! loop index
  integer       :: N;                   ! number of cells
  real          :: u0, u1;              ! cell averages
  real          :: h0, h1;              ! corresponding cell widths
  real          :: slope;               ! retained PLM slope
  real          :: a, b;                ! auxiliary variables
  real          :: u0_l, u0_r;          ! edge values

  N = grid%nb_cells

  ! -----------------------------------------
  ! Left edge value in the left boundary cell
  ! -----------------------------------------
  h0 = grid%h(1)
  h1 = grid%h(2)

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

  if ( (u1 - u0) * (ppoly%E(2,1) - u0_r) .LT. 0.0 ) then
    slope = 2.0 * ( ppoly%E(2,1) - u0 )
  end if
  
  ! Using the limited slope, the left edge value is reevaluated and 
  ! the interpolant coefficients recomputed
  if ( h0 .NE. 0.0 ) then
    ppoly%E(1,1) = u0 - 0.5 * slope
  else
    ppoly%E(1,1) = u0
  end if    
  
  ppoly%coefficients(1,1) = ppoly%E(1,1)
  ppoly%coefficients(1,2) = ppoly%E(1,2) - ppoly%E(1,1)
  
  ! ------------------------------------------
  ! Right edge value in the left boundary cell
  ! ------------------------------------------
  h0 = grid%h(N-1)
  h1 = grid%h(N)

  u0 = u(N-1)
  u1 = u(N)

  slope = 2.0 * ( u1 - u0 )
  
  u0_l = u1 - 0.5 * slope

  if ( (u1 - u0) * (u0_l - ppoly%E(N-1,2)) .LT. 0.0 ) then
    slope = 2.0 * ( u1 - ppoly%E(N-1,2) )
  end if
  
  if ( h1 .NE. 0.0 ) then
    ppoly%E(N,2) = u1 + 0.5 * slope
  else
    ppoly%E(N,2) = u1
  end if    
  
  ppoly%coefficients(N,1) = ppoly%E(N,1)
  ppoly%coefficients(N,2) = ppoly%E(N,2) - ppoly%E(N,1)

end subroutine p1m_boundary_extrapolation

end module regrid_p1m
