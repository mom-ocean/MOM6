module regrid_edge_values
!==============================================================================
!
! This file is part of MOM.
!
! Date of creation: 2008.06.09
! L. White
!
! This module contains routines that estimate edge values to be used in
! high-order reconstruction schemes.
!
!==============================================================================
use regrid_grid1d_class, only : grid1D_t
use regrid_solvers, only : solve_linear_system, solve_tridiagonal_system
use polynomial_functions, only : evaluation_polynomial

implicit none ; private

! -----------------------------------------------------------------------------
! Private variables used only in this module
! -----------------------------------------------------------------------------
type, public :: edgeValueArrays
  private
  real, dimension(:), allocatable :: tri_l    ! trid. system (lower diagonal)
  real, dimension(:), allocatable :: tri_d    ! trid. system (middle diagonal)
  real, dimension(:), allocatable :: tri_u    ! trid. system (upper diagonal)
  real, dimension(:), allocatable :: tri_x    ! trid. system (unknowns vector)
  real, dimension(:), allocatable :: tri_b    ! trid. system (rhs)
end type edgeValueArrays

! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public bound_edge_values
public average_discontinuous_edge_values
public check_discontinuous_edge_values
public edge_values_explicit_h2
public edge_values_explicit_h4
public edge_values_implicit_h4
public edge_values_implicit_h6
public triDiagEdgeWorkAllocate
public triDiagEdgeWorkDeallocate

contains

!------------------------------------------------------------------------------
! Bound edge values by neighboring cell averages
!------------------------------------------------------------------------------
subroutine bound_edge_values( grid, u, edge_values )
! ------------------------------------------------------------------------------
! In this routine, we loop on all cells to bound their left and right
! edge values by the cell averages. That is, the left edge value must lie
! between the left cell average and the central cell average. A similar
! reasoning applies to the right edge values. 
!
! Both boundary edge values are set equal to the boundary cell averages.
! Any extrapolation scheme is applied after this routine has been called.
! Therefore, boundary cells are treated as if they were local extrama.
! ------------------------------------------------------------------------------

  ! Arguments
  type(grid1D_t), intent(in)            :: grid             
  real, dimension(:), intent(in)        :: u            ! cell averages
  real, dimension(:,:), intent(inout)   :: edge_values      

  ! Local variables
  integer       :: k            ! loop index
  integer       :: N            ! number of cells
  integer       :: k0, k1, k2
  real          :: h_l, h_c, h_r
  real          :: u_l, u_c, u_r
  real          :: u0_l, u0_r
  real          :: sigma_l, sigma_c, sigma_r    ! left, center and right 
                                                ! van Leer slopes   
  real          :: slope                ! retained PLM slope

  N = grid%nb_cells

  ! Loop on cells to bound edge value
  do k = 1,N
    
    ! For the sake of bounding boundary edge values, the left neighbor
    ! of the left boundary cell is assumed to be the same as the left
    ! boundary cell and the right neighbor of the right boundary cell
    ! is assumed to be the same as the right boundary cell. This 
    ! effectively makes boundary cells look like extrema.
    if ( k .EQ. 1 ) then
      k0 = 1
      k1 = 1
      k2 = 2
    else if ( k .EQ. N ) then
      k0 = N-1
      k1 = N
      k2 = N
    else
      k0 = k-1
      k1 = k
      k2 = k+1
    end if

    ! All cells can now be treated equally
    h_l = grid%h(k0)
    h_c = grid%h(k1)
    h_r = grid%h(k2)
    
    u_l = u(k0)
    u_c = u(k1)
    u_r = u(k2)

    u0_l = edge_values(k,1)
    u0_r = edge_values(k,2)

    sigma_l = 2.0 * ( u_c - u_l ) / h_c
    sigma_c = 2.0 * ( u_r - u_l ) / ( h_l + 2.0*h_c + h_r )
    sigma_r = 2.0 * ( u_r - u_c ) / h_c
    
    if ( (sigma_l * sigma_r) .GT. 0.0 ) then
      slope = sign( min(abs(sigma_l),abs(sigma_c),abs(sigma_r)), sigma_c )
    else
      slope = 0.0
    end if
    
    ! The limiter must be used in the local coordinate system to each cell.
    ! Hence, we must multiply the slope by h1. The multiplication by 0.5 is
    ! simply a way to make it useable in the limiter (cfr White and Adcroft 
    ! JCP 2008 Eqs 19 and 20)
    slope = slope * h_c * 0.5
    
    if ( (u_l-u0_l)*(u0_l-u_c) .LT. 0.0 ) then
      u0_l = u_c - sign( min( abs(slope), abs(u0_l-u_c) ), slope )
    end if
    
    if ( (u_r-u0_r)*(u0_r-u_c) .LT. 0.0 ) then
      u0_r = u_c + sign( min( abs(slope), abs(u0_r-u_c) ), slope )
    end if

    ! Store edge values
    edge_values(k,1) = u0_l
    edge_values(k,2) = u0_r

  end do ! loop on interior edges

end subroutine bound_edge_values


!------------------------------------------------------------------------------
! Average discontinuous edge values (systematically)
!------------------------------------------------------------------------------
subroutine average_discontinuous_edge_values( grid, u, edge_values )
! ------------------------------------------------------------------------------
! For each interior edge, check whether the edge values are discontinuous.
! If so, compute the average and replace the edge values by the average.!
! ------------------------------------------------------------------------------
  
  ! Arguments
  type(grid1D_t), intent(in)            :: grid             
  real, dimension(:), intent(in)        :: u            ! cell averages
  real, dimension(:,:), intent(inout)   :: edge_values      

  ! Local variables
  integer       :: k            ! loop index
  integer       :: N            ! number of cells
  real          :: u0_minus     ! left value at given edge
  real          :: u0_plus      ! right value at given edge
  real          :: u0_avg       ! avg value at given edge

  N = grid%nb_cells

  ! Loop on interior edges
  do k = 1,N-1
  
    ! Edge value on the left of the edge
    u0_minus = edge_values(k,2)

    ! Edge value on the right of the edge
    u0_plus  = edge_values(k+1,1)

    if ( u0_minus .NE. u0_plus ) then
      u0_avg = 0.5 * ( u0_minus + u0_plus )
      edge_values(k,2) = u0_avg
      edge_values(k+1,1) = u0_avg
    end if
    
  end do ! end loop on interior edges

end subroutine average_discontinuous_edge_values


!------------------------------------------------------------------------------
! Check discontinuous edge values and take average is not monotonic
!------------------------------------------------------------------------------
subroutine check_discontinuous_edge_values( grid, u, edge_values )
! ------------------------------------------------------------------------------
! For each interior edge, check whether the edge values are discontinuous.
! If so and if they are not monotonic, replace each edge value by their average.
! ------------------------------------------------------------------------------
  
  ! Arguments
  type(grid1D_t), intent(in)            :: grid             
  real, dimension(:), intent(in)        :: u            ! cell averages
  real, dimension(:,:), intent(inout)   :: edge_values      

  ! Local variables
  integer       :: k            ! loop index
  integer       :: N            ! number of cells
  real          :: u0_minus     ! left value at given edge
  real          :: u0_plus      ! right value at given edge
  real          :: um_minus     ! left cell average
  real          :: um_plus      ! right cell average
  real          :: u0_avg       ! avg value at given edge

  N = grid%nb_cells

  ! Loop on interior cells
  do k = 1,N-1
  
    ! Edge value on the left of the edge
    u0_minus = edge_values(k,2)

    ! Edge value on the right of the edge
    u0_plus  = edge_values(k+1,1)

    ! Left cell average
    um_minus = u(k)

    ! Right cell average
    um_plus = u(k+1)

    if ( (u0_plus - u0_minus)*(um_plus - um_minus) .LT. 0.0 ) then
      u0_avg = 0.5 * ( u0_minus + u0_plus )
      edge_values(k,2) = u0_avg
      edge_values(k+1,1) = u0_avg
    end if
    
  end do ! end loop on interior edges

end subroutine check_discontinuous_edge_values


!------------------------------------------------------------------------------
! Compute h2 edge values (explicit second order accurate)
!------------------------------------------------------------------------------
subroutine edge_values_explicit_h2( grid, u, edge_values )
! ------------------------------------------------------------------------------
! Compute edge values based on second-order explicit estimates.
! These estimates are based on a straight line spanning two cells and evaluated
! at the location of the middle edge. An interpolant spanning cells
! k-1 and k is evaluated at edge k-1/2. The estimate for each edge is unique.
!
!       k-1     k   
! ..--o------o------o--..
!          k-1/2
! 
! Boundary edge values are set to be equal to the boundary cell averages.
! ------------------------------------------------------------------------------

  ! Arguments
  type(grid1D_t), intent(in)            :: grid             
  real, dimension(:), intent(in)        :: u            ! cell averages
  real, dimension(:,:), intent(inout)   :: edge_values      

  ! Local variables
  integer   :: k        ! loop index
  integer   :: N        ! number of cells
  real      :: h0, h1   ! cell widths
  real      :: u0, u1   ! cell averages
  
  N = grid%nb_cells
  
  ! Loop on interior cells
  do k = 2,N
  
    h0 = grid%h(k-1)
    h1 = grid%h(k)

    u0 = u(k-1)
    u1 = u(k)

    ! Compute left edge value
    edge_values(k,1) = ( u0*h1 + u1*h0 ) / ( h0 + h1 )

    ! Left edge value of the current cell is equal to right edge 
    ! value of left cell
    edge_values(k-1,2) = edge_values(k,1)
  
  end do ! end loop on interior cells
  
  ! Boundary edge values are simply equal to the boundary cell averages
  edge_values(1,1) = u(1)
  edge_values(N,2) = u(N)

end subroutine edge_values_explicit_h2


!------------------------------------------------------------------------------
! Compute h4 edge values (explicit fourth order accurate)
!------------------------------------------------------------------------------
subroutine edge_values_explicit_h4( grid, u, edge_values )
! -----------------------------------------------------------------------------
! Compute edge values based on fourth-order explicit estimates.
! These estimates are based on a cubic interpolant spanning four cells 
! and evaluated at the location of the middle edge. An interpolant spanning 
! cells i-2, i-1, i and i+1 is evaluated at edge i-1/2. The estimate for 
! each edge is unique.
!
!       i-2    i-1     i     i+1
! ..--o------o------o------o------o--..
!                 i-1/2
! 
! The first two edge values are estimated by evaluating the first available 
! cubic interpolant, i.e., the interpolant spanning cells 1, 2, 3 and 4. 
! Similarly, the last two edge values are estimated by evaluating the last 
! available interpolant.
!
! For this fourth-order scheme, at least four cells must exist.
! -----------------------------------------------------------------------------

  ! Arguments
  type(grid1D_t), intent(in)            :: grid             
  real, dimension(:), intent(in)        :: u            ! cell averages
  real, dimension(:,:), intent(inout)   :: edge_values      

  ! Local variables
  integer               :: i, j
  integer               :: N                ! number of grid cells
  real                  :: u0, u1, u2, u3
  real                  :: h0, h1, h2, h3
  real                  :: href
  real                  :: f1, f2, f3       ! auxiliary variables
  real                  :: e                ! edge value    
  real, dimension(5)    :: x                ! used to compute edge
  real, dimension(4,4)  :: A                ! values near the boundaries
  real, dimension(4)    :: B, C
  
  N = grid%nb_cells
 
  ! Loop on interior cells
  do i = 3,N-1
    
    h0 = grid%h(i-2)
    h1 = grid%h(i-1)
    h2 = grid%h(i)
    h3 = grid%h(i+1)
    
    u0 = u(i-2)
    u1 = u(i-1)
    u2 = u(i)
    u3 = u(i+1)
    
    f1 = (h0+h1) * (h2+h3) / (h1+h2)
    f2 = u1 * h2 + u2 * h1
    f3 = 1.0 / (h0+h1+h2) + 1.0 / (h1+h2+h3)
    
    e = f1 * f2 * f3

    f1 = h2 * (h2+h3) / ( (h0+h1+h2)*(h0+h1) )
    f2 = u1*(h0+2.0*h1) - u0*h1

    e = e + f1*f2
    
    f1 = h1 * (h0+h1) / ( (h1+h2+h3)*(h2+h3) )
    f2 = u2*(2.0*h2+h3) - u3*h2
    
    e = e + f1*f2
    
    e = e / ( h0 + h1 + h2 + h3)
    
    edge_values(i,1) = e
    edge_values(i-1,2) = e

  end do ! end loop on interior cells

  ! Determine first two edge values
  x(1) = 0.0
  do i = 2,5
    x(i) = x(i-1) + grid%h(i-1)
  end do

  do i = 1,4
    
    do j = 1,4
      A(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / j
    end do
    
    B(i) = u(i) * ( grid%h(i) )
    
  end do    

  call solve_linear_system( A, B, C, 4 )

  ! First edge value
  edge_values(1,1) = evaluation_polynomial( C, 4, x(1) )
  
  ! Second edge value
  edge_values(1,2) = evaluation_polynomial( C, 4, x(2) )
  edge_values(2,1) = edge_values(1,2)

  ! Determine last two edge values
  x(1) = 0.0
  do i = 2,5
    x(i) = x(i-1) + grid%h(N-5+i)
  end do

  do i = 1,4
    
    do j = 1,4
      A(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / j
    end do
    
    B(i) = u(N-4+i) * ( grid%h(N-4+i) )
    
  end do    

  call solve_linear_system( A, B, C, 4 )
  
  ! Last edge value
  edge_values(N,2) = evaluation_polynomial( C, 4, x(5) )
  
  ! Second to last edge value
  edge_values(N,1) = evaluation_polynomial( C, 4, x(4) )
  
  edge_values(N-1,2) = edge_values(N,1)

end subroutine edge_values_explicit_h4


!------------------------------------------------------------------------------
! Compute ih4 edge values (implicit fourth order accurate)
!------------------------------------------------------------------------------
subroutine edge_values_implicit_h4( grid, work, u, edge_values )
! -----------------------------------------------------------------------------
! Compute edge values based on fourth-order implicit estimates.
!
! Fourth-order implicit estimates of edge values are based on a two-cell 
! stencil. A tridiagonal system is set up and is based on expressing the 
! edge values in terms of neighboring cell averages. The generic 
! relationship is
!
! \alpha u_{i-1/2} + u_{i+1/2} + \beta u_{i+3/2} = a \bar{u}_i + b \bar{u}_{i+1}
!
! and the stencil looks like this
!
!          i     i+1
!   ..--o------o------o--..
!     i-1/2  i+1/2  i+3/2
!
! In this routine, the coefficients \alpha, \beta, a and b are computed, 
! the tridiagonal system is built, boundary conditions are prescribed and 
! the system is solved to yield edge-value estimates. 
! 
! There are N+1 unknowns and we are able to write N-1 equations. The 
! boundary conditions close the system.
! -----------------------------------------------------------------------------

  ! Arguments
  type(grid1D_t), intent(in)            :: grid             
  type(edgeValueArrays), intent(inout)  :: work         ! Work space
  real, dimension(:), intent(in)        :: u            ! cell averages
  real, dimension(:,:), intent(inout)   :: edge_values      

  ! Local variables
  integer               :: i, j                 ! loop indexes
  integer               :: N                    ! number of cells
  real                  :: h0, h1               ! cell widths
  real                  :: h0_2, h1_2, h0h1
  real                  :: d2, d4
  real                  :: alpha, beta          ! stencil coefficients
  real                  :: a, b
  real, dimension(5)    :: x                    ! system used to enforce
  real, dimension(4,4)  :: Asys                 ! boundary conditions
  real, dimension(4)    :: Bsys, Csys

  ! Get number of cells (there are N+1 edge values to estimate)
  N = grid%nb_cells

  ! Loop on cells (except last one)
  do i = 1,N-1

    ! Get cell widths
    h0 = grid%h(i)
    h1 = grid%h(i+1)

    ! Auxiliary calculations
    d2 = (h0 + h1) ** 2
    d4 = d2 ** 2
    h0h1 = h0 * h1
    h0_2 = h0 * h0
    h1_2 = h1 * h1

    ! Coefficients
    alpha = h1_2 / d2
    beta = h0_2 / d2
    a = 2.0 * h1_2 * ( h1_2 + 2.0 * h0_2 + 3.0 * h0h1 ) / d4
    b = 2.0 * h0_2 * ( h0_2 + 2.0 * h1_2 + 3.0 * h0h1 ) / d4
    
    work%tri_l(i+1) = alpha
    work%tri_d(i+1) = 1.0
    work%tri_u(i+1) = beta
    
    work%tri_b(i+1) = a * u(i) + b * u(i+1)
    
  end do ! end loop on cells

  ! Boundary conditions: left boundary
  x(1) = 0.0
  do i = 2,5
    x(i) = x(i-1) + grid%h(i-1)
  end do

  do i = 1,4
    
    do j = 1,4
      Asys(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / j
    end do
    
    Bsys(i) = u(i) * ( grid%h(i) )
    
  end do    

  call solve_linear_system( Asys, Bsys, Csys, 4 )
  
  work%tri_d(1) = 1.0
  work%tri_u(1) = 0.0
  work%tri_b(1) = evaluation_polynomial( Csys, 4, x(1) )        ! first edge value
  
  ! Boundary conditions: right boundary
  x(1) = 0.0
  do i = 2,5
    x(i) = x(i-1) + grid%h(N-5+i)
  end do

  do i = 1,4
    
    do j = 1,4
      Asys(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / j
    end do
    
    Bsys(i) = u(N-4+i) * ( grid%h(N-4+i) )
    
  end do    

  call solve_linear_system( Asys, Bsys, Csys, 4 )
  
  work%tri_l(N+1) = 0.0
  work%tri_d(N+1) = 1.0
  work%tri_b(N+1) = evaluation_polynomial( Csys, 4, x(5) )      ! last edge value

  ! Solve tridiagonal system and assign edge values
  call solve_tridiagonal_system( work%tri_l, work%tri_d, work%tri_u, work%tri_b, work%tri_x, N+1 )

  do i = 2,N
    edge_values(i,1)   = work%tri_x(i)
    edge_values(i-1,2) = work%tri_x(i)
  end do
  edge_values(1,1) = work%tri_x(1)
  edge_values(N,2) = work%tri_x(N+1)

end subroutine edge_values_implicit_h4


!------------------------------------------------------------------------------
! Compute ih6 edge values (implicit sixth order accurate)
!------------------------------------------------------------------------------
subroutine edge_values_implicit_h6( grid, work, u, edge_values )
! -----------------------------------------------------------------------------
! Sixth-order implicit estimates of edge values are based on a four-cell, 
! three-edge stencil. A tridiagonal system is set up and is based on 
! expressing the edge values in terms of neighboring cell averages. 
! 
! The generic relationship is
!
! \alpha u_{i-1/2} + u_{i+1/2} + \beta u_{i+3/2} = 
! a \bar{u}_{i-1} + b \bar{u}_i + c \bar{u}_{i+1} + d \bar{u}_{i+2}
!
! and the stencil looks like this
!
!         i-1     i     i+1    i+2
!   ..--o------o------o------o------o--..
!            i-1/2  i+1/2  i+3/2
!
! In this routine, the coefficients \alpha, \beta, a, b, c and d are 
! computed, the tridiagonal system is built, boundary conditions are 
! prescribed and the system is solved to yield edge-value estimates. 
!
! Note that the centered stencil only applies to edges 3 to N-1 (edges are
! numbered 1 to n+1), which yields N-3 equations for N+1 unknowns. Two other
! equations are written by using a right-biased stencil for edge 2 and a
! left-biased stencil for edge N. The prescription of boundary conditions
! (using sixth-order polynomials) closes the system.
!
! CAUTION: For each edge, in order to determine the coefficients of the
!          implicit expression, a 6x6 linear system is solved. This may
!          become computationally expensive if regridding is carried out
!          often. Figuring out closed-form expressions for these coefficients
!          on nonuniform meshes turned out to be intractable.
! -----------------------------------------------------------------------------

  ! Arguments
  type(grid1D_t), intent(in)            :: grid             
  type(edgeValueArrays), intent(inout)  :: work         ! Work space
  real, dimension(:), intent(in)        :: u            ! cell averages
  real, dimension(:,:), intent(inout)   :: edge_values      

  ! Local variables
  integer               :: i, j, k              ! loop indexes
  integer               :: N                    ! number of cells
  real                  :: h0, h1, h2, h3       ! cell widths
  real                  :: g, g_2, g_3          ! the following are 
  real                  :: g_4, g_5, g_6        ! auxiliary variables
  real                  :: d2, d3, d4, d5, d6   ! to set up the systems
  real                  :: n2, n3, n4, n5, n6   ! used to compute the
  real                  :: h1_2, h2_2           ! the coefficients of the
  real                  :: h1_3, h2_3           ! tridiagonal system
  real                  :: h1_4, h2_4           ! ...
  real                  :: h1_5, h2_5           ! ...
  real                  :: h1_6, h2_6           ! ...
  real                  :: h0ph1, h0ph1_2       ! ...
  real                  :: h0ph1_3, h0ph1_4     ! ...
  real                  :: h2ph3, h2ph3_2       ! ...
  real                  :: h2ph3_3, h2ph3_4     ! ...
  real                  :: h0ph1_5, h2ph3_5     ! ...
  real                  :: alpha, beta          ! stencil coefficients
  real                  :: a, b, c, d           ! "
  real, dimension(7)    :: x                    ! system used to enforce
  real, dimension(6,6)  :: Asys                 ! boundary conditions
  real, dimension(6)    :: Bsys, Csys           ! ...

  ! Get number of cells (there are N+1 edge values to estimate)
  N = grid%nb_cells

  ! Loop on cells (except last one)
  do k = 2,N-2

    ! Cell widths
    h0 = grid%h(k-1)
    h1 = grid%h(k+0)
    h2 = grid%h(k+1)
    h3 = grid%h(k+2)

    ! Auxiliary calculations
    h1_2 = h1 * h1
    h1_3 = h1_2 * h1 
    h1_4 = h1_2 * h1_2
    h1_5 = h1_3 * h1_2
    h1_6 = h1_3 * h1_3
    
    h2_2 = h2 * h2
    h2_3 = h2_2 * h2
    h2_4 = h2_2 * h2_2
    h2_5 = h2_3 * h2_2
    h2_6 = h2_3 * h2_3
    
    g   = h0 + h1
    g_2 = g * g
    g_3 = g * g_2
    g_4 = g_2 * g_2
    g_5 = g_4 * g
    g_6 = g_3 * g_3

    d2 = ( h1_2 - g_2 ) / h0
    d3 = ( h1_3 - g_3 ) / h0
    d4 = ( h1_4 - g_4 ) / h0
    d5 = ( h1_5 - g_5 ) / h0
    d6 = ( h1_6 - g_6 ) / h0
    
    g   = h2 + h3
    g_2 = g * g
    g_3 = g * g_2
    g_4 = g_2 * g_2
    g_5 = g_4 * g
    g_6 = g_3 * g_3
    
    n2 = ( g_2 - h2_2 ) / h3
    n3 = ( g_3 - h2_3 ) / h3
    n4 = ( g_4 - h2_4 ) / h3
    n5 = ( g_5 - h2_5 ) / h3
    n6 = ( g_6 - h2_6 ) / h3

    ! Compute matrix entries
    Asys(1,1) = 1.0
    Asys(1,2) = 1.0
    Asys(1,3) = -1.0
    Asys(1,4) = -1.0
    Asys(1,5) = -1.0
    Asys(1,6) = -1.0
    
    Asys(2,1) = - h1
    Asys(2,2) = h2
    Asys(2,3) = -0.5 * d2
    Asys(2,4) = 0.5 * h1
    Asys(2,5) = -0.5 * h2
    Asys(2,6) = -0.5 * n2
    
    Asys(3,1) = 0.5 * h1_2
    Asys(3,2) = 0.5 * h2_2
    Asys(3,3) = d3 / 6.0
    Asys(3,4) = - h1_2 / 6.0
    Asys(3,5) = - h2_2 / 6.0
    Asys(3,6) = - n3 / 6.0
    
    Asys(4,1) = - h1_3 / 6.0
    Asys(4,2) = h2_3 / 6.0
    Asys(4,3) = - d4 / 24.0
    Asys(4,4) = h1_3 / 24.0
    Asys(4,5) = - h2_3 / 24.0
    Asys(4,6) = - n4 / 24.0
    
    Asys(5,1) = h1_4 / 24.0
    Asys(5,2) = h2_4 / 24.0
    Asys(5,3) = d5 / 120.0
    Asys(5,4) = - h1_4 / 120.0
    Asys(5,5) = - h2_4 / 120.0
    Asys(5,6) = - n5 / 120.0
    
    Asys(6,1) = - h1_5 / 120.0
    Asys(6,2) = h2_5 / 120.0
    Asys(6,3) = - d6 / 720.0
    Asys(6,4) = h1_5 / 720.0
    Asys(6,5) = - h2_5 / 720.0
    Asys(6,6) = - n6 / 720.0
    
    Bsys(:) = (/ -1.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
    call solve_linear_system( Asys, Bsys, Csys, 6 )

    alpha = Csys(1)
    beta  = Csys(2)
    a = Csys(3)
    b = Csys(4)
    c = Csys(5)
    d = Csys(6)

    work%tri_l(k+1) = alpha
    work%tri_d(k+1) = 1.0
    work%tri_u(k+1) = beta
    work%tri_b(k+1) = a * u(k-1) + b * u(k) + c * u(k+1) + d * u(k+2)
    
  end do ! end loop on cells
  
  ! Use a right-biased stencil for the second row
  
  ! Cell widths
  h0 = grid%h(1)
  h1 = grid%h(2)
  h2 = grid%h(3)
  h3 = grid%h(4)
    
  ! Auxiliary calculations
  h1_2 = h1 * h1
  h1_3 = h1_2 * h1
  h1_4 = h1_2 * h1_2
  h1_5 = h1_3 * h1_2
  h1_6 = h1_3 * h1_3
    
  h2_2 = h2 * h2
  h2_3 = h2_2 * h2
  h2_4 = h2_2 * h2_2
  h2_5 = h2_3 * h2_2
  h2_6 = h2_3 * h2_3
    
  g   = h0 + h1
  g_2 = g * g
  g_3 = g * g_2
  g_4 = g_2 * g_2
  g_5 = g_4 * g
  g_6 = g_3 * g_3

  h0ph1   = h0 + h1
  h0ph1_2 = h0ph1 * h0ph1
  h0ph1_3 = h0ph1_2 * h0ph1
  h0ph1_4 = h0ph1_2 * h0ph1_2
  h0ph1_5 = h0ph1_3 * h0ph1_2

  d2 = ( h1_2 - g_2 ) / h0
  d3 = ( h1_3 - g_3 ) / h0
  d4 = ( h1_4 - g_4 ) / h0
  d5 = ( h1_5 - g_5 ) / h0
  d6 = ( h1_6 - g_6 ) / h0
    
  g   = h2 + h3
  g_2 = g * g
  g_3 = g * g_2
  g_4 = g_2 * g_2
  g_5 = g_4 * g
  g_6 = g_3 * g_3
    
  n2 = ( g_2 - h2_2 ) / h3
  n3 = ( g_3 - h2_3 ) / h3
  n4 = ( g_4 - h2_4 ) / h3
  n5 = ( g_5 - h2_5 ) / h3
  n6 = ( g_6 - h2_6 ) / h3
    
  ! Compute matrix entries
  Asys(1,1) = 1.0
  Asys(1,2) = 1.0
  Asys(1,3) = -1.0
  Asys(1,4) = -1.0
  Asys(1,5) = -1.0
  Asys(1,6) = -1.0
    
  Asys(2,1) = - h0ph1
  Asys(2,2) = 0.0
  Asys(2,3) = -0.5 * d2
  Asys(2,4) = 0.5 * h1
  Asys(2,5) = -0.5 * h2
  Asys(2,6) = -0.5 * n2
    
  Asys(3,1) = 0.5 * h0ph1_2
  Asys(3,2) = 0.0
  Asys(3,3) = d3 / 6.0
  Asys(3,4) = - h1_2 / 6.0
  Asys(3,5) = - h2_2 / 6.0
  Asys(3,6) = - n3 / 6.0
    
  Asys(4,1) = - h0ph1_3 / 6.0
  Asys(4,2) = 0.0
  Asys(4,3) = - d4 / 24.0
  Asys(4,4) = h1_3 / 24.0
  Asys(4,5) = - h2_3 / 24.0
  Asys(4,6) = - n4 / 24.0
    
  Asys(5,1) = h0ph1_4 / 24.0
  Asys(5,2) = 0.0
  Asys(5,3) = d5 / 120.0
  Asys(5,4) = - h1_4 / 120.0
  Asys(5,5) = - h2_4 / 120.0
  Asys(5,6) = - n5 / 120.0
    
  Asys(6,1) = - h0ph1_5 / 120.0
  Asys(6,2) = 0.0
  Asys(6,3) = - d6 / 720.0
  Asys(6,4) = h1_5 / 720.0
  Asys(6,5) = - h2_5 / 720.0
  Asys(6,6) = - n6 / 720.0
    
  Bsys(:) = (/ -1.0, h1, -0.5*h1_2, h1_3/6.0, -h1_4/24.0, h1_5/120.0 /)

  call solve_linear_system( Asys, Bsys, Csys, 6 )

  alpha = Csys(1)
  beta  = Csys(2)
  a = Csys(3)
  b = Csys(4)
  c = Csys(5)
  d = Csys(6)

  work%tri_l(2) = alpha
  work%tri_d(2) = 1.0
  work%tri_u(2) = beta
  work%tri_b(2) = a * u(1) + b * u(2) + c * u(3) + d * u(4)
  
  ! Boundary conditions: left boundary
  x(1) = 0.0
  do i = 2,7
    x(i) = x(i-1) + grid%h(i-1)
  end do

  do i = 1,6
    
    do j = 1,6
      Asys(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / j
    end do
    
    Bsys(i) = u(i) * grid%h(i)
    
  end do    

  call solve_linear_system( Asys, Bsys, Csys, 6 )
  
  work%tri_l(1) = 0.0
  work%tri_d(1) = 1.0
  work%tri_u(1) = 0.0
  work%tri_b(1) = evaluation_polynomial( Csys, 6, x(1) )        ! first edge value
  
  ! Use a left-biased stencil for the second to last row
  
  ! Cell widths
  h0 = grid%h(N-3)
  h1 = grid%h(N-2)
  h2 = grid%h(N-1)
  h3 = grid%h(N)
    
  ! Auxiliary calculations
  h1_2 = h1 * h1
  h1_3 = h1_2 * h1
  h1_4 = h1_2 * h1_2
  h1_5 = h1_3 * h1_2
  h1_6 = h1_3 * h1_3
    
  h2_2 = h2 * h2
  h2_3 = h2_2 * h2
  h2_4 = h2_2 * h2_2
  h2_5 = h2_3 * h2_2
  h2_6 = h2_3 * h2_3
    
  g   = h0 + h1
  g_2 = g * g
  g_3 = g * g_2
  g_4 = g_2 * g_2
  g_5 = g_4 * g
  g_6 = g_3 * g_3

  h2ph3   = h2 + h3
  h2ph3_2 = h2ph3 * h2ph3
  h2ph3_3 = h2ph3_2 * h2ph3
  h2ph3_4 = h2ph3_2 * h2ph3_2
  h2ph3_5 = h2ph3_3 * h2ph3_2

  d2 = ( h1_2 - g_2 ) / h0
  d3 = ( h1_3 - g_3 ) / h0
  d4 = ( h1_4 - g_4 ) / h0
  d5 = ( h1_5 - g_5 ) / h0
  d6 = ( h1_6 - g_6 ) / h0
    
  g   = h2 + h3
  g_2 = g * g
  g_3 = g * g_2
  g_4 = g_2 * g_2
  g_5 = g_4 * g
  g_6 = g_3 * g_3
    
  n2 = ( g_2 - h2_2 ) / h3
  n3 = ( g_3 - h2_3 ) / h3
  n4 = ( g_4 - h2_4 ) / h3
  n5 = ( g_5 - h2_5 ) / h3
  n6 = ( g_6 - h2_6 ) / h3

  ! Compute matrix entries
  Asys(1,1) = 1.0
  Asys(1,2) = 1.0
  Asys(1,3) = -1.0
  Asys(1,4) = -1.0
  Asys(1,5) = -1.0
  Asys(1,6) = -1.0
    
  Asys(2,1) =   0.0
  Asys(2,2) = h2ph3
  Asys(2,3) =   -0.5 * d2
  Asys(2,4) =   0.5 * h1
  Asys(2,5) =   -0.5 * h2
  Asys(2,6) =   -0.5 * n2
    
  Asys(3,1) =   0.0
  Asys(3,2) = 0.5 * h2ph3_2
  Asys(3,3) =   d3 / 6.0
  Asys(3,4) =   - h1_2 / 6.0
  Asys(3,5) =   - h2_2 / 6.0
  Asys(3,6) =   - n3 / 6.0
    
  Asys(4,1) =   0.0
  Asys(4,2) = h2ph3_3 / 6.0
  Asys(4,3) =   - d4 / 24.0
  Asys(4,4) =   h1_3 / 24.0
  Asys(4,5) =   - h2_3 / 24.0
  Asys(4,6) =   - n4 / 24.0
    
  Asys(5,1) =   0.0
  Asys(5,2) = h2ph3_4 / 24.0
  Asys(5,3) =   d5 / 120.0
  Asys(5,4) =   - h1_4 / 120.0
  Asys(5,5) =   - h2_4 / 120.0
  Asys(5,6) =   - n5 / 120.0
    
  Asys(6,1) =   0.0
  Asys(6,2) = h2ph3_5 / 120.0
  Asys(6,3) =   - d6 / 720.0
  Asys(6,4) =   h1_5 / 720.0
  Asys(6,5) =   - h2_5 / 720.0
  Asys(6,6) =   - n6 / 720.0

  Bsys(:) = (/ -1.0, -h2, -0.5*h2_2, -h2_3/6.0, -h2_4/24.0, -h2_5/120.0 /)

  call solve_linear_system( Asys, Bsys, Csys, 6 )

  alpha = Csys(1)
  beta  = Csys(2)
  a = Csys(3)
  b = Csys(4)
  c = Csys(5)
  d = Csys(6)

  work%tri_l(N) = alpha
  work%tri_d(N) = 1.0
  work%tri_u(N) = beta
  work%tri_b(N) = a * u(N-3) + b * u(N-2) + c * u(N-1) + d * u(N)

  ! Boundary conditions: right boundary
  x(1) = 0.0
  do i = 2,7
    x(i) = x(i-1) + grid%h(N-7+i)
  end do

  do i = 1,6
    
    do j = 1,6
      Asys(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / j
    end do
    
    Bsys(i) = u(N-6+i) * grid%h(N-6+i)
    
  end do    

  call solve_linear_system( Asys, Bsys, Csys, 6 )
  
  work%tri_l(N+1) = 0.0
  work%tri_d(N+1) = 1.0
  work%tri_u(N+1) = 0.0
  work%tri_b(N+1) = evaluation_polynomial( Csys, 6, x(7) )      ! last edge value
  
  ! Solve tridiagonal system and assign edge values
  call solve_tridiagonal_system( work%tri_l, work%tri_d, work%tri_u, work%tri_b, work%tri_x, N+1 )

  do i = 2,N
    edge_values(i,1)   = work%tri_x(i)
    edge_values(i-1,2) = work%tri_x(i)
  end do
  edge_values(1,1) = work%tri_x(1)
  edge_values(N,2) = work%tri_x(N+1)

end subroutine edge_values_implicit_h6


!------------------------------------------------------------------------------
! Allocate memory for tridiagonal system used to compute edge values
!------------------------------------------------------------------------------
subroutine triDiagEdgeWorkAllocate( N, work )
!------------------------------------------------------------------------------
! In this routine, we allocate memory for the tridiagonal system that will
! be used to compute implicit edge-value estimates. The argument 'N' is the
! number of layers and therefore, the size of the system is N+1.
!------------------------------------------------------------------------------

  ! Argument
  integer, intent(in)   :: N
  type(edgeValueArrays), intent(inout)  :: work          ! Work space

  allocate( work%tri_l(N+1) )
  allocate( work%tri_d(N+1) )
  allocate( work%tri_u(N+1) )
  allocate( work%tri_b(N+1) )
  allocate( work%tri_x(N+1) )

end subroutine triDiagEdgeWorkAllocate

!------------------------------------------------------------------------------
! Deallocate memory for tridiagonal system used to compute edge values
!------------------------------------------------------------------------------
subroutine triDiagEdgeWorkDeallocate( work )
  type(edgeValueArrays), intent(inout)  :: work          ! Work space

  deallocate( work%tri_l )
  deallocate( work%tri_d )
  deallocate( work%tri_u )
  deallocate( work%tri_b )
  deallocate( work%tri_x )

end subroutine triDiagEdgeWorkDeallocate


end module regrid_edge_values
