module regrid_ppm
!==============================================================================
!
! This file is part of MOM.
!
! Date of creation: 2008.06.06
! L. White
!
! This module contains routines that handle one-dimensionnal finite volume
! reconstruction using the piecewise parabolic method (PPM).
!
!==============================================================================
use regrid_grid1d_class ! see 'regrid_grid1d_class.F90'
use regrid_ppoly_class  ! see 'regrid_ppoly.F90'
use regrid_edge_values  ! see 'regrid_edge_values.F90'

implicit none ; private

public ppm_reconstruction
public ppm_boundary_extrapolation

contains

!------------------------------------------------------------------------------
! ppm_reconstruction
! -----------------------------------------------------------------------------
subroutine ppm_reconstruction ( grid, ppoly, u )
!------------------------------------------------------------------------------
! Reconstruction by quadratic polynomials within each cell.
!
! The edge values MUST have been estimated prior to calling this routine
!
! grid:  one-dimensional grid (see grid.F90)
! ppoly: piecewise quadratic polynomial to be reconstructed (see ppoly.F90)
! u:     cell averages
!
! It is assumed that the dimension of 'u' is equal to the number of cells
! defining 'grid' and 'ppoly'. No consistency check is performed.
!------------------------------------------------------------------------------

  ! Arguments
  type(grid1d_t), intent(in)     :: grid
  type(ppoly_t), intent(inout)   :: ppoly
  real, dimension(:), intent(in) :: u

  ! Local variables
  integer   :: k;           ! loop index
  integer   :: N;           ! number of cells
  real      :: u0_l, u0_r;  ! edge values (left and right)
  real      :: a, b, c;     ! parabola coefficients
  
  N = grid%nb_cells

  ! PPM limiter
  call ppm_limiter_standard ( grid, ppoly, u )

  ! Loop on cells to construct the parabola within each cell
  do k = 1,N
  
    u0_l = ppoly%E(k,1)
    u0_r = ppoly%E(k,2)

    a = u0_l
    b = 6.0 * u(k) - 4.0 * u0_l - 2.0 * u0_r
    c = 3.0 * ( u0_r + u0_l - 2.0 * u(k) )
    
    ! Store coefficients
    ppoly%coefficients(k,1) = a
    ppoly%coefficients(k,2) = b
    ppoly%coefficients(k,3) = c
  
  end do ! end loop on interior cells

end subroutine ppm_reconstruction


!------------------------------------------------------------------------------
! Limit ppm
! -----------------------------------------------------------------------------
subroutine ppm_limiter_standard ( grid, ppoly, u )
!------------------------------------------------------------------------------
! Standard PPM limiter (Colella & Woodward, JCP 1984).
!
! grid:  one-dimensional grid (see grid.F90)
! ppoly: piecewise quadratic polynomial to be reconstructed (see ppoly.F90)
! u:     cell averages
!
! It is assumed that the dimension of 'u' is equal to the number of cells
! defining 'grid' and 'ppoly'. No consistency check is performed.
!------------------------------------------------------------------------------

  ! Arguments
  type(grid1d_t), intent(in)     :: grid
  type(ppoly_t), intent(inout)   :: ppoly
  real, dimension(:), intent(in) :: u

  ! Local variables
  integer   :: k;               ! loop index
  integer   :: N;               ! number of cells
  real      :: u_l, u_c, u_r;   ! cell averages (left, center and right)
  real      :: u0_l, u0_r;      ! edge values (left and right)
  real      :: expr1, expr2
  
  N = grid%nb_cells

  ! Bound edge values
  call bound_edge_values ( grid, u, ppoly%E )

  ! Make discontinuous edge values monotonic
  call check_discontinuous_edge_values ( grid, u, ppoly%E )

  ! Loop on interior cells to apply the standard 
  ! PPM limiter (Colella & Woodward, JCP 84)
  do k = 2,N-1
    
    ! Get cell averages
    u_l = u(k-1)
    u_c = u(k)
    u_r = u(k+1)
  
    u0_l = ppoly%E(k,1)
    u0_r = ppoly%E(k,2)
    
    ! Auxiliary variables
    expr1 = (u0_r - u0_l) * (u_c - 0.5*(u0_l+u0_r))
    expr2 = (u0_r - u0_l) * (u0_r - u0_l) / 6.0
    
    ! Flatten extremum
    if ( (u_r - u_c)*(u_c - u_l) .LE. 0.0) then
      u0_l = u_c
      u0_r = u_c
    end if

    if ( expr1 .GT. expr2 ) then  
      u0_l = 3.0 * u_c - 2.0 * u0_r
    end if  
    
    if ( expr1 .LT. -expr2 ) then  
      u0_r = 3.0 * u_c - 2.0 * u0_l
    end if  
    
    ppoly%E(k,1) = u0_l
    ppoly%E(k,2) = u0_r
  
  end do ! end loop on interior cells

  ! PCM within boundary cells
  ppoly%E(1,:) = u(1)
  ppoly%E(N,:) = u(N)

end subroutine ppm_limiter_standard


!------------------------------------------------------------------------------
! ppm boundary extrapolation
! -----------------------------------------------------------------------------
subroutine ppm_boundary_extrapolation ( grid, ppoly, u )
!------------------------------------------------------------------------------
! Reconstruction by parabolas within boundary cells.
!
! The following explanations apply to the left boundary cell. The same
! reasoning holds for the right boundary cell.
!
! A parabola needs to be built in the cell and requires three degrees of
! freedom, which are the right edge value and slope and the cell average.
! The right edge values and slopes are taken to be that of the neighboring
! cell (i.e., the left edge value and slope of the neighboring cell).
! The resulting parabola is not necessarily monotonic and the traditional
! PPM limiter is used to modify one of the edge values in order to yield
! a monotonic parabola.
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
  integer       :: k;       ! loop index
  integer       :: N;       ! number of cells
  integer       :: i0, i1
  integer       :: monotonic
  real          :: u0, u1
  real          :: h0, h1
  real          :: a, b, c, d, e
  real          :: u0_l, u0_r
  real          :: u1_l, u1_r
  real          :: slope
  real          :: exp1, exp2

  N = grid%nb_cells
  
  ! ----- Left boundary -----
  i0 = 1
  i1 = 2
  h0 = grid%h(i0)
  h1 = grid%h(i1)
  u0 = u(i0)
  u1 = u(i1)

  ! Compute the left edge slope in neighboring cell and express it in
  ! the global coordinate system
  b = ppoly%coefficients(i1,2)
  u1_r = b *(h0/h1);    ! derivative evaluated at xi = 0.0, 
                        ! expressed w.r.t. xi (local coord. system)
  
  ! Limit the right slope by the PLM limited slope
  slope = 2.0 * ( u1 - u0 )
  if ( abs(u1_r) .GT. abs(slope) ) then
    u1_r = slope
  end if

  ! The right edge value in the boundary cell is taken to be the left
  ! edge value in the neighboring cell
  u0_r = ppoly%E(i1,1)

  ! Given the right edge value and slope, we determine the left
  ! edge value and slope by computing the parabola as determined by
  ! the right edge value and slope and the boundary cell average
  u0_l = 3.0 * u0 + 0.5 * u1_r - 2.0 * u0_r

  ! Apply the traditional PPM limiter
  exp1 = (u0_r - u0_l) * (u0 - 0.5*(u0_l+u0_r))
  exp2 = (u0_r - u0_l) * (u0_r - u0_l) / 6.0
    
  if ( exp1 .GT. exp2 ) then  
    u0_l = 3.0 * u0 - 2.0 * u0_r
  end if  

  if ( exp1 .LT. -exp2 ) then  
    u0_r = 3.0 * u0 - 2.0 * u0_l
  end if  

  ppoly%E(i0,1) = u0_l
  ppoly%E(i0,2) = u0_r
    
  a = u0_l
  b = 6.0 * u0 - 4.0 * u0_l - 2.0 * u0_r
  c = 3.0 * ( u0_r + u0_l - 2.0 * u0 )

  ppoly%coefficients(i0,1) = a
  ppoly%coefficients(i0,2) = b
  ppoly%coefficients(i0,3) = c
  
  ! ----- Right boundary -----
  i0 = N-1
  i1 = N
  h0 = grid%h(i0)
  h1 = grid%h(i1)
  u0 = u(i0)
  u1 = u(i1)

  ! Compute the right edge slope in neighboring cell and express it in
  ! the global coordinate system
  b = ppoly%coefficients(i0,2)
  c = ppoly%coefficients(i0,3)
  u1_l = (b + 2*c);                 ! derivative evaluated at xi = 1.0
  u1_l = u1_l * (h1/h0)
  
  ! Limit the left slope by the PLM limited slope
  slope = 2.0 * ( u1 - u0 )
  if ( abs(u1_l) .GT. abs(slope) ) then
    u1_l = slope
  end if

  ! The left edge value in the boundary cell is taken to be the right
  ! edge value in the neighboring cell
  u0_l = ppoly%E(i0,2)

  ! Given the left edge value and slope, we determine the right
  ! edge value and slope by computing the parabola as determined by
  ! the left edge value and slope and the boundary cell average
  u0_r = 3.0 * u1 - 0.5 * u1_l - 2.0 * u0_l
  
  ! Apply the traditional PPM limiter
  exp1 = (u0_r - u0_l) * (u1 - 0.5*(u0_l+u0_r))
  exp2 = (u0_r - u0_l) * (u0_r - u0_l) / 6.0
    
  if ( exp1 .GT. exp2 ) then  
    u0_l = 3.0 * u1 - 2.0 * u0_r
  end if  

  if ( exp1 .LT. -exp2 ) then  
    u0_r = 3.0 * u1 - 2.0 * u0_l
  end if  

  ppoly%E(i1,1) = u0_l
  ppoly%E(i1,2) = u0_r
    
  a = u0_l
  b = 6.0 * u1 - 4.0 * u0_l - 2.0 * u0_r
  c = 3.0 * ( u0_r + u0_l - 2.0 * u1 )

  ppoly%coefficients(i1,1) = a
  ppoly%coefficients(i1,2) = b
  ppoly%coefficients(i1,3) = c

end subroutine ppm_boundary_extrapolation

end module regrid_ppm
