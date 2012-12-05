module regrid_plm
!==============================================================================
!
! This file is part of MOM.
!
! Date of creation: 2008.06.06
! L. White
!
! This module contains routines that handle one-dimensionnal finite volume
! reconstruction using the piecewise linear method (PLM).
!
!==============================================================================
use regrid_grid1d_class    ! see 'regrid_grid1d_class.F90'
use regrid_ppoly_class     ! see 'regrid_ppoly.F90'

implicit none ; private

public plm_reconstruction
public plm_boundary_extrapolation

contains

!------------------------------------------------------------------------------
! plm_reconstruction
! -----------------------------------------------------------------------------
subroutine plm_reconstruction ( grid, ppoly, u )
!------------------------------------------------------------------------------
! Reconstruction by linear polynomials within each cell.
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
  real          :: u_l, u_c, u_r;       ! left, center and right cell averages
  real          :: h_l, h_c, h_r;       ! left, center and right cell widths
  real          :: sigma_l, sigma_c, sigma_r;   ! left, center and right 
                                                ! van Leer slopes   
  real          :: slope;               ! retained PLM slope
  real          :: a, b;                ! auxiliary variables
  real, dimension(:,:), allocatable :: E_old

  N = grid%nb_cells

  ! Loop on interior cells
  do k = 2,N-1
    
    ! Get cell averages
    u_l = u(k-1)
    u_c = u(k)
    u_r = u(k+1)

    ! Get cell widths
    h_l = grid%h(k-1)
    h_c = grid%h(k)
    h_r = grid%h(k+1)

    ! Compute limited slope
    sigma_l = 2.0 * ( u_c - u_l ) / h_c
    sigma_c = 2.0 * ( u_r - u_l ) / ( h_l + 2.0*h_c + h_r )
    sigma_r = 2.0 * ( u_r - u_c ) / h_c

    if ( (sigma_l * sigma_r) .GT. 0.0 ) then
      slope = sign ( min (abs(sigma_l),abs(sigma_c),abs(sigma_r)), sigma_c )
    else
      slope = 0.0
    end if

    ! Determine coefficients of straight line.
    ! CAUTION: the slope 'slope' is computed with respect to the global 
    ! coordinate x. The slope with respect to the normalized local 
    ! coordinate is sigma * h_c (d_P/d_xi = d_P/d_x * d_x/d_xi = d_P/d_x * h_c)
    slope = slope * h_c

    a = u_c - 0.5 * slope
    b = slope

    ppoly%coefficients(k,1) = a
    ppoly%coefficients(k,2) = b
    ppoly%E(k,1) = a
    ppoly%E(k,2) = a+b
        
  end do ! end loop on interior cells

  ! In both boundary cells, a piecewise constant approximation is used.
  ! Extrapolation -- if any -- to increase the accuracy is considered
  ! in another routine, namely 'plm_boundary_extrapolation'
  
  ! Left (top) boundary cell
  slope = u(2) - u(1)
  slope = 0.0
  a = u(1) - 0.5 * slope
  b = slope
  ppoly%coefficients(1,1) = a
  ppoly%coefficients(1,2) = b
  ppoly%E(1,1) = a
  ppoly%E(1,2) = a+b
  
  ! Right (bottom) boundary cell
  slope = u(N) - u(N-1)
  slope = 0.0
  a = u(N) - 0.5 * slope
  b = slope
  ppoly%coefficients(N,1) = a
  ppoly%coefficients(N,2) = b
  ppoly%E(N,1) = a
  ppoly%E(N,2) = a+b

  ! Second pass: we need to check for nonmonotonic discontinuous edge values.
  ! When this occurs, the PLM slope is redefined so as to ensure monotonic edge
  ! values across edges.
  allocate ( E_old(N,2) )
  E_old = ppoly%E
  
  do k = 2,N-1

    ! By default, the right and left slopes are both equal to the original slope
    slope = ppoly%coefficients(k,2)
    sigma_l = slope
    sigma_r = slope
    
    ! If the edge values across the left edge are nonmonotonic (relative to the
    ! cell averages), the PLM slope is redefined in terms of the edge value
    ! lying in the neighboring cell (the left cell).
    if ( (u(k)-u(k-1)) * (E_old(k,1)-E_old(k-1,2)) .LT. 0.0 ) then
      sigma_l = 2.0 * ( u(k) - E_old(k-1,2) )
    end if
    
    ! If the edge values across the right edge are nonmonotonic (relative to the
    ! cell averages), the PLM slope is redefined in terms of the edge value
    ! lying in the neighboring cell (the right cell).
    if ( (u(k+1)-u(k)) * (E_old(k+1,1)-E_old(k,2)) .LT. 0.0 ) then
      sigma_r = 2.0 * ( E_old(k+1,1) - u(k) )
    end if

    ! Take the minimum of both new slopes. If all edge values are 
    ! monotonic, the new slope is simply equal to the old one.
    slope = sign ( min( abs(sigma_l), abs(sigma_r) ), slope )
    
    a = u(k) - 0.5 * slope
    b = slope

    ppoly%coefficients(k,1) = a
    ppoly%coefficients(k,2) = b
    ppoly%E(k,1) = a
    ppoly%E(k,2) = a+b
    
  end do ! end loop on interior cells

  deallocate ( E_old )
  
end subroutine plm_reconstruction


!------------------------------------------------------------------------------
! plm boundary extrapolation
! -----------------------------------------------------------------------------
subroutine plm_boundary_extrapolation ( grid, ppoly, u )
!------------------------------------------------------------------------------
! Reconstruction by linear polynomials within boundary cells.
! The left and right edge values in the left and right boundary cells,
! respectively, are estimated using a linear extrapolation within the cells.
! 
! This extrapolation is EXACT when the underlying profile is linear.
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
  real          :: u0_l, u0_r;          ! edge values

  N = grid%nb_cells

  ! -----------------------------------------
  ! Left edge value in the left boundary cell
  ! -----------------------------------------
  h0 = grid%h(1)
  h1 = grid%h(2)

  u0 = u(1)
  u1 = u(2)

  ! The h2 scheme is used to compute the right edge value
  ppoly%E(1,2) = (u0*h1 + u1*h0) / (h0 + h1)

  ! The standard PLM slope is computed as a first estimate for the
  ! reconstruction within the cell
  slope = 2.0 * ( ppoly%E(1,2) - u0 )

  ppoly%E(1,1) = u0 - 0.5 * slope
  ppoly%E(1,2) = u0 + 0.5 * slope
  
  ppoly%coefficients(1,1) = ppoly%E(1,1)
  ppoly%coefficients(1,2) = ppoly%E(1,2) - ppoly%E(1,1)
  
  ! ------------------------------------------
  ! Right edge value in the left boundary cell
  ! ------------------------------------------
  h0 = grid%h(N-1)
  h1 = grid%h(N)

  u0 = u(N-1)
  u1 = u(N)
  
  ! The h2 scheme is used to compute the right edge value
  ppoly%E(N,1) = (u0*h1 + u1*h0) / (h0 + h1)

  ! The standard PLM slope is computed as a first estimate for the
  ! reconstruction within the cell
  slope = 2.0 * ( u1 - ppoly%E(N,1) )
  
  ppoly%E(N,1) = u1 - 0.5 * slope
  ppoly%E(N,2) = u1 + 0.5 * slope
  
  ppoly%coefficients(N,1) = ppoly%E(N,1)
  ppoly%coefficients(N,2) = ppoly%E(N,2) - ppoly%E(N,1)

end subroutine plm_boundary_extrapolation

end module regrid_plm
