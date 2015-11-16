module PLM_functions
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

implicit none ; private

public PLM_reconstruction, PLM_boundary_extrapolation

real, parameter :: h_neglect = 1.E-30

contains

!------------------------------------------------------------------------------
! PLM_reconstruction
! -----------------------------------------------------------------------------
subroutine PLM_reconstruction( N, h, u, ppoly_E, ppoly_coefficients )
!------------------------------------------------------------------------------
! Reconstruction by linear polynomials within each cell.
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
  integer,              intent(in)    :: N ! Number of cells
  real, dimension(:),   intent(in)    :: h ! cell widths (size N)
  real, dimension(:),   intent(in)    :: u ! cell averages (size N)
  real, dimension(:,:), intent(inout) :: ppoly_E
  real, dimension(:,:), intent(inout) :: ppoly_coefficients

  ! Local variables
  integer       :: k                    ! loop index
  real          :: u_l, u_c, u_r        ! left, center and right cell averages
  real          :: h_l, h_c, h_r        ! left, center and right cell widths
  real          :: sigma_l, sigma_c, sigma_r    ! left, center and right 
                                                ! van Leer slopes   
  real          :: slope                ! retained PLM slope
  real          :: a, b                 ! auxiliary variables
  real, dimension(N,2) :: E_old

  ! Loop on interior cells
  do k = 2,N-1
    
    ! Get cell averages
    u_l = u(k-1)
    u_c = u(k)
    u_r = u(k+1)

    ! Get cell widths
    h_l = h(k-1)
    h_c = h(k)
    h_r = h(k+1)

    ! Compute limited slope
    sigma_l = 2.0 * ( u_c - u_l ) / ( h_c + h_neglect)
    sigma_c = 2.0 * ( u_r - u_l ) / ( h_l + 2.0*h_c + h_r + h_neglect)
    sigma_r = 2.0 * ( u_r - u_c ) / ( h_c + h_neglect)

    if ( (sigma_l * sigma_r) .GT. 0.0 ) then
      slope = sign( min(abs(sigma_l),abs(sigma_c),abs(sigma_r)), sigma_c )
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

    ppoly_coefficients(k,1) = a
    ppoly_coefficients(k,2) = b
    ppoly_E(k,1) = a
    ppoly_E(k,2) = a+b
        
  end do ! end loop on interior cells

  ! In both boundary cells, a piecewise constant approximation is used.
  ! Extrapolation -- if any -- to increase the accuracy is considered
  ! in another routine, namely 'PLM_boundary_extrapolation'
  
  ! Left (top) boundary cell
  slope = u(2) - u(1)
  slope = 0.0
  a = u(1) - 0.5 * slope
  b = slope
  ppoly_coefficients(1,1) = a
  ppoly_coefficients(1,2) = b
  ppoly_E(1,1) = a
  ppoly_E(1,2) = a+b
  
  ! Right (bottom) boundary cell
  slope = u(N) - u(N-1)
  slope = 0.0
  a = u(N) - 0.5 * slope
  b = slope
  ppoly_coefficients(N,1) = a
  ppoly_coefficients(N,2) = b
  ppoly_E(N,1) = a
  ppoly_E(N,2) = a+b

  ! Second pass: we need to check for nonmonotonic discontinuous edge values.
  ! When this occurs, the PLM slope is redefined so as to ensure monotonic edge
  ! values across edges.
  E_old(1:N,1:2) = ppoly_E(1:N,1:2)
  
  do k = 2,N-1

    ! By default, the right and left slopes are both equal to the original slope
    slope = ppoly_coefficients(k,2)
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
    slope = sign( min( abs(sigma_l), abs(sigma_r) ), slope )
    
    a = u(k) - 0.5 * slope
    b = slope

    ppoly_coefficients(k,1) = a
    ppoly_coefficients(k,2) = b
    ppoly_E(k,1) = a
    ppoly_E(k,2) = a+b
    
  end do ! end loop on interior cells

end subroutine PLM_reconstruction


!------------------------------------------------------------------------------
! plm boundary extrapolation
! -----------------------------------------------------------------------------
subroutine PLM_boundary_extrapolation( N, h, u, ppoly_E, ppoly_coefficients )
!------------------------------------------------------------------------------
! Reconstruction by linear polynomials within boundary cells.
! The left and right edge values in the left and right boundary cells,
! respectively, are estimated using a linear extrapolation within the cells.
! 
! This extrapolation is EXACT when the underlying profile is linear.
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
  integer,              intent(in)    :: N ! Number of cells
  real, dimension(:),   intent(in)    :: h ! cell widths (size N)
  real, dimension(:),   intent(in)    :: u ! cell averages (size N)
  real, dimension(:,:), intent(inout) :: ppoly_E
  real, dimension(:,:), intent(inout) :: ppoly_coefficients

  ! Local variables
  real          :: u0, u1               ! cell averages
  real          :: h0, h1               ! corresponding cell widths
  real          :: slope                ! retained PLM slope

  ! -----------------------------------------
  ! Left edge value in the left boundary cell
  ! -----------------------------------------
  h0 = h(1) + h_neglect
  h1 = h(2) + h_neglect

  u0 = u(1)
  u1 = u(2)

  ! The h2 scheme is used to compute the right edge value
  ppoly_E(1,2) = (u0*h1 + u1*h0) / (h0 + h1)

  ! The standard PLM slope is computed as a first estimate for the
  ! reconstruction within the cell
  slope = 2.0 * ( ppoly_E(1,2) - u0 )

  ppoly_E(1,1) = u0 - 0.5 * slope
  ppoly_E(1,2) = u0 + 0.5 * slope
  
  ppoly_coefficients(1,1) = ppoly_E(1,1)
  ppoly_coefficients(1,2) = ppoly_E(1,2) - ppoly_E(1,1)
  
  ! ------------------------------------------
  ! Right edge value in the left boundary cell
  ! ------------------------------------------
  h0 = h(N-1) + h_neglect
  h1 = h(N) + h_neglect

  u0 = u(N-1)
  u1 = u(N)
  
  ! The h2 scheme is used to compute the right edge value
  ppoly_E(N,1) = (u0*h1 + u1*h0) / (h0 + h1)

  ! The standard PLM slope is computed as a first estimate for the
  ! reconstruction within the cell
  slope = 2.0 * ( u1 - ppoly_E(N,1) )
  
  ppoly_E(N,1) = u1 - 0.5 * slope
  ppoly_E(N,2) = u1 + 0.5 * slope
  
  ppoly_coefficients(N,1) = ppoly_E(N,1)
  ppoly_coefficients(N,2) = ppoly_E(N,2) - ppoly_E(N,1)

end subroutine PLM_boundary_extrapolation

end module PLM_functions
