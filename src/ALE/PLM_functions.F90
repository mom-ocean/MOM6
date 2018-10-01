!> Piecewise linear reconstruction functions
module PLM_functions

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

public PLM_reconstruction, PLM_boundary_extrapolation

real, parameter :: hNeglect_dflt = 1.E-30 !< Default negligible cell thickness

contains

!> Reconstruction by linear polynomials within each cell
!!
!! It is assumed that the size of the array 'u' is equal to the number of cells
!! defining 'grid' and 'ppoly'. No consistency check is performed here.
subroutine PLM_reconstruction( N, h, u, ppoly_E, ppoly_coef, h_neglect )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N)
  real, dimension(:),   intent(in)    :: u !< cell averages (size N)
  real, dimension(:,:), intent(inout) :: ppoly_E !< edge values of piecewise polynomials,
                                           !! with the same units as u.
  real, dimension(:,:), intent(inout) :: ppoly_coef !< coefficients of piecewise polynomials, mainly
                                           !! with the same units as u.
  real,       optional, intent(in)    :: h_neglect !< A negligibly small width for
                                           !! the purpose of cell reconstructions
                                           !! in the same units as h

  ! Local variables
  integer       :: k                    ! loop index
  real          :: u_l, u_c, u_r        ! left, center and right cell averages
  real          :: h_l, h_c, h_r, h_cn  ! left, center and right cell widths
  real          :: sigma_l, sigma_c, sigma_r    ! left, center and right
                                                ! van Leer slopes
  real          :: slope                ! retained PLM slope
  real          :: a, b                 ! auxiliary variables
  real          :: u_min, u_max, e_l, e_r, edge
  real          :: almost_one, almost_two
  real, dimension(N) :: slp, mslp
  real    :: hNeglect

  hNeglect = hNeglect_dflt ; if (present(h_neglect)) hNeglect = h_neglect

  almost_one = 1. - epsilon(slope)
  almost_two = 2. * almost_one

  ! Loop on interior cells
  do k = 2,N-1

    ! Get cell averages
    u_l = u(k-1) ; u_c = u(k) ; u_r = u(k+1)

    ! Get cell widths
    h_l = h(k-1) ; h_c = h(k) ; h_r = h(k+1)
    h_cn = max( h_c, hNeglect ) ! To avoid division by zero

    ! Side differences
    sigma_r = u_r - u_c
    sigma_l = u_c - u_l

    ! This is the second order slope given by equation 1.7 of
    ! Piecewise Parabolic Method, Colella and Woodward (1984),
    ! http://dx.doi.org/10.1016/0021-991(84)90143-8.
    ! For uniform resolution it simplifies to ( u_r - u_l )/2 .
 !  sigma_c = ( h_c / ( h_cn + ( h_l + h_r ) ) ) * ( &
 !                ( 2.*h_l + h_c ) / ( h_r + h_cn ) * sigma_r &
 !              + ( 2.*h_r + h_c ) / ( h_l + h_cn ) * sigma_l )

    ! This is the original estimate of the second order slope from Laurent
    ! but multiplied by h_c
    sigma_c = 2.0 * ( u_r - u_l ) * ( h_c / ( h_l + 2.0*h_c + h_r + hNeglect) )

    if ( (sigma_l * sigma_r) > 0.0 ) then
      ! This limits the slope so that the edge values are bounded by the
      ! two cell averages spanning the edge.
      u_min = min( u_l, u_c, u_r )
      u_max = max( u_l, u_c, u_r )
      slope = sign( min( abs(sigma_c), 2.*min( u_c - u_min, u_max - u_c ) ), sigma_c )
    else
      ! Extrema in the mean values require a PCM reconstruction avoid generating
      ! larger extreme values.
      slope = 0.0
    endif

    ! This block tests to see if roundoff causes edge values to be out of bounds
    u_min = min( u_l, u_c, u_r )
    u_max = max( u_l, u_c, u_r )
    if (u_c - 0.5*abs(slope) < u_min .or.  u_c + 0.5*abs(slope) > u_max) then
      slope = slope * almost_one
    endif

    ! An attempt to avoid inconsistency when the values become unrepresentable.
    if (abs(slope) < 1.E-140) slope = 0.

    ! Safety check - this block really should not be needed ...
!   if (u_c - 0.5*abs(slope) < u_min .or.  u_c + 0.5*abs(slope) > u_max) then
!     write(0,*) 'l,c,r=',u_l,u_c,u_r
!     write(0,*) 'min,max=',u_min,u_max
!     write(0,*) 'slp=',slope
!     sigma_l = u_c-0.5*abs(slope)
!     sigma_r = u_c+0.5*abs(slope)
!     write(0,*) 'lo,hi=',sigma_l,sigma_r
!     write(0,*) 'elo,ehi=',sigma_l-u_min,sigma_r-u_max
!     stop 'Limiter failed!'
!   endif

    slp(k) = slope
    ppoly_E(k,1) = u_c - 0.5 * slope
    ppoly_E(k,2) = u_c + 0.5 * slope

  enddo ! end loop on interior cells

  ! Boundary cells use PCM. Extrapolation is handled in a later routine.
  slp(1) = 0.
  ppoly_E(1,2) = u(1)
  slp(N) = 0.
  ppoly_E(N,1) = u(N)

  ! This loop adjusts the slope so that edge values are monotonic.
  do K = 2, N-1
    u_l = u(k-1) ; u_c = u(k) ; u_r = u(k+1)
    e_r = ppoly_E(k-1,2) ! Right edge from cell k-1
    e_l = ppoly_E(k+1,1) ! Left edge from cell k
    mslp(k) = abs(slp(k))
    u_min = min(e_r, u_c)
    u_max = max(e_r, u_c)
    edge = u_c - 0.5 * slp(k)
    if ( ( edge - e_r ) * ( u_c - edge ) < 0. ) then
      edge = 0.5 * ( edge + e_r ) ! * almost_one?
      mslp(k) = min( mslp(k), abs( edge - u_c ) * almost_two )
    endif
    edge = u_c + 0.5 * slp(k)
    if ( ( edge - u_c ) * ( e_l - edge ) < 0. ) then
      edge = 0.5 * ( edge + e_l ) ! * almost_one?
      mslp(k) = min( mslp(k), abs( edge - u_c ) * almost_two )
    endif
  enddo ! end loop on interior cells
  mslp(1) = 0.
  mslp(N) = 0.

  ! Check that the above adjustment worked
! do K = 2, N-1
!   u_r = u(k-1) + 0.5 * sign( mslp(k-1), slp(k-1) ) ! Right edge from cell k-1
!   u_l = u(k) - 0.5 * sign( mslp(k), slp(k) ) ! Left edge from cell k
!   if ( (u(k)-u(k-1)) * (u_l-u_r) < 0. ) then
!     stop 'Adjustment failed!'
!   endif
! enddo ! end loop on interior cells

  ! Store and return edge values and polynomial coefficients.
  ppoly_E(1,1) = u(1)
  ppoly_E(1,2) = u(1)
  ppoly_coef(1,1) = u(1)
  ppoly_coef(1,2) = 0.
  do k = 2, N-1
    slope = sign( mslp(k), slp(k) )
    u_l = u(k) - 0.5 * slope ! Left edge value of cell k
    u_r = u(k) + 0.5 * slope ! Right edge value of cell k

    ! Check that final edge values are bounded
    u_min = min( u(k-1), u(k) )
    u_max = max( u(k-1), u(k) )
    if (u_l<u_min .or. u_l>u_max) then
      write(0,*) 'u(k-1)=',u(k-1),'u(k)=',u(k),'slp=',slp(k),'u_l=',u_l
      stop 'Left edge out of bounds'
    endif
    u_min = min( u(k+1), u(k) )
    u_max = max( u(k+1), u(k) )
    if (u_r<u_min .or. u_r>u_max) then
      write(0,*) 'u(k)=',u(k),'u(k+1)=',u(k+1),'slp=',slp(k),'u_r=',u_r
      stop 'Right edge out of bounds'
    endif

    ppoly_E(k,1) = u_l
    ppoly_E(k,2) = u_r
    ppoly_coef(k,1) = u_l
    ppoly_coef(k,2) = ( u_r - u_l )
    ! Check to see if this evaluation of the polynomial at x=1 would be
    ! monotonic w.r.t. the next cell's edge value. If not, scale back!
    edge = ppoly_coef(k,2) + ppoly_coef(k,1)
    e_r = u(k+1) - 0.5 * sign( mslp(k+1), slp(k+1) )
    if ( (edge-u(k))*(e_r-edge)<0.) then
      ppoly_coef(k,2) = ppoly_coef(k,2) * almost_one
    endif
  enddo
  ppoly_E(N,1) = u(N)
  ppoly_E(N,2) = u(N)
  ppoly_coef(N,1) = u(N)
  ppoly_coef(N,2) = 0.

end subroutine PLM_reconstruction


!> Reconstruction by linear polynomials within boundary cells
!!
!! The left and right edge values in the left and right boundary cells,
!! respectively, are estimated using a linear extrapolation within the cells.
!!
!! This extrapolation is EXACT when the underlying profile is linear.
!!
!! It is assumed that the size of the array 'u' is equal to the number of cells
!! defining 'grid' and 'ppoly'. No consistency check is performed here.

subroutine PLM_boundary_extrapolation( N, h, u, ppoly_E, ppoly_coef, h_neglect )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N)
  real, dimension(:),   intent(in)    :: u !< cell averages (size N)
  real, dimension(:,:), intent(inout) :: ppoly_E !< edge values of piecewise polynomials,
                                           !! with the same units as u.
  real, dimension(:,:), intent(inout) :: ppoly_coef !< coefficients of piecewise polynomials, mainly
                                           !! with the same units as u.
  real,       optional, intent(in)    :: h_neglect !< A negligibly small width for
                                           !! the purpose of cell reconstructions
                                           !! in the same units as h

  ! Local variables
  real    :: u0, u1               ! cell averages
  real    :: h0, h1               ! corresponding cell widths
  real    :: slope                ! retained PLM slope
  real    :: hNeglect

  hNeglect = hNeglect_dflt ; if (present(h_neglect)) hNeglect = h_neglect

  ! -----------------------------------------
  ! Left edge value in the left boundary cell
  ! -----------------------------------------
  h0 = h(1) + hNeglect
  h1 = h(2) + hNeglect

  u0 = u(1)
  u1 = u(2)

  ! The h2 scheme is used to compute the right edge value
  ppoly_E(1,2) = (u0*h1 + u1*h0) / (h0 + h1)

  ! The standard PLM slope is computed as a first estimate for the
  ! reconstruction within the cell
  slope = 2.0 * ( ppoly_E(1,2) - u0 )

  ppoly_E(1,1) = u0 - 0.5 * slope
  ppoly_E(1,2) = u0 + 0.5 * slope

  ppoly_coef(1,1) = ppoly_E(1,1)
  ppoly_coef(1,2) = ppoly_E(1,2) - ppoly_E(1,1)

  ! ------------------------------------------
  ! Right edge value in the left boundary cell
  ! ------------------------------------------
  h0 = h(N-1) + hNeglect
  h1 = h(N) + hNeglect

  u0 = u(N-1)
  u1 = u(N)

  ! The h2 scheme is used to compute the right edge value
  ppoly_E(N,1) = (u0*h1 + u1*h0) / (h0 + h1)

  ! The standard PLM slope is computed as a first estimate for the
  ! reconstruction within the cell
  slope = 2.0 * ( u1 - ppoly_E(N,1) )

  ppoly_E(N,1) = u1 - 0.5 * slope
  ppoly_E(N,2) = u1 + 0.5 * slope

  ppoly_coef(N,1) = ppoly_E(N,1)
  ppoly_coef(N,2) = ppoly_E(N,2) - ppoly_E(N,1)

end subroutine PLM_boundary_extrapolation

!> \namespace plm_functions
!!
!! Date of creation: 2008.06.06
!! L. White
!!
!! This module contains routines that handle one-dimensionnal finite volume
!! reconstruction using the piecewise linear method (PLM).

end module PLM_functions
