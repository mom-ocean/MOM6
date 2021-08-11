!> Piecewise linear reconstruction functions
module PLM_functions

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

public PLM_boundary_extrapolation
public PLM_extrapolate_slope
public PLM_monotonized_slope
public PLM_reconstruction
public PLM_slope_wa
public PLM_slope_cw

real, parameter :: hNeglect_dflt = 1.E-30 !< Default negligible cell thickness

contains

!> Returns a limited PLM slope following White and Adcroft, 2008. [units of u]
!! Note that this is not the same as the Colella and Woodward method.
real elemental pure function PLM_slope_wa(h_l, h_c, h_r, h_neglect, u_l, u_c, u_r)
  real, intent(in) :: h_l !< Thickness of left cell [units of grid thickness]
  real, intent(in) :: h_c !< Thickness of center cell [units of grid thickness]
  real, intent(in) :: h_r !< Thickness of right cell [units of grid thickness]
  real, intent(in) :: h_neglect !< A negligible thickness [units of grid thickness]
  real, intent(in) :: u_l !< Value of left cell [units of u]
  real, intent(in) :: u_c !< Value of center cell [units of u]
  real, intent(in) :: u_r !< Value of right cell [units of u]
  ! Local variables
  real :: sigma_l, sigma_c, sigma_r ! Left, central and right slope estimates as
                                    ! differences across the cell [units of u]
  real :: u_min, u_max ! Minimum and maximum value across cell [units of u]

  ! Side differences
  sigma_r = u_r - u_c
  sigma_l = u_c - u_l

  ! Quasi-second order difference
  sigma_c = 2.0 * ( u_r - u_l ) * ( h_c / ( h_l + 2.0*h_c + h_r + h_neglect) )

  ! Limit slope so that reconstructions are bounded by neighbors
  u_min = min( u_l, u_c, u_r )
  u_max = max( u_l, u_c, u_r )
  if ( (sigma_l * sigma_r) > 0.0 ) then
    ! This limits the slope so that the edge values are bounded by the
    ! two cell averages spanning the edge.
    PLM_slope_wa = sign( min( abs(sigma_c), 2.*min( u_c - u_min, u_max - u_c ) ), sigma_c )
  else
    ! Extrema in the mean values require a PCM reconstruction avoid generating
    ! larger extreme values.
    PLM_slope_wa = 0.0
  endif

  ! This block tests to see if roundoff causes edge values to be out of bounds
  if (u_c - 0.5*abs(PLM_slope_wa) < u_min .or.  u_c + 0.5*abs(PLM_slope_wa) > u_max) then
    PLM_slope_wa = PLM_slope_wa * ( 1. - epsilon(PLM_slope_wa) )
  endif

  ! An attempt to avoid inconsistency when the values become unrepresentable.
  ! ### The following 1.E-140 is dimensionally inconsistent. A newer version of
  ! PLM is progress that will avoid the need for such rounding.
  if (abs(PLM_slope_wa) < 1.E-140) PLM_slope_wa = 0.

end function PLM_slope_wa

!> Returns a limited PLM slope following Colella and Woodward 1984.
real elemental pure function PLM_slope_cw(h_l, h_c, h_r, h_neglect, u_l, u_c, u_r)
  real, intent(in) :: h_l !< Thickness of left cell [units of grid thickness]
  real, intent(in) :: h_c !< Thickness of center cell [units of grid thickness]
  real, intent(in) :: h_r !< Thickness of right cell [units of grid thickness]
  real, intent(in) :: h_neglect !< A negligible thickness [units of grid thickness]
  real, intent(in) :: u_l !< Value of left cell [units of u]
  real, intent(in) :: u_c !< Value of center cell [units of u]
  real, intent(in) :: u_r !< Value of right cell [units of u]
  ! Local variables
  real :: sigma_l, sigma_c, sigma_r ! Left, central and right slope estimates as
                                    ! differences across the cell [units of u]
  real :: u_min, u_max ! Minimum and maximum value across cell [units of u]
  real :: h_cn ! Thickness of center cell [units of grid thickness]

  h_cn = h_c + h_neglect

  ! Side differences
  sigma_r = u_r - u_c
  sigma_l = u_c - u_l

  ! This is the second order slope given by equation 1.7 of
  ! Piecewise Parabolic Method, Colella and Woodward (1984),
  ! http://dx.doi.org/10.1016/0021-991(84)90143-8.
  ! For uniform resolution it simplifies to ( u_r - u_l )/2 .
  sigma_c = ( h_c / ( h_cn + ( h_l + h_r ) ) ) * ( &
                ( 2.*h_l + h_c ) / ( h_r + h_cn ) * sigma_r &
              + ( 2.*h_r + h_c ) / ( h_l + h_cn ) * sigma_l )

  ! Limit slope so that reconstructions are bounded by neighbors
  u_min = min( u_l, u_c, u_r )
  u_max = max( u_l, u_c, u_r )
  if ( (sigma_l * sigma_r) > 0.0 ) then
    ! This limits the slope so that the edge values are bounded by the
    ! two cell averages spanning the edge.
    PLM_slope_cw = sign( min( abs(sigma_c), 2.*min( u_c - u_min, u_max - u_c ) ), sigma_c )
  else
    ! Extrema in the mean values require a PCM reconstruction avoid generating
    ! larger extreme values.
    PLM_slope_cw = 0.0
  endif

  ! This block tests to see if roundoff causes edge values to be out of bounds
  if (u_c - 0.5*abs(PLM_slope_cw) < u_min .or.  u_c + 0.5*abs(PLM_slope_cw) > u_max) then
    PLM_slope_cw = PLM_slope_cw * ( 1. - epsilon(PLM_slope_cw) )
  endif

  ! An attempt to avoid inconsistency when the values become unrepresentable.
  ! ### The following 1.E-140 is dimensionally inconsistent. A newer version of
  ! PLM is progress that will avoid the need for such rounding.
  if (abs(PLM_slope_cw) < 1.E-140) PLM_slope_cw = 0.

end function PLM_slope_cw

!> Returns a limited PLM slope following Colella and Woodward 1984.
real elemental pure function PLM_monotonized_slope(u_l, u_c, u_r, s_l, s_c, s_r)
  real, intent(in) :: u_l !< Value of left cell [units of u]
  real, intent(in) :: u_c !< Value of center cell [units of u]
  real, intent(in) :: u_r !< Value of right cell [units of u]
  real, intent(in) :: s_l !< PLM slope of left cell [units of u]
  real, intent(in) :: s_c !< PLM slope of center cell [units of u]
  real, intent(in) :: s_r !< PLM slope of right cell [units of u]
  ! Local variables
  real :: e_r, e_l, edge ! Right, left and temporary edge values [units of u]
  real :: almost_two ! The number 2, almost.
  real :: slp ! Magnitude of PLM central slope [units of u]

  almost_two = 2. * ( 1. - epsilon(s_c) )

  ! Edge values of neighbors abutting this cell
  e_r = u_l + 0.5*s_l
  e_l = u_r - 0.5*s_r
  slp = abs(s_c)

  ! Check that left edge is between right edge of cell to the left and this cell mean
  edge = u_c - 0.5 * s_c
  if ( ( edge - e_r ) * ( u_c - edge ) < 0. ) then
    edge = 0.5 * ( edge + e_r )
    slp = min( slp, abs( edge - u_c ) * almost_two )
  endif

  ! Check that right edge is between left edge of cell to the right and this cell mean
  edge = u_c + 0.5 * s_c
  if ( ( edge - u_c ) * ( e_l - edge ) < 0. ) then
    edge = 0.5 * ( edge + e_l )
    slp = min( slp, abs( edge - u_c ) * almost_two )
  endif

  PLM_monotonized_slope = sign( slp, s_c )

end function PLM_monotonized_slope

!> Returns a PLM slope using h2 extrapolation from a cell to the left.
!! Use the negative to extrapolate from the a cell to the right.
real elemental pure function PLM_extrapolate_slope(h_l, h_c, h_neglect, u_l, u_c)
  real, intent(in) :: h_l !< Thickness of left cell [units of grid thickness]
  real, intent(in) :: h_c !< Thickness of center cell [units of grid thickness]
  real, intent(in) :: h_neglect !< A negligible thickness [units of grid thickness]
  real, intent(in) :: u_l !< Value of left cell [units of u]
  real, intent(in) :: u_c !< Value of center cell [units of u]
  ! Local variables
  real :: left_edge ! Left edge value [units of u]
  real :: hl, hc ! Left and central cell thicknesses [units of grid thickness]

  ! Avoid division by zero for vanished cells
  hl = h_l + h_neglect
  hc = h_c + h_neglect

  ! The h2 scheme is used to compute the left edge value
  left_edge = (u_l*hc + u_c*hl) / (hl + hc)

  PLM_extrapolate_slope = 2.0 * ( u_c - left_edge )

end function PLM_extrapolate_slope


!> Reconstruction by linear polynomials within each cell
!!
!! It is assumed that the size of the array 'u' is equal to the number of cells
!! defining 'grid' and 'ppoly'. No consistency check is performed here.
subroutine PLM_reconstruction( N, h, u, edge_values, ppoly_coef, h_neglect )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N)
  real, dimension(:),   intent(in)    :: u !< cell averages (size N)
  real, dimension(:,:), intent(inout) :: edge_values !< edge values of piecewise polynomials,
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
  real          :: slope                ! retained PLM slope
  real          :: a, b                 ! auxiliary variables
  real          :: u_min, u_max, e_l, e_r, edge
  real          :: almost_one
  real, dimension(N) :: slp, mslp
  real    :: hNeglect

  hNeglect = hNeglect_dflt ; if (present(h_neglect)) hNeglect = h_neglect

  almost_one = 1. - epsilon(slope)

  ! Loop on interior cells
  do k = 2,N-1
    slp(k) = PLM_slope_wa(h(k-1), h(k), h(k+1), hNeglect, u(k-1), u(k), u(k+1))
  enddo ! end loop on interior cells

  ! Boundary cells use PCM. Extrapolation is handled after monotonization.
  slp(1) = 0.
  slp(N) = 0.

  ! This loop adjusts the slope so that edge values are monotonic.
  do K = 2, N-1
    mslp(k) = PLM_monotonized_slope( u(k-1), u(k), u(k+1), slp(k-1), slp(k), slp(k+1) )
  enddo ! end loop on interior cells
  mslp(1) = 0.
  mslp(N) = 0.

  ! Store and return edge values and polynomial coefficients.
  edge_values(1,1) = u(1)
  edge_values(1,2) = u(1)
  ppoly_coef(1,1) = u(1)
  ppoly_coef(1,2) = 0.
  do k = 2, N-1
    slope = mslp(k)
    u_l = u(k) - 0.5 * slope ! Left edge value of cell k
    u_r = u(k) + 0.5 * slope ! Right edge value of cell k

    edge_values(k,1) = u_l
    edge_values(k,2) = u_r
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
  edge_values(N,1) = u(N)
  edge_values(N,2) = u(N)
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
subroutine PLM_boundary_extrapolation( N, h, u, edge_values, ppoly_coef, h_neglect )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N)
  real, dimension(:),   intent(in)    :: u !< cell averages (size N)
  real, dimension(:,:), intent(inout) :: edge_values !< edge values of piecewise polynomials,
                                           !! with the same units as u.
  real, dimension(:,:), intent(inout) :: ppoly_coef !< coefficients of piecewise polynomials, mainly
                                           !! with the same units as u.
  real,       optional, intent(in)    :: h_neglect !< A negligibly small width for
                                           !! the purpose of cell reconstructions
                                           !! in the same units as h
  ! Local variables
  real    :: slope                ! retained PLM slope
  real    :: hNeglect

  hNeglect = hNeglect_dflt ; if (present(h_neglect)) hNeglect = h_neglect

  ! Extrapolate from 2 to 1 to estimate slope
  slope = - PLM_extrapolate_slope( h(2), h(1), hNeglect, u(2), u(1) )

  edge_values(1,1) = u(1) - 0.5 * slope
  edge_values(1,2) = u(1) + 0.5 * slope

  ppoly_coef(1,1) = edge_values(1,1)
  ppoly_coef(1,2) = edge_values(1,2) - edge_values(1,1)

  ! Extrapolate from N-1 to N to estimate slope
  slope = PLM_extrapolate_slope( h(N-1), h(N), hNeglect, u(N-1), u(N) )

  edge_values(N,1) = u(N) - 0.5 * slope
  edge_values(N,2) = u(N) + 0.5 * slope

  ppoly_coef(N,1) = edge_values(N,1)
  ppoly_coef(N,2) = edge_values(N,2) - edge_values(N,1)

end subroutine PLM_boundary_extrapolation

!> \namespace plm_functions
!!
!! Date of creation: 2008.06.06
!! L. White
!!
!! This module contains routines that handle one-dimensionnal finite volume
!! reconstruction using the piecewise linear method (PLM).

end module PLM_functions
