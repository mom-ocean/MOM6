!> Cubic interpolation functions
module P3M_functions

! This file is part of MOM6. See LICENSE.md for the license.

use regrid_edge_values, only : bound_edge_values, average_discontinuous_edge_values

implicit none ; private

public P3M_interpolation
public P3M_boundary_extrapolation

real, parameter :: hNeglect_dflt = 1.E-30 !< Default value of a negligible cell thickness
real, parameter :: hNeglect_edge_dflt = 1.E-10 !< Default value of a negligible edge thickness

contains

!> Set up a piecewise cubic interpolation from cell averages and estimated
!! edge slopes and values
!!
!! Cubic interpolation between edges.
!!
!! The edge values and slopes MUST have been estimated prior to calling
!! this routine.
!!
!! It is assumed that the size of the array 'u' is equal to the number of cells
!! defining 'grid' and 'ppoly'. No consistency check is performed here.
subroutine P3M_interpolation( N, h, u, edge_values, ppoly_S, ppoly_coef, h_neglect, answers_2018 )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell averages (size N) in arbitrary units [A]
  real, dimension(:,:), intent(inout) :: edge_values   !< Edge value of polynomial [A]
  real, dimension(:,:), intent(inout) :: ppoly_S   !< Edge slope of polynomial [A H-1].
  real, dimension(:,:), intent(inout) :: ppoly_coef !< Coefficients of polynomial [A]
  real,       optional, intent(in)    :: h_neglect !< A negligibly small width for the
                                          !! purpose of cell reconstructions [H]
  logical,    optional, intent(in)    :: answers_2018 !< If true use older, less acccurate expressions.

  ! Call the limiter for p3m, which takes care of everything from
  ! computing the coefficients of the cubic to monotonizing it.
  ! This routine could be called directly instead of having to call
  ! 'P3M_interpolation' first but we do that to provide an homogeneous
  ! interface.
  call P3M_limiter( N, h, u, edge_values, ppoly_S, ppoly_coef, h_neglect, answers_2018 )

end subroutine P3M_interpolation

!> Adust a piecewise cubic reconstruction with a limiter that adjusts the edge
!! values and slopes
!!
!! The p3m limiter operates as follows:
!!
!! 1. Edge values are bounded
!! 2. Discontinuous edge values are systematically averaged
!! 3. Loop on cells and do the following
!!    a. Build cubic curve
!!    b. Check if cubic curve is monotonic
!!    c. If not, monotonize cubic curve and rebuild it
!!
!! Step 3 of the monotonization process leaves all edge values unchanged.
subroutine P3M_limiter( N, h, u, edge_values, ppoly_S, ppoly_coef, h_neglect, answers_2018 )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell averages (size N) in arbitrary units [A]
  real, dimension(:,:), intent(inout) :: edge_values !< Edge value of polynomial [A]
  real, dimension(:,:), intent(inout) :: ppoly_S  !< Edge slope of polynomial [A H-1]
  real, dimension(:,:), intent(inout) :: ppoly_coef !< Coefficients of polynomial [A]
  real,       optional, intent(in)    :: h_neglect !< A negligibly small width for
                                           !! the purpose of cell reconstructions [H]
  logical,    optional, intent(in)    :: answers_2018 !< If true use older, less acccurate expressions.

  ! Local variables
  integer :: k            ! loop index
  logical :: monotonic    ! boolean indicating whether the cubic is monotonic
  real    :: u0_l, u0_r   ! edge values [A]
  real    :: u1_l, u1_r   ! edge slopes [A H-1]
  real    :: u_l, u_c, u_r        ! left, center and right cell averages [A]
  real    :: h_l, h_c, h_r        ! left, center and right cell widths [H]
  real    :: sigma_l, sigma_c, sigma_r  ! left, center and right van Leer slopes [A H-1]
  real    :: slope        ! retained PLM slope [A H-1]
  real    :: eps
  real    :: hNeglect     ! A negligibly small thickness [H]

  hNeglect = hNeglect_dflt ; if (present(h_neglect)) hNeglect = h_neglect

  eps = 1e-10

  ! 1. Bound edge values (boundary cells are assumed to be local extrema)
  call bound_edge_values( N, h, u, edge_values, hNeglect, answers_2018 )

  ! 2. Systematically average discontinuous edge values
  call average_discontinuous_edge_values( N, edge_values )


  ! 3. Loop on cells and do the following
  !     (a) Build cubic curve
  !     (b) Check if cubic curve is monotonic
  !     (c) If not, monotonize cubic curve and rebuild it
  do k = 1,N

    ! Get edge values, edge slopes and cell width
    u0_l = edge_values(k,1)
    u0_r = edge_values(k,2)
    u1_l = ppoly_S(k,1)
    u1_r = ppoly_S(k,2)

    ! Get cell widths and cell averages (boundary cells are assumed to
    ! be local extrema for the sake of slopes)
    u_c = u(k)
    h_c = h(k)

    if ( k == 1 ) then
      h_l = h(k)
      u_l = u(k)
    else
      h_l = h(k-1)
      u_l = u(k-1)
    endif

    if ( k == N ) then
      h_r = h(k)
      u_r = u(k)
    else
      h_r = h(k+1)
      u_r = u(k+1)
    endif

    ! Compute limited slope
    sigma_l = 2.0 * ( u_c - u_l ) / ( h_c + hNeglect )
    sigma_c = 2.0 * ( u_r - u_l ) / ( h_l + 2.0*h_c + h_r + hNeglect )
    sigma_r = 2.0 * ( u_r - u_c ) / ( h_c + hNeglect )

    if ( (sigma_l * sigma_r) > 0.0 ) then
      slope = sign( min(abs(sigma_l),abs(sigma_c),abs(sigma_r)), sigma_c )
    else
      slope = 0.0
    endif

    ! If the slopes are small, set them to zero to prevent asymmetric representation near extrema.
    if ( abs(u1_l*h_c) < epsilon(u_c)*abs(u_c) ) u1_l = 0.0
    if ( abs(u1_r*h_c) < epsilon(u_c)*abs(u_c) ) u1_r = 0.0

    ! The edge slopes are limited from above by the respective
    ! one-sided slopes
    if ( abs(u1_l) > abs(sigma_l) ) then
      u1_l = sigma_l
    endif

    if ( abs(u1_r) > abs(sigma_r) ) then
      u1_r = sigma_r
    endif

    ! Build cubic interpolant (compute the coefficients)
    call build_cubic_interpolant( h, k, edge_values, ppoly_S, ppoly_coef )

    ! Check whether cubic is monotonic
    monotonic = is_cubic_monotonic( ppoly_coef, k )

    ! If cubic is not monotonic, monotonize it by modifiying the
    ! edge slopes, store the new edge slopes and recompute the
    ! cubic coefficients
    if ( .not.monotonic ) then
      call monotonize_cubic( h_c, u0_l, u0_r, sigma_l, sigma_r, slope, u1_l, u1_r )
    endif

    ! Store edge slopes
    ppoly_S(k,1) = u1_l
    ppoly_S(k,2) = u1_r

    ! Recompute coefficients of cubic
    call build_cubic_interpolant( h, k, edge_values, ppoly_S, ppoly_coef )

  enddo ! loop on cells

end subroutine P3M_limiter


!> Calculate the edge values and slopes at boundary cells as part of building a
!! piecewise cubic sub-grid scale profiles
!!
!! The following explanations apply to the left boundary cell. The same
!! reasoning holds for the right boundary cell.
!!
!! A cubic needs to be built in the cell and requires four degrees of freedom,
!! which are the edge values and slopes. The right edge values and slopes are
!! taken to be that of the neighboring cell (i.e., the left edge value and slope
!! of the neighboring cell). The left edge value and slope are determined by
!! computing the parabola based on the cell average and the right edge value
!! and slope. The resulting cubic is not necessarily monotonic and the slopes
!! are subsequently modified to yield a monotonic cubic.
subroutine P3M_boundary_extrapolation( N, h, u, edge_values, ppoly_S, ppoly_coef, &
                                       h_neglect, h_neglect_edge )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell averages (size N) in arbitrary units [A]
  real, dimension(:,:), intent(inout) :: edge_values !< Edge value of polynomial [A]
  real, dimension(:,:), intent(inout) :: ppoly_S !< Edge slope of polynomial [A H-1]
  real, dimension(:,:), intent(inout) :: ppoly_coef !< Coefficients of polynomial [A]
  real,       optional, intent(in)    :: h_neglect !< A negligibly small width for the
                                          !! purpose of cell reconstructions [H]
  real,       optional, intent(in)    :: h_neglect_edge !< A negligibly small width
                                          !! for the purpose of finding edge values [H]
  ! Local variables
  integer :: i0, i1
  logical :: monotonic    ! boolean indicating whether the cubic is monotonic
  real    :: u0, u1  ! Values of u in two adjacent cells [A]
  real    :: h0, h1  ! Values of h in two adjacent cells, plus a smal increment [H]
  real    :: b, c, d ! Temporary variables [A]
  real    :: u0_l, u0_r ! Left and right edge values [A]
  real    :: u1_l, u1_r ! Left and right edge slopes [A H-1]
  real    :: slope   ! The cell center slope [A H-1]
  real    :: hNeglect, hNeglect_edge ! Negligibly small thickness [H]

  hNeglect = hNeglect_dflt ; if (present(h_neglect)) hNeglect = h_neglect
  hNeglect_edge = hNeglect_edge_dflt ; if (present(h_neglect_edge)) hNeglect_edge = h_neglect_edge

  ! ----- Left boundary -----
  i0 = 1
  i1 = 2
  h0 = h(i0) + hNeglect_edge
  h1 = h(i1) + hNeglect_edge
  u0 = u(i0)
  u1 = u(i1)

  ! Compute the left edge slope in neighboring cell and express it in
  ! the global coordinate system
  b = ppoly_coef(i1,2)
  u1_r = b / h1     ! derivative evaluated at xi = 0.0, expressed w.r.t. x

  ! Limit the right slope by the PLM limited slope
  slope = 2.0 * ( u1 - u0 ) / ( h0 + hNeglect )
  if ( abs(u1_r) > abs(slope) ) then
    u1_r = slope
  endif

  ! The right edge value in the boundary cell is taken to be the left
  ! edge value in the neighboring cell
  u0_r = edge_values(i1,1)

  ! Given the right edge value and slope, we determine the left
  ! edge value and slope by computing the parabola as determined by
  ! the right edge value and slope and the boundary cell average
  u0_l = 3.0 * u0 + 0.5 * h0*u1_r - 2.0 * u0_r
  u1_l = ( - 6.0 * u0 - 2.0 * h0*u1_r + 6.0 * u0_r) / ( h0 + hNeglect )

  ! Check whether the edge values are monotonic. For example, if the left edge
  ! value is larger than the right edge value while the slope is positive, the
  ! edge values are inconsistent and we need to modify the left edge value
  if ( (u0_r-u0_l) * slope < 0.0 ) then
    u0_l = u0_r
    u1_l = 0.0
    u1_r = 0.0
  endif

  ! Store edge values and slope, build cubic and check monotonicity
  edge_values(i0,1) = u0_l
  edge_values(i0,2) = u0_r
  ppoly_S(i0,1) = u1_l
  ppoly_S(i0,2) = u1_r

  ! Store edge values and slope, build cubic and check monotonicity
  call build_cubic_interpolant( h, i0, edge_values, ppoly_S, ppoly_coef )
  monotonic = is_cubic_monotonic( ppoly_coef, i0 )

  if ( .not.monotonic ) then
    call monotonize_cubic( h0, u0_l, u0_r, 0.0, slope, slope, u1_l, u1_r )

    ! Rebuild cubic after monotonization
    ppoly_S(i0,1) = u1_l
    ppoly_S(i0,2) = u1_r
    call build_cubic_interpolant( h, i0, edge_values, ppoly_S, ppoly_coef )

  endif

  ! ----- Right boundary -----
  i0 = N-1
  i1 = N
  h0 = h(i0) + hNeglect_edge
  h1 = h(i1) + hNeglect_edge
  u0 = u(i0)
  u1 = u(i1)

  ! Compute the right edge slope in neighboring cell and express it in
  ! the global coordinate system
  b = ppoly_coef(i0,2)
  c = ppoly_coef(i0,3)
  d = ppoly_coef(i0,4)
  u1_l = (b + 2*c + 3*d) / ( h0 + hNeglect ) ! derivative evaluated at xi = 1.0

  ! Limit the left slope by the PLM limited slope
  slope = 2.0 * ( u1 - u0 ) / ( h1 + hNeglect )
  if ( abs(u1_l) > abs(slope) ) then
    u1_l = slope
  endif

  ! The left edge value in the boundary cell is taken to be the right
  ! edge value in the neighboring cell
  u0_l = edge_values(i0,2)

  ! Given the left edge value and slope, we determine the right
  ! edge value and slope by computing the parabola as determined by
  ! the left edge value and slope and the boundary cell average
  u0_r = 3.0 * u1 - 0.5 * h1*u1_l - 2.0 * u0_l
  u1_r = ( 6.0 * u1 - 2.0 * h1*u1_l - 6.0 * u0_l) / ( h1 + hNeglect )

  ! Check whether the edge values are monotonic. For example, if the right edge
  ! value is smaller than the left edge value while the slope is positive, the
  ! edge values are inconsistent and we need to modify the right edge value
  if ( (u0_r-u0_l) * slope < 0.0 ) then
    u0_r = u0_l
    u1_l = 0.0
    u1_r = 0.0
  endif

  ! Store edge values and slope, build cubic and check monotonicity
  edge_values(i1,1) = u0_l
  edge_values(i1,2) = u0_r
  ppoly_S(i1,1) = u1_l
  ppoly_S(i1,2) = u1_r

  call build_cubic_interpolant( h, i1, edge_values, ppoly_S, ppoly_coef )
  monotonic = is_cubic_monotonic( ppoly_coef, i1 )

  if ( .not.monotonic ) then
    call monotonize_cubic( h1, u0_l, u0_r, slope, 0.0, slope, u1_l, u1_r )

    ! Rebuild cubic after monotonization
    ppoly_S(i1,1) = u1_l
    ppoly_S(i1,2) = u1_r
    call build_cubic_interpolant( h, i1, edge_values, ppoly_S, ppoly_coef )

  endif

end subroutine P3M_boundary_extrapolation


!> Build cubic interpolant in cell k
!!
!! Given edge values and edge slopes, compute coefficients of cubic in cell k.
!!
!! NOTE: edge values and slopes MUST have been properly calculated prior to
!! calling this routine.
subroutine build_cubic_interpolant( h, k, edge_values, ppoly_S, ppoly_coef )
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  integer,              intent(in)    :: k !< The index of the cell to work on
  real, dimension(:,:), intent(in)    :: edge_values !< Edge value of polynomial in arbitrary units [A]
  real, dimension(:,:), intent(in)    :: ppoly_S    !< Edge slope of polynomial [A H-1]
  real, dimension(:,:), intent(inout) :: ppoly_coef !< Coefficients of polynomial [A]

  ! Local variables
  real          :: u0_l, u0_r       ! edge values [A]
  real          :: u1_l, u1_r       ! edge slopes times the cell width [A]
  real          :: h_c              ! cell width  [H]
  real          :: a0, a1, a2, a3   ! cubic coefficients [A]

  h_c = h(k)

  u0_l = edge_values(k,1)
  u0_r = edge_values(k,2)

  u1_l = ppoly_S(k,1) * h_c
  u1_r = ppoly_S(k,2) * h_c

  a0 = u0_l
  a1 = u1_l
  a2 = 3.0 * ( u0_r - u0_l ) - u1_r - 2.0 * u1_l
  a3 = u1_r + u1_l + 2.0 * ( u0_l - u0_r )

  ppoly_coef(k,1) = a0
  ppoly_coef(k,2) = a1
  ppoly_coef(k,3) = a2
  ppoly_coef(k,4) = a3

end subroutine build_cubic_interpolant


!> Check whether the cubic reconstruction in cell k is monotonic
!!
!! This function checks whether the cubic curve in cell k is monotonic.
!! If so, returns 1. Otherwise, returns 0.
!!
!! The cubic is monotonic if the first derivative is single-signed in (0,1).
!! Hence, we check whether the roots (if any) lie inside this interval. If there
!! is no root or if both roots lie outside this interval, the cubic is monotonic.
logical function is_cubic_monotonic( ppoly_coef, k )
  real, dimension(:,:), intent(in) :: ppoly_coef !< Coefficients of cubic polynomial in arbitary units [A]
  integer,              intent(in) :: k  !< The index of the cell to work on
  ! Local variables
  real :: a, b, c   ! Coefficients of the first derivative of the cubic [A]

  a = ppoly_coef(k,2)
  b = 2.0 * ppoly_coef(k,3)
  c = 3.0 * ppoly_coef(k,4)

  ! Look for real roots of the quadratic derivative equation, c*x**2 + b*x + a = 0, in (0, 1)
  if (b*b - 4.0*a*c <= 0.0) then  ! The cubic is monotonic everywhere.
    is_cubic_monotonic = .true.
  elseif (a * (a + (b + c)) < 0.0) then ! The derivative changes sign between the endpoints of (0, 1)
    is_cubic_monotonic = .false.
  elseif (b * (b + 2.0*c) < 0.0) then ! The second derivative changes sign inside of (0, 1)
    is_cubic_monotonic = .false.
  else
    is_cubic_monotonic = .true.
  endif

end function is_cubic_monotonic

!> Monotonize a cubic curve by modifying the edge slopes.
!!
!! This routine takes care of monotonizing a cubic on [0,1] by modifying the
!! edge slopes. The edge values are NOT modified. The cubic is entirely
!! determined by the four degrees of freedom u0_l, u0_r, u1_l and u1_r.
!!
!! u1_l and u1_r are the edge slopes expressed in the GLOBAL coordinate system.
!!
!! The monotonization occurs as follows.
!
!! 1. The edge slopes are set to 0 if they are inconsistent with the limited
!!    PLM slope
!! 2. We check whether we can find an inflexion point in [0,1]. At most one
!!    inflexion point may exist.
!!    a. If there is no inflexion point, the cubic is monotonic.
!!    b. If there is one inflexion point and it lies outside [0,1], the
!!       cubic is monotonic.
!!    c. If there is one inflexion point and it lies in [0,1] and the slope
!!       at the location of the inflexion point is consistent, the cubic
!!       is monotonic.
!!    d. If the inflexion point lies in [0,1] but the slope is inconsistent,
!!       we go to (3) to shift the location of the inflexion point to the left
!!       or to the right. To the left when the 2nd-order left slope is smaller
!!       than the 2nd order right slope.
!! 3. Edge slopes are modified to shift the inflexion point, either onto the left
!!    edge or onto the right edge.

subroutine monotonize_cubic( h, u0_l, u0_r, sigma_l, sigma_r, slope, u1_l, u1_r )
  real, intent(in)      :: h       !< cell width [H]
  real, intent(in)      :: u0_l    !< left edge value in arbitrary units [A]
  real, intent(in)      :: u0_r    !< right edge value [A]
  real, intent(in)      :: sigma_l !< left 2nd-order slopes [A H-1]
  real, intent(in)      :: sigma_r !< right 2nd-order slopes [A H-1]
  real, intent(in)      :: slope   !< limited PLM slope [A H-1]
  real, intent(inout)   :: u1_l    !< left edge slopes [A H-1]
  real, intent(inout)   :: u1_r    !< right edge slopes [A H-1]
  ! Local variables
  logical       :: found_ip
  logical       :: inflexion_l  ! bool telling if inflex. pt must be on left
  logical       :: inflexion_r  ! bool telling if inflex. pt must be on right
  real          :: a1, a2, a3   ! Temporary slopes times the cell width [A]
  real          :: u1_l_tmp     ! trial left edge slope [A H-1]
  real          :: u1_r_tmp     ! trial right edge slope [A H-1]
  real          :: xi_ip        ! location of inflexion point in cell coordinates (0,1) [nondim]
  real          :: slope_ip     ! slope at inflexion point times cell width [A]

  found_ip = .false.
  inflexion_l = .false.
  inflexion_r = .false.

  ! If the edge slopes are inconsistent w.r.t. the limited PLM slope,
  ! set them to zero
  if ( u1_l*slope <= 0.0 ) then
    u1_l = 0.0
  endif

  if ( u1_r*slope <= 0.0 ) then
    u1_r = 0.0
  endif

  ! Compute the location of the inflexion point, which is the root
  ! of the second derivative
  a1 = h * u1_l
  a2 = 3.0 * ( u0_r - u0_l ) - h*(u1_r + 2.0*u1_l)
  a3 = h*(u1_r + u1_l) + 2.0*(u0_l - u0_r)

  ! There is a possible root (and inflexion point) only if a3 is nonzero.
  ! When a3 is zero, the second derivative of the cubic is constant (the
  ! cubic degenerates into a parabola) and no inflexion point exists.
  if ( a3 /= 0.0 ) then
    ! Location of inflexion point
    xi_ip = - a2 / (3.0 * a3)

    ! If the inflexion point lies in [0,1], change boolean value
    if ( (xi_ip >= 0.0) .AND. (xi_ip <= 1.0) ) then
      found_ip = .true.
    endif
  endif

  ! When there is an inflexion point within [0,1], check the slope
  ! to see if it is consistent with the limited PLM slope. If not,
  ! decide on which side we want to collapse the inflexion point.
  ! If the inflexion point lies on one of the edges, the cubic is
  ! guaranteed to be monotonic
  if ( found_ip ) then
    slope_ip = a1 + 2.0*a2*xi_ip + 3.0*a3*xi_ip*xi_ip

    ! Check whether slope is consistent
    if ( slope_ip*slope < 0.0 ) then
      if ( abs(sigma_l) < abs(sigma_r)  ) then
        inflexion_l = .true.
      else
        inflexion_r = .true.
      endif
    endif
  endif ! found_ip

  ! At this point, if the cubic is not monotonic, we know where the
  ! inflexion point should lie. When the cubic is monotonic, both
  ! 'inflexion_l' and 'inflexion_r' are false and nothing is to be done.

  ! Move inflexion point on the left
  if ( inflexion_l ) then

    u1_l_tmp = 1.5*(u0_r-u0_l)/h - 0.5*u1_r
    u1_r_tmp = 3.0*(u0_r-u0_l)/h - 2.0*u1_l

    if ( (u1_l_tmp*slope < 0.0) .AND. (u1_r_tmp*slope < 0.0) ) then

      u1_l = 0.0
      u1_r = 3.0 * (u0_r - u0_l) / h

    elseif (u1_l_tmp*slope < 0.0) then

      u1_r = u1_r_tmp
      u1_l = 1.5*(u0_r - u0_l)/h - 0.5*u1_r

    elseif (u1_r_tmp*slope < 0.0) then

      u1_l = u1_l_tmp
      u1_r = 3.0*(u0_r - u0_l)/h - 2.0*u1_l

    else

      u1_l = u1_l_tmp
      u1_r = u1_r_tmp

    endif

  endif ! end treating case with inflexion point on the left

  ! Move inflexion point on the right
  if ( inflexion_r ) then

    u1_l_tmp = 3.0*(u0_r-u0_l)/h - 2.0*u1_r
    u1_r_tmp = 1.5*(u0_r-u0_l)/h - 0.5*u1_l

    if ( (u1_l_tmp*slope < 0.0) .AND. (u1_r_tmp*slope < 0.0) ) then

      u1_l = 3.0 * (u0_r - u0_l) / h
      u1_r = 0.0

    elseif (u1_l_tmp*slope < 0.0) then

      u1_r = u1_r_tmp
      u1_l = 3.0*(u0_r - u0_l)/h - 2.0*u1_r

    elseif (u1_r_tmp*slope < 0.0) then

      u1_l = u1_l_tmp
      u1_r = 1.5*(u0_r - u0_l)/h - 0.5*u1_l

    else

      u1_l = u1_l_tmp
      u1_r = u1_r_tmp

    endif

  endif ! end treating case with inflexion point on the right

  ! Zero out negligibly small slopes.
  if ( abs(u1_l*h) < epsilon(u0_l) * (abs(u0_l) + abs(u0_r)) ) u1_l = 0.0
  if ( abs(u1_r*h) < epsilon(u0_l) * (abs(u0_l) + abs(u0_r)) ) u1_r = 0.0

end subroutine monotonize_cubic

!> \namespace p3m_functions
!!
!! Date of creation: 2008.06.09
!! L. White
!!
!! This module contains p3m interpolation routines.
!!
!! p3m interpolation is performed by estimating the edge values and slopes
!! and constructing a cubic polynomial. We then make sure that the edge values
!! are bounded and continuous and we then modify the slopes to get a monotonic
!! cubic curve.

end module P3M_functions
