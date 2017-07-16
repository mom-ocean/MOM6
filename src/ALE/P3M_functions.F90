module P3M_functions
!==============================================================================
!
! This file is part of MOM.
!
! Date of creation: 2008.06.09
! L. White
!
! This module contains p3m interpolation routines.
!
! p3m interpolation is performed by estimating the edge values and slopes
! and constructing a cubic polynomial. We then make sure that the edge values
! are bounded and continuous and we then modify the slopes to get a monotonic
! cubic curve.
!
!==============================================================================
use regrid_edge_values, only : bound_edge_values, average_discontinuous_edge_values

implicit none ; private

public P3M_interpolation
public P3M_boundary_extrapolation

real, parameter :: h_neglect = 1.E-30

contains

!------------------------------------------------------------------------------
! p3m interpolation
! -----------------------------------------------------------------------------
subroutine P3M_interpolation( N, h, u, ppoly_E, ppoly_S, ppoly_coefficients )
!------------------------------------------------------------------------------
! Cubic interpolation between edges.
!
! The edge values and slopes MUST have been estimated prior to calling
! this routine.
!
! It is assumed that the size of the array 'u' is equal to the number of cells
! defining 'grid' and 'ppoly'. No consistency check is performed here.
!------------------------------------------------------------------------------

  ! Arguments
  integer,              intent(in)    :: N ! Number of cells
  real, dimension(:),   intent(in)    :: h ! cell widths (size N)
  real, dimension(:),   intent(in)    :: u ! cell averages (size N)
  real, dimension(:,:), intent(inout) :: ppoly_E            !Edge value of polynomial
  real, dimension(:,:), intent(inout) :: ppoly_S            !Edge slope of polynomial
  real, dimension(:,:), intent(inout) :: ppoly_coefficients !Coefficients of polynomial


  ! Call the limiter for p3m, which takes care of everything from
  ! computing the coefficients of the cubic to monotonizing it.
  ! This routine could be called directly instead of having to call
  ! 'P3M_interpolation' first but we do that to provide an homogeneous
  ! interface.

  call P3M_limiter( N, h, u, ppoly_E, ppoly_S, ppoly_coefficients )

end subroutine P3M_interpolation


!------------------------------------------------------------------------------
! p3m limiter
! -----------------------------------------------------------------------------
subroutine P3M_limiter( N, h, u, ppoly_E, ppoly_S, ppoly_coefficients )
!------------------------------------------------------------------------------
! The p3m limiter operates as follows:
!
! (1) Edge values are bounded
! (2) Discontinuous edge values are systematically averaged
! (3) Loop on cells and do the following
!       (a) Build cubic curve
!       (b) Check if cubic curve is monotonic
!       (c) If not, monotonize cubic curve and rebuild it
!
! Step (3) of the monotonization process leaves all edge values unchanged.
!------------------------------------------------------------------------------

  ! Arguments
  integer,              intent(in)    :: N ! Number of cells
  real, dimension(:),   intent(in)    :: h ! cell widths (size N)
  real, dimension(:),   intent(in)    :: u ! cell averages (size N)
  real, dimension(:,:), intent(inout) :: ppoly_E            !Edge value of polynomial
  real, dimension(:,:), intent(inout) :: ppoly_S            !Edge slope of polynomial
  real, dimension(:,:), intent(inout) :: ppoly_coefficients !Coefficients of polynomial

!  real, dimension(:,:), intent(inout) :: ppoly_coefficients

  ! Local variables
  integer   :: k            ! loop index
  integer   :: monotonic    ! boolean indicating whether the cubic is monotonic
  real      :: u0_l, u0_r   ! edge values
  real      :: u1_l, u1_r   ! edge slopes
  real      :: u_l, u_c, u_r        ! left, center and right cell averages
  real      :: h_l, h_c, h_r        ! left, center and right cell widths
  real      :: sigma_l, sigma_c, sigma_r    ! left, center and right
                                            ! van Leer slopes
  real      :: slope        ! retained PLM slope
  real      :: eps

  eps = 1e-10

  ! 1. Bound edge values (boundary cells are assumed to be local extrema)
  call bound_edge_values( N, h, u, ppoly_E )

  ! 2. Systematically average discontinuous edge values
  call average_discontinuous_edge_values( N, ppoly_E )


  ! 3. Loop on cells and do the following
  !     (a) Build cubic curve
  !     (b) Check if cubic curve is monotonic
  !     (c) If not, monotonize cubic curve and rebuild it
  do k = 1,N

    ! Get edge values, edge slopes and cell width
    u0_l = ppoly_E(k,1)
    u0_r = ppoly_E(k,2)
    u1_l = ppoly_S(k,1)
    u1_r = ppoly_S(k,2)

    ! Get cell widths and cell averages (boundary cells are assumed to
    ! be local extrema for the sake of slopes)
    u_c = u(k)
    h_c = h(k)

    if ( k .EQ. 1 ) then
      h_l = h(k)
      u_l = u(k)
    else
      h_l = h(k-1)
      u_l = u(k-1)
    end if

    if ( k .EQ. N ) then
      h_r = h(k)
      u_r = u(k)
    else
      h_r = h(k+1)
      u_r = u(k+1)
    end if

    ! Compute limited slope
    sigma_l = 2.0 * ( u_c - u_l ) / ( h_c + h_neglect )
    sigma_c = 2.0 * ( u_r - u_l ) / ( h_l + 2.0*h_c + h_r + h_neglect )
    sigma_r = 2.0 * ( u_r - u_c ) / ( h_c + h_neglect )

    if ( (sigma_l * sigma_r) .GT. 0.0 ) then
      slope = sign( min(abs(sigma_l),abs(sigma_c),abs(sigma_r)), sigma_c )
    else
      slope = 0.0
    end if

    ! If the slopes are close to zero in machine precision and in absolute
    ! value, we set the slope to zero. This prevents asymmetric representation
    ! near extrema.
    if ( abs(u1_l*h_c) .LT. eps ) then
      u1_l = 0.0
    end if

    if ( abs(u1_r*h_c) .LT. eps ) then
      u1_r = 0.0
    end if

    ! The edge slopes are limited from above by the respective
    ! one-sided slopes
    if ( abs(u1_l) .GT. abs(sigma_l) ) then
      u1_l = sigma_l
    end if

    if ( abs(u1_r) .GT. abs(sigma_r) ) then
      u1_r = sigma_r
    end if

    ! Build cubic interpolant (compute the coefficients)
    call build_cubic_interpolant( h, k, ppoly_E, ppoly_S, ppoly_coefficients )

    ! Check whether cubic is monotonic
    monotonic = is_cubic_monotonic( ppoly_coefficients, k )

    ! If cubic is not monotonic, monotonize it by modifiying the
    ! edge slopes, store the new edge slopes and recompute the
    ! cubic coefficients
    if ( monotonic .EQ. 0 ) then
      call monotonize_cubic( h_c, u0_l, u0_r, sigma_l, sigma_r, slope, u1_l, u1_r )
    end if

    ! Store edge slopes
    ppoly_S(k,1) = u1_l
    ppoly_S(k,2) = u1_r

    ! Recompute coefficients of cubic
    call build_cubic_interpolant( h, k, ppoly_E, ppoly_S, ppoly_coefficients )

  end do ! loop on cells

end subroutine P3M_limiter


!------------------------------------------------------------------------------
! p3m boundary extrapolation
! -----------------------------------------------------------------------------
subroutine P3M_boundary_extrapolation( N, h, u, ppoly_E, ppoly_S, ppoly_coefficients )
!------------------------------------------------------------------------------
! The following explanations apply to the left boundary cell. The same
! reasoning holds for the right boundary cell.
!
! A cubic needs to be built in the cell and requires four degrees of freedom,
! which are the edge values and slopes. The right edge values and slopes are
! taken to be that of the neighboring cell (i.e., the left edge value and slope
! of the neighboring cell). The left edge value and slope are determined by
! computing the parabola based on the cell average and the right edge value
! and slope. The resulting cubic is not necessarily monotonic and the slopes
! are subsequently modified to yield a monotonic cubic.
!------------------------------------------------------------------------------

  ! Arguments
  integer,              intent(in)    :: N ! Number of cells
  real, dimension(:),   intent(in)    :: h ! cell widths (size N)
  real, dimension(:),   intent(in)    :: u ! cell averages (size N)
  real, dimension(:,:), intent(inout) :: ppoly_E            !Edge value of polynomial
  real, dimension(:,:), intent(inout) :: ppoly_S            !Edge slope of polynomial
  real, dimension(:,:), intent(inout) :: ppoly_coefficients !Coefficients of polynomial

  ! Local variables
  integer       :: i0, i1
  integer       :: monotonic
  real          :: u0, u1
  real          :: h0, h1
  real          :: b, c, d
  real          :: u0_l, u0_r
  real          :: u1_l, u1_r
  real          :: eps
  real          :: slope

  eps = 1e-10

  ! ----- Left boundary -----
  i0 = 1
  i1 = 2
  h0 = h(i0) + eps
  h1 = h(i1) + eps
  u0 = u(i0)
  u1 = u(i1)

  ! Compute the left edge slope in neighboring cell and express it in
  ! the global coordinate system
  b = ppoly_coefficients(i1,2)
  u1_r = b / h1     ! derivative evaluated at xi = 0.0, expressed w.r.t. x

  ! Limit the right slope by the PLM limited slope
  slope = 2.0 * ( u1 - u0 ) / ( h0 + h_neglect )
  if ( abs(u1_r) .GT. abs(slope) ) then
    u1_r = slope
  end if

  ! The right edge value in the boundary cell is taken to be the left
  ! edge value in the neighboring cell
  u0_r = ppoly_E(i1,1)

  ! Given the right edge value and slope, we determine the left
  ! edge value and slope by computing the parabola as determined by
  ! the right edge value and slope and the boundary cell average
  u0_l = 3.0 * u0 + 0.5 * h0*u1_r - 2.0 * u0_r
  u1_l = ( - 6.0 * u0 - 2.0 * h0*u1_r + 6.0 * u0_r) / ( h0 + h_neglect )

  ! Check whether the edge values are monotonic. For example, if the left edge
  ! value is larger than the right edge value while the slope is positive, the
  ! edge values are inconsistent and we need to modify the left edge value
  if ( (u0_r-u0_l) * slope .LT. 0.0 ) then
    u0_l = u0_r
    u1_l = 0.0
    u1_r = 0.0
  end if

  ! Store edge values and slope, build cubic and check monotonicity
  ppoly_E(i0,1) = u0_l
  ppoly_E(i0,2) = u0_r
  ppoly_S(i0,1) = u1_l
  ppoly_S(i0,2) = u1_r

  ! Store edge values and slope, build cubic and check monotonicity
  call build_cubic_interpolant( h, i0, ppoly_E, ppoly_S, ppoly_coefficients )
  monotonic = is_cubic_monotonic( ppoly_coefficients, i0 )

  if ( monotonic .EQ. 0 ) then
    call monotonize_cubic( h0, u0_l, u0_r, 0.0, slope, slope, u1_l, u1_r )

    ! Rebuild cubic after monotonization
    ppoly_S(i0,1) = u1_l
    ppoly_S(i0,2) = u1_r
    call build_cubic_interpolant( h, i0, ppoly_E, ppoly_S, ppoly_coefficients )

  end if

  ! ----- Right boundary -----
  i0 = N-1
  i1 = N
  h0 = h(i0) + eps
  h1 = h(i1) + eps
  u0 = u(i0)
  u1 = u(i1)

  ! Compute the right edge slope in neighboring cell and express it in
  ! the global coordinate system
  b = ppoly_coefficients(i0,2)
  c = ppoly_coefficients(i0,3)
  d = ppoly_coefficients(i0,4)
  u1_l = (b + 2*c + 3*d) / ( h0 + h_neglect ) ! derivative evaluated at xi = 1.0

  ! Limit the left slope by the PLM limited slope
  slope = 2.0 * ( u1 - u0 ) / ( h1 + h_neglect )
  if ( abs(u1_l) .GT. abs(slope) ) then
    u1_l = slope
  end if

  ! The left edge value in the boundary cell is taken to be the right
  ! edge value in the neighboring cell
  u0_l = ppoly_E(i0,2)

  ! Given the left edge value and slope, we determine the right
  ! edge value and slope by computing the parabola as determined by
  ! the left edge value and slope and the boundary cell average
  u0_r = 3.0 * u1 - 0.5 * h1*u1_l - 2.0 * u0_l
  u1_r = ( 6.0 * u1 - 2.0 * h1*u1_l - 6.0 * u0_l) / ( h1 + h_neglect )

  ! Check whether the edge values are monotonic. For example, if the right edge
  ! value is smaller than the left edge value while the slope is positive, the
  ! edge values are inconsistent and we need to modify the right edge value
  if ( (u0_r-u0_l) * slope .LT. 0.0 ) then
    u0_r = u0_l
    u1_l = 0.0
    u1_r = 0.0
  end if

  ! Store edge values and slope, build cubic and check monotonicity
  ppoly_E(i1,1) = u0_l
  ppoly_E(i1,2) = u0_r
  ppoly_S(i1,1) = u1_l
  ppoly_S(i1,2) = u1_r

  call build_cubic_interpolant( h, i1, ppoly_E, ppoly_S, ppoly_coefficients )
  monotonic = is_cubic_monotonic( ppoly_coefficients, i1 )

  if ( monotonic .EQ. 0 ) then
    call monotonize_cubic( h1, u0_l, u0_r, slope, 0.0, slope, u1_l, u1_r )

    ! Rebuild cubic after monotonization
    ppoly_S(i1,1) = u1_l
    ppoly_S(i1,2) = u1_r
    call build_cubic_interpolant( h, i1, ppoly_E, ppoly_S, ppoly_coefficients )

  end if

end subroutine P3M_boundary_extrapolation


!------------------------------------------------------------------------------
! Build cubic interpolant in cell k
! -----------------------------------------------------------------------------
subroutine build_cubic_interpolant( h, k, ppoly_E, ppoly_S, ppoly_coefficients )
!------------------------------------------------------------------------------
! Given edge values and edge slopes, compute coefficients of cubic in cell k.
!
! NOTE: edge values and slopes MUST have been properly calculated prior to
! calling this routine.
!------------------------------------------------------------------------------

  ! Arguments
  real, dimension(:),   intent(in)    :: h ! cell widths (size N)
  integer,              intent(in)    :: k
  real, dimension(:,:), intent(in)    :: ppoly_E            !Edge value of polynomial
  real, dimension(:,:), intent(in)    :: ppoly_S            !Edge slope of polynomial
  real, dimension(:,:), intent(inout) :: ppoly_coefficients !Coefficients of polynomial

  ! Local variables
  real          :: u0_l, u0_r       ! edge values
  real          :: u1_l, u1_r       ! edge slopes
  real          :: h_c              ! cell width
  real          :: a0, a1, a2, a3   ! cubic coefficients

  h_c = h(k)

  u0_l = ppoly_E(k,1)
  u0_r = ppoly_E(k,2)

  u1_l = ppoly_S(k,1) * h_c
  u1_r = ppoly_S(k,2) * h_c

  a0 = u0_l
  a1 = u1_l
  a2 = 3.0 * ( u0_r - u0_l ) - u1_r - 2.0 * u1_l
  a3 = u1_r + u1_l + 2.0 * ( u0_l - u0_r )

  ppoly_coefficients(k,1) = a0
  ppoly_coefficients(k,2) = a1
  ppoly_coefficients(k,3) = a2
  ppoly_coefficients(k,4) = a3

end subroutine build_cubic_interpolant


!------------------------------------------------------------------------------
! Check whether cubic is monotonic
! -----------------------------------------------------------------------------
integer function is_cubic_monotonic( ppoly_coefficients, k )
!------------------------------------------------------------------------------
! This function checks whether the cubic curve in cell k is monotonic.
! If so, returns 1. Otherwise, returns 0.
!
! The cubic is monotonic if the first derivative is single-signed in [0,1].
! Hence, we check whether the roots (if any) lie inside this interval. If there
! is no root or if both roots lie outside this interval, the cubic is monotnic.
!------------------------------------------------------------------------------

  ! Arguments
  real, dimension(:,:), intent(in) :: ppoly_coefficients
  integer, intent(in)              :: k

  ! Local variables
  integer       :: monotonic        ! boolean indicating if monotonic or not
  real          :: a0, a1, a2, a3   ! cubic coefficients
  real          :: a, b, c          ! coefficients of first derivative
  real          :: xi_0, xi_1       ! roots of first derivative (if any !)
  real          :: rho
  real          :: eps

  ! Define the radius of the ball around 0 and 1 in which all values are assumed
  ! to be equal to 0 or 1, respectively
  eps = 1e-14

  a0 = ppoly_coefficients(k,1)
  a1 = ppoly_coefficients(k,2)
  a2 = ppoly_coefficients(k,3)
  a3 = ppoly_coefficients(k,4)

  a = a1
  b = 2.0 * a2
  c = 3.0 * a3

  xi_0 = -1.0
  xi_1 = -1.0

  rho = b*b - 4.0*a*c

  if ( rho .GE. 0.0 ) then
    if ( abs(c) .GT. 1e-15 ) then
      xi_0 = 0.5 * ( -b - sqrt( rho ) ) / c
      xi_1 = 0.5 * ( -b + sqrt( rho ) ) / c
    else if ( abs(b) .GT. 1e-15 ) then
      xi_0 = - a / b
      xi_1 = - a / b
    end if

    ! If one of the roots of the first derivative lies in (0,1),
    ! the cubic is not monotonic.
    if ( ( (xi_0 .GT. eps) .AND. (xi_0 .LT. 1.0-eps) ) .OR. &
         ( (xi_1 .GT. eps) .AND. (xi_1 .LT. 1.0-eps) ) ) then
      monotonic = 0
    else
      monotonic = 1
    end if

  else ! there are no real roots --> cubic is monotonic
    monotonic = 1
  end if

  ! Set the return value
  is_cubic_monotonic = monotonic

end function is_cubic_monotonic


!------------------------------------------------------------------------------
! Monotonize cubic curve
! -----------------------------------------------------------------------------
subroutine monotonize_cubic( h, u0_l, u0_r, sigma_l, sigma_r, slope, u1_l, u1_r )
!------------------------------------------------------------------------------
! This routine takes care of monotonizing a cubic on [0,1] by modifying the
! edge slopes. The edge values are NOT modified. The cubic is entirely
! determined by the four degrees of freedom u0_l, u0_r, u1_l and u1_r.
!
! u1_l and u1_r are the edge slopes expressed in the GLOBAL coordinate system.
!
! The monotonization occurs as follows.

! 1. The edge slopes are set to 0 if they are inconsistent with the limited
!    PLM slope
! 2. We check whether we can find an inflexion point in [0,1]. At most one
!    inflexion point may exist.
!    (a) If there is no inflexion point, the cubic is monotonic.
!    (b) If there is one inflexion point and it lies outside [0,1], the
!        cubic is monotonic.
!    (c) If there is one inflexion point and it lies in [0,1] and the slope
!        at the location of the inflexion point is consistent, the cubic
!        is monotonic.
!    (d) If the inflexion point lies in [0,1] but the slope is inconsistent,
!        we go to (3) to shift the location of the inflexion point to the left
!        or to the right. To the left when the 2nd-order left slope is smaller
!        than the 2nd order right slope.
! 3. Edge slopes are modified to shift the inflexion point, either onto the left
!    edge or onto the right edge.
!
!------------------------------------------------------------------------------

  ! Arguments
  real, intent(in)      :: h                ! cell width
  real, intent(in)      :: u0_l, u0_r       ! edge values
  real, intent(in)      :: sigma_l, sigma_r ! left and right 2nd-order slopes
  real, intent(in)      :: slope            ! limited PLM slope
  real, intent(inout)   :: u1_l, u1_r       ! edge slopes

  ! Local variables
  integer       :: found_ip
  integer       :: inflexion_l  ! bool telling if inflex. pt must be on left
  integer       :: inflexion_r  ! bool telling if inflex. pt must be on right
  real          :: eps
  real          :: a1, a2, a3
  real          :: u1_l_tmp     ! trial left edge slope
  real          :: u1_r_tmp     ! trial right edge slope
  real          :: xi_ip        ! location of inflexion point
  real          :: slope_ip     ! slope at inflexion point

  eps = 1e-14

  found_ip = 0
  inflexion_l = 0
  inflexion_r = 0

  ! If the edge slopes are inconsistent w.r.t. the limited PLM slope,
  ! set them to zero
  if ( u1_l*slope .LE. 0.0 ) then
    u1_l = 0.0
  end if

  if ( u1_r*slope .LE. 0.0 ) then
    u1_r = 0.0
  end if

  ! Compute the location of the inflexion point, which is the root
  ! of the second derivative
  a1 = h * u1_l
  a2 = 3.0 * ( u0_r - u0_l ) - h*(u1_r + 2.0*u1_l)
  a3 = h*(u1_r + u1_l) + 2.0*(u0_l - u0_r)

  ! There is a possible root (and inflexion point) only if a3 is nonzero.
  ! When a3 is zero, the second derivative of the cubic is constant (the
  ! cubic degenerates into a parabola) and no inflexion point exists.
  if ( a3 .NE. 0.0 ) then
    ! Location of inflexion point
    xi_ip = - a2 / (3.0 * a3)

    ! If the inflexion point lies in [0,1], change boolean value
    if ( (xi_ip .GE. 0.0) .AND. (xi_ip .LE. 1.0) ) then
      found_ip = 1
    end if
  end if

  ! When there is an inflexion point within [0,1], check the slope
  ! to see if it is consistent with the limited PLM slope. If not,
  ! decide on which side we want to collapse the inflexion point.
  ! If the inflexion point lies on one of the edges, the cubic is
  ! guaranteed to be monotonic
  if ( found_ip .EQ. 1 ) then
    slope_ip = a1 + 2.0*a2*xi_ip + 3.0*a3*xi_ip*xi_ip

    ! Check whether slope is consistent
    if ( slope_ip*slope .LT. 0.0 ) then
      if ( abs(sigma_l) .LT. abs(sigma_r)  ) then
        inflexion_l = 1
      else
        inflexion_r = 1
      end if
    end if
  end if ! found_ip

  ! At this point, if the cubic is not monotonic, we know where the
  ! inflexion point should lie. When the cubic is monotonic, both
  ! 'inflexion_l' and 'inflexion_r' are set to 0 and nothing is to be done.

  ! Move inflexion point on the left
  if ( inflexion_l .EQ. 1 ) then

    u1_l_tmp = 1.5*(u0_r-u0_l)/h - 0.5*u1_r
    u1_r_tmp = 3.0*(u0_r-u0_l)/h - 2.0*u1_l

    if ( (u1_l_tmp*slope .LT. 0.0) .AND. (u1_r_tmp*slope .LT. 0.0) ) then

      u1_l = 0.0
      u1_r = 3.0 * (u0_r - u0_l) / h

    else if (u1_l_tmp*slope .LT. 0.0) then

      u1_r = u1_r_tmp
      u1_l = 1.5*(u0_r - u0_l)/h - 0.5*u1_r

    else if (u1_r_tmp*slope .LT. 0.0) then

      u1_l = u1_l_tmp
      u1_r = 3.0*(u0_r - u0_l)/h - 2.0*u1_l

    else

      u1_l = u1_l_tmp
      u1_r = u1_r_tmp

    end if

  end if ! end treating case with inflexion point on the left

  ! Move inflexion point on the right
  if ( inflexion_r .EQ. 1 ) then

    u1_l_tmp = 3.0*(u0_r-u0_l)/h - 2.0*u1_r
    u1_r_tmp = 1.5*(u0_r-u0_l)/h - 0.5*u1_l

    if ( (u1_l_tmp*slope .LT. 0.0) .AND. (u1_r_tmp*slope .LT. 0.0) ) then

      u1_l = 3.0 * (u0_r - u0_l) / h
      u1_r = 0.0

    else if (u1_l_tmp*slope .LT. 0.0) then

      u1_r = u1_r_tmp
      u1_l = 3.0*(u0_r - u0_l)/h - 2.0*u1_r

    else if (u1_r_tmp*slope .LT. 0.0) then

      u1_l = u1_l_tmp
      u1_r = 1.5*(u0_r - u0_l)/h - 0.5*u1_l

    else

      u1_l = u1_l_tmp
      u1_r = u1_r_tmp

    end if

  end if ! end treating case with inflexion point on the right

  if ( abs(u1_l*h) .LT. eps ) then
    u1_l = 0.0
  end if

  if ( abs(u1_r*h) .LT. eps ) then
    u1_r = 0.0
  end if

end subroutine monotonize_cubic

end module P3M_functions
