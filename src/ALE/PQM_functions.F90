module PQM_functions
!==============================================================================
!
! This file is part of MOM.
!
! Date of creation: 2008.06.06
! L. White
!
! This module contains routines that handle one-dimensionnal finite volume
! reconstruction using the piecewise quartic method (PQM).
!
!==============================================================================
use regrid_edge_values, only : bound_edge_values, check_discontinuous_edge_values

implicit none ; private

public PQM_reconstruction, PQM_boundary_extrapolation, PQM_boundary_extrapolation_v1

real, parameter :: h_neglect = 1.E-30

contains

!------------------------------------------------------------------------------
! PQM_reconstruction
! -----------------------------------------------------------------------------
subroutine PQM_reconstruction( N, h, u, ppoly_E, ppoly_S, ppoly_coefficients )
!------------------------------------------------------------------------------
! Reconstruction by quartic polynomials within each cell.
!
! grid:  one-dimensional grid (see grid.F90)
! ppoly: piecewise quartic polynomial to be reconstructed (see ppoly.F90)
! u:     cell averages
!
! It is assumed that the dimension of 'u' is equal to the number of cells
! defining 'grid' and 'ppoly'. No consistency check is performed.
!------------------------------------------------------------------------------

  ! Arguments
  integer,              intent(in)    :: N ! Number of cells
  real, dimension(:),   intent(in)    :: h ! cell widths (size N)
  real, dimension(:),   intent(in)    :: u ! cell averages (size N)
  real, dimension(:,:), intent(inout) :: ppoly_E            !Edge value of polynomial
  real, dimension(:,:), intent(inout) :: ppoly_S            !Edge slope of polynomial
  real, dimension(:,:), intent(inout) :: ppoly_coefficients !Coefficients of polynomial

  ! Local variables
  integer   :: k                ! loop index
  real      :: h_c              ! cell width
  real      :: u0_l, u0_r       ! edge values (left and right)
  real      :: u1_l, u1_r       ! edge slopes (left and right)
  real      :: a, b, c, d, e    ! parabola coefficients

  ! PQM limiter
  call PQM_limiter( N, h, u, ppoly_E, ppoly_S )

  ! Loop on cells to construct the cubic within each cell
  do k = 1,N

    u0_l = ppoly_E(k,1)
    u0_r = ppoly_E(k,2)

    u1_l = ppoly_S(k,1)
    u1_r = ppoly_S(k,2)

    h_c = h(k)

    a = u0_l
    b = h_c * u1_l
    c = 30.0 * u(k) - 12.0*u0_r - 18.0*u0_l + 1.5*h_c*(u1_r - 3.0*u1_l)
    d = -60.0 * u(k) + h_c *(6.0*u1_l - 4.0*u1_r) + 28.0*u0_r + 32.0*u0_l
    e = 30.0 * u(k) + 2.5*h_c*(u1_r - u1_l) - 15.0*(u0_l + u0_r)

    ! Store coefficients
    ppoly_coefficients(k,1) = a
    ppoly_coefficients(k,2) = b
    ppoly_coefficients(k,3) = c
    ppoly_coefficients(k,4) = d
    ppoly_coefficients(k,5) = e

  end do ! end loop on cells

end subroutine PQM_reconstruction


!------------------------------------------------------------------------------
! Limit pqm
! -----------------------------------------------------------------------------
subroutine PQM_limiter( N, h, u, ppoly_E, ppoly_S )
!------------------------------------------------------------------------------
! Standard PQM limiter (White & Adcroft, JCP 2008).
!
! grid:  one-dimensional grid (see grid.F90)
! ppoly: piecewise quadratic polynomial to be reconstructed (see ppoly.F90)
! u:     cell averages
!
! It is assumed that the dimension of 'u' is equal to the number of cells
! defining 'grid' and 'ppoly'. No consistency check is performed.
!------------------------------------------------------------------------------

  ! Arguments
  integer,              intent(in)    :: N ! Number of cells
  real, dimension(:),   intent(in)    :: h ! cell widths (size N)
  real, dimension(:),   intent(in)    :: u ! cell averages (size N)
  real, dimension(:,:), intent(inout) :: ppoly_E            !Edge value of polynomial
  real, dimension(:,:), intent(inout) :: ppoly_S            !Edge slope of polynomial

  ! Local variables
  integer   :: k            ! loop index
  integer   :: inflexion_l
  integer   :: inflexion_r
  real      :: u0_l, u0_r   ! edge values
  real      :: u1_l, u1_r   ! edge slopes
  real      :: u_l, u_c, u_r        ! left, center and right cell averages
  real      :: h_l, h_c, h_r        ! left, center and right cell widths
  real      :: sigma_l, sigma_c, sigma_r    ! left, center and right
                                            ! van Leer slopes
  real      :: slope        ! retained PLM slope
  real      :: a, b, c, d, e
  real      :: alpha1, alpha2, alpha3
  real      :: rho, sqrt_rho
  real      :: gradient1, gradient2
  real      :: x1, x2

  ! Bound edge values
  call bound_edge_values( N, h, u, ppoly_E )

  ! Make discontinuous edge values monotonic (thru averaging)
  call check_discontinuous_edge_values( N, u, ppoly_E )

  ! Loop on interior cells to apply the PQM limiter
  do k = 2,N-1

    !if ( h(k) .lt. 1.0 ) cycle

    inflexion_l = 0
    inflexion_r = 0

    ! Get edge values, edge slopes and cell width
    u0_l = ppoly_E(k,1)
    u0_r = ppoly_E(k,2)
    u1_l = ppoly_S(k,1)
    u1_r = ppoly_S(k,2)

    ! Get cell widths and cell averages (boundary cells are assumed to
    ! be local extrema for the sake of slopes)
    h_l = h(k-1)
    h_c = h(k)
    h_r = h(k+1)
    u_l = u(k-1)
    u_c = u(k)
    u_r = u(k+1)

    ! Compute limited slope
    sigma_l = 2.0 * ( u_c - u_l ) / ( h_c + h_neglect )
    sigma_c = 2.0 * ( u_r - u_l ) / ( h_l + 2.0*h_c + h_r + h_neglect )
    sigma_r = 2.0 * ( u_r - u_c ) / ( h_c + h_neglect )

    if ( (sigma_l * sigma_r) .GT. 0.0 ) then
      slope = sign( min(abs(sigma_l),abs(sigma_c),abs(sigma_r)), sigma_c )
    else
      slope = 0.0
    end if

    ! If one of the slopes has the wrong sign compared with the
    ! limited PLM slope, it is set equal to the limited PLM slope
    if ( u1_l*slope .le. 0.0 ) u1_l = slope
    if ( u1_r*slope .le. 0.0 ) u1_r = slope

    ! Local extremum --> flatten
    if ( (u0_r - u_c) * (u_c - u0_l) .le. 0.0) then
      u0_l = u_c
      u0_r = u_c
      u1_l = 0.0
      u1_r = 0.0
      inflexion_l = -1
      inflexion_r = -1
    end if

    ! Edge values are bounded and averaged when discontinuous and not
    ! monotonic, edge slopes are consistent and the cell is not an extremum.
    ! We now need to check and encorce the monotonicity of the quartic within
    ! the cell
    if ( (inflexion_l .EQ. 0) .AND. (inflexion_r .EQ. 0) ) then

      a = u0_l
      b = h_c * u1_l
      c = 30.0 * u(k) - 12.0*u0_r - 18.0*u0_l + 1.5*h_c*(u1_r - 3.0*u1_l)
      d = -60.0 * u(k) + h_c *(6.0*u1_l - 4.0*u1_r) + 28.0*u0_r + 32.0*u0_l
      e = 30.0 * u(k) + 2.5*h_c*(u1_r - u1_l) - 15.0*(u0_l + u0_r)

      ! Determine the coefficients of the second derivative
      ! alpha1 xi^2 + alpha2 xi + alpha3
      alpha1 = 6*e
      alpha2 = 3*d
      alpha3 = c

      rho = alpha2 * alpha2 - 4.0 * alpha1 * alpha3

      ! Check whether inflexion points exist
      if (( alpha1 .ne. 0.0 ) .and. ( rho .ge. 0.0 )) then

        sqrt_rho = sqrt( rho )

        x1 = 0.5 * ( - alpha2 - sqrt_rho ) / alpha1
        x2 = 0.5 * ( - alpha2 + sqrt_rho ) / alpha1

        ! Check whether both inflexion points lie in [0,1]
        if ( (x1 .GE. 0.0) .AND. (x1 .LE. 1.0) .AND. &
             (x2 .GE. 0.0) .AND. (x2 .LE. 1.0) ) then

          gradient1 = 4.0 * e * (x1**3) + 3.0 * d * (x1**2) + 2.0 * c * x1 + b
          gradient2 = 4.0 * e * (x2**3) + 3.0 * d * (x2**2) + 2.0 * c * x2 + b

          ! Check whether one of the gradients is inconsistent
          if ( (gradient1 * slope .LT. 0.0) .OR. &
               (gradient2 * slope .LT. 0.0) ) then
            ! Decide where to collapse inflexion points
            ! (depends on one-sided slopes)
            if ( abs(sigma_l) .LT. abs(sigma_r) ) then
              inflexion_l = 1
            else
              inflexion_r = 1
            end if
          end if

        ! If both x1 and x2 do not lie in [0,1], check whether
        ! only x1 lies in [0,1]
        else if ( (x1 .GE. 0.0) .AND. (x1 .LE. 1.0) ) then

          gradient1 = 4.0 * e * (x1**3) + 3.0 * d * (x1**2) + 2.0 * c * x1 + b

          ! Check whether the gradient is inconsistent
          if ( gradient1 * slope .LT. 0.0 ) then
            ! Decide where to collapse inflexion points
            ! (depends on one-sided slopes)
            if ( abs(sigma_l) .LT. abs(sigma_r) ) then
              inflexion_l = 1
            else
              inflexion_r = 1
            end if
          end if

        ! If x1 does not lie in [0,1], check whether x2 lies in [0,1]
        else if ( (x2 .GE. 0.0) .AND. (x2 .LE. 1.0) ) then

          gradient2 = 4.0 * e * (x2**3) + 3.0 * d * (x2**2) + 2.0 * c * x2 + b

          ! Check whether the gradient is inconsistent
          if ( gradient2 * slope .LT. 0.0 ) then
            ! Decide where to collapse inflexion points
            ! (depends on one-sided slopes)
            if ( abs(sigma_l) .LT. abs(sigma_r) ) then
              inflexion_l = 1
            else
              inflexion_r = 1
            end if
          end if

        end if ! end checking where the inflexion points lie

      end if ! end checking if alpha1 != 0 AND rho >= 0

      ! If alpha1 is zero, the second derivative of the quartic reduces
      ! to a straight line
      if (( alpha1 .eq. 0.0 ) .and. ( alpha2 .ne. 0.0 )) then

          x1 = - alpha3 / alpha2
          if ( (x1 .ge. 0.0) .AND. (x1 .le. 1.0) ) then

            gradient1 = 4.0 * e * (x1**3) + 3.0 * d * (x1**2) + 2.0 * c * x1 + b

            ! Check whether the gradient is inconsistent
            if ( gradient1 * slope .LT. 0.0 ) then
              ! Decide where to collapse inflexion points
              ! (depends on one-sided slopes)
              if ( abs(sigma_l) .LT. abs(sigma_r) ) then
                inflexion_l = 1
              else
                inflexion_r = 1
              end if
            end if ! check slope consistency

          end if

      end if ! end check whether we can find the root of the straight line

    end if ! end checking whether to shift inflexion points

    ! At this point, we know onto which edge to shift inflexion points
    if ( inflexion_l .EQ. 1 ) then

      ! We modify the edge slopes so that both inflexion points
      ! collapse onto the left edge
      u1_l = ( 10.0 * u_c - 2.0 * u0_r - 8.0 * u0_l ) / (3.0*h_c + h_neglect )
      u1_r = ( -10.0 * u_c + 6.0 * u0_r + 4.0 * u0_l ) / ( h_c + h_neglect )

      ! One of the modified slopes might be inconsistent. When that happens,
      ! the inconsistent slope is set equal to zero and the opposite edge value
      ! and edge slope are modified in compliance with the fact that both
      ! inflexion points must still be located on the left edge
      if ( u1_l * slope .LT. 0.0 ) then

        u1_l = 0.0
        u0_r = 5.0 * u_c - 4.0 * u0_l
        u1_r = 20.0 * (u_c - u0_l) / ( h_c + h_neglect )

      else if ( u1_r * slope .LT. 0.0 ) then

        u1_r = 0.0
        u0_l = (5.0*u_c - 3.0*u0_r) / 2.0
        u1_l = 10.0 * (-u_c + u0_r) / (3.0 * h_c + h_neglect)

      end if

    else if ( inflexion_r .EQ. 1 ) then

      ! We modify the edge slopes so that both inflexion points
      ! collapse onto the right edge
      u1_r = ( -10.0 * u_c + 8.0 * u0_r + 2.0 * u0_l ) / (3.0 * h_c + h_neglect)
      u1_l = ( 10.0 * u_c - 4.0 * u0_r - 6.0 * u0_l ) / (h_c + h_neglect)

      ! One of the modified slopes might be inconsistent. When that happens,
      ! the inconsistent slope is set equal to zero and the opposite edge value
      ! and edge slope are modified in compliance with the fact that both
      ! inflexion points must still be located on the right edge
      if ( u1_l * slope .LT. 0.0 ) then

        u1_l = 0.0
        u0_r = ( 5.0 * u_c - 3.0 * u0_l ) / 2.0
        u1_r = 10.0 * (u_c - u0_l) / (3.0 * h_c + h_neglect)

      else if ( u1_r * slope .LT. 0.0 ) then

        u1_r = 0.0
        u0_l = 5.0 * u_c - 4.0 * u0_r
        u1_l = 20.0 * ( -u_c + u0_r ) / (h_c + h_neglect)

      end if

    end if ! clause to check where to collapse inflexion points

    ! Save edge values and edge slopes for reconstruction
    ppoly_E(k,1) = u0_l
    ppoly_E(k,2) = u0_r
    ppoly_S(k,1) = u1_l
    ppoly_S(k,2) = u1_r

  end do ! end loop on interior cells

  ! Constant reconstruction within boundary cells
  ppoly_E(1,:) = u(1)
  ppoly_S(1,:) = 0.0

  ppoly_E(N,:) = u(N)
  ppoly_S(N,:) = 0.0

end subroutine PQM_limiter


!------------------------------------------------------------------------------
! pqm boundary extrapolation
! -----------------------------------------------------------------------------
subroutine PQM_boundary_extrapolation( N, h, u, ppoly_E, ppoly_coefficients )
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
  integer,              intent(in)    :: N ! Number of cells
  real, dimension(:),   intent(in)    :: h ! cell widths (size N)
  real, dimension(:),   intent(in)    :: u ! cell averages (size N)
  real, dimension(:,:), intent(inout) :: ppoly_E            !Edge value of polynomial
  real, dimension(:,:), intent(inout) :: ppoly_coefficients !Coefficients of polynomial

  ! Local variables
  integer       :: i0, i1
  real          :: u0, u1
  real          :: h0, h1
  real          :: a, b, c, d, e
  real          :: u0_l, u0_r
  real          :: u1_l, u1_r
  real          :: slope
  real          :: exp1, exp2

  ! ----- Left boundary -----
  i0 = 1
  i1 = 2
  h0 = h(i0)
  h1 = h(i1)
  u0 = u(i0)
  u1 = u(i1)

  ! Compute the left edge slope in neighboring cell and express it in
  ! the global coordinate system
  b = ppoly_coefficients(i1,2)
  u1_r = b *(h0/h1)     ! derivative evaluated at xi = 0.0,
                        ! expressed w.r.t. xi (local coord. system)

  ! Limit the right slope by the PLM limited slope
  slope = 2.0 * ( u1 - u0 )
  if ( abs(u1_r) .GT. abs(slope) ) then
    u1_r = slope
  end if

  ! The right edge value in the boundary cell is taken to be the left
  ! edge value in the neighboring cell
  u0_r = ppoly_E(i1,1)

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

  ppoly_E(i0,1) = u0_l
  ppoly_E(i0,2) = u0_r

  a = u0_l
  b = 6.0 * u0 - 4.0 * u0_l - 2.0 * u0_r
  c = 3.0 * ( u0_r + u0_l - 2.0 * u0 )

  ! The quartic is reduced to a parabola in the boundary cell
  ppoly_coefficients(i0,1) = a
  ppoly_coefficients(i0,2) = b
  ppoly_coefficients(i0,3) = c
  ppoly_coefficients(i0,4) = 0.0
  ppoly_coefficients(i0,5) = 0.0

  ! ----- Right boundary -----
  i0 = N-1
  i1 = N
  h0 = h(i0)
  h1 = h(i1)
  u0 = u(i0)
  u1 = u(i1)

  ! Compute the right edge slope in neighboring cell and express it in
  ! the global coordinate system
  b = ppoly_coefficients(i0,2)
  c = ppoly_coefficients(i0,3)
  d = ppoly_coefficients(i0,4)
  e = ppoly_coefficients(i0,5)
  u1_l = (b + 2*c + 3*d + 4*e)      ! derivative evaluated at xi = 1.0
  u1_l = u1_l * (h1/h0)

  ! Limit the left slope by the PLM limited slope
  slope = 2.0 * ( u1 - u0 )
  if ( abs(u1_l) .GT. abs(slope) ) then
    u1_l = slope
  end if

  ! The left edge value in the boundary cell is taken to be the right
  ! edge value in the neighboring cell
  u0_l = ppoly_E(i0,2)

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

  ppoly_E(i1,1) = u0_l
  ppoly_E(i1,2) = u0_r

  a = u0_l
  b = 6.0 * u1 - 4.0 * u0_l - 2.0 * u0_r
  c = 3.0 * ( u0_r + u0_l - 2.0 * u1 )

  ! The quartic is reduced to a parabola in the boundary cell
  ppoly_coefficients(i1,1) = a
  ppoly_coefficients(i1,2) = b
  ppoly_coefficients(i1,3) = c
  ppoly_coefficients(i1,4) = 0.0
  ppoly_coefficients(i1,5) = 0.0

end subroutine PQM_boundary_extrapolation


!------------------------------------------------------------------------------
! pqm boundary extrapolation using rational function
! -----------------------------------------------------------------------------
subroutine PQM_boundary_extrapolation_v1( N, h, u, ppoly_E, ppoly_S, ppoly_coefficients )
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
  integer,              intent(in)    :: N ! Number of cells
  real, dimension(:),   intent(in)    :: h ! cell widths (size N)
  real, dimension(:),   intent(in)    :: u ! cell averages (size N)
  real, dimension(:,:), intent(inout) :: ppoly_E            !Edge value of polynomial
  real, dimension(:,:), intent(inout) :: ppoly_S            !Edge slope of polynomial
  real, dimension(:,:), intent(inout) :: ppoly_coefficients !Coefficients of polynomial

  ! Local variables
  integer       :: i0, i1
  integer       :: inflexion_l
  integer       :: inflexion_r
  real          :: u0, u1, um
  real          :: h0, h1
  real          :: a, b, c, d, e
  real          :: ar, br, beta
  real          :: u0_l, u0_r
  real          :: u1_l, u1_r
  real          :: u_plm
  real          :: slope
  real          :: alpha1, alpha2, alpha3
  real          :: rho, sqrt_rho
  real          :: gradient1, gradient2
  real          :: x1, x2

  ! ----- Left boundary (TOP) -----
  i0 = 1
  i1 = 2
  h0 = h(i0)
  h1 = h(i1)
  u0 = u(i0)
  u1 = u(i1)
  um = u0

  ! Compute real slope and express it w.r.t. local coordinate system
  ! within boundary cell
  slope = 2.0 * ( u1 - u0 ) / ( ( h0 + h1 ) + h_neglect )
  slope = slope * h0

  ! The right edge value and slope of the boundary cell are taken to be the
  ! left edge value and slope of the adjacent cell
  a = ppoly_coefficients(i1,1)
  b = ppoly_coefficients(i1,2)

  u0_r = a          ! edge value
  u1_r = b / (h1 + h_neglect) ! edge slope (w.r.t. global coord.)

  ! Compute coefficient for rational function based on mean and right
  ! edge value and slope
  if (u1_r.ne.0.) then ! HACK by AJA
    beta = 2.0 * ( u0_r - um ) / ( (h0 + h_neglect)*u1_r) - 1.0
  else
    beta = 0.
  endif ! HACK by AJA
  br = u0_r + beta*u0_r - um
  ar = um + beta*um - br

  ! Left edge value estimate based on rational function
  u0_l = ar

  ! Edge value estimate based on PLM
  u_plm = um - 0.5 * slope

  ! Check whether the left edge value is bounded by the mean and
  ! the PLM edge value. If so, keep it and compute left edge slope
  ! based on the rational function. If not, keep the PLM edge value and
  ! compute corresponding slope.
  if ( abs(um-u0_l) .lt. abs(um-u_plm) ) then
    u1_l = 2.0 * ( br - ar*beta)
    u1_l = u1_l / (h0 + h_neglect)
  else
    u0_l = u_plm
    u1_l = slope / (h0 + h_neglect)
  end if

  ! Monotonize quartic
  inflexion_l = 0

  a = u0_l
  b = h0 * u1_l
  c = 30.0 * um - 12.0*u0_r - 18.0*u0_l + 1.5*h0*(u1_r - 3.0*u1_l)
  d = -60.0 * um + h0 *(6.0*u1_l - 4.0*u1_r) + 28.0*u0_r + 32.0*u0_l
  e = 30.0 * um + 2.5*h0*(u1_r - u1_l) - 15.0*(u0_l + u0_r)

  alpha1 = 6*e
  alpha2 = 3*d
  alpha3 = c

  rho = alpha2 * alpha2 - 4.0 * alpha1 * alpha3

  ! Check whether inflexion points exist. If so, transform the quartic
  ! so that both inflexion points coalesce on the left edge.
  if (( alpha1 .ne. 0.0 ) .and. ( rho .ge. 0.0 )) then

    sqrt_rho = sqrt( rho )

    x1 = 0.5 * ( - alpha2 - sqrt_rho ) / alpha1
    if ( (x1 .gt. 0.0) .and. (x1 .lt. 1.0) ) then
      gradient1 = 4.0 * e * (x1**3) + 3.0 * d * (x1**2) + 2.0 * c * x1 + b
      if ( gradient1 * slope .lt. 0.0 ) then
        inflexion_l = 1
      end if
    end if

    x2 = 0.5 * ( - alpha2 + sqrt_rho ) / alpha1
    if ( (x2 .gt. 0.0) .and. (x2 .lt. 1.0) ) then
      gradient2 = 4.0 * e * (x2**3) + 3.0 * d * (x2**2) + 2.0 * c * x2 + b
      if ( gradient2 * slope .lt. 0.0 ) then
        inflexion_l = 1
      end if
    end if

  end if

  if (( alpha1 .eq. 0.0 ) .and. ( alpha2 .ne. 0.0 )) then

    x1 = - alpha3 / alpha2
    if ( (x1 .ge. 0.0) .and. (x1 .le. 1.0) ) then
      gradient1 = 3.0 * d * (x1**2) + 2.0 * c * x1 + b
      if ( gradient1 * slope .lt. 0.0 ) then
        inflexion_l = 1
      end if
    end if

  end if

  if ( inflexion_l .eq. 1 ) then

    ! We modify the edge slopes so that both inflexion points
    ! collapse onto the left edge
    u1_l = ( 10.0 * um - 2.0 * u0_r - 8.0 * u0_l ) / (3.0*h0 + h_neglect)
    u1_r = ( -10.0 * um + 6.0 * u0_r + 4.0 * u0_l ) / (h0 + h_neglect)

    ! One of the modified slopes might be inconsistent. When that happens,
    ! the inconsistent slope is set equal to zero and the opposite edge value
    ! and edge slope are modified in compliance with the fact that both
    ! inflexion points must still be located on the left edge
    if ( u1_l * slope .LT. 0.0 ) then

      u1_l = 0.0
      u0_r = 5.0 * um - 4.0 * u0_l
      u1_r = 20.0 * (um - u0_l) / ( h0 + h_neglect )

    else if ( u1_r * slope .LT. 0.0 ) then

      u1_r = 0.0
      u0_l = (5.0*um - 3.0*u0_r) / 2.0
      u1_l = 10.0 * (-um + u0_r) / (3.0 * h0 + h_neglect )

    end if

  end if

  ! Store edge values, edge slopes and coefficients
  ppoly_E(i0,1) = u0_l
  ppoly_E(i0,2) = u0_r
  ppoly_S(i0,1) = u1_l
  ppoly_S(i0,2) = u1_r

  a = u0_l
  b = h0 * u1_l
  c = 30.0 * um - 12.0*u0_r - 18.0*u0_l + 1.5*h0*(u1_r - 3.0*u1_l)
  d = -60.0 * um + h0 *(6.0*u1_l - 4.0*u1_r) + 28.0*u0_r + 32.0*u0_l
  e = 30.0 * um + 2.5*h0*(u1_r - u1_l) - 15.0*(u0_l + u0_r)

    ! Store coefficients
  ppoly_coefficients(i0,1) = a
  ppoly_coefficients(i0,2) = b
  ppoly_coefficients(i0,3) = c
  ppoly_coefficients(i0,4) = d
  ppoly_coefficients(i0,5) = e

  ! ----- Right boundary (BOTTOM) -----
  i0 = N-1
  i1 = N
  h0 = h(i0)
  h1 = h(i1)
  u0 = u(i0)
  u1 = u(i1)
  um = u1

  ! Compute real slope and express it w.r.t. local coordinate system
  ! within boundary cell
  slope = 2.0 * ( u1 - u0 ) / ( h0 + h1 )
  slope = slope * h1

  ! The left edge value and slope of the boundary cell are taken to be the
  ! right edge value and slope of the adjacent cell
  a = ppoly_coefficients(i0,1)
  b = ppoly_coefficients(i0,2)
  c = ppoly_coefficients(i0,3)
  d = ppoly_coefficients(i0,4)
  e = ppoly_coefficients(i0,5)
  u0_l = a + b + c + d + e                  ! edge value
  u1_l = (b + 2*c + 3*d + 4*e) / h0         ! edge slope (w.r.t. global coord.)

  ! Compute coefficient for rational function based on mean and left
  ! edge value and slope
  if (um-u0_l.ne.0.) then ! HACK by AJA
    beta = 0.5*h1*u1_l / (um-u0_l) - 1.0
  else
    beta = 0.
  endif ! HACK by AJA
  br = beta*um + um - u0_l
  ar = u0_l

  ! Right edge value estimate based on rational function
  if (1+beta.ne.0.) then ! HACK by AJA
    u0_r = (ar + 2*br + beta*br ) / ((1+beta)*(1+beta))
  else
    u0_r = um + 0.5 * slope ! PLM
  endif ! HACK by AJA

  ! Right edge value estimate based on PLM
  u_plm = um + 0.5 * slope

  ! Check whether the right edge value is bounded by the mean and
  ! the PLM edge value. If so, keep it and compute right edge slope
  ! based on the rational function. If not, keep the PLM edge value and
  ! compute corresponding slope.
  if ( abs(um-u0_r) .lt. abs(um-u_plm) ) then
    u1_r = 2.0 * ( br - ar*beta ) / ( (1+beta)*(1+beta)*(1+beta) )
    u1_r = u1_r / h1
  else
    u0_r = u_plm
    u1_r = slope / h1
  end if

  ! Monotonize quartic
  inflexion_r = 0

  a = u0_l
  b = h1 * u1_l
  c = 30.0 * um - 12.0*u0_r - 18.0*u0_l + 1.5*h1*(u1_r - 3.0*u1_l)
  d = -60.0 * um + h1*(6.0*u1_l - 4.0*u1_r) + 28.0*u0_r + 32.0*u0_l
  e = 30.0 * um + 2.5*h1*(u1_r - u1_l) - 15.0*(u0_l + u0_r)

  alpha1 = 6*e
  alpha2 = 3*d
  alpha3 = c

  rho = alpha2 * alpha2 - 4.0 * alpha1 * alpha3

  ! Check whether inflexion points exist. If so, transform the quartic
  ! so that both inflexion points coalesce on the right edge.
  if (( alpha1 .ne. 0.0 ) .and. ( rho .ge. 0.0 )) then

    sqrt_rho = sqrt( rho )

    x1 = 0.5 * ( - alpha2 - sqrt_rho ) / alpha1
    if ( (x1 .gt. 0.0) .and. (x1 .lt. 1.0) ) then
      gradient1 = 4.0 * e * (x1**3) + 3.0 * d * (x1**2) + 2.0 * c * x1 + b
      if ( gradient1 * slope .lt. 0.0 ) then
        inflexion_r = 1
      end if
    end if

    x2 = 0.5 * ( - alpha2 + sqrt_rho ) / alpha1
    if ( (x2 .gt. 0.0) .and. (x2 .lt. 1.0) ) then
      gradient2 = 4.0 * e * (x2**3) + 3.0 * d * (x2**2) + 2.0 * c * x2 + b
      if ( gradient2 * slope .lt. 0.0 ) then
        inflexion_r = 1
      end if
    end if

  end if

  if (( alpha1 .eq. 0.0 ) .and. ( alpha2 .ne. 0.0 )) then

    x1 = - alpha3 / alpha2
    if ( (x1 .ge. 0.0) .and. (x1 .le. 1.0) ) then
      gradient1 = 3.0 * d * (x1**2) + 2.0 * c * x1 + b
      if ( gradient1 * slope .lt. 0.0 ) then
        inflexion_r = 1
      end if
    end if

  end if

  if ( inflexion_r .eq. 1 ) then

    ! We modify the edge slopes so that both inflexion points
    ! collapse onto the right edge
    u1_r = ( -10.0 * um + 8.0 * u0_r + 2.0 * u0_l ) / (3.0 * h1)
    u1_l = ( 10.0 * um - 4.0 * u0_r - 6.0 * u0_l ) / h1

    ! One of the modified slopes might be inconsistent. When that happens,
    ! the inconsistent slope is set equal to zero and the opposite edge value
    ! and edge slope are modified in compliance with the fact that both
    ! inflexion points must still be located on the right edge
    if ( u1_l * slope .lt. 0.0 ) then

      u1_l = 0.0
      u0_r = ( 5.0 * um - 3.0 * u0_l ) / 2.0
      u1_r = 10.0 * (um - u0_l) / (3.0 * h1)

    else if ( u1_r * slope .lt. 0.0 ) then

      u1_r = 0.0
      u0_l = 5.0 * um - 4.0 * u0_r
      u1_l = 20.0 * ( -um + u0_r ) / h1

    end if

  end if

  ! Store edge values, edge slopes and coefficients
  ppoly_E(i1,1) = u0_l
  ppoly_E(i1,2) = u0_r
  ppoly_S(i1,1) = u1_l
  ppoly_S(i1,2) = u1_r

  a = u0_l
  b = h1 * u1_l
  c = 30.0 * um - 12.0*u0_r - 18.0*u0_l + 1.5*h1*(u1_r - 3.0*u1_l)
  d = -60.0 * um + h1 *(6.0*u1_l - 4.0*u1_r) + 28.0*u0_r + 32.0*u0_l
  e = 30.0 * um + 2.5*h1*(u1_r - u1_l) - 15.0*(u0_l + u0_r)

  ppoly_coefficients(i1,1) = a
  ppoly_coefficients(i1,2) = b
  ppoly_coefficients(i1,3) = c
  ppoly_coefficients(i1,4) = d
  ppoly_coefficients(i1,5) = e

end subroutine PQM_boundary_extrapolation_v1

end module PQM_functions
