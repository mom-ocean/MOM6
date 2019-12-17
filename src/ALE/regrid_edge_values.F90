!> Edge value estimation for high-order resconstruction
module regrid_edge_values

! This file is part of MOM6. See LICENSE.md for the license.

use regrid_solvers, only : solve_linear_system, solve_tridiagonal_system, solve_diag_dominant_tridiag
use polynomial_functions, only : evaluation_polynomial

implicit none ; private

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

! The following parameters are used to avoid singular matrices for boundary
! extrapolation. The are needed only in the case where thicknesses vanish
! to a small enough values such that the eigenvalues of the matrix can not
! be separated.
!   Specifying a dimensional parameter value, as is done here, is a terrible idea.
real, parameter :: hNeglect_edge_dflt = 1.e-10 !< The default value for cut-off minimum
                                          !! thickness for sum(h) in edge value inversions
real, parameter :: hNeglect_dflt = 1.e-30 !< The default value for cut-off minimum
                                          !! thickness for sum(h) in other calculations
real, parameter :: hMinFrac      = 1.e-5  !< A minimum fraction for min(h)/sum(h)

contains

!> Bound edge values by neighboring cell averages
!!
!! In this routine, we loop on all cells to bound their left and right
!! edge values by the cell averages. That is, the left edge value must lie
!! between the left cell average and the central cell average. A similar
!! reasoning applies to the right edge values.
!!
!! Both boundary edge values are set equal to the boundary cell averages.
!! Any extrapolation scheme is applied after this routine has been called.
!! Therefore, boundary cells are treated as if they were local extrama.
subroutine bound_edge_values( N, h, u, edge_val, h_neglect, answers_2018 )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< cell widths [H]
  real, dimension(N),   intent(in)    :: u !< cell average properties in arbitrary units [A]
  real, dimension(N,2), intent(inout) :: edge_val !< Potentially modified edge values [A]; the
                                           !! second index is for the two edges of each cell.
  real,       optional, intent(in)    :: h_neglect !< A negligibly small width [H]
  logical,    optional, intent(in)    :: answers_2018 !< If true use older, less acccurate expressions.
  ! Local variables
  real    :: sigma_l, sigma_c, sigma_r    ! left, center and right van Leer slopes [A H-1] or [A]
  real    :: slope_x_h     ! retained PLM slope times  half grid step [A]
  real    :: hNeglect      ! A negligible thickness [H].
  logical :: use_2018_answers  ! If true use older, less acccurate expressions.
  integer :: k, km1, kp1   ! Loop index and the values to either side.

  use_2018_answers = .true. ; if (present(answers_2018)) use_2018_answers = answers_2018
  if (use_2018_answers) then
    hNeglect = hNeglect_dflt ; if (present(h_neglect)) hNeglect = h_neglect
  endif

  ! Loop on cells to bound edge value
  do k = 1,N

    ! For the sake of bounding boundary edge values, the left neighbor of the left boundary cell
    ! is assumed to be the same as the left boundary cell and the right neighbor of the right
    ! boundary cell is assumed to be the same as the right boundary cell. This effectively makes
    ! boundary cells look like extrema.
    km1 = max(1,k-1) ; kp1 = min(k+1,N)

    slope_x_h = 0.0
    if (use_2018_answers) then
      sigma_l = 2.0 * ( u(k) - u(km1) ) / ( h(k) + hNeglect )
      sigma_c = 2.0 * ( u(kp1) - u(km1) ) / ( h(km1) + 2.0*h(k) + h(kp1) + hNeglect )
      sigma_r = 2.0 * ( u(kp1) - u(k) ) / ( h(k) + hNeglect )

      ! The limiter is used in the local coordinate system to each cell, so for convenience store
      ! the slope times a half grid spacing.  (See White and Adcroft JCP 2008 Eqs 19 and 20)
      if ( (sigma_l * sigma_r) > 0.0 ) &
        slope_x_h = 0.5 * h(k) * sign( min(abs(sigma_l),abs(sigma_c),abs(sigma_r)), sigma_c )
    elseif ( ((h(km1) + h(kp1)) + 2.0*h(k)) > 0.0 ) then
      sigma_l = ( u(k) - u(km1) )
      sigma_c = ( u(kp1) - u(km1) ) * ( h(k) / ((h(km1) + h(kp1)) + 2.0*h(k)) )
      sigma_r = ( u(kp1) - u(k) )

      ! The limiter is used in the local coordinate system to each cell, so for convenience store
      ! the slope times a half grid spacing.  (See White and Adcroft JCP 2008 Eqs 19 and 20)
      if ( (sigma_l * sigma_r) > 0.0 ) &
        slope_x_h = sign( min(abs(sigma_l),abs(sigma_c),abs(sigma_r)), sigma_c )
    endif

    ! Limit the edge values
    if ( (u(km1)-edge_val(k,1)) * (edge_val(k,1)-u(k)) < 0.0 ) then
      edge_val(k,1) = u(k) - sign( min( abs(slope_x_h), abs(edge_val(k,1)-u(k)) ), slope_x_h )
    endif

    if ( (u(kp1)-edge_val(k,2)) * (edge_val(k,2)-u(k)) < 0.0 ) then
      edge_val(k,2) = u(k) + sign( min( abs(slope_x_h), abs(edge_val(k,2)-u(k)) ), slope_x_h )
    endif

    ! Finally bound by neighboring cell means in case of roundoff
    edge_val(k,1) = max( min( edge_val(k,1), max(u(km1), u(k)) ), min(u(km1), u(k)) )
    edge_val(k,2) = max( min( edge_val(k,2), max(u(kp1), u(k)) ), min(u(kp1), u(k)) )

  enddo ! loop on interior edges

end subroutine bound_edge_values

!> Replace discontinuous collocated edge values with their average
!!
!! For each interior edge, check whether the edge values are discontinuous.
!! If so, compute the average and replace the edge values by the average.
subroutine average_discontinuous_edge_values( N, edge_val )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N,2), intent(inout) :: edge_val !< Edge values that may be modified [A]; the
                                           !! second index is for the two edges of each cell.
  ! Local variables
  integer       :: k            ! loop index
  real          :: u0_avg       ! avg value at given edge

  ! Loop on interior edges
  do k = 1,N-1
    ! Compare edge values on the right and left sides of the edge
    if ( edge_val(k,2) /= edge_val(k+1,1) ) then
      u0_avg = 0.5 * ( edge_val(k,2) + edge_val(k+1,1) )
      edge_val(k,2) = u0_avg
      edge_val(k+1,1) = u0_avg
    endif

  enddo ! end loop on interior edges

end subroutine average_discontinuous_edge_values

!> Check discontinuous edge values and replace them with their average if not monotonic
!!
!! For each interior edge, check whether the edge values are discontinuous.
!! If so and if they are not monotonic, replace each edge value by their average.
subroutine check_discontinuous_edge_values( N, u, edge_val )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: u !< cell averages in arbitrary units [A]
  real, dimension(N,2), intent(inout) :: edge_val !< Cell edge values [A]; the
                                           !! second index is for the two edges of each cell.
  ! Local variables
  integer       :: k            ! loop index
  real          :: u0_avg       ! avg value at given edge [A]

  do k = 1,N-1
    if ( (edge_val(k+1,1) - edge_val(k,2)) * (u(k+1) - u(k)) < 0.0 ) then
      u0_avg = 0.5 * ( edge_val(k,2) + edge_val(k+1,1) )
      u0_avg = max( min( u0_avg, max(u(k), u(k+1)) ), min(u(k), u(k+1)) )
      edge_val(k,2) = u0_avg
      edge_val(k+1,1) = u0_avg
    endif
  enddo ! end loop on interior edges

end subroutine check_discontinuous_edge_values


!> Compute h2 edge values (explicit second order accurate)
!! in the same units as h.
!
!! Compute edge values based on second-order explicit estimates.
!! These estimates are based on a straight line spanning two cells and evaluated
!! at the location of the middle edge. An interpolant spanning cells
!! k-1 and k is evaluated at edge k-1/2. The estimate for each edge is unique.
!!
!!       k-1     k
!! ..--o------o------o--..
!!          k-1/2
!!
!! Boundary edge values are set to be equal to the boundary cell averages.
subroutine edge_values_explicit_h2( N, h, u, edge_val )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< cell widths [H]
  real, dimension(N),   intent(in)    :: u !< cell average properties in arbitrary units [A]
  real, dimension(N,2), intent(inout) :: edge_val !< Returned edge values [A]; the
                                           !! second index is for the two edges of each cell.

  ! Local variables
  integer   :: k        ! loop index

  ! Boundary edge values are simply equal to the boundary cell averages
  edge_val(1,1) = u(1)
  edge_val(N,2) = u(N)

  do k = 2,N
    ! Compute left edge value
    if (h(k-1) + h(k) == 0.0) then    ! Avoid singularities when h0+h1=0
      edge_val(k,1) = 0.5 * (u(k-1) + u(k))
    else
      edge_val(k,1) = ( u(k-1)*h(k) + u(k)*h(k-1) ) / ( h(k-1) + h(k) )
    endif

    ! Left edge value of the current cell is equal to right edge value of left cell
    edge_val(k-1,2) = edge_val(k,1)
  enddo

end subroutine edge_values_explicit_h2

!> Compute h4 edge values (explicit fourth order accurate)
!! in the same units as h.
!!
!! Compute edge values based on fourth-order explicit estimates.
!! These estimates are based on a cubic interpolant spanning four cells
!! and evaluated at the location of the middle edge. An interpolant spanning
!! cells i-2, i-1, i and i+1 is evaluated at edge i-1/2. The estimate for
!! each edge is unique.
!!
!!       i-2    i-1     i     i+1
!! ..--o------o------o------o------o--..
!!                 i-1/2
!!
!! The first two edge values are estimated by evaluating the first available
!! cubic interpolant, i.e., the interpolant spanning cells 1, 2, 3 and 4.
!! Similarly, the last two edge values are estimated by evaluating the last
!! available interpolant.
!!
!! For this fourth-order scheme, at least four cells must exist.
subroutine edge_values_explicit_h4( N, h, u, edge_val, h_neglect, answers_2018 )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< cell widths [H]
  real, dimension(N),   intent(in)    :: u !< cell average properties in arbitrary units [A]
  real, dimension(N,2), intent(inout) :: edge_val !< Returned edge values [A]; the second index
                                           !! is for the two edges of each cell.
  real,       optional, intent(in)    :: h_neglect !< A negligibly small width [H]
  logical,    optional, intent(in)    :: answers_2018 !< If true use older, less acccurate expressions.

  ! Local variables
  integer               :: i, j
  real                  :: h0, h1, h2, h3   ! temporary thicknesses [H]
  real                  :: h_sum            ! A sum of adjacent thicknesses [H]
  real                  :: h_min            ! A minimal cell width [H]
  real                  :: f1, f2, f3       ! auxiliary variables with various units
  real                  :: et1, et2, et3    ! terms the expresson for edge values [A H]
  real, dimension(5)    :: x          ! Coordinate system with 0 at edges [H]
  real, parameter       :: C1_12 = 1.0 / 12.0
  real                  :: dx, xavg   ! Differences and averages of successive values of x [same units as h]
  real, dimension(4,4)  :: A                ! values near the boundaries
  real, dimension(4)    :: B, C
  real      :: hNeglect ! A negligible thickness in the same units as h.
  logical   :: use_2018_answers  ! If true use older, less acccurate expressions.

  use_2018_answers = .true. ; if (present(answers_2018)) use_2018_answers = answers_2018
  hNeglect = hNeglect_edge_dflt ; if (present(h_neglect)) hNeglect = h_neglect

  ! Loop on interior cells
  do i = 3,N-1

    h0 = h(i-2)
    h1 = h(i-1)
    h2 = h(i)
    h3 = h(i+1)

    ! Avoid singularities when consecutive pairs of h vanish
    if (h0+h1==0.0 .or. h1+h2==0.0 .or. h2+h3==0.0) then
      if (use_2018_answers) then
        h_min = hMinFrac*max( hNeglect, h0+h1+h2+h3 )
      else
        h_min = hMinFrac*max( hNeglect, (h0+h1)+(h2+h3) )
      endif
      h0 = max( h_min, h(i-2) )
      h1 = max( h_min, h(i-1) )
      h2 = max( h_min, h(i) )
      h3 = max( h_min, h(i+1) )
    endif

    if (use_2018_answers) then
      f1 = (h0+h1) * (h2+h3) / (h1+h2)
      f2 = h2 * u(i-1) + h1 * u(i)
      f3 = 1.0 / (h0+h1+h2) + 1.0 / (h1+h2+h3)
      et1 = f1 * f2 * f3
    else
      et1 = ( (h0+h1) * (h2+h3) * ((h1+h2+h3) + (h0+h1+h2)) / &
              (((h1+h2) * ((h0+h1+h2) * (h1+h2+h3)))) ) * &
            (h2 * u(i-1) + h1 * u(i))
    endif

    et2 = ( h2 * (h2+h3) / ( (h0+h1+h2)*(h0+h1) ) ) * &
          ((h0+2.0*h1) * u(i-1) - h1 * u(i-2))

    et3 = ( h1 * (h0+h1) / ( (h1+h2+h3)*(h2+h3) ) ) * &
          ((2.0*h2+h3) * u(i) - h2 * u(i+1))

    if (use_2018_answers) then
      edge_val(i,1) = (et1 + et2 + et3) / ( h0 + h1 + h2 + h3)
    else
      edge_val(i,1) = (et1 + (et2 + et3)) / ((h0 + h1) + (h2 + h3))
    endif
    edge_val(i-1,2) = edge_val(i,1)

  enddo ! end loop on interior cells

  ! Determine first two edge values
  if (use_2018_answers) then
    h_min = max( hNeglect, hMinFrac*sum(h(1:4)) )
    x(1) = 0.0
    do i = 1,4
      dx = max(h_min, h(i) )
      x(i+1) = x(i) + dx
      do j = 1,4 ; A(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / real(j) ; enddo
      B(i) = u(i) * dx
    enddo

    call solve_linear_system( A, B, C, 4 )

    ! Set the edge values of the first cell
    edge_val(1,1) = evaluation_polynomial( C, 4, x(1) )
    edge_val(1,2) = evaluation_polynomial( C, 4, x(2) )
  else  ! Use expressions with less sensitivity to roundoff
    h_min = hMinFrac*((h(1) + h(2)) + (h(3) + h(4)))
    if (h_min == 0.0) h_min = 1.0  ! Handle the case of all massless layers.
    x(1) = 0.0
    do i = 1,4
      dx = max(h_min, h(i) )
      x(i+1) = x(i) + dx
      xavg = 0.5 * (x(i+1) + x(i))
      A(i,1) = dx
      A(i,2) = dx * xavg
      A(i,3) = dx * (xavg**2 + C1_12*dx**2)
      A(i,4) = dx * xavg * (xavg**2 + 0.25*dx**2)
      B(i) = u(i) * dx
    enddo

    call solve_linear_system( A, B, C, 4 )

    ! Set the edge values of the first cell
    edge_val(1,1) = C(1) ! x(1) = 0 so ignore + x(1)*(C(2) + x(1)*(C(3) + x(1)*C(4)))
    edge_val(1,2) = C(1) + x(2)*(C(2) + x(2)*(C(3) + x(2)*C(4)))
  endif
  edge_val(2,1) = edge_val(1,2)

  ! Determine two edge values of the last cell
  if (use_2018_answers) then
    h_min = max( hNeglect, hMinFrac*sum(h(N-3:N)) )

    x(1) = 0.0
    do i = 1,4
      dx = max(h_min, h(N-4+i) )
      x(i+1) = x(i) + dx
      do j = 1,4 ; A(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / real(j) ; enddo
      B(i) = u(N-4+i) * dx
    enddo

    call solve_linear_system( A, B, C, 4 )

    ! Set the last and second to last edge values
    edge_val(N,2) = evaluation_polynomial( C, 4, x(5) )
    edge_val(N,1) = evaluation_polynomial( C, 4, x(4) )
  else
    ! Use expressions with less sensitivity to roundoff, including using a coordinate
    ! system that sets the origin at the last interface in the domain.
    h_min = hMinFrac * ((h(N-3) + h(N-2)) + (h(N-1) + h(N)))
    if (h_min == 0.0) h_min = 1.0  ! Handle the case of all massless layers.

    x(1) = 0.0

    do i=1,4
      dx = max(h_min, h(N+1-i) )
      x(i+1) = x(i) + dx
      xavg = x(i) + 0.5*dx

      A(i,1) = dx
      A(i,2) = dx * xavg
      A(i,3) = dx * (xavg**2 + C1_12*dx**2)
      A(i,4) = dx * xavg * (xavg**2 + 0.25*dx**2)

      B(i) = u(N+1-i) * dx
    enddo

    call solve_linear_system( A, B, C, 4 )

    ! Set the last and second to last edge values
    edge_val(N,2) = C(1)
    edge_val(N,1) = C(1) + x(2)*(C(2) + x(2)*(C(3) + x(2)*C(4)))
  endif
  edge_val(N-1,2) = edge_val(N,1)

end subroutine edge_values_explicit_h4

!> Compute ih4 edge values (implicit fourth order accurate)
!! in the same units as h.
!!
!! Compute edge values based on fourth-order implicit estimates.
!!
!! Fourth-order implicit estimates of edge values are based on a two-cell
!! stencil. A tridiagonal system is set up and is based on expressing the
!! edge values in terms of neighboring cell averages. The generic
!! relationship is
!!
!! \f[
!! \alpha u_{i-1/2} + u_{i+1/2} + \beta u_{i+3/2} = a \bar{u}_i + b \bar{u}_{i+1}
!! \f]
!!
!! and the stencil looks like this
!!
!!          i     i+1
!!   ..--o------o------o--..
!!     i-1/2  i+1/2  i+3/2
!!
!! In this routine, the coefficients \f$\alpha\f$, \f$\beta\f$, \f$a\f$ and \f$b\f$ are
!! computed, the tridiagonal system is built, boundary conditions are prescribed and
!! the system is solved to yield edge-value estimates.
!!
!! There are N+1 unknowns and we are able to write N-1 equations. The
!! boundary conditions close the system.
subroutine edge_values_implicit_h4( N, h, u, edge_val, h_neglect, answers_2018 )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< cell widths [H]
  real, dimension(N),   intent(in)    :: u !< cell average properties in arbitrary units [A]
  real, dimension(N,2), intent(inout) :: edge_val !< Returned edge values [A]; the second index
                                           !! is for the two edges of each cell.
  real,       optional, intent(in)    :: h_neglect !< A negligibly small width [H]
  logical,    optional, intent(in)    :: answers_2018 !< If true use older, less acccurate expressions.

  ! Local variables
  integer               :: i, j                 ! loop indexes
  real                  :: h0, h1               ! cell widths [H]
  real                  :: h_min                ! A minimal cell width [H]
  real                  :: h_sum                ! A sum of adjacent thicknesses [H]
  real                  :: h0_2, h1_2, h0h1
  real                  :: d2, d4
  real                  :: alpha, beta          ! stencil coefficients [nondim]
  real                  :: I_h2, abmix          ! stencil coefficients [nondim]
  real                  :: a, b
  real, dimension(5)    :: x                    ! Coordinate system with 0 at edges [H]
  real, parameter       :: C1_12 = 1.0 / 12.0
  real                  :: dx, xavg             ! Differences and averages of successive values of x [H]
  real, dimension(4,4)  :: Asys                 ! boundary conditions
  real, dimension(4)    :: Bsys, Csys
  real, dimension(N+1)  :: tri_l, &     ! tridiagonal system (lower diagonal) [nondim]
                           tri_d, &     ! tridiagonal system (middle diagonal) [nondim]
                           tri_c, &     ! tridiagonal system central value, with tri_d = tri_c+tri_l+tri_u
                           tri_u, &     ! tridiagonal system (upper diagonal) [nondim]
                           tri_b, &     ! tridiagonal system (right hand side) [A]
                           tri_x        ! tridiagonal system (solution vector) [A]
  real      :: hNeglect          ! A negligible thickness [H]
  logical   :: use_2018_answers  ! If true use older, less acccurate expressions.

  use_2018_answers = .true. ; if (present(answers_2018)) use_2018_answers = answers_2018
  hNeglect = hNeglect_edge_dflt ; if (present(h_neglect)) hNeglect = h_neglect

  ! Loop on cells (except last one)
  do i = 1,N-1

    ! Get cell widths
    h0 = h(i)
    h1 = h(i+1)

    if (use_2018_answers) then
      ! Avoid singularities when h0+h1=0
      if (h0+h1==0.) then
        h0 = hNeglect
        h1 = hNeglect
      endif

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

      tri_d(i+1) = 1.0
    else  ! Use expressions with less sensitivity to roundoff
      if (h0+h1==0.) then  ! Avoid singularities when h0+h1=0
        alpha = 0.25 ; beta = 0.25 ; abmix = 0.25
      else
        ! The 1e-12 here attempts to balance truncation errors from the differences of
        ! large numbers against errors from approximating thin layers as non-vanishing.
        if (abs(h0) < 1.0e-12*abs(h1)) h0 = 1.0e-12*h1
        if (abs(h1) < 1.0e-12*abs(h0)) h1 = 1.0e-12*h0
        I_h2 = 1.0 / ((h0 + h1)**2)
        alpha = (h1 * h1) * I_h2
        beta = (h0 * h0) * I_h2
        abmix = (h0 * h1) * I_h2
      endif
      a = 2.0 * alpha * ( alpha + 2.0 * beta + 3.0 * abmix )
      b = 2.0 * beta * ( beta + 2.0 * alpha + 3.0 * abmix )

      tri_c(i+1) = 2.0*abmix  ! = 1.0 - alpha - beta
    endif

    tri_l(i+1) = alpha
    tri_u(i+1) = beta

    tri_b(i+1) = a * u(i) + b * u(i+1)

  enddo ! end loop on cells

  ! Boundary conditions: set the first boundary value
  if (use_2018_answers) then
    h_min = max( hNeglect, hMinFrac*sum(h(1:4)) )
    x(1) = 0.0
    do i = 1,4
      dx = max(h_min, h(i) )
      x(i+1) = x(i) + dx
      do j = 1,4 ; Asys(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / j ; enddo
      Bsys(i) = u(i) * dx
    enddo

    call solve_linear_system( Asys, Bsys, Csys, 4 )

    tri_b(1) = evaluation_polynomial( Csys, 4, x(1) )  ! Set the first edge value
    tri_d(1) = 1.0
  else ! Use expressions with less sensitivity to roundoff
    h_min = max( hNeglect, hMinFrac * ((h(1) + h(2)) + (h(3) + h(4))) )
    x(1) = 0.0
    do i = 1,4
      dx = max(h_min, h(i) )
      x(i+1) = x(i) + dx
      xavg = x(i) + 0.5*dx
      Asys(i,1) = dx
      Asys(i,2) = dx * xavg
      Asys(i,3) = dx * (xavg**2 + C1_12*dx**2)
      Asys(i,4) = dx * xavg * (xavg**2 + 0.25*dx**2)
      Bsys(i) = u(i) * dx
    enddo

    call solve_linear_system( Asys, Bsys, Csys, 4 )

    tri_b(1) = Csys(1)  ! Set the first edge value, using the fact that x(1) = 0.
    tri_c(1) = 1.0
  endif
  tri_u(1) = 0.0 ! tri_l(1) = 0.0

  ! Boundary conditions: set the last boundary value
  if (use_2018_answers) then
    h_min = max( hNeglect, hMinFrac*sum(h(N-3:N)) )
    x(1) = 0.0
    do i=1,4
      dx = max(h_min, h(N-4+i) )
      x(i+1) = x(i) + dx
      do j = 1,4 ; Asys(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / j ; enddo
      Bsys(i) = u(N-4+i) * dx
    enddo

    call solve_linear_system( Asys, Bsys, Csys, 4 )

    ! Set the last edge value
    tri_b(N+1) = evaluation_polynomial( Csys, 4, x(5) )
    tri_d(N+1) = 1.0

  else
    ! Use expressions with less sensitivity to roundoff, including using a coordinate
    ! system that sets the origin at the last interface in the domain.
    h_min = max( hNeglect, hMinFrac * ((h(N-3) + h(N-2)) + (h(N-1) + h(N))) )
    x(1) = 0.0
    do i=1,4
      dx = max(h_min, h(N+1-i) )
      x(i+1) = x(i) + dx
      xavg = x(i) + 0.5*dx

      Asys(i,1) = dx
      Asys(i,2) = dx * xavg
      Asys(i,3) = dx * (xavg**2 + C1_12*dx**2)
      Asys(i,4) = dx * xavg * (xavg**2 + 0.25*dx**2)

      Bsys(i) = u(N+1-i) * dx
    enddo

    call solve_linear_system( Asys, Bsys, Csys, 4 )

    ! Set the last edge value
    tri_b(N+1) = Csys(1)
    tri_c(N+1) = 1.0
  endif
  tri_l(N+1) = 0.0 ! tri_u(N+1) = 0.0

  ! Solve tridiagonal system and assign edge values
  if (use_2018_answers) then
    call solve_tridiagonal_system( tri_l, tri_d, tri_u, tri_b, tri_x, N+1 )
  else
    call solve_diag_dominant_tridiag( tri_l, tri_c, tri_u, tri_b, tri_x, N+1 )
  endif

  edge_val(1,1) = tri_x(1)
  do i=2,N
    edge_val(i,1)   = tri_x(i)
    edge_val(i-1,2) = tri_x(i)
  enddo
  edge_val(N,2) = tri_x(N+1)

end subroutine edge_values_implicit_h4

!> Compute ih6 edge values (implicit sixth order accurate)
                                           !! in the same units as h.
!!
!! Sixth-order implicit estimates of edge values are based on a four-cell,
!! three-edge stencil. A tridiagonal system is set up and is based on
!! expressing the edge values in terms of neighboring cell averages.
!!
!! The generic relationship is
!!
!! \f[
!! \alpha u_{i-1/2} + u_{i+1/2} + \beta u_{i+3/2} =
!! a \bar{u}_{i-1} + b \bar{u}_i + c \bar{u}_{i+1} + d \bar{u}_{i+2}
!! \f]
!!
!! and the stencil looks like this
!!
!!         i-1     i     i+1    i+2
!!   ..--o------o------o------o------o--..
!!            i-1/2  i+1/2  i+3/2
!!
!! In this routine, the coefficients \f$\alpha\f$, \f$\beta\f$, a, b, c and d are
!! computed, the tridiagonal system is built, boundary conditions are
!! prescribed and the system is solved to yield edge-value estimates.
!!
!! Note that the centered stencil only applies to edges 3 to N-1 (edges are
!! numbered 1 to n+1), which yields N-3 equations for N+1 unknowns. Two other
!! equations are written by using a right-biased stencil for edge 2 and a
!! left-biased stencil for edge N. The prescription of boundary conditions
!! (using sixth-order polynomials) closes the system.
!!
!! CAUTION: For each edge, in order to determine the coefficients of the
!!          implicit expression, a 6x6 linear system is solved. This may
!!          become computationally expensive if regridding is carried out
!!          often. Figuring out closed-form expressions for these coefficients
!!          on nonuniform meshes turned out to be intractable.
subroutine edge_values_implicit_h6( N, h, u, edge_val, h_neglect, answers_2018 )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< cell widths [H]
  real, dimension(N),   intent(in)    :: u !< cell average properties (size N) in arbitrary units [A]
  real, dimension(N,2), intent(inout) :: edge_val  !< Returned edge values [A]; the second index
                                           !! is for the two edges of each cell.
  real,       optional, intent(in)    :: h_neglect !< A negligibly small width [H]
  logical,    optional, intent(in)    :: answers_2018 !< If true use older, less acccurate expressions.

  ! Local variables
  integer               :: i, j, k              ! loop indexes
  real                  :: h0, h1, h2, h3       ! cell widths [H]
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
  real, dimension(7)    :: x          ! Coordinate system with 0 at edges [same units as h]
  real, parameter       :: C1_12 = 1.0 / 12.0
  real, parameter       :: C5_6 = 5.0 / 6.0
  real                  :: dx, xavg   ! Differences and averages of successive values of x [same units as h]
  real, dimension(6,6)  :: Asys                 ! boundary conditions
  real, dimension(6)    :: Bsys, Csys           ! ...
  real, dimension(N+1)  :: tri_l, &             ! trid. system (lower diagonal)
                           tri_d, &             ! trid. system (middle diagonal)
                           tri_u, &             ! trid. system (upper diagonal)
                           tri_b, &             ! trid. system (unknowns vector)
                           tri_x                ! trid. system (rhs)
  real      :: hNeglect          ! A negligible thickness [H].
  logical   :: use_2018_answers  ! If true use older, less acccurate expressions.

  use_2018_answers = .true. ; if (present(answers_2018)) use_2018_answers = answers_2018
  hNeglect = hNeglect_edge_dflt ; if (present(h_neglect)) hNeglect = h_neglect

  ! Loop on cells (except last one)
  do k = 2,N-2

    ! Cell widths
    h0 = h(k-1)
    h1 = h(k+0)
    h2 = h(k+1)
    h3 = h(k+2)

    ! Avoid singularities when h0=0 or h3=0
    if (h0*h3==0.) then
      g = max( hNeglect, h0+h1+h2+h3 )
      h0 = max( hMinFrac*g, h0 )
      h1 = max( hMinFrac*g, h1 )
      h2 = max( hMinFrac*g, h2 )
      h3 = max( hMinFrac*g, h3 )
    endif

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

    tri_l(k+1) = alpha
    tri_d(k+1) = 1.0
    tri_u(k+1) = beta
    tri_b(k+1) = a * u(k-1) + b * u(k) + c * u(k+1) + d * u(k+2)

  enddo ! end loop on cells

  ! Use a right-biased stencil for the second row

  ! Cell widths
  h0 = h(1)
  h1 = h(2)
  h2 = h(3)
  h3 = h(4)

  ! Avoid singularities when h0=0 or h3=0
  if (h0*h3==0.) then
    g = max( hNeglect, h0+h1+h2+h3 )
    h0 = max( hMinFrac*g, h0 )
    h1 = max( hMinFrac*g, h1 )
    h2 = max( hMinFrac*g, h2 )
    h3 = max( hMinFrac*g, h3 )
  endif

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

  tri_l(2) = alpha
  tri_d(2) = 1.0
  tri_u(2) = beta
  tri_b(2) = a * u(1) + b * u(2) + c * u(3) + d * u(4)

  ! Boundary conditions: left boundary
!  h_sum = (h(1) + h(2)) + (h(5) + h(6)) + (h(3) + h(4))
  g = max( hNeglect, hMinFrac*sum(h(1:6)) )
  x(1) = 0.0
  do i = 2,7
    x(i) = x(i-1) + max( g, h(i-1) )
  enddo

  do i = 1,6
    dx = max( g, h(i) )
    if (use_2018_answers) then
      do j = 1,6 ; Asys(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / j ; enddo
    else  ! Use expressions with less sensitivity to roundoff
      xavg = 0.5 * (x(i+1) + x(i))
      Asys(i,1) = dx
      Asys(i,2) = dx * xavg
      Asys(i,3) = dx * (xavg**2 + C1_12*dx**2)
      Asys(i,4) = dx * xavg * (xavg**2 + 0.25*dx**2)
      Asys(i,5) = dx * (xavg**4 + 0.5*xavg**2*dx**2 + 0.0125*dx**4)
      Asys(i,6) = dx * xavg * (xavg**4 + C5_6*xavg**2*dx**2 + 0.0625*dx**4)
    endif
    Bsys(i) = u(i) * dx

  enddo

  call solve_linear_system( Asys, Bsys, Csys, 6 )

  tri_l(1) = 0.0
  tri_d(1) = 1.0
  tri_u(1) = 0.0
  tri_b(1) = evaluation_polynomial( Csys, 6, x(1) )        ! first edge value

  ! Use a left-biased stencil for the second to last row

  ! Cell widths
  h0 = h(N-3)
  h1 = h(N-2)
  h2 = h(N-1)
  h3 = h(N)

  ! Avoid singularities when h0=0 or h3=0
  if (h0*h3==0.) then
    g = max( hNeglect, h0+h1+h2+h3 )
    h0 = max( hMinFrac*g, h0 )
    h1 = max( hMinFrac*g, h1 )
    h2 = max( hMinFrac*g, h2 )
    h3 = max( hMinFrac*g, h3 )
  endif

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

  tri_l(N) = alpha
  tri_d(N) = 1.0
  tri_u(N) = beta
  tri_b(N) = a * u(N-3) + b * u(N-2) + c * u(N-1) + d * u(N)

  ! Boundary conditions: right boundary
!  h_sum = (h(N-3) + h(N-2)) + ((h(N-1) + h(N)) + (h(N-5) + h(N-4)))
  g = max( hNeglect, hMinFrac*sum(h(N-5:N)) )
  x(1) = 0.0
  do i = 2,7
    x(i) = x(i-1) + max( g, h(N-7+i) )
  enddo

  do i = 1,6
    dx = max( g, h(N-6+i) )
    if (use_2018_answers) then
      do j = 1,6 ; Asys(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / j ; enddo
    else  ! Use expressions with less sensitivity to roundoff
      xavg = 0.5 * (x(i+1) + x(i))
      Asys(i,1) = dx
      Asys(i,2) = dx * xavg
      Asys(i,3) = dx * (xavg**2 + C1_12*dx**2)
      Asys(i,4) = dx * xavg * (xavg**2 + 0.25*dx**2)
      Asys(i,5) = dx * (xavg**4 + 0.5*xavg**2*dx**2 + 0.0125*dx**4)
      Asys(i,6) = dx * xavg * (xavg**4 + C5_6*xavg**2*dx**2 + 0.0625*dx**4)
    endif
    Bsys(i) = u(N-6+i) * dx

  enddo

  call solve_linear_system( Asys, Bsys, Csys, 6 )

  tri_l(N+1) = 0.0
  tri_d(N+1) = 1.0
  tri_u(N+1) = 0.0
  tri_b(N+1) = evaluation_polynomial( Csys, 6, x(7) )      ! last edge value

  ! Solve tridiagonal system and assign edge values
  call solve_tridiagonal_system( tri_l, tri_d, tri_u, tri_b, tri_x, N+1 )

  do i = 2,N
    edge_val(i,1)   = tri_x(i)
    edge_val(i-1,2) = tri_x(i)
  enddo
  edge_val(1,1) = tri_x(1)
  edge_val(N,2) = tri_x(N+1)

end subroutine edge_values_implicit_h6

end module regrid_edge_values
