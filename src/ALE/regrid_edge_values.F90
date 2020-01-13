!> Edge value estimation for high-order resconstruction
module regrid_edge_values

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL
use regrid_solvers, only : solve_linear_system, solve_tridiagonal_system
use regrid_solvers, only : solve_diag_dominant_tridiag, linear_solver
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
  real, dimension(4)    :: dz               ! A temporary array of limited layer thicknesses [H]
  real, dimension(4)    :: u_tmp            ! A temporary array of cell average properties [A]
  real, parameter       :: C1_12 = 1.0 / 12.0
  real                  :: dx, xavg   ! Differences and averages of successive values of x [same units as h]
  real, dimension(4,4)  :: A                ! values near the boundaries
  real, dimension(4)    :: B, C
  real      :: hNeglect ! A negligible thickness in the same units as h.
  logical   :: use_2018_answers  ! If true use older, less acccurate expressions.

  use_2018_answers = .true. ; if (present(answers_2018)) use_2018_answers = answers_2018
  if (use_2018_answers) then
    hNeglect = hNeglect_edge_dflt ; if (present(h_neglect)) hNeglect = h_neglect
  else
    hNeglect = hNeglect_dflt ; if (present(h_neglect)) hNeglect = h_neglect
  endif

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

    do i=1,4
      dz(i) = max(h_min, h(i) )
      u_tmp(i) = u(i)
    enddo
    call end_value_h4(dz, u_tmp, C)

    ! Set the edge values of the first cell
    edge_val(1,1) = C(1)
    edge_val(1,2) = C(1) + dz(1)*(C(2) + dz(1)*(C(3) + dz(1)*C(4)))
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
    do i=1,4
      dz(i) = max(h_min, h(N+1-i) )
      u_tmp(i) = u(N+1-i)
    enddo
    call end_value_h4(dz, u_tmp, C)

    ! Set the last and second to last edge values
    edge_val(N,2) = C(1)
    edge_val(N,1) = C(1) + dz(1)*(C(2) + dz(1)*(C(3) + dz(1)*C(4)))
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
  real                  :: h0, h1, h2           ! cell widths [H]
  real                  :: h_min                ! A minimal cell width [H]
  real                  :: h_sum                ! A sum of adjacent thicknesses [H]
  real                  :: h0_2, h1_2, h0h1
  real                  :: h0ph1_2, h0ph1_4
  real                  :: alpha, beta          ! stencil coefficients [nondim]
  real                  :: I_h2, abmix          ! stencil coefficients [nondim]
  real                  :: a, b
  real, dimension(5)    :: x                    ! Coordinate system with 0 at edges [H]
  real, parameter       :: C1_12 = 1.0 / 12.0
  real, parameter       :: C1_3 = 1.0 / 3.0
  real, dimension(4)    :: dz                   ! A temporary array of limited layer thicknesses [H]
  real, dimension(4)    :: u_tmp                ! A temporary array of cell average properties [A]
  real                  :: dx, xavg             ! Differences and averages of successive values of x [H]
  real, dimension(4,4)  :: Asys                 ! boundary conditions
  real, dimension(4)    :: Bsys, Csys
  real, dimension(4,4)  :: Asys_orig            ! boundary conditions
  real, dimension(4)    :: Bsys_orig
  real, dimension(N+1)  :: tri_l, &     ! tridiagonal system (lower diagonal) [nondim]
                           tri_d, &     ! tridiagonal system (middle diagonal) [nondim]
                           tri_c, &     ! tridiagonal system central value, with tri_d = tri_c+tri_l+tri_u
                           tri_u, &     ! tridiagonal system (upper diagonal) [nondim]
                           tri_b, &     ! tridiagonal system (right hand side) [A]
                           tri_x        ! tridiagonal system (solution vector) [A]
  real      :: hNeglect          ! A negligible thickness [H]
  logical   :: use_2018_answers  ! If true use older, less acccurate expressions.

  use_2018_answers = .true. ; if (present(answers_2018)) use_2018_answers = answers_2018
  if (use_2018_answers) then
    hNeglect = hNeglect_edge_dflt ; if (present(h_neglect)) hNeglect = h_neglect
  else
    hNeglect = hNeglect_dflt ; if (present(h_neglect)) hNeglect = h_neglect
  endif

  ! Loop on cells (except last one)
  do i = 1,N-1
    if (use_2018_answers) then
      ! Get cell widths
      h0 = h(i)
      h1 = h(i+1)
      ! Avoid singularities when h0+h1=0
      if (h0+h1==0.) then
        h0 = hNeglect
        h1 = hNeglect
      endif

      ! Auxiliary calculations
      h0ph1_2 = (h0 + h1)**2
      h0ph1_4 = h0ph1_2**2
      h0h1 = h0 * h1
      h0_2 = h0 * h0
      h1_2 = h1 * h1

      ! Coefficients
      alpha = h1_2 / h0ph1_2
      beta = h0_2 / h0ph1_2
      a = 2.0 * h1_2 * ( h1_2 + 2.0 * h0_2 + 3.0 * h0h1 ) / h0ph1_4
      b = 2.0 * h0_2 * ( h0_2 + 2.0 * h1_2 + 3.0 * h0h1 ) / h0ph1_4

      tri_d(i+1) = 1.0
    else  ! Use expressions with less sensitivity to roundoff
      ! Get cell widths
      h0 = max(h(i), hNeglect)
      h1 = max(h(i+1), hNeglect)
      ! The 1e-12 here attempts to balance truncation errors from the differences of
      ! large numbers against errors from approximating thin layers as non-vanishing.
      if (abs(h0) < 1.0e-12*abs(h1)) h0 = 1.0e-12*h1
      if (abs(h1) < 1.0e-12*abs(h0)) h1 = 1.0e-12*h0
      I_h2 = 1.0 / ((h0 + h1)**2)
      alpha = (h1 * h1) * I_h2
      beta = (h0 * h0) * I_h2
      abmix = (h0 * h1) * I_h2
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
    do i=1,4
      dz(i) = max(h_min, h(i) )
      u_tmp(i) = u(i)
    enddo
    call end_value_h4(dz, u_tmp, Csys)

    tri_b(1) = Csys(1)  ! Set the first edge value.

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
    do i=1,4
      dz(i) = max(h_min, h(N+1-i) )
      u_tmp(i) = u(N+1-i)
    enddo

    call end_value_h4(dz, u_tmp, Csys)

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

!> Determine a one-sided 4th order polynomial fit of u to the data points for the purposes of specifying
!! edge values, as described in the appendix of White and Adcroft JCP 2008.
subroutine end_value_h4(dz, u, Csys)
  real, dimension(4), intent(in)  :: dz    !< The thicknesses of 4 layers, starting at the edge [H].
                                           !! The values of dz must be positive.
  real, dimension(4), intent(in)  :: u     !< The average properties of 4 layers, starting at the edge [A]
  real, dimension(4), intent(out) :: Csys  !< The four coefficients of a 4th order polynomial fit
                                           !! of u as a function of z [A H-(n-1)]

  ! Local variables
  real :: Wt(3,4)         ! The weights of successive u differences in the 4 closed form expressions.
                          ! The units of Wt vary with the second index as [H-(n-1)].
  real :: h1, h2, h3, h4  ! Copies of the layer thicknesses [H]
  real :: h12, h23, h34   ! Sums of two successive thicknesses [H]
  real :: h123, h234      ! Sums of three successive thicknesses [H]
  real :: h1234           ! Sums of all four thicknesses [H]
  ! real :: I_h1          ! The inverse of the a thickness [H-1]
  real :: I_h12, I_h23, I_h34 ! The inverses of sums of two thicknesses [H-1]
  real :: I_h123, I_h234  ! The inverse of the sum of three thicknesses [H-1]
  real :: I_h1234         ! The inverse of the sum of all four thicknesses [H-1]
  real :: I_denom         ! The inverse of the denominator some expressions [H-3]
  real :: I_denB3         ! The inverse of the product of three sums of thicknesses [H-3]
  real, parameter :: C1_3 = 1.0 / 3.0
  integer :: i, j, k

  ! These are only used for code verification
  real, dimension(4) :: Atest  ! The  coefficients of an expression that is being tested.
  real :: zavg, u_mag, c_mag
  character(len=128) :: mesg
  real, parameter :: C1_12 = 1.0 / 12.0

 ! if ((dz(1) == dz(2)) .and. (dz(1) == dz(3)) .and. (dz(1) == dz(4))) then
 !   ! There are simple closed-form expressions in this case
 !   I_h1 = 0.0 ; if (dz(1) > 0.0) I_h1 = 1.0 / dz(1)
 !   Csys(1) = u(1) + (-13.0 * (u(2)-u(1)) + 10.0 * (u(3)-u(2)) - 3.0 * (u(4)-u(3))) * (0.25*C1_3)
 !   Csys(2) = (35.0 * (u(2)-u(1)) - 34.0 * (u(3)-u(2)) + 11.0 * (u(4)-u(3))) * (0.25*C1_3 * I_h1)
 !   Csys(3) = (-5.0 * (u(2)-u(1)) + 8.0 * (u(3)-u(2)) - 3.0 * (u(4)-u(3))) * (0.25 * I_h1**2)
 !   Csys(4) = ((u(2)-u(1)) - 2.0 * (u(3)-u(2)) + (u(4)-u(3))) * (0.5*C1_3)
 ! else

  ! Express the coefficients as sums of the differences between properties of succesive layers.

  h1 = dz(1) ; h2 = dz(2) ; h3 = dz(3) ; h4 = dz(4)
  h12 = h1+h2 ; h23 = h2+h3 ; h34 = h3+h4
  h123 = h12 + h3 ; h234 = h2 + h34 ; h1234 = h12 + h34
  ! Find 3 reciprocals with a single division for efficiency.
  I_denB3 = 1.0 / (h123 * h12 * h23)
  I_h12 = (h123 * h23) * I_denB3
  I_h23 = (h12 * h123) * I_denB3
  I_h123 = (h12 * h23) * I_denB3
  I_denom = 1.0 / ( h1234 * (h234 * h34) )
  I_h34 = (h1234 * h234) * I_denom
  I_h234 = (h1234 * h34) * I_denom
  I_h1234 = (h234 * h34) * I_denom

  ! Calculation coefficients in the four equations

  ! The expressions for Csys(3) and Csys(4) come from reducing the 4x4 matrix problem into the following 2x2
  ! matrix problem, then manipulating the analytic solution to avoid any subtraction and simplifying.
  !  (C1_3 * h123 * h23) * Csys(3) + (0.25 * h123 * h23 * (h3 + 2.0*h2 + 3.0*h1)) * Csys(4) =
  !            (u(3)-u(1)) - (u(2)-u(1)) * (h12 + h23) * I_h12
  !  (C1_3 * ((h23 + h34) * h1234 + h23 * h3)) * Csys(3) +
  !  (0.25 * ((h1234 + h123 + h12 + h1) * h23 * h3 + (h1234 + h12 + h1) * (h23 + h34) * h1234)) * Csys(4) =
  !            (u(4)-u(1)) - (u(2)-u(1)) * (h123 + h234) * I_h12
  ! The final expressions for Csys(1) and Csys(2) were derived by algebraically manipulating the following expressions:
  !  Csys(1) = (C1_3 * h1 * h12 * Csys(3) + 0.25 * h1 * h12 * (2.0*h1+h2) * Csys(4)) - &
  !            (h1*I_h12)*(u(2)-u(1)) + u(1)
  !  Csys(2) = (-2.0*C1_3 * (2.0*h1+h2) * Csys(3) - 0.5 * (h1**2 + h12 * (2.0*h1+h2)) * Csys(4)) + &
  !            2.0*I_h12 * (u(2)-u(1))
  ! These expressions are typically evaluated at x=0 and x=h1, so it is important that these are well behaved
  ! for these values, suggesting that h1/h23 and h1/h34 should not be allowed to be too large.

  Wt(1,1) = -h1 * (I_h1234 + I_h123 + I_h12)                              ! > -3
  Wt(2,1) =  h1 * h12 * ( I_h234 * I_h1234 + I_h23 * (I_h234 + I_h123) )  ! < (h1/h234) + (h1/h23)*(2+(h1/h234))
  Wt(3,1) = -h1 * h12 * h123 * I_denom                                    ! > -(h1/h34)*(1+(h1/h234))

  Wt(1,2) =  2.0 * (I_h12*(1.0 + (h1+h12) * (I_h1234 + I_h123)) + h1 * I_h1234*I_h123) ! < 10/h12
  Wt(2,2) = -2.0 * ((h1 * h12 * I_h1234) *       (I_h23 * (I_h234 + I_h123)) + &       ! > -(10+6*(h1/h234))/h23
                    (h1+h12) * ( I_h1234*I_h234 + I_h23 * (I_h234 + I_h123) ) )
  Wt(3,2) =  2.0 * ((h1+h12) * h123 + h1*h12 ) * I_denom                               ! < (2+(6*h1/h234)) / h34

  Wt(1,3) = -3.0 * I_h12 * I_h123* ( 1.0 + I_h1234 * ((h1+h12)+h123) )                 ! > -12 / (h12*h123)
  Wt(2,3) =  3.0 * I_h23 * ( I_h123 + I_h1234 * ((h1+h12)+h123) * (I_h123 + I_h234) )  ! < 12 / (h23^2)
  Wt(3,3) = -3.0 * ((h1+h12)+h123) * I_denom                                           ! > -9 / (h234*h23)

  Wt(1,4) =  4.0 * I_h1234 * I_h123 * I_h12                          ! Wt*h1^3 < 4
  Wt(2,4) = -4.0 * I_h1234 * (I_h23 * (I_h123 + I_h234))             ! Wt*h1^3 > -4* (h1/h23)*(1+h1/h234)
  Wt(3,4) =  4.0 * I_denom  ! = 4.0*I_h1234 * I_h234 * I_h34         ! Wt*h1^3 < 4 * (h1/h234)*(h1/h34)

  Csys(1) = ((u(1) + Wt(1,1) * (u(2)-u(1))) + Wt(2,1) * (u(3)-u(2))) + Wt(3,1) * (u(4)-u(3))
  Csys(2) = (Wt(1,2) * (u(2)-u(1)) + Wt(2,2) * (u(3)-u(2))) + Wt(3,2) * (u(4)-u(3))
  Csys(3) = (Wt(1,3) * (u(2)-u(1)) + Wt(2,3) * (u(3)-u(2))) + Wt(3,3) * (u(4)-u(3))
  Csys(4) = (Wt(1,4) * (u(2)-u(1)) + Wt(2,4) * (u(3)-u(2))) + Wt(3,4) * (u(4)-u(3))

  ! endif ! End of non-uniform layer thickness branch.

  ! To verify that these answers are correct, uncomment the following:
!  u_mag = 0.0 ; do i=1,4 ; u_mag = max(u_mag, abs(u(i))) ; enddo
!  do i = 1,4
!    if (i==1) then ; zavg = 0.5*dz(i) ; else ; zavg = zavg + 0.5*(dz(i-1)+dz(i)) ; endif
!    Atest(1) = 1.0
!    Atest(2) = zavg                             ! = ( (z(i+1)**2) - (z(i)**2) ) / (2*dz(i))
!    Atest(3) = (zavg**2 + 0.25*C1_3*dz(i)**2)   ! = ( (z(i+1)**3) - (z(i)**3) ) / (3*dz(i))
!    Atest(4) = zavg * (zavg**2 + 0.25*dz(i)**2) ! = ( (z(i+1)**4) - (z(i)**4) ) / (4*dz(i))
!    c_mag = 1.0 ; do k=0,3 ; do j=1,3 ; c_mag = c_mag + abs(Wt(j,k+1) * zavg**k) ; enddo ; enddo
!    write(mesg, '("end_value_h4 line ", i2, " c_mag = ", es10.2, " u_mag = ", es10.2)') i, c_mag, u_mag
!    call test_line(mesg, 4, Atest, Csys, u(i), u_mag*c_mag, tolerance=1.0e-15)
!  enddo

end subroutine end_value_h4


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
!! computed, the tridiagonal system is built, boundary conditions are prescribed and
!! the system is solved to yield edge-value estimates.  This scheme is described in detail
!! by White and Adcroft, 2009, J. Comp. Phys, https://doi.org/10.1016/j.jcp.2008.04.026
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
  real :: h0, h1, h2, h3       ! cell widths [H]
  real :: hMin                 ! The minimum thickness used in these calculations [H]
  real :: h01, h01_2, h01_3, h01_4, h01_5, h01_6  ! Summed thicknesses to various powers [H^n ~> m^n or kg^n m-2n]
  real :: h23, h23_2, h23_3, h23_4, h23_5, h23_6  ! Summed thicknesses to various powers [H^n ~> m^n or kg^n m-2n]
  real :: hNeglect             ! A negligible thickness [H].
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
  integer :: i, j, k              ! loop indexes

  hNeglect = hNeglect_edge_dflt ; if (present(h_neglect)) hNeglect = h_neglect

  ! Loop on interior cells
  do k = 2,N-2
    ! Store temporary cell widths, avoiding singularities from zero thicknesses or extreme changes.
    hMin = max(hNeglect, hMinFrac*((h(k-1) + h(k)) + (h(k+1) + h(k+2))))
    h0 = max(h(k-1), hMin) ; h1 = max(h(k), hMin)
    h2 = max(h(k+1), hMin) ; h3 = max(h(k+2), hMin)

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

    h01   = h0 + h1
    h01_2 = h01 * h01
    h01_3 = h01 * h01_2
    h01_4 = h01_2 * h01_2
    h01_5 = h01_4 * h01
    h01_6 = h01_3 * h01_3

    d2 = ( h01_2 - h1_2 ) / h0 ! = (2.0*h1 + h0)
    d3 = ( h01_3 - h1_3 ) / h0 ! = (3.0*h1_2 + h0*(3.0*h1 + h0))
    d4 = ( h01_4 - h1_4 ) / h0 ! = (4.0*h1_3 + h0*(6.0*h1_2 + h0*(4.0*h1 + h0)))
    d5 = ( h01_5 - h1_5 ) / h0 ! = (5.0*h1_4 + h0*(10.0*h1_3 + h0*(10.0*h1_2 + h0*(5.0*h1 + h0))))
    d6 = ( h01_6 - h1_6 ) / h0 ! = (6.0*h1_5 + h0*(15.0*h1_4 + h0*(20.0*h1_3 + h0*(15.0*h1_2 + h0*(6.0*h1 + h0)))))

    h23   = h2 + h3
    h23_2 = h23 * h23
    h23_3 = h23 * h23_2
    h23_4 = h23_2 * h23_2
    h23_5 = h23_4 * h23
    h23_6 = h23_3 * h23_3

    n2 = ( h23_2 - h2_2 ) / h3  !  = 2.0*h2   + h3
    n3 = ( h23_3 - h2_3 ) / h3  !  = 3.0*h2_2 + h3*(3.0*h2   + h3)
    n4 = ( h23_4 - h2_4 ) / h3  !  = 4.0*h2_3 + h3*(6.0*h2_2 + h3*(4.0*h2 + h3))
    n5 = ( h23_5 - h2_5 ) / h3  !  = 5.0*h2_4 + h3*(10.0*h2_3 + h3*(10.0*h2_2 + h3*(5.0*h2 + h3)))
    n6 = ( h23_6 - h2_6 ) / h3  !  = 6.0*h2_5 + h3*(15.0*h2_4 + h3*(20.0*h2_3 + h3*(15.0*h2_2 + h3*(6.0*h2 + h3))))

    ! Compute matrix entries as described in Eq. (48) of White and Adcroft (2009)
    Asys(1,1) = 1.0
    Asys(1,2) = 1.0
    Asys(1,3) = -1.0
    Asys(1,4) = -1.0
    Asys(1,5) = -1.0
    Asys(1,6) = -1.0

    Asys(2,1) = - h1
    Asys(2,2) = h2
    Asys(2,3) = 0.5 * d2
    Asys(2,4) = 0.5 * h1
    Asys(2,5) = -0.5 * h2
    Asys(2,6) = -0.5 * n2

    Asys(3,1) = 0.5 * h1_2
    Asys(3,2) = 0.5 * h2_2
    Asys(3,3) = -d3 / 6.0
    Asys(3,4) = - h1_2 / 6.0
    Asys(3,5) = - h2_2 / 6.0
    Asys(3,6) = - n3 / 6.0

    Asys(4,1) = - h1_3 / 6.0
    Asys(4,2) = h2_3 / 6.0
    Asys(4,3) = d4 / 24.0
    Asys(4,4) = h1_3 / 24.0
    Asys(4,5) = - h2_3 / 24.0
    Asys(4,6) = - n4 / 24.0

    Asys(5,1) = h1_4 / 24.0
    Asys(5,2) = h2_4 / 24.0
    Asys(5,3) = -d5 / 120.0
    Asys(5,4) = - h1_4 / 120.0
    Asys(5,5) = - h2_4 / 120.0
    Asys(5,6) = - n5 / 120.0

    Asys(6,1) = - h1_5 / 120.0
    Asys(6,2) = h2_5 / 120.0
    Asys(6,3) = d6 / 720.0
    Asys(6,4) = h1_5 / 720.0
    Asys(6,5) = - h2_5 / 720.0
    Asys(6,6) = - n6 / 720.0

    Bsys(:) = (/ -1.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

    call solve_linear_system( Asys, Bsys, Csys, 6, .false. )

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

  ! Use a right-biased stencil for the second row, as described in Eq. (49) of White and Adcroft (2009).

  ! Store temporary cell widths, avoiding singularities from zero thicknesses or extreme changes.
  hMin = max(hNeglect, hMinFrac*((h(1) + h(2)) + (h(3) + h(4))))
  h0 = max(h(1), hMin) ; h1 = max(h(2), hMin)
  h2 = max(h(3), hMin) ; h3 = max(h(4), hMin)

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

  h01   = h0 + h1 ;  h01_2 = h01 * h01 ;  h01_3 = h01 * h01_2
  h01_4 = h01_2 * h01_2 ; h01_5 = h01_4 * h01 ; h01_6 = h01_3 * h01_3

  d2 = ( h01_2 - h1_2 ) / h0 ! = (2.0*h1 + h0)
  d3 = ( h01_3 - h1_3 ) / h0 ! = (3.0*h1_2 + h0*(3.0*h1 + h0))
  d4 = ( h01_4 - h1_4 ) / h0 ! = (4.0*h1_3 + h0*(6.0*h1_2 + h0*(4.0*h1 + h0)))
  d5 = ( h01_5 - h1_5 ) / h0 ! = (5.0*h1_4 + h0*(10.0*h1_3 + h0*(10.0*h1_2 + h0*(5.0*h1 + h0))))
  d6 = ( h01_6 - h1_6 ) / h0 ! = (6.0*h1_5 + h0*(15.0*h1_4 + h0*(20.0*h1_3 + h0*(15.0*h1_2 + h0*(6.0*h1 + h0)))))

  h23   = h2 + h3
  h23_2 = h23 * h23
  h23_3 = h23 * h23_2
  h23_4 = h23_2 * h23_2
  h23_5 = h23_4 * h23
  h23_6 = h23_3 * h23_3

  n2 = ( h23_2 - h2_2 ) / h3  !  = 2.0*h2   + h3
  n3 = ( h23_3 - h2_3 ) / h3  !  = 3.0*h2_2 + h3*(3.0*h2   + h3)
  n4 = ( h23_4 - h2_4 ) / h3  !  = 4.0*h2_3 + h3*(6.0*h2_2 + h3*(4.0*h2 + h3))
  n5 = ( h23_5 - h2_5 ) / h3  !  = 5.0*h2_4 + h3*(10.0*h2_3 + h3*(10.0*h2_2 + h3*(5.0*h2 + h3)))
  n6 = ( h23_6 - h2_6 ) / h3  !  = 6.0*h2_5 + h3*(15.0*h2_4 + h3*(20.0*h2_3 + h3*(15.0*h2_2 + h3*(6.0*h2 + h3))))

  h0ph1   = h0 + h1
  h0ph1_2 = h0ph1 * h0ph1
  h0ph1_3 = h0ph1_2 * h0ph1
  h0ph1_4 = h0ph1_2 * h0ph1_2
  h0ph1_5 = h0ph1_3 * h0ph1_2

  ! Compute matrix entries
  Asys(1,1) = 1.0
  Asys(1,2) = 1.0
  Asys(1,3) = -1.0
  Asys(1,4) = -1.0
  Asys(1,5) = -1.0
  Asys(1,6) = -1.0

  Asys(2,1) = - h0ph1
  Asys(2,2) = 0.0
  Asys(2,3) = 0.5 * d2
  Asys(2,4) = 0.5 * h1
  Asys(2,5) = -0.5 * h2
  Asys(2,6) = -0.5 * n2

  Asys(3,1) = 0.5 * h0ph1_2
  Asys(3,2) = 0.0
  Asys(3,3) = -d3 / 6.0
  Asys(3,4) = - h1_2 / 6.0
  Asys(3,5) = - h2_2 / 6.0
  Asys(3,6) = - n3 / 6.0

  Asys(4,1) = - h0ph1_3 / 6.0
  Asys(4,2) = 0.0
  Asys(4,3) = d4 / 24.0
  Asys(4,4) = h1_3 / 24.0
  Asys(4,5) = - h2_3 / 24.0
  Asys(4,6) = - n4 / 24.0

  Asys(5,1) = h0ph1_4 / 24.0
  Asys(5,2) = 0.0
  Asys(5,3) = -d5 / 120.0
  Asys(5,4) = - h1_4 / 120.0
  Asys(5,5) = - h2_4 / 120.0
  Asys(5,6) = - n5 / 120.0

  Asys(6,1) = - h0ph1_5 / 120.0
  Asys(6,2) = 0.0
  Asys(6,3) = d6 / 720.0
  Asys(6,4) = h1_5 / 720.0
  Asys(6,5) = - h2_5 / 720.0
  Asys(6,6) = - n6 / 720.0

  Bsys(:) = (/ -1.0, h1, -0.5*h1_2, h1_3/6.0, -h1_4/24.0, h1_5/120.0 /)

  call solve_linear_system( Asys, Bsys, Csys, 6, .false. )

  alpha = Csys(1)
  beta  = Csys(2)

  tri_l(2) = alpha
  tri_d(2) = 1.0
  tri_u(2) = beta
  tri_b(2) = Csys(3) * u(1) + Csys(4) * u(2) + Csys(5) * u(3) + Csys(6) * u(4)

  ! Boundary conditions: left boundary
  hMin = max( hNeglect, hMinFrac*((h(1)+h(2)) + (h(5)+h(6)) + (h(3)+h(4))) )
  x(1) = 0.0
  do i = 1,6
    dx = max( hMin, h(i) )
    xavg = x(i) + 0.5*dx
    Asys(i,1) = 1.0
    Asys(i,2) = xavg
    Asys(i,3) = (xavg**2 + C1_12*dx**2)
    Asys(i,4) = xavg * (xavg**2 + 0.25*dx**2)
    Asys(i,5) = (xavg**4 + 0.5*xavg**2*dx**2 + 0.0125*dx**4)
    Asys(i,6) = xavg * (xavg**4 + C5_6*xavg**2*dx**2 + 0.0625*dx**4)
    Bsys(i) = u(i)
    x(i+1) = x(i) + dx
  enddo

  call solve_linear_system( Asys, Bsys, Csys, 6, .false. )

  tri_l(1) = 0.0
  tri_d(1) = 1.0
  tri_u(1) = 0.0
  tri_b(1) = evaluation_polynomial( Csys, 6, x(1) )        ! first edge value

  ! Use a left-biased stencil for the second to last row

  ! Store temporary cell widths, avoiding singularities from zero thicknesses or extreme changes.
  hMin = max(hNeglect, hMinFrac*((h(N-3) + h(N-2)) + (h(N-1) + h(N))))
  h0 = max(h(N-3), hMin) ; h1 = max(h(N-2), hMin)
  h2 = max(h(N-1), hMin) ; h3 = max(h(N), hMin)

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

  h01   = h0 + h1
  h01_2 = h01 * h01
  h01_3 = h01 * h01_2
  h01_4 = h01_2 * h01_2
  h01_5 = h01_4 * h01
  h01_6 = h01_3 * h01_3

  d2 = ( h01_2 - h1_2 ) / h0 ! = (2.0*h1 + h0)
  d3 = ( h01_3 - h1_3 ) / h0 ! = (3.0*h1_2 + h0*(3.0*h1 + h0))
  d4 = ( h01_4 - h1_4 ) / h0 ! = (4.0*h1_3 + h0*(6.0*h1_2 + h0*(4.0*h1 + h0)))
  d5 = ( h01_5 - h1_5 ) / h0 ! = (5.0*h1_4 + h0*(10.0*h1_3 + h0*(10.0*h1_2 + h0*(5.0*h1 + h0))))
  d6 = ( h01_6 - h1_6 ) / h0 ! = (6.0*h1_5 + h0*(15.0*h1_4 + h0*(20.0*h1_3 + h0*(15.0*h1_2 + h0*(6.0*h1 + h0)))))

  h23   = h2 + h3
  h23_2 = h23 * h23
  h23_3 = h23 * h23_2
  h23_4 = h23_2 * h23_2
  h23_5 = h23_4 * h23
  h23_6 = h23_3 * h23_3

  n2 = ( h23_2 - h2_2 ) / h3  !  = 2.0*h2   + h3
  n3 = ( h23_3 - h2_3 ) / h3  !  = 3.0*h2_2 + h3*(3.0*h2   + h3)
  n4 = ( h23_4 - h2_4 ) / h3  !  = 4.0*h2_3 + h3*(6.0*h2_2 + h3*(4.0*h2 + h3))
  n5 = ( h23_5 - h2_5 ) / h3  !  = 5.0*h2_4 + h3*(10.0*h2_3 + h3*(10.0*h2_2 + h3*(5.0*h2 + h3)))
  n6 = ( h23_6 - h2_6 ) / h3  !  = 6.0*h2_5 + h3*(15.0*h2_4 + h3*(20.0*h2_3 + h3*(15.0*h2_2 + h3*(6.0*h2 + h3))))

  h2ph3   = h2 + h3
  h2ph3_2 = h2ph3 * h2ph3
  h2ph3_3 = h2ph3_2 * h2ph3
  h2ph3_4 = h2ph3_2 * h2ph3_2
  h2ph3_5 = h2ph3_3 * h2ph3_2

  ! Compute matrix entries
  Asys(1,1) = 1.0
  Asys(1,2) = 1.0
  Asys(1,3) = -1.0
  Asys(1,4) = -1.0
  Asys(1,5) = -1.0
  Asys(1,6) = -1.0

  Asys(2,1) =   0.0
  Asys(2,2) = h2ph3
  Asys(2,3) =   0.5 * d2
  Asys(2,4) =   0.5 * h1
  Asys(2,5) =   -0.5 * h2
  Asys(2,6) =   -0.5 * n2

  Asys(3,1) =   0.0
  Asys(3,2) = 0.5 * h2ph3_2
  Asys(3,3) =   -d3 / 6.0
  Asys(3,4) =   - h1_2 / 6.0
  Asys(3,5) =   - h2_2 / 6.0
  Asys(3,6) =   - n3 / 6.0

  Asys(4,1) =   0.0
  Asys(4,2) = h2ph3_3 / 6.0
  Asys(4,3) =   d4 / 24.0
  Asys(4,4) =   h1_3 / 24.0
  Asys(4,5) =   - h2_3 / 24.0
  Asys(4,6) =   - n4 / 24.0

  Asys(5,1) =   0.0
  Asys(5,2) = h2ph3_4 / 24.0
  Asys(5,3) =   - d5 / 120.0
  Asys(5,4) =   - h1_4 / 120.0
  Asys(5,5) =   - h2_4 / 120.0
  Asys(5,6) =   - n5 / 120.0

  Asys(6,1) =   0.0
  Asys(6,2) = h2ph3_5 / 120.0
  Asys(6,3) =   d6 / 720.0
  Asys(6,4) =   h1_5 / 720.0
  Asys(6,5) =   - h2_5 / 720.0
  Asys(6,6) =   - n6 / 720.0

  Bsys(:) = (/ -1.0, -h2, -0.5*h2_2, -h2_3/6.0, -h2_4/24.0, -h2_5/120.0 /)

  call solve_linear_system( Asys, Bsys, Csys, 6, .false. )

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
  hMin = max( hNeglect, hMinFrac*(h(N-3) + h(N-2)) + ((h(N-1) + h(N)) + (h(N-5) + h(N-4))) )
  x(1) = 0.0
  do i = 1,6
    dx = max( hMin, h(N-6+i) )
    xavg = x(i) + 0.5 * dx
    Asys(i,1) = 1.0
    Asys(i,2) = xavg
    Asys(i,3) = (xavg**2 + C1_12*dx**2)
    Asys(i,4) = xavg * (xavg**2 + 0.25*dx**2)
    Asys(i,5) = (xavg**4 + 0.5*xavg**2*dx**2 + 0.0125*dx**4)
    Asys(i,6) = xavg * (xavg**4 + C5_6*xavg**2*dx**2 + 0.0625*dx**4)
    Bsys(i) = u(N-6+i)
    x(i+1) = x(i) + dx
  enddo

  call solve_linear_system( Asys, Bsys, Csys, 6, .false. )

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


! Verify that A*C = R to within roundoff.
subroutine test_line(msg, N, A, C, R, mag, tolerance)
  integer,            intent(in) :: N
  real, dimension(4), intent(in) :: A
  real, dimension(4), intent(in) :: C
  real,               intent(in) :: R
  real,               intent(in) :: mag  !< The magnitude of leading order terms in this line
  real, optional,     intent(in) :: tolerance
  character(len=*) :: msg

  real :: sum, sum_mag
  real :: tol
  character(len=128) :: mesg2
  integer :: i

  tol = 1.0e-12 ; if (present(tolerance)) tol = tolerance

  sum = 0.0 ; sum_mag = max(0.0,mag)

  do i=1,N
    sum = sum + A(i) * C(i)
    sum_mag = sum_mag + abs(A(i) * C(i))
  enddo

  if (abs(sum - R) > tol * (sum_mag + abs(R))) then
    write(mesg2, '(", Fractional error = ", es12.4,", sum = ", es12.4)') (sum - R) / (sum_mag + abs(R)), sum
    call MOM_error(FATAL, "Failed line test: "//trim(msg)//trim(mesg2))
  endif

end subroutine test_line


end module regrid_edge_values
