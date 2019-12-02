!> Edge value estimation for high-order resconstruction
module regrid_edge_values

! This file is part of MOM6. See LICENSE.md for the license.

use regrid_solvers, only : solve_linear_system, solve_tridiagonal_system
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

#undef __DO_SAFETY_CHECKS__

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
subroutine bound_edge_values( N, h, u, edge_val, h_neglect )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell average properties (size N) in arbitrary units [A]
  real, dimension(:,:), intent(inout) :: edge_val !< Potentially modified edge values [A]
  real,       optional, intent(in)    :: h_neglect !< A negligibly small width [H]
  ! Local variables
  integer       :: k            ! loop index
  integer       :: k0, k1, k2
  real          :: h_l, h_c, h_r ! Layer thicknesses [H]
  real          :: u_l, u_c, u_r ! Cell average properties [A]
  real          :: u0_l, u0_r    ! Edge values of properties [A]
  real          :: sigma_l, sigma_c, sigma_r    ! left, center and right
                                                ! van Leer slopes [A H-1]
  real          :: slope         ! retained PLM slope [A H-1]
  real          :: hNeglect      ! A negligible thickness [H].

  hNeglect = hNeglect_dflt ; if (present(h_neglect)) hNeglect = h_neglect

  ! Loop on cells to bound edge value
  do k = 1,N

    ! For the sake of bounding boundary edge values, the left neighbor
    ! of the left boundary cell is assumed to be the same as the left
    ! boundary cell and the right neighbor of the right boundary cell
    ! is assumed to be the same as the right boundary cell. This
    ! effectively makes boundary cells look like extrema.
    if ( k == 1 ) then
      k0 = 1
      k1 = 1
      k2 = 2
    elseif ( k == N ) then
      k0 = N-1
      k1 = N
      k2 = N
    else
      k0 = k-1
      k1 = k
      k2 = k+1
    endif

    ! All cells can now be treated equally
    h_l = h(k0)
    h_c = h(k1)
    h_r = h(k2)

    u_l = u(k0)
    u_c = u(k1)
    u_r = u(k2)

    u0_l = edge_val(k,1)
    u0_r = edge_val(k,2)

    sigma_l = 2.0 * ( u_c - u_l ) / ( h_c + hNeglect )
    sigma_c = 2.0 * ( u_r - u_l ) / ( h_l + 2.0*h_c + h_r + hNeglect )
    sigma_r = 2.0 * ( u_r - u_c ) / ( h_c + hNeglect )

    if ( (sigma_l * sigma_r) > 0.0 ) then
      slope = sign( min(abs(sigma_l),abs(sigma_c),abs(sigma_r)), sigma_c )
    else
      slope = 0.0
    endif

    ! The limiter must be used in the local coordinate system to each cell.
    ! Hence, we must multiply the slope by h1. The multiplication by 0.5 is
    ! simply a way to make it useable in the limiter (cfr White and Adcroft
    ! JCP 2008 Eqs 19 and 20)
    slope = slope * h_c * 0.5

    if ( (u_l-u0_l)*(u0_l-u_c) < 0.0 ) then
      u0_l = u_c - sign( min( abs(slope), abs(u0_l-u_c) ), slope )
    endif

    if ( (u_r-u0_r)*(u0_r-u_c) < 0.0 ) then
      u0_r = u_c + sign( min( abs(slope), abs(u0_r-u_c) ), slope )
    endif

    ! Finally bound by neighboring cell means in case of round off
    u0_l = max( min( u0_l, max(u_l, u_c) ), min(u_l, u_c) )
    u0_r = max( min( u0_r, max(u_r, u_c) ), min(u_r, u_c) )

    ! Store edge values
    edge_val(k,1) = u0_l
    edge_val(k,2) = u0_r

  enddo ! loop on interior edges

end subroutine bound_edge_values

!> Replace discontinuous collocated edge values with their average
!!
!! For each interior edge, check whether the edge values are discontinuous.
!! If so, compute the average and replace the edge values by the average.
subroutine average_discontinuous_edge_values( N, edge_val )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:,:), intent(inout) :: edge_val !< Edge values that may be modified
                                           !! the second index size is 2.
  ! Local variables
  integer       :: k            ! loop index
  real          :: u0_minus     ! left value at given edge
  real          :: u0_plus      ! right value at given edge
  real          :: u0_avg       ! avg value at given edge

  ! Loop on interior edges
  do k = 1,N-1

    ! Edge value on the left of the edge
    u0_minus = edge_val(k,2)

    ! Edge value on the right of the edge
    u0_plus  = edge_val(k+1,1)

    if ( u0_minus /= u0_plus ) then
      u0_avg = 0.5 * ( u0_minus + u0_plus )
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
  real, dimension(:),   intent(in)    :: u !< cell averages (size N) in arbitrary units [A]
  real, dimension(:,:), intent(inout) :: edge_val !< Cell edge values [A].
  ! Local variables
  integer       :: k            ! loop index
  real          :: u0_minus     ! left value at given edge [A]
  real          :: u0_plus      ! right value at given edge [A]
  real          :: um_minus     ! left cell average [A]
  real          :: um_plus      ! right cell average [A]
  real          :: u0_avg       ! avg value at given edge [A]

  ! Loop on interior cells
  do k = 1,N-1

    ! Edge value on the left of the edge
    u0_minus = edge_val(k,2)

    ! Edge value on the right of the edge
    u0_plus  = edge_val(k+1,1)

    ! Left cell average
    um_minus = u(k)

    ! Right cell average
    um_plus = u(k+1)

    if ( (u0_plus - u0_minus)*(um_plus - um_minus) < 0.0 ) then
      u0_avg = 0.5 * ( u0_minus + u0_plus )
      u0_avg = max( min( u0_avg, max(um_minus, um_plus) ), min(um_minus, um_plus) )
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
subroutine edge_values_explicit_h2( N, h, u, edge_val, h_neglect )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell average properties (size N) in arbitrary units [A]
  real, dimension(:,:), intent(inout) :: edge_val !< Returned edge values [A]; the second index size is 2.
  real,       optional, intent(in)    :: h_neglect !< A negligibly small width [H]
  ! Local variables
  integer   :: k        ! loop index
  real      :: h0, h1   ! cell widths [H]
  real      :: u0, u1   ! cell averages [A]
  real      :: hNeglect ! A negligible thickness [H]

  hNeglect = hNeglect_edge_dflt ; if (present(h_neglect)) hNeglect = h_neglect

  ! Loop on interior cells
  do k = 2,N

    h0 = h(k-1)
    h1 = h(k)

    ! Avoid singularities when h0+h1=0
    if (h0+h1==0.) then
      h0 = hNeglect
      h1 = hNeglect
    endif

    u0 = u(k-1)
    u1 = u(k)

    ! Compute left edge value
    edge_val(k,1) = ( u0*h1 + u1*h0 ) / ( h0 + h1 )

    ! Left edge value of the current cell is equal to right edge
    ! value of left cell
    edge_val(k-1,2) = edge_val(k,1)

  enddo ! end loop on interior cells

  ! Boundary edge values are simply equal to the boundary cell averages
  edge_val(1,1) = u(1)
  edge_val(N,2) = u(N)

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
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell average properties (size N) in arbitrary units [A]
  real, dimension(:,:), intent(inout) :: edge_val !< Returned edge values [A]; the second index size is 2.
  real,       optional, intent(in)    :: h_neglect !< A negligibly small width [H]
  logical,    optional, intent(in)    :: answers_2018 !< If true use older, less acccurate expressions.

  ! Local variables
  integer               :: i, j
  real                  :: u0, u1, u2, u3   ! temporary properties [A]
  real                  :: h0, h1, h2, h3   ! temporary thicknesses [H]
  real                  :: f1, f2, f3       ! auxiliary variables with various units
  real                  :: e                ! edge value
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
    if (h0+h1==0. .or. h1+h2==0. .or. h2+h3==0.) then
      f1 = max( hNeglect, h0+h1+h2+h3 )
      h0 = max( hMinFrac*f1, h(i-2) )
      h1 = max( hMinFrac*f1, h(i-1) )
      h2 = max( hMinFrac*f1, h(i) )
      h3 = max( hMinFrac*f1, h(i+1) )
    endif

    u0 = u(i-2)
    u1 = u(i-1)
    u2 = u(i)
    u3 = u(i+1)

    f1 = (h0+h1) * (h2+h3) / (h1+h2)
    f2 = u1 * h2 + u2 * h1
    f3 = 1.0 / (h0+h1+h2) + 1.0 / (h1+h2+h3)

    e = f1 * f2 * f3

    f1 = h2 * (h2+h3) / ( (h0+h1+h2)*(h0+h1) )
    f2 = u1*(h0+2.0*h1) - u0*h1

    e = e + f1*f2

    f1 = h1 * (h0+h1) / ( (h1+h2+h3)*(h2+h3) )
    f2 = u2*(2.0*h2+h3) - u3*h2

    e = e + f1*f2

    e = e / ( h0 + h1 + h2 + h3)

    edge_val(i,1) = e
    edge_val(i-1,2) = e

#ifdef __DO_SAFETY_CHECKS__
    if (e /= e) then
      write(0,*) 'NaN in explicit_edge_h4 at k=',i
      write(0,*) 'u0-u3=',u0,u1,u2,u3
      write(0,*) 'h0-h3=',h0,h1,h2,h3
      write(0,*) 'f1-f3=',f1,f2,f3
      stop 'Nan during edge_values_explicit_h4'
    endif
#endif

  enddo ! end loop on interior cells

  ! Determine first two edge values
  f1 = max( hNeglect, hMinFrac*sum(h(1:4)) )
  x(1) = 0.0
  do i = 2,5
    x(i) = x(i-1) + max(f1, h(i-1))
  enddo

  do i = 1,4
    dx = max(f1, h(i) )
    if (use_2018_answers) then
      do j = 1,4 ; A(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / real(j) ; enddo
    else  ! Use expressions with less sensitivity to roundoff
      xavg = 0.5 * (x(i+1) + x(i))
      A(i,1) = dx
      A(i,2) = dx * xavg
      A(i,3) = dx * (xavg**2 + C1_12*dx**2)
      A(i,4) = dx * xavg * (xavg**2 + 0.25*dx**2)
    endif

    B(i) = u(i) * dx

  enddo

  call solve_linear_system( A, B, C, 4 )

  ! First edge value
  edge_val(1,1) = evaluation_polynomial( C, 4, x(1) )

  ! Second edge value
  edge_val(1,2) = evaluation_polynomial( C, 4, x(2) )
  edge_val(2,1) = edge_val(1,2)

#ifdef __DO_SAFETY_CHECKS__
  if (edge_val(1,1) /= edge_val(1,1) .or. edge_val(1,2) /= edge_val(1,2)) then
    write(0,*) 'NaN in explicit_edge_h4 at k=',1
    write(0,*) 'A=',A
    write(0,*) 'B=',B
    write(0,*) 'C=',C
    write(0,*) 'h(1:4)=',h(1:4)
    write(0,*) 'x=',x
    stop 'Nan during edge_values_explicit_h4'
  endif
#endif

  ! Determine last two edge values
  f1 = max( hNeglect, hMinFrac*sum(h(N-3:N)) )
  x(1) = 0.0
  do i = 2,5
    x(i) = x(i-1) + max(f1, h(N-5+i))
  enddo

  do i = 1,4
    dx = max(f1, h(N-4+i) )
    if (use_2018_answers) then
      do j = 1,4 ; A(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / real(j) ; enddo
    else  ! Use expressions with less sensitivity to roundoff
      xavg = 0.5 * (x(i+1) + x(i))
      A(i,1) = dx
      A(i,2) = dx * xavg
      A(i,3) = dx * (xavg**2 + C1_12*dx**2)
      A(i,4) = dx * xavg * (xavg**2 + 0.25*dx**2)
    endif

    B(i) = u(N-4+i) * dx

  enddo

  call solve_linear_system( A, B, C, 4 )

  ! Last edge value
  edge_val(N,2) = evaluation_polynomial( C, 4, x(5) )

  ! Second to last edge value
  edge_val(N,1) = evaluation_polynomial( C, 4, x(4) )
  edge_val(N-1,2) = edge_val(N,1)

#ifdef __DO_SAFETY_CHECKS__
  if (edge_val(N,1) /= edge_val(N,1) .or. edge_val(N,2) /= edge_val(N,2)) then
    write(0,*) 'NaN in explicit_edge_h4 at k=',N
    write(0,*) 'A='
    do i = 1,4
      do j = 1,4
        A(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / real(j)
      enddo
      write(0,*) A(i,:)
      B(i) = u(N-4+i) * ( h(N-4+i) )
    enddo
    write(0,*) 'B=',B
    write(0,*) 'C=',C
    write(0,*) 'h(:N)=',h(N-3:N)
    write(0,*) 'x=',x
    stop 'Nan during edge_values_explicit_h4'
  endif
#endif

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
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell average properties (size N) in arbitrary units [A]
  real, dimension(:,:), intent(inout) :: edge_val !< Returned edge values [A]; the second index size is 2.
  real,       optional, intent(in)    :: h_neglect !< A negligibly small width [H]
  logical,    optional, intent(in)    :: answers_2018 !< If true use older, less acccurate expressions.

  ! Local variables
  integer               :: i, j                 ! loop indexes
  real                  :: h0, h1               ! cell widths [H]
  real                  :: h0_2, h1_2, h0h1
  real                  :: d2, d4
  real                  :: alpha, beta          ! stencil coefficients
  real                  :: a, b
  real, dimension(5)    :: x                    ! Coordinate system with 0 at edges [H]
  real, parameter       :: C1_12 = 1.0 / 12.0
  real                  :: dx, xavg             ! Differences and averages of successive values of x [H]
  real, dimension(4,4)  :: Asys                 ! boundary conditions
  real, dimension(4)    :: Bsys, Csys
  real, dimension(N+1)  :: tri_l, &             ! trid. system (lower diagonal)
                           tri_d, &             ! trid. system (middle diagonal)
                           tri_u, &             ! trid. system (upper diagonal)
                           tri_b, &             ! trid. system (unknowns vector)
                           tri_x                ! trid. system (rhs)
  real      :: hNeglect          ! A negligible thickness [H]
  logical   :: use_2018_answers  ! If true use older, less acccurate expressions.

  use_2018_answers = .true. ; if (present(answers_2018)) use_2018_answers = answers_2018
  hNeglect = hNeglect_edge_dflt ; if (present(h_neglect)) hNeglect = h_neglect

  ! Loop on cells (except last one)
  do i = 1,N-1

    ! Get cell widths
    h0 = h(i)
    h1 = h(i+1)

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

    tri_l(i+1) = alpha
    tri_d(i+1) = 1.0
    tri_u(i+1) = beta

    tri_b(i+1) = a * u(i) + b * u(i+1)

  enddo ! end loop on cells

  ! Boundary conditions: left boundary
  h0 = max( hNeglect, hMinFrac*sum(h(1:4)) )
  x(1) = 0.0
  do i = 2,5
    x(i) = x(i-1) + max( h0, h(i-1) )
  enddo

  do i = 1,4
    dx = max(h0, h(i) )
    if (use_2018_answers) then
      do j = 1,4 ; Asys(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / j ; enddo
    else  ! Use expressions with less sensitivity to roundoff
      xavg = 0.5 * (x(i+1) + x(i))
      Asys(i,1) = dx
      Asys(i,2) = dx * xavg
      Asys(i,3) = dx * (xavg**2 + C1_12*dx**2)
      Asys(i,4) = dx * xavg * (xavg**2 + 0.25*dx**2)
    endif

    Bsys(i) = u(i) * dx

  enddo

  call solve_linear_system( Asys, Bsys, Csys, 4 )

  tri_d(1) = 1.0
  tri_u(1) = 0.0
  tri_b(1) = evaluation_polynomial( Csys, 4, x(1) )        ! first edge value

  ! Boundary conditions: right boundary
  h0 = max( hNeglect, hMinFrac*sum(h(N-3:N)) )
  x(1) = 0.0
  do i = 2,5
    x(i) = x(i-1) + max( h0, h(N-5+i) )
  enddo

  do i = 1,4
    dx = max(h0, h(N-4+i) )
    if (use_2018_answers) then
      do j = 1,4 ; Asys(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / j ; enddo
    else  ! Use expressions with less sensitivity to roundoff
      xavg = 0.5 * (x(i+1) + x(i))
      Asys(i,1) = dx
      Asys(i,2) = dx * xavg
      Asys(i,3) = dx * (xavg**2 + C1_12*dx**2)
      Asys(i,4) = dx * xavg * (xavg**2 + 0.25*dx**2)
    endif
    Bsys(i) = u(N-4+i) * dx

  enddo

  call solve_linear_system( Asys, Bsys, Csys, 4 )

  tri_l(N+1) = 0.0
  tri_d(N+1) = 1.0
  tri_b(N+1) = evaluation_polynomial( Csys, 4, x(5) )      ! last edge value

  ! Solve tridiagonal system and assign edge values
  call solve_tridiagonal_system( tri_l, tri_d, tri_u, tri_b, tri_x, N+1 )

  do i = 2,N
    edge_val(i,1)   = tri_x(i)
    edge_val(i-1,2) = tri_x(i)
  enddo
  edge_val(1,1) = tri_x(1)
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
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell average properties (size N) in arbitrary units [A]
  real, dimension(:,:), intent(inout) :: edge_val  !< Returned edge values [A]; the second index size is 2.
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
