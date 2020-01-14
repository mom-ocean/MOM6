!> Routines that estimate edge slopes to be used in
!! high-order reconstruction schemes.
module regrid_edge_slopes

! This file is part of MOM6. See LICENSE.md for the license.

use regrid_solvers, only : solve_linear_system, solve_tridiagonal_system
use regrid_solvers, only : solve_diag_dominant_tridiag, linear_solver
use polynomial_functions, only : evaluation_polynomial

implicit none ; private

public edge_slopes_implicit_h3
public edge_slopes_implicit_h5

! Specifying a dimensional parameter value, as is done here, is a terrible idea.
real, parameter :: hNeglect_dflt = 1.E-30 !< Default negligible cell thickness
real, parameter :: hMinFrac      = 1.e-5  !< A minimum fraction for min(h)/sum(h)

contains

!------------------------------------------------------------------------------
!> Compute ih3 edge slopes (implicit third order accurate)
!! in the same units as h.
!!
!! Compute edge slopes based on third-order implicit estimates. Note that
!! the estimates are fourth-order accurate on uniform grids
!!
!! Third-order implicit estimates of edge slopes are based on a two-cell
!! stencil. A tridiagonal system is set up and is based on expressing the
!! edge slopes in terms of neighboring cell averages. The generic
!! relationship is
!!
!! \f[
!! \alpha u'_{i-1/2} + u'_{i+1/2} + \beta u'_{i+3/2} =
!! a \bar{u}_i + b \bar{u}_{i+1}
!! \f]
!!
!! and the stencil looks like this
!!
!!          i     i+1
!!   ..--o------o------o--..
!!     i-1/2  i+1/2  i+3/2
!!
!! In this routine, the coefficients \f$\alpha\f$, \f$\beta\f$, a and b are computed,
!! the tridiagonal system is built, boundary conditions are prescribed and
!! the system is solved to yield edge-slope estimates.
!!
!! There are N+1 unknowns and we are able to write N-1 equations. The
!! boundary conditions close the system.
subroutine edge_slopes_implicit_h3( N, h, u, edge_slopes, h_neglect, answers_2018 )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< cell widths [H]
  real, dimension(N),   intent(in)    :: u !< cell average properties in arbitrary units [A]
  real, dimension(N,2), intent(inout) :: edge_slopes !< Returned edge slopes [A H-1]; the
                                           !! second index is for the two edges of each cell.
  real,       optional, intent(in)    :: h_neglect !< A negligibly small width [H]
  logical,    optional, intent(in)    :: answers_2018 !< If true use older, less acccurate expressions.
  ! Local variables
  integer               :: i, j                 ! loop indexes
  real                  :: h0, h1               ! cell widths [H or nondim]
  real                  :: h0_2, h1_2, h0h1     ! products of cell widths [H2 or nondim]
  real                  :: h0_3, h1_3           ! products of three cell widths [H3 or nondim]
  real                  :: h_min                ! A minimal cell width [H]
  real                  :: d                    ! A temporary variable [H3]
  real                  :: I_d                  ! A temporary variable [nondim]
  real                  :: I_h                  ! Inverses of thicknesses [H-1]
  real                  :: alpha, beta          ! stencil coefficients [nondim]
  real                  :: a, b                 ! weights of cells [H-1]
  real, parameter       :: C1_12 = 1.0 / 12.0
  real, dimension(5)    :: x          ! Coordinate system with 0 at edges [H]
  real                  :: dx, xavg   ! Differences and averages of successive values of x [H]
  real, dimension(4,4)  :: Asys       ! matrix used to find boundary conditions
  real, dimension(4)    :: Bsys, Csys
  real, dimension(3)    :: Dsys
  real, dimension(N+1)  :: tri_l, &     ! tridiagonal system (lower diagonal) [nondim]
                           tri_d, &     ! tridiagonal system (middle diagonal) [nondim]
                           tri_c, &     ! tridiagonal system central value, with tri_d = tri_c+tri_l+tri_u
                           tri_u, &     ! tridiagonal system (upper diagonal) [nondim]
                           tri_b, &     ! tridiagonal system (right hand side) [A H-1]
                           tri_x        ! tridiagonal system (solution vector) [A H-1]
  real      :: hNeglect  ! A negligible thickness [H].
  real      :: hNeglect3 ! hNeglect^3 [H3].
  logical   :: use_2018_answers  ! If true use older, less acccurate expressions.

  hNeglect = hNeglect_dflt ; if (present(h_neglect))  hNeglect = h_neglect
  hNeglect3 = hNeglect**3
  use_2018_answers = .true. ; if (present(answers_2018)) use_2018_answers = answers_2018

  ! Loop on cells (except last one)
  do i = 1,N-1

    if (use_2018_answers) then
      ! Get cell widths
      h0 = h(i)
      h1 = h(i+1)

      ! Auxiliary calculations
      h0h1 = h0 * h1
      h0_2 = h0 * h0
      h1_2 = h1 * h1
      h0_3 = h0_2 * h0
      h1_3 = h1_2 * h1

      d = 4.0 * h0h1 * ( h0 + h1 ) + h1_3 + h0_3

      ! Coefficients
      alpha = h1 * (h0_2 + h0h1 - h1_2) / ( d + hNeglect3 )
      beta  = h0 * (h1_2 + h0h1 - h0_2) / ( d + hNeglect3 )
      a = -12.0 * h0h1 / ( d + hNeglect3 )
      b = -a

      tri_l(i+1) = alpha
      tri_d(i+1) = 1.0
      tri_u(i+1) = beta

      tri_b(i+1) = a * u(i) + b * u(i+1)
    else
      ! Get cell widths
      h0 = max(h(i), hNeglect)
      h1 = max(h(i+1), hNeglect)

      I_h = 1.0 / (h0 + h1)
      h0 = h0 * I_h ; h1 = h1 * I_h

      h0h1 = h0 * h1 ; h0_2 = h0 * h0 ; h1_2 = h1 * h1
      h0_3 = h0_2 * h0 ; h1_3 = h1_2 * h1

      I_d = 1.0 / (4.0 * h0h1 * ( h0 + h1 ) + h1_3 + h0_3) ! = 1 / ((h0 + h1)**3 + h0*h1*(h0 + h1))

      ! Set the tridiagonal coefficients
      tri_l(i+1) = (h1 * ((h0_2 + h0h1) - h1_2)) * I_d
      ! tri_d(i+1) = 1.0
      tri_c(i+1) = 2.0 * ((h0_2 + h1_2) * (h0 + h1)) * I_d
      tri_u(i+1) = (h0 * ((h1_2 + h0h1) - h0_2)) * I_d
      ! The following expressions have been simplified using the nondimensionalization above:
      ! I_d = 1.0 / (1.0 + h0h1)
      ! tri_l(i+1) = (h0h1 - h1_3) * I_d
      ! tri_c(i+1) = 2.0 * (h0_2 + h1_2) * I_d
      ! tri_u(i+1) = (h0h1 - h0_3) * I_d

      tri_b(i+1) = 12.0 * (h0h1 * I_d) * ((u(i+1) - u(i)) * I_h)
    endif

  enddo ! end loop on cells

  ! Boundary conditions: set the first edge slope
  if (use_2018_answers) then
    x(1) = 0.0
    do i = 1,4
      dx = h(i)
      x(i+1) = x(i) + dx
      do j = 1,4 ; Asys(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / j ; enddo
      Bsys(i) = u(i) * dx
    enddo

    call solve_linear_system( Asys, Bsys, Csys, 4 )

    Dsys(1) = Csys(2) ; Dsys(2) = 2.0 * Csys(3) ; Dsys(3) = 3.0 * Csys(4)
    tri_b(1) = evaluation_polynomial( Dsys, 3, x(1) )  ! Set the first edge slope
    tri_d(1) = 1.0
  else ! Use expressions with less sensitivity to roundoff
    h_min = max( hNeglect, hMinFrac * ((h(1) + h(2)) + (h(3) + h(4))) )
    x(1) = 0.0
    do i = 1,4
      dx = max(h_min, h(i) )
      x(i+1) = x(i) + dx
      xavg = x(i) + 0.5*dx
      Asys(1,i) = 1.0
      Asys(2,i) = xavg
      Asys(3,i) = (xavg**2 + C1_12*dx**2)
      Asys(4,i) = xavg * (xavg**2 + 0.25*dx**2)
      Bsys(i) = u(i)
    enddo

    call linear_solver( 4, Asys, Bsys, Csys )

    ! Set the first edge slope
    tri_b(1) = Csys(2) ! + x(1)*(2.0*Csys(3) + x(1)*(3.0*Csys(4)))
    tri_c(1) = 1.0
  endif
  tri_u(1) = 0.0 ! tri_l(1) = 0.0

  ! Boundary conditions: set the last edge slope
  if (use_2018_answers) then
    x(1) = 0.0
    do i = 1,4
      dx = h(N-4+i)
      x(i+1) = x(i) + dx
      do j = 1,4 ; Asys(i,j) = ( (x(i+1)**j) - (x(i)**j) ) / j ; enddo
      Bsys(i) = u(N-4+i) * dx
    enddo

    call solve_linear_system( Asys, Bsys, Csys, 4 )

    Dsys(1) = Csys(2) ; Dsys(2) = 2.0 * Csys(3) ; Dsys(3) = 3.0 * Csys(4)
    ! Set the last edge slope
    tri_b(N+1) = evaluation_polynomial( Dsys, 3, x(5) )
    tri_d(N+1) = 1.0
  else
    ! Use expressions with less sensitivity to roundoff, including using a coordinate
    ! system that sets the origin at the last interface in the domain.
    h_min = max( hNeglect, hMinFrac * ((h(N-3) + h(N-2)) + (h(N-1) + h(N))) )
    x(1) = 0.0
    do i = 1,4
      dx = max(h_min, h(N+1-i) )
      x(i+1) = x(i) + dx
      xavg = x(i) + 0.5*dx
      Asys(1,i) = 1.0
      Asys(2,i) = xavg
      Asys(3,i) = (xavg**2 + C1_12*dx**2)
      Asys(4,i) = xavg * (xavg**2 + 0.25*dx**2)
      Bsys(i) = u(N+1-i)
    enddo

    call linear_solver( 4, Asys, Bsys, Csys )

    ! Set the last edge slope
    tri_b(N+1) = Csys(2)
    tri_c(N+1) = 1.0
  endif
  tri_l(N+1) = 0.0 ! tri_u(N+1) = 0.0

  ! Solve tridiagonal system and assign edge slopes
  if (use_2018_answers) then
    call solve_tridiagonal_system( tri_l, tri_d, tri_u, tri_b, tri_x, N+1 )
  else
    call solve_diag_dominant_tridiag( tri_l, tri_c, tri_u, tri_b, tri_x, N+1 )
  endif

  do i = 2,N
    edge_slopes(i,1)   = tri_x(i)
    edge_slopes(i-1,2) = tri_x(i)
  enddo
  edge_slopes(1,1) = tri_x(1)
  edge_slopes(N,2) = tri_x(N+1)

end subroutine edge_slopes_implicit_h3


!------------------------------------------------------------------------------
!> Compute ih5 edge values (implicit fifth order accurate)
subroutine edge_slopes_implicit_h5( N, h, u, edge_slopes, h_neglect, answers_2018 )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< cell widths [H]
  real, dimension(N),   intent(in)    :: u !< cell average properties in arbitrary units [A]
  real, dimension(N,2), intent(inout) :: edge_slopes !< Returned edge slopes [A H-1]; the
                                           !! second index is for the two edges of each cell.
  real, optional,       intent(in)    :: h_neglect !< A negligibly small width [H]
  logical,    optional, intent(in)    :: answers_2018 !< If true use older, less acccurate expressions.
! -----------------------------------------------------------------------------
! Fifth-order implicit estimates of edge values are based on a four-cell,
! three-edge stencil. A tridiagonal system is set up and is based on
! expressing the edge slopes in terms of neighboring cell averages.
!
! The generic relationship is
!
! \alpha u'_{i-1/2} + u'_{i+1/2} + \beta u'_{i+3/2} =
! a \bar{u}_{i-1} + b \bar{u}_i + c \bar{u}_{i+1} + d \bar{u}_{i+2}
!
! and the stencil looks like this
!
!         i-1     i     i+1    i+2
!   ..--o------o------o------o------o--..
!            i-1/2  i+1/2  i+3/2
!
! In this routine, the coefficients \alpha, \beta, a, b, c and d are
! computed, the tridiagonal system is built, boundary conditions are
! prescribed and the system is solved to yield edge-value estimates.
!
! Note that the centered stencil only applies to edges 3 to N-1 (edges are
! numbered 1 to n+1), which yields N-3 equations for N+1 unknowns. Two other
! equations are written by using a right-biased stencil for edge 2 and a
! left-biased stencil for edge N. The prescription of boundary conditions
! (using sixth-order polynomials) closes the system.
!
! CAUTION: For each edge, in order to determine the coefficients of the
!          implicit expression, a 6x6 linear system is solved. This may
!          become computationally expensive if regridding is carried out
!          often. Figuring out closed-form expressions for these coefficients
!          on nonuniform meshes turned out to be intractable.
! -----------------------------------------------------------------------------

  ! Local variables
  real :: h0, h1, h2, h3       ! cell widths [H]
  real :: hMin                 ! The minimum thickness used in these calculations [H]
  real :: h01, h01_2           ! Summed thicknesses to various powers [H^n ~> m^n or kg^n m-2n]
  real :: h23, h23_2           ! Summed thicknesses to various powers [H^n ~> m^n or kg^n m-2n]
  real :: hNeglect             ! A negligible thickness [H].
  real                  :: h1_2, h2_2           ! the coefficients of the
  real                  :: h1_3, h2_3           ! tridiagonal system
  real                  :: h1_4, h2_4           ! ...
  real                  :: h1_5, h2_5           ! ...
  real                  :: alpha, beta          ! stencil coefficients
  real, dimension(7)    :: x                    ! Coordinate system with 0 at edges [same units as h]
  real, parameter       :: C1_12 = 1.0 / 12.0
  real, parameter       :: C5_6 = 5.0 / 6.0
  real                  :: dx, xavg             ! Differences and averages of successive values of x [same units as h]
  real, dimension(6,6)  :: Asys                 ! matrix used to find  boundary conditions
  real, dimension(6)    :: Bsys, Csys           ! ...
  real, dimension(5)    :: Dsys                 ! derivative
  real, dimension(N+1)  :: tri_l, &             ! trid. system (lower diagonal)
                           tri_d, &             ! trid. system (middle diagonal)
                           tri_u, &             ! trid. system (upper diagonal)
                           tri_b, &             ! trid. system (unknowns vector)
                           tri_x                ! trid. system (rhs)
  real :: h_Min_Frac = 1.0e-4
  integer :: i, j, k              ! loop indexes

  hNeglect = hNeglect_dflt ; if (present(h_neglect)) hNeglect = h_neglect

  ! Loop on cells (except the first and last ones)
  do k = 2,N-2
    ! Store temporary cell widths, avoiding singularities from zero thicknesses or extreme changes.
    hMin = max(hNeglect, h_Min_Frac*((h(k-1) + h(k)) + (h(k+1) + h(k+2))))
    h0 = max(h(k-1), hMin) ; h1 = max(h(k), hMin)
    h2 = max(h(k+1), hMin) ; h3 = max(h(k+2), hMin)

    ! Auxiliary calculations
    h1_2 = h1 * h1 ; h1_3 = h1_2 * h1 ; h1_4 = h1_2 * h1_2 ; h1_5 = h1_3 * h1_2
    h2_2 = h2 * h2 ; h2_3 = h2_2 * h2 ; h2_4 = h2_2 * h2_2 ; h2_5 = h2_3 * h2_2

    ! Compute matrix entries as described in Eq. (52) of White and Adcroft (2009)
    Asys(1,1) = 0.0
    Asys(1,2) = 0.0
    Asys(1,3) = 1.0
    Asys(1,4) = 1.0
    Asys(1,5) = 1.0
    Asys(1,6) = 1.0

    Asys(2,1) = 2.0
    Asys(2,2) = 2.0
    Asys(2,3) = (2.0*h1 + h0)
    Asys(2,4) =  h1
    Asys(2,5) = -h2
    Asys(2,6) = -(2.0*h2   + h3)

    Asys(3,1) = 6.0*h1
    Asys(3,2) = -6.0* h2
    Asys(3,3) = (3.0*h1_2 + h0*(3.0*h1 + h0)) ! = ((h0+h1)**3 - h1**3) / h0
    Asys(3,4) = h1_2
    Asys(3,5) = h2_2
    Asys(3,6) = (3.0*h2_2 + h3*(3.0*h2 + h3)) ! = ((h2+h3)**3 - h2**3) / h3

    Asys(4,1) = -12.0* h1_2
    Asys(4,2) = -12.0* h2_2
    Asys(4,3) = -(4.0*h1_3 + h0*(6.0*h1_2 + h0*(4.0*h1 + h0))) ! = -((h0+h1)**4 - h1**4) / h0
    Asys(4,4) = - h1_3
    Asys(4,5) = h2_3
    Asys(4,6) = (4.0*h2_3 + h3*(6.0*h2_2 + h3*(4.0*h2 + h3))) ! = ((h2+h3)**4 - h2**4)/ h3

    Asys(5,1) = 20.0*h1_3
    Asys(5,2) = -20.0* h2_3
    Asys(5,3) = (5.0*h1_4 + h0*(10.0*h1_3 + h0*(10.0*h1_2 + h0*(5.0*h1 + h0))))
    Asys(5,4) = h1_4
    Asys(5,5) = h2_4
    Asys(5,6) = (5.0*h2_4 + h3*(10.0*h2_3 + h3*(10.0*h2_2 + h3*(5.0*h2 + h3))))

    Asys(6,1) = -30.0*h1_4
    Asys(6,2) = -30.0*h2_4
    Asys(6,3) = -(6.0*h1_5 + h0*(15.0*h1_4 + h0*(20.0*h1_3 + h0*(15.0*h1_2 + h0*(6.0*h1 + h0)))))
    Asys(6,4) = -h1_5
    Asys(6,5) = h2_5
    Asys(6,6) = (6.0*h2_5 + h3*(15.0*h2_4 + h3*(20.0*h2_3 + h3*(15.0*h2_2 + h3*(6.0*h2 + h3)))))

    Bsys(:) = (/ 0.0, -2.0, 0.0, 0.0, 0.0, 0.0 /)

    call solve_linear_system( Asys, Bsys, Csys, 6, .false. )

    alpha = Csys(1)
    beta  = Csys(2)

    tri_l(k+1) = alpha
    tri_d(k+1) = 1.0
    tri_u(k+1) = beta
    tri_b(k+1) = Csys(3) * u(k-1) + Csys(4) * u(k) + Csys(5) * u(k+1) + Csys(6) * u(k+2)

  enddo ! end loop on cells

  ! Use a right-biased stencil for the second row, as described in Eq. (53) of White and Adcroft (2009).

  ! Store temporary cell widths, avoiding singularities from zero thicknesses or extreme changes.
  hMin = max(hNeglect, h_Min_Frac*((h(1) + h(2)) + (h(3) + h(4))))
  h0 = max(h(1), hMin) ; h1 = max(h(2), hMin)
  h2 = max(h(3), hMin) ; h3 = max(h(4), hMin)

  ! Auxiliary calculations
  h1_2 = h1 * h1 ; h1_3 = h1_2 * h1 ; h1_4 = h1_2 * h1_2 ; h1_5 = h1_3 * h1_2
  h2_2 = h2 * h2 ; h2_3 = h2_2 * h2 ; h2_4 = h2_2 * h2_2 ; h2_5 = h2_3 * h2_2
  h01 = h0 + h1 ; h01_2 = h01 * h01

  ! Compute matrix entries
  Asys(1,1) = 0.0
  Asys(1,2) = 0.0
  Asys(1,3) = 1.0
  Asys(1,4) = 1.0
  Asys(1,5) = 1.0
  Asys(1,6) = 1.0

  Asys(2,1) = 2.0
  Asys(2,2) = 2.0
  Asys(2,3) = (2.0*h1 + h0)
  Asys(2,4) = h1
  Asys(2,5) = -h2
  Asys(2,6) = -(2.0*h2 + h3)

  Asys(3,1) = 6.0*h01
  Asys(3,2) = 0.0
  Asys(3,3) = (3.0*h1_2 + h0*(3.0*h1 + h0))
  Asys(3,4) = h1_2
  Asys(3,5) = h2_2
  Asys(3,6) = 3.0*h2_2 + h3*(3.0*h2 + h3)

  Asys(4,1) = -12.0*h01_2
  Asys(4,2) = 0.0
  Asys(4,3) = -(4.0*h1_3 + h0*(6.0*h1_2 + h0*(4.0*h1 + h0)))
  Asys(4,4) = -h1_3
  Asys(4,5) = h2_3
  Asys(4,6) = 4.0*h2_3 + h3*(6.0*h2_2 + h3*(4.0*h2 + h3))

  Asys(5,1) = 20.0*(h01*h01_2)
  Asys(5,2) = 0.0
  Asys(5,3) = (5.0*h1_4 + h0*(10.0*h1_3 + h0*(10.0*h1_2 + h0*(5.0*h1 + h0))))
  Asys(5,4) = h1_4
  Asys(5,5) = h2_4
  Asys(5,6) = 5.0*h2_4 + h3*(10.0*h2_3 + h3*(10.0*h2_2 + h3*(5.0*h2 + h3)))

  Asys(6,1) = -30.0*(h01_2*h01_2)
  Asys(6,2) = 0.0
  Asys(6,3) = -(6.0*h1_5 + h0*(15.0*h1_4 + h0*(20.0*h1_3 + h0*(15.0*h1_2 + h0*(6.0*h1 + h0)))))
  Asys(6,4) = -h1_5
  Asys(6,5) = h2_5
  Asys(6,6) = 6.0*h2_5 + h3*(15.0*h2_4 + h3*(20.0*h2_3 + h3*(15.0*h2_2 + h3*(6.0*h2 + h3))))

  Bsys(:) = (/ 0.0, -2.0, -6.0*h1, 12.0*h1_2, -20.0*h1_3, 30.0*h1_4 /)

  call solve_linear_system( Asys, Bsys, Csys, 6, .false. )

  alpha = Csys(1)
  beta  = Csys(2)

  tri_l(2) = alpha
  tri_d(2) = 1.0
  tri_u(2) = beta
  tri_b(2) = Csys(3) * u(1) + Csys(4) * u(2) + Csys(5) * u(3) + Csys(6) * u(4)

  ! Boundary conditions: left boundary
  x(1) = 0.0
  do i = 1,6
    dx = h(i)
    xavg = x(i) + 0.5 * dx
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

  Dsys(1) = Csys(2)
  Dsys(2) = 2.0 * Csys(3)
  Dsys(3) = 3.0 * Csys(4)
  Dsys(4) = 4.0 * Csys(5)
  Dsys(5) = 5.0 * Csys(6)

  tri_d(1) = 0.0
  tri_d(1) = 1.0
  tri_u(1) = 0.0
  tri_b(1) = evaluation_polynomial( Dsys, 5, x(1) )        ! first edge value

  ! Use a left-biased stencil for the second to last row, as described in Eq. (54) of White and Adcroft (2009).

  ! Store temporary cell widths, avoiding singularities from zero thicknesses or extreme changes.
  hMin = max(hNeglect, h_Min_Frac*((h(N-3) + h(N-2)) + (h(N-1) + h(N))))
  h0 = max(h(N-3), hMin) ; h1 = max(h(N-2), hMin)
  h2 = max(h(N-1), hMin) ; h3 = max(h(N), hMin)

  ! Auxiliary calculations
  h1_2 = h1 * h1 ; h1_3 = h1_2 * h1 ; h1_4 = h1_2 * h1_2 ; h1_5 = h1_3 * h1_2
  h2_2 = h2 * h2 ; h2_3 = h2_2 * h2 ; h2_4 = h2_2 * h2_2 ; h2_5 = h2_3 * h2_2

  h23 = h2 + h3 ; h23_2 = h23 * h23

  ! Compute matrix entries
  Asys(1,1) = 0.0
  Asys(1,2) = 0.0
  Asys(1,3) = 1.0
  Asys(1,4) = 1.0
  Asys(1,5) = 1.0
  Asys(1,6) = 1.0

  Asys(2,1) = 2.0
  Asys(2,2) = 2.0
  Asys(2,3) = (2.0*h1 + h0)
  Asys(2,4) = h1
  Asys(2,5) = -h2
  Asys(2,6) = -(2.0*h2 + h3)

  Asys(3,1) = 0.0
  Asys(3,2) = -6.0*h23
  Asys(3,3) = (3.0*h1_2 + h0*(3.0*h1 + h0))
  Asys(3,4) = h1_2
  Asys(3,5) = h2_2
  Asys(3,6) = 3.0*h2_2 + h3*(3.0*h2 + h3)

  Asys(4,1) = 0.0
  Asys(4,2) = -12.0*h23_2
  Asys(4,3) = -(4.0*h1_3 + h0*(6.0*h1_2 + h0*(4.0*h1 + h0)))
  Asys(4,4) = -h1_3
  Asys(4,5) = h2_3
  Asys(4,6) = 4.0*h2_3 + h3*(6.0*h2_2 + h3*(4.0*h2 + h3))

  Asys(5,1) = 0.0
  Asys(5,2) = -20.0*(h23*h23_2)
  Asys(5,3) = (5.0*h1_4 + h0*(10.0*h1_3 + h0*(10.0*h1_2 + h0*(5.0*h1 + h0))))
  Asys(5,4) = h1_4
  Asys(5,5) = h2_4
  Asys(5,6) = 5.0*h2_4 + h3*(10.0*h2_3 + h3*(10.0*h2_2 + h3*(5.0*h2 + h3)))

  Asys(6,1) = 0.0
  Asys(6,2) = -30.0*(h23_2*h23_2)
  Asys(6,3) = -(6.0*h1_5 + h0*(15.0*h1_4 + h0*(20.0*h1_3 + h0*(15.0*h1_2 + h0*(6.0*h1 + h0)))))
  Asys(6,4) = -h1_5
  Asys(6,5) = h2_5
  Asys(6,6) = 6.0*h2_5 + h3*(15.0*h2_4 + h3*(20.0*h2_3 + h3*(15.0*h2_2 + h3*(6.0*h2 + h3))))

  Bsys(:) = (/ 0.0, -2.0, 6.0*h2, 12.0*h2_2, 20.0*h2_3, 30.0*h2_4 /)

  call solve_linear_system( Asys, Bsys, Csys, 6, .false. )

  alpha = Csys(1)
  beta  = Csys(2)

  tri_l(N) = alpha
  tri_d(N) = 1.0
  tri_u(N) = beta
  tri_b(N) = Csys(3) * u(N-3) + Csys(4) * u(N-2) + Csys(5) * u(N-1) + Csys(6) * u(N)

  ! Boundary conditions: right boundary
  x(1) = 0.0
  do i = 1,6
    dx = h(N-6+i)
    xavg = x(i) + 0.5*dx
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

  Dsys(1) = Csys(2)
  Dsys(2) = 2.0 * Csys(3)
  Dsys(3) = 3.0 * Csys(4)
  Dsys(4) = 4.0 * Csys(5)
  Dsys(5) = 5.0 * Csys(6)

  tri_l(N+1) = 0.0
  tri_d(N+1) = 1.0
  tri_u(N+1) = 0.0
  tri_b(N+1) = evaluation_polynomial( Dsys, 5, x(7) )      ! last edge value

  ! Solve tridiagonal system and assign edge values
  call solve_tridiagonal_system( tri_l, tri_d, tri_u, tri_b, tri_x, N+1 )

  do i = 2,N
    edge_slopes(i,1)   = tri_x(i)
    edge_slopes(i-1,2) = tri_x(i)
  enddo
  edge_slopes(1,1) = tri_x(1)
  edge_slopes(N,2) = tri_x(N+1)

end subroutine edge_slopes_implicit_h5

end module regrid_edge_slopes
