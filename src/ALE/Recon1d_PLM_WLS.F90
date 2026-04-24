!> Piecewise Linear Method using Weighted Conservative Least Squares 1D reconstruction
module Recon1d_PLM_WLS

! This file is part of MOM6. See LICENSE.md for the license.

use Recon1d_type, only : Recon1d, testing

implicit none ; private

public PLM_WLS, testing

!> PLM reconstruction using Weighted Least Squares constrained to conserve for central cell
!!
!! The source for the methods ultimately used by this class are:
!! - init()                    *locally defined
!! - reconstruct()             *locally defined
!! - average()                 *locally defined
!! - f()                       *locally defined
!! - dfdx()                    *locally defined
!! - x()                    -> recon1d_type.x()
!! - check_reconstruction()    *locally defined
!! - unit_tests()              *locally defined
!! - destroy()                 *locally defined
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> init()
!! - reconstruct_parent()   -> reconstruct()
type, extends (Recon1d) :: PLM_WLS

  real, allocatable :: ul(:) !< Left edge value [A]
  real, allocatable :: ur(:) !< Right edge value [A]
  real, allocatable, private :: slp(:) !< Difference across cell, ur - ul [A].
                              !! This is redundant with ul and ur and not used
                              !! in any evaluations, but is needed for testing.

contains
  !> Implementation of the PLM_WLS initialization
  procedure :: init => init
  !> Implementation of the PLM_WLS reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of the PLM_WLS average over an interval [A]
  procedure :: average => average
  !> Implementation of evaluating the PLM_WLS reconstruction at a point [A]
  procedure :: f => f
  !> Implementation of the derivative of the PLM_WLS reconstruction at a point [A]
  procedure :: dfdx => dfdx
  !> Implementation of deallocation for PLM_WLS
  procedure :: destroy => destroy
  !> Implementation of check reconstruction for the PLM_WLS reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the PLM_WLS reconstruction
  procedure :: unit_tests => unit_tests

  !> Duplicate interface to init()
  procedure :: init_parent => init
  !> Duplicate interface to reconstruct()
  procedure :: reconstruct_parent => reconstruct

end type PLM_WLS

contains

!> Initialize a 1D PLM reconstruction for n cells
subroutine init(this, n, h_neglect, check)
  class(PLM_WLS),    intent(out) :: this      !< This reconstruction
  integer,           intent(in)  :: n         !< Number of cells in this column
  real, optional,    intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H]
  logical, optional, intent(in)  :: check     !< If true, enable some consistency checking

  this%n = n

  allocate( this%u_mean(n) )
  allocate( this%ul(n) )
  allocate( this%ur(n) )
  allocate( this%slp(n) )

  this%h_neglect = tiny( this%u_mean(1) )
  if (present(h_neglect)) this%h_neglect = h_neglect
  this%check = .false.
  if (present(check)) this%check = check

end subroutine init

!> Calculate a 1D PLM_WLS reconstruction based on h(:) and u(:)
subroutine reconstruct(this, h, u)
  class(PLM_WLS), intent(inout) :: this !< This reconstruction
  real,           intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,           intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
  real :: slp ! The PLM slopes (difference across cell) [A]
  real :: u_l, u_r, u_c ! Left, right, and center values [A]
  real :: h_l, h_c, h_r ! Thickness of left, center and right cells [H]
  real :: h_l0, h_r0 ! Thickness of left and right cells with h_neglect added [H]
  real :: hx2l, hx2r ! Contributions to denominator, <h x^2> [H3]
  real :: hxyl, hxyr ! Contributions to numerator, <h x y> [H2 A]
  integer :: n, km1, k, kp1

  n = this%n

  ! Loop over all cells
  do k = 1, n
    km1 = max(1, k-1)
    kp1 = min(n, k+1)
    u_l = u(km1)
    u_c = u(k)
    u_r = u(kp1)

    h_l = h(km1) * real( k - km1 ) ! This zeroes h_l at k==1
    h_c = h(k)
    h_r = h(kp1) * real( kp1 - k ) ! This zeroes h_r at k==n

    ! This is the slope that minimizes the error
    !  sum_l={-1,1} h(k+l) * [ u(k+l) - u(k) + slp * ( z(k+l) - z(k) ) ]
    ! i.e. volume weighted least squares
    h_l0 = h_l + this%h_neglect
    h_r0 = h_r + this%h_neglect
    hxyl = ( h_l * ( h_c + h_l ) ) * ( u_c - u_l )
    hxyr = ( h_r * ( h_c + h_r ) ) * ( u_r - u_c )
    hx2l = h_l0 * ( h_c + h_l0 )**2
    hx2r = h_r0 * ( h_c + h_r0 )**2
    slp = 2. * h_c * ( hxyr + hxyl ) / ( hx2l + hx2r )

    ! Mean value
    this%u_mean(k) = u_c

    ! Left edge
    this%ul(k) = u_c - 0.5 * slp

    ! Right edge
    this%ur(k) = u_c + 0.5 * slp

    ! Store slope
    this%slp(k) = slp
  enddo

end subroutine reconstruct

!> Value of PLM_WLS reconstruction at a point in cell k [A]
real function f(this, k, x)
  class(PLM_WLS), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: x    !< Non-dimensional position within element [nondim]
  real :: du ! Difference across cell [A]

  du = this%ur(k) - this%ul(k)

  ! This expression might be used beyond the element to evaluate
  ! LS errors. In other PLM implementations x is bounded to the
  ! element and the expressions are constructed to not exceed
  ! bounds. There are no such constraints for PLM_WLS.
  f = this%u_mean(k) + du * ( x - 0.5)
  !f = this%u_mean(k) + this%slp(k) * ( x - 0.5)

end function f

!> Derivative of PLM_WLS reconstruction at a point in cell k [A]
real function dfdx(this, k, x)
  class(PLM_WLS), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: x    !< Non-dimensional position within element [nondim]

  dfdx = this%ur(k) - this%ul(k)

end function dfdx

!> Average between xa and xb for cell k of a 1D PLM reconstruction [A]
real function average(this, k, xa, xb)
  class(PLM_WLS), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: xa   !< Start of averaging interval on element (0 to 1)
  real,           intent(in) :: xb   !< End of averaging interval on element (0 to 1)
  real :: xmab ! Mid-point between xa and xb (0 to 1)
  real :: u_a, u_b ! Values at xa and xb [A]

  ! Mid-point between xa and xb
  xmab = 0.5 * ( xa + xb )

  ! This expression for u_a can overshoot u_r but is good for xmab<<1
  u_a = this%ul(k) + ( this%ur(k)  - this%ul(k) ) * xmab
  ! This expression for u_b can overshoot u_l but is good for 1-xmab<<1
  u_b = this%ur(k) + ( this%ul(k)  - this%ur(k) ) * ( 1. - xmab )

  ! Since u_a and u_b are both bounded, this will perserve uniformity but will the
  ! sum be bounded? Emperically it seems to work...
  average = 0.5 * ( u_a + u_b )

end function average

!> Deallocate the PLM reconstruction
subroutine destroy(this)
  class(PLM_WLS), intent(inout) :: this !< This reconstruction

  deallocate( this%u_mean, this%ul, this%ur )

end subroutine destroy

!> Checks the PLM_WLS reconstruction for consistency
logical function check_reconstruction(this, h, u)
  class(PLM_WLS), intent(in) :: this !< This reconstruction
  real,           intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,           intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables
  integer :: k
  real :: slp ! Cell slope [A]
  type(PLM_WLS) :: perturbed !< A perturbed reconstruction
  real :: u_l, u_r, u_c ! Left, right, and center values [A]
  real :: h_l, h_c, h_r ! Thickness of left, center and right cells [H]
  real :: h_l0, h_r0, h_c0 ! Thickness of left, right, center cells with h_neglect added [H]
  real :: x_l, x_r ! Positions of left and right cells [H]
  real :: hx2l, hx2r ! Contributions to denominator, <h x^2> [H3]
  real :: hxyl, hxyr ! Contributions to numerator, <h x y> [H2 A]
  real :: hy2l, hy2r ! Contributions to error, <h y^2> [H3]
  real :: y_l, y_r ! Left, right, value differencess [A]
  real :: b_h, bp_h ! slp / h_c [A H-1]
  integer :: km1, kp1

  check_reconstruction = .false.

  do k = 1, this%n
    if ( abs( this%u_mean(k) - u(k) ) > 0. ) check_reconstruction = .true.
  enddo

  ! Check the cell reconstruction is monotonic within each cell (it should be as a straight line)
  do k = 1, this%n
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ur(k) - this%u_mean(k) ) < 0. ) check_reconstruction = .true.
  enddo

  ! Check the cell is a straight line (to within machine precision)
  do k = 1, this%n
    if ( abs(2. * this%u_mean(k) - ( this%ul(k) + this%ur(k) )) > epsilon(this%u_mean(1)) * &
         max(abs(2. * this%u_mean(k)), abs(this%ul(k)), abs(this%ur(k))) ) check_reconstruction = .true.
  enddo

  ! Create a perturbable reconstruction
  call perturbed%init( this%n, h_neglect=this%h_neglect )
  call perturbed%reconstruct( h, u ) ! Should reproduce "this"
  ! Check the copy is identical
  do k = 1, this%n
    if ( abs( perturbed%u_mean(k) - this%u_mean(k) ) > 0. ) check_reconstruction = .true.
    if ( abs( perturbed%ul(k) - this%ul(k) ) > 0. ) check_reconstruction = .true.
    if ( abs( perturbed%ur(k) - this%ur(k) ) > 0. ) check_reconstruction = .true.
    if ( abs( perturbed%slp(k) - this%slp(k) ) > 0. ) check_reconstruction = .true.
  enddo
  ! Now perturb the slope. The local error should not decrease.
  do k = 1, this%n
    slp = this%slp(k) * ( 1.0 + 1. * epsilon(slp) )
    perturbed%slp(k) = slp
    perturbed%ul(k) = u(k) - 0.5 * slp
    perturbed%ur(k) = u(k) + 0.5 * slp
    if ( LS_error(perturbed, k, h, u) < LS_error(this, k, h, u) ) check_reconstruction = .true.

    slp = this%slp(k) * ( 1.0 - 1. * epsilon(slp) )
    perturbed%slp(k) = slp
    perturbed%ul(k) = u(k) - 0.5 * slp
    perturbed%ur(k) = u(k) + 0.5 * slp
    if ( LS_error(perturbed, k, h, u) < LS_error(this, k, h, u) ) check_reconstruction = .true.
  enddo

end function check_reconstruction

!> Returns local least squares error for a particular cell
!!
!! Note that this is the error relative to the minimum of the loss function so that at the
!! true solution this function returns zero. See module documentation.
real function LS_error(this, k, h, u)
  type(PLM_WLS), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,           intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables
  real :: u_l, u_r, u_c ! Left, right, and center values [A]
  real :: h_l, h_c, h_r ! Thickness of left, center and right cells [H]
  real :: h_l0, h_r0, hc0 ! Thickness of left, right, center cells with h_neglect added [H]
  real :: hx2l, hx2r ! Contributions to denominator, <h x^2> [H3]
  real :: hxyl, hxyr ! Contributions to numerator, <h x y> [H2 A]
  integer :: km1, kp1

  km1 = max(1, k-1)
  kp1 = min(this%n, k+1)
  u_l = u(km1)
  u_c = u(k)
  u_r = u(kp1)

  h_l = h(km1) * real( k - km1 ) ! This zeroes h_l at k==1
  h_r = h(kp1) * real( kp1 - k ) ! This zeroes h_r at k==n
  h_c = h(k)
  hc0 = h_c + this%h_neglect

  h_l0 = h_l + this%h_neglect
  h_r0 = h_r + this%h_neglect
  hxyl = ( h_l * 0.5 * ( h_c + h_l ) ) * ( u_c - u_l )
  hxyr = ( h_r * 0.5 * ( h_c + h_r ) ) * ( u_r - u_c )
  hx2l = h_l0 * 0.25 * ( h_c + h_l0 )**2
  hx2r = h_r0 * 0.25 * ( h_c + h_r0 )**2
  LS_error = h_c * ( ( hx2l + hx2r ) * this%slp(k) - h(k) * ( hxyl + hxyr ) )**2
  LS_error = LS_error / ( hc0 * ( hx2l + hx2r ) )
end function LS_error

!> Runs PLM_WLS reconstruction unit tests and returns True for any fails, False otherwise
logical function unit_tests(this, verbose, stdout, stderr)
  class(PLM_WLS), intent(inout) :: this    !< This reconstruction
  logical,        intent(in)    :: verbose !< True, if verbose
  integer,        intent(in)    :: stdout  !< I/O channel for stdout
  integer,        intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables
  real, allocatable :: ul(:), ur(:), um(:) ! test values [A]
  real, allocatable :: ull(:), urr(:) ! test values [A]
  type(testing) :: test ! convenience functions
  integer :: k

  call test%set( stdout=stdout ) ! Sets the stdout channel in test
  call test%set( stderr=stderr ) ! Sets the stderr channel in test
  call test%set( verbose=verbose ) ! Sets the verbosity flag in test

  call this%init(3, h_neglect=1.e-20)
  call test%test( this%n /= 3, "Setting number of levels")
  allocate( um(3), ul(3), ur(3), ull(3), urr(3) )

  call this%reconstruct( (/1.,1.,1./), (/-1.,0.,2./) )
  call test%real_arr(3, this%slp, (/1.,1.5,2./), "(1,1,1)(-1,0,2) slope")

  do k = 1, 3
    um(k) = LS_error(this, k, (/1.,1.,1./), (/-1.,0.,2./) )
  enddo
  call test%real_arr(3, um, (/0.,0.,0./), "(1,1,1)(-1,0,2) LS' rel error")

  call this%reconstruct( (/0.,1.,1./), (/-1.,0.,2./) )
  call test%real_arr(3, this%slp, (/0.,2.,2./), "(0,1,1)(-1,0,2) slope")

  do k = 1, 3
    um(k) = LS_error(this, k, (/0.,1.,1./), (/-1.,0.,2./) )
  enddo
  call test%real_arr(3, um, (/0.,0.,0./), "(0,1,1)(-1,0,2) LS' rel error")

  call this%reconstruct( (/1.,1.,1./), (/-2.,0.,1./) )
  call test%real_arr(3, this%slp, (/2.,1.5,1./), "(1,1,1)(-2,0,1) slope")

  call this%reconstruct( (/1.,1.,0./), (/-2.,0.,1./) )
  call test%real_arr(3, this%slp, (/2.,2.,0./), "(1,1,0)(-2,0,1) slope")

  call this%destroy()
  call this%init(3) ! Reset to defaults

  ! Straight line data on uniform grid
  call this%reconstruct( (/2.,2.,2./), (/1.,3.,5./) )
  call test%real_arr(3, this%u_mean, (/1.,3.,5./), "Straight line data")

  do k = 1, 3
    ul(k) = this%f(k, 0.)
    um(k) = this%f(k, 0.5)
    ur(k) = this%f(k, 1.)
  enddo
  call test%real_arr(3, ul, (/0.,2.,4./), "Evaluation on left edge")
  call test%real_arr(3, um, (/1.,3.,5./), "Evaluation in center")
  call test%real_arr(3, ur, (/2.,4.,6./), "Evaluation on right edge")

  do k = 1, 3
    ul(k) = this%dfdx(k, 0.)
    um(k) = this%dfdx(k, 0.5)
    ur(k) = this%dfdx(k, 1.)
  enddo
  call test%real_arr(3, ul, (/2.,2.,2./), "dfdx on left edge")
  call test%real_arr(3, um, (/2.,2.,2./), "dfdx in center")
  call test%real_arr(3, ur, (/2.,2.,2./), "dfdx on right edge")

  do k = 1, 3
    um(k) = LS_error(this, k, (/2.,2.,2./), (/1.,3.,5./) )
  enddo
  call test%real_arr(3, um, (/0.,0.,0./), "Rel error is 0")

  do k = 1, 3
    um(k) = this%average(k, 0.5, 0.75) ! Average from x=0.5 to 0.75 in each cell
  enddo
  call test%real_arr(3, um, (/1.25,3.25,5.25/), "Return interval average")

  call test%real_scalar( this%x(1,0.), 0., 'f-1(1,0)=0')
  call test%real_scalar( this%x(1,1.), 0.5, 'f-1(1,1)=0.5')
  call test%real_scalar( this%x(1,3.), 1., 'f-1(1,3)=1')
  call test%real_scalar( this%x(2,1.), 0., 'f-1(2,1)=0')
  call test%real_scalar( this%x(2,3.), 0.5, 'f-1(2,3)=0.5')
  call test%real_scalar( this%x(2,5.), 1., 'f-1(2,5)=1')
  call test%real_scalar( this%x(3,3.), 0., 'f-1(3,3)=0')
  call test%real_scalar( this%x(3,5.), 0.5, 'f-1(3,5)=0.5')
  call test%real_scalar( this%x(3,7.), 1., 'f-1(3,7)=1')

  call this%destroy()
  deallocate( um, ul, ur, ull, urr )

  allocate( um(4), ul(4), ur(4) )
  call this%init(4)

  deallocate( um, ul, ur )

  unit_tests = test%summarize("PLM_WLS:unit_tests")

end function unit_tests

!> \namespace recon1d_plm_wls
!!
!! This implementation of PLM fits the slope using least squares, but retains conservation
!! for the central cell by passing through the central value.
!! Cell-wise reconstructions are NOT limited by neighbours.
!! Thus, this reconstruction does not yield monotonic profiles needed for the general remapping problem.
!!
!! The algorithm solves the least squares problem of fitting a straight line through
!! the neighboring data. The line is constained to pass through the center cell,
!! \f$ (x_{k}, y_{k}) \f$, so that the construction is conservative. The more general
!! function \f$ f(x) = a_{k} + b_{k} x \f$ would not conserve for arbitrary data.
!!
!! The unknown parameter \f$ b_{k} \f$ in the line
!! \f[
!!  f(x) = y_{k} + b_{k} ( x - x_{k} ) / h_{k}
!! \f]
!! is fit to neighbors \f$ x_{k-1}, y_{k-1} \f$ and \f$ x_{k+1}, y_{k+1} \f$.
!!
!! Denoting \f$ y'_{k+j} = y_{k+j} - y_{k} \f$ and \f$ x'_{k+j} = x_{k+j} - x_{k} \f$
!! the local error is
!! \f{align}{
!!  e_{k+j} &= b_k \frac{ x_{k+j} - x_{k} }{ h_{k} } + y_{k} - y_{k+j} \\\\
!!          &= b_k \frac{ x'_{k+j} }{ h_{k} } - y'_{k+j}
!! \;\; . \f}
!!
!! We use volume weighting in the loss
!! \f[
!!  G(b) = h_{k-1} e_{k-1}^2 + h_{k+1} e_{k+1}^2
!! \;\; . \f]
!!
!! When solving for \f$ b_k \f$, we solve \f$ dG/db = 0 \f$ where
!! \f{align}{
!!  dG/db &= 2 h_{k-1} e_{k-1} \frac{ de_{k-1} }{db} + 2 h_{k+1} e_{k+1} \frac{ de_{k+1} }{db} \\\\
!!        &= 2 h_{k-1} ( b_k \frac{ x'_{k-1} }{ h_{k} } - \frac{ y'_{k-1} ) x'_{k-1} }{ h_{k} } +
!!           2 h_{k+1} ( b_k \frac{ x'_{k+1} }{ h_{k} } - \frac{ y'_{k+1} ) x'_{k+1} }{ h_{k} } \\\\
!!        &= 4 b_k \frac{ < h x'^2 > }{ h_{k}^2 } - 4 \frac{ < h x' y' > }{ h_{k} }
!! \f}
!! and where \f$ < a > = \frac{1}{2} ( a_{k-1} + a_{k+1} ) \f$.
!! Thus
!! \f[
!!  b_k = \frac{ h_{k} < h x' y' > }{ < h x'^2 > } \;\; .
!! \f]
!!
!! When evaluating the loss, \f$ G \f$, some rearrangement is necessary to reduce truncation
!! errors. Since
!! \f{align}{
!!  e_{k+j}^2 &= \left( b \frac{ x'_{k+j} }{ h_{k} } - y'_{k+j} \right)^2 \\\\
!!            &= b^2 \frac{ {x'}_{k+j}^2 }{ h_{k}^2 } - 2 b \frac{ x'_{k+j} y'_{k+j} }{ h_{k} } + {y'}_{k+j}^2
!! \f}
!! then
!! \f{align}{
!!  G(b) &= 2 < h e^2 > \\\\
!!       &= 2 b^2 \frac{ < h {x'}^2 > }{ h_{k}^2 } - 4 b \frac{ < h x' y' > }{ h_{k} } + 2 < h' {y'}^2 >
!! \;\; .
!! \f}
!!
!! If we denote the value of b that yields the minimum value as \f$ b^* \f$ then
!! \f[
!!  G(b^*) = < h {y'}^2 > - \frac{ < h x' y' >^2 }{ < h {x'}^2 > }
!! \;\; .
!! \f]
!!
!! Let
!! \f{align}{
!!  G''(b) &= G(b) - G(b^*) \\\\
!!         &= b^2 \frac{ < h {x'}^2 > }{ h_{k}^2 } - 2 b \frac{ < h x' y' > }{ h_{k} }
!!            + \frac{ < h x' y' > }{ < h {x'}^2 > } \\\\
!!         &= \frac{ \left( b < h {x'}^2 > - h_{k} < h x' y' > \right)^2 }{ h_{k} < h {x'}^2 > }
!! \;\; .
!! \f}
!! Minimizing \f$ G''(b) \f$ is equivalent to minimizing \f$ G(b) \f$ for the same data.
!! \f$ G''(b^*)=0 \f$ so evaluation with the last form, in the vicinity of \f$ b^* \f$, avoids
!! large cancelling terms.

end module Recon1d_PLM_WLS
