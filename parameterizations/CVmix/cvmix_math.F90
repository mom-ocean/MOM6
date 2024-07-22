!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module cvmix_math

!BOP
!\newpage
! !MODULE: cvmix_math
!
! !AUTHOR:
!  Michael N. Levy, NCAR (mlevy@ucar.edu)
!
! !DESCRIPTION:
!  This module contains routines to compute polynomial interpolations (linear,
!  quadratic, or cubic spline), evaluate  third-order polynomials and their
!  derivatives at specific values, and compute roots of these polynomials.
!\\
!\\
!
! !REVISION HISTORY:
!  $Id$
!  $URL$

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                                 &
                                    cvmix_one

!EOP

  implicit none
  private
  save

!BOP

! !DEFINED PARAMETERS:
  integer, parameter, public :: CVMIX_MATH_INTERP_LINEAR      = 1
  integer, parameter, public :: CVMIX_MATH_INTERP_QUAD        = 2
  integer, parameter, public :: CVMIX_MATH_INTERP_CUBE_SPLINE = 3

  real(cvmix_r8), parameter :: CVMIX_MATH_NEWTON_TOL       = 1.0e-12_cvmix_r8
  integer,        parameter :: CVMIX_MATH_MAX_NEWTON_ITERS = 100

! !PUBLIC MEMBER FUNCTIONS:

  public :: cvmix_math_poly_interp
  public :: cvmix_math_cubic_root_find
  public :: cvmix_math_evaluate_cubic

!EOP

  contains

!BOP

! !IROUTINE: cvmix_math_poly_interp
! !INTERFACE:

  subroutine cvmix_math_poly_interp(coeffs, interp_type, x, y, x0, y0)

! !DESCRIPTION:
!  Given (x(1), y(1)), (x(2), y(2)), and possibly (x0, y0), compute coeffs =
!  $(/a_0, a_1, a_2, a_3/)$ such that, for $f(x) = \sum a_nx^n$, the following
!  hold: $f(x(1)) = y(1)$ and $f(x(2)) = y(2)$. For both quadratic and cubic
!  interpolation, $f'(x(1)) = (y(1)-y0)/(x(1)-x0)$ as well, and for cubic splines
!  $f'(x(2)) = (y(2) - y(1))/(x(2) - x(1))$.
!  \\
!  \\

! !INPUT PARAMETERS:
    integer,                      intent(in)    :: interp_type
    real(cvmix_r8), dimension(2), intent(in)    :: x, y
    real(cvmix_r8), optional,     intent(in)    :: x0, y0
! !OUTPUT PARAMETERS:
    real(cvmix_r8), dimension(4), intent(inout) :: coeffs

!EOP
!BOC

    ! Local variables
    real(cvmix_r8) :: det
    integer        :: k, k2
    real(kind=cvmix_r8), dimension(4,4) :: Minv
    real(kind=cvmix_r8), dimension(4)   :: rhs

    ! All interpolation assumes form of
    ! y = dx^3 + cx^2 + bx + a
    ! linear => c = d = 0
    ! quad   => d = 0
    coeffs(1:4) = 0.0_cvmix_r8
    select case (interp_type)
      case (CVMIX_MATH_INTERP_LINEAR)
        ! Match y(1) and y(2)
!        print*, "Linear interpolation"
        coeffs(2) = (y(2)-y(1))/(x(2)-x(1))
        coeffs(1) = y(1)-coeffs(2)*x(1)
      case (CVMIX_MATH_INTERP_QUAD)
        ! Match y(1), y(2), and y'(1) [requires x(0)]
!        print*, "Quadratic interpolation"
        ! [ x2^2 x2 1 ][ c ]   [    y2 ]
        ! [ x1^2 x1 1 ][ b ] = [    y1 ]
        ! [  2x1  1 0 ][ a ]   [ slope ]
        !      ^^^
        !       M
        det = -((x(2)-x(1))**2)
        ! only using 3x3 block of Minv and first 3 elements of rhs
        rhs(1) = y(2)
        rhs(2) = y(1)
        if (present(x0).and.present(y0)) then
          rhs(3) = (y(1)-y0)/(x(1)-x0)
        else
          rhs(3) = 0.0_cvmix_r8
        end if

        Minv(1,1) = -cvmix_one/det
        Minv(1,2) = cvmix_one/det
        Minv(1,3) = -cvmix_one/(x(2)-x(1))
        Minv(2,1) = real(2, cvmix_r8)*x(1)/det
        Minv(2,2) = -real(2, cvmix_r8)*x(1)/det
        Minv(2,3) = (x(2)+x(1))/(x(2)-x(1))
        Minv(3,1) = -(x(1)**2)/det
        Minv(3,2) = x(2)*(real(2, cvmix_r8)*x(1)-x(2))/det
        Minv(3,3) = -x(2)*x(1)/(x(2)-x(1))

        do k=1,3
          do k2=1,3
            ! Note: weird "4-k2" term is used because I switched from
            ! y= 0x^3 + bx^2 + cx + d to
            ! y = a + bx + cx^2 + 0x^3
            coeffs(k2) = coeffs(k2)+Minv(4-k2,k)*rhs(k)
          end do
        end do
      case (CVMIX_MATH_INTERP_CUBE_SPLINE)
        ! Match y(1), y(2), y'(1), and y'(2)
!        print*, "Cubic spline interpolation"
        ! [ x2^3 x2^2 x2 1 ][ d ]   [     y2 ]
        ! [ x1^3 x1^2 x1 1 ][ c ] = [     y1 ]
        ! [  3x1  2x1  1 0 ][ b ]   [ slope1 ]
        ! [  3x2  2x2  1 0 ][ a ]   [ slope2 ]
        !      ^^^
        !       M
        det = -((x(2)-x(1))**3)
        rhs(1) = y(2)
        rhs(2) = y(1)
        if (present(x0).and.present(y0)) then
          rhs(3) = (y(1)-y0)/(x(1)-x0)
        else
          rhs(3) = 0.0_cvmix_r8
        end if
        rhs(4) = (y(2)-y(1))/(x(2)-x(1))

        Minv(1,1) = real(2, cvmix_r8)/det
        Minv(1,2) = -real(2, cvmix_r8)/det
        Minv(1,3) = (x(1)-x(2))/det
        Minv(1,4) = (x(1)-x(2))/det
        Minv(2,1) = -real(3, cvmix_r8)*(x(2)+x(1))/det
        Minv(2,2) = real(3, cvmix_r8)*(x(2)+x(1))/det
        Minv(2,3) = (x(2)-x(1))*(real(2, cvmix_r8)*x(2)+x(1))/det
        Minv(2,4) = (x(2)-x(1))*(real(2, cvmix_r8)*x(1)+x(2))/det
        Minv(3,1) = real(6, cvmix_r8)*x(2)*x(1)/det
        Minv(3,2) = -real(6, cvmix_r8)*x(2)*x(1)/det
        Minv(3,3) = -x(2)*(x(2)-x(1))*(real(2, cvmix_r8)*x(1)+x(2))/det
        Minv(3,4) = -x(1)*(x(2)-x(1))*(real(2, cvmix_r8)*x(2)+x(1))/det
        Minv(4,1) = -(x(1)**2)*(real(3, cvmix_r8)*x(2)-x(1))/det
        Minv(4,2) = -(x(2)**2)*(-real(3, cvmix_r8)*x(1)+x(2))/det
        Minv(4,3) = x(1)*(x(2)**2)*(x(2)-x(1))/det
        Minv(4,4) = x(2)*(x(1)**2)*(x(2)-x(1))/det

        do k=1,4
          do k2=1,4
            ! Note: weird "5-k2" term is used because I switched from
            ! y = a + bx + cx^2 + dx^3 to
            ! y= ax^3 + bx^2 + cx + d
            coeffs(k2) = coeffs(k2)+Minv(5-k2,k)*rhs(k)
          end do
        end do
    end select

!EOC

  end subroutine cvmix_math_poly_interp

  function cvmix_math_cubic_root_find(coeffs, x0)

    real(cvmix_r8), dimension(4), intent(in) :: coeffs
    real(cvmix_r8),               intent(in) :: x0

    real(cvmix_r8) :: cvmix_math_cubic_root_find
    real(cvmix_r8) :: fun_val, root, slope
    integer :: it_cnt

    root = x0
    fun_val = coeffs(4)*(root**3)+coeffs(3)*(root**2)+coeffs(2)*root+coeffs(1)
    do it_cnt = 1, CVMIX_MATH_MAX_NEWTON_ITERS
      if (abs(fun_val).lt.CVMIX_MATH_NEWTON_TOL) &
        exit
      slope = 3.0_cvmix_r8*coeffs(4)*(root**2)+2.0_cvmix_r8*coeffs(3)*root+coeffs(2)
      root = root - fun_val/slope
      fun_val = coeffs(4)*(root**3)+coeffs(3)*(root**2)+coeffs(2)*root+coeffs(1)
    end do
    cvmix_math_cubic_root_find = root

  end function cvmix_math_cubic_root_find

!BOP

! !IROUTINE: cvmix_math_evaluate_cubic
! !INTERFACE:

  function cvmix_math_evaluate_cubic(coeffs, x_in, fprime)

! !DESCRIPTION:
!  Computes $f(x) = a_0 + a_1x + a_2x^2 + a_3x^3$ at $x = $\verb|x_in|, where
!  \verb|coeffs|$ = (/a_0, a_1, a_2, a_3/)$. If requested, can also return
!  $f'(x)$
!\\
!\\

! !USES:
!  Only those used by entire module.

! !INPUT PARAMETERS:
    real(cvmix_r8), dimension(4), intent(in) :: coeffs
    real(cvmix_r8),               intent(in) :: x_in

! !OUTPUT PARAMETERS:
    real(cvmix_r8) :: cvmix_math_evaluate_cubic
    real(cvmix_r8), optional, intent(out) :: fprime

!EOP
!BOC

    ! Local Variables
    integer :: i

    ! Initialize both the cubic and its derivative to its constant term and
    ! then add the powers of x_in via a do-loop. This both reduces the number
    ! of arithmetic steps in the algorithm and avoids possible compiler issues
    ! if x_in = 0 (because 0*0 is undefined in some compilers)
    cvmix_math_evaluate_cubic = coeffs(1)
    if (present(fprime)) &
      fprime = coeffs(2)
    do i=2,4
      cvmix_math_evaluate_cubic = cvmix_math_evaluate_cubic +                 &
                                  coeffs(i)*(x_in**(i-1))
      if (present(fprime).and.(i.gt.2)) &
        fprime = fprime + coeffs(i)*real(i-1,cvmix_r8)*(x_in**(i-2))
    end do

  end function cvmix_math_evaluate_cubic

end module cvmix_math
