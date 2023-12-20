!> A module with intrinsic functions that are used by MOM but are not supported
!!  by some compilers.
module MOM_intrinsic_functions

! This file is part of MOM6. See LICENSE.md for the license.

use iso_fortran_env, only : stdout => output_unit, stderr => error_unit

implicit none ; private

public :: invcosh, cuberoot
public :: intrinsic_functions_unit_tests

contains

!> Evaluate the inverse cosh, either using a math library or an
!! equivalent expression
function invcosh(x)
  real, intent(in) :: x !< The argument of the inverse of cosh [nondim].  NaNs will
                        !! occur if x<1, but there is no error checking
  real :: invcosh  ! The inverse of cosh of x [nondim]

#ifdef __INTEL_COMPILER
  invcosh = acosh(x)
#else
  invcosh = log(x+sqrt(x*x-1))
#endif

end function invcosh

!> Returns the cube root of a real argument at roundoff accuracy, in a form that works properly with
!! rescaling of the argument by integer powers of 8.  If the argument is a NaN, a NaN is returned.
elemental function cuberoot(x) result(root)
  real, intent(in) :: x !< The argument of cuberoot in arbitrary units cubed [A3]
  real :: root !< The real cube root of x in arbitrary units [A]

  real :: asx ! The absolute value of x rescaled by an integer power of 8 to put it into
              ! the range from 0.125 < asx <= 1.0, in ambiguous units cubed [B3]
  real :: root_asx ! The cube root of asx [B]
  real :: num ! The numerator of an expression for the evolving estimate of the cube root of asx
              ! in arbitrary units that can grow or shrink with each iteration [B C]
  real :: den ! The denominator of an expression for the evolving estimate of the cube root of asx
              ! in arbitrary units that can grow or shrink with each iteration [C]
  real :: num_prev ! The numerator of an expression for the previous iteration of the evolving estimate
              ! of the cube root of asx in arbitrary units that can grow or shrink with each iteration [B D]
  real :: den_prev ! The denominator of an expression for the previous iteration of the evolving estimate of
              ! the cube root of asx in arbitrary units that can grow or shrink with each iteration [D]
  real, parameter :: den_min = 2.**(minexponent(1.) / 4 + 4)  ! A value of den that triggers rescaling [C]
  real, parameter :: den_max = 2.**(maxexponent(1.) / 4 - 2)  ! A value of den that triggers rescaling [C]
  integer :: ex_3 ! One third of the exponent part of x, used to rescale x to get a.
  integer :: itt

  if ((x >= 0.0) .eqv. (x <= 0.0)) then
    ! Return 0 for an input of 0, or NaN for a NaN input.
    root = x
  else
    ex_3 = ceiling(exponent(x) / 3.)
    ! Here asx is in the range of 0.125 <= asx < 1.0
    asx = scale(abs(x), -3*ex_3)

    ! This first estimate is one iteration of Newton's method with a starting guess of 1.  It is
    ! always an over-estimate of the true solution, but it is a good approximation for asx near 1.
    num = 2.0 + asx
    den = 3.0
    ! Iteratively determine Root = asx**1/3 using Newton's method, noting that in this case Newton's
    ! method converges monotonically from above and needs no bounding.  For the range of asx from
    ! 0.125 to 1.0 with the first guess used above, 6 iterations suffice to converge to roundoff.
    do itt=1,9
      ! Newton's method iterates estimates as Root = Root - (Root**3 - asx) / (3.0 * Root**2), or
      ! equivalently as Root = (2.0*Root**2 + asx) / (3.0 * Root**2).
      ! Keeping the estimates in a fractional form Root = num / den allows this calculation with
      ! fewer (or no) real divisions during the iterations before doing a single real division
      ! at the end, and it is therefore more computationally efficient.

      num_prev = num ; den_prev = den
      num = 2.0 * num_prev**3 + asx * den_prev**3
      den = 3.0 * (den_prev * num_prev**2)

      if ((num * den_prev == num_prev * den) .or. (itt == 9)) then
        !   If successive estimates of root are identical, this is a converged solution.
        root_asx = num / den
        exit
      elseif (num * den_prev > num_prev * den) then
        !   If the estimates are increasing, this also indicates convergence, but for a more subtle
        ! reason.  Because Newton's method converges monotonically from above (at least for infinite
        ! precision math), the only reason why this estimate could increase is if the iterations
        ! have converged to a roundoff-level limit cycle around an irrational or otherwise
        ! unrepresentable solution, with values only changing in the last bit or two.  If so, we
        ! should stop iterating and accept the one of the current or previous solutions, both of
        ! which will be within numerical roundoff of the true solution.
        root_asx = num / den
        ! Pick the more accurate of the last two iterations.
        ! Given that both of the two previous iterations are within roundoff of the true
        ! solution, this next step might be overkill.
        if ( abs(den_prev**3*root_asx**3 - den_prev**3*asx) > abs(num_prev**3 - den_prev**3*asx) ) then
          ! The previous iteration was slightly more accurate, so use that for root_asx.
          root_asx = num_prev / den_prev
        endif
        exit
      endif

      ! Because successive estimates of the numerator and denominator tend to be the cube of their
      ! predecessors, the numerator and denominator need to be rescaled by division when they get
      ! too large or small to avoid overflow or underflow in the convergence test below.
      if ((den > den_max) .or. (den < den_min)) then
        num = scale(num, -exponent(den))
        den = scale(den, -exponent(den))
      endif

    enddo

    root = sign(scale(root_asx, ex_3), x)
  endif

end function cuberoot

!> Returns true if any unit test of intrinsic_functions fails, or false if they all pass.
logical function intrinsic_functions_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout

  ! Local variables
  real :: testval  ! A test value for self-consistency testing [nondim]
  logical :: fail, v

  fail = .false.
  v = verbose
  write(stdout,*) '==== MOM_intrinsic_functions: intrinsic_functions_unit_tests ==='

  fail = fail .or. Test_cuberoot(v, 1.2345678901234e9)
  fail = fail .or. Test_cuberoot(v, -9.8765432109876e-21)
  fail = fail .or. Test_cuberoot(v, 64.0)
  fail = fail .or. Test_cuberoot(v, -0.5000000000001)
  fail = fail .or. Test_cuberoot(v, 0.0)
  fail = fail .or. Test_cuberoot(v, 1.0)
  fail = fail .or. Test_cuberoot(v, 0.125)
  fail = fail .or. Test_cuberoot(v, 0.965)

end function intrinsic_functions_unit_tests

!> True if the cube of cuberoot(val) does not closely match val. False otherwise.
logical function Test_cuberoot(verbose, val)
  logical, intent(in) :: verbose !< If true, write results to stdout
  real, intent(in) :: val  !< The real value to test, in arbitrary units [A]
  ! Local variables
  real :: diff ! The difference between val and the cube root of its cube.

  diff = val - cuberoot(val**3)
  Test_cuberoot = (abs(diff) > 2.0e-15*abs(val))

  if (Test_cuberoot) then
    write(stdout, '("For val = ",ES22.15,", (val - cuberoot(val**3))) = ",ES9.2," <-- FAIL")') val, diff
  elseif (verbose) then
    write(stdout, '("For val = ",ES22.15,", (val - cuberoot(val**3))) = ",ES9.2)') val, diff

  endif
end function Test_cuberoot

end module MOM_intrinsic_functions
