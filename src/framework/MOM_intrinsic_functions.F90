!> A module with intrinsic functions that are used by MOM but are not supported
!!  by some compilers.
module MOM_intrinsic_functions

! This file is part of MOM6. See LICENSE.md for the license.

use iso_fortran_env, only : stdout => output_unit, stderr => error_unit
use iso_fortran_env, only : int64, real64

implicit none ; private

public :: invcosh, cuberoot
public :: intrinsic_functions_unit_tests

! Floating point model, if bit layout from high to low is (sign, exp, frac)

integer, parameter :: bias = maxexponent(1.) - 1
  !< The double precision exponent offset
integer, parameter :: signbit = storage_size(1.) - 1
  !< Position of sign bit
integer, parameter :: explen = 1 + ceiling(log(real(bias))/log(2.))
  !< Bit size of exponent
integer, parameter :: expbit = signbit - explen
  !< Position of lowest exponent bit
integer, parameter :: fraclen = expbit
  !< Length of fractional part

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
  real :: ra_3 ! root_asx cubed [B3]
  real :: num ! The numerator of an expression for the evolving estimate of the cube root of asx
              ! in arbitrary units that can grow or shrink with each iteration [B C]
  real :: den ! The denominator of an expression for the evolving estimate of the cube root of asx
              ! in arbitrary units that can grow or shrink with each iteration [C]
  real :: num_prev ! The numerator of an expression for the previous iteration of the evolving estimate
              ! of the cube root of asx in arbitrary units that can grow or shrink with each iteration [B D]
  real :: np_3 ! num_prev cubed  [B3 D3]
  real :: den_prev ! The denominator of an expression for the previous iteration of the evolving estimate of
              ! the cube root of asx in arbitrary units that can grow or shrink with each iteration [D]
  real :: dp_3 ! den_prev cubed  [C3]
  real :: r0  ! Initial value of the iterative solver. [B C]
  real :: r0_3 ! r0 cubed [B3 C3]
  integer :: itt

  integer(kind=int64) :: e_x, s_x

  if ((x >= 0.0) .eqv. (x <= 0.0)) then
    ! Return 0 for an input of 0, or NaN for a NaN input.
    root = x
  else
    call rescale_cbrt(x, asx, e_x, s_x)

    !   Iteratively determine root_asx = asx**1/3 using Halley's method and then Newton's method,
    ! noting that Halley's method onverges monotonically and needs no bounding.  Halley's method is
    ! slightly more complicated that Newton's method, but converges in a third fewer iterations.
    !   Keeping the estimates in a fractional form Root = num / den allows this calculation with
    ! no real divisions during the iterations before doing a single real division at the end,
    ! and it is therefore more computationally efficient.

    ! This first estimate gives the same magnitude of errors for 0.125 and 1.0 after two iterations.
    ! The first iteration is applied explicitly.
    r0 = 0.707106
    r0_3 = r0 * r0 * r0
    num = r0 * (r0_3 + 2.0 * asx)
    den = 2.0 * r0_3 + asx

    do itt=1,2
      ! Halley's method iterates estimates as Root = Root * (Root**3 + 2.*asx) / (2.*Root**3 + asx).
      num_prev = num ; den_prev = den

      ! Pre-compute these as integer powers, to avoid `pow()`-like intrinsics.
      np_3 = num_prev * num_prev * num_prev
      dp_3 = den_prev * den_prev * den_prev

      num = num_prev * (np_3 + 2.0 * asx * dp_3)
      den = den_prev * (2.0 * np_3 + asx * dp_3)
      ! Equivalent to:  root_asx = root_asx * (root_asx**3 + 2.*asx) / (2.*root_asx**3 + asx)
    enddo
    ! At this point the error in root_asx is better than 1 part in 3e14.
    root_asx = num / den

    ! One final iteration with Newton's method polishes up the root and gives a solution
    ! that is within the last bit of the true solution.
    ra_3 = root_asx * root_asx * root_asx
    root_asx = root_asx - (ra_3 - asx) / (3.0 * (root_asx * root_asx))

    root = descale(root_asx, e_x, s_x)
  endif
end function cuberoot


!> Rescale `a` to the range [0.125, 1) and compute its cube-root exponent.
pure subroutine rescale_cbrt(a, x, e_r, s_a)
  real, intent(in) :: a
    !< The real parameter to be rescaled for cube root in abitrary units cubed [A3]
  real, intent(out) :: x
    !< The rescaled value of a in the range from 0.125 < asx <= 1.0, in ambiguous units cubed [B3]
  integer(kind=int64), intent(out) :: e_r
    !< Cube root of the exponent of the rescaling of `a`
  integer(kind=int64), intent(out) :: s_a
    !< The sign bit of a

  integer(kind=int64) :: xb
    ! Floating point value of a, bit-packed as an integer
  integer(kind=int64) :: e_a
    ! Unscaled exponent of a
  integer(kind=int64) :: e_x
    ! Exponent of x
  integer(kind=int64) :: e_div, e_mod
    ! Quotient and remainder of e in e = 3*(e/3) + modulo(e,3).

  ! Pack bits of a into xb and extract its exponent and sign.
  xb = transfer(a, 1_int64)
  s_a = ibits(xb, signbit, 1)
  e_a = ibits(xb, expbit, explen) - bias

  ! Compute terms of exponent decomposition e = 3*(e/3) + modulo(e,3).
  ! (Fortran division is round-to-zero, so we must emulate floor division.)
  e_mod = modulo(e_a, 3_int64)
  e_div = (e_a - e_mod)/3

  ! Our scaling decomposes e_a into e = {3*(e/3) + 3} + {modulo(e,3) - 3}.

  ! The first term is a perfect cube, whose cube root is computed below.
  e_r = e_div + 1

  ! The second term ensures that x is shifted to [0.125, 1).
  e_x = e_mod - 3

  ! Insert the new 11-bit exponent into xb and write to x and extend the
  ! bitcount to 12, so that the sign bit is zero and x is always positive.
  call mvbits(e_x + bias, 0, explen + 1, xb, fraclen)
  x = transfer(xb, 1.)
end subroutine rescale_cbrt


!> Undo the rescaling of a real number back to its original base.
pure function descale(x, e_a, s_a) result(a)
  real, intent(in) :: x
    !< The rescaled value which is to be restored in ambiguous units [B]
  integer(kind=int64), intent(in) :: e_a
    !< Exponent of the unscaled value
  integer(kind=int64), intent(in) :: s_a
    !< Sign bit of the unscaled value
  real :: a
    !< Restored value with the corrected exponent and sign in abitrary units [A]

  integer(kind=int64) :: xb
    ! Bit-packed real number into integer form
  integer(kind=int64) :: e_x
    ! Biased exponent of x

  ! Apply the corrected exponent and sign to x.
  xb = transfer(x, 1_int64)
  e_x = ibits(xb, expbit, explen)
  call mvbits(e_a + e_x, 0, explen, xb, expbit)
  call mvbits(s_a, 0, 1, xb, signbit)
  a = transfer(xb, 1.)
end function descale


!> Returns true if any unit test of intrinsic_functions fails, or false if they all pass.
function intrinsic_functions_unit_tests(verbose) result(fail)
  logical, intent(in) :: verbose !< If true, write results to stdout
  logical :: fail !< True if any of the unit tests fail

  ! Local variables
  real :: testval  ! A test value for self-consistency testing [nondim]
  logical :: v
  integer :: n

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
  fail = fail .or. Test_cuberoot(v, 1.0 - epsilon(1.0))
  fail = fail .or. Test_cuberoot(v, 1.0 - 0.5*epsilon(1.0))

  testval = 1.0e-99
  v = .false.
  do n=-160,160
    fail = fail .or. Test_cuberoot(v, testval)
    testval = (-2.908 * (1.414213562373 + 1.2345678901234e-5*n)) * testval
  enddo
end function intrinsic_functions_unit_tests

!> True if the cube of cuberoot(val) does not closely match val. False otherwise.
logical function Test_cuberoot(verbose, val)
  logical, intent(in) :: verbose !< If true, write results to stdout
  real, intent(in) :: val  !< The real value to test, in arbitrary units [A]
  ! Local variables
  real :: diff ! The difference between val and the cube root of its cube [A].

  diff = val - cuberoot(val)**3
  Test_cuberoot = (abs(diff) > 2.0e-15*abs(val))

  if (Test_cuberoot) then
    write(stdout, '("For val = ",ES22.15,", (val - cuberoot(val**3))) = ",ES9.2," <-- FAIL")') val, diff
  elseif (verbose) then
    write(stdout, '("For val = ",ES22.15,", (val - cuberoot(val**3))) = ",ES9.2)') val, diff

  endif
end function Test_cuberoot

end module MOM_intrinsic_functions
