!> A module with intrinsic functions that are used by MOM but are not supported
!!  by some compilers.
module MOM_intrinsic_functions

! This file is part of MOM6. See LICENSE.md for the license.

use iso_fortran_env, only : stdout => output_unit, stderr => error_unit
use iso_fortran_env, only : int64, real64

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
  integer :: itt

  integer(kind=int64) :: e_x, s_x

  if ((x >= 0.0) .eqv. (x <= 0.0)) then
    ! Return 0 for an input of 0, or NaN for a NaN input.
    root = x
  else
    call rescale_exp(x, asx, e_x, s_x)

    !   Iteratively determine root_asx = asx**1/3 using Halley's method and then Newton's method,
    ! noting that Halley's method onverges monotonically and needs no bounding.  Halley's method is
    ! slightly more complicated that Newton's method, but converges in a third fewer iterations.
    !   Keeping the estimates in a fractional form Root = num / den allows this calculation with
    ! no real divisions during the iterations before doing a single real division at the end,
    ! and it is therefore more computationally efficient.

    ! This first estimate gives the same magnitude of errors for 0.125 and 1.0 after two iterations.
    ! The first iteration is applied explicitly.
    num = 0.707106 * (0.707106**3 + 2.0 * asx)
    den = 2.0 * (0.707106**3) + asx

    do itt=1,2
      ! Halley's method iterates estimates as Root = Root * (Root**3 + 2.*asx) / (2.*Root**3 + asx).
      num_prev = num ; den_prev = den
      num = num_prev * (num_prev**3 + 2.0 * asx * (den_prev**3))
      den = den_prev * (2.0 * num_prev**3 + asx * (den_prev**3))
      ! Equivalent to:  root_asx = root_asx * (root_asx**3 + 2.*asx) / (2.*root_asx**3 + asx)
    enddo
    ! At this point the error in root_asx is better than 1 part in 3e14.
    root_asx = num / den

    ! One final iteration with Newton's method polishes up the root and gives a solution
    ! that is within the last bit of the true solution.
    root_asx = root_asx - (root_asx**3 - asx) / (3.0 * (root_asx**2))

    root = descale_cbrt(root_asx, e_x, s_x)
  endif
end function cuberoot


!> Rescale `a` to the range [0.125, 1) while preserving its fractional term.
pure subroutine rescale_exp(a, x, e_a, s_a)
  real, intent(in) :: a
    !< The value to be rescaled
  real, intent(out) :: x
    !< The rescaled value of `a`
  integer(kind=int64), intent(out) :: e_a
    !< The biased exponent of `a`
  integer(kind=int64), intent(out) :: s_a
    !< The sign bit of `a`

  ! Floating point model, if format is (sign, exp, frac)
  integer, parameter :: bias = maxexponent(1.) - 1
    !< The double precision exponent offset (assuming a balanced range)
  integer, parameter :: signbit = storage_size(1.) - 1
    !< Position of sign bit
  integer, parameter :: explen = 1 + ceiling(log(real(bias))/log(2.))
    !< Bit size of exponent
  integer, parameter :: expbit = signbit - explen
    !< Position of lowest exponent bit
  integer, parameter :: fraclen = expbit
    !< Length of fractional part

  integer(kind=int64) :: xb
    !< A floating point number, bit-packed as an integer
  integer(kind=int64) :: e_scaled
    !< The new rescaled exponent of `a` (i.e. the exponent of `x`)

  ! Pack bits of `a` into `xb` and extract its exponent and sign
  xb = transfer(a, 1_int64)
  s_a = ibits(xb, signbit, 1)
  e_a = ibits(xb, expbit, explen)

  ! Decompose the exponent as `e = modulo(e,3) + 3*(e/3)` and extract the
  ! rescaled exponent, now in {-3,-2,-1}
  e_scaled = modulo(e_a, 3) - 3 + bias

  ! Insert the new 11-bit exponent into `xb`, while also setting the sign bit
  ! to zero, ensuring that `xb` is always positive.
  call mvbits(e_scaled, 0, explen + 1, xb, fraclen)

  ! Transfer the final modified value to `x`
  x = transfer(xb, 1.)
end subroutine rescale_exp


!> Descale a real number to its original base, and apply the cube root to the
!! remaining exponent.
pure function descale_cbrt(x, e_a, s_a) result(r)
  real, intent(in) :: x
    !< Cube root of the rescaled value, which was rescaled to [0.125, 1.0)
  integer(kind=int64), intent(in) :: e_a
    !< Exponent of the original value to be cube rooted
  integer(kind=int64), intent(in) :: s_a
    !< Sign bit of the original value to be cube rooted
  real :: r
    !< Restored value with the cube root applied to its exponent

  ! Floating point model, if format is (sign, exp, frac)
  integer, parameter :: bias = maxexponent(1.) - 1
    !< The double precision exponent offset (assuming a balanced range)
  integer, parameter :: signbit = storage_size(1.) - 1
    !< Position of sign bit
  integer, parameter :: explen = 1 + ceiling(log(real(bias))/log(2.))
    !< Bit size of exponent
  integer, parameter :: expbit = signbit - explen
    !< Position of lowest exponent bit
  integer, parameter :: fraclen = expbit
    !< Length of fractional part

  integer(kind=int64) :: xb
    ! Bit-packed real number into integer form
  integer(kind=int64) :: e_r
    ! Exponent of the descaled value

  ! Extract the exponent of the rescaled value, in {-3, -2, -1}
  xb = transfer(x, 1_8)
  e_r = ibits(xb, expbit, explen)

  ! Apply the cube root to the old exponent (after removing its bias) and add
  ! to the rescaled exponent.  Correct the previous -3 with a +1.
  e_r = e_r + (e_a/3 - bias/3 + 1)

  ! Apply the corrected exponent and sign and convert back to real
  call mvbits(e_r, 0, explen, xb, expbit)
  call mvbits(s_a, 0, 1, xb, signbit)
  r = transfer(xb, 1.)
end function descale_cbrt



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
