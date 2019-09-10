!> Polynomial functions
module polynomial_functions

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

public :: evaluation_polynomial, integration_polynomial, first_derivative_polynomial

contains

!> Pointwise evaluation of a polynomial at x
!!
!! The polynomial is defined by the coefficients contained in the
!! array of the same name, as follows: C(1) + C(2)x + C(3)x^2 + C(4)x^3 + ...
!! where C refers to the array 'coeff'.
!! The number of coefficients is given by ncoef and x
!! is the coordinate where the polynomial is to be evaluated.
real function evaluation_polynomial( coeff, ncoef, x )
  real, dimension(:), intent(in) :: coeff !< The coefficients of the polynomial
  integer,            intent(in) :: ncoef !< The number of polynomial coefficients
  real,               intent(in) :: x     !< The position at which to evaluate the polynomial
  ! Local variables
  integer :: k
  real    :: f    ! value of polynomial at x

  f = 0.0
  do k = 1,ncoef
    f = f + coeff(k) * ( x**(k-1) )
  enddo

  evaluation_polynomial = f

end function evaluation_polynomial

!> Calculates the first derivative of a polynomial evaluated at a point x
!!
!! The polynomial is defined by the coefficients contained in the
!! array of the same name, as follows: C(1) + C(2)x + C(3)x^2 + C(4)x^3 + ...
!! where C refers to the array 'coeff'.
!! The number of coefficients is given by ncoef and x
!! is the coordinate where the polynomial's derivative is to be evaluated.
real function first_derivative_polynomial( coeff, ncoef, x )
  real, dimension(:), intent(in) :: coeff !< The coefficients of the polynomial
  integer,            intent(in) :: ncoef !< The number of polynomial coefficients
  real, intent(in)               :: x     !< The position at which to evaluate the derivative
  ! Local variables
  integer                               :: k
  real                                  :: f    ! value of polynomial at x

  f = 0.0
  do k = 2,ncoef
    f = f + REAL(k-1)*coeff(k) * ( x**(k-2) )
  enddo

  first_derivative_polynomial = f

end function first_derivative_polynomial

!> Exact integration of polynomial of degree npoly
!!
!! The array of coefficients (Coeff) must be of size npoly+1.
real function integration_polynomial( xi0, xi1, Coeff, npoly )
  real,               intent(in) :: xi0   !< The lower bound of the integral
  real,               intent(in) :: xi1   !< The lower bound of the integral
  real, dimension(:), intent(in) :: Coeff !< The coefficients of the polynomial
  integer,            intent(in) :: npoly !< The degree of the polynomial
  ! Local variables
  integer                           :: k
  real                              :: integral

  integral = 0.0

  do k = 1,npoly+1
    integral = integral + Coeff(k) * (xi1**k - xi0**k) / real(k)
  enddo
!
!One non-answer-changing way of unrolling the above is:
!  k=1
!  integral = integral + Coeff(k) * (xi1**k - xi0**k) / real(k)
!  if (npoly>=1) then
!    k=2
!    integral = integral + Coeff(k) * (xi1**k - xi0**k) / real(k)
!  endif
!  if (npoly>=2) then
!    k=3
!    integral = integral + Coeff(k) * (xi1**k - xi0**k) / real(k)
!  endif
!  if (npoly>=3) then
!    k=4
!    integral = integral + Coeff(k) * (xi1**k - xi0**k) / real(k)
!  endif
!  if (npoly>=4) then
!    k=5
!    integral = integral + Coeff(k) * (xi1**k - xi0**k) / real(k)
!  endif
!
  integration_polynomial = integral

end function integration_polynomial

!> \namespace polynomial_functions
!!
!! Date of creation: 2008.06.12
!! L. White
!!
!! This module contains routines that handle polynomials.

end module polynomial_functions
