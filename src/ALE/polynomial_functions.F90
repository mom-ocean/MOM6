module polynomial_functions
!==============================================================================
!
! This file is part of MOM.
!
! Date of creation: 2008.06.12
! L. White
!
! This module contains routines that handle polynomials.
!
!==============================================================================

implicit none ; private

public :: evaluation_polynomial, integration_polynomial

! -----------------------------------------------------------------------------
! This module contains the following routines
! -----------------------------------------------------------------------------
contains

! -----------------------------------------------------------------------------
! Pointwise evaluation of a polynomial
! -----------------------------------------------------------------------------
real function evaluation_polynomial( coefficients, nb_coefficients, x )
! -----------------------------------------------------------------------------
! The polynomial is defined by the coefficients contained in the
! array of the same name, as follows: C(1) + C(2)x + C(3)x^2 + C(4)x^3 + ...
! where C refers to the array 'coefficients'.
! The number of coefficients is given by nb_coefficients and x
! is the coordinate where the polynomial is to be evaluated.
!
! The function returns the value of the polynomial at x.
! -----------------------------------------------------------------------------

  ! Arguments
  real, dimension(:), intent(in)        :: coefficients
  integer, intent(in)                   :: nb_coefficients
  real, intent(in)                      :: x

  ! Local variables
  integer                               :: k
  real                                  :: f    ! value of polynomial at x

  f = 0.0
  do k = 1,nb_coefficients
    f = f + coefficients(k) * ( x**(k-1) )
  end do

  evaluation_polynomial = f

end function evaluation_polynomial

! -----------------------------------------------------------------------------
! Exact integration of polynomial of degree n
! -----------------------------------------------------------------------------
real function integration_polynomial( xi0, xi1, C, n )
! -----------------------------------------------------------------------------
! Exact integration of a polynomial of degree n over the interval [xi0,xi1].
! The array of coefficients (C) must be of size n+1, where n is the degree of
! the polynomial to integrate.
! -----------------------------------------------------------------------------

  ! Arguments
  real, intent(in)                  :: xi0, xi1
  real, dimension(:), intent(in)    :: C
  integer, intent(in)               :: n

  ! Local variables
  integer                           :: k
  real                              :: integral

  integral = 0.0

  do k = 1,(n+1)
    integral = integral + C(k) * (xi1**k - xi0**k) / real(k)
  end do
!
!One non-answer-changing way of unrolling the above is:
!  k=1
!  integral = integral + C(k) * (xi1**k - xi0**k) / real(k)
!  if (n>=1) then
!    k=2
!    integral = integral + C(k) * (xi1**k - xi0**k) / real(k)
!  endif
!  if (n>=2) then
!    k=3
!    integral = integral + C(k) * (xi1**k - xi0**k) / real(k)
!  endif
!  if (n>=3) then
!    k=4
!    integral = integral + C(k) * (xi1**k - xi0**k) / real(k)
!  endif
!  if (n>=4) then
!    k=5
!    integral = integral + C(k) * (xi1**k - xi0**k) / real(k)
!  endif
!
  integration_polynomial = integral

end function integration_polynomial

end module polynomial_functions
