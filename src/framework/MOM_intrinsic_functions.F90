module MOM_intrinsic_functions

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*   This file is a part of MOM.  See MOM.F90 for licensing.           *
!*                                                                     *
!*   This module holds intrinsic functions which are used by MOM but   *
!* are not supported  by some compilers.                               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

  implicit none
  private

  public :: invcosh

  contains

  function invcosh(x)
    real, intent(in) :: x
    real :: invcosh

#ifdef __INTEL_COMPILER
    invcosh=acosh(x)
#else
    invcosh=log(x+sqrt(x*x-1))
#endif

  end function invcosh

end module MOM_intrinsic_functions
