dnl MOM6_FC_DO_CONCURRENT_LOCAL
dnl
dnl Determine if the compiler support do-concurrent locality specifiers.
dnl Currently only LOCAL is tested, but this should also include LOCAL_INIT,
dnl SHARED, and DEFAULT(NONE).
dnl
dnl Results are cached in the `mom6_cv_fc_do_concurrent_local` variable.
dnl
AC_DEFUN([MOM6_FC_DO_CONCURRENT_LOCAL], [
  AC_REQUIRE([AC_PROG_FC])
  AC_LANG_PUSH([Fortran])
  AC_CACHE_CHECK([whether $FC supports do concurrent locality specifiers],
    [mom6_cv_fc_do_concurrent_local], [
      mom6_cv_fc_do_concurrent_local=no
      AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM([], [
dnl ---
      do concurrent(i=1:2) local(a,b)
      end do
dnl ---
        ])
      ], [
        mom6_cv_fc_do_concurrent_local=yes
      ])
    ]
  )
  AS_IF([test "x$mom6_cv_fc_do_concurrent_local" = "xyes"], [
    AC_DEFINE([HAVE_FC_DO_CONCURRENT_LOCAL], [1],
      [Define if do concurrent supports locality specifiers.])
  ])
  AC_LANG_POP([Fortran])
])
