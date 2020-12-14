dnl AX_FC_CRAY_POINTER([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE])
dnl
dnl This macro tests if any flags are required to enable Cray pointers.
dnl
dnl Cray pointers provided a means for more direct access to memory.  Since
dnl such references can potentially violate certain requirements of the
dnl language standard, they are typically considered a vendor extension.
dnl
dnl Most compilers provide these in some form.  A partial list of supported
dnl flags are shown below, but additional feedback is required for other
dnl compilers.
dnl
dnl The known flags are:
dnl     GCC               -fcray-pointer
dnl     Intel Fortran     none
dnl     PGI Fortran       -Mcray=pointer
dnl     Cray Fortran      none
dnl
AC_DEFUN([AX_FC_CRAY_POINTER], [
  AC_LANG_ASSERT([Fortran])
  AC_MSG_CHECKING([for $FC option to support Cray pointers])
  AC_CACHE_VAL([ac_cv_fc_cray_ptr], [
    ac_cv_fc_cray_ptr='unknown'
    ac_save_FCFLAGS=$FCFLAGS
    for ac_option in none -fcray-pointer -Mcray=pointer; do
      test "$ac_option" != none && FCFLAGS="$ac_save_FCFLAGS $ac_option"
      AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM([], [
      integer aptr(2)
      pointer (iptr, aptr)
        ])],
        [ac_cv_fc_cray_ptr=$ac_option],
      )
      FCFLAGS=$ac_save_FCFLAGS
      AS_IF([test "$ac_cv_fc_cray_ptr" != unknown], [break])
    done
  ])
  AS_CASE([ac_cv_fc_cray_ptr],
    [none], [AC_MSG_RESULT([none_needed])],
    [unknown], [AC_MSG_RESULT([unsupported])],
    [AC_MSG_RESULT([$ac_cv_fc_cray_ptr])]
  )
  AS_IF([test "$ac_cv_fc_cray_ptr" != unknown], [
    m4_default([$1], [
      AS_IF([test "$ac_cv_fc_cray_ptr" != none],
        [FCFLAGS="$FCFLAGS $ac_cv_fc_cray_ptr"]
      )
    ])],
    [m4_default([$2], [AC_MSG_ERROR(["$FC does not support Cray pointers"])])]
  )
])
