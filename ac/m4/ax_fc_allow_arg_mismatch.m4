dnl Test if mismatched function arguments are permitted.
dnl
dnl This macro tests if a flag is required to enable mismatched functions in
dnl a single translation unit (aka file).
dnl
dnl If a compiler encounters two undefined programs with different input
dnl argument types, then it may regard this as a mismatch which requires action
dnl from the user.  A common example is a procedure which may be called with
dnl a variable of either an integer or a real type.
dnl
dnl This can happen, for example, if one is relying on an interface to resolve
dnl such differences, but one is also relying on a legacy header interface via
dnl `#include` rather than an explicit module which includes the complete
dnl interface specification.
dnl
dnl No modern project is expected to see these issues, but this is helpful for
dnl older projects which used legacy headers.
dnl
dnl Flags:
dnl     GNU: -fallow-argument-mismatch
dnl
AC_DEFUN([AX_FC_ALLOW_ARG_MISMATCH],
  [ALLOW_ARG_MISMATCH_FCFLAGS=
  AC_CACHE_CHECK(
    [for $FC option to support mismatched procedure arguments],
    [ac_cv_prog_fc_arg_mismatch],
    [AC_COMPILE_IFELSE(
      [AC_LANG_PROGRAM([], [
      call f(1)
      call f(1.0)
      ])],
      [ac_cv_prog_fc_arg_mismatch='none needed'],
      [ac_cv_prog_fc_arg_mismatch='unsupported'
      for ac_option in -fallow-argument-mismatch; do
        ac_save_FCFLAGS=$FCFLAGS
        FCFLAGS="$FCFLAGS $ac_option"
        AC_COMPILE_IFELSE(
          [AC_LANG_PROGRAM([], [
      call f(1)
      call f(1.0)
          ])],
          [ac_cv_prog_fc_arg_mismatch=$ac_option]
        )
        FCFLAGS=$ac_save_FCFLAGS
        if test "$ac_cv_prog_fc_arg_mismatch" != unsupported; then
          break
        fi
      done])
    ]
  )
  case $ac_cv_prog_fc_arg_mismatch in #(
    "none needed" | unsupported)
  ;; #(
    *)
  ALLOW_ARG_MISMATCH_FCFLAGS=$ac_cv_prog_fc_arg_mismatch ;;
  esac
  AC_SUBST(ALLOW_ARG_MISMATCH_FCFLAGS)
])
