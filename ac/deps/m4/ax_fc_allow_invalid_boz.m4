dnl Test if BOZ literal assignment is supported.
dnl
dnl This macro tests if a flag is required to enable BOZ literal assignments
dnl for variables.
dnl
dnl BOZ literals (e.g. Z'FFFF') are typeless, and formally cannot be assigned
dnl to typed variables.  Nonetheless, few compilers forbid such operations,
dnl despite the potential pitfalls around interpreting such values.
dnl
dnl As of version 10.1, gfortran now forbids such assignments and requires a
dnl flag to convert the raised errors into warnings.
dnl
dnl While the best solution is to replace such assignments with proper
dnl conversion functions, this test is useful to accommodate older projects.
dnl
dnl Flags:
dnl     GNU: -fallow-invalid-boz
AC_DEFUN([AX_FC_ALLOW_INVALID_BOZ],
  [ALLOW_INVALID_BOZ_FCFLAGS=
  AC_CACHE_CHECK(
    [for $FC option to support invalid BOZ assignment],
    [ac_cv_prog_fc_invalid_boz],
    [AC_COMPILE_IFELSE(
      [AC_LANG_PROGRAM([], [
      integer n
      n = z'ff'
      ])],
      [ac_cv_prog_fc_invalid_boz='none needed'],
      [ac_cv_prog_fc_invalid_boz='unsupported'
      for ac_option in -fallow-invalid-boz; do
        ac_save_FCFLAGS=$FCFLAGS
        FCFLAGS="$FCFLAGS $ac_option"
        AC_LINK_IFELSE(
          [AC_LANG_PROGRAM([], [
      integer n
      n = z'ff'
          ])],
          [ac_cv_prog_fc_invalid_boz=$ac_option]
        )
        FCFLAGS=$ac_save_FCFLAGS
        if test "$ac_cv_prog_fc_invalid_boz" != unsupported; then
          break
        fi
      done])
    ]
  )
  case $ac_cv_prog_fc_invalid_boz in #(
    "none needed" | unsupported)
  ;; #(
    *)
  ALLOW_INVALID_BOZ_FCFLAGS=$ac_cv_prog_fc_invalid_boz ;;
  esac
  AC_SUBST(ALLOW_INVALID_BOZ_FCFLAGS)]
)
