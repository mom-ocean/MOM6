dnl AX_FC_CHECK_MODULE(MODULE,
dnl                    [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
dnl                    [OTHER-FCFLAGS])
dnl
dnl This macro checks if a Fortran module is available to the compiler.
dnl
dnl The fourth argument (optional) allows for specification of supplemental
dnl FCFLAGS arguments.  This would primarily be used to test additional
dnl paths (typically using -I) for the module file.
dnl
dnl Results are cached in the ax_fc_cv_mod_MODULE variable.
dnl
AC_DEFUN([AX_FC_CHECK_MODULE],
[
  AS_VAR_PUSHDEF([ax_fc_Module], [ax_fc_cv_mod_$1])
  AC_CACHE_CHECK([if $FC can use module $1$ax_fc_mod_msg_FCFLAGS], [ax_fc_cv_mod_$1],[
    ax_fc_chk_mod_save_FCFLAGS=$FCFLAGS
    FCFLAGS="$4 $FCFLAGS"
    AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM([],[use $1])],
        [AS_VAR_SET([ax_fc_Module], [yes])],
        [AS_VAR_SET([ax_fc_Module], [no])]
    )
    FCFLAGS=$ax_fc_chk_mod_save_FCFLAGS
  ])
  AS_VAR_IF([ax_fc_Module], [yes], [$2], [$3])
  AS_VAR_POPDEF([ax_fc_Module])
])
