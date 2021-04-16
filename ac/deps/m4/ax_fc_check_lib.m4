dnl AX_FC_CHECK_LIB(LIBRARY, FUNCTION,
dnl                 [MODULE], [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
dnl                 [OTHER-LDFLAGS], [OTHER-LIBS])
dnl
dnl This macro checks if a Fortran library containing a designated function
dnl is available to the compiler.  For the most part, this macro should behave
dnl like the Autoconf AC_CHECK_LIB macro.
dnl
dnl This macro differs somewhat from AC_CHECK_LIB, since it includes two
dnl additional features:
dnl
dnl 1. The third argument (optional) allows us to specify a Fortran module,
dnl    which may be required to access the library's functions.
dnl
dnl 2. The sixth argument (optional) allows specification of supplemental
dnl    LDFLAGS arguments.  This can be used, for example, to test for the
dnl    library with different -L flags, or perhaps other ld configurations.
dnl
dnl Results are cached in the ax_fc_cv_lib_LIBRARY_FUNCTION variable.
dnl
AC_DEFUN([AX_FC_CHECK_LIB],[dnl
  AS_VAR_PUSHDEF([ax_fc_Lib], [ax_fc_cv_lib_$1_$2])
  m4_ifval([$6],
    [ax_fc_lib_msg_LDFLAGS=" with $6"],
    [ax_fc_lib_msg_LDFLAGS=""]
  )
  AC_CACHE_CHECK([for $2 in -l$1$ax_fc_lib_msg_LDFLAGS], [ax_fc_cv_lib_$1_$2],[
    ax_fc_check_lib_save_LDFLAGS=$LDFLAGS
    LDFLAGS="$6 $LDFLAGS"
    ax_fc_check_lib_save_LIBS=$LIBS
    LIBS="-l$1 $7 $LIBS"
    AS_IF([test -n $3],
      [ax_fc_use_mod="use $3"],
      [ax_fc_use_mod=""])
    AC_LINK_IFELSE([
        AC_LANG_PROGRAM([], [dnl
        $ax_fc_use_mod
        call $2]dnl
        )
      ],
      [AS_VAR_SET([ax_fc_Lib], [yes])],
      [AS_VAR_SET([ax_fc_Lib], [no])]
    )
    LIBS=$ax_fc_check_lib_save_LIBS
    LDFLAGS=$ax_fc_check_lib_save_LDFLAGS
  ])
  AS_VAR_IF([ax_fc_Lib], [yes],
    [m4_default([$4], [LIBS="-l$1 $LIBS"])],
    [$5]
  )
  AS_VAR_POPDEF([ax_fc_Lib])
])
