dnl MOM6_FC_CHECK_LIB(LIBRARY, PROCEDURE,
dnl                   [MODULE], [ARGS], [FUNC-RESULT], [DECLS],
dnl                   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
dnl                   [OTHER-LDFLAGS], [OTHER-LIBS])
dnl
dnl This macro checks if a Fortran library containing a designated function
dnl is available to the compiler.  For the most part, this macro should behave
dnl like the Autoconf AC_CHECK_LIB macro.
dnl
dnl This macro differs from AC_CHECK_LIB, since it includes several additional
dnl arguments.  Although the next four arguments are optional, they are
dnl required for many function tests.
dnl
dnl - MODULE specifies the Fortran module containing the procedure.
dnl
dnl - ARGS is used to specify any arguments of the procedure.
dnl
dnl - FUNC-RESULT, if set, identifies the procedure as a function rather than
dnl   a subroutine, and specifies the function test result.
dnl
dnl - DECLS is used as a code block to explicitly declare variables, when
dnl   implicit typing is not sufficient.
dnl
dnl The following argument has also been added.
dnl
dnl - OTHER-LDFLAGS allows specification of supplemental LDFLAGS arguments.
dnl   This can be used, for example, to test for the library with different
dnl   -L flags, or perhaps other ld configurations.
dnl
dnl Results are cached in the mom6_fc_cv_lib_LIBRARY_PROCEDURE variable.
dnl
AC_DEFUN([MOM6_FC_CHECK_LIB],[
  AS_VAR_PUSHDEF([mom6_fc_Lib], [mom6_fc_cv_lib_$1_$2])
  m4_ifval([$9],
    [mom6_fc_lib_msg_LDFLAGS=" with $9"],
    [mom6_fc_lib_msg_LDFLAGS=""]
  )
  AC_CACHE_CHECK(
    [for $2 in -l$1$mom6_fc_lib_msg_LDFLAGS],
    [mom6_fc_Lib],[
      mom6_fc_check_lib_save_LDFLAGS=$LDFLAGS
      LDFLAGS="$9 $LDFLAGS"
      mom6_fc_check_lib_save_LIBS=$LIBS
      LIBS="-l$1 $10 $LIBS"
      AS_IF([test -n "$3"],
        [mom6_fc_use_mod="use $3"],
        [mom6_fc_use_mod=""]
      )
      AS_IF([test -n "$5"],
        [mom6_fc_proc="$5 = $2"],
        [mom6_fc_proc="call $2"]
      )
      AS_IF([test -n "$4"],
        [mom6_fc_proc="${mom6_fc_proc}($4)"]
      )
      AS_IF([test -n "$6"],
        [mom6_fc_decls="$6"],
        [mom6_fc_decls=""]
      )
      AC_LANG_PUSH([Fortran])
      AC_LINK_IFELSE([dnl
dnl Begin 7-column code block
AC_LANG_PROGRAM([], [dnl
        $mom6_fc_use_mod
        $mom6_fc_decls
        $mom6_fc_proc])dnl
dnl End code block
        ],
        [AS_VAR_SET([mom6_fc_Lib], [yes])],
        [AS_VAR_SET([mom6_fc_Lib], [no])]
      )
      AC_LANG_POP([Fortran])
      LIBS=$mom6_fc_check_lib_save_LIBS
      LDFLAGS=$mom6_fc_check_lib_save_LDFLAGS
    ]
  )
  AS_VAR_IF([mom6_fc_Lib], [yes],
    [m4_default([$7], [LIBS="-l$1 $LIBS"])],
    [$8]
  )
  AS_VAR_POPDEF([mom6_fc_Lib])
])
