dnl AX_FC_CHECK_C_LIB(LIBRARY, FUNCTION,
dnl                   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
dnl                   [OTHER-LDFLAGS], [OTHER-LIBS])
dnl
dnl This macro checks if a C library can be referenced by a Fortran compiler.
dnl
dnl Results are cached in `ax_fc_cv_c_lib_LIBRARY_FUNCTION`.
dnl
dnl NOTE: Might be possible to rewrite this to use `AX_FC_CHECK_BIND_C`.
dnl
AC_DEFUN([AX_FC_CHECK_C_LIB], [
  AS_VAR_PUSHDEF([ax_fc_C_Lib], [ax_fc_cv_c_lib_$1_$2])
  m4_ifval([$5],
    [ax_fc_c_lib_msg_LDFLAGS=" with $5"],
    [ax_fc_c_lib_msg_LDFLAGS=""]
  )
  AC_CACHE_CHECK(
    [for $2 in -l$1$ax_fc_c_lib_msg_LDFLAGS], [ax_fc_cv_c_lib_$1_$2], [
      ax_fc_check_c_lib_save_LDFLAGS=$LDFLAGS
      LDFLAGS="$6 $LDFLAGS"
      ax_fc_check_c_lib_save_LIBS=$LIBS
      LIBS="-l$1 $7 $LIBS"
      AC_LINK_IFELSE(
        [AC_LANG_PROGRAM([],[dnl
dnl begin code block
        interface
        subroutine test() bind(c, name="$2")
        end subroutine test
        end interface
        call test])
dnl end code block
        ],
        [AS_VAR_SET([ax_fc_C_Lib], [yes])],
        [AS_VAR_SET([ax_fc_C_Lib], [no])]
      )
      LDFLAGS=$ax_fc_check_c_lib_save_LDFLAGS
      LIBS=$ax_fc_check_c_lib_save_LIBS
    ]
  )
  AS_VAR_IF([ax_fc_C_Lib], [yes],
    [m4_default([$3], [LIBS="-l$1 $LIBS"])],
    [$4]
  )
  AS_VAR_POPDEF([ax_fc_C_Lib])
])
