dnl AX_FC_CHECK_C_LIB(FUNCTION,
dnl                   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
dnl                   [OTHER-LDFLAGS], [OTHER-LIBS])
dnl
dnl This macro checks if a C binding is available to the compiler.
dnl
dnl Equivalently, it checks if the Fortran compiler can see a C function.
dnl
dnl Results are cached in `ax_fc_cv_bind_c_FUNCTION`.
dnl
AC_DEFUN([AX_FC_CHECK_BIND_C], [
  AS_VAR_PUSHDEF([ax_fc_Bind_C], [ax_fc_cv_bind_c_$1])
  m4_ifval([$4],
    [ax_fc_bind_c_msg_LDFLAGS=" with $4"],
    [ax_fc_bind_c_msg_LDFLAGS=""]
  )
  AC_CACHE_CHECK(
    [if $FC can bind $1$ax_fc_bind_c_msg_LDFLAGS], [ax_fc_cv_bind_c_$1], [
      ax_fc_check_bind_c_save_LDFLAGS=$LDFLAGS
      LDFLAGS="$4 $LDFLAGS"
      ax_fc_check_bind_c_save_LIBS=$LIBS
      LIBS="$5 $LIBS"
      AC_LINK_IFELSE(
        [AC_LANG_PROGRAM([],[dnl
dnl begin code block
        interface
        subroutine test() bind(c, name="$1")
        end subroutine test
        end interface
        call test])
dnl end code block
        ],
        [AS_VAR_SET([ax_fc_Bind_C], [yes])],
        [AS_VAR_SET([ax_fc_Bind_C], [no])]
      )
      LDFLAGS=$ax_fc_check_bind_c_save_LDFLAGS
      LIBS=$ax_fc_check_bind_c_save_LIBS
    ]
  )
  AS_VAR_IF([ax_fc_Bind_C], [yes], [$2], [$3])
  AS_VAR_POPDEF([ax_fc_Bind_C])
])
