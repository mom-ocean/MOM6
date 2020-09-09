dnl Test if Cray pointers are supported.
dnl
dnl This macro tests if any flags are required to enable Cray pointers.
dnl
dnl Cray pointers provided a means for more direct access to memory.  Since
dnl such references can potentially violate certain requirements of the
dnl language standard, they are typically considered a vendor extension.
dnl
dnl Most compilers provide these in some form.  As far as I can tell, only GNU
dnl explicitly requires a flag.  Known tests are shown below, but additional
dnl feedback is required to fill this out.
dnl
dnl Flags
dnl     GCC               -fcray-pointer
dnl     Intel Fortran     none
dnl     PGI Fortran       none
dnl     Cray Fortran      none
dnl
AC_DEFUN([AX_FC_CRAY_POINTER],
  [CRAY_POINTER_FCFLAGS=
  AC_CACHE_CHECK([for $FC option to support Cray pointers],
    [ac_cv_prog_fc_cray_ptr],
    [AC_COMPILE_IFELSE(
      [AC_LANG_PROGRAM([], [
      integer aptr(2)
      pointer (iptr, aptr)
      ])],
      [ac_cv_prog_fc_cray_ptr='none needed'],
      [ac_cv_prog_fc_cray_ptr='unsupported'
      for ac_option in -fcray-pointer; do
        ac_save_FCFLAGS=$FCFLAGS
        FCFLAGS="$FCFLAGS $ac_option"
        AC_LINK_IFELSE(
          [AC_LANG_PROGRAM([], [
          integer aptr(2)
          pointer (iptr, aptr)
          ])],
          [ac_cv_prog_fc_cray_ptr=$ac_option]
        )
        FCFLAGS=$ac_save_FCFLAGS
        if test "$ac_cv_prog_fc_cray_ptr" != unsupported; then
          break
        fi
      done])
  ])
  case $ac_cv_prog_fc_cray_ptr in #(
    "none needed" | unsupported)
  ;; #(
    *)
  CRAY_POINTER_FCFLAGS=$ac_cv_prog_fc_cray_ptr ;;
  esac
  AC_SUBST(CRAY_POINTER_FCFLAGS)
])
