dnl Determine the flag required to force 64-bit reals.
dnl
dnl Many applications do not specify the kind of its real variables, even
dnl though the code may intrinsically require double-precision.  Most compilers
dnl will also default to using single-precision (32-bit) reals.
dnl
dnl This test determines the flag required to set reals without explcit kind to
dnl 64-bit double precision floats.  Ideally, we also desire to leave any
dnl `DOUBLE PRECISION` variable as 64-bit.  But this does not appear to always
dnl be possible, such as in NAG Fortran (see below).
dnl
dnl This does not test if the behavior of integers is changed; for example,
dnl Cray's Fortran wrapper's -default will double both.  This is addressed by
dnl avoiding any flags with affect integers, but this should still be used with
dnl some care.
dnl
dnl   GCC               -fdefault-real-8, -fdefault-double-8
dnl   [Common alias]    -r8
dnl   Intel Fortran     -real-kind 64
dnl   PGI Fortran       -Mr8
dnl   Cray Fortran      -s real64
dnl   NAG               -double
dnl
dnl NOTE:
dnl   - Many compilers accept -r8 for real and double precision sizes, but
dnl     several compiler-specific options are also provided.
dnl
dnl   - -r8 in NAG will attempt to also set double precision to 16 bytes if
dnl     available, which is generally undesired.
dnl
dnl     Additionally, the -double flag, which doubles *all* types, appears to
dnl     be the preferred flag here.
dnl
dnl     Neither flag describes what we actually want, but we include it here
dnl     as a last resort.
dnl
AC_DEFUN([AX_FC_REAL8],
[
  REAL8_FCFLAGS=
  AC_ARG_ENABLE([real8],
    [AS_HELP_STRING([--disable-real-8], [do not force 8-byte reals])])
  if test "$enable_real8" != no; then
    AC_CACHE_CHECK([for $FC option to force 8-byte reals],
      [ac_cv_prog_fc_real8],
      [AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM([], [
      real :: x(4)
      double precision :: y(4)
      integer, parameter :: &
        m = merge(1, 0, kind(x(1)) == selected_real_kind(15, 307)), &
        n = merge(1, 0, kind(y(1)) == selected_real_kind(15, 307))
      print *, x(::m)
      print *, y(::n)
        ])],
        [ac_cv_prog_fc_real8='none needed'],
        [ac_cv_prog_fc_real8='unsupported'
        for ac_option in "-fdefault-real-8 -fdefault-double-8" -r8 "-real-kind 64" -Mr8 "-s real64" -double; do
          ac_save_FCFLAGS=$FCFLAGS
          FCFLAGS="$FCFLAGS $ac_option"
          AC_LINK_IFELSE(
            [AC_LANG_PROGRAM([], [
      real :: x(4)
      double precision :: y(4)
      integer, parameter :: &
        m = merge(1, 0, kind(x(1)) == selected_real_kind(15, 307)), &
        n = merge(1, 0, kind(y(1)) == selected_real_kind(15, 307))
      print *, x(::m)
      print *, y(::n)
            ])],
            [ac_cv_prog_fc_real8=$ac_option]
          )
          FCFLAGS=$ac_save_FCFLAGS
          if test "$ac_cv_prog_fc_real8" != unsupported; then
            break
          fi
        done])
      ])
      case $ac_cv_prog_fc_real8 in #(
        "none needed" | unsupported)
      ;; #(
        *)
      REAL8_FCFLAGS=$ac_cv_prog_fc_real8 ;;
      esac
  fi
  AC_SUBST(REAL8_FCFLAGS)
])
