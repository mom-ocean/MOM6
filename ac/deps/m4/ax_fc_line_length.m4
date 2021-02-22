# AX_FC_LINE_LENGTH([LENGTH], [ACTION-IF-SUCCESS],
#		    [ACTION-IF-FAILURE = FAILURE])
# ------------------------------------------------
# This is a backport of the AC_FC_LINE_LENGTH macro in Autoconf 2.67 and newer.
# Comments below are from the Autoconf 2.69 implementation.
#
# Look for a compiler flag to make the Fortran (FC) compiler accept long lines
# in the current (free- or fixed-format) source code, and adds it to FCFLAGS.
# The optional LENGTH may be 80, 132 (default), or `unlimited' for longer
# lines.  Note that line lengths above 250 columns are not portable, and some
# compilers (hello ifort) do not accept more than 132 columns at least for
# fixed format.  Call ACTION-IF-SUCCESS (defaults to nothing) if successful
# (i.e. can compile code using new extension) and ACTION-IF-FAILURE (defaults
# to failing with an error message) if not.  (Defined via DEFUN_ONCE to
# prevent flag from being added to FCFLAGS multiple times.)
# You should call AC_FC_FREEFORM or AC_FC_FIXEDFORM to set the desired format
# prior to using this macro.
#
# The known flags are:
# -f{free,fixed}-line-length-N with N 72, 80, 132, or 0 or none for none.
# -ffree-line-length-none: GNU gfortran
# -ffree-line-length-huge: g95 (also -ffixed-line-length-N as above)
#       -qfixed=132 80 72: IBM compiler (xlf)
#                -Mextend: Cray
#            -132 -80 -72: Intel compiler (ifort)
#                          Needs to come before -extend_source because ifort
#                          accepts that as well with an optional parameter and
#                          doesn't fail but only warns about unknown arguments.
#          -extend_source: SGI compiler
#  -W, -WNN (132, 80, 72): Absoft Fortran
#     +es, +extend_source: HP Fortran (254 in either form, default is 72 fixed,
#			   132 free)
#            -w, (-)-wide: Lahey/Fujitsu Fortran (255 cols in fixed form)
#                      -e: Sun Fortran compiler (132 characters)
#                    -132: NAGWare
#         -72, -f, -Wf,-f: f2c (a weak form of "free-form" and long lines).
#                  /XLine: Open Watcom

AC_DEFUN_ONCE([AX_FC_LINE_LENGTH], [
  AC_LANG_ASSERT([Fortran])
  m4_case(m4_default([$1], [132]),
    [unlimited], [
      ac_fc_line_len_string=unlimited
      ac_fc_line_len=0
      ac_fc_line_length_test='
      subroutine longer_than_132(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,'\
'arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19)'
    ],
    [132], [
      ac_fc_line_len=132
      ac_fc_line_length_test='
      subroutine longer_than_80(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,'\
'arg10)'
    ],
    [80], [
      ac_fc_line_len=80
      ac_fc_line_length_test='
      subroutine longer_than_72(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)'
    ],
    [m4_warning([Invalid length argument `$1'])]
  )
  : ${ac_fc_line_len_string=$ac_fc_line_len}
  AC_MSG_CHECKING([for Fortran flag needed to accept $ac_fc_line_len_string column source lines])
  AC_CACHE_VAL([ac_cv_fc_line_length], [
    ac_cv_fc_line_length=unknown
    ac_save_FCFLAGS=$FCFLAGS
    for ac_flag in none \
      -ffree-line-length-none \
      -ffixed-line-length-none \
      -ffree-line-length-huge \
      -ffree-line-length-$ac_fc_line_len \
      -ffixed-line-length-$ac_fc_line_len \
      -qfixed=$ac_fc_line_len \
      -Mextend \
      -$ac_fc_line_len \
      -extend_source \
      -W$ac_fc_line_len \
      -W +extend_source +es -wide --wide -w -e -f -Wf,-f -xline
    do
      test "$ac_flag" != none && FCFLAGS="$ac_save_FCFLAGS $ac_flag"
      AC_COMPILE_IFELSE([$ac_fc_line_length_test
      end subroutine
        ], [ac_cv_fc_line_length=$ac_flag]
      )
      FCFLAGS=$ac_save_FCFLAGS
      dnl TODO: Remove conftest.{err,$ac_objext,$ac_ext} ??
      AS_IF([test "$ac_cv_fc_line_length" != unknown], [break])
    done
  ])
  AC_MSG_RESULT([$ac_cv_fc_line_length])
  AS_IF([test "$ac_cv_fc_line_length" != unknown], [
    m4_default([$2], [
      AS_IF([test "$ac_cv_fc_line_length" != none], [
        FCFLAGS="$FCFLAGS $ac_cv_fc_line_length"
      ])
    ])], [
    m4_default([$3], [
      AC_MSG_ERROR([Fortran does not accept long source lines], 77)
    ])
  ])
])
