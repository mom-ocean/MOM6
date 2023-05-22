# Autoconf configuration

# NOTE:
# - We currently do not use a MOM6 version tag, but this would be one option in
#   the future:
#     [m4_esyscmd_s([git describe])]
# - Another option is `git rev-parse HEAD` for the full hash.
# - We would probably run this inside of a script to avoid the explicit
#   dependency on git.

AC_PREREQ([2.63])

AC_INIT(
  [MOM6],
  [ ],
  [https://github.com/NOAA-GFDL/MOM6/issues],
  [],
  [https://github.com/NOAA-GFDL/MOM6])

#---
# NOTE: For the autoconf-adverse, the configuration files and autoreconf output
#   are kept in the `ac` directory.
#
# This breaks the convention where configure.ac resides in the top directory.
#
# As a result, $srcdir initially points to the `ac` directory, rather than the
# top directory of the codebase.
#
# In order to balance this, we up-path (../) srcdir and point AC_CONFIG_SRCDIR
# to srcdir and point AC_CONFIG_SRCDIR to the parent directory.
#
# Someday we may revert this and work from the top-level directory.  But for
# now we will isolate autoconf to a subdirectory.
#---

# Validate srdcir and configure input
AC_CONFIG_SRCDIR([../src/core/MOM.F90])
AC_CONFIG_MACRO_DIR([m4])
srcdir=$srcdir/..


# Configure the memory layout header

AC_ARG_VAR([MOM_MEMORY],
  [Path to MOM_memory.h header, describing the field memory layout: dynamic
   symmetric (default), dynamic asymmetric, or static.]
)

AS_VAR_IF([MOM_MEMORY], [],
  [MOM_MEMORY=${srcdir}/config_src/memory/dynamic_symmetric/MOM_memory.h]
)

# Confirm that MOM_MEMORY is named 'MOM_memory.h'
AS_IF([test $(basename "${MOM_MEMORY}") == "MOM_memory.h"], [],
 [AC_MSG_ERROR([MOM_MEMORY header ${MOM_MEMORY} must be named 'MOM_memory.h'])]
)

# Confirm that the file exists
AC_CHECK_FILE(["$MOM_MEMORY"], [],
  [AC_MSG_ERROR([MOM_MEMORY header ${MOM_MEMORY} not found.])]
)

MOM_MEMORY_DIR=$(AS_DIRNAME(["${MOM_MEMORY}"]))
AC_SUBST([MOM_MEMORY_DIR])


# Default to solo_driver
DRIVER_DIR=${srcdir}/config_src/drivers/solo_driver
AC_ARG_WITH([driver],
  AS_HELP_STRING(
    [--with-driver=coupled_driver|solo_driver|unit_tests],
    [Select directory for driver source code]
  )
)
AS_IF([test "x$with_driver" != "x"],
  [DRIVER_DIR=${srcdir}/config_src/drivers/${with_driver}])

# TODO: Rather than point to a pre-configured header file, autoconf could be
# used to configure a header based on a template.
#AC_CONFIG_HEADERS(["$MEM_LAYOUT/MOM_memory.h"])

# Select the model framework (default: FMS1)
# NOTE: We can phase this out after the FMS1 I/O has been removed from FMS and
#   replace with a detection test.  For now, it is a user-defined switch.
MODEL_FRAMEWORK=${srcdir}/config_src/infra/FMS1
AC_ARG_WITH([framework],
  AS_HELP_STRING([--with-framework=fms1|fms2], [Select the model framework]))
AS_CASE(["$with_framework"],
  [fms1], [MODEL_FRAMEWORK=${srcdir}/config_src/infra/FMS1],
  [fms2], [MODEL_FRAMEWORK=${srcdir}/config_src/infra/FMS2],
  [MODEL_FRAMEWORK=${srcdir}/config_src/infra/FMS1]
)


# Explicitly assume free-form Fortran
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)


# Determine MPI compiler wrappers
# NOTE:
# - AX_MPI invokes AC_PROG_FC, often with gfortran, even if the MPI launcher
#   does not use gfortran.
# - This can cause standard AC_PROG_FC tests to fail if FCFLAGS is configured
#   with flags from another compiler.
# - I do not yet know how to resolve this possible issue.
AX_MPI([], [
  AC_MSG_ERROR([Could not find MPI launcher.])
])


# Explicitly replace FC and LD with MPI wrappers
# NOTE: This is yet another attempt to manage the potential mismatches between
#   FC and MPIFC.  Without this step, the tests below would not use MPIFC.
AC_SUBST(FC, $MPIFC)
AC_SUBST(LD, $MPIFC)

# Confirm that FC can see the Fortran 90 MPI module.
AX_FC_CHECK_MODULE([mpi],
  [], [AC_MSG_ERROR([Could not find MPI Fortran module.])])


# netCDF configuration

# Search for the Fortran netCDF module.
AX_FC_CHECK_MODULE([netcdf], [], [
  AS_UNSET([ax_fc_cv_mod_netcdf])
  AC_PATH_PROG([NF_CONFIG], [nf-config])
  AS_IF([test -n "$NF_CONFIG"], [
    AC_SUBST([FCFLAGS], ["$FCFLAGS -I$($NF_CONFIG --includedir)"])
    ], [AC_MSG_ERROR([Could not find nf-config.])]
  )
  AX_FC_CHECK_MODULE([netcdf], [], [
    AC_MSG_ERROR([Could not find netcdf module.])
  ])
])

# Confirm that the Fortran compiler can link the netCDF C library
AX_FC_CHECK_C_LIB([netcdf], [nc_create], [], [
  AS_UNSET([ax_fc_cv_c_lib_netcdf_nc_create])
  AC_PATH_PROG([NC_CONFIG], [nc-config])
  AS_IF([test -n "$NC_CONFIG"], [
    AC_SUBST([LDFLAGS], ["$LDFLAGS -L$($NC_CONFIG --libdir)"])
  ], [
    AC_MSG_ERROR([Could not find nc-config.])
  ])
  AX_FC_CHECK_C_LIB([netcdf], [nc_create], [], [
    AC_MSG_ERROR([Could not find netCDF C library.])
  ])
])

# Confirm that the Fortran compiler can link to the netCDF Fortran library.
# NOTE:
#   - We test nf_create, rather than nf90_create, since AX_FC_CHECK_LIB can
#     not currently probe the Fortran 90 interfaces.
#   - nf-config does not have --libdir, so we parse the --flibs output.
AX_FC_CHECK_LIB([netcdff], [nf_create], [], [], [
  AS_UNSET([ax_fc_cv_lib_netcdff_nf_create])
  AC_PATH_PROG([NF_CONFIG], [nf-config])
  AS_IF([test -n "$NF_CONFIG"], [
    AC_SUBST([LDFLAGS],
      ["$LDFLAGS $($NF_CONFIG --flibs | xargs -n1 | grep "^-L" | sort -u | xargs)"]
    )
  ], [
    AC_MSG_ERROR([Could not find nf-config.])
  ])
  AX_FC_CHECK_LIB([netcdff], [nf_create], [], [], [
    AC_MSG_ERROR([Could not find netCDF Fortran library.])
  ])
])


# Force 8-byte reals
AX_FC_REAL8
AS_IF(
  [test "$enable_real8" != no],
  [FCFLAGS="$FCFLAGS $REAL8_FCFLAGS"])


# OpenMP configuration

# NOTE: AC_OPENMP fails on `Fortran` for Autoconf <2.69 due to a m4 bug.
#   For older versions, we test against CC and use the result for FC.
m4_version_prereq([2.69], [AC_OPENMP], [
  AC_LANG_PUSH([C])
  AC_OPENMP
  AC_LANG_POP([C])
  OPENMP_FCFLAGS="$OPENMP_CFLAGS"
])

# NOTE: Only apply OpenMP flags if explicitly enabled.
AS_IF(
  [test "$enable_openmp" = yes], [
  FCFLAGS="$FCFLAGS $OPENMP_FCFLAGS"
  LDFLAGS="$LDFLAGS $OPENMP_FCFLAGS"
])


# FMS support

# Test for fms_mod to verify FMS module access
AX_FC_CHECK_MODULE([fms_mod], [], [
  AS_UNSET([ax_fc_cv_mod_fms_mod])
  AX_FC_CHECK_MODULE([fms_mod],
    [AC_SUBST([FCFLAGS], ["-I${srcdir}/ac/deps/include $FCFLAGS"])],
    [AC_MSG_ERROR([Could not find fms_mod Fortran module.])],
    [-I${srcdir}/ac/deps/include])
])

# Test for fms_init to verify FMS library linking
AX_FC_CHECK_LIB([FMS], [fms_init], [fms_mod],
  [], [
    AS_UNSET([ax_fc_cv_lib_FMS_fms_init])
    AX_FC_CHECK_LIB([FMS], [fms_init], [fms_mod], [
      AC_SUBST([LDFLAGS], ["-L${srcdir}/ac/deps/lib $LDFLAGS"])
      AC_SUBST([LIBS], ["-lFMS $LIBS"])
    ],
    [AC_MSG_ERROR([Could not find FMS library.])],
    [-L${srcdir}/ac/deps/lib])
  ]
)


# Verify that FMS is at least 2019.01.02
# NOTE: 2019.01.02 introduced two changes:
#   - diag_axis_init supports an optional domain_position argument
#   - position values NORTH, EAST, CENTER were added to diag_axis_mod
# For our versioning test, we check the second feature.
AC_MSG_CHECKING([if diag_axis_mod supports domain positions])
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([], [use diag_axis_mod, only: NORTH, EAST, CENTER])],
  [AC_MSG_RESULT([yes])],
  [
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([diag_axis_mod in MOM6 requires FMS 2019.01.02 or newer.])
  ]
)

# Python interpreter test

AC_ARG_VAR([PYTHON], [Python interpreter command])

AS_VAR_SET_IF([PYTHON], [
  AC_PATH_PROGS([PYTHON], ["$PYTHON"], [none])
], [
  AC_PATH_PROGS([PYTHON], [python python3 python2], [none])
])
AS_VAR_IF([PYTHON], [none], [
  AC_MSG_ERROR([Python interpreter not found.])
])


# Makedep test
AC_PATH_PROG([MAKEDEP], [makedep], [${srcdir}/ac/makedep])
AC_SUBST([MAKEDEP])


# Generate source list and configure dependency command
AC_SUBST([SRC_DIRS], ["\\
  ${srcdir}/src \\
  ${MODEL_FRAMEWORK} \\
  ${srcdir}/config_src/external \\
  ${DRIVER_DIR} \\
  ${MOM_MEMORY_DIR}"]
)
AC_CONFIG_COMMANDS(Makefile.dep, [make depend])


# POSIX verification tests

# Symbols in <setjmp.h> may be defined as macros, making them inaccessible by
# Fortran C bindings.  `sigsetjmp` is known to have an internal symbol in
# glibc, so we check for this possibility.  For the others, we only check for
# existence.

# If the need arises, we may want to define these under a standalone macro.

# Validate the setjmp symbol
AX_FC_CHECK_BIND_C([setjmp],
  [SETJMP="setjmp"], [SETJMP="setjmp_missing"]
)
AC_DEFINE_UNQUOTED([SETJMP_NAME], ["${SETJMP}"])

# Validate the longjmp symbol
AX_FC_CHECK_BIND_C([longjmp],
  [LONGJMP="longjmp"], [LONGJMP="longjmp_missing"]
)
AC_DEFINE_UNQUOTED([LONGJMP_NAME], ["${LONGJMP}"])

# Determine the sigsetjmp symbol.  If missing, then point to sigsetjmp_missing.
#
# Supported symbols:
#   sigsetjmp       POSIX, BSD libc (MacOS)
#   __sigsetjmp     glibc (Linux)
SIGSETJMP="sigsetjmp_missing"
for sigsetjmp_fn in sigsetjmp __sigsetjmp; do
  AX_FC_CHECK_BIND_C([${sigsetjmp_fn}], [
    SIGSETJMP=${sigsetjmp_fn}
    break
  ])
done
AC_DEFINE_UNQUOTED([SIGSETJMP_NAME], ["${SIGSETJMP}"])

# Validate the siglongjmp symbol
AX_FC_CHECK_BIND_C([siglongjmp],
  [SIGLONGJMP="siglongjmp"], [SETJMP="siglongjmp_missing"]
)
AC_DEFINE_UNQUOTED([SIGLONGJMP_NAME], ["${SIGLONGJMP}"])


# Verify the size of nonlocal jump buffer structs
# NOTE: This requires C compiler, but can it be done with a Fortran compiler?
AC_LANG_PUSH([C])

AX_MPI([], [AC_MSG_ERROR([Could not find MPI launcher.])])
AC_SUBST([CC], [$MPICC])

AC_CHECK_SIZEOF([jmp_buf], [], [#include <setjmp.h>])
AC_CHECK_SIZEOF([sigjmp_buf], [], [#include <setjmp.h>])

AC_LANG_POP([C])


# Prepare output
AC_SUBST([CPPFLAGS])
AC_CONFIG_FILES([Makefile:${srcdir}/ac/Makefile.in])
AC_OUTPUT
