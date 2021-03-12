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


# Default to symmetric grid
# NOTE: --enable is more properly used to add a feature, rather than to select
#   a compile-time mode, so this is not exactly being used as intended.
MEM_LAYOUT=${srcdir}/config_src/memory/dynamic_symmetric
AC_ARG_ENABLE([asymmetric],
  AS_HELP_STRING([--enable-asymmetric], [Use the asymmetric grid]))
AS_IF([test "$enable_asymmetric" = yes],
  [MEM_LAYOUT=${srcdir}/config_src/memory/dynamic_nonsymmetric])

# Default to solo_driver
DRIVER_DIR=${srcdir}/config_src/drivers/solo_driver
AC_ARG_WITH([driver],
  AS_HELP_STRING([--with-driver=coupled_driver|solo_driver], [Select directory for driver source code]))
AS_IF([test "x$with_driver" != "x"],
  [DRIVER_DIR=${srcdir}/config_src/drivers/${with_driver}])

# TODO: Rather than point to a pre-configured header file, autoconf could be
# used to configure a header based on a template.
#AC_CONFIG_HEADERS(["$MEM_LAYOUT/MOM_memory.h"])


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
AX_MPI([],
  [AC_MSG_ERROR([Could not find MPI launcher.])])


# Explicitly replace FC and LD with MPI wrappers
# NOTE: This is yet another attempt to manage the potential mismatches between
#   FC and MPIFC.  Without this step, the tests below would not use MPIFC.
AC_SUBST(FC, $MPIFC)
AC_SUBST(LD, $MPIFC)

# Confirm that FC can see the Fortran 90 MPI module.
AX_FC_CHECK_MODULE([mpi],
  [], [AC_MSG_ERROR([Could not find MPI Fortran module.])])


# netCDF configuration

# Search for the Fortran netCDF module, fallback to nf-config.
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

# FMS may invoke netCDF C calls, so we link to libnetcdf.
AC_LANG_PUSH([C])
AC_CHECK_LIB([netcdf], [nc_create], [], [
  AS_UNSET([ac_cv_lib_netcdf_nc_create])
  AC_PATH_PROG([NC_CONFIG], [nc-config])
  AS_IF([test -n "$NC_CONFIG"], [
    AC_SUBST([LDFLAGS],
      ["$LDFLAGS -L$($NC_CONFIG --libdir)"]
    )
    ], [AC_MSG_ERROR([Could not find nc-config.])]
  )
  AC_CHECK_LIB([netcdf], [nc_create], [], [
    AC_MSG_ERROR([Could not find libnetcdf.])
  ])
])
AC_LANG_POP([C])

# NOTE: We test for nf_create, rather than nf90_create, because AX_FC_CHECK_LIB
#   is currently not yet able to properly probe inside modules.
#   Testing of the nf90_* functions will require a macro update.
# NOTE: nf-config does not have a --libdir flag, so we use --prefix and assume
#   that libraries are in the $prefix/lib directory.

# Link to Fortran netCDF library, netcdff
AX_FC_CHECK_LIB([netcdff], [nf_create], [], [], [
  AS_UNSET([ax_fc_cv_lib_netcdff_nf_create])
  AC_PATH_PROG([NF_CONFIG], [nf-config])
  AS_IF([test -n "$NF_CONFIG"], [
    AC_SUBST([LDFLAGS],
      ["$LDFLAGS -L$($NF_CONFIG --prefix)/lib"]
    )
    ], [AC_MSG_ERROR([Could not find nf-config.])]
  )
  AX_FC_CHECK_LIB([netcdff], [nf_create], [], [], [
    AC_MSG_ERROR([Could not find libnetcdff.])
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


# Search for mkmf build tools
AC_PATH_PROG([LIST_PATHS], [list_paths])
AS_IF([test -z "$LIST_PATHS"], [
  AC_PATH_PROG([LIST_PATHS], [list_paths], [], ["$PATH:${srcdir}/ac/deps/bin"])
  AS_IF([test -z "$LIST_PATHS"],
    [AC_MSG_ERROR([Could not find list_paths.])],
    [AC_SUBST(PATH, ["$PATH:${srcdir}/ac/deps/bin"])])
  ]
)

AC_PATH_PROG([MKMF], [mkmf])
AS_IF([test -z "$MKMF"], [
  AC_PATH_PROG([MKMF], [mkmf], [], ["$PATH:${srcdir}/ac/deps/bin"])
  AS_IF([test -z "$MKMF"],
    [AC_MSG_ERROR([Could not find mkmf.])],
    [AC_SUBST(PATH, ["$PATH:${srcdir}/ac/deps/bin"])])
  ]
)


# NOTE: MEM_LAYOUT unneeded if we shift to MOM_memory.h.in template
AC_CONFIG_COMMANDS([path_names],
  [list_paths -l \
    ${srcdir}/src \
    ${srcdir}/config_src/infra/FMS1 \
    ${srcdir}/config_src/ext* \
    ${DRIVER_DIR} \
    ${MEM_LAYOUT}
], [MEM_LAYOUT=$MEM_LAYOUT DRIVER_DIR=$DRIVER_DIR])


AC_CONFIG_COMMANDS([Makefile.mkmf],
  [mkmf -p MOM6 -m Makefile.mkmf path_names])


# Prepare output
AC_SUBST(CPPFLAGS)
AC_CONFIG_FILES([Makefile:${srcdir}/ac/Makefile.in])
AC_OUTPUT
