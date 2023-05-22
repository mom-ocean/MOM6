# tc4 preprocessor configuration
AC_PREREQ([2.63])
AC_INIT([], [])

# Validate srdcir and configure input
AC_CONFIG_SRCDIR([gen_grid.F90])
AC_CONFIG_MACRO_DIR([../../ac/m4])


# Explicitly assume free-form Fortran
AC_LANG([Fortran])
AC_FC_SRCEXT([f90])

# We do not need MPI, but we want to emulate the executable used in MOM6
AX_MPI([], [AC_MSG_ERROR([Could not find MPI launcher.])])
AC_SUBST([FC], [$MPIFC])
AC_SUBST([LD], [$MPILD])


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


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
