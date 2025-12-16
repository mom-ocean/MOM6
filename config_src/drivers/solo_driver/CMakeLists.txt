# Copyright ACCESS-NRI and contributors. See the top-level LICENSE file for details.

add_fortran_target(MOM6 EXE mod)

target_link_libraries(MOM6 PRIVATE
  mom6shared
)

target_sources(MOM6 PRIVATE
  MOM_surface_forcing.F90
  MESO_surface_forcing.F90
  user_surface_forcing.F90
  MOM_driver.F90
)

set_target_properties(MOM6 PROPERTIES
  LINKER_LANGUAGE Fortran
  OUTPUT_NAME mom6-solo
)

# Install
install(TARGETS MOM6
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)