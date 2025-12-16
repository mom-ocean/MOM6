# Copyright ACCESS-NRI and contributors. See the top-level LICENSE file for details.

add_fortran_target(mom6lib LIB mod STATIC)

target_link_libraries(mom6lib PRIVATE
  mom6shared
  Access3::nuopc_cap_share
  Access3::share
  Access3::timing
  Access3::cdeps-common
  ESMF::ESMF
)

target_compile_definitions(mom6lib PRIVATE CESMCOUPLED)

target_sources(mom6lib PRIVATE
  mom_cap_time.F90
  mom_surface_forcing_nuopc.F90
  ocn_comp_NUOPC.F90
  time_utils.F90
  mom_cap.F90
  mom_cap_methods.F90
  mom_cap_gtracer_flux.F90
  mom_ocean_model_nuopc.F90
  mom_cap_outputlog.F90
  mom_cap_profiling.F90
  mom_inline_mod.F90
)

# Install and export
set_target_properties(mom6lib PROPERTIES
  OUTPUT_NAME access-mom6lib
  EXPORT_NAME mom6lib
)
install(TARGETS mom6lib mom6shared
  EXPORT Mom6libTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT AccessMOM6_runtime NAMELINK_COMPONENT AccessMOM6_Development
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT AccessMOM6_Development
)
# Fortran module files are a special case, as currently there is no standard
# way of handling them in CMake
target_include_directories(mom6lib PUBLIC "$<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}/Mom6lib>")
  get_target_property(mom_moddir mom6lib Fortran_MODULE_DIRECTORY)
  install(FILES ${mom_moddir}/ocn_comp_nuopc.mod ${mom_moddir}/mom_cap_mod.mod
  DESTINATION ${CMAKE_INSTALL_MODULEDIR}/Mom6lib
  COMPONENT AccessMOM6_Development
)
install(EXPORT Mom6libTargets
  FILE Mom6libTargets.cmake
  NAMESPACE Access3::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mom6lib
)

configure_package_config_file(
  cmake/Mom6libConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/Mom6libConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mom6lib
)

install(FILES ${CMAKE_SOURCE_DIR}/cmake/FindNetCDF.cmake ${CMAKE_CURRENT_BINARY_DIR}/Mom6libConfig.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Mom6lib
  COMPONENT AccessMOM6_Development
)
