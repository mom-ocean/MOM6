# Copyright ACCESS-NRI and contributors. See the top-level LICENSE file for details.

function(add_fortran_target TARGET_NAME TYPE MOD_DIR)
  if(${TYPE} STREQUAL "LIB")
      add_library(${TARGET_NAME} ${ARGN})
  elseif(${TYPE} STREQUAL "EXE")
      add_executable(${TARGET_NAME} ${ARGN})
  else()
      message(FATAL_ERROR "Invalid TYPE '${TYPE}' for add_fortran_target. Use LIB or EXE.")
  endif()

  # Set Fortran module directory
  get_target_property(TGT_BUILD_DIR ${TARGET_NAME} BINARY_DIR)
  set_target_properties(${TARGET_NAME} PROPERTIES
      Fortran_MODULE_DIRECTORY "${TGT_BUILD_DIR}/${MOD_DIR}"
  )

  # Add the module directory to the include path
  target_include_directories(${TARGET_NAME} INTERFACE "$<BUILD_INTERFACE:${TGT_BUILD_DIR}/${MOD_DIR}>")
endfunction()
