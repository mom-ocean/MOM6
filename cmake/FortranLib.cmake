# Copyright ACCESS-NRI and contributors. See the top-level LICENSE file for details.

function(add_fortran_library LIB MOD_DIR)
    add_library(${LIB} ${ARGN})

    get_target_property(LIB_DIR ${LIB} BINARY_DIR)
    set_target_properties(${LIB} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_DIR}/${MOD_DIR})

    target_include_directories(${LIB} INTERFACE "$<BUILD_INTERFACE:${LIB_DIR}/${MOD_DIR}>")
  endfunction(add_fortran_library)
