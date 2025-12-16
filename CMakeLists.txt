# Copyright ACCESS-NRI and contributors. See the top-level LICENSE file for details.

cmake_minimum_required(VERSION 3.18)

#[==============================================================================[
#                         Basic project definition                              #
#]==============================================================================]
project(MOM6
        DESCRIPTION "Modular Ocean Model 6"
        HOMEPAGE_URL https://github.com/ACCESS-NRI/MOM6
        LANGUAGES C Fortran)

#[==============================================================================[
#                           Project configuration                               #
#]==============================================================================]
message(STATUS "Build options")

option(MOM6_ACCESS3    "Build MOM6 library for use with ACCESS3"  OFF)
option(MOM6_SOLO       "Build MOM6 solo executable"               ON)
option(MOM6_OPENMP     "Enable OpenMP threading"                  OFF)
option(MOM6_ASYMMETRIC "Use MOM asymmetric memory"                OFF)
message(STATUS "  - MOM6_ACCESS3     ${MOM6_ACCESS3}")
message(STATUS "  - MOM6_SOLO        ${MOM6_SOLO}")
message(STATUS "  - MOM6_OPENMP      ${MOM6_OPENMP}")
message(STATUS "  - MOM6_ASYMMETRIC  ${MOM6_ASYMMETRIC}")

if(NOT MOM6_ACCESS3 AND NOT MOM6_SOLO)
  message(FATAL_ERROR "One of MOM6_ACCESS3 or MOM6_SOLO must be enabled")
endif()

#[==============================================================================[
#                           Project configuration                               #
#]==============================================================================]

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(GNUInstallDirs)
include(FortranLib)
include(CMakePackageConfigHelpers)

# Fortran compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbacktrace -fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none -fdefault-real-8 -fdefault-double-8 -g -grecord-gcc-switches")
  if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif()
  set(CMAKE_Fortran_FLAGS_RELEASE "-O")
  set(CMAKE_Fortran_FLAGS_DEBUG "-Wall -Og -ffpe-trap=zero,overflow -fcheck=bounds")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qno-opt-dynamic-align  -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -fp-model precise -r8 -g -grecord-gcc-switches")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check noarg_temp_created -check all") #  -fpe0 see https://github.com/ACCESS-NRI/MOM6/issues/39
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

# C compiler flags
if(CMAKE_C_COMPILER_ID MATCHES "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -g -grecord-gcc-switches")
  set(CMAKE_C_FLAGS_RELEASE "-O")
  set(CMAKE_C_FLAGS_DEBUG "-Wall -Og -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=bounds")
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -traceback -qno-opt-dynamic-align -fp-model precise -std=gnu99 -g -grecord-gcc-switches")
  set(CMAKE_C_FLAGS_RELEASE "-O2")
  set(CMAKE_C_FLAGS_DEBUG "-O0")
else()
  message(WARNING "C compiler with ID ${CMAKE_C_COMPILER_ID} will be used with CMake default options")
endif()

# See https://github.com/ACCESS-NRI/GFDL-generic-tracers/issues/29
add_compile_definitions(_USE_GENERIC_TRACER) # _USE_MOM6_DIAG)

## Fortran modules path; currently this is simply set to be the include dir
set(CMAKE_INSTALL_MODULEDIR ${CMAKE_INSTALL_INCLUDEDIR}
  CACHE STRING
  "Fortran module installation path (Not a cmake native variable)"
)

#[==============================================================================[
#                              External packages                                #
#]==============================================================================]

find_package(fms COMPONENTS R8 REQUIRED)
find_package(GFDLGTracers REQUIRED)

if(NOT TARGET NetCDF::NetCDF_Fortran OR NOT TARGET NetCDF::NetCDF_C) #fms provides NETCDF inteface
  message(FATAL_ERROR "NetCDF interface missing from FMS package")
endif()

if(MOM6_ACCESS3)
  find_package(Access3Share REQUIRED cdeps timing share nuopc_cap_share)
  if(NOT TARGET ESMF::ESMF)
    message(FATAL_ERROR "ESMF interface missing from Access3Share package")
  endif()
endif()

if(MOM6_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

#[==============================================================================[
#                             Main definitions                                  #
#]==============================================================================]

## Build MOM6 shared library
add_subdirectory(src)

## Build MOM6 ACCESS3 library and/or MOM6 solo executable
add_subdirectory(config_src)
