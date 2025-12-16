# Copyright ACCESS-NRI and contributors. See the top-level LICENSE file for details.

if(MOM6_ACCESS3)
  ## Build MOM6 library for use in ACCESS3
  add_subdirectory(nuopc_cap)
endif()

if(MOM6_SOLO)
  ## Build MOM6 solo executable
  add_subdirectory(solo_driver)
endif()