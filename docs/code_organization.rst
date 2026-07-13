Organization of the code
========================

The MOM6 source code is divided into a tree of directories (folders) to group
related code (e.g. `src/core`) or similar modules (e.g.
`src/parameterizations/vertical`).

The highest level of directories are:

`src/`
  Code underneath `src/` is always required and compiled.

`config_src/`
  Under `config_src` are various drivers and memory configuration sources that
  can only be compiled in limited configurations. See :ref:`config_src`

`pkg/`
  Packages (source code) from other sources/parties only some of which might
  be used. We include the entire package as a sub-module but use
  symbolic-links to extract the parts the MOM6 uses.

`docs/`
  The directory that contains this documentation, namely that beyond the
  in-code API documentation. Some of the files are configuration files
  needed for running doxygen and sphinx. Most documentation in this folder
  is in the form of reStructuredText (.rst) files.

`.testing/`
  A directory for running tests on MOM6. These are some of the
  smaller/simpler examples that MOM6 can run out of the box, without
  large netCDF files.

The directory tree is::

  MOM6
  +-- config_src
  |   +-- drivers
  |   |   +-- FMS_cap                  # GFDL coupler interface
  |   |   +-- ice_solo_driver          # Ice-only standalone
  |   |   +-- nuopc_cap                # NUOPC/CESM coupling
  |   |   +-- solo_driver              # Ocean-only standalone
  |   |   +-- timing_tests             # Performance benchmarks
  |   |   +-- unit_tests               # Unit test executables
  |   +-- external                     # Null hooks for optional components
  |   +-- infra                        # Framework interface (FMS1/FMS2 wrappers)
  |   |   +-- FMS1                     # FMS1 wrappers
  |   |   +-- FMS2                     # FMS2 wrappers
  |   +-- memory
  |       +-- dynamic_nonsymmetric
  |       +-- dynamic_symmetric
  +-- docs
  +-- pkg
  |   +-- CVMix-src                    # Community Vertical Mixing
  |   +-- GSW-Fortran                  # TEOS-10 Gibbs Seawater
  +-- src
  |   +-- ALE                          # Vertical remapping/regridding
  |   +-- core                         # Dynamical core
  |   +-- diagnostics                  # Diagnostic calculations
  |   +-- equation_of_state            # EOS implementations
  |   +-- framework                    # Infrastructure (diagnostics, I/O, parsing, domains)
  |   +-- ice_shelf                    # Ice shelf dynamics
  |   +-- initialization               # Grid/state initialization
  |   +-- ocean_data_assim             # Data assimilation
  |   +-- parameterizations
  |   |   +-- lateral                  # Lateral parameterizations
  |   |   +-- vertical                 # Vertical mixing
  |   +-- tracer                       # Tracer transport and specific tracers
  |   +-- user                         # Idealized configuration initialization
  +-- .testing
      +-- tc0                          # Unit tests
      +-- tc1 / tc1.a / tc1.b          # Benchmark configurations
      +-- tc2 / tc2.a                  # ALE with tides / sigma-coordinate
      +-- tc3                          # Open boundary conditions
      +-- tc4                          # Sponges and I/O initialization

.. _config_src:

`config_src/`
-------------

`memory/dynamic_nonsymmetric/`, `memory/dynamic_symmetric/`
  One or none of `config_src/memory/dynamic_nonsymmetric/` or
  `config_src/memory/dynamic_symmetric/` can be included at compile time. If neither
  is used then a `MOM_memory.h` file specific to the model configuration must be
  present - this is known as a "static" compile with fixed layout and domain shape.

`external/`
  Contains "null" modules providing the API to optional components to use
  with MOM6. Currently available are ocean data assimilation (`ODA_hooks`) and
  the GFDL ocean bio-geochemistry model (`GFDL_ocean_BGC`). When building
  MOM6 in stand-alone ocean-only mode these modules should be compiled in.
  To use the actual ODA or BGC, add the appropriate source to the search
  paths .

`infra/`
  Contains MOM6-specific thin wrappers to all of the FMS types and routines that
  are used by MOM6.  The code in this directory should only be called by the
  infrastructure-agnostic code in src/framework.

`drivers/ice_solo_driver/`
  This driver produces a stand-alone ice-shelf executable that steps the
  ice-shelf model without any ocean dynamics.

`drivers/solo_driver/`
  This driver produces an ocean-only executable with no other coupled
  components (no sea-ice, no atmosphere, etc.). It is the simplest
  configuration and fastest to compile and thus used for a lot of testing.

`drivers/FMS_cap/`
  This driver provides an interface for the GFDL coupler to call. It requires
  compiling MOM6 along with at least a sea-ice model and possibly all other
  components in a coupled model.

`drivers/nuopc_cap/`
  This driver provides a NUOPC-compliant interface for coupling MOM6 within
  CESM or other NUOPC-based coupled systems.

`drivers/unit_tests/`
  Unit test executables for testing individual MOM6 components in isolation.

`drivers/timing_tests/`
  Performance benchmark executables for profiling MOM6 routines.

.. _src:

`src/`
------

`core/`
  The dynamical core modules (except for the ALE remapping/regridding).

`ALE/`
  Functions for remapping from between arbitrary vertical grids
  and generating grids.

`diagnostics/`
  Some diagnostic calculations

`equation_of_state/`
  Various equations of state (linear; Wright, 1997; TEOS-10; ...).

`framework/`
  Low-level wrappers for communication, diagnostics management, parsing
  of input parameters, time management, CPU clocks.

`ice_shelf/`
  Ice shelf dynamics and thermodynamics.

`initialization/`
  Initialization of the horizontal grid, vertical coordinate, and the state.

`ocean_data_assim/`
  Data assimilation interfaces.

`parameterizations/lateral`
  Sub-grid scale parameterization with fluxes primarily oriented in the
  lateral direction.

`parameterizations/vertical`
  Sub-grid scale parameterization with fluxes primarily oriented in the
  vertical direction, including the top and bottom boundary layer schemes.

`tracer/`
  Everything to do with tracers, including advection and isopycnal stirring.

`user/`
  Initialization and forcing for specific (coded) configurations.
