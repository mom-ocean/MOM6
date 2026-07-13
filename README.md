[![Read The Docs Status](https://readthedocs.org/projects/mom6/badge/?version=main)](https://mom6.readthedocs.io/en/main/?badge=main)
[![codecov](https://codecov.io/gh/NOAA-GFDL/MOM6/branch/dev/gfdl/graph/badge.svg?token=uF8SVydCdp)](https://codecov.io/gh/NOAA-GFDL/MOM6)

# MOM6

MOM6 (Modular Ocean Model, version 6) is a next-generation open-development-source ocean model developed by a consortium of institutions. It uses a modern Fortran codebase solving the primitive equations for ocean dynamics on an Arakawa C-grid.

- Arbitrary Lagrangian-Eulerian (ALE) vertical coordinate
- Boussinesq and non-Boussinesq modes
- Flexible equation of state (Wright, TEOS-10, linear, Roquet (TEOS-10), ...)
- Comprehensive parameterization library (ePBL, KPP, lateral mixing, tidal forcing)
- Coupled to SIS2 or CICE (sea ice), ice shelves, and Earth system models via the FMS or NUOPC couplers, or run in stand-alone ocean-only configurations
- Dimensional unit scaling for consistency testing


# Where to find information

- [MOM6-examples wiki](https://github.com/NOAA-GFDL/MOM6-examples/wiki) -- installation instructions and tutorials
- [Source code documentation](https://mom6.readthedocs.io/en/main/) -- hosted on Read the Docs
- [Developers wiki](https://github.com/NOAA-GFDL/MOM6/wiki) -- developer guides and conventions
- [Developer workflow](https://github.com/NOAA-GFDL/MOM6/wiki/Developer-workflow)
- [Runtime parameter system](https://github.com/NOAA-GFDL/MOM6/wiki/MOM6-run-time-parameter-system)
- [Repository policies](https://github.com/NOAA-GFDL/MOM6-examples/wiki/MOM6-repository-policies)
- [MOM6 forum](https://bb.cgd.ucar.edu/cesm/forums/mom6.148/)
- [CVMix](https://github.com/CVMix/CVMix-src) -- Community Vertical Mixing
- [TEOS-10 (GSW)](http://www.teos-10.org/) -- Gibbs Seawater
- [GOTM](https://gotm.net/) -- General Ocean Turbulence Model


# What files are what

The top level directory structure groups source code and input files as follow:

| File/directory    | Purpose |
| --------------    | ------- |
| ```LICENSE.md```  | A copy of the Gnu lesser general public license, version 3. |
| ```README.md```   | This file with basic pointers to more information. |
| ```src/```        | Contains the source code for MOM6 that is always compiled. |
| ```config_src/``` | Contains optional source code depending on mode and configuration such as dynamic-memory versus static, ocean-only versus coupled. |
| ```pkg/```        | Contains third party (non-MOM6 or FMS) code that is compiled into MOM6. |
| ```docs/```       | Workspace for generated documentation.  See [docs/README.md](docs/README.md) |
| ```.testing/```   | Contains the verification test suite.  See [.testing/README.md](.testing/README.md) |
| ```ac/```         | Contains the autoconf build configuration files. See [ac/README.md](ac/README.md) |


# Quick start guide

To quickly get started and build an ocean-only MOM6 executable, see the
[autoconf README](ac/README.md).

For setting up an experiment, or building an executable for coupled modeling,
consult the [MOM6-examples wiki](https://github.com/NOAA-GFDL/MOM6-examples/wiki).


# Disclaimer

The United States Department of Commerce (DOC) GitHub project code is provided
on an "as is" basis and the user assumes responsibility for its use. DOC has
relinquished control of the information and no longer has responsibility to
protect the integrity, confidentiality, or availability of the information. Any
claims against the Department of Commerce stemming from the use of its GitHub
project will be governed by all applicable Federal law. Any reference to
specific commercial products, processes, or services by service mark,
trademark, manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce. The
Department of Commerce seal and logo, or the seal and logo of a DOC bureau,
shall not be used in any manner to imply endorsement of any commercial product
or activity by DOC or the United States Government.

This project code is made available through GitHub but is managed by NOAA-GFDL
at https://gitlab.gfdl.noaa.gov.
