[![Build Status](https://travis-ci.org/NOAA-GFDL/MOM6.svg?branch=dev/master)](https://travis-ci.org/NOAA-GFDL/MOM6)
[![Read The Docs Status](https://readthedocs.org/projects/mom6/badge/?badge=latest)](http://mom6.readthedocs.io/)
[![codecov](https://codecov.io/gh/NOAA-GFDL/MOM6/branch/dev%2Fmaster/graph/badge.svg)](https://codecov.io/gh/NOAA-GFDL/MOM6)

# MOM6

This is the MOM6 source code.


# Where to find information

Start at the [MOM6-examples wiki](https://github.com/NOAA-GFDL/MOM6-examples/wiki) which has installation instructions.

[Source code documentation](http://mom6.readthedocs.io/) is hosted on read the docs.


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
