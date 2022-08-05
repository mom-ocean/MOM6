# Overview
This module is designed to be used in conjunction with the SmartSim and
SmartRedis libraries found at https://github.com/CrayLabs/. These
libraries are used to perform machine-learning inference and online
analysis using a Redis-based database.

An earlier implementation of these routines was used in Partee et al. [2022]:
"Using Machine Learning at scale in numerical simulations with SmartSim:
An application to ocean climate modeling" (doi.org/10.1016/j.jocs.2022.101707)
to predict eddy kinetic energy for use in the MEKE module. The additional
scripts and installation instructions for compiling MOM6 for this case can
be found at: https://github.com/CrayLabs/NCAR_ML_EKE/. The substantive
code in the new implementation is part of `MOM_MEKE.F90`.

# File description

- `MOM_database_comms` contains just method signatures and elements of the
  control structure that are imported elsewhere within the primary MOM6
  code. This includes: `dbcomms_CS_type`, `dbclient_type`, and `database_comms_init`

- `database_client_interface.F90` contains the methods for a communication client
  to transfer data and/or commands between MOM6 and a remote database. This is
  roughly based on the SmartRedis library, though only the methods that are most
  likely to be used with MOM6 are retained. This is to ensure that the API can be
  tested without requiring MOM6 users to compile in the the full library.
