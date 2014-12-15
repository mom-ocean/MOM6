# MOM6

This is the MOM6 source code.

# Where to find information

Start at the [MOM6-examples wiki](https://github.com/CommerceGov/NOAA-GFDL-MOM6-examples/wiki) which has installation instructions.

# What files are what

The top level directory structure groups source code and input files as follow:

| File/directory    | Purpose |
| --------------    | ------- |
| ```LICENSE.md```  | a copy of the Gnu general public license, version 3. |
| ```README.md```   | this file with basic pointers to more information |
| ```src/```        | contains the source code for MOM6 that is always compiled |
| ```config_src/``` | contains optional source code depending on mode and configuration such as dynamic-memory versus static, ocean-only versus coupled. |
| ```pkg/```        | contains third party (non-MOM6 or FMS) code that is compiled into MOM6 |
