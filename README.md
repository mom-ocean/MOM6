# MOM6 documentation

# Where to find information

During development, the MOM6 wiki is the primary place to find more information:
  http://wiki.gfdl.noaa.gov/index.php/MOM6

In particular, to setup your development-mode working directory there are extensive instructions at:
  http://wiki.gfdl.noaa.gov/index.php/MOM6_setup_instructions

If you are working from outside of GFDL, we have a minimal quick-start guide at the GitHub site
which you should be able to see if you have access to the GitHub repository:
  https://github.com/CommerceGov/NOAA-GFDL-MOM6/wiki

# What files are what

The top level directory structure groups source code and input files as follow:

| File/directory    | Purpose |
| --------------    | ------- |
| ```LICENSE.md```  | a copy of the Gnu general public license, version 3. |
| ```README.md```   | this file with basic pointers to more information |
| ```src/```        | contains the source code for MOM6 that is always compiled |
| ```config_src/``` | contains optional source code depending on mode and configuration such as dynamic-memory versus static, ocean-only versus coupled. |
| ```examples/```   | contains parameters, input data, paths to data, and some source code for static compiles. examples/ is sub-divided into four directories named for the style of compiled executable |
| ```pkg/```        | contains third party (non-MOM6 or FMS) code that is compiled into MOM6 |
| ```tools/```      | tools for working with MOM6 (not source code and not necessarily supported) |

## Sub-directories of ```examples/```

The examples are grouped by class of experiment for which dynamic executables can be shared:

| Directory                            | Experment class |
| ---------                            | --------------- |
| ```examples/solo_ocean```            | uses just MOM6 code |
| ```examples/ocean_SIS```             | uses just MOM6 and SIS code in coupled mode |
| ```examples/ocean_SIS2```            | uses just MOM6 and SIS2 code in coupled mode |
| ```examples/coupled_AM2_SIS```       | uses MOM6, SIS, LM2 and AM2 code ie. fully coupled |
| ```examples/coupled_AM2_LM3_SIS/```  | uses MOM6, SIS, LM3 and AM2 code |
| ```examples/coupled_AM2_LM3_SIS2/``` | uses MOM6, SIS2, ML3 and AM2 code |
