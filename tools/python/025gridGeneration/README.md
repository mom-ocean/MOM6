The Makefile in this directory manages the complete work flow to create
the grids, topography and some oceanographic fields needed to run OM4_025.

The environment must be setup outside of make. At GFDL do:

  source $MODULESHOME/init/csh
  module use -a /home/fms/local/modulefiles
  module load netcdf/4.2 intel_compilers
  module load nco/4.3.1
  module load python
