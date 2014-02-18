#!/bin/csh -f
#------------------------------------
#PBS -N ocean_ts_annual
#PBS -l size=1
#PBS -l walltime=04:00:00
#PBS -r y
#PBS -j oe
#PBS -o
#PBS -q batch 
#----------------------------------

echo $0 $*

# variables set by frepp
set in_data_dir
set in_data_file
set descriptor
set out_dir
set yr1
set yr2
set databegyr
set dataendyr
set datachunk
set fremodule
set freanalysismodule = fre-analysis/test

# make sure valid platform and required modules are loaded
if (`gfdl_platform` == "hpcs-csc") then
   source $MODULESHOME/init/csh
   module use -a /home/fms/local/modulefiles /usr/local/paida/Modules
   module purge
   module load $fremodule
   module load $freanalysismodule
   module load gcc
   module load netcdf/4.2
   module load python/2.7.3
else
   echo "ERROR: invalid platform"
   exit 1
endif

# check again?
if (! $?FRE_ANALYSIS_HOME) then
   echo "ERROR: environment variable FRE_ANALYSIS_HOME not set."
   exit 1
endif

# Run script
set src_dir = ${out_dir}/mom6/tools/analysis
set script_bash = ${src_dir}/MOM6_ts_annual.bash
$script_bash -s $src_dir -i $in_data_dir -o $out_dir -d $descriptor -y $yr1,$yr2,$databegyr,$dataendyr,$datachunk

