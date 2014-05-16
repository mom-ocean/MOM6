#!/bin/csh
#------------------------------------------------------------------------------
#  MOM6_refineDiag.csh
#
#  DESCRIPTION: This is a script that is inteded to drive all the 
#               pre-postprocessing stages of data manipulations on analysis nodes.
#               It is intended to be called by frepp at the "refineDiag" stage 
#               which happens just before the components are post processed by frepp.
#               To make this happens a path to the (would be) script should appear 
#               in the <refineDiag> tag of the xmls, e.g., 
#      <refineDiag script="/nbhome/$USER/$(FRE_STEM)$(DEBUGLEVEL)/mom6/tools/analysis/untar_history_files.csh"/>
#               Note that the above script should exist when frepp is called.
#               This could be achieved by cloning the mom6 git repo in the <csh> section of the setup block 
#               in the corresponding the gfdl platfrom. E.g., 
#        <csh><![CDATA[
#           source $MODULESHOME/init/csh
#           module use -a /home/John.Krasting/local/modulefiles
#           module purge
#           module load jpk-analysis/0.0.4
#           module load $(FRE_VERSION)
###         The following clones the mom6 git repo which should contain all the pp scripts  
#           setenv NBROOT /nbhome/$USER/$(FRE_STEM)$(DEBUGLEVEL)
#           mkdir -p $NBROOT
#           cd $NBROOT
#           git clone /home/fms/git/ocean/mom6
#         ]]></csh>
#
#------------------------------------------------------------------------------
echo ""
echo "  -- begin MOM6_refineDiag.csh --  "
echo ""

#Unpack the history files and save them in the archive
echo ""
echo "  ---------- begin untar_history_files ----------  "
echo ""
     if (! -d $histDir/unpack) then
       mkdir -p $histDir/unpack
     endif
 
     foreach f (*.nc)
      if(! -e $histDir/unpack/${f} ) then
      gcp -v ${f} gfdl:${histDir}/unpack/
      endif
     end
echo "  ---------- end untar_history_files ----------  "
echo ""
echo "  ---------- begin yearly analysis ----------  "
echo ""
#
#Generate this year's analysis figures based on the unpacked history files
#
#Niki: Note that here we do not have any FRE environment to /nbhome
#      NBROOT should be set as an environ vatiable at the setup in the xml
#      setenv NBROOT /nbhome/$USER/$(FRE_STEM)$(DEBUGLEVEL)
set src_dir=$NBROOT/mom6/tools/analysis 
# The following variables are set by frepp, but frepp is not called at refineDiag stage, so wee need to set them here
set in_data_dir = $histDir/unpack
set descriptor = $name   
set out_dir = $NBROOT                  #Niki: How can we set this to frepp analysisdir /nbhome/...
set yr1 = $oname 
set yr2 = $oname 
set databegyr = $oname 
set dataendyr = $oname 
set datachunk = 1
set fremodule = fre            #Niki: How to set this to the caller version?
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

# List of variables
set varlist='temp salt'
# Determine available input variables and generate list of input files
set filelist="$yr1.ocean_static.nc $yr1.ocean_annual.nc" 

# Prepare workdir
set workdir=$TMPDIR/nnz_ocean_ts_annual
rm -rf $workdir
mkdir -p $workdir
# dmget and copy files
cd ${in_data_dir}
dmget ${filelist}
gcp ${filelist} $workdir
cd $workdir

# Generate include file common to all scripts
cat > files.txt << EOF
   descriptor = "${descriptor}"
   years = "${yr1}"
EOF

# Generate plots using python
cp ${src_dir}/MOM6_annual_analysis.py .
python MOM6_annual_analysis.py

# Copy output files to their destination
set dest_dir=${out_dir}/ocean_${yr1}/MOM6_test.ts/
if( ! -d $dest_dir ) mkdir -p $dest_dir
gcp -cd *.ps ${dest_dir}


echo "  ---------- end yearly analysis ----------  "

echo "  -- end   MOM6_refineDiag.csh --  "

exit
