#!/bin/bash

# Environment variables set by gitlab (the CI environment)
if [ -z $JOB_DIR ]; then
  echo Environment variable "$"JOB_DIR should be defined to point to a unique directory for these scripts to use.
  echo '$JOB_DIR is derived from $CI_PIPELINE_ID in MOM6/.gitlab-ci.yml'
  echo 'To use interactively try:'
  echo '  JOB_DIR=tmp' $0 $@
  exit 911
fi
if [ -z $CI_PROJECT_DIR ]; then
  echo Environment variable "$"CI_PROJECT_DIR should be defined and point to where gitlab has cloned the MOM6 repository for this pipeline.
  echo 'To use interactively try:'
  echo '  CI_PROJECT_DIR=.' $0 $@
  exit 911
else
  CI_PROJECT_DIR=`realpath $CI_PROJECT_DIR`
fi
if [ -z $CI_COMMIT_SHA ]; then
  echo Environment variable "$"CI_COMMIT_SHA should be defined and indicate the MOM6 commit to used in this pipeline.
  echo 'To use interactively try:'
  echo '  CI_COMMIT_SHA=`git rev-parse HEAD`' $0 $@
  exit 911
fi

# Use CI=true to enable the gitlab folding

set -e # Stop if we encounter an error

# Environment variables that can be set outside
STATS_REPO_URL="${STATS_REPO_URL:-https://gitlab.gfdl.noaa.gov/ogrp/Gaea-stats-MOM6-examples.git}"
STATS_REPO_BRANCH="${STATS_REPO_BRANCH:-dev/gfdl}"
CONFIGS_DIR="${CONFIGS_DIR:-MOM6-examples}"
CONFIGS_REPO_BRANCH="${CONFIGS_REPO_BRANCH:-$STATS_REPO_BRANCH}"

# Global variables derived from the above
DRYRUN=
STATS_REPO=$(basename $STATS_REPO_URL)
STATS_REPO_DIR=$(basename $STATS_REPO .git)

# Static variables
RED=$'\033[1;31m'
GRN=$'\033[1;32m'
OFF=$'\e[m'

# Print the start of a fold in the log
section-start () {
  echo -e "\e[0Ksection_start:`date +%s`:$1[collapsed=true]\r\e[0K$2"
}

# Print the start of a fold in the log but not collapsed
section-start-open () {
  echo -e "\e[0Ksection_start:`date +%s`:$1[collapsed=false]\r\e[0K$2"
}

# Print the end of a fold in the log
section-end () {
  echo -e "\e[0Ksection_end:`date +%s`:$1\r\e[0K"
}

# Create $JOB_DIR and clean out any prior work-spaces
# Location: run in MOM6 directory
clean-job-dir () {
  section-start clean-job-dir "Cleaning $JOB_DIR directory"
  if [ ! $DRYRUN ] ; then
    #NOT USED? cd $CI_PROJECT_DIR
    #NOT USED? git submodule init ; git submodule update
    echo Job directory set to $JOB_DIR
    mkdir -p $JOB_DIR
    cd $JOB_DIR
    test -d $STATS_REPO_DIR && rm -rf $STATS_REPO_DIR # In case we are re-running this stage
  fi
  section-end clean-job-dir
}

# Create the full work space starting at the regression repository (usually Gaea-stats-MOM6-examples)
# Location: run in MOM6 directory
create-job-dir () {
  section-start create-job-dir "Creating and populating $JOB_DIR"
  if [ ! $DRYRUN ] ; then
    mkdir -p $JOB_DIR
    cd $JOB_DIR
    git clone $STATS_REPO_URL $STATS_REPO_DIR
    cd $STATS_REPO_DIR
    git checkout $STATS_REPO_BRANCH
    git submodule update --init
    cd $CONFIGS_DIR
    git checkout $CONFIGS_REPO_BRANCH
    git submodule init
    git submodule set-url src/MOM6 $CI_PROJECT_DIR/.git
    git submodule update --recursive --jobs 8
    (cd src/MOM6 ; git checkout $CI_COMMIT_SHA) # Get commit to be tested
    (cd src/MOM6 ; git submodule update --recursive --init)
    make -f tools/MRS/Makefile.clone clone_gfdl -j # Extras and link to datasets
    bash tools/MRS/generate_manifest.sh . tools/MRS/excluded-expts.txt > manifest.mk
    mkdir -p results
  fi
  section-end create-job-dir
}

# Create a copy of the configurations working directory
# Location: run in MOM6 directory
copy-test-space () {
  if [ -z $1 ]; then echo "copy-test-space needs an argument" ; exit 911 ; fi
  section-start copy-test-space-$1 "Copying $CONFIGS_DIR for $1"
  if [ ! $DRYRUN ] ; then
    COPIED_DIR=tmp-$CONFIGS_DIR-$1
    cd $JOB_DIR/$STATS_REPO_DIR
    git clone -s $CONFIGS_DIR/.git $COPIED_DIR
    cd $COPIED_DIR
    ln -s ../$CONFIGS_DIR/{build,results,.datasets} .
    cp ../$CONFIGS_DIR/manifest.mk .
  fi
  section-end copy-test-space-$1
}

# Build a group of executables using the tools/MRS/Makefile.build template
# Location: run in MOM6 directory
mrs-compile () {
  if [ -z $1 ]; then echo "mrs-compile needs an argument" ; exit 911 ; fi
  section-start mrs-compile-$1 "Compiling target $1"
  if [ ! $DRYRUN ] ; then
    cd $JOB_DIR/$STATS_REPO_DIR/$CONFIGS_DIR
    time make -f tools/MRS/Makefile.build $1 -s -j
  fi
  section-end mrs-compile-$1
}

# Build an ocean-only executable without intermediate libraries
# Location: run in MOM6 directory
nolibs-ocean-only-compile () {
  if [ -z $1 ]; then echo "nolibs-ocean-only-compile needs an argument" ; exit 911 ; fi
  section-start nolibs-ocean-only-compile-$1 "Compiling ocean-only $1 executable"
  if [ ! $DRYRUN ] ; then
    cd $JOB_DIR/$STATS_REPO_DIR/$CONFIGS_DIR
    mkdir -p build-ocean-only-nolibs-$1
    cd build-ocean-only-nolibs-$1
    make -f ../tools/MRS/Makefile.build ./$1/env BUILD=. ENVIRON=../../environ -s
    ../src/mkmf/bin/list_paths -l ../src/MOM6/config_src/{drivers/solo_driver,memory/dynamic_symmetric,infra/FMS1,ext*} ../src/MOM6/src ../src/FMS1
    sed -i '/FMS1\/.*\/test_/d' path_names
    ../src/mkmf/bin/mkmf -t ../src/mkmf/templates/ncrc5-$1.mk -p MOM6 -c"-Duse_libMPI -Duse_netCDF" path_names
    (source $1/env ; make NETCDF=3 REPRO=1 MOM6 -s -j)
  fi
  section-end nolibs-ocean-only-compile-$1
}

# Build an ocean-ice executable without intermediate libraries
# Location: run in MOM6 directory
nolibs-ocean-ice-compile () {
  if [ -z $1 ]; then echo "nolibs-ocean-ice-compile needs an argument" ; exit 911 ; fi
  section-start nolibs-ocean-ice-compile-$1 "Compiling ocean-ice $1 executable"
  if [ ! $DRYRUN ] ; then
    cd $JOB_DIR/$STATS_REPO_DIR/$CONFIGS_DIR
    mkdir -p build-ocean-ice-nolibs-$1
    cd build-ocean-ice-nolibs-$1
    make -f ../tools/MRS/Makefile.build ./$1/env BUILD=. ENVIRON=../../environ -s
    ../src/mkmf/bin/list_paths -l ../src/MOM6/config_src/{drivers/FMS_cap,memory/dynamic_symmetric,infra/FMS1,ext*} ../src/MOM6/src ../src/SIS2/*src ../src/icebergs/src ../src/{FMS1,coupler,ice_param,land_null,atmos_null}
    sed -i '/FMS1\/.*\/test_/d' path_names
    ../src/mkmf/bin/mkmf -t ../src/mkmf/templates/ncrc5-$1.mk -p MOM6 -c"-Duse_libMPI -Duse_netCDF -D_USE_LEGACY_LAND_ -Duse_AM3_physics" path_names
    (source $1/env ; make NETCDF=3 REPRO=1 MOM6 -s -j)
  fi
  section-end nolibs-ocean-ice-compile-$1
}

# Internal function to clean up stats files
# Args: list of top level directories to scan
clean-stats () {
  find $@ -name "*.stats.*[a-z][a-z][a-z]" -delete
}

# Internal function to clean up param files
# Args: list of top level directories to scan
clean-params () {
  find $@ -name "*_parameter_doc.*" -delete
  find $@ -name "*available_diags*" -delete
}

# Internal function to check for core files
# Args: list of top level directories to scan
check-for-core-files () {
  EXIT_CODE=0
  find $@ -name core -type f | grep . && EXIT_CODE=1
  if [[ $EXIT_CODE -gt 0 ]]
  then
    echo "Error: core files found!"
    exit 911
  fi
}

# Internal function to clean up core files (needed for re-running)
# Args: list of top level directories to scan
clean-core-files () {
  find $@ -name core -type f -delete
}

# Internal function to run a sub-suite and copy results to storage
# Args:
# $1 is compiler (gnu, intel, pgi, ...)
# $2 is sub-suite (_all, _ocean_only, _static_ocean_only, ...)
# $3 is MEMORY macro (dynamic_symmetric, dynamic_nonsymmetric, static)
# $4 is MODE macro (repro, debug)
# $5 is LAYOUT macro (def, alt)
mrs-run-sub-suite () {
  if [ "$#" -ne 5 ]; then echo "mrs-run-sub-suite needs 5 arguments" ; exit 911 ; fi
  section-start mrs-run-sub-suite-$1-$2-$3-$4-$5 "Running target $1-$2-$3-$4-$5"
  EXP_GROUPS=`grep / manifest.mk | sed 's:/.*::' | uniq`
  clean-stats $EXP_GROUPS
  clean-params $EXP_GROUPS
  clean-core-files $EXP_GROUPS
  if [[ "$3" == *"_nonsym"* ]]; then
    set -e
    time make -f tools/MRS/Makefile.run ocean_only/circle_obcs/ocean.stats.$1 MEMORY=${3/_nonsym/_sym} MODE=$4 LAYOUT=$5 -s -j
  fi
  set -e
  time make -f tools/MRS/Makefile.run $1_$2 MEMORY=$3 MODE=$4 LAYOUT=$5 -s -j
  tar cf - `find $EXP_GROUPS -name "*.stats.*[a-z][a-z][a-z]"` | tar --one-top-level=results/$1-$2-$3-$4-$5-stats -xf -
  tar cf - `find $EXP_GROUPS -name "*_parameter_doc.*" -o -name "*available_diags*"` | tar --one-top-level=results/$1-$2-$3-$4-$5-params -xf -
  check-for-core-files $EXP_GROUPS
  section-end mrs-run-sub-suite-$1-$2-$3-$4-$5
}

# Internal function to run restarts on a sub-suite and copy results to storage
# Args:
#   $1 is compiler (gnu, intel, pgi, ...)
#   $2 is sub-suite (_all, _ocean_only, _static_ocean_only, ...)
#   $3 is MEMORY macro (dynamic_symmetric, dynamic_nonsymmetric, static)
#   $4 is MODE macro (repro, debug)
#   $5 is LAYOUT macro (def, alt)
mrs-run-restarts-sub-suite () {
  if [ "$#" -ne 5 ]; then echo "mrs-run-restarts-sub-suite needs 5 arguments" ; exit 911 ; fi
  section-start mrs-run-restarts-sub-suite-$1-$2-$3-$4-$5 "Running target $1-$2-$3-$4-$5"
  clean-stats $2
  clean-core-files $2
  time make -f tools/MRS/Makefile.restart $1_$2 MEMORY=$3 MODE=$4 LAYOUT=$5 -s -j RESTART_STAGE=01
  check-for-core-files $2
  time make -f tools/MRS/Makefile.restart $1_$2 MEMORY=$3 MODE=$4 LAYOUT=$5 -s -j RESTART_STAGE=02
  check-for-core-files $2
  time make -f tools/MRS/Makefile.restart $1_$2 MEMORY=$3 MODE=$4 LAYOUT=$5 -s -j RESTART_STAGE=12
  check-for-core-files $2
  section-end mrs-run-restarts-sub-suite-$1-$2-$3-$4-$5
}

# Run a suite of experiments
#   $1 - compiler brand
#   $2 - any combination of "SNLDTR"
#        S = symmetric
#        N = non-symmetric
#        L = layout
#        D = debug
#        R = restarts
run-suite () {
  if [ "$#" -ne 2 ]; then echo "run-suite needs 2 arguments" ; exit 911 ; fi
  section-start run-suite-$1-$2 "Running suite for $1-$2"
  WORK_DIR=tmp-$CONFIGS_DIR-$1
  rm -f $JOB_DIR/CI-BATCH-SUCCESS-$1-$2
  set -e
  set -v

  pushd $JOB_DIR/$STATS_REPO_DIR/$WORK_DIR > /dev/null
  if [[ "$2" =~ "S" ]]; then # Symmetric
    mrs-run-sub-suite $1 all dynamic_symmetric repro def
  fi
  if [[ "$2" =~ "N" ]]; then # Non-symmetric
    mrs-run-sub-suite $1 all dynamic_nonsymmetric repro def
  fi
  if [[ "$2" =~ "L" ]]; then # Layout
    mrs-run-sub-suite $1 all dynamic_symmetric repro alt
  fi
  if [[ "$2" =~ "D" ]]; then # Debug
    mrs-run-sub-suite $1 ocean_only dynamic_symmetric debug def
  fi
  if [[ "$2" =~ "T" ]]; then # sTatic
    mrs-run-sub-suite $1 static_ocean_only static repro def
  fi
  popd > /dev/null
  if [[ "$2" =~ "R" ]]; then # Restarts
    pushd $JOB_DIR/$STATS_REPO_DIR/$WORK_DIR-rst > /dev/null
    mrs-run-restarts-sub-suite $1 ocean_only dynamic_symmetric repro def
    mrs-run-restarts-sub-suite $1 ice_ocean_SIS2 dynamic_symmetric repro def
    popd > /dev/null
  fi

  # Indicate all went well
  touch $JOB_DIR/CI-BATCH-SUCCESS-$1-$2

  section-end run-suite-$1-$2
}

# Test the value of stats files. All files in results/ are checked for in regressions/. It is assumed
# missing files are intended and failed runs were caught earlier in the CI process.
# Args:
#   $1 is path of results to test (relative to $STATS_REPO_DIR)
#   $2 is path of correct results to test against (relative to $STATS_REPO_DIR)
compare-stats () {
  if [ "$#" -ne 2 ]; then echo "compare-stats needs 2 arguments" ; exit 911 ; fi
  section-start-open compare-stats-$1-$2-$3-$4-$5 "Checking stats for '$1' against '$2'"
  # This checks that any file in the results directory is exactly the same as in regressions/
  ( cd $JOB_DIR/$STATS_REPO_DIR/$1 ; md5sum `find * -type f` ) | ( cd $JOB_DIR/$STATS_REPO_DIR/$2 ; md5sum -c ) 2>&1 | sed "s/ OK/$GRN&$OFF/;s/ FAILED/$RED&$OFF/;s/WARNING/$RED&$OFF/"
  FAIL=${PIPESTATUS[1]}
  if [ ! $FAIL == 0 ]; then
    exit 911
  fi
  section-end compare-stats-$1-$2-$3-$4-$5
}

# Test the value of stats files for a class of run
#   $1 - compiler brand
#   $2 - any combination of "SNLDTR"
#        S = symmetric
#        N = non-symmetric
#        L = layout
#        D = debug
#        T = static
#        R = restarts
#
# Many tests are tested against the "dynamic_symmetric repro" suite which must also have been run.
# The "dynamic_symmetric repro" tests alone are checked against the regressions. This is so that
# the pipelines might separate errors that are internally inconsistent.
check-stats () {
  if [ "$#" -ne 2 ]; then echo "check-stats needs 2 arguments" ; exit 911 ; fi

  if [[ "$2" =~ "S" ]]; then # Symmetric
    compare-stats $CONFIGS_DIR/results/$1-all-dynamic_symmetric-repro-def-stats regressions
  fi
  if [[ "$2" =~ "N" ]]; then # Non-symmetric
    compare-stats $CONFIGS_DIR/results/$1-all-dynamic_nonsymmetric-repro-def-stats $CONFIGS_DIR/results/$1-all-dynamic_symmetric-repro-def-stats
  fi
  if [[ "$2" =~ "L" ]]; then # Layout
    compare-stats $CONFIGS_DIR/results/$1-all-dynamic_symmetric-repro-alt-stats $CONFIGS_DIR/results/$1-all-dynamic_symmetric-repro-def-stats
  fi
  if [[ "$2" =~ "D" ]]; then # Debug
    compare-stats $CONFIGS_DIR/results/$1-ocean_only-dynamic_symmetric-debug-def-stats $CONFIGS_DIR/results/$1-all-dynamic_symmetric-repro-def-stats
  fi
  if [[ "$2" =~ "T" ]]; then # sTatic
    compare-stats $CONFIGS_DIR/results/$1-static_ocean_only-static-repro-def-stats $CONFIGS_DIR/results/$1-all-dynamic_symmetric-repro-def-stats
  fi
  if [[ "$2" =~ "R" ]]; then # Restarts
    pushd $JOB_DIR/$STATS_REPO_DIR/tmp-$CONFIGS_DIR-$1-rst > /dev/null
    make -f tools/MRS/Makefile.restart restart_$1_ocean_only restart_$1_ice_ocean_SIS2 -s -k
    popd > /dev/null
  fi

}

# Helper function to compare two files
# Args:
#   $1 is source directory
#   $2 is target directory
#   $3- are file names that should exist relative to both $1 and $2
#
# Operations for `compare-files src/ tgt/ file1 file2 file3`:
#   1. create the md5sum of file1, file2, and file3, in src/ and then run `md5sum-c` in tgt/
#   2. if differences are detected,
#      a. report the "OK" results first, then the "FAILED".
#      b. report the "FAILED".
#      c. for each failed file, show the `diff src/$f tgt/$f`
#   3. if no differences are detected, show `md5sum -c` output so the log lists all files that were checked
compare-files () {
  SRC=$1
  TGT=$2
  shift; shift
  FILES=$@
  ( cd $SRC ; md5sum $FILES ) | ( cd $TGT ; md5sum -c ) | sed -r "s/([A-Za-z0-9_\.\/\-]*): ([A-Z]*)/\2 \1/;s/OK /${GRN}PASS$OFF /;s/FAILED /${RED}FAILED$OFF /"
  FAIL=${PIPESTATUS[1]}
  if [ ! $FAIL == 0 ]; then
    echo Differences follow:
    # All is not well so re-order md5sum to summarize status
    DFILES=$( ( cd $SRC ; md5sum $FILES ) | ( cd $TGT ; md5sum -c 2> /dev/null ) | grep ": FAILED" | sed 's/:.*//')
    for f in $DFILES; do
      echo diff $SRC/$f $TGT/$f | sed "s:$JOB_DIR/$STATS_REPO_DIR/::g;s:$CONFIGS_DIR/results/::"
      diff $SRC/$f $TGT/$f || true
    done
    echo Files $DFILES had differences
    exit 911
  fi
}

# Test the value of param files. All files generated in results/ are looked for $CONFIGS_DIR
# Args:
#   $1 is compiler (gnu, intel, pgi, ...)
check-params () {
  if [ "$#" -ne 1 ]; then echo "check-params needs 1 argument" ; exit 911 ; fi
  section-start-open check-params-$1 "Checking params for $1"
  SRC=$JOB_DIR/$STATS_REPO_DIR/$CONFIGS_DIR/results/$1-all-dynamic_symmetric-repro-def-params
  FILES=$( cd $SRC ; find * -name "*parameter_doc*" -type f )
  compare-files $SRC $JOB_DIR/$STATS_REPO_DIR/$CONFIGS_DIR $FILES
  section-end check-params-$1
}

# Test the value of available_diag files. Only those recorded in $CONFIGS_DIR are checked.
# Args:
#   $1 is compiler (gnu, intel, pgi, ...)
check-diags () {
  if [ "$#" -ne 1 ]; then echo "check-diags needs 1 argument" ; exit 911 ; fi
  section-start-open check-diags-$1 "Checking diagnostics for $1"
  # This checks that any file in the results directory is exactly the same as in regressions/
  SRC=$JOB_DIR/$STATS_REPO_DIR/$CONFIGS_DIR/results/$1-all-dynamic_symmetric-repro-def-params
  TGT=$JOB_DIR/$STATS_REPO_DIR/$CONFIGS_DIR
  EXP_GROUPS=`grep / $TGT/manifest.mk | sed 's:/.*::' | uniq`
  #FILES=$( cd $TGT ; find $EXP_GROUPS -name "*available_diags*" -type f )
  # The following option finds the intersection between all available_diags in both $TGT and $SRC because
  # $SRC contains more than are recorded in $TGT but $TGT might have some that we no longer monitor
  FILES=$( comm -12 <(cd $SRC; find $EXP_GROUPS -name '*available_diags*' -type f | sort) <(cd $TGT; find $EXP_GROUPS -name '*available_diags*' -type f | sort) )
  compare-files $SRC $TGT $FILES
  section-end check-diags-$1
}

# Process command line
START_DIR=`pwd`
while [[ $# -gt 0 ]]; do # Loop through arguments
  cd $START_DIR
  arg=$1
  shift
  case "$arg" in
    -n | --norun)
      DRYRUN=1; echo Dry-run enabled; continue ;;
    +n | ++norun)
      DRYRUN=; echo Dry-run disabled; continue ;;
    -x)
      set -x; continue ;;
    +x)
      set +x; continue ;;
    clean-job-dir)
      clean-job-dir; continue ;;
    create-job-dir)
      create-job-dir https://gitlab.gfdl.noaa.gov/ogrp/Gaea-stats-MOM6-examples.git dev/gfdl; continue ;;
    copy-test-space)
      copy-test-space $1; shift; continue ;;
    mrs-compile)
      mrs-compile $1; shift; continue ;;
    nolibs-ocean-only-compile)
      nolibs-ocean-only-compile $1; shift; continue ;;
    nolibs-ocean-ice-compile)
      nolibs-ocean-ice-compile $1; shift; continue ;;
    run-suite)
      run-suite $1 $2; shift; shift; continue ;;
    check-stats)
      check-stats $1 $2; shift; shift; continue ;;
    check-params)
      check-params $1; shift; continue ;;
    check-diags)
      check-diags $1; shift; continue ;;
    *)
      echo \"$arg\" is not a recognized argument! ; exit 9 ;;
  esac
done
