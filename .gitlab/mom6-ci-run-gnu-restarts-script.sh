#!/bin/bash

sect=none
clean_stats () { # fn to clean up stats files
  find [oicl]* -name "*.stats.*[a-z][a-z][a-z]" -delete
}
section_start () { # fn to print fold-able banner in CI
  echo -e "\e[0Ksection_start:`date +%s`:$1[collapsed=true]\r\e[0K$2"
  sect=$1
}
section_end () { # fn to close fold-able banner in CI and clean up stats
  echo -e "\e[0Ksection_end:`date +%s`:$sect\r\e[0K"
  clean_stats
}
check_for_core_files () {
  EXIT_CODE=0
  find [oilc]* -name core | grep . && EXIT_CODE=1
  if [[ $EXIT_CODE -gt 0 ]]
  then
    echo "Error: core files found!"
    exit 911
  fi
}

# Make sure we have a clean start
clean_stats
find [oilc]* -name core -delete
rm -f .CI-GNU-RESTARTS-BATCH-SUCCESS

set -e
set -v

# Run symmetric gnu restart tests
section_start gnu_restarts "Running symmetric gnu restart tests"
time make -f tools/MRS/Makefile.restart gnu_ocean_only -s -j RESTART_STAGE=01
time make -f tools/MRS/Makefile.restart gnu_ice_ocean_SIS2 -s -j RESTART_STAGE=01
time make -f tools/MRS/Makefile.restart gnu_ocean_only -s -j RESTART_STAGE=02
time make -f tools/MRS/Makefile.restart gnu_ice_ocean_SIS2 -s -j RESTART_STAGE=02
time make -f tools/MRS/Makefile.restart gnu_ocean_only -s -j RESTART_STAGE=12
time make -f tools/MRS/Makefile.restart gnu_ice_ocean_SIS2 -s -j RESTART_STAGE=12
tar cf - `find [oilc]*/ -path "*/??.ignore/*" -name "ocean.stats.*[a-z][a-z][a-z]"` | tar --one-top-level=results/gnu_restarts -xf -
check_for_core_files
find [oilc]* -name "*.ignore" -type d -prune -exec rm -rf {} \;
section_end

# Indicate all went well
touch .CI-GNU-RESTARTS-BATCH-SUCCESS
