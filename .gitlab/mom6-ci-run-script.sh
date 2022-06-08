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
rm -f .CI-BATCH-SUCCESS

set -e
set -v

# Run symmetric gnu regressions
section_start gnu_all_sym "Running symmetric gnu"
time make -f tools/MRS/Makefile.run gnu_all -s -j
tar cf gnu_all_sym.tar `find [oicl]* -name "*.stats.*[a-z][a-z][a-z]"`
tar cf gnu_params.tar `find [oicl]* -name "*_parameter_doc.*"`
check_for_core_files
section_end

# Run non-symmetric gnu regressions
section_start gnu_all_nonsym "Running nonsymmetric gnu"
time make -f tools/MRS/Makefile.run ocean_only/circle_obcs/ocean.stats.gnu -s # work around
time make -f tools/MRS/Makefile.run gnu_all -s -j MEMORY=dynamic_nonsymmetric
tar cf gnu_all_nonsym.tar `find [oicl]* -name "*.stats.*[a-z][a-z][a-z]"`
check_for_core_files
section_end

# Run symmetric gnu regressions with alternate layout
section_start gnu_all_layout "Running symmetric gnu with alternate layouts"
time make -f tools/MRS/Makefile.run gnu_all -s -j LAYOUT=alt
tar cf gnu_all_layout.tar `find [oicl]* -name "*.stats.*[a-z][a-z][a-z]"`
check_for_core_files
section_end

# Run symmetric gnu regressions with debug executable
section_start gnu_ocean_only_debug "Running symmetric gnu_ocean_only with debug executable"
time make -f tools/MRS/Makefile.run gnu_ocean_only -s -j MODE=debug
tar cf gnu_ocean_only_debug.tar `find [oicl]* -name "*.stats.*[a-z][a-z][a-z]"`
check_for_core_files
section_end

# Run symmetric static gnu regressions
section_start gnu_all_static "Running symmetric gnu with static executable"
time make -f tools/MRS/Makefile.run gnu_static_ocean_only MEMORY=static -s -j
tar cf gnu_all_static.tar `find [oicl]* -name "*.stats.*[a-z][a-z][a-z]"`
check_for_core_files
section_end

section_start gnu_restarts "Running symmetric gnu restart tests"
time make -f tools/MRS/Makefile.restart gnu_ocean_only -s -j RESTART_STAGE=01
time make -f tools/MRS/Makefile.restart gnu_ice_ocean_SIS2 -s -j RESTART_STAGE=01
time make -f tools/MRS/Makefile.restart gnu_ocean_only -s -j RESTART_STAGE=02
time make -f tools/MRS/Makefile.restart gnu_ice_ocean_SIS2 -s -j RESTART_STAGE=02
time make -f tools/MRS/Makefile.restart gnu_ocean_only -s -j RESTART_STAGE=12
time make -f tools/MRS/Makefile.restart gnu_ice_ocean_SIS2 -s -j RESTART_STAGE=12
tar cf gnu_restarts.tar `find [oilc]*/ -path "*/??.ignore/*" -name "ocean.stats.*[a-z][a-z][a-z]"`
check_for_core_files
find [oilc]* -name "*.ignore" -type d -prune -exec rm -rf {} \;
section_end

# Run symmetric intel regressions
section_start intel_all_sym "Running symmetric intel"
time make -f tools/MRS/Makefile.run intel_all -s -j
tar cf intel_all_sym.tar `find [a-z]* -name "*.stats.*[a-z][a-z][a-z]"`
tar cf intel_params.tar `find [a-z]* -name "*_parameter_doc.*"`
check_for_core_files
section_end

# Run non-symmetric intel regressions
section_start intel_all_nonsym "Running nonsymmetric intel"
time make -f tools/MRS/Makefile.run ocean_only/circle_obcs/ocean.stats.intel -s # work around
time make -f tools/MRS/Makefile.run intel_all -s -j MEMORY=dynamic_nonsymmetric
tar cf intel_all_nonsym.tar `find [oicl]* -name "*.stats.*[a-z][a-z][a-z]"`
check_for_core_files
section_end

# Run symmetric intel regressions with alternate layout
section_start intel_all_layout "Running symmetric intel with alternate layouts"
time make -f tools/MRS/Makefile.run intel_all -s -j LAYOUT=alt
tar cf intel_all_layout.tar `find [oicl]* -name "*.stats.*[a-z][a-z][a-z]"`
check_for_core_files
section_end

# Run symmetric pgi regressions
section_start pgi_all_sym "Running symmetric pgi"
time make -f tools/MRS/Makefile.run pgi_all -s -j
tar cf pgi_all_sym.tar `find [a-z]* -name "*.stats.*[a-z][a-z][a-z]"`
tar cf pgi_params.tar `find [a-z]* -name "*_parameter_doc.*"`
check_for_core_files
section_end

# Run non-symmetric pgi regressions
section_start pgi_all_nonsym "Running nonsymmetric pgi"
time make -f tools/MRS/Makefile.run ocean_only/circle_obcs/ocean.stats.pgi -s # work around
time make -f tools/MRS/Makefile.run pgi_all -s -j MEMORY=dynamic_nonsymmetric
tar cf pgi_all_nonsym.tar `find [oicl]* -name "*.stats.*[a-z][a-z][a-z]"`
check_for_core_files
section_end

# Run symmetric pgi regressions with alternate layout
section_start pgi_all_layout "Running symmetric pgi with alternate layouts"
time make -f tools/MRS/Makefile.run gnu_all -s -j LAYOUT=alt
tar cf pgi_all_layout.tar `find [oicl]* -name "*.stats.*[a-z][a-z][a-z]"`
check_for_core_files
section_end

# Indicate all went well
touch .CI-BATCH-SUCCESS
