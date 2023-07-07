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
rm -f .CI-GNU-BATCH-SUCCESS

set -e
set -v

# Run symmetric gnu regressions
section_start gnu_all_sym "Running symmetric gnu"
time make -f tools/MRS/Makefile.run gnu_all -s -j
tar cf - `find [oicl]* -name "*.stats.*[a-z][a-z][a-z]"` | tar --one-top-level=results/gnu_all_sym -xf -
tar cf - `find [oicl]* -name "*_parameter_doc.*" -o -name "*available_diags*"` | tar --one-top-level=results/gnu_params -xf -
check_for_core_files
section_end

# Run non-symmetric gnu regressions
section_start gnu_all_nonsym "Running nonsymmetric gnu"
time make -f tools/MRS/Makefile.run ocean_only/circle_obcs/ocean.stats.gnu
time make -f tools/MRS/Makefile.run gnu_all -s -j MEMORY=dynamic_nonsymmetric
tar cf - `find [oicl]* -name "*.stats.*[a-z][a-z][a-z]"` | tar --one-top-level=results/gnu_all_nonsym -xf -
check_for_core_files
section_end

# Run symmetric gnu regressions with alternate layout
section_start gnu_all_layout "Running symmetric gnu with alternate layouts"
time make -f tools/MRS/Makefile.run gnu_all -s -j LAYOUT=alt
tar cf - `find [oicl]* -name "*.stats.*[a-z][a-z][a-z]"` | tar --one-top-level=results/gnu_all_layout -xf -
check_for_core_files
section_end

# Run symmetric gnu regressions with debug executable
section_start gnu_ocean_only_debug "Running symmetric gnu_ocean_only with debug executable"
time make -f tools/MRS/Makefile.run gnu_ocean_only -s -j MODE=debug
tar cf - `find [oicl]* -name "*.stats.*[a-z][a-z][a-z]"` | tar --one-top-level=results/gnu_ocean_only_debug -xf -
check_for_core_files
section_end

# Run symmetric static gnu regressions
section_start gnu_all_static "Running symmetric gnu with static executable"
time make -f tools/MRS/Makefile.run gnu_static_ocean_only MEMORY=static -s -j
tar cf - `find [oicl]* -name "*.stats.*[a-z][a-z][a-z]"` | tar --one-top-level=results/gnu_all_static -xf -
check_for_core_files
section_end

# Indicate all went well
touch .CI-GNU-BATCH-SUCCESS
