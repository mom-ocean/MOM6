# MOM6 Test suite Makefile
#
# Usage:
#   make -j
#       Build the FMS library and test executables
#
#   make -j test
#       Run the test suite, defined in the `tc` directores.
#
#   make clean
#       Delete the MOM6 test executables and dependency builds (FMS)
#
#   make clean.build
#       Delete only the MOM6 test executables
#
#
# Configuration:
#   These settings can be provided as either command-line flags, or saved in a
#   `config.mk` file.
#
# General test configuration:
#   MPIRUN                  MPI job launcher (mpirun, srun, etc)
#   DO_REPRO_TESTS          Enable production ("repro") testing equivalence
#   DO_REGRESSION_TESTS     Enable regression tests (usually dev/gfdl)
#   DO_COVERAGE             Enable code coverage and generate .gcov reports
#   DO_PROFILE              Enable performance profiler comparison tests
#   REQUIRE_CODECOV_UPLOAD  Abort as error if upload to codecov.io fails.
#
# Compiler configuration:
#   CC      C compiler
#   MPICC   MPI C compiler
#   FC      Fortran compiler
#   MPIFC   MPI Fortran compiler
# (NOTE: These are environment variables and may be inherited from a shell.)
#
# Build configuration:
#   FCFLAGS_DEBUG       Testing ("debug") compiler flags
#   FCFLAGS_REPRO       Production ("repro") compiler flags
#   FCFLAGS_OPT         Aggressive optimization compiler flags
#   FCFLAGS_INIT        Variable initialization flags
#   FCFLAGS_COVERAGE    Code coverage flags
#   FCFLAGS_FMS         FMS build flags (default: FCFLAGS_DEBUG)
#
#   LDFLAGS_COVERAGE    Linker coverage flags
#   LDFLAGS_USER        User-defined linker flags (used for all MOM/FMS builds)
#
# Experiment Configuration:
#   EXECS   Executables to be built by `make` or `make all`
#   CONFIGS Model configurations to test (default: `tc*`)
#   DIMS    Dimensional scaling tests
#   TESTS   Tests to run
#
# Regression repository ("target") configuration:
#   MOM_TARGET_SLUG             URL slug (minus domain) of the target repo
#   MOM_TARGET_URL              Full URL of the target repo
#   MOM_TARGET_LOCAL_BRANCH     Target branch name
# (NOTE: These would typically be configured by a CI.)
#
# Output paths:
#   BUILD   Compiled executables and libraries
#   DEPS    Compiled dependencies
#   WORK    Test model output

# TODO: POSIX shell compatibility
SHELL = bash

# No implicit rules, suffixes, or variables
MAKEFLAGS += -rR

# Determine the MOM6 autoconf srcdir
AC_SRCDIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))../ac

# User-defined configuration
-include config.mk

# Set the FMS library
FMS_COMMIT ?= 2023.03
FMS_URL ?= https://github.com/NOAA-GFDL/FMS.git
export FMS_COMMIT
export FMS_URL

# Set the MPI launcher here
# TODO: This needs more automated configuration
MPIRUN ?= mpirun

# Generic compiler variables are passed through to the builds
export CC
export MPICC
export FC
export MPIFC

# Builds are distinguished by FCFLAGS
FCFLAGS ?= -g -O0

FCFLAGS_DEBUG ?= $(FCFLAGS)
FCFLAGS_REPRO ?= -g -O2
FCFLAGS_OPT ?= -g -O3 -mavx -fno-omit-frame-pointer
FCFLAGS_INIT ?=
FCFLAGS_COVERAGE ?= -g -O0 -fbacktrace --coverage
FCFLAGS_FMS ?= $(FCFLAGS)
# Additional notes:
# - These default values are simple, minimalist flags, supported by nearly all
#   compilers, and are somewhat analogous to GFDL's DEBUG and REPRO builds.
#
# - These flags can be configured outside of the Makefile, either with
#   config.mk or as environment variables.

LDFLAGS_COVERAGE ?= --coverage
LDFLAGS_USER ?=

# Set to verify identical DEBUG and REPRO results
DO_REPRO_TESTS ?=

# Enable profiling
DO_PROFILE ?=

# Enable code coverage runs
DO_COVERAGE ?=

# Enable unit tests
DO_UNIT_TESTS ?=

# Check for regressions with target branch
DO_REGRESSION_TESTS ?=

# Report failure if coverage report is not uploaded
REQUIRE_COVERAGE_UPLOAD ?=

# Print logs if an error is encountered
REPORT_ERROR_LOGS ?=

# Time measurement (configurable by the CI)
TIME ?= time

# Legacy external work directory
#WORKSPACE ?=
WORKSPACE ?= .

# Set directories for build/ and work/
BUILD ?= $(WORKSPACE)/build
DEPS ?= $(BUILD)/deps
WORK ?= $(WORKSPACE)/work

# Experiment configuration
EXECS ?= symmetric/MOM6 asymmetric/MOM6 openmp/MOM6
CONFIGS ?= $(wildcard tc*)
DIMS ?= t l h z q r
TESTS ?= grid layout rotate restart openmp nan $(foreach d,$(DIMS),dim.$(d))

# Unit test executables
UNIT_EXECS ?= \
  $(basename $(notdir $(wildcard ../config_src/drivers/unit_tests/*.F90)))

# Timing test executables
TIMING_EXECS ?= \
  $(basename $(notdir $(wildcard ../config_src/drivers/timing_tests/*.F90)))


#---
# Test configuration

# Set if either DO_COVERAGE or DO_UNIT_TESTS is set
run_unit_tests =

# REPRO and DEBUG equivalence
ifdef DO_REPRO_TESTS
  EXECS += repro/MOM6
  TESTS += repro
endif

# Profiling
ifdef DO_PROFILE
  EXECS += opt/MOM6 opt_target/MOM6
endif

# Coverage
ifdef DO_COVERAGE
  EXECS += cov/MOM6
  run_unit_execs = yes
endif

# Unit test executables
ifdef DO_UNIT_TESTS
  run_unit_tests = yes
endif

# If either coverage or unit tests are enabled, build the unit test execs
ifdef run_unit_tests
  EXECS += $(foreach e, $(UNIT_EXECS), unit/$(e))
endif

# Regression testing
ifdef DO_REGRESSION_TESTS
  EXECS += target/MOM6
  TESTS += regression

  MOM_TARGET_SLUG ?= NOAA-GFDL/MOM6
  MOM_TARGET_URL ?= https://github.com/$(MOM_TARGET_SLUG)

  MOM_TARGET_LOCAL_BRANCH ?= dev/gfdl
  MOM_TARGET_BRANCH := origin/$(MOM_TARGET_LOCAL_BRANCH)

  TARGET_CODEBASE = $(BUILD)/target_codebase
else
  MOM_TARGET_URL =
  MOM_TARGET_BRANCH =
  TARGET_CODEBASE =
endif


## Rules

.PHONY: all build.regressions build.prof
all: $(foreach b,$(EXECS),$(BUILD)/$(b))
build.regressions: $(foreach b,symmetric target,$(BUILD)/$(b)/MOM6)
build.prof: $(foreach b,opt opt_target,$(BUILD)/$(b)/MOM6)

# Executable
.PRECIOUS: $(foreach b,$(EXECS),$(BUILD)/$(b))


# Compiler flags

# .testing dependencies
FCFLAGS_DEPS = -I$(abspath $(DEPS)/include)
LDFLAGS_DEPS = -L$(abspath $(DEPS)/lib)
PATH_DEPS = PATH="${PATH}:$(abspath $(DEPS)/bin)"


# Define the build targets in terms of the traditional DEBUG/REPRO/etc labels
SYMMETRIC_FCFLAGS := FCFLAGS="$(FCFLAGS_DEBUG) $(FCFLAGS_INIT) $(FCFLAGS_DEPS)"
ASYMMETRIC_FCFLAGS := FCFLAGS="$(FCFLAGS_DEBUG) $(FCFLAGS_INIT) $(FCFLAGS_DEPS)"
REPRO_FCFLAGS := FCFLAGS="$(FCFLAGS_REPRO) $(FCFLAGS_DEPS)"
OPT_FCFLAGS := FCFLAGS="$(FCFLAGS_OPT) $(FCFLAGS_DEPS)"
OPENMP_FCFLAGS := FCFLAGS="$(FCFLAGS_DEBUG) $(FCFLAGS_INIT) $(FCFLAGS_DEPS)"
TARGET_FCFLAGS := FCFLAGS="$(FCFLAGS_DEBUG) $(FCFLAGS_INIT) $(FCFLAGS_DEPS)"
COV_FCFLAGS := FCFLAGS="$(FCFLAGS_COVERAGE) $(FCFLAGS_DEPS)"

MOM_LDFLAGS := LDFLAGS="$(LDFLAGS_DEPS) $(LDFLAGS_USER)"
COV_LDFLAGS := LDFLAGS="$(LDFLAGS_COVERAGE) $(LDFLAGS_DEPS) $(LDFLAGS_USER)"


# Environment variable configuration
MOM_ENV := $(PATH_FMS)
$(BUILD)/symmetric/Makefile: MOM_ENV += $(SYMMETRIC_FCFLAGS) $(MOM_LDFLAGS)
$(BUILD)/asymmetric/Makefile: MOM_ENV += $(ASYMMETRIC_FCFLAGS) $(MOM_LDFLAGS) \
  MOM_MEMORY=$(AC_SRCDIR)/../config_src/memory/dynamic_nonsymmetric/MOM_memory.h
$(BUILD)/repro/Makefile: MOM_ENV += $(REPRO_FCFLAGS) $(MOM_LDFLAGS)
$(BUILD)/openmp/Makefile: MOM_ENV += $(OPENMP_FCFLAGS) $(MOM_LDFLAGS)
$(BUILD)/target/Makefile: MOM_ENV += $(TARGET_FCFLAGS) $(MOM_LDFLAGS)
$(BUILD)/opt/Makefile: MOM_ENV += $(OPT_FCFLAGS) $(MOM_LDFLAGS)
$(BUILD)/opt_target/Makefile: MOM_ENV += $(OPT_FCFLAGS) $(MOM_LDFLAGS)
$(BUILD)/coupled/Makefile: MOM_ENV += $(SYMMETRIC_FCFLAGS) $(MOM_LDFLAGS)
$(BUILD)/nuopc/Makefile: MOM_ENV += $(SYMMETRIC_FCFLAGS) $(MOM_LDFLAGS)
$(BUILD)/cov/Makefile: MOM_ENV += $(COV_FCFLAGS) $(COV_LDFLAGS)
$(BUILD)/unit/Makefile: MOM_ENV += $(COV_FCFLAGS) $(COV_LDFLAGS)
$(BUILD)/timing/Makefile: MOM_ENV += $(OPT_FCFLAGS) $(MOM_LDFLAGS)

# Configure script flags
$(BUILD)/openmp/Makefile: MOM_ACFLAGS += --enable-openmp
$(BUILD)/coupled/Makefile: MOM_ACFLAGS += --with-driver=FMS_cap
$(BUILD)/nuopc/Makefile: MOM_ACFLAGS += --with-driver=nuopc_cap
$(BUILD)/unit/Makefile: MOM_ACFLAGS += --with-driver=unit_tests
$(BUILD)/timing/Makefile: MOM_ACFLAGS += --with-driver=timing_tests


# Build executables
.NOTPARALLEL:$(foreach e,$(UNIT_EXECS),$(BUILD)/unit/$(e))
$(BUILD)/unit/test_%: $(BUILD)/unit/Makefile FORCE
	cd $(@D) && $(TIME) $(MAKE) $(@F) -j
$(BUILD)/unit/Makefile: $(foreach e,$(UNIT_EXECS),../config_src/drivers/unit_tests/$(e).F90)

.NOTPARALLEL:$(foreach e,$(TIMING_EXECS),$(BUILD)/timing/$(e))
$(BUILD)/timing/time_%: $(BUILD)/timing/Makefile FORCE
	cd $(@D) && $(TIME) $(MAKE) $(@F) -j
$(BUILD)/timing/Makefile: $(foreach e,$(TIMING_EXECS),../config_src/drivers/timing_tests/$(e).F90)

$(BUILD)/%/MOM6: $(BUILD)/%/Makefile FORCE
	cd $(@D) && $(TIME) $(MAKE) $(@F) -j

# Target codebase should use its own build system
$(BUILD)/target/MOM6: $(BUILD)/target FORCE | $(TARGET_CODEBASE)
	$(MAKE) -C $(TARGET_CODEBASE)/.testing BUILD=build build/symmetric/MOM6

$(BUILD)/target: | $(TARGET_CODEBASE)
	ln -s $(abspath $(TARGET_CODEBASE))/.testing/build/symmetric $@

FORCE:


## Use autoconf to construct the Makefile for each target
# TODO: This could all be moved to a top-level MOM6 Makefile
.PRECIOUS: $(BUILD)/%/Makefile
.PRECIOUS: $(BUILD)/%/Makefile.in
.PRECIOUS: $(BUILD)/%/configure
.PRECIOUS: $(BUILD)/%/config.status
.PRECIOUS: $(BUILD)/%/configure.ac
.PRECIOUS: $(BUILD)/%/m4/

$(BUILD)/%/Makefile: $(BUILD)/%/Makefile.in $(BUILD)/%/config.status
	cd $(@D) && ./config.status

$(BUILD)/%/config.status: $(BUILD)/%/configure $(DEPS)/lib/libFMS.a
	cd $(@D) && $(MOM_ENV) ./configure -n --srcdir=$(AC_SRCDIR) $(MOM_ACFLAGS) \
	 || (cat config.log && false)

$(BUILD)/%/Makefile.in: ../ac/Makefile.in | $(BUILD)/%/
	cp ../ac/Makefile.in $(@D)

$(BUILD)/%/configure: $(BUILD)/%/configure.ac $(BUILD)/%/m4/
	autoreconf -if $(@D)

$(BUILD)/%/configure.ac: ../ac/configure.ac | $(BUILD)/%/
	cp ../ac/configure.ac $(@D)

$(BUILD)/%/m4/: ../ac/m4/ | $(BUILD)/%/
	cp -r ../ac/m4 $(@D)

ALL_EXECS = symmetric asymmetric repro openmp opt opt_target coupled nuopc \
  cov unit timing
$(foreach b,$(ALL_EXECS),$(BUILD)/$(b)/):
	mkdir -p $@

# Fetch the regression target codebase
$(TARGET_CODEBASE):
	git clone --recursive $(MOM_TARGET_URL) $@
	cd $@ && git checkout --recurse-submodules $(MOM_TARGET_BRANCH)


## FMS

# Set up the FMS build environment variables
FMS_ENV = \
  PATH="${PATH}:$(realpath ../ac)" \
  FCFLAGS="$(FCFLAGS_FMS)" \
  REPORT_ERROR_LOGS="$(REPORT_ERROR_LOGS)"

$(DEPS)/lib/libFMS.a: $(DEPS)/Makefile $(DEPS)/Makefile.fms.in $(DEPS)/configure.fms.ac $(DEPS)/m4
	$(FMS_ENV) $(MAKE) -C $(DEPS) lib/libFMS.a

$(DEPS)/Makefile: ../ac/deps/Makefile | $(DEPS)
	cp ../ac/deps/Makefile $(DEPS)/Makefile

$(DEPS)/Makefile.fms.in: ../ac/deps/Makefile.fms.in | $(DEPS)
	cp ../ac/deps/Makefile.fms.in $(DEPS)/Makefile.fms.in

$(DEPS)/configure.fms.ac: ../ac/deps/configure.fms.ac | $(DEPS)
	cp ../ac/deps/configure.fms.ac $(DEPS)/configure.fms.ac

$(DEPS)/m4: ../ac/deps/m4 | $(DEPS)
	cp -r ../ac/deps/m4 $(DEPS)/

$(DEPS):
	mkdir -p $(DEPS)

#---
# Verify that the coupled model drivers can be compiled.  This does not verify
# that they can be run, since it would require external submodels.

# NUOPC driver
$(BUILD)/nuopc/mom_ocean_model_nuopc.o: $(BUILD)/nuopc/Makefile
	cd $(@D) && make $(@F)
check_mom6_api_nuopc: $(BUILD)/nuopc/mom_ocean_model_nuopc.o

# GFDL coupled driver
$(BUILD)/coupled/ocean_model_MOM.o: $(BUILD)/coupled/Makefile
	cd $(@D) && make $(@F)
check_mom6_api_coupled: $(BUILD)/coupled/ocean_model_MOM.o


## Testing

.PHONY: test
test: $(foreach t,$(TESTS),test.$(t))

# NOTE: We remove tc3 (OBC) from grid test since it cannot run asymmetric grids

# NOTE: rotation diag chksum disabled since we cannot yet compare rotationally
#       equivalent diagnostics

# TODO: restart checksum comparison is not yet implemented

.PHONY: $(foreach t,$(TESTS),test.$(t))
test.grid: $(foreach c,$(filter-out tc3,$(CONFIGS)),$(c).grid $(c).grid.diag)
test.layout: $(foreach c,$(CONFIGS),$(c).layout $(c).layout.diag)
test.rotate: $(foreach c,$(CONFIGS),$(c).rotate)
test.restart: $(foreach c,$(CONFIGS),$(c).restart)
test.repro: $(foreach c,$(CONFIGS),$(c).repro $(c).repro.diag)
test.openmp: $(foreach c,$(CONFIGS),$(c).openmp $(c).openmp.diag)
test.nan: $(foreach c,$(CONFIGS),$(c).nan $(c).nan.diag)
test.regression: $(foreach c,$(CONFIGS),$(c).regression $(c).regression.diag)
test.dim: $(foreach d,$(DIMS),test.dim.$(d))
define TEST_DIM_RULE
test.dim.$(1): $(foreach c,$(CONFIGS),$(c).dim.$(1) $(c).dim.$(1).diag)
endef
$(foreach d,$(DIMS),$(eval $(call TEST_DIM_RULE,$(d))))

.PHONY: run.symmetric run.asymmetric run.nans run.openmp run.cov
run.symmetric: $(foreach c,$(CONFIGS),$(WORK)/$(c)/symmetric/ocean.stats)
run.asymmetric: $(foreach c,$(filter-out tc3,$(CONFIGS)),$(CONFIGS),$(WORK)/$(c)/asymmetric/ocean.stats)
run.nan: $(foreach c,$(CONFIGS),$(WORK)/$(c)/nan/ocean.stats)
run.openmp: $(foreach c,$(CONFIGS),$(WORK)/$(c)/openmp/ocean.stats)
run.cov: $(foreach c,$(CONFIGS),$(WORK)/$(c)/cov/ocean.stats)

# Configuration test rules
# $(1): Configuration name (tc1, tc2, &c.)
# $(2): Excluded tests
.PRECIOUS: $(foreach c,$(CONFIGS),$(c))
define CONFIG_RULE
$(1): \
  $(foreach t,$(filter-out $(2),$(TESTS)),$(1).$(t)) \
  $(foreach t,$(filter-out $(2) rotate restart,$(TESTS)),$(1).$(t).diag)
endef
$(foreach c,$(filter-out tc3,$(CONFIGS)),$(eval $(call CONFIG_RULE,$(c),)))
# NOTE: tc3 uses OBCs and does not support asymmetric grid
$(eval $(call CONFIG_RULE,tc3,grid))

# Color highlights for test results
RED = \033[0;31m
GREEN = \033[0;32m
YELLOW = \033[0;33m
MAGENTA = \033[0;35m
RESET = \033[0m

DONE = ${GREEN}DONE${RESET}
PASS = ${GREEN}PASS${RESET}
WARN = ${YELLOW}WARN${RESET}
FAIL = ${RED}FAIL${RESET}

# Comparison rules
# $(1): Configuration (tc1, tc2, &c.)
# $(2): Test type (grid, layout, &c.)
# $(3): Comparison targets (symmetric asymmetric, symmetric layout, &c.)
define CMP_RULE
.PRECIOUS: $(foreach b,$(3),$(WORK)/$(1)/$(b)/ocean.stats)
$(1).$(2): $(foreach b,$(3),$(WORK)/$(1)/$(b)/ocean.stats)
	@test "$$(shell ls -A $(WORK)/results/$(1) 2>/dev/null)" || rm -rf $(WORK)/results/$(1)
	@cmp $$^ || !( \
	  mkdir -p $(WORK)/results/$(1); \
	  (diff $$^ | tee $(WORK)/results/$(1)/ocean.stats.$(2).diff | head -n 20) ; \
	  echo -e "$(FAIL): Solutions $(1).$(2) have changed." \
	)
	@echo -e "$(PASS): Solutions $(1).$(2) agree."

.PRECIOUS: $(foreach b,$(3),$(WORK)/$(1)/$(b)/chksum_diag)
$(1).$(2).diag: $(foreach b,$(3),$(WORK)/$(1)/$(b)/chksum_diag)
	@cmp $$^ || !( \
	  mkdir -p $(WORK)/results/$(1); \
	  (diff $$^ | tee $(WORK)/results/$(1)/chksum_diag.$(2).diff | head -n 20) ; \
	  echo -e "$(FAIL): Diagnostics $(1).$(2).diag have changed." \
	)
	@echo -e "$(PASS): Diagnostics $(1).$(2).diag agree."
endef

$(foreach c,$(CONFIGS),$(eval $(call CMP_RULE,$(c),grid,symmetric asymmetric)))
$(foreach c,$(CONFIGS),$(eval $(call CMP_RULE,$(c),layout,symmetric layout)))
$(foreach c,$(CONFIGS),$(eval $(call CMP_RULE,$(c),rotate,symmetric rotate)))
$(foreach c,$(CONFIGS),$(eval $(call CMP_RULE,$(c),repro,symmetric repro)))
$(foreach c,$(CONFIGS),$(eval $(call CMP_RULE,$(c),openmp,symmetric openmp)))
$(foreach c,$(CONFIGS),$(eval $(call CMP_RULE,$(c),nan,symmetric nan)))
define CONFIG_DIM_RULE
$(1).dim: $(foreach d,$(DIMS),$(1).dim.$(d))
$(foreach d,$(DIMS),$(eval $(call CMP_RULE,$(1),dim.$(d),symmetric dim.$(d))))
endef
$(foreach c,$(CONFIGS),$(eval $(call CONFIG_DIM_RULE,$(c))))


# Custom comparison rules

# Restart tests only compare the final stat record
.PRECIOUS: $(foreach b,symmetric restart target,$(WORK)/%/$(b)/ocean.stats)
%.restart: $(foreach b,symmetric restart,$(WORK)/%/$(b)/ocean.stats)
	@test "$(shell ls -A $(WORK)/results/$* 2>/dev/null)" || rm -rf $(WORK)/results/$*
	@cmp $(foreach f,$^,<(tr -s ' ' < $(f) | cut -d ' ' -f3- | tail -n 1)) \
	  || !( \
	    mkdir -p $(WORK)/results/$*; \
	    (diff $^ | tee $(WORK)/results/$*/chksum_diag.restart.diff | head -n 20) ; \
	    echo -e "$(FAIL): Solutions $*.restart have changed." \
	  )
	@echo -e "$(PASS): Solutions $*.restart agree."

# TODO: chksum_diag parsing of restart files

# stats rule is unchanged, but we cannot use CMP_RULE to generate it.
%.regression: $(foreach b,symmetric target,$(WORK)/%/$(b)/ocean.stats)
	@test "$(shell ls -A $(WORK)/results/$* 2>/dev/null)" || rm -rf $(WORK)/results/$*
	@cmp $^ || !( \
	  mkdir -p $(WORK)/results/$*; \
	  (diff $^ | tee $(WORK)/results/$*/ocean.stats.regression.diff | head -n 20) ; \
	  echo -e "$(FAIL): Solutions $*.regression have changed." \
	)
	@echo -e "$(PASS): Solutions $*.regression agree."

# Regression testing only checks for changes in existing diagnostics
.PRECIOUS: $(WORK)/%/target/chksum_diag
%.regression.diag: $(foreach b,symmetric target,$(WORK)/%/$(b)/chksum_diag)
	@! diff $^ | grep "^[<>]" | grep "^>" > /dev/null \
	  || ! (\
	    mkdir -p $(WORK)/results/$*; \
	    (diff $^ | tee $(WORK)/results/$*/chksum_diag.regression.diff | head -n 20) ; \
	    echo -e "$(FAIL): Diagnostics $*.regression.diag have changed." \
	  )
	@cmp $^ || ( \
	  diff $^ | head -n 20; \
	  echo -e "$(WARN): New diagnostics in $<" \
	)
	@echo -e "$(PASS): Diagnostics $*.regression.diag agree."


#---
# Preprocessing
# NOTE: This only support tc4, but can be generalized over all tests.
.PHONY: preproc
preproc: tc4/Makefile
	cd tc4 && $(MAKE) LAUNCHER="$(MPIRUN)"
preproc-compile: tc4/Makefile
	cd tc4 && $(MAKE) executables

tc4/Makefile: tc4/configure tc4/Makefile.in
	cd $(@D) && ./configure || (cat config.log && false)

tc4/configure: tc4/configure.ac
	cd $(@D) && autoreconf -if


#---
# Test run output files

# Rule to build $(WORK)/<tc>/{ocean.stats,chksum_diag}.<tag>
# $(1): Test configuration name <tag>
# $(2): Executable type
# $(3): Enable coverage flag
# $(4): MOM_override configuration
# $(5): Environment variables
# $(6): Number of MPI ranks

define STAT_RULE
$(WORK)/%/$(1)/ocean.stats $(WORK)/%/$(1)/chksum_diag: $(BUILD)/$(2)/MOM6 | preproc
	@echo "Running test $$*.$(1)..."
	mkdir -p $$(@D)
	cp -RL $$*/* $$(@D)
	echo -e "$(4)" > $$(@D)/MOM_override
	rm -f $(WORK)/results/$$*/std.$(1).{out,err}
	cd $$(@D) \
	  && $(TIME) $(5) $(MPIRUN) -n $(6) $$(abspath $$<) 2> std.err > std.out \
	  || !( \
	    mkdir -p ../../../results/$$*/ ; \
	    cat std.out | tee ../../../results/$$*/std.$(1).out | tail -n 40 ; \
	    cat std.err | tee ../../../results/$$*/std.$(1).err | tail -n 40 ; \
	    rm ocean.stats chksum_diag ; \
	    echo -e "$(FAIL): $$*.$(1) failed at runtime." \
	  )
	@echo -e "$(DONE): $$*.$(1); no runtime errors."
	if [ $(3) ]; then \
	  mkdir -p $(WORK)/results/$$* ; \
	  cd $(BUILD)/$(2) ; \
	  gcov -b *.gcda > gcov.$$*.$(1).out ; \
	  find -name "*.gcov" -exec sed -i -r 's/^( *[0-9]*)\*:/ \1:/g' {} \; ; \
	fi
endef


# Upload coverage reports
CODECOV_UPLOADER_URL ?= https://uploader.codecov.io/latest/linux/codecov
CODECOV_TOKEN ?=

ifdef CODECOV_TOKEN
  CODECOV_TOKEN_ARG = -t $(CODECOV_TOKEN)
else
  CODECOV_TOKEN_ARG =
endif

codecov:
	curl -s $(CODECOV_UPLOADER_URL) -o $@
	chmod +x codecov

.PHONY: report.cov
report.cov: run.cov codecov
	./codecov $(CODECOV_TOKEN_ARG) -R $(BUILD)/cov -Z -f "*.gcov" \
	  > $(BUILD)/cov/codecov.out \
	  2> $(BUILD)/cov/codecov.err \
	  && echo -e "${MAGENTA}Report uploaded to codecov.${RESET}" \
	  || { \
	    cat $(BUILD)/cov/codecov.err ; \
	    echo -e "${RED}Failed to upload report.${RESET}" ; \
	    if [ "$(REQUIRE_COVERAGE_UPLOAD)" = true ] ; then false ; fi ; \
	  }

# Define $(,) as comma escape character
, := ,

$(eval $(call STAT_RULE,symmetric,symmetric,,,,1))
$(eval $(call STAT_RULE,asymmetric,asymmetric,,,,1))
$(eval $(call STAT_RULE,target,target,,,,1))
$(eval $(call STAT_RULE,repro,repro,,,,1))
$(eval $(call STAT_RULE,openmp,openmp,,,GOMP_CPU_AFFINITY=0,1))
$(eval $(call STAT_RULE,layout,symmetric,,LAYOUT=2$(,)1,,2))
$(eval $(call STAT_RULE,rotate,symmetric,,ROTATE_INDEX=True\nINDEX_TURNS=1,,1))
$(eval $(call STAT_RULE,nan,symmetric,,,MALLOC_PERTURB_=1,1))
$(eval $(call STAT_RULE,dim.t,symmetric,,T_RESCALE_POWER=11,,1))
$(eval $(call STAT_RULE,dim.l,symmetric,,L_RESCALE_POWER=11,,1))
$(eval $(call STAT_RULE,dim.h,symmetric,,H_RESCALE_POWER=11,,1))
$(eval $(call STAT_RULE,dim.z,symmetric,,Z_RESCALE_POWER=11,,1))
$(eval $(call STAT_RULE,dim.q,symmetric,,Q_RESCALE_POWER=11,,1))
$(eval $(call STAT_RULE,dim.r,symmetric,,R_RESCALE_POWER=11,,1))

$(eval $(call STAT_RULE,cov,cov,true,,,1))

# Generate the half-period input namelist as follows:
#  1. Fetch DAYMAX and TIMEUNIT from MOM_input
#  2. Convert DAYMAX from TIMEUNIT to seconds
#  3. Apply seconds to `ocean_solo_nml` inside input.nml.
# NOTE: Assumes that runtime set by DAYMAX, will fail if set by input.nml
$(WORK)/%/restart/ocean.stats: $(BUILD)/symmetric/MOM6 | preproc
	rm -rf $(@D)
	mkdir -p $(@D)
	cp -RL $*/* $(@D)
	mkdir -p $(@D)/RESTART
	# Set the half-period
	cd $(@D) \
	  && daymax=$$(grep DAYMAX MOM_input | cut -d '!' -f 1 | cut -d '=' -f 2 | xargs) \
	  && timeunit=$$(grep TIMEUNIT MOM_input | cut -d '!' -f 1 | cut -d '=' -f 2 | xargs) \
	  && if [ -z "$${timeunit}" ]; then timeunit="8.64e4"; fi \
	  && printf -v timeunit_int "%.f" "$${timeunit}" \
	  && halfperiod=$$(awk -v t=$${daymax} -v dt=$${timeunit} 'BEGIN {printf "%.f", 0.5*t*dt}') \
	  && printf "\n&ocean_solo_nml\n    seconds = $${halfperiod}\n/\n" >> input.nml
	# Remove any previous archived output
	rm -f $(WORK)/results/$*/std.restart{1,2}.{out,err}
	# Run the first half-period
	cd $(@D) && $(TIME) $(MPIRUN) -n 1 $(abspath $<) 2> std1.err > std1.out \
	  || !( \
	    cat std1.out | tee ../../../results/$*/std.restart1.out | tail -n 40 ; \
	    cat std1.err | tee ../../../results/$*/std.restart1.err | tail -n 40 ; \
	    echo -e "$(FAIL): $*.restart failed at runtime." \
	  )
	# Setup the next inputs
	cd $(@D) && rm -rf INPUT && mv RESTART INPUT
	mkdir $(@D)/RESTART
	cd $(@D) && sed -i -e "s/input_filename *= *'n'/input_filename = 'r'/g" input.nml
	# Run the second half-period
	cd $(@D) && $(TIME) $(MPIRUN) -n 1 $(abspath $<) 2> std2.err > std2.out \
	  || !( \
	    cat std2.out | tee ../../../results/$*/std.restart2.out | tail -n 40 ; \
	    cat std2.err | tee ../../../results/$*/std.restart2.err | tail -n 40 ; \
	    echo -e "$(FAIL): $*.restart failed at runtime." \
	  )

# TODO: Restart checksum diagnostics

#---
# Not a true rule; only call this after `make test` to summarize test results.
.PHONY: test.summary
test.summary:
	./tools/report_test_results.sh $(WORK)/results


#---
# Unit test

# NOTE: Using file parser gcov report as a proxy for test completion
.PHONY: run.cov.unit
run.cov.unit: $(foreach t,$(UNIT_EXECS),$(BUILD)/unit/$(t).F90.gcov)

.PHONY: build.unit
build.unit: $(foreach f, $(UNIT_EXECS), $(BUILD)/unit/$(f))
.PHONY: run.unit
run.unit: $(foreach f, $(UNIT_EXECS), work/unit/$(f).out)
.PHONY: build.timing
build.timing: $(foreach f, $(TIMING_EXECS), $(BUILD)/timing/$(f))
.PHONY: run.timing
run.timing: $(foreach f, $(TIMING_EXECS), work/timing/$(f).out)
.PHONY: show.timing
show.timing: $(foreach f, $(TIMING_EXECS), work/timing/$(f).show)
$(WORK)/timing/%.show:
	./tools/disp_timing.py $(@:.show=.out)


# Invoke the above unit/timing rules for a "target" code
# Invoke with appropriate macros defines, i.e.
#   make build.timing_target MOM_TARGET_URL=... MOM_TARGET_BRANCH=... TARGET_CODEBASE=$(BUILD)/target_codebase
#   make run.timing_target TARGET_CODEBASE=$(BUILD)/target_codebase

TIMING_TARGET_EXECS ?= $(basename $(notdir $(wildcard $(TARGET_CODEBASE)/config_src/drivers/timing_tests/*.F90) ) )

.PHONY: build.timing_target
build.timing_target: $(foreach f, $(TIMING_TARGET_EXECS), $(TARGET_CODEBASE)/.testing/$(BUILD)/timing/$(f))
.PHONY: run.timing_target
run.timing_target: $(foreach f, $(TIMING_TARGET_EXECS), $(TARGET_CODEBASE)/.testing/work/timing/$(f).out)
.PHONY: compare.timing
compare.timing: $(foreach f, $(filter $(TIMING_EXECS),$(TIMING_TARGET_EXECS)), work/timing/$(f).compare)
$(WORK)/timing/%.compare: $(TARGET_CODEBASE)
	./tools/disp_timing.py -r $(TARGET_CODEBASE)/.testing/$(@:.compare=.out) $(@:.compare=.out)
$(TARGET_CODEBASE)/.testing/%: | $(TARGET_CODEBASE)
	cd $(TARGET_CODEBASE)/.testing && make $*


# General rule to run a unit test executable
# Pattern is to run $(BUILD)/unit/executable and direct output to executable.out
$(WORK)/unit/%.out: $(BUILD)/unit/%
	@mkdir -p $(@D)
	cd $(@D) ; $(TIME) $(MPIRUN) -n 1 $(abspath $<) 2> >(tee $*.err) > $*.out

# The file parser uses a separate rule to support two-core tests.
$(WORK)/unit/test_MOM_file_parser.out: $(BUILD)/unit/test_MOM_file_parser
	if [ $(REPORT_COVERAGE) ]; then \
	  find $(BUILD)/unit -name *.gcda -exec rm -f '{}' \; ; \
	fi
	mkdir -p $(@D)
	cd $(@D) \
	  && rm -f input.nml logfile.0000*.out *_input MOM_parameter_doc.* \
	  && $(TIME) $(MPIRUN) -n 1 $(abspath $<) 2> test_MOM_file_parser.err > test_MOM_file_parser.out \
	  || !( \
	    cat test_MOM_file_parser.out | tail -n 100 ; \
	    cat test_MOM_file_parser.err | tail -n 100 ; \
	  )
	cd $(@D) \
	  && $(TIME) $(MPIRUN) -n 2 $(abspath $<) 2> p2.test_MOM_file_parser.err > p2.test_MOM_file_parser.out \
	  || !( \
	    cat p2.test_MOM_file_parser.out | tail -n 100 ; \
	    cat p2.test_MOM_file_parser.err | tail -n 100 ; \
	  )

$(BUILD)/unit/test_%.F90.gcov: $(WORK)/unit/test_%.out
	cd $(@D) \
	  && gcov -b *.gcda > gcov.unit.out
	find $(@D) -name "*.gcov" -exec sed -i -r 's/^( *[0-9]*)\*:/ \1:/g' {} \;

.PHONY: report.cov.unit
report.cov.unit: $(foreach t,$(UNIT_EXECS),$(BUILD)/unit/$(t).F90.gcov) codecov
	./codecov $(CODECOV_TOKEN_ARG) -R $(BUILD)/unit -f "*.gcov" -Z -n "Unit tests" \
	    > $(BUILD)/unit/codecov.out \
	    2> $(BUILD)/unit/codecov.err \
	  && echo -e "${MAGENTA}Report uploaded to codecov.${RESET}" \
	  || { \
	    cat $(BUILD)/unit/codecov.err ; \
	    echo -e "${RED}Failed to upload report.${RESET}" ; \
	    if [ "$(REQUIRE_COVERAGE_UPLOAD)" = true ] ; then false ; fi ; \
	  }

$(WORK)/timing/%.out: $(BUILD)/timing/% FORCE
	@mkdir -p $(@D)
	@echo Running $< in $(@D)
	@cd $(@D) ; $(TIME) $(MPIRUN) -n 1 $(abspath $<) 2> $*.err > $*.out


## Profiling based on FMS clocks

PCONFIGS = p0

.PHONY: profile
profile: $(foreach p,$(PCONFIGS), prof.$(p))

.PHONY: prof.p0
prof.p0: $(WORK)/p0/opt/clocks.json $(WORK)/p0/opt_target/clocks.json
	python tools/compare_clocks.py $^

$(WORK)/p0/%/clocks.json: $(WORK)/p0/%/std.out
	python tools/parse_fms_clocks.py -d $(@D) $^ > $@ \
	  || !( rm $@ )

$(WORK)/p0/opt/std.out: $(BUILD)/opt/MOM6
$(WORK)/p0/opt_target/std.out: $(BUILD)/opt_target/MOM6

$(WORK)/p0/%/std.out:
	mkdir -p $(@D)
	cp -RL p0/* $(@D)
	mkdir -p $(@D)/RESTART
	echo -e "" > $(@D)/MOM_override
	cd $(@D) \
	  && $(MPIRUN) -n 1 $(abspath $<) 2> std.err > std.out


## Profiling based on perf output

# TODO: This expects the -e flag, can I handle it in the command?
PERF_EVENTS ?=

.PHONY: perf
perf: $(foreach p,$(PCONFIGS), perf.$(p))

.PHONY: prof.p0
perf.p0: $(WORK)/p0/opt/profile.json $(WORK)/p0/opt_target/profile.json
	python tools/compare_perf.py $^

$(WORK)/p0/%/profile.json: $(WORK)/p0/%/perf.data
	python tools/parse_perf.py -f $< > $@

$(WORK)/p0/opt/perf.data: $(BUILD)/opt/MOM6
$(WORK)/p0/opt_target/perf.data: $(BUILD)/opt_target/MOM6

$(WORK)/p0/%/perf.data:
	mkdir -p $(@D)
	cp -RL p0/* $(@D)
	mkdir -p $(@D)/RESTART
	echo -e "" > $(@D)/MOM_override
	cd $(@D) \
	  && perf record \
	    -F 3999 \
	    ${PERF_EVENTS} \
	    ../../../$< 2> std.perf.err > std.perf.out \
	  || cat std.perf.err


## Cleanup
# NOTE: These tests assert that we are in the .testing directory.

.PHONY: clean
clean: clean.build clean.stats
	rm -rf $(BUILD)


.PHONY: clean.build
clean.build:
	@[ $$(basename $$(pwd)) = .testing ]
	for b in $(ALL_EXECS); do \
	  rm -rf $(BUILD)/$${b}; \
	done


.PHONY: clean.stats
clean.stats:
	@[ $$(basename $$(pwd)) = .testing ]
	rm -rf $(WORK)


.PHONY: clean.preproc
clean.preproc:
	@if [ -f tc4/Makefile ] ; then \
	  cd tc4 && make clean ; \
	fi
