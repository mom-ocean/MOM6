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
#   	Delete only the MOM6 test executables
#
#
# Configuration:
#   These settings can be provided as either command-line flags, or saved in a
#   `config.mk` file.
#
# General test configuration:
#   MPIRUN                  MPI job launcher (mpirun, srun, etc)
#   FRAMEWORK               Model framework (fms1 or fms2)
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
#   BUILDS  Executables to be built by `make` or `make all`
#   CONFIGS Model configurations to test (default: `tc*`)
#   TESTS   Tests to run
#   DIMS    Dimensional scaling tests
#
# Regression repository ("target") configuration:
#   MOM_TARGET_SLUG             URL slug (minus domain) of the target repo
#   MOM_TARGET_URL              Full URL of the target repo
#   MOM_TARGET_LOCAL_BRANCH     Target branch name
# (NOTE: These would typically be configured by a CI.)
#
#----

# TODO: POSIX shell compatibility
SHELL = bash

# No implicit rules
.SUFFIXES:

# No implicit variables
MAKEFLAGS += -R

# User-defined configuration
-include config.mk

# Set the infra framework
FRAMEWORK ?= fms1

# Set the MPI launcher here
# TODO: This needs more automated configuration
MPIRUN ?= mpirun

# Generic compiler variables are pass through to the builds
export CC
export MPICC
export FC
export MPIFC

# Builds are distinguished by FCFLAGS
FCFLAGS_DEBUG ?= -g -O0
FCFLAGS_REPRO ?= -g -O2
FCFLAGS_OPT ?= -g -O3 -mavx -fno-omit-frame-pointer
FCFLAGS_INIT ?=
FCFLAGS_COVERAGE ?= -g -O0 -fbacktrace --coverage
FCFLAGS_FMS ?= $(FCFLAGS_DEBUG)
# Additional notes:
# - These default values are simple, minimalist flags, supported by nearly all
#   compilers, and are somewhat analogous to GFDL's DEBUG and REPRO builds.
#
# - These flags can be configured outside of the Makefile, either with
#   config.mk or as environment variables.

LDFLAGS_COVERAGE ?= --coverage
LDFLAGS_USER ?=

# Set to `true` to require identical results from DEBUG and REPRO builds
# NOTE: Many compilers (Intel, GCC on ARM64) do not produce identical results
#   across DEBUG and REPRO builds (as defined below), so we disable on
#   default.
DO_REPRO_TESTS ?=

# Enable profiling
DO_PROFILE ?=

# Enable code coverage runs
DO_COVERAGE ?=

# Report failure if coverage report is not uploaded
REQUIRE_COVERAGE_UPLOAD ?=

# Print logs if an error is encountered
REPORT_ERROR_LOGS ?=

# Time measurement (configurable by the CI)
TIME ?= time


# Experiment configuration
BUILDS ?= symmetric/MOM6 asymmetric/MOM6 openmp/MOM6
CONFIGS ?= $(wildcard tc*)
TESTS ?= grid layout rotate restart openmp nan $(foreach d,$(DIMS),dim.$(d))
DIMS ?= t l h z q r


#---
# Test configuration

# REPRO tests enable reproducibility with optimization, and often do not match
# the DEBUG results in older GCCs and vendor compilers, so we can optionally
# disable them.
ifeq ($(DO_REPRO_TESTS), true)
  BUILDS += repro/MOM6
  TESTS += repro
endif

# Profiling
ifeq ($(DO_PROFILE), true)
  BUILDS += opt/MOM6 opt_target/MOM6
endif

# Unit testing
UNIT_EXECS ?= MOM_unit_tests
ifeq ($(DO_COVERAGE), true)
  BUILDS += cov/MOM6 $(foreach e, $(UNIT_EXECS), unit/$(e))
endif

ifeq ($(DO_PROFILE), false)
  BUILDS += opt/MOM6 opt_target/MOM6
endif


DO_REGRESSION_TESTS ?=
ifeq ($(DO_REGRESSION_TESTS), true)
  BUILDS += target/MOM6
  TESTS += regression

  MOM_TARGET_SLUG ?= NOAA-GFDL/MOM6
  MOM_TARGET_URL ?= https://github.com/$(MOM_TARGET_SLUG)

  MOM_TARGET_LOCAL_BRANCH ?= dev/gfdl
  MOM_TARGET_BRANCH := origin/$(MOM_TARGET_LOCAL_BRANCH)

  TARGET_CODEBASE = build/target_codebase
else
  MOM_TARGET_URL =
  MOM_TARGET_BRANCH =
  TARGET_CODEBASE =
endif

# List of source files to link this Makefile's dependencies to model Makefiles
# Assumes a depth of two, and the following extensions: F90 inc c h
# (1): Root directory
# NOTE: extensions could be a second variable
SOURCE = \
  $(foreach ext,F90 inc c h,$(wildcard $(1)/*/*.$(ext) $(1)/*/*/*.$(ext)))

MOM_SOURCE = \
  $(call SOURCE,../src) \
  $(wildcard ../config_src/drivers/solo_driver/*.F90) \
  $(wildcard ../config_src/ext*/*/*.F90)

TARGET_SOURCE = \
  $(call SOURCE,build/target_codebase/src) \
  $(wildcard build/target_codebase/config_src/drivers/solo_driver/*.F90) \
  $(wildcard build/target_codebase/config_src/ext*/*.F90)

# NOTE: Current default framework is FMS1, but this could change.
ifeq ($(FRAMEWORK), fms2)
  MOM_SOURCE +=$(wildcard ../config_src/infra/FMS2/*.F90)
  TARGET_SOURCE += $(wildcard build/target_codebase/config_src/infra/FMS2/*.F90)
else
  MOM_SOURCE += $(wildcard ../config_src/infra/FMS1/*.F90)
  TARGET_SOURCE += $(wildcard build/target_codebase/config_src/infra/FMS1/*.F90)
endif

FMS_SOURCE = $(call SOURCE,deps/fms/src)


#---
# Python preprocessing environment configuration

HAS_NUMPY = $(shell python -c "import numpy" 2> /dev/null && echo "yes")
HAS_NETCDF4 = $(shell python -c "import netCDF4" 2> /dev/null && echo "yes")

USE_VENV =
ifneq ($(HAS_NUMPY), yes)
  USE_VENV = yes
endif
ifneq ($(HAS_NETCDF4), yes)
  USE_VENV = yes
endif

# When disabled, activation is a null operation (`true`)
VENV_PATH =
VENV_ACTIVATE = true
ifeq ($(USE_VENV), yes)
  VENV_PATH = work/local-env
  VENV_ACTIVATE = . $(VENV_PATH)/bin/activate
endif


#---
# Rules

.PHONY: all build.regressions build.prof
all: $(foreach b,$(BUILDS),build/$(b)) $(VENV_PATH)
build.regressions: $(foreach b,symmetric target,build/$(b)/MOM6)
build.prof: $(foreach b,opt opt_target,build/$(b)/MOM6)

# Executable
.PRECIOUS: $(foreach b,$(BUILDS),build/$(b))


# Compiler flags

# .testing dependencies
# TODO: We should probably build TARGET with the FMS that it was configured
#   to use.  But for now we use the same FMS over all builds.
FCFLAGS_DEPS = -I../../deps/include
LDFLAGS_DEPS = -L../../deps/lib
PATH_DEPS = PATH="${PATH}:../../deps/bin"


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
build/symmetric/Makefile: MOM_ENV=$(PATH_FMS) $(SYMMETRIC_FCFLAGS) $(MOM_LDFLAGS)
build/asymmetric/Makefile: MOM_ENV=$(PATH_FMS) $(ASYMMETRIC_FCFLAGS) $(MOM_LDFLAGS)
build/repro/Makefile: MOM_ENV=$(PATH_FMS) $(REPRO_FCFLAGS) $(MOM_LDFLAGS)
build/openmp/Makefile: MOM_ENV=$(PATH_FMS) $(OPENMP_FCFLAGS) $(MOM_LDFLAGS)
build/target/Makefile: MOM_ENV=$(PATH_FMS) $(TARGET_FCFLAGS) $(MOM_LDFLAGS)
build/opt/Makefile: MOM_ENV=$(PATH_FMS) $(OPT_FCFLAGS) $(MOM_LDFLAGS)
build/opt_target/Makefile: MOM_ENV=$(PATH_FMS) $(OPT_FCFLAGS) $(MOM_LDFLAGS)
build/coupled/Makefile: MOM_ENV=$(PATH_FMS) $(SYMMETRIC_FCFLAGS) $(MOM_LDFLAGS)
build/nuopc/Makefile: MOM_ENV=$(PATH_FMS) $(SYMMETRIC_FCFLAGS) $(MOM_LDFLAGS)
build/mct/Makefile: MOM_ENV=$(PATH_FMS) $(SYMMETRIC_FCFLAGS) $(MOM_LDFLAGS)
build/cov/Makefile: MOM_ENV=$(PATH_FMS) $(COV_FCFLAGS) $(COV_LDFLAGS)
build/unit/Makefile: MOM_ENV=$(PATH_FMS) $(COV_FCFLAGS) $(COV_LDFLAGS)

# Configure script flags
build/symmetric/Makefile: MOM_ACFLAGS=
build/asymmetric/Makefile: MOM_ACFLAGS=--enable-asymmetric
build/repro/Makefile: MOM_ACFLAGS=
build/openmp/Makefile: MOM_ACFLAGS=--enable-openmp
build/target/Makefile: MOM_ACFLAGS=
build/opt/Makefile: MOM_ACFLAGS=
build/opt_target/Makefile: MOM_ACFLAGS=
build/coupled/Makefile: MOM_ACFLAGS=--with-driver=FMS_cap
build/nuopc/Makefile: MOM_ACFLAGS=--with-driver=nuopc_cap
build/mct/Makefile: MOM_ACFLAGS=--with-driver=mct_cap
build/cov/Makefile: MOM_ACFLAGS=
build/unit/Makefile: MOM_ACFLAGS=--with-driver=unit_tests

# Fetch regression target source code
build/target/Makefile: | $(TARGET_CODEBASE)
build/opt_target/Makefile: | $(TARGET_CODEBASE)


# Define source code dependencies
# NOTE: ./configure is too much, but Makefile is not enough!
#   Ideally we only want to re-run both Makefile and mkmf, but the mkmf call
#   is inside ./configure, so we must re-run ./configure as well.
build/target_codebase/configure: $(TARGET_SOURCE)


# Build executables
$(foreach e,$(UNIT_EXECS),build/unit/$(e)): build/unit/Makefile $(MOM_SOURCE)
	cd $(@D) && $(TIME) $(MAKE) -j
build/%/MOM6: build/%/Makefile $(MOM_SOURCE)
	cd $(@D) && $(TIME) $(MAKE) -j


# Use autoconf to construct the Makefile for each target
.PRECIOUS: build/%/Makefile
build/%/Makefile: ../ac/configure ../ac/Makefile.in deps/lib/libFMS.a
	mkdir -p $(@D)
	cd $(@D) \
	  && $(MOM_ENV) ../../../ac/configure $(MOM_ACFLAGS) --with-framework=$(FRAMEWORK) \
	  || (cat config.log && false)


../ac/configure: ../ac/configure.ac ../ac/m4
	autoreconf -i $<


# Fetch the regression target codebase
build/target/Makefile build/opt_target/Makefile: \
  $(TARGET_CODEBASE)/ac/configure deps/lib/libFMS.a
	mkdir -p $(@D)
	cd $(@D) \
	  && $(MOM_ENV) ../../$(TARGET_CODEBASE)/ac/configure $(MOM_ACFLAGS) \
	  || (cat config.log && false)


$(TARGET_CODEBASE)/ac/configure: $(TARGET_CODEBASE)
	autoreconf -i $</ac


$(TARGET_CODEBASE):
	git clone --recursive $(MOM_TARGET_URL) $@
	cd $@ && git checkout --recurse-submodules $(MOM_TARGET_BRANCH)


#---
# FMS

# Set up the FMS build environment variables
FMS_ENV = \
  PATH="${PATH}:$(realpath ../ac)" \
  FCFLAGS="$(FCFLAGS_FMS)" \
  REPORT_ERROR_LOGS="$(REPORT_ERROR_LOGS)"

deps/lib/libFMS.a: deps/fms/build/libFMS.a
	$(MAKE) -C deps lib/libFMS.a

deps/fms/build/libFMS.a: deps/fms/build/Makefile
	$(MAKE) -C deps fms/build/libFMS.a

deps/fms/build/Makefile: deps/fms/src/configure deps/Makefile.fms.in
	$(FMS_ENV) $(MAKE) -C deps fms/build/Makefile

deps/Makefile.fms.in: ../ac/deps/Makefile.fms.in deps/Makefile
	cp $< deps

# TODO: m4 dependencies?
deps/fms/src/configure: ../ac/deps/configure.fms.ac deps/Makefile $(FMS_SOURCE) | deps/fms/src
	cp ../ac/deps/configure.fms.ac deps
	cp -r ../ac/deps/m4 deps
	$(MAKE) -C deps fms/src/configure

deps/fms/src: deps/Makefile
	make -C deps fms/src

# Dependency init
deps/Makefile: ../ac/deps/Makefile
	mkdir -p $(@D)
	cp $< $@


#---
# The following block does a non-library build of a coupled driver interface to
# MOM, along with everything below it.  This simply checks that we have not
# broken the ability to compile.  This is not a means to build a complete
# coupled executable.
# TODO:
# 	- Avoid re-building FMS and MOM6 src by re-using existing object/mod files
# 	- Use autoconf rather than mkmf templates
MK_TEMPLATE ?= ../../deps/mkmf/templates/ncrc-gnu.mk

# NUOPC driver
build/nuopc/mom_ocean_model_nuopc.o: build/nuopc/Makefile
	cd $(@D) && make $(@F)
check_mom6_api_nuopc: build/nuopc/mom_ocean_model_nuopc.o

# GFDL coupled driver
build/coupled/ocean_model_MOM.o: build/coupled/Makefile
	cd $(@D) && make $(@F)
check_mom6_api_coupled: build/coupled/ocean_model_MOM.o

# MCT driver
build/mct/mom_ocean_model_mct.o: build/mct/Makefile
	cd $(@D) && make $(@F)
check_mom6_api_mct: build/mct/mom_ocean_model_mct.o


#---
# Python preprocessing

# NOTE: Some less mature environments (e.g. Arm64 Ubuntu) require explicit
#   installation of numpy before netCDF4, as well as wheel and cython support.
work/local-env:
	python3 -m venv $@
	. $@/bin/activate \
	  && python3 -m pip install --upgrade pip \
	  && pip3 install wheel \
	  && pip3 install cython \
	  && pip3 install numpy \
	  && pip3 install netCDF4


#---
# Testing

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
run.symmetric: $(foreach c,$(CONFIGS),work/$(c)/symmetric/ocean.stats)
run.asymmetric: $(foreach c,$(filter-out tc3,$(CONFIGS)),$(CONFIGS),work/$(c)/asymmetric/ocean.stats)
run.nan: $(foreach c,$(CONFIGS),work/$(c)/nan/ocean.stats)
run.openmp: $(foreach c,$(CONFIGS),work/$(c)/openmp/ocean.stats)
run.cov: $(foreach c,$(CONFIGS),work/$(c)/cov/ocean.stats)

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
.PRECIOUS: $(foreach b,$(3),work/$(1)/$(b)/ocean.stats)
$(1).$(2): $(foreach b,$(3),work/$(1)/$(b)/ocean.stats)
	@test "$$(shell ls -A results/$(1) 2>/dev/null)" || rm -rf results/$(1)
	@cmp $$^ || !( \
	  mkdir -p results/$(1); \
	  (diff $$^ | tee results/$(1)/ocean.stats.$(2).diff | head -n 20) ; \
	  echo -e "$(FAIL): Solutions $(1).$(2) have changed." \
	)
	@echo -e "$(PASS): Solutions $(1).$(2) agree."

.PRECIOUS: $(foreach b,$(3),work/$(1)/$(b)/chksum_diag)
$(1).$(2).diag: $(foreach b,$(3),work/$(1)/$(b)/chksum_diag)
	@cmp $$^ || !( \
	  mkdir -p results/$(1); \
	  (diff $$^ | tee results/$(1)/chksum_diag.$(2).diff | head -n 20) ; \
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
.PRECIOUS: $(foreach b,symmetric restart target,work/%/$(b)/ocean.stats)
%.restart: $(foreach b,symmetric restart,work/%/$(b)/ocean.stats)
	@test "$(shell ls -A results/$* 2>/dev/null)" || rm -rf results/$*
	@cmp $(foreach f,$^,<(tr -s ' ' < $(f) | cut -d ' ' -f3- | tail -n 1)) \
	  || !( \
	    mkdir -p results/$*; \
	    (diff $^ | tee results/$*/chksum_diag.restart.diff | head -n 20) ; \
	    echo -e "$(FAIL): Solutions $*.restart have changed." \
	  )
	@echo -e "$(PASS): Solutions $*.restart agree."

# TODO: chksum_diag parsing of restart files

# stats rule is unchanged, but we cannot use CMP_RULE to generate it.
%.regression: $(foreach b,symmetric target,work/%/$(b)/ocean.stats)
	@test "$(shell ls -A results/$* 2>/dev/null)" || rm -rf results/$*
	@cmp $^ || !( \
	  mkdir -p results/$*; \
	  (diff $^ | tee results/$*/ocean.stats.regression.diff | head -n 20) ; \
	  echo -e "$(FAIL): Solutions $*.regression have changed." \
	)
	@echo -e "$(PASS): Solutions $*.regression agree."

# Regression testing only checks for changes in existing diagnostics
%.regression.diag: $(foreach b,symmetric target,work/%/$(b)/chksum_diag)
	@! diff $^ | grep "^[<>]" | grep "^>" > /dev/null \
	  || ! (\
	    mkdir -p results/$*; \
	    (diff $^ | tee results/$*/chksum_diag.regression.diff | head -n 20) ; \
	    echo -e "$(FAIL): Diagnostics $*.regression.diag have changed." \
	  )
	@cmp $^ || ( \
	  diff $^ | head -n 20; \
	  echo -e "$(WARN): New diagnostics in $<" \
	)
	@echo -e "$(PASS): Diagnostics $*.regression.diag agree."


#---
# Test run output files

# Rule to build work/<tc>/{ocean.stats,chksum_diag}.<tag>
# $(1): Test configuration name <tag>
# $(2): Executable type
# $(3): Enable coverage flag
# $(4): MOM_override configuration
# $(5): Environment variables
# $(6): Number of MPI ranks

define STAT_RULE
work/%/$(1)/ocean.stats work/%/$(1)/chksum_diag: build/$(2)/MOM6 $(VENV_PATH)
	@echo "Running test $$*.$(1)..."
	mkdir -p $$(@D)
	cp -RL $$*/* $$(@D)
	if [ -f $$(@D)/Makefile ]; then \
	  $$(VENV_ACTIVATE) \
	    && cd $$(@D) \
	    && $(MAKE); \
	else \
	  cd $$(@D); \
	fi
	mkdir -p $$(@D)/RESTART
	echo -e "$(4)" > $$(@D)/MOM_override
	rm -f results/$$*/std.$(1).{out,err}
	cd $$(@D) \
	  && $(TIME) $(5) $(MPIRUN) -n $(6) ../../../$$< 2> std.err > std.out \
	  || !( \
	    mkdir -p ../../../results/$$*/ ; \
	    cat std.out | tee ../../../results/$$*/std.$(1).out | tail -n 20 ; \
	    cat std.err | tee ../../../results/$$*/std.$(1).err | tail -n 20 ; \
	    rm ocean.stats chksum_diag ; \
	    echo -e "$(FAIL): $$*.$(1) failed at runtime." \
	  )
	@echo -e "$(DONE): $$*.$(1); no runtime errors."
	if [ $(3) ]; then \
	  mkdir -p results/$$* ; \
	  cd build/$(2) ; \
	  gcov -b *.gcda > gcov.$$*.$(1).out ; \
	  find -name "*.gcov" -exec sed -i -r 's/^( *[0-9]*)\*:/ \1:/g' {} \; ; \
	fi
endef


# Upload coverage reports
CODECOV_UPLOADER_URL ?= https://uploader.codecov.io/latest/linux/codecov
codecov:
	curl -s $(CODECOV_UPLOADER_URL) -o $@
	chmod +x codecov

.PHONY: report.cov
report.cov: run.cov codecov
	./codecov -R build/cov -Z -f "*.gcov" \
	  > build/cov/codecov.out \
	  2> build/cov/codecov.err \
	  && echo -e "${MAGENTA}Report uploaded to codecov.${RESET}" \
	  || { \
	    cat build/cov/codecov.err ; \
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
work/%/restart/ocean.stats: build/symmetric/MOM6 $(VENV_PATH)
	rm -rf $(@D)
	mkdir -p $(@D)
	cp -RL $*/* $(@D)
	if [ -f $(@D)/Makefile ]; then \
	  $(VENV_ACTIVATE) \
	    && cd work/$*/restart \
	    && $(MAKE); \
	else \
	  cd work/$*/restart; \
	fi
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
	rm -f results/$*/std.restart{1,2}.{out,err}
	# Run the first half-period
	cd $(@D) && $(TIME) $(MPIRUN) -n 1 ../../../$< 2> std1.err > std1.out \
	  || !( \
	    cat std1.out | tee ../../../results/$*/std.restart1.out | tail -n 20 ; \
	    cat std1.err | tee ../../../results/$*/std.restart1.err | tail -n 20 ; \
	    echo -e "$(FAIL): $*.restart failed at runtime." \
	  )
	# Setup the next inputs
	cd $(@D) && rm -rf INPUT && mv RESTART INPUT
	mkdir $(@D)/RESTART
	cd $(@D) && sed -i -e "s/input_filename *= *'n'/input_filename = 'r'/g" input.nml
	# Run the second half-period
	cd $(@D) && $(TIME) $(MPIRUN) -n 1 ../../../$< 2> std2.err > std2.out \
	  || !( \
	    cat std2.out | tee ../../../results/$*/std.restart2.out | tail -n 20 ; \
	    cat std2.err | tee ../../../results/$*/std.restart2.err | tail -n 20 ; \
	    echo -e "$(FAIL): $*.restart failed at runtime." \
	  )

# TODO: Restart checksum diagnostics

#---
# Not a true rule; only call this after `make test` to summarize test results.
.PHONY: test.summary
test.summary:
	@if ls results/*/* &> /dev/null; then \
	  if ls results/*/std.*.err &> /dev/null; then \
	    echo "The following tests failed to complete:" ; \
	    ls results/*/std.*.out \
	      | awk '{split($$0,a,"/"); split(a[3],t,"."); v=t[2]; if(length(t)>3) v=v"."t[3]; print a[2],":",v}'; \
	  fi; \
	  if ls results/*/ocean.stats.*.diff &> /dev/null; then \
	    echo "The following tests report solution regressions:" ; \
	    ls results/*/ocean.stats.*.diff \
	      | awk '{split($$0,a,"/"); split(a[3],t,"."); v=t[3]; if(length(t)>4) v=v"."t[4]; print a[2],":",v}'; \
	  fi; \
	  if ls results/*/chksum_diag.*.diff &> /dev/null; then \
	    echo "The following tests report diagnostic regressions:" ; \
	    ls results/*/chksum_diag.*.diff \
	      | awk '{split($$0,a,"/"); split(a[3],t,"."); v=t[2]; if(length(t)>3) v=v"."t[3]; print a[2],":",v}'; \
	  fi; \
	  false ; \
	else \
	  echo -e "$(PASS): All tests passed!"; \
	fi


#---
# Unit test

# NOTE: Using file parser gcov report as a proxy for test completion
.PHONY: run.cov.unit
run.cov.unit: build/unit/MOM_file_parser_tests.F90.gcov

work/unit/std.out: build/unit/MOM_unit_tests
	if [ $(REPORT_COVERAGE) ]; then \
	  find build/unit -name *.gcda -exec rm -f '{}' \; ; \
	fi
	rm -rf $(@D)
	mkdir -p $(@D)
	cd $(@D) \
	  && $(TIME) $(MPIRUN) -n 1 ../../$< 2> std.err > std.out \
	  || !( \
	    cat std.out | tail -n 100 ; \
	    cat std.err | tail -n 100 ; \
	  )
	cd $(@D) \
	  && $(TIME) $(MPIRUN) -n 2 ../../$< 2> p2.std.err > p2.std.out \
	  || !( \
	    cat p2.std.out | tail -n 100 ; \
	    cat p2.std.err | tail -n 100 ; \
	  )

# NOTE: .gcov actually depends on .gcda, but .gcda is produced with std.out
# TODO: Replace work/unit/std.out with *.gcda?
build/unit/MOM_file_parser_tests.F90.gcov: work/unit/std.out
	cd $(@D) \
	  && gcov -b *.gcda > gcov.unit.out
	find $(@D) -name "*.gcov" -exec sed -i -r 's/^( *[0-9]*)\*:/ \1:/g' {} \;

.PHONY: report.cov.unit
report.cov.unit: build/unit/MOM_file_parser_tests.F90.gcov codecov
	./codecov -R build/unit -f "*.gcov" -Z -n "Unit tests" \
	    > build/unit/codecov.out \
	    2> build/unit/codecov.err \
	  && echo -e "${MAGENTA}Report uploaded to codecov.${RESET}" \
	  || { \
	    cat build/unit/codecov.err ; \
	    echo -e "${RED}Failed to upload report.${RESET}" ; \
		if [ "$(REQUIRE_COVERAGE_UPLOAD)" = true ] ; then false ; fi ; \
	  }


#---
# Profiling based on FMS clocks

PCONFIGS = p0

.PHONY: profile
profile: $(foreach p,$(PCONFIGS), prof.$(p))

.PHONY: prof.p0
prof.p0: work/p0/opt/clocks.json work/p0/opt_target/clocks.json
	python tools/compare_clocks.py $^

work/p0/%/clocks.json: work/p0/%/std.out
	python tools/parse_fms_clocks.py -d $(@D) $^ > $@

work/p0/opt/std.out: build/opt/MOM6
work/p0/opt_target/std.out: build/opt_target/MOM6

work/p0/%/std.out:
	mkdir -p $(@D)
	cp -RL p0/* $(@D)
	mkdir -p $(@D)/RESTART
	echo -e "" > $(@D)/MOM_override
	cd $(@D) \
	  && $(MPIRUN) -n 1 ../../../$< 2> std.err > std.out


#---
# Profiling based on perf output

# TODO: This expects the -e flag, can I handle it in the command?
PERF_EVENTS ?=

.PHONY: perf
perf: $(foreach p,$(PCONFIGS), perf.$(p))

.PHONY: prof.p0
perf.p0: work/p0/opt/profile.json work/p0/opt_target/profile.json
	python tools/compare_perf.py $^

work/p0/%/profile.json: work/p0/%/perf.data
	python tools/parse_perf.py -f $< > $@

work/p0/opt/perf.data: build/opt/MOM6
work/p0/opt_target/perf.data: build/opt_target/MOM6

work/p0/%/perf.data:
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


#----
# NOTE: These tests assert that we are in the .testing directory.

.PHONY: clean
clean: clean.build clean.stats
	@[ $$(basename $$(pwd)) = .testing ]
	rm -rf deps


.PHONY: clean.build
clean.build:
	@[ $$(basename $$(pwd)) = .testing ]
	rm -rf build


.PHONY: clean.stats
clean.stats:
	@[ $$(basename $$(pwd)) = .testing ]
	rm -rf work results
