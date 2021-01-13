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
#       Wipe the MOM6 test executables
#       (NOTE: This does not delete FMS in the `deps`)
#
#
# Configuration:
#   These settings can be provided as either command-line flags, or saved in a
#   `config.mk` file.
#
# Test suite configuration:
#
#   MPIRUN                  MPI job launcher (mpirun, srun, etc)
#   DO_REPRO_TESTS          Enable production ("repro") testing equivalence
#   DO_REGRESSION_TESTS:    Enable regression tests (usually dev/gfdl)
#   REPORT_COVERAGE         Enable code coverage and report to codecov
#
# Compiler configuration:
#   (NOTE: These are environment variables and may be inherited from a shell.)
#
#   CC      C compiler
#   MPICC   MPI C compiler
#   FC      Fortran compiler
#   MPIFC   MPI Fortran compiler
#
# Build configuration:
#
#   FCFLAGS_DEBUG       Testing ("debug") compiler flags
#   FCFLAGS_REPRO       Production ("repro") compiler flags
#   FCFLAGS_INIT        Variable initialization flags
#   FCFLAGS_COVERAGE    Code coverage flags
#
# Regression repository ("target") configuration:
#   (NOTE: These would typically be configured by a CI such as Travis.)
#
#   MOM_TARGET_SLUG             URL slug (minus domain) of the target repo
#   MOM_TARGET_URL              Full URL of the target repo
#   MOM_TARGET_LOCAL_BRANCH     Target branch name
#
#----

# TODO: Bourne shell compatibility
SHELL = bash

# No implicit rules
.SUFFIXES:

# No implicit variables
MAKEFLAGS += -R

# User-defined configuration
-include config.mk

# Set the MPI launcher here
# TODO: This needs more automated configuration
MPIRUN ?= mpirun

# Generic compiler variables are pass through to the builds
export CC
export MPICC
export FC
export MPIFC

# Builds are distinguished by FCFLAGS
# NOTE: FMS will be built using FCFLAGS_DEBUG
FCFLAGS_DEBUG ?= -g -O0
FCFLAGS_REPRO ?= -g -O2
FCFLAGS_INIT ?=
FCFLAGS_COVERAGE ?=
# Additional notes:
#
# - The default values are simple, minimalist flags, supported by nearly all
#   compilers which are comparable to GFDL's canonical DEBUG and REPRO builds.
#
# - These flags should be configured outside of the Makefile, either with
#   config.mk or as environment variables.
#
# - FMS cannot be built with the same aggressive initialization flags as MOM6,
#   so FCFLAGS_INIT is used to provide additional MOM6 configuration.

# Set to `true` to require identical results from DEBUG and REPRO builds
DO_REPRO_TESTS ?=

TIME ?= time

# Many compilers (Intel, GCC on ARM64) do not yet produce identical results
# across DEBUG and REPRO builds (as defined below), so we disable on default.

#---
# Dependencies
DEPS = deps

# mkmf, list_paths (GFDL build toolchain)
LIST_PATHS := $(DEPS)/bin/list_paths
MKMF := $(DEPS)/bin/mkmf


#---
# Test configuration

# Executables
BUILDS = symmetric asymmetric repro openmp
CONFIGS := $(wildcard tc*)
TESTS = grids layouts restarts nans dims openmps rotations
DIMS = t l h z q r

# REPRO tests enable reproducibility with optimization, and often do not match
# the DEBUG results in older GCCs and vendor compilers, so we can optionally
# disable them.
ifeq ($(DO_REPRO_TESTS), true)
  BUILDS += repro
  TESTS += repros
endif

# The following variables are configured by Travis:
#   DO_REGRESSION_TESTS: true if $(TRAVIS_PULL_REQUEST) is a PR number
#   MOM_TARGET_SLUG: TRAVIS_REPO_SLUG
#   MOM_TARGET_LOCAL_BRANCH: TRAVIS_BRANCH
# These are set to true by our Travis configuration if testing a pull request
DO_REGRESSION_TESTS ?=
REPORT_COVERAGE ?=

ifeq ($(DO_REGRESSION_TESTS), true)
  BUILDS += target
  TESTS += regressions

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

MOM_SOURCE = $(call SOURCE,../src) \
  $(wildcard ../config_src/solo_driver/*.F90) \
  $(wildcard ../config_src/ext*/*/*.F90)
TARGET_SOURCE = $(call SOURCE,build/target_codebase/src) \
  $(wildcard build/target_codebase/config_src/solo_driver/*.F90) \
  $(wildcard build/target_codebase/config_src/ext*/*.F90)
FMS_SOURCE = $(call SOURCE,$(DEPS)/fms/src)


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

.PHONY: all build.regressions
all: $(foreach b,$(BUILDS),build/$(b)/MOM6)
build.regressions: $(foreach b,symmetric target,build/$(b)/MOM6)

# Executable
BUILD_TARGETS = MOM6 Makefile path_names
.PRECIOUS: $(foreach b,$(BUILDS),$(foreach f,$(BUILD_TARGETS),build/$(b)/$(f)))

# Compiler flags

# Conditionally build symmetric with coverage support
COVERAGE=$(if $(REPORT_COVERAGE),$(FCFLAGS_COVERAGE),)

# .testing dependencies
# TODO: We should probably build TARGET with the FMS that it was configured
#   to use.  But for now we use the same FMS over all builds.
FCFLAGS_FMS = -I../../$(DEPS)/include
LDFLAGS_FMS = -L../../$(DEPS)/lib
PATH_FMS = PATH="${PATH}:../../$(DEPS)/bin"


# Define the build targets in terms of the traditional DEBUG/REPRO/etc labels
SYMMETRIC_FCFLAGS := FCFLAGS="$(FCFLAGS_DEBUG) $(FCFLAGS_INIT) $(COVERAGE) $(FCFLAGS_FMS)"
ASYMMETRIC_FCFLAGS := FCFLAGS="$(FCFLAGS_DEBUG) $(FCFLAGS_INIT) $(FCFLAGS_FMS)"
REPRO_FCFLAGS := FCFLAGS="$(FCFLAGS_REPRO) $(FCFLAGS_FMS)"
OPENMP_FCFLAGS := FCFLAGS="$(FCFLAGS_DEBUG) $(FCFLAGS_INIT) $(FCFLAGS_FMS)"
TARGET_FCFLAGS := FCFLAGS="$(FCFLAGS_DEBUG) $(FCFLAGS_INIT) $(FCFLAGS_FMS)"

MOM_LDFLAGS := LDFLAGS="$(LDFLAGS_FMS)"
SYMMETRIC_LDFLAGS := LDFLAGS="$(COVERAGE) $(LDFLAGS_FMS)"


# Environment variable configuration
build/symmetric/Makefile: MOM_ENV=$(PATH_FMS) $(SYMMETRIC_FCFLAGS) $(SYMMETRIC_LDFLAGS)
build/asymmetric/Makefile: MOM_ENV=$(PATH_FMS) $(ASYMMETRIC_FCFLAGS) $(MOM_LDFLAGS)
build/repro/Makefile: MOM_ENV=$(PATH_FMS) $(REPRO_FCFLAGS) $(MOM_LDFLAGS)
build/openmp/Makefile: MOM_ENV=$(PATH_FMS) $(OPENMP_FCFLAGS) $(MOM_LDFLAGS)
build/target/Makefile: MOM_ENV=$(PATH_FMS) $(TARGET_FCFLAGS) $(MOM_LDFLAGS)
build/coupled/Makefile: MOM_ENV=$(PATH_FMS) $(SYMMETRIC_FCFLAGS) $(SYMMETRIC_LDFLAGS)
build/nuopc/Makefile: MOM_ENV=$(PATH_FMS) $(SYMMETRIC_FCFLAGS) $(SYMMETRIC_LDFLAGS)
build/mct/Makefile: MOM_ENV=$(PATH_FMS) $(SYMMETRIC_FCFLAGS) $(SYMMETRIC_LDFLAGS)

# Configure script flags
build/symmetric/Makefile: MOM_ACFLAGS=
build/asymmetric/Makefile: MOM_ACFLAGS=--enable-asymmetric
build/repro/Makefile: MOM_ACFLAGS=
build/openmp/Makefile: MOM_ACFLAGS=--enable-openmp
build/target/Makefile: MOM_ACFLAGS=
build/coupled/Makefile: MOM_ACFLAGS=--with-driver=coupled_driver
build/nuopc/Makefile: MOM_ACFLAGS=--with-driver=nuopc_driver
build/mct/Makefile: MOM_ACFLAGS=--with-driver=mct_driver

# Fetch regression target source code
build/target/Makefile: | $(TARGET_CODEBASE)


# Define source code dependencies
# NOTE: ./configure is too much, but Makefile is not enough!
#   Ideally we would want to re-run both Makefile and mkmf, but our mkmf call
#   is inside ./configure, so we must re-run ./configure as well.
$(foreach b,$(filter-out target,$(BUILDS)),build/$(b)/Makefile): $(MOM_SOURCE)
build/target/configure: $(TARGET_SOURCE)


# Build MOM6
.PRECIOUS: $(foreach b,$(BUILDS),build/$(b)/MOM6)
build/%/MOM6: build/%/Makefile
	cd $(@D) && $(TIME) $(MAKE) -j


# Use autoconf to construct the Makefile for each target
.PRECIOUS: $(foreach b,$(BUILDS),build/$(b)/Makefile)
build/%/Makefile: ../ac/configure ../ac/Makefile.in $(DEPS)/lib/libFMS.a $(MKMF) $(LIST_PATHS)
	mkdir -p $(@D)
	cd $(@D) \
	  && $(MOM_ENV) ../../../ac/configure $(MOM_ACFLAGS) \
	  || (cat config.log && false)


../ac/configure: ../ac/configure.ac ../ac/m4
	autoreconf -i $<


# Fetch the regression target codebase
build/target/Makefile: $(TARGET_CODEBASE)/ac/configure $(DEPS)/lib/libFMS.a $(MKMF) $(LIST_PATHS)
	mkdir -p $(@D)
	cd $(@D) \
	  && $(MOM_ENV) ../../$(TARGET_CODEBASE)/ac/configure $(MOM_ACFLAGS) \
	  || (cat config.log && false)


$(TARGET_CODEBASE)/ac/configure: $(TARGET_CODEBASE)
	autoreconf -i $</ac


$(TARGET_CODEBASE):
	git clone --recursive $(MOM_TARGET_URL) $@
	cd $@ && git checkout $(MOM_TARGET_BRANCH)
	# Copy modern autoconf files to target?
	mkdir -p $(TARGET_CODEBASE)/ac
	cp -r ../ac/{configure.ac,Makefile.in,m4} $(TARGET_CODEBASE)/ac


#---
# FMS

# Set up the FMS build environment variables
FMS_ENV = PATH="${PATH}:../../bin" FCFLAGS="$(FCFLAGS_DEBUG)"

# TODO: *.mod dependencies?
$(DEPS)/lib/libFMS.a: $(DEPS)/fms/build/libFMS.a
	$(MAKE) -C $(DEPS) lib/libFMS.a

$(DEPS)/fms/build/libFMS.a: $(DEPS)/fms/build/Makefile
	$(MAKE) -C $(DEPS) fms/build/libFMS.a

$(DEPS)/fms/build/Makefile: $(DEPS)/fms/src/configure $(DEPS)/Makefile.fms.in $(MKMF) $(LIST_PATHS)
	$(FMS_ENV) $(MAKE) -C $(DEPS) fms/build/Makefile
	$(MAKE) -C $(DEPS) fms/build/Makefile

$(DEPS)/Makefile.fms.in: ../ac/deps/Makefile.fms.in $(DEPS)/Makefile
	cp $< $(DEPS)

# TODO: m4 dependencies?
$(DEPS)/fms/src/configure: ../ac/deps/configure.fms.ac $(DEPS)/Makefile $(FMS_SOURCE) | $(DEPS)/fms/src
	cp ../ac/deps/configure.fms.ac $(DEPS)
	cp -r ../ac/deps/m4 $(DEPS)
	$(MAKE) -C $(DEPS) fms/src/configure

$(DEPS)/fms/src: $(DEPS)/Makefile
	make -C $(DEPS) fms/src

# mkmf
$(MKMF) $(LIST_PATHS): $(DEPS)/mkmf
	$(MAKE) -C $(DEPS) bin/$(@F)

$(DEPS)/mkmf: $(DEPS)/Makefile
	$(MAKE) -C $(DEPS) mkmf

# Dependency init
$(DEPS)/Makefile: ../ac/deps/Makefile
	mkdir -p $(@D)
	cp $< $@

#---
# The following block does a non-library build of a coupled driver interface to MOM, along with everything below it.
# This simply checks that we have not broken the ability to compile. This is not a means to build a complete coupled executable.
# Todo: 
# - avoid re-building FMS and MOM6 src by re-using existing object/mod files
# - use autoconf rather than mkmf templates
MK_TEMPLATE ?= ../../$(DEPS)/mkmf/templates/ncrc-gnu.mk
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
test.grids: $(foreach c,$(filter-out tc3,$(CONFIGS)),$(c).grid $(c).grid.diag)
test.layouts: $(foreach c,$(CONFIGS),$(c).layout $(c).layout.diag)
test.rotations: $(foreach c,$(CONFIGS),$(c).rotate)
test.restarts: $(foreach c,$(CONFIGS),$(c).restart)
test.repros: $(foreach c,$(CONFIGS),$(c).repro $(c).repro.diag)
test.openmps: $(foreach c,$(CONFIGS),$(c).openmp $(c).openmp.diag)
test.nans: $(foreach c,$(CONFIGS),$(c).nan $(c).nan.diag)
test.dims: $(foreach c,$(CONFIGS),$(foreach d,$(DIMS),$(c).dim.$(d) $(c).dim.$(d).diag))
test.regressions: $(foreach c,$(CONFIGS),$(c).regression $(c).regression.diag)

.PHONY: run.symmetric run.asymmetric run.nans run.openmp
run.symmetric: $(foreach c,$(CONFIGS),work/$(c)/symmetric/ocean.stats)
run.asymmetric: $(foreach c,$(filter-out tc3,$(CONFIGS)),$(CONFIGS),work/$(c)/asymmetric/ocean.stats)
run.nans: $(foreach c,$(CONFIGS),work/$(c)/nan/ocean.stats)
run.openmp: $(foreach c,$(CONFIGS),work/$(c)/openmp/ocean.stats)

# Color highlights for test results
RED = \033[0;31m
YELLOW = \033[0;33m
GREEN = \033[0;32m
MAGENTA = \033[0;35m
RESET = \033[0m

DONE = ${GREEN}DONE${RESET}
PASS = ${GREEN}PASS${RESET}
WARN = ${YELLOW}WARN${RESET}
FAIL = ${RED}FAIL${RESET}

# Comparison rules
# $(1): Test type (grid, layout, &c.)
# $(2): Comparison targets (symmetric asymmetric, symmetric layout, &c.)
define CMP_RULE
.PRECIOUS: $(foreach b,$(2),work/%/$(b)/ocean.stats)
%.$(1): $(foreach b,$(2),work/%/$(b)/ocean.stats)
	@test "$$(shell ls -A results/$$* 2>/dev/null)" || rm -rf results/$$*
	@cmp $$^ || !( \
	  mkdir -p results/$$*; \
	  (diff $$^ | tee results/$$*/ocean.stats.$(1).diff | head -n 20) ; \
	  echo -e "$(FAIL): Solutions $$*.$(1) have changed." \
	)
	@echo -e "$(PASS): Solutions $$*.$(1) agree."

.PRECIOUS: $(foreach b,$(2),work/%/$(b)/chksum_diag)
%.$(1).diag: $(foreach b,$(2),work/%/$(b)/chksum_diag)
	@cmp $$^ || !( \
	  mkdir -p results/$$*; \
	  (diff $$^ | tee results/$$*/chksum_diag.$(1).diff | head -n 20) ; \
	  echo -e "$(FAIL): Diagnostics $$*.$(1).diag have changed." \
	)
	@echo -e "$(PASS): Diagnostics $$*.$(1).diag agree."
endef

$(eval $(call CMP_RULE,grid,symmetric asymmetric))
$(eval $(call CMP_RULE,layout,symmetric layout))
$(eval $(call CMP_RULE,rotate,symmetric rotate))
$(eval $(call CMP_RULE,repro,symmetric repro))
$(eval $(call CMP_RULE,openmp,symmetric openmp))
$(eval $(call CMP_RULE,nan,symmetric nan))
$(foreach d,$(DIMS),$(eval $(call CMP_RULE,dim.$(d),symmetric dim.$(d))))

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
	if [ $(3) ]; then find build/$(2) -name *.gcda -exec rm -f '{}' \; ; fi
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
	  bash <(curl -s https://codecov.io/bash) -n $$@ \
	    > work/$$*/codecov.$(1).out \
	    2> work/$$*/codecov.$(1).err \
	    && echo -e "${MAGENTA}Report uploaded to codecov.${RESET}"; \
	fi
endef


# Define $(,) as comma escape character
, := ,

$(eval $(call STAT_RULE,symmetric,symmetric,$(REPORT_COVERAGE),,,1))
$(eval $(call STAT_RULE,asymmetric,asymmetric,,,,1))
$(eval $(call STAT_RULE,target,target,,,,1))
$(eval $(call STAT_RULE,repro,repro,,,,1))
$(eval $(call STAT_RULE,openmp,openmp,,,GOMP_CPU_AFFINITY=0,1))
$(eval $(call STAT_RULE,layout,symmetric,,LAYOUT=2$(,)1,,2))
$(eval $(call STAT_RULE,rotate,symmetric,,ROTATE_INDEX=True\nINDEX_TURNS=1,,1))
$(eval $(call STAT_RULE,nan,symmetric,,,MALLOC_PERTURB_=256,1))
$(eval $(call STAT_RULE,dim.t,symmetric,,T_RESCALE_POWER=11,,1))
$(eval $(call STAT_RULE,dim.l,symmetric,,L_RESCALE_POWER=11,,1))
$(eval $(call STAT_RULE,dim.h,symmetric,,H_RESCALE_POWER=11,,1))
$(eval $(call STAT_RULE,dim.z,symmetric,,Z_RESCALE_POWER=11,,1))
$(eval $(call STAT_RULE,dim.q,symmetric,,Q_RESCALE_POWER=11,,1))
$(eval $(call STAT_RULE,dim.r,symmetric,,R_RESCALE_POWER=11,,1))


# Restart tests require significant preprocessing, and are handled separately.
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
	# Generate the half-period input namelist
	# TODO: Assumes that runtime set by DAYMAX, will fail if set by input.nml
	cd $(@D) \
	  && daymax=$$(grep DAYMAX MOM_input | cut -d '!' -f 1 | cut -d '=' -f 2 | xargs) \
	  && timeunit=$$(grep TIMEUNIT MOM_input | cut -d '!' -f 1 | cut -d '=' -f 2 | xargs) \
	  && if [ -z "$${timeunit}" ]; then timeunit="8.64e4"; fi \
	  && printf -v timeunit_int "%.f" "$${timeunit}" \
	  && halfperiod=$$(printf "%.f" $$(bc <<< "scale=10; 0.5 * $${daymax} * $${timeunit_int}")) \
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


#----
# NOTE: These tests assert that we are in the .testing directory.

.PHONY: clean
clean: clean.stats
	@[ $$(basename $$(pwd)) = .testing ]
	rm -rf build


.PHONY: clean.stats
clean.stats:
	@[ $$(basename $$(pwd)) = .testing ]
	rm -rf work results
