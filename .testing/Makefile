SHELL = bash
.SUFFIXES:

# User-defined configuration
-include config.mk

# Default configurations
MPIRUN ?= mpirun
DO_REPRO_TESTS ?= true

#---
# Dependencies
DEPS = deps

# mkmf, list_paths (GFDL build toolchain)
MKMF_URL ?= https://github.com/NOAA-GFDL/mkmf.git
MKMF_COMMIT ?= master
LIST_PATHS := $(abspath $(DEPS)/mkmf/bin/list_paths)
MKMF := $(abspath $(DEPS)/mkmf/bin/mkmf)

# FMS framework
FMS_URL ?= https://github.com/NOAA-GFDL/FMS.git
FMS_COMMIT ?= 2019.01.01
FMS := $(DEPS)/fms

#---
# Build configuration

# Build settings
MKMF_CPP = "-Duse_libMPI -Duse_netCDF -DSPMD"

# Environment
# TODO: This info ought to be determined by CMake, automake, etc.
#MKMF_TEMPLATE ?= linux-ubuntu-xenial-gnu.mk
MKMF_TEMPLATE ?= deps/mkmf/templates/ncrc-gnu.mk
#MKMF_TEMPLATE ?= deps/mkmf/templates/ncrc-intel.mk

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

MOM_SOURCE = $(call SOURCE,../src) $(wildcard ../config_src/solo_driver/*.F90)
TARGET_SOURCE = $(call SOURCE,build/target_codebase/src) \
  $(wildcard build/target_codebase/config_src/solo_driver/*.F90)
FMS_SOURCE = $(call SOURCE,$(DEPS)/fms/src)

#---
# Rules

.PHONY: all build.regressions
all: $(foreach b,$(BUILDS),build/$(b)/MOM6)
build.regressions: $(foreach b,symmetric target,build/$(b)/MOM6)

# Executable
BUILD_TARGETS = MOM6 Makefile path_names
.PRECIOUS: $(foreach b,$(BUILDS),$(foreach f,$(BUILD_TARGETS),build/$(b)/$(f)))

# Conditionally build symmetric with coverage support
COVFLAG=$(if $(REPORT_COVERAGE),COVERAGE=1,)

build/target/MOM6: MOMFLAGS=NETCDF=3 DEBUG=1 INIT=1
build/symmetric/MOM6: MOMFLAGS=NETCDF=3 DEBUG=1 INIT=1 $(COVFLAG)
build/asymmetric/MOM6: MOMFLAGS=NETCDF=3 DEBUG=1 INIT=1
build/repro/MOM6: MOMFLAGS=NETCDF=3 REPRO=1
build/openmp/MOM6: MOMFLAGS=NETCDF=3 DEBUG=1 OPENMP=1 INIT=1

build/asymmetric/path_names: GRID_SRC=config_src/dynamic
build/%/path_names: GRID_SRC=config_src/dynamic_symmetric

build/%/MOM6: build/%/Makefile $(FMS)/lib/libfms.a
	make -C $(@D) $(MOMFLAGS) $(@F)

build/%/Makefile: build/%/path_names
	cp $(MKMF_TEMPLATE) $(@D)
	cd $(@D) && $(MKMF) \
	  -t $(notdir $(MKMF_TEMPLATE)) \
	  -o '-I ../../$(DEPS)/fms/build' \
	  -p MOM6 \
	  -l '../../$(DEPS)/fms/lib/libfms.a' \
	  -c $(MKMF_CPP) \
	  path_names

# NOTE: These path_names rules could be merged

build/target/path_names: $(LIST_PATHS) $(TARGET_CODEBASE) $(TARGET_SOURCE)
	mkdir -p $(@D)
	cd $(@D) && $(LIST_PATHS) -l \
	  ../../$(TARGET_CODEBASE)/src \
	  ../../$(TARGET_CODEBASE)/config_src/solo_driver \
	  ../../$(TARGET_CODEBASE)/$(GRID_SRC)

build/%/path_names: $(LIST_PATHS) $(MOM_SOURCE)
	mkdir -p $(@D)
	cd $(@D) && $(LIST_PATHS) -l \
	  ../../../src \
	  ../../../config_src/solo_driver \
	  ../../../$(GRID_SRC)

# Target repository for regression tests
$(TARGET_CODEBASE):
	git clone --recursive $(MOM_TARGET_URL) $@
	cd $@ && git checkout $(MOM_TARGET_BRANCH)


#----
# FMS build

$(FMS)/lib/libfms.a: $(FMS)/build/Makefile
	mkdir -p $(FMS)/lib
	cd $(FMS)/build && make NETCDF=3 DEBUG=1 ../lib/libfms.a

$(FMS)/build/Makefile: $(FMS)/build/path_names
	cp $(MKMF_TEMPLATE) $(@D)
	cd $(@D) && $(MKMF) \
	  -t $(notdir $(MKMF_TEMPLATE)) \
	  -p ../lib/libfms.a \
	  -c $(MKMF_CPP) \
	  path_names

$(FMS)/build/path_names: $(LIST_PATHS) $(FMS)/src $(FMS_SOURCE)
	mkdir -p $(@D)
	cd $(@D) && $(LIST_PATHS) -l ../src

$(FMS)/src:
	git clone $(FMS_URL) $@
	cd $@; git checkout $(FMS_COMMIT)


#---
# Build Toolchain

$(LIST_PATHS) $(MKMF):
	git clone $(MKMF_URL) $(DEPS)/mkmf
	cd $(DEPS)/mkmf; git checkout $(MKMF_COMMIT)


#----
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

# Color highlights for test results
RED=\033[0;31m
GREEN=\033[0;32m
RESET=\033[0m

DONE=${GREEN}DONE${RESET}
PASS=${GREEN}PASS${RESET}
FAIL=${RED}FAIL${RESET}

# Comparison rules
# $(1): Test type (grid, layout, &c.)
# $(2): Comparison targets (symmetric asymmetric, symmetric layout, &c.)
define CMP_RULE
.PRECIOUS: $(foreach b,$(2),work/%/$(b)/ocean.stats)
%.$(1): $(foreach b,$(2),work/%/$(b)/ocean.stats)
	@cmp $$^ || !( \
	  mkdir -p results/$$*; \
	  (diff $$^ | tee results/$$*/ocean.stats.$(1).diff | head) ; \
	  echo -e "${FAIL}: Solutions $$*.$(1) have changed." \
	)
	@echo -e "${PASS}: Solutions $$*.$(1) agree."

.PRECIOUS: $(foreach b,$(2),work/%/$(b)/chksum_diag)
%.$(1).diag: $(foreach b,$(2),work/%/$(b)/chksum_diag)
	@cmp $$^ || !( \
	  mkdir -p results/$$*; \
	  (diff $$^ | tee results/$$*/chksum_diag.$(1).diff | head) ; \
	  echo -e "${FAIL}: Diagnostics $$*.$(1).diag have changed." \
	)
	@echo -e "${PASS}: Diagnostics $$*.$(1).diag agree."
endef

$(eval $(call CMP_RULE,grid,symmetric asymmetric))
$(eval $(call CMP_RULE,layout,symmetric layout))
$(eval $(call CMP_RULE,rotate,symmetric rotate))
$(eval $(call CMP_RULE,repro,symmetric repro))
$(eval $(call CMP_RULE,openmp,symmetric openmp))
$(eval $(call CMP_RULE,nan,symmetric nan))
$(foreach d,$(DIMS),$(eval $(call CMP_RULE,dim.$(d),symmetric dim.$(d))))
$(eval $(call CMP_RULE,regression,symmetric target))

# Custom comparison rules

# Restart tests only compare the final stat record
.PRECIOUS: $(foreach b,symmetric restart target,work/%/$(b)/ocean.stats)
%.restart: $(foreach b,symmetric restart,work/%/$(b)/ocean.stats)
	#cmp $(foreach f,$^,<(tr -s ' ' < $(f) | cut -d ' ' -f3- | tail -n 1)) \
	#  || diff $^
	@cmp $(foreach f,$^,<(tr -s ' ' < $(f) | cut -d ' ' -f3- | tail -n 1)) \
	  || !( \
	    mkdir -p results/$*; \
	    (diff $$^ | tee results/$*/chksum_diag.restart.diff | head) ; \
	    echo -e "${FAIL}: Diagnostics $*.restart.diag have changed." \
	  )
	@echo -e "${PASS}: Diagnostics $*.restart.diag agree."

# TODO: chksum_diag parsing of restart files


#---
# Test run output files

# Generalized MPI environment variable support
# XXX: Using `-env` in the MPICH test can erroneously producing an `nv` file.
# $(1): Environment variables
ifeq ($(shell $(MPIRUN) -x tmp=1 true 2> /dev/null ; echo $$?), 0)
  MPIRUN_CMD=$(MPIRUN) $(if $(1),-x $(1),)
else ifeq ($(shell $(MPIRUN) -env tmp=1 true 2> /dev/null ; echo $$?), 0)
  MPIRUN_CMD=$(MPIRUN) $(if $(1),-env $(1),)
else
  MPIRUN_CMD=$(1) $(MPIRUN)
endif


# Rule to build work/<tc>/{ocean.stats,chksum_diag}.<tag>
# $(1): Test configuration name <tag>
# $(2): Executable type
# $(3): Enable coverage flag
# $(4): MOM_override configuration
# $(5): Environment variables
# $(6): Number of MPI ranks
define STAT_RULE
work/%/$(1)/ocean.stats work/%/$(1)/chksum_diag: build/$(2)/MOM6
	@echo "Running test $$*.$(1)..."
	if [ $(3) ]; then find build/$(2) -name *.gcda -exec rm -f '{}' \; ; fi
	mkdir -p $$(@D)
	cp -rL $$*/* $$(@D)
	cd $$(@D) && if [ -f Makefile ]; then make; fi
	mkdir -p $$(@D)/RESTART
	echo -e "$(4)" > $$(@D)/MOM_override
	cd $$(@D) \
	  && $$(call MPIRUN_CMD,$(5)) -n $(6) ../../../$$< 2> std.err > std.out \
	  || !( \
	    mkdir -p ../../../results/$$*/ ; \
	    cat std.out | tee ../../../results/$$*/std.$(1).out | tail ; \
		cat std.err | tee ../../../results/$$*/std.$(1).err | tail ; \
		rm ocean.stats chksum_diag ; \
		echo -e "${FAIL}: $$*.$(1) failed at runtime." \
	  )
	@echo -e "${DONE}: $$*.$(1); no runtime errors."
	if [ $(3) ]; then \
	  mkdir -p results/$$* ; \
	  bash <(curl -s https://codecov.io/bash) -n $$@ \
	    > work/$$*/codecov.$(1).out \
		2> work/$$*/codecov.$(1).err ; \
	fi
endef


# Define $(,) as comma escape character
, := ,

$(eval $(call STAT_RULE,symmetric,symmetric,$(REPORT_COVERAGE),,,1))
$(eval $(call STAT_RULE,asymmetric,asymmetric,,,,1))
$(eval $(call STAT_RULE,target,target,,,,1))
$(eval $(call STAT_RULE,repro,repro,,,,1))
$(eval $(call STAT_RULE,openmp,openmp,,,,1))
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
work/%/restart/ocean.stats: build/symmetric/MOM6
	rm -rf $(@D)
	mkdir -p $(@D)
	cp -rL $*/* $(@D)
	cd work/$*/restart && if [ -f Makefile ]; then make; fi
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
	# Run the first half-period
	cd $(@D) && $(MPIRUN) -n 1 ../../../$< 2> std1.err > std1.out \
	  || !( \
	    cat std1.out | tee ../../../results/$*/std.restart1.out | tail ; \
		cat std1.err | tee ../../../results/$*/std.restart1.err | tail ; \
		echo -e "${FAIL}: $*.restart failed at runtime." \
	  )
	# Setup the next inputs
	cd $(@D) && rm -rf INPUT && mv RESTART INPUT
	mkdir $(@D)/RESTART
	cd $(@D) && sed -i -e "s/input_filename *= *'n'/input_filename = 'r'/g" input.nml
	# Run the second half-period
	cd $(@D) && $(MPIRUN) -n 1 ../../../$< 2> std2.err > std2.out \
	  || !( \
	    cat std2.out | tee ../../../results/$*/std.restart2.out | tail ; \
		cat std2.err | tee ../../../results/$*/std.restart2.err | tail ; \
		echo -e "${FAIL}: $*.restart failed at runtime." \
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
	  echo -e "${PASS}: All tests passed!"; \
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
