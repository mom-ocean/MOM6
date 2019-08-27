SHELL = bash
MPIRUN ?= mpirun

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
FMS_COMMIT ?= f2e2c86f6c0eb6d389a20509a8a60fa22924e16b
FMS := $(DEPS)/fms

#---
# Build configuration

# Build settings
MKMF_CPP = "-Duse_libMPI -Duse_netCDF -DSPMD"

# Environment
# TODO: This info ought to be determined by CMake, automake, etc.
#MKMF_TEMPLATE ?= .testing/linux-ubuntu-xenial-gnu.mk
MKMF_TEMPLATE ?= $(DEPS)/mkmf/templates/ncrc-gnu.mk
#MKMF_TEMPLATE ?= $(DEPS)/mkmf/templates/ncrc-intel.mk

#---
# Test configuration

# Executables
BUILDS = symmetric asymmetric repro
CONFIGS := $(foreach n,$(shell seq 0 3),tc$(n))
TESTS = grids layouts restarts repros nans dims

# The following variables are configured by Travis:
#   DO_REGRESSION_TESTS: true if $(TRAVIS_PULL_REQUEST) is a PR number
#   MOM_TARGET_SLUG: TRAVIS_REPO_SLUG
#   MOM_TARGET_LOCAL_BRANCH: TRAVIS_BRANCH
#
# MOM_TARGET_REMOTE_NAME and MOM_MERGED_BRANCH are arbitrary labels which may
# cause namespace conflicts, and should probably be renamed or phased out.

# These are set to true by Travis if testing a pull request
DO_REGRESSION_TESTS ?=
REPORT_COVERAGE ?=

ifeq ($(DO_REGRESSION_TESTS), true)
	BUILDS += target
	TEST += regressions

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


#---
# Rules

.PHONY: all
all: $(foreach b,$(BUILDS),build/$(b)/MOM6)

# Executable
BUILD_TARGETS = MOM6 Makefile path_names
.PRECIOUS: $(foreach b,$(BUILDS),$(foreach f,$(BUILD_TARGETS),build/$(b)/$(f)))

build/target/MOM6: MOMFLAGS=NETCDF=3 DEBUG=1
build/symmetric/MOM6: MOMFLAGS=NETCDF=3 DEBUG=1 COVERAGE=1
build/asymmetric/MOM6: MOMFLAGS=NETCDF=3 DEBUG=1
build/repro/MOM6: MOMFLAGS=NETCDF=3 REPRO=1

build/asymmetric/path_names: GRID_SRC=config_src/dynamic
build/%/path_names: GRID_SRC=config_src/dynamic_symmetric

build/%/MOM6: build/%/Makefile $(FMS)/lib/libfms.a
	make -C $(@D) $(MOMFLAGS) $(@F)

build/%/Makefile: build/%/path_names
	cp $(MKMF_TEMPLATE) $(@D)
	cd $(@D) && $(MKMF) \
		-t $(notdir $(MKMF_TEMPLATE)) \
		-o '-I ../../$(FMS)/build' \
		-p MOM6 \
		-l '../../$(FMS)/lib/libfms.a' \
		-c $(MKMF_CPP) \
		path_names

# NOTE: These path_names rules could be merged

build/target/path_names: $(LIST_PATHS) $(TARGET_CODEBASE)
	mkdir -p $(@D)
	cd $(@D) && $(LIST_PATHS) -l \
		../../$(TARGET_CODEBASE)/src \
		../../$(TARGET_CODEBASE)/config_src/solo_driver \
		../../$(TARGET_CODEBASE)/$(GRID_SRC)

build/%/path_names: $(LIST_PATHS)
	mkdir -p $(@D)
	cd $(@D) && $(LIST_PATHS) -l \
		../../src \
		../../config_src/solo_driver \
		../../$(GRID_SRC)

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

$(FMS)/build/path_names: $(FMS)/src $(FMS_FILES) $(LIST_PATHS)
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

# NOTE: We remove tc3 (OBC) from grid test

.PHONY: $(foreach t,$(TESTS),test.$(t))
test.regressions: $(foreach c,$(CONFIGS),$(c).regression $(c).regression.diag)
test.grids: $(foreach c,$(filter-out tc3,$(CONFIGS)),$(c).grid $(c).grid.diag)
test.layouts: $(foreach c,$(CONFIGS),$(c).layout $(c).layout.diag)
test.restarts: $(foreach c,$(CONFIGS),$(c).restart)
test.repros: $(foreach c,$(CONFIGS),$(c).repro $(c).repro.diag)
test.nans: $(foreach c,$(CONFIGS),$(c).nan $(c).nan.diag)
test.dims: $(foreach c,$(CONFIGS),$(foreach d,t l h z,$(c).dim.$(d) $(c).dim.$(d).diag))

# NOTE: chksum_diag return code of cmp is currently ignored since many fail!
define CMP_RULE
%.$(1): $(foreach b,$(2),.testing/%/ocean.stats.$(b))
	cmp $$^

%.$(1).diag: $(foreach b,$(2),.testing/%/chksum_diag.$(b))
	-cmp $$^
endef

$(eval $(call CMP_RULE,regression,symmetric target))
$(eval $(call CMP_RULE,grid,symmetric asymmetric))
$(eval $(call CMP_RULE,layout,symmetric layout))
$(eval $(call CMP_RULE,repro,symmetric repro))
$(eval $(call CMP_RULE,nan,symmetric nan))
$(foreach d,t l h z,$(eval $(call CMP_RULE,dim.$(d),symmetric dim.$(d))))

# Restart tests only compare the final stat record
%.restart: $(foreach b,symmetric restart,.testing/%/ocean.stats.$(b))
	cmp $(foreach f,$^,<(tr -s ' ' < $(f) | cut -d ' ' -f3- | tail -n 1))

# TODO: chksum_diag parsing of restart files


#---
# Test run output files

#(1): Configuration name
#(2): Executable type
#(3): Enable coverage flag
#(4): MOM_override configuration
#(5): Environment variables
#(6): Number of MPI ranks

# Simple environment variable support for Slurm and OpenMPI
# This will probably need to be generalized in the future.
ifeq ($(MPIRUN), srun)
MPIRUN_CMD=$(1) $(MPIRUN)
else
MPIRUN_CMD=$(MPIRUN) $(if $(1),-x $(1),)
endif

define STAT_RULE
.testing/%/ocean.stats.$(1): build/$(2)/MOM6
	if [ $(3) ]; then find build -name *.gcda -exec rm -f '{}' \; ; fi
	mkdir -p $$(@D)/RESTART
	echo $(4) > $$(@D)/MOM_override
	cd $$(@D) && $$(call MPIRUN_CMD,$(5)) -n $(6) ../../$$< 2> debug.out
	cp $$(@D)/ocean.stats $$@
	> $$(@D)/MOM_override
	if [ $(3) ]; then bash <(curl -s https://codecov.io/bash) -n $$@; fi

.testing/%/chksum_diag.$(1): .testing/%/ocean.stats.$(1)
	cp $$(@D)/chksum_diag $$@
endef

# Define $(,) as comma escape character
, := ,

$(eval $(call STAT_RULE,symmetric,symmetric,$(REPORT_COVERAGE),,,1))
$(eval $(call STAT_RULE,asymmetric,asymmetric,,,,1))
$(eval $(call STAT_RULE,target,target,,,,1))
$(eval $(call STAT_RULE,repro,repro,,,,1))
$(eval $(call STAT_RULE,layout,symmetric,,LAYOUT=2$(,)1,,2))
$(eval $(call STAT_RULE,nan,symmetric,,,MALLOC_PERTURB_=256,1))
$(eval $(call STAT_RULE,dim.t,symmetric,,T_RESCALE_POWER=11,,1))
$(eval $(call STAT_RULE,dim.l,symmetric,,L_RESCALE_POWER=11,,1))
$(eval $(call STAT_RULE,dim.h,symmetric,,H_RESCALE_POWER=11,,1))
$(eval $(call STAT_RULE,dim.z,symmetric,,Z_RESCALE_POWER=11,,1))

# Restart tests require signicant preprocessing, and are handled separately.
.testing/%/ocean.stats.restart: build/symmetric/MOM6
	# Cleanup
	mkdir -p $(@D)/RESTART
	git checkout $(@D)/input.nml
	> $(@D)/MOM_override
	# Generate the half-period input namelist
	# TODO: Assumes runtime set by DAYMAX, will fail if set by input.nml
	cd $(@D) \
		&& daymax=$$(grep DAYMAX MOM_input | cut -d '!' -f 1 | cut -d '=' -f 2 | xargs) \
		&& timeunit=$$(grep TIMEUNIT MOM_input | cut -d '!' -f 1 | cut -d '=' -f 2 | xargs) \
		&& if [ -z "$${timeunit}" ]; then timeunit="8.64e4"; fi \
		&& printf -v timeunit_int "%.f" "$${timeunit}" \
		&& halfperiod=$$(printf "%.f" $$(bc <<< "scale=10; 0.5 * $${daymax} * $${timeunit_int}")) \
		&& f90nml -g ocean_solo_nml -v seconds=$${halfperiod} input.nml > tmp.nml \
		&& mv tmp.nml input.nml \
		&& echo $${daymax} $${timeunit}
	# Run the first half-period
	cd $(@D) && $(MPIRUN) -n 1 ../../$< 2> debug.out
	# Setup the next inputs
	rm -rf $(@D)/INPUT && mv $(@D)/RESTART $(@D)/INPUT
	mkdir $(@D)/RESTART
	cd $(@D) && f90nml -g mom_input_nml -v input_filename='r' input.nml > tmp.nml
	mv $(@D)/tmp.nml $(@D)/input.nml
	# Run the second half-period
	cd $(@D) && $(MPIRUN) -n 1 ../../$< 2> debug.out
	# Archive the results and cleanup
	cp $(@D)/ocean.stats $@
	rm -rf $(@D)/INPUT
	git checkout $(@D)/input.nml

# TODO: Restart checksum diagnostics


#----
.PHONY: clean
clean: clean.stats
	rm -rf build

.PHONY: clean.stats
clean.stats:
	find .testing -name ocean.stats* -exec rm {} \;
	find .testing -name chksum_diag* -exec rm {} \;
