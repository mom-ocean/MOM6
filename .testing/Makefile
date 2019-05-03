# Makefile steps to run on Travis-CI
# e.g. make MEMORY_SHAPE=dynamic_symmetric REPRO=1 OPENMP=1

# Versions to use
FMS_COMMIT ?= xanadu
MKMF_COMMIT ?= master

# Where to clone from
FMS_URL ?= https://github.com/NOAA-GFDL/FMS.git
MKMF_URL ?= https://github.com/NOAA-GFDL/mkmf.git
CONFIGS_URL ?= https://github.com/NOAA-GFDL/MOM6-examples.git
REGRESSIONS_URL ?= https://github.com/adcroft/Gaea-stats-MOM6-examples

# Experiments to run
ifeq ($(MEMORY_SHAPE),"dynamic_symmetric")
EXPERIMENTS ?= unit_tests double_gyre flow_downslope/z CVmix_SCM_tests/cooling_only/EPBL circle_obcs
else
EXPERIMENTS ?= unit_tests double_gyre flow_downslope/z CVmix_SCM_tests/cooling_only/EPBL
endif

FMS_PACKAGES ?= platform,include,memutils,constants,mpp,fms,time_manager,diag_manager,data_override,coupler/coupler_types.F90,coupler/ensemble_manager.F90,axis_utils,horiz_interp,time_interp,astronomy,mosaic,random_numbers
TEMPLATE ?= .testing/linux-ubuntu-xenial-gnu.mk
MPIRUN ?= mpirun

# MEMORY_SHAPE must be defined for this Makefile to work
MEMORY_SHAPE ?= dynamic_symmetric

# Everything above is above is "configurable" with environment variables
SHELL = bash

# Path where executable will be built
BUILD_PATH = build
###/$(MEMORY_SHAPE)-$(EXEC_MODE)
# Root of configurations (MOM6-examples)
EXPERIMENTS_ROOT = experiments
# Regression results
REGRESSIONS_ROOT = answers

.PRECIOUS: %/ocean.stats

run: $(foreach e,$(EXPERIMENTS),$(EXPERIMENTS_ROOT)/ocean_only/$(e)/ocean.stats)

test: $(foreach e,$(EXPERIMENTS),$(REGRESSIONS_ROOT)/regressions/ocean_only/$(e)/ocean.stats.gnu)

compile: $(BUILD_PATH)/MOM6

$(BUILD_PATH)/MOM6: FMS mkmf
	mkdir -p $(@D)
	cd $(@D); \
	../mkmf/bin/list_paths -l ../FMS/{$(FMS_PACKAGES)} ../config_src/{$(MEMORY_SHAPE),solo_driver} ../src \
	&& ../mkmf/bin/mkmf -t ../$(TEMPLATE) -c '-Duse_libMPI -Duse_netCDF -DSPMD -DUSE_LOG_DIAG_FIELD_INFO -DMAXFIELDMETHODS_=500' -p $(@F) path_names \
	&& make -j NETCDF=3 $(@F)

$(EXPERIMENTS_ROOT)/%/ocean.stats: $(EXPERIMENTS_ROOT)
	mkdir -p $(@D)/RESTART
	cd $(@D) ; $(MPIRUN) -n 1 $(PWD)/$(BUILD_PATH)/MOM6

$(REGRESSIONS_ROOT)/regressions/%/ocean.stats.gnu: $(EXPERIMENTS_ROOT)/%/ocean.stats $(REGRESSIONS_ROOT)
	cp $< $@
	cd $(@D) ; git status --porcelain $(@F)

# Targets to clone repositories needed to build
FMS:
	git clone -q $(FMS_URL)
	cd $@ ; git checkout -q $(FMS_COMMIT)

mkmf:
	git clone -q $(MKMF_URL)
	cd $@ ; git checkout -q $(MKMF_COMMIT)

$(EXPERIMENTS_ROOT):
	mkdir -p $(@D)
	cd $(@D) ; git clone --depth 1 $(CONFIGS_URL) experiments

$(REGRESSIONS_ROOT):
	mkdir -p $(@D)
	cd $(@D) ; git clone --depth 1 $(REGRESSIONS_URL) answers
