DEPS = deps
SHELL = bash
MPIRUN ?= mpirun # srun in Slurm, aprun in Cray-land, etc.

# GFDL build toolchain
MKMF_URL ?= https://github.com/NOAA-GFDL/mkmf.git
MKMF_COMMIT ?= master

LIST_PATHS = $(abspath $(DEPS)/mkmf/bin/list_paths)
MKMF = $(abspath $(DEPS)/mkmf/bin/mkmf)

# FMS framework
FMS_URL ?= https://github.com/NOAA-GFDL/FMS.git
FMS_COMMIT ?= f2e2c86f6c0eb6d389a20509a8a60fa22924e16b
FMS := $(DEPS)/fms

# Function to get list of files
# TODO: .h and .c files
#FMS_FILES = $(sort $(shell find $(FMS)/src -name '*.F90'))
#MOM_FILES = $(sort $(shell find src config_src/solo_driver -name '*.F90'))

# Build settings
MKMF_CPP = "-Duse_libMPI -Duse_netCDF -DSPMD"

# Environment
# TODO: This info ought to be determined by CMake, automake, etc.
MKMF_TEMPLATE = "linux-ubuntu-xenial-gnu.mk"
#MKMF_TEMPLATE = "ncrc-gnu.mk"
#MKMF_TEMPLATE = "ncrc-intel.mk"

#-------------------

.PHONY: all
all: MOM6

# TODO: Split into libmom6 and executable?
MOM6: build/Makefile $(FMS)/lib/libfms.a
	make -C build NETCDF=3 DEBUG=1 COVERAGE=1 ../$@

build/Makefile: build/path_names
	cp .testing/$(MKMF_TEMPLATE) $(@D)
	cd $(@D) && $(MKMF) \
		-t $(MKMF_TEMPLATE) \
		-o '-I ../$(FMS)/build' \
		-p ../MOM6 \
		-l '../$(FMS)/lib/libfms.a' \
		-c $(MKMF_CPP) \
		$(notdir $<)

build/path_names: src $(MOM_FILES) $(LIST_PATHS)
	mkdir -p $(@D)
	cd $(@D) && $(LIST_PATHS) -l \
		../src \
		../config_src/solo_driver \
		../config_src/dynamic_symmetric

$(FMS)/lib/libfms.a: $(FMS)/build/Makefile
	mkdir -p $(FMS)/lib
	cd $(FMS)/build && make NETCDF=3 DEBUG=1 ../lib/libfms.a

$(FMS)/build/Makefile: $(FMS)/build/path_names
	cp .testing/$(MKMF_TEMPLATE) $(@D)
	cd $(@D) && $(MKMF) \
		-t $(MKMF_TEMPLATE) \
		-p ../lib/libfms.a \
		-c $(MKMF_CPP) \
		$(notdir $<)

$(FMS)/build/path_names: $(FMS)/src $(FMS_FILES) $(LIST_PATHS)
	mkdir -p $(@D)
	cd $(@D) && $(LIST_PATHS) -l ../src

$(FMS)/src:
	git clone $(FMS_URL) $@
	cd $@; git checkout $(FMS_COMMIT)

$(LIST_PATHS) $(MKMF):
	git clone $(MKMF_URL) $(DEPS)/mkmf
	cd $(DEPS)/mkmf; git checkout $(MKMF_COMMIT)

#----
TESTS = benchmark unit_tests circle_obcs

.PHONY: test $(TESTS)
test: $(TESTS)

$(TESTS):
	find $(BUILD_PATH) -name *.gcda -exec rm -f '{}' \;
	mkdir -p .testing/$@/RESTART
	cd .testing/$@ && $(MPIRUN) -n 1 ../../MOM6
	bash <(curl -s https://codecov.io/bash) -n $@
