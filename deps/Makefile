SHELL = bash
.SUFFIXES:

# FMS build configuration
# TODO: Include FC?
FCFLAGS_FMS ?=

# Only set FCFLAGS if an argument is provided.
FCFLAGS_ENV =
ifneq ($(FCFLAGS_FMS),)
  FCFLAGS_ENV := FCFLAGS="$(FCFLAGS_FMS)"
endif


# mkmf, list_paths (GFDL build toolchain)
MKMF_URL ?= https://github.com/NOAA-GFDL/mkmf.git
MKMF_COMMIT ?= master

# FMS framework
FMS_URL ?= https://github.com/NOAA-GFDL/FMS.git
FMS_COMMIT ?= 2019.01.03


# List of source files to link this Makefile's dependencies to model Makefiles
# Assumes a depth of two, and the following extensions: F90 inc c h
# (1): Root directory
# NOTE: extensions could be a second variable
SOURCE = \
  $(foreach ext,F90 inc c h,$(wildcard $(1)/*/*.$(ext) $(1)/*/*/*.$(ext)))

FMS_SOURCE = $(call SOURCE,fms/src)


#---
# Rules

.PHONY: all
all: bin/mkmf bin/list_paths lib/libFMS.a


#---
# mkmf checkout

bin/mkmf bin/list_paths: mkmf
	mkdir -p $(@D)
	cp $^/$@ $@

mkmf:
	git clone $(MKMF_URL) $@
	git -C $@ checkout $(MKMF_COMMIT)


#---
# FMS build

# NOTE: We emulate the automake `make install` stage by storing libFMS.a to
#   ${srcdir}/deps/lib and copying module files to ${srcdir}/deps/include.
#   This is a flawed approach, since module files are untracked and could be
#   handled more safely, but this is adequate for now.

# TODO: track *.mod copy?
lib/libFMS.a: fms/build/Makefile
	mkdir -p {lib,include}
	cd fms/build && $(MAKE) -j ../../lib/libFMS.a
	cp fms/build/*.mod include


# TODO: Include FC, CC, CFLAGS?
fms/build/Makefile: FC_ENV=$(FCFLAGS_ENV)

fms/build/Makefile: fms/src/configure bin/mkmf bin/list_paths
	mkdir -p fms/build
	cd $(@D) && $(FC_ENV) ../src/configure --srcdir=../src


# TODO: Track m4 macros?
fms/src/configure: fms/src configure.fms.ac Makefile.fms.in $(FMS_SOURCE)
	cp configure.fms.ac $</configure.ac
	cp Makefile.fms.in $</Makefile.in
	cp -r ../m4 $(@D)
	cd $(@D) && autoreconf -i


fms/src:
	git clone $(FMS_URL) $@
	git -C $@ checkout $(FMS_COMMIT)


.PHONY: clean
clean:
	rm -rf fms/build lib include bin

distclean: clean
	rm -rf fms mkmf
