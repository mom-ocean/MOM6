SHELL = bash
.SUFFIXES:

# FMS build configuration
PATH_ENV ?=
FCFLAGS_ENV ?=

# Only set FCFLAGS if an argument is provided.
FMS_FCFLAGS =
ifneq ($(FCFLAGS_ENV),)
  FMS_FCFLAGS := FCFLAGS="$(FCFLAGS_ENV)"
endif


# Ditto for path
FMS_PATH =
ifneq ($(PATH_ENV),)
  FMS_PATH := PATH="$(PATH_ENV)"
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
lib/libFMS.a: fms/build/libFMS.a fms/build/Makefile
	mkdir -p {lib,include}
	cp fms/build/libFMS.a lib/libFMS.a
	cp fms/build/*.mod include


fms/build/libFMS.a: fms/build/Makefile
	make -C fms/build libFMS.a


# TODO: Include FC, CC, CFLAGS?
fms/build/Makefile: FMS_ENV=$(FMS_PATH) $(FMS_FCFLAGS)

fms/build/Makefile: Makefile.fms.in fms/src/configure bin/mkmf bin/list_paths
	mkdir -p fms/build
	cp Makefile.fms.in fms/src/Makefile.in
	cd $(@D) && $(FMS_ENV) ../src/configure --srcdir=../src


# TODO: Track m4 macros?
fms/src/configure: configure.fms.ac $(FMS_SOURCE) | fms/src
	cp configure.fms.ac fms/src/configure.ac
	cp -r m4 $(@D)
	cd $(@D) && autoreconf -i

fms/src:
	git clone $(FMS_URL) $@
	git -C $@ checkout $(FMS_COMMIT)


.PHONY: clean
clean:
	rm -rf fms/build lib include bin


.PHONY: distclean
distclean: clean
	rm -rf fms mkmf
