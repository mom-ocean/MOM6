SHELL = bash

# Disable implicit rules
.SUFFIXES:

# Disable implicit variables
MAKEFLAGS += -R

# FMS framework
FMS_URL ?= https://github.com/NOAA-GFDL/FMS.git
FMS_COMMIT ?= 2019.01.03


# List of source files to link this Makefile's dependencies to model Makefiles
# Assumes a depth of two, and the following extensions: F90 inc c h
# (1): Root directory
SOURCE = \
  $(foreach ext,F90 inc c h,$(wildcard $(1)/*/*.$(ext) $(1)/*/*/*.$(ext)))

FMS_SOURCE = $(call SOURCE,fms/src)


# If `true`, print logs if an error is encountered.
REPORT_ERROR_LOGS ?=


#---
# Rules

.PHONY: all
all: lib/libFMS.a


#---
# FMS build

# NOTE: We emulate the automake `make install` stage by storing libFMS.a to
#   ${srcdir}/deps/lib and copying module files to ${srcdir}/deps/include.
lib/libFMS.a: fms/build/libFMS.a
	mkdir -p lib include
	cp fms/build/libFMS.a lib/libFMS.a
	cp fms/build/*.mod include

fms/build/libFMS.a: fms/build/Makefile
	$(MAKE) -C fms/build libFMS.a

fms/build/Makefile: fms/build/Makefile.in fms/build/configure
	cd $(@D) && { \
	  ./configure --srcdir=../src \
	  || { \
	    if [ "${REPORT_ERROR_LOGS}" = true ]; then cat config.log ; fi ; \
	    false; \
	  } \
	}

fms/build/Makefile.in: Makefile.fms.in | fms/build
	cp Makefile.fms.in fms/build/Makefile.in

fms/build/configure: fms/build/configure.ac $(FMS_SOURCE) | fms/src
	autoreconf fms/build

fms/build/configure.ac: configure.fms.ac m4 | fms/build
	cp configure.fms.ac fms/build/configure.ac
	cp -r m4 fms/build

fms/build:
	mkdir -p fms/build

fms/src:
	git clone $(FMS_URL) $@
	git -C $@ checkout $(FMS_COMMIT)

# Cleanup

.PHONY: clean
clean:
	rm -rf fms/build lib include

.PHONY: distclean
distclean: clean
	rm -rf fms
