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
	make -C fms/build libFMS.a


fms/build/Makefile: Makefile.fms.in fms/src/configure
	mkdir -p fms/build
	cp Makefile.fms.in fms/src/Makefile.in
	cd $(@D) && { \
	  ../src/configure --srcdir=../src \
	  || { \
	    if [ "${REPORT_ERROR_LOGS}" = true ]; then cat config.log ; fi ; \
	    false; \
	  } \
	}


fms/src/configure: configure.fms.ac $(FMS_SOURCE) | fms/src
	cp configure.fms.ac fms/src/configure.ac
	cp -r m4 $(@D)
	cd $(@D) && autoreconf -i


fms/src:
	git clone $(FMS_URL) $@
	git -C $@ checkout $(FMS_COMMIT)


.PHONY: clean
clean:
	rm -rf fms/build lib include

.PHONY: distclean
distclean: clean
	rm -rf fms
