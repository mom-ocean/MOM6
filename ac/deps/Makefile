# This file is part of MOM6, the Modular Ocean Model version 6.
# See the LICENSE file for licensing information.
# SPDX-License-Identifier: Apache-2.0

SHELL = bash

# Disable implicit rules
.SUFFIXES:

# Disable implicit variables
MAKEFLAGS += -R

# FMS framework
FMS_URL ?= https://github.com/NOAA-GFDL/FMS.git
FMS_COMMIT ?= 2023.03

GSW_URL ?= https://github.com/mom-ocean/GSW-Fortran.git
GSW_COMMIT ?= 29e64d652786e1d076a05128c920f394202bfe10

CVMIX_URL ?= https://github.com/mom-ocean/CVMix-src.git
CVMIX_COMMIT ?= 65ef5c73bc7f5663d5688f75c3855d431da4baea


# List of source files to link this Makefile's dependencies to model Makefiles
# Assumes a depth of two, and the following extensions: F90 inc c h
# (1): Root directory
SOURCE = \
  $(foreach ext,F90 inc c h,$(wildcard $(1)/*/*.$(ext) $(1)/*/*/*.$(ext)))

FMS_SOURCE = $(call SOURCE,fms/src)
GSW_SOURCE = $(call SOURCE,gsw/src)
CVMIX_SOURCE = $(call SOURCE,CVMix-src/src/shared)


# If `true`, print logs if an error is encountered.
REPORT_ERROR_LOGS ?=


# If set, use the submodule repositories in pkg/
PKG ?= $(abspath ../../pkg)
MAKEDEP ?= $(abspath ../makedep)

#---
# Rules

.PHONY: all
all: lib/libFMS.a
all: lib/libgsw.a
all: lib/libcvmix.a


# Library build rules template
#
# $(1): target library
# $(2): dependency label
# $(3): library source files

define LIB_RULES
lib/$(1): $(2)/build/$(1)
	mkdir -p $$(@D) include/
	cp $$< $$@
	cp $$(dir $$<)/*.mod include/

$(2)/build/$(1): $(2)/build/Makefile
	$$(MAKE) -C $$(@D) $(1)

$(2)/build/Makefile: $(2)/build/Makefile.in $(2)/build/configure
	cd $$(@D) && { \
	  MAKEDEP=$$(MAKEDEP) \
	  ./configure --srcdir=../src \
	  || { \
	    if [ "$${REPORT_ERROR_LOGS}" = true ]; then cat config.log ; fi ; \
	    false; \
	  } \
	}

$(2)/build/Makefile.in: Makefile.$(2).in | $(2)/build
	cp $$< $$@

$(2)/build/configure: $(2)/build/configure.ac $(3) | $(2)/src
	autoreconf $$(@D)

$(2)/build/configure.ac: configure.$(2).ac m4 | $(2)/build
	cp $$< $$@
	cp -r m4 $$(@D)

$(2)/build:
	mkdir -p $$@
endef

$(eval $(call LIB_RULES,libFMS.a,fms,$(FMS_SOURCE)))
$(eval $(call LIB_RULES,libgsw.a,gsw,$(GSW_SOURCE)))
$(eval $(call LIB_RULES,libcvmix.a,cvmix,$(CVMIX_SOURCE)))


# Dependency source

fms/src:
	git clone $(FMS_URL) $@
	git -C $@ checkout $(FMS_COMMIT)


ifdef PKG
gsw/src: | gsw/build
	ln -s $(PKG)/GSW-Fortran gsw/src

cvmix/src: | cvmix/build
	ln -s $(PKG)/CVMix-src cvmix/src

else
gsw/src:
	git clone $(GSW_URL) $@
	git -C $@ checkout $(GSW_COMMIT)

cvmix/src:
	git clone $(CVMIX_URL) $@
	git -C $@ checkout $(CVMIX_COMMIT)
endif


# Cleanup

.PHONY: clean
clean:
	rm -rf fms/build gsw/build cvmix/build lib include

.PHONY: distclean
distclean: clean
	rm -rf fms gsw cvmix
