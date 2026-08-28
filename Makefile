# This file is part of MOM6, the Modular Ocean Model version 6.
# See the LICENSE file for licensing information.
# SPDX-License-Identifier: Apache-2.0

BUILD ?= build
FMS_BUILD ?= ac/deps/fms/build
MOM_MEMORY ?=

.PHONY: all
all: $(BUILD)/MOM6

$(BUILD)/MOM6: $(BUILD)/Makefile $(FMS_BUILD)/libFMS.a
	if test $(FMS_BUILD)/libFMS.a -nt $@ ; then \
	  $(MAKE) -C $(BUILD) clean ; \
	fi
	$(MAKE) -C $(BUILD) MOM6


# Makefile setup

$(BUILD)/Makefile: $(BUILD)/config.status ac/Makefile.in
	cd $(BUILD) && ./config.status

$(BUILD)/config.status: configure $(FMS_BUILD)/libFMS.a | $(BUILD)
	cd $(BUILD) && \
	PATH="${PATH}:$(CURDIR)/ac" \
	$(CURDIR)/configure -n $(CONFIG_FLAGS)


# configure setup

configure: ac/configure.ac ac/aclocal.m4 | ac/m4
	cd ac && autoconf -o $(abspath $@)

ac/aclocal.m4: ac/configure.ac | ac/m4
	cd ac && aclocal

$(BUILD):
	mkdir -p $@


# Dependencies

# NOTE: If libFMS has changed, then we completely rebuild MOM6
$(FMS_BUILD)/libFMS.a: FORCE
	$(MAKE) -C ac/deps \
	  BUILD=$(abspath $(FMS_BUILD)) \
	  CODEBASE=$(abspath $(FMS_CODEBASE))

FORCE:


# Cleanup

# Remove build output
.PHONY: clean
clean:
	rm -rf $(BUILD)
	$(MAKE) -C ac/deps clean

# Remove generated autoconf output
.PHONY: ac-clean
ac-clean:
	rm -f configure
	rm -f configure~
	rm -f ac/aclocal.m4
	rm -rf ac/autom4te.cache/

# Remove all build products
.PHONY: maintainer-clean
maintainer-clean: clean ac-clean
