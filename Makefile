#
# SEATREE top-level build
#
# This builds the compiled components that the mantle convection
# modules need:
#
#   make hc        HC mantle flow code (git submodule modules/mc/hc)
#   make conman    ConMan 2-D convection code (git submodule modules/mc/ConMan)
#   make all       both of the above (default)
#   make test      run the HC regression tests and a short ConMan run
#   make clean     remove the objects and binaries of both codes
#
# The seismology modules (syn2d, larry, larry3d, nonlinloc, EQdyna)
# are still built by the python installer, see install.seatree.sh.
#
# HC is built without its GMT bindings (NO_GMT). SEATREE only uses the
# ASCII output of the HC tools and converts it with the system GMT,
# so no GMT 4 development files are needed. To build the GMT 4 based
# HC tools as well, run make inside modules/mc/hc with GMT4HOME and
# NETCDFHOME set, see modules/mc/hc/README.md.
#
HC_DIR = modules/mc/hc
CONMAN_DIR = modules/mc/ConMan

# HC places binaries in bin/$(ARCH) if ARCH is set. We build into
# plain bin/ so that the python side finds the binaries regardless of
# the run-time environment. With NO_GMT, the HC Makefile (7bb6d03 or
# later) includes sh_syn in the all_no_gmt target.
HC_MAKE = $(MAKE) -C $(HC_DIR) ARCH=
HC_FLAGS = ADD_FLAGS=-DNO_GMT

all: hc conman

#
# fetch the submodules if the checkout did not use --recurse-submodules
#
submodules:
	@if [ ! -f $(HC_DIR)/Makefile ] || [ ! -f $(CONMAN_DIR)/src/Makefile ]; then \
		echo "fetching git submodules"; \
		git submodule update --init; \
	fi

hc: submodules
	$(HC_MAKE) $(HC_FLAGS) all_no_gmt
	@echo "HC binaries are in $(HC_DIR)/bin"

conman: submodules
	$(MAKE) -C $(CONMAN_DIR)/src
	@echo "ConMan binary is $(CONMAN_DIR)/conman"

test: hc-test conman-test

hc-test: hc
	cd $(HC_DIR)/test && ./run_test.sh ../bin compare 1e-6 "-solver 0"
	cd $(HC_DIR)/test && ./run_test.sh ../bin compare 1e-6 "-solver 1"

conman-test: conman
	./install/test.conman.sh

clean:
	-$(HC_MAKE) clean
	-$(MAKE) -C $(CONMAN_DIR)/src clean

.PHONY: all submodules hc conman test hc-test conman-test clean
