#
# SEATREE top-level build
#
# Builds the compiled components of the SEATREE modules:
#
#   make all        hc conman syn2d larry larry3d nll (default)
#   make optional   eqdyna hypoinverse (need MPI / are not used by the GUI)
#   make test       HC regression tests and a short ConMan run
#   make clean      remove objects and binaries of everything above
#
# Individual targets: hc, conman, syn2d, larry, larry3d, nll, eqdyna,
# hypoinverse. The python installer (install/configure.python3.gtk4)
# calls these targets and then writes the configuration files.
#
# HC is built without its GMT bindings (NO_GMT). SEATREE only uses the
# ASCII output of the HC tools and converts it with the system GMT,
# so no GMT 4 development files are needed. To build the GMT 4 based
# HC tools as well, run make inside modules/mc/hc with GMT4HOME and
# NETCDFHOME set, see modules/mc/hc/README.md.
#
HC_DIR = modules/mc/hc
CONMAN_DIR = modules/mc/ConMan
SYN2D_DIR = modules/seismo/syn2d
LARRY_DIR = modules/seismo/larry
LARRY3D_DIR = modules/seismo/larry3d
NLL_DIR = modules/seismo/nonlinloc/nll
EQDYNA_DIR = modules/seismo/EQdyna
HYPO_DIR = modules/seismo/hypoinverse/src_hyp2000g

#
# compilers. The module makefiles use $(CC), $(FC), $(F77) and $(F90);
# all are passed down explicitly so that they can be set on the
# command line, e.g.
#
#   make FC=/usr/bin/gfortran all
#
# make's built-in defaults (cc, f77) are replaced by the full paths of
# gcc and gfortran unless CC or FC were given on the command line or
# in the environment. The full path is used because GNU make 4.3
# fails with "Permission denied" when it runs a compiler directly and
# a directory of the same name sits in an earlier PATH entry (for
# instance with "." in the path), whereas the shell's command -v skips
# such entries.
#
ifeq ($(origin CC),default)
CC = $(shell command -v gcc 2>/dev/null || echo gcc)
endif
ifeq ($(origin FC),default)
FC = $(shell command -v gfortran 2>/dev/null || echo gfortran)
endif
# the seismology makefiles distinguish F77 and F90; both default to FC
F77 ?= $(FC)
F90 ?= $(FC)

# flags for the fixed-form Fortran 77 sources of the seismology
# modules (gfortran). F_EXT_SOURCE_FLAG allows 132 column lines.
CFLAGS ?= -O2 -fcommon
FFLAGS ?= -O2 -fallow-argument-mismatch
F_EXT_SOURCE_FLAG ?= -ffixed-line-length-132

# larry and larry3d put their binaries in a subdirectory named after
# the architecture; the python side (conf/larry*/..Conf.xml) uses the
# same name
ARCH ?= $(shell uname -m)

# variables handed to the seismology makefiles
SEISMO_VARS = CC="$(CC)" F77="$(F77)" F90="$(F90)" CFLAGS="$(CFLAGS)" \
	FFLAGS="$(FFLAGS)" F_EXT_SOURCE_FLAG="$(F_EXT_SOURCE_FLAG)" ARCH="$(ARCH)"

# HC places binaries in bin/$(ARCH) if ARCH is set. We build into
# plain bin/ so that the python side finds the binaries regardless of
# the run-time environment. With NO_GMT, the HC Makefile (7bb6d03 or
# later) includes sh_syn in the all_no_gmt target.
HC_MAKE = $(MAKE) -C $(HC_DIR) ARCH= CC=$(CC)
HC_FLAGS = ADD_FLAGS=-DNO_GMT

# EQdyna's makefile selects compiler and NetCDF locations by MACHINE
UNAME_S = $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
EQDYNA_MACHINE = macos
else
EQDYNA_MACHINE = ubuntu
endif

all: hc conman syn2d larry larry3d nll

optional: eqdyna hypoinverse

#
# fetch the submodules if the checkout did not use --recurse-submodules
#
submodules:
	@if [ ! -f $(HC_DIR)/Makefile ] || [ ! -f $(CONMAN_DIR)/src/Makefile ] || \
	    [ ! -f $(EQDYNA_DIR)/src/makefile ]; then \
		echo "fetching git submodules"; \
		git submodule update --init; \
	fi

#
# mantle convection
#
hc: submodules
	$(HC_MAKE) $(HC_FLAGS) all_no_gmt
	@echo "HC binaries are in $(HC_DIR)/bin"

conman: submodules
	$(MAKE) -C $(CONMAN_DIR)/src FC=$(FC)
	@echo "ConMan binary is $(CONMAN_DIR)/conman"

#
# seismology
#
syn2d:
	mkdir -p $(SYN2D_DIR)/makemodel/bin $(SYN2D_DIR)/makedata/bin $(SYN2D_DIR)/inversion/bin
	$(MAKE) -C $(SYN2D_DIR) $(SEISMO_VARS) build_only
	@echo "syn2d binaries are in $(SYN2D_DIR)/{makemodel,makedata,inversion}/bin"

larry:
	$(MAKE) -C $(LARRY_DIR) $(SEISMO_VARS)
	@echo "larry binaries are in $(LARRY_DIR)/$(ARCH)"

larry3d:
	$(MAKE) -C $(LARRY3D_DIR) $(SEISMO_VARS)
	@echo "larry3d binaries are in $(LARRY3D_DIR)/$(ARCH)"

nll:
	mkdir -p $(NLL_DIR)/bin
	$(MAKE) -C $(NLL_DIR)/src CC="$(CC)" CFLAGS="$(CFLAGS)"
	@echo "NonLinLoc binaries are in $(NLL_DIR)/bin"

#
# optional: EQdyna needs an MPI Fortran compiler (mpif90) and the
# NetCDF Fortran library; its makefile picks them by MACHINE
# (ubuntu or macos). hypoinverse is not used by the GUI.
#
eqdyna: submodules
	$(MAKE) -C $(EQDYNA_DIR)/src MACHINE=$(EQDYNA_MACHINE)
	mkdir -p $(EQDYNA_DIR)/bin
	mv $(EQDYNA_DIR)/src/eqdyna $(EQDYNA_DIR)/bin/
	@echo "EQdyna binary is $(EQDYNA_DIR)/bin/eqdyna"

hypoinverse:
	mkdir -p $(HYPO_DIR)/../bin
	$(MAKE) -C $(HYPO_DIR) F77="$(F77)" FFLAGS="$(FFLAGS) $(F_EXT_SOURCE_FLAG)"
	@echo "hypoinverse binary is $(HYPO_DIR)/../bin/hyp"

#
# tests
#
test: hc-test conman-test

hc-test: hc
	cd $(HC_DIR)/test && ./run_test.sh ../bin compare 1e-6 "-solver 0"
	cd $(HC_DIR)/test && ./run_test.sh ../bin compare 1e-6 "-solver 1"

conman-test: conman
	./install/test.conman.sh

clean:
	-$(HC_MAKE) clean
	-$(MAKE) -C $(CONMAN_DIR)/src FC=$(FC) clean
	-$(MAKE) -C $(SYN2D_DIR) $(SEISMO_VARS) clean_all
	-rm -rf $(LARRY_DIR)/$(ARCH) $(LARRY3D_DIR)/$(ARCH)
	-rm -f $(NLL_DIR)/src/*.o $(NLL_DIR)/bin/*
	-$(MAKE) -C $(EQDYNA_DIR)/src MACHINE=$(EQDYNA_MACHINE) clean
	-rm -rf $(EQDYNA_DIR)/bin
	-$(MAKE) -C $(HYPO_DIR) clean

.PHONY: all optional submodules hc conman syn2d larry larry3d nll eqdyna hypoinverse \
	test hc-test conman-test clean
